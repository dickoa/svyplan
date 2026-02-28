#' Estimate Variance Components
#'
#' Estimate between- and within-stage variance components using nested
#' ANOVA decomposition. Supports SRS and PPS first-stage designs.
#'
#' @param x A formula, numeric vector, or survey design object
#'   (see Details).
#' @param ... Additional arguments passed to methods.
#'
#' @return A `svyplan_varcomp` object with components:
#' \describe{
#'   \item{`varb`}{Between-PSU variance (scalar).}
#'   \item{`varw`}{Within-PSU variance. Scalar for 2-stage, length-2
#'     vector (`varw_psu`, `varw_ssu`) for 3-stage.}
#'   \item{`delta`}{Measure of homogeneity. Length 1 for 2-stage,
#'     length 2 (`delta_psu`, `delta_ssu`) for 3-stage.}
#'   \item{`k`}{Ratio parameter(s), same length as `delta`
#'     (`k_psu`, `k_ssu` for 3-stage).}
#'   \item{`rel_var`}{Unit relvariance (scalar).}
#'   \item{`stages`}{Number of stages (2 or 3).}
#' }
#'
#' @details
#' The interface is determined by the class of `x`:
#'
#' - **Formula**: `varcomp(income ~ district, data = frame)`. The LHS is
#'   the analysis variable, RHS terms are stage IDs (outermost first).
#' - **Numeric vector**: `varcomp(y, stage_id = list(cluster_ids))`.
#' - **survey.design**: `varcomp(design, ~y)`. Cluster structure is
#'   extracted from the design object. Requires the survey package.
#'
#' When `prob` is `NULL`, SRS first-stage is assumed. When provided, PPS
#' variance estimation is used.
#'
#' The returned `delta` is the measure of homogeneity
#' \eqn{\delta = V_b / (V_b + V_w)}{delta = Vb / (Vb + Vw)} following
#' Valliant, Dever, and Kreuter (2018, Ch. 9). Unlike the traditional ANOVA
#' intraclass correlation coefficient, `delta` is constrained to \eqn{[0, 1]}
#' and should not be compared directly to mixed-model ICCs (e.g. from lme4)
#' which can be negative.
#'
#' Clusters containing a single observation have undefined within-cluster
#' variance. In this case, the within-cluster variance is imputed as the
#' mean variance of the remaining clusters.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 9.
#'
#' Hansen, M. H., Hurwitz, W. N., and Madow, W. G. (1953).
#' *Sample Survey Methods and Theory* (Vol. I). Wiley.
#'
#' @seealso [n_cluster()] which accepts a `svyplan_varcomp` as `delta`.
#'
#' @examples
#' # 2-stage SRS using formula
#' set.seed(42)
#' frame <- data.frame(
#'   income = rnorm(200, 50000, 10000),
#'   district = rep(1:20, each = 10)
#' )
#' vc <- varcomp(income ~ district, data = frame)
#' vc
#'
#' # Feed into n_cluster
#' n_cluster(cost = c(500, 50), delta = vc, budget = 100000)
#'
#' @export
varcomp <- function(x, ...) {
  UseMethod("varcomp")
}

#' @describeIn varcomp Method for formula interface.
#'
#' @param data A data frame (required for formula interface).
#' @param prob First-stage selection probabilities. A one-sided formula
#'   (e.g., `~pp`) when using the formula interface, or a numeric vector.
#'   `NULL` (default) assumes SRS.
#'
#' @export
varcomp.formula <- function(x, ..., data = NULL, prob = NULL) {
  .varcomp_formula(x, data = data, prob = prob)
}

#' @describeIn varcomp Default method for numeric vectors.
#'
#' @param stage_id A list of stage-ID vectors (required for vector interface).
#'   Length determines the number of stage boundaries (stages - 1).
#'
#' @export
varcomp.default <- function(x, ..., stage_id = NULL, prob = NULL) {
  if (is.numeric(x)) {
    .varcomp_vector(x, stage_id = stage_id, prob = prob)
  } else {
    stop("'x' must be a formula, numeric vector, or survey design object",
         call. = FALSE)
  }
}

#' @describeIn varcomp Method for survey design objects. Pass a one-sided
#'   formula (e.g., `~y`) to specify the outcome variable. Cluster
#'   structure is extracted from the design.
#'
#' @export
varcomp.survey.design <- function(x, ..., prob = NULL) {
  formula <- NULL
  for (a in list(...)) {
    if (inherits(a, "formula")) {
      formula <- a
      break
    }
  }
  if (is.null(formula)) {
    stop("a one-sided formula specifying the outcome is required (e.g., ~y)",
         call. = FALSE)
  }

  y_name <- all.vars(formula)
  if (length(y_name) != 1L) {
    stop("formula must reference exactly one variable", call. = FALSE)
  }

  y <- x$variables[[y_name]]
  if (is.null(y)) {
    stop(sprintf("variable '%s' not found in design variables", y_name),
         call. = FALSE)
  }

  cl <- x$cluster
  n_stages <- ncol(cl)

  if (n_stages == 1L && length(unique(cl[[1L]])) == nrow(cl)) {
    stop("design has no clusters (ids = ~1); varcomp requires a clustered design",
         call. = FALSE)
  }

  stage_id <- lapply(seq_len(n_stages), function(j) cl[[j]])
  .varcomp_dispatch(y, stage_id, prob)
}

#' Parse formula interface and dispatch
#' @keywords internal
#' @noRd
.varcomp_formula <- function(formula, data, prob) {
  if (is.null(data)) {
    stop("'data' is required for the formula interface", call. = FALSE)
  }

  vars <- all.vars(formula)
  if (length(vars) < 2L) {
    stop("formula must have a response and at least one stage ID (e.g., y ~ cluster)",
         call. = FALSE)
  }

  y_name <- vars[1L]
  stage_names <- vars[-1L]

  y <- data[[y_name]]
  if (is.null(y)) {
    stop(sprintf("variable '%s' not found in 'data'", y_name), call. = FALSE)
  }

  stage_id <- lapply(stage_names, function(nm) {
    col <- data[[nm]]
    if (is.null(col)) {
      stop(sprintf("variable '%s' not found in 'data'", nm), call. = FALSE)
    }
    col
  })

  pp <- NULL
  if (!is.null(prob)) {
    if (inherits(prob, "formula")) {
      prob_name <- all.vars(prob)
      if (length(prob_name) != 1L) {
        stop("'prob' formula must reference exactly one variable", call. = FALSE)
      }
      pp <- data[[prob_name]]
      if (is.null(pp)) {
        stop(sprintf("variable '%s' not found in 'data'", prob_name),
             call. = FALSE)
      }
    } else {
      pp <- prob
    }
  }

  .varcomp_dispatch(y, stage_id, pp)
}

#' Vector interface
#' @keywords internal
#' @noRd
.varcomp_vector <- function(x, stage_id, prob) {
  if (is.null(stage_id) || !is.list(stage_id)) {
    stop("'stage_id' must be a list of ID vectors", call. = FALSE)
  }
  .varcomp_dispatch(x, stage_id, prob)
}

#' Dispatch to correct variance component estimator
#' @keywords internal
#' @noRd
.varcomp_dispatch <- function(y, stage_id, prob) {
  if (!is.numeric(y) || length(y) == 0L) {
    stop("'y' must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("outcome vector must not contain NA values", call. = FALSE)
  }
  if (length(stage_id) == 0L) {
    stop("'stage_id' must not be empty", call. = FALSE)
  }
  for (i in seq_along(stage_id)) {
    if (length(stage_id[[i]]) != length(y)) {
      stop(sprintf("'stage_id[[%d]]' must have length %d (same as outcome vector)",
                   i, length(y)), call. = FALSE)
    }
  }

  n_boundaries <- length(stage_id)
  stages <- n_boundaries + 1L

  if (stages > 3L) {
    stop("4+ stage variance components are not yet supported", call. = FALSE)
  }

  has_prob <- !is.null(prob)

  if (stages == 2L && !has_prob) {
    .varcomp_2stage_srs(y, stage_id[[1L]])
  } else if (stages == 2L && has_prob) {
    .varcomp_2stage_pps(y, stage_id[[1L]], prob)
  } else if (stages == 3L) {
    if (!has_prob) {
      stop("3-stage varcomp requires 'prob' (PPS probabilities)", call. = FALSE)
    }
    .varcomp_3stage_pps(y, stage_id[[1L]], stage_id[[2L]], prob)
  }
}

#' Impute singleton cluster variances (NA from var() on length-1 groups)
#' @keywords internal
#' @noRd
.impute_singleton_var <- function(x) {
  lonely <- is.na(x)
  if (all(lonely)) {
    warning("all clusters are singletons; within-cluster variance set to 0",
            call. = FALSE)
    x[] <- 0
  } else if (any(lonely)) {
    x[lonely] <- mean(x[!lonely])
  }
  x
}

#' 2-stage SRS variance components
#' @keywords internal
#' @noRd
.varcomp_2stage_srs <- function(y, psu_id) {
  uid <- unique(psu_id)
  M <- length(uid)
  idx <- match(psu_id, uid)
  Ni <- tabulate(idx, nbins = M)

  grp <- split(y, idx)
  ti <- vapply(grp, sum, numeric(1L))
  S2Ui <- vapply(grp, var, numeric(1L))

  S2Ui <- .impute_singleton_var(S2Ui)

  tbarU <- mean(ti)
  tU <- M * tbarU
  S2U1 <- var(ti)
  vb <- S2U1 / tbarU^2

  ybarU <- mean(y)
  S2U <- var(y)

  vw <- M * sum(Ni^2 * S2Ui) / tU^2

  eps <- sqrt(.Machine$double.eps)
  rel_var <- if (abs(ybarU) < eps && S2U < eps) 0 else S2U / ybarU^2

  total_v <- vb + vw
  if (total_v < eps || rel_var < eps || !is.finite(rel_var)) {
    warning("outcome variance is approximately zero; delta set to 0 by convention",
            call. = FALSE)
    delta <- 0
    k <- 1
  } else {
    delta <- vb / total_v
    k <- total_v / rel_var
  }

  .new_svyplan_varcomp(
    varb    = vb,
    varw    = vw,
    delta   = delta,
    k       = k,
    rel_var = rel_var,
    stages  = 2L
  )
}

#' 2-stage PPS variance components
#' @keywords internal
#' @noRd
.varcomp_2stage_pps <- function(y, psu_id, pp) {
  unique_psu <- sort(unique(psu_id))
  M <- length(unique_psu)
  idx <- match(psu_id, unique_psu)
  Ni <- tabulate(idx, nbins = M)

  # If pp has one value per element, extract per-PSU values
  if (length(pp) == length(y)) {
    pp_psu <- pp[match(unique_psu, psu_id)]
  } else if (length(pp) == M) {
    pp_psu <- pp
  } else {
    stop("'prob' must have length equal to number of observations or number of PSUs",
         call. = FALSE)
  }

  if (abs(sum(pp_psu) - 1) >= 1e-3) {
    stop("'prob' values must sum to 1", call. = FALSE)
  }

  grp <- split(y, idx)
  cl_tots <- vapply(grp, sum, numeric(1L))
  cl_vars <- vapply(grp, var, numeric(1L))

  cl_vars <- .impute_singleton_var(cl_vars)

  tU <- sum(cl_tots)
  S2U1 <- sum(pp_psu * (cl_tots / pp_psu - tU)^2)
  vb <- S2U1 / tU^2

  ybarU <- mean(y)
  vw <- sum(Ni^2 * cl_vars / pp_psu) / tU^2
  S2U <- var(y)

  eps <- sqrt(.Machine$double.eps)
  rel_var <- if (abs(ybarU) < eps && S2U < eps) 0 else S2U / ybarU^2

  total_v <- vb + vw
  if (total_v < eps || rel_var < eps || !is.finite(rel_var)) {
    warning("outcome variance is approximately zero; delta set to 0 by convention",
            call. = FALSE)
    delta <- 0
    k <- 1
  } else {
    delta <- vb / total_v
    k <- total_v / rel_var
  }

  .new_svyplan_varcomp(
    varb    = vb,
    varw    = vw,
    delta   = delta,
    k       = k,
    rel_var = rel_var,
    stages  = 2L
  )
}

#' 3-stage PPS variance components
#' @keywords internal
#' @noRd
.varcomp_3stage_pps <- function(y, psu_id, ssu_id, pp) {
  unique_psu <- sort(unique(psu_id))
  M <- length(unique_psu)
  psu_idx <- match(psu_id, unique_psu)

  # Map pp to per-PSU
  if (length(pp) == length(y)) {
    pp_psu <- pp[match(unique_psu, psu_id)]
  } else if (length(pp) == M) {
    pp_psu <- pp
  } else {
    stop("'prob' must have length equal to number of observations or number of PSUs",
         call. = FALSE)
  }

  if (abs(sum(pp_psu) - 1) >= 1e-3) {
    stop("'prob' values must sum to 1", call. = FALSE)
  }

  # PSU totals
  grp_psu <- split(y, psu_idx)
  tUi <- vapply(grp_psu, sum, numeric(1L))
  tU <- sum(tUi)

  # Between-PSU variance
  S2U1pwr <- sum(pp_psu * (tUi / pp_psu - tU)^2)
  B <- S2U1pwr / tU^2

  # Nest SSU IDs within PSU (handles non-unique SSU IDs across PSUs)
  ssu_nested <- interaction(psu_id, ssu_id, drop = TRUE)
  unique_ssu <- unique(ssu_nested)
  first_psu_per_ssu <- psu_id[match(unique_ssu, ssu_nested)]
  psu_of_ssu_idx <- match(first_psu_per_ssu, unique_psu)
  Ni <- tabulate(psu_of_ssu_idx, nbins = M)

  # SSU totals and their variances within PSU
  ssu_idx <- match(ssu_nested, unique_ssu)
  grp_ssu <- split(y, ssu_idx)
  tij <- vapply(grp_ssu, sum, numeric(1L))
  grp_tij <- split(tij, psu_of_ssu_idx)
  S2U2i <- vapply(grp_tij, var, numeric(1L))

  S2U2i <- .impute_singleton_var(S2U2i)
  vw2 <- sum(Ni^2 * S2U2i / pp_psu) / tU^2

  # Element-level variance within PSU (for delta1 = B/(B+W))
  Qi <- tabulate(psu_idx, nbins = M)
  S2U3i <- vapply(grp_psu, var, numeric(1L))
  S2U3i <- .impute_singleton_var(S2U3i)
  W <- sum(Qi^2 * S2U3i / pp_psu) / tU^2

  # Element-level variance within SSU
  n_ssu <- length(unique_ssu)
  Qij <- tabulate(ssu_idx, nbins = n_ssu)
  S2U3ij <- vapply(grp_ssu, var, numeric(1L))
  S2U3ij <- .impute_singleton_var(S2U3ij)

  # Replicate pp and Ni to SSU level
  pp_ssu <- pp_psu[psu_of_ssu_idx]
  Ni_ssu <- Ni[psu_of_ssu_idx]

  vw3 <- sum(Ni_ssu * Qij^2 * S2U3ij / pp_ssu) / tU^2

  eps <- sqrt(.Machine$double.eps)
  ybarU <- mean(y)
  S2U <- var(y)
  V <- if (abs(ybarU) < eps && S2U < eps) 0 else S2U / ybarU^2

  warn <- FALSE
  bw <- B + W
  if (bw < eps || V < eps || !is.finite(V)) {
    warn <- TRUE
    delta1 <- 0
    k1 <- 1
  } else {
    delta1 <- B / bw
    k1 <- bw / V
  }
  ww <- vw2 + vw3
  if (ww < eps || V < eps || !is.finite(V)) {
    warn <- TRUE
    delta2 <- 0
    k2 <- 1
  } else {
    delta2 <- vw2 / ww
    k2 <- ww / V
  }
  if (warn) {
    warning("outcome variance is approximately zero; delta set to 0 by convention",
            call. = FALSE)
  }

  .new_svyplan_varcomp(
    varb    = B,
    varw    = c(varw_psu = vw2, varw_ssu = vw3),
    delta   = c(delta_psu = delta1, delta_ssu = delta2),
    k       = c(k_psu = k1, k_ssu = k2),
    rel_var = V,
    stages  = 3L
  )
}
