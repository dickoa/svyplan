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
#'   \item{`strata`}{Per-stratum component table when `strata` is
#'     supplied, otherwise `NULL` (see Details).}
#' }
#'
#' @details
#' The interface is determined by the class of `x`:
#'
#' - **Formula**: `varcomp(income ~ district/village, data = frame)`. The
#'   LHS is the analysis variable, RHS terms express nesting using
#'   `/` or `%in%` (see [formula]): `y ~ psu/ssu` or equivalently
#'   `y ~ ssu %in% psu` (SSUs nested within PSUs). With `/` the
#'   outermost stage comes first; with `%in%` the innermost comes
#'   first. Be careful not to reverse the order.
#' - **Numeric vector**: `varcomp(y, stage_id = list(cluster_ids))`.
#' - **survey.design**: `varcomp(design, ~y)`. Cluster structure and
#'   design weights are extracted from the design object. Requires the
#'   survey package. Weights are treated as inverse inclusion
#'   probabilities: cluster sizes and totals are estimated by summed
#'   weights, and the estimation variance of the weighted cluster
#'   totals is subtracted from the between-stage variance terms, so
#'   unequal-probability samples from a previous round give
#'   approximately design-unbiased components. Unit weights recover the
#'   frame formulas exactly. The weight *scale* matters: within each
#'   cluster the summed weights should estimate the cluster population
#'   size, since the implied sampling fraction drives the correction.
#'   The correction assumes noninformative (SRS-like) subsampling
#'   within clusters; informative within-cluster sampling remains
#'   approximate. For a PPS first stage, also pass the PSU selection
#'   probabilities via `prob`.
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
#' With `strata`, components are estimated separately within each
#' stratum and returned as a per-stratum table in `$strata` (also via
#' `as.data.frame()`); the pooled fields (`varb`, `delta`, ...) are not
#' filled. The table's columns (`sd`, `mean`, `delta_psu`, `k_psu`)
#' match the [n_alloc()] frame contract, so after adding stratum `N`
#' it feeds a stratified two-stage allocation directly. When `prob` is
#' combined with `strata`, supply one value per observation, summing
#' to 1 within each stratum.
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
#' # 2-stage SRS using formula (PSU = district)
#' set.seed(42)
#' frame2 <- data.frame(
#'   income = rnorm(200, 50000, 10000),
#'   district = rep(1:20, each = 10)
#' )
#' vc2 <- varcomp(income ~ district, data = frame2)
#' vc2
#'
#' # Feed into n_cluster
#' n_cluster(stage_cost = c(500, 50), delta = vc2, budget = 100000)
#'
#' # Per-stratum components for a stratified two-stage plan
#' set.seed(42)
#' frame_s <- data.frame(
#'   income = rnorm(400, 50000, 10000),
#'   district = rep(1:40, each = 10),
#'   region = rep(c("North", "South"), each = 200)
#' )
#' varcomp(income ~ district, data = frame_s, strata = ~region)
#'
#' # 3-stage SRS using formula: villages nested within districts
#' # "/" expresses nesting (outermost stage first, see ?formula)
#' set.seed(42)
#' frame3 <- data.frame(
#'   income = rnorm(400, 50000, 10000),
#'   district = rep(1:20, each = 20),
#'   village = rep(1:100, each = 4)
#' )
#' vc3 <- varcomp(income ~ district/village, data = frame3)
#' vc3
#'
#' # 3-stage PPS (explicit first-stage probabilities)
#' frame3$pp <- rep(1 / 20, 400)
#' vc3_pps <- varcomp(income ~ district/village, data = frame3, prob = ~pp)
#'
#' # Vector (list) interface
#' varcomp(frame3$income,
#'         stage_id = list(frame3$district, frame3$village))
#'
#' @export
varcomp <- function(x, ...) {
  UseMethod("varcomp")
}

#' @describeIn varcomp Method for formula interface.
#'
#' @param data A data frame (required for formula interface).
#' @param prob First-stage selection probabilities (PPS). A one-sided
#'   formula (e.g., `~pp`) when using the formula interface, or a numeric
#'   vector: either one value per observation (constant within each PSU)
#'   or one value per PSU. A one-per-PSU vector is matched by name when
#'   named; unnamed values are taken in sorted order of the unique PSU
#'   identifiers. Values must be strictly between 0 and 1 and sum to 1
#'   across PSUs (tolerance 1e-3). `NULL` (default) assumes SRS.
#' @param strata Optional stratification: a one-sided formula (formula
#'   and survey.design interfaces) or a vector (default interface).
#'   Components are then estimated per stratum; see Details.
#'
#' @export
varcomp.formula <- function(x, ..., data = NULL, prob = NULL, strata = NULL) {
  if (inherits(strata, "formula")) {
    strata_name <- all.vars(strata)
    if (length(strata_name) != 1L) {
      stop("'strata' formula must reference exactly one variable",
           call. = FALSE)
    }
    strata <- data[[strata_name]]
    if (is.null(strata)) {
      stop(sprintf("variable '%s' not found in 'data'", strata_name),
           call. = FALSE)
    }
  }
  .varcomp_formula(x, data = data, prob = prob, strata = strata)
}

#' @describeIn varcomp Default method for numeric vectors.
#'
#' @param stage_id A list of stage-ID vectors (required for vector interface).
#'   Length determines the number of stage boundaries (stages - 1).
#'
#' @export
varcomp.default <- function(x, ..., stage_id = NULL, prob = NULL,
                            strata = NULL) {
  if (is.numeric(x)) {
    .varcomp_vector(x, stage_id = stage_id, prob = prob, strata = strata)
  } else {
    stop("'x' must be a formula, numeric vector, or survey design object",
         call. = FALSE)
  }
}

#' @describeIn varcomp Method for survey design objects. Pass a one-sided
#'   formula (e.g., `~y`) to specify the outcome variable. Cluster
#'   structure and design weights are extracted from the design.
#'
#' @export
varcomp.survey.design <- function(x, ..., prob = NULL, strata = NULL) {
  if (inherits(strata, "formula")) {
    strata_name <- all.vars(strata)
    if (length(strata_name) != 1L) {
      stop("'strata' formula must reference exactly one variable",
           call. = FALSE)
    }
    strata <- x$variables[[strata_name]]
    if (is.null(strata)) {
      stop(sprintf("variable '%s' not found in design variables", strata_name),
           call. = FALSE)
    }
  }
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

  w <- as.numeric(stats::weights(x))
  if (anyNA(w) || any(!is.finite(w)) || any(w <= 0)) {
    stop("design weights must be positive and finite", call. = FALSE)
  }
  rng <- range(w)
  if ((rng[2L] - rng[1L]) / rng[2L] <= 1e-8 && abs(rng[1L] - 1) <= 1e-8) {
    w <- NULL
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
  .varcomp_dispatch(y, stage_id, prob, w, strata)
}

#' Parse formula interface and dispatch
#' @keywords internal
#' @noRd
.varcomp_formula <- function(formula, data, prob, strata = NULL) {
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

  if (length(stage_names) > 1L) {
    tl <- attr(terms(formula), "term.labels")
    n_stg <- length(stage_names)
    slash_ok <- length(tl) == n_stg && !grepl(":", tl[1L], fixed = TRUE)
    if (slash_ok) {
      for (i in seq_along(tl)[-1L]) {
        if (!startsWith(tl[i], paste0(tl[i - 1L], ":"))) {
          slash_ok <- FALSE
          break
        }
      }
    }
    in_ok <- !slash_ok && length(tl) == 1L &&
      length(strsplit(tl, ":")[[1L]]) == n_stg
    if (slash_ok) {
      stage_names <- tl[1L]
      for (i in seq_along(tl)[-1L]) {
        stage_names <- c(stage_names,
                         sub(paste0(tl[i - 1L], ":"), "", tl[i], fixed = TRUE))
      }
    } else if (in_ok) {
      stage_names <- rev(strsplit(tl, ":")[[1L]])
    } else {
      stop(
        "multi-stage formula must express nesting ",
        "(e.g., y ~ psu/ssu or y ~ ssu %in% psu); see ?formula",
        call. = FALSE
      )
    }
  }

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

  .varcomp_dispatch(y, stage_id, pp, strata = strata)
}

#' Vector interface
#' @keywords internal
#' @noRd
.varcomp_vector <- function(x, stage_id, prob, strata = NULL) {
  if (is.null(stage_id) || !is.list(stage_id)) {
    stop("'stage_id' must be a list of ID vectors", call. = FALSE)
  }
  .varcomp_dispatch(x, stage_id, prob, strata = strata)
}

#' Dispatch to correct variance component estimator
#' @keywords internal
#' @noRd
.varcomp_dispatch <- function(y, stage_id, prob, w = NULL, strata = NULL) {
  if (!is.numeric(y) || length(y) == 0L) {
    stop("'y' must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("outcome vector must not contain NA values", call. = FALSE)
  }
  if (any(!is.finite(y))) {
    stop("outcome vector must contain only finite values", call. = FALSE)
  }
  if (length(stage_id) == 0L) {
    stop("'stage_id' must not be empty", call. = FALSE)
  }
  for (i in seq_along(stage_id)) {
    if (length(stage_id[[i]]) != length(y)) {
      stop(sprintf("'stage_id[[%d]]' must have length %d (same as outcome vector)",
                   i, length(y)), call. = FALSE)
    }
    if (anyNA(stage_id[[i]])) {
      stop(sprintf("'stage_id[[%d]]' must not contain NA values", i),
           call. = FALSE)
    }
  }

  n_boundaries <- length(stage_id)
  stages <- n_boundaries + 1L

  if (stages > 3L) {
    stop(
      "4+ stage variance components are not supported; estimate the top three stages and fold deeper stages (e.g. persons within households) into the 'deff' passed to n_prop() or n_mean()",
      call. = FALSE
    )
  }

  if (!is.null(strata)) {
    return(.varcomp_strata(y, stage_id, prob, w, strata))
  }

  has_prob <- !is.null(prob)

  if (stages == 2L && !has_prob) {
    .varcomp_2stage_srs(y, stage_id[[1L]], w)
  } else if (stages == 2L && has_prob) {
    .varcomp_2stage_pps(y, stage_id[[1L]], prob, w)
  } else if (stages == 3L && !has_prob) {
    M <- length(unique(stage_id[[1L]]))
    .varcomp_3stage_pps(y, stage_id[[1L]], stage_id[[2L]], rep(1 / M, M), w)
  } else {
    .varcomp_3stage_pps(y, stage_id[[1L]], stage_id[[2L]], prob, w)
  }
}

#' Per-stratum variance components
#'
#' Splits the data by stratum and runs the requested estimator within
#' each. Column names of the result (`sd`, `mean`, `delta_psu`, `k_psu`,
#' ...) deliberately match the n_alloc() frame contract so the table can
#' be merged into an allocation frame directly.
#' @keywords internal
#' @noRd
.varcomp_strata <- function(y, stage_id, prob, w, strata) {
  if (length(strata) != length(y)) {
    stop("'strata' must have the same length as the outcome", call. = FALSE)
  }
  if (anyNA(strata)) {
    stop("'strata' must not contain NA values", call. = FALSE)
  }
  if (!is.null(prob) && length(prob) != length(y)) {
    stop(
      "with 'strata', 'prob' must have one value per observation (summing to 1 within each stratum)",
      call. = FALSE
    )
  }

  strata <- as.character(strata)
  lev <- unique(strata)
  rows <- lapply(lev, function(s) {
    idx <- which(strata == s)
    vc <- tryCatch(
      .varcomp_dispatch(
        y[idx],
        lapply(stage_id, function(id) id[idx]),
        if (!is.null(prob)) prob[idx],
        if (!is.null(w)) w[idx]
      ),
      error = function(e) {
        stop(sprintf("stratum '%s': %s", s, conditionMessage(e)),
             call. = FALSE)
      }
    )
    row <- if (is.null(w)) {
      data.frame(stratum = s, sd = sd(y[idx]), mean = mean(y[idx]))
    } else {
      data.frame(
        stratum = s,
        sd = sqrt(.varcomp_wvar(y[idx], w[idx])),
        mean = sum(w[idx] * y[idx]) / sum(w[idx])
      )
    }
    if (vc$stages == 2L) {
      row$delta_psu <- vc$delta
      row$k_psu <- vc$k
      row$varb <- vc$varb
      row$varw <- vc$varw
    } else {
      row$delta_psu <- vc$delta[["delta_psu"]]
      row$delta_ssu <- vc$delta[["delta_ssu"]]
      row$k_psu <- vc$k[["k_psu"]]
      row$k_ssu <- vc$k[["k_ssu"]]
      row$varb <- vc$varb
      row$varw_psu <- vc$varw[["varw_psu"]]
      row$varw_ssu <- vc$varw[["varw_ssu"]]
    }
    row$rel_var <- vc$rel_var
    row
  })

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  .new_svyplan_varcomp(
    varb    = NULL,
    varw    = NULL,
    delta   = NULL,
    k       = NULL,
    rel_var = NULL,
    stages  = length(stage_id) + 1L,
    strata  = tab
  )
}

#' Map and validate PPS probabilities to per-PSU values
#'
#' Accepts one value per observation (must be constant within PSU) or one
#' per PSU. A one-per-PSU vector is matched by name when named; unnamed
#' values are taken in sorted order of the unique PSU identifiers.
#' @keywords internal
#' @noRd
.varcomp_map_pp <- function(pp, y, psu_id, unique_psu) {
  M <- length(unique_psu)
  if (!is.numeric(pp) || anyNA(pp) || any(!is.finite(pp)) ||
      any(pp <= 0) || any(pp >= 1)) {
    stop("'prob' values must be strictly between 0 and 1", call. = FALSE)
  }
  if (length(pp) == length(y)) {
    grp_pp <- split(pp, match(psu_id, unique_psu))
    spread <- vapply(grp_pp, function(v) max(v) - min(v), numeric(1L))
    if (any(spread > 1e-12)) {
      stop("'prob' must be constant within each PSU when given per observation",
           call. = FALSE)
    }
    pp_psu <- pp[match(unique_psu, psu_id)]
  } else if (length(pp) == M) {
    if (!is.null(names(pp))) {
      if (anyDuplicated(names(pp))) {
        stop("names of 'prob' must be unique", call. = FALSE)
      }
      m <- match(as.character(unique_psu), names(pp))
      if (anyNA(m)) {
        stop("names of 'prob' must match the PSU identifiers", call. = FALSE)
      }
      pp_psu <- as.numeric(pp[m])
    } else {
      pp_psu <- pp
    }
  } else {
    stop("'prob' must have length equal to number of observations or number of PSUs",
         call. = FALSE)
  }
  if (abs(sum(pp_psu) - 1) >= 1e-3) {
    stop("'prob' values must sum to 1", call. = FALSE)
  }
  pp_psu
}

#' Require at least two PSUs for a between-PSU variance
#' @keywords internal
#' @noRd
.check_min_psu <- function(M) {
  if (M < 2L) {
    stop(
      "at least two PSUs are required to estimate between-PSU variance; collapse single-PSU strata with a neighbour",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Group-safe weighted variance (NA for singletons)
#' @keywords internal
#' @noRd
.varcomp_wvar <- function(x, w) {
  if (length(x) < 2L) {
    return(NA_real_)
  }
  .wtdvar(x, w)
}

#' Estimation variance of a weighted group total under SRSWOR within
#' the group; zero when the group's weights are all 1 (Nhat = n)
#' @keywords internal
#' @noRd
.varcomp_total_var <- function(n_i, Nhat_i, S2_i) {
  Nhat_i^2 * pmax(1 - n_i / Nhat_i, 0) * S2_i / n_i
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
.varcomp_2stage_srs <- function(y, psu_id, w = NULL) {
  uid <- unique(psu_id)
  M <- length(uid)
  .check_min_psu(M)
  idx <- match(psu_id, uid)
  grp <- split(y, idx)

  if (is.null(w)) {
    Ni <- tabulate(idx, nbins = M)
    ti <- vapply(grp, sum, numeric(1L))
    S2Ui <- .impute_singleton_var(vapply(grp, var, numeric(1L)))
    S2U1 <- var(ti)
    ybarU <- mean(y)
    S2U <- var(y)
  } else {
    grp_w <- split(w, idx)
    Ni <- vapply(grp_w, sum, numeric(1L))
    ti <- vapply(seq_len(M), function(i) sum(grp_w[[i]] * grp[[i]]),
                 numeric(1L))
    S2Ui <- .impute_singleton_var(
      vapply(seq_len(M), function(i) .varcomp_wvar(grp[[i]], grp_w[[i]]),
             numeric(1L))
    )
    ni <- tabulate(idx, nbins = M)
    Vti <- .varcomp_total_var(ni, Ni, S2Ui)
    S2U1 <- max(var(ti) - mean(Vti), 0)
    ybarU <- sum(w * y) / sum(w)
    S2U <- .varcomp_wvar(y, w)
  }

  tbarU <- mean(ti)
  tU <- M * tbarU
  vb <- S2U1 / tbarU^2
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
.varcomp_2stage_pps <- function(y, psu_id, pp, w = NULL) {
  unique_psu <- sort(unique(psu_id))
  M <- length(unique_psu)
  .check_min_psu(M)
  idx <- match(psu_id, unique_psu)

  pp_psu <- .varcomp_map_pp(pp, y, psu_id, unique_psu)

  grp <- split(y, idx)

  if (is.null(w)) {
    Ni <- tabulate(idx, nbins = M)
    cl_tots <- vapply(grp, sum, numeric(1L))
    cl_vars <- .impute_singleton_var(vapply(grp, var, numeric(1L)))
    tU <- sum(cl_tots)
    S2U1 <- sum(pp_psu * (cl_tots / pp_psu - tU)^2)
    ybarU <- mean(y)
    S2U <- var(y)
  } else {
    grp_w <- split(w, idx)
    Ni <- vapply(grp_w, sum, numeric(1L))
    cl_tots <- vapply(seq_len(M), function(i) sum(grp_w[[i]] * grp[[i]]),
                      numeric(1L))
    cl_vars <- .impute_singleton_var(
      vapply(seq_len(M), function(i) .varcomp_wvar(grp[[i]], grp_w[[i]]),
             numeric(1L))
    )
    tU <- sum(cl_tots)
    ni <- tabulate(idx, nbins = M)
    Vti <- .varcomp_total_var(ni, Ni, cl_vars)
    S2U1 <- max(
      sum(pp_psu * (cl_tots / pp_psu - tU)^2) -
        sum(Vti * (1 - pp_psu) / pp_psu),
      0
    )
    ybarU <- sum(w * y) / sum(w)
    S2U <- .varcomp_wvar(y, w)
  }

  vb <- S2U1 / tU^2
  vw <- sum(Ni^2 * cl_vars / pp_psu) / tU^2

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
.varcomp_3stage_pps <- function(y, psu_id, ssu_id, pp, w = NULL) {
  unique_psu <- sort(unique(psu_id))
  M <- length(unique_psu)
  .check_min_psu(M)
  psu_idx <- match(psu_id, unique_psu)

  pp_psu <- .varcomp_map_pp(pp, y, psu_id, unique_psu)

  weighted <- !is.null(w)

  # PSU totals
  grp_psu <- split(y, psu_idx)
  if (weighted) {
    grp_psu_w <- split(w, psu_idx)
    tUi <- vapply(seq_len(M), function(i) sum(grp_psu_w[[i]] * grp_psu[[i]]),
                  numeric(1L))
  } else {
    tUi <- vapply(grp_psu, sum, numeric(1L))
  }
  tU <- sum(tUi)

  # Between-PSU variance (bias-corrected below for weighted samples)
  S2U1pwr <- sum(pp_psu * (tUi / pp_psu - tU)^2)

  # Nest SSU IDs within PSU (handles non-unique SSU IDs across PSUs)
  ssu_nested <- interaction(psu_id, ssu_id, drop = TRUE)
  unique_ssu <- unique(ssu_nested)
  first_psu_per_ssu <- psu_id[match(unique_ssu, ssu_nested)]
  psu_of_ssu_idx <- match(first_psu_per_ssu, unique_psu)
  Ni <- tabulate(psu_of_ssu_idx, nbins = M)

  # SSU totals and their variances within PSU
  ssu_idx <- match(ssu_nested, unique_ssu)
  n_ssu <- length(unique_ssu)
  grp_ssu <- split(y, ssu_idx)
  if (weighted) {
    grp_ssu_w <- split(w, ssu_idx)
    tij <- vapply(seq_len(n_ssu),
                  function(j) sum(grp_ssu_w[[j]] * grp_ssu[[j]]),
                  numeric(1L))
  } else {
    tij <- vapply(grp_ssu, sum, numeric(1L))
  }
  grp_tij <- split(tij, psu_of_ssu_idx)
  S2U2i <- .impute_singleton_var(vapply(grp_tij, var, numeric(1L)))

  # Element-level variance within PSU (for delta1 = B/(B+W))
  if (weighted) {
    Qi <- vapply(grp_psu_w, sum, numeric(1L))
    S2U3i <- vapply(seq_len(M), function(i) {
      .varcomp_wvar(grp_psu[[i]], grp_psu_w[[i]])
    }, numeric(1L))
  } else {
    Qi <- tabulate(psu_idx, nbins = M)
    S2U3i <- vapply(grp_psu, var, numeric(1L))
  }
  S2U3i <- .impute_singleton_var(S2U3i)
  W <- sum(Qi^2 * S2U3i / pp_psu) / tU^2

  # Element-level variance within SSU
  if (weighted) {
    Qij <- vapply(grp_ssu_w, sum, numeric(1L))
    S2U3ij <- vapply(seq_len(n_ssu), function(j) {
      .varcomp_wvar(grp_ssu[[j]], grp_ssu_w[[j]])
    }, numeric(1L))
  } else {
    Qij <- tabulate(ssu_idx, nbins = n_ssu)
    S2U3ij <- vapply(grp_ssu, var, numeric(1L))
  }
  S2U3ij <- .impute_singleton_var(S2U3ij)

  if (weighted) {
    # Estimation variance of SSU and PSU totals (zero for w = 1)
    qij <- tabulate(ssu_idx, nbins = n_ssu)
    Vtij <- .varcomp_total_var(qij, Qij, S2U3ij)
    Vtij_psu <- split(Vtij, psu_of_ssu_idx)
    S2U2i <- pmax(S2U2i - vapply(Vtij_psu, mean, numeric(1L)), 0)
    Vti <- vapply(Vtij_psu, sum, numeric(1L))
    S2U1pwr <- max(S2U1pwr - sum(Vti * (1 - pp_psu) / pp_psu), 0)
  }

  B <- S2U1pwr / tU^2
  vw2 <- sum(Ni^2 * S2U2i / pp_psu) / tU^2

  # Replicate pp and Ni to SSU level
  pp_ssu <- pp_psu[psu_of_ssu_idx]
  Ni_ssu <- Ni[psu_of_ssu_idx]

  vw3 <- sum(Ni_ssu * Qij^2 * S2U3ij / pp_ssu) / tU^2

  eps <- sqrt(.Machine$double.eps)
  if (weighted) {
    ybarU <- sum(w * y) / sum(w)
    S2U <- .varcomp_wvar(y, w)
  } else {
    ybarU <- mean(y)
    S2U <- var(y)
  }
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
