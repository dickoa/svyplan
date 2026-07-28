#' Estimate Variance Components
#'
#' Estimate between- and within-stage variance components using nested
#' ANOVA decomposition. Supports SRS and PPS first-stage designs.
#'
#' @param x A formula, numeric vector, or survey design object
#'   (see Details).
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#'
#' @return A `svyplan_varcomp` object with components:
#' \describe{
#'   \item{`varb`}{Between-PSU variance (scalar).}
#'   \item{`varw`}{Within-PSU variance. Scalar for 2-stage, length-2
#'     vector (`varw_psu`, `varw_ssu`) for 3-stage.}
#'   \item{`icc`}{Measure of homogeneity. Length 1 for 2-stage,
#'     length 2 (`icc_psu`, `icc_ssu`) for 3-stage.}
#'   \item{`var_ratio`}{Ratio parameter(s), same length as `icc`
#'     (`var_ratio_psu`, `var_ratio_ssu` for 3-stage).}
#'   \item{`unit_relvar`}{Unit relvariance (scalar).}
#'   \item{`stages`}{Number of stages (2 or 3).}
#'   \item{`strata`}{Per-stratum component table when `strata` is
#'     supplied, otherwise `NULL` (see Details).}
#' }
#' The result can be exported with [as.data.frame()]. Unstratified results
#' become a one-row two- or three-stage component table. Stratified results
#' return the table stored in `$strata`.
#'
#' @details
#' The interface is determined by the class of `x`:
#'
#' - **Formula**: `varcomp(income ~ district/village, data = frame)`. The
#'   LHS is the analysis variable, RHS terms express nesting using
#'   `/` or `%in%` (see [formula]): `y ~ psu/ssu` or equivalently
#'   `y ~ ssu %in% psu` (SSUs nested within PSUs). With `/` the
#'   outermost stage comes first, whereas with `%in%` the innermost comes
#'   first. Be careful not to reverse the order.
#' - **Numeric vector**: `varcomp(y, stage_id = list(cluster_ids))`.
#' - **survey.design**: `varcomp(design, ~y)`. Cluster structure and
#'   design weights are extracted from the design object. Requires the
#'   survey package. For a PPS first stage, also pass the PSU selection
#'   probabilities via `prob`.
#'
#' The formula and vector interfaces compute frame components: `data`
#' is treated as a complete population. To estimate components from a
#' sample instead, supply `weights` (formula and vector interfaces) or
#' use the survey.design method. Weights are treated as inverse
#' inclusion probabilities: cluster sizes and totals are estimated by
#' summed weights, and the estimation variance of the weighted cluster
#' totals is subtracted from the between-stage variance terms, so
#' unequal-probability samples from a previous round give approximately
#' design-unbiased components. Unit weights recover the frame formulas
#' exactly. The weight *scale* matters: within each cluster the summed
#' weights should estimate the cluster population size, since the
#' implied sampling fraction drives the correction. For a multi-stage
#' sample this means the *within-cluster* weights (the product of the
#' stage-2 and later weights), not the full design weight, whose
#' stage-1 factor would overstate every cluster size. The correction
#' assumes noninformative (SRS-like) subsampling within clusters.
#' Informative within-cluster sampling remains approximate. In the
#' 3-stage weighted case the between-PSU correction removes
#' element-stage estimation noise but not SSU-stage subsampling noise,
#' so when few SSUs are sampled per PSU the between-PSU component (and
#' `icc_psu`) is conservative: biased upward, never downward.
#'
#' When `prob` is `NULL`, SRS first-stage is assumed. When provided, PPS
#' variance estimation is used. `prob` must sum to 1 across the PSUs
#' present in the data. On a complete frame these are the one-draw
#' selection probabilities themselves. In a sample containing only some
#' of the PSUs, renormalize over the sampled PSUs:
#' `prob = pp[sampled] / sum(pp[sampled])`. For a fixed-size PPS design
#' the stage-1 inclusion probabilities are proportional to the one-draw
#' probabilities, so shares derived from inverse stage-1 weights,
#' `pi1 / sum(pi1)`, are identical and can be used when the frame is
#' not at hand. Certainty (take-all) PSUs have no place in the PPS
#' path: they contribute no between-PSU variance and should be treated
#' as separate strata, with `prob` renormalized over the remaining
#' PSUs.
#'
#' The returned `icc` is the design-based measure of homogeneity
#' \eqn{\delta = V_b / (V_b + V_w)}{icc = Vb / (Vb + Vw)}, written
#' \eqn{\delta} in Valliant, Dever, and Kreuter (2018, Ch. 9).
#' Unlike the traditional ANOVA
#' intraclass correlation coefficient, `icc` is constrained to \eqn{[0, 1]}
#' and should not be compared directly to mixed-model ICCs (e.g. from lme4)
#' which can be negative.
#'
#' ## Component estimands
#'
#' For two-stage SRS, let \eqn{M} be the number of PSUs, \eqn{N_i} the
#' ultimate-unit count in PSU \eqn{i}, \eqn{t_i} its outcome total, and
#' \eqn{S_i^2} its within-PSU variance. With
#' \eqn{t_U=M\bar t}, the complete-frame components returned are
#' \deqn{V_b=s^2(t_i)/\bar t^2, \qquad
#'       V_w=M\sum_iN_i^2S_i^2/t_U^2.}
#' For PPS with one-draw probabilities \eqn{p_i}, they are
#' \deqn{V_b=\sum_i p_i(t_i/p_i-t_U)^2/t_U^2, \qquad
#'       V_w=\sum_iN_i^2S_i^2/(p_it_U^2).}
#' In both cases `unit_relvar` is the ultimate-unit variance divided by the
#' squared ultimate-unit mean, `icc = varb / (varb + varw)`, and
#' `var_ratio = (varb + varw) / unit_relvar`.
#'
#' For three stages, `varb` uses the PPS between-PSU expression above.
#' `varw_psu` applies the same PPS scaling to the variance of SSU totals
#' within each PSU, and `varw_ssu` aggregates the within-SSU ultimate-unit
#' variances. The returned PSU homogeneity compares the between-PSU component
#' with the element-level variance aggregated within PSUs; the SSU
#' homogeneity is `varw_psu / (varw_psu + varw_ssu)`. This distinction is why
#' the first `icc` is not generally `varb / sum(c(varb, varw))` for a
#' three-stage result.
#'
#' The two `var_ratio` values are estimated from their own stage decompositions
#' rather than by imposing the identity `var_ratio_ssu = var_ratio_psu * (1 - icc_psu)`
#' that the planning formula uses (see [design_effect()]). On small clusters
#' the two can differ by several percent. Passing a `svyplan_varcomp` to a
#' planning function uses the estimated pair as given; omitting `var_ratio_ssu` there
#' applies the identity instead.
#'
#' When weights are supplied, \eqn{N_i}, totals, means, and variances in these
#' expressions are their weighted estimates. The code subtracts the stated
#' SRS-without-replacement estimation variance of estimated cluster totals and
#' truncates negative corrected components to zero. These corrections are
#' planning approximations, not exact variance estimators for arbitrary
#' informative multistage samples.
#'
#' Clusters containing a single observation have undefined within-cluster
#' variance. In this case, the within-cluster variance is imputed as the
#' mean variance of the remaining clusters.
#'
#' With `strata`, components are estimated separately within each
#' stratum and returned as a per-stratum table in `$strata` (also via
#' `as.data.frame()`). The pooled fields (`varb`, `icc`, ...) are not
#' filled. The table's columns (`sd`, `mean`, `icc_psu`, `var_ratio_psu`)
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
#' @seealso [n_cluster()] which accepts a `svyplan_varcomp` as `icc`.
#'
#' @examples
#' # 2-stage SRS using formula (PSU = district)
#' set.seed(314)
#' frame2 <- data.frame(
#'   income = rnorm(200, 50000, 10000),
#'   district = rep(1:20, each = 10)
#' )
#' vc2 <- varcomp(income ~ district, data = frame2)
#' vc2
#' as.data.frame(vc2)
#'
#' # Feed into n_cluster
#' n_cluster(stage_cost = c(500, 50), icc = vc2, budget = 100000)
#'
#' # Estimate components from a two-stage sample: within-cluster
#' # weights (summing to each cluster's population size) and, for a
#' # PPS first stage, one-draw probabilities renormalized over the
#' # sampled clusters
#' sampled <- unlist(lapply(split(seq_len(200), frame2$district)[1:8],
#'                          sample, size = 4))
#' samp <- frame2[sampled, ]
#' samp$w <- 10 / 4
#' pp <- rep(1 / 20, 8)
#' varcomp(income ~ district, data = samp, weights = ~w,
#'         prob = pp / sum(pp))
#'
#' # Per-stratum components for a stratified two-stage plan
#' set.seed(2718)
#' frame_s <- data.frame(
#'   income = rnorm(400, 50000, 10000),
#'   district = rep(1:40, each = 10),
#'   region = rep(c("North", "South"), each = 200)
#' )
#' varcomp(income ~ district, data = frame_s, strata = ~region)
#'
#' # 3-stage SRS using formula: villages nested within districts
#' # "/" expresses nesting (outermost stage first, see ?formula)
#' set.seed(1618)
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
#'   named. Unnamed values are taken in sorted order of the unique PSU
#'   identifiers. Values must be strictly between 0 and 1 and sum to 1
#'   across PSUs (tolerance 1e-3). `NULL` (default) assumes SRS.
#' @param strata Optional stratification: a one-sided formula (formula
#'   and survey.design interfaces) or a vector (default interface).
#'   Components are then estimated per stratum. See Details.
#' @param weights Optional sampling weights for estimating components
#'   from a sample rather than a frame: a one-sided formula (e.g.,
#'   `~w`) when using the formula interface, or a numeric vector with
#'   one positive value per observation. Within each cluster the
#'   summed weights must estimate the cluster population size. For a
#'   multi-stage sample use the within-cluster weights, not the full
#'   design weight (see Details). `NULL` (default) computes frame
#'   components.
#'
#' @export
varcomp.formula <- function(x, ..., data = NULL, prob = NULL, strata = NULL,
                            weights = NULL) {
  .check_unused_dots(...)
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
  if (inherits(weights, "formula")) {
    w_name <- all.vars(weights)
    if (length(w_name) != 1L) {
      stop("'weights' formula must reference exactly one variable",
           call. = FALSE)
    }
    weights <- data[[w_name]]
    if (is.null(weights)) {
      stop(sprintf("variable '%s' not found in 'data'", w_name),
           call. = FALSE)
    }
  }
  if (is.null(weights) && inherits(data, "tbl_sample")) {
    warning(
      "'data' is a samplyr sample but the formula interface treats it as a complete population frame. Pass within-cluster 'weights' (and 'prob' for a PPS first stage) for design-based components",
      call. = FALSE
    )
  }
  .varcomp_formula(x, data = data, prob = prob, strata = strata, w = weights)
}

#' @describeIn varcomp Default method for numeric vectors.
#'
#' @param stage_id A list of stage-ID vectors (required for vector interface).
#'   Length determines the number of stage boundaries (stages - 1).
#'
#' @export
varcomp.default <- function(x, ..., stage_id = NULL, prob = NULL,
                            strata = NULL, weights = NULL) {
  .check_unused_dots(...)
  if (is.numeric(x)) {
    .varcomp_vector(x, stage_id = stage_id, prob = prob, strata = strata,
                    w = weights)
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
  dots <- list(...)
  is_formula <- vapply(dots, inherits, logical(1L), "formula")
  formula_idx <- which(is_formula)
  if (length(formula_idx) == 0L) {
    stop("a one-sided formula specifying the outcome is required (e.g., ~y)",
         call. = FALSE)
  }
  if (length(formula_idx) > 1L) {
    stop("supply exactly one outcome formula (e.g., ~y)", call. = FALSE)
  }
  .stop_unused_dots(names(dots), setdiff(seq_along(dots), formula_idx))

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
  formula <- dots[[formula_idx]]

  w <- .varcomp_check_w(as.numeric(stats::weights(x)), "design weights")

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
    stop("design has no clusters (ids = ~1). varcomp requires a clustered design",
         call. = FALSE)
  }

  stage_id <- lapply(seq_len(n_stages), function(j) cl[[j]])
  .varcomp_dispatch(y, stage_id, prob, w, strata)
}

#' Validate sampling weights. Unit weights collapse to NULL so the
#' exact frame formulas apply
#' @keywords internal
#' @noRd
.varcomp_check_w <- function(w, what = "'weights'") {
  if (is.null(w)) {
    return(NULL)
  }
  if (!is.numeric(w) || length(w) == 0L || anyNA(w) ||
      any(!is.finite(w)) || any(w <= 0)) {
    stop(sprintf("%s must be positive and finite", what), call. = FALSE)
  }
  rng <- range(w)
  if ((rng[2L] - rng[1L]) / rng[2L] <= 1e-8 && abs(rng[1L] - 1) <= 1e-8) {
    return(NULL)
  }
  w
}

#' Parse formula interface and dispatch
#' @keywords internal
#' @noRd
.varcomp_formula <- function(formula, data, prob, strata = NULL, w = NULL) {
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
        "(e.g., y ~ psu/ssu or y ~ ssu %in% psu). See ?formula",
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

  .varcomp_dispatch(y, stage_id, pp, .varcomp_check_w(w), strata)
}

#' Vector interface
#' @keywords internal
#' @noRd
.varcomp_vector <- function(x, stage_id, prob, strata = NULL, w = NULL) {
  if (is.null(stage_id) || !is.list(stage_id)) {
    stop("'stage_id' must be a list of ID vectors", call. = FALSE)
  }
  .varcomp_dispatch(x, stage_id, prob, .varcomp_check_w(w), strata)
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
  if (!is.null(w) && length(w) != length(y)) {
    stop("'weights' must have the same length as the outcome vector",
         call. = FALSE)
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
      "4+ stage variance components are not supported. Estimate the top three stages and fold deeper stages (e.g. persons within households) into the 'deff' passed to n_prop() or n_mean()",
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
#' each. Column names of the result (`sd`, `mean`, `icc_psu`, `var_ratio_psu`,
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
        sd = sqrt(.vc_group_var(y[idx], w[idx])),
        mean = sum(w[idx] * y[idx]) / sum(w[idx])
      )
    }
    if (vc$stages == 2L) {
      row$icc_psu <- vc$icc
      row$var_ratio_psu <- vc$var_ratio
      row$varb <- vc$varb
      row$varw <- vc$varw
    } else {
      row$icc_psu <- vc$icc[["icc_psu"]]
      row$icc_ssu <- vc$icc[["icc_ssu"]]
      row$var_ratio_psu <- vc$var_ratio[["var_ratio_psu"]]
      row$var_ratio_ssu <- vc$var_ratio[["var_ratio_ssu"]]
      row$varb <- vc$varb
      row$varw_psu <- vc$varw[["varw_psu"]]
      row$varw_ssu <- vc$varw[["varw_ssu"]]
    }
    row$unit_relvar <- vc$unit_relvar
    row
  })

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  .new_svyplan_varcomp(
    varb    = NULL,
    varw    = NULL,
    icc   = NULL,
    var_ratio       = NULL,
    unit_relvar = NULL,
    stages  = length(stage_id) + 1L,
    strata  = tab
  )
}

#' Map and validate PPS probabilities to per-PSU values
#'
#' Accepts one value per observation (must be constant within PSU) or one
#' per PSU. A one-per-PSU vector is matched by name when named. Unnamed
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
.check_min_psu <- function(psu_count) {
  if (psu_count < 2L) {
    stop(
      "at least two PSUs are required to estimate between-PSU variance. Collapse single-PSU strata with a neighbour",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Tolerance below which a variance counts as zero
#' @keywords internal
#' @noRd
.vc_eps <- function() {
  sqrt(.Machine$double.eps)
}

#' Split observations into groups, preserving the caller's level order
#'
#' Returns the per-observation group index and the number of groups, so
#' every downstream summary can be built with a single split.
#' @keywords internal
#' @noRd
.vc_group <- function(id, levels) {
  list(index = match(id, levels), count = length(levels))
}

#' Per-group size, total, and unit variance
#'
#' Size is the count of observations without weights and the summed weight
#' with them, so it estimates the group's population count either way.
#' Totals and variances follow the same substitution. Groups with a single
#' observation have no variance estimate and are filled in by the caller.
#' @keywords internal
#' @noRd
.vc_summarise <- function(y, w, group) {
  parts <- split(y, group$index)
  observed <- tabulate(group$index, nbins = group$count)
  if (is.null(w)) {
    list(
      size = as.numeric(observed),
      total = vapply(parts, sum, numeric(1L), USE.NAMES = FALSE),
      var = vapply(parts, var, numeric(1L), USE.NAMES = FALSE),
      observed = observed
    )
  } else {
    weight_parts <- split(w, group$index)
    list(
      size = vapply(weight_parts, sum, numeric(1L), USE.NAMES = FALSE),
      total = vapply(
        seq_len(group$count),
        function(g) sum(weight_parts[[g]] * parts[[g]]),
        numeric(1L)
      ),
      var = vapply(
        seq_len(group$count),
        function(g) .vc_group_var(parts[[g]], weight_parts[[g]]),
        numeric(1L)
      ),
      observed = observed
    )
  }
}

#' Weighted variance of one group, NA when it holds a single observation
#' @keywords internal
#' @noRd
.vc_group_var <- function(y, w) {
  if (length(y) < 2L) {
    return(NA_real_)
  }
  .wtdvar(y, w)
}

#' Fill in the variance of groups holding a single observation
#'
#' A singleton group carries no information about within-group spread, so
#' it borrows the average of the groups that do.
#' @keywords internal
#' @noRd
.vc_fill_singletons <- function(group_var) {
  alone <- is.na(group_var)
  if (all(alone)) {
    warning("all clusters are singletons. Within-cluster variance set to 0",
            call. = FALSE)
    group_var[] <- 0
    return(group_var)
  }
  if (any(alone)) {
    group_var[alone] <- mean(group_var[!alone])
  }
  group_var
}

#' Sampling variance of an estimated group total under SRSWOR in the group
#'
#' The usual N^2 (1 - f) S^2 / n, with the group's population size
#' estimated by its summed weights. Zero when the weights are all 1,
#' because the group is then fully enumerated.
#' @keywords internal
#' @noRd
.vc_total_var <- function(observed, size, group_var) {
  size^2 * pmax(1 - observed / size, 0) * group_var / observed
}

#' One variance component, as a relvariance and as a bare numerator
#'
#' Every component here is a variance divided by the square of a total or a
#' mean, and every one of those denominators is the same multiple of the
#' population mean. Carrying the numerator on that common scale alongside
#' the relvariance is what lets `icc` and `var_ratio` be formed when the
#' mean is zero and each relvariance on its own is infinite.
#' @keywords internal
#' @noRd
.vc_component <- function(num, den) {
  list(
    num = num,
    rel = if (den > 0) num / den else if (num > 0) Inf else 0
  )
}

#' Relvariance of the between-PSU term, SRS first stage
#'
#' The frame variance of the PSU totals relative to their mean.
#' @keywords internal
#' @noRd
.vc_between_srs <- function(psu_total, adjustment = 0) {
  spread <- max(var(psu_total) - adjustment, 0)
  n <- length(psu_total)
  .vc_component(spread * n^2, sum(psu_total)^2)
}

#' Relvariance of the between-PSU term, PPS first stage
#'
#' The with-replacement (Hansen-Hurwitz) dispersion of the inflated PSU
#' totals about the grand total, relative to its square.
#' @keywords internal
#' @noRd
.vc_between_pps <- function(psu_total, draw_prob, adjustment = 0) {
  grand_total <- sum(psu_total)
  spread <- sum(draw_prob * (psu_total / draw_prob - grand_total)^2)
  .vc_component(max(spread - adjustment, 0), grand_total^2)
}

#' Relvariance contributed by the spread inside each group
#'
#' Aggregates group-level variances up to the estimator scale by weighting
#' each group by the square of its size and the inverse of its first-stage
#' draw probability.
#' @keywords internal
#' @noRd
.vc_within <- function(size, group_var, draw_prob, grand_total) {
  .vc_component(sum(size^2 * group_var / draw_prob), grand_total^2)
}

#' Unit relvariance of the analysis variable
#' @keywords internal
#' @noRd
.vc_unit_relvar <- function(y, w) {
  centre <- if (is.null(w)) mean(y) else sum(w * y) / sum(w)
  spread <- if (is.null(w)) var(y) else .wtdvar(y, w)
  total_w <- if (is.null(w)) length(y) else sum(w)
  eps <- .vc_eps()
  if (abs(centre) < eps && spread < eps) {
    return(.vc_component(0, 1))
  }
  .vc_component(spread * total_w^2, (centre * total_w)^2)
}

#' Homogeneity and ratio parameter for one pair of components
#'
#' `icc` is the share of the pair carried by the first component and
#' `var_ratio` rescales the pair to unit relvariance. Both are ratios of
#' components sharing a denominator, so both are formed from the
#' numerators and stay defined wherever the outcome has variance to split,
#' including at a mean of zero. Only an outcome with no variance leaves
#' them at their neutral values.
#' @keywords internal
#' @noRd
.vc_ratios <- function(first, second, unit) {
  pair <- first$num + second$num
  if (pair <= 0 || unit$num <= 0) {
    return(list(icc = 0, var_ratio = 1, degenerate = TRUE))
  }
  list(icc = first$num / pair, var_ratio = pair / unit$num, degenerate = FALSE)
}

#' Warn once when a component pair could not be identified
#' @keywords internal
#' @noRd
.vc_warn_degenerate <- function(degenerate) {
  if (any(degenerate)) {
    warning("the outcome has no variance to split. 'icc' set to 0 by convention",
            call. = FALSE)
  }
  invisible(NULL)
}

#' Warn when the relvariances are undefined but the ratios are not
#'
#' Relvariance is variance over the square of the mean, so an outcome
#' centred on zero has none, however variable it is. That is a different
#' condition from a constant outcome and it leaves `icc` and `var_ratio`
#' perfectly well defined, since neither depends on the mean.
#' @keywords internal
#' @noRd
.vc_warn_zero_mean <- function(unit) {
  if (unit$num > 0 && !is.finite(unit$rel)) {
    warning(paste("the outcome mean is approximately zero, so the relvariances",
                  "'varb', 'varw' and 'unit_relvar' are infinite; 'icc' and",
                  "'var_ratio' do not depend on the mean and are reported"),
            call. = FALSE)
  }
  invisible(NULL)
}

#' Two-stage variance components
#'
#' Shared by the SRS and PPS first-stage paths, which differ only in how
#' the between-PSU dispersion of the PSU totals is measured and in the
#' first-stage draw probabilities used to scale the within-PSU term.
#' @keywords internal
#' @noRd
.varcomp_2stage <- function(y, psu_id, draw_prob, w, levels) {
  psu <- .vc_group(psu_id, levels)
  .check_min_psu(psu$count)
  summary <- .vc_summarise(y, w, psu)
  summary$var <- .vc_fill_singletons(summary$var)

  correction <- 0
  if (!is.null(w)) {
    correction <- .vc_total_var(summary$observed, summary$size, summary$var)
  }

  if (is.null(draw_prob)) {
    between <- .vc_between_srs(
      summary$total,
      adjustment = if (is.null(w)) 0 else mean(correction)
    )
    scale <- rep(1 / psu$count, psu$count)
  } else {
    between <- .vc_between_pps(
      summary$total, draw_prob,
      adjustment = if (is.null(w)) {
        0
      } else {
        sum(correction * (1 - draw_prob) / draw_prob)
      }
    )
    scale <- draw_prob
  }

  within <- .vc_within(
    summary$size, summary$var, scale, sum(summary$total)
  )
  unit <- .vc_unit_relvar(y, w)
  ratios <- .vc_ratios(between, within, unit)
  .vc_warn_degenerate(ratios$degenerate)
  .vc_warn_zero_mean(unit)

  .new_svyplan_varcomp(
    varb    = between$rel,
    varw    = within$rel,
    icc   = ratios$icc,
    var_ratio       = ratios$var_ratio,
    unit_relvar = unit$rel,
    stages  = 2L
  )
}

#' 2-stage SRS variance components
#' @keywords internal
#' @noRd
.varcomp_2stage_srs <- function(y, psu_id, w = NULL) {
  .varcomp_2stage(y, psu_id, NULL, w, unique(psu_id))
}

#' 2-stage PPS variance components
#' @keywords internal
#' @noRd
.varcomp_2stage_pps <- function(y, psu_id, pp, w = NULL) {
  levels <- sort(unique(psu_id))
  draw_prob <- .varcomp_map_pp(pp, y, psu_id, levels)
  .varcomp_2stage(y, psu_id, draw_prob, w, levels)
}

#' 3-stage PPS variance components
#'
#' Four terms are built from two nested groupings. At the PSU level the
#' between-PSU dispersion of PSU totals gives `varb` and the spread of
#' elements within a PSU gives the companion term the PSU homogeneity is
#' measured against. At the SSU level the spread of SSU totals within a
#' PSU gives `varw_psu` and the spread of elements within an SSU gives
#' `varw_ssu`.
#' @keywords internal
#' @noRd
.varcomp_3stage_pps <- function(y, psu_id, ssu_id, pp, w = NULL) {
  psu_levels <- sort(unique(psu_id))
  psu <- .vc_group(psu_id, psu_levels)
  .check_min_psu(psu$count)
  draw_prob <- .varcomp_map_pp(pp, y, psu_id, psu_levels)

  # SSU labels need only be unique inside a PSU, so nest before grouping.
  nested <- interaction(psu_id, ssu_id, drop = TRUE)
  ssu_levels <- unique(nested)
  ssu <- .vc_group(nested, ssu_levels)
  psu_of_ssu <- match(psu_id[match(ssu_levels, nested)], psu_levels)
  ssu_per_psu <- tabulate(psu_of_ssu, nbins = psu$count)

  by_psu <- .vc_summarise(y, w, psu)
  by_ssu <- .vc_summarise(y, w, ssu)
  by_psu$var <- .vc_fill_singletons(by_psu$var)
  by_ssu$var <- .vc_fill_singletons(by_ssu$var)

  ssu_total_var <- .vc_fill_singletons(
    vapply(split(by_ssu$total, psu_of_ssu), var, numeric(1L), USE.NAMES = FALSE)
  )

  psu_adjustment <- 0
  if (!is.null(w)) {
    # Estimated SSU totals carry element-stage sampling noise. Remove it
    # from both terms built on those totals.
    ssu_noise <- .vc_total_var(by_ssu$observed, by_ssu$size, by_ssu$var)
    grouped_noise <- split(ssu_noise, psu_of_ssu)
    ssu_total_var <- pmax(
      ssu_total_var - vapply(grouped_noise, mean, numeric(1L),
                             USE.NAMES = FALSE),
      0
    )
    psu_noise <- vapply(grouped_noise, sum, numeric(1L), USE.NAMES = FALSE)
    psu_adjustment <- sum(psu_noise * (1 - draw_prob) / draw_prob)
  }

  grand_total <- sum(by_psu$total)
  between_psu <- .vc_between_pps(by_psu$total, draw_prob,
                                 adjustment = psu_adjustment)
  element_in_psu <- .vc_within(by_psu$size, by_psu$var, draw_prob, grand_total)
  between_ssu <- .vc_within(ssu_per_psu, ssu_total_var, draw_prob, grand_total)
  element_in_ssu <- .vc_component(
    sum(
      ssu_per_psu[psu_of_ssu] * by_ssu$size^2 * by_ssu$var /
        draw_prob[psu_of_ssu]
    ),
    grand_total^2
  )

  unit <- .vc_unit_relvar(y, w)
  psu_ratios <- .vc_ratios(between_psu, element_in_psu, unit)
  ssu_ratios <- .vc_ratios(between_ssu, element_in_ssu, unit)
  .vc_warn_degenerate(c(psu_ratios$degenerate, ssu_ratios$degenerate))
  .vc_warn_zero_mean(unit)

  .new_svyplan_varcomp(
    varb    = between_psu$rel,
    varw    = c(varw_psu = between_ssu$rel, varw_ssu = element_in_ssu$rel),
    icc   = c(icc_psu = psu_ratios$icc, icc_ssu = ssu_ratios$icc),
    var_ratio       = c(var_ratio_psu = psu_ratios$var_ratio, var_ratio_ssu = ssu_ratios$var_ratio),
    unit_relvar = unit$rel,
    stages  = 3L
  )
}
