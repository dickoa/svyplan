#' Sampling precision for a multistage cluster allocation
#'
#' Compute the sampling error (SE, MOE, CV) for a given S-stage cluster sample
#' allocation (S >= 2). This is the inverse of [n_cluster()].
#'
#' @param x First argument used for S3 dispatch. Supply a numeric vector of
#'   per-stage sample sizes (same as `n`), a data frame of indicators (one row
#'   per indicator carrying stage sizes `n_1`, `take_2`, ..., `take_S` and
#'   homogeneities `icc_1`, ..., `icc_{S-1}`), a `svyplan_cluster` allocation
#'   from [n_cluster()], or a [svyplan()] profile. Leave `NULL` when supplying
#'   `n` by name.
#' @param n Numeric vector of per-stage sample sizes of length S >= 2
#'   (outermost stage 1 to innermost stage S). Optional unique names are
#'   preserved. Optional when passed as `x` in the default method.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param icc Numeric vector of stage homogeneities of length `stages - 1`, or
#'   a `svyplan_varcomp` object.
#' @param unit_relvar Unit relvariance, default 1. Required beside a
#'   `svyplan_varcomp` backed out of a design effect and refused beside one
#'   estimated from data. See [n_cluster()].
#' @param var_ratio Ratio of stage component unit variance to the analysis
#'   variable's unit variance, default 1. A scalar supplies stage 1 and derives
#'   deeper ratios via `k_s = k_{s-1} * (1 - icc_{s-1})`. In a vector, an `NA`
#'   after the first entry is derived the same way. Refused when `icc` is a
#'   `svyplan_varcomp`.
#' @param resp_rate Expected response rate(s) in `(0, 1]`, default 1. A scalar
#'   applies to the final stage (ultimate units) with earlier stages set to 1.
#'   A vector of length `stages` supplies stage-1 through stage-S rates. Its
#'   names, if given, place each rate, as in [n_cluster()].
#' @param eligibility_rate Share of issued final-stage units in scope, in
#'   `(0, 1]`. Default 1.
#' @param eligible_per_unit Eligible units per final-stage unit, a positive
#'   number, default 1. For the `data.frame` method, an `eligible_per_unit`
#'   column sets it per row with `NA` falling back to this scalar.
#' @param N Optional numeric vector of length `stages` giving aggregate
#'   population counts from stage 1 to stage S (outermost to innermost). Its
#'   names, if given, place each count, as in [n_cluster()].
#' @param fpc Which finite population correction the variance carries.
#'   `"auto"` (default) applies stage-by-stage correction (`"stage"`) whenever
#'   `N` (or population columns) is supplied, and `"none"` otherwise.
#'   `"stage"` requires population counts. Supply population counts only
#'   with `"auto"` or `"stage"`.
#' @param domains Optional character vector naming domain columns in `x`
#'   (data frame method).
#' @param stage_cost Optional per-stage costs (length `stages`), recorded on
#'   the result so a later round trip to [n_cluster()] can re-solve the
#'   design. Costs do not enter the precision calculation.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "cluster"` and components
#'   `$se`, `$moe`, and `$cv`. In the single-indicator method, only `$cv` is
#'   computable (`$se` and `$moe` are `NA`). With an indicator data frame,
#'   per-indicator precision is in `$detail` (`.se`, `.moe`, `.rmoe`, `.cv`).
#'   `$params` stores the stage response rates as `resp_rate_1`, ...,
#'   `resp_rate_{S-1}` and `resp_rate`, the names the allocation frames use.
#'
#' @details
#' For an \eqn{S}-stage design with gross stage sizes \eqn{n_1, \dots, n_S},
#' effective stage response multipliers \eqn{r_1, \dots, r_S}, stage
#' variance components \eqn{R_s = V k_s \delta_s} (\eqn{s < S}) and
#' \eqn{R_S = V k_{S-1}(1 - \delta_{S-1})}{R_S = V k_(S-1) (1 - delta_(S-1))},
#' and per-parent population capacities \eqn{M_1 = N_1},
#' \eqn{M_s = N_s / N_{s-1}}{M_s = N_s / N_(s-1)} for \eqn{1 < s < S}{1 < s < S}
#' and \eqn{M_S = N e / N_{S-1}}{M_S = N e / N_(S-1)}, where \eqn{e} is
#' `eligibility_rate * eligible_per_unit` (all infinite under
#' `fpc = "none"`),
#' \deqn{CV = \sqrt{\sum_{s=1}^{S} \frac{R_s}{\prod_{j=1}^{s} n_j r_j}
#'   \left(1 - \frac{n_s r_s}{M_s}\right)}.}{CV = sqrt(sum_{s=1}^S (R_s / prod_{j=1}^s (n_j r_j)) (1 - n_s r_s / M_s)).}
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer. Chapter 9.
#'
#' @family cluster design functions
#' @seealso [n_cluster()] for the inverse operation, [varcomp()] for
#'   estimating variance components, and [prec_multi()] for single-stage
#'   multi-indicator precision.
#'
#' @examples
#' # Direct usage: 2-stage, 3-stage, and 4-stage
#' prec_cluster(c(50, 12), icc = 0.05)
#' prec_cluster(c(50, 12, 8), icc = c(0.01, 0.05))
#' prec_cluster(c(15, 5, 4, 8), icc = c(0.01, 0.03, 0.06))
#'
#' # Round-trip from n_cluster
#' res <- n_cluster(c(500, 50), icc = 0.05, cv = 0.05)
#' prec_cluster(res)
#' prec_cluster(res$operational$n, icc = 0.05)
#'
#' # Stage finite population correction (auto-enabled when N is supplied)
#' prec_cluster(c(12, 20), icc = 0.05, N = c(16, 2400))
#'
#' # Several indicators under one 2-stage allocation
#' indicators <- data.frame(
#'   name   = c("stunting", "anemia"),
#'   p      = c(0.30, 0.10),
#'   n_1    = c(60, 60),
#'   take_2 = c(12, 12),
#'   icc_1  = c(0.02, 0.05)
#' )
#' prec_cluster(indicators)
#'
#' @export
prec_cluster <- function(x = NULL, ...) {
  if ("indicators" %in% ...names()) {
    stop(
      "pass the indicator table as the first argument: prec_cluster(indicators, ...)",
      call. = FALSE
    )
  }
  if (!missing(x)) {
    .res <- .dispatch_plan(x, "x", prec_cluster, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_cluster")
}

#' @rdname prec_cluster
#' @export
prec_cluster.default <- function(
  x = NULL,
  ...,
  n = NULL,
  icc = NULL,
  unit_relvar = NULL,
  var_ratio = NULL,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  N = NULL,
  fpc = c("auto", "none", "stage"),
  stage_cost = NULL,
  plan = NULL
) {
  dots <- list(...)
  if ("indicators" %in% names(dots)) {
    stop(
      "pass the indicator table as the first argument: prec_cluster(indicators, ...)",
      call. = FALSE
    )
  }
  if ("domains" %in% names(dots)) {
    stop(
      "'domains' applies only when the first argument to prec_cluster() is an indicator data frame",
      call. = FALSE
    )
  }
  if (!is.null(x)) {
    if (!is.null(n)) {
      stop(
        "supply per-stage sample sizes either as the first argument or via 'n', not both",
        call. = FALSE
      )
    }
    n <- x
    x <- NULL
  }
  fpc <- match.arg(fpc)
  n_len <- if (is.numeric(n)) length(n) else NULL
  plan <- .plan_for_cluster_mode(plan, FALSE, fpc, stages = n_len)
  plan <- .check_varcomp_components(icc, unit_relvar, var_ratio, plan,
                                    "prec_cluster()")
  .check_relvar_identified(
    icc,
    !is.null(unit_relvar) || "unit_relvar" %in% names(plan$defaults),
    "prec_cluster()"
  )
  .plan <- .merge_plan_args(plan, prec_cluster.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_cluster.default, c(.plan, dots)))
  .check_unused_dots(...)

  if (!is.numeric(n) || length(n) < 2L) {
    stop("'n' must be a numeric vector of length >= 2", call. = FALSE)
  }
  if (anyNA(n)) stop("'n' must not contain NA values", call. = FALSE)
  if (any(!is.finite(n))) stop("'n' must contain only finite values", call. = FALSE)
  if (any(n <= 0)) stop("all elements of 'n' must be positive", call. = FALSE)
  n <- .reorder_n_vec(n)
  stages <- length(n)
  stage_nms <- .stage_names(stages, names(n))
  names(n) <- stage_nms

  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)
  resp_vec <- .expand_stage_resp_rate(resp_rate, stages, stage_nms)
  fpc_resolved <- .resolve_stage_fpc(fpc, N)
  yield <- eligibility_rate * eligible_per_unit
  pop <- if (identical(fpc_resolved, "stage")) {
    .stage_population(N, stages, yield = yield, stage_nms = stage_nms)
  } else {
    NULL
  }

  unit_relvar <- unit_relvar %||% 1
  var_ratio <- var_ratio %||% 1
  if (is.null(icc)) {
    stop("'icc' is required (directly or via plan)", call. = FALSE)
  }
  if (inherits(icc, "svyplan_varcomp")) {
    vc <- icc
    if (!is.null(vc$strata)) {
      stop(
        "stratified varcomp: merge its $strata columns into an n_alloc() frame, or pass one stratum's values",
        call. = FALSE
      )
    }
    if (!is.null(pop)) .check_stage_fpc_varcomp(vc)
    icc <- vc$icc
    var_ratio <- vc$var_ratio
    if (!is.na(vc$unit_relvar)) unit_relvar <- vc$unit_relvar
  }

  icc <- .reorder_stage_vec(icc, "icc")
  var_ratio <- .reorder_stage_vec(var_ratio, "var_ratio")
  check_icc(icc, expected_length = stages - 1L)
  check_scalar(unit_relvar, "unit_relvar")
  .check_var_ratio(var_ratio)
  var_ratio_input <- .var_ratio_given(var_ratio, length(icc))
  var_ratio <- .expand_var_ratio(var_ratio, icc)

  if (!is.null(stage_cost)) {
    check_stage_cost(stage_cost)
    stage_cost <- .reorder_stage_cost(stage_cost)
    if (length(stage_cost) != stages) {
      stop(sprintf("'stage_cost' must have length %d", stages), call. = FALSE)
    }
    names(stage_cost) <- stage_nms
  }

  rates <- resp_vec
  rates[stages] <- rates[stages] * yield
  rate_label <- paste(
    c("resp_rate", if (eligibility_rate != 1) "eligibility_rate",
      if (eligible_per_unit != 1) "eligible_per_unit"),
    collapse = " * "
  )
  for (s in 2L:stages) {
    .check_expected_take(
      n[[s]], rates[s], stage_nms[s],
      if (s == stages) rate_label else sprintf("resp_rate[%d]", s)
    )
  }

  if (!is.null(pop)) .check_stage_sizes(n, pop, stage_nms)

  R_mat <- .stage_R_matrix(icc, unit_relvar, var_ratio)
  rates_mat <- matrix(rates, nrow = 1L)
  M_mat <- matrix(if (!is.null(pop)) pop$M else rep(Inf, stages), nrow = 1L)
  cv_val <- sqrt(.stage_cv2(unname(n), R_mat, rates_mat, M_mat)[1L])

  names(icc) <- paste0("icc_", seq_len(stages - 1L))
  names(var_ratio) <- paste0("var_ratio_", seq_len(stages - 1L))
  params <- list(
    n = n,
    icc = icc,
    unit_relvar = unit_relvar,
    var_ratio = var_ratio,
    resp_rate = .name_stage_resp_rate(resp_vec),
    stages = stages
  )
  if (eligibility_rate != 1) params$eligibility_rate <- eligibility_rate
  if (eligible_per_unit != 1) params$eligible_per_unit <- eligible_per_unit
  if (!is.null(stage_cost)) params$stage_cost <- stage_cost
  if (!is.null(pop)) {
    params$fpc <- "stage"
    params$N <- pop$N
  }
  if (anyNA(var_ratio_input)) params$var_ratio_input <- var_ratio_input

  .new_svyplan_prec(
    se = NA_real_,
    moe = NA_real_,
    cv = cv_val,
    type = "cluster",
    params = params
  )
}

#' @rdname prec_cluster
#' @export
prec_cluster.svyplan_cluster <- function(x, ...) {
  if (!is.null(x[["indicators"]])) {
    return(.prec_cluster_from_indicator_fit(x, list(...)))
  }
  .prec_cluster_from_fit(x, list(...))
}

#' Read a single-indicator cluster allocation back as precision
#' @keywords internal
#' @noRd
.prec_cluster_from_fit <- function(x, dots, field = FALSE) {
  p <- x$params
  args <- list(
    n = if (field) x$operational$n else x$n,
    icc = p$icc,
    unit_relvar = p$unit_relvar,
    var_ratio = p[["var_ratio_input"]] %||% p$var_ratio,
    resp_rate = p[["resp_rate"]] %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    stage_cost = p$stage_cost,
    N = p[["N"]],
    fpc = p[["fpc"]] %||% "auto"
  )
  out <- do.call(
    prec_cluster.default,
    .roundtrip_args(args, dots, prec_cluster.default)
  )
  if (!is.null(p$budget)) out$params$budget <- p$budget
  if (!is.null(p$stage_n)) out$params$stage_n <- p$stage_n
  if (!is.null(p$fixed_cost)) out$params$fixed_cost <- p$fixed_cost
  if (!is.null(p[["control"]])) out$params$control <- p[["control"]]
  if (field) out else .mark_continuous_optimum(out, dots)
}

#' Refuse scalar inputs that the indicator frame carries per row
#' @keywords internal
#' @noRd
.check_prec_cluster_indicator_args <- function(n, icc, unit_relvar, var_ratio) {
  cols <- c(
    n = "n_1", icc = "icc_1", unit_relvar = "unit_relvar",
    var_ratio = "var_ratio_1"
  )
  supplied <- !vapply(
    list(n = n, icc = icc, unit_relvar = unit_relvar, var_ratio = var_ratio),
    is.null,
    logical(1L)
  )
  if (any(supplied)) .stop_carried_by_column(names(cols)[supplied], cols)
  invisible(TRUE)
}

#' Read a several-indicators cluster allocation back as precision
#' @keywords internal
#' @noRd
.prec_cluster_from_indicator_fit <- function(x, dots, field = FALSE) {
  tgt <- x$indicators
  stages <- x$stages
  stage_nms <- names(x$n) %||% paste0("stage_", seq_len(stages))
  if ("resp_rate" %in% names(dots)) {
    for (s in seq_len(stages - 1L)) {
      col <- paste0("resp_rate_", s)
      if (col %in% names(tgt)) tgt[[col]] <- NA_real_
    }
    if ("resp_rate" %in% names(tgt)) tgt$resp_rate <- NA_real_
  }
  if ("eligible_per_unit" %in% names(dots) &&
      "eligible_per_unit" %in% names(tgt)) {
    tgt$eligible_per_unit <- NA_real_
  }
  tgt$cv <- NULL
  tgt$moe <- NULL
  stage_n <- if (field) x$operational$n else x$n
  tgt$n_1 <- unname(stage_n[1L])
  for (s in 2L:stages) {
    tgt[[paste0("take_", s)]] <- unname(stage_n[s])
  }

  dom_cols <- x$params$domain_cols %||% character(0)
  if (!is.null(x$domains) && length(dom_cols) > 0L) {
    dom <- if (field) x$operational$domains else x$domains
    tgt_key <- .domain_key(tgt, dom_cols)
    dom_key <- .domain_key(dom, dom_cols)
    dom_idx <- match(tgt_key, dom_key)
    tgt$n_1 <- dom[[stage_nms[1L]]][dom_idx]
    for (s in 2L:stages) {
      tgt[[paste0("take_", s)]] <- dom[[stage_nms[s]]][dom_idx]
    }
  }

  args <- list(
    indicators = tgt,
    stages = stages,
    stage_cost = x$params$stage_cost,
    domain_cols = dom_cols,
    fpc = x$params[["fpc"]] %||% "auto",
    N_arg = x$params[["N"]],
    resp_rate = x$params[["resp_rate"]] %||% 1,
    eligibility_rate = x$params[["eligibility_rate"]] %||% 1,
    eligible_per_unit = x$params[["eligible_per_unit"]] %||% 1
  )
  allowed <- c(
    "stage_cost", "domains", "resp_rate", "eligibility_rate",
    "eligible_per_unit"
  )
  dot_names <- names(dots) %||% rep("", length(dots))
  unknown <- setdiff(dot_names, allowed)
  if (length(unknown) > 0L || any(!nzchar(dot_names))) {
    do.call(.check_unused_dots, dots)
  }
  if ("stage_cost" %in% names(dots)) {
    check_stage_cost(dots$stage_cost)
    override_cost <- .reorder_stage_cost(dots$stage_cost)
    if (length(override_cost) != stages) {
      stop(sprintf("'stage_cost' must have length %d", stages), call. = FALSE)
    }
    args$stage_cost <- override_cost
  }
  if ("domains" %in% names(dots)) {
    if (!is.null(dots$domains) &&
        (!is.character(dots$domains) || anyNA(dots$domains) ||
         any(!dots$domains %in% names(tgt)))) {
      stop("'domains' must name columns in indicators", call. = FALSE)
    }
    args$domain_cols <- dots$domains %||% character(0)
  }
  if ("resp_rate" %in% names(dots)) {
    args$resp_rate <- .expand_stage_resp_rate(dots$resp_rate, stages)
  }
  if ("eligibility_rate" %in% names(dots)) {
    check_resp_rate(dots$eligibility_rate, "eligibility_rate")
    args$eligibility_rate <- dots$eligibility_rate
  }
  if ("eligible_per_unit" %in% names(dots)) {
    check_eligible_per_unit(dots$eligible_per_unit)
    args$eligible_per_unit <- dots$eligible_per_unit
  }
  res <- do.call(.prec_multi_cluster, args)
  for (p in c("budget", "stage_n", "joint", "allocation", "fixed_cost",
              "min_n_domain", "mode", "control")) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res$params$cv_target <- x$indicators$cv
  res
}
