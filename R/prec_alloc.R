#' Precision of a stratified allocation
#'
#' Compute the expected sampling precision of a stratum allocation over the
#' whole population and by domain. This is the inverse of [n_alloc()].
#'
#' @param frame For the default method, a stratum-level data frame as in
#'   [n_alloc()]. For `svyplan_alloc` objects, a result from [n_alloc()].
#'   Assess a joint allocation from [n_alloc_joint()] with [prec_alloc_joint()].
#' @param n Numeric vector of stratum sample sizes in frame row order,
#'   counting issued units. No entry may exceed its stratum's `N`, or, when
#'   the frame gives `N_1`, `N_1` times the product of the stratum's takes.
#'   This holds whatever `fpc` is.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param domains Optional character vector of domain column names in `frame`.
#' @param alpha Significance level, default `0.05`.
#' @param deff Design effect multiplier (`> 0`), scalar or length `nrow(frame)`.
#' @param resp_rate Expected response rate among eligible units in `(0, 1]`.
#' @param eligibility_rate Share of issued units in scope in `(0, 1]`.
#' @param eligible_per_unit Eligible units per sampled unit (`> 0`). These
#'   four combine with frame columns of the same name as in [n_alloc()].
#' @param df Optional degrees of freedom for interval quantiles.
#' @param unit_cost Optional scalar or length-`nrow(frame)` per-unit cost.
#' @param fpc Finite population correction (`"unit"`, `"stage"`, or `"none"`).
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "alloc"`, overall `se`, `moe`,
#'   `cv`, and `rmoe`, stratum table `$detail`, and optional `$domains`.
#'
#' @details
#' Passing a fitted allocation assesses its continuous counts and takes by
#' default. For a clustered field design, supply `n = fit$detail$n_int`
#' and any changed whole takes as named `take_2`, ..., `take_S` arguments
#' in `...`. The whole take columns are `take_2_int`, ..., `take_S_int`.
#'
#' @family stratified design functions
#' @seealso [n_alloc()] for the inverse allocator and [prec_alloc_joint()] for
#'   joint multivariate allocation precision.
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("A", "B", "C"),
#'   N = c(4000, 3000, 3000),
#'   sd = c(10, 15, 8),
#'   mean = c(50, 60, 55)
#' )
#' fit <- n_alloc(frame, n = 600)
#' prec_alloc(fit)
#' prec_alloc(frame, n = c(200, 200, 200))
#'
#' cluster <- transform(frame, N_1 = N / 20, take_2 = 10, icc_1 = 0.05,
#'                      cost_1 = 50, cost_2 = 1)
#' fit <- n_alloc(cluster, n = 600)
#' prec_alloc(fit, n = fit$detail$n_int, take_2 = fit$detail$take_2_int)
#'
#' @export
prec_alloc <- function(frame, ...) {
  if (!missing(frame)) {
    .res <- .dispatch_plan(frame, "frame", prec_alloc.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_alloc")
}

#' @rdname prec_alloc
#' @export
prec_alloc.default <- function(
  frame,
  n = NULL,
  ...,
  domains = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  df = NULL,
  unit_cost = NULL,
  fpc = c("unit", "stage", "none"),
  plan = NULL
) {
  .plan <- .merge_plan_args(
    plan,
    prec_alloc.default,
    match.call(),
    environment()
  )
  if (!is.null(.plan)) {
    return(do.call(prec_alloc.default, c(.plan, list(...))))
  }
  .check_alloc_joint_args(list(...), "prec_alloc_joint()")
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(frame, "frame")
  .refuse_argument_beside_column(
    frame, resp_rate = resp_rate, deff = deff,
    eligibility_rate = eligibility_rate, eligible_per_unit = eligible_per_unit
  )
  fpc <- match.arg(fpc)
  check_alpha(alpha)

  prep <- .alloc_prepare_frame(
    frame,
    domains = domains,
    unit_cost = unit_cost,
    deff = deff,
    resp_rate = resp_rate,
    fpc = fpc,
    eligibility_rate = eligibility_rate,
    eligible_per_unit = eligible_per_unit
  )
  deff_arg <- deff
  resp_rate_arg <- resp_rate
  deff <- prep$deff_h
  resp_rate <- prep$resp_rate_h
  H <- length(prep$N_h)

  if (!is.numeric(n) || anyNA(n) || any(!is.finite(n)) || length(n) != H) {
    stop(
      "'n' must be a finite numeric vector with length nrow(frame)",
      call. = FALSE
    )
  }
  if (any(n <= 0)) {
    stop("all 'n' elements must be positive", call. = FALSE)
  }
  over <- n - prep$N_h > pmin(1e-9 * prep$N_h, 1e-6)
  if (any(over)) {
    stop(sprintf("'n' exceeds the population bound for stratum: %s",
                 paste(prep$stratum[over], collapse = ", ")), call. = FALSE)
  }
  prep <- .cluster_take_rise(prep, n)
  if (isTRUE(prep$cluster) && !is.null(prep$N_psu_h)) {
    cap <- prep$N_psu_h * prep$take_prod_h
    over <- !is.na(cap) & n - cap > pmin(1e-9 * cap, 1e-6)
    if (any(over)) {
      stop(sprintf(
        "'n' needs more PSUs than 'N_1' holds at its takes for stratum: %s",
        paste(prep$stratum[over], collapse = ", ")), call. = FALSE)
    }
  }

  metrics <- .alloc_metrics(
    N_h = prep$N_h,
    N_pop = .alloc_scope_N(prep),
    S_h = prep$S_h,
    mean_h = prep$mean_h,
    n_h = n,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = prep$cost_h,
    df = df,
    N_fpc = prep$N_fpc_h
  )

  detail <- .alloc_detail(
    prep = prep,
    n_h = n,
    m_h = rep(NA_real_, H),
    M_h = prep$N_h,
    resp_rate = resp_rate,
    deff = deff,
    metrics = metrics
  )
  .alloc_warn_psu_room(prep, detail)
  domain_summary <- .alloc_domain_summary(
    prep = prep,
    n_h = n,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  carried <- .alloc_carry_labels(detail, frame)
  detail <- carried$detail

  .new_svyplan_prec(
    se = metrics$se,
    moe = metrics$moe,
    cv = metrics$cv,
    rmoe = metrics$rmoe,
    type = "alloc",
    params = list(
      frame = frame,
      n = n,
      alpha = alpha,
      fpc = fpc,
      deff = deff_arg,
      resp_rate = resp_rate_arg,
      eligibility_rate = if (!identical(eligibility_rate, 1)) eligibility_rate,
      eligible_per_unit = if (!identical(eligible_per_unit, 1)) eligible_per_unit,
      df = df,
      cost_h = prep$cost_h,
      domain_cols = prep$domain_cols,
      domain_idx = prep$domain_idx,
      strata_keys = carried$keys,
      achieved = list(n = sum(n), cv = metrics$cv, cost = metrics$cost)
    ),
    detail = detail,
    domains = domain_summary
  )
}

#' @rdname prec_alloc
#' @export
prec_alloc.svyplan_alloc <- function(frame, ...) {
  obj <- frame
  p <- obj$params
  n_h <- p$n_h
  if (is.null(n_h) && !is.null(obj$detail) && "n" %in% names(obj$detail)) {
    n_h <- obj$detail$n
  }
  if (is.null(n_h)) {
    stop(
      "allocation detail does not contain stratum sample sizes",
      call. = FALSE
    )
  }

  ## Assess at the searched takes, and keep the sizing frame for inversion.
  take <- .alloc_take_override(p$frame, list(...), .alloc_is_cluster(p$frame))
  sizing <- take$frame
  searched <- p$searched_takes
  if (!is.null(searched)) {
    searched <- searched[, setdiff(colnames(searched),
                                   grep("^take_[0-9]+$", names(list(...)),
                                        value = TRUE)), drop = FALSE]
  }
  args <- list(
    frame = .alloc_frame_at_takes(sizing, searched),
    n = n_h,
    domains = p$domain_cols,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    df = p$df
  )
  if (!.alloc_is_cluster(p$frame)) {
    args$unit_cost <- p$cost_h
  }
  args$fpc <- p$fpc %||% "unit"
  out <- do.call(
    prec_alloc.default,
    .roundtrip_args(args, take$dots, prec_alloc.default)
  )
  out$params$alloc <- p$alloc
  out$params$alloc_q <- p$alloc_q
  out$params$min_n_stratum <- p$min_n_stratum
  out$params$min_n_domain <- p[["min_n_domain"]]
  out$params$control <- p[["control"]]
  if (!is.null(searched) && ncol(searched) > 0L) {
    out$params$frame <- sizing
    out$params$searched_takes <- searched
  }
  dots_in <- list(...)
  take_nms <- grep("^take_[0-9]+$", names(dots_in), value = TRUE)
  .mark_continuous_optimum(out, dots_in, c("n", take_nms))
}

#' @rdname prec_alloc
#' @export
prec_alloc.svyplan_alloc_joint <- function(frame, ...) {
  stop(
    "use prec_alloc_joint() for a svyplan_alloc_joint object from n_alloc_joint()",
    call. = FALSE
  )
}
