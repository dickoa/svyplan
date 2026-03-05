#' Constrained Stratified Allocation
#'
#' Compute a constrained stratum allocation under a fixed total sample size,
#' target CV, or budget.
#'
#' @param frame For the default method: data frame with one row per stratum.
#'   Must include `N_h` and either `S_h` or `var`.
#'   Optional columns: `mean_h` (or `p_h`), `cost_h`, `max_weight`, `take_all`.
#'   Any unrecognized columns are treated as domain identifiers.
#'   For `svyplan_prec` objects: a precision result from [prec_alloc()].
#' @param ... Additional arguments passed to methods.
#' @param n Total sample size. Specify exactly one of `n`, `cv`, or `budget`.
#' @param cv Target coefficient of variation (requires `mean_h` or `p_h`
#'   in `frame`). When domain columns are present, this target is enforced
#'   in each domain. Specify exactly one of `n`, `cv`, or `budget`.
#' @param budget Total field budget. Specify exactly one of `n`, `cv`,
#'   or `budget`.
#' @param alloc Allocation rule: `"neyman"` (default), `"optimal"`,
#'   `"proportional"`, or `"power"`.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of
#'   per-stratum unit costs, overriding `frame$cost_h`.
#' @param alpha Significance level, default 0.05.
#' @param deff Design effect multiplier (> 0).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1.
#' @param min_n Optional minimum sample size per stratum.
#' @param power_q Bankier power parameter from 0 to 1, used when
#'   `alloc = "power"`.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "alloc"` and a stratum-level
#'   allocation table in `$detail`.
#'
#' @details
#' Allocation is controlled by the `alloc` parameter (same methods as
#' [strata_bound()]):
#' - **proportional**: \eqn{n_h \propto N_h}
#' - **neyman**: \eqn{n_h \propto N_h S_h}
#' - **optimal**: \eqn{n_h \propto N_h S_h / \sqrt{c_h}}
#' - **power**: Bankier (1988), \eqn{n_h \propto S_h N_h^{power\_q}}
#'
#' Stratum allocations are rounded to integers using the ORIC method
#' (Cont and Heidari, 2015). Constraints (`min_n`, `max_weight`, `take_all`)
#' are enforced via recursive Neyman allocation (RNA, Wesolowski et al., 2021).
#'
#' When `cv` is specified with domain columns, the algorithm finds the
#' minimum total `n` such that every domain achieves the target CV.
#'
#' When `budget` is specified, the algorithm finds the maximum affordable
#' allocation under unit costs.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapter 5.
#'
#' Bankier, M. D. (1988). Power allocations: determining sample sizes for
#'   subnational areas. *The American Statistician*, 42(3), 174--177.
#'
#' @seealso [prec_alloc()], [strata_bound()].
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("A", "B", "C"),
#'   N_h = c(4000, 3000, 3000),
#'   S_h = c(10, 15, 8),
#'   mean_h = c(50, 60, 55),
#'   cost_h = c(1, 1.5, 1)
#' )
#'
#' n_alloc(frame, n = 600)
#' n_alloc(frame, cv = 0.03)
#'
#' @export
n_alloc <- function(frame, ...) {
  if (!missing(frame)) {
    .res <- .dispatch_plan(frame, "frame", n_alloc.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_alloc")
}

#' @rdname n_alloc
#' @export
n_alloc.default <- function(
  frame,
  n = NULL,
  cv = NULL,
  budget = NULL,
  alloc = c("neyman", "optimal", "proportional", "power"),
  unit_cost = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  min_n = NULL,
  power_q = 0.5,
  plan = NULL,
  ...
) {
  .plan <- .merge_plan_args(plan, n_alloc.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_alloc.default, c(.plan, list(...))))
  alloc <- match.arg(alloc)
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)

  if (!is.null(min_n)) check_scalar(min_n, "min_n")
  if (alloc == "power") {
    if (!is.numeric(power_q) || length(power_q) != 1L || is.na(power_q) ||
        power_q < 0 || power_q > 1) {
      stop("'power_q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }

  mode_count <- (!is.null(n)) + (!is.null(cv)) + (!is.null(budget))
  if (mode_count != 1L) {
    stop("specify exactly one of 'n', 'cv', or 'budget'", call. = FALSE)
  }

  prep <- .alloc_prepare_frame(frame, unit_cost = unit_cost)
  N_h <- prep$N_h
  S_h <- prep$S_h
  mean_h <- prep$mean_h
  cost_h <- prep$cost_h

  bounds <- .alloc_bounds(
    N_h = N_h,
    max_weight = prep$max_weight,
    take_all = prep$take_all,
    min_n = min_n
  )
  m_h <- bounds$m_h
  M_h <- bounds$M_h

  lo <- sum(m_h)
  hi <- sum(M_h)
  tol <- 1e-8

  a_h <- .alloc_weights(alloc, power_q, N_h, S_h, cost_h)
  if (!is.finite(sum(a_h)) || sum(a_h) <= 0) a_h <- N_h

  mode <- if (!is.null(n)) "n" else if (!is.null(cv)) "cv" else "budget"
  target_total <- NA_real_

  if (mode == "n") {
    check_scalar(n, "n")
    if (n < lo - tol) {
      stop("'n' is below the minimum feasible total under constraints",
           call. = FALSE)
    }
    if (n > hi + tol) {
      stop("'n' exceeds the maximum feasible total (census bound)",
           call. = FALSE)
    }
    target_total <- n
  } else if (mode == "cv") {
    check_scalar(cv, "cv")
    if (anyNA(mean_h)) {
      stop("'mean_h' (or 'p_h') is required in frame when solving for 'cv'",
           call. = FALSE)
    }

    cv_for_total <- function(n_total) {
      n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
      if (length(prep$domain_idx) == 0L) {
        .alloc_metrics(
          N_h = N_h, S_h = S_h, mean_h = mean_h, n_h = n_h,
          alpha = alpha, deff = deff, resp_rate = resp_rate,
          cost_h = cost_h
        )$cv
      } else {
        .alloc_domain_cv_max(
          prep = prep, n_h = n_h,
          alpha = alpha, deff = deff, resp_rate = resp_rate
        )
      }
    }

    cv_lo <- cv_for_total(lo)
    if (cv <= 0) stop("'cv' must be positive", call. = FALSE)
    if (cv_lo <= cv + tol) {
      target_total <- lo
    } else {
      cv_hi <- cv_for_total(hi)
      if (!is.finite(cv_hi) || cv_hi > cv + tol) {
        stop("target 'cv' is unattainable under current constraints",
             call. = FALSE)
      }
      target_total <- uniroot(
        function(x) cv_for_total(x) - cv,
        interval = c(lo, hi), tol = 1e-8
      )$root
    }
  } else {
    check_scalar(budget, "budget")
    cost_lo <- sum(m_h * cost_h)
    cost_hi <- sum(M_h * cost_h)

    if (budget < cost_lo - tol) {
      stop("'budget' is below the minimum feasible cost under constraints",
           call. = FALSE)
    }
    if (budget > cost_hi + tol) {
      stop("'budget' exceeds the maximum feasible cost (census bound)",
           call. = FALSE)
    }

    if (abs(budget - cost_lo) <= tol) {
      target_total <- lo
    } else {
      cost_for_total <- function(n_total) {
        n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
        sum(n_h * cost_h)
      }
      target_total <- uniroot(
        function(x) cost_for_total(x) - budget,
        interval = c(lo, hi), tol = 1e-8
      )$root
    }
  }

  n_h <- .rna_alloc(a_h, target_total, m_h, M_h)
  metrics <- .alloc_metrics(
    N_h = N_h, S_h = S_h, mean_h = mean_h, n_h = n_h,
    alpha = alpha, deff = deff, resp_rate = resp_rate, cost_h = cost_h
  )

  detail <- .alloc_detail(
    prep = prep, n_h = n_h, m_h = m_h, M_h = M_h,
    resp_rate = resp_rate, deff = deff
  )
  domains <- .alloc_domain_summary(
    prep = prep, n_h = n_h,
    alpha = alpha, deff = deff, resp_rate = resp_rate
  )

  params <- list(
    frame = frame,
    alloc = alloc,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = cost_h,
    min_n = min_n,
    domain_cols = prep$domain_cols,
    power_q = if (alloc == "power") power_q else NULL,
    n_h = n_h,
    achieved = list(n = sum(n_h), cv = metrics$cv, cost = metrics$cost)
  )
  params[[mode]] <- switch(mode, n = n, cv = cv, budget = budget)

  obj <- .new_svyplan_n(
    n = sum(n_h),
    type = "alloc",
    method = alloc,
    params = params,
    detail = detail,
    binding = .alloc_binding_label(n_h, m_h, M_h),
    domains = domains
  )
  obj$se <- metrics$se
  obj$moe <- metrics$moe
  obj$cv <- metrics$cv

  obj
}

#' @rdname n_alloc
#' @export
n_alloc.svyplan_prec <- function(
  frame,
  n = NULL,
  cv = NULL,
  budget = NULL,
  ...
) {
  x <- frame
  if (x$type != "alloc") {
    stop("n_alloc requires a svyplan_prec of type 'alloc'", call. = FALSE)
  }
  p <- x$params

  if (is.null(n) && is.null(cv) && is.null(budget)) {
    n <- p$achieved$n
  }

  n_alloc.default(
    frame = p$frame,
    n = n, cv = cv, budget = budget,
    alloc = p$alloc %||% "neyman",
    unit_cost = p$cost_h,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    min_n = p$min_n,
    power_q = p$power_q %||% 0.5
  )
}

#' Precision for a Constrained Allocation
#'
#' Compute aggregate precision for a stratum allocation.
#'
#' @param x For the default method: a frame data frame with stratum metadata.
#'   For `svyplan_n` objects: an allocation result from [n_alloc()].
#' @param ... Additional arguments passed to methods.
#' @param n Stratum sample sizes, length `nrow(x)` (default method only).
#' @param alpha Significance level, default 0.05.
#' @param deff Design effect multiplier (> 0).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1.
#' @param unit_cost Optional scalar or length-`nrow(x)` vector of
#'   per-stratum unit costs, overriding `x$cost_h`.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "alloc"`.
#'
#' @seealso [n_alloc()].
#'
#' @examples
#' frame <- data.frame(
#'   N_h = c(4000, 3000, 3000),
#'   S_h = c(10, 15, 8),
#'   mean_h = c(50, 60, 55)
#' )
#' res <- n_alloc(frame, n = 600)
#' prec_alloc(res)
#'
#' @export
prec_alloc <- function(x, ...) {
  if (!missing(x)) {
    .res <- .dispatch_plan(x, "x", prec_alloc.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_alloc")
}

#' @rdname prec_alloc
#' @export
prec_alloc.default <- function(
  x,
  n,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  unit_cost = NULL,
  plan = NULL,
  ...
) {
  .plan <- .merge_plan_args(plan, prec_alloc.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_alloc.default, c(.plan, list(...))))
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)

  prep <- .alloc_prepare_frame(x, unit_cost = unit_cost)
  H <- length(prep$N_h)

  if (!is.numeric(n) || anyNA(n) || any(!is.finite(n)) || length(n) != H) {
    stop("'n' must be a finite numeric vector with length nrow(x)",
         call. = FALSE)
  }
  if (any(n <= 0)) {
    stop("all 'n' elements must be positive", call. = FALSE)
  }

  metrics <- .alloc_metrics(
    N_h = prep$N_h, S_h = prep$S_h, mean_h = prep$mean_h, n_h = n,
    alpha = alpha, deff = deff, resp_rate = resp_rate,
    cost_h = prep$cost_h
  )

  detail <- .alloc_detail(
    prep = prep, n_h = n,
    m_h = rep(NA_real_, H), M_h = prep$N_h,
    resp_rate = resp_rate, deff = deff
  )

  .new_svyplan_prec(
    se = metrics$se,
    moe = metrics$moe,
    cv = metrics$cv,
    type = "alloc",
    params = list(
      frame = x, n = n, alpha = alpha, deff = deff,
      resp_rate = resp_rate, cost_h = prep$cost_h,
      achieved = list(n = sum(n), cv = metrics$cv, cost = metrics$cost)
    ),
    detail = detail
  )
}

#' @rdname prec_alloc
#' @export
prec_alloc.svyplan_n <- function(x, ...) {
  obj <- x
  if (obj$type != "alloc") {
    stop("prec_alloc requires a svyplan_n of type 'alloc'", call. = FALSE)
  }
  p <- obj$params
  n_h <- p$n_h
  if (is.null(n_h) && !is.null(obj$detail) && "n_h" %in% names(obj$detail)) {
    n_h <- obj$detail$n_h
  }
  if (is.null(n_h)) {
    stop("allocation detail does not contain stratum sample sizes",
         call. = FALSE)
  }

  prec_alloc.default(
    x = p$frame,
    n = n_h,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    unit_cost = p$cost_h
  )
}

#' @keywords internal
#' @noRd
.alloc_prepare_frame <- function(frame, unit_cost = NULL) {
  if (!is.data.frame(frame) || nrow(frame) == 0L) {
    stop("'frame' must be a non-empty data frame", call. = FALSE)
  }

  if (!"N_h" %in% names(frame)) {
    stop("'frame' must contain an 'N_h' column", call. = FALSE)
  }
  N_h <- frame$N_h
  if (!is.numeric(N_h) || anyNA(N_h) || any(!is.finite(N_h)) || any(N_h <= 0)) {
    stop("'N_h' must contain positive finite values", call. = FALSE)
  }

  if ("S_h" %in% names(frame)) {
    S_h <- frame$S_h
    if (!is.numeric(S_h) || anyNA(S_h) || any(!is.finite(S_h)) || any(S_h < 0)) {
      stop("'S_h' must contain non-negative finite values", call. = FALSE)
    }
  } else if ("var" %in% names(frame)) {
    vv <- frame$var
    if (!is.numeric(vv) || anyNA(vv) || any(!is.finite(vv)) || any(vv < 0)) {
      stop("'var' must contain non-negative finite values", call. = FALSE)
    }
    S_h <- sqrt(vv)
  } else {
    stop("'frame' must contain either 'S_h' or 'var'", call. = FALSE)
  }

  mean_h <- rep(NA_real_, nrow(frame))
  if ("mean_h" %in% names(frame)) {
    mean_h <- frame$mean_h
  } else if ("p_h" %in% names(frame)) {
    mean_h <- frame$p_h
  }
  if (any(!is.na(mean_h) & !is.finite(mean_h))) {
    stop("'mean_h' (or 'p_h') must contain finite values", call. = FALSE)
  }

  if (!is.null(unit_cost)) {
    if (!is.numeric(unit_cost) || anyNA(unit_cost) || any(!is.finite(unit_cost)) || any(unit_cost <= 0)) {
      stop("'unit_cost' must contain positive finite values", call. = FALSE)
    }
    if (length(unit_cost) == 1L) {
      cost_h <- rep(unit_cost, nrow(frame))
    } else if (length(unit_cost) == nrow(frame)) {
      cost_h <- unit_cost
    } else {
      stop("'unit_cost' must have length 1 or nrow(frame)", call. = FALSE)
    }
  } else if ("cost_h" %in% names(frame)) {
    cost_h <- frame$cost_h
    if (!is.numeric(cost_h) || anyNA(cost_h) || any(!is.finite(cost_h)) ||
        any(cost_h <= 0)) {
      stop("'cost_h' must contain positive finite values", call. = FALSE)
    }
  } else {
    cost_h <- rep(1, nrow(frame))
  }

  max_weight <- rep(NA_real_, nrow(frame))
  if ("max_weight" %in% names(frame)) {
    max_weight <- frame$max_weight
    bad <- !is.na(max_weight) & (!is.finite(max_weight) | max_weight <= 0)
    if (any(bad)) {
      stop("'max_weight' must be positive and finite when provided",
           call. = FALSE)
    }
  }

  take_all <- rep(FALSE, nrow(frame))
  if ("take_all" %in% names(frame)) {
    ta <- frame$take_all
    if (is.numeric(ta)) ta <- ta != 0
    if (!is.logical(ta) || anyNA(ta)) {
      stop("'take_all' must be logical (or 0/1)", call. = FALSE)
    }
    take_all <- ta
  }

  stratum <- if ("stratum" %in% names(frame)) {
    as.character(frame$stratum)
  } else {
    as.character(seq_len(nrow(frame)))
  }

  known <- c(
    "stratum", "N_h", "S_h", "var", "mean_h", "p_h",
    "cost_h", "max_weight", "take_all"
  )
  domain_cols <- setdiff(names(frame), known)
  domain_idx <- list()
  domain_values <- NULL
  if (length(domain_cols) > 0L) {
    key <- if (length(domain_cols) == 1L) {
      as.character(frame[[domain_cols]])
    } else {
      as.character(interaction(frame[domain_cols], drop = TRUE, sep = ":"))
    }
    lev <- unique(key)
    domain_idx <- setNames(lapply(lev, function(k) which(key == k)), lev)
    domain_values <- frame[match(lev, key), domain_cols, drop = FALSE]
    rownames(domain_values) <- NULL
  }

  list(
    frame = frame,
    N_h = as.numeric(N_h),
    S_h = as.numeric(S_h),
    mean_h = as.numeric(mean_h),
    cost_h = as.numeric(cost_h),
    max_weight = as.numeric(max_weight),
    take_all = as.logical(take_all),
    stratum = stratum,
    domain_cols = domain_cols,
    domain_idx = domain_idx,
    domain_values = domain_values
  )
}

#' @keywords internal
#' @noRd
.alloc_bounds <- function(N_h, max_weight, take_all, min_n = NULL) {
  H <- length(N_h)
  m_h <- pmin(rep(1, H), N_h)
  M_h <- as.numeric(N_h)

  if (!is.null(min_n)) m_h <- pmax(m_h, min_n)

  has_wmax <- !is.na(max_weight)
  if (any(has_wmax)) {
    m_h[has_wmax] <- pmax(m_h[has_wmax], N_h[has_wmax] / max_weight[has_wmax])
  }

  if (any(take_all)) {
    m_h[take_all] <- N_h[take_all]
    M_h[take_all] <- N_h[take_all]
  }

  if (any(m_h > M_h + 1e-8)) {
    stop("constraints are infeasible: lower bounds exceed stratum population",
         call. = FALSE)
  }

  list(m_h = m_h, M_h = M_h)
}

#' @keywords internal
#' @noRd
.alloc_metrics <- function(N_h, S_h, mean_h, n_h, alpha, deff,
                           resp_rate, cost_h) {
  n_eff <- n_h * resp_rate / deff
  W_h <- N_h / sum(N_h)

  term <- numeric(length(n_h))
  good <- n_eff > 0
  fpc <- pmax(0, 1 - n_eff / N_h)
  term[good] <- W_h[good]^2 * S_h[good]^2 * fpc[good] / n_eff[good]

  zero_ok <- !good & S_h == 0
  bad <- !good & !zero_ok
  if (any(bad)) {
    stop(
      "allocation implies zero effective sample in a stratum with positive variability",
      call. = FALSE
    )
  }

  V <- sum(term)
  se <- sqrt(max(V, 0))
  moe <- qnorm(1 - alpha / 2) * se

  cv <- NA_real_
  if (!all(is.na(mean_h))) {
    if (anyNA(mean_h)) {
      stop("'mean_h' (or 'p_h') must be complete to compute aggregate CV",
           call. = FALSE)
    }
    ybar <- sum(W_h * mean_h)
    cv <- if (ybar == 0) Inf else se / abs(ybar)
  }

  list(se = se, moe = moe, cv = cv, cost = sum(n_h * cost_h))
}

#' @keywords internal
#' @noRd
.alloc_detail <- function(prep, n_h, m_h, M_h, resp_rate, deff) {
  out <- data.frame(
    stratum = prep$stratum,
    N_h = prep$N_h,
    S_h = prep$S_h,
    cost_h = prep$cost_h,
    n_h = n_h,
    n_h_int = .round_oric(n_h),
    weight = prep$N_h / n_h,
    n_eff = n_h * resp_rate / deff,
    .lower = m_h,
    .upper = M_h,
    .binding = abs(n_h - m_h) < 1e-6 | abs(n_h - M_h) < 1e-6
  )
  if (any(prep$take_all)) out$take_all <- prep$take_all
  if (!all(is.na(prep$mean_h))) out$mean_h <- prep$mean_h
  out
}

#' @keywords internal
#' @noRd
.alloc_binding_label <- function(n_h, m_h, M_h) {
  on_lower <- abs(n_h - m_h) < 1e-6
  on_upper <- abs(n_h - M_h) < 1e-6
  if (any(on_lower)) return("lower_bound")
  if (any(on_upper)) return("upper_bound")
  "none"
}

#' @keywords internal
#' @noRd
.alloc_domain_cv_max <- function(prep, n_h, alpha, deff, resp_rate) {
  cvs <- vapply(
    prep$domain_idx,
    function(idx) {
      .alloc_metrics(
        N_h = prep$N_h[idx], S_h = prep$S_h[idx],
        mean_h = prep$mean_h[idx], n_h = n_h[idx],
        alpha = alpha, deff = deff, resp_rate = resp_rate,
        cost_h = prep$cost_h[idx]
      )$cv
    },
    numeric(1)
  )
  max(cvs)
}

#' @keywords internal
#' @noRd
.alloc_domain_summary <- function(prep, n_h, alpha, deff, resp_rate) {
  if (length(prep$domain_idx) == 0L) return(NULL)

  out <- vector("list", length(prep$domain_idx))
  dom_names <- names(prep$domain_idx)
  for (i in seq_along(prep$domain_idx)) {
    idx <- prep$domain_idx[[i]]
    met <- .alloc_metrics(
      N_h = prep$N_h[idx], S_h = prep$S_h[idx],
      mean_h = prep$mean_h[idx], n_h = n_h[idx],
      alpha = alpha, deff = deff, resp_rate = resp_rate,
      cost_h = prep$cost_h[idx]
    )
    row <- if (!is.null(prep$domain_values)) {
      prep$domain_values[i, , drop = FALSE]
    } else {
      data.frame()
    }
    row$.domain <- dom_names[i]
    row$.n <- sum(n_h[idx])
    row$.se <- met$se
    row$.moe <- met$moe
    row$.cv <- met$cv
    row$.cost <- met$cost
    out[[i]] <- row
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
