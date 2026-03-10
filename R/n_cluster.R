#' Optimal Multistage Cluster Allocation
#'
#' Compute optimal per-stage sample sizes for a multistage cluster
#' design, minimizing cost for a given precision or minimizing
#' variance for a given budget.
#'
#' @param stage_cost For the default method: numeric vector of per-stage
#'   costs. Length determines the number of stages (2 or 3). Named vectors
#'   are accepted with stage names `cost_psu`, `cost_ssu`, `cost_tsu`
#'   (`cost_tsu` aliases `cost_ssu` in 2-stage).
#'   For `svyplan_prec` objects: a precision result from [prec_cluster()].
#' @param ... Additional arguments passed to methods.
#' @param delta Numeric vector of homogeneity measures (length = stages - 1),
#'   or a `svyplan_varcomp` object.
#' @param rel_var Unit relvariance (default 1).
#' @param k Ratio parameter(s). Scalar for 2-stage, length-2 vector for
#'   3-stage (default 1).
#' @param cv Target coefficient of variation. Specify exactly one of `cv`
#'   or `budget`.
#' @param budget Total budget. Specify exactly one of `cv` or `budget`.
#' @param n_psu Fixed number of PSUs (stage-1 sample size). `NULL` (default)
#'   means optimize. For 2-stage, at most one of `n_psu` or `psu_size` may be
#'   specified. For 3-stage, up to two of `n_psu`, `psu_size`, `ssu_size` may
#'   be fixed.
#' @param psu_size Fixed cluster size (stage-2 sample size per PSU). `NULL`
#'   (default) means optimize. This is the typical MICS/DHS parameterization
#'   where the number of households per cluster is fixed.
#' @param ssu_size Fixed SSU take size (stage-3 sample size per SSU). `NULL`
#'   (default) means optimize. Only valid for 3-stage designs.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The stage-1 sample size is inflated by `1 / resp_rate`.
#' @param fixed_cost Fixed overhead cost (C0). Default 0.
#'   The total cost model becomes
#'   `C = C0 + c1*n_psu + c2*n_psu*psu_size [+ c3*n_psu*psu_size*ssu_size]`.
#'   In budget mode, only `budget - fixed_cost` is available for variable costs;
#'   in CV mode, `fixed_cost` is added to the variable cost.
#' @param plan Optional [svyplan()] object providing design defaults
#'   (including `stage_cost`, `delta`, `rel_var`, `k`, `resp_rate`, `fixed_cost`).
#'
#' @return A `svyplan_cluster` object with components:
#' \describe{
#'   \item{`n`}{Named numeric vector of continuous per-stage sample sizes
#'     (e.g. `c(n_psu = 84.1, psu_size = 13.8)`). Use `ceiling()` for
#'     operational (integer) values.}
#'   \item{`stages`}{Number of stages (2 or 3).}
#'   \item{`total_n`}{Continuous total sample size (`prod(n)`). Use
#'     `as.integer()` for the operational total (product of ceiled stages),
#'     or `as.double()` for this continuous value.}
#'   \item{`cv`}{Achieved coefficient of variation (based on continuous
#'     optimum).}
#'   \item{`cost`}{Total cost.}
#'   \item{`params`}{List of input parameters.}
#' }
#'
#' @details
#' Stage count is determined by `length(stage_cost)`. Two dispatch dimensions:
#'
#' - **2-stage vs 3-stage** (vector length)
#' - **budget vs cv** mode (which is non-NULL)
#'
#' One or more stage sizes can be fixed, leaving the remaining stage(s) to be
#' optimized or derived from the constraint. For 2-stage designs, at most one
#' stage may be fixed. For 3-stage designs, up to two stages may be fixed;
#' the remaining free stage is derived from the budget or CV constraint.
#'
#' If `delta` is a `svyplan_varcomp` object, `delta`, `rel_var`, and `k`
#' are extracted automatically.
#'
#' Boundary and near-boundary homogeneity values are not supported by the
#' analytical optimum used here. When `delta` is near 0, most variability is
#' within PSUs, so the closed-form optimum collapses toward taking many units
#' in very few PSUs. When `delta` is near 1, most variability is between PSUs,
#' so the optimum collapses toward taking very few units in many PSUs. In both
#' cases the analytical allocation becomes degenerate, so `n_cluster()`
#' rejects values numerically too close to 0 or 1.
#'
#' These functions assume sampling fractions are negligible at each stage
#' (equivalent to sampling with replacement). No finite population correction
#' is applied. This is standard for multistage planning when cluster
#' populations are large relative to the sample.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 9.
#'
#' @seealso [prec_cluster()] for the inverse, [varcomp()] for estimating
#'   variance components.
#'
#' @examples
#' # 2-stage, budget mode
#' n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000)
#'
#' # 2-stage, CV mode
#' n_cluster(stage_cost = c(500, 50), delta = 0.05, cv = 0.05)
#'
#' # 2-stage, fixed n_psu
#' n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000, n_psu = 40)
#'
#' # 2-stage, fixed psu_size (MICS/DHS style: 20 households per cluster)
#' n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000, psu_size = 20)
#'
#' # 3-stage
#' n_cluster(stage_cost = c(500, 100, 50), delta = c(0.01, 0.05), cv = 0.05)
#'
#' # 3-stage, fixed n_psu + ssu_size (solve for psu_size)
#' n_cluster(
#'   stage_cost = c(500, 100, 50), delta = c(0.01, 0.05),
#'   budget = 500000, n_psu = 50, ssu_size = 8
#' )
#'
#' # With fixed overhead cost
#' n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000, fixed_cost = 5000)
#'
#' @export
n_cluster <- function(stage_cost, ...) {
  if (!missing(stage_cost)) {
    .res <- .dispatch_plan(stage_cost, "stage_cost", n_cluster.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_cluster")
}

#' @rdname n_cluster
#' @export
n_cluster.default <- function(
  stage_cost = NULL,
  delta = NULL,
  rel_var = 1,
  k = 1,
  cv = NULL,
  budget = NULL,
  n_psu = NULL,
  psu_size = NULL,
  ssu_size = NULL,
  resp_rate = 1,
  fixed_cost = 0,
  plan = NULL,
  ...
) {
  .plan <- .merge_plan_args(plan, n_cluster.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_cluster.default, c(.plan, list(...))))
  if (is.null(stage_cost))
    stop("'stage_cost' is required (directly or via plan)", call. = FALSE)
  if (is.null(delta))
    stop("'delta' is required (directly or via plan)", call. = FALSE)
  check_stage_cost(stage_cost)
  stage_cost <- .reorder_stage_cost(stage_cost)
  check_resp_rate(resp_rate)

  if (inherits(delta, "svyplan_varcomp")) {
    vc <- delta
    delta <- vc$delta
    rel_var <- vc$rel_var
    k <- vc$k
  }

  check_scalar(rel_var, "rel_var")
  if (
    !is.numeric(k) ||
      length(k) == 0L ||
      anyNA(k) ||
      any(k <= 0) ||
      any(!is.finite(k))
  ) {
    stop("'k' must contain positive finite values", call. = FALSE)
  }

  stages <- length(stage_cost)
  delta <- .reorder_stage_vec(delta, "delta")
  k <- .reorder_stage_vec(k, "k")
  check_delta(delta, expected_length = stages - 1L)
  .check_cluster_delta_open(delta, context = "n_cluster()")

  has_cv <- !is.null(cv)
  has_budget <- !is.null(budget)
  if (has_cv == has_budget) {
    stop("specify exactly one of 'cv' or 'budget'", call. = FALSE)
  }
  if (has_cv) {
    check_scalar(cv, "cv")
  }
  if (has_budget) {
    check_scalar(budget, "budget")
  }
  if (!is.null(n_psu)) check_scalar(n_psu, "n_psu")
  if (!is.null(psu_size)) check_scalar(psu_size, "psu_size")
  if (!is.null(ssu_size)) check_scalar(ssu_size, "ssu_size")
  if (stages == 2L && !is.null(ssu_size)) {
    stop("'ssu_size' is not applicable for 2-stage designs", call. = FALSE)
  }
  n_fixed <- sum(!is.null(n_psu), !is.null(psu_size), !is.null(ssu_size))
  if (n_fixed >= stages) {
    stop("cannot fix all stages; use prec_cluster() instead", call. = FALSE)
  }
  check_fixed_cost(fixed_cost, budget)

  res <- if (stages == 2L) {
    .n_cluster_2stage(
      stage_cost,
      delta,
      rel_var,
      k,
      cv,
      budget,
      n_psu,
      psu_size,
      resp_rate,
      fixed_cost
    )
  } else {
    k <- rep_len(k, 2L)
    .n_cluster_3stage(
      stage_cost,
      delta,
      rel_var,
      k,
      cv,
      budget,
      n_psu,
      psu_size,
      ssu_size,
      resp_rate,
      fixed_cost
    )
  }
  res
}

#' @rdname n_cluster
#' @export
n_cluster.svyplan_prec <- function(stage_cost, cv = NULL, budget = NULL, ...) {
  x <- stage_cost
  if (x$type != "cluster") {
    stop("n_cluster requires a svyplan_prec of type 'cluster'", call. = FALSE)
  }
  p <- x$params
  if (is.null(cv) && is.null(budget)) {
    cv <- x$cv
  }
  n_cluster.default(
    stage_cost = p$stage_cost,
    delta = p$delta,
    rel_var = p$rel_var,
    k = p$k,
    cv = cv,
    budget = budget,
    n_psu = p$n_psu,
    psu_size = p$psu_size,
    ssu_size = p$ssu_size,
    resp_rate = p$resp_rate %||% 1,
    fixed_cost = p$fixed_cost %||% 0
  )
}

#' @keywords internal
#' @noRd
.n_cluster_2stage <- function(
  stage_cost,
  delta,
  rel_var,
  k,
  cv,
  budget,
  n_psu,
  psu_size,
  resp_rate,
  fixed_cost = 0
) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  var_budget <- if (!is.null(budget)) budget - fixed_cost else NULL

  if (is.null(n_psu) && is.null(psu_size)) {
    n2_opt <- sqrt(C1 / C2 * (1 - delta) / delta)

    if (!is.null(budget)) {
      n1_opt <- var_budget / (C1 + C2 * n2_opt)
      n1_eff <- n1_opt * resp_rate
      cv_achieved <- sqrt(
        rel_var / (n1_eff * n2_opt) * k * (1 + delta * (n2_opt - 1))
      )
      total_cost <- budget
    } else {
      n1_eff_needed <- rel_var *
        k *
        (1 + delta * (n2_opt - 1)) /
        (n2_opt * cv^2)
      n1_opt <- n1_eff_needed / resp_rate
      total_cost <- fixed_cost + C1 * n1_opt + C2 * n1_opt * n2_opt
      cv_achieved <- cv
    }
  } else if (!is.null(n_psu)) {
    n1_opt <- n_psu
    n1_eff <- n_psu * resp_rate
    if (!is.null(budget)) {
      n2_opt <- (var_budget - C1 * n_psu) / (C2 * n_psu)
      if (n2_opt <= 0) {
        stop(
          "budget is too small for the given fixed stage-1 size",
          call. = FALSE
        )
      }
      cv_achieved <- sqrt(
        rel_var * k / (n1_eff * n2_opt) * (1 + delta * (n2_opt - 1))
      )
      total_cost <- budget
    } else {
      n2_opt <- (1 - delta) / (cv^2 * n1_eff / (rel_var * k) - delta)
      if (n2_opt <= 0) {
        stop(
          "target CV is too small for the given fixed stage-1 size and parameters",
          call. = FALSE
        )
      }
      total_cost <- fixed_cost + C1 * n_psu + C2 * n_psu * n2_opt
      cv_achieved <- cv
    }
  } else {
    n2_opt <- psu_size
    if (!is.null(budget)) {
      n1_opt <- var_budget / (C1 + C2 * psu_size)
      n1_eff <- n1_opt * resp_rate
      cv_achieved <- sqrt(
        rel_var * k / (n1_eff * n2_opt) * (1 + delta * (n2_opt - 1))
      )
      total_cost <- budget
    } else {
      n1_eff_needed <- rel_var *
        k *
        (1 + delta * (n2_opt - 1)) /
        (n2_opt * cv^2)
      n1_opt <- n1_eff_needed / resp_rate
      total_cost <- fixed_cost + C1 * n1_opt + C2 * n1_opt * n2_opt
      cv_achieved <- cv
    }
  }

  n_vec <- c(n_psu = n1_opt, psu_size = n2_opt)
  total_n <- prod(n_vec)

  params <- list(
    stage_cost = c(cost_psu = stage_cost[1L], cost_ssu = stage_cost[2L]),
    delta = c(delta_psu = delta[1L]),
    rel_var = rel_var,
    k = c(k_psu = k[1L]),
    resp_rate = resp_rate
  )
  if (!is.null(cv)) {
    params$cv <- cv
  }
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) {
    params$n_psu <- n_psu
  }
  if (!is.null(psu_size)) {
    params$psu_size <- psu_size
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 2L,
    total_n = total_n,
    cv = cv_achieved,
    cost = total_cost,
    params = params
  )
}

#' @keywords internal
#' @noRd
.n_cluster_3stage <- function(
  stage_cost,
  delta,
  rel_var,
  k,
  cv,
  budget,
  n_psu,
  psu_size,
  ssu_size,
  resp_rate,
  fixed_cost = 0
) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  delta1 <- delta[1L]
  delta2 <- delta[2L]
  k1 <- k[1L]
  k2 <- k[2L]
  var_budget <- if (!is.null(budget)) budget - fixed_cost else NULL

  .cv3 <- function(n1e, n2v, n3v) {
    sqrt(
      rel_var / (n1e * n2v * n3v) *
        (k1 * delta1 * n2v * n3v + k2 * (1 + delta2 * (n3v - 1)))
    )
  }

  solve_for <- if (!is.null(n_psu) && !is.null(psu_size)) {
    "n3"
  } else if (!is.null(n_psu)) {
    "n2"
  } else {
    "n1"
  }

  n3 <- if (!is.null(ssu_size)) {
    ssu_size
  } else if (solve_for != "n3") {
    sqrt((1 - delta2) / delta2 * C2 / C3)
  }

  n2 <- if (!is.null(psu_size)) {
    psu_size
  } else if (solve_for != "n2") {
    if (!is.null(ssu_size)) {
      sqrt(
        k2 * (1 + delta2 * (n3 - 1)) * C1 /
          (n3 * k1 * delta1 * (C2 + C3 * n3))
      )
    } else {
      1 / n3 * sqrt((1 - delta2) / delta1 * C1 / C3 * k2 / k1)
    }
  }

  n1 <- n_psu

  if (solve_for == "n1") {
    if (!is.null(budget)) {
      n1 <- var_budget / (C1 + C2 * n2 + C3 * n2 * n3)
      cv_achieved <- .cv3(n1 * resp_rate, n2, n3)
      total_cost <- budget
    } else {
      n1_eff <- rel_var / (cv^2 * n2 * n3) *
        (k1 * delta1 * n2 * n3 + k2 * (1 + delta2 * (n3 - 1)))
      n1 <- n1_eff / resp_rate
      total_cost <- fixed_cost +
        C1 * n1 + C2 * n1 * n2 + C3 * n1 * n2 * n3
      cv_achieved <- cv
    }
  } else if (solve_for == "n2") {
    n1_eff <- n_psu * resp_rate
    if (!is.null(budget)) {
      n2 <- (var_budget / n_psu - C1) / (C2 + C3 * n3)
      if (n2 <= 0) {
        stop(
          "budget is too small for the given fixed stage sizes",
          call. = FALSE
        )
      }
      cv_achieved <- .cv3(n1_eff, n2, n3)
      total_cost <- budget
    } else {
      n2 <- k2 * (1 + delta2 * (n3 - 1)) /
        (n3 * (cv^2 * n1_eff / rel_var - k1 * delta1))
      if (n2 <= 0) {
        stop(
          "target CV is too small for the given fixed stage sizes and parameters",
          call. = FALSE
        )
      }
      total_cost <- fixed_cost +
        C1 * n_psu + C2 * n_psu * n2 + C3 * n_psu * n2 * n3
      cv_achieved <- cv
    }
  } else {
    n1_eff <- n_psu * resp_rate
    if (!is.null(budget)) {
      n3 <- (var_budget - C1 * n_psu - C2 * n_psu * n2) / (C3 * n_psu * n2)
      if (n3 <= 0) {
        stop(
          "budget is too small for the given fixed stage sizes",
          call. = FALSE
        )
      }
      cv_achieved <- .cv3(n1_eff, n2, n3)
      total_cost <- budget
    } else {
      denom <- n2 * (cv^2 * n1_eff / rel_var - k1 * delta1) - k2 * delta2
      n3 <- k2 * (1 - delta2) / denom
      if (n3 <= 0) {
        stop(
          "target CV is too small for the given fixed stage sizes and parameters",
          call. = FALSE
        )
      }
      total_cost <- fixed_cost +
        C1 * n_psu + C2 * n_psu * n2 + C3 * n_psu * n2 * n3
      cv_achieved <- cv
    }
  }

  n_vec <- c(n_psu = n1, psu_size = n2, ssu_size = n3)
  total_n <- prod(n_vec)

  params <- list(
    stage_cost = c(cost_psu = C1, cost_ssu = C2, cost_tsu = C3),
    delta = c(delta_psu = delta1, delta_ssu = delta2),
    rel_var = rel_var,
    k = c(k_psu = k1, k_ssu = k2),
    resp_rate = resp_rate
  )
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(psu_size)) params$psu_size <- psu_size
  if (!is.null(ssu_size)) params$ssu_size <- ssu_size
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost

  .new_svyplan_cluster(
    n = n_vec,
    stages = 3L,
    total_n = total_n,
    cv = cv_achieved,
    cost = total_cost,
    params = params
  )
}
