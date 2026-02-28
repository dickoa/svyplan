#' Optimal Multistage Cluster Allocation
#'
#' Compute optimal per-stage sample sizes for a multistage cluster
#' design, minimizing cost for a given precision or minimizing
#' variance for a given budget.
#'
#' @param cost For the default method: numeric vector of per-stage costs.
#'   Length determines the number of stages (2 or 3).
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
#'   means optimize all stages.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The stage-1 sample size is inflated by `1 / resp_rate`.
#' @param fixed_cost Fixed overhead cost (C0). Default 0.
#'   The total cost model becomes
#'   `C = C0 + c1*n_psu + c2*n_psu*psu_size [+ c3*n_psu*psu_size*ssu_size]`.
#'   In budget mode, only `budget - fixed_cost` is available for variable costs;
#'   in CV mode, `fixed_cost` is added to the variable cost.
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
#' Stage count is determined by `length(cost)`. Two dispatch dimensions:
#'
#' - **2-stage vs 3-stage** (vector length)
#' - **budget vs cv** mode (which is non-NULL)
#'
#' When `n_psu` is specified, stage 1 is fixed and only stage 2+ are optimized.
#'
#' If `delta` is a `svyplan_varcomp` object, `delta`, `rel_var`, and `k`
#' are extracted automatically.
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
#' n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
#'
#' # 2-stage, CV mode
#' n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
#'
#' # 2-stage, fixed n_psu
#' n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000, n_psu = 40)
#'
#' # 3-stage
#' n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05), cv = 0.05)
#'
#' # With fixed overhead cost
#' n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000, fixed_cost = 5000)
#'
#' @export
n_cluster <- function(cost, ...) {
  UseMethod("n_cluster")
}

#' @rdname n_cluster
#' @export
n_cluster.default <- function(
  cost,
  delta,
  rel_var = 1,
  k = 1,
  cv = NULL,
  budget = NULL,
  n_psu = NULL,
  resp_rate = 1,
  fixed_cost = 0,
  ...
) {
  check_cost(cost)
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

  stages <- length(cost)
  delta <- .reorder_stage_vec(delta, "delta")
  k <- .reorder_stage_vec(k, "k")
  check_delta(delta, expected_length = stages - 1L)

  if (any(delta == 0) || any(delta == 1)) {
    stop(
      "'delta' must be in (0, 1) for n_cluster(); boundary values make cluster optimization degenerate",
      call. = FALSE
    )
  }

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
  if (!is.null(n_psu)) {
    check_scalar(n_psu, "n_psu")
  }
  check_fixed_cost(fixed_cost, budget)

  if (stages == 2L) {
    .n_cluster_2stage(
      cost,
      delta,
      rel_var,
      k,
      cv,
      budget,
      n_psu,
      resp_rate,
      fixed_cost
    )
  } else {
    k <- rep_len(k, 2L)
    .n_cluster_3stage(
      cost,
      delta,
      rel_var,
      k,
      cv,
      budget,
      n_psu,
      resp_rate,
      fixed_cost
    )
  }
}

#' @rdname n_cluster
#' @export
n_cluster.svyplan_prec <- function(cost, cv = NULL, budget = NULL, ...) {
  x <- cost
  if (x$type != "cluster") {
    stop("n_cluster requires a svyplan_prec of type 'cluster'", call. = FALSE)
  }
  p <- x$params
  if (is.null(cv) && is.null(budget)) {
    cv <- x$cv
  }
  n_cluster.default(
    cost = p$cost,
    delta = p$delta,
    rel_var = p$rel_var,
    k = p$k,
    cv = cv,
    budget = budget,
    n_psu = p$n_psu,
    resp_rate = p$resp_rate %||% 1,
    fixed_cost = p$fixed_cost %||% 0
  )
}

#' @keywords internal
#' @noRd
.n_cluster_2stage <- function(
  cost,
  delta,
  rel_var,
  k,
  cv,
  budget,
  n_psu,
  resp_rate,
  fixed_cost = 0
) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  var_budget <- if (!is.null(budget)) budget - fixed_cost else NULL

  if (is.null(n_psu)) {
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
  } else {
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
  }

  n_vec <- c(n_psu = n1_opt, psu_size = n2_opt)
  total_n <- prod(n_vec)

  params <- list(
    cost = cost,
    delta = delta,
    rel_var = rel_var,
    k = k,
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
  cost,
  delta,
  rel_var,
  k,
  cv,
  budget,
  n_psu,
  resp_rate,
  fixed_cost = 0
) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  C3 <- cost[3L]
  delta1 <- delta[1L]
  delta2 <- delta[2L]
  k1 <- k[1L]
  k2 <- k[2L]
  var_budget <- if (!is.null(budget)) budget - fixed_cost else NULL

  n3_opt <- sqrt((1 - delta2) / delta2 * C2 / C3)

  if (is.null(n_psu)) {
    n2_opt <- 1 / n3_opt * sqrt((1 - delta2) / delta1 * C1 / C3 * k2 / k1)

    if (!is.null(budget)) {
      n1_opt <- var_budget / (C1 + C2 * n2_opt + C3 * n2_opt * n3_opt)
      n1_eff <- n1_opt * resp_rate
      cv_achieved <- sqrt(
        rel_var /
          (n1_eff * n2_opt * n3_opt) *
          (k1 * delta1 * n2_opt * n3_opt + k2 * (1 + delta2 * (n3_opt - 1)))
      )
      total_cost <- budget
    } else {
      n1_eff_needed <- rel_var /
        (cv^2 * n2_opt * n3_opt) *
        (k1 * delta1 * n2_opt * n3_opt + k2 * (1 + delta2 * (n3_opt - 1)))
      n1_opt <- n1_eff_needed / resp_rate
      total_cost <- fixed_cost +
        C1 * n1_opt +
        C2 * n1_opt * n2_opt +
        C3 * n1_opt * n2_opt * n3_opt
      cv_achieved <- cv
    }
  } else {
    n1_opt <- n_psu
    n1_eff <- n_psu * resp_rate
    if (!is.null(budget)) {
      C_prime <- var_budget / n_psu - C1
      n2_opt <- C_prime / (C2 + C3 * n3_opt)
      if (n2_opt <= 0) {
        stop(
          "budget is too small for the given fixed stage-1 size",
          call. = FALSE
        )
      }
      cv_achieved <- sqrt(
        rel_var /
          (n1_eff * n2_opt * n3_opt) *
          (k1 * delta1 * n2_opt * n3_opt + k2 * (1 + delta2 * (n3_opt - 1)))
      )
      total_cost <- C1 *
        n_psu +
        C2 * n_psu * n2_opt +
        C3 * n_psu * n2_opt * n3_opt +
        fixed_cost
    } else {
      n2_opt <- k2 *
        (1 + delta2 * (n3_opt - 1)) /
        (n3_opt * (cv^2 * n1_eff / rel_var - k1 * delta1))
      if (n2_opt <= 0) {
        stop(
          "target CV is too small for the given fixed stage-1 size and parameters",
          call. = FALSE
        )
      }
      total_cost <- fixed_cost +
        C1 * n_psu +
        C2 * n_psu * n2_opt +
        C3 * n_psu * n2_opt * n3_opt
      cv_achieved <- cv
    }
  }

  n_vec <- c(n_psu = n1_opt, psu_size = n2_opt, ssu_size = n3_opt)
  total_n <- prod(n_vec)

  params <- list(
    cost = cost,
    delta = delta,
    rel_var = rel_var,
    k = k,
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
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 3L,
    total_n = total_n,
    cv = cv_achieved,
    cost = total_cost,
    params = params
  )
}
