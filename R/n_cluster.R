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
#' @param m Fixed stage-1 sample size. `NULL` (default) means optimize all
#'   stages.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The stage-1 sample size is inflated by `1 / resp_rate`.
#'
#' @return A `svyplan_cluster` object with components:
#' \describe{
#'   \item{`n`}{Named numeric vector of per-stage sample sizes
#'     (e.g. `c(n1 = 80, n2 = 12)`).}
#'   \item{`stages`}{Number of stages (2 or 3).}
#'   \item{`total_n`}{Total sample size (`prod(n)`).}
#'   \item{`cv`}{Achieved coefficient of variation.}
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
#' When `m` is specified, stage 1 is fixed and only stage 2+ are optimized.
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
#' # 2-stage, fixed m
#' n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000, m = 40)
#'
#' # 3-stage
#' n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05), cv = 0.05)
#'
#' @export
n_cluster <- function(cost, ...) UseMethod("n_cluster")

#' @rdname n_cluster
#' @export
n_cluster.default <- function(cost, delta, rel_var = 1, k = 1,
                              cv = NULL, budget = NULL, m = NULL,
                              resp_rate = 1, ...) {
  check_cost(cost)
  check_resp_rate(resp_rate)

  if (inherits(delta, "svyplan_varcomp")) {
    vc <- delta
    delta <- vc$delta
    rel_var <- vc$rel_var
    k <- vc$k
  }

  check_scalar(rel_var, "rel_var")
  if (!is.numeric(k) || length(k) == 0L || anyNA(k) || any(k <= 0))
    stop("'k' must contain positive values", call. = FALSE)

  stages <- length(cost)
  check_delta(delta, expected_length = stages - 1L)

  if (any(delta == 0) || any(delta == 1))
    stop("'delta' must be in (0, 1) for n_cluster(); boundary values make cluster optimization degenerate",
         call. = FALSE)

  has_cv <- !is.null(cv)
  has_budget <- !is.null(budget)
  if (has_cv == has_budget) {
    stop("specify exactly one of 'cv' or 'budget'", call. = FALSE)
  }
  if (has_cv) check_scalar(cv, "cv")
  if (has_budget) check_scalar(budget, "budget")
  if (!is.null(m)) check_scalar(m, "m")

  if (stages == 2L) {
    .n_cluster_2stage(cost, delta, rel_var, k, cv, budget, m, resp_rate)
  } else {
    k <- rep_len(k, 2L)
    .n_cluster_3stage(cost, delta, rel_var, k, cv, budget, m, resp_rate)
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
    cost     = p$cost,
    delta    = p$delta,
    rel_var  = p$rel_var,
    k        = p$k,
    cv       = cv,
    budget   = budget,
    m        = p$m,
    resp_rate = p$resp_rate %||% 1
  )
}

#' @keywords internal
#' @noRd
.n_cluster_2stage <- function(cost, delta, rel_var, k, cv, budget, m,
                              resp_rate) {
  C1 <- cost[1L]
  C2 <- cost[2L]

  if (is.null(m)) {
    n2_opt <- sqrt(C1 / C2 * (1 - delta) / delta)

    if (!is.null(budget)) {
      n1_opt <- budget / (C1 + C2 * n2_opt)
      n1_eff <- n1_opt * resp_rate
      cv_achieved <- sqrt(rel_var / (n1_eff * n2_opt) * k *
                            (1 + delta * (n2_opt - 1)))
      total_cost <- budget
    } else {
      n1_eff_needed <- rel_var * k * (1 + delta * (n2_opt - 1)) / (n2_opt * cv^2)
      n1_opt <- n1_eff_needed / resp_rate
      total_cost <- C1 * n1_opt + C2 * n1_opt * n2_opt
      cv_achieved <- cv
    }
  } else {
    n1_opt <- m
    n1_eff <- m * resp_rate
    if (!is.null(budget)) {
      n2_opt <- (budget - C1 * m) / (C2 * m)
      if (n2_opt <= 0) {
        stop("budget is too small for the given fixed stage-1 size",
             call. = FALSE)
      }
      cv_achieved <- sqrt(rel_var * k / (n1_eff * n2_opt) *
                            (1 + delta * (n2_opt - 1)))
      total_cost <- budget
    } else {
      n2_opt <- (1 - delta) / (cv^2 * n1_eff / (rel_var * k) - delta)
      if (n2_opt <= 0) {
        stop("target CV is too small for the given fixed stage-1 size and parameters", call. = FALSE)
      }
      total_cost <- C1 * m + C2 * m * n2_opt
      cv_achieved <- cv
    }
  }

  n_vec <- c(n1 = n1_opt, n2 = n2_opt)
  total_n <- prod(n_vec)

  params <- list(cost = cost, delta = delta, rel_var = rel_var, k = k,
                 resp_rate = resp_rate)
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(m)) params$m <- m

  .new_svyplan_cluster(
    n       = n_vec,
    stages  = 2L,
    total_n = total_n,
    cv      = cv_achieved,
    cost    = total_cost,
    params  = params
  )
}

#' @keywords internal
#' @noRd
.n_cluster_3stage <- function(cost, delta, rel_var, k, cv, budget, m,
                              resp_rate) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  C3 <- cost[3L]
  delta1 <- delta[1L]
  delta2 <- delta[2L]
  k1 <- k[1L]
  k2 <- k[2L]

  n3_opt <- sqrt((1 - delta2) / delta2 * C2 / C3)

  if (is.null(m)) {
    n2_opt <- 1 / n3_opt * sqrt((1 - delta2) / delta1 * C1 / C3 * k2 / k1)

    if (!is.null(budget)) {
      n1_opt <- budget / (C1 + C2 * n2_opt + C3 * n2_opt * n3_opt)
      n1_eff <- n1_opt * resp_rate
      cv_achieved <- sqrt(rel_var / (n1_eff * n2_opt * n3_opt) *
                            (k1 * delta1 * n2_opt * n3_opt +
                               k2 * (1 + delta2 * (n3_opt - 1))))
      total_cost <- budget
    } else {
      n1_eff_needed <- rel_var / (cv^2 * n2_opt * n3_opt) *
        (k1 * delta1 * n2_opt * n3_opt + k2 * (1 + delta2 * (n3_opt - 1)))
      n1_opt <- n1_eff_needed / resp_rate
      total_cost <- C1 * n1_opt + C2 * n1_opt * n2_opt +
        C3 * n1_opt * n2_opt * n3_opt
      cv_achieved <- cv
    }
  } else {
    n1_opt <- m
    n1_eff <- m * resp_rate
    if (!is.null(budget)) {
      C_prime <- budget / m - C1
      n2_opt <- C_prime / (C2 + C3 * n3_opt)
      if (n2_opt <= 0) {
        stop("budget is too small for the given fixed stage-1 size",
             call. = FALSE)
      }
      cv_achieved <- sqrt(rel_var / (n1_eff * n2_opt * n3_opt) *
                            (k1 * delta1 * n2_opt * n3_opt +
                               k2 * (1 + delta2 * (n3_opt - 1))))
      total_cost <- C1 * m + C2 * m * n2_opt + C3 * m * n2_opt * n3_opt
    } else {
      n2_opt <- k2 * (1 + delta2 * (n3_opt - 1)) /
        (n3_opt * (cv^2 * n1_eff / rel_var - k1 * delta1))
      if (n2_opt <= 0) {
        stop("target CV is too small for the given fixed stage-1 size and parameters", call. = FALSE)
      }
      total_cost <- C1 * m + C2 * m * n2_opt + C3 * m * n2_opt * n3_opt
      cv_achieved <- cv
    }
  }

  n_vec <- c(n1 = n1_opt, n2 = n2_opt, n3 = n3_opt)
  total_n <- prod(n_vec)

  params <- list(cost = cost, delta = delta, rel_var = rel_var, k = k,
                 resp_rate = resp_rate)
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(m)) params$m <- m

  .new_svyplan_cluster(
    n       = n_vec,
    stages  = 3L,
    total_n = total_n,
    cv      = cv_achieved,
    cost    = total_cost,
    params  = params
  )
}
