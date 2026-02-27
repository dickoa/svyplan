#' Check that a value is a positive scalar
#' @keywords internal
#' @noRd
check_scalar <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop(sprintf("'%s' must be a numeric scalar", name), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("'%s' must not be NA", name), call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(sprintf("'%s' must be positive", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that a value is in (0, 1)
#' @keywords internal
#' @noRd
check_proportion <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop(sprintf("'%s' must be a numeric scalar", name), call. = FALSE)
  }
  if (anyNA(x) || x <= 0 || x >= 1) {
    stop(sprintf("'%s' must be in (0, 1)", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that exactly one of moe/cv is specified
#' @keywords internal
#' @noRd
check_precision <- function(moe, cv) {
  has_moe <- !is.null(moe)
  has_cv <- !is.null(cv)
  if (has_moe == has_cv) {
    stop("specify exactly one of 'moe' or 'cv'", call. = FALSE)
  }
  if (has_moe) {
    check_scalar(moe, "moe")
  }
  if (has_cv) {
    check_scalar(cv, "cv")
  }
  invisible(TRUE)
}

#' Check alpha in (0, 1)
#' @keywords internal
#' @noRd
check_alpha <- function(alpha) {
  check_proportion(alpha, "alpha")
}

#' Check deff > 0
#' @keywords internal
#' @noRd
check_deff <- function(deff) {
  check_scalar(deff, "deff")
  invisible(TRUE)
}

#' Check population size N > 1 or Inf
#' @keywords internal
#' @noRd
check_population_size <- function(N) {
  if (!is.numeric(N) || length(N) != 1L || anyNA(N) || N <= 1) {
    stop("'N' must be greater than 1 (or Inf)", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check weights vector: numeric, positive, non-empty
#' @keywords internal
#' @noRd
check_weights <- function(w, name = "x") {
  if (!is.numeric(w) || length(w) == 0L) {
    stop(sprintf("'%s' must be a non-empty numeric vector", name), call. = FALSE)
  }
  if (anyNA(w)) {
    stop(sprintf("'%s' must not contain NA values", name), call. = FALSE)
  }
  if (any(w <= 0)) {
    stop(sprintf("'%s' must contain only positive values", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check cost vector
#' @keywords internal
#' @noRd
check_cost <- function(cost) {
  if (!is.numeric(cost) || length(cost) < 2L) {
    stop("'cost' must be a numeric vector of length >= 2", call. = FALSE)
  }
  if (anyNA(cost) || any(cost <= 0)) {
    stop("'cost' must contain positive values only", call. = FALSE)
  }
  if (length(cost) > 3L) {
    stop("4+ stage optimization is not yet supported", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check fixed cost (C0)
#' @keywords internal
#' @noRd
check_fixed_cost <- function(fixed_cost, budget = NULL) {
  if (!is.numeric(fixed_cost) || length(fixed_cost) != 1L ||
      is.na(fixed_cost) || fixed_cost < 0) {
    stop("'fixed_cost' must be a non-negative numeric scalar", call. = FALSE)
  }
  if (!is.null(budget) && fixed_cost >= budget) {
    stop("'fixed_cost' must be less than 'budget'", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check delta vector
#' @keywords internal
#' @noRd
check_delta <- function(delta, expected_length = NULL) {
  if (!is.numeric(delta) || length(delta) == 0L) {
    stop("'delta' must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(delta) || any(delta < 0) || any(delta > 1)) {
    stop("'delta' values must be in [0, 1]", call. = FALSE)
  }
  if (!is.null(expected_length) && length(delta) != expected_length) {
    stop(sprintf("'delta' must have length %d (stages - 1)", expected_length),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check sides parameter (1 or 2)
#' @keywords internal
#' @noRd
check_sides <- function(sides) {
  if (!identical(sides, 1L) && !identical(sides, 1) &&
      !identical(sides, 2L) && !identical(sides, 2)) {
    stop("'sides' must be 1 or 2", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check overlap fraction in \[0, 1\]
#' @keywords internal
#' @noRd
check_overlap <- function(overlap) {
  if (!is.numeric(overlap) || length(overlap) != 1L ||
      anyNA(overlap) || overlap < 0 || overlap > 1) {
    stop("'overlap' must be a number in [0, 1]", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check correlation coefficient in \[0, 1\]
#' @keywords internal
#' @noRd
check_rho <- function(rho) {
  if (!is.numeric(rho) || length(rho) != 1L ||
      anyNA(rho) || rho < 0 || rho > 1) {
    stop("'rho' must be a number in [0, 1]", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check response rate in (0, 1]
#' @keywords internal
#' @noRd
check_resp_rate <- function(resp_rate) {
  if (!is.numeric(resp_rate) || length(resp_rate) != 1L ||
      anyNA(resp_rate) || resp_rate <= 0 || resp_rate > 1) {
    stop("'resp_rate' must be a number in (0, 1]", call. = FALSE)
  }
  invisible(TRUE)
}

#' Apply finite population correction
#' @keywords internal
#' @noRd
.apply_fpc <- function(n0, N) {
  if (is.infinite(N)) n0 else n0 / (1 + n0 / N)
}

#' Apply response rate adjustment (inflate n)
#' @keywords internal
#' @noRd
.apply_resp_rate <- function(n, resp_rate) {
  n / resp_rate
}

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Clamp FPC to 0 when n_eff >= N (census)
#' @keywords internal
#' @noRd
.clamp_fpc <- function(fpc, n_eff, N) {
  if (is.infinite(N)) return(fpc)
  if (n_eff >= N) {
    warning("effective sample size (", round(n_eff, 1),
            ") >= population size (", N,
            "); FPC set to 0 (census)", call. = FALSE)
    return(0)
  }
  fpc
}

#' Validate common optional columns in multi-indicator targets
#' @keywords internal
#' @noRd
.validate_common_columns <- function(targets) {
  if ("alpha" %in% names(targets)) {
    vals <- targets$alpha[!is.na(targets$alpha)]
    if (any(vals <= 0 | vals >= 1))
      stop("'alpha' values must be in (0, 1)", call. = FALSE)
  }
  if ("deff" %in% names(targets)) {
    vals <- targets$deff[!is.na(targets$deff)]
    if (any(vals <= 0))
      stop("'deff' values must be positive", call. = FALSE)
  }
  if ("N" %in% names(targets)) {
    vals <- targets$N[!is.na(targets$N)]
    if (any(vals <= 1))
      stop("'N' values must be greater than 1 (or Inf)", call. = FALSE)
  }
  if ("resp_rate" %in% names(targets)) {
    vals <- targets$resp_rate[!is.na(targets$resp_rate)]
    if (any(vals <= 0 | vals > 1))
      stop("'resp_rate' values must be in (0, 1]", call. = FALSE)
  }
  if ("n" %in% names(targets)) {
    if (anyNA(targets$n) || any(targets$n <= 0))
      stop("'n' values must be positive and non-NA", call. = FALSE)
  }
  invisible(TRUE)
}

#' Weighted variance
#' @keywords internal
#' @noRd
.wtdvar <- function(x, w) {
  n <- length(w)
  sw <- sum(w)
  xbarw <- sum(w * x) / sw
  n / (n - 1) * sum(w * (x - xbarw)^2) / sw
}
