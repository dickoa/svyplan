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

#' Check deff >= 1
#' @keywords internal
#' @noRd
check_deff <- function(deff) {
  check_scalar(deff, "deff")
  if (deff < 1) {
    stop("'deff' must be >= 1", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check population size N > 0 or Inf
#' @keywords internal
#' @noRd
check_population_size <- function(N) {
  if (!is.numeric(N) || length(N) != 1L || anyNA(N) || N <= 0) {
    stop("'N' must be a positive number or Inf", call. = FALSE)
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
    warning(sprintf("some values in '%s' are <= 0", name), call. = FALSE)
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

#' Apply finite population correction
#' @keywords internal
#' @noRd
.apply_fpc <- function(n0, N) {
  if (is.infinite(N)) n0 else n0 / (1 + n0 / N)
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
