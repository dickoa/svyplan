#' Effective Sample Size
#'
#' Compute the effective sample size, adjusting for design effects.
#' This is an S3 generic that mirrors [design_effect()].
#'
#' @param x A numeric vector of survey weights (for diagnostic methods),
#'   or `NULL` (for the `"cluster"` planning method).
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric scalar.
#'
#' @seealso [design_effect()] for the underlying design effect computation.
#'
#' @examples
#' # Effective sample size from weights (Kish)
#' w <- runif(100, 1, 5)
#' effective_n(w, method = "kish")
#'
#' # Planning: effective n given cluster design
#' effective_n(delta = 0.05, m = 25, n = 800, method = "cluster")
#'
#' @name effective_n
#' @export
effective_n <- function(x = NULL, ...) {
  if (is.null(x)) {
    effective_n.default(x, ...)
  } else {
    UseMethod("effective_n")
  }
}

#' @describeIn effective_n Method for numeric weights vector.
#'
#' @param y Outcome variable (for `"henry"`, `"spencer"`, `"cr"`).
#' @param x_cal Calibration covariate (for `"henry"`).
#' @param p 1-draw selection probabilities (for `"spencer"`).
#' @param strvar Stratum IDs (for `"cr"`).
#' @param clvar Cluster IDs (for `"cr"`).
#' @param stages Integer vector of stages per stratum (for `"cr"`).
#' @param method One of `"kish"`, `"henry"`, `"spencer"`, or `"cr"`.
#'
#' @export
effective_n.numeric <- function(x, ..., y = NULL, x_cal = NULL, p = NULL,
                                strvar = NULL, clvar = NULL, stages = NULL,
                                method = "kish") {
  method <- match.arg(method, c("kish", "henry", "spencer", "cr"))
  check_weights(x)

  if (method == "kish") {
    # Direct formula: sum(w)^2 / sum(w^2)
    sum(x)^2 / sum(x^2)
  } else {
    deff <- design_effect.numeric(x, y = y, x_cal = x_cal, p = p,
                                  strvar = strvar, clvar = clvar,
                                  stages = stages, method = method)
    n <- length(x)
    if (is.list(deff)) {
      n / deff$overall
    } else {
      n / deff
    }
  }
}

#' @describeIn effective_n Planning method (no weights needed).
#'
#' @param delta ICC / homogeneity measure, scalar or `svyplan_varcomp`.
#' @param m Mean cluster size (scalar).
#' @param n Total sample size (required for the cluster method).
#'
#' @export
effective_n.default <- function(x = NULL, ..., delta = NULL, m = NULL,
                                n = NULL, method = "cluster") {
  method <- match.arg(method, "cluster")

  if (is.null(n)) {
    stop("'n' is required for the cluster method", call. = FALSE)
  }
  check_scalar(n, "n")

  deff <- design_effect.default(delta = delta, m = m, method = method)
  n / deff
}
