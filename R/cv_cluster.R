#' CV for a Multistage Cluster Allocation
#'
#' Compute the coefficient of variation for a given multistage sample
#' allocation. This is the inverse of [n_cluster()].
#'
#' @param n Numeric vector of per-stage sample sizes: `c(n1, n2)` for
#'   2-stage or `c(n1, n2, n3)` for 3-stage.
#' @param delta Numeric vector of homogeneity measures (length = stages - 1),
#'   or a `svyplan_varcomp` object.
#' @param rel_var Unit relvariance (default 1).
#' @param k Ratio parameter(s). Scalar for 2-stage, length-2 vector for
#'   3-stage (default 1).
#'
#' @return A numeric scalar: the coefficient of variation.
#'
#' @details
#' Stage count is determined by `length(n)`.
#'
#' **2-stage** (Valliant et al., 2018, Eq. 9.2.23):
#' \deqn{CV = \sqrt{\frac{V \cdot k}{n_1 \cdot n_2} (1 + \delta (n_2 - 1))}}
#'
#' **3-stage**:
#' \deqn{CV = \sqrt{\frac{V}{n_1 \cdot n_2 \cdot n_3} (k_1 \delta_1 n_2 n_3
#'   + k_2 (1 + \delta_2 (n_3 - 1)))}}
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 9.
#'
#' @seealso [n_cluster()] for the inverse operation.
#'
#' @examples
#' cv_cluster(n = c(50, 12), delta = 0.05)
#' cv_cluster(n = c(50, 12, 8), delta = c(0.01, 0.05))
#'
#' @export
cv_cluster <- function(n, delta, rel_var = 1, k = 1) {
  if (inherits(delta, "svyplan_varcomp")) {
    vc <- delta
    delta <- vc$delta
    rel_var <- vc$rel_var
    k <- vc$k
  }

  if (!is.numeric(n) || length(n) < 2L) {
    stop("'n' must be a numeric vector of length >= 2", call. = FALSE)
  }
  if (length(n) > 3L) {
    stop("4+ stage CV calculation is not yet supported", call. = FALSE)
  }

  stages <- length(n)
  check_delta(delta, expected_length = stages - 1L)

  if (any(n <= 0)) {
    stop("all elements of 'n' must be positive", call. = FALSE)
  }

  if (stages == 2L) {
    .cv_cluster_2stage(n, delta, rel_var, k)
  } else {
    k <- rep_len(k, 2L)
    .cv_cluster_3stage(n, delta, rel_var, k)
  }
}

#' @keywords internal
#' @noRd
.cv_cluster_2stage <- function(n, delta, rel_var, k) {
  n1 <- n[1L]
  n2 <- n[2L]
  sqrt(rel_var / (n1 * n2) * k * (1 + delta * (n2 - 1)))
}

#' @keywords internal
#' @noRd
.cv_cluster_3stage <- function(n, delta, rel_var, k) {
  n1 <- n[1L]
  n2 <- n[2L]
  n3 <- n[3L]
  delta1 <- delta[1L]
  delta2 <- delta[2L]
  k1 <- k[1L]
  k2 <- k[2L]
  sqrt(rel_var / (n1 * n2 * n3) * (k1 * delta1 * n2 * n3 +
                                     k2 * (1 + delta2 * (n3 - 1))))
}
