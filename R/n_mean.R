#' Sample Size for a Mean
#'
#' Compute the required sample size for estimating a population mean
#' with a specified margin of error or coefficient of variation.
#'
#' @param var Population variance \eqn{S^2} (required).
#' @param mu Population mean. Required when `cv` is specified.
#' @param moe Desired margin of error. Specify exactly one of `moe` or `cv`.
#' @param cv Target coefficient of variation. Specify exactly one of `moe`
#'   or `cv`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (>= 1).
#'
#' @return A `svyplan_n` object.
#'
#' @details
#' Two modes:
#'
#' - **MOE mode**: `n = z^2 * var / (moe^2 + z^2 * var / N)`, then
#'   multiplied by `deff`.
#' - **CV mode**: Computes `CVpop = sqrt(var) / mu`, then
#'   `n = CVpop^2 / (cv^2 + CVpop^2 / N)`, multiplied by `deff`.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @seealso [n_prop()] for proportions, [n_cluster()] for multistage designs.
#'
#' @examples
#' # MOE mode
#' n_mean(var = 100, moe = 2)
#'
#' # CV mode
#' n_mean(var = 100, mu = 50, cv = 0.05)
#'
#' # With FPC and design effect
#' n_mean(var = 100, moe = 2, N = 5000, deff = 1.5)
#'
#' @export
n_mean <- function(var, mu = NULL, moe = NULL, cv = NULL, alpha = 0.05,
                   N = Inf, deff = 1) {
  check_scalar(var, "var")
  check_precision(moe, cv)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)

  if (!is.null(cv) && is.null(mu)) {
    stop("'mu' is required when 'cv' is specified", call. = FALSE)
  }
  if (!is.null(mu)) {
    check_scalar(mu, "mu")
  }

  z <- qnorm(1 - alpha / 2)

  if (!is.null(moe)) {
    # MOE mode: n = z^2 * S^2 / (moe^2 + z^2 * S^2 / N)
    n <- z^2 * var / (moe^2 + z^2 * var / N)
  } else {
    # CV mode: CVpop = sqrt(S^2) / mu, n = CVpop^2 / (cv^2 + CVpop^2 / N)
    CVpop <- sqrt(var) / mu
    n <- CVpop^2 / (cv^2 + CVpop^2 / N)
  }

  n <- n * deff

  params <- list(var = var, alpha = alpha, N = N, deff = deff)
  if (!is.null(mu)) params$mu <- mu
  if (!is.null(moe)) params$moe <- moe else params$cv <- cv

  .new_svyplan_n(
    n      = n,
    type   = "mean",
    params = params
  )
}
