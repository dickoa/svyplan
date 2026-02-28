#' Sample Size for a Mean
#'
#' Compute the required sample size for estimating a population mean
#' with a specified margin of error or coefficient of variation.
#'
#' @param var For the default method: population variance \eqn{S^2}.
#'   For `svyplan_prec` objects: a precision result from [prec_mean()].
#' @param ... Additional arguments passed to methods.
#' @param mu Population mean magnitude (positive). Required when `cv` is specified.
#' @param moe Desired margin of error. Specify exactly one of `moe` or `cv`.
#' @param cv Target coefficient of variation. Specify exactly one of `moe`
#'   or `cv`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The sample size is inflated by `1 / resp_rate`.
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
#' The FPC is the standard Cochran (1977) one-step form. Unlike [n_prop()],
#' no `N/(N-1)` adjustment is needed because `var` is already defined on
#' `N-1` degrees of freedom.
#'
#' All methods use the normal (z) quantile. This is standard for survey
#' sampling where the sample size is large enough for the CLT to apply.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @seealso [n_prop()] for proportions, [n_cluster()] for multistage designs,
#'   [n_multi()] for multiple indicators, [prec_mean()] for the inverse.
#'
#' @examples
#' # MOE mode
#' n_mean(var = 100, moe = 2)
#'
#' # CV mode
#' n_mean(var = 100, mu = 50, cv = 0.05)
#'
#' # With FPC, design effect, and response rate
#' n_mean(var = 100, moe = 2, N = 5000, deff = 1.5, resp_rate = 0.8)
#'
#' @export
n_mean <- function(var, ...) UseMethod("n_mean")

#' @rdname n_mean
#' @export
n_mean.default <- function(var, mu = NULL, moe = NULL, cv = NULL, alpha = 0.05,
                           N = Inf, deff = 1, resp_rate = 1, ...) {
  check_scalar(var, "var")
  check_precision(moe, cv)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)

  if (!is.null(cv) && is.null(mu)) {
    stop("'mu' is required when 'cv' is specified", call. = FALSE)
  }
  if (!is.null(mu)) {
    check_scalar(mu, "mu")
  }

  z <- qnorm(1 - alpha / 2)

  if (!is.null(moe)) {
    n <- z^2 * var / (moe^2 + z^2 * var / N)
  } else {
    CVpop <- sqrt(var) / mu
    n <- CVpop^2 / (cv^2 + CVpop^2 / N)
  }

  n <- n * deff
  n <- .apply_resp_rate(n, resp_rate)

  if (!is.infinite(N) && n > N)
    warning("Calculated sample size exceeds population size N.", call. = FALSE)

  params <- list(var = var, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate)
  if (!is.null(mu)) params$mu <- mu
  if (!is.null(moe)) params$moe <- moe else params$cv <- cv

  .new_svyplan_n(
    n      = n,
    type   = "mean",
    params = params
  )
}

#' @rdname n_mean
#' @export
n_mean.svyplan_prec <- function(var, moe = NULL, cv = NULL, ...) {
  x <- var
  if (x$type != "mean") {
    stop("n_mean requires a svyplan_prec of type 'mean'", call. = FALSE)
  }
  par <- x$params
  if (is.null(moe) && is.null(cv)) {
    moe <- x$moe
  }
  n_mean.default(
    var       = par$var,
    mu        = par$mu,
    moe       = moe,
    cv        = cv,
    alpha     = par$alpha,
    N         = par$N,
    deff      = par$deff,
    resp_rate = par$resp_rate
  )
}
