#' Sample Size for a Mean
#'
#' Compute the required sample size for estimating a population mean
#' with a specified margin of error or coefficient of variation.
#'
#' @param var For the default method: population variance \eqn{S^2}.
#'   Estimate from a pilot study, a previous survey, or published data
#'   for a similar population. When uncertain, use a conservative
#'   (larger) estimate to avoid under-sizing.
#'   For `svyplan_prec` objects: a precision result from [prec_mean()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param sd Population standard deviation, an alternative spelling of
#'   `var`. Supply exactly one of `var` or `sd`. Stratum frames and
#'   published survey reports usually quote standard deviations.
#' @param mu Population mean. Required when `cv` is
#'   specified, because CV is defined as SE / mean.
#' @param moe Desired margin of error, the half-width of the confidence
#'   interval, in the same units as the variable. For example, if
#'   measuring income in dollars, `moe = 50` means the 95 percent CI
#'   should be no wider than +/- $50. Specify exactly one of `moe`
#'   or `cv`.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.05` means the standard error should be at
#'   most 5 percent of the estimate. Use `cv` when you want precision
#'   to scale with the estimate (common in economic surveys). Use `moe`
#'   when you want a fixed absolute precision. Requires `mu`. Specify
#'   exactly one of `moe` or `cv`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction. Setting a finite `N` reduces the required sample size
#'   when the sampling fraction is non-negligible (rule of thumb:
#'   matters when n/N > 5 percent).
#' @param deff Design effect multiplier (> 0). See [n_prop()] for
#'   guidance on estimating DEFF. Values < 1 are valid for efficient
#'   designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The required sample size is inflated by `1 / resp_rate`.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "mean"`:
#' \describe{
#'   \item{`n`}{Required sample size, continuous and gross. It already
#'     carries `deff` and the `1 / resp_rate` inflation, so it counts the
#'     units to release, not the completed interviews. `$n` and
#'     `as.double()` keep the unrounded value, which is what makes the
#'     round trip through [prec_mean()] exact; `print()` and
#'     `as.integer()` round it up to the whole units you would field.
#'     Take the field figure from `as.integer()` rather than from `$n`.}
#'   \item{`se`, `moe`, `cv`}{Precision the design achieves at that `n`,
#'     the same values [prec_mean()] reports for the same inputs. `moe`
#'     is `qnorm(1 - alpha / 2) * se`, and the interval is symmetric
#'     about the mean. `cv` is `NA` unless `mu` was supplied, since a
#'     relative standard error needs a mean to be relative to.}
#'   \item{`params`}{The validated inputs (`var`, `alpha`, `N`, `deff`,
#'     `resp_rate`, `mu` when given, and whichever of `moe` or `cv` was
#'     the target). Dispersion is always stored as `var`, including when
#'     you supplied `sd`. [predict()], [confint()] and the [prec_mean()]
#'     round trip read the design back from here.}
#' }
#'
#' @details
#' Two modes:
#'
#' - **MOE mode**: `n = deff * z^2 * var / (moe^2 + deff * z^2 * var / N)`.
#' - **CV mode**: Computes `CVpop = sqrt(var) / abs(mu)`, then
#'   `n = deff * CVpop^2 / (cv^2 + deff * CVpop^2 / N)`.
#'
#' `deff` appears in the denominator as well as the numerator: it inflates
#' the variance the finite population correction is then applied to, rather
#' than scaling a size already corrected. The two coincide only at infinite
#' `N`. For `var = 100`, `moe = 2`, `N = 100` and `deff = 2`, inflating
#' afterwards would give 97.98 against the correct 65.76.
#'
#' ## Finite population correction
#'
#' Setting `N` to a finite value reduces the required sample size when
#' the sampling fraction (n/N) is non-negligible (rule of thumb: matters
#' when n/N > 5 percent). Unlike [n_prop()], no `N/(N-1)` adjustment is
#' needed because `var` is already defined on `N-1` degrees of freedom.
#' See [n_prop()] for a fuller explanation of FPC.
#'
#' All methods use the normal (z) quantile. This is standard for survey
#' sampling where the sample size is large enough for the CLT to apply.
#'
#' ## Sample size for a total
#'
#' A separate `n_total()` function is not needed because the sample size
#' for a population total \eqn{\hat{Y} = N \bar{y}} is identical to
#' the sample size for the mean. The two are related by a factor of
#' \eqn{N}:
#'
#' - **CV mode**: \eqn{CV(\hat{Y}) = CV(\bar{y})}, so the required
#'   sample size is the same. Use `n_mean(var, mu, cv)` directly.
#' - **MOE mode**: \eqn{MOE(\hat{Y}) = N \times MOE(\bar{y})}, so
#'   divide the target MOE for the total by \eqn{N}:
#'   `n_mean(var, moe = moe_total / N, N = N)`.
#'
#' See Examples below.
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
#' # --- Sample size for a total ---
#' # Target: estimate total income (N = 10000) with MOE of 500000
#' n_mean(var = 2500, moe = 500000 / 10000, N = 10000)
#'
#' # CV mode: identical for means and totals
#' n_mean(var = 2500, mu = 300, cv = 0.05, N = 10000)
#'
#' @export
n_mean <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", n_mean.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_mean")
}

#' @rdname n_mean
#' @export
n_mean.default <- function(
  var = NULL,
  ...,
  sd = NULL,
  mu = NULL,
  moe = NULL,
  cv = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_mean.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_mean.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  var <- .resolve_var(var, sd)
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
    check_mu(mu)
  }

  z <- qnorm(1 - alpha / 2)

  if (!is.null(moe)) {
    n <- z^2 * deff * var / (moe^2 + z^2 * deff * var / N)
  } else {
    CVpop2 <- var / mu^2
    n <- deff * CVpop2 / (cv^2 + deff * CVpop2 / N)
  }

  n <- .apply_resp_rate(n, resp_rate)
  .check_attainable(n, N, resp_rate)

  params <- list(
    var = var,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate
  )

  if (!is.null(mu)) {
    params$mu <- mu
  }
  if (!is.null(moe)) {
    params$moe <- moe
  } else {
    params$cv <- cv
  }

  .new_svyplan_n(
    n = n,
    type = "mean",
    params = params
  )
}

#' @rdname n_mean
#' @export
n_mean.svyplan_prec <- function(var, ..., moe = NULL, cv = NULL) {
  x <- var
  if (x$type != "mean") {
    stop("n_mean requires a svyplan_prec of type 'mean'", call. = FALSE)
  }
  par <- x$params
  if (is.null(moe) && is.null(cv)) {
    moe <- x$moe
  }
  args <- list(
    var = par$var,
    mu = par$mu,
    moe = moe,
    cv = cv,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate
  )
  do.call(n_mean.default, .roundtrip_args(args, list(...), n_mean.default))
}
