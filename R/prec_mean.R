#' Sampling precision for a mean
#'
#' Compute the sampling error (SE, margin of error, CV) for estimating a
#' population mean given a sample size. This is the inverse of [n_mean()].
#'
#' @param var For the default method: population variance \eqn{S^2}.
#'   For `svyplan_n` objects: a sample size result from [n_mean()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param sd Population standard deviation, an alternative spelling of
#'   `var`. Supply exactly one of `var` or `sd`. Stratum frames and
#'   published survey reports usually quote standard deviations.
#' @param n Sample size, measured as gross units drawn and bounded by a
#'   finite `N`.
#' @param mu Population mean. Required for the CV component.
#' @param cv Target relative standard error, supplied instead of `mu` to solve
#'   for the smallest mean the design measures that precisely. Supply at most
#'   one of `mu` or `cv`. Omitting both leaves `cv` undefined, as before.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The effective sample size is `n * resp_rate`.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()]. It switches the interval quantile from normal to t.
#'   `NULL` (default) applies no adjustment. See [n_prop()] for the full
#'   account.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "mean"`:
#' \describe{
#'   \item{`se`}{Standard error of the planned estimate, computed on the
#'     net sample `n * resp_rate`.}
#'   \item{`moe`}{Margin of error, `q * se`, with `q` defined in Details. The
#'     interval is symmetric, so the limits are `mu - moe` and
#'     `mu + moe`.}
#'   \item{`cv`}{Relative standard error, `se / abs(mu)`. `NA` when `mu` is
#'     not supplied, since a relative standard error needs a mean to be
#'     relative to.}
#'   \item{`solved`}{`"mu"` when the mean was solved for, and absent
#'     otherwise. The solved value is in `params$mu`.}
#'   \item{`params`}{The validated inputs (`var`, `n`, `alpha`, `N`,
#'     `deff`, `resp_rate`, `df`, and `mu` when given or solved for). Dispersion
#'     is always stored as `var`, including when you supplied `sd`.
#'     [predict()], [confint()] and the [n_mean()] round trip read the
#'     design back from here.}
#' }
#'
#' Nothing here is rounded: `n` is taken as given, so passing a
#' continuous `n` back from [n_mean()] reproduces its `se`, `moe` and
#' `cv` exactly.
#'
#' @details
#' Computes the standard error for the given sample size and design
#' parameters, then derives the margin of error and coefficient of
#' variation. The effective sample size is `n * resp_rate / deff`, with
#' optional finite population correction.
#' The quantile `q` is `qnorm(1 - alpha / 2)` by default and
#' `qt(1 - alpha / 2, df)` when `df` is supplied.
#'
#' Supplying `cv` in place of `mu` solves the same equation in the remaining
#' direction, returning the smallest mean the design measures that precisely.
#' The standard error of a mean does not involve the mean, so this is a
#' division rather than a search, and only the magnitude is recoverable: `cv`
#' is defined against `abs(mu)`, and the positive root is returned. See
#' [prec_prop()] for the proportion case, where the same reading answers which
#' estimates a fielded design can carry.
#'
#' ## Round-trip with n_mean
#'
#' `prec_mean()` is the inverse of [n_mean()]: if you compute
#' `res <- n_mean(var = 100, moe = 2)` and then call `prec_mean(res)`,
#' you will recover `moe = 2`. You can also pass an `svyplan_n` object
#' directly: `prec_mean(res)`.
#'
#' ## Precision of a total
#'
#' A separate `prec_total()` is not needed. The precision of the
#' estimated total \eqn{\hat{Y} = N \bar{y}}{Yhat = N ybar} is a rescaling of the
#' mean precision:
#'
#' - \eqn{SE(\hat{Y}) = N \times SE(\bar{y})}{SE(Yhat) = N * SE(ybar)}
#' - \eqn{MOE(\hat{Y}) = N \times MOE(\bar{y})}{MOE(Yhat) = N * MOE(ybar)}
#' - \eqn{CV(\hat{Y}) = CV(\bar{y})}{CV(Yhat) = CV(ybar)}
#'
#' Call `prec_mean()` and multiply `$se` and `$moe` by \eqn{N} for
#' the total. The `$cv` component is identical.
#' See Examples below.
#'
#' @family sample size and precision functions
#' @seealso [n_mean()] for the inverse (compute n from a precision target),
#'   [prec_prop()] for proportions.
#'
#' @examples
#' # Precision with n = 400
#' prec_mean(var = 100, n = 400, mu = 50)
#'
#' # Without mu (CV will be NA)
#' prec_mean(var = 100, n = 400)
#'
#' # Smallest mean 400 units can report at a 5 percent CV
#' prec_mean(var = 100, n = 400, cv = 0.05)$params$mu
#'
#' # Round-trip from n_mean
#' res <- n_mean(var = 100, moe = 2)
#' prec_mean(res)
#'
#' ## Precision of a total
#' # Precision of a mean, then scale to total
#' N <- 10000
#' p <- prec_mean(var = 2500, n = 400, mu = 300, N = N)
#' p$moe * N   # MOE for the estimated total
#' p$se * N    # SE for the estimated total
#' p$cv        # CV is the same for mean and total
#'
#' @export
prec_mean <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", prec_mean.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_mean")
}

#' @rdname prec_mean
#' @export
prec_mean.default <- function(
  var = NULL,
  n,
  ...,
  sd = NULL,
  mu = NULL,
  cv = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, prec_mean.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_mean.default, c(.plan, list(...))))
  .check_unused_dots(...)
  var <- .resolve_var(var, sd)
  check_scalar(var, "var")
  check_scalar(n, "n")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  .check_gross_n(n, N)
  if (!is.null(mu) && !is.null(cv)) {
    stop("specify at most one of 'mu' or 'cv'", call. = FALSE)
  }
  if (!is.null(mu)) {
    check_mu(mu)
  }

  if (!is.null(df)) check_df(df)

  solved <- if (!is.null(cv)) "mu" else NULL
  if (!is.null(solved)) {
    check_scalar(cv, "cv")
    mu <- .prec_solve_mean(cv, var, n, alpha, N, deff, resp_rate, df)
  }

  prec <- .prec_engine_mean(var, mu, n, alpha, N, deff, resp_rate, df)
  se <- prec$se
  moe <- prec$moe
  cv_val <- prec$cv

  params <- list(
    var = var,
    n = n,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  if (!is.null(mu)) {
    params$mu <- mu
  }

  .new_svyplan_prec(
    se = se,
    moe = moe,
    cv = cv_val,
    type = "mean",
    params = params,
    solved = solved
  )
}

#' @rdname prec_mean
#' @export
prec_mean.svyplan_n <- function(var, ...) {
  x <- var
  if (x$type != "mean") {
    stop("prec_mean requires a svyplan_n of type 'mean'", call. = FALSE)
  }
  par <- x$params
  args <- list(
    var = par$var,
    n = x$n,
    mu = par$mu,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate %||% 1,
    df = par$df
  )
  do.call(prec_mean.default, .roundtrip_args(args, list(...), prec_mean.default))
}
