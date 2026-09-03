#' Precision of a ratio of two totals at a given sample size
#'
#' Evaluate the standard error, margin of error, and coefficient of variation
#' a planned design achieves for a ratio of two population totals.
#'
#' @section Why there is no solve-for-level mode:
#'
#' [prec_mean()] accepts a target `cv` and solves for the smallest detectable
#' `mu`, because the standard error of a mean does not involve the mean. That
#' inverse does not exist here. The relative standard error of a ratio,
#'
#' \deqn{CV(\hat{R}) = \sqrt{\mathrm{deff} \cdot L_R (1/n_{net} - 1/N)}}{CV(Rhat) = sqrt(deff * L_R * (1/n_net - 1/N))}
#'
#' does not involve the magnitude of `R` at all, so either every non-zero
#' ratio with those component moments meets a `cv` target or none does.
#'
#' A `moe` target is algebraically invertible, giving
#' \eqn{|R| = \mathrm{moe} / (q \sqrt{L_R \cdot fpc / n_{eff}})}{abs(R) = moe / (q * sqrt(L_R * fpc / n_eff))}, and is still
#' not offered. It recovers only the magnitude, so it cannot return the
#' estimand that was asked for, and the question it answers, the largest
#' ratio whose absolute margin of error stays inside a bound at a fixed
#' sample size, is not one survey planning asks.
#'
#' @param r For the default method: the anticipated ratio, `mean(y) / mean(x)`.
#'   May be negative, but not zero.
#'   For `svyplan_n` objects: a sample size result from [n_ratio()].
#' @param n Sample size to evaluate, gross of nonresponse. The responding
#'   sample is `n * resp_rate`.
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param cv_num Coefficient of variation of the numerator variable,
#'   strictly positive.
#' @param cv_den Coefficient of variation of the denominator variable,
#'   strictly positive.
#' @param component_cor Correlation between numerator and denominator across
#'   units, in \[-1, 1\].
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect of the ratio estimator (> 0). See [n_ratio()]
#'   for what this is the design effect of.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1.
#' @param df Degrees of freedom of the variance estimator, switching the
#'   interval quantile from normal to t. `NULL` (default) applies none.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "ratio"` and
#'   `method = "linearization"`:
#' \describe{
#'   \item{`se`}{Standard error of the estimated ratio, on the ratio scale.}
#'   \item{`moe`}{Half-width of the confidence interval, `q * se`.}
#'   \item{`cv`}{Relative standard error, `se / abs(r)`.}
#'   \item{`rmoe`}{Margin of error relative to `r`, `moe / abs(r)`.}
#'   \item{`params`}{The inputs, plus `unit_relvar` and `n`.}
#' }
#'
#' @details
#'
#' At a gross sample `n`, with `n_net = n * resp_rate` and
#' `n_eff = n_net / deff`,
#'
#' \deqn{SE(\hat{R}) = |R| \sqrt{L_R (1 - n_{net}/N) / n_{eff}}}{SE(Rhat) = abs(R) * sqrt(L_R * (1 - n_net/N) / n_eff)}
#'
#' where \eqn{L_R} is the unit relative variance defined in [n_ratio()]. The
#' interval is the symmetric first-order one and is not a Fieller interval.
#' The method omits the ratio estimator's bias and assumes a denominator
#' safely away from zero.
#'
#' @family proportion, mean and ratio functions
#' @seealso [n_ratio()] for the inverse, [prec_mean()] for a mean,
#'   [prec_cluster()] for a multistage design.
#'
#' @examples
#' # Precision of an issued sample of 1200
#' prec_ratio(r = 420, n = 1200, cv_num = 1.20, cv_den = 0.45,
#'            component_cor = 0.65)
#'
#' # With a design effect and nonresponse
#' prec_ratio(r = 420, n = 1200, cv_num = 1.20, cv_den = 0.45,
#'            component_cor = 0.65, deff = 1.3, resp_rate = 0.85)
#'
#' # Round trip: a size and the precision it buys agree exactly
#' size <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
#'                 cv = 0.05)
#' prec_ratio(size)$cv
#'
#' @export
prec_ratio <- function(r = NULL, ...) {
  if (!is.null(r)) {
    .res <- .dispatch_plan(r, "r", prec_ratio.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_ratio")
}

#' @rdname prec_ratio
#' @export
prec_ratio.default <- function(
  r = NULL,
  n = NULL,
  ...,
  cv_num = NULL,
  cv_den = NULL,
  component_cor = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, prec_ratio.default, match.call(),
                            environment())
  if (!is.null(.plan)) {
    return(do.call(prec_ratio.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  .check_ratio_moments(r, cv_num, cv_den, component_cor)
  check_scalar(n, "n")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  .check_gross_n(n, N)
  if (!is.null(df)) {
    check_df(df)
  }

  unit_relvar <- .ratio_unit_relvar(r, cv_num, cv_den, component_cor)
  .check_ratio_relvar(unit_relvar, cv_num, cv_den, component_cor)

  prec <- .prec_engine_ratio(r, unit_relvar, n, alpha, N, deff, resp_rate, df)

  params <- list(
    r = r,
    cv_num = cv_num,
    cv_den = cv_den,
    component_cor = component_cor,
    unit_relvar = unit_relvar,
    n = n,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  .new_svyplan_prec(
    se = prec$se,
    moe = prec$moe,
    cv = prec$cv,
    type = "ratio",
    method = "linearization",
    params = params
  )
}

#' @rdname prec_ratio
#' @export
prec_ratio.svyplan_n <- function(r, ...) {
  x <- r
  if (x$type != "ratio") {
    stop("prec_ratio requires a svyplan_n of type 'ratio'", call. = FALSE)
  }
  par <- x$params
  args <- list(
    r = par$r,
    n = x$n,
    cv_num = par$cv_num,
    cv_den = par$cv_den,
    component_cor = par$component_cor,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate %||% 1,
    df = par$df
  )
  do.call(prec_ratio.default,
          .roundtrip_args(args, list(...), prec_ratio.default))
}
