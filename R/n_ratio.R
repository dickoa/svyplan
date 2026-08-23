#' Sample size for a ratio of two totals
#'
#' Compute the required sample size for estimating a ratio of two population
#' totals, such as consumption per person or yield per hectare, with a
#' specified margin of error or coefficient of variation.
#'
#' @section What this ratio is:
#'
#' The estimand is \eqn{R = T_y / T_x}, the ratio of two totals from one
#' finite population, equivalently the ratio of two means observed on the
#' same sampled units. The numerator and denominator must be measured on the
#' same units, and `resp_rate` is the expected rate of records usable for
#' both.
#'
#' It is not the population mean of the unit ratios \eqn{y_i / x_i}, which is
#' an ordinary mean and belongs in [n_mean()]. It is not a ratio whose parts
#' come from different samples, nor a comparison of two ratios, nor a
#' regression slope.
#'
#' A rate whose denominator is a subpopulation count, such as a coverage rate
#' among an eligible group, is usually easier to plan as a domain proportion
#' with [n_prop()]. This function earns its place when the denominator is a
#' continuous per-unit quantity, where there is no domain share to inflate by.
#'
#' @section Design effect and clustering:
#'
#' `deff` is the design effect of the ratio estimator, or first-order
#' equivalently of the linearized variable \eqn{e = y - Rx}. It is not the
#' design effect of the numerator, nor of the denominator, and it cannot be
#' derived from their design effects. The clustering loss on a ratio is
#' routinely an order of magnitude away from the loss on either component,
#' because the correlation that makes a ratio precise also removes much of
#' the between-cluster variation.
#'
#' For a multistage plan, compute the linearized variable from pilot data and
#' pass its variance components to [n_cluster()] rather than passing a `deff`
#' here. See [varcomp()] and the example below.
#'
#' @param r For the default method: the anticipated ratio, `mean(y) / mean(x)`.
#'   It may be negative but not zero, since a relative precision has no scale
#'   at zero. Named `r` rather than `ratio`, which already means an allocation
#'   ratio in [svyplan()] profiles.
#'   For `svyplan_prec` objects: a precision result from [prec_ratio()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param cv_num Coefficient of variation of the numerator variable,
#'   `sd(y) / abs(mean(y))`, taken on the magnitude so it carries no sign.
#'   It must be strictly positive: a component with no variation at all is
#'   not a sampling problem, and at `cv_num = 0` the ratio's own variance
#'   comes entirely from the denominator.
#' @param cv_den Coefficient of variation of the denominator variable,
#'   `sd(x) / abs(mean(x))`.
#' @param component_cor Correlation between the numerator and denominator
#'   variables across units, in \[-1, 1\]. Distinct from `overlap_cor`, which
#'   is a correlation between occasions. What makes a ratio precise is
#'   `component_cor * sign(r)` being close to 1 while the two coefficients of
#'   variation are similar, so for a positive ratio that means a strongly
#'   positive correlation and for a negative ratio a strongly negative one.
#' @param moe Desired margin of error on the ratio scale, the half-width of
#'   the confidence interval in the units of `r`. Specify exactly one of
#'   `moe`, `cv`, or `rmoe`.
#' @param cv Target coefficient of variation of the estimated ratio. Specify
#'   exactly one of `moe`, `cv`, or `rmoe`.
#' @param rmoe Target margin of error relative to `r`, so `rmoe = 0.05` asks
#'   for an interval whose half-width is 5 percent of the ratio. Specify
#'   exactly one of `moe`, `cv`, or `rmoe`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size, in the units the ratio is observed on. `Inf`
#'   (default) means no finite population correction.
#' @param deff Design effect of the ratio estimator (> 0). See the section
#'   above, which is the one place this differs from [n_mean()].
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1. The rate
#'   is for records usable on both components, not on either alone.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, available from [design_df()]. It switches the interval
#'   quantile from normal to t. `NULL` (default) applies no adjustment.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "ratio"` and
#'   `method = "linearization"`:
#' \describe{
#'   \item{`n`}{Required sample size, continuous and gross. It carries `deff`
#'     and the `1 / resp_rate` inflation, so it counts units to release
#'     rather than completed interviews. `print()` and `as.integer()` round
#'     it up.}
#'   \item{`se`, `moe`, `cv`, `rmoe`}{Precision achieved at that `n`, on the
#'     ratio scale.}
#'   \item{`params`}{The inputs, plus `unit_relvar`, the coefficient
#'     \eqn{L_R} defined below.}
#' }
#'
#' @details
#'
#' Write \eqn{CV_y} and \eqn{CV_x} for the component coefficients of
#' variation and \eqn{\rho} for their correlation. The unit relative variance
#' of the ratio is
#'
#' \deqn{L_R = CV_y^2 + CV_x^2 - 2 \rho \, CV_y CV_x \, \mathrm{sign}(R)}{L_R = CV_y^2 + CV_x^2 - 2 rho CV_y CV_x sign(R)}
#'
#' which equals \eqn{S_e^2 / \mu_y^2} for \eqn{e = y - Rx}. The `sign(R)`
#' factor matters only when the two means differ in sign, and it keeps the
#' result unchanged when either component is negated.
#'
#' The required responding sample for a target relative standard error
#' \eqn{c} is
#'
#' \deqn{n = \frac{\mathrm{deff} \cdot L_R}{c^2 + \mathrm{deff} \cdot L_R / N}}{n = deff * L_R / (c^2 + deff * L_R / N)}
#'
#' and the gross sample is that divided by `resp_rate`. A margin of error
#' converts to \eqn{c} exactly, as
#' \eqn{c = \mathrm{moe} / (q |R|)}{c = moe / (q * abs(R))}.
#'
#' The method is first order. It omits the ratio estimator's
#' \eqn{O(1/n)} bias and assumes the denominator is far enough from zero
#' that its sampling distribution does not approach it. A large `cv_den`, a
#' skewed denominator, or a small effective sample all call for simulation
#' rather than this approximation. The interval is the symmetric
#' \eqn{\hat{R} \pm q \cdot SE}{Rhat +/- q * SE}, not a Fieller interval.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @family sample size functions
#' @seealso [prec_ratio()] for the inverse, [n_mean()] for a mean,
#'   [n_cluster()] for a multistage design, [varcomp()] for the linearized
#'   variable's components.
#'
#' @examples
#' # Per-capita household consumption: total consumption over total members
#' n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45, component_cor = 0.65,
#'         cv = 0.05)
#'
#' # Margin of error on the ratio scale
#' n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45, component_cor = 0.65,
#'         moe = 20)
#'
#' # Relative margin of error, with a design effect and nonresponse
#' n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45, component_cor = 0.65,
#'         rmoe = 0.05, deff = 1.3, resp_rate = 0.85)
#'
#' # A high correlation between the parts makes the ratio cheap to estimate
#' n_ratio(r = 420, cv_num = 1.20, cv_den = 1.15, component_cor = 0.95,
#'         cv = 0.05)
#'
#' ## Planning a clustered ratio from pilot data
#' # The homogeneity that matters belongs to e = y - R * x, not to y or x.
#' set.seed(1)
#' psu <- rep(1:30, each = 20)
#' x <- stats::rlnorm(600, 1, 0.4)
#' y <- 3 * x + stats::rnorm(600, 0, 0.8)
#' R <- mean(y) / mean(x)
#' e <- y - R * x
#' vc <- suppressWarnings(varcomp(e, stage_id = list(psu)))
#' # varcomp reports an infinite unit_relvar because e has mean zero, so take
#' # the ratios from it and supply the ratio's own coefficient separately.
#' n_cluster(stage_cost = c(500, 50), icc = vc$icc,
#'           unit_relvar = stats::var(e) / mean(y)^2,
#'           var_ratio = vc$var_ratio, cv = 0.05)
#'
#' @export
n_ratio <- function(r = NULL, ...) {
  if (!is.null(r)) {
    .res <- .dispatch_plan(r, "r", n_ratio.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_ratio")
}

#' @rdname n_ratio
#' @export
n_ratio.default <- function(
  r = NULL,
  ...,
  cv_num = NULL,
  cv_den = NULL,
  component_cor = NULL,
  moe = NULL,
  cv = NULL,
  rmoe = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_ratio.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_ratio.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  .check_ratio_moments(r, cv_num, cv_den, component_cor)
  check_precision(moe, cv, rmoe)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  if (!is.null(df)) {
    check_df(df)
  }

  unit_relvar <- .ratio_unit_relvar(r, cv_num, cv_den, component_cor)
  .check_ratio_relvar(unit_relvar, cv_num, cv_den, component_cor)

  n <- .n_ratio_from_target(r, unit_relvar, moe, cv, rmoe, alpha, N, deff,
                            resp_rate, df)

  params <- list(
    r = r,
    cv_num = cv_num,
    cv_den = cv_den,
    component_cor = component_cor,
    unit_relvar = unit_relvar,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  if (!is.null(moe)) {
    params$moe <- moe
  } else if (!is.null(rmoe)) {
    params$rmoe <- rmoe
  } else {
    params$cv <- cv
  }

  .new_svyplan_n(
    n = n,
    type = "ratio",
    method = "linearization",
    params = params
  )
}

#' @rdname n_ratio
#' @export
n_ratio.svyplan_prec <- function(r, ..., moe = NULL, cv = NULL, rmoe = NULL) {
  x <- r
  if (x$type != "ratio") {
    stop("n_ratio requires a svyplan_prec of type 'ratio'", call. = FALSE)
  }
  par <- x$params
  if (is.null(moe) && is.null(cv) && is.null(rmoe)) {
    moe <- x$moe
  }
  args <- list(
    r = par$r,
    cv_num = par$cv_num,
    cv_den = par$cv_den,
    component_cor = par$component_cor,
    rmoe = rmoe,
    moe = moe,
    cv = cv,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    df = par$df
  )
  do.call(n_ratio.default, .roundtrip_args(args, list(...), n_ratio.default))
}
