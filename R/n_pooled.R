#' Sample size for the average of several occasions of a repeated survey
#'
#' Compute the sample size per occasion required to estimate the average of
#' several occasions of a repeated survey, the equal-weight mean of the
#' occasion estimates, with a specified margin of error or coefficient of
#' variation, given how far the occasions overlap. An annual average built
#' from quarterly rounds is the ordinary case, and that average is what
#' "pooled" names here. This is the inverse of [prec_pooled()].
#'
#' @param var For the default method: the population variance \eqn{S^2} on
#'   one occasion, taken to be the same on each. Estimate from a pilot
#'   study, a previous round, or published data. Supply this, `p` or `r`.
#'   For `svyplan_prec` objects: a precision result from [prec_pooled()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param occasions Number of occasions entering the average, at least 2.
#' @param sd Population standard deviation on one occasion, an alternative
#'   spelling of `var`. Supply exactly one of `var` or `sd`.
#' @param p The proportion being averaged, as an alternative to `var`. The
#'   occasion variance is then \eqn{Np(1-p)/(N-1)}, the same finite
#'   population variance [n_prop()] uses, and `mu` is determined and must not
#'   be supplied.
#' @param r The ratio of totals being averaged, `mean(y) / mean(x)` as in
#'   [n_ratio()], taken to be the same on each occasion, as an alternative
#'   to `var` and `p`. `mu` is then `r` and must not be supplied. See the
#'   section on a ratio of totals.
#' @param cv_num Coefficient of variation of the ratio's numerator
#'   variable. Used with `r` only.
#' @param cv_den Coefficient of variation of the ratio's denominator
#'   variable. Used with `r` only.
#' @param component_cor Correlation between the ratio's numerator and
#'   denominator on one unit, in \eqn{[-1, 1]}. Used with `r` only. It is
#'   not `overlap_cor`, which correlates the occasions.
#' @param mu The level being averaged, on the `var` scale. Required when `cv`
#'   or `rmoe` is specified, both being defined against it. Determined by `p`
#'   on the proportion scale.
#' @param moe Desired margin of error on the pooled estimate, the half-width
#'   of its confidence interval. Specify exactly one of `moe`, `cv`, or
#'   `rmoe`.
#' @param cv Target standard error relative to the level. Requires `mu`.
#'   Specify exactly one of `moe`, `cv`, or `rmoe`.
#' @param rmoe Target margin of error relative to the level. It is
#'   `moe / abs(mu)` and requires `mu`. Specify exactly one of `moe`, `cv`,
#'   or `rmoe`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0), applied to the variance of the
#'   pooled estimate rather than to any one occasion.
#' @param resp_rate Expected response rate among eligible units, in (0, 1\].
#'   Default 1 (no adjustment). Each occasion's size is inflated by
#'   `1 / resp_rate`. It is a single round's response, not attrition across a
#'   panel. See the nonresponse section of [svyplan-package] for what the rate
#'   counts and what this adjustment does and does not claim.
#' @param eligibility_rate Share of issued units in scope, such as dwellings
#'   found occupied, in (0, 1\]. Default 1. It applies to every occasion, and
#'   with a finite `N` the population in scope is `N * eligibility_rate`. See
#'   the nonresponse section of [svyplan-package].
#' @param overlap Fraction of one occasion's **responding** sample carried
#'   into a later one, in \eqn{[0, 1]}. One number, meaning the same overlap at
#'   every lag, or one per lag. A [design_overlap()] result is accepted
#'   directly at `resp_rate = 1`. See [prec_pooled()] on why not below it.
#'   `0` (default) makes the occasions independent.
#' @param overlap_cor Correlation between two occasions among the units they
#'   share, in \eqn{[0, 1]}. One number for every lag, or one per lag.
#'   Default 0, which makes overlap irrelevant. Supply this or `cor_decay`,
#'   not both. For a ratio it is the correlation of the ratio's linearized
#'   variable.
#' @param cor_decay Correlation at lag \eqn{m} taken as \eqn{\rho^m}{rho^m},
#'   an alternative to stating `overlap_cor` lag by lag.
#' @param overlap_draw How the shared units are drawn, `"fixed"` (default)
#'   for a fixed count drawn with the rest of each occasion's sample, or
#'   `"independent"` for a retained panel topped up on each occasion by a
#'   fresh sample drawn independently of it from the whole frame. The
#'   second is defined for two occasions, or for more when no lag overlaps.
#'   See [prec_pooled()].
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, available from [design_df()]. It switches the interval
#'   quantile from normal to t. `NULL` (default) applies no adjustment.
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param .overlap_basis Internal. Records whether a stored `overlap`
#'   profile was resolved from a [design_overlap()] object (`"issued"`) or
#'   supplied as respondent overlap (`"respondent"`), so that a round trip
#'   or a grid meets the same refusal the first call would have. Set from
#'   the result being re-read, and there is no reason to pass it by hand.
#'
#' @return A `svyplan_n` object with `type = "pooled"`:
#' \describe{
#'   \item{`n`}{Required size per occasion, continuous and gross. It already
#'     carries `deff` and the `1 / resp_rate` inflation, so it counts the
#'     units to release on each occasion, not the completed interviews. `$n`
#'     and `as.double()` keep the unrounded value, which is what makes the
#'     round trip through [prec_pooled()] exact, whereas `print()` rounds
#'     up. This
#'     is not a count of distinct frame units: `overlap` is measured among
#'     respondents, so it says how often a respondent is measured again and
#'     leaves the issued sample's own overlap unstated. At full response the
#'     two coincide and the series consumes fewer than `occasions * n`
#'     distinct units. Below it, how many is a question this function has
#'     not been told enough to answer.}
#'   \item{`se`, `moe`, `cv`, `rmoe`}{Precision the design achieves at that
#'     size, the same values [prec_pooled()] reports for the same inputs.
#'     `cv` and `rmoe` are `NA` unless the level is known.}
#'   \item{`params`}{The validated inputs, including whichever of `moe`,
#'     `cv`, or `rmoe` was the target, `overlap` and `overlap_cor` as
#'     resolved vectors of one entry per lag, and the `overlap_basis` they
#'     were resolved under. A ratio keeps `r`, its three moments and its
#'     `unit_relvar`, `var` then being the variance of its linearized
#'     variable.}
#' }
#'
#' @details
#' A required size must fit the finite frame at each positive-overlap lag:
#' `(2 - overlap[lag]) * n * resp_rate <= N`. A target violating a pairwise
#' bound is rejected, as is a size producing an invalid covariance matrix.
#' See [prec_pooled()] for these necessary feasibility conditions.
#'
#' The pooled variance is \eqn{deff\{A/(n r) + B\}}{deff (A/(n r) + B)} in the gross size per
#' occasion, with \eqn{r} the response rate times the eligibility rate,
#' \eqn{N} the population in scope, `N * eligibility_rate`, \eqn{T} the
#' number of occasions and
#'
#' \deqn{A = \frac{S^2}{T^2}\Big[T + 2\sum_{m=1}^{T-1}(T-m)\rho_m o_m\Big],
#'       \qquad
#'       B = -\frac{S^2}{NT^2}\Big[T + 2\sum_{m=1}^{T-1}(T-m)\rho_m
#'       \mathbb{1}\{o_m > 0\}\Big],}{A = (S^2/T^2) [T + 2 sum_m (T - m) rho_m o_m],   B = -(S^2/(N T^2)) [T + 2 sum_m (T - m) rho_m 1\{o_m > 0\}],}
#'
#' so the target inverts in closed form,
#'
#' \deqn{n = \frac{deff\,A}{r\,(se^2 - deff\,B)}.}{n = deff A / (r (se^2 - deff B)).}
#'
#' \eqn{B \le 0} always under the nonnegative correlation contract, since its
#' bracket is at least \eqn{T}, so the divisor cannot vanish and no target is
#' out of reach for want of precision however much the occasions overlap.
#' Overlap raises the size a target needs, but it does not put a floor under
#' the precision. The size can still exceed \eqn{N}, which is the ordinary
#' boundary and is reported as unattainable there as everywhere else.
#'
#' ## A ratio of totals
#'
#' A ratio \eqn{R = \bar Y / \bar X}{R = Ybar / Xbar} enters through its
#' linearized variable \eqn{z = (y - Rx)/\bar X}{z = (y - R x) / Xbar},
#' whose variance is \eqn{R^2 L}, with \eqn{L} the unit relative variance
#' [n_ratio()] reports. That variance takes the place of `var` here and is
#' what `params$var` records, and `overlap_cor` and `cor_decay` describe the
#' correlation of \eqn{z} between occasions. A ratio whose level moves
#' between occasions is outside the model, as a moving mean is.
#'
#' What is sized is the equal-weight average of the occasion ratios. A ratio
#' of the pooled totals, \eqn{\sum_t \hat Y_t / \sum_t \hat X_t}{sum_t Yhat_t / sum_t Xhat_t},
#' is a different estimator even at one level. To first order its error is
#' \eqn{\sum_t w_t(\hat R_t - R)}{sum_t w_t (Rhat_t - R)}, with
#' \eqn{w_t} the occasion's share of the denominator total, so it has the
#' variance computed here only when the denominator totals are equal across
#' occasions. Otherwise its variance differs, and over independent
#' occasions it is larger, equal weights being the ones that minimize it
#' there. Over two independent occasions whose denominator totals stand 1
#' to 3, it is 25 percent larger than the equal-weight average's. Sizing an annual rate defined on pooled
#' totals with this function treats the denominator totals as equal, which
#' is close for a stable population and is an approximation to state where
#' it is not.
#'
#' ## Sizing for a pooled estimate or for a change
#'
#' The mean number of visits per distinct unit, which [design_overlap()]
#' reports for a window, measures burden and is not the variance multiplier
#' of a pooled average. That multiplier follows from the overlap at every
#' lag. Under `"2-2-2"` over four quarters with `overlap_cor = 1` and no
#' finite population correction it is 15/8, against 16/9 visits per unit,
#' since a unit interviewed twice carries twice the weight.
#'
#' The two arms pull in opposite directions against one design lever, and
#' which arm gains depends on the sign of the covariance the overlap
#' induces. Above a sampling fraction of \eqn{o_m > n/N}{o_m > n/N}, the
#' ordinary case and the only one without a finite population correction,
#' overlap improves a change and raises the size a pooled target needs.
#' Below it the occasions share fewer units than chance would give them and
#' the directions reverse. [prec_pooled()] works the boundary. The level at
#' a single occasion is unaffected either way. A design serving a pooled
#' target and a change target is sized by taking the larger of `n_pooled()`
#' and [n_change()], since neither dominates.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley. Chapter 12.
#'
#' @family repeated survey functions
#' @seealso [prec_pooled()] for the inverse (compute precision from a size)
#'   and for what the covariance assumes, [n_change()] for the other arm of
#'   the trade-off, [design_overlap()] for the overlap a rotation schedule
#'   gives, [n_mean()] and [n_prop()] for a single occasion.
#'
#' @examples
#' # Four independent quarterly rounds averaged into an annual figure
#' n_pooled(var = 100, moe = 1, occasions = 4)
#'
#' # The same target from a rotating panel needs more per occasion
#' n_pooled(var = 100, moe = 1, occasions = 4, overlap = 0.75,
#'          cor_decay = 0.8)
#'
#' # An annual average of a proportion, to within one point
#' n_pooled(p = 0.3, moe = 0.01, occasions = 12, overlap = 0.75,
#'          cor_decay = 0.9)
#'
#' # An annual unemployment rate from four quarters, to a cv of 3 percent
#' n_pooled(r = 0.0615, cv_num = 4.9, cv_den = 0.73, component_cor = 0.15,
#'          cv = 0.03, occasions = 4, overlap = 0.75, cor_decay = 0.8)
#'
#' # Relative target: the average to within 5 percent of itself
#' n_pooled(var = 100, mu = 20, rmoe = 0.05, occasions = 4)
#'
#' # With FPC, design effect, and response
#' n_pooled(var = 100, moe = 1, occasions = 4, N = 50000, deff = 1.5,
#'          resp_rate = 0.8)
#'
#' # Two occasions on a small frame, a retained panel plus independent draws
#' n_pooled(var = 100, moe = 2.5, occasions = 2, N = 400, overlap = 0.5,
#'          overlap_cor = 0.8, overlap_draw = "independent")
#'
#' @export
n_pooled <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", n_pooled.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_pooled")
}

#' @rdname n_pooled
#' @export
n_pooled.default <- function(
  var = NULL,
  occasions,
  ...,
  sd = NULL,
  p = NULL,
  r = NULL,
  cv_num = NULL,
  cv_den = NULL,
  component_cor = NULL,
  mu = NULL,
  moe = NULL,
  cv = NULL,
  rmoe = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  eligibility_rate = 1,
  overlap = 0,
  overlap_cor = NULL,
  cor_decay = NULL,
  overlap_draw = c("fixed", "independent"),
  df = NULL,
  plan = NULL,
  .overlap_basis = NULL
) {
  .plan <- .merge_plan_args(plan, n_pooled.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_pooled.default, c(.plan, list(...))))
  .check_unused_dots(...)
  check_population_size(N)
  check_resp_rate(eligibility_rate, "eligibility_rate")
  N_elig <- N * eligibility_rate
  .check_occasions(occasions)
  inputs <- .pooled_inputs(var, sd, p, mu, N_elig, r, cv_num, cv_den,
                           component_cor)
  s <- inputs$scale
  check_precision(moe, cv, rmoe)
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)
  rate <- resp_rate * eligibility_rate
  if (!is.null(df)) check_df(df)
  ov <- .resolve_lag_overlap(overlap, occasions, resp_rate, .overlap_basis)
  rho <- .resolve_lag_cor(overlap_cor, cor_decay, occasions - 1L)
  overlap_draw <- match.arg(overlap_draw)
  .check_pooled_draw(overlap_draw, ov$overlap, occasions)

  q <- .q_alpha(alpha, df)
  if (!is.null(cv)) {
    if (is.null(inputs$mu)) {
      stop("'mu' is required when 'cv' is specified", call. = FALSE)
    }
    if (inputs$mu == 0) {
      stop("'cv' is undefined at 'mu' = 0", call. = FALSE)
    }
    target_se <- cv * abs(inputs$mu)
  } else {
    # On the ratio scale the engine works in units of `s`, so an absolute
    # margin is divided by it. A relative one is already scale-free.
    moe_used <- if (is.null(rmoe)) {
      moe / s
    } else {
      .moe_from_rmoe(rmoe, inputs$mu, "mu")
    }
    target_se <- moe_used / q
  }

  n <- .n_pooled_from_se(
    target_se, inputs$var, N_elig, deff, rate, occasions, ov$overlap, rho,
    overlap_draw
  )
  n <- .census_size(n, N)
  .check_attainable(n, N, resp_rate, eligibility_rate)

  prec <- .prec_engine_pooled(
    n, inputs$var, inputs$mu, N_elig, deff, rate, occasions, ov$overlap, rho,
    overlap_draw, alpha, df
  )
  # The inversion is closed form on n <= N, so a mismatch here means the
  # linear form and the variance it inverts have drifted apart, not that
  # the inputs were unusual.
  if (abs(prec$se - target_se) > 1e-6 * max(1, target_se)) {
    stop(
      "internal error: solved size does not reproduce the target standard error",
      call. = FALSE
    )
  }

  params <- list(
    var = inputs$var_stored,
    occasions = occasions,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    overlap = ov$overlap,
    overlap_cor = rho,
    overlap_basis = ov$basis,
    overlap_draw = overlap_draw,
    df = df
  )
  if (eligibility_rate != 1) params$eligibility_rate <- eligibility_rate
  if (!is.null(inputs$p)) params$p <- inputs$p
  params <- .add_ratio_params(params, inputs$ratio)
  if (!is.null(inputs$mu)) params$mu <- inputs$mu_stored
  if (!is.null(moe)) {
    params$moe <- moe
  } else if (!is.null(rmoe)) {
    params$rmoe <- rmoe
  } else {
    params$cv <- cv
  }

  .new_svyplan_n(
    n = n,
    type = "pooled",
    params = params,
    se = prec$se * s,
    moe = prec$moe * s,
    cv = prec$cv
  )
}

#' @rdname n_pooled
#' @export
n_pooled.svyplan_prec <- function(var, ..., moe = NULL, cv = NULL,
                                  rmoe = NULL) {
  x <- var
  if (x$type != "pooled") {
    stop("n_pooled requires a svyplan_prec of type 'pooled'", call. = FALSE)
  }
  par <- x$params
  # The achieved margin of error is the implied target, but only when the
  # caller named none of the three: restoring it alongside an override
  # would send two targets into a function that takes one.
  if (is.null(moe) && is.null(cv) && is.null(rmoe)) {
    moe <- x$moe
  }
  args <- c(.repeated_scale_args(par, "mu"), list(
    occasions = par$occasions,
    moe = moe,
    cv = cv,
    rmoe = rmoe,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    eligibility_rate = par[["eligibility_rate"]] %||% 1,
    overlap = par$overlap,
    overlap_cor = par$overlap_cor,
    overlap_draw = par$overlap_draw,
    df = par$df,
    .overlap_basis = par$overlap_basis
  ))
  do.call(n_pooled.default, .roundtrip_args(args, list(...), n_pooled.default))
}
