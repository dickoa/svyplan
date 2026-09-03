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
#'   study, a previous round, or published data. Supply this or `p`, not
#'   both. For `svyplan_prec` objects: a precision result from
#'   [prec_pooled()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param occasions Number of occasions entering the average, at least 2.
#' @param sd Population standard deviation on one occasion, an alternative
#'   spelling of `var`. Supply exactly one of `var` or `sd`.
#' @param p The proportion being averaged, as an alternative to `var`. The
#'   occasion variance is then \eqn{Np(1-p)/(N-1)}, the same finite
#'   population variance [n_prop()] uses, and `mu` is determined and must not
#'   be supplied.
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
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). Each occasion's size is inflated by `1 / resp_rate`. It is
#'   a single round's response, not attrition across a panel. See the
#'   nonresponse section of [svyplan-package] for what this adjustment does
#'   and does not claim.
#' @param overlap Fraction of one occasion's **responding** sample carried
#'   into a later one, in \[0, 1\]. One number, meaning the same overlap at
#'   every lag, or one per lag. A [design_overlap()] result is accepted
#'   directly at `resp_rate = 1`. See [prec_pooled()] on why not below it.
#'   `0` (default) makes the occasions independent.
#' @param overlap_cor Correlation between two occasions among the units they
#'   share, in \[0, 1\]. One number for every lag, or one per lag. Default 0,
#'   which makes overlap irrelevant. Supply this or `cor_decay`, not both.
#' @param cor_decay Correlation at lag \eqn{m} taken as \eqn{\rho^m}{rho^m},
#'   an alternative to stating `overlap_cor` lag by lag.
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
#'     were resolved under.}
#' }
#'
#' @details
#' The pooled variance is \eqn{deff\{A/(n r) + B\}}{deff (A/(n r) + B)} in the gross size per
#' occasion, with \eqn{r} the response rate, \eqn{T} the number of occasions
#' and
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
#' ## Sizing for a pooled estimate or for a change
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
#' # Relative target: the average to within 5 percent of itself
#' n_pooled(var = 100, mu = 20, rmoe = 0.05, occasions = 4)
#'
#' # With FPC, design effect, and response
#' n_pooled(var = 100, moe = 1, occasions = 4, N = 50000, deff = 1.5,
#'          resp_rate = 0.8)
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
  mu = NULL,
  moe = NULL,
  cv = NULL,
  rmoe = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  overlap = 0,
  overlap_cor = NULL,
  cor_decay = NULL,
  df = NULL,
  plan = NULL,
  .overlap_basis = NULL
) {
  .plan <- .merge_plan_args(plan, n_pooled.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_pooled.default, c(.plan, list(...))))
  .check_unused_dots(...)
  check_population_size(N)
  .check_occasions(occasions)
  inputs <- .pooled_inputs(var, sd, p, mu, N)
  check_precision(moe, cv, rmoe)
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)
  if (!is.null(df)) check_df(df)
  ov <- .resolve_lag_overlap(overlap, occasions, resp_rate, .overlap_basis)
  rho <- .resolve_lag_cor(overlap_cor, cor_decay, occasions - 1L)

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
    moe_used <- if (is.null(rmoe)) {
      moe
    } else {
      .moe_from_rmoe(rmoe, inputs$mu, "mu")
    }
    target_se <- moe_used / q
  }

  n <- .n_pooled_from_se(
    target_se, inputs$var, N, deff, resp_rate, occasions, ov$overlap, rho
  )
  .check_attainable(n, N, resp_rate)

  prec <- .prec_engine_pooled(
    n, inputs$var, inputs$mu, N, deff, resp_rate, occasions, ov$overlap, rho,
    alpha, df
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
    var = inputs$var,
    occasions = occasions,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    overlap = ov$overlap,
    overlap_cor = rho,
    overlap_basis = ov$basis,
    df = df
  )
  if (!is.null(inputs$p)) params$p <- inputs$p
  if (!is.null(inputs$mu)) params$mu <- inputs$mu
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
    se = prec$se,
    moe = prec$moe,
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
  args <- list(
    var = if (is.null(par$p)) par$var else NULL,
    p = par$p,
    occasions = par$occasions,
    mu = if (is.null(par$p)) par$mu else NULL,
    moe = moe,
    cv = cv,
    rmoe = rmoe,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    overlap = par$overlap,
    overlap_cor = par$overlap_cor,
    df = par$df,
    .overlap_basis = par$overlap_basis
  )
  do.call(n_pooled.default, .roundtrip_args(args, list(...), n_pooled.default))
}
