#' Sampling precision for a change between two occasions
#'
#' Compute the sampling error (SE, margin of error, CV) for the change in a
#' mean or a proportion between two occasions of the same population, given
#' the size of each occasion and how far the two samples overlap. This is
#' the inverse of [n_change()].
#'
#' @param var For the default method: the population variance \eqn{S^2} on
#'   each occasion, one value for both or one per occasion. Supply this or
#'   `p`, not both. For `svyplan_n` objects: a sample size result from
#'   [n_change()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param n Sample size, gross units drawn, either one size for both
#'   occasions or one per occasion. Two sizes make `overlap` directional,
#'   since it counts the shared units against the first occasion.
#' @param sd Population standard deviation on each occasion, an alternative
#'   spelling of `var`. Supply exactly one of `var` or `sd`.
#' @param p The two occasion proportions, `c(p1, p2)`, as an alternative to
#'   `var`. The occasion variances are then \eqn{Np(1-p)/(N-1)}, the same
#'   finite population variance [prec_prop()] uses, and the change is
#'   \eqn{p_2 - p_1}, so `change` is determined and must not be supplied.
#' @param change Expected change, \eqn{\mu_2 - \mu_1}, on the `var` scale.
#'   It may be negative. Required only for the relative measures: `cv` and
#'   `rmoe` are defined against it and are `NA` without it. Determined by
#'   `p` on the proportion scale.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction. One value covers both occasions; two are accepted only at
#'   `overlap = 0`, where the occasions are independent and may legitimately
#'   be different populations. A positive `overlap` requires a single `N`,
#'   since units can only be shared by samples drawn from one population.
#' @param deff Design effect multiplier (> 0), applied to the variance of
#'   the change rather than to either occasion separately.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). It applies to both occasions, and nets each size down to
#'   `n * resp_rate` before the variance is formed. It is a single round's
#'   response, not attrition across a panel.
#' @param overlap Fraction of the first occasion's **responding** sample
#'   carried into the second, in \[0, 1\], measured against the first
#'   occasion as `n12 / n1`, so it cannot exceed `n[2] / n[1]`. `0` (default)
#'   drops the covariance entirely and reproduces the two-independent-samples
#'   result. `1` means every responding unit of the first occasion is
#'   measured again, which is a full panel when the occasions are the same
#'   size and a subsample of the second when it is larger. See Details on
#'   what a response rate below 1 does to this reading.
#' @param overlap_cor Correlation between the two occasions among the
#'   overlapping units, in \[0, 1\]. Default 0, which makes overlap
#'   worthless, since it is the product `overlap * overlap_cor` that buys precision
#'   on a change, so a full panel of uncorrelated measurements buys nothing.
#'   On the proportion scale two Bernoulli marginals bound the correlation
#'   they can have, and a value above that bound is rejected.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()]. It switches the interval quantile from normal to t.
#'   `NULL` (default) applies no adjustment.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "change"`:
#' \describe{
#'   \item{`se`}{Standard error of the estimated change, computed on the net
#'     sizes `n * resp_rate`.}
#'   \item{`moe`}{Margin of error, `qnorm(1 - alpha / 2) * se`. The interval
#'     is symmetric about the change, so the limits are `change - moe` and
#'     `change + moe`.}
#'   \item{`cv`}{Standard error relative to the change, `se / abs(change)`.
#'     `NA` when `change` is neither supplied nor determined by `p`. A
#'     change near zero makes it large by construction, which is a statement
#'     about the estimand and not about the design.}
#'   \item{`params`}{The validated inputs. Dispersion is always stored as
#'     `var`, a pair, including when you supplied `sd` or `p`; `p` is kept as
#'     well when the proportion scale was used, and it is what the
#'     [n_change()] round trip reads the scale back from.}
#' }
#'
#' Nothing here is rounded, so passing a continuous `n` back from
#' [n_change()] reproduces its `se`, `moe` and `cv` exactly.
#'
#' @details
#' The estimand is the change in one population measured twice, not a
#' difference between two populations. Writing \eqn{v_1, v_2} for the
#' occasion variances, \eqn{n_1, n_2} for the net sizes and \eqn{k} for the
#' shared units, the expression is piecewise in `overlap`. At
#' `overlap = 0` the two occasions are treated as independent, which is what
#' also lets them come from different populations:
#'
#' \deqn{V = \frac{v_1}{n_1}\Big(1 - \frac{n_1}{N_1}\Big)
#'         + \frac{v_2}{n_2}\Big(1 - \frac{n_2}{N_2}\Big).}{V = v_1/n_1 (1 - n_1/N_1 ) + v_2/n_2 (1 - n_2/N_2 ).}
#'
#' Above zero the occasions share one population of size \eqn{N} and the
#' covariance enters:
#'
#' \deqn{V = \frac{v_1}{n_1} + \frac{v_2}{n_2}
#'       - \frac{2\rho\,\mathrm{overlap}\sqrt{v_1v_2}}{n_2}
#'       - \frac{v_1 + v_2 - 2\rho\sqrt{v_1v_2}}{N},}{V = v_1/n_1 + v_2/n_2 - (2 rho overlap sqrt(v_1v_2))/n_2 - (v_1 + v_2 - 2 rho sqrt(v_1v_2))/N,}
#'
#' a per-unit part less a population part. The marginal terms carry their
#' own finite population correction; the overlap covariance does not, since
#' \eqn{Cov(\bar y_1, \bar y_2) = \rho S_1 S_2 \{k/(n_1n_2) - 1/N\}}{Cov(ybar_1, ybar_2) = rho S_1 S_2 \{k/(n_1n_2) - 1/N\}} enters
#' the population once as \eqn{1/N}. The two expressions agree at
#' \eqn{\rho = 0} and differ otherwise, because the second keeps the
#' population term the first drops along with the rest of the covariance.
#' At `overlap = 1` and `overlap_cor = 1` with equal sizes and variances the
#' population terms cancel exactly and the whole reduces to
#' \eqn{2S^2(1 - \mathrm{overlap})/n}{2S^2(1 - overlap)/n}, which is zero, the same units being
#' measured twice and the change is observed rather than estimated.
#'
#' `deff` multiplies the assembled variance. It is the design effect of the
#' change, which is not in general either occasion's design effect, and a
#' clustered design measuring the same clusters twice will have a smaller
#' one than a design that reclusters.
#'
#' ## What the overlap covariance assumes
#'
#' The covariance is a model form, a correlation `overlap_cor` between
#' occasions among the units the two occasions share, and simple random
#' sampling otherwise. It is not the design-based covariance of two waves of
#' a complex design, which carries the pairwise inclusion probabilities of
#' the master sample the waves were coordinated from. Under a clustered or
#' unequal-probability master the two can differ materially, and this
#' function plans from the schedule rather than measuring a realization.
#' `overlap` and `overlap_cor` are asserted here rather than derived, and nothing checks
#' them against a rotation pattern.
#'
#' The sizes enter the covariance after `resp_rate` has netted them down, so
#' `overlap` is the overlap between the two **responding** samples, not
#' between the issued ones. The two coincide at `resp_rate = 1`. Below it,
#' converting an issued-sample overlap into this one needs an assumption
#' about how response at the second occasion depends on response at the
#' first, which this function does not make. Under independent response the
#' responding overlap is the lower quantity, about `resp_rate` times the
#' issued one. Supply the figure you expect among respondents.
#' [design_overlap()] computes the issued figure a rotation schedule gives
#' and states the conversion.
#'
#' ## Round-trip with n_change
#'
#' `prec_change()` is the inverse of [n_change()]: computing
#' `res <- n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6)`
#' and then `prec_change(res)` recovers `moe = 2`.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley. Chapter 12.
#'
#' @family precision functions
#' @seealso [n_change()] for the inverse (compute n from a precision
#'   target), [prec_mean()] and [prec_prop()] for a single occasion,
#'   [power_mean()] and [power_prop()] to frame the same overlap as a
#'   hypothesis test.
#'
#' @examples
#' # Two independent occasions of 500: the flat two-sample result
#' prec_change(var = 100, n = 500)
#'
#' # A half-overlapping panel, correlated 0.6, buys precision
#' prec_change(var = 100, n = 500, overlap = 0.5, overlap_cor = 0.6)
#'
#' # A change in a proportion: variances and the change come from p
#' prec_change(p = c(0.30, 0.36), n = 1200, overlap = 0.75, overlap_cor = 0.5)
#'
#' # Relative measures need the change on the var scale
#' prec_change(var = 100, n = 500, change = 5)$cv
#'
#' # Unequal occasions: overlap counts against the first
#' prec_change(var = 100, n = c(800, 400), overlap = 0.5, overlap_cor = 0.6)
#'
#' # Round-trip from n_change
#' res <- n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6)
#' prec_change(res)
#'
#' @export
prec_change <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", prec_change.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_change")
}

#' @rdname prec_change
#' @export
prec_change.default <- function(
  var = NULL,
  n,
  ...,
  sd = NULL,
  p = NULL,
  change = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  overlap = 0,
  overlap_cor = 0,
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, prec_change.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_change.default, c(.plan, list(...))))
  .check_unused_dots(...)
  N_pair <- .check_power_N(N)
  inputs <- .change_inputs(var, sd, p, change, N_pair)
  n <- .check_change_n(n)
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)
  check_overlap(overlap)
  check_overlap_cor(overlap_cor)
  .check_overlap_N(overlap, N_pair)
  .check_overlap_n(overlap, n = n)
  .check_bernoulli_cor(inputs$p, overlap_cor, overlap)
  .check_gross_n(
    if (length(n) == 1L) c(n, n) else n,
    N_pair,
    label = c("occasion 1", "occasion 2")
  )
  if (!is.null(df)) check_df(df)

  prec <- .prec_engine_change(
    n, inputs$var_pair, inputs$change, N_pair, deff, resp_rate,
    overlap, overlap_cor, alpha, df
  )

  params <- list(
    var = inputs$var_pair,
    n = n,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    overlap = overlap,
    overlap_cor = overlap_cor,
    df = df
  )
  if (!is.null(inputs$p)) params$p <- inputs$p
  if (!is.null(inputs$change)) params$change <- inputs$change

  .new_svyplan_prec(
    se = prec$se,
    moe = prec$moe,
    cv = prec$cv,
    type = "change",
    params = params
  )
}

#' @rdname prec_change
#' @export
prec_change.svyplan_n <- function(var, ...) {
  x <- var
  if (x$type != "change") {
    stop("prec_change requires a svyplan_n of type 'change'", call. = FALSE)
  }
  par <- x$params
  args <- list(
    var = if (is.null(par$p)) par$var else NULL,
    p = par$p,
    n = x$n,
    change = if (is.null(par$p)) par$change else NULL,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate %||% 1,
    overlap = par$overlap,
    overlap_cor = par$overlap_cor,
    df = par$df
  )
  do.call(prec_change.default, .roundtrip_args(args, list(...), prec_change.default))
}
