#' Sampling precision for an estimate pooled across occasions
#'
#' Compute the sampling error (SE, margin of error, CV) for the equal-weight
#' mean of the occasion estimates of a repeated survey, given the size of
#' each occasion and how far consecutive occasions overlap. An annual
#' average from quarterly rounds is the ordinary case. This is the inverse
#' of [n_pooled()].
#'
#' @param var For the default method: the population variance \eqn{S^2} on
#'   one occasion, taken to be the same on each. Supply this or `p`, not
#'   both. For `svyplan_n` objects: a sample size result from [n_pooled()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param n Sample size per occasion, gross units drawn.
#' @param occasions Number of occasions entering the average, at least 2.
#' @param sd Population standard deviation on one occasion, an alternative
#'   spelling of `var`. Supply exactly one of `var` or `sd`.
#' @param p The proportion being averaged, as an alternative to `var`. The
#'   occasion variance is then \eqn{Np(1-p)/(N-1)}, the same finite
#'   population variance [prec_prop()] uses, and `mu` is determined and must
#'   not be supplied.
#' @param mu The level being averaged, on the `var` scale. Required only for
#'   the relative measures: `cv` and `rmoe` are defined against it and are
#'   `NA` without it. Determined by `p` on the proportion scale.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0), applied to the variance of
#'   the pooled estimate rather than to any one occasion. See Details.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). It applies to every occasion and nets each size down to
#'   `n * resp_rate` before the covariance is formed. It is a single round's
#'   response, not attrition across a panel.
#' @param overlap Fraction of one occasion's **responding** sample carried
#'   into a later one, in \[0, 1\]. One number, meaning the same overlap at
#'   every lag, or one per lag, so `overlap[1]` is the consecutive-occasion
#'   figure. A [design_overlap()] result is accepted directly at
#'   `resp_rate = 1`; see Details on why not below it. `0` (default) makes
#'   the occasions independent, and the pooled variance is then a single
#'   occasion's divided by `occasions`.
#' @param overlap_cor Correlation between two occasions among the units they
#'   share, in \[0, 1\]. One number for every lag, or one per lag. Default
#'   0, which makes overlap irrelevant, since it is the product
#'   `overlap * overlap_cor` that moves the variance. Supply this or
#'   `cor_decay`, not both.
#' @param cor_decay Correlation at lag \eqn{m} taken as
#'   \eqn{\rho^m}{rho^m}, an alternative to stating `overlap_cor` lag by
#'   lag. It is the shape most panels have, a correlation falling away with
#'   distance, and a single number is rarely the right one over a long
#'   horizon.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()]. It switches the interval quantile from normal to t.
#'   `NULL` (default) applies no adjustment.
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param .overlap_basis Internal. Records whether a stored `overlap`
#'   profile was resolved from a [design_overlap()] object (`"issued"`) or
#'   supplied as respondent overlap (`"respondent"`), so that a round trip
#'   or a grid meets the same refusal the first call would have. Set from
#'   the result being re-read; there is no reason to pass it by hand.
#'
#' @return A `svyplan_prec` object with `type = "pooled"`:
#' \describe{
#'   \item{`se`}{Standard error of the pooled estimate, computed on the net
#'     sizes `n * resp_rate`.}
#'   \item{`moe`}{Margin of error, `qnorm(1 - alpha / 2) * se`.}
#'   \item{`cv`}{Standard error relative to the level, `se / abs(mu)`. `NA`
#'     when `mu` is neither supplied nor determined by `p`.}
#'   \item{`params`}{The validated inputs. Dispersion is always stored as
#'     `var`, and `overlap` and `overlap_cor` as resolved vectors of one
#'     entry per lag, whatever form they were supplied in, alongside the
#'     `overlap_basis` they were resolved under.}
#' }
#'
#' Nothing here is rounded, so passing a continuous `n` back from
#' [n_pooled()] reproduces its `se`, `moe` and `cv` exactly.
#'
#' @details
#' The estimand is the equal-weight mean of the occasion estimates,
#' \eqn{\bar y = T^{-1}\sum_t \bar y_t}{ybar = (1/T) sum_t ybar_t}. Three
#' things share that name and this is none of the other two: it is not a
#' pooled variance, and it is not the pooling of the rotating cohorts that
#' make up one occasion, which [n_panel()] describes.
#'
#' Writing \eqn{T} for `occasions`, \eqn{n} for the net size,
#' \eqn{o_m} for the overlap at lag \eqn{m} and \eqn{\rho_m} for the
#' correlation there, the occasions carry
#'
#' \deqn{Var(\bar y_t) = \frac{S^2}{n}\Big(1 - \frac{n}{N}\Big), \qquad
#'       Cov(\bar y_t, \bar y_s) = \rho_m S^2
#'       \Big(\frac{o_m}{n} - \frac{1}{N}\Big),}{Var(ybar_t) = S^2/n (1 - n/N),   Cov(ybar_t, ybar_s) = rho_m S^2 (o_m/n - 1/N),}
#'
#' the same kernel [prec_change()] takes a difference on, and the average of
#' them has
#'
#' \deqn{V = \frac{1}{T^2}\Big[\,T\,Var(\bar y_t)
#'       + 2\sum_{m=1}^{T-1}(T-m)\,Cov_m\Big].}{V = (1/T^2) [ T Var(ybar_t) + 2 sum_m (T - m) Cov_m ].}
#'
#' A lag whose overlap is zero contributes exactly zero, which is the rule
#' [prec_change()] applies at `overlap = 0` read lag by lag. Two boundaries
#' fix the whole: a fresh sample each occasion pools to
#' \eqn{Var(\bar y_t)/T}, independent averaging, and a full panel measured
#' with correlation 1 pools to \eqn{Var(\bar y_t)}, the same units every
#' time and averaging buys nothing.
#'
#' `deff` multiplies the assembled variance. It is the design effect of the
#' pooled estimate, which is not in general any one occasion's.
#'
#' ## Overlap works against a pooled estimate, when it works at all
#'
#' A covariance is subtracted in a difference and added in a sum, so
#' whichever sign it carries, it moves a change and a pooled average in
#' opposite directions. **When the covariance is positive**, which is the
#' ordinary case, overlap improves a change and inflates a pooled average,
#' and planners routinely get that backwards.
#'
#' The covariance is positive exactly when
#' \eqn{o_m > n/N}{o_m > n/N}: shared units above what independent draws
#' from the same finite population would already give the two occasions.
#' Without a finite population correction any positive overlap qualifies,
#' so the ordinary case is the only case there. With one, an overlap
#' *below* the sampling fraction means the occasions are negatively
#' coordinated, sharing fewer units than chance, and **the directions
#' reverse**: at `N = 1000`, `n = 500`, `overlap = 0.25` and
#' `overlap_cor = 0.8`, two occasions have a pooled variance of 0.030
#' against 0.050 for independent ones, while the change variance rises from
#' 0.200 to 0.280. The kernel is a valid covariance throughout.
#'
#' The level at a single occasion is unaffected either way, being a function
#' of that occasion's size alone, so it is a reference line rather than a
#' third position: what trades off is the change against the pooled
#' estimate. `?prec_change` and this page size the two arms.
#'
#' ## Issued overlap and respondent overlap
#'
#' A bare number is the overlap between the **responding** samples, matching
#' [prec_change()]. [design_overlap()] reports the overlap between the
#' **issued** ones, and the two coincide only at full response, so its
#' result is accepted here at `resp_rate = 1` and refused below it. The
#' conversion needs an assumption about how response at one occasion depends
#' on response at another, which this function does not make: under
#' independent response the responding overlap is about `resp_rate` times
#' the issued one, and stating it is the planner's decision rather than this
#' function's.
#'
#' ## When the inputs do not describe a covariance
#'
#' A lag's covariance changes sign once its overlap falls below the sampling
#' fraction `n / N`, and a covariance whose sign varies across lags need not
#' be a valid covariance at all. This is checked on the assembled matrix,
#' because a kernel that fails it can still return a positive pooled
#' variance, and a design refused here is one whose stated overlap,
#' correlation and sampling fraction cannot hold together. High sampling
#' fractions are where it bites, from about `n / N = 0.6` upwards over a
#' long horizon, and much of that region is also a rotation that would
#' exhaust its own population; the message says so when it does.
#'
#' ## What the overlap covariance assumes
#'
#' As in [prec_change()], the covariance is a model form, a correlation
#' between occasions among the units they share and simple random sampling
#' otherwise. It is not the design-based covariance of a complex master
#' sample, which carries its own pairwise inclusion probabilities.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley. Chapter 12.
#'
#' @family precision functions
#' @seealso [n_pooled()] for the inverse (compute n from a precision
#'   target), [prec_change()] for the other arm of the trade-off,
#'   [design_overlap()] for the overlap a rotation schedule gives,
#'   [prec_mean()] and [prec_prop()] for a single occasion.
#'
#' @examples
#' # Four independent quarterly rounds averaged into an annual figure
#' prec_pooled(var = 100, n = 500, occasions = 4)
#'
#' # The same rounds from a rotating panel: overlap inflates the average
#' prec_pooled(var = 100, n = 500, occasions = 4, overlap = 0.75,
#'             cor_decay = 0.8)
#'
#' # An annual average of a proportion measured monthly
#' prec_pooled(p = 0.3, n = 1200, occasions = 12, overlap = 0.75,
#'             cor_decay = 0.9)
#'
#' # Overlap read off a rotation schedule, at full response
#' prec_pooled(var = 100, n = 500, occasions = 8,
#'             overlap = design_overlap("4", max_lag = 7), cor_decay = 0.8)
#'
#' # Round-trip from n_pooled
#' res <- n_pooled(var = 100, moe = 1, occasions = 4, overlap = 0.5,
#'                 cor_decay = 0.7)
#' prec_pooled(res)
#'
#' @export
prec_pooled <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", prec_pooled.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_pooled")
}

#' @rdname prec_pooled
#' @export
prec_pooled.default <- function(
  var = NULL,
  n,
  occasions,
  ...,
  sd = NULL,
  p = NULL,
  mu = NULL,
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
  .plan <- .merge_plan_args(plan, prec_pooled.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_pooled.default, c(.plan, list(...))))
  .check_unused_dots(...)
  check_population_size(N)
  .check_occasions(occasions)
  inputs <- .pooled_inputs(var, sd, p, mu, N)
  check_scalar(n, "n")
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)
  if (!is.null(df)) check_df(df)
  ov <- .resolve_lag_overlap(overlap, occasions, resp_rate, .overlap_basis)
  rho <- .resolve_lag_cor(overlap_cor, cor_decay, occasions - 1L)
  .check_gross_n(n, N, label = "each occasion")

  prec <- .prec_engine_pooled(
    n, inputs$var, inputs$mu, N, deff, resp_rate, occasions, ov$overlap, rho,
    alpha, df
  )

  params <- list(
    var = inputs$var,
    n = n,
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

  .new_svyplan_prec(
    se = prec$se,
    moe = prec$moe,
    cv = prec$cv,
    type = "pooled",
    params = params
  )
}

#' @rdname prec_pooled
#' @export
prec_pooled.svyplan_n <- function(var, ...) {
  x <- var
  if (x$type != "pooled") {
    stop("prec_pooled requires a svyplan_n of type 'pooled'", call. = FALSE)
  }
  par <- x$params
  args <- list(
    var = if (is.null(par$p)) par$var else NULL,
    p = par$p,
    n = x$n,
    occasions = par$occasions,
    mu = if (is.null(par$p)) par$mu else NULL,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate %||% 1,
    overlap = par$overlap,
    overlap_cor = par$overlap_cor,
    df = par$df,
    .overlap_basis = par$overlap_basis
  )
  do.call(prec_pooled.default, .roundtrip_args(args, list(...), prec_pooled.default))
}
