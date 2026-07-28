#' Sample Size for a Proportion
#'
#' Compute the required sample size for estimating a population proportion
#' with a specified margin of error or coefficient of variation.
#'
#' @param p For the default method: expected proportion, in (0, 1).
#'   For `svyplan_prec` objects: a precision result from [prec_prop()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param moe Desired margin of error, the half-width of the confidence
#'   interval on the proportion scale. For example, `moe = 0.05` means
#'   the 95 percent CI should be no wider than +/- 5 percentage points.
#'   Specify exactly one of `moe` or `cv`.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.10` means the standard error should be at
#'   most 10 percent of the estimate. Use `cv` when you want precision
#'   to scale with the estimate (common in economic surveys). Use `moe`
#'   when you want a fixed absolute precision (common in health/DHS
#'   surveys). Specify exactly one of `moe` or `cv`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Accounts for the loss of
#'   precision from a complex design (clustering, unequal weights)
#'   compared to simple random sampling. A DEFF of 1.5 means 50 percent
#'   more interviews are needed for the same precision. Estimate from a
#'   previous survey, use [design_effect()] to compute it, or apply a
#'   rule of thumb (1.5--2.0 for typical cluster designs). Values < 1
#'   are valid for efficient designs (e.g., stratified sampling with
#'   Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The required sample size is inflated by `1 / resp_rate`.
#'   Estimate from response rates observed in similar surveys in the same
#'   population.
#' @param method One of `"wald"` (default), `"wilson"`, `"logodds"`, or
#'   `"beta"`.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, typically sampled PSUs minus strata. Used by
#'   `method = "beta"` only, where it widens the interval for a variance
#'   estimated from few clusters (Korn and Graubard, 1998, eq. 2.2). `NULL`
#'   (default) applies no adjustment, which is equivalent to `df = n - 1`,
#'   the value a simple random sample would have.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "proportion"`:
#' \describe{
#'   \item{`n`}{Required sample size, continuous and gross. It already
#'     carries `deff` and the `1 / resp_rate` inflation, so it counts the
#'     units to release, not the completed interviews. `$n` and
#'     `as.double()` keep the unrounded value, which is what makes the
#'     round trip through [prec_prop()] exact; `print()` and
#'     `as.integer()` round it up to the whole units you would field.
#'     Take the field figure from `as.integer()` rather than from `$n`.}
#'   \item{`se`, `moe`, `cv`}{Precision the design achieves at that `n`,
#'     the same values [prec_prop()] reports for the same inputs. `moe`
#'     is `qnorm(1 - alpha / 2) * se`; for the three asymmetric methods
#'     it is half the interval's length rather than an offset from either
#'     limit, so read the limits with [confint()].}
#'   \item{`method`}{The interval method used.}
#'   \item{`params`}{The validated inputs (`p`, `alpha`, `N`, `deff`,
#'     `resp_rate`, `df`, and whichever of `moe` or `cv` was the target).
#'     [predict()], [confint()] and the [prec_prop()] round trip read the
#'     design back from here.}
#' }
#'
#' @details
#' Four confidence interval methods are available:
#'
#' - **Wald** (`"wald"`): Standard normal approximation
#'   (Cochran, 1977, Ch. 3). Supports both `moe` and `cv` modes,
#'   with optional finite population correction.
#' - **Wilson** (`"wilson"`): Wilson (1927) score interval. Only `moe`
#'   mode, with optional finite population correction.
#' - **Log-odds** (`"logodds"`): Log-odds (logit) transform interval.
#'   Only `moe` mode, with optional finite population correction.
#' - **Beta** (`"beta"`): Korn-Graubard (1998) interval, the Clopper-Pearson
#'   limits evaluated at the effective sample size, optionally widened for
#'   the degrees of freedom of the variance estimator via `df`. Only `moe`
#'   mode, with optional finite population correction. This is the method
#'   `survey::svyciprop(method = "beta")` reports, and the two agree exactly
#'   on a design where they see the same effective size.
#'
#' All four read one variance,
#' \eqn{\mathrm{deff}\,\frac{N}{N-1}\,p(1-p)(1/n-1/N)}. They differ only in
#' the interval built around it, so they agree closely whenever the margin
#' of error is small and diverge only where the normal approximation itself
#' is doubtful. In practice the choice matters for a rare or near-universal
#' outcome and is immaterial otherwise; see the section below.
#'
#' The design effect and the finite population correction enter every method
#' through the effective sample size
#' \eqn{n_\mathrm{eff}=n_\mathrm{net}/(\mathrm{deff}\cdot\mathrm{fpc})}, the
#' size at which an infinite-population simple random sample would carry the
#' same variance. For Wald this reproduces the usual closed form exactly; for
#' the other three it is the approximation that keeps `deff` and `N` acting
#' on the interval as they act on the variance. A census therefore yields a
#' zero margin of error under all four.
#'
#' The Wilson, log-odds, and beta intervals are not symmetric about `p`, so
#' the reported `moe` is half the interval's length rather than an offset
#' either limit sits at. For the limits themselves, use [confint()], which
#' returns the interval the chosen method actually produces. The gap matters
#' most for the beta method at a rare outcome, where the upper arm can be
#' several times the lower one.
#'
#' ## Choosing a method
#'
#' Reach for Wald unless you have a reason not to. It is the only method
#' that also solves a `cv` target, it is what every survey text and every
#' sample size table reports, and for a proportion between roughly 0.1 and
#' 0.9 at any usable sample size the four methods differ by less than a
#' percent of the margin of error.
#'
#' The exception is a rare or near-universal outcome. As `p` approaches 0 or
#' 1 the Wald interval loses coverage and can extend past 0 or 1, while the
#' score interval stays inside the parameter space and keeps its nominal
#' coverage far better. Planning a survey for a 2 percent prevalence is the
#' case that justifies `method = "wilson"`.
#'
#' Log-odds is the narrower case again: it respects the parameter space like
#' Wilson but keeps the estimate at the centre of the interval on the logit
#' scale, which matters when the plan will be reported as an odds ratio or
#' fed into a logistic model. If you are not doing either, prefer Wilson.
#'
#' Reach for `"beta"` when the expected number of positive cases is small in
#' absolute terms, not merely when `p` is small: a few dozen cases or fewer,
#' which is the regime Korn and Graubard wrote for and the one where the
#' normality of the estimated proportion breaks down however large the
#' sample is. Rare-outcome domain estimates in a clustered survey are the
#' standard case. It is the most conservative of the four, it is the only
#' one whose interval is guaranteed to stay inside \eqn{[0, 1]} by
#' construction rather than by clamping, and it is the only one that can
#' account for a variance estimated from few clusters, through `df`. Its
#' cost is a larger planned sample.
#'
#' A note on where the methods bind. For *sizing*, the choice is close to
#' immaterial: across the usual range the four sizes differ by a few
#' percent, far less than the uncertainty in the assumed `p`, `deff`, and
#' response rate. For *assessing* an achieved design the choice can dominate,
#' because that is where small realized samples of rare outcomes appear. If
#' you are unsure, plan with Wald and report with `"beta"`.
#'
#' ## Finite population correction
#'
#' Setting `N` to a finite value reduces the required sample size when
#' the sampling fraction (n/N) is non-negligible. As a rule of thumb,
#' FPC has little effect when n/N < 5 percent. The Wald FPC uses the
#' Cochran (1977, Ch. 3) form with an `N/(N-1)` factor to account for
#' the Bernoulli finite-population variance. This differs from
#' [n_mean()], where no `N/(N-1)` adjustment is needed because the
#' variance is already defined on `N-1` degrees of freedom.
#'
#' All methods use the normal (z) quantile. This is standard for survey
#' sampling where the sample size is large enough for the CLT to apply.
#'
#' When called on a `svyplan_prec` object, parameters are extracted from the
#' stored result. Any argument of the default method (e.g. `method`, `deff`,
#' `N`) can be overridden through `...`. Unknown argument names are an
#' error. Passing a different `method` evaluates the stored precision
#' target under that formula. The round-trip will not be exact because the
#' precision was computed under the original method.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession,
#' and statistical inference. *Journal of the American Statistical
#' Association*, 22(158), 209--212.
#'
#' Brown, L. D., Cai, T. T. and DasGupta, A. (2001). Interval estimation
#' for a binomial proportion. *Statistical Science*, 16(2), 101--133.
#' The coverage comparison that motivates preferring the score interval
#' for a rare or near-universal outcome.
#'
#' Korn, E. L. and Graubard, B. I. (1998). Confidence intervals for
#' proportions with small expected number of positive counts estimated from
#' survey data. *Survey Methodology*, 24(2), 193--201.
#'
#' Clopper, C. J. and Pearson, E. S. (1934). The use of confidence or
#' fiducial limits illustrated in the case of the binomial.
#' *Biometrika*, 26(4), 404--413.
#'
#' @seealso [n_mean()] for continuous variables, [n_cluster()] for
#'   multistage designs, [n_multi()] for multiple indicators,
#'   [prec_prop()] for the inverse.
#'
#' @examples
#' # Wald, absolute margin of error
#' n_prop(p = 0.3, moe = 0.05)
#'
#' # Wald, target CV with finite population
#' n_prop(p = 0.5, cv = 0.10, N = 10000)
#'
#' # Wilson score interval
#' n_prop(p = 0.1, moe = 0.03, method = "wilson")
#'
#' # Korn-Graubard interval for a rare outcome, and the asymmetric limits
#' # it implies
#' rare <- n_prop(p = 0.02, moe = 0.01, method = "beta")
#' confint(rare)
#'
#' # The same, allowing for a variance estimated from 30 PSUs in 5 strata
#' n_prop(p = 0.02, moe = 0.01, method = "beta", df = 25)
#'
#' # With design effect and response rate
#' n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8)
#'
#' # MICS/DHS-style relative margin of error (RME)
#' # RME = moe / p, so moe = RME * p
#' p <- 0.2
#' n_prop(p = p, moe = 0.12 * p, deff = 1.5, resp_rate = 0.9)
#'
#' @export
n_prop <- function(p, ...) {
  if (!missing(p)) {
    .res <- .dispatch_plan(p, "p", n_prop.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_prop")
}

#' @rdname n_prop
#' @export
n_prop.default <- function(
  p,
  ...,
  moe = NULL,
  cv = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  method = c("wald", "wilson", "logodds", "beta"),
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_prop.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_prop.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  check_proportion(p, "p")
  check_precision(moe, cv)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  method <- match.arg(method)

  if (!is.null(df) && method != "beta") {
    stop("'df' applies to method = 'beta' only", call. = FALSE)
  }

  n <- switch(
    method,
    wald = .n_prop_wald(p, moe, cv, alpha, N, deff),
    wilson = .n_prop_wilson(p, moe, cv, alpha, N, deff),
    logodds = .n_prop_logodds(p, moe, cv, alpha, N, deff),
    beta = .n_prop_beta(p, moe, cv, alpha, N, deff, df)
  )

  n <- .apply_resp_rate(n, resp_rate)
  .check_attainable(n, N, resp_rate)

  params <- list(
    p = p,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  if (!is.null(moe)) {
    params$moe <- moe
  } else {
    params$cv <- cv
  }

  .new_svyplan_n(
    n = n,
    type = "proportion",
    method = method,
    params = params
  )
}

#' @keywords internal
#' @noRd
.n_prop_wald <- function(p, moe, cv, alpha, N, deff = 1) {
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p
  a <- ifelse(is.infinite(N), 1, N / (N - 1))

  if (!is.null(moe)) {
    n0 <- a * z^2 * deff * p * q / (moe^2 + z^2 * deff * p * q / (N - 1))
  } else {
    n0 <- a * deff * q / p / (cv^2 + deff * q / p / (N - 1))
  }

  n0
}

#' Invert a margin of error into a sample size
#'
#' Both non-Wald methods define the margin of error as the half-width of a
#' confidence interval that shrinks monotonically in the sample size. The
#' sample size is therefore recovered by solving `half_width(n) = moe`
#' rather than by carrying a separate closed form, which keeps each `n_*`
#' method an exact inverse of the corresponding `prec_*` method.
#'
#' `upper` bounds the search: the population size when it is finite,
#' otherwise `NULL` to expand the bracket until the half-width falls below
#' the target.
#' @keywords internal
#' @noRd
.solve_n_from_moe <- function(half_width, moe, upper = NULL, what) {
  gap <- function(n) half_width(n) - moe
  lower <- 1e-8
  if (gap(lower) <= 0) {
    stop(
      sprintf("the %s margin of error is unattainable: it never exceeds %.4g",
              what, moe),
      call. = FALSE
    )
  }
  if (is.null(upper)) {
    upper <- 1
    for (i in seq_len(200L)) {
      if (gap(upper) < 0) break
      upper <- upper * 4
    }
    if (gap(upper) >= 0) {
      stop(sprintf("the %s sample size search did not converge", what),
           call. = FALSE)
    }
  }
  uniroot(gap, c(lower, upper), tol = .Machine$double.eps^0.75)$root
}

#' @keywords internal
#' @noRd
.n_prop_wilson <- function(p, moe, cv, alpha, N, deff = 1) {
  if (is.null(moe)) {
    stop("Wilson method requires 'moe' (not 'cv')", call. = FALSE)
  }
  z <- qnorm(1 - alpha / 2)
  # The Wilson half-width tends to 1/2 as n tends to 0, so a wider margin
  # is never achieved by any sample size.
  if (moe >= 0.5) {
    stop("Wilson method requires 'moe' < 0.5", call. = FALSE)
  }
  # Solve on the effective scale, where the score interval is defined, then
  # convert back through the shared variance so that 'deff' and 'N' enter
  # exactly as they do for the Wald and log-odds methods.
  n_eff <- .solve_n_from_moe(
    function(n) .wilson_moe(p, n, z), moe, what = "Wilson"
  )
  .n_from_effective(n_eff, N, deff)
}

#' Korn-Graubard sample size
#'
#' Unlike the other methods, the degrees-of-freedom adjustment depends on the
#' net size itself, so the half-width is inverted directly on `n_net` rather
#' than on the effective scale.
#' @keywords internal
#' @noRd
.n_prop_beta <- function(p, moe, cv, alpha, N, deff = 1, df = NULL) {
  if (is.null(moe)) {
    stop("beta method requires 'moe' (not 'cv')", call. = FALSE)
  }
  if (moe >= 0.5) {
    stop("beta method requires 'moe' < 0.5", call. = FALSE)
  }
  .solve_n_from_moe(
    function(n_net) {
      .beta_moe(
        p,
        .kg_effective(.effective_from_n(n_net, N, deff), n_net, alpha, df),
        alpha
      )
    },
    moe,
    upper = if (is.infinite(N)) NULL else N * (1 - 1e-12),
    what = "Korn-Graubard"
  )
}

#' @keywords internal
#' @noRd
.n_prop_logodds <- function(p, moe, cv, alpha, N, deff = 1) {
  if (is.null(moe)) {
    stop("Log-odds method requires 'moe' (not 'cv')", call. = FALSE)
  }
  if (moe >= 0.5) {
    stop("log-odds method requires 'moe' < 0.5", call. = FALSE)
  }
  .n_prop_logodds_raw(p, moe, alpha, N, deff)
}

#' Core log-odds n solver (no cv check).
#' Used by both .n_prop_logodds() and prec_prop()'s round trip.
#' @keywords internal
#' @noRd
.n_prop_logodds_raw <- function(p, e, alpha, N, deff = 1) {
  .solve_n_from_moe(
    function(n) .logodds_moe(p, n, alpha, N, deff),
    e,
    upper = if (is.infinite(N)) NULL else N * (1 - 1e-12),
    what = "log-odds"
  )
}

#' @rdname n_prop
#' @export
n_prop.svyplan_prec <- function(p, ..., moe = NULL, cv = NULL) {
  x <- p
  if (x$type != "proportion") {
    stop("n_prop requires a svyplan_prec of type 'proportion'", call. = FALSE)
  }
  par <- x$params
  if (is.null(moe) && is.null(cv)) {
    moe <- x$moe
  }
  args <- list(
    p = par$p,
    moe = moe,
    cv = cv,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    method = x$method %||% "wald",
    df = par$df
  )
  do.call(n_prop.default, .roundtrip_args(args, list(...), n_prop.default))
}
