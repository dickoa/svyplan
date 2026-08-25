#' Sample size for a proportion
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
#'   Specify exactly one of `moe`, `cv`, or `rmoe`.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.10` means the standard error should be at
#'   most 10 percent of the estimate. Use `cv` when you want precision
#'   to scale with the estimate (common in economic surveys). Use `moe`
#'   when you want a fixed absolute precision (common in health surveys).
#'   Specify exactly one of `moe`, `cv`, or `rmoe`.
#' @param rmoe Target margin of error relative to `p`, so `rmoe = 0.12`
#'   asks for a 95 percent interval whose half-width is 12 percent of the
#'   proportion. It is `moe / p`, and therefore fixes the same interval `moe` does
#'   while scaling with the estimate the way `cv` does. See the precision
#'   quantities section of [prec_prop()]. Specify exactly one of `moe`,
#'   `cv`, or `rmoe`.
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
#'   population. See the nonresponse section of [svyplan-package] for what
#'   this adjustment does and does not claim.
#' @param method One of `"wald"` (default), `"wilson"`, `"logodds"`, or
#'   `"beta"`.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()] where the plan is in hand. It widens the interval for a
#'   variance estimated from few clusters. `NULL` (default) applies no
#'   adjustment, treating the variance as known. See Details, including for
#'   what an unset `df` means under `method = "beta"`, where the reference
#'   is the value a simple random sample would have rather than infinity.
#' @param min_cases Minimum expected number of positive cases the sample
#'   must yield, an alternative constraint to precision for a rare outcome.
#'   The returned `n` satisfies both. `NULL` (default) sizes on precision
#'   alone. See Details.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "proportion"`:
#' \describe{
#'   \item{`n`}{Required sample size, continuous and gross. It already
#'     carries `deff` and the `1 / resp_rate` inflation, so it counts the
#'     units to release, not the completed interviews. `$n` and
#'     `as.double()` keep the unrounded value, which is what makes the
#'     round trip through [prec_prop()] exact, whereas `print()` and
#'     `as.integer()` round it up to the whole units you would field.
#'     Take the field figure from `as.integer()` rather than from `$n`.}
#'   \item{`se`, `moe`, `cv`, `rmoe`}{Precision the design achieves at that `n`,
#'     the same values [prec_prop()] reports for the same inputs. `se` is
#'     the sampling standard error and `cv` is `se / p`, so both are the
#'     same under all four methods. `moe` is half the length of the
#'     interval the chosen `method` builds, which equals
#'     `q * se` under `"wald"` alone, with `q` defined in Details. Where
#'     [confint()] truncates a Wald interval at 0 or 1, `moe` stays
#'     `q * se` and exceeds half the reported width (see Details). For the
#'     three asymmetric methods it is not an offset from either limit, so read
#'     the limits with [confint()]. `rmoe` is `moe / p`, so a target
#'     stated as a relative margin of error reads back in the units it
#'     was stated in.}
#'   \item{`method`}{The interval method used.}
#'   \item{`expected_cases`}{Positive cases the design expects to yield,
#'     `n * resp_rate * p`.}
#'   \item{`binding`}{Which constraint set the size, `"precision"` or
#'     `"min_cases"`. `NULL` when no `min_cases` was given, there being
#'     nothing for precision to bind against.}
#'   \item{`params`}{The validated inputs (`p`, `alpha`, `N`, `deff`,
#'     `resp_rate`, `df`, `min_cases` when given, and whichever of `moe`,
#'     `cv`, or `rmoe` was the target). [predict()], [confint()] and the
#'     [prec_prop()] round trip read the design back from here.}
#' }
#'
#' @details
#' Four confidence interval methods are available:
#'
#' - **Wald** (`"wald"`): Standard normal approximation
#'   (Cochran, 1977, Ch. 3). Supports `moe`, `rmoe`, and `cv` targets,
#'   with optional finite population correction.
#' - **Wilson** (`"wilson"`): Wilson (1927) score interval. `moe` and
#'   `rmoe` targets only, with optional finite population correction.
#' - **Log-odds** (`"logodds"`): Log-odds (logit) transform interval.
#'   `moe` and `rmoe` targets only, with optional finite population
#'   correction.
#' - **Beta** (`"beta"`): Korn-Graubard (1998) interval, the Clopper-Pearson
#'   limits evaluated at the effective sample size, optionally widened for
#'   the degrees of freedom of the variance estimator via `df`. `moe` and
#'   `rmoe` targets only, with optional finite population correction. This
#'   is the method
#'   `survey::svyciprop(method = "beta")` reports, and the two agree exactly
#'   on a design where they see the same effective size.
#'
#' All four read one variance,
#' \eqn{\mathrm{deff}\,\frac{N}{N-1}\,p(1-p)(1/n-1/N)}{deff N/(N-1) p(1-p)(1/n-1/N)}. They differ only in
#' the interval built around it, so they agree closely whenever the margin
#' of error is small and diverge only where the normal approximation itself
#' is doubtful. In practice the choice matters for a rare or near-universal
#' outcome and is immaterial otherwise. See the section below.
#'
#' The design effect and the finite population correction enter every method
#' through the effective sample size
#' \eqn{n_\mathrm{eff}=n_\mathrm{net}/(\mathrm{deff}\cdot\mathrm{fpc})}{n_eff=n_net/(deff * fpc)}, the
#' size at which an infinite-population simple random sample would carry the
#' same variance. For Wald this reproduces the usual closed form exactly, and
#' for the other three it is the approximation that keeps `deff` and `N`
#' acting on the interval as they act on the variance. A census therefore
#' yields a zero margin of error under all four.
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
#' that also solves a `cv` target, and it is what survey texts and sample
#' size tables conventionally report. `moe` and `rmoe` targets are available
#' under all four.
#'
#' How much the choice costs depends on the sample size, and the four
#' methods converge on each other slowly. With `alpha = 0.05`, `deff = 1`,
#' `N = Inf`, `resp_rate = 1`, no supplied `df`, and `p` from 0.1 to 0.9,
#' the widest and narrowest margins of error differ by 8.2 percent at
#' `n = 100`, 3.8 percent at 500, 2.7 percent at 1000, and 1.2 percent at
#' 5000. Below a few thousand the choice is worth a moment. At larger survey
#' sizes its effect on required size is often modest over this range.
#'
#' The case that decides it is a rare or near-universal outcome. As `p`
#' approaches 0 or 1 the Wald interval loses coverage and can extend past 0
#' or 1. [confint()] truncates it at the boundary when it does, and `moe`
#' is then no longer half the reported interval's length. Wilson and
#' log-odds lie strictly inside \eqn{(0, 1)} and beta inside \eqn{[0, 1]},
#' all three by construction, so none of them needs that truncation.
#' Planning a survey for a 2 percent prevalence is the case that justifies
#' `method = "wilson"`, since the score interval keeps its nominal coverage
#' far better there.
#'
#' Log-odds is the narrower case again. It respects the parameter space like
#' Wilson but keeps the estimate at the center of the interval on the logit
#' scale. Back-transforming its limits gives an interval for the proportion,
#' or equivalently for its odds. It is not an interval for an odds ratio and
#' does not size a logistic regression. Prefer Wilson unless reporting on the
#' logit or odds scale is itself useful.
#'
#' Reach for `"beta"` when the expected number of outcomes in the category
#' nearest a boundary is small in absolute terms, not merely when `p` is
#' numerically close to a boundary. Rare-outcome domain estimates in a
#' clustered survey are the standard case. Near `p = 1`, apply the same
#' reasoning to expected negative outcomes. Beta is the widest over most
#' of the range, and so the most demanding to plan for, but not
#' everywhere: for a very
#' rare outcome at a small sample, log-odds is wider still. Nor is it the
#' only method that answers to a variance estimated from few clusters. All
#' four read `df` through the same t quantile and widen by it, so `df`
#' alone is not a reason to choose beta.
#'
#' A note on where the methods bind. For *sizing*, the sizes differ by the
#' measured few-percent amounts over the range above. That does not make the
#' interval methods inferentially interchangeable. When the reported
#' interval's width is the binding requirement, size under the method you
#' will report. Compare methods with separate [n_prop()] or [prec_prop()]
#' calls, or with separately constructed objects. [predict()] varies the
#' permitted planning parameters while keeping an object's method fixed. For
#' *assessing* an achieved design, method choice can dominate when realized
#' counts are sparse or the estimate lies near a boundary.
#'
#' ## A minimum expected number of cases
#'
#' For a rare outcome in a small domain the constraint that actually binds
#' is often not a margin of error but an expected count. An indicator may
#' require 30 positive cases on average, or a subgroup analysis may need a
#' planning floor for its expected events. `min_cases` states that requirement
#' directly. The size it demands is
#' \deqn{n_\mathrm{cases} = \mathrm{min\_cases} / (p \cdot
#'       \mathrm{resp\_rate}),}{n_cases = min_cases / (p * resp_rate),}
#' and the result is the larger of that and the size precision asks for.
#' Thus the count floor holds in expectation. The realized count is random,
#' and the calculation does not control the probability of observing at
#' least `min_cases`. Such an assurance constraint is outside this function's
#' scope. `binding` says which planning constraint decided, and
#' `expected_cases` reports the expected positive count either way.
#'
#' `deff` does not enter this size, and that is deliberate rather than an
#' omission. A design effect describes how precisely the proportion is
#' estimated and the number of positive cases that turn up in a sample of
#' \eqn{n} is a property of the sample size and the prevalence alone. The
#' response rate does enter, because the cases are counted among
#' respondents and the returned `n` is gross, on the same footing as every
#' other size the package reports.
#'
#' `expected_cases` is reported on every proportion result, with or without
#' `min_cases`. It is the expected number of positive respondents after the
#' response-rate adjustment, not an observed count or an assurance
#' probability. It does not incorporate `deff` and is not an effective event
#' count. It helps interpret a rare positive outcome, but it is not a complete
#' interval diagnostic. Near `p = 1`, inspect the expected negative count
#' instead. The `min_cases` constraint acts only on expected positive cases.
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
#' ## Degrees of freedom
#'
#' Write \eqn{q} for the interval quantile. By default
#' \eqn{q = \Phi^{-1}(1-\alpha/2)}{q = qnorm(1 - alpha/2)}, which treats the
#' variance as known. Supplying `df` instead gives
#' \eqn{q = t_{1-\alpha/2}(\mathrm{df})}{q = qt(1 - alpha/2, df)}.
#'
#' Supplying `df` says the variance will be estimated from a design with
#' that many degrees of freedom, and switches the quantile to
#' \eqn{t_{1-\alpha/2}(\mathrm{df})}{t_(1-alpha/2)(df)}. It is opt-in rather than derived: in
#' [n_cluster()] the df depends on the PSU count being solved for, so an
#' automatic version would need a fixed point. Read it off a plan you
#' already have with [design_df()], or pass a count directly.
#'
#' The three interval methods and the Korn-Graubard one apply that quantile
#' widening by different routes. `"wald"`, `"wilson"` and `"logodds"`
#' substitute the quantile in the half-width. `"beta"` instead scales the
#' effective sample size by the squared ratio of the two t quantiles
#' (Korn and Graubard, 1998, eq. 2.2), which is the same widening expressed
#' on the sample size rather than on the interval.
#'
#' The two routes differ in what an *unset* `df` means, because they differ
#' in what they measure it against. For the three substituting methods
#' `NULL` and `Inf` are the same statement, the normal quantile. For
#' `"beta"` the reference is the value a simple random sample of the same
#' size would have, so `NULL` is equivalent to `df = n - 1`, and `df = Inf`
#' claims more degrees of freedom than an SRS has and narrows the interval
#' accordingly.
#'
#' The reported `se` is the sampling standard error from the common variance
#' above, and `cv = se / p`. Neither depends on `method`, `alpha`, or `df`.
#' The margin of error reads the interval construction and is therefore
#' method-specific. Under `"wald"`, `moe = q * se`. The other three methods
#' are asymmetric, so their half-length has no universal `q * se` identity.
#'
#' `df` is deliberately absent from the power functions. There the
#' quantile is a normal deviate for an alternative rather than an interval
#' half-width, and a t-based power calculation is a different procedure,
#' not a substituted quantile.
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
#' @family sample size functions
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
#' # At least 30 expected cases of a rare outcome, whatever precision asks
#' rare <- n_prop(p = 0.02, moe = 0.02, min_cases = 30)
#' rare$binding
#' rare$expected_cases
#'
#' # MICS-style relative margin of error, 12% of the proportion
#' n_prop(p = 0.2, rmoe = 0.12, deff = 1.5, resp_rate = 0.9)
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
  rmoe = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  method = c("wald", "wilson", "logodds", "beta"),
  df = NULL,
  min_cases = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_prop.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_prop.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  check_proportion(p, "p")
  check_precision(moe, cv, rmoe)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  method <- match.arg(method)

  if (!is.null(df)) {
    check_df(df)
  }

  # Normalized once, here, so no solver or engine below ever sees 'rmoe'.
  moe_used <- if (is.null(rmoe)) moe else .moe_from_rmoe(rmoe, p, "p")

  n <- switch(
    method,
    wald = .n_prop_wald(p, moe_used, cv, alpha, N, deff, df),
    wilson = .n_prop_wilson(p, moe_used, cv, alpha, N, deff, df),
    logodds = .n_prop_logodds(p, moe_used, cv, alpha, N, deff, df),
    beta = .n_prop_beta(p, moe_used, cv, alpha, N, deff, df)
  )

  n <- .apply_resp_rate(n, resp_rate)

  binding <- NULL
  if (!is.null(min_cases)) {
    n_cases <- .n_from_min_cases(min_cases, p, resp_rate)
    binding <- if (n_cases > n) "min_cases" else "precision"
    n <- max(n, n_cases)
  }
  .check_attainable(n, N, resp_rate)

  params <- list(
    p = p,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  # Exactly one of the three, as supplied: a round trip through predict()
  # rebuilds the call from these, and two targets at once would fail the gate.
  if (!is.null(moe)) {
    params$moe <- moe
  } else if (!is.null(rmoe)) {
    params$rmoe <- rmoe
  } else {
    params$cv <- cv
  }
  if (!is.null(min_cases)) {
    params$min_cases <- min_cases
  }

  .new_svyplan_n(
    n = n,
    type = "proportion",
    method = method,
    params = params,
    binding = binding
  )
}

#' Effective SRS size a method needs to reach a margin of error
#'
#' The interval-method-specific content of a `moe` target, expressed as the
#' infinite-population simple random sample size that delivers it with no
#' design effect and full response. It is what lets a margin of error be
#' restated as a sampling CV without assuming the Wald half-width, since
#' `sqrt(p (1 - p) / n_eff)` is the standard error at which the method's own
#' interval closes to `moe`. Under `"wald"` it returns `z^2 p q / moe^2`, so
#' the restatement reduces to `moe / (z p)` exactly.
#' @keywords internal
#' @noRd
.n_prop_effective <- function(p, moe, alpha, method = "wald", df = NULL) {
  switch(
    method,
    wald = .n_prop_wald(p, moe, NULL, alpha, Inf, 1, df),
    wilson = .n_prop_wilson(p, moe, NULL, alpha, Inf, 1, df),
    logodds = .n_prop_logodds(p, moe, NULL, alpha, Inf, 1, df),
    beta = .n_prop_beta(p, moe, NULL, alpha, Inf, 1, df)
  )
}

#' @keywords internal
#' @noRd
.n_prop_wald <- function(p, moe, cv, alpha, N, deff = 1, df = NULL) {
  z <- .q_alpha(alpha, df)
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
      sprintf(
        "the %s margin of error is unattainable: it never exceeds %.4g",
        what,
        moe
      ),
      call. = FALSE
    )
  }
  if (is.null(upper)) {
    upper <- 1
    for (i in seq_len(200L)) {
      if (gap(upper) < 0) {
        break
      }
      upper <- upper * 4
    }
    if (gap(upper) >= 0) {
      stop(
        sprintf("the %s sample size search did not converge", what),
        call. = FALSE
      )
    }
  }
  uniroot(gap, c(lower, upper), tol = .Machine$double.eps^0.75)$root
}

#' @keywords internal
#' @noRd
.n_prop_wilson <- function(p, moe, cv, alpha, N, deff = 1, df = NULL) {
  if (is.null(moe)) {
    stop("Wilson method requires 'moe' (not 'cv')", call. = FALSE)
  }
  z <- .q_alpha(alpha, df)
  # The Wilson half-width tends to 1/2 as n tends to 0, so a wider margin
  # is never achieved by any sample size.
  if (moe >= 0.5) {
    stop("Wilson method requires 'moe' < 0.5", call. = FALSE)
  }
  # Solve on the effective scale, where the score interval is defined, then
  # convert back through the shared variance so that 'deff' and 'N' enter
  # exactly as they do for the Wald and log-odds methods.
  n_eff <- .solve_n_from_moe(
    function(n) .wilson_moe(p, n, z),
    moe,
    what = "Wilson"
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
.n_prop_logodds <- function(p, moe, cv, alpha, N, deff = 1, df = NULL) {
  if (is.null(moe)) {
    stop("Log-odds method requires 'moe' (not 'cv')", call. = FALSE)
  }
  if (moe >= 0.5) {
    stop("log-odds method requires 'moe' < 0.5", call. = FALSE)
  }
  .n_prop_logodds_raw(p, moe, alpha, N, deff, df)
}

#' Core log-odds n solver (no cv check).
#' Used by both .n_prop_logodds() and prec_prop()'s round trip.
#' @keywords internal
#' @noRd
.n_prop_logodds_raw <- function(p, e, alpha, N, deff = 1, df = NULL) {
  .solve_n_from_moe(
    function(n) .logodds_moe(p, n, alpha, N, deff, df),
    e,
    upper = if (is.infinite(N)) NULL else N * (1 - 1e-12),
    what = "log-odds"
  )
}

#' @rdname n_prop
#' @export
n_prop.svyplan_prec <- function(p, ..., moe = NULL, cv = NULL, rmoe = NULL) {
  x <- p
  if (x$type != "proportion") {
    stop("n_prop requires a svyplan_prec of type 'proportion'", call. = FALSE)
  }
  par <- x$params
  # The achieved margin of error is the implied target, but only when the
  # caller named none of the three: restoring it alongside an override
  # would send two targets into a function that takes one.
  if (is.null(moe) && is.null(cv) && is.null(rmoe)) {
    moe <- x$moe
  }
  args <- list(
    p = par$p,
    moe = moe,
    cv = cv,
    rmoe = rmoe,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    method = x$method %||% "wald",
    df = par$df,
    min_cases = par$min_cases
  )
  do.call(n_prop.default, .roundtrip_args(args, list(...), n_prop.default))
}
