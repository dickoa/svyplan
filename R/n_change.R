#' Sample size for a change between two occasions
#'
#' Compute the sample size per occasion required to estimate the change in a
#' mean or a proportion between two occasions of the same population with a
#' specified margin of error or coefficient of variation, given how far the
#' two samples overlap.
#'
#' @param var For the default method: the population variance \eqn{S^2} on
#'   each occasion, one value for both or one per occasion. Estimate from a
#'   pilot study, a previous round, or published data. Supply this or `p`,
#'   not both. For `svyplan_prec` objects: a precision result from
#'   [prec_change()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param sd Population standard deviation on each occasion, an alternative
#'   spelling of `var`. Supply exactly one of `var` or `sd`.
#' @param p The two occasion proportions, `c(p1, p2)`, as an alternative to
#'   `var`. The occasion variances are then \eqn{Np(1-p)/(N-1)}, the same
#'   finite population variance [n_prop()] uses, and the change is
#'   \eqn{p_2 - p_1}, so `change` is determined and must not be supplied.
#'   Sizing for a change in a proportion is the ordinary case for a repeated
#'   household survey tracking an indicator.
#' @param change Expected change, \eqn{\mu_2 - \mu_1}, on the `var` scale. It
#'   may be negative. Required when `cv` or `rmoe` is specified, both being
#'   defined against it. Determined by `p` on the proportion scale.
#' @param moe Desired margin of error on the change, the half-width of its
#'   confidence interval, in the units the change is measured in. Specify
#'   exactly one of `moe`, `cv`, or `rmoe`.
#' @param cv Target standard error relative to the change, so `cv = 0.10`
#'   asks for a standard error one tenth of the change being measured.
#'   Requires `change`. Specify exactly one of `moe`, `cv`, or `rmoe`.
#' @param rmoe Target margin of error relative to the change, so
#'   `rmoe = 0.25` asks for a 95 percent interval whose half-width is a
#'   quarter of the change. It is `moe / abs(change)` and requires `change`.
#'   Specify exactly one of `moe`, `cv`, or `rmoe`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction. One value covers both occasions; two are accepted only at
#'   `overlap = 0`, where the occasions are independent and may legitimately
#'   be different populations. A positive `overlap` requires a single `N`,
#'   since units can only be shared by samples drawn from one population.
#' @param deff Design effect multiplier (> 0), applied to the variance of
#'   the change rather than to either occasion separately. See Details.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). Each occasion's size is inflated by `1 / resp_rate`. It is
#'   a single round's response, not attrition across a panel.
#' @param ratio Size of the first occasion relative to the second, so
#'   `ratio = 2` sizes a large baseline against a smaller follow-up. Default
#'   1, equal occasions, which returns a single size. Because `overlap` is
#'   measured against the first occasion it cannot exceed `1 / ratio` once
#'   the first occasion is the larger.
#' @param overlap Fraction of the first occasion's **responding** sample
#'   carried into the second, in \[0, 1\]. `0` (default) drops the covariance
#'   entirely and reproduces the two-independent-samples size. `1` means
#'   every responding unit of the first occasion is measured again, which is
#'   a full panel when the occasions are the same size. A rotation pattern
#'   determines the issued figure, `(D - 1) / D` between consecutive
#'   occasions of an in-for-`D` schedule, and [design_overlap()] computes it
#'   at every lag. The two are the same number at full response only. See
#'   [prec_change()] on what a response rate below 1 does to this reading.
#' @param overlap_cor Correlation between the two occasions among the
#'   overlapping units, in \[0, 1\]. Default 0, which makes overlap
#'   worthless: it is the product `overlap * overlap_cor` that reduces the
#'   size, so a full panel of uncorrelated measurements saves nothing. On the
#'   proportion scale two Bernoulli marginals bound the correlation they can
#'   have, and a value above that bound is rejected.
#' @param df Degrees of freedom of the variance estimator the planned design
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()]. It switches the interval quantile from normal to t.
#'   `NULL` (default) applies no adjustment.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "change"`:
#' \describe{
#'   \item{`n`}{Required size per occasion, continuous and gross. It already
#'     carries `deff` and the `1 / resp_rate` inflation, so it counts the
#'     units to release on each occasion, not the completed interviews. A
#'     single value when `ratio = 1`, otherwise one per occasion. `$n` and
#'     `as.double()` keep the unrounded value, which is what makes the round
#'     trip through [prec_change()] exact; `print()` rounds up to the whole
#'     units you would field. A positive overlap means the occasions are not
#'     disjoint, so the number of distinct units sampled is less than the
#'     sum; the number of interviews is not, since a shared unit is
#'     interviewed on both occasions.}
#'   \item{`se`, `moe`, `cv`, `rmoe`}{Precision the design achieves at that
#'     size, the same values [prec_change()] reports for the same inputs.
#'     `cv` and `rmoe` are `NA` unless the change is known.}
#'   \item{`params`}{The validated inputs, including whichever of `moe`,
#'     `cv`, or `rmoe` was the target. Dispersion is always stored as `var`,
#'     a pair; `p` is kept as well when the proportion scale was used.}
#' }
#'
#' @details
#' The change variance is \eqn{A/n_2 + B} in the net second-occasion size,
#' where
#'
#' \deqn{A = v_1/r + v_2 - 2\,\mathrm{overlap}\,\rho\sqrt{v_1v_2}, \qquad
#'       B = -\frac{v_1 + v_2 - 2\rho\sqrt{v_1v_2}}{N},}{A = v_1/r + v_2 - 2 overlap rho sqrt(v_1v_2), B = -(v_1 + v_2 - 2 rho sqrt(v_1v_2))/N,}
#'
#' with \eqn{r} the `ratio` and \eqn{\rho} the `overlap_cor`. The finite
#' population terms collect into \eqn{B}, a constant no sample size can
#' move, so the target inverts in closed form rather than by search. Both
#' coefficients are written above for a positive `overlap`; at
#' `overlap = 0` the covariance is dropped entirely, so
#' \eqn{A = v_1/r + v_2} and \eqn{B = -v_1/N_1 - v_2/N_2}, which is what
#' lets the two occasions come from different populations. See
#' [prec_change()] for the variance itself and for the piecewise statement.
#'
#' Overlap and correlation only ever act together. At `overlap_cor = 0` the
#' size is what two independent samples would need whatever the overlap, and
#' the two arguments are worth setting from the same evidence: a rotation
#' pattern fixes `overlap`, and a previous round of the same survey is what
#' identifies `overlap_cor`.
#'
#' The design effect is that of the change, not of either occasion. A
#' clustered design revisiting the same clusters has a smaller one than a
#' design that reclusters between occasions, because the cluster effects
#' partly cancel in the difference. Passing an occasion's design effect
#' here therefore over-sizes a design that keeps its clusters.
#'
#' ## When no size is enough, and when any size is
#'
#' Two boundaries are worth knowing before reading a result. Full overlap
#' with unit correlation at equal sizes and variances leaves the change with
#' no sampling variance at all: the same units are measured twice, so no
#' size is identified and the function says so rather than returning zero.
#'
#' At a finite `N` with equal occasions and full response, the requirement
#' rises towards `N` but never past it, because two censuses of one
#' population measure the change exactly. Every target is attainable there,
#' and a demanding one simply returns a near-census. Unequal occasions or a
#' response rate below 1 remove that guarantee, since the size released is
#' then larger than the sample carrying the precision, and a target beyond
#' reach is reported as unattainable.
#'
#' ## Sizing for a change, a level, or an average
#'
#' Overlap improves a change and does nothing for the level at a single
#' occasion, which is a function of that occasion's size alone. A fresh
#' sample each round is not the reverse of this: it does nothing for the
#' level either, and the two designs give the same single-occasion
#' precision at the same size. What overlap does work against is an
#' estimate pooled across occasions, an annual average of quarterly rounds
#' say, because the covariance it induces is subtracted in a difference and
#' added in a sum. That holds while the covariance is positive, which needs
#' the overlap to exceed the sampling fraction when `N` is finite; see
#' [prec_pooled()] for the boundary. A design serving more than one of the
#' three is sized by
#' taking the largest of `n_change()`, [n_pooled()] and the corresponding
#' [n_mean()] or [n_prop()] requirement, since none dominates.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley. Chapter 12.
#'
#' @family sample size functions
#' @seealso [prec_change()] for the inverse (compute precision from a size),
#'   [n_mean()] and [n_prop()] for a single occasion, [power_mean()] and
#'   [power_prop()] to frame the same overlap as a hypothesis test.
#'
#' @examples
#' # Two independent occasions
#' n_change(var = 100, moe = 2)
#'
#' # A half-overlapping panel correlated 0.6 needs fewer units per occasion
#' n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6)
#'
#' # A 6-point rise in a proportion, resolved to within 2 points
#' n_change(p = c(0.30, 0.36), moe = 0.02)
#'
#' # The same target under a consecutive-month CPS-style overlap of 0.75
#' n_change(p = c(0.30, 0.36), moe = 0.02, overlap = 0.75, overlap_cor = 0.5)
#'
#' # Relative target: resolve the change to within a quarter of itself
#' n_change(var = 100, change = 5, rmoe = 0.25)
#'
#' # A large baseline against a smaller follow-up
#' n_change(var = 100, moe = 2, ratio = 2)
#'
#' # With FPC, design effect, and response
#' n_change(var = 100, moe = 2, N = 20000, deff = 1.5, resp_rate = 0.8)
#'
#' @export
n_change <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", n_change.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_change")
}

#' @rdname n_change
#' @export
n_change.default <- function(
  var = NULL,
  ...,
  sd = NULL,
  p = NULL,
  change = NULL,
  moe = NULL,
  cv = NULL,
  rmoe = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  ratio = 1,
  overlap = 0,
  overlap_cor = 0,
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_change.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_change.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  N_pair <- .check_power_N(N)
  inputs <- .change_inputs(var, sd, p, change, N_pair)
  check_precision(moe, cv, rmoe)
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)
  check_scalar(ratio, "ratio")
  check_overlap(overlap)
  check_overlap_cor(overlap_cor)
  .check_overlap_N(overlap, N_pair)
  .check_overlap_n(overlap, ratio = ratio)
  .check_bernoulli_cor(inputs$p, overlap_cor, overlap)
  if (!is.null(df)) check_df(df)

  q <- .q_alpha(alpha, df)
  if (!is.null(cv)) {
    if (is.null(inputs$change)) {
      stop("'change' is required when 'cv' is specified", call. = FALSE)
    }
    if (inputs$change == 0) {
      stop("'cv' is undefined at 'change' = 0", call. = FALSE)
    }
    target_se <- cv * abs(inputs$change)
  } else {
    moe_used <- if (is.null(rmoe)) {
      moe
    } else {
      .moe_from_rmoe(rmoe, inputs$change, "change")
    }
    target_se <- moe_used / q
  }

  n <- .n_change_from_se(
    target_se, inputs$var_pair, N_pair, deff, resp_rate,
    ratio, overlap, overlap_cor
  )

  n_vec <- if (length(n) == 1L) c(n, n) else n
  .check_attainable(n_vec[1L], N_pair[1L], resp_rate)
  .check_attainable(n_vec[2L], N_pair[2L], resp_rate)

  prec <- .prec_engine_change(
    n, inputs$var_pair, inputs$change, N_pair, deff, resp_rate,
    overlap, overlap_cor, alpha, df
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
    var = inputs$var_pair,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    ratio = ratio,
    overlap = overlap,
    overlap_cor = overlap_cor,
    df = df
  )
  if (!is.null(inputs$p)) params$p <- inputs$p
  if (!is.null(inputs$change)) params$change <- inputs$change
  if (!is.null(moe)) {
    params$moe <- moe
  } else if (!is.null(rmoe)) {
    params$rmoe <- rmoe
  } else {
    params$cv <- cv
  }

  .new_svyplan_n(
    n = n,
    type = "change",
    params = params,
    se = prec$se,
    moe = prec$moe,
    cv = prec$cv
  )
}

#' @rdname n_change
#' @export
n_change.svyplan_prec <- function(var, ..., moe = NULL, cv = NULL,
                                  rmoe = NULL) {
  x <- var
  if (x$type != "change") {
    stop("n_change requires a svyplan_prec of type 'change'", call. = FALSE)
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
    change = if (is.null(par$p)) par$change else NULL,
    moe = moe,
    cv = cv,
    rmoe = rmoe,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    ratio = if (length(par$n) == 2L) par$n[1L] / par$n[2L] else 1,
    overlap = par$overlap,
    overlap_cor = par$overlap_cor,
    df = par$df
  )
  do.call(n_change.default, .roundtrip_args(args, list(...), n_change.default))
}
