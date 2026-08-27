#' Power analysis for proportions
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-sample test of proportions. Leave exactly one of `n`, `power`, or
#' `p2` as `NULL` to solve for that quantity.
#'
#' @param p1 Baseline proportion, in (0, 1).
#' @param p2 Alternative proportion, in (0, 1). Leave `NULL` to solve for
#'   the minimum detectable effect (MDE). The solver searches both above and
#'   below `p1` and returns the alternative closest to `p1` that achieves the
#'   target power. When `p1` is near 0 or 1, the MDE may only be detectable
#'   in one direction.
#' @param n Per-group sample size, measured as gross units drawn and bounded
#'   by the corresponding finite `N`. Scalar (equal groups) or length-2 vector
#'   `c(n1, n2)` for unequal groups. Leave `NULL` to solve for sample size.
#' @param power Target power, in (0, 1). When solving for sample size or MDE,
#'   it must exceed `alpha`, the power at zero effect. Leave `NULL` to solve
#'   for power.
#' @param alpha Significance level, default 0.05.
#' @param N Population size for finite-population correction. A scalar applies
#'   to both groups. A length-2 vector `c(N1, N2)` sets group-specific
#'   population sizes. `Inf` disables FPC for the corresponding group.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The required sample size is inflated by `1 / resp_rate`.
#' @param alternative Character: `"two.sided"` (default) or `"one.sided"`.
#'   Use `"two.sided"` when either an increase or decrease matters (the
#'   usual default). Use `"one.sided"` when only one direction is of
#'   interest (e.g. "has the intervention *reduced* stunting?"), which
#'   requires a smaller sample for the same power.
#' @param ratio Allocation ratio n1/n2 (default 1). Only used when solving
#'   for n (`n = NULL`). For example, `ratio = 2` means group 1 gets twice
#'   the sample of group 2.
#' @param overlap Panel overlap fraction in \[0, 1\], for repeated surveys.
#'   Defined as the fraction of group 1 that also appears in group 2
#'   (`overlap = n12 / n1`). Only supported with `method = "wald"`.
#' @param overlap_cor Correlation between occasions in \[0, 1\].
#' @param method Variance method: `"wald"` (default), `"arcsine"`, or
#'   `"logodds"`. The arcsine square-root transform is variance-stabilizing.
#'   The log-odds transform is not, because its variance still depends on the
#'   proportion. Both provide alternatives to the Wald calculation near a
#'   boundary (Valliant, 2018, sections 4.3.4--4.3.5).
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#'
#' @return A `svyplan_power` object with components:
#' \describe{
#'   \item{`n`}{Per-group sample size (scalar or length-2 for unequal groups).}
#'   \item{`power`}{Achieved power.}
#'   \item{`effect`}{Difference in proportions (`abs(p2 - p1)`).}
#'   \item{`solved`}{Which quantity was solved for
#'     (`"n"`, `"power"`, or `"mde"`).}
#'   \item{`params`}{List of input parameters.}
#' }
#'
#' @details
#' ## Choosing a method
#'
#' \describe{
#'   \item{`"wald"` (default)}{Standard Wald test. Appropriate for
#'     proportions in the 0.2--0.8 range. The only method that
#'     supports panel `overlap`.}
#'   \item{`"arcsine"`}{Arcsine-transformed test. The arcsine variance
#'     \eqn{1/(4n)} is approximately constant in \eqn{p}, making it
#'     more accurate when either proportion is below 0.15 or above 0.85.}
#'   \item{`"logodds"`}{Log-odds (logit) test with separate null and
#'     alternative variances. Also suitable for rare or extreme
#'     proportions. Uses Valliant Eq (4.23)/(4.24).}
#' }
#'
#' For proportions in the 0.2--0.8 range, all three methods give similar
#' results. Near a boundary, the arcsine and log-odds calculations avoid
#' relying on the untransformed Wald scale. Their operating characteristics
#' still depend on the design and the planning values supplied.
#'
#' ## Null variance convention
#'
#' `method = "wald"` uses the unpooled variance
#' \eqn{p_1 q_1 / n_1 + p_2 q_2 / n_2} for both the critical value and the
#' power shift. [stats::power.prop.test()] instead evaluates the critical
#' value under the null, using the pooled variance
#' \eqn{\bar{p} \bar{q} (1/n_1 + 1/n_2)}{pbar qbar (1/n_1 + 1/n_2)} with
#' \eqn{\bar{p} = (p_1 + p_2) / 2}{pbar = (p_1 + p_2) / 2}. Both conventions are standard, and they
#' differ by a few tenths of a percent in the resulting size:
#'
#' ```
#' power_prop(p1 = 0.3, p2 = 0.4, power = 0.8)$n   # 353.2
#' power.prop.test(p1 = 0.3, p2 = 0.4, power = 0.8)$n   # 355.9
#' ```
#'
#' The unpooled form is used here because it is the variance the package
#' reports everywhere else. It is the same expression that `deff`, the
#' finite population correction, and `resp_rate` act on in [n_prop()] and
#' [prec_prop()], so a power calculation and a precision calculation for the
#' same design stay on one scale. The pooled form has no finite-population
#' analogue that keeps that correspondence.
#'
#' `method = "logodds"` does carry a pooled null, because the statistic it
#' powers has one. It refers the log-odds difference to a critical value
#' computed under the null and a power shift computed under the alternative,
#' rejecting when
#'
#' \deqn{|\mathrm{logit}(\hat p_1) - \mathrm{logit}(\hat p_2)|
#'   > z_{\alpha} \sqrt{V_0},}{|logit(phat_1) - logit(phat_2)| > z_alpha sqrt(V_0),}
#'
#' with the two variances
#'
#' \deqn{V_0 = d \left( \frac{f_1}{n_1 \bar{p} \bar{q}}
#'                    + \frac{f_2}{n_2 \bar{p} \bar{q}} \right), \qquad
#'       V_A = d \left( \frac{f_1}{n_1 p_1 q_1}
#'                    + \frac{f_2}{n_2 p_2 q_2} \right),}{V_0 = d ( f_1/(n_1 pbar qbar) + f_2/(n_2 pbar qbar) ), V_A = d ( f_1/(n_1 p_1 q_1) + f_2/(n_2 p_2 q_2) ),}
#'
#' where \eqn{d} is `deff` and \eqn{f_i} the finite population correction of
#' group \eqn{i}. Under the null the two groups share one proportion, and the
#' estimator of it is the pooled one, so \eqn{\bar{p}}{pbar} is the
#' **sample-size-weighted** mean
#'
#' \deqn{\bar{p} = \frac{n_1 p_1 + n_2 p_2}{n_1 + n_2}
#'               = \frac{r p_1 + p_2}{r + 1},}{pbar = (n_1 p_1 + n_2 p_2)/(n_1 + n_2) = (r p_1 + p_2)/(r + 1),}
#'
#' with \eqn{r} the allocation `ratio`. At `ratio = 1` this is
#' \eqn{(p_1 + p_2) / 2}. Away from it the weighted and unweighted nulls give
#' materially different sizes, so the weighting is not a refinement. Sizing
#' `p1 = 0.1` against `p2 = 0.2` at `ratio = 4` needs 523 and 131, against
#' 457 and 115 for an unweighted null: 14 percent more fieldwork, because the
#' larger group is the one with the smaller proportion and pulls the pooled
#' null toward it. Both the size-solving and the power-computing path use this
#' \eqn{\bar{p}}{pbar}, so `power_prop()` inverts itself under any allocation.
#'
#' ## Normal approximation
#'
#' All three methods compute critical values and power from the standard
#' normal, with no degrees-of-freedom correction and the variance treated
#' as known. This is the convention for survey-scale samples and is what
#' lets `deff`, a finite `N`, and `resp_rate` enter the variance directly.
#' The choice of `method` addresses accuracy in \eqn{p}, not in \eqn{n}:
#' `"arcsine"` and `"logodds"` improve the normal approximation for an
#' extreme proportion, but none of the three is exact at small `n`.
#'
#' The `df` argument that [n_prop()], [n_mean()] and [n_alloc()] accept has
#' no counterpart here, and its absence is a decision rather than an
#' omission. There the quantile is the half-width of a confidence interval
#' and a t quantile substitutes for a normal one directly. Here it is a
#' normal deviate for an alternative, and a t-based power calculation is a
#' different procedure. Passing `df` is an error that says so.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapter 4.
#'
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' @seealso [power_mean()] for continuous outcomes, [power_did()] for
#'   difference-in-differences, [n_prop()] for estimation precision.
#'   With `overlap` set, this is the two-occasion change of one population:
#'   [n_change()] and [prec_change()] size and evaluate the same quantity as
#'   an estimate with a margin of error rather than as a test, and
#'   [design_overlap()] derives the overlap from a rotation schedule.
#'
#' @examples
#' # Sample size to detect a 5pp change from 30%
#' power_prop(p1 = 0.30, p2 = 0.35)
#'
#' # Power given n = 500
#' power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
#'
#' # MDE with n = 1000
#' power_prop(p1 = 0.30, n = 1000)
#'
#' # Arcsine transform for rare proportions
#' power_prop(p1 = 0.15, p2 = 0.18, alternative = "one.sided",
#'            method = "arcsine")
#'
#' # Log-odds transform
#' power_prop(p1 = 0.15, p2 = 0.18, alternative = "one.sided",
#'            method = "logodds")
#'
#' # Allocation ratio 2:1
#' power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
#'
#' @export
power_prop <- function(p1, ...) {
  if (!missing(p1)) {
    .res <- .dispatch_plan(p1, "p1", power_prop.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("power_prop")
}

#' @rdname power_prop
#' @export
power_prop.default <- function(p1, ..., p2 = NULL, n = NULL, power = 0.80,
                       alpha = 0.05, N = Inf, deff = 1,
                       resp_rate = 1,
                       alternative = c("two.sided", "one.sided"),
                       ratio = 1,
                       overlap = 0, overlap_cor = 0,
                       method = c("wald", "arcsine", "logodds"),
                       plan = NULL) {
  .plan <- .merge_plan_args(plan, power_prop.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(power_prop.default, c(.plan, list(...))))
  .stop_power_df(...)
  .check_unused_dots(...)
  check_proportion(p1, "p1")
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  check_alpha(alpha)
  N_pair <- .check_power_N(N)
  check_deff(deff)
  check_overlap(overlap)
  check_overlap_cor(overlap_cor)
  .check_overlap_N(overlap, N_pair)
  check_resp_rate(resp_rate)
  ratio <- .resolve_ratio(n, ratio)

  if (method != "wald" && overlap > 0)
    stop("overlap is only supported with method = 'wald'", call. = FALSE)

  null_count <- is.null(p2) + is.null(n) + is.null(power)
  if (null_count != 1L)
    stop("leave exactly one of 'n', 'power', or 'p2' as NULL", call. = FALSE)

  if (!is.null(p2)) check_proportion(p2, "p2")
  if (!is.null(n)) {
    n <- .check_power_n(n)
    .check_overlap_n(overlap, n = n)
  } else {
    .check_overlap_n(overlap, ratio = ratio)
  }
  if (!is.null(power)) check_proportion(power, "power")
  if (!is.null(power) && (is.null(n) || is.null(p2))) {
    .check_power_target(power, alpha)
  }
  if (method == "wald" && !is.null(p2)) {
    .check_bernoulli_cor(c(p1, p2), overlap_cor, overlap)
  }

  params <- list(p1 = p1, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate, alternative = alternative,
                 ratio = ratio, overlap = overlap, overlap_cor = overlap_cor,
                 method = method)

  if (is.null(n)) {
    if (p1 == p2) stop("'p1' and 'p2' must differ", call. = FALSE)
    params$p2 <- p2
    params$power <- power

    res <- switch(method,
      wald    = .power_prop_n_wald(p1, p2, power, alpha, N_pair, deff,
                                  alternative, overlap, overlap_cor, ratio, resp_rate),
      arcsine = .power_prop_n_arcsine(p1, p2, power, alpha, N_pair, deff,
                                      alternative, ratio, resp_rate),
      logodds = .power_prop_n_logodds(p1, p2, power, alpha, N_pair, deff,
                                      alternative, ratio, resp_rate))

    n_eff <- res * resp_rate
    achieved_power <- switch(method,
      wald = .power_prop_power_wald(
        p1, p2, n_eff, alpha, N_pair, deff,
        alternative, overlap, overlap_cor
      ),
      arcsine = .power_prop_power_arcsine(
        p1, p2, n_eff, alpha, N_pair, deff, alternative
      ),
      logodds = .power_prop_power_logodds(
        p1, p2, n_eff, alpha, N_pair, deff, alternative
      )
    )

    .new_svyplan_power(n = res, power = achieved_power, effect = abs(p1 - p2),
                       type = "proportion", solved = "n", params = params)

  } else if (is.null(power)) {
    if (p1 == p2) stop("'p1' and 'p2' must differ", call. = FALSE)
    params$p2 <- p2
    params$n <- n
    .check_gross_n(if (length(n) == 1L) c(n, n) else n, N_pair,
                   label = c("group 1", "group 2"))
    n_eff <- n * resp_rate

    res <- switch(method,
      wald    = .power_prop_power_wald(p1, p2, n_eff, alpha, N_pair, deff,
                                      alternative, overlap, overlap_cor),
      arcsine = .power_prop_power_arcsine(p1, p2, n_eff, alpha, N_pair, deff,
                                          alternative),
      logodds = .power_prop_power_logodds(p1, p2, n_eff, alpha, N_pair, deff,
                                          alternative))

    .new_svyplan_power(n = n, power = res, effect = abs(p1 - p2),
                       type = "proportion", solved = "power", params = params)

  } else {
    params$n <- n
    params$power <- power
    .check_gross_n(if (length(n) == 1L) c(n, n) else n, N_pair,
                   label = c("group 1", "group 2"))
    n_eff <- n * resp_rate

    res <- switch(method,
      wald    = .power_prop_mde_wald(p1, n_eff, power, alpha, N_pair, deff,
                                    alternative, overlap, overlap_cor),
      arcsine = .power_prop_mde_arcsine(p1, n_eff, power, alpha, N_pair, deff,
                                        alternative),
      logodds = .power_prop_mde_logodds(p1, n_eff, power, alpha, N_pair, deff,
                                        alternative))

    .new_svyplan_power(n = n, power = power, effect = abs(p1 - res),
                       type = "proportion", solved = "mde",
                       params = c(params, list(p2 = res)))
  }
}

#' @keywords internal
#' @noRd
.prop_var <- function(p1, p2, overlap, overlap_cor) {
  p1 * (1 - p1) + p2 * (1 - p2) -
    2 * overlap * overlap_cor * sqrt(p1 * (1 - p1) * p2 * (1 - p2))
}

## Wald internals

.power_prop_n_wald <- function(p1, p2, power, alpha, N_pair, deff,
                               alternative, overlap, overlap_cor, ratio, resp_rate) {
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  icc <- abs(p1 - p2)
  q1 <- 1 - p1; q2 <- 1 - p2
  ov_term <- 2 * overlap * overlap_cor * sqrt(p1 * q1 * p2 * q2)

  r <- ratio
  power_n2 <- function(n2) {
    n_eff <- if (r == 1) c(n2, n2) else c(r * n2, n2)
    n_eff <- n_eff * resp_rate
    .power_prop_power_wald(p1, p2, n_eff, alpha, N_pair, deff,
                           alternative, overlap, overlap_cor)
  }
  if (alternative == "one.sided" && all(is.infinite(N_pair))) {
    V_r <- p1 * q1 / r + p2 * q2 - ov_term
    n2 <- (z_a + z_b)^2 * V_r * deff / icc^2 / resp_rate
    n2 <- max(n2, 2, 2 / r)
  } else {
    n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
  }
  if (r == 1) n2 else c(r * n2, n2)
}

.power_prop_power_wald <- function(p1, p2, n_eff, alpha, N_pair, deff,
                                   alternative, overlap, overlap_cor) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff

  z_a <- .z_alpha(alpha, alternative)
  icc <- abs(p1 - p2)

  V_d <- .diff_var_fpc(n_vec, .bernoulli_var(c(p1, p2), N_pair), N_pair, deff,
                       overlap, overlap_cor)
  V_d <- .safe_variance(V_d, "difference variance")
  if (V_d == 0) return(1)

  se <- sqrt(V_d)
  pw <- pnorm(icc / se - z_a)
  if (alternative == "two.sided")
    pw <- pw + pnorm(-icc / se - z_a)
  min(pw, 1)
}

.power_prop_mde_wald <- function(p1, n_eff, power, alpha, N_pair, deff,
                                  alternative, overlap, overlap_cor) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  if (all(vapply(seq_len(2L), function(i) {
    .fpc_factor_prop(n_vec[i], N_pair[i]) == 0
  }, logical(1L)))) .stop_census_mde()

  target_fn <- function(p2) {
    .power_prop_power_wald(p1, p2, n_eff, alpha, N_pair, deff,
                           alternative, overlap, overlap_cor) - power
  }

  eps <- 1e-8
  roots <- list()
  feasible <- .bernoulli_cor_interval(p1, overlap_cor, overlap)

  up_lo <- p1 + eps; up_hi <- min(1 - eps, feasible[2L])
  if (up_lo < up_hi) {
    f_lo <- target_fn(up_lo); f_hi <- target_fn(up_hi)
    if (is.finite(f_lo) && is.finite(f_hi)) {
      if (f_lo == 0) roots <- c(roots, list(up_lo))
      else if (f_hi == 0) roots <- c(roots, list(up_hi))
      else if (sign(f_lo) != sign(f_hi))
        roots <- c(roots, list(uniroot(target_fn, c(up_lo, up_hi), tol = eps)$root))
    }
  }

  dn_lo <- max(eps, feasible[1L]); dn_hi <- p1 - eps
  if (dn_lo < dn_hi) {
    f_lo <- target_fn(dn_lo); f_hi <- target_fn(dn_hi)
    if (is.finite(f_lo) && is.finite(f_hi)) {
      if (f_hi == 0) roots <- c(roots, list(dn_hi))
      else if (f_lo == 0) roots <- c(roots, list(dn_lo))
      else if (sign(f_lo) != sign(f_hi))
        roots <- c(roots, list(uniroot(target_fn, c(dn_lo, dn_hi), tol = eps)$root))
    }
  }

  if (length(roots) == 0L)
    stop("no detectable alternative exists for the given n and power", call. = FALSE)

  dists <- vapply(roots, function(r) abs(r - p1), numeric(1L))
  roots[[which.min(dists)]]
}

## Arcsine internals

.power_prop_n_arcsine <- function(p1, p2, power, alpha, N_pair, deff,
                                   alternative, ratio, resp_rate) {
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  diff_phi <- asin(sqrt(p1)) - asin(sqrt(p2))

  r <- ratio
  power_n2 <- function(n2) {
    n_eff <- if (r == 1) c(n2, n2) else c(r * n2, n2)
    n_eff <- n_eff * resp_rate
    .power_prop_power_arcsine(p1, p2, n_eff, alpha, N_pair, deff,
                              alternative)
  }
  if (alternative == "one.sided" && all(is.infinite(N_pair))) {
    n2 <- ((z_a + z_b) / abs(diff_phi))^2 * deff * (1 / r + 1) / 4
    n2 <- n2 / resp_rate
    n2 <- max(n2, 2, 2 / r)
  } else {
    n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
  }
  if (r == 1) n2 else c(r * n2, n2)
}

.power_prop_power_arcsine <- function(p1, p2, n_eff, alpha, N_pair, deff,
                                      alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff

  z_a <- .z_alpha(alpha, alternative)
  diff_phi <- asin(sqrt(p1)) - asin(sqrt(p2))

  fpc1 <- .fpc_factor_prop(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor_prop(n_vec[2], N_pair[2])
  se_phi <- sqrt(deff * (fpc1 / (4 * n_vec[1]) + fpc2 / (4 * n_vec[2])))

  if (se_phi == 0) return(1)

  z_test <- abs(diff_phi) / se_phi - z_a
  pw <- pnorm(z_test)
  if (alternative == "two.sided")
    pw <- pw + pnorm(-abs(diff_phi) / se_phi - z_a)
  min(pw, 1)
}

.power_prop_mde_arcsine <- function(p1, n_eff, power, alpha, N_pair, deff,
                                     alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  fpc1 <- .fpc_factor_prop(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor_prop(n_vec[2], N_pair[2])
  se_phi <- sqrt(deff * (fpc1 / (4 * n_vec[1]) + fpc2 / (4 * n_vec[2])))
  if (se_phi == 0) .stop_census_mde()
  mde_phi <- .solve_normal_mde(power, se_phi, alpha, alternative)

  phi1 <- asin(sqrt(p1))

  roots <- list()
  if (phi1 - mde_phi >= 0) {
    p2_dn <- sin(phi1 - mde_phi)^2
    if (p2_dn > 0 && p2_dn < 1 && p2_dn != p1) roots <- c(roots, list(p2_dn))
  }
  if (phi1 + mde_phi <= pi / 2) {
    p2_up <- sin(phi1 + mde_phi)^2
    if (p2_up > 0 && p2_up < 1 && p2_up != p1) roots <- c(roots, list(p2_up))
  }

  if (length(roots) == 0L)
    stop("no detectable alternative exists for the given n and power", call. = FALSE)

  dists <- vapply(roots, function(r) abs(r - p1), numeric(1L))
  roots[[which.min(dists)]]
}

## Log-odds internals

.power_prop_n_logodds <- function(p1, p2, power, alpha, N_pair, deff,
                                   alternative, ratio, resp_rate) {
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  q1 <- 1 - p1; q2 <- 1 - p2
  diff_phi <- log(p1 / q1) - log(p2 / q2)
  # Allocation-weighted pooled null. resp_rate is common and cancels.
  p_bar <- (ratio * p1 + p2) / (ratio + 1)
  q_bar <- 1 - p_bar

  r <- ratio
  power_n2 <- function(n2) {
    n_eff <- if (r == 1) c(n2, n2) else c(r * n2, n2)
    n_eff <- n_eff * resp_rate
    .power_prop_power_logodds(p1, p2, n_eff, alpha, N_pair, deff,
                              alternative)
  }
  if (alternative == "one.sided" && all(is.infinite(N_pair))) {
    V0_coeff <- (1 / r + 1) / (p_bar * q_bar)
    VA_coeff <- 1 / (r * p1 * q1) + 1 / (p2 * q2)
    n2 <- ((z_a * sqrt(V0_coeff) + z_b * sqrt(VA_coeff)) /
      abs(diff_phi))^2 * deff / resp_rate
    n2 <- max(n2, 2, 2 / r)
  } else {
    n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
  }
  if (r == 1) n2 else c(r * n2, n2)
}

.power_prop_power_logodds <- function(p1, p2, n_eff, alpha, N_pair, deff,
                                      alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff

  z_a <- .z_alpha(alpha, alternative)
  q1 <- 1 - p1; q2 <- 1 - p2
  diff_phi <- log(p1 / q1) - log(p2 / q2)
  # Pooled null proportion, weighted by the sizes actually realized here. The
  # size-solving path weights by 'ratio', which is these two sizes in the same
  # proportion, so the two paths agree.
  p_bar <- (n_vec[1] * p1 + n_vec[2] * p2) / (n_vec[1] + n_vec[2])
  q_bar <- 1 - p_bar

  fpc1 <- .fpc_factor_prop(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor_prop(n_vec[2], N_pair[2])

  V0 <- deff * (fpc1 / (n_vec[1] * p_bar * q_bar) +
                fpc2 / (n_vec[2] * p_bar * q_bar))
  VA <- deff * (fpc1 / (n_vec[1] * p1 * q1) +
                fpc2 / (n_vec[2] * p2 * q2))

  if (VA == 0) return(1)

  z_test <- (abs(diff_phi) - z_a * sqrt(V0)) / sqrt(VA)
  pw <- pnorm(z_test)
  if (alternative == "two.sided")
    pw <- pw + pnorm((-abs(diff_phi) - z_a * sqrt(V0)) / sqrt(VA))
  min(pw, 1)
}

.power_prop_mde_logodds <- function(p1, n_eff, power, alpha, N_pair, deff,
                                     alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  if (all(vapply(seq_len(2L), function(i) {
    .fpc_factor_prop(n_vec[i], N_pair[i]) == 0
  }, logical(1L)))) .stop_census_mde()

  target_fn <- function(p2) {
    .power_prop_power_logodds(p1, p2, n_eff, alpha, N_pair, deff,
                              alternative) - power
  }

  eps <- 1e-8
  roots <- list()

  up_lo <- p1 + eps; up_hi <- 1 - eps
  if (up_lo < up_hi) {
    f_lo <- target_fn(up_lo); f_hi <- target_fn(up_hi)
    if (is.finite(f_lo) && is.finite(f_hi) && sign(f_lo) != sign(f_hi))
      roots <- c(roots, list(uniroot(target_fn, c(up_lo, up_hi), tol = eps)$root))
  }

  dn_lo <- eps; dn_hi <- p1 - eps
  if (dn_lo < dn_hi) {
    f_lo <- target_fn(dn_lo); f_hi <- target_fn(dn_hi)
    if (is.finite(f_lo) && is.finite(f_hi) && sign(f_lo) != sign(f_hi))
      roots <- c(roots, list(uniroot(target_fn, c(dn_lo, dn_hi), tol = eps)$root))
  }

  if (length(roots) == 0L)
    stop("no detectable alternative exists for the given n and power", call. = FALSE)

  dists <- vapply(roots, function(r) abs(r - p1), numeric(1L))
  roots[[which.min(dists)]]
}
