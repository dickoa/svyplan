#' Power analysis for means
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-sample test of means. Leave exactly one of `n`, `power`, or `effect`
#' as `NULL` to solve for that quantity.
#'
#' @param var Within-group variance. Scalar (equal variances in both groups)
#'   or length-2 vector `c(var1, var2)` for unequal group variances.
#' @param sd Population standard deviation, an alternative spelling of
#'   `var`. Supply exactly one of `var` or `sd`. Stratum frames and
#'   published survey reports usually quote standard deviations.
#' @param effect Absolute difference in means (effect-size magnitude, positive).
#'   Leave `NULL` to solve for MDE.
#' @param n Per-group sample size. Scalar (equal groups) or length-2 vector
#'   `c(n1, n2)` for unequal groups. Leave `NULL` to solve for sample size.
#' @param power Target power, in (0, 1). Leave `NULL` to solve for power.
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
#'   interest, which requires a smaller sample for the same power.
#' @param ratio Allocation ratio n1/n2 (default 1). Only used when solving
#'   for n (`n = NULL`). For example, `ratio = 2` means group 1 gets twice
#'   the sample of group 2.
#' @param overlap Panel overlap fraction in \[0, 1\], for repeated surveys.
#'   Defined as the fraction of group 1 that also appears in group 2
#'   (`overlap = n12 / n1`). A positive value describes a coordinated
#'   design with that many units deliberately held in common, which with a
#'   finite `N` requires both occasions to sample one population. The
#'   default 0 is the ordinary two-group comparison, where the groups are
#'   independent and may be different populations of different sizes. It is
#'   not the same model as a deliberately disjoint pair, which is why the
#'   two need not agree in the limit when `N` is small.
#' @param overlap_cor Correlation between occasions in \[0, 1\].
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#'
#' @return A `svyplan_power` object with components:
#' \describe{
#'   \item{`n`}{Per-group sample size (scalar or length-2 for unequal groups).}
#'   \item{`power`}{Achieved power.}
#'   \item{`effect`}{Effect size (difference in means).}
#'   \item{`solved`}{Which quantity was solved for
#'     (`"n"`, `"power"`, or `"mde"`).}
#'   \item{`params`}{List of input parameters.}
#' }
#'
#' @details
#' The `var` argument is the within-group population variance. Estimate it
#' from a pilot study, a previous survey, or published data for a similar
#' population. When uncertain, use a conservative (larger) estimate. This
#' inflates the sample size and reduces the risk of an underpowered design.
#'
#' To specify the effect in terms of Cohen's d (standardized effect size),
#' convert via `effect = d * sqrt(mean(var))`, where `d` follows Cohen's
#' conventions: 0.2 (small), 0.5 (medium), 0.8 (large).
#'
#' When `var` is a length-2 vector, the variance of the difference is:
#'
#' \deqn{V = \sigma^2_1 / r + \sigma^2_2 - 2 \cdot \text{overlap} \cdot
#'   \rho \cdot \sigma_1 \sigma_2}{V = sigma^2_1 / r + sigma^2_2 - 2 * overlap * rho * sigma_1 sigma_2}
#'
#' where `r` is the allocation ratio n1/n2 (default 1). When `var` is
#' scalar and `ratio = 1`, this simplifies to the familiar
#' `V = 2 * var * (1 - overlap * overlap_cor)`.
#'
#' With a finite `N` the correction applies to the marginal terms but not
#' to the overlap covariance, which carries a single \eqn{1/N}: for two
#' SRSWOR samples sharing \eqn{k = overlap \cdot n_1}{k = overlap * n_1} units,
#' \eqn{Cov(\bar y_1, \bar y_2) = \rho S_1S_2\{k/(n_1n_2) - 1/N\}}{Cov(ybar_1, ybar_2) = rho S_1S_2\{k/(n_1n_2) - 1/N\}}. At
#' `overlap_cor = 1` with equal sizes and variances the population terms
#' cancel exactly and the difference variance is
#' \eqn{2S^2(1 - overlap)/n}, free of `N`.
#'
#' ## Normal approximation
#'
#' Critical values and power are computed from the standard normal, not
#' from a t distribution with an estimated denominator: no degrees of
#' freedom enter, and the variance is treated as known. This is the
#' convention for survey-scale samples, where the two agree closely, and
#' it is what makes `deff` and a finite `N` insertable directly into the
#' variance. At small `n` the sizes are correspondingly smaller than
#' [stats::power.t.test()], which uses a noncentral t. Use that function
#' instead when the sample is small enough for the difference to matter.
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
#' @seealso [power_prop()] for proportions, [power_did()] for
#'   difference-in-differences, [n_mean()] for estimation precision.
#'   With `overlap` set, this is the two-occasion change of one population:
#'   [n_change()] and [prec_change()] size and evaluate the same quantity as
#'   an estimate with a margin of error rather than as a test, and
#'   [design_overlap()] derives the overlap from a rotation schedule.
#'
#' @examples
#' # Sample size to detect a difference of 5 with variance 100
#' power_mean(100, effect = 5)
#'
#' # Power given n = 200
#' power_mean(100, effect = 5, n = 200, power = NULL)
#'
#' # MDE with n = 500
#' power_mean(100, n = 500)
#'
#' # With design effect
#' power_mean(100, effect = 5, deff = 1.5)
#'
#' # Unequal group variances
#' power_mean(c(80, 120), effect = 5)
#'
#' # Allocation ratio 2:1
#' power_mean(100, effect = 5, ratio = 2)
#'
#' @export
power_mean <- function(var = NULL, ...) {
  if (!is.null(var)) {
    .res <- .dispatch_plan(var, "var", power_mean.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("power_mean")
}

#' @rdname power_mean
#' @export
power_mean.default <- function(var = NULL, ..., sd = NULL, effect = NULL, n = NULL, power = 0.80,
                       alpha = 0.05, N = Inf, deff = 1,
                       resp_rate = 1,
                       alternative = c("two.sided", "one.sided"),
                       ratio = 1,
                       overlap = 0, overlap_cor = 0,
                       plan = NULL) {
  .plan <- .merge_plan_args(plan, power_mean.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(power_mean.default, c(.plan, list(...))))
  .stop_power_df(...)
  .check_unused_dots(...)
  alternative <- match.arg(alternative)
  var <- .resolve_var(var, sd)
  var_pair <- .as_pair(var, "var")
  check_alpha(alpha)
  N_pair <- .check_power_N(N)
  check_deff(deff)
  check_overlap(overlap)
  check_overlap_cor(overlap_cor)
  .check_overlap_N(overlap, N_pair)
  check_resp_rate(resp_rate)
  ratio <- .resolve_ratio(n, ratio)

  null_count <- is.null(effect) + is.null(n) + is.null(power)
  if (null_count != 1L)
    stop("leave exactly one of 'n', 'power', or 'effect' as NULL", call. = FALSE)

  if (!is.null(effect)) check_scalar(effect, "effect")
  if (!is.null(n)) {
    n <- .check_power_n(n)
    .check_overlap_n(overlap, n = n)
  } else {
    .check_overlap_n(overlap, ratio = ratio)
  }
  if (!is.null(power)) check_proportion(power, "power")

  z_a <- .z_alpha(alpha, alternative)

  params <- list(var = var, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate, alternative = alternative,
                 ratio = ratio, overlap = overlap, overlap_cor = overlap_cor)


  ov_term <- 2 * overlap * overlap_cor * sqrt(var_pair[1] * var_pair[2])

  if (is.null(n)) {
    params$effect <- effect
    params$power <- power
    z_b <- qnorm(power)

    if (all(is.infinite(N_pair))) {
      if (ratio == 1) {
        V <- var_pair[1] + var_pair[2] - ov_term
        n2_0 <- (z_a + z_b)^2 * V * deff / effect^2
        n0 <- n2_0 / resp_rate
      } else {
        r <- ratio
        V_r <- var_pair[1] / r + var_pair[2] - ov_term
        n2 <- (z_a + z_b)^2 * V_r * deff / effect^2
        n2 <- n2 / resp_rate
        n0 <- c(r * n2, n2)
      }
    } else {
      r <- ratio
      power_n2 <- function(n2) {
        n_vec <- if (r == 1) c(n2, n2) else c(r * n2, n2)
        n_eff <- n_vec * resp_rate

        V_d <- .diff_var_fpc(n_eff, var_pair, N_pair, deff, overlap, overlap_cor)
        V_d <- .safe_variance(V_d, "difference variance")
        if (V_d == 0) return(1)
        se <- sqrt(V_d)
        pw <- pnorm(abs(effect) / se - z_a)
        if (alternative == "two.sided")
          pw <- pw + pnorm(-abs(effect) / se - z_a)
        min(pw, 1)
      }

      n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
      n0 <- if (r == 1) n2 else c(r * n2, n2)
    }

    .new_svyplan_power(n = n0, power = power, effect = effect,
                       type = "mean", solved = "n", params = params)

  } else if (is.null(power)) {
    params$effect <- effect
    params$n <- n
    n_vec <- if (length(n) == 1L) c(n, n) else n
    .check_gross_n(n_vec, N_pair, label = c("group 1", "group 2"))
    n_eff <- n_vec * resp_rate

    V_d <- .diff_var_fpc(n_eff, var_pair, N_pair, deff, overlap, overlap_cor)
    V_d <- .safe_variance(V_d, "difference variance")
    if (V_d == 0) {
      pw <- 1
    } else {
      se <- sqrt(V_d)
      pw <- pnorm(abs(effect) / se - z_a)
      if (alternative == "two.sided")
        pw <- pw + pnorm(-abs(effect) / se - z_a)
      pw <- min(pw, 1)
    }

    .new_svyplan_power(n = n, power = pw, effect = effect,
                       type = "mean", solved = "power", params = params)

  } else {
    params$n <- n
    params$power <- power
    z_b <- qnorm(power)
    n_vec <- if (length(n) == 1L) c(n, n) else n
    .check_gross_n(n_vec, N_pair, label = c("group 1", "group 2"))
    n_eff <- n_vec * resp_rate

    V_d <- .diff_var_fpc(n_eff, var_pair, N_pair, deff, overlap, overlap_cor)
    V_d <- .safe_variance(V_d, "difference variance")
    se <- sqrt(V_d)
    mde <- (z_a + z_b) * se

    .new_svyplan_power(n = n, power = power, effect = mde,
                       type = "mean", solved = "mde", params = params)
  }
}
