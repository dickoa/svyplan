#' Power Analysis for Means
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-sample test of means. Leave exactly one of `n`, `power`, or `effect`
#' as `NULL` to solve for that quantity.
#'
#' @param effect Absolute difference in means (effect-size magnitude, positive).
#'   Leave `NULL` to solve for MDE.
#' @param var Within-group variance. Scalar (equal variances in both groups)
#'   or length-2 vector `c(var1, var2)` for unequal group variances.
#' @param n Per-group sample size. Scalar (equal groups) or length-2 vector
#'   `c(n1, n2)` for unequal groups. Leave `NULL` to solve for sample size.
#' @param power Target power, in (0, 1). Leave `NULL` to solve for power.
#' @param alpha Significance level, default 0.05.
#' @param N Population size for finite-population correction. A scalar applies
#'   to both groups; a length-2 vector `c(N1, N2)` sets group-specific
#'   population sizes. `Inf` disables FPC for the corresponding group.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The sample size is inflated by `1 / resp_rate`.
#' @param alternative Character: `"two.sided"` (default) or `"one.sided"`.
#' @param ratio Allocation ratio n1/n2 (default 1). Only used when solving
#'   for n (`n = NULL`). For example, `ratio = 2` means group 1 gets twice
#'   the sample of group 2.
#' @param overlap Panel overlap fraction in \[0, 1\], for repeated surveys.
#'   Defined as the fraction of group 1 that also appears in group 2
#'   (`overlap = n12 / n1`).
#' @param rho Correlation between occasions in \[0, 1\].
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param ... Additional arguments passed to methods.
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
#' To specify the effect in terms of Cohen's d (standardized effect size),
#' convert via `effect = d * sqrt(mean(var))`, where `d` follows Cohen's
#' conventions: 0.2 (small), 0.5 (medium), 0.8 (large).
#'
#' When `var` is a length-2 vector, the variance of the difference is:
#'
#' \deqn{V = \sigma^2_1 / r + \sigma^2_2 - 2 \cdot \text{overlap} \cdot
#'   \rho \cdot \sigma_1 \sigma_2}
#'
#' where `r` is the allocation ratio n1/n2 (default 1). When `var` is
#' scalar and `ratio = 1`, this simplifies to the familiar
#' `V = 2 * var * (1 - overlap * rho)`.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapter 4.
#'
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' @seealso [power_prop()] for proportions, [n_mean()] for estimation
#'   precision.
#'
#' @examples
#' # Sample size to detect a difference of 5 with variance 100
#' power_mean(effect = 5, var = 100)
#'
#' # Power given n = 200
#' power_mean(effect = 5, var = 100, n = 200, power = NULL)
#'
#' # MDE with n = 500
#' power_mean(var = 100, n = 500)
#'
#' # With design effect
#' power_mean(effect = 5, var = 100, deff = 1.5)
#'
#' # Unequal group variances
#' power_mean(effect = 5, var = c(80, 120))
#'
#' # Allocation ratio 2:1
#' power_mean(effect = 5, var = 100, ratio = 2)
#'
#' @export
power_mean <- function(effect, ...) {
  if (!missing(effect)) {
    .res <- .dispatch_plan(effect, "effect", power_mean.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("power_mean")
}

#' @rdname power_mean
#' @export
power_mean.default <- function(effect = NULL, var, n = NULL, power = 0.80,
                       alpha = 0.05, N = Inf, deff = 1,
                       resp_rate = 1,
                       alternative = c("two.sided", "one.sided"),
                       ratio = 1,
                       overlap = 0, rho = 0,
                       plan = NULL, ...) {
  .plan <- .merge_plan_args(plan, power_mean.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(power_mean.default, c(.plan, list(...))))
  alternative <- match.arg(alternative)
  var_pair <- .as_pair(var, "var")
  check_alpha(alpha)
  N_pair <- .check_power_N(N)
  check_deff(deff)
  check_overlap(overlap)
  check_rho(rho)
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
                 ratio = ratio, overlap = overlap, rho = rho)


  ov_term <- 2 * overlap * rho * sqrt(var_pair[1] * var_pair[2])

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

        if (any(!is.infinite(N_pair) & n_eff >= N_pair)) return(1)

        fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
        fpc2 <- .fpc_factor(n_eff[2], N_pair[2])
        V_d <- deff * (var_pair[1] * fpc1 / n_eff[1] +
                       var_pair[2] * fpc2 / n_eff[2])
        if (overlap > 0) {
          fpc_ov <- .fpc_factor(n_eff[2], N_pair[2])
          V_d <- V_d - 2 * overlap * rho * sqrt(var_pair[1] * var_pair[2]) *
                 deff * fpc_ov / n_eff[2]
        }
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
    n_eff <- n_vec * resp_rate

    if (any(!is.infinite(N_pair) & n_eff >= N_pair)) {
      warning("effective sample size >= population size; returning power = 1 (census)",
              call. = FALSE)
      pw <- 1
    } else {
      fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
      fpc2 <- .fpc_factor(n_eff[2], N_pair[2])
      V_d <- deff * (var_pair[1] * fpc1 / n_eff[1] +
                     var_pair[2] * fpc2 / n_eff[2])
      if (overlap > 0) {
        fpc_ov <- .fpc_factor(n_eff[2], N_pair[2])
        V_d <- V_d - 2 * overlap * rho * sqrt(var_pair[1] * var_pair[2]) *
               deff * fpc_ov / n_eff[2]
      }
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
    }

    .new_svyplan_power(n = n, power = pw, effect = effect,
                       type = "mean", solved = "power", params = params)

  } else {
    params$n <- n
    params$power <- power
    z_b <- qnorm(power)
    n_vec <- if (length(n) == 1L) c(n, n) else n
    n_eff <- n_vec * resp_rate

    fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
    fpc2 <- .fpc_factor(n_eff[2], N_pair[2])
    V_d <- deff * (var_pair[1] * fpc1 / n_eff[1] +
                   var_pair[2] * fpc2 / n_eff[2])
    if (overlap > 0) {
      fpc_ov <- .fpc_factor(n_eff[2], N_pair[2])
      V_d <- V_d - 2 * overlap * rho * sqrt(var_pair[1] * var_pair[2]) *
             deff * fpc_ov / n_eff[2]
    }
    V_d <- .safe_variance(V_d, "difference variance")
    se <- sqrt(V_d)
    mde <- (z_a + z_b) * se

    .new_svyplan_power(n = n, power = power, effect = mde,
                       type = "mean", solved = "mde", params = params)
  }
}
