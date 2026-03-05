#' Power Analysis for Difference-in-Differences Designs
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-group, two-period difference-in-differences (DiD) contrast.
#' Leave exactly one of `n`, `power`, or `effect` as `NULL`.
#'
#' @param treat Numeric length-2 vector for treated group outcomes:
#'   `c(baseline, endline)`.
#' @param control Numeric length-2 vector for control group outcomes:
#'   `c(baseline, endline)`.
#' @param outcome Outcome scale: `"mean"` (default) or `"prop"`.
#' @param var Outcome variance for `outcome = "mean"`.
#'   Length 1: common variance for all four cells.
#'   Length 2: group-specific variances `c(var_treat, var_control)`,
#'   assumed equal across waves.
#'   Length 4: cell-specific variances in order
#'   `c(var_treat_baseline, var_treat_endline, var_control_baseline,
#'   var_control_endline)`.
#' @param effect Absolute DiD effect size to detect (> 0).
#'   Leave `NULL` to solve for MDE.
#' @param n Per-arm sample size per wave. Scalar (equal treated/control)
#'   or length-2 vector `c(n_treat, n_control)`. Leave `NULL` to solve
#'   for n.
#' @param power Target power, in (0, 1). Leave `NULL` to solve for power.
#' @param alpha Significance level, default 0.05.
#' @param N Population size for finite-population correction.
#'   Scalar applies to both arms; length-2 vector `c(N_treat, N_control)`.
#'   `Inf` (default) disables FPC.
#' @param deff Design effect multiplier (> 0).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1.
#' @param alternative Character: `"two.sided"` (default) or `"one.sided"`.
#' @param ratio Allocation ratio `n_treat / n_control` (default 1).
#'   Used only when solving for `n` (`n = NULL`).
#' @param overlap Panel overlap fraction in \[0, 1\] within each arm
#'   across baseline and endline.
#' @param rho Correlation between baseline and endline outcomes within
#'   overlapping units, in \[0, 1\].
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `svyplan_power` object with components:
#' \describe{
#'   \item{`n`}{Per-arm sample size (scalar or length-2).}
#'   \item{`power`}{Achieved power.}
#'   \item{`effect`}{DiD effect size.}
#'   \item{`type`}{`"did_prop"` or `"did_mean"`.}
#'   \item{`solved`}{Which quantity was solved for
#'     (`"n"`, `"power"`, or `"mde"`).}
#'   \item{`params`}{List of input parameters.}
#' }
#'
#' @details
#' The DiD estimand is
#' `(treat_endline - treat_baseline) - (control_endline - control_baseline)`.
#'
#' Per-arm variance of the before-after change accounts for panel overlap:
#'
#' \deqn{V_{\text{arm}} = V_{\text{pre}} + V_{\text{post}}
#'   - 2 \cdot \text{overlap} \cdot \rho \cdot
#'   \sqrt{V_{\text{pre}} \cdot V_{\text{post}}}}
#'
#' The DiD test statistic variance is then:
#'
#' \deqn{V_d = \text{deff} \left(
#'   \frac{V_{\text{trt}} \cdot \text{fpc}_t}{n_t} +
#'   \frac{V_{\text{ctrl}} \cdot \text{fpc}_c}{n_c}
#' \right)}
#'
#' When `overlap = 0`, this reduces to the classical flat-variance formula.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapter 4.
#'
#' @seealso [power_prop()] for two-sample proportions, [power_mean()] for
#'   two-sample means.
#'
#' @examples
#' # DiD sample size for means
#' power_did(
#'   treat = c(50, 55), control = c(50, 52),
#'   outcome = "mean", var = 100, effect = 3
#' )
#'
#' # DiD power for proportions
#' power_did(
#'   treat = c(0.30, 0.36), control = c(0.30, 0.33),
#'   outcome = "prop", effect = 0.03, n = 800, power = NULL
#' )
#'
#' # MDE for means with n = 500
#' power_did(
#'   treat = c(50, 55), control = c(50, 52),
#'   outcome = "mean", var = 100, effect = NULL, n = 500
#' )
#'
#' # Panel overlap reduces required n
#' power_did(
#'   treat = c(0.50, 0.55), control = c(0.50, 0.48),
#'   outcome = "prop", effect = 0.07, overlap = 0.5, rho = 0.6
#' )
#'
#' @export
power_did <- function(treat, ...) {
  if (!missing(treat)) {
    .res <- .dispatch_plan(treat, "treat", power_did.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("power_did")
}

#' @rdname power_did
#' @export
power_did.default <- function(
  treat,
  control,
  outcome = c("mean", "prop"),
  var = NULL,
  effect = NULL,
  n = NULL,
  power = 0.80,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  alternative = c("two.sided", "one.sided"),
  ratio = 1,
  overlap = 0,
  rho = 0,
  plan = NULL,
  ...
) {
  .plan <- .merge_plan_args(plan, power_did.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(power_did.default, c(.plan, list(...))))
  outcome <- match.arg(outcome)
  alternative <- match.arg(alternative)

  treat <- .check_did_path(treat, "treat", outcome)
  control <- .check_did_path(control, "control", outcome)

  check_alpha(alpha)
  N_pair <- .check_power_N(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  check_overlap(overlap)
  check_rho(rho)

  ratio <- .resolve_ratio(n, ratio)

  null_count <- is.null(n) + is.null(power) + is.null(effect)
  if (null_count != 1L) {
    stop("leave exactly one of 'n', 'power', or 'effect' as NULL",
         call. = FALSE)
  }

  if (!is.null(n)) {
    n <- .check_power_n(n)
    .check_overlap_n(overlap, n = n)
  } else {
    .check_overlap_n(overlap, ratio = ratio)
  }
  if (!is.null(power)) check_proportion(power, "power")
  if (!is.null(effect)) check_scalar(effect, "effect")

  if (outcome == "mean") {
    var_parts <- .did_var_parts(var)
    var_terms <- .did_var_terms_mean(var_parts, overlap, rho)
  } else {
    var_parts <- NULL
    var_terms <- .did_var_terms_prop(treat, control, overlap, rho)
  }

  type <- if (outcome == "prop") "did_prop" else "did_mean"

  params <- list(
    treat = treat,
    control = control,
    outcome = outcome,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    alternative = alternative,
    ratio = ratio,
    overlap = overlap,
    rho = rho
  )

  if (!is.null(var_parts)) params$var <- var_parts

  if (is.null(n)) {
    params$effect <- effect
    params$power <- power

    z_a <- .z_alpha(alpha, alternative)
    z_b <- qnorm(power)

    if (all(is.infinite(N_pair))) {
      if (ratio == 1) {
        n0 <- (z_a + z_b)^2 * deff * sum(var_terms) / effect^2
        n0 <- n0 / resp_rate
      } else {
        n_c <- (z_a + z_b)^2 * deff *
          (var_terms[1] / ratio + var_terms[2]) / effect^2
        n_c <- n_c / resp_rate
        n0 <- c(ratio * n_c, n_c)
      }
    } else {
      r <- ratio
      power_n2 <- function(n2) {
        n_vec <- if (r == 1) c(n2, n2) else c(r * n2, n2)
        n_eff <- n_vec * resp_rate
        .power_did_power_core(
          effect = effect,
          n_eff = n_eff,
          var_terms = var_terms,
          alpha = alpha,
          N_pair = N_pair,
          deff = deff,
          alternative = alternative
        )
      }
      n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
      n0 <- if (r == 1) n2 else c(r * n2, n2)
    }

    .new_svyplan_power(
      n = n0, power = power, effect = effect,
      type = type, solved = "n", params = params
    )

  } else if (is.null(power)) {
    params$effect <- effect
    params$n <- n
    n_vec <- if (length(n) == 1L) c(n, n) else n
    n_eff <- n_vec * resp_rate

    pw <- .power_did_power_core(
      effect = effect,
      n_eff = n_eff,
      var_terms = var_terms,
      alpha = alpha,
      N_pair = N_pair,
      deff = deff,
      alternative = alternative
    )

    .new_svyplan_power(
      n = n, power = pw, effect = effect,
      type = type, solved = "power", params = params
    )

  } else {
    params$n <- n
    params$power <- power
    n_vec <- if (length(n) == 1L) c(n, n) else n
    n_eff <- n_vec * resp_rate

    mde <- .power_did_mde_core(
      n_eff = n_eff,
      power = power,
      alpha = alpha,
      N_pair = N_pair,
      deff = deff,
      alternative = alternative,
      var_terms = var_terms
    )

    .new_svyplan_power(
      n = n, power = power, effect = mde,
      type = type, solved = "mde", params = params
    )
  }
}

#' @keywords internal
#' @noRd
.check_did_path <- function(x, name, outcome) {
  if (!is.numeric(x) || length(x) != 2L || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("'%s' must be a finite numeric vector of length 2", name),
         call. = FALSE)
  }
  if (outcome == "prop") {
    if (any(x <= 0 | x >= 1)) {
      stop(sprintf("all '%s' values must be in (0, 1)", name), call. = FALSE)
    }
  }
  x
}

#' @keywords internal
#' @noRd
.did_var_parts <- function(var) {
  if (is.null(var)) {
    stop("'var' is required when outcome = 'mean'", call. = FALSE)
  }
  if (!is.numeric(var) || anyNA(var) || any(!is.finite(var))) {
    stop("'var' must be finite numeric", call. = FALSE)
  }
  if (!length(var) %in% c(1L, 2L, 4L)) {
    stop("'var' must have length 1, 2, or 4", call. = FALSE)
  }
  if (any(var <= 0)) {
    stop("all 'var' values must be positive", call. = FALSE)
  }
  if (length(var) == 1L) return(rep(var, 4L))
  if (length(var) == 2L) return(c(var[1], var[1], var[2], var[2]))
  var
}

#' @keywords internal
#' @noRd
.did_var_terms_prop <- function(treat, control, overlap, rho) {
  t0 <- treat[1]; t1 <- treat[2]
  c0 <- control[1]; c1 <- control[2]

  vt <- t0 * (1 - t0) + t1 * (1 - t1) -
    2 * overlap * rho * sqrt(t0 * (1 - t0) * t1 * (1 - t1))
  vc <- c0 * (1 - c0) + c1 * (1 - c1) -
    2 * overlap * rho * sqrt(c0 * (1 - c0) * c1 * (1 - c1))

  c(
    .safe_variance(vt, "treated change variance"),
    .safe_variance(vc, "control change variance")
  )
}

#' @keywords internal
#' @noRd
.did_var_terms_mean <- function(var_parts, overlap, rho) {
  vt0 <- var_parts[1]; vt1 <- var_parts[2]
  vc0 <- var_parts[3]; vc1 <- var_parts[4]

  vt <- vt0 + vt1 - 2 * overlap * rho * sqrt(vt0 * vt1)
  vc <- vc0 + vc1 - 2 * overlap * rho * sqrt(vc0 * vc1)

  c(
    .safe_variance(vt, "treated change variance"),
    .safe_variance(vc, "control change variance")
  )
}

#' @keywords internal
#' @noRd
.power_did_power_core <- function(
  effect, n_eff, var_terms, alpha, N_pair, deff, alternative
) {
  if (length(n_eff) == 1L) n_eff <- c(n_eff, n_eff)
  if (any(!is.infinite(N_pair) & n_eff >= N_pair)) {
    warning(
      "effective sample size >= population size; returning power = 1 (census)",
      call. = FALSE
    )
    return(1)
  }

  z_a <- .z_alpha(alpha, alternative)
  fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
  fpc2 <- .fpc_factor(n_eff[2], N_pair[2])

  V_d <- deff * (
    var_terms[1] * fpc1 / n_eff[1] +
    var_terms[2] * fpc2 / n_eff[2]
  )
  V_d <- .safe_variance(V_d, "DiD variance")
  if (V_d == 0) return(1)

  se <- sqrt(V_d)
  delta <- abs(effect)
  pw <- pnorm(delta / se - z_a)
  if (alternative == "two.sided") {
    pw <- pw + pnorm(-delta / se - z_a)
  }
  min(pw, 1)
}

#' @keywords internal
#' @noRd
.power_did_mde_core <- function(
  n_eff, power, alpha, N_pair, deff, alternative, var_terms
) {
  if (length(n_eff) == 1L) n_eff <- c(n_eff, n_eff)
  fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
  fpc2 <- .fpc_factor(n_eff[2], N_pair[2])
  V_d <- deff * (
    var_terms[1] * fpc1 / n_eff[1] +
    var_terms[2] * fpc2 / n_eff[2]
  )
  V_d <- .safe_variance(V_d, "DiD variance")

  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  (z_a + z_b) * sqrt(V_d)
}
