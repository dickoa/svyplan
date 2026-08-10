#' Power analysis for difference-in-differences designs
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-group, two-period difference-in-differences (DiD) contrast.
#'
#' `treat` and `control` already determine the contrast, so `effect` is
#' optional: leave `n` or `power` as `NULL` to solve for it. Leaving both
#' `n` and `power` supplied while `effect` is `NULL` solves for the MDE
#' instead. Exactly one quantity must remain unknown.
#'
#' @param treat Numeric length-2 vector for treated group outcomes:
#'   `c(baseline, endline)`.
#' @param control Numeric length-2 vector for control group outcomes:
#'   `c(baseline, endline)`.
#' @param outcome Outcome scale: `"mean"` (default) or `"prop"`.
#' @param var Outcome variance. Applies to `outcome = "mean"` only; under
#'   `outcome = "prop"` the cell variances follow from `treat` and `control`,
#'   so supplying it is an error rather than a silent no-op.
#'   Length 1: common variance for all four cells.
#'   Length 2: group-specific variances `c(var_treat, var_control)`,
#'   assumed equal across waves.
#'   Length 4: cell-specific variances in order
#'   `c(var_treat_baseline, var_treat_endline, var_control_baseline,
#'   var_control_endline)`.
#' @param sd Outcome standard deviation, an alternative spelling of `var`
#'   taking the same lengths. Supply at most one of `var` or `sd`.
#' @param effect Absolute DiD effect size to detect (> 0). Defaults to the
#'   contrast implied by `treat` and `control`,
#'   `|(treat[2] - treat[1]) - (control[2] - control[1])|`, so it only needs
#'   to be supplied to plan for an effect other than the one those paths
#'   describe. Supplying a value that disagrees with them warns. Leave `NULL`
#'   with both `n` and `power` given to solve for the MDE.
#' @param n Per-arm sample size per wave. Scalar (equal treated/control)
#'   or length-2 vector `c(n_treat, n_control)`. Leave `NULL` to solve
#'   for n.
#' @param power Target power, in (0, 1). Leave `NULL` to solve for power.
#' @param alpha Significance level, default 0.05.
#' @param N Population size for finite-population correction.
#'   A scalar applies to both arms. Use a length-2 vector
#'   `c(N_treat, N_control)` for arm-specific population sizes.
#'   `Inf` (default) disables FPC.
#' @param deff Design effect multiplier (> 0).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The required sample size is inflated by `1 / resp_rate`.
#' @param alternative Character: `"two.sided"` (default) or `"one.sided"`.
#'   Use `"two.sided"` when the intervention could plausibly move the
#'   outcome in either direction (the usual default). Use `"one.sided"`
#'   when only one direction is of interest, which requires a smaller
#'   sample for the same power.
#' @param ratio Allocation ratio `n_treat / n_control` (default 1).
#'   Used only when solving for `n` (`n = NULL`).
#' @param overlap Panel overlap fraction in \[0, 1\] within each arm
#'   across baseline and endline.
#' @param overlap_cor Correlation between baseline and endline outcomes within
#'   overlapping units, in \[0, 1\].
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
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
#'   \sqrt{V_{\text{pre}} \cdot V_{\text{post}}}}{V_arm = V_pre + V_post - 2 * overlap * rho * sqrt(V_pre * V_post)}
#'
#' The DiD test statistic variance is then:
#'
#' \deqn{V_d = \text{deff} \left(
#'   \frac{V_{\text{trt}} \cdot \text{fpc}_t}{n_t} +
#'   \frac{V_{\text{ctrl}} \cdot \text{fpc}_c}{n_c}
#' \right)}{V_d = deff ( (V_trt * fpc_t)/n_t + (V_ctrl * fpc_c)/n_c )}
#'
#' When `overlap = 0`, this reduces to the classical flat-variance formula.
#'
#' ## What `treat` and `control` are used for
#'
#' Both paths always supply the contrast. Whether they also supply the
#' variance depends on `outcome`: for `"prop"` the four cell variances are
#' \eqn{p(1-p)} at each path value, while for `"mean"` the variance comes
#' entirely from `var` and the paths are used for the contrast alone.
#'
#' Because the contrast is derived, `effect` is redundant with the paths and
#' is only needed to plan for a different effect than they describe, for
#' example a conservative target below the change a pilot observed. Doing so
#' warns, since the printed result then shows an `effect` that the displayed
#' paths do not produce.
#'
#' ## Normal approximation
#'
#' Critical values and power come from the standard normal, with the
#' variance treated as known and no degrees-of-freedom correction, as in
#' [power_mean()] and [power_prop()]. The cluster count, not the unit
#' count, is what governs the accuracy of that approximation for a
#' clustered DiD design: `deff` inflates the variance but does not model
#' the loss of degrees of freedom, so a design with few clusters per arm
#' is optimistic here by more than the unit count suggests.
#'
#' The `df` argument that [n_prop()], [n_mean()] and [n_alloc()] accept has
#' no counterpart here, and its absence is a decision rather than an
#' omission. There the quantile is the half-width of a confidence interval
#' and a t quantile substitutes for a normal one directly; here it is a
#' normal deviate for an alternative, and a t-based power calculation is a
#' different procedure. Passing `df` is an error that says so.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapter 4.
#'
#' @seealso [power_prop()] for two-sample proportions, [power_mean()] for
#'   two-sample means.
#'
#' @examples
#' # DiD sample size for means. The effect is the contrast the two paths
#' # describe: (55 - 50) - (52 - 50) = 3.
#' power_did(
#'   treat = c(50, 55), control = c(50, 52),
#'   outcome = "mean", var = 100
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
#'   outcome = "prop", effect = 0.07, overlap = 0.5, overlap_cor = 0.6
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
  ...,
  outcome = c("mean", "prop"),
  var = NULL,
  sd = NULL,
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
  overlap_cor = 0,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, power_did.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(power_did.default, c(.plan, list(...))))
  .stop_power_df(...)
  .check_unused_dots(...)
  outcome <- match.arg(outcome)
  if (outcome == "prop" && (!is.null(var) || !is.null(sd))) {
    stop(
      "'var' and 'sd' apply only to outcome = \"mean\"; under outcome = \"prop\" the cell variances follow from 'treat' and 'control'",
      call. = FALSE
    )
  }
  if (!is.null(sd)) var <- .resolve_var(var, sd)
  alternative <- match.arg(alternative)

  treat <- .check_did_path(treat, "treat", outcome)
  control <- .check_did_path(control, "control", outcome)

  check_alpha(alpha)
  N_pair <- .check_power_N(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  check_overlap(overlap)
  check_overlap_cor(overlap_cor)

  ratio <- .resolve_ratio(n, ratio)

  effect_implied <- (treat[2L] - treat[1L]) - (control[2L] - control[1L])

  null_count <- is.null(n) + is.null(power) + is.null(effect)
  # 'treat' and 'control' already determine the contrast, so an omitted
  # 'effect' is filled from them rather than being a second unknown. Only
  # when it is the second unknown: with both 'n' and 'power' supplied, a
  # NULL 'effect' still means "solve for the MDE".
  if (null_count == 2L && is.null(effect)) {
    if (abs(effect_implied) <= 0) {
      stop(
        "'treat' and 'control' imply a difference-in-differences of 0; supply 'effect', or leave both 'n' and 'power' to solve for the minimum detectable effect",
        call. = FALSE
      )
    }
    effect <- abs(effect_implied)
    null_count <- 1L
  }
  if (null_count != 1L) {
    stop("leave exactly one of 'n', 'power', or 'effect' as NULL",
         call. = FALSE)
  }
  if (!is.null(effect) &&
      abs(abs(effect) - abs(effect_implied)) >
        1e-8 * max(1, abs(effect_implied))) {
    warning(
      sprintf(
        "'effect' (%g) differs from the difference-in-differences implied by 'treat' and 'control' (%g); the supplied value is used for the contrast and 'treat'/'control' only for the variance",
        effect, abs(effect_implied)
      ),
      call. = FALSE
    )
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
    var_terms <- .did_var_terms_mean(var_parts, overlap, overlap_cor)
  } else {
    var_parts <- NULL
    var_terms <- .did_var_terms_prop(treat, control, overlap, overlap_cor)
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
    overlap_cor = overlap_cor
  )

  if (!is.null(var_parts)) params$var <- var_parts

  if (is.null(n)) {
    params$effect <- effect
    params$power <- power

    z_a <- .z_alpha(alpha, alternative)
    z_b <- qnorm(power)

    if (all(is.infinite(N_pair))) {
      if (ratio == 1) {
        n0 <- (z_a + z_b)^2 * deff * sum(var_terms$change) / effect^2
        n0 <- n0 / resp_rate
      } else {
        n_c <- (z_a + z_b)^2 * deff *
          (var_terms$change[1] / ratio + var_terms$change[2]) / effect^2
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
    .check_gross_n(n_vec, N_pair,
                   label = c("the treatment group", "the control group"))
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
    .check_gross_n(n_vec, N_pair,
                   label = c("the treatment group", "the control group"))
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

#' Per-arm variance of the before-after change
#'
#' Returns the per-unit variance and, alongside it, the census term the
#' finite population correction subtracts. The overlap covariance carries a
#' single \eqn{1/N} rather than the arm's marginal factor, so the change
#' variance splits as \eqn{v/n - v^{census}/N}{v/n - v^census/N} with the census term the
#' same expression at full overlap. See `.diff_var_fpc()`.
#' @keywords internal
#' @noRd
.did_var_pair <- function(v0, v1, overlap, overlap_cor, what) {
  cross <- overlap_cor * sqrt(v0 * v1)
  list(
    change = .safe_variance(v0 + v1 - 2 * overlap * cross, what),
    census = .safe_variance(v0 + v1 - 2 * cross, what)
  )
}

#' @keywords internal
#' @noRd
.did_var_terms_prop <- function(treat, control, overlap, overlap_cor) {
  t0 <- treat[1]; t1 <- treat[2]
  c0 <- control[1]; c1 <- control[2]

  vt <- .did_var_pair(t0 * (1 - t0), t1 * (1 - t1), overlap, overlap_cor,
                      "treated change variance")
  vc <- .did_var_pair(c0 * (1 - c0), c1 * (1 - c1), overlap, overlap_cor,
                      "control change variance")

  list(change = c(vt$change, vc$change), census = c(vt$census, vc$census))
}

#' @keywords internal
#' @noRd
.did_var_terms_mean <- function(var_parts, overlap, overlap_cor) {
  vt <- .did_var_pair(var_parts[1], var_parts[2], overlap, overlap_cor,
                      "treated change variance")
  vc <- .did_var_pair(var_parts[3], var_parts[4], overlap, overlap_cor,
                      "control change variance")

  list(change = c(vt$change, vc$change), census = c(vt$census, vc$census))
}

#' DiD variance with the correction applied to each arm
#' @keywords internal
#' @noRd
.did_var_d <- function(n_eff, var_terms, N_pair, deff) {
  if (length(n_eff) == 1L) n_eff <- c(n_eff, n_eff)
  per_unit <- var_terms$change / n_eff
  census <- ifelse(is.infinite(N_pair), 0, var_terms$census / N_pair)
  deff * sum(per_unit - census)
}

#' @keywords internal
#' @noRd
.power_did_power_core <- function(
  effect, n_eff, var_terms, alpha, N_pair, deff, alternative
) {
  z_a <- .z_alpha(alpha, alternative)
  V_d <- .did_var_d(n_eff, var_terms, N_pair, deff)
  V_d <- .safe_variance(V_d, "DiD variance")
  if (V_d == 0) return(1)

  se <- sqrt(V_d)
  icc <- abs(effect)
  pw <- pnorm(icc / se - z_a)
  if (alternative == "two.sided") {
    pw <- pw + pnorm(-icc / se - z_a)
  }
  min(pw, 1)
}

#' @keywords internal
#' @noRd
.power_did_mde_core <- function(
  n_eff, power, alpha, N_pair, deff, alternative, var_terms
) {
  V_d <- .did_var_d(n_eff, var_terms, N_pair, deff)
  V_d <- .safe_variance(V_d, "DiD variance")

  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  (z_a + z_b) * sqrt(V_d)
}
