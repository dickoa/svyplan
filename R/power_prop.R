#' Power Analysis for Proportions
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
#'   (`overlap = n12 / n1`). Only supported with `method = "wald"`.
#' @param rho Correlation between occasions in \[0, 1\].
#' @param method Variance method: `"wald"` (default), `"arcsine"`, or
#'   `"logodds"`. Arcsine and log-odds transforms are variance-stabilizing
#'   and perform better for rare or extreme proportions (Valliant,
#'   2018, \ifelse{html}{\out{&sect;}}{\enc{§}{S}}4.3.4--4.3.5).
#' @param plan Optional [svyplan()] object providing design defaults.
#' @param ... Additional arguments passed to methods.
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
#' Three methods are available:
#'
#' \describe{
#'   \item{`"wald"`}{Standard Wald test. Supports panel overlap.}
#'   \item{`"arcsine"`}{Arcsine-transformed test. The arcsine variance
#'     \eqn{1/(4n)} is approximately constant in \eqn{p}, making it
#'     more accurate for rare proportions.}
#'   \item{`"logodds"`}{Log-odds (logit) test with separate null and
#'     alternative variances. Uses Valliant Eq (4.23)/(4.24).}
#' }
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapter 4.
#'
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' @seealso [power_mean()] for continuous outcomes, [n_prop()] for
#'   estimation precision.
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
power_prop.default <- function(p1, p2 = NULL, n = NULL, power = 0.80,
                       alpha = 0.05, N = Inf, deff = 1,
                       resp_rate = 1,
                       alternative = c("two.sided", "one.sided"),
                       ratio = 1,
                       overlap = 0, rho = 0,
                       method = c("wald", "arcsine", "logodds"),
                       plan = NULL, ...) {
  .plan <- .merge_plan_args(plan, power_prop.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(power_prop.default, c(.plan, list(...))))
  check_proportion(p1, "p1")
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  check_alpha(alpha)
  N_pair <- .check_power_N(N)
  check_deff(deff)
  check_overlap(overlap)
  check_rho(rho)
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

  params <- list(p1 = p1, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate, alternative = alternative,
                 ratio = ratio, overlap = overlap, rho = rho,
                 method = method)

  if (is.null(n)) {
    if (p1 == p2) stop("'p1' and 'p2' must differ", call. = FALSE)
    params$p2 <- p2
    params$power <- power

    res <- switch(method,
      wald    = .power_prop_n_wald(p1, p2, power, alpha, N_pair, deff,
                                  alternative, overlap, rho, ratio, resp_rate),
      arcsine = .power_prop_n_arcsine(p1, p2, power, alpha, N_pair, deff,
                                      alternative, ratio, resp_rate),
      logodds = .power_prop_n_logodds(p1, p2, power, alpha, N_pair, deff,
                                      alternative, ratio, resp_rate))

    .new_svyplan_power(n = res, power = power, effect = abs(p1 - p2),
                       type = "proportion", solved = "n", params = params)

  } else if (is.null(power)) {
    if (p1 == p2) stop("'p1' and 'p2' must differ", call. = FALSE)
    params$p2 <- p2
    params$n <- n
    n_eff <- n * resp_rate

    res <- switch(method,
      wald    = .power_prop_power_wald(p1, p2, n_eff, alpha, N_pair, deff,
                                      alternative, overlap, rho),
      arcsine = .power_prop_power_arcsine(p1, p2, n_eff, alpha, N_pair, deff,
                                          alternative),
      logodds = .power_prop_power_logodds(p1, p2, n_eff, alpha, N_pair, deff,
                                          alternative))

    .new_svyplan_power(n = n, power = res, effect = abs(p1 - p2),
                       type = "proportion", solved = "power", params = params)

  } else {
    params$n <- n
    params$power <- power
    n_eff <- n * resp_rate

    res <- switch(method,
      wald    = .power_prop_mde_wald(p1, n_eff, power, alpha, N_pair, deff,
                                    alternative, overlap, rho),
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
.prop_var <- function(p1, p2, overlap, rho) {
  p1 * (1 - p1) + p2 * (1 - p2) -
    2 * overlap * rho * sqrt(p1 * (1 - p1) * p2 * (1 - p2))
}

# --- Wald internals ---

.power_prop_n_wald <- function(p1, p2, power, alpha, N_pair, deff,
                               alternative, overlap, rho, ratio, resp_rate) {
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  delta <- abs(p1 - p2)
  q1 <- 1 - p1; q2 <- 1 - p2
  ov_term <- 2 * overlap * rho * sqrt(p1 * q1 * p2 * q2)

  if (all(is.infinite(N_pair))) {
    if (ratio == 1) {
      V <- p1 * q1 + p2 * q2 - ov_term
      n2_0 <- (z_a + z_b)^2 * V * deff / delta^2
      n0 <- n2_0 / resp_rate
    } else {
      r <- ratio
      V_r <- p1 * q1 / r + p2 * q2 - ov_term
      n2 <- (z_a + z_b)^2 * V_r * deff / delta^2
      n2 <- n2 / resp_rate
      n0 <- c(r * n2, n2)
    }
  } else {
    r <- ratio
    power_n2 <- function(n2) {
      n_eff <- if (r == 1) c(n2, n2) else c(r * n2, n2)
      n_eff <- n_eff * resp_rate
      .power_prop_power_wald(p1, p2, n_eff, alpha, N_pair, deff,
                             alternative, overlap, rho)
    }
    n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
    n0 <- if (r == 1) n2 else c(r * n2, n2)
  }
  n0
}

.power_prop_power_wald <- function(p1, p2, n_eff, alpha, N_pair, deff,
                                   alternative, overlap, rho) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  if (any(!is.infinite(N_pair) & n_vec >= N_pair)) {
    warning("effective sample size >= population size; returning power = 1 (census)",
            call. = FALSE)
    return(1)
  }

  z_a <- .z_alpha(alpha, alternative)
  q1 <- 1 - p1; q2 <- 1 - p2
  delta <- abs(p1 - p2)

  fpc1 <- .fpc_factor(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor(n_vec[2], N_pair[2])
  V_d <- deff * (p1 * q1 * fpc1 / n_vec[1] + p2 * q2 * fpc2 / n_vec[2])
  if (overlap > 0) {
    fpc_ov <- .fpc_factor(n_vec[2], N_pair[2])
    V_d <- V_d - 2 * overlap * rho * sqrt(p1 * q1 * p2 * q2) *
           deff * fpc_ov / n_vec[2]
  }
  V_d <- .safe_variance(V_d, "difference variance")
  if (V_d == 0) return(1)

  se <- sqrt(V_d)
  pw <- pnorm(delta / se - z_a)
  if (alternative == "two.sided")
    pw <- pw + pnorm(-delta / se - z_a)
  min(pw, 1)
}

.power_prop_mde_wald <- function(p1, n_eff, power, alpha, N_pair, deff,
                                  alternative, overlap, rho) {
  target_fn <- function(p2) {
    .power_prop_power_wald(p1, p2, n_eff, alpha, N_pair, deff,
                           alternative, overlap, rho) - power
  }

  eps <- 1e-8
  roots <- list()

  up_lo <- p1 + eps; up_hi <- 1 - eps
  if (up_lo < up_hi) {
    f_lo <- target_fn(up_lo); f_hi <- target_fn(up_hi)
    if (is.finite(f_lo) && is.finite(f_hi)) {
      if (f_lo == 0) roots <- c(roots, list(up_lo))
      else if (f_hi == 0) roots <- c(roots, list(up_hi))
      else if (sign(f_lo) != sign(f_hi))
        roots <- c(roots, list(uniroot(target_fn, c(up_lo, up_hi), tol = eps)$root))
    }
  }

  dn_lo <- eps; dn_hi <- p1 - eps
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

# --- Arcsine internals ---

.power_prop_n_arcsine <- function(p1, p2, power, alpha, N_pair, deff,
                                   alternative, ratio, resp_rate) {
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  delta_phi <- asin(sqrt(p1)) - asin(sqrt(p2))

  if (all(is.infinite(N_pair))) {
    if (ratio == 1) {
      n2_0 <- ((z_a + z_b) / (sqrt(2) * abs(delta_phi)))^2 * deff
      n0 <- n2_0 / resp_rate
    } else {
      r <- ratio
      n2 <- ((z_a + z_b) / abs(delta_phi))^2 * deff * (1 / r + 1) / 4
      n2 <- n2 / resp_rate
      n0 <- c(r * n2, n2)
    }
  } else {
    r <- ratio
    power_n2 <- function(n2) {
      n_eff <- if (r == 1) c(n2, n2) else c(r * n2, n2)
      n_eff <- n_eff * resp_rate
      .power_prop_power_arcsine(p1, p2, n_eff, alpha, N_pair, deff,
                                alternative)
    }
    n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
    n0 <- if (r == 1) n2 else c(r * n2, n2)
  }
  n0
}

.power_prop_power_arcsine <- function(p1, p2, n_eff, alpha, N_pair, deff,
                                      alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  if (any(!is.infinite(N_pair) & n_vec >= N_pair)) {
    warning("effective sample size >= population size; returning power = 1 (census)",
            call. = FALSE)
    return(1)
  }

  z_a <- .z_alpha(alpha, alternative)
  delta_phi <- asin(sqrt(p1)) - asin(sqrt(p2))

  fpc1 <- .fpc_factor(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor(n_vec[2], N_pair[2])
  se_phi <- sqrt(deff * (fpc1 / (4 * n_vec[1]) + fpc2 / (4 * n_vec[2])))

  z_test <- abs(delta_phi) / se_phi - z_a
  pw <- pnorm(z_test)
  if (alternative == "two.sided")
    pw <- pw + pnorm(-abs(delta_phi) / se_phi - z_a)
  min(pw, 1)
}

.power_prop_mde_arcsine <- function(p1, n_eff, power, alpha, N_pair, deff,
                                     alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)

  fpc1 <- .fpc_factor(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor(n_vec[2], N_pair[2])
  se_phi <- sqrt(deff * (fpc1 / (4 * n_vec[1]) + fpc2 / (4 * n_vec[2])))
  mde_phi <- (z_a + z_b) * se_phi

  phi1 <- asin(sqrt(p1))
  p2_up <- sin(phi1 - mde_phi)^2
  p2_dn <- sin(phi1 + mde_phi)^2

  roots <- list()
  if (p2_up > 0 && p2_up < 1 && p2_up != p1) roots <- c(roots, list(p2_up))
  if (p2_dn > 0 && p2_dn < 1 && p2_dn != p1) roots <- c(roots, list(p2_dn))

  if (length(roots) == 0L)
    stop("no detectable alternative exists for the given n and power", call. = FALSE)

  dists <- vapply(roots, function(r) abs(r - p1), numeric(1L))
  roots[[which.min(dists)]]
}

# --- Log-odds internals ---

.power_prop_n_logodds <- function(p1, p2, power, alpha, N_pair, deff,
                                   alternative, ratio, resp_rate) {
  z_a <- .z_alpha(alpha, alternative)
  z_b <- qnorm(power)
  q1 <- 1 - p1; q2 <- 1 - p2
  delta_phi <- log(p1 / q1) - log(p2 / q2)
  p_bar <- (p1 + p2) / 2
  q_bar <- 1 - p_bar

  if (all(is.infinite(N_pair))) {
    if (ratio == 1) {
      V0 <- 2 / (p_bar * q_bar)
      VA <- 1 / (p1 * q1) + 1 / (p2 * q2)
      n2_0 <- ((z_a * sqrt(V0) + z_b * sqrt(VA)) / abs(delta_phi))^2 * deff
      n0 <- n2_0 / resp_rate
    } else {
      r <- ratio
      V0_coeff <- (1 / r + 1) / (p_bar * q_bar)
      VA_coeff <- 1 / (r * p1 * q1) + 1 / (p2 * q2)
      n2 <- ((z_a * sqrt(V0_coeff) + z_b * sqrt(VA_coeff)) / abs(delta_phi))^2 * deff
      n2 <- n2 / resp_rate
      n0 <- c(r * n2, n2)
    }
  } else {
    r <- ratio
    power_n2 <- function(n2) {
      n_eff <- if (r == 1) c(n2, n2) else c(r * n2, n2)
      n_eff <- n_eff * resp_rate
      .power_prop_power_logodds(p1, p2, n_eff, alpha, N_pair, deff,
                                alternative)
    }
    n2 <- .solve_n2_from_power(power, power_n2, N_pair, r, resp_rate)
    n0 <- if (r == 1) n2 else c(r * n2, n2)
  }
  n0
}

.power_prop_power_logodds <- function(p1, p2, n_eff, alpha, N_pair, deff,
                                      alternative) {
  n_vec <- if (length(n_eff) == 1L) c(n_eff, n_eff) else n_eff
  if (any(!is.infinite(N_pair) & n_vec >= N_pair)) {
    warning("effective sample size >= population size; returning power = 1 (census)",
            call. = FALSE)
    return(1)
  }

  z_a <- .z_alpha(alpha, alternative)
  q1 <- 1 - p1; q2 <- 1 - p2
  delta_phi <- log(p1 / q1) - log(p2 / q2)
  p_bar <- (p1 + p2) / 2
  q_bar <- 1 - p_bar

  fpc1 <- .fpc_factor(n_vec[1], N_pair[1])
  fpc2 <- .fpc_factor(n_vec[2], N_pair[2])

  V0 <- deff * (fpc1 / (n_vec[1] * p_bar * q_bar) +
                fpc2 / (n_vec[2] * p_bar * q_bar))
  VA <- deff * (fpc1 / (n_vec[1] * p1 * q1) +
                fpc2 / (n_vec[2] * p2 * q2))

  z_test <- (abs(delta_phi) - z_a * sqrt(V0)) / sqrt(VA)
  pw <- pnorm(z_test)
  if (alternative == "two.sided")
    pw <- pw + pnorm((-abs(delta_phi) - z_a * sqrt(V0)) / sqrt(VA))
  min(pw, 1)
}

.power_prop_mde_logodds <- function(p1, n_eff, power, alpha, N_pair, deff,
                                     alternative) {
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
