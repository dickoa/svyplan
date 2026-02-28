#' Power Analysis for Proportions
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-sample test of proportions. Leave exactly one of `n`, `power`, or
#' `p2` as `NULL` to solve for that quantity.
#'
#' @param p1 Baseline proportion, in (0, 1).
#' @param p2 Alternative proportion, in (0, 1). Leave `NULL` to solve for MDE.
#' @param n Per-group sample size. Leave `NULL` to solve for sample size.
#' @param power Target power, in (0, 1). Leave `NULL` to solve for power.
#' @param alpha Significance level, default 0.05.
#' @param N Per-group population size. `Inf` (default) means no finite
#'   population correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param sides `1` for one-sided or `2` (default) for two-sided test.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The sample size is inflated by `1 / resp_rate`.
#' @param overlap Panel overlap fraction in \[0, 1\], for repeated surveys.
#' @param rho Correlation between occasions in \[0, 1\].
#'
#' @return A `svyplan_power` object with components:
#' \describe{
#'   \item{`n`}{Per-group sample size.}
#'   \item{`power`}{Achieved power.}
#'   \item{`effect`}{Difference in proportions (`abs(p2 - p1)`).}
#'   \item{`solved`}{Which quantity was solved for
#'     (`"n"`, `"power"`, or `"mde"`).}
#'   \item{`params`}{List of input parameters.}
#' }
#'
#' @details
#' The effective variance for the difference in proportions is:
#'
#' \deqn{V = p_1(1-p_1) + p_2(1-p_2) - 2 \cdot \text{overlap} \cdot \rho
#'   \sqrt{p_1(1-p_1) \cdot p_2(1-p_2)}}
#'
#' The panel overlap term reduces the variance when the same units are
#' observed at both occasions, following the approach in Kish (1965, Ch. 11).
#'
#' - **Solve n**: \eqn{n_0 = (z_{\alpha/s} + z_\beta)^2 \cdot V \cdot
#'   \text{deff} / \delta^2}, then finite population correction.
#' - **Solve power**: Compute SE, then
#'   \eqn{\text{power} = \Phi(\delta / \text{SE} - z_{\alpha/s})}.
#'   For two-sided tests both tails are included.
#' - **Solve MDE**: `uniroot` search for the `p2` closest to `p1` that
#'   achieves the target power. Both directions (`p2 > p1` and `p2 < p1`)
#'   are searched; the closer alternative is returned.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Kish, L. (1965). *Survey Sampling*. Wiley.
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
#' # MDE with n = 1000 (searches both directions, returns closest p2)
#' power_prop(p1 = 0.30, n = 1000)
#'
#' # MDE near boundary (downward alternative found automatically)
#' power_prop(p1 = 0.999, n = 100)
#'
#' # With design effect
#' power_prop(p1 = 0.30, p2 = 0.35, deff = 1.5)
#'
#' # Panel survey with 50% overlap
#' power_prop(p1 = 0.30, p2 = 0.35, overlap = 0.5, rho = 0.6)
#'
#' @export
power_prop <- function(
  p1,
  p2 = NULL,
  n = NULL,
  power = 0.80,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  sides = 2,
  overlap = 0,
  rho = 0
) {
  check_proportion(p1, "p1")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_sides(sides)
  check_overlap(overlap)
  check_rho(rho)
  check_resp_rate(resp_rate)

  null_count <- is.null(p2) + is.null(n) + is.null(power)
  if (null_count != 1L) {
    stop("leave exactly one of 'n', 'power', or 'p2' as NULL", call. = FALSE)
  }

  if (!is.null(p2)) {
    check_proportion(p2, "p2")
  }
  if (!is.null(n)) {
    check_scalar(n, "n")
    if (n < 2) stop("'n' must be >= 2", call. = FALSE)
  }
  if (!is.null(power)) {
    check_proportion(power, "power")
  }

  params <- list(
    p1 = p1,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    sides = sides,
    overlap = overlap,
    rho = rho
  )

  if (is.null(n)) {
    if (p1 == p2) {
      stop("'p1' and 'p2' must differ", call. = FALSE)
    }
    params$p2 <- p2
    params$power <- power
    res <- .power_prop_n(p1, p2, power, alpha, N, deff, sides, overlap, rho)
    res <- .apply_resp_rate(res, resp_rate)
    .new_svyplan_power(
      n = res,
      power = power,
      effect = abs(p1 - p2),
      type = "proportion",
      solved = "n",
      params = params
    )
  } else if (is.null(power)) {
    if (p1 == p2) {
      stop("'p1' and 'p2' must differ", call. = FALSE)
    }
    params$p2 <- p2
    params$n <- n
    n_eff <- n * resp_rate
    res <- .power_prop_power(p1, p2, n_eff, alpha, N, deff, sides, overlap, rho)
    .new_svyplan_power(
      n = n,
      power = res,
      effect = abs(p1 - p2),
      type = "proportion",
      solved = "power",
      params = params
    )
  } else {
    params$n <- n
    params$power <- power
    n_eff <- n * resp_rate
    res <- .power_prop_mde(
      p1,
      n_eff,
      power,
      alpha,
      N,
      deff,
      sides,
      overlap,
      rho
    )
    .new_svyplan_power(
      n = n,
      power = power,
      effect = abs(p1 - res),
      type = "proportion",
      solved = "mde",
      params = c(params, list(p2 = res))
    )
  }
}

#' @keywords internal
#' @noRd
.prop_var <- function(p1, p2, overlap, rho) {
  p1 *
    (1 - p1) +
    p2 * (1 - p2) -
    2 * overlap * rho * sqrt(p1 * (1 - p1) * p2 * (1 - p2))
}

#' @keywords internal
#' @noRd
.power_prop_n <- function(p1, p2, power, alpha, N, deff, sides, overlap, rho) {
  z_alpha <- qnorm(1 - alpha / sides)
  z_beta <- qnorm(power)
  V <- .prop_var(p1, p2, overlap, rho)
  delta <- abs(p1 - p2)
  n0 <- (z_alpha + z_beta)^2 * V * deff / delta^2
  .apply_fpc(n0, N)
}

#' @keywords internal
#' @noRd
.power_prop_power <- function(p1, p2, n, alpha, N, deff, sides, overlap, rho) {
  if (!is.infinite(N) && n >= N) {
    warning(
      "effective sample size (",
      round(n, 1),
      ") >= population size (",
      N,
      "); returning power = 1 (census)",
      call. = FALSE
    )
    return(1)
  }
  z_alpha <- qnorm(1 - alpha / sides)
  V <- .prop_var(p1, p2, overlap, rho)
  delta <- abs(p1 - p2)
  f <- if (is.infinite(N)) 0 else n / N
  se <- sqrt(V * deff * (1 - f) / n)
  pw <- pnorm(delta / se - z_alpha)
  if (sides == 2L || sides == 2) {
    pw <- pw + pnorm(-delta / se - z_alpha)
  }
  min(pw, 1)
}

#' @keywords internal
#' @noRd
.power_prop_mde <- function(p1, n, power, alpha, N, deff, sides, overlap, rho) {
  target_fn <- function(p2) {
    .power_prop_power(p1, p2, n, alpha, N, deff, sides, overlap, rho) - power
  }

  eps <- 1e-8
  roots <- list()

  # Try upward: p2 in (p1 + eps, 1 - eps)
  up_lo <- p1 + eps
  up_hi <- 1 - eps
  if (up_lo < up_hi) {
    f_lo <- target_fn(up_lo)
    f_hi <- target_fn(up_hi)
    if (is.finite(f_lo) && is.finite(f_hi)) {
      if (f_lo == 0) {
        roots <- c(roots, list(up_lo))
      } else if (f_hi == 0) {
        roots <- c(roots, list(up_hi))
      } else if (sign(f_lo) != sign(f_hi)) {
        r <- uniroot(target_fn, c(up_lo, up_hi), tol = eps)
        roots <- c(roots, list(r$root))
      }
    }
  }

  # Try downward: p2 in (eps, p1 - eps)
  dn_lo <- eps
  dn_hi <- p1 - eps
  if (dn_lo < dn_hi) {
    f_lo <- target_fn(dn_lo)
    f_hi <- target_fn(dn_hi)
    if (is.finite(f_lo) && is.finite(f_hi)) {
      if (f_hi == 0) {
        roots <- c(roots, list(dn_hi))
      } else if (f_lo == 0) {
        roots <- c(roots, list(dn_lo))
      } else if (sign(f_lo) != sign(f_hi)) {
        r <- uniroot(target_fn, c(dn_lo, dn_hi), tol = eps)
        roots <- c(roots, list(r$root))
      }
    }
  }

  if (length(roots) == 0L) {
    stop(
      "no detectable alternative exists for the given n and power",
      call. = FALSE
    )
  }

  # Return the p2 closest to p1 (smallest MDE)
  dists <- vapply(roots, function(r) abs(r - p1), numeric(1L))
  roots[[which.min(dists)]]
}
