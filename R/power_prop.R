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
#' @param deff Design effect multiplier (>= 1).
#' @param sides `1` for one-sided or `2` (default) for two-sided test.
#' @param overlap Panel overlap fraction in \[0, 1\], for repeated surveys.
#' @param rho Correlation between occasions in \[0, 1\].
#'
#' @return A `svyplan_power` object.
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
#' - **Solve MDE**: `uniroot` search over `p2`.
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
#' # MDE with n = 1000
#' power_prop(p1 = 0.30, n = 1000)
#'
#' # With design effect
#' power_prop(p1 = 0.30, p2 = 0.35, deff = 1.5)
#'
#' # Panel survey with 50% overlap
#' power_prop(p1 = 0.30, p2 = 0.35, overlap = 0.5, rho = 0.6)
#'
#' @export
power_prop <- function(p1, p2 = NULL, n = NULL, power = 0.80,
                       alpha = 0.05, N = Inf, deff = 1,
                       sides = 2, overlap = 0, rho = 0) {
  check_proportion(p1, "p1")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_sides(sides)
  check_overlap(overlap)
  check_rho(rho)

  null_count <- is.null(p2) + is.null(n) + is.null(power)
  if (null_count != 1L) {
    stop("leave exactly one of 'n', 'power', or 'p2' as NULL", call. = FALSE)
  }

  if (!is.null(p2)) check_proportion(p2, "p2")
  if (!is.null(n)) {
    check_scalar(n, "n")
    if (n < 2) stop("'n' must be >= 2", call. = FALSE)
  }
  if (!is.null(power)) check_proportion(power, "power")

  params <- list(p1 = p1, alpha = alpha, N = N, deff = deff,
                 sides = sides, overlap = overlap, rho = rho)

  if (is.null(n)) {
    if (p1 == p2) stop("'p1' and 'p2' must differ", call. = FALSE)
    params$p2 <- p2
    params$power <- power
    res <- .power_prop_n(p1, p2, power, alpha, N, deff, sides, overlap, rho)
    .new_svyplan_power(n = res, power = power,
                       delta = abs(p1 - p2), type = "proportion",
                       solved = "n", params = params)

  } else if (is.null(power)) {
    if (p1 == p2) stop("'p1' and 'p2' must differ", call. = FALSE)
    params$p2 <- p2
    params$n <- n
    res <- .power_prop_power(p1, p2, n, alpha, N, deff, sides, overlap, rho)
    .new_svyplan_power(n = n, power = res,
                       delta = abs(p1 - p2), type = "proportion",
                       solved = "power", params = params)

  } else {
    params$n <- n
    params$power <- power
    res <- .power_prop_mde(p1, n, power, alpha, N, deff, sides, overlap, rho)
    .new_svyplan_power(n = n, power = power,
                       delta = abs(p1 - res), type = "proportion",
                       solved = "mde", params = c(params, list(p2 = res)))
  }
}

.prop_var <- function(p1, p2, overlap, rho) {
  p1 * (1 - p1) + p2 * (1 - p2) -
    2 * overlap * rho * sqrt(p1 * (1 - p1) * p2 * (1 - p2))
}

.power_prop_n <- function(p1, p2, power, alpha, N, deff, sides, overlap, rho) {
  z_alpha <- qnorm(1 - alpha / sides)
  z_beta <- qnorm(power)
  V <- .prop_var(p1, p2, overlap, rho)
  delta <- abs(p1 - p2)
  n0 <- (z_alpha + z_beta)^2 * V * deff / delta^2
  .apply_fpc(n0, N)
}

.power_prop_power <- function(p1, p2, n, alpha, N, deff, sides, overlap, rho) {
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

.power_prop_mde <- function(p1, n, power, alpha, N, deff, sides, overlap, rho) {
  target_fn <- function(p2) {
    .power_prop_power(p1, p2, n, alpha, N, deff, sides, overlap, rho) - power
  }

  upper <- min(1 - 1e-8, p1 + 0.5)
  lower <- p1 + 1e-8

  if (target_fn(upper) < 0) {
    stop("no detectable alternative exists for the given n and power",
         call. = FALSE)
  }

  res <- tryCatch(
    uniroot(target_fn, interval = c(lower, upper), tol = 1e-8),
    error = function(e) {
      stop("no detectable alternative exists for the given n and power",
           call. = FALSE)
    }
  )
  res$root
}
