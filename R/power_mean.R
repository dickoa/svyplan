#' Power Analysis for Means
#'
#' Compute sample size, power, or minimum detectable effect (MDE) for a
#' two-sample test of means. Leave exactly one of `n`, `power`, or `delta`
#' as `NULL` to solve for that quantity.
#'
#' @param delta Effect size (difference in means). Leave `NULL` to solve
#'   for MDE.
#' @param var Within-group variance (required).
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
#'   \item{`delta`}{Effect size (difference in means).}
#'   \item{`solved`}{Which quantity was solved for
#'     (`"n"`, `"power"`, or `"mde"`).}
#'   \item{`params`}{List of input parameters.}
#' }
#'
#' @details
#' The effective variance for the difference in means is:
#'
#' \deqn{V = 2 \cdot \text{var} \cdot (1 - \text{overlap} \cdot \rho)}
#'
#' The factor of 2 accounts for the two independent groups, reduced by
#' the panel overlap term (Kish, 1965, Ch. 11).
#'
#' - **Solve n**: \eqn{n_0 = (z_{\alpha/s} + z_\beta)^2 \cdot V \cdot
#'   \text{deff} / \delta^2}, then finite population correction.
#' - **Solve power**: Compute SE, then
#'   \eqn{\text{power} = \Phi(\delta / \text{SE} - z_{\alpha/s})}.
#' - **Solve MDE**: Analytical formula,
#'   \eqn{\delta = (z_{\alpha/s} + z_\beta) \sqrt{V \cdot \text{deff}
#'   \cdot (1-f) / n}}.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' @seealso [power_prop()] for proportions, [n_mean()] for estimation
#'   precision.
#'
#' @examples
#' # Sample size to detect a difference of 5 with variance 100
#' power_mean(delta = 5, var = 100)
#'
#' # Power given n = 200
#' power_mean(delta = 5, var = 100, n = 200, power = NULL)
#'
#' # MDE with n = 500
#' power_mean(var = 100, n = 500)
#'
#' # With design effect
#' power_mean(delta = 5, var = 100, deff = 1.5)
#'
#' @export
power_mean <- function(delta = NULL, var, n = NULL, power = 0.80,
                       alpha = 0.05, N = Inf, deff = 1,
                       resp_rate = 1,
                       sides = 2, overlap = 0, rho = 0) {
  check_scalar(var, "var")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_sides(sides)
  check_overlap(overlap)
  check_rho(rho)
  check_resp_rate(resp_rate)

  null_count <- is.null(delta) + is.null(n) + is.null(power)
  if (null_count != 1L) {
    stop("leave exactly one of 'n', 'power', or 'delta' as NULL", call. = FALSE)
  }

  if (!is.null(delta)) {
    check_scalar(delta, "delta")
  }
  if (!is.null(n)) {
    check_scalar(n, "n")
    if (n < 2) stop("'n' must be >= 2", call. = FALSE)
  }
  if (!is.null(power)) check_proportion(power, "power")

  V <- 2 * var * (1 - overlap * rho)

  params <- list(var = var, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate,
                 sides = sides, overlap = overlap, rho = rho)

  if (is.null(n)) {
    params$delta <- delta
    params$power <- power
    res <- .power_mean_n(delta, V, power, alpha, N, deff, sides)
    res <- .apply_resp_rate(res, resp_rate)
    .new_svyplan_power(n = res, power = power,
                       delta = delta, type = "mean",
                       solved = "n", params = params)

  } else if (is.null(power)) {
    params$delta <- delta
    params$n <- n
    n_eff <- n * resp_rate
    res <- .power_mean_power(delta, V, n_eff, alpha, N, deff, sides)
    .new_svyplan_power(n = n, power = res,
                       delta = delta, type = "mean",
                       solved = "power", params = params)

  } else {
    params$n <- n
    params$power <- power
    n_eff <- n * resp_rate
    res <- .power_mean_mde(V, n_eff, power, alpha, N, deff, sides)
    .new_svyplan_power(n = n, power = power,
                       delta = res, type = "mean",
                       solved = "mde", params = params)
  }
}

.power_mean_n <- function(delta, V, power, alpha, N, deff, sides) {
  z_alpha <- qnorm(1 - alpha / sides)
  z_beta <- qnorm(power)
  n0 <- (z_alpha + z_beta)^2 * V * deff / delta^2
  .apply_fpc(n0, N)
}

.power_mean_power <- function(delta, V, n, alpha, N, deff, sides) {
  if (!is.infinite(N) && n >= N) {
    warning("effective sample size (", round(n, 1),
            ") >= population size (", N,
            "); returning power = 1 (census)", call. = FALSE)
    return(1)
  }
  z_alpha <- qnorm(1 - alpha / sides)
  f <- if (is.infinite(N)) 0 else n / N
  se <- sqrt(V * deff * (1 - f) / n)
  pw <- pnorm(abs(delta) / se - z_alpha)
  if (sides == 2L || sides == 2) {
    pw <- pw + pnorm(-abs(delta) / se - z_alpha)
  }
  min(pw, 1)
}

.power_mean_mde <- function(V, n, power, alpha, N, deff, sides) {
  z_alpha <- qnorm(1 - alpha / sides)
  z_beta <- qnorm(power)
  f <- if (is.infinite(N)) 0 else n / N
  (z_alpha + z_beta) * sqrt(V * deff * (1 - f) / n)
}
