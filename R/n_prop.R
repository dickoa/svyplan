#' Sample Size for a Proportion
#'
#' Compute the required sample size for estimating a population proportion
#' with a specified margin of error or coefficient of variation.
#'
#' @param p For the default method: expected proportion, in (0, 1).
#'   For `svyplan_prec` objects: a precision result from [prec_prop()].
#' @param ... Additional arguments passed to methods.
#' @param moe Desired margin of error. Specify exactly one of `moe` or `cv`.
#' @param cv Target coefficient of variation. Specify exactly one of `moe`
#'   or `cv`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The sample size is inflated by `1 / resp_rate`.
#' @param method One of `"wald"` (default), `"wilson"`, or `"logodds"`.
#'
#' @return A `svyplan_n` object.
#'
#' @details
#' Three confidence interval methods are available:
#'
#' - **Wald** (`"wald"`): Standard normal approximation
#'   (Cochran, 1977, Ch. 3). Supports both `moe` and `cv` modes,
#'   with optional finite population correction.
#' - **Wilson** (`"wilson"`): Wilson (1927) score interval. Only `moe`
#'   mode, no FPC.
#' - **Log-odds** (`"logodds"`): Log-odds (logit) transform interval.
#'   Only `moe` mode, with optional FPC.
#'
#' For proportions near 0 or 1 (below 0.1 or above 0.9), the Wald interval
#' has poor coverage; `method = "wilson"` is recommended in those cases.
#'
#' For the Wilson and log-odds methods, the design effect is applied as a
#' multiplicative factor to the final SRS sample size, which is an
#' approximation.
#'
#' The Wald FPC uses the Cochran (1977, Ch. 3) form with an `N/(N-1)` factor
#' to account for the Bernoulli finite-population variance.
#'
#' All methods use the normal (z) quantile. This is standard for survey
#' sampling where the sample size is large enough for the CLT to apply.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession,
#' and statistical inference. *Journal of the American Statistical
#' Association*, 22(158), 209--212.
#'
#' @seealso [n_mean()] for continuous variables, [n_cluster()] for
#'   multistage designs, [n_multi()] for multiple indicators,
#'   [prec_prop()] for the inverse.
#'
#' @examples
#' # Wald, absolute margin of error
#' n_prop(p = 0.3, moe = 0.05)
#'
#' # Wald, target CV with finite population
#' n_prop(p = 0.5, cv = 0.10, N = 10000)
#'
#' # Wilson score interval
#' n_prop(p = 0.1, moe = 0.03, method = "wilson")
#'
#' # With design effect and response rate
#' n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8)
#'
#' @export
n_prop <- function(p, ...) UseMethod("n_prop")

#' @rdname n_prop
#' @export
n_prop.default <- function(p, moe = NULL, cv = NULL, alpha = 0.05,
                           N = Inf, deff = 1, resp_rate = 1,
                           method = "wald", ...) {
  check_proportion(p, "p")
  check_precision(moe, cv)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  method <- match.arg(method, c("wald", "wilson", "logodds"))

  n <- switch(method,
    wald    = .n_prop_wald(p, moe, cv, alpha, N),
    wilson  = .n_prop_wilson(p, moe, cv, alpha),
    logodds = .n_prop_logodds(p, moe, cv, alpha, N)
  )

  n <- n * deff
  n <- .apply_resp_rate(n, resp_rate)

  if (!is.infinite(N) && n > N)
    warning("Calculated sample size exceeds population size N.", call. = FALSE)

  params <- list(p = p, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate)
  if (!is.null(moe)) params$moe <- moe else params$cv <- cv

  .new_svyplan_n(
    n      = n,
    type   = "proportion",
    method = method,
    params = params
  )
}

#' @keywords internal
#' @noRd
.n_prop_wald <- function(p, moe, cv, alpha, N) {
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p

  if (!is.null(moe)) {
    a <- if (is.infinite(N)) 1 else N / (N - 1)
    n0 <- a * z^2 * p * q / (moe^2 + z^2 * p * q / (N - 1))
  } else {
    a <- if (is.infinite(N)) 1 else N / (N - 1)
    n0 <- a * q / p / (cv^2 + q / p / (N - 1))
  }

  n0
}

#' @keywords internal
#' @noRd
.n_prop_wilson <- function(p, moe, cv, alpha) {
  if (is.null(moe)) {
    stop("Wilson method requires 'moe' (not 'cv')", call. = FALSE)
  }
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p
  e <- moe

  rad <- e^2 - p * q * (4 * e^2 - p * q)
  if (rad < 0) {
    stop("Wilson formula radicand is negative; check parameter values",
         call. = FALSE)
  }

  (p * q - 2 * e^2 + sqrt(rad)) * (z / e)^2 / 2
}

#' @keywords internal
#' @noRd
.n_prop_logodds <- function(p, moe, cv, alpha, N) {
  if (is.null(moe)) {
    stop("Log-odds method requires 'moe' (not 'cv')", call. = FALSE)
  }
  .n_prop_logodds_raw(p, moe, alpha, N)
}

#' Core log-odds n formula (no cv check).
#' Used by both .n_prop_logodds() and .logodds_moe().
#' @keywords internal
#' @noRd
.n_prop_logodds_raw <- function(p, e, alpha, N) {
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p
  kk <- q / p
  a <- if (is.infinite(N)) 1 else N / (N - 1)

  rad <- e^2 * (1 + kk^2)^2 + kk^2 * (1 - 2 * e) * (1 + 2 * e)
  x <- e * (1 + kk^2) + sqrt(rad)
  x <- x / (kk * (1 - 2 * e))

  inv_n <- a * (sqrt(p * q) / z * log(x))^2 + 1 / N
  1 / inv_n
}

#' @rdname n_prop
#' @export
n_prop.svyplan_prec <- function(p, moe = NULL, cv = NULL, ...) {
  x <- p
  if (x$type != "proportion") {
    stop("n_prop requires a svyplan_prec of type 'proportion'", call. = FALSE)
  }
  par <- x$params
  if (is.null(moe) && is.null(cv)) {
    moe <- x$moe
  }
  n_prop.default(
    p         = par$p,
    moe       = moe,
    cv        = cv,
    alpha     = par$alpha,
    N         = par$N,
    deff      = par$deff,
    resp_rate = par$resp_rate,
    method    = x$method %||% "wald"
  )
}
