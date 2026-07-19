#' Sample Size for a Proportion
#'
#' Compute the required sample size for estimating a population proportion
#' with a specified margin of error or coefficient of variation.
#'
#' @param p For the default method: expected proportion, in (0, 1).
#'   For `svyplan_prec` objects: a precision result from [prec_prop()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param moe Desired margin of error, the half-width of the confidence
#'   interval on the proportion scale. For example, `moe = 0.05` means
#'   the 95 percent CI should be no wider than +/- 5 percentage points.
#'   Specify exactly one of `moe` or `cv`.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.10` means the standard error should be at
#'   most 10 percent of the estimate. Use `cv` when you want precision
#'   to scale with the estimate (common in economic surveys). Use `moe`
#'   when you want a fixed absolute precision (common in health/DHS
#'   surveys). Specify exactly one of `moe` or `cv`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Accounts for the loss of
#'   precision from a complex design (clustering, unequal weights)
#'   compared to simple random sampling. A DEFF of 1.5 means 50 percent
#'   more interviews are needed for the same precision. Estimate from a
#'   previous survey, use [design_effect()] to compute it, or apply a
#'   rule of thumb (1.5--2.0 for typical cluster designs). Values < 1
#'   are valid for efficient designs (e.g., stratified sampling with
#'   Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The required sample size is inflated by `1 / resp_rate`.
#'   Estimate from response rates observed in similar surveys in the same
#'   population.
#' @param method One of `"wald"` (default), `"wilson"`, or `"logodds"`.
#' @param plan Optional [svyplan()] object providing design defaults.
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
#' has poor coverage. The recommended choice in those cases is
#' `method = "wilson"`.
#'
#' For the Wilson and log-odds methods, the design effect is applied as a
#' multiplicative factor to the final SRS sample size, which is an
#' approximation.
#'
#' ## Finite population correction
#'
#' Setting `N` to a finite value reduces the required sample size when
#' the sampling fraction (n/N) is non-negligible. As a rule of thumb,
#' FPC has little effect when n/N < 5 percent. The Wald FPC uses the
#' Cochran (1977, Ch. 3) form with an `N/(N-1)` factor to account for
#' the Bernoulli finite-population variance. This differs from
#' [n_mean()], where no `N/(N-1)` adjustment is needed because the
#' variance is already defined on `N-1` degrees of freedom.
#'
#' All methods use the normal (z) quantile. This is standard for survey
#' sampling where the sample size is large enough for the CLT to apply.
#'
#' When called on a `svyplan_prec` object, parameters are extracted from the
#' stored result. Any argument of the default method (e.g. `method`, `deff`,
#' `N`) can be overridden through `...`. Unknown argument names are an
#' error. Passing a different `method` evaluates the stored precision
#' target under that formula. The round-trip will not be exact because the
#' precision was computed under the original method.
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
#' # MICS/DHS-style relative margin of error (RME)
#' # RME = moe / p, so moe = RME * p
#' p <- 0.2
#' n_prop(p = p, moe = 0.12 * p, deff = 1.5, resp_rate = 0.9)
#'
#' @export
n_prop <- function(p, ...) {
  if (!missing(p)) {
    .res <- .dispatch_plan(p, "p", n_prop.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_prop")
}

#' @rdname n_prop
#' @export
n_prop.default <- function(
  p,
  ...,
  moe = NULL,
  cv = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  method = c("wald", "wilson", "logodds"),
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_prop.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_prop.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  check_proportion(p, "p")
  check_precision(moe, cv)
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  method <- match.arg(method)

  n <- switch(
    method,
    wald = .n_prop_wald(p, moe, cv, alpha, N, deff),
    wilson = .n_prop_wilson(p, moe, cv, alpha) * deff,
    logodds = .n_prop_logodds(p, moe, cv, alpha, N, deff)
  )

  n <- .apply_resp_rate(n, resp_rate)
  .check_attainable(n, N, resp_rate)

  params <- list(
    p = p,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate
  )

  if (!is.null(moe)) {
    params$moe <- moe
  } else {
    params$cv <- cv
  }

  .new_svyplan_n(
    n = n,
    type = "proportion",
    method = method,
    params = params
  )
}

#' @keywords internal
#' @noRd
.n_prop_wald <- function(p, moe, cv, alpha, N, deff = 1) {
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p
  a <- ifelse(is.infinite(N), 1, N / (N - 1))

  if (!is.null(moe)) {
    n0 <- a * z^2 * deff * p * q / (moe^2 + z^2 * deff * p * q / (N - 1))
  } else {
    n0 <- a * deff * q / p / (cv^2 + deff * q / p / (N - 1))
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
    stop(
      "Wilson formula radicand is negative; check parameter values",
      call. = FALSE
    )
  }

  (p * q - 2 * e^2 + sqrt(rad)) * (z / e)^2 / 2
}

#' @keywords internal
#' @noRd
.n_prop_logodds <- function(p, moe, cv, alpha, N, deff = 1) {
  if (is.null(moe)) {
    stop("Log-odds method requires 'moe' (not 'cv')", call. = FALSE)
  }
  if (moe >= 0.5) {
    stop("log-odds method requires 'moe' < 0.5", call. = FALSE)
  }
  .n_prop_logodds_raw(p, moe, alpha, N, deff)
}

#' Core log-odds n formula (no cv check).
#' Used by both .n_prop_logodds() and .logodds_moe().
#' @keywords internal
#' @noRd
.n_prop_logodds_raw <- function(p, e, alpha, N, deff = 1) {
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p
  kk <- q / p
  a <- ifelse(is.infinite(N), 1, N / (N - 1))

  rad <- e^2 * (1 + kk^2)^2 + kk^2 * (1 - 2 * e) * (1 + 2 * e)
  x <- e * (1 + kk^2) + sqrt(rad)
  x <- x / (kk * (1 - 2 * e))

  inv_n <- a * (sqrt(p * q) / z * log(x))^2 / deff + 1 / N
  1 / inv_n
}

#' @rdname n_prop
#' @export
n_prop.svyplan_prec <- function(p, ..., moe = NULL, cv = NULL) {
  x <- p
  if (x$type != "proportion") {
    stop("n_prop requires a svyplan_prec of type 'proportion'", call. = FALSE)
  }
  par <- x$params
  if (is.null(moe) && is.null(cv)) {
    moe <- x$moe
  }
  args <- list(
    p = par$p,
    moe = moe,
    cv = cv,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate,
    method = x$method %||% "wald"
  )
  do.call(n_prop.default, .roundtrip_args(args, list(...), n_prop.default))
}
