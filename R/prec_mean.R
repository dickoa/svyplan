#' Sampling Precision for a Mean
#'
#' Compute the sampling error (SE, margin of error, CV) for estimating a
#' population mean given a sample size. This is the inverse of [n_mean()].
#'
#' @param var For the default method: population variance \eqn{S^2}.
#'   For `svyplan_n` objects: a sample size result from [n_mean()].
#' @param ... Additional arguments passed to methods.
#' @param n Sample size.
#' @param mu Population mean. Required for the CV component.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The effective sample size is deflated by `resp_rate`.
#'
#' @return A `svyplan_prec` object with components `$se`, `$moe`, and `$cv`.
#'   `$cv` is `NA` when `mu` is not provided.
#'
#' @details
#' Computes the standard error for the given sample size and design
#' parameters, then derives the margin of error and coefficient of
#' variation. The effective sample size is `n * resp_rate / deff`, with
#' optional finite population correction.
#'
#' @seealso [n_mean()] for the inverse (compute n from a precision target),
#'   [prec_prop()] for proportions.
#'
#' @examples
#' # Precision with n = 400
#' prec_mean(var = 100, n = 400, mu = 50)
#'
#' # Without mu (CV will be NA)
#' prec_mean(var = 100, n = 400)
#'
#' @export
prec_mean <- function(var, ...) UseMethod("prec_mean")

#' @rdname prec_mean
#' @export
prec_mean.default <- function(var, n, mu = NULL, alpha = 0.05, N = Inf,
                              deff = 1, resp_rate = 1, ...) {
  check_scalar(var, "var")
  check_scalar(n, "n")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  if (!is.null(mu)) check_scalar(mu, "mu")

  n_eff <- n * resp_rate / deff
  z <- qnorm(1 - alpha / 2)

  fpc <- if (is.infinite(N)) 1 else 1 - n_eff / N
  fpc <- .clamp_fpc(fpc, n_eff, N)

  se <- sqrt(var * fpc / n_eff)
  moe <- z * se
  cv_val <- if (!is.null(mu)) se / mu else NA_real_

  params <- list(var = var, n = n, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate)
  if (!is.null(mu)) params$mu <- mu

  .new_svyplan_prec(
    se     = se,
    moe    = moe,
    cv     = cv_val,
    type   = "mean",
    params = params
  )
}

#' @rdname prec_mean
#' @export
prec_mean.svyplan_n <- function(var, ...) {
  x <- var
  if (x$type != "mean") {
    stop("prec_mean requires a svyplan_n of type 'mean'", call. = FALSE)
  }
  par <- x$params
  prec_mean.default(
    var       = par$var,
    n         = x$n,
    mu        = par$mu,
    alpha     = par$alpha,
    N         = par$N,
    deff      = par$deff,
    resp_rate = par$resp_rate %||% 1
  )
}
