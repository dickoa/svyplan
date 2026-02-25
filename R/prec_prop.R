#' Sampling Precision for a Proportion
#'
#' Compute the sampling error (SE, margin of error, CV) for estimating a
#' population proportion given a sample size. This is the inverse of
#' [n_prop()].
#'
#' @param p For the default method: expected proportion, in (0, 1).
#'   For `svyplan_n` objects: a sample size result from [n_prop()].
#' @param ... Additional arguments passed to methods.
#' @param n Sample size.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The effective sample size is deflated by `resp_rate`.
#' @param method One of `"wald"` (default), `"wilson"`, or `"logodds"`.
#'
#' @return A `svyplan_prec` object with components `$se`, `$moe`, and `$cv`.
#'
#' @details
#' Computes the standard error for the given sample size and design
#' parameters, then derives the margin of error and coefficient of
#' variation. The effective sample size is `n * resp_rate / deff`, with
#' optional finite population correction.
#'
#' The Wald FPC uses the Cochran (1977, Ch. 3) form: the finite-population
#' correction for a Bernoulli proportion is `(N - n_eff) / (N - 1)`, not
#' the simpler `1 - n_eff / N` used for means.
#'
#' @seealso [n_prop()] for the inverse (compute n from a precision target),
#'   [prec_mean()] for continuous variables.
#'
#' @examples
#' # Precision with n = 400
#' prec_prop(p = 0.3, n = 400)
#'
#' # With design effect and response rate
#' prec_prop(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.8)
#'
#' @export
prec_prop <- function(p, ...) UseMethod("prec_prop")

#' @rdname prec_prop
#' @export
prec_prop.default <- function(p, n, alpha = 0.05, N = Inf, deff = 1,
                              resp_rate = 1, method = "wald", ...) {
  check_proportion(p, "p")
  check_scalar(n, "n")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  method <- match.arg(method, c("wald", "wilson", "logodds"))

  n_eff <- n * resp_rate / deff
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p

  if (method == "wald") {
    fpc <- if (is.infinite(N)) 1 else (N - n_eff) / (N - 1)
    fpc <- .clamp_fpc(fpc, n_eff, N)
    se <- sqrt(p * q * fpc / n_eff)
    moe <- z * se
  } else if (method == "wilson") {
    moe <- z * sqrt(p * q / n_eff + z^2 / (4 * n_eff^2)) / (1 + z^2 / n_eff)
    se <- moe / z
  } else {
    moe <- .logodds_moe(p, n_eff, alpha, N)
    se <- moe / z
  }
  cv_val <- se / p

  params <- list(p = p, n = n, alpha = alpha, N = N, deff = deff,
                 resp_rate = resp_rate)

  .new_svyplan_prec(
    se     = se,
    moe    = moe,
    cv     = cv_val,
    type   = "proportion",
    method = method,
    params = params
  )
}

#' @rdname prec_prop
#' @export
prec_prop.svyplan_n <- function(p, ...) {
  x <- p
  if (x$type != "proportion") {
    stop("prec_prop requires a svyplan_n of type 'proportion'", call. = FALSE)
  }
  par <- x$params
  prec_prop.default(
    p         = par$p,
    n         = x$n,
    alpha     = par$alpha,
    N         = par$N,
    deff      = par$deff,
    resp_rate = par$resp_rate %||% 1,
    method    = x$method %||% "wald"
  )
}

#' Invert log-odds n formula to get MOE for a given n_eff.
#' @keywords internal
#' @noRd
.logodds_moe <- function(p, n_eff, alpha, N) {
  uniroot(
    function(e) .n_prop_logodds_raw(p, e, alpha, N) - n_eff,
    interval = c(1e-6, p * (1 - p)),
    tol = .Machine$double.eps^0.5
  )$root
}
