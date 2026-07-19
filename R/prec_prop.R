#' Sampling Precision for a Proportion
#'
#' Compute the sampling error (SE, margin of error, CV) for estimating a
#' population proportion given a sample size. This is the inverse of
#' [n_prop()].
#'
#' @param p For the default method: expected proportion, in (0, 1).
#'   For `svyplan_n` objects: a sample size result from [n_prop()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param n Sample size, measured as gross units drawn and bounded by a
#'   finite `N`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The effective sample size is `n * resp_rate`.
#' @param method One of `"wald"` (default), `"wilson"`, or `"logodds"`.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with components `$se`, `$moe`, and `$cv`.
#'
#' @details
#' Computes the standard error for the given sample size and design
#' parameters, then derives the margin of error and coefficient of
#' variation. The variance equation is
#' `se^2 = deff * p * (1 - p) * fpc(n_net) / n_net`, where
#' `n_net = n * resp_rate` is the expected number of responding units:
#' `deff` multiplies the SRS variance at the realized sample size, and
#' the finite population correction uses the actual sampling fraction
#' `n_net / N`. Equivalently, `se^2 = p * (1 - p) * fpc(n_net) / n_eff`
#' with the effective sample size `n_eff = n_net / deff`.
#'
#' The Wald FPC uses the Cochran (1977, Ch. 3) form: the finite-population
#' correction for a Bernoulli proportion is `(N - n_net) / (N - 1)`, not
#' the simpler `1 - n_net / N` used for means. The Wilson method has no
#' finite-population form. It is evaluated at `n_eff` and ignores `N`.
#'
#' When called on a `svyplan_n` object, parameters are extracted from the
#' stored result. Any argument of the default method (e.g. `method`, `deff`,
#' `N`) can be overridden through `...`. Unknown argument names are an
#' error. Passing a different `method` evaluates the stored sample
#' size under that formula. The round-trip will not be exact because `n`
#' was determined under the original method.
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
prec_prop <- function(p, ...) {
  if (!missing(p)) {
    .res <- .dispatch_plan(p, "p", prec_prop.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_prop")
}

#' @rdname prec_prop
#' @export
prec_prop.default <- function(
  p,
  n,
  ...,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  method = c("wald", "wilson", "logodds"),
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, prec_prop.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_prop.default, c(.plan, list(...))))
  .check_unused_dots(...)
  check_proportion(p, "p")
  check_scalar(n, "n")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  .check_gross_n(n, N)
  method <- match.arg(method)

  prec <- .prec_engine_prop(p, n, alpha, N, deff, resp_rate, method)
  se <- prec$se
  moe <- prec$moe
  cv_val <- prec$cv

  params <- list(
    p = p,
    n = n,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate
  )


  .new_svyplan_prec(
    se = se,
    moe = moe,
    cv = cv_val,
    type = "proportion",
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
  args <- list(
    p = par$p,
    n = x$n,
    alpha = par$alpha,
    N = par$N,
    deff = par$deff,
    resp_rate = par$resp_rate %||% 1,
    method = x$method %||% "wald"
  )
  do.call(prec_prop.default, .roundtrip_args(args, list(...), prec_prop.default))
}

#' Invert log-odds n formula to get MOE for a given net sample size.
#' @keywords internal
#' @noRd
.logodds_moe <- function(p, n_net, alpha, N, deff = 1) {
  if (!is.infinite(N) && n_net >= N) {
    warning("net sample size >= population size; moe is 0", call. = FALSE)
    return(0)
  }
  f <- function(e) .n_prop_logodds_raw(p, e, alpha, N, deff) - n_net
  lo <- 1e-6
  hi <- 0.5 - 1e-6
  f_lo <- f(lo)
  f_hi <- f(hi)
  if (sign(f_lo) == sign(f_hi)) {
    stop(
      "log-odds MOE inversion failed: no sign change on [",
      lo,
      ", ",
      hi,
      "] for p = ",
      p,
      ", n_net = ",
      round(n_net, 2),
      ", N = ",
      N,
      call. = FALSE
    )
  }
  tryCatch(
    uniroot(
      f,
      interval = c(lo, hi),
      f.lower = f_lo,
      f.upper = f_hi,
      tol = .Machine$double.eps^0.5
    )$root,
    error = function(e) {
      stop(
        "log-odds MOE inversion failed for p = ",
        p,
        ", n_net = ",
        round(n_net, 2),
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}
