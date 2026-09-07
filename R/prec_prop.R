#' Sampling precision for a proportion
#'
#' Compute the sampling error (SE, margin of error, CV) for estimating a
#' population proportion given a sample size. This is the inverse of
#' [n_prop()].
#'
#' @param p For the default method: expected proportion, in (0, 1). Leave
#'   it `NULL` and give `cv` or `rmoe` to solve for it instead.
#'   For `svyplan_n` objects: a sample size result from [n_prop()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param n Sample size, measured as gross units drawn and bounded by a
#'   finite `N`.
#' @param cv Target relative standard error, supplied instead of `p` to solve
#'   for the smallest proportion the design measures that precisely. Specify
#'   exactly one of `p`, `cv`, or `rmoe`.
#' @param rmoe Target margin of error relative to the proportion, `moe / p`,
#'   supplied instead of `p` to solve for the smallest proportion the design
#'   can report at that relative half-width. Unlike the `cv` solve this one
#'   reads the interval rather than the variance and so depends on `method`.
#'   Specify exactly one of `p`, `cv`, or `rmoe`.
#' @param alpha Significance level, default 0.05.
#' @param N Population size. `Inf` (default) means no finite population
#'   correction.
#' @param deff Design effect multiplier (> 0). Values < 1 are valid for
#'   efficient designs (e.g., stratified sampling with Neyman allocation).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1 (no
#'   adjustment). The effective sample size is `n * resp_rate`.
#' @param method One of `"wald"` (default), `"wilson"`, `"logodds"`, or
#'   `"beta"` (Korn-Graubard). See [n_prop()] for how to choose.
#' @param df Degrees of freedom of the variance estimator, typically sampled
#'   PSUs minus strata, and available from [design_df()]. It switches the
#'   interval quantile from normal to t under every method. `NULL`
#'   (default) applies no adjustment. See [n_prop()].
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "proportion"`:
#' \describe{
#'   \item{`se`}{Sampling standard error of the planned estimate, computed
#'     on the net sample `n * resp_rate`. It describes the estimator, so it
#'     is the same under all four methods and moves only with `p`, `n`,
#'     `deff`, `N`, and `resp_rate`.}
#'   \item{`moe`}{Half the length of the interval the chosen `method`
#'     builds at that standard error. Only under `"wald"` is this
#'     `q * se`, with `q` defined in [n_prop()]. Where [confint()]
#'     truncates a Wald interval at 0 or 1, `moe` stays `q * se` and
#'     exceeds half the reported width. The other three close at a
#'     different half-width, which is the whole of what distinguishes
#'     them. None of the three is symmetric about `p`, so `p - moe` and
#'     `p + moe` are not the endpoints. Use [confint()] for those.}
#'   \item{`cv`}{Relative standard error, `se / p`, and therefore also the
#'     same under all four methods.}
#'   \item{`rmoe`}{Margin of error relative to the proportion, `moe / p`.
#'     It reads the interval, so like `moe` it is method-specific.}
#'   \item{`method`}{The interval method used.}
#'   \item{`expected_cases`}{Positive cases the design expects to yield,
#'     `n * resp_rate * p`, the same quantity [n_prop()] reports and the
#'     one its `min_cases` constraint acts on.}
#'   \item{`solved`}{`"p"` when the proportion was solved for, and absent
#'     otherwise. The solved value is in `params$p`, whichever of `cv` or
#'     `rmoe` asked for it.}
#'   \item{`params`}{The validated inputs (`p`, `n`, `alpha`, `N`,
#'     `deff`, `resp_rate`, `df`, and `rmoe` when the proportion was solved
#'     from one). [predict()], [confint()] and the
#'     [n_prop()] round trip read the design back from here.}
#' }
#'
#' Nothing here is rounded: `n` is taken as given, so passing a
#' continuous `n` back from [n_prop()] reproduces its `se`, `moe`, `cv`
#' and `rmoe` exactly.
#'
#' @details
#' Computes the standard error for the given sample size and design
#' parameters, then derives the margin of error and coefficient of
#' variation. The variance equation is
#' `se^2 = deff * p * (1 - p) * fpc(n_net) / n_net`, where
#' `n_net = n * resp_rate` is the expected number of responding units:
#' `deff` multiplies the SRS variance at the realized sample size, and
#' the finite population correction uses the actual sampling fraction
#' `n_net / N`. Equivalently, `se^2 = p * (1 - p) / n_eff` with the
#' effective sample size `n_eff = n_net / (deff * fpc(n_net))`, the size at
#' which an infinite-population simple random sample would carry the same
#' variance.
#'
#' The FPC uses the Cochran (1977, Ch. 3) form. The finite-population
#' correction for a Bernoulli proportion is `(N - n_net) / (N - 1)`, not
#' the simpler `1 - n_net / N` used for means. All four methods read the
#' same variance through `n_eff`, so all four respond to `deff` and `N`
#' and all four return zero at a census. They differ only in the interval
#' built around that variance. See [n_prop()] for how to choose.
#'
#' That division of labour is worth stating plainly, because it decides what
#' each returned quantity answers. `se` and `cv` describe the estimator and
#' are identical under all four methods, and changing `method` on the same
#' design does not move them. `moe` and `rmoe` describe the interval and are
#' the only place the methods differ. Reading `se` back out of a half-width as
#' `moe / q` would make the reported standard error of the estimate change
#' when nothing but the interval construction had.
#'
#' ## The four precision quantities
#'
#' The package reports the same precision four ways, along two axes. One
#' axis is what is being described, the estimator's sampling variance or
#' the interval drawn around it. The other is whether it is stated in the
#' estimand's units or as a fraction of it.
#'
#' |            | Absolute | Relative |
#' |------------|----------|----------|
#' | Sampling   | `se`     | `cv`     |
#' | Interval   | `moe`    | `rmoe`   |
#'
#' The sampling row is free of `alpha` and of `method`. The interval row
#' depends on both: `alpha` fixes the confidence level and `method` fixes
#' the construction, and `rmoe = moe / p` inherits each. `rmoe = q * cv`
#' holds under `"wald"` alone, where `q` is the interval quantile the
#' design uses, normal by default and `t` on `df` degrees of freedom when
#' `df` is supplied. Under the other three the interval does not close at
#' `q * se`, and reading a relative margin of error off the CV that way
#' overstates the precision by 4 to 10 percent at `p = 0.02`, `n = 900`,
#' `deff = 2`.
#'
#' [n_prop()] takes any of `moe`, `rmoe`, and `cv` as the target, `cv`
#' under `"wald"` alone for the reason above, and all four quantities are
#' reported back by both functions.
#'
#' `$moe` is half the length of that interval, which equals `q * se` under
#' `"wald"` alone. Only the Wald interval is symmetric about `p`, so for the
#' other three `p - moe` and `p + moe` are not the limits: use [confint()]
#' for those. The `"beta"` method makes this most visible, since its interval
#' is deliberately asymmetric for a rare outcome.
#'
#' Supplying `cv` in place of `p` solves the same equation in the remaining
#' direction. [n_prop()] fixes the proportion and the precision and returns
#' the sample size. This function fixes the proportion and the sample size and
#' returns the precision. Giving `n` and `cv` returns the proportion. That
#' third reading answers which estimates a fielded design can carry. Because
#' `se(p) / p` falls as `p` rises, the solution is a floor, and every larger
#' proportion meets the target. It is the planning form of a publication
#' threshold. `expected_cases` is reported alongside it, so the precision
#' floor and the `min_cases` floor [n_prop()] applies can be read together.
#'
#' The solution is closed form, `p = 1 / (1 + n_eff * cv^2)`, and because
#' `cv` reads the sampling variance rather than an interval, it is the same
#' under all four methods: `method` changes the interval this design would
#' report, not the smallest proportion it can measure to a given relative
#' standard error. Every positive `cv` has a root in (0, 1), so no target is
#' unattainable. A census is an error, having no sampling variance for any
#' proportion to be measured against.
#'
#' Giving `n` and `rmoe` asks the same question of the interval instead of
#' the variance, and the answer is method-specific. It is the smallest
#' proportion the design can report at a relative margin of error of `rmoe`
#' under the chosen construction. There is no closed form, and two properties the
#' `cv` solve does not share bound it. The Wilson and Korn-Graubard
#' half-widths do not vanish as `p` approaches 1, so their relative margin
#' of error has a positive floor and a target below it is refused with the
#' floor named. The back-transformed log-odds half-width turns upward once
#' the logit spread outgrows `logit(p)`, near `p = 0.999` at `n = 1500` but
#' as low as `p = 0.97` at `n = 30`, so under `"logodds"` the target is met
#' on an interval rather than on every larger proportion, and the reported
#' root is its lower end.
#'
#' When called on a `svyplan_n` object, parameters are extracted from the
#' stored result. Any argument of the default method (e.g. `method`, `deff`,
#' `N`) can be overridden through `...`. Unknown argument names are an
#' error. Passing a different `method` evaluates the stored sample
#' size under that formula. The round-trip will not be exact because `n`
#' was determined under the original method.
#'
#' @family sample size and precision functions
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
#' # Korn-Graubard interval for a rare outcome in a clustered design
#' rare <- prec_prop(p = 0.02, n = 900, deff = 2, method = "beta", df = 25)
#' confint(rare)
#'
#' # Smallest proportion 1500 units can report at a 10 percent CV
#' floor <- prec_prop(n = 1500, cv = 0.10, N = 2e6)
#' floor$params$p
#'
#' # The same question of the interval: a 20 percent relative margin of
#' # error under the Korn-Graubard construction
#' prec_prop(n = 1500, rmoe = 0.20, method = "beta", N = 2e6)$params$p
#'
#' # How that floor moves with the sample size
#' predict(floor, expand.grid(n = c(500, 1500, 5000)))
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
  p = NULL,
  n,
  ...,
  cv = NULL,
  rmoe = NULL,
  alpha = 0.05,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  method = c("wald", "wilson", "logodds", "beta"),
  df = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, prec_prop.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_prop.default, c(.plan, list(...))))
  .check_unused_dots(...)
  if (is.null(p) + is.null(cv) + is.null(rmoe) != 2L) {
    stop("specify exactly one of 'p', 'cv', or 'rmoe'", call. = FALSE)
  }
  check_scalar(n, "n")
  check_alpha(alpha)
  check_population_size(N)
  check_deff(deff)
  check_resp_rate(resp_rate)
  .check_gross_n(n, N)
  method <- match.arg(method)
  if (!is.null(df)) check_df(df)

  solved <- if (is.null(p)) "p" else NULL
  if (!is.null(solved)) {
    if (!is.null(cv)) {
      check_scalar(cv, "cv")
      p <- .prec_solve_prop(cv, n, N, deff, resp_rate)
    } else {
      check_scalar(rmoe, "rmoe")
      p <- .prec_solve_prop_rmoe(rmoe, n, alpha, N, deff, resp_rate, method, df)
    }
  }
  check_proportion(p, "p")

  prec <- .prec_engine_prop(p, n, alpha, N, deff, resp_rate, method, df)
  se <- prec$se
  moe <- prec$moe
  cv_val <- prec$cv

  params <- list(
    p = p,
    n = n,
    alpha = alpha,
    N = N,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )
  # Which target drove the solve, since 'cv' and 'rmoe' both report
  # solved = "p" but re-solve differently: the CV root is the same under
  # every method and the relative-margin root is not.
  if (!is.null(rmoe)) {
    params$rmoe <- rmoe
  }

  .new_svyplan_prec(
    se = se,
    moe = moe,
    cv = cv_val,
    type = "proportion",
    method = method,
    params = params,
    solved = solved
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
    method = x$method %||% "wald",
    df = par$df
  )
  do.call(prec_prop.default, .roundtrip_args(args, list(...), prec_prop.default))
}
