#' Effective Sample Size
#'
#' Convert a planned sample size into the simple-random-sample size that
#' would give the same precision, `n * resp_rate / deff`. This is the
#' mirror of [design_effect()] and takes the same arguments.
#'
#' @param x A `svyplan_deff` from [design_effect()], a [n_cluster()] or
#'   [prec_cluster()] allocation, a [varcomp()] estimate, an [n_alloc()]
#'   allocation, or `NULL` (default) to build the design effect from the
#'   component arguments. A [n_prop()] or [n_mean()] result is not one of
#'   these: its design effect is the `deff` you supplied rather than
#'   something to be derived, so pass the two directly as
#'   `effective_n(n = x$n, deff = x$params$deff)`.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#'
#' @return A numeric scalar: the effective sample size,
#'   `n * resp_rate / deff`. It is not rounded, and it can exceed `n` when
#'   the design effect is below 1 (a stratification gain, for instance).
#'
#' @details
#' The design effect is built exactly as in [design_effect()], from
#' whichever of the clustering, weighting, and stratification components
#' you supply. `n` is the gross planned sample size; it is taken from the
#' arguments when it can be (the length of `weights`, the total of
#' `strata$n`, or the total of the plan you pass as `x`), and must be given
#' explicitly otherwise.
#'
#' Sizes are counted as units **issued**, so nonresponse has to be taken
#' off before the design effect is applied: a design that issues `n` and
#' analyses `n * resp_rate` of them carries the information of
#' `n * resp_rate / deff` simple random draws. This is the identity the
#' allocation and precision functions plan on, and the one reported in the
#' `n_eff` column of an [n_alloc()] table. Passing a plan as `x` picks up
#' the response rate it was built with, so `effective_n(plan)` and that
#' column agree; supply `resp_rate` yourself when you pass a bare `n` that
#' has not already been netted down.
#'
#' The design effect it divides by is a without-FPC planning quantity, so
#' for an [n_alloc()] result the effective size inherits that scale: it
#' answers how many simple random draws carry the same information under
#' the planning model, not how many the finite-population variance from
#' [prec_alloc()] would imply once sampling fractions are material. A
#' generalized allocation has no single design effect and is refused here
#' for the same reason it is in [design_effect()].
#'
#' `effective_n()` is a planning tool. To compute the effective sample size
#' realized by collected data, use `survey::svymean(..., deff = TRUE)` with
#' the actual weights, strata, and clusters.
#'
#' @seealso [design_effect()] for the design effect itself.
#'
#' @examples
#' # From a design effect you already built
#' effective_n(design_effect(icc = 0.05, n_per_psu = 25), n = 1200)
#'
#' # Straight from the components
#' effective_n(n = 1200, icc = 0.05, n_per_psu = 25)
#'
#' # From planned weights, where n is the number of units
#' effective_n(weights = rep(c(1, 4), c(300, 100)))
#'
#' # From a planned allocation
#' effective_n(strata = data.frame(N = c(50000, 120000), n = c(600, 400)))
#'
#' # From an optimized cluster plan
#' effective_n(n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05))
#'
#' @name effective_n
#' @export
effective_n <- function(x = NULL, ...) {
  # See design_effect(): dispatch must not fall through to the first ...
  # element when no object is supplied.
  if (missing(x) || is.null(x)) {
    return(effective_n.default(NULL, ...))
  }
  UseMethod("effective_n")
}

#' @describeIn effective_n Build the design effect from planning components.
#'
#' @param n Gross planned sample size, counted as units issued. Required
#'   unless it can be derived from `weights`, `strata$n`, or `x`.
#' @param deff Design effect to apply directly, instead of building one
#'   from components.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1. It nets
#'   `n` down to the units the design expects to analyse before the design
#'   effect is applied. When `x` is a plan, its own rate is used unless you
#'   override it here.
#' @param icc,n_per_psu,n_per_ssu,var_ratio,weights,strata Design components, with
#'   the same meaning as in [design_effect()].
#'
#' @export
effective_n.default <- function(
  x = NULL,
  ...,
  n = NULL,
  deff = NULL,
  resp_rate = 1,
  icc = NULL,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  var_ratio = 1,
  weights = NULL,
  strata = NULL
) {
  .check_unused_dots(...)
  if (!is.null(x)) {
    stop(
      "'x' must be a svyplan result or NULL. Pass planning components by name (n, deff, icc, n_per_psu, weights, strata)",
      call. = FALSE
    )
  }
  if (!is.null(deff)) {
    if (!is.null(icc) || !is.null(n_per_psu) || !is.null(n_per_ssu) ||
        !is.null(weights) || !is.null(strata)) {
      stop(
        "supply either 'deff' or the design components it would be built from",
        call. = FALSE
      )
    }
    check_deff(as.double(deff))
    value <- as.double(deff)
  } else {
    value <- as.double(design_effect(
      icc = icc, n_per_psu = n_per_psu, n_per_ssu = n_per_ssu, var_ratio = var_ratio,
      weights = weights, strata = strata
    ))
  }
  .effective_n(.effective_n_size(n, weights, strata), value, resp_rate)
}

#' @describeIn effective_n Apply a design effect already built by
#'   [design_effect()].
#' @export
effective_n.svyplan_deff <- function(x, ..., n = NULL, resp_rate = 1) {
  .check_unused_dots(...)
  .effective_n(.effective_n_size(n, NULL, NULL), as.double(x), resp_rate)
}

#' @describeIn effective_n Use the total sample size and design features of
#'   a [n_cluster()] or [prec_cluster()] allocation.
#' @export
effective_n.svyplan_cluster <- function(x, ..., n = NULL, resp_rate = NULL,
                                        weights = NULL, strata = NULL) {
  .check_unused_dots(...)
  deff <- design_effect(x, weights = weights, strata = strata)
  .effective_n(n %||% x$total_n, as.double(deff),
               .effective_n_resp(resp_rate, x))
}

#' @describeIn effective_n Use the total and design features of a
#'   [prec_cluster()] result. Other `svyplan_prec` types carry the `deff`
#'   you supplied rather than one to be derived, and are rejected.
#' @export
effective_n.svyplan_prec <- function(x, ..., n = NULL, resp_rate = NULL,
                                     weights = NULL, strata = NULL) {
  .check_unused_dots(...)
  deff <- design_effect(x, weights = weights, strata = strata)
  .effective_n(n %||% prod(x$params$n), as.double(deff),
               .effective_n_resp(resp_rate, x))
}

#' @describeIn effective_n Use the design features of a [varcomp()]
#'   estimate. `n` and the stage takes are your design choice.
#' @export
effective_n.svyplan_varcomp <- function(x, ..., n = NULL, resp_rate = 1,
                                        n_per_psu = NULL,
                                        n_per_ssu = NULL, weights = NULL,
                                        strata = NULL) {
  .check_unused_dots(...)
  deff <- design_effect(
    x, n_per_psu = n_per_psu, n_per_ssu = n_per_ssu,
    weights = weights, strata = strata
  )
  .effective_n(.effective_n_size(n, weights, strata), as.double(deff),
               resp_rate)
}

#' @describeIn effective_n Use the total and design features of an
#'   [n_alloc()] allocation, dividing by the allocation's own variance
#'   ratio. The frame needs `mean` alongside `N` and `sd`, since without
#'   the stratum means the population variance that ratio divides by is
#'   not identified.
#' @export
effective_n.svyplan_n <- function(x, ..., n = NULL, resp_rate = NULL,
                                  weights = NULL) {
  .check_unused_dots(...)
  deff <- design_effect(x, weights = weights)
  .effective_n(n %||% x$n, as.double(deff), .effective_n_resp(resp_rate, x))
}

#' Net a validated sample size down for response, then divide by a design
#' effect
#' @keywords internal
#' @noRd
.effective_n <- function(n, deff, resp_rate = 1) {
  check_scalar(n, "n")
  check_resp_rate(resp_rate)
  n * resp_rate / deff
}

#' Response rate for a plan, explicit argument first
#' @keywords internal
#' @noRd
.effective_n_resp <- function(resp_rate, x) {
  resp_rate %||% x$params$resp_rate %||% 1
}

#' Recover the gross sample size from whichever component supplies it
#' @keywords internal
#' @noRd
.effective_n_size <- function(n, weights, strata) {
  if (!is.null(n)) {
    return(n)
  }
  if (!is.null(weights)) {
    return(length(weights))
  }
  if (is.data.frame(strata) && !is.null(strata[["n"]])) {
    return(sum(strata$n))
  }
  stop(
    "'n' is required. It can only be derived from 'weights', from 'strata$n', or from a plan passed as the first argument",
    call. = FALSE
  )
}
