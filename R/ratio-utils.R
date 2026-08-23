#' Check the correlation between a ratio's numerator and denominator
#'
#' Distinct from `check_overlap_cor()`, which bounds a repeated-occasion
#' correlation at zero. Two variables measured on the same unit may be
#' negatively correlated, so the admissible range is the full one.
#' @keywords internal
#' @noRd
check_component_cor <- function(component_cor, name = "component_cor") {
  if (
    !is.numeric(component_cor) ||
      length(component_cor) != 1L ||
      anyNA(component_cor) ||
      component_cor < -1 ||
      component_cor > 1
  ) {
    stop(sprintf("'%s' must be a number in [-1, 1]", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the four moments a ratio estimand is planned from
#' @keywords internal
#' @noRd
.check_ratio_moments <- function(r, cv_num, cv_den, component_cor) {
  if (is.null(r) || is.null(cv_num) || is.null(cv_den) ||
        is.null(component_cor)) {
    stop(
      "a ratio needs 'r', 'cv_num', 'cv_den', and 'component_cor'",
      call. = FALSE
    )
  }
  check_mu(r, "r")
  check_scalar(cv_num, "cv_num")
  check_scalar(cv_den, "cv_den")
  check_component_cor(component_cor)
  invisible(TRUE)
}

#' Unit relative variance of a ratio of two totals
#'
#' The coefficient `L_R` of Var(R_hat)/R^2, equal to
#' `var(y - R * x) / mean(y)^2`. `sign(r)` carries the sign of
#' `mu_y * mu_x`, without which the covariance term is subtracted in the
#' wrong direction whenever the two means differ in sign.
#'
#' Cancellation is possible when the components are near-proportional, so a
#' result inside a relative tolerance of the terms that formed it is returned
#' as an exact zero. A value below its negative is unreachable from valid
#' moments, since `L_R >= (cv_num - cv_den)^2`, and is an internal error
#' rather than an input one.
#' @keywords internal
#' @noRd
.ratio_unit_relvar <- function(r, cv_num, cv_den, component_cor) {
  cross <- cv_num * cv_den
  value <- cv_num^2 + cv_den^2 - 2 * component_cor * cross * sign(r)
  tol <- 8 * .Machine$double.eps *
    pmax(1, cv_num^2 + cv_den^2 + 2 * abs(component_cor) * cross)
  # A coefficient of variation large enough to overflow when squared makes the
  # tolerance infinite, and `abs(value) <= Inf` would then clamp a genuine Inf
  # to zero. Only a finite comparison can call a value cancellation rather
  # than signal, so a non-finite result is passed through for validation.
  decidable <- is.finite(value) & is.finite(tol)
  if (any(decidable & value < -tol)) {
    stop(
      "internal error: negative unit relative variance for a ratio",
      call. = FALSE
    )
  }
  ifelse(decidable & abs(value) <= tol, 0, value)
}

#' Refuse a ratio with no sampling variance, and flag one that nearly has none
#'
#' At `L_R = 0` the first-order model says `y` is proportional to `x`, every
#' positive sample estimates the ratio exactly, and no size is identified. Just
#' above zero the size is real but is dominated by the last digits of the
#' correlation, which is worth saying out loud rather than returning in
#' silence.
#' @keywords internal
#' @noRd
.check_ratio_relvar <- function(unit_relvar, cv_num, cv_den, component_cor) {
  if (!is.finite(unit_relvar)) {
    stop(
      "the ratio's unit relative variance is not finite at these moments: ",
      "'cv_num' and 'cv_den' are too large to square without overflow",
      call. = FALSE
    )
  }
  if (unit_relvar == 0) {
    stop(
      "the ratio has no sampling variance at these moments: 'cv_num' and ",
      "'cv_den' are equal and 'component_cor' is at its bound, which makes ",
      "the numerator proportional to the denominator and the ratio constant. ",
      "No sample size is identified",
      call. = FALSE
    )
  }
  if (unit_relvar < 1e-3 * max(cv_num^2, cv_den^2)) {
    warning(
      "the ratio's unit relative variance (",
      signif(unit_relvar, 3),
      ") is far below its components': the plan rests on 'component_cor' = ",
      signif(component_cor, 4),
      " being known to that precision. Vary it with predict() before ",
      "committing to this size",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Precision of a ratio estimate at a gross sample size
#'
#' Numerically the mean engine at `var = r^2 * unit_relvar` and `mu = r`, so
#' the finite population correction, response inflation, census clamp and
#' interval quantile are the ones every other estimand uses. Only ratio-scale
#' quantities are returned: the equivalent variance is a step in the
#' calculation and is not a variance of either observed component.
#' @keywords internal
#' @noRd
.prec_engine_ratio <- function(r, unit_relvar, n, alpha, N, deff, resp_rate,
                               df = NULL) {
  # The relative standard error first, then the scale. Forming r^2 overflows
  # for a large finite ratio and underflows for a small one, in both cases
  # losing a standard error the ratio scale represents perfectly well: at
  # r = 1e200 the answer is about 4e198, not Inf.
  rel <- .prec_engine_mean(unit_relvar, 1, n, alpha, N, deff, resp_rate, df)
  list(se = abs(r) * rel$se, moe = abs(r) * rel$moe, cv = rel$cv)
}

#' The relative standard error a ratio target implies
#'
#' All three targets reduce to one, exactly rather than approximately: a
#' margin of error on the ratio scale is `q * abs(r)` times its relative
#' standard error.
#' @keywords internal
#' @noRd
.ratio_target_cv <- function(r, moe, cv, rmoe, q) {
  moe_used <- if (is.null(rmoe)) moe else .moe_from_rmoe(rmoe, r, "r")
  if (is.null(moe_used)) cv else moe_used / (q * abs(r))
}

#' Rows of an indicator table that declare a ratio estimand
#'
#' A non-missing `r` is the marker. `unit_relvar` cannot serve as one: it is
#' already a valid column for cluster planning of means and proportions and
#' says nothing about the estimand's scale.
#' @keywords internal
#' @noRd
.is_ratio_row <- function(indicators) {
  if (!"r" %in% names(indicators)) {
    return(rep(FALSE, nrow(indicators)))
  }
  !is.na(indicators$r)
}

#' Validate the moment quartet on every row that declares a ratio
#'
#' Row by row, because a mixed table leaves the ratio columns missing on its
#' proportion and mean rows, and that is not an error there.
#' @keywords internal
#' @noRd
.check_ratio_rows <- function(indicators) {
  is_ratio <- .is_ratio_row(indicators)
  if (!any(is_ratio)) {
    return(invisible(TRUE))
  }
  absent <- setdiff(c("cv_num", "cv_den", "component_cor"), names(indicators))
  if (length(absent) > 0L) {
    stop(
      sprintf(
        "a ratio row needs 'r', 'cv_num', 'cv_den', and 'component_cor'; column(s) %s are absent",
        paste(sQuote(absent), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  for (i in which(is_ratio)) {
    moments <- c(
      cv_num = indicators$cv_num[i],
      cv_den = indicators$cv_den[i],
      component_cor = indicators$component_cor[i]
    )
    if (anyNA(moments)) {
      stop(
        sprintf(
          "row %d sets 'r' but leaves %s missing",
          i, paste(sQuote(names(moments)[is.na(moments)]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    .check_ratio_moments(indicators$r[i], moments[["cv_num"]],
                         moments[["cv_den"]], moments[["component_cor"]])
    relvar <- .ratio_unit_relvar(indicators$r[i], moments[["cv_num"]],
                                 moments[["cv_den"]],
                                 moments[["component_cor"]])
    .check_ratio_relvar(relvar, moments[["cv_num"]], moments[["cv_den"]],
                        moments[["component_cor"]])
  }
  if ("unit_relvar" %in% names(indicators)) {
    # A ratio row derives this from its moments, so a value alongside them can
    # only agree or contradict. It agrees whenever the table came back from a
    # result, which stores the derived coefficient, so the check is on the
    # value rather than on its presence.
    for (i in which(is_ratio & !is.na(indicators$unit_relvar))) {
      derived <- .ratio_unit_relvar(
        indicators$r[i], indicators$cv_num[i], indicators$cv_den[i],
        indicators$component_cor[i]
      )
      supplied <- indicators$unit_relvar[i]
      if (!isTRUE(all.equal(supplied, derived, tolerance = 1e-8))) {
        stop(
          sprintf(
            paste0(
              "row %d gives 'unit_relvar' = %s, but its moments imply %s. A ",
              "ratio row derives its unit relative variance from 'cv_num', ",
              "'cv_den' and 'component_cor', so drop 'unit_relvar' on that ",
              "row or correct the moments"
            ),
            i, format(supplied), format(derived)
          ),
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

#' Gross sample size a ratio target needs
#' @keywords internal
#' @noRd
.n_ratio_from_target <- function(r, unit_relvar, moe, cv, rmoe, alpha, N,
                                 deff, resp_rate, df = NULL) {
  q <- .q_alpha(alpha, df)
  target_cv <- .ratio_target_cv(r, moe, cv, rmoe, q)
  n <- deff * unit_relvar / (target_cv^2 + deff * unit_relvar / N)
  n <- .apply_resp_rate(n, resp_rate)
  .check_attainable(n, N, resp_rate)
  n
}
