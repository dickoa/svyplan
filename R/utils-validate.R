#' Reject arguments that were not consumed from dots
#'
#' Uses the dots metadata primitives so argument promises are not evaluated
#' merely to construct an error message.
#' @keywords internal
#' @noRd
.check_unused_dots <- function(...) {
  n <- ...length()
  if (n == 0L) {
    return(invisible(NULL))
  }
  .stop_unused_dots(...names(), seq_len(n))
}

#' Report selected unused dots
#' @keywords internal
#' @noRd
.stop_unused_dots <- function(names, positions) {
  if (length(positions) == 0L) {
    return(invisible(NULL))
  }
  if (is.null(names)) {
    names <- rep("", max(positions))
  }
  labels <- ifelse(
    nzchar(names[positions]),
    paste0("'", names[positions], "'"),
    paste0("..", positions)
  )
  stop(
    sprintf(
      "unused argument%s: %s",
      if (length(positions) == 1L) "" else "s",
      paste(labels, collapse = ", ")
    ),
    call. = FALSE
  )
}

#' Validate and return a whole-number count
#' @keywords internal
#' @noRd
check_count <- function(x, name, minimum = 1L) {
  valid <- is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x == floor(x) &&
    x >= minimum &&
    x <= .Machine$integer.max
  if (!valid) {
    stop(
      sprintf("'%s' must be an integer >= %d", name, minimum),
      call. = FALSE
    )
  }
  as.integer(x)
}

#' Check that a value is a positive scalar
#' @keywords internal
#' @noRd
check_scalar <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop(sprintf("'%s' must be a numeric scalar", name), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("'%s' must not be NA", name), call. = FALSE)
  }
  if (is.infinite(x)) {
    stop(sprintf("'%s' must be finite", name), call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(sprintf("'%s' must be positive", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Aggregate mean of a stratified frame, with a zero that survives rounding
#'
#' A CV divides by this mean, so whether it is zero decides between a
#' number and `Inf`. Stratum means that cancel do not cancel to exactly
#' zero: what is left depends on the order the sum was taken in and on the
#' width of the accumulator, which differs between platforms. A total below
#' the scale of the terms that formed it is that zero, and is returned as
#' an exact one so every caller can keep testing for it directly. The
#' tolerance is the same relative one [varcomp()] uses to decide that an
#' outcome is centered on zero.
#' @keywords internal
#' @noRd
.aggregate_mean <- function(share, mean) {
  ybar <- sum(share * mean)
  if (abs(ybar) <= sqrt(.Machine$double.eps) * sum(share * abs(mean))) {
    return(0)
  }
  ybar
}

#' Check a population mean used as the denominator of a CV
#'
#' A coefficient of variation is a magnitude, `SE / abs(mu)`, and the
#' sample-size formulas use the square of the mean, so the sign of a mean
#' is immaterial to them: a quantity that is negative on average, like a
#' net change, can be planned for exactly as one that is positive. Zero is
#' the value with no relative scale, and that is what is rejected.
#' @keywords internal
#' @noRd
check_mu <- function(mu, name = "mu", zero_ok = FALSE) {
  if (!is.numeric(mu) || length(mu) != 1L) {
    stop(sprintf("'%s' must be a numeric scalar", name), call. = FALSE)
  }
  if (anyNA(mu)) {
    stop(sprintf("'%s' must not be NA", name), call. = FALSE)
  }
  if (!is.finite(mu)) {
    stop(sprintf("'%s' must be finite", name), call. = FALSE)
  }
  if (mu == 0 && !zero_ok) {
    stop(sprintf(
      "'%s' must not be zero: a coefficient of variation has no scale there",
      name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check that a value is in (0, 1)
#' @keywords internal
#' @noRd
check_proportion <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop(sprintf("'%s' must be a numeric scalar", name), call. = FALSE)
  }
  if (anyNA(x) || x <= 0 || x >= 1) {
    stop(sprintf("'%s' must be in (0, 1)", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Default `var_ratio_2` implied by the three-stage variance decomposition
#'
#' The three-stage multiplier is
#' `D = var_ratio_1 icc_1 m q + var_ratio_2 (1 + icc_2 (q - 1))`, where
#' `var_ratio_1 = S^2 / unit_relvar` rescales the components' unit variance to the
#' analysis variable and `var_ratio_2` does the same for the within-PSU part,
#' `S_w^2 = S^2 (1 - icc_1)` when SSU size plays no part in the within-PSU
#' estimator (SSUs of equal size, or drawn proportional to size). Then
#' `var_ratio_2` is not free once `var_ratio_1` and `icc_1` are known:
#' \deqn{k_{ssu} = k_{psu}(1 - \delta_{psu}).}{k_ssu = k_psu(1 - delta_psu).}
#' SRS of SSUs of unequal size gives a larger value, which only `varcomp()` or
#' an explicit pair can supply. On a census-tract frame it reached three
#' times this one.
#' Equivalently `var_ratio_1 icc_1 + var_ratio_2 = var_ratio_1`, which is what makes
#' `D` collapse to `var_ratio_1` when `m = q = 1` and there is no clustering left
#' to inflate anything. Defaulting `var_ratio_2` to 1 asserts that the within-PSU
#' variance is the whole variance, contradicting any positive `icc_1` in
#' the same expression, and inflates the design effect by
#' `var_ratio_1 icc_1 (1 + icc_2 (q - 1))`.
#' @keywords internal
#' @noRd
.var_ratio_ssu_default <- function(var_ratio_psu, icc_psu) {
  value <- var_ratio_psu * (1 - icc_psu)
  if (any(!is.finite(value)) || any(value <= 0)) {
    stop(
      "'icc_1' of 1 leaves no within-PSU variance, so a three-stage design is not identified. Supply 'var_ratio_2' explicitly or plan two stages",
      call. = FALSE
    )
  }
  value
}

#' Expand a three-stage `var_ratio` to the (var_ratio_1, var_ratio_2) pair
#'
#' A scalar `var_ratio` supplies `var_ratio_1` only, and `var_ratio_2`
#' then follows from the decomposition. A length-2 `var_ratio` is taken as
#' supplied and checked for consistency.
#' @keywords internal
#' @noRd
.stage_k_pair <- function(var_ratio, icc) {
  if (length(var_ratio) == 1L) {
    return(c(var_ratio, .var_ratio_ssu_default(var_ratio, icc[1L])))
  }
  rep_len(var_ratio, 2L)
}

#' Check that exactly one of moe/cv/rmoe is specified
#'
#' `rmoe` is the margin of error relative to the estimand, so it competes
#' with both of the others: it fixes the same interval half-width `moe`
#' does, and it states it in the relative units `cv` uses.
#' @keywords internal
#' @noRd
check_precision <- function(moe, cv, rmoe = NULL) {
  has_moe <- !is.null(moe)
  has_cv <- !is.null(cv)
  has_rmoe <- !is.null(rmoe)
  if (has_moe + has_cv + has_rmoe != 1L) {
    stop("specify exactly one of 'moe', 'cv', or 'rmoe'", call. = FALSE)
  }
  if (has_moe) {
    check_scalar(moe, "moe")
  }
  if (has_cv) {
    check_scalar(cv, "cv")
  }
  if (has_rmoe) {
    check_scalar(rmoe, "rmoe")
  }
  invisible(TRUE)
}

#' Turn a relative margin of error into an absolute one
#'
#' The single conversion point for the scalar interfaces. `abs()` on the
#' estimand keeps a negative mean giving a positive margin of error, the
#' same convention `.convert_moe_to_cv()` uses in the other direction.
#' Neither current caller can observe the sign, `p` being positive and the
#' mean formulas reading `moe` squared, so `abs()` is what keeps the
#' helper right for a caller that reads it unsquared, as
#' `.indicators_moe_from_rmoe()` does on the table side.
#' @keywords internal
#' @noRd
.moe_from_rmoe <- function(rmoe, estimand, estimand_name) {
  if (is.null(estimand) || anyNA(estimand)) {
    stop(
      sprintf("'%s' is required when 'rmoe' is specified", estimand_name),
      call. = FALSE
    )
  }
  if (any(estimand == 0)) {
    stop(
      sprintf("'rmoe' is undefined at '%s' = 0", estimand_name),
      call. = FALSE
    )
  }
  rmoe * abs(estimand)
}

#' Express a margin of error relative to the estimand it bounds
#'
#' It answers where `cv` answers, and in the same two edge states: `NA`
#' when the estimand is unknown, which is what an allocation over a frame
#' carrying no `mean` or `p` reports, and `Inf` at an estimand of zero,
#' which no relative quantity is defined against. Both fall out of the
#' division, so neither is special-cased.
#' @keywords internal
#' @noRd
.rmoe_from_moe <- function(moe, estimand, at_zero = Inf) {
  if (is.null(estimand) || length(estimand) == 0L) {
    return(rep(NA_real_, length(moe)))
  }
  out <- moe / abs(estimand)
  # Allocations keep Inf at a zero stratum mean, as their .cv does. The
  # single-estimand families report NA there, a level with no scale.
  out[!is.na(estimand) & estimand == 0] <- at_zero
  out
}

#' Standard error relative to a level, NA where there is no level
#'
#' Absent or zero, the level gives a relative error no scale, so the answer
#' is `NA` rather than `Inf`, the state `.rmoe_from_moe()` returns for the
#' same families.
#' @keywords internal
#' @noRd
.cv_of <- function(se, level) {
  if (is.null(level) || is.na(level) || level == 0) {
    return(NA_real_)
  }
  se / abs(level)
}

#' Refuse an expected take the planning approximation cannot carry
#'
#' Response enters these formulas as a deterministic expected take, `take *
#' rate`. Below one expected respondent that substitution stops describing
#' the design: the within-unit variance term it feeds is built for a take
#' that yields at least one observation. The design itself is perfectly
#' valid, so the message says the approximation is unsupported there rather
#' than that the design is impossible.
#' @keywords internal
#' @noRd
.check_expected_take <- function(take, rate, take_name, rate_name,
                                 where = "") {
  if (!is.finite(take) || !is.finite(rate)) {
    return(invisible(TRUE))
  }
  if (take * rate >= 1 - 1e-9) {
    return(invisible(TRUE))
  }
  stop(
    sprintf(
      "'%s' = %.4g at '%s' = %.4g%s expects %.4g responses per unit, and this package's expected-take planning approximation is unsupported below one. Raise '%s', raise '%s', or plan the stage with a design that does not rely on it",
      take_name, take, rate_name, rate, where, take * rate, take_name,
      rate_name
    ),
    call. = FALSE
  )
}

#' The smallest whole gross take with at least one expected respondent
#'
#' The whole-unit counterpart of `.check_expected_take()`. Rows sharing a
#' take need it for the lowest rate among them. A rate above one, which
#' `eligible_per_unit` allows, gives one.
#' @keywords internal
#' @noRd
.take_floor <- function(rate) {
  max(1L, as.integer(ceiling(.take_bound(rate) - 1e-9)))
}

#' The whole PSU count a fixed `n_psu` fields
#' @keywords internal
#' @noRd
.whole_psu <- function(n_psu) {
  max(1, round(n_psu))
}

#' The whole take a fixed take fields
#'
#' Capped at the whole population cap first, then rounded, and never below
#' the take floor of the response rates.
#' @keywords internal
#' @noRd
.whole_take <- function(fixed, floor, cap = Inf) {
  max(floor, round(min(fixed, cap)))
}

#' The smallest continuous gross take with one expected respondent
#' @keywords internal
#' @noRd
.take_bound <- function(rate) {
  max(1, 1 / min(rate))
}

#' `optimize()` over a take from its supported bound
#'
#' The search runs over `c(1, upper)` first, so a design the bound does not
#' touch is found exactly as before it existed. Only an optimum below the
#' bound is searched again from the bound. `optimize()` never evaluates the
#' ends of its interval, so a cost still falling at the bound would return a
#' point just inside it. The bound itself is scored and kept when it is no
#' worse.
#' @keywords internal
#' @noRd
.optimize_from <- function(f, lower, upper, ...) {
  if (upper > 1) {
    opt <- optimize(f, interval = c(1, upper), ...)
    if (opt$minimum >= lower) {
      return(opt)
    }
  }
  at_lower <- list(minimum = lower, objective = f(lower))
  if (upper <= lower) {
    return(at_lower)
  }
  opt <- optimize(f, interval = c(lower, upper), ...)
  if (is.finite(at_lower$objective) && at_lower$objective <= opt$objective) {
    at_lower
  } else {
    opt
  }
}

#' Check alpha in (0, 1)
#' @keywords internal
#' @noRd
check_alpha <- function(alpha) {
  check_proportion(alpha, "alpha")
}

#' Check deff > 0
#' @keywords internal
#' @noRd
check_deff <- function(deff) {
  check_scalar(deff, "deff")
  invisible(TRUE)
}

#' Resolve a design parameter that may vary by stratum
#'
#' Three sources, in decreasing precedence: a length-H vector argument, a
#' frame column (NA falling back to the scalar), and the scalar argument.
#' This is the `unit_cost` precedent extended to parameters whose scalar
#' default is not `NULL`, so "supplied" cannot be detected by missingness.
#' The `measures` tables of joint mode read `NA` the same way.
#' @keywords internal
#' @noRd
.alloc_resolve_h <- function(arg, frame, name, H, validate) {
  if (length(arg) > 1L) {
    if (length(arg) != H) {
      stop(sprintf("'%s' must have length 1 or nrow(frame)", name),
           call. = FALSE)
    }
    out <- arg
  } else {
    check_scalar(arg, name, positive = FALSE)
    out <- if (name %in% names(frame)) {
      col <- frame[[name]]
      if (!is.numeric(col)) {
        stop(sprintf("'frame$%s' must be numeric", name), call. = FALSE)
      }
      ifelse(is.na(col), arg, col)
    } else {
      rep(arg, H)
    }
  }
  validate(out, name)
  as.numeric(out)
}

#' Refuse a scalar argument that the joint columns override on every row
#'
#' A joint fit reads a measure row's value first, then its stratum's frame
#' column, then the scalar, so the scalar is ignored when every measure row
#' that a target or the objective uses finds a value in one of the two
#' columns.
#' @keywords internal
#' @noRd
.refuse_joint_argument_beside_columns <- function(frame, measures, targets,
                                                  objective, alpha, df, ...) {
  ## A malformed frame is named by the frame check that follows.
  if (!is.data.frame(frame)) {
    return(invisible(NULL))
  }
  .refuse_argument_beside_column(frame, ...)
  vals <- list(...)
  check <- names(vals)[vapply(vals, function(v) {
    length(v) == 1L && !isTRUE(v == 1)
  }, logical(1))]
  check <- intersect(check, names(measures))
  if (length(check) == 0L || !is.data.frame(measures) ||
        nrow(measures) == 0L || is.null(measures[["stratum"]]) ||
        is.null(frame[["stratum"]]) || (is.null(targets) && is.null(objective))) {
    return(invisible(NULL))
  }
  used <- .bethel_required_measures(
    frame, measures, targets, objective, alpha, df
  )$measures
  f_stratum <- as.character(frame[["stratum"]])
  for (nm in check) {
    value <- used[[nm]]
    if (nm %in% names(frame)) {
      fallback <- frame[[nm]][match(as.character(used[["stratum"]]), f_stratum)]
      value[is.na(value)] <- fallback[is.na(value)]
    }
    if (!anyNA(value)) {
      stop(
        sprintf(
          "every row of 'measures' that a target or the objective uses takes '%s' from its own column or the frame column, so the scalar argument would be ignored. Give it as a column or as the argument, not both",
          nm
        ),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

#' Subset a design parameter that may be scalar or per stratum
#' @keywords internal
#' @noRd
.subset_h <- function(x, idx) {
  if (length(x) == 1L) x else x[idx]
}

#' @keywords internal
#' @noRd
.check_deff_h <- function(x, name = "deff") {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
    stop(sprintf("'%s' must contain positive finite values", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.check_resp_rate_h <- function(x, name = "resp_rate") {
  if (!is.numeric(x) || anyNA(x) || any(x <= 0) || any(x > 1)) {
    stop(sprintf("'%s' must contain values in (0, 1]", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check population size N > 1 or Inf
#' @keywords internal
#' @noRd
check_population_size <- function(N) {
  if (!is.numeric(N) || length(N) != 1L || anyNA(N) || N <= 1) {
    stop("'N' must be greater than 1 (or Inf)", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check weights vector: numeric, positive, non-empty
#' @keywords internal
#' @noRd
check_weights <- function(w, name = "x") {
  if (!is.numeric(w) || length(w) == 0L) {
    stop(
      sprintf("'%s' must be a non-empty numeric vector", name),
      call. = FALSE
    )
  }
  if (anyNA(w)) {
    stop(sprintf("'%s' must not contain NA values", name), call. = FALSE)
  }
  if (any(!is.finite(w))) {
    stop(sprintf("'%s' must contain only finite values", name), call. = FALSE)
  }
  if (any(w <= 0)) {
    stop(sprintf("'%s' must contain only positive values", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check stage cost vector
#' @keywords internal
#' @noRd
check_stage_cost <- function(stage_cost) {
  if (!is.numeric(stage_cost) || length(stage_cost) < 2L) {
    stop("'stage_cost' must be a numeric vector of length >= 2", call. = FALSE)
  }
  if (anyNA(stage_cost) || any(stage_cost <= 0)) {
    stop("'stage_cost' must contain positive values only", call. = FALSE)
  }
  if (any(is.infinite(stage_cost))) {
    stop("'stage_cost' must contain finite values only", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check fixed cost (C0)
#' @keywords internal
#' @noRd
check_fixed_cost <- function(fixed_cost, budget = NULL) {
  if (
    !is.numeric(fixed_cost) ||
      length(fixed_cost) != 1L ||
      is.na(fixed_cost) ||
      !is.finite(fixed_cost) ||
      fixed_cost < 0
  ) {
    stop("'fixed_cost' must be a non-negative numeric scalar", call. = FALSE)
  }
  if (!is.null(budget) && fixed_cost >= budget) {
    stop("'fixed_cost' must be less than 'budget'", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check icc vector
#' @keywords internal
#' @noRd
check_icc <- function(icc, expected_length = NULL) {
  if (!is.numeric(icc) || length(icc) == 0L) {
    stop("'icc' must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(icc) || any(icc < 0) || any(icc > 1)) {
    stop("'icc' values must be in [0, 1]", call. = FALSE)
  }
  if (!is.null(expected_length) && length(icc) != expected_length) {
    stop(
      sprintf("'icc' must have length %d (stages - 1)", expected_length),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Require a unit relvariance a design effect could not supply
#'
#' A `varcomp()` estimated from data carries its own `unit_relvar`, and
#' `.check_varcomp_components()` refuses another. One backed out of a design
#' effect has none to give: a scalar design effect fixes the product
#' `var_ratio (1 + icc (b - 1))` and says nothing about the scale of the
#' outcome. Leaving the default in place would price a CV off a
#' relvariance of 1 without ever saying so.
#' @keywords internal
#' @noRd
.check_relvar_identified <- function(icc, supplied, context) {
  if (inherits(icc, "svyplan_varcomp") && identical(icc$source, "deff") &&
      !supplied) {
    stop(
      sprintf(
        "'unit_relvar' is not identified by a design effect, which fixes 'icc' and 'var_ratio' only. Supply 'unit_relvar' directly to %s",
        context
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Refuse a component a `varcomp()` passed as `icc` already sets
#'
#' Runs before the profile merge, the only point at which an explicit
#' argument is distinguishable from a profile default. An explicit value is
#' refused, and a profile value is dropped, since an explicit `icc` overrides
#' the profile. A back-out sets no `unit_relvar`, so that one stays.
#' @keywords internal
#' @noRd
.check_varcomp_components <- function(icc, unit_relvar, var_ratio, plan,
                                      context) {
  vc <- icc %||% if (inherits(plan, "svyplan")) plan$defaults$icc
  if (!inherits(vc, "svyplan_varcomp")) return(plan)
  backout <- identical(vc$source, "deff")
  if (!is.null(var_ratio)) {
    stop(
      sprintf(
        "'var_ratio' is set by the varcomp() result passed as 'icc'. Leave it unset in %s, %s",
        context,
        if (backout) {
          "and give it to varcomp(deff = , var_ratio = ) instead"
        } else {
          "or pass 'icc' as numbers to use another value"
        }
      ),
      call. = FALSE
    )
  }
  if (!backout && !is.null(unit_relvar)) {
    stop(
      sprintf(
        "'unit_relvar' is estimated by the varcomp() result passed as 'icc'. Leave it unset in %s, or pass 'icc' as numbers to use another value",
        context
      ),
      call. = FALSE
    )
  }
  if (inherits(plan, "svyplan")) {
    drop <- c("var_ratio", if (!backout) "unit_relvar")
    plan$defaults <- plan$defaults[!names(plan$defaults) %in% drop]
  }
  plan
}

#' Restrict a profile to the defaults the chosen cluster entry mode accepts
#'
#' A profile describes a design, not the way a design is entered. Its scalar
#' homogeneity is refused by the several-indicators mode, where the frame
#' carries one value per row, and its domain floor means nothing to a single
#' indicator. Its population is read only by the stage correction. Dropping
#' the inapplicable defaults keeps one profile usable from every mode.
#' @keywords internal
#' @noRd
.plan_for_cluster_mode <- function(plan, has_indicators, fpc = "auto",
                                   stages = NULL) {
  if (!inherits(plan, "svyplan") || length(plan$defaults) == 0L) {
    return(plan)
  }
  drop <- if (has_indicators) {
    c("icc", "unit_relvar", "var_ratio")
  } else {
    "min_n_domain"
  }
  # Drop scalar N from a shared 1-stage profile unless fpc = "stage" explicitly;
  # drop any N when fpc = "none".
  if (identical(fpc, "none")) {
    drop <- c(drop, "N")
  } else if (identical(fpc, "auto") && "N" %in% names(plan$defaults)) {
    n_len <- length(plan$defaults[["N"]])
    if (n_len < 2L || (!is.null(stages) && n_len != stages)) {
      drop <- c(drop, "N")
    }
  }
  plan$defaults <- plan$defaults[!names(plan$defaults) %in% drop]
  plan
}

#' Refuse an argument the indicator frame carries per row
#' @keywords internal
#' @noRd
.stop_carried_by_column <- function(given, cols) {
  stop(
    sprintf(
      "%s %s carried by the 'indicators' column%s %s, one value per row",
      .quote_names(given),
      if (length(given) > 1L) "are" else "is",
      if (length(given) > 1L) "s" else "",
      .quote_names(unname(cols[given]))
    ),
    call. = FALSE
  )
}

#' Refuse a table of indicators handed to the numeric first argument
#'
#' The two entry modes take their first argument from different vocabularies,
#' and `indicators` is never positional, so a frame arriving here is a
#' misplaced table rather than a malformed vector.
#' @keywords internal
#' @noRd
.stop_frame_in_first_slot <- function(arg, what) {
  stop(
    sprintf(
      "'%s' is a numeric vector of %s. A table of indicators goes to 'indicators'",
      arg, what
    ),
    call. = FALSE
  )
}

#' Refuse an argument that belongs to the several-indicators mode
#' @keywords internal
#' @noRd
.stop_applies_to_indicators <- function(given) {
  stop(
    sprintf(
      "%s appl%s to 'indicators', a data frame with one row per indicator",
      .quote_names(given),
      if (length(given) > 1L) "y" else "ies"
    ),
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.quote_names <- function(x) paste(sprintf("'%s'", x), collapse = ", ")

#' Tolerance for detecting numerically degenerate cluster homogeneity
#' @keywords internal
#' @noRd
.cluster_icc_tol <- function() {
  1e-8
}

#' Reject cluster homogeneity values that are too close to 0 or 1
#' @keywords internal
#' @noRd
.check_cluster_icc_open <- function(icc, context = "n_cluster()") {
  tol <- .cluster_icc_tol()
  bad <- icc <= tol | icc >= 1 - tol
  if (any(bad)) {
    stop(
      sprintf(
        "'icc' must stay away from 0 and 1 for %s; values <= %.0e or >= %.8f make cluster optimization degenerate",
        context,
        tol,
        1 - tol
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate and preserve or reorder a named stage cost vector
#'
#' Accepts numbered names (`cost_1..cost_S`) in any order, or arbitrary unique
#' user stage names (e.g. `c(district = 500, ea = 100, hh = 50)`), keeping the
#' caller's stage names when provided.
#' @keywords internal
#' @noRd
.reorder_stage_cost <- function(stage_cost) {
  nms <- names(stage_cost)
  if (is.null(nms)) {
    return(unname(stage_cost))
  }
  if (anyNA(nms) || any(!nzchar(nms))) {
    stop("'stage_cost' names must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("duplicate names in 'stage_cost'", call. = FALSE)
  }
  stages <- length(stage_cost)
  numbered <- paste0("cost_", seq_len(stages))
  if (all(nms %in% numbered)) {
    out <- unname(stage_cost[match(numbered, nms)])
    names(out) <- paste0("stage_", seq_len(stages))
    return(out)
  }
  stage_numbered <- paste0("stage_", seq_len(stages))
  if (all(nms %in% stage_numbered)) {
    out <- unname(stage_cost[match(stage_numbered, nms)])
    names(out) <- stage_numbered
    return(out)
  }
  stage_cost
}

#' Validate and name a stage sample-size vector of length S >= 2
#' @keywords internal
#' @noRd
.reorder_n_vec <- function(n) {
  stages <- length(n)
  nms <- names(n)
  default_nms <- paste0("stage_", seq_len(stages))
  if (is.null(nms)) {
    names(n) <- default_nms
    return(n)
  }
  if (anyNA(nms) || any(!nzchar(nms))) {
    stop("'n' names must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("duplicate names in 'n'", call. = FALSE)
  }
  if (all(nms %in% default_nms)) {
    out <- unname(n[match(default_nms, nms)])
    names(out) <- default_nms
    return(out)
  }
  n
}

#' Map stage cost columns in predict(newdata) for arbitrary S >= 2
#' @keywords internal
#' @noRd
.cluster_cost_col_map <- function(cols, stages, stage_names = NULL) {
  canonical <- paste0("cost_", seq_len(stages))
  allowed <- canonical
  dict <- setNames(canonical, canonical)
  used <- intersect(cols, allowed)
  if (length(used) == 0L) {
    return(list(allowed = allowed, map = setNames(character(0), character(0)),
                canonical = canonical))
  }
  mapped <- unname(dict[used])
  list(allowed = allowed, map = setNames(mapped, used), canonical = canonical)
}

#' Map numbered stage parameter columns in predict(newdata) for arbitrary S >= 2
#' @keywords internal
#' @noRd
.cluster_stage_col_map <- function(
  cols,
  name,
  stage_count,
  allow_scalar_alias = FALSE
) {
  canonical <- paste0(name, "_", seq_len(stage_count))
  allowed <- canonical
  dict <- setNames(canonical, canonical)
  if (allow_scalar_alias && stage_count == 1L) {
    allowed <- c(allowed, name)
    dict <- c(dict, setNames(canonical[1L], name))
  }
  used <- intersect(cols, allowed)
  if (length(used) == 0L) {
    return(list(
      allowed = allowed,
      map = setNames(character(0), character(0)),
      canonical = canonical
    ))
  }
  mapped <- unname(dict[used])
  if (anyDuplicated(mapped)) {
    dup_stage <- mapped[duplicated(mapped)][1L]
    offenders <- used[mapped == dup_stage]
    stop(
      sprintf(
        "newdata cannot contain multiple columns for %s: %s",
        dup_stage,
        paste(sQuote(offenders), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  list(allowed = allowed, map = setNames(mapped, used), canonical = canonical)
}

#' Apply mapped cost columns from predict row params to base cost
#' @keywords internal
#' @noRd
.apply_cluster_cost_cols <- function(base_cost, params, cost_col_map) {
  if (length(cost_col_map) == 0L) {
    return(base_cost)
  }
  out <- base_cost
  canonical <- paste0("cost_", seq_along(base_cost))
  for (col in names(cost_col_map)) {
    idx <- match(cost_col_map[[col]], canonical)
    out[idx] <- params[[col]]
  }
  out
}

#' Apply mapped stage columns from predict row params to a base stage vector
#' @keywords internal
#' @noRd
.apply_cluster_stage_cols <- function(
  base_stage,
  params,
  stage_col_map,
  canonical
) {
  out <- unname(base_stage)
  if (length(stage_col_map) == 0L) {
    return(out)
  }
  for (col in names(stage_col_map)) {
    idx <- match(stage_col_map[[col]], canonical)
    out[idx] <- params[[col]]
  }
  out
}

#' Reorder a named stage-indexed vector (numbered `name_1..name_k`)
#' @keywords internal
#' @noRd
.reorder_stage_vec <- function(x, name) {
  nms <- names(x)
  if (is.null(nms)) {
    return(x)
  }
  if (anyNA(nms) || any(!nzchar(nms))) {
    stop(sprintf("'%s' names must be non-empty", name), call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop(sprintf("duplicate names in '%s'", name), call. = FALSE)
  }
  k <- length(x)
  canonical <- paste0(name, "_", seq_len(k))
  unknown <- setdiff(nms, canonical)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "unrecognized names in '%s': %s. Expected %s",
        name,
        paste(sQuote(unknown), collapse = ", "),
        paste(sQuote(canonical), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unname(x[match(canonical, nms)])
}

#' Check overlap fraction in \eqn{[0, 1]}
#' @keywords internal
#' @noRd
check_overlap <- function(overlap) {
  # The lag is named where it is extracted, never defaulted here.
  if (inherits(overlap, "svyplan_overlap")) {
    stop(
      "'overlap' is one lag, and a rotation overlap covers every lag its schedule reaches; name the one this change spans, as in design_overlap(\"4-8-4\")[1] between consecutive occasions or [12] a year apart on monthly ones",
      call. = FALSE
    )
  }
  if (
    !is.numeric(overlap) ||
      length(overlap) != 1L ||
      anyNA(overlap) ||
      overlap < 0 ||
      overlap > 1
  ) {
    stop("'overlap' must be a number in [0, 1]", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check correlation coefficient in \eqn{[0, 1]}
#' @keywords internal
#' @noRd
check_overlap_cor <- function(overlap_cor) {
  if (
    !is.numeric(overlap_cor) || length(overlap_cor) != 1L || anyNA(overlap_cor) || overlap_cor < 0 || overlap_cor > 1
  ) {
    stop("'overlap_cor' must be a number in [0, 1]", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check response rate in (0, 1]
#' @keywords internal
#' @noRd
check_eligible_per_unit <- function(eligible_per_unit, name = "eligible_per_unit") {
  if (
    !is.numeric(eligible_per_unit) ||
      length(eligible_per_unit) != 1L ||
      anyNA(eligible_per_unit) ||
      !is.finite(eligible_per_unit) ||
      eligible_per_unit <= 0
  ) {
    stop(sprintf("'%s' must be a positive finite number", name),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the eligible population a proportion's FPC reads
#'
#' The Bernoulli correction divides by `N * eligible_per_unit - 1`, so a finite
#' eligible population must exceed one unit.
#' @keywords internal
#' @noRd
check_eligible_population <- function(N, eligible_per_unit,
                                      eligibility_rate = 1, label = NULL) {
  N_elig <- N * eligibility_rate * eligible_per_unit
  if (any(is.finite(N_elig) & N_elig <= 1)) {
    stop(
      sprintf(
        "'%s', the eligible population, must be greater than 1",
        label %||% if (all(eligibility_rate == 1)) "N * eligible_per_unit" else
          "N * eligibility_rate * eligible_per_unit"
      ),
      call. = FALSE
    )
  }
  invisible(N_elig)
}

#' Eligible units per issued unit stored in a result's parameters
#' @keywords internal
#' @noRd
.params_yield <- function(params) {
  (params[["eligibility_rate"]] %||% 1) * (params[["eligible_per_unit"]] %||% 1)
}

#' Refuse an eligible-units rate where the allocation cannot carry it
#' @keywords internal
#' @noRd
.check_eligible_scope <- function(eligible_per_unit, frame, measures, psu) {
  given <- !identical(eligible_per_unit, 1) ||
    "eligible_per_unit" %in% c(names(frame), names(measures))
  if (!given) return(invisible(TRUE))
  if (!is.null(psu)) {
    stop(
      "'eligible_per_unit' is not available with a PSU register. Fold a rate at most 1 into 'resp_rate', or drop 'psu' and give 'N_1'",
      call. = FALSE
    )
  }
  if (is.null(measures)) {
    stop(
      "'eligible_per_unit' applies to joint allocation. Give 'measures' and 'targets', one row each for a single indicator",
      call. = FALSE
    )
  }
  check_eligible_per_unit(eligible_per_unit)
  invisible(TRUE)
}

#' Stratum floors from the argument or the frame column of the same name
#'
#' One value applies to every stratum, and one per stratum, as a vector or
#' as a `frame` column, sets each floor apart. The argument overrides the
#' column. `NA` sets no floor for that stratum.
#' @keywords internal
#' @noRd
.resolve_min_n_stratum <- function(min_n_stratum, frame) {
  value <- min_n_stratum %||% if (is.data.frame(frame)) {
    frame[["min_n_stratum"]]
  }
  if (is.null(value)) {
    return(NULL)
  }
  if (length(value) == 1L && is.na(value) && is.data.frame(frame)) {
    value <- rep(NA_real_, nrow(frame))
  }
  if (length(value) == 1L && !is.na(value) && !is.null(min_n_stratum)) {
    check_scalar(value, "min_n_stratum")
    return(value)
  }
  if (!is.numeric(value) && !all(is.na(value))) {
    stop("'min_n_stratum' must be numeric", call. = FALSE)
  }
  if (!is.data.frame(frame) || length(value) != nrow(frame)) {
    stop(
      "'min_n_stratum' must be one number, or one per 'frame' row",
      call. = FALSE
    )
  }
  value <- as.numeric(value)
  if (any(!is.na(value) & (!is.finite(value) | value <= 0))) {
    stop(
      "'min_n_stratum' values must be positive and finite, or NA for no floor",
      call. = FALSE
    )
  }
  value
}

check_resp_rate <- function(resp_rate, name = "resp_rate") {
  if (
    !is.numeric(resp_rate) ||
      length(resp_rate) != 1L ||
      anyNA(resp_rate) ||
      resp_rate <= 0 ||
      resp_rate > 1
  ) {
    stop(sprintf("'%s' must be a number in (0, 1]", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Apply finite population correction
#' @keywords internal
#' @noRd
.apply_fpc <- function(n0, N) {
  if (is.infinite(N)) n0 else n0 / (1 + n0 / N)
}

#' Apply response rate adjustment (inflate n)
#' @keywords internal
#' @noRd
.apply_resp_rate <- function(n, resp_rate) {
  n / resp_rate
}

#' Stage column names before stages were numbered, and their numbered forms
#' @keywords internal
#' @noRd
.unnumbered_stage_columns <- c(
  icc_psu = "icc_1", icc_ssu = "icc_2",
  var_ratio_psu = "var_ratio_1", var_ratio_ssu = "var_ratio_2",
  cost_psu = "cost_1", cost_ssu = "cost_2", cost_tsu = "cost_3",
  N_psu = "N_1", N_ssu = "N_2",
  n_per_psu = "take_2", n_per_ssu = "take_3", n_psu = "n_1",
  resp_rate_psu = "resp_rate_1", resp_rate_ssu = "resp_rate_2"
)

#' Refuse stage columns named by unit rather than numbered by stage
#' @keywords internal
#' @noRd
.refuse_unnumbered_stage_columns <- function(df, arg) {
  if (!is.data.frame(df)) return(invisible(NULL))
  hit <- intersect(names(df), names(.unnumbered_stage_columns))
  if (length(hit) == 0L) return(invisible(NULL))
  stop(
    sprintf(
      "'%s' has %s. Stage columns are numbered by stage: %s",
      arg,
      paste(sQuote(hit), collapse = ", "),
      paste(sprintf("%s for %s", sQuote(.unnumbered_stage_columns[hit]),
                    sQuote(hit)), collapse = ", ")
    ),
    call. = FALSE
  )
}

#' Refuse a scalar argument that a complete frame column overrides
#' @keywords internal
#' @noRd
.refuse_argument_beside_column <- function(frame, ...) {
  vals <- list(...)
  for (nm in names(vals)) {
    v <- vals[[nm]]
    if (nm %in% names(frame) && length(v) == 1L && !isTRUE(v == 1) &&
          !anyNA(frame[[nm]])) {
      stop(
        sprintf(
          "the frame column '%s' has a value on every row, so the scalar argument would be ignored. Give it as the column or as the argument, not both",
          nm
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}
