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
check_mu <- function(mu, name = "mu") {
  if (!is.numeric(mu) || length(mu) != 1L) {
    stop(sprintf("'%s' must be a numeric scalar", name), call. = FALSE)
  }
  if (anyNA(mu)) {
    stop(sprintf("'%s' must not be NA", name), call. = FALSE)
  }
  if (!is.finite(mu)) {
    stop(sprintf("'%s' must be finite", name), call. = FALSE)
  }
  if (mu == 0) {
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

#' Default `var_ratio_ssu` implied by the three-stage variance decomposition
#'
#' The three-stage multiplier is
#' `D = var_ratio_psu icc_psu m q + var_ratio_ssu (1 + icc_ssu (q - 1))`, where
#' `var_ratio_psu = S^2 / unit_relvar` rescales the components' unit variance to the
#' analysis variable and `var_ratio_ssu` does the same for the within-PSU part,
#' `S_w^2 = S^2 (1 - icc_psu)`. So `var_ratio_ssu` is not free once `var_ratio_psu` and
#' `icc_psu` are known:
#' \deqn{k_{ssu} = k_{psu}(1 - \delta_{psu}).}{k_ssu = k_psu(1 - delta_psu).}
#' Equivalently `var_ratio_psu icc_psu + var_ratio_ssu = var_ratio_psu`, which is what makes
#' `D` collapse to `var_ratio_psu` when `m = q = 1` and there is no clustering left
#' to inflate anything. Defaulting `var_ratio_ssu` to 1 asserts that the within-PSU
#' variance is the whole variance, contradicting any positive `icc_psu` in
#' the same expression, and inflates the design effect by
#' `var_ratio_psu icc_psu (1 + icc_ssu (q - 1))`.
#' @keywords internal
#' @noRd
.var_ratio_ssu_default <- function(var_ratio_psu, icc_psu) {
  value <- var_ratio_psu * (1 - icc_psu)
  if (any(!is.finite(value)) || any(value <= 0)) {
    stop(
      "'icc_psu' of 1 leaves no within-PSU variance, so a three-stage design is not identified; supply 'var_ratio_ssu' explicitly or plan two stages",
      call. = FALSE
    )
  }
  value
}

#' Expand a three-stage `var_ratio` to the (var_ratio_psu, var_ratio_ssu) pair
#'
#' A scalar `var_ratio` supplies `var_ratio_psu` only, and `var_ratio_ssu`
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
.rmoe_from_moe <- function(moe, estimand) {
  if (is.null(estimand) || length(estimand) == 0L) {
    return(rep(NA_real_, length(moe)))
  }
  moe / abs(estimand)
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
.check_expected_take <- function(take, rate, take_name, rate_name) {
  if (!is.finite(take) || !is.finite(rate)) {
    return(invisible(TRUE))
  }
  if (take * rate >= 1 - 1e-9) {
    return(invisible(TRUE))
  }
  stop(
    sprintf(
      "'%s' = %.4g at '%s' = %.4g expects %.4g responses per unit, and this package's expected-take planning approximation is unsupported below one. Raise '%s', raise '%s', or plan the stage with a design that does not rely on it",
      take_name, take, rate_name, rate, take * rate, take_name, rate_name
    ),
    call. = FALSE
  )
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

#' Validate a numeric covariate vector against weights
#' @keywords internal
#' @noRd
check_covariate <- function(x, n, name) {
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be numeric", name), call. = FALSE)
  }
  if (length(x) != n) {
    stop(
      sprintf("'%s' must have the same length as weights", name),
      call. = FALSE
    )
  }
  if (anyNA(x)) {
    stop(sprintf("'%s' must not contain NA values", name), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(
      sprintf("'%s' must contain only finite values", name),
      call. = FALSE
    )
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
  if (length(stage_cost) > 3L) {
    stop(
      "4+ stage optimization is not supported; plan the first three stages and fold deeper stages (e.g. persons within households) into the 'deff' passed to n_prop() or n_mean()",
      call. = FALSE
    )
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
#' A `varcomp()` estimated from data carries its own `unit_relvar` and
#' overrides whatever the caller passed. One backed out of a design effect
#' has none to give: a scalar design effect fixes the product
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

#' Restrict a profile to the defaults the chosen cluster entry mode accepts
#'
#' A profile describes a design, not the way a design is entered. Its scalar
#' homogeneity is refused by the several-indicators mode, where the frame
#' carries one value per row, and its domain floor means nothing to a single
#' indicator. Dropping the inapplicable defaults keeps one profile usable
#' from both modes.
#' @keywords internal
#' @noRd
.plan_for_cluster_mode <- function(plan, has_indicators) {
  if (!inherits(plan, "svyplan") || length(plan$defaults) == 0L) {
    return(plan)
  }
  drop <- if (has_indicators) {
    c("icc", "unit_relvar", "var_ratio")
  } else {
    "min_n_domain"
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

#' Refuse the arguments the other cluster entry mode owns
#'
#' Read the mode from `indicators` alone, so a scalar input that the frame
#' already carries per row is an error rather than a value the solver never
#' looks at.
#' @keywords internal
#' @noRd
.check_cluster_mode_args <- function(indicators, icc, unit_relvar, var_ratio,
                                     cv, domains, allocation, domain_sampling,
                                     min_n_domain) {
  if (is.null(indicators)) {
    given <- c(
      if (!is.null(domains)) "domains",
      if (length(allocation) == 1L) "allocation",
      if (length(domain_sampling) == 1L) "domain_sampling",
      if (!is.null(min_n_domain)) "min_n_domain"
    )
    if (length(given) > 0L) .stop_applies_to_indicators(given)
    return(invisible(TRUE))
  }
  cols <- c(
    icc = "icc_psu", unit_relvar = "unit_relvar",
    var_ratio = "var_ratio_psu", cv = "cv"
  )
  supplied <- !vapply(
    list(icc = icc, unit_relvar = unit_relvar, var_ratio = var_ratio, cv = cv),
    is.null,
    logical(1L)
  )
  if (any(supplied)) .stop_carried_by_column(names(cols)[supplied], cols)
  invisible(TRUE)
}

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

#' Asymptotic CV floor for a multistage cluster design at fixed `n_psu`
#'
#' As the free stage sizes grow, CV approaches an irreducible between-PSU
#' variance floor. When both `n_psu` and `n_per_psu` are fixed, the floor
#' is tighter (includes the between-SSU component scaled by 1/n_per_psu).
#' @keywords internal
#' @noRd
.multistage_cv_floor <- function(unit_relvar, var_ratio_psu, icc_psu,
                                 var_ratio_ssu = NULL, icc_ssu = NULL,
                                 rr, n_psu, n_per_psu = NULL) {
  comp <- var_ratio_psu * icc_psu
  if (!is.null(n_per_psu) && !is.null(var_ratio_ssu) && !is.null(icc_ssu)) {
    comp <- comp + var_ratio_ssu * icc_ssu / n_per_psu
  }
  sqrt(unit_relvar * comp / (n_psu * rr))
}

#' Check target CV against the achievable floor and error with a diagnostic
#'
#' @keywords internal
#' @noRd
.check_multistage_feasibility <- function(cv_t, cv_floor, n_psu,
                                          unit_relvar, var_ratio_psu, icc_psu,
                                          var_ratio_ssu = NULL, icc_ssu = NULL,
                                          n_per_psu = NULL,
                                          rr,
                                          labels = NULL,
                                          context = "n_multi()") {
  bad <- cv_floor >= cv_t
  if (!any(bad)) {
    return(invisible(TRUE))
  }

  if (is.null(labels)) {
    labels <- as.character(seq_along(cv_t))
  }
  idx <- which(bad)[1L]
  comp_j <- var_ratio_psu[idx] * icc_psu[idx]
  if (!is.null(n_per_psu) && !is.null(var_ratio_ssu) && !is.null(icc_ssu)) {
    comp_j <- comp_j + var_ratio_ssu[idx] * icc_ssu[idx] / n_per_psu
  }
  required_n_psu <- ceiling(
    unit_relvar[idx] * comp_j / (cv_t[idx]^2 * rr[idx])
  )
  more <- if (sum(bad) > 1L) {
    sprintf(" (%d other indicator(s) also infeasible)", sum(bad) - 1L)
  } else {
    ""
  }
  stop(
    sprintf(
      "%s: target CV %.4g for indicator '%s' is below the achievable floor %.4g at n_psu = %d; increase n_psu to at least %d, or relax target CV above %.4g%s",
      context,
      cv_t[idx],
      labels[idx],
      cv_floor[idx],
      n_psu,
      required_n_psu,
      cv_floor[idx],
      more
    ),
    call. = FALSE
  )
}

#' Reorder a named stage-indexed vector
#'
#'
#' If `x` is unnamed, return as-is. If named, validate and reorder
#' to canonical `c(..._psu, ..._ssu)` order.
#' @param x Numeric vector (length 1 or 2).
#' @param name Parameter name for error messages (e.g., `"icc"`, `"var_ratio"`).
#' @keywords internal
#' @noRd
.reorder_named_vec <- function(x, name, canonical, aliases = character(0)) {
  nms <- names(x)
  if (is.null(nms)) {
    return(unname(x))
  }
  if (anyNA(nms) || any(nms == "")) {
    stop(sprintf("'%s' names must be non-empty", name), call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop(sprintf("duplicate names in '%s'", name), call. = FALSE)
  }

  dict <- setNames(canonical, canonical)
  if (length(aliases) > 0L) {
    dict <- c(dict, aliases)
  }

  unknown <- setdiff(nms, names(dict))
  if (length(unknown) > 0L) {
    expected <- unique(c(canonical, names(aliases)))
    stop(
      sprintf(
        "unrecognized names in '%s': %s; expected %s",
        name,
        paste(sQuote(unknown), collapse = ", "),
        paste(sQuote(expected), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  target <- unname(dict[nms])
  if (anyDuplicated(target)) {
    stop(
      sprintf("named '%s' has overlapping aliases for the same stage", name),
      call. = FALSE
    )
  }

  missing <- setdiff(canonical, target)
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "named '%s' must include all stages: %s",
        name,
        paste(sQuote(canonical), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  idx <- match(canonical, target)
  unname(x[idx])
}

#' Reorder named stage cost vector
#'
#' Supports stage names `cost_psu`, `cost_ssu`, `cost_tsu`.
#' For 2-stage designs, `cost_tsu` is accepted as an alias of `cost_ssu`.
#' @keywords internal
#' @noRd
.reorder_stage_cost <- function(stage_cost) {
  nms <- names(stage_cost)
  if (is.null(nms)) {
    return(unname(stage_cost))
  }

  stages <- length(stage_cost)
  canonical <- if (stages == 2L) {
    c("cost_psu", "cost_ssu")
  } else {
    c("cost_psu", "cost_ssu", "cost_tsu")
  }
  aliases <- if (stages == 2L) {
    c(cost_tsu = "cost_ssu")
  } else {
    character(0)
  }
  .reorder_named_vec(stage_cost, "stage_cost", canonical, aliases)
}

#' Reorder named stage sample-size vector
#'
#' Supports stage names `n_psu`, `n_per_psu`, `n_per_ssu`.
#' @keywords internal
#' @noRd
.reorder_n_vec <- function(n) {
  stages <- length(n)
  out_names <- if (stages == 2L) {
    c("n_psu", "n_per_psu")
  } else {
    c("n_psu", "n_per_psu", "n_per_ssu")
  }

  if (is.null(names(n))) {
    names(n) <- out_names
    return(n)
  }

  n <- .reorder_named_vec(n, "n", out_names)
  names(n) <- out_names
  n
}

#' Reorder named stage cost columns in predict(newdata)
#'
#' @keywords internal
#' @noRd
.cluster_cost_col_map <- function(cols, stages) {
  if (stages == 2L) {
    allowed <- c("cost_psu", "cost_ssu", "cost_tsu")
    dict <- c(
      cost_psu = "cost_psu",
      cost_ssu = "cost_ssu",
      cost_tsu = "cost_ssu"
    )
  } else {
    allowed <- c("cost_psu", "cost_ssu", "cost_tsu")
    dict <- c(
      cost_psu = "cost_psu",
      cost_ssu = "cost_ssu",
      cost_tsu = "cost_tsu"
    )
  }

  used <- intersect(cols, allowed)
  if (length(used) == 0L) {
    return(list(allowed = allowed, map = setNames(character(0), character(0))))
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

  list(allowed = allowed, map = setNames(mapped, used))
}

#' Reorder named stage parameter columns in predict(newdata)
#'
#' Supports stage names like `icc_psu` / `icc_ssu` (or `var_ratio_psu` / `var_ratio_ssu`).
#' Optionally supports a scalar alias (`icc` or `var_ratio`) for 2-stage.
#' @keywords internal
#' @noRd
.cluster_stage_col_map <- function(
  cols,
  name,
  stage_count,
  allow_scalar_alias = FALSE
) {
  canonical <- if (stage_count == 1L) {
    paste0(name, "_psu")
  } else {
    c(paste0(name, "_psu"), paste0(name, "_ssu"))
  }
  allowed <- canonical
  dict <- setNames(canonical, canonical)
  if (allow_scalar_alias) {
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
  out <- unname(base_cost)
  for (col in names(cost_col_map)) {
    idx <- match(cost_col_map[[col]], c("cost_psu", "cost_ssu", "cost_tsu"))
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

#' Reorder a named stage-indexed vector
#'
#' If `x` is unnamed, return as-is. If named, validate and reorder
#' to canonical stage order.
#' @param x Numeric vector (length 1 or 2).
#' @param name Parameter name for error messages (e.g., `"icc"`, `"var_ratio"`).
#' @keywords internal
#' @noRd
.reorder_stage_vec <- function(x, name) {
  nms <- names(x)
  if (is.null(nms)) {
    return(x)
  }
  psu_nm <- paste0(name, "_psu")
  ssu_nm <- paste0(name, "_ssu")
  if (length(x) == 1L) {
    return(.reorder_named_vec(x, name, psu_nm))
  }
  .reorder_named_vec(x, name, c(psu_nm, ssu_nm))
}

#' Check overlap fraction in \[0, 1\]
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

#' Check correlation coefficient in \[0, 1\]
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

#' Normalize scalar or length-2 input to always length-2
#' @keywords internal
#' @noRd
.as_pair <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || anyNA(x) || !all(is.finite(x))) {
    stop(sprintf("'%s' must be finite numeric", name), call. = FALSE)
  }
  if (length(x) == 1L) {
    if (positive && x <= 0) {
      stop(sprintf("'%s' must be positive", name), call. = FALSE)
    }
    return(c(x, x))
  }
  if (length(x) == 2L) {
    if (positive && any(x <= 0)) {
      stop(sprintf("all '%s' elements must be positive", name), call. = FALSE)
    }
    return(x)
  }
  stop(sprintf("'%s' must be length 1 or 2", name), call. = FALSE)
}

#' Validate and normalize n for power functions
#' @keywords internal
#' @noRd
.check_power_n <- function(n) {
  if (!is.numeric(n) || anyNA(n) || !all(is.finite(n))) {
    stop("'n' must be finite numeric", call. = FALSE)
  }
  if (!length(n) %in% c(1L, 2L)) {
    stop(
      "'n' must be length 1 (equal groups) or 2 (unequal groups)",
      call. = FALSE
    )
  }
  if (any(n < 2)) {
    stop("'n' must be >= 2", call. = FALSE)
  }
  n
}

#' Resolve n/ratio when solving for n
#' @keywords internal
#' @noRd
.resolve_ratio <- function(n, ratio) {
  if (is.null(ratio)) {
    return(1)
  }
  if (
    !is.numeric(ratio) ||
      length(ratio) != 1L ||
      is.na(ratio) ||
      ratio <= 0 ||
      !is.finite(ratio)
  ) {
    stop("'ratio' must be a positive finite scalar", call. = FALSE)
  }
  if (!is.null(n) && ratio != 1) {
    stop("'ratio' cannot be used when 'n' is provided", call. = FALSE)
  }
  ratio
}

#' Two-sided interval quantile, normal or t
#'
#' Every confidence interval in the package is a half-width `q * se`, and
#' `q` is the normal quantile when the variance is treated as known and the
#' t quantile on `df` when it is estimated from a design with that many
#' degrees of freedom. Routing every interval site through one function is
#' what keeps a target and the sensitivity of that target to it reading the
#' same quantile: they are the same expression differentiated, and a
#' mismatch between them misreports without failing.
#'
#' `.z_alpha()` below is deliberately separate and stays normal-only. It
#' serves the power functions, where the quantile is a normal deviate for
#' an alternative rather than an interval half-width, and a t-based power
#' calculation is a different procedure than a quantile substitution.
#' @keywords internal
#' @noRd
.q_alpha <- function(alpha, df = NULL) {
  q <- qnorm(1 - alpha / 2)
  if (is.null(df) || all(is.na(df))) {
    return(q)
  }
  df <- rep_len(as.double(df), length(alpha))
  use <- !is.na(df)
  check_df(df[use])
  q[use] <- stats::qt(1 - alpha[use] / 2, df[use])
  q
}

#' Check degrees of freedom for an interval quantile
#'
#' Below one degree of freedom the t quantile is undefined and no interval
#' is identified from the design at all.
#' @keywords internal
#' @noRd
check_df <- function(df, name = "df") {
  # Inf is admitted and is the identity: qt at infinite df is the normal
  # quantile, which is the same statement as treating the variance as known.
  if (!is.numeric(df) || anyNA(df) || any(df < 1)) {
    stop(sprintf("'%s' must be a number >= 1", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Refuse degrees of freedom where the quantile is not an interval
#'
#' A power calculation's quantile is a normal deviate for an alternative,
#' not the half-width of a confidence interval, so substituting a t
#' quantile for it does not produce a t-based power calculation. That is a
#' different procedure. The argument is refused with the reason rather
#' than reported as an unknown name.
#' @keywords internal
#' @noRd
.stop_power_df <- function(...) {
  if ("df" %in% (...names() %||% character(0))) {
    stop(
      "'df' is not accepted by the power functions: their quantile is a normal deviate for an alternative, not an interval half-width, and a t-based power calculation is a different procedure rather than a substituted quantile. It applies to n_prop(), n_mean(), n_alloc() and their precision counterparts",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Compute z_alpha for power functions
#' @keywords internal
#' @noRd
.z_alpha <- function(alpha, alternative) {
  qnorm(1 - alpha / if (alternative == "two.sided") 2 else 1)
}

#' Check power population size N (scalar or length-2)
#' @keywords internal
#' @noRd
.check_power_N <- function(N) {
  if (!is.numeric(N) || anyNA(N) || !all(is.finite(N) | is.infinite(N))) {
    stop("'N' must be numeric, finite, or Inf", call. = FALSE)
  }
  if (!length(N) %in% c(1L, 2L)) {
    stop(
      "'N' must be length 1 (equal groups) or 2 (group-specific)",
      call. = FALSE
    )
  }
  if (any(N <= 1 & !is.infinite(N))) {
    stop("'N' values must be greater than 1 (or Inf)", call. = FALSE)
  }
  if (length(N) == 1L) c(N, N) else N
}

#' Effective N for n2 when ratio constraint links n1 = ratio * n2
#' @keywords internal
#' @noRd
.N2_for_ratio <- function(N_pair, ratio) {
  min(N_pair[2], N_pair[1] / ratio)
}

#' Upper bound for gross n2 from frame and respondent-overlap constraints
#' @keywords internal
#' @noRd
.n2_upper_bound <- function(N_pair, ratio, resp_rate, overlap = 0,
                            within_arm = FALSE) {
  b1 <- if (is.infinite(N_pair[1])) Inf else N_pair[1] / ratio
  b2 <- N_pair[2]
  cap <- min(b1, b2)
  if (overlap > 0) {
    overlap_cap <- if (within_arm) {
      cap / (resp_rate * (2 - overlap))
    } else {
      N_pair[1L] / (resp_rate * (1 + ratio * (1 - overlap)))
    }
    cap <- min(cap, overlap_cap)
  }
  cap
}

#' Solve n2 by inverting a monotone power function
#' @keywords internal
#' @noRd
.solve_n2_from_power <- function(
  target_power,
  power_fn,
  N_pair,
  ratio,
  resp_rate,
  tol = 1e-8,
  overlap = 0,
  within_arm = FALSE
) {
  lo <- max(2, 2 / ratio)
  hi_cap <- .n2_upper_bound(N_pair, ratio, resp_rate, overlap, within_arm)
  if (hi_cap < lo) {
    stop("target power is unattainable under finite population and overlap constraints",
         call. = FALSE)
  }
  p_lo <- suppressWarnings(power_fn(lo))
  if (!is.finite(p_lo)) {
    p_lo <- 0
  }
  if (target_power <= p_lo) {
    return(lo)
  }

  if (is.finite(hi_cap)) {
    hi <- hi_cap
    if (hi <= lo) {
      stop(
        "target power is unattainable under finite population and overlap constraints",
        call. = FALSE
      )
    }
    p_hi <- suppressWarnings(power_fn(hi))
    if (!is.finite(p_hi) || p_hi < target_power - 32 * .Machine$double.eps) {
      stop(
        "target power is unattainable under finite population and overlap constraints",
        call. = FALSE
      )
    }
    if (abs(p_hi - target_power) <= 32 * .Machine$double.eps) return(hi)
  } else {
    hi <- max(2, lo * 2)
    p_hi <- suppressWarnings(power_fn(hi))
    iter <- 0L
    while ((is.na(p_hi) || p_hi < target_power) && iter < 100L) {
      hi <- hi * 2
      p_hi <- suppressWarnings(power_fn(hi))
      iter <- iter + 1L
    }
    if (!is.finite(p_hi) || p_hi < target_power) {
      stop("could not bracket sample size for target power", call. = FALSE)
    }
  }

  uniroot(
    function(n2) power_fn(n2) - target_power,
    interval = c(lo, hi),
    tol = tol
  )$root
}

#' Validate a power target used for inversion
#' @keywords internal
#' @noRd
.check_power_target <- function(power, alpha) {
  if (power <= alpha) {
    stop(
      sprintf(
        "'power' must be greater than 'alpha' (%.4g) when solving for sample size or minimum detectable effect; alpha is the power at zero effect",
        alpha
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Normal-approximation power for a standardized absolute effect
#' @keywords internal
#' @noRd
.normal_power <- function(effect, se, alpha, alternative) {
  if (se == 0) {
    return(if (effect == 0) alpha else 1)
  }
  z_a <- .z_alpha(alpha, alternative)
  z <- abs(effect) / se
  out <- pnorm(z - z_a)
  if (alternative == "two.sided") {
    out <- out + pnorm(-z - z_a)
  }
  min(out, 1)
}

#' Refuse a minimum detectable effect that a census does not have
#'
#' Enumerating the population leaves no sampling variance, so power is
#' `alpha` at exactly zero effect and 1 at every nonzero effect. No effect
#' has any intermediate power, so the requested target has no solution and
#' reporting zero alongside it would state a power the design never achieves.
#' @keywords internal
#' @noRd
.stop_census_mde <- function() {
  stop(
    "no minimum detectable effect exists: the design enumerates its population, so the sampling variance is zero and power is 'alpha' at zero effect and 1 at any nonzero effect",
    call. = FALSE
  )
}

#' Invert normal-approximation power for an absolute effect
#' @keywords internal
#' @noRd
.solve_normal_mde <- function(power, se, alpha, alternative, tol = 1e-10) {
  if (se == 0) {
    .stop_census_mde()
  }
  if (alternative == "one.sided") {
    return((.z_alpha(alpha, alternative) + qnorm(power)) * se)
  }
  hi <- (.z_alpha(alpha, alternative) + qnorm(power)) * se
  hi <- max(hi, se * sqrt(.Machine$double.eps))
  while (.normal_power(hi, se, alpha, alternative) < power) {
    hi <- hi * 2
  }
  uniroot(
    function(effect) .normal_power(effect, se, alpha, alternative) - power,
    c(0, hi),
    tol = tol
  )$root
}

#' Per-group FPC factor (returns 1 for Inf, 0-clamped)
#' @keywords internal
#' @noRd
.fpc_factor <- function(n, N) {
  if (is.infinite(N)) {
    return(1)
  }
  max(0, 1 - n / N)
}

#' Finite population correction on N - 1 degrees of freedom
#'
#' The transformed proportion scales carry their dispersion through
#' \eqn{Var(\hat p) = S^2(1/n - 1/N)} with the Bernoulli population variance
#' \eqn{S^2 = Np(1-p)/(N-1)}, so the factor that multiplies the per-unit term
#' is \eqn{(N-n)/(N-1)}, not \eqn{1 - n/N}. The delta-method variances on the
#' arcsine and log-odds scales inherit it unchanged, because both divide
#' \eqn{Var(\hat p)} by a function of \eqn{p} alone. Matches
#' `.bernoulli_var()` and `.prec_engine_prop()`, and reduces to
#' `.fpc_factor()` as \eqn{N} grows.
#' @keywords internal
#' @noRd
.fpc_factor_prop <- function(n, N) {
  if (is.infinite(N)) {
    return(1)
  }
  max(0, (N - n) / (N - 1))
}

#' Smallest proportion a design measures to a given CV
#'
#' The relative standard error se(p)/p decreases in p, so the CV target has
#' one root over the range that matters and the equality solution is also the
#' threshold: every larger p meets the target. `se` is the sampling standard
#' error under every interval method, so the root is the same under all four
#' and inverts in closed form from se^2 = p q / n_eff. The interval method
#' governs `moe` rather than `se`, so it does not enter here and the
#' signature does not carry it.
#' @keywords internal
#' @noRd
.prec_solve_prop <- function(cv, n, N, deff, resp_rate) {
  n_net <- n * resp_rate
  if (!is.infinite(N) && n_net >= N) {
    stop(
      "a census carries no sampling variance, so every proportion meets any 'cv'; supply 'p' instead",
      call. = FALSE
    )
  }

  # Every positive cv has a root in (0, 1) here, so there is no unattainable
  # case to report: p -> 1 as cv -> 0 and p -> 0 as cv grows.
  1 / (1 + .effective_from_n(n_net, N, deff) * cv^2)
}

#' Smallest proportion a design measures to a given relative margin of error
#'
#' The counterpart of `.prec_solve_prop()` on the interval row of the
#' package's precision 2x2. `moe(p) / p` reads the interval the chosen
#' method builds rather than the sampling variance, so unlike the `cv`
#' solve this one is genuinely method-specific and has no closed form.
#'
#' Two properties the `cv` solve does not share govern the search:
#'
#' - Wilson and Korn-Graubard half-widths do not vanish as `p` approaches
#'   1, so their relative margin of error has a positive floor and a target
#'   below it is unattainable rather than merely extreme.
#' - The back-transformed log-odds half-width turns upward once the logit
#'   spread outgrows `logit(p)`, near `p = 0.999` at `n = 1500` but as low
#'   as `p = 0.97` at `n = 30`. `moe(p) / p` is therefore decreasing then
#'   increasing, and the bracket ends at that turn: the root below it is
#'   the smallest proportion meeting the target, which is what this solve
#'   reports.
#' @keywords internal
#' @noRd
.prec_solve_prop_rmoe <- function(rmoe, n, alpha, N, deff, resp_rate, method,
                                  df = NULL) {
  n_net <- n * resp_rate
  if (!is.infinite(N) && n_net >= N) {
    stop(
      "a census carries no sampling variance, so every proportion meets any 'rmoe'; supply 'p' instead",
      call. = FALSE
    )
  }
  relative <- function(p) {
    .prec_engine_prop(p, n, alpha, N, deff, resp_rate, method, df)$moe / p
  }
  lower <- 1e-12
  cap <- 1 - 1e-9
  upper <- if (method == "logodds") {
    exp(stats::optimize(
      function(l) relative(exp(l)), c(log(lower), log(cap))
    )$minimum)
  } else {
    cap
  }
  attainable <- relative(upper)
  if (rmoe <= attainable) {
    stop(
      sprintf(
        "'rmoe' = %.4g is unattainable: at n = %.4g the '%s' interval's relative margin of error never falls below %.4g",
        rmoe, n, method, attainable
      ),
      call. = FALSE
    )
  }
  if (relative(lower) <= rmoe) {
    return(lower)
  }
  uniroot(
    function(p) relative(p) - rmoe, c(lower, upper),
    tol = .Machine$double.eps^0.75
  )$root
}

#' Shared precision engine for proportions
#'
#' One documented variance equation: n_net = n * resp_rate responding
#' units drive both the leading term and the FPC's sampling fraction, and
#' deff multiplies the SRSWOR variance at n_net, so
#' se^2 = deff * p * q * fpc(n_net) / n_net with fpc(n) = (N - n)/(N - 1).
#' All four methods read that variance through one quantity, the effective
#' size n_eff = n_net / (deff * fpc) at which an infinite-population SRS
#' would reproduce it. A census drives fpc to 0, n_eff to infinity, and
#' every method's margin of error to 0.
#'
#' `se` is that sampling standard error and `cv` is se/p, so both are the
#' same under all four methods: they describe the estimator, not the
#' interval drawn around it. The methods differ only in `moe`, the half-width
#' of the interval they construct at that variance. Deriving `se` as moe/q
#' instead would make the reported standard error move when only the interval
#' construction changed.
#' @keywords internal
#' @noRd
.prec_engine_prop <- function(p, n, alpha, N, deff, resp_rate, method,
                              df = NULL) {
  q <- .q_alpha(alpha, df)
  q_complement <- 1 - p
  n_net <- n * resp_rate
  fpc <- if (is.infinite(N)) 1 else (N - n_net) / (N - 1)
  fpc <- .clamp_fpc(fpc, n_net, N)
  n_eff <- .effective_from_n(n_net, N, deff)
  # A census has no sampling variance under any method. Returning here keeps
  # the single census warning raised by .clamp_fpc() above.
  if (is.infinite(n_eff)) {
    return(list(se = 0, moe = 0, cv = 0))
  }

  se <- sqrt(p * q_complement / n_eff)

  moe <- if (method == "wald") {
    q * se
  } else if (method == "wilson") {
    .wilson_moe(p, n_eff, q)
  } else if (method == "logodds") {
    .logodds_moe(p, n_net, alpha, N, deff, df)
  } else {
    .beta_moe(p, .kg_effective(n_eff, n_net, alpha, df), alpha)
  }

  list(se = se, moe = moe, cv = se / p)
}

#' Convert a net sample size to its infinite-population equivalent
#'
#' `n_eff` is the size at which an infinite-population SRS has the same
#' variance as `n_net` responding units under `deff` and the finite
#' population correction, so `p q / n_eff == deff p q fpc / n_net`. It is
#' infinite for a census, where the design carries no sampling variance.
#' @keywords internal
#' @noRd
.effective_from_n <- function(n_net, N, deff) {
  if (is.infinite(N)) {
    return(n_net / deff)
  }
  if (n_net >= N) {
    return(Inf)
  }
  n_net * (N - 1) / (deff * (N - n_net))
}

#' Invert `.effective_from_n()`
#'
#' Solves `n_eff = n_net (N - 1) / (deff (N - n_net))` for `n_net`. The
#' result approaches but never reaches `N`, so a finite frame can always
#' meet a positive margin of error.
#' @keywords internal
#' @noRd
.n_from_effective <- function(n_eff, N, deff) {
  if (is.infinite(N)) {
    return(n_eff * deff)
  }
  n_eff * deff * N / (N - 1 + n_eff * deff)
}

#' Half-width of the Wilson score interval
#'
#' Centered at `(n p + z^2/2) / (n + z^2)` with half-width
#' `z sqrt(n p q + z^2/4) / (n + z^2)`, evaluated at the effective sample
#' size `n_eff` so that the design effect and the finite population
#' correction enter through the same variance the other methods use. The
#' half-width tends to 1/2 as `n_eff` tends to 0 and to 0 as it grows.
#' @keywords internal
#' @noRd
.wilson_moe <- function(p, n_eff, z) {
  z * sqrt(p * (1 - p) / n_eff + z^2 / (4 * n_eff^2)) / (1 + z^2 / n_eff)
}

#' Degrees-of-freedom adjusted effective sample size (Korn-Graubard 2.2)
#'
#' A variance estimated from `df` degrees of freedom is less stable than one
#' from a simple random sample of `n_net` units, and the Korn-Graubard
#' interval widens for it by scaling the effective size by
#' `(t_{n_net-1}(1 - alpha/2) / t_df(1 - alpha/2))^2`. Leaving `df` `NULL`
#' means "as stable as a simple random sample of the same size" and applies
#' no adjustment.
#' @keywords internal
#' @noRd
.kg_effective <- function(n_eff, n_net, alpha, df = NULL) {
  if (is.null(df) || is.infinite(n_eff)) {
    return(n_eff)
  }
  check_df(df)
  reference <- max(n_net - 1, 1)
  n_eff * (stats::qt(1 - alpha / 2, reference) /
             stats::qt(1 - alpha / 2, df))^2
}

#' Korn-Graubard (1998) confidence limits for a proportion
#'
#' The Clopper-Pearson limits evaluated at the effective sample size, with
#' `x = p * n_eff` playing the role of the observed positive count. Written
#' in the equivalent beta form of Korn and Graubard's equation (1.2).
#' @keywords internal
#' @noRd
.beta_limits <- function(p, n_eff, alpha) {
  if (is.infinite(n_eff)) {
    return(c(p, p))
  }
  x <- p * n_eff
  c(
    stats::qbeta(alpha / 2, x, n_eff - x + 1),
    stats::qbeta(1 - alpha / 2, x + 1, n_eff - x)
  )
}

#' Half-width of the Korn-Graubard interval
#'
#' The interval is asymmetric about `p`, so this half-width is a summary of
#' its length rather than an offset either limit sits at. Use [confint()] for
#' the limits themselves. It decreases monotonically in `n_eff`, which is
#' what lets `n_prop()` invert it.
#' @keywords internal
#' @noRd
.beta_moe <- function(p, n_eff, alpha) {
  limits <- .beta_limits(p, n_eff, alpha)
  (limits[2L] - limits[1L]) / 2
}

#' Half-width of the back-transformed log-odds interval
#'
#' Builds a symmetric interval on the logit scale from the delta-method
#' standard error of `logit(p_hat)`, then maps both endpoints back to the
#' probability scale and halves their distance. The variance of `p_hat`
#' follows the package convention shared with the Wald method,
#' `deff * N / (N - 1) * p q (1/n_net - 1/N)`, so the two agree as the
#' margin of error shrinks.
#' @keywords internal
#' @noRd
.logodds_moe <- function(p, n_net, alpha, N, deff = 1, df = NULL) {
  if (!is.infinite(N) && n_net >= N) {
    warning("net sample size >= population size; moe is 0", call. = FALSE)
    return(0)
  }
  bernoulli <- if (is.infinite(N)) 1 else N / (N - 1)
  fraction <- if (is.infinite(N)) 0 else 1 / N
  var_p <- deff * bernoulli * p * (1 - p) * (1 / n_net - fraction)
  spread <- .q_alpha(alpha, df) * sqrt(var_p) / (p * (1 - p))
  center <- qlogis(p)
  (plogis(center + spread) - plogis(center - spread)) / 2
}

#' Shared precision engine for means
#'
#' Same convention as .prec_engine_prop(), with fpc(n) = 1 - n / N.
#' Smallest mean a design measures to a given CV
#'
#' The standard error of a mean does not involve the mean, so the CV target
#' inverts directly. Only the magnitude is recoverable, `cv` being defined
#' against `abs(mu)`, and the positive root is returned.
#' @keywords internal
#' @noRd
.prec_solve_mean <- function(cv, var, n, alpha, N, deff, resp_rate,
                             df = NULL) {
  se <- .prec_engine_mean(var, NULL, n, alpha, N, deff, resp_rate, df)$se
  if (se == 0) {
    stop(
      "a census carries no sampling variance, so every mean meets any 'cv'; supply 'mu' instead",
      call. = FALSE
    )
  }
  se / cv
}

#' @keywords internal
#' @noRd
.prec_engine_mean <- function(var, mu, n, alpha, N, deff, resp_rate,
                              df = NULL) {
  q <- .q_alpha(alpha, df)
  n_net <- n * resp_rate
  n_eff <- n_net / deff
  fpc <- if (is.infinite(N)) 1 else 1 - n_net / N
  fpc <- .clamp_fpc(fpc, n_net, N)
  se <- sqrt(var * fpc / n_eff)
  list(se = se, moe = q * se,
       cv = if (!is.null(mu)) se / abs(mu) else NA_real_)
}

#' Collision-free key for domain grouping and matching
#'
#' Length-prefixed encoding of each column value, so pasted keys are
#' injective regardless of separators inside the values. Keys built from
#' different data frames are comparable. Missing domain values error.
#' @keywords internal
#' @noRd
.domain_key <- function(df, cols) {
  vals <- lapply(cols, function(col) {
    v <- df[[col]]
    if (anyNA(v)) {
      stop(
        sprintf("domain column '%s' must not contain missing values", col),
        call. = FALSE
      )
    }
    v <- as.character(v)
    paste0(nchar(v), "_", v)
  })
  do.call(paste, c(vals, list(sep = ":")))
}

#' Error when a supplied gross sample exceeds a finite frame
#'
#' Vectorized over groups and indicators. `label` names the offender.
#' @keywords internal
#' @noRd
.check_gross_n <- function(n, N, label = NULL) {
  len <- max(length(n), length(N))
  n <- rep_len(n, len)
  N <- rep_len(N, len)
  bad <- is.finite(N) & n > N * (1 + 1e-9)
  if (any(bad)) {
    i <- which(bad)[1L]
    where <- if (!is.null(label)) {
      sprintf(" for %s", rep_len(label, len)[i])
    } else {
      ""
    }
    stop(
      sprintf(
        "'n' exceeds the population%s: cannot draw %s units without replacement from N = %s",
        where, round(n[i], 1), N[i]
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Error when a gross sample cannot be drawn from a finite frame
#' @keywords internal
#' @noRd
.check_attainable <- function(n, N, resp_rate) {
  if (!is.infinite(N) && n > N * (1 + 1e-9)) {
    stop(
      sprintf(
        "required sample (%s units drawn) exceeds the population (N = %s): even a census yields about %s respondents at resp_rate = %s; the target is unattainable",
        round(n, 1), N, round(N * resp_rate), resp_rate
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Expected number of positive cases a proportion design will yield
#'
#' A count, not a variance: `deff` has no part in it, since a design effect
#' says how precisely the proportion is estimated and nothing about how
#' many cases turn up. The sample is gross, so the response rate nets it
#' down first. `NULL` for anything but a proportion, which is what makes
#' the field conditional on the constructors.
#' @keywords internal
#' @noRd
.expected_cases <- function(n, params, type) {
  if (!identical(type, "proportion") || is.null(params$p) || is.null(n)) {
    return(NULL)
  }
  n * (params$resp_rate %||% 1) * params$p
}

#' Sample size the minimum expected count of positive cases demands
#'
#' Gross units, on the same footing as every other size the package
#' returns: the count is realized among respondents, so the drawn sample
#' carries the `1 / resp_rate` inflation.
#' @keywords internal
#' @noRd
.n_from_min_cases <- function(min_cases, p, resp_rate) {
  check_scalar(min_cases, "min_cases")
  min_cases / (p * resp_rate)
}

#' Smallest issued size whose respondent count clears a target with
#' probability at least `level`
#'
#' Planning at the expected respondent count leaves about half the designs
#' short, so an assurance level converts a required number of respondents
#' into an issued number. Shared by `n_twophase()` and `n_panel()`, whose
#' assurance vocabulary is one thing.
#'
#' The tail is non-decreasing in the issued size, so the answer is bracketed
#' and bisected rather than walked. **Walking up from the expected count
#' returns a size that is not the smallest**, in two directions: below a level
#' of a half the answer lies under the expected count, and at a high response
#' rate the expected count itself overshoots, `need = 2` at a rate of 0.99
#' clearing 0.8 assurance with 2 issued where `ceiling(2 / 0.99)` is 3. The
#' bracket is anchored at `need - 1`, which clears no level at all, fewer
#' issued than needed carrying no chance of reaching it. The expected count is
#' the upper end of the first bracket, which is a starting point rather than a
#' contract: correctness rests on that anchor and on the monotonicity, so any
#' other start converges on the same answer and only the count of
#' `pbinom()` calls changes.
#' @keywords internal
#' @noRd
.assure_size <- function(m, r, level) {
  r <- rep_len(r, length(m))
  vapply(seq_along(m), function(i) {
    need <- ceiling(m[i])
    if (need <= 0) {
      return(0)
    }
    clears <- function(g) {
      stats::pbinom(need - 1L, g, r[i], lower.tail = FALSE) >= level
    }
    lo <- need - 1L
    # A rate is at most 1, so this is at least `need` and no separate floor is
    # needed.
    hi <- ceiling(need / r[i])
    step <- max(1, hi - lo)
    while (!clears(hi)) {
      lo <- hi
      hi <- hi + step
      step <- step * 2
      if (hi > 1e9) {
        stop("the assurance search did not converge", call. = FALSE)
      }
    }
    while (hi - lo > 1) {
      mid <- lo + (hi - lo) %/% 2
      if (clears(mid)) hi <- mid else lo <- mid
    }
    as.double(hi)
  }, numeric(1))
}

#' Validate computed variance term before sqrt
#' @keywords internal
#' @noRd
.safe_variance <- function(V, what = "variance") {
  if (!is.finite(V)) {
    stop(
      sprintf("%s is not finite; check input parameters", what),
      call. = FALSE
    )
  }
  tol <- sqrt(.Machine$double.eps)
  if (V < -tol) {
    stop(
      sprintf("%s is negative; reduce overlap or overlap_cor, or adjust inputs", what),
      call. = FALSE
    )
  }
  if (V < 0) 0 else V
}

#' Validate overlap against ratio or explicit n
#' @keywords internal
#' @noRd
.check_overlap_n <- function(overlap, n = NULL, ratio = NULL) {
  if (overlap == 0) {
    return(invisible(NULL))
  }
  if (!is.null(ratio) && ratio > 1 && overlap > 1 / ratio) {
    stop(
      sprintf("overlap (%.3g) must be <= 1/ratio (%.3g)", overlap, 1 / ratio),
      call. = FALSE
    )
  }
  if (!is.null(n) && length(n) == 2L && overlap > n[2] / n[1]) {
    stop(
      sprintf(
        "overlap (%.3g) must be <= n[2]/n[1] (%.3g)",
        overlap,
        n[2] / n[1]
      ),
      call. = FALSE
    )
  }
}

#' Positive-overlap samples must fit their shared finite frame
#'
#' Sizes here count respondents. Zero overlap retains the public independent-
#' samples convention. Vector inputs check several pairs (DiD arms or pooled
#' lags), not the union of an entire rotation schedule.
#' @keywords internal
#' @noRd
.check_overlap_frame <- function(n1_net, n2_net, N, overlap) {
  distinct <- n2_net + (1 - overlap) * n1_net
  bad <- which(overlap > 0 & is.finite(N) &
                 distinct > N + 32 * .Machine$double.eps * pmax(1, N))
  if (length(bad)) {
    j <- bad[1L]
    # Recycle scalar arguments explicitly for an informative vector diagnostic.
    counts <- distinct[j]
    population <- rep_len(N, length(distinct))[j]
    stop(sprintf(
      "finite population overlap is infeasible: the two responding samples require %.10g distinct units but 'N' is %.10g; increase 'overlap', reduce 'n', or relax the sizing target",
      counts, population
    ), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate and normalize a take_all column
#'
#' Documented as logical or 0/1, so any other number is a mistake rather
#' than a truthy value. `as.logical()` on its own would take 2 for `TRUE`
#' and pin a stratum the caller never meant to pin.
#' @keywords internal
#' @noRd
.check_take_all <- function(take_all, n) {
  if (is.null(take_all)) {
    return(rep(FALSE, n))
  }
  if (is.numeric(take_all)) {
    if (anyNA(take_all) || any(!take_all %in% c(0, 1))) {
      stop("'take_all' must be logical (or 0/1)", call. = FALSE)
    }
    return(take_all != 0)
  }
  if (!is.logical(take_all) || anyNA(take_all)) {
    stop("'take_all' must be logical (or 0/1)", call. = FALSE)
  }
  take_all
}

#' Overlapping samples must share one population
#'
#' Two occasions can only have units in common if they are drawn from the
#' same frame, so a single `N` is the only coherent input once `overlap` is
#' positive. Without overlap the two groups are independent and may
#' legitimately be different populations of different sizes.
#' @keywords internal
#' @noRd
.check_overlap_N <- function(overlap, N_pair) {
  if (overlap > 0 && !isTRUE(all.equal(N_pair[1], N_pair[2]))) {
    stop(
      paste("overlapping samples are drawn from one population, so 'N' must",
            "be a single size or two equal ones"),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Design variance of a difference between two overlapping occasion means
#'
#' The marginal terms carry their own finite population correction, whereas
#' the overlap covariance does not. For two SRSWOR samples drawn from one
#' population of size \eqn{N} and sharing \eqn{k = overlap \cdot n_1}{k = overlap * n_1} units,
#' \deqn{Cov(\bar y_1, \bar y_2) = \rho S_1S_2\{k/(n_1n_2) - 1/N\},}{Cov(ybar_1, ybar_2) = rho S_1S_2\{k/(n_1n_2) - 1/N\},}
#' so the population term enters once as \eqn{1/N} rather than through
#' sample 2's marginal factor. Collecting terms,
#' \deqn{V = \frac{v_1}{n_1} + \frac{v_2}{n_2}
#'       - \frac{2\rho\,overlap\sqrt{v_1v_2}}{n_2}
#'       - \frac{v_1 + v_2 - 2\rho\sqrt{v_1v_2}}{N},}{V = v_1/n_1 + v_2/n_2 - (2 rho overlap sqrt(v_1v_2))/n_2 - (v_1 + v_2 - 2 rho sqrt(v_1v_2))/N,}
#' a per-unit part less a census part. At \eqn{\rho = 1} with equal sizes
#' and variances the census part is zero, the population terms cancel
#' exactly, and the whole reduces to \eqn{2S^2(1 - overlap)/n}.
#' @keywords internal
#' @noRd
.diff_var_fpc <- function(n_eff, var_pair, N_pair, deff, overlap, overlap_cor) {
  if (length(n_eff) == 1L) n_eff <- c(n_eff, n_eff)
  .check_overlap_frame(n_eff[1L], n_eff[2L], N_pair[1L], overlap)
  fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
  fpc2 <- .fpc_factor(n_eff[2], N_pair[2])
  V <- var_pair[1] * fpc1 / n_eff[1] + var_pair[2] * fpc2 / n_eff[2]
  if (overlap > 0) {
    cross <- overlap_cor * sqrt(var_pair[1] * var_pair[2])
    V <- V - 2 * overlap * cross / n_eff[2] +
      if (is.infinite(N_pair[1])) 0 else 2 * cross / N_pair[1]
  }
  deff * V
}

#' Finite-population variance of a Bernoulli variable
#'
#' The variance consumed by the package's without-replacement kernels is the
#' population variance on `N - 1` degrees of freedom. For a population
#' proportion `p`, that is `N * p * (1 - p) / (N - 1)`. The infinite-
#' population limit is the usual `p * (1 - p)`.
#' @keywords internal
#' @noRd
.bernoulli_var <- function(p, N) {
  adj <- ifelse(is.infinite(N), 1, N / (N - 1))
  p * (1 - p) * adj
}

#' Resolve the two occasions' variances and the change they bound
#'
#' A change is planned on one of two scales. On the mean scale the
#' dispersion is supplied directly and the expected change is optional,
#' needed only by the relative measures. On the proportion scale the two
#' occasion proportions supply both: the variances are \eqn{p(1-p)} and the
#' change is \eqn{p_2 - p_1}, so supplying `change` as well would let the
#' two disagree and is refused rather than silently preferred.
#'
#' The proportion scale carries the same \eqn{N/(N-1)} adjustment
#' `.prec_engine_prop()` applies, because the variance the change formula
#' wants is the population variance \eqn{S^2} on \eqn{N-1} degrees of
#' freedom and a Bernoulli population's is \eqn{Np(1-p)/(N-1)}, not
#' \eqn{p(1-p)}. Without it a finite `N` would make one occasion of
#' `prec_change(p = )` disagree with `prec_prop()` on the same design. The
#' mean scale needs no adjustment: `var` is already defined on \eqn{N-1}.
#' @keywords internal
#' @noRd
.change_inputs <- function(var, sd, p, change, N_pair) {
  has_var <- !is.null(var) || !is.null(sd)
  has_p <- !is.null(p)
  if (has_var + has_p != 1L) {
    stop("specify exactly one of 'var' (or 'sd') or 'p'", call. = FALSE)
  }
  if (has_p) {
    if (!is.null(change)) {
      stop(
        "'change' is p[2] - p[1] on the proportion scale; do not supply it",
        call. = FALSE
      )
    }
    if (!is.numeric(p) || length(p) != 2L) {
      stop("'p' must be two proportions, one per occasion", call. = FALSE)
    }
    check_proportion(p[1L], "p[1]")
    check_proportion(p[2L], "p[2]")
    return(list(
      var_pair = .bernoulli_var(p, N_pair),
      change = p[2L] - p[1L],
      p = p
    ))
  }
  var <- .resolve_var(var, sd)
  if (!is.null(change)) check_scalar(change, "change", positive = FALSE)
  list(var_pair = .as_pair(var, "var"), change = change, p = NULL)
}

#' Refuse a between-occasion correlation two proportions cannot have
#'
#' Two Bernoulli variables with means \eqn{p_1} and \eqn{p_2} admit
#' correlations only up to the Frechet-Hoeffding bound
#' \eqn{(\min(p_1,p_2) - p_1p_2)/\sqrt{p_1q_1p_2q_2}}{(min(p_1,p_2) - p_1p_2)/sqrt(p_1q_1p_2q_2)}. Above it no joint
#' distribution exists, so the covariance the change formula subtracts
#' describes nothing. Checked only where the marginals are known, which is
#' the proportion scale, and only where the correlation is used at all.
#' @keywords internal
#' @noRd
.check_bernoulli_cor <- function(p, overlap_cor, overlap) {
  if (is.null(p) || overlap == 0 || overlap_cor == 0) {
    return(invisible(TRUE))
  }
  bound <- (min(p) - p[1L] * p[2L]) /
    sqrt(p[1L] * (1 - p[1L]) * p[2L] * (1 - p[2L]))
  if (overlap_cor > bound * (1 + 1e-9)) {
    stop(
      sprintf(
        "'overlap_cor' (%.4g) exceeds the largest correlation two proportions of %.4g and %.4g can have (%.4g)",
        overlap_cor, p[1L], p[2L], bound
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Feasible second marginal at a requested Bernoulli correlation
#' @keywords internal
#' @noRd
.bernoulli_cor_interval <- function(p1, overlap_cor, overlap) {
  if (overlap == 0 || overlap_cor == 0) {
    return(c(0, 1))
  }
  rho2 <- overlap_cor^2
  q1 <- 1 - p1
  c(
    rho2 * p1 / (q1 + rho2 * p1),
    p1 / (p1 + rho2 * q1)
  )
}

#' Check a size supplied per occasion
#'
#' One size means both occasions are the same size, which is the ordinary
#' repeated survey. Two mean the second occasion was resized, which is what
#' makes `overlap` directional.
#' @keywords internal
#' @noRd
.check_change_n <- function(n) {
  if (!is.numeric(n) || anyNA(n) || !all(is.finite(n))) {
    stop("'n' must be finite numeric", call. = FALSE)
  }
  if (!length(n) %in% c(1L, 2L)) {
    stop(
      "'n' must be length 1 (equal occasions) or 2 (one size per occasion)",
      call. = FALSE
    )
  }
  if (any(n <= 0)) {
    stop("'n' must be positive", call. = FALSE)
  }
  n
}

#' Sampling variance of a change between two occasions of one population
#'
#' A thin wrapper over `.diff_var_fpc()` that nets the gross sizes down and
#' pairs a scalar `n`, so the size the user states is the size every entry
#' point in the package states: units drawn, not interviews completed.
#' @keywords internal
#' @noRd
.change_var <- function(n, var_pair, N_pair, deff, resp_rate, overlap,
                        overlap_cor) {
  n_vec <- if (length(n) == 1L) c(n, n) else n
  V <- .diff_var_fpc(
    n_vec * resp_rate, var_pair, N_pair, deff, overlap, overlap_cor
  )
  .safe_variance(V, "change variance")
}

#' Precision of a change at a given size
#' @keywords internal
#' @noRd
.prec_engine_change <- function(n, var_pair, change, N_pair, deff, resp_rate,
                                overlap, overlap_cor, alpha, df = NULL) {
  se <- sqrt(
    .change_var(n, var_pair, N_pair, deff, resp_rate, overlap, overlap_cor)
  )
  list(
    se = se,
    moe = .q_alpha(alpha, df) * se,
    cv = if (!is.null(change)) se / abs(change) else NA_real_
  )
}

#' Second-occasion size that reaches a target standard error on the change
#'
#' The change variance is \eqn{A/n_2 + B} in the net second-occasion size,
#' with the finite population terms collected into a constant \eqn{B} that
#' no sample size can move. Both coefficients are exact on \eqn{n \le N},
#' where the correction factor is still \eqn{1 - n/N}, so the inversion is a
#' division rather than a search. \eqn{B \le 0} always, since the population
#' term it carries is \eqn{-(v_1 + v_2 - 2\rho\sqrt{v_1v_2})/N}{-(v_1 + v_2 - 2 rho sqrt(v_1v_2))/N} at the common
#' \eqn{N} that a positive overlap requires, so the divisor cannot vanish.
#'
#' \eqn{A} can vanish. At equal sizes, equal variances, full overlap and unit
#' correlation the two occasions share every unit and the change is measured
#' without sampling error at any size, so no size is identified and the
#' caller is told which input to relax rather than handed a zero.
#' @keywords internal
#' @noRd
.n_change_from_se <- function(target_se, var_pair, N_pair, deff, resp_rate,
                              ratio, overlap, overlap_cor) {
  cross <- if (overlap > 0) overlap_cor * sqrt(var_pair[1] * var_pair[2]) else 0
  A <- var_pair[1] / ratio + var_pair[2] - 2 * overlap * cross
  B <- -var_pair[1] / N_pair[1] - var_pair[2] / N_pair[2] +
    if (is.infinite(N_pair[1])) 0 else 2 * cross / N_pair[1]
  if (A <= sqrt(.Machine$double.eps) * sum(var_pair)) {
    stop(
      "the change carries no sampling variance at any size (overlap and overlap_cor leave nothing to sample); reduce either, or size the occasions unequally",
      call. = FALSE
    )
  }
  n2_net <- A / (target_se^2 / deff - B)
  n2 <- n2_net / resp_rate
  if (ratio == 1) n2 else c(ratio * n2, n2)
}

#' Number of occasions entering a pooled estimate
#' @keywords internal
#' @noRd
.MAX_OCCASIONS <- 1000L

.check_occasions <- function(occasions) {
  if (
    !is.numeric(occasions) || length(occasions) != 1L || anyNA(occasions) ||
      !is.finite(occasions) || occasions < 2 ||
      occasions != round(occasions)
  ) {
    stop("'occasions' must be a whole number of at least 2", call. = FALSE)
  }
  # Dense assembly is cubic, so a mistyped figure fails rather than allocates.
  if (occasions > .MAX_OCCASIONS) {
    stop(
      sprintf(
        "'occasions' must be at most %d; the covariance is assembled densely, so a larger horizon is a typing error more often than a design",
        .MAX_OCCASIONS
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve the overlap at every lag a pooled estimate spans
#'
#' A bare number is respondent overlap, the reading every other entry point
#' in the package uses. A `svyplan_overlap` is issued overlap, which is a
#' different quantity below full response, so it is accepted only where the
#' two coincide and refused with the conversion named otherwise. Beyond a
#' schedule's life nothing is shared, which is why an object shorter than
#' the horizon is padded with zeros rather than rejected.
#' @keywords internal
#' @noRd
.resolve_lag_overlap <- function(overlap, occasions, resp_rate,
                                 basis = NULL) {
  max_lag <- occasions - 1L
  if (inherits(overlap, "svyplan_overlap")) {
    .check_overlap_basis("issued", resp_rate)
    ov <- as.numeric(overlap)
    return(list(
      overlap = c(ov, rep(0, max(0L, max_lag - length(ov))))[seq_len(max_lag)],
      basis = "issued"
    ))
  }
  if (
    !is.numeric(overlap) || anyNA(overlap) || any(overlap < 0) ||
      any(overlap > 1)
  ) {
    stop("'overlap' must be numbers in [0, 1]", call. = FALSE)
  }
  # A stored profile carries the basis it was resolved under, so a round
  # trip that lowers the response rate meets the same refusal the first
  # call would have.
  basis <- basis %||% "respondent"
  .check_overlap_basis(basis, resp_rate)
  if (length(overlap) == 1L) {
    return(list(overlap = rep(overlap, max_lag), basis = basis))
  }
  if (length(overlap) != max_lag) {
    stop(
      sprintf(
        "'overlap' must be one number or one per lag (%d for %d occasions)",
        max_lag, occasions
      ),
      call. = FALSE
    )
  }
  list(overlap = overlap, basis = basis)
}

#' An issued overlap is only a respondent overlap at full response
#'
#' The refusal names the conversion rather than performing it, since which
#' assumption to make about response persisting across occasions is the
#' planner's to state. Checked wherever a profile is resolved, so a stored
#' one meets it again on a round trip or in a grid.
#' @keywords internal
#' @noRd
.check_overlap_basis <- function(basis, resp_rate) {
  if (identical(basis, "issued") && !isTRUE(all.equal(resp_rate, 1))) {
    stop(
      "design_overlap() reports the overlap between issued samples, and this variance is formed on respondents; the two coincide only at resp_rate = 1. Supply the respondent overlap you expect, which under independent response is resp_rate times the issued figure",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve the correlation at every lag
#'
#' A panel's correlation falls away with distance, so a single number is the
#' wrong shape for anything but a first look and `cor_decay` is the form most
#' panels have. The two are alternatives rather than a value and a modifier,
#' so supplying both is refused.
#' @keywords internal
#' @noRd
.resolve_lag_cor <- function(overlap_cor, cor_decay, max_lag) {
  if (!is.null(cor_decay)) {
    if (!is.null(overlap_cor)) {
      stop("supply exactly one of 'overlap_cor' or 'cor_decay'", call. = FALSE)
    }
    if (
      !is.numeric(cor_decay) || length(cor_decay) != 1L || anyNA(cor_decay) ||
        cor_decay < 0 || cor_decay > 1
    ) {
      stop("'cor_decay' must be a number in [0, 1]", call. = FALSE)
    }
    return(cor_decay^seq_len(max_lag))
  }
  cor <- overlap_cor %||% 0
  if (!is.numeric(cor) || anyNA(cor) || any(cor < 0) || any(cor > 1)) {
    stop("'overlap_cor' must be numbers in [0, 1]", call. = FALSE)
  }
  if (length(cor) == 1L) {
    return(rep(cor, max_lag))
  }
  if (length(cor) != max_lag) {
    stop(
      sprintf(
        "'overlap_cor' must be one number or one per lag (%d for %d occasions)",
        max_lag, max_lag + 1L
      ),
      call. = FALSE
    )
  }
  cor
}

#' Covariance between two occasion means at each lag
#'
#' The kernel `.diff_var_fpc()` already commits the package to, read at a
#' general lag. A lag whose overlap is zero contributes exactly zero,
#' population term included, which is the same piecewise rule
#' `prec_change()` applies at `overlap = 0` and the only one under which a
#' fresh sample each round pools to `V_level / occasions` exactly.
#' @keywords internal
#' @noRd
.pooled_lag_cov <- function(var, n_net, N, ov, rho) {
  pop <- if (is.infinite(N)) 0 else 1 / N
  out <- rho * var * (ov / n_net - pop)
  out[ov <= 0] <- 0
  out
}

#' The occasion-by-occasion covariance a rotation implies
#' @keywords internal
#' @noRd
.pooled_kernel <- function(var, n_net, N, occasions, ov, rho) {
  c0 <- var * .fpc_factor(n_net, N) / n_net
  lag <- abs(outer(seq_len(occasions), seq_len(occasions), "-"))
  matrix(
    c(c0, .pooled_lag_cov(var, n_net, N, ov, rho))[lag + 1L],
    occasions, occasions
  )
}

#' The assembled covariance must be a covariance
#'
#' Checked on the matrix and not on the quantity reported from it. A kernel
#' that is not positive semidefinite can still give a positive pooled
#' variance, so `.safe_variance()` on the result passes cases this rejects.
#' The failure is ordinary rather than exotic: a lag's covariance changes
#' sign once its overlap falls below the sampling fraction, and a Toeplitz
#' kernel whose signs vary across lags need not be positive semidefinite.
#' An AR(1) correlation does not escape it, since the assembled matrix is a
#' difference of two positive semidefinite parts and not the Schur product
#' of the overlap and correlation kernels alone.
#' @keywords internal
#' @noRd
.check_lag_psd <- function(K, n_gross, n_net, N, ov, rho, occasions) {
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  diag_term <- K[1L, 1L]
  # Relative to the spectrum: an absolute floor would tie the verdict to the
  # outcome's units. The zero matrix has scale zero and passes on equality.
  scale <- max(abs(ev))
  if (min(ev) >= -sqrt(.Machine$double.eps) * scale) {
    return(invisible(NULL))
  }
  cov_m <- .pooled_lag_cov(1, n_net, N, ov, rho)
  neg <- which(cov_m < 0)
  msg <- sprintf(
    "the assembled between-occasion covariance is not a valid covariance (minimum eigenvalue %.3g against a diagonal of %.3g)",
    min(ev), diag_term
  )
  if (length(neg) && !is.infinite(N)) {
    msg <- paste0(
      msg,
      sprintf(
        "; at a sampling fraction of %.3g the covariance at lag %d is negative, which happens once a lag's overlap falls below n/N",
        n_net / N, neg[1L]
      )
    )
  }
  # Interviews per cohort, not calendar span. Incomplete inside the horizon,
  # where the diagnostic is dropped.
  complete <- length(ov) > 0L && ov[length(ov)] == 0
  interviews <- 1 + 2 * sum(ov)
  if (complete && !is.infinite(N) && interviews > 0 &&
        occasions * n_gross / interviews > N) {
    msg <- paste0(
      msg,
      sprintf(
        ". Over %d occasions a rotation interviewing each cohort %.3g times at this size draws %.2f times the population, so the schedule could not be fielded either",
        occasions, interviews, occasions * n_gross / (interviews * N)
      )
    )
  }
  stop(
    paste0(
      msg,
      ". Reduce the sampling fraction, shorten the horizon, or supply an overlap and correlation that hold together across lags"
    ),
    call. = FALSE
  )
}

#' Variance of the equal-weight mean of the occasion estimates
#' @keywords internal
#' @noRd
.pooled_var <- function(n, var, N, deff, resp_rate, occasions, ov, rho) {
  n_net <- n * resp_rate
  K <- .pooled_kernel(var, n_net, N, occasions, ov, rho)
  .check_lag_psd(K, n, n_net, N, ov, rho, occasions)
  .check_overlap_frame(n_net, n_net, N, ov)
  .safe_variance(deff * sum(K) / occasions^2, "pooled variance")
}

#' Precision of a pooled estimate at a given size
#' @keywords internal
#' @noRd
.prec_engine_pooled <- function(n, var, mu, N, deff, resp_rate, occasions,
                                ov, rho, alpha, df = NULL) {
  se <- sqrt(.pooled_var(n, var, N, deff, resp_rate, occasions, ov, rho))
  list(
    se = se,
    moe = .q_alpha(alpha, df) * se,
    cv = if (!is.null(mu)) se / abs(mu) else NA_real_
  )
}

#' Gross size per occasion that reaches a target standard error
#'
#' The pooled variance is \eqn{deff\{A/(n r) + B\}} in the gross size, with
#' \eqn{r} the response rate and the finite population terms collected into
#' \eqn{B}. \eqn{B \le 0} always, since its bracket is at least
#' `occasions` under the nonnegative correlation contract, so the divisor
#' cannot vanish and no target is unattainable for want of precision. The
#' size can still exceed \eqn{N}, which is the ordinary boundary and is
#' checked where every other size is.
#' @keywords internal
#' @noRd
.n_pooled_from_se <- function(target_se, var, N, deff, resp_rate, occasions,
                              ov, rho) {
  m <- seq_len(occasions - 1L)
  w <- 2 * (occasions - m)
  pos <- ov > 0
  A <- var * (occasions + sum(w[pos] * rho[pos] * ov[pos])) / occasions^2
  B <- if (is.infinite(N)) {
    0
  } else {
    -var * (occasions + sum(w[pos] * rho[pos])) / (N * occasions^2)
  }
  deff * A / (resp_rate * (target_se^2 - deff * B))
}

#' Resolve the dispersion and level a pooled estimate is planned on
#'
#' The proportion scale carries the same \eqn{N/(N-1)} adjustment
#' `.change_inputs()` applies and for the same reason: the variance the
#' engine wants is \eqn{S^2} on \eqn{N-1} degrees of freedom.
#' @keywords internal
#' @noRd
.pooled_inputs <- function(var, sd, p, mu, N) {
  has_var <- !is.null(var) || !is.null(sd)
  has_p <- !is.null(p)
  if (has_var + has_p != 1L) {
    stop("specify exactly one of 'var' (or 'sd') or 'p'", call. = FALSE)
  }
  if (has_p) {
    if (!is.null(mu)) {
      stop(
        "'mu' is 'p' on the proportion scale; do not supply it",
        call. = FALSE
      )
    }
    check_proportion(p, "p")
    adj <- if (is.infinite(N)) 1 else N / (N - 1)
    return(list(var = p * (1 - p) * adj, mu = p, p = p))
  }
  var <- .resolve_var(var, sd)
  check_scalar(var, "var")
  if (!is.null(mu)) check_scalar(mu, "mu", positive = FALSE)
  list(var = var, mu = mu, p = NULL)
}

#' Null-coalescing operator
#' Internal fallback for R versions before base::`%||%` became available.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Dispatch helper for svyplan pipe support
#'
#' Intercepts `plan |> fn(arg = val)` calls where R's named-argument
#' matching pushes the svyplan object into `...` instead of dispatching
#' on it.
#' @return The function result, or `NULL` if no interception needed.
#' @keywords internal
#' @noRd
.dispatch_plan <- function(first, first_name, fn, ...) {
  dots <- list(...)
  nms <- names(dots) %||% rep("", length(dots))
  has_named_plan <- any(nms == "plan") &&
    inherits(dots[[match("plan", nms)]], "svyplan")
  unnamed <- which(!nzchar(nms))
  n_unnamed_plan <- sum(
    vapply(dots[unnamed], inherits, logical(1), "svyplan")
  )
  n_total <- inherits(first, "svyplan") + n_unnamed_plan + has_named_plan
  if (n_total > 1L) {
    stop(
      "multiple svyplan objects supplied; use plan= to specify one",
      call. = FALSE
    )
  }
  if (inherits(first, "svyplan")) {
    return(do.call(fn, c(list(plan = first), dots)))
  }
  if (n_unnamed_plan == 1L) {
    plan_idx <- unnamed[vapply(dots[unnamed], inherits, logical(1), "svyplan")]
    plan <- dots[[plan_idx]]
    dots <- dots[-plan_idx]
    return(do.call(fn, c(setNames(list(first), first_name), list(plan = plan), dots)))
  }
  NULL
}

#' Merge plan defaults into a function call
#'
#' Pure function that returns a merged argument list when plan defaults
#' need to be applied, or `NULL` when no merging is needed. Uses
#' `formals()` introspection on the target function to determine which
#' plan defaults are applicable.
#'
#' @param plan A `svyplan` object or `NULL`.
#' @param fn The target function (e.g., `n_prop.default`).
#' @param mc The result of `match.call()` from the caller.
#' @param env The calling function's `environment()`.
#' @return A named list of merged arguments, or `NULL`.
#' @keywords internal
#' @noRd
.merge_plan_args <- function(plan, fn, mc, env) {
  if (is.null(plan)) return(NULL)
  if (!inherits(plan, "svyplan"))
    stop("'plan' must be a svyplan object", call. = FALSE)

  defaults <- plan$defaults
  if (length(defaults) == 0L) return(NULL)

  target_fmls <- setdiff(names(formals(fn)), c("...", "plan"))
  if (!is.null(defaults$prop_method) && "method" %in% target_fmls &&
      !"prop_method" %in% target_fmls) {
    choices <- tryCatch(eval(formals(fn)$method), error = function(e) NULL)
    has_choices <- is.character(choices) && length(choices) > 1L
    if (!has_choices || defaults$prop_method %in% choices) {
      defaults$method <- defaults$prop_method
    }
  }
  applicable <- defaults[names(defaults) %in% target_fmls]
  if (length(applicable) == 0L) return(NULL)

  explicit <- names(mc)[-1L]
  fill <- applicable[!names(applicable) %in% explicit]
  if (length(fill) == 0L) return(NULL)

  args <- mget(target_fmls, envir = env)
  args[names(fill)] <- fill
  args
}

#' Merge round-trip overrides from ... into stored arguments
#'
#' Named dots override the stored values, and a NULL value unsets the
#' stored one (mode switching). Unknown names error instead of being dropped.
#' @keywords internal
#' @noRd
.roundtrip_args <- function(args, dots, fn) {
  if (length(dots) == 0L) {
    return(args)
  }
  nms <- names(dots) %||% rep("", length(dots))
  if (any(!nzchar(nms))) {
    stop("arguments passed via ... must be named", call. = FALSE)
  }
  unknown <- setdiff(nms, setdiff(names(formals(fn)), "..."))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "unused argument%s: %s",
        if (length(unknown) > 1L) "s" else "",
        paste0("'", unknown, "'", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  args[nms] <- dots
  args
}

#' Clamp FPC to 0 when the net sample reaches N (census)
#' @keywords internal
#' @noRd
.clamp_fpc <- function(fpc, n_net, N) {
  if (is.infinite(N)) {
    return(fpc)
  }
  if (n_net >= N) {
    warning(
      "net sample size (",
      round(n_net, 1),
      ") >= population size (",
      N,
      "); FPC set to 0 (census)",
      call. = FALSE
    )
    return(0)
  }
  fpc
}

#' Column names that are near misses for a recognized indicator column
#'
#' Extra columns are allowed, as in the `n_alloc()` frame, so this rejects
#' only names that the rest of the package would lead a reader to expect
#' here: the spelling the `n_alloc()` frame uses for the same quantity,
#' the argument name a neighboring function takes, or the singular of
#' the argument this table is passed as.
#' @keywords internal
#' @noRd
.indicator_column_aliases <- function() {
  c(
    indicator = "name", indicators = "name", label = "name",
    mean = "mu", method = "prop_method",
    icc = "icc_psu", var_ratio = "var_ratio_psu",
    rme = "rmoe", RMoE = "rmoe"
  )
}

#' Take the dispersion of a continuous indicator in either spelling
#'
#' The scalar methods accept `var` or `sd`, so a table written for one
#' should not have to be rewritten for the other. They are different
#' quantities rather than synonyms, so supplying both is an error rather
#' than a preference, and the rest of the code sees only `var`.
#' @keywords internal
#' @noRd
.indicators_var_from_sd <- function(indicators, domains = NULL) {
  nms <- setdiff(names(indicators), domains)
  if (!"sd" %in% nms) {
    return(indicators)
  }
  if ("var" %in% nms) {
    stop("supply the dispersion as 'sd' or as 'var', not both", call. = FALSE)
  }
  values <- indicators$sd
  present <- !is.na(values)
  if (!is.numeric(values) || any(values[present] <= 0) ||
      any(!is.finite(values[present]))) {
    stop("all 'sd' values must be positive and finite", call. = FALSE)
  }
  indicators$var <- values^2
  indicators$sd <- NULL
  indicators
}

#' The estimand each indicator row declares
#'
#' Three markers, one per estimand: a non-missing `p`, `var`, or `r`. Returned
#' as a logical matrix so callers can both count them and name the offenders.
#' @keywords internal
#' @noRd
.estimand_markers <- function(indicators) {
  n <- nrow(indicators)
  absent <- rep(FALSE, n)
  cbind(
    p = if ("p" %in% names(indicators)) !is.na(indicators$p) else absent,
    var = if ("var" %in% names(indicators)) !is.na(indicators$var) else absent,
    r = .is_ratio_row(indicators)
  )
}

#' Require exactly one estimand per indicator row
#'
#' Shared by the sizing and precision paths so a table that `n_multi()`
#' refuses cannot be evaluated by `prec_multi()`.
#' @keywords internal
#' @noRd
.check_estimand_markers <- function(indicators) {
  markers <- .estimand_markers(indicators)
  n_markers <- rowSums(markers)
  if (any(n_markers == 0L)) {
    stop(
      sprintf(
        "row(s) %s must have a non-NA 'p', 'var', or 'r' value",
        paste(which(n_markers == 0L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(n_markers > 1L)) {
    i <- which(n_markers > 1L)[1L]
    stop(
      sprintf(
        "row %d sets %s; each row must have only one of 'p', 'var', or 'r'",
        i, paste(sQuote(colnames(markers)[markers[i, ]]), collapse = " and ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The value a row's relative quantities are relative to
#'
#' Chosen by the row's estimand marker, not by taking the first non-missing
#' of `p`, `mu`, `r`. A ratio row carrying an incidental `mu`, which is
#' legitimate in a mixed table where some other row is a mean, would otherwise
#' have its `rmoe` scaled by that mean.
#' @keywords internal
#' @noRd
.indicator_scale <- function(indicators) {
  out <- rep(NA_real_, nrow(indicators))
  # Precedence r, then p, then mu, each filling only what is still unset. The
  # `is.na(out)` guard is what stops a mean's `mu`, legitimate in a mixed
  # table, from rescaling a ratio row. Validation rejects a row carrying two
  # markers before any scale is reported, so the precedence between them is
  # not observable and no test can pin it.
  for (col in c("r", "p", "mu")) {
    if (!col %in% names(indicators)) {
      next
    }
    take <- is.na(out) & !is.na(indicators[[col]])
    out[take] <- as.numeric(indicators[[col]][take])
  }
  out
}

#' Take a relative margin of error target in a table and make it absolute
#'
#' Follows `.indicators_var_from_sd()`: accept the alternate spelling,
#' refuse it alongside the one it competes with, convert once, and let
#' every path downstream see only `moe`. A row's scale comes from its own
#' estimand marker via `.indicator_scale()`, so an `rmoe` row without one is
#' an error rather than a row that quietly drops out of the target set.
#' @keywords internal
#' @noRd
.indicators_moe_from_rmoe <- function(indicators, domains = NULL) {
  nms <- setdiff(names(indicators), domains)
  if (!"rmoe" %in% nms) {
    return(indicators)
  }
  values <- indicators$rmoe
  if (!is.numeric(values)) {
    stop("'rmoe' values must be numeric", call. = FALSE)
  }
  set <- !is.na(values)
  if (any(values[set] <= 0) || any(!is.finite(values[set]))) {
    stop("'rmoe' values must be positive and finite", call. = FALSE)
  }
  for (other in c("moe", "cv")) {
    if (other %in% nms && any(set & !is.na(indicators[[other]]))) {
      stop(
        sprintf("each row must have only one of 'rmoe' or '%s'", other),
        call. = FALSE
      )
    }
  }
  scale <- .indicator_scale(indicators[, nms, drop = FALSE])
  if (any(set & (is.na(scale) | scale == 0))) {
    stop(
      "'rmoe' rows need the estimand it is relative to: 'p' for a proportion, 'mu' for a mean, or 'r' for a ratio, and not zero",
      call. = FALSE
    )
  }
  moe <- if ("moe" %in% nms) indicators$moe else rep(NA_real_, nrow(indicators))
  moe[set] <- values[set] * abs(scale[set])
  indicators$moe <- moe
  indicators$rmoe <- NULL
  indicators
}

#' Reject indicator columns that would be silently ignored
#'
#' `domains` names are user-chosen, so they are exempt: a domain may
#' legitimately be called `label` or `method`.
#' @keywords internal
#' @noRd
.check_indicator_columns <- function(indicators, domains = NULL) {
  nms <- setdiff(names(indicators), domains)
  aliases <- .indicator_column_aliases()
  hit <- intersect(nms, names(aliases))
  hit <- hit[!aliases[hit] %in% nms]
  if (length(hit) > 0L) {
    stop(
      sprintf(
        "unrecognized indicator column(s): %s. Did you mean %s?",
        paste(sQuote(hit), collapse = ", "),
        paste(sQuote(unname(aliases[hit])), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Per-row degrees of freedom, or NULL when the row carries none
#'
#' `NULL` is what the single-indicator engines take to mean "no adjustment",
#' and an all-NA column is the common case, so an absent or NA entry has to
#' arrive there as `NULL` rather than as `NA`.
#' @keywords internal
#' @noRd
.row_df <- function(indicators, i) {
  if (!"df" %in% names(indicators)) {
    return(NULL)
  }
  value <- indicators$df[i]
  if (is.na(value)) NULL else value
}

#' Per-row minimum expected case counts, as a floor on each row's size
#'
#' Returns one gross size per row, `NA` where the row sets no floor, so the
#' caller can take an elementwise maximum against the sizes precision asks
#' for. A count of positive cases is defined for a proportion only, so a
#' floor on a mean row is refused rather than ignored.
#' @keywords internal
#' @noRd
.multi_min_cases_n <- function(indicators) {
  if (!"min_cases" %in% names(indicators)) {
    return(NULL)
  }
  set <- !is.na(indicators$min_cases)
  if (!any(set)) {
    return(NULL)
  }
  is_prop <- if ("p" %in% names(indicators)) {
    !is.na(indicators$p)
  } else {
    rep(FALSE, nrow(indicators))
  }
  bad <- which(set & !is_prop)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "'min_cases' counts positive cases and applies to proportion rows only. Row(s) %s carry a %s",
        paste(bad, collapse = ", "),
        if (any(.is_ratio_row(indicators)[bad])) "ratio" else "mean"
      ),
      call. = FALSE
    )
  }
  out <- rep(NA_real_, nrow(indicators))
  out[set] <- vapply(
    which(set),
    function(i) {
      .n_from_min_cases(indicators$min_cases[i], indicators$p[i],
                        indicators$resp_rate[i])
    },
    numeric(1L)
  )
  out
}

#' Refuse a sizing-only case floor where there is no size to set
#' @keywords internal
#' @noRd
.stop_min_cases_column <- function(indicators, context) {
  if ("min_cases" %in% names(indicators) && any(!is.na(indicators$min_cases))) {
    stop(
      sprintf("'min_cases' sizes a single-stage indicator table and %s", context),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate common optional columns in a multi-indicator table
#' @keywords internal
#' @noRd
.validate_common_columns <- function(indicators) {
  if ("min_cases" %in% names(indicators)) {
    vals <- indicators$min_cases[!is.na(indicators$min_cases)]
    if (length(vals) > 0L &&
        (!is.numeric(vals) || any(vals <= 0) || any(!is.finite(vals)))) {
      stop("'min_cases' values must be positive and finite", call. = FALSE)
    }
  }
  if ("df" %in% names(indicators)) {
    vals <- indicators$df[!is.na(indicators$df)]
    if (length(vals) > 0L && (!is.numeric(vals) || any(vals <= 0))) {
      stop("'df' values must be positive", call. = FALSE)
    }
  }
  if ("alpha" %in% names(indicators)) {
    vals <- indicators$alpha[!is.na(indicators$alpha)]
    if (any(vals <= 0 | vals >= 1)) {
      stop("'alpha' values must be in (0, 1)", call. = FALSE)
    }
  }
  if ("deff" %in% names(indicators)) {
    vals <- indicators$deff[!is.na(indicators$deff)]
    if (any(vals <= 0) || any(!is.finite(vals))) {
      stop("'deff' values must be positive and finite", call. = FALSE)
    }
  }
  if ("N" %in% names(indicators)) {
    vals <- indicators$N[!is.na(indicators$N)]
    if (any(vals <= 1)) {
      stop("'N' values must be greater than 1 (or Inf)", call. = FALSE)
    }
  }
  for (rate in intersect(
    c("resp_rate_psu", "resp_rate_ssu", "resp_rate"),
    names(indicators)
  )) {
    vals <- indicators[[rate]][!is.na(indicators[[rate]])]
    if (!is.numeric(vals) || any(vals <= 0 | vals > 1)) {
      stop(sprintf("'%s' values must be in (0, 1]", rate), call. = FALSE)
    }
  }
  if ("n" %in% names(indicators)) {
    if (anyNA(indicators$n) || any(indicators$n <= 0) || any(!is.finite(indicators$n))) {
      stop("'n' values must be positive, finite, and non-NA", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Weighted variance
#'
#' Normalizes the weights to sum to 1, then takes the weighted mean of the
#' squared deviations from the weighted center and rescales by n / (n - 1)
#' so unit weights reproduce [var()].
#' @keywords internal
#' @noRd
.wtdvar <- function(x, w) {
  if (!is.numeric(x) || !is.numeric(w) || length(x) != length(w)) {
    stop("'x' and 'w' must be numeric vectors of equal length", call. = FALSE)
  }
  if (length(w) < 2L) {
    stop("'x' and 'w' must have length >= 2", call. = FALSE)
  }
  if (anyNA(x) || anyNA(w) || any(!is.finite(x)) || any(!is.finite(w))) {
    stop(
      "'x' and 'w' must contain only finite, non-missing values",
      call. = FALSE
    )
  }
  n <- length(w)
  total <- sum(w)
  if (total <= 0) {
    stop("sum of weights must be positive", call. = FALSE)
  }
  share <- w / total
  center <- sum(share * x)
  n / (n - 1) * sum(share * (x - center)^2)
}

#' Resolve the dispersion input of the mean-based functions
#'
#' `var` and `sd` are two spellings of the same input, and exactly one is
#' required. `sd` is accepted because stratum frames and published survey
#' reports quote standard deviations rather than variances.
#'
#' @param var Population variance, or `NULL`.
#' @param sd Population standard deviation, or `NULL`.
#' @return The variance.
#' @keywords internal
#' @noRd
.resolve_var <- function(var, sd) {
  if (!is.null(sd)) {
    if (!is.null(var)) {
      stop("supply exactly one of 'var' or 'sd'", call. = FALSE)
    }
    if (!is.numeric(sd) || anyNA(sd) || any(!is.finite(sd)) || any(sd <= 0)) {
      stop("'sd' must contain positive finite values", call. = FALSE)
    }
    return(sd^2)
  }
  if (is.null(var)) {
    stop("supply exactly one of 'var' or 'sd'", call. = FALSE)
  }
  var
}
