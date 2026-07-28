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
#' outcome is centred on zero.
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
#' \deqn{k_{ssu} = k_{psu}(1 - \delta_{psu}).}
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
#' A scalar `var_ratio` supplies `var_ratio_psu` only; `var_ratio_ssu` then follows from the
#' decomposition. A length-2 `var_ratio` is taken as supplied and checked for
#' consistency.
#' @keywords internal
#' @noRd
.stage_k_pair <- function(var_ratio, icc) {
  if (length(var_ratio) == 1L) {
    return(c(var_ratio, .var_ratio_ssu_default(var_ratio, icc[1L])))
  }
  rep_len(var_ratio, 2L)
}

#' Check that exactly one of moe/cv is specified
#' @keywords internal
#' @noRd
check_precision <- function(moe, cv) {
  has_moe <- !is.null(moe)
  has_cv <- !is.null(cv)
  if (has_moe == has_cv) {
    stop("specify exactly one of 'moe' or 'cv'", call. = FALSE)
  }
  if (has_moe) {
    check_scalar(moe, "moe")
  }
  if (has_cv) {
    check_scalar(cv, "cv")
  }
  invisible(TRUE)
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

#' Check target CV against the achievable floor; error with diagnostic
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
check_resp_rate <- function(resp_rate) {
  if (
    !is.numeric(resp_rate) ||
      length(resp_rate) != 1L ||
      anyNA(resp_rate) ||
      resp_rate <= 0 ||
      resp_rate > 1
  ) {
    stop("'resp_rate' must be a number in (0, 1]", call. = FALSE)
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

#' Upper bound for n2 from finite-population constraints on effective n
#' @keywords internal
#' @noRd
.n2_upper_bound <- function(N_pair, ratio, resp_rate) {
  b1 <- if (is.infinite(N_pair[1])) Inf else N_pair[1] / ratio
  b2 <- N_pair[2]
  min(b1, b2)
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
  tol = 1e-8
) {
  lo <- sqrt(.Machine$double.eps)
  p_lo <- suppressWarnings(power_fn(lo))
  if (!is.finite(p_lo)) {
    p_lo <- 0
  }
  if (target_power <= p_lo) {
    return(lo)
  }

  hi_cap <- .n2_upper_bound(N_pair, ratio, resp_rate)
  if (is.finite(hi_cap)) {
    hi <- hi_cap * (1 - 1e-7)
    if (hi <= lo) {
      stop(
        "target power is unattainable under finite population constraints",
        call. = FALSE
      )
    }
    p_hi <- suppressWarnings(power_fn(hi))
    if (!is.finite(p_hi) || p_hi < target_power - tol) {
      stop(
        "target power is unattainable under finite population constraints",
        call. = FALSE
      )
    }
  } else {
    hi <- max(2, lo * 2)
    p_hi <- suppressWarnings(power_fn(hi))
    iter <- 0L
    while ((is.na(p_hi) || p_hi < target_power) && hi < 1e12 && iter < 100L) {
      hi <- hi * 2
      p_hi <- suppressWarnings(power_fn(hi))
      iter <- iter + 1L
    }
    if (!is.finite(p_hi) || p_hi < target_power - tol) {
      stop("could not bracket sample size for target power", call. = FALSE)
    }
  }

  uniroot(
    function(n2) power_fn(n2) - target_power,
    interval = c(lo, hi),
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

#' Shared precision engine for proportions
#'
#' One documented variance equation: n_net = n * resp_rate responding
#' units drive both the leading term and the FPC's sampling fraction;
#' deff multiplies the SRSWOR variance at n_net, so
#' se^2 = deff * p * q * fpc(n_net) / n_net with fpc(n) = (N - n)/(N - 1).
#' All three methods read that variance through one quantity, the effective
#' size n_eff = n_net / (deff * fpc) at which an infinite-population SRS
#' would reproduce it. A census drives fpc to 0, n_eff to infinity, and
#' every method's margin of error to 0.
#' @keywords internal
#' @noRd
.prec_engine_prop <- function(p, n, alpha, N, deff, resp_rate, method,
                              df = NULL) {
  z <- qnorm(1 - alpha / 2)
  q <- 1 - p
  n_net <- n * resp_rate
  fpc <- if (is.infinite(N)) 1 else (N - n_net) / (N - 1)
  fpc <- .clamp_fpc(fpc, n_net, N)
  n_eff <- .effective_from_n(n_net, N, deff)
  # A census has no sampling variance under any method. Returning here keeps
  # the single census warning raised by .clamp_fpc() above.
  if (is.infinite(n_eff)) {
    return(list(se = 0, moe = 0, cv = 0))
  }

  if (method == "wald") {
    se <- sqrt(p * q / n_eff)
    moe <- z * se
  } else if (method == "wilson") {
    moe <- .wilson_moe(p, n_eff, z)
    se <- moe / z
  } else if (method == "logodds") {
    moe <- .logodds_moe(p, n_net, alpha, N, deff)
    se <- moe / z
  } else {
    moe <- .beta_moe(p, .kg_effective(n_eff, n_net, alpha, df), alpha)
    se <- moe / z
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
#' Centred at `(n p + z^2/2) / (n + z^2)` with half-width
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
  # With fewer than one degree of freedom the t quantile is undefined and no
  # interval is identified from the design.
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df < 1) {
    stop("'df' must be a number >= 1", call. = FALSE)
  }
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
#' Builds a symmetric interval on the logit scale from the icc-method
#' standard error of `logit(p_hat)`, then maps both endpoints back to the
#' probability scale and halves their distance. The variance of `p_hat`
#' follows the package convention shared with the Wald method,
#' `deff * N / (N - 1) * p q (1/n_net - 1/N)`, so the two agree as the
#' margin of error shrinks.
#' @keywords internal
#' @noRd
.logodds_moe <- function(p, n_net, alpha, N, deff = 1) {
  if (!is.infinite(N) && n_net >= N) {
    warning("net sample size >= population size; moe is 0", call. = FALSE)
    return(0)
  }
  bernoulli <- if (is.infinite(N)) 1 else N / (N - 1)
  fraction <- if (is.infinite(N)) 0 else 1 / N
  var_p <- deff * bernoulli * p * (1 - p) * (1 / n_net - fraction)
  spread <- qnorm(1 - alpha / 2) * sqrt(var_p) / (p * (1 - p))
  centre <- qlogis(p)
  (plogis(centre + spread) - plogis(centre - spread)) / 2
}

#' Shared precision engine for means
#'
#' Same convention as .prec_engine_prop(), with fpc(n) = 1 - n / N.
#' @keywords internal
#' @noRd
.prec_engine_mean <- function(var, mu, n, alpha, N, deff, resp_rate) {
  z <- qnorm(1 - alpha / 2)
  n_net <- n * resp_rate
  n_eff <- n_net / deff
  fpc <- if (is.infinite(N)) 1 else 1 - n_net / N
  fpc <- .clamp_fpc(fpc, n_net, N)
  se <- sqrt(var * fpc / n_eff)
  list(se = se, moe = z * se,
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
#' Vectorized over groups/indicators; `label` names the offender.
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
#' The marginal terms carry their own finite population correction; the
#' overlap covariance does not. For two SRSWOR samples drawn from one
#' population of size \eqn{N} and sharing \eqn{k = overlap \cdot n_1} units,
#' \deqn{Cov(\bar y_1, \bar y_2) = \rho S_1S_2\{k/(n_1n_2) - 1/N\},}
#' so the population term enters once as \eqn{1/N} rather than through
#' sample 2's marginal factor. Collecting terms,
#' \deqn{V = \frac{v_1}{n_1} + \frac{v_2}{n_2}
#'       - \frac{2\rho\,overlap\sqrt{v_1v_2}}{n_2}
#'       - \frac{v_1 + v_2 - 2\rho\sqrt{v_1v_2}}{N},}
#' a per-unit part less a census part. At \eqn{\rho = 1} with equal sizes
#' and variances the census part is zero, the population terms cancel
#' exactly, and the whole reduces to \eqn{2S^2(1 - overlap)/n}.
#' @keywords internal
#' @noRd
.diff_var_fpc <- function(n_eff, var_pair, N_pair, deff, overlap, overlap_cor) {
  if (length(n_eff) == 1L) n_eff <- c(n_eff, n_eff)
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
#' Named dots override the stored values; a NULL value unsets the stored
#' one (mode switching). Unknown names error instead of being dropped.
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
#' the argument name a neighbouring function takes, or the singular of
#' the argument this table is passed as.
#' @keywords internal
#' @noRd
.indicator_column_aliases <- function() {
  c(
    indicator = "name", indicators = "name", label = "name",
    mean = "mu", method = "prop_method",
    icc = "icc_psu", var_ratio = "var_ratio_psu"
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

#' Validate common optional columns in a multi-indicator table
#' @keywords internal
#' @noRd
.validate_common_columns <- function(indicators) {
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
  if ("resp_rate" %in% names(indicators)) {
    vals <- indicators$resp_rate[!is.na(indicators$resp_rate)]
    if (any(vals <= 0 | vals > 1)) {
      stop("'resp_rate' values must be in (0, 1]", call. = FALSE)
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
#' squared deviations from the weighted centre and rescales by n / (n - 1)
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
  centre <- sum(share * x)
  n / (n - 1) * sum(share * (x - centre)^2)
}

#' Resolve the dispersion input of the mean-based functions
#'
#' `var` and `sd` are two spellings of the same input; exactly one is
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
