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
    stop("4+ stage optimization is not yet supported", call. = FALSE)
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

#' Check delta vector
#' @keywords internal
#' @noRd
check_delta <- function(delta, expected_length = NULL) {
  if (!is.numeric(delta) || length(delta) == 0L) {
    stop("'delta' must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(delta) || any(delta < 0) || any(delta > 1)) {
    stop("'delta' values must be in [0, 1]", call. = FALSE)
  }
  if (!is.null(expected_length) && length(delta) != expected_length) {
    stop(
      sprintf("'delta' must have length %d (stages - 1)", expected_length),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Tolerance for detecting numerically degenerate cluster homogeneity
#' @keywords internal
#' @noRd
.cluster_delta_tol <- function() {
  1e-8
}

#' Reject cluster homogeneity values that are too close to 0 or 1
#' @keywords internal
#' @noRd
.check_cluster_delta_open <- function(delta, context = "n_cluster()") {
  tol <- .cluster_delta_tol()
  bad <- delta <= tol | delta >= 1 - tol
  if (any(bad)) {
    stop(
      sprintf(
        "'delta' must stay away from 0 and 1 for %s; values <= %.0e or >= %.8f make cluster optimization degenerate",
        context,
        tol,
        1 - tol
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Reorder a named stage-indexed vector
#'
#'
#' If `x` is unnamed, return as-is. If named, validate and reorder
#' to canonical `c(..._psu, ..._ssu)` order.
#' @param x Numeric vector (length 1 or 2).
#' @param name Parameter name for error messages (e.g., `"delta"`, `"k"`).
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
#' Supports stage names `n_psu`, `psu_size`, `ssu_size`.
#' @keywords internal
#' @noRd
.reorder_n_vec <- function(n) {
  stages <- length(n)
  out_names <- if (stages == 2L) {
    c("n_psu", "psu_size")
  } else {
    c("n_psu", "psu_size", "ssu_size")
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
#' Supports stage names like `delta_psu` / `delta_ssu` (or `k_psu` / `k_ssu`).
#' Optionally supports a scalar alias (`delta` or `k`) for 2-stage.
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
#' @param name Parameter name for error messages (e.g., `"delta"`, `"k"`).
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
check_rho <- function(rho) {
  if (
    !is.numeric(rho) || length(rho) != 1L || anyNA(rho) || rho < 0 || rho > 1
  ) {
    stop("'rho' must be a number in [0, 1]", call. = FALSE)
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
  b1 <- if (is.infinite(N_pair[1])) Inf else N_pair[1] / (ratio * resp_rate)
  b2 <- if (is.infinite(N_pair[2])) Inf else N_pair[2] / resp_rate
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
      sprintf("%s is negative; reduce overlap or rho, or adjust inputs", what),
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
  applicable <- defaults[names(defaults) %in% target_fmls]
  if (length(applicable) == 0L) return(NULL)

  explicit <- names(mc)[-1L]
  fill <- applicable[!names(applicable) %in% explicit]
  if (length(fill) == 0L) return(NULL)

  args <- mget(target_fmls, envir = env)
  args[names(fill)] <- fill
  args
}

#' Clamp FPC to 0 when n_eff >= N (census)
#' @keywords internal
#' @noRd
.clamp_fpc <- function(fpc, n_eff, N) {
  if (is.infinite(N)) {
    return(fpc)
  }
  if (n_eff >= N) {
    warning(
      "effective sample size (",
      round(n_eff, 1),
      ") >= population size (",
      N,
      "); FPC set to 0 (census)",
      call. = FALSE
    )
    return(0)
  }
  fpc
}

#' Validate common optional columns in multi-indicator targets
#' @keywords internal
#' @noRd
.validate_common_columns <- function(targets) {
  if ("alpha" %in% names(targets)) {
    vals <- targets$alpha[!is.na(targets$alpha)]
    if (any(vals <= 0 | vals >= 1)) {
      stop("'alpha' values must be in (0, 1)", call. = FALSE)
    }
  }
  if ("deff" %in% names(targets)) {
    vals <- targets$deff[!is.na(targets$deff)]
    if (any(vals <= 0) || any(!is.finite(vals))) {
      stop("'deff' values must be positive and finite", call. = FALSE)
    }
  }
  if ("N" %in% names(targets)) {
    vals <- targets$N[!is.na(targets$N)]
    if (any(vals <= 1)) {
      stop("'N' values must be greater than 1 (or Inf)", call. = FALSE)
    }
  }
  if ("resp_rate" %in% names(targets)) {
    vals <- targets$resp_rate[!is.na(targets$resp_rate)]
    if (any(vals <= 0 | vals > 1)) {
      stop("'resp_rate' values must be in (0, 1]", call. = FALSE)
    }
  }
  if ("n" %in% names(targets)) {
    if (anyNA(targets$n) || any(targets$n <= 0) || any(!is.finite(targets$n))) {
      stop("'n' values must be positive, finite, and non-NA", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Weighted variance
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
  sw <- sum(w)
  if (sw <= 0) {
    stop("sum of weights must be positive", call. = FALSE)
  }
  xbarw <- sum(w * x) / sw
  n / (n - 1) * sum(w * (x - xbarw)^2) / sw
}
