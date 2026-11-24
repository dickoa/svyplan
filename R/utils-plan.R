#' Refuse a sizing result where a sizing function takes its inputs
#'
#' A result would otherwise reach the default method as if it were a frame
#' or an indicator table, and fail on a missing input.
#' @keywords internal
#' @noRd
.refuse_result_as_input <- function(x, arg, fn) {
  results <- c("svyplan_n", "svyplan_alloc", "svyplan_cluster",
               "svyplan_twophase", "svyplan_power")
  if (!inherits(x, results)) {
    return(invisible(NULL))
  }
  scenarios <- if (inherits(x, "svyplan_twophase")) {
    ""
  } else {
    ", or vary them with predict() on the result"
  }
  stop(
    sprintf(
      "'%s' is a sizing result, not an input. To size again, call %s() with the original inputs and the values to change%s",
      arg, fn, scenarios
    ),
    call. = FALSE
  )
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

#' Name of the function a profile is merged into
#'
#' A profile passed first or by position reaches the method through
#' `do.call()`, so the call carries the function object, not its name.
#' @keywords internal
#' @noRd
.plan_target_name <- function(fn, mc) {
  name <- if (is.function(mc[[1L]])) {
    ns <- asNamespace("svyplan")
    Find(function(nm) identical(ns[[nm]], fn), ls(ns)) %||% "this function"
  } else {
    deparse(mc[[1L]])
  }
  sub("\\.default$", "", name)
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
  if (!is.null(defaults[["eligibility_rate"]]) &&
      defaults[["eligibility_rate"]] != 1 &&
      !"eligibility_rate" %in% target_fmls) {
    stop(
      sprintf(
        "the profile sets 'eligibility_rate', which %s() does not take. Drop it from the profile for this call",
        .plan_target_name(fn, mc)
      ),
      call. = FALSE
    )
  }
  if (!is.null(defaults$eligible_per_unit) && defaults$eligible_per_unit != 1 &&
      !"eligible_per_unit" %in% target_fmls) {
    stop(
      sprintf(
        "the profile sets 'eligible_per_unit', which %s() does not take. Drop it from the profile for this call, or fold it into 'resp_rate' when it is at most 1",
        .plan_target_name(fn, mc)
      ),
      call. = FALSE
    )
  }
  explicit <- names(mc)[-1L]
  by_stage <- c("resp_rate", "N", "stage_n")
  long <- by_stage[vapply(by_stage, function(k) length(defaults[[k]]) > 1L,
                          logical(1))]
  long <- setdiff(intersect(long, target_fmls), explicit)
  if (length(long) > 0L && !"stage_cost" %in% target_fmls) {
    stop(
      sprintf(
        "the profile gives %s one value per stage, and %s() does not read stages that way. Set it in the call, or use a profile without it%s",
        paste(sQuote(long), collapse = " and "),
        .plan_target_name(fn, mc),
        if ("frame" %in% target_fmls && "resp_rate" %in% long) {
          ". A cluster frame takes stage response rates as 'resp_rate_<k>' columns"
        } else {
          ""
        }
      ),
      call. = FALSE
    )
  }
  applicable <- defaults[names(defaults) %in% target_fmls]
  if (length(applicable) == 0L) return(NULL)

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

#' Mark a precision result read back from a fit at its own sizes
#'
#' The print then says the sizes are the fit's continuous optimum. A size
#' passed through `...` makes them typed input, so the mark is dropped.
#' @keywords internal
#' @noRd
.mark_continuous_optimum <- function(out, dots, sizes = "n") {
  if (!any(sizes %in% names(dots))) {
    out$params$continuous_optimum <- TRUE
  }
  out
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
    mu = "mean", method = "prop_method",
    icc = "icc_1", var_ratio = "var_ratio_1",
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
  for (col in c("r", "p", "mean")) {
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
      "'rmoe' rows need the estimand it is relative to: 'p' for a proportion, 'mean' for a mean, or 'r' for a ratio, and not zero",
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
.check_indicator_columns <- function(indicators, domains = NULL,
                                     eligible = FALSE, stages = NULL,
                                     context = "n_multi()") {
  domains <- as.character(domains)
  reserved <- domains[domains == "name" | startsWith(domains, ".")]
  if (length(reserved) > 0L) {
    stop(
      sprintf(
        "domain column(s) %s are named like a result column. Rename them, for example without the leading dot",
        paste(sQuote(reserved), collapse = ", ")
      ),
      call. = FALSE
    )
  }
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
  if (!eligible && "eligibility_rate" %in% nms) {
    stop(
      "an 'eligibility_rate' column applies to n_multi() and prec_multi(). With cluster stages, give the argument 'eligibility_rate', which applies to every row",
      call. = FALSE
    )
  }
  .check_stage_columns(nms, stages, context)
  .check_unique_indicators(indicators, domains)
  invisible(TRUE)
}

#' Refuse indicator columns that the design's stages cannot use
#' @keywords internal
#' @noRd
.check_stage_columns <- function(nms, stages, context) {
  m <- regmatches(nms, regexec("^(icc|var_ratio|N|take|n|cost)_([0-9]+)$",
                               nms))
  hit <- lengths(m) == 3L
  if (is.null(stages)) {
    if (any(hit)) {
      stop(
        sprintf(
          "column(s) %s describe cluster stages, which %s does not read. Size a cluster design with n_cluster(<table>, stage_cost = ), or assess one with prec_cluster(<table>)",
          paste(sQuote(nms[hit]), collapse = ", "), context
        ),
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }
  stem <- vapply(m[hit], `[`, "", 2L)
  k <- as.integer(vapply(m[hit], `[`, "", 3L))
  past <- ifelse(stem == "take", k > stages | k < 2L,
                 ifelse(stem == "n", k != 1L,
                        ifelse(stem == "cost", k > stages, k >= stages)))
  if (any(past)) {
    stop(
      sprintf(
        "column(s) %s are past the last stage of this %d-stage design. Drop them, or give 'stage_cost' one value per stage",
        paste(sQuote(nms[hit][past]), collapse = ", "), stages
      ),
      call. = FALSE
    )
  }
  if ("deff" %in% nms) {
    stop(
      "a 'deff' column is not read by a cluster design, which builds its design effect from icc. Drop it",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Refuse a missing indicator name, or one that labels two rows of a domain
#'
#' Results name rows by `name`, as the binding indicator and as the key of
#' every per-row table, so a missing name or one repeated within a domain
#' makes them ambiguous. Domain columns that are absent are left to their
#' own check.
#' @keywords internal
#' @noRd
.check_unique_indicators <- function(indicators, domains = NULL) {
  if (!"name" %in% names(indicators)) {
    return(invisible(TRUE))
  }
  name <- as.character(indicators$name)
  missing_name <- which(is.na(name) | !nzchar(trimws(name)))
  if (length(missing_name) > 0L) {
    stop(
      sprintf(
        "indicator name missing on row(s) %s. Give each row a 'name', or drop the column to number rows",
        paste(missing_name, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  domain_cols <- intersect(as.character(domains), names(indicators))
  key <- indicators[c("name", domain_cols)]
  repeated <- unique(as.character(key$name[duplicated(key)]))
  if (length(repeated) > 0L) {
    stop(
      sprintf(
        "indicator name(s) %s label more than one row%s. Give each row its own 'name'",
        paste(sQuote(repeated), collapse = ", "),
        if (length(domain_cols) > 0L) " of the same domain" else ""
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
  elig <- .indicator_rate(indicators, "eligible_per_unit", nrow(indicators)) *
    .indicator_rate(indicators, "eligibility_rate", nrow(indicators))
  out <- rep(NA_real_, nrow(indicators))
  out[set] <- vapply(
    which(set),
    function(i) {
      .n_from_min_cases(indicators$min_cases[i], indicators$p[i],
                        indicators$resp_rate[i] * elig[i])
    },
    numeric(1L)
  )
  out
}

#' Refuse a sizing-only case floor where there is no size to set
#' @keywords internal
#' @noRd
.stop_min_cases_column <- function(indicators, context) {
  for (col in c("min_cases", "min_respondents")) {
    if (col %in% names(indicators) && any(!is.na(indicators[[col]]))) {
      stop(
        sprintf("'%s' sizes a single-stage indicator table and %s", col, context),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Rows that carry only a respondent floor
#'
#' Such a row stands for an estimate whose precision is not modeled, such as
#' a fertility rate sized by a count of interviews, so it has no estimand
#' and no target and is set aside by every precision step.
#' @keywords internal
#' @noRd
.floor_only_rows <- function(indicators) {
  if (!"min_respondents" %in% names(indicators)) {
    return(rep(FALSE, nrow(indicators)))
  }
  set <- function(col) {
    if (col %in% names(indicators)) !is.na(indicators[[col]]) else FALSE
  }
  estimand <- set("p") | set("var") | set("sd") | .is_ratio_row(indicators)
  !is.na(indicators$min_respondents) & !estimand
}

#' Validate the respondent floors of an indicator table
#' @keywords internal
#' @noRd
.check_respondent_floors <- function(indicators, floor_only, domains = NULL) {
  if (!"min_respondents" %in% names(indicators)) {
    return(invisible(TRUE))
  }
  vals <- indicators$min_respondents[!is.na(indicators$min_respondents)]
  if (length(vals) > 0L &&
      (!is.numeric(vals) || any(vals <= 0) || any(!is.finite(vals)))) {
    stop("'min_respondents' values must be positive and finite", call. = FALSE)
  }
  if (all(floor_only)) {
    stop(
      "'indicators' needs at least one row with an estimate and a target. A row with only 'min_respondents' raises the size those rows set",
      call. = FALSE
    )
  }
  if (length(domains) > 0L && all(domains %in% names(indicators))) {
    key <- .domain_key(indicators, domains)
    bare_keys <- names(which(tapply(floor_only, key, all)))
    bare <- vapply(bare_keys, function(k) {
      vals <- vapply(indicators[match(k, key), domains, drop = FALSE],
                     as.character, character(1L))
      if (length(domains) == 1L) {
        vals
      } else {
        paste(paste0(domains, " = ", vals), collapse = ", ")
      }
    }, character(1L))
    if (length(bare) > 0L) {
      stop(
        sprintf(
          "domain(s) %s carry only 'min_respondents'. Every domain needs a row with an estimate and a target, which the floor then raises",
          paste(sQuote(bare), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  targets <- intersect(c("cv", "moe", "rmoe"), names(indicators))
  for (col in targets) {
    hit <- which(floor_only & !is.na(indicators[[col]]))
    if (length(hit) > 0L) {
      stop(
        sprintf(
          "row(s) %s carry '%s' but no 'p', 'var' or 'r'. A row with only 'min_respondents' has no estimate to target",
          paste(hit, collapse = ", "), col
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Size each row's respondent floor needs
#'
#' The floor counts eligible respondents, `n * resp_rate * eligible_per_unit`,
#' with no proportion in it, which is what sets it apart from `min_cases`.
#' @keywords internal
#' @noRd
.multi_min_respondents_n <- function(indicators) {
  if (!"min_respondents" %in% names(indicators) ||
      all(is.na(indicators$min_respondents))) {
    return(NULL)
  }
  elig <- .indicator_rate(indicators, "eligible_per_unit", nrow(indicators)) *
    .indicator_rate(indicators, "eligibility_rate", nrow(indicators))
  indicators$min_respondents / (indicators$resp_rate * elig)
}

#' Validate common optional columns in a multi-indicator table
#' @keywords internal
#' @noRd
.validate_common_columns <- function(indicators) {
  ## A column that is all NA is absent, whatever type data.frame() gave it.
  ## Any other column must be numeric before its values are checked.
  values <- function(name) {
    col <- indicators[[name]]
    if (is.null(col) || all(is.na(col))) {
      return(NULL)
    }
    if (!is.numeric(col)) {
      stop(sprintf("'%s' must be a numeric column", name), call. = FALSE)
    }
    col[!is.na(col)]
  }
  vals <- values("min_cases")
  if (length(vals) > 0L && (any(vals <= 0) || any(!is.finite(vals)))) {
    stop("'min_cases' values must be positive and finite", call. = FALSE)
  }
  vals <- values("df")
  if (length(vals) > 0L && any(vals < 1)) {
    stop("'df' values must be numbers >= 1", call. = FALSE)
  }
  vals <- values("alpha")
  if (any(vals <= 0 | vals >= 1)) {
    stop("'alpha' values must be in (0, 1)", call. = FALSE)
  }
  vals <- values("deff")
  if (any(vals <= 0) || any(!is.finite(vals))) {
    stop("'deff' values must be positive and finite", call. = FALSE)
  }
  vals <- values("N")
  if (any(vals <= 1)) {
    stop("'N' values must be greater than 1 (or Inf)", call. = FALSE)
  }
  rate_cols <- c(
    grep("^resp_rate_[0-9]+$", names(indicators), value = TRUE),
    intersect(c("resp_rate", "eligibility_rate"), names(indicators))
  )
  for (rate in rate_cols) {
    vals <- values(rate)
    if (any(vals <= 0 | vals > 1)) {
      stop(sprintf("'%s' values must be in (0, 1]", rate), call. = FALSE)
    }
  }
  if ("n" %in% names(indicators)) {
    n <- indicators$n
    if (!is.numeric(n) || anyNA(n) || any(n <= 0) || any(!is.finite(n))) {
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
