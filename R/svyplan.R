#' Survey plan profile
#'
#' Create a reusable profile that captures shared design defaults
#' for survey sample size and power calculations. A plan can be passed
#' to functions like [n_prop()], [n_mean()], [power_prop()], or
#' [n_cluster()] either through the `plan` argument or by piping the plan
#' in as the first argument, as the examples show.
#'
#' @param ... Named defaults to reuse across calls. Allowed names:
#'   `alpha`, `df`, `N`, `deff`, `resp_rate`, `eligibility_rate`,
#'   `eligible_per_unit`,
#'   `prop_method`,
#'   `stage_cost`, `stage_n`, `icc`, `unit_relvar`, `var_ratio`, `fixed_cost`,
#'   `unit_cost`,
#'   `alternative`, `ratio`, `overlap`, `overlap_cor`, `overlap_draw`,
#'   `alloc`, `min_n_stratum`, `min_n_domain`, `alloc_q`, `phase1_cost`.
#'
#' @return A `svyplan` object.
#'
#' @details
#' The profile stores design parameters that are shared across multiple
#' functions. Each default fills in only arguments not explicitly provided
#' by the caller. Explicit arguments always take precedence:
#'
#' ```
#' plan <- svyplan(deff = 1.8, resp_rate = 0.85)
#' n_prop(p = 0.3, moe = 0.05, plan = plan)           # uses plan defaults
#' n_prop(p = 0.3, moe = 0.05, plan = plan, deff = 2)  # deff = 2 wins
#' ```
#'
#' Defaults are applied only when their names match the called function's
#' formals. Irrelevant defaults are silently ignored, except
#' `eligibility_rate` and `eligible_per_unit`. A function without
#' `eligibility_rate` would size from the units issued as if all were in
#' scope, and one without `eligible_per_unit` would count eligible persons
#' where the profile counts households, so either stops instead.
#' `prop_method` must
#' be `"wald"`, `"wilson"`, `"logodds"`, or `"beta"` (validated at
#' construction)
#' and also fills the `method` argument of [n_prop()], [prec_prop()], and
#' [power_prop()] when its value is valid for that function (so
#' `svyplan(prop_method = "wilson")` applies to `n_prop()` but is ignored
#' by `power_prop()`, which has no Wilson method).
#'
#' `resp_rate`, `N` and `stage_n` may hold one value per stage, for
#' [n_cluster()] and [prec_cluster()]. Only those two read them by stage.
#' Any other function given such a profile stops rather than read the
#' values as something else, such as one response rate per stratum in
#' [n_alloc()].
#'
#' All stored defaults are validated when the profile is created or updated.
#' When both `stage_cost` and `icc` are supplied, their lengths must describe
#' the same number of stages. Length checks that depend on call-specific data,
#' such as matching `unit_cost` to an allocation frame, occur when the plan is
#' used. `unit_cost` is ordered by allocation-frame row, whereas
#' [strata_bound()] orders costs from the lowest to the highest stratum, so
#' only a scalar `unit_cost` reaches that function from a profile. A vector
#' one is rejected there rather than applied to the wrong strata.
#'
#' Estimand-specific values (`p`, `var`, `mu`, `moe`, `cv`, `n`, `power`,
#' `effect`) should be passed directly to each function, not stored in
#' the plan.
#'
#' When piping, the plan is detected automatically. The function's
#' original first argument (e.g., `p` in [n_prop()]) can be passed
#' either positionally or by name.
#'
#' @family design input functions
#' @seealso [n_prop()], [n_mean()], [n_cluster()], [power_prop()],
#'   [power_mean()], [power_did()].
#'
#' @examples
#' # Create a plan with common design parameters
#' plan <- svyplan(deff = 1.8, resp_rate = 0.85, N = 50000)
#'
#' # Use via named argument
#' n_prop(p = 0.3, moe = 0.05, plan = plan)
#' n_mean(var = 100, mu = 50, cv = 0.05, plan = plan)
#' power_prop(p1 = 0.30, p2 = 0.35, plan = plan)
#' n_multi(data.frame(p = 0.05, moe = 0.02), plan = svyplan(prop_method = "wilson"))
#'
#' # Use via pipe
#' plan |> n_prop(0.3, moe = 0.05)
#' plan |> n_prop(p = 0.3, moe = 0.05)
#' plan |> n_mean(100, mu = 50, cv = 0.05)
#' plan |> power_prop(0.30, p2 = 0.35)
#'
#' # Cluster context
#' cl_plan <- svyplan(stage_cost = c(500, 50), icc = 0.05, resp_rate = 0.85)
#' cl_plan |> n_cluster(cv = 0.05)
#'
#' # Override a plan default
#' n_prop(p = 0.3, moe = 0.05, plan = plan, deff = 2.0)
#'
#' @export
svyplan <- function(...) {
  defaults <- list(...)
  .validate_plan_defaults(defaults)
  ## An explicit NULL unsets a default, as it does in update().
  defaults <- defaults[!vapply(defaults, is.null, logical(1))]
  structure(list(defaults = defaults), class = "svyplan")
}

#' @export
print.svyplan <- function(x, ...) {
  .check_unused_dots(...)
  cat("svyplan profile\n")
  if (length(x$defaults) == 0L) {
    cat("  (no defaults)\n")
  } else {
    nms <- names(x$defaults)
    max_w <- max(nchar(nms))
    for (nm in nms) {
      val <- x$defaults[[nm]]
      val_str <- if (is.object(val) && !is.null(utils::getS3method(
        "format", class(val)[1L], optional = TRUE))) {
        paste(format(val), collapse = " ")
      } else if (length(val) > 1L) {
        paste0("c(", paste(val, collapse = ", "), ")")
      } else {
        as.character(val)
      }
      cat(sprintf("  %-*s = %s\n", max_w, nm, val_str))
    }
  }
  invisible(x)
}

#' @rdname svyplan
#' @param object A `svyplan` object to update.
#' @export
update.svyplan <- function(object, ...) {
  dots <- list(...)
  nms <- names(dots)
  if (length(dots) > 0L &&
        (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms))) {
    stop("arguments to update() must be named, each name once",
         call. = FALSE)
  }
  new_defaults <- modifyList(object$defaults, dots)
  do.call(svyplan, new_defaults)
}

#' @keywords internal
#' @noRd
.svyplan_allowed_defaults <- function() {
  c(
    "alpha", "N", "deff", "resp_rate", "eligibility_rate",
    "eligible_per_unit", "df",
    "prop_method",
    "stage_cost", "stage_n", "icc", "unit_relvar", "var_ratio", "fixed_cost",
    "unit_cost",
    "alternative", "ratio", "overlap", "overlap_cor", "overlap_draw",
    "alloc", "min_n_stratum", "min_n_domain", "alloc_q",
    "phase1_cost"
  )
}

#' @keywords internal
#' @noRd
.validate_plan_defaults <- function(defaults) {
  if (length(defaults) == 0L) return(invisible(NULL))

  nms <- names(defaults)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("all arguments to svyplan() must be named", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("default names must be unique", call. = FALSE)
  }

  allowed <- .svyplan_allowed_defaults()
  bad <- setdiff(nms, allowed)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "unknown default(s): %s. Allowed: %s",
        paste(bad, collapse = ", "),
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if ("prop_method" %in% nms) {
    pm <- defaults$prop_method
    if (!is.character(pm) || length(pm) != 1L ||
        !pm %in% c("wald", "wilson", "logodds", "beta")) {
      stop(
        "'prop_method' must be one of \"wald\", \"wilson\", \"logodds\", \"beta\"",
        call. = FALSE
      )
    }
  }
  if ("alpha" %in% nms) check_alpha(defaults$alpha)
  if ("N" %in% nms && !is.null(defaults$N)) {
    if (length(defaults$N) == 0L) {
      stop("'N' must have at least one value", call. = FALSE)
    }
    if (length(defaults$N) == 1L) {
      check_population_size(defaults$N)
    } else {
      .stage_population(defaults$N, length(defaults$N))
    }
  }
  if ("deff" %in% nms) check_deff(defaults$deff)
  if ("resp_rate" %in% nms) {
    rr <- defaults$resp_rate
    if (!is.numeric(rr) || length(rr) == 0L || anyNA(rr) ||
        any(!is.finite(rr)) || any(rr <= 0) || any(rr > 1)) {
      stop("'resp_rate' must contain values in (0, 1]", call. = FALSE)
    }
  }
  if ("eligibility_rate" %in% nms)
    check_resp_rate(defaults$eligibility_rate, "eligibility_rate")
  if ("eligible_per_unit" %in% nms) check_eligible_per_unit(defaults$eligible_per_unit)
  # design_df() returns a classed numeric, so a plan can store the count a
  # plan object already knows: svyplan(df = design_df(alloc)).
  if ("df" %in% nms && !is.null(defaults$df)) {
    if (!is.numeric(defaults$df)) {
      stop("'df' must be a single number >= 1", call. = FALSE)
    }
    check_df(as.double(defaults$df))
  }

  stage_cost <- NULL
  if ("stage_cost" %in% nms) {
    check_stage_cost(defaults$stage_cost)
    stage_cost <- .reorder_stage_cost(defaults$stage_cost)
  }

  if ("stage_n" %in% nms) {
    sn <- defaults$stage_n
    if (!is.numeric(sn) || length(sn) < 2L ||
        any(!is.na(sn) & (!is.finite(sn) | sn < 1))) {
      stop("'stage_n' must be a numeric vector of length >= 2 with NA or values >= 1",
           call. = FALSE)
    }
    if (!is.null(stage_cost) && length(sn) != length(stage_cost)) {
      stop("'stage_n' must have the same length as 'stage_cost'", call. = FALSE)
    }
  }

  if ("resp_rate" %in% nms && length(defaults$resp_rate) > 1L &&
      !is.null(stage_cost) && length(defaults$resp_rate) != length(stage_cost)) {
    stop("'resp_rate' must be a scalar or have the same length as 'stage_cost'",
         call. = FALSE)
  }

  if ("N" %in% nms && length(defaults$N) > 1L && !is.null(stage_cost) &&
      length(defaults$N) != length(stage_cost)) {
    stop("'N' must have the same length as 'stage_cost'", call. = FALSE)
  }

  icc <- NULL
  if ("icc" %in% nms) {
    icc <- defaults$icc
    if (inherits(icc, "svyplan_varcomp")) {
      if (!is.null(icc$strata)) {
        stop(
          "stratified varcomp cannot be stored as a cluster plan default",
          call. = FALSE
        )
      }
      icc <- icc$icc
    }
    icc <- .reorder_stage_vec(icc, "icc")
    check_icc(icc)
    .check_cluster_icc_open(icc, context = "svyplan()")
    if (!is.null(stage_cost)) {
      check_icc(icc, expected_length = length(stage_cost) - 1L)
    }
  }

  if ("unit_relvar" %in% nms) check_scalar(defaults$unit_relvar, "unit_relvar")
  if ("var_ratio" %in% nms) {
    var_ratio <- .reorder_stage_vec(defaults$var_ratio, "var_ratio")
    if (!is.numeric(var_ratio) || length(var_ratio) == 0L || anyNA(var_ratio) ||
        any(!is.finite(var_ratio)) || any(var_ratio <= 0)) {
      stop("'var_ratio' must contain positive finite values", call. = FALSE)
    }
    if (!is.null(stage_cost) && !length(var_ratio) %in% c(1L, length(stage_cost) - 1L)) {
      stop(
        sprintf(
          "'var_ratio' must have length 1 or %d for a %d-stage plan",
          length(stage_cost) - 1L, length(stage_cost)
        ),
        call. = FALSE
      )
    }
  }
  if ("fixed_cost" %in% nms) check_fixed_cost(defaults$fixed_cost)
  if ("unit_cost" %in% nms) check_weights(defaults$unit_cost, "unit_cost")

  if ("alternative" %in% nms) {
    alternative <- defaults$alternative
    if (!is.character(alternative) || length(alternative) != 1L ||
        is.na(alternative) ||
        !alternative %in% c("two.sided", "one.sided")) {
      stop("'alternative' must be 'two.sided' or 'one.sided'", call. = FALSE)
    }
  }
  if ("ratio" %in% nms) check_scalar(defaults$ratio, "ratio")
  if ("overlap" %in% nms) check_overlap(defaults$overlap)
  if ("overlap_cor" %in% nms) check_overlap_cor(defaults$overlap_cor)
  if ("overlap_draw" %in% nms) {
    overlap_draw <- defaults$overlap_draw
    if (!is.character(overlap_draw) || length(overlap_draw) != 1L ||
        is.na(overlap_draw) ||
        !overlap_draw %in% c("fixed", "independent")) {
      stop("'overlap_draw' must be 'fixed' or 'independent'", call. = FALSE)
    }
  }

  if ("alloc" %in% nms) {
    alloc <- defaults$alloc
    if (!is.character(alloc) || length(alloc) != 1L || is.na(alloc) ||
        !alloc %in% c("neyman", "optimal", "proportional", "power")) {
      stop(
        "'alloc' must be 'neyman', 'optimal', 'proportional', or 'power'",
        call. = FALSE
      )
    }
  }
  if ("min_n_stratum" %in% nms) check_scalar(defaults$min_n_stratum, "min_n_stratum")
  if ("min_n_domain" %in% nms) {
    check_scalar(defaults$min_n_domain, "min_n_domain")
  }
  if ("alloc_q" %in% nms) {
    alloc_q <- defaults$alloc_q
    if (!is.numeric(alloc_q) || length(alloc_q) != 1L || is.na(alloc_q) ||
        !is.finite(alloc_q) || alloc_q < 0 || alloc_q > 1) {
      stop("'alloc_q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }
  if ("phase1_cost" %in% nms) check_scalar(defaults$phase1_cost, "phase1_cost")
  invisible(NULL)
}
