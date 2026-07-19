#' Survey Plan Profile
#'
#' Create a reusable profile that captures shared design defaults
#' for survey sample size and power calculations. A plan can be passed
#' to functions like [n_prop()], [n_mean()], [power_prop()], or
#' [n_cluster()] either via the `plan` named argument or by piping:
#'
#' ```
#' plan <- svyplan(deff = 1.8, resp_rate = 0.85, N = 50000)
#'
#' # Named argument
#' n_prop(p = 0.3, moe = 0.05, plan = plan)
#'
#' # Pipe
#' plan |> n_prop(0.3, moe = 0.05)
#' plan |> n_prop(p = 0.3, moe = 0.05)
#' plan |> n_cluster(cv = 0.05)
#' ```
#'
#' @param ... Named defaults to reuse across calls. Allowed names:
#'   `alpha`, `N`, `deff`, `resp_rate`,
#'   `prop_method`,
#'   `stage_cost`, `delta`, `rel_var`, `k`, `fixed_cost`,
#'   `unit_cost`,
#'   `alternative`, `ratio`, `overlap`, `rho`,
#'   `alloc`, `min_n`, `power_q`.
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
#' formals. Irrelevant defaults are silently ignored. `prop_method` must
#' be `"wald"`, `"wilson"`, or `"logodds"` (validated at construction)
#' and also fills the `method` argument of [n_prop()], [prec_prop()], and
#' [power_prop()] when its value is valid for that function (so
#' `svyplan(prop_method = "wilson")` applies to `n_prop()` but is ignored
#' by `power_prop()`, which has no Wilson method).
#'
#' All stored defaults are validated when the profile is created or updated.
#' When both `stage_cost` and `delta` are supplied, their lengths must describe
#' the same number of stages. Length checks that depend on call-specific data,
#' such as matching `unit_cost` to an allocation frame, occur when the plan is
#' used.
#'
#' Estimand-specific values (`p`, `var`, `mu`, `moe`, `cv`, `n`, `power`,
#' `effect`) should be passed directly to each function, not stored in
#' the plan.
#'
#' When piping, the plan is detected automatically. The function's
#' original first argument (e.g., `p` in [n_prop()]) can be passed
#' either positionally or by name.
#'
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
#' cl_plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.85)
#' cl_plan |> n_cluster(cv = 0.05)
#'
#' # Override a plan default
#' n_prop(p = 0.3, moe = 0.05, plan = plan, deff = 2.0)
#'
#' @export
svyplan <- function(...) {
  defaults <- list(...)
  .validate_plan_defaults(defaults)
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
      val_str <- if (length(val) > 1L) {
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
  new_defaults <- modifyList(object$defaults, dots)
  do.call(svyplan, new_defaults)
}

#' @keywords internal
#' @noRd
.svyplan_allowed_defaults <- function() {
  c(
    "alpha", "N", "deff", "resp_rate",
    "prop_method",
    "stage_cost", "delta", "rel_var", "k", "fixed_cost",
    "unit_cost",
    "alternative", "ratio", "overlap", "rho",
    "alloc", "min_n", "power_q"
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
        !pm %in% c("wald", "wilson", "logodds")) {
      stop("'prop_method' must be one of \"wald\", \"wilson\", \"logodds\"",
           call. = FALSE)
    }
  }
  if ("alpha" %in% nms) check_alpha(defaults$alpha)
  if ("N" %in% nms) check_population_size(defaults$N)
  if ("deff" %in% nms) check_deff(defaults$deff)
  if ("resp_rate" %in% nms) check_resp_rate(defaults$resp_rate)

  stage_cost <- NULL
  if ("stage_cost" %in% nms) {
    check_stage_cost(defaults$stage_cost)
    stage_cost <- .reorder_stage_cost(defaults$stage_cost)
  }

  delta <- NULL
  if ("delta" %in% nms) {
    delta <- defaults$delta
    if (inherits(delta, "svyplan_varcomp")) {
      if (!is.null(delta$strata)) {
        stop(
          "stratified varcomp cannot be stored as a cluster plan default",
          call. = FALSE
        )
      }
      delta <- delta$delta
    }
    delta <- .reorder_stage_vec(delta, "delta")
    if (!length(delta) %in% 1:2) {
      stop("'delta' must have length 1 or 2", call. = FALSE)
    }
    check_delta(delta)
    .check_cluster_delta_open(delta, context = "svyplan()")
    if (!is.null(stage_cost)) {
      check_delta(delta, expected_length = length(stage_cost) - 1L)
    }
  }

  if ("rel_var" %in% nms) check_scalar(defaults$rel_var, "rel_var")
  if ("k" %in% nms) {
    k <- .reorder_stage_vec(defaults$k, "k")
    if (!is.numeric(k) || !length(k) %in% 1:2 || anyNA(k) ||
        any(!is.finite(k)) || any(k <= 0)) {
      stop("'k' must contain one or two positive finite values", call. = FALSE)
    }
    if (!is.null(stage_cost) && length(stage_cost) == 2L && length(k) != 1L) {
      stop("'k' must have length 1 for a 2-stage plan", call. = FALSE)
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
  if ("rho" %in% nms) check_rho(defaults$rho)

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
  if ("min_n" %in% nms) check_scalar(defaults$min_n, "min_n")
  if ("power_q" %in% nms) {
    power_q <- defaults$power_q
    if (!is.numeric(power_q) || length(power_q) != 1L || is.na(power_q) ||
        !is.finite(power_q) || power_q < 0 || power_q > 1) {
      stop("'power_q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }
  invisible(NULL)
}
