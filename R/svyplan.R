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
#'   `stage_cost`, `delta`, `rel_var`, `k`, `fixed_cost`,
#'   `unit_cost`,
#'   `alternative`, `method`, `ratio`, `overlap`, `rho`,
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
#' formals. Irrelevant defaults are silently ignored.
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
    "stage_cost", "delta", "rel_var", "k", "fixed_cost",
    "unit_cost",
    "alternative", "method", "ratio", "overlap", "rho",
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

  if ("alpha" %in% nms) check_alpha(defaults$alpha)
  if ("N" %in% nms && any(is.finite(defaults$N))) check_population_size(defaults$N)
  if ("deff" %in% nms) check_deff(defaults$deff)
  if ("resp_rate" %in% nms) check_resp_rate(defaults$resp_rate)
  if ("stage_cost" %in% nms) check_stage_cost(defaults$stage_cost)
  invisible(NULL)
}
