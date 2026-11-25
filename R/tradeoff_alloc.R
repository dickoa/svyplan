#' Trade the precision of one objective component for another's
#'
#' Tighten the CV of one component of a budget fit and see what a second
#' component gives up, at the same budget, with every other component held
#' at its baseline CV.
#'
#' Each row refits `object` with the CV of `improve` capped at `cv_mult`
#' times its baseline, the CV of every component other than `allow` capped
#' at its baseline, and the relative variance of `allow` as the objective.
#' The fit's targets, stratum bounds, sample floors and budget all hold. A
#' component whose baseline CV is 0, a census, keeps its strata at their
#' upper counts. The refit is a sum-rule budget fit, so its continuous and
#' whole designs, statuses and search bound mean what they mean in
#' [n_alloc_joint()].
#'
#' With `baseline = "continuous"` and slopes the fit reports as
#' `"available"`, `estimated_allow_cv` is the first-order estimate from the
#' frontier's normal: the priorities under the sum rule, and
#' `$optimization$objective_weights` under the worst rule. It is `NA`
#' otherwise.
#'
#' @param object A budget fit from [n_alloc_joint()] on an element frame.
#' @param improve,allow Two different components of `object$objective$component`.
#' @param cv_mult Positive multipliers of the baseline CV of `improve`, each
#'   at most 1. One row per value.
#' @param baseline `"field"` (default) compares with the whole design of
#'   `object`, `"continuous"` with its continuous design.
#' @param control Search settings as in [n_alloc_joint()]. The default is the
#'   fit's own.
#' @return A data frame, one row per `cv_mult`, with `requested_cv` for
#'   `improve` and `baseline_allow_cv`, `estimated_allow_cv`, the continuous
#'   design's `continuous_improve_cv`, `continuous_allow_cv`,
#'   `continuous_cost` and `continuous_status`, the whole design's
#'   `field_improve_cv`, `field_allow_cv`, `field_cost`, `field_status` and
#'   `field_allow_lower_bound`, a CV no whole design within the budget beats,
#'   and the list columns `components` (each component's baseline, cap,
#'   continuous and field CV) and `n_continuous` and `n_field` (the
#'   allocations, `NULL` when none). A status is `"optimal"`, `"feasible"`,
#'   `"limit"` or `"gap"` as in [n_alloc_joint()], `"infeasible"` when no
#'   design meets the caps within the budget, or `"unknown"` when the solve
#'   could not tell.
#' @seealso [n_alloc_joint()], [prec_alloc_joint()]
#' @examples
#' frame <- data.frame(stratum = c("A", "B", "C"), N = c(4000, 3000, 2000),
#'                     unit_cost = c(1, 2, 3))
#' measures <- data.frame(
#'   stratum = rep(frame$stratum, 2),
#'   name = rep(c("income", "coverage"), each = 3),
#'   mean = c(50, 60, 70, NA, NA, NA), sd = c(20, 25, 40, NA, NA, NA),
#'   p = c(NA, NA, NA, 0.3, 0.5, 0.2)
#' )
#' fit <- n_alloc_joint(frame, measures, objective = c("income", "coverage"),
#'                      budget = 900)
#' tradeoff_alloc(fit, improve = "income@.overall",
#'                allow = "coverage@.overall", cv_mult = c(0.98, 0.95))[
#'   c("cv_mult", "requested_cv", "field_improve_cv", "field_allow_cv",
#'     "field_cost", "field_status")]
#' @export
tradeoff_alloc <- function(object, improve, allow, cv_mult = 0.99,
                           baseline = c("field", "continuous"),
                           control = NULL) {
  baseline <- match.arg(baseline)
  p <- object$params
  if (!inherits(object, "svyplan_alloc_joint") ||
        !identical(p$mode, "budget_objective") ||
        !identical(object$stages, 1L) || !is.null(p$psu) ||
        isTRUE(p$free_take)) {
    stop("'object' must be a budget fit from n_alloc_joint() on an element frame",
         call. = FALSE)
  }
  ids <- object$objective$component
  if (!is.character(improve) || !is.character(allow) ||
        length(improve) != 1L || length(allow) != 1L ||
        !improve %in% ids || !allow %in% ids || improve == allow) {
    stop("'improve' and 'allow' must be two different components of object$objective$component",
         call. = FALSE)
  }
  if (!is.numeric(cv_mult) || length(cv_mult) == 0L ||
        any(!is.finite(cv_mult) | cv_mult <= 0 | cv_mult > 1)) {
    stop("'cv_mult' must hold positive multipliers of at most 1",
         call. = FALSE)
  }
  i <- match(improve, ids)
  j <- match(allow, ids)
  base <- if (baseline == "field") object$operational$objective else
    object$objective
  base_cv <- base$.cv[match(ids, base$component)]
  rate <- .tradeoff_rate(object, i, j, baseline)
  problem <- p$problem
  rows <- lapply(cv_mult, function(m) {
    requested <- base_cv
    requested[i] <- base_cv[i] * m
    out <- .tradeoff_refit(object, requested, j, control)
    cv_at <- function(x) {
      if (is.null(x)) return(rep(NA_real_, length(ids)))
      .bethel_objective_value(problem, x)$components$.cv
    }
    cont <- cv_at(out$x)
    field <- cv_at(out$whole)
    components <- object$objective[c("component", "name", "domain", "level")]
    components$baseline_cv <- base_cv
    components$cap_cv <- requested
    components$cap_cv[j] <- NA_real_
    components$continuous_cv <- cont
    components$field_cv <- field
    data.frame(
      cv_mult = m,
      requested_cv = requested[i],
      baseline_allow_cv = base_cv[j],
      estimated_allow_cv = sqrt(base_cv[j]^2 +
                                  rate * (base_cv[i]^2 - requested[i]^2)),
      continuous_improve_cv = cont[i],
      continuous_allow_cv = cont[j],
      continuous_cost = out$cost %||% NA_real_,
      continuous_status = out$continuous_status,
      field_improve_cv = field[i],
      field_allow_cv = field[j],
      field_cost = out$field_cost %||% NA_real_,
      field_status = out$field_status,
      field_allow_lower_bound = out$bound %||% NA_real_,
      components = I(list(components)),
      n_continuous = I(list(out$x)),
      n_field = I(list(out$whole)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' The local rate at which `allow`'s relative variance buys `improve`'s
#'
#' Along the frontier at a fixed budget, with the other components held,
#' `normal_i d relvar_i + normal_j d relvar_j = 0`. `NA` unless the
#' continuous baseline has slopes and both normals are positive.
#' @keywords internal
#' @noRd
.tradeoff_rate <- function(object, i, j, baseline) {
  if (!identical(baseline, "continuous") ||
        !identical(object$optimization$slopes, "available")) {
    return(NA_real_)
  }
  normal <- if (identical(object$params$objective_rule %||% "sum", "worst")) {
    unname(object$optimization$objective_weights)
  } else {
    object$objective$priority
  }
  if (length(normal) < max(i, j) || !all(is.finite(normal[c(i, j)])) ||
        !all(normal[c(i, j)] > 0)) {
    return(NA_real_)
  }
  normal[i] / normal[j]
}

#' One trade-off refit through the joint engine
#'
#' The caps join the fit's targets as cv rows, tightening a cv target the
#' fit already holds on the same estimate. A component at CV 0 holds its
#' strata at their upper counts. The continuous design is solved first, so
#' a whole design that the budget cannot buy leaves it in place.
#' @keywords internal
#' @noRd
.tradeoff_refit <- function(object, requested, allow, control) {
  p <- object$params
  obj <- object$objective
  problem <- p$problem
  targets <- p$targets
  if (is.null(targets) || nrow(targets) == 0L) {
    targets <- data.frame(name = character(0), domain = character(0),
                          level = character(0), cv = numeric(0),
                          constraint = character(0), stringsAsFactors = FALSE)
  }
  if (!"cv" %in% names(targets)) targets$cv <- NA_real_
  min_n <- .resolve_min_n_stratum(p$min_n_stratum, p$frame)
  lower <- rep_len(if (is.null(min_n)) NA_real_ else min_n,
                   length(problem$cost))
  added <- list()
  for (k in setdiff(seq_along(requested), allow)) {
    if (!(requested[k] > 0)) {
      held <- problem$objective$A[, k] > 0
      lower[held] <- problem$upper[held] * problem$stage$take[held]
      next
    }
    same <- which(targets$name == obj$name[k] &
                    targets$domain == obj$domain[k] &
                    (is.na(targets$level) & is.na(obj$level[k]) |
                       targets$level %in% obj$level[k]) &
                    !is.na(targets$cv))
    if (length(same)) {
      targets$cv[same] <- pmin(targets$cv[same], requested[k])
    } else {
      added[[length(added) + 1L]] <- data.frame(
        name = obj$name[k], domain = obj$domain[k], level = obj$level[k],
        cv = requested[k], constraint = paste0(".tradeoff:", obj$component[k]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(added)) {
    extra <- do.call(rbind, added)
    for (col in setdiff(names(targets), names(extra))) extra[[col]] <- NA
    for (col in setdiff(names(extra), names(targets))) targets[[col]] <- NA
    targets <- rbind(targets, extra[names(targets)])
  }
  if (nrow(targets) == 0L) targets <- NULL
  objective <- data.frame(name = obj$name[allow], domain = obj$domain[allow],
                          level = obj$level[allow], priority = 1,
                          stringsAsFactors = FALSE)
  sc <- .search_control(control %||% p[["control"]])
  run <- function(whole) {
    tryCatch(
      .n_alloc_bethel(
        frame = p$frame, measures = p$measures, targets = targets,
        unit_cost = p$unit_cost, alpha = p$alpha, deff = p$deff,
        resp_rate = p$resp_rate,
        min_n_stratum = if (any(!is.na(lower))) lower else NULL,
        objective = objective, budget = p$budget, df = p$df,
        fpc = p$fpc %||% "unit", .finite_prop_var = TRUE,
        eligible_per_unit = p$eligible_per_unit %||% 1,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        max_nodes = sc$max_nodes, gap = sc$gap, .whole = whole,
        objective_rule = "sum", min_n_domain = p$min_n_domain
      ),
      error = function(e) e
    )
  }
  status_of <- function(e) {
    if (inherits(e, "svyplan_plan_unknown")) return("unknown")
    msg <- conditionMessage(e)
    proven <- c("infeasible", "cannot fund", "no integer allocation meets",
                "no feasible integer allocation", "exceeds")
    if (any(vapply(proven, grepl, logical(1), x = msg, fixed = TRUE))) {
      "infeasible"
    } else {
      "unknown"
    }
  }
  out <- list(continuous_status = "unknown", field_status = "unknown")
  cont <- run(FALSE)
  if (inherits(cont, "error")) {
    out$continuous_status <- status_of(cont)
    out$field_status <- out$continuous_status
    return(out)
  }
  out$x <- stats::setNames(cont$detail$n, cont$detail$stratum)
  out$cost <- cont$params$achieved$cost
  out$continuous_status <- cont$optimization$classification
  fit <- run(TRUE)
  if (inherits(fit, "error")) {
    out$field_status <- status_of(fit)
    return(out)
  }
  out$whole <- stats::setNames(as.numeric(fit$detail$n_int),
                               fit$detail$stratum)
  out$field_cost <- fit$operational$cost
  out$field_status <- fit$operational$status
  out$bound <- sqrt(max(0, fit$optimization$search$bound))
  out
}
