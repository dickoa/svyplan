#' Summarize stratified allocation results
#'
#' `summary()` separates a classic [n_alloc()] or [prec_alloc()] result into
#' the overall answer, the per-stratum allocation, achieved precision,
#' allocation bounds, domains, and assumptions. The design summary evaluates
#' the whole-unit field allocation. The precision summary evaluates the
#' supplied allocation exactly. Joint results from [n_alloc_joint()] or
#' [prec_alloc_joint()] have a diagnostic summary with constraint,
#' objective, optimization, and bound tables.
#'
#' @param object A classic result from [n_alloc()] or [prec_alloc()], or a
#'   joint result from [n_alloc_joint()] or [prec_alloc_joint()].
#' @param x A `summary.svyplan_alloc` or `summary.svyplan_bethel` object.
#' @param ... Additional arguments are not supported.
#'
#' @return For a classic allocation, `summary()` returns a
#'   `summary.svyplan_alloc` object with fields `kind`, `question`, `method`,
#'   `mode`, `overall`, `continuous`, `allocation`, `precision`, `bounds`,
#'   `domains`, and `assumptions`. For a joint allocation or assessment it
#'   returns a `summary.svyplan_bethel` object as described in Details.
#'   Numeric quantities remain unformatted. Formatting is applied by the summary print
#'   methods. `print()` returns the summary invisibly. The summaries of
#'   other `svyplan_n` results are described in [summary.svyplan_n]. Any
#'   other `svyplan_prec` result has no structured summary yet, and
#'   `summary()` returns the result itself, which prints as `print()` does.
#'
#' @details
#' For `n_alloc()`, `overall`, `allocation`, `precision`, and `domains`
#' describe `$detail$n_int` and the operational cluster takes when present.
#' The mathematical optimum remains in `continuous`. For `prec_alloc()`, the
#' same fields describe `$detail$n`, the allocation actually supplied to the
#' precision calculation. `continuous` and `bounds` are `NULL`.
#'
#' Stratum standard errors, margins of error, and CVs do not add to their
#' overall counterparts. `variance_share` is the stratum's contribution to
#' the overall variance and sums to one when that variance is positive.
#' Domain rows are separate subpopulation assessments, not contributions to
#' add to the overall row.
#'
#' A joint allocation summary has fields `kind`, `question`, `mode`,
#' `stages`, `status`, `overall`, `continuous`, `allocation`, `constraints`,
#' `operational_constraints`, `objective`, `operational_objective`, `bounds`,
#' `optimization`, and `assumptions`, plus `sensitivity_note` in a budget
#' fit and in a fit whose slopes are withheld. For a fitted design the headline and
#' allocation describe the integer field recommendation, while `continuous`
#' and `constraints` describe the continuous design. For a precision
#' assessment they describe the supplied allocation exactly.
#'
#' Constraint rows are simultaneous requirements, not an additive
#' decomposition. `.ratio` is achieved divided by target, `.residual` is that
#' ratio minus one, and `.sensitivity` is the local change in minimum cost per
#' unit relaxation of the stated target, so `-.sensitivity * 0.1 * .target`
#' is the first-order saving from relaxing it by 10 percent. In a budget
#' fit the objective is held fixed for `.sensitivity`, and
#' `.objective_sensitivity` is the local change in the objective per unit
#' relaxation with the budget held fixed. `sensitivity_note` states both,
#' and why slopes are missing when they are. A row with no slope shares its
#' price with another active row. Objective
#' `.contribution` and
#' `.share` columns are genuinely additive across components.
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("Urban", "Rural", "Remote"),
#'   N = c(12000, 30000, 8000),
#'   sd = c(5, 8, 12),
#'   mean = c(20, 18, 15),
#'   unit_cost = c(1, 1.4, 3)
#' )
#' allocation <- n_alloc(frame, n = 1200)
#' summary(allocation)
#' summary(prec_alloc(allocation, n = allocation$detail$n_int))
#'
#' @seealso [n_alloc()] and [prec_alloc()] for classic allocation,
#'   [n_alloc_joint()] and [prec_alloc_joint()] for joint allocation,
#'   [strata_bound()] for constructing the strata to allocate over, and
#'   [design_df()] for the degrees of freedom the allocation leaves.
#'
#' @name summary.svyplan_alloc
NULL

#' @rdname summary.svyplan_alloc
#' @export
summary.svyplan_prec <- function(object, ...) {
  if (identical(object$type, "alloc_joint") ||
      (identical(object$type, "alloc") && identical(object$method, "bethel"))) {
    .check_unused_dots(...)
    return(.summary_bethel_alloc(object, kind = "assessment"))
  }
  if (identical(object$type, "alloc")) {
    .check_unused_dots(...)
    return(.summary_classic_alloc(object, kind = "assessment"))
  }
  .summary_as_print(object, ...)
}

#' Stages of the planned design, not of the solve
#'
#' A register fit's continuous relaxation is solved as one stage, its
#' clustering handed over as a design effect, but the plan it describes has
#' two.
#' @keywords internal
#' @noRd
.bethel_plan_stages <- function(object, problem) {
  if (is.null(object$params$psu)) problem$stages else 2L
}

#' Build a structured diagnostic summary for a generalized allocation
#' @keywords internal
#' @noRd
.summary_bethel_alloc <- function(object, kind) {
  p <- object$params
  problem <- p$problem
  budget <- p[["budget"]]
  budget_residual <- p[["budget_residual"]]
  if (is.null(problem) || is.null(problem$stratum_ids)) {
    stop(
      "generalized allocation result is missing its planning problem",
      call. = FALSE
    )
  }
  design <- identical(kind, "design")
  constraints <- if (design) object$constraints else object$detail
  operational_constraints <- if (design) {
    object$operational$constraints
  } else {
    NULL
  }
  mode <- if (design) {
    p$mode
  } else if (!is.null(object$objective) && !is.null(budget)) {
    "budget_objective"
  } else if (!is.null(object$objective)) {
    "objective"
  } else {
    "targets"
  }
  allocation <- .bethel_summary_allocation(object, design)
  bounds <- .bethel_summary_bounds(object, allocation, design)
  objective <- object$objective
  operational_objective <- if (design) object$operational$objective else NULL

  primary_constraints <- if (design) operational_constraints else constraints
  primary_objective <- if (design) operational_objective else objective
  floors <- object$domain_floors
  operational_floors <- if (design) object$operational$domain_floors
  overall <- .bethel_summary_overall(
    n = if (design) object$operational$n else p$achieved$n,
    cost = if (design) object$operational$cost else p$achieved$cost,
    constraints = primary_constraints,
    floors = if (design) operational_floors else floors,
    bounds = bounds,
    objective = primary_objective,
    objective_value = if (design) {
      object$operational$objective_value
    } else {
      object$objective_value
    },
    budget = budget,
    budget_residual = if (design) {
      object$operational$budget_residual
    } else {
      budget_residual
    }
  )
  continuous <- if (design) {
    .bethel_summary_overall(
      n = object$n,
      cost = p$achieved$cost,
      constraints = constraints,
      floors = floors,
      bounds = bounds,
      objective = objective,
      objective_value = object$objective_value,
      budget = budget,
      budget_residual = object$optimization$budget_residual
    )
  } else {
    NULL
  }

  # `certainty` is attached only when there is one, so a plan with no register
  # keeps the schema it has always had.
  out <- list(
    kind = kind,
    question = .bethel_summary_question(mode, kind, budget,
                                        p$objective_rule %||% "sum"),
    mode = mode,
    stages = .bethel_plan_stages(object, problem),
    status = if (design) object$optimization$classification else "assessed",
    overall = overall,
    continuous = continuous,
    allocation = allocation,
    constraints = constraints,
    operational_constraints = operational_constraints,
    objective = objective,
    operational_objective = operational_objective,
    bounds = bounds,
    optimization = if (design) {
      .bethel_summary_optimization(object)
    } else {
      NULL
    },
    assumptions = .bethel_summary_assumptions(object)
  )
  if (design) {
    note <- .bethel_sensitivity_note(
      object$optimization$slopes,
      if (identical(mode, "budget_objective")) budget,
      p$objective_rule %||% "sum"
    )
    if (length(note) > 0L) out$sensitivity_note <- note
  }
  if (!is.null(floors)) {
    out$domain_floors <- floors
    out$operational_domain_floors <- operational_floors
  }
  cert <- .bethel_summary_certainty(object)
  if (!is.null(cert)) {
    out$certainty <- cert
    sources <- c("threshold", "operational", "supplied", "planned")
    out$certainty_sources <- vapply(
      sources,
      function(src) sum(object$psu$.certainty_source %in% src),
      integer(1)
    )
  }
  structure(out, class = c("summary.svyplan_bethel", "list"))
}

#' The certainty split by stratum, which the printed line points at
#'
#' Counts are PSUs available in each part. The threshold and the design effect
#' are read from the plan as returned, so the table describes the answer.
#' @keywords internal
#' @noRd
.bethel_summary_certainty <- function(object) {
  if (is.null(object$params$psu) || is.null(object$psu)) {
    return(NULL)
  }
  detail <- object$detail
  out <- data.frame(
    stratum = detail$stratum,
    certainty = detail$n_psu_certain,
    remainder = detail$n_psu_rest,
    threshold = detail$threshold,
    stringsAsFactors = FALSE
  )
  # The threshold is the cutoff's, so a relaxed rule is shown beside it.
  cutoff <- .psu_cutoff(object$params$certainty_cutoff, object$params$frame)
  if (any(cutoff < 1)) out$cutoff <- cutoff
  if (!is.null(detail$n_zone)) out$zones <- detail$n_zone
  # The design effects the split produced already have a home, one row per
  # stratum and constraint in the resolved planning inputs below, so they are
  # not repeated here.
  out
}

#' Headline values for one numerical version of a generalized allocation
#' @keywords internal
#' @noRd
.bethel_summary_overall <- function(
  n,
  cost,
  constraints,
  bounds,
  floors = NULL,
  objective,
  objective_value,
  budget = NULL,
  budget_residual = NULL
) {
  n_constraints <- if (is.null(constraints)) 0L else nrow(constraints)
  n_pass <- if (n_constraints == 0L) 0L else sum(constraints$.pass)
  n_bounds_bad <- if (is.null(bounds)) {
    0L
  } else {
    status <- if ("field_status" %in% names(bounds)) {
      bounds$field_status
    } else {
      bounds$status
    }
    sum(grepl("violation", status), na.rm = TRUE)
  }
  list(
    n = n,
    cost = cost,
    n_constraints = n_constraints,
    n_pass = n_pass,
    n_violated = n_constraints - n_pass,
    all_pass = (n_constraints == 0L || n_pass == n_constraints) &&
      (is.null(floors) || all(floors$.pass)),
    n_floors = if (is.null(floors)) 0L else nrow(floors),
    n_floors_pass = if (is.null(floors)) 0L else sum(floors$.pass),
    n_bound_violations = n_bounds_bad,
    objective_value = objective_value,
    objective_components = if (is.null(objective)) 0L else nrow(objective),
    objective_worst = !is.null(objective) && ".binding" %in% names(objective),
    budget = budget,
    budget_residual = budget_residual
  )
}

#' Say what each constraint sensitivity holds fixed in a budget fit
#' @keywords internal
#' @noRd
.bethel_sensitivity_note <- function(status, budget, rule = "sum") {
  if (is.null(status)) {
    return(NULL)
  }
  budget_mode <- !is.null(budget)
  held <- if (budget_mode) {
    sprintf(
      paste(
        "Sens. is the change in minimum cost per unit relaxation of a target,",
        "the objective held fixed. Obj. sens. is the change in the %s",
        "relative variance per unit relaxation, the budget held at %.6g."
      ),
      if (identical(rule, "worst")) "worst weighted" else "weighted",
      budget
    )
  }
  refit <- "Refit with the changed target to measure its effect."
  slopes <- if (budget_mode) "The slopes are" else "Sens. is"
  reason <- switch(
    status,
    available = NULL,
    budget_not_binding = paste(
      "The budget does not bind, so relaxing a target leaves the objective",
      "unchanged."
    ),
    budget_at_minimum = paste(
      "The slopes are not defined for active targets because the budget",
      "equals the cheapest cost that meets them.", refit
    ),
    objective_multiplier_not_positive = paste(
      "Obj. sens. is not defined for active targets because the objective",
      "has no positive price at this budget.", refit
    ),
    not_identifiable = paste(
      "Obj. sens. is not defined for active targets because the objective",
      "shares its price with them.", refit
    ),
    register = paste(
      slopes, "not reported for a PSU register, whose prices hold for one",
      "certainty configuration only.", refit
    ),
    not_optimal = paste(
      slopes, "not reported because the continuous design is feasible but",
      "not certified optimal.", refit
    ),
    multipliers_not_verified = paste(
      slopes, "not reported because no multipliers price the full",
      "free-take problem at this design.", refit
    ),
    NULL
  )
  c(held, reason)
}

#' State which generalized allocation question was answered
#' @keywords internal
#' @noRd
.bethel_summary_question <- function(mode, kind, budget, rule = "sum") {
  if (identical(kind, "assessment")) {
    if (identical(rule, "worst")) {
      return(
        "assess a supplied joint allocation's precision and worst weighted objective"
      )
    }
    if (identical(mode, "budget_objective")) {
      return(
        "assess a supplied joint allocation against its targets, objective, and budget"
      )
    }
    if (identical(mode, "objective")) {
      return(
        "assess a supplied joint allocation's precision and weighted objective"
      )
    }
    return("assess a supplied joint allocation against every precision target")
  }
  if (identical(mode, "budget_objective") && identical(rule, "worst")) {
    return(sprintf(
      "minimize the worst weighted relative variance within budget = %.6g",
      budget
    ))
  }
  if (identical(mode, "budget_objective")) {
    return(sprintf(
      "find the best joint allocation affordable within budget = %.6g",
      budget
    ))
  }
  "find the cheapest joint allocation meeting every precision target"
}

#' Per-stratum decisions on their public ultimate-unit scale
#' @keywords internal
#' @noRd
.bethel_summary_allocation <- function(object, design) {
  p <- object$params
  problem <- p$problem
  take <- problem$stage$take
  if (design) {
    d <- object$detail
    n_continuous <- d$n
    n_primary <- d$n_int
  } else {
    n_continuous <- NULL
    n_primary <- as.numeric(p$n)
  }
  decision <- n_primary / take
  out <- data.frame(
    stratum = problem$stratum_ids,
    population = problem$population_N,
    stringsAsFactors = FALSE
  )
  if (design) {
    out$n_continuous <- n_continuous
  }
  out[[if (design) "n_field" else "n_supplied"]] <- n_primary
  out$weight <- problem$population_N / n_primary
  out$cost <- problem$cost * decision
  # A register's solver prices interviews at the marginal cost of its held
  # classification, so its certainty visits are added here, and the field
  # design is priced as it is drawn.
  if (!is.null(p$psu) && all(c("cost_1", "cost_2") %in% names(p$frame))) {
    cost_1 <- rep_len(as.numeric(p$frame$cost_1), nrow(out))
    cost_2 <- rep_len(as.numeric(p$frame$cost_2), nrow(out))
    held <- if (design) object$detail$n_psu_certain else
      p$certainty$n_psu_certain
    draw <- if (design) object$detail$n_psu_draw else p$certainty$n_psu_draw
    # A supplied allocation has no draw, so only its certainty visits are
    # added to the solver's price.
    out$cost <- if (is.null(draw)) {
      held * cost_1 + out$cost
    } else {
      (held + draw) * cost_1 + n_primary * cost_2
    }
  }
  if (problem$stages > 1L) {
    S <- problem$stages
    for (s in seq_len(S - 1L)) {
      out[[paste0("N_", s)]] <- problem$stage$N_stages[, s]
    }
    if (design) {
      out$n_1_continuous <- if ("n_1" %in% names(d)) d$n_1 else n_continuous / take
    }
    out[[if (design) "n_1_field" else "n_1_supplied"]] <-
      if (design && "n_1_int" %in% names(d)) d$n_1_int else decision
    free <- design && isTRUE(p$free_take) && "take_2" %in% names(d)
    for (s in 2L:S) {
      tc <- paste0("take_", s)
      if (free) {
        out[[paste0(tc, "_continuous")]] <- d[[tc]]
        out[[paste0(tc, "_field")]] <- problem$stage$takes[, s - 1L]
      } else {
        out[[tc]] <- problem$stage$takes[, s - 1L]
      }
    }
  }
  out
}

#' Bound positions and violations on the ultimate-unit scale
#' @keywords internal
#' @noRd
.bethel_summary_bounds <- function(object, allocation, design) {
  p <- object$params
  problem <- p$problem
  lower <- problem$lower * problem$stage$take
  upper <- problem$upper * problem$stage$take
  if (design) {
    continuous <- object$detail$n
    field <- object$detail$n_int
    data.frame(
      stratum = problem$stratum_ids,
      lower = lower,
      continuous = continuous,
      field = field,
      upper = upper,
      continuous_status = .bethel_bound_status(continuous, lower, upper),
      field_status = .bethel_bound_status(field, lower, upper),
      stringsAsFactors = FALSE
    )
  } else {
    supplied <- allocation$n_supplied
    data.frame(
      stratum = problem$stratum_ids,
      lower = lower,
      supplied = supplied,
      upper = upper,
      lower_violation = supplied < lower - p$feasibility_tolerance,
      upper_violation = supplied > upper + p$feasibility_tolerance,
      status = .bethel_bound_status(supplied, lower, upper),
      pass = !(supplied < lower - p$feasibility_tolerance |
        supplied > upper + p$feasibility_tolerance),
      stringsAsFactors = FALSE
    )
  }
}

#' Label one allocation's relation to each pair of bounds
#' @keywords internal
#' @noRd
.bethel_bound_status <- function(value, lower, upper, tolerance = 1e-6) {
  tol <- tolerance * pmax(1, abs(value), abs(lower), abs(upper))
  out <- rep("interior", length(value))
  out[value < lower - tol] <- "lower violation"
  out[value > upper + tol] <- "upper violation"
  fixed <- abs(upper - lower) <= tol
  out[fixed & abs(value - lower) <= tol] <- "fixed"
  out[!fixed & abs(value - lower) <= tol] <- "lower"
  out[!fixed & abs(value - upper) <= tol] <- "upper"
  out
}

#' Make solver diagnostics readable without discarding their raw values
#' @keywords internal
#' @noRd
.bethel_summary_optimization <- function(object) {
  opt <- object$optimization
  ids <- object$params$problem$constraint_ids
  solver_ids <- .bethel_solver_ids(
    object, length(opt$constraint_scale %||% numeric(0))
  )
  active <- opt$active_precision %||% integer(0)
  active_precision <- .bethel_solver_ids(object, max(c(0, active)))[active]
  diagnostic_names <- c(
    "primal_residual",
    "projected_dual_residual",
    "stationarity_residual",
    "complementarity_residual",
    "relative_duality_gap"
  )
  present <- diagnostic_names[diagnostic_names %in% names(opt)]
  diagnostics <- data.frame(
    diagnostic = present,
    value = unname(unlist(opt[present])),
    tolerance = rep(opt$kkt_tolerance %||% NA_real_, length(present)),
    stringsAsFactors = FALSE
  )
  diagnostics$pass <- diagnostics$value <= diagnostics$tolerance
  multiplier_ids <- .bethel_solver_ids(
    object, length(opt$multiplier_identifiable %||% logical(0))
  )
  identifiable <- opt$multiplier_identifiable %||% logical(0)
  if (length(identifiable) > 0L) {
    names(identifiable) <- multiplier_ids
  }
  constraint_scale <- opt$constraint_scale %||% numeric(0)
  constraint_residual <- opt$constraint_residual %||% numeric(0)
  if (length(constraint_scale) > 0L) {
    names(constraint_scale) <- solver_ids
  }
  if (length(constraint_residual) > 0L) {
    names(constraint_residual) <- solver_ids
  }
  list(
    classification = opt$classification,
    converged = opt$converged,
    search = opt$search,
    message = opt$message,
    iterations = opt$iterations,
    convergence_code = opt$convergence_code,
    feasibility_tolerance = opt$feasibility_tolerance,
    kkt_tolerance = opt$kkt_tolerance,
    diagnostics = diagnostics,
    active_precision = unname(active_precision),
    active_lower = object$params$problem$stratum_ids[
      opt$active_lower %||% integer(0)
    ],
    active_upper = object$params$problem$stratum_ids[
      opt$active_upper %||% integer(0)
    ],
    multiplier_identifiable = identifiable,
    continuous_cost = opt$cost,
    cost_scale = opt$cost_scale,
    constraint_scale = constraint_scale,
    constraint_residual = constraint_residual,
    dual_value = opt$dual_value,
    budget = opt$budget,
    budget_binding = opt$budget_binding,
    budget_residual = opt$budget_residual,
    budget_sensitivity = opt$budget_sensitivity,
    objective_bound = opt$objective_bound,
    objective_multiplier = opt$objective_multiplier,
    root_iterations = opt$root_iterations,
    objective_weights = opt$objective_weights,
    second_stage = opt$second_stage,
    slopes = opt$slopes
  )
}

#' Names of the solver rows: the targets, then the objective or its caps
#'
#' A sum-rule budget solve appends one objective row. A worst-rule solve
#' appends one cap per component with a variance term.
#' @keywords internal
#' @noRd
.bethel_solver_ids <- function(object, n) {
  ids <- object$params$problem$constraint_ids
  if (n <= length(ids)) {
    return(ids[seq_len(n)])
  }
  obj <- object$params$problem$objective
  extra <- if (identical(obj$rule, "worst")) {
    paste0("worst: ", obj$component[colSums(obj$A) > 0])
  } else {
    "objective"
  }
  c(ids, extra)[seq_len(n)]
}

#' Label a component by its name, or its full key inside a domain
#' @keywords internal
#' @noRd
.bethel_component_label <- function(tab) {
  ifelse(tab$domain == ".overall", tab$name, tab$component)
}

#' Resolved numerical assumptions behind each target-stratum contribution
#' @keywords internal
#' @noRd
.bethel_summary_assumptions <- function(object) {
  p <- object$params
  problem <- p$problem
  K <- length(problem$constraint_ids)
  H <- length(problem$stratum_ids)
  model <- NULL
  if (K > 0L) {
    grid <- expand.grid(
      stratum_index = seq_len(H),
      constraint_index = seq_len(K),
      KEEP.OUT.ATTRS = FALSE
    )
    keep <- problem$membership[cbind(grid$stratum_index, grid$constraint_index)]
    grid <- grid[keep, , drop = FALSE]
    at <- cbind(grid$stratum_index, grid$constraint_index)
    model <- data.frame(
      stratum = problem$stratum_ids[grid$stratum_index],
      constraint = problem$constraint_ids[grid$constraint_index],
      mean = problem$mean_hk[at],
      variance = problem$var_hk[at],
      deff = problem$deff_hk[at],
      response_rate = problem$resp_hk[at],
      stringsAsFactors = FALSE
    )
    scope <- problem$scope_hk
    if (!is.null(scope) && any(scope[at] != 1)) {
      model$eligibility_rate <- scope[at]
    }
    elig <- problem$elig_hk
    if (!is.null(elig) && any(elig[at] != 1)) {
      model$eligible_per_unit <- elig[at]
    }
  }
  list(
    alpha = p$alpha,
    df = problem$df,
    stages = .bethel_plan_stages(object, problem),
    # The solver saw a one-stage problem because a register hands it the
    # clustering as a design effect, so the model it recorded describes the
    # solve rather than the plan.
    variance_model = if (is.null(p$psu)) {
      problem$variance_model
    } else {
      "two_part_certainty_wald"
    },
    feasibility_tolerance = p$feasibility_tolerance,
    fpc = p$fpc %||% "unit",
    min_n_stratum = p$min_n_stratum,
    free_take = isTRUE(p$free_take),
    fixed_takes = if (problem$stages > 1L) {
      ft <- data.frame(
        stratum = problem$stratum_ids,
        stringsAsFactors = FALSE
      )
      for (s in 2L:problem$stages) {
        ft[[paste0("take_", s)]] <- problem$stage$takes[, s - 1L]
      }
      ft
    } else {
      NULL
    },
    model = model
  )
}

#' Build one structured summary for design and assessment directions
#' @keywords internal
#' @noRd
.summary_classic_alloc <- function(object, kind) {
  detail <- object$detail
  if (is.null(detail) || !all(c("stratum", "N", "n") %in% names(detail))) {
    stop(
      "classic allocation result is missing its per-stratum detail",
      call. = FALSE
    )
  }

  design <- identical(kind, "design")
  n_used <- if (design) detail$n_int else detail$n
  state <- .alloc_summary_state(object, n_used, operational = design)
  overall <- .alloc_summary_overall(state, n_used)
  if (design) {
    df <- tryCatch(
      suppressWarnings(as.double(design_df(object))),
      error = function(e) NA_real_
    )
    overall$design_df <- df
  }

  continuous <- NULL
  if (design) {
    continuous <- list(
      n = object$n,
      expected_respondents = sum(detail$n * state$response_rate),
      cost = if (all(is.na(state$cost_h))) {
        NA_real_
      } else {
        object$params$achieved$cost
      },
      se = object$se,
      moe = object$moe,
      rmoe = object$rmoe,
      cv = object$cv
    )
  }

  domains <- .alloc_summary_domains(object, state, n_used)
  if (!is.null(domains)) {
    overall$worst_domain_cv <- max(domains$cv)
    if (design && !is.null(object$domains)) {
      continuous$worst_domain_cv <- max(object$domains$.cv)
    }
  }

  mode <- if (design) object$params$mode else "assessment"
  method <- if (design) object$method else NULL
  out <- structure(
    list(
      kind = kind,
      question = .alloc_summary_question(object, kind),
      method = method,
      mode = mode,
      overall = overall,
      continuous = continuous,
      allocation = .alloc_summary_allocation(
        object,
        state,
        n_used,
        kind = kind
      ),
      precision = if (!.alloc_size_only(object)) {
        .alloc_summary_precision(state, n_used)
      },
      bounds = if (design) .alloc_summary_bounds(detail, n_used) else NULL,
      domains = domains,
      assumptions = list(
        alpha = object$params$alpha,
        deff = state$deff,
        response_rate = state$response_rate / state$eligibility_rate,
        eligibility_rate = if (any(state$scope != 1)) state$scope,
        eligible_per_unit = if (any(state$per_unit != 1)) state$per_unit,
        unit_response_rate = state$unit_response_rate,
        psu_response_rate = state$psu_response_rate,
        df = object$params$df,
        fpc = object$params$fpc %||% "unit",
        min_n_stratum = object$params$min_n_stratum,
        alloc_q = object$params$alloc_q
      )
    ),
    class = c("summary.svyplan_alloc", "list")
  )
  if (design && !is.null(object$optimization)) out$search <- object$optimization
  out
}

#' Resolve the variance, response, and cost model at the displayed allocation
#' @keywords internal
#' @noRd
.alloc_summary_state <- function(object, n_h, operational) {
  p <- object$params
  d <- object$detail
  frame <- p$frame
  H <- nrow(d)
  deff <- .alloc_resolve_h(p$deff, frame, "deff", H, .check_deff_h)
  rates <- .alloc_summary_response_rates(object)
  unit_response <- rates$unit
  cluster <- .alloc_is_cluster(frame)
  psu_response <- rates$psu
  response <- rates$combined
  S_h <- d$sd
  cost_h <- d$unit_cost
  N_scope <- d$N * rates$eligibility
  N_fpc <- if (identical(p$fpc, "none")) rep(Inf, H) else N_scope

  if (cluster) {
    S <- .alloc_detect_stages(frame)
    takes_mat <- matrix(0, nrow = H, ncol = S - 1L)
    for (s in 2L:S) {
      col_int <- paste0("take_", s, "_int")
      col_cont <- paste0("take_", s)
      takes_mat[, s - 1L] <- if (operational && col_int %in% names(d)) {
        d[[col_int]]
      } else {
        d[[col_cont]]
      }
    }
    icc_mat <- matrix(0, nrow = H, ncol = S - 1L)
    vr_mat <- matrix(0, nrow = H, ncol = S - 1L)
    for (s in seq_len(S - 1L)) {
      icc_mat[, s] <- frame[[paste0("icc_", s)]]
      vr_col <- paste0("var_ratio_", s)
      if (vr_col %in% names(frame)) {
        vr_mat[, s] <- frame[[vr_col]]
      } else if (s == 1L) {
        vr_mat[, 1L] <- 1
      } else {
        vr_mat[, s] <- vr_mat[, s - 1L] * (1 - icc_mat[, s - 1L])
      }
    }
    fpc_mode <- p$fpc %||% "unit"
    rr_mat <- rates$stage_mat
    inflation <- numeric(H)
    for (h in seq_len(H)) {
      eff_t <- takes_mat[h, ] * rr_mat[h, 2L:S]
      Q <- numeric(S)
      Q[S] <- 1
      for (s in (S - 1L):1L) Q[s] <- Q[s + 1L] * eff_t[s]
      R_h <- numeric(S)
      for (s in seq_len(S - 1L)) R_h[s] <- vr_mat[h, s] * icc_mat[h, s]
      R_h[S] <- vr_mat[h, S - 1L] * (1 - icc_mat[h, S - 1L])
      f_stage <- numeric(S)
      if (identical(fpc_mode, "stage")) {
        for (s in 2L:S) {
          M_s <- if (s < S) {
            frame[[paste0("N_", s)]][h] / frame[[paste0("N_", s - 1L)]][h]
          } else {
            N_scope[h] / frame[[paste0("N_", S - 1L)]][h]
          }
          f_stage[s] <- min(1, eff_t[s - 1L] / M_s)
        }
      }
      inflation[h] <- R_h[1L] * Q[1L] +
        sum((1 - f_stage[2L:S]) * R_h[2L:S] * Q[2L:S])
    }
    R1_vec <- vr_mat[, 1L] * icc_mat[, 1L]
    S_h <- d$sd * sqrt(inflation)
    N_fpc <- switch(
      fpc_mode,
      unit = N_scope,
      none = rep(Inf, H),
      stage = ifelse(R1_vec > 0, frame$N_1 * inflation / R1_vec, Inf)
    )
    cost_cols <- paste0("cost_", seq_len(S))
    if (all(cost_cols %in% names(frame))) {
      cost_h <- numeric(H)
      for (h in seq_len(H)) {
        cum_t <- c(1, cumprod(takes_mat[h, ]))
        c_vec <- vapply(cost_cols, function(col) frame[[col]][h], numeric(1L))
        cost_h[h] <- sum(c_vec * cum_t) / prod(takes_mat[h, ])
      }
    } else {
      cost_h <- rep(NA_real_, H)
    }
  }

  mean_h <- if ("mean" %in% names(d)) d$mean else rep(NA_real_, H)
  metrics <- .alloc_metrics(
    N_h = d$N,
    N_pop = N_scope,
    S_h = S_h,
    mean_h = mean_h,
    n_h = n_h,
    alpha = p$alpha,
    deff = deff,
    resp_rate = response,
    cost_h = cost_h,
    df = p$df,
    N_fpc = N_fpc
  )
  list(
    metrics = metrics,
    stratum = d$stratum,
    population = d$N,
    population_in_scope = N_scope,
    N_fpc = N_fpc,
    sd = S_h,
    mean = mean_h,
    deff = deff,
    unit_response_rate = unit_response,
    psu_response_rate = psu_response,
    response_rate = response,
    eligibility_rate = rates$eligibility,
    scope = rates$scope,
    per_unit = rates$per_unit,
    cost_h = cost_h,
    cluster = cluster
  )
}

#' Resolve the stage rates whose product leaves responding ultimate units
#' @keywords internal
#' @noRd
.alloc_summary_response_rates <- function(object) {
  p <- object$params
  frame <- p$frame
  H <- nrow(object$detail)
  unit_raw <- .alloc_resolve_h(
    p$resp_rate,
    frame,
    "resp_rate",
    H,
    .check_resp_rate_h
  )
  S <- .alloc_detect_stages(frame)
  stage_mat <- matrix(1, nrow = H, ncol = S)
  if (S >= 2L) {
    for (s in seq_len(S - 1L)) {
      stage_mat[, s] <- .alloc_resolve_h(
        1, frame, paste0("resp_rate_", s), H, .check_resp_rate_h
      )
    }
  }
  psu <- if (S >= 2L) stage_mat[, 1L] else rep(1, H)
  eligibility <- .alloc_resolve_h(
    p[["eligibility_rate"]] %||% 1,
    frame,
    "eligibility_rate",
    H,
    .check_resp_rate_h
  )
  eligible_per_unit <- .alloc_resolve_h(
    p[["eligible_per_unit"]] %||% 1,
    frame,
    "eligible_per_unit",
    H,
    .check_deff_h
  )
  yield <- eligibility * eligible_per_unit
  stage_mat[, S] <- unit_raw * yield
  combined <- apply(stage_mat, 1L, prod)
  response <- combined / yield
  list(unit = unit_raw * yield, psu = psu, stage_mat = stage_mat,
       combined = combined, eligibility = yield, response = response,
       scope = eligibility, per_unit = eligible_per_unit)
}

#' Overall fields on the same numerical basis as the stratum tables
#' @keywords internal
#' @noRd
.alloc_summary_overall <- function(state, n_h) {
  m <- state$metrics
  list(
    n = sum(n_h),
    expected_respondents = sum(n_h * state$response_rate),
    cost = m$cost,
    se = m$se,
    moe = m$moe,
    rmoe = m$rmoe,
    cv = m$cv
  )
}

#' State the allocation question rather than only its algorithm
#' @keywords internal
#' @noRd
.alloc_summary_question <- function(object, kind) {
  if (identical(kind, "assessment")) {
    return("assess a supplied stratified allocation")
  }
  p <- object$params
  method <- paste0(
    toupper(substr(object$method, 1L, 1L)),
    substring(object$method, 2L)
  )
  switch(
    p$mode,
    n = sprintf(
      "distribute a fixed sample of %g using %s allocation",
      p$n,
      method
    ),
    cv = sprintf(
      "find the smallest %s allocation attaining cv = %.4g",
      method,
      p$cv
    ),
    moe = sprintf(
      "find the smallest %s allocation attaining moe = %.4g",
      method,
      p$moe
    ),
    rmoe = sprintf(
      "find the smallest %s allocation attaining rmoe = %.4g",
      method,
      p$rmoe
    ),
    budget = sprintf(
      "find the best %s allocation within budget = %.6g",
      method,
      p$budget
    ),
    sprintf("construct a %s allocation", method)
  )
}

#' Field or supplied allocation table
#' @keywords internal
#' @noRd
.alloc_summary_allocation <- function(object, state, n_h, kind) {
  d <- object$detail
  out <- data.frame(
    stratum = d$stratum,
    population = d$N,
    response_rate = state$response_rate / state$eligibility_rate,
    unit_cost = state$cost_h,
    stringsAsFactors = FALSE
  )
  extra <- list(
    eligibility_rate = if (any(state$scope != 1)) state$scope,
    eligible_per_unit = if (any(state$per_unit != 1)) state$per_unit
  )
  extra <- Filter(Negate(is.null), extra)
  if (length(extra) > 0L) {
    out <- cbind(out[1:3], as.data.frame(extra), out[-(1:3)])
  }
  if (identical(kind, "design")) {
    out$n_continuous <- d$n
    out$n_field <- n_h
  } else {
    out$n_supplied <- n_h
  }
  if (state$cluster) {
    take_cols <- grep("^take_[0-9]+$", names(d), value = TRUE)
    if (identical(kind, "design")) {
      out$n_1_continuous <- d$n_1
      out$n_1_field <- d$n_1_int
      for (tc in take_cols) {
        out[[paste0(tc, "_continuous")]] <- d[[tc]]
        tc_int <- paste0(tc, "_int")
        if (tc_int %in% names(d)) {
          out[[paste0(tc, "_field")]] <- d[[tc_int]]
        }
      }
    } else {
      out$n_1 <- d$n_1
      for (tc in take_cols) {
        out[[tc]] <- d[[tc]]
      }
    }
  }
  out$expected_respondents <- n_h * state$response_rate
  out$weight <- d$N / n_h
  out$cost <- n_h * state$cost_h
  if (all(is.na(out$unit_cost))) {
    out$unit_cost <- NULL
  }
  if (all(is.na(out$cost))) {
    out$cost <- NULL
  }
  out
}

#' Per-stratum precision at the field or supplied allocation
#' @keywords internal
#' @noRd
.alloc_summary_precision <- function(state, n_h) {
  m <- state$metrics
  data.frame(
    stratum = state$stratum,
    effective_n = n_h * state$response_rate / state$deff,
    se = m$se_h,
    moe = m$moe_h,
    rmoe = m$rmoe_h,
    cv = m$cv_h,
    variance_share = m$share_h,
    stringsAsFactors = FALSE
  )
}

#' Bound table for a fitted classic allocation
#' @keywords internal
#' @noRd
.alloc_summary_bounds <- function(detail, n_h) {
  data.frame(
    stratum = detail$stratum,
    lower = detail$.lower,
    field = n_h,
    upper = detail$.upper,
    binding = detail$.binding,
    source = detail$.bound_source,
    stringsAsFactors = FALSE
  )
}

#' Recompute domains on the same basis as the headline
#' @keywords internal
#' @noRd
.alloc_summary_domains <- function(object, state, n_h) {
  p <- object$params
  idx_list <- p$domain_idx
  if (is.null(idx_list) || length(idx_list) == 0L) {
    return(NULL)
  }
  frame <- p$frame
  ids <- p$domain_cols %||% character(0)
  rows <- lapply(seq_along(idx_list), function(i) {
    idx <- idx_list[[i]]
    m <- .alloc_metrics(
      N_h = state$population[idx],
      N_pop = state$population_in_scope[idx],
      S_h = state$sd[idx],
      mean_h = state$mean[idx],
      n_h = n_h[idx],
      alpha = p$alpha,
      deff = state$deff[idx],
      resp_rate = state$response_rate[idx],
      cost_h = state$cost_h[idx],
      df = p$df,
      N_fpc = state$N_fpc[idx]
    )
    id <- if (length(ids) > 0L) {
      frame[idx[1L], ids, drop = FALSE]
    } else {
      data.frame()
    }
    id$.domain <- names(idx_list)[i]
    id$population <- sum(state$population[idx])
    id$n <- sum(n_h[idx])
    id$expected_respondents <- sum(n_h[idx] * state$response_rate[idx])
    id$se <- m$se
    id$moe <- m$moe
    id$rmoe <- m$rmoe
    id$cv <- m$cv
    id$cost <- m$cost
    id
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (all(is.na(out$cost))) {
    out$cost <- NULL
  }
  out
}

#' @rdname summary.svyplan_alloc
#' @export
print.summary.svyplan_bethel <- function(x, ...) {
  .check_unused_dots(...)
  design <- identical(x$kind, "design")
  cat(
    if (design) {
      "Generalized Bethel allocation summary\n\n"
    } else {
      "Generalized Bethel precision summary\n\n"
    }
  )
  cat(sprintf("Question: %s\n", x$question))
  cat(sprintf("Status: %s\n", x$status))
  cat(sprintf("Stages: %d\n\n", x$stages))

  .print_bethel_overall(
    x$overall,
    if (design) "Operational field design" else "Supplied allocation"
  )
  if (!is.null(x$continuous)) {
    diag <- x$optimization$diagnostics
    gap_note <- .bethel_gap_note(
      x$status,
      diag$value[diag$diagnostic %in% "relative_duality_gap"]
    )
    cat(if (is.null(gap_note)) "\nContinuous optimum\n" else
      sprintf("\nContinuous design, feasible, %s\n", gap_note))
    .print_bethel_overall(x$continuous, title = NULL)
  }

  cat("\nAllocation by stratum\n\n")
  main <- intersect(
    c(
      "stratum",
      "population",
      "n_continuous",
      "n_field",
      "n_supplied",
      "weight",
      "cost"
    ),
    names(x$allocation)
  )
  .print_bethel_table(x$allocation[main])
  stage_pat <- paste0(
    "^(N_[0-9]+|n_1_(continuous|field|supplied)|",
    "take_[0-9]+(_continuous|_field)?)$"
  )
  stages <- c("stratum", grep(stage_pat, names(x$allocation), value = TRUE))
  if (length(stages) > 1L) {
    cat("\nStage decisions by stratum\n\n")
    .print_bethel_table(x$allocation[stages])
  }

  .print_bethel_constraints(
    x$constraints,
    if (design) {
      "Continuous precision constraints"
    } else {
      "Assessed precision constraints"
    },
    note = x$sensitivity_note
  )
  if (design) {
    .print_bethel_constraints(
      x$operational_constraints,
      "Operational precision constraints"
    )
  }

  .print_bethel_floors(
    x$domain_floors,
    if (design) "Continuous sample floors" else "Assessed sample floors"
  )
  if (design) {
    .print_bethel_floors(x$operational_domain_floors,
                         "Operational sample floors")
  }

  .print_bethel_objective(
    x$objective,
    if (design) {
      "Continuous objective components"
    } else {
      "Assessed objective components"
    }
  )
  if (design) {
    .print_bethel_objective(
      x$operational_objective,
      "Operational objective components"
    )
  }

  cat("\nAllocation bounds\n\n")
  .print_bethel_table(x$bounds)
  if (x$overall$n_bound_violations == 0L) {
    cat("\nNo allocation bounds are violated.\n")
  }

  if (!is.null(x$optimization)) {
    .print_bethel_optimization(x$optimization)
  }

  cat("\nAssumptions\n")
  cat(sprintf("  variance model: %s\n", x$assumptions$variance_model))
  cat(sprintf("  alpha: %s\n", .fmt_prob(x$assumptions$alpha)))
  cat(sprintf(
    "  feasibility tolerance: %.3g\n",
    x$assumptions$feasibility_tolerance
  ))
  .print_alloc_fpc(x$assumptions$fpc)
  if (!is.null(x$assumptions$df) && any(is.finite(x$assumptions$df))) {
    cat(sprintf(
      "  interval df: %s\n",
      paste(
        unique(x$assumptions$df[is.finite(x$assumptions$df)]),
        collapse = ", "
      )
    ))
  }
  if (!is.null(x$assumptions$fixed_takes)) {
    cat(if (isTRUE(x$assumptions$free_take)) {
      "\nOptimized later-stage takes\n\n"
    } else {
      "\nFixed later-stage takes\n\n"
    })
    .print_bethel_table(x$assumptions$fixed_takes)
  }
  if (!is.null(x$certainty)) {
    cat("\nCertainty split by stratum\n")
    by_source <- x$certainty_sources[x$certainty_sources > 0L]
    if (length(by_source)) {
      cat(sprintf(
        "certainty PSUs by source: %s\n",
        paste(names(by_source), by_source, collapse = ", ")
      ))
    }
    tab <- x$certainty
    if (nrow(tab) > 0L) {
      cat("\n")
      if (is.numeric(tab$threshold)) tab$threshold <- round(tab$threshold)
      .print_bethel_table(tab)
    }
  }
  if (!is.null(x$assumptions$model) && nrow(x$assumptions$model) > 0L) {
    cat("\nResolved target-stratum planning inputs\n\n")
    .print_bethel_table(x$assumptions$model)
  }
  invisible(x)
}

#' Print generalized allocation headline values
#' @keywords internal
#' @noRd
.print_bethel_overall <- function(values, title) {
  if (!is.null(title)) {
    cat(title, "\n", sep = "")
  }
  rows <- list(
    Sample = .fmt_continuous_n(values$n),
    Cost = sprintf("%.2f", values$cost)
  )
  rows[["Hard targets"]] <- if (values$n_constraints == 0L) {
    "none"
  } else {
    sprintf(
      "%d of %d pass%s",
      values$n_pass,
      values$n_constraints,
      if (values$n_violated == 0L) {
        ""
      } else {
        sprintf(", %d violated", values$n_violated)
      }
    )
  }
  if (isTRUE(values$n_floors > 0L)) {
    rows[["Sample floors"]] <- sprintf("%d of %d met", values$n_floors_pass,
                                       values$n_floors)
  }
  if (!is.null(values$objective_value)) {
    rows[[if (isTRUE(values$objective_worst)) {
      "Worst weighted rel. var."
    } else {
      "Weighted rel. variance"
    }]] <- sprintf(
      "%.6g (%d component%s)",
      values$objective_value,
      values$objective_components,
      if (values$objective_components == 1L) "" else "s"
    )
  }
  if (!is.null(values$budget)) {
    rows[["Budget"]] <- sprintf("%.6g", values$budget)
    if (!is.null(values$budget_residual)) {
      rows[["Residual"]] <- sprintf("%.6g", values$budget_residual)
    }
  }
  width <- max(nchar(names(rows)))
  for (name in names(rows)) {
    cat(sprintf("  %-*s  %s\n", width, name, rows[[name]]))
  }
  invisible(NULL)
}

#' Print the regional sample floors, required against achieved
#' @keywords internal
#' @noRd
.print_bethel_floors <- function(tab, title) {
  if (is.null(tab) || nrow(tab) == 0L) return(invisible(NULL))
  cat("\n", title, "\n\n", sep = "")
  .print_bethel_table(tab[c("domain", "level", "required", "achieved",
                            ".pass")])
  invisible(NULL)
}

#' Print both the result and numerical diagnostics for every hard target
#' @keywords internal
#' @noRd
.print_bethel_constraints <- function(tab, title, note = NULL) {
  cat("\n", title, "\n\n", sep = "")
  if (is.null(tab) || nrow(tab) == 0L) {
    cat("No hard precision constraints.\n")
    return(invisible(NULL))
  }
  result <- intersect(
    c(
      "constraint",
      ".metric",
      ".target",
      ".achieved",
      ".ratio",
      ".pass",
      ".binding"
    ),
    names(tab)
  )
  .print_bethel_table(tab[result])
  diagnostics <- intersect(
    c("constraint", ".residual", ".tolerance", ".multiplier", ".sensitivity",
      ".objective_sensitivity"),
    names(tab)
  )
  if (length(diagnostics) > 1L) {
    cat("\nConstraint diagnostics\n\n")
    .print_bethel_table(tab[diagnostics])
    if (length(note) > 0L) {
      cat("\n")
      cat(strwrap(paste(note, collapse = " "), width = 78), sep = "\n")
    }
  }
  invisible(NULL)
}

#' Print a genuinely additive generalized-allocation objective
#' @keywords internal
#' @noRd
.print_bethel_objective <- function(tab, title) {
  if (is.null(tab) || nrow(tab) == 0L) {
    return(invisible(NULL))
  }
  worst <- ".binding" %in% names(tab)
  cat("\n", title, if (worst) {
    " (the largest contribution is the objective)"
  } else {
    " (contributions are additive)"
  }, "\n\n", sep = "")
  cols <- intersect(
    c("component", "priority", ".relvar", ".cv", ".contribution", ".share",
      ".binding"),
    names(tab)
  )
  .print_bethel_table(tab[cols])
  invisible(NULL)
}

#' Print solver certification and local sensitivity
#' @keywords internal
#' @noRd
.print_bethel_optimization <- function(opt) {
  cat("\nOptimization diagnostics\n")
  cat(sprintf("  classification: %s\n", opt$classification))
  cat(sprintf("  converged: %s\n", if (isTRUE(opt$converged)) "yes" else "no"))
  if (!is.null(opt$iterations)) {
    cat(sprintf("  solver iterations: %d\n", opt$iterations))
  }
  if (length(opt$active_precision) > 0L) {
    cat(sprintf(
      "  active precision: %s\n",
      paste(opt$active_precision, collapse = ", ")
    ))
  }
  cat(sprintf(
    "  active bounds: %d lower, %d upper\n",
    length(opt$active_lower),
    length(opt$active_upper)
  ))
  if (!is.null(opt$search)) {
    sr <- opt$search
    cat(sprintf(
      "  %s: %s after %s node%s of %s, bound %s (gap %.3g%%)\n",
      if (identical(sr$family, "whole units")) "whole-unit search" else
        "field design search",
      sr$status,
      format(sr$nodes, big.mark = ","),
      if (sr$nodes == 1) "" else "s",
      format(sr$max_nodes, big.mark = ","),
      format(signif(sr$bound, 6), big.mark = ","),
      100 * sr$gap
    ))
    if (!is.na(sr$states %||% NA)) {
      cat(sprintf(
        "  states: %s listed, %s kept by the bound\n",
        format(sr$states, big.mark = ","),
        format(sr$kept, big.mark = ",")
      ))
    }
    if (!is.null(sr$trials)) {
      cat(sprintf(
        "  worst-level bisection: %d trial%s, second stage %s\n",
        sr$trials, if (sr$trials == 1L) "" else "s", sr$second_stage
      ))
    }
  }
  if (!is.null(opt$second_stage)) {
    cat(sprintf("  continuous second stage: %s\n", opt$second_stage))
  }
  if (!is.null(opt$objective_weights) && !anyNA(opt$objective_weights)) {
    cat("  sum-rule priorities giving the same continuous design:\n")
    cat(sprintf("    %s %.4g\n", names(opt$objective_weights),
                opt$objective_weights), sep = "")
  }
  if (!is.null(opt$budget_binding)) {
    cat(sprintf(
      "  budget binding: %s\n",
      if (isTRUE(opt$budget_binding)) "yes" else "no"
    ))
  }
  if (!is.null(opt$budget_sensitivity) && is.finite(opt$budget_sensitivity)) {
    cat(sprintf("  budget sensitivity: %.6g\n", opt$budget_sensitivity))
  } else if (isTRUE(opt$budget_binding) && !is.null(opt$slopes)) {
    cat(sprintf("  budget sensitivity: not reported (%s)\n", opt$slopes))
  }
  if (nrow(opt$diagnostics) > 0L) {
    cat("\nKKT certification\n\n")
    .print_bethel_table(opt$diagnostics)
  }
  invisible(NULL)
}

#' Format and print one bounded-width diagnostic table
#' @keywords internal
#' @noRd
.print_bethel_table <- function(tab, max_rows = 20L) {
  total <- nrow(tab)
  shown <- utils::head(tab, max_rows)
  print(
    .fmt_bethel_summary_table(shown),
    row.names = FALSE,
    quote = FALSE,
    right = TRUE
  )
  if (total > nrow(shown)) {
    cat(sprintf(
      "... %d of %d rows shown (all rows retained in the object)\n",
      nrow(shown),
      total
    ))
  }
  invisible(NULL)
}

#' Format generalized diagnostic tables without changing stored values
#' @keywords internal
#' @noRd
.fmt_bethel_summary_table <- function(tab) {
  out <- tab
  id_cols <- intersect(c("constraint", "component"), names(out))
  for (name in id_cols) {
    out[[name]] <- vapply(out[[name]], .abbrev_bethel_id, character(1L))
  }
  count_cols <- union(
    intersect(c("population", "n_field", "field"), names(out)),
    grep("^N_[0-9]+$", names(out), value = TRUE)
  )
  for (name in count_cols) {
    out[[name]] <- vapply(out[[name]], .fmt_count_n, character(1L))
  }
  two_cols <- intersect(
    c(
      "n_continuous",
      "n_supplied",
      "continuous",
      "supplied",
      "lower",
      "upper",
      "n_1_continuous",
      "n_1_supplied",
      "cost"
    ),
    names(out)
  )
  for (name in two_cols) {
    out[[name]] <- sprintf("%.2f", out[[name]])
  }
  if ("weight" %in% names(out)) {
    out$weight <- sprintf("%.3f", out$weight)
  }
  number_cols <- intersect(
    c(
      "mean",
      "variance",
      "deff",
      "response_rate",
      "priority",
      ".metric",
      ".target",
      ".achieved",
      ".ratio",
      ".residual",
      ".tolerance",
      ".multiplier",
      ".sensitivity",
      ".objective_sensitivity",
      ".relvar",
      ".cv",
      ".contribution",
      ".share",
      "value",
      "tolerance"
    ),
    names(out)
  )
  number_cols <- setdiff(number_cols, ".metric")
  for (name in number_cols) {
    out[[name]] <- ifelse(is.na(out[[name]]), "", sprintf("%.4g", out[[name]]))
  }
  logical_cols <- names(out)[vapply(out, is.logical, logical(1L))]
  for (name in logical_cols) {
    out[[name]] <- ifelse(
      is.na(out[[name]]),
      "",
      ifelse(out[[name]], "yes", "no")
    )
  }
  labels <- c(
    stratum = "Stratum",
    population = "Pop.",
    n_continuous = "Cont.",
    n_field = "Field",
    n_supplied = "Supplied",
    weight = "Weight",
    cost = "Cost",
    n_1_continuous = "n_1 cont.",
    n_1_field = "n_1 field",
    n_1_supplied = "n_1 supplied",
    constraint = "Constraint",
    component = "Component",
    .metric = "Metric",
    .target = "Target",
    .achieved = "Ach.",
    .ratio = "Ratio",
    .pass = "Pass",
    .binding = "Bind",
    .residual = "Residual",
    .tolerance = "Tol.",
    .multiplier = "Mult.",
    .sensitivity = "Sens.",
    .objective_sensitivity = "Obj. sens.",
    priority = "Priority",
    .relvar = "Rel. var.",
    .cv = "CV",
    .contribution = "Contribution",
    .share = "Share",
    lower = "Lower",
    continuous = "Cont.",
    field = "Field",
    supplied = "Supplied",
    upper = "Upper",
    continuous_status = "Cont. status",
    field_status = "Field status",
    lower_violation = "Lower bad",
    upper_violation = "Upper bad",
    status = "Status",
    pass = "Pass",
    mean = "Mean",
    variance = "Variance",
    deff = "Deff",
    response_rate = "Resp.",
    eligibility_rate = "Scope",
    eligible_per_unit = "Elig.",
    diagnostic = "Diagnostic",
    value = "Value",
    tolerance = "Tolerance"
  )
  hit <- names(out) %in% names(labels)
  names(out)[hit] <- unname(labels[names(out)[hit]])
  names(out) <- .alloc_take_labels(names(out))
  out
}

#' Labels for planned and fielded takes at every stage
#' @keywords internal
#' @noRd
.alloc_take_labels <- function(cols) {
  cols <- sub("^take_([0-9]+)_continuous$", "take_\\1 cont.", cols)
  sub("^take_([0-9]+)_field$", "take_\\1 field", cols)
}

#' Keep identifier columns from determining console width
#' @keywords internal
#' @noRd
.abbrev_bethel_id <- function(x, width = 20L) {
  x <- as.character(x)
  if (is.na(x) || nchar(x) <= width) {
    return(x)
  }
  paste0(substr(x, 1L, width - 3L), "...")
}

#' @rdname summary.svyplan_alloc
#' @export
print.summary.svyplan_alloc <- function(x, ...) {
  .check_unused_dots(...)
  cluster <- any(grepl("^n_1($|_)", names(x$allocation)))
  cat(
    if (identical(x$kind, "design")) {
      if (cluster) {
        "Stratified cluster allocation summary\n\n"
      } else {
        "Stratified allocation summary\n\n"
      }
    } else {
      if (cluster) {
        "Cluster allocation precision summary\n\n"
      } else {
        "Allocation precision summary\n\n"
      }
    }
  )
  cat(sprintf("Question: %s\n\n", x$question))

  .print_alloc_summary_overall(
    x$overall,
    if (identical(x$kind, "design")) {
      "Overall field design"
    } else {
      "Overall precision"
    },
    show_response = any(x$assumptions$response_rate < 1) ||
      !is.null(x$assumptions$eligibility_rate)
  )
  .print_search_skipped(x$search, indent = "  ")
  if (!is.null(x$continuous)) {
    cat("\nContinuous design\n")
    .print_alloc_summary_values(
      x$continuous,
      show_response = any(x$assumptions$response_rate < 1) ||
      !is.null(x$assumptions$eligibility_rate),
      indent = "  "
    )
  }

  cat("\nAllocation by stratum\n\n")
  alloc_main <- intersect(
    c(
      "stratum",
      "population",
      "response_rate",
      "eligibility_rate",
      "n_continuous",
      "n_field",
      "n_supplied",
      "expected_respondents"
    ),
    names(x$allocation)
  )
  print(
    .fmt_alloc_summary_table(x$allocation[alloc_main]),
    row.names = FALSE,
    quote = FALSE,
    right = TRUE
  )

  stage_cols <- c(
    "stratum",
    grep("^(n_1|take_[0-9]+)($|_continuous|_field)", names(x$allocation),
         value = TRUE)
  )
  if (length(stage_cols) > 1L) {
    cat("\nCluster stages by stratum\n\n")
    print(
      .fmt_alloc_summary_table(x$allocation[stage_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  fieldwork_cols <- intersect(
    c("stratum", "unit_cost", "weight", "cost"),
    names(x$allocation)
  )
  if (length(fieldwork_cols) > 1L) {
    cat("\nCost and weights by stratum\n\n")
    print(
      .fmt_alloc_summary_table(x$allocation[fieldwork_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  if (is.null(x$precision)) {
    cat("\nPrecision not assessed (no 'sd', 'var' or 'p' in frame)\n")
  } else {
    cat("\nAchieved precision by stratum\n\n")
    print(
      .fmt_alloc_summary_table(x$precision),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  if (!is.null(x$bounds)) {
    cat("\nAllocation bounds\n\n")
    print(
      .fmt_alloc_summary_table(x$bounds),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
    active <- sum(x$bounds$binding, na.rm = TRUE)
    cat(sprintf(
      "\n%s\n",
      if (active == 0L) {
        "No allocation bounds are active."
      } else {
        sprintf(
          "%d allocation bound%s active.",
          active,
          if (active == 1L) " is" else "s are"
        )
      }
    ))
  }

  if (!is.null(x$domains)) {
    measures <- c(
      ".domain",
      "population",
      "n",
      "expected_respondents",
      "se",
      "moe",
      "rmoe",
      "cv",
      "cost"
    )
    ids <- setdiff(names(x$domains), measures)
    key <- if (length(ids) > 0L) ids else ".domain"
    allocation_cols <- intersect(
      c(key, "population", "n", "expected_respondents", "cost"),
      names(x$domains)
    )
    precision_cols <- intersect(
      c(key, "se", "moe", "rmoe", "cv"),
      names(x$domains)
    )
    cat("\nDomain allocation (not additive to the overall row)\n\n")
    print(
      .fmt_alloc_summary_table(x$domains[allocation_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
    cat("\nDomain precision\n\n")
    print(
      .fmt_alloc_summary_table(x$domains[precision_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  cat("\nAssumptions\n")
  cat(sprintf("  alpha: %s\n", .fmt_prob(x$assumptions$alpha)))
  cat(sprintf(
    "  design effect: %s\n",
    .fmt_by_stratum(x$assumptions$deff, function(v) sprintf("%.2f", v))
  ))
  cat(sprintf(
    "  response rate: %s\n",
    .fmt_by_stratum(x$assumptions$response_rate, .fmt_prob)
  ))
  if (!is.null(x$assumptions$eligibility_rate)) {
    cat(sprintf(
      "  eligibility rate: %s\n",
      .fmt_by_stratum(x$assumptions$eligibility_rate, .fmt_prob)
    ))
  }
  if (!is.null(x$assumptions$eligible_per_unit)) {
    cat(sprintf(
      "  eligible per unit: %s\n",
      .fmt_by_stratum(x$assumptions$eligible_per_unit,
                      function(v) sprintf("%.4g", v))
    ))
  }
  if (!is.null(x$assumptions$df)) {
    cat(sprintf("  interval df: %g\n", x$assumptions$df))
  }
  .print_alloc_fpc(x$assumptions$fpc)
  invisible(x)
}

#' The finite population correction line of an allocation summary
#'
#' Printed only away from the default `"unit"`.
#' @keywords internal
#' @noRd
.print_alloc_fpc <- function(fpc) {
  if (!is.null(fpc) && !identical(fpc, "unit")) {
    cat(sprintf(
      "  finite population correction: %s\n",
      if (identical(fpc, "none")) "none" else "by stage"
    ))
  }
}

#' Print named overall allocation values
#' @keywords internal
#' @noRd
.print_alloc_summary_overall <- function(values, title, show_response) {
  cat(title, "\n", sep = "")
  .print_alloc_summary_values(values, show_response, indent = "  ")
}

#' @keywords internal
#' @noRd
.print_alloc_summary_values <- function(values, show_response, indent) {
  labels <- c(
    n = "Sample",
    expected_respondents = "Expected respondents",
    cost = "Cost",
    se = "SE",
    moe = "MOE",
    rmoe = "Relative MOE",
    cv = "CV",
    worst_domain_cv = "Worst domain CV",
    design_df = "Design df"
  )
  keep <- names(values)
  keep <- keep[keep %in% names(labels)]
  if (!show_response) {
    keep <- setdiff(keep, "expected_respondents")
  }
  keep <- keep[vapply(
    values[keep],
    function(v) {
      length(v) == 1L && !is.null(v) && !is.na(v)
    },
    logical(1L)
  )]
  width <- max(nchar(labels[keep]))
  for (name in keep) {
    value <- values[[name]]
    shown <- switch(
      name,
      n = .fmt_continuous_n(value),
      expected_respondents = sprintf("%.1f", value),
      cost = sprintf("%.1f", value),
      design_df = sprintf("%g", value),
      sprintf("%.4f", value)
    )
    cat(sprintf("%s%-*s  %s\n", indent, width, labels[[name]], shown))
  }
  invisible(NULL)
}

#' Format allocation summary tables without changing stored values
#' @keywords internal
#' @noRd
.fmt_alloc_summary_table <- function(tab) {
  out <- tab
  count_cols <- union(
    intersect(c("population", "n_field", "n_1_field", "field"), names(out)),
    grep("^take_[0-9]+_field$", names(out), value = TRUE)
  )
  for (name in count_cols) {
    out[[name]] <- vapply(out[[name]], .fmt_count_n, character(1L))
  }
  two_cols <- union(
    intersect(
      c(
        "unit_cost",
        "n_continuous",
        "n_supplied",
        "n_1_continuous",
        "n_1",
        "expected_respondents",
        "weight",
        "cost",
        "n"
      ),
      names(out)
    ),
    grep("^take_[0-9]+($|_continuous$)", names(out), value = TRUE)
  )
  for (name in two_cols) {
    out[[name]] <- ifelse(is.na(out[[name]]), "", sprintf("%.2f", out[[name]]))
  }
  for (name in intersect(c("lower", "upper"), names(out))) {
    out[[name]] <- ifelse(
      is.na(out[[name]]),
      "",
      vapply(out[[name]], .fmt_continuous_n, character(1L))
    )
  }
  if ("response_rate" %in% names(out)) {
    out$response_rate <- sprintf("%.3g", out$response_rate)
  }
  four_cols <- intersect(
    c("effective_n", "se", "moe", "rmoe", "cv"),
    names(out)
  )
  for (name in four_cols) {
    out[[name]] <- ifelse(is.na(out[[name]]), "", sprintf("%.4f", out[[name]]))
  }
  if ("variance_share" %in% names(out)) {
    out$variance_share <- ifelse(
      is.na(out$variance_share),
      "",
      sprintf("%.3f", out$variance_share)
    )
  }
  if ("binding" %in% names(out)) {
    out$binding <- ifelse(is.na(out$binding), "", ifelse(out$binding, "*", ""))
  }
  if ("source" %in% names(out)) {
    out$source[is.na(out$source)] <- ""
  }

  labels <- c(
    stratum = "Stratum",
    population = "Pop.",
    response_rate = "Resp.",
    eligibility_rate = "Scope",
    eligible_per_unit = "Elig.",
    unit_cost = "Cost/unit",
    n_continuous = "Cont.",
    n_field = "Field",
    n_supplied = "Supplied",
    n_1_continuous = "n_1 cont.",
    n_1_field = "n_1 field",
    n_1 = "n_1",
    expected_respondents = "Exp. resp.",
    weight = "Weight",
    cost = "Cost",
    effective_n = "Eff. n",
    se = "SE",
    moe = "MOE",
    rmoe = "Rel. MOE",
    cv = "CV",
    variance_share = "Var. share",
    lower = "Lower",
    field = "Field",
    upper = "Upper",
    binding = "Binding",
    source = "Source",
    .domain = "Domain",
    n = "Sample"
  )
  hit <- names(out) %in% names(labels)
  names(out)[hit] <- unname(labels[names(out)[hit]])
  names(out) <- .alloc_take_labels(names(out))
  out
}

