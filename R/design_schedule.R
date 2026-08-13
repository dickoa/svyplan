#' Construct an operational schedule from a rotating panel plan
#'
#' Translate the launch and recruitment quantities in a rotating
#' [svyplan_panel][n_panel()] into named cohort components, a dense
#' panel-by-occasion schedule, and any interview commitments beyond the
#' planning horizon. The result is a planning object: it does not inspect a
#' frame, execute a sample, or construct analysis weights.
#'
#' @param panel_plan A rotating `svyplan_panel` from [n_panel()] or
#'   [prec_panel()] with `start = "immediate"` or `start = "gradual"`.
#' @param life An unbroken equal-take `svyplan_overlap`, such as
#'   `design_overlap("4")`. Its length must match the panel life.
#' @param horizon A positive whole number of program occasions in the finite
#'   planning window.
#' @param horizon_policy One of `"continuing"`, `"truncate_lives"`, or
#'   `"close_intake"`. Under `"continuing"`, interviews owed after `horizon`
#'   are returned separately in `tail_commitments`. `"truncate_lives"` cuts
#'   those interviews at the horizon. `"close_intake"` stops recruitment soon
#'   enough for every recruited cohort to finish by the horizon.
#' @param refreshment Whether later cohorts are drawn from an
#'   `"entrant_register"` or a `"whole_vintage"`. This records a frame role;
#'   it does not establish disjointness or combine weights.
#' @param frame_vintage `NULL`, or a named character vector supplying frame
#'   identifiers for any subset of the generated cohort names. Unnamed or
#'   unknown entries are refused rather than matched by position.
#' @param rounding The explicit operational rounding rule. The first schema
#'   supports `"ceiling"`, applied to each startup panel and intake cohort
#'   before totals are formed. Operational counts are stored as whole-valued
#'   numerics so plans above R's 32-bit integer range remain representable.
#'
#' @return An `svyplan_schedule` list with scalar planning metadata and four
#'   data frames:
#' \describe{
#'   \item{`components`}{One row per execution: startup followed by the intake
#'     cohorts permitted by the horizon policy. Continuous and operational
#'     issue remain separate.}
#'   \item{`schedule`}{The dense component-panel grid over occasions 1 through
#'     `horizon`. Inactive rows have `active = FALSE` and `life_stage = NA`.}
#'   \item{`tail_commitments`}{A sparse table of interviews after `horizon`
#'     owed to components entering within a `"continuing"` window. It has a
#'     fixed empty schema under the other policies.}
#'   \item{`issue`}{The continuous and operational recruitment issued at each
#'     in-horizon occasion. [as.data.frame()] returns this view.}
#' }
#'
#' The generated activity is checked against `panel_plan$launch_waves` over
#' every shared period that the horizon policy has not changed. Its realized
#' issued overlap is also checked against `life` at every mature membership
#' origin for which the destination occasion exists.
#'
#' @examples
#' panel <- n_panel(
#'   n_prop(p = 0.5, moe = 0.03),
#'   retention = c(0.9, 0.9, 0.9),
#'   resp_rate = 0.75,
#'   design = "rotating",
#'   start = "immediate"
#' )
#' schedule <- design_schedule(
#'   panel,
#'   design_overlap("4"),
#'   horizon = 6,
#'   horizon_policy = "continuing",
#'   refreshment = "entrant_register",
#'   rounding = "ceiling"
#' )
#' schedule
#' schedule$components
#'
#' @family repeated survey planning
#' @export
design_schedule <- function(
  panel_plan,
  life,
  horizon,
  horizon_policy,
  refreshment,
  frame_vintage = NULL,
  rounding
) {
  if (!inherits(panel_plan, "svyplan_panel") ||
      !identical(panel_plan$design, "rotating")) {
    stop("'panel_plan' must be a rotating svyplan_panel", call. = FALSE)
  }
  if (is.null(panel_plan$start) || is.null(panel_plan$launch) ||
      is.null(panel_plan$launch_waves)) {
    stop(
      "'panel_plan' must carry a launch; call n_panel() or prec_panel() with start = \"immediate\" or \"gradual\"",
      call. = FALSE
    )
  }
  if (!inherits(life, "svyplan_overlap")) {
    stop("'life' must be a svyplan_overlap from design_overlap()", call. = FALSE)
  }

  take <- attr(life, "schedule", exact = TRUE)
  if (any(take == 0)) {
    stop("'life' must be unbroken in the first schedule schema", call. = FALSE)
  }
  if (length(unique(signif(take, 14L))) != 1L) {
    stop("'life' must be equal-take in the first schedule schema", call. = FALSE)
  }
  life_length <- length(take)
  if (life_length != nrow(panel_plan$waves)) {
    stop(
      sprintf(
        "the overlap life length is %d but 'panel_plan' has a %d-wave life",
        life_length, nrow(panel_plan$waves)
      ),
      call. = FALSE
    )
  }

  horizon <- .schedule_whole(horizon, "horizon")
  if (missing(horizon_policy)) {
    stop("'horizon_policy' must be declared explicitly", call. = FALSE)
  }
  horizon_policy <- .schedule_choice(
    horizon_policy,
    c("continuing", "truncate_lives", "close_intake"),
    "horizon_policy"
  )
  if (identical(horizon_policy, "close_intake") && horizon < life_length) {
    stop(
      "'close_intake' requires a horizon no shorter than the full cohort life",
      call. = FALSE
    )
  }
  if (missing(refreshment)) {
    stop("'refreshment' must be declared explicitly", call. = FALSE)
  }
  refreshment <- .schedule_choice(
    refreshment,
    c("entrant_register", "whole_vintage"),
    "refreshment"
  )
  if (missing(rounding)) {
    stop("'rounding' must be declared explicitly", call. = FALSE)
  }
  rounding <- .schedule_choice(rounding, "ceiling", "rounding")

  entrant_issue <- unname(panel_plan$n_entrants)
  # Whole-valued doubles remain exact far beyond R's 32-bit integer range.
  # Operational counts are counts by value, not storage-mode integers.
  panel_issue <- ceiling(entrant_issue)
  max_entry <- if (identical(horizon_policy, "close_intake")) {
    horizon - life_length + 1L
  } else {
    horizon
  }

  component_names <- c(
    "startup",
    if (max_entry >= 2L) paste0("intake_", 2:max_entry)
  )
  entry_wave <- c(1L, if (max_entry >= 2L) 2:max_entry)
  startup_panels <- if (identical(panel_plan$start, "immediate")) {
    life_length
  } else {
    1L
  }
  panels <- c(startup_panels, rep(1L, length(component_names) - 1L))
  frame_role <- c("startup", rep(refreshment, length(component_names) - 1L))
  vintages <- .schedule_frame_vintages(frame_vintage, component_names)

  components <- data.frame(
    cohort = component_names,
    entry_wave = as.integer(entry_wave),
    frame_role = frame_role,
    planned_issue = entrant_issue * panels,
    operational_issue = panel_issue * panels,
    panels = as.integer(panels),
    panel_issue = rep(panel_issue, length(component_names)),
    frame_vintage = vintages,
    status = rep("planned", length(component_names)),
    stringsAsFactors = FALSE
  )

  catalog <- .schedule_panel_catalog(
    components, panel_plan$start, life_length, horizon, horizon_policy
  )
  schedule <- .schedule_dense_grid(catalog, horizon, life_length)
  tail_commitments <- .schedule_tail(catalog, horizon, horizon_policy)
  steady <- .schedule_steady(schedule, life_length)
  schedule$steady_state <- steady[schedule$wave]
  schedule <- schedule[c(
    "wave", "cohort", "panel", "life_stage", "life_length", "active",
    "steady_state", "relative_take"
  )]

  steady_waves <- which(steady)
  steady_state_from <- if (length(steady_waves) == 0L) {
    NA_integer_
  } else {
    as.integer(steady_waves[[1L]])
  }

  issue <- .schedule_issue_profile(components, horizon, steady)
  panel_parameters <- .schedule_panel_parameters(
    schedule, tail_commitments, startup_panels
  )

  out <- .new_svyplan_schedule(
    schema_version = 1L,
    life_spec = as.numeric(take),
    life_notation = "resolved_per_occasion_take",
    life_length = life_length,
    horizon = horizon,
    horizon_policy = horizon_policy,
    launch_policy = panel_plan$start,
    target_estimand = list(type = panel_plan$type, method = panel_plan$method),
    n_target = unname(panel_plan$n_target),
    response_model = list(
      resp_rate = panel_plan$params$resp_rate,
      retention = panel_plan$params$retention,
      q = panel_plan$waves$q
    ),
    rounding = list(rule = rounding, level = "panel_and_cohort"),
    steady_state_from = steady_state_from,
    overlap = as.data.frame(life),
    refreshment = refreshment,
    requires_disjoint_entrants = identical(refreshment, "entrant_register"),
    n_entrants = entrant_issue,
    n_in_sample = unname(panel_plan$n_in_sample),
    n_cohorts = as.integer(panel_plan$n_cohorts),
    panel_parameters = panel_parameters,
    components = components,
    schedule = schedule,
    tail_commitments = tail_commitments,
    issue = issue
  )

  .check_schedule_launch(out, panel_plan$launch_waves, horizon_policy, max_entry)
  .check_schedule_overlap(out, as.double(life))
  out
}

#' @keywords internal
#' @noRd
.schedule_whole <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != trunc(x)) {
    stop(sprintf("'%s' must be a positive whole number", name), call. = FALSE)
  }
  as.integer(x)
}

#' @keywords internal
#' @noRd
.schedule_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    stop(
      sprintf("'%s' must be one of %s", name, paste(sprintf("\"%s\"", choices), collapse = ", ")),
      call. = FALSE
    )
  }
  x
}

#' @keywords internal
#' @noRd
.schedule_frame_vintages <- function(x, cohorts) {
  out <- rep(NA_character_, length(cohorts))
  if (is.null(x)) {
    return(out)
  }
  ok <- is.character(x) && !anyNA(x) && length(x) > 0L &&
    !is.null(names(x)) && all(nzchar(names(x))) && anyDuplicated(names(x)) == 0L
  if (!ok || !all(names(x) %in% cohorts)) {
    stop(
      "'frame_vintage' must be NULL or a named character vector whose names are generated cohorts",
      call. = FALSE
    )
  }
  out[match(names(x), cohorts)] <- unname(x)
  out
}

#' @keywords internal
#' @noRd
.schedule_panel_catalog <- function(components, start, life_length, horizon,
                                    horizon_policy) {
  rows <- vector("list", nrow(components))
  for (i in seq_len(nrow(components))) {
    component <- components[i, ]
    k <- component$panels
    lengths <- if (i == 1L && identical(start, "immediate")) {
      rev(seq_len(life_length))
    } else {
      rep(life_length, k)
    }
    if (identical(horizon_policy, "truncate_lives")) {
      lengths <- pmin(lengths, horizon - component$entry_wave + 1L)
    }
    rows[[i]] <- data.frame(
      cohort = rep(component$cohort, k),
      panel = seq_len(k),
      entry_wave = rep(component$entry_wave, k),
      life_length = as.integer(lengths),
      relative_take = rep(1, k),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.schedule_dense_grid <- function(catalog, horizon, life_length) {
  rows <- lapply(seq_len(horizon), function(wave) {
    out <- catalog
    out$wave <- wave
    stage <- wave - out$entry_wave + 1L
    out$active <- stage >= 1L & stage <= out$life_length
    out$life_stage <- ifelse(out$active, stage, NA_integer_)
    out[c(
      "wave", "cohort", "panel", "life_stage", "life_length", "active",
      "relative_take"
    )]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$wave <- as.integer(out$wave)
  out$panel <- as.integer(out$panel)
  out$life_stage <- as.integer(out$life_stage)
  out
}

#' @keywords internal
#' @noRd
.empty_schedule_tail <- function() {
  data.frame(
    wave = integer(),
    cohort = character(),
    panel = integer(),
    life_stage = integer(),
    life_length = integer(),
    relative_take = numeric(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.schedule_tail <- function(catalog, horizon, policy) {
  if (!identical(policy, "continuing")) {
    return(.empty_schedule_tail())
  }
  rows <- lapply(seq_len(nrow(catalog)), function(i) {
    z <- catalog[i, ]
    last <- z$entry_wave + z$life_length - 1L
    if (last <= horizon) {
      return(NULL)
    }
    waves <- seq.int(horizon + 1L, last)
    data.frame(
      wave = as.integer(waves),
      cohort = rep(z$cohort, length(waves)),
      panel = rep(as.integer(z$panel), length(waves)),
      life_stage = as.integer(waves - z$entry_wave + 1L),
      life_length = rep(as.integer(z$life_length), length(waves)),
      relative_take = rep(z$relative_take, length(waves)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(.empty_schedule_tail())
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$wave, match(out$cohort, unique(catalog$cohort)), out$panel), ]
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.schedule_steady <- function(schedule, life_length) {
  vapply(seq_len(max(schedule$wave)), function(wave) {
    at <- schedule[
      schedule$wave == wave & schedule$active,
      c("life_stage", "relative_take"),
      drop = FALSE
    ]
    totals <- numeric(life_length)
    if (nrow(at) > 0L) {
      by_stage <- tapply(at$relative_take, at$life_stage, sum)
      totals[as.integer(names(by_stage))] <- as.numeric(by_stage)
    }
    isTRUE(all.equal(totals, rep(1, life_length), tolerance = 1e-12))
  }, logical(1L))
}

#' @keywords internal
#' @noRd
.schedule_issue_profile <- function(components, horizon, steady) {
  rows <- lapply(seq_len(horizon), function(wave) {
    at <- components[components$entry_wave == wave, , drop = FALSE]
    data.frame(
      wave = wave,
      cohort = if (nrow(at) == 0L) NA_character_ else at$cohort[[1L]],
      planned_issue = if (nrow(at) == 0L) 0 else at$planned_issue[[1L]],
      operational_issue = if (nrow(at) == 0L) 0 else at$operational_issue[[1L]],
      steady_state = steady[[wave]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$wave <- as.integer(out$wave)
  out
}

#' @keywords internal
#' @noRd
.schedule_panel_parameters <- function(schedule, tail, k) {
  if (k < 2L) {
    return(list(k = as.integer(k), r_min = 1L, block_width = NA_integer_))
  }
  startup <- schedule[schedule$cohort == "startup" & schedule$active, ]
  if (nrow(tail) > 0L) {
    startup <- rbind(
      startup[c("wave", "cohort", "panel", "life_stage", "life_length",
                "relative_take")],
      tail[tail$cohort == "startup", , drop = FALSE]
    )
  }
  per_wave <- table(startup$wave)
  r_min <- as.integer(min(per_wave))
  list(
    k = as.integer(k),
    r_min = r_min,
    block_width = as.integer(k * ceiling(2 / r_min))
  )
}

#' @keywords internal
#' @noRd
.new_svyplan_schedule <- function(...) {
  out <- list(...)
  .validate_schedule_schema(out)
  structure(out, class = c("svyplan_schedule", "list"))
}

#' Validate the package-neutral schedule schema
#' @keywords internal
#' @noRd
.validate_schedule_schema <- function(x) {
  fail <- function(detail) {
    stop(
      sprintf("internal schedule defect: %s", detail),
      call. = FALSE
    )
  }
  same <- function(a, b) {
    isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = TRUE))
  }
  scalar_number <- function(z, lower = -Inf, upper = Inf) {
    is.numeric(z) && length(z) == 1L && !is.na(z) && is.finite(z) &&
      z >= lower && z <= upper
  }
  whole_counts <- function(z, lower = 0) {
    is.numeric(z) && !anyNA(z) && all(is.finite(z)) &&
      all(z >= lower) && all(z == floor(z))
  }

  required <- c(
    "schema_version", "life_spec", "life_notation", "life_length",
    "horizon", "horizon_policy", "launch_policy", "target_estimand",
    "n_target", "response_model", "rounding", "steady_state_from",
    "overlap", "refreshment", "requires_disjoint_entrants", "n_entrants",
    "n_in_sample", "n_cohorts", "panel_parameters", "components",
    "schedule", "tail_commitments", "issue"
  )
  if (!is.list(x) || !identical(names(x), required)) {
    fail("top-level fields do not match schema version 1")
  }
  if (!identical(x$schema_version, 1L) ||
      !is.integer(x$life_length) || length(x$life_length) != 1L ||
      is.na(x$life_length) || x$life_length < 2L ||
      !is.integer(x$horizon) || length(x$horizon) != 1L ||
      is.na(x$horizon) || x$horizon < 1L ||
      !identical(x$life_notation, "resolved_per_occasion_take") ||
      !is.character(x$horizon_policy) || length(x$horizon_policy) != 1L ||
      is.na(x$horizon_policy) ||
      !x$horizon_policy %in% c("continuing", "truncate_lives", "close_intake") ||
      !is.character(x$launch_policy) || length(x$launch_policy) != 1L ||
      is.na(x$launch_policy) ||
      !x$launch_policy %in% c("immediate", "gradual")) {
    fail("scalar metadata does not match schema version 1")
  }
  if (!is.numeric(x$life_spec) || length(x$life_spec) != x$life_length ||
      anyNA(x$life_spec) || any(!is.finite(x$life_spec)) ||
      any(x$life_spec <= 0) ||
      length(unique(signif(x$life_spec, 14L))) != 1L) {
    fail("life_spec must be an unbroken equal-take life")
  }
  if (identical(x$horizon_policy, "close_intake") &&
      x$horizon < x$life_length) {
    fail("close_intake has a horizon shorter than the cohort life")
  }
  if (!is.list(x$target_estimand) ||
      !identical(names(x$target_estimand), c("type", "method")) ||
      !is.character(x$target_estimand$type) ||
      length(x$target_estimand$type) != 1L || is.na(x$target_estimand$type) ||
      !(is.null(x$target_estimand$method) ||
        (is.character(x$target_estimand$method) &&
           length(x$target_estimand$method) == 1L &&
           !is.na(x$target_estimand$method))) ||
      !scalar_number(x$n_target, lower = 0) || x$n_target == 0) {
    fail("target metadata is inconsistent")
  }
  response <- x$response_model
  if (!is.list(response) ||
      !identical(names(response), c("resp_rate", "retention", "q")) ||
      !scalar_number(response$resp_rate, lower = 0, upper = 1) ||
      response$resp_rate == 0 ||
      !is.numeric(response$retention) ||
      length(response$retention) != x$life_length - 1L ||
      anyNA(response$retention) || any(!is.finite(response$retention)) ||
      any(response$retention <= 0 | response$retention > 1) ||
      !is.numeric(response$q) || length(response$q) != x$life_length ||
      anyNA(response$q) || any(!is.finite(response$q))) {
    fail("response metadata is inconsistent")
  }
  expected_q <- response$resp_rate * cumprod(c(1, response$retention))
  if (!same(as.numeric(response$q), expected_q)) {
    fail("response probabilities do not follow resp_rate and retention")
  }
  if (!is.list(x$rounding) ||
      !identical(x$rounding, list(
        rule = "ceiling", level = "panel_and_cohort"
      ))) {
    fail("rounding metadata does not match schema version 1")
  }
  if (!is.character(x$refreshment) || length(x$refreshment) != 1L ||
      is.na(x$refreshment) ||
      !x$refreshment %in% c("entrant_register", "whole_vintage") ||
      !is.logical(x$requires_disjoint_entrants) ||
      length(x$requires_disjoint_entrants) != 1L ||
      is.na(x$requires_disjoint_entrants) ||
      !identical(
        x$requires_disjoint_entrants,
        identical(x$refreshment, "entrant_register")
      )) {
    fail("refreshment metadata is inconsistent")
  }
  if (!scalar_number(x$n_entrants, lower = 0) || x$n_entrants == 0 ||
      !scalar_number(x$n_in_sample, lower = 0) || x$n_in_sample == 0 ||
      !identical(x$n_cohorts, x$life_length) ||
      !same(x$n_in_sample, x$n_cohorts * x$n_entrants)) {
    fail("panel counts are inconsistent")
  }

  component_cols <- c(
    "cohort", "entry_wave", "frame_role", "planned_issue",
    "operational_issue", "panels", "panel_issue", "frame_vintage", "status"
  )
  schedule_cols <- c(
    "wave", "cohort", "panel", "life_stage", "life_length", "active",
    "steady_state", "relative_take"
  )
  tail_cols <- c(
    "wave", "cohort", "panel", "life_stage", "life_length", "relative_take"
  )
  issue_cols <- c(
    "wave", "cohort", "planned_issue", "operational_issue", "steady_state"
  )
  exact_table <- function(tab, cols) {
    is.data.frame(tab) && identical(names(tab), cols)
  }
  if (!exact_table(x$components, component_cols) ||
      !exact_table(x$schedule, schedule_cols) ||
      !exact_table(x$tail_commitments, tail_cols) ||
      !exact_table(x$issue, issue_cols)) {
    fail("one or more tables do not have the fixed schema-version-1 columns")
  }

  components <- x$components
  max_entry <- if (identical(x$horizon_policy, "close_intake")) {
    x$horizon - x$life_length + 1L
  } else {
    x$horizon
  }
  expected_cohorts <- c(
    "startup", if (max_entry >= 2L) paste0("intake_", 2:max_entry)
  )
  startup_panels <- if (identical(x$launch_policy, "immediate")) {
    x$life_length
  } else {
    1L
  }
  expected_panels <- c(
    startup_panels, rep(1L, length(expected_cohorts) - 1L)
  )
  expected_roles <- c(
    "startup", rep(x$refreshment, length(expected_cohorts) - 1L)
  )
  if (!identical(components$cohort, expected_cohorts) ||
      anyNA(components$cohort) || anyDuplicated(components$cohort) ||
      !is.integer(components$entry_wave) || anyNA(components$entry_wave) ||
      !identical(components$entry_wave, seq_len(nrow(components))) ||
      !is.numeric(components$planned_issue) ||
      any(!is.finite(components$planned_issue) | components$planned_issue <= 0) ||
      !whole_counts(components$operational_issue, lower = 1) ||
      any(components$operational_issue < components$planned_issue) ||
      !identical(components$panels, expected_panels) ||
      !whole_counts(components$panel_issue, lower = 1) ||
      !is.character(components$frame_vintage) ||
      length(components$frame_vintage) != nrow(components) ||
      any(!is.na(components$frame_vintage) &
            !nzchar(components$frame_vintage)) ||
      !identical(components$frame_role, expected_roles) ||
      !identical(components$status, rep("planned", nrow(components)))) {
    fail("the component manifest is inconsistent")
  }
  expected_panel_issue <- rep(ceiling(x$n_entrants), nrow(components))
  if (!same(components$planned_issue, x$n_entrants * expected_panels) ||
      !same(components$panel_issue, expected_panel_issue) ||
      !same(
        components$operational_issue,
        expected_panel_issue * expected_panels
      )) {
    fail("component issue does not follow panel-level ceiling")
  }

  schedule <- x$schedule
  n_panels <- sum(components$panels)
  if (!is.integer(schedule$wave) || !is.integer(schedule$panel) ||
      !is.integer(schedule$life_stage) || !is.integer(schedule$life_length) ||
      !is.logical(schedule$active) || !is.logical(schedule$steady_state) ||
      !is.numeric(schedule$relative_take) ||
      nrow(schedule) != x$horizon * n_panels ||
      anyNA(schedule[c("wave", "cohort", "panel", "life_length", "active",
                       "steady_state", "relative_take")]) ||
      any(!schedule$cohort %in% components$cohort) ||
      any(schedule$wave < 1L | schedule$wave > x$horizon) ||
      any(schedule$panel < 1L) ||
      any(schedule$life_length < 1L | schedule$life_length > x$life_length) ||
      any(schedule$relative_take != 1) ||
      any(is.na(schedule$life_stage) != !schedule$active) ||
      any(schedule$life_stage[schedule$active] < 1L |
            schedule$life_stage[schedule$active] >
              schedule$life_length[schedule$active])) {
    fail("the dense schedule grid is inconsistent")
  }
  grid_key <- paste(schedule$wave, schedule$cohort, schedule$panel, sep = "|")
  if (anyDuplicated(grid_key)) {
    fail("the dense schedule grid contains duplicate component-panel rows")
  }

  catalog <- .schedule_panel_catalog(
    components, x$launch_policy, x$life_length, x$horizon, x$horizon_policy
  )
  expected_schedule <- .schedule_dense_grid(
    catalog, x$horizon, x$life_length
  )
  expected_steady <- .schedule_steady(expected_schedule, x$life_length)
  expected_schedule$steady_state <- expected_steady[expected_schedule$wave]
  expected_schedule <- expected_schedule[names(schedule)]
  if (!identical(schedule, expected_schedule)) {
    fail("the dense schedule does not match the component manifest")
  }

  tail <- x$tail_commitments
  if (!is.integer(tail$wave) || !is.integer(tail$panel) ||
      !is.integer(tail$life_stage) || !is.integer(tail$life_length) ||
      !is.numeric(tail$relative_take) || anyNA(tail) ||
      any(tail$wave <= x$horizon) ||
      any(!tail$cohort %in% components$cohort) ||
      any(tail$panel < 1L) ||
      any(tail$life_stage < 1L | tail$life_stage > tail$life_length) ||
      any(tail$life_length < 1L | tail$life_length > x$life_length) ||
      any(tail$relative_take != 1) ||
      (!identical(x$horizon_policy, "continuing") && nrow(tail) > 0L)) {
    fail("the tail-commitment table is inconsistent with the horizon policy")
  }
  expected_tail <- .schedule_tail(
    catalog, x$horizon, x$horizon_policy
  )
  if (!identical(tail, expected_tail)) {
    fail("tail commitments do not match the component manifest")
  }

  issue <- x$issue
  if (!is.integer(issue$wave) ||
      !identical(issue$wave, seq_len(x$horizon)) ||
      !is.numeric(issue$planned_issue) ||
      !whole_counts(issue$operational_issue) ||
      !is.logical(issue$steady_state) || anyNA(issue$steady_state) ||
      any(issue$planned_issue < 0) ||
      any(issue$operational_issue < issue$planned_issue)) {
    fail("the occasion issue profile is inconsistent")
  }
  expected_issue <- .schedule_issue_profile(
    components, x$horizon, expected_steady
  )
  if (!identical(issue, expected_issue)) {
    fail("the occasion issue profile does not match the component manifest")
  }

  expected_steady_waves <- which(expected_steady)
  expected_steady_from <- if (length(expected_steady_waves) == 0L) {
    NA_integer_
  } else {
    as.integer(expected_steady_waves[[1L]])
  }
  if (!is.integer(x$steady_state_from) ||
      length(x$steady_state_from) != 1L ||
      !identical(x$steady_state_from, expected_steady_from)) {
    fail("steady-state metadata does not match the schedule")
  }

  expected_overlap <- as.data.frame(design_overlap(x$life_spec))
  if (!identical(x$overlap, expected_overlap)) {
    fail("overlap metadata does not match life_spec")
  }
  expected_parameters <- .schedule_panel_parameters(
    expected_schedule, expected_tail, startup_panels
  )
  if (!identical(x$panel_parameters, expected_parameters)) {
    fail("panel parameters do not match the schedule")
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.schedule_active_aggregate <- function(x) {
  tab <- x$schedule[x$schedule$active, , drop = FALSE]
  q <- x$response_model$q
  out <- data.frame(
    period = tab$wave,
    wave = tab$life_stage,
    n_issued = x$n_entrants * tab$relative_take,
    n_resp = x$n_entrants * tab$relative_take * q[tab$life_stage]
  )
  out <- stats::aggregate(
    cbind(n_issued, n_resp) ~ period + wave,
    data = out,
    FUN = sum
  )
  out <- out[order(out$period, out$wave), ]
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.check_schedule_launch <- function(x, launch_waves, policy, max_entry) {
  top <- min(x$horizon, x$life_length + 1L)
  if (identical(policy, "close_intake")) {
    top <- min(top, max_entry)
  }
  got <- .schedule_active_aggregate(x)
  got <- got[got$period <= top, , drop = FALSE]
  want <- launch_waves[launch_waves$period <= top, , drop = FALSE]
  rownames(want) <- NULL
  if (!isTRUE(all.equal(got, want, tolerance = 1e-10, check.attributes = FALSE))) {
    stop(
      "internal schedule defect: generated activity does not reconcile with panel_plan$launch_waves",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.check_schedule_overlap <- function(x, expected) {
  active <- x$schedule[x$schedule$active, , drop = FALSE]
  key <- paste(active$cohort, active$panel, sep = "|")
  take <- vapply(seq_len(x$horizon), function(w) {
    sum(active$relative_take[active$wave == w])
  }, numeric(1L))
  for (lag in seq_len(min(length(expected), x$horizon - 1L))) {
    origins <- which(
      seq_len(x$horizon) + lag <= x$horizon &
        abs(take - x$life_length) <= 1e-10
    )
    for (origin in origins) {
      first <- unique(key[active$wave == origin])
      second <- unique(key[active$wave == origin + lag])
      overlap <- length(intersect(first, second)) / length(first)
      if (!isTRUE(all.equal(overlap, expected[[lag]], tolerance = 1e-10))) {
        stop(
          sprintf(
            "internal schedule defect: overlap at origin %d and lag %d does not match the declared life",
            origin, lag
          ),
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

#' Print, summarise and coerce a longitudinal design schedule
#'
#' `print()` states the issue profile as a run of takes rather than as one
#' row per occasion, because a schedule's occasions are mostly identical and
#' its length is set by the reporting horizon rather than by the design.
#' `summary()` gives the occasion-by-occasion tables: the issue profile in
#' full, the component activity, the overlap the rotation produces, and the
#' interviews owed after the horizon. [as.data.frame()] returns the issue
#' profile, and every table stays reachable as a field.
#'
#' @param x An `svyplan_schedule`, or the `summary.svyplan_schedule` object
#'   `summary()` returns.
#' @param object An `svyplan_schedule`.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   [as.data.frame()] arguments.
#' @param name,i,value Standard replacement arguments for `$<-`, `[[<-`, and
#'   `[<-`. These operations are refused because they would break
#'   reconciliation within the schedule object.
#' @param ... Additional arguments are not supported and produce an error.
#' @return `print()` returns its argument invisibly, `summary()` an object of
#'   class `summary.svyplan_schedule` carrying the schedule and its four
#'   tables, `format()` one descriptive string, and `as.data.frame()` the
#'   issue profile.
#'
#' Direct field and table replacement with `$<-`, `[[<-`, or `[<-` is refused
#' because the metadata and tables describe one reconciled design. Recompute
#' with [design_schedule()] to change the plan, or extract a table first to
#' modify a plain data frame.
#' @name print.svyplan_schedule
NULL

#' @rdname print.svyplan_schedule
#' @export
print.svyplan_schedule <- function(x, ...) {
  .check_unused_dots(...)
  cat(sprintf(
    "Longitudinal design schedule (%s launch, %s)\n",
    x$launch_policy, x$horizon_policy
  ))
  cat(.fmt_schedule_life(x))
  cat(sprintf(
    "rounding: %s at %s level\n",
    x$rounding$rule, gsub("_", " ", x$rounding$level, fixed = TRUE)
  ))
  cat(sprintf("issue: %s\n", .fmt_schedule_issue(x)))
  cat(.fmt_schedule_tail(x))
  cat("# summary() for the occasion-by-occasion profile and the overlap\n")
  invisible(x)
}

#' The life, and where the response composition settles
#' @keywords internal
#' @noRd
.fmt_schedule_life <- function(x) {
  sprintf(
    "life: %d stages over %d occasions, %s\n",
    x$life_length, x$horizon,
    if (is.na(x$steady_state_from)) {
      "steady state not reached in the window"
    } else {
      sprintf("steady from occasion %d", x$steady_state_from)
    }
  )
}

#' The issue profile as the runs it is made of
#'
#' A schedule's length is set by the reporting horizon, not by the design, so
#' a row per occasion grows without bound while saying the same thing: a
#' startup take, then one intake repeated. Naming the runs states the profile
#' at a size that does not move when the horizon does. A run of zeros is
#' intake having closed, which `horizon_policy = "close_intake"` produces and
#' which is a fact about the schedule rather than a gap in it.
#' @keywords internal
#' @noRd
.fmt_schedule_issue <- function(x) {
  v <- x$issue$operational_issue
  parts <- sprintf("%s at startup", format(v[[1L]]))
  if (length(v) > 1L) {
    runs <- rle(v[-1L])
    at <- 2L
    for (j in seq_along(runs$lengths)) {
      lo <- at
      hi <- at + runs$lengths[[j]] - 1L
      at <- hi + 1L
      parts <- c(parts, if (runs$values[[j]] == 0) {
        sprintf("intake closed from %d", lo)
      } else {
        sprintf(
          "%s per occasion (%s)", format(runs$values[[j]]),
          if (lo == hi) sprintf("occasion %d", lo) else {
            sprintf("occasions %d-%d", lo, hi)
          }
        )
      })
    }
  }
  paste(parts, collapse = ", ")
}

#' @keywords internal
#' @noRd
.fmt_schedule_tail <- function(x) {
  if (nrow(x$tail_commitments) == 0L) {
    return("tail commitments: none\n")
  }
  sprintf(
    "tail commitments: %d panel-interviews after occasion %d\n",
    nrow(x$tail_commitments), x$horizon
  )
}

#' Drop the component columns that carry no distinction
#'
#' `frame_vintage` is absent unless one was declared and `status` is one
#' value for a plan that has not been executed against, so on an ordinary
#' schedule both are a column of the same entry repeated.
#' @keywords internal
#' @noRd
.fmt_schedule_components <- function(d) {
  if (all(is.na(d$frame_vintage))) {
    d$frame_vintage <- NULL
  }
  if (length(unique(d$status)) == 1L) {
    d$status <- NULL
  }
  d
}

#' @rdname print.svyplan_schedule
#' @export
summary.svyplan_schedule <- function(object, ...) {
  .check_unused_dots(...)
  structure(
    list(
      schedule = object,
      issue = object$issue,
      components = object$components,
      overlap = object$overlap,
      tail_commitments = object$tail_commitments,
      activity = object$schedule
    ),
    class = "summary.svyplan_schedule"
  )
}

#' @rdname print.svyplan_schedule
#' @export
print.summary.svyplan_schedule <- function(x, ...) {
  .check_unused_dots(...)
  s <- x$schedule
  cat(sprintf(
    "Analysis of a longitudinal design schedule (%s launch, %s)\n\n",
    s$launch_policy, s$horizon_policy
  ))
  cat(.fmt_schedule_life(s))
  cat(sprintf(
    "rounding: %s at %s level\n",
    s$rounding$rule, gsub("_", " ", s$rounding$level, fixed = TRUE)
  ))
  # Whole units, and the same ones the issue profile shows: the standing
  # sample is the cohorts times the take each one was rounded to, not the
  # continuous total rounded once.
  entrants <- ceiling(s$n_entrants)
  cat(sprintf(
    "entrants: %s per occasion, %s in sample across %d cohorts\n",
    format(entrants), format(entrants * s$n_cohorts), s$n_cohorts
  ))
  cat(sprintf("refreshment: %s\n", s$refreshment))

  cat("\nIssue profile\n")
  print(x$issue, row.names = FALSE, right = TRUE)

  cat("\nComponents\n")
  print(.fmt_schedule_components(x$components), row.names = FALSE,
        right = TRUE)

  cat("\nOverlap the rotation produces\n")
  print(as.data.frame(x$overlap), row.names = FALSE, right = TRUE)

  cat(.fmt_schedule_tail(s))
  cat("# $activity for the panel-level schedule, $tail_commitments for the tail\n")
  invisible(x)
}

#' @rdname print.svyplan_schedule
#' @export
format.svyplan_schedule <- function(x, ...) {
  .check_unused_dots(...)
  sprintf(
    "svyplan_schedule [%s, %s, %d-stage life, horizon %d]",
    x$launch_policy, x$horizon_policy, x$life_length, x$horizon
  )
}

#' @rdname print.svyplan_schedule
#' @export
as.data.frame.svyplan_schedule <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  out <- x$issue
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}

#' Refuse mutation of a reconciled schedule object
#' @keywords internal
#' @noRd
.no_schedule_replacement <- function() {
  stop(
    paste0(
      "a longitudinal design schedule cannot be modified in place: its ",
      "metadata and tables describe one reconciled design; recompute with ",
      "design_schedule(), or extract the table you want to modify"
    ),
    call. = FALSE
  )
}

#' @rdname print.svyplan_schedule
#' @export
`$<-.svyplan_schedule` <- function(x, name, value) {
  .no_schedule_replacement()
}

#' @rdname print.svyplan_schedule
#' @export
`[<-.svyplan_schedule` <- function(x, i, value) {
  .no_schedule_replacement()
}

#' @rdname print.svyplan_schedule
#' @export
`[[<-.svyplan_schedule` <- function(x, i, value) {
  .no_schedule_replacement()
}
