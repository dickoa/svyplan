#' Construct an operational schedule from a rotating panel plan
#'
#' Translate the launch and recruitment quantities in a rotating
#' [svyplan_panel][n_panel()] into named cohort components, a dense
#' panel-by-occasion schedule, and any interview commitments beyond the
#' planning horizon. The result is a planning object. It does not inspect
#' a frame, execute a sample, or construct analysis weights.
#'
#' @param panel_plan A rotating `svyplan_panel` from [n_panel()] or
#'   [prec_panel()] with `start = "immediate"` or `start = "gradual"`. Its
#'   `rotation` is the life the schedule lays out, and the overlap the
#'   manifest records is computed from it, so neither is supplied here.
#' @param horizon A positive whole number of program occasions in the finite
#'   planning window.
#' @param horizon_policy One of `"continuing"`, `"truncate_lives"`, or
#'   `"close_intake"`. Under `"continuing"`, interviews owed after `horizon`
#'   are returned separately in `tail_commitments`. `"truncate_lives"` cuts
#'   those interviews at the horizon, and a life it cuts ends at its last
#'   interview inside the window. `"close_intake"` stops recruitment soon
#'   enough for every recruited cohort to finish by the horizon.
#' @param refreshment Whether later cohorts are drawn from an
#'   `"entrant_register"` or a `"whole_vintage"`. This records a frame role
#'   and does not establish disjointness or combine weights.
#' @param frame_vintage `NULL`, or a named character vector supplying frame
#'   identifiers for any subset of the generated cohort names. Unnamed or
#'   unknown entries are refused rather than matched by position.
#' @param rounding The explicit operational rounding rule. The schema
#'   supports `"ceiling"`, applied to each startup panel and intake cohort
#'   before totals are formed. Operational counts are stored as whole-valued
#'   numerics so plans above R's 32-bit integer range remain representable.
#' @param take `NULL` (default), or a positive whole number of units a PSU
#'   takes per tenure. With it the assignment unit is the PSU. Each panel
#'   issues `ceiling(n_entrants / take)` PSUs, rounded at the PSU, and the
#'   units it stands for are that number times `take`. The unit quantities
#'   stay with the panel plan, whose target carries the clustering through
#'   its `deff`.
#' @param tenures `1` (default) or `"open"`. With `1` a PSU leaves with its
#'   units at the end of their life. With `"open"` the PSUs stay through the
#'   horizon and only the units rotate. Each panel is then a PSU group whose
#'   life repeats the rotation back to back, and `visit` restarts at 1 at
#'   each tenure, where the group's new units are drawn. A finite number of
#'   tenures above 1 is not available.
#' @param subperiods A positive whole number, `1` by default. The rotation,
#'   the panel plan and `horizon` are stated per reporting period, and the
#'   schedule is written in `horizon * subperiods` fieldwork waves. See the
#'   section on sub-periods.
#'
#' @return An `svyplan_schedule` list with scalar planning metadata and four
#'   data frames. `schema_version` is 2. The metadata carries `n_entrants`,
#'   `n_in_sample` and `n_held` from the panel plan, and `n_cohorts`, the
#'   cohorts in sample at an occasion, which is the number of visits.
#'   `life_length` is the life in occasions, gaps included. `take`,
#'   `tenures` (`1`, or `Inf` for `"open"`) and `subperiods` record the
#'   arguments.
#' \describe{
#'   \item{`components`}{One row per execution: startup followed by the intake
#'     cohorts permitted by the horizon policy. Continuous and operational
#'     issue remain separate.}
#'   \item{`schedule`}{The dense component-panel grid over occasions 1 through
#'     `horizon`. `period` names each wave's reporting period and equals
#'     `wave`. A row is in one of three states. Before a panel's entry or
#'     after its life it has `active = FALSE` and `life_stage` and `visit`
#'     `NA`. Held out of sample, in a gap of its life, it has
#'     `active = FALSE`, a `life_stage` and `visit = NA`. Interviewed, it has
#'     `active = TRUE`, a `life_stage` and a `visit`.}
#'   \item{`tail_commitments`}{A sparse table of interviews after `horizon`
#'     owed to components entering within a `"continuing"` window, with the
#'     same `period`, `life_stage` and `visit`. It lists interviews only and
#'     has a fixed empty schema under the other policies.}
#'   \item{`issue`}{The continuous and operational recruitment issued at each
#'     in-horizon occasion, in assignment units, and `tenure_issue`, the units
#'     first interviewed there, which is where a tenure's units are drawn.
#'     [as.data.frame()] returns this view.}
#' }
#'
#' The generated activity is checked against `panel_plan$launch_waves` over
#' every shared period that the horizon policy has not changed. Its realized
#' issued overlap is also checked against the overlap the rotation
#' produces, at every mature membership origin for which the destination
#' occasion exists.
#'
#' @details
#' `life_stage` counts the occasions of a panel's own life and `visit` its
#' interviews, both from 1. They are equal for a life without a break. Under
#' `"2-2-2"` an intake cohort has `life_stage` 1 to 6 and `visit` 1, 2,
#' `NA`, `NA`, 3, 4.
#'
#' An immediate start covers every stage of the life at the first occasion,
#' one startup panel per stage, as [n_panel()] plans it. A panel covering a
#' stage out of sample is selected at the first occasion and held until the
#' rest of its life returns to sample. Its first interview, where `visit` is
#' 1, comes later. Under `"2-2-2"` the panels covering stages 3 and 4 are
#' first interviewed at occasions 3 and 2. A panel numbers its own stages
#' and interviews, so every panel has a `visit` of 1 at its first interview.
#'
#' ## PSUs and tenures
#'
#' With `take`, the schedule counts PSUs. Each panel is a group of
#' `ceiling(n_entrants / take)` PSUs, and a group takes `take` units in each
#' PSU at every tenure. The PSU allocation then reads its total from the
#' schedule, as in `n_alloc(frame, n = schedule$components$operational_issue[1],
#' n_step = G)` with `G` the number of startup panels, so that the master has
#' the PSU count the schedule promises.
#'
#' Under `tenures = "open"` no PSU is replaced inside the horizon. An
#' immediate start selects every group at the first occasion, each group's
#' first tenure being the rest of the life from the stage it covers, and a
#' gradual start brings in one group an occasion until the life is covered.
#' A new tenure is a new draw of units, so the overlap of the units is still
#' the one [design_overlap()] gives for the rotation. The horizon policies
#' apply to tenures. `"continuing"` lists the rest of every tenure opened in
#' the window, `"truncate_lives"` cuts tenures at the horizon, and
#' `"close_intake"` opens no tenure that cannot finish by it.
#'
#' ## Sub-periods
#'
#' Monthly fieldwork for quarterly estimates is a schedule over months. With
#' `subperiods = s`, each reporting panel `g` becomes `s` panels, numbered
#' `(g - 1) * s + j` for `j` in 1 to `s`, so the group of a panel is
#' `ceiling(panel / s)`. Panel `j` is interviewed at wave
#' `(period - 1) * s + j` of every period its group is interviewed in, and is
#' held in the other waves of the group's life. A unit keeps its sub-period at
#' every visit and across tenures. `period` names the reporting period of
#' each wave, and `life_stage`, `life_length` and `visit` count periods and
#' interviews, not waves. `horizon`, `wave`, `entry_wave` and
#' `steady_state_from` count waves, and a cohort entering in period `t`
#' enters at wave `(t - 1) * s + 1`, every sub-period panel being selected
#' then. Each sub-period panel takes `ceiling(n_entrants / s)` units, or that
#' over `take` in PSUs, so a period issues up to `s - 1` more per group than
#' the continuous plan. The launch and overlap checks run on the reporting
#' grid, and `$overlap` is the reporting-period profile.
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
#'   horizon = 6,
#'   horizon_policy = "continuing",
#'   refreshment = "entrant_register",
#'   rounding = "ceiling"
#' )
#' schedule
#' schedule$components
#'
#' # A 2-2-2 rotation: six startup panels, four interviewed at each occasion
#' quarterly <- n_panel(
#'   n_prop(p = 0.5, moe = 0.03),
#'   retention = c(0.9, 0.85, 0.9),
#'   resp_rate = 0.75,
#'   design = "rotating",
#'   start = "immediate",
#'   rotation = "2-2-2"
#' )
#' gapped <- design_schedule(
#'   quarterly,
#'   horizon = 8,
#'   horizon_policy = "continuing",
#'   refreshment = "entrant_register",
#'   rounding = "ceiling"
#' )
#' gapped
#' startup <- gapped$schedule[gapped$schedule$cohort == "startup", ]
#' xtabs(visit ~ panel + wave, startup, addNA = TRUE)
#'
#' # Quarterly reporting on a life of five quarters, 14 units per PSU, PSUs
#' # kept and monthly fieldwork: 15 panels, each in its own month
#' monthly <- design_schedule(
#'   n_panel(n_prop(p = 0.5, moe = 0.02, deff = 1.8),
#'           retention = rep(0.9, 4), resp_rate = 0.8,
#'           design = "rotating", start = "immediate"),
#'   horizon = 4,
#'   horizon_policy = "continuing",
#'   refreshment = "entrant_register",
#'   rounding = "ceiling",
#'   take = 14,
#'   tenures = "open",
#'   subperiods = 3
#' )
#' monthly
#' xtabs(visit ~ panel + wave, monthly$schedule, addNA = TRUE)
#'
#' @seealso [design_rotation()] for the pattern this schedules, [n_panel()]
#'   and [prec_panel()] for the rotating plan, [plot.svyplan_overlap()] for
#'   the rotation chart, and [print.svyplan_schedule()] for the printed form
#'   and its `summary()`.
#'
#' @family repeated survey functions
#' @export
design_schedule <- function(
  panel_plan,
  horizon,
  horizon_policy,
  refreshment,
  frame_vintage = NULL,
  rounding,
  take = NULL,
  tenures = 1,
  subperiods = 1
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
  rotation <- panel_plan$rotation
  if (!inherits(rotation, "svyplan_rotation")) {
    stop(
      "'panel_plan' must carry its rotation; plan it with n_panel() or prec_panel()",
      call. = FALSE
    )
  }
  w <- as.double(rotation)
  life_length <- attr(rotation, "life", exact = TRUE)

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
  if (identical(horizon_policy, "truncate_lives") &&
      identical(panel_plan$start, "immediate")) {
    first_interview <- max(vapply(seq_len(life_length), function(s) {
      which(w[s:life_length] > 0)[[1L]]
    }, integer(1L)))
    if (horizon < first_interview) {
      stop(
        sprintf(
          "under \"truncate_lives\" the horizon has to reach the first interview of every startup panel, which the last of them has at occasion %d",
          first_interview
        ),
        call. = FALSE
      )
    }
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
  take <- .schedule_take(take)
  tenures <- .schedule_tenures(tenures)
  subperiods <- .schedule_whole(subperiods, "subperiods")

  entrant_issue <- unname(panel_plan$n_entrants)
  plan <- .schedule_plan(
    panel_plan$start, life_length, horizon, horizon_policy, tenures
  )
  per_panel <- entrant_issue / (.schedule_unit(take) * subperiods)
  # Whole-valued doubles remain exact far beyond R's 32-bit integer range.
  # Operational counts are counts by value, not storage-mode integers.
  panel_issue <- ceiling(per_panel)
  panels <- plan$panels * subperiods
  components <- data.frame(
    cohort = plan$cohort,
    entry_wave = as.integer((plan$entry_wave - 1L) * subperiods + 1L),
    frame_role = c("startup", rep(refreshment, nrow(plan) - 1L)),
    planned_issue = per_panel * panels,
    operational_issue = panel_issue * panels,
    panels = as.integer(panels),
    panel_issue = rep(panel_issue, nrow(plan)),
    frame_vintage = .schedule_frame_vintages(frame_vintage, plan$cohort),
    status = rep("planned", nrow(plan)),
    stringsAsFactors = FALSE
  )

  tables <- .schedule_tables(
    plan, panel_plan$start, w, horizon, horizon_policy, tenures, subperiods
  )
  overlap <- design_overlap(rotation)
  reporting <- list(
    schedule = tables$reporting,
    n_entrants = entrant_issue,
    response_model = list(q = panel_plan$waves$q),
    horizon = horizon,
    life_length = life_length,
    life_spec = w
  )
  last_open <- if (identical(horizon_policy, "close_intake")) {
    horizon - life_length + 1L
  } else {
    horizon
  }
  .check_schedule_launch(
    reporting, panel_plan$launch_waves, horizon_policy, last_open
  )
  .check_schedule_overlap(reporting, as.double(overlap))

  schedule <- tables$schedule
  steady_waves <- which(tables$steady)
  steady_state_from <- if (length(steady_waves) == 0L) {
    NA_integer_
  } else {
    as.integer(steady_waves[[1L]])
  }

  .new_svyplan_schedule(
    schema_version = 2L,
    life_spec = as.numeric(w),
    life_notation = "resolved_per_occasion_take",
    life_length = life_length,
    horizon = horizon * subperiods,
    horizon_policy = horizon_policy,
    launch_policy = panel_plan$start,
    target_estimand = list(type = panel_plan$type, method = panel_plan$method),
    n_target = unname(panel_plan$n_target),
    response_model = list(
      resp_rate = panel_plan$params$resp_rate,
      retention = panel_plan$params$retention,
      q = panel_plan$waves$q
    ),
    rounding = list(
      rule = rounding, level = .schedule_rounding_level(take, subperiods)
    ),
    take = take,
    tenures = tenures,
    subperiods = subperiods,
    steady_state_from = steady_state_from,
    overlap = as.data.frame(overlap),
    refreshment = refreshment,
    requires_disjoint_entrants = identical(refreshment, "entrant_register"),
    n_entrants = entrant_issue,
    n_in_sample = unname(panel_plan$n_in_sample),
    n_held = unname(panel_plan$n_held),
    n_cohorts = as.integer(panel_plan$n_cohorts),
    panel_parameters = .schedule_panel_parameters(
      schedule, tables$tail, components$panels[[1L]]
    ),
    components = components,
    schedule = schedule,
    tail_commitments = tables$tail,
    issue = .schedule_issue_profile(
      components, schedule, tables$steady, .schedule_unit(take)
    )
  )
}

#' @keywords internal
#' @noRd
.schedule_take <- function(take) {
  if (is.null(take)) {
    return(NA_real_)
  }
  if (!is.numeric(take) || length(take) != 1L || is.na(take) ||
      !is.finite(take) || take < 1 || take != trunc(take)) {
    stop(
      "'take' must be NULL or a positive whole number of units per PSU per tenure",
      call. = FALSE
    )
  }
  as.double(take)
}

#' @keywords internal
#' @noRd
.schedule_tenures <- function(tenures) {
  if (identical(tenures, "open")) {
    return(Inf)
  }
  if (is.numeric(tenures) && length(tenures) == 1L && !is.na(tenures) &&
      tenures == 1) {
    return(1)
  }
  stop(
    "'tenures' must be 1, a PSU leaving with its units, or \"open\", PSUs kept through the horizon; a finite number of tenures above 1 is not described yet",
    call. = FALSE
  )
}

#' Units an assignment unit stands for, 1 without a take
#' @keywords internal
#' @noRd
.schedule_unit <- function(take) {
  if (is.na(take)) 1 else take
}

#' @keywords internal
#' @noRd
.schedule_rounding_level <- function(take, subperiods) {
  parts <- c(if (!is.na(take)) "psu", if (subperiods > 1L) "subperiod")
  if (length(parts) == 0L) "panel_and_cohort" else paste(parts, collapse = "_and_")
}

#' The executions of a schedule on the reporting grid
#' @keywords internal
#' @noRd
.schedule_plan <- function(start, life_length, horizon, policy, tenures) {
  max_entry <- if (identical(policy, "close_intake")) {
    horizon - life_length + 1L
  } else {
    horizon
  }
  if (is.infinite(tenures)) {
    max_entry <- if (identical(start, "immediate")) {
      1L
    } else {
      min(life_length, max_entry)
    }
  }
  max_entry <- as.integer(max_entry)
  data.frame(
    cohort = c("startup", if (max_entry >= 2L) paste0("intake_", 2:max_entry)),
    entry_wave = seq_len(max_entry),
    panels = c(
      if (identical(start, "immediate")) life_length else 1L,
      rep(1L, max_entry - 1L)
    ),
    stringsAsFactors = FALSE
  )
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
.schedule_panel_catalog <- function(plan, start, w, horizon,
                                    horizon_policy, tenures) {
  life <- length(w)
  launch <- .rotation_cohorts(w, start, 1L)
  startup_first <- rev(life - lengths(lapply(launch, `[[`, "w")) + 1L)
  rows <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    k <- plan$panels[[i]]
    first <- if (i == 1L) startup_first else rep(1L, k)
    h <- horizon - plan$entry_wave[[i]] + 1L
    len <- vapply(first, function(f) {
      .schedule_life_length(w, f, h, horizon_policy, tenures)
    }, integer(1L))
    rows[[i]] <- data.frame(
      cohort = rep(plan$cohort[[i]], k),
      panel = seq_len(k),
      entry_wave = rep(plan$entry_wave[[i]], k),
      first_stage = as.integer(first),
      life_length = len,
      relative_take = rep(1, k),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Where a panel's own stage sits in the life, tenures repeating it
#' @keywords internal
#' @noRd
.schedule_position <- function(w, first, stage) {
  life <- length(w)
  first_tenure <- life - first + 1L
  later <- stage > first_tenure
  in_tenure <- ifelse(later, (stage - first_tenure - 1L) %% life + 1L, stage)
  list(
    pos = ifelse(later, in_tenure, first + stage - 1L),
    start = ifelse(later, 1L, first),
    in_tenure = in_tenure
  )
}

#' Occasions of a panel's own life inside and beyond the window
#' @keywords internal
#' @noRd
.schedule_life_length <- function(w, first, h, policy, tenures) {
  life <- length(w)
  first_tenure <- life - first + 1L
  if (identical(policy, "truncate_lives")) {
    cut <- if (is.finite(tenures)) min(first_tenure, h) else h
    on <- w[.schedule_position(w, first, seq_len(cut))$pos] > 0
    return(as.integer(max(which(on))))
  }
  if (is.finite(tenures) || h <= first_tenure) {
    return(as.integer(first_tenure))
  }
  more <- if (identical(policy, "close_intake")) {
    (h - first_tenure) %/% life
  } else {
    ceiling((h - first_tenure) / life)
  }
  as.integer(first_tenure + more * life)
}

#' Column order of the dense schedule
#' @keywords internal
#' @noRd
.schedule_columns <- function() {
  c("wave", "period", "cohort", "panel", "life_stage", "visit",
    "life_length", "active", "steady_state", "relative_take")
}

#' Every catalog panel at one wave: held, interviewed, and its visit
#'
#' `.tstage` is the stage within the current tenure, which the steady test
#' counts and the stored tables drop.
#' @keywords internal
#' @noRd
.schedule_rows <- function(catalog, wave, w) {
  seen <- c(0L, cumsum(w > 0))
  stage <- wave - catalog$entry_wave + 1L
  held <- stage >= 1L & stage <= catalog$life_length
  at <- .schedule_position(w, catalog$first_stage, ifelse(held, stage, 1L))
  active <- held & w[at$pos] > 0
  data.frame(
    wave = rep(as.integer(wave), nrow(catalog)),
    period = rep(as.integer(wave), nrow(catalog)),
    cohort = catalog$cohort,
    panel = as.integer(catalog$panel),
    life_stage = ifelse(held, as.integer(stage), NA_integer_),
    visit = ifelse(
      active, as.integer(seen[at$pos + 1L] - seen[at$start]), NA_integer_
    ),
    life_length = as.integer(catalog$life_length),
    active = active,
    relative_take = catalog$relative_take,
    .tstage = ifelse(held, as.integer(at$in_tenure), NA_integer_),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.schedule_dense_grid <- function(catalog, horizon, w) {
  out <- do.call(rbind, lapply(seq_len(horizon), function(wave) {
    .schedule_rows(catalog, wave, w)
  }))
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.empty_schedule_tail <- function() {
  data.frame(
    wave = integer(),
    period = integer(),
    cohort = character(),
    panel = integer(),
    life_stage = integer(),
    visit = integer(),
    life_length = integer(),
    relative_take = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Interviews owed after the horizon, and only interviews
#' @keywords internal
#' @noRd
.schedule_tail <- function(catalog, horizon, policy, w) {
  last <- max(catalog$entry_wave + catalog$life_length - 1L)
  if (!identical(policy, "continuing") || last <= horizon) {
    return(.empty_schedule_tail())
  }
  out <- do.call(rbind, lapply(seq.int(horizon + 1L, last), function(wave) {
    .schedule_rows(catalog, wave, w)
  }))
  out <- out[out$active, names(.empty_schedule_tail()), drop = FALSE]
  if (nrow(out) == 0L) {
    return(.empty_schedule_tail())
  }
  rownames(out) <- NULL
  out
}

#' Occasions where every stage of a tenure is held exactly once
#' @keywords internal
#' @noRd
.schedule_steady <- function(schedule, life_length) {
  vapply(seq_len(max(schedule$wave)), function(wave) {
    at <- schedule[
      schedule$wave == wave & !is.na(schedule$.tstage),
      c(".tstage", "relative_take"),
      drop = FALSE
    ]
    totals <- numeric(life_length)
    if (nrow(at) > 0L) {
      by_stage <- tapply(at$relative_take, at$.tstage, sum)
      totals[as.integer(names(by_stage))] <- as.numeric(by_stage)
    }
    isTRUE(all.equal(totals, rep(1, life_length), tolerance = 1e-12))
  }, logical(1L))
}

#' The reporting grid, and the fieldwork tables split from it
#' @keywords internal
#' @noRd
.schedule_tables <- function(plan, start, w, horizon, policy, tenures,
                             subperiods) {
  catalog <- .schedule_panel_catalog(plan, start, w, horizon, policy, tenures)
  grid <- .schedule_dense_grid(catalog, horizon, w)
  tail <- .schedule_tail(catalog, horizon, policy, w)
  steady <- .schedule_steady(grid, length(w))
  grid$steady_state <- steady[grid$wave]
  grid <- grid[.schedule_columns()]
  if (subperiods == 1L) {
    return(list(reporting = grid, schedule = grid, tail = tail,
                steady = steady))
  }
  list(
    reporting = grid,
    schedule = .schedule_split(grid, subperiods, dense = TRUE),
    tail = .schedule_split(tail, subperiods, dense = FALSE),
    steady = rep(steady, each = subperiods)
  )
}

#' Split each reporting panel into its sub-period panels
#'
#' Panel `(g - 1) s + j` is interviewed at wave `(t - 1) s + j` of every
#' period `t` its group is interviewed in, and held in the others.
#' @keywords internal
#' @noRd
.schedule_split <- function(tab, s, dense) {
  if (nrow(tab) == 0L) {
    return(tab)
  }
  rank <- match(tab$cohort, unique(tab$cohort))
  combos <- if (dense) {
    expand.grid(j = seq_len(s), m = seq_len(s))
  } else {
    data.frame(j = seq_len(s), m = seq_len(s))
  }
  row <- rep(seq_len(nrow(tab)), each = nrow(combos))
  j <- rep(combos$j, times = nrow(tab))
  m <- rep(combos$m, times = nrow(tab))
  out <- tab[row, , drop = FALSE]
  out$wave <- as.integer((tab$period[row] - 1L) * s + m)
  out$panel <- as.integer((tab$panel[row] - 1L) * s + j)
  if (dense) {
    out$active <- tab$active[row] & m == j
    out$visit <- ifelse(out$active, tab$visit[row], NA_integer_)
  }
  out <- out[order(out$wave, rank[row], out$panel), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.schedule_issue_profile <- function(components, schedule, steady, unit) {
  opening <- schedule$active & !is.na(schedule$visit) & schedule$visit == 1L
  opened <- tabulate(schedule$wave[opening], length(steady))
  rows <- lapply(seq_along(steady), function(wave) {
    at <- components[components$entry_wave == wave, , drop = FALSE]
    data.frame(
      wave = wave,
      cohort = if (nrow(at) == 0L) NA_character_ else at$cohort[[1L]],
      planned_issue = if (nrow(at) == 0L) 0 else at$planned_issue[[1L]],
      operational_issue = if (nrow(at) == 0L) 0 else at$operational_issue[[1L]],
      tenure_issue = opened[[wave]] * components$panel_issue[[1L]] * unit,
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
  per_wave <- table(c(startup$wave, tail$wave[tail$cohort == "startup"]))
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
    "n_target", "response_model", "rounding", "take", "tenures",
    "subperiods", "steady_state_from",
    "overlap", "refreshment", "requires_disjoint_entrants", "n_entrants",
    "n_in_sample", "n_held", "n_cohorts", "panel_parameters", "components",
    "schedule", "tail_commitments", "issue"
  )
  if (!is.list(x) || !identical(names(x), required)) {
    fail("top-level fields do not match schema version 2")
  }
  if (!identical(x$schema_version, 2L) ||
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
    fail("scalar metadata does not match schema version 2")
  }
  if (!is.numeric(x$life_spec) || length(x$life_spec) != x$life_length ||
      anyNA(x$life_spec) || any(!is.finite(x$life_spec)) ||
      any(x$life_spec < 0) || x$life_spec[[1L]] <= 0 ||
      x$life_spec[[x$life_length]] <= 0 ||
      any(x$life_spec[x$life_spec > 0] != 1)) {
    fail("life_spec must be an equal-take life that starts and ends in sample")
  }
  visits <- as.integer(sum(x$life_spec > 0))
  if (identical(x$horizon_policy, "close_intake") &&
      x$horizon %/% x$subperiods < x$life_length) {
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
      length(response$retention) != visits - 1L ||
      anyNA(response$retention) || any(!is.finite(response$retention)) ||
      any(response$retention <= 0 | response$retention > 1) ||
      !is.numeric(response$q) || length(response$q) != visits ||
      anyNA(response$q) || any(!is.finite(response$q))) {
    fail("response metadata is inconsistent")
  }
  expected_q <- response$resp_rate * cumprod(c(1, response$retention))
  if (!same(as.numeric(response$q), expected_q)) {
    fail("response probabilities do not follow resp_rate and retention")
  }
  if (!is.double(x$take) || length(x$take) != 1L ||
      !(is.na(x$take) || (is.finite(x$take) && x$take >= 1 &&
                            x$take == trunc(x$take))) ||
      !identical(x$tenures, 1) && !identical(x$tenures, Inf) ||
      !is.integer(x$subperiods) || length(x$subperiods) != 1L ||
      is.na(x$subperiods) || x$subperiods < 1L ||
      x$horizon %% x$subperiods != 0L) {
    fail("take, tenures or subperiods metadata is inconsistent")
  }
  s <- x$subperiods
  horizon <- x$horizon %/% s
  unit <- .schedule_unit(x$take)
  if (!is.list(x$rounding) ||
      !identical(x$rounding, list(
        rule = "ceiling", level = .schedule_rounding_level(x$take, s)
      ))) {
    fail("rounding metadata does not match schema version 2")
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
      !scalar_number(x$n_held, lower = 0) ||
      !identical(x$n_cohorts, visits) ||
      !same(x$n_in_sample, x$n_cohorts * x$n_entrants) ||
      !same(x$n_held, x$life_length * x$n_entrants)) {
    fail("panel counts are inconsistent")
  }

  component_cols <- c(
    "cohort", "entry_wave", "frame_role", "planned_issue",
    "operational_issue", "panels", "panel_issue", "frame_vintage", "status"
  )
  schedule_cols <- .schedule_columns()
  tail_cols <- names(.empty_schedule_tail())
  issue_cols <- c(
    "wave", "cohort", "planned_issue", "operational_issue", "tenure_issue",
    "steady_state"
  )
  exact_table <- function(tab, cols) {
    is.data.frame(tab) && identical(names(tab), cols)
  }
  if (!exact_table(x$components, component_cols) ||
      !exact_table(x$schedule, schedule_cols) ||
      !exact_table(x$tail_commitments, tail_cols) ||
      !exact_table(x$issue, issue_cols)) {
    fail("one or more tables do not have the fixed schema-version-2 columns")
  }

  components <- x$components
  plan <- .schedule_plan(
    x$launch_policy, x$life_length, horizon, x$horizon_policy, x$tenures
  )
  expected_panels <- as.integer(plan$panels * s)
  expected_roles <- c("startup", rep(x$refreshment, nrow(plan) - 1L))
  if (!identical(components$cohort, plan$cohort) ||
      anyNA(components$cohort) || anyDuplicated(components$cohort) ||
      !is.integer(components$entry_wave) || anyNA(components$entry_wave) ||
      !identical(components$entry_wave,
                 as.integer((plan$entry_wave - 1L) * s + 1L)) ||
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
  per_panel <- x$n_entrants / (unit * s)
  expected_panel_issue <- rep(ceiling(per_panel), nrow(components))
  if (!same(components$planned_issue, per_panel * expected_panels) ||
      !same(components$panel_issue, expected_panel_issue) ||
      !same(
        components$operational_issue,
        expected_panel_issue * expected_panels
      )) {
    fail("component issue does not follow panel-level ceiling")
  }

  schedule <- x$schedule
  n_panels <- sum(components$panels)
  bounded <- is.finite(x$tenures)
  if (!is.integer(schedule$wave) || !is.integer(schedule$period) ||
      !is.integer(schedule$panel) ||
      !is.integer(schedule$life_stage) || !is.integer(schedule$visit) ||
      !is.integer(schedule$life_length) ||
      !is.logical(schedule$active) || !is.logical(schedule$steady_state) ||
      !is.numeric(schedule$relative_take) ||
      nrow(schedule) != x$horizon * n_panels ||
      anyNA(schedule[c("wave", "period", "cohort", "panel", "life_length",
                       "active", "steady_state", "relative_take")]) ||
      any(schedule$period != (schedule$wave - 1L) %/% s + 1L) ||
      any(!schedule$cohort %in% components$cohort) ||
      any(schedule$wave < 1L | schedule$wave > x$horizon) ||
      any(schedule$panel < 1L) ||
      any(schedule$life_length < 1L) ||
      bounded && any(schedule$life_length > x$life_length) ||
      any(schedule$relative_take != 1) ||
      any(is.na(schedule$visit) != !schedule$active) ||
      any(is.na(schedule$life_stage) & schedule$active) ||
      any(!is.na(schedule$life_stage) &
            (schedule$life_stage < 1L |
               schedule$life_stage > schedule$life_length))) {
    fail("the dense schedule grid is inconsistent")
  }
  grid_key <- paste(schedule$wave, schedule$cohort, schedule$panel, sep = "|")
  if (anyDuplicated(grid_key)) {
    fail("the dense schedule grid contains duplicate component-panel rows")
  }
  interviews <- rbind(
    schedule[schedule$active, c("cohort", "panel", "wave", "visit")],
    x$tail_commitments[c("cohort", "panel", "wave", "visit")]
  )
  interviews <- interviews[order(interviews$cohort, interviews$panel,
                                 interviews$wave), ]
  if (anyNA(interviews$visit)) {
    fail("an interview has no visit")
  }
  step <- stats::ave(interviews$visit, interviews$cohort, interviews$panel,
                     FUN = function(v) c(v[[1L]], diff(v)))
  restart <- interviews$visit == 1L
  first <- !duplicated(interviews[c("cohort", "panel")])
  if (any(first & !restart) ||
      any(!first & step != 1L & !(restart & !bounded))) {
    fail("a visit does not count the interviews of its tenure")
  }

  tables <- .schedule_tables(
    plan, x$launch_policy, x$life_spec, horizon, x$horizon_policy,
    x$tenures, s
  )
  if (!identical(schedule, tables$schedule)) {
    fail("the dense schedule does not match the component manifest")
  }

  tail <- x$tail_commitments
  if (!is.integer(tail$wave) || !is.integer(tail$period) ||
      !is.integer(tail$panel) ||
      !is.integer(tail$life_stage) || !is.integer(tail$visit) ||
      !is.integer(tail$life_length) ||
      any(tail$period != (tail$wave - 1L) %/% s + 1L) ||
      !is.numeric(tail$relative_take) || anyNA(tail) ||
      any(tail$wave <= x$horizon) ||
      any(!tail$cohort %in% components$cohort) ||
      any(tail$panel < 1L) ||
      any(tail$life_stage < 1L | tail$life_stage > tail$life_length) ||
      any(tail$life_length < 1L) ||
      bounded && any(tail$life_length > x$life_length) ||
      any(tail$relative_take != 1) ||
      (!identical(x$horizon_policy, "continuing") && nrow(tail) > 0L)) {
    fail("the tail-commitment table is inconsistent with the horizon policy")
  }
  if (!identical(tail, tables$tail)) {
    fail("tail commitments do not match the component manifest")
  }

  issue <- x$issue
  if (!is.integer(issue$wave) ||
      !identical(issue$wave, seq_len(x$horizon)) ||
      !is.numeric(issue$planned_issue) ||
      !whole_counts(issue$operational_issue) ||
      !whole_counts(issue$tenure_issue) ||
      !is.logical(issue$steady_state) || anyNA(issue$steady_state) ||
      any(issue$planned_issue < 0) ||
      any(issue$operational_issue < issue$planned_issue)) {
    fail("the occasion issue profile is inconsistent")
  }
  expected_issue <- .schedule_issue_profile(
    components, tables$schedule, tables$steady, unit
  )
  if (!identical(issue, expected_issue)) {
    fail("the occasion issue profile does not match the component manifest")
  }

  expected_steady_waves <- which(tables$steady)
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
    tables$schedule, tables$tail, expected_panels[[1L]]
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
    wave = tab$visit,
    n_issued = x$n_entrants * tab$relative_take,
    n_resp = x$n_entrants * tab$relative_take * q[tab$visit]
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
  active <- active[order(active$cohort, active$panel, active$wave), ]
  tenure <- stats::ave(as.integer(active$visit == 1L), active$cohort,
                       active$panel, FUN = cumsum)
  key <- paste(active$cohort, active$panel, tenure, sep = "|")
  take <- vapply(seq_len(x$horizon), function(w) {
    sum(active$relative_take[active$wave == w])
  }, numeric(1L))
  for (lag in seq_len(min(length(expected), x$horizon - 1L))) {
    origins <- which(
      seq_len(x$horizon) + lag <= x$horizon &
        abs(take - sum(x$life_spec > 0)) <= 1e-10
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

#' Print, summarize and coerce a longitudinal design schedule
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
#'
#' @seealso [design_schedule()], which builds these objects, [n_panel()] for
#'   the rotating plan they schedule, and [print.svyplan_overlap()] for the
#'   rotation the life describes.
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
#'   horizon = 6,
#'   horizon_policy = "continuing",
#'   refreshment = "entrant_register",
#'   rounding = "ceiling"
#' )
#' schedule
#'
#' # the occasion-by-occasion tables the print block summarizes
#' s <- summary(schedule)
#' s
#'
#' # every table stays reachable as a field
#' schedule$components
#' schedule$tail_commitments
#'
#' # as.data.frame() returns the issue profile
#' as.data.frame(schedule)
#'
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
  cat(.fmt_schedule_rotation(x))
  cat(.fmt_schedule_psu(x))
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
  steady <- if (is.na(x$steady_state_from)) {
    "steady state not reached in the window"
  } else {
    sprintf("steady from %s %d", .fmt_schedule_occasion(x),
            x$steady_state_from)
  }
  if (x$subperiods == 1L) {
    return(sprintf(
      "life: %d stages over %d occasions, %s\n", x$life_length, x$horizon,
      steady
    ))
  }
  sprintf(
    "life: %d periods over %d periods of %d waves, %s\n", x$life_length,
    x$horizon %/% x$subperiods, x$subperiods, steady
  )
}

#' What one row of the schedule's grid is called
#' @keywords internal
#' @noRd
.fmt_schedule_occasion <- function(x) {
  if (x$subperiods == 1L) "occasion" else "wave"
}

#' The PSU take and the PSU stay, where either was declared
#' @keywords internal
#' @noRd
.fmt_schedule_psu <- function(x) {
  paste0(
    if (!is.na(x$take)) {
      sprintf("take: %s units per PSU per tenure\n", format(x$take))
    },
    if (is.infinite(x$tenures)) "tenures: open, PSUs kept through the horizon\n"
  )
}

#' The rotation of a life with a gap
#' @keywords internal
#' @noRd
.fmt_schedule_rotation <- function(x) {
  if (all(x$life_spec > 0)) {
    return("")
  }
  sprintf("rotation: %s\n", .fmt_rotation(x$life_spec))
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
  unit <- "occasion"
  if (x$subperiods > 1L) {
    v <- as.vector(tapply(v, (x$issue$wave - 1L) %/% x$subperiods, sum))
    unit <- "period"
  }
  psu <- if (is.na(x$take)) "" else " PSUs"
  closed <- if (is.infinite(x$tenures)) "no PSU intake from" else {
    "intake closed from"
  }
  parts <- sprintf("%s%s at startup", format(v[[1L]]), psu)
  if (length(v) > 1L) {
    runs <- rle(v[-1L])
    at <- 2L
    for (j in seq_along(runs$lengths)) {
      lo <- at
      hi <- at + runs$lengths[[j]] - 1L
      at <- hi + 1L
      parts <- c(parts, if (runs$values[[j]] == 0) {
        if (x$subperiods == 1L) sprintf("%s %d", closed, lo) else {
          sprintf("%s %s %d", closed, unit, lo)
        }
      } else {
        sprintf(
          "%s%s per %s (%s)", format(runs$values[[j]]), psu, unit,
          if (lo == hi) sprintf("%s %d", unit, lo) else {
            sprintf("%ss %d-%d", unit, lo, hi)
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
    "tail commitments: %d panel-interviews after %s %d\n",
    nrow(x$tail_commitments), .fmt_schedule_occasion(x), x$horizon
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
  cat(.fmt_schedule_rotation(s))
  cat(.fmt_schedule_psu(s))
  cat(sprintf(
    "rounding: %s at %s level\n",
    s$rounding$rule, gsub("_", " ", s$rounding$level, fixed = TRUE)
  ))
  # From the rounded panels, which ceiling(n_entrants) understates.
  entrants <- s$components$panel_issue[[1L]] * .schedule_unit(s$take) *
    s$subperiods
  cat(sprintf(
    "entrants: %s per %s, %s in sample across %d cohorts\n",
    format(entrants), if (s$subperiods == 1L) "occasion" else "period",
    format(entrants * s$n_cohorts), s$n_cohorts
  ))
  if (s$life_length > s$n_cohorts) {
    cat(sprintf(
      "held: %s across %d cohorts, %d of them out of sample\n",
      format(entrants * s$life_length), s$life_length,
      s$life_length - s$n_cohorts
    ))
  }
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
