#' Construct an operational schedule from a rotating panel plan
#'
#' Translate the start-up and recruitment quantities in a rotating
#' [svyplan_panel][n_panel()] into named selections, a dense panel-by-period
#' schedule, and any interview commitments beyond the planning horizon. The
#' result is a planning object. It does not inspect a frame, execute a
#' sample, or construct analysis weights.
#'
#' @param panel_plan A rotating `svyplan_panel` from [n_panel()] or
#'   [prec_panel()] with `start = "immediate"` or `start = "gradual"`. Its
#'   `rotation` is the life the schedule lays out, and the overlap the
#'   manifest records is computed from it, so neither is supplied here.
#' @param horizon A positive whole number of reporting periods in the finite
#'   planning window.
#' @param horizon_policy One of `"continuing"`, `"truncate_lives"`, or
#'   `"close_intake"`. Under `"continuing"`, interviews owed after `horizon`
#'   are returned separately in `tail_commitments`. `"truncate_lives"` cuts
#'   those interviews at the horizon, and a life it cuts ends at its last
#'   interview inside the window. `"close_intake"` stops recruitment soon
#'   enough for every recruited cohort to finish by the horizon.
#' @param refreshment Whether later selections are drawn from an
#'   `"entrant_register"` or a `"whole_vintage"`. This records a frame role
#'   and does not establish disjointness or combine weights.
#' @param frame_vintage `NULL`, or a named character vector supplying frame
#'   identifiers for any subset of the generated selection names. Unnamed or
#'   unknown entries are refused rather than matched by position.
#' @param rounding The operational rounding rule, `"ceiling"` (default and
#'   the one rule supported), applied to each panel before totals are formed
#'   and recorded in `$rounding`. Operational counts are stored as
#'   whole-valued numerics so plans above R's 32-bit integer range remain
#'   representable. With a finite population the rounded schedule must not
#'   hold more frame units at any period than the target's `N`, out-of-sample
#'   spells included, and one that does is refused.
#' @param take `NULL` (default), or a positive whole number of units a PSU
#'   takes per cohort. With it the assignment unit is the PSU. Each panel
#'   issues `ceiling(n_entrants / take)` PSUs, rounded at the PSU, and the
#'   units it stands for are that number times `take`. The unit quantities
#'   stay with the panel plan, whose target carries the clustering through
#'   its `deff`.
#' @param retain_psus `FALSE` (default) or `TRUE`, and `TRUE` needs `take`.
#'   With `FALSE` a PSU leaves with its cohort at the end of the cohort's
#'   life. With `TRUE` the PSUs are retained through the horizon and only
#'   their households rotate (Kalton, 2009, p. 95, "may retain the same set
#'   of PSUs"). Each panel is then a PSU group that takes a new cohort each
#'   time its life returns to the start of the pattern. See the section on
#'   retained PSUs.
#' @param entrants Which entrant count to field when the plan carries an
#'   assurance, `"expected"` for `$n_entrants` or `"assured"` for
#'   `$n_assured`. It must be declared for such a plan, since the expected
#'   count does not keep the assurance. Without an assurance it is
#'   `"expected"` and may be left `NULL`. Recorded as `$entrant_basis`.
#' @param subrounds A positive whole number, `1` by default, the subrounds
#'   each period's fieldwork is divided into. The rotation, the panel plan
#'   and `horizon` are stated per reporting period. See the section on
#'   subrounds.
#'
#' @return An `svyplan_schedule` list with scalar planning metadata and four
#'   data frames. The metadata carries `start`, the
#'   panel plan's start, `n_entrants`, `n_in_sample` and `n_held` from the
#'   panel plan, at its assured count when `entrant_basis` is `"assured"`,
#'   and `n_cohorts`, the cohorts in sample at a period, which is
#'   the number of visits. `life_length` is the life in periods, gaps
#'   included, and `horizon` counts periods. `equilibrium_from` is the first
#'   period at equilibrium, `NA` when the window ends before it. `take`,
#'   `retain_psus` and `subrounds` record the arguments.
#' \describe{
#'   \item{`selections`}{One row per selection from the frame: `startup` at
#'     period 1, then one named `period_k` for each later period the horizon
#'     policy allows, each with its `entry_period`, frame role, number of
#'     panels and issue. Continuous and operational issue remain separate.}
#'   \item{`schedule`}{The dense panel grid over periods 1 through `horizon`,
#'     one row per panel at each `period` and `subround`, with the columns
#'     described under Details. A row is in one of three states. Before a
#'     panel's entry or after its life it has `active = FALSE` and `cohort`,
#'     `wave`, `visit` and `life_stage` `NA`. Held out of sample, in a gap of
#'     its life, it has `active = FALSE`, a `cohort` and a `life_stage`, and
#'     `wave` and `visit` `NA`. Interviewed, it has `active = TRUE` and all
#'     four.}
#'   \item{`tail_commitments`}{A sparse table of interviews after `horizon`
#'     owed to panels selected within a `"continuing"` window, with the
#'     columns of `schedule` other than `active` and `equilibrium`. It lists
#'     interviews only and has fixed empty columns under the other
#'     policies.}
#'   \item{`issue`}{One row per subround of each period: the `selection`
#'     made there, its continuous and operational issue in assignment units,
#'     `cohort_issue`, the units of the cohorts first interviewed there,
#'     which is where their units are drawn, and `equilibrium`.
#'     [as.data.frame()] returns this view.}
#' }
#'
#' The generated activity is checked against `panel_plan$period_visits`
#' over every shared period that the horizon policy has not changed. Its
#' realized issued overlap is also checked against the overlap the rotation
#' produces, at every period with the full sample in the field for which
#' the destination period exists.
#'
#' @details
#' The schedule describes each interview in the terms the cell table of
#' [plot.svyplan_rotation()] uses, and for the same rotation and start the
#' two agree row for row.
#'
#' - `selection` names the selection from the frame a panel comes from. A
#'   selection is one draw of first-stage units, or of sample units when
#'   there is no `take`, made at one period. The households a retained PSU
#'   takes later are a new cohort, not a new selection.
#' - `cohort` says whose interview the row is. A cohort is the units at the
#'   same wave at every period, the units whose life in the pattern starts
#'   at the same period, as in the cohorts of the UK LFS (Office for
#'   National Statistics, 2026, section 3). A start-up cohort's start is
#'   counted back before period 1, so the six start-up cohorts of `"2-2-2"`
#'   are distinct, each at its own place in the pattern. Cohorts are
#'   numbered in the chart's order.
#' - `panel` is the unit of assignment to the rotation, the sample units or,
#'   with `take`, a PSU group. Panels are numbered across the whole schedule
#'   as the chart numbers them, the start-up panels of an immediate start
#'   from the one with a single period left to the one with the whole life
#'   (Steel and McLaren, 2008).
#' - `rotation_group` is the slot a new panel takes over when the panel
#'   before it rotates out, assigned by the rule on
#'   [plot.svyplan_rotation()].
#' - `wave` is the place of the interview in the rotation pattern, the
#'   meaning of Regulation (EU) 2019/2240 and of the UK LFS. `visit` is the
#'   interview count of the row's cohort, which is what attrition follows.
#'   They differ only for the start-up cohorts of an immediate start, whose
#'   first interview can come at a later wave.
#' - `life_stage` is the stage of the pattern the cohort is at,
#'   out-of-sample stages counted. `life_length` is the number of periods
#'   the panel is held. A life is what Kalton (2009) calls a panel's
#'   duration.
#' - `equilibrium` is `TRUE` at a period holding one cohort at every
#'   interview of the life, after which response and precision stop
#'   changing (Lynn, 2012). Under an immediate start the full sample is in
#'   the field from period 1, and equilibrium comes once the start-up
#'   cohorts have aged through, the start-up period of Kalton (2009, p. 97).
#'
#' Under `"2-2-2"` a cohort selected at a later period has `life_stage` 1 to
#' 6, `wave` and `visit` 1, 2, `NA`, `NA`, 3, 4. The start-up cohort with two
#' periods left has `life_stage` 5 and 6 and `wave` 3 and 4, at `visit` 1
#' and 2.
#'
#' An immediate start covers every stage of the life at the first period,
#' one start-up cohort per stage, as [n_panel()] plans it. A cohort covering
#' a stage out of sample is selected at the first period and held until the
#' rest of its life returns to sample. Its first interview, where `visit` is
#' 1, comes later. Under `"2-2-2"` the cohorts covering stages 3 and 4 are
#' first interviewed at periods 3 and 2.
#'
#' ## PSUs, and retained PSUs
#'
#' With `take`, the schedule counts PSUs. Each panel is a group of
#' `ceiling(n_entrants / take)` PSUs, and a group takes `take` units in each
#' PSU for every cohort. The PSU allocation then reads its total from the
#' schedule, as in `n_alloc(frame, n = schedule$selections$operational_issue[1],
#' n_step = G)` with `G` the number of start-up panels, so that the master
#' has the PSU count the schedule promises.
#'
#' Under `retain_psus = TRUE` no PSU is replaced inside the horizon. An
#' immediate start selects every group at the first period, each group's
#' first cohort covering the rest of the life from the stage it starts at,
#' and a gradual start brings in one group a period until the life is
#' covered. When a group's cohort finishes, the group takes a new cohort at
#' wave 1, as the PSUs of a Statistics Canada rotation group take new
#' dwellings (Statistics Canada, 2017). Every panel is then its own rotation
#' group. A new cohort is a new draw of units, so the overlap of the units
#' is still the one [design_overlap()] gives for the rotation, and `visit`
#' restarts at 1 with each cohort. The horizon policies apply to cohorts.
#' `"continuing"` lists the rest of every cohort taken in the window,
#' `"truncate_lives"` cuts cohorts at the horizon, and `"close_intake"`
#' takes no cohort that cannot finish by it.
#'
#' ## Subrounds
#'
#' Monthly fieldwork for quarterly estimates divides each quarterly round
#' into monthly subrounds, as the ILO manual describes it (Hussmanns,
#' Mehran and Verma, 1990, p. 225). With `subrounds = s`, each reporting
#' panel `g` becomes `s` panels, numbered `(g - 1) * s + j` for `j` in 1 to
#' `s`, so the reporting panel of a panel is `ceiling(panel / s)`. Panel `j`
#' is interviewed in subround `j` of every period its reporting panel is
#' interviewed in, and is held in the other subrounds. A cohort is then
#' `s` panels, one per subround. A unit keeps its subround at every visit,
#' and a split panel keeps the rotation group and the cohorts of its
#' reporting panel. Every count in the schedule is in periods, and a cohort
#' entering in period `t` has all its subround panels selected then. Each
#' subround panel takes `ceiling(n_entrants / s)` units, or that over `take`
#' in PSUs, so a period issues up to `s - 1` more per group than the
#' continuous plan. The start-up and overlap checks run on the reporting
#' grid, and `$overlap` is the reporting-period profile.
#'
#' ## Keeping a schedule
#'
#' A schedule is rebuilt from its panel plan, not reread. Keep the panel
#' plan's call and the call to `design_schedule()`, and run them again with
#' the installed svyplan. A schedule saved with [saveRDS()] under another
#' version of svyplan may carry columns whose meaning has since changed.
#'
#' @references
#' Commission Implementing Regulation (EU) 2019/2240 of 16 December 2019
#' specifying the technical items of the data set, establishing the
#' technical formats for transmission of information and specifying the
#' detailed arrangements and content of the quality reports on the
#' organisation of a sample survey in the labour force domain in accordance
#' with Regulation (EU) 2019/1700 of the European Parliament and of the
#' Council. Article 2.
#'
#' Hussmanns, R., Mehran, F., and Verma, V. (1990). *Surveys of Economically
#' Active Population, Employment, Unemployment and Underemployment: An ILO
#' Manual on Concepts and Methods*. International Labour Office, Geneva.
#' Chapter 11.
#'
#' Kalton, G. (2009). Designs for surveys over time. In D. Pfeffermann and
#' C. R. Rao (eds.), *Handbook of Statistics 29A, Sample Surveys: Design,
#' Methods and Applications*, 89--108. Elsevier.
#'
#' Lynn, P. (2012). *Longitudinal Survey Methods for the Household Finances
#' and Consumption Survey*. Report prepared for the European Central Bank.
#'
#' Office for National Statistics (2026). *Labour Force Survey User Guide,
#' Volume 1: LFS Background and Methodology*. Section 3.
#'
#' Statistics Canada (2017). *Methodology of the Canadian Labour Force
#' Survey*. Catalogue no. 71-526-X.
#'
#' Steel, D. and McLaren, C. (2008). *Design and Analysis of Repeated
#' Surveys*. Working Paper 11-08, Centre for Statistical and Survey
#' Methodology, University of Wollongong.
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
#' schedule$selections
#'
#' # A 2-2-2 rotation: six start-up cohorts, four interviewed at each period
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
#' startup <- gapped$schedule[gapped$schedule$selection == "startup", ]
#' xtabs(visit ~ panel + period, startup, addNA = TRUE)
#' xtabs(wave ~ panel + period, startup, addNA = TRUE)
#'
#' # Quarterly reporting on a life of five quarters, 14 units per PSU, PSUs
#' # retained and monthly fieldwork: 15 panels, each in its own month
#' monthly <- design_schedule(
#'   n_panel(n_prop(p = 0.5, moe = 0.02, deff = 1.8),
#'           retention = rep(0.9, 4), resp_rate = 0.8,
#'           design = "rotating", start = "immediate"),
#'   horizon = 4,
#'   horizon_policy = "continuing",
#'   refreshment = "entrant_register",
#'   rounding = "ceiling",
#'   take = 14,
#'   retain_psus = TRUE,
#'   subrounds = 3
#' )
#' monthly
#' xtabs(visit ~ panel + subround, monthly$schedule[monthly$schedule$period == 2, ])
#'
#' @seealso [design_rotation()] for the pattern this schedules, [n_panel()]
#'   and [prec_panel()] for the rotating plan, [plot.svyplan_rotation()] for
#'   the rotation chart and its cell table, and [print.svyplan_schedule()]
#'   for the printed form and its `summary()`.
#'
#' @family repeated survey functions
#' @export
design_schedule <- function(
  panel_plan,
  horizon,
  horizon_policy,
  refreshment,
  frame_vintage = NULL,
  rounding = "ceiling",
  take = NULL,
  retain_psus = FALSE,
  subrounds = 1,
  entrants = NULL
) {
  if (!inherits(panel_plan, "svyplan_panel") ||
      !identical(panel_plan$design, "rotating")) {
    stop("'panel_plan' must be a rotating svyplan_panel", call. = FALSE)
  }
  if (is.null(panel_plan$start) || is.null(panel_plan$periods) ||
      is.null(panel_plan$period_visits)) {
    stop(
      "'panel_plan' must carry a start; call n_panel() or prec_panel() with start = \"immediate\" or \"gradual\"",
      call. = FALSE
    )
  }
  entrants <- .schedule_entrants(entrants, panel_plan)
  if (identical(entrants, "assured")) {
    # Every reconciliation below, the start-up visits, the standing sample
    # and the held units, then follows from the count fielded.
    panel_plan <- .panel_at(panel_plan, panel_plan$n_assured)
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
          "under \"truncate_lives\" the horizon has to reach the first interview of every startup panel, which the last of them has at period %d",
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
  rounding <- .schedule_choice(rounding, "ceiling", "rounding")
  take <- .schedule_take(take)
  retain_psus <- .schedule_retain(retain_psus, take)
  subrounds <- .schedule_whole(subrounds, "subrounds")

  entrant_issue <- unname(panel_plan$n_entrants)
  plan <- .schedule_plan(
    panel_plan$start, life_length, horizon, horizon_policy, retain_psus
  )
  per_panel <- entrant_issue / (.schedule_unit(take) * subrounds)
  # Whole-valued doubles remain exact far beyond R's 32-bit integer range.
  # Operational counts are counts by value, not storage-mode integers.
  panel_issue <- ceiling(per_panel)
  panels <- plan$panels * subrounds
  selections <- data.frame(
    selection = plan$selection,
    entry_period = as.integer(plan$entry_period),
    frame_role = c("startup", rep(refreshment, nrow(plan) - 1L)),
    planned_issue = per_panel * panels,
    operational_issue = panel_issue * panels,
    panels = as.integer(panels),
    panel_issue = rep(panel_issue, nrow(plan)),
    frame_vintage = .schedule_frame_vintages(frame_vintage, plan$selection),
    status = rep("planned", nrow(plan)),
    stringsAsFactors = FALSE
  )

  tables <- .schedule_tables(
    plan, panel_plan$start, w, horizon, horizon_policy, retain_psus,
    subrounds
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
  .check_schedule_startup(
    reporting, panel_plan$period_visits, horizon_policy, last_open
  )
  .check_schedule_overlap(reporting, as.double(overlap))
  .check_schedule_population(
    tables$schedule, panel_issue, take,
    panel_plan$target$params$N %||% Inf, panel_plan$n_held
  )

  schedule <- tables$schedule
  settled <- which(tables$equilibrium)
  equilibrium_from <- if (length(settled) == 0L) {
    NA_integer_
  } else {
    as.integer(settled[[1L]])
  }

  .new_svyplan_schedule(
    life_spec = as.numeric(w),
    life_notation = "resolved_per_period_take",
    life_length = life_length,
    horizon = horizon,
    horizon_policy = horizon_policy,
    start = panel_plan$start,
    target_estimand = list(type = panel_plan$type, method = panel_plan$method),
    n_target = unname(panel_plan$n_target),
    response_model = list(
      resp_rate = panel_plan$params$resp_rate,
      eligibility_rate = panel_plan$params[["eligibility_rate"]] %||% 1,
      retention = panel_plan$params$retention,
      q = panel_plan$waves$q
    ),
    rounding = list(
      rule = rounding, level = .schedule_rounding_level(take, subrounds)
    ),
    take = take,
    retain_psus = retain_psus,
    subrounds = subrounds,
    equilibrium_from = equilibrium_from,
    overlap = as.data.frame(overlap),
    refreshment = refreshment,
    requires_disjoint_entrants = identical(refreshment, "entrant_register"),
    n_entrants = entrant_issue,
    entrant_basis = entrants,
    n_in_sample = unname(panel_plan$n_in_sample),
    n_held = unname(panel_plan$n_held),
    n_cohorts = as.integer(panel_plan$n_cohorts),
    panel_parameters = .schedule_panel_parameters(
      schedule, tables$tail, selections$panels[[1L]]
    ),
    selections = selections,
    schedule = schedule,
    tail_commitments = tables$tail,
    issue = .schedule_issue_profile(
      selections, schedule, tables$equilibrium, subrounds,
      .schedule_unit(take)
    )
  )
}

#' Which entrant count the schedule fields
#'
#' A plan with an assurance carries two counts, and fielding the expected
#' one silently gives up the assurance, so the choice is required there.
#' @keywords internal
#' @noRd
.schedule_entrants <- function(entrants, panel_plan) {
  assured <- panel_plan$n_assured
  if (is.null(entrants)) {
    if (!is.null(assured)) {
      stop(
        sprintf(
          "the plan carries an assured entrant count, %s per period against %s expected; declare entrants = \"assured\" to field it, or entrants = \"expected\" to field the expected count and give up the assurance",
          .fmt_count_n(assured), format(signif(panel_plan$n_entrants, 4))
        ),
        call. = FALSE
      )
    }
    return("expected")
  }
  entrants <- .schedule_choice(entrants, c("expected", "assured"), "entrants")
  if (identical(entrants, "assured") && is.null(assured)) {
    stop(
      "entrants = \"assured\" needs a plan with an assurance; plan it with n_panel(..., assurance = )",
      call. = FALSE
    )
  }
  entrants
}

#' Refuse a rounded schedule that holds more units than the population
#'
#' The plan passed its own frame check on continuous counts. Rounding each
#' panel up, by subround and by PSU, can push the units held at once past
#' `N`. Held counts every selected panel whose life covers the period, its
#' out-of-sample spells included, since its units stay reserved.
#' @keywords internal
#' @noRd
.check_schedule_population <- function(schedule, panel_issue, take, N,
                                       planned_held) {
  if (!is.finite(N)) {
    return(invisible(TRUE))
  }
  units_per_panel <- panel_issue * .schedule_unit(take)
  size <- if (is.na(take)) {
    sprintf("%s units", .fmt_count_n(units_per_panel))
  } else {
    sprintf("%s PSUs of %s units", .fmt_count_n(panel_issue),
            .fmt_count_n(take))
  }
  held <- schedule[!is.na(schedule$life_stage), c("period", "panel")]
  panels <- tapply(held$panel, held$period, function(z) length(unique(z)))
  units <- panels * units_per_panel
  worst <- which.max(units)
  if (units[[worst]] > N) {
    stop(
      sprintf(
        "rounding holds %s frame units at period %s (%s panels of %s), more than the population of %s; the plan holds %s before rounding. Use a smaller 'take' or fewer 'subrounds', or plan from a smaller target",
        .fmt_count_n(units[[worst]]), names(units)[[worst]],
        .fmt_count_n(panels[[worst]]), size,
        .fmt_count_n(N), format(signif(planned_held, 4))
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
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
      "'take' must be NULL or a positive whole number of units per PSU per cohort",
      call. = FALSE
    )
  }
  as.double(take)
}

#' @keywords internal
#' @noRd
.schedule_retain <- function(retain_psus, take) {
  if (!is.logical(retain_psus) || length(retain_psus) != 1L ||
      is.na(retain_psus)) {
    stop("'retain_psus' must be TRUE or FALSE", call. = FALSE)
  }
  if (retain_psus && is.na(take)) {
    stop(
      "'retain_psus' needs 'take': the schedule counts PSUs only when it knows how many units each takes per cohort",
      call. = FALSE
    )
  }
  retain_psus
}

#' Units an assignment unit stands for, 1 without a take
#' @keywords internal
#' @noRd
.schedule_unit <- function(take) {
  if (is.na(take)) 1 else take
}

#' @keywords internal
#' @noRd
.schedule_rounding_level <- function(take, subrounds) {
  parts <- c(if (!is.na(take)) "psu", if (subrounds > 1L) "subround")
  if (length(parts) == 0L) "panel" else paste(parts, collapse = "_and_")
}

#' The selections of a schedule on the reporting grid
#' @keywords internal
#' @noRd
.schedule_plan <- function(start, life_length, horizon, policy, retain) {
  max_entry <- if (identical(policy, "close_intake")) {
    horizon - life_length + 1L
  } else {
    horizon
  }
  if (retain) {
    max_entry <- if (identical(start, "immediate")) {
      1L
    } else {
      min(life_length, max_entry)
    }
  }
  max_entry <- as.integer(max_entry)
  data.frame(
    selection = c("startup", if (max_entry >= 2L) paste0("period_", 2:max_entry)),
    entry_period = seq_len(max_entry),
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
.schedule_frame_vintages <- function(x, selections) {
  out <- rep(NA_character_, length(selections))
  if (is.null(x)) {
    return(out)
  }
  ok <- is.character(x) && !anyNA(x) && length(x) > 0L &&
    !is.null(names(x)) && all(nzchar(names(x))) && anyDuplicated(names(x)) == 0L
  if (!ok || !all(names(x) %in% selections)) {
    stop(
      "'frame_vintage' must be NULL or a named character vector whose names are generated selections",
      call. = FALSE
    )
  }
  out[match(names(x), selections)] <- unname(x)
  out
}

#' The panels a schedule fields, in the chart's order
#'
#' Panels are numbered across the whole schedule as the cell table numbers
#' them: the start-up panels from the one-period tail to the full life, as
#' [plot.svyplan_rotation()] and Lynn's Figure 5 number them, then one per
#' later selection. A panel's rotation group is read from its first cohort's
#' whole life, untouched by the horizon policy, so a schedule and the cell
#' table of the same rotation assign the same groups. With PSUs retained the
#' group is read from every cohort the panel takes. The cohort offset makes
#' the first start-up cohort number 1.
#' @keywords internal
#' @noRd
.schedule_panel_catalog <- function(plan, start, w, horizon,
                                    horizon_policy, retain) {
  life <- length(w)
  startup <- .rotation_cohorts(w, start, 1L)
  startup_first <- life - lengths(lapply(startup, `[[`, "w")) + 1L
  rows <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    k <- plan$panels[[i]]
    first <- if (i == 1L) startup_first else rep(1L, k)
    h <- horizon - plan$entry_period[[i]] + 1L
    len <- vapply(first, function(f) {
      .schedule_life_length(w, f, h, horizon_policy, retain)
    }, integer(1L))
    rows[[i]] <- data.frame(
      selection = rep(plan$selection[[i]], k),
      entry_period = rep(plan$entry_period[[i]], k),
      first_stage = as.integer(first),
      life_length = len,
      relative_take = rep(1, k),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$panel <- seq_len(nrow(out))
  out$rotation_group <- .rotation_groups(lapply(seq_len(nrow(out)), function(i) {
    first <- out$first_stage[[i]]
    stage <- seq_len(
      if (!retain) life - first + 1L else out$life_length[[i]]
    )
    at <- .schedule_position(w, first, stage)
    (out$entry_period[[i]] + stage - 1L)[w[at$pos] > 0]
  }))
  attr(out, "cohort_offset") <- if (identical(start, "immediate")) {
    life - 1L
  } else {
    0L
  }
  out
}

#' Where a panel's own stage sits in the life, each new cohort repeating it
#' @keywords internal
#' @noRd
.schedule_position <- function(w, first, stage) {
  life <- length(w)
  first_cohort <- life - first + 1L
  later <- stage > first_cohort
  in_cohort <- ifelse(later, (stage - first_cohort - 1L) %% life + 1L, stage)
  list(
    pos = ifelse(later, in_cohort, first + stage - 1L),
    start = ifelse(later, 1L, first),
    in_cohort = in_cohort
  )
}

#' Periods of a panel's own life inside and beyond the window
#' @keywords internal
#' @noRd
.schedule_life_length <- function(w, first, h, policy, retain) {
  life <- length(w)
  first_cohort <- life - first + 1L
  if (identical(policy, "truncate_lives")) {
    cut <- if (!retain) min(first_cohort, h) else h
    on <- w[.schedule_position(w, first, seq_len(cut))$pos] > 0
    return(as.integer(max(which(on))))
  }
  if (!retain || h <= first_cohort) {
    return(as.integer(first_cohort))
  }
  more <- if (identical(policy, "close_intake")) {
    (h - first_cohort) %/% life
  } else {
    ceiling((h - first_cohort) / life)
  }
  as.integer(first_cohort + more * life)
}

#' Column order of the dense schedule
#' @keywords internal
#' @noRd
.schedule_columns <- function() {
  c("period", "subround", "selection", "cohort", "panel", "rotation_group",
    "wave", "visit", "life_stage", "life_length", "active", "equilibrium",
    "relative_take")
}

#' Every catalog panel at one period: held, interviewed, and its wave
#'
#' `wave` is the place of the interview in the rotation pattern, `visit`
#' the interview count of the row's cohort, and `life_stage` the stage of
#' the pattern the cohort is at. The cohort is numbered from the period its
#' life starts, counted back before period 1 for a start-up cohort.
#' `.cstage` is the stage within the current cohort's life, which the
#' equilibrium test counts and the stored tables drop.
#' @keywords internal
#' @noRd
.schedule_rows <- function(catalog, period, w) {
  seen <- c(0L, cumsum(w > 0))
  stage <- period - catalog$entry_period + 1L
  held <- stage >= 1L & stage <= catalog$life_length
  at <- .schedule_position(w, catalog$first_stage, ifelse(held, stage, 1L))
  active <- held & w[at$pos] > 0
  offset <- attr(catalog, "cohort_offset", exact = TRUE)
  data.frame(
    period = rep(as.integer(period), nrow(catalog)),
    subround = rep(1L, nrow(catalog)),
    selection = catalog$selection,
    cohort = ifelse(
      held, as.integer(period - at$pos + 1L + offset), NA_integer_
    ),
    panel = as.integer(catalog$panel),
    rotation_group = as.integer(catalog$rotation_group),
    wave = ifelse(active, as.integer(seen[at$pos + 1L]), NA_integer_),
    visit = ifelse(
      active, as.integer(seen[at$pos + 1L] - seen[at$start]), NA_integer_
    ),
    life_stage = ifelse(held, as.integer(at$pos), NA_integer_),
    life_length = as.integer(catalog$life_length),
    active = active,
    relative_take = catalog$relative_take,
    .cstage = ifelse(held, as.integer(at$in_cohort), NA_integer_),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.schedule_dense_grid <- function(catalog, horizon, w) {
  out <- do.call(rbind, lapply(seq_len(horizon), function(period) {
    .schedule_rows(catalog, period, w)
  }))
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.empty_schedule_tail <- function() {
  data.frame(
    period = integer(),
    subround = integer(),
    selection = character(),
    cohort = integer(),
    panel = integer(),
    rotation_group = integer(),
    wave = integer(),
    visit = integer(),
    life_stage = integer(),
    life_length = integer(),
    relative_take = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Interviews owed after the horizon, and only interviews
#' @keywords internal
#' @noRd
.schedule_tail <- function(catalog, horizon, policy, w) {
  last <- max(catalog$entry_period + catalog$life_length - 1L)
  if (!identical(policy, "continuing") || last <= horizon) {
    return(.empty_schedule_tail())
  }
  out <- do.call(rbind, lapply(seq.int(horizon + 1L, last), function(period) {
    .schedule_rows(catalog, period, w)
  }))
  out <- out[out$active, names(.empty_schedule_tail()), drop = FALSE]
  if (nrow(out) == 0L) {
    return(.empty_schedule_tail())
  }
  rownames(out) <- NULL
  out
}

#' Periods where every stage of a life is held exactly once
#'
#' That is equilibrium in the sense of Lynn (2012): one panel at each
#' interview of the life, so response and precision stop changing.
#' @keywords internal
#' @noRd
.schedule_equilibrium <- function(schedule, life_length) {
  vapply(seq_len(max(schedule$period)), function(period) {
    at <- schedule[
      schedule$period == period & !is.na(schedule$.cstage),
      c(".cstage", "relative_take"),
      drop = FALSE
    ]
    totals <- numeric(life_length)
    if (nrow(at) > 0L) {
      by_stage <- tapply(at$relative_take, at$.cstage, sum)
      totals[as.integer(names(by_stage))] <- as.numeric(by_stage)
    }
    isTRUE(all.equal(totals, rep(1, life_length), tolerance = 1e-12))
  }, logical(1L))
}

#' The reporting grid, and the fieldwork tables split from it
#' @keywords internal
#' @noRd
.schedule_tables <- function(plan, start, w, horizon, policy, retain,
                             subrounds) {
  catalog <- .schedule_panel_catalog(plan, start, w, horizon, policy, retain)
  grid <- .schedule_dense_grid(catalog, horizon, w)
  tail <- .schedule_tail(catalog, horizon, policy, w)
  equilibrium <- .schedule_equilibrium(grid, length(w))
  grid$equilibrium <- equilibrium[grid$period]
  grid <- grid[.schedule_columns()]
  if (subrounds == 1L) {
    return(list(reporting = grid, schedule = grid, tail = tail,
                equilibrium = equilibrium))
  }
  list(
    reporting = grid,
    schedule = .schedule_split(grid, subrounds, dense = TRUE),
    tail = .schedule_split(tail, subrounds, dense = FALSE),
    equilibrium = equilibrium
  )
}

#' Split each reporting panel into its subround panels
#'
#' Panel `(g - 1) s + j` is interviewed in subround `j` of every period its
#' group is interviewed in, and held in the other subrounds.
#' @keywords internal
#' @noRd
.schedule_split <- function(tab, s, dense) {
  if (nrow(tab) == 0L) {
    return(tab)
  }
  combos <- if (dense) {
    expand.grid(j = seq_len(s), m = seq_len(s))
  } else {
    data.frame(j = seq_len(s), m = seq_len(s))
  }
  row <- rep(seq_len(nrow(tab)), each = nrow(combos))
  j <- rep(combos$j, times = nrow(tab))
  m <- rep(combos$m, times = nrow(tab))
  out <- tab[row, , drop = FALSE]
  out$subround <- as.integer(m)
  out$panel <- as.integer((tab$panel[row] - 1L) * s + j)
  if (dense) {
    out$active <- tab$active[row] & m == j
    out$visit <- ifelse(out$active, tab$visit[row], NA_integer_)
    out$wave <- ifelse(out$active, tab$wave[row], NA_integer_)
  }
  out <- out[order(out$period, out$subround, out$panel), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Issue at each subround of each period
#' @keywords internal
#' @noRd
.schedule_issue_profile <- function(selections, schedule, equilibrium, s,
                                    unit) {
  horizon <- length(equilibrium)
  opening <- schedule$active & !is.na(schedule$visit) & schedule$visit == 1L
  key <- (schedule$period - 1L) * s + schedule$subround
  opened <- tabulate(key[opening], horizon * s)
  rows <- lapply(seq_len(horizon * s), function(r) {
    period <- (r - 1L) %/% s + 1L
    subround <- (r - 1L) %% s + 1L
    at <- selections[
      selections$entry_period == period & subround == 1L, , drop = FALSE
    ]
    data.frame(
      period = period,
      subround = subround,
      selection = if (nrow(at) == 0L) NA_character_ else at$selection[[1L]],
      planned_issue = if (nrow(at) == 0L) 0 else at$planned_issue[[1L]],
      operational_issue = if (nrow(at) == 0L) 0 else at$operational_issue[[1L]],
      cohort_issue = opened[[r]] * selections$panel_issue[[1L]] * unit,
      equilibrium = equilibrium[[period]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$period <- as.integer(out$period)
  out$subround <- as.integer(out$subround)
  out
}

#' @keywords internal
#' @noRd
.schedule_panel_parameters <- function(schedule, tail, k) {
  if (k < 2L) {
    return(list(k = as.integer(k), r_min = 1L, block_width = NA_integer_))
  }
  startup <- schedule[schedule$selection == "startup" & schedule$active, ]
  later <- tail[tail$selection == "startup", ]
  per_period <- table(c(paste(startup$period, startup$subround),
                          paste(later$period, later$subround)))
  r_min <- as.integer(min(per_period))
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

#' Validate the package-neutral schedule layout
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
    "life_spec", "life_notation", "life_length",
    "horizon", "horizon_policy", "start", "target_estimand",
    "n_target", "response_model", "rounding", "take", "retain_psus",
    "subrounds", "equilibrium_from",
    "overlap", "refreshment", "requires_disjoint_entrants", "n_entrants",
    "entrant_basis", "n_in_sample", "n_held", "n_cohorts", "panel_parameters", "selections",
    "schedule", "tail_commitments", "issue"
  )
  if (!is.list(x) || !identical(names(x), required)) {
    fail("top-level fields are not those of a schedule")
  }
  if (!is.integer(x$life_length) || length(x$life_length) != 1L ||
      is.na(x$life_length) || x$life_length < 2L ||
      !is.integer(x$horizon) || length(x$horizon) != 1L ||
      is.na(x$horizon) || x$horizon < 1L ||
      !identical(x$life_notation, "resolved_per_period_take") ||
      !is.character(x$horizon_policy) || length(x$horizon_policy) != 1L ||
      is.na(x$horizon_policy) ||
      !x$horizon_policy %in% c("continuing", "truncate_lives", "close_intake") ||
      !is.character(x$start) || length(x$start) != 1L ||
      is.na(x$start) ||
      !x$start %in% c("immediate", "gradual")) {
    fail("scalar metadata is inconsistent")
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
      !identical(names(response),
                 c("resp_rate", "eligibility_rate", "retention", "q")) ||
      !scalar_number(response$resp_rate, lower = 0, upper = 1) ||
      response$resp_rate == 0 ||
      !scalar_number(response$eligibility_rate, lower = 0, upper = 1) ||
      response$eligibility_rate == 0 ||
      !is.numeric(response$retention) ||
      length(response$retention) != visits - 1L ||
      anyNA(response$retention) || any(!is.finite(response$retention)) ||
      any(response$retention <= 0 | response$retention > 1) ||
      !is.numeric(response$q) || length(response$q) != visits ||
      anyNA(response$q) || any(!is.finite(response$q))) {
    fail("response metadata is inconsistent")
  }
  expected_q <- response$resp_rate * response$eligibility_rate *
    cumprod(c(1, response$retention))
  if (!same(as.numeric(response$q), expected_q)) {
    fail("response probabilities do not follow resp_rate, eligibility_rate and retention")
  }
  if (!is.double(x$take) || length(x$take) != 1L ||
      !(is.na(x$take) || (is.finite(x$take) && x$take >= 1 &&
                            x$take == trunc(x$take))) ||
      !isTRUE(x$retain_psus) && !isFALSE(x$retain_psus) ||
      x$retain_psus && is.na(x$take) ||
      !is.integer(x$subrounds) || length(x$subrounds) != 1L ||
      is.na(x$subrounds) || x$subrounds < 1L) {
    fail("take, retain_psus or subrounds metadata is inconsistent")
  }
  s <- x$subrounds
  horizon <- x$horizon
  unit <- .schedule_unit(x$take)
  if (!is.list(x$rounding) ||
      !identical(x$rounding, list(
        rule = "ceiling", level = .schedule_rounding_level(x$take, s)
      ))) {
    fail("rounding metadata is inconsistent")
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
  if (!is.character(x$entrant_basis) || length(x$entrant_basis) != 1L ||
      !x$entrant_basis %in% c("expected", "assured")) {
    fail("entrant_basis must be \"expected\" or \"assured\"")
  }
  if (!scalar_number(x$n_entrants, lower = 0) || x$n_entrants == 0 ||
      !scalar_number(x$n_in_sample, lower = 0) || x$n_in_sample == 0 ||
      !scalar_number(x$n_held, lower = 0) ||
      !identical(x$n_cohorts, visits) ||
      !same(x$n_in_sample, x$n_cohorts * x$n_entrants) ||
      !same(x$n_held, x$life_length * x$n_entrants)) {
    fail("panel counts are inconsistent")
  }

  selection_cols <- c(
    "selection", "entry_period", "frame_role", "planned_issue",
    "operational_issue", "panels", "panel_issue", "frame_vintage", "status"
  )
  schedule_cols <- .schedule_columns()
  tail_cols <- names(.empty_schedule_tail())
  issue_cols <- c(
    "period", "subround", "selection", "planned_issue", "operational_issue",
    "cohort_issue", "equilibrium"
  )
  exact_table <- function(tab, cols) {
    is.data.frame(tab) && identical(names(tab), cols)
  }
  if (!exact_table(x$selections, selection_cols) ||
      !exact_table(x$schedule, schedule_cols) ||
      !exact_table(x$tail_commitments, tail_cols) ||
      !exact_table(x$issue, issue_cols)) {
    fail("one or more tables do not have the schedule's columns")
  }

  selections <- x$selections
  plan <- .schedule_plan(
    x$start, x$life_length, horizon, x$horizon_policy, x$retain_psus
  )
  expected_panels <- as.integer(plan$panels * s)
  expected_roles <- c("startup", rep(x$refreshment, nrow(plan) - 1L))
  if (!identical(selections$selection, plan$selection) ||
      anyNA(selections$selection) || anyDuplicated(selections$selection) ||
      !is.integer(selections$entry_period) ||
      anyNA(selections$entry_period) ||
      !identical(selections$entry_period, as.integer(plan$entry_period)) ||
      !is.numeric(selections$planned_issue) ||
      any(!is.finite(selections$planned_issue) | selections$planned_issue <= 0) ||
      !whole_counts(selections$operational_issue, lower = 1) ||
      any(selections$operational_issue < selections$planned_issue) ||
      !identical(selections$panels, expected_panels) ||
      !whole_counts(selections$panel_issue, lower = 1) ||
      !is.character(selections$frame_vintage) ||
      length(selections$frame_vintage) != nrow(selections) ||
      any(!is.na(selections$frame_vintage) &
            !nzchar(selections$frame_vintage)) ||
      !identical(selections$frame_role, expected_roles) ||
      !identical(selections$status, rep("planned", nrow(selections)))) {
    fail("the selection manifest is inconsistent")
  }
  per_panel <- x$n_entrants / (unit * s)
  expected_panel_issue <- rep(ceiling(per_panel), nrow(selections))
  if (!same(selections$planned_issue, per_panel * expected_panels) ||
      !same(selections$panel_issue, expected_panel_issue) ||
      !same(
        selections$operational_issue,
        expected_panel_issue * expected_panels
      )) {
    fail("selection issue does not follow panel-level ceiling")
  }

  schedule <- x$schedule
  n_panels <- sum(selections$panels)
  bounded <- !x$retain_psus
  int_cols <- c("period", "subround", "cohort", "panel", "rotation_group",
                "wave", "visit", "life_stage", "life_length")
  if (!all(vapply(schedule[int_cols], is.integer, logical(1L))) ||
      !is.logical(schedule$active) || !is.logical(schedule$equilibrium) ||
      !is.numeric(schedule$relative_take) ||
      nrow(schedule) != horizon * s * n_panels ||
      anyNA(schedule[c("period", "subround", "selection", "panel",
                       "rotation_group", "life_length", "active",
                       "equilibrium", "relative_take")]) ||
      any(schedule$period < 1L | schedule$period > horizon) ||
      any(schedule$subround < 1L | schedule$subround > s) ||
      any(!schedule$selection %in% selections$selection) ||
      any(schedule$panel < 1L) ||
      any(is.na(schedule$cohort) != is.na(schedule$life_stage)) ||
      any(!is.na(schedule$cohort) & schedule$cohort < 1L) ||
      any(schedule$rotation_group < 1L) ||
      any(schedule$life_length < 1L) ||
      bounded && any(schedule$life_length > x$life_length) ||
      any(schedule$relative_take != 1) ||
      any(is.na(schedule$visit) != !schedule$active) ||
      any(is.na(schedule$wave) != !schedule$active) ||
      any(is.na(schedule$life_stage) & schedule$active) ||
      any(!is.na(schedule$life_stage) &
            (schedule$life_stage < 1L |
               schedule$life_stage > x$life_length)) ||
      any(!is.na(schedule$wave) &
            (schedule$wave < 1L | schedule$wave > visits))) {
    fail("the dense schedule grid is inconsistent")
  }
  grid_key <- paste(schedule$period, schedule$subround, schedule$panel,
                    sep = "|")
  if (anyDuplicated(grid_key)) {
    fail("the dense schedule grid contains duplicate panel rows")
  }
  if (any(tapply(schedule$selection, schedule$panel,
                 function(z) length(unique(z))) != 1L)) {
    fail("a panel belongs to more than one selection")
  }
  period_cols <- c("panel", "cohort", "period", "subround", "visit")
  interviews <- rbind(
    schedule[schedule$active, period_cols],
    x$tail_commitments[period_cols]
  )
  interviews <- interviews[order(interviews$panel, interviews$cohort,
                                 interviews$period, interviews$subround), ]
  if (anyNA(interviews$visit) || anyNA(interviews$cohort)) {
    fail("an interview has no visit or no cohort")
  }
  counted <- stats::ave(interviews$visit, interviews$panel, interviews$cohort,
                        FUN = seq_along)
  if (!identical(as.integer(counted), interviews$visit)) {
    fail("a visit does not count the interviews of its cohort")
  }

  tables <- .schedule_tables(
    plan, x$start, x$life_spec, horizon, x$horizon_policy,
    x$retain_psus, s
  )
  if (!identical(schedule, tables$schedule)) {
    fail("the dense schedule does not match the selection manifest")
  }

  tail <- x$tail_commitments
  if (!all(vapply(tail[int_cols], is.integer, logical(1L))) ||
      !is.numeric(tail$relative_take) || anyNA(tail) ||
      any(tail$period <= horizon) ||
      any(tail$subround < 1L | tail$subround > s) ||
      any(!tail$selection %in% selections$selection) ||
      any(tail$panel < 1L) ||
      any(tail$life_stage < 1L | tail$life_stage > x$life_length) ||
      any(tail$wave < 1L | tail$wave > visits) ||
      any(tail$life_length < 1L) ||
      bounded && any(tail$life_length > x$life_length) ||
      any(tail$relative_take != 1) ||
      (!identical(x$horizon_policy, "continuing") && nrow(tail) > 0L)) {
    fail("the tail-commitment table is inconsistent with the horizon policy")
  }
  if (!identical(tail, tables$tail)) {
    fail("tail commitments do not match the selection manifest")
  }

  issue <- x$issue
  if (!is.integer(issue$period) || !is.integer(issue$subround) ||
      !identical(issue$period, rep(seq_len(horizon), each = s)) ||
      !identical(issue$subround, rep(seq_len(s), times = horizon)) ||
      !is.numeric(issue$planned_issue) ||
      !whole_counts(issue$operational_issue) ||
      !whole_counts(issue$cohort_issue) ||
      !is.logical(issue$equilibrium) || anyNA(issue$equilibrium) ||
      any(issue$planned_issue < 0) ||
      any(issue$operational_issue < issue$planned_issue)) {
    fail("the issue profile is inconsistent")
  }
  expected_issue <- .schedule_issue_profile(
    selections, tables$schedule, tables$equilibrium, s, unit
  )
  if (!identical(issue, expected_issue)) {
    fail("the issue profile does not match the selection manifest")
  }

  settled <- which(tables$equilibrium)
  expected_from <- if (length(settled) == 0L) {
    NA_integer_
  } else {
    as.integer(settled[[1L]])
  }
  if (!is.integer(x$equilibrium_from) ||
      length(x$equilibrium_from) != 1L ||
      !identical(x$equilibrium_from, expected_from)) {
    fail("equilibrium metadata does not match the schedule")
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
    period = tab$period,
    visit = tab$visit,
    n_issued = x$n_entrants * tab$relative_take,
    n_resp = x$n_entrants * tab$relative_take * q[tab$visit]
  )
  out <- stats::aggregate(
    cbind(n_issued, n_resp) ~ period + visit,
    data = out,
    FUN = sum
  )
  out <- out[order(out$period, out$visit), ]
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.check_schedule_startup <- function(x, period_visits, policy, max_entry) {
  top <- min(x$horizon, x$life_length + 1L)
  if (identical(policy, "close_intake")) {
    top <- min(top, max_entry)
  }
  got <- .schedule_active_aggregate(x)
  got <- got[got$period <= top, , drop = FALSE]
  want <- period_visits[period_visits$period <= top, , drop = FALSE]
  rownames(want) <- NULL
  if (!isTRUE(all.equal(got, want, tolerance = 1e-10, check.attributes = FALSE))) {
    stop(
      "internal schedule defect: generated activity does not reconcile with panel_plan$period_visits",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.check_schedule_overlap <- function(x, expected) {
  active <- x$schedule[x$schedule$active, , drop = FALSE]
  key <- paste(active$panel, active$cohort, sep = "|")
  take <- vapply(seq_len(x$horizon), function(p) {
    sum(active$relative_take[active$period == p])
  }, numeric(1L))
  for (lag in seq_len(min(length(expected), x$horizon - 1L))) {
    origins <- which(
      seq_len(x$horizon) + lag <= x$horizon &
        abs(take - sum(x$life_spec > 0)) <= 1e-10
    )
    for (origin in origins) {
      first <- unique(key[active$period == origin])
      second <- unique(key[active$period == origin + lag])
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
#' row per period, because a schedule's periods are mostly identical and
#' its length is set by the reporting horizon rather than by the design.
#' `summary()` gives the period-by-period tables: the issue profile in
#' full, the selections, the overlap the rotation produces, and the
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
#' # the period-by-period tables the print block summarizes
#' s <- summary(schedule)
#' s
#'
#' # every table stays reachable as a field
#' schedule$selections
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
    "Longitudinal design schedule (%s start, %s)\n",
    x$start, x$horizon_policy
  ))
  cat(.fmt_schedule_life(x))
  cat(.fmt_schedule_psu(x))
  cat(sprintf(
    "rounding: %s at %s level\n",
    x$rounding$rule, gsub("_", " ", x$rounding$level, fixed = TRUE)
  ))
  cat(sprintf("issue: %s\n", .fmt_schedule_issue(x)))
  cat(.fmt_schedule_tail(x))
  cat("# summary() for the period-by-period profile and the overlap\n")
  invisible(x)
}

#' The life, the window, and the period the interview mix settles
#' @keywords internal
#' @noRd
.fmt_schedule_life <- function(x) {
  settled <- if (is.na(x$equilibrium_from)) {
    "equilibrium not reached in the window"
  } else {
    sprintf("equilibrium from period %d", x$equilibrium_from)
  }
  paste0(
    sprintf("life: %d periods, %s in sample each period\n", x$life_length,
            format(sum(x$life_spec))),
    .fmt_schedule_rotation(x),
    sprintf("horizon: %d periods, %s\n", x$horizon, settled),
    if (x$subrounds > 1L) sprintf("subrounds: %d per period\n", x$subrounds)
  )
}

#' The PSU take and the PSU stay, where either was declared
#' @keywords internal
#' @noRd
.fmt_schedule_psu <- function(x) {
  paste0(
    if (!is.na(x$take)) {
      sprintf("take: %s units per PSU per cohort\n", format(x$take))
    },
    if (x$retain_psus) {
      "PSUs: retained through the horizon, a new cohort at each return\n"
    }
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
#' a row per period grows without bound while saying the same thing: a
#' start-up take, then one selection repeated. Naming the runs states the profile
#' at a size that does not move when the horizon does. A run of zeros is
#' intake having closed, which `horizon_policy = "close_intake"` produces and
#' which is a fact about the schedule rather than a gap in it.
#' @keywords internal
#' @noRd
.fmt_schedule_issue <- function(x) {
  v <- as.vector(tapply(x$issue$operational_issue, x$issue$period, sum))
  unit <- "period"
  psu <- if (is.na(x$take)) "" else " PSUs"
  closed <- if (x$retain_psus) "no new PSUs from" else {
    "intake closed from"
  }
  parts <- sprintf("%s%s at start-up", format(v[[1L]]), psu)
  if (length(v) > 1L) {
    runs <- rle(v[-1L])
    at <- 2L
    for (j in seq_along(runs$lengths)) {
      lo <- at
      hi <- at + runs$lengths[[j]] - 1L
      at <- hi + 1L
      parts <- c(parts, if (runs$values[[j]] == 0) {
        sprintf("%s %s %d", closed, unit, lo)
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
    "tail commitments: %d panel-interviews after period %d\n",
    nrow(x$tail_commitments), x$horizon
  )
}

#' Drop the selection columns that carry no distinction
#'
#' `frame_vintage` is absent unless one was declared and `status` is one
#' value for a plan that has not been executed against, so on an ordinary
#' schedule both are a column of the same entry repeated.
#' @keywords internal
#' @noRd
.fmt_schedule_selections <- function(d) {
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
      selections = object$selections,
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
    "Analysis of a longitudinal design schedule (%s start, %s)\n\n",
    s$start, s$horizon_policy
  ))
  cat(.fmt_schedule_life(s))
  cat(.fmt_schedule_psu(s))
  # From the rounded panels, which ceiling(n_entrants) understates.
  entrants <- s$selections$panel_issue[[1L]] * .schedule_unit(s$take) *
    s$subrounds
  cat(sprintf(
    "rounding: %s at %s level, %s entrants per period against the plan's %s (+%.1f%%)\n",
    s$rounding$rule, gsub("_", " ", s$rounding$level, fixed = TRUE),
    format(entrants), format(signif(s$n_entrants, 4)),
    100 * (entrants / s$n_entrants - 1)
  ))
  cat(sprintf(
    "entrants: %s per period%s, %s in sample across %d cohorts\n",
    format(entrants),
    if (identical(s$entrant_basis, "assured")) " (assured count)" else "",
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

  cat("\nSelections\n")
  print(.fmt_schedule_selections(x$selections), row.names = FALSE,
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
    "svyplan_schedule [%s start, %s, %d-period life, horizon %d periods]",
    x$start, x$horizon_policy, x$life_length, x$horizon
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
