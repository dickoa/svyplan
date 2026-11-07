schedule_panel_plan <- function(start = "immediate", life = 4L) {
  n_panel(
    n_prop(p = 0.5, moe = 0.03),
    retention = rep(0.9, life - 1L),
    resp_rate = 0.75,
    design = "rotating",
    start = start
  )
}

schedule_active_counts <- function(x) {
  tab <- x$schedule[x$schedule$active, , drop = FALSE]
  e <- x$n_entrants
  q <- x$response_model$q
  out <- data.frame(
    period = tab$period,
    visit = tab$visit,
    n_issued = e * tab$relative_take,
    n_resp = e * tab$relative_take * q[tab$visit]
  )
  out <- aggregate(
    cbind(n_issued, n_resp) ~ period + visit,
    data = out,
    FUN = sum
  )
  out[order(out$period, out$visit), ]
}

schedule_realized_overlap <- function(origin, x, lag) {
  active <- x$schedule[x$schedule$active, , drop = FALSE]
  key <- paste(active$selection, active$panel, sep = "|")
  first <- unique(key[active$period == origin])
  second <- unique(key[active$period == origin + lag])
  length(intersect(first, second)) / length(first)
}

test_that("an immediate schedule has stable components and a dense grid", {
  plan <- schedule_panel_plan()

  x <- design_schedule(
    plan,
    horizon = 6,
    horizon_policy = "continuing",
    refreshment = "entrant_register",
    frame_vintage = c(startup = "2026Q1", period_2 = "2026Q2"),
    rounding = "ceiling"
  )

  expect_s3_class(x, "svyplan_schedule")
  expect_identical(x$start, "immediate")
  expect_identical(x$horizon_policy, "continuing")
  expect_identical(x$equilibrium_from, 4L)
  expect_identical(x$selections$selection, c("startup", paste0("period_", 2:6)))
  expect_identical(x$selections$panels, c(4L, rep(1L, 5L)))
  expect_true(all(x$selections$status == "planned"))
  expect_identical(x$selections$frame_vintage[1:2], c("2026Q1", "2026Q2"))
  expect_true(all(is.na(x$selections$frame_vintage[-(1:2)])))

  period_n <- ceiling(plan$n_entrants)
  expect_equal(x$selections$planned_issue[[1L]], plan$n_in_sample)
  expect_equal(x$selections$operational_issue[[1L]], 4L * period_n)
  expect_equal(x$selections$panel_issue[[1L]], period_n)

  expect_identical(nrow(x$schedule), 6L * (4L + 5L))
  expect_true(all(is.na(x$schedule$life_stage[!x$schedule$active])))
  expect_false(anyNA(x$schedule$life_stage[x$schedule$active]))
  expect_true(all(x$schedule$equilibrium[x$schedule$period >= 4L]))
  expect_false(any(x$schedule$equilibrium[x$schedule$period < 4L]))

  # startup panels run from the one-period tail to the full life, as the
  # chart numbers them, each at its first interview and its own wave
  startup_1 <- x$schedule[
    x$schedule$selection == "startup" & x$schedule$period == 1L,
  ]
  expect_true(all(startup_1$active))
  expect_identical(startup_1$panel, 1:4)
  expect_identical(startup_1$life_length, 1:4)
  expect_identical(startup_1$life_stage, 4:1)
  expect_identical(startup_1$wave, 4:1)
  expect_identical(startup_1$visit, rep(1L, 4L))
  full_life <- x$schedule[
    x$schedule$selection == "startup" & x$schedule$panel == 4L & x$schedule$active,
  ]
  expect_identical(full_life$life_stage, 1:4)
  expect_identical(full_life$wave, 1:4)
  period_2 <- x$schedule[
    x$schedule$selection == "period_2" & x$schedule$active,
  ]
  expect_identical(period_2$period, 2:5)
  expect_identical(period_2$life_stage, 1:4)
  expect_identical(period_2$wave, 1:4)

  expect_gt(nrow(x$tail_commitments), 0L)
  expect_gt(min(x$tail_commitments$period), x$horizon)
  expect_lte(max(x$tail_commitments$period), x$horizon + x$life_length - 1L)
  expect_identical(x$panel_parameters$k, 4L)
  expect_identical(x$panel_parameters$r_min, 1L)
  expect_identical(x$panel_parameters$block_width, 8L)
})

test_that("generated schedules reconcile with the period-visit table where applicable", {
  for (start in c("immediate", "gradual")) {
    plan <- schedule_panel_plan(start)
    for (policy in c("continuing", "truncate_lives")) {
      x <- design_schedule(
        plan,
        6,
        policy,
        refreshment = "entrant_register",
        rounding = "ceiling"
      )
      got <- schedule_active_counts(x)
      got <- got[got$period <= 5L, ]
      want <- plan$period_visits[plan$period_visits$period <= 5L, ]
      rownames(got) <- NULL
      rownames(want) <- NULL
      expect_equal(got, want, tolerance = 1e-12, info = paste(start, policy))

      later <- schedule_active_counts(x)
      later <- later[later$period == 6L, -1L, drop = FALSE]
      settled <- plan$period_visits[
        plan$period_visits$period == 5L,
        -1L,
        drop = FALSE
      ]
      rownames(later) <- rownames(settled) <- NULL
      expect_equal(
        later,
        settled,
        tolerance = 1e-12,
        info = paste(start, policy, "after the start-up")
      )
    }
  }

  plan <- schedule_panel_plan("immediate")
  closed <- design_schedule(
    plan,
    6,
    "complete_lives",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )
  got <- schedule_active_counts(closed)
  got <- got[got$period <= 3L, ]
  want <- plan$period_visits[plan$period_visits$period <= 3L, ]
  rownames(got) <- rownames(want) <- NULL
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("overlap matches the life only from mature membership origins", {
  expected <- as.double(design_overlap("4"))

  immediate <- design_schedule(
    schedule_panel_plan("immediate"),
    8,
    "continuing",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )
  for (lag in seq_along(expected)) {
    origins <- seq_len(immediate$horizon - lag)
    got <- vapply(
      origins,
      schedule_realized_overlap,
      numeric(1L),
      x = immediate,
      lag = lag
    )
    expect_equal(got, rep(expected[[lag]], length(got)), tolerance = 1e-12)
  }

  gradual <- design_schedule(
    schedule_panel_plan("gradual"),
    8,
    "continuing",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )
  expect_equal(schedule_realized_overlap(1L, gradual, 1L), 1)
  for (lag in seq_along(expected)) {
    origins <- 4:(gradual$horizon - lag)
    if (length(origins) > 0L && all(origins <= gradual$horizon - lag)) {
      got <- vapply(
        origins,
        schedule_realized_overlap,
        numeric(1L),
        x = gradual,
        lag = lag
      )
      expect_equal(got, rep(expected[[lag]], length(got)), tolerance = 1e-12)
    }
  }
})

test_that("horizon policies separate issue, commitments and tail composition", {
  plan <- schedule_panel_plan()
  make <- function(policy) {
    design_schedule(
      plan,
      6,
      policy,
      refreshment = "entrant_register",
      rounding = "ceiling"
    )
  }
  continuing <- make("continuing")
  truncated <- make("truncate_lives")
  closed <- make("complete_lives")

  expect_equal(
    sum(continuing$issue$operational_issue),
    sum(truncated$issue$operational_issue)
  )
  expect_gt(
    sum(continuing$issue$operational_issue),
    sum(closed$issue$operational_issue)
  )
  expect_gt(nrow(continuing$tail_commitments), 0L)
  expect_identical(nrow(truncated$tail_commitments), 0L)
  expect_identical(nrow(closed$tail_commitments), 0L)

  cont_active <- continuing$schedule$active
  trunc_active <- truncated$schedule$active
  expect_identical(cont_active, trunc_active)
  expect_lt(
    truncated$schedule$life_length[
      truncated$schedule$selection == "period_6" &
        truncated$schedule$period == 6L
    ][[1L]],
    continuing$life_length
  )
  active_by_period <- tapply(closed$schedule$active, closed$schedule$period,
                             sum)
  expect_true(all(diff(tail(active_by_period, 4L)) < 0L))
  expect_true(is.na(closed$equilibrium_from))
})

test_that("a gradual start has one start-up panel and reaches composition at L", {
  x <- design_schedule(
    schedule_panel_plan("gradual"),
    6,
    "continuing",
    refreshment = "whole_vintage",
    rounding = "ceiling"
  )

  expect_identical(x$selections$panels[[1L]], 1L)
  expect_true(all(x$selections$frame_role[-1L] == "whole_vintage"))
  expect_identical(x$equilibrium_from, 4L)
  expect_true(is.na(x$panel_parameters$block_width))
})

test_that("inputs that cannot be scheduled are refused", {
  rotating <- schedule_panel_plan()
  fixed <- n_panel(
    n_prop(p = 0.5, moe = 0.03),
    retention = rep(0.9, 3),
    resp_rate = 0.75
  )
  startless <- n_panel(
    n_prop(p = 0.5, moe = 0.03),
    retention = rep(0.9, 3),
    resp_rate = 0.75,
    design = "rotating"
  )

  expect_error(
    design_schedule(
      fixed,
      6,
      "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "rotating"
  )
  expect_error(
    design_schedule(
      startless,
      6,
      "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "start"
  )
  expect_error(
    design_schedule(
      rotating,
      0,
      "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "horizon"
  )
  expect_error(
    design_schedule(
      rotating,
      3,
      "complete_lives",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "shorter"
  )
  expect_error(
    design_schedule(
      rotating,
      6,
      "continuing",
      refreshment = "entrant_register",
      rounding = "nearest"
    ),
    "rounding"
  )
  expect_error(
    design_schedule(
      rotating,
      6,
      "later",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "horizon_policy"
  )
  expect_error(
    design_schedule(
      rotating,
      6,
      "continuing",
      refreshment = "other",
      rounding = "ceiling"
    ),
    "refreshment"
  )
  expect_error(
    design_schedule(
      rotating,
      6,
      "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling",
      frame_vintage = c(nope = "x")
    ),
    "frame_vintage"
  )
})

test_that("horizon policy is explicit and rounding defaults to ceiling", {
  plan <- schedule_panel_plan()
  expect_error(
    design_schedule(
      plan,
      horizon_policy = "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "horizon"
  )
  expect_error(
    design_schedule(
      plan,
      6,
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "horizon_policy"
  )
  defaulted <- design_schedule(
    plan,
    6,
    "continuing",
    refreshment = "entrant_register"
  )
  expect_identical(defaulted$rounding$rule, "ceiling")
  expect_identical(
    defaulted,
    design_schedule(plan, 6, "continuing", refreshment = "entrant_register",
                    rounding = "ceiling")
  )
  expect_error(
    design_schedule(
      plan,
      6,
      "continuing",
      rounding = "ceiling"
    ),
    "refreshment"
  )
})

test_that("the schedule is validated at construction", {
  x <- design_schedule(
    schedule_panel_plan(),
    6,
    "continuing",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )

  expect_invisible(svyplan:::.validate_schedule_schema(unclass(x)))

  broken <- unclass(x)
  broken$schedule$life_stage[!broken$schedule$active][[1L]] <- 1L
  expect_error(
    svyplan:::.validate_schedule_schema(broken),
    "dense schedule grid is inconsistent"
  )
  broken <- unclass(x)
  broken$schedule$cohort[broken$schedule$active][[1L]] <- 99L
  expect_error(
    svyplan:::.validate_schedule_schema(broken),
    "does not match the selection manifest|count the interviews of its cohort"
  )

  broken <- unclass(x)
  broken$tail_commitments$active <- TRUE
  expect_error(
    svyplan:::.validate_schedule_schema(broken),
    "do not have the schedule's columns"
  )

  mutations <- list(
    panel_number = function(z) {
      z$schedule$panel[[1L]] <- 99L
      z
    },
    rounding_rule = function(z) {
      z$rounding$rule <- "floor"
      z
    },
    disjointness = function(z) {
      z$requires_disjoint_entrants <- NA
      z
    },
    tail_stage = function(z) {
      z$tail_commitments$life_stage[[1L]] <- 99L
      z
    },
    in_sample = function(z) {
      z$n_in_sample <- 1
      z
    },
    relative_take = function(z) {
      z$schedule$relative_take[] <- 0.5
      z
    },
    response = function(z) {
      z$response_model$q[[2L]] <- z$response_model$q[[1L]]
      z
    },
    eligibility = function(z) {
      z$response_model$eligibility_rate <- 0.5
      z
    },
    component_issue = function(z) {
      z$selections$planned_issue[[1L]] <- 1
      z
    },
    issue_profile = function(z) {
      z$issue$selection[[2L]] <- "startup"
      z
    },
    overlap = function(z) {
      z$overlap$overlap[[1L]] <- 0
      z
    },
    equilibrium = function(z) {
      z$equilibrium_from <- 1L
      z
    },
    wave = function(z) {
      z$schedule$wave[z$schedule$active][[1L]] <- 9L
      z
    },
    rotation_group = function(z) {
      z$schedule$rotation_group[[1L]] <- 0L
      z
    },
    panel_parameters = function(z) {
      z$panel_parameters$r_min <- 2L
      z
    }
  )
  for (nm in names(mutations)) {
    expect_error(
      svyplan:::.validate_schedule_schema(mutations[[nm]](unclass(x))),
      "internal schedule defect",
      info = nm
    )
  }
})

test_that("a panel with eligibility below one gives a valid schedule", {
  plan <- n_panel(
    n_prop(p = 0.5, moe = 0.03),
    retention = c(0.9, 0.9, 0.9),
    resp_rate = 0.8,
    eligibility_rate = 0.9,
    design = "rotating",
    start = "immediate"
  )
  x <- design_schedule(plan, horizon = 6, horizon_policy = "continuing",
                       refreshment = "entrant_register", rounding = "ceiling")
  expect_equal(x$response_model$q, plan$waves$q)
  expect_equal(x$response_model$eligibility_rate, 0.9)
  expect_invisible(svyplan:::.validate_schedule_schema(unclass(x)))
})

test_that("operational counts remain whole-valued above the integer range", {
  large <- as.double(.Machine$integer.max) + 1000
  plan <- prec_panel(
    large,
    n_prop(p = 0.5, moe = 0.03),
    retention = rep(0.9, 3),
    design = "rotating",
    start = "immediate"
  )

  expect_warning(
    x <- design_schedule(
      plan,
      6,
      "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    NA
  )
  expect_type(x$selections$panel_issue, "double")
  expect_equal(x$selections$panel_issue[[1L]], large)
  expect_equal(x$selections$operational_issue[[1L]], 4 * large)
  expect_true(all(
    x$selections$operational_issue == floor(x$selections$operational_issue)
  ))
  expect_invisible(svyplan:::.validate_schedule_schema(unclass(x)))
})

test_that("schedule replacement cannot create cross-field contradictions", {
  x <- design_schedule(
    schedule_panel_plan(),
    6,
    "continuing",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )

  expect_error(
    {
      x$rounding <- list(rule = "floor")
    },
    "cannot be modified in place"
  )
  expect_error(
    {
      x[["n_in_sample"]] <- 1
    },
    "cannot be modified in place"
  )
  expect_error(
    {
      x["schedule"] <- list(NULL)
    },
    "cannot be modified in place"
  )
  expect_error(
    {
      x$schedule$panel[[1L]] <- 99L
    },
    "cannot be modified in place"
  )
})

test_that("schedule methods lead with the issue profile", {
  x <- design_schedule(
    schedule_panel_plan(),
    6,
    "continuing",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )

  expect_identical(as.data.frame(x), x$issue)
  expect_match(format(x), "svyplan_schedule.*immediate.*continuing")
  out <- capture.output(print(x))
  expect_match(out[[1L]], "Longitudinal design schedule")
  expect_true(any(grepl("tail commitments", out, fixed = TRUE)))
  expect_error(as.data.frame(x, view = "schedule"), "unused argument")
})

test_that("a schedule survives an RDS round trip", {
  x <- design_schedule(
    schedule_panel_plan(),
    6,
    "continuing",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))

  saveRDS(x, path)
  restored <- readRDS(path)

  expect_identical(restored, x)
  expect_s3_class(restored, "svyplan_schedule")
  expect_invisible(svyplan:::.validate_schedule_schema(unclass(restored)))
  expect_identical(as.data.frame(restored), restored$issue)
  expect_match(format(restored), "svyplan_schedule.*immediate.*continuing")
})

## S13. The printed block does not grow with the reporting horizon

test_that("print states the issue profile as runs, not one row per period", {
  mk <- function(h, policy = "continuing") {
    design_schedule(schedule_panel_plan(), h, policy,
                    refreshment = "entrant_register", rounding = "ceiling")
  }
  # A schedule's length is set by the reporting window, not by the design, so
  # the block must be the same size at every horizon.
  for (h in c(6L, 12L, 24L, 52L)) {
    out <- capture.output(print(mk(h)))
    expect_length(out, 7L)
    expect_lt(max(nchar(out)), 80L)
    expect_true(any(grepl(sprintf("periods 2-%d\\)", h), out)))
  }
  # The periods are gone from print and whole in summary.
  out <- capture.output(print(mk(24L)))
  expect_false(any(grepl("period_7", out, fixed = TRUE)))
  expect_identical(nrow(summary(mk(24L))$issue), 24L)
})

test_that("each horizon policy states its own profile", {
  mk <- function(h, policy) {
    design_schedule(schedule_panel_plan(), h, policy,
                    refreshment = "entrant_register", rounding = "ceiling")
  }
  # Taking no new cohort is a fact about the schedule, not a gap in the report.
  closed <- capture.output(print(mk(10L, "complete_lives")))
  expect_true(any(grepl("no new cohort from period 8", closed, fixed = TRUE)))
  expect_true(any(grepl("periods 2-7", closed, fixed = TRUE)))
  expect_length(closed, 7L)
  for (policy in c("continuing", "truncate_lives")) {
    out <- capture.output(print(mk(12L, policy)))
    expect_false(any(grepl("no new cohort", out, fixed = TRUE)))
    expect_length(out, 7L)
  }
})

test_that("the life line keeps the branch where equilibrium is not reached", {
  short <- design_schedule(schedule_panel_plan(), 2L,
                           "continuing", refreshment = "entrant_register",
                           rounding = "ceiling")
  expect_true(is.na(short$equilibrium_from))
  out <- capture.output(print(short))
  expect_identical(out[2L], "life: 4 periods, 4 in sample each period")
  expect_identical(out[3L],
                   "horizon: 2 periods, equilibrium not reached in the window")
})

test_that("summary carries every table print no longer shows", {
  x <- design_schedule(schedule_panel_plan(), 6L,
                       "continuing", refreshment = "entrant_register",
                       rounding = "ceiling")
  sm <- summary(x)
  expect_s3_class(sm, "summary.svyplan_schedule")
  expect_identical(sm$issue, x$issue)
  expect_identical(sm$selections, x$selections)
  expect_identical(sm$activity, x$schedule)
  expect_identical(sm$tail_commitments, x$tail_commitments)
  out <- capture.output(print(sm))
  expect_true(any(grepl("Issue profile", out, fixed = TRUE)))
  expect_true(any(grepl("Overlap the rotation produces", out, fixed = TRUE)))
  # Whole units, and the same ones the issue profile shows: the standing
  # sample is the cohorts times the take each was rounded to.
  expect_true(any(grepl(
    sprintf("in sample across %d cohorts", x$n_cohorts), out, fixed = TRUE
  )))
  shown <- ceiling(x$n_entrants) * x$n_cohorts
  expect_true(any(grepl(sprintf("%d in sample", shown), out)))
  expect_identical(shown, x$issue$operational_issue[[1L]])
})

## A rotation, and nothing that merely describes one

test_that("the schedule lays out the rotation the panel was planned on", {
  plan <- n_panel(n_prop(p = 0.5, moe = 0.03), retention = c(0.9, 0.85, 0.9),
                  resp_rate = 0.75, design = "rotating", start = "gradual",
                  rotation = "2-2-2")
  x <- design_schedule(plan, 6, "continuing",
                       refreshment = "entrant_register", rounding = "ceiling")
  expect_identical(x$life_spec, c(1, 1, 0, 0, 1, 1))
  expect_identical(x$life_length, 6L)
  expect_identical(x$n_cohorts, 4L)
  expect_error(
    design_schedule(plan, design_rotation("2-2-2"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "horizon"
  )
  stripped <- plan
  stripped$rotation <- NULL
  expect_error(
    design_schedule(stripped, 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "must carry its rotation"
  )
})

test_that("the manifest computes its own overlap from the rotation", {
  out <- design_schedule(
    schedule_panel_plan(), 6, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
  )
  expect_identical(out$overlap, as.data.frame(design_overlap("4")))
})

## A life with a gap

.gapped_schedule <- function(start = "immediate", horizon = 8,
                             policy = "continuing") {
  design_schedule(
    n_panel(n_prop(p = 0.5, moe = 0.03), retention = c(0.9, 0.85, 0.9),
            resp_rate = 0.75, design = "rotating", start = start,
            rotation = "2-2-2"),
    horizon, policy, refreshment = "entrant_register", rounding = "ceiling"
  )
}

test_that("a 2-2-2 startup holds six panels and interviews four", {
  x <- .gapped_schedule()
  s <- x$schedule
  expect_identical(x$panel_parameters$k, 6L)
  expect_equal(as.vector(tapply(s$active, s$period, sum)), rep(4L, 8L))
  expect_equal(as.vector(tapply(!is.na(s$life_stage), s$period, sum)),
               rep(6L, 8L))
  expect_equal(x$n_held, 6 * x$n_entrants)
  st <- s[s$selection == "startup", ]
  first <- vapply(1:6, function(j) min(st$period[st$panel == j & st$active]),
                  numeric(1L))
  expect_identical(first, c(1, 1, 2, 3, 1, 1))
  expect_identical(st$visit[st$panel == 6L], c(1L, 2L, NA, NA, 3L, 4L, NA, NA))
  # NISR (2024), Annex 1: the startup panels open at waves D, C, C, C, B, A
  expect_identical(
    vapply(1:6, function(j) st$wave[st$panel == j & st$active][[1L]],
           integer(1L)),
    c(4L, 3L, 3L, 3L, 2L, 1L)
  )
  # the panels covering stages 3 and 4 are held, not absent, before
  # their first interview, at the stages of the pattern they cover
  expect_identical(st$life_stage[st$panel == 4L][1:2], 3:4)
  expect_false(any(st$active[st$panel == 4L][1:2]))
  intake <- s[s$selection == "period_2", ]
  expect_identical(intake$visit, c(NA, 1L, 2L, NA, NA, 3L, 4L, NA))
  expect_identical(intake$wave, c(NA, 1L, 2L, NA, NA, 3L, 4L, NA))
  expect_identical(x$equilibrium_from, 6L)
  expect_true(all(!is.na(x$tail_commitments$visit)))
  expect_true(all(x$schedule$subround == 1L))
})

test_that("a gapped schedule matches an enumeration of its grid", {
  specs <- c("3", "5", "2-2-2", "4-8-4", "2-1-3", "1-0-1", "1-1-0-1")
  for (spec in specs) for (start in c("gradual", "immediate")) {
    w <- as.double(design_rotation(spec))
    L <- length(w)
    v <- sum(w > 0)
    plan <- n_panel(n_prop(p = 0.5, moe = 0.03),
                    retention = seq(0.9, 0.8, length.out = v - 1L),
                    resp_rate = 0.8, design = "rotating", start = start,
                    rotation = spec)
    for (h in unique(c(2L, L, L + 3L))) {
      for (policy in c("continuing", "truncate_lives", "complete_lives")) {
        if (policy == "complete_lives" && h < L) next
        first <- max(vapply(seq_len(L), function(j) {
          which(w[j:L] > 0)[[1L]]
        }, integer(1L)))
        if (policy == "truncate_lives" && start == "immediate" && h < first) {
          next
        }
        lab <- paste(spec, start, h, policy)
        x <- design_schedule(plan, h, policy,
                             refreshment = "entrant_register",
                             rounding = "ceiling")
        last_entry <- if (policy == "complete_lives") h - L + 1L else h
        o <- startup_oracle(w, start, last_entry + L)
        o$entry <- ifelse(substr(o$id, 1, 1) == "s", 1L,
                          as.integer(substring(o$id, 2)))
        o <- o[o$entry <= last_entry, ]
        if (policy == "truncate_lives") o <- o[o$period <= h, ]
        inside <- o[o$period <= h, ]
        got <- x$schedule[x$schedule$active, ]
        expect_identical(
          as.vector(table(factor(got$period, levels = seq_len(h)))),
          as.vector(table(factor(inside$period, levels = seq_len(h)))),
          label = lab
        )
        expect_identical(sort(got$visit), sort(inside$visit), label = lab)
        tail_o <- o[o$period > h, ]
        expect_identical(nrow(x$tail_commitments),
                         if (policy == "continuing") nrow(tail_o) else 0L,
                         label = lab)
        if (policy == "truncate_lives") {
          held <- x$schedule[!is.na(x$schedule$life_stage), ]
          ends <- stats::aggregate(period ~ selection + panel, held, max)
          last_int <- stats::aggregate(period ~ selection + panel,
                                       x$schedule[x$schedule$active, ], max)
          m <- merge(ends, last_int, by = c("selection", "panel"))
          expect_identical(m$period.x, m$period.y, label = lab)
        }
      }
    }
  }
})

test_that("truncate_lives refuses a horizon before a startup panel's first interview", {
  expect_error(.gapped_schedule(horizon = 2, policy = "truncate_lives"),
               "which the last of them has at period 3")
  expect_s3_class(.gapped_schedule(horizon = 3, policy = "truncate_lives"),
                  "svyplan_schedule")
  expect_s3_class(.gapped_schedule(start = "gradual", horizon = 2,
                                   policy = "truncate_lives"),
                  "svyplan_schedule")
})

test_that("the version-3 validator refuses each broken rule", {
  x <- unclass(.gapped_schedule())
  held <- which(!x$schedule$active & !is.na(x$schedule$life_stage))[[1L]]
  live <- which(x$schedule$active)
  mutations <- list(
    unknown_field = function(z) {
      z$schema_version <- 4L
      z
    },
    visit_skips = function(z) {
      z$schedule$visit[live[[2L]]] <- z$schedule$visit[live[[2L]]] + 1L
      z
    },
    visit_on_inactive = function(z) {
      z$schedule$visit[held] <- 1L
      z
    },
    held_without_stage = function(z) {
      z$schedule$life_stage[held] <- NA_integer_
      z
    },
    tail_not_interview = function(z) {
      z$tail_commitments$visit[[1L]] <- NA_integer_
      z
    },
    period = function(z) {
      z$schedule$period[[1L]] <- 2L
      z
    },
    subround = function(z) {
      z$schedule$subround[[1L]] <- 2L
      z
    },
    wave_on_held = function(z) {
      z$schedule$wave[held] <- 1L
      z
    },
    n_held = function(z) {
      z$n_held <- z$n_in_sample
      z
    },
    n_cohorts = function(z) {
      z$n_cohorts <- z$life_length
      z
    }
  )
  for (nm in names(mutations)) {
    expect_error(
      svyplan:::.validate_schedule_schema(mutations[[nm]](x)),
      "internal schedule defect",
      info = nm
    )
  }
  expect_invisible(svyplan:::.validate_schedule_schema(x))
})

test_that("print names a gapped rotation and summary the units held", {
  x <- .gapped_schedule()
  out <- capture.output(print(x))
  expect_true(any(out == "rotation: 2 in, 2 out, 2 in"))
  s <- capture.output(print(summary(x)))
  expect_true(any(grepl(
    sprintf("held: %d across 6 cohorts, 2 of them out of sample",
            6L * x$selections$panel_issue[[1L]]), s, fixed = TRUE
  )))
  plain <- design_schedule(schedule_panel_plan(), 6, "continuing",
                           refreshment = "entrant_register",
                           rounding = "ceiling")
  expect_false(any(grepl("rotation:|held:",
                         capture.output(print(summary(plain))))))
})

test_that("summary counts the units the rounded panels issue", {
  plan <- n_panel(n_prop(p = 0.5, moe = 0.02, deff = 1.8),
                  retention = rep(0.9, 4), resp_rate = 0.8,
                  design = "rotating", start = "immediate")
  gapped <- n_panel(n_prop(p = 0.5, moe = 0.02, deff = 1.8),
                    retention = rep(0.9, 3), resp_rate = 0.8,
                    design = "rotating", start = "immediate",
                    rotation = "2-2-2")
  mk <- function(p, ...) {
    design_schedule(p, horizon = 6, horizon_policy = "continuing",
                    refreshment = "entrant_register", rounding = "ceiling",
                    ...)
  }
  cases <- list(
    take = mk(plan, take = 14),
    subrounds = mk(plan, subrounds = 7),
    both = mk(plan, take = 14, retain_psus = TRUE, subrounds = 3),
    gapped = mk(gapped, take = 14, subrounds = 3)
  )
  for (nm in names(cases)) {
    x <- cases[[nm]]
    e <- x$selections$panel_issue[[1L]] * c(x$take, 1)[[1L + is.na(x$take)]] *
      x$subrounds
    expect_false(e == ceiling(x$n_entrants), info = nm)
    s <- capture.output(print(summary(x)))
    expect_true(any(s == sprintf(
      "entrants: %d per period, %d in sample across %d cohorts",
      e, e * x$n_cohorts, x$n_cohorts
    )), info = nm)
    if (x$life_length > x$n_cohorts) {
      expect_true(any(startsWith(s, sprintf("held: %d across", e *
                                              x$life_length))), info = nm)
    }
  }
  s <- capture.output(print(summary(cases$both)))
  expect_true(any(s == paste("entrants: 1344 per period, 6720 in sample",
                             "across 5 cohorts")))
})

test_that("the overlap check runs at the mature origins of a gapped life", {
  x <- .gapped_schedule(horizon = 12)
  want <- as.double(design_overlap("2-2-2"))
  expect_invisible(svyplan:::.check_schedule_overlap(x, want))
  expect_error(svyplan:::.check_schedule_overlap(x, rev(want)),
               "overlap at origin 1 and lag 1")
})

## PSUs, retained PSUs and subrounds

.c_plan <- function(spec = "5", start = "immediate") {
  w <- as.double(design_rotation(spec))
  n_panel(n_prop(p = 0.5, moe = 0.02, deff = 1.8),
          retention = seq(0.9, 0.8, length.out = sum(w > 0) - 1L),
          resp_rate = 0.8, design = "rotating", start = start,
          rotation = spec)
}

.c_schedule <- function(plan, horizon, policy = "continuing", ...) {
  design_schedule(plan, horizon, policy, refreshment = "entrant_register",
                  rounding = "ceiling", ...)
}

test_that("take counts PSUs and the units drawn for each cohort", {
  plan <- .c_plan()
  x <- .c_schedule(plan, 8, take = 14)
  e <- plan$n_entrants
  psu <- x$selections$panel_issue[[1L]]
  expect_identical(psu, ceiling(e / 14))
  expect_gte(psu * 14, e)
  expect_lt(psu * 14 - e, 14)
  expect_equal(x$selections$planned_issue, e / 14 * x$selections$panels)
  expect_identical(x$selections$operational_issue, psu * x$selections$panels)
  expect_identical(x$rounding$level, "psu")
  expect_identical(x$take, 14)
  opened <- tapply(x$schedule$active & x$schedule$visit %in% 1L,
                   x$schedule$period, sum)
  expect_equal(x$issue$cohort_issue, as.vector(opened) * psu * 14)
  expect_error(.c_schedule(plan, 8, take = 2.5), "'take' must be NULL")
  expect_error(.c_schedule(plan, 8, take = 14, retain_psus = "open"),
               "TRUE or FALSE")
  expect_error(.c_schedule(plan, 8, retain_psus = TRUE), "needs 'take'")
})

test_that("retained PSUs keep their groups and take new cohorts", {
  x <- .c_schedule(.c_plan(), 12, take = 1, retain_psus = TRUE)
  expect_identical(x$selections$selection, "startup")
  expect_identical(x$selections$panels, 5L)
  expect_true(x$retain_psus)
  per_period <- tapply(x$schedule$active, x$schedule$period, sum)
  expect_equal(as.vector(per_period), rep(5, 12))
  opening <- x$issue$cohort_issue / x$selections$panel_issue[[1L]]
  expect_identical(opening, c(5, rep(1, 11)))
  grad <- .c_schedule(.c_plan(start = "gradual"), 12, take = 1, retain_psus = TRUE)
  expect_identical(grad$selections$selection, c("startup", paste0("period_", 2:5)))
  expect_identical(grad$issue$cohort_issue / grad$selections$panel_issue[[1L]],
                   rep(1, 12))
  gap <- .c_schedule(.c_plan("2-2-2"), 10, take = 1, retain_psus = TRUE)
  expect_identical(gap$selections$panels, 6L)
  expect_equal(as.vector(tapply(gap$schedule$active, gap$schedule$period,
                                sum)), rep(4, 10))
  expect_equal(as.vector(tapply(!is.na(gap$schedule$life_stage),
                                gap$schedule$period, sum)), rep(6, 10))
  full <- gap$schedule$panel == 6L & gap$schedule$active
  expect_true(any(diff(gap$schedule$visit[full]) < 0))
  # a new cohort restarts the visits and keeps the waves of the pattern
  expect_identical(gap$schedule$wave[full], gap$schedule$visit[full])
})

test_that("retained and replaced PSUs field the same units", {
  specs <- c("3", "5", "2-2-2", "4-8-4", "2-1-3", "1-0-1", "1-1-0-1")
  for (spec in specs) for (start in c("gradual", "immediate")) {
    plan <- .c_plan(spec, start)
    L <- length(as.double(plan$rotation))
    for (h in unique(c(L, L + 3L, 2L * L + 1L))) {
      for (policy in c("continuing", "truncate_lives", "complete_lives")) {
        lab <- paste(spec, start, h, policy)
        one <- .c_schedule(plan, h, policy)
        open <- .c_schedule(plan, h, policy, take = 1, retain_psus = TRUE)
        a <- one$schedule[one$schedule$active, ]
        b <- open$schedule[open$schedule$active, ]
        expect_identical(tabulate(b$period, h), tabulate(a$period, h),
                         label = lab)
        expect_identical(sort(b$visit), sort(a$visit), label = lab)
        expect_identical(open$issue$cohort_issue, one$issue$cohort_issue,
                         label = lab)
        expect_identical(nrow(open$tail_commitments),
                         nrow(one$tail_commitments), label = lab)
        expect_identical(open$equilibrium_from, one$equilibrium_from,
                         label = lab)
        expect_identical(open$issue$equilibrium, one$issue$equilibrium,
                         label = lab)
        expect_lte(sum(open$selections$panels), L)
      }
    }
  }
})

test_that("subrounds split each group into its own months", {
  for (start in c("immediate", "gradual")) for (retain in c(FALSE, TRUE)) {
    plan <- .c_plan("5", start)
    r <- .c_schedule(plan, 4, take = 1, retain_psus = retain)
    x <- .c_schedule(plan, 4, take = 1, retain_psus = retain, subrounds = 3)
    lab <- paste(start, retain)
    s <- x$schedule
    expect_identical(x$horizon, 4L, label = lab)
    expect_identical(x$subrounds, 3L, label = lab)
    expect_identical(x$panel_parameters$k,
                     if (start == "immediate") 15L else 3L, label = lab)
    expect_true(all(s$subround %in% 1:3), label = lab)
    expect_identical(x$selections$entry_period, r$selections$entry_period,
                     label = lab)
    expect_identical(x$selections$panels, r$selections$panels * 3L,
                     label = lab)
    on <- s[s$active, ]
    expect_true(all(on$subround - 1L == (on$panel - 1L) %% 3L),
                label = lab)
    # collapsing the subround panels of a group gives back the grid
    s$group <- (s$panel - 1L) %/% 3L + 1L
    back <- stats::aggregate(active ~ period + selection + group, s, any)
    want <- stats::aggregate(active ~ period + selection + panel, r$schedule,
                             any)
    m <- merge(back, want, by.x = c("period", "selection", "group"),
               by.y = c("period", "selection", "panel"))
    expect_identical(nrow(m), nrow(want), label = lab)
    expect_identical(m$active.x, m$active.y, label = lab)
    stage <- unique(s[!is.na(s$life_stage),
                      c("period", "selection", "group", "life_stage")])
    expect_identical(
      unique(s[c("selection", "group", "rotation_group")])$rotation_group,
      unique(r$schedule[c("selection", "panel", "rotation_group")])$rotation_group,
      label = lab
    )
    expect_identical(anyDuplicated(stage[c("period", "selection", "group")]), 0L,
                     label = lab)
    expect_identical(sort(on$visit),
                     sort(rep(r$schedule$visit[r$schedule$active], each = 3L)),
                     label = lab)
    expect_identical(x$equilibrium_from, r$equilibrium_from, label = lab)
  }
})

test_that("a quarterly design with monthly fieldwork and kept PSUs runs as in samplyr's probe", {
  x <- .c_schedule(.c_plan("5"), 4, take = 14, retain_psus = TRUE,
                   subrounds = 3)
  s <- x$schedule
  expect_identical(s$active, s$subround - 1L == (s$panel - 1L) %% 3L)
  opening <- s[s$active & s$visit == 1L, ]
  expect_true(all(opening$subround - 1L == (opening$panel - 1L) %% 3L))
  expect_identical(x$rounding$level, "psu_and_subround")
  e <- x$n_entrants
  psu <- x$selections$panel_issue[[1L]]
  expect_identical(psu, ceiling(e / (14 * 3)))
  expect_lt(psu * 14 * 3 - e, 14 * 3)
  out <- capture.output(print(x))
  expect_true(any(out == "take: 14 units per PSU per cohort"))
  expect_true(any(
    out == "PSUs: retained through the horizon, a new cohort at each return"
  ))
  expect_true(any(out == "life: 5 periods, 5 in sample each period"))
  expect_true(any(grepl("^horizon: 4 periods, ", out)))
  expect_true(any(out == "subrounds: 3 per period"))
  expect_true(any(grepl("no new PSUs from period 2", out, fixed = TRUE)))
})

test_that("subrounds split any rotation and keep its reporting grid", {
  specs <- c("3", "2-2-2", "4-8-4", "1-0-1")
  for (spec in specs) for (start in c("gradual", "immediate")) {
    for (s in c(1L, 3L, 13L)) for (retain in c(FALSE, TRUE)) {
      plan <- .c_plan(spec, start)
      L <- length(as.double(plan$rotation))
      lab <- paste(spec, start, s, retain)
      r <- .c_schedule(plan, L + 2L, take = 1, retain_psus = retain)
      x <- .c_schedule(plan, L + 2L, take = 1, retain_psus = retain,
                       subrounds = s)
      if (s == 1L) {
        expect_identical(x, r, label = lab)
        next
      }
      on <- x$schedule[x$schedule$active, ]
      want <- r$schedule[r$schedule$active, ]
      expect_identical(nrow(on), nrow(want) * s, label = lab)
      expect_equal(tabulate(on$period, L + 2L) / s,
                   tabulate(want$period, L + 2L), label = lab)
      expect_identical(nrow(x$tail_commitments),
                       nrow(r$tail_commitments) * s, label = lab)
    }
  }
})

test_that("the validator refuses broken PSU and subround rules", {
  x <- unclass(.c_schedule(.c_plan("5"), 4, take = 14, retain_psus = TRUE,
                           subrounds = 3))
  one <- unclass(.c_schedule(.c_plan("5"), 8))
  mutations <- list(
    retain_psus = function(z) {
      z$retain_psus <- FALSE
      z
    },
    subrounds = function(z) {
      z$subrounds <- 2L
      z
    },
    take = function(z) {
      z$take <- 10
      z
    },
    level = function(z) {
      z$rounding$level <- "psu"
      z
    },
    period = function(z) {
      z$schedule$period <- z$schedule$subround
      z
    },
    subround = function(z) {
      z$schedule$subround <- 1L
      z
    }
  )
  for (nm in names(mutations)) {
    expect_error(
      svyplan:::.validate_schedule_schema(mutations[[nm]](x)),
      "internal schedule defect", info = nm
    )
  }
  restart <- one
  at <- which(restart$schedule$selection == "startup" &
                restart$schedule$panel == 5L & restart$schedule$active)
  restart$schedule$visit[at[[3L]]] <- 1L
  expect_error(svyplan:::.validate_schedule_schema(restart),
               "does not count the interviews of its cohort")
})

## The schedule and the cell table describe each interview alike

test_that("the schedule's interviews are the cell table's, column for column", {
  for (spec in c("6", "2-2-2", "4-8-4")) {
    for (start in c("gradual", "immediate")) for (s in c(1L, 3L)) {
      rot <- design_rotation(spec)
      L <- length(rot)
      plan <- n_panel(n_prop(p = 0.5, moe = 0.05),
                      retention = rep(1, sum(rot > 0) - 1L),
                      design = "rotating", start = start, rotation = rot)
      x <- design_schedule(plan, horizon = 14, horizon_policy = "continuing",
                           refreshment = "entrant_register",
                           rounding = "ceiling", subrounds = s)
      a <- x$schedule[x$schedule$active, ]
      a <- a[(a$panel - 1L) %% s == 0L, ]
      reporting <- (a$panel - 1L) %/% s + 1L
      d <- as.data.frame(rot, start = start, horizon = 14)
      lab <- paste(spec, start, s)
      # every chart panel is one cohort, and a cohort holds s panels
      expect_identical(a$cohort, reporting, label = lab)
      expect_identical(
        sort(paste(a$period, reporting, a$rotation_group, a$wave, a$visit,
                   a$life_stage)),
        sort(paste(d$period, d$panel, d$rotation_group, d$wave, d$visit,
                   d$life_stage)),
        label = lab
      )
    }
  }
})

## Selections, cohorts and retained PSUs

test_that("a cohort is the units at the same wave at every period", {
  for (spec in c("6", "2-2-2", "4-8-4")) {
    for (start in c("gradual", "immediate")) for (s in c(1L, 3L)) {
      rot <- design_rotation(spec)
      L <- length(rot)
      plan <- n_panel(n_prop(p = 0.5, moe = 0.05),
                      retention = rep(1, sum(rot > 0) - 1L),
                      design = "rotating", start = start, rotation = rot)
      x <- design_schedule(plan, horizon = 12, horizon_policy = "continuing",
                           refreshment = "entrant_register",
                           rounding = "ceiling", subrounds = s)
      lab <- paste(spec, start, s)
      held <- x$schedule[!is.na(x$schedule$life_stage), ]
      offset <- if (start == "immediate") L - 1L else 0L
      expect_identical(held$cohort,
                       held$period - held$life_stage + 1L + offset,
                       label = lab)
      # a cohort keeps one wave per period, and holds s panels
      w <- unique(held[held$active, c("period", "cohort", "wave")])
      expect_identical(anyDuplicated(w[c("period", "cohort")]), 0L,
                       label = lab)
      per <- tapply(held$panel, held$cohort, function(p) length(unique(p)))
      expect_true(all(per == s), label = lab)
      start_up <- unique(held$cohort[held$selection == "startup"])
      expect_identical(length(start_up),
                       if (start == "immediate") L else 1L, label = lab)
      expect_identical(anyDuplicated(x$schedule[c("period", "subround",
                                                  "panel")]), 0L,
                       label = lab)
    }
  }
})

test_that("retained PSUs take a new cohort at wave 1 on the Rwanda plan", {
  plan <- prec_panel(1656, target = n_prop(p = 0.15, cv = 0.12),
                     retention = rep(1, 3), resp_rate = 0.9,
                     design = "rotating", start = "immediate",
                     rotation = "2-2-2")
  mk <- function(retain) {
    design_schedule(plan, horizon = 8, horizon_policy = "continuing",
                    refreshment = "whole_vintage", rounding = "ceiling",
                    take = 12, retain_psus = retain)
  }
  open <- mk(TRUE)
  one <- mk(FALSE)
  expect_identical(open$selections$selection, "startup")
  p6 <- open$schedule[open$schedule$panel == 6L, ]
  expect_identical(p6$wave, c(1L, 2L, NA, NA, 3L, 4L, 1L, 2L))
  expect_identical(p6$visit, c(1L, 2L, NA, NA, 3L, 4L, 1L, 2L))
  expect_identical(p6$cohort, c(rep(6L, 6L), 12L, 12L))
  expect_identical(open$schedule$rotation_group, open$schedule$panel)
  expect_identical(open$issue$cohort_issue, one$issue$cohort_issue)
  # the cohorts are the ones a schedule replacing its PSUs would field
  key <- function(z) {
    a <- z$schedule[z$schedule$active, ]
    sort(paste(a$period, a$cohort, a$wave, a$visit))
  }
  expect_identical(key(open), key(one))
})

test_that("a schedule names selections, cohorts and retained PSUs", {
  x <- .c_schedule(.c_plan("2-2-2"), 8, take = 14, retain_psus = TRUE)
  for (nm in c("start-up_policy", "tenures", "components")) {
    expect_false(nm %in% names(unclass(x)), label = nm)
  }
  expect_false("tenure_issue" %in% names(x$issue))
  expect_true("in_sample" %in% names(x$overlap))
  expect_false(any(grepl("intake_", x$selections$selection, fixed = TRUE)))
  expect_error(
    .c_schedule(.c_plan("2-2-2"), 8, tenures = "open"),
    "unused argument"
  )
})

## Rounding must not hold more units than the population

test_that("a rounded schedule beyond the finite population is refused", {
  panel <- n_panel(prec_prop(n = 90, p = 0.5, N = 100), retention = 1,
                   design = "rotating", start = "immediate")
  expect_error(
    design_schedule(panel, horizon = 3, horizon_policy = "continuing",
                    refreshment = "entrant_register", take = 30),
    "holds 120 frame units at period 1 \\(2 panels of 2 PSUs of 30 units\\), more than the population of 100"
  )
  expect_error(
    design_schedule(panel, horizon = 3, horizon_policy = "continuing",
                    refreshment = "entrant_register", take = 30,
                    retain_psus = TRUE),
    "more than the population of 100"
  )
  sub <- n_panel(prec_prop(n = 79.51, p = 0.5, N = 100), retention = 1,
                 design = "rotating", start = "immediate")
  expect_error(
    design_schedule(sub, horizon = 3, horizon_policy = "continuing",
                    refreshment = "entrant_register", subrounds = 30),
    "holds 120 frame units"
  )
  # a gap holds units too: 2-2-2 holds six cohorts, four of them in sample
  gap <- n_panel(prec_prop(n = 64, p = 0.5, N = 100), retention = c(1, 1, 1),
                 design = "rotating", rotation = "2-2-2", start = "immediate")
  expect_error(
    design_schedule(gap, horizon = 8, horizon_policy = "continuing",
                    refreshment = "entrant_register", take = 17),
    "holds 102 frame units"
  )
  # control: rounding up to the population itself is allowed
  ok <- n_panel(prec_prop(n = 96, p = 0.5, N = 100), retention = 1,
                design = "rotating", start = "immediate")
  x <- design_schedule(ok, horizon = 3, horizon_policy = "continuing",
                       refreshment = "entrant_register", take = 25)
  expect_equal(x$selections$operational_issue[[1L]] * 25, 100)
  # infinite N never triggers
  inf <- n_panel(prec_prop(n = 90, p = 0.5), retention = 1,
                 design = "rotating", start = "immediate")
  expect_s3_class(
    design_schedule(inf, horizon = 3, horizon_policy = "continuing",
                    refreshment = "entrant_register", take = 30),
    "svyplan_schedule"
  )
})

test_that("summary() reports the inflation rounding adds", {
  panel <- n_panel(prec_prop(n = 79.51, p = 0.5, N = 1000), retention = 1,
                   design = "rotating", start = "immediate")
  x <- design_schedule(panel, horizon = 3, horizon_policy = "continuing",
                       refreshment = "entrant_register", subrounds = 30)
  out <- capture.output(print(summary(x)))
  expect_true(any(grepl(
    "^rounding: ceiling at subround level, 60 entrants per period against the plan's 39.76 \\(\\+50.9%\\)$",
    out
  )))
})

## Expected or assured entrants

test_that("a plan with an assurance declares which count it fields", {
  plan <- n_panel(n_mean(var = 1, moe = 0.196), retention = 1,
                  resp_rate = 0.8, design = "rotating", rotation = "2",
                  assurance = 0.95, start = "immediate")
  expect_error(
    design_schedule(plan, 4, "continuing", "entrant_register"),
    "68 per period against 62.5 expected; declare entrants"
  )
  assured <- design_schedule(plan, 4, "continuing", "entrant_register",
                             entrants = "assured")
  expected <- design_schedule(plan, 4, "continuing", "entrant_register",
                              entrants = "expected")
  expect_identical(assured$entrant_basis, "assured")
  expect_identical(expected$entrant_basis, "expected")
  expect_equal(assured$n_entrants, 68)
  expect_equal(assured$selections$operational_issue[[2L]], 68)
  expect_equal(expected$selections$operational_issue[[2L]], 63)
  expect_equal(assured$n_in_sample, 136)
  expect_true(any(grepl("^entrants: 68 per period \\(assured count\\)",
                        capture.output(print(summary(assured))))))
  # the assured count reaches 100 respondents with probability 0.95
  expect_gte(1 - stats::pbinom(99, 2 * 68, 0.8), 0.95)
  expect_lt(1 - stats::pbinom(99, 2 * 63, 0.8), 0.95)
  no_assurance <- n_panel(n_mean(var = 1, moe = 0.196), retention = 1,
                          design = "rotating", rotation = "2",
                          start = "immediate")
  expect_error(
    design_schedule(no_assurance, 4, "continuing", "entrant_register",
                    entrants = "assured"),
    "needs a plan with an assurance"
  )
  expect_identical(
    design_schedule(no_assurance, 4, "continuing",
                    "entrant_register")$entrant_basis,
    "expected"
  )
})
