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
    period = tab$wave,
    wave = tab$visit,
    n_issued = e * tab$relative_take,
    n_resp = e * tab$relative_take * q[tab$visit]
  )
  out <- aggregate(
    cbind(n_issued, n_resp) ~ period + wave,
    data = out,
    FUN = sum
  )
  out[order(out$period, out$wave), ]
}

schedule_realized_overlap <- function(origin, x, lag) {
  active <- x$schedule[x$schedule$active, , drop = FALSE]
  key <- paste(active$cohort, active$panel, sep = "|")
  first <- unique(key[active$wave == origin])
  second <- unique(key[active$wave == origin + lag])
  length(intersect(first, second)) / length(first)
}

test_that("an immediate schedule has stable components and a dense grid", {
  plan <- schedule_panel_plan()

  x <- design_schedule(
    plan,
    horizon = 6,
    horizon_policy = "continuing",
    refreshment = "entrant_register",
    frame_vintage = c(startup = "2026Q1", intake_2 = "2026Q2"),
    rounding = "ceiling"
  )

  expect_s3_class(x, "svyplan_schedule")
  expect_identical(x$schema_version, 2L)
  expect_identical(x$launch_policy, "immediate")
  expect_identical(x$horizon_policy, "continuing")
  expect_identical(x$steady_state_from, 4L)
  expect_identical(x$components$cohort, c("startup", paste0("intake_", 2:6)))
  expect_identical(x$components$panels, c(4L, rep(1L, 5L)))
  expect_true(all(x$components$status == "planned"))
  expect_identical(x$components$frame_vintage[1:2], c("2026Q1", "2026Q2"))
  expect_true(all(is.na(x$components$frame_vintage[-(1:2)])))

  intake_n <- ceiling(plan$n_entrants)
  expect_equal(x$components$planned_issue[[1L]], plan$n_in_sample)
  expect_equal(x$components$operational_issue[[1L]], 4L * intake_n)
  expect_equal(x$components$panel_issue[[1L]], intake_n)

  expect_identical(nrow(x$schedule), 6L * (4L + 5L))
  expect_true(all(is.na(x$schedule$life_stage[!x$schedule$active])))
  expect_false(anyNA(x$schedule$life_stage[x$schedule$active]))
  expect_true(all(x$schedule$steady_state[x$schedule$wave >= 4L]))
  expect_false(any(x$schedule$steady_state[x$schedule$wave < 4L]))

  startup_1 <- x$schedule[
    x$schedule$cohort == "startup" & x$schedule$wave == 1L,
  ]
  expect_true(all(startup_1$active))
  expect_identical(startup_1$life_stage, rep(1L, 4L))
  expect_identical(startup_1$life_length, 4:1)
  startup_panel_1 <- x$schedule[
    x$schedule$cohort == "startup" & x$schedule$panel == 1L & x$schedule$active,
  ]
  expect_identical(startup_panel_1$life_stage, 1:4)
  intake_2 <- x$schedule[
    x$schedule$cohort == "intake_2" & x$schedule$active,
  ]
  expect_identical(intake_2$wave, 2:5)
  expect_identical(intake_2$life_stage, 1:4)

  expect_gt(nrow(x$tail_commitments), 0L)
  expect_gt(min(x$tail_commitments$wave), x$horizon)
  expect_lte(max(x$tail_commitments$wave), x$horizon + x$life_length - 1L)
  expect_identical(x$panel_parameters$k, 4L)
  expect_identical(x$panel_parameters$r_min, 1L)
  expect_identical(x$panel_parameters$block_width, 8L)
})

test_that("generated schedules reconcile with the launch table where applicable", {
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
      want <- plan$launch_waves[plan$launch_waves$period <= 5L, ]
      rownames(got) <- NULL
      rownames(want) <- NULL
      expect_equal(got, want, tolerance = 1e-12, info = paste(start, policy))

      later <- schedule_active_counts(x)
      later <- later[later$period == 6L, -1L, drop = FALSE]
      settled <- plan$launch_waves[
        plan$launch_waves$period == 5L,
        -1L,
        drop = FALSE
      ]
      rownames(later) <- rownames(settled) <- NULL
      expect_equal(
        later,
        settled,
        tolerance = 1e-12,
        info = paste(start, policy, "post-launch")
      )
    }
  }

  plan <- schedule_panel_plan("immediate")
  closed <- design_schedule(
    plan,
    6,
    "close_intake",
    refreshment = "entrant_register",
    rounding = "ceiling"
  )
  got <- schedule_active_counts(closed)
  got <- got[got$period <= 3L, ]
  want <- plan$launch_waves[plan$launch_waves$period <= 3L, ]
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
  closed <- make("close_intake")

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
      truncated$schedule$cohort == "intake_6" &
        truncated$schedule$wave == 6L
    ][[1L]],
    continuing$life_length
  )
  active_by_wave <- tapply(closed$schedule$active, closed$schedule$wave, sum)
  expect_true(all(diff(tail(active_by_wave, 4L)) < 0L))
  expect_true(is.na(closed$steady_state_from))
})

test_that("a gradual launch has one startup panel and reaches composition at L", {
  x <- design_schedule(
    schedule_panel_plan("gradual"),
    6,
    "continuing",
    refreshment = "whole_vintage",
    rounding = "ceiling"
  )

  expect_identical(x$components$panels[[1L]], 1L)
  expect_true(all(x$components$frame_role[-1L] == "whole_vintage"))
  expect_identical(x$steady_state_from, 4L)
  expect_true(is.na(x$panel_parameters$block_width))
})

test_that("inputs that do not define the first schema are refused", {
  rotating <- schedule_panel_plan()
  fixed <- n_panel(
    n_prop(p = 0.5, moe = 0.03),
    retention = rep(0.9, 3),
    resp_rate = 0.75
  )
  launchless <- n_panel(
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
      launchless,
      6,
      "continuing",
      refreshment = "entrant_register",
      rounding = "ceiling"
    ),
    "launch"
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
      "close_intake",
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

test_that("horizon policy and rounding are explicit", {
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
  expect_error(
    design_schedule(
      plan,
      6,
      "continuing",
      refreshment = "entrant_register"
    ),
    "rounding"
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

test_that("the versioned schedule schema is validated at construction", {
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
    "does not match the component manifest"
  )

  broken <- unclass(x)
  broken$tail_commitments$active <- TRUE
  expect_error(
    svyplan:::.validate_schedule_schema(broken),
    "fixed schema-version-2 columns"
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
    component_issue = function(z) {
      z$components$planned_issue[[1L]] <- 1
      z
    },
    issue_profile = function(z) {
      z$issue$cohort[[2L]] <- "startup"
      z
    },
    overlap = function(z) {
      z$overlap$overlap[[1L]] <- 0
      z
    },
    steady_state = function(z) {
      z$steady_state_from <- 1L
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
  expect_type(x$components$panel_issue, "double")
  expect_equal(x$components$panel_issue[[1L]], large)
  expect_equal(x$components$operational_issue[[1L]], 4 * large)
  expect_true(all(
    x$components$operational_issue == floor(x$components$operational_issue)
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
  expect_identical(restored$schema_version, 2L)
  expect_invisible(svyplan:::.validate_schedule_schema(unclass(restored)))
  expect_identical(as.data.frame(restored), restored$issue)
  expect_match(format(restored), "svyplan_schedule.*immediate.*continuing")
})

## S13. The printed block does not grow with the reporting horizon

test_that("print states the issue profile as runs, not one row per occasion", {
  mk <- function(h, policy = "continuing") {
    design_schedule(schedule_panel_plan(), h, policy,
                    refreshment = "entrant_register", rounding = "ceiling")
  }
  # A schedule's length is set by the reporting window, not by the design, so
  # the block must be the same size at every horizon.
  for (h in c(6L, 12L, 24L, 52L)) {
    out <- capture.output(print(mk(h)))
    expect_length(out, 6L)
    expect_lt(max(nchar(out)), 80L)
    expect_true(any(grepl(sprintf("occasions 2-%d\\)", h), out)))
  }
  # The occasions are gone from print and whole in summary.
  out <- capture.output(print(mk(24L)))
  expect_false(any(grepl("intake_7", out, fixed = TRUE)))
  expect_identical(nrow(summary(mk(24L))$issue), 24L)
})

test_that("each horizon policy states its own profile", {
  mk <- function(h, policy) {
    design_schedule(schedule_panel_plan(), h, policy,
                    refreshment = "entrant_register", rounding = "ceiling")
  }
  # Intake stopping is a fact about the schedule, not a gap in the report.
  closed <- capture.output(print(mk(10L, "close_intake")))
  expect_true(any(grepl("intake closed from 8", closed, fixed = TRUE)))
  expect_true(any(grepl("occasions 2-7", closed, fixed = TRUE)))
  expect_length(closed, 6L)
  for (policy in c("continuing", "truncate_lives")) {
    out <- capture.output(print(mk(12L, policy)))
    expect_false(any(grepl("intake closed", out, fixed = TRUE)))
    expect_length(out, 6L)
  }
})

test_that("the life line keeps the branch where no steady state is reached", {
  short <- design_schedule(schedule_panel_plan(), 2L,
                           "continuing", refreshment = "entrant_register",
                           rounding = "ceiling")
  expect_true(is.na(short$steady_state_from))
  out <- capture.output(print(short))
  expect_match(out[2L], "^life: 4 stages over 2 occasions, ")
  expect_match(out[2L], "steady state not reached in the window$")
})

test_that("summary carries every table print no longer shows", {
  x <- design_schedule(schedule_panel_plan(), 6L,
                       "continuing", refreshment = "entrant_register",
                       rounding = "ceiling")
  sm <- summary(x)
  expect_s3_class(sm, "summary.svyplan_schedule")
  expect_identical(sm$issue, x$issue)
  expect_identical(sm$components, x$components)
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

## Schema version 2 on a life with a gap

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
  expect_equal(as.vector(tapply(s$active, s$wave, sum)), rep(4L, 8L))
  expect_equal(as.vector(tapply(!is.na(s$life_stage), s$wave, sum)),
               rep(6L, 8L))
  expect_equal(x$n_held, 6 * x$n_entrants)
  st <- s[s$cohort == "startup", ]
  first <- vapply(1:6, function(j) min(st$wave[st$panel == j & st$active]),
                  numeric(1L))
  expect_identical(first, c(1, 1, 3, 2, 1, 1))
  expect_identical(st$visit[st$panel == 1L], c(1L, 2L, NA, NA, 3L, 4L, NA, NA))
  # the panels covering stages 3 and 4 are held, not absent, before
  # their first interview
  expect_identical(st$life_stage[st$panel == 3L][1:2], 1:2)
  expect_false(any(st$active[st$panel == 3L][1:2]))
  intake <- s[s$cohort == "intake_2", ]
  expect_identical(intake$visit, c(NA, 1L, 2L, NA, NA, 3L, 4L, NA))
  expect_identical(x$steady_state_from, 6L)
  expect_true(all(!is.na(x$tail_commitments$visit)))
  expect_identical(x$schedule$period, x$schedule$wave)
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
      for (policy in c("continuing", "truncate_lives", "close_intake")) {
        if (policy == "close_intake" && h < L) next
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
        last_entry <- if (policy == "close_intake") h - L + 1L else h
        o <- launch_oracle(w, start, last_entry + L)
        o$entry <- ifelse(substr(o$id, 1, 1) == "s", 1L,
                          as.integer(substring(o$id, 2)))
        o <- o[o$entry <= last_entry, ]
        if (policy == "truncate_lives") o <- o[o$period <= h, ]
        inside <- o[o$period <= h, ]
        got <- x$schedule[x$schedule$active, ]
        expect_identical(
          as.vector(table(factor(got$wave, levels = seq_len(h)))),
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
          ends <- stats::aggregate(wave ~ cohort + panel, held, max)
          last_int <- stats::aggregate(wave ~ cohort + panel,
                                       x$schedule[x$schedule$active, ], max)
          m <- merge(ends, last_int, by = c("cohort", "panel"))
          expect_identical(m$wave.x, m$wave.y, label = lab)
        }
      }
    }
  }
})

test_that("truncate_lives refuses a horizon before a startup panel's first interview", {
  expect_error(.gapped_schedule(horizon = 2, policy = "truncate_lives"),
               "which the last of them has at occasion 3")
  expect_s3_class(.gapped_schedule(horizon = 3, policy = "truncate_lives"),
                  "svyplan_schedule")
  expect_s3_class(.gapped_schedule(start = "gradual", horizon = 2,
                                   policy = "truncate_lives"),
                  "svyplan_schedule")
})

test_that("the version-2 validator refuses each broken rule", {
  x <- unclass(.gapped_schedule())
  held <- which(!x$schedule$active & !is.na(x$schedule$life_stage))[[1L]]
  live <- which(x$schedule$active)
  mutations <- list(
    version = function(z) {
      z$schema_version <- 1L
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
            6L * x$components$panel_issue[[1L]]), s, fixed = TRUE
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
    subperiods = mk(plan, subperiods = 7),
    both = mk(plan, take = 14, tenures = "open", subperiods = 3),
    gapped = mk(gapped, take = 14, subperiods = 3)
  )
  for (nm in names(cases)) {
    x <- cases[[nm]]
    e <- x$components$panel_issue[[1L]] * c(x$take, 1)[[1L + is.na(x$take)]] *
      x$subperiods
    expect_false(e == ceiling(x$n_entrants), info = nm)
    s <- capture.output(print(summary(x)))
    per <- if (x$subperiods == 1L) "occasion" else "period"
    expect_true(any(s == sprintf(
      "entrants: %d per %s, %d in sample across %d cohorts",
      e, per, e * x$n_cohorts, x$n_cohorts
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

## PSUs, tenures and sub-periods

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

test_that("take counts PSUs and the units drawn at each tenure", {
  plan <- .c_plan()
  x <- .c_schedule(plan, 8, take = 14)
  e <- plan$n_entrants
  psu <- x$components$panel_issue[[1L]]
  expect_identical(psu, ceiling(e / 14))
  expect_gte(psu * 14, e)
  expect_lt(psu * 14 - e, 14)
  expect_equal(x$components$planned_issue, e / 14 * x$components$panels)
  expect_identical(x$components$operational_issue, psu * x$components$panels)
  expect_identical(x$rounding$level, "psu")
  expect_identical(x$take, 14)
  opened <- tapply(x$schedule$active & x$schedule$visit %in% 1L,
                   x$schedule$wave, sum)
  expect_equal(x$issue$tenure_issue, as.vector(opened) * psu * 14)
  expect_error(.c_schedule(plan, 8, take = 2.5), "'take' must be NULL")
  expect_error(.c_schedule(plan, 8, tenures = 3), "not described yet")
})

test_that("open tenures keep the PSUs and rotate the units", {
  x <- .c_schedule(.c_plan(), 12, tenures = "open")
  expect_identical(x$components$cohort, "startup")
  expect_identical(x$components$panels, 5L)
  expect_identical(x$tenures, Inf)
  per_wave <- tapply(x$schedule$active, x$schedule$wave, sum)
  expect_equal(as.vector(per_wave), rep(5, 12))
  opening <- x$issue$tenure_issue / x$components$panel_issue[[1L]]
  expect_identical(opening, c(5, rep(1, 11)))
  grad <- .c_schedule(.c_plan(start = "gradual"), 12, tenures = "open")
  expect_identical(grad$components$cohort, c("startup", paste0("intake_", 2:5)))
  expect_identical(grad$issue$tenure_issue / grad$components$panel_issue[[1L]],
                   rep(1, 12))
  gap <- .c_schedule(.c_plan("2-2-2"), 10, tenures = "open")
  expect_identical(gap$components$panels, 6L)
  expect_equal(as.vector(tapply(gap$schedule$active, gap$schedule$wave, sum)),
               rep(4, 10))
  expect_equal(as.vector(tapply(!is.na(gap$schedule$life_stage),
                                gap$schedule$wave, sum)), rep(6, 10))
  expect_true(any(diff(gap$schedule$visit[gap$schedule$panel == 1L &
                                            gap$schedule$active]) < 0))
})

test_that("open tenures and one tenure field the same units", {
  specs <- c("3", "5", "2-2-2", "4-8-4", "2-1-3", "1-0-1", "1-1-0-1")
  for (spec in specs) for (start in c("gradual", "immediate")) {
    plan <- .c_plan(spec, start)
    L <- length(as.double(plan$rotation))
    for (h in unique(c(L, L + 3L, 2L * L + 1L))) {
      for (policy in c("continuing", "truncate_lives", "close_intake")) {
        lab <- paste(spec, start, h, policy)
        one <- .c_schedule(plan, h, policy)
        open <- .c_schedule(plan, h, policy, tenures = "open")
        a <- one$schedule[one$schedule$active, ]
        b <- open$schedule[open$schedule$active, ]
        expect_identical(tabulate(b$wave, h), tabulate(a$wave, h), label = lab)
        expect_identical(sort(b$visit), sort(a$visit), label = lab)
        expect_identical(open$issue$tenure_issue, one$issue$tenure_issue,
                         label = lab)
        expect_identical(nrow(open$tail_commitments),
                         nrow(one$tail_commitments), label = lab)
        expect_identical(open$steady_state_from, one$steady_state_from,
                         label = lab)
        expect_identical(open$issue$steady_state, one$issue$steady_state,
                         label = lab)
        expect_lte(sum(open$components$panels), L)
      }
    }
  }
})

test_that("sub-periods split each group into its own months", {
  for (start in c("immediate", "gradual")) for (tenures in list(1, "open")) {
    plan <- .c_plan("5", start)
    r <- .c_schedule(plan, 4, tenures = tenures)
    x <- .c_schedule(plan, 4, tenures = tenures, subperiods = 3)
    lab <- paste(start, tenures)
    s <- x$schedule
    expect_identical(x$horizon, 12L, label = lab)
    expect_identical(x$subperiods, 3L, label = lab)
    expect_identical(x$panel_parameters$k,
                     if (start == "immediate") 15L else 3L, label = lab)
    expect_identical(s$period, (s$wave - 1L) %/% 3L + 1L, label = lab)
    expect_identical(x$components$entry_wave,
                     (r$components$entry_wave - 1L) * 3L + 1L, label = lab)
    expect_identical(x$components$panels, r$components$panels * 3L,
                     label = lab)
    on <- s[s$active, ]
    expect_true(all((on$wave - 1L) %% 3L == (on$panel - 1L) %% 3L),
                label = lab)
    # collapsing the sub-period panels of a group gives back the grid
    s$group <- (s$panel - 1L) %/% 3L + 1L
    back <- stats::aggregate(active ~ period + cohort + group, s, any)
    want <- stats::aggregate(active ~ wave + cohort + panel, r$schedule, any)
    m <- merge(back, want, by.x = c("period", "cohort", "group"),
               by.y = c("wave", "cohort", "panel"))
    expect_identical(nrow(m), nrow(want), label = lab)
    expect_identical(m$active.x, m$active.y, label = lab)
    stage <- unique(s[!is.na(s$life_stage),
                      c("period", "cohort", "group", "life_stage")])
    expect_identical(anyDuplicated(stage[c("period", "cohort", "group")]), 0L,
                     label = lab)
    expect_identical(sort(on$visit),
                     sort(rep(r$schedule$visit[r$schedule$active], each = 3L)),
                     label = lab)
    if (!is.na(r$steady_state_from)) {
      expect_identical(x$steady_state_from,
                       (r$steady_state_from - 1L) * 3L + 1L, label = lab)
    }
  }
})

test_that("a quarterly design with monthly fieldwork and kept PSUs runs as in samplyr's probe", {
  x <- .c_schedule(.c_plan("5"), 4, take = 14, tenures = "open",
                   subperiods = 3)
  s <- x$schedule
  expect_identical(s$active, (s$wave - 1L) %% 3L == (s$panel - 1L) %% 3L)
  opening <- s[s$active & s$visit == 1L, ]
  expect_true(all((opening$wave - 1L) %% 3L == (opening$panel - 1L) %% 3L))
  expect_identical(x$rounding$level, "psu_and_subperiod")
  e <- x$n_entrants
  psu <- x$components$panel_issue[[1L]]
  expect_identical(psu, ceiling(e / (14 * 3)))
  expect_lt(psu * 14 * 3 - e, 14 * 3)
  out <- capture.output(print(x))
  expect_true(any(out == "take: 14 units per PSU per tenure"))
  expect_true(any(out == "tenures: open, PSUs kept through the horizon"))
  expect_true(any(grepl("^life: 5 periods over 4 periods of 3 waves", out)))
  expect_true(any(grepl("no PSU intake from period 2", out, fixed = TRUE)))
})

test_that("sub-periods split any rotation and keep its reporting grid", {
  specs <- c("3", "2-2-2", "4-8-4", "1-0-1")
  for (spec in specs) for (start in c("gradual", "immediate")) {
    for (s in c(1L, 3L, 13L)) for (tenures in list(1, "open")) {
      plan <- .c_plan(spec, start)
      L <- length(as.double(plan$rotation))
      lab <- paste(spec, start, s, tenures)
      r <- .c_schedule(plan, L + 2L, tenures = tenures)
      x <- .c_schedule(plan, L + 2L, tenures = tenures, subperiods = s)
      if (s == 1L) {
        expect_identical(x, r, label = lab)
        next
      }
      on <- x$schedule[x$schedule$active, ]
      want <- r$schedule[r$schedule$active, ]
      expect_identical(nrow(on), nrow(want) * s, label = lab)
      expect_equal(tabulate(on$period, L + 2L) / s,
                   tabulate(want$wave, L + 2L), label = lab)
      expect_identical(nrow(x$tail_commitments),
                       nrow(r$tail_commitments) * s, label = lab)
    }
  }
})

test_that("the version-2 validator refuses broken tenure and sub-period rules", {
  x <- unclass(.c_schedule(.c_plan("5"), 4, take = 14, tenures = "open",
                           subperiods = 3))
  one <- unclass(.c_schedule(.c_plan("5"), 8))
  mutations <- list(
    tenures = function(z) {
      z$tenures <- 1
      z
    },
    subperiods = function(z) {
      z$subperiods <- 2L
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
      z$schedule$period <- z$schedule$wave
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
  at <- which(restart$schedule$cohort == "startup" &
                restart$schedule$panel == 1L & restart$schedule$active)
  restart$schedule$visit[at[[3L]]] <- 1L
  expect_error(svyplan:::.validate_schedule_schema(restart),
               "does not count the interviews of its tenure")
})
