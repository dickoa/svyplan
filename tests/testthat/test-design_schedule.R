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
    wave = tab$life_stage,
    n_issued = e * tab$relative_take,
    n_resp = e * tab$relative_take * q[tab$life_stage]
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
    design_overlap("4"),
    horizon = 6,
    horizon_policy = "continuing",
    refreshment = "entrant_register",
    frame_vintage = c(startup = "2026Q1", intake_2 = "2026Q2"),
    rounding = "ceiling"
  )

  expect_s3_class(x, "svyplan_schedule")
  expect_identical(x$schema_version, 1L)
  expect_identical(x$launch_policy, "immediate")
  expect_identical(x$horizon_policy, "continuing")
  expect_identical(x$steady_state_from, 4L)
  expect_identical(x$components$cohort,
                   c("startup", paste0("intake_", 2:6)))
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
    x$schedule$cohort == "startup" & x$schedule$panel == 1L &
      x$schedule$active,
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
        plan, design_overlap("4"), 6, policy,
        refreshment = "entrant_register", rounding = "ceiling"
      )
      got <- schedule_active_counts(x)
      got <- got[got$period <= 5L, ]
      want <- plan$launch_waves[plan$launch_waves$period <= 5L, ]
      rownames(got) <- NULL
      rownames(want) <- NULL
      expect_equal(got, want, tolerance = 1e-12,
                   info = paste(start, policy))

      later <- schedule_active_counts(x)
      later <- later[later$period == 6L, -1L, drop = FALSE]
      settled <- plan$launch_waves[
        plan$launch_waves$period == 5L, -1L, drop = FALSE
      ]
      rownames(later) <- rownames(settled) <- NULL
      expect_equal(later, settled, tolerance = 1e-12,
                   info = paste(start, policy, "post-launch"))
    }
  }

  plan <- schedule_panel_plan("immediate")
  closed <- design_schedule(
    plan, design_overlap("4"), 6, "close_intake",
    refreshment = "entrant_register", rounding = "ceiling"
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
    schedule_panel_plan("immediate"), design_overlap("4"), 8, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
  )
  for (lag in seq_along(expected)) {
    origins <- seq_len(immediate$horizon - lag)
    got <- vapply(origins, schedule_realized_overlap, numeric(1L),
                  x = immediate, lag = lag)
    expect_equal(got, rep(expected[[lag]], length(got)), tolerance = 1e-12)
  }

  gradual <- design_schedule(
    schedule_panel_plan("gradual"), design_overlap("4"), 8, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
  )
  expect_equal(schedule_realized_overlap(1L, gradual, 1L), 1)
  for (lag in seq_along(expected)) {
    origins <- 4:(gradual$horizon - lag)
    if (length(origins) > 0L && all(origins <= gradual$horizon - lag)) {
      got <- vapply(origins, schedule_realized_overlap, numeric(1L),
                    x = gradual, lag = lag)
      expect_equal(got, rep(expected[[lag]], length(got)), tolerance = 1e-12)
    }
  }
})

test_that("horizon policies separate issue, commitments and tail composition", {
  plan <- schedule_panel_plan()
  make <- function(policy) design_schedule(
    plan, design_overlap("4"), 6, policy,
    refreshment = "entrant_register", rounding = "ceiling"
  )
  continuing <- make("continuing")
  truncated <- make("truncate_lives")
  closed <- make("close_intake")

  expect_equal(sum(continuing$issue$operational_issue),
               sum(truncated$issue$operational_issue))
  expect_gt(sum(continuing$issue$operational_issue),
            sum(closed$issue$operational_issue))
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
    schedule_panel_plan("gradual"), design_overlap("4"), 6, "continuing",
    refreshment = "whole_vintage", rounding = "ceiling"
  )

  expect_identical(x$components$panels[[1L]], 1L)
  expect_true(all(x$components$frame_role[-1L] == "whole_vintage"))
  expect_identical(x$steady_state_from, 4L)
  expect_true(is.na(x$panel_parameters$block_width))
})

test_that("inputs that do not define the first schema are refused", {
  rotating <- schedule_panel_plan()
  fixed <- n_panel(
    n_prop(p = 0.5, moe = 0.03), retention = rep(0.9, 3),
    resp_rate = 0.75
  )
  launchless <- n_panel(
    n_prop(p = 0.5, moe = 0.03), retention = rep(0.9, 3),
    resp_rate = 0.75, design = "rotating"
  )

  expect_error(
    design_schedule(fixed, design_overlap("4"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "rotating"
  )
  expect_error(
    design_schedule(launchless, design_overlap("4"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "launch"
  )
  expect_error(
    design_schedule(rotating, design_overlap("3"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "life length"
  )
  expect_error(
    design_schedule(rotating, design_overlap("1-0-1"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "unbroken"
  )
  expect_error(
    design_schedule(rotating, design_overlap(c(1, 0.5, 0.5, 0.5)), 6,
                    "continuing", refreshment = "entrant_register",
                    rounding = "ceiling"),
    "equal-take"
  )
  expect_error(
    design_schedule(rotating, design_overlap("4"), 0, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "horizon"
  )
  expect_error(
    design_schedule(rotating, design_overlap("4"), 3, "close_intake",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "shorter"
  )
  expect_error(
    design_schedule(rotating, design_overlap("4"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "nearest"),
    "rounding"
  )
  expect_error(
    design_schedule(rotating, design_overlap("4"), 6, "later",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "horizon_policy"
  )
  expect_error(
    design_schedule(rotating, design_overlap("4"), 6, "continuing",
                    refreshment = "other", rounding = "ceiling"),
    "refreshment"
  )
  expect_error(
    design_schedule(rotating, design_overlap("4"), 6, "continuing",
                    refreshment = "entrant_register", rounding = "ceiling",
                    frame_vintage = c(nope = "x")),
    "frame_vintage"
  )
})

test_that("horizon policy and rounding are explicit", {
  plan <- schedule_panel_plan()
  expect_error(
    design_schedule(plan, design_overlap("4"),
                    horizon_policy = "continuing",
                    refreshment = "entrant_register", rounding = "ceiling"),
    "horizon"
  )
  expect_error(
    design_schedule(plan, design_overlap("4"), 6,
                    refreshment = "entrant_register", rounding = "ceiling"),
    "horizon_policy"
  )
  expect_error(
    design_schedule(plan, design_overlap("4"), 6, "continuing",
                    refreshment = "entrant_register"),
    "rounding"
  )
  expect_error(
    design_schedule(plan, design_overlap("4"), 6, "continuing",
                    rounding = "ceiling"),
    "refreshment"
  )
})

test_that("the versioned schedule schema is validated at construction", {
  x <- design_schedule(
    schedule_panel_plan(), design_overlap("4"), 6, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
  )

  expect_invisible(svyplan:::.validate_schedule_schema(unclass(x)))

  broken <- unclass(x)
  broken$schedule$life_stage[!broken$schedule$active][[1L]] <- 1L
  expect_error(
    svyplan:::.validate_schedule_schema(broken),
    "dense schedule grid"
  )

  broken <- unclass(x)
  broken$tail_commitments$active <- TRUE
  expect_error(
    svyplan:::.validate_schedule_schema(broken),
    "fixed schema-version-1 columns"
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
      plan, design_overlap("4"), 6, "continuing",
      refreshment = "entrant_register", rounding = "ceiling"
    ),
    NA
  )
  expect_type(x$components$panel_issue, "double")
  expect_equal(x$components$panel_issue[[1L]], large)
  expect_equal(x$components$operational_issue[[1L]], 4 * large)
  expect_true(all(x$components$operational_issue ==
                    floor(x$components$operational_issue)))
  expect_invisible(svyplan:::.validate_schedule_schema(unclass(x)))
})

test_that("schedule replacement cannot create cross-field contradictions", {
  x <- design_schedule(
    schedule_panel_plan(), design_overlap("4"), 6, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
  )

  expect_error({x$rounding <- list(rule = "floor")},
               "cannot be modified in place")
  expect_error({x[["n_in_sample"]] <- 1}, "cannot be modified in place")
  expect_error({x["schedule"] <- list(NULL)}, "cannot be modified in place")
  expect_error({x$schedule$panel[[1L]] <- 99L},
               "cannot be modified in place")
})

test_that("schedule methods lead with the issue profile", {
  x <- design_schedule(
    schedule_panel_plan(), design_overlap("4"), 6, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
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
    schedule_panel_plan(), design_overlap("4"), 6, "continuing",
    refreshment = "entrant_register", rounding = "ceiling"
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))

  saveRDS(x, path)
  restored <- readRDS(path)

  expect_identical(restored, x)
  expect_s3_class(restored, "svyplan_schedule")
  expect_identical(restored$schema_version, 1L)
  expect_invisible(svyplan:::.validate_schedule_schema(unclass(restored)))
  expect_identical(as.data.frame(restored), restored$issue)
  expect_match(format(restored), "svyplan_schedule.*immediate.*continuing")
})
