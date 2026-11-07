test_that("plot.svyplan_strata runs without error", {
  set.seed(1)
  x <- rlnorm(500, 6, 1)
  sb <- strata_bound(x, n_strata = 3, n = 100, method = "cumrootf")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(sb))
})

test_that("plot.svyplan_strata returns invisible(x)", {
  set.seed(1)
  x <- rlnorm(500, 6, 1)
  sb <- strata_bound(x, n_strata = 3, n = 100, method = "cumrootf")
  pdf(tempfile())
  on.exit(dev.off())
  out <- plot(sb)
  expect_identical(out, sb)
})

test_that("plot.svyplan_power works (solved for n)", {
  pw <- power_prop(p1 = 0.30, p2 = 0.40, power = 0.80)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(pw))
})

test_that("plot.svyplan_power works (solved for power)", {
  pw <- power_prop(p1 = 0.30, p2 = 0.40, n = 500, power = NULL)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(pw))
})

test_that("plot.svyplan_power works (solved for mde)", {
  pw <- power_prop(p1 = 0.30, n = 500)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(pw))
})

test_that("plot.svyplan_strata accepts user overrides", {
  set.seed(1)
  x <- rlnorm(500, 6, 1)
  sb <- strata_bound(x, n_strata = 3, n = 100, method = "cumrootf")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(sb, main = "Custom", col = "steelblue", ylab = "f_h"))
})

test_that("plot.svyplan_power accepts user overrides", {
  pw <- power_prop(p1 = 0.30, p2 = 0.40, power = 0.80)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(pw, main = "Custom", col = "red", lwd = 2))
})

test_that("plot.svyplan_power returns invisible(x)", {
  pw <- power_mean(effect = 5, var = 100)
  pdf(tempfile())
  on.exit(dev.off())
  out <- plot(pw)
  expect_identical(out, pw)
})

test_that("plot.svyplan_power errors for vector n", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  expect_error(plot(pw), "does not support power objects with unequal-group n")
})

test_that("plot.svyplan_n draws the budget frontier", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 700)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(fit, npoints = 6L))
  expect_identical(plot(fit, npoints = 6L), fit)
})

test_that("the default frontier grid contains no infeasible budgets", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 700)
  grid <- .budget_grid(fit, 8L)
  expect_length(grid, 8L)
  expect_true(all(is.finite(grid)))

  # every point solves, so plotting emits no infeasibility warning
  fr <- predict(fit, data.frame(budget = grid))
  expect_true(all(fr$.feasible))
})

test_that("plot.svyplan_n rejects results with no frontier", {
  expect_error(plot(n_prop(p = 0.3, moe = 0.05)),
               "fixed-budget joint allocations")

  z <- .bethel_fixture()
  min_cost <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  expect_error(plot(min_cost), "fixed-budget joint allocations")
})

test_that("plot.svyplan_n accepts an explicit budget grid", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 700)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(fit, newdata = data.frame(budget = c(600, 700, 800))))
  expect_error(
    plot(fit, newdata = data.frame(budget = 700)),
    "too few feasible budgets"
  )
})

test_that("the frontier title names one component and counts several", {
  one <- data.frame(component = "income@.overall")
  many <- data.frame(component = c("income@.overall", "vaccination@region"))
  expect_identical(.frontier_title(one), "Budget frontier for income@.overall")
  expect_identical(.frontier_title(many),
                   "Budget frontier (2 objective components)")
  expect_identical(.frontier_title(NULL), "Budget frontier")
})

## P1. The rotation chart and its cell table

.cell_totals <- function(d) {
  vapply(seq_len(attr(d, "horizon")), function(t) {
    sum(d$take[d$period == t])
  }, numeric(1L))
}

test_that("each object draws its own chart and returns invisible(x)", {
  rot <- design_rotation("1-1-0-0-1-1")
  ov <- design_overlap(rot)
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(rot))
  expect_silent(plot(ov))
  expect_identical(plot(rot), rot)
  expect_identical(plot(ov), ov)
})

test_that("the overlap chart points a chart argument to the rotation", {
  ov <- design_overlap("6")
  pdf(tempfile())
  on.exit(dev.off())
  for (arg in c("type", "start", "horizon", "panel_plan", "label", "fill")) {
    args <- setNames(list(ov, "x"), c("x", arg))
    expect_error(do.call(plot, args), "drawn from the rotation", fixed = TRUE)
  }
  expect_silent(plot(ov, col = "steelblue", main = "CPS"))
})

test_that("the chart draws patterns of every shape", {
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(design_rotation("4-8-4")))
  expect_silent(plot(design_rotation("5")))
  expect_silent(plot(design_rotation(c(1, 1, 0.5, 0.5))))
  expect_silent(plot(design_rotation("2"), horizon = 3L))
  expect_silent(plot(design_rotation("3-1-2"), rows = "rotation_group",
                     fill = "panel", horizon = 14L))
})

test_that("horizon is validated", {
  x <- design_rotation("5")
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(x, horizon = 0), "whole number >= 1")
  expect_error(plot(x, horizon = 2.5), "whole number >= 1")
  expect_error(plot(x, horizon = c(4, 5)), "whole number >= 1")
  expect_error(as.data.frame(x, horizon = 0), "whole number >= 1")
  expect_error(plot(x, rows = "nope"), "should be one of")
  expect_error(plot(x, bogus = 1), "unused argument")
})

test_that("a rotating panel labels the chart in units", {
  target <- n_prop(p = 0.5, moe = 0.031)
  rot <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                 resp_rate = 0.728, design = "rotating")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(design_rotation("5"), panel_plan = rot))
  d <- as.data.frame(design_rotation("5"), panel_plan = rot)
  expect_equal(d$n, d$take * ceiling(rot$n_entrants))
  expect_false("n" %in% names(as.data.frame(design_rotation("5"))))
})

test_that("the chart's entrant count is the one print reports", {
  # the same object must not report two recruitments, so the chart issues
  # whole units the way print.svyplan_panel() does
  target <- n_prop(p = 0.5, moe = 0.031)
  rot <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                 resp_rate = 0.728, design = "rotating")
  take <- .chart_take(rot, rep(1, 5))
  expect_equal(take, ceiling(rot$n_entrants))
  out <- capture.output(print(rot))
  expect_true(any(grepl(sprintf("^entrants: %d per period", take), out)))
  # and the full-sample total the chart draws is the one summary() names
  expect_true(any(grepl(sprintf("^in sample: %d ", take * rot$n_cohorts),
                        capture.output(print(summary(rot))))))
})

test_that("a panel that does not describe the rotation is refused", {
  target <- n_prop(p = 0.5, moe = 0.031)
  fixed <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                   resp_rate = 0.728)
  rot <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                 resp_rate = 0.728, design = "rotating")
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(design_rotation("5"), panel_plan = fixed), "fixed panel")
  expect_error(plot(design_rotation("4"), panel_plan = rot),
               "follows the rotation 5 in and the chart draws 4 in")
  expect_error(plot(design_rotation("5"), panel_plan = target),
               "must be an n_panel\\(\\) or prec_panel\\(\\) result")
})

test_that("the cell table has one row per interview, in a fixed shape", {
  d <- as.data.frame(design_rotation("2-2-2"), start = "immediate",
                     horizon = 8)
  expect_named(d, c("period", "panel", "entry", "startup", "rotation_group",
                    "wave", "visit", "life_stage", "take"))
  for (col in c("period", "panel", "entry", "rotation_group", "wave",
                "visit", "life_stage")) {
    expect_type(d[[col]], "integer")
  }
  expect_type(d$startup, "logical")
  expect_type(d$take, "double")
  expect_identical(order(d$panel, d$period), seq_len(nrow(d)))
  expect_false(anyDuplicated(d[c("panel", "period")]) > 0L)
  expect_true(all(d$period >= 1L & d$period <= 8L))
  expect_identical(attr(d, "start"), "immediate")
  expect_identical(attr(d, "horizon"), 8L)
  expect_identical(attr(d, "full_sample_from"), 1L)
  expect_s3_class(attr(d, "rotation"), "svyplan_rotation")
  g <- as.data.frame(design_rotation("2-2-2"))
  expect_identical(attr(g, "full_sample_from"), 6L)
  expect_identical(attr(g, "horizon"), 10L)
  expect_false(any(g$startup))
})

test_that("waves count the pattern and visits count interviews", {
  d <- as.data.frame(design_rotation("6"), start = "immediate",
                     horizon = 12)
  first <- d[d$period == 1L, ]
  # Statistics Canada, Figure 2.1: January holds waves 6 down to 1, every
  # one of them a first interview
  expect_identical(first$panel, 1:6)
  expect_identical(first$wave, 6:1)
  expect_identical(first$visit, rep(1L, 6L))
  expect_identical(d$wave[!d$startup], d$visit[!d$startup])
  expect_true(any(d$wave[d$startup] != d$visit[d$startup]))
})

test_that("an immediate 2-2-2 start reproduces the Rwanda LFS annex", {
  # NISR (2024), Annex 1, rows 1 to 12 over 2024 Q1 to 2026 Q4. A to D are
  # the first to fourth appearance in the 2-2-2 pattern
  d <- as.data.frame(design_rotation("2-2-2"), start = "immediate",
                     horizon = 12)
  d <- d[d$panel <= 12L, ]
  got <- sort(paste0(d$panel, LETTERS[d$wave], "@", d$period))
  startup <- c("1D@1", "2C@1", "2D@2", "3C@2", "3D@3", "4C@3", "4D@4",
              "5B@1", "5C@4", "5D@5")
  later <- unlist(lapply(6:12, function(k) {
    paste0(k, c("A", "B", "C", "D"), "@", k - c(5L, 4L, 1L, 0L))
  }))
  expect_identical(got, sort(c(startup, later)))
})

test_that("rotation groups are the slots the published charts show", {
  g <- as.data.frame(design_rotation("2-2-2"), horizon = 30)
  by_panel <- tapply(g$rotation_group, g$panel, `[`, 1L)
  # NISR (2024), Figure 1: rotation group A holds the panels entering at
  # quarters 1, 3, 9 and 11
  expect_identical(unname(which(by_panel == 1L))[1:4], c(1L, 3L, 9L, 11L))
  full <- g[g$period >= 6L, ]
  expect_false(anyDuplicated(full[c("period", "rotation_group")]) > 0L)
  expect_identical(max(g$rotation_group), 4L)
  for (spec in c("6", "4-8-4")) {
    d <- as.data.frame(design_rotation(spec), horizon = 40)
    life <- length(design_rotation(spec))
    held <- d[d$period >= life, ]
    expect_identical(max(d$rotation_group),
                     as.integer(sum(design_rotation(spec) > 0)))
    expect_false(anyDuplicated(held[c("period", "rotation_group")]) > 0L)
  }
  # 3-1-2 cannot fill every group at every period: six groups for five
  # panels in sample
  d <- as.data.frame(design_rotation("3-1-2"), horizon = 40)
  expect_identical(max(d$rotation_group), 6L)
  expect_true(all(table(d$period[d$period >= 6L]) == 5L))
})

test_that("widening the window never moves a panel to another group", {
  for (start in c("gradual", "immediate")) {
    a <- as.data.frame(design_rotation("2-2-2"), start = start, horizon = 8)
    b <- as.data.frame(design_rotation("2-2-2"), start = start, horizon = 30)
    ga <- tapply(a$rotation_group, a$panel, `[`, 1L)
    gb <- tapply(b$rotation_group, b$panel, `[`, 1L)
    expect_identical(ga, gb[names(ga)])
  }
})

test_that("the cell table holds the interviews the schedule fields", {
  for (spec in c("6", "2-2-2", "4-8-4")) {
    for (start in c("gradual", "immediate")) {
      rot <- design_rotation(spec)
      plan <- n_panel(n_prop(p = 0.5, moe = 0.05),
                      retention = rep(1, sum(rot > 0) - 1L),
                      design = "rotating", start = start, rotation = rot)
      s <- design_schedule(plan, horizon = 20, horizon_policy = "continuing",
                           refreshment = "entrant_register",
                           rounding = "ceiling")
      a <- s$schedule[s$schedule$active, ]
      d <- as.data.frame(rot, start = start, horizon = 20)
      expect_identical(sort(paste(d$period, d$visit)),
                       sort(paste(a$period, a$visit)))
    }
  }
})

test_that("the total climbs to the full sample and holds", {
  # the overlap describes the full sample, so the chart has to reach it:
  # from the period that spans the life, the total is design_overlap()'s
  # own in_sample and stays there
  for (spec in c("4-8-4", "5", "1-1-0-0-1-1")) {
    rot <- design_rotation(spec)
    life <- length(rot)
    d <- as.data.frame(rot, horizon = life + 4L)
    totals <- .cell_totals(d)
    expect_equal(totals[life:(life + 4L)],
                 rep(attr(rot, "in_sample"), 5L), tolerance = 1e-12)
    expect_lt(totals[1L], attr(rot, "in_sample"))
    expect_false(is.unsorted(totals))
    expect_identical(attr(d, "full_sample_from"), life)
  }
  # a varying take counts the units, not the periods
  d <- as.data.frame(design_rotation(c(1, 1, 0.5, 0.5)), horizon = 6)
  expect_equal(.cell_totals(d), c(1, 2, 2.5, 3, 3, 3), tolerance = 1e-12)
})

test_that("the chart draws an immediate start, full from the first period", {
  x <- design_rotation("6")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(x, start = "immediate"))
  expect_identical(plot(x, start = "immediate"), x)
  expect_error(plot(x, start = "later"), "should be one of")

  # Lynn Figure 5: start-up cohorts remain for one through six interviews, so
  # the sample is full from period 1 and never climbs
  co <- .rotation_cohorts(as.double(x), "immediate", 10L)
  expect_length(co, length(x) + 10L - 1L)
  expect_identical(
    vapply(co[seq_along(x)], function(z) length(z$w), integer(1L)),
    seq_along(x)
  )
  expect_identical(
    vapply(co[seq_along(x)], `[[`, integer(1L), "entry"),
    rep(1L, length(x))
  )
  expect_identical(co[[length(x) + 1L]]$entry, 2L)
  d <- as.data.frame(x, start = "immediate", horizon = 10)
  expect_equal(.cell_totals(d), rep(attr(x, "in_sample"), 10L),
               tolerance = 1e-12)
  grad <- .cell_totals(as.data.frame(x, horizon = 10))
  expect_lt(grad[[1L]], attr(x, "in_sample"))
  expect_equal(grad[length(x):10L],
               rep(attr(x, "in_sample"), 10L - length(x) + 1L),
               tolerance = 1e-12)
})

test_that("an immediate start of a gapped life is drawn as n_panel plans it", {
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(design_rotation("1-1-0-0-1-1"), start = "immediate"))
  expect_silent(plot(design_rotation("4-8-4"), start = "immediate"))
  gap <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(0.9, 0.85, 0.9),
                 resp_rate = 0.8, design = "rotating", start = "immediate",
                 rotation = "2-2-2")
  expect_silent(plot(gap$rotation, panel_plan = gap))
  n <- length(gap$rotation) + 1L
  d <- as.data.frame(gap$rotation, panel_plan = gap, horizon = n)
  expect_equal(.cell_totals(d) * gap$n_entrants, gap$periods$n_in_sample,
               tolerance = 1e-12)
})

test_that("a panel may only scale the rotation it was planned on", {
  target <- n_prop(p = 0.5, moe = 0.031)
  six <- n_panel(target, retention = rep(0.9, 5), resp_rate = 0.8,
                 design = "rotating")
  four <- n_panel(target, retention = rep(0.9, 3), resp_rate = 0.8,
                  design = "rotating")
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(design_rotation("1-1-0-0-1-1"), panel_plan = six),
               "follows the rotation 6 in and the chart draws 2 in, 2 out, 2 in")
  expect_error(plot(design_rotation(c(1, 1, 0.5, 0.5)), panel_plan = four),
               "must describe the same design")
  expect_error(plot(design_rotation("4"), panel_plan = six),
               "must describe the same design")
  expect_silent(plot(design_rotation("6"), panel_plan = six))
  expect_silent(plot(design_rotation("4"), panel_plan = four))
})

test_that("the chart draws the start-up its panel was planned with", {
  target <- n_prop(p = 0.5, moe = 0.031)
  ret <- c(0.9, 0.9)
  imm <- n_panel(target, retention = ret, resp_rate = 0.8,
                 design = "rotating", start = "immediate")
  grad <- n_panel(target, retention = ret, resp_rate = 0.8,
                  design = "rotating", start = "gradual")
  plain <- n_panel(target, retention = ret, resp_rate = 0.8,
                   design = "rotating")
  # inherited when not stated
  expect_identical(.check_chart_start(NULL, imm), "immediate")
  expect_identical(.check_chart_start(NULL, grad), "gradual")
  # and "gradual" where the panel carries no start-up, or there is no panel
  expect_identical(.check_chart_start(NULL, plain), "gradual")
  expect_identical(.check_chart_start(NULL, NULL), "gradual")
  # stating it wins, which is how one plan's two starts are compared
  expect_identical(.check_chart_start("gradual", imm), "gradual")
  expect_identical(.check_chart_start("immediate", grad), "immediate")
  expect_identical(attr(as.data.frame(imm$rotation, panel_plan = imm), "start"),
                   "immediate")
  expect_identical(
    attr(as.data.frame(imm$rotation, panel_plan = imm, start = "gradual"),
         "start"),
    "gradual"
  )

  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(imm$rotation, panel_plan = imm))
  expect_silent(plot(imm$rotation, panel_plan = imm, start = "gradual"))
})

test_that("labels, fills and period labels take presets, columns or functions", {
  x <- design_rotation("2-2-2")
  d <- as.data.frame(x, start = "immediate", horizon = 12)
  expect_identical(.chart_labels(d, "wave"), paste0("W", d$wave))
  expect_null(.chart_labels(d, "none"))
  expect_identical(.chart_labels(d, "visit"), as.character(d$visit))
  expect_identical(.chart_labels(d, function(z) paste0(z$panel, "x")),
                   paste0(d$panel, "x"))
  expect_error(.chart_labels(d, "nope"), "a column of the cell table")
  expect_error(.chart_labels(d, function(z) 1:3), "returned 3 values")
  expect_error(.chart_labels(d, c("wave", "visit")),
               "a column of the cell table")

  expect_identical(unique(.chart_fill(d, "none", NULL)), "grey75")
  expect_error(.chart_fill(d, "none", c("red", "blue")), "one color")
  f <- .chart_fill(d, "wave", c("red", "blue", "green", "black"))
  expect_identical(as.vector(f), c("red", "blue", "green", "black")[d$wave])
  expect_identical(attr(f, "legend"),
                   list(values = 1:4, col = c("red", "blue", "green", "black")))
  expect_true(.chart_needs_legend(f, "wave", "visit"))
  expect_false(.chart_needs_legend(f, "wave", "wave"))
  expect_false(.chart_needs_legend(.chart_fill(d, "none", NULL), "none",
                                   "visit"))
  expect_length(unique(.chart_fill(d, "rotation_group", NULL)), 4L)
  expect_error(.chart_fill(d, "wave", "red"), "1 colors for 4 fill values")
  expect_error(.chart_fill(d, "nope", NULL), "a column of the cell table")
  expect_error(.chart_fill(d, function(z) 1, NULL), "returned 1 values")
  expect_identical(
    as.vector(.chart_fill(d, function(z) z$wave > 2, c("a", "b"))),
    c("a", "b")[(d$wave > 2) + 1L]
  )

  expect_identical(.chart_period_labels(NULL, 3L), 1:3)
  expect_identical(.chart_period_labels(c("Q1", "Q2"), 2L), c("Q1", "Q2"))
  expect_error(.chart_period_labels(month.abb, 10L), "12 labels for 10")

  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(x, start = "immediate", horizon = 12,
                     label = function(z) paste0(z$panel, LETTERS[z$wave]),
                     fill = "wave", period_labels = paste0("Q", 1:12)))
  expect_silent(plot(x, rows = "rotation_group", label = "panel",
                     marker = FALSE))
  expect_error(plot(x, marker = NA), "TRUE or FALSE")
})

## P2. Curves and grids read the fit's rates

test_that("a power curve passes through the fitted point", {
  e <- list(resp_rate = 0.8, eligibility_rate = 0.5)
  fits <- list(
    do.call(power_prop, c(list(p1 = 0.3, p2 = 0.4, N = 5000), e)),
    do.call(power_mean, c(list(var = 100, effect = 2, N = 5000), e)),
    do.call(power_did, c(list(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                              outcome = "prop", N = 6000, overlap = 0.5,
                              overlap_cor = 0.4), e)),
    do.call(power_prop, c(list(p1 = 0.3, p2 = 0.4, N = 5000,
                               eligible_per_unit = 2.6), e)),
    do.call(power_mean, c(list(var = 100, effect = 2, N = 5000,
                               eligible_per_unit = 2.6), e))
  )
  for (x in fits) {
    expect_equal(svyplan:::.power_curve_grid(x, x$n), x$power,
                 info = x$type)
  }
})

test_that("a budget grid starts at the cheapest whole design in scope", {
  frame <- data.frame(stratum = c("North urban", "North rural", "South urban"),
                      region = c("North", "North", "South"),
                      N = c(1000, 1800, 1200))
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("coverage", "income"), each = 3),
    p = c(0.5, 0.4, 0.6, NA, NA, NA),
    mean = c(NA, NA, NA, 50, 55, 60),
    sd = c(NA, NA, NA, 10, 12, 15)
  )
  targets <- data.frame(name = c("coverage", "income"),
                        domain = c(".overall", "region"),
                        level = c(NA, "North"), cv = c(0.06, NA),
                        moe = c(NA, 3))
  fit <- n_alloc(frame, measures = measures, targets = targets,
                 objective = "income", budget = 1000,
                 eligibility_rate = 0.7)
  floor_fit <- n_alloc(frame, measures = measures, targets = targets,
                       eligibility_rate = 0.7)
  expect_equal(svyplan:::.budget_grid(fit, 5)[1L],
               floor_fit$operational$cost)
})

test_that("the chart and the cell table take the schedule's argument names", {
  rot <- design_rotation("2-2-2")
  expect_error(as.data.frame(rot, n_period = 8), "n_period")
  d <- as.data.frame(rot, start = "immediate", horizon = 8)
  expect_identical(attr(d, "horizon"), 8L)
  expect_false("launch" %in% names(d))
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(rot, n_period = 8), "n_period")
})
