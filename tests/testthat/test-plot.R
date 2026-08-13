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

## P1. The rotation chart

test_that("plot.svyplan_overlap draws both views and returns invisible(x)", {
  x <- design_overlap("1-1-0-0-1-1")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(x))
  expect_silent(plot(x, type = "overlap"))
  expect_identical(plot(x), x)
  expect_identical(plot(x, type = "overlap"), x)
})

test_that("the chart draws schedules of every shape", {
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(design_overlap("4-8-4")))          # a gap
  expect_silent(plot(design_overlap("5")))              # no break
  expect_silent(plot(design_overlap(c(1, 1, 0.5, 0.5)))) # a varying take
  expect_silent(plot(design_overlap("2"), n_period = 3L))
})

test_that("n_period is validated", {
  x <- design_overlap("5")
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(x, n_period = 0), "whole number >= 1")
  expect_error(plot(x, n_period = 2.5), "whole number >= 1")
  expect_error(plot(x, n_period = c(4, 5)), "whole number >= 1")
  expect_error(plot(x, type = "nope"), "should be one of")
})

test_that("a rotating panel labels the chart in units", {
  target <- n_prop(p = 0.5, moe = 0.031)
  rot <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                 resp_rate = 0.728, design = "rotating")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(design_overlap("5"), panel = rot))
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
  expect_true(any(grepl(sprintf("^entrants: %d per occasion", take), out)))
  # and the steady-state total the chart draws is the one summary() names
  expect_true(any(grepl(sprintf("^in sample: %d ", take * rot$n_cohorts),
                        capture.output(print(summary(rot))))))
})

test_that("a panel that does not describe the schedule is refused", {
  target <- n_prop(p = 0.5, moe = 0.031)
  fixed <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                   resp_rate = 0.728)
  rot <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
                 resp_rate = 0.728, design = "rotating")
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(design_overlap("5"), panel = fixed), "fixed panel")
  expect_error(plot(design_overlap("4"), panel = rot),
               "5 live cohorts and the schedule has a life of 4")
  expect_error(plot(design_overlap("5"), panel = target),
               "must be an n_panel\\(\\) or prec_panel\\(\\) result")
})

test_that("the chart's total row climbs to the steady state and holds", {
  # the overlap describes a steady state, so the chart has to reach one:
  # from the period that spans the life, the total is design_overlap()'s
  # own n_occasion and stays there
  for (spec in c("4-8-4", "5", "1-1-0-0-1-1")) {
    x <- design_overlap(spec)
    life <- x$life
    totals <- .chart_totals(.chart_cohorts(x$schedule, "gradual", life + 4L),
                            life + 4L)
    expect_equal(totals[life:(life + 4L)], rep(x$n_occasion, 5L),
                 tolerance = 1e-12)
    expect_lt(totals[1L], x$n_occasion)
    expect_false(is.unsorted(totals))
  }
  # a varying take counts the units, not the occasions
  expect_equal(
    .chart_totals(.chart_cohorts(c(1, 1, 0.5, 0.5), "gradual", 6L), 6L),
    c(1, 2, 2.5, 3, 3, 3), tolerance = 1e-12
  )
})

test_that("the chart draws an immediate launch, full from the first period", {
  x <- design_overlap("6")
  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(x, start = "immediate"))
  expect_identical(plot(x, start = "immediate"), x)
  expect_error(plot(x, start = "later"), "should be one of")

  # Lynn Figure 5: launch cohorts remain for one through six interviews, so
  # the sample is full from period 1 and never climbs
  co <- .chart_cohorts(x$schedule, "immediate", 10L)
  expect_length(co, x$life + 10L - 1L)
  # Lynn orders the six launch samples by remaining life: Sample 1 appears
  # once, Sample 2 twice, ..., and Sample 6 for all six waves. Sample 7 is
  # the first ordinary entrant, at period 2.
  expect_identical(
    vapply(co[seq_len(x$life)], function(z) length(z$w), integer(1L)),
    seq_len(x$life)
  )
  expect_identical(
    vapply(co[seq_len(x$life)], `[[`, integer(1L), "entry"),
    rep(1L, x$life)
  )
  expect_identical(co[[x$life + 1L]]$entry, 2L)
  expect_identical(length(co[[x$life + 1L]]$w), x$life)
  totals <- .chart_totals(co, 10L)
  expect_equal(totals, rep(x$n_occasion, 10L), tolerance = 1e-12)

  # against the gradual launch, which reaches the same figure at the life
  grad <- .chart_totals(.chart_cohorts(x$schedule, "gradual", 10L), 10L)
  expect_lt(grad[[1L]], x$n_occasion)
  expect_equal(grad[x$life:10L], rep(x$n_occasion, 10L - x$life + 1L),
               tolerance = 1e-12)
})

test_that("an immediate launch of a gapped life is refused, as in n_panel", {
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(design_overlap("1-1-0-0-1-1"), start = "immediate"),
               "selected before their first interview")
  expect_error(plot(design_overlap("4-8-4"), start = "immediate"),
               "selected before their first interview")
  expect_silent(plot(design_overlap("1-1-0-0-1-1")))
})

test_that("a panel may only scale a schedule it represents", {
  # n_panel() models equal cohorts interviewed at every wave, so its
  # n_in_sample is n_cohorts * entrants. A gapped or subsampled schedule
  # holds fewer at an occasion, and scaling one by the other would put two
  # designs on one chart: 6 * 267 = 1602 against the chart's 4 * 267 = 1068.
  target <- n_prop(p = 0.5, moe = 0.031)
  six <- n_panel(target, retention = rep(0.9, 5), resp_rate = 0.8,
                 design = "rotating")
  four <- n_panel(target, retention = rep(0.9, 3), resp_rate = 0.8,
                  design = "rotating")
  pdf(tempfile())
  on.exit(dev.off())
  expect_error(plot(design_overlap("1-1-0-0-1-1"), panel = six),
               "a life with a break in it")
  expect_error(plot(design_overlap(c(1, 1, 0.5, 0.5)), panel = four),
               "subsamples a later wave")
  # the life still has to match, and a schedule the panel does represent draws
  expect_error(plot(design_overlap("4"), panel = six), "6 live cohorts")
  expect_silent(plot(design_overlap("6"), panel = six))
  expect_silent(plot(design_overlap("4"), panel = four))
})

test_that("the chart draws the launch its panel was planned with", {
  target <- n_prop(p = 0.5, moe = 0.031)
  ret <- c(0.9, 0.9)
  imm <- n_panel(target, retention = ret, resp_rate = 0.8,
                 design = "rotating", start = "immediate")
  grad <- n_panel(target, retention = ret, resp_rate = 0.8,
                  design = "rotating", start = "gradual")
  plain <- n_panel(target, retention = ret, resp_rate = 0.8,
                   design = "rotating")
  w <- design_overlap("3")$schedule
  # inherited when not stated
  expect_identical(.check_chart_start(NULL, w, imm), "immediate")
  expect_identical(.check_chart_start(NULL, w, grad), "gradual")
  # and "gradual" where the panel carries no launch, or there is no panel
  expect_identical(.check_chart_start(NULL, w, plain), "gradual")
  expect_identical(.check_chart_start(NULL, w, NULL), "gradual")
  # stating it wins, which is how one plan's two launches are compared
  expect_identical(.check_chart_start("gradual", w, imm), "gradual")
  expect_identical(.check_chart_start("immediate", w, grad), "immediate")

  pdf(tempfile())
  on.exit(dev.off())
  expect_silent(plot(design_overlap("3"), panel = imm))
  expect_silent(plot(design_overlap("3"), panel = imm, start = "gradual"))
})
