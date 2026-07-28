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
