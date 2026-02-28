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
