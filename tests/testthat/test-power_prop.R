test_that("power_prop solve-n matches formula", {
  p1 <- 0.30; p2 <- 0.35; power <- 0.80; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  V <- p1 * (1 - p1) + p2 * (1 - p2)
  delta <- abs(p1 - p2)
  expected <- (z_a + z_b)^2 * V / delta^2

  res <- power_prop(p1 = p1, p2 = p2)
  expect_s3_class(res, "svyplan_power")
  expect_equal(res$n, expected, tolerance = 1e-6)
  expect_equal(res$solved, "n")
  expect_equal(res$type, "proportion")
  expect_equal(res$power, 0.80)
  expect_equal(res$delta, 0.05)
})

test_that("power_prop solve-power matches formula", {
  p1 <- 0.30; p2 <- 0.35; n <- 1500; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  V <- p1 * (1 - p1) + p2 * (1 - p2)
  delta <- abs(p1 - p2)
  se <- sqrt(V / n)
  expected <- pnorm(delta / se - z_a) + pnorm(-delta / se - z_a)

  res <- power_prop(p1 = p1, p2 = p2, n = n, power = NULL)
  expect_equal(res$power, expected, tolerance = 1e-6)
  expect_equal(res$solved, "power")
  expect_equal(res$n, 1500)
})

test_that("power_prop solve-mde returns valid p2", {
  res <- power_prop(p1 = 0.30, n = 1000)
  expect_equal(res$solved, "mde")
  expect_true(res$params$p2 > 0.30)
  expect_true(res$params$p2 < 1)
  expect_equal(res$delta, res$params$p2 - 0.30, tolerance = 1e-8)
})

test_that("power_prop deff multiplier", {
  base <- power_prop(p1 = 0.30, p2 = 0.35)
  deff2 <- power_prop(p1 = 0.30, p2 = 0.35, deff = 2)
  expect_equal(deff2$n, base$n * 2, tolerance = 1e-6)
})

test_that("power_prop FPC reduces n", {
  base <- power_prop(p1 = 0.30, p2 = 0.35)
  fpc <- power_prop(p1 = 0.30, p2 = 0.35, N = 5000)
  expect_true(fpc$n < base$n)
})

test_that("power_prop panel overlap reduces n", {
  base <- power_prop(p1 = 0.30, p2 = 0.35)
  panel <- power_prop(p1 = 0.30, p2 = 0.35, overlap = 0.5, rho = 0.6)
  expect_true(panel$n < base$n)
})

test_that("power_prop one-sided needs smaller n", {
  two <- power_prop(p1 = 0.30, p2 = 0.35, sides = 2)
  one <- power_prop(p1 = 0.30, p2 = 0.35, sides = 1)
  expect_true(one$n < two$n)
})

test_that("power_prop round-trip n -> power", {
  res_n <- power_prop(p1 = 0.30, p2 = 0.35, power = 0.80)
  res_pw <- power_prop(p1 = 0.30, p2 = 0.35, n = res_n$n, power = NULL)
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("power_prop round-trip n -> mde", {
  res_n <- power_prop(p1 = 0.30, p2 = 0.35, power = 0.80)
  res_mde <- power_prop(p1 = 0.30, n = res_n$n, power = 0.80)
  expect_equal(res_mde$params$p2, 0.35, tolerance = 1e-3)
})

test_that("power_prop coercion methods", {
  res <- power_prop(p1 = 0.30, p2 = 0.35)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
  expect_equal(as.double(res), res$n)
})

test_that("power_prop validates inputs", {
  expect_error(power_prop(p1 = 0), "must be in \\(0, 1\\)")
  expect_error(power_prop(p1 = 0.3, p2 = 0.4, n = 100), "exactly one")
  expect_error(power_prop(p1 = 0.3, p2 = 0.3), "must differ")
  expect_error(power_prop(p1 = 0.3, p2 = 0.4, n = 100, power = NULL, sides = 3),
               "sides")
  expect_error(power_prop(p1 = 0.3, p2 = 0.4, n = 100, power = NULL, overlap = -1),
               "overlap")
  expect_error(power_prop(p1 = 0.3, p2 = 0.4, n = 100, power = NULL, rho = 2),
               "rho")
})

test_that("power_prop format and print", {
  res <- power_prop(p1 = 0.30, p2 = 0.35)
  expect_match(format(res), "svyplan_power")
  expect_output(print(res), "Power analysis for proportions")
  expect_output(print(res), "solved for sample size")
})
