test_that("power_mean solve-n matches formula", {
  delta <- 5; var <- 100; power <- 0.80; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  V <- 2 * var
  expected <- (z_a + z_b)^2 * V / delta^2

  res <- power_mean(delta = delta, var = var)
  expect_s3_class(res, "svyplan_power")
  expect_equal(res$n, expected, tolerance = 1e-6)
  expect_equal(res$solved, "n")
  expect_equal(res$type, "mean")
  expect_equal(res$delta, 5)
})

test_that("power_mean solve-power matches formula", {
  delta <- 5; var <- 100; n <- 200; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  V <- 2 * var
  se <- sqrt(V / n)
  expected <- pnorm(delta / se - z_a) + pnorm(-delta / se - z_a)

  res <- power_mean(delta = delta, var = var, n = n, power = NULL)
  expect_equal(res$power, expected, tolerance = 1e-6)
  expect_equal(res$solved, "power")
})

test_that("power_mean solve-mde matches analytical formula", {
  var <- 100; n <- 500; power <- 0.80; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  V <- 2 * var
  expected <- (z_a + z_b) * sqrt(V / n)

  res <- power_mean(var = var, n = n)
  expect_equal(res$delta, expected, tolerance = 1e-6)
  expect_equal(res$solved, "mde")
})

test_that("power_mean deff multiplier", {
  base <- power_mean(delta = 5, var = 100)
  deff2 <- power_mean(delta = 5, var = 100, deff = 2)
  expect_equal(deff2$n, base$n * 2, tolerance = 1e-6)
})

test_that("power_mean FPC reduces n", {
  base <- power_mean(delta = 5, var = 100)
  fpc <- power_mean(delta = 5, var = 100, N = 2000)
  expect_true(fpc$n < base$n)
})

test_that("power_mean panel overlap reduces n", {
  base <- power_mean(delta = 5, var = 100)
  panel <- power_mean(delta = 5, var = 100, overlap = 0.5, rho = 0.6)
  expect_true(panel$n < base$n)
})

test_that("power_mean one-sided needs smaller n", {
  two <- power_mean(delta = 5, var = 100, sides = 2)
  one <- power_mean(delta = 5, var = 100, sides = 1)
  expect_true(one$n < two$n)
})

test_that("power_mean round-trip n -> power", {
  res_n <- power_mean(delta = 5, var = 100, power = 0.80)
  res_pw <- power_mean(delta = 5, var = 100, n = res_n$n, power = NULL)
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("power_mean round-trip n -> mde", {
  res_n <- power_mean(delta = 5, var = 100, power = 0.80)
  res_mde <- power_mean(var = 100, n = res_n$n, power = 0.80)
  expect_equal(res_mde$delta, 5, tolerance = 1e-3)
})

test_that("power_mean coercion methods", {
  res <- power_mean(delta = 5, var = 100)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
  expect_equal(as.double(res), res$n)
})

test_that("power_mean validates inputs", {
  expect_error(power_mean(var = 100), "exactly one")
  expect_error(power_mean(delta = 5, var = 100, n = 200), "exactly one")
  expect_error(power_mean(delta = 5, var = 100, sides = 3), "sides")
  expect_error(power_mean(delta = 5, var = 100, overlap = -1), "overlap")
})

test_that("power_mean format and print", {
  res <- power_mean(delta = 5, var = 100)
  expect_match(format(res), "svyplan_power")
  expect_output(print(res), "Power analysis for means")
  expect_output(print(res), "solved for sample size")
})
