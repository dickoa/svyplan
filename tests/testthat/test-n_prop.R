test_that("n_prop wald MOE mode", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_s3_class(result, "svyplan_n")
  expect_equal(result$type, "proportion")
  expect_equal(result$method, "wald")

  # Infinite pop: a=1, z=1.96, n = 1*1.96^2*0.3*0.7 / (0.05^2 + 1.96^2*0.3*0.7/Inf)
  z <- qnorm(0.975)
  expected <- z^2 * 0.3 * 0.7 / 0.05^2
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop wald CV mode", {
  result <- n_prop(p = 0.5, cv = 0.10)
  # a=1, q/p = 1, n = 1 * 1 / (0.01 + 1/Inf) = 100
  expect_equal(result$n, 100, tolerance = 1e-6)
})

test_that("n_prop wald MOE with FPC", {
  result <- n_prop(p = 0.3, moe = 0.05, N = 1000)
  z <- qnorm(0.975)
  a <- 1000 / 999
  expected <- a * z^2 * 0.3 * 0.7 / (0.05^2 + z^2 * 0.3 * 0.7 / 999)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop wald CV with FPC", {
  result <- n_prop(p = 0.5, cv = 0.10, N = 5000)
  a <- 5000 / 4999
  expected <- a * 1 / (0.01 + 1 / 4999)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop wilson method", {
  result <- n_prop(p = 0.3, moe = 0.05, method = "wilson")
  z <- qnorm(0.975)
  q <- 0.7
  e <- 0.05
  rad <- e^2 - 0.3 * q * (4 * e^2 - 0.3 * q)
  expected <- (0.3 * q - 2 * e^2 + sqrt(rad)) * (z / e)^2 / 2
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop logodds method", {
  result <- n_prop(p = 0.3, moe = 0.05, method = "logodds")
  z <- qnorm(0.975)
  q <- 0.7
  e <- 0.05
  kk <- q / 0.3
  rad <- e^2 * (1 + kk^2)^2 + kk^2 * (1 - 2 * e) * (1 + 2 * e)
  x <- (e * (1 + kk^2) + sqrt(rad)) / (kk * (1 - 2 * e))
  expected <- 1 / ((sqrt(0.3 * q) / z * log(x))^2)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop deff multiplier works", {
  base <- n_prop(p = 0.3, moe = 0.05)
  with_deff <- n_prop(p = 0.3, moe = 0.05, deff = 2)
  expect_equal(with_deff$n, base$n * 2, tolerance = 1e-6)
})

test_that("n_prop returns correct S3 class", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_true(inherits(result, "svyplan_n"))
  expect_true(is.list(result))
  expect_equal(as.integer(result), ceiling(result$n))
  expect_equal(as.double(result), result$n)
})

test_that("n_prop validates inputs", {
  expect_error(n_prop(p = 0, moe = 0.05), "must be in \\(0, 1\\)")
  expect_error(n_prop(p = 1, moe = 0.05), "must be in \\(0, 1\\)")
  expect_error(n_prop(p = 0.3), "specify exactly one")
  expect_error(n_prop(p = 0.3, moe = 0.05, cv = 0.10), "specify exactly one")
  expect_s3_class(n_prop(p = 0.3, moe = 0.05, deff = 0.5), "svyplan_n")
  expect_error(n_prop(p = 0.3, moe = 0.05, deff = 0), "must be positive")
  expect_error(n_prop(p = 0.3, moe = 0.05, deff = -1), "must be positive")
  expect_error(n_prop(p = 0.3, moe = 0.05, N = -1), "positive")
  expect_error(n_prop(p = 0.3, cv = 0.10, method = "wilson"), "Wilson.*moe")
  expect_error(n_prop(p = 0.3, cv = 0.10, method = "logodds"), "Log-odds.*moe")
})

test_that("svyplan_n has se/moe/cv fields for proportion", {
  res <- n_prop(p = 0.3, moe = 0.05)
  expect_true(!is.null(res$se))
  expect_true(!is.null(res$moe))
  expect_true(!is.null(res$cv))
  expect_true(is.numeric(res$se))
  expect_true(res$se > 0)
  expect_true(res$moe > 0)
  expect_true(res$cv > 0)
})

test_that("se/moe/cv use Cochran FPC for proportion (infinite pop)", {
  res <- n_prop(p = 0.3, moe = 0.05)
  z <- qnorm(0.975)
  n_eff <- res$n
  se_expected <- sqrt(0.3 * 0.7 / n_eff)
  expect_equal(res$se, se_expected, tolerance = 1e-6)
  expect_equal(res$moe, z * se_expected, tolerance = 1e-6)
  expect_equal(res$cv, se_expected / 0.3, tolerance = 1e-6)
})

test_that("se/moe/cv use Cochran FPC for proportion (finite pop)", {
  res <- n_prop(p = 0.3, moe = 0.05, N = 500)
  z <- qnorm(0.975)
  n_eff <- res$n
  fpc <- (500 - n_eff) / (500 - 1)
  se_expected <- sqrt(0.3 * 0.7 * fpc / n_eff)
  expect_equal(res$se, se_expected, tolerance = 1e-6)
  expect_equal(res$moe, z * se_expected, tolerance = 1e-6)
})

test_that("proportion FPC round-trip: moe from n matches target", {
  res <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$moe, 0.05, tolerance = 1e-6)
})

test_that("proportion FPC round-trip with finite N", {
  res <- n_prop(p = 0.3, moe = 0.05, N = 1000)
  expect_equal(res$moe, 0.05, tolerance = 1e-6)
})
