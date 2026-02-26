test_that("n_mean MOE mode", {
  result <- n_mean(var = 100, moe = 2)
  z <- qnorm(0.975)
  expected <- z^2 * 100 / (4 + z^2 * 100 / Inf)
  expect_equal(result$n, expected, tolerance = 1e-6)
  expect_s3_class(result, "svyplan_n")
  expect_equal(result$type, "mean")
})

test_that("n_mean CV mode", {
  result <- n_mean(var = 100, mu = 50, cv = 0.05)
  CVpop <- sqrt(100) / 50  # 0.2
  expected <- CVpop^2 / (0.05^2 + CVpop^2 / Inf)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_mean MOE with FPC", {
  result <- n_mean(var = 100, moe = 2, N = 5000)
  z <- qnorm(0.975)
  expected <- z^2 * 100 / (4 + z^2 * 100 / 5000)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_mean CV with FPC", {
  result <- n_mean(var = 100, mu = 50, cv = 0.05, N = 5000)
  CVpop <- sqrt(100) / 50
  expected <- CVpop^2 / (0.05^2 + CVpop^2 / 5000)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_mean deff multiplier works", {
  base <- n_mean(var = 100, moe = 2)
  with_deff <- n_mean(var = 100, moe = 2, deff = 1.5)
  expect_equal(with_deff$n, base$n * 1.5, tolerance = 1e-6)
})

test_that("n_mean returns correct S3 class", {
  result <- n_mean(var = 100, moe = 2)
  expect_true(inherits(result, "svyplan_n"))
  expect_equal(as.integer(result), ceiling(result$n))
  expect_equal(as.double(result), result$n)
})

test_that("n_mean validates inputs", {
  expect_error(n_mean(var = -1, moe = 2), "positive")
  expect_error(n_mean(var = 100), "specify exactly one")
  expect_error(n_mean(var = 100, moe = 2, cv = 0.05), "specify exactly one")
  expect_error(n_mean(var = 100, cv = 0.05), "'mu' is required")
  expect_s3_class(n_mean(var = 100, moe = 2, deff = 0.5), "svyplan_n")
  expect_error(n_mean(var = 100, moe = 2, deff = 0), "must be positive")
  expect_error(n_mean(var = 100, moe = 2, deff = -1), "must be positive")
})

test_that("svyplan_n has se/moe/cv fields for mean", {
  res <- n_mean(var = 100, moe = 2, mu = 50)
  expect_true(is.numeric(res$se))
  expect_true(res$se > 0)
  expect_true(res$moe > 0)
  expect_true(res$cv > 0)
})

test_that("mean moe round-trip: moe from n matches target", {
  res <- n_mean(var = 100, moe = 2)
  expect_equal(res$moe, 2, tolerance = 1e-6)
})

test_that("mean moe round-trip with finite N", {
  res <- n_mean(var = 100, moe = 2, N = 5000)
  expect_equal(res$moe, 2, tolerance = 1e-6)
})

test_that("mean cv is NA when mu not provided", {
  res <- n_mean(var = 100, moe = 2)
  expect_true(is.na(res$cv))
})

test_that("n_mean rejects N = 1", {
  expect_error(n_mean(var = 100, moe = 2, N = 1), "greater than 1")
})
