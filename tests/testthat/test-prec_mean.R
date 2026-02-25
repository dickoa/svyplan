test_that("prec_mean computes SE, MOE, CV", {
  result <- prec_mean(var = 100, n = 400, mu = 50)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "mean")

  se_exp <- sqrt(100 / 400)
  z <- qnorm(0.975)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
  expect_equal(result$moe, z * se_exp, tolerance = 1e-6)
  expect_equal(result$cv, se_exp / 50, tolerance = 1e-6)
})

test_that("prec_mean without mu gives NA cv", {
  result <- prec_mean(var = 100, n = 400)
  expect_false(is.na(result$se))
  expect_false(is.na(result$moe))
  expect_true(is.na(result$cv))
})

test_that("prec_mean with FPC", {
  result <- prec_mean(var = 100, n = 400, N = 5000)
  f <- 400 / 5000
  se_exp <- sqrt(100 * (1 - f) / 400)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

test_that("prec_mean with deff and resp_rate", {
  result <- prec_mean(var = 100, n = 400, deff = 2, resp_rate = 0.8)
  n_eff <- 400 * 0.8 / 2
  se_exp <- sqrt(100 / n_eff)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

test_that("prec_mean validates inputs", {
  expect_error(prec_mean(var = -1, n = 400), "must be positive")
  expect_error(prec_mean(var = 100, n = -1), "must be positive")
})
