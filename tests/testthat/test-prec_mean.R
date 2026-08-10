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

test_that("prec_mean with combined deff, FPC, resp_rate", {
  result <- prec_mean(var = 100, n = 400, mu = 50,
                      deff = 1.5, N = 5000, resp_rate = 0.9)
  n_net <- 400 * 0.9
  fpc <- 1 - n_net / 5000
  se_exp <- sqrt(1.5 * 100 * fpc / n_net)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
  expect_equal(result$moe, qnorm(0.975) * se_exp, tolerance = 1e-6)
  expect_equal(result$cv, se_exp / 50, tolerance = 1e-6)
})

test_that("prec_mean with small N (FPC dominant)", {
  result <- prec_mean(var = 100, n = 50, N = 100)
  n_eff <- 50
  fpc <- 1 - 50 / 100
  se_exp <- sqrt(100 * fpc / n_eff)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

test_that("prec_mean with alpha = 0.10", {
  result <- prec_mean(var = 100, n = 400, alpha = 0.10)
  z <- qnorm(0.95)
  se_exp <- sqrt(100 / 400)
  expect_equal(result$moe, z * se_exp, tolerance = 1e-6)
})

test_that("prec_mean with deff < 1 (efficient design)", {
  base <- prec_mean(var = 100, n = 400)
  eff <- prec_mean(var = 100, n = 400, deff = 0.7)
  expect_true(eff$se < base$se)
})

test_that("prec_mean round-trip with n_mean moe mode", {
  s <- n_mean(var = 100, moe = 2)
  p <- prec_mean(var = 100, n = s$n)
  expect_equal(p$moe, 2, tolerance = 1e-6)
})

test_that("prec_mean round-trip with n_mean cv mode", {
  s <- n_mean(var = 100, mu = 50, cv = 0.05)
  p <- prec_mean(var = 100, n = s$n, mu = 50)
  expect_equal(p$cv, 0.05, tolerance = 1e-6)
})

test_that("prec_mean round-trip with FPC and deff", {
  s <- n_mean(var = 100, moe = 2, N = 5000, deff = 1.5)
  p <- prec_mean(var = 100, n = s$n, N = 5000, deff = 1.5)
  expect_equal(p$moe, 2, tolerance = 1e-4)
})

test_that("prec_mean round-trip via S3 dispatch", {
  s <- n_mean(var = 100, mu = 50, moe = 2)
  p <- prec_mean(s)
  expect_equal(p$moe, 2, tolerance = 1e-6)
  s2 <- n_mean(p)
  expect_equal(s2$n, s$n, tolerance = 1e-6)
})

test_that("prec_mean large n approaches zero SE", {
  result <- prec_mean(var = 100, n = 1e10)
  expect_true(result$se < 0.001)
})

test_that("prec_mean census gives zero SE", {
  result <- suppressWarnings(prec_mean(var = 100, n = 500, N = 500))
  expect_equal(result$se, 0)
  expect_equal(result$moe, 0)
})

test_that("prec_mean validates inputs", {
  expect_error(prec_mean(var = -1, n = 400), "must be positive")
  expect_error(prec_mean(var = 100, n = -1), "must be positive")
  expect_error(prec_mean(var = 100, n = 400, alpha = 0), "alpha")
  expect_error(prec_mean(var = 100, n = 400, alpha = 1), "alpha")
  expect_error(prec_mean(var = 100, n = 400, deff = 0), "deff")
  expect_error(prec_mean(var = 100, n = 400, deff = -1), "deff")
  expect_error(prec_mean(var = 100, n = 400, resp_rate = 0), "resp_rate")
  expect_error(prec_mean(var = 100, n = 400, resp_rate = 1.5), "resp_rate")
  expect_error(prec_mean(var = 100, n = 400, N = -1), "must be greater than 1")
  expect_error(prec_mean(var = 100, n = 400, mu = 0), "must not be zero")
})

## Solving for the mean

test_that("solving for mu inverts the precision it reports", {
  r <- prec_mean(var = 25, n = 1500, cv = 0.10, N = 2e6)
  expect_identical(r$solved, "mu")
  expect_equal(r$cv, 0.10)
  expect_equal(r$params$mu, r$se / 0.10)
  expect_equal(prec_mean(var = 25, n = 1500, mu = r$params$mu, N = 2e6)$cv, 0.10)
})

test_that("mu and cv are mutually exclusive but both stay optional", {
  expect_error(prec_mean(var = 25, n = 400, mu = 5, cv = 0.1),
               "at most one of 'mu' or 'cv'")
  bare <- prec_mean(var = 25, n = 400)
  expect_true(is.na(bare$cv))
  expect_null(bare$solved)
})

test_that("predict on a solved mean varies cv, not mu", {
  r <- prec_mean(var = 25, n = 1500, cv = 0.10)
  grid <- predict(r, expand.grid(n = c(500, 1500)))
  expect_true("mu" %in% names(grid))
  expect_equal(grid$cv, rep(0.10, 2))
  expect_equal(grid$mu[2], r$params$mu)
})
