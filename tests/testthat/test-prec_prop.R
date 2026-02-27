test_that("prec_prop computes SE, MOE, CV for wald method", {
  result <- prec_prop(p = 0.3, n = 400)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "proportion")
  expect_equal(result$method, "wald")

  z <- qnorm(0.975)
  se_exp <- sqrt(0.3 * 0.7 / 400)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
  expect_equal(result$moe, z * se_exp, tolerance = 1e-6)
  expect_equal(result$cv, se_exp / 0.3, tolerance = 1e-6)
})

test_that("prec_prop with FPC uses Cochran form", {
  result <- prec_prop(p = 0.3, n = 400, N = 5000)
  n_eff <- 400
  fpc <- (5000 - n_eff) / (5000 - 1)
  se_exp <- sqrt(0.3 * 0.7 * fpc / n_eff)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

test_that("prec_prop with deff", {
  base <- prec_prop(p = 0.3, n = 400)
  deff2 <- prec_prop(p = 0.3, n = 400, deff = 2)
  expect_equal(deff2$se, base$se * sqrt(2), tolerance = 1e-6)
})

test_that("prec_prop with resp_rate", {
  base <- prec_prop(p = 0.3, n = 400)
  rr <- prec_prop(p = 0.3, n = 400, resp_rate = 0.8)
  n_eff <- 400 * 0.8
  se_exp <- sqrt(0.3 * 0.7 / n_eff)
  expect_equal(rr$se, se_exp, tolerance = 1e-6)
})

test_that("prec_prop validates inputs", {
  expect_error(prec_prop(p = 0, n = 400), "must be in \\(0, 1\\)")
  expect_error(prec_prop(p = 1, n = 400), "must be in \\(0, 1\\)")
  expect_error(prec_prop(p = 0.3, n = -1), "must be positive")
  expect_s3_class(prec_prop(p = 0.3, n = 400, deff = 0.5), "svyplan_prec")
  expect_error(prec_prop(p = 0.3, n = 400, deff = 0), "must be positive")
  expect_error(prec_prop(p = 0.3, n = 400, deff = -1), "must be positive")
  expect_error(prec_prop(p = 0.3, n = 400, resp_rate = 0), "resp_rate")
  expect_error(prec_prop(p = 0.3, n = 400, resp_rate = 1.5), "resp_rate")
})

test_that("prec_prop prints correctly", {
  result <- prec_prop(p = 0.3, n = 400)
  out <- capture.output(print(result))
  expect_match(out[1], "Sampling precision for proportion")
  expect_match(out[2], "n = 400")
  expect_match(out[3], "se =")
})

test_that("prec_prop format returns string", {
  result <- prec_prop(p = 0.3, n = 400)
  expect_true(is.character(format(result)))
  expect_match(format(result), "svyplan_prec")
})

test_that("prec_prop logodds works for extreme p", {
  res <- prec_prop(p = 0.9, n = 20, method = "logodds")
  expect_true(res$moe > 0 && res$moe < 0.5)
  expect_true(res$se > 0)

  res2 <- prec_prop(p = 0.05, n = 100, method = "logodds")
  expect_true(res2$moe > 0 && res2$moe < 0.5)
})

test_that("prec_prop logodds extreme p round-trip", {
  for (pp in c(0.05, 0.95)) {
    s1 <- n_prop(p = pp, moe = 0.04, method = "logodds")
    p1 <- prec_prop(s1)
    s2 <- n_prop(p1)
    expect_equal(s2$n, s1$n, tolerance = 1)
  }
})

test_that("prec_prop logodds census returns moe = 0 with warning", {
  expect_warning(
    res <- prec_prop(p = 0.5, n = 50, N = 50, method = "logodds"),
    "effective sample size >= population size"
  )
  expect_equal(res$moe, 0)
})
