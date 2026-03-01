test_that("confint.svyplan_n works for proportion", {
  result <- n_prop(p = 0.3, moe = 0.05)
  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1L)
  expect_equal(ncol(ci), 2L)
  mid <- (ci[1, 1] + ci[1, 2]) / 2
  expect_equal(mid, 0.3, tolerance = 1e-6)
})

test_that("confint.svyplan_n works for mean", {
  result <- n_mean(var = 100, moe = 2, mu = 50)
  ci <- confint(result)
  mid <- (ci[1, 1] + ci[1, 2]) / 2
  expect_equal(mid, 50, tolerance = 1e-6)
})

test_that("confint.svyplan_n fails for mean without mu", {
  result <- n_mean(var = 100, moe = 2)
  expect_error(confint(result), "mu.*required")
})

test_that("confint.svyplan_n fails for multi type", {
  tgt <- data.frame(p = c(0.3, 0.5), moe = c(0.05, 0.05))
  result <- n_multi(tgt)
  expect_error(confint(result), "single-indicator")
})

test_that("confint.svyplan_n respects level argument", {
  result <- n_prop(p = 0.3, moe = 0.05)
  ci95 <- confint(result, level = 0.95)
  ci99 <- confint(result, level = 0.99)
  width95 <- ci95[1, 2] - ci95[1, 1]
  width99 <- ci99[1, 2] - ci99[1, 1]
  expect_true(width99 > width95)
})

test_that("confint.svyplan_prec works for proportion", {
  result <- prec_prop(p = 0.3, n = 400)
  ci <- confint(result)
  expect_true(is.matrix(ci))
  mid <- (ci[1, 1] + ci[1, 2]) / 2
  expect_equal(mid, 0.3, tolerance = 1e-6)
})

test_that("confint.svyplan_prec uses Wilson endpoints for method='wilson'", {
  p <- 0.1
  n <- 100
  alpha <- 0.05
  result <- prec_prop(p = p, n = n, alpha = alpha, method = "wilson")
  ci <- confint(result)

  z <- qnorm(1 - alpha / 2)
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  expected <- c(center - half, center + half)

  expect_equal(as.numeric(ci[1, ]), expected, tolerance = 1e-10)
})

test_that("confint.svyplan_n uses Wilson endpoints for method='wilson'", {
  p <- 0.1
  alpha <- 0.05
  result <- n_prop(p = p, moe = 0.06, alpha = alpha, method = "wilson")
  ci <- confint(result)

  n_eff <- result$n
  z <- qnorm(1 - alpha / 2)
  center <- (p + z^2 / (2 * n_eff)) / (1 + z^2 / n_eff)
  half <- z * sqrt(p * (1 - p) / n_eff + z^2 / (4 * n_eff^2)) / (1 + z^2 / n_eff)
  expected <- c(center - half, center + half)

  expect_equal(as.numeric(ci[1, ]), expected, tolerance = 1e-10)
})

test_that("confint.svyplan_prec logodds interval is asymmetric around p", {
  p <- 0.1
  result <- prec_prop(p = p, n = 120, method = "logodds")
  ci <- confint(result)
  lower_gap <- p - ci[1, 1]
  upper_gap <- ci[1, 2] - p
  expect_gt(abs(lower_gap - upper_gap), 1e-4)
})

test_that("confint.svyplan_prec works for mean", {
  result <- prec_mean(var = 100, n = 400, mu = 50)
  ci <- confint(result)
  mid <- (ci[1, 1] + ci[1, 2]) / 2
  expect_equal(mid, 50, tolerance = 1e-6)
})

test_that("confint.svyplan_prec fails for mean without mu", {
  result <- prec_mean(var = 100, n = 400)
  expect_error(confint(result), "mu.*required")
})

test_that("confint.svyplan_prec fails for cluster type", {
  result <- prec_cluster(n = c(50, 12), delta = 0.05)
  expect_error(confint(result), "not supported")
})

test_that("confint.svyplan_n rejects invalid level", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_error(confint(result, level = 1.2), "level.*\\(0, 1\\)")
  expect_error(confint(result, level = 0), "level.*\\(0, 1\\)")
  expect_error(confint(result, level = -1), "level.*\\(0, 1\\)")
})

test_that("confint.svyplan_prec rejects invalid level", {
  result <- prec_prop(p = 0.3, n = 400)
  expect_error(confint(result, level = 1.2), "level.*\\(0, 1\\)")
  expect_error(confint(result, level = 0), "level.*\\(0, 1\\)")
  expect_error(confint(result, level = -1), "level.*\\(0, 1\\)")
})

test_that("Wilson/logodds confint widens with higher level", {
  res_w <- n_prop(p = 0.3, moe = 0.05, method = "wilson")
  ci90 <- confint(res_w, level = 0.90)
  ci99 <- confint(res_w, level = 0.99)
  expect_true((ci99[1, 2] - ci99[1, 1]) > (ci90[1, 2] - ci90[1, 1]))

  res_l <- n_prop(p = 0.2, moe = 0.04, method = "logodds")
  ci90_l <- confint(res_l, level = 0.90)
  ci99_l <- confint(res_l, level = 0.99)
  expect_true((ci99_l[1, 2] - ci99_l[1, 1]) > (ci90_l[1, 2] - ci90_l[1, 1]))
})

test_that("logodds confint: finite N narrows interval vs N=Inf", {
  res_inf <- prec_prop(p = 0.2, n = 500, method = "logodds")
  res_fin <- prec_prop(p = 0.2, n = 500, method = "logodds", N = 5000)
  ci_inf <- confint(res_inf)
  ci_fin <- confint(res_fin)
  expect_true((ci_fin[1, 2] - ci_fin[1, 1]) < (ci_inf[1, 2] - ci_inf[1, 1]))
})

test_that("Wilson confint: higher deff widens CI", {
  res1 <- prec_prop(p = 0.3, n = 400, method = "wilson", deff = 1)
  res2 <- prec_prop(p = 0.3, n = 400, method = "wilson", deff = 2)
  ci1 <- confint(res1)
  ci2 <- confint(res2)
  expect_true((ci2[1, 2] - ci2[1, 1]) > (ci1[1, 2] - ci1[1, 1]))
})
