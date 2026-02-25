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
