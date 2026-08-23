## Regression for a pre-existing defect: predict() dropped a 'df' the object
## itself was built with, so a predicted row reported a normal-quantile
## interval for an object whose own interval used a t quantile. The proportion
## branches always forwarded it, which is why it stayed invisible.
##
## Both directions are covered here because they are separate branches:
## predict.svyplan_n() and predict.svyplan_prec().

test_that("predict on a sized mean keeps the object's own df", {
  obj <- n_mean(var = 4, mu = 10, cv = 0.05, df = 12)
  grid <- predict(obj, data.frame(deff = c(1, 1.5)))

  # At deff = 1 the row is the object, so its interval must be the object's.
  expect_equal(grid$n[1L], obj$n, tolerance = 1e-10)
  expect_equal(grid$moe[1L], obj$moe, tolerance = 1e-10)
  expect_equal(grid$rmoe[1L], obj$rmoe, tolerance = 1e-10)
  # And that interval is wider than the normal one it used to report.
  expect_gt(grid$moe[1L], qnorm(0.975) * grid$se[1L])
  expect_equal(grid$moe[1L], qt(0.975, 12) * grid$se[1L], tolerance = 1e-10)
})

test_that("predict on a mean precision result keeps the object's own df", {
  obj <- prec_mean(var = 4, mu = 10, n = 200, df = 12)
  grid <- predict(obj, data.frame(n = c(200, 400)))

  expect_equal(grid$se[1L], obj$se, tolerance = 1e-10)
  expect_equal(grid$moe[1L], obj$moe, tolerance = 1e-10)
  expect_equal(grid$moe[1L], qt(0.975, 12) * grid$se[1L], tolerance = 1e-10)
  expect_gt(grid$moe[1L], qnorm(0.975) * grid$se[1L])
})

test_that("df can be varied in a mean grid, both directions", {
  sized <- n_mean(var = 4, mu = 10, moe = 0.5)
  by_df <- predict(sized, data.frame(df = c(3, 10, 1e6)))
  expect_true(all(diff(by_df$n) < 0))

  prec <- prec_mean(var = 4, mu = 10, n = 200)
  prec_df <- predict(prec, data.frame(df = c(3, 10, 1e6)))
  expect_true(all(diff(prec_df$moe) < 0))
  expect_equal(prec_df$se[1L], prec_df$se[3L])
})

test_that("a mean without df is unaffected", {
  obj <- n_mean(var = 4, mu = 10, cv = 0.05)
  grid <- predict(obj, data.frame(deff = c(1, 1.5)))

  expect_equal(grid$moe[1L], obj$moe, tolerance = 1e-10)
  expect_equal(grid$moe[1L], qnorm(0.975) * grid$se[1L], tolerance = 1e-10)
})

test_that("predict on a ratio keeps the object's own df, both directions", {
  sized <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                   moe = 0.2, df = 12)
  grid <- predict(sized, data.frame(deff = c(1, 1.5)))
  expect_equal(grid$n[1L], sized$n, tolerance = 1e-10)
  expect_equal(grid$moe[1L], qt(0.975, 12) * grid$se[1L], tolerance = 1e-10)

  prec <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                     component_cor = 0.7, df = 12)
  prec_grid <- predict(prec, data.frame(n = c(400, 900)))
  expect_equal(prec_grid$moe[1L], prec$moe, tolerance = 1e-10)
  expect_equal(prec_grid$moe[1L], qt(0.975, 12) * prec_grid$se[1L],
               tolerance = 1e-10)
})
