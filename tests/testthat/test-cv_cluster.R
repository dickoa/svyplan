test_that("cv_cluster 2-stage formula", {
  result <- cv_cluster(n = c(50, 12), delta = 0.05)
  expected <- sqrt(1 / (50 * 12) * 1 * (1 + 0.05 * (12 - 1)))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cv_cluster 3-stage formula", {
  result <- cv_cluster(n = c(50, 12, 8), delta = c(0.01, 0.05))
  expected <- sqrt(1 / (50 * 12 * 8) * (1 * 0.01 * 12 * 8 +
                                          1 * (1 + 0.05 * (8 - 1))))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cv_cluster with custom rel_var and k", {
  result <- cv_cluster(n = c(50, 12), delta = 0.05, rel_var = 1.5, k = 1.2)
  expected <- sqrt(1.5 / (50 * 12) * 1.2 * (1 + 0.05 * (12 - 1)))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cv_cluster accepts svyplan_varcomp", {
  vc <- .new_svyplan_varcomp(
    var_between = 0.01, var_within = 1.0, delta = 0.0099, k = 1.01,
    rel_var = 1.5, stages = 2L
  )
  result <- cv_cluster(n = c(50, 12), delta = vc)
  expected <- sqrt(1.5 / (50 * 12) * 1.01 * (1 + 0.0099 * (12 - 1)))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cv_cluster validates inputs", {
  expect_error(cv_cluster(n = 10, delta = 0.05), "length >= 2")
  expect_error(cv_cluster(n = c(50, 12, 8, 5), delta = c(0.01, 0.02, 0.03)),
               "not yet supported")
  expect_error(cv_cluster(n = c(50, 12), delta = c(0.05, 0.10)),
               "must have length")
  expect_error(cv_cluster(n = c(50, -1), delta = 0.05),
               "must be positive")
})
