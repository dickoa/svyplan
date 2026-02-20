test_that("n_cluster 2-stage budget mode", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_s3_class(result, "svyplan_cluster")

  n2_opt <- sqrt(500 / 50 * (1 - 0.05) / 0.05)
  n1_opt <- 100000 / (500 + 50 * n2_opt)
  cv_expected <- sqrt(1 / (n1_opt * n2_opt) * 1 * (1 + 0.05 * (n2_opt - 1)))

  expect_equal(result$n[["n2"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n1"]], n1_opt, tolerance = 1e-6)
  expect_equal(result$cv, cv_expected, tolerance = 1e-6)
  expect_equal(result$cost, 100000)
  expect_equal(result$stages, 2L)
})

test_that("n_cluster 2-stage CV mode", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)

  n2_opt <- sqrt(500 / 50 * (1 - 0.05) / 0.05)
  n1_opt <- 1 * 1 * (1 + 0.05 * (n2_opt - 1)) / (n2_opt * 0.05^2)
  cost_expected <- 500 * n1_opt + 50 * n1_opt * n2_opt

  expect_equal(result$n[["n2"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n1"]], n1_opt, tolerance = 1e-6)
  expect_equal(result$cost, cost_expected, tolerance = 1e-4)
})

test_that("n_cluster 2-stage fixed m budget mode", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05,
                      budget = 100000, m = 40)

  n2_expected <- (100000 - 500 * 40) / (50 * 40)
  cv_expected <- sqrt(1 * 1 / (40 * n2_expected) *
                        (1 + 0.05 * (n2_expected - 1)))

  expect_equal(result$n[["n1"]], 40, tolerance = 1e-6)
  expect_equal(result$n[["n2"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$cv, cv_expected, tolerance = 1e-6)
})

test_that("n_cluster 2-stage fixed m CV mode", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05, m = 40)

  n2_expected <- (1 - 0.05) / (0.05^2 * 40 / (1 * 1) - 0.05)
  cost_expected <- 500 * 40 + 50 * 40 * n2_expected

  expect_equal(result$n[["n1"]], 40, tolerance = 1e-6)
  expect_equal(result$n[["n2"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$cost, cost_expected, tolerance = 1e-4)
})

test_that("n_cluster 3-stage budget mode", {
  result <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                      budget = 500000)

  n3_opt <- sqrt((1 - 0.05) / 0.05 * 100 / 50)
  n2_opt <- 1 / n3_opt * sqrt((1 - 0.05) / 0.01 * 500 / 50 * 1 / 1)
  n1_opt <- 500000 / (500 + 100 * n2_opt + 50 * n2_opt * n3_opt)

  expect_equal(result$n[["n3"]], n3_opt, tolerance = 1e-6)
  expect_equal(result$n[["n2"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n1"]], n1_opt, tolerance = 1e-6)
  expect_equal(result$stages, 3L)
})

test_that("n_cluster 3-stage CV mode", {
  result <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                      cv = 0.05)

  n3_opt <- sqrt((1 - 0.05) / 0.05 * 100 / 50)
  n2_opt <- 1 / n3_opt * sqrt((1 - 0.05) / 0.01 * 500 / 50 * 1 / 1)
  n1_opt <- 1 / (0.05^2 * n2_opt * n3_opt) *
    (1 * 0.01 * n2_opt * n3_opt + 1 * (1 + 0.05 * (n3_opt - 1)))

  expect_equal(result$n[["n3"]], n3_opt, tolerance = 1e-6)
  expect_equal(result$n[["n2"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n1"]], n1_opt, tolerance = 1e-6)
})

test_that("n_cluster accepts svyplan_varcomp", {
  vc <- .new_svyplan_varcomp(
    var_between = 0.01, var_within = 1.0, delta = 0.05, k = 1.0,
    rel_var = 1.0, stages = 2L
  )
  result <- n_cluster(cost = c(500, 50), delta = vc, budget = 100000)
  expect_s3_class(result, "svyplan_cluster")
  expect_equal(result$stages, 2L)
})

test_that("n_cluster round-trips with cv_cluster", {
  plan <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  cv_check <- cv_cluster(n = unname(plan$n), delta = 0.05)
  expect_equal(cv_check, plan$cv, tolerance = 1e-6)
})

test_that("n_cluster rejects boundary delta values", {
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0, budget = 1e5),
    "\\(0, 1\\)"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 1, cv = 0.05),
    "\\(0, 1\\)"
  )
})

test_that("n_cluster validates inputs", {
  expect_error(n_cluster(cost = 500, delta = 0.05, budget = 100000),
               "length >= 2")
  expect_error(n_cluster(cost = c(500, 50), delta = 0.05),
               "specify exactly one")
  expect_error(n_cluster(cost = c(500, 50), delta = 0.05,
                         cv = 0.05, budget = 100000),
               "specify exactly one")
  expect_error(n_cluster(cost = c(500, 50, 20, 10), delta = c(0.01, 0.02, 0.03),
                         budget = 100000),
               "not yet supported")
})
