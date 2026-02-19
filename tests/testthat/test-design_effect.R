test_that("design_effect kish formula", {
  set.seed(42)
  w <- runif(100, 1, 5)
  result <- design_effect(w, method = "kish")

  # Kish: 1 + sum((w - mean(w))^2) / n / mean(w)^2
  n <- length(w)
  expected <- 1 + sum((w - mean(w))^2) / n / mean(w)^2

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("design_effect cluster planning formula", {
  result <- design_effect(delta = 0.05, m = 25, method = "cluster")
  expect_equal(result, 1 + (25 - 1) * 0.05)
  expect_equal(result, 2.2)
})

test_that("design_effect cluster with svyplan_varcomp", {
  vc <- .new_svyplan_varcomp(
    var_between = 0.01, var_within = 1.0, delta = 0.05, k = 1.0,
    rel_var = 1.0, stages = 2L
  )
  result <- design_effect(delta = vc, m = 25, method = "cluster")
  expect_equal(result, 1 + (25 - 1) * 0.05)
})

test_that("design_effect henry method", {
  set.seed(42)
  n <- 100
  x_cal <- runif(n, 10, 100)
  y <- 5 + 0.3 * x_cal + rnorm(n, 0, 2)
  w <- runif(n, 1, 5)

  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1L)
})

test_that("design_effect spencer method", {
  set.seed(42)
  n <- 100
  p_sel <- runif(n, 0.01, 0.2)
  y <- 50 + rnorm(n, 0, 10)
  w <- 1 / p_sel

  result <- design_effect(w, y = y, p = p_sel, method = "spencer")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1L)
})

test_that("design_effect validates inputs", {
  expect_error(design_effect(delta = 0.05, method = "cluster"),
               "'delta' and 'm' are required")
  expect_error(design_effect(m = 25, method = "cluster"),
               "'delta' and 'm' are required")
  w <- runif(10, 1, 5)
  expect_error(design_effect(w, method = "henry"),
               "'y' and 'x_cal' are required")
  expect_error(design_effect(w, method = "spencer"),
               "'y' and 'p' are required")
  expect_error(design_effect(w, y = rnorm(10), method = "cr"),
               "'strvar' or 'clvar'")
})

test_that("design_effect cr requires survey package", {
  skip_if_not_installed("survey")
  set.seed(42)
  n <- 100
  w <- rep(50, n)  # Must sum to more than n for CR
  y <- rnorm(n, 50, 10)
  strvar <- rep(1:5, each = 20)
  clvar <- rep(1:25, each = 4)
  stages <- rep(2L, 5)

  result <- design_effect(w, y = y, strvar = strvar, clvar = clvar,
                          stages = stages, method = "cr")
  expect_true(is.list(result))
  expect_true("strata" %in% names(result))
  expect_true("overall" %in% names(result))
  expect_true(is.numeric(result$overall))
})
