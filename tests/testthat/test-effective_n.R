test_that("effective_n kish uses direct formula", {
  set.seed(42)
  w <- runif(100, 1, 5)
  result <- effective_n(w, method = "kish")

  # Direct: sum(w)^2 / sum(w^2)
  expected <- sum(w)^2 / sum(w^2)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("effective_n kish consistent with design_effect", {
  set.seed(42)
  w <- runif(100, 1, 5)
  n_eff <- effective_n(w, method = "kish")
  deff <- design_effect(w, method = "kish")

  # n_eff should equal n / deff
  expect_equal(n_eff, length(w) / deff, tolerance = 1e-6)
})

test_that("effective_n cluster planning", {
  result <- effective_n(delta = 0.05, psu_size =25, n = 800, method = "cluster")
  deff <- 1 + (25 - 1) * 0.05  # 2.2
  expect_equal(result, 800 / deff, tolerance = 1e-10)
})

test_that("effective_n cluster validates n argument", {
  expect_error(effective_n(delta = 0.05, psu_size =25, method = "cluster"),
               "'n' is required")
})

test_that("effective_n with equal weights equals n", {
  w <- rep(1, 50)
  result <- effective_n(w, method = "kish")
  expect_equal(result, 50, tolerance = 1e-10)
})

test_that("effective_n rejects Inf weights", {
  expect_error(effective_n(c(1, Inf), method = "kish"), "finite")
})

test_that("effective_n rejects -Inf weights", {
  expect_error(effective_n(c(1, -Inf), method = "kish"), "finite")
})

test_that("effective_n kish with highly variable weights", {
  w <- c(rep(1, 99), 100)
  result <- effective_n(w, method = "kish")
  expected <- sum(w)^2 / sum(w^2)
  expect_equal(result, expected, tolerance = 1e-10)
  expect_true(result < length(w))
})

test_that("effective_n kish with length-1 weight", {
  result <- effective_n(5, method = "kish")
  expect_equal(result, 1)
})

test_that("effective_n cluster with small psu_size", {
  result <- effective_n(delta = 0.05, psu_size = 2, n = 100, method = "cluster")
  deff <- 1 + (2 - 1) * 0.05
  expect_equal(result, 100 / deff, tolerance = 1e-10)
})

test_that("effective_n cluster with delta = 0 equals n", {
  result <- effective_n(delta = 0, psu_size = 25, n = 800, method = "cluster")
  expect_equal(result, 800, tolerance = 1e-10)
})

test_that("effective_n cluster with delta = 1 (max homogeneity)", {
  result <- effective_n(delta = 1, psu_size = 25, n = 800, method = "cluster")
  deff <- 1 + (25 - 1) * 1
  expect_equal(result, 800 / deff, tolerance = 1e-10)
  expect_equal(result, 32, tolerance = 1e-10)
})

test_that("effective_n cluster validates psu_size", {
  expect_error(effective_n(delta = 0.05, n = 800, method = "cluster"),
               "psu_size")
})

test_that("effective_n cluster validates delta", {
  expect_error(effective_n(delta = -0.1, psu_size = 25, n = 800, method = "cluster"))
  expect_error(effective_n(delta = 1.5, psu_size = 25, n = 800, method = "cluster"))
})

test_that("effective_n is inverse of design_effect * n", {
  set.seed(123)
  w <- runif(200, 1, 10)
  for (m in c("kish")) {
    n_eff <- effective_n(w, method = m)
    deff <- design_effect(w, method = m)
    expect_equal(n_eff, length(w) / deff, tolerance = 1e-6)
  }
})
