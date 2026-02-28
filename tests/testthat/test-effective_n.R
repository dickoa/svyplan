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
