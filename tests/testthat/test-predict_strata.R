test_that("predict returns factor with correct levels", {
  set.seed(42)
  x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 4, n = 100)
  f <- predict(sb, x)
  expect_s3_class(f, "factor")
  expect_equal(nlevels(f), 4L)
  expect_equal(length(f), length(x))
  expect_false(anyNA(f))
})

test_that("predict assigns all observations", {
  set.seed(42)
  x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 4, n = 100)
  f <- predict(sb, x)
  expect_equal(as.integer(table(f)), sb$strata$N_h)
})

test_that("predict accepts custom labels", {
  set.seed(42)
  x <- rlnorm(200, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 3, n = 60)
  labs <- c("Low", "Mid", "High")
  f <- predict(sb, x, labels = labs)
  expect_equal(levels(f), labs)
})

test_that("predict returns NA for out-of-range values", {
  set.seed(42)
  x <- rlnorm(200, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 3, n = 60)
  newx <- c(-1, x, max(x) + 1e6)
  f <- predict(sb, newx)
  expect_true(is.na(f[1L]))
  expect_true(is.na(f[length(f)]))
  expect_false(anyNA(f[2L:(length(f) - 1L)]))
})

test_that("predict works with new data of different length", {
  set.seed(42)
  x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 4, n = 100)
  newx <- rlnorm(50, meanlog = 6, sdlog = 1.5)
  f <- predict(sb, newx)
  expect_equal(length(f), 50L)
  expect_equal(nlevels(f), 4L)
})

test_that("predict works with certain stratum", {
  set.seed(42)
  x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 3, n = 80, certain = quantile(x, 0.95))
  f <- predict(sb, x)
  expect_equal(nlevels(f), 3L)
  expect_equal(length(f), length(x))
})

test_that("predict validates inputs", {
  set.seed(42)
  x <- rlnorm(200, meanlog = 6, sdlog = 1.5)
  sb <- strata_bound(x, n_strata = 3, n = 60)
  expect_error(predict(sb, "text"), "numeric vector")
  expect_error(predict(sb, x, labels = c("A", "B")), "length 3")
})
