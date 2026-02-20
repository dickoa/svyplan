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

test_that("design_effect henry returns 1 for equal weights", {
  set.seed(42)
  n <- 100
  w <- rep(5, n)
  x_cal <- runif(n, 10, 100)
  y <- 5 + 0.3 * x_cal + rnorm(n, 0, 2)
  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_equal(result, 1.0)
})

test_that("design_effect spencer returns 1 for equal weights", {
  set.seed(42)
  n <- 100
  w <- rep(5, n)
  p <- rep(0.1, n)
  y <- 50 + rnorm(n, 0, 10)
  result <- design_effect(w, y = y, p = p, method = "spencer")
  expect_equal(result, 1.0)
})

test_that("design_effect henry returns 1 for constant y", {
  set.seed(42)
  n <- 100
  w <- runif(n, 1, 5)
  x_cal <- runif(n, 10, 100)
  y <- rep(50, n)
  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_equal(result, 1.0)
})

test_that("design_effect spencer returns 1 for constant y", {
  set.seed(42)
  n <- 100
  w <- runif(n, 1, 5)
  p <- runif(n, 0.01, 0.2)
  y <- rep(50, n)
  result <- design_effect(w, y = y, p = p, method = "spencer")
  expect_equal(result, 1.0)
})

test_that("design_effect spencer handles equal probabilities", {
  set.seed(42)
  n <- 100
  w <- runif(n, 1, 5)
  p <- rep(0.1, n)
  y <- 50 + rnorm(n, 0, 10)
  result <- design_effect(w, y = y, p = p, method = "spencer")
  expect_true(is.numeric(result))
  expect_false(is.nan(result))
})

test_that("design_effect cluster rejects invalid delta", {
  expect_error(
    design_effect(delta = -0.2, m = 10, method = "cluster"),
    "\\[0, 1\\]"
  )
  expect_error(
    design_effect(delta = 1.5, m = 10, method = "cluster"),
    "\\[0, 1\\]"
  )
})

test_that("design_effect rejects zero weights", {
  expect_error(
    design_effect(c(1, 0, 2), method = "kish"),
    "positive"
  )
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

test_that("design_effect cr stratified + clustered", {
  set.seed(42)
  n <- 200
  strvar <- rep(1:5, each = 40)
  clvar <- rep(1:50, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- rep(2L, 5)

  result <- design_effect(w, y = y, strvar = strvar, clvar = clvar,
                          stages = stages, method = "cr")
  expect_true(is.list(result))
  expect_true("strata" %in% names(result))
  expect_true("overall" %in% names(result))
  expect_true(is.numeric(result$overall))
  expect_equal(nrow(result$strata), 5L)
  expect_true(all(c("deff_w", "deff_c", "deff_s") %in% names(result$strata)))
  expect_equal(result$overall, sum(result$strata$deff_w *
                                     result$strata$deff_c *
                                     result$strata$deff_s))
})

test_that("design_effect cr stratified, no clusters", {
  set.seed(42)
  n <- 100
  strvar <- rep(1:5, each = 20)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, strvar = strvar, method = "cr")
  expect_true(is.list(result))
  expect_equal(result$overall, sum(result$strata$deff_w * result$strata$deff_s))
})

test_that("design_effect cr unstratified + clustered", {
  set.seed(42)
  n <- 100
  clvar <- rep(1:25, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, clvar = clvar, method = "cr")
  expect_true(is.list(result))
  expect_true("rho" %in% names(result$strata))
  expect_equal(result$overall, result$strata$deff_w * result$strata$deff_c)
})

test_that("design_effect cr mixed stages", {
  set.seed(42)
  n <- 160
  strvar <- rep(1:4, each = 40)
  clvar <- rep(1:40, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- c(1L, 2L, 2L, 1L)

  result <- design_effect(w, y = y, strvar = strvar, clvar = clvar,
                          stages = stages, method = "cr")
  expect_equal(result$strata$deff_c[1], 1)
  expect_equal(result$strata$deff_c[4], 1)
})

test_that("design_effect cr equal weights gives deff_w near 1", {
  set.seed(42)
  n <- 100
  clvar <- rep(1:25, each = 4)
  w <- rep(50, n)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, clvar = clvar, method = "cr")
  expect_equal(result$strata$deff_w, 1, tolerance = 1e-10)
})

test_that("design_effect cr agrees with survey package", {
  skip_if_not_installed("survey")
  set.seed(42)
  n <- 200
  strvar <- rep(1:4, each = 50)
  clvar <- rep(1:50, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- rep(2L, 4)

  our <- design_effect(w, y = y, strvar = strvar, clvar = clvar,
                       stages = stages, method = "cr")

  dsgn <- survey::svydesign(ids = ~clvar, strata = ~strvar,
                            data = data.frame(y = y), weights = w,
                            nest = TRUE)
  mn <- survey::svymean(~y, design = dsgn, deff = TRUE)
  survey_deff <- as.numeric(survey::deff(mn))

  expect_equal(our$overall, survey_deff, tolerance = 0.05)
})
