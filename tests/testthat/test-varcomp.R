test_that("varcomp 2-stage SRS formula interface works", {
  set.seed(42)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  result <- varcomp(income ~ district, data = frame)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)
  expect_true(result$delta >= 0 && result$delta <= 1)
  expect_true(result$k > 0)
  expect_true(result$rel_var > 0)
})

test_that("varcomp 2-stage SRS vector interface matches formula", {
  set.seed(42)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  result_formula <- varcomp(income ~ district, data = frame)
  result_vector <- varcomp(frame$income, stage_id = list(frame$district))

  expect_equal(result_vector$var_between, result_formula$var_between, tolerance = 1e-10)
  expect_equal(result_vector$var_within, result_formula$var_within, tolerance = 1e-10)
  expect_equal(result_vector$delta, result_formula$delta, tolerance = 1e-10)
})

test_that("varcomp 2-stage SRS known values", {
  # Manually construct a known case
  set.seed(123)
  psu_id <- rep(1:10, each = 5)
  y <- rnorm(50, 100, 20)

  result <- varcomp(y, stage_id = list(psu_id))

  # Replicate ANOVA decomposition
  M <- 10L
  Ni <- as.numeric(table(psu_id))
  ti <- as.numeric(by(y, INDICES = psu_id, FUN = sum))
  S2Ui <- as.numeric(by(y, INDICES = psu_id, FUN = var))
  tbarU <- mean(ti)
  tU <- M * tbarU
  S2U1 <- var(ti)
  B2 <- S2U1 / tbarU^2
  W2 <- M * sum(Ni^2 * S2Ui) / tU^2

  expect_equal(result$var_between, B2, tolerance = 1e-10)
  expect_equal(result$var_within, W2, tolerance = 1e-10)
  expect_equal(result$delta, B2 / (B2 + W2), tolerance = 1e-10)
})

test_that("varcomp 2-stage PPS known values", {
  set.seed(456)
  psu_id <- rep(1:10, each = 5)
  y <- rnorm(50, 100, 20)
  pp <- rep(1 / 10, 10)

  result <- varcomp(y, stage_id = list(psu_id), prob = pp)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)

  # Replicate PPS decomposition
  Ni <- as.numeric(table(psu_id))
  cl_tots <- as.numeric(by(y, INDICES = psu_id, FUN = sum))
  cl_vars <- as.numeric(by(y, INDICES = psu_id, FUN = var))
  tU <- sum(cl_tots)
  S2U1 <- sum(pp * (cl_tots / pp - tU)^2)
  B2 <- S2U1 / tU^2
  W2 <- sum(Ni^2 * cl_vars / pp) / tU^2

  expect_equal(result$var_between, B2, tolerance = 1e-10)
  expect_equal(result$var_within, W2, tolerance = 1e-10)
})

test_that("varcomp 2-stage PPS with formula prob", {
  set.seed(789)
  frame <- data.frame(
    income = rnorm(50, 50000, 10000),
    district = rep(1:10, each = 5),
    pp = rep(1 / 10, 50)
  )
  result <- varcomp(income ~ district, data = frame, prob = ~pp)
  expect_s3_class(result, "svyplan_varcomp")
})

test_that("varcomp integrates with n_cluster", {
  set.seed(42)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  vc <- varcomp(income ~ district, data = frame)
  plan <- n_cluster(cost = c(500, 50), delta = vc, budget = 100000)
  expect_s3_class(plan, "svyplan_cluster")
})

test_that("varcomp validates inputs", {
  expect_error(varcomp("not_numeric"), "must be a formula or a numeric")
  expect_error(varcomp(1:10), "'stage_id' must be a list")
  expect_error(varcomp(~ x, data = data.frame(x = 1:10)),
               "must have a response")
  expect_error(
    varcomp(income ~ district, data = data.frame(x = 1:10)),
    "not found"
  )
})

test_that("varcomp handles lonely SSUs", {
  # One cluster with a single element
  psu_id <- c(1, 2, 2, 3, 3, 3)
  y <- c(10, 20, 25, 30, 35, 40)
  result <- varcomp(y, stage_id = list(psu_id))
  expect_s3_class(result, "svyplan_varcomp")
  expect_false(is.na(result$delta))
})
