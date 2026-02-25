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

  expect_equal(
    result_vector$var_between,
    result_formula$var_between,
    tolerance = 1e-10
  )
  expect_equal(
    result_vector$var_within,
    result_formula$var_within,
    tolerance = 1e-10
  )
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
  expect_error(varcomp("not_numeric"), "must be a formula")
  expect_error(varcomp(1:10), "'stage_id' must be a list")
  expect_error(varcomp(~x, data = data.frame(x = 1:10)), "must have a response")
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

test_that("varcomp.survey.design matches formula interface", {
  skip_if_not_installed("survey")
  set.seed(42)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  ref <- varcomp(income ~ district, data = frame)

  dsgn <- survey::svydesign(
    ids = ~district,
    data = frame,
    weights = rep(1, 200)
  )
  result <- varcomp(dsgn, ~income)

  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$var_between, ref$var_between, tolerance = 1e-10)
  expect_equal(result$var_within, ref$var_within, tolerance = 1e-10)
  expect_equal(result$delta, ref$delta, tolerance = 1e-10)
})

test_that("varcomp.survey.design works with strata", {
  skip_if_not_installed("survey")
  set.seed(1)
  frame <- data.frame(
    y = rnorm(200, 50, 10),
    cluster = rep(1:20, each = 10),
    stratum = rep(1:4, each = 50)
  )
  dsgn <- survey::svydesign(
    ids = ~cluster,
    strata = ~stratum,
    data = frame,
    weights = rep(1, 200),
    nest = TRUE
  )
  result <- varcomp(dsgn, ~y)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)
  expect_true(result$delta >= 0 && result$delta <= 1)
})

test_that("varcomp.survey.design feeds into n_cluster", {
  skip_if_not_installed("survey")
  set.seed(2)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  dsgn <- survey::svydesign(
    ids = ~district,
    data = frame,
    weights = rep(1, 200)
  )
  vc <- varcomp(dsgn, ~income)
  plan <- n_cluster(cost = c(500, 50), delta = vc, budget = 100000)
  expect_s3_class(plan, "svyplan_cluster")
})

test_that("3-stage varcomp with non-nested SSU IDs matches nested", {
  set.seed(3)
  psu_id <- rep(1:5, each = 10)
  ssu_id_non_nested <- rep(rep(1:2, each = 5), 5)
  ssu_id_nested <- interaction(psu_id, ssu_id_non_nested, drop = TRUE)
  y <- rnorm(50, 100, 20)
  pp <- rep(0.2, 5)

  res_non <- varcomp(y, stage_id = list(psu_id, ssu_id_non_nested), prob = pp)
  res_nested <- varcomp(y, stage_id = list(psu_id, ssu_id_nested), prob = pp)

  expect_equal(res_non$delta, res_nested$delta, tolerance = 1e-10)
  expect_equal(res_non$var_between, res_nested$var_between, tolerance = 1e-10)
  expect_equal(res_non$var_within, res_nested$var_within, tolerance = 1e-10)
})

test_that("all-singleton 2-stage returns delta = 1 with warning", {
  psu_id <- 1:10
  y <- rnorm(10, 100, 20)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "all clusters are singletons"
  )
  expect_equal(res$delta, 1)
  expect_false(is.nan(res$delta))
  expect_true(is.finite(res$k))
})

test_that("varcomp.survey.design validates inputs", {
  skip_if_not_installed("survey")
  frame <- data.frame(y = rnorm(50), cl = rep(1:10, each = 5))
  dsgn <- survey::svydesign(ids = ~cl, data = frame, weights = rep(1, 50))

  expect_error(varcomp(dsgn), "one-sided formula")
  expect_error(varcomp(dsgn, ~nonexistent), "not found")

  dsgn_nocl <- survey::svydesign(ids = ~1, data = frame, weights = rep(1, 50))
  expect_error(varcomp(dsgn_nocl, ~y), "no clusters")
})
