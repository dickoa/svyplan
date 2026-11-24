test_that("simple and cluster multi-indicator APIs have invariant classes", {
  simple_targets <- data.frame(p = 0.30, moe = 0.05)
  cluster_targets <- data.frame(
    p = 0.30,
    cv = 0.10,
    icc_1 = 0.05
  )

  simple <- n_multi(simple_targets)
  cluster <- n_cluster(
    cluster_targets,
    stage_cost = c(500, 50)
  )

  expect_s3_class(simple, "svyplan_n")
  expect_s3_class(cluster, "svyplan_cluster")
  expect_error(
    n_multi(cluster_targets, stage_cost = c(500, 50)),
    "belongs to n_cluster"
  )
})

test_that("simple API identity is not changed by cluster plan defaults", {
  plan <- svyplan(stage_cost = c(500, 50))
  result <- n_multi(data.frame(p = 0.30, moe = 0.05), plan = plan)

  expect_s3_class(result, "svyplan_n")
})

test_that("cluster precision does not require costs", {
  targets <- data.frame(
    p = 0.30,
    n_1 = 60,
    take_2 = 12,
    icc_1 = 0.05
  )

  precision <- prec_cluster(targets)

  expect_s3_class(precision, "svyplan_prec")
  expect_identical(precision$type, "cluster")
  expect_null(precision$params$stage_cost)
  expect_error(n_cluster(precision), "'stage_cost' is required")
  expect_s3_class(
    n_cluster(prec_cluster(targets, stage_cost = c(500, 50))),
    "svyplan_cluster"
  )
})

test_that("simple and cluster round trips cannot be mixed", {
  cluster <- n_cluster(
    data.frame(p = 0.30, cv = 0.10, icc_1 = 0.05),
    stage_cost = c(500, 50)
  )
  precision <- prec_cluster(cluster)

  expect_error(prec_multi(cluster), "prec_cluster")
  expect_error(n_multi(precision), "n_cluster")
})

## A one-domain joint budget fit, and non-Wald round trips

test_that("a one-domain joint budget fit holds min_n_domain", {
  tab <- data.frame(name = "y", p = 0.3, cv = 0.1, icc_1 = 0.05, region = "A",
                    N = 1e6, N_1 = 100)
  fit <- function(budget, ...) {
    n_cluster(tab, stage_cost = c(100, 10), domains = "region",
              allocation = "joint", budget = budget, fpc = "stage", ...)
  }
  expect_error(fit(2450, stage_n = c(2, NA), min_n_domain = 500),
               "min_n_domain = 500 not achievable")
  held <- expect_no_warning(fit(8000, stage_n = c(2, NA), min_n_domain = 500))
  expect_gte(held$domains$.total_n, 500 * (1 - 1e-8))
  free <- fit(2450)
  expect_equal(free$domains$.cost, 2450, tolerance = 1e-6)
})

test_that("a non-Wald precision result sizes again in moe mode", {
  fit <- n_multi(data.frame(name = "a", p = 0.3, cv = 0.05))
  for (m in c("wilson", "logodds", "beta")) {
    pr <- prec_multi(fit, prop_method = m)
    again <- n_multi(pr)
    expect_equal(again$n, fit$n, tolerance = 1e-8)
    expect_equal(prec_multi(again)$detail$.moe, pr$detail$.moe,
                 tolerance = 1e-8)
  }
})

test_that("a field size rounded up to the census does not warn", {
  expect_no_warning(x <- n_multi(data.frame(name = "a", p = 0.5, moe = 0.005,
                                            N = 100)))
  expect_lt(x$n, 100)
  expect_equal(x$detail$.cv_field, 0)
})

test_that("indicator table messages name the mean column", {
  expect_error(n_multi(data.frame(name = "a", var = 4, rmoe = 0.1)),
               "'mean' for a mean")
})
