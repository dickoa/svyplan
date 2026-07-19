test_that("simple and cluster multi-indicator APIs have invariant classes", {
  simple_targets <- data.frame(p = 0.30, moe = 0.05)
  cluster_targets <- data.frame(
    p = 0.30,
    cv = 0.10,
    delta_psu = 0.05
  )

  simple <- n_multi(simple_targets)
  cluster <- n_multi_cluster(
    cluster_targets,
    stage_cost = c(500, 50)
  )

  expect_s3_class(simple, "svyplan_n")
  expect_s3_class(cluster, "svyplan_cluster")
  expect_error(
    n_multi(cluster_targets, stage_cost = c(500, 50)),
    "moved to n_multi_cluster"
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
    n = 60,
    psu_size = 12,
    delta_psu = 0.05
  )

  precision <- prec_multi_cluster(targets)

  expect_s3_class(precision, "svyplan_prec")
  expect_identical(precision$params$design, "cluster")
  expect_null(precision$params$stage_cost)
  expect_error(n_multi_cluster(precision), "'stage_cost' is required")
  expect_s3_class(
    n_multi_cluster(precision, stage_cost = c(500, 50)),
    "svyplan_cluster"
  )
})

test_that("simple and cluster round trips cannot be mixed", {
  cluster <- n_multi_cluster(
    data.frame(p = 0.30, cv = 0.10, delta_psu = 0.05),
    stage_cost = c(500, 50)
  )
  precision <- prec_multi_cluster(cluster)

  expect_error(prec_multi(cluster), "prec_multi_cluster")
  expect_error(n_multi(precision), "n_multi_cluster")
})
