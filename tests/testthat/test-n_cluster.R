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

test_that("n_cluster round-trips with prec_cluster", {
  plan <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  cv_check <- prec_cluster(n = unname(plan$n), delta = 0.05)$cv
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

test_that("n_cluster is an S3 generic", {
  expect_true(is.function(n_cluster))
  # UseMethod body check
  expect_true(grepl("UseMethod", deparse(body(n_cluster))[1]))
  # Verify dispatch works on numeric (default method)
  res <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_s3_class(res, "svyplan_cluster")
})

test_that("n_cluster params store cv in CV mode", {
  res <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  expect_equal(res$params$cv, 0.05)
  expect_null(res$params$budget)
  expect_null(res$params$m)
})

test_that("n_cluster params store budget and m", {
  res <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000, m = 40)
  expect_equal(res$params$budget, 100000)
  expect_equal(res$params$m, 40)
  expect_null(res$params$cv)
})

test_that("n_cluster params store budget without m", {
  res <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_equal(res$params$budget, 100000)
  expect_null(res$params$m)
})

test_that("n_cluster 3-stage params store cv", {
  res <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05), cv = 0.05)
  expect_equal(res$params$cv, 0.05)
  expect_null(res$params$budget)
})

test_that("n_cluster 3-stage params store budget and m", {
  res <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                   budget = 500000, m = 50)
  expect_equal(res$params$budget, 500000)
  expect_equal(res$params$m, 50)
})

test_that("svyplan_cluster has se/moe/cv fields", {
  res <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_true("se" %in% names(res))
  expect_true("moe" %in% names(res))
  expect_true("cv" %in% names(res))
  expect_true(is.na(res$se))
  expect_true(is.na(res$moe))
  expect_true(res$cv > 0)
})

test_that("n_cluster rejects invalid rel_var", {
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, rel_var = -1, cv = 0.05),
    "must be positive"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, rel_var = 0, cv = 0.05),
    "must be positive"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, rel_var = NA_real_, cv = 0.05),
    "must not be NA"
  )
})

test_that("n_cluster rejects invalid k", {
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, k = -1, cv = 0.05),
    "positive"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, k = 0, cv = 0.05),
    "positive"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, k = NA, cv = 0.05),
    "positive"
  )
  expect_error(
    n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
              k = c(1, -1), cv = 0.05),
    "positive"
  )
})

test_that("n_cluster fixed-m CV error says 'too small'", {
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.001, m = 5),
    "too small"
  )
  expect_error(
    n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
              cv = 0.001, m = 5),
    "too small"
  )
})

test_that("fixed_cost = 0 is backward compatible", {
  r1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  r2 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000,
                   fixed_cost = 0)
  expect_equal(r1$n, r2$n)
  expect_equal(r1$cost, r2$cost)
  expect_null(r2$params$fixed_cost)
})

test_that("fixed_cost 2-stage budget mode reduces n1, n2 unchanged", {
  base <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  fc <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000,
                   fixed_cost = 5000)
  expect_equal(fc$n[["n2"]], base$n[["n2"]], tolerance = 1e-10)
  expect_true(fc$n[["n1"]] < base$n[["n1"]])
  expect_equal(fc$cost, 100000)
  expect_equal(fc$params$fixed_cost, 5000)
})

test_that("fixed_cost 2-stage CV mode leaves n unchanged, adds to cost", {
  base <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  fc <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05,
                   fixed_cost = 5000)
  expect_equal(fc$n, base$n, tolerance = 1e-10)
  expect_equal(fc$cost, base$cost + 5000, tolerance = 1e-6)
})

test_that("fixed_cost 3-stage budget mode reduces n1", {
  base <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                    budget = 500000)
  fc <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                   budget = 500000, fixed_cost = 10000)
  expect_equal(fc$n[["n2"]], base$n[["n2"]], tolerance = 1e-10)
  expect_equal(fc$n[["n3"]], base$n[["n3"]], tolerance = 1e-10)
  expect_true(fc$n[["n1"]] < base$n[["n1"]])
  expect_equal(fc$cost, 500000)
})

test_that("fixed_cost 3-stage CV mode adds to cost", {
  base <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                    cv = 0.05)
  fc <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05),
                   cv = 0.05, fixed_cost = 10000)
  expect_equal(fc$n, base$n, tolerance = 1e-10)
  expect_equal(fc$cost, base$cost + 10000, tolerance = 1e-6)
})

test_that("fixed_cost validation rejects bad inputs", {
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05, fixed_cost = -1),
    "non-negative"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05, fixed_cost = NA),
    "non-negative"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05,
              fixed_cost = c(1, 2)),
    "non-negative"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000,
              fixed_cost = 100000),
    "less than"
  )
  expect_error(
    n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000,
              fixed_cost = 200000),
    "less than"
  )
})

test_that("fixed_cost round-trip n_cluster -> prec_cluster -> n_cluster", {
  orig <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05,
                    fixed_cost = 5000)
  prec <- prec_cluster(orig)
  expect_equal(prec$params$fixed_cost, 5000)
  back <- n_cluster(prec)
  expect_equal(unname(back$n), unname(orig$n), tolerance = 1e-4)
  expect_equal(back$params$fixed_cost, 5000)
})
