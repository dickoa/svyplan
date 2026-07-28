test_that("prec_cluster computes CV for 2-stage", {
  result <- prec_cluster(n = c(50, 12), icc = 0.05)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "cluster")

  cv_exp <- sqrt(1 * 1 / (50 * 12) * (1 + 0.05 * (12 - 1)))
  expect_equal(result$cv, cv_exp, tolerance = 1e-6)
  expect_true(is.na(result$se))
  expect_true(is.na(result$moe))
})

test_that("prec_cluster computes CV for 3-stage", {
  result <- prec_cluster(n = c(50, 12, 8), icc = c(0.01, 0.05))
  expect_s3_class(result, "svyplan_prec")

  cv_exp <- sqrt(1 / (50 * 12 * 8) *
    (1 * 0.01 * 12 * 8 + (1 - 0.01) * (1 + 0.05 * (8 - 1))))
  expect_equal(result$cv, cv_exp, tolerance = 1e-6)
})

test_that("prec_cluster with resp_rate deflates stage-1", {
  base <- prec_cluster(n = c(50, 12), icc = 0.05)
  rr <- prec_cluster(n = c(50, 12), icc = 0.05, resp_rate = 0.8)
  cv_exp <- sqrt(1 / (40 * 12) * (1 + 0.05 * 11))
  expect_equal(rr$cv, cv_exp, tolerance = 1e-6)
  expect_true(rr$cv > base$cv)
})

test_that("prec_cluster.svyplan_cluster carries stage_cost metadata", {
  s1 <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  expect_equal(unname(p1$params$stage_cost), c(500, 50))
  expect_equal(names(p1$params$stage_cost), c("cost_psu", "cost_ssu"))
  expect_equal(p1$params$budget, 100000)
})

test_that("prec_cluster validates inputs", {
  expect_error(prec_cluster(n = 10, icc = 0.05), "length >= 2")
  expect_error(prec_cluster(n = c(50, -1), icc = 0.05), "positive")
  expect_error(prec_cluster(n = c(50, 12), icc = c(0.05, 0.1)),
               "must have length")
})

test_that("prec_cluster rejects NA in n", {
  expect_error(prec_cluster(n = c(50, NA), icc = 0.05), "NA")
})

test_that("prec_cluster rejects infinite n", {
  expect_error(
    prec_cluster(n = c(Inf, 10), icc = 0.05),
    "finite"
  )
})

test_that("prec_cluster prints cluster format", {
  result <- prec_cluster(n = c(50, 12), icc = 0.05)
  out <- capture.output(print(result))
  expect_match(out[1], "Sampling precision for 2-stage cluster")
  expect_match(out[2], "n_psu = 50")
  expect_match(out[2], "n_per_psu = 12")
  expect_match(out[3], "cv =")
})

test_that("prec_cluster rejects negative unit_relvar", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, unit_relvar = -1),
    "unit_relvar.*positive"
  )
})

test_that("prec_cluster rejects zero unit_relvar", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, unit_relvar = 0),
    "unit_relvar.*positive"
  )
})

test_that("prec_cluster rejects negative var_ratio", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, var_ratio = -1),
    "var_ratio.*positive"
  )
})

test_that("prec_cluster rejects NA var_ratio", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, var_ratio = NA),
    "var_ratio.*positive"
  )
})

test_that("prec_cluster accepts named n vectors and reorders", {
  ref2 <- prec_cluster(n = c(50, 12), icc = 0.05)
  named2 <- prec_cluster(n = c(n_per_psu = 12, n_psu = 50), icc = 0.05)
  alias2 <- prec_cluster(n = c(n_per_psu = 12, n_psu = 50), icc = 0.05)
  expect_equal(named2$cv, ref2$cv, tolerance = 1e-10)
  expect_equal(alias2$cv, ref2$cv, tolerance = 1e-10)
  expect_equal(names(named2$params$n), c("n_psu", "n_per_psu"))

  ref3 <- prec_cluster(n = c(50, 12, 8), icc = c(0.01, 0.05))
  named3 <- prec_cluster(
    n = c(n_per_ssu = 8, n_psu = 50, n_per_psu = 12),
    icc = c(0.01, 0.05)
  )
  alias3 <- prec_cluster(
    n = c(n_per_ssu = 8, n_psu = 50, n_per_psu = 12),
    icc = c(0.01, 0.05)
  )
  expect_equal(named3$cv, ref3$cv, tolerance = 1e-10)
  expect_equal(alias3$cv, ref3$cv, tolerance = 1e-10)
  expect_equal(names(named3$params$n), c("n_psu", "n_per_psu", "n_per_ssu"))

  expect_error(
    prec_cluster(n = c(n1 = 50, n2 = 12), icc = 0.05),
    "unrecognized names"
  )
})
