test_that("prec_cluster computes CV for 2-stage", {
  result <- prec_cluster(n = c(50, 12), delta = 0.05)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "cluster")

  cv_exp <- sqrt(1 * 1 / (50 * 12) * (1 + 0.05 * (12 - 1)))
  expect_equal(result$cv, cv_exp, tolerance = 1e-6)
  expect_true(is.na(result$se))
  expect_true(is.na(result$moe))
})

test_that("prec_cluster computes CV for 3-stage", {
  result <- prec_cluster(n = c(50, 12, 8), delta = c(0.01, 0.05))
  expect_s3_class(result, "svyplan_prec")

  cv_exp <- sqrt(1 / (50 * 12 * 8) *
    (1 * 0.01 * 12 * 8 + 1 * (1 + 0.05 * (8 - 1))))
  expect_equal(result$cv, cv_exp, tolerance = 1e-6)
})

test_that("prec_cluster with resp_rate deflates stage-1", {
  base <- prec_cluster(n = c(50, 12), delta = 0.05)
  rr <- prec_cluster(n = c(50, 12), delta = 0.05, resp_rate = 0.8)
  cv_exp <- sqrt(1 / (40 * 12) * (1 + 0.05 * 11))
  expect_equal(rr$cv, cv_exp, tolerance = 1e-6)
  expect_true(rr$cv > base$cv)
})

test_that("prec_cluster.svyplan_cluster carries cost metadata", {
  s1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  expect_equal(p1$params$cost, c(500, 50))
  expect_equal(p1$params$budget, 100000)
})

test_that("prec_cluster validates inputs", {
  expect_error(prec_cluster(n = 10, delta = 0.05), "length >= 2")
  expect_error(prec_cluster(n = c(50, -1), delta = 0.05), "positive")
  expect_error(prec_cluster(n = c(50, 12), delta = c(0.05, 0.1)),
               "must have length")
})

test_that("prec_cluster prints cluster format", {
  result <- prec_cluster(n = c(50, 12), delta = 0.05)
  out <- capture.output(print(result))
  expect_match(out[1], "Sampling precision for 2-stage cluster")
  expect_match(out[2], "n1 = 50")
  expect_match(out[2], "n2 = 12")
  expect_match(out[3], "cv =")
})
