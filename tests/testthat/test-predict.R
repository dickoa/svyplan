test_that("predict.svyplan_n returns correct shape for proportions", {
  x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
  nd <- expand.grid(deff = c(1, 1.5, 2), resp_rate = c(0.8, 0.9))
  res <- predict(x, nd)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 6L)
  expect_true(all(c("deff", "resp_rate", "n", "se", "moe", "cv") %in% names(res)))
})

test_that("higher deff produces larger n for proportions", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(deff = c(1, 2, 3))
  res <- predict(x, nd)

  expect_true(all(diff(res$n) > 0))
})

test_that("single-row identity for n_prop", {
  x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.9)
  nd <- data.frame(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.9)
  res <- predict(x, nd)

  expect_equal(res$n, x$n, tolerance = 1e-8)
  expect_equal(res$se, x$se, tolerance = 1e-8)
})

test_that("predict.svyplan_n switches from moe to cv mode", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(cv = c(0.05, 0.10))
  res <- predict(x, nd)

  ref1 <- n_prop(p = 0.3, cv = 0.05)
  ref2 <- n_prop(p = 0.3, cv = 0.10)
  expect_equal(res$n[1], ref1$n, tolerance = 1e-8)
  expect_equal(res$n[2], ref2$n, tolerance = 1e-8)
})

test_that("predict.svyplan_n works for means", {
  x <- n_mean(var = 100, moe = 2, deff = 1.5)
  nd <- data.frame(deff = c(1, 2))
  res <- predict(x, nd)

  expect_equal(nrow(res), 2L)
  expect_true(res$n[2] > res$n[1])
})

test_that("single-row identity for n_mean", {
  x <- n_mean(var = 100, mu = 50, cv = 0.05, N = 5000)
  nd <- data.frame(var = 100, mu = 50, cv = 0.05, N = 5000)
  res <- predict(x, nd)

  expect_equal(res$n, x$n, tolerance = 1e-8)
})

test_that("no duplicate columns when newdata overlaps result names", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(moe = c(0.03, 0.05))
  res <- predict(x, nd)
  expect_equal(sum(names(res) == "moe"), 1L)
  expect_equal(ncol(res), 4L)

  nd2 <- data.frame(cv = c(0.05, 0.10))
  res2 <- predict(x, nd2)
  expect_equal(sum(names(res2) == "cv"), 1L)

  cl <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  nd3 <- data.frame(cv = c(0.03, 0.05))
  res3 <- predict(cl, nd3)
  expect_equal(sum(names(res3) == "cv"), 1L)

  pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
  nd4 <- data.frame(n = c(200, 500))
  res4 <- predict(pw, nd4)
  expect_equal(sum(names(res4) == "n"), 1L)
})

test_that("predict.svyplan_n errors on both moe and cv in newdata", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(moe = 0.05, cv = 0.10)
  expect_error(predict(x, nd), "cannot contain both")
})

test_that("predict.svyplan_n errors for multi-indicator results", {
  targets <- data.frame(
    name = c("a", "b"),
    p    = c(0.3, 0.5),
    moe  = c(0.05, 0.05)
  )
  x <- n_multi(targets)
  expect_error(predict(x, data.frame(deff = 1)), "multi-indicator")
})

test_that("predict.svyplan_cluster varies budget", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  nd <- data.frame(budget = c(50000, 100000, 200000))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(all(c("n_psu", "psu_size", "total_n", "cv", "cost") %in% names(res)))
  expect_true(all(diff(res$cv) < 0))
})

test_that("predict.svyplan_cluster varies cv", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  nd <- data.frame(cv = c(0.03, 0.05, 0.10))
  res <- predict(x, nd)

  expect_true(all(diff(res$cost) < 0))
})

test_that("predict.svyplan_cluster uses original mode when no cv/budget in newdata", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  nd <- data.frame(resp_rate = c(0.8, 0.9, 1.0))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(res$n_psu[1] > res$n_psu[3])
})

test_that("predict.svyplan_cluster errors on both cv and budget", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  nd <- data.frame(cv = 0.05, budget = 100000)
  expect_error(predict(x, nd), "cannot contain both")
})

test_that("predict.svyplan_cluster rejects multi-indicator results", {
  targets <- data.frame(
    name   = c("a", "b"),
    p      = c(0.3, 0.1),
    cv     = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  x <- n_multi(targets, cost = c(500, 50))
  expect_error(predict(x, data.frame(cv = 0.05)), "multi-indicator")
})

test_that("predict.svyplan_power varies n (solved for power)", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
  nd <- data.frame(n = seq(200, 800, 200))
  res <- predict(pw, nd)

  expect_equal(nrow(res), 4L)
  expect_true(all(c("n", "power", "effect") %in% names(res)))
  expect_true(all(diff(res$power) > 0))
})

test_that("predict.svyplan_power varies power (solved for n)", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35)
  nd <- data.frame(power = c(0.70, 0.80, 0.90))
  res <- predict(pw, nd)

  expect_true(all(diff(res$n) > 0))
})

test_that("predict.svyplan_power works for means", {
  pw <- power_mean(effect = 5, var = 100, n = 200, power = NULL)
  nd <- data.frame(n = c(100, 200, 400))
  res <- predict(pw, nd)

  expect_true(all(diff(res$power) > 0))
})

test_that("predict.svyplan_power errors for solved-for param in newdata", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35)
  expect_error(predict(pw, data.frame(n = 500)), "unknown parameter")
})

test_that("predict.svyplan_prec varies n for proportions", {
  x <- prec_prop(p = 0.3, n = 400)
  nd <- data.frame(n = c(100, 400, 1600))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(all(c("se", "moe", "cv") %in% names(res)))
  expect_true(all(diff(res$se) < 0))
})

test_that("predict.svyplan_prec varies n for means", {
  x <- prec_mean(var = 100, n = 400, mu = 50)
  nd <- data.frame(n = c(100, 400, 1600))
  res <- predict(x, nd)

  expect_true(all(diff(res$se) < 0))
})

test_that("single-row identity for prec_prop", {
  x <- prec_prop(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.9)
  nd <- data.frame(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.9)
  res <- predict(x, nd)

  expect_equal(res$se, x$se, tolerance = 1e-8)
  expect_equal(res$moe, x$moe, tolerance = 1e-8)
})

test_that("predict.svyplan_prec errors for cluster type", {
  x <- prec_cluster(n = c(50, 12), delta = 0.05)
  expect_error(predict(x, data.frame(n = 100)), "not supported")
})

test_that("predict errors for non-data.frame newdata", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(predict(x, list(deff = 1)), "must be a data frame")
})

test_that("predict errors for empty newdata", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(predict(x, data.frame()), "at least one")
})

test_that("predict errors for unknown param names", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(predict(x, data.frame(bogus = 1)), "unknown parameter")
})

test_that("predict errors for non-numeric columns", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(
    predict(x, data.frame(deff = "high", stringsAsFactors = FALSE)),
    "non-numeric"
  )
})

test_that("predict produces NA + warning on row failure", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(p = c(0.3, 0, 0.5))
  expect_warning(
    res <- predict(x, nd),
    "predict row 2 failed"
  )
  expect_true(is.na(res$n[2]))
  expect_false(is.na(res$n[1]))
  expect_false(is.na(res$n[3]))
})

test_that("predict.svyplan_cluster supports fixed_cost in newdata", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05,
                  fixed_cost = 5000)
  nd <- data.frame(fixed_cost = c(0, 5000, 10000))
  res <- predict(x, nd)
  expect_equal(nrow(res), 3L)
  expect_true("fixed_cost" %in% names(res))
  expect_true(res$cost[3] > res$cost[1])
  expect_equal(res$n_psu[1], res$n_psu[2], tolerance = 1e-8)
})
