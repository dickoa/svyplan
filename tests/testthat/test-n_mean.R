test_that("n_mean MOE mode", {
  result <- n_mean(var = 100, moe = 2)
  z <- qnorm(0.975)
  expected <- z^2 * 100 / (4 + z^2 * 100 / Inf)
  expect_equal(result$n, expected, tolerance = 1e-6)
  expect_s3_class(result, "svyplan_n")
  expect_equal(result$type, "mean")
})

test_that("n_mean CV mode", {
  result <- n_mean(var = 100, mu = 50, cv = 0.05)
  CVpop <- sqrt(100) / 50  # 0.2
  expected <- CVpop^2 / (0.05^2 + CVpop^2 / Inf)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_mean MOE with FPC", {
  result <- n_mean(var = 100, moe = 2, N = 5000)
  z <- qnorm(0.975)
  expected <- z^2 * 100 / (4 + z^2 * 100 / 5000)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_mean CV with FPC", {
  result <- n_mean(var = 100, mu = 50, cv = 0.05, N = 5000)
  CVpop <- sqrt(100) / 50
  expected <- CVpop^2 / (0.05^2 + CVpop^2 / 5000)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_mean deff multiplier works", {
  base <- n_mean(var = 100, moe = 2)
  with_deff <- n_mean(var = 100, moe = 2, deff = 1.5)
  expect_equal(with_deff$n, base$n * 1.5, tolerance = 1e-6)
})

test_that("n_mean returns correct S3 class", {
  result <- n_mean(var = 100, moe = 2)
  expect_true(inherits(result, "svyplan_n"))
  expect_equal(as.integer(result), ceiling(result$n))
  expect_equal(as.double(result), result$n)
})

test_that("n_mean validates inputs", {
  expect_error(n_mean(var = -1, moe = 2), "positive")
  expect_error(n_mean(var = 100), "specify exactly one")
  expect_error(n_mean(var = 100, moe = 2, cv = 0.05), "specify exactly one")
  expect_error(n_mean(var = 100, cv = 0.05), "'mu' is required")
  expect_s3_class(n_mean(var = 100, moe = 2, deff = 0.5), "svyplan_n")
  expect_error(n_mean(var = 100, moe = 2, deff = 0), "must be positive")
  expect_error(n_mean(var = 100, moe = 2, deff = -1), "must be positive")
})

test_that("svyplan_n has se/moe/cv fields for mean", {
  res <- n_mean(var = 100, moe = 2, mu = 50)
  expect_true(is.numeric(res$se))
  expect_true(res$se > 0)
  expect_true(res$moe > 0)
  expect_true(res$cv > 0)
})

test_that("mean moe round-trip: moe from n matches target", {
  res <- n_mean(var = 100, moe = 2)
  expect_equal(res$moe, 2, tolerance = 1e-6)
})

test_that("mean moe round-trip with finite N", {
  res <- n_mean(var = 100, moe = 2, N = 5000)
  expect_equal(res$moe, 2, tolerance = 1e-6)
})

test_that("mean cv is NA when mu not provided", {
  res <- n_mean(var = 100, moe = 2)
  expect_true(is.na(res$cv))
})

test_that("n_mean rejects N = 1", {
  expect_error(n_mean(var = 100, moe = 2, N = 1), "greater than 1")
})

test_that("deff multiplies the SRS variance at the actual net n (finite N)", {
  z <- qnorm(0.975)
  for (deff in c(0.8, 2)) {
    got <- n_mean(var = 100, moe = 2, N = 500, deff = deff)$n
    expected <- z^2 * deff * 100 / (2^2 + z^2 * deff * 100 / 500)
    expect_equal(got, expected, tolerance = 1e-10)
    pm <- prec_mean(var = 100, n = got, N = 500, deff = deff)
    expect_equal(pm$moe, 2, tolerance = 1e-8)
  }
})

test_that("a CV target is planned on the magnitude of the mean", {
  # the formula uses mu^2, so the sign cannot matter
  expect_equal(n_mean(sd = 10, mu = -20, cv = 0.1)$n,
               n_mean(sd = 10, mu = 20, cv = 0.1)$n)
  expect_equal(prec_mean(sd = 10, mu = -20, n = 100)$cv,
               prec_mean(sd = 10, mu = 20, n = 100)$cv)
  expect_gt(prec_mean(sd = 10, mu = -20, n = 100)$cv, 0)

  # zero is the value with no relative scale
  expect_error(n_mean(sd = 10, mu = 0, cv = 0.1), "'cv' is undefined at 'mu' = 0")
  expect_true(is.na(prec_mean(sd = 10, mu = 0, n = 100)$cv))

  # and the frame paths agree with the scalar ones
  frame <- data.frame(stratum = c("A", "B"), N = c(100, 100), sd = c(1, 2),
                      mean = c(-10, -20))
  flipped <- transform(frame, mean = c(10, 20))
  expect_equal(n_alloc(frame, n = 100, alloc = "neyman")$cv,
               n_alloc(flipped, n = 100, alloc = "neyman")$cv)

  two_phase <- data.frame(stratum = c("A", "B"), N = c(600, 400),
                          sd = c(10, 20), mean = c(-50, -70), unit_cost = 1)
  expect_equal(
    n_twophase(two_phase, phase1_cost = 1, budget = 5000)$cv,
    n_twophase(transform(two_phase, mean = c(50, 70)),
               phase1_cost = 1, budget = 5000)$cv
  )
})

test_that("a tiny target gives the census rather than a size that misses it", {
  for (t in 10^seq(-6, -10, by = -0.25)) {
    n <- suppressWarnings(n_mean(var = 4, moe = t, N = 1000, deff = 2))$n
    got <- suppressWarnings(prec_mean(var = 4, n = n, N = 1000, deff = 2))$moe
    expect_lte(got, t * (1 + 1e-8), label = t)
  }
  fit <- n_mean(var = 4, moe = 1e-3, N = 1000, deff = 2)
  expect_lt(fit$n, 1000)
  expect_equal(prec_mean(var = 4, n = fit$n, N = 1000, deff = 2)$moe, 1e-3,
               tolerance = 1e-8)
})

## A zero mean serves an absolute target

test_that("a zero mean is accepted with an absolute target, its relative measures NA", {
  # Closed form: the size and se do not involve mu at all.
  plain <- n_mean(var = 4, moe = 1)
  zero <- n_mean(var = 4, mu = 0, moe = 1)
  expect_equal(zero$n, plain$n)
  expect_true(is.na(zero$cv))
  expect_true(is.na(zero$rmoe))
  p <- prec_mean(var = 4, mu = 0, n = 100)
  expect_equal(p$se, sqrt(4 / 100))
  expect_true(is.na(p$cv))
  expect_true(is.na(p$rmoe))
  expect_true(is.na(prec_mean(zero)$cv))
  expect_error(n_mean(var = 4, mu = 0, cv = 0.1),
               "'cv' is undefined at 'mu' = 0")
  expect_error(n_mean(var = 4, mu = 0, rmoe = 0.1),
               "'rmoe' is undefined at 'mu' = 0")
})

test_that("pooled and change results give NA, not Inf, at a zero level", {
  checked <- character(0)
  for (case in list(
    list("pooled", prec_pooled(var = 4, n = 100, occasions = 2, mu = 0)),
    list("pooled_n", n_pooled(var = 4, moe = 1, occasions = 2, mu = 0)),
    list("change", prec_change(var = 4, n = 100, change = 0))
  )) {
    expect_true(is.na(case[[2]]$cv), info = case[[1]])
    expect_true(is.na(case[[2]]$rmoe), info = case[[1]])
    checked <- c(checked, case[[1]])
  }
  expect_setequal(checked, c("pooled", "pooled_n", "change"))
})

test_that("n_multi takes a zero mean on a moe row and refuses it on a cv row", {
  m <- n_multi(data.frame(name = c("a", "b"), var = 4, mean = c(0, 5),
                          moe = 1))
  expect_equal(m$detail$.n[1], n_mean(var = 4, moe = 1)$n)
  expect_true(is.na(m$detail$.cv_achieved[1]))
  expect_error(n_multi(data.frame(var = 4, mean = 0, cv = 0.1)),
               "must be non-zero on a row with a 'cv' target")
  expect_error(
    n_cluster(data.frame(var = 4, mean = 0, moe = 1,
                                      icc_1 = 0.05),
              stage_cost = c(500, 50)),
    "multistage mode"
  )
})
