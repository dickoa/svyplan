## Relative margin of error (rmoe)

## The scalar equivalences

test_that("n_prop(rmoe) equals n_prop(moe = rmoe * p) under all four methods", {
  p <- 0.2
  r <- 0.12
  for (m in c("wald", "wilson", "logodds", "beta")) {
    by_rmoe <- n_prop(p = p, rmoe = r, method = m, deff = 1.5,
                      resp_rate = 0.9, df = 25)
    by_moe <- n_prop(p = p, moe = r * p, method = m, deff = 1.5,
                     resp_rate = 0.9, df = 25)
    expect_identical(by_rmoe$n, by_moe$n)
    expect_equal(by_rmoe$rmoe, r, tolerance = 1e-10)
  }
})

test_that("n_mean(rmoe) uses the magnitude of mu", {
  for (mu in c(50, -50)) {
    by_rmoe <- n_mean(var = 2500, mu = mu, rmoe = 0.05, N = 10000)
    by_moe <- n_mean(var = 2500, mu = mu, moe = 0.05 * abs(mu), N = 10000)
    expect_identical(by_rmoe$n, by_moe$n)
    expect_equal(by_rmoe$rmoe, 0.05, tolerance = 1e-10)
  }
  expect_identical(
    n_mean(var = 2500, mu = 50, rmoe = 0.05)$n,
    n_mean(var = 2500, mu = -50, rmoe = 0.05)$n
  )
})

test_that("rmoe = q * cv holds under wald alone", {
  # The engine's quantile, which df = 25 makes a t rather than a normal.
  q <- stats::qt(0.975, 25)
  p <- 0.02
  args <- list(p = p, n = 900, deff = 2, df = 25)
  gaps <- vapply(c("wald", "wilson", "logodds", "beta"), function(m) {
    res <- do.call(prec_prop, c(args, list(method = m)))
    res$rmoe / (q * res$cv) - 1
  }, numeric(1L))
  expect_equal(unname(gaps[["wald"]]), 0, tolerance = 1e-12)
  expect_true(all(gaps[c("wilson", "logodds", "beta")] > 0.04))
})

## Reported back in the units it was stated in

test_that("an rmoe target round trips through prec_prop", {
  for (m in c("wald", "wilson", "logodds", "beta")) {
    sized <- n_prop(p = 0.2, rmoe = 0.12, method = m)
    expect_equal(prec_prop(sized)$rmoe, 0.12, tolerance = 1e-9)
  }
})

test_that("params records the target as supplied", {
  expect_equal(n_prop(p = 0.2, rmoe = 0.12)$params$rmoe, 0.12)
  expect_null(n_prop(p = 0.2, rmoe = 0.12)$params$moe)
  expect_null(n_prop(p = 0.2, rmoe = 0.12)$params$cv)
  expect_null(n_prop(p = 0.2, moe = 0.05)$params$rmoe)
})

test_that("rmoe is NA when the estimand carries no scale", {
  expect_true(is.na(n_mean(var = 100, moe = 2)$rmoe))
  expect_false(is.na(n_mean(var = 100, mu = 50, moe = 2)$rmoe))
})

## prec_prop solves for the level from an rmoe target

test_that("the rmoe solve reduces to the cv solve under wald", {
  q <- stats::qnorm(0.975)
  by_cv <- prec_prop(n = 1500, cv = 0.10, N = 2e6)
  by_rmoe <- prec_prop(n = 1500, rmoe = q * 0.10, method = "wald", N = 2e6)
  expect_equal(by_rmoe$params$p, by_cv$params$p, tolerance = 1e-8)
  expect_identical(by_rmoe$solved, "p")
})

test_that("the rmoe solve is method-specific and meets its target", {
  q <- stats::qnorm(0.975)
  levels <- vapply(c("wald", "wilson", "logodds", "beta"), function(m) {
    res <- prec_prop(n = 1500, rmoe = q * 0.10, method = m, N = 2e6)
    expect_equal(res$rmoe, q * 0.10, tolerance = 1e-8)
    res$params$p
  }, numeric(1L))
  # The interval methods order as they do everywhere else: a wider interval
  # needs a larger proportion to reach the same relative half-width.
  expect_equal(unname(levels), sort(unname(levels)))
  expect_gt(levels[["beta"]], levels[["wald"]])
  expect_false(isTRUE(all.equal(levels[["wilson"]], levels[["wald"]])))
})

test_that("an rmoe below the method's floor is refused with the floor named", {
  expect_error(prec_prop(n = 1500, rmoe = 1e-5, method = "wilson"),
               "unattainable")
  expect_error(prec_prop(n = 1500, rmoe = 1e-5, method = "beta"),
               "never falls below")
  expect_error(prec_prop(n = 1500, rmoe = 1e-5, method = "logodds"),
               "unattainable")
  # Wald's relative half-width does go to zero, so the same target is met.
  expect_lt(prec_prop(n = 1500, rmoe = 1e-5, method = "wald")$rmoe, 1.1e-5)
})

test_that("a census has no proportion to solve for from rmoe", {
  expect_error(prec_prop(n = 500, rmoe = 0.1, N = 500), "census")
})

test_that("predict on an rmoe solve keeps the rmoe metric", {
  fitted <- prec_prop(n = 1500, rmoe = 0.20, method = "beta", N = 2e6)
  grid <- predict(fitted, expand.grid(n = c(500, 1500, 5000)))
  expect_true(all(c("p", "rmoe") %in% names(grid)))
  expect_equal(grid$rmoe, rep(0.20, 3L), tolerance = 1e-8)
  # Larger designs report the same relative half-width at a smaller level.
  expect_equal(grid$p, sort(grid$p, decreasing = TRUE))
  expect_error(predict(fitted, data.frame(cv = 0.1)),
               "unknown parameter")
})

## Indicator tables

test_that("an rmoe column is normalized to moe at ingestion", {
  by_rmoe <- data.frame(
    name = c("stunting", "anaemia"), p = c(0.25, 0.40),
    rmoe = c(0.12, 0.10), deff = 1.5
  )
  by_moe <- data.frame(
    name = c("stunting", "anaemia"), p = c(0.25, 0.40),
    moe = c(0.12 * 0.25, 0.10 * 0.40), deff = 1.5
  )
  expect_identical(n_multi(by_rmoe)$n, n_multi(by_moe)$n)
  expect_equal(n_multi(by_rmoe)$detail, n_multi(by_moe)$detail)
})

test_that("an rmoe column works in the cluster path and reads the row method", {
  by_rmoe <- data.frame(
    name = c("a", "b"), p = c(0.20, 0.35), rmoe = c(0.10, 0.12),
    icc_psu = 0.05, prop_method = "wilson"
  )
  by_moe <- data.frame(
    name = c("a", "b"), p = c(0.20, 0.35),
    moe = c(0.10 * 0.20, 0.12 * 0.35),
    icc_psu = 0.05, prop_method = "wilson"
  )
  expect_identical(
    n_cluster(indicators = by_rmoe, stage_cost = c(500, 50))$n,
    n_cluster(indicators = by_moe, stage_cost = c(500, 50))$n
  )
})

test_that("an rmoe row takes its scale from mu for a mean indicator", {
  by_rmoe <- data.frame(name = "spend", var = 2500, mu = -40, rmoe = 0.05)
  by_moe <- data.frame(name = "spend", var = 2500, mu = -40, moe = 0.05 * 40)
  expect_identical(n_multi(by_rmoe)$n, n_multi(by_moe)$n)
})

test_that("prec_multi reports .rmoe against each row's own estimand", {
  indicators <- data.frame(
    name = c("stunting", "spend"), p = c(0.25, NA), var = c(NA, 2500),
    mu = c(NA, -40), n = 800
  )
  res <- prec_multi(indicators)
  expect_true(".rmoe" %in% names(res$detail))
  expect_equal(res$detail$.rmoe, res$detail$.moe / c(0.25, 40),
               tolerance = 1e-12)
})

## Generalized allocation

test_that("an rmoe target allocates exactly as cv = rmoe / q", {
  frame <- data.frame(stratum = paste0("s", 1:6),
                      N = c(1200, 800, 2400, 1500, 900, 3000))
  measures <- data.frame(
    stratum = frame$stratum, name = "A",
    p = c(0.12, 0.20, 0.08, 0.30, 0.15, 0.22)
  )
  q <- stats::qnorm(0.975)
  r <- 0.06
  by_rmoe <- n_alloc(frame, measures = measures,
                     targets = data.frame(name = "A", rmoe = r))
  by_cv <- n_alloc(frame, measures = measures,
                   targets = data.frame(name = "A", cv = r / q))
  expect_equal(by_rmoe$detail$n, by_cv$detail$n, tolerance = 1e-10)
  expect_equal(by_rmoe$n, by_cv$n, tolerance = 1e-10)

  # Reported in the units the target was stated in, not in its cv form.
  expect_identical(by_rmoe$constraints$.metric, "rmoe")
  expect_equal(by_rmoe$constraints$.achieved, r, tolerance = 1e-8)
  expect_equal(by_rmoe$constraints$.ratio, 1, tolerance = 1e-8)
  expect_equal(by_rmoe$constraints$.achieved,
               q * by_cv$constraints$.achieved, tolerance = 1e-8)

  # Sensitivity is d(Vmax)/dr scaled by the multiplier: the cv branch over q^2.
  total <- sum(frame$N * measures$p)
  expect_equal(
    by_rmoe$constraints$.sensitivity,
    -2 * by_rmoe$constraints$.multiplier * r * total^2 / q^2,
    tolerance = 1e-6
  )
})

test_that("an rmoe allocation round trips through prec_alloc", {
  frame <- data.frame(stratum = paste0("s", 1:5),
                      N = c(1200, 800, 2400, 1500, 900))
  measures <- data.frame(stratum = frame$stratum, name = "A",
                         p = c(0.12, 0.20, 0.08, 0.30, 0.15))
  fitted <- n_alloc(frame, measures = measures,
                    targets = data.frame(name = "A", rmoe = 0.06))
  assessed <- prec_alloc(fitted)
  expect_identical(assessed$detail$.metric, "rmoe")
  refitted <- n_alloc(assessed)
  expect_identical(refitted$constraints$.metric, "rmoe")
  expect_equal(refitted$detail$n, fitted$detail$n, tolerance = 1e-8)
})

test_that("stratum and domain tables carry .rmoe beside .moe", {
  frame <- data.frame(stratum = paste0("s", 1:3), N = c(1000, 2000, 1500),
                      sd = c(5, 8, 6), mean = c(20, 30, 25))
  res <- n_alloc(frame, n = 600)
  expect_true(".rmoe" %in% names(res$detail))
  expect_equal(res$detail$.rmoe, res$detail$.moe / frame$mean,
               tolerance = 1e-12)

  # No mean on the frame leaves nothing for a relative quantity to measure
  # against, exactly as it leaves .cv NA.
  bare <- n_alloc(data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                             sd = c(5, 8)), n = 300)
  expect_true(all(is.na(bare$detail$.rmoe)))
  expect_true(all(is.na(bare$detail$.cv)))

  # A stratum mean of zero is known, not unknown: no relative quantity is
  # defined against it, and .rmoe says so the same way .cv does.
  zero <- n_alloc(
    data.frame(stratum = c("a", "b"), N = c(1000, 2000), sd = c(5, 8),
               mean = c(0, 30)),
    n = 300
  )
  expect_identical(zero$detail$.rmoe[1L], Inf)
  expect_identical(zero$detail$.cv[1L], Inf)
  expect_true(is.finite(zero$detail$.rmoe[2L]))
})

test_that("an allocation reports rmoe against its population mean", {
  frame <- data.frame(stratum = paste0("s", 1:3), N = c(1000, 2000, 1500),
                      sd = c(5, 8, 6), mean = c(20, 30, 25))
  fitted <- n_alloc(frame, n = 600)
  ybar <- stats::weighted.mean(frame$mean, frame$N)
  expect_equal(fitted$rmoe, fitted$moe / ybar, tolerance = 1e-12)
  expect_equal(prec_alloc(fitted)$rmoe, fitted$rmoe, tolerance = 1e-12)

  # A proportion frame measures against the same weighted estimand.
  prop_frame <- data.frame(stratum = c("a", "b"), N = c(1000, 3000),
                           sd = c(0.4, 0.5), p = c(0.2, 0.5))
  pf <- n_alloc(prop_frame, n = 400)
  expect_equal(pf$rmoe,
               pf$moe / stats::weighted.mean(prop_frame$p, prop_frame$N),
               tolerance = 1e-12)

  # A frame with no estimand leaves it NA, exactly as it leaves cv NA.
  bare <- n_alloc(data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                             sd = c(5, 8)), n = 300)
  expect_true(is.na(bare$rmoe))
  expect_true(is.na(bare$cv))
})

## Errors

test_that("rmoe needs the estimand it is relative to", {
  expect_error(n_mean(var = 100, rmoe = 0.1),
               "'mu' is required when 'rmoe' is specified")
  expect_error(n_multi(data.frame(name = "x", var = 100, rmoe = 0.1)),
               "'rmoe' rows need the estimand")
})

test_that("rmoe is exclusive with moe and with cv", {
  expect_error(n_prop(p = 0.2, rmoe = 0.1, moe = 0.02),
               "exactly one of 'moe', 'cv', or 'rmoe'")
  expect_error(n_prop(p = 0.2, rmoe = 0.1, cv = 0.05),
               "exactly one of 'moe', 'cv', or 'rmoe'")
  expect_error(n_mean(var = 100, mu = 50, rmoe = 0.1, cv = 0.05),
               "exactly one of 'moe', 'cv', or 'rmoe'")
  expect_error(prec_prop(n = 400, p = 0.3, rmoe = 0.1),
               "exactly one of 'p', 'cv', or 'rmoe'")
  expect_error(
    n_multi(data.frame(name = "x", p = 0.2, rmoe = 0.1, moe = 0.02)),
    "only one of 'rmoe' or 'moe'"
  )
  expect_error(
    n_multi(data.frame(name = "x", p = 0.2, rmoe = 0.1, cv = 0.05)),
    "only one of 'rmoe' or 'cv'"
  )
})

test_that("a target table row states exactly one of cv, moe, and rmoe", {
  frame <- data.frame(stratum = "a", N = 100)
  measures <- data.frame(stratum = "a", name = "A", p = 0.2)
  expect_error(
    n_alloc(frame, measures = measures,
            targets = data.frame(name = "A", cv = 0.1, rmoe = 0.1)),
    "exactly one positive finite 'cv', 'moe', or 'rmoe'"
  )
})

test_that("an objective cannot carry a precision requirement", {
  frame <- data.frame(stratum = c("a", "b"), N = c(100, 200))
  measures <- data.frame(stratum = c("a", "b"), name = "A", p = c(0.2, 0.3))
  expect_error(
    n_alloc(frame, measures = measures, budget = 1000,
            objective = data.frame(name = "A", rmoe = 0.1)),
    "put 'cv', 'moe', or 'rmoe' requirements in 'targets'"
  )
})

test_that("the MICS and surveyplanning spellings point at rmoe", {
  expect_error(n_multi(data.frame(name = "x", p = 0.2, rme = 0.12)),
               "Did you mean 'rmoe'")
  expect_error(n_multi(data.frame(name = "x", p = 0.2, RMoE = 0.12)),
               "Did you mean 'rmoe'")
})

test_that("rmoe must be positive and finite", {
  expect_error(n_prop(p = 0.2, rmoe = 0), "'rmoe' must be positive")
  expect_error(n_prop(p = 0.2, rmoe = Inf), "'rmoe' must be finite")
  expect_error(n_multi(data.frame(name = "x", p = 0.2, rmoe = -0.1)),
               "'rmoe' values must be positive and finite")
})

## predict on a sizing result

test_that("predict varies rmoe as a target and reports it back", {
  sized <- n_prop(p = 0.3, rmoe = 0.10)
  grid <- predict(sized, data.frame(rmoe = c(0.05, 0.10, 0.20)))
  expect_true("rmoe" %in% names(grid))
  expect_equal(grid$rmoe, c(0.05, 0.10, 0.20), tolerance = 1e-10)
  expect_equal(grid$n, sort(grid$n, decreasing = TRUE))
  expect_equal(grid$n[2L], sized$n, tolerance = 1e-10)
  expect_error(predict(sized, data.frame(rmoe = 0.1, cv = 0.05)),
               "cannot contain more than one of")
})

## rmoe as a round-trip override

test_that("the prec round trip accepts an rmoe override", {
  # The achieved margin of error is the implied target only when the caller
  # names none of the three; restoring it alongside an override would send
  # two targets into a function that takes one.
  m <- n_mean(prec_mean(var = 100, n = 400, mu = 50), rmoe = 0.05)
  expect_equal(m$rmoe, 0.05, tolerance = 1e-10)
  expect_equal(m$params$rmoe, 0.05)
  expect_null(m$params$moe)

  p <- n_prop(prec_prop(p = 0.3, n = 400), rmoe = 0.1)
  expect_equal(p$rmoe, 0.1, tolerance = 1e-10)
  expect_equal(p$params$rmoe, 0.1)
  expect_null(p$params$moe)
})

test_that("naming no target still implies the achieved moe", {
  expect_equal(n_mean(prec_mean(var = 100, n = 400, mu = 50))$n, 400,
               tolerance = 1e-9)
  expect_equal(n_prop(prec_prop(p = 0.3, n = 400))$n, 400, tolerance = 1e-9)
  expect_equal(n_mean(prec_mean(var = 100, n = 400, mu = 50), cv = 0.05)$cv,
               0.05, tolerance = 1e-10)
})
