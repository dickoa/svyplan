test_that("power inversion targets must exceed alpha", {
  for (alternative in c("one.sided", "two.sided")) {
    for (target in c(0.05, 0.025)) {
      expect_error(
        power_mean(var = 100, effect = 5, power = target,
                   alternative = alternative),
        "greater than 'alpha'"
      )
      expect_error(
        power_mean(var = 100, n = 100, effect = NULL, power = target,
                   alternative = alternative),
        "greater than 'alpha'"
      )
      expect_error(
        power_prop(0.3, p2 = 0.4, power = target,
                   alternative = alternative),
        "greater than 'alpha'"
      )
      expect_error(
        power_prop(0.3, n = 100, p2 = NULL, power = target,
                   alternative = alternative),
        "greater than 'alpha'"
      )
      expect_error(
        power_did(
          c(50, 55), control = c(50, 50), outcome = "mean", var = 100,
          effect = 5, power = target, alternative = alternative
        ),
        "greater than 'alpha'"
      )
      expect_error(
        power_did(
          c(50, 55), control = c(50, 50), outcome = "mean", var = 100,
          effect = NULL, n = 100, power = target,
          alternative = alternative
        ),
        "greater than 'alpha'"
      )
    }
  }
})

test_that("solved power sizes respect two units in both groups", {
  mean_fit <- power_mean(var = 1, effect = 10, power = 0.8)
  expect_equal(mean_fit$n, 2)
  expect_gt(mean_fit$power, 0.8)

  prop_fit <- power_prop(0.01, p2 = 0.99, power = 0.8)
  expect_equal(prop_fit$n, 2)
  expect_gt(prop_fit$power, 0.8)

  did_fit <- power_did(
    c(0, 10), control = c(0, 0), outcome = "mean", var = 1,
    effect = 10, power = 0.8, ratio = 0.1
  )
  expect_equal(unname(did_fit$n), c(2, 20))
  expect_gt(did_fit$power, 0.8)
})

test_that("two-sided sizing stores power achieved by its returned n", {
  for (N in c(Inf, 1e6)) {
    for (target in c(0.06, 0.999)) {
      mean_fit <- power_mean(var = 1, effect = 0.2, power = target, N = N)
      mean_back <- power_mean(
        var = 1, effect = 0.2, n = mean_fit$n, power = NULL, N = N
      )
      expect_equal(mean_fit$power, mean_back$power, tolerance = 1e-8)
      expect_equal(mean_fit$power, target, tolerance = 1e-8)

      prop_fit <- power_prop(0.3, p2 = 0.4, power = target, N = N)
      prop_back <- power_prop(
        0.3, p2 = 0.4, n = prop_fit$n, power = NULL, N = N
      )
      expect_equal(prop_fit$power, prop_back$power, tolerance = 1e-8)
      expect_equal(prop_fit$power, target, tolerance = 1e-8)

      did_fit <- power_did(
        c(0, 0.2), control = c(0, 0), outcome = "mean", var = 1,
        effect = 0.2, power = target, N = N
      )
      did_back <- power_did(
        c(0, 0.2), control = c(0, 0), outcome = "mean", var = 1,
        effect = 0.2, n = did_fit$n, power = NULL, N = N
      )
      expect_equal(did_fit$power, did_back$power, tolerance = 1e-8)
      expect_equal(did_fit$power, target, tolerance = 1e-8)
    }
  }
})

## Every proportion scale carries the finite-population Bernoulli variance,
## so the transformed scales use (N - n)/(N - 1), not 1 - n/N.

test_that("finite-population proportion power matches a hand oracle on every scale", {
  p1 <- 0.2
  p2 <- 0.6
  N <- 10
  n <- 5
  z_a <- qnorm(0.975)
  fpc <- (N - n) / (N - 1)

  two_tail <- function(effect, se) {
    pnorm(abs(effect) / se - z_a) + pnorm(-abs(effect) / se - z_a)
  }

  # p(1-p) with (N - n)/(N - 1), which is how N/(N - 1) times 1 - n/N lands.
  se_wald <- sqrt(sum(c(p1 * (1 - p1), p2 * (1 - p2)) * fpc / n))
  expect_equal(
    power_prop(p1, p2 = p2, n = n, power = NULL, N = N, method = "wald")$power,
    two_tail(p1 - p2, se_wald),
    tolerance = 1e-10
  )
  expect_equal(
    .bernoulli_var(p1, N) * .fpc_factor(n, N),
    p1 * (1 - p1) * .fpc_factor_prop(n, N),
    tolerance = 1e-12
  )

  phi <- asin(sqrt(p1)) - asin(sqrt(p2))
  expect_equal(
    power_prop(p1, p2 = p2, n = n, power = NULL, N = N,
               method = "arcsine")$power,
    two_tail(phi, sqrt(2 * fpc / (4 * n))),
    tolerance = 1e-10
  )

  q1 <- 1 - p1
  q2 <- 1 - p2
  p_bar <- (p1 + p2) / 2
  q_bar <- 1 - p_bar
  lambda <- log(p1 / q1) - log(p2 / q2)
  V0 <- 2 * fpc / (n * p_bar * q_bar)
  VA <- fpc / (n * p1 * q1) + fpc / (n * p2 * q2)
  expect_equal(
    power_prop(p1, p2 = p2, n = n, power = NULL, N = N,
               method = "logodds")$power,
    pnorm((abs(lambda) - z_a * sqrt(V0)) / sqrt(VA)) +
      pnorm((-abs(lambda) - z_a * sqrt(V0)) / sqrt(VA)),
    tolerance = 1e-10
  )
})

test_that("a large finite population reproduces the infinite-population power", {
  for (method in c("wald", "arcsine", "logodds")) {
    expect_equal(
      power_prop(0.2, p2 = 0.6, n = 50, power = NULL, N = 1e9,
                 method = method)$power,
      power_prop(0.2, p2 = 0.6, n = 50, power = NULL, method = method)$power,
      tolerance = 1e-6,
      info = method
    )
  }
})

## A census has no sampling variance, so no effect has the requested power.

test_that("a census refuses to report a minimum detectable effect", {
  expect_error(
    power_mean(var = 100, n = 100, N = 100, power = 0.8, effect = NULL),
    "no minimum detectable effect exists"
  )
  for (method in c("wald", "arcsine", "logodds")) {
    expect_error(
      power_prop(0.3, n = 100, N = 100, power = 0.8, p2 = NULL,
                 method = method),
      "no minimum detectable effect exists",
      info = method
    )
  }
  expect_error(
    power_did(
      treat = c(0.3, 0.4), control = c(0.3, 0.3), outcome = "prop",
      n = 100, N = 100, power = 0.8, effect = NULL
    ),
    "no minimum detectable effect exists"
  )

  # One unit short of a census still has variance and still solves.
  expect_gt(
    power_mean(var = 100, n = 99, N = 100, power = 0.8, effect = NULL)$effect,
    0
  )
})

## The DiD finite-population correction subtracts a census term that differs
## from the per-unit numerator at partial overlap.

test_that("the documented DiD variance decomposition is the implemented one", {
  n_eff <- c(300, 300)
  N_pair <- c(1000, 1000)
  terms <- .did_var_terms_prop(
    c(0.3, 0.4), c(0.3, 0.32), 0.5, 0.6, N_pair
  )

  expect_equal(
    .did_var_d(n_eff, terms, N_pair, 1),
    sum(terms$change / n_eff - terms$census / N_pair),
    tolerance = 1e-12
  )

  # The factor form is correct only at the two extremes.
  for (overlap in c(0, 1)) {
    ends <- .did_var_terms_prop(
      c(0.3, 0.4), c(0.3, 0.32), overlap, 0.6, N_pair
    )
    expect_equal(
      .did_var_d(n_eff, ends, N_pair, 1),
      sum(ends$change * (1 - n_eff / N_pair) / n_eff),
      tolerance = 1e-12,
      info = as.character(overlap)
    )
  }
  expect_false(isTRUE(all.equal(
    .did_var_d(n_eff, terms, N_pair, 1),
    sum(terms$change * (1 - n_eff / N_pair) / n_eff)
  )))
})
