test_that("power_did prop: solve for n matches the closed-form starting value", {
  treat <- c(0.50, 0.55)
  control <- c(0.50, 0.48)
  effect <- abs((treat[2] - treat[1]) - (control[2] - control[1]))
  # Per-arm variance (no overlap): V_arm = p0(1-p0) + p1(1-p1)
  V_trt <- treat[1] * (1 - treat[1]) + treat[2] * (1 - treat[2])
  V_ctrl <- control[1] * (1 - control[1]) + control[2] * (1 - control[2])
  V <- V_trt + V_ctrl
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  n_expected <- (z_a + z_b)^2 * V / effect^2

  res <- power_did(treat = treat, control = control, outcome = "prop",
                    effect = effect)
  expect_equal(res$n, n_expected, tolerance = 1e-4)
  expect_equal(res$effect, effect)
  expect_equal(res$power, 0.80)
  expect_equal(res$solved, "n")
  expect_equal(res$type, "did_prop")
})

test_that("power_did mean: solve for n matches the closed-form starting value", {
  var <- 100
  effect <- 5
  V <- 4 * var  # 4 cells equal variance, no overlap
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  n_expected <- (z_a + z_b)^2 * V / effect^2

  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = var, effect = effect
  )
  expect_equal(res$n, n_expected, tolerance = 1e-4)
  expect_equal(res$power, 0.8, tolerance = 1e-8)
  expect_equal(res$type, "did_mean")
  expect_equal(res$solved, "n")
})

test_that("power_did mean: length-2 var", {
  var2 <- c(80, 120)
  effect <- 5
  # var_parts = c(80, 80, 120, 120), V_trt = 160, V_ctrl = 240, sum = 400
  V <- 2 * var2[1] + 2 * var2[2]
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  n_expected <- (z_a + z_b)^2 * V / effect^2

  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = var2, effect = effect
  )
  expect_equal(res$n, n_expected, tolerance = 1e-4)
})

test_that("power_did mean: length-4 var", {
  var4 <- c(80, 120, 90, 110)
  effect <- 5
  V <- (var4[1] + var4[2]) + (var4[3] + var4[4])
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  n_expected <- (z_a + z_b)^2 * V / effect^2

  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = var4, effect = effect
  )
  expect_equal(res$n, n_expected, tolerance = 1e-4)
})

test_that("power_did: solve for power", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, power = NULL, n = 500
  )
  expect_equal(res$solved, "power")
  expect_true(res$power > 0 && res$power < 1)
  expect_equal(res$n, 500)
})

test_that("power_did: solve for MDE (prop)", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = NULL, n = 500, power = 0.80
  )
  expect_equal(res$solved, "mde")
  expect_true(res$effect > 0)
  expect_equal(res$n, 500)
  expect_equal(res$power, 0.80)
})

test_that("power_did: solve for MDE (mean)", {
  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = NULL, n = 500, power = 0.80
  )
  expect_equal(res$solved, "mde")
  expect_true(res$effect > 0)

  # Verify against manual formula
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  V <- 4 * 100
  se <- sqrt(V / 500)
  mde_expected <- (z_a + z_b) * se
  expect_equal(res$effect, mde_expected, tolerance = 1e-4)
  expect_equal(
    power_did(
      treat = c(50, 50 + res$effect), control = c(50, 50), outcome = "mean",
      var = 100, effect = res$effect, n = 500, power = NULL
    )$power,
    0.8,
    tolerance = 1e-8
  )
})

test_that("power_did round-trip: n -> power -> n", {
  res1 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, power = 0.90
  )
  res2 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, power = NULL, n = res1$n
  )
  expect_equal(res2$power, 0.90, tolerance = 1e-4)
})

test_that("power_did round-trip: n -> MDE -> n (mean)", {
  var4 <- c(80, 120, 90, 110)
  effect <- 5
  res1 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = var4, effect = effect, power = 0.80
  )
  res2 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = var4, effect = NULL, n = res1$n, power = 0.80
  )
  expect_equal(res2$effect, effect, tolerance = 1e-4)
})

test_that("power_did: deff multiplier increases n", {
  res1 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  res2 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, deff = 2
  )
  expect_equal(res2$n, res1$n * 2, tolerance = 1e-6)
})

test_that("power_did: overlap reduces n", {
  res1 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  res2 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, overlap = 0.5, overlap_cor = 0.6
  )
  expect_true(res2$n < res1$n)
})

test_that("power_did: per-arm overlap formula gives same as flat when overlap=0", {
  # With overlap=0, per-arm formula reduces to V_arm = V_pre + V_post
  # which is the same as the flat formula
  res_ov0 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, overlap = 0, overlap_cor = 0.5
  )
  res_no_ov <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  expect_equal(res_ov0$n, res_no_ov$n, tolerance = 1e-10)
})

test_that("power_did finite-N overlap zero ignores overlap_cor", {
  power_at <- function(rho) {
    power_did(
      treat = c(50, 55), control = c(50, 50), outcome = "mean",
      var = 100, effect = 5, n = 100, power = NULL, N = 1000,
      overlap = 0, overlap_cor = rho
    )$power
  }
  expect_equal(power_at(0), power_at(0.3), tolerance = 1e-12)
  expect_equal(power_at(0), power_at(0.9), tolerance = 1e-12)
})

test_that("power_did proportions use finite-population Bernoulli variances", {
  treat <- c(0.2, 0.6)
  control <- c(0.3, 0.4)
  N <- c(10, 12)
  n <- c(5, 6)
  effect <- abs(diff(treat) - diff(control))
  S2_t <- N[1] / (N[1] - 1) * treat * (1 - treat)
  S2_c <- N[2] / (N[2] - 1) * control * (1 - control)
  V <- sum(S2_t) * (1 / n[1] - 1 / N[1]) +
    sum(S2_c) * (1 / n[2] - 1 / N[2])
  z_a <- qnorm(0.975)
  expected <- pnorm(effect / sqrt(V) - z_a) +
    pnorm(-effect / sqrt(V) - z_a)

  got <- power_did(
    treat, control = control, outcome = "prop", effect = effect,
    n = n, power = NULL, N = N
  )
  expect_equal(got$power, expected, tolerance = 1e-12)
})

test_that("power_did proportions check each arm's attainable correlation", {
  expect_error(
    power_did(
      treat = c(0.01, 0.9), control = c(0.4, 0.4), outcome = "prop",
      effect = 0.89, n = 100, power = NULL,
      overlap = 1, overlap_cor = 0.9
    ),
    "0.0335"
  )
  expect_error(
    power_did(
      treat = c(0.4, 0.4), control = c(0.01, 0.9), outcome = "prop",
      effect = 0.89, n = 100, power = NULL,
      overlap = 1, overlap_cor = 0.9
    ),
    "0.0335"
  )
})

test_that("power_did: one-sided gives smaller n", {
  res1 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, alternative = "two.sided"
  )
  res2 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, alternative = "one.sided"
  )
  expect_true(res2$n < res1$n)
})

test_that("power_did: FPC reduces n (finite N)", {
  res_inf <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5, N = Inf
  )
  res_fin <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5, N = 5000
  )
  expect_true(res_fin$n < res_inf$n)
})

test_that("power_did: FPC round-trip", {
  res1 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5, N = 2000, power = 0.80
  )
  res2 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5, power = NULL, n = res1$n, N = 2000
  )
  expect_equal(res2$power, 0.80, tolerance = 1e-3)
})

test_that("power_did: ratio (unequal groups)", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, ratio = 2
  )
  expect_equal(length(res$n), 2L)
  expect_equal(res$n[1] / res$n[2], 2, tolerance = 1e-6)
})

test_that("power_did: ratio round-trip", {
  res1 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, ratio = 2, power = 0.80
  )
  res2 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, power = NULL, n = res1$n
  )
  expect_equal(res2$power, 0.80, tolerance = 1e-4)
})

test_that("power_did: resp_rate inflates n", {
  res1 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  res2 <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, resp_rate = 0.8
  )
  expect_equal(res2$n, res1$n / 0.8, tolerance = 1e-6)
})

test_that("power_did: census guard with finite N", {
  full <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5, power = NULL, n = 100, N = 100
  )
  expect_equal(full$power, 1)

  partial <- power_did(
    treat = c(50, 55), control = c(50, 50), outcome = "mean",
    var = 100, effect = 5, power = NULL,
    n = c(100, 20), N = c(100, 100)
  )
  expect_lt(partial$power, 1)

  expect_error(
    power_did(
      treat = c(50, 55), control = c(50, 50),
      outcome = "mean", var = 100, effect = 5, power = NULL, n = 5000, N = 100
    ),
    "cannot draw"
  )
})

test_that("power_did: print output", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  out <- capture.output(print(res))
  expect_true(any(grepl("DiD proportions", out)))
  expect_true(any(grepl("sample size", out)))
  expect_true(any(grepl("treat", out)))
  expect_true(any(grepl("control", out)))
})

test_that("power_did: format output", {
  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5
  )
  fmt <- format(res)
  expect_true(grepl("did_mean", fmt))
})

test_that("power_did: plot runs without error", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  pdf(tempfile())
  on.exit(dev.off())
  expect_no_error(plot(res))
})

test_that("power_did: predict grid (power solve)", {
  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5, power = NULL, n = 500
  )
  grid <- predict(res, data.frame(n = c(300, 500, 700)))
  expect_equal(nrow(grid), 3L)
  expect_true("power" %in% names(grid))
  expect_true(grid$power[3] > grid$power[1])
})

test_that("power_did: predict for MDE", {
  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = NULL, n = 500, power = 0.80
  )
  grid <- predict(res, data.frame(n = c(300, 500, 700)))
  expect_equal(nrow(grid), 3L)
  expect_true("effect" %in% names(grid))
  expect_true(grid$effect[1] > grid$effect[3])
})

test_that("power_did: predict for n solve", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, power = 0.80
  )
  grid <- predict(res, data.frame(power = c(0.70, 0.80, 0.90)))
  expect_equal(nrow(grid), 3L)
  expect_true(grid$n[3] > grid$n[1])
})

test_that("power_did: coercion", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07
  )
  expect_equal(as.double(res), res$n)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
})

test_that("power_did: validation errors", {
  expect_error(
    power_did(treat = c(0.5), control = c(0.5, 0.5), outcome = "prop"),
    "length 2"
  )
  expect_error(
    power_did(treat = c(0.5, 1.5), control = c(0.5, 0.5), outcome = "prop"),
    "\\(0, 1\\)"
  )
  expect_error(
    power_did(treat = c(50, 55), control = c(50, 50), outcome = "mean",
              var = c(1, 2, 3), effect = 5),
    "length 1, 2, or 4"
  )
  expect_error(
    power_did(treat = c(50, 55), control = c(50, 50), outcome = "mean",
              var = -1, effect = 5),
    "positive"
  )
  expect_error(
    power_did(treat = c(50, 55), control = c(50, 50), outcome = "mean",
              effect = 5),
    "required"
  )
  expect_error(
    power_did(treat = c(0.5, 0.55), control = c(0.5, 0.48),
              outcome = "prop", n = NULL, power = NULL),
    "exactly one"
  )
  expect_error(
    power_did(treat = c(0.5, 0.55), control = c(0.5, 0.48),
              outcome = "prop", power = NULL, n = 100, ratio = 2),
    "ratio"
  )
})

test_that("power_did: scalar var equals length-4 equal", {
  res1 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5
  )
  res2 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = c(100, 100, 100, 100), effect = 5
  )
  expect_equal(res1$n, res2$n, tolerance = 1e-10)
})

test_that("power_did: unequal var gives different n than equal var", {
  res1 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5
  )
  res2 <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = c(80, 120, 90, 150), effect = 5
  )
  expect_false(abs(res1$n - res2$n) < 1e-6)
})

test_that("power_did: explicit effect overrides derived, with a warning", {
  expect_warning(
    res <- power_did(
      treat = c(0.50, 0.55), control = c(0.50, 0.48),
      outcome = "prop", effect = 0.10, power = NULL, n = 300
    ),
    "differs from the difference-in-differences implied"
  )
  expect_equal(res$effect, 0.10)
})

test_that("power_did derives effect from treat and control when omitted", {
  for (spec in list(
    list(outcome = "prop", treat = c(0.50, 0.55), control = c(0.50, 0.48),
         implied = 0.07),
    list(outcome = "mean", treat = c(50, 55), control = c(50, 50),
         var = 100, implied = 5)
  )) {
    args <- spec[c("treat", "control", "outcome")]
    if (!is.null(spec$var)) args$var <- spec$var

    solved_n <- do.call(power_did, args)
    expect_equal(solved_n$effect, spec$implied)
    expect_equal(solved_n$solved, "n")
    expect_equal(
      solved_n$n,
      do.call(power_did, c(args, list(effect = spec$implied)))$n
    )

    solved_power <- do.call(power_did, c(args, list(n = 400, power = NULL)))
    expect_equal(solved_power$effect, spec$implied)
    expect_equal(solved_power$solved, "power")

    # Both n and power supplied still means "solve for the MDE".
    mde <- do.call(power_did, c(args, list(n = 400, power = 0.80)))
    expect_equal(mde$solved, "mde")
  }
})

test_that("power_did rejects a zero implied effect when effect is omitted", {
  expect_error(
    power_did(treat = c(50, 55), control = c(50, 55),
              outcome = "mean", var = 100),
    "imply a difference-in-differences of 0"
  )
})

test_that("power_did still requires exactly one unknown", {
  expect_error(
    power_did(treat = c(50, 55), control = c(50, 50),
              outcome = "mean", var = 100, n = NULL, power = NULL,
              effect = 5),
    "leave exactly one"
  )
})

test_that("power_did: plot for did_mean", {
  res <- power_did(
    treat = c(50, 55), control = c(50, 50),
    outcome = "mean", var = 100, effect = 5
  )
  pdf(tempfile())
  on.exit(dev.off())
  expect_no_error(plot(res))
})

test_that("power_did: predict for did_prop", {
  res <- power_did(
    treat = c(0.50, 0.55), control = c(0.50, 0.48),
    outcome = "prop", effect = 0.07, power = NULL, n = 500
  )
  grid <- predict(res, data.frame(n = c(300, 500, 700)))
  expect_equal(nrow(grid), 3L)
})

test_that("power_did: params stored correctly", {
  treat <- c(0.50, 0.55)
  control <- c(0.50, 0.48)
  res <- power_did(treat = treat, control = control, outcome = "prop",
                    effect = 0.07)
  expect_equal(res$params$treat, treat)
  expect_equal(res$params$control, control)
  expect_equal(res$params$outcome, "prop")
  expect_null(res$params$var)
})

test_that("var and sd are rejected under outcome = 'prop'", {
  expect_error(
    power_did(treat = c(0.2, 0.3), control = c(0.2, 0.25),
              outcome = "prop", var = 100),
    "apply only to outcome"
  )
  expect_error(
    power_did(treat = c(0.2, 0.3), control = c(0.2, 0.25),
              outcome = "prop", sd = 10),
    "apply only to outcome"
  )
  expect_s3_class(
    power_did(treat = c(0.2, 0.3), control = c(0.2, 0.25), outcome = "prop"),
    "svyplan_power"
  )
})

test_that("sd is an alternative spelling of var under outcome = 'mean'", {
  a <- power_did(treat = c(50, 55), control = c(50, 52),
                 outcome = "mean", var = 100)
  b <- power_did(treat = c(50, 55), control = c(50, 52),
                 outcome = "mean", sd = 10)
  expect_equal(a$n, b$n)
})

test_that("the within-arm overlap correction uses that arm's own N", {
  # each arm is its own population, and inside it the two occasions share
  # units, so the covariance carries 1/N for that arm
  za <- qnorm(0.975)
  vp <- c(4, 4, 9, 9); nn <- c(50, 70); NN <- c(300, 500)
  for (rho in c(0.5, 1)) {
    for (ov in c(0.4, 0.9)) {
      arm <- function(v0, v1, n, N) {
        v0 * (1 / n - 1 / N) + v1 * (1 / n - 1 / N) -
          2 * rho * sqrt(v0 * v1) * (ov / n - 1 / N)
      }
      exact <- arm(vp[1], vp[2], nn[1], NN[1]) + arm(vp[3], vp[4], nn[2], NN[2])
      got <- power_did(treat = c(10, 11), control = c(10, 10), outcome = "mean",
                       var = vp, effect = 1, n = nn, power = NULL, N = NN,
                       overlap = ov, overlap_cor = rho)
      se <- 1 / (qnorm(got$power - pnorm(-1 / sqrt(exact) - za)) + za)
      expect_equal(se^2, exact, tolerance = 1e-7)
    }
  }
})
