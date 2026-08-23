## Kernel for ratio-of-totals planning. Expected values are derived here from
## the residual definition or by hand, never read back from the implementation.
## dev/ratio-probe.R is the wider reference these fixtures were taken from.

## A deterministic population, so the identities hold without an RNG.
ratio_pop <- function() {
  x <- c(2, 3, 5, 7, 11, 4, 6, 9, 8, 13, 1.5, 2.5)
  y <- c(7, 8, 16, 20, 34, 11, 19, 28, 23, 40, 5, 9)
  list(x = x, y = y)
}

test_that(".ratio_unit_relvar equals var(y - R*x) / mean(y)^2", {
  pop <- ratio_pop()
  r <- mean(pop$y) / mean(pop$x)
  e <- pop$y - r * pop$x

  value <- .ratio_unit_relvar(
    r = r,
    cv_num = sd(pop$y) / abs(mean(pop$y)),
    cv_den = sd(pop$x) / abs(mean(pop$x)),
    component_cor = cor(pop$y, pop$x)
  )

  expect_equal(value, var(e) / mean(pop$y)^2, tolerance = 1e-12)
})

test_that("r^2 times the unit relative variance is var(e) / mean(x)^2", {
  pop <- ratio_pop()
  r <- mean(pop$y) / mean(pop$x)
  e <- pop$y - r * pop$x

  value <- .ratio_unit_relvar(
    r = r,
    cv_num = sd(pop$y) / abs(mean(pop$y)),
    cv_den = sd(pop$x) / abs(mean(pop$x)),
    component_cor = cor(pop$y, pop$x)
  )

  expect_equal(r^2 * value, var(e) / mean(pop$x)^2, tolerance = 1e-12)
})

test_that("the sign of r flips the covariance term", {
  positive <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  negative <- .ratio_unit_relvar(-2, 1.1, 0.6, 0.7)

  expect_equal(positive, 1.21 + 0.36 - 2 * 0.7 * 1.1 * 0.6)
  expect_equal(negative, 1.21 + 0.36 + 2 * 0.7 * 1.1 * 0.6)
  expect_equal(positive, 0.646)
})

test_that("a mixed-sign population needs the signed form", {
  pop <- ratio_pop()
  y <- -pop$y
  x <- pop$x
  r <- mean(y) / mean(x)
  truth <- var(y - r * x) / mean(y)^2

  signed <- .ratio_unit_relvar(
    r = r,
    cv_num = sd(y) / abs(mean(y)),
    cv_den = sd(x) / abs(mean(x)),
    component_cor = cor(y, x)
  )
  unsigned <- (sd(y) / abs(mean(y)))^2 + (sd(x) / abs(mean(x)))^2 -
    2 * cor(y, x) * (sd(y) / abs(mean(y))) * (sd(x) / abs(mean(x)))

  expect_equal(signed, truth, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unsigned, truth)))
})

test_that("the unit relative variance is invariant to rescaling either part", {
  pop <- ratio_pop()
  relvar_of <- function(y, x) {
    .ratio_unit_relvar(
      r = mean(y) / mean(x),
      cv_num = sd(y) / abs(mean(y)),
      cv_den = sd(x) / abs(mean(x)),
      component_cor = cor(y, x)
    )
  }
  base <- relvar_of(pop$y, pop$x)

  expect_equal(relvar_of(3 * pop$y, pop$x), base, tolerance = 1e-12)
  expect_equal(relvar_of(pop$y, 5 * pop$x), base, tolerance = 1e-12)
  expect_equal(relvar_of(-pop$y, pop$x), base, tolerance = 1e-12)
  expect_equal(relvar_of(pop$y, -pop$x), base, tolerance = 1e-12)
})

test_that("proportional components give exactly zero, not a residue", {
  expect_identical(.ratio_unit_relvar(2, 0.8, 0.8, 1), 0)
  expect_identical(.ratio_unit_relvar(-2, 0.8, 0.8, -1), 0)
})

test_that("a negative unit relative variance is an internal error", {
  # Unreachable from valid moments, since L_R >= (cv_num - cv_den)^2. Only a
  # correlation outside its own bound can produce it, which is what the
  # validator exists to stop.
  expect_error(
    .ratio_unit_relvar(2, 1, 1, 1.5),
    "internal error"
  )
})

test_that("the cancellation clamp is tighter than any real difference", {
  # Just outside the tolerance the value survives rather than being zeroed.
  value <- .ratio_unit_relvar(2, 1, 1, 1 - 1e-9)
  expect_gt(value, 0)
  expect_equal(value, 2e-9, tolerance = 1e-6)
})

test_that(".check_ratio_relvar refuses a ratio with no sampling variance", {
  expect_error(
    .check_ratio_relvar(0, 0.8, 0.8, 1),
    "no sampling variance"
  )
})

test_that(".check_ratio_relvar warns just below the threshold and not above", {
  # cv_num = cv_den = 1 makes the threshold 1e-3 exactly.
  expect_warning(
    .check_ratio_relvar(.ratio_unit_relvar(2, 1, 1, 0.9996), 1, 1, 0.9996),
    "known to that precision"
  )
  expect_silent(
    .check_ratio_relvar(.ratio_unit_relvar(2, 1, 1, 0.999), 1, 1, 0.999)
  )
})

test_that("check_component_cor admits negatives and rejects out of range", {
  expect_silent(check_component_cor(-1))
  expect_silent(check_component_cor(0))
  expect_silent(check_component_cor(1))
  expect_error(check_component_cor(1.01), "must be a number in \\[-1, 1\\]")
  expect_error(check_component_cor(-1.01), "must be a number in \\[-1, 1\\]")
  expect_error(check_component_cor(NA_real_), "must be a number in \\[-1, 1\\]")
  expect_error(check_component_cor(c(0.2, 0.3)), "must be a number in \\[-1, 1\\]")
  expect_error(check_component_cor("0.5"), "must be a number in \\[-1, 1\\]")
})

test_that(".check_ratio_moments names a missing moment", {
  expect_error(.check_ratio_moments(NULL, 1.1, 0.6, 0.7), "needs 'r'")
  expect_error(.check_ratio_moments(2, NULL, 0.6, 0.7), "needs 'r'")
  expect_error(.check_ratio_moments(2, 1.1, NULL, 0.7), "needs 'r'")
  expect_error(.check_ratio_moments(2, 1.1, 0.6, NULL), "needs 'r'")
})

test_that(".check_ratio_moments rejects a zero ratio and a negative CV", {
  expect_error(.check_ratio_moments(0, 1.1, 0.6, 0.7), "must not be zero")
  expect_error(.check_ratio_moments(2, -1.1, 0.6, 0.7), "'cv_num' must be positive")
  expect_error(.check_ratio_moments(2, 1.1, 0, 0.7), "'cv_den' must be positive")
  expect_error(.check_ratio_moments(2, Inf, 0.6, 0.7), "'cv_num' must be finite")
  expect_silent(.check_ratio_moments(-2, 1.1, 0.6, -0.7))
})

test_that(".prec_engine_ratio matches the closed form across the design grid", {
  r <- 2
  relvar <- .ratio_unit_relvar(r, 1.1, 0.6, 0.7)
  grid <- expand.grid(
    n = c(80, 900),
    N = c(Inf, 20000),
    deff = c(1, 1.6),
    resp_rate = c(1, 0.75),
    alpha = c(0.05, 0.1),
    df = c(NA, 30)
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    dfv <- if (is.na(g$df)) NULL else g$df
    n_net <- g$n * g$resp_rate
    n_eff <- n_net / g$deff
    fpc <- if (is.infinite(g$N)) 1 else 1 - n_net / g$N
    q <- if (is.null(dfv)) {
      qnorm(1 - g$alpha / 2)
    } else {
      qt(1 - g$alpha / 2, dfv)
    }
    expected_se <- abs(r) * sqrt(relvar * fpc / n_eff)

    prec <- .prec_engine_ratio(r, relvar, g$n, g$alpha, g$N, g$deff,
                               g$resp_rate, dfv)

    expect_equal(prec$se, expected_se, tolerance = 1e-12)
    expect_equal(prec$moe, q * expected_se, tolerance = 1e-12)
    expect_equal(prec$cv, expected_se / abs(r), tolerance = 1e-12)
  }
})

test_that(".prec_engine_ratio reports on the ratio scale, not the y scale", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  prec <- .prec_engine_ratio(2, relvar, 500, 0.05, Inf, 1, 1)

  expect_named(prec, c("se", "moe", "cv"))
  expect_equal(prec$cv, prec$se / 2)
  expect_equal(prec$se, 2 * sqrt(relvar / 500), tolerance = 1e-12)
})

test_that("a negative ratio gives the same relative precision as its mirror", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  positive <- .prec_engine_ratio(2, relvar, 500, 0.05, Inf, 1, 1)
  negative <- .prec_engine_ratio(-2, relvar, 500, 0.05, Inf, 1, 1)

  expect_equal(positive$cv, negative$cv)
  expect_equal(positive$se, negative$se)
})

test_that(".n_ratio_from_target agrees across the three target modes", {
  r <- 2
  relvar <- .ratio_unit_relvar(r, 1.1, 0.6, 0.7)
  q <- qnorm(0.975)
  target_cv <- 0.05

  by_cv <- .n_ratio_from_target(r, relvar, NULL, target_cv, NULL,
                                0.05, Inf, 1, 1)
  by_moe <- .n_ratio_from_target(r, relvar, q * target_cv * abs(r), NULL, NULL,
                                 0.05, Inf, 1, 1)
  by_rmoe <- .n_ratio_from_target(r, relvar, NULL, NULL, q * target_cv,
                                  0.05, Inf, 1, 1)

  expect_equal(by_cv, relvar / target_cv^2, tolerance = 1e-10)
  expect_equal(by_moe, by_cv, tolerance = 1e-10)
  expect_equal(by_rmoe, by_cv, tolerance = 1e-10)
})

test_that(".n_ratio_from_target puts deff inside the FPC, not after it", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  deff <- 1.8
  N <- 4000

  got <- .n_ratio_from_target(2, relvar, NULL, 0.05, NULL, 0.05, N, deff, 1)

  expect_equal(got, deff * relvar / (0.05^2 + deff * relvar / N),
               tolerance = 1e-10)
  # Applying deff to a size already corrected for the frame is the wrong order
  # and gives a different answer, so the test can tell them apart.
  wrong <- deff * (relvar / (0.05^2 + relvar / N))
  expect_false(isTRUE(all.equal(got, wrong)))
})

test_that(".n_ratio_from_target inflates for response, and n is gross", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  net <- .n_ratio_from_target(2, relvar, NULL, 0.05, NULL, 0.05, Inf, 1, 1)
  gross <- .n_ratio_from_target(2, relvar, NULL, 0.05, NULL, 0.05, Inf, 1, 0.8)

  expect_equal(gross, net / 0.8, tolerance = 1e-10)
})

test_that(".n_ratio_from_target uses a t quantile when df is supplied", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  normal <- .n_ratio_from_target(2, relvar, 0.2, NULL, NULL, 0.05, Inf, 1, 1)
  student <- .n_ratio_from_target(2, relvar, 0.2, NULL, NULL, 0.05, Inf, 1, 1,
                                  df = 8)

  expect_gt(student, normal)
  expect_equal(
    student,
    relvar / (0.2 / (qt(0.975, 8) * 2))^2,
    tolerance = 1e-10
  )
})

test_that(".n_ratio_from_target refuses a target the frame cannot reach", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  expect_error(
    .n_ratio_from_target(2, relvar, NULL, 0.001, NULL, 0.05, 500, 1, 0.5),
    "unattainable"
  )
})

test_that("an rmoe target needs the ratio it is relative to", {
  relvar <- .ratio_unit_relvar(2, 1.1, 0.6, 0.7)
  expect_error(
    .ratio_target_cv(NULL, NULL, NULL, 0.1, qnorm(0.975)),
    "'r' is required when 'rmoe' is specified"
  )
})

test_that("the residual route and the moment route size identically", {
  # VDK Example 3.15 sizes the mean of y under a ratio model from the residual
  # variance. Section 12.2 of dev/PLAN-RATIO-ESTIMANDS-20260823.md claims that
  # is the same arithmetic as sizing the ratio itself in CV mode.
  pop <- ratio_pop()
  r <- mean(pop$y) / mean(pop$x)
  e <- pop$y - r * pop$x
  N <- length(pop$y)
  target_cv <- 0.05

  relvar <- .ratio_unit_relvar(
    r = r,
    cv_num = sd(pop$y) / abs(mean(pop$y)),
    cv_den = sd(pop$x) / abs(mean(pop$x)),
    component_cor = cor(pop$y, pop$x)
  )
  via_moments <- .n_ratio_from_target(r, relvar, NULL, target_cv, NULL,
                                      0.05, N, 1, 1)

  residual_relvar <- var(e) / mean(pop$y)^2
  via_residual <- residual_relvar / (target_cv^2 + residual_relvar / N)

  expect_equal(via_moments, via_residual, tolerance = 1e-10)
})
