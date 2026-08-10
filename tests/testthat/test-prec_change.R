## P1. The variance the two occasions actually have

test_that("prec_change at overlap 0 is two independent samples", {
  res <- prec_change(var = 100, n = 500)
  expect_s3_class(res, "svyplan_prec")
  expect_equal(res$type, "change")
  expect_equal(res$se, sqrt(100 / 500 + 100 / 500), tolerance = 1e-12)
  expect_equal(res$moe, qnorm(0.975) * res$se, tolerance = 1e-12)
})

test_that("prec_change subtracts the overlap covariance", {
  res <- prec_change(var = 100, n = 500, overlap = 0.5, overlap_cor = 0.6)
  v <- 100 / 500 + 100 / 500 - 2 * 0.5 * 0.6 * 100 / 500
  expect_equal(res$se, sqrt(v), tolerance = 1e-12)
})

test_that("overlap and overlap_cor only ever act together", {
  base <- prec_change(var = 100, n = 900)$se
  expect_equal(prec_change(var = 100, n = 900, overlap = 1, overlap_cor = 0)$se,
               base, tolerance = 1e-12)
  expect_equal(prec_change(var = 100, n = 900, overlap = 0, overlap_cor = 1)$se,
               base, tolerance = 1e-12)
})

test_that("a full correlated panel measures the change without error", {
  res <- prec_change(var = 100, n = 500, overlap = 1, overlap_cor = 1)
  expect_equal(res$se, 0)
  expect_equal(res$moe, 0)
})

test_that("two censuses of one population measure the change exactly", {
  expect_equal(prec_change(var = 100, n = 5000, N = 5000)$se, 0)
})

test_that("the overlap covariance carries a single 1/N", {
  # The marginal terms take their own fpc; the covariance enters the
  # population once, which is what makes the census cancellation exact.
  n <- 400
  N <- 10000
  rho <- 0.7
  ov <- 0.5
  v <- 100
  ne <- n
  expected <- v * (1 - ne / N) / ne + v * (1 - ne / N) / ne -
    2 * ov * rho * v / ne + 2 * rho * v / N
  res <- prec_change(var = v, n = n, N = N, overlap = ov, overlap_cor = rho)
  expect_equal(res$se^2, expected, tolerance = 1e-12)
})

## P2. The two scales

test_that("the proportion scale is p(1 - p) on each occasion", {
  a <- prec_change(p = c(0.3, 0.36), n = 900)
  b <- prec_change(var = c(0.3 * 0.7, 0.36 * 0.64), n = 900)
  expect_equal(a$se, b$se, tolerance = 1e-12)
})

test_that("the proportion scale determines the change", {
  res <- prec_change(p = c(0.3, 0.36), n = 900)
  expect_equal(res$params$change, 0.06, tolerance = 1e-12)
  expect_equal(res$cv, res$se / 0.06, tolerance = 1e-12)
  expect_equal(res$rmoe, res$moe / 0.06, tolerance = 1e-12)
})

test_that("a negative change keeps positive relative measures", {
  res <- prec_change(p = c(0.36, 0.30), n = 900)
  expect_equal(res$params$change, -0.06, tolerance = 1e-12)
  expect_gt(res$cv, 0)
  expect_gt(res$rmoe, 0)
})

test_that("relative measures are NA without a known change", {
  res <- prec_change(var = 100, n = 500)
  expect_false(is.na(res$se))
  expect_true(is.na(res$cv))
  expect_true(is.na(res$rmoe))
})

test_that("sd is an alternative spelling of var", {
  expect_equal(prec_change(sd = 10, n = 500)$se,
               prec_change(var = 100, n = 500)$se, tolerance = 1e-12)
  expect_equal(prec_change(sd = c(10, 12), n = 500)$se,
               prec_change(var = c(100, 144), n = 500)$se, tolerance = 1e-12)
})

## P3. Unequal occasions, where overlap becomes directional

test_that("overlap counts against the first occasion", {
  res <- prec_change(var = 100, n = c(800, 400), overlap = 0.5,
                     overlap_cor = 0.6)
  v <- 100 / 800 + 100 / 400 - 2 * 0.5 * 0.6 * 100 / 400
  expect_equal(res$se, sqrt(v), tolerance = 1e-12)
})

test_that("overlap cannot exceed the share the second occasion can hold", {
  expect_error(
    prec_change(var = 100, n = c(1000, 400), overlap = 0.9, overlap_cor = 0.5),
    "must be <= n\\[2\\]/n\\[1\\]"
  )
  expect_silent(
    prec_change(var = 100, n = c(1000, 400), overlap = 0.4, overlap_cor = 0.5)
  )
})

## P4. deff, response, and the interval quantile

test_that("deff multiplies the assembled change variance", {
  a <- prec_change(var = 100, n = 500, overlap = 0.5, overlap_cor = 0.6)
  b <- prec_change(var = 100, n = 500, overlap = 0.5, overlap_cor = 0.6,
                   deff = 2)
  expect_equal(b$se^2, 2 * a$se^2, tolerance = 1e-12)
})

test_that("resp_rate nets both occasions down before the variance forms", {
  a <- prec_change(var = 100, n = 1000, resp_rate = 0.5)
  b <- prec_change(var = 100, n = 500)
  expect_equal(a$se, b$se, tolerance = 1e-12)
})

test_that("df switches the interval quantile from normal to t", {
  a <- prec_change(var = 100, n = 500)
  b <- prec_change(var = 100, n = 500, df = 30)
  expect_equal(a$se, b$se, tolerance = 1e-12)
  expect_equal(b$moe, qt(0.975, 30) * b$se, tolerance = 1e-12)
  expect_gt(b$moe, a$moe)
})

## P5. Agreement with the power family, which shares the variance

test_that("prec_change moe equals the power_mean MDE at power 0.5", {
  cfgs <- list(
    list(N = Inf, ov = 0, rho = 0, n = 900, rr = 1, deff = 1),
    list(N = Inf, ov = 0.5, rho = 0.6, n = 900, rr = 1, deff = 1),
    list(N = 20000, ov = 0.75, rho = 0.5, n = 900, rr = 0.8, deff = 1.5),
    list(N = Inf, ov = 0.4, rho = 0.9, n = c(1800, 900), rr = 1, deff = 1),
    list(N = 50000, ov = 0, rho = 0, n = c(1350, 900), rr = 0.7, deff = 2)
  )
  for (cfg in cfgs) {
    pw <- power_mean(var = 100, n = cfg$n, power = 0.5, N = cfg$N,
                     deff = cfg$deff, resp_rate = cfg$rr,
                     overlap = cfg$ov, overlap_cor = cfg$rho)
    pc <- prec_change(var = 100, n = cfg$n, N = cfg$N, deff = cfg$deff,
                      resp_rate = cfg$rr, overlap = cfg$ov,
                      overlap_cor = cfg$rho)
    expect_equal(pw$effect, pc$moe, tolerance = 1e-10)
  }
})

test_that("prec_change se reproduces the power_prop power", {
  for (ov in c(0, 0.5, 0.75)) {
    pc <- prec_change(p = c(0.30, 0.36), n = 2000, overlap = ov,
                      overlap_cor = 0.5)
    pw <- power_prop(p1 = 0.30, p2 = 0.36, n = 2000, power = NULL,
                     overlap = ov, overlap_cor = 0.5)
    z <- qnorm(0.975)
    manual <- pnorm(0.06 / pc$se - z) + pnorm(-0.06 / pc$se - z)
    expect_equal(pw$power, manual, tolerance = 1e-10)
  }
})

## P6. Rejected inputs

test_that("prec_change requires exactly one dispersion scale", {
  expect_error(prec_change(n = 500), "exactly one of 'var'")
  expect_error(prec_change(var = 100, p = c(0.3, 0.36), n = 500),
               "exactly one of 'var'")
})

test_that("change is refused on the proportion scale", {
  expect_error(prec_change(p = c(0.3, 0.36), change = 0.06, n = 500),
               "do not supply it")
})

test_that("p must be two proportions", {
  expect_error(prec_change(p = 0.3, n = 500), "two proportions")
  expect_error(prec_change(p = c(0.3, 0.36, 0.4), n = 500), "two proportions")
  expect_error(prec_change(p = c(0, 0.36), n = 500), "must be in \\(0, 1\\)")
})

test_that("overlapping occasions must share one population", {
  expect_error(
    prec_change(var = 100, n = 500, N = c(10000, 20000), overlap = 0.5,
                overlap_cor = 0.5),
    "drawn from one population"
  )
  expect_silent(prec_change(var = 100, n = 500, N = c(10000, 20000)))
})

test_that("n is one size or one per occasion", {
  expect_error(prec_change(var = 100, n = c(500, 400, 300)),
               "length 1 \\(equal occasions\\)")
  expect_error(prec_change(var = 100, n = -5), "must be positive")
})

test_that("prec_change rejects unused arguments", {
  expect_error(prec_change(var = 100, n = 500, nope = 1), "unused argument")
})

## P7. Object surface

test_that("params keep the pair and the scale they were built from", {
  res <- prec_change(p = c(0.3, 0.36), n = 900, overlap = 0.5,
                     overlap_cor = 0.6)
  expect_equal(res$params$p, c(0.3, 0.36))
  expect_equal(res$params$var, c(0.3 * 0.7, 0.36 * 0.64), tolerance = 1e-12)
  expect_equal(res$params$overlap, 0.5)
  expect_equal(res$params$overlap_cor, 0.6)

  res_mean <- prec_change(var = 100, n = 900)
  expect_null(res_mean$params$p)
  expect_equal(res_mean$params$var, c(100, 100))
})

test_that("confint centres the interval on the change", {
  res <- prec_change(p = c(0.3, 0.36), n = 2000)
  ci <- confint(res)
  expect_equal(unname(ci[1L]), 0.06 - res$moe, tolerance = 1e-12)
  expect_equal(unname(ci[2L]), 0.06 + res$moe, tolerance = 1e-12)
})

test_that("confint needs a change to centre on", {
  expect_error(confint(prec_change(var = 100, n = 500)),
               "required to compute a confidence interval for a change")
})

test_that("print reports the scale, the size and the overlap", {
  out <- capture.output(print(prec_change(p = c(0.3, 0.36), n = 2000,
                                          overlap = 0.75, overlap_cor = 0.5)))
  expect_match(out[1L], "proportion scale")
  expect_match(out[2L], "n = 2000 per occasion")
  expect_true(any(grepl("overlap = 0.75", out)))

  flat <- capture.output(print(prec_change(var = 100, n = 500)))
  expect_match(flat[1L], "mean scale")
  expect_true(any(grepl("No between-occasion covariance", flat)))

  pair <- capture.output(print(prec_change(var = 100, n = c(800, 400))))
  expect_match(pair[2L], "n = 800 then 400")
})

test_that("a plan supplies overlap defaults", {
  pl <- svyplan(overlap = 0.5, overlap_cor = 0.6, alpha = 0.10)
  expect_equal(
    prec_change(var = 100, n = 500, plan = pl)$moe,
    prec_change(var = 100, n = 500, overlap = 0.5, overlap_cor = 0.6,
                alpha = 0.10)$moe,
    tolerance = 1e-12
  )
})

## P8. Agreement with prec_prop at a finite N

test_that("the proportion scale uses the same finite population variance", {
  # p(1-p) is the Bernoulli variance; the change formula wants the
  # population variance on N-1 degrees of freedom, N/(N-1) times larger.
  # Without the adjustment one occasion would disagree with prec_prop().
  for (cfg in list(list(p = 0.5, N = 5, n = 2), list(p = 0.3, N = 800, n = 200),
                   list(p = 0.12, N = 50000, n = 4000))) {
    pc <- prec_change(p = rep(cfg$p, 2L), n = cfg$n, N = cfg$N)
    pp <- prec_prop(p = cfg$p, n = cfg$n, N = cfg$N)
    expect_equal(pc$se^2 / 2, pp$se^2, tolerance = 1e-12)
  }
})

test_that("at infinite N the two scales are unadjusted", {
  pc <- prec_change(p = c(0.3, 0.3), n = 400)
  pp <- prec_prop(p = 0.3, n = 400)
  expect_equal(pc$se^2 / 2, pp$se^2, tolerance = 1e-12)
  expect_equal(pc$params$var, c(0.21, 0.21), tolerance = 1e-12)
})

test_that("the adjustment is per occasion when the two N differ", {
  res <- prec_change(p = c(0.3, 0.4), n = 100, N = c(500, 2000))
  expect_equal(
    res$params$var,
    c(0.3 * 0.7 * 500 / 499, 0.4 * 0.6 * 2000 / 1999),
    tolerance = 1e-12
  )
})

## P9. The correlation two proportions can actually have

test_that("overlap_cor is bounded by the Bernoulli marginals", {
  # max corr = (min(p1, p2) - p1 p2) / sqrt(p1 q1 p2 q2)
  bound <- (0.30 - 0.30 * 0.36) / sqrt(0.30 * 0.70 * 0.36 * 0.64)
  expect_error(
    prec_change(p = c(0.30, 0.36), n = 900, overlap = 0.5,
                overlap_cor = bound + 0.01),
    "exceeds the largest correlation"
  )
  expect_silent(
    prec_change(p = c(0.30, 0.36), n = 900, overlap = 0.5,
                overlap_cor = bound - 0.01)
  )
})

test_that("the bound is not checked where the correlation is unused", {
  expect_silent(prec_change(p = c(0.30, 0.36), n = 900, overlap_cor = 1))
  expect_silent(prec_change(p = c(0.30, 0.36), n = 900, overlap = 0.5,
                            overlap_cor = 0))
})

test_that("the mean scale has no such bound", {
  expect_silent(prec_change(var = 100, n = 900, overlap = 0.5,
                            overlap_cor = 1))
})

## P10. The variance is piecewise in overlap

test_that("overlap 0 drops the whole covariance, population term included", {
  # Not the same as evaluating the overlap > 0 expression at overlap = 0,
  # which would keep +2 rho sqrt(v1 v2) / N.
  res <- prec_change(var = 100, n = 500, N = 10000, overlap = 0,
                     overlap_cor = 0.8)
  independent <- 100 / 500 + 100 / 500 - 100 / 10000 - 100 / 10000
  expect_equal(res$se^2, independent, tolerance = 1e-12)
  expect_equal(res$se, prec_change(var = 100, n = 500, N = 10000)$se,
               tolerance = 1e-12)
})

test_that("print names the absent covariance without calling it disjoint", {
  out <- capture.output(print(prec_change(var = 100, n = 500, overlap = 1,
                                          overlap_cor = 0)))
  expect_true(any(grepl("No between-occasion covariance", out)))
  expect_false(any(grepl("independent", out)))
})
