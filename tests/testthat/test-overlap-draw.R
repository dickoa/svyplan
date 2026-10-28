## How the shared units of two occasions are drawn
##
## "fixed" draws the shared units as a fixed count with the rest of each
## occasion's sample. "independent" draws a retained panel once and tops
## each occasion up with a fresh sample from the whole frame. The exact
## variances are checked by enumeration in test-finite-population.R. These
## tests cover the boundaries, the frame check, the families that share the
## kernel, and the round trips.

two_sided_power <- function(effect, se, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  pnorm(abs(effect) / se - z) + pnorm(-abs(effect) / se - z)
}

## T1. Where the two draws agree, and where they part

test_that("the default is the fixed draw", {
  args <- list(var = 100, n = 40, N = 400, overlap = 0.5, overlap_cor = 0.8)
  expect_identical(
    do.call(prec_change, args),
    do.call(prec_change, c(args, overlap_draw = "fixed"))
  )
  expect_identical(
    do.call(prec_change, args)$params$overlap_draw, "fixed"
  )
})

test_that("the draws agree at full overlap and in an infinite population", {
  full <- list(var = c(100, 120), n = 40, N = 400, overlap = 1,
               overlap_cor = 0.7)
  expect_equal(
    do.call(prec_change, c(full, overlap_draw = "independent"))$se,
    do.call(prec_change, full)$se,
    tolerance = 1e-12
  )
  inf <- list(var = c(100, 120), n = 40, overlap = 0.5, overlap_cor = 0.7)
  expect_equal(
    do.call(prec_change, c(inf, overlap_draw = "independent"))$se,
    do.call(prec_change, inf)$se,
    tolerance = 1e-12
  )
})

test_that("only the independent draw is continuous at zero overlap", {
  base <- prec_change(var = 100, n = 40, N = 400)$se
  near <- function(draw) {
    prec_change(var = 100, n = 40, N = 400, overlap = 1e-9,
                overlap_cor = 0.5, overlap_draw = draw)$se
  }
  expect_equal(near("independent"), base, tolerance = 1e-8)
  expect_gt(abs(near("fixed") / base - 1), 0.01)
})

test_that("the two-occasion average parts from the fixed draw as documented", {
  gap <- function(f) {
    args <- list(var = 1, n = 100, occasions = 2, N = 100 / f,
                 overlap = 0.5, overlap_cor = 0.8)
    indep <- do.call(prec_pooled, c(args, overlap_draw = "independent"))
    100 * (indep$se^2 / do.call(prec_pooled, args)$se^2 - 1)
  }
  expect_equal(round(gap(0.05), 1), 4.2)
  expect_equal(round(gap(0.10), 1), 9.0)
})

test_that("the change parts from the fixed draw as documented", {
  gap <- function(f) {
    args <- list(var = 1, n = 100, N = 100 / f, overlap = 0.5,
                 overlap_cor = 0.8)
    indep <- do.call(prec_change, c(args, overlap_draw = "independent"))
    100 * (indep$se^2 / do.call(prec_change, args)$se^2 - 1)
  }
  expect_equal(round(gap(0.05), 1), -0.8)
  expect_equal(round(gap(0.10), 1), -1.7)
})

## T2. Each part fits the frame on its own

test_that("independent draws lift the distinct-units bound", {
  expect_error(
    prec_change(var = 100, n = 60, N = 100, overlap = 0.2,
                overlap_cor = 0.5),
    "108 distinct units"
  )
  expect_no_error(
    prec_change(var = 100, n = 60, N = 100, overlap = 0.2,
                overlap_cor = 0.5, overlap_draw = "independent")
  )
  expect_error(
    prec_pooled(var = 100, n = 60, occasions = 2, N = 100, overlap = 0.2,
                overlap_cor = 0.5),
    "distinct units"
  )
  expect_no_error(
    prec_pooled(var = 100, n = 60, occasions = 2, N = 100, overlap = 0.2,
                overlap_cor = 0.5, overlap_draw = "independent")
  )
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.4, n = 60, N = 100, power = NULL,
               overlap = 0.2, overlap_cor = 0.3),
    "distinct units"
  )
  expect_no_error(
    power_prop(p1 = 0.3, p2 = 0.4, n = 60, N = 100, power = NULL,
               overlap = 0.2, overlap_cor = 0.3,
               overlap_draw = "independent")
  )
  expect_error(
    power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32), outcome = "prop",
              n = 60, N = 100, power = NULL, overlap = 0.2,
              overlap_cor = 0.3),
    "distinct units"
  )
  expect_no_error(
    power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32), outcome = "prop",
              n = 60, N = 100, power = NULL, overlap = 0.2,
              overlap_cor = 0.3, overlap_draw = "independent")
  )
})

test_that("a power solve under independent draws reaches past the fixed cap", {
  args <- list(var = 100, effect = 3, N = 100, overlap = 0.2,
               overlap_cor = 0.5)
  target <- do.call(power_mean, c(args, n = 70, power = list(NULL),
                                  overlap_draw = "independent"))$power
  solved <- do.call(power_mean, c(args, power = target,
                                  overlap_draw = "independent"))
  expect_equal(solved$n, 70, tolerance = 1e-6)
  expect_error(do.call(power_mean, c(args, power = target)), "unattainable")
})

## T3. The families that share the kernel

test_that("power_mean's MDE at power 0.5 tracks prec_change under both draws", {
  for (draw in c("fixed", "independent")) {
    args <- list(var = 100, n = 40, N = 400, overlap = 0.5,
                 overlap_cor = 0.78, overlap_draw = draw)
    moe <- do.call(prec_change, args)$moe
    mde <- do.call(power_mean, c(args, power = 0.5))$effect
    expect_equal(mde, moe, tolerance = 1e-4, info = draw)
  }
  indep <- power_mean(var = 100, n = 40, N = 400, power = 0.5, overlap = 0.5,
                      overlap_cor = 0.78, overlap_draw = "independent")
  fixed <- power_mean(var = 100, n = 40, N = 400, power = 0.5, overlap = 0.5,
                      overlap_cor = 0.78)
  expect_gt(abs(indep$effect / fixed$effect - 1), 1e-3)
})

test_that("power_prop under independent draws is the change's normal power", {
  args <- list(n = 60, N = 300, overlap = 0.5, overlap_cor = 0.4,
               overlap_draw = "independent")
  se <- do.call(prec_change, c(args, list(p = c(0.3, 0.4))))$se
  pw <- do.call(power_prop, c(args, p1 = 0.3, p2 = 0.4,
                              power = list(NULL)))$power
  expect_equal(pw, two_sided_power(0.1, se), tolerance = 1e-12)
})

test_that("power_did under independent draws sums the two arms' changes", {
  v <- c(100, 120, 90, 110)
  arm <- function(var) {
    prec_change(var = var, n = 80, N = 500, overlap = 0.5, overlap_cor = 0.6,
                overlap_draw = "independent")$se^2
  }
  se <- sqrt(arm(v[1:2]) + arm(v[3:4]))
  pw <- power_did(treat = c(10, 14), control = c(10, 11), outcome = "mean",
                  var = v, n = 80, N = 500, power = NULL, overlap = 0.5,
                  overlap_cor = 0.6, overlap_draw = "independent")$power
  expect_equal(pw, two_sided_power(3, se), tolerance = 1e-12)
})

test_that("power_did draws agree at full overlap", {
  args <- list(treat = c(0.3, 0.4), control = c(0.3, 0.32), outcome = "prop",
               n = 100, N = 400, power = NULL, overlap = 1, overlap_cor = 0.6)
  expect_equal(
    do.call(power_did, c(args, overlap_draw = "independent"))$power,
    do.call(power_did, args)$power,
    tolerance = 1e-12
  )
})

## T4. Unequal occasions

test_that("unequal occasions use each occasion's retained share", {
  res <- n_change(var = 100, moe = 4, N = 400, ratio = 2, overlap = 0.4,
                  overlap_cor = 0.6, overlap_draw = "independent")
  n <- res$n * c(1, 1)
  w <- c(0.4, 0.4 * n[1] / n[2])
  cross <- 0.6 * 100
  hand <- sum(100 * (1 / n - (w^2 + (1 - w)^2) / 400)) -
    2 * cross * (0.4 / n[2] - w[1] * w[2] / 400)
  expect_equal(res$se^2, hand, tolerance = 1e-12)
  expect_equal(res$moe, 4, tolerance = 1e-9)
})

## T5. Pooled estimates beyond two occasions

test_that("independent draws pool two occasions only when units are shared", {
  expect_error(
    prec_pooled(var = 100, n = 40, occasions = 3, N = 400, overlap = 0.5,
                overlap_cor = 0.5, overlap_draw = "independent"),
    "two occasions only"
  )
  expect_error(
    n_pooled(var = 100, moe = 3, occasions = 4, N = 400,
             overlap = design_overlap("2-2-2"), cor_decay = 0.5,
             overlap_draw = "independent"),
    "two occasions only"
  )
  fresh <- list(var = 100, n = 40, occasions = 3, N = 400)
  expect_equal(
    do.call(prec_pooled, c(fresh, overlap_draw = "independent"))$se,
    do.call(prec_pooled, fresh)$se,
    tolerance = 1e-12
  )
})

## T6. Round trips, predict, plan profiles

test_that("round trips keep the independent draw", {
  res <- n_change(var = 100, moe = 4, N = 400, overlap = 0.5,
                  overlap_cor = 0.8, overlap_draw = "independent")
  back <- prec_change(res)
  expect_identical(back$params$overlap_draw, "independent")
  expect_equal(back$moe, 4, tolerance = 1e-9)
  expect_equal(n_change(back)$n, res$n, tolerance = 1e-9)

  pool <- n_pooled(var = 100, moe = 2.5, occasions = 2, N = 400,
                   overlap = 0.5, overlap_cor = 0.8,
                   overlap_draw = "independent")
  pback <- prec_pooled(pool)
  expect_identical(pback$params$overlap_draw, "independent")
  expect_equal(pback$moe, 2.5, tolerance = 1e-9)
  expect_equal(n_pooled(pback)$n, pool$n, tolerance = 1e-9)
})

test_that("predict reproduces direct calls under the independent draw", {
  res <- n_change(var = 100, moe = 4, N = 400, overlap = 0.5,
                  overlap_cor = 0.8, overlap_draw = "independent")
  got <- predict(res, data.frame(N = c(300, 800)))
  for (i in 1:2) {
    direct <- n_change(var = 100, moe = 4, N = got$N[i], overlap = 0.5,
                       overlap_cor = 0.8, overlap_draw = "independent")
    expect_equal(got$n1[i], direct$n, tolerance = 1e-12)
  }

  pool <- n_pooled(var = 100, moe = 2.5, occasions = 2, N = 400,
                   overlap = 0.5, overlap_cor = 0.8,
                   overlap_draw = "independent")
  pgot <- predict(pool, data.frame(N = 800))
  expect_equal(
    pgot$n,
    n_pooled(var = 100, moe = 2.5, occasions = 2, N = 800, overlap = 0.5,
             overlap_cor = 0.8, overlap_draw = "independent")$n,
    tolerance = 1e-12
  )

  pm <- power_mean(var = 100, effect = 4, N = 400, overlap = 0.5,
                   overlap_cor = 0.8, overlap_draw = "independent")
  expect_equal(
    predict(pm, data.frame(N = 800))$n,
    power_mean(var = 100, effect = 4, N = 800, overlap = 0.5,
               overlap_cor = 0.8, overlap_draw = "independent")$n,
    tolerance = 1e-12
  )
  pp <- power_prop(p1 = 0.3, p2 = 0.4, N = 600, overlap = 0.5,
                   overlap_cor = 0.4, overlap_draw = "independent")
  expect_equal(
    predict(pp, data.frame(N = 900))$n,
    power_prop(p1 = 0.3, p2 = 0.4, N = 900, overlap = 0.5,
               overlap_cor = 0.4, overlap_draw = "independent")$n,
    tolerance = 1e-12
  )
  pd <- power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                  outcome = "prop", n = 100, N = 400, power = NULL,
                  overlap = 0.5, overlap_cor = 0.6,
                  overlap_draw = "independent")
  expect_equal(
    predict(pd, data.frame(N = 800))$power,
    power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32), outcome = "prop",
              n = 100, N = 800, power = NULL, overlap = 0.5,
              overlap_cor = 0.6, overlap_draw = "independent")$power,
    tolerance = 1e-12
  )
})

test_that("a plan profile carries the draw and validates it", {
  plan <- svyplan(overlap_draw = "independent")
  expect_equal(
    prec_change(var = 100, n = 40, N = 400, overlap = 0.5,
                overlap_cor = 0.8, plan = plan)$se,
    prec_change(var = 100, n = 40, N = 400, overlap = 0.5,
                overlap_cor = 0.8, overlap_draw = "independent")$se,
    tolerance = 1e-12
  )
  expect_error(svyplan(overlap_draw = "random"), "'fixed' or 'independent'")
  expect_error(
    prec_change(var = 100, n = 40, overlap_draw = "random"),
    "should be one of"
  )
})

## T7. Printing

test_that("print names the independent draw only where it moves the result", {
  shows <- function(x) any(grepl("independent draws", capture.output(print(x))))
  args <- list(var = 100, n = 40, overlap = 0.5, overlap_cor = 0.8)
  expect_true(shows(do.call(prec_change,
                            c(args, N = 400, overlap_draw = "independent"))))
  expect_false(shows(do.call(prec_change, c(args, N = 400))))
  expect_false(shows(do.call(prec_change,
                             c(args, overlap_draw = "independent"))))
  expect_false(shows(prec_change(var = 100, n = 40, N = 400,
                                 overlap_draw = "independent")))
  expect_true(shows(prec_pooled(var = 100, n = 40, occasions = 2, N = 400,
                                overlap = 0.5, overlap_cor = 0.8,
                                overlap_draw = "independent")))
  expect_true(shows(power_mean(var = 100, effect = 4, N = 400, overlap = 0.5,
                               overlap_cor = 0.8,
                               overlap_draw = "independent")))
})
