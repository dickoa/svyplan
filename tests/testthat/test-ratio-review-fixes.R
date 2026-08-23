## Regressions for six defects found in review after the ratio work landed.
## Each is checked on the simple and the multistage path, because the row
## classification is shared and five of the six reached both.

## L_R for r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7.
RV <- 0.646

test_that("an incidental mu does not rescale a ratio row's rmoe", {
  # A mixed table legitimately carries 'mu' for its mean rows. Selecting the
  # scale as the first non-missing of p, mu, r gave the ratio row the mean's
  # scale and undersized it by a factor of 50.
  bare <- n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                             component_cor = 0.7, rmoe = 0.1))
  with_mu <- n_multi(data.frame(r = 2, mu = 100, cv_num = 1.1, cv_den = 0.6,
                                component_cor = 0.7, rmoe = 0.1))

  expect_equal(with_mu$n, bare$n, tolerance = 1e-10)
  expect_equal(bare$n, RV / (0.1 / qnorm(0.975))^2, tolerance = 1e-8)
})

test_that("an incidental mu does not rescale a clustered ratio row", {
  bare <- n_multi_cluster(
    data.frame(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
               rmoe = 0.1, icc_psu = 0.05),
    stage_cost = c(500, 50)
  )
  with_mu <- n_multi_cluster(
    data.frame(r = 2, mu = 100, cv_num = 1.1, cv_den = 0.6,
               component_cor = 0.7, rmoe = 0.1, icc_psu = 0.05),
    stage_cost = c(500, 50)
  )

  expect_equal(with_mu$total_n, bare$total_n, tolerance = 1e-8)
})

test_that("a mean row still takes its scale from mu", {
  by_rmoe <- n_multi(data.frame(var = 100, mu = 50, rmoe = 0.1))
  by_moe <- n_multi(data.frame(var = 100, mu = 50, moe = 5))
  expect_equal(by_rmoe$n, by_moe$n, tolerance = 1e-10)
})

test_that("a proportion row still takes its scale from p", {
  by_rmoe <- n_multi(data.frame(p = 0.3, rmoe = 0.1))
  by_moe <- n_multi(data.frame(p = 0.3, moe = 0.03))
  expect_equal(by_rmoe$n, by_moe$n, tolerance = 1e-10)
})

test_that("a unit_relvar contradicting the ratio moments is refused", {
  # Sizing recomputed the coefficient from the moments while achieved
  # precision read the supplied value, so one result reported two designs.
  expect_error(
    n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7, cv = 0.05, unit_relvar = 100)),
    "but its moments imply"
  )
  expect_error(
    n_multi_cluster(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                               component_cor = 0.7, cv = 0.05,
                               icc_psu = 0.05, unit_relvar = 100),
                    stage_cost = c(500, 50)),
    "but its moments imply"
  )
})

test_that("a unit_relvar agreeing with the moments is accepted", {
  # This is what a result stores, so rejecting on presence rather than on
  # value broke every S3 round trip through a ratio row.
  expect_s3_class(
    n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7, cv = 0.05, unit_relvar = RV)),
    "svyplan_n"
  )
})

test_that("a non-ratio row may still supply unit_relvar", {
  expect_s3_class(
    n_multi_cluster(data.frame(var = 100, mu = 50, cv = 0.05, icc_psu = 0.05,
                               unit_relvar = 0.04),
                    stage_cost = c(500, 50)),
    "svyplan_cluster"
  )
})

test_that("the precision paths require exactly one estimand per row", {
  two <- data.frame(p = 0.3, r = 2, cv_num = 1.1, cv_den = 0.6,
                    component_cor = 0.7)
  expect_error(
    prec_multi(cbind(two, n = 500)),
    "only one of 'p', 'var', or 'r'"
  )
  expect_error(
    prec_multi_cluster(cbind(two, n = 50, n_per_psu = 10, icc_psu = 0.05),
                       stage_cost = c(500, 50)),
    "only one of 'p', 'var', or 'r'"
  )
  expect_error(
    prec_multi(data.frame(var = 100, mu = 50, r = 2, cv_num = 1.1,
                          cv_den = 0.6, component_cor = 0.7, n = 500)),
    "only one of 'p', 'var', or 'r'"
  )
})

test_that("prec_multi_cluster validates ratio moments", {
  base <- list(n = 50, n_per_psu = 10, icc_psu = 0.05)
  expect_error(
    do.call(prec_multi_cluster,
            list(data.frame(r = 0, cv_num = 1.1, cv_den = 0.6,
                            component_cor = 0.7, n = 50, n_per_psu = 10,
                            icc_psu = 0.05),
                 stage_cost = c(500, 50))),
    "must not be zero"
  )
  expect_error(
    prec_multi_cluster(data.frame(r = 2, cv_num = 1.1, n = 50,
                                  n_per_psu = 10, icc_psu = 0.05),
                       stage_cost = c(500, 50)),
    "column\\(s\\) .*'cv_den'.*are absent"
  )
  expect_error(
    prec_multi_cluster(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                                  component_cor = 5, n = 50, n_per_psu = 10,
                                  icc_psu = 0.05),
                       stage_cost = c(500, 50)),
    "must be a number in \\[-1, 1\\]"
  )
})

test_that("an extreme but finite ratio keeps its standard error", {
  # se = abs(r) * sqrt(L_R / n). Forming r^2 first overflows at r = 1e200 and
  # underflows at r = 1e-200, in both cases destroying a representable answer.
  for (r in c(1e200, 1e-200)) {
    res <- prec_ratio(r = r, n = 500, cv_num = 1.1, cv_den = 0.6,
                      component_cor = 0.7)
    expect_equal(res$se, abs(r) * sqrt(RV / 500), tolerance = 1e-10)
    expect_true(is.finite(res$se))
    expect_gt(res$se, 0)
    # The relative measure is free of the scale entirely.
    expect_equal(res$cv, sqrt(RV / 500), tolerance = 1e-12)
  }
})

test_that("an extreme ratio sizes as its relative precision demands", {
  small <- n_ratio(r = 1e-200, cv_num = 1.1, cv_den = 0.6,
                   component_cor = 0.7, cv = 0.05)
  ordinary <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6,
                      component_cor = 0.7, cv = 0.05)
  expect_equal(small$n, ordinary$n, tolerance = 1e-10)
})

test_that("a coefficient of variation too large to square is refused", {
  # The cancellation tolerance is built from the same squared terms, so an
  # overflow there made it infinite and clamped every value to zero.
  expect_error(
    prec_ratio(r = 2, n = 500, cv_num = 1e200, cv_den = 0.6,
               component_cor = 0.7),
    "not finite at these moments"
  )
  expect_error(
    n_ratio(r = 2, cv_num = 1e200, cv_den = 0.6, component_cor = 0.7,
            cv = 0.05),
    "not finite at these moments"
  )
})

test_that("the clamp still zeroes genuine cancellation", {
  expect_identical(.ratio_unit_relvar(2, 0.8, 0.8, 1), 0)
  expect_gt(.ratio_unit_relvar(2, 1, 1, 1 - 1e-9), 0)
})

test_that("confint on a ratio uses the t quantile when df is set", {
  # The shared help page said means use a z interval even when df selects t.
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05, df = 7)
  ci <- confint(res)
  half <- (ci[2L] - ci[1L]) / 2

  expect_equal(half, qt(0.975, 7) * res$se, tolerance = 1e-10)
  expect_gt(half, qnorm(0.975) * res$se)
})

test_that("a negative ratio is precise at a strongly negative correlation", {
  # The help said a high correlation makes a ratio precise, which holds only
  # for a positive ratio. The quantity is component_cor * sign(r).
  negative_cor <- .ratio_unit_relvar(-2, 1.1, 0.6, -0.9)
  positive_cor <- .ratio_unit_relvar(-2, 1.1, 0.6, 0.9)

  expect_lt(negative_cor, positive_cor)
  expect_equal(negative_cor, .ratio_unit_relvar(2, 1.1, 0.6, 0.9))
})

## The S3 round trips, exercised through object dispatch rather than a
## hand-rebuilt indicator table. Rebuilding by hand is what let the stored
## unit_relvar go unnoticed: a result stores the derived coefficient, and
## handing that table straight back looked like a user supplying it.

test_that("a ratio result round-trips through prec_multi() and back", {
  sized <- n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                              component_cor = 0.7, cv = 0.05))
  prec <- prec_multi(sized)
  back <- n_multi(prec)

  expect_s3_class(prec, "svyplan_prec")
  expect_equal(prec$detail$.cv, 0.05, tolerance = 1e-10)
  expect_equal(back$n, sized$n, tolerance = 1e-10)
})

test_that("a clustered ratio result round-trips in both directions", {
  sized <- n_multi_cluster(
    data.frame(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
               cv = 0.05, icc_psu = 0.05),
    stage_cost = c(500, 50)
  )
  prec <- prec_multi_cluster(sized)
  back <- n_multi_cluster(prec)

  expect_equal(prec$detail$.cv, 0.05, tolerance = 1e-6)
  expect_equal(unname(back$n), unname(sized$n), tolerance = 1e-6)
})

test_that("a mixed table round-trips with its ratio row intact", {
  tbl <- data.frame(
    name = c("literacy", "consumption"),
    p = c(0.35, NA), r = c(NA, 420),
    cv_num = c(NA, 1.20), cv_den = c(NA, 0.45),
    component_cor = c(NA, 0.65), cv = 0.05
  )
  sized <- n_multi(tbl)
  back <- n_multi(prec_multi(sized))

  expect_equal(back$n, sized$n, tolerance = 1e-10)
  expect_equal(back$indicators$r, sized$indicators$r)
})

test_that("the derived coefficient is stored and survives the round trip", {
  # It is derived from the moments on the way in and carried on the way back,
  # rather than blanked, so a ratio row is not asymmetric with a mean row,
  # which keeps its own. Blanking it was an abandoned fix for the round-trip
  # regression; consistency-checking the value is what shipped.
  sized <- n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                              component_cor = 0.7, cv = 0.05))
  expect_equal(sized$indicators$unit_relvar, RV)

  prec <- prec_multi(sized)
  expect_equal(prec$params$indicators$unit_relvar, RV)
})

test_that("a contradictory unit_relvar is still refused after the fix", {
  expect_error(
    n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7, cv = 0.05, unit_relvar = 100)),
    "but its moments imply"
  )
})

test_that("multistage cv mode reports .cv alone for every estimand", {
  # Not a ratio limitation: se and moe are NA in cv mode for proportions and
  # means too, which test-n_multi.R pins by name. Recorded here so the ratio
  # row is not later "fixed" on its own.
  args <- list(n = 45, n_per_psu = 14, icc_psu = 0.05)
  rows <- list(
    prop = data.frame(p = 0.3),
    mean = data.frame(var = 2500, mu = 300),
    ratio = data.frame(r = 420, cv_num = 1.2, cv_den = 0.45,
                       component_cor = 0.65)
  )
  for (nm in names(rows)) {
    res <- prec_multi_cluster(cbind(rows[[nm]], as.data.frame(args)),
                              stage_cost = c(500, 50))
    expect_true(is.na(res$detail$.se), info = nm)
    expect_true(is.finite(res$detail$.cv), info = nm)
  }
})
