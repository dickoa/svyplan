## The row classification is shared, so each defect is checked on the simple
## and the multistage path.

## L_R for r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7.
RV <- 0.646

test_that("an incidental mu does not rescale a ratio row's rmoe", {
  # A mixed table carries 'mu' for its mean rows, and taking the first
  # non-missing of p, mu, r gave the ratio row that scale.
  bare <- n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                             component_cor = 0.7, rmoe = 0.1))
  with_mu <- n_multi(data.frame(r = 2, mu = 100, cv_num = 1.1, cv_den = 0.6,
                                component_cor = 0.7, rmoe = 0.1))

  expect_equal(with_mu$n, bare$n, tolerance = 1e-10)
  expect_equal(bare$n, RV / (0.1 / qnorm(0.975))^2, tolerance = 1e-8)
})

test_that("an incidental mu does not rescale a clustered ratio row", {
  bare <- n_cluster(
    indicators = data.frame(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
               rmoe = 0.1, icc_psu = 0.05),
    stage_cost = c(500, 50)
  )
  with_mu <- n_cluster(
    indicators = data.frame(r = 2, mu = 100, cv_num = 1.1, cv_den = 0.6,
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
  # Sizing recomputed it from the moments while achieved precision read the
  # supplied value, so one result reported two designs.
  expect_error(
    n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7, cv = 0.05, unit_relvar = 100)),
    "but its moments imply"
  )
  expect_error(
    n_cluster(indicators = data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                               component_cor = 0.7, cv = 0.05,
                               icc_psu = 0.05, unit_relvar = 100),
                    stage_cost = c(500, 50)),
    "but its moments imply"
  )
})

test_that("a unit_relvar agreeing with the moments is accepted", {
  # A result stores this, so rejecting on presence broke every round trip.
  expect_s3_class(
    n_multi(data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7, cv = 0.05, unit_relvar = RV)),
    "svyplan_n"
  )
})

test_that("a non-ratio row may still supply unit_relvar", {
  expect_s3_class(
    n_cluster(indicators = data.frame(var = 100, mu = 50, cv = 0.05, icc_psu = 0.05,
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
    prec_cluster(indicators = cbind(two, n = 50, n_per_psu = 10, icc_psu = 0.05),
                       stage_cost = c(500, 50)),
    "only one of 'p', 'var', or 'r'"
  )
  expect_error(
    prec_multi(data.frame(var = 100, mu = 50, r = 2, cv_num = 1.1,
                          cv_den = 0.6, component_cor = 0.7, n = 500)),
    "only one of 'p', 'var', or 'r'"
  )
})

test_that("prec_cluster(indicators = ) validates ratio moments", {
  base <- list(n = 50, n_per_psu = 10, icc_psu = 0.05)
  expect_error(
    do.call(prec_cluster,
            list(indicators = data.frame(r = 0, cv_num = 1.1, cv_den = 0.6,
                                         component_cor = 0.7, n = 50,
                                         n_per_psu = 10, icc_psu = 0.05),
                 stage_cost = c(500, 50))),
    "must not be zero"
  )
  expect_error(
    prec_cluster(indicators = data.frame(r = 2, cv_num = 1.1, n = 50,
                                  n_per_psu = 10, icc_psu = 0.05),
                       stage_cost = c(500, 50)),
    "column\\(s\\) .*'cv_den'.*are absent"
  )
  expect_error(
    prec_cluster(indicators = data.frame(r = 2, cv_num = 1.1, cv_den = 0.6,
                                  component_cor = 5, n = 50, n_per_psu = 10,
                                  icc_psu = 0.05),
                       stage_cost = c(500, 50)),
    "must be a number in \\[-1, 1\\]"
  )
})

test_that("an extreme but finite ratio keeps its standard error", {
  # Forming r^2 first overflows at 1e200 and underflows at 1e-200, both
  # destroying a representable answer.
  for (r in c(1e200, 1e-200)) {
    res <- prec_ratio(r = r, n = 500, cv_num = 1.1, cv_den = 0.6,
                      component_cor = 0.7)
    expect_equal(res$se, abs(r) * sqrt(RV / 500), tolerance = 1e-10)
    expect_true(is.finite(res$se))
    expect_gt(res$se, 0)
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
  # The tolerance is built from the same squared terms, so an overflow made
  # it infinite and clamped every value to zero.
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
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05, df = 7)
  ci <- confint(res)
  half <- (ci[2L] - ci[1L]) / 2

  expect_equal(half, qt(0.975, 7) * res$se, tolerance = 1e-10)
  expect_gt(half, qnorm(0.975) * res$se)
})

test_that("a negative ratio is precise at a strongly negative correlation", {
  # The quantity is component_cor * sign(r), not component_cor.
  negative_cor <- .ratio_unit_relvar(-2, 1.1, 0.6, -0.9)
  positive_cor <- .ratio_unit_relvar(-2, 1.1, 0.6, 0.9)

  expect_lt(negative_cor, positive_cor)
  expect_equal(negative_cor, .ratio_unit_relvar(2, 1.1, 0.6, 0.9))
})

## Round trips through object dispatch, not a hand-rebuilt indicator table.
## Rebuilding by hand is what let the stored unit_relvar go unnoticed.

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
  sized <- n_cluster(
    indicators = data.frame(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
               cv = 0.05, icc_psu = 0.05),
    stage_cost = c(500, 50)
  )
  prec <- prec_cluster(sized)
  back <- n_cluster(prec)

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
  # Carried on the way back rather than blanked, so a ratio row is not
  # asymmetric with a mean row, which keeps its own.
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

test_that("multistage cv mode scales .se by the estimand of every row type", {
  # Every estimand fills .se the same way, off its own scale. Checked across
  # all three so the ratio row is never given a branch of its own.
  args <- list(n = 45, n_per_psu = 14, icc_psu = 0.05)
  rows <- list(
    prop = data.frame(p = 0.3),
    mean = data.frame(var = 2500, mu = 300),
    ratio = data.frame(r = 420, cv_num = 1.2, cv_den = 0.45,
                       component_cor = 0.65)
  )
  estimand <- c(prop = 0.3, mean = 300, ratio = 420)
  for (nm in names(rows)) {
    res <- prec_cluster(indicators = cbind(rows[[nm]], as.data.frame(args)),
                              stage_cost = c(500, 50))
    expect_true(is.finite(res$detail$.cv), info = nm)
    expect_true(is.finite(res$detail$.se), info = nm)
    expect_equal(res$detail$.se, res$detail$.cv * estimand[[nm]],
                 tolerance = 1e-12, info = nm)
    expect_equal(res$detail$.moe, qnorm(0.975) * res$detail$.se,
                 tolerance = 1e-10, info = nm)
  }
})
