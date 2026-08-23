## Ratio rows in the multi-indicator grammar. The reference for every ratio
## quantity is the scalar pair, which is itself checked against
## dev/ratio-probe.R, so these assert agreement rather than re-deriving.

ratio_row <- function(...) {
  base <- data.frame(r = 420, cv_num = 1.20, cv_den = 0.45,
                     component_cor = 0.65)
  extra <- list(...)
  for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  base
}

test_that("a one-row ratio table matches n_ratio()", {
  for (target in list(list(cv = 0.05), list(moe = 20), list(rmoe = 0.05))) {
    tbl <- do.call(ratio_row, target)
    solo <- do.call(
      n_ratio,
      c(list(r = 420, cv_num = 1.20, cv_den = 0.45, component_cor = 0.65),
        target)
    )
    expect_equal(n_multi(tbl)$n, solo$n, tolerance = 1e-10,
                 info = names(target))
  }
})

test_that("a one-row ratio precision table matches prec_ratio()", {
  tbl <- ratio_row(n = 1200)
  solo <- prec_ratio(r = 420, n = 1200, cv_num = 1.20, cv_den = 0.45,
                     component_cor = 0.65)
  res <- prec_multi(tbl)

  expect_equal(res$se, solo$se, tolerance = 1e-10)
  expect_equal(res$moe, solo$moe, tolerance = 1e-10)
  expect_equal(res$cv, solo$cv, tolerance = 1e-10)
  expect_equal(res$detail$.rmoe, solo$rmoe, tolerance = 1e-10)
})

test_that("design columns reach a ratio row", {
  tbl <- ratio_row(cv = 0.05, N = 9000, deff = 1.4, resp_rate = 0.8,
                   alpha = 0.10)
  solo <- n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45,
                  component_cor = 0.65, cv = 0.05, N = 9000, deff = 1.4,
                  resp_rate = 0.8, alpha = 0.10)

  expect_equal(n_multi(tbl)$n, solo$n, tolerance = 1e-10)
})

test_that("a per-row df reaches a ratio row", {
  tbl <- ratio_row(moe = 20, df = 8)
  solo <- n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45,
                  component_cor = 0.65, moe = 20, df = 8)

  expect_equal(n_multi(tbl)$n, solo$n, tolerance = 1e-10)
  expect_gt(n_multi(tbl)$n, n_multi(ratio_row(moe = 20))$n)
})

test_that("a mixed table sizes proportions, means, and ratios together", {
  tbl <- data.frame(
    name = c("literacy", "income", "consumption"),
    p = c(0.35, NA, NA),
    var = c(NA, 2500, NA),
    mu = c(NA, 300, NA),
    r = c(NA, NA, 420),
    cv_num = c(NA, NA, 1.20),
    cv_den = c(NA, NA, 0.45),
    component_cor = c(NA, NA, 0.65),
    cv = 0.05
  )
  res <- n_multi(tbl)

  expect_equal(nrow(res$detail), 3L)
  expect_equal(res$binding, "literacy")
  expect_equal(
    res$detail$.n[3L],
    n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45, component_cor = 0.65,
            cv = 0.05)$n,
    tolerance = 1e-10
  )
  # Every row reports the CV it achieves at the winning size.
  expect_equal(res$detail$.cv_achieved[1L], 0.05, tolerance = 1e-10)
  expect_true(all(res$detail$.cv_achieved[2:3] < 0.05))
})

test_that("a ratio row can be the binding indicator", {
  tbl <- data.frame(
    name = c("literacy", "consumption"),
    p = c(0.35, NA),
    r = c(NA, 420),
    cv_num = c(NA, 1.20),
    cv_den = c(NA, 0.45),
    component_cor = c(NA, 0.65),
    cv = c(0.20, 0.02)
  )
  res <- n_multi(tbl)

  expect_equal(res$binding, "consumption")
  expect_true(res$detail$.binding[2L])
})

test_that("rmoe on a ratio row is relative to r", {
  by_rmoe <- n_multi(ratio_row(rmoe = 0.05))
  by_moe <- n_multi(ratio_row(moe = 0.05 * 420))

  expect_equal(by_rmoe$n, by_moe$n, tolerance = 1e-10)
})

test_that("the reported rmoe of a ratio row uses r, not a mean", {
  res <- prec_multi(ratio_row(n = 1200))
  expect_equal(res$detail$.rmoe, res$detail$.moe / 420, tolerance = 1e-10)
})

test_that("the ratio inputs survive a prec_multi round trip", {
  sized <- n_multi(ratio_row(cv = 0.05))
  prec <- prec_multi(ratio_row(n = sized$n))
  back <- n_multi(prec)

  expect_equal(back$n, sized$n, tolerance = 1e-10)
  for (nm in c("r", "cv_num", "cv_den", "component_cor")) {
    expect_equal(back$indicators[[nm]], sized$indicators[[nm]], info = nm)
  }
})

test_that("domain ratios are sized from their own moments", {
  tbl <- data.frame(
    region = c("north", "south"),
    r = c(420, 380),
    cv_num = c(1.20, 0.80),
    cv_den = c(0.45, 0.45),
    component_cor = 0.65,
    cv = 0.05
  )
  res <- n_multi(tbl, domains = "region")

  expect_equal(nrow(res$domains), 2L)
  expect_equal(
    res$domains$.n[1L],
    n_ratio(r = 420, cv_num = 1.20, cv_den = 0.45, component_cor = 0.65,
            cv = 0.05)$n,
    tolerance = 1e-10
  )
  expect_equal(
    res$domains$.n[2L],
    n_ratio(r = 380, cv_num = 0.80, cv_den = 0.45, component_cor = 0.65,
            cv = 0.05)$n,
    tolerance = 1e-10
  )
})

test_that("a CV target does not depend on the ratio's magnitude", {
  # CV(R_hat) does not involve R, so two domains with the same component
  # moments and different ratios need the same size. Documented, not a bug.
  tbl <- data.frame(
    region = c("north", "south"),
    r = c(420, 38),
    cv_num = 1.20, cv_den = 0.45, component_cor = 0.65,
    cv = 0.05
  )
  res <- n_multi(tbl, domains = "region")
  expect_equal(res$domains$.n[1L], res$domains$.n[2L], tolerance = 1e-10)
})

test_that("a row may declare only one estimand", {
  expect_error(
    n_multi(data.frame(p = 0.3, r = 420, cv_num = 1.2, cv_den = 0.45,
                       component_cor = 0.65, cv = 0.05)),
    "only one of 'p', 'var', or 'r'"
  )
  expect_error(
    n_multi(data.frame(var = 100, mu = 50, r = 420, cv_num = 1.2,
                       cv_den = 0.45, component_cor = 0.65, cv = 0.05)),
    "only one of 'p', 'var', or 'r'"
  )
})

test_that("an incomplete ratio quartet is refused", {
  expect_error(
    n_multi(data.frame(r = 420, cv_num = 1.2, cv = 0.05)),
    "column\\(s\\) .*'cv_den'.*are absent"
  )
  expect_error(
    n_multi(data.frame(r = 420, cv_num = NA, cv_den = 0.45,
                       component_cor = 0.65, cv = 0.05)),
    "row 1 sets 'r' but leaves .*cv_num.* missing"
  )
})

test_that("ratio moments are validated row by row", {
  expect_error(
    n_multi(ratio_row(component_cor = 1.5, cv = 0.05)),
    "must be a number in \\[-1, 1\\]"
  )
  expect_error(
    n_multi(data.frame(r = 2, cv_num = 0.8, cv_den = 0.8, component_cor = 1,
                       cv = 0.05)),
    "no sampling variance"
  )
})

test_that("unit_relvar alone does not make a row a ratio", {
  # It is already a valid column for cluster planning of means and
  # proportions, so it cannot mark an estimand scale.
  expect_error(
    n_multi(data.frame(unit_relvar = 0.6, cv = 0.05)),
    "must contain a 'p', 'var', or 'r' column"
  )
})

test_that("min_cases is refused on a ratio row and names it", {
  expect_error(
    n_multi(ratio_row(cv = 0.05, min_cases = 50)),
    "Row\\(s\\) 1 carry a ratio"
  )
})

test_that("the multistage paths take ratio rows", {
  # They refused them until the Step 5 gate passed, which required a one-row
  # ratio table to reproduce the direct n_cluster() call. See
  # test-multi-cluster-ratio.R for that comparison and the stage semantics.
  sized <- n_multi_cluster(ratio_row(cv = 0.05, icc_psu = 0.05),
                           stage_cost = c(500, 50))
  expect_s3_class(sized, "svyplan_cluster")
  expect_length(sized$n, 2L)

  prec <- prec_multi_cluster(
    ratio_row(n = sized$n[1L], n_per_psu = sized$n[2L], icc_psu = 0.05),
    stage_cost = c(500, 50)
  )
  expect_equal(prec$detail$.cv, 0.05, tolerance = 1e-6)
})

test_that("existing indicator tables are unaffected", {
  prop_only <- data.frame(p = c(0.3, 0.5), moe = 0.05)
  mean_only <- data.frame(var = c(100, 200), mu = c(50, 60), cv = 0.05)
  mixed <- data.frame(p = c(0.3, NA), var = c(NA, 100), mu = c(NA, 50),
                      cv = 0.05)

  expect_s3_class(n_multi(prop_only), "svyplan_n")
  expect_s3_class(n_multi(mean_only), "svyplan_n")
  expect_s3_class(n_multi(mixed), "svyplan_n")
  expect_false("r" %in% names(n_multi(prop_only)$indicators))
  expect_false("cv_num" %in% names(n_multi(mean_only)$indicators))
})
