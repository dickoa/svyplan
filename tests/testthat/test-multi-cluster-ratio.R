## A ratio row reaches the optimizer as a unit relative variance and a CV
## target, exactly as a mean row does, so these benchmark against the direct
## n_cluster()/prec_cluster() call under the same stage inputs.

## L_R = 1.20^2 + 0.45^2 - 2 * 0.65 * 1.20 * 0.45 = 0.9405 exactly.
CL_RELVAR <- 0.9405

## The two routes agree to the optimizer's convergence, ~1.6e-9, not bit for
## bit. Not ratio-specific: a mean row shows 1.61e-9. Asserting tighter would
## test the search's arithmetic rather than the claim.
OPT_TOL <- 1e-7

## Budget mode converges more loosely, ~1.46e-6, also across estimands:
## ratio 1.456819e-06, mean 1.456818e-06, proportion 1.456821e-06.
BUDGET_TOL <- 1e-5

cluster_row <- function(...) {
  base <- data.frame(r = 420, cv_num = 1.20, cv_den = 0.45,
                     component_cor = 0.65)
  extra <- list(...)
  for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  base
}

test_that("the derived coefficient is the one the optimizer receives", {
  expect_equal(.ratio_unit_relvar(420, 1.20, 0.45, 0.65), CL_RELVAR)
})

test_that("a one-row two-stage ratio table matches n_cluster()", {
  tbl <- cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1)
  via_table <- n_cluster(indicators = tbl, stage_cost = c(500, 50))
  direct <- n_cluster(stage_cost = c(500, 50), icc = 0.05,
                      unit_relvar = CL_RELVAR, var_ratio = 1, cv = 0.05)

  expect_equal(unname(via_table$n), unname(direct$n), tolerance = OPT_TOL)
  expect_equal(via_table$total_n, direct$total_n, tolerance = OPT_TOL)
  expect_equal(via_table$cost, direct$cost, tolerance = OPT_TOL)
})

test_that("a one-row three-stage ratio table matches n_cluster()", {
  tbl <- cluster_row(cv = 0.05, icc_psu = 0.04, icc_ssu = 0.10,
                     var_ratio_psu = 1)
  via_table <- n_cluster(indicators = tbl, stage_cost = c(800, 120, 40))
  direct <- n_cluster(stage_cost = c(800, 120, 40),
                      icc = c(0.04, 0.10), unit_relvar = CL_RELVAR,
                      var_ratio = 1, cv = 0.05)

  expect_equal(unname(via_table$n), unname(direct$n), tolerance = OPT_TOL)
  expect_equal(via_table$total_n, direct$total_n, tolerance = OPT_TOL)
})

test_that("a moe target converts exactly on the ratio scale", {
  by_moe <- n_cluster(
    indicators = cluster_row(moe = 20, icc_psu = 0.05, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )
  by_cv <- n_cluster(
    indicators = cluster_row(cv = 20 / (qnorm(0.975) * 420), icc_psu = 0.05,
                var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )

  expect_equal(unname(by_moe$n), unname(by_cv$n), tolerance = OPT_TOL)
})

test_that("budget mode routes a ratio row to the same allocation", {
  tbl <- cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1)
  via_table <- n_cluster(indicators = tbl, stage_cost = c(500, 50), budget = 80000)
  direct <- n_cluster(stage_cost = c(500, 50), icc = 0.05,
                      unit_relvar = CL_RELVAR, var_ratio = 1, budget = 80000)

  expect_equal(unname(via_table$n), unname(direct$n), tolerance = BUDGET_TOL)
})

test_that("a fixed stage size is honoured on a ratio row", {
  tbl <- cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1)
  via_table <- n_cluster(indicators = tbl, stage_cost = c(500, 50), n_per_psu = 12)
  direct <- n_cluster(stage_cost = c(500, 50), icc = 0.05,
                      unit_relvar = CL_RELVAR, var_ratio = 1, cv = 0.05,
                      n_per_psu = 12)

  expect_equal(unname(via_table$n), unname(direct$n), tolerance = OPT_TOL)
  expect_equal(unname(via_table$n)[2L], 12)
})

test_that("a ratio row of indicators matches the scalar prec_cluster()", {
  via_table <- prec_cluster(
    indicators = cluster_row(n = 45, n_per_psu = 14, icc_psu = 0.05, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )
  direct <- prec_cluster(n = c(45, 14), icc = 0.05,
                         unit_relvar = CL_RELVAR, var_ratio = 1)

  expect_equal(via_table$detail$.cv, direct$cv, tolerance = OPT_TOL)
})

test_that("a moe-sized ratio row reports se and moe on the ratio scale", {
  # The reconstruction runs only in moe mode, reached through the round trip
  # from a moe-sized result. In cv mode every estimand reports .cv alone.
  sized <- n_cluster(
    indicators = cluster_row(moe = 20, icc_psu = 0.05, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )
  res <- prec_cluster(sized)

  expect_equal(res$detail$.se, res$detail$.cv * 420, tolerance = 1e-10)
  expect_equal(res$detail$.moe, qnorm(0.975) * res$detail$.se,
               tolerance = 1e-10)
  expect_equal(res$detail$.rmoe, res$detail$.moe / 420, tolerance = 1e-10)
  expect_equal(res$detail$.moe, 20, tolerance = 1e-5)
})

test_that("the multistage round trip returns the same allocation", {
  tbl <- cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1)
  sized <- n_cluster(indicators = tbl, stage_cost = c(500, 50))
  prec <- prec_cluster(
    indicators = cluster_row(n = sized$n[1L], n_per_psu = sized$n[2L], icc_psu = 0.05,
                var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )

  expect_equal(prec$detail$.cv, 0.05, tolerance = 1e-6)
})

test_that("a mixed cluster table sizes all three estimands", {
  tbl <- data.frame(
    name = c("literacy", "income", "consumption"),
    p = c(0.35, NA, NA),
    var = c(NA, 2500, NA),
    mu = c(NA, 300, NA),
    r = c(NA, NA, 420),
    cv_num = c(NA, NA, 1.20),
    cv_den = c(NA, NA, 0.45),
    component_cor = c(NA, NA, 0.65),
    cv = 0.05,
    icc_psu = 0.05
  )
  res <- n_cluster(indicators = tbl, stage_cost = c(500, 50))

  expect_equal(nrow(res$detail), 3L)
  expect_true(all(is.finite(res$detail$.cv_target)))
  expect_equal(sum(res$detail$.binding), 1L)
})

test_that("a component's icc cannot stand in for the linearized one", {
  # VDK Example 9.3: the clustering effect on a ratio and on a component
  # total differ by an order of magnitude in the same population. Planning
  # with the wrong one is a different design, not a conservative one.
  ratio_icc <- n_cluster(
    indicators = cluster_row(cv = 0.05, icc_psu = 0.00088, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )
  component_icc <- n_cluster(
    indicators = cluster_row(cv = 0.05, icc_psu = 0.02251, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )

  expect_false(isTRUE(all.equal(unname(ratio_icc$n),
                                unname(component_icc$n))))
  expect_gt(component_icc$total_n, ratio_icc$total_n)
})

test_that("stage response rates reach a ratio row", {
  with_rate <- n_cluster(
    indicators = cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1,
                resp_rate_psu = 0.8),
    stage_cost = c(500, 50)
  )
  without <- n_cluster(
    indicators = cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )

  expect_gt(with_rate$n[1L], without$n[1L])
})

test_that("the operational design is whole units", {
  res <- n_cluster(
    indicators = cluster_row(cv = 0.05, icc_psu = 0.05, var_ratio_psu = 1),
    stage_cost = c(500, 50)
  )
  op <- res$operational$n

  expect_true(all(op == floor(op)))
  expect_true(all(op > 0))
})

test_that("an incomplete ratio quartet is still refused in cluster mode", {
  expect_error(
    n_cluster(indicators = data.frame(r = 420, cv_num = 1.2, cv = 0.05,
                               icc_psu = 0.05),
                    stage_cost = c(500, 50)),
    "column\\(s\\) .*'cv_den'.*are absent"
  )
})

test_that("a cluster row may still declare only one estimand", {
  expect_error(
    n_cluster(indicators = data.frame(p = 0.3, r = 420, cv_num = 1.2, cv_den = 0.45,
                               component_cor = 0.65, cv = 0.05,
                               icc_psu = 0.05),
                    stage_cost = c(500, 50)),
    "only one of 'p', 'var', or 'r'"
  )
})
