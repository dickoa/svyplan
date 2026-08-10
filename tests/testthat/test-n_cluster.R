test_that("n_cluster 2-stage budget mode", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_s3_class(result, "svyplan_cluster")

  n2_opt <- sqrt(500 / 50 * (1 - 0.05) / 0.05)
  n1_opt <- 100000 / (500 + 50 * n2_opt)
  cv_expected <- sqrt(1 / (n1_opt * n2_opt) * 1 * (1 + 0.05 * (n2_opt - 1)))

  expect_equal(result$n[["n_per_psu"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_opt, tolerance = 1e-6)
  expect_equal(result$cv, cv_expected, tolerance = 1e-6)
  expect_equal(result$cost, 100000)
  expect_equal(result$stages, 2L)
})

test_that("n_cluster 2-stage CV mode", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)

  n2_opt <- sqrt(500 / 50 * (1 - 0.05) / 0.05)
  n1_opt <- 1 * 1 * (1 + 0.05 * (n2_opt - 1)) / (n2_opt * 0.05^2)
  cost_expected <- 500 * n1_opt + 50 * n1_opt * n2_opt

  expect_equal(result$n[["n_per_psu"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_opt, tolerance = 1e-6)
  expect_equal(result$cost, cost_expected, tolerance = 1e-4)
})

test_that("n_cluster 2-stage fixed m budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    budget = 100000,
    n_psu = 40
  )

  n2_expected <- (100000 - 500 * 40) / (50 * 40)
  cv_expected <- sqrt(
    1 * 1 / (40 * n2_expected) * (1 + 0.05 * (n2_expected - 1))
  )

  expect_equal(result$n[["n_psu"]], 40, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$cv, cv_expected, tolerance = 1e-6)
})

test_that("n_cluster 2-stage fixed m CV mode", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, n_psu = 40)

  n2_expected <- (1 - 0.05) / (0.05^2 * 40 / (1 * 1) - 0.05)
  cost_expected <- 500 * 40 + 50 * 40 * n2_expected

  expect_equal(result$n[["n_psu"]], 40, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$cost, cost_expected, tolerance = 1e-4)
})

test_that("n_cluster 3-stage budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    budget = 500000
  )

  k2 <- 1 - 0.01
  n3_opt <- sqrt((1 - 0.05) / 0.05 * 100 / 50)
  n2_opt <- 1 / n3_opt * sqrt((1 - 0.05) / 0.01 * 500 / 50 * k2 / 1)
  n1_opt <- 500000 / (500 + 100 * n2_opt + 50 * n2_opt * n3_opt)

  expect_equal(result$n[["n_per_ssu"]], n3_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_opt, tolerance = 1e-6)
  expect_equal(result$stages, 3L)
})

test_that("n_cluster 3-stage CV mode", {
  result <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)

  k2 <- 1 - 0.01
  n3_opt <- sqrt((1 - 0.05) / 0.05 * 100 / 50)
  n2_opt <- 1 / n3_opt * sqrt((1 - 0.05) / 0.01 * 500 / 50 * k2 / 1)
  n1_opt <- 1 /
    (0.05^2 * n2_opt * n3_opt) *
    (1 * 0.01 * n2_opt * n3_opt + k2 * (1 + 0.05 * (n3_opt - 1)))

  expect_equal(result$n[["n_per_ssu"]], n3_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_opt, tolerance = 1e-6)
})

test_that("n_cluster accepts svyplan_varcomp", {
  vc <- .new_svyplan_varcomp(
    varb = 0.01,
    varw = 1.0,
    icc = 0.05,
    var_ratio = 1.0,
    unit_relvar = 1.0,
    stages = 2L
  )
  result <- n_cluster(stage_cost = c(500, 50), icc = vc, budget = 100000)
  expect_s3_class(result, "svyplan_cluster")
  expect_equal(result$stages, 2L)
})

test_that("n_cluster round-trips with prec_cluster", {
  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  cv_check <- prec_cluster(n = unname(plan$n), icc = 0.05)$cv
  expect_equal(cv_check, plan$cv, tolerance = 1e-6)
})

test_that("n_cluster rejects boundary icc values", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0, budget = 1e5),
    "stay away from 0 and 1"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 1, cv = 0.05),
    "stay away from 0 and 1"
  )
})

test_that("n_cluster rejects near-boundary icc values", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 1e-12, budget = 1e5),
    "stay away from 0 and 1"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 1 - 1e-12, cv = 0.05),
    "stay away from 0 and 1"
  )
})

test_that("n_cluster rejects svyplan_varcomp with near-zero icc", {
  vc <- .new_svyplan_varcomp(
    varb = 1e-12,
    varw = 1,
    icc = 1e-12,
    var_ratio = 1,
    unit_relvar = 1,
    stages = 2L
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = vc, cv = 0.1),
    "stay away from 0 and 1"
  )
})

test_that("n_cluster validates inputs", {
  expect_error(
    n_cluster(stage_cost = 500, icc = 0.05, budget = 100000),
    "length >= 2"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05),
    "specify exactly one"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, budget = 100000),
    "specify exactly one"
  )
  expect_error(
    n_cluster(
      stage_cost = c(500, 50, 20, 10),
      icc = c(0.01, 0.02, 0.03),
      budget = 100000
    ),
    "fold deeper stages"
  )
})

test_that("n_cluster is an S3 generic", {
  expect_true(is.function(n_cluster))
  # Verify dispatch works on numeric (default method)
  res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_s3_class(res, "svyplan_cluster")
})

test_that("n_cluster params store cv in CV mode", {
  res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  expect_equal(res$params$cv, 0.05)
  expect_null(res$params$budget)
  expect_null(res$params$n_psu)
})

test_that("n_cluster params store budget and m", {
  res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000, n_psu = 40)
  expect_equal(res$params$budget, 100000)
  expect_equal(res$params$n_psu, 40)
  expect_null(res$params$cv)
})

test_that("n_cluster params store budget without m", {
  res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_equal(res$params$budget, 100000)
  expect_null(res$params$n_psu)
})

test_that("n_cluster 3-stage params store cv", {
  res <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)
  expect_equal(res$params$cv, 0.05)
  expect_null(res$params$budget)
})

test_that("n_cluster 3-stage params store budget and m", {
  res <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    budget = 500000,
    n_psu = 50
  )
  expect_equal(res$params$budget, 500000)
  expect_equal(res$params$n_psu, 50)
})

test_that("svyplan_cluster has se/moe/cv fields", {
  res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_true("se" %in% names(res))
  expect_true("moe" %in% names(res))
  expect_true("cv" %in% names(res))
  expect_true(is.na(res$se))
  expect_true(is.na(res$moe))
  expect_true(res$cv > 0)
})

test_that("n_cluster rejects invalid unit_relvar", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, unit_relvar = -1, cv = 0.05),
    "must be positive"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, unit_relvar = 0, cv = 0.05),
    "must be positive"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, unit_relvar = NA_real_, cv = 0.05),
    "must not be NA"
  )
})

test_that("n_cluster rejects invalid var_ratio", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = -1, cv = 0.05),
    "positive"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = 0, cv = 0.05),
    "positive"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = NA, cv = 0.05),
    "positive"
  )
  expect_error(
    n_cluster(
      stage_cost = c(500, 100, 50),
      icc = c(0.01, 0.05),
      var_ratio = c(1, -1),
      cv = 0.05
    ),
    "positive"
  )
})

test_that("n_cluster fixed-m CV error reports actionable floor diagnostic", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.001, n_psu = 5),
    "below the achievable floor"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.001, n_psu = 5),
    "increase n_psu to at least"
  )
  expect_error(
    n_cluster(
      stage_cost = c(500, 100, 50),
      icc = c(0.01, 0.05),
      cv = 0.001,
      n_psu = 5
    ),
    "below the achievable floor"
  )
})

test_that("fixed_cost = 0 is backward compatible", {
  r1 <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  r2 <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    budget = 100000,
    fixed_cost = 0
  )
  expect_equal(r1$n, r2$n)
  expect_equal(r1$cost, r2$cost)
  expect_null(r2$params$fixed_cost)
})

test_that("fixed_cost 2-stage budget mode reduces n1, n2 unchanged", {
  base <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  fc <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    budget = 100000,
    fixed_cost = 5000
  )
  expect_equal(fc$n[["n_per_psu"]], base$n[["n_per_psu"]], tolerance = 1e-10)
  expect_true(fc$n[["n_psu"]] < base$n[["n_psu"]])
  expect_equal(fc$cost, 100000)
  expect_equal(fc$params$fixed_cost, 5000)
})

test_that("fixed_cost 2-stage CV mode leaves n unchanged, adds to cost", {
  base <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  fc <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, fixed_cost = 5000)
  expect_equal(fc$n, base$n, tolerance = 1e-10)
  expect_equal(fc$cost, base$cost + 5000, tolerance = 1e-6)
})

test_that("fixed_cost 3-stage budget mode reduces n1", {
  base <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    budget = 500000
  )
  fc <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    budget = 500000,
    fixed_cost = 10000
  )
  expect_equal(fc$n[["n_per_psu"]], base$n[["n_per_psu"]], tolerance = 1e-10)
  expect_equal(fc$n[["n_per_ssu"]], base$n[["n_per_ssu"]], tolerance = 1e-10)
  expect_true(fc$n[["n_psu"]] < base$n[["n_psu"]])
  expect_equal(fc$cost, 500000)
})

test_that("fixed_cost 3-stage CV mode adds to cost", {
  base <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)
  fc <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    cv = 0.05,
    fixed_cost = 10000
  )
  expect_equal(fc$n, base$n, tolerance = 1e-10)
  expect_equal(fc$cost, base$cost + 10000, tolerance = 1e-6)
})

test_that("fixed_cost validation rejects bad inputs", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, fixed_cost = -1),
    "non-negative"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, fixed_cost = NA),
    "non-negative"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, fixed_cost = c(1, 2)),
    "non-negative"
  )
  expect_error(
    n_cluster(
      stage_cost = c(500, 50),
      icc = 0.05,
      budget = 100000,
      fixed_cost = 100000
    ),
    "less than"
  )
  expect_error(
    n_cluster(
      stage_cost = c(500, 50),
      icc = 0.05,
      budget = 100000,
      fixed_cost = 200000
    ),
    "less than"
  )
})

test_that("as.integer returns the operational stage vector", {
  res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_identical(as.integer(res), as.integer(res$operational$n))
  expect_named(res$operational$n, c("n_psu", "n_per_psu"))
})

test_that("fixed_cost round-trip n_cluster -> prec_cluster -> n_cluster", {
  orig <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    cv = 0.05,
    fixed_cost = 5000
  )
  prec <- prec_cluster(orig)
  expect_equal(prec$params$fixed_cost, 5000)
  back <- n_cluster(prec)
  expect_equal(unname(back$n), unname(orig$n), tolerance = 1e-4)
  expect_equal(back$params$fixed_cost, 5000)
})

test_that("cluster display totals come from the operational design", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  op <- x$operational
  expect_equal(op$total_n, prod(op$n))
  expect_match(format(x), as.character(op$total_n))
  out <- capture.output(print(x))
  expect_true(any(grepl(sprintf("total n = %d", op$total_n), out)))
})

test_that("n_cluster accepts named icc vector in correct order", {
  ref <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)
  named <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(icc_psu = 0.01, icc_ssu = 0.05),
    cv = 0.05
  )
  expect_equal(named$n, ref$n)
  expect_equal(named$cv, ref$cv)
})

test_that("n_cluster reorders named icc vector", {
  ref <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)
  swapped <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(icc_ssu = 0.05, icc_psu = 0.01),
    cv = 0.05
  )
  expect_equal(swapped$n, ref$n)
  expect_equal(swapped$cv, ref$cv)
})

test_that("n_cluster rejects bad names in icc", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = c(foo = 0.05), cv = 0.05),
    "unrecognized names"
  )
})

test_that("n_cluster accepts named var_ratio vector", {
  ref <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    var_ratio = c(1.2, 1.2 * 0.99),
    cv = 0.05
  )
  swapped <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    var_ratio = c(var_ratio_ssu = 1.2 * 0.99, var_ratio_psu = 1.2),
    cv = 0.05
  )
  expect_equal(swapped$n, ref$n)
})

test_that("n_cluster accepts named cost vectors and reorders", {
  ref2 <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  named2 <- n_cluster(
    stage_cost = c(cost_ssu = 50, cost_psu = 500),
    icc = 0.05,
    cv = 0.05
  )
  alias2 <- n_cluster(
    stage_cost = c(cost_tsu = 50, cost_psu = 500),
    icc = 0.05,
    cv = 0.05
  )
  expect_equal(named2$n, ref2$n)
  expect_equal(alias2$n, ref2$n)

  ref3 <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)
  named3 <- n_cluster(
    stage_cost = c(cost_tsu = 50, cost_psu = 500, cost_ssu = 100),
    icc = c(0.01, 0.05),
    cv = 0.05
  )
  expect_equal(named3$n, ref3$n)
})

test_that("n_cluster rejects bad names in cost", {
  expect_error(
    n_cluster(stage_cost = c(foo = 500, bar = 50), icc = 0.05, cv = 0.05),
    "unrecognized names"
  )
  expect_error(
    n_cluster(stage_cost = c(cost1 = 500, cost2 = 50), icc = 0.05, cv = 0.05),
    "unrecognized names"
  )
})

test_that("n_cluster stores stage-named cost in params (2-stage)", {
  x <- n_cluster(stage_cost = c(750, 100), icc = 0.05, budget = 1e5)
  expect_equal(names(x$params$stage_cost), c("cost_psu", "cost_ssu"))
  expect_equal(unname(x$params$stage_cost), c(750, 100))
})

test_that("n_cluster stores stage-named cost/icc/var_ratio in params (3-stage)", {
  x <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    var_ratio = c(1.2, 1.2 * 0.99), cv = 0.05
  )
  expect_equal(names(x$params$stage_cost), c("cost_psu", "cost_ssu", "cost_tsu"))
  expect_equal(names(x$params$icc), c("icc_psu", "icc_ssu"))
  expect_equal(names(x$params$var_ratio), c("var_ratio_psu", "var_ratio_ssu"))
})

test_that("cost names survive prec_cluster round-trip", {
  x <- n_cluster(stage_cost = c(750, 100), icc = 0.05, budget = 1e5)
  pr <- prec_cluster(x)
  x2 <- n_cluster(pr, budget = 1e5)
  expect_equal(names(x2$params$stage_cost), c("cost_psu", "cost_ssu"))
  expect_equal(x2$n, x$n, tolerance = 1e-6)
})

test_that("n_cluster 2-stage fixed n_per_psu budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    budget = 100000,
    n_per_psu = 20
  )

  n1_expected <- 100000 / (500 + 50 * 20)
  n1_eff <- n1_expected * 1
  cv_expected <- sqrt(1 * 1 / (n1_eff * 20) * (1 + 0.05 * (20 - 1)))

  expect_equal(result$n[["n_per_psu"]], 20, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cv, cv_expected, tolerance = 1e-6)
  expect_equal(result$cost, 100000)
})

test_that("n_cluster 2-stage fixed n_per_psu CV mode", {
  result <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    cv = 0.05,
    n_per_psu = 20
  )

  n1_expected <- 1 * 1 * (1 + 0.05 * (20 - 1)) / (20 * 0.05^2)
  cost_expected <- 500 * n1_expected + 50 * n1_expected * 20

  expect_equal(result$n[["n_per_psu"]], 20, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cv, 0.05, tolerance = 1e-6)
  expect_equal(result$cost, cost_expected, tolerance = 1e-4)
})

# Conditional on a fixed middle take, the final-stage optimum is not the
# unrestricted sqrt((1 - icc_ssu) / icc_ssu * C2 / C3): that take ignores both
# the fixed n_per_psu and the per-PSU cost C1. These two tests previously
# asserted the unrestricted formula, which is the defect rather than the
# contract.
n3_conditional <- function(n2, icc_psu, icc_ssu, stage_cost) {
  k2 <- 1 - icc_psu
  a <- icc_psu + k2 * icc_ssu / n2
  g <- k2 * (1 - icc_ssu) / n2
  sqrt(g * (stage_cost[1L] + stage_cost[2L] * n2) /
         (a * stage_cost[3L] * n2))
}

test_that("n_cluster 3-stage fixed n_per_psu budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    budget = 500000,
    n_per_psu = 10
  )

  n3_opt <- n3_conditional(10, 0.01, 0.05, c(500, 100, 50))
  n1_expected <- 500000 / (500 + 100 * 10 + 50 * 10 * n3_opt)

  expect_equal(result$n[["n_per_psu"]], 10, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], n3_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cost, 500000)

  # The property the formula exists for: no other final take buys more
  # precision from the same budget.
  cv_at <- function(q) {
    n1 <- 500000 / (500 + 100 * 10 + 50 * 10 * q)
    prec_cluster(n = c(n1, 10, q), icc = c(0.01, 0.05))$cv
  }
  expect_equal(result$n[["n_per_ssu"]],
               stats::optimize(cv_at, c(1.0001, 200), tol = 1e-10)$minimum,
               tolerance = 1e-5)
})

test_that("n_cluster 3-stage fixed n_per_psu CV mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    cv = 0.05,
    n_per_psu = 10
  )

  k2 <- 1 - 0.01
  n3_opt <- n3_conditional(10, 0.01, 0.05, c(500, 100, 50))
  n1_expected <- 1 /
    (0.05^2 * 10 * n3_opt) *
    (1 * 0.01 * 10 * n3_opt + k2 * (1 + 0.05 * (n3_opt - 1)))

  expect_equal(result$n[["n_per_psu"]], 10, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], n3_opt, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cv, 0.05, tolerance = 1e-6)

  # And no other final take reaches that CV more cheaply.
  cost_at <- function(q) {
    p <- prec_cluster(n = c(1, 10, q), icc = c(0.01, 0.05))
    (p$cv / 0.05)^2 * (500 + 100 * 10 + 50 * 10 * q)
  }
  expect_equal(result$n[["n_per_ssu"]],
               stats::optimize(cost_at, c(1.0001, 200), tol = 1e-10)$minimum,
               tolerance = 1e-5)
  expect_lt(result$cost,
            local({
              q <- sqrt((1 - 0.05) / 0.05 * 100 / 50)
              cost_at(q)
            }))
})

test_that("n_cluster rejects all stages fixed (2-stage)", {
  expect_error(
    n_cluster(
      stage_cost = c(500, 50), icc = 0.05, budget = 100000,
      n_psu = 40, n_per_psu = 20
    ),
    "cannot fix all stages"
  )
})

test_that("n_cluster rejects n_per_ssu for 2-stage", {
  expect_error(
    n_cluster(
      stage_cost = c(500, 50), icc = 0.05, budget = 100000,
      n_per_ssu = 8
    ),
    "not applicable"
  )
})

test_that("n_cluster rejects all stages fixed (3-stage)", {
  expect_error(
    n_cluster(
      stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), budget = 500000,
      n_psu = 50, n_per_psu = 10, n_per_ssu = 8
    ),
    "cannot fix all stages"
  )
})

test_that("n_cluster params store n_per_psu", {
  res <- n_cluster(
    stage_cost = c(500, 50), icc = 0.05, budget = 100000, n_per_psu = 20
  )
  expect_equal(res$params$n_per_psu, 20)
  expect_null(res$params$n_psu)
})

test_that("n_per_psu round-trip n_cluster -> prec_cluster -> n_cluster", {
  orig <- n_cluster(
    stage_cost = c(500, 50),
    icc = 0.05,
    cv = 0.05,
    n_per_psu = 20
  )
  prec <- prec_cluster(orig)
  expect_equal(prec$params$n_per_psu, 20)
  back <- n_cluster(prec)
  expect_equal(unname(back$n), unname(orig$n), tolerance = 1e-4)
  expect_equal(back$params$n_per_psu, 20)
})

test_that("fixed_cost with n_per_psu 2-stage budget mode reduces n1", {
  base <- n_cluster(
    stage_cost = c(500, 50), icc = 0.05, budget = 100000, n_per_psu = 20
  )
  fc <- n_cluster(
    stage_cost = c(500, 50), icc = 0.05, budget = 100000,
    n_per_psu = 20, fixed_cost = 5000
  )
  expect_equal(fc$n[["n_per_psu"]], 20, tolerance = 1e-10)
  expect_true(fc$n[["n_psu"]] < base$n[["n_psu"]])
  expect_equal(fc$cost, 100000)
})

test_that("n_cluster 3-stage n_per_ssu only budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_per_ssu = 8
  )
  k2 <- 1 - 0.01
  n3 <- 8
  n2_expected <- sqrt(
    k2 * (1 + 0.05 * (n3 - 1)) * 500 /
      (n3 * 1 * 0.01 * (100 + 50 * n3))
  )
  n1_expected <- 500000 / (500 + 100 * n2_expected + 50 * n2_expected * n3)

  expect_equal(result$n[["n_per_ssu"]], 8, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cost, 500000)
})

test_that("n_cluster 3-stage n_per_ssu only CV mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    cv = 0.05, n_per_ssu = 8
  )
  k2 <- 1 - 0.01
  n3 <- 8
  n2 <- sqrt(
    k2 * (1 + 0.05 * (n3 - 1)) * 500 /
      (n3 * 1 * 0.01 * (100 + 50 * n3))
  )
  n1_expected <- 1 / (0.05^2 * n2 * n3) *
    (1 * 0.01 * n2 * n3 + k2 * (1 + 0.05 * (n3 - 1)))

  expect_equal(result$n[["n_per_ssu"]], 8, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cv, 0.05, tolerance = 1e-6)
})

test_that("n_cluster 3-stage n_psu + n_per_ssu budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_psu = 50, n_per_ssu = 8
  )
  n3 <- 8
  n2_expected <- (500000 / 50 - 500) / (100 + 50 * n3)

  expect_equal(result$n[["n_psu"]], 50, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], 8, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$cost, 500000)
})

test_that("n_cluster 3-stage n_psu + n_per_ssu CV mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    cv = 0.05, n_psu = 50, n_per_ssu = 8
  )
  k2 <- 1 - 0.01
  n3 <- 8
  n1_eff <- 50
  n2_expected <- k2 * (1 + 0.05 * (n3 - 1)) /
    (n3 * (0.05^2 * n1_eff / 1 - 1 * 0.01))

  expect_equal(result$n[["n_psu"]], 50, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], 8, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], n2_expected, tolerance = 1e-6)
  expect_equal(result$cv, 0.05, tolerance = 1e-6)
})

test_that("n_cluster 3-stage n_per_psu + n_per_ssu budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_per_psu = 10, n_per_ssu = 8
  )
  n1_expected <- 500000 / (500 + 100 * 10 + 50 * 10 * 8)

  expect_equal(result$n[["n_per_psu"]], 10, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], 8, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_expected, tolerance = 1e-6)
  expect_equal(result$cost, 500000)
})

test_that("n_cluster 3-stage n_per_psu + n_per_ssu CV mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    cv = 0.05, n_per_psu = 10, n_per_ssu = 8
  )
  k2 <- 1 - 0.01
  n1_eff_expected <- 1 / (0.05^2 * 10 * 8) *
    (1 * 0.01 * 10 * 8 + k2 * (1 + 0.05 * (8 - 1)))

  expect_equal(result$n[["n_per_psu"]], 10, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], 8, tolerance = 1e-6)
  expect_equal(result$n[["n_psu"]], n1_eff_expected, tolerance = 1e-6)
  expect_equal(result$cv, 0.05, tolerance = 1e-6)
})

test_that("n_cluster 3-stage n_psu + n_per_psu budget mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_psu = 50, n_per_psu = 10
  )
  n3_expected <- (500000 - 500 * 50 - 100 * 50 * 10) / (50 * 50 * 10)

  expect_equal(result$n[["n_psu"]], 50, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], 10, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], n3_expected, tolerance = 1e-6)
  expect_equal(result$cost, 500000)
})

test_that("n_cluster 3-stage n_psu + n_per_psu CV mode", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    cv = 0.0465, n_psu = 50, n_per_psu = 10
  )
  k2 <- 1 - 0.01
  n1_eff <- 50
  denom <- 10 * (0.0465^2 * n1_eff / 1 - 1 * 0.01) - k2 * 0.05
  n3_expected <- k2 * (1 - 0.05) / denom

  expect_equal(result$n[["n_psu"]], 50, tolerance = 1e-6)
  expect_equal(result$n[["n_per_psu"]], 10, tolerance = 1e-6)
  expect_equal(result$n[["n_per_ssu"]], n3_expected, tolerance = 1e-6)
  expect_gte(result$n[["n_per_ssu"]], 1)
  expect_equal(result$cv, 0.0465, tolerance = 1e-6)
})

test_that("n_cluster clamps sub-unit stage sizes to 1 and reports achieved cv", {
  result <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    cv = 0.05, n_psu = 50, n_per_psu = 10
  )
  expect_equal(result$n[["n_per_ssu"]], 1)
  cv_expected <- sqrt(1 / (50 * 10) * (0.01 * 10 + (1 - 0.01)))
  expect_equal(result$cv, cv_expected, tolerance = 1e-8)
  expect_lt(result$cv, 0.05)
})

test_that("n_cluster params store n_per_ssu", {
  res <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_per_ssu = 8
  )
  expect_equal(res$params$n_per_ssu, 8)
  expect_null(res$params$n_psu)
  expect_null(res$params$n_per_psu)
})

test_that("n_per_ssu round-trip n_cluster -> prec_cluster -> n_cluster", {
  orig <- n_cluster(
    stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    cv = 0.05, n_psu = 50, n_per_ssu = 8
  )
  prec <- prec_cluster(orig)
  expect_equal(prec$params$n_per_ssu, 8)
  expect_equal(prec$params$n_psu, 50)
  back <- n_cluster(prec)
  expect_equal(unname(back$n), unname(orig$n), tolerance = 1e-4)
})

test_that("n_cluster 3-stage all combos verify against prec_cluster", {
  verify <- function(...) {
    res <- n_cluster(...)
    cv_check <- prec_cluster(n = unname(res$n), icc = c(0.01, 0.05))$cv
    expect_equal(cv_check, res$cv, tolerance = 1e-6)
  }
  verify(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_per_ssu = 8)
  verify(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_psu = 50, n_per_ssu = 8)
  verify(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_per_psu = 10, n_per_ssu = 8)
  verify(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
    budget = 500000, n_psu = 50, n_per_psu = 10)
})

test_that("prec_cluster accepts named icc and reorders", {
  ref <- prec_cluster(n = c(50, 12, 8), icc = c(0.01, 0.05))
  swapped <- prec_cluster(
    n = c(50, 12, 8),
    icc = c(icc_ssu = 0.05, icc_psu = 0.01)
  )
  expect_equal(swapped$cv, ref$cv)
})

test_that("budget mode rejects budgets below one realizable PSU", {
  expect_error(
    suppressWarnings(n_cluster(stage_cost = c(1, 1000), icc = 0.9, budget = 100)),
    "too small for any realizable design"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), budget = 600),
    "too small for any realizable design"
  )
})

test_that("cost-optimal n_per_psu below 1 is clamped with a warning", {
  expect_warning(
    res <- n_cluster(stage_cost = c(1, 1000), icc = 0.9, budget = 5000),
    "clamped to 1"
  )
  expect_equal(res$n[["n_per_psu"]], 1)
  expect_gte(res$n[["n_psu"]], 1)
})

test_that("fixed stage sizes below 1 are rejected", {
  expect_error(n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                         n_per_psu = 0.5),
               "'n_per_psu' must be at least 1")
})

test_that("operational budget designs never exceed the budget", {
  grid <- expand.grid(C1 = c(50, 500, 2000), C2 = c(10, 50),
                      icc = c(0.02, 0.3), budget = c(1500, 25000))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    x <- suppressWarnings(try(
      n_cluster(stage_cost = c(g$C1, g$C2), icc = g$icc, budget = g$budget),
      silent = TRUE
    ))
    if (inherits(x, "try-error")) next
    op <- x$operational
    expect_lte(op$cost, g$budget + 1e-8)
    expect_true(all(op$n >= 1))
    expect_equal(op$cv, prec_cluster(n = op$n, icc = g$icc)$cv,
                 tolerance = 1e-10)
  }
})

test_that("operational cv designs meet the target at whole units", {
  for (target in c(0.03, 0.05, 0.1)) {
    x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = target)
    op <- x$operational
    expect_lte(op$cv, target + 1e-10)
    expect_true(all(op$n == round(op$n)))
    x3 <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
                    cv = target)
    op3 <- x3$operational
    expect_lte(op3$cv, target + 1e-10)
    expect_gte(x3$cost, 0.99 * op3$cost - 2 * (500 + 100 + 50))
  }
})

test_that("3-stage operational budget design fits and beats naive ceiling", {
  x <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
                 budget = 20000)
  op <- x$operational
  expect_lte(op$cost, 20000)
  expect_equal(op$total_n, prod(op$n))
  naive <- ceiling(x$n)
  naive_cost <- naive[[1]] * (500 + 100 * naive[[2]] + 50 * naive[[2]] * naive[[3]])
  expect_true(naive_cost > 20000 || op$cv <= x$cv * 1.5)
})

## Three-stage whole-unit search over the ultimate take

# Reference design, written straight from the documented three-stage cv so it
# shares no algebra with the search it checks. Enumerates every whole design in
# the box and returns the cheapest that meets the target, or the most precise
# that fits the budget.
op3_reference <- function(stage_cost, icc, var_ratio, unit_relvar = 1,
                          cv = NULL, budget = NULL, n_psu = NULL,
                          n_per_psu = NULL, n_per_ssu = NULL,
                          resp_rate_psu = 1, resp_rate_ssu = 1, resp_rate = 1,
                          m_max = 60L, q_max = 20000L) {
  ms <- if (is.null(n_per_psu)) seq_len(m_max) else as.integer(n_per_psu)
  qs <- if (is.null(n_per_ssu)) seq_len(q_max) else as.integer(n_per_ssu)
  m <- rep(ms, times = length(qs))
  q <- rep(qs, each = length(ms))
  mr <- m * resp_rate_ssu
  qr <- q * resp_rate
  numer <- unit_relvar / (resp_rate_psu * mr * qr) *
    (var_ratio[1L] * icc[1L] * mr * qr +
       var_ratio[2L] * (1 + icc[2L] * (qr - 1)))
  per_psu <- stage_cost[1L] + stage_cost[2L] * m + stage_cost[3L] * m * q
  if (is.null(budget)) {
    a <- if (is.null(n_psu)) pmax(1, ceiling(numer / cv^2 - 1e-9)) else
      rep(as.integer(n_psu), length(m))
    key <- a * per_psu
    key[numer / a > cv^2 + 1e-12] <- Inf
  } else {
    a <- if (is.null(n_psu)) floor(budget / per_psu) else
      rep(as.integer(n_psu), length(m))
    key <- numer / a
    key[a < 1 | a * per_psu > budget + 1e-8] <- Inf
  }
  i <- which.min(key)
  list(n = c(a[i], m[i], q[i]), cost = a[i] * per_psu[i],
       cv = sqrt(numer[i] / a[i]))
}

test_that("the three-stage search finds the whole-unit optimum a q window misses", {
  x <- suppressWarnings(n_cluster(
    cv = 0.05, stage_cost = c(1000, 200, 0.02),
    icc = c(0.02, 0.0005), unit_relvar = 1.5
  ))
  op <- x$operational
  expect_identical(as.integer(op$n), c(13L, 1L, 833L))
  expect_equal(op$cost, 15816.58, tolerance = 1e-6)
  expect_lte(op$cv, 0.05)
  # The optimum sits where the PSU count is a tight ceiling, five times below
  # the continuous take and past any window anchored on it.
  expect_gt(x$n[["n_per_ssu"]], 4000)
  expect_gt(op$n[["n_per_ssu"]], 400)
})

test_that("three-stage cv designs match an exhaustive whole-unit search", {
  grid <- expand.grid(
    C3 = c(0.02, 5, 60), icc2 = c(0.0005, 0.05), cv = c(0.05, 0.1),
    V = c(1, 1.5)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    icc <- c(0.02, g$icc2)
    vr <- c(1, 1 - 0.02)
    x <- suppressWarnings(n_cluster(
      cv = g$cv, stage_cost = c(1000, 200, g$C3), icc = icc,
      unit_relvar = g$V
    ))
    ref <- op3_reference(c(1000, 200, g$C3), icc, vr, g$V, cv = g$cv)
    expect_lte(x$operational$cv, g$cv + 1e-12)
    expect_equal(x$operational$cost, ref$cost, tolerance = 1e-9)
  }
})

test_that("three-stage budget designs match an exhaustive whole-unit search", {
  grid <- expand.grid(C3 = c(1, 20, 50), icc2 = c(0.005, 0.05),
                      budget = c(20000, 250000))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    icc <- c(0.01, g$icc2)
    vr <- c(1, 1 - 0.01)
    x <- suppressWarnings(n_cluster(
      budget = g$budget, stage_cost = c(500, 100, g$C3), icc = icc
    ))
    ref <- op3_reference(c(500, 100, g$C3), icc, vr, budget = g$budget)
    expect_lte(x$operational$cost, g$budget + 1e-8)
    expect_lte(x$operational$cv, ref$cv * (1 + 1e-9))
  }
})

test_that("a cheap final stage buys a take past any window on the continuous one", {
  # The most precise affordable design spends the residual budget on ultimate
  # units, so the take runs well past where a window around the continuous
  # optimum would stop.
  x <- suppressWarnings(n_cluster(budget = 5e4, stage_cost = c(1000, 200, 0.02),
                                  icc = c(0.02, 0.0005)))
  op <- x$operational
  expect_gt(op$n[["n_per_ssu"]], 400)
  expect_lte(op$cost, 5e4 + 1e-8)
  ref <- op3_reference(c(1000, 200, 0.02), c(0.02, 0.0005), c(1, 1 - 0.02),
                       budget = 5e4, q_max = 6000L)
  expect_lte(op$cv, ref$cv * (1 + 1e-9))
})

test_that("the three-stage search is exact with stages fixed", {
  cost <- c(500, 100, 20)
  icc <- c(0.01, 0.03)
  vr <- c(1, 1 - 0.01)
  combos <- list(
    list(n_psu = 40), list(n_per_psu = 6), list(n_per_ssu = 4),
    list(n_psu = 40, n_per_psu = 6), list(n_psu = 40, n_per_ssu = 4),
    list(n_per_psu = 6, n_per_ssu = 4)
  )
  for (fix in combos) {
    for (mode in list(list(cv = 0.05), list(budget = 3e5))) {
      args <- c(list(stage_cost = cost, icc = icc), fix, mode)
      x <- suppressWarnings(do.call(n_cluster, args))
      ref <- do.call(op3_reference, c(
        list(cost, icc, vr), fix, mode, list(q_max = 5000L)
      ))
      if (is.null(mode$budget)) {
        expect_lte(x$operational$cv, 0.05 + 1e-12)
        expect_equal(x$operational$cost, ref$cost, tolerance = 1e-9)
      } else {
        expect_lte(x$operational$cost, 3e5 + 1e-8)
        expect_lte(x$operational$cv, ref$cv * (1 + 1e-9))
      }
    }
  }
})

test_that("stage response rates keep the three-stage search exact", {
  cost <- c(800, 150, 3)
  icc <- c(0.02, 0.01)
  vr <- c(1, 1 - 0.02)
  x <- suppressWarnings(n_cluster(
    cv = 0.06, stage_cost = cost, icc = icc, resp_rate_psu = 0.9,
    resp_rate_ssu = 0.85, resp_rate = 0.8
  ))
  ref <- op3_reference(cost, icc, vr, cv = 0.06, resp_rate_psu = 0.9,
                       resp_rate_ssu = 0.85, resp_rate = 0.8)
  expect_lte(x$operational$cv, 0.06 + 1e-12)
  expect_equal(x$operational$cost, ref$cost, tolerance = 1e-9)
})

test_that("a loose target settles on a single PSU rather than below one", {
  x <- suppressWarnings(n_cluster(cv = 0.4, stage_cost = c(500, 100, 20),
                                  icc = c(0.005, 0.02)))
  op <- x$operational
  expect_identical(op$n[["n_psu"]], 1L)
  expect_lte(op$cv, 0.4 + 1e-12)
  ref <- op3_reference(c(500, 100, 20), c(0.005, 0.02), c(1, 1 - 0.005),
                       cv = 0.4)
  expect_equal(op$cost, ref$cost, tolerance = 1e-9)
})

test_that("n_multi_cluster reaches the same three-stage design as n_cluster", {
  ind <- data.frame(name = "x", var = 1.5, mu = 1, cv = 0.05,
                    icc_psu = 0.02, icc_ssu = 0.0005)
  nc <- suppressWarnings(n_cluster(cv = 0.05, icc = c(0.02, 0.0005),
                                   unit_relvar = 1.5,
                                   stage_cost = c(1000, 200, 0.02)))
  nm <- suppressWarnings(n_multi_cluster(ind, stage_cost = c(1000, 200, 0.02)))
  expect_identical(as.integer(nm$operational$n), as.integer(nc$operational$n))
  expect_equal(nm$operational$cost, nc$operational$cost)

  ncb <- suppressWarnings(n_cluster(budget = 2e5, icc = c(0.02, 0.005),
                                    stage_cost = c(500, 100, 20)))
  nmb <- suppressWarnings(n_multi_cluster(
    data.frame(name = "x", var = 1, mu = 1, cv = 0.05,
               icc_psu = 0.02, icc_ssu = 0.005),
    stage_cost = c(500, 100, 20), budget = 2e5
  ))
  expect_identical(as.integer(nmb$operational$n), as.integer(ncb$operational$n))
})

test_that("the binding indicator drives the multi-indicator three-stage search", {
  ind <- data.frame(
    name = c("loose", "tight"), var = c(1, 1.5), mu = c(1, 1),
    cv = c(0.2, 0.05), icc_psu = c(0.01, 0.02), icc_ssu = c(0.02, 0.0005)
  )
  nm <- suppressWarnings(n_multi_cluster(ind, stage_cost = c(1000, 200, 0.02)))
  tight <- suppressWarnings(n_multi_cluster(ind[2, ],
                                            stage_cost = c(1000, 200, 0.02)))
  expect_identical(as.integer(nm$operational$n),
                   as.integer(tight$operational$n))
  expect_true(all(nm$operational$cv_by_target <= ind$cv + 1e-12))
})
