## Contracts at the edges of the feasible region
##
## Optima inside the feasible region are checked across the suite. These
## pin what happens at its edges, where the 2026-10-27 reviews found every
## numerical defect: the cheapest affordable design, a census, a narrow
## feasible budget split, and a rounded take that moves the stage fractions.

## B1. A budget equal to the cheapest design returns that design

test_that("a budget for one PSU at the smallest take returns one PSU", {
  cases <- list(
    list(stage_cost = c(500, 50), icc = 0.05, resp_rate = 1),
    list(stage_cost = c(20, 2), icc = 0.3, resp_rate = 0.5),
    list(stage_cost = c(300, 50, 10), icc = c(0.04, 0.1), resp_rate = 0.8)
  )
  pops <- list(list(N = 10000, N_psu = 100), list(N = 40, N_psu = 5),
               list(N = 12000, N_psu = 40, N_ssu = 400))
  for (i in seq_along(cases)) {
    cs <- cases[[i]]
    r <- cs$resp_rate
    smallest <- max(1, ceiling(1 / r - 1e-9))
    stages <- length(cs$stage_cost)
    cheapest <- if (stages == 2L) {
      cs$stage_cost[1] + cs$stage_cost[2] * smallest
    } else {
      sum(cs$stage_cost * c(1, 1, smallest))
    }
    for (fpc in c("none", "stage")) {
      args <- c(cs, list(unit_relvar = 2, budget = cheapest, fpc = fpc))
      if (fpc == "stage") args <- c(args, pops[[i]])
      fit <- suppressWarnings(do.call(n_cluster, args))
      expect_gte(fit$n[["n_psu"]], 1 - 1e-9)
      expect_lte(fit$cost, cheapest + 1e-8)
      expect_lte(fit$operational$cost, cheapest + 1e-8)
      expect_equal(fit$operational$n[["n_psu"]], 1)
    }
  }
})

test_that("the review's one-PSU designs under the stage correction", {
  fit <- suppressWarnings(n_cluster(stage_cost = c(500, 50), icc = 0.05,
                                    budget = 600, fpc = "stage", N = 10000,
                                    N_psu = 100))
  expect_equal(unname(fit$n), c(1, 2))
  expect_equal(fit$cv, 0.717635, tolerance = 1e-6)
  low <- suppressWarnings(n_cluster(stage_cost = c(20, 2), icc = 0.3,
                                    unit_relvar = 4, budget = 30,
                                    resp_rate = 0.5, fpc = "stage", N = 40,
                                    N_psu = 5))
  expect_equal(unname(low$operational$n), c(1, 5))
  expect_equal(low$cv, 1.315295, tolerance = 1e-6)
  expect_error(
    n_cluster(stage_cost = c(20, 2), icc = 0.3, unit_relvar = 4, budget = 23,
              resp_rate = 0.5, fpc = "stage", N = 40, N_psu = 5),
    "one PSU with 2 units, costs 24"
  )
})

## B2. A census converts back through precision

test_that("every census survives the round trip through precision", {
  two <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 1e6,
                   fpc = "stage", N = 10000, N_psu = 100)
  three <- n_cluster(stage_cost = c(300, 50, 10), icc = c(0.04, 0.1),
                     budget = 1e7, fpc = "stage", N = 1200, N_psu = 10,
                     N_ssu = 40)
  rows <- n_cluster(
    indicators = data.frame(name = c("y", "z"), p = c(0.3, 0.2),
                            cv = c(0.1, 0.12), icc_psu = c(0.05, 0.1)),
    stage_cost = c(500, 50), budget = 1e6, fpc = "stage", N = 10000,
    N_psu = 100
  )
  for (fit in list(two, three, rows)) {
    expect_equal(fit$cv, 0)
    expect_equal(n_cluster(prec_cluster(fit))$n, fit$n)
  }
  census <- prec_cluster(n = c(100, 100), icc = 0.05, fpc = "stage",
                         N = 10000, N_psu = 100, stage_cost = c(500, 50))
  expect_error(n_cluster(census), "cv is 0 and states no target")
  expect_equal(n_cluster(census, budget = 1e6)$n, two$n)
})

## B3. A feasible joint budget split is found however narrow

test_that("a joint split inside a narrow feasible interval is found", {
  base <- data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
                     region = c("a", "b"))
  for (stage in c(FALSE, TRUE)) {
    tab <- if (stage) transform(base, N = 10000, N_psu = 100) else base
    args <- list(indicators = tab, stage_cost = c(100, 10),
                 domains = "region", n_psu = 2, budget = 450,
                 allocation = "joint")
    if (stage) args$fpc <- "stage"
    fit <- do.call(n_cluster, args)
    expect_equal(fit$domains$.cost, c(225, 225), tolerance = 1e-6)
    expect_equal(fit$domains$n_per_psu, c(1.25, 1.25), tolerance = 1e-6)
    args$budget <- 430
    expect_error(do.call(n_cluster, args),
                 "every domain together costs 440")
  }
})

## B4. A rounded take that moves both stage fractions

test_that("the classic whole-cluster repair meets the stage-corrected target", {
  target <- prec_cluster(n = c(5.5, 2.1), icc = 0.1, unit_relvar = 1,
                         fpc = "stage", N = 30, N_psu = 10)$cv
  frame <- data.frame(stratum = "S", N = 30, sd = 1, mean = 1, icc_psu = 0.1,
                      n_per_psu = 2.1, N_psu = 10, cost_psu = 1, cost_ssu = 10)
  fit <- n_alloc(frame, cv = target, fpc = "stage", alloc = "proportional")
  expect_equal(fit$detail$n_psu_int, 7L)
  expect_equal(fit$detail$n_per_psu_int, 2)
  expect_lte(fit$operational$cv, target)
  expect_equal(
    fit$operational$cv,
    prec_cluster(n = c(7, 2), icc = 0.1, unit_relvar = 1, fpc = "stage",
                 N = 30, N_psu = 10)$cv,
    tolerance = 1e-10
  )
})
