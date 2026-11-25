## Oracle for every budget case below: a hand count. Stratum A takes its
## only unit, which leaves room for exactly one unit of cost 1 in B. A
## relative slack of 1e-9 on a budget of 1e10 would buy nine more.

budget <- 1e10 + 1

test_that("a joint budget fit buys no unit beyond the budget", {
  frame <- data.frame(stratum = c("A", "B"), N = c(1, 10),
                      unit_cost = c(1e10, 1))
  measures <- data.frame(stratum = frame$stratum, name = "y",
                         mean = 1, sd = c(0, 1))
  for (rule in c("sum", "worst")) {
    fit <- n_alloc_joint(frame, measures, objective = "y", budget = budget,
                         fpc = "none", objective_rule = rule)
    expect_equal(fit$detail$n, c(1, 1), tolerance = 1e-9, label = rule)
    expect_identical(as.numeric(fit$detail$n_int), c(1, 1), label = rule)
    expect_lte(fit$params$achieved$cost, budget)
    expect_lte(fit$operational$cost, budget)
  }
})

test_that("an element and a cluster allocation buy no unit beyond the budget", {
  frame <- data.frame(stratum = c("A", "B"), N = c(1, 10), sd = c(0, 1),
                      mean = 1, unit_cost = c(1e10, 1))
  fit <- n_alloc(frame, budget = budget, alloc = "optimal")
  expect_identical(as.numeric(fit$detail$n_int), c(1, 1))

  frame <- data.frame(stratum = c("A", "B"), N_1 = c(1, 10), N = c(5, 50),
                      take_2 = c(1, 1), cost_1 = c(1e10 - 1, 0.5),
                      cost_2 = c(1, 0.5), sd = c(0, 1), mean = 1,
                      icc_1 = 0.1)
  checked <- character(0)
  for (take in c("fixed", "free")) {
    for (fpc in c("unit", "stage")) {
      f <- if (take == "fixed") frame else frame[names(frame) != "take_2"]
      fit <- suppressWarnings(n_alloc(f, budget = budget, alloc = "optimal",
                                      fpc = fpc))
      label <- paste(take, fpc)
      d <- fit$detail
      expect_identical(as.numeric(d$n_1_int), c(1, 1), label = label)
      expect_identical(as.numeric(d$take_2_int), c(1, 1), label = label)
      checked <- c(checked, label)
    }
  }
  expect_setequal(checked, c("fixed unit", "fixed stage", "free unit",
                             "free stage"))
})

test_that("a cluster field design buys no unit beyond the budget", {
  ## One PSU at 1e10 - 4 leaves room for five second-stage units.
  fit <- suppressWarnings(n_cluster(stage_cost = c(1e10 - 4, 1), icc = 0.1,
                                    budget = budget, N = c(1, 30)))
  expect_identical(as.numeric(fit$operational$n), c(1, 5))
  expect_lte(fit$operational$cost, budget)
})

test_that("the free-take worst floor reads a near-census remainder exactly", {
  ## Oracle: with every PSU at take 5 of 20, the stage-corrected variance of
  ## the mean is 3 delta / 200 in closed form.
  frame <- data.frame(stratum = "A", N = 200, N_1 = 10, cost_1 = 5,
                      cost_2 = 1)
  for (k in 46:48) {
    delta <- 2^-k
    measures <- data.frame(stratum = "A", name = "y", mean = 1, sd = 1,
                           icc_1 = 1 - delta)
    fit <- n_alloc_joint(frame, measures, objective = "y", budget = 100,
                         fpc = "stage", objective_rule = "worst",
                         take_2_min = 5, take_2_max = 5)
    label <- sprintf("delta = 2^-%d", k)
    search <- fit$optimization$search
    expected <- 3 * delta / 200
    expect_equal(fit$operational$objective_value / expected, 1,
                 tolerance = 1e-12, label = label)
    expect_identical(fit$optimization$classification, "optimal",
                     label = label)
    expect_identical(search$status, "optimal", label = label)
    expect_equal(search$bound / expected, 1, tolerance = 1e-12,
                 label = label)
    expect_identical(search$trials, 0L, label = label)
  }
})

test_that("a rule allocation rounds within the budget", {
  frame <- data.frame(stratum = c("A", "B"), N = c(1, 10), sd = 1, mean = 1,
                      unit_cost = c(1e10, 1))
  checked <- character(0)
  for (alloc in c("proportional", "neyman", "power")) {
    fit <- suppressWarnings(n_alloc(frame, budget = budget, alloc = alloc))
    expect_identical(as.numeric(fit$detail$n_int), c(1, 1), label = alloc)
    checked <- c(checked, alloc)
  }
  frame <- data.frame(stratum = c("A", "B"), N_1 = c(1, 10), N = c(5, 50),
                      take_2 = 1, cost_1 = c(1e10 - 1, 0.5),
                      cost_2 = c(1, 0.5), sd = 1, mean = 1, icc_1 = 0.1)
  for (alloc in c("proportional", "neyman")) {
    fit <- suppressWarnings(n_alloc(frame, budget = budget, alloc = alloc))
    label <- paste("cluster", alloc)
    expect_identical(as.numeric(fit$detail$n_1_int), c(1, 1), label = label)
    checked <- c(checked, label)
  }
  ## Three stages, one PSU of 1e10 in A and PSUs of cost 1 in B.
  frame <- data.frame(stratum = c("A", "B"), N_1 = c(1, 10),
                      N_2 = c(1, 10), N = c(1, 10), take_2 = 1, take_3 = 1,
                      cost_1 = c(1e10 - 2, 0.5), cost_2 = c(1, 0.25),
                      cost_3 = c(1, 0.25), sd = 1, mean = 1, icc_1 = 0.1,
                      icc_2 = 0.1)
  fit <- suppressWarnings(n_alloc(frame, budget = budget,
                                  alloc = "proportional"))
  expect_identical(as.numeric(fit$detail$n_1_int), c(1, 1))
  expect_lte(fit$operational$cost, budget)
  checked <- c(checked, "three stages")
  expect_setequal(checked, c("proportional", "neyman", "power",
                             "cluster proportional", "cluster neyman",
                             "three stages"))
})

test_that("a domain floor is met within the budget or refused", {
  ## A takes 1e10 and D at least one unit, so B + C = 3 needs 1e10 + 4.
  frame <- data.frame(stratum = c("A", "D", "B", "C"),
                      dom = c("d1", "d1", "d2", "d2"), N = c(1, 10, 10, 10),
                      sd = c(1, 5, 1, 1), mean = 1,
                      unit_cost = c(1e10, 1, 1, 1))
  floor <- data.frame(domain = "dom", level = "d2", min_n = 3)
  run <- function(b) {
    suppressWarnings(n_alloc(frame, budget = b, alloc = "optimal",
                             domains = "dom", min_n_domain = floor))
  }
  expect_error(run(1e10 + 3), "within the budget")
  fit <- run(1e10 + 4)
  n <- as.numeric(fit$detail$n_int)
  expect_identical(n[1:2], c(1, 1))
  expect_identical(sum(n[3:4]), 3)
})

test_that("a joint cluster split buys no unit beyond the budget", {
  ## Two PSUs of 5e9 - 2 leave 6 units, three per symmetric domain.
  tab <- data.frame(name = "y", p = 0.3, cv = 0.1, icc_1 = 0.05,
                    region = c("A", "B"), N = 30, N_1 = 1)
  fit <- suppressWarnings(n_cluster(tab, stage_cost = c(5e9 - 2, 1),
                                    domains = "region", allocation = "joint",
                                    budget = 1e10 + 2, fpc = "stage"))
  expect_identical(as.numeric(fit$operational$domains$.total_n), c(3, 3))
  expect_lte(fit$operational$cost, 1e10 + 2)
})

test_that("a two-phase design buys no unit beyond the budget", {
  frame <- data.frame(stratum = c("A", "B"), N = 500, sd = c(0.01, 1),
                      mean = 10, unit_cost = c(1e10, 1))
  for (b in c(1e10 + 5, 1e10 + 20)) {
    fit <- suppressWarnings(n_twophase(frame, phase1_cost = 0.25,
                                       budget = b))
    expect_lte(fit$operational$cost, b)
  }
})

test_that("a register plan buys no unit beyond the budget", {
  ## A certainty PSU of 1e10 - 4 beside PSUs of cost 1 in B.
  psu <- data.frame(stratum = c("A", rep("B", 12)), N = 5)
  frame <- data.frame(stratum = c("A", "B"), N = c(5, 60), take_2 = c(5, 1),
                      cost_1 = c(1e10 - 5, 0.5), cost_2 = c(1, 0.5))
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 10,
                         sd = c(1, 3), icc_1 = 0.1)
  run <- function(b) {
    suppressWarnings(n_alloc_joint(frame, measures, objective = "y",
                                   budget = b, psu = psu))
  }
  for (b in c(1e10 + 2, 1e10 + 4)) {
    fit <- run(b)
    expect_lte(fit$operational$cost, b)
    expect_lte(fit$params$achieved$cost, b)
  }
  ## The cheapest plan, one unit in A and one PSU in B, costs 1e10 - 3.
  expect_error(run(1e10 - 4), "below the cheapest field design")
  expect_lte(run(1e10 - 3)$operational$cost, 1e10 - 3)
})

test_that("an exact cheapest register plan proves a budget below it short", {
  ## The cheapest plan meeting the target is found exactly, so a budget a
  ## relative 5e-10 below its cost is refused, not left unknown.
  psu <- data.frame(stratum = c("A", rep("B", 12)), N = 5)
  frame <- data.frame(stratum = c("A", "B"), N = c(5, 60), take_2 = c(5, 1),
                      cost_1 = c(95, 0.5), cost_2 = c(1, 0.5))
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 10,
                         sd = c(1, 3), icc_1 = 0.1)
  targets <- data.frame(name = "y", cv = 0.03)
  cheapest <- suppressWarnings(n_alloc_joint(frame, measures,
                                             targets = targets, psu = psu))
  expect_identical(cheapest$optimization$search$status, "optimal")
  cost <- cheapest$operational$cost
  run <- function(b) {
    suppressWarnings(n_alloc_joint(frame, measures, targets = targets,
                                   objective = "y", budget = b, psu = psu))
  }
  expect_error(run(cost * (1 - 5e-10)), "cannot fund the precision targets")
  expect_equal(run(cost)$operational$cost, cost)
})
