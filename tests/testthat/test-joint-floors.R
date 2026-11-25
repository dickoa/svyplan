## Regional sample floors in joint mode

## Whole designs of a small element frame, every one listed.
floor_states <- function(frame) {
  grid <- expand.grid(lapply(frame$N, seq_len))
  names(grid) <- frame$stratum
  as.matrix(grid)
}

floor_relvars <- function(states, frame, measures, names) {
  total <- sum(frame$N)
  vapply(names, function(nm) {
    sd <- measures$sd[measures$name == nm]
    mu <- sum(frame$N * measures$mean[measures$name == nm]) / total
    drop((1 / states) %*% (frame$N * sd)^2) / (total * mu)^2
  }, numeric(nrow(states)))
}

test_that("a regional floor reproduces the education survey's allocation", {
  ## Oracle: Willems and Muennich (2026), Section 6, the budget 241 fit with
  ## a regional cv of 0.04, stratum floors of 1% of N and a floor of 61
  ## units in each region.
  z <- bethel_reference()
  N <- z$N
  frame <- z$frame
  regional <- z$measures[z$measures$name == "Floorspace", ]
  regional$name <- "Regional"
  regional$mean <- sum(N) / stats::ave(N, frame$region, FUN = sum)
  floor <- data.frame(domain = "region", level = c("r1", "r2", "r3"),
                      min_n = 61)
  targets <- data.frame(name = "Regional", domain = "region", cv = 0.04)
  fit <- n_alloc_joint(frame, rbind(z$measures, regional), targets = targets,
                       objective = rownames(z$d_qh), budget = 241,
                       fpc = "none", min_n_stratum = 0.01 * N,
                       min_n_domain = floor)
  expect_equal(fit$detail$n, c(85.36, 30.58, 25.31, 38.76, 37.08, 23.92),
               tolerance = 3e-4)
  expect_equal(fit$domain_floors$achieved[3], 61, tolerance = 1e-8)
  expect_true(all(fit$domain_floors$.pass))
  expect_true(all(fit$operational$domain_floors$.pass))
  expect_true(all(fit$operational$constraints$.pass))
  expect_lte(fit$operational$cost, 241)
  assessed <- prec_alloc_joint(fit, n = fit$detail$n_int)
  expect_equal(assessed$domain_floors$achieved,
               fit$operational$domain_floors$achieved)
})

test_that("a fractional floor buys the best whole design that meets it", {
  ## Oracle: every whole design of three strata, under both rules.
  frame <- data.frame(stratum = letters[1:3], N = c(9, 8, 9),
                      region = c("N", "N", "S"), unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 3), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1))
  floors <- data.frame(domain = "region", level = "N", min_n = 8.2)
  states <- floor_states(frame)
  keep <- rowSums(states[, 1:2]) >= 9 & drop(states %*% frame$unit_cost) <= 21
  relvar <- floor_relvars(states[keep, ], frame, measures, c("y", "z"))
  checked <- character(0)
  for (rule in c("sum", "worst")) {
    fit <- n_alloc_joint(frame, measures, objective = c("y", "z"),
                         budget = 21, objective_rule = rule,
                         min_n_domain = floors, fpc = "none")
    best <- if (rule == "sum") min(rowSums(relvar)) else
      min(apply(relvar, 1L, max))
    expect_equal(fit$operational$objective_value, best, tolerance = 1e-12,
                 label = rule)
    expect_identical(fit$operational$domain_floors$required, 9)
    expect_identical(fit$domain_floors$required, 8.2)
    expect_gte(fit$domain_floors$achieved, 8.2 * (1 - 1e-12))
    again <- n_alloc_joint(prec_alloc_joint(fit, n = fit$detail$n_int))
    expect_identical(again$params$min_n_domain$min_n, 8.2)
    expect_false(summary(prec_alloc_joint(fit, n = c(2, 2, 2)))$overall$all_pass)
    checked <- c(checked, rule)
  }
  expect_setequal(checked, c("sum", "worst"))
})

test_that("a whole design exactly on its floor survives budget pruning", {
  ## Oracle: the enumerated best design at budget 37, which sits on the
  ## floor of ceiling(9.2) = 10 units in region N.
  frame <- data.frame(stratum = letters[1:3], N = 9:7,
                      region = c("N", "N", "S"), unit_cost = c(6, 2, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 3), mean = 1,
                         sd = c(1.90730886, 0.63679961, 1.90721897,
                                1.95284597, 1.97822312, 1.39572087))
  floors <- data.frame(domain = "region", level = "N", min_n = 9.2)
  fit <- n_alloc_joint(frame, measures, objective = c("y", "z"), budget = 37,
                       min_n_domain = floors, fpc = "none")
  states <- floor_states(frame)
  keep <- rowSums(states[, 1:2]) >= 10 & drop(states %*% frame$unit_cost) <= 37
  relvar <- rowSums(floor_relvars(states[keep, ], frame, measures,
                                  c("y", "z")))
  expect_identical(fit$operational$status, "optimal")
  expect_identical(as.numeric(fit$detail$n_int),
                   as.numeric(states[keep, ][which.min(relvar), ]))
  expect_identical(fit$operational$domain_floors$achieved, 10)
})

test_that("whole floored designs match enumeration under both rules", {
  ## Oracle: every whole design of small frames with one or two floored
  ## regions, in cost mode and budget mode.
  set.seed(20261125)
  checked <- character(0)
  for (i in 1:12) {
    H <- 4
    frame <- data.frame(stratum = letters[1:H], N = sample(3:6, H, TRUE),
                        region = rep(c("N", "S"), each = 2),
                        unit_cost = sample(1:4, H, TRUE))
    measures <- data.frame(stratum = rep(frame$stratum, 2),
                           name = rep(c("y", "z"), each = H), mean = 1,
                           sd = round(stats::runif(2 * H, 0.5, 3), 2))
    levels <- if (i %% 2 == 0) c("N", "S") else "N"
    cap <- tapply(frame$N, frame$region, sum)[levels]
    floors <- data.frame(domain = "region", level = levels,
                         min_n = pmax(2, round(0.6 * cap) - 0.5))
    states <- floor_states(frame)
    held <- rep(TRUE, nrow(states))
    for (k in seq_len(nrow(floors))) {
      cols <- frame$region == floors$level[k]
      held <- held & rowSums(states[, cols, drop = FALSE]) >=
        ceiling(floors$min_n[k])
    }
    cost <- drop(states %*% frame$unit_cost)
    budget <- floor(stats::quantile(cost[held], 0.3))
    keep <- held & cost <= budget
    relvar <- floor_relvars(states[keep, , drop = FALSE], frame, measures,
                            c("y", "z"))
    for (rule in c("sum", "worst")) {
      fit <- n_alloc_joint(frame, measures, objective = c("y", "z"),
                           budget = budget, objective_rule = rule,
                           min_n_domain = floors, fpc = "none")
      best <- if (rule == "sum") min(rowSums(relvar)) else
        min(apply(relvar, 1L, max))
      label <- sprintf("frame %d %s", i, rule)
      expect_identical(fit$operational$status, "optimal", label = label)
      expect_equal(fit$operational$objective_value, best, tolerance = 1e-12,
                   label = label)
      expect_true(all(fit$operational$domain_floors$.pass), label = label)
      checked <- c(checked, label)
    }
  }
  expect_length(checked, 24L)
})

test_that("a take-all stratum counts toward its region's floor", {
  ## Oracle: A is a census of 4 units, so the floor of 6 in region N
  ## leaves 2 for B. The continuous cost minimum is then closed form.
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(4, 20, 20),
                      region = c("N", "N", "S"), unit_cost = 1,
                      take_all = c(TRUE, FALSE, FALSE))
  measures <- data.frame(stratum = frame$stratum, name = "y", mean = 1,
                         sd = c(1, 1e-3, 1))
  floors <- data.frame(domain = "region", level = "N", min_n = 6)
  fit <- n_alloc_joint(frame, measures, objective = "y", budget = 12,
                       min_n_domain = floors, fpc = "none")
  expect_equal(fit$detail$n[1:2], c(4, 2), tolerance = 1e-8)
  expect_equal(fit$detail$n[3], 6, tolerance = 1e-8)
  expect_identical(as.numeric(fit$detail$n_int), c(4, 2, 6))
})

test_that("a floored budget is refused below the cheapest floored design", {
  ## A takes its only unit at 1e10, and the floor of 3 in region S costs 3.
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(1, 10, 10),
                      region = c("N", "S", "S"), unit_cost = c(1e10, 1, 1))
  measures <- data.frame(stratum = frame$stratum, name = "y", mean = 1,
                         sd = 1)
  floors <- data.frame(domain = "region", level = "S", min_n = 3)
  run <- function(b) {
    n_alloc_joint(frame, measures, objective = "y", budget = b,
                  min_n_domain = floors, fpc = "none")
  }
  ## The cheapest design is accepted up to a relative 1e-6 over the
  ## budget, so the refusal starts below that.
  expect_error(run(1e10 * (1 - 2e-6)), "regional sample floors")
  expect_equal(run(1e10 + 2)$operational$cost, 1e10 + 3)
  fit <- run(1e10 + 3)
  n <- as.numeric(fit$detail$n_int)
  expect_identical(n[1], 1)
  expect_identical(sum(n[2:3]), 3)
  expect_lte(fit$operational$cost, 1e10 + 3)
})

test_that("joint floors refuse what they cannot hold", {
  frame <- data.frame(stratum = letters[1:3], N = c(9, 8, 9),
                      region = c("N", "N", "S"), unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = frame$stratum, name = "y", mean = 1,
                         sd = 1)
  floor <- function(min_n, level = "N") {
    data.frame(domain = "region", level = level, min_n = min_n)
  }
  run <- function(...) {
    n_alloc_joint(frame, measures, objective = "y", budget = 40,
                  fpc = "none", ...)
  }
  expect_error(run(min_n_domain = floor(18)), "capacity")
  expect_error(run(min_n_domain = floor(3, "W")), "empty")
  expect_error(run(min_n_domain = data.frame(level = "N", min_n = 3)),
               "domain")
  expect_error(run(min_n_domain = floor(-1)), "min_n")
  expect_error(run(min_n_domain = 3), "table")
})

test_that("a binding floor's slope matches a refit on both sides", {
  ## Oracle: central differences of the minimum cost, and of the objective
  ## at a fixed budget, in the floor.
  frame <- data.frame(stratum = letters[1:4], N = c(400, 300, 500, 200),
                      region = c("N", "N", "S", "S"),
                      unit_cost = c(2, 1, 3, 1))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 4), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1, 2, 1))
  floor <- function(f) data.frame(domain = "region", level = "N", min_n = f)
  cost_fit <- function(f) {
    n_alloc_joint(frame, measures,
                  targets = data.frame(name = c("y", "z"), cv = 0.1),
                  min_n_domain = floor(f), fpc = "none")
  }
  f0 <- 1.3 * cost_fit(0)$domain_floors$achieved
  fit <- cost_fit(f0)
  h <- 1e-3 * f0
  expect_identical(fit$optimization$slopes, "available")
  expect_equal(fit$domain_floors$.sensitivity,
               (cost_fit(f0 + h)$params$achieved$cost -
                  cost_fit(f0 - h)$params$achieved$cost) / (2 * h),
               tolerance = 1e-5)
  budget_fit <- function(f) {
    n_alloc_joint(frame, measures, objective = c("y", "z"), budget = 900,
                  min_n_domain = floor(f), fpc = "none")
  }
  f1 <- 1.3 * budget_fit(0)$domain_floors$achieved
  fit <- budget_fit(f1)
  h <- 1e-3 * f1
  expect_identical(fit$optimization$slopes, "available")
  expect_equal(fit$domain_floors$.objective_sensitivity,
               (budget_fit(f1 + h)$objective_value -
                  budget_fit(f1 - h)$objective_value) / (2 * h),
               tolerance = 1e-5)
  expect_identical(budget_fit(0)$domain_floors$.objective_sensitivity, 0)
})

test_that("joint floors are refused off the element route", {
  floor <- data.frame(domain = "region", level = "N", min_n = 4)
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 10,
                         sd = 3, icc_1 = 0.1)
  fixed <- data.frame(stratum = c("A", "B"), region = c("N", "S"),
                      N = c(400, 400), N_1 = 20, take_2 = 5, cost_1 = 10,
                      cost_2 = 1)
  free <- fixed[names(fixed) != "take_2"]
  psu <- data.frame(stratum = rep(c("A", "B"), each = 20), N = 20)
  run <- function(frame, ...) {
    n_alloc_joint(frame, measures, objective = "y", budget = 300,
                  min_n_domain = floor, ...)
  }
  expect_error(run(fixed), "element frames")
  expect_error(run(free), "element frames")
  expect_error(run(fixed, psu = psu), "element frames")
})

## predict() grids against direct calls

test_that("a budget grid reads the field objective of each direct fit", {
  z <- bethel_reference()
  objective <- data.frame(name = rownames(z$d_qh), priority = c(1, 2, 1, 3))
  checked <- character(0)
  for (rule in c("sum", "worst")) {
    fit <- n_alloc_joint(z$frame, z$measures, objective = objective,
                         budget = 241, objective_rule = rule, fpc = "none")
    grid <- predict(fit, data.frame(budget = c(230, 260)))
    for (i in 1:2) {
      direct <- n_alloc_joint(z$frame, z$measures, objective = objective,
                              budget = grid$budget[i], objective_rule = rule,
                              fpc = "none")
      label <- sprintf("%s %g", rule, grid$budget[i])
      expect_equal(grid$objective_value[i], direct$objective_value,
                   tolerance = 1e-10, label = label)
      expect_equal(grid$objective_value_int[i],
                   direct$operational$objective_value, tolerance = 1e-12,
                   label = label)
      expect_equal(grid$cost_int[i], direct$operational$cost, label = label)
      checked <- c(checked, label)
    }
  }
  expect_length(checked, 4L)
})

test_that("a floored grid keeps the floors of each direct fit", {
  frame <- data.frame(stratum = letters[1:3], N = c(9, 8, 9),
                      region = c("N", "N", "S"), unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 3), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1))
  floors <- data.frame(domain = "region", level = "N", min_n = 8.2)
  fit <- n_alloc_joint(frame, measures, objective = c("y", "z"), budget = 21,
                       min_n_domain = floors, fpc = "none")
  grid <- predict(fit, data.frame(budget = c(21, 30)))
  for (i in 1:2) {
    direct <- n_alloc_joint(frame, measures, objective = c("y", "z"),
                            budget = grid$budget[i], min_n_domain = floors,
                            fpc = "none")
    expect_equal(grid$n_int[i], direct$operational$n)
    expect_equal(grid$objective_value_int[i],
                 direct$operational$objective_value, tolerance = 1e-12)
  }
  assessed <- prec_alloc_joint(fit, n = c(2, 2, 2))
  expect_false(all(assessed$domain_floors$.pass))
  grid <- predict(assessed, data.frame(deff = c(1, 1.1)))
  expect_identical(grid$.feasible, c(FALSE, FALSE))
  met <- predict(prec_alloc_joint(fit, n = fit$detail$n_int),
                 data.frame(deff = 1))
  expect_identical(met$.feasible, TRUE)
})

## The whole search with falling loads

test_that("the kept sizes of a falling load match a listing of every size", {
  ## Oracle: every size of every branch, kept when its reduced cost is
  ## within the slack of its stratum's least.
  set.seed(11)
  checked <- 0L
  for (i in 1:200) {
    H <- 3
    lo <- rep(1, H)
    hi <- sample(5:40, H, TRUE)
    q <- svyplan:::.whole_problem(
      o1 = runif(H, 0, 2), o2 = runif(H, 0, 30),
      g1 = cbind(runif(H, 0, 0.05), -runif(H, 0, 0.1)),
      g2 = cbind(runif(H, 0, 1), 0), lo = lo, hi = hi,
      g0 = cbind(0, runif(H, 0, 2))
    )
    lambda <- runif(2, 0, 20)
    slack <- runif(1, 0, 5)
    iv <- svyplan:::.whole_intervals(q, lambda, slack)
    abc <- svyplan:::.whole_abc(q, lambda)
    for (b in seq_len(q$B)) {
      x <- q$lo[b]:q$hi[b]
      r <- abc$gamma[b] + abc$alpha[b] * x + abc$beta[b] / x
      rmin <- min(vapply(which(q$owner == q$owner[b]), function(k) {
        z <- q$lo[k]:q$hi[k]
        min(abc$gamma[k] + abc$alpha[k] * z + abc$beta[k] / z)
      }, 1))
      kept <- x[r <= rmin + slack + 1e-12 * abs(rmin + slack)]
      listed <- if (iv$lo[b] <= iv$hi[b]) iv$lo[b]:iv$hi[b] else integer(0)
      expect_true(all(kept %in% listed), label = sprintf("case %d", i))
      checked <- checked + (abc$alpha[b] < 0)
    }
  }
  expect_gt(checked, 50L)
})

test_that("a floor leaves the whole search as strong as without it", {
  ## Oracles: a slack floor leaves the optimum unchanged, and a binding one
  ## is certified within the gap the search reports.
  set.seed(3)
  H <- 10
  frame <- data.frame(stratum = sprintf("s%02d", 1:H),
                      N = round(stats::runif(H, 2e4, 8e4)),
                      region = rep(c("N", "S"), each = 5),
                      unit_cost = round(stats::runif(H, 1, 5), 1))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = H), mean = 10,
                         sd = stats::runif(2 * H, 1, 10))
  run <- function(floor, rule) {
    n_alloc_joint(frame, measures, objective = c("y", "z"), budget = 8000,
                  objective_rule = rule, fpc = "unit",
                  min_n_domain = if (floor > 0) {
                    data.frame(domain = "region", level = "S",
                               min_n = floor + 0.5)
                  })
  }
  free <- run(0, "sum")
  slack <- run(1500, "sum")
  expect_gt(sum(free$detail$n_int[6:10]), 1500)
  expect_identical(slack$operational$status, "optimal")
  expect_equal(slack$operational$objective_value,
               free$operational$objective_value, tolerance = 1e-12)
  tight <- run(2500, "sum")
  expect_identical(tight$operational$status, "optimal")
  expect_identical(sum(tight$detail$n_int[6:10]), 2501)
  worst <- run(2500, "worst")
  expect_lte(worst$optimization$search$gap, 0.01)
  expect_true(all(worst$operational$domain_floors$.pass))
  ## The cost search starts from the continuous design rounded up and its
  ## multipliers: without them the gap reads 0.11 and 33468 states.
  cost <- n_alloc_joint(frame, measures,
                        targets = data.frame(name = c("y", "z"), cv = 0.02),
                        fpc = "unit",
                        min_n_domain = data.frame(domain = "region",
                                                  level = "S", min_n = 900.5))
  expect_lte(cost$optimization$search$gap, 0.01)
  expect_lt(cost$optimization$search$states, 5000)
  expect_true(all(cost$operational$domain_floors$.pass))
})

test_that("a target and a floor on one stratum share a price", {
  ## The floor sits exactly at the size the target alone needs, so both
  ## bind and neither has a slope of its own.
  frame <- data.frame(stratum = c("A", "B"), N = 500, region = c("N", "S"),
                      unit_cost = c(1, 2))
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 1,
                         sd = c(2, 1))
  targets <- data.frame(name = "y", domain = "region", level = "N", cv = 0.1)
  alone <- n_alloc_joint(frame, measures, targets = targets, fpc = "none")
  floor <- data.frame(domain = "region", level = "N",
                      min_n = alone$detail$n[1])
  fit <- n_alloc_joint(frame, measures, targets = targets, fpc = "none",
                       min_n_domain = floor)
  expect_true(fit$constraints$.binding)
  expect_true(is.na(fit$constraints$.sensitivity))
  expect_true(is.na(fit$domain_floors$.sensitivity))
})

test_that("a budget at the cheapest floored cost withholds its slopes", {
  frame <- data.frame(stratum = letters[1:3], N = c(9, 8, 9),
                      region = c("N", "N", "S"), unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 3), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1))
  ## The floor of 8 in N is met most cheaply by B at 7 with A at 1, and C
  ## stays at 1: 2 + 7 + 3 = 12, in continuous and whole units alike.
  fit <- n_alloc_joint(frame, measures, objective = c("y", "z"),
                       budget = 12, fpc = "none",
                       min_n_domain = data.frame(domain = "region",
                                                 level = "N", min_n = 8))
  expect_identical(as.numeric(fit$detail$n_int), c(1, 7, 1))
  expect_identical(fit$optimization$slopes, "budget_at_minimum")
})

test_that("a floor met at a cost the multiplier reaches still certifies", {
  ## At the top of the budget bracket the objective's multiplier is tiny, so
  ## the strata of the floored level carry almost no weight and the floor
  ## stops at the cheapest one's cost. Its size must fill the floor, not
  ## jump to its upper bound.
  frame <- data.frame(stratum = letters[1:4], N = c(400, 300, 500, 200),
                      region = c("N", "N", "S", "S"),
                      unit_cost = c(2, 1, 3, 1))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 4), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1, 2, 1))
  fit <- n_alloc_joint(frame, measures,
                       targets = data.frame(name = "y", domain = "region",
                                            level = "N", cv = 0.12),
                       objective = c("y", "z"), budget = 900, fpc = "none",
                       min_n_domain = data.frame(domain = "region",
                                                 level = "S", min_n = 150))
  expect_identical(fit$optimization$classification, "optimal")
  expect_true(all(fit$operational$domain_floors$.pass))
  expect_lte(fit$objective_value,
             fit$operational$objective_value * (1 + 1e-12))
  cheap <- svyplan:::.bethel_floor_min(c(3, 1), c(1e-30, 1e-30), c(1, 1),
                                       c(500, 200),
                                       list(spec = data.frame(min_n = 150),
                                            membership = matrix(1, 2, 1),
                                            indices = list(1:2)))
  expect_equal(sum(cheap$x), 150)
  expect_equal(cheap$x, c(1, 149))
})
