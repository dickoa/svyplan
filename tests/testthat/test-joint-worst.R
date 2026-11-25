# The worst rule of n_alloc_joint(): the least largest weighted relative
# variance within a budget

worst_fit <- function(objective, budget = 241, ...) {
  ref <- bethel_reference()
  n_alloc_joint(ref$frame, ref$measures, objective = objective,
                budget = budget, fpc = "none", objective_rule = "worst", ...)
}

## Refusals

test_that("the worst rule is refused where it is not defined", {
  ref <- bethel_reference()
  expect_error(
    n_alloc_joint(ref$frame, ref$measures,
                  targets = data.frame(name = "Age", cv = 0.06),
                  fpc = "none", objective_rule = "worst"),
    "needs 'objective'"
  )
  expect_error(
    worst_fit(data.frame(name = c("Age", "Oil"), priority = c(1, 0))),
    "every objective 'priority' must be positive"
  )
  expect_error(worst_fit("Age", control = list(gap = 0.01)),
               "'control\\$gap' is refused")
  expect_error(
    n_alloc_joint(ref$frame, ref$measures, objective = "Age", budget = 241,
                  objective_rule = "max"),
    "should be one of"
  )
  register <- data.frame(stratum = c("A", "B"), N = c(2000, 1600),
                         cost_1 = 200, cost_2 = 5, take_2 = c(10, 8))
  psu <- data.frame(stratum = rep(register$stratum, each = 8),
                    N = c(600, rep(200, 7), 450, 220, 165, rep(153, 5)))
  y <- data.frame(stratum = register$stratum, name = "y",
                  p = c(0.5, 0.4), icc_1 = 0.05)
  expect_error(
    n_alloc_joint(register, y, objective = "y", budget = 3000, psu = psu,
                  objective_rule = "worst"),
    "not available with a PSU register"
  )
  expect_error(
    prec_alloc_joint(ref$frame, n = rep(40, 6), measures = ref$measures,
                     targets = data.frame(name = "Age", cv = 0.06),
                     fpc = "none", objective_rule = "worst"),
    "needs 'objective'"
  )
})

## Continuous design

test_that("the continuous worst rule reproduces the published allocations", {
  ## Willems and Muennich (2026, Section 6): n* and the CVs to the digits
  ## printed, lambda* and mu* to their four, and the worst CV against an
  ## independent conic solve of the same problem.
  ref <- bethel_reference()
  components <- rownames(ref$d_qh)
  fit <- worst_fit(components)
  expect_identical(fit$optimization$classification, "optimal")
  expect_equal(round(fit$detail$n, 2),
               c(89.98, 28.83, 26.79, 42.90, 34.42, 18.08))
  expect_equal(round(fit$objective$.cv, 4), c(0.0600, 0.0600, 0.0482, 0.0503))
  expect_equal(sqrt(fit$objective_value), 0.06001741, tolerance = 1e-7)
  expect_equal(fit$params$achieved$cost, 241, tolerance = 1e-9)
  expect_equal(unname(fit$optimization$objective_weights),
               c(0.6669, 0.3331, 0, 0), tolerance = 1e-3)
  expect_equal(-1e5 * fit$optimization$budget_sensitivity, 1.4947,
               tolerance = 1e-3)
  expect_identical(fit$objective$.binding, c(TRUE, TRUE, FALSE, FALSE))

  low <- worst_fit(components, budget = 217)
  expect_equal(sqrt(low$objective_value), 0.06324933, tolerance = 1e-7)

  ## The trade-off weights of the paper, w* scaled by (0.95^-2, 1.0933^-2,
  ## 1, 1).
  trade <- worst_fit(data.frame(
    name = components,
    priority = 1 / (c(0.95^2, 1.0933^2, 1, 1) * fit$objective$.relvar)
  ))
  expect_equal(round(trade$detail$n, 2),
               c(75.21, 34.58, 29.54, 45.81, 36.60, 19.27))
  expect_equal(round(trade$objective$.cv, 4),
               c(0.0587, 0.0642, 0.0483, 0.0517))
  expect_equal(trade$objective_value / trade$objective$priority[1],
               0.003442, tolerance = 1e-3)
})

test_that("the dual weights are sum-rule priorities giving the same design", {
  ref <- bethel_reference()
  components <- rownames(ref$d_qh)
  for (priority in list(c(1, 1, 1, 1), c(2, 1, 3, 1))) {
    worst <- worst_fit(data.frame(name = components, priority = priority))
    weights <- worst$optimization$objective_weights
    expect_equal(sum(weights), 1, tolerance = 1e-12)
    sum_rule <- n_alloc_joint(
      ref$frame, ref$measures,
      objective = data.frame(name = components, priority = unname(weights)),
      budget = 241, fpc = "none"
    )
    expect_equal(sum_rule$detail$n, worst$detail$n, tolerance = 1e-6)
  }
})

test_that("the second stage spends what the worst level leaves", {
  ## The worst of 2 / x1 and 1 / x2 with each count in [1, 2] at unit cost.
  ## Level 1 needs only (2, 1), so a budget of 3.5 has 0.5 left, which the
  ## second stage gives to the second stratum.
  frame <- data.frame(stratum = c("A", "B"), N = 2, unit_cost = 1)
  measures <- data.frame(stratum = c("A", "B", "A", "B"),
                         name = c("y1", "y1", "y2", "y2"), mean = 1,
                         sd = c(sqrt(2), 0, 0, 1))
  objective <- data.frame(name = c("y1", "y2"), domain = "stratum",
                          level = c("A", "B"))
  run <- function(budget) {
    n_alloc_joint(frame, measures, objective = objective, budget = budget,
                  fpc = "none", objective_rule = "worst")
  }
  fit <- run(3.5)
  expect_equal(fit$detail$n, c(2, 1.5), tolerance = 1e-8)
  expect_equal(fit$objective_value, 1, tolerance = 1e-8)
  expect_identical(fit$optimization$second_stage, "optimal")
  expect_equal(run(4)$detail$n, c(2, 2))

  ## Whole counts in [1, 4] at budget 7. Level 0.5 is reached at (4, 2),
  ## cost 6, and at (4, 3), cost 7. The second stage takes (4, 3).
  frame$N <- 4
  whole <- run(7)
  expect_identical(whole$detail$n_int, c(4, 3))
  expect_equal(whole$operational$objective_value, 0.5, tolerance = 1e-12)
  expect_identical(whole$optimization$search$second_stage, "optimal")
})

test_that("the worst-rule slope matches centered refits", {
  ## Oracle: the worst level refitted on both sides of a binding target.
  ref <- bethel_reference()
  objective <- c("Floorspace", "Employees", "Oil")
  run <- function(cv) {
    n_alloc_joint(ref$frame, ref$measures,
                  targets = data.frame(name = "Age", cv = cv),
                  objective = objective, budget = 241, fpc = "none",
                  objective_rule = "worst")
  }
  fit <- run(0.061)
  expect_true(fit$constraints$.binding)
  expect_identical(fit$optimization$slopes, "available")
  h <- 1e-6
  expect_equal((run(0.061 + h)$objective_value -
                  run(0.061 - h)$objective_value) / (2 * h),
               fit$constraints$.objective_sensitivity, tolerance = 1e-5)
  expect_match(
    gsub("\\s+", " ", paste(capture.output(print(summary(fit))),
                            collapse = " ")),
    "change in the worst weighted relative variance", fixed = TRUE
  )
})

## Whole design

test_that("the whole worst design is the reference design with a closed bracket", {
  fit <- worst_fit(rownames(bethel_reference()$d_qh))
  search <- fit$optimization$search
  expect_identical(fit$operational$status, "optimal")
  expect_identical(fit$detail$n_int, c(90, 29, 26, 44, 34, 18))
  expect_equal(sqrt(fit$operational$objective_value), 0.06003494,
               tolerance = 1e-7)
  expect_lte(search$gap, 1e-9)
  expect_lte(search$bound, fit$operational$objective_value)
  expect_lte(fit$operational$cost, 241)
})

test_that("the whole worst design is the best of every affordable design", {
  ## Oracle: every whole design within the budget, scored with the
  ## stratified variance sum_h W_h^2 S_h^2 (1 / n_h - 1 / N_h).
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(12, 8, 5),
                      unit_cost = c(1, 2.5, 4), region = c("n", "n", "s"))
  means <- list(y = c(10, 12, 9), z = c(5, 4, 6))
  sds <- list(y = c(4, 6, 5), z = c(2, 3, 2))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 3),
                         mean = unlist(means), sd = unlist(sds))
  objective <- data.frame(name = c("y", "z", "y"),
                          domain = c(".overall", ".overall", "region"),
                          level = c(NA, NA, "n"), priority = c(1, 2, 0.5))
  relvar <- function(n, h, s, m) {
    W <- frame$N[h] / sum(frame$N[h])
    sum(W^2 * s[h]^2 * (1 / n[h] - 1 / frame$N[h])) / sum(W * m[h])^2
  }
  contributions <- function(n) {
    c(relvar(n, 1:3, sds$y, means$y), 2 * relvar(n, 1:3, sds$z, means$z),
      0.5 * relvar(n, 1:2, sds$y, means$y))
  }
  worst_level <- function(n) max(contributions(n))
  grid <- expand.grid(a = 1:12, b = 1:8, c = 1:5)
  grid$cost <- grid$a + 2.5 * grid$b + 4 * grid$c
  for (budget in c(18, 27.5, 40)) {
    within <- grid[grid$cost <= budget, ]
    best <- min(apply(within[1:3], 1, worst_level))
    fit <- n_alloc_joint(frame, measures, objective = objective,
                         budget = budget, fpc = "unit",
                         objective_rule = "worst")
    label <- sprintf("budget %g", budget)
    expect_equal(fit$operational$objective_value, best, tolerance = 1e-10,
                 label = label)
    expect_equal(worst_level(fit$detail$n_int), best, tolerance = 1e-10,
                 label = label)
    expect_lte(fit$optimization$search$bound, best * (1 + 1e-12))
    expect_identical(fit$operational$status, "optimal")
    ## Among the designs at the least level, the second stage keeps the one
    ## with the least weighted sum. At budget 18 two designs tie.
    tied <- within[apply(within[1:3], 1, worst_level) <= best * (1 + 1e-12),
                   1:3]
    expect_equal(sum(contributions(fit$detail$n_int)),
                 min(apply(tied, 1, function(n) sum(contributions(n)))),
                 tolerance = 1e-10, label = label)
  }

  capped <- n_alloc_joint(frame, measures, objective = objective,
                          budget = 27.5, fpc = "unit",
                          objective_rule = "worst",
                          control = list(max_nodes = 3))
  search <- capped$optimization$search
  expect_identical(capped$operational$status, "limit")
  expect_lte(search$nodes, 3)
  expect_lte(search$bound, capped$operational$objective_value)
  expect_lte(search$bound, 0.0144431220 * (1 + 1e-9))
  expect_true(capped$operational$all_pass)
  expect_lte(capped$operational$cost, 27.5)
})

test_that("the whole worst search keeps the nominal budget of the sum rule", {
  ## Both whole searches buy better designs only within the budget and its
  ## rounding room, so a budget of 2e6 unit-cost units buys 2000000, not the
  ## 2000002 a tolerance of 1e-6 would allow. The bisection's bound stays
  ## at or below the design's level.
  frame <- data.frame(stratum = "A", N = 4e6, unit_cost = 1)
  measures <- data.frame(stratum = "A", name = "y", mean = 1, sd = 1)
  for (rule in c("sum", "worst")) {
    fit <- n_alloc_joint(frame, measures, objective = "y", budget = 2e6,
                         fpc = "none", objective_rule = rule)
    expect_identical(fit$detail$n_int, 2e6, label = rule)
  }
  expect_lte(fit$optimization$search$bound,
             fit$operational$objective_value)
})

test_that("a larger census frame keeps its continuous and whole designs", {
  ## Under fpc = "unit" the floor is a census, and the first trial above it
  ## stops short of certification. Its dual bound still proves the cost
  ## above the budget. The whole bisection starts from the rounded
  ## continuous design, which keeps the gap small within the default nodes.
  set.seed(11)
  H <- 30
  frame <- data.frame(
    stratum = sprintf("s%02d", 1:H),
    N = round(exp(runif(H, log(500), log(20000)))),
    unit_cost = round(runif(H, 1, 4), 2),
    region = rep(sprintf("r%d", 1:5), each = 6)
  )
  names <- c("a", "b", "c", "d")
  measures <- data.frame(stratum = rep(frame$stratum, 4),
                         name = rep(names, each = H),
                         mean = runif(H * 4, 5, 50), sd = runif(H * 4, 2, 30))
  objective <- rbind(
    data.frame(name = names, domain = ".overall", level = NA),
    data.frame(name = "a", domain = "region", level = sprintf("r%d", 1:5))
  )
  for (budget in c(3000, 8000)) {
    fit <- n_alloc_joint(frame, measures, objective = objective,
                         budget = budget, fpc = "unit",
                         objective_rule = "worst")
    search <- fit$optimization$search
    label <- sprintf("budget %g", budget)
    expect_identical(fit$optimization$classification, "optimal", label = label)
    expect_equal(fit$params$achieved$cost, budget, tolerance = 1e-8,
                 label = label)
    expect_lte(fit$operational$cost, budget * (1 + 1e-9))
    expect_lte(search$bound, fit$operational$objective_value)
    expect_lte(search$gap, 0.01)
  }
})

test_that("the whole bisection alone finds the design and spends the leftover", {
  ## Without the rounded start, the bisection and the second stage must find
  ## the designs themselves.
  whole_alone <- function(frame, measures, objective, budget) {
    problem <- svyplan:::.build_bethel_problem(
      frame = frame, measures = measures, targets = NULL, unit_cost = NULL,
      deff = 1, resp_rate = 1, eligible_per_unit = 1, eligibility_rate = 1,
      alpha = 0.05, df = NULL, min_n_stratum = NULL, objective = objective,
      fpc = "none", .finite_prop_var = TRUE
    )
    problem$objective$rule <- "worst"
    svyplan:::.bethel_whole_worst(problem, budget)
  }
  ## A budget of 2e6 buys 2000000 within the nominal budget, and the bound
  ## stays below the design's level.
  one <- whole_alone(data.frame(stratum = "A", N = 4e6, unit_cost = 1),
                     data.frame(stratum = "A", name = "y", mean = 1, sd = 1),
                     "y", 2e6)
  expect_identical(one$allocation, 2e6)
  expect_lte(one$search$bound, one$objective_value)
  ## Level 0.5 at (4, 2), cost 6, leaves one unit that the second stage
  ## spends on the second stratum.
  two <- whole_alone(
    data.frame(stratum = c("A", "B"), N = 4, unit_cost = 1),
    data.frame(stratum = c("A", "B", "A", "B"),
               name = c("y1", "y1", "y2", "y2"), mean = 1,
               sd = c(sqrt(2), 0, 0, 1)),
    data.frame(name = c("y1", "y2"), domain = "stratum", level = c("A", "B")),
    7
  )
  expect_identical(two$allocation, c(4, 3))
  expect_identical(two$search$second_stage, "optimal")
})

## Invariants

test_that("the worst rule keeps the invariants of its definition", {
  ref <- bethel_reference()
  components <- rownames(ref$d_qh)
  one <- worst_fit("Age")
  sum_one <- n_alloc_joint(ref$frame, ref$measures, objective = "Age",
                           budget = 241, fpc = "none")
  expect_equal(one$detail$n, sum_one$detail$n, tolerance = 1e-10)
  expect_equal(one$objective_value, sum_one$objective_value,
               tolerance = 1e-10)

  base <- worst_fit(components)
  scaled <- worst_fit(data.frame(name = components, priority = 3))
  expect_equal(scaled$objective_value, 3 * base$objective_value,
               tolerance = 1e-10)
  expect_equal(scaled$detail$n, base$detail$n, tolerance = 1e-10)
  expect_identical(scaled$detail$n_int, base$detail$n_int)

  reversed <- worst_fit(rev(components))
  expect_equal(reversed$objective_value, base$objective_value,
               tolerance = 1e-12)
  expect_equal(reversed$detail$n, base$detail$n, tolerance = 1e-10)

  ## At the cost of the equal-CV minimum-cost design, that design is the
  ## least worst one.
  cheapest <- n_alloc_joint(ref$frame, ref$measures,
                            targets = data.frame(name = components, cv = 0.06),
                            fpc = "none")
  at_cost <- worst_fit(components, budget = cheapest$params$achieved$cost)
  expect_equal(at_cost$objective_value, 0.06^2, tolerance = 1e-9)
  expect_equal(at_cost$detail$n, cheapest$detail$n, tolerance = 1e-8)
})

test_that("finite population cases keep the worst level exact", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(400, 120, 15),
                      unit_cost = c(1, 2, 3), region = c("n", "n", "s"))
  measures <- data.frame(
    stratum = rep(frame$stratum, 3),
    name = rep(c("y", "z", "y2"), each = 3),
    mean = c(10, 12, 9, 5, 4, 6, 20, 24, 18),
    sd = c(4, 6, 5, 2, 3, 2, 8, 12, 10)
  )
  ## y2 is y doubled, so its cap duplicates y's. The regional component
  ## overlaps the overall one.
  objective <- data.frame(name = c("y", "z", "y", "y2"),
                          domain = c(".overall", ".overall", "region",
                                     ".overall"),
                          level = c(NA, NA, "n", NA))
  run <- function(budget) {
    n_alloc_joint(frame, measures, objective = objective, budget = budget,
                  fpc = "unit", objective_rule = "worst")
  }
  fit <- run(300)
  expect_identical(fit$optimization$classification, "optimal")
  expect_equal(fit$params$achieved$cost, 300, tolerance = 1e-9)
  expect_equal(fit$objective$.relvar[1], fit$objective$.relvar[4],
               tolerance = 1e-12)
  expect_equal(fit$objective_value, max(fit$objective$.contribution))
  expect_lte(fit$operational$cost, 300)
  expect_lte(fit$optimization$search$gap, 1e-9)

  census <- run(700)
  expect_identical(census$detail$n, frame$N)
  expect_identical(census$objective_value, 0)
  expect_identical(census$operational$objective_value, 0)
})

test_that("the worst-rule slope is missing at the cheapest budget", {
  ref <- bethel_reference()
  targets <- data.frame(name = "Age", cv = 0.06)
  continuous <- function(budget) {
    svyplan:::.n_alloc_bethel(
      frame = ref$frame, measures = ref$measures, targets = targets,
      unit_cost = NULL, alpha = 0.05, deff = 1, resp_rate = 1,
      eligible_per_unit = 1, eligibility_rate = 1, df = NULL,
      min_n_stratum = NULL, objective = c("Floorspace", "Oil"),
      budget = budget, fpc = "none", .finite_prop_var = TRUE,
      .whole = FALSE, objective_rule = "worst"
    )
  }
  cheapest <- n_alloc_joint(ref$frame, ref$measures, targets = targets,
                            fpc = "none")
  fit <- continuous(cheapest$params$achieved$cost)
  expect_identical(fit$optimization$slopes,
                   "budget_at_minimum")
  expect_equal(fit$detail$n, cheapest$detail$n, tolerance = 1e-5)
  expect_true(is.na(fit$constraints$.objective_sensitivity))
})

## Outputs, assessment and round trip

test_that("a worst-rule fit assesses, round trips and predicts under its rule", {
  ref <- bethel_reference()
  fit <- worst_fit(rownames(ref$d_qh))
  expect_identical(fit$params$objective_rule, "worst")
  expect_identical(names(fit$objective),
                   c("component", "name", "domain", "level", "priority",
                     ".relvar", ".cv", ".contribution", ".binding"))

  continuous <- prec_alloc_joint(fit)
  expect_equal(continuous$objective_value, fit$objective_value,
               tolerance = 1e-12)
  expect_identical(continuous$params$objective_rule, "worst")
  whole <- prec_alloc_joint(fit, n = fit$detail$n_int)
  expect_equal(whole$objective_value, fit$operational$objective_value,
               tolerance = 1e-12)

  again <- n_alloc_joint(continuous)
  expect_identical(again$params$objective_rule, "worst")
  expect_equal(again$detail$n, fit$detail$n, tolerance = 1e-8)

  grid <- predict(fit, data.frame(budget = c(217, 241)))
  expect_equal(grid$objective_value[2], fit$objective_value,
               tolerance = 1e-10)
  expect_equal(grid$cv, sqrt(grid$objective_value))
  expect_equal(grid$cv[1], 0.06324933, tolerance = 1e-7)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(plot(fit, newdata = data.frame(budget = c(217, 260))),
                   fit)

  sum_rule <- n_alloc_joint(ref$frame, ref$measures, objective = "Age",
                            budget = 241, fpc = "none")
  expect_identical(sum_rule$params$objective_rule, "sum")
  expect_true(".share" %in% names(sum_rule$objective))
})

test_that("a worst-rule fit prints its rule at width 80", {
  old <- options(width = 80)
  on.exit(options(old))
  fit <- worst_fit(rownames(bethel_reference()$d_qh))
  printed <- capture.output(print(fit))
  expect_true(any(grepl(
    "question: least worst weighted relative variance within a budget of 241",
    printed, fixed = TRUE
  )))
  expect_true(any(grepl("^worst weighted rel. variance: ", printed)))
  expect_true(any(grepl("^worst cv = 0\\.0600 \\(", printed)))
  expect_true(any(grepl("\\*", printed)))
  expect_false(any(grepl(".share", printed, fixed = TRUE)))
  summarized <- capture.output(print(summary(fit)))
  expect_true(any(grepl(
    "Question: minimize the worst weighted relative variance within budget = 241",
    summarized, fixed = TRUE
  )))
  expect_true(any(grepl("the largest contribution is the objective",
                        summarized, fixed = TRUE)))
  expect_true(any(grepl("worst: Floorspace@.overall", summarized,
                        fixed = TRUE)))
  expect_true(any(grepl("sum-rule priorities giving the same continuous",
                        summarized, fixed = TRUE)))
  long <- c(printed, summarized)
  long <- long[nchar(long) > 80]
  expect_true(all(grepl("whole-unit search", long, fixed = TRUE)))
})

test_that("an affordable domain census reads exactly zero", {
  ## Random frames whose budget buys a census of the objective's domain
  ## and the floors elsewhere. Oracle: a census has no variance under unit
  ## FPC, so the level, its whole counterpart and the gap are all zero.
  set.seed(3)
  for (i in 1:40) {
    H <- sample(3:5, 1)
    frame <- data.frame(stratum = paste0("s", 1:H),
                        N = sample(3:9, H, TRUE),
                        unit_cost = round(runif(H, 0.3, 3.7), 3),
                        region = c("in", sample(c("in", "out"), H - 1, TRUE)))
    measures <- data.frame(stratum = frame$stratum, name = "y",
                           mean = runif(H, 1, 9), sd = runif(H, 0.3, 7))
    inside <- frame$region == "in"
    budget <- sum(frame$N[inside] * frame$unit_cost[inside]) +
      sum(frame$unit_cost[!inside]) * 1.0001
    fit <- n_alloc_joint(frame, measures,
                         objective = data.frame(name = "y", domain = "region",
                                                level = "in"),
                         budget = budget, fpc = "unit",
                         objective_rule = "worst")
    label <- sprintf("case %d", i)
    expect_identical(fit$objective_value, 0, label = label)
    expect_identical(fit$operational$objective_value, 0, label = label)
    expect_identical(fit$operational$status, "optimal", label = label)
    expect_identical(fit$optimization$search$gap, 0, label = label)
  }
})

test_that("the whole worst level is bought within the nominal budget", {
  ## Stratum A costs 1e7, so a budget tolerance of 1e-6 would buy ten units
  ## of the cheap strata. Oracle: every allocation of B and C within the
  ## budget and its rounding room, scored with the closed form below. At
  ## the cheapest cost only (1, 1, 1) fits.
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(1, 20, 20),
                      unit_cost = c(1e7, 6, 5))
  sd <- function(v) sqrt(v) * 41 / c(1, 20, 20)
  cases <- list(
    one = list(measures = data.frame(stratum = frame$stratum, name = "a",
                                     mean = 1, sd = sd(c(0, 1.02, 1))),
               objective = "a",
               level = function(B, C) 1.02 / B + 1 / C),
    two = list(measures = data.frame(stratum = rep(frame$stratum, 2),
                                     name = rep(c("a", "b"), each = 3),
                                     mean = 1,
                                     sd = c(sd(c(0, 1.02, 1)), sd(c(0, 1.35, 0)))),
               objective = c("a", "b"),
               level = function(B, C) pmax(1.02 / B + 1 / C, 1.35 / B))
  )
  for (budget in c(10000011, 10000017, 10000027)) {
    grid <- expand.grid(B = 1:20, C = 1:20)
    grid <- grid[1e7 + 6 * grid$B + 5 * grid$C <= budget * (1 + 1e-9), ]
    for (name in names(cases)) {
      z <- cases[[name]]
      label <- paste(name, format(budget))
      fit <- n_alloc_joint(frame, z$measures, objective = z$objective,
                           budget = budget, fpc = "none",
                           objective_rule = "worst")
      expect_equal(fit$operational$objective_value,
                   min(z$level(grid$B, grid$C)), tolerance = 1e-10,
                   label = label)
      expect_lte(fit$operational$cost, budget * (1 + 1e-9))
      expect_lte(fit$params$achieved$cost, budget * (1 + 1e-9))
      expect_lte(fit$optimization$search$bound,
                 fit$operational$objective_value)
      expect_identical(fit$operational$status, "optimal", label = label)
      if (budget == 10000011) {
        expect_equal(fit$detail$n_int, c(1, 1, 1), label = label)
      }
    }
  }
})

test_that("a kink in the worst level withholds the budget slope", {
  ## Both components bind with stratum A at its upper bound, so the one
  ## free stratum cannot separate their prices. Oracle: refits on each side
  ## of the budget give different slopes.
  W <- c(2, 98) / 100
  frame <- data.frame(stratum = c("A", "B"), N = c(2, 98), unit_cost = 1)
  measures <- data.frame(stratum = c("A", "B", "A", "B"),
                         name = c("y1", "y1", "y2", "y2"), mean = 1,
                         sd = c(1, 2, sqrt(3), sqrt(2)) / rep(W, 2))
  fit <- function(budget) {
    n_alloc_joint(frame, measures, objective = c("y1", "y2"), budget = budget,
                  fpc = "none", objective_rule = "worst")
  }
  at <- fit(4)
  h <- 1e-5
  left <- (at$objective_value - fit(4 - h)$objective_value) / h
  right <- (fit(4 + h)$objective_value - at$objective_value) / h
  expect_gt(abs(left - right), 0.1)
  expect_identical(at$objective$.binding, c(TRUE, TRUE))
  expect_identical(at$optimization$slopes, "not_identifiable")
  expect_true(is.na(at$optimization$budget_sensitivity))
})

test_that("a census beside a low-variance stratum keeps the small variance", {
  ## Oracle: the relative variance of the sampled stratum,
  ## W_B^2 S_B^2 (1 / 10 - 1 / 100), set to 2.25e-17. The census stratum
  ## contributes nothing.
  sd_b <- sqrt(2.25e-17 / ((100 / 110)^2 * (1 / 10 - 1 / 100)))
  frame <- data.frame(stratum = c("A", "B"), N = c(10, 100))
  ## The census stratum's large SD puts an aggregate cancellation threshold
  ## near 1e-11, which would erase the sampled stratum's variance.
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 1,
                         sd = c(1000, sd_b))
  for (rule in c("sum", "worst")) {
    assessed <- prec_alloc_joint(frame, n = c(10, 10), measures = measures,
                                 objective = "y", fpc = "unit",
                                 objective_rule = rule)
    expect_equal(1e17 * assessed$objective$.relvar, 2.25,
                 tolerance = 1e-9, label = rule)
    expect_gt(assessed$objective$.relvar, 0)
  }
})

test_that("a trial with a large multiplier still certifies", {
  ## Two components, one stratum at its population size at the optimum.
  ## The binding cap carries a multiplier near 2400, so the gradient test
  ## once stopped the dual search short of the duality-gap tolerance.
  ## Oracle: the continuous level from the one-variable profile along the
  ## budget line, and the whole level by enumerating every design.
  fixtures <- list(
    list(N = c(6, 17), cost = c(7, 7), budget = 136.782517130766,
         A = matrix(c(1.78879471676428, 0.042740524327607,
                      25.1005075579283, 0.0493505668191266), 2)),
    list(N = c(7, 15), cost = c(5, 8), budget = 93.0452435655519,
         A = matrix(c(5.99395408198966, 0.313897359540511,
                      46.2855485857171, 0.0339708626482544), 2))
  )
  for (x in fixtures) {
    frame <- data.frame(stratum = c("A", "B"), N = x$N, unit_cost = x$cost)
    measures <- data.frame(stratum = rep(frame$stratum, 2),
                           name = rep(c("a", "b"), each = 2), mean = 1,
                           sd = sqrt(c(x$A)) * sum(x$N) / x$N)
    fit <- n_alloc_joint(frame, measures, objective = c("a", "b"),
                         objective_rule = "worst", budget = x$budget,
                         fpc = "none")
    value <- function(n1, n2) max(x$A[1, ] / n1 + x$A[2, ] / n2)
    profile <- function(n1) {
      value(n1, min(x$N[2], (x$budget - x$cost[1] * n1) / x$cost[2]))
    }
    lo <- max(1, (x$budget - x$N[2] * x$cost[2]) / x$cost[1])
    hi <- min(x$N[1], (x$budget - x$cost[2]) / x$cost[1])
    continuous <- min(profile(lo), profile(hi),
                      stats::optimize(profile, c(lo, hi), tol = 1e-11)$objective)
    grid <- expand.grid(n1 = seq_len(x$N[1]), n2 = seq_len(x$N[2]))
    grid <- grid[grid$n1 * x$cost[1] + grid$n2 * x$cost[2] <=
                   x$budget * (1 + 1e-9), ]
    whole <- min(mapply(value, grid$n1, grid$n2))
    expect_identical(fit$optimization$classification, "optimal")
    expect_equal(fit$objective_value, continuous, tolerance = 1e-6)
    expect_identical(fit$operational$status, "optimal")
    expect_equal(fit$operational$objective_value, whole, tolerance = 1e-9)
  }
})

test_that("a budget just above the cheapest cost buys the worst optimum", {
  ## Stratum A has one unit at cost 1e7, so the budget tolerance covers ten
  ## units of the other strata. Oracle: component a binds alone, so the
  ## optimum is its own allocation, n_B / n_C = sqrt((1.02 / 6) / (1 / 5)),
  ## along the budget line, where component b stays below it.
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(1, 20, 20),
                      unit_cost = c(1e7, 6, 5))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("a", "b"), each = 3), mean = 1,
                         sd = sqrt(c(0, 1.02, 1, 0, 1.35, 0)) * 41 /
                           rep(frame$N, 2))
  r <- sqrt((1.02 / 6) / (1 / 5))
  for (budget in c(10000017, 10000021)) {
    n_c <- (budget - 1e7) / (6 * r + 5)
    n_b <- r * n_c
    expect_lt(1.35 / n_b, 1.02 / n_b + 1 / n_c)
    fit <- n_alloc_joint(frame, measures, objective = c("a", "b"),
                         budget = budget, fpc = "none",
                         objective_rule = "worst")
    expect_equal(fit$objective_value, 1.02 / n_b + 1 / n_c,
                 tolerance = 1e-5, label = format(budget))
    expect_lte(fit$params$achieved$cost, budget * (1 + 1e-9))
  }
})

test_that("a trial that fits the budget but fails certification keeps the search going", {
  ## Every worst-level trial whose design costs at most 0.9 of the budget
  ## is reported uncertified, with its design kept. Such a design meets the
  ## caps and fits, so it proves its level affordable. Oracle: the
  ## continuous level from the one-variable profile along the budget line.
  N <- c(6, 17)
  cost <- c(7, 7)
  budget <- 136.782517130766
  A <- matrix(c(1.78879471676428, 0.042740524327607,
                25.1005075579283, 0.0493505668191266), 2)
  frame <- data.frame(stratum = c("A", "B"), N = N, unit_cost = cost)
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("a", "b"), each = 2), mean = 1,
                         sd = sqrt(c(A)) * sum(N) / N)
  real <- .bethel_solve
  flipped <- 0L
  local_mocked_bindings(.bethel_solve = function(A, bound, cost, ...,
                                                 constraint_ids) {
    s <- real(A, bound, cost, ..., constraint_ids = constraint_ids)
    trial <- any(startsWith(constraint_ids, ".worst:")) &&
      !any(constraint_ids == ".objective")
    if (trial && identical(s$classification, "optimal") &&
        sum(cost * s$allocation) <= 0.9 * budget) {
      flipped <<- flipped + 1L
      s$classification <- "failed"
    }
    s
  })
  fit <- n_alloc_joint(frame, measures, objective = c("a", "b"),
                       objective_rule = "worst", budget = budget,
                       fpc = "none")
  profile <- function(n1) {
    max(A[1, ] / n1 + A[2, ] / min(N[2], (budget - cost[1] * n1) / cost[2]))
  }
  lo <- max(1, (budget - N[2] * cost[2]) / cost[1])
  hi <- min(N[1], (budget - cost[2]) / cost[1])
  continuous <- min(profile(lo), profile(hi),
                    stats::optimize(profile, c(lo, hi), tol = 1e-11)$objective)
  expect_gt(flipped, 0L)
  expect_equal(fit$objective_value, continuous, tolerance = 1e-6)
})
