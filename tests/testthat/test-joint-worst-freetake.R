# The worst rule on the free-take route of n_alloc_joint()

## Two strata whose whole designs can all be listed: PSU counts 1 to 6,
## takes 1 to 3 and 1 to 4, 432 designs in all.
freetake_fixture <- function() {
  list(
    frame = data.frame(stratum = c("A", "B"), N = c(18, 24), N_1 = c(6, 6),
                       cost_1 = c(12, 8), cost_2 = c(2, 3)),
    measures = data.frame(
      stratum = rep(c("A", "B"), 2), name = rep(c("cov", "inc"), each = 2),
      p = c(0.35, 0.45, NA, NA), mean = c(NA, NA, 100, 120),
      sd = c(NA, NA, 25, 30), icc_1 = c(0.34, 0.08, 0.04, 0.42)
    ),
    objective = data.frame(name = c("cov", "inc"), priority = c(1, 28))
  )
}

freetake_fit <- function(rule = "worst", cv = 0.26, budget = 190,
                         fpc = "none", ...) {
  z <- freetake_fixture()
  n_alloc_joint(z$frame, z$measures,
                targets = data.frame(name = "cov", cv = cv),
                objective = z$objective, budget = budget, fpc = fpc,
                take_2_min = 1, take_2_max = c(3, 4),
                objective_rule = rule, ...)
}

test_that("the free-take worst rule matches the enumerated whole designs", {
  ## Oracle: the 432 whole designs scored with the closed-form stratified
  ## two-stage relative variance sum_h W_h^2 S_h^2 (1 + (m_h - 1) icc_h) /
  ## (x_h m_h) / Ybar^2, under fpc = "none".
  z <- freetake_fixture()
  W <- z$frame$N / sum(z$frame$N)
  S2 <- list(cov = c(0.35, 0.45) * (1 - c(0.35, 0.45)), inc = c(25, 30)^2)
  mu <- list(cov = c(0.35, 0.45), inc = c(100, 120))
  icc <- list(cov = c(0.34, 0.08), inc = c(0.04, 0.42))
  relvar <- function(j, x, m) {
    sum(W^2 * S2[[j]] * (1 + (m - 1) * icc[[j]]) / (x * m)) /
      sum(W * mu[[j]])^2
  }
  grid <- expand.grid(x1 = 1:6, x2 = 1:6, m1 = 1:3, m2 = 1:4)
  grid$cost <- with(grid, 12 * x1 + 8 * x2 + 2 * x1 * m1 + 3 * x2 * m2)
  grid$cov <- with(grid, mapply(function(a, b, c, d) relvar("cov", c(a, b), c(c, d)),
                                x1, x2, m1, m2))
  grid$inc <- with(grid, mapply(function(a, b, c, d) relvar("inc", c(a, b), c(c, d)),
                                x1, x2, m1, m2))
  ok <- grid$cost <= 190 & sqrt(grid$cov) <= 0.26
  level <- pmax(grid$cov, 28 * grid$inc)[ok]
  best <- grid[ok, ][which.min(level), ]

  fit <- freetake_fit("worst")
  expect_equal(fit$detail$n_1_int, c(best$x1, best$x2))
  expect_equal(fit$detail$take_2_int, c(best$m1, best$m2))
  expect_identical(fit$detail$n_1_int, c(4, 6))
  expect_identical(fit$detail$take_2_int, c(3, 3))
  expect_identical(fit$operational$cost, 174)
  expect_equal(fit$operational$objective_value, min(level), tolerance = 1e-10)
  expect_equal(fit$operational$objective_value, 0.09104372, tolerance = 1e-7)
  expect_identical(fit$operational$status, "optimal")
  expect_lte(fit$optimization$search$bound, fit$operational$objective_value)

  ## The sum rule keeps its own design.
  summed <- freetake_fit("sum")
  expect_identical(summed$detail$n_1_int, c(5, 5))
  expect_identical(summed$detail$take_2_int, c(3, 4))
})

test_that("the continuous free-take level is certified by its dual bound", {
  fit <- freetake_fit("worst")
  opt <- fit$optimization
  expect_identical(opt$classification, "optimal")
  expect_equal(fit$objective_value, 0.0859138222, tolerance = 1e-8)
  expect_gte(opt$objective_lower_bound, 0.0859137)
  expect_lte(opt$objective_lower_bound, fit$objective_value)
  expect_equal(fit$params$achieved$cost, 190, tolerance = 1e-8)
  expect_identical(fit$objective$.binding, c(FALSE, TRUE))
})

test_that("free takes reproduce fixed takes and beat every fixed take grid", {
  z <- freetake_fixture()
  worst_at <- function(frame, ...) {
    n_alloc_joint(frame, z$measures, objective = z$objective, budget = 190,
                  objective_rule = "worst", ...)$objective_value
  }
  for (fpc in c("none", "stage", "unit")) {
    held <- worst_at(z$frame, fpc = fpc, take_2_min = c(3, 4),
                     take_2_max = c(3, 4))
    fixed <- worst_at(transform(z$frame, take_2 = c(3, 4)), fpc = fpc)
    expect_equal(held, fixed, tolerance = 1e-6, label = fpc)

    ## Oracle: the continuous fit at every whole take pair. Free takes may
    ## exceed the best of them only within the certified gap of 1e-6.
    free <- worst_at(z$frame, fpc = fpc, take_2_min = 1, take_2_max = c(3, 4))
    grid <- expand.grid(m1 = 1:3, m2 = 1:4)
    each <- mapply(function(a, b) {
      worst_at(transform(z$frame, take_2 = c(a, b)), fpc = fpc)
    }, grid$m1, grid$m2)
    expect_lte(free, min(each) * (1 + 1e-6), label = fpc)
  }
})

test_that("free-take slopes need an optimal fit and verified multipliers", {
  ## D4: a slope is reported only for an "optimal" fit whose multipliers
  ## price the full free-take problem. Oracle for the reported slope:
  ## centered refits of the worst level.
  fit <- freetake_fit("worst", cv = 0.245)
  expect_identical(fit$optimization$classification, "optimal")
  expect_true(fit$constraints$.binding)
  expect_identical(fit$optimization$slopes, "available")
  h <- 1e-7
  moved <- (freetake_fit("worst", cv = 0.245 + h)$objective_value -
              freetake_fit("worst", cv = 0.245 - h)$objective_value) / (2 * h)
  expect_equal(fit$constraints$.objective_sensitivity, moved, tolerance = 1e-4)
  up <- freetake_fit("worst", cv = 0.245, budget = 190 + 1e-4)
  expect_equal(fit$optimization$budget_sensitivity,
               (up$objective_value - fit$objective_value) / 1e-4,
               tolerance = 1e-3)

  for (rule in c("sum", "worst")) {
    feasible <- freetake_fit(rule, cv = 0.243)
    expect_identical(feasible$optimization$classification, "feasible",
                     label = rule)
    expect_identical(feasible$optimization$slopes,
                     "not_optimal", label = rule)
    expect_true(is.na(feasible$constraints$.objective_sensitivity))
    expect_true(is.na(feasible$optimization$budget_sensitivity))
    expect_match(
      gsub("\\s+", " ", paste(capture.output(print(summary(feasible))),
                              collapse = " ")),
      "feasible but not certified optimal", fixed = TRUE
    )
  }
  slack <- freetake_fit("worst", cv = 0.26)
  expect_identical(slack$constraints$.objective_sensitivity, 0)
})

test_that("free-take slopes come only from multipliers the full dual verifies", {
  ## A fixed-take repair can replace the solver's multipliers with ones that
  ## price only the repaired takes. The check takes the dual's vector when
  ## the solver's fails, and reports none when both fail.
  fit <- freetake_fit("worst", cv = 0.245)
  z <- freetake_fixture()
  spec <- svyplan:::.resolve_freetake_spec(z$frame, z$measures,
                                           take_2_min = 1,
                                           take_2_max = c(3, 4))
  prob <- svyplan:::.freetake_build_problem(
    frame = z$frame, measures = z$measures,
    targets = data.frame(name = "cov", cv = 0.245), ft_spec = spec,
    objective = z$objective, fpc = "none"
  )
  control <- svyplan:::.bethel_control()
  solved <- svyplan:::.freetake_worst_solve(prob, 190, control)
  caps <- svyplan:::.freetake_worst_caps(prob)
  cols <- list(a = caps$a, b = caps$b, d = caps$d, e = caps$e,
               vmax = rep(solved$objective_bound, caps$J))
  slopes <- function(s) {
    svyplan:::.freetake_slopes(prob, s, fit$constraints, cols,
                                         control)
  }
  reference <- slopes(solved)
  expect_identical(reference$status, "available")
  repaired <- solved
  repaired$solved$lambda <- 1.1 * repaired$solved$lambda
  expect_identical(slopes(repaired)$status, "available")
  expect_equal(slopes(repaired)$objective, reference$objective,
               tolerance = 1e-10)
  repaired$solved$dual_lambda <- 1.1 * repaired$solved$dual_lambda
  expect_identical(slopes(repaired)$status, "multipliers_not_verified")
  expect_true(is.na(slopes(repaired)$objective))
  expect_true(is.na(slopes(repaired)$sensitivity))
})

test_that("the free-take dual weights are sum-rule priorities giving the same design", {
  ## At priorities 1 and 18.5 both components sit at the level, so both
  ## weights are positive.
  z <- freetake_fixture()
  run <- function(priority, rule) {
    n_alloc_joint(z$frame, z$measures,
                  objective = data.frame(name = c("cov", "inc"),
                                         priority = priority),
                  budget = 190, fpc = "none", take_2_min = 1,
                  take_2_max = c(3, 4), objective_rule = rule)
  }
  fit <- run(c(1, 18.5), "worst")
  expect_identical(fit$objective$.binding, c(TRUE, TRUE))
  weights <- fit$optimization$objective_weights
  expect_equal(sum(weights), 1, tolerance = 1e-12)
  expect_true(all(weights > 0.1))
  summed <- run(unname(weights), "sum")
  expect_equal(summed$detail$n_1, fit$detail$n_1, tolerance = 1e-4)
  expect_equal(summed$detail$take_2, fit$detail$take_2, tolerance = 1e-4)
})

test_that("a census floor under unit FPC keeps the free-take search going", {
  ## The first trial above the census floor has no sign at this cap, so the
  ## lower end is raised rather than taken as the answer.
  fit <- freetake_fit("worst", cv = 0.27, fpc = "unit")
  expect_identical(fit$optimization$classification, "optimal")
  expect_equal(fit$params$achieved$cost, 190, tolerance = 1e-6)
})

test_that("the free-take slope is missing at the cheapest budget", {
  z <- freetake_fixture()
  targets <- data.frame(name = "cov", cv = 0.24)
  spec <- svyplan:::.resolve_freetake_spec(z$frame, z$measures,
                                           take_2_min = 1,
                                           take_2_max = c(3, 4))
  prob <- svyplan:::.freetake_build_problem(
    frame = z$frame, measures = z$measures, targets = targets,
    ft_spec = spec, objective = z$objective, fpc = "none"
  )
  cheapest <- svyplan:::.freetake_budget_open(prob, 1e9,
                                              svyplan:::.bethel_control())
  fit <- svyplan:::.freetake_worst_solve(prob, cheapest$cost_min)
  expect_true(fit$at_minimum)
  caps <- svyplan:::.freetake_worst_caps(prob)
  constraints <- data.frame(.binding = TRUE, .metric = "cv", .target = 0.24)
  slopes <- svyplan:::.freetake_slopes(
    prob, fit, constraints,
    list(a = caps$a, b = caps$b, d = caps$d, e = caps$e,
         vmax = rep(fit$objective_bound, caps$J)),
    svyplan:::.bethel_control()
  )
  expect_identical(slopes$status, "budget_at_minimum")
  expect_true(is.na(slopes$objective))
  expect_true(is.na(slopes$sensitivity))
})

test_that("free-take rows lose their own prices only in a dependency", {
  shared <- svyplan:::.freetake_shared_rows
  prob <- list(H = 2L, N_1 = c(10, 10), N = c(100, 100), L_ult = c(1, 1),
               ft_spec = list(take_2_min = c(1, 1), take_2_max_cont = c(9, 9)))
  a <- matrix(c(1, 2, 3, 1), 2)
  b <- matrix(c(4, 1, 1, 5), 2)
  d <- matrix(0, 2, 2)
  expect_identical(shared(prob, a, b, d, x = c(4, 5), y = c(20, 25)),
                   c(FALSE, FALSE))
  expect_identical(shared(prob, cbind(a, a[, 1] * 3), cbind(b, b[, 1] * 3),
                          cbind(d, 0), x = c(4, 5), y = c(20, 25)),
                   c(TRUE, FALSE, TRUE))
  ## Every x at its upper bound and every y at its take ceiling leaves no
  ## free direction.
  expect_true(all(shared(prob, a, b, d, x = c(10, 10), y = c(90, 90))))
})

test_that("the free-take second stage spends what the worst level leaves", {
  ## Stratum A's component reaches its floor, the whole of A, for less than
  ## the budget, so the level binds there and B gets the rest.
  frame <- data.frame(stratum = c("A", "B"), N = c(6, 400), N_1 = c(2, 40),
                      cost_1 = 10, cost_2 = 1)
  measures <- data.frame(stratum = c("A", "B", "A", "B"),
                         name = c("y1", "y1", "y2", "y2"),
                         mean = 1, sd = c(1, 0, 0, 0.5), icc_1 = 0.1)
  objective <- data.frame(name = c("y1", "y2"), domain = "stratum",
                          level = c("A", "B"))
  fit <- n_alloc_joint(frame, measures, objective = objective, budget = 300,
                       fpc = "none", take_2_min = 1, take_2_max = c(3, 10),
                       objective_rule = "worst")
  expect_identical(fit$optimization$budget_binding, FALSE)
  expect_true(fit$optimization$second_stage %in% c("optimal", "feasible"))
  expect_equal(fit$params$achieved$cost, 300, tolerance = 1e-6)
  expect_identical(fit$objective$.binding, c(TRUE, FALSE))
  expect_lt(fit$objective$.contribution[2], fit$objective$.contribution[1])
})

test_that("the free-take whole second stage spends the leftover", {
  ## Fixed takes of 1, so the counts alone set the components 2 / x1 and
  ## 1 / x2 with x in [1, 4] at unit cost. Level 0.5 is reached at (4, 2),
  ## cost 6, and at (4, 3), cost 7. Called without the rounded start, the
  ## bisection returns the cheaper one and the second stage takes (4, 3).
  frame <- data.frame(stratum = c("A", "B"), N = 4, N_1 = 4, cost_1 = 0.5,
                      cost_2 = 0.5)
  measures <- data.frame(stratum = c("A", "B", "A", "B"),
                         name = c("y1", "y1", "y2", "y2"), mean = 1,
                         sd = c(sqrt(2), 0, 0, 1), icc_1 = 0.1)
  objective <- data.frame(name = c("y1", "y2"), domain = "stratum",
                          level = c("A", "B"))
  spec <- svyplan:::.resolve_freetake_spec(frame, measures, take_2_min = 1,
                                           take_2_max = 1)
  prob <- svyplan:::.freetake_build_problem(
    frame = frame, measures = measures, targets = NULL, ft_spec = spec,
    objective = objective, fpc = "none"
  )
  whole <- svyplan:::.freetake_whole_worst(prob, 7,
                                           svyplan:::.bethel_control(), 5000)
  expect_equal(whole$n_1_int, c(4, 3))
  expect_identical(whole$search$second_stage, "optimal")
})

test_that("mixed, stepped and rated free takes match their enumeration", {
  ## Stratum A has free takes 1 to 3, stratum B takes 2, 4 or 6, with a PSU
  ## response rate and an objective whose domain overlaps the overall one.
  ## Oracle: every whole design within the budget, each assessed at its
  ## fixed takes.
  frame <- data.frame(stratum = c("A", "B"), N = c(30, 72), N_1 = c(5, 6),
                      cost_1 = c(15, 10), cost_2 = c(2, 3),
                      resp_rate_1 = c(0.9, 0.8), region = c("n", "s"))
  measures <- data.frame(stratum = rep(c("A", "B"), 2),
                         name = rep(c("cov", "inc"), each = 2),
                         p = c(0.3, 0.5, NA, NA), mean = c(NA, NA, 80, 95),
                         sd = c(NA, NA, 20, 35), icc_1 = c(0.2, 0.1, 0.05, 0.3))
  objective <- data.frame(name = c("cov", "inc", "inc"),
                          domain = c(".overall", ".overall", "region"),
                          level = c(NA, NA, "n"), priority = c(1, 20, 5))
  budget <- 160
  for (fpc in c("none", "unit")) {
    fit <- n_alloc_joint(frame, measures, objective = objective,
                         budget = budget, fpc = fpc, take_2_min = c(1, 2),
                         take_2_max = c(3, 6), take_2_step = c(1, 2),
                         objective_rule = "worst")
    grid <- expand.grid(x1 = 1:5, x2 = 1:6, m1 = 1:3, m2 = c(2, 4, 6))
    grid <- grid[with(grid, (15 + 2 * m1) * x1 + (10 + 3 * m2) * x2) <= budget &
                   grid$x1 * grid$m1 <= 30 & grid$x2 * grid$m2 <= 72, ]
    level <- mapply(function(a, b, c, d) {
      prec_alloc_joint(transform(frame, take_2 = c(c, d)), n = c(a * c, b * d),
                       measures = measures, objective = objective,
                       fpc = fpc, objective_rule = "worst")$objective_value
    }, grid$x1, grid$x2, grid$m1, grid$m2)
    expect_equal(fit$operational$objective_value, min(level),
                 tolerance = 1e-10, label = fpc)
    expect_lte(fit$optimization$search$bound, min(level) * (1 + 1e-12))
  }
})

test_that("a free-take domain census is certified at exact zero", {
  ## The objective covers region A only, whose census is affordable. The
  ## census cancels the variance only to rounding, and the level must still
  ## read zero with a closed certificate.
  frame <- data.frame(stratum = c("A", "B"), region = c("A", "B"),
                      N = c(20, 100), N_1 = c(4, 10),
                      cost_1 = c(10, 100), cost_2 = c(1, 10))
  measures <- data.frame(stratum = "A", name = "y", mean = 1, sd = 1,
                         icc_1 = 0.2)
  fit <- n_alloc_joint(frame, measures,
                       objective = data.frame(name = "y", domain = "region",
                                              level = "A"),
                       budget = 175, take_2_max = 5, objective_rule = "worst")
  expect_identical(fit$optimization$classification, "optimal")
  expect_identical(fit$operational$status, "optimal")
  expect_identical(fit$objective_value, 0)
  expect_identical(fit$operational$objective_value, 0)
  expect_equal(fit$operational$cost, 170)
})

test_that("a free-take census and a node limit keep their records", {
  z <- freetake_fixture()
  census <- n_alloc_joint(z$frame, z$measures, objective = z$objective,
                          budget = 240, fpc = "unit", take_2_min = 1,
                          take_2_max = c(3, 4), objective_rule = "worst")
  expect_identical(census$objective_value, 0)
  expect_identical(census$operational$objective_value, 0)
  expect_identical(census$detail$n_1_int * census$detail$take_2_int,
                   z$frame$N)

  capped <- freetake_fit("worst", control = list(max_nodes = 2))
  search <- capped$optimization$search
  expect_identical(capped$operational$status, "limit")
  expect_lte(search$nodes, 2)
  expect_lte(search$bound, capped$operational$objective_value)
  expect_true(capped$operational$all_pass)
  expect_lte(capped$operational$cost, 190 * (1 + 1e-6))
})

test_that("a larger free-take frame solves within the default nodes", {
  set.seed(5)
  H <- 12
  N_1 <- sample(20:60, H, replace = TRUE)
  frame <- data.frame(stratum = sprintf("s%02d", 1:H), N_1 = N_1,
                      N = N_1 * sample(40:120, H, replace = TRUE),
                      cost_1 = round(runif(H, 40, 150)),
                      cost_2 = round(runif(H, 3, 12)),
                      region = rep(c("r1", "r2", "r3"), each = 4))
  names <- c("a", "b", "c")
  measures <- data.frame(stratum = rep(frame$stratum, 3),
                         name = rep(names, each = H),
                         mean = runif(3 * H, 5, 50), sd = runif(3 * H, 2, 30),
                         icc_1 = runif(3 * H, 0.01, 0.2))
  objective <- rbind(
    data.frame(name = names, domain = ".overall", level = NA),
    data.frame(name = "a", domain = "region", level = c("r1", "r2", "r3"))
  )
  fit <- n_alloc_joint(frame, measures, objective = objective, budget = 8000,
                       fpc = "unit", take_2_min = 2, take_2_max = 25,
                       objective_rule = "worst")
  search <- fit$optimization$search
  expect_identical(fit$optimization$classification, "optimal")
  expect_lte(search$bound, fit$operational$objective_value)
  expect_lte(search$gap, 0.015)
  expect_lte(fit$operational$cost, 8000 * (1 + 1e-6))
  ## Regression values: the rounded start and the half-nodes trials each
  ## lower the whole level here, by 7.5e-4 and 3.7e-3 at 8000 and by 1.3e-3
  ## and 3.3e-3 at 20000 (relative).
  for (case in list(c(8000, 0.0038508423), c(20000, 0.0013382616))) {
    staged <- n_alloc_joint(frame, measures, objective = objective,
                            budget = case[1], fpc = "stage", take_2_min = 2,
                            take_2_max = 25, objective_rule = "worst")
    expect_lte(staged$operational$objective_value, case[2] * (1 + 1e-6))
  }
})

test_that("a free-take worst fit assesses, round trips and predicts", {
  fit <- freetake_fit("worst")
  continuous <- prec_alloc_joint(fit)
  expect_identical(continuous$params$objective_rule, "worst")
  expect_equal(continuous$objective_value, fit$objective_value,
               tolerance = 1e-10)
  whole <- prec_alloc_joint(fit, n = fit$detail$n_int)
  expect_equal(whole$objective_value, fit$operational$objective_value,
               tolerance = 1e-12)
  expect_equal(whole$params$achieved$cost, fit$operational$cost)
  again <- n_alloc_joint(continuous)
  expect_identical(again$params$objective_rule, "worst")
  expect_equal(again$objective_value, fit$objective_value, tolerance = 1e-6)
  grid <- predict(fit, data.frame(budget = c(180, 190)))
  expect_equal(grid$objective_value[2], fit$objective_value,
               tolerance = 1e-8)
  old <- options(width = 80)
  on.exit(options(old))
  printed <- capture.output(print(fit), print(summary(fit)))
  expect_true(any(grepl("least worst weighted relative variance", printed,
                        fixed = TRUE)))
  expect_true(all(nchar(printed) <= 80 |
                    grepl("whole-unit search|field design search", printed)))
})

test_that("a withheld free-take fit reports no slope on a slack row either", {
  ## The returned design is only feasible, so a target slack there may bind
  ## at the optimum. Every row is NA, zero included.
  z <- freetake_fixture()
  for (rule in c("sum", "worst")) {
    fit <- n_alloc_joint(
      z$frame, z$measures,
      targets = data.frame(name = c("cov", "inc"), cv = c(0.243, 0.2)),
      objective = z$objective, budget = 190, fpc = "none",
      take_2_min = 1, take_2_max = c(3, 4), objective_rule = rule
    )
    expect_identical(fit$optimization$slopes, "not_optimal", label = rule)
    expect_identical(fit$constraints$.binding, c(TRUE, FALSE), label = rule)
    expect_true(all(is.na(fit$constraints$.sensitivity)), label = rule)
    expect_true(all(is.na(fit$constraints$.objective_sensitivity)),
                label = rule)
  }
})

test_that("a budget just above the cheapest free-take cost buys the optimum", {
  ## Stratum C is one PSU of one unit with no variance, so its cost only
  ## shifts the budget. At cost 1e7 the budget tolerance covers ten units
  ## of the other strata. Oracle: the same spend above the cheapest cost
  ## buys the same level when C costs 1.
  z <- freetake_fixture()
  frame <- rbind(z$frame, data.frame(stratum = "C", N = 1, N_1 = 1,
                                     cost_1 = 1, cost_2 = 1))
  measures <- rbind(z$measures, data.frame(
    stratum = "C", name = c("cov", "inc"), p = c(0, NA), mean = c(NA, 100),
    sd = c(NA, 0), icc_1 = 0.1
  ))
  targets <- data.frame(name = "cov", cv = 0.26)
  cheapest <- n_alloc_joint(frame, measures, targets = targets,
                            take_2_min = 1, take_2_max = c(3, 4, 1),
                            fpc = "none")$params$achieved$cost
  fit <- function(cost_c, budget) {
    frame$cost_1[3] <- cost_c
    n_alloc_joint(frame, measures, targets = targets,
                  objective = z$objective, budget = budget, fpc = "none",
                  take_2_min = 1, take_2_max = c(3, 4, 1),
                  objective_rule = "worst")
  }
  cheap <- fit(1, cheapest + 9)
  dear <- fit(1e7, cheapest + 9 - 1 + 1e7)
  expect_identical(cheap$optimization$classification, "optimal")
  expect_identical(dear$optimization$classification, "optimal")
  expect_equal(dear$objective_value, cheap$objective_value, tolerance = 1e-5)
})
