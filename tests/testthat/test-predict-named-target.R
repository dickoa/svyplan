## predict() grids that name one target per row

named_inputs <- function() {
  frame <- data.frame(
    stratum = c("North urban", "North rural", "South urban"),
    region = c("North", "North", "South"),
    N = c(1000, 1800, 1200)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("coverage", "income"), each = 3),
    p = c(0.5, 0.4, 0.6, NA, NA, NA),
    mean = c(NA, NA, NA, 50, 55, 60),
    sd = c(NA, NA, NA, 10, 12, 15)
  )
  targets <- data.frame(name = c("coverage", "income"),
                        domain = c(".overall", "region"),
                        cv = c(0.06, NA), moe = c(NA, 3))
  list(frame = frame, measures = measures, targets = targets)
}

## Oracle for each row: n_alloc_joint() refitted with that one target scaled.
expect_named_rows <- function(fit, grid, refit) {
  checked <- character(0)
  for (i in seq_len(nrow(grid))) {
    id <- grid$target_constraint[i]
    label <- sprintf("row %d %s", i, id %||% "all")
    direct <- refit(id, grid$target_mult[i])
    expect_equal(grid$n[i], direct$n, tolerance = 1e-8, label = label)
    expect_equal(grid$cost_int[i], direct$operational$cost, label = label)
    expect_equal(grid$n_int[i], direct$operational$n, label = label)
    expect_identical(grid$.status[i], direct$operational$status, label = label)
    if (!is.na(id)) {
      k <- match(id, direct$constraints$constraint)
      expect_equal(grid$.target[i], direct$constraints$.target[k],
                   label = label)
      expect_equal(grid$.achieved[i], direct$constraints$.achieved[k],
                   tolerance = 1e-10, label = label)
      expect_equal(grid$.achieved_int[i],
                   direct$operational$constraints$.achieved[k],
                   tolerance = 1e-10, label = label)
    }
    checked <- c(checked, label)
  }
  checked
}

test_that("a named row equals a refit with that target scaled", {
  z <- named_inputs()
  checked <- character(0)
  for (rule in c("sum", "worst")) {
    objective <- if (rule == "sum") "income" else c("coverage", "income")
    fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                         objective = objective, budget = 400,
                         objective_rule = rule)
    newdata <- data.frame(
      target_constraint = c("coverage@.overall:cv", "income@region=South:moe",
                            NA),
      target_mult = c(0.95, 1.2, 0.98)
    )
    grid <- predict(fit, newdata)
    refit <- function(id, m) {
      targets <- fit$params$targets
      if (is.na(id)) {
        targets$cv <- targets$cv * m
        targets$moe <- targets$moe * m
      } else {
        k <- match(id, targets$constraint)
        col <- if (is.na(targets$cv[k])) "moe" else "cv"
        targets[[col]][k] <- targets[[col]][k] * m
      }
      n_alloc_joint(z$frame, z$measures, targets = targets,
                    objective = objective, budget = 400,
                    objective_rule = rule)
    }
    expect_equal(grid$.target_base[1:2], c(0.06, 3))
    expect_equal(grid$.target[1:2], c(0.057, 3.6))
    expect_true(all(is.na(grid[3, c(".target_base", ".target")])))
    expect_equal(grid$objective_value_int,
                 vapply(seq_len(3), function(i) {
                   refit(newdata$target_constraint[i],
                         newdata$target_mult[i])$operational$objective_value
                 }, 1), tolerance = 1e-12)
    checked <- c(checked, paste(rule, expect_named_rows(fit, grid, refit)))
  }
  expect_length(checked, 6L)
  ## A row without a name is the target_mult row of a grid without names.
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  all_rows <- predict(fit, data.frame(target_mult = 0.98))
  mixed <- predict(fit, data.frame(target_constraint = NA_character_,
                                   target_mult = 0.98))
  shared <- intersect(names(all_rows), names(mixed))
  expect_equal(mixed[shared], all_rows[shared])
})

test_that("a named row on a cluster and a free-take frame matches its refit", {
  cluster <- data.frame(stratum = c("A", "B"), N = c(20000, 8000),
                        N_1 = c(400, 160), cost_1 = 400, cost_2 = 40,
                        take_2 = c(12, 8), region = c("N", "S"))
  coverage <- data.frame(stratum = c("A", "B"), name = "coverage",
                         p = c(0.5, 0.4), icc_1 = 0.05)
  goal <- data.frame(name = "coverage", domain = c(".overall", "region"),
                     cv = c(0.10, 0.2), constraint = c("national", NA))
  goal <- rbind(goal[1, ], data.frame(name = "coverage", domain = "region",
                                      cv = 0.2, constraint = c("north", "south")))
  goal$level <- c(NA, "N", "S")
  checked <- character(0)
  for (free in c(FALSE, TRUE)) {
    frame <- if (free) cluster[names(cluster) != "take_2"] else cluster
    args <- list(frame, coverage, targets = goal)
    if (free) args <- c(args, list(take_2_min = 2, take_2_max = 20))
    fit <- do.call(n_alloc_joint, args)
    grid <- predict(fit, data.frame(target_constraint = c("north", "national"),
                                    target_mult = c(0.8, 0.9)))
    refit <- function(id, m) {
      targets <- goal
      targets$cv[targets$constraint == id] <- targets$cv[targets$constraint == id] * m
      do.call(n_alloc_joint, c(list(frame, coverage, targets = targets),
                               if (free) list(take_2_min = 2, take_2_max = 20)))
    }
    checked <- c(checked, paste(free, expect_named_rows(fit, grid, refit)))
  }
  expect_length(checked, 4L)
})

test_that("a named grid takes a budget column too", {
  z <- named_inputs()
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  grid <- predict(fit, data.frame(budget = c(400, 500),
                                  target_constraint = "coverage@.overall:cv",
                                  target_mult = 0.9))
  for (i in 1:2) {
    targets <- fit$params$targets
    targets$cv[1] <- 0.054
    direct <- n_alloc_joint(z$frame, z$measures, targets = targets,
                            objective = "income", budget = grid$budget[i])
    expect_equal(grid$objective_value[i], direct$objective_value,
                 tolerance = 1e-10)
    expect_equal(grid$.achieved_int[i],
                 direct$operational$constraints$.achieved[1],
                 tolerance = 1e-10)
  }
})

test_that("every joint grid reports its search status", {
  ## A node-limited search keeps its design and reads "limit".
  set.seed(3)
  H <- 10
  frame <- data.frame(stratum = sprintf("s%02d", 1:H),
                      N = round(stats::runif(H, 2e4, 8e4)),
                      unit_cost = round(stats::runif(H, 1, 5), 1))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = H), mean = 10,
                         sd = stats::runif(2 * H, 1, 10))
  fit <- n_alloc_joint(frame, measures, objective = c("y", "z"),
                       budget = 8000, objective_rule = "worst", fpc = "unit")
  grid <- predict(fit, data.frame(budget = 8000))
  expect_identical(grid$.status, fit$operational$status)
  expect_identical(grid$.status, "limit")
  expect_true(grid$.feasible)
  expect_equal(grid$objective_value_int, fit$operational$objective_value)

  ## A refit that proves nothing reads "unknown" with .feasible NA, and a
  ## refused one "infeasible" with .feasible FALSE.
  z <- named_inputs()
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  real <- n_alloc_joint.default
  local_mocked_bindings(n_alloc_joint.default = function(...) {
    args <- list(...)
    if (isTRUE(args$budget == 450)) {
      stop(structure(list(message = "no whole design was found", call = NULL),
                     class = c("svyplan_plan_unknown", "error", "condition")))
    }
    real(...)
  })
  grid <- suppressWarnings(predict(fit, data.frame(budget = c(450, 250))))
  expect_identical(grid$.status, c("unknown", "infeasible"))
  expect_identical(grid$.feasible, c(NA, FALSE))
  expect_warning(predict(fit, data.frame(budget = 450)), "unknown")
})

test_that("a named grid refuses what it cannot scale", {
  z <- named_inputs()
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  expect_error(predict(fit, data.frame(target_constraint = "missing",
                                       target_mult = 0.9)),
               "names no target of the fit")
  expect_error(predict(fit, data.frame(target_constraint =
                                         "coverage@.overall:cv")),
               "needs a 'target_mult'")
  plain <- n_alloc_joint(z$frame, z$measures, objective = "income",
                         budget = 400)
  expect_error(predict(plain, data.frame(target_constraint =
                                           "coverage@.overall:cv",
                                         target_mult = 0.9)),
               "unknown parameter")
})

test_that("a register grid carries a plan only to the same scenario", {
  ## A row stopped at the node limit takes the best plan of a lower budget
  ## in the same scenario. The id "NA" and a missing selector are different
  ## scenarios, as are multipliers equal to 7 digits.
  z <- .psu_fixture()
  z$frame$cost_1 <- 400
  z$frame$cost_2 <- 20
  z$targets$constraint <- c("NA", "income")
  z$targets$cv <- c(0.045, 0.019)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets,
                       objective = "literacy", budget = 40000, psu = z$psu)
  real <- n_alloc_joint.default
  local_mocked_bindings(n_alloc_joint.default = function(...) {
    if (isTRUE(list(...)$budget == 50000)) {
      stop(structure(list(message = "no plan found within the node limit",
                          call = NULL),
                     class = c("svyplan_plan_unknown", "error", "condition")))
    }
    real(...)
  })
  apart <- suppressWarnings(predict(fit, data.frame(
    budget = c(40000, 50000), target_constraint = c("NA", NA),
    target_mult = 0.9
  )))
  expect_identical(apart$.status, c("optimal", "unknown"))
  expect_identical(apart$.feasible, c(TRUE, NA))
  close <- suppressWarnings(predict(fit, data.frame(
    budget = c(40000, 50000), target_constraint = "NA",
    target_mult = c(0.9, 0.9 + 1e-9)
  )))
  expect_identical(close$.status, c("optimal", "unknown"))
  ## The same scenario at a higher budget takes the plan, and the named
  ## target's whole precision with it.
  same <- suppressWarnings(predict(fit, data.frame(
    budget = c(40000, 50000), target_constraint = "NA", target_mult = 0.9
  )))
  expect_identical(same$.status, c("optimal", "limit"))
  expect_identical(same$cost_int[2], same$cost_int[1])
  expect_identical(same$.achieved_int[2], same$.achieved_int[1])
  expect_false(is.na(same$.achieved_int[2]))
})

test_that("a named rmoe target on fixed and free takes matches its refit", {
  cluster <- data.frame(stratum = c("A", "B"), N = c(20000, 8000),
                        N_1 = c(400, 160), cost_1 = 400, cost_2 = 40,
                        take_2 = c(12, 8))
  measures <- data.frame(stratum = rep(c("A", "B"), 2),
                         name = rep(c("coverage", "income"), each = 2),
                         p = c(0.5, 0.4, NA, NA), mean = c(NA, NA, 50, 60),
                         sd = c(NA, NA, 20, 30), icc_1 = 0.05)
  targets <- data.frame(name = "coverage", rmoe = 0.2)
  checked <- character(0)
  for (free in c(FALSE, TRUE)) {
    frame <- if (free) cluster[names(cluster) != "take_2"] else cluster
    takes <- if (free) list(take_2_min = 2, take_2_max = 20)
    for (rule in c("sum", "worst")) {
      fit_with <- function(t) {
        do.call(n_alloc_joint, c(list(frame, measures, targets = t,
                                      objective = c("coverage", "income"),
                                      budget = 60000, objective_rule = rule),
                                 takes))
      }
      fit <- fit_with(targets)
      grid <- predict(fit, data.frame(target_constraint = "coverage@.overall:rmoe",
                                      target_mult = 0.8))
      direct <- fit_with(transform(targets, rmoe = 0.2 * 0.8))
      label <- paste(if (free) "free" else "fixed", rule)
      expect_equal(grid$n, direct$n, tolerance = 1e-8, label = label)
      expect_equal(grid$cost_int, direct$operational$cost, label = label)
      expect_equal(grid$objective_value_int,
                   direct$operational$objective_value, tolerance = 1e-12,
                   label = label)
      expect_identical(grid$.status, direct$operational$status, label = label)
      expect_equal(grid$.target, 0.16, tolerance = 1e-15, label = label)
      expect_equal(grid$.achieved, direct$constraints$.achieved,
                   tolerance = 1e-10, label = label)
      expect_equal(grid$.achieved_int,
                   direct$operational$constraints$.achieved,
                   tolerance = 1e-10, label = label)
      checked <- c(checked, label)
    }
  }
  expect_length(checked, 4L)
})

test_that("a fit saved without objective_rule predicts under the sum rule", {
  z <- named_inputs()
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  legacy <- fit
  legacy$params$objective_rule <- NULL
  newdata <- data.frame(target_constraint = "coverage@.overall:cv",
                        target_mult = 0.95, budget = c(400, 450))
  expect_equal(predict(legacy, newdata), predict(fit, newdata))
})

test_that("a grid reads a bare NA selector as every target", {
  z <- named_inputs()
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  bare <- predict(fit, data.frame(target_constraint = NA, target_mult = 0.98))
  typed <- predict(fit, data.frame(target_constraint = NA_character_,
                                   target_mult = 0.98))
  expect_equal(bare, typed)
})

test_that("a grid refuses a column given twice and a take beside its bounds", {
  z <- named_inputs()
  fit <- n_alloc_joint(z$frame, z$measures, targets = z$targets,
                       objective = "income", budget = 400)
  twice <- data.frame(target_mult = 1, target_mult = 0.5, check.names = FALSE)
  expect_error(predict(fit, twice), "more than once")
  expect_error(predict(fit, 1:3), "must be a data frame")
  expect_error(predict(fit, matrix(1, 2, 2)), "must be a data frame")
  expect_error(predict(n_prop(p = 0.3, moe = 0.05),
                       data.frame(deff = 1, deff = 2, check.names = FALSE)),
               "more than once")
  cluster <- data.frame(stratum = c("A", "B"), N = c(20000, 8000),
                        N_1 = c(400, 160), cost_1 = 400, cost_2 = 40)
  coverage <- data.frame(stratum = c("A", "B"), name = "coverage",
                         p = c(0.5, 0.4), icc_1 = 0.05)
  free <- n_alloc_joint(cluster, coverage,
                        targets = data.frame(name = "coverage", cv = 0.10),
                        take_2_min = 2, take_2_max = 20)
  expect_error(predict(free, data.frame(take_2 = 12, take_2_max = 2)),
               "beside 'take_2_min' or 'take_2_max'")
})

test_that("a fractional floor's frontier starts at the cheapest whole design", {
  ## Oracle: every whole design of three strata. The floor of 8.2 in N is
  ## 9 units in whole numbers, which cost 13 at the least.
  frame <- data.frame(stratum = letters[1:3], N = c(9, 8, 9),
                      region = c("N", "N", "S"), unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 3), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1))
  fit <- n_alloc_joint(frame, measures, objective = c("y", "z"), budget = 21,
                       min_n_domain = data.frame(domain = "region",
                                                 level = "N", min_n = 8.2),
                       fpc = "none")
  states <- as.matrix(expand.grid(lapply(frame$N, seq_len)))
  cheapest <- min(states[rowSums(states[, 1:2]) >= 9, ] %*% frame$unit_cost)
  expect_identical(svyplan:::.budget_grid(fit, 2)[1], cheapest)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, npoints = 2))
})
