test_that("fixed-budget objective mode has stable public result semantics", {
  z <- .bethel_fixture()
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = "income", budget = 5000
  )

  expect_s3_class(fit, "svyplan_n")
  expect_identical(fit$method, "bethel")
  expect_identical(fit$params$mode, "budget_objective")
  expect_equal(fit$params$budget, 5000)
  expect_s3_class(fit$objective, "data.frame")
  expect_true(all(c("component", "name", "domain", "level", "priority",
                    ".relvar", ".cv", ".contribution", ".share") %in%
                    names(fit$objective)))
  expect_equal(fit$objective_value, sum(fit$objective$.contribution))
  expect_equal(fit$objective$.cv, sqrt(fit$objective$.relvar))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
  expect_lte(fit$operational$cost, 5000)
  expect_equal(fit$operational$budget_residual, 5000 - fit$operational$cost)
  expect_true(fit$optimization$budget_binding)
})

test_that("the mode table rejects incomplete objective specifications", {
  z <- .bethel_fixture()
  expect_error(
    n_alloc(z$frame, measures = z$measures, objective = "income"),
    "'objective' requires 'budget'"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            objective = "income"),
    "'objective' requires 'budget'"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            budget = 5000),
    "'budget' requires 'objective'"
  )
  expect_error(
    n_alloc(z$frame, objective = "income", budget = 5000),
    "requires 'measures'"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, objective = "income",
            budget = 5000, n = 100),
    "cannot be combined"
  )
})

test_that("an objective needs no targets", {
  z <- .bethel_fixture()
  fit <- n_alloc(
    z$frame, measures = z$measures, objective = "vaccination",
    budget = 4000, min_n_stratum = 5
  )
  expect_identical(nrow(fit$constraints), 0L)
  expect_true(fit$operational$all_pass)
  expect_lte(fit$operational$cost, 4000)
  expect_equal(fit$params$achieved$cost, 4000)
})

test_that("objective specifications are validated", {
  z <- .bethel_fixture()
  bad <- function(objective, ...) {
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            objective = objective, budget = 5000, ...)
  }
  expect_error(
    bad(data.frame(name = "income", weight = 2)),
    "'priority'"
  )
  expect_error(
    bad(data.frame(name = "income", cv = 0.05)),
    "relative-variance units"
  )
  expect_error(
    bad(data.frame(name = "income", priority = -1)),
    "non-negative"
  )
  expect_error(
    bad(data.frame(name = c("income", "vaccination"), priority = c(0, 0))),
    "at least one objective 'priority' must be positive"
  )
  expect_error(
    bad(data.frame(name = c("income", "income"), priority = c(1, 1))),
    "duplicate indicator-domain objective components"
  )
  expect_error(
    bad(data.frame(name = "income", domain = "region", level = "Nowhere",
                   priority = 1)),
    "objective component domain 'region=Nowhere' is empty"
  )
  expect_error(
    bad(data.frame(name = "income", domain = "region", priority = 1)),
    "must be supplied together"
  )
  expect_error(
    bad(data.frame(name = "absent", priority = 1)),
    "missing measures for indicator 'absent'"
  )
})

test_that("an objective on a negligible total is rejected like a CV target", {
  z <- .bethel_fixture()
  # N = c(1000, 2000, 1500, 1000), so these means total exactly zero
  z$measures$mean[z$measures$name == "income"] <- c(2, -1, 0, 0)
  expect_error(
    n_alloc(z$frame, measures = z$measures,
            objective = "income", budget = 5000),
    "objective is undefined for component 'income@.overall'"
  )
})

test_that("priorities are scale free and zero priorities are inert", {
  z <- .bethel_fixture()
  base <- data.frame(
    name = c("vaccination", "income"), domain = ".overall",
    level = NA_character_, priority = c(2, 1)
  )
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = base, budget = 6000)

  scaled <- base
  scaled$priority <- scaled$priority * 37.5
  rescaled <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                      objective = scaled, budget = 6000)
  expect_equal(rescaled$detail$n, fit$detail$n, tolerance = 1e-8)
  expect_equal(rescaled$objective_value, 37.5 * fit$objective_value,
               tolerance = 1e-8)

  inert <- rbind(base, data.frame(
    name = "income", domain = "residence", level = "Urban", priority = 0
  ))
  ignored <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     objective = inert, budget = 6000)
  expect_equal(ignored$detail$n, fit$detail$n, tolerance = 1e-8)
  expect_equal(ignored$objective_value, fit$objective_value, tolerance = 1e-10)
})

test_that("the character shorthand matches the equal-priority long form", {
  z <- .bethel_fixture()
  short <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   objective = c("vaccination", "income"), budget = 6000)
  long <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = data.frame(
      name = c("vaccination", "income"), domain = ".overall",
      level = NA_character_, priority = c(1, 1)
    ),
    budget = 6000
  )
  expect_equal(short$detail$n, long$detail$n)
  expect_equal(short$objective_value, long$objective_value)
})

test_that("the optimum is monotone in the budget and in target tightness", {
  z <- .bethel_fixture()
  objective <- "income"
  budgets <- c(3000, 4000, 5000, 6000)
  values <- vapply(budgets, function(b) {
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            objective = objective, budget = b)$objective_value
  }, numeric(1L))
  expect_true(all(diff(values) <= 1e-12))

  loose <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   objective = objective, budget = 5000)
  tight <- z$targets
  tight$cv[1L] <- tight$cv[1L] / 2
  tightened <- n_alloc(z$frame, measures = z$measures, targets = tight,
                       objective = objective, budget = 5000)
  expect_gte(tightened$objective_value, loose$objective_value - 1e-12)
})

test_that("non-binding targets leave the budget-only solution unchanged", {
  z <- .bethel_fixture()
  bare <- n_alloc(z$frame, measures = z$measures, objective = "income",
                  budget = 4000, min_n_stratum = 5)
  slack <- data.frame(
    name = "vaccination", domain = ".overall", level = NA_character_, cv = 5
  )
  padded <- n_alloc(z$frame, measures = z$measures, targets = slack,
                    objective = "income", budget = 4000, min_n_stratum = 5)
  expect_equal(padded$detail$n, bare$detail$n, tolerance = 1e-6)
  expect_equal(padded$objective_value, bare$objective_value, tolerance = 1e-10)
  expect_false(padded$constraints$.binding)
})

test_that("budget mode is the dual of minimum-cost mode", {
  z <- .bethel_fixture()
  budgeted <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                      objective = "income", budget = 6000)
  pinned <- rbind(z$targets, data.frame(
    name = "income", domain = ".overall", level = NA_character_,
    cv = sqrt(budgeted$objective_value), moe = NA_real_
  ))
  cheapest <- n_alloc(z$frame, measures = z$measures, targets = pinned)

  expect_equal(cheapest$params$achieved$cost, 6000, tolerance = 1e-6)
  expect_equal(cheapest$detail$n, budgeted$detail$n, tolerance = 1e-5)
})

test_that("the continuous solution spends the budget unless bounds intervene", {
  z <- .bethel_fixture()
  spent <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   objective = "income", budget = 6000)
  expect_equal(spent$params$achieved$cost, 6000, tolerance = 1e-6)
  expect_true(spent$optimization$budget_binding)
  expect_true(is.finite(spent$optimization$budget_sensitivity))
  expect_lt(spent$optimization$budget_sensitivity, 0)

  census_cost <- sum(z$frame$N * z$frame$unit_cost)
  capped <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                    objective = "income", budget = census_cost * 2)
  expect_false(capped$optimization$budget_binding)
  expect_equal(capped$detail$n, z$frame$N)
  expect_equal(capped$params$achieved$cost, census_cost)
  expect_gt(capped$optimization$budget_residual, 0)
})

test_that("the three infeasibility cases are distinguished", {
  z <- .bethel_fixture()
  unattainable <- .bethel_multistage_fixture(2L)
  unattainable$targets$cv[1L] <- 1e-6
  expect_error(
    n_alloc(unattainable$frame, measures = unattainable$measures,
            targets = unattainable$targets, objective = "vaccination",
            budget = 1e7),
    "precision targets are infeasible"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            objective = "income", budget = 100),
    "cannot fund the precision targets.*short by"
  )
  # Continuous-feasible but integer-infeasible: the fractional lower bounds fit
  # the budget while their whole-unit counterparts do not.
  expect_error(
    n_alloc(z$frame, measures = z$measures, objective = "income",
            budget = 473, min_n_stratum = 100.5),
    "no integer allocation fits the budget"
  )
})

test_that("the operational design is feasible, affordable, and whole", {
  z <- .bethel_fixture()
  for (budget in c(3500, 4200, 5000, 6000)) {
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   objective = "income", budget = budget)
    expect_lte(fit$operational$cost, budget)
    expect_true(all(fit$operational$constraints$.pass))
    expect_true(all(fit$detail$n_int == floor(fit$detail$n_int)))
    expect_true(all(fit$detail$n_int >= ceiling(fit$detail$.lower - 1e-9)))
    expect_true(all(fit$detail$n_int <= floor(fit$detail$.upper + 1e-9)))
    expect_gte(fit$operational$objective_value, fit$objective_value - 1e-12)
  }
})

test_that("fixed-take multistage budget objectives keep the take identities", {
  for (stages in c(2L, 3L)) {
    z <- .bethel_multistage_fixture(stages)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   objective = "vaccination", budget = 400000)
    take <- if (stages == 2L) fit$detail$n_per_psu else
      fit$detail$n_per_psu * fit$detail$n_per_ssu
    expect_equal(fit$detail$n, fit$detail$n_psu * take)
    expect_equal(fit$detail$n_int, fit$detail$n_psu_int * take)
    expect_equal(fit$detail$n_psu_int, floor(fit$detail$n_psu_int))
    expect_lte(fit$operational$cost, 400000)
    expect_true(all(fit$operational$constraints$.pass))
  }
})

test_that("budget-objective results round trip through prec_alloc", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 6000)
  assessed <- prec_alloc(fit)
  expect_equal(assessed$objective_value, fit$objective_value)
  expect_equal(assessed$objective, fit$objective)
  expect_equal(assessed$params$budget, 6000)

  again <- n_alloc(assessed)
  expect_equal(again$detail$n, fit$detail$n, tolerance = 1e-6)
  expect_equal(again$objective_value, fit$objective_value, tolerance = 1e-9)

  bare <- n_alloc(z$frame, measures = z$measures, objective = "vaccination",
                  budget = 4000, min_n_stratum = 5)
  expect_equal(n_alloc(prec_alloc(bare))$detail$n, bare$detail$n,
               tolerance = 1e-6)
})

test_that("direct joint assessment reports objective components", {
  z <- .bethel_fixture()
  assessed <- prec_alloc(
    z$frame, n = c(400, 700, 600, 400), measures = z$measures,
    targets = z$targets, objective = "income", budget = 3000
  )
  expect_equal(assessed$objective_value, sum(assessed$objective$.contribution))
  expect_equal(assessed$params$budget_residual,
               3000 - assessed$params$achieved$cost)
  expect_error(
    prec_alloc(z$frame, n = c(400, 700, 600, 400), measures = z$measures,
               targets = z$targets, budget = 3000),
    "'budget' requires 'objective'"
  )
})

test_that("predict returns the budget frontier and flags infeasible points", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 5000)
  frontier <- predict(fit, data.frame(budget = c(4000, 5000, 6000)))

  expect_equal(frontier$budget, c(4000, 5000, 6000))
  expect_true(all(frontier$.feasible))
  expect_true(all(diff(frontier$objective_value) <= 1e-12))
  expect_equal(frontier$cv, sqrt(frontier$objective_value))
  expect_true(all(frontier$cost_int <= frontier$budget))
  expect_equal(frontier$objective_value[2L], fit$objective_value)

  expect_warning(
    partial <- predict(fit, data.frame(budget = c(100, 5000))),
    "infeasible"
  )
  expect_false(partial$.feasible[1L])
  expect_true(is.na(partial$objective_value[1L]))
  expect_true(partial$.feasible[2L])

  expect_error(
    predict(fit, data.frame(cv = 0.1)),
    "unknown parameter"
  )
  targets_only <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  expect_error(
    predict(targets_only, data.frame(budget = 5000)),
    "modify 'targets'"
  )
})

test_that("print names the question each mode answers", {
  z <- .bethel_fixture()
  budgeted <- capture.output(print(n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = "income", budget = 5000
  )))
  expect_true(any(grepl("best design affordable within a budget", budgeted)))
  expect_true(any(grepl("weighted relative variance", budgeted)))

  cheapest <- capture.output(print(
    n_alloc(z$frame, measures = z$measures, targets = z$targets)
  ))
  expect_true(any(grepl("cheapest design meeting every precision target",
                        cheapest)))
  expect_false(any(grepl("weighted relative variance", cheapest)))
})

test_that("an unaffordable budget names a budget that would work", {
  z <- .bethel_fixture()
  base <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  cost0 <- base$params$achieved$cost

  err <- tryCatch(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            objective = "income", budget = cost0),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("Raise 'budget' to at least", err, fixed = TRUE))

  suggested <- as.numeric(sub(
    ".*Raise 'budget' to at least ([0-9.]+).*", "\\1",
    gsub("\n", " ", err)
  ))
  expect_true(is.finite(suggested))
  expect_gt(suggested, cost0)

  # the number the message names must itself be affordable
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = suggested)
  expect_lte(fit$operational$cost, suggested)
})
