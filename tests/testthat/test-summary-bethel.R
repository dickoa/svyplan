test_that("minimum-cost Bethel summary separates field and continuous designs", {
  z <- .bethel_fixture()
  x <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  s <- summary(x)

  expect_s3_class(s, "summary.svyplan_bethel")
  expect_identical(
    names(s),
    c("kind", "question", "mode", "stages", "status", "overall",
      "continuous", "allocation", "constraints",
      "operational_constraints", "objective", "operational_objective",
      "bounds", "optimization", "assumptions")
  )
  expect_identical(s$kind, "design")
  expect_identical(s$mode, "targets")
  expect_identical(s$status, "optimal")
  expect_match(s$question, "cheapest joint allocation")

  expect_equal(s$overall$n, x$operational$n)
  expect_equal(s$overall$cost, x$operational$cost)
  expect_true(s$overall$all_pass)
  expect_equal(s$continuous$n, x$n)
  expect_equal(s$continuous$cost, x$params$achieved$cost)
  expect_equal(s$allocation$n_field, x$detail$n_int)
  expect_equal(s$allocation$n_continuous, x$detail$n)
  expect_equal(sum(s$allocation$cost), x$operational$cost)
  expect_equal(s$constraints, x$constraints)
  expect_equal(s$operational_constraints, x$operational$constraints)
  expect_true(any(is.finite(s$constraints$.sensitivity)))
  expect_true(all(s$optimization$diagnostics$pass))
  expect_identical(s$optimization$active_precision, x$binding)
  expect_identical(nrow(s$bounds), nrow(z$frame))
  expect_identical(nrow(s$assumptions$model), 8L)
  expect_error(summary(x, digits = 3), "unused argument")
})

test_that("budget-objective summary retains both additive decompositions", {
  z <- .bethel_fixture()
  x <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = c("vaccination", "income"), budget = 5000
  )
  s <- summary(x)

  expect_identical(s$mode, "budget_objective")
  expect_equal(s$overall$objective_value, x$operational$objective_value)
  expect_equal(s$continuous$objective_value, x$objective_value)
  expect_equal(s$objective, x$objective)
  expect_equal(s$operational_objective, x$operational$objective)
  expect_equal(sum(s$objective$.contribution), x$objective_value)
  expect_equal(sum(s$objective$.share), 1)
  expect_equal(
    sum(s$operational_objective$.contribution),
    x$operational$objective_value
  )
  expect_equal(s$overall$budget_residual, x$operational$budget_residual)
  expect_true(s$optimization$budget_binding)
  expect_equal(
    s$optimization$budget_sensitivity,
    x$optimization$budget_sensitivity
  )
  expect_identical(s$optimization$active_precision, "objective")
})

test_that("Bethel precision summary assesses the supplied allocation exactly", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  supplied <- fit$detail$n_int
  supplied[1L] <- 0.5
  x <- prec_alloc(fit, n = supplied)
  s <- summary(x)

  expect_identical(s$kind, "assessment")
  expect_identical(s$status, "assessed")
  expect_null(s$continuous)
  expect_null(s$operational_constraints)
  expect_null(s$optimization)
  expect_equal(s$overall$n, sum(supplied))
  expect_equal(s$overall$cost, x$params$achieved$cost)
  expect_equal(s$allocation$n_supplied, supplied)
  expect_equal(s$constraints, x$detail)
  expect_true(s$overall$n_violated > 0L)
  expect_equal(s$overall$n_bound_violations, 1L)
  expect_identical(s$bounds$status[1L], "lower violation")
  expect_false(s$bounds$pass[1L])

  objective_only <- prec_alloc(
    z$frame, n = fit$detail$n_int, measures = z$measures,
    targets = z$targets, objective = "income"
  )
  objective_summary <- summary(objective_only)
  expect_identical(objective_summary$mode, "objective")
  expect_null(objective_summary$overall$budget)
  expect_match(objective_summary$question, "weighted objective")
  expect_false(any(grepl("budget:", capture.output(print(objective_only)),
                         fixed = TRUE)))
})

test_that("multistage Bethel summary exposes PSU decisions and fixed takes", {
  z <- .bethel_multistage_fixture(3L)
  x <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  s <- summary(x)

  expect_identical(s$stages, 3L)
  expect_equal(s$allocation$psu_continuous, x$detail$n_psu)
  expect_equal(s$allocation$psu_field, x$detail$n_psu_int)
  expect_equal(s$allocation$take_per_psu, z$frame$n_per_psu)
  expect_equal(s$allocation$take_per_ssu, z$frame$n_per_ssu)
  expect_equal(sum(s$allocation$cost), x$operational$cost)
  expect_equal(s$assumptions$fixed_takes$n_per_psu, z$frame$n_per_psu)
  expect_equal(s$assumptions$fixed_takes$n_per_ssu, z$frame$n_per_ssu)
})

test_that("Bethel summary print is sectioned, bounded, and explicit", {
  z <- .bethel_fixture()
  x <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = "income", budget = 5000
  )
  out <- capture.output(expect_invisible(print(summary(x))))

  expect_match(out[1L], "Generalized Bethel allocation summary")
  expect_true(any(grepl("Operational field design", out, fixed = TRUE)))
  expect_true(any(grepl("Continuous precision constraints", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Operational precision constraints", out,
                        fixed = TRUE)))
  expect_true(any(grepl("contributions are additive", out, fixed = TRUE)))
  expect_true(any(grepl("Allocation bounds", out, fixed = TRUE)))
  expect_true(any(grepl("Optimization diagnostics", out, fixed = TRUE)))
  expect_true(any(grepl("KKT certification", out, fixed = TRUE)))
  expect_true(any(grepl("Resolved target-stratum", out, fixed = TRUE)))
  expect_true(all(nchar(out) <= 80L))
})
