test_that("two-stage generalized allocation exposes both public units", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc_joint(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  d <- fit$detail

  expect_s3_class(fit, "svyplan_alloc_joint")
  expect_s3_class(fit, "svyplan_alloc")
  expect_identical(fit$method, "bethel")
  expect_identical(fit$params$stages, 2L)
  expect_true(all(c("n_1", "n_1_int", "take_2") %in% names(d)))
  expect_equal(d$n, d$n_1 * d$take_2)
  expect_equal(d$n_int, d$n_1_int * d$take_2)
  expect_type(d$n_int, "double")
  expect_type(d$n_1_int, "double")
  expect_true(all(d$n_1_int == floor(d$n_1_int)))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
  expect_equal(
    fit$params$achieved$cost,
    sum(d$n_1 * (d$cost_1 + d$cost_2 * d$take_2))
  )
})

test_that("three-stage generalized allocation preserves fixed-take identities", {
  z <- .bethel_multistage_fixture(3L)
  fit <- n_alloc_joint(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  d <- fit$detail

  expect_identical(fit$params$stages, 3L)
  expect_true(all(c("N_2", "take_3", "cost_3") %in% names(d)))
  expect_equal(d$n, d$n_1 * d$take_2 * d$take_3)
  expect_equal(d$n_int, d$n_1_int * d$take_2 * d$take_3)
  expect_type(d$n_int, "double")
  expect_type(d$n_1_int, "double")
  expect_equal(fit$n, sum(d$n))
  expect_equal(fit$operational$n, sum(d$n_int))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
})

test_that("the take sweep reproduces a refit at every grid point", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets)
  takes <- c(6, 8, 10, 14, 20)

  sweep <- predict(fit, data.frame(take_2 = takes))

  manual <- do.call(rbind, lapply(takes, function(m) {
    frame <- z$frame
    frame$take_2 <- m
    refit <- n_alloc_joint(frame, measures = z$measures, targets = z$targets)
    data.frame(
      n_1 = sum(refit$detail$n_1),
      n = refit$n,
      cost = refit$params$achieved$cost,
      n_1_int = sum(refit$detail$n_1_int),
      n_int = as.numeric(refit$operational$n),
      cost_int = refit$operational$cost
    )
  }))

  expect_equal(sweep$take_2, takes)
  expect_true(all(sweep$.feasible))
  expect_equal(sweep[names(manual)], manual, ignore_attr = TRUE)
  # The ultimate sample rises with the take while the PSU count falls, which
  # is what makes the cost curve over the take U-shaped.
  expect_true(all(diff(sweep$n) > 0))
  expect_true(all(diff(sweep$n_1) < 0))
})

test_that("the take sweep varies every stratum and holds the rest of the plan", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets,
                 min_n_stratum = 20)

  sweep <- predict(fit, data.frame(take_2 = 9))
  refit <- n_alloc_joint(
    transform(z$frame, take_2 = 9), measures = z$measures,
    targets = z$targets, min_n_stratum = 20
  )

  expect_equal(sweep$n, refit$n)
  expect_equal(sweep$n_1, sum(refit$detail$n_1))
  # A scalar take applies to every stratum, so the per-stratum takes the fit
  # was built from (8, 10, 12, 7) are replaced rather than scaled.
  expect_true(all(refit$detail$take_2 == 9))
})

test_that("a three-stage fit sweeps either take", {
  z <- .bethel_multistage_fixture(3L)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets)

  ssu <- predict(fit, data.frame(take_3 = c(2, 3, 5)))
  expect_equal(ssu$take_3, c(2, 3, 5))
  expect_true(all(ssu$.feasible))
  expect_true(all(diff(ssu$n_1) < 0))

  crossed <- predict(fit, expand.grid(take_2 = c(4, 8), take_3 = c(2, 4)))
  expect_equal(nrow(crossed), 4L)
  expect_true(all(crossed$.feasible))
  expect_equal(
    crossed$n[crossed$take_2 == 4 & crossed$take_3 == 2],
    predict(fit, data.frame(take_2 = 4, take_3 = 2))$n
  )
})

test_that("an infeasible take gives an NA row and keeps the rest of the grid", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets)

  expect_warning(
    sweep <- predict(fit, data.frame(take_2 = c(10, 100000))),
    "infeasible"
  )
  expect_true(sweep$.feasible[1L])
  expect_false(sweep$.feasible[2L])
  expect_false(is.na(sweep$n[1L]))
  expect_true(all(is.na(sweep[2L, c("n_1", "n", "cost", "n_int", "cost_int")])))
  expect_type(sweep$.feasible, "logical")
})

test_that("the take sweep refuses grids it cannot honour", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets)

  # A fractional take is a malformed grid, not an infeasible design point.
  expect_error(
    predict(fit, data.frame(take_2 = c(6, 7.5))),
    "must contain positive whole numbers"
  )
  expect_error(
    predict(fit, data.frame(take_2 = 0)),
    "must contain positive whole numbers"
  )
  # A two-stage fit has no third stage to sweep.
  expect_error(
    predict(fit, data.frame(take_3 = 3)),
    "unknown parameter"
  )
  # Targets stay out of reach: the answer is to edit them and rerun.
  expect_error(
    predict(fit, data.frame(cv = 0.05)),
    "unknown parameter"
  )
  # A one-stage joint fit has no take to vary at all.
  z1 <- .bethel_fixture()
  flat <- n_alloc_joint(z1$frame, measures = z1$measures, targets = z1$targets)
  expect_error(
    predict(flat, data.frame(take_2 = 10)),
    "unknown parameter"
  )
})

test_that("a budget-objective fit sweeps takes alongside the budget", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 60000)

  crossed <- predict(fit, expand.grid(take_2 = c(8, 14),
                                      budget = c(50000, 60000)))
  expect_equal(nrow(crossed), 4L)
  expect_true(all(crossed$.feasible))
  expect_true(all(crossed$cost_int <= crossed$budget))
  expect_equal(crossed$cv, sqrt(crossed$objective_value))
  # A budget-mode grid reports the PSU count a multistage design turns on.
  expect_true(all(c("n_1", "n_1_int") %in% names(crossed)))

  # Takes alone hold the budget at the value the fit was built with.
  takes_only <- predict(fit, data.frame(take_2 = c(8, 14)))
  expect_equal(takes_only$n, crossed$n[crossed$budget == 60000])
})
