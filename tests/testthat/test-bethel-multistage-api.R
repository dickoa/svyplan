test_that("two-stage generalized allocation exposes both public units", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  d <- fit$detail

  expect_s3_class(fit, "svyplan_n")
  expect_identical(fit$method, "bethel")
  expect_identical(fit$params$stages, 2L)
  expect_true(all(c("n_psu", "n_psu_int", "n_per_psu") %in% names(d)))
  expect_equal(d$n, d$n_psu * d$n_per_psu)
  expect_equal(d$n_int, d$n_psu_int * d$n_per_psu)
  expect_type(d$n_int, "double")
  expect_type(d$n_psu_int, "double")
  expect_true(all(d$n_psu_int == floor(d$n_psu_int)))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
  expect_equal(
    fit$params$achieved$cost,
    sum(d$n_psu * (d$cost_psu + d$cost_ssu * d$n_per_psu))
  )
})

test_that("three-stage generalized allocation preserves fixed-take identities", {
  z <- .bethel_multistage_fixture(3L)
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  d <- fit$detail

  expect_identical(fit$params$stages, 3L)
  expect_true(all(c("N_ssu", "n_per_ssu", "cost_tsu") %in% names(d)))
  expect_equal(d$n, d$n_psu * d$n_per_psu * d$n_per_ssu)
  expect_equal(d$n_int, d$n_psu_int * d$n_per_psu * d$n_per_ssu)
  expect_type(d$n_int, "double")
  expect_type(d$n_psu_int, "double")
  expect_equal(fit$n, sum(d$n))
  expect_equal(fit$operational$n, sum(d$n_int))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
})

test_that("the take sweep reproduces a refit at every grid point", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  takes <- c(6, 8, 10, 14, 20)

  sweep <- predict(fit, data.frame(n_per_psu = takes))

  manual <- do.call(rbind, lapply(takes, function(m) {
    frame <- z$frame
    frame$n_per_psu <- m
    refit <- n_alloc(frame, measures = z$measures, targets = z$targets)
    data.frame(
      n_psu = sum(refit$detail$n_psu),
      n = refit$n,
      cost = refit$params$achieved$cost,
      n_psu_int = sum(refit$detail$n_psu_int),
      n_int = as.numeric(refit$operational$n),
      cost_int = refit$operational$cost
    )
  }))

  expect_equal(sweep$n_per_psu, takes)
  expect_true(all(sweep$.feasible))
  expect_equal(sweep[names(manual)], manual, ignore_attr = TRUE)
  # The ultimate sample rises with the take while the PSU count falls, which
  # is what makes the cost curve over the take U-shaped.
  expect_true(all(diff(sweep$n) > 0))
  expect_true(all(diff(sweep$n_psu) < 0))
})

test_that("the take sweep varies every stratum and holds the rest of the plan", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 min_n_stratum = 20)

  sweep <- predict(fit, data.frame(n_per_psu = 9))
  refit <- n_alloc(
    transform(z$frame, n_per_psu = 9), measures = z$measures,
    targets = z$targets, min_n_stratum = 20
  )

  expect_equal(sweep$n, refit$n)
  expect_equal(sweep$n_psu, sum(refit$detail$n_psu))
  # A scalar take applies to every stratum, so the per-stratum takes the fit
  # was built from (8, 10, 12, 7) are replaced rather than scaled.
  expect_true(all(refit$detail$n_per_psu == 9))
})

test_that("a three-stage fit sweeps either take", {
  z <- .bethel_multistage_fixture(3L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)

  ssu <- predict(fit, data.frame(n_per_ssu = c(2, 3, 5)))
  expect_equal(ssu$n_per_ssu, c(2, 3, 5))
  expect_true(all(ssu$.feasible))
  expect_true(all(diff(ssu$n_psu) < 0))

  crossed <- predict(fit, expand.grid(n_per_psu = c(4, 8), n_per_ssu = c(2, 4)))
  expect_equal(nrow(crossed), 4L)
  expect_true(all(crossed$.feasible))
  expect_equal(
    crossed$n[crossed$n_per_psu == 4 & crossed$n_per_ssu == 2],
    predict(fit, data.frame(n_per_psu = 4, n_per_ssu = 2))$n
  )
})

test_that("an infeasible take gives an NA row and keeps the rest of the grid", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)

  expect_warning(
    sweep <- predict(fit, data.frame(n_per_psu = c(10, 100000))),
    "infeasible"
  )
  expect_true(sweep$.feasible[1L])
  expect_false(sweep$.feasible[2L])
  expect_false(is.na(sweep$n[1L]))
  expect_true(all(is.na(sweep[2L, c("n_psu", "n", "cost", "n_int", "cost_int")])))
  expect_type(sweep$.feasible, "logical")
})

test_that("the take sweep refuses grids it cannot honour", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)

  # A fractional take is a malformed grid, not an infeasible design point.
  expect_error(
    predict(fit, data.frame(n_per_psu = c(6, 7.5))),
    "must contain positive whole numbers"
  )
  expect_error(
    predict(fit, data.frame(n_per_psu = 0)),
    "must contain positive whole numbers"
  )
  # A two-stage fit has no third stage to sweep.
  expect_error(
    predict(fit, data.frame(n_per_ssu = 3)),
    "unknown parameter"
  )
  # Targets stay out of reach: the answer is to edit them and rerun.
  expect_error(
    predict(fit, data.frame(cv = 0.05)),
    "unknown parameter"
  )
  # A one-stage joint fit has no take to vary at all.
  z1 <- .bethel_fixture()
  flat <- n_alloc(z1$frame, measures = z1$measures, targets = z1$targets)
  expect_error(
    predict(flat, data.frame(n_per_psu = 10)),
    "modify 'targets'"
  )
})

test_that("a budget-objective fit sweeps takes alongside the budget", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 objective = "income", budget = 60000)

  crossed <- predict(fit, expand.grid(n_per_psu = c(8, 14),
                                      budget = c(50000, 60000)))
  expect_equal(nrow(crossed), 4L)
  expect_true(all(crossed$.feasible))
  expect_true(all(crossed$cost_int <= crossed$budget))
  expect_equal(crossed$cv, sqrt(crossed$objective_value))
  # A budget-mode grid reports the PSU count a multistage design turns on.
  expect_true(all(c("n_psu", "n_psu_int") %in% names(crossed)))

  # Takes alone hold the budget at the value the fit was built with.
  takes_only <- predict(fit, data.frame(n_per_psu = c(8, 14)))
  expect_equal(takes_only$n, crossed$n[crossed$budget == 60000])
})
