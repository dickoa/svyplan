## The finite population correction a cluster allocation carries
##
## `fpc` picks which correction the variance carries, and the three choices are
## nested rather than arbitrary. Writing f1 for the PSU fraction and f2 for the
## within-PSU fraction, the exact equal-take two-stage without-replacement
## variance corrects each component by its own stage:
##
##   V = (1 - f1) S1^2 / a + (1 - f2) S2^2 / (a m)
##
## "none" applies no correction and matches n_cluster() at any fraction.
## "unit" applies 1 - n/N to the whole inflated variance, and since
## n/N = f1 f2 that factor exceeds both exact ones, so it overstates both
## components. "stage" is the exact form above. The ordering
## none >= unit >= stage therefore holds at every sampling fraction, and the
## three coincide as the fractions vanish.

.fpc_p <- 0.30
.fpc_relvar <- (1 - .fpc_p) / .fpc_p

.fpc_frame <- function(N, N_psu, m, icc = 0.05) {
  data.frame(
    stratum = "S", N = N, mean = .fpc_p, sd = sqrt(.fpc_p * (1 - .fpc_p)),
    icc_psu = icc, n_per_psu = m, N_psu = N_psu,
    cost_psu = 300, cost_ssu = 25
  )
}

## T1. The three corrections are ordered, at every sampling fraction

test_that("none >= unit >= stage at every sampling fraction", {
  m <- 20
  for (N in c(2e6, 2e5, 2e4, 6e3)) {
    frame <- .fpc_frame(N, N / 40, m)
    got <- vapply(
      c("none", "unit", "stage"),
      function(f) n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = f)$detail$n,
      numeric(1)
    )
    expect_gte(got[["none"]], got[["unit"]] - 1e-8)
    expect_gte(got[["unit"]], got[["stage"]] - 1e-8)
  }
})

test_that("the three corrections converge as the fractions vanish", {
  # Both fractions have to be made small, and they are controlled by
  # different quantities. f1 falls with N_psu large relative to the PSUs
  # taken, but f2 = m * r / M is set by the mean PSU size M = N / N_psu and
  # does not move with N at all. A large N over few large PSUs leaves f2
  # wherever it was.
  frame <- .fpc_frame(2e9, 2e5, 20)
  got <- vapply(
    c("none", "unit", "stage"),
    function(f) n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = f)$detail$n,
    numeric(1)
  )
  expect_equal(got[["unit"]], got[["none"]], tolerance = 2e-3)
  expect_equal(got[["stage"]], got[["none"]], tolerance = 2e-3)
})

test_that("growing N alone does not close the gap, because f2 does not move", {
  # The same take in the same size of PSU keeps the same within-PSU
  # fraction however large the population is. This is the reason "stage" is
  # not a small correction to "none" in general.
  wide <- vapply(
    c(2e5, 2e7, 2e9),
    function(N) {
      frame <- .fpc_frame(N, N / 40, 20)
      n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = "stage")$detail$n /
        n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = "none")$detail$n
    },
    numeric(1)
  )
  expect_equal(wide[2], wide[3], tolerance = 1e-3)
  expect_lt(wide[3], 0.9)
})

## T2. "none" is n_cluster()'s model, and matches it away from negligible f too

test_that("fpc = none reproduces n_cluster at an appreciable fraction", {
  icc <- 0.05
  m <- 20
  by_cluster <- n_cluster(
    icc = icc, n_per_psu = m, cv = 0.10, resp_rate = 0.5,
    unit_relvar = .fpc_relvar, stage_cost = c(300, 25)
  )$n[["n_psu"]]
  # N chosen so that n/N is around 12 percent, where "unit" visibly departs.
  frame <- .fpc_frame(5e3, 250, m)
  by_alloc <- n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = "none")$detail$n
  expect_equal(by_alloc / m, by_cluster, tolerance = 1e-6)

  unit <- n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = "unit")$detail$n
  expect_lt(unit, by_alloc - 1)
})

## T3. "stage" is the exact component-wise variance, not an approximation of it

test_that("fpc = stage reproduces the exact two-stage WOR variance", {
  icc <- 0.05
  m <- 10
  r <- 0.6
  N <- 2e4
  N_psu <- 500
  M <- N / N_psu
  frame <- .fpc_frame(N, N_psu, m, icc)

  for (n_psu in c(20, 80, 300, 500)) {
    n <- n_psu * m
    got <- prec_alloc(frame, n = n, resp_rate = r, fpc = "stage")$se^2
    S2 <- .fpc_p * (1 - .fpc_p)
    f1 <- n_psu / N_psu
    f2 <- m * r / M
    exact <- (1 - f1) * S2 * icc / n_psu +
      (1 - f2) * S2 * (1 - icc) / (n_psu * m * r)
    expect_equal(got, exact, tolerance = 1e-10)
  }
})

test_that("fpc = stage reproduces the exact three-stage variance", {
  d1 <- 0.05
  d2 <- 0.10
  m <- 5
  q <- 4
  r <- 0.6
  N <- 2e5
  N_psu <- 2000
  N_ssu <- 20000
  frame <- data.frame(
    stratum = "S", N = N, N_psu = N_psu, N_ssu = N_ssu,
    n_per_psu = m, n_per_ssu = q,
    cost_psu = 300, cost_ssu = 50, cost_tsu = 25
  )
  measures <- data.frame(
    stratum = "S", name = "y", p = .fpc_p, icc_psu = d1, icc_ssu = d2,
    var_ratio_psu = 1, resp_rate = r
  )
  targets <- data.frame(
    name = "y", domain = ".overall", level = NA, cv = 0.08
  )
  fit <- n_alloc(frame, measures = measures, targets = targets, fpc = "stage")

  a <- fit$detail$n / (m * q)
  S2 <- .fpc_p * (1 - .fpc_p)
  k2 <- 1 - d1
  qr <- q * r
  f1 <- a / N_psu
  f2 <- m / (N_ssu / N_psu)
  f3 <- qr / (N / N_ssu)
  exact <- (1 - f1) * S2 * d1 / a +
    (1 - f2) * S2 * k2 * d2 / (a * m) +
    (1 - f3) * S2 * k2 * (1 - d2) / (a * m * qr)
  expect_equal(sqrt(exact) / .fpc_p, 0.08, tolerance = 1e-6)
})

## T4. Zero variance only where the relevant stages are actually enumerated

test_that("a full census is zero variance under unit and stage", {
  N <- 2000
  N_psu <- 100
  M <- N / N_psu
  frame <- .fpc_frame(N, N_psu, M)
  expect_equal(prec_alloc(frame, n = N, fpc = "unit")$se, 0)
  expect_equal(prec_alloc(frame, n = N, fpc = "stage")$se, 0)
  expect_gt(prec_alloc(frame, n = N, fpc = "none")$se, 0)
})

test_that("taking every PSU is not a census while a take is in force", {
  N <- 2000
  N_psu <- 100
  M <- N / N_psu
  frame <- .fpc_frame(N, N_psu, M / 2)
  got <- prec_alloc(frame, n = N_psu * M / 2, fpc = "stage")$se
  # Every PSU is in, so the between-PSU component is gone, but half of each
  # PSU is unobserved and that component is not.
  S2 <- .fpc_p * (1 - .fpc_p)
  exact <- (1 - 0.5) * S2 * (1 - 0.05) / (N_psu * M / 2)
  expect_equal(got, sqrt(exact), tolerance = 1e-10)
  expect_gt(got, 0)
})

## T5. The two allocators agree under every correction

test_that("ordinary and generalized cluster paths agree under every fpc", {
  icc <- 0.05
  m <- 20
  r <- 0.5
  N <- 2e5
  N_psu <- N / 40
  ordinary <- .fpc_frame(N, N_psu, m, icc)
  frame <- data.frame(
    stratum = "S", N = N, N_psu = N_psu, n_per_psu = m,
    cost_psu = 300, cost_ssu = 25
  )
  measures <- data.frame(
    stratum = "S", name = "y", p = .fpc_p, icc_psu = icc,
    var_ratio_psu = 1, resp_rate = r
  )
  targets <- data.frame(
    name = "y", domain = ".overall", level = NA, cv = 0.10
  )
  for (mode in c("none", "unit", "stage")) {
    expect_equal(
      n_alloc(frame, measures = measures, targets = targets, fpc = mode)$detail$n,
      n_alloc(ordinary, cv = 0.10, resp_rate = r, fpc = mode)$detail$n,
      tolerance = 1e-6
    )
  }
})

## T6. Round trips, including through the operational take

test_that("prec_alloc reproduces the target under every fpc", {
  frame <- .fpc_frame(2e5, 5e3, 20)
  for (mode in c("none", "unit", "stage")) {
    fit <- n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = mode)
    back <- prec_alloc(frame, n = fit$detail$n, resp_rate = 0.5, fpc = mode)
    expect_equal(back$cv, 0.10, tolerance = 1e-8)
  }
})

test_that("the operational block reads the take it actually fields", {
  # No 'n_per_psu', so the take is cost-optimal and the whole-unit search
  # picks its own integer take. Both halves of the stage correction have to
  # be rebuilt on that take rather than on the continuous one.
  frame <- data.frame(
    stratum = "S", N = 2e4, mean = .fpc_p, sd = sqrt(.fpc_p * (1 - .fpc_p)),
    icc_psu = 0.05, N_psu = 500, cost_psu = 300, cost_ssu = 25
  )
  fit <- n_alloc(frame, cv = 0.10, resp_rate = 0.6, fpc = "stage")
  a <- fit$detail$n_psu_int
  b <- fit$detail$n_per_psu_int
  S2 <- .fpc_p * (1 - .fpc_p)
  f1 <- a / 500
  f2 <- b * 0.6 / (2e4 / 500)
  exact <- (1 - f1) * S2 * 0.05 / a +
    (1 - f2) * S2 * (1 - 0.05) / (a * b * 0.6)
  expect_equal(fit$operational$cv, sqrt(exact) / .fpc_p, tolerance = 1e-8)
})

## T7. The argument is refused where it has no stages to choose between

test_that("fpc is refused outside cluster allocation", {
  flat <- data.frame(
    stratum = c("a", "b"), N = c(1000, 2000), mean = c(10, 12),
    sd = c(3, 4)
  )
  expect_error(n_alloc(flat, cv = 0.05, fpc = "stage"), "cluster allocation only")
  expect_error(n_alloc(flat, cv = 0.05, fpc = "none"), "cluster allocation only")
  expect_silent(n_alloc(flat, cv = 0.05))

  measures <- data.frame(
    stratum = c("a", "b"), name = "y", p = c(0.3, 0.4)
  )
  targets <- data.frame(
    name = "y", domain = ".overall", level = NA, cv = 0.05
  )
  expect_error(
    n_alloc(flat[, c("stratum", "N")], measures = measures,
            targets = targets, fpc = "stage"),
    "cluster allocation only"
  )
})

test_that("fpc = stage needs the PSU population", {
  frame <- data.frame(
    stratum = "S", N = 2e5, mean = .fpc_p, sd = sqrt(.fpc_p * (1 - .fpc_p)),
    icc_psu = 0.05, n_per_psu = 20, cost_psu = 300, cost_ssu = 25
  )
  expect_error(n_alloc(frame, cv = 0.10, fpc = "stage"), "N_psu")
  expect_silent(n_alloc(frame, cv = 0.10, fpc = "none"))
})

## T8. The chosen correction survives every way of getting back to the fit
##
## A mode that is honoured on the way in and dropped on the way out reports a
## precision the design does not have, and silently, since every value stays
## plausible. These pin the propagation rather than the arithmetic.

test_that("fpc survives the round trip through prec_alloc", {
  frame <- .fpc_frame(2e5, 5e3, 20)
  for (mode in c("none", "unit", "stage")) {
    fit <- n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = mode)
    expect_identical(fit$params$fpc, mode)
    expect_equal(prec_alloc(fit)$cv, fit$cv, tolerance = 1e-8)
  }
})

test_that("fpc survives the round trip back to n_alloc", {
  frame <- .fpc_frame(2e5, 5e3, 20)
  for (mode in c("none", "unit", "stage")) {
    fit <- n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = mode)
    back <- n_alloc(prec_alloc(fit))
    expect_equal(back$detail$n, fit$detail$n, tolerance = 1e-6)
  }
})

test_that("fpc survives the generalized round trip", {
  frame <- data.frame(
    stratum = "S", N = 2e5, N_psu = 5e3, n_per_psu = 20,
    cost_psu = 300, cost_ssu = 25
  )
  measures <- data.frame(
    stratum = "S", name = "y", p = .fpc_p, icc_psu = 0.05,
    var_ratio_psu = 1, resp_rate = 0.5
  )
  targets <- data.frame(
    name = "y", domain = ".overall", level = NA, cv = 0.10
  )
  for (mode in c("none", "unit", "stage")) {
    fit <- n_alloc(frame, measures = measures, targets = targets, fpc = mode)
    expect_identical(fit$params$fpc, mode)
    back <- prec_alloc(fit)
    expect_equal(
      back$detail$.achieved[1L], fit$constraints$.achieved[1L],
      tolerance = 1e-6
    )
  }
})

test_that("predict() sweeps under the fit's own correction", {
  frame <- .fpc_frame(2e5, 5e3, 20)
  for (mode in c("none", "unit", "stage")) {
    fit <- n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = mode)
    swept <- predict(fit, data.frame(cv = c(0.10, 0.08)))
    expect_equal(swept$n[1L], fit$detail$n, tolerance = 1e-6)
  }
  # The modes must stay distinguishable through predict(), or a dropped
  # mode would look like agreement rather than like a bug.
  n_at <- vapply(
    c("none", "unit", "stage"),
    function(m) {
      predict(
        n_alloc(frame, cv = 0.10, resp_rate = 0.5, fpc = m),
        data.frame(cv = 0.08)
      )$n
    },
    numeric(1)
  )
  expect_gt(n_at[["none"]], n_at[["stage"]])
})

## T9. A register is a different design, and says so rather than ignoring it

test_that("fpc is refused with a PSU register", {
  frame <- data.frame(
    stratum = "S", N = 1000, n_per_psu = 10, cost_psu = 300, cost_ssu = 25
  )
  psu <- data.frame(stratum = "S", N = rep(20, 50))
  measures <- data.frame(stratum = "S", name = "y", p = .fpc_p, icc_psu = 0.05)
  targets <- data.frame(
    name = "y", domain = ".overall", level = NA, cv = 0.10
  )
  for (mode in c("none", "stage")) {
    expect_error(
      n_alloc(frame, psu = psu, measures = measures, targets = targets,
              fpc = mode),
      "PSU register"
    )
  }
  expect_s3_class(
    n_alloc(frame, psu = psu, measures = measures, targets = targets),
    "svyplan_n"
  )
})

## T10. The generalized sweep and the frontier floor rebuild the fit too
##
## predict() on a generalized fit re-solves n_alloc() per grid row, and the
## budget-frontier plot re-solves it again to find where the frontier starts.
## Both are rebuilds, so both drop the correction unless handed it.

.fpc_joint <- function() {
  list(
    frame = data.frame(
      stratum = "S", N = 2e5, N_psu = 5e3, n_per_psu = 20,
      cost_psu = 300, cost_ssu = 25
    ),
    measures = data.frame(
      stratum = "S", name = "y", p = .fpc_p, icc_psu = 0.05,
      var_ratio_psu = 1, resp_rate = 0.5
    ),
    targets = data.frame(
      name = "y", domain = ".overall", level = NA, cv = 0.10
    )
  )
}

test_that("predict() on a generalized fit keeps the fit's correction", {
  z <- .fpc_joint()
  got <- vapply(
    c("none", "unit", "stage"),
    function(mode) {
      fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     fpc = mode)
      pred <- predict(fit, data.frame(n_per_psu = 20))
      # Re-solving at the fitted take must return the fitted design.
      expect_equal(pred$n, fit$detail$n, tolerance = 1e-6)
      pred$n
    },
    numeric(1)
  )
  # And the three must stay distinguishable, or a dropped mode reads as
  # agreement rather than as a bug.
  expect_gt(got[["none"]], got[["unit"]])
  expect_gt(got[["unit"]], got[["stage"]])
})

test_that("the generalized take sweep keeps the correction across the grid", {
  z <- .fpc_joint()
  swept <- lapply(
    c("none", "unit", "stage"),
    function(mode) {
      fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     fpc = mode)
      predict(fit, data.frame(n_per_psu = c(10, 20, 30)))$n
    }
  )
  # Every take in the grid, not only the fitted one.
  expect_true(all(swept[[1]] > swept[[3]]))
  expect_true(all(swept[[2]] > swept[[3]]))
})

test_that("the budget frontier floor is found under the fit's correction", {
  z <- .fpc_joint()
  floors <- vapply(
    c("none", "unit", "stage"),
    function(mode) {
      fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     objective = "y", budget = 4e5, fpc = mode)
      grid <- .budget_grid(fit, npoints = 5L)
      grid[1L]
    },
    numeric(1)
  )
  # A cheaper correction starts the frontier lower, because the cheapest
  # design meeting the targets is itself cheaper.
  expect_gt(floors[["none"]], floors[["stage"]])
  expect_gt(floors[["unit"]], floors[["stage"]])
})

## T11. domain_sampling, the same class found by dev/api-gate.R
##
## Not an fpc test. It lives here because it is the same defect: a mode stored
## on the fit and dropped by the round trip, so n_multi() silently re-solved
## under "separate". prec_multi() carried a fixed list of parameters forward
## and this one was not on it.

test_that("domain_sampling survives the n_multi round trip", {
  indicators <- data.frame(
    name = c("a", "b"), p = c(0.3, 0.5), cv = c(0.08, 0.06),
    domain = c("north", "south"), share = c(0.4, 0.6)
  )
  sizes <- vapply(
    c("separate", "natural"),
    function(mode) {
      fit <- n_multi(indicators, domains = "domain", domain_sampling = mode)
      expect_identical(fit$params$domain_sampling, mode)
      expect_identical(prec_multi(fit)$params$domain_sampling, mode)
      expect_equal(n_multi(prec_multi(fit))$n, fit$n, tolerance = 1e-8)
      fit$n
    },
    numeric(1)
  )
  # The two modes must stay apart, or the round-trip check above passes on a
  # design that has quietly reverted.
  expect_gt(sizes[["natural"]], sizes[["separate"]])
})
