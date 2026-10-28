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

## T11. domain_sampling, the same defect class as fpc
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

## T12. A free take under "stage" is the optimum of the corrected variance
##
## With D = icc - (1 - icc) / M, the corrected variance per unit of budget is
## minimized at sqrt(c1 (1 - icc) / (c2 r D)), and at the whole PSU when
## D <= 0. Taking the uncorrected take instead loses up to a fifth of the
## variance when PSUs are small.

.stage_take_frame <- function(icc = 5 / 54, N = 4900, N_psu = 100) {
  data.frame(
    N = N, N_psu = N_psu, sd = sqrt(54), mean = 10,
    icc_psu = icc, cost_psu = 12, cost_ssu = 3
  )
}

.stage_take_var <- function(b, budget, frame, r = 1) {
  M <- frame$N / frame$N_psu
  S2 <- frame$sd^2
  Sb <- S2 * frame$icc_psu
  Sw <- S2 * (1 - frame$icc_psu)
  a <- budget / (frame$cost_psu + frame$cost_ssu * b)
  (1 - a / frame$N_psu) * Sb / a + (1 - b * r / M) * Sw / (a * b * r)
}

test_that("the stage take minimizes the corrected variance at fixed budget", {
  frame <- .stage_take_frame()
  fit <- n_alloc(frame, budget = 1000, fpc = "stage")
  best <- optimize(
    .stage_take_var, c(1, 49), budget = 1000, frame = frame, tol = 1e-10
  )$minimum
  expect_equal(fit$detail$n_per_psu, 7, tolerance = 1e-10)
  expect_equal(fit$detail$n_per_psu, best, tolerance = 1e-6)
  expect_equal(fit$detail$n_psu, 1000 / 33, tolerance = 1e-10)
  expect_equal(fit$se^2, .stage_take_var(7, 1000, frame), tolerance = 1e-10)
})

test_that("the stage take reads the unit response rate", {
  frame <- .stage_take_frame()
  r <- 0.8
  fit <- n_alloc(frame, budget = 1000, resp_rate = r, fpc = "stage")
  D <- 5 / 54 - (49 / 54) / 49
  expect_equal(
    fit$detail$n_per_psu,
    sqrt(12 * (49 / 54) / (3 * r * D)),
    tolerance = 1e-10
  )
})

test_that("below 1 / (M + 1) the stage take is the whole PSU", {
  frame <- .stage_take_frame(icc = 0.015)
  expect_no_warning(fit <- n_alloc(frame, budget = 1000, fpc = "stage"))
  expect_equal(fit$detail$n_per_psu, 49)
})

test_that("unit and none keep the uncorrected take", {
  frame <- .stage_take_frame()
  uncorrected <- sqrt(12 / 3 * (49 / 54) / (5 / 54))
  for (f in c("unit", "none")) {
    fit <- n_alloc(frame, budget = 1000, fpc = f)
    expect_equal(fit$detail$n_per_psu, uncorrected, tolerance = 1e-10)
  }
})

test_that("the stage take survives the round trip through prec_alloc", {
  fit <- n_alloc(.stage_take_frame(), budget = 1000, fpc = "stage")
  back <- n_alloc(prec_alloc(fit))
  expect_equal(back$detail$n_per_psu, fit$detail$n_per_psu, tolerance = 1e-10)
  expect_equal(back$detail$n, fit$detail$n, tolerance = 1e-8)
})

## T13. n_cluster() and prec_cluster() under the stage correction
##
## Gross PSU count a and take b give
##   cv^2 = [B - W/M + W/(b r)] / (a r1) - B/A,  B = U k icc, W = U k (1 - icc)
## A design never takes more than N_psu PSUs or M = N / N_psu units per PSU.

.st_args <- list(
  stage_cost = c(12, 3), icc = 5 / 54, unit_relvar = 54 / 100,
  fpc = "stage", N = 4900, N_psu = 100
)

.st_cv2 <- function(a, b, icc = 5 / 54, U = 0.54, A = 100, M = 49, r1 = 1,
                    r = 1) {
  U * icc * (1 - a * r1 / A) / (a * r1) +
    U * (1 - icc) * (1 - b * r / M) / (a * r1 * b * r)
}

test_that("prec_cluster applies each stage's own correction", {
  got <- prec_cluster(
    n = c(20, 10), icc = 0.1, unit_relvar = 2, fpc = "stage",
    N = 4900, N_psu = 100, resp_rate = 0.8, resp_rate_psu = 0.9
  )$cv
  expect_equal(got^2, .st_cv2(20, 10, icc = 0.1, U = 2, r1 = 0.9, r = 0.8),
               tolerance = 1e-12)
})

test_that("the stage take minimizes the corrected variance in n_cluster", {
  fit <- do.call(n_cluster, c(.st_args, list(budget = 1000)))
  expect_equal(fit$n[["n_per_psu"]], 7, tolerance = 1e-10)
  expect_equal(fit$n[["n_psu"]], 1000 / 33, tolerance = 1e-10)
  expect_equal(fit$cv^2, .st_cv2(1000 / 33, 7), tolerance = 1e-12)
})

test_that("below 1 / (M + 1) n_cluster takes whole PSUs", {
  fit <- do.call(
    n_cluster, modifyList(.st_args, list(icc = 0.015, budget = 1000))
  )
  expect_equal(fit$n[["n_per_psu"]], 49)
})

test_that("past N_psu PSUs the budget goes to the take", {
  fit <- do.call(n_cluster, c(.st_args, list(budget = 6000)))
  expect_equal(fit$n[["n_psu"]], 100)
  expect_equal(fit$n[["n_per_psu"]], (6000 / 100 - 12) / 3, tolerance = 1e-10)
  expect_lte(fit$operational$n[["n_psu"]], 100)
  expect_lte(fit$operational$n[["n_per_psu"]], 49)
})

test_that("past N_psu PSUs the precision target goes to the take", {
  fit <- do.call(n_cluster, c(.st_args, list(cv = 0.01)))
  expect_equal(fit$n[["n_psu"]], 100)
  expect_equal(fit$cv, 0.01, tolerance = 1e-10)
  expect_lte(fit$n[["n_per_psu"]], 49)
  expect_lte(fit$operational$cv, 0.01 + 1e-12)
  expect_lte(fit$operational$n[["n_psu"]], 100)
})

test_that("a fixed take that needs more than N_psu PSUs is refused", {
  expect_error(
    do.call(n_cluster, c(.st_args, list(cv = 0.01, n_per_psu = 5))),
    "A larger take may reach it"
  )
})

test_that("population bounds are checked on fixed sizes", {
  expect_error(
    do.call(n_cluster, c(.st_args, list(budget = 1e5, n_psu = 120))),
    "exceeds 'N_psu'"
  )
  expect_error(
    prec_cluster(n = c(20, 60), icc = 0.1, fpc = "stage", N = 4900,
                 N_psu = 100),
    "exceeds the units per PSU"
  )
})

test_that("the stage correction needs its population and only it reads it", {
  expect_error(
    n_cluster(stage_cost = c(12, 3), icc = 0.1, budget = 1000, N = 4900),
    "fpc = \"stage\" corrects for"
  )
  expect_error(
    n_cluster(stage_cost = c(12, 3), icc = 0.1, budget = 1000, fpc = "stage"),
    "needs the population"
  )
  plan <- svyplan(N = 4900)
  expect_no_error(
    n_cluster(stage_cost = c(12, 3), icc = 0.1, budget = 1000, plan = plan)
  )
  expect_equal(
    n_cluster(stage_cost = c(12, 3), icc = 0.1, budget = 1000, plan = plan,
              fpc = "stage", N_psu = 100)$n,
    n_cluster(stage_cost = c(12, 3), icc = 0.1, budget = 1000,
              fpc = "stage", N = 4900, N_psu = 100)$n
  )
})

test_that("the stage correction is refused where its model does not hold", {
  fr <- data.frame(
    y = c(5, 10, 12, 3, 6, 15, 13, 6, 11, 4, 8, 9),
    psu = rep(1:4, each = 3),
    share = rep(c(0.1, 0.2, 0.3, 0.4), each = 3)
  )
  pps <- varcomp(y ~ psu, data = fr, prob = ~share)
  expect_error(
    prec_cluster(n = c(2, 2), icc = pps, fpc = "stage", N = 12, N_psu = 4),
    "size-proportional first stage"
  )
  fr3 <- data.frame(
    y = c(5, 10, 12, 3, 6, 15, 13, 6, 11, 4, 8, 9, 7, 2, 5, 9),
    psu = rep(1:2, each = 8),
    ssu = rep(1:4, each = 4),
    share = rep(c(0.4, 0.6, 0.3, 0.7), each = 4)
  )
  pps_ssu <- varcomp(fr3$y, stage_id = list(fr3$psu, fr3$ssu),
                     prob = list(NULL, fr3$share))
  expect_error(
    prec_cluster(n = c(2, 2, 2), icc = pps_ssu, fpc = "stage", N = 16,
                 N_psu = 2, N_ssu = 4),
    "size-proportional second stage"
  )
  expect_error(
    n_cluster(
      indicators = data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
                              region = c("a", "b")),
      domains = "region", stage_cost = c(12, 3), fpc = "stage", N = 4900,
      N_psu = 100
    ),
    "reads each domain's population from 'indicators' columns"
  )
})

test_that("three stages read N_ssu and two stages refuse it", {
  expect_error(
    n_cluster(stage_cost = c(12, 3, 1), icc = c(0.1, 0.1), budget = 1000,
              fpc = "stage", N = 4900, N_psu = 100),
    "also needs 'N_ssu'"
  )
  expect_error(
    n_cluster(stage_cost = c(12, 3), icc = 0.1, budget = 1000,
              fpc = "stage", N = 4900, N_psu = 100, N_ssu = 700),
    "applies to three-stage designs"
  )
  expect_error(
    prec_cluster(n = c(10, 2, 3), icc = c(0.1, 0.1), N_ssu = 700),
    "fpc = \"stage\" corrects for"
  )
  expect_error(
    prec_cluster(n = c(10, 2, 3), icc = c(0.1, 0.1), fpc = "stage",
                 N = 4900, N_psu = 100, N_ssu = 50),
    "between 'N_psu' and 'N'"
  )
})

test_that("the stage correction survives every round trip", {
  fit <- do.call(n_cluster, c(.st_args, list(budget = 1000)))
  prec <- prec_cluster(fit)
  expect_equal(prec$cv, fit$cv, tolerance = 1e-12)
  expect_identical(prec$params$fpc, "stage")
  back <- n_cluster(prec)
  expect_equal(back$n, fit$n, tolerance = 1e-10)
  swept <- predict(fit, data.frame(budget = c(800, 1000)))
  expect_equal(swept$n_per_psu[2], 7, tolerance = 1e-10)
  expect_equal(swept$cv[2], fit$cv, tolerance = 1e-12)
})

## T14. Three stages under the stage correction
##
## With R1 = U k1 d1, R2 = U k2 d2, R3 = U k2 (1 - d2) and responding sizes,
##   cv^2 = R1 (1 - a/A)/a + R2 (1 - b/Q)/(ab) + R3 (1 - c/C)/(abc)
## for A PSUs of Q SSUs of C units. The SSU fraction lowers the PSU term to
## R1 - R2/Q and the unit fraction lowers the SSU term to R2 - R3/C.

.st3_args <- list(
  stage_cost = c(300, 50, 10), icc = c(0.04, 0.10), unit_relvar = 7 / 3,
  fpc = "stage", N = 12000, N_psu = 40, N_ssu = 400
)

.st3_cv2 <- function(a, b, cc, icc = c(0.04, 0.10), U = 7 / 3,
                     k = c(1, 0.96), A = 40, Q = 10, C = 30) {
  U * k[1] * icc[1] * (1 - a / A) / a +
    U * k[2] * icc[2] * (1 - b / Q) / (a * b) +
    U * k[2] * (1 - icc[2]) * (1 - cc / C) / (a * b * cc)
}

## The cheapest design for the target over a grid of takes, within the bounds.
.st3_grid_cost <- function(cv, icc = c(0.04, 0.10), Q = 10, C = 30,
                           cost = c(300, 50, 10), A = 40, k = c(1, 0.96)) {
  g <- expand.grid(b = seq(1, Q, length.out = 200),
                   cc = seq(1, C, length.out = 200))
  U <- 7 / 3
  R <- c(U * k[1] * icc[1], U * k[2] * icc[2], U * k[2] * (1 - icc[2]))
  G <- R[1] - R[2] / Q + (R[2] - R[3] / C) / g$b + R[3] / (g$b * g$cc)
  a <- G / (cv^2 + R[1] / A)
  ok <- a <= A
  min((a * (cost[1] + cost[2] * g$b + cost[3] * g$b * g$cc))[ok])
}

test_that("prec_cluster corrects each of three stages by its own fraction", {
  got <- do.call(prec_cluster, c(
    .st3_args[c("icc", "unit_relvar", "fpc", "N", "N_psu", "N_ssu")],
    list(n = c(30, 4, 8))
  ))$cv
  expect_equal(got^2, .st3_cv2(30, 4, 8), tolerance = 1e-12)
})

test_that("the three-stage stage design is no costlier than any on a grid", {
  for (cv in c(0.05, 0.08)) {
    fit <- do.call(n_cluster, c(.st3_args, list(cv = cv)))
    expect_equal(fit$cv, cv, tolerance = 1e-10)
    expect_lte(fit$cost, .st3_grid_cost(cv))
    expect_equal(fit$cv^2, .st3_cv2(fit$n[[1]], fit$n[[2]], fit$n[[3]]),
                 tolerance = 1e-12)
  }
})

## The best design over a grid of takes within the population, the cheapest
## for a target or the most precise for a budget.
.st3_grid_best <- function(icc, A, Q, C, cv = NULL, budget = NULL,
                           cost = c(300, 50, 10), U = 7 / 3) {
  g <- expand.grid(b = seq(1, Q, length.out = 300),
                   cc = seq(1, C, length.out = 300))
  k2 <- 1 - icc[1]
  R <- c(U * icc[1], U * k2 * icc[2], U * k2 * (1 - icc[2]))
  G <- R[1] - R[2] / Q + (R[2] - R[3] / C) / g$b + R[3] / (g$b * g$cc)
  H <- cost[1] + cost[2] * g$b + cost[3] * g$b * g$cc
  if (!is.null(cv)) {
    a <- G / (cv^2 + R[1] / A)
    min((a * H)[a <= A])
  } else {
    a <- pmin(A, budget / H)
    sqrt(min((G / a - R[1] / A)[a >= 1]))
  }
}

.st3_corner <- function(icc, N, N_ssu, ...) {
  n_cluster(stage_cost = c(300, 50, 10), icc = icc, unit_relvar = 7 / 3,
            fpc = "stage", N = N, N_psu = 40, N_ssu = N_ssu, ...)
}

test_that("a unit fraction at or past R2 / R3 = 1 / C takes whole SSUs", {
  ## icc_ssu = 0.1 is below 1 / (C + 1) = 1 / 6, so R2 - R3 / C < 0.
  icc <- c(0.04, 0.10)
  by_budget <- .st3_corner(icc, 2000, 400, budget = 20000)
  expect_equal(by_budget$n[["n_per_ssu"]], 5)
  expect_lte(by_budget$cv,
             .st3_grid_best(icc, 40, 10, 5, budget = 20000) + 1e-12)
  by_cv <- .st3_corner(icc, 2000, 400, cv = 0.05)
  expect_equal(by_cv$n[["n_per_ssu"]], 5)
  expect_lte(by_cv$cost, .st3_grid_best(icc, 40, 10, 5, cv = 0.05) + 1e-8)
  for (fit in list(by_budget, by_cv)) {
    expect_lte(fit$operational$n[["n_per_ssu"]], 5)
  }
})

test_that("an SSU fraction at or past R1 = R2 / Q takes every SSU of a PSU", {
  ## R1 - R2 / Q = U (0.05 - 0.95 * 0.5 / 4) < 0.
  icc <- c(0.05, 0.5)
  by_budget <- .st3_corner(icc, 4800, 160, budget = 20000)
  expect_equal(by_budget$n[["n_per_psu"]], 4)
  expect_lte(by_budget$cv,
             .st3_grid_best(icc, 40, 4, 30, budget = 20000) + 1e-12)
  by_cv <- .st3_corner(icc, 4800, 160, cv = 0.08)
  expect_equal(by_cv$n[["n_per_psu"]], 4)
  expect_lte(by_cv$cost, .st3_grid_best(icc, 40, 4, 30, cv = 0.08) + 1e-8)
  for (fit in list(by_budget, by_cv)) {
    expect_lte(fit$operational$n[["n_per_psu"]], 4)
  }
})

test_that("past N_psu PSUs the three-stage budget and target go to the takes", {
  fit <- do.call(n_cluster, c(.st3_args, list(budget = 1e5)))
  expect_equal(fit$n[["n_psu"]], 40)
  expect_equal(fit$cost, 1e5, tolerance = 1e-10)
  expect_lte(fit$n[["n_per_psu"]], 10)
  expect_lte(fit$n[["n_per_ssu"]], 30)
  op <- fit$operational$n
  expect_true(op[["n_psu"]] <= 40 && op[["n_per_psu"]] <= 10 &&
                op[["n_per_ssu"]] <= 30)
  expect_lte(fit$operational$cost, 1e5)

  tight <- do.call(n_cluster, c(.st3_args, list(cv = 0.02)))
  expect_equal(tight$n[["n_psu"]], 40)
  expect_equal(tight$cv, 0.02, tolerance = 1e-10)
  expect_lte(tight$operational$cv, 0.02 + 1e-12)
  op <- tight$operational$n
  expect_true(op[["n_psu"]] <= 40 && op[["n_per_psu"]] <= 10 &&
                op[["n_per_ssu"]] <= 30)
  ## At full response every unit taken is a census, so only nonresponse
  ## leaves a floor above zero.
  expect_error(
    do.call(n_cluster, c(.st3_args, list(cv = 0.001, resp_rate = 0.8))),
    "out of reach"
  )
})

test_that("three-stage stage sizes are bounded by the population", {
  expect_error(
    do.call(n_cluster, c(.st3_args, list(cv = 0.05, n_per_psu = 11))),
    "exceeds the SSUs per PSU"
  )
  expect_error(
    do.call(n_cluster, c(.st3_args, list(cv = 0.05, n_per_ssu = 31))),
    "exceeds the units per SSU"
  )
  expect_error(
    do.call(n_cluster, c(.st3_args, list(budget = 1e5, n_psu = 41))),
    "exceeds 'N_psu'"
  )
})

test_that("the three-stage stage correction survives every round trip", {
  fit <- do.call(n_cluster, c(.st3_args, list(cv = 0.06, resp_rate = 0.8)))
  prec <- prec_cluster(fit)
  expect_equal(prec$cv, fit$cv, tolerance = 1e-12)
  expect_identical(prec$params$N_ssu, 400)
  back <- n_cluster(prec)
  expect_equal(back$n, fit$n, tolerance = 1e-10)
  swept <- predict(fit, data.frame(N_ssu = c(400, 800), cv = 0.06))
  expect_equal(swept$n_psu[1], fit$n[["n_psu"]], tolerance = 1e-10)
  expect_false(isTRUE(all.equal(swept$n_psu[2], swept$n_psu[1])))
})

## T15. N counts final-stage units, and eligible units are N * eligible_rate
##
## The take is bounded by the N / N_psu final-stage units a PSU holds, and
## the within-PSU fraction is the responding share of its eligible units.

test_that("the take counts final-stage units and the variance eligible ones", {
  args <- list(stage_cost = c(12, 3), icc = 0.015, unit_relvar = 0.54,
               budget = 1000, fpc = "stage", N = 4900, N_psu = 100)
  ## 24.5 eligible units per PSU, and 0.015 < 1 / 25.5, so whole PSUs.
  half <- do.call(n_cluster, c(args, list(eligible_rate = 0.5)))
  expect_equal(half$n[["n_per_psu"]], 49)
  expect_lte(half$operational$n[["n_per_psu"]], 49)
  expect_identical(half$params$N, 4900)
  ## 98 eligible units per PSU, and 0.015 > 1 / 99, so the corrected take.
  double <- do.call(n_cluster, c(args, list(eligible_rate = 2)))
  D <- 0.015 - 0.985 / 98
  expect_equal(double$n[["n_per_psu"]], sqrt(12 * 0.985 / (3 * 2 * D)),
               tolerance = 1e-10)
  expect_equal(
    prec_cluster(n = c(10, 30), icc = 0.1, unit_relvar = 2, fpc = "stage",
                 N = 4900, N_psu = 100, eligible_rate = 2)$cv^2,
    .st_cv2(10, 30, icc = 0.1, U = 2, M = 98, r = 2),
    tolerance = 1e-12
  )
  expect_error(
    prec_cluster(n = c(10, 60), icc = 0.1, fpc = "stage", N = 4900,
                 N_psu = 100, eligible_rate = 0.5),
    "exceeds the units per PSU"
  )
})

## T16. Several indicators under the stage correction
##
## Every row reads the one population. A one-row table returns the scalar
## design, which test-cluster-interface-equivalence.R pins. Here, rows that
## differ meet their own corrected targets.

.st_rows <- data.frame(
  name = c("stunting", "anemia"), p = c(0.3, 0.1), cv = c(0.08, 0.12),
  icc_psu = c(0.03, 0.12)
)

test_that("every row meets its corrected target and one binds", {
  fit <- n_cluster(indicators = .st_rows, stage_cost = c(120, 10),
                   fpc = "stage", N = 3000, N_psu = 60)
  each <- vapply(seq_len(2), function(j) {
    prec_cluster(n = fit$n, icc = .st_rows$icc_psu[j],
                 unit_relvar = (1 - .st_rows$p[j]) / .st_rows$p[j],
                 fpc = "stage", N = 3000, N_psu = 60)$cv
  }, numeric(1))
  expect_equal(fit$detail$.cv_achieved, each, tolerance = 1e-12)
  expect_true(all(each <= .st_rows$cv * (1 + 1e-8)))
  expect_equal(max(each / .st_rows$cv), 1, tolerance = 1e-6)
  expect_lte(fit$n[["n_psu"]], 60)
  expect_true(all(fit$operational$cv_by_target <= .st_rows$cv + 1e-12))
})

test_that("prec_cluster reads the indicators under the stage correction", {
  fit <- n_cluster(indicators = .st_rows, stage_cost = c(120, 10),
                   fpc = "stage", N = 3000, N_psu = 60, resp_rate = 0.8)
  prec <- prec_cluster(fit)
  expect_identical(prec$params$fpc, "stage")
  expect_equal(prec$detail$.cv, fit$detail$.cv_achieved, tolerance = 1e-10)
  back <- n_cluster(prec)
  expect_equal(back$n, fit$n, tolerance = 1e-7)
  expect_identical(back$params$N_psu, 60)
  expect_output(print(fit), "corrected by stage for N_psu = 60 PSUs")
  expect_output(print(prec), "corrected by stage for N_psu = 60 PSUs")
  three <- do.call(n_cluster, c(.st3_args, list(cv = 0.06)))
  expect_output(print(three), "N_ssu = 400 SSUs of N = 12000 units")
})

test_that("three-stage indicators stay within the population", {
  rows <- transform(.st_rows, icc_ssu = c(0.08, 0.2))
  fit <- n_cluster(indicators = rows, stage_cost = c(300, 50, 10),
                   cv = NULL, budget = 60000, fpc = "stage",
                   N = 12000, N_psu = 40, N_ssu = 400)
  expect_lte(fit$n[["n_psu"]], 40)
  expect_lte(fit$n[["n_per_psu"]], 10)
  expect_lte(fit$n[["n_per_ssu"]], 30)
  op <- fit$operational
  expect_true(op$n[["n_psu"]] <= 40 && op$n[["n_per_psu"]] <= 10 &&
                op$n[["n_per_ssu"]] <= 30)
  expect_lte(op$cost, 60000)
})

test_that("the whole-unit three-stage indicators design is the best affordable", {
  ## Ten PSUs, so each row's correction R1 / A is large enough to reorder
  ## the designs the two targets prefer.
  rows <- data.frame(name = c("a", "b"), p = c(0.16, 0.48), cv = c(0.12, 0.06),
                     icc_psu = c(0.15, 0.06), icc_ssu = c(0.06, 0.27))
  fit <- n_cluster(indicators = rows, stage_cost = c(300, 50, 10),
                   budget = 8100, fpc = "stage", N = 720, N_psu = 10,
                   N_ssu = 60)
  g <- expand.grid(a = 1:10, m = 1:6, q = 1:12)
  g <- g[g$a * (300 + 50 * g$m + 10 * g$m * g$q) <= 8100, ]
  worst <- rep(0, nrow(g))
  for (j in 1:2) {
    cv <- vapply(seq_len(nrow(g)), function(i) {
      prec_cluster(n = c(g$a[i], g$m[i], g$q[i]),
                   icc = c(rows$icc_psu[j], rows$icc_ssu[j]),
                   unit_relvar = (1 - rows$p[j]) / rows$p[j], fpc = "stage",
                   N = 720, N_psu = 10, N_ssu = 60)$cv
    }, numeric(1))
    worst <- pmax(worst, cv / rows$cv[j])
  }
  expect_equal(max(fit$operational$cv_by_target / rows$cv), min(worst),
               tolerance = 1e-12)
})

## T17. Domains under the stage correction
##
## Each domain is a population of its own, given by N, N_psu and N_ssu
## columns constant within it. Separate domains are the single-population
## design per domain, and a joint budget equalizes the worst ratio.

.dm_rows <- data.frame(name = c("stunting", "anemia"), p = c(0.30, 0.10),
                       cv = c(0.08, 0.12), icc_psu = c(0.03, 0.12))
.dm_pops <- list(capital = c(N = 4500, N_psu = 30),
                 rural = c(N = 135000, N_psu = 900))
.dm_table <- function(rows = .dm_rows, pops = .dm_pops) {
  do.call(rbind, lapply(names(pops), function(d) {
    cbind(rows, region = d, as.data.frame(as.list(pops[[d]])))
  }))
}
.dm_one <- function(d, rows = .dm_rows, pops = .dm_pops, ...) {
  args <- c(list(indicators = rows, fpc = "stage"), as.list(pops[[d]]),
            list(...))
  do.call(n_cluster, args)
}

test_that("separate domains are each domain's own design", {
  for (mode in list(list(), list(budget = 20000))) {
    fit <- do.call(n_cluster, c(list(indicators = .dm_table(),
                                     stage_cost = c(120, 10),
                                     domains = "region", fpc = "stage"),
                                mode))
    for (i in 1:2) {
      one <- do.call(.dm_one, c(list(names(.dm_pops)[i],
                                     stage_cost = c(120, 10)), mode))
      expect_identical(unname(unlist(fit$domains[i, c("n_psu", "n_per_psu")])),
                       unname(one$n))
      expect_identical(fit$domains$.cv[i], one$cv)
    }
  }
  ## The capital of 30 EAs takes all of them, where the uncorrected plan
  ## asks for 134.
  capital <- .dm_one("capital", stage_cost = c(120, 10))
  expect_equal(capital$n[["n_psu"]], 30)
  expect_equal(capital$cost, 8501, tolerance = 1e-3)
})

test_that("three-stage domains read N_ssu columns", {
  rows <- transform(.dm_rows, icc_ssu = c(0.08, 0.2))
  pops <- list(a = c(N = 12000, N_psu = 40, N_ssu = 400),
               b = c(N = 60000, N_psu = 200, N_ssu = 2000))
  fit <- n_cluster(indicators = .dm_table(rows, pops),
                   stage_cost = c(300, 50, 10), domains = "region",
                   fpc = "stage")
  one <- .dm_one("a", rows, pops, stage_cost = c(300, 50, 10))
  expect_identical(unname(unlist(fit$domains[1, c("n_psu", "n_per_psu",
                                                  "n_per_ssu")])),
                   unname(one$n))
  expect_error(
    n_cluster(indicators = .dm_table(rows, lapply(pops, `[`, 1:2)),
              stage_cost = c(300, 50, 10), domains = "region",
              fpc = "stage"),
    "also needs 'N_ssu'"
  )
})

test_that("a joint budget equalizes the corrected worst ratio", {
  fit <- n_cluster(indicators = .dm_table(), stage_cost = c(120, 10),
                   domains = "region", fpc = "stage", budget = 25000,
                   allocation = "joint")
  ratio <- vapply(1:2, function(i) {
    max(prec_cluster(
      indicators = cbind(.dm_rows, n = fit$domains$n_psu[i],
                         n_per_psu = fit$domains$n_per_psu[i]),
      fpc = "stage", N = .dm_pops[[i]][["N"]], N_psu = .dm_pops[[i]][["N_psu"]]
    )$cv / .dm_rows$cv)
  }, numeric(1))
  expect_equal(ratio[1], ratio[2], tolerance = 1e-3)
  expect_equal(fit$cost, 25000, tolerance = 1e-8)
  expect_lte(fit$domains$n_psu[1], 30)
  ## No split of the budget does better on the worst domain.
  grid <- vapply(seq(0.05, 0.95, by = 0.01), function(w) {
    max(vapply(1:2, function(i) {
      share <- if (i == 1) w else 1 - w
      one <- .dm_one(names(.dm_pops)[i], stage_cost = c(120, 10),
                     budget = share * 25000)
      max(one$detail$.cv_achieved / .dm_rows$cv)
    }, numeric(1)))
  }, numeric(1))
  expect_lte(max(ratio), min(grid) + 1e-6)

  single <- .dm_table(pops = .dm_pops["capital"])
  by_joint <- n_cluster(indicators = single, stage_cost = c(120, 10),
                        domains = "region", fpc = "stage", budget = 8000,
                        allocation = "joint")
  by_separate <- n_cluster(indicators = single, stage_cost = c(120, 10),
                           domains = "region", fpc = "stage", budget = 8000)
  expect_equal(by_joint$domains$.cv, by_separate$domains$.cv,
               tolerance = 1e-4)
})

test_that("precision with domains reads each domain's population", {
  fit <- n_cluster(indicators = .dm_table(), stage_cost = c(120, 10),
                   domains = "region", fpc = "stage")
  prec <- prec_cluster(fit)
  expect_identical(prec$params$fpc, "stage")
  for (i in 1:2) {
    scalar <- vapply(1:2, function(j) {
      prec_cluster(
        n = unlist(fit$domains[i, c("n_psu", "n_per_psu")]),
        icc = .dm_rows$icc_psu[j],
        unit_relvar = (1 - .dm_rows$p[j]) / .dm_rows$p[j], fpc = "stage",
        N = .dm_pops[[i]][["N"]], N_psu = .dm_pops[[i]][["N_psu"]]
      )$cv
    }, numeric(1))
    expect_equal(prec$detail$.cv[prec$detail$region == names(.dm_pops)[i]],
                 scalar, tolerance = 1e-12)
  }
  back <- n_cluster(prec)
  expect_equal(back$domains, fit$domains, tolerance = 1e-8)
  expect_output(print(fit),
                "corrected by stage for each domain's population")
  budgeted <- n_cluster(indicators = .dm_table(), stage_cost = c(120, 10),
                        domains = "region", fpc = "stage", budget = 20000)
  expect_output(print(budgeted), "budget = 20000 for each domain")
  expect_equal(budgeted$cost, 40000)
})

test_that("population columns are read only under the stage correction", {
  expect_error(
    n_cluster(indicators = transform(.dm_rows, N = 500, N_psu = 20),
              stage_cost = c(120, 10)),
    "Add fpc = \"stage\", or drop them"
  )
  expect_error(
    n_cluster(indicators = .dm_table(), stage_cost = c(120, 10),
              domains = "region"),
    "Add fpc = \"stage\""
  )
  expect_error(
    prec_cluster(indicators = transform(.dm_rows, n = 20, n_per_psu = 10,
                                        N_psu = 30)),
    "Add fpc = \"stage\", or drop it"
  )
  expect_error(
    n_cluster(indicators = .dm_table(), stage_cost = c(120, 10),
              domains = "region", fpc = "stage", N = 4500, N_psu = 30),
    "given twice"
  )
  varying <- .dm_table()
  varying$N_psu[1] <- 31
  expect_error(
    n_cluster(indicators = varying, stage_cost = c(120, 10),
              domains = "region", fpc = "stage"),
    "constant within each domain"
  )
  whole <- n_cluster(indicators = transform(.dm_rows, N = 4500, N_psu = 30),
                     stage_cost = c(120, 10), fpc = "stage")
  by_args <- n_cluster(indicators = .dm_rows, stage_cost = c(120, 10),
                       fpc = "stage", N = 4500, N_psu = 30)
  expect_identical(whole$n, by_args$n)
  expect_null(whole$params$N_psu)
  expect_identical(prec_cluster(whole)$detail$.cv,
                   prec_cluster(by_args)$detail$.cv)
})
