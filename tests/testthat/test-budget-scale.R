## Budgets at any monetary scale

scale_families <- list(
  `free take` = function(k) {
    n_alloc_joint(
      data.frame(stratum = c("A", "B"), N = c(621, 400), N_1 = c(27, 20),
                 cost_1 = c(17, 20) * k, cost_2 = c(2.56, 3) * k),
      data.frame(stratum = c("A", "B"), name = "y", p = c(.37, .2),
                 icc_1 = c(.31, .1)),
      objective = "y", budget = 500 * k
    )
  },
  `fixed take` = function(k) {
    n_alloc_joint(
      data.frame(stratum = c("A", "B"), N = c(621, 400), N_1 = c(27, 20),
                 take_2 = c(5, 4), cost_1 = c(17, 20) * k,
                 cost_2 = c(2.56, 3) * k),
      data.frame(stratum = c("A", "B"), name = "y", p = c(.37, .2),
                 icc_1 = c(.31, .1)),
      objective = "y", budget = 500 * k
    )
  },
  `two stages` = function(k) {
    n_alloc(data.frame(stratum = c("A", "B"), N = c(6000, 4000),
                       N_1 = c(60, 40), mean = c(10, 12), sd = c(3, 5),
                       icc_1 = c(.1, .2), cost_1 = c(50, 70) * k,
                       cost_2 = c(4, 6) * k),
            budget = 3000 * k, alloc = "optimal")
  },
  element = function(k) {
    n_alloc(data.frame(stratum = c("A", "B", "C"), N = c(6000, 4000, 2000),
                       mean = c(10, 12, 9), sd = c(3, 5, 4),
                       unit_cost = c(4, 6, 9) * k),
            budget = 1500 * k, alloc = "optimal")
  },
  `three stages` = function(k) {
    n_alloc(data.frame(stratum = c("A", "B"), N = c(20000, 30000),
                       N_1 = c(200, 300), N_2 = c(2000, 3000),
                       mean = c(40, 50), sd = c(10, 14), icc_1 = c(.05, .1),
                       icc_2 = c(.15, .2), cost_1 = c(400, 600) * k,
                       cost_2 = c(50, 40) * k, cost_3 = c(10, 12) * k),
            budget = 20000 * k, alloc = "optimal")
  },
  `two phases` = function(k) {
    n_twophase(data.frame(stratum = c("A", "B", "C", "D"),
                          N = c(3500, 2500, 2500, 1500),
                          sd = c(12, 25, 8, 40), mean = c(40, 70, 35, 90),
                          unit_cost = c(2, 5, 1, 9) * k),
               phase1_cost = k, budget = 50000 * k)
  }
)

budget_of <- c(`free take` = 500, `fixed take` = 500, `two stages` = 3000,
               element = 1500, `three stages` = 20000, `two phases` = 50000)

test_that("a common cost scale changes neither the design nor its spending", {
  checked <- character(0)
  for (nm in names(scale_families)) {
    at <- lapply(c(1e-12, 1, 1e12), function(k) {
      fit <- suppressWarnings(scale_families[[nm]](k))
      cont <- fit$params$achieved$cost %||% fit$cost
      c(whole = fit$operational$cost / k, cont = cont / k,
        objective = fit$operational$cv %||% fit$objective_value)
    })
    for (a in at) {
      expect_lte(a[["whole"]], budget_of[[nm]] * (1 + 1e-6))
      expect_lte(a[["cont"]], budget_of[[nm]] * (1 + 1e-6))
      expect_equal(a[["objective"]], at[[2]][["objective"]],
                   tolerance = 1e-9, info = nm)
    }
    checked <- c(checked, nm)
  }
  expect_setequal(checked, names(budget_of))
})

test_that("a continuous free-take design stays within its budget", {
  fb <- data.frame(stratum = c("A", "B"), N = c(600, 400), N_1 = c(20, 10),
                   cost_1 = c(10, 12), cost_2 = c(1, 2))
  for (icc in c(.9, .99, .999999)) {
    mb <- data.frame(stratum = c("A", "B"), name = "y", p = .37,
                     icc_1 = icc)
    fit <- n_alloc_joint(fb, mb, objective = "y", budget = 600)
    expect_lte(fit$params$achieved$cost, 600 * (1 + 1e-6))
    expect_lte(fit$operational$cost, 600)
    expect_identical(fit$optimization$classification, "feasible")
  }
})

test_that("a budget fit is optimal only within its precision bound", {
  f <- data.frame(stratum = "A", N = 200, N_1 = 20, cost_1 = 30, cost_2 = 5)
  d <- data.frame(stratum = "A", name = "y", p = .3, icc_1 = .2)
  fit <- n_alloc_joint(f, d, objective = "y", budget = 900, fpc = "unit",
                       resp_rate = .8, take_2_min = 2, take_2_max = 10)
  opt <- fit$optimization
  expect_identical(opt$classification, "optimal")
  expect_lte(opt$objective_lower_bound,
             fit$params$achieved$objective_value * (1 + 1e-12))
  expect_lte(opt$objective_gap, opt$kkt_tolerance)
})

## The cheapest budget and the order of the precision bound

test_that("a budget that only the cheapest design fits returns that design", {
  f <- data.frame(stratum = "A", N = 100, N_1 = 10, cost_1 = 10, cost_2 = 1)
  m <- data.frame(stratum = "A", name = "y", p = .37, icc_1 = .31)
  fit <- n_alloc_joint(f, m, objective = "y", budget = 11)
  expect_identical(fit$optimization$classification, "feasible")
  expect_true(is.na(fit$optimization$objective_lower_bound))
  expect_equal(fit$params$achieved$cost, 11)
  expect_equal(fit$operational$cost, 11)
  expect_equal(c(fit$detail$n_1_int, fit$detail$take_2_int), c(1, 1))
  expect_error(n_alloc_joint(f, m, objective = "y", budget = 10.9),
               "cannot fund the precision targets")
})

test_that("the precision bound never exceeds the objective it certifies", {
  f <- data.frame(stratum = "A", N = 100, N_1 = 10, cost_1 = 10, cost_2 = 1)
  m <- data.frame(stratum = "A", name = "y", p = .37, icc_1 = .31)
  for (b in c(11.00001, 11.5, 150, 199)) {
    opt <- n_alloc_joint(f, m, objective = "y", budget = b)$optimization
    expect_gte(opt$objective_gap, 0)
    expect_identical(opt$classification, "optimal")
  }
  ## At budget 150 the optimum is a = 7.5 PSUs of take 10, from the
  ## closed form of one stratum and one target.
  fit <- n_alloc_joint(f, m, objective = "y", budget = 150)
  v <- function(a) (.63 / .37) * (.31 / a + .69 / (10 * a)) * (1 - 10 * a / 100)
  expect_lte(fit$optimization$objective_lower_bound, v(7.5) * (1 + 1e-12))
  expect_lte(fit$optimization$objective_lower_bound,
             fit$params$achieved$objective_value * (1 + 1e-12))
})

## One budget line for every whole search

test_that("a budget just below a design's cost does not buy it", {
  frames <- list(
    three = data.frame(stratum = "A", N = 100, N_1 = 10, N_2 = 20, mean = 10,
                       sd = 2, icc_1 = .1, icc_2 = .1, take_2 = 1, take_3 = 1,
                       cost_1 = 1, cost_2 = 1, cost_3 = 1),
    two = data.frame(stratum = "A", N = 100, N_1 = 10, mean = 10, sd = 2,
                     icc_1 = .1, take_2 = 2, cost_1 = 1, cost_2 = 1),
    element = data.frame(stratum = c("A", "B"), N = 100, mean = 10, sd = 2,
                         unit_cost = 1)
  )
  ## Each design of cost 6 is one cost step above the next cheaper one.
  below <- c(three = 3, two = 3, element = 5)
  cost_cols <- c("cost_1", "cost_2", "cost_3", "unit_cost")
  checked <- character(0)
  for (nm in names(frames)) {
    for (scale in c(1e-12, 1, 1e12)) {
      fr <- frames[[nm]]
      cols <- intersect(cost_cols, names(fr))
      fr[cols] <- lapply(fr[cols], function(v) v * scale)
      over <- suppressWarnings(n_alloc(fr, budget = 6 * scale * (1 - 1e-13),
                                       alloc = "optimal", fpc = "none"))
      expect_equal(over$operational$cost / scale, below[[nm]])
      fits <- suppressWarnings(n_alloc(fr, budget = 6 * scale,
                                       alloc = "optimal", fpc = "none"))
      expect_equal(fits$operational$cost / scale, 6)
    }
    checked <- c(checked, nm)
  }
  expect_setequal(checked, names(frames))
})

## Affordability at the accepted targets

test_that("a budget that funds a design within the target tolerance is not refused", {
  f <- data.frame(stratum = "A", N = 5000, N_1 = 50, cost_1 = 100, cost_2 = 1)
  m <- data.frame(stratum = "A", name = "y", p = .5, icc_1 = .5)
  ## 10 PSUs of 10 cost 1100 and reach cv^2 = .5 / 10 + .5 / 100 = .055,
  ## inside the relative 1e-6 tolerance of each target below.
  for (e in c(.9e-6, .99e-6)) {
    targets <- data.frame(name = "y", cv = sqrt(.055) / (1 + e))
    fit <- n_alloc_joint(f, m, targets = targets, objective = "y",
                         budget = 1100, fpc = "none")
    expect_lte(fit$params$achieved$cost, 1100 * (1 + 1e-6))
    expect_equal(fit$operational$cost, 1100)
    expect_true(all(prec_alloc_joint(fit)$detail$.pass))
    expect_error(n_alloc_joint(f, m, targets = targets, objective = "y",
                               budget = 1099.9, fpc = "none"),
                 "cannot fund the precision targets")
  }
})
