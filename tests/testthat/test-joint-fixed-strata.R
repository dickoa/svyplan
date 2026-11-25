# Strata fixed at their size, and small variances beside them

## A take-all stratum A of 50 units beside B and C, 1000 units each with SD
## 1 and 2. A census under unit FPC adds no variance, so the free units
## follow Neyman allocation whatever A's SD.
takeall_fixture <- function(sd_a) {
  list(
    frame = data.frame(stratum = c("A", "B", "C"), N = c(50, 1000, 1000)),
    measures = data.frame(stratum = c("A", "B", "C"), name = "y", mean = 1,
                          sd = c(sd_a, 1, 2))
  )
}

test_that("a take-all stratum leaves the budget allocation of the others exact", {
  ## Oracle: 200 free units split as N_h S_h, 1:2, and the relvar
  ## sum_h N_h^2 S_h^2 (1 / n_h - 1 / N_h) / 2050^2.
  n_free <- 200 * c(1, 2) / 3
  relvar <- sum(1000^2 * c(1, 4) * (1 / n_free - 1 / 1000)) / 2050^2
  for (sd_a in c(1e6, 1e7)) {
    z <- takeall_fixture(sd_a)
    for (rule in c("sum", "worst")) {
      fit <- n_alloc_joint(z$frame, z$measures, objective = "y",
                           budget = 250, min_n_stratum = c(50, 1, 1),
                           objective_rule = rule)
      label <- paste(rule, sd_a)
      expect_equal(fit$detail$n, c(50, n_free), tolerance = 1e-6,
                   label = label)
      expect_equal(fit$objective_value, relvar, tolerance = 1e-6,
                   label = label)
      expect_identical(fit$optimization$slopes, "available", label = label)
    }
  }
})

test_that("a take-all stratum leaves the cost-minimum allocation of the others exact", {
  ## Oracle: minimizing n_B + n_C subject to
  ## sum_h N_h^2 S_h^2 / n_h - sum_h N_h S_h^2 <= V gives
  ## n_h = N_h S_h sum(N S) / (V + sum N S^2), with V = (0.05 * 2050)^2.
  NS <- 1000 * c(1, 2)
  n_free <- NS * sum(NS) / ((0.05 * 2050)^2 + sum(1000 * c(1, 4)))
  for (sd_a in c(1e6, 1e7)) {
    z <- takeall_fixture(sd_a)
    fit <- n_alloc_joint(z$frame, z$measures,
                         targets = data.frame(name = "y", cv = 0.05),
                         min_n_stratum = c(50, 1, 1))
    expect_equal(fit$detail$n, c(50, n_free), tolerance = 1e-6,
                 label = format(sd_a))
    expect_true(all(fit$operational$constraints$.pass), label = format(sd_a))
  }
})

## Free takes: stratum C is every one of its 10 PSUs and every unit in them,
## so it adds no variance and its SD cannot move the design. Its floor of 40
## units at a largest take of 4 needs all 10 PSUs, so a smallest take of 1
## leaves the same single design as a take fixed at 4.
freetake_takeall <- function(sd_c, mode, c_min = 4) {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(180, 240, 40),
                      N_1 = c(60, 60, 10), cost_1 = c(12, 8, 5),
                      cost_2 = c(2, 3, 1))
  measures <- data.frame(stratum = c("A", "B", "C"), name = "inc",
                         mean = 100, sd = c(25, 30, sd_c),
                         icc_1 = c(0.04, 0.42, 0.1))
  args <- list(frame, measures, fpc = "unit", take_2_min = c(1, 1, c_min),
               take_2_max = c(3, 4, 4), min_n_stratum = c(1, 1, 40))
  if (identical(mode, "budget")) {
    do.call(n_alloc_joint, c(args, list(objective = "inc", budget = 1000)))
  } else {
    do.call(n_alloc_joint,
            c(args, list(targets = data.frame(name = "inc", cv = 0.004))))
  }
}

test_that("a free-take census stratum does not move the design of the others", {
  ## Oracle: the fit at SD 1e3 with the take fixed at 4, where nothing
  ## cancels badly.
  for (mode in c("budget", "target")) {
    base <- freetake_takeall(1e3, mode)
    expect_identical(base$optimization$classification, "optimal")
    for (case in list(c(1e6, 4), c(1e7, 4), c(1e6, 1), c(1e7, 1))) {
      sd_c <- case[1]
      label <- paste(mode, sd_c, "smallest take", case[2])
      fit <- expect_no_warning(freetake_takeall(sd_c, mode, c_min = case[2]))
      expect_identical(fit$optimization$classification, "optimal",
                       label = label)
      expect_equal(fit$detail$n_1, base$detail$n_1, tolerance = 1e-6,
                   label = label)
      expect_equal(fit$detail$take_2, base$detail$take_2, tolerance = 1e-6,
                   label = label)
      expect_identical(fit$detail$n_1_int, base$detail$n_1_int,
                       label = label)
      expect_identical(fit$detail$take_2_int, base$detail$take_2_int,
                       label = label)
    }
  }
})

test_that("a PSU census keeps the within-PSU variance however small", {
  ## Every PSU is taken, 5 of the 20 units in each, ICC 1 - delta with
  ## unit variance. The between-PSU term vanishes, and the total's
  ## variance is 10 * 20^2 (1 - 5 / 20) delta / 5 = 600 delta, so the
  ## relvar is 600 delta / 200^2 = 3 delta / 200.
  delta <- 2^-48
  measures <- data.frame(stratum = "A", name = "y", mean = 1, sd = 1,
                         icc_1 = 1 - delta)
  frame <- data.frame(stratum = "A", N = 200, N_1 = 10, cost_1 = 5,
                      cost_2 = 1, take_2 = 5)
  for (rule in c("sum", "worst")) {
    assessed <- prec_alloc_joint(frame, n = 50, measures = measures,
                                 objective = "y", fpc = "stage",
                                 objective_rule = rule)
    expect_equal(assessed$objective_value / (3 * delta / 200), 1,
                 tolerance = 1e-10, label = paste("assessed", rule))
    fit <- n_alloc_joint(frame[setdiff(names(frame), "take_2")], measures,
                         objective = "y", budget = 100, fpc = "stage",
                         take_2_min = 5, take_2_max = 5,
                         objective_rule = rule)
    expect_equal(fit$objective_value / (3 * delta / 200), 1,
                 tolerance = 1e-10, label = paste("fitted", rule))
  }
})

test_that("a census beside a low-variance stratum keeps the target's small variance", {
  ## Oracle: the relvar of the sampled stratum,
  ## W_B^2 S_B^2 (1 / 10 - 1 / 100), set to 2.25e-17, so the CV is
  ## sqrt(2.25e-17). The census stratum contributes nothing.
  sd_b <- sqrt(2.25e-17 / ((100 / 110)^2 * (1 / 10 - 1 / 100)))
  frame <- data.frame(stratum = c("A", "B"), N = c(10, 100))
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 1,
                         sd = c(1000, sd_b))
  assessed <- prec_alloc_joint(frame, n = c(10, 10), measures = measures,
                               targets = data.frame(name = "y", cv = 0.1),
                               fpc = "unit")
  expect_equal(unname(assessed$cv) / sqrt(2.25e-17), 1, tolerance = 1e-9)
})

test_that("a census beside a low-variance stratum still spends the budget", {
  ## Stratum A is fixed at its census, stratum B's relvar falls as n_B
  ## grows, and the budget of 110 buys n_B = 10 at unit cost.
  frame <- data.frame(stratum = c("A", "B"), N = 100)
  measures <- data.frame(stratum = c("A", "B"), name = "a", mean = 1,
                         sd = c(1, sqrt(1e-15)))
  for (rule in c("sum", "worst")) {
    fit <- n_alloc_joint(frame, measures, objective = "a", budget = 110,
                         min_n_stratum = c(100, 1), objective_rule = rule)
    expect_equal(fit$detail$n, c(100, 10), tolerance = 1e-6, label = rule)
    expect_identical(fit$optimization$slopes, "available", label = rule)
  }
})

test_that("the share a correction leaves is exact one unit short of a census", {
  ## Oracle: (P - n) / P is 1 / P exactly when n = P - 1, which 1 - n / P
  ## gives only to about 1e-7 relative at P = 1e9. The free-take between-PSU
  ## term reads this share at the PSU count itself.
  P <- 1e9
  expect_equal(.fpc_remaining(P - 1, P, 1) * P, 1, tolerance = 1e-12)
  expect_equal(.fpc_remaining(P, P, 0.8), 0.2, tolerance = 1e-15)
  expect_identical(.fpc_remaining(P, P, 1), 0)
})

test_that("an assessment below a stratum's floor reads the design given", {
  ## min_n_stratum = 100 fixes stratum A at its census for the solvers, but
  ## the design assessed takes 50 of its 100 units. Oracle: the relvar
  ## sum_h W_h^2 S^2 (1 / n_h - 1 / N_h) with W_h = 1/2, S = 1 and
  ## n_h = 50, which is 0.005, CV 0.0707, above the target of 0.06.
  frame <- data.frame(stratum = c("A", "B"), N = 100)
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 1, sd = 1)
  relvar <- 2 * 0.25 * (1 / 50 - 1 / 100)
  assessed <- prec_alloc_joint(frame, n = c(50, 50), measures = measures,
                               targets = data.frame(name = "y", cv = 0.06),
                               min_n_stratum = c(100, 1), fpc = "unit")
  expect_equal(unname(assessed$cv), sqrt(relvar), tolerance = 1e-12)
  expect_false(assessed$detail$.pass)
  expect_identical(assessed$bounds$.lower_violation, c(TRUE, FALSE))
  for (rule in c("sum", "worst")) {
    assessed <- prec_alloc_joint(frame, n = c(50, 50), measures = measures,
                                 objective = "y", min_n_stratum = c(100, 1),
                                 fpc = "unit", objective_rule = rule)
    expect_equal(assessed$objective_value, relvar, tolerance = 1e-12,
                 label = rule)
  }
})

test_that("a free-take assessment below a stratum's floor reads the design given", {
  ## Take 5 and min_n_stratum = 100 fix stratum A at all 20 PSUs for the
  ## solvers, but the design assessed takes 10. Oracle under fpc = "none":
  ## sum_h W_h^2 S^2 (1 + (m - 1) icc) / n_h with W_h = 1/2, m = 5,
  ## icc = 0.1 and n_h = 50, which is 0.014.
  frame <- data.frame(stratum = c("A", "B"), N = 100, N_1 = 20, cost_1 = 5,
                      cost_2 = 1)
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 1,
                         sd = 1, icc_1 = 0.1)
  assessed <- .freetake_assess(
    frame, measures, targets = NULL, objective = "y", n_1 = c(10, 10),
    take_2 = c(5, 5), alpha = 0.05, deff = 1, resp_rate = 1,
    eligible_per_unit = 1, eligibility_rate = 1, df = NULL,
    min_n_stratum = c(100, 1), fpc = "none"
  )
  expect_equal(assessed$objective$components$.relvar,
               2 * 0.25 * (1 + 4 * 0.1) / 50, tolerance = 1e-12)
})
