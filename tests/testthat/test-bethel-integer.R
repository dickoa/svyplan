## The joint allocation's whole design against every whole allocation

joint_frame <- function(seed) {
  set.seed(seed)
  H <- sample(2:3, 1)
  frame <- data.frame(stratum = LETTERS[seq_len(H)], N = sample(30:90, H),
                      unit_cost = round(stats::runif(H, 1, 10), 1))
  list(
    frame = frame,
    measures = data.frame(
      stratum = rep(frame$stratum, 2), name = rep(c("y", "x"), each = H),
      p = c(stats::runif(H, 0.2, 0.6), rep(NA, H)),
      mean = c(rep(NA, H), stats::runif(H, 10, 20)),
      sd = c(rep(NA, H), stats::runif(H, 3, 8))
    ),
    targets = data.frame(name = c("y", "x"),
                         cv = c(stats::runif(1, 0.08, 0.15),
                                stats::runif(1, 0.04, 0.08)))
  )
}

# Every whole decision of a fit's problem, with its cost, its objective
# and whether it meets every target and the budget at the tolerances the
# assessment accepts.
every_plan <- function(problem, budget = NULL) {
  grid <- as.matrix(expand.grid(lapply(
    seq_along(problem$cost),
    function(h) problem$lower_int[h]:problem$upper_int[h]
  )))
  accept <- 1 + .bethel_canonical_tolerance(problem, 1e-6)
  load <- (1 / grid) %*% problem$A
  ok <- rowSums(load <= rep(problem$bound * accept, each = nrow(grid))) ==
    ncol(problem$A)
  cost <- drop(grid %*% problem$cost)
  if (!is.null(budget)) ok <- ok & cost <= budget + 1e-6 * max(1, budget)
  objective <- if (!is.null(problem$objective)) {
    problem$objective$Q0 + drop((1 / grid) %*% problem$objective$q_h)
  }
  list(cost = cost, objective = objective, ok = ok)
}

test_that("a joint allocation is the cheapest whole allocation", {
  for (seed in 1:20) {
    z <- joint_frame(seed)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    plans <- every_plan(fit$params$problem)
    label <- paste("seed", seed)
    expect_equal(fit$operational$cost, min(plans$cost[plans$ok]),
                 tolerance = 1e-9, label = label)
    expect_identical(fit$operational$status, "optimal", label = label)
    expect_identical(fit$optimization$search$family, "whole units")
    expect_true(all(fit$operational$constraints$.pass), label = label)
  }
})

test_that("a joint budget buys the best whole allocation", {
  for (seed in 1:12) {
    z <- joint_frame(seed)
    base <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    for (hard in c(FALSE, TRUE)) {
      budget <- (if (hard) 1.2 else 0.8) * base$operational$cost
      fit <- n_alloc(z$frame, measures = z$measures,
                     targets = if (hard) z$targets,
                     objective = c("y", "x"), budget = budget)
      plans <- every_plan(fit$params$problem, budget)
      label <- paste("seed", seed, "hard", hard)
      expect_equal(fit$operational$objective_value,
                   min(plans$objective[plans$ok]), tolerance = 1e-9,
                   label = label)
      expect_identical(fit$operational$status, "optimal", label = label)
      expect_equal(fit$optimization$search$bound,
                   fit$operational$objective_value, tolerance = 1e-9)
      expect_lte(fit$operational$cost, budget * (1 + 1e-6))
      expect_equal(fit$operational$budget_residual,
                   budget - fit$operational$cost)
    }
  }
})

test_that("a two-stage joint allocation is the cheapest in whole PSUs", {
  for (seed in 1:8) {
    set.seed(100 + seed)
    H <- sample(2:3, 1)
    frame <- data.frame(stratum = LETTERS[seq_len(H)],
                        N = sample(2000:6000, H), N_psu = sample(15:40, H),
                        n_per_psu = sample(c(5, 8, 10), H, replace = TRUE),
                        cost_psu = sample(c(50, 120, 300), H, replace = TRUE),
                        cost_ssu = sample(c(5, 10, 20), H, replace = TRUE))
    measures <- data.frame(stratum = frame$stratum, name = "y",
                           p = stats::runif(H, 0.2, 0.6),
                           icc_psu = stats::runif(H, 0.02, 0.1))
    fit <- n_alloc(frame, measures = measures,
                   targets = data.frame(name = "y", cv = 0.08))
    plans <- every_plan(fit$params$problem)
    expect_equal(fit$operational$cost, min(plans$cost[plans$ok]),
                 tolerance = 1e-9, label = paste("seed", seed))
    expect_identical(fit$detail$n_psu_int, round(fit$detail$n_psu_int))
    expect_identical(fit$operational$status, "optimal")
  }
})

test_that("a search stopped at its node limit reports a valid bound", {
  z <- joint_frame(20)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 control = list(max_nodes = 1))
  best <- min(every_plan(fit$params$problem)$cost[
    every_plan(fit$params$problem)$ok])
  sr <- fit$optimization$search
  expect_identical(fit$operational$status, "limit")
  expect_identical(sr$max_nodes, 1)
  expect_gte(fit$operational$cost, best - 1e-9)
  expect_lte(sr$bound, best + 1e-9)
  expect_true(all(fit$operational$constraints$.pass))
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("^search: stopped at the node limit, gap ", out)))
  expect_true(any(grepl("whole-unit search: limit after",
                        utils::capture.output(summary(fit)))))
  expect_false(any(grepl("^search:", utils::capture.output(print(
    n_alloc(z$frame, measures = z$measures, targets = z$targets)
  )))))
  # The limit survives the round trip through the precision result.
  again <- n_alloc(prec_alloc(fit))
  expect_identical(again$optimization$search$max_nodes, 1)
  expect_identical(again$params$control, list(max_nodes = 1))
})

test_that("an unaffordable joint budget names the cheapest whole cost", {
  z <- joint_frame(5)
  base <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  # Between the continuous optimum's cost and the whole design's, the
  # budget funds the targets in the continuous model only.
  expect_lt(base$params$achieved$cost, base$operational$cost)
  err <- tryCatch(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            objective = "y", budget = base$params$achieved$cost),
    error = conditionMessage
  )
  expect_match(err, sprintf("costs %.6g against a budget",
                            base$operational$cost), fixed = TRUE)
  expect_match(err, "Raise 'budget' to at least", fixed = TRUE)
})

test_that("a target within tolerance below the census floor is met there", {
  fr <- data.frame(stratum = c("A", "B"), N = c(40, 60), unit_cost = c(1, 2))
  ms <- data.frame(stratum = c("A", "B"), name = "x", mean = 10, sd = c(3, 5))
  floor_cv <- prec_alloc(data.frame(N = c(40, 60), sd = c(3, 5), mean = 10),
                         n = c(40, 60), resp_rate = 0.8)$cv
  # The public tolerance is 1e-6 on the cv ratio.
  for (k in c(1e-12, 1e-9, 1e-7)) {
    fit <- n_alloc(fr, measures = ms,
                   targets = data.frame(name = "x", cv = floor_cv * (1 - k)),
                   resp_rate = 0.8)
    expect_identical(as.integer(fit$detail$n_int), c(40L, 60L), label = k)
    expect_identical(fit$operational$status, "optimal")
  }
  # Past it the gate refuses, rather than the whole-unit step failing on a
  # target the continuous solve admitted on its own scale.
  expect_error(
    n_alloc(fr, measures = ms,
            targets = data.frame(name = "x", cv = floor_cv * (1 - 2e-6)),
            resp_rate = 0.8),
    "precision targets are infeasible"
  )
})

test_that("the starting multipliers keep the top when the residual stays positive", {
  # The constraint at the upper bounds exceeds its bound by less than the
  # tolerance, so the bracket would hold two positive ends.
  A <- matrix(0.5, 2, 1)
  upper <- rep(1 - 2e-9, 2)
  expect_gt(sum(A / upper) - 1, 0)
  lambda <- .bethel_initial_lambda(A, c(1, 1), c(0.1, 0.1), upper)
  expect_true(is.finite(lambda) && lambda > 0)
})

test_that("a budget at the floor solves cleanly when a target parallels the objective", {
  # With one stratum, and with the objective the targeted indicator, the
  # target row is the objective's at another bound. Near the floor budget
  # the two coincided and the inner solve failed inside the root search.
  fr <- data.frame(stratum = "A", N = 100000)
  ms <- data.frame(stratum = "A", name = "y", p = 0.3)
  tg <- data.frame(name = "y", cv = 0.05)
  floor_cost <- n_alloc(fr, measures = ms, targets = tg,
                        eligible_rate = 0.1)$operational$cost
  for (b in c(floor_cost, floor_cost + 1)) {
    fit <- expect_no_warning(
      n_alloc(fr, measures = ms, targets = tg, objective = "y", budget = b,
              eligible_rate = 0.1)
    )
    expect_identical(as.integer(fit$detail$n_int), as.integer(b))
    expect_true(fit$operational$all_pass)
    opt <- fit$optimization
    expect_equal(opt$constraint_scale[1], fit$params$problem$bound[1])
    expect_length(opt$constraint_residual, 2L)
  }
  expect_equal(.bethel_parallel_scale(cbind(c(2, 4), c(1, 3)), c(1, 2)),
               c(2, NA))
})
