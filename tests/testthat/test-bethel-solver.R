test_that("Bethel solver reproduces the one-constraint analytic solution", {
  A <- matrix(c(4, 9, 16), ncol = 1)
  cost <- c(1, 2, 4)
  bound <- 1
  expected <- sqrt(A[, 1] / cost) *
    sum(sqrt(A[, 1] * cost)) / bound

  out <- .bethel_solve(
    A, bound, cost,
    lower = rep(0.1, 3), upper = rep(100, 3)
  )

  expect_identical(out$classification, "optimal")
  expect_equal(out$allocation, expected, tolerance = 1e-5)
  expect_lte(out$primal_residual, 1e-6)
  expect_lte(out$projected_dual_residual, 1e-6)
  expect_lte(out$stationarity_residual, 1e-6)
  expect_lte(out$relative_duality_gap, 1e-6)
  expect_identical(out$feasibility_tolerance, .bethel_control()$tolerance)
  expect_lte(out$kkt_tolerance, out$feasibility_tolerance)
})

test_that("solver certification and public feasibility share one contract", {
  control <- .bethel_control()
  expect_lte(control$kkt_tolerance, control$tolerance)
  expect_error(
    .bethel_control(tolerance = 1e-8, kkt_tolerance = 1e-6),
    "must not exceed"
  )

  problem <- .new_bethel_problem(
    A = matrix(1, 1, 1), B = 0, bound = 1,
    cost = 1, lower = 0.5, upper = 2,
    stratum_ids = "a", constraint_ids = "target",
    constraint_meta = data.frame(
      constraint = "target", name = "target", domain = ".overall",
      level = NA_character_, .metric = "cv", .target = 1
    ),
    total = 1, domain_N = 1, alpha = 0.05
  )
  # A ratio excess inside the KKT/public feasibility contract must not be
  # reported as a failure by the public precision evaluator.
  allocation <- 1 / (1 + control$tolerance / 2)^2
  evaluated <- .precision_from_allocation(problem, allocation)
  expect_true(evaluated$all_pass)
  expect_equal(evaluated$constraints$.tolerance, control$tolerance)
  expect_gt(evaluated$constraints$.residual, 0)
})

test_that("Bethel solver handles lower, upper, and infeasible shortcuts", {
  lower <- .bethel_solve(matrix(1, 1, 1), 1, 1, 1, 10)
  upper <- .bethel_solve(matrix(1, 1, 1), 0.1, 1, 1, 10)
  bad <- .bethel_solve(
    matrix(1, 1, 1), 0.09, 1, 1, 10,
    constraint_ids = "hard"
  )

  expect_equal(lower$allocation, 1)
  expect_equal(lower$iterations, 0L)
  expect_equal(upper$allocation, 10)
  expect_identical(upper$classification, "optimal")
  expect_identical(bad$classification, "infeasible")
  expect_identical(bad$infeasible, "hard")
})

test_that("Bethel scaling and permutation preserve the solution", {
  A <- matrix(c(4, 1, 1, 9, 3, 2), nrow = 3)
  bound <- c(1, 1.2)
  cost <- c(1, 2, 1.5)
  lower <- rep(1, 3)
  upper <- rep(30, 3)
  x <- .bethel_solve(A, bound, cost, lower, upper)
  y <- .bethel_solve(
    A[, 2:1, drop = FALSE] * rep(c(10, 0.2), each = 3),
    bound[2:1] * c(10, 0.2),
    cost * 7, lower, upper
  )

  expect_identical(x$classification, "optimal")
  expect_identical(y$classification, "optimal")
  expect_equal(x$allocation, y$allocation, tolerance = 1e-5)
})

test_that("Bethel solver certifies deterministic random feasible problems", {
  set.seed(412)
  for (i in seq_len(20)) {
    H <- sample(3:8, 1)
    K <- sample(1:4, 1)
    A <- matrix(rexp(H * K), H, K)
    cost <- runif(H, 0.5, 3)
    lower <- runif(H, 0.5, 2)
    upper <- lower + runif(H, 5, 20)
    witness <- lower + runif(H, 0.2, 0.8) * (upper - lower)
    bound <- drop(crossprod(A, 1 / witness)) * runif(K, 1, 1.2)

    out <- .bethel_solve(A, bound, cost, lower, upper)
    expect_identical(out$classification, "optimal")
    expect_lte(out$primal_residual, out$feasibility_tolerance)
    expect_true(all(drop(crossprod(A, 1 / out$allocation)) <=
                    bound * (1 + 1e-6)))
  }
})

test_that("dual multiplier matches the finite-difference cost sensitivity", {
  A <- matrix(c(4, 9, 16), ncol = 1)
  bound <- 1
  args <- list(
    A = A, cost = c(1, 2, 4),
    lower = rep(0.1, 3), upper = rep(100, 3)
  )
  fit <- do.call(.bethel_solve, c(args, list(bound = bound)))
  step <- 1e-5
  plus <- do.call(.bethel_solve, c(args, list(bound = bound + step)))
  minus <- do.call(.bethel_solve, c(args, list(bound = bound - step)))
  derivative <- (plus$cost - minus$cost) / (2 * step)

  expect_equal(derivative, -fit$lambda, tolerance = 1e-4)
})

test_that("duplicate active constraints do not invent multiplier splits", {
  A <- cbind(c(4, 9, 16), c(8, 18, 32))
  fit <- .bethel_solve(
    A, bound = c(1, 2), cost = c(1, 2, 4),
    lower = rep(0.1, 3), upper = rep(100, 3)
  )

  expect_identical(fit$classification, "optimal")
  expect_true(all(!fit$multiplier_identifiable))
  expect_true(all(is.na(fit$lambda)))
})

test_that("integer cleanup leaves no feasible one-unit deletion", {
  A <- matrix(c(20, 15, 8, 10, 5, 18), nrow = 3)
  problem <- .new_bethel_problem(
    A = A, B = c(0, 0), bound = c(1.5, 1.2),
    cost = c(1, 1.5, 2), lower = rep(1, 3), upper = rep(50, 3),
    stratum_ids = c("a", "b", "c"), constraint_ids = c("x", "y"),
    constraint_meta = data.frame(
      constraint = c("x", "y"), name = c("x", "y"),
      domain = ".overall", level = NA_character_,
      .metric = "moe", .target = qnorm(0.975) * sqrt(c(1.5, 1.2))
    ),
    total = c(1, 1), domain_N = c(1, 1), alpha = c(0.05, 0.05)
  )
  solved <- .bethel_solve(A, problem$bound, problem$cost,
                          problem$lower, problem$upper)
  integer <- .integerize_bethel(solved$allocation, problem)
  normalized <- sweep(A, 2, problem$bound, "/")

  for (h in which(integer$allocation > problem$lower_int)) {
    candidate <- integer$allocation
    candidate[h] <- candidate[h] - 1L
    expect_true(any(.bethel_constraint_value(normalized, candidate) > 1 + 1e-8))
  }
})
