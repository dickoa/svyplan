test_that("generalized allocation is invariant to table permutations", {
  for (stages in 1:3) {
    z <- if (stages == 1L) .bethel_fixture() else
      .bethel_multistage_fixture(stages)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    permuted <- n_alloc(
      z$frame[nrow(z$frame):1, , drop = FALSE],
      measures = z$measures[nrow(z$measures):1, , drop = FALSE],
      targets = z$targets[nrow(z$targets):1, , drop = FALSE]
    )

    x <- setNames(fit$detail$n, fit$detail$stratum)
    y <- setNames(permuted$detail$n, permuted$detail$stratum)
    expect_equal(x[sort(names(x))], y[sort(names(y))], tolerance = 1e-5)
    xc <- fit$constraints[order(fit$constraints$constraint), ]
    yc <- permuted$constraints[order(permuted$constraints$constraint), ]
    expect_equal(xc$.achieved, yc$.achieved, tolerance = 1e-7)
  }
})

test_that("rescaling every field cost preserves the continuous allocation", {
  scale <- 17
  for (stages in 1:3) {
    z <- if (stages == 1L) .bethel_fixture() else
      .bethel_multistage_fixture(stages)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    scaled_frame <- z$frame
    if (stages == 1L) {
      scaled_frame$unit_cost <- scaled_frame$unit_cost * scale
    } else {
      cost_names <- c("cost_psu", "cost_ssu", "cost_tsu")
      for (nm in intersect(cost_names, names(scaled_frame))) {
        scaled_frame[[nm]] <- scaled_frame[[nm]] * scale
      }
    }
    scaled <- n_alloc(
      scaled_frame, measures = z$measures, targets = z$targets
    )

    expect_equal(scaled$detail$n, fit$detail$n, tolerance = 1e-5)
    expect_equal(
      scaled$params$achieved$cost,
      scale * fit$params$achieved$cost,
      tolerance = 1e-6
    )
  }
})

test_that("tightening all targets cannot reduce minimum continuous cost", {
  for (stages in 1:3) {
    z <- if (stages == 1L) .bethel_fixture() else
      .bethel_multistage_fixture(stages)
    loose <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    tight_targets <- z$targets
    if ("cv" %in% names(tight_targets)) {
      tight_targets$cv <- tight_targets$cv * 0.9
    }
    if ("moe" %in% names(tight_targets)) {
      tight_targets$moe <- tight_targets$moe * 0.9
    }
    tight <- n_alloc(
      z$frame, measures = z$measures, targets = tight_targets
    )
    expect_gte(
      tight$params$achieved$cost,
      loose$params$achieved$cost * (1 - 1e-7)
    )
  }
})

test_that("joint result schemas are stable across one to three stages", {
  constraint_names <- c(
    "constraint", "name", "domain", "level", ".metric", ".target",
    ".achieved", ".ratio", ".residual", ".tolerance", ".se", ".cv",
    ".moe", ".pass", ".binding", ".multiplier", ".sensitivity"
  )
  bounds_names <- c(
    "stratum", "n", ".lower", ".upper", ".lower_violation",
    ".upper_violation", ".pass"
  )
  for (stages in 1:3) {
    z <- if (stages == 1L) .bethel_fixture() else
      .bethel_multistage_fixture(stages)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    assessed <- prec_alloc(fit)

    expect_identical(names(fit$constraints), constraint_names)
    expect_identical(names(fit$operational), c(
      "n", "cost", "constraints", "all_pass", "repair_iterations"
    ))
    expect_identical(names(assessed$bounds), bounds_names)
    expect_identical(fit$params$stages, stages)
    expect_identical(assessed$params$stages, stages)
    expect_identical(
      unique(fit$constraints$.tolerance),
      fit$params$feasibility_tolerance
    )
    expect_true(is.double(fit$detail$n))
    expect_true(is.double(fit$detail$n_int))
    expect_true(is.logical(fit$constraints$.pass))
    expect_true(is.logical(assessed$bounds$.pass))
  }
})

test_that("integer repair is feasible and near the exhaustive tiny optimum", {
  set.seed(919)
  for (case in seq_len(40)) {
    H <- sample(2:4, 1)
    K <- sample(1:3, 1)
    A <- matrix(runif(H * K, 0.2, 5), H, K)
    lower <- rep(1, H)
    upper <- sample(3:7, H, replace = TRUE)
    grid <- expand.grid(lapply(upper, seq_len))
    candidates <- as.matrix(grid)
    witness <- candidates[sample(nrow(candidates), 1), ]
    bound <- drop(crossprod(A, 1 / witness)) * runif(K, 1, 1.3)
    cost <- runif(H, 0.5, 3)
    feasible <- apply(candidates, 1, function(n) {
      all(drop(crossprod(A, 1 / n)) <= bound * (1 + 1e-6))
    })
    expect_true(any(feasible), info = paste("case", case))

    ids <- paste0("c", seq_len(K))
    meta <- data.frame(
      constraint = ids, name = ids, domain = ".overall",
      level = NA_character_, .metric = "moe",
      .target = stats::qnorm(0.975) * sqrt(bound)
    )
    problem <- .new_bethel_problem(
      A, rep(0, K), bound, cost, lower, upper,
      paste0("h", seq_len(H)), ids, meta,
      total = rep(1, K), domain_N = rep(1, K), alpha = rep(0.05, K)
    )
    solved <- .bethel_solve(
      A, bound, cost, lower, upper, constraint_ids = ids
    )
    expect_identical(solved$classification, "optimal", info = paste("case", case))
    integer <- .integerize_bethel(solved$allocation, problem)
    exact_cost <- min(drop(candidates[feasible, , drop = FALSE] %*% cost))

    expect_true(integer$all_pass, info = paste("case", case))
    # Integer cleanup is deliberately local rather than an exact
    # combinatorial optimizer. This deterministic audit set caps its observed
    # cost gap while the exact optimum remains an independent feasibility oracle.
    expect_true(
      integer$cost <= 1.20 * exact_cost,
      info = paste("case", case)
    )
  }
})
