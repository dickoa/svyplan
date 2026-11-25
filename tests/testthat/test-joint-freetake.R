test_that("a small relative variance changes the scale of the budget objective only", {
  frame <- data.frame(stratum = "A", N = 10000, N_1 = 100,
                      cost_1 = 1000, cost_2 = 100)
  measures <- data.frame(stratum = "A", name = "y", mean = 1, sd = 1,
                         icc_1 = 0.05)
  for (fpc in c("none", "stage", "unit")) {
    solve <- function(scale) {
      measures$mean <- scale
      n_alloc_joint(frame, measures, targets = data.frame(name = "y", cv = 1 / scale),
                    budget = 5500, objective = "y", fpc = fpc)
    }
    reference <- solve(1)
    fit <- solve(1e9)
    expect_equal(fit$params$achieved$cost, 5500, tolerance = 1e-8)
    expect_equal(fit$objective_value * 1e18, reference$objective_value,
                 tolerance = 1e-8)
    expect_lte(fit$objective_value / fit$operational$objective_value, 1 + 1e-8)
    expect_identical(fit$optimization$classification, "optimal")
  }
})

test_that("the free-take objective minimum uses cost only on objective ties", {
  ## Large objectives at the left edge cannot make the two right-hand
  ## corners equivalent. With A = B = D = 0 the lowest-cost corner wins.
  pt <- .freetake_stratum_obj_min(0.1, 1, 1, 0.01, 0, 1e13, 1e13, 1, 1, 10)
  expect_equal(pt, list(x = 1e13, y = 1e13))
  flat <- .freetake_stratum_obj_min(2, 3, 0, 0, 0, 10, 100, 12, 2, 5)
  expect_equal(flat, list(x = 12 / 5, y = 12))
  ## With A = D * N the objective is flat along y = N, so the bend at
  ## N / m_max ties the corner at N / m_min and costs less.
  bend <- .freetake_stratum_obj_min(1, 1, 16, 1, 0.25, 32, 64, 1, 1, 4)
  value <- function(x) 16 / x + 1 / 64 - 0.25 * 64 / x
  expect_identical(value(16), value(32))
  expect_equal(bend, list(x = 16, y = 64))
})

test_that("a nonbinding budget buys only the required take outside the objective domain", {
  frame <- data.frame(stratum = c("A", "B"), region = c("a", "b"),
                      N = 100, N_1 = 10, cost_1 = 10, cost_2 = 1)
  measures <- data.frame(stratum = c("A", "B"), name = "y", mean = 1, sd = 1,
                         icc_1 = 0.05)
  objective <- data.frame(name = "y", domain = "region", level = "a")
  fit <- n_alloc_joint(frame, measures, objective = objective, budget = 300,
                       fpc = "none")
  expect_equal(fit$detail$n_1, c(10, 1))
  expect_equal(fit$detail$take_2, c(10, 1))
  expect_identical(fit$optimization$classification, "optimal")
  expect_false(fit$optimization$budget_binding)
})

test_that("free-take 2-stage n_alloc_joint matches n_cluster on 1-stratum 1-indicator", {
  f <- data.frame(
    stratum = "A", N = 50000, N_1 = 500,
    cost_1 = 200, cost_2 = 15, take_2 = NA_real_
  )
  m <- data.frame(stratum = "A", name = "y", mean = 100, sd = 20, icc_1 = 0.05)
  t <- data.frame(name = "y", cv = 0.03)

  # fpc = "none"
  j_none <- n_alloc_joint(f, m, targets = t, fpc = "none")
  c_none <- n_cluster(
    stage_cost = c(200, 15), icc = 0.05, cv = 0.03,
    unit_relvar = (20 / 100)^2
  )
  expect_equal(j_none$detail$n_1, unname(c_none$n[1]), tolerance = 1e-5)
  expect_equal(j_none$detail$take_2, unname(c_none$n[2]), tolerance = 1e-5)
  expect_equal(j_none$params$achieved$cost, c_none$cost, tolerance = 1e-5)
  expect_equal(j_none$operational$status, "optimal")
  expect_true(j_none$operational$all_pass)

  # fpc = "stage"
  j_stage <- n_alloc_joint(f, m, targets = t, fpc = "stage")
  c_stage <- n_cluster(
    stage_cost = c(200, 15), icc = 0.05, cv = 0.03,
    unit_relvar = (20 / 100)^2, N = c(500, 50000), fpc = "stage"
  )
  expect_equal(j_stage$detail$n_1, unname(c_stage$n[1]), tolerance = 1e-5)
  expect_equal(j_stage$detail$take_2, unname(c_stage$n[2]), tolerance = 1e-5)
  expect_equal(j_stage$params$achieved$cost, c_stage$cost, tolerance = 1e-5)
  expect_equal(j_stage$operational$status, "optimal")
})

test_that("free-take 2-stage integer optimum matches brute-force enumeration across FPC modes", {
  frame <- data.frame(
    stratum = c("U", "R"),
    region = c("North", "South"),
    N = c(1200, 1600),
    N_1 = c(30, 40),
    cost_1 = c(100, 180),
    cost_2 = c(15, 10)
  )
  measures <- data.frame(
    stratum = rep(c("U", "R"), 2),
    name = rep(c("prev", "inc"), each = 2),
    p = c(0.4, 0.3, NA, NA),
    mean = c(NA, NA, 50, 60),
    sd = c(NA, NA, 12, 18),
    icc_1 = c(0.08, 0.03, 0.02, 0.09)
  )
  targets <- data.frame(
    name = c("prev", "inc"),
    cv = c(0.07, 0.045)
  )

  # Brute-force reference over all integer (take_U, take_R) in [4, 20]
  brute_min_cost <- function(fpc_mode, t_vals_U, t_vals_R) {
    best_cost <- Inf
    best_pair <- NULL
    for (mU in t_vals_U) {
      for (mR in t_vals_R) {
        f_fix <- frame
        f_fix$take_2 <- c(mU, mR)
        res <- tryCatch(
          n_alloc_joint(f_fix, measures, targets = targets, fpc = fpc_mode),
          error = function(e) NULL
        )
        if (!is.null(res) && res$operational$all_pass &&
            res$operational$cost < best_cost) {
          best_cost <- res$operational$cost
          best_pair <- c(mU, mR)
        }
      }
    }
    list(cost = best_cost, take = best_pair)
  }

  for (fpc_mode in c("unit", "stage", "none")) {
    fit <- n_alloc_joint(
      frame, measures, targets = targets,
      take_2_min = 4, take_2_max = 20, take_2_step = 2, fpc = fpc_mode
    )
    bf <- brute_min_cost(fpc_mode, seq(4, 20, by = 2), seq(4, 20, by = 2))
    expect_equal(fit$operational$status, "optimal")
    expect_true(fit$operational$all_pass)
    expect_equal(fit$operational$cost, bf$cost, tolerance = 1e-8)
    expect_lte(fit$params$achieved$cost, fit$operational$cost + 1e-6)
    expect_true(all((fit$detail$take_2_int - 4) %% 2 == 0))

    # Round-trip prec_alloc_joint for both continuous and integer
    pc <- prec_alloc_joint(fit)
    pi <- prec_alloc_joint(fit, n = fit$detail$n_int)
    expect_equal(pc$detail$.achieved, fit$constraints$.achieved, tolerance = 1e-10)
    expect_equal(pi$detail$.achieved, fit$operational$constraints$.achieved,
                 tolerance = 1e-10)
  }
})

test_that("free-take 2-stage handles small ICC under fpc = 'stage' (negative A_hk)", {
  # N / N_1 = 40, so icc_1 = 0.005 < 1 / 41 gives R1 - R2 / M_bar < 0
  frame <- data.frame(
    stratum = c("A", "B"),
    N = c(2000, 2000),
    N_1 = c(50, 50),
    cost_1 = c(120, 120),
    cost_2 = c(10, 10)
  )
  measures <- data.frame(
    stratum = c("A", "B"),
    name = c("y", "y"),
    mean = c(50, 60),
    sd = c(15, 18),
    icc_1 = c(0.005, 0.005)
  )
  targets <- data.frame(name = "y", cv = 0.01)

  fit <- n_alloc_joint(frame, measures, targets = targets,
                       take_2_min = 5, take_2_max = 40, fpc = "stage")
  expect_equal(fit$operational$status, "optimal")
  expect_true(fit$operational$all_pass)
  # Continuous optimum hits take = 40 in both strata; integer optimum can trim
  # take_2 slightly after rounding n_1 up and must beat fixed take_2 = 40.
  expect_equal(fit$detail$take_2, c(40, 40), tolerance = 1e-5)
  f40 <- frame
  f40$take_2 <- c(40, 40)
  fix40 <- n_alloc_joint(f40, measures, targets = targets, fpc = "stage")
  expect_lte(fit$operational$cost, fix40$operational$cost)
})

test_that("free-take 2-stage supports mixed fixed/free strata and frame bound columns", {
  frame <- data.frame(
    stratum = c("U", "R"),
    N = c(2000, 3000),
    N_1 = c(40, 60),
    cost_1 = c(120, 200),
    cost_2 = c(15, 12),
    take_2 = c(12, NA_real_),
    take_2_min = c(NA_real_, 6),
    take_2_max = c(NA_real_, 24),
    take_2_step = c(NA_real_, 3)
  )
  measures <- data.frame(
    stratum = rep(c("U", "R"), 2),
    name = rep(c("a", "b"), each = 2),
    p = c(0.4, 0.35, NA, NA),
    mean = c(NA, NA, 80, 90),
    sd = c(NA, NA, 20, 25),
    icc_1 = c(0.06, 0.04, 0.03, 0.08)
  )
  targets <- data.frame(name = c("a", "b"), cv = c(0.06, 0.04))

  fit <- n_alloc_joint(frame, measures, targets = targets)
  expect_equal(fit$operational$status, "optimal")
  expect_equal(fit$detail$take_2[1], 12)
  expect_equal(fit$detail$take_2_int[1], 12)
  expect_true(fit$detail$take_2_int[2] %in% seq(6, 24, by = 3))
  expect_true(fit$operational$all_pass)
})

test_that("free-take 2-stage budget-objective mode matches brute-force enumeration", {
  frame <- data.frame(
    stratum = c("U", "R"),
    N = c(1600, 2000),
    N_1 = c(40, 50),
    cost_1 = c(110, 190),
    cost_2 = c(14, 10)
  )
  measures <- data.frame(
    stratum = rep(c("U", "R"), 2),
    name = rep(c("p1", "y1"), each = 2),
    p = c(0.35, 0.45, NA, NA),
    mean = c(NA, NA, 100, 120),
    sd = c(NA, NA, 25, 30),
    icc_1 = c(0.07, 0.04, 0.03, 0.08)
  )
  targets <- data.frame(name = "p1", cv = 0.08)
  budget <- 8000

  fit <- n_alloc_joint(
    frame, measures, targets = targets,
    objective = "y1", budget = budget,
    take_2_min = 6, take_2_max = 18, take_2_step = 3
  )
  expect_equal(fit$operational$status, "optimal")
  expect_true(fit$operational$all_pass)
  expect_lte(fit$operational$cost, budget + 1e-6)
  expect_lte(fit$objective_value, fit$operational$objective_value + 1e-8)

  # Brute-force check over all (take_U, take_R) in {6, 9, 12, 15, 18}
  best_obj <- Inf
  for (mU in seq(6, 18, by = 3)) {
    for (mR in seq(6, 18, by = 3)) {
      f_fix <- frame
      f_fix$take_2 <- c(mU, mR)
      res <- tryCatch(
        n_alloc_joint(f_fix, measures, targets = targets,
                      objective = "y1", budget = budget),
        error = function(e) NULL
      )
      if (!is.null(res) && res$operational$all_pass &&
          res$operational$cost <= budget + 1e-6 &&
          res$operational$objective_value < best_obj) {
        best_obj <- res$operational$objective_value
      }
    }
  }
  expect_equal(fit$operational$objective_value, best_obj, tolerance = 1e-10)
})

test_that("free-take 2-stage with PSU register matches brute-force over fixed-take runs", {
  psu <- data.frame(
    psu_id = sprintf("P%02d", 1:20),
    stratum = rep(c("U", "R"), each = 10),
    N = c(c(400, rep(100, 9)), c(350, rep(90, 9)))
  )
  frame <- data.frame(
    stratum = c("U", "R"),
    N = c(sum(psu$N[psu$stratum == "U"]), sum(psu$N[psu$stratum == "R"])),
    cost_1 = c(180, 240),
    cost_2 = c(12, 10)
  )
  measures <- data.frame(
    stratum = rep(c("U", "R"), 2),
    name = rep(c("p1", "y1"), each = 2),
    p = c(0.4, 0.35, NA, NA),
    mean = c(NA, NA, 50, 60),
    sd = c(NA, NA, 12, 15),
    icc_1 = c(0.06, 0.04, 0.03, 0.07)
  )
  targets <- data.frame(name = c("p1", "y1"), cv = c(0.08, 0.055))

  # Target mode with certainty_cutoff + n_psu_per_zone = 2
  fit <- n_alloc_joint(
    frame, measures, targets = targets, psu = psu,
    certainty_cutoff = 0.85, n_psu_per_zone = 2,
    take_2_min = 6, take_2_max = 18, take_2_step = 4
  )
  expect_equal(fit$operational$status, "optimal")
  expect_true(fit$operational$all_pass)
  expect_true(all(c("take_2", "take_2_int") %in% names(fit$detail)))
  expect_true(all(fit$detail$take_2_int %in% c(6, 10, 14, 18)))

  # Non-certainty PSUs in $psu$take_2 match optimal take_2_int
  for (h in seq_len(nrow(frame))) {
    st <- frame$stratum[h]
    rem <- fit$psu$stratum == st & !fit$psu$certainty
    if (any(rem)) {
      expect_true(all(fit$psu$take_2[rem] == fit$detail$take_2_int[h]))
    }
  }

  # Brute-force check against fixed-take PSU runs across every (mU, mR)
  best_cost <- Inf
  for (mU in c(6, 10, 14, 18)) {
    for (mR in c(6, 10, 14, 18)) {
      f_fix <- frame
      f_fix$take_2 <- c(mU, mR)
      res <- n_alloc_joint(
        f_fix, measures, targets = targets, psu = psu,
        certainty_cutoff = 0.85, n_psu_per_zone = 2
      )
      if (res$operational$cost < best_cost) {
        best_cost <- res$operational$cost
      }
    }
  }
  expect_equal(fit$operational$cost, best_cost, tolerance = 1e-8)

  # Round-trip prec_alloc_joint and n_alloc_joint
  pr <- prec_alloc_joint(fit)
  expect_equal(pr$detail$.achieved, fit$operational$constraints$.achieved,
               tolerance = 1e-10)
  expect_equal(pr$params$achieved$cost, fit$operational$cost, tolerance = 1e-8)

  # Budget mode with PSU register free-take
  budget <- fit$operational$cost * 1.15
  fit_b <- n_alloc_joint(
    frame, measures, targets = targets,
    objective = "y1", budget = budget, psu = psu,
    take_2_min = 6, take_2_max = 18, take_2_step = 4
  )
  expect_equal(fit_b$operational$status, "optimal")
  expect_true(fit_b$operational$all_pass)
  expect_lte(fit_b$operational$cost, budget + 1e-6)

  best_obj <- Inf
  for (mU in c(6, 10, 14, 18)) {
    for (mR in c(6, 10, 14, 18)) {
      f_fix <- frame
      f_fix$take_2 <- c(mU, mR)
      res <- tryCatch(
        n_alloc_joint(f_fix, measures, targets = targets,
                      objective = "y1", budget = budget, psu = psu),
        error = function(e) NULL
      )
      if (!is.null(res) && res$operational$objective_value < best_obj) {
        best_obj <- res$operational$objective_value
      }
    }
  }
  expect_equal(fit_b$operational$objective_value, best_obj, tolerance = 1e-10)
})

test_that("free-take 2-stage prints, summarizes, and exports design_df cleanly", {
  frame <- data.frame(
    stratum = c("A", "B"), N = c(2000, 3000), N_1 = c(50, 60),
    cost_1 = c(100, 150), cost_2 = c(10, 10)
  )
  measures <- data.frame(
    stratum = c("A", "B"), name = c("y", "y"),
    mean = c(50, 60), sd = c(15, 20), icc_1 = c(0.05, 0.08)
  )
  targets <- data.frame(name = "y", cv = 0.05)

  fit <- n_alloc_joint(frame, measures, targets = targets, take_2_min = 4, take_2_max = 25)
  expect_no_error(capture.output(print(fit)))
  expect_no_error(capture.output(summary(fit)))
  df <- design_df(fit)
  expect_s3_class(df, "svyplan_df")
  expect_equal(nrow(df$strata), 2L)
})

test_that("free-take 2-stage rejects invalid inputs clearly", {
  # 1-stage with take_2_min
  f1 <- data.frame(stratum = "A", N = 1000)
  m1 <- data.frame(stratum = "A", name = "y", mean = 50, sd = 10)
  t1 <- data.frame(name = "y", cv = 0.05)
  expect_error(
    n_alloc_joint(f1, m1, targets = t1, take_2_min = 5),
    "two-stage"
  )

  # 3-stage with NA take_2
  f3 <- data.frame(
    stratum = "A", N = 10000, N_1 = 50, N_2 = 500,
    cost_1 = 200, cost_2 = 50, cost_3 = 10,
    take_2 = NA_real_, take_3 = 5
  )
  m3 <- data.frame(
    stratum = "A", name = "y", mean = 50, sd = 10,
    icc_1 = 0.05, icc_2 = 0.1
  )
  expect_error(
    n_alloc_joint(f3, m3, targets = t1),
    "two-stage"
  )

  # Fixed take_2 with take_2_min
  f2_fix <- data.frame(
    stratum = "A", N = 2000, N_1 = 50,
    cost_1 = 100, cost_2 = 10, take_2 = 10
  )
  m2 <- data.frame(stratum = "A", name = "y", mean = 50, sd = 10, icc_1 = 0.05)
  expect_error(
    n_alloc_joint(f2_fix, m2, targets = t1, take_2_min = 5),
    "free stratum"
  )

  # Missing stage costs with free take
  f2_nocost <- data.frame(stratum = "A", N = 2000, N_1 = 50, take_2 = NA_real_)
  expect_error(
    n_alloc_joint(f2_nocost, m2, targets = t1),
    "cost_1"
  )

  # take_2_max > N / N_1 under fpc = "stage"
  f2_free <- data.frame(
    stratum = "A", N = 1000, N_1 = 50, cost_1 = 100, cost_2 = 10
  )
  expect_error(
    n_alloc_joint(f2_free, m2, targets = t1, take_2_max = 30, fpc = "stage"),
    "average PSU size"
  )

  # PSU register with PSU smaller than take_2_min
  psu <- data.frame(
    psu_id = c("P1", "P2", "P3"),
    stratum = rep("A", 3),
    N = c(5, 100, 100)
  )
  f_psu <- data.frame(stratum = "A", N = 205, cost_1 = 100, cost_2 = 10)
  expect_error(
    n_alloc_joint(f_psu, m2, targets = t1, psu = psu, take_2_min = 10),
    "take_2_min"
  )
})

## Continuous certification of the free-take route

freetake_textbook_cost <- function(frame, measures, targets, fpc, lower = 1,
                                   take_2_min = 1,
                                   take_2_max = frame$N / frame$N_1) {
  N <- frame$N
  N1 <- frame$N_1
  M <- N / N1
  rows <- measures[match(targets$name, measures$name), ]
  v <- rows$p * (1 - rows$p)
  icc <- rows$icc_1
  v_max <- (targets$cv * N * rows$p)^2
  cost <- function(m) {
    deff <- 1 + icc * (m - 1)
    a <- switch(fpc,
      unit = , none = N^2 * v * deff / m,
      stage = N^2 * v * (icc + (1 - icc) * (1 / m - 1 / M))
    )
    shift <- switch(fpc, unit = N * v * deff, none = 0, stage = N * M * v * icc)
    x <- max(a / (v_max + shift), 1, lower / m)
    if (x > min(N1, N / m) * (1 + 1e-9)) return(1e300)
    (frame$cost_1 + frame$cost_2 * m) * x
  }
  grid <- seq(take_2_min, take_2_max, length.out = 4001L)
  vals <- vapply(grid, cost, 1)
  j <- which.min(vals)
  o <- stats::optimize(
    cost, grid[c(max(1L, j - 1L), min(length(grid), j + 1L))],
    tol = 1e-12
  )
  min(vals[j], o$objective)
}

test_that("the free-take continuous optimum at a PSU-count bound is found and certified", {
  frame <- data.frame(stratum = "S1", N = 621, N_1 = 27,
                      cost_1 = 0.17, cost_2 = 2.56)
  measures <- data.frame(stratum = "S1", name = "y", p = 0.37, icc_1 = 0.31)
  targets <- data.frame(name = "y", cv = 0.01)
  fit <- n_alloc_joint(frame, measures, targets = targets)
  oracle <- freetake_textbook_cost(frame, measures, targets, "unit",
                                   lower = fit$detail$.lower)
  opt <- fit$optimization
  expect_lte(fit$params$achieved$cost, oracle * (1 + 1e-6))
  expect_equal(fit$detail$n_1, 27, tolerance = 1e-6)
  expect_equal(fit$detail$take_2, 22.89278, tolerance = 1e-5)
  expect_identical(opt$classification, "optimal")
  expect_lte(opt$relative_duality_gap, opt$kkt_tolerance)
  expect_lte(opt$dual_value, oracle * (1 + 1e-9))
  expect_lte(fit$params$achieved$cost, fit$operational$cost)
})

test_that("the free-take route never certifies a continuous design above the optimum", {
  set.seed(1)
  seen <- character(0)
  for (i in seq_len(120)) {
    fpc <- c("unit", "unit", "stage", "none")[i %% 4 + 1]
    N1 <- sample(3:30, 1)
    frame <- data.frame(
      stratum = "A", N = N1 * sample(3:30, 1), N_1 = N1,
      cost_1 = 10^runif(1, -1, 4), cost_2 = 10^runif(1, -1, 2)
    )
    measures <- data.frame(
      stratum = "A", name = paste0("y", 1:4), p = runif(4, 0.01, 0.95),
      icc_1 = 10^runif(4, -3, -0.05)
    )
    targets <- data.frame(name = measures$name, cv = 10^runif(4, -2.4, -0.1))
    info <- sprintf("%s instance %d", fpc, i)
    fit <- tryCatch(n_alloc_joint(frame, measures, targets = targets, fpc = fpc),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      expect_match(conditionMessage(fit), "infeasible", info = info)
      next
    }
    opt <- fit$optimization
    oracle <- freetake_textbook_cost(frame, measures, targets, fpc,
                                     lower = fit$detail$.lower)
    seen <- c(seen, opt$classification)
    expect_gte(fit$params$achieved$cost, oracle * (1 - 1e-6))
    expect_lte(opt$dual_value, oracle * (1 + 1e-9))
    expect_lte(fit$params$achieved$cost, fit$operational$cost * (1 + 1e-9))
    if (identical(opt$classification, "optimal")) {
      expect_lte(fit$params$achieved$cost, oracle * (1 + 1e-6))
      expect_lte(opt$relative_duality_gap, opt$kkt_tolerance)
      expect_gte(opt$relative_duality_gap, -1e-6)
    } else {
      expect_identical(opt$classification, "feasible", info = info)
      expect_gt(opt$relative_duality_gap, opt$kkt_tolerance)
    }
  }
  expect_true("optimal" %in% seen)
})

test_that("a two-stratum free-take design reaches the continuous optimum", {
  frame <- data.frame(stratum = c("A", "B"), region = c("N", "S"),
                      N = c(325, 450), N_1 = c(13, 10),
                      cost_1 = c(290, 41), cost_2 = c(16, 16))
  targets <- data.frame(name = c("a", "b", "a"),
                        domain = c(".overall", ".overall", "region"),
                        level = c(NA, NA, "N"))
  exact <- data.frame(
    stratum = c("A", "B", "A", "B"), name = c("a", "a", "b", "b"), mean = 50,
    sd = c(16.02248745, 21.14480039, 7.528977685, 22.6367637),
    icc_1 = c(0.8086754819, 0.7207273967, 0.03605451067, 0.04637340541),
    resp_rate = c(0.7764507945, 0.5294326289, 0.5780274237, 0.9691422998)
  )
  targets$cv <- c(0.05776353821, 0.02820768642, 0.06520689642)
  fit <- n_alloc_joint(frame, exact, targets = targets)
  ## 12198.1: a grid over both takes, each a convex fixed-take solve.
  expect_lte(fit$params$achieved$cost, 12198.1 * (1 + 1e-6))
  expect_lte(fit$optimization$dual_value, 12198.1 * (1 + 1e-6))
  expect_lte(fit$params$achieved$cost, fit$operational$cost)
  expect_true(fit$optimization$classification %in% c("optimal", "feasible"))

  rounded <- exact
  rounded$sd <- c(16.022487, 21.1448, 7.528978, 22.636764)
  rounded$icc_1 <- c(0.80867548, 0.7207274, 0.03605451, 0.04637341)
  rounded$resp_rate <- c(0.7764508, 0.5294326, 0.5780274, 0.9691423)
  targets$cv <- c(0.05776354, 0.02820769, 0.0652069)
  fit_r <- n_alloc_joint(frame, rounded, targets = targets)
  expect_equal(fit_r$operational$cost, 12276)
  expect_lte(fit_r$params$achieved$cost, fit_r$operational$cost)
})

test_that("a free-take design whose first dual bound is loose is still certified", {
  frame <- data.frame(stratum = c("A", "B"), N = c(425, 513), N_1 = c(17, 27),
                      cost_1 = c(44.3397741145977, 3.81739465041128),
                      cost_2 = c(10.8789360625289, 0.405358116889432))
  measures <- data.frame(
    stratum = rep(c("A", "B"), 4), name = rep(paste0("y", 1:4), each = 2),
    p = c(0.202813018038869, 0.0260903034312651, 0.480407773791812,
          0.466324549131095, 0.473794465600513, 0.0521620801975951,
          0.250823933295906, 0.944637100710534),
    icc_1 = c(0.0225403309994826, 0.301306604317985, 0.871869574703354,
              0.00108353452458534, 0.525482796871637, 0.0315158155058987,
              0.00520101613670893, 0.101420281632981)
  )
  targets <- data.frame(name = paste0("y", 1:4),
                        cv = c(0.121803356582652, 0.0519542479980103,
                               0.011264489453123, 0.0127033760692222))
  fit <- n_alloc_joint(frame, measures, targets = targets, take_2_max = c(51, 54))
  opt <- fit$optimization
  ## 5215.5653: multi-start Nelder-Mead over the takes, each a convex solve.
  expect_lte(fit$params$achieved$cost, 5215.5653 * (1 + 1e-6))
  expect_identical(opt$classification, "optimal")
  expect_lte(abs(opt$relative_duality_gap), opt$kkt_tolerance)
})

test_that("a free-take optimum on the edge of the feasible takes is reached", {
  ## Exact bits matter: the optimum sits where the PSU cap binds.
  frame <- data.frame(stratum = "A", N = 112, N_1 = 4,
                      cost_1 = 0x1.3143d22e5c3b6p-3,
                      cost_2 = 0x1.3c003c007fb9cp-3)
  measures <- data.frame(
    stratum = "A", name = paste0("y", 1:4),
    p = c(0x1.890e546deb852p-2, 0x1.e1ebf881ae147p-2, 0x1.34d4c5abf5c29p-1,
          0x1.6147afe05c28fp-1),
    icc_1 = c(0x1.b991c6d903b02p-10, 0x1.1eded88195764p-6,
              0x1.7dbc8b16a5dd3p-1, 0x1.01a5ecbd79fd8p-2)
  )
  targets <- data.frame(
    name = measures$name,
    cv = c(0x1.aab17402427efp-7, 0x1.b44fe51dcdf25p-5, 0x1.150d7d3a6064bp-8,
           0x1.0091d24752fb4p-6)
  )
  fit <- n_alloc_joint(frame, measures, targets = targets, fpc = "stage")
  oracle <- freetake_textbook_cost(frame, measures, targets, "stage",
                                   lower = fit$detail$.lower)
  expect_lte(fit$params$achieved$cost, oracle * (1 + 1e-7))
  expect_identical(fit$optimization$classification, "optimal")
  expect_true(all(fit$constraints$.ratio <= 1 + 1e-6))
})

test_that("a free-take census design is reached under the unit correction", {
  frame <- data.frame(stratum = "A", N = 120, N_1 = 6,
                      cost_1 = 0.147357733818111, cost_2 = 18.7709820518239)
  measures <- data.frame(
    stratum = "A", name = paste0("y", 1:4),
    p = c(0.0160983291547745, 0.731206055628136, 0.86502676170785,
          0.717094809561968),
    icc_1 = c(0.0227719753833856, 0.0427764141909372, 0.00697165759685054,
              0.0218147980390309)
  )
  targets <- data.frame(name = measures$name,
                        cv = c(0.00457457541857384, 0.0235637937935708,
                               0.0268362286729755, 0.708292137835873))
  fit <- n_alloc_joint(frame, measures, targets = targets, take_2_max = 20)
  oracle <- freetake_textbook_cost(frame, measures, targets, "unit",
                                   lower = fit$detail$.lower, take_2_max = 20)
  expect_lte(fit$params$achieved$cost, oracle * (1 + 1e-7))
  expect_identical(fit$optimization$classification, "optimal")
  expect_true(all(fit$constraints$.ratio <= 1 + 1e-6))
})

test_that("a free-take budget is refused only below the dual bound", {
  frame <- data.frame(stratum = c("A", "B"), N = c(532, 30), N_1 = c(28, 5),
                      cost_1 = c(7.07063587717158, 0.578738214079003),
                      cost_2 = c(22.5576906511801, 16.9300337906728))
  measures <- data.frame(stratum = c("A", "B"), name = "y1",
                         p = c(0.432623439654708, 0.424562202217057),
                         icc_1 = c(0.0477240467650706, 0.793064531126775))
  targets <- data.frame(name = "y1", cv = 0.00733288432015429)
  take_2_max <- c(84, 10)
  base <- n_alloc_joint(frame, measures, targets = targets, take_2_max = take_2_max)
  opt <- base$optimization
  expect_identical(opt$classification, "feasible")
  low <- opt$dual_value
  high <- base$params$achieved$cost
  expect_lt(low, high)
  expect_output(print(base), "continuous: feasible, within 0.582% of the dual bound")
  expect_output(print(base), "continuous design: n = ")
  expect_output(print(summary(base)), "Continuous design, feasible, within")

  run <- function(budget) {
    n_alloc_joint(frame, measures, targets = targets, objective = "y1",
                  budget = budget, take_2_max = take_2_max)
  }
  expect_error(run(low * 0.999), "costs at least")
  expect_error(run(0.5 * (low + high)), class = "svyplan_plan_unknown")
  fit <- run(high * 1.002)
  expect_identical(fit$optimization$classification, "feasible")
  expect_lte(fit$params$achieved$cost, high * 1.002 * (1 + 1e-6))
  expect_lte(fit$operational$cost, high * 1.002 * (1 + 1e-6))
  expect_true(fit$operational$all_pass)

  objective_only <- function(budget) {
    n_alloc_joint(frame, measures, objective = "y1", budget = budget,
                  take_2_max = take_2_max)$optimization
  }
  expect_identical(objective_only(8000)$classification, "optimal")
  expect_identical(objective_only(high * 1.002)$classification, "feasible")
})

test_that("a looser target never worsens the free-take budget objective", {
  ## Invariant oracle: a slack target leaves the optimum unchanged. Under
  ## "unit" and "stage" a census gives a zero objective floor, and the
  ## first root-search trial sits just above it.
  frame <- data.frame(stratum = c("A", "B"), N = c(18, 24), N_1 = c(6, 6),
                      cost_1 = c(12, 8), cost_2 = c(2, 3))
  measures <- data.frame(
    stratum = rep(c("A", "B"), 2), name = rep(c("cov", "inc"), each = 2),
    p = c(0.35, 0.45, NA, NA), mean = c(NA, NA, 100, 120),
    sd = c(NA, NA, 25, 30), icc_1 = c(0.34, 0.08, 0.04, 0.42)
  )
  objective <- data.frame(name = c("cov", "inc"), priority = c(1, 28))
  run <- function(cv, fpc) {
    n_alloc_joint(frame, measures, targets = data.frame(name = "cov", cv = cv),
                  objective = objective, budget = 190, fpc = fpc,
                  take_2_min = 1, take_2_max = c(3, 4))
  }
  caps <- list(
    unit = c(0.2, 0.22, 0.23, 0.235, 0.23501, 0.25, 0.27, 0.29, 1),
    stage = c(0.2, 0.205, 0.22, 0.25, 1)
  )
  for (fpc in names(caps)) {
    floor_fit <- run(1, fpc)
    expect_identical(floor_fit$optimization$classification, "optimal")
    for (cv in caps[[fpc]]) {
      fit <- run(cv, fpc)
      label <- sprintf("%s cv = %s", fpc, cv)
      expect_lt(fit$constraints$.cv, cv, label = label)
      expect_equal(fit$objective_value, floor_fit$objective_value,
                   tolerance = 1e-8, label = label)
      expect_equal(fit$params$achieved$cost, 190, tolerance = 1e-6,
                   label = label)
      expect_identical(fit$optimization$classification, "optimal",
                       label = label)
    }
  }
})

test_that("the stratum minimum of the free-take dual is global in every FPC mode", {
  stratum_min <- svyplan:::.freetake_stratum_min
  set.seed(7)
  for (fpc in c("none", "stage", "unit")) {
    for (i in seq_len(50)) {
      N1 <- sample(2:60, 1)
      N <- round(N1 * exp(runif(1, log(3), log(80))))
      M <- N / N1
      v <- exp(runif(1, log(0.05), log(500)))
      icc <- runif(1, 0, 0.9)
      r2 <- runif(1, 0.3, 1)
      a <- switch(fpc, stage = N^2 * v * (icc - (1 - icc) / M), N^2 * v * icc)
      b <- N^2 * v * (1 - icc) / r2
      d <- if (fpc == "unit") N * v * icc * r2 else 0
      lam <- exp(runif(1, log(1e-8), log(1e2))) / (N^2 * v)
      m_lo <- if (runif(1) < 0.5) 1 else runif(1, 1, 5)
      m_hi <- if (runif(1) < 0.3) m_lo else min(M, m_lo + exp(runif(1, 0, log(80))))
      lo <- if (runif(1) < 0.5) 1 else runif(1, 1, min(N, N1 * m_hi))
      c1 <- exp(runif(1, log(0.1), log(500)))
      c2 <- exp(runif(1, log(0.1), log(50)))
      s <- stratum_min(c1, c2, lam * a, lam * b, lam * d, N1, N, lo, m_lo, m_hi)
      best <- Inf
      for (x in seq(1, N1, length.out = 200)) {
        y_lo <- max(lo, m_lo * x)
        y_hi <- min(N, m_hi * x)
        if (y_lo > y_hi) next
        y <- seq(y_lo, y_hi, length.out = 200)
        best <- min(best, c1 * x + c2 * y + lam * (a / x + b / y - d * y / x))
      }
      scale <- c1 * s$x + c2 * s$y + lam * (abs(a) / s$x + b / s$y + d * s$y / s$x)
      expect_lte(s$val, best + 1e-9 * scale, label = sprintf("%s %d", fpc, i))
    }
  }
})

## Numbered take bounds

bound_fixture <- function() {
  list(
    frame = data.frame(stratum = c("U", "R"), N = c(1600, 2000),
                       N_1 = c(40, 50), cost_1 = c(110, 190),
                       cost_2 = c(14, 10)),
    measures = data.frame(stratum = rep(c("U", "R"), 2),
                          name = rep(c("p1", "y1"), each = 2),
                          p = c(0.35, 0.45, NA, NA),
                          mean = c(NA, NA, 100, 120),
                          sd = c(NA, NA, 25, 30),
                          icc_1 = c(0.07, 0.04, 0.03, 0.08)),
    targets = data.frame(name = c("p1", "y1"), cv = c(0.08, 0.04))
  )
}

register_bound_fixture <- function() {
  frame <- data.frame(stratum = c("A", "B"), N = c(2000, 1600),
                      cost_1 = 200, cost_2 = 5)
  list(
    frame = frame,
    psu = data.frame(stratum = rep(frame$stratum, each = 8),
                     N = c(600, rep(200, 7), 450, 220, 165, rep(153, 5))),
    measures = data.frame(stratum = frame$stratum, name = "y",
                          p = c(0.5, 0.4), icc_1 = 0.05),
    targets = data.frame(name = "y", cv = 0.1)
  )
}

field_cols <- c("n_1", "take_2", "take_2_int", "n_int")

test_that("take bounds give one fit as arguments or frame columns, and arguments win", {
  fx <- bound_fixture()
  fit <- function(frame, ...) {
    n_alloc_joint(frame, fx$measures, targets = fx$targets, ...)
  }
  open <- fit(fx$frame)
  by_arg <- fit(fx$frame, take_2_min = 6, take_2_max = 18, take_2_step = 4)
  by_col <- fit(transform(fx$frame, take_2_min = 6, take_2_max = 18,
                          take_2_step = 4))
  both <- fit(transform(fx$frame, take_2_min = 1, take_2_max = 30,
                        take_2_step = 1),
              take_2_min = 6, take_2_max = 18, take_2_step = 4)

  expect_false(all(open$detail$take_2_int %in% seq(6, 18, by = 4)))
  expect_true(all(by_arg$detail$take_2_int %in% seq(6, 18, by = 4)))
  expect_equal(by_col$detail[field_cols], by_arg$detail[field_cols])
  expect_equal(both$detail[field_cols], by_arg$detail[field_cols])
  expect_identical(by_arg$params$take_2_step, 4)
})

test_that("a non-default take step survives size to precision to size", {
  fx <- bound_fixture()
  fit <- n_alloc_joint(fx$frame, fx$measures, targets = fx$targets,
                       take_2_min = 6, take_2_max = 18, take_2_step = 4)
  again <- n_alloc_joint(prec_alloc_joint(fit))
  expect_identical(again$params$take_2_step, 4)
  expect_equal(again$detail[field_cols], fit$detail[field_cols])

  rx <- register_bound_fixture()
  reg <- n_alloc_joint(rx$frame, rx$measures, targets = rx$targets,
                       psu = rx$psu, take_2_min = 2, take_2_max = 17,
                       take_2_step = 6)
  expect_true(all(reg$detail$take_2_int %in% c(2, 8, 14)))
  reg_again <- n_alloc_joint(prec_alloc_joint(reg))
  expect_identical(reg_again$params$take_2_step, 6)
  expect_equal(reg_again$detail[c("take_2_int", "n_int")],
               reg$detail[c("take_2_int", "n_int")])
})

test_that("a bound sweep row equals a direct refit at the fitted step", {
  fx <- bound_fixture()
  fit <- function(...) {
    n_alloc_joint(fx$frame, fx$measures, targets = fx$targets,
                  take_2_step = 4, ...)
  }
  base <- fit(take_2_min = 6, take_2_max = 18)
  rows <- data.frame(take_2_min = c(6, 10), take_2_max = c(18, 15))
  sweep <- predict(base, rows)
  rx <- register_bound_fixture()
  reg_fit <- function(...) {
    n_alloc_joint(rx$frame, rx$measures, targets = rx$targets,
                  psu = rx$psu, take_2_step = 6, ...)
  }
  reg_base <- reg_fit(take_2_min = 2, take_2_max = 20)
  reg_rows <- data.frame(take_2_min = c(2, 8), take_2_max = c(17, 20))
  reg_sweep <- predict(reg_base, reg_rows)

  checked <- character()
  for (case in list(list("plain", fit, rows, sweep),
                    list("register", reg_fit, reg_rows, reg_sweep))) {
    for (i in seq_len(nrow(case[[3]]))) {
      direct <- case[[2]](take_2_min = case[[3]]$take_2_min[i],
                          take_2_max = case[[3]]$take_2_max[i])
      row <- case[[4]][i, ]
      label <- sprintf("%s row %d", case[[1]], i)
      expect_equal(row$n, direct$n, tolerance = 1e-8, label = label)
      expect_equal(row$n_int, sum(direct$detail$n_int), label = label)
      expect_equal(row$cost_int, direct$operational$cost, label = label)
      checked <- c(checked, label)
    }
  }
  expect_setequal(checked, c("plain row 1", "plain row 2",
                             "register row 1", "register row 2"))
})

test_that("free-take cost slopes need a certified fit with verified multipliers", {
  frame <- data.frame(stratum = "A", N = 1000, N_1 = 100, cost_1 = 10,
                      cost_2 = 1)
  measures <- data.frame(stratum = "A", name = c("early", "late"), mean = 1,
                         sd = 1, icc_1 = c(0.8, 0.2))
  cost_fit <- function(cv) {
    n_alloc_joint(frame, measures,
                  targets = data.frame(name = "early", cv = cv), fpc = "none")
  }
  ## Oracle: the minimum cost refitted on both sides of the target.
  fit <- cost_fit(0.28)
  expect_identical(fit$optimization$classification, "optimal")
  expect_identical(fit$optimization$slopes, "available")
  expect_equal(fit$constraints$.sensitivity,
               (cost_fit(0.28001)$params$achieved$cost -
                  cost_fit(0.27999)$params$achieved$cost) / 2e-5,
               tolerance = 1e-5)

  ## A design certified only as feasible reports no slope.
  frame <- data.frame(stratum = c("A", "B"), N = c(532, 30), N_1 = c(28, 5),
                      cost_1 = c(7.07063587717158, 0.578738214079003),
                      cost_2 = c(22.5576906511801, 16.9300337906728))
  measures <- data.frame(stratum = c("A", "B"), name = "y1",
                         p = c(0.432623439654708, 0.424562202217057),
                         icc_1 = c(0.0477240467650706, 0.793064531126775))
  feasible <- n_alloc_joint(frame, measures,
                            targets = data.frame(name = "y1",
                                                 cv = 0.00733288432015429),
                            take_2_max = c(84, 10))
  expect_identical(feasible$optimization$classification, "feasible")
  expect_identical(feasible$optimization$slopes, "not_optimal")
  expect_true(all(is.na(feasible$constraints$.sensitivity)))
  expect_output(print(summary(feasible)), "Sens. is not reported")
})
