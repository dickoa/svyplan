# Joint allocation against published and independently solved results
## Bethel (1989) educational survey, as set up by Willems and Muennich (2026)

test_that("the Bethel minimum cost at a CV of 0.06 matches an independent solve", {
  ## Reference values from an independent conic solve of the same problem.
  ## The CVs and the multiplier shares agree with the paper's f(n*) and
  ## lambda* to the digits it prints (Section 6).
  ref <- bethel_reference()
  fit <- n_alloc_joint(ref$frame, ref$measures,
                       targets = data.frame(name = rownames(ref$d_qh), cv = 0.06),
                       fpc = "none")
  expect_identical(fit$optimization$classification, "optimal")
  expect_equal(fit$params$achieved$cost, 241.13987, tolerance = 1e-6)
  expect_equal(fit$constraints$.cv, c(0.06, 0.06, 0.04820675, 0.05025695),
               tolerance = 1e-6)
  ## Every target shares one total and one CV, so the multipliers normalized
  ## to sum 1 are the paper's lambda*.
  m <- fit$constraints$.multiplier
  expect_equal(m / sum(m), c(0.6669, 0.3331, 0, 0), tolerance = 1e-3)
  expect_identical(m[3:4], c(0, 0))
  ## Unit costs make every whole cost an integer, so 242 is the least whole
  ## cost above 241.14.
  expect_identical(fit$operational$cost, 242)
  expect_true(fit$operational$all_pass)
})

test_that("the equal-priority sum within a budget matches its closed form", {
  ## Closed form: with unit costs and no binding bound or target, the sum of
  ## relative variances sum_h D_h / n_h is least at n_h proportional to
  ## sqrt(D_h), where D_h = sum_q d_qh.
  ref <- bethel_reference()
  fit <- n_alloc_joint(ref$frame, ref$measures,
                       objective = rownames(ref$d_qh), budget = 241,
                       fpc = "none")
  D <- colSums(ref$d_qh)
  n <- 241 * sqrt(D) / sum(sqrt(D))
  relvar <- drop(ref$d_qh %*% (1 / n))
  expect_identical(fit$optimization$classification, "optimal")
  expect_equal(fit$detail$n, unname(n), tolerance = 1e-8)
  expect_equal(fit$objective$.relvar, unname(relvar), tolerance = 1e-8)
  expect_equal(fit$objective_value, sum(relvar), tolerance = 1e-8)
  expect_equal(max(sqrt(fit$objective$.relvar)), 0.0645398, tolerance = 1e-6)
})

test_that("regional CV caps within a budget match an independent solve", {
  ## The paper's regional Floorspace caps (Section 6), without its regional
  ## minimum sample sizes. Reference values from an independent conic solve
  ## with stratum bounds 1 and N_h. The paper bounds each region's
  ## sum_h d_1h / n_h by 0.04^2, a relative variance on the national total,
  ## which is a regional CV of 0.04 * sum(N) / N_r.
  ref <- bethel_reference()
  N_r <- tapply(ref$N, ref$frame$region, sum)
  targets <- data.frame(name = "Floorspace", domain = "region",
                        level = names(N_r),
                        cv = 0.04 * sum(ref$N) / as.numeric(N_r))
  fit <- n_alloc_joint(ref$frame, ref$measures, targets = targets,
                       objective = rownames(ref$d_qh), budget = 241,
                       fpc = "none")
  opt <- fit$optimization
  expect_identical(opt$classification, "optimal")
  expect_equal(fit$params$achieved$cost, 241, tolerance = 1e-8)
  expect_equal(fit$objective_value, 0.0115320035, tolerance = 1e-6)
  ## The tolerance is relative only for values above it, hence the scale.
  expect_equal(-1e5 * opt$budget_sensitivity, 5.240969, tolerance = 1e-3)
  ## The multiplier prices a cap on the regional variance in cost units.
  ## Times sum(N)^2 and over the objective multiplier, it is the objective's
  ## price per unit of the paper's normalized bound.
  price <- fit$constraints$.multiplier * sum(ref$N)^2 /
    opt$objective_multiplier
  expect_identical(price[1], 0)
  expect_equal(price[2:3], c(0.59553, 0.091177), tolerance = 1e-3)
})

## Objective change per unit of target at a held budget

budget_slope_fit <- function(targets, budget = 241, objective = "Floorspace",
                             ...) {
  ref <- bethel_reference()
  n_alloc_joint(ref$frame, ref$measures, targets = targets,
                objective = objective, budget = budget, fpc = "none", ...)
}

summary_text <- function(fit) {
  gsub("\\s+", " ", paste(capture.output(print(summary(fit))), collapse = " "))
}

centered_slope <- function(make, value, h, ...) {
  up <- budget_slope_fit(make(value + h), ...)$objective_value
  down <- budget_slope_fit(make(value - h), ...)$objective_value
  (up - down) / (2 * h)
}

test_that("the objective slope at a held budget matches centered refits", {
  ## Oracle: the optimal objective refitted on both sides of the target,
  ## which does not read the multipliers.
  age <- function(cv) data.frame(name = "Age", cv = cv)
  fit <- budget_slope_fit(age(0.061))
  slope <- fit$constraints$.objective_sensitivity
  expect_identical(fit$optimization$slopes, "available")
  expect_equal(slope, -0.05035477, tolerance = 1e-6)
  expect_equal(slope, fit$constraints$.sensitivity /
                 fit$optimization$objective_multiplier, tolerance = 1e-12)
  for (h in c(1e-4, 1e-5, 1e-6)) {
    expect_equal(centered_slope(age, 0.061, h), slope, tolerance = 1e-4,
                 label = sprintf("h = %g", h))
  }
  for (metric in c("moe", "rmoe")) {
    target <- function(value) {
      out <- data.frame(name = "Age")
      out[[metric]] <- value
      out
    }
    value <- if (metric == "moe") 0.12 else 0.13
    fit <- budget_slope_fit(target(value), df = 12)
    expect_true(fit$constraints$.binding, label = metric)
    expect_equal(centered_slope(target, value, 1e-6, df = 12),
                 fit$constraints$.objective_sensitivity, tolerance = 1e-5,
                 label = metric)
  }
})

test_that("a slack target has a zero slope and an unidentified one is missing", {
  fit <- budget_slope_fit(data.frame(name = c("Age", "Oil"),
                                     cv = c(0.061, 0.2)))
  expect_identical(fit$constraints$.binding, c(TRUE, FALSE))
  expect_identical(fit$constraints$.objective_sensitivity[2], 0)
  expect_false(is.na(fit$constraints$.objective_sensitivity[1]))

  ## A CV target and the RMOE target it implies are the same constraint.
  twin <- data.frame(name = c("Age", "Age"), cv = c(0.061, NA),
                     rmoe = c(NA, 0.061 * qnorm(0.975)))
  ## The two rows share one price, so neither has a slope of its own, while
  ## the objective's price stays its own.
  fit <- budget_slope_fit(twin)
  expect_identical(fit$optimization$slopes, "available")
  expect_identical(fit$constraints$.objective_sensitivity,
                   c(NA_real_, NA_real_))
  expect_identical(fit$constraints$.sensitivity, c(NA_real_, NA_real_))

  loose <- budget_slope_fit(data.frame(name = "Age", cv = 0.061),
                            budget = 1e6)
  expect_identical(loose$optimization$slopes,
                   "budget_not_binding")
  expect_identical(loose$constraints$.objective_sensitivity, 0)
  expect_match(summary_text(loose), "The budget does not bind", fixed = TRUE)

  ## The objective covers region r1 only, so it reaches its floor below the
  ## budget while the r3 target still binds at its cheapest design.
  ref <- bethel_reference()
  r3 <- function(cv) {
    data.frame(name = "Floorspace", domain = "region", level = "r3", cv = cv)
  }
  floor_fit <- function(cv) {
    n_alloc_joint(ref$frame, ref$measures, targets = r3(cv),
                  objective = data.frame(name = "Floorspace",
                                         domain = "region", level = "r1"),
                  budget = 7990, fpc = "none")
  }
  fit <- floor_fit(0.8)
  expect_identical(fit$optimization$slopes,
                   "budget_not_binding")
  expect_true(fit$constraints$.binding)
  expect_identical(fit$constraints$.objective_sensitivity, 0)
  moved <- floor_fit(0.8 + 1e-4)$objective_value -
    floor_fit(0.8 - 1e-4)$objective_value
  expect_lt(abs(moved), 1e-8 * fit$objective_value)
})

test_that("the slope is missing at the cheapest budget and exact above it", {
  ## The whole search cannot fund a continuous minimum cost, so the
  ## continuous engine is called as the register route calls it.
  ref <- bethel_reference()
  targets <- data.frame(name = rownames(ref$d_qh), cv = 0.06)
  continuous <- function(targets, budget) {
    svyplan:::.n_alloc_bethel(
      frame = ref$frame, measures = ref$measures, targets = targets,
      unit_cost = NULL, alpha = 0.05, deff = 1, resp_rate = 1,
      eligible_per_unit = 1, eligibility_rate = 1, df = NULL,
      min_n_stratum = NULL, objective = rownames(ref$d_qh),
      budget = budget, fpc = "none", .finite_prop_var = TRUE,
      .whole = FALSE
    )
  }
  cheapest <- n_alloc_joint(ref$frame, ref$measures, targets = targets,
                            fpc = "none")$params$achieved$cost
  at_min <- continuous(targets, cheapest)
  expect_identical(at_min$optimization$slopes,
                   "budget_at_minimum")
  expect_identical(at_min$constraints$.objective_sensitivity,
                   c(NA, NA, 0, 0))
  expect_true(is.na(at_min$optimization$budget_sensitivity))
  expect_true(all(is.na(at_min$constraints$.sensitivity[1:2])))
  budget <- 1.01 * cheapest
  above <- continuous(targets, budget)
  expect_identical(above$optimization$slopes, "available")
  h <- 1e-6
  relaxed <- function(cv) {
    out <- targets
    out$cv[1] <- cv
    continuous(out, budget)$objective_value
  }
  expect_equal((relaxed(0.06 + h) - relaxed(0.06 - h)) / (2 * h),
               above$constraints$.objective_sensitivity[1],
               tolerance = 1e-5)
})

test_that("only the rows in a dependency lose their own prices", {
  ## Oracle: a row's price is its own exactly when every null vector of the
  ## gradients vanishes on it, which these small matrices show by hand.
  shared <- svyplan:::.ambiguous_columns
  A <- cbind(c(1, 2, 3, 0), c(3, 1, 1, 0), c(0, 0, 1, 5))
  expect_identical(shared(A), c(FALSE, FALSE, FALSE))
  ## A doubled column trades with its original, and the third keeps its own.
  expect_identical(shared(cbind(A[, 1], A[, 3], 2 * A[, 1])),
                   c(TRUE, FALSE, TRUE))
  ## A sum of two columns trades with both of them.
  expect_identical(shared(cbind(A, A[, 1] + A[, 2])),
                   c(TRUE, TRUE, FALSE, TRUE))
  ## A zero column has no price, and more columns than directions share.
  expect_identical(shared(cbind(A[, 1], 0)), c(FALSE, TRUE))
  expect_true(all(shared(cbind(A[1:2, 1:2], c(1, 1)))))
})

test_that("a register budget fit leaves the objective slope missing", {
  frame <- data.frame(stratum = c("A", "B"), N = c(2000, 1600),
                      cost_1 = 200, cost_2 = 5, take_2 = c(10, 8))
  psu <- data.frame(stratum = rep(frame$stratum, each = 8),
                    N = c(600, rep(200, 7), 450, 220, 165, rep(153, 5)))
  measures <- data.frame(stratum = frame$stratum, name = "y",
                         p = c(0.5, 0.4), icc_1 = 0.05)
  fit <- n_alloc_joint(frame, measures,
                       targets = data.frame(name = "y", cv = 0.2),
                       objective = "y", budget = 3000, psu = psu)
  expect_identical(fit$optimization$slopes, "register")
  expect_true(all(is.na(fit$constraints$.objective_sensitivity)))
  expect_true(all(is.na(fit$constraints$.sensitivity)))
  expect_true(is.na(fit$optimization$budget_sensitivity))
  cheapest <- n_alloc_joint(frame, measures,
                            targets = data.frame(name = "y", cv = 0.1),
                            psu = psu)
  expect_identical(cheapest$optimization$slopes, "register")
  expect_true(all(is.na(cheapest$constraints$.sensitivity)))
  expect_match(summary_text(fit), "not reported for a PSU register",
               fixed = TRUE)
})

test_that("a budget fit saved before the objective slope still prints", {
  fit <- budget_slope_fit(data.frame(name = "Age", cv = 0.061))
  expect_match(summary_text(fit), "the budget held at 241.", fixed = TRUE)
  fit$constraints$.objective_sensitivity <- NULL
  fit$optimization$slopes <- NULL
  expect_output(print(fit), "question")
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("Obj. sens.", out, fixed = TRUE)))
})

## Slope status across modes and routes

test_that("a shared price leaves the other targets their own slopes", {
  ## The two Age rows are one constraint, and Oil binds on its own. Oracle
  ## for Oil's slope: centered refits.
  ref <- bethel_reference()
  run <- function(oil) {
    n_alloc_joint(ref$frame, ref$measures,
                  targets = data.frame(name = c("Age", "Age", "Oil"),
                                       cv = c(0.061, NA, oil),
                                       rmoe = c(NA, 0.061 * qnorm(0.975), NA)),
                  objective = "Floorspace", budget = 241, fpc = "none")
  }
  fit <- run(0.05)
  expect_identical(fit$constraints$.binding, c(TRUE, TRUE, TRUE))
  expect_identical(fit$optimization$slopes, "available")
  expect_true(all(is.na(fit$constraints$.objective_sensitivity[1:2])))
  expect_true(all(is.na(fit$constraints$.sensitivity[1:2])))
  h <- 1e-6
  expect_equal(fit$constraints$.objective_sensitivity[3],
               (run(0.05 + h)$objective_value -
                  run(0.05 - h)$objective_value) / (2 * h),
               tolerance = 1e-5)
  expect_true(is.finite(fit$optimization$budget_sensitivity))
})

test_that("coincident targets in a minimum-cost fit have no slopes", {
  z <- .bethel_fixture()
  copy <- transform(z$measures[z$measures$name == "vaccination", ],
                    name = "other")
  fit <- n_alloc_joint(z$frame, rbind(z$measures, copy),
                       targets = data.frame(name = c("vaccination", "other"),
                                            cv = 0.05))
  expect_equal(fit$constraints$.achieved, c(0.05, 0.05), tolerance = 1e-6)
  expect_true(all(is.na(fit$constraints$.sensitivity)))
  expect_identical(fit$optimization$slopes, "available")
})

test_that("every slope is withheld at a whole budget floor", {
  ## One stratum whose cheapest design meeting the target costs exactly 100.
  fit <- n_alloc_joint(data.frame(stratum = "A", N = 1000),
                       data.frame(stratum = "A", name = "y", mean = 1, sd = 1),
                       targets = data.frame(name = "y", cv = 0.1),
                       objective = "y", budget = 100, fpc = "none")
  expect_identical(fit$operational$cost, 100)
  expect_identical(fit$optimization$slopes, "budget_at_minimum")
  expect_true(is.na(fit$optimization$budget_sensitivity))
  expect_true(is.na(fit$constraints$.objective_sensitivity))
  expect_true(is.na(fit$constraints$.sensitivity))
  expect_match(summary_text(fit), "budget sensitivity: not reported",
               fixed = TRUE)
})

test_that("the worst rule ties beyond an expensive stratum and reaches a census", {
  ## Stratum c lies outside the objective and costs 100 per unit, so the
  ## budget of 104 buys both domains whole and one unit of c.
  frame <- data.frame(stratum = letters[1:3], district = letters[1:3],
                      N = 2, unit_cost = c(1, 1, 100))
  measures <- data.frame(stratum = letters[1:2], name = "y", mean = 1,
                         sd = c(sqrt(2), 1))
  objective <- data.frame(name = "y", domain = "district", level = letters[1:2])
  for (fpc in c("none", "unit")) {
    fit <- n_alloc_joint(frame, measures, objective = objective, budget = 104,
                         objective_rule = "worst", fpc = fpc)
    expect_equal(fit$detail$n_int, c(2, 2, 1), label = fpc)
    expect_identical(fit$operational$status, "optimal", label = fpc)
    expect_equal(fit$objective_value, if (fpc == "unit") 0 else 1,
                 tolerance = 1e-12, label = fpc)
  }
})
