## tradeoff_alloc()

test_that("a trade-off on the education survey prices the local rate", {
  ## Oracle: the first-order rate of the frontier at the continuous worst-
  ## rule optimum, checked against the refit at a multiplier near 1.
  z <- bethel_reference()
  fit <- n_alloc_joint(z$frame, z$measures, objective = rownames(z$d_qh),
                       budget = 241, objective_rule = "worst", fpc = "none")
  res <- tradeoff_alloc(fit, "Floorspace@.overall", "Age@.overall",
                        c(0.999, 0.95), baseline = "continuous")
  expect_identical(res$continuous_status, c("optimal", "infeasible"))
  expect_identical(res$field_status[2], "infeasible")
  expect_equal(res$estimated_allow_cv[1], res$continuous_allow_cv[1],
               tolerance = 2e-4)
  comp <- res$components[[1]]
  held <- !is.na(comp$cap_cv)
  expect_true(all(comp$continuous_cv[held] <= comp$cap_cv[held] * (1 + 1e-6)))
  expect_true(all(comp$field_cv[held] <= comp$cap_cv[held] * (1 + 1e-6)))
  expect_lte(res$field_cost[1], 241)
  expect_lte(res$field_allow_lower_bound[1],
             res$field_allow_cv[1] * (1 + 1e-9))
  expect_null(res$n_field[[2]])
})

test_that("a trade-off's field design is the best whole design under the caps", {
  ## Oracle: every whole design of three strata within the budget, with the
  ## other components at their baseline field CVs and the improved one at
  ## its requested CV.
  frame <- data.frame(stratum = letters[1:3], N = c(8, 7, 9),
                      unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 3),
                         name = rep(c("x", "y", "z"), each = 3), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1, 2, 1, 3))
  fit <- n_alloc_joint(frame, measures, objective = c("x", "y", "z"),
                       budget = 24, fpc = "none")
  res <- tradeoff_alloc(fit, "x@.overall", "y@.overall", c(1, 0.97, 0.9))
  states <- as.matrix(expand.grid(lapply(frame$N, seq_len)))
  total <- sum(frame$N)
  cv <- vapply(c("x", "y", "z"), function(nm) {
    sd <- measures$sd[measures$name == nm]
    sqrt(drop((1 / states) %*% (frame$N * sd)^2)) / total
  }, numeric(nrow(states)))
  base <- fit$operational$objective$.cv
  within <- drop(states %*% frame$unit_cost) <= 24
  for (r in seq_len(nrow(res))) {
    ok <- within & cv[, "x"] <= res$requested_cv[r] * (1 + 1e-9) &
      cv[, "z"] <= base[3] * (1 + 1e-9)
    label <- sprintf("cv_mult %g", res$cv_mult[r])
    if (!any(ok)) {
      expect_identical(res$field_status[r], "infeasible", label = label)
      next
    }
    expect_identical(res$field_status[r], "optimal", label = label)
    expect_equal(res$field_allow_cv[r], min(cv[ok, "y"]), tolerance = 1e-10,
                 label = label)
  }
  ## The original design meets every cap at cv_mult = 1.
  expect_lte(res$field_allow_cv[1], base[2] * (1 + 1e-12))
})

test_that("a trade-off keeps a continuous result and protects a census", {
  frame <- data.frame(stratum = LETTERS[1:3], region = LETTERS[1:3], N = 30)
  measures <- data.frame(stratum = frame$stratum, name = "y", mean = 1,
                         sd = 1)
  objective <- data.frame(name = "y", domain = "region",
                          level = frame$region,
                          component = c("Improve", "Allow", "Protect"))
  fit <- n_alloc_joint(frame, measures, objective = objective, budget = 10,
                       min_n_stratum = c(1, 2, 1), fpc = "none")
  res <- tradeoff_alloc(fit, "Improve", "Allow", 0.88,
                        baseline = "continuous")
  expect_identical(res$continuous_status, "optimal")
  expect_identical(res$field_status, "infeasible")
  expect_lte(res$continuous_improve_cv, res$requested_cv * (1 + 1e-6))
  expect_null(res$n_field[[1]])

  frame$N[3] <- 3
  frame$take_all <- c(FALSE, FALSE, TRUE)
  fit <- n_alloc_joint(frame, measures, objective = objective, budget = 13)
  res <- tradeoff_alloc(fit, "Improve", "Allow", 0.99)
  expect_identical(res$field_status, "optimal")
  expect_identical(res$components[[1]]$field_cv[3], 0)
  expect_equal(unname(res$n_field[[1]][3]), 3)
})

test_that("a trade-off keeps the fit's floors and targets", {
  frame <- data.frame(stratum = letters[1:4], N = c(400, 300, 500, 200),
                      region = c("N", "N", "S", "S"),
                      unit_cost = c(2, 1, 3, 1))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("y", "z"), each = 4), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1, 2, 1))
  floor <- data.frame(domain = "region", level = "S", min_n = 150)
  targets <- data.frame(name = "y", domain = "region", level = "N",
                        cv = 0.12)
  fit <- n_alloc_joint(frame, measures, targets = targets,
                       objective = c("y", "z"), budget = 900, fpc = "none",
                       min_n_domain = floor)
  res <- tradeoff_alloc(fit, "z@.overall", "y@.overall", 0.98)
  expect_identical(res$field_status, "optimal")
  n <- res$n_field[[1]]
  expect_gte(sum(n[3:4]), 150)
  assessed <- prec_alloc_joint(fit, n = unname(n))
  expect_true(all(assessed$detail$.pass))
  expect_true(all(assessed$domain_floors$.pass))
  expect_lte(sum(n * frame$unit_cost), 900)
})

test_that("tradeoff_alloc() refuses what it cannot compare", {
  z <- bethel_reference()
  fit <- n_alloc_joint(z$frame, z$measures, objective = rownames(z$d_qh),
                       budget = 241, fpc = "none")
  cost <- n_alloc_joint(z$frame, z$measures,
                        targets = data.frame(name = "Age", cv = 0.08),
                        fpc = "none")
  expect_error(tradeoff_alloc(cost, "Age@.overall", "Oil@.overall"),
               "budget fit")
  expect_error(tradeoff_alloc(fit, "Age@.overall", "Age@.overall"),
               "two different components")
  expect_error(tradeoff_alloc(fit, "Age@.overall", "Missing"),
               "two different components")
  expect_error(tradeoff_alloc(fit, "Age@.overall", "Oil@.overall", 1.1),
               "at most 1")
})

test_that("a trade-off holds a census the budget bought", {
  ## C is not take-all, but the fit takes all 3 of its units, so its CV is
  ## 0 and must stay 0.
  frame <- data.frame(stratum = LETTERS[1:3], region = LETTERS[1:3],
                      N = c(30, 30, 3))
  measures <- data.frame(stratum = frame$stratum, name = "y", mean = 1,
                         sd = c(1, 1, 5))
  objective <- data.frame(name = "y", domain = "region",
                          level = frame$region,
                          component = c("Improve", "Allow", "Protect"))
  fit <- n_alloc_joint(frame, measures, objective = objective, budget = 15)
  expect_identical(as.numeric(fit$detail$n_int[3]), 3)
  res <- tradeoff_alloc(fit, "Improve", "Allow", 0.95)
  expect_identical(res$field_status, "optimal")
  expect_identical(res$components[[1]]$field_cv[3], 0)
})

test_that("a trade-off tightens a target the fit holds on the same estimate", {
  frame <- data.frame(stratum = letters[1:3], N = c(8, 7, 9),
                      unit_cost = c(2, 1, 3))
  measures <- data.frame(stratum = rep(frame$stratum, 2),
                         name = rep(c("x", "y"), each = 3), mean = 1,
                         sd = c(1, 2, 2, 3, 1, 1))
  fit <- n_alloc_joint(frame, measures,
                       targets = data.frame(name = "x", cv = 0.5),
                       objective = c("x", "y"), budget = 40, fpc = "none")
  res <- tradeoff_alloc(fit, "x@.overall", "y@.overall", 0.97)
  expect_identical(res$field_status, "optimal")
  expect_lte(res$field_improve_cv, res$requested_cv * (1 + 1e-9))
})
