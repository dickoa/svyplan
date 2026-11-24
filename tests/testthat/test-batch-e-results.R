## Results, displays and messages

test_that("a cluster domain design reports its national cv in $operational", {
  frame <- data.frame(region = c("N", "N", "S", "S"),
                      stratum = c("a", "b", "c", "d"), N = 2000,
                      sd = c(14, 16, 15, 13), mean = c(58, 62, 60, 59),
                      icc_1 = c(0.04, 0.06, 0.05, 0.03),
                      cost_1 = c(400, 550, 450, 500),
                      cost_2 = c(45, 60, 50, 55))
  two <- suppressWarnings(n_alloc(frame, cv = 0.05, alloc = "optimal",
                                  domains = "region"))
  op <- two$operational
  expect_equal(op$cv, op$rmoe / qnorm(0.975), tolerance = 1e-12)
  expect_equal(op$cv, op$se / sum(frame$N * frame$mean / sum(frame$N)),
               tolerance = 1e-12)
  frame3 <- transform(frame, icc_2 = 0.1, cost_3 = 8)
  three <- suppressWarnings(n_alloc(frame3, cv = 0.05, alloc = "optimal",
                                    domains = "region"))
  op3 <- three$operational
  expect_equal(op3$cv, op3$rmoe / qnorm(0.975), tolerance = 1e-12)
})

test_that("$detail flags the continuous binding row, as $domains names it", {
  tab <- data.frame(region = rep(c("N", "S"), each = 2),
                    name = rep(c("a", "b"), 2), p = 0.3,
                    cv = c(0.05, 0.0685, 0.05, 0.06),
                    icc_1 = c(0.02, 0.2, 0.02, 0.2))
  fit <- suppressWarnings(n_cluster(tab, stage_cost = c(500, 50),
                                    domains = "region"))
  expect_false(identical(fit$domains$.binding, fit$operational$domains$.binding))
  for (lev in c("N", "S")) {
    det <- fit$detail[fit$detail$region == lev, ]
    expect_identical(det$name[det$.binding],
                     fit$domains$.binding[fit$domains$region == lev])
  }
})

test_that("the binding domain is the one furthest from its target", {
  tab <- data.frame(region = c("N", "S"), name = c("a", "b"), p = 0.3,
                    cv = c(0.08, 0.03), icc_1 = 0.05)
  fit <- suppressWarnings(n_cluster(tab, stage_cost = c(500, 50),
                                    domains = "region", budget = 40000))
  expect_equal(fit$domains$.cv[1], fit$domains$.cv[2])
  expect_identical(fit$binding, "b")
  expect_identical(fit$operational$binding, "b")
})

test_that("the two-phase print counts phase-1 respondents among eligible units", {
  frame <- data.frame(stratum = paste0("s", 1:4), N = c(4000, 3000, 2000, 1000),
                      mean = c(50, 60, 70, 80), sd = c(10, 12, 15, 20),
                      unit_cost = c(20, 25, 30, 40))
  fit <- n_twophase(frame, phase1_cost = 1, budget = 50000, resp_rate = 0.9,
                    eligibility_rate = 0.8)
  rate <- fit$responding[["n_phase1"]] / fit$n[["n_phase1"]]
  expect_equal(rate, 0.72)
  expected <- round(fit$operational$n[["n_phase1"]] * rate)
  expect_output(print(fit), sprintf("expected responding: n_phase1 = %d ",
                                    expected))
})

test_that("a cluster n_alloc() print names its stages and its PSU total", {
  frame <- data.frame(stratum = c("A", "B"), N = c(6000, 9000), sd = c(10, 14),
                      mean = c(40, 50), icc_1 = c(0.05, 0.1),
                      cost_1 = c(400, 600), cost_2 = c(50, 40))
  two <- n_alloc(frame, cv = 0.03, alloc = "optimal")
  expect_output(print(two), "Stratum allocation \\(optimal, two-stage, 2 strata\\)")
  expect_output(print(two), sprintf("n_1 = %d", sum(two$detail$n_1_int)))
  three <- n_alloc(transform(frame, icc_2 = 0.1, cost_3 = 5), cv = 0.03,
                   alloc = "optimal")
  expect_output(print(three), "three-stage, 2 strata")
  expect_output(print(three), sprintf("n_1 = %d", sum(three$detail$n_1_int)))
})

test_that("prec_cluster() on an indicator table prints the precision columns", {
  tab <- data.frame(region = c("N", "S"), name = c("a", "b"), p = c(0.3, 0.4),
                    cv = 0.05, icc_1 = c(0.05, 0.1))
  fit <- n_cluster(tab, stage_cost = c(500, 50), domains = "region")
  p <- prec_cluster(fit)
  out <- capture.output(print(p))
  header <- out[grepl("\\.se", out)]
  expect_length(header, 1L)
  expect_identical(strsplit(trimws(header), " +")[[1]],
                   c("region", "name", ".se", ".moe", ".rmoe", ".cv"))
  expect_true(all(c("icc_1", "n_1", "take_2") %in% names(p$detail)))
})

test_that("a free-take summary shows both takes and calls them optimized", {
  frame <- data.frame(stratum = c("U", "R"), N = c(1600, 2000),
                      N_1 = c(40, 50), cost_1 = c(110, 190),
                      cost_2 = c(14, 10))
  measures <- data.frame(stratum = c("U", "R"), name = "y", p = c(0.35, 0.45),
                         icc_1 = c(0.07, 0.04))
  fit <- n_alloc_joint(frame, measures, targets = data.frame(name = "y",
                                                              cv = 0.06))
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Optimized later-stage takes", out)))
  expect_false(any(grepl("Fixed later-stage takes", out)))
  stage <- out[grep("Stage decisions", out) + 2L]
  expect_match(stage, "take_2 cont\\.")
  expect_match(stage, "take_2 field")
  fixed <- transform(frame, take_2 = c(12, 10))
  fit_f <- n_alloc_joint(fixed, measures, targets = data.frame(name = "y",
                                                                cv = 0.06))
  expect_true(any(grepl("Fixed later-stage takes",
                        capture.output(print(summary(fit_f))))))
})

test_that("prec_cluster() of a fit at an explicit n is not called the continuous optimum", {
  fit <- n_cluster(c(500, 50), icc = 0.05, cv = 0.05)
  at_n <- prec_cluster(fit, n = c(100, 10))
  expect_null(at_n$params$continuous_optimum)
  expect_false(any(grepl("continuous optimum", capture.output(print(at_n)))))
  expect_true(prec_cluster(fit)$params$continuous_optimum)
})

test_that("stage labels read as words and population sizes carry stage names", {
  fit <- n_cluster(c(ea = 500, hh = 50), icc = 0.05, cv = 0.05,
                   N = c(200, 20000))
  expect_output(print(fit), "corrected by stage for ea = 200 \\| hh = 20000")
  expect_output(print(summary(fit)), "population size: ea = 200, hh = 20000")
  plain <- n_cluster(c(500, 50), icc = 0.05, cv = 0.05, N = c(200, 20000))
  expect_output(print(plain),
                "corrected by stage for stage_1 = 200 \\| stage_2 = 20000")
  set.seed(1)
  dat <- data.frame(cl = rep(1:40, each = 6), y = rnorm(240))
  dat$y <- dat$y + rnorm(40)[dat$cl]
  vc <- varcomp(y ~ cl, data = dat)
  expect_output(print(vc), "SRS stage 1")
  expect_false(any(grepl("stage_1", capture.output(print(vc)))))
})

test_that("the budget frontier of a bounded free take starts at a feasible budget", {
  frame <- data.frame(stratum = c("U", "R"), N = c(1600, 2000),
                      N_1 = c(40, 50), cost_1 = c(110, 190),
                      cost_2 = c(14, 10))
  measures <- data.frame(stratum = rep(c("U", "R"), 2),
                         name = rep(c("p1", "y1"), each = 2),
                         p = c(0.35, 0.45, NA, NA), mean = c(NA, NA, 100, 120),
                         sd = c(NA, NA, 25, 30),
                         icc_1 = c(0.07, 0.04, 0.03, 0.08))
  targets <- data.frame(name = "p1", cv = 0.08)
  fit <- n_alloc_joint(frame, measures, targets = targets, objective = "y1",
                       budget = 8000, take_2_min = 6, take_2_max = 18,
                       take_2_step = 3)
  grid <- svyplan:::.budget_grid(fit, 5)
  floor_fit <- n_alloc_joint(frame, measures, targets = targets,
                             take_2_min = 6, take_2_max = 18, take_2_step = 3)
  expect_equal(grid[1], floor_fit$operational$cost)
  expect_no_error(expect_no_warning(
    n_alloc_joint(frame, measures, targets = targets, objective = "y1",
                  budget = grid[1], take_2_min = 6, take_2_max = 18,
                  take_2_step = 3)
  ))
})

test_that("a register budget that buys a census gives objective 0", {
  frame <- data.frame(stratum = c("A", "B"), N = c(231, 309),
                      cost_1 = c(329, 156), cost_2 = c(7, 3), take_2 = c(2, 2))
  psu <- data.frame(stratum = rep(c("A", "B"), c(6, 5)),
                    N = c(11, 14, 13, 44, 101, 48, 40, 25, 15, 98, 131))
  measures <- data.frame(
    stratum = c("A", "B", "A", "B"), name = c("a", "a", "b", "b"),
    p = c(0.269824, 0.332814, NA, NA), mean = c(NA, NA, 76.9329, 39.4353),
    sd = c(NA, NA, 22.9506, 29.0158),
    icc_1 = c(0.213607, 0.124009, 0.186257, 0.280169)
  )
  fit <- n_alloc_joint(frame, measures, objective = c("a", "b"),
                       budget = 7000, psu = psu)
  expect_equal(fit$operational$objective_value, 0, tolerance = 1e-12)
  expect_lte(fit$operational$cost, 7000)
})

test_that("the summary of a register precision assessment prints", {
  frame <- data.frame(stratum = c("A", "B"), N = c(231, 309),
                      cost_1 = c(329, 156), cost_2 = c(7, 3), take_2 = c(2, 2))
  psu <- data.frame(stratum = rep(c("A", "B"), c(6, 5)),
                    N = c(11, 14, 13, 44, 101, 48, 40, 25, 15, 98, 131))
  measures <- data.frame(
    stratum = c("A", "B", "A", "B"), name = c("a", "a", "b", "b"),
    p = c(0.269824, 0.332814, NA, NA), mean = c(NA, NA, 76.9329, 39.4353),
    sd = c(NA, NA, 22.9506, 29.0158),
    icc_1 = c(0.213607, 0.124009, 0.186257, 0.280169)
  )
  expect_warning(
    fit <- n_alloc_joint(frame, measures, targets = data.frame(name = "a", cv = 0.2),
                         psu = psu),
    class = "svyplan_single_draw"
  )
  out <- expect_no_error(capture.output(print(summary(prec_alloc_joint(fit)))))
  expect_true(any(grepl("Certainty split by stratum", out)))
})

test_that("the precision functions and confint() refuse an allocation", {
  frame <- data.frame(stratum = c("A", "B"), N = c(6000, 9000), sd = c(10, 14),
                      mean = c(40, 50))
  fit <- n_alloc(frame, cv = 0.03)
  types <- c(prec_prop = "proportion", prec_mean = "mean",
             prec_ratio = "ratio", prec_multi = "multi",
             prec_change = "change", prec_pooled = "pooled")
  for (fn in names(types)) {
    expect_error(
      get(fn)(fit),
      sprintf("%s requires a svyplan_n of type '%s'.*prec_alloc\\(\\)",
              fn, types[[fn]]),
      info = fn
    )
  }
  expect_error(confint(fit), "confint not supported for this type.*prec_alloc\\(\\)")
})

test_that("effective_n() takes aux_r2 wherever it builds a design effect", {
  expect_equal(effective_n(n = 1200, icc = 0.05, take = 25, aux_r2 = 0.3),
               1200 / (2.2 * 0.7))
  expect_error(effective_n(n = 1200, deff = 2, aux_r2 = 0.3),
               "either 'deff' or the design components")
  fit <- n_cluster(c(500, 50), icc = 0.05, cv = 0.05)
  expect_equal(effective_n(fit, aux_r2 = 0.3), effective_n(fit) / 0.7,
               tolerance = 1e-12)
  expect_equal(effective_n(prec_cluster(fit), aux_r2 = 0.3),
               effective_n(prec_cluster(fit)) / 0.7, tolerance = 1e-12)
  frame <- data.frame(stratum = c("A", "B"), N = c(6000, 9000),
                      sd = c(10, 14), mean = c(40, 50))
  alloc <- n_alloc(frame, cv = 0.03)
  expect_equal(effective_n(alloc, aux_r2 = 0.3), effective_n(alloc) / 0.7,
               tolerance = 1e-12)
})

test_that("a domain summary flags the field binding beside the field precision", {
  tab <- data.frame(region = rep(c("N", "S"), each = 2),
                    name = rep(c("a", "b"), 2), p = 0.3,
                    cv = c(0.05, 0.0685, 0.05, 0.06),
                    icc_1 = c(0.02, 0.2, 0.02, 0.2))
  fit <- suppressWarnings(n_cluster(tab, stage_cost = c(500, 50),
                                    domains = "region"))
  s <- summary(fit)
  field <- fit$operational$domains$.binding
  expect_identical(s$indicators$name[s$indicators$binding], field)
})
