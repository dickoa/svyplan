test_that("classic scenarios replace eligibility columns only when varied", {
  frame <- data.frame(stratum = c("A", "B"), N = c(1000, 2000),
                      mean = c(40, 50), sd = c(8, 15),
                      eligibility_rate = c(0.8, 0.7),
                      eligible_per_unit = c(1.5, 2))
  fit <- n_alloc(frame, cv = 0.05)
  assessment <- prec_alloc(fit)
  scenarios <- list(eligibility_rate = c(0.6, 1), eligible_per_unit = c(1, 2.5))
  for (key in names(scenarios)) {
    values <- scenarios[[key]]
    grid <- setNames(data.frame(values), key)
    sized <- predict(fit, grid)
    assessed <- predict(assessment, grid)
    edited <- frame
    edited[[key]] <- NULL
    for (i in seq_along(values)) {
      arg <- setNames(list(values[i]), key)
      size <- do.call(n_alloc, c(list(frame = edited, cv = 0.05), arg))
      precision <- do.call(
        prec_alloc, c(list(frame = edited, n = assessment$params$n), arg)
      )
      expect_equal(sized$n[i], size$n)
      expect_equal(assessed$cv[i], precision$cv)
    }
  }
  expect_equal(predict(fit, data.frame(cv = 0.04))$n,
               n_alloc(frame, cv = 0.04)$n)
  expect_equal(predict(assessment, data.frame(alpha = 0.05))$cv, assessment$cv)
})

test_that("predict.svyplan_alloc_joint supports global re-optimization scenarios", {
  frame <- data.frame(
    stratum = c("U", "R"),
    N = c(4000, 6000),
    N_1 = c(80, 100),
    cost_1 = c(120, 200),
    cost_2 = c(15, 12),
    take_2 = c(10, 12)
  )
  measures <- data.frame(
    stratum = rep(c("U", "R"), 2),
    name = rep(c("p1", "y1"), each = 2),
    p = c(0.4, 0.3, NA, NA),
    mean = c(NA, NA, 50, 60),
    sd = c(NA, NA, 15, 20),
    icc_1 = c(0.06, 0.04, 0.03, 0.08)
  )
  targets <- data.frame(name = c("p1", "y1"), cv = c(0.06, 0.045))

  fit <- n_alloc_joint(frame, measures, targets = targets)

  # target_mult: tighter target -> higher cost and n
  g_t <- predict(fit, data.frame(target_mult = c(0.85, 1.0, 1.2)))
  expect_equal(nrow(g_t), 3L)
  expect_equal(g_t$cost[2], fit$params$achieved$cost, tolerance = 1e-8)
  expect_gt(g_t$cost[1], g_t$cost[2])
  expect_lt(g_t$cost[3], g_t$cost[2])

  # cost_mult: scales cost linearly, leaves n unchanged in minimum-cost mode
  g_c <- predict(fit, data.frame(cost_mult = c(1.0, 1.5)))
  expect_equal(g_c$n[1], g_c$n[2], tolerance = 1e-6)
  expect_equal(g_c$cost[2], 1.5 * g_c$cost[1], tolerance = 1e-6)
  expect_equal(g_c$cost_int[2], 1.5 * g_c$cost_int[1], tolerance = 1e-6)

  # icc_mult: higher ICC -> higher sample size and cost
  g_i <- predict(fit, data.frame(icc_mult = c(0.5, 1.0, 1.5)))
  expect_lt(g_i$n[1], g_i$n[2])
  expect_gt(g_i$n[3], g_i$n[2])

  # resp_rate and deff overrides
  g_r <- predict(fit, data.frame(resp_rate = c(0.7, 1.0), deff = c(1.2, 1.0)))
  expect_gt(g_r$n[1], g_r$n[2])
  expect_equal(g_r$n[2], fit$n, tolerance = 1e-6)

  # Free-take fit varies take_2_min / take_2_max
  frame_ft <- frame
  frame_ft$take_2 <- NA_real_
  fit_ft <- n_alloc_joint(frame_ft, measures, targets = targets,
                          take_2_min = 4, take_2_max = 25)
  g_ft <- predict(fit_ft, data.frame(take_2_min = c(4, 15), take_2_max = c(25, 25)))
  expect_equal(nrow(g_ft), 2L)
  expect_true(all(g_ft$.feasible))
  expect_lte(g_ft$cost[1], g_ft$cost[2] + 1e-6)
})

test_that("predict.svyplan_alloc_joint works on 1-stage element fits", {
  frame <- data.frame(stratum = c("A", "B"), N = c(1000, 2000), unit_cost = c(10, 15))
  measures <- data.frame(
    stratum = rep(c("A", "B"), 2),
    name = rep(c("p1", "y1"), each = 2),
    p = c(0.3, 0.5, NA, NA),
    mean = c(NA, NA, 50, 60),
    sd = c(NA, NA, 10, 14)
  )
  targets <- data.frame(name = c("p1", "y1"), cv = c(0.08, 0.04))
  fit <- n_alloc_joint(frame, measures, targets = targets)

  res <- predict(fit, data.frame(target_mult = c(0.9, 1.0, 1.1)))
  expect_equal(nrow(res), 3L)
  expect_equal(res$n[2], fit$n, tolerance = 1e-6)
  expect_gt(res$n[1], res$n[2])
  expect_lt(res$n[3], res$n[2])
})

test_that("predict.svyplan_cluster supports multi-indicator scenario grids", {
  ind <- data.frame(
    name = c("prev", "inc"),
    p = c(0.3, NA),
    mean = c(NA, 100),
    sd = c(NA, 25),
    cv = c(0.08, 0.05),
    icc_1 = c(0.03, 0.06)
  )
  fit <- n_cluster(ind, stage_cost = c(400, 40))

  g <- predict(fit, data.frame(
    target_mult = c(0.9, 1.0, 1.2),
    icc_mult = c(1.0, 1.0, 1.5),
    cost_mult = c(1.0, 1.0, 1.2)
  ))
  expect_equal(nrow(g), 3L)
  expect_equal(g$total_n[2], fit$total_n, tolerance = 1e-6)
  expect_equal(g$cost[2], fit$cost, tolerance = 1e-6)
  expect_gt(g$total_n[1], g$total_n[2])
  expect_true("binding" %in% names(g))
})

test_that("predict() on prec_alloc_joint() runs fixed-design robustness sweeps (including PSU register)", {
  frame <- data.frame(
    stratum = c("U", "R"),
    N = c(4000, 6000),
    N_1 = c(80, 100),
    cost_1 = c(120, 200),
    cost_2 = c(15, 12),
    take_2 = c(10, 12)
  )
  measures <- data.frame(
    stratum = rep(c("U", "R"), 2),
    name = rep(c("p1", "y1"), each = 2),
    p = c(0.4, 0.3, NA, NA),
    mean = c(NA, NA, 50, 60),
    sd = c(NA, NA, 15, 20),
    icc_1 = c(0.06, 0.04, 0.03, 0.08)
  )
  targets <- data.frame(name = c("p1", "y1"), cv = c(0.06, 0.045))
  fit <- n_alloc_joint(frame, measures, targets = targets)
  pr <- prec_alloc_joint(fit, n = fit$detail$n_int)

  rob <- predict(pr, data.frame(
    icc_mult = c(1.0, 1.5, 2.0),
    var_mult = c(1.0, 1.0, 1.25),
    cost_mult = c(1.0, 1.1, 1.2)
  ))
  expect_equal(nrow(rob), 3L)
  # Sample size is held FIXED across robustness scenarios
  expect_equal(rob$n, rep(sum(fit$detail$n_int), 3L))
  expect_true(rob$.feasible[1])
  expect_gt(rob$.max_ratio[2], rob$.max_ratio[1])
  expect_gt(rob$.max_ratio[3], rob$.max_ratio[2])
  expect_equal(rob$cost, pr$params$achieved$cost * c(1.0, 1.1, 1.2),
               tolerance = 1e-6)

  psu <- data.frame(
    psu_id = sprintf("P%02d", 1:20),
    stratum = rep(c("U", "R"), each = 10),
    N = c(c(400, rep(100, 9)), c(350, rep(90, 9)))
  )
  f_psu <- data.frame(
    stratum = c("U", "R"),
    N = c(sum(psu$N[psu$stratum == "U"]), sum(psu$N[psu$stratum == "R"])),
    cost_1 = c(150, 200),
    cost_2 = c(12, 10),
    take_2 = c(10, 10)
  )
  fit_psu <- n_alloc_joint(f_psu, measures, targets = targets, psu = psu)
  pr_psu <- prec_alloc_joint(fit_psu)

  rob_psu <- predict(pr_psu, data.frame(
    icc_mult = c(1.0, 1.5),
    var_mult = c(1.0, 1.2),
    cost_mult = c(1.0, 1.25),
    resp_rate = c(1.0, 0.85)
  ))
  expect_equal(nrow(rob_psu), 2L)
  expect_equal(rob_psu$n[1], rob_psu$n[2])
  expect_true(rob_psu$.feasible[1])
  expect_gt(rob_psu$.max_ratio[2], rob_psu$.max_ratio[1])
  expect_equal(rob_psu$cost[2], 1.25 * rob_psu$cost[1], tolerance = 1e-6)

  # Varying n on a PSU register precision result is refused clearly
  expect_error(
    predict(pr_psu, data.frame(n = c(300, 400))),
    "cannot vary 'n'"
  )
})

test_that("predict() validates icc_mult overflow cleanly", {
  frame <- data.frame(
    stratum = "A", N = 2000, N_1 = 50,
    cost_1 = 100, cost_2 = 10, take_2 = 10
  )
  measures <- data.frame(stratum = "A", name = "y", mean = 50, sd = 10, icc_1 = 0.6)
  targets <- data.frame(name = "y", cv = 0.05)
  fit <- n_alloc_joint(frame, measures, targets = targets)
  pr <- prec_alloc_joint(fit, n = fit$detail$n_int)

  expect_error(predict(pr, data.frame(icc_mult = 2)), "exceeds 1")
  expect_error(predict(fit, data.frame(icc_mult = 2)), "exceeds 1")
})

## Domain fits and stated targets reach the refit

# The table a user would pass for a target multiple of `m`.
scale_targets <- function(tbl, m) {
  for (col in intersect(c("cv", "moe", "rmoe"), names(tbl))) {
    tbl[[col]] <- tbl[[col]] * m
  }
  tbl
}

multi_domain_tables <- function() {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.1), moe = c(0.05, 0.03))
  rep_names <- rbind(transform(ind, region = "N"),
                     transform(ind, region = "S", p = c(0.4, 0.2)))
  uniq <- rep_names
  uniq$name <- paste0(uniq$name, uniq$region)
  two <- rbind(transform(rep_names, sex = "F"),
               transform(rep_names, sex = "M", p = rep_names$p + 0.05))
  list(rep_names = rep_names, uniq = uniq, two = two)
}

test_that("predict() on an n_multi() domain fit keeps its domains", {
  tbls <- multi_domain_tables()
  cases <- list(
    list(indicators = tbls$rep_names, domains = "region"),
    list(indicators = tbls$uniq, domains = "region"),
    list(indicators = tbls$rep_names, domains = "region", min_n_domain = 800),
    list(indicators = tbls$two, domains = c("region", "sex")),
    list(indicators = transform(tbls$rep_names, share = rep(c(0.3, 0.7), each = 2)),
         domains = "region", domain_sampling = "natural")
  )
  for (cs in cases) {
    fit <- do.call(n_multi, cs)
    g <- predict(fit, data.frame(target_mult = c(1, 1.2)))
    expect_equal(g$n[1], fit$n, tolerance = 1e-10)
    cs$indicators <- scale_targets(cs$indicators, 1.2)
    expect_equal(g$n[2], do.call(n_multi, cs)$n, tolerance = 1e-10)
  }
})

test_that("predict() on an n_cluster() domain fit keeps its domains", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.1), cv = c(0.1, 0.15),
                    icc_1 = c(0.02, 0.05))
  rep_names <- rbind(transform(ind, region = "N"),
                     transform(ind, region = "S", p = c(0.4, 0.2)))
  uniq <- rep_names
  uniq$name <- paste0(uniq$name, uniq$region)
  pops <- transform(rep_names, N_1 = rep(c(60, 90), each = 2),
                    N = rep(c(6000, 9000), each = 2))
  cases <- list(
    list(x = rep_names),
    list(x = uniq),
    list(x = rep_names, budget = 60000),
    list(x = rep_names, allocation = "joint", budget = 60000),
    list(x = pops),
    list(x = rbind(transform(rep_names, sex = "F"),
                   transform(rep_names, sex = "M")),
         domains = c("region", "sex"))
  )
  for (cs in cases) {
    cs$stage_cost <- c(500, 50)
    cs$domains <- cs$domains %||% "region"
    fit <- do.call(n_cluster, cs)
    g <- predict(fit, data.frame(target_mult = c(1, 1.2)))
    expect_equal(g$total_n[1], fit$total_n, tolerance = 1e-10)
    expect_equal(g$cost[1], fit$cost, tolerance = 1e-10)
    cs$x <- scale_targets(cs$x, 1.2)
    ref <- do.call(n_cluster, cs)
    expect_equal(g$total_n[2], ref$total_n, tolerance = 1e-10)
    expect_equal(g$cost[2], ref$cost, tolerance = 1e-10)
  }
})

test_that("predict() on n_cluster() refits from the targets as stated", {
  base <- data.frame(name = c("a", "b"), p = c(0.3, 0.1),
                     icc_1 = c(0.02, 0.05))
  tbls <- list(
    transform(base, moe = c(0.05, 0.03)),
    transform(base, rmoe = c(0.2, 0.35)),
    transform(base, moe = c(0.05, NA), cv = c(NA, 0.15)),
    transform(base, moe = c(0.05, 0.03), resp_rate = c(0.8, 0.9))
  )
  for (tbl in tbls) {
    for (er in c(1, 0.9)) {
      fit <- n_cluster(tbl, stage_cost = c(500, 50), eligibility_rate = er)
      g <- predict(fit, data.frame(target_mult = c(1, 1.5)))
      expect_equal(g$total_n[1], fit$total_n, tolerance = 1e-10)
      ref <- n_cluster(scale_targets(tbl, 1.5), stage_cost = c(500, 50),
                       eligibility_rate = er)
      expect_equal(g$total_n[2], ref$total_n, tolerance = 1e-10)
      expect_equal(g$cost[2], ref$cost, tolerance = 1e-10)
    }
  }
})

test_that("predict() refuses a domain fit whose domain columns are lost", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.1), moe = c(0.05, 0.03))
  tbl <- rbind(transform(ind, region = "N"), transform(ind, region = "S"))
  fit <- n_multi(tbl, domains = "region")
  fit$params$domain_cols <- character(0)
  expect_error(predict(fit, data.frame(target_mult = 1)),
               "records no domain columns")
  cl <- n_cluster(transform(tbl, icc_1 = 0.03), stage_cost = c(500, 50),
                  domains = "region")
  cl$params$domain_cols <- NULL
  expect_error(predict(cl, data.frame(target_mult = 1)),
               "records no domain columns")
})

## A continuous joint precision result predicts at its own sizes

joint_fixed_take <- function(stages) {
  if (stages == 2L) {
    frame <- data.frame(stratum = c("A", "B"), N = c(4000, 6000),
                        N_1 = c(80, 100), cost_1 = c(120, 200),
                        cost_2 = c(15, 12), take_2 = c(10, 12))
    measures <- data.frame(stratum = c("A", "B"), name = "y", p = c(0.4, 0.3),
                           icc_1 = c(0.06, 0.04))
  } else {
    frame <- data.frame(stratum = c("A", "B"), N = c(40000, 60000),
                        N_1 = c(80, 100), N_2 = c(800, 1000),
                        cost_1 = c(300, 400), cost_2 = c(60, 50),
                        cost_3 = c(15, 12), take_2 = c(4, 5),
                        take_3 = c(10, 12))
    measures <- data.frame(stratum = c("A", "B"), name = "y", p = c(0.4, 0.3),
                           icc_1 = c(0.06, 0.04), icc_2 = c(0.1, 0.08))
  }
  n_alloc_joint(frame, measures, targets = data.frame(name = "y", cv = 0.06))
}

test_that("predict() on a continuous fixed-take joint precision result", {
  for (stages in 2:3) {
    fit <- joint_fixed_take(stages)
    pr <- prec_alloc_joint(fit)
    expect_false(all(pr$params$n == round(pr$params$n)))
    g <- predict(pr, data.frame(var_mult = c(1, 1), icc_mult = c(1, 1.2)))
    expect_equal(g$n[1], pr$params$achieved$n, tolerance = 1e-10)
    expect_equal(g$cost[1], pr$params$achieved$cost, tolerance = 1e-10)
    expect_equal(g$.max_ratio[1], max(pr$detail$.ratio), tolerance = 1e-10)
    scaled <- fit
    for (col in grep("^icc_", names(scaled$params$measures), value = TRUE)) {
      scaled$params$measures[[col]] <- scaled$params$measures[[col]] * 1.2
    }
    pr_s <- prec_alloc_joint(scaled)
    expect_equal(g$n[2], g$n[1], tolerance = 1e-10)
    expect_equal(g$.max_ratio[2], max(pr_s$detail$.ratio), tolerance = 1e-10)
    expect_gt(g$.max_ratio[2], g$.max_ratio[1])
  }
})

test_that("an explicit fractional PSU count is still refused", {
  fit <- joint_fixed_take(2L)
  pr <- prec_alloc_joint(fit, n = fit$detail$n_int)
  g <- predict(pr, data.frame(var_mult = 1))
  expect_equal(g$n, sum(fit$detail$n_int))
  expect_warning(
    expect_error(predict(pr, data.frame(n = I(list(c(100.5, 240))))),
                 "all rows failed"),
    "whole number of PSUs"
  )
  expect_error(
    prec_alloc_joint(fit$params$frame, n = c(100.5, 240),
                     measures = fit$params$measures,
                     targets = fit$params$targets),
    "whole number of PSUs"
  )
})

## A scenario icc moves the var_ratio it implies

test_that("predict() re-derives a default var_ratio from the scenario icc", {
  sc3 <- c(500, 100, 20)
  sc4 <- c(800, 200, 50, 10)
  check <- function(fit, row, ...) {
    g <- predict(fit, row)
    ref <- n_cluster(...)
    expect_equal(g$total_n, ref$total_n, tolerance = 1e-10)
    expect_equal(g$cost, ref$cost, tolerance = 1e-10)
    expect_equal(g$cv, ref$cv, tolerance = 1e-10)
  }
  c3 <- n_cluster(stage_cost = sc3, icc = c(0.05, 0.10), cv = 0.05)
  check(c3, data.frame(icc_1 = 0.08),
        stage_cost = sc3, icc = c(0.08, 0.10), cv = 0.05)
  check(c3, data.frame(var_ratio_2 = 0.7),
        stage_cost = sc3, icc = c(0.05, 0.10), var_ratio = c(1, 0.7),
        cv = 0.05)
  c4 <- n_cluster(stage_cost = sc4, icc = c(0.05, 0.10, 0.2), cv = 0.05)
  check(c4, data.frame(icc_1 = 0.08, icc_2 = 0.15),
        stage_cost = sc4, icc = c(0.08, 0.15, 0.2), cv = 0.05)
  c3s <- n_cluster(stage_cost = sc3, icc = c(0.05, 0.10), var_ratio = 1.3,
                   cv = 0.05)
  check(c3s, data.frame(icc_1 = 0.08),
        stage_cost = sc3, icc = c(0.08, 0.10), var_ratio = 1.3, cv = 0.05)
})

test_that("predict() keeps a var_ratio the user gave", {
  sc3 <- c(500, 100, 20)
  fit <- n_cluster(stage_cost = sc3, icc = c(0.05, 0.10),
                   var_ratio = c(1, 0.7), cv = 0.05)
  g <- predict(fit, data.frame(icc_1 = 0.08))
  ref <- n_cluster(stage_cost = sc3, icc = c(0.08, 0.10),
                   var_ratio = c(1, 0.7), cv = 0.05)
  expect_equal(g$total_n, ref$total_n, tolerance = 1e-10)
})

test_that("icc_mult on a three-stage indicator fit equals the direct refit", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.1), cv = c(0.1, 0.15),
                    icc_1 = c(0.02, 0.05), icc_2 = c(0.1, 0.08))
  fit <- n_cluster(ind, stage_cost = c(500, 100, 20))
  g <- predict(fit, data.frame(icc_mult = c(1, 1.5)))
  expect_equal(g$total_n[1], fit$total_n, tolerance = 1e-10)
  scaled <- transform(ind, icc_1 = icc_1 * 1.5, icc_2 = icc_2 * 1.5)
  ref <- n_cluster(scaled, stage_cost = c(500, 100, 20))
  expect_equal(g$total_n[2], ref$total_n, tolerance = 1e-10)
  expect_equal(g$cost[2], ref$cost, tolerance = 1e-10)
})

test_that("a refit at another icc from a precision result re-derives var_ratio", {
  sc3 <- c(500, 100, 20)
  c3 <- n_cluster(stage_cost = sc3, icc = c(0.05, 0.10), cv = 0.05)
  ref <- n_cluster(stage_cost = sc3, icc = c(0.08, 0.10), cv = 0.05)
  expect_equal(n_cluster(prec_cluster(c3), icc = c(0.08, 0.10))$total_n,
               ref$total_n, tolerance = 1e-10)
  direct <- prec_cluster(n = c3$n, stage_cost = sc3, icc = c(0.08, 0.10))
  expect_equal(prec_cluster(c3, icc = c(0.08, 0.10))$cv, direct$cv,
               tolerance = 1e-12)
  kept <- n_cluster(stage_cost = sc3, icc = c(0.05, 0.10),
                    var_ratio = c(1, 0.7), cv = 0.05)
  ref_k <- n_cluster(stage_cost = sc3, icc = c(0.08, 0.10),
                     var_ratio = c(1, 0.7), cv = 0.05)
  expect_equal(n_cluster(prec_cluster(kept), icc = c(0.08, 0.10))$total_n,
               ref_k$total_n, tolerance = 1e-10)
})

test_that("an NA var_ratio past the first entry is derived", {
  sc4 <- c(800, 200, 50, 10)
  icc <- c(0.05, 0.10, 0.2)
  k <- c(1.2, 1.2 * 0.95, 0.7)
  a <- n_cluster(stage_cost = sc4, icc = icc, var_ratio = c(1.2, NA, 0.7),
                 cv = 0.05)
  b <- n_cluster(stage_cost = sc4, icc = icc, var_ratio = k, cv = 0.05)
  expect_equal(a$total_n, b$total_n, tolerance = 1e-12)
  expect_equal(unname(a$params$var_ratio), k)
  pa <- prec_cluster(n = b$n, icc = icc, var_ratio = c(1.2, NA, 0.7))
  expect_equal(pa$cv, prec_cluster(n = b$n, icc = icc, var_ratio = k)$cv)
  expect_error(n_cluster(stage_cost = sc4, icc = icc,
                         var_ratio = c(NA, 1, 0.7), cv = 0.05),
               "NA only past the first")
  expect_error(prec_cluster(n = b$n, icc = icc, var_ratio = c(1, -1, NA)),
               "positive finite")
})

## Scenario grids keep earlier-stage rates, column types and per-row checks

test_that("a varied resp_rate replaces only the final stage of a multi-indicator fit", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.4), cv = 0.1, icc_1 = 0.05)
  fit <- n_cluster(ind, stage_cost = c(100, 10), resp_rate = c(0.95, 0.80))
  g <- predict(fit, data.frame(resp_rate = c(0.80, 0.70)))
  expect_equal(g$total_n[1], fit$total_n)
  direct <- n_cluster(ind, stage_cost = c(100, 10), resp_rate = c(0.95, 0.70))
  expect_equal(g$total_n[2], direct$total_n)
  expect_equal(g$cost[2], direct$cost)
})

test_that("a failed row keeps the column types of the good rows", {
  fit <- n_alloc_joint(data.frame(stratum = c("A", "B"), N = c(1000, 2000)),
                       measures = data.frame(stratum = c("A", "B"), name = "y",
                                             mean = 10, var = c(4, 9)),
                       targets = data.frame(name = "y", cv = 0.05))
  pr <- prec_alloc_joint(fit)
  for (rr in list(c(0.8, -1), c(-1, 0.8))) {
    g <- suppressWarnings(predict(pr, data.frame(resp_rate = rr)))
    expect_type(g$.feasible, "logical")
    expect_type(g$.n_fail, "integer")
    expect_true(is.na(g$.feasible[rr < 0]))
  }
  expect_error(suppressWarnings(predict(pr, data.frame(resp_rate = c(-1, -2)))),
               "all rows failed")
})

test_that("an issued overlap fails per row and a bad column type is refused first", {
  x <- n_pooled(var = 100, moe = 1, occasions = 8,
                overlap = design_overlap("4", max_lag = 7), cor_decay = 0.8)
  expect_warning(g <- predict(x, data.frame(resp_rate = c(1, 0.8))),
                 "predict row 2 failed: design_overlap\\(\\) of a rotation")
  expect_equal(g$n[1], x$n)
  expect_true(is.na(g$n[2]))
  expect_error(predict(x, data.frame(resp_rate = c("a", "b"))),
               "non-numeric columns in newdata")
  pr <- prec_pooled(x)
  expect_warning(g <- predict(pr, data.frame(resp_rate = c(1, 0.8))),
                 "predict row 2 failed")
  expect_equal(g$se[1], pr$se)
})

test_that("a joint fit varies eligible_per_unit", {
  fr <- data.frame(stratum = c("A", "B"), N = c(1000, 2000))
  me <- data.frame(stratum = c("A", "B"), name = "y", p = c(0.2, 0.4))
  tg <- data.frame(name = "y", cv = 0.05)
  g <- predict(n_alloc_joint(fr, me, targets = tg),
               data.frame(eligible_per_unit = c(1, 2)))
  expect_equal(g$n[2], sum(n_alloc_joint(fr, me, targets = tg,
                                         eligible_per_unit = 2)$detail$n))
})
