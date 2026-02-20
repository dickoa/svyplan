test_that("n_multi rejects non-data-frame", {
  expect_error(n_multi(list(p = 0.3)), "non-empty data frame")
  expect_error(n_multi(data.frame()), "non-empty data frame")
})

test_that("n_multi requires p or var column", {
  expect_error(n_multi(data.frame(moe = 0.05)), "must contain 'p' or 'var'")
})

test_that("n_multi requires moe or cv column", {
  expect_error(n_multi(data.frame(p = 0.3)), "must contain 'moe' or 'cv'")
})

test_that("n_multi rejects p outside (0,1)", {
  expect_error(n_multi(data.frame(p = 0, moe = 0.05)), "must be in \\(0, 1\\)")
  expect_error(n_multi(data.frame(p = 1, moe = 0.05)), "must be in \\(0, 1\\)")
})

test_that("n_multi rejects rows with both p and var set", {
  df <- data.frame(p = 0.3, var = 10, moe = 0.05)
  expect_error(n_multi(df), "only one of")
})

test_that("n_multi rejects budget without cost", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(n_multi(df, budget = 10000), "'budget' requires 'cost'")
})

test_that("n_multi rejects m without cost", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(n_multi(df, m = 50), "'m' requires 'cost'")
})

test_that("multistage requires cv column", {
  df <- data.frame(p = 0.3, moe = 0.05, delta1 = 0.02)
  expect_error(n_multi(df, cost = c(500, 50)), "requires 'cv'")
})

test_that("multistage requires delta1 column", {
  df <- data.frame(p = 0.3, cv = 0.10)
  expect_error(n_multi(df, cost = c(500, 50)), "requires 'delta1'")
})

test_that("var + cv requires mu", {
  df <- data.frame(var = 100, cv = 0.05, delta1 = 0.02)
  expect_error(n_multi(df, cost = c(500, 50)), "'mu' is required")
})

test_that("n_multi rejects negative cv in multistage", {
  expect_error(
    n_multi(data.frame(p = 0.3, cv = -0.1, delta1 = 0.05), cost = c(500, 50)),
    "positive"
  )
})

test_that("n_multi rejects delta1 outside [0, 1]", {
  expect_error(
    n_multi(data.frame(p = 0.3, cv = 0.1, delta1 = 1.5), cost = c(500, 50)),
    "\\[0, 1\\]"
  )
})

test_that("n_multi rejects negative moe", {
  expect_error(
    n_multi(data.frame(p = 0.3, moe = -0.05)),
    "positive"
  )
})

test_that("single proportion indicator matches n_prop", {
  df <- data.frame(p = 0.3, moe = 0.05)
  res <- n_multi(df)
  ref <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$type, "multi")
})

test_that("single mean indicator matches n_mean", {
  df <- data.frame(var = 100, moe = 2)
  res <- n_multi(df)
  ref <- n_mean(var = 100, moe = 2)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("single proportion CV mode matches n_prop", {
  df <- data.frame(p = 0.5, cv = 0.10)
  res <- n_multi(df)
  ref <- n_prop(p = 0.5, cv = 0.10)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("single mean CV mode matches n_mean", {
  df <- data.frame(var = 100, mu = 50, cv = 0.05)
  res <- n_multi(df)
  ref <- n_mean(var = 100, mu = 50, cv = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("multiple proportions: max is selected", {
  df <- data.frame(
    name = c("stunting", "vaccination", "anemia"),
    p = c(0.30, 0.70, 0.10),
    moe = c(0.05, 0.05, 0.03)
  )
  res <- n_multi(df)

  # Verify each individual n
  n1 <- n_prop(p = 0.30, moe = 0.05)$n
  n2 <- n_prop(p = 0.70, moe = 0.05)$n
  n3 <- n_prop(p = 0.10, moe = 0.03)$n

  expect_equal(res$n, max(n1, n2, n3), tolerance = 1e-6)
  expect_equal(res$binding, "anemia")
  expect_true(res$detail$.binding[3])
  expect_false(res$detail$.binding[1])
})

test_that("mixed proportion + mean indicators work", {
  df <- data.frame(
    name = c("prop_ind", "mean_ind"),
    p = c(0.5, NA),
    var = c(NA, 100),
    moe = c(0.05, 2)
  )
  res <- n_multi(df)

  n1 <- n_prop(p = 0.5, moe = 0.05)$n
  n2 <- n_mean(var = 100, moe = 2)$n

  expect_equal(res$n, max(n1, n2), tolerance = 1e-6)
})

test_that("deff is applied per indicator", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.5, 0.3),
    moe = c(0.05, 0.05),
    deff = c(1, 2)
  )
  res <- n_multi(df)
  ref_a <- n_prop(p = 0.5, moe = 0.05, deff = 1)$n
  ref_b <- n_prop(p = 0.3, moe = 0.05, deff = 2)$n
  expect_equal(res$n, max(ref_a, ref_b), tolerance = 1e-6)
})

test_that("FPC is applied per indicator", {
  df <- data.frame(
    p = c(0.5, 0.5),
    moe = c(0.05, 0.05),
    N = c(Inf, 500)
  )
  res <- n_multi(df)
  n_inf <- n_prop(p = 0.5, moe = 0.05)$n
  n_500 <- n_prop(p = 0.5, moe = 0.05, N = 500)$n
  expect_equal(res$n, max(n_inf, n_500), tolerance = 1e-6)
})

test_that("alpha defaults to 0.05 when not specified", {
  df <- data.frame(p = 0.3, moe = 0.05)
  res <- n_multi(df)
  ref <- n_prop(p = 0.3, moe = 0.05, alpha = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("custom alpha is respected", {
  df <- data.frame(p = 0.3, moe = 0.05, alpha = 0.10)
  res <- n_multi(df)
  ref <- n_prop(p = 0.3, moe = 0.05, alpha = 0.10)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("detail dataframe has correct structure", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05)
  )
  res <- n_multi(df)
  expect_true(is.data.frame(res$detail))
  expect_true(".n" %in% names(res$detail))
  expect_true(".binding" %in% names(res$detail))
  expect_equal(nrow(res$detail), 2L)
  expect_equal(sum(res$detail$.binding), 1L)
})

test_that("domain column is auto-detected", {
  df <- data.frame(
    p = c(0.3, 0.3),
    moe = c(0.05, 0.05),
    region = c("North", "South")
  )
  res <- n_multi(df)
  expect_s3_class(res, "svyplan_n")
  expect_true(!is.null(res$domains))
  expect_equal(nrow(res$domains), 2L)
  expect_true("region" %in% names(res$domains))
})

test_that("per-domain results differ when parameters differ", {
  df <- data.frame(
    name = c("anemia", "anemia"),
    p = c(0.10, 0.30),
    moe = c(0.03, 0.05),
    region = c("North", "South")
  )
  res <- n_multi(df)
  expect_true(res$domains$.n[1] != res$domains$.n[2])
})

test_that("multi-indicator multi-domain takes max across domains", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.3, 0.4, 0.1, 0.2),
    moe = rep(0.05, 4),
    region = rep(c("R1", "R2"), 2)
  )
  res <- n_multi(df)
  expect_equal(res$n, max(res$domains$.n), tolerance = 1e-6)
})

test_that("domains with two grouping columns work", {
  df <- data.frame(
    p = rep(0.3, 4),
    moe = rep(0.05, 4),
    region = c("N", "N", "S", "S"),
    urban = c("U", "R", "U", "R")
  )
  res <- n_multi(df)
  expect_equal(nrow(res$domains), 4L)
})

test_that("single indicator 2-stage matches n_cluster", {
  df <- data.frame(
    p = 0.3,
    cv = 0.05,
    delta1 = 0.05
  )
  res <- n_multi(df, cost = c(500, 50))
  ref <- n_cluster(
    cost = c(500, 50),
    delta = 0.05,
    rel_var = (1 - 0.3) / 0.3,
    cv = 0.05
  )
  expect_equal(res$n[["n2"]], ref$n[["n2"]], tolerance = 0.5)
  expect_equal(res$cv, ref$cv, tolerance = 0.01)
  expect_s3_class(res, "svyplan_cluster")
})

test_that("multi-indicator 2-stage has correct structure", {
  df <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50))
  expect_s3_class(res, "svyplan_cluster")
  expect_equal(res$stages, 2L)
  expect_true(!is.null(res$detail))
  expect_true(!is.null(res$binding))
  expect_equal(nrow(res$detail), 2L)
  expect_true(".cv_achieved" %in% names(res$detail))
  expect_true(".cv_target" %in% names(res$detail))
})

test_that("multi-indicator 2-stage achieves CV for binding indicator", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50))
  # Binding indicator should approximately meet its target
  binding_row <- which(res$detail$.binding)
  expect_equal(
    res$detail$.cv_achieved[binding_row],
    res$detail$.cv_target[binding_row],
    tolerance = 0.01
  )
})

test_that("round-trip: cv_cluster confirms achieved CVs", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50))

  for (j in 1:2) {
    rv <- (1 - df$p[j]) / df$p[j]
    cv_check <- unname(cv_cluster(
      n = res$n,
      delta = df$delta1[j],
      rel_var = rv
    ))
    expect_equal(cv_check, res$detail$.cv_achieved[j], tolerance = 1e-4)
  }
})

test_that("2-stage budget mode respects budget", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000)
  expect_equal(res$cost, 100000, tolerance = 1)
})

test_that("2-stage budget mode returns achieved CVs", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000)
  expect_true(all(res$detail$.cv_achieved > 0))
})

test_that("2-stage budget + fixed m works", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000, m = 80)
  expect_equal(res$n[["n1"]], 80)
  expect_true(res$cost <= 100000 + 1)
})

test_that("2-stage budget too small raises error", {
  df <- data.frame(
    p = 0.3,
    cv = 0.10,
    delta1 = 0.05
  )
  expect_error(
    n_multi(df, cost = c(500, 50), budget = 100, m = 80),
    "budget is too small"
  )
})

test_that("2-stage with domains produces domain results", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.3, 0.4, 0.1, 0.2),
    cv = rep(0.10, 4),
    delta1 = rep(0.02, 4),
    region = rep(c("North", "South"), 2)
  )
  res <- n_multi(df, cost = c(500, 50))
  expect_s3_class(res, "svyplan_cluster")
  expect_true(!is.null(res$domains))
  expect_equal(nrow(res$domains), 2L)
  expect_true("n1" %in% names(res$domains))
  expect_true("n2" %in% names(res$domains))
  expect_true(".binding" %in% names(res$domains))
})

test_that("3-stage multi-indicator has correct structure", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.01, 0.02),
    delta2 = c(0.05, 0.08)
  )
  res <- n_multi(df, cost = c(500, 100, 50))
  expect_s3_class(res, "svyplan_cluster")
  expect_equal(res$stages, 3L)
  expect_equal(length(res$n), 3L)
  expect_true(!is.null(res$detail))
})

test_that("3-stage budget mode works", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.01, 0.02),
    delta2 = c(0.05, 0.08)
  )
  res <- n_multi(df, cost = c(500, 100, 50), budget = 200000)
  expect_equal(res$cost, 200000, tolerance = 1)
  expect_true(all(res$detail$.cv_achieved > 0))
})

test_that("identical targets yield same per-indicator n", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.3),
    moe = c(0.05, 0.05)
  )
  res <- n_multi(df)
  expect_equal(res$detail$.n[1], res$detail$.n[2], tolerance = 1e-6)
})

test_that("name column is optional", {
  df <- data.frame(p = c(0.3, 0.5), moe = c(0.05, 0.05))
  res <- n_multi(df)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$detail$name, c(1L, 2L))
})

test_that("rel_var can be supplied directly", {
  df <- data.frame(
    p = 0.3,
    cv = 0.10,
    delta1 = 0.05,
    rel_var = 3.0
  )
  res <- n_multi(df, cost = c(500, 50))
  # Should use supplied rel_var, not derived (1-0.3)/0.3 = 2.33
  ref <- n_cluster(cost = c(500, 50), delta = 0.05, rel_var = 3.0, cv = 0.10)
  expect_equal(res$n[["n2"]], ref$n[["n2"]], tolerance = 0.5)
})

test_that("print.svyplan_n works for multi type", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05)
  )
  res <- n_multi(df)
  expect_output(print(res), "Multi-indicator sample size")
  expect_output(print(res), "binding")
})

test_that("print.svyplan_n works for multi type with domains", {
  df <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05),
    region = c("A", "B")
  )
  res <- n_multi(df)
  expect_output(print(res), "Multi-indicator sample size")
  expect_output(print(res), "domains")
})

test_that("print.svyplan_cluster works for multi type", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- n_multi(df, cost = c(500, 50))
  expect_output(print(res), "Multi-indicator optimal allocation")
  expect_output(print(res), "binding")
})

test_that("print.svyplan_cluster works for multi type with domains", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- n_multi(df, cost = c(500, 50))
  expect_output(print(res), "Multi-indicator optimal allocation")
  expect_output(print(res), "domains")
})


test_that("invalid joint value raises error", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(n_multi(df, joint = "yes"), "TRUE or FALSE")
  expect_error(n_multi(df, joint = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(n_multi(df, joint = NA), "TRUE or FALSE")
})

test_that("joint budget: worst CV ratio <= equal-split (asymmetric domains)", {
  df_hard <- data.frame(
    name = c("a", "b"),
    p = c(0.05, 0.10),
    cv = c(0.20, 0.15),
    delta1 = c(0.08, 0.05)
  )
  df_easy <- data.frame(
    name = c("a", "b"),
    p = c(0.40, 0.30),
    cv = c(0.10, 0.10),
    delta1 = c(0.02, 0.02)
  )
  budget <- 100000
  cost <- c(500, 50)

  res_h_eq <- n_multi(df_hard, cost = cost, budget = budget / 2)
  res_e_eq <- n_multi(df_easy, cost = cost, budget = budget / 2)
  worst_equal <- max(
    max(res_h_eq$detail$.cv_achieved / res_h_eq$detail$.cv_target),
    max(res_e_eq$detail$.cv_achieved / res_e_eq$detail$.cv_target)
  )

  df <- rbind(
    cbind(df_hard, region = "Hard"),
    cbind(df_easy, region = "Easy")
  )
  res_jnt <- n_multi(df, cost = cost, budget = budget, joint = TRUE)
  hard_b <- res_jnt$domains$.cost[res_jnt$domains$region == "Hard"]
  easy_b <- res_jnt$domains$.cost[res_jnt$domains$region == "Easy"]

  res_h_jnt <- n_multi(df_hard, cost = cost, budget = hard_b)
  res_e_jnt <- n_multi(df_easy, cost = cost, budget = easy_b)
  worst_joint <- max(
    max(res_h_jnt$detail$.cv_achieved / res_h_jnt$detail$.cv_target),
    max(res_e_jnt$detail$.cv_achieved / res_e_jnt$detail$.cv_target)
  )

  expect_true(worst_joint <= worst_equal + 1e-4)
})

test_that("joint budget: single domain identical to non-joint", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("R1", "R1")
  )
  res_ind <- n_multi(df, cost = c(500, 50), budget = 80000, joint = FALSE)
  res_jnt <- n_multi(df, cost = c(500, 50), budget = 80000, joint = TRUE)
  expect_equal(res_ind$domains$.cv, res_jnt$domains$.cv, tolerance = 1e-4)
  expect_equal(res_ind$cost, res_jnt$cost, tolerance = 1)
})

test_that("joint budget: equal domains get approximately equal budgets", {
  df <- data.frame(
    p = c(0.3, 0.3),
    cv = c(0.10, 0.10),
    delta1 = c(0.03, 0.03),
    region = c("A", "B")
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  budgets <- res$domains$.cost
  expect_equal(budgets[1], budgets[2], tolerance = budgets[1] * 0.05)
})

test_that("joint budget: doubling budget improves worst CV", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.30, 0.10, 0.10, 0.30),
    cv = c(0.10, 0.15, 0.15, 0.10),
    delta1 = c(0.02, 0.05, 0.05, 0.02),
    region = rep(c("R1", "R2"), 2)
  )
  res1 <- n_multi(df, cost = c(500, 50), budget = 50000, joint = TRUE)
  res2 <- n_multi(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  expect_true(max(res2$domains$.cv) < max(res1$domains$.cv))
})

test_that("joint budget: total cost equals budget", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  expect_equal(sum(res$domains$.cost), 100000, tolerance = 1)
})

test_that("joint budget: budget shifts toward harder domain", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.05, 0.40, 0.10, 0.30),
    cv = c(0.20, 0.10, 0.15, 0.10),
    delta1 = c(0.08, 0.02, 0.05, 0.02),
    region = rep(c("Hard", "Easy"), 2)
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  hard_budget <- res$domains$.cost[res$domains$region == "Hard"]
  easy_budget <- res$domains$.cost[res$domains$region == "Easy"]
  expect_true(hard_budget > easy_budget)
})

test_that("joint budget: 2-domain grid search confirms optimizer", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.10, 0.30, 0.20, 0.40),
    cv = c(0.12, 0.10, 0.10, 0.08),
    delta1 = c(0.05, 0.03, 0.04, 0.02),
    region = rep(c("R1", "R2"), 2)
  )
  budget <- 80000
  cost <- c(500, 50)

  res_jnt <- n_multi(df, cost = cost, budget = budget, joint = TRUE)
  worst_jnt <- max(
    res_jnt$domains$.cv /
      c(
        max(df$cv[df$region == "R1"]),
        max(df$cv[df$region == "R2"])
      )
  )

  grid_fracs <- seq(0.05, 0.95, by = 0.01)
  grid_worst <- vapply(
    grid_fracs,
    function(w1) {
      r1 <- n_multi(
        df[df$region == "R1", names(df) != "region"],
        cost = cost,
        budget = w1 * budget
      )
      r2 <- n_multi(
        df[df$region == "R2", names(df) != "region"],
        cost = cost,
        budget = (1 - w1) * budget
      )
      max(
        max(r1$detail$.cv_achieved / r1$detail$.cv_target),
        max(r2$detail$.cv_achieved / r2$detail$.cv_target)
      )
    },
    numeric(1L)
  )
  grid_best <- min(grid_worst)

  expect_equal(worst_jnt, grid_best, tolerance = 0.01)
})

test_that("joint budget: output has correct structure", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  expect_s3_class(res, "svyplan_cluster")
  expect_true(!is.null(res$domains))
  expect_true("region" %in% names(res$domains))
  expect_true(".cost" %in% names(res$domains))
  expect_true(".cv" %in% names(res$domains))
  expect_true(isTRUE(res$params$joint))
})

test_that("joint = TRUE without budget is same as independent (CV mode)", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res_ind <- n_multi(df, cost = c(500, 50), joint = FALSE)
  res_jnt <- n_multi(df, cost = c(500, 50), joint = TRUE)
  expect_equal(res_ind$domains$.cv, res_jnt$domains$.cv, tolerance = 1e-4)
})

test_that("joint = TRUE without domains is same as joint = FALSE", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res_ind <- n_multi(df, cost = c(500, 50), budget = 80000, joint = FALSE)
  res_jnt <- n_multi(df, cost = c(500, 50), budget = 80000, joint = TRUE)
  expect_equal(res_ind$cv, res_jnt$cv, tolerance = 1e-4)
  expect_equal(res_ind$cost, res_jnt$cost, tolerance = 1)
})

test_that("joint = TRUE in simple mode is ignored", {
  df <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05)
  )
  res_ind <- n_multi(df, joint = FALSE)
  res_jnt <- n_multi(df, joint = TRUE)
  expect_equal(res_ind$n, res_jnt$n)
})

test_that("joint budget: 3+ domains work", {
  df <- data.frame(
    p = c(0.3, 0.1, 0.2),
    cv = c(0.10, 0.15, 0.12),
    delta1 = c(0.02, 0.05, 0.03),
    region = c("A", "B", "C")
  )
  res <- n_multi(df, cost = c(500, 50), budget = 150000, joint = TRUE)
  expect_s3_class(res, "svyplan_cluster")
  expect_equal(nrow(res$domains), 3L)
  expect_equal(sum(res$domains$.cost), 150000, tolerance = 1)
  expect_true(isTRUE(res$params$joint))
})

test_that("joint budget: 3-stage design works", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.01, 0.02),
    delta2 = c(0.05, 0.08),
    region = c("A", "B")
  )
  res <- n_multi(df, cost = c(500, 100, 50), budget = 200000, joint = TRUE)
  expect_s3_class(res, "svyplan_cluster")
  expect_equal(res$stages, 3L)
  expect_equal(sum(res$domains$.cost), 200000, tolerance = 1)
  expect_true(isTRUE(res$params$joint))
})

test_that("joint budget: multiple indicators per domain", {
  df <- data.frame(
    name = rep(c("stunting", "anemia", "vaccination"), each = 2),
    p = c(0.30, 0.25, 0.10, 0.15, 0.70, 0.60),
    cv = rep(c(0.10, 0.15, 0.08), each = 2),
    delta1 = rep(c(0.02, 0.05, 0.03), each = 2),
    region = rep(c("Urban", "Rural"), 3)
  )
  res <- n_multi(df, cost = c(500, 50), budget = 120000, joint = TRUE)
  expect_equal(nrow(res$domains), 2L)
  expect_equal(sum(res$domains$.cost), 120000, tolerance = 1)
})

test_that("joint budget: print shows joint label", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- n_multi(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  expect_output(print(res), "joint")
})
