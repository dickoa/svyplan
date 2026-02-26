nm <- function(...) suppressMessages(n_multi(...))

test_that("n_multi rejects non-data-frame", {
  expect_error(nm(list(p = 0.3)), "non-empty data frame")
  expect_error(nm(data.frame()), "non-empty data frame")
})

test_that("n_multi requires p or var column", {
  expect_error(nm(data.frame(moe = 0.05)), "must contain 'p' or 'var'")
})

test_that("n_multi requires moe or cv column", {
  expect_error(nm(data.frame(p = 0.3)), "must contain 'moe' or 'cv'")
})

test_that("n_multi rejects p outside (0,1)", {
  expect_error(nm(data.frame(p = 0, moe = 0.05)), "must be in \\(0, 1\\)")
  expect_error(nm(data.frame(p = 1, moe = 0.05)), "must be in \\(0, 1\\)")
})

test_that("n_multi rejects rows with both p and var set", {
  df <- data.frame(p = 0.3, var = 10, moe = 0.05)
  expect_error(nm(df), "only one of")
})

test_that("n_multi rejects budget without cost", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, budget = 10000), "'budget' requires 'cost'")
})

test_that("n_multi rejects m without cost", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, m = 50), "'m' requires 'cost'")
})

test_that("multistage requires cv column", {
  df <- data.frame(p = 0.3, moe = 0.05, delta1 = 0.02)
  expect_error(nm(df, cost = c(500, 50)), "requires 'cv'")
})

test_that("multistage requires delta1 column", {
  df <- data.frame(p = 0.3, cv = 0.10)
  expect_error(nm(df, cost = c(500, 50)), "requires 'delta1'")
})

test_that("var + cv requires mu", {
  df <- data.frame(var = 100, cv = 0.05, delta1 = 0.02)
  expect_error(nm(df, cost = c(500, 50)), "'mu' is required")
})

test_that("n_multi rejects negative cv in multistage", {
  expect_error(
    nm(data.frame(p = 0.3, cv = -0.1, delta1 = 0.05), cost = c(500, 50)),
    "positive"
  )
})

test_that("n_multi rejects delta1 outside [0, 1]", {
  expect_error(
    nm(data.frame(p = 0.3, cv = 0.1, delta1 = 1.5), cost = c(500, 50)),
    "\\[0, 1\\]"
  )
})

test_that("n_multi rejects negative moe", {
  expect_error(
    nm(data.frame(p = 0.3, moe = -0.05)),
    "positive"
  )
})

test_that("single proportion indicator matches n_prop", {
  df <- data.frame(p = 0.3, moe = 0.05)
  res <- nm(df)
  ref <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$type, "multi")
})

test_that("single mean indicator matches n_mean", {
  df <- data.frame(var = 100, moe = 2)
  res <- nm(df)
  ref <- n_mean(var = 100, moe = 2)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("single proportion CV mode matches n_prop", {
  df <- data.frame(p = 0.5, cv = 0.10)
  res <- nm(df)
  ref <- n_prop(p = 0.5, cv = 0.10)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("single mean CV mode matches n_mean", {
  df <- data.frame(var = 100, mu = 50, cv = 0.05)
  res <- nm(df)
  ref <- n_mean(var = 100, mu = 50, cv = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("multiple proportions: max is selected", {
  df <- data.frame(
    name = c("stunting", "vaccination", "anemia"),
    p = c(0.30, 0.70, 0.10),
    moe = c(0.05, 0.05, 0.03)
  )
  res <- nm(df)

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
  res <- nm(df)

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
  res <- nm(df)
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
  res <- nm(df)
  n_inf <- n_prop(p = 0.5, moe = 0.05)$n
  n_500 <- n_prop(p = 0.5, moe = 0.05, N = 500)$n
  expect_equal(res$n, max(n_inf, n_500), tolerance = 1e-6)
})

test_that("alpha defaults to 0.05 when not specified", {
  df <- data.frame(p = 0.3, moe = 0.05)
  res <- nm(df)
  ref <- n_prop(p = 0.3, moe = 0.05, alpha = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("custom alpha is respected", {
  df <- data.frame(p = 0.3, moe = 0.05, alpha = 0.10)
  res <- nm(df)
  ref <- n_prop(p = 0.3, moe = 0.05, alpha = 0.10)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("detail dataframe has correct structure", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05)
  )
  res <- nm(df)
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
  res <- nm(df)
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
  res <- nm(df)
  expect_true(res$domains$.n[1] != res$domains$.n[2])
})

test_that("multi-indicator multi-domain takes max across domains", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.3, 0.4, 0.1, 0.2),
    moe = rep(0.05, 4),
    region = rep(c("R1", "R2"), 2)
  )
  res <- nm(df)
  expect_equal(res$n, max(res$domains$.n), tolerance = 1e-6)
})

test_that("domains with two grouping columns work", {
  df <- data.frame(
    p = rep(0.3, 4),
    moe = rep(0.05, 4),
    region = c("N", "N", "S", "S"),
    urban = c("U", "R", "U", "R")
  )
  res <- nm(df)
  expect_equal(nrow(res$domains), 4L)
})

test_that("single indicator 2-stage matches n_cluster", {
  df <- data.frame(
    p = 0.3,
    cv = 0.05,
    delta1 = 0.05
  )
  res <- nm(df, cost = c(500, 50))
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
  res <- nm(df, cost = c(500, 50))
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
  res <- nm(df, cost = c(500, 50))
  # Binding indicator should approximately meet its target
  binding_row <- which(res$detail$.binding)
  expect_equal(
    res$detail$.cv_achieved[binding_row],
    res$detail$.cv_target[binding_row],
    tolerance = 0.01
  )
})

test_that("round-trip: prec_cluster confirms achieved CVs", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- nm(df, cost = c(500, 50))

  for (j in 1:2) {
    rv <- (1 - df$p[j]) / df$p[j]
    cv_check <- unname(
      prec_cluster(
        n = res$n,
        delta = df$delta1[j],
        rel_var = rv
      )$cv
    )
    expect_equal(cv_check, res$detail$.cv_achieved[j], tolerance = 1e-4)
  }
})

test_that("2-stage budget mode respects budget", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- nm(df, cost = c(500, 50), budget = 100000)
  expect_equal(res$cost, 100000, tolerance = 1)
})

test_that("2-stage budget mode returns achieved CVs", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- nm(df, cost = c(500, 50), budget = 100000)
  expect_true(all(res$detail$.cv_achieved > 0))
})

test_that("2-stage budget + fixed m works", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res <- nm(df, cost = c(500, 50), budget = 100000, m = 80)
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
    nm(df, cost = c(500, 50), budget = 100, m = 80),
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
  res <- nm(df, cost = c(500, 50))
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
  res <- nm(df, cost = c(500, 100, 50))
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
  res <- nm(df, cost = c(500, 100, 50), budget = 200000)
  expect_equal(res$cost, 200000, tolerance = 1)
  expect_true(all(res$detail$.cv_achieved > 0))
})

test_that("identical targets yield same per-indicator n", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.3),
    moe = c(0.05, 0.05)
  )
  res <- nm(df)
  expect_equal(res$detail$.n[1], res$detail$.n[2], tolerance = 1e-6)
})

test_that("name column is optional", {
  df <- data.frame(p = c(0.3, 0.5), moe = c(0.05, 0.05))
  res <- nm(df)
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
  res <- nm(df, cost = c(500, 50))
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
  res <- nm(df)
  expect_output(print(res), "Multi-indicator sample size")
  expect_output(print(res), "binding")
})

test_that("print.svyplan_n works for multi type with domains", {
  df <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05),
    region = c("A", "B")
  )
  res <- nm(df)
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
  res <- nm(df, cost = c(500, 50))
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
  res <- nm(df, cost = c(500, 50))
  expect_output(print(res), "Multi-indicator optimal allocation")
  expect_output(print(res), "domains")
})


test_that("invalid joint value raises error", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, joint = "yes"), "TRUE or FALSE")
  expect_error(nm(df, joint = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(nm(df, joint = NA), "TRUE or FALSE")
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

  res_h_eq <- nm(df_hard, cost = cost, budget = budget / 2)
  res_e_eq <- nm(df_easy, cost = cost, budget = budget / 2)
  worst_equal <- max(
    max(res_h_eq$detail$.cv_achieved / res_h_eq$detail$.cv_target),
    max(res_e_eq$detail$.cv_achieved / res_e_eq$detail$.cv_target)
  )

  df <- rbind(
    cbind(df_hard, region = "Hard"),
    cbind(df_easy, region = "Easy")
  )
  res_jnt <- nm(df, cost = cost, budget = budget, joint = TRUE)
  hard_b <- res_jnt$domains$.cost[res_jnt$domains$region == "Hard"]
  easy_b <- res_jnt$domains$.cost[res_jnt$domains$region == "Easy"]

  res_h_jnt <- nm(df_hard, cost = cost, budget = hard_b)
  res_e_jnt <- nm(df_easy, cost = cost, budget = easy_b)
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
  res_ind <- nm(df, cost = c(500, 50), budget = 80000, joint = FALSE)
  res_jnt <- nm(df, cost = c(500, 50), budget = 80000, joint = TRUE)
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
  res <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
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
  res1 <- nm(df, cost = c(500, 50), budget = 50000, joint = TRUE)
  res2 <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  expect_true(max(res2$domains$.cv) < max(res1$domains$.cv))
})

test_that("joint budget: total cost equals budget", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
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
  res <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
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

  res_jnt <- nm(df, cost = cost, budget = budget, joint = TRUE)
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
      r1 <- nm(
        df[df$region == "R1", names(df) != "region"],
        cost = cost,
        budget = w1 * budget
      )
      r2 <- nm(
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
  res <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
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
  res_ind <- nm(df, cost = c(500, 50), joint = FALSE)
  res_jnt <- nm(df, cost = c(500, 50), joint = TRUE)
  expect_equal(res_ind$domains$.cv, res_jnt$domains$.cv, tolerance = 1e-4)
})

test_that("joint = TRUE without domains is same as joint = FALSE", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res_ind <- nm(df, cost = c(500, 50), budget = 80000, joint = FALSE)
  res_jnt <- nm(df, cost = c(500, 50), budget = 80000, joint = TRUE)
  expect_equal(res_ind$cv, res_jnt$cv, tolerance = 1e-4)
  expect_equal(res_ind$cost, res_jnt$cost, tolerance = 1)
})

test_that("joint = TRUE in simple mode is ignored", {
  df <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05)
  )
  res_ind <- nm(df, joint = FALSE)
  res_jnt <- nm(df, joint = TRUE)
  expect_equal(res_ind$n, res_jnt$n)
})

test_that("joint budget: 3+ domains work", {
  df <- data.frame(
    p = c(0.3, 0.1, 0.2),
    cv = c(0.10, 0.15, 0.12),
    delta1 = c(0.02, 0.05, 0.03),
    region = c("A", "B", "C")
  )
  res <- suppressWarnings(nm(df, cost = c(500, 50), budget = 150000, joint = TRUE))
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
  res <- nm(df, cost = c(500, 100, 50), budget = 200000, joint = TRUE)
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
  res <- nm(df, cost = c(500, 50), budget = 120000, joint = TRUE)
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
  res <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  expect_output(print(res), "joint")
})

test_that("resp_rate defaults to 1 when not specified", {
  df <- data.frame(p = 0.3, moe = 0.05)
  res <- nm(df)
  ref <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$n, ref$n, tolerance = 1e-6)
})

test_that("simple mode resp_rate inflates n per indicator", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05),
    resp_rate = c(0.8, 0.9)
  )
  res <- nm(df)
  n_a <- n_prop(p = 0.3, moe = 0.05)$n / 0.8
  n_b <- n_prop(p = 0.5, moe = 0.05)$n / 0.9
  expect_equal(res$n, max(n_a, n_b), tolerance = 1e-6)
})

test_that("resp_rate is not misdetected as a domain column", {
  df <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05),
    resp_rate = c(0.8, 0.9)
  )
  res <- nm(df)
  expect_null(res$domains)
})

test_that("2-stage resp_rate: per-indicator inflation in optimization", {
  df_equal <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    resp_rate = c(1, 1)
  )
  df_diff <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    resp_rate = c(0.5, 0.9)
  )
  res_eq <- nm(df_equal, cost = c(500, 50))
  res_diff <- nm(df_diff, cost = c(500, 50))
  expect_true(res_diff$n[["n1"]] > res_eq$n[["n1"]])
})

test_that("2-stage resp_rate: per-indicator vs global min gives different results", {
  # When indicators have very different resp_rates, per-indicator inflation

  # should give a smaller n1 than using the global minimum resp_rate.
  # Indicator A: needs n1_eff=80, resp_rate=0.5 -> n1_actual=160
  # Indicator B: needs n1_eff=100, resp_rate=0.9 -> n1_actual=111
  # Correct: max(160, 111) = 160
  # Wrong (global min): max(80, 100)/0.5 = 200
  df <- data.frame(
    name = c("easy", "hard"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.08),
    delta1 = c(0.02, 0.05),
    resp_rate = c(0.5, 0.9)
  )
  res <- nm(df, cost = c(500, 50))
  # Verify the resp_rate-adjusted CV still meets target for binding indicator
  binding_row <- which(res$detail$.binding)
  expect_lte(
    res$detail$.cv_achieved[binding_row],
    res$detail$.cv_target[binding_row] * 1.01
  )
})

test_that("2-stage budget mode: resp_rate affects CV achieved", {
  df_no_rr <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  df_rr <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    resp_rate = c(0.8, 0.7)
  )
  res_no <- nm(df_no_rr, cost = c(500, 50), budget = 100000)
  res_rr <- nm(df_rr, cost = c(500, 50), budget = 100000)
  # Lower resp_rate means worse precision for same budget
  expect_true(max(res_rr$detail$.cv_achieved) > max(res_no$detail$.cv_achieved))
})

test_that("3-stage resp_rate: per-indicator inflation works", {
  df <- data.frame(
    name = c("a", "b"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.01, 0.02),
    delta2 = c(0.05, 0.08),
    resp_rate = c(0.8, 0.7)
  )
  res <- nm(df, cost = c(500, 100, 50))
  expect_s3_class(res, "svyplan_cluster")
  expect_equal(res$stages, 3L)
  binding_row <- which(res$detail$.binding)
  expect_lte(
    res$detail$.cv_achieved[binding_row],
    res$detail$.cv_target[binding_row] * 1.01
  )
})

test_that("multistage params store budget and m", {
  df <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  res_cv <- nm(df, cost = c(500, 50))
  expect_equal(res_cv$params$cost, c(500, 50))
  expect_null(res_cv$params$budget)
  expect_null(res_cv$params$m)

  res_b <- nm(df, cost = c(500, 50), budget = 100000)
  expect_equal(res_b$params$budget, 100000)

  res_bm <- nm(df, cost = c(500, 50), budget = 100000, m = 80)
  expect_equal(res_bm$params$budget, 100000)
  expect_equal(res_bm$params$m, 80)
})

test_that("domain multistage params store budget", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- nm(df, cost = c(500, 50), budget = 100000)
  expect_equal(res$params$budget, 100000)
})

test_that("uniform resp_rate matches no resp_rate", {
  df1 <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  df2 <- data.frame(
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    resp_rate = c(1, 1)
  )
  res1 <- nm(df1, cost = c(500, 50))
  res2 <- nm(df2, cost = c(500, 50))
  expect_equal(res1$n, res2$n, tolerance = 1e-6)
  expect_equal(res1$cv, res2$cv, tolerance = 1e-6)
})

test_that("min_n validation: rejects non-numeric", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, min_n = "100"), "positive numeric scalar")
})

test_that("min_n validation: rejects negative", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, min_n = -10), "positive numeric scalar")
})

test_that("min_n validation: rejects NA", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, min_n = NA_real_), "positive numeric scalar")
})

test_that("min_n validation: rejects length > 1", {
  df <- data.frame(p = 0.3, moe = 0.05)
  expect_error(nm(df, min_n = c(100, 200)), "positive numeric scalar")
})

test_that("min_n: no domains -> silently ignored", {
  df <- data.frame(p = 0.3, moe = 0.05)
  res_no <- nm(df)
  res_mn <- nm(df, min_n = 9999)
  expect_equal(res_no$n, res_mn$n, tolerance = 1e-6)
})

test_that("min_n: simple + domains floor applied", {
  df <- data.frame(
    name = rep("ind1", 2),
    p = c(0.50, 0.30),
    moe = c(0.05, 0.05),
    region = c("Easy", "Hard")
  )
  res_no <- nm(df)
  easy_n <- res_no$domains$.n[res_no$domains$region == "Easy"]
  hard_n <- res_no$domains$.n[res_no$domains$region == "Hard"]

  floor_val <- ceiling(max(easy_n, hard_n)) + 100
  res_mn <- nm(df, min_n = floor_val)
  expect_true(all(res_mn$domains$.n >= floor_val))
})

test_that("min_n: simple + domains .binding updated for floored domains", {
  df <- data.frame(
    p = c(0.50, 0.05),
    moe = c(0.05, 0.03),
    region = c("A", "B")
  )
  res_no <- nm(df)
  n_B <- res_no$domains$.n[res_no$domains$region == "B"]
  n_A <- res_no$domains$.n[res_no$domains$region == "A"]
  # Floor should be > A but < B so only A is floored
  floor_val <- ceiling(max(n_A, n_B)) + 50
  res_mn <- nm(df, min_n = floor_val)
  expect_true(all(res_mn$domains$.binding == "(min_n)"))
})

test_that("min_n: simple + domains, floor below all domains -> no effect", {
  df <- data.frame(
    p = c(0.30, 0.10),
    moe = c(0.05, 0.03),
    region = c("A", "B")
  )
  res_no <- nm(df)
  res_mn <- nm(df, min_n = 1)
  expect_equal(res_no$domains$.n, res_mn$domains$.n)
  expect_true(all(res_mn$domains$.binding != "(min_n)"))
})

test_that("min_n: simple + domains stores min_n in params", {
  df <- data.frame(
    p = c(0.30, 0.10),
    moe = c(0.05, 0.03),
    region = c("A", "B")
  )
  res <- nm(df, min_n = 500)
  expect_equal(res$params$min_n, 500)
})

test_that("min_n: joint budget ensures all domains >= min_n", {
  df <- data.frame(
    name = rep(c("stunting", "anemia"), each = 2),
    p = c(0.30, 0.25, 0.10, 0.15),
    cv = c(0.10, 0.10, 0.15, 0.15),
    delta1 = c(0.02, 0.03, 0.05, 0.04),
    region = rep(c("Urban", "Rural"), 2)
  )
  res_no <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE)
  min_total <- min(res_no$domains$.total_n)
  floor_val <- ceiling(min_total * 0.5)

  res_mn <- nm(
    df,
    cost = c(500, 50),
    budget = 100000,
    joint = TRUE,
    min_n = floor_val
  )
  expect_true(all(res_mn$domains$.total_n >= floor_val - 1))
  expect_equal(res_mn$params$min_n, floor_val)
})

test_that("min_n: joint monotonicity -> larger min_n -> worse overall CV", {
  df <- data.frame(
    name = rep(c("a", "b"), each = 2),
    p = c(0.05, 0.40, 0.10, 0.30),
    cv = c(0.20, 0.10, 0.15, 0.10),
    delta1 = c(0.08, 0.02, 0.05, 0.02),
    region = rep(c("Hard", "Easy"), 2)
  )
  res1 <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE, min_n = 50)
  res2 <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE, min_n = 500)
  expect_true(max(res2$domains$.cv) >= max(res1$domains$.cv) - 1e-4)
})

test_that("min_n: joint feasibility error when single domain impossible", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  expect_error(
    nm(df, cost = c(500, 50), budget = 5000, joint = TRUE, min_n = 999999),
    "not achievable for domain"
  )
})

test_that("min_n: joint feasibility error when budget too small for all domains", {
  df <- data.frame(
    p = c(0.3, 0.1, 0.2),
    cv = c(0.10, 0.15, 0.12),
    delta1 = c(0.02, 0.05, 0.03),
    region = c("A", "B", "C")
  )
  res_full <- suppressWarnings(nm(
    df,
    cost = c(500, 50),
    budget = 150000,
    joint = TRUE
  ))
  big_floor <- ceiling(max(res_full$domains$.total_n) * 0.9)
  expect_error(
    nm(df, cost = c(500, 50), budget = 150000, joint = TRUE, min_n = big_floor),
    "not achievable"
  )
})

test_that("min_n: non-joint multistage warns when total_n < min_n", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- nm(df, cost = c(500, 50))
  big_floor <- ceiling(max(res$domains$.total_n)) + 1000
  expect_warning(
    nm(df, cost = c(500, 50), min_n = big_floor),
    "total_n below min_n"
  )
})

test_that("min_n: non-joint multistage no warning when all above floor", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  expect_no_warning(suppressMessages(nm(df, cost = c(500, 50), min_n = 1)))
})

test_that("min_n: multistage stores min_n in params", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- suppressWarnings(nm(df, cost = c(500, 50), min_n = 500))
  expect_equal(res$params$min_n, 500)
})

test_that("min_n: print shows min_n for simple domains", {
  df <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.05),
    region = c("A", "B")
  )
  res <- nm(df, min_n = 500)
  expect_output(print(res), "min_n = 500")
})

test_that("min_n: print shows min_n for cluster domains", {
  df <- data.frame(
    p = c(0.3, 0.1),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05),
    region = c("A", "B")
  )
  res <- nm(df, cost = c(500, 50), budget = 100000, joint = TRUE, min_n = 50)
  expect_output(print(res), "min_n = 50")
})

test_that("n_multi rejects negative var", {
  df <- data.frame(var = -10, mu = 5, moe = 1)
  expect_error(nm(df), "var.*positive")
  df2 <- data.frame(var = 0, mu = 5, moe = 1)
  expect_error(nm(df2), "var.*positive")
})

test_that("n_multi rejects non-positive mu", {
  df <- data.frame(var = 10, mu = 0, cv = 0.1)
  expect_error(nm(df), "mu.*positive")
  df2 <- data.frame(var = 10, mu = -5, cv = 0.1)
  expect_error(nm(df2), "mu.*positive")
})

test_that("n_multi rejects N = 1 in simple mode", {
  df <- data.frame(p = 0.3, moe = 0.05, N = 1)
  expect_error(nm(df), "greater than 1")
})
