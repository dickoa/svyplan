test_that("planning and precision functions reject unused arguments", {
  frame <- data.frame(N = c(100, 200), sd = c(10, 15))
  targets_n <- data.frame(p = 0.3, moe = 0.05)
  targets_p <- data.frame(p = 0.3, n = 100)

  expect_error(n_prop(0.3, moe = 0.05, defff = 2), "unused argument.*defff")
  expect_error(n_mean(100, moe = 2, defff = 2), "unused argument.*defff")
  expect_error(
    n_cluster(c(500, 50), delta = 0.05, budget = 1e5, resp_rte = 0.8),
    "unused argument.*resp_rte"
  )
  expect_error(n_multi(targets_n, prop_methd = "wald"),
               "unused argument.*prop_methd")
  expect_error(n_alloc(frame, n = 50, aloc = "optimal"),
               "unused argument.*aloc")

  expect_error(prec_prop(0.3, n = 100, alpa = 0.1),
               "unused argument.*alpa")
  expect_error(prec_mean(100, n = 100, alpa = 0.1),
               "unused argument.*alpa")
  expect_error(
    prec_cluster(c(20, 10), delta = 0.05, resp_rte = 0.8),
    "unused argument.*resp_rte"
  )
  expect_error(prec_multi(targets_p, prop_methd = "wald"),
               "unused argument.*prop_methd")
  expect_error(prec_alloc(frame, n = c(20, 30), resp_rte = 0.8),
               "unused argument.*resp_rte")
})

test_that("power and design-effect functions reject unused arguments", {
  expect_error(power_prop(0.3, p2 = 0.4, powr = 0.9),
               "unused argument.*powr")
  expect_error(power_mean(1, effect = 0.1, powr = 0.9),
               "unused argument.*powr")
  expect_error(
    power_did(c(1, 2), c(1, 1), outcome = "mean", var = 1,
              effect = 0.1, powr = 0.9),
    "unused argument.*powr"
  )
  expect_error(design_effect(1:5, methd = "kish"),
               "unused argument.*methd")
  expect_error(effective_n(1:5, methd = "kish"),
               "unused argument.*methd")
})

test_that("varcomp rejects unused arguments but retains its survey formula", {
  d <- data.frame(y = 1:12, psu = rep(1:3, each = 4))

  expect_error(varcomp(d$y, stage_id = list(d$psu), weigths = rep(1, 12)),
               "unused argument.*weigths")
  expect_error(varcomp(y ~ psu, data = d, weigths = ~w),
               "unused argument.*weigths")

  skip_if_not_installed("survey")
  design <- survey::svydesign(ids = ~psu, data = d, weights = ~1)
  expect_s3_class(varcomp(design, ~y), "svyplan_varcomp")
  expect_error(varcomp(design, ~y, mystery = 1),
               "unused argument.*mystery")
  expect_error(varcomp(design, ~y, ~psu), "exactly one outcome formula")
})

test_that("plan dispatch does not drop unused arguments", {
  plan <- svyplan(alpha = 0.1, prop_method = "wilson")
  targets <- data.frame(p = 0.3, n = 100)

  expect_error(n_prop(0.3, moe = 0.05, plan = plan, defff = 2),
               "unused argument.*defff")
  expect_error(prec_multi(targets, plan = plan, prop_methd = "wald"),
               "unused argument.*prop_methd")
})

test_that("round-trip overrides remain supported and validated", {
  n_result <- n_prop(0.3, moe = 0.05)

  expect_s3_class(prec_prop(n_result, alpha = 0.1), "svyplan_prec")
  expect_error(prec_prop(n_result, alpa = 0.1), "unused argument.*alpa")

  targets <- data.frame(p = 0.3, cv = 0.1, delta_psu = 0.05)
  cluster_result <- n_multi_cluster(targets, stage_cost = c(500, 50))
  expect_error(prec_multi_cluster(cluster_result, alpa = 0.1),
               "unused argument.*alpa")
})

test_that("prediction and display methods reject unused arguments", {
  result <- n_prop(0.3, moe = 0.05)

  expect_error(predict(result, data.frame(p = 0.2), na_rm = TRUE),
               "unused argument.*na_rm")
  expect_error(format(result, digts = 3), "unused argument.*digts")
  expect_error(as.integer(result, roundng = "up"),
               "unused argument.*roundng")
  expect_s3_class(data.frame(result), "data.frame")
  expect_s3_class(as.data.frame(result, optional = TRUE), "data.frame")
})
