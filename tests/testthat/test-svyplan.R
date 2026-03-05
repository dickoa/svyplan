test_that("svyplan constructor creates valid object", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  expect_s3_class(plan, "svyplan")
  expect_equal(plan$defaults$deff, 1.5)
  expect_equal(plan$defaults$resp_rate, 0.9)
})

test_that("svyplan empty constructor", {
  plan <- svyplan()
  expect_s3_class(plan, "svyplan")
  expect_length(plan$defaults, 0L)
})

test_that("svyplan with all core params", {
  plan <- svyplan(alpha = 0.10, N = 50000, deff = 1.8, resp_rate = 0.85)
  expect_equal(plan$defaults$alpha, 0.10)
  expect_equal(plan$defaults$N, 50000)
  expect_equal(plan$defaults$deff, 1.8)
  expect_equal(plan$defaults$resp_rate, 0.85)
})

test_that("svyplan with cluster context", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.85)
  expect_equal(plan$defaults$stage_cost, c(500, 50))
  expect_equal(plan$defaults$delta, 0.05)
  expect_equal(plan$defaults$resp_rate, 0.85)
})

test_that("svyplan with extended defaults", {
  plan <- svyplan(alternative = "one.sided", method = "arcsine",
                  overlap = 0.5, rho = 0.6)
  expect_equal(plan$defaults$alternative, "one.sided")
  expect_equal(plan$defaults$method, "arcsine")
  expect_equal(plan$defaults$overlap, 0.5)
  expect_equal(plan$defaults$rho, 0.6)
})

test_that("svyplan validates core typed params", {
  expect_error(svyplan(alpha = 2), "alpha")
  expect_error(svyplan(alpha = 0), "alpha")
  expect_error(svyplan(deff = -1), "deff")
  expect_error(svyplan(resp_rate = 0), "resp_rate")
  expect_error(svyplan(resp_rate = 1.5), "resp_rate")
  expect_error(svyplan(N = 0), "'N' must be greater than 1")
  expect_error(svyplan(stage_cost = 10), "'stage_cost' must be a numeric vector")
})

test_that("svyplan rejects unknown defaults", {
  expect_error(svyplan(p = 0.3), "unknown default")
  expect_error(svyplan(var = 100), "unknown default")
  expect_error(svyplan(n = 400), "unknown default")
  expect_error(svyplan(effect = 5), "unknown default")
})

test_that("svyplan rejects unnamed arguments", {
  expect_error(svyplan(0.05), "all arguments.*must be named")
})

test_that("svyplan rejects duplicate names", {
  expect_error(do.call(svyplan, list(deff = 1, deff = 2)), "unique")
})

test_that("print.svyplan produces output", {
  plan <- svyplan(deff = 1.5, N = 10000)
  out <- capture.output(print(plan))
  expect_true(any(grepl("svyplan profile", out)))
  expect_true(any(grepl("deff", out)))
  expect_true(any(grepl("10000", out)))
})

test_that("print.svyplan empty plan", {
  plan <- svyplan()
  out <- capture.output(print(plan))
  expect_true(any(grepl("no defaults", out)))
})

test_that("print.svyplan formats vector values", {
  plan <- svyplan(stage_cost = c(500, 50))
  out <- capture.output(print(plan))
  expect_true(any(grepl("c\\(500, 50\\)", out)))
})

test_that("update.svyplan modifies parameters", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  plan2 <- update(plan, deff = 2.0)
  expect_equal(plan2$defaults$deff, 2.0)
  expect_equal(plan2$defaults$resp_rate, 0.9)
})

test_that("update.svyplan adds new parameters", {
  plan <- svyplan(deff = 1.5)
  plan2 <- update(plan, alpha = 0.10)
  expect_equal(plan2$defaults$deff, 1.5)
  expect_equal(plan2$defaults$alpha, 0.10)
})

test_that("update.svyplan validates modified params", {
  plan <- svyplan()
  expect_error(update(plan, alpha = 5), "alpha")
})

test_that("n_prop uses plan defaults", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9, N = 10000)
  res <- n_prop(p = 0.3, moe = 0.05, plan = plan)
  ref <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.9, N = 10000)
  expect_equal(res$n, ref$n)
  expect_equal(res$params$deff, 1.5)
  expect_equal(res$params$resp_rate, 0.9)
  expect_equal(res$params$N, 10000)
})

test_that("n_prop explicit override wins over plan", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- n_prop(p = 0.3, moe = 0.05, plan = plan, deff = 2.0)
  ref <- n_prop(p = 0.3, moe = 0.05, deff = 2.0, resp_rate = 0.9)
  expect_equal(res$n, ref$n)
  expect_equal(res$params$deff, 2.0)
})

test_that("n_prop without plan is unchanged", {
  res <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$params$deff, 1)
  expect_equal(res$params$resp_rate, 1)
  expect_equal(res$params$N, Inf)
})

test_that("n_prop ignores irrelevant plan defaults", {
  plan <- svyplan(deff = 1.5, stage_cost = c(500, 50), delta = 0.05)
  res <- n_prop(p = 0.3, moe = 0.05, plan = plan)
  ref <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
  expect_equal(res$n, ref$n)
})

test_that("n_mean uses plan defaults", {
  plan <- svyplan(deff = 1.8, N = 50000)
  res <- n_mean(var = 100, moe = 2, plan = plan)
  ref <- n_mean(var = 100, moe = 2, deff = 1.8, N = 50000)
  expect_equal(res$n, ref$n)
})

test_that("n_mean explicit override wins over plan", {
  plan <- svyplan(deff = 1.8)
  res <- n_mean(var = 100, moe = 2, plan = plan, deff = 1.0)
  ref <- n_mean(var = 100, moe = 2, deff = 1.0)
  expect_equal(res$n, ref$n)
})

test_that("prec_prop uses plan defaults", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- prec_prop(p = 0.3, n = 400, plan = plan)
  ref <- prec_prop(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$se, ref$se)
  expect_equal(res$moe, ref$moe)
})

test_that("prec_mean uses plan defaults", {
  plan <- svyplan(deff = 1.5, N = 5000)
  res <- prec_mean(var = 100, n = 400, mu = 50, plan = plan)
  ref <- prec_mean(var = 100, n = 400, mu = 50, deff = 1.5, N = 5000)
  expect_equal(res$se, ref$se)
  expect_equal(res$cv, ref$cv)
})

test_that("power_prop uses plan defaults", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- power_prop(p1 = 0.30, p2 = 0.35, plan = plan)
  ref <- power_prop(p1 = 0.30, p2 = 0.35, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$n, ref$n)
})

test_that("power_prop explicit alpha wins over plan", {
  plan <- svyplan(alpha = 0.10)
  res <- power_prop(p1 = 0.30, p2 = 0.35, plan = plan, alpha = 0.05)
  ref <- power_prop(p1 = 0.30, p2 = 0.35, alpha = 0.05)
  expect_equal(res$n, ref$n)
})

test_that("power_mean uses plan defaults", {
  plan <- svyplan(deff = 1.5, alpha = 0.10)
  res <- power_mean(effect = 5, var = 100, plan = plan)
  ref <- power_mean(effect = 5, var = 100, deff = 1.5, alpha = 0.10)
  expect_equal(res$n, ref$n)
})

test_that("power_did uses plan defaults", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- power_did(
    treat = c(50, 55), control = c(50, 52),
    outcome = "mean", var = 100, effect = 3, plan = plan
  )
  ref <- power_did(
    treat = c(50, 55), control = c(50, 52),
    outcome = "mean", var = 100, effect = 3,
    deff = 1.5, resp_rate = 0.9
  )
  expect_equal(res$n, ref$n)
})

test_that("n_cluster uses plan for stage_cost/delta", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.85)
  res <- n_cluster(cv = 0.05, plan = plan)
  ref <- n_cluster(stage_cost = c(500, 50), delta = 0.05, cv = 0.05,
                   resp_rate = 0.85)
  expect_equal(res$n, ref$n)
  expect_equal(res$cv, ref$cv)
})

test_that("n_cluster with budget from plan", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05)
  res <- n_cluster(budget = 100000, plan = plan)
  ref <- n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_equal(res$n, ref$n)
})

test_that("n_cluster explicit stage_cost overrides plan", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05)
  res <- n_cluster(stage_cost = c(800, 80), cv = 0.05, plan = plan)
  ref <- n_cluster(stage_cost = c(800, 80), delta = 0.05, cv = 0.05)
  expect_equal(res$n, ref$n)
})

test_that("n_cluster requires stage_cost and delta", {
  expect_error(n_cluster(cv = 0.05), "'stage_cost' is required")
  expect_error(n_cluster(stage_cost = c(500, 50), cv = 0.05),
               "'delta' is required")
})

test_that("prec_cluster uses plan defaults", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.9)
  res <- prec_cluster(n = c(50, 12), plan = plan)
  ref <- prec_cluster(n = c(50, 12), delta = 0.05, resp_rate = 0.9)
  expect_equal(res$cv, ref$cv)
})

test_that("prec_cluster requires delta", {
  expect_error(prec_cluster(n = c(50, 12)), "'delta' is required")
})

test_that("n_multi uses plan for stage_cost", {
  targets <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  plan <- svyplan(stage_cost = c(500, 50))
  res <- suppressMessages(n_multi(targets, plan = plan))
  ref <- suppressMessages(n_multi(targets, stage_cost = c(500, 50)))
  expect_equal(res$n, ref$n)
})

test_that("n_alloc uses plan defaults", {
  frame <- data.frame(
    N_h = c(4000, 3000, 3000),
    S_h = c(10, 15, 8),
    mean_h = c(50, 60, 55)
  )
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- n_alloc(frame, n = 600, plan = plan)
  ref <- n_alloc(frame, n = 600, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$n, ref$n)
})

test_that("prec_alloc uses plan defaults", {
  frame <- data.frame(
    N_h = c(4000, 3000, 3000),
    S_h = c(10, 15, 8),
    mean_h = c(50, 60, 55)
  )
  plan <- svyplan(deff = 1.5, alpha = 0.10)
  res <- prec_alloc(frame, n = c(200, 200, 200), plan = plan)
  ref <- prec_alloc(frame, n = c(200, 200, 200), deff = 1.5, alpha = 0.10)
  expect_equal(res$se, ref$se)
  expect_equal(res$cv, ref$cv)
})

test_that("plan = NULL gives identical results to no plan", {
  res1 <- n_prop(p = 0.3, moe = 0.05, plan = NULL)
  res2 <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res1$n, res2$n)
})

test_that("empty plan is no-op", {
  plan <- svyplan()
  res <- n_prop(p = 0.3, moe = 0.05, plan = plan)
  ref <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$n, ref$n)
})

test_that("invalid plan argument errors", {
  expect_error(n_prop(p = 0.3, moe = 0.05, plan = list(deff = 1.5)),
               "'plan' must be a svyplan object")
})

test_that("plan with explicit default value (deff=1) -- explicit wins", {
  plan <- svyplan(deff = 1.8)
  res <- n_prop(p = 0.3, moe = 0.05, plan = plan, deff = 1)
  ref <- n_prop(p = 0.3, moe = 0.05, deff = 1)
  expect_equal(res$n, ref$n)
  expect_equal(res$params$deff, 1)
})

test_that("pipe: plan |> n_prop", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- plan |> n_prop(0.3, moe = 0.05)
  ref <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> n_mean", {
  plan <- svyplan(deff = 1.8, N = 50000)
  res <- plan |> n_mean(100, moe = 2)
  ref <- n_mean(var = 100, moe = 2, deff = 1.8, N = 50000)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> prec_prop", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- plan |> prec_prop(0.3, n = 400)
  ref <- prec_prop(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$se, ref$se)
})

test_that("pipe: plan |> prec_mean", {
  plan <- svyplan(deff = 1.5)
  res <- plan |> prec_mean(100, n = 400, mu = 50)
  ref <- prec_mean(var = 100, n = 400, mu = 50, deff = 1.5)
  expect_equal(res$cv, ref$cv)
})

test_that("pipe: plan |> power_prop", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- plan |> power_prop(0.30, p2 = 0.35)
  ref <- power_prop(p1 = 0.30, p2 = 0.35, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> power_mean", {
  plan <- svyplan(deff = 1.5)
  res <- plan |> power_mean(5, var = 100)
  ref <- power_mean(effect = 5, var = 100, deff = 1.5)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> power_did", {
  plan <- svyplan(deff = 1.5)
  res <- plan |> power_did(c(50, 55), control = c(50, 52),
                           outcome = "mean", var = 100, effect = 3)
  ref <- power_did(treat = c(50, 55), control = c(50, 52),
                   outcome = "mean", var = 100, effect = 3, deff = 1.5)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> n_cluster (stage_cost/delta from plan)", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.85)
  res <- plan |> n_cluster(cv = 0.05)
  ref <- n_cluster(stage_cost = c(500, 50), delta = 0.05, cv = 0.05,
                   resp_rate = 0.85)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> prec_cluster", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05)
  res <- plan |> prec_cluster(c(50, 12))
  ref <- prec_cluster(n = c(50, 12), delta = 0.05)
  expect_equal(res$cv, ref$cv)
})

test_that("pipe: plan |> n_multi", {
  targets <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  plan <- svyplan(stage_cost = c(500, 50))
  res <- suppressMessages(plan |> n_multi(targets))
  ref <- suppressMessages(n_multi(targets, stage_cost = c(500, 50)))
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> n_alloc", {
  frame <- data.frame(
    N_h = c(4000, 3000, 3000),
    S_h = c(10, 15, 8),
    mean_h = c(50, 60, 55)
  )
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res <- plan |> n_alloc(frame, n = 600)
  ref <- n_alloc(frame, n = 600, deff = 1.5, resp_rate = 0.9)
  expect_equal(res$n, ref$n)
})

test_that("pipe: plan |> prec_alloc", {
  frame <- data.frame(
    N_h = c(4000, 3000, 3000),
    S_h = c(10, 15, 8),
    mean_h = c(50, 60, 55)
  )
  plan <- svyplan(deff = 1.5, alpha = 0.10)
  res <- plan |> prec_alloc(frame, n = c(200, 200, 200))
  ref <- prec_alloc(frame, n = c(200, 200, 200), deff = 1.5, alpha = 0.10)
  expect_equal(res$se, ref$se)
})

test_that("pipe result equals plan= result for n_prop", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9, N = 10000)
  res_pipe <- plan |> n_prop(0.3, moe = 0.05)
  res_plan <- n_prop(p = 0.3, moe = 0.05, plan = plan)
  expect_equal(res_pipe$n, res_plan$n)
})

test_that("pipe result equals plan= result for n_cluster", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.85)
  res_pipe <- plan |> n_cluster(cv = 0.05)
  res_plan <- n_cluster(cv = 0.05, plan = plan)
  expect_equal(res_pipe$n, res_plan$n)
})

test_that("named pipe: plan |> n_prop(p = ...) matches all styles", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res_named_plan <- n_prop(p = 0.3, moe = 0.05, plan = plan)
  res_pos_pipe <- plan |> n_prop(0.3, moe = 0.05)
  res_named_pipe <- plan |> n_prop(p = 0.3, moe = 0.05)
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> n_mean(var = ...) matches all styles", {
  plan <- svyplan(deff = 1.8, N = 50000)
  res_named_plan <- n_mean(var = 100, moe = 2, plan = plan)
  res_pos_pipe <- plan |> n_mean(100, moe = 2)
  res_named_pipe <- plan |> n_mean(var = 100, moe = 2)
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> prec_prop(p = ...) matches all styles", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res_named_plan <- prec_prop(p = 0.3, n = 400, plan = plan)
  res_pos_pipe <- plan |> prec_prop(0.3, n = 400)
  res_named_pipe <- plan |> prec_prop(p = 0.3, n = 400)
  expect_equal(res_named_pipe$se, res_named_plan$se)
  expect_equal(res_named_pipe$se, res_pos_pipe$se)
})

test_that("named pipe: plan |> prec_mean(var = ...) matches all styles", {
  plan <- svyplan(deff = 1.5)
  res_named_plan <- prec_mean(var = 100, n = 400, mu = 50, plan = plan)
  res_pos_pipe <- plan |> prec_mean(100, n = 400, mu = 50)
  res_named_pipe <- plan |> prec_mean(var = 100, n = 400, mu = 50)
  expect_equal(res_named_pipe$cv, res_named_plan$cv)
  expect_equal(res_named_pipe$cv, res_pos_pipe$cv)
})

test_that("named pipe: plan |> power_prop(p1 = ...) matches all styles", {
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res_named_plan <- power_prop(p1 = 0.30, p2 = 0.35, plan = plan)
  res_pos_pipe <- plan |> power_prop(0.30, p2 = 0.35)
  res_named_pipe <- plan |> power_prop(p1 = 0.30, p2 = 0.35)
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> power_mean(effect = ...) matches all styles", {
  plan <- svyplan(deff = 1.5)
  res_named_plan <- power_mean(effect = 5, var = 100, plan = plan)
  res_pos_pipe <- plan |> power_mean(5, var = 100)
  res_named_pipe <- plan |> power_mean(effect = 5, var = 100)
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> power_did(treat = ...) matches all styles", {
  plan <- svyplan(deff = 1.5)
  args <- list(control = c(50, 52), outcome = "mean", var = 100, effect = 3)
  res_named_plan <- do.call(power_did, c(list(treat = c(50, 55), plan = plan), args))
  res_pos_pipe <- do.call(power_did, c(list(plan, c(50, 55)), args))
  res_named_pipe <- plan |> power_did(treat = c(50, 55), control = c(50, 52),
                                       outcome = "mean", var = 100, effect = 3)
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> n_cluster(cv = ...) matches all styles", {
  plan <- svyplan(stage_cost = c(500, 50), delta = 0.05, resp_rate = 0.85)
  res_named_plan <- n_cluster(cv = 0.05, plan = plan)
  res_pos_pipe <- plan |> n_cluster(cv = 0.05)
  expect_equal(res_pos_pipe$n, res_named_plan$n)
})

test_that("named pipe: plan |> prec_cluster(n = ...) matches all styles", {
  plan <- svyplan(delta = 0.05)
  res_named_plan <- prec_cluster(n = c(50, 12), plan = plan)
  res_pos_pipe <- plan |> prec_cluster(c(50, 12))
  res_named_pipe <- plan |> prec_cluster(n = c(50, 12))
  expect_equal(res_named_pipe$cv, res_named_plan$cv)
  expect_equal(res_named_pipe$cv, res_pos_pipe$cv)
})

test_that("named pipe: plan |> n_multi(targets = ...) matches all styles", {
  targets <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  plan <- svyplan(stage_cost = c(500, 50))
  res_named_plan <- suppressMessages(n_multi(targets, plan = plan))
  res_pos_pipe <- suppressMessages(plan |> n_multi(targets))
  res_named_pipe <- suppressMessages(plan |> n_multi(targets = targets))
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> n_alloc(frame = ...) matches all styles", {
  frame <- data.frame(
    N_h = c(4000, 3000, 3000),
    S_h = c(10, 15, 8),
    mean_h = c(50, 60, 55)
  )
  plan <- svyplan(deff = 1.5, resp_rate = 0.9)
  res_named_plan <- n_alloc(frame, n = 600, plan = plan)
  res_pos_pipe <- plan |> n_alloc(frame, n = 600)
  res_named_pipe <- plan |> n_alloc(frame = frame, n = 600)
  expect_equal(res_named_pipe$n, res_named_plan$n)
  expect_equal(res_named_pipe$n, res_pos_pipe$n)
})

test_that("named pipe: plan |> prec_alloc(x = ...) matches all styles", {
  frame <- data.frame(
    N_h = c(4000, 3000, 3000),
    S_h = c(10, 15, 8),
    mean_h = c(50, 60, 55)
  )
  plan <- svyplan(deff = 1.5)
  res_named_plan <- prec_alloc(frame, n = c(200, 200, 200), plan = plan)
  res_pos_pipe <- plan |> prec_alloc(frame, n = c(200, 200, 200))
  res_named_pipe <- plan |> prec_alloc(x = frame, n = c(200, 200, 200))
  expect_equal(res_named_pipe$se, res_named_plan$se)
  expect_equal(res_named_pipe$se, res_pos_pipe$se)
})

test_that("named pipe: plan |> prec_multi(targets = ...) with stage_cost", {
  targets <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    n = c(50, 50),
    psu_size = c(12, 12),
    delta_psu = c(0.02, 0.05)
  )
  plan <- svyplan(stage_cost = c(500, 50))
  res_explicit <- prec_multi(targets, stage_cost = c(500, 50))
  res_named_plan <- prec_multi(targets, plan = plan)
  res_named_pipe <- plan |> prec_multi(targets = targets)
  expect_equal(res_named_plan$cv, res_explicit$cv)
  expect_equal(res_named_pipe$cv, res_explicit$cv)
})

test_that("dispatch does not intercept svyplan_prec objects", {
  res_n <- n_prop(p = 0.3, moe = 0.05)
  res_prec <- prec_prop(res_n)
  expect_s3_class(res_prec, "svyplan_prec")
})

test_that("dispatch does not intercept svyplan_n objects", {
  res_prec <- prec_prop(p = 0.3, n = 400)
  res_n <- n_prop(res_prec)
  expect_s3_class(res_n, "svyplan_n")
})

test_that("multiple unnamed svyplan objects in ... errors", {
  plan1 <- svyplan(deff = 1.5)
  plan2 <- svyplan(deff = 2.0)
  expect_error(
    n_prop(0.3, plan1, plan2, moe = 0.05),
    "multiple svyplan objects"
  )
})

test_that("unnamed plan + named plan= errors", {
  p1 <- svyplan(deff = 1.5)
  p2 <- svyplan(deff = 2.0)
  expect_error(
    n_prop(0.3, p1, moe = 0.05, plan = p2),
    "multiple svyplan objects"
  )
})

test_that("first-arg plan + named plan= errors", {
  p1 <- svyplan(deff = 1.5)
  p2 <- svyplan(deff = 2.0)
  expect_error(
    n_prop(p1, p = 0.3, moe = 0.05, plan = p2),
    "multiple svyplan objects"
  )
})
