test_that("predict.svyplan_n returns correct shape for proportions", {
  x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
  nd <- expand.grid(deff = c(1, 1.5, 2), resp_rate = c(0.8, 0.9))
  res <- predict(x, nd)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 6L)
  expect_true(all(c("deff", "resp_rate", "n", "se", "moe", "cv") %in% names(res)))
})

test_that("higher deff produces larger n for proportions", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(deff = c(1, 2, 3))
  res <- predict(x, nd)

  expect_true(all(diff(res$n) > 0))
})

test_that("single-row identity for n_prop", {
  x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.9)
  nd <- data.frame(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.9)
  res <- predict(x, nd)

  expect_equal(res$n, x$n, tolerance = 1e-8)
  expect_equal(res$se, x$se, tolerance = 1e-8)
})

test_that("predict.svyplan_n switches from moe to cv mode", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(cv = c(0.05, 0.10))
  res <- predict(x, nd)

  ref1 <- n_prop(p = 0.3, cv = 0.05)
  ref2 <- n_prop(p = 0.3, cv = 0.10)
  expect_equal(res$n[1], ref1$n, tolerance = 1e-8)
  expect_equal(res$n[2], ref2$n, tolerance = 1e-8)
})

test_that("predict.svyplan_n works for means", {
  x <- n_mean(var = 100, moe = 2, deff = 1.5)
  nd <- data.frame(deff = c(1, 2))
  res <- predict(x, nd)

  expect_equal(nrow(res), 2L)
  expect_true(res$n[2] > res$n[1])
})

test_that("single-row identity for n_mean", {
  x <- n_mean(var = 100, mu = 50, cv = 0.05, N = 5000)
  nd <- data.frame(var = 100, mu = 50, cv = 0.05, N = 5000)
  res <- predict(x, nd)

  expect_equal(res$n, x$n, tolerance = 1e-8)
})

test_that("no duplicate columns when newdata overlaps result names", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(moe = c(0.03, 0.05))
  res <- predict(x, nd)
  expect_equal(sum(names(res) == "moe"), 1L)
  expect_equal(ncol(res), 5L)

  nd2 <- data.frame(cv = c(0.05, 0.10))
  res2 <- predict(x, nd2)
  expect_equal(sum(names(res2) == "cv"), 1L)

  cl <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  nd3 <- data.frame(cv = c(0.03, 0.05))
  res3 <- predict(cl, nd3)
  expect_equal(sum(names(res3) == "cv"), 1L)

  pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
  nd4 <- data.frame(n = c(200, 500))
  res4 <- predict(pw, nd4)
  expect_equal(sum(names(res4) == "n"), 1L)
})

test_that("predict.svyplan_n errors on more than one precision target in newdata", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(moe = 0.05, cv = 0.10)
  expect_error(predict(x, nd), "cannot contain more than one of")
})

test_that("predict.svyplan_n errors for multi-indicator results", {
  targets <- data.frame(
    name = c("a", "b"),
    p    = c(0.3, 0.5),
    moe  = c(0.05, 0.05)
  )
  x <- n_multi(targets)
  expect_error(predict(x, data.frame(deff = 1)), "multi-indicator")
})

test_that("predict.svyplan_cluster varies budget", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  nd <- data.frame(budget = c(50000, 100000, 200000))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(all(c("n_psu", "n_per_psu", "total_n", "cv", "cost") %in% names(res)))
  expect_true(all(diff(res$cv) < 0))
})

test_that("predict.svyplan_cluster varies cv", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  nd <- data.frame(cv = c(0.03, 0.05, 0.10))
  res <- predict(x, nd)

  expect_true(all(diff(res$cost) < 0))
})

test_that("predict.svyplan_cluster supports varying icc only", {
  cl <- n_cluster(stage_cost = c(750, 100), icc = 0.05, unit_relvar = 1, var_ratio = 1,
                  budget = 1e5)
  nd <- data.frame(icc = c(0.01, 0.05, 0.10, 0.20))
  res <- predict(cl, nd)

  expect_equal(nrow(res), 4L)
  expect_equal(res$icc, nd$icc)

  ref_cv <- vapply(nd$icc, function(d) {
    n_cluster(
      stage_cost = c(750, 100), icc = d, unit_relvar = 1, var_ratio = 1, budget = 1e5
    )$cv
  }, numeric(1))
  expect_equal(res$cv, ref_cv, tolerance = 1e-8)
})

test_that("predict.svyplan_cluster supports varying var_ratio in 2-stage", {
  cl <- n_cluster(stage_cost = c(750, 100), icc = 0.05, unit_relvar = 1, var_ratio = 1,
                  budget = 1e5)
  nd <- data.frame(var_ratio_psu = c(0.8, 1.0, 1.2))
  res <- predict(cl, nd)

  ref_cv <- vapply(nd$var_ratio_psu, function(kval) {
    n_cluster(
      stage_cost = c(750, 100), icc = 0.05, unit_relvar = 1, var_ratio = kval, budget = 1e5
    )$cv
  }, numeric(1))
  expect_equal(res$cv, ref_cv, tolerance = 1e-8)

  nd_alias <- data.frame(var_ratio = c(0.8, 1.0, 1.2))
  res_alias <- predict(cl, nd_alias)
  expect_equal(res_alias$cv, ref_cv, tolerance = 1e-8)
})

test_that("predict.svyplan_cluster varies unit_relvar", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, unit_relvar = 1, cv = 0.05)
  nd <- data.frame(unit_relvar = c(0.5, 1, 2))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(all(diff(res$total_n) > 0))
})

test_that("predict.svyplan_cluster uses original mode when no cv/budget in newdata", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  nd <- data.frame(resp_rate_psu = c(0.8, 0.9, 1.0))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(res$n_psu[1] > res$n_psu[3])
})

test_that("predict.svyplan_cluster supports varying stage costs", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  nd <- data.frame(cost_psu = c(500, 800), cost_ssu = c(50, 120))
  res <- predict(x, nd)
  ref <- vapply(seq_len(nrow(nd)), function(i) {
    n_cluster(
      stage_cost = c(nd$cost_psu[i], nd$cost_ssu[i]),
      icc = 0.05,
      budget = 100000
    )$cv
  }, numeric(1))
  expect_equal(res$cv, ref, tolerance = 1e-8)

  nd_alias <- data.frame(cost_psu = c(500, 800), cost_tsu = c(50, 120))
  res_alias <- predict(x, nd_alias)
  expect_equal(res_alias$cv, ref, tolerance = 1e-8)
})

test_that("predict.svyplan_cluster supports stage-wise icc and var_ratio in 3-stage", {
  x <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    var_ratio = c(1.0, 1.0),
    budget = 200000
  )
  nd <- data.frame(
    icc_psu = c(0.01, 0.02),
    icc_ssu = c(0.05, 0.08),
    var_ratio_psu = c(1.0, 1.2),
    var_ratio_ssu = c(1.0 * (1 - 0.01), 1.2 * (1 - 0.02))
  )
  res <- predict(x, nd)

  ref <- vapply(seq_len(nrow(nd)), function(i) {
    n_cluster(
      stage_cost = c(500, 100, 50),
      icc = c(nd$icc_psu[i], nd$icc_ssu[i]),
      var_ratio = c(nd$var_ratio_psu[i], nd$var_ratio_ssu[i]),
      budget = 200000
    )$cv
  }, numeric(1))
  expect_equal(res$cv, ref, tolerance = 1e-8)
})

test_that("predict.svyplan_cluster rejects overlapping cost aliases", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  nd <- data.frame(cost_ssu = 50, cost_tsu = 50)
  expect_error(predict(x, nd), "multiple columns")
})

test_that("predict.svyplan_cluster rejects overlapping icc/var_ratio aliases", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = 1, budget = 100000)
  expect_error(
    predict(x, data.frame(icc = 0.05, icc_psu = 0.05)),
    "multiple columns"
  )
  expect_error(
    predict(x, data.frame(var_ratio = 1, var_ratio_psu = 1)),
    "multiple columns"
  )
})

test_that("predict.svyplan_cluster rejects scalar icc/var_ratio for 3-stage", {
  x <- n_cluster(
    stage_cost = c(500, 100, 50),
    icc = c(0.01, 0.05),
    var_ratio = c(1, 1),
    budget = 200000
  )
  expect_error(
    predict(x, data.frame(icc = 0.02)),
    "icc_psu'.*icc_ssu"
  )
  expect_error(
    predict(x, data.frame(var_ratio = 1.1)),
    "var_ratio_psu'.*var_ratio_ssu"
  )
})

test_that("predict.svyplan_cluster errors on both cv and budget", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  nd <- data.frame(cv = 0.05, budget = 100000)
  expect_error(predict(x, nd), "cannot contain both")
})

test_that("predict.svyplan_cluster rejects a several-indicators result", {
  targets <- data.frame(
    name   = c("a", "b"),
    p      = c(0.3, 0.1),
    cv     = c(0.10, 0.15),
    icc_psu = c(0.02, 0.05)
  )
  x <- n_cluster(indicators = targets, stage_cost = c(500, 50))
  expect_error(predict(x, data.frame(cv = 0.05)), "several-indicators")
})

test_that("predict.svyplan_power varies n (solved for power)", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
  nd <- data.frame(n = seq(200, 800, 200))
  res <- predict(pw, nd)

  expect_equal(nrow(res), 4L)
  expect_true(all(c("n", "power", "effect") %in% names(res)))
  expect_true(all(diff(res$power) > 0))
})

test_that("predict.svyplan_power varies power (solved for n)", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35)
  nd <- data.frame(power = c(0.70, 0.80, 0.90))
  res <- predict(pw, nd)

  expect_true(all(diff(res$n) > 0))
})

test_that("predict.svyplan_power works for means", {
  pw <- power_mean(effect = 5, var = 100, n = 200, power = NULL)
  nd <- data.frame(n = c(100, 200, 400))
  res <- predict(pw, nd)

  expect_true(all(diff(res$power) > 0))
})

test_that("predict.svyplan_power errors for solved-for param in newdata", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35)
  expect_error(predict(pw, data.frame(n = 500)), "unknown parameter")
})

test_that("predict.svyplan_prec varies n for proportions", {
  x <- prec_prop(p = 0.3, n = 400)
  nd <- data.frame(n = c(100, 400, 1600))
  res <- predict(x, nd)

  expect_equal(nrow(res), 3L)
  expect_true(all(c("se", "moe", "cv") %in% names(res)))
  expect_true(all(diff(res$se) < 0))
})

test_that("predict.svyplan_prec varies n for means", {
  x <- prec_mean(var = 100, n = 400, mu = 50)
  nd <- data.frame(n = c(100, 400, 1600))
  res <- predict(x, nd)

  expect_true(all(diff(res$se) < 0))
})

test_that("single-row identity for prec_prop", {
  x <- prec_prop(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.9)
  nd <- data.frame(p = 0.3, n = 400, deff = 1.5, resp_rate = 0.9)
  res <- predict(x, nd)

  expect_equal(res$se, x$se, tolerance = 1e-8)
  expect_equal(res$moe, x$moe, tolerance = 1e-8)
})

test_that("predict.svyplan_prec errors for cluster type", {
  x <- prec_cluster(n = c(50, 12), icc = 0.05)
  expect_error(predict(x, data.frame(n = 100)), "not supported")
})

test_that("predict errors for non-data.frame newdata", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(predict(x, list(deff = 1)), "must be a data frame")
})

test_that("predict errors for empty newdata", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(predict(x, data.frame()), "at least one")
})

test_that("predict errors for unknown param names", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(predict(x, data.frame(bogus = 1)), "unknown parameter")
})

test_that("predict errors for non-numeric columns", {
  x <- n_prop(p = 0.3, moe = 0.05)
  expect_error(
    predict(x, data.frame(deff = "high", stringsAsFactors = FALSE)),
    "non-numeric"
  )
})

test_that("predict produces NA + warning on row failure", {
  x <- n_prop(p = 0.3, moe = 0.05)
  nd <- data.frame(p = c(0.3, 0, 0.5))
  expect_warning(
    res <- predict(x, nd),
    "predict row 2 failed"
  )
  expect_true(is.na(res$n[2]))
  expect_false(is.na(res$n[1]))
  expect_false(is.na(res$n[3]))
})

test_that("predict.svyplan_cluster supports fixed_cost in newdata", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                  fixed_cost = 5000)
  nd <- data.frame(fixed_cost = c(0, 5000, 10000))
  res <- predict(x, nd)
  expect_equal(nrow(res), 3L)
  expect_true("fixed_cost" %in% names(res))
  expect_true(res$cost[3] > res$cost[1])
  expect_equal(res$n_psu[1], res$n_psu[2], tolerance = 1e-8)
})

test_that("predict.svyplan_power errors for vector n object", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  expect_error(
    predict(pw, data.frame(power = 0.90)),
    "does not support power objects with unequal-group n"
  )
})

test_that("predict.svyplan_power passes method through for power_prop", {
  pw <- power_prop(p1 = 0.15, p2 = 0.18, n = 2000, power = NULL,
                    alternative = "one.sided", method = "arcsine")
  nd <- data.frame(n = c(1000, 2000, 3000))
  res <- predict(pw, nd)
  expect_equal(nrow(res), 3L)
  expect_true(all(diff(res$power) > 0))
})

test_that("predict.svyplan_power supports varying alternative", {
  pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
  nd <- data.frame(
    alternative = c("two.sided", "one.sided"),
    stringsAsFactors = FALSE
  )
  res <- predict(pw, nd)
  expect_equal(nrow(res), 2L)
  expect_true(res$power[2] > res$power[1])
})

test_that("predict works for cluster-mode n_alloc objects", {
  fr <- data.frame(stratum = c("A", "B"), N = c(50000, 150000),
                   sd = c(0.45, 0.48), mean = c(0.35, 0.25),
                   icc_psu = c(0.03, 0.08), n_per_psu = c(12, 12))
  x <- n_alloc(fr, n = 1000)
  p <- predict(x, data.frame(n = c(1000, 1100)))
  expect_equal(nrow(p), 2L)
  expect_equal(p$cv[1], x$cv, tolerance = 1e-10)
})

test_that("predict keeps the beta interval's df", {
  x <- n_prop(0.02, moe = 0.01, method = "beta", df = 25)
  # deff = 1 is what the object already has, so the size must not move
  expect_equal(predict(x, data.frame(deff = 1))$n, x$n)

  p <- prec_prop(p = 0.02, n = 900, deff = 2, method = "beta", df = 25)
  expect_equal(predict(p, data.frame(deff = 2))$moe, p$moe)
})

test_that("df can be varied on a beta grid", {
  x <- n_prop(0.02, moe = 0.01, method = "beta", df = 25)
  g <- predict(x, data.frame(df = c(10, 25, 100)))

  expect_equal(g$n[2], x$n)
  # fewer degrees of freedom widen the interval, so the size rises
  expect_true(all(diff(g$n) < 0))
})

test_that("documented sizing parameters are accepted in compatible grids", {
  cases <- list(
    proportion = list(
      object = n_prop(p = 0.3, moe = 0.05),
      grids = list(
        p = data.frame(p = 0.35), moe = data.frame(moe = 0.04),
        rmoe = data.frame(rmoe = 0.15), cv = data.frame(cv = 0.1),
        alpha = data.frame(alpha = 0.1), N = data.frame(N = 5000),
        deff = data.frame(deff = 1.2), resp_rate = data.frame(resp_rate = 0.9),
        df = data.frame(df = 30), min_cases = data.frame(min_cases = 20)
      )
    ),
    mean = list(
      object = n_mean(var = 100, mu = 50, moe = 2),
      grids = list(
        var = data.frame(var = 120), mu = data.frame(mu = 45),
        moe = data.frame(moe = 1.8), rmoe = data.frame(rmoe = 0.04),
        cv = data.frame(cv = 0.04), alpha = data.frame(alpha = 0.1),
        N = data.frame(N = 5000), deff = data.frame(deff = 1.2),
        resp_rate = data.frame(resp_rate = 0.9), df = data.frame(df = 30)
      )
    ),
    ratio = list(
      object = n_ratio(r = 2, cv_num = 0.5, cv_den = 0.3,
                       component_cor = 0.4, cv = 0.1),
      grids = list(
        r = data.frame(r = 2.2), cv_num = data.frame(cv_num = 0.55),
        cv_den = data.frame(cv_den = 0.32),
        component_cor = data.frame(component_cor = 0.3),
        moe = data.frame(moe = 0.2), rmoe = data.frame(rmoe = 0.1),
        cv = data.frame(cv = 0.08), alpha = data.frame(alpha = 0.1),
        N = data.frame(N = 5000), deff = data.frame(deff = 1.2),
        resp_rate = data.frame(resp_rate = 0.9), df = data.frame(df = 30)
      )
    )
  )

  for (case in cases) {
    for (grid in case$grids) expect_s3_class(predict(case$object, grid), "data.frame")
    expect_error(predict(case$object, data.frame(undocumented = 1)), "unknown parameter")
  }
})

test_that("documented precision parameters are accepted in compatible grids", {
  cases <- list(
    proportion = list(
      object = prec_prop(p = 0.3, n = 400),
      grids = list(
        p = data.frame(p = 0.35), n = data.frame(n = 500),
        alpha = data.frame(alpha = 0.1), N = data.frame(N = 5000),
        deff = data.frame(deff = 1.2), resp_rate = data.frame(resp_rate = 0.9),
        df = data.frame(df = 30)
      )
    ),
    mean = list(
      object = prec_mean(var = 100, n = 400, mu = 50),
      grids = list(
        var = data.frame(var = 120), n = data.frame(n = 500),
        mu = data.frame(mu = 45), alpha = data.frame(alpha = 0.1),
        N = data.frame(N = 5000), deff = data.frame(deff = 1.2),
        resp_rate = data.frame(resp_rate = 0.9), df = data.frame(df = 30)
      )
    ),
    ratio = list(
      object = prec_ratio(r = 2, n = 400, cv_num = 0.5, cv_den = 0.3,
                          component_cor = 0.4),
      grids = list(
        r = data.frame(r = 2.2), n = data.frame(n = 500),
        cv_num = data.frame(cv_num = 0.55), cv_den = data.frame(cv_den = 0.32),
        component_cor = data.frame(component_cor = 0.3),
        alpha = data.frame(alpha = 0.1), N = data.frame(N = 5000),
        deff = data.frame(deff = 1.2), resp_rate = data.frame(resp_rate = 0.9),
        df = data.frame(df = 30)
      )
    )
  )

  for (case in cases) {
    for (grid in case$grids) expect_s3_class(predict(case$object, grid), "data.frame")
    expect_error(predict(case$object, data.frame(undocumented = 1)), "unknown parameter")
  }

  solved_p <- prec_prop(p = NULL, n = 400, cv = 0.1)
  expect_s3_class(predict(solved_p, data.frame(cv = 0.08)), "data.frame")
  solved_rmoe <- prec_prop(p = NULL, n = 400, rmoe = 0.3, method = "wilson")
  expect_s3_class(predict(solved_rmoe, data.frame(rmoe = 0.25)), "data.frame")
  solved_mu <- prec_mean(var = 100, n = 400, cv = 0.05)
  expect_s3_class(predict(solved_mu, data.frame(cv = 0.04)), "data.frame")
})
