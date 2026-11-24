## Eligible units per sampled unit: the proportion family

.wald_cv_n <- function(p, cv, N, deff) {
  q <- 1 - p
  a <- if (is.infinite(N)) 1 else N / (N - 1)
  a * deff * q / p / (cv^2 + deff * q / p / (N - 1))
}

test_that("at N = Inf a rate at most 1 equals folding it into resp_rate", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    targets <- list(list(moe = 0.03), list(rmoe = 0.15))
    if (method == "wald") targets <- c(targets, list(list(cv = 0.08)))
    for (target in targets) {
      elig <- do.call(n_prop, c(list(p = 0.28, deff = 2.5, resp_rate = 0.9,
                                     eligible_per_unit = 0.41, method = method),
                                target))
      folded <- do.call(n_prop, c(list(p = 0.28, deff = 2.5,
                                       resp_rate = 0.9 * 0.41,
                                       method = method), target))
      expect_identical(elig$n, folded$n)
    }
  }
})

test_that("a finite N puts the eligible population in the correction", {
  x <- n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9,
              eligible_per_unit = 0.41, N = 20000)
  expect_equal(x$n, .wald_cv_n(0.28, 0.05, 20000 * 0.41, 2.5) / (0.9 * 0.41))
  folded <- n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9 * 0.41,
                   N = 20000)
  expect_gt(folded$n, x$n)
})

test_that("a rate above 1 is accepted and solved on the same scale", {
  x <- n_prop(p = 0.2, cv = 0.08, eligible_per_unit = 2.1, N = 5000)
  expect_equal(x$n, .wald_cv_n(0.2, 0.08, 5000 * 2.1, 1) / 2.1)
  expect_lt(x$n, n_prop(p = 0.2, cv = 0.08, N = 5000)$n)
})

test_that("cases count eligible respondents", {
  x <- n_prop(p = 0.02, moe = 0.02, resp_rate = 0.9, eligible_per_unit = 0.6,
              min_cases = 40)
  expect_identical(x$binding, "min_cases")
  expect_equal(x$n, 40 / (0.02 * 0.9 * 0.6))
  expect_equal(x$expected_cases, 40)
})

test_that("prec_prop forms the same precision as n_prop targets", {
  ev <- prec_prop(p = 0.28, n = 6000, deff = 2.5, resp_rate = 0.9,
                  eligible_per_unit = 0.41, N = 20000)
  n_net <- 6000 * 0.9 * 0.41
  N_elig <- 20000 * 0.41
  n_eff <- n_net * (N_elig - 1) / (2.5 * (N_elig - n_net))
  expect_equal(ev$se, sqrt(0.28 * 0.72 / n_eff))

  floor_p <- prec_prop(n = 6000, cv = 0.1, resp_rate = 0.9,
                       eligible_per_unit = 0.41, N = 20000)
  expect_equal(prec_prop(p = floor_p$params$p, n = 6000, resp_rate = 0.9,
                         eligible_per_unit = 0.41, N = 20000)$cv, 0.1)
})

test_that("the rate survives every round trip", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    x <- n_prop(p = 0.28, moe = 0.03, deff = 2.5, resp_rate = 0.9,
                eligible_per_unit = 0.41, N = 20000, method = method)
    expect_equal(x$moe, 0.03)
    ev <- prec_prop(x)
    expect_equal(ev$moe, 0.03)
    expect_identical(ev$params$eligible_per_unit, 0.41)
    expect_equal(n_prop(ev)$n, x$n)
  }
  x <- n_prop(p = 0.28, cv = 0.05, eligible_per_unit = 1.6, N = 8000)
  swept <- predict(x, data.frame(eligible_per_unit = c(0.4, 1.6, 2.5)))
  expect_equal(swept$cv, rep(0.05, 3))
  expect_equal(swept$n[2], x$n)
  ev <- prec_prop(x)
  expect_equal(predict(ev, data.frame(eligible_per_unit = 1.6))$cv, 0.05)
  expect_equal(predict(x, data.frame(cv = c(0.05, 0.04)))$n[1], x$n)
  expect_equal(predict(ev, data.frame(n = x$n))$cv, 0.05)
})

test_that("a profile carries the rate and refuses a function without it", {
  pl <- svyplan(eligible_per_unit = 0.41, resp_rate = 0.9)
  expect_identical(
    n_prop(p = 0.28, cv = 0.05, plan = pl)$n,
    n_prop(p = 0.28, cv = 0.05, resp_rate = 0.9, eligible_per_unit = 0.41)$n
  )
  expect_identical(
    prec_prop(p = 0.28, n = 6000, plan = pl)$params$eligible_per_unit, 0.41
  )
  expect_identical(
    n_mean(var = 4, mu = 10, cv = 0.05, plan = pl)$n,
    n_mean(var = 4, mu = 10, cv = 0.05, resp_rate = 0.9,
           eligible_per_unit = 0.41)$n
  )
  expect_equal(
    n_multi(data.frame(name = "a", p = 0.3, moe = 0.05), plan = pl)$n,
    n_prop(p = 0.3, moe = 0.05, resp_rate = 0.9, eligible_per_unit = 0.41)$n
  )
  piped <- tryCatch(
    pl |> n_ratio(r = 0.5, var_y = 1, var_x = 1, cov_xy = 0.5, mu_x = 2,
                  cv = 0.05),
    error = conditionMessage
  )
  expect_identical(
    piped,
    paste("the profile sets 'eligible_per_unit', which n_ratio() does not take.",
          "Drop it from the profile for this call, or fold it into",
          "'resp_rate' when it is at most 1")
  )
  expect_identical((pl |> prec_mean(var = 4, n = 500))$params$eligible_per_unit,
                   0.41)
  expect_no_error(n_ratio(r = 0.5, cv_num = 0.3, cv_den = 0.2,
                          component_cor = 0.5, cv = 0.05,
                          plan = svyplan(eligible_per_unit = 1)))
  expect_error(svyplan(eligible_per_unit = 0), "positive finite")
})

test_that("the rate is validated and stays out of a default result", {
  for (bad in list(0, -1, Inf, NA_real_, c(0.4, 0.5), "0.4")) {
    expect_error(n_prop(p = 0.3, moe = 0.05, eligible_per_unit = bad),
                 "'eligible_per_unit' must be a positive finite number")
    expect_error(prec_prop(p = 0.3, n = 400, eligible_per_unit = bad),
                 "'eligible_per_unit' must be a positive finite number")
  }
  expect_null(n_prop(p = 0.3, moe = 0.05)$params$eligible_per_unit)
  expect_null(prec_prop(p = 0.3, n = 400)$params$eligible_per_unit)
  expect_error(n_prop(p = 0.3, cv = 0.01, resp_rate = 0.5,
                      eligible_per_unit = 0.2, N = 1000),
               "exceeds the population")
})

test_that("the print counts eligible respondents", {
  x <- n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9,
              eligible_per_unit = 0.41)
  out <- capture.output(print(x))
  expect_match(out[2], "eligible respondents: 2572")
  expect_match(out[2], "eligible_per_unit = 0.41")
  expect_match(capture.output(print(prec_prop(x)))[2],
               "eligible respondents: 2571.429, eligible_per_unit = 0.41")
  expect_false(any(grepl("eligible",
                         capture.output(print(n_prop(p = 0.3, moe = 0.05))))))
})

test_that("the interval is as wide as the reported margin", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    for (rate in c(0.25, 2.6)) {
      x <- n_prop(p = 0.3, moe = 0.04, eligible_per_unit = rate, N = 20000,
                  deff = 1.5, resp_rate = 0.8, method = method)
      prec <- prec_prop(p = 0.3, n = x$n, eligible_per_unit = rate, N = 20000,
                        deff = 1.5, resp_rate = 0.8, method = method)
      label <- paste(method, rate)
      expect_equal(diff(confint(x)[1, ]) / 2, 0.04, ignore_attr = TRUE,
                   tolerance = 1e-8, label = label)
      expect_equal(diff(confint(prec)[1, ]) / 2, prec$moe,
                   ignore_attr = TRUE, tolerance = 1e-8, label = label)
    }
  }
})

## Eligible units per sampled unit: the cluster family

test_that("cluster sizes equal folding a rate at most 1 into resp_rate", {
  calls <- list(
    list(stage_cost = c(500, 50), icc = 0.05, cv = 0.05),
    list(stage_cost = c(500, 50), icc = 0.05, budget = 1e5),
    list(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, stage_n = c(80, NA)),
    list(stage_cost = c(500, 50), icc = 0.05, budget = 1e5, stage_n = c(NA, 25)),
    list(stage_cost = c(500, 100, 20), icc = c(0.03, 0.1), cv = 0.05)
  )
  for (args in calls) {
    elig <- do.call(n_cluster, c(args, resp_rate = 0.9, eligible_per_unit = 0.41))
    folded <- do.call(n_cluster, c(args, resp_rate = 0.9 * 0.41))
    expect_equal(elig$n, folded$n)
    expect_equal(elig$operational$n, folded$operational$n)
    S <- length(args$stage_cost)
    expect_identical(
      elig$params$resp_rate,
      stats::setNames(c(rep(1, S - 1L), 0.9),
                      c(paste0("resp_rate_", seq_len(S - 1L)), "resp_rate"))
    )
  }
})

test_that("the household take is the eligible-person optimum converted", {
  for (e in c(0.41, 2.5)) {
    x <- n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05,
                   resp_rate = 0.9, eligible_per_unit = e)
    cost_per_eligible <- 50 / (e * 0.9)
    persons <- sqrt(500 / cost_per_eligible * (1 - 0.07) / 0.07)
    expect_equal(unname(x$n["stage_2"]), persons / (e * 0.9))
  }
})

test_that("prec_cluster reads the eligible respondents per cluster", {
  ev <- prec_cluster(n = c(600, 20), icc = 0.07, resp_rate = 0.9,
                     eligible_per_unit = 0.41)
  b <- 20 * 0.9 * 0.41
  expect_equal(ev$cv, sqrt((1 + 0.07 * (b - 1)) / (600 * b)))
  expect_identical(ev$params$eligible_per_unit, 0.41)
  expect_identical(ev$params$resp_rate, c(resp_rate_1 = 1, resp_rate = 0.9))
})

test_that("the cluster rate survives every round trip", {
  for (args in list(
    list(stage_cost = c(500, 50), icc = 0.07, cv = 0.05),
    list(stage_cost = c(500, 100, 20), icc = c(0.03, 0.1), budget = 1e5)
  )) {
    x <- do.call(n_cluster, c(args, resp_rate = 0.9, eligible_per_unit = 1.8))
    ev <- prec_cluster(x)
    expect_equal(ev$cv, x$cv)
    expect_identical(ev$params$eligible_per_unit, 1.8)
    expect_equal(n_cluster(ev)$n, x$n)
  }
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05,
                 eligible_per_unit = 0.41)
  expect_equal(n_cluster(prec_cluster(x))$n, x$n)
  swept <- predict(x, data.frame(eligible_per_unit = c(0.41, 1.2)))
  expect_equal(swept$cv, c(0.05, 0.05))
  expect_equal(swept$stage_1[1], unname(x$n["stage_1"]))
  by_cv <- predict(x, data.frame(cv = c(0.05, 0.04)))
  expect_equal(by_cv$stage_2[1], unname(x$n["stage_2"]))
})

test_that("an indicator table sizes each row with its own eligible_per_unit", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.5), moe = 0.05,
                    N = 20000, resp_rate = 0.9, eligible_per_unit = c(0.41, 2.3))
  fit <- n_multi(ind)
  expect_equal(
    fit$detail$.n,
    c(n_prop(0.3, moe = 0.05, N = 20000, resp_rate = 0.9,
             eligible_per_unit = 0.41)$n,
      n_prop(0.5, moe = 0.05, N = 20000, resp_rate = 0.9,
             eligible_per_unit = 2.3)$n)
  )
  expect_equal(
    fit$detail$.cv_achieved[2],
    prec_prop(0.5, n = fit$n, N = 20000, resp_rate = 0.9,
              eligible_per_unit = 2.3)$cv
  )
  expect_equal(n_multi(prec_multi(fit))$n, fit$n)

  scalar <- n_multi(transform(ind, eligible_per_unit = NULL), eligible_per_unit = 2.3)
  expect_equal(scalar$detail$.n[2], fit$detail$.n[2])
  expect_equal(
    n_multi(transform(ind, eligible_per_unit = c(NA, 2.3)), eligible_per_unit = 0.41)$n,
    fit$n
  )

  doms <- transform(rbind(ind, ind), dom = rep(c("A", "B"), each = 2))
  expect_equal(
    n_multi(doms, domains = "dom")$domains$.n,
    rep(fit$n, 2)
  )

  pr <- prec_multi(transform(ind, moe = NULL, n = 800))
  expect_equal(
    pr$cv[1],
    prec_prop(0.3, n = 800, N = 20000, resp_rate = 0.9,
              eligible_per_unit = 0.41)$cv
  )
  expect_null(n_multi(ind[, 1:3])$indicators$eligible_per_unit)
})

test_that("a mean row takes the rate, a ratio row refuses it, and the cluster mode reads it", {
  mixed <- data.frame(name = c("a", "m"), p = c(0.3, NA), var = c(NA, 4),
                      mean = c(NA, 10), cv = c(0.1, 0.02))
  fit <- n_multi(mixed, eligible_per_unit = 2)
  expect_equal(fit$detail$.n[2L],
               n_mean(var = 4, mu = 10, cv = 0.02, eligible_per_unit = 2)$n)
  expect_equal(fit$detail$.cv_achieved[2L],
               prec_mean(var = 4, mu = 10, n = max(fit$detail$.n),
                         eligible_per_unit = 2)$cv)
  expect_equal(fit$detail$.n[2L],
               n_mean(var = 4, mu = 10, cv = 0.02)$n / 2)
  expect_equal(
    prec_multi(transform(mixed, n = 300), eligible_per_unit = 2)$detail$.se[2L],
    prec_mean(var = 4, mu = 10, n = 300, eligible_per_unit = 2)$se
  )
  ratio <- data.frame(name = "r", r = 420, cv_num = 1.2, cv_den = 0.45,
                      component_cor = 0.65, cv = 0.05)
  expect_error(n_multi(ratio, eligible_per_unit = 2),
               "proportion and mean rows.*Row\\(s\\) 1 carry a ratio")
  expect_error(
    n_multi(data.frame(name = "a", p = 0.3, cv = 0.1, eligible_per_unit = 0)),
    "positive and finite"
  )
  ind <- data.frame(name = "a", p = 0.3, moe = 0.05, icc_1 = 0.05,
                    eligible_per_unit = 0.41)
  cl_col <- n_cluster(ind, stage_cost = c(500, 50))
  cl_arg <- n_cluster(transform(ind, eligible_per_unit = NULL),
                      stage_cost = c(500, 50), eligible_per_unit = 0.41)
  expect_equal(cl_col$n, cl_arg$n)
  expect_equal(prec_cluster(cl_col)$detail$.cv, cl_col$cv)
  expect_equal(
    prec_cluster(transform(ind, moe = NULL, n_1 = 50, take_2 = 10))$detail$.cv,
    prec_cluster(n = c(50, 10), icc = 0.05, unit_relvar = 0.7 / 0.3,
                 eligible_per_unit = 0.41)$cv
  )
})

test_that("cluster indicator mode supports eligible_per_unit and the take check names it", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.2), cv = c(0.1, 0.12),
                    icc_1 = c(0.05, 0.08))
  fit_arg <- n_cluster(ind, stage_cost = c(500, 50), eligible_per_unit = 0.4)
  fit_col <- n_cluster(transform(ind, eligible_per_unit = c(0.4, 0.4)),
                       stage_cost = c(500, 50))
  expect_equal(fit_arg$n, fit_col$n)
  # Per-row eligible_per_unit including > 1 and stage FPC round-trip
  ind_var <- transform(ind, eligible_per_unit = c(0.5, 2.2))
  fit_fpc <- n_cluster(ind_var, stage_cost = c(500, 50), fpc = "stage",
                       N = c(500, 50000))
  pr_fpc <- prec_cluster(fit_fpc)
  expect_identical(pr_fpc$type, "cluster")
  expect_equal(n_cluster(pr_fpc)$n, fit_fpc$n, tolerance = 1e-5)
  expect_equal(
    prec_cluster(transform(ind, n_1 = 50, take_2 = 10),
                 eligible_per_unit = 0.4)$detail$.cv[1L],
    prec_cluster(n = c(50, 10), icc = 0.05, unit_relvar = 0.7 / 0.3,
                 eligible_per_unit = 0.4)$cv
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
              stage_n = c(NA, 2), resp_rate = 0.9, eligible_per_unit = 0.4),
    "at 'resp_rate \\* eligible_per_unit' = 0.36"
  )
  expect_error(
    prec_cluster(n = c(50, 2), icc = 0.05, resp_rate = 0.9,
                 eligible_per_unit = 0.4),
    "at 'resp_rate \\* eligible_per_unit' = 0.36"
  )
})

test_that("a profile carries the rate to the cluster functions", {
  pl <- svyplan(eligible_per_unit = 0.41, resp_rate = 0.9)
  expect_equal(
    n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05, plan = pl)$n,
    n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05,
              resp_rate = 0.9, eligible_per_unit = 0.41)$n
  )
  expect_identical(
    prec_cluster(n = c(60, 20), icc = 0.07, plan = pl)$params$eligible_per_unit,
    0.41
  )
})

test_that("the cluster print counts eligible respondents", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                 resp_rate = 0.9, eligible_per_unit = 0.41)
  out <- capture.output(print(x))
  expect_match(out[2], "eligible respondents: 558")
  expect_match(out[3], "resp_rate = 0.90, eligible_per_unit = 0.41")
  expect_match(capture.output(print(prec_cluster(x)))[2],
               "eligible respondents: 547.4634, eligible_per_unit = 0.41")
  plain <- capture.output(print(n_cluster(stage_cost = c(500, 50),
                                          icc = 0.05, cv = 0.05)))
  expect_false(any(grepl("eligible", plain)))
  expect_null(n_cluster(stage_cost = c(500, 50), icc = 0.05,
                        cv = 0.05)$params$eligible_per_unit)
})

## Eligible units per sampled unit: joint allocation

.elig_joint <- function() {
  frame <- data.frame(
    stratum = c("A-u", "A-r", "B-u", "B-r"),
    region = rep(c("A", "B"), each = 2),
    N = c(4e5, 6e5, 3e5, 5e5), take_2 = 20, cost_1 = 400,
    cost_2 = 30, N_1 = c(2000, 3000, 1500, 2500)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("mcpr", "stunting"), each = 4),
    p = c(0.30, 0.24, 0.28, 0.22, 0.20, 0.30, 0.22, 0.33),
    icc_1 = rep(c(0.07, 0.05), each = 4), resp_rate = 0.9,
    eligible_per_unit = rep(c(0.41, 0.27), each = 4)
  )
  targets <- data.frame(name = c("mcpr", "stunting"),
                        domain = "region", cv = c(0.08, 0.10))
  list(frame = frame, measures = measures, targets = targets)
}

.fold_rate <- function(measures) {
  measures$resp_rate <- measures$resp_rate * measures$eligible_per_unit
  measures$eligible_per_unit <- NULL
  measures
}

test_that("without a correction the joint plan equals the folded rate", {
  z <- .elig_joint()
  elig <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets,
                  fpc = "none")
  folded <- n_alloc_joint(z$frame, measures = .fold_rate(z$measures),
                    targets = z$targets, fpc = "none")
  expect_equal(elig$detail$n, folded$detail$n)
  expect_identical(elig$detail$n_1_int, folded$detail$n_1_int)

  default <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets)
  folded_default <- n_alloc_joint(z$frame, measures = .fold_rate(z$measures),
                            targets = z$targets)
  expect_lt(sum(default$detail$n), sum(folded_default$detail$n))
  expect_equal(default$detail$n, folded_default$detail$n, tolerance = 0.01)
})

test_that("every variance term reads the eligible population", {
  fr <- data.frame(stratum = c("a", "b", "c"), N = c(5000, 8000, 3000))
  ms <- data.frame(stratum = fr$stratum, name = "y",
                   p = c(0.3, 0.2, 0.4), resp_rate = 0.8,
                   eligible_per_unit = c(0.5, 1.2, 0.9))
  tg <- data.frame(name = "y", cv = 0.05)
  elig <- n_alloc_joint(fr, measures = ms, targets = tg)
  scaled <- n_alloc_joint(
    transform(fr, N = N * ms$eligible_per_unit),
    measures = transform(ms, resp_rate = resp_rate * eligible_per_unit,
                         eligible_per_unit = NULL),
    targets = tg
  )
  expect_equal(elig$detail$n, scaled$detail$n)
  expect_identical(elig$detail$N, fr$N)
  expect_equal(elig$detail$weight, fr$N / elig$detail$n)
})

test_that("a measures row beats the frame column, which beats the scalar", {
  z <- .elig_joint()
  no_col <- transform(z$measures, eligible_per_unit = NULL)
  by_scalar <- n_alloc_joint(z$frame, measures = no_col, targets = z$targets,
                       eligible_per_unit = 0.41)
  by_frame <- n_alloc_joint(
    transform(z$frame, eligible_per_unit = c(0.41, 0.41, 0.41, NA)),
    measures = no_col, targets = z$targets, eligible_per_unit = 2
  )
  frame_spelled <- n_alloc_joint(
    transform(z$frame, eligible_per_unit = c(0.41, 0.41, 0.41, 2)),
    measures = no_col, targets = z$targets
  )
  expect_equal(by_frame$detail$n, frame_spelled$detail$n)
  expect_false(isTRUE(all.equal(by_frame$detail$n, by_scalar$detail$n)))
  by_row <- n_alloc_joint(
    transform(z$frame, eligible_per_unit = 2),
    measures = transform(no_col, eligible_per_unit = c(rep(0.41, 7), NA)),
    targets = z$targets
  )
  row_spelled <- n_alloc_joint(
    z$frame, measures = transform(no_col, eligible_per_unit = c(rep(0.41, 7), 2)),
    targets = z$targets
  )
  expect_equal(by_row$detail$n, row_spelled$detail$n)
  expect_error(
    n_alloc_joint(transform(z$frame, eligible_per_unit = 0.41),
                  measures = no_col, targets = z$targets,
                  eligible_per_unit = 2),
    "scalar argument would be ignored"
  )
  expect_error(
    n_alloc_joint(z$frame,
                  measures = transform(no_col, eligible_per_unit = 0.41),
                  targets = z$targets, eligible_per_unit = 3),
    "scalar argument would be ignored"
  )
  expect_identical(by_scalar$params$eligible_per_unit, 0.41)
  expect_null(n_alloc_joint(z$frame, measures = no_col,
                      targets = z$targets)$params$eligible_per_unit)
  per_indicator <- n_alloc_joint(z$frame, measures = z$measures,
                           targets = z$targets)
  common <- n_alloc_joint(z$frame, measures = no_col, targets = z$targets,
                    eligible_per_unit = 0.27)
  expect_false(isTRUE(all.equal(per_indicator$detail$n, common$detail$n)))
  expect_no_error(n_alloc_joint(z$frame,
                          measures = transform(z$measures, eligible_per_unit = 3.2),
                          targets = z$targets))
  expect_error(n_alloc_joint(z$frame,
                       measures = transform(z$measures, eligible_per_unit = 0),
                       targets = z$targets),
               "'measures\\$eligible_per_unit' must contain positive finite")
})

test_that("joint round trips, predict() and the summary keep the rate", {
  z <- .elig_joint()
  x <- n_alloc_joint(z$frame, measures = z$measures, targets = z$targets)
  ev <- prec_alloc_joint(x)
  expect_equal(ev$detail$.achieved[ev$detail$.binding],
               ev$detail$.target[ev$detail$.binding])
  expect_equal(n_alloc_joint(ev)$detail$n, x$detail$n, tolerance = 1e-5)

  xs <- n_alloc_joint(z$frame, measures = transform(z$measures, eligible_per_unit = NULL),
                targets = z$targets, eligible_per_unit = 0.41)
  expect_identical(prec_alloc_joint(xs)$params$eligible_per_unit, 0.41)
  expect_equal(n_alloc_joint(prec_alloc_joint(xs))$detail$n, xs$detail$n,
               tolerance = 1e-5)
  swept <- predict(xs, data.frame(take_2 = 20))
  expect_equal(swept$n, sum(xs$detail$n_int), tolerance = 0.02)

  summarized <- capture.output(summary(x))
  expect_true(any(grepl("Elig\\.", summarized)))
  plain <- n_alloc_joint(z$frame,
                   measures = transform(z$measures, eligible_per_unit = NULL),
                   targets = z$targets)
  expect_false(any(grepl("Elig\\.", capture.output(summary(plain)))))
})

test_that("classic allocation supports eligible_per_unit and a register refuses it", {
  classic <- data.frame(N = c(100, 200), sd = c(2, 3), mean = c(10, 12))
  a_arg <- n_alloc(classic, n = 50, eligible_per_unit = 0.4)
  a_col <- n_alloc(transform(classic, eligible_per_unit = 0.4), n = 50)
  a_plan <- n_alloc(classic, n = 50, plan = svyplan(eligible_per_unit = 0.4))
  expect_equal(a_arg$se, a_col$se)
  expect_equal(a_arg$se, a_plan$se)
  expect_equal(prec_alloc(a_arg)$se, a_arg$se)
  expect_equal(
    prec_alloc(classic, n = c(20, 30), eligible_per_unit = 0.4)$se,
    prec_alloc(transform(classic, N = N * 0.4), n = c(20, 30), resp_rate = 0.4)$se
  )
  z <- .elig_joint()
  reg <- data.frame(stratum = rep(z$frame$stratum, each = 2),
                    N = rep(z$frame$N / 2, each = 2))
  expect_error(
    n_alloc_joint(transform(z$frame, N_1 = NULL), measures = z$measures,
            targets = z$targets, psu = reg),
    "not available with a PSU register"
  )
})

test_that("the substitution holds for moe targets and the staged correction", {
  fr <- data.frame(stratum = c("a", "b", "c"), N = c(5000, 8000, 3000))
  ms <- data.frame(stratum = fr$stratum, name = "y",
                   p = c(0.3, 0.2, 0.4), resp_rate = 0.8,
                   eligible_per_unit = c(0.5, 1.2, 0.9))
  tg <- data.frame(name = "y", moe = 0.03)
  scaled_fr <- transform(fr, N = N * ms$eligible_per_unit)
  folded <- transform(ms, resp_rate = resp_rate * eligible_per_unit,
                      eligible_per_unit = NULL)
  elig <- n_alloc_joint(fr, measures = ms, targets = tg)
  scaled <- n_alloc_joint(scaled_fr, measures = folded, targets = tg)
  expect_equal(elig$detail$n, scaled$detail$n)
  expect_equal(prec_alloc_joint(elig)$detail$.se, prec_alloc_joint(scaled)$detail$.se)

  cl <- transform(fr, N_1 = c(100, 160, 60), take_2 = 12,
                  cost_1 = 300, cost_2 = 20)
  cms <- transform(ms, resp_rate = 0.9, eligible_per_unit = c(0.6, 1.1, 0.8),
                   icc_1 = 0.06)
  ctg <- data.frame(name = "y", cv = 0.06)
  staged <- n_alloc_joint(cl, measures = cms, targets = ctg, fpc = "stage")
  staged_scaled <- n_alloc_joint(
    transform(cl, N = N * cms$eligible_per_unit),
    measures = transform(cms, resp_rate = resp_rate * eligible_per_unit,
                         eligible_per_unit = NULL),
    targets = ctg, fpc = "stage"
  )
  expect_equal(staged$detail$n, staged_scaled$detail$n)
})

test_that("a direct joint assessment takes the scalar rate", {
  z <- .elig_joint()
  no_col <- transform(z$measures, eligible_per_unit = NULL)
  n_h <- c(700, 1000, 740, 1160)
  scalar <- prec_alloc_joint(z$frame, n = n_h, measures = no_col,
                       targets = z$targets, eligible_per_unit = 0.41)
  column <- prec_alloc_joint(z$frame, n = n_h,
                       measures = transform(no_col, eligible_per_unit = 0.41),
                       targets = z$targets)
  none <- prec_alloc_joint(z$frame, n = n_h, measures = no_col,
                     targets = z$targets)
  expect_equal(scalar$detail$.achieved, column$detail$.achieved)
  expect_false(isTRUE(all.equal(scalar$detail$.achieved,
                                none$detail$.achieved)))
  expect_identical(scalar$params$eligible_per_unit, 0.41)
})

test_that("the substitution holds through three stages", {
  fr <- data.frame(
    stratum = c("a", "b"), N = c(60000, 90000), N_1 = c(60, 90),
    N_2 = c(6000, 9000), take_2 = 4, take_3 = 3,
    cost_1 = 500, cost_2 = 60, cost_3 = 15
  )
  ms <- data.frame(stratum = fr$stratum, name = "y", p = c(0.3, 0.25),
                   icc_1 = 0.05, icc_2 = 0.2, resp_rate = 0.85,
                   eligible_per_unit = c(0.7, 1.15))
  tg <- data.frame(name = "y", cv = 0.08)
  for (fpc in c("unit", "stage")) {
    elig <- n_alloc_joint(fr, measures = ms, targets = tg, fpc = fpc)
    scaled <- n_alloc_joint(
      transform(fr, N = N * ms$eligible_per_unit),
      measures = transform(ms, resp_rate = resp_rate * eligible_per_unit,
                           eligible_per_unit = NULL),
      targets = tg, fpc = fpc
    )
    expect_equal(elig$detail$n, scaled$detail$n)
  }
})

test_that("an eligible population of one unit or less is refused", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    expect_error(
      n_prop(0.3, cv = 0.1, N = 10, eligible_per_unit = 0.1, method = method),
      "'N * eligible_per_unit', the eligible population, must be greater than 1",
      fixed = TRUE
    )
    expect_error(
      prec_prop(0.3, n = 5, N = 10, eligible_per_unit = 0.05, method = method),
      "must be greater than 1", fixed = TRUE
    )
  }
  expect_error(prec_prop(n = 5, cv = 0.2, N = 10, eligible_per_unit = 0.1),
               "must be greater than 1", fixed = TRUE)
  expect_true(is.finite(n_prop(0.3, cv = 0.1, N = 10, eligible_per_unit = 0.2)$n))
  expect_true(is.finite(n_prop(0.3, cv = 0.1, eligible_per_unit = 0.01)$n))
  fr <- data.frame(stratum = c("A", "B"), N = c(10, 2000))
  expect_error(
    n_alloc_joint(fr, measures = data.frame(stratum = fr$stratum, name = "y",
                                      p = c(0.3, 0.5)),
            targets = data.frame(name = "y", cv = 0.05), eligible_per_unit = 0.05),
    "'N * eligible_per_unit', the eligible population", fixed = TRUE
  )
})

test_that("the default frontier starts at the eligible floor", {
  fr <- data.frame(stratum = "A", N = 100000)
  ms <- data.frame(stratum = "A", name = "y", p = 0.3)
  tg <- data.frame(name = "y", cv = 0.05)
  floor_cost <- n_alloc_joint(fr, measures = ms, targets = tg,
                        eligible_per_unit = 0.1)$operational$cost
  fit <- n_alloc_joint(fr, measures = ms, targets = tg, objective = "y",
                 budget = 15000, eligible_per_unit = 0.1)
  grid <- svyplan:::.budget_grid(fit, 5)
  expect_equal(grid[1], floor_cost)
  rows <- expect_no_warning(predict(fit, data.frame(budget = grid)))
  expect_true(all(rows$.feasible))
})

## Means and the power family

test_that("a mean is sized on eligible respondents in the eligible population", {
  z <- qnorm(0.975)
  for (N in c(Inf, 5000)) {
    fit <- n_mean(var = 4, moe = 0.2, N = N, deff = 1.5, resp_rate = 0.9,
                  eligibility_rate = 0.95, eligible_per_unit = 2.6)
    N_elig <- N * 0.95 * 2.6
    net <- z^2 * 1.5 * 4 / (0.2^2 + z^2 * 1.5 * 4 / N_elig)
    expect_equal(fit$n, net / (0.9 * 0.95 * 2.6), info = format(N))
    back <- prec_mean(var = 4, n = fit$n, N = N, deff = 1.5, resp_rate = 0.9,
                      eligibility_rate = 0.95, eligible_per_unit = 2.6)
    expect_equal(back$moe, 0.2, info = format(N))
    expect_equal(prec_mean(fit)$moe, 0.2, info = format(N))
    expect_equal(n_mean(back)$n, fit$n, info = format(N))
  }
  expect_identical(fit$params$eligible_per_unit, 2.6)
  expect_null(n_mean(var = 4, moe = 0.2)$params$eligible_per_unit)
  expect_error(n_mean(var = 4, moe = 0.2, eligible_per_unit = 0),
               "'eligible_per_unit' must be a positive finite number")
  expect_error(prec_mean(var = 4, n = 1, N = 2, eligible_per_unit = 0.4),
               "the eligible population, must be greater than 1")
})

test_that("a power solve counts eligible respondents and round-trips", {
  cases <- list(
    function(...) power_prop(p1 = 0.3, p2 = 0.35, N = 4000, overlap = 0.5,
                             overlap_cor = 0.4, ...),
    function(...) power_mean(var = 100, effect = 2, N = 4000, overlap = 0.5,
                             overlap_cor = 0.4, ...),
    function(...) power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                            outcome = "prop", N = 6000, overlap = 0.5,
                            overlap_cor = 0.4, ...)
  )
  for (call in cases) {
    fit <- call(resp_rate = 0.9, eligible_per_unit = 2.6)
    expect_identical(fit$params$eligible_per_unit, 2.6)
    expect_equal(
      call(resp_rate = 0.9, eligible_per_unit = 2.6, n = fit$n,
           power = NULL)$power,
      0.8, tolerance = 1e-6, info = fit$type
    )
  }
  # Without a finite population the count only divides the issued size.
  expect_equal(power_prop(p1 = 0.3, p2 = 0.35, eligible_per_unit = 2.6)$n,
               power_prop(p1 = 0.3, p2 = 0.35)$n / 2.6)
  expect_equal(power_mean(var = 100, effect = 2, eligible_per_unit = 0.45)$n,
               power_mean(var = 100, effect = 2)$n / 0.45)
  pl <- svyplan(resp_rate = 0.9, eligible_per_unit = 0.45)
  expect_equal((pl |> power_prop(p1 = 0.3, p2 = 0.35))$n,
               power_prop(p1 = 0.3, p2 = 0.35, resp_rate = 0.9,
                          eligible_per_unit = 0.45)$n)
})

test_that("a count below 1 scales the population as eligibility does", {
  # Both multiply the rate and the population in scope, so at 0.5 they agree
  # exactly, finite population and overlap included.
  sized <- list(
    function(...) n_mean(var = 4, moe = 0.2, N = 3000, resp_rate = 0.9, ...),
    function(...) power_prop(p1 = 0.3, p2 = 0.35, N = 3000, overlap = 0.5,
                             overlap_cor = 0.4, resp_rate = 0.9, ...),
    function(...) power_mean(var = 100, effect = 2, N = 3000, overlap = 0.5,
                             overlap_cor = 0.4, resp_rate = 0.9, ...),
    function(...) power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                            outcome = "prop", N = 3000, overlap = 0.5,
                            overlap_cor = 0.4, resp_rate = 0.9, ...)
  )
  for (call in sized) {
    expect_equal(call(eligible_per_unit = 0.5)$n,
                 call(eligibility_rate = 0.5)$n)
  }
  expect_equal(
    prec_mean(var = 4, n = 900, N = 3000, eligible_per_unit = 0.5)$se,
    prec_mean(var = 4, n = 900, N = 3000, eligibility_rate = 0.5)$se
  )
})

## The repeated-survey and ratio functions check the eligible population

test_that("repeated-survey and ratio functions refuse one eligible unit or less", {
  msg <- "'N * eligibility_rate', the eligible population, must be greater than 1"
  expect_error(prec_pooled(p = 0.3, n = 5, occasions = 4, N = 10,
                           eligibility_rate = 0.1), msg, fixed = TRUE)
  expect_error(prec_pooled(var = 1, n = 5, occasions = 4, N = 10,
                           eligibility_rate = 0.1), msg, fixed = TRUE)
  expect_error(n_pooled(p = 0.3, moe = 0.1, occasions = 4, N = 10,
                        eligibility_rate = 0.05), msg, fixed = TRUE)
  expect_error(prec_change(p = c(0.3, 0.4), n = 5, N = 10,
                           eligibility_rate = 0.05), msg, fixed = TRUE)
  expect_error(prec_change(p = c(0.3, 0.4), n = 5, N = c(100, 5),
                           eligibility_rate = 0.2), msg, fixed = TRUE)
  expect_error(n_change(p = c(0.3, 0.4), moe = 0.1, N = c(100, 4),
                        eligibility_rate = 0.25), msg, fixed = TRUE)
  ratio <- list(r = 1, cv_num = 0.5, cv_den = 0.4, component_cor = 0.5)
  expect_error(do.call(n_ratio, c(ratio, cv = 0.1, N = 2,
                                  eligibility_rate = 0.2)), msg, fixed = TRUE)
  expect_error(do.call(prec_ratio, c(ratio, n = 1, N = 2,
                                     eligibility_rate = 0.2)), msg, fixed = TRUE)
  expect_true(is.finite(n_change(p = c(0.3, 0.4), moe = 0.1,
                                 N = c(100, 200), eligibility_rate = 0.5)$n))
  expect_true(is.finite(n_pooled(p = 0.3, moe = 0.05, occasions = 3, N = 1e4,
                                 eligibility_rate = 0.5)$n))
  expect_error(prec_mean(var = 1, n = 5, N = 10, eligibility_rate = 0.1),
               "the eligible population, must be greater than 1")
})

test_that("an unattainable size names the eligible units per unit", {
  expect_error(n_prop(p = 0.3, moe = 0.001, N = 100, eligible_per_unit = 2,
                      resp_rate = 0.4),
               "about 80 respondents at resp_rate = 0.4 with eligible_per_unit = 2")
})
