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
                                     eligible_rate = 0.41, method = method),
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
              eligible_rate = 0.41, N = 20000)
  expect_equal(x$n, .wald_cv_n(0.28, 0.05, 20000 * 0.41, 2.5) / (0.9 * 0.41))
  folded <- n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9 * 0.41,
                   N = 20000)
  expect_gt(folded$n, x$n)
})

test_that("a rate above 1 is accepted and solved on the same scale", {
  x <- n_prop(p = 0.2, cv = 0.08, eligible_rate = 2.1, N = 5000)
  expect_equal(x$n, .wald_cv_n(0.2, 0.08, 5000 * 2.1, 1) / 2.1)
  expect_lt(x$n, n_prop(p = 0.2, cv = 0.08, N = 5000)$n)
})

test_that("cases count eligible respondents", {
  x <- n_prop(p = 0.02, moe = 0.02, resp_rate = 0.9, eligible_rate = 0.6,
              min_cases = 40)
  expect_identical(x$binding, "min_cases")
  expect_equal(x$n, 40 / (0.02 * 0.9 * 0.6))
  expect_equal(x$expected_cases, 40)
})

test_that("prec_prop forms the same precision as n_prop targets", {
  ev <- prec_prop(p = 0.28, n = 6000, deff = 2.5, resp_rate = 0.9,
                  eligible_rate = 0.41, N = 20000)
  n_net <- 6000 * 0.9 * 0.41
  N_elig <- 20000 * 0.41
  n_eff <- n_net * (N_elig - 1) / (2.5 * (N_elig - n_net))
  expect_equal(ev$se, sqrt(0.28 * 0.72 / n_eff))

  floor_p <- prec_prop(n = 6000, cv = 0.1, resp_rate = 0.9,
                       eligible_rate = 0.41, N = 20000)
  expect_equal(prec_prop(p = floor_p$params$p, n = 6000, resp_rate = 0.9,
                         eligible_rate = 0.41, N = 20000)$cv, 0.1)
})

test_that("the rate survives every round trip", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    x <- n_prop(p = 0.28, moe = 0.03, deff = 2.5, resp_rate = 0.9,
                eligible_rate = 0.41, N = 20000, method = method)
    expect_equal(x$moe, 0.03)
    ev <- prec_prop(x)
    expect_equal(ev$moe, 0.03)
    expect_identical(ev$params$eligible_rate, 0.41)
    expect_equal(n_prop(ev)$n, x$n)
  }
  x <- n_prop(p = 0.28, cv = 0.05, eligible_rate = 1.6, N = 8000)
  swept <- predict(x, data.frame(eligible_rate = c(0.4, 1.6, 2.5)))
  expect_equal(swept$cv, rep(0.05, 3))
  expect_equal(swept$n[2], x$n)
  ev <- prec_prop(x)
  expect_equal(predict(ev, data.frame(eligible_rate = 1.6))$cv, 0.05)
  expect_equal(predict(x, data.frame(cv = c(0.05, 0.04)))$n[1], x$n)
  expect_equal(predict(ev, data.frame(n = x$n))$cv, 0.05)
})

test_that("a profile carries the rate and refuses a function without it", {
  pl <- svyplan(eligible_rate = 0.41, resp_rate = 0.9)
  expect_identical(
    n_prop(p = 0.28, cv = 0.05, plan = pl)$n,
    n_prop(p = 0.28, cv = 0.05, resp_rate = 0.9, eligible_rate = 0.41)$n
  )
  expect_identical(
    prec_prop(p = 0.28, n = 6000, plan = pl)$params$eligible_rate, 0.41
  )
  expect_error(n_mean(var = 4, mu = 10, cv = 0.05, plan = pl),
               "profile sets 'eligible_rate', which n_mean\\(\\) does not take")
  expect_error(n_multi(data.frame(name = "a", p = 0.3, moe = 0.05), plan = pl),
               "which n_multi\\(\\) does not take")
  piped <- tryCatch(
    pl |> n_multi(data.frame(name = "a", p = 0.3, moe = 0.05)),
    error = conditionMessage
  )
  expect_identical(
    piped,
    paste("the profile sets 'eligible_rate', which n_multi() does not take.",
          "Drop it from the profile for this call, or fold it into",
          "'resp_rate' when it is at most 1")
  )
  expect_error(pl |> n_mean(var = 4, mu = 10, cv = 0.05),
               "which n_mean\\(\\) does not take")
  expect_no_error(n_mean(var = 4, mu = 10, cv = 0.05,
                         plan = svyplan(eligible_rate = 1)))
  expect_error(svyplan(eligible_rate = 0), "positive finite")
})

test_that("the rate is validated and stays out of a default result", {
  for (bad in list(0, -1, Inf, NA_real_, c(0.4, 0.5), "0.4")) {
    expect_error(n_prop(p = 0.3, moe = 0.05, eligible_rate = bad),
                 "'eligible_rate' must be a positive finite number")
    expect_error(prec_prop(p = 0.3, n = 400, eligible_rate = bad),
                 "'eligible_rate' must be a positive finite number")
  }
  expect_null(n_prop(p = 0.3, moe = 0.05)$params$eligible_rate)
  expect_null(prec_prop(p = 0.3, n = 400)$params$eligible_rate)
  expect_error(n_prop(p = 0.3, cv = 0.01, resp_rate = 0.5,
                      eligible_rate = 0.2, N = 1000),
               "exceeds the population")
})

test_that("the print counts eligible respondents", {
  x <- n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9,
              eligible_rate = 0.41)
  out <- capture.output(print(x))
  expect_match(out[2], "eligible respondents: 2572")
  expect_match(out[2], "eligible_rate = 0.41")
  expect_match(capture.output(print(prec_prop(x)))[2],
               "eligible respondents: 2572, eligible_rate = 0.41")
  expect_false(any(grepl("eligible",
                         capture.output(print(n_prop(p = 0.3, moe = 0.05))))))
})

test_that("the interval is as wide as the reported margin", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    for (rate in c(0.25, 2.6)) {
      x <- n_prop(p = 0.3, moe = 0.04, eligible_rate = rate, N = 20000,
                  deff = 1.5, resp_rate = 0.8, method = method)
      prec <- prec_prop(p = 0.3, n = x$n, eligible_rate = rate, N = 20000,
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
    list(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, n_psu = 80),
    list(stage_cost = c(500, 50), icc = 0.05, budget = 1e5, n_per_psu = 25),
    list(stage_cost = c(500, 100, 20), icc = c(0.03, 0.1), cv = 0.05)
  )
  for (args in calls) {
    elig <- do.call(n_cluster, c(args, resp_rate = 0.9, eligible_rate = 0.41))
    folded <- do.call(n_cluster, c(args, resp_rate = 0.9 * 0.41))
    expect_equal(elig$n, folded$n)
    expect_equal(elig$operational$n, folded$operational$n)
    expect_identical(elig$params$resp_rate, 0.9)
  }
})

test_that("the household take is the eligible-person optimum converted", {
  for (e in c(0.41, 2.5)) {
    x <- n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05,
                   resp_rate = 0.9, eligible_rate = e)
    cost_per_eligible <- 50 / (e * 0.9)
    persons <- sqrt(500 / cost_per_eligible * (1 - 0.07) / 0.07)
    expect_equal(unname(x$n["n_per_psu"]), persons / (e * 0.9))
  }
})

test_that("prec_cluster reads the eligible respondents per cluster", {
  ev <- prec_cluster(n = c(600, 20), icc = 0.07, resp_rate = 0.9,
                     eligible_rate = 0.41)
  b <- 20 * 0.9 * 0.41
  expect_equal(ev$cv, sqrt((1 + 0.07 * (b - 1)) / (600 * b)))
  expect_identical(ev$params$eligible_rate, 0.41)
  expect_identical(ev$params$resp_rate, 0.9)
})

test_that("the cluster rate survives every round trip", {
  for (args in list(
    list(stage_cost = c(500, 50), icc = 0.07, cv = 0.05),
    list(stage_cost = c(500, 100, 20), icc = c(0.03, 0.1), budget = 1e5)
  )) {
    x <- do.call(n_cluster, c(args, resp_rate = 0.9, eligible_rate = 1.8))
    ev <- prec_cluster(x)
    expect_equal(ev$cv, x$cv)
    expect_identical(ev$params$eligible_rate, 1.8)
    expect_equal(n_cluster(ev)$n, x$n)
  }
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05,
                 eligible_rate = 0.41)
  expect_equal(n_cluster(prec_cluster(x))$n, x$n)
  swept <- predict(x, data.frame(eligible_rate = c(0.41, 1.2)))
  expect_equal(swept$cv, c(0.05, 0.05))
  expect_equal(swept$n_psu[1], unname(x$n["n_psu"]))
  by_cv <- predict(x, data.frame(cv = c(0.05, 0.04)))
  expect_equal(by_cv$n_per_psu[1], unname(x$n["n_per_psu"]))
})

test_that("an indicators table refuses an eligible_rate column", {
  msg <- "does not take an 'eligible_rate' column"
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.5), moe = 0.05,
                    eligible_rate = 0.41)
  expect_error(n_multi(ind), msg)
  expect_error(
    n_multi(transform(rbind(ind, ind), dom = rep(1:2, each = 2)),
            domains = "dom"),
    msg
  )
  expect_error(prec_multi(transform(ind, moe = NULL, n = 500)), msg)
  expect_error(
    n_cluster(indicators = transform(ind, icc_psu = 0.05),
              stage_cost = c(500, 50)),
    msg
  )
  expect_error(
    prec_cluster(indicators = transform(ind, moe = NULL, icc_psu = 0.05,
                                        n = 50, n_per_psu = 10)),
    msg
  )
  expect_no_error(
    n_multi(transform(ind, eligible_rate = NULL, grp = "x"))
  )
})

test_that("indicators mode refuses the rate and the take check names it", {
  ind <- data.frame(name = "a", p = 0.3, cv = 0.1, icc_psu = 0.05)
  expect_error(
    n_cluster(indicators = ind, stage_cost = c(500, 50), eligible_rate = 0.4),
    "applies to the single-indicator call"
  )
  expect_error(
    prec_cluster(indicators = transform(ind, n = 50, n_per_psu = 10),
                 eligible_rate = 0.4),
    "applies to the single-indicator call"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05, n_per_psu = 2,
              resp_rate = 0.9, eligible_rate = 0.4),
    "at 'resp_rate \\* eligible_rate' = 0.36"
  )
  expect_error(
    prec_cluster(n = c(50, 2), icc = 0.05, resp_rate = 0.9,
                 eligible_rate = 0.4),
    "at 'resp_rate \\* eligible_rate' = 0.36"
  )
})

test_that("a profile carries the rate to the cluster functions", {
  pl <- svyplan(eligible_rate = 0.41, resp_rate = 0.9)
  expect_equal(
    n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05, plan = pl)$n,
    n_cluster(stage_cost = c(500, 50), icc = 0.07, cv = 0.05,
              resp_rate = 0.9, eligible_rate = 0.41)$n
  )
  expect_identical(
    prec_cluster(n = c(60, 20), icc = 0.07, plan = pl)$params$eligible_rate,
    0.41
  )
})

test_that("the cluster print counts eligible respondents", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                 resp_rate = 0.9, eligible_rate = 0.41)
  out <- capture.output(print(x))
  expect_match(out[2], "eligible respondents: 558")
  expect_match(out[3], "resp_rate = 0.90, eligible_rate = 0.41")
  expect_match(capture.output(print(prec_cluster(x)))[2],
               "eligible respondents: 561, eligible_rate = 0.41")
  plain <- capture.output(print(n_cluster(stage_cost = c(500, 50),
                                          icc = 0.05, cv = 0.05)))
  expect_false(any(grepl("eligible", plain)))
  expect_null(n_cluster(stage_cost = c(500, 50), icc = 0.05,
                        cv = 0.05)$params$eligible_rate)
})

## Eligible units per sampled unit: joint allocation

.elig_joint <- function() {
  frame <- data.frame(
    stratum = c("A-u", "A-r", "B-u", "B-r"),
    region = rep(c("A", "B"), each = 2),
    N = c(4e5, 6e5, 3e5, 5e5), n_per_psu = 20, cost_psu = 400,
    cost_ssu = 30, N_psu = c(2000, 3000, 1500, 2500)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("mcpr", "stunting"), each = 4),
    p = c(0.30, 0.24, 0.28, 0.22, 0.20, 0.30, 0.22, 0.33),
    icc_psu = rep(c(0.07, 0.05), each = 4), resp_rate = 0.9,
    eligible_rate = rep(c(0.41, 0.27), each = 4)
  )
  targets <- data.frame(name = c("mcpr", "stunting"),
                        domain = "region", cv = c(0.08, 0.10))
  list(frame = frame, measures = measures, targets = targets)
}

.fold_rate <- function(measures) {
  measures$resp_rate <- measures$resp_rate * measures$eligible_rate
  measures$eligible_rate <- NULL
  measures
}

test_that("without a correction the joint plan equals the folded rate", {
  z <- .elig_joint()
  elig <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  fpc = "none")
  folded <- n_alloc(z$frame, measures = .fold_rate(z$measures),
                    targets = z$targets, fpc = "none")
  expect_equal(elig$detail$n, folded$detail$n)
  expect_identical(elig$detail$n_psu_int, folded$detail$n_psu_int)

  default <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  folded_default <- n_alloc(z$frame, measures = .fold_rate(z$measures),
                            targets = z$targets)
  expect_lt(sum(default$detail$n), sum(folded_default$detail$n))
  expect_equal(default$detail$n, folded_default$detail$n, tolerance = 0.01)
})

test_that("every variance term reads the eligible population", {
  fr <- data.frame(stratum = c("a", "b", "c"), N = c(5000, 8000, 3000))
  ms <- data.frame(stratum = fr$stratum, name = "y",
                   p = c(0.3, 0.2, 0.4), resp_rate = 0.8,
                   eligible_rate = c(0.5, 1.2, 0.9))
  tg <- data.frame(name = "y", cv = 0.05)
  elig <- n_alloc(fr, measures = ms, targets = tg)
  scaled <- n_alloc(
    transform(fr, N = N * ms$eligible_rate),
    measures = transform(ms, resp_rate = resp_rate * eligible_rate,
                         eligible_rate = NULL),
    targets = tg
  )
  expect_equal(elig$detail$n, scaled$detail$n)
  expect_identical(elig$detail$N, fr$N)
  expect_equal(elig$detail$weight, fr$N / elig$detail$n)
})

test_that("a measures row beats the frame column, which beats the scalar", {
  z <- .elig_joint()
  no_col <- transform(z$measures, eligible_rate = NULL)
  by_scalar <- n_alloc(z$frame, measures = no_col, targets = z$targets,
                       eligible_rate = 0.41)
  by_frame <- n_alloc(transform(z$frame, eligible_rate = 0.41),
                      measures = no_col, targets = z$targets,
                      eligible_rate = 2)
  by_row <- n_alloc(transform(z$frame, eligible_rate = 2),
                    measures = transform(no_col, eligible_rate = 0.41),
                    targets = z$targets, eligible_rate = 3)
  expect_equal(by_frame$detail$n, by_scalar$detail$n)
  expect_equal(by_row$detail$n, by_scalar$detail$n)
  expect_identical(by_scalar$params$eligible_rate, 0.41)
  expect_null(n_alloc(z$frame, measures = no_col,
                      targets = z$targets)$params$eligible_rate)
  per_indicator <- n_alloc(z$frame, measures = z$measures,
                           targets = z$targets)
  common <- n_alloc(z$frame, measures = no_col, targets = z$targets,
                    eligible_rate = 0.27)
  expect_false(isTRUE(all.equal(per_indicator$detail$n, common$detail$n)))
  expect_no_error(n_alloc(z$frame,
                          measures = transform(z$measures, eligible_rate = 3.2),
                          targets = z$targets))
  expect_error(n_alloc(z$frame,
                       measures = transform(z$measures, eligible_rate = 0),
                       targets = z$targets),
               "'measures\\$eligible_rate' must contain positive finite")
})

test_that("joint round trips, predict() and the summary keep the rate", {
  z <- .elig_joint()
  x <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  ev <- prec_alloc(x)
  expect_equal(ev$detail$.achieved[ev$detail$.binding],
               ev$detail$.target[ev$detail$.binding])
  expect_equal(n_alloc(ev)$detail$n, x$detail$n, tolerance = 1e-5)

  xs <- n_alloc(z$frame, measures = transform(z$measures, eligible_rate = NULL),
                targets = z$targets, eligible_rate = 0.41)
  expect_identical(prec_alloc(xs)$params$eligible_rate, 0.41)
  expect_equal(n_alloc(prec_alloc(xs))$detail$n, xs$detail$n,
               tolerance = 1e-5)
  swept <- predict(xs, data.frame(n_per_psu = 20))
  expect_equal(swept$n, sum(xs$detail$n_int), tolerance = 0.02)

  summarized <- capture.output(summary(x))
  expect_true(any(grepl("Elig\\.", summarized)))
  plain <- n_alloc(z$frame,
                   measures = transform(z$measures, eligible_rate = NULL),
                   targets = z$targets)
  expect_false(any(grepl("Elig\\.", capture.output(summary(plain)))))
})

test_that("classic allocation and a register refuse the rate", {
  classic <- data.frame(N = c(100, 200), sd = 1)
  expect_error(n_alloc(classic, n = 50, eligible_rate = 0.4),
               "applies to joint allocation")
  expect_error(n_alloc(transform(classic, eligible_rate = 0.4), n = 50),
               "applies to joint allocation")
  expect_error(prec_alloc(classic, n = c(20, 30), eligible_rate = 0.4),
               "applies to joint allocation")
  expect_error(n_alloc(classic, n = 50,
                       plan = svyplan(eligible_rate = 0.4)),
               "applies to joint allocation")
  z <- .elig_joint()
  reg <- data.frame(stratum = rep(z$frame$stratum, each = 2),
                    N = rep(z$frame$N / 2, each = 2))
  expect_error(
    n_alloc(transform(z$frame, N_psu = NULL), measures = z$measures,
            targets = z$targets, psu = reg),
    "not available with a PSU register"
  )
})

test_that("the substitution holds for moe targets and the staged correction", {
  fr <- data.frame(stratum = c("a", "b", "c"), N = c(5000, 8000, 3000))
  ms <- data.frame(stratum = fr$stratum, name = "y",
                   p = c(0.3, 0.2, 0.4), resp_rate = 0.8,
                   eligible_rate = c(0.5, 1.2, 0.9))
  tg <- data.frame(name = "y", moe = 0.03)
  scaled_fr <- transform(fr, N = N * ms$eligible_rate)
  folded <- transform(ms, resp_rate = resp_rate * eligible_rate,
                      eligible_rate = NULL)
  elig <- n_alloc(fr, measures = ms, targets = tg)
  scaled <- n_alloc(scaled_fr, measures = folded, targets = tg)
  expect_equal(elig$detail$n, scaled$detail$n)
  expect_equal(prec_alloc(elig)$detail$.se, prec_alloc(scaled)$detail$.se)

  cl <- transform(fr, N_psu = c(100, 160, 60), n_per_psu = 12,
                  cost_psu = 300, cost_ssu = 20)
  cms <- transform(ms, resp_rate = 0.9, eligible_rate = c(0.6, 1.1, 0.8),
                   icc_psu = 0.06)
  ctg <- data.frame(name = "y", cv = 0.06)
  staged <- n_alloc(cl, measures = cms, targets = ctg, fpc = "stage")
  staged_scaled <- n_alloc(
    transform(cl, N = N * cms$eligible_rate),
    measures = transform(cms, resp_rate = resp_rate * eligible_rate,
                         eligible_rate = NULL),
    targets = ctg, fpc = "stage"
  )
  expect_equal(staged$detail$n, staged_scaled$detail$n)
})

test_that("a direct joint assessment takes the scalar rate", {
  z <- .elig_joint()
  no_col <- transform(z$measures, eligible_rate = NULL)
  n_h <- c(700, 1000, 740, 1160)
  scalar <- prec_alloc(z$frame, n = n_h, measures = no_col,
                       targets = z$targets, eligible_rate = 0.41)
  column <- prec_alloc(z$frame, n = n_h,
                       measures = transform(no_col, eligible_rate = 0.41),
                       targets = z$targets)
  none <- prec_alloc(z$frame, n = n_h, measures = no_col,
                     targets = z$targets)
  expect_equal(scalar$detail$.achieved, column$detail$.achieved)
  expect_false(isTRUE(all.equal(scalar$detail$.achieved,
                                none$detail$.achieved)))
  expect_identical(scalar$params$eligible_rate, 0.41)
})

test_that("the substitution holds through three stages", {
  fr <- data.frame(
    stratum = c("a", "b"), N = c(60000, 90000), N_psu = c(60, 90),
    N_ssu = c(6000, 9000), n_per_psu = 4, n_per_ssu = 3,
    cost_psu = 500, cost_ssu = 60, cost_tsu = 15
  )
  ms <- data.frame(stratum = fr$stratum, name = "y", p = c(0.3, 0.25),
                   icc_psu = 0.05, icc_ssu = 0.2, resp_rate = 0.85,
                   eligible_rate = c(0.7, 1.15))
  tg <- data.frame(name = "y", cv = 0.08)
  for (fpc in c("unit", "stage")) {
    elig <- n_alloc(fr, measures = ms, targets = tg, fpc = fpc)
    scaled <- n_alloc(
      transform(fr, N = N * ms$eligible_rate),
      measures = transform(ms, resp_rate = resp_rate * eligible_rate,
                           eligible_rate = NULL),
      targets = tg, fpc = fpc
    )
    expect_equal(elig$detail$n, scaled$detail$n)
  }
})

test_that("an eligible population of one unit or less is refused", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    expect_error(
      n_prop(0.3, cv = 0.1, N = 10, eligible_rate = 0.1, method = method),
      "'N * eligible_rate', the eligible population, must be greater than 1",
      fixed = TRUE
    )
    expect_error(
      prec_prop(0.3, n = 5, N = 10, eligible_rate = 0.05, method = method),
      "must be greater than 1", fixed = TRUE
    )
  }
  expect_error(prec_prop(n = 5, cv = 0.2, N = 10, eligible_rate = 0.1),
               "must be greater than 1", fixed = TRUE)
  expect_true(is.finite(n_prop(0.3, cv = 0.1, N = 10, eligible_rate = 0.2)$n))
  expect_true(is.finite(n_prop(0.3, cv = 0.1, eligible_rate = 0.01)$n))
  fr <- data.frame(stratum = c("A", "B"), N = c(10, 2000))
  expect_error(
    n_alloc(fr, measures = data.frame(stratum = fr$stratum, name = "y",
                                      p = c(0.3, 0.5)),
            targets = data.frame(name = "y", cv = 0.05), eligible_rate = 0.05),
    "'N * eligible_rate', the eligible population", fixed = TRUE
  )
})

test_that("the default frontier starts at the eligible floor", {
  fr <- data.frame(stratum = "A", N = 100000)
  ms <- data.frame(stratum = "A", name = "y", p = 0.3)
  tg <- data.frame(name = "y", cv = 0.05)
  floor_cost <- n_alloc(fr, measures = ms, targets = tg,
                        eligible_rate = 0.1)$operational$cost
  fit <- n_alloc(fr, measures = ms, targets = tg, objective = "y",
                 budget = 15000, eligible_rate = 0.1)
  grid <- svyplan:::.budget_grid(fit, 5)
  expect_equal(grid[1], floor_cost)
  rows <- expect_no_warning(predict(fit, data.frame(budget = grid)))
  expect_true(all(rows$.feasible))
})
