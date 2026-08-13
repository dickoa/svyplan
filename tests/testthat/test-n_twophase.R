frame4 <- data.frame(
  stratum   = c("A", "B", "C", "D"),
  N         = c(3500, 2500, 2500, 1500),
  sd        = c(12, 25, 8, 40),
  mean      = c(40, 70, 35, 90),
  unit_cost = c(2, 5, 1, 9)
)

# brute-force the constrained optimum of the VC product
vc_opt <- function(W, S, cost, c_a, A, pinned = rep(FALSE, length(W))) {
  VC <- function(nu) (A + sum(W * S^2 / nu)) * (c_a + sum(cost * W * nu))
  free <- !pinned
  p <- stats::optim(rep(0.5, sum(free)),
                    function(p) { nu <- rep(1, length(W)); nu[free] <- p; VC(nu) },
                    method = "L-BFGS-B", lower = 1e-9, upper = 1)$par
  nu <- rep(1, length(W)); nu[free] <- p; nu
}

test_that("the closed form attains the constrained optimum", {
  W <- frame4$N / sum(frame4$N)
  A <- sum(W * (frame4$mean - sum(W * frame4$mean))^2)
  fit <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  ref <- vc_opt(W, frame4$sd, frame4$unit_cost, 1, A)
  expect_equal(fit$detail$nu, ref, tolerance = 1e-5)
})

test_that("weak stratification truncates strata one at a time", {
  # a batch-pinning implementation returns all ones here and is strictly worse
  flat <- transform(frame4, mean = c(55, 55.4, 55, 55.4))
  W <- flat$N / sum(flat$N)
  A <- sum(W * (flat$mean - sum(W * flat$mean))^2)
  fit <- n_twophase(flat, phase1_cost = 1, budget = 50000)
  ref <- vc_opt(W, flat$sd, flat$unit_cost, 1, A)
  expect_equal(fit$detail$nu, ref, tolerance = 1e-5)
  expect_false(all(fit$detail$nu == 1))
  expect_true(any(fit$detail$take_all))
})

test_that("no between-stratum component still gives an interior optimum", {
  # A = 0 makes the multiplier infinite on the first pass; the tie-break must
  # rank on sd/sqrt(unit_cost), not on nu
  flat <- frame4[, setdiff(names(frame4), "mean")]
  W <- flat$N / sum(flat$N)
  fit <- n_twophase(flat, phase1_cost = 1, budget = 50000, mu = 50)
  ref <- vc_opt(W, flat$sd, flat$unit_cost, 1, 0)
  expect_equal(fit$detail$nu, ref, tolerance = 1e-5)
  expect_false(all(fit$detail$nu == 1))
})

test_that("nonresponse follow-up reproduces the standard optimum", {
  for (p in list(list(c1 = 50, c2 = 200, th = 0.50, n1 = 828, n2 = 293, cv = 0.0382),
                 list(c1 = 75, c2 = 150, th = 0.70, n1 = 949, n2 = 241, cv = 0.10))) {
    nrfu <- data.frame(
      stratum   = c("respondents", "nonrespondents"),
      N         = c(p$th, 1 - p$th),
      sd        = c(1, 1),
      unit_cost = c(0, p$c2),
      take_all  = c(TRUE, FALSE)
    )
    fit <- if (p$cv == 0.10) {
      n_twophase(nrfu, phase1_cost = p$c1, cv = 0.10, mu = 1 / 3)
    } else {
      n_twophase(nrfu, phase1_cost = p$c1, budget = 100000, mu = 1)
    }
    expect_equal(fit$detail$nu[2], sqrt(p$c1 / (p$c2 * p$th)), tolerance = 1e-8)
    expect_equal(round(fit$n[["n_phase1"]]), p$n1, tolerance = 1)
    expect_equal(round(fit$detail$n_issued[2]), p$n2, tolerance = 1)
  }
})

test_that("cv and budget modes are inverses of each other", {
  a <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  b <- n_twophase(frame4, phase1_cost = 1, cv = a$cv)
  expect_equal(b$n[["n_phase1"]], a$n[["n_phase1"]], tolerance = 1e-6)
  expect_equal(b$cost, a$cost, tolerance = 1e-6)
  expect_equal(b$detail$nu, a$detail$nu, tolerance = 1e-8)
})

test_that("prec_twophase inverts n_twophase", {
  fit <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  pr <- prec_twophase(fit)
  expect_s3_class(pr, "svyplan_prec")
  expect_equal(pr$cv, fit$cv, tolerance = 1e-8)
})

test_that("the single-phase comparison is reported and can win", {
  fit <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  expect_true(is.finite(fit$single_phase$n))
  expect_true(is.logical(fit$single_phase$better))
  # an expensive screener with a useless stratifier cannot pay for itself
  useless <- transform(frame4, mean = rep(55, 4))
  lose <- n_twophase(useless, phase1_cost = 50, budget = 50000, mu = 55)
  expect_true(lose$single_phase$better)
  # a cheap screener with strong stratification does
  strong <- transform(frame4, mean = c(10, 200, 5, 400))
  win <- n_twophase(strong, phase1_cost = 0.01, budget = 50000)
  expect_false(win$single_phase$better)
})

test_that("take_all strata are pinned and cost nothing extra when free", {
  pinned <- transform(frame4, take_all = c(FALSE, TRUE, FALSE, FALSE))
  fit <- n_twophase(pinned, phase1_cost = 1, budget = 50000)
  expect_equal(fit$detail$nu[2], 1)
  expect_true(fit$detail$take_all[2])
})

test_that("degenerate strata behave", {
  zero_sd <- transform(frame4, sd = c(12, 0, 8, 40))
  expect_equal(n_twophase(zero_sd, phase1_cost = 1, budget = 50000)$detail$nu[2], 0)
  free_zero_cost <- transform(frame4, unit_cost = c(2, 0, 1, 9))
  expect_equal(n_twophase(free_zero_cost, phase1_cost = 1, budget = 50000)$detail$nu[2], 1)
})

test_that("input validation is specific", {
  expect_error(n_twophase(frame4[0, ], phase1_cost = 1, budget = 1), "non-empty")
  expect_error(n_twophase(frame4, budget = 1), "'phase1_cost' is required")
  expect_error(n_twophase(frame4, phase1_cost = 1), "exactly one of 'cv' or 'budget'")
  expect_error(n_twophase(frame4, phase1_cost = 1, cv = 0.05, budget = 1),
               "exactly one of 'cv' or 'budget'")
  expect_error(n_twophase(frame4[, c("N", "mean")], phase1_cost = 1, budget = 1),
               "must contain a 'sd' column")
  expect_error(n_twophase(transform(frame4, cost = 1), phase1_cost = 1, budget = 1),
               "per-stratum cost column is 'unit_cost'")
  expect_error(n_twophase(transform(frame4, take_all = TRUE), phase1_cost = 1, budget = 1),
               "nothing is left to subsample")
  expect_error(n_twophase(transform(frame4, N = c(-1, 1, 1, 1)), phase1_cost = 1, budget = 1),
               "'N' must contain positive finite values")
  no_mean <- frame4[, setdiff(names(frame4), "mean")]
  expect_error(n_twophase(no_mean, phase1_cost = 1, cv = 0.05), "'mu' is required")
})

test_that("a plan supplies defaults", {
  p <- svyplan(N = 1e6)
  fit <- n_twophase(frame4, phase1_cost = 1, budget = 50000, plan = p)
  expect_equal(fit$params$N, 1e6)
})

test_that("the two design effects enter their own components", {
  # phase1_deff multiplies the between component only
  base <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  hi1 <- n_twophase(frame4, phase1_cost = 1, budget = 50000, phase1_deff = 4)
  # a larger phase-1 deff makes the between term dearer, so subsample harder
  expect_true(all(hi1$detail$nu <= base$detail$nu + 1e-9))
  expect_gt(hi1$cv, base$cv)

  # a per-stratum deff column raises that stratum's fraction
  d2 <- transform(frame4, deff = c(1, 3, 1, 1))
  hi2 <- n_twophase(d2, phase1_cost = 1, budget = 50000)
  expect_gt(hi2$detail$nu[2] / hi2$detail$nu[1],
            base$detail$nu[2] / base$detail$nu[1])
})

test_that("the closed form with deffs attains the constrained optimum", {
  fr <- transform(frame4, deff = c(1.2, 2.5, 0.8, 1.7))
  W <- fr$N / sum(fr$N)
  A <- sum(W * (fr$mean - sum(W * fr$mean))^2)
  d1 <- 2.2
  fit <- n_twophase(fr, phase1_cost = 1, budget = 50000, phase1_deff = d1)
  VC <- function(nu) (d1 * A + sum(fr$deff * W * fr$sd^2 / nu)) *
    (1 + sum(fr$unit_cost * W * nu))
  ref <- stats::optim(rep(0.5, 4), VC, method = "L-BFGS-B",
                      lower = 1e-9, upper = 1)$par
  expect_equal(fit$detail$nu, ref, tolerance = 1e-5)
})

test_that("single_deff moves only the comparator", {
  a <- n_twophase(frame4, phase1_cost = 1, budget = 50000, single_deff = 1)
  b <- n_twophase(frame4, phase1_cost = 1, budget = 50000, single_deff = 3)
  expect_equal(a$detail$nu, b$detail$nu)
  expect_equal(a$cv, b$cv)
  expect_gt(b$single_phase$cv, a$single_phase$cv)
  # a badly clustered comparator can flip the verdict
  expect_true(a$single_phase$better)
  expect_false(b$single_phase$better)
})

test_that("deffs default to 1 and reproduce the SRS result", {
  a <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  b <- n_twophase(transform(frame4, deff = 1), phase1_cost = 1,
                  budget = 50000, phase1_deff = 1, single_deff = 1)
  expect_equal(a$detail$nu, b$detail$nu)
  expect_equal(a$cv, b$cv)
})

test_that("design effects are validated", {
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1, phase1_deff = 0),
               "'phase1_deff' must be a positive finite scalar")
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1, single_deff = -1),
               "'single_deff' must be a positive finite scalar")
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1,
                          phase1_deff = c(1, 2)),
               "'phase1_deff' must be a positive finite scalar")
  expect_error(n_twophase(transform(frame4, deff = c(1, 0, 1, 1)),
                          phase1_cost = 1, budget = 1),
               "'deff' column must contain positive finite values")
})

test_that("response rates divide their own component", {
  base <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  # phase-1 response inflates the between component and costs precision
  p1 <- n_twophase(frame4, phase1_cost = 1, budget = 50000, resp_rate = 0.8)
  expect_gt(p1$cv, base$cv)
  # a per-stratum phase-2 rate raises that stratum's fraction, like a deff
  r2 <- transform(frame4, resp_rate = c(1, 0.5, 1, 1))
  hi <- n_twophase(r2, phase1_cost = 1, budget = 50000)
  expect_gt(hi$detail$nu[2] / hi$detail$nu[1],
            base$detail$nu[2] / base$detail$nu[1])
  # d/r is what enters, so a deff of 2 and a response of 0.5 coincide
  a <- n_twophase(transform(frame4, deff = c(1, 2, 1, 1)),
                  phase1_cost = 1, budget = 50000)
  b <- n_twophase(transform(frame4, resp_rate = c(1, 0.5, 1, 1)),
                  phase1_cost = 1, budget = 50000)
  expect_equal(a$detail$nu, b$detail$nu)
})

test_that("a common response rate cancels from the relative allocation", {
  base <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  same <- n_twophase(transform(frame4, resp_rate = 0.7),
                     phase1_cost = 1, budget = 50000)
  expect_equal(same$detail$nu / sum(same$detail$nu),
               base$detail$nu / sum(base$detail$nu), tolerance = 1e-8)
})

test_that("phase 2 cannot outdraw the pool phase 1 classified", {
  flat <- transform(frame4, mean = c(55, 55.2, 55, 55.2))
  for (r in c(0.4, 0.6, 0.9)) {
    fit <- n_twophase(flat, phase1_cost = 1, budget = 50000, resp_rate = r)
    expect_lte(max(fit$detail$nu), r + 1e-9)
  }
})

test_that("issued and responding are reported as separate quantities", {
  fr <- transform(frame4, resp_rate = c(0.9, 0.7, 0.85, 0.6))
  fit <- n_twophase(fr, phase1_cost = 1, budget = 50000, resp_rate = 0.8)
  expect_equal(fit$responding[["n_phase1"]], 0.8 * fit$n[["n_phase1"]])
  expect_equal(fit$detail$n_resp, fr$resp_rate * fit$detail$n_issued)
  expect_equal(fit$responding[["n_phase2"]], sum(fit$detail$n_resp))
  expect_true(all(fit$responding <= fit$n))
})

test_that("prec_twophase honours response and round trips", {
  fr <- transform(frame4, resp_rate = c(0.9, 0.7, 0.85, 0.6),
                  deff = c(1.2, 2.5, 0.8, 1.7))
  fit <- n_twophase(fr, phase1_cost = 1, budget = 50000,
                    resp_rate = 0.8, phase1_deff = 2.2)
  pr <- prec_twophase(fit)
  expect_equal(pr$cv, fit$cv, tolerance = 1e-8)
  back <- n_twophase(pr)
  expect_equal(back$cv, fit$cv, tolerance = 1e-6)
  # ignoring response would give a different answer
  no_r <- prec_twophase(
    transform(fit$detail[, c("stratum", "N", "sd", "unit_cost", "deff", "nu")],
              resp_rate = 1),
    n_phase1 = fit$n[["n_phase1"]], phase1_deff = 2.2
  )
  expect_false(isTRUE(all.equal(no_r$cv, fit$cv)))
})

test_that("the single-phase comparator carries its own response rate", {
  a <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  b <- n_twophase(frame4, phase1_cost = 1, budget = 50000,
                  single_resp_rate = 0.5)
  expect_equal(a$detail$nu, b$detail$nu)
  expect_gt(b$single_phase$cv, a$single_phase$cv)
})

test_that("response rates are validated", {
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1, resp_rate = 0),
               "resp_rate")
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1, resp_rate = 1.5),
               "resp_rate")
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1,
                          single_resp_rate = -1), "resp_rate")
  expect_error(n_twophase(transform(frame4, resp_rate = c(1, 0, 1, 1)),
                          phase1_cost = 1, budget = 1),
               "'resp_rate' column must contain values in \\(0, 1\\]")
})

test_that("nonresponse follow-up is unaffected by the defaults", {
  theta <- 0.5
  nrfu <- data.frame(
    stratum   = c("respondents", "nonrespondents"),
    N         = c(theta, 1 - theta),
    sd        = c(1, 1),
    unit_cost = c(0, 200),
    take_all  = c(TRUE, FALSE)
  )
  fit <- n_twophase(nrfu, phase1_cost = 50, budget = 100000, mu = 1)
  expect_equal(fit$detail$nu[2], sqrt(50 / (200 * theta)), tolerance = 1e-8)
  expect_equal(unname(fit$responding), unname(fit$n))
})

test_that("fixing n_phase1 keeps the shares and moves only the scale", {
  free <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  for (na in c(8000, 20000, 30000)) {
    fx <- n_twophase(frame4, phase1_cost = 1, budget = 50000, n_phase1 = na)
    expect_equal(fx$n[["n_phase1"]], na)
    if (!all(fx$detail$nu == 1)) {
      expect_equal(fx$detail$nu / sum(fx$detail$nu),
                   free$detail$nu / sum(free$detail$nu), tolerance = 1e-8)
    }
  }
})

test_that("the free optimum is the best member of the fixed family", {
  free <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  for (na in c(0.6, 0.85, 1.2, 1.6) * free$n[["n_phase1"]]) {
    fx <- n_twophase(frame4, phase1_cost = 1, budget = 50000,
                     n_phase1 = round(na))
    expect_gte(fx$cv, free$cv - 1e-9)
  }
})

test_that("fixed n_phase1 with a cv target meets it and reports the cost", {
  fx <- n_twophase(frame4, phase1_cost = 1, cv = 0.006, n_phase1 = 20000)
  expect_equal(fx$cv, 0.006, tolerance = 1e-8)
  expect_equal(fx$n[["n_phase1"]], 20000)
  expect_true(is.finite(fx$cost))
})

test_that("fixed n_phase1 has its own two feasibility errors", {
  expect_error(
    n_twophase(frame4, phase1_cost = 1, budget = 50000, n_phase1 = 60000),
    "budget does not reach phase 2"
  )
  expect_error(
    n_twophase(frame4, phase1_cost = 1, cv = 0.0001, n_phase1 = 8000),
    "target 'cv' is out of reach"
  )
  expect_error(
    n_twophase(frame4, phase1_cost = 1, budget = 1, n_phase1 = 1e9, N = 1e6),
    "exceeds the population size"
  )
})

test_that("the whole-unit design respects budget, pool and target", {
  b <- n_twophase(frame4, phase1_cost = 1, budget = 50000)
  expect_true(all(b$operational$n_int == floor(b$operational$n_int)))
  expect_lte(b$operational$cost, 50000)
  expect_lte(b$operational$n[["n_phase2"]], b$operational$n[["n_phase1"]])
  # cv mode rounds up, so the integer design must not miss the target
  cvfit <- n_twophase(frame4, phase1_cost = 1, cv = 0.006)
  expect_lte(cvfit$operational$cv, 0.006 + 1e-9)
  # the pool bounds every stratum
  r <- n_twophase(frame4, phase1_cost = 1, budget = 50000, resp_rate = 0.6)
  expect_true(all(r$operational$n_int <=
                    floor(0.6 * r$detail$share * r$operational$n[["n_phase1"]])))
})

test_that("the assurance issue is the smallest that clears the level", {
  # The grid reaches below a level of a half and up to a response rate of
  # 0.99, the two places where the expected count is not the right place to
  # start looking: below a half the answer sits under it, and at a high rate
  # the expected count itself overshoots.
  for (r in c(0.05, 0.5, 0.7, 0.9, 0.99, 1)) {
    for (m in c(1, 2, 20, 100, 500)) {
      for (lvl in c(0.01, 0.1, 0.49, 0.5, 0.51, 0.8, 0.95, 0.999)) {
        g <- svyplan:::.assure_size(m, r, lvl)
        expect_gte(stats::pbinom(m - 1, g, r, lower.tail = FALSE), lvl)
        expect_lt(stats::pbinom(m - 1, g - 1, r, lower.tail = FALSE), lvl)
        # Never fewer issued than the respondents required.
        expect_gte(g, m)
      }
    }
  }
  expect_equal(svyplan:::.assure_size(0, 0.5, 0.9), 0)
  # Two cases the walk from the expected count got wrong.
  expect_equal(svyplan:::.assure_size(100, 0.64, 0.1), 144)
  expect_equal(svyplan:::.assure_size(2, 0.99, 0.8), 2)
  # A rate shorter than the counts is recycled rather than read as NA.
  expect_equal(
    svyplan:::.assure_size(c(10, 20), 0.5, 0.9),
    c(svyplan:::.assure_size(10, 0.5, 0.9), svyplan:::.assure_size(20, 0.5, 0.9))
  )
})

test_that("assurance inflates issue over the expectation and costs more", {
  fr <- transform(frame4, resp_rate = c(0.6, 0.8, 0.7, 0.9))
  fit <- n_twophase(fr, phase1_cost = 1, budget = 50000,
                    resp_rate = 0.75, assurance = 0.9)
  o <- fit$operational
  expect_true(all(o$assured >= o$n_int))
  expect_gt(o$assured_cost, o$cost)
  # and it really does clear the level for the respondents planned on
  need <- ceiling(o$n_int * fr$resp_rate)
  expect_true(all(stats::pbinom(need - 1, o$assured, fr$resp_rate,
                                lower.tail = FALSE) >= 0.9))
  # absent by default
  expect_null(n_twophase(fr, phase1_cost = 1, budget = 50000)$operational$assured)
})

test_that("assurance is validated", {
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1e5, assurance = 0),
               "'assurance' must be a probability")
  expect_error(n_twophase(frame4, phase1_cost = 1, budget = 1e5, assurance = 1),
               "'assurance' must be a probability")
})

test_that("the phase-1 correction applies component by component", {
  fr <- data.frame(
    stratum   = c("A", "B"),
    N         = c(600, 400),
    sd        = c(10, 20),
    mean      = c(50, 70),
    unit_cost = 1,
    nu        = c(0.3, 0.5)
  )
  W <- fr$N / sum(fr$N)
  A <- sum(W * (fr$mean - sum(W * fr$mean))^2)
  mu <- sum(W * fr$mean)
  for (N in c(1000, 2000, 10000)) {
    for (n1 in c(200, 500, 900)) {
      ref <- A * (1 / n1 - 1 / N) +
        sum(W * fr$sd^2 * (1 / (fr$nu * n1) - 1 / N))
      got <- prec_twophase(fr, n_phase1 = n1, N = N)
      expect_equal(got$se, sqrt(ref), tolerance = 1e-10)
      expect_equal(got$cv, sqrt(ref) / mu, tolerance = 1e-10)
    }
  }
})

test_that("a phase-1 census keeps the variance phase 2 still has to measure", {
  # phase 1 fixes the stratum weights, phase 2 still sees half of each stratum
  fr <- data.frame(
    stratum = c("A", "B"), N = c(60, 40), sd = c(10, 20),
    mean = c(50, 70), unit_cost = 1, nu = c(0.5, 0.5)
  )
  W <- fr$N / sum(fr$N)
  n2 <- fr$nu * fr$N
  ref <- sum(W^2 * (1 - n2 / fr$N) * fr$sd^2 / n2)
  got <- prec_twophase(fr, n_phase1 = 100, N = 100)
  expect_equal(got$se, sqrt(ref), tolerance = 1e-10)
  expect_gt(got$cv, 0)
  # and a census that also measures everything is exactly zero-variance
  full <- transform(fr, nu = c(1, 1))
  expect_equal(prec_twophase(full, n_phase1 = 100, N = 100)$se, 0)
})

test_that("fixed n_phase1 in budget mode pays for the take-all strata first", {
  fr <- data.frame(
    stratum   = c("pin", "free"),
    N         = c(50, 50),
    sd        = c(1, 1),
    mean      = c(1, 2),
    unit_cost = c(100, 1),
    take_all  = c(TRUE, FALSE)
  )
  expect_error(
    n_twophase(fr, phase1_cost = 1, budget = 200, n_phase1 = 100),
    "measuring the take-all strata"
  )
  # the budget that clears the mandatory cost does return a usable design
  fit <- n_twophase(fr, phase1_cost = 1, budget = 6000, n_phase1 = 100)
  expect_true(all(fit$detail$nu > 0))
  expect_true(all(fit$detail$n_issued > 0))
  expect_true(all(fit$operational$n_int >= 0))
  expect_gt(fit$n[["n_phase2"]], 0)
})

test_that("fixed n_phase1 in cv mode counts the residual of pinned strata", {
  fr <- data.frame(
    stratum   = c("pin", "free"),
    N         = c(50, 50),
    sd        = c(100, 1),
    mean      = c(1, 1),
    unit_cost = 1,
    take_all  = c(TRUE, FALSE)
  )
  # cv 0.1 is unreachable: the pinned stratum alone leaves cv near 7
  expect_error(
    n_twophase(fr, phase1_cost = 1, cv = 0.1, n_phase1 = 100, mu = 1),
    "target 'cv' is out of reach"
  )
  # a reachable target is still met
  fit <- n_twophase(fr, phase1_cost = 1, cv = 8, n_phase1 = 100, mu = 1)
  expect_lte(fit$cv, 8 * (1 + 1e-6))
})

test_that("integerization does not strand a stratum it can afford", {
  fr <- data.frame(
    stratum   = c("bulk", "mid", "thin"),
    N         = c(600, 380, 20),
    sd        = c(30, 15, 0.6),
    mean      = c(50, 52, 60),
    unit_cost = c(1, 1, 3)
  )
  # the buy-back scores a first unit against an infinite term, not a zero gain
  fit <- n_twophase(fr, phase1_cost = 1, budget = 204)
  expect_true(all(fit$operational$n_int > 0))
  expect_true(is.finite(fit$operational$cv))
  expect_lte(fit$operational$cost, 204)
})

test_that("a stratum the phase-1 pool cannot reach is reported, not hidden", {
  fr <- data.frame(
    stratum   = c("large", "tiny"),
    N         = c(999, 1),
    sd        = c(10, 10),
    mean      = c(50, 80),
    unit_cost = 1
  )
  expect_warning(
    fit <- n_twophase(fr, phase1_cost = 1, budget = 100),
    "rounds to zero units in 'tiny'"
  )
  expect_equal(fit$operational$n_int[2], 0)
  expect_true(is.infinite(fit$operational$cv))
})

test_that("a stratum with no variance needs no phase-2 sample", {
  fr <- data.frame(
    stratum   = c("a", "b"),
    N         = c(999, 1),
    sd        = c(10, 0),
    mean      = c(50, 80),
    unit_cost = 1
  )
  fit <- n_twophase(fr, phase1_cost = 1, budget = 100)
  expect_equal(fit$operational$n_int[2], 0)
  expect_true(is.finite(fit$operational$cv))
  expect_true(is.finite(fit$cv))
})

test_that("take_all rejects numbers that are not 0 or 1", {
  legacy <- data.frame(stratum = c("A", "B"), N = c(100, 100), sd = c(1, 2),
                       take_all = c(2, 0))
  expect_error(n_alloc(legacy, n = 150, alloc = "neyman"),
               "must be logical")
  two_phase <- data.frame(stratum = c("A", "B"), N = c(100, 100),
                          sd = c(1, 2), mean = c(10, 10), unit_cost = 1,
                          take_all = c(2, 0))
  expect_error(n_twophase(two_phase, phase1_cost = 1, budget = 300),
               "must be logical")
  # 0/1 and logical both still work
  ok <- transform(legacy, take_all = c(1, 0))
  expect_true(n_alloc(ok, n = 150, alloc = "neyman")$detail$take_all[1])
  ok2 <- transform(legacy, take_all = c(TRUE, FALSE))
  expect_true(n_alloc(ok2, n = 150, alloc = "neyman")$detail$take_all[1])
})

test_that("the precision inverse rebuilds the same problem", {
  frame <- data.frame(
    stratum = c("A", "B"), N = c(600, 400), sd = c(10, 20),
    mean = c(50, 70), unit_cost = c(2, 3)
  )

  # a fixed phase-1 size survives the round trip
  fixed <- n_twophase(frame, phase1_cost = 1, cv = 0.08, n_phase1 = 500,
                      mu = 58)
  back <- n_twophase(prec_twophase(fixed))
  expect_equal(back$n, fixed$n)
  expect_equal(back$cv, fixed$cv)
  expect_equal(back$cost, fixed$cost)
  expect_equal(back$detail$nu, fixed$detail$nu)
  expect_equal(back$params$n_phase1_fixed, 500)

  # so do the take-all set, the assurance level and the comparator
  pinned <- transform(frame, take_all = c(TRUE, FALSE))
  full <- suppressWarnings(
    n_twophase(pinned, phase1_cost = 1, budget = 20000, assurance = 0.9,
               single_deff = 2.5, resp_rate = 0.9)
  )
  inverse <- suppressWarnings(n_twophase(prec_twophase(full)))
  expect_equal(inverse$params$take_all, c(TRUE, FALSE))
  expect_equal(inverse$params$assurance, 0.9)
  expect_equal(inverse$params$single_deff, 2.5)
  expect_equal(inverse$detail$nu, full$detail$nu)
  expect_equal(inverse$n, full$n)

  # and the caller can still ask for the other mode
  by_budget <- suppressWarnings(
    n_twophase(prec_twophase(full), budget = 20000)
  )
  expect_equal(by_budget$n, full$n)
})

test_that("a fully pinned frame is priced but not allocated", {
  frame <- data.frame(
    stratum = c("A", "B"), N = c(600, 400), sd = c(10, 20),
    mean = c(50, 70), unit_cost = c(2, 3), take_all = c(TRUE, TRUE)
  )
  expect_error(n_twophase(frame, phase1_cost = 1, budget = 1e4),
               "nothing is left to subsample")
  # evaluating one is well defined
  priced <- prec_twophase(transform(frame, nu = c(1, 1)), n_phase1 = 200)
  expect_true(is.finite(priced$cv))
  expect_gt(priced$cv, 0)
})

test_that("assurance reports a phase-2 issue the phase-1 pool cannot supply", {
  frame <- data.frame(
    stratum = LETTERS[1:3], N = c(100, 100, 100), sd = c(10, 12, 14),
    mean = c(50, 55, 60), unit_cost = 2, resp_rate = 0.8
  )
  expect_warning(
    fit <- n_twophase(frame, phase1_cost = 1, budget = 500, resp_rate = 0.9,
                      assurance = 0.95),
    "exceeds the pool phase 1 supplies"
  )
  expect_true(any(fit$operational$assured > fit$operational$pool))

  # a design that leaves headroom in the pool does not warn
  costly <- transform(frame, unit_cost = 50)
  expect_silent(
    ok <- n_twophase(costly, phase1_cost = 1, budget = 5000, resp_rate = 0.9,
                     assurance = 0.95)
  )
  expect_true(all(ok$operational$assured <= ok$operational$pool))
})

test_that("svyplan can carry the two-phase screening cost", {
  plan <- svyplan(phase1_cost = 1, N = 1e5)
  frame <- data.frame(
    stratum = c("A", "B"), N = c(600, 400), sd = c(10, 20),
    mean = c(50, 70), unit_cost = c(2, 3)
  )
  fit <- n_twophase(frame, budget = 5000, plan = plan)
  expect_equal(fit$params$phase1_cost, 1)
  expect_error(svyplan(phase1_cost = -1), "phase1_cost")
})

test_that("no predict method is promised for the two-phase class", {
  expect_false(exists("predict.svyplan_twophase", mode = "function"))
})

test_that("the census term carries design effects but not response", {
  # issuing phase 2 to everyone still leaves the variance nonresponse costs
  fr <- data.frame(
    stratum = c("A", "B"), N = c(60, 40), sd = c(10, 20),
    mean = c(50, 70), unit_cost = 1, nu = 1, resp_rate = 0.8
  )
  W <- fr$N / sum(fr$N)
  mu <- sum(W * fr$mean)
  reference <- sum(W * fr$sd^2 * (1 / (fr$resp_rate * fr$nu) - 1)) / sum(fr$N)

  got <- prec_twophase(fr, n_phase1 = 100, N = 100)
  expect_equal(got$se, sqrt(reference), tolerance = 1e-10)
  expect_equal(got$cv, sqrt(reference) / mu, tolerance = 1e-10)

  # complete response and complete measurement is still exactly zero
  expect_equal(prec_twophase(transform(fr, resp_rate = 1),
                             n_phase1 = 100, N = 100)$se, 0)

  # a phase-1 census that classified only 80% of units cannot carry more than
  # 80% into phase 2, and what it misses stays in the variance
  partial <- transform(fr[, setdiff(names(fr), "resp_rate")], nu = 0.8)
  total <- sum(W * (fr$mean - mu)^2) + sum(W * fr$sd^2)
  expect_equal(
    prec_twophase(partial, n_phase1 = 100, N = 100, resp_rate = 0.8)$se,
    sqrt(total * (1 / 0.8 - 1) / 100), tolerance = 1e-10
  )

  # a target under the response-limited floor is refused
  alloc <- fr[, setdiff(names(fr), "nu")]
  expect_error(
    n_twophase(alloc, phase1_cost = 1, cv = 0.01, n_phase1 = 100, N = 100),
    "out of reach"
  )
})

test_that("the single-phase comparator uses the same census correction", {
  fr <- data.frame(
    stratum = c("A", "B"), N = c(60, 40), sd = c(10, 20),
    mean = c(50, 70), unit_cost = 1
  )
  W <- fr$N / sum(fr$N)
  mu <- sum(W * fr$mean)
  total <- sum(W * (fr$mean - mu)^2) + sum(W * fr$sd^2)

  # issued to the whole frame at 80% response: not a census of the outcome
  single <- n_twophase(fr, phase1_cost = 1, budget = 1e6, N = 100,
                       single_resp_rate = 0.8)$single_phase
  expect_equal(single$n, 100)
  expect_equal(single$cv, sqrt(total * (1 / 0.8 - 1) / 100) / mu,
               tolerance = 1e-10)

  # a comparator that cannot reach the target is capped at N, reports the
  # cv a census of the frame would give, and does not win on cost alone
  hard <- n_twophase(fr, phase1_cost = 5, cv = 0.012, N = 1000,
                     single_resp_rate = 0.3, single_cost = 0.001)
  expect_equal(hard$single_phase$n, 1000)
  expect_false(hard$single_phase$reaches_target)
  expect_gt(hard$single_phase$cv, 0.012)
  expect_lt(hard$single_phase$cost, hard$cost)
  expect_false(hard$single_phase$better)
})

test_that("the corrected comparator can change the design verdict", {
  fr <- data.frame(
    stratum = c("A", "B"), N = c(600, 400), sd = c(10, 20),
    mean = c(50, 52), unit_cost = 1
  )
  fit <- n_twophase(fr, phase1_cost = 0.5, budget = 600, N = 1000,
                    single_resp_rate = 0.5)
  # reusing var_unit as the census term reported cv 0.01068 here and handed
  # the verdict to the single-phase design; the correct figure is above the
  # two-phase design's own
  expect_gt(fit$single_phase$cv, fit$cv)
  expect_false(fit$single_phase$better)
})

test_that("the two-phase variance is non-negative before it is clamped", {
  grid <- expand.grid(
    N = c(Inf, 400, 1000), n1 = c(50, 200, 400),
    r1 = c(1, 0.7), r2 = c(1, 0.6), nu = c(0.3, 1), d1 = c(1, 2), d2 = c(1, 1.5)
  )
  grid <- grid[grid$n1 <= grid$N & grid$nu <= grid$r1, ]
  expect_gt(nrow(grid), 20)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    fr <- data.frame(
      stratum = c("A", "B"), N = c(60, 40), sd = c(10, 20),
      mean = c(50, 70), unit_cost = 1, nu = g$nu, deff = g$d2,
      resp_rate = g$r2
    )
    got <- prec_twophase(fr, n_phase1 = g$n1, N = g$N, resp_rate = g$r1,
                         phase1_deff = g$d1)
    W <- fr$N / sum(fr$N)
    A <- sum(W * (fr$mean - sum(W * fr$mean))^2)
    var_unit <- g$d1 * A / g$r1 + sum(g$d2 / g$r2 * W * fr$sd^2 / g$nu)
    var_pop <- g$d1 * A + sum(g$d2 * W * fr$sd^2)
    raw <- var_unit / g$n1 - if (is.infinite(g$N)) 0 else var_pop / g$N
    expect_gte(raw, -1e-9)
    expect_equal(got$se, sqrt(max(raw, 0)), tolerance = 1e-10)
  }
})

test_that("prec_twophase refuses a phase-2 issue the pool cannot supply", {
  fr <- data.frame(
    stratum = c("A", "B"), N = c(60, 40), sd = c(10, 20),
    mean = c(50, 70), unit_cost = 1, nu = c(0.9, 0.9)
  )
  expect_error(
    prec_twophase(fr, n_phase1 = 100, N = 1000, resp_rate = 0.8),
    "cannot exceed the classified fraction"
  )
  # exactly at the classification rate is allowed, as is the default
  expect_s3_class(
    prec_twophase(transform(fr, nu = c(0.8, 0.8)), n_phase1 = 100, N = 1000,
                  resp_rate = 0.8),
    "svyplan_prec"
  )
  expect_s3_class(prec_twophase(fr, n_phase1 = 100, N = 1000), "svyplan_prec")
})

test_that("take_all pins at the classification rate, not at one", {
  fr <- data.frame(
    stratum = c("pin", "free"), N = c(60, 40), sd = c(10, 1),
    mean = c(50, 70), unit_cost = c(1, 100), take_all = c(TRUE, FALSE)
  )
  fit <- n_twophase(fr, phase1_cost = 1, budget = 1000, resp_rate = 0.8)
  expect_equal(fit$detail$nu[1], 0.8)
  # with full classification the two readings coincide
  full <- n_twophase(fr, phase1_cost = 1, budget = 1000)
  expect_equal(full$detail$nu[1], 1)
})

test_that("two-phase stratum labels identify one stratum each", {
  fr <- data.frame(
    stratum = c("A", "A"), N = c(60, 40), sd = c(10, 20),
    mean = c(50, 70), unit_cost = 1
  )
  expect_error(n_twophase(fr, phase1_cost = 1, budget = 1000), "must be unique")
  expect_error(
    n_twophase(transform(fr, stratum = c("A", NA)), phase1_cost = 1,
               budget = 1000),
    "missing or empty"
  )
  expect_error(
    n_twophase(transform(fr, stratum = c("A", " ")), phase1_cost = 1,
               budget = 1000),
    "missing or empty"
  )
  # omitted labels fall back to unambiguous row numbers
  bare <- n_twophase(fr[, setdiff(names(fr), "stratum")], phase1_cost = 1,
                     budget = 1000)
  expect_equal(bare$detail$stratum, c("1", "2"))
})

test_that("stratum means that cancel do not pass for a usable 'mu'", {
  f <- data.frame(
    stratum = c("A", "B"), N = c(5000, 5000), sd = c(2, 2),
    mean = c(0.1, -0.1), unit_cost = c(2, 2)
  )
  expect_error(
    n_twophase(f, phase1_cost = 1, cv = 0.05),
    "'mu' is required in cv mode"
  )
  # an explicit 'mu' is the user's own number and is used as given
  expect_s3_class(n_twophase(f, phase1_cost = 1, cv = 0.05, mu = 3),
                  "svyplan_twophase")
})

test_that("a two-phase assurance level prints as itself", {
  fr <- data.frame(N = c(1000, 2000), sd = c(10, 15), mean = c(5, 8),
                   unit_cost = c(2, 3))
  out <- capture.output(print(
    n_twophase(fr, phase1_cost = 1, cv = 0.05, assurance = 0.999)
  ))
  expect_true(any(grepl("assured (0.999)", out, fixed = TRUE)))
  expect_false(any(grepl("assured (1.00)", out, fixed = TRUE)))
})

## TP-print. The printed block is the design as it would be fielded

.tp_frame <- function() {
  data.frame(
    stratum   = c("A", "B", "C", "D"),
    N         = c(3500, 2500, 2500, 1500),
    sd        = c(12, 25, 8, 40),
    mean      = c(40, 70, 35, 90),
    unit_cost = c(2, 5, 1, 9)
  )
}

test_that("print carries one reading of the design, not two", {
  plan <- n_twophase(.tp_frame(), phase1_cost = 1, budget = 50000)
  out <- capture.output(print(plan))
  expect_length(out, 10L)
  expect_lt(max(nchar(out)), 80L)
  # The continuous optimum and the fielded design differ by a unit or two and
  # report the same cv, so only the fielded one is printed.
  expect_false(any(grepl("^issued:", out)))
  expect_false(any(grepl("n_issued", out, fixed = TRUE)))
  expect_false(any(grepl("^---$", out)))
  expect_match(out[2L], "^field design: n_phase1 = [0-9]+ \\| n_phase2 = [0-9]+$")
  # The comparator is a verdict, so it is one line.
  expect_length(grep("single-phase", out), 1L)
  expect_identical(out[length(out)],
                   "# summary() for the continuous optimum and the comparator")
})

test_that("every printed count is on the fielded path", {
  frame <- transform(.tp_frame(), resp_rate = c(0.8, 0.75, 0.9, 0.7))
  plan <- n_twophase(frame, phase1_cost = 1, budget = 50000)
  shown <- svyplan:::.twophase_shown(plan)
  expect_identical(shown$n, plan$operational$n)
  expect_identical(shown$n_int, plan$operational$n_int)
  # The header sizes, the stratum takes and the responding counts must all
  # reconcile: a responding figure read off the continuous design is a number
  # nothing else in the block adds up to.
  expect_identical(sum(shown$n_int), unname(shown$n[["n_phase2"]]))
  brief <- svyplan:::.fmt_twophase_detail(plan, brief = TRUE)
  expect_identical(brief$n_int, format(shown$n_int))
  expect_identical(brief$n_resp,
                   format(round(shown$n_int * plan$detail$resp_rate)))
  resp_line <- svyplan:::.fmt_twophase_responding(plan)
  expect_match(resp_line, sprintf(
    "n_phase2 = %s", format(round(sum(shown$n_int * plan$detail$resp_rate)))
  ), fixed = TRUE)
  # The stored fields stay continuous, which is what keeps the round trip
  # exact.
  expect_false(isTRUE(all.equal(plan$detail$n_issued, plan$operational$n_int)))
})

test_that("the conditional blocks all survive the cut", {
  frame <- transform(.tp_frame(), resp_rate = c(0.8, 0.75, 0.9, 0.7),
                     deff = c(1.2, 1.1, 1.3, 1.5))
  out <- capture.output(print(
    n_twophase(frame, phase1_cost = 1, budget = 50000, assurance = 0.9,
               fixed_cost = 5000)
  ))
  expect_true(any(grepl("^expected responding:", out)))
  expect_true(any(grepl("(fixed: 5000)", out, fixed = TRUE)))
  expect_true(any(grepl("^assured \\(0.90\\):", out)))
  expect_match(out[grep("^ stratum", out)], "deff")
  expect_match(out[grep("^ stratum", out)], "resp")
  expect_match(out[grep("^ stratum", out)], "n_resp")
  expect_lt(max(nchar(out)), 80L)
})

test_that("summary carries the continuous optimum and the comparator", {
  plan <- n_twophase(.tp_frame(), phase1_cost = 1, budget = 50000)
  sm <- summary(plan)
  expect_s3_class(sm, "summary.svyplan_twophase")
  expect_identical(sm$plan, plan)
  expect_identical(sm$n, plan$n)
  expect_identical(sm$single_phase, plan$single_phase)
  # The full table is the brief one plus the continuous issue.
  brief <- svyplan:::.fmt_twophase_detail(plan, brief = TRUE)
  expect_identical(sm$detail[names(brief)], brief)
  expect_identical(sm$detail$n_issued, format(round(plan$detail$n_issued)))
  out <- capture.output(print(sm))
  expect_true(any(grepl("^continuous optimum: ", out)))
  expect_true(any(grepl("^Single-phase comparator: ", out)))
  expect_true(any(grepl("reaches the target", out, fixed = TRUE)))
})
