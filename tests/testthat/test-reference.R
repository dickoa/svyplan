# Comparative tests against published results
# Cochran (1977) Chapter 3: Proportions
test_that("n_prop matches Cochran (1977) SRS formula", {
  # Cochran Eq 3.22: n0 = z^2 * p*q / e^2
  # with FPC: n = n0 / (1 + n0/N)
  p <- 0.5
  e <- 0.05
  z <- qnorm(0.975)
  n0 <- z^2 * p * (1 - p) / e^2
  result <- n_prop(p = p, moe = e)
  expect_equal(result$n, n0, tolerance = 1e-6)

  # with finite population
  N <- 1000
  n_fpc <- n0 / (1 + (n0 - 1) / N)
  result_fpc <- n_prop(p = p, moe = e, N = N)
  expect_equal(result_fpc$n, n_fpc, tolerance = 1e-2)
})

test_that("n_prop p=0.5 gives largest sample (conservative)", {
  # Cochran (1977): p=0.5 maximizes p*q, so gives largest n
  n_05 <- n_prop(p = 0.5, moe = 0.05)$n
  n_03 <- n_prop(p = 0.3, moe = 0.05)$n
  n_01 <- n_prop(p = 0.1, moe = 0.05)$n
  n_09 <- n_prop(p = 0.9, moe = 0.05)$n
  expect_true(n_05 >= n_03)
  expect_true(n_05 >= n_01)
  expect_true(n_05 >= n_09)
})

# Cochran (1977) Chapter 5: Means
test_that("n_mean matches Cochran (1977) SRS formula", {
  # Cochran Eq 5.3: n0 = z^2 * S^2 / e^2
  # with FPC: n = n0 / (1 + n0/N)
  S2 <- 100
  e <- 2
  z <- qnorm(0.975)
  n0 <- z^2 * S2 / e^2
  result <- n_mean(var = S2, moe = e)
  expect_equal(result$n, n0, tolerance = 1e-6)

  # with finite population
  N <- 5000
  n_fpc <- n0 / (1 + n0 / N)
  result_fpc <- n_mean(var = S2, moe = e, N = N)
  expect_equal(result_fpc$n, n_fpc, tolerance = 1e-2)
})

test_that("n_mean CV mode matches Cochran formula", {
  # CV = SE/mu, so SE = cv*mu, and n = S^2/(mu^2 * cv^2) for SRS
  S2 <- 100
  mu <- 50
  cv <- 0.05
  CVpop <- sqrt(S2) / mu
  n0 <- CVpop^2 / cv^2
  result <- n_mean(var = S2, mu = mu, cv = cv)
  expect_equal(result$n, n0, tolerance = 1e-6)
})

# VDK (2018) Chapter 3: Design Effect
test_that("design_effect weighting matches VDK Eq 3.5", {
  # VDK Eq 3.5: deff_w = n * sum(w^2) / sum(w)^2
  # equivalently: 1 + CV_pop(w)^2
  w <- c(1, 1, 1, 1, 5)
  n <- length(w)
  deff_exp <- n * sum(w^2) / sum(w)^2
  result <- design_effect(weights = w)
  expect_equal(as.double(result), deff_exp, tolerance = 1e-6)
})

test_that("design_effect weighting = 1 for equal weights (VDK)", {
  # VDK: equal weights => deff_w = 1
  w <- rep(3, 100)
  result <- design_effect(weights = w)
  expect_equal(as.double(result), 1, tolerance = 1e-10)
})

test_that("design_effect clustering matches VDK Eq 3.18", {
  # VDK Eq 3.18: deff_c = 1 + (b_bar - 1) * icc
  # where b_bar = avg cluster size, icc = ICC
  icc <- 0.05
  b_bar <- 20
  deff_exp <- 1 + (b_bar - 1) * icc
  result <- design_effect(icc = icc, n_per_psu = b_bar)
  expect_equal(as.double(result), deff_exp, tolerance = 1e-10)
})

# VDK (2018) Chapter 4: Power Analysis
test_that("power_prop matches VDK Eq 4.7 (Wald, SRS)", {
  # VDK Eq 4.7: n = (z_a + z_b)^2 * (p1*q1 + p2*q2) / icc^2
  p1 <- 0.30
  p2 <- 0.35
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  icc <- abs(p1 - p2)
  V <- p1 * (1 - p1) + p2 * (1 - p2)
  n_exp <- (z_a + z_b)^2 * V / icc^2
  result <- power_prop(p1 = p1, p2 = p2)
  expect_equal(result$n, n_exp, tolerance = 1e-4)
})

test_that("power_mean matches VDK Eq 4.14 (SRS)", {
  # VDK Eq 4.14: n = (z_a + z_b)^2 * 2 * sigma^2 / icc^2
  sigma2 <- 100
  icc <- 5
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  n_exp <- (z_a + z_b)^2 * 2 * sigma2 / icc^2
  result <- power_mean(effect = icc, var = sigma2)
  expect_equal(result$n, n_exp, tolerance = 1e-4)
})

test_that("n_prop MICS-style with resp_rate and deff", {
  # Full MICS formula: n = z^2 * p*(1-p) * deff / (RME*p)^2 / resp_rate
  # RME = moe / p
  p <- 0.15
  RME <- 0.15
  deff <- 2.0
  rr <- 0.90
  moe <- RME * p
  z <- qnorm(0.975)
  n0 <- z^2 * p * (1 - p) / moe^2
  n_exp <- n0 * deff / rr
  result <- n_prop(p = p, moe = moe, deff = deff, resp_rate = rr)
  expect_equal(result$n, n_exp, tolerance = 1e-4)
})

# Cluster design: optimal allocation (VDK Ch 9)
test_that("n_cluster optimal allocation matches VDK formula", {
  # VDK Eq 9.14: m_opt = sqrt(c1/c2 * (1-icc)/icc)
  c1 <- 500
  c2 <- 50
  icc <- 0.05
  m_opt <- sqrt(c1 / c2 * (1 - icc) / icc)
  result <- n_cluster(stage_cost = c(c1, c2), icc = icc, cv = 0.05)
  expect_equal(unname(result$n[2]), max(2, ceiling(m_opt)), tolerance = 1)
})

# Precision round-trips preserve VDK relationships
test_that("prec_prop SE matches VDK formula", {
  # SE(p_hat) = sqrt(p*q / n_eff), where n_eff = n * resp_rate / deff
  p <- 0.3
  n <- 400
  deff <- 1.5
  rr <- 0.9
  n_eff <- n * rr / deff
  se_exp <- sqrt(p * (1 - p) / n_eff)
  result <- prec_prop(p = p, n = n, deff = deff, resp_rate = rr)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

test_that("prec_mean SE with FPC matches Cochran formula", {
  # SE = sqrt(S^2 * (1 - n/N) / n)
  S2 <- 100
  n <- 200
  N <- 1000
  fpc <- 1 - n / N
  se_exp <- sqrt(S2 * fpc / n)
  result <- prec_mean(var = S2, n = n, N = N)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

# VDK (2018) Example 5.2: minimum relvariance of estimated total revenue under
# a fixed budget, book pp. 133-140 with the comparison in Table 5.4 p. 162
.vdk_example_5_2 <- function() {
  sector <- c("Manufacturing", "Retail", "Wholesale", "Service", "Finance")
  N <- c(6221, 11738, 4333, 22809, 5467)
  # VDK compute proportion SDs with the finite-population factor
  sd_prop <- function(p) sqrt(p * (1 - p) * N / (N - 1))
  list(
    frame = data.frame(
      stratum = sector, N = N, unit_cost = c(120, 80, 80, 90, 150),
      stringsAsFactors = FALSE
    ),
    measures = data.frame(
      stratum = rep(sector, 4),
      name = rep(c("revenue", "employees", "research", "offshore"), each = 5),
      mean = c(
        85, 11, 23, 17, 126,
        511, 21, 70, 32, 157,
        0.8, 0.2, 0.5, 0.3, 0.9,
        0.06, 0.03, 0.03, 0.21, 0.77
      ),
      sd = c(
        170.0, 8.8, 23.0, 25.5, 315.0,
        255.50, 5.25, 35.00, 32.00, 471.00,
        sd_prop(c(0.8, 0.2, 0.5, 0.3, 0.9)),
        sd_prop(c(0.06, 0.03, 0.03, 0.21, 0.77))
      ),
      stringsAsFactors = FALSE
    ),
    targets = data.frame(
      name = c("employees", "research", "offshore"),
      domain = ".overall", level = NA_character_,
      cv = c(0.05, 0.03, 0.03),
      stringsAsFactors = FALSE
    )
  )
}

test_that("n_alloc reproduces VDK Example 5.2 under a fixed budget", {
  z <- .vdk_example_5_2()
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = "revenue", budget = 300000, min_n_stratum = 100
  )

  # Table 5.4: proc nlp and nloptr both land here
  expect_equal(round(fit$detail$n), c(413, 318, 124, 1397, 596))
  expect_equal(sum(fit$detail$n), 2847.56, tolerance = 1e-4)
  expect_equal(fit$params$achieved$cost, 300000, tolerance = 1e-6)

  # proc nlp reports 0.0021705237; svyplan agrees to nine significant figures
  expect_equal(fit$objective_value, 0.002170523648, tolerance = 1e-9)
  expect_equal(sqrt(fit$objective_value), 0.04658888, tolerance = 1e-7)

  achieved <- setNames(fit$constraints$.achieved, fit$constraints$name)
  expect_equal(unname(achieved["employees"]), 0.0239, tolerance = 1e-3)
  expect_equal(unname(achieved["research"]), 0.0208, tolerance = 1e-3)
  # the offshore constraint is the binding one, at exactly its 3% limit
  expect_equal(unname(achieved["offshore"]), 0.03, tolerance = 1e-9)
  expect_identical(fit$binding, "offshore@.overall:cv")
})

test_that("the VDK Example 5.2 rerun at 350000 improves on the published CV", {
  z <- .vdk_example_5_2()
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    objective = "revenue", budget = 350000, min_n_stratum = 100
  )
  # VDK report CV(revenue) = 0.0409 at this budget (p. 140). Their Solver run
  # starts from the 300000 solution and stops at a local point; the
  # epsilon-constraint route finds a strictly feasible design at 0.0388.
  expect_equal(sqrt(fit$objective_value), 0.03875365, tolerance = 1e-6)
  expect_lt(sqrt(fit$objective_value), 0.0409)
  expect_equal(fit$params$achieved$cost, 350000, tolerance = 1e-6)
  expect_true(all(fit$constraints$.pass))
  expect_lte(fit$operational$cost, 350000)
})

test_that("VDK Example 5.2 at 250000 cannot fund its targets", {
  z <- .vdk_example_5_2()
  expect_error(
    n_alloc(
      z$frame, measures = z$measures, targets = z$targets,
      objective = "revenue", budget = 250000, min_n_stratum = 100
    ),
    "cheapest target-feasible design costs 265192, short by 15191"
  )
})

## surveyplanning (Breidaks, Liberts, Jukams; CSB Latvia) v4.0
# Values produced by their functions and pasted in as literals, so the
# comparison does not depend on that package being installed.

test_that("prec_alloc stratum precision matches surveyplanning expvar", {
  # expvar(Yh = "Yh", H = "H", s2h = "s2h", nh = "nh", poph = "poph",
  #        Dom = "dom", dataset = data.table(
  #          H = 1:4, Yh = c(4000, 9000, 12000, 900), s2h = c(9, 16, 25, 4),
  #          poph = c(800, 1200, 1500, 300), nh = c(100, 150, 200, 60),
  #          dom = c("N", "N", "S", "S")))
  f <- data.frame(
    N = c(800, 1200, 1500, 300),
    sd = sqrt(c(9, 16, 25, 4)),
    mean = c(4000, 9000, 12000, 900) / c(800, 1200, 1500, 300),
    dom = c("N", "N", "S", "S")
  )
  p <- prec_alloc(f, n = c(100, 150, 200, 60), domains = "dom")

  # their var and se are on the stratum-total scale, ours on the mean scale
  expect_equal(p$detail$.se * f$N,
               c(224.49944, 366.60606, 493.71044, 69.28203),
               tolerance = 1e-7)
  expect_equal((p$detail$.se * f$N)^2,
               c(50400, 134400, 243750, 4800), tolerance = 1e-7)
  # cv is scale free, so it compares directly (theirs is a percentage)
  expect_equal(p$detail$.cv * 100,
               c(5.612486, 4.073401, 4.114254, 7.698004), tolerance = 1e-6)
  expect_equal(p$domains$.cv * 100, c(3.306798, 3.864712), tolerance = 1e-6)
  expect_equal(p$cv * 100, 2.541673, tolerance = 1e-6)
})

test_that("n_alloc reproduces surveyplanning optsize, including take-all", {
  # optsize(H = "H", n = 600, poph = "poph", s2h = "s2h", dataset = ...)
  f <- data.frame(N = c(800, 1200, 1500, 300), sd = sqrt(c(9, 16, 25, 4)))
  expect_equal(n_alloc(f, n = 600, alloc = "neyman")$detail$n,
               c(94.11765, 188.23529, 294.11765, 23.52941), tolerance = 1e-6)

  # with fullsampleh = c(0, 0, 0, 1)
  expect_equal(
    n_alloc(transform(f, take_all = c(FALSE, FALSE, FALSE, TRUE)),
            n = 600, alloc = "neyman")$detail$n,
    c(48.97959, 97.95918, 153.06122, 300), tolerance = 1e-6
  )
})

test_that("the solved proportion matches surveyplanning min_prop", {
  # min_prop(n = 15e3, pop = 2e6, RMoE = 0.1, R = 0.75, deff_sam = 1.4)
  # their RMoE is moe/p, so the equivalent cv target is RMoE / qnorm(0.975)
  r <- prec_prop(n = 15e3, cv = 0.1 / qnorm(0.975), N = 2e6,
                 resp_rate = 0.75, deff = 1.4)
  expect_equal(r$params$p, 0.0453788177, tolerance = 1e-9)
  # min_count is pop * min_prop
  expect_equal(r$params$p * 2e6, 90757.6353, tolerance = 1e-6)
})

test_that("n_alloc reproduces optsize with per-stratum Rh and deffh", {
  # optsize(H = "H", n = 600, poph = "poph", s2h = "s2h", Rh = "Rh",
  #         deffh = "deffh", dataset = data.table(
  #           H = 1:4, s2h = c(9, 16, 25, 4), poph = c(800, 1200, 1500, 300),
  #           Rh = c(.9, .8, .85, 1), deffh = c(1.5, 1.2, 2, 1)))
  f <- data.frame(N = c(800, 1200, 1500, 300), sd = sqrt(c(9, 16, 25, 4)))
  d <- c(1.5, 1.2, 2, 1)
  r <- c(0.9, 0.8, 0.85, 1)

  expect_equal(n_alloc(f, n = 600, alloc = "neyman", deff = d, resp_rate = r)$detail$n,
               c(88.18253, 167.31458, 327.42642, 17.07647), tolerance = 1e-7)

  # with fullsampleh = c(0, 0, 0, 1)
  expect_equal(
    n_alloc(transform(f, take_all = c(FALSE, FALSE, FALSE, TRUE)),
            n = 600, alloc = "neyman", deff = d, resp_rate = r)$detail$n,
    c(45.38290, 86.10799, 168.50911, 300), tolerance = 1e-7
  )
})
