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
test_that("design_effect kish matches VDK Eq 3.5", {
  # VDK Eq 3.5: deff_w = n * sum(w^2) / sum(w)^2
  # equivalently: 1 + CV_pop(w)^2
  w <- c(1, 1, 1, 1, 5)
  n <- length(w)
  deff_exp <- n * sum(w^2) / sum(w)^2
  result <- design_effect(w, method = "kish")
  expect_equal(as.double(result), deff_exp, tolerance = 1e-6)
})

test_that("design_effect kish = 1 for equal weights (VDK)", {
  # VDK: equal weights => deff_w = 1
  w <- rep(3, 100)
  result <- design_effect(w, method = "kish")
  expect_equal(as.double(result), 1, tolerance = 1e-10)
})

test_that("design_effect cluster matches VDK Eq 3.18", {
  # VDK Eq 3.18: deff_c = 1 + (b_bar - 1) * delta
  # where b_bar = avg cluster size, delta = ICC
  delta <- 0.05
  b_bar <- 20
  deff_exp <- 1 + (b_bar - 1) * delta
  result <- design_effect(delta = delta, psu_size = b_bar, method = "cluster")
  expect_equal(as.double(result), deff_exp, tolerance = 1e-10)
})

# VDK (2018) Chapter 4: Power Analysis
test_that("power_prop matches VDK Eq 4.7 (Wald, SRS)", {
  # VDK Eq 4.7: n = (z_a + z_b)^2 * (p1*q1 + p2*q2) / delta^2
  p1 <- 0.30
  p2 <- 0.35
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  delta <- abs(p1 - p2)
  V <- p1 * (1 - p1) + p2 * (1 - p2)
  n_exp <- (z_a + z_b)^2 * V / delta^2
  result <- power_prop(p1 = p1, p2 = p2)
  expect_equal(result$n, n_exp, tolerance = 1e-4)
})

test_that("power_mean matches VDK Eq 4.14 (SRS)", {
  # VDK Eq 4.14: n = (z_a + z_b)^2 * 2 * sigma^2 / delta^2
  sigma2 <- 100
  delta <- 5
  z_a <- qnorm(0.975)
  z_b <- qnorm(0.80)
  n_exp <- (z_a + z_b)^2 * 2 * sigma2 / delta^2
  result <- power_mean(effect = delta, var = sigma2)
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
  # VDK Eq 9.14: m_opt = sqrt(c1/c2 * (1-delta)/delta)
  c1 <- 500
  c2 <- 50
  delta <- 0.05
  m_opt <- sqrt(c1 / c2 * (1 - delta) / delta)
  result <- n_cluster(stage_cost = c(c1, c2), delta = delta, cv = 0.05)
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
