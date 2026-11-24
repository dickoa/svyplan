# Comparative tests against published results
## Cochran (1977, 3rd ed.), Chapter 4: the estimation of sample size
test_that("n_prop matches Cochran (1977) SRS formula", {
  # Cochran (1977), Sect. 4.4, pp. 75 to 76, Eqs. (4.1) to (4.3):
  # n0 = z^2 * p*q / e^2, and with the fpc n = n0 / (1 + (n0 - 1)/N), (4.3)
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

test_that("n_mean matches Cochran (1977) SRS formula", {
  # Cochran (1977), Sect. 4.6, p. 78: for an absolute error d,
  # n0 = t^2 * S^2 / d^2, and with the fpc n = n0 / (1 + n0/N) as in (4.3)
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
  # Cochran (1977), Eq. (4.5): n0 = (S/Ybar)^2 / C, C the desired cv^2
  S2 <- 100
  mu <- 50
  cv <- 0.05
  CVpop <- sqrt(S2) / mu
  n0 <- CVpop^2 / cv^2
  result <- n_mean(var = S2, mu = mu, cv = cv)
  expect_equal(result$n, n0, tolerance = 1e-6)
})

## Valliant, Dever and Kreuter (2018, 2nd ed.)
test_that("design_effect weighting matches VDK Eq 14.8", {
  # VDK (2018), Sect. 14.4.1, p. 396, Eq. (14.8): Kish's
  # deff_K = 1 + relvar(w), which equals n * sum(w^2) / sum(w)^2
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

test_that("design_effect clustering matches VDK Eq 9.5", {
  # VDK (2018), Sect. 9.2.1, p. 213, Eqs. (9.5) and (9.6) at k = 1:
  # deff_c = 1 + (b_bar - 1) * delta, b_bar the elements per PSU
  icc <- 0.05
  b_bar <- 20
  deff_exp <- 1 + (b_bar - 1) * icc
  result <- design_effect(icc = icc, take = b_bar)
  expect_equal(as.double(result), deff_exp, tolerance = 1e-10)
})

test_that("power_prop matches VDK Eq 4.18 (Wald, SRS)", {
  # VDK (2018), Sect. 4.3.3, Eq. (4.18), with the pooled null term Delta1
  # replaced by the unpooled Delta2, as the Wald method does, and the
  # quantile at 1 - alpha/2 for a two-sided test:
  # n = (z_a + z_b)^2 * (p1*q1 + p2*q2) / delta^2
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

test_that("power_mean matches VDK Eq 4.15 (SRS)", {
  # VDK (2018), Sect. 4.3.2, Eq. (4.15) at r = 1, which is (4.14) with no
  # overlap, with the quantile at 1 - alpha/2 for a two-sided test:
  # n = (z_a + z_b)^2 * 2 * sigma^2 / delta^2
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

test_that("n_cluster optimal allocation matches VDK formula", {
  # VDK (2018), Sect. 9.3.1, p. 232, Eq. (9.24):
  # m_opt = sqrt(c1/c2 * (1 - delta)/delta)
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
  fit <- n_alloc_joint(
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
  fit <- n_alloc_joint(
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
    n_alloc_joint(
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

## Further worked examples, at their printed precision

test_that("VDK Example 3.9: stratified allocation of the SMHO frame", {
  smho <- data.frame(
    stratum = as.character(1:6), N = c(215, 65, 252, 50, 149, 144),
    mean = c(21240408, 10024876, 4913008, 11927573, 6118415, 15567731),
    sd = c(26787207, 10645109, 6909676, 11085034, 9817762, 44553355),
    unit_cost = c(1400, 200, 300, 600, 450, 1000)
  )
  share <- function(fit) round(fit$detail$n / sum(fit$detail$n), 3)
  expect_equal(share(n_alloc(smho, n = 100, alloc = "neyman")),
               c(0.346, 0.042, 0.105, 0.033, 0.088, 0.386))
  b100 <- n_alloc(smho, budget = 100000, alloc = "optimal")
  expect_equal(share(b100), c(0.257, 0.082, 0.168, 0.038, 0.115, 0.339))
  expect_equal(round(b100$n), 119)
  expect_equal(round(n_alloc(smho, budget = 200000, alloc = "optimal")$n), 238)
  expect_equal(round(n_alloc(smho, cv = 0.05, alloc = "optimal")$n), 406)
  expect_equal(round(n_alloc(smho, cv = 0.10, alloc = "optimal")$n), 198)
})

test_that("VDK Example 9.7: two-stage optimum under a budget", {
  fits <- lapply(c(0.01, 0.05, 0.10, 0.20), function(d) {
    n_cluster(stage_cost = c(750, 100), icc = d, unit_relvar = 1,
              budget = 100000)
  })
  expect_equal(vapply(fits, function(f) round(f$n[[1]], 1), 1),
               c(28.8, 51.4, 63.6, 77.1))
  expect_equal(vapply(fits, function(f) round(f$n[[2]], 1), 1),
               c(27.2, 11.9, 8.2, 5.5))
  expect_equal(vapply(fits, function(f) round(f$cv, 4), 1),
               c(0.0401, 0.0502, 0.0574, 0.0670))
})

test_that("VDK Example 9.8: three-stage optimum under a budget", {
  one <- n_cluster(stage_cost = c(500, 100, 120), icc = c(0.01, 0.10),
                   unit_relvar = 1, var_ratio = c(1, 1), budget = 100000)
  expect_equal(round(unname(one$n), 1), c(28.3, 7.1, 2.7))
  expect_equal(round(one$cv, 4), 0.0499)
  fits <- lapply(c(0.01, 0.05, 0.10), function(d1) {
    n_cluster(stage_cost = c(500, 100, 120), icc = c(d1, 0.10),
              unit_relvar = 2, var_ratio = c(1, 1), budget = 100000)
  })
  sizes <- t(vapply(fits, function(f) round(unname(f$n), 1), numeric(3)))
  expect_equal(sizes[, 1], c(28.3, 53.9, 68.6))
  expect_equal(sizes[, 2], c(7.1, 3.2, 2.2))
  expect_equal(sizes[, 3], c(2.7, 2.7, 2.7))
  expect_equal(vapply(fits, function(f) round(f$cv, 4), 1),
               c(0.0706, 0.0830, 0.0922))
})

test_that("VDK Example 17.13: double sampling for stratification", {
  frame <- data.frame(stratum = as.character(1:4), N = rep(0.25, 4),
                      p = c(0.02, 0.12, 0.37, 0.54), unit_cost = 50)
  fit <- n_twophase(frame, phase1_cost = 10, budget = 20000)
  expect_equal(round((fit$cv / fit$single_phase$cv)^2, 2), 1.06)
})

test_that("VDK Example 17.14: nonresponse follow-up under a budget", {
  frame <- data.frame(stratum = c("resp", "nonresp"), N = c(0.5, 0.5),
                      mean = 1, sd = 1, unit_cost = c(0, 200),
                      take_all = c(TRUE, FALSE))
  fit <- n_twophase(frame, phase1_cost = 50, budget = 100000,
                    single_cost = 50, single_resp_rate = 0.5)
  expect_equal(round(fit$n[["n_phase1"]]), 828)
  expect_equal(round(fit$detail$n_issued[2]), 293)
  expect_equal(round(fit$detail$nu[2], 4), 0.7071)
  expect_equal(round(fit$cv, 4), 0.0382)
})

## Hankin, Mohr and Newman (2019), Examples 5.7, 7.1, 9.1.3, 10.3, Table 13.1

test_that("Satterthwaite degrees of freedom of six stratified allocations", {
  frame <- data.frame(stratum = c("1", "2", "3"), N = c(100, 200, 50),
                      var = c(200, 500, 4000))
  allocations <- list(c(20, 18, 2), c(25, 10, 5), c(20, 10, 10),
                      c(10, 10, 20), c(5, 5, 30), c(11, 23, 6))
  df <- vapply(allocations, function(a) {
    floor(as.double(design_df(prec_alloc(frame, n = a),
                              method = "satterthwaite")))
  }, 1)
  expect_equal(df, c(1, 11, 16, 13, 5, 12))
  expect_equal(n_alloc(frame, n = 40, alloc = "proportional")$detail$n_int,
               c(11, 23, 6))
  ney <- n_alloc(frame, n = 40, alloc = "neyman")
  expect_equal(ney$detail$n_int, c(6, 20, 14))
  expect_equal(floor(as.double(design_df(ney, method = "satterthwaite"))), 36)
})

test_that("the delta-method variance of a ratio from a five-unit population", {
  x <- c(10, 20, 30, 40, 50)
  y <- c(2, 5, 7, 14, 20)
  r <- mean(y) / mean(x)
  args <- list(r = r, cv_num = sd(y) / mean(y), cv_den = sd(x) / mean(x),
               component_cor = cor(x, y), N = 5)
  p <- do.call(prec_ratio, c(args, n = 3))
  expect_equal(round((mean(x) * p$se)^2, 2), 0.92)
  expect_equal(do.call(n_ratio, c(args, moe = qnorm(0.975) * p$se))$n, 3,
               tolerance = 1e-9)
})

test_that("the variance of a two-stage mean from a nine-unit population", {
  dat <- data.frame(psu = rep(1:3, each = 3),
                    y = c(5, 10, 12, 3, 15, 6, 13, 6, 11))
  p <- prec_cluster(c(2, 2), icc = varcomp(y ~ psu, data = dat),
                    N = c(3, 9), fpc = "stage")
  expect_equal(round((p$cv * mean(dat$y))^2, 3), 1.972)
})

test_that("the efficiency of two-phase sampling for stratification", {
  frame <- data.frame(stratum = c("A", "B"), N = 5000, sd = sqrt(5),
                      mean = c(100 - sqrt(95), 100 + sqrt(95)),
                      unit_cost = 10)
  fit <- n_twophase(frame, phase1_cost = 1, budget = 10000)
  expect_equal(fit$detail$nu, rep(sqrt(0.1 * 0.05 / 0.95), 2),
               tolerance = 1e-6)
  expect_equal(round((fit$single_phase[["cv"]] / fit$cv)^2, 1), 3.5)
})

test_that("the variance of a change under retention, replacement and partial overlap", {
  N <- 400
  n <- 40
  v <- c(100, 164)
  covar <- 100
  r <- covar / sqrt(prod(v))
  change_var <- function(...) {
    (N * prec_change(var = v, n = n, N = N, ...)$se)^2
  }
  expect_equal(change_var(overlap = 1, overlap_cor = r),
               N^2 * (N - n) / N * (sum(v) - 2 * covar) / n)
  expect_equal(change_var(overlap = 0), N^2 * (N - n) / N * sum(v) / n)
  m <- n / 2
  partial <- 0.25 * N^2 * (N - m) / N * (sum(v) - 2 * covar) / m +
    0.25 * N^2 * (N - m) / N * sum(v) / m
  expect_equal(change_var(overlap = 0.5, overlap_cor = r,
                          overlap_draw = "independent"), partial)
})
