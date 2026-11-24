# Analytical and finite-population Monte Carlo verification for arbitrary
# S-stage designs (S = 2, 3, 4, 5) across varcomp(), n_cluster(),
# prec_cluster(), design_effect(), effective_n(), n_alloc(), and n_alloc_joint().

test_that("4-stage and 5-stage n_cluster / prec_cluster round-trip exactly", {
  # 4-stage CV mode with stage FPC, vector resp_rate, and custom stage names
  sc4 <- c(province = 1500, district = 400, ea = 100, hh = 20)
  icc4 <- c(0.02, 0.05, 0.10)
  rr4 <- c(0.95, 0.92, 0.90, 0.85)
  N4 <- c(40, 240, 1200, 12000)

  plan4 <- n_cluster(
    stage_cost = sc4, icc = icc4, unit_relvar = 1.6,
    cv = 0.05, resp_rate = rr4, N = N4
  )
  expect_s3_class(plan4, "svyplan_cluster")
  expect_identical(plan4$stages, 4L)
  expect_identical(names(plan4$n), names(sc4))
  expect_identical(names(plan4$operational$n), names(sc4))
  expect_equal(plan4$cv, 0.05, tolerance = 1e-6)
  expect_lte(plan4$operational$cv, 0.05 + 1e-8)

  p4 <- prec_cluster(plan4)
  expect_equal(p4$cv, plan4$cv, tolerance = 1e-7)
  rt4 <- n_cluster(p4)
  expect_equal(unname(rt4$n), unname(plan4$n), tolerance = 1e-4)

  # 4-stage partial fixed stage_n
  fix4 <- n_cluster(
    stage_cost = sc4, icc = icc4, unit_relvar = 1.6,
    cv = 0.05, stage_n = c(NA, 3, NA, 8)
  )
  expect_equal(unname(fix4$n[c(2, 4)]), c(3, 8))
  expect_equal(prec_cluster(fix4)$cv, 0.05, tolerance = 1e-6)

  # 5-stage budget mode
  sc5 <- c(1500, 600, 250, 80, 15)
  icc5 <- c(0.01, 0.03, 0.05, 0.10)
  plan5 <- n_cluster(
    stage_cost = sc5, icc = icc5, unit_relvar = 1.2, budget = 120000
  )
  expect_identical(plan5$stages, 5L)
  expect_identical(names(plan5$n), paste0("stage_", 1:5))
  expect_equal(plan5$cost, 120000, tolerance = 1e-4)
  expect_lte(plan5$operational$cost, 120000 + 1e-8)
  expect_equal(prec_cluster(plan5)$cv, plan5$cv, tolerance = 1e-7)
})

test_that("full census gives zero CV and design_effect identity holds for S = 4, 5", {
  # 4-stage census: M = c(10, 5, 4, 6) => N = c(10, 50, 200, 1200)
  pc4 <- prec_cluster(
    n = c(10, 5, 4, 6), icc = c(0.02, 0.05, 0.10),
    unit_relvar = 1.5, N = c(10, 50, 200, 1200)
  )
  expect_equal(pc4$cv, 0, tolerance = 1e-12)

  # 5-stage census: M = c(6, 4, 3, 3, 5) => N = c(6, 24, 72, 216, 1080)
  pc5 <- prec_cluster(
    n = c(6, 4, 3, 3, 5), icc = c(0.01, 0.03, 0.06, 0.10),
    unit_relvar = 1.2, N = c(6, 24, 72, 216, 1080)
  )
  expect_equal(pc5$cv, 0, tolerance = 1e-12)

  # design_effect() at unit takes equals var_ratio[1] for S = 4, 5
  for (vr1 in c(1, 1.35)) {
    expect_equal(
      as.double(design_effect(icc = c(0.02, 0.05, 0.10), take = c(1, 1, 1),
                              var_ratio = vr1)),
      vr1,
      tolerance = 1e-12
    )
    expect_equal(
      as.double(design_effect(icc = c(0.01, 0.03, 0.06, 0.10),
                              take = c(1, 1, 1, 1), var_ratio = vr1)),
      vr1,
      tolerance = 1e-12
    )
  }

  # effective_n() == 1 / cv^2 (with unit_relvar = 1, no FPC) for S = 4, 5
  p4_nf <- prec_cluster(
    n = c(20, 4, 3, 8), icc = c(0.02, 0.05, 0.10),
    unit_relvar = 1, resp_rate = c(0.95, 0.9, 0.85, 0.8)
  )
  expect_equal(effective_n(p4_nf), 1 / p4_nf$cv^2, tolerance = 1e-10)

  p5_nf <- prec_cluster(
    n = c(15, 3, 2, 3, 5), icc = c(0.01, 0.03, 0.05, 0.10),
    unit_relvar = 1, resp_rate = c(0.95, 0.95, 0.9, 0.9, 0.8)
  )
  expect_equal(effective_n(p5_nf), 1 / p5_nf$cv^2, tolerance = 1e-10)
})

test_that("C branch-and-bound finds exact global integer optimum on 4-stage problems", {
  # Exhaustive brute-force check on 4-stage box (1..15 x 1..5 x 1..4 x 1..15)
  sc <- c(1000, 250, 80, 15)
  icc <- c(0.03, 0.06, 0.12)
  vr <- 1.5
  N_vec <- c(15, 75, 300, 4500) # M = c(15, 5, 4, 15)

  grid <- expand.grid(n1 = 1:15, n2 = 1:5, n3 = 1:4, n4 = 1:15)
  M <- c(15, 5, 4, 15)
  k <- c(1, 1 - icc[1], (1 - icc[1]) * (1 - icc[2]))
  R <- vr * c(k[1] * icc[1], k[2] * icc[2], k[3] * icc[3], k[3] * (1 - icc[3]))
  cv2_all <- R[1] * (1 - grid$n1 / M[1]) / grid$n1 +
    R[2] * (1 - grid$n2 / M[2]) / (grid$n1 * grid$n2) +
    R[3] * (1 - grid$n3 / M[3]) / (grid$n1 * grid$n2 * grid$n3) +
    R[4] * (1 - grid$n4 / M[4]) / (grid$n1 * grid$n2 * grid$n3 * grid$n4)
  cost_all <- sc[1] * grid$n1 + sc[2] * grid$n1 * grid$n2 +
    sc[3] * grid$n1 * grid$n2 * grid$n3 +
    sc[4] * grid$n1 * grid$n2 * grid$n3 * grid$n4

  # CV mode
  target_cv <- 0.08
  fit_cv <- n_cluster(
    stage_cost = sc, icc = icc, unit_relvar = vr, cv = target_cv, N = N_vec
  )
  feas_cv <- which(sqrt(pmax(cv2_all, 0)) <= target_cv + 1e-9)
  best_cost <- min(cost_all[feas_cv])
  expect_equal(fit_cv$operational$cost, best_cost, tolerance = 1e-8)
  expect_lte(fit_cv$operational$cv, target_cv + 1e-8)

  # Budget mode
  budget <- 16000
  fit_b <- n_cluster(
    stage_cost = sc, icc = icc, unit_relvar = vr, budget = budget, N = N_vec
  )
  feas_b <- which(cost_all <= budget + 1e-9)
  best_cv <- sqrt(min(pmax(cv2_all[feas_b], 0)))
  expect_equal(fit_b$operational$cv, best_cv, tolerance = 1e-8)
  expect_lte(fit_b$operational$cost, budget + 1e-8)
})

test_that("4-stage domain, classic and joint designs reassess at their whole sizes", {
  ind4 <- data.frame(
    region = c("N", "N", "S", "S"),
    name = c("y1", "y2", "y1", "y2"),
    p = c(0.3, 0.5, 0.25, 0.45),
    icc_1 = c(0.02, 0.03, 0.03, 0.04),
    icc_2 = c(0.05, 0.04, 0.06, 0.05),
    icc_3 = c(0.08, 0.10, 0.09, 0.11),
    cv = c(0.06, 0.05, 0.06, 0.05)
  )
  sc <- c(1000, 300, 80, 15)
  cl4 <- n_cluster(ind4, stage_cost = sc, domains = "region")
  op <- cl4$operational$domains
  for (lev in c("N", "S")) {
    rows <- ind4[ind4$region == lev, setdiff(names(ind4), c("region", "cv"))]
    n_d <- unlist(op[op$region == lev, paste0("stage_", 1:4)])
    rows$n_1 <- n_d[1]
    rows$take_2 <- n_d[2]
    rows$take_3 <- n_d[3]
    rows$take_4 <- n_d[4]
    assessed <- prec_cluster(rows, stage_cost = sc)
    expect_true(all(assessed$detail$.cv <= ind4$cv[ind4$region == lev] + 1e-12))
    expect_equal(op$.cost[op$region == lev], sum(sc * cumprod(n_d)))
  }

  cl4_b <- n_cluster(ind4, stage_cost = sc, domains = "region",
                     allocation = "joint", budget = 150000)
  expect_equal(cl4_b$operational$cost, sum(cl4_b$operational$domains$.cost))
  expect_lte(cl4_b$operational$cost, 150000)

  fr4 <- data.frame(
    stratum = c("Urban", "Rural"), N = c(120000, 180000), sd = c(12, 16),
    mean = c(50, 60), icc_1 = c(0.02, 0.03), icc_2 = c(0.05, 0.06),
    icc_3 = c(0.10, 0.12), cost_1 = c(1000, 1200), cost_2 = c(300, 350),
    cost_3 = c(80, 90), cost_4 = c(15, 20), N_1 = c(30, 40),
    N_2 = c(300, 400), N_3 = c(3000, 4000)
  )
  a4 <- n_alloc(fr4, cv = 0.025, alloc = "optimal", fpc = "stage")
  d <- a4$detail
  whole <- transform(fr4, take_2 = d$take_2_int, take_3 = d$take_3_int,
                     take_4 = d$take_4_int)
  expect_equal(prec_alloc(whole, n = d$n_int, fpc = "stage")$cv,
               a4$operational$cv, tolerance = 1e-12)
  cum <- cbind(1, d$take_2_int, d$take_2_int * d$take_3_int,
               d$take_2_int * d$take_3_int * d$take_4_int)
  costs <- as.matrix(fr4[paste0("cost_", 1:4)])
  expect_equal(a4$operational$cost, sum(d$n_1_int * rowSums(costs * cum)))
  expect_lte(a4$operational$cv, 0.025)

  jf4 <- data.frame(
    stratum = c("Urban", "Rural"), N = c(120000, 180000), N_1 = c(40, 50),
    N_2 = c(400, 500), N_3 = c(4000, 5000), take_2 = c(3, 3),
    take_3 = c(2, 2), take_4 = c(5, 5), cost_1 = c(900, 1100),
    cost_2 = c(220, 260), cost_3 = c(60, 70), cost_4 = c(15, 18)
  )
  jm4 <- data.frame(
    stratum = rep(c("Urban", "Rural"), 2),
    name = rep(c("prev", "income"), each = 2),
    p = c(0.35, 0.45, NA, NA), mean = c(NA, NA, 500, 420),
    sd = c(NA, NA, 180, 160), icc_1 = c(0.02, 0.03, 0.03, 0.04),
    icc_2 = c(0.05, 0.06, 0.04, 0.05), icc_3 = c(0.09, 0.11, 0.08, 0.10)
  )
  jt4 <- data.frame(name = c("prev", "income"), cv = c(0.05, 0.04))
  j4 <- n_alloc_joint(jf4, measures = jm4, targets = jt4, fpc = "stage")
  field <- prec_alloc_joint(j4, n = j4$detail$n_int)
  expect_equal(field$detail$.achieved, j4$operational$constraints$.achieved,
               tolerance = 1e-12)
  pred4 <- predict(j4, data.frame(take_4 = c(4, 5, 6)))
  expect_equal(pred4$cost_int[pred4$take_4 == 5], j4$operational$cost)
})

srs_stage_var <- function(y, M, n) {
  S <- length(M)
  arr <- array(y, dim = rev(M))
  V <- 0
  for (s in seq_len(S)) {
    m <- apply(arr, rev(seq_len(S))[seq_len(s)], mean)
    S2 <- if (s == 1L) var(as.vector(m)) else mean(apply(m, seq_len(s - 1L), var))
    V <- V + (1 - n[s] / M[s]) * S2 / prod(n[seq_len(s)])
  }
  V
}

nested_population <- function(M, sds, mu) {
  S <- length(M)
  y <- rep(mu, prod(M))
  ids <- list()
  for (s in seq_len(S - 1L)) {
    id <- rep(seq_len(prod(M[seq_len(s)])), each = prod(M[(s + 1L):S]))
    ids[[s]] <- id
    y <- y + rnorm(prod(M[seq_len(s)]), 0, sds[s])[id]
  }
  y <- y + rnorm(prod(M), 0, sds[S])
  c(list(y = y), stats::setNames(ids, paste0("s", seq_len(S - 1L))))
}

test_that("4-stage SRSWOR precision equals the stage-wise variance of the population", {
  set.seed(202604)
  M <- c(12L, 6L, 5L, 8L)
  pop <- as.data.frame(nested_population(M, c(5, 6, 8, 12), 100))
  vc <- varcomp(y ~ s1 / s2 / s3, data = pop)
  n_draw <- c(5L, 3L, 2L, 4L)
  pred <- prec_cluster(n = n_draw, icc = vc, N = cumprod(M))
  closed <- sqrt(srs_stage_var(pop$y, M, n_draw)) / mean(pop$y)
  expect_equal(pred$cv, closed, tolerance = 1e-12)

  y_arr <- array(pop$y, dim = rev(M))
  reps <- 2500L
  est <- numeric(reps)
  for (r in seq_len(reps)) {
    tot <- 0
    for (a in sample.int(M[1], n_draw[1])) {
      for (b in sample.int(M[2], n_draw[2])) {
        for (cc in sample.int(M[3], n_draw[3])) {
          tot <- tot + sum(y_arr[sample.int(M[4], n_draw[4]), cc, b, a])
        }
      }
    }
    est[r] <- tot / prod(n_draw)
  }
  ## the relative standard error of a standard deviation from reps draws
  expect_equal(sd(est) / mean(pop$y), closed,
               tolerance = 4 / sqrt(2 * (reps - 1)))
})

test_that("5-stage SRSWOR precision equals the stage-wise variance of the population", {
  set.seed(202606)
  M <- c(8L, 4L, 4L, 3L, 6L)
  pop <- as.data.frame(nested_population(M, c(5, 5, 6, 7, 10), 120))
  vc <- varcomp(y ~ s1 / s2 / s3 / s4, data = pop)
  n_draw <- c(4L, 2L, 2L, 2L, 3L)
  pred <- prec_cluster(n = n_draw, icc = vc, N = cumprod(M))
  expect_equal(pred$cv, sqrt(srs_stage_var(pop$y, M, n_draw)) / mean(pop$y),
               tolerance = 1e-12)
})

test_that("stratified 4-stage SRSWOR precision equals the stage-wise variance", {
  set.seed(202607)
  MA <- c(10L, 5L, 4L, 6L)
  MB <- c(12L, 4L, 4L, 6L)
  dA <- as.data.frame(nested_population(MA, c(5, 6, 7, 10), 100))
  dB <- as.data.frame(nested_population(MB, c(7, 5, 8, 12), 140))
  pop <- rbind(cbind(stratum = "A", dA), cbind(stratum = "B", dB))
  vc <- varcomp(y ~ s1 / s2 / s3, data = pop, strata = ~stratum)
  st <- vc$strata
  n1_h <- c(4L, 5L)
  takes <- c(2L, 2L, 3L)
  frame <- data.frame(
    stratum = c("A", "B"), N = c(nrow(dA), nrow(dB)),
    N_1 = c(MA[1], MB[1]), N_2 = c(prod(MA[1:2]), prod(MB[1:2])),
    N_3 = c(prod(MA[1:3]), prod(MB[1:3])),
    take_2 = takes[1], take_3 = takes[2], take_4 = takes[3],
    cost_1 = c(500, 600), cost_2 = c(150, 180), cost_3 = c(50, 60),
    cost_4 = c(10, 12)
  )
  measures <- data.frame(
    stratum = st$stratum, name = "y", mean = st$mean, sd = st$sd,
    icc_1 = st$icc_1, icc_2 = st$icc_2, icc_3 = st$icc_3,
    var_ratio_1 = st$var_ratio_1, var_ratio_2 = st$var_ratio_2,
    var_ratio_3 = st$var_ratio_3
  )
  prec <- prec_alloc_joint(
    frame, n = n1_h * prod(takes), measures = measures,
    targets = data.frame(name = "y", cv = 0.05), fpc = "stage"
  )
  W <- c(nrow(dA), nrow(dB)) / nrow(pop)
  V <- W[1]^2 * srs_stage_var(dA$y, MA, c(n1_h[1], takes)) +
    W[2]^2 * srs_stage_var(dB$y, MB, c(n1_h[2], takes))
  expect_equal(prec$detail$.achieved[1L], sqrt(V) / mean(pop$y),
               tolerance = 1e-12)
})
