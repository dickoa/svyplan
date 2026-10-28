## Small finite populations, reproduced exactly
##
## Each case is a population small enough to enumerate or to work by hand,
## so the stage-by-stage correction can be checked against the variance of
## the actual sampling design rather than against another formula.

## A two-stage population of three PSUs of three units

.fp_y <- matrix(c(5, 10, 12, 3, 6, 15, 13, 6, 11), nrow = 3, byrow = TRUE)

test_that("the stage correction equals the enumerated two-stage variance", {
  A <- 3
  B <- 3
  a <- 2
  b <- 2
  mu <- mean(.fp_y)
  est <- c()
  for (ps in combn(A, a, simplify = FALSE)) {
    for (u in combn(B, b, simplify = FALSE)) {
      for (v in combn(B, b, simplify = FALSE)) {
        est <- c(est, mean(c(.fp_y[ps[1], u], .fp_y[ps[2], v])))
      }
    }
  }
  enumerated <- mean((est - mu)^2)
  vc <- varcomp(as.vector(t(.fp_y)), stage_id = list(rep(1:3, each = 3)))
  by_cluster <- prec_cluster(n = c(a, b), icc = vc, fpc = "stage",
                             N = 9, N_psu = 3)
  expect_equal((by_cluster$cv * mu)^2, enumerated, tolerance = 1e-12)
  frame <- data.frame(
    N = 9, N_psu = 3, n_per_psu = b, sd = sd(as.vector(.fp_y)), mean = mu,
    icc_psu = vc$icc, var_ratio_psu = vc$var_ratio
  )
  by_alloc <- prec_alloc(frame, n = a * b, fpc = "stage")
  expect_equal(by_alloc$se^2, enumerated, tolerance = 1e-12)
})

## A three-stage population of three PSUs of three SSUs of three units
##
## The between-PSU component varcomp() returns here has the divisor A - 1,
## as at two stages, which is what makes the PSU term exact.

.fp_y3 <- array(
  c(5, 9, 2, 7, 4, 8, 6, 3, 10, 12, 8, 14, 9, 11, 7, 13, 10, 6,
    3, 5, 4, 8, 2, 6, 7, 9, 1),
  dim = c(3, 3, 3)
)

test_that("the stage correction equals the enumerated three-stage variance", {
  pairs <- combn(3, 2, simplify = FALSE)
  psu_means <- lapply(1:3, function(i) {
    out <- c()
    for (s in pairs) for (u in pairs) for (v in pairs) {
      out <- c(out, mean(c(.fp_y3[i, s[1], u], .fp_y3[i, s[2], v])))
    }
    out
  })
  est <- c()
  for (ps in pairs) {
    est <- c(est, as.vector(outer(psu_means[[ps[1]]], psu_means[[ps[2]]],
                                  "+")) / 2)
  }
  mu <- mean(.fp_y3)
  enumerated <- mean((est - mu)^2)
  psu <- rep(1:3, 9)
  ssu <- paste(psu, rep(rep(1:3, each = 3), 3))
  vc <- varcomp(as.vector(.fp_y3), stage_id = list(psu, ssu))
  by_cluster <- prec_cluster(n = c(2, 2, 2), icc = vc, fpc = "stage",
                             N = 27, N_psu = 3, N_ssu = 9)
  expect_equal((by_cluster$cv * mu)^2, enumerated, tolerance = 1e-12)
})

## A forest of 100 plots of 49 trees, between and within variances 5 and 49

test_that("the cost-optimal take keeps the within-PSU correction", {
  fit <- n_cluster(stage_cost = c(12, 3), icc = 5 / 54, unit_relvar = 0.54,
                   budget = 1000, fpc = "stage", N = 4900, N_psu = 100)
  expect_equal(fit$n[["n_per_psu"]], sqrt(12 * 49 / (3 * (5 - 49 / 49))))
  expect_equal(fit$n[["n_per_psu"]], 7)
  expect_equal(fit$n[["n_psu"]], 1000 / 33)
})

## A single-stage cluster population of four clusters of three units

.fp_z <- matrix(c(1, 4, 7, 3, 10, 5, 2, 8, 2, 4, 9, 2), nrow = 4, byrow = TRUE)

test_that("whole clusters beat SRS exactly below icc = 1 / (M + 1)", {
  vc <- varcomp(as.vector(t(.fp_z)), stage_id = list(rep(1:4, each = 3)))
  expect_lt(vc$icc, 1 / (3 + 1))
  cluster <- (1 - 2 / 4) * var(rowMeans(.fp_z)) / 2
  srs <- (1 - 6 / 12) * var(as.vector(.fp_z)) / 6
  expect_lt(cluster, srs)
  by_cluster <- prec_cluster(n = c(2, 3), icc = vc, fpc = "stage",
                             N = 12, N_psu = 4)
  expect_equal((by_cluster$cv * mean(.fp_z))^2, cluster, tolerance = 1e-12)
  fit <- n_cluster(stage_cost = c(5, 1), icc = vc, budget = 16,
                   fpc = "stage", N = 12, N_psu = 4)
  expect_equal(fit$n[["n_per_psu"]], 3)
})

test_that("a fixed whole-cluster take sizes the PSU count exactly", {
  vc <- varcomp(as.vector(t(.fp_z)), stage_id = list(rep(1:4, each = 3)))
  mu <- mean(.fp_z)
  S_B2 <- var(rowMeans(.fp_z))
  target <- 0.1
  fit <- n_cluster(stage_cost = c(5, 1), icc = vc, cv = target / mu,
                   n_per_psu = 3, fpc = "stage", N = 12, N_psu = 4)
  expect_equal(fit$n[["n_psu"]], S_B2 / (target^2 + S_B2 / 4),
               tolerance = 1e-10)
})

## Change between two occasions on a population of three units

test_that("full and independent panels give the enumerated change variances", {
  y1 <- c(12, 20, 40)
  y2 <- c(10, 20, 32)
  args <- list(var = c(var(y1), var(y2)), n = c(2, 2), N = 3)
  panel <- do.call(prec_change, c(args, list(overlap = 1,
                                             overlap_cor = cor(y1, y2))))
  indep <- do.call(prec_change, args)
  expect_equal(panel$se^2 * 9, 26, tolerance = 1e-12)
  expect_equal(indep$se^2 * 9, 494, tolerance = 1e-12)
})

## A retained panel and independent fresh draws on a population of four
##
## Each occasion pools the retained units with its own fresh units, drawn
## from the whole frame independently of the panel, so a unit can be drawn
## twice on one occasion. Every combination of draws is equally likely.

enumerate_draws <- function(y1, y2, k, f1, f2) {
  N <- length(y1)
  pick <- function(m) combn(N, m, simplify = FALSE)
  est <- list()
  for (r in pick(k)) for (i in pick(f1)) for (j in pick(f2)) {
    est[[length(est) + 1L]] <- c(mean(y1[c(r, i)]), mean(y2[c(r, j)]))
  }
  est <- do.call(rbind, est)
  pv <- function(x) mean((x - mean(x))^2)
  c(change = pv(est[, 2] - est[, 1]), average = pv(rowMeans(est)))
}

test_that("independent draws give the enumerated change variance", {
  y1 <- c(1, 2, 4, 9)
  y2 <- c(3, 5, 4, 12)
  for (cfg in list(c(1, 1, 1), c(2, 1, 1), c(1, 1, 2), c(3, 1, 1))) {
    k <- cfg[1]
    n <- k + cfg[2:3]
    got <- prec_change(var = c(var(y1), var(y2)), n = n, N = 4,
                       overlap = k / n[1], overlap_cor = cor(y1, y2),
                       overlap_draw = "independent")
    expect_equal(got$se^2, enumerate_draws(y1, y2, k, cfg[2], cfg[3])[[1]],
                 tolerance = 1e-12, info = paste(cfg, collapse = "-"))
  }
})

test_that("independent draws give the enumerated two-occasion average", {
  y1 <- c(1, 2, 4, 9)
  y2 <- c(2, 1, 4, 9)
  for (k in 1:3) {
    got <- prec_pooled(var = var(y1), n = k + 1, occasions = 2, N = 4,
                       overlap = k / (k + 1), overlap_cor = cor(y1, y2),
                       overlap_draw = "independent")
    expect_equal(got$se^2, enumerate_draws(y1, y2, k, 1, 1)[[2]],
                 tolerance = 1e-12, info = k)
  }
})

## Half of 40 units retained on a frame of 400, variances 100 and 164 and a
## covariance of 100, on the total scale

test_that("independent draws reproduce the partial-retention change", {
  args <- list(var = c(100, 164), n = 40, N = 400, overlap = 0.5,
               overlap_cor = 100 / sqrt(100 * 164))
  indep <- do.call(prec_change, c(args, overlap_draw = "independent"))
  fixed <- do.call(prec_change, args)
  expect_equal(indep$se^2 * 400^2, 623200, tolerance = 1e-12)
  expect_equal(fixed$se^2 * 400^2, 630400, tolerance = 1e-12)
})

## Two-phase sampling on ten units in two strata of five

test_that("an exact between component gives the finite SRS variance", {
  small <- data.frame(N = c(5, 5), sd = c(2, 2), mean = c(10, 20),
                      nu = c(1, 1))
  S2 <- (sum(4 * small$sd^2) + sum(5 * (small$mean - 15)^2)) / 9
  got <- prec_twophase(small, n_phase1 = 4, N = 10,
                       between = S2 - sum(0.5 * small$sd^2))
  expect_equal(got$se^2, S2 * (1 / 4 - 1 / 10), tolerance = 1e-12)
})
