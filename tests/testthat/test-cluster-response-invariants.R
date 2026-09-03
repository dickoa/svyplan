## Stage-specific response in the clustered planners
##
## One indicator is one design, so n_cluster() must return what
## n_cluster() returns for the same problem. That invariant is what catches a
## branch that carries the response rates into the variance it reports but not
## into the optimization that chose the design.

## Helpers

mc_vs_cluster <- function(mode = c("cv", "budget"), rr = 1, rp = 1,
                          icc = 0.05, p = 0.3, cost = c(500, 50),
                          target = 0.05, budget = 1e5, ...) {
  mode <- match.arg(mode)
  fixed <- list(...)
  ind <- data.frame(name = "a", p = p, icc_psu = icc, cv = target,
                    resp_rate = rr, resp_rate_psu = rp)
  multi_args <- c(list(indicators = ind, stage_cost = cost), fixed)
  cluster_args <- c(
    list(stage_cost = cost, icc = icc, unit_relvar = (1 - p) / p,
         resp_rate = rr, resp_rate_psu = rp),
    fixed
  )
  if (mode == "cv") {
    cluster_args$cv <- target
  } else {
    multi_args$budget <- budget
    cluster_args$budget <- budget
  }
  list(multi = do.call(n_cluster, multi_args),
       cluster = do.call(n_cluster, cluster_args))
}

## The two-stage optimum

test_that("two-stage CV mode matches n_cluster at every response rate", {
  for (rr in c(1, 0.85, 0.5)) {
    for (rp in c(1, 0.9)) {
      got <- mc_vs_cluster("cv", rr = rr, rp = rp)
      expect_equal(got$multi$n[[2L]], got$cluster$n[[2L]], tolerance = 1e-6,
                   info = sprintf("take, resp_rate=%.2f resp_rate_psu=%.2f",
                                  rr, rp))
      expect_equal(got$multi$n[[1L]], got$cluster$n[[1L]], tolerance = 1e-6,
                   info = sprintf("n_psu, resp_rate=%.2f", rr))
      expect_equal(got$multi$cost, got$cluster$cost, tolerance = 1e-5)
    }
  }
})

test_that("the two-stage take follows the response-aware closed form", {
  b_star <- function(c1, c2, icc, r) sqrt(c1 * (1 - icc) / (c2 * icc * r))
  # 0.2 and 0.1 matter on their own: the optimum grows as 1 / sqrt(r), so
  # below r = 0.25 it escapes any bracket drawn at the gross scale, however
  # generous the safety margin on it.
  for (rr in c(1, 0.85, 0.5, 0.2, 0.1)) {
    fit <- n_cluster(
      indicators = data.frame(name = "a", p = 0.3, cv = 0.05, icc_psu = 0.05,
                 resp_rate = rr),
      stage_cost = c(500, 50)
    )
    expect_equal(fit$n[[2L]], b_star(500, 50, 0.05, rr), tolerance = 1e-6,
                 info = sprintf("resp_rate=%.2f", rr))
  }
})

test_that("two-stage budget mode matches n_cluster at every response rate", {
  for (rr in c(1, 0.85, 0.5)) {
    got <- mc_vs_cluster("budget", rr = rr)
    expect_equal(got$multi$n[[2L]], got$cluster$n[[2L]], tolerance = 1e-4,
                 info = sprintf("take, resp_rate=%.2f", rr))
    expect_equal(got$multi$n[[1L]], got$cluster$n[[1L]], tolerance = 1e-4)
  }
})

## The reported CV must be the one the precision function confirms

test_that("budget-mode CV is the CV prec_cluster() recomputes", {
  for (stages in 2:3) {
    for (rr in c(1, 0.85, 0.5)) {
      ind <- data.frame(name = "a", p = 0.3, cv = 0.10, icc_psu = 0.04,
                        resp_rate = rr)
      cost <- if (stages == 2L) c(500, 50) else c(500, 100, 50)
      if (stages == 3L) ind$icc_ssu <- 0.06
      fit <- n_cluster(indicators = ind, stage_cost = cost, budget = 1e5)
      expect_equal(fit$cv, prec_cluster(fit)$cv, tolerance = 1e-8,
                   info = sprintf("%d-stage, resp_rate=%.2f", stages, rr))
    }
  }
})

test_that("CV-mode results already agree, and stay agreeing", {
  for (rr in c(1, 0.5)) {
    ind <- data.frame(name = "a", p = 0.3, cv = 0.10, icc_psu = 0.04,
                      resp_rate = rr)
    fit <- n_cluster(indicators = ind, stage_cost = c(500, 50))
    expect_equal(fit$cv, prec_cluster(fit)$cv, tolerance = 1e-8)
  }
})

## Fixed stage sizes

test_that("a fixed n_psu is respected by the continuous two-stage solution", {
  fit <- n_cluster(
    indicators = data.frame(name = "a", p = 0.3, cv = 0.08, icc_psu = 0.05),
    stage_cost = c(500, 50), n_psu = 30
  )
  expect_equal(fit$n[[1L]], 30, tolerance = 1e-8)
  expect_equal(fit$operational$n[[1L]], 30)
  # The take must be the one that reaches the target at that PSU count, not
  # the unconstrained cost optimum.
  expect_equal(prec_cluster(fit)$cv, 0.08, tolerance = 1e-6)
})

test_that("a fixed n_psu matches n_cluster, response or not", {
  for (rr in c(1, 0.6)) {
    got <- mc_vs_cluster("cv", rr = rr, target = 0.10, n_psu = 40)
    expect_equal(got$multi$n[[1L]], 40, tolerance = 1e-8)
    expect_equal(got$multi$n[[2L]], got$cluster$n[[2L]], tolerance = 1e-5,
                 info = sprintf("resp_rate=%.2f", rr))
  }
})

## Joint-domain budget allocation

test_that("joint-domain budget allocation reads the ultimate response rate", {
  mk <- function(rr) {
    data.frame(name = rep("a", 2L), dom = c("A", "B"), p = 0.3,
               cv = 0.10, icc_psu = 0.05, resp_rate = rr)
  }
  full <- n_cluster(indicators = mk(1), stage_cost = c(500, 50), domains = "dom",
                          budget = 2e5, allocation = "joint")
  half <- n_cluster(indicators = mk(0.5), stage_cost = c(500, 50), domains = "dom",
                          budget = 2e5, allocation = "joint")
  # Half the ultimate units respond, so the same budget buys less precision.
  expect_gt(max(half$domains$.cv), max(full$domains$.cv))
})

## The classic clustered allocator's operational search

test_that("clustered n_alloc integerizes on the responding take", {
  frame <- data.frame(stratum = "s1", N = 1e6, sd = 1, mean = 1,
                      icc_psu = 0.05, n_per_psu = 20, resp_rate = 0.5)
  fit <- n_alloc(frame, cv = 0.05)
  # 58 PSUs of 20 issued units, half responding, meet the target exactly under
  # the documented expected-take model.
  expect_equal(fit$detail$n_psu_int, 58)
  responding <- 20 * 0.5
  expect_equal(
    sqrt((1 + 0.05 * (responding - 1)) / (fit$detail$n_psu_int * responding)),
    0.05,
    tolerance = 1e-6
  )
})

test_that("the operational cluster CV is the one prec_alloc confirms", {
  frame <- data.frame(stratum = c("s1", "s2"), N = c(4e5, 6e5),
                      sd = c(1, 1.4), mean = c(1, 1.2), icc_psu = 0.05,
                      n_per_psu = c(20, 15), resp_rate = c(0.5, 0.8))
  fit <- n_alloc(frame, cv = 0.06)
  expect_true(all(fit$detail$n_psu_int <= fit$detail$n_psu + 1))
  expect_lte(prec_alloc(fit, n = fit$detail$n_int)$cv, 0.06 * (1 + 1e-6))
})

## Three-stage fixed-stage inversions
##
## Every allowed combination, not only the two-stage fixed-PSU case. Each
## inversion is read off search_coef, which already carries every stage's
## response rate, so a target that is attainable at the fixed count is met
## exactly rather than approached.

test_that("every three-stage fixed-stage combination meets an attainable target", {
  base <- data.frame(name = "a", p = 0.3, cv = 0.10, icc_psu = 0.02,
                     icc_ssu = 0.05, resp_rate_psu = 0.9,
                     resp_rate_ssu = 0.8, resp_rate = 0.7)
  specs <- list(
    list(n_psu = 20),
    list(n_psu = 20, n_per_ssu = 5),
    list(n_psu = 20, n_per_psu = 10),
    list(n_per_psu = 8),
    list(n_per_ssu = 4),
    list(n_per_psu = 8, n_per_ssu = 4)
  )
  for (spec in specs) {
    fit <- do.call(n_cluster,
                   c(list(indicators = base, stage_cost = c(500, 100, 50)),
                     spec))
    label <- paste(names(spec), collapse = "+")
    expect_equal(fit$cv, 0.10, tolerance = 1e-6, info = label)
    expect_equal(prec_cluster(fit)$cv, 0.10, tolerance = 1e-6,
                 info = label)
    for (nm in names(spec)) {
      expect_equal(fit$n[[nm]], spec[[nm]], tolerance = 1e-8, info = label)
    }
  }
})

test_that("a fixed three-stage count below the PSU floor is refused", {
  base <- data.frame(name = "a", p = 0.3, cv = 0.02, icc_psu = 0.10,
                     icc_ssu = 0.05, resp_rate = 0.7)
  expect_error(
    n_cluster(indicators = base, stage_cost = c(500, 100, 50), n_psu = 3),
    "floor|below|achievable"
  )
})

## PSU response must not enter the clustering bracket

test_that("PSU response leaves the integer take alone", {
  takes <- vapply(c(1, 0.75, 0.5), function(rp) {
    frame <- data.frame(stratum = "s1", N = 1e6, sd = 1, mean = 1,
                        icc_psu = 0.05, cost_psu = 500, cost_ssu = 50,
                        resp_rate = 0.8, resp_rate_psu = rp)
    fit <- n_alloc(frame, cv = 0.05)
    c(fit$detail$n_per_psu, fit$detail$n_per_psu_int)
  }, numeric(2L))
  # The continuous optimum does not move with PSU response, and neither may
  # the whole take the operational search rounds it to.
  expect_equal(diff(range(takes[1L, ])), 0, tolerance = 1e-9)
  expect_equal(diff(range(takes[2L, ])), 0)
})

test_that("the operational cluster CV reads the ultimate rate in the bracket", {
  frame <- data.frame(stratum = "s1", N = 1e6, sd = 1, mean = 1,
                      icc_psu = 0.05, n_per_psu = 20,
                      resp_rate = 0.8, resp_rate_psu = 0.5)
  fit <- n_alloc(frame, cv = 0.05)
  responding <- 20 * 0.8
  n_net <- fit$detail$n_psu_int * responding * 0.5
  expect_equal(
    fit$operational$cv,
    sqrt((1 + 0.05 * (responding - 1)) / n_net *
           (1 - fit$detail$n_psu_int * 20 * 0.8 * 0.5 / 1e6)),
    tolerance = 1e-4
  )
})

## n_cluster's conditional final-stage optimum

test_that("a fixed middle take gets the conditional final-stage optimum", {
  for (ps in c(5, 10, 20)) {
    fit <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05),
                     cv = 0.10, n_per_psu = ps)
    cost_at <- function(q) {
      p <- prec_cluster(n = c(1, ps, q), icc = c(0.02, 0.05))
      (p$cv / 0.10)^2 * (500 + 100 * ps + 50 * ps * q)
    }
    expect_equal(fit$n[["n_per_ssu"]],
                 stats::optimize(cost_at, c(1.0001, 200), tol = 1e-10)$minimum,
                 tolerance = 1e-5, info = sprintf("n_per_psu=%d", ps))
  }
})

test_that("both stages free keeps the unrestricted optimum", {
  fit <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05),
                   cv = 0.10)
  expect_equal(fit$n[["n_per_ssu"]], sqrt((1 - 0.05) / 0.05 * 100 / 50),
               tolerance = 1e-6)
})
