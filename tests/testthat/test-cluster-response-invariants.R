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
  ind <- data.frame(name = "a", p = p, icc_1 = icc, cv = target,
                    resp_rate = rr, resp_rate_1 = rp)
  multi_args <- c(list(x = ind, stage_cost = cost), fixed)
  cluster_args <- c(
    list(stage_cost = cost, icc = icc, unit_relvar = (1 - p) / p,
         resp_rate = c(rp, rr)),
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
      data.frame(name = "a", p = 0.3, cv = 0.05, icc_1 = 0.05,
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
      ind <- data.frame(name = "a", p = 0.3, cv = 0.10, icc_1 = 0.04,
                        resp_rate = rr)
      cost <- if (stages == 2L) c(500, 50) else c(500, 100, 50)
      if (stages == 3L) ind$icc_2 <- 0.06
      fit <- n_cluster(ind, stage_cost = cost, budget = 1e5)
      expect_equal(fit$cv, prec_cluster(fit)$cv, tolerance = 1e-8,
                   info = sprintf("%d-stage, resp_rate=%.2f", stages, rr))
    }
  }
})

test_that("CV-mode results already agree, and stay agreeing", {
  for (rr in c(1, 0.5)) {
    ind <- data.frame(name = "a", p = 0.3, cv = 0.10, icc_1 = 0.04,
                      resp_rate = rr)
    fit <- n_cluster(ind, stage_cost = c(500, 50))
    expect_equal(fit$cv, prec_cluster(fit)$cv, tolerance = 1e-8)
  }
})

## Fixed stage sizes

test_that("a fixed stage_1 is respected by the continuous two-stage solution", {
  fit <- n_cluster(
    data.frame(name = "a", p = 0.3, cv = 0.08, icc_1 = 0.05),
    stage_cost = c(500, 50), stage_n = c(30, NA)
  )
  expect_equal(fit$n[[1L]], 30, tolerance = 1e-8)
  expect_equal(fit$operational$n[[1L]], 30)
  # The take must be the one that reaches the target at that PSU count, not
  # the unconstrained cost optimum.
  expect_equal(prec_cluster(fit)$cv, 0.08, tolerance = 1e-6)
})

test_that("a fixed stage_1 matches n_cluster, response or not", {
  for (rr in c(1, 0.6)) {
    got <- mc_vs_cluster("cv", rr = rr, target = 0.10, stage_n = c(40, NA))
    expect_equal(got$multi$n[[1L]], 40, tolerance = 1e-8)
    expect_equal(got$multi$n[[2L]], got$cluster$n[[2L]], tolerance = 1e-5,
                 info = sprintf("resp_rate=%.2f", rr))
  }
})

## Joint-domain budget allocation

test_that("joint-domain budget allocation reads the ultimate response rate", {
  mk <- function(rr) {
    data.frame(name = rep("a", 2L), dom = c("A", "B"), p = 0.3,
               cv = 0.10, icc_1 = 0.05, resp_rate = rr)
  }
  full <- n_cluster(mk(1), stage_cost = c(500, 50), domains = "dom",
                          budget = 2e5, allocation = "joint")
  half <- n_cluster(mk(0.5), stage_cost = c(500, 50), domains = "dom",
                          budget = 2e5, allocation = "joint")
  # Half the ultimate units respond, so the same budget buys less precision.
  expect_gt(max(half$domains$.cv), max(full$domains$.cv))
})

## The classic clustered allocator's operational search

test_that("clustered n_alloc integerizes on the responding take", {
  frame <- data.frame(stratum = "s1", N = 1e6, sd = 1, mean = 1,
                      icc_1 = 0.05, take_2 = 20, resp_rate = 0.5)
  fit <- n_alloc(frame, cv = 0.05)
  # 58 PSUs of 20 issued units, half responding, meet the target exactly under
  # the documented expected-take model.
  expect_equal(fit$detail$n_1_int, 58)
  responding <- 20 * 0.5
  expect_equal(
    sqrt((1 + 0.05 * (responding - 1)) / (fit$detail$n_1_int * responding)),
    0.05,
    tolerance = 1e-6
  )
})

test_that("the operational cluster CV is the one prec_alloc confirms", {
  frame <- data.frame(stratum = c("s1", "s2"), N = c(4e5, 6e5),
                      sd = c(1, 1.4), mean = c(1, 1.2), icc_1 = 0.05,
                      take_2 = c(20, 15), resp_rate = c(0.5, 0.8))
  fit <- n_alloc(frame, cv = 0.06)
  expect_true(all(fit$detail$n_1_int <= fit$detail$n_1 + 1))
  expect_lte(prec_alloc(fit, n = fit$detail$n_int)$cv, 0.06 * (1 + 1e-6))
})

## Three-stage fixed-stage inversions
##
## Every allowed combination, not only the two-stage fixed-PSU case. Each
## inversion is read off search_coef, which already carries every stage's
## response rate, so a target that is attainable at the fixed count is met
## exactly rather than approached.

test_that("every three-stage fixed-stage combination meets an attainable target", {
  base <- data.frame(name = "a", p = 0.3, cv = 0.10, icc_1 = 0.02,
                     icc_2 = 0.05, resp_rate_1 = 0.9,
                     resp_rate_2 = 0.8, resp_rate = 0.7)
  specs <- list(
    c(20, NA, NA),
    c(20, NA, 5),
    c(20, 10, NA),
    c(NA, 8, NA),
    c(NA, NA, 4),
    c(NA, 8, 4)
  )
  for (sn in specs) {
    fit <- n_cluster(base, stage_cost = c(500, 100, 50), stage_n = sn)
    label <- paste(sn, collapse = ",")
    expect_equal(fit$cv, 0.10, tolerance = 1e-6, info = label)
    expect_equal(prec_cluster(fit)$cv, 0.10, tolerance = 1e-6,
                 info = label)
    for (s in which(!is.na(sn))) {
      expect_equal(unname(fit$n[s]), sn[s], tolerance = 1e-8, info = label)
    }
  }
})

test_that("a fixed three-stage count below the PSU floor is refused", {
  base <- data.frame(name = "a", p = 0.3, cv = 0.02, icc_1 = 0.10,
                     icc_2 = 0.05, resp_rate = 0.7)
  expect_error(
    n_cluster(base, stage_cost = c(500, 100, 50), stage_n = c(3, NA, NA)),
    "floor|below|achievable"
  )
})

## PSU response must not enter the clustering bracket

test_that("PSU response leaves the integer take alone", {
  takes <- vapply(c(1, 0.75, 0.5), function(rp) {
    frame <- data.frame(stratum = "s1", N = 1e6, sd = 1, mean = 1,
                        icc_1 = 0.05, cost_1 = 500, cost_2 = 50,
                        resp_rate = 0.8, resp_rate_1 = rp)
    fit <- n_alloc(frame, cv = 0.05)
    c(fit$detail$take_2, fit$detail$take_2_int)
  }, numeric(2L))
  # The continuous optimum does not move with PSU response, and neither may
  # the whole take the operational search rounds it to.
  expect_equal(diff(range(takes[1L, ])), 0, tolerance = 1e-9)
  expect_equal(diff(range(takes[2L, ])), 0)
})

test_that("the operational cluster CV reads the ultimate rate in the bracket", {
  frame <- data.frame(stratum = "s1", N = 1e6, sd = 1, mean = 1,
                      icc_1 = 0.05, take_2 = 20,
                      resp_rate = 0.8, resp_rate_1 = 0.5)
  fit <- n_alloc(frame, cv = 0.05)
  responding <- 20 * 0.8
  n_net <- fit$detail$n_1_int * responding * 0.5
  expect_equal(
    fit$operational$cv,
    sqrt((1 + 0.05 * (responding - 1)) / n_net *
           (1 - fit$detail$n_1_int * 20 * 0.8 * 0.5 / 1e6)),
    tolerance = 1e-4
  )
})

## n_cluster's conditional final-stage optimum

test_that("a fixed middle take gets the conditional final-stage optimum", {
  for (ps in c(5, 10, 20)) {
    fit <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05),
                     cv = 0.10, stage_n = c(NA, ps, NA))
    cost_at <- function(q) {
      p <- prec_cluster(n = c(1, ps, q), icc = c(0.02, 0.05))
      (p$cv / 0.10)^2 * (500 + 100 * ps + 50 * ps * q)
    }
    expect_equal(fit$n[["stage_3"]],
                 stats::optimize(cost_at, c(1.0001, 200), tol = 1e-10)$minimum,
                 tolerance = 1e-5, info = sprintf("take_2=%d", ps))
  }
})

test_that("both stages free keeps the unrestricted optimum", {
  fit <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05),
                   cv = 0.10)
  expect_equal(fit$n[["stage_3"]], sqrt((1 - 0.05) / 0.05 * 100 / 50),
               tolerance = 1e-6)
})

## Stored stage response rates

test_that("cluster plans store stage response rates by name", {
  two <- n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = 1,
                   cv = 0.05, resp_rate = c(0.9, 0.8))
  three <- n_cluster(stage_cost = c(ea = 500, hh = 50, p = 5),
                     icc = c(0.05, 0.1), var_ratio = c(1, 1), cv = 0.05,
                     resp_rate = c(0.9, 0.95, 0.8))
  scalar <- n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = 1,
                      cv = 0.05, resp_rate = 0.8)
  multi <- n_cluster(
    data.frame(name = c("a", "b"), p = c(0.3, 0.5), icc_1 = c(0.05, 0.08),
               cv = c(0.08, 0.06)),
    stage_cost = c(500, 50), resp_rate = c(0.9, 0.85)
  )
  prec <- prec_cluster(n = c(40, 12), stage_cost = c(500, 50), icc = 0.05,
                       var_ratio = 1, unit_relvar = 1,
                       resp_rate = c(0.9, 0.8))
  expect_identical(two$params$resp_rate, c(resp_rate_1 = 0.9, resp_rate = 0.8))
  expect_identical(three$params$resp_rate,
                   c(resp_rate_1 = 0.9, resp_rate_2 = 0.95, resp_rate = 0.8))
  expect_identical(scalar$params$resp_rate,
                   c(resp_rate_1 = 1, resp_rate = 0.8))
  expect_identical(multi$params$resp_rate,
                   c(resp_rate_1 = 0.9, resp_rate = 0.85))
  expect_identical(prec$params$resp_rate,
                   c(resp_rate_1 = 0.9, resp_rate = 0.8))
})

test_that("stored rate names keep each rate on its stage when read back", {
  ## Swapped rates give another design, so a rate read back on the wrong
  ## stage changes every comparison below.
  fit_with <- function(rr, ...) {
    n_cluster(stage_cost = c(ea = 500, hh = 50, p = 5), icc = c(0.05, 0.1),
              var_ratio = c(1, 1), cv = 0.05, resp_rate = rr, ...)
  }
  fit <- fit_with(c(0.6, 0.95, 0.8))
  expect_false(isTRUE(all.equal(fit$n, fit_with(c(0.8, 0.95, 0.6))$n)))

  again <- n_cluster(prec_cluster(fit))
  expect_equal(again$n, fit$n)
  expect_identical(again$params$resp_rate, fit$params$resp_rate)

  grid <- predict(fit, data.frame(resp_rate = 0.7))
  expect_equal(unlist(grid[1, c("ea", "hh", "p")]),
               fit_with(c(0.6, 0.95, 0.7))$n)

  direct <- prec_cluster(n = fit$n, stage_cost = c(ea = 500, hh = 50, p = 5),
                         icc = c(0.05, 0.1), var_ratio = c(1, 1),
                         unit_relvar = 1, resp_rate = c(0.6, 0.95, 0.7))
  expect_equal(as.double(effective_n(fit, resp_rate = 0.7)),
               as.double(effective_n(direct)))
  expect_false(isTRUE(all.equal(as.double(effective_n(fit, resp_rate = 0.7)),
                                as.double(effective_n(fit)))))
})

test_that("a cluster summary prints with named stage rates", {
  fit <- n_cluster(stage_cost = c(500, 50), icc = 0.05, var_ratio = 1,
                   cv = 0.05, resp_rate = c(0.9, 0.8))
  s <- summary(fit)
  expect_identical(names(s$field),
                   c("respondents", "net_take", "deff", "effective", "df"))
  expect_equal(s$field[["net_take"]], 14 * 0.8)
  out <- capture.output(print(s))
  expect_true(any(grepl("response rate: 0.9, 0.8", out, fixed = TRUE)))
})
