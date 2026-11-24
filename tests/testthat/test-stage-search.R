## The compiled stage search against enumeration of small boxes

# Every whole design of a two-stage box, its cost and its relvariance.
stage_grid <- function(R, cost, upper) {
  g <- expand.grid(n1 = seq_len(upper[1]), n2 = seq_len(upper[2]))
  g$cost <- cost[1] * g$n1 + cost[2] * g$n1 * g$n2
  g$cv2 <- R[1] / g$n1 + R[2] / (g$n1 * g$n2)
  g
}

test_that("a missing warm start value leaves the budget search exact", {
  R <- matrix(c(0.04, 0.9), 1L)
  cost <- c(40, 5)
  upper <- c(30, 20)
  g <- stage_grid(R, cost, upper)
  budget <- 900
  want <- min(g$cv2[g$cost <= budget])
  for (n_init in list(c(NA_real_, 1), c(NaN, NaN), c(Inf, 3))) {
    r <- .stage_solve_discrete(R, matrix(1, 1, 2), matrix(Inf, 1, 2), 1,
                               cost, c(1, 1), upper, budget_var = budget,
                               n_init = n_init)
    expect_true(r$found)
    expect_false(anyNA(r$n))
    n <- as.numeric(r$n)
    expect_lte(cost[1] * n[1] + cost[2] * prod(n), budget)
    expect_equal(R[1] / n[1] + R[2] / prod(n), want, tolerance = 1e-12)
  }
})

test_that("a missing warm start value never makes a budget feasible", {
  r <- .stage_solve_discrete(matrix(c(0.04, 0.9), 1L), matrix(1, 1, 2),
                             matrix(Inf, 1, 2), 1, c(40, 5), c(1, 1),
                             c(30, 20), budget_var = 10,
                             n_init = c(NA_real_, 1))
  expect_false(r$found)
})

test_that("a missing warm start value leaves the cv search exact", {
  R <- matrix(c(0.04, 0.9), 1L)
  cost <- c(40, 5)
  upper <- c(30, 20)
  g <- stage_grid(R, cost, upper)
  want <- min(g$cost[g$cv2 <= 0.06^2 * (1 + 1e-9)])
  r <- .stage_solve_discrete(R, matrix(1, 1, 2), matrix(Inf, 1, 2), 0.06^2,
                             cost, c(1, 1), upper, n_init = c(NA_real_, 1))
  expect_true(r$found)
  n <- as.numeric(r$n)
  expect_equal(cost[1] * n[1] + cost[2] * prod(n), want)
})

test_that("the stage search refuses arguments of the wrong shape", {
  R <- matrix(c(0.04, 0.9), 1L)
  one <- matrix(1, 1, 2)
  inf <- matrix(Inf, 1, 2)
  call <- function(...) {
    a <- list(R, one, inf, 1, c(40, 5), c(1, 1), c(30, 20), 900, 0, NULL,
              2e7)
    b <- list(...)
    a[as.integer(names(b))] <- b
    do.call(.Call, c(list(C_stage_search), a))
  }
  expect_error(call(`2` = matrix(1, 1, 3)), "'rates' must be a 1 x 2")
  expect_error(call(`3` = matrix(1L, 1, 2)), "'M' must be a 1 x 2")
  expect_error(call(`4` = c(1, 1)), "'cv2_target' must be a double vector")
  expect_error(call(`5` = 40), "'cost' must be a double vector")
  expect_error(call(`7` = c(30, 20, 1)), "'upper' must be a double vector")
  expect_error(call(`10` = c(1, 2, 3)), "'n_init' must be a double vector")
  expect_error(call(`6` = c(NA_real_, 1)), "'lower' must be finite")
  expect_error(call(`6` = c(0, 1)), "'lower' must be finite")
  expect_error(call(`7` = c(NaN, 20)), "'upper' not NaN")
})

## Scale-free search: enumeration with the documented order of ties

# Every whole design of a box, with its cost, relvariances and tie-break keys.
stage_box <- function(R, rates, M, cost, upper) {
  S <- length(cost)
  g <- as.matrix(expand.grid(lapply(upper, seq_len)))
  P <- t(apply(g, 1L, cumprod))
  cost_all <- drop(P %*% cost)
  rbar <- t(apply(rates, 1L, cumprod))
  cv2 <- matrix(0, nrow(g), nrow(R))
  for (j in seq_len(nrow(R))) {
    for (s in seq_len(S)) {
      f <- if (is.finite(M[j, s])) pmax(0, 1 - g[, s] * rates[j, s] / M[j, s]) else 1
      cv2[, j] <- cv2[, j] + R[j, s] * f / (P[, s] * rbar[j, s])
    }
  }
  take <- if (S > 1L) rowSums(g[, -1L, drop = FALSE]) else 0
  list(n = g, cost = cost_all, cv2 = cv2, tn = P[, S], take = take)
}

# The design the documented order picks, or NULL when none is feasible.
stage_pick <- function(box, cv2_target = NULL, budget = NULL) {
  if (is.null(budget)) {
    ok <- apply(box$cv2, 1L, function(v) all(v <= cv2_target * (1 + 1e-9)))
    prim <- box$cost
  } else {
    ok <- box$cost <= budget * (1 + 1e-9)
    prim <- apply(box$cv2, 1L, max)
  }
  if (!any(ok)) return(NULL)
  i <- which(ok)
  best <- min(prim[i])
  i <- i[prim[i] <= best + 1e-12 * abs(best)]
  many <- length(i) > 1L
  if (is.null(budget)) i <- i[box$tn[i] == min(box$tn[i])]
  i <- i[box$take[i] == min(box$take[i])]
  list(n = box$n[i[1L], ], prim = best, many = many)
}

stage_random <- function(seed) {
  set.seed(seed)
  S <- sample(2:4, 1)
  J <- sample(1:2, 1)
  upper <- c(8, 6, 5, 4)[seq_len(S)]
  scale <- 10^stats::runif(1, -10, 3)
  cost <- sample(1:12, S, replace = TRUE) * scale
  R <- matrix(stats::runif(J * S, 0.05, 2), J) * 10^stats::runif(1, -10, 0)
  rates <- matrix(sample(c(1, 0.8, 0.5), J * S, replace = TRUE), J)
  M <- matrix(Inf, J, S)
  if (stats::runif(1) < 0.4) M[, S] <- 12
  list(R = R, rates = rates, M = M, cost = cost, upper = upper, S = S)
}

test_that("the stage search matches enumeration at any cost and variance scale", {
  for (seed in 1:40) {
    x <- stage_random(seed)
    box <- stage_box(x$R, x$rates, x$M, x$cost, x$upper)
    lower <- rep(1, x$S)
    tgt <- apply(box$cv2, 2L, stats::quantile, 0.3)
    want <- stage_pick(box, cv2_target = tgt)
    r <- .stage_solve_discrete(x$R, x$rates, x$M, tgt, x$cost, lower, x$upper)
    expect_identical(r$found, !is.null(want))
    if (!is.null(want)) {
      expect_identical(as.numeric(r$n), as.numeric(want$n), info = paste("seed", seed))
    }
    budget <- stats::quantile(box$cost, 0.4)
    want <- stage_pick(box, budget = budget)
    r <- .stage_solve_discrete(x$R, x$rates, x$M, rep(1, nrow(x$R)), x$cost,
                               lower, x$upper, budget_var = budget)
    expect_identical(r$found, !is.null(want))
    if (!is.null(want)) {
      expect_identical(as.numeric(r$n), as.numeric(want$n), info = paste("seed", seed))
    }
  }
})

test_that("a tie on cost goes to the smaller total size at every scale", {
  # (6, 1) and (4, 2) both cost 12, with total sizes 6 and 8.
  R <- matrix(c(1, 1), 1L)
  box <- stage_box(R, matrix(1, 1, 2), matrix(Inf, 1, 2), c(1, 1), c(8, 8))
  want <- stage_pick(box, cv2_target = 0.375)
  expect_true(want$many)
  expect_identical(as.numeric(want$n), c(6, 1))
  for (m in c(1e-10, 1, 1e6)) {
    r <- .stage_solve_discrete(R, matrix(1, 1, 2), matrix(Inf, 1, 2), 0.375,
                               c(1, 1) * m, c(1, 1), c(8, 8))
    expect_identical(as.numeric(r$n), c(6, 1))
  }
})

test_that("a cluster design does not move when costs or variances are rescaled", {
  for (sc in list(c(500, 50), c(800, 120, 20))) {
    icc <- if (length(sc) == 2L) 0.05 else c(0.05, 0.1)
    cv <- n_cluster(sc, icc = icc, cv = 0.05)$operational$n
    bud <- n_cluster(sc, icc = icc, budget = 1e5)$operational$n
    for (m in c(1e-10, 10)) {
      expect_identical(n_cluster(sc * m, icc = icc, cv = 0.05)$operational$n,
                       cv)
      expect_identical(
        n_cluster(sc * m, icc = icc, budget = 1e5 * m)$operational$n, bud
      )
    }
    expect_identical(
      n_cluster(sc, icc = icc, unit_relvar = 1e-10, budget = 1e5)$operational$n,
      bud
    )
  }
})

test_that("small cost scales neither overrun a budget nor pick a dearer design", {
  a <- n_cluster(stage_cost = c(1e-8, 1e-9), icc = 0.5, budget = 5e-8)
  expect_identical(unname(a$operational$n), c(4L, 2L))
  expect_lte(a$operational$cost, 5e-8)
  b <- n_cluster(stage_cost = c(1e-8, 1e-10), icc = 0.5, cv = 0.05)
  expect_identical(unname(b$operational$n), c(220L, 10L))
  expect_identical(
    unname(n_cluster(c(500, 50), icc = 0.05, unit_relvar = 1e-10,
                     budget = 1e5)$operational$n),
    c(80L, 15L)
  )
})

## The search's node limit is reported, never hidden

test_that("a fit records whether its whole-unit search finished", {
  f <- n_cluster(c(500, 50), icc = 0.05, cv = 0.05)
  expect_identical(f$operational$search$status, "optimal")
  expect_gt(f$operational$search$nodes, 0)
  lim <- n_cluster(c(800, 120, 20), icc = c(0.05, 0.1), cv = 0.05,
                   control = list(max_nodes = 5))
  expect_identical(lim$operational$search$status, "limit")
  expect_lte(lim$operational$search$nodes, 5)
  expect_output(print(lim), "stopped at the node limit of 5")
  expect_false(any(grepl("node limit", capture.output(print(f)))))
  again <- n_cluster(prec_cluster(lim))
  expect_identical(again$operational$search$status, "limit")
  g <- predict(lim, data.frame(cv = c(0.05, 0.06)))
  expect_equal(g$total_n[1], lim$total_n)
  expect_error(n_cluster(c(500, 50), icc = 0.05, cv = 0.05,
                         control = list(gap = 0.1)),
               "'control\\$gap' does not apply")
  expect_error(n_cluster(c(500, 50), icc = 0.05, cv = 0.05,
                         control = list(max_nodes = 0)),
               "max_nodes")
})

test_that("a search stopped before any design says feasibility is unknown", {
  R <- matrix(c(0.04, 0.3, 0.9), 1L)
  e <- tryCatch(
    .stage_solve_discrete(R, matrix(1, 1, 3), matrix(Inf, 1, 3), 1,
                          c(40, 10, 5), c(1, 1, 1), c(30, 20, 10),
                          budget_var = 2000, max_nodes = 1),
    error = identity
  )
  expect_s3_class(e, "svyplan_plan_unknown")
  expect_match(conditionMessage(e), "whether one exists is unknown")
})

test_that("a domain fit reports the limit of any of its domains", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.1), cv = c(0.1, 0.15),
                    icc_1 = c(0.02, 0.05), icc_2 = c(0.1, 0.08))
  tbl <- rbind(transform(ind, region = "N"),
               transform(ind, region = "S", p = c(0.4, 0.2)))
  ok <- n_cluster(tbl, stage_cost = c(800, 120, 20), domains = "region")
  expect_identical(ok$operational$search$status, "optimal")
  lim <- n_cluster(tbl, stage_cost = c(800, 120, 20), domains = "region",
                   control = list(max_nodes = 5))
  expect_identical(lim$operational$search$status, "limit")
  expect_named(lim$operational$search$by_domain, c("N", "S"))
  expect_output(print(lim), "node limit of 5 in")
  expect_identical(
    n_cluster(prec_cluster(lim))$operational$search$status, "limit"
  )
  expect_identical(
    n_cluster(ind, stage_cost = c(800, 120, 20),
              control = list(max_nodes = 5))$operational$search$status,
    "limit"
  )
})

test_that("a target just below a fixed stage's floor gets the floor error", {
  expect_error(
    n_cluster(c(500, 50), icc = 0.04, cv = 0.05 * (1 - 1e-13),
              stage_n = c(16, NA)),
    "below the achievable floor"
  )
})

## Counts past the integer range

test_that("a joint allocation past the integer range keeps its counts", {
  j <- expect_no_warning(n_alloc_joint(
    data.frame(stratum = "A", N = 1e10),
    measures = data.frame(stratum = "A", name = "y", p = 0.3),
    targets = data.frame(name = "y", cv = 0.1), min_n_stratum = 3e9
  ))
  expect_identical(j$operational$n, 3e9)
  expect_false(anyNA(j$detail$n_int))
  expect_true(j$operational$all_pass)
  pr <- prec_alloc_joint(j, n = j$detail$n_int)
  expect_equal(pr$params$achieved$n, 3e9)
  expect_output(print(j))
  below <- n_alloc_joint(
    data.frame(stratum = "A", N = 1e10),
    measures = data.frame(stratum = "A", name = "y", p = 0.3),
    targets = data.frame(name = "y", cv = 0.1),
    min_n_stratum = .Machine$integer.max - 10
  )
  expect_identical(below$operational$n, .Machine$integer.max - 10)
  expect_error(.whole_storage(2^53 + 2), "exceeds 2\\^53")
})

test_that("searches that count in integers refuse a larger design", {
  expect_error(n_cluster(c(1, 1e-6), icc = 0.5, cv = 1e-5),
               "counts at most 2,147,483,647 units at stage 1")
  expect_error(n_cluster(c(100, 1), icc = 0.05, budget = 1e12),
               "counts at most 2,147,483,647 units at stage 1")
  expect_error(n_alloc(data.frame(stratum = "A", N = 1e10, sd = 1), n = 3e9),
               "counts at most 2,147,483,647 units in total")
  expect_s3_class(n_cluster(c(100, 1), icc = 0.05, budget = 1e8),
                  "svyplan_cluster")
})

## Two whole-unit search details

test_that("opposing constraints still give the subgradient a plan", {
  # x must lie in [5, 7], and doubled equal multipliers stall near 4.18.
  q <- .whole_problem(o1 = 1, o2 = 0, g1 = matrix(c(0, 1 / 7), 1),
                      g2 = matrix(c(2.5, 0), 1), lo = 1, hi = 10,
                      g0 = matrix(c(0.5, 0), 1))
  expect_null(.whole_plan_from(q, c(0.5, 0.5)))
  r <- .whole_search(q)
  expect_identical(r$status, "optimal")
  expect_equal(r$value, 5)
})

test_that("the smaller root of a reduced-cost interval keeps its digits", {
  # Reference by bisection on a x^2 - top x + b over [0, top / (2 a)].
  ref <- function(a, b, top) {
    f <- function(x) a * x^2 - top * x + b
    lo <- 0
    hi <- top / (2 * a)
    for (i in 1:200) {
      mid <- (lo + hi) / 2
      if (f(mid) > 0) lo <- mid else hi <- mid
    }
    (lo + hi) / 2
  }
  for (cs in list(c(1e-14, 100, 1), c(1e-10, 5, 2), c(0.5, 0.3, 4))) {
    a <- cs[1]
    b <- cs[2]
    top <- cs[3]
    d <- sqrt(top^2 - 4 * a * b)
    expect_equal(.whole_small_root(a, b, top, d), ref(a, b, top),
                 tolerance = 1e-12)
  }
})

## Free takes past 1e5 in the multi-indicator solver

test_that("an indicator table holding one indicator twice gives the single design", {
  one <- n_cluster(stage_cost = c(1e8, 1), icc = 1e-6, budget = 1e10,
                   unit_relvar = 1)
  tab <- data.frame(name = c("a", "b"), mean = 1, sd = 1, cv = 0.1,
                    icc_1 = 1e-6)
  two <- n_cluster(tab, stage_cost = c(1e8, 1), budget = 1e10)
  expect_equal(unname(two$n), unname(one$n), tolerance = 1e-6)
  expect_equal(unname(two$n[2]), 9999995, tolerance = 1e-6)
  expect_equal(two$cv, one$cv, tolerance = 1e-6)
})

test_that("a supplied take cap above 1e5 binds in the multi-indicator solver", {
  tab <- data.frame(name = c("a", "b"), mean = 1, sd = 1, cv = 0.1,
                    icc_1 = 1e-6)
  capped <- n_cluster(tab, stage_cost = c(1e8, 1), budget = 1e10,
                      N = c(100, 5e8))
  single <- n_cluster(stage_cost = c(1e8, 1), icc = 1e-6, budget = 1e10,
                      unit_relvar = 1, N = c(100, 5e8))
  expect_equal(unname(capped$n[2]), 5e6, tolerance = 1e-9)
  expect_equal(unname(capped$n), unname(single$n), tolerance = 1e-6)
  loose <- n_cluster(tab, stage_cost = c(1e8, 1), budget = 1e10,
                     N = c(100, 2e9))
  expect_gt(unname(loose$n[2]), 1e7)
})

test_that("the multi-indicator cv mode finds takes past 1e5", {
  tab <- data.frame(name = c("a", "b"), mean = 1, sd = 1, cv = c(1e-3, 2e-3),
                    icc_1 = 1e-7)
  fit <- suppressWarnings(n_cluster(tab, stage_cost = c(1e9, 1)))
  single <- suppressWarnings(n_cluster(stage_cost = c(1e9, 1), icc = 1e-7,
                                       cv = 1e-3, unit_relvar = 1))
  expect_equal(unname(fit$n), unname(single$n), tolerance = 1e-6)
  expect_gt(unname(fit$n[2]), 1e5)
})

## What a stopped search reports

test_that("a stopped cluster search does not claim optimality anywhere", {
  f <- n_cluster(c(1, 1), icc = .5, budget = 1e8,
                 control = list(max_nodes = 1))
  out <- capture.output(print(f))
  expect_identical(out[1], "2-stage allocation")
  expect_match(paste(trimws(out), collapse = " "), "not proven most precise")
  expect_true(all(nchar(out) <= 80))
  sm <- capture.output(print(summary(f)))
  expect_identical(sm[1], "2-stage allocation")
  expect_true(any(grepl("search: stopped at the node limit", sm)))
  grid <- predict(f, data.frame(budget = c(1e3, 1e8)))
  expect_identical(grid$.status[2], "limit")
  expect_equal(grid$.nodes[2], 1)

  g <- n_cluster(stage_cost = c(500, 50), icc = .05, cv = .05)
  expect_identical(capture.output(print(g))[1], "Optimal 2-stage allocation")
  expect_identical(capture.output(print(summary(g)))[1],
                   "Optimal 2-stage allocation")
  expect_identical(predict(g, data.frame(cv = c(.05, .1)))$.status,
                   c("optimal", "optimal"))
})

test_that("a target search stopped at its limit is not proven cheapest", {
  tab <- data.frame(name = c("y", "z"), p = c(.3, .2), cv = c(.05, .08),
                    icc_1 = .05, icc_2 = .1)
  fit <- n_cluster(tab, stage_cost = c(200, 40, 10),
                   control = list(max_nodes = 1))
  out <- capture.output(print(fit))
  expect_identical(fit$operational$search$status, "limit")
  expect_match(out[1], "^Multi-indicator allocation")
  expect_match(paste(trimws(out), collapse = " "), "not proven cheapest")
})

test_that("a scenario whose search finds nothing reads unknown", {
  g <- n_cluster(stage_cost = c(500, 50), icc = .05, cv = .05)
  local_mocked_bindings(n_cluster.default = function(...) {
    stop(structure(list(message = "no design", call = NULL),
                   class = c("svyplan_plan_unknown", "error", "condition")))
  })
  grid <- predict(g, data.frame(cv = c(.05, .1)))
  expect_identical(grid$.status, c("unknown", "unknown"))
  expect_true(all(is.na(grid$total_n)))
})

test_that("a joint search record counts each domain's nodes", {
  di <- data.frame(region = c("A", "B"), name = "y", p = c(.2, .4),
                   cv = c(.3, .2), icc_1 = c(.1, .2))
  fit <- suppressWarnings(n_cluster(di, stage_cost = c(4, 1),
                                    domains = "region", allocation = "joint",
                                    budget = 12))
  expect_named(fit$operational$search$nodes_by_domain, c("A", "B"))
  expect_true(all(fit$operational$search$nodes_by_domain >= 0))
})

test_that("a several-indicator scenario whose search finds nothing reads unknown", {
  tab <- data.frame(name = c("y", "z"), p = c(.3, .2), cv = c(.1, .15),
                    icc_1 = .05, icc_2 = .1)
  m <- n_cluster(tab, stage_cost = c(200, 40, 10), budget = 30000)
  local_mocked_bindings(n_cluster.data.frame = function(...) {
    stop(structure(list(message = "no design", call = NULL),
                   class = c("svyplan_plan_unknown", "error", "condition")))
  })
  grid <- predict(m, data.frame(budget = c(20000, 30000)))
  expect_identical(grid$.status, c("unknown", "unknown"))
  expect_true(all(is.na(grid$binding)))
})
