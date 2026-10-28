## The menu search against every plan of small problems

# The cheapest plan meeting every column, over every plan.
menu_brute <- function(p) {
  grid <- as.matrix(expand.grid(p$groups))
  cost <- rowSums(matrix(p$cost[grid], nrow(grid)))
  ok <- rep(TRUE, nrow(grid))
  for (k in seq_len(ncol(p$C))) {
    ok <- ok & rowSums(matrix(p$C[grid, k], nrow(grid))) <= 1 + 1e-9
  }
  if (any(ok)) min(cost[ok]) else Inf
}

# A dearer state loads every column less, the shape of a size menu, and
# each column skips some groups, as a domain target does. With `shaped =
# FALSE` loads are drawn apart from cost, so many states are dominated.
menu_random <- function(seed, H = sample(2:5, 1), size = 7L,
                        integer = TRUE, shaped = TRUE) {
  set.seed(seed)
  K <- sample(1:3, 1)
  owner <- rep(seq_len(H), sample(seq_len(size), H, replace = TRUE))
  x <- stats::runif(length(owner), 0.2, 3)
  cost <- if (integer) round(10 * x) else 10 * x
  if (!shaped) x <- stats::runif(length(owner), 0.2, 3)
  C <- matrix(stats::runif(length(owner) * K, 0.5, 1.5), ncol = K) /
    (x * H) * stats::runif(1, 0.3, 1.6)
  for (k in seq_len(K)) {
    if (stats::runif(1) < 0.4) C[owner == sample(H, 1), k] <- 0
  }
  .menu_problem(cost, C, owner)
}

check_plan <- function(p, r, want) {
  expect_equal(r$status, "optimal")
  expect_equal(r$value, want, tolerance = 1e-9)
  expect_equal(sum(p$cost[r$pick]), r$value, tolerance = 1e-12)
  expect_true(.menu_feasible(p, r$pick))
  expect_identical(p$owner[r$pick], seq_len(p$H))
}

test_that("the search proves the optimum of every small problem", {
  infeasible <- 0
  for (seed in 1:150) {
    integer <- seed %% 3 != 0
    p <- menu_random(seed, integer = integer)
    want <- menu_brute(p)
    sg <- .menu_subgradient(p)
    r <- .menu_search(p, sg$lambda, sg$pick, grain = if (integer) 1 else 0)
    if (is.finite(want)) {
      check_plan(p, r, want)
    } else {
      infeasible <- infeasible + 1
      expect_equal(r$status, "infeasible")
      expect_null(r$pick)
    }
  }
  expect_gt(infeasible, 0)
})

test_that("the search proves the optimum where it has to branch", {
  deep <- 0
  for (seed in 1:80) {
    integer <- seed %% 2 == 0
    p <- menu_random(seed, H = 8L, size = 4L, integer = integer)
    want <- menu_brute(p)
    sg <- .menu_subgradient(p)
    r <- .menu_search(p, sg$lambda, sg$pick, grain = if (integer) 1 else 0)
    if (is.finite(want)) check_plan(p, r, want) else
      expect_equal(r$status, "infeasible")
    deep <- deep + (r$nodes > 5)
  }
  expect_gt(deep, 5)
})

test_that("the search finds the optimum with no plan to start from", {
  # Without an incumbent, a node or state cut too early loses the optimum.
  for (seed in 1:80) {
    integer <- seed %% 2 == 0
    p <- menu_random(seed, H = if (seed > 40) 8L else sample(2:5, 1),
                     size = if (seed > 40) 4L else 7L, integer = integer)
    want <- menu_brute(p)
    if (!is.finite(want)) next
    grain <- if (integer) 1 else 0
    lambda <- .menu_subgradient(p)$lambda
    check_plan(p, .menu_search(p, lambda, NULL, grain = grain), want)
    check_plan(p, .menu_search(p, rep(0, ncol(p$C)), NULL, grain = grain),
               want)
  }
})

test_that("the search proves the optimum where many states are dominated", {
  for (seed in 1:80) {
    integer <- seed %% 2 == 0
    p <- menu_random(seed, H = 8L, size = 5L, integer = integer,
                     shaped = FALSE)
    want <- menu_brute(p)
    if (!is.finite(want)) next
    grain <- if (integer) 1 else 0
    lambda <- .menu_subgradient(p)$lambda
    check_plan(p, .menu_search(p, lambda, NULL, grain = grain), want)
  }
})

test_that("the multipliers never return a plan dearer than their start", {
  for (seed in 1:40) {
    p <- menu_random(seed, H = 8L, size = 4L)
    want <- menu_brute(p)
    if (!is.finite(want)) next
    best <- .menu_search(p, .menu_subgradient(p)$lambda, NULL, grain = 1)
    for (warm in c(FALSE, TRUE)) {
      sg <- .menu_subgradient(p, lambda = if (warm) rep(1, ncol(p$C)),
                              start = best$pick)
      expect_equal(sg$ub, want)
      expect_equal(sum(p$cost[sg$pick]), sg$ub)
      expect_lte(sg$L, want + 1e-9 * want)
    }
  }
})

test_that("past the open-list cap the search goes on depth-first", {
  for (seed in 1:40) {
    p <- menu_random(seed, H = 6L, size = 6L)
    want <- menu_brute(p)
    if (!is.finite(want)) next
    sg <- .menu_subgradient(p)
    check_plan(p, .menu_search(p, sg$lambda, sg$pick, grain = 1,
                               max_open = 1), want)
  }
})

test_that("a node limit returns a plan with a bound below the optimum", {
  hit <- 0
  for (seed in 1:40) {
    p <- menu_random(seed, H = 6L, size = 6L, integer = FALSE)
    want <- menu_brute(p)
    if (!is.finite(want)) next
    sg <- .menu_subgradient(p)
    r <- .menu_search(p, sg$lambda, sg$pick, max_nodes = 3)
    if (r$status == "limit") hit <- hit + 1
    expect_true(r$status %in% c("optimal", "limit"))
    expect_true(.menu_feasible(p, r$pick))
    expect_equal(sum(p$cost[r$pick]), r$value, tolerance = 1e-12)
    expect_lte(r$bound, want + 1e-9 * want)
    expect_gte(r$value, want - 1e-9 * want)
    if (r$status == "optimal") expect_equal(r$value, want, tolerance = 1e-9)
  }
  expect_gt(hit, 0)
})

test_that("groups that share no column are solved apart", {
  for (seed in 1:30) {
    a <- menu_random(seed, H = 3L)
    b <- menu_random(seed + 1000, H = 3L)
    C <- rbind(cbind(a$C, matrix(0, nrow(a$C), ncol(b$C))),
               cbind(matrix(0, nrow(b$C), ncol(a$C)), b$C))
    p <- .menu_problem(c(a$cost, b$cost), C, c(a$owner, b$owner + 3L))
    want <- menu_brute(a) + menu_brute(b)
    sg <- .menu_subgradient(p)
    r <- .menu_search(p, sg$lambda, sg$pick, grain = 1)
    if (is.finite(want)) {
      check_plan(p, r, want)
    } else {
      expect_equal(r$status, "infeasible")
    }
  }
})

test_that("the tolerance follows the plan, not the dearest states", {
  # A budget-mode objective near 1e-4 beside states costing 1e3. A tolerance
  # read off those would hide an improvement of 1e-8.
  p <- .menu_problem(c(1e-4, 1e-4 - 1e-8, 1e3, 1e-4, 1e3),
                     matrix(c(0.5, 0.4, 0, 0.5, 0)), c(1, 1, 1, 2, 2))
  r <- .menu_search(p, 0, pick = c(1L, 4L))
  expect_identical(r$pick, c(2L, 4L))
  expect_equal(r$value, 2e-4 - 1e-8, tolerance = 1e-12)
  expect_identical(r$status, "optimal")
})

test_that("a problem with no columns takes each group's cheapest state", {
  p <- .menu_problem(c(3, 1, 2, 5, 4), matrix(0, 5, 0), c(1, 1, 1, 2, 2))
  r <- .menu_search(p, numeric(0))
  expect_identical(r$pick, c(2L, 5L))
  expect_equal(r$value, 5)
  expect_equal(r$status, "optimal")
})

test_that("states in any order give the same optimum", {
  for (seed in 1:30) {
    p <- menu_random(seed, H = 6L, size = 5L)
    want <- menu_brute(p)
    if (!is.finite(want)) next
    o <- sample(length(p$cost))
    q <- .menu_problem(p$cost[o], p$C[o, , drop = FALSE], p$owner[o])
    expect_null(q$ends)
    sg <- .menu_subgradient(q)
    check_plan(q, .menu_search(q, sg$lambda, NULL, grain = 1), want)
  }
})

test_that("group minima match a loop over the runs", {
  set.seed(3)
  x <- matrix(round(stats::runif(60), 1), 20, 3)
  ends <- c(1L, 4L, 9L, 10L, 20L)
  got <- .menu_group_min(x, ends)
  runs <- split(seq_len(20), rep(seq_along(ends), diff(c(0L, ends))))
  for (k in 1:3) {
    for (g in seq_along(runs)) {
      v <- x[runs[[g]], k]
      expect_equal(got$min[g, k], min(v))
      expect_equal(got$arg[g, k], runs[[g]][which.min(v)])
      expect_equal(got$second[g, k], if (length(v) > 1) sort(v)[2] else Inf)
    }
  }
  one <- .menu_group_min(x[, 1], ends)
  expect_null(dim(one$min))
  expect_equal(one$min, got$min[, 1])
  expect_error(.menu_group_min(x, c(4L, 2L)), "increase")
  expect_error(.menu_group_min(x, c(4L, 30L)), "increase")
})

test_that("the dominance pass matches a loop over the states", {
  loop <- function(C, ord) {
    kept <- integer(0)
    vapply(ord, function(i) {
      dominated <- length(kept) &&
        any(rowSums(C[kept, , drop = FALSE] <=
                      rep(C[i, ], each = length(kept)) + 1e-15) == ncol(C))
      if (!dominated) kept <<- c(kept, i)
      !dominated
    }, logical(1))
  }
  for (seed in 1:40) {
    set.seed(seed)
    K <- sample(1:4, 1)
    n <- sample(1:60, 1)
    # Coarse values give ties and repeated rows, and a load 1e-16 above
    # another sits inside the tolerance.
    C <- matrix(round(stats::runif(n * K), 1), n, K) * sample(c(1, 100), 1)
    if (n > 2) C[2, ] <- C[1, ] + 1e-16
    ord <- sample(n, sample(n, 1))
    expect_identical(.menu_dominated(C, ord), loop(C, ord))
  }
  C <- matrix(c(1, 2, 2, 1, 2, 1, 2, 2), 4, 2)
  expect_identical(.menu_dominated(C, 1:4), c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(.menu_dominated(C, integer(0)), logical(0))
  # A kept load 2e-16 above is within the tolerance. At 100 the tolerance
  # rounds away, so a repeated row is dropped by the equality alone.
  near <- matrix(c(0.5 + 2e-16, 0.5, 0.5, 0.5), 2, 2)
  expect_identical(.menu_dominated(near, 1:2), c(TRUE, FALSE))
  expect_identical(.menu_dominated(matrix(100, 2, 2), 1:2), c(TRUE, FALSE))
  expect_error(.menu_dominated(C, c(1L, 5L)), "index rows")
  expect_error(.menu_dominated(C, c(0L, 1L)), "index rows")
  expect_error(.menu_dominated(1:4, 1L), "double matrix")
})

## A dive before best-first search
##
## Where many plans cost within a grain or two, best-first can spread over
## them without finding the cheapest, which a depth-first dive reaches and
## which then proves itself against the bound.

.dive_menu <- function(seed, H = 24, S = 60, K = 6) {
  set.seed(seed)
  owner <- rep(seq_len(H), each = S)
  base <- rep(runif(H, 100, 600), each = S)
  size <- runif(H * S, 0.2, 3)
  cost <- round(base * size / 2)
  C <- matrix(0, H * S, K)
  for (k in seq_len(K)) {
    w <- rep(runif(H, 0.5, 2), each = S) * (k %% 3 + 1)
    C[, k] <- w / (size * H * runif(1, 0.6, 1.1)) * runif(H * S, 0.9, 1.1) / 3
  }
  .menu_problem(cost, C, owner)
}

test_that("the dive proves a menu best-first leaves open", {
  p <- .dive_menu(348)
  sg <- .menu_subgradient(p)
  grain <- .menu_grain(p$cost)
  tree <- .menu_tree(p, sg$lambda, sg$pick, grain, 2000, 1e4)
  both <- .menu_component(p, sg$lambda, sg$pick, grain, 2000, 1e4)
  expect_identical(tree$status, "limit")
  expect_identical(both$status, "optimal")
  expect_lte(both$value, tree$value)
  expect_lt(both$nodes, 2000)
  ## Too few nodes to prove it: the status says so, and a proof would
  ## carry a bound equal to its value.
  short <- .menu_component(p, sg$lambda, sg$pick, grain, 200, 1e4)
  expect_false(identical(short$status, "optimal"))
  expect_lt(short$bound, short$value)
})

test_that("the dive and the search agree with enumeration on small menus", {
  for (seed in 1:30) {
    set.seed(seed)
    H <- 4
    S <- 5
    owner <- rep(seq_len(H), each = S)
    cost <- round(runif(H * S, 1, 40))
    C <- matrix(runif(H * S * 2, 0, 0.45), H * S, 2)
    p <- .menu_problem(cost, C, owner)
    plans <- as.matrix(expand.grid(lapply(p$groups, identity)))
    load_ok <- apply(plans, 1, function(j) all(colSums(C[j, ]) <= 1 + 1e-9))
    sg <- .menu_subgradient(p)
    if (!any(load_ok)) next
    best <- min(apply(plans[load_ok, , drop = FALSE], 1, function(j) sum(cost[j])))
    for (max_nodes in c(5, 50, 5000)) {
      r <- .menu_component(p, sg$lambda, sg$pick, .menu_grain(cost), max_nodes,
                           1e4)
      if (is.null(r$pick)) next
      expect_lte(r$bound, best + 1e-9)
      if (r$status == "optimal") {
        expect_equal(r$value, best)
        expect_equal(r$bound, r$value)
      }
      expect_gte(r$value, best - 1e-9)
    }
  }
})
