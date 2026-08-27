test_that("bounded allocation respects random feasible totals and bounds", {
  set.seed(20260827)
  for (i in seq_len(100L)) {
    L <- sample(2:12, 1L)
    m_h <- sample(0:5, L, replace = TRUE)
    M_h <- m_h + sample(0:30, L, replace = TRUE)
    a_h <- rexp(L)
    a_h[sample(c(TRUE, FALSE), L, replace = TRUE, prob = c(0.2, 0.8))] <- 0
    n_total <- runif(1L, sum(m_h), sum(M_h))

    n_h <- svyplan:::.rna_alloc(a_h, n_total, m_h, M_h)

    expect_equal(sum(n_h), n_total, tolerance = 1e-8)
    expect_true(all(n_h >= m_h - 1e-8))
    expect_true(all(n_h <= M_h + 1e-8))
  }
})

test_that("bounded allocation handles simultaneous lower and upper pressure", {
  n_h <- svyplan:::.rna_alloc(
    a_h = c(1295187.5, 631803.5, 61558864.5, 0),
    n_total = 100,
    m_h = c(2, 2, 2, 1),
    M_h = c(354, 58, 87, 1)
  )

  expect_equal(sum(n_h), 100, tolerance = 1e-8)
  expect_equal(n_h[3:4], c(87, 1), tolerance = 1e-8)
  expect_true(all(n_h >= c(2, 2, 2, 1)))
  expect_true(all(n_h <= c(354, 58, 87, 1)))
})

test_that("bounded allocation rejects infeasible totals", {
  expect_error(
    svyplan:::.rna_alloc(c(1, 2), 2.9, c(2, 1), c(5, 5)),
    "feasible range is 3 to 10",
    class = "svyplan_alloc_infeasible"
  )
  expect_error(
    svyplan:::.rna_alloc(c(1, 2), 10.1, c(2, 1), c(5, 5)),
    "feasible range is 3 to 10",
    class = "svyplan_alloc_infeasible"
  )
})

## An infeasible total must be refused rather than met by collapsing a stratum
## to a single unit, which lowers its own minimum from two to one.

test_that("a boundary search refuses to reach a total through a singleton stratum", {
  set.seed(1)
  x <- c(rlnorm(160), rlnorm(40, 4, 0.2))
  threshold <- as.numeric(quantile(x, 0.8))

  for (method in c("lh", "kozak")) {
    args <- list(
      x = x, n_strata = 3, take_all_above = threshold,
      method = method, max_iter = 50L
    )
    if (method == "kozak") args$n_restart <- 5L

    for (n in 42:43) {
      expect_error(
        do.call(strata_bound, c(args, list(n = n))),
        "infeasible",
        info = paste(method, n)
      )
    }

    res <- do.call(strata_bound, c(args, list(n = 44)))
    expect_equal(sum(res$strata$n), 44, info = method)
    expect_true(all(res$strata$N >= 2), info = method)
  }
})

test_that("no search method returns a stratum it cannot sample", {
  set.seed(1)
  x <- c(rlnorm(160), rlnorm(40, 4, 0.2))

  for (method in c("lh", "kozak")) {
    args <- list(x = x, n_strata = 4, method = method, max_iter = 50L)
    if (method == "kozak") args$n_restart <- 10L

    expect_error(
      do.call(strata_bound, c(args, list(n = 7))),
      "infeasible",
      info = method
    )

    for (target in c(0.3, 0.1)) {
      res <- do.call(strata_bound, c(args, list(cv = target)))
      expect_true(all(res$strata$N >= 2), info = paste(method, target))
      expect_true(all(res$strata$n >= 2), info = paste(method, target))
    }
  }
})

test_that("a take-all stratum of one unit stays legal", {
  expect_false(svyplan:::.strata_degenerate(c(5, 1), take_all_idx = 2L))
  expect_true(svyplan:::.strata_degenerate(c(5, 1), take_all_idx = NULL))
  expect_true(svyplan:::.strata_degenerate(c(0, 5), take_all_idx = 2L))
})

## Ties can leave no valid stratification at all. Filtering candidates is not
## enough there, because the search has nothing to fall back to.

test_that("data admitting no valid stratification is refused, not approximated", {
  x <- c(rep(1, 10), 2, rep(3, 10), 4)

  for (method in c("lh", "kozak", "cumrootf", "geo")) {
    expect_error(
      suppressWarnings(
        strata_bound(x, n_strata = 4, n = 8, method = method)
      ),
      "at least two population units",
      info = method
    )
  }

  expect_error(
    strata_bound(x, n_strata = 3, n = 8, method = "lh"),
    "at least two population units"
  )
  expect_error(
    strata_bound(x, n_strata = 4, cv = 0.05, method = "lh"),
    "at least two population units"
  )
  expect_error(
    suppressWarnings(strata_bound(x, n_strata = 4, method = "cumrootf")),
    "at least two population units"
  )

  res <- strata_bound(x, n_strata = 2, n = 8, method = "lh")
  expect_equal(res$strata$N, c(11L, 11L))
  expect_equal(sum(res$strata$n), 8)
})

test_that("the fixed-n floor is refused before any boundary is searched", {
  set.seed(7)
  x <- c(runif(60, 1, 10), 100)

  expect_error(
    strata_bound(x, n_strata = 3, n = 4, take_all_above = 50, method = "lh"),
    "minimum feasible is 5"
  )
  expect_error(
    strata_bound(x, n_strata = 4, n = 7, method = "lh"),
    "minimum feasible is 8"
  )

  res <- strata_bound(x, n_strata = 3, n = 5, take_all_above = 50,
                      method = "lh")
  expect_equal(sum(res$strata$n), 5)
})

test_that("an explicit one-unit take-all stratum is accepted", {
  set.seed(7)
  x <- c(runif(60, 1, 10), 100)

  res <- strata_bound(x, n_strata = 3, n = 20, take_all_above = 50,
                      method = "lh")

  expect_equal(res$strata$N[3], 1L)
  expect_equal(res$strata$n[3], 1)
  expect_true(res$strata$take_all[3])
  expect_true(all(res$strata$N[1:2] >= 2))
  expect_equal(sum(res$strata$n), 20)
})

test_that("cv mode never returns a stratum it cannot sample", {
  set.seed(1)
  x <- c(rlnorm(160), rlnorm(40, 4, 0.2))
  threshold <- as.numeric(quantile(x, 0.8))

  for (target in c(0.2, 0.1, 0.02)) {
    for (method in c("lh", "kozak")) {
      res <- strata_bound(x, n_strata = 3, cv = target,
                          take_all_above = threshold, method = method,
                          max_iter = 40L)
      sampled <- !res$strata$take_all
      expect_true(all(res$strata$N[sampled] >= 2),
                  info = paste(method, target))
      expect_true(all(res$strata$n[sampled] >= 2),
                  info = paste(method, target))
    }
  }
})


## The unconstrained ratio allocation is the bounded solution whenever it
## already satisfies every bound, so taking it directly is exact.

test_that("the unconstrained fast path equals an independent KKT oracle", {
  # Solve sum(clamp(lambda * a, m, M)) = total directly. This states the
  # bounded ratio rule rather than reproducing the breakpoint implementation.
  kkt <- function(a, total, m, M) {
    S <- function(l) sum(pmin(pmax(l * a, m), M))
    hi <- 1
    while (S(hi) < total) hi <- hi * 2
    lambda <- uniroot(function(l) S(l) - total, c(0, hi), tol = 1e-14)$root
    pmin(pmax(lambda * a, m), M)
  }

  set.seed(4242)
  for (i in seq_len(200L)) {
    L <- sample(2:6, 1L)
    a_h <- runif(L, 0.1, 5)
    N_h <- sample(20:200, L, replace = TRUE)
    m_h <- rep(2, L)
    M_h <- N_h
    total <- runif(1, sum(m_h), sum(M_h))
    if (total > sum(M_h)) next

    expect_equal(
      svyplan:::.rna_alloc(a_h, total, m_h, M_h),
      kkt(a_h, total, m_h, M_h),
      tolerance = 1e-8
    )
  }
})

test_that("the fast path is declined exactly when a bound would bind", {
  taken <- function(a_h, total, m_h, M_h) {
    raw <- total * a_h / sum(a_h)
    all(raw >= m_h) && all(raw <= M_h)
  }

  # Interior: the raw ratio solution is strictly inside every bound.
  a <- c(4, 3, 2, 1)
  m <- c(2, 2, 2, 2)
  M <- c(500, 300, 150, 60)
  expect_true(taken(a, 100, m, M))
  expect_equal(svyplan:::.rna_alloc(a, 100, m, M), 100 * a / sum(a))

  # Exactly on a lower bound stays on the fast path and returns the ratio.
  a2 <- c(48, 1)
  m2 <- c(2, 2)
  M2 <- c(500, 300)
  expect_true(taken(a2, 98, m2, M2))
  expect_equal(svyplan:::.rna_alloc(a2, 98, m2, M2), c(96, 2))

  # Exactly on an upper bound likewise.
  a3 <- c(1, 1)
  expect_true(taken(a3, 100, c(2, 2), c(50, 300)))
  expect_equal(svyplan:::.rna_alloc(a3, 100, c(2, 2), c(50, 300)), c(50, 50))

  # A zero weight whose lower bound is zero may use the fast path.
  expect_true(taken(c(3, 1, 0), 40, c(0, 0, 0), c(50, 50, 50)))
  expect_equal(svyplan:::.rna_alloc(c(3, 1, 0), 40, c(0, 0, 0), c(50, 50, 50)),
               c(30, 10, 0))

  # A zero weight with a positive minimum, a take-all component included,
  # must decline and take the general path.
  expect_false(taken(c(3, 1, 0), 60, c(2, 2, 20), c(50, 50, 20)))
  res <- svyplan:::.rna_alloc(c(3, 1, 0), 60, c(2, 2, 20), c(50, 50, 20))
  expect_equal(sum(res), 60)
  expect_equal(res[3], 20)

  # Simultaneous lower and upper pressure must decline.
  expect_false(taken(c(1295187.5, 631803.5, 61558864.5, 0), 100,
                     c(2, 2, 2, 1), c(354, 58, 87, 1)))
})

test_that("every route out of the allocator satisfies the post-condition", {
  expect_error(
    svyplan:::.rna_check(c(1, 1), 5, c(0, 0), c(10, 10), 1e-10),
    "violated its total or bounds"
  )
  expect_error(
    svyplan:::.rna_check(c(1, 9), 10, c(2, 2), c(10, 10), 1e-10),
    "violated its total or bounds"
  )
  expect_equal(
    svyplan:::.rna_check(c(4, 6), 10, c(2, 2), c(10, 10), 1e-10),
    c(4, 6)
  )
})
