## One sizing solver serves both search methods, and the ordinary case is
## solved in closed form rather than by bisection.

fixture <- function() {
  set.seed(3)
  x <- rlnorm(2000, 6, 1)
  bk <- unname(quantile(x, c(0.25, 0.5, 0.75)))
  pre <- svyplan:::.strata_precompute(sort(x))
  idx <- svyplan:::.bk_to_idx(sort(x), bk)
  st <- svyplan:::.strata_stats_from_prefix(pre, idx)
  list(x = x, bk = bk, N_h = st$N_h, W_h = st$W_h, S_h = st$S_h,
       mean_h = st$mean_h)
}

test_that("the analytic size matches a hand-written inversion", {
  f <- fixture()

  for (alloc in c("proportional", "neyman", "optimal", "power")) {
    for (deff in c(1, 1.8)) {
      for (resp_rate in c(1, 0.85)) {
        a_h <- svyplan:::.alloc_weights(alloc, 0.5, f$N_h, f$S_h,
                                        rep(1, length(f$N_h)))
        m_h <- pmin(rep(2, length(f$N_h)), f$N_h)
        M_h <- f$N_h
        target_V <- (0.05 * abs(sum(f$W_h * f$mean_h)))^2

        got <- svyplan:::.strata_size_for_target(
          f$W_h, f$S_h, f$N_h, a_h, m_h, M_h, target_V, deff, resp_rate
        )

        # V(n) = A/n - B while the allocation stays proportional.
        c_h <- a_h / sum(a_h)
        ws <- f$W_h^2 * f$S_h^2
        A <- deff / resp_rate * sum(ws / c_h)
        B <- deff * sum(ws / f$N_h)
        want <- A / (target_V + B)

        lab <- paste(alloc, deff, resp_rate)
        expect_equal(got, want, tolerance = 1e-8, info = lab)

        # And the size it returns actually reaches the target.
        n_h <- svyplan:::.rna_alloc(a_h, got, m_h, M_h)
        expect_lte(
          svyplan:::.strata_variance(f$W_h, f$S_h, n_h, f$N_h, deff,
                                     resp_rate),
          target_V * (1 + 1e-8)
        )
      }
    }
  }
})

test_that("a bound-violating candidate falls back and still meets the target", {
  f <- fixture()
  a_h <- svyplan:::.alloc_weights("neyman", 0.5, f$N_h, f$S_h,
                                  rep(1, length(f$N_h)))
  # Cap the heaviest stratum so the proportional candidate breaches its upper
  # bound, while leaving the target reachable through the other strata.
  m_h <- pmin(rep(2, length(f$N_h)), f$N_h)
  M_h <- f$N_h
  j <- which.max(a_h)
  M_h[j] <- max(m_h[j] + 1, floor(f$N_h[j] / 4))
  target_V <- (0.08 * abs(sum(f$W_h * f$mean_h)))^2

  got <- svyplan:::.strata_size_for_target(
    f$W_h, f$S_h, f$N_h, a_h, m_h, M_h, target_V
  )
  expect_true(is.finite(got))

  n_h <- svyplan:::.rna_alloc(a_h, got, m_h, M_h)
  expect_lte(svyplan:::.strata_variance(f$W_h, f$S_h, n_h, f$N_h), target_V)
})

test_that("the returned size is the smallest that meets the target", {
  f <- fixture()
  a_h <- svyplan:::.alloc_weights("neyman", 0.5, f$N_h, f$S_h,
                                  rep(1, length(f$N_h)))
  m_h <- pmin(rep(2, length(f$N_h)), f$N_h)
  M_h <- f$N_h
  target_V <- (0.05 * abs(sum(f$W_h * f$mean_h)))^2

  got <- svyplan:::.strata_size_for_target(
    f$W_h, f$S_h, f$N_h, a_h, m_h, M_h, target_V
  )
  var_at <- function(n) {
    svyplan:::.strata_variance(
      f$W_h, f$S_h, svyplan:::.rna_alloc(a_h, n, m_h, M_h), f$N_h
    )
  }
  expect_lte(var_at(got), target_V * (1 + 1e-9))
  # A materially smaller size must miss it, so the answer is not merely
  # somewhere in a flat region.
  expect_gt(var_at(got * 0.99), target_V)
})

test_that("an unattainable target is reported rather than approximated", {
  f <- fixture()
  a_h <- svyplan:::.alloc_weights("neyman", 0.5, f$N_h, f$S_h,
                                  rep(1, length(f$N_h)))
  m_h <- pmin(rep(2, length(f$N_h)), f$N_h)

  # A full census drives the variance to zero, so a target is unattainable
  # only when the design cannot reach a census. Cap the takes to create that.
  capped <- pmax(m_h, floor(f$N_h / 10))
  expect_equal(
    svyplan:::.strata_size_for_target(f$W_h, f$S_h, f$N_h, a_h, m_h, capped,
                                      target_V = 1e-12),
    Inf
  )

  # A target already met by the minimum returns the minimum.
  expect_equal(
    svyplan:::.strata_size_for_target(f$W_h, f$S_h, f$N_h, a_h, m_h, f$N_h,
                                      target_V = 1e12),
    sum(m_h)
  )
})

## The objective reads the partition, so equivalent boundaries are one point.

test_that("boundaries inside one partition give one objective value", {
  set.seed(3)
  x_sort <- sort(rlnorm(500, 6, 1))
  pre <- svyplan:::.strata_precompute(x_sort)
  bk <- unname(quantile(x_sort, c(0.3, 0.6)))
  idx <- findInterval(bk, x_sort)

  # Two numerically different vectors that cut the population identically.
  alt <- c(
    (x_sort[idx[1]] + x_sort[idx[1] + 1L]) / 2,
    (x_sort[idx[2]] + x_sort[idx[2] + 1L]) / 2
  )
  expect_equal(findInterval(alt, x_sort), idx)

  obj <- function(b) {
    svyplan:::.strata_obj(x_sort, b, 200, "neyman", 0.5, rep(1, 3),
                          .pre = pre)
  }
  expect_equal(obj(bk), obj(alt), tolerance = 1e-12)
})

## LH stops on the partition, not on numeric drift inside one.

test_that("the search stops on its own rather than spending the whole budget", {
  set.seed(3)
  x <- rlnorm(2000, 6, 1)

  # Stopping is expressed on the cut indices, so a fit that has settled gives
  # the same answer whether it is allowed 25 sweeps or the full 200. A search
  # still grinding through the budget would not. Some allocations settle on a
  # fixed point and report converged, others end on a repeated partition and
  # honestly report that they did not; both must terminate.
  for (alloc in c("proportional", "neyman", "optimal", "power")) {
    short <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh",
                          alloc = alloc, max_iter = 25L)
    full <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh",
                         alloc = alloc, max_iter = 200L)

    expect_equal(full$boundaries, short$boundaries, info = alloc)
    expect_equal(full$n, short$n, info = alloc)
    expect_lte(full$cv, 0.05 + 1e-8)
    expect_equal(sum(full$strata$n), full$n, info = alloc)
  }
})

test_that("convergence survives ties and near-coincident values", {
  # Ties at plausible cut points.
  tied <- rep(c(1, 5, 20, 60, 200), times = c(40, 30, 25, 20, 15))
  res <- strata_bound(tied, n_strata = 3, n = 60, method = "lh")
  expect_equal(sum(res$strata$n), 60)
  expect_true(all(res$strata$N >= 2))

  # Adjacent unique values far closer than diff(range(x)) * 1e-6.
  near <- c(seq(1, 2, length.out = 60), 1e6 + seq(0, 1e-4, length.out = 40))
  res2 <- strata_bound(near, n_strata = 3, n = 50, method = "lh")
  expect_equal(sum(res2$strata$n), 50)
  expect_true(all(res2$strata$N >= 2))
})

test_that("LH and Kozak agree that a design meets its target", {
  set.seed(3)
  x <- rlnorm(1500, 6, 1)
  lh <- strata_bound(x, n_strata = 4, cv = 0.06, method = "lh")
  set.seed(5)
  kz <- strata_bound(x, n_strata = 4, cv = 0.06, method = "kozak",
                     n_restart = 5L, max_iter = 50L)

  expect_lte(lh$cv, 0.06 + 1e-8)
  expect_lte(kz$cv, 0.06 + 1e-8)
  # Both search the same objective, so neither should be far off the other.
  expect_lt(abs(lh$n - kz$n) / lh$n, 0.25)
})

## `converged` is a claim about the design returned, so it has to survive one
## more sweep started from exactly those boundaries. Comparing two runs of the
## same search cannot catch a wrong claim, because both take the same exit.

test_that("a result labelled converged is stable under one further sweep", {
  # One canonical coordinate sweep, written here rather than reused from the
  # package so the test cannot inherit the defect it is checking for.
  resweep <- function(x, boundaries, alloc, mode, target) {
    x_sort <- sort(x)
    L <- length(boundaries) + 1L
    pre <- svyplan:::.strata_precompute(x_sort)
    obj <- if (mode == "cv") {
      function(b) {
        svyplan:::.strata_n_for_cv(x_sort, b, target, alloc, 0.5,
                                   rep(1, L), NULL, .pre = pre)
      }
    } else {
      function(b) {
        svyplan:::.strata_obj(x_sort, b, target, alloc, 0.5,
                              rep(1, L), NULL, .pre = pre)
      }
    }
    canon <- function(b) {
      p <- findInterval(b, x_sort)
      cand <- x_sort[pmax(p, 1L)]
      ok <- p >= 1L & findInterval(cand, x_sort) == p
      b[ok] <- cand[ok]
      b
    }
    x_uniq <- sort(unique(x_sort))
    nu <- length(x_uniq)
    tol <- diff(range(x_sort)) * 1e-6
    bk <- boundaries
    for (h in seq_len(L - 1L)) {
      lo <- if (h == 1L) x_uniq[1L] + tol else bk[h - 1L] + tol
      hi <- if (h == L - 1L) x_uniq[nu] - tol else bk[h + 1L] - tol
      if (lo >= hi) next
      o <- suppressWarnings(optimize(
        function(b) {
          bt <- bk
          bt[h] <- b
          obj(bt)
        },
        interval = c(lo, hi)
      ))
      cand <- bk
      cand[h] <- o$minimum
      cur <- obj(bk)
      if (!is.finite(cur) || obj(cand) <= cur) bk <- cand
    }
    findInterval(canon(bk), x_sort)
  }

  set.seed(99)
  checked <- 0L
  for (i in seq_len(6L)) {
    n <- sample(c(300, 500, 800), 1L)
    x <- rlnorm(n, 6, 1)
    for (alloc in c("proportional", "neyman", "optimal", "power")) {
      mode <- if (i %% 2L == 0L) "cv" else "n"
      target <- if (mode == "cv") 0.05 else max(50, round(n / 10))
      res <- if (mode == "cv") {
        strata_bound(x, n_strata = 4, cv = target, method = "lh",
                     alloc = alloc)
      } else {
        strata_bound(x, n_strata = 4, n = target, method = "lh",
                     alloc = alloc)
      }
      if (!isTRUE(res$converged)) next
      checked <- checked + 1L
      expect_equal(
        resweep(x, res$boundaries, alloc, mode, target),
        findInterval(res$boundaries, sort(x)),
        info = paste(alloc, mode, n)
      )
    }
  }
  # The assertion above is only meaningful if convergence is actually reported.
  expect_gt(checked, 10L)
})

test_that("a converged result is the best partition the search found", {
  set.seed(1)
  x <- rlnorm(500, 6, 1)
  # The case that used to report convergence at its untouched starting
  # quantiles, which one further sweep left immediately.
  res <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh",
                      alloc = "proportional")
  start <- unname(quantile(sort(x), c(0.25, 0.5, 0.75)))
  expect_true(res$converged)
  expect_false(identical(
    findInterval(res$boundaries, sort(x)),
    findInterval(start, sort(x))
  ))
  expect_true(is.finite(res$n))
})
