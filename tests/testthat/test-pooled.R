## T1. The two boundary identities that fix the estimand

test_that("a fresh sample each occasion pools to the level over occasions", {
  for (Tocc in c(2, 3, 8, 20)) {
    for (N in c(Inf, 1e5)) {
      lvl <- prec_mean(var = 100, n = 1000, N = N)$se^2
      pooled <- prec_pooled(var = 100, n = 1000, occasions = Tocc, N = N)$se^2
      expect_equal(pooled, lvl / Tocc)
    }
  }
})

test_that("a correlation without overlap pools as if independent", {
  a <- prec_pooled(var = 100, n = 1000, occasions = 6, N = 1e5)
  b <- prec_pooled(var = 100, n = 1000, occasions = 6, N = 1e5,
                   overlap = 0, overlap_cor = 0.9)
  expect_equal(a$se, b$se)
})

test_that("a full panel at unit correlation pools to a single occasion", {
  for (Tocc in c(2, 5, 12)) {
    for (N in c(Inf, 1e5)) {
      lvl <- prec_mean(var = 100, n = 1000, N = N)$se^2
      pooled <- prec_pooled(var = 100, n = 1000, occasions = Tocc, N = N,
                            overlap = 1, overlap_cor = 1)$se^2
      expect_equal(pooled, lvl)
    }
  }
})

## T2. The identity that ties the pooled engine to the shipped change engine

test_that("V_pooled = V_level - V_change/4 at two occasions", {
  grid <- expand.grid(
    n = c(100, 1000, 30000),
    ov = c(0, 0.25, 0.5, 1),
    rho = c(0, 0.3, 0.9, 1),
    N = c(Inf, 1e5)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    lvl <- prec_mean(var = 100, n = g$n, N = g$N)$se^2
    chg <- prec_change(var = 100, n = g$n, N = g$N,
                       overlap = g$ov, overlap_cor = g$rho)$se^2
    pool <- prec_pooled(var = 100, n = g$n, occasions = 2, N = g$N,
                        overlap = g$ov, overlap_cor = g$rho)$se^2
    expect_equal(pool, lvl - chg / 4,
                 info = sprintf("n=%g ov=%g rho=%g N=%g",
                                g$n, g$ov, g$rho, g$N))
  }
})

test_that("the identity survives deff and a response rate", {
  for (de in c(1, 1.4)) {
    for (rr in c(1, 0.7)) {
      lvl <- prec_mean(var = 100, n = 800, N = 1e5, deff = de,
                       resp_rate = rr)$se^2
      chg <- prec_change(var = 100, n = 800, N = 1e5, deff = de,
                         resp_rate = rr, overlap = 0.5,
                         overlap_cor = 0.8)$se^2
      pool <- prec_pooled(var = 100, n = 800, occasions = 2, N = 1e5,
                          deff = de, resp_rate = rr, overlap = 0.5,
                          overlap_cor = 0.8)$se^2
      expect_equal(pool, lvl - chg / 4)
    }
  }
})

## T3. The direction of the trade-off

test_that("overlap inflates a pooled estimate while it induces a positive covariance", {
  # n/N = 0.01 here, so every overlap tried is above the sampling fraction.
  se <- vapply(c(0.02, 0.25, 0.5, 0.75, 1), function(o) {
    prec_pooled(var = 100, n = 1000, occasions = 8, N = 1e5,
                overlap = o, cor_decay = 0.8)$se
  }, numeric(1))
  expect_true(all(diff(se) > 0))
})

test_that("without an FPC any positive overlap inflates a pooled estimate", {
  se <- vapply(c(0, 0.25, 0.5, 0.75, 1), function(o) {
    prec_pooled(var = 100, n = 1000, occasions = 8, overlap = o,
                cor_decay = 0.8)$se
  }, numeric(1))
  expect_true(all(diff(se) > 0))
})

test_that("an overlap below the sampling fraction reverses both directions", {
  # The covariance is rho * S^2 * (overlap/n - 1/N), so it is negative when
  # the occasions share fewer units than independent draws would give them.
  # The kernel stays a valid covariance and the arms swap roles.
  S2 <- 100; N <- 1000; n <- 500; rho <- 0.8
  indep_p <- prec_pooled(var = S2, n = n, occasions = 2, N = N)$se^2
  indep_c <- prec_change(var = S2, n = n, N = N)$se^2
  below_p <- prec_pooled(var = S2, n = n, occasions = 2, N = N,
                         overlap = 0.25, overlap_cor = rho)$se^2
  below_c <- prec_change(var = S2, n = n, N = N,
                         overlap = 0.25, overlap_cor = rho)$se^2
  expect_lt(below_p, indep_p)
  expect_gt(below_c, indep_c)
  # At exactly the sampling fraction the covariance vanishes and both
  # return to the independent values.
  at_p <- prec_pooled(var = S2, n = n, occasions = 2, N = N,
                      overlap = n / N, overlap_cor = rho)$se^2
  at_c <- prec_change(var = S2, n = n, N = N,
                      overlap = n / N, overlap_cor = rho)$se^2
  expect_equal(at_p, indep_p)
  expect_equal(at_c, indep_c)
  # Above it the ordinary direction holds.
  above_p <- prec_pooled(var = S2, n = n, occasions = 2, N = N,
                         overlap = 0.75, overlap_cor = rho)$se^2
  above_c <- prec_change(var = S2, n = n, N = N,
                         overlap = 0.75, overlap_cor = rho)$se^2
  expect_gt(above_p, indep_p)
  expect_lt(above_c, indep_c)
})

test_that("the level is flat in the schedule while the arms move apart", {
  lvl <- prec_mean(var = 100, n = 1000, N = 1e5)$se
  arms <- lapply(c(2, 4, 6, 8), function(L) {
    ov <- as.numeric(design_overlap(as.character(L), max_lag = 7))
    ov <- c(ov, rep(0, 7 - length(ov)))
    list(
      change = prec_change(var = 100, n = 1000, N = 1e5,
                           overlap = ov[1], overlap_cor = 0.8)$se,
      pooled = prec_pooled(var = 100, n = 1000, occasions = 8, N = 1e5,
                           overlap = ov, cor_decay = 0.8)$se
    )
  })
  chg <- vapply(arms, function(a) a$change, numeric(1))
  pool <- vapply(arms, function(a) a$pooled, numeric(1))
  expect_true(all(diff(chg) < 0))
  expect_true(all(diff(pool) > 0))
  # The reference line does not move with the schedule the arms are read from.
  expect_equal(lvl, prec_mean(var = 100, n = 1000, N = 1e5)$se)
})

## T4. The per-lag zero rule

test_that("a lag beyond the life contributes exactly zero, not a population term", {
  # Life 2 over 4 occasions shares nothing at lags 2 and 3. Padding the
  # profile with zeros must agree with stating only the lags that share.
  a <- prec_pooled(var = 100, n = 1000, occasions = 4, N = 1e5,
                   overlap = c(0.5, 0, 0), overlap_cor = 0.8)
  b <- prec_pooled(var = 100, n = 1000, occasions = 4, N = 1e5,
                   overlap = c(0.5, 0, 0), overlap_cor = c(0.8, 0, 0))
  expect_equal(a$se, b$se)
})

test_that("zero overlap is not the limit of a vanishing one", {
  # prec_change() drops the population term at overlap = 0, so the pooled
  # engine reading the same kernel must have the same discontinuity.
  at_zero <- prec_pooled(var = 100, n = 1000, occasions = 4, N = 1e5,
                         overlap = 0, overlap_cor = 0.8)$se^2
  near_zero <- prec_pooled(var = 100, n = 1000, occasions = 4, N = 1e5,
                           overlap = 1e-12, overlap_cor = 0.8)$se^2
  expect_true(near_zero < at_zero)
  # The dropped population term is one 2 (T - m) rho S^2 / N per lag, over T^2.
  expect_equal(at_zero - near_zero, 2 * sum(4 - (1:3)) * 0.8 * 100 / (1e5 * 16),
               tolerance = 1e-6)
})

## T5. The PSD guard, on the assembled matrix and not on the reported number

test_that("an AR(1) correlation is rejected once the sampling fraction breaks it", {
  N <- 1e5
  ov <- as.numeric(design_overlap("4", max_lag = 7))
  ok <- prec_pooled(var = 100, n = 0.5 * N, occasions = 8, N = N,
                    overlap = ov, cor_decay = 0.8)
  expect_true(ok$se > 0)
  expect_error(
    prec_pooled(var = 100, n = 0.67 * N, occasions = 8, N = N,
                overlap = ov, cor_decay = 0.8),
    "not a valid covariance"
  )
})

test_that("the guard names the sampling fraction and the lag that turns negative", {
  N <- 1e5
  ov <- as.numeric(design_overlap("4", max_lag = 7))
  expect_error(
    prec_pooled(var = 100, n = 0.67 * N, occasions = 8, N = N,
                overlap = ov, cor_decay = 0.8),
    "sampling fraction of 0.67"
  )
  expect_error(
    prec_pooled(var = 100, n = 0.67 * N, occasions = 8, N = N,
                overlap = ov, cor_decay = 0.8),
    "lag 2 is negative"
  )
})

test_that("the guard catches a kernel whose pooled variance is still positive", {
  # This is the case a .safe_variance() check on the result would pass: the
  # reported quantity is positive while the matrix is not a covariance.
  N <- 1e5
  ov <- as.numeric(design_overlap("4", max_lag = 7))
  n_net <- 0.67 * N
  K <- svyplan:::.pooled_kernel(100, n_net, N, 8, ov, 0.8^seq_len(7))
  expect_true(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) < 0)
  expect_true(sum(K) / 64 > 0)
  expect_error(
    prec_pooled(var = 100, n = n_net, occasions = 8, N = N,
                overlap = ov, cor_decay = 0.8),
    "not a valid covariance"
  )
})

test_that("a supplied correlation vector that is not PSD is rejected", {
  comb <- rep(c(1, 0), length.out = 7)
  expect_error(
    prec_pooled(var = 100, n = 1000, occasions = 8, N = 1e5,
                overlap = 1, overlap_cor = comb),
    "not a valid covariance"
  )
})

test_that("an infinite population never trips the guard", {
  # Without an FPC no lag covariance can turn negative, so every profile in
  # the contract's range must be accepted.
  for (o in c(0, 0.3, 1)) {
    for (r in c(0, 0.5, 1)) {
      expect_silent(
        prec_pooled(var = 100, n = 1000, occasions = 10, overlap = o,
                    overlap_cor = r)
      )
    }
  }
})

test_that("the sampling fraction is swept across the range the contract claims", {
  N <- 1e5
  for (f in c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95)) {
    res <- try(
      prec_pooled(var = 100, n = f * N, occasions = 4, N = N,
                  overlap = 0.75, cor_decay = 0.9),
      silent = TRUE
    )
    if (!inherits(res, "try-error")) expect_true(res$se >= 0)
  }
  # A full panel is safe at every fraction, its covariance never changing
  # sign, so the sweep must not reject it anywhere.
  for (f in c(0.01, 0.5, 0.9, 0.99)) {
    expect_silent(
      prec_pooled(var = 100, n = f * N, occasions = 6, N = N,
                  overlap = 1, overlap_cor = 1)
    )
  }
})

## T6. Issued versus respondent overlap

test_that("a design_overlap object is accepted at full response", {
  a <- prec_pooled(var = 100, n = 500, occasions = 8,
                   overlap = design_overlap("4", max_lag = 7), cor_decay = 0.8)
  b <- prec_pooled(var = 100, n = 500, occasions = 8,
                   overlap = c(0.75, 0.5, 0.25, 0, 0, 0, 0), cor_decay = 0.8)
  expect_equal(a$se, b$se)
})

test_that("a design_overlap object is refused below full response", {
  expect_error(
    prec_pooled(var = 100, n = 500, occasions = 8, resp_rate = 0.8,
                overlap = design_overlap("4", max_lag = 7)),
    "issued samples"
  )
  expect_error(
    n_pooled(var = 100, moe = 1, occasions = 8, resp_rate = 0.8,
             overlap = design_overlap("4", max_lag = 7)),
    "resp_rate = 1"
  )
})

test_that("a bare number is respondent overlap at any response rate", {
  expect_silent(
    prec_pooled(var = 100, n = 500, occasions = 4, resp_rate = 0.8,
                overlap = 0.75, overlap_cor = 0.5)
  )
})

test_that("an object shorter than the horizon is padded, not rejected", {
  # design_overlap() defaults to max_lag = life - 1, and beyond a life
  # nothing is shared.
  a <- prec_pooled(var = 100, n = 500, occasions = 8,
                   overlap = design_overlap("4"), cor_decay = 0.8)
  b <- prec_pooled(var = 100, n = 500, occasions = 8,
                   overlap = design_overlap("4", max_lag = 7), cor_decay = 0.8)
  expect_equal(a$se, b$se)
})

## T7. The lag correlation model

test_that("cor_decay gives rho^m", {
  a <- prec_pooled(var = 100, n = 500, occasions = 4, overlap = 1,
                   cor_decay = 0.8)
  b <- prec_pooled(var = 100, n = 500, occasions = 4, overlap = 1,
                   overlap_cor = 0.8^(1:3))
  expect_equal(a$se, b$se)
})

test_that("overlap_cor and cor_decay are alternatives", {
  expect_error(
    prec_pooled(var = 100, n = 500, occasions = 4, overlap_cor = 0.5,
                cor_decay = 0.8),
    "exactly one"
  )
})

test_that("a profile of the wrong length is refused", {
  expect_error(
    prec_pooled(var = 100, n = 500, occasions = 4, overlap = c(0.5, 0.4)),
    "one per lag"
  )
  expect_error(
    prec_pooled(var = 100, n = 500, occasions = 4, overlap = 0.5,
                overlap_cor = c(0.5, 0.4)),
    "one per lag"
  )
})

test_that("out-of-range profiles are refused", {
  expect_error(prec_pooled(var = 100, n = 500, occasions = 4, overlap = 1.5),
               "\\[0, 1\\]")
  expect_error(prec_pooled(var = 100, n = 500, occasions = 4,
                           overlap = c(0.5, -0.1, 0.2)), "\\[0, 1\\]")
  expect_error(prec_pooled(var = 100, n = 500, occasions = 4, overlap = 1,
                           cor_decay = 1.2), "\\[0, 1\\]")
})

## T8. Sizing

test_that("n_pooled inverts prec_pooled exactly", {
  grid <- expand.grid(
    Tocc = c(2, 4, 12),
    ov = c(0, 0.5, 1),
    N = c(Inf, 1e6),
    deff = c(1, 1.6),
    rr = c(1, 0.75)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    res <- n_pooled(var = 100, moe = 1, occasions = g$Tocc, N = g$N,
                    deff = g$deff, resp_rate = g$rr, overlap = g$ov,
                    cor_decay = 0.8)
    back <- prec_pooled(res)
    expect_equal(back$moe, 1,
                 info = sprintf("T=%g ov=%g N=%g deff=%g rr=%g",
                                g$Tocc, g$ov, g$N, g$deff, g$rr))
    expect_equal(back$se, res$se)
  }
})

test_that("the size is gross, carrying deff and the response inflation", {
  base <- n_pooled(var = 100, moe = 1, occasions = 4)
  expect_equal(
    n_pooled(var = 100, moe = 1, occasions = 4, resp_rate = 0.5)$n,
    base$n * 2
  )
  expect_equal(
    n_pooled(var = 100, moe = 1, occasions = 4, deff = 2)$n,
    base$n * 2
  )
})

test_that("overlap raises the size a pooled target needs", {
  n <- vapply(c(0, 0.5, 1), function(o) {
    n_pooled(var = 100, moe = 1, occasions = 8, overlap = o,
             cor_decay = 0.8)$n
  }, numeric(1))
  expect_true(all(diff(n) > 0))
})

test_that("no pooled target is unattainable for want of precision", {
  # B <= 0 under the nonnegative correlation contract, so however far the
  # occasions overlap a finite size reaches any target that fits in N.
  for (o in c(0, 0.5, 1)) {
    res <- n_pooled(var = 100, moe = 0.05, occasions = 12, N = 1e9,
                    overlap = o, overlap_cor = 1)
    expect_true(is.finite(res$n) && res$n > 0)
  }
})

test_that("every target is attainable at full response with no overlap", {
  # The same boundary n_change() documents: disjoint occasions of one
  # population pool to zero variance at a census, so the size rises towards
  # N and never past it.
  res <- n_pooled(var = 100, moe = 0.01, occasions = 4, N = 500)
  expect_true(res$n < 500)
  expect_equal(prec_pooled(res)$moe, 0.01)
})

test_that("a size beyond the population is reported as unattainable", {
  # A response rate below 1 removes the guarantee, the size released being
  # larger than the sample carrying the precision.
  expect_error(
    n_pooled(var = 100, moe = 0.01, occasions = 4, N = 500, resp_rate = 0.5),
    "exceeds the population"
  )
})

test_that("relative targets round-trip on both scales", {
  a <- n_pooled(var = 100, mu = 20, rmoe = 0.05, occasions = 4)
  expect_equal(prec_pooled(a)$rmoe, 0.05)
  b <- n_pooled(var = 100, mu = 20, cv = 0.02, occasions = 4)
  expect_equal(prec_pooled(b)$cv, 0.02)
  d <- n_pooled(p = 0.3, rmoe = 0.05, occasions = 4)
  expect_equal(prec_pooled(d)$rmoe, 0.05)
})

test_that("cv and rmoe need a level", {
  expect_error(n_pooled(var = 100, cv = 0.05, occasions = 4), "'mu' is required")
  expect_error(n_pooled(var = 100, rmoe = 0.05, occasions = 4), "mu")
})

## T9. Scales and their guards

test_that("the proportion scale carries the N/(N-1) adjustment", {
  p <- 0.3
  N <- 1000
  a <- prec_pooled(p = p, n = 100, occasions = 2, N = N, overlap = 0)
  b <- prec_pooled(var = p * (1 - p) * N / (N - 1), n = 100, occasions = 2,
                   N = N, overlap = 0)
  expect_equal(a$se, b$se)
})

test_that("one occasion of a pooled proportion agrees with prec_prop", {
  a <- prec_pooled(p = 0.3, n = 500, occasions = 2, N = 1e4)$se
  b <- prec_prop(p = 0.3, n = 500, N = 1e4)$se
  expect_equal(a, b / sqrt(2))
})

test_that("the two scales are exclusive and mu is refused beside p", {
  expect_error(prec_pooled(var = 100, p = 0.3, n = 500, occasions = 4),
               "exactly one")
  expect_error(prec_pooled(n = 500, occasions = 4), "exactly one")
  expect_error(prec_pooled(p = 0.3, mu = 0.3, n = 500, occasions = 4),
               "do not supply it")
})

test_that("sd is an alternative spelling of var", {
  expect_equal(
    prec_pooled(sd = 10, n = 500, occasions = 4)$se,
    prec_pooled(var = 100, n = 500, occasions = 4)$se
  )
})

## T10. occasions

test_that("occasions must be a whole number of at least 2", {
  for (bad in list(1, 0, -3, 2.5, NA_real_, Inf, c(2, 3), "4")) {
    expect_error(prec_pooled(var = 100, n = 500, occasions = bad),
                 "whole number")
  }
})

## T11. Methods

test_that("print reports the estimand, the horizon and the profile", {
  res <- prec_pooled(var = 100, n = 500, occasions = 4, overlap = 0.75,
                     cor_decay = 0.8)
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "pooled estimate")
  expect_match(out, "4 occasions")
  expect_match(out, "shared out to lag 3")
  flat <- prec_pooled(var = 100, n = 500, occasions = 4)
  expect_match(paste(capture.output(print(flat)), collapse = "\n"),
               "No between-occasion covariance")
})

test_that("print reports the net size when response is below full", {
  res <- n_pooled(var = 100, moe = 1, occasions = 4, resp_rate = 0.5)
  expect_match(paste(capture.output(print(res)), collapse = "\n"), "net:")
})

test_that("confint centres on the level and refuses without one", {
  res <- n_pooled(var = 100, mu = 20, moe = 1, occasions = 4)
  ci <- confint(res)
  expect_equal(unname(ci[1, 1]), 20 - res$moe)
  expect_equal(unname(ci[1, 2]), 20 + res$moe)
  bare <- n_pooled(var = 100, moe = 1, occasions = 4)
  expect_error(confint(bare), "'mu'")
  expect_equal(unname(confint(prec_pooled(res))[1, 1]), 20 - res$moe)
})

test_that("predict varies a pooled target over a grid", {
  res <- n_pooled(var = 100, moe = 1, occasions = 4, overlap = 0.5,
                  overlap_cor = 0.8)
  out <- predict(res, data.frame(moe = c(0.5, 1, 2)))
  expect_equal(nrow(out), 3L)
  expect_true(all(diff(out$n) < 0))
  expect_equal(out$moe, c(0.5, 1, 2))
})

test_that("predict varies occasions only where the profile is flat", {
  flat <- n_pooled(var = 100, moe = 1, occasions = 4, overlap = 0.5,
                   overlap_cor = 0.8)
  out <- predict(flat, data.frame(occasions = c(2, 4, 8)))
  expect_equal(nrow(out), 3L)
  shaped <- n_pooled(var = 100, moe = 1, occasions = 4,
                     overlap = c(0.75, 0.5, 0.25), overlap_cor = 0.8)
  expect_error(predict(shaped, data.frame(occasions = c(2, 4))),
               "unknown parameter")
})

test_that("predict on a pooled precision object is refused", {
  res <- prec_pooled(var = 100, n = 500, occasions = 4)
  expect_error(predict(res, data.frame(n = c(100, 200))), "not supported")
})

test_that("the round trip refuses a result of the wrong type", {
  expect_error(prec_pooled(n_change(var = 100, moe = 2)), "type 'pooled'")
  expect_error(n_pooled(prec_change(var = 100, n = 500)), "type 'pooled'")
})

test_that("unused arguments are rejected", {
  expect_error(prec_pooled(var = 100, n = 500, occasions = 4, nonsense = 1))
  expect_error(n_pooled(var = 100, moe = 1, occasions = 4, nonsense = 1))
})

test_that("a svyplan object supplies design defaults", {
  plan <- svyplan(deff = 2, resp_rate = 0.8)
  a <- prec_pooled(var = 100, n = 500, occasions = 4, plan = plan)
  b <- prec_pooled(var = 100, n = 500, occasions = 4, deff = 2,
                   resp_rate = 0.8)
  expect_equal(a$se, b$se)
})

## T12. The stored profile is the resolved one

test_that("params carry one overlap and one correlation per lag", {
  res <- prec_pooled(var = 100, n = 500, occasions = 5, overlap = 0.5,
                     cor_decay = 0.8)
  expect_length(res$params$overlap, 4L)
  expect_length(res$params$overlap_cor, 4L)
  expect_equal(res$params$overlap, rep(0.5, 4))
  expect_equal(res$params$overlap_cor, 0.8^(1:4))
})

## T13. The PSD guard does not depend on the outcome's units

test_that("the PSD verdict is invariant to rescaling the outcome", {
  # Positive semidefiniteness is a property of the correlation structure.
  # A tolerance with an absolute floor rejects a design at one scale and
  # accepts the same design rescaled, which is the defect this pins.
  N <- 1e5
  ov <- as.numeric(design_overlap("4", max_lag = 7))
  for (v in c(1e6, 100, 1, 1e-4, 1e-8)) {
    expect_error(
      prec_pooled(var = v, n = 0.67 * N, occasions = 8, N = N,
                  overlap = ov, cor_decay = 0.8),
      "not a valid covariance",
      info = sprintf("var = %g", v)
    )
  }
})

test_that("a valid design stays accepted at every scale, and se scales with sd", {
  N <- 1e5
  ov <- as.numeric(design_overlap("4", max_lag = 7))
  base <- prec_pooled(var = 1, n = 0.5 * N, occasions = 8, N = N,
                      overlap = ov, cor_decay = 0.8)$se
  for (v in c(1e6, 100, 1e-4, 1e-8)) {
    got <- prec_pooled(var = v, n = 0.5 * N, occasions = 8, N = N,
                       overlap = ov, cor_decay = 0.8)$se
    expect_equal(got, base * sqrt(v))
  }
})

test_that("a census at full overlap gives the zero matrix and no variance", {
  # The scale is zero there, so the tolerance must pass on the equality
  # rather than dividing by it.
  res <- prec_pooled(var = 100, n = 1000, occasions = 4, N = 1000,
                     overlap = 1, overlap_cor = 1)
  expect_equal(res$se, 0)
})

## T14. Issued provenance survives the round trip and the grid

test_that("a stored issued profile is refused when a round trip lowers response", {
  x <- n_pooled(var = 100, moe = 1, occasions = 8,
                overlap = design_overlap("4", max_lag = 7), cor_decay = 0.8)
  expect_equal(x$params$overlap_basis, "issued")
  expect_error(prec_pooled(x, resp_rate = 0.8), "issued samples")
  expect_error(predict(x, data.frame(resp_rate = 0.8)), "issued samples")
  # Re-reading it unchanged is fine, response still being full.
  expect_equal(prec_pooled(x)$moe, 1)
  expect_equal(nrow(predict(x, data.frame(moe = c(0.5, 1)))), 2L)
})

test_that("a respondent profile round-trips at any response rate", {
  x <- n_pooled(var = 100, moe = 1, occasions = 8, overlap = 0.75,
                cor_decay = 0.8)
  expect_equal(x$params$overlap_basis, "respondent")
  expect_silent(prec_pooled(x, resp_rate = 0.8))
  expect_equal(nrow(predict(x, data.frame(resp_rate = c(0.7, 0.9)))), 2L)
})

test_that("the basis survives an n_pooled round trip from a precision object", {
  x <- prec_pooled(var = 100, n = 500, occasions = 8,
                   overlap = design_overlap("4", max_lag = 7), cor_decay = 0.8)
  expect_equal(x$params$overlap_basis, "issued")
  expect_error(n_pooled(x, resp_rate = 0.8), "issued samples")
})

## T15. The exhaustion diagnostic counts interviews, not calendar span

test_that("interviews per cohort are recovered from the overlap profile", {
  # A cohort interviewed k times contributes k(k-1)/2 shared pairs against
  # k occasions of membership, so 1 + 2 * sum(overlap) is k whether or not
  # the schedule has gaps. Reading the last positive lag gives the span.
  for (spec in c("4", "6", "4-8-4", "1-1-0-0-1-1", "2-2-2")) {
    s <- design_overlap(spec, max_lag = 30)
    expect_equal(1 + 2 * sum(as.numeric(s)), sum(as.double(s$rotation)),
                 info = spec)
  }
})

test_that("the diagnostic names interviews, not the span, for a gapped schedule", {
  # "4-8-4" spans 16 occasions and interviews 8. A size that exhausts the
  # population must be reported against the 8: reading the last positive
  # lag would call it 16 and understate what the design consumes, enough
  # here to drop the sentence entirely.
  N <- 1e5
  ov <- as.numeric(design_overlap("4-8-4", max_lag = 19))[1:19]
  msg <- tryCatch(
    prec_pooled(var = 100, n = 0.8 * N, occasions = 20, N = N,
                overlap = ov, cor_decay = 0.9),
    error = function(e) conditionMessage(e)
  )
  expect_true(is.character(msg))
  expect_match(msg, "not a valid covariance")
  expect_match(msg, "each cohort 8 times")
  expect_match(msg, "draws 2.00 times the population")
})

test_that("the diagnostic is dropped when the profile is truncated", {
  # Shared units at the largest lag mean the life runs past the horizon and
  # sum(overlap) is incomplete, so no interview count can be claimed. These
  # profiles do fail the guard, so the branch is reached rather than skipped.
  N <- 1e5
  # The first three would trip the size condition if `complete` were
  # ignored, so the branch is load-bearing rather than incidentally quiet.
  for (ov in list(c(0.9, 0.05, 0.05), c(0.95, 0.02, 0.03), c(0.8, 0.05, 0.1),
                  c(0.9, 0.5, 0.3), c(0.99, 0.1, 0.99),
                  c(0.9, 0.1, 0.8, 0.1, 0.7))) {
    msg <- tryCatch(
      prec_pooled(var = 100, n = 0.8 * N, occasions = length(ov) + 1L,
                  N = N, overlap = ov, overlap_cor = 0.9),
      error = function(e) conditionMessage(e)
    )
    expect_true(is.character(msg), info = paste(ov, collapse = ","))
    expect_match(msg, "not a valid covariance")
    expect_false(grepl("each cohort", msg), info = paste(ov, collapse = ","))
  }
})

## T16. occasions is bounded

test_that("occasions is capped, the covariance being assembled densely", {
  expect_error(prec_pooled(var = 100, n = 500, occasions = 1001),
               "at most 1000")
  expect_error(n_pooled(var = 100, moe = 1, occasions = 1e6),
               "at most 1000")
  expect_silent(prec_pooled(var = 100, n = 500, occasions = 1000))
})
