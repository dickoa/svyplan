test_that("n_prop wald MOE mode", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_s3_class(result, "svyplan_n")
  expect_equal(result$type, "proportion")
  expect_equal(result$method, "wald")

  # Infinite pop: a=1, z=1.96, n = 1*1.96^2*0.3*0.7 / (0.05^2 + 1.96^2*0.3*0.7/Inf)
  z <- qnorm(0.975)
  expected <- z^2 * 0.3 * 0.7 / 0.05^2
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop wald CV mode", {
  result <- n_prop(p = 0.5, cv = 0.10)
  # a=1, q/p = 1, n = 1 * 1 / (0.01 + 1/Inf) = 100
  expect_equal(result$n, 100, tolerance = 1e-6)
})

test_that("n_prop wald MOE with FPC", {
  result <- n_prop(p = 0.3, moe = 0.05, N = 1000)
  z <- qnorm(0.975)
  a <- 1000 / 999
  expected <- a * z^2 * 0.3 * 0.7 / (0.05^2 + z^2 * 0.3 * 0.7 / 999)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop wald CV with FPC", {
  result <- n_prop(p = 0.5, cv = 0.10, N = 5000)
  a <- 5000 / 4999
  expected <- a * 1 / (0.01 + 1 / 4999)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop wilson method", {
  result <- n_prop(p = 0.3, moe = 0.05, method = "wilson")
  z <- qnorm(0.975)
  q <- 0.7
  e <- 0.05
  rad <- e^2 - 0.3 * q * (4 * e^2 - 0.3 * q)
  expected <- (0.3 * q - 2 * e^2 + sqrt(rad)) * (z / e)^2 / 2
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop logodds method", {
  result <- n_prop(p = 0.3, moe = 0.05, method = "logodds")
  z <- qnorm(0.975)
  q <- 0.7
  e <- 0.05
  kk <- q / 0.3
  rad <- e^2 * (1 + kk^2)^2 + kk^2 * (1 - 2 * e) * (1 + 2 * e)
  x <- (e * (1 + kk^2) + sqrt(rad)) / (kk * (1 - 2 * e))
  expected <- 1 / ((sqrt(0.3 * q) / z * log(x))^2)
  expect_equal(result$n, expected, tolerance = 1e-6)
})

test_that("n_prop deff multiplier works", {
  base <- n_prop(p = 0.3, moe = 0.05)
  with_deff <- n_prop(p = 0.3, moe = 0.05, deff = 2)
  expect_equal(with_deff$n, base$n * 2, tolerance = 1e-6)
})

test_that("n_prop returns correct S3 class", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_true(inherits(result, "svyplan_n"))
  expect_true(is.list(result))
  expect_equal(as.integer(result), ceiling(result$n))
  expect_equal(as.double(result), result$n)
})

test_that("n_prop validates inputs", {
  expect_error(n_prop(p = 0, moe = 0.05), "must be in \\(0, 1\\)")
  expect_error(n_prop(p = 1, moe = 0.05), "must be in \\(0, 1\\)")
  expect_error(n_prop(p = 0.3), "specify exactly one")
  expect_error(n_prop(p = 0.3, moe = 0.05, cv = 0.10), "specify exactly one")
  expect_s3_class(n_prop(p = 0.3, moe = 0.05, deff = 0.5), "svyplan_n")
  expect_error(n_prop(p = 0.3, moe = 0.05, deff = 0), "must be positive")
  expect_error(n_prop(p = 0.3, moe = 0.05, deff = -1), "must be positive")
  expect_error(n_prop(p = 0.3, moe = 0.05, N = -1), "greater than 1")
  expect_error(n_prop(p = 0.3, moe = 0.05, N = 1), "greater than 1")
  expect_error(n_prop(p = 0.3, cv = 0.10, method = "wilson"), "Wilson.*moe")
  expect_error(n_prop(p = 0.3, cv = 0.10, method = "logodds"), "Log-odds.*moe")
})

test_that("svyplan_n has se/moe/cv fields for proportion", {
  res <- n_prop(p = 0.3, moe = 0.05)
  expect_true(!is.null(res$se))
  expect_true(!is.null(res$moe))
  expect_true(!is.null(res$cv))
  expect_true(is.numeric(res$se))
  expect_true(res$se > 0)
  expect_true(res$moe > 0)
  expect_true(res$cv > 0)
})

test_that("se/moe/cv use Cochran FPC for proportion (infinite pop)", {
  res <- n_prop(p = 0.3, moe = 0.05)
  z <- qnorm(0.975)
  n_eff <- res$n
  se_expected <- sqrt(0.3 * 0.7 / n_eff)
  expect_equal(res$se, se_expected, tolerance = 1e-6)
  expect_equal(res$moe, z * se_expected, tolerance = 1e-6)
  expect_equal(res$cv, se_expected / 0.3, tolerance = 1e-6)
})

test_that("se/moe/cv use Cochran FPC for proportion (finite pop)", {
  res <- n_prop(p = 0.3, moe = 0.05, N = 500)
  z <- qnorm(0.975)
  n_eff <- res$n
  fpc <- (500 - n_eff) / (500 - 1)
  se_expected <- sqrt(0.3 * 0.7 * fpc / n_eff)
  expect_equal(res$se, se_expected, tolerance = 1e-6)
  expect_equal(res$moe, z * se_expected, tolerance = 1e-6)
})

test_that("proportion FPC round-trip: moe from n matches target", {
  res <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(res$moe, 0.05, tolerance = 1e-6)
})

test_that("proportion FPC round-trip with finite N", {
  res <- n_prop(p = 0.3, moe = 0.05, N = 1000)
  expect_equal(res$moe, 0.05, tolerance = 1e-6)
})

test_that("n_prop logodds rejects moe >= 0.5", {
  expect_error(n_prop(p = 0.3, moe = 0.5, method = "logodds"),
               "moe.*< 0.5")
  expect_error(n_prop(p = 0.3, moe = 0.6, method = "logodds"),
               "moe.*< 0.5")
})

test_that("deff multiplies the SRS variance at the actual net n (finite N)", {
  z <- qnorm(0.975)
  for (deff in c(0.8, 2)) {
    got <- n_prop(p = 0.5, moe = 0.1, N = 100, deff = deff)$n
    expected <- z^2 * deff * 0.25 * 100 / (0.1^2 * 99 + z^2 * deff * 0.25)
    expect_equal(got, expected, tolerance = 1e-10)
    pr <- prec_prop(p = 0.5, n = got, N = 100, deff = deff)
    expect_equal(pr$moe, 0.1, tolerance = 1e-8)
  }
})

test_that("a census has zero sampling variance regardless of deff", {
  expect_warning(pr <- prec_prop(p = 0.5, n = 100, N = 100, deff = 2),
                 "census")
  expect_equal(pr$se, 0)
  expect_warning(pm <- prec_mean(var = 100, n = 100, N = 100, deff = 2),
                 "census")
  expect_equal(pm$se, 0)
})

test_that("unattainable finite-population targets error", {
  expect_error(n_mean(var = 100, moe = 1, N = 100, deff = 2, resp_rate = 0.8),
               "unattainable")
  expect_error(n_prop(p = 0.5, moe = 0.02, N = 100, resp_rate = 0.5),
               "unattainable")
})

test_that("Wilson and log-odds sizes invert their precision counterparts", {
  for (method in c("wilson", "logodds")) {
    for (p in c(0.05, 0.3, 0.5, 0.9)) {
      for (moe in c(0.01, 0.03, 0.08)) {
        n <- n_prop(p = p, moe = moe, method = method)$n
        back <- prec_prop(p = p, n = n, method = method)$moe
        expect_equal(back, moe, tolerance = 1e-8)
      }
    }
  }
})

test_that("every method inverts its precision counterpart under deff, N, resp_rate", {
  for (method in c("wald", "wilson", "logodds")) {
    for (N in c(Inf, 10000, 2000, 600)) {
      for (deff in c(1, 2)) {
        for (resp_rate in c(1, 0.8)) {
          n <- n_prop(p = 0.2, moe = 0.04, method = method, N = N,
                      deff = deff, resp_rate = resp_rate)$n
          back <- prec_prop(p = 0.2, n = n, method = method, N = N,
                            deff = deff, resp_rate = resp_rate)$moe
          expect_equal(back, 0.04, tolerance = 1e-8)
        }
      }
    }
  }
})

test_that("the Wilson method applies the finite population correction", {
  # Sizes shrink as the frame shrinks, and the infinite-population default
  # is unchanged.
  sizes <- vapply(
    c(Inf, 1e5, 5000, 1000, 500),
    function(N) n_prop(p = 0.2, moe = 0.04, method = "wilson", N = N)$n,
    numeric(1L)
  )
  expect_true(all(diff(sizes) < 0))
  expect_equal(sizes[1L], 382.4532, tolerance = 1e-4)

  # A finite frame never needs more than a census, whatever the target.
  expect_lt(n_prop(p = 0.5, moe = 0.001, method = "wilson", N = 800)$n, 800)

  # The three methods agree closely at the same allocation, because they
  # now share one variance.
  moes <- vapply(
    c("wald", "wilson", "logodds"),
    function(m) prec_prop(p = 0.2, n = 400, N = 500, method = m)$moe,
    numeric(1L)
  )
  expect_equal(max(moes) / min(moes), 1, tolerance = 0.01)
})

test_that("Wilson and log-odds inversion holds with deff and finite N", {
  n <- n_prop(p = 0.25, moe = 0.04, method = "wilson", deff = 1.8)$n
  expect_equal(prec_prop(p = 0.25, n = n, method = "wilson", deff = 1.8)$moe,
               0.04, tolerance = 1e-8)
  n <- n_prop(p = 0.25, moe = 0.04, method = "logodds", deff = 1.8,
              N = 5000)$n
  expect_equal(
    prec_prop(p = 0.25, n = n, method = "logodds", deff = 1.8, N = 5000)$moe,
    0.04, tolerance = 1e-8
  )
})

test_that("log-odds converges to Wald as the margin of error shrinks", {
  # Both use Var(p_hat) = deff * N/(N-1) * p q (1/n - 1/N); the log-odds
  # interval is asymmetric, so they agree only in the limit.
  for (N in c(2000, Inf)) {
    ratios <- vapply(c(0.02, 0.005, 0.001, 0.0002), function(moe) {
      n_prop(p = 0.3, moe = moe, method = "logodds", N = N)$n /
        n_prop(p = 0.3, moe = moe, method = "wald", N = N)$n
    }, numeric(1L))
    expect_true(all(diff(abs(ratios - 1)) < 0))
    expect_equal(ratios[length(ratios)], 1, tolerance = 1e-6)
  }
})

test_that("margins of error at or above 1/2 are rejected, not silently wrong", {
  expect_error(n_prop(p = 0.5, moe = 0.5, method = "wilson"),
               "requires 'moe' < 0.5")
  expect_error(n_prop(p = 0.5, moe = 0.6, method = "wilson"),
               "requires 'moe' < 0.5")
  expect_error(n_prop(p = 0.5, moe = 0.5, method = "logodds"),
               "requires 'moe' < 0.5")
})

test_that("Wilson sizes fall with the margin of error and rise toward p = 0.5", {
  sizes <- vapply(c(0.02, 0.03, 0.05, 0.1), function(moe) {
    n_prop(p = 0.3, moe = moe, method = "wilson")$n
  }, numeric(1L))
  expect_true(all(diff(sizes) < 0))
  expect_gt(n_prop(p = 0.5, moe = 0.03, method = "wilson")$n,
            n_prop(p = 0.1, moe = 0.03, method = "wilson")$n)
})

test_that("the beta method reproduces survey::svyciprop(method = 'beta')", {
  skip_if_not_installed("survey")
  # An SRS design has deff 1 and no FPC, so the effective size svyciprop
  # derives from its variance estimate is the one svyplan uses, and the two
  # interval constructions must agree exactly.
  for (spec in list(c(n = 200, k = 6), c(n = 500, k = 25), c(n = 123, k = 13))) {
    y <- c(rep(1, spec[["k"]]), rep(0, spec[["n"]] - spec[["k"]]))
    design <- suppressWarnings(
      survey::svydesign(id = ~1, data = data.frame(y = y))
    )
    fit <- survey::svyciprop(~y, design, method = "beta", level = 0.95)
    p_hat <- as.numeric(coef(fit))
    n_eff <- p_hat * (1 - p_hat) / as.numeric(vcov(fit))
    expect_equal(
      as.numeric(.beta_limits(p_hat, n_eff, 0.05)),
      as.numeric(attr(fit, "ci")),
      tolerance = 1e-10
    )
  }
})

test_that("the degrees-of-freedom adjustment follows Korn-Graubard eq (2.2)", {
  # Korn and Graubard (1998), Table 7: the HHANES cocaine example has
  # n = 123 with 8 degrees of freedom, for which the paper states that the
  # second factor of (2.2) is 0.794.
  n_star <- 1000
  adjusted <- .kg_effective(n_star, n_net = 123, alpha = 0.10, df = 8)
  expect_equal(adjusted / n_star, 0.794, tolerance = 1e-3)

  # Fewer degrees of freedom widen the interval and so raise the size.
  expect_gt(
    n_prop(p = 0.05, moe = 0.02, method = "beta", df = 20)$n,
    n_prop(p = 0.05, moe = 0.02, method = "beta")$n
  )
  # Omitting 'df' is the same as claiming the variance is as stable as one
  # from a simple random sample of the same size, i.e. df = n_net - 1.
  expect_equal(
    prec_prop(p = 0.05, n = 400, method = "beta", df = 399)$moe,
    prec_prop(p = 0.05, n = 400, method = "beta")$moe
  )
  # A variance estimated with more degrees of freedom than that is more
  # stable, so the interval narrows.
  expect_lt(
    prec_prop(p = 0.05, n = 400, method = "beta", df = Inf)$moe,
    prec_prop(p = 0.05, n = 400, method = "beta")$moe
  )
})

test_that("the beta interval stays inside (0, 1) and is asymmetric", {
  fit <- prec_prop(p = 0.02, n = 150, method = "beta")
  ci <- confint(fit)
  expect_gt(ci[1L], 0)
  expect_lt(ci[2L], 1)
  # The upper arm is longer than the lower one for a rare outcome, which is
  # the behaviour Wald cannot reproduce: its lower limit is negative here.
  expect_gt(ci[2L] - 0.02, 0.02 - ci[1L])
  expect_lt(
    confint(prec_prop(p = 0.02, n = 150, method = "wald"))[1L] |> unname(),
    1e-12
  )
})

test_that("the beta method honours deff, N, resp_rate and df", {
  expect_lt(
    prec_prop(p = 0.08, n = 400, method = "beta", N = 600)$moe,
    prec_prop(p = 0.08, n = 400, method = "beta")$moe
  )
  expect_gt(
    prec_prop(p = 0.08, n = 400, method = "beta", deff = 2)$moe,
    prec_prop(p = 0.08, n = 400, method = "beta")$moe
  )
  expect_warning(
    res <- prec_prop(p = 0.08, n = 500, N = 500, method = "beta"),
    "net sample size .* >= population size"
  )
  expect_equal(res$moe, 0)

  expect_error(n_prop(p = 0.3, cv = 0.1, method = "beta"), "requires 'moe'")
  expect_error(n_prop(p = 0.3, moe = 0.05, method = "wilson", df = 10),
               "applies to method = 'beta' only")
  expect_error(prec_prop(p = 0.3, n = 100, method = "beta", df = 0.5),
               "'df' must be a number >= 1")
})

test_that("reported moe equals the confint half-width for every method", {
  for (method in c("wald", "wilson", "logodds", "beta")) {
    for (N in c(Inf, 5000, 600)) {
      for (deff in c(1, 2)) {
        fit <- prec_prop(p = 0.08, n = 400, method = method, N = N,
                         deff = deff)
        ci <- confint(fit)
        expect_equal(as.numeric((ci[2L] - ci[1L]) / 2), fit$moe,
                     tolerance = 1e-10)
      }
    }
  }
})
