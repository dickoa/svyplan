test_that("prec_prop computes SE, MOE, CV for wald method", {
  result <- prec_prop(p = 0.3, n = 400)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "proportion")
  expect_equal(result$method, "wald")

  z <- qnorm(0.975)
  se_exp <- sqrt(0.3 * 0.7 / 400)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
  expect_equal(result$moe, z * se_exp, tolerance = 1e-6)
  expect_equal(result$cv, se_exp / 0.3, tolerance = 1e-6)
})

test_that("prec_prop with FPC uses Cochran form", {
  result <- prec_prop(p = 0.3, n = 400, N = 5000)
  n_eff <- 400
  fpc <- (5000 - n_eff) / (5000 - 1)
  se_exp <- sqrt(0.3 * 0.7 * fpc / n_eff)
  expect_equal(result$se, se_exp, tolerance = 1e-6)
})

test_that("prec_prop with deff", {
  base <- prec_prop(p = 0.3, n = 400)
  deff2 <- prec_prop(p = 0.3, n = 400, deff = 2)
  expect_equal(deff2$se, base$se * sqrt(2), tolerance = 1e-6)
})

test_that("prec_prop with resp_rate", {
  base <- prec_prop(p = 0.3, n = 400)
  rr <- prec_prop(p = 0.3, n = 400, resp_rate = 0.8)
  n_eff <- 400 * 0.8
  se_exp <- sqrt(0.3 * 0.7 / n_eff)
  expect_equal(rr$se, se_exp, tolerance = 1e-6)
})

test_that("prec_prop validates inputs", {
  expect_error(prec_prop(p = 0, n = 400), "must be in \\(0, 1\\)")
  expect_error(prec_prop(p = 1, n = 400), "must be in \\(0, 1\\)")
  expect_error(prec_prop(p = 0.3, n = -1), "must be positive")
  expect_s3_class(prec_prop(p = 0.3, n = 400, deff = 0.5), "svyplan_prec")
  expect_error(prec_prop(p = 0.3, n = 400, deff = 0), "must be positive")
  expect_error(prec_prop(p = 0.3, n = 400, deff = -1), "must be positive")
  expect_error(prec_prop(p = 0.3, n = 400, resp_rate = 0), "resp_rate")
  expect_error(prec_prop(p = 0.3, n = 400, resp_rate = 1.5), "resp_rate")
})

test_that("prec_prop prints correctly", {
  result <- prec_prop(p = 0.3, n = 400)
  out <- capture.output(print(result))
  expect_match(out[1], "Sampling precision for proportion")
  expect_match(out[2], "n = 400")
  expect_match(out[3], "se =")
})

test_that("prec_prop format returns string", {
  result <- prec_prop(p = 0.3, n = 400)
  expect_true(is.character(format(result)))
  expect_match(format(result), "svyplan_prec")
})

test_that("prec_prop logodds works for extreme p", {
  res <- prec_prop(p = 0.9, n = 20, method = "logodds")
  expect_true(res$moe > 0 && res$moe < 0.5)
  expect_true(res$se > 0)

  res2 <- prec_prop(p = 0.05, n = 100, method = "logodds")
  expect_true(res2$moe > 0 && res2$moe < 0.5)
})

test_that("prec_prop logodds extreme p round-trip", {
  for (pp in c(0.05, 0.95)) {
    s1 <- n_prop(p = pp, moe = 0.04, method = "logodds")
    p1 <- prec_prop(s1)
    s2 <- n_prop(p1)
    expect_equal(s2$n, s1$n, tolerance = 1)
  }
})

test_that("a census returns moe = 0 with one warning under every method", {
  for (m in c("wald", "wilson", "logodds")) {
    expect_warning(
      res <- prec_prop(p = 0.5, n = 50, N = 50, method = m),
      "net sample size .* >= population size"
    )
    expect_equal(res$moe, 0)
    expect_equal(res$se, 0)
  }
})

test_that("supplied gross n above a finite N is rejected", {
  expect_error(prec_prop(p = 0.5, n = 120, N = 100, resp_rate = 0.8),
               "cannot draw 120 units")
  expect_error(prec_mean(var = 1, n = 120, N = 100, resp_rate = 0.8),
               "cannot draw 120 units")
  expect_warning(res <- prec_prop(p = 0.5, n = 100, N = 100), "census")
  expect_equal(res$se, 0)
})

## Solving for the proportion

test_that("solving for p inverts the precision it reports", {
  # the solve runs on the log scale, so accuracy is relative to the root and
  # holds for a threshold of 1e-5 as well as one of 0.5
  for (m in c("wald", "wilson", "logodds", "beta")) {
    for (n in c(200, 1500, 1e5)) {
      for (target in c(0.02, 0.10, 0.50)) {
        r <- prec_prop(n = n, cv = target, N = 2e6, method = m)
        expect_identical(r$solved, "p")
        expect_equal(r$cv, target, tolerance = 1e-10)
        # the reported precision is what the solved p actually achieves
        expect_equal(prec_prop(p = r$params$p, n = n, N = 2e6, method = m)$cv,
                     target, tolerance = 1e-10)
      }
    }
  }
})

test_that("the Wald solution is the closed form", {
  n <- 1500
  cv <- 0.1
  n_eff <- n / (1 + 0)  # infinite N, no deff, full response
  expect_equal(prec_prop(n = n, cv = cv)$params$p, 1 / (1 + n_eff * cv^2))
})

test_that("a solved p is a floor: every larger proportion meets the target", {
  r <- prec_prop(n = 1500, cv = 0.10, N = 2e6, method = "wilson")
  p0 <- r$params$p
  for (mult in c(1.01, 1.5, 4)) {
    expect_lt(prec_prop(p = min(p0 * mult, 0.99), n = 1500, N = 2e6,
                        method = "wilson")$cv, 0.10)
  }
  expect_gt(prec_prop(p = p0 * 0.99, n = 1500, N = 2e6, method = "wilson")$cv,
            0.10)
})

test_that("solve mode round trips through n_prop", {
  r <- prec_prop(n = 1500, cv = 0.10, N = 2e6, deff = 1.5, resp_rate = 0.8)
  expect_equal(n_prop(r)$n, 1500)
})

test_that("p, cv, and rmoe are mutually exclusive and one is required", {
  expect_error(prec_prop(n = 400), "exactly one of 'p', 'cv', or 'rmoe'")
  expect_error(prec_prop(p = 0.3, cv = 0.1, n = 400),
               "exactly one of 'p', 'cv', or 'rmoe'")
})

test_that("a census has no proportion to solve for", {
  expect_error(prec_prop(n = 500, cv = 0.1, N = 500),
               "census|exceeds the population")
})

test_that("every positive cv target has a root, under every method", {
  # The sampling CV falls monotonically in p and covers (0, Inf) over (0, 1),
  # so there is no unattainable case. The target below was previously refused
  # because the logit-scale half-width it was read through turns back upward
  # above about 0.999.
  r <- prec_prop(n = 50, cv = 0.02, method = "logodds")
  expect_equal(r$solved, "p")
  expect_equal(r$params$p, 1 / 1.02, tolerance = 1e-10)
  expect_equal(r$cv, 0.02, tolerance = 1e-10)
})

test_that("se is the sampling standard error under every method", {
  n_eff <- 900 / 2
  target <- sqrt(0.02 * 0.98 / n_eff)
  for (m in c("wald", "wilson", "logodds", "beta")) {
    res <- prec_prop(p = 0.02, n = 900, deff = 2, method = m, df = 25)
    expect_equal(res$se, target, tolerance = 1e-12)
    expect_equal(res$cv, target / 0.02, tolerance = 1e-12)
  }
})

test_that("the interval still depends on the method it was built with", {
  moes <- vapply(c("wald", "wilson", "logodds", "beta"), function(m) {
    prec_prop(p = 0.02, n = 900, deff = 2, method = m, df = 25)$moe
  }, numeric(1L))
  expect_equal(length(unique(round(moes, 10))), 4L)

  widths <- vapply(c("wald", "wilson", "logodds", "beta"), function(m) {
    ci <- confint(prec_prop(p = 0.02, n = 900, deff = 2, method = m, df = 25))
    ci[2L] - ci[1L]
  }, numeric(1L))
  expect_equal(length(unique(round(widths, 10))), 4L)
})

test_that("solving p from cv is the same under every method", {
  ps <- vapply(c("wald", "wilson", "logodds", "beta"), function(m) {
    prec_prop(n = 1500, cv = 0.10, N = 2e6, method = m)$params$p
  }, numeric(1L))
  expect_equal(length(unique(round(ps, 12))), 1L)
})

test_that("the default direction keeps its schema", {
  r <- prec_prop(p = 0.3, n = 400)
  expect_null(r$solved)
})

test_that("predict on a solved result varies cv, not p", {
  r <- prec_prop(n = 1500, cv = 0.10, N = 2e6)
  grid <- predict(r, expand.grid(n = c(500, 1500, 5000)))
  expect_true("p" %in% names(grid))
  expect_equal(grid$cv, rep(0.10, 3))
  expect_true(all(diff(grid$p) < 0))
  expect_equal(grid$p[2], r$params$p)
  expect_error(predict(r, expand.grid(p = c(0.1, 0.2))), "unknown parameter")
})

## Claims the "Choosing a method" section of ?n_prop makes

# The prose there was wrong on four counts before it was measured. These pin
# the measurements, so a change in behaviour fails a test rather than turning
# the help page into fiction.

prop_methods <- c("wald", "wilson", "logodds", "beta")

test_that("df widens the margin of error under every method", {
  for (m in prop_methods) {
    wide <- prec_prop(p = 0.1, n = 100, method = m, df = 10)$moe
    plain <- prec_prop(p = 0.1, n = 100, method = m)$moe
    expect_gt(wide, plain)
  }
  # and it is the same t quantile doing it in each case
  ratio <- vapply(prop_methods, function(m) {
    prec_prop(p = 0.5, n = 400, method = m, df = 12)$moe /
      prec_prop(p = 0.5, n = 400, method = m)$moe
  }, numeric(1L))
  expect_equal(unname(ratio), rep(qt(0.975, 12) / qnorm(0.975), 4),
               tolerance = 1e-3)
})

test_that("beta is the widest method over most of the range but not all", {
  grid <- expand.grid(p = c(0.01, 0.02, 0.05, 0.1, 0.3, 0.5),
                      n = c(100, 250, 500, 1000, 5000), deff = c(1, 2))
  widest <- vapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    w <- vapply(prop_methods, function(m) {
      prec_prop(p = g$p, n = g$n, deff = g$deff, method = m)$moe
    }, numeric(1L))
    prop_methods[which.max(w)]
  }, character(1L))
  expect_gt(mean(widest == "beta"), 0.7)
  expect_true(any(widest != "beta"))
  # the exceptions are rare outcomes at small samples, and log-odds owns them
  expect_setequal(unique(widest[widest != "beta"]), "logodds")
  expect_lte(max(grid$p[widest != "beta"]), 0.05)
  expect_lte(max(grid$n[widest != "beta"]), 500)
})

test_that("only Wald leaves the parameter space and needs truncating", {
  grid <- expand.grid(p = c(0.001, 0.005, 0.01, 0.05, 0.5, 0.95, 0.99),
                      n = c(20, 30, 60, 100, 500), deff = c(1, 2.5))
  hits <- function(m) {
    vapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      x <- prec_prop(p = g$p, n = g$n, deff = g$deff, method = m)
      ci <- confint(x)
      # a truncated interval is one whose half-length no longer equals moe
      abs((ci[2] - ci[1]) / 2 - x$moe) > 1e-9
    }, logical(1L))
  }
  expect_true(any(hits("wald")))
  for (m in c("wilson", "logodds", "beta")) expect_false(any(hits(m)))

  # Wilson and log-odds stay strictly inside (0, 1); beta inside [0, 1]
  for (m in c("wilson", "logodds")) {
    ci <- confint(prec_prop(p = 0.005, n = 30, method = m))
    expect_gt(ci[1], 0)
    expect_lt(ci[2], 1)
  }
  ci_b <- confint(prec_prop(p = 0.005, n = 30, method = "beta"))
  expect_gte(ci_b[1], 0)
  expect_lte(ci_b[2], 1)
})

test_that("the four methods converge slowly in the sample size", {
  spread <- function(n) {
    ps <- seq(0.1, 0.9, by = 0.05)
    max(vapply(ps, function(p) {
      w <- vapply(prop_methods, function(m) {
        prec_prop(p = p, n = n, method = m)$moe
      }, numeric(1L))
      max(w) / min(w) - 1
    }, numeric(1L)))
  }
  # the figures quoted in ?n_prop, to the precision it quotes them
  expect_equal(100 * spread(100), 8.2, tolerance = 0.05)
  expect_equal(100 * spread(500), 3.8, tolerance = 0.05)
  expect_equal(100 * spread(1000), 2.7, tolerance = 0.05)
  expect_equal(100 * spread(5000), 1.2, tolerance = 0.05)
  # "less than a percent" is not true at any of them
  expect_gt(spread(5000), 0.01)
})

test_that("Wald is the only method that solves a cv target", {
  expect_s3_class(n_prop(p = 0.1, cv = 0.1, method = "wald"), "svyplan_n")
  for (m in c("wilson", "logodds", "beta")) {
    expect_error(n_prop(p = 0.1, cv = 0.1, method = m), "'moe'")
  }
})
