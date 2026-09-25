bankier_frame <- function() {
  data.frame(stratum = c("A", "B"), N = c(1000, 1000),
             sd = c(2, 2), mean = c(10, 2), alloc_measure = c(10000, 2000))
}

test_that("Bankier uses population CV times the chosen measure to a power", {
  f <- bankier_frame()
  for (q in c(0, .5, 1)) {
    fit <- n_alloc(f, n = 60, alloc = "power", alloc_q = q)
    weights <- c(.2, 1) * c(10000, 2000)^q
    expect_equal(fit$detail$n, 60 * weights / sum(weights))
    expect_equal(sum(fit$detail$n_int), 60)
  }
  expect_equal(n_alloc(f, n = 60, alloc = "power", alloc_q = 0)$detail$n,
               c(10, 50))
  expect_equal(n_alloc(f, n = 60, alloc = "power", alloc_q = 1)$detail$n,
               c(30, 30))
})

test_that("default totals give Neyman at q one and population CVs at q zero", {
  f <- bankier_frame()
  f$alloc_measure <- NULL
  expect_equal(n_alloc(f, n = 60, alloc = "power", alloc_q = 1)$detail$n,
               n_alloc(f, n = 60, alloc = "neyman")$detail$n)
  expect_equal(n_alloc(f, n = 60, alloc = "power", alloc_q = 0)$detail$n,
               c(10, 50))
  f$alloc_measure <- f$N
  expect_equal(n_alloc(f, n = 60, alloc = "power", alloc_q = 1)$detail$n,
               c(10, 50))
})

test_that("Bankier supports signed means and is invariant to common scales", {
  f <- bankier_frame()
  base <- n_alloc(f, n = 60, alloc = "power")$detail$n
  f$alloc_measure <- f$alloc_measure * 1e200
  f$mean <- -f$mean
  expect_equal(n_alloc(f, n = 60, alloc = "power")$detail$n, base)
  f$alloc_measure <- NULL
  base <- n_alloc(f, n = 60, alloc = "power")$detail$n
  f$mean <- f$mean * 1e100
  f$sd <- f$sd * 1e100
  expect_equal(n_alloc(f, n = 60, alloc = "power")$detail$n, base)
})

test_that("Bankier validates means and measures and retains zero-variance bounds", {
  f <- bankier_frame()
  for (bad in c(0, NA_real_, Inf)) {
    f$mean[1] <- bad
    expect_error(n_alloc(f, n = 60, alloc = "power"), "mean|finite")
  }
  f <- bankier_frame()
  f$mean <- NULL
  expect_error(n_alloc(f, n = 60, alloc = "power"), "nonzero.*mean")
  for (bad in c(0, -1, NA_real_, Inf)) {
    f <- bankier_frame()
    f$alloc_measure[1] <- bad
    expect_error(n_alloc(f, n = 60, alloc = "power"), "alloc_measure.*positive finite")
  }
  f <- bankier_frame()
  f$sd[1] <- 0
  expect_equal(n_alloc(f, n = 60, alloc = "power", min_n_stratum = 4)$detail$n,
               c(4, 56))
  f$sd[] <- 0
  expect_warning(fit <- n_alloc(f, n = 60, alloc = "power"), "all 'sd'")
  expect_equal(fit$detail$n, c(30, 30))
})

test_that("custom Bankier measures survive prediction and precision round trips", {
  f <- bankier_frame()
  f$alloc_measure <- c(500, 8000)
  fit <- n_alloc(f, n = 60, alloc = "power", alloc_q = .3)
  precision <- prec_alloc(fit)
  expect_equal(n_alloc(precision, n = 60)$detail$n, fit$detail$n)
  expect_equal(precision$params$frame$alloc_measure, f$alloc_measure)
  grid <- predict(fit, data.frame(alloc_q = c(0, .5, 1)))
  for (i in 1:3) {
    direct <- n_alloc(f, n = 60, alloc = "power", alloc_q = c(0, .5, 1)[i])
    expect_equal(grid$cv[i], direct$cv)
  }
  sized <- n_alloc(f, cv = .1, alloc = "power")
  expect_equal(prec_alloc(sized)$cv, .1, tolerance = 1e-6)
  f$unit_cost <- c(1, 3)
  budget <- n_alloc(f, budget = 120, alloc = "power")
  weights <- f$sd / abs(f$mean) * sqrt(f$alloc_measure)
  expect_equal(budget$detail$n, 120 * weights / sum(weights * f$unit_cost))
})

test_that("boundary methods recompute Bankier population CVs and totals", {
  x <- exp(seq(0, 4, length.out = 200))
  for (method in c("cumrootf", "geo", "lh", "kozak")) {
    set.seed(928)
    fit <- strata_bound(x, n_strata = 3, n = 60, method = method,
                         alloc = "power", alloc_q = .5)
    f <- fit$strata
    weights <- f$sd / abs(f$mean) * sqrt(f$N * abs(f$mean))
    continuous <- .rna_alloc(weights, 60, pmin(2, f$N), f$N)
    expected <- .round_oric_bounded(continuous, as.integer(pmin(2, f$N)),
                                    as.integer(f$N))
    expect_equal(f$n, expected)
    expect_equal(prec_alloc(f, n = f$n)$cv, fit$cv)
    target <- strata_bound(x, n_strata = 3, cv = .08, method = method,
                            alloc = "power", alloc_q = .5)
    expect_lte(target$cv, .08 + 1e-8)
  }
})

test_that("alloc_measure is refused where the allocation cannot use it", {
  fr <- data.frame(
    stratum = c("a", "b", "c"), N = c(1000, 2000, 500),
    sd = c(10, 20, 5), mean = c(50, 80, 20),
    alloc_measure = c(3, 1, 2)
  )
  expect_error(n_alloc(fr, n = 300, alloc = "neyman"), "alloc = \"power\"")
  expect_error(n_alloc(fr, n = 300, alloc = "proportional"), "not alloc = \"proportional\"")
  expect_no_error(n_alloc(fr, n = 300, alloc = "power"))

  measures <- data.frame(stratum = fr$stratum, name = "y",
                         mean = fr$mean, sd = fr$sd)
  targets <- data.frame(name = "y", cv = 0.05)
  expect_error(
    n_alloc(fr[c("stratum", "N", "alloc_measure")], measures = measures,
            targets = targets),
    "not used for joint"
  )
})

test_that("Bankier stratification refuses a variable taking both signs", {
  set.seed(1)
  expect_error(
    strata_bound(rnorm(500), n_strata = 3, n = 100, alloc = "power"),
    "both signs"
  )
  expect_no_error(
    strata_bound(rlnorm(500), n_strata = 3, n = 100, alloc = "power")
  )
  expect_no_error(
    strata_bound(rnorm(500), n_strata = 3, n = 100, alloc = "neyman")
  )
})
