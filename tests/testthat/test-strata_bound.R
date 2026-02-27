set.seed(123)
x_lnorm <- rlnorm(1000, meanlog = 6, sdlog = 1.5)
x_unif <- runif(500, 10, 100)

test_that("validates x input", {

  expect_error(strata_bound("abc"), "numeric vector")
  expect_error(strata_bound(1), "at least 2")
  expect_error(strata_bound(c(1, NA, 3)), "NA")
})

test_that("validates n_strata", {
  expect_error(strata_bound(x_unif, n_strata = 1), "integer >= 2")
  expect_error(strata_bound(x_unif, n_strata = NA), "integer >= 2")
  expect_error(strata_bound(c(1, 1, 1), n_strata = 4, n = 2, method = "cumrootf"),
               "fewer unique")
})

test_that("validates method", {
  expect_error(strata_bound(x_unif, method = "invalid"), "'arg' should be one of")
})

test_that("lh and kozak require n or cv", {
  expect_error(strata_bound(x_unif, method = "lh"), "requires")
  expect_error(strata_bound(x_unif, method = "kozak"), "requires")
})

test_that("cumrootf and geo work without n or cv", {
  res <- strata_bound(x_unif, n_strata = 3, method = "cumrootf")
  expect_s3_class(res, "svyplan_strata")
  res2 <- strata_bound(x_lnorm, n_strata = 3, method = "geo")
  expect_s3_class(res2, "svyplan_strata")
})

test_that("cannot specify both n and cv", {
  expect_error(strata_bound(x_unif, n = 100, cv = 0.05), "at most one")
})

test_that("validates alloc", {
  expect_error(strata_bound(x_unif, n = 100, alloc = 42),
               "must be one of")
  expect_error(strata_bound(x_unif, n = 100, alloc = "invalid"),
               "'arg' should be one of")
  expect_error(strata_bound(x_unif, n = 100, alloc = list(q1 = 1)),
               "must be one of")
  expect_error(strata_bound(x_unif, n = 100, alloc = "power", q = 2),
               "numeric scalar in \\[0, 1\\]")
  expect_error(strata_bound(x_unif, n = 100, alloc = "power", q = -0.1),
               "numeric scalar in \\[0, 1\\]")
  expect_error(strata_bound(x_unif, n = 100, alloc = "power", q = "a"),
               "numeric scalar in \\[0, 1\\]")
})

test_that("validates cost", {
  expect_error(strata_bound(x_unif, n = 100, cost = -1), "positive")
  expect_error(strata_bound(x_unif, n = 100, cost = c(1, NA)), "positive")
})

test_that("validates certain", {
  expect_error(strata_bound(x_unif, n = 100, certain = "a"), "numeric scalar")
  expect_error(strata_bound(x_unif, n = 100, certain = max(x_unif) + 1),
               "no units")
  expect_error(strata_bound(x_unif, n = 100, certain = min(x_unif) - 1),
               "all units")
})

test_that("certain marks only the last stratum as take-all", {
  set.seed(1)
  x <- c(runif(95, 1, 100), runif(5, 200, 300))
  thr <- as.numeric(quantile(x, 0.95))
  res <- strata_bound(x, n_strata = 3, n = 30, certain = thr, method = "cumrootf")
  expect_equal(which(res$strata$certain), 3L)
  expect_equal(res$strata$n_h[3], res$strata$N_h[3])
})

test_that("cumrootf: uniform data yields reasonable strata", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_s3_class(res, "svyplan_strata")
  expect_equal(res$n_strata, 3L)
  expect_length(res$boundaries, 2L)
  expect_true(all(res$boundaries > min(x_unif)))
  expect_true(all(res$boundaries < max(x_unif)))
  expect_true(all(diff(res$boundaries) > 0))
})

test_that("cumrootf: custom nclass respected", {
  res1 <- strata_bound(x_lnorm, n_strata = 3, n = 100,
                        method = "cumrootf", nclass = 50)
  res2 <- strata_bound(x_lnorm, n_strata = 3, n = 100,
                        method = "cumrootf", nclass = 200)
  expect_s3_class(res1, "svyplan_strata")
  expect_s3_class(res2, "svyplan_strata")
})

test_that("cumrootf: 4 strata yields 3 boundaries", {
  res <- strata_bound(x_lnorm, n_strata = 4, n = 150, method = "cumrootf")
  expect_length(res$boundaries, 3L)
  expect_equal(nrow(res$strata), 4L)
})

test_that("geo: boundaries form geometric progression", {
  x_pos <- x_lnorm
  res <- strata_bound(x_pos, n_strata = 4, n = 200, method = "geo")
  bk <- c(min(x_pos), res$boundaries, max(x_pos))
  ratios <- bk[-1] / bk[-length(bk)]
  expect_true(max(abs(diff(ratios))) < 0.01)
})

test_that("geo: errors on non-positive x", {
  x_neg <- c(-1, 1:10)
  expect_error(strata_bound(x_neg, n_strata = 3, method = "geo"), "positive")
})

test_that("geo: correct boundary count", {
  res <- strata_bound(x_lnorm, n_strata = 5, n = 200, method = "geo")
  expect_length(res$boundaries, 4L)
})

test_that("lh: converges on lognormal data", {
  skip_on_cran()
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "lh")
  expect_s3_class(res, "svyplan_strata")
  expect_true(res$converged)
  expect_equal(res$method, "lh")
})

test_that("kozak: converges on lognormal data", {
  skip_on_cran()
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "kozak",
                       niter = 5L, maxiter = 50L)
  expect_s3_class(res, "svyplan_strata")
  expect_equal(res$method, "kozak")
})

test_that("kozak cv matches target approximately", {
  skip_on_cran()
  target_cv <- 0.10
  res <- strata_bound(x_lnorm, n_strata = 4, cv = target_cv, method = "kozak",
                       niter = 10L, maxiter = 100L)
  expect_true(res$cv <= target_cv * 1.5)
})

test_that("kozak n matches target approximately", {
  skip_on_cran()
  target_n <- 200
  res <- strata_bound(x_lnorm, n_strata = 3, n = target_n, method = "kozak",
                       niter = 5L, maxiter = 50L)
  expect_true(abs(res$n - target_n) / target_n < 0.5)
})

test_that("proportional allocation: n_h proportional to N_h", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                       alloc = "proportional")
  df <- res$strata
  prop_alloc <- df$n_h / sum(df$n_h)
  prop_pop <- df$N_h / sum(df$N_h)
  expect_true(max(abs(prop_alloc - prop_pop)) < 0.15)
})

test_that("neyman allocation: n_h proportional to N_h * S_h", {
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                       alloc = "neyman")
  df <- res$strata
  expected_prop <- df$N_h * df$S_h
  expected_prop <- expected_prop / sum(expected_prop)
  actual_prop <- df$n_h / sum(df$n_h)
  expect_true(cor(actual_prop, expected_prop) > 0.8)
})

test_that("power allocation works", {
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                       alloc = "power", q = 0.5)
  expect_s3_class(res, "svyplan_strata")
  expect_equal(res$alloc, "power")
  expect_equal(res$params$q, 0.5)
})

test_that("n_h >= 2 per stratum", {
  res <- strata_bound(x_lnorm, n_strata = 5, n = 50, method = "cumrootf")
  expect_true(all(res$strata$n_h >= 2))
})

test_that("take-all stratum works", {
  skip_on_cran()
  thresh <- quantile(x_lnorm, 0.90)
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, certain = thresh)
  expect_s3_class(res, "svyplan_strata")
  expect_true(any(res$strata$certain))
  certain_row <- res$strata[res$strata$certain, ]
  expect_equal(certain_row$n_h, certain_row$N_h)
})

test_that("output class is svyplan_strata", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_s3_class(res, "svyplan_strata")
  expect_true(inherits(res, "list"))
})

test_that("strata df has correct columns", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expected_cols <- c("stratum", "lower", "upper", "N_h", "W_h", "S_h",
                     "n_h", "certain")
  expect_equal(names(res$strata), expected_cols)
})

test_that("as.integer returns ceiling(n)", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
})

test_that("as.double returns boundaries", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_equal(as.double(res), res$boundaries)
})

test_that("as.data.frame returns strata df", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_identical(as.data.frame(res), res$strata)
})

test_that("print returns invisible(x)", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  out <- capture.output(val <- print(res))
  expect_identical(val, res)
  expect_true(length(out) > 0)
})

test_that("format returns informative string", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  fmt <- format(res)
  expect_true(grepl("svyplan_strata", fmt))
  expect_true(grepl("cumrootf", fmt))
  expect_true(grepl("3 strata", fmt))
})

test_that("boundaries partition x into exactly n_strata groups", {
  res <- strata_bound(x_lnorm, n_strata = 4, n = 200, method = "cumrootf")
  bins <- findInterval(x_lnorm, res$boundaries, left.open = TRUE) + 1L
  expect_equal(length(unique(bins)), 4L)
  expect_equal(sum(res$strata$N_h), length(x_lnorm))
})

test_that("sum(n_h) equals n", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_equal(sum(res$strata$n_h), res$n)
})

test_that("works on simulated lognormal (realistic skewed data)", {
  skip_on_cran()
  set.seed(999)
  x_skew <- rlnorm(2000, meanlog = 8, sdlog = 2)
  res <- strata_bound(x_skew, n_strata = 5, n = 500, method = "kozak",
                       niter = 5L, maxiter = 50L)
  expect_s3_class(res, "svyplan_strata")
  expect_equal(nrow(res$strata), 5L)
  expect_true(all(res$strata$N_h > 0))
})

test_that("lh with cv mode works", {
  skip_on_cran()
  res <- strata_bound(x_lnorm, n_strata = 3, cv = 0.10, method = "lh")
  expect_s3_class(res, "svyplan_strata")
  expect_true(res$n > 0)
})

test_that("kozak outperforms or matches cumrootf", {
  skip_on_cran()
  res_cr <- strata_bound(x_lnorm, n_strata = 4, n = 200, method = "cumrootf")
  res_kz <- strata_bound(x_lnorm, n_strata = 4, n = 200, method = "kozak",
                          niter = 10L, maxiter = 100L)
  expect_true(res_kz$cv <= res_cr$cv * 1.2)
})

test_that("cost parameter is scalar-recycled", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                       cost = 5)
  expect_s3_class(res, "svyplan_strata")
})

test_that("cost parameter with per-stratum vector", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                       cost = c(1, 2, 5))
  expect_s3_class(res, "svyplan_strata")
})

test_that("params captures expected fields", {
  skip_on_cran()
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "kozak",
                       niter = 15L, maxiter = 50L)
  expect_equal(res$params$N, length(x_lnorm))
  expect_equal(res$params$maxiter, 50L)
  expect_equal(res$params$niter, 15L)
})

test_that("2 strata works (single boundary)", {
  res <- strata_bound(x_unif, n_strata = 2, n = 80, method = "cumrootf")
  expect_length(res$boundaries, 1L)
  expect_equal(nrow(res$strata), 2L)
})

test_that("lh handles skewed data without NA crash", {
  skip_on_cran()
  set.seed(5)
  x_skew <- rlnorm(500, meanlog = 8, sdlog = 3)
  res <- strata_bound(x_skew, n_strata = 4, n = 100, method = "lh")
  expect_s3_class(res, "svyplan_strata")
})

test_that("lh maxiter > 1 improves over maxiter = 1 on skewed data", {
  skip_on_cran()
  set.seed(42)
  x_skew <- rlnorm(800, meanlog = 6, sdlog = 2)
  res1 <- strata_bound(x_skew, n_strata = 4, n = 200, method = "lh",
                        maxiter = 1L)
  res200 <- strata_bound(x_skew, n_strata = 4, n = 200, method = "lh",
                          maxiter = 200L)
  expect_true(res200$cv <= res1$cv)
})

test_that("power alloc q = 1 matches neyman", {
  res_ney <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                           alloc = "neyman")
  res_pow <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                           alloc = "power", q = 1)
  expect_equal(res_pow$strata$n_h, res_ney$strata$n_h)
})

test_that("power alloc q = 0 differs from neyman on skewed data", {
  set.seed(7)
  x_skew <- rlnorm(2000, meanlog = 6, sdlog = 2)
  res_ney <- strata_bound(x_skew, n_strata = 4, n = 400, method = "cumrootf",
                           alloc = "neyman")
  res_pow <- strata_bound(x_skew, n_strata = 4, n = 400, method = "cumrootf",
                           alloc = "power", q = 0)
  expect_false(identical(res_pow$strata$n_h, res_ney$strata$n_h))
  expect_equal(res_pow$alloc, "power")
})

test_that("power alloc uses default q = 0.5", {
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                       alloc = "power")
  expect_equal(res$alloc, "power")
  expect_equal(res$params$q, 0.5)
})

test_that("print shows allocation label", {
  res_ney <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                           alloc = "neyman")
  out_ney <- capture.output(print(res_ney))
  expect_true(any(grepl("Allocation: neyman", out_ney)))

  res_pow <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                           alloc = "power", q = 0.3)
  out_pow <- capture.output(print(res_pow))
  expect_true(any(grepl("power \\(q = 0\\.30\\)", out_pow)))
})

test_that("alloc field stores method name for all methods", {
  for (a in c("proportional", "neyman", "optimal")) {
    res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                         alloc = a)
    expect_equal(res$alloc, a)
    expect_null(res$params$q)
  }
})
