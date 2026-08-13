set.seed(123)
x_lnorm <- rlnorm(1000, meanlog = 6, sdlog = 1.5)
x_unif <- runif(500, 10, 100)

test_that("validates x input", {

  expect_error(strata_bound("abc"), "numeric vector")
  expect_error(strata_bound(1), "at least 2")
  expect_error(strata_bound(c(1, NA, 3)), "NA")
  expect_error(strata_bound(c(1, Inf, 3)), "finite")
  expect_error(strata_bound(c(1, -Inf, 3)), "finite")
})

test_that("validates n_strata", {
  expect_error(strata_bound(x_unif), "'n_strata' is required")
  expect_error(strata_bound(x_unif, n_strata = 1), "integer >= 2")
  expect_error(strata_bound(x_unif, n_strata = NA), "integer >= 2")
  expect_error(strata_bound(x_unif, n_strata = 3.5), "integer >= 2")
  expect_error(strata_bound(x_unif, n_strata = "3"), "integer >= 2")
  expect_error(strata_bound(c(1, 1, 1), n_strata = 4, n = 2, method = "cumrootf"),
               "fewer unique")
})

test_that("validates whole-number controls without rejecting whole doubles", {
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100.5, method = "cumrootf"),
               "'n' must be an integer")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf", n_class = 2.5),
               "'n_class' must be an integer")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, method = "lh", max_iter = 0),
               "'max_iter' must be an integer")
  expect_error(
    strata_bound(x_unif, n_strata = 3, n = 100, method = "kozak", n_restart = 1.5),
    "'n_restart' must be an integer"
  )
  expect_error(
    strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                 n_class = .Machine$integer.max + 1),
    "'n_class' must be an integer"
  )

  expect_no_error(
    strata_bound(x_unif, n_strata = 3.0, n = 100.0, method = "cumrootf",
                 n_class = 50.0)
  )
  expect_no_error(
    strata_bound(x_unif, n_strata = 3.0, n = 100.0, method = "kozak",
                 max_iter = 200.0, n_restart = 30.0)
  )
})

test_that("validates finite thresholds and costs", {
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                            take_all_above = Inf),
               "finite numeric scalar")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                            unit_cost = Inf),
               "positive and finite")
})

test_that("validates method", {
  expect_error(strata_bound(x_unif, n_strata = 3, method = "invalid"), "'arg' should be one of")
})

test_that("lh and kozak require n or cv", {
  expect_error(strata_bound(x_unif, n_strata = 3, method = "lh"), "requires")
  expect_error(strata_bound(x_unif, n_strata = 3, method = "kozak"), "requires")
})

test_that("cumrootf and geo work without n or cv", {
  res <- strata_bound(x_unif, n_strata = 3, method = "cumrootf")
  expect_s3_class(res, "svyplan_strata")
  res2 <- strata_bound(x_lnorm, n_strata = 3, method = "geo")
  expect_s3_class(res2, "svyplan_strata")
})

test_that("cannot specify both n and cv", {
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, cv = 0.05), "at most one")
})

test_that("errors when requested n is below feasible minimum", {
  expect_error(
    strata_bound(x_unif, n_strata = 4, n = 3, method = "cumrootf"),
    "minimum feasible"
  )
})

test_that("errors when requested n is above feasible maximum", {
  expect_error(
    strata_bound(x_unif, n_strata = 3, n = length(x_unif) + 1, method = "cumrootf"),
    "maximum feasible"
  )
})

test_that("validates alloc", {
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, alloc = 42),
               "must be one of")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, alloc = "invalid"),
               "'arg' should be one of")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, alloc = list(q1 = 1)),
               "must be one of")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, alloc = "power", alloc_q =2),
               "numeric scalar in \\[0, 1\\]")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, alloc = "power", alloc_q =-0.1),
               "numeric scalar in \\[0, 1\\]")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, alloc = "power", alloc_q ="a"),
               "numeric scalar in \\[0, 1\\]")
})

test_that("validates unit_cost", {
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, unit_cost = -1), "positive")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, unit_cost = c(1, NA)), "positive")
})

test_that("validates take_all", {
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, take_all_above = "a"), "numeric scalar")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, take_all_above = max(x_unif) + 1),
               "no units")
  expect_error(strata_bound(x_unif, n_strata = 3, n = 100, take_all_above = min(x_unif) - 1),
               "all units")
})

test_that("take_all marks only the last stratum as take-all", {
  set.seed(1)
  x <- c(runif(95, 1, 100), runif(5, 200, 300))
  thr <- as.numeric(quantile(x, 0.95))
  res <- strata_bound(x, n_strata = 3, n = 30, take_all_above = thr, method = "cumrootf")
  expect_equal(which(res$strata$take_all), 3L)
  expect_equal(res$strata$n[3], res$strata$N[3])
})

test_that("take_all includes equality and CV sizing accounts for the census stratum", {
  x <- c(1:89, rep(90, 5), 91:100)
  res <- strata_bound(
    x, n_strata = 3, cv = 0.05, method = "cumrootf", take_all_above = 90
  )

  expect_equal(res$strata$N[3], sum(x >= 90))
  expect_equal(res$strata$n[3], res$strata$N[3])
  expect_lte(res$cv, 0.05 + 1e-10)
})

test_that("take_all CV regression meets its requested precision", {
  res <- strata_bound(
    1:100, n_strata = 3, cv = 0.05,
    method = "cumrootf", take_all_above = 90
  )
  expect_equal(res$strata$N[3], 11L)
  expect_lte(res$cv, 0.05 + 1e-10)
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

test_that("cumrootf: custom n_class respected", {
  res1 <- strata_bound(x_lnorm, n_strata = 3, n = 100,
                        method = "cumrootf", n_class = 50)
  res2 <- strata_bound(x_lnorm, n_strata = 3, n = 100,
                        method = "cumrootf", n_class = 200)
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
                       n_restart = 5L, max_iter = 50L)
  expect_s3_class(res, "svyplan_strata")
  expect_equal(res$method, "kozak")
})

test_that("kozak cv matches target approximately", {
  skip_on_cran()
  target_cv <- 0.10
  res <- strata_bound(x_lnorm, n_strata = 4, cv = target_cv, method = "kozak",
                       n_restart = 10L, max_iter = 100L)
  expect_true(res$cv <= target_cv * 1.5)
})

test_that("kozak n matches target approximately", {
  skip_on_cran()
  target_n <- 200
  res <- strata_bound(x_lnorm, n_strata = 3, n = target_n, method = "kozak",
                       n_restart = 5L, max_iter = 50L)
  expect_true(abs(res$n - target_n) / target_n < 0.5)
})

test_that("proportional allocation: n_h proportional to N", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                       alloc = "proportional")
  df <- res$strata
  prop_alloc <- df$n / sum(df$n)
  prop_pop <- df$N / sum(df$N)
  expect_true(max(abs(prop_alloc - prop_pop)) < 0.15)
})

test_that("neyman allocation: n_h proportional to N * sd", {
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                       alloc = "neyman")
  df <- res$strata
  expected_prop <- df$N * df$sd
  expected_prop <- expected_prop / sum(expected_prop)
  actual_prop <- df$n / sum(df$n)
  expect_true(cor(actual_prop, expected_prop) > 0.8)
})

test_that("power allocation works", {
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                       alloc = "power", alloc_q =0.5)
  expect_s3_class(res, "svyplan_strata")
  expect_equal(res$alloc, "power")
  expect_equal(res$params$alloc_q, 0.5)
})

test_that("n_h >= 2 per stratum", {
  res <- strata_bound(x_lnorm, n_strata = 5, n = 50, method = "cumrootf")
  expect_true(all(res$strata$n >= 2))
})

test_that("take-all stratum works", {
  skip_on_cran()
  thresh <- quantile(x_lnorm, 0.90)
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, take_all_above = thresh)
  expect_s3_class(res, "svyplan_strata")
  expect_true(any(res$strata$take_all))
  take_all_row <- res$strata[res$strata$take_all, ]
  expect_equal(take_all_row$n, take_all_row$N)
})

test_that("output class is svyplan_strata", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_s3_class(res, "svyplan_strata")
  expect_true(inherits(res, "list"))
})

test_that("strata df has correct columns", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expected_cols <- c("stratum", "lower", "upper", "N", "share", "sd",
                     "mean", "n", "take_all")
  expect_equal(names(res$strata), expected_cols)
})

test_that("as.integer returns ceiling(n)", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
})

test_that("as.double returns total sample size", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_identical(as.double(res), as.double(res$n))
  expect_length(res$boundaries, res$n_strata - 1L)
})

test_that("as.data.frame returns strata df", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_identical(as.data.frame(res), res$strata)
})

test_that("strata summaries retain full numerical precision", {
  res <- strata_bound(x_lnorm, n_strata = 4, n = 100, method = "cumrootf")
  bins <- findInterval(x_lnorm, res$boundaries, left.open = TRUE) + 1L
  groups <- split(x_lnorm, factor(bins, levels = seq_len(res$n_strata)))
  expected_sd <- vapply(groups, function(x) if (length(x) < 2L) 0 else sd(x),
                        numeric(1L))
  expected_share <- tabulate(bins, nbins = res$n_strata) / length(x_lnorm)

  expect_equal(res$strata$sd, unname(expected_sd), tolerance = 1e-12)
  expect_identical(res$strata$share, expected_share)
  expect_false(identical(res$strata$sd, round(res$strata$sd, 4L)))
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
  expect_equal(sum(res$strata$N), length(x_lnorm))
})

test_that("sum(n_h) equals n", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf")
  expect_equal(sum(res$strata$n), res$n)
})

test_that("works on simulated lognormal (realistic skewed data)", {
  skip_on_cran()
  set.seed(999)
  x_skew <- rlnorm(2000, meanlog = 8, sdlog = 2)
  res <- strata_bound(x_skew, n_strata = 5, n = 500, method = "kozak",
                       n_restart = 5L, max_iter = 50L)
  expect_s3_class(res, "svyplan_strata")
  expect_equal(nrow(res$strata), 5L)
  expect_true(all(res$strata$N > 0))
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
                          n_restart = 10L, max_iter = 100L)
  expect_true(res_kz$cv <= res_cr$cv * 1.2)
})

test_that("unit_cost parameter is scalar-recycled", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                       unit_cost = 5)
  expect_s3_class(res, "svyplan_strata")
})

test_that("unit_cost parameter with per-stratum vector", {
  res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                       unit_cost = c(1, 2, 5))
  expect_s3_class(res, "svyplan_strata")
})

test_that("params captures expected fields", {
  skip_on_cran()
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "kozak",
                       n_restart = 15L, max_iter = 50L)
  expect_equal(res$params$N, length(x_lnorm))
  expect_equal(res$params$max_iter, 50L)
  expect_equal(res$params$n_restart, 15L)
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

test_that("lh max_iter > 1 improves over max_iter = 1 on skewed data", {
  skip_on_cran()
  set.seed(1021)
  x_skew <- rlnorm(800, meanlog = 6, sdlog = 2)
  res1 <- strata_bound(x_skew, n_strata = 4, n = 200, method = "lh",
                        max_iter = 1L)
  res200 <- strata_bound(x_skew, n_strata = 4, n = 200, method = "lh",
                          max_iter = 200L)
  expect_true(res200$cv <= res1$cv)
})

test_that("power alloc q = 1 matches neyman", {
  res_ney <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                           alloc = "neyman")
  res_pow <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                           alloc = "power", alloc_q =1)
  expect_equal(res_pow$strata$n, res_ney$strata$n)
})

test_that("power alloc q = 0 differs from neyman on skewed data", {
  set.seed(7)
  x_skew <- rlnorm(2000, meanlog = 6, sdlog = 2)
  res_ney <- strata_bound(x_skew, n_strata = 4, n = 400, method = "cumrootf",
                           alloc = "neyman")
  res_pow <- strata_bound(x_skew, n_strata = 4, n = 400, method = "cumrootf",
                           alloc = "power", alloc_q =0)
  expect_false(identical(res_pow$strata$n, res_ney$strata$n))
  expect_equal(res_pow$alloc, "power")
})

test_that("power alloc uses default q = 0.5", {
  res <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                       alloc = "power")
  expect_equal(res$alloc, "power")
  expect_equal(res$params$alloc_q, 0.5)
})

test_that("print shows allocation label", {
  res_ney <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                           alloc = "neyman")
  out_ney <- capture.output(print(res_ney))
  expect_true(any(grepl("allocation: neyman", out_ney)))

  res_pow <- strata_bound(x_lnorm, n_strata = 3, n = 200, method = "cumrootf",
                           alloc = "power", alloc_q =0.3)
  out_pow <- capture.output(print(res_pow))
  expect_true(any(grepl("power \\(alloc_q = 0\\.30\\)", out_pow)))
})

test_that("print is the header the table does not already carry", {
  res <- strata_bound(x_lnorm, n_strata = 4, n = 100)
  out <- capture.output(print(res))
  # Two header lines, the blank that sets the table off, and one row per
  # stratum, plus the column names.
  expect_length(out, 2L + 1L + 1L + 4L)
  expect_lt(max(nchar(out)), 80L)
  # The cut points are the lower/upper columns; naming them again above the
  # table invites a reader to look for a difference that is not there.
  expect_false(any(grepl("^Boundaries:", out)))
  expect_false(any(grepl("^---$", out)))
  # Convergence reads as part of the method that searched for the boundaries.
  expect_match(out[1L], "coordinate search, 4 strata, converged\\)$")
  expect_false(any(grepl("^Converged:", out)))
})

test_that("a non-iterative method claims no convergence either way", {
  # `converged` is NA for cumrootf, and an unmeasured fact must not print as
  # "not converged": isFALSE(NA) is FALSE, which is what keeps it silent.
  res <- strata_bound(x_lnorm, n_strata = 4, method = "cumrootf", n = 100)
  expect_true(is.na(res$converged))
  out <- capture.output(print(res))
  expect_match(out[1L], "\\(Dalenius-Hodges, 4 strata\\)$")
  expect_false(any(grepl("converged", out, ignore.case = TRUE)))
})

test_that("boundaries print at reading precision, never in exponent form", {
  # The frame reaches 3e4 here, which is where a %g format would switch.
  res <- strata_bound(x_lnorm, n_strata = 4, n = 100)
  out <- capture.output(print(res))
  expect_false(any(grepl("e\\+", out)))
  expect_identical(.fmt_boundary(c(8.634646, 428.413168, 30184.6939)),
                   c("8.6346", "428.41", "30185"))
  # Formatted one at a time: a shared format pads every cut point to the
  # widest one's decimals.
  expect_identical(.fmt_boundary(c(1.5, 1000)), c("1.5", "1000"))
})

test_that("alloc field stores method name for all methods", {
  for (a in c("proportional", "neyman", "optimal")) {
    res <- strata_bound(x_unif, n_strata = 3, n = 100, method = "cumrootf",
                         alloc = a)
    expect_equal(res$alloc, a)
    expect_null(res$params$alloc_q)
  }
})

test_that("cumrootf falls back to distinct-value boundaries on discrete data", {
  expect_warning(
    x <- strata_bound(rep(1:4, c(100, 1, 1, 1)), n_strata = 4, n = 20,
                      method = "cumrootf"),
    "adjacent distinct values"
  )
  d <- x$strata
  expect_equal(nrow(d), 4L)
  expect_true(all(d$N >= 1))
  expect_equal(sum(d$N), 103)
  expect_length(x$boundaries, 3L)
  expect_true(all(diff(x$boundaries) > 0))
  expect_equal(sum(d$n), 20L)
})

test_that("unit_cost length must be 1 or n_strata", {
  set.seed(1)
  x <- rlnorm(200)
  expect_error(
    strata_bound(x, n_strata = 3, n = 30, method = "cumrootf",
                 alloc = "optimal", unit_cost = c(1, 2)),
    "length 1 or 3"
  )
})

test_that("cv-mode integer allocation meets the cv target", {
  set.seed(123)
  aux <- rlnorm(100)
  x <- strata_bound(aux, n_strata = 2, cv = 0.211, method = "cumrootf")
  d <- x$strata
  W <- d$N / sum(d$N)
  cv_int <- sqrt(sum(W^2 * d$sd^2 * (1 - d$n / d$N) / d$n)) / mean(aux)
  expect_lte(cv_int, 0.211 * 1.01)
})

test_that("strata cv describes the integer allocation", {
  set.seed(123)
  aux <- rlnorm(100)
  x <- strata_bound(aux, n_strata = 2, cv = 0.211, method = "cumrootf")
  d <- x$strata
  W <- d$N / sum(d$N)
  cv_chk <- sqrt(sum(W^2 * d$sd^2 * (1 - d$n / d$N) / d$n)) / mean(aux)
  expect_equal(x$cv, cv_chk, tolerance = 1e-4)
  expect_lte(x$cv, 0.211 + 1e-10)
  expect_equal(x$params$cv_target, 0.211)
  y <- strata_bound(aux, n_strata = 3, n = 30, method = "cumrootf")
  expect_equal(sum(y$strata$n), 30L)
  expect_true(all(y$strata$n >= 2))
})

test_that("kozak handles highly discrete data without NA crashes", {
  set.seed(7)
  x <- strata_bound(rep(1:6, c(200, 1, 1, 1, 1, 1)), n_strata = 4, n = 30,
                    method = "kozak")
  expect_equal(nrow(x$strata), 4L)
  expect_true(all(x$strata$N >= 1))
  expect_equal(sum(x$strata$n), 30L)
})

test_that("cumrootf errors when distinct values are fewer than strata", {
  expect_error(
    strata_bound(rep(1:3, c(50, 3, 2)), n_strata = 4, n = 20,
                 method = "cumrootf"),
    "unique values"
  )
})

test_that("the strata table is a valid n_alloc frame", {
  set.seed(3)
  x <- rlnorm(4000, 6, 1)
  sb <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh")

  expect_true(all(c("N", "sd", "mean") %in% names(sb$strata)))
  # The means are the stratum means of x, recomputed from the boundaries.
  bins <- cut(x, c(-Inf, sb$boundaries, Inf), labels = FALSE)
  expect_equal(sb$strata$mean, as.numeric(tapply(x, bins, mean)))

  # The handoff runs without the caller reconstructing anything.
  expect_s3_class(n_alloc(sb$strata, cv = 0.05), "svyplan_n")
})

test_that("strata_bound and n_alloc agree on the same design", {
  set.seed(3)
  x <- rlnorm(4000, 6, 1)
  for (spec in list(c(1, 1), c(1.8, 0.85), c(2.5, 0.7))) {
    deff <- spec[1L]
    resp_rate <- spec[2L]
    sb <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh",
                       deff = deff, resp_rate = resp_rate)
    continuous <- .strata_n_for_cv(
      x, sb$boundaries, 0.05, "neyman", 0.5, rep(1, 4),
      deff = deff, resp_rate = resp_rate
    )
    expect_equal(
      continuous,
      n_alloc(sb$strata, cv = 0.05, deff = deff, resp_rate = resp_rate)$n,
      tolerance = 1e-8
    )
  }
})

test_that("deff and resp_rate default to the identity", {
  set.seed(3)
  x <- rlnorm(4000, 6, 1)
  base <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh")
  expect_equal(base$n, 46)
  expect_equal(base$cv, 0.04863648, tolerance = 1e-6)
  expect_equal(
    base$n,
    strata_bound(x, n_strata = 4, cv = 0.05, method = "lh",
                 deff = 1, resp_rate = 1)$n
  )

  # A design effect raises the required sample; a response rate raises it
  # further, because 'n' is what gets fielded.
  expect_gt(strata_bound(x, n_strata = 4, cv = 0.05, deff = 1.8)$n, base$n)
  expect_gt(
    strata_bound(x, n_strata = 4, cv = 0.05, deff = 1.8, resp_rate = 0.85)$n,
    strata_bound(x, n_strata = 4, cv = 0.05, deff = 1.8)$n
  )
  # At a fixed n, they worsen the achieved cv instead.
  expect_gt(
    strata_bound(x, n_strata = 4, n = 400, deff = 2)$cv,
    strata_bound(x, n_strata = 4, n = 400)$cv
  )
})

test_that("a scalar deff leaves the boundaries where they were", {
  set.seed(3)
  x <- rlnorm(4000, 6, 1)
  plain <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh")$boundaries
  scaled <- strata_bound(x, n_strata = 4, cv = 0.05, method = "lh",
                         deff = 1.8, resp_rate = 0.85)$boundaries
  expect_equal(scaled, plain, tolerance = 0.01)
})

test_that("strata_bound accepts a plan", {
  set.seed(3)
  x <- rlnorm(4000, 6, 1)
  plan <- svyplan(deff = 1.8, resp_rate = 0.85, alloc = "power",
                  alloc_q = 0.3)

  from_plan <- strata_bound(x, n_strata = 4, cv = 0.05, plan = plan)
  explicit <- strata_bound(x, n_strata = 4, cv = 0.05, deff = 1.8,
                           resp_rate = 0.85, alloc = "power", alloc_q = 0.3)
  expect_equal(from_plan$boundaries, explicit$boundaries)
  expect_equal(from_plan$n, explicit$n)

  # An explicit argument still wins.
  expect_equal(
    strata_bound(x, n_strata = 4, cv = 0.05, plan = plan, deff = 1)$n,
    strata_bound(x, n_strata = 4, cv = 0.05, resp_rate = 0.85,
                 alloc = "power", alloc_q = 0.3)$n
  )

  # A scalar unit_cost is a design constant and passes through; a vector one
  # is frame-ordered and is rejected rather than silently misapplied.
  expect_s3_class(
    strata_bound(x, n_strata = 4, cv = 0.05, plan = svyplan(unit_cost = 2)),
    "svyplan_strata"
  )
  expect_error(
    strata_bound(x, n_strata = 4, cv = 0.05,
                 plan = svyplan(unit_cost = c(1, 2, 3, 4))),
    "orders costs from the lowest to the highest stratum"
  )
  expect_error(strata_bound(x, n_strata = 4, cv = 0.05, plan = "nope"),
               "must be a svyplan object")
})

test_that("controls the chosen method cannot use are rejected", {
  x <- rlnorm(500, meanlog = 5, sdlog = 1)
  expect_error(
    strata_bound(x, n_strata = 3, n = 100, method = "geo", n_restart = 99),
    "'n_restart' applies only to method"
  )
  expect_error(
    strata_bound(x, n_strata = 3, n = 100, method = "lh", n_class = 20),
    "'n_class' applies only to method"
  )
  expect_error(
    strata_bound(x, n_strata = 3, n = 100, method = "cumrootf", max_iter = 50),
    "'max_iter' applies only to method"
  )
  expect_s3_class(
    strata_bound(x, n_strata = 3, n = 100, method = "kozak",
                 n_restart = 5, max_iter = 50),
    "svyplan_strata"
  )
  expect_s3_class(
    strata_bound(x, n_strata = 3, n = 100, method = "cumrootf", n_class = 20),
    "svyplan_strata"
  )
})

test_that("max_iter defaults to 200 for the iterative methods", {
  x <- rlnorm(300, meanlog = 5, sdlog = 1)
  res <- strata_bound(x, n_strata = 3, n = 100, method = "lh")
  expect_equal(res$params$max_iter, 200)
})
