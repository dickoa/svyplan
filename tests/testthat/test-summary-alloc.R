summary_alloc_frame <- function() {
  data.frame(
    stratum = c("Urban", "Rural", "Remote"),
    N = c(12000, 30000, 8000),
    sd = c(5, 8, 12),
    mean = c(20, 18, 15),
    unit_cost = c(1, 1.4, 3),
    resp_rate = c(0.90, 0.82, 0.70)
  )
}

summary_alloc_cluster_frame <- function(domains = FALSE) {
  out <- data.frame(
    stratum = c("Urban", "Rural"),
    N = c(50000, 150000),
    sd = c(0.45, 0.48),
    mean = c(0.35, 0.25),
    icc_psu = c(0.03, 0.08),
    cost_psu = c(300, 600),
    cost_ssu = c(40, 60),
    resp_rate = c(0.90, 0.80),
    resp_rate_psu = c(0.95, 0.90),
    N_psu = c(1000, 2000)
  )
  if (domains) out$region <- c("N", "S")
  out
}

test_that("n_alloc summary evaluates the operational field design", {
  x <- n_alloc(summary_alloc_frame(), n = 1200, deff = 1.4)
  s <- summary(x)

  expect_s3_class(s, "summary.svyplan_alloc")
  expect_identical(
    names(s),
    c("kind", "question", "method", "mode", "overall", "continuous",
      "allocation", "precision", "bounds", "domains", "assumptions")
  )
  expect_identical(s$kind, "design")
  expect_identical(s$mode, "n")
  expect_match(s$question, "fixed sample of 1200.*Neyman")
  expect_equal(s$overall$n, x$operational$n)
  expect_equal(s$overall$cost, x$operational$cost)
  expect_equal(s$overall$se, x$operational$se)
  expect_equal(s$overall$moe, x$operational$moe)
  expect_equal(s$overall$cv, x$operational$cv)
  expect_equal(s$continuous$n, x$n)
  expect_equal(s$continuous$se, x$se)
  expect_equal(s$allocation$n_field, x$detail$n_int)
  expect_equal(s$allocation$n_continuous, x$detail$n)
  expect_equal(
    s$overall$expected_respondents,
    sum(x$detail$n_int * c(0.90, 0.82, 0.70))
  )
  expect_equal(sum(s$precision$variance_share), 1)
  expect_equal(s$overall$design_df, as.double(design_df(x)))
  expect_identical(nrow(s$bounds), nrow(x$detail))
  expect_null(s$domains)
  expect_true(all(vapply(s$allocation, function(z) {
    is.numeric(z) || is.character(z)
  }, logical(1L))))
  expect_error(summary(x, digits = 3), "unused argument")

  compact <- capture.output(print(x))
  expect_true(any(grepl("expected respondents = 960.4", compact,
                        fixed = TRUE)))
})

test_that("n_alloc summary print is sectioned and fits an ordinary console", {
  x <- n_alloc(summary_alloc_frame(), n = 1200, deff = 1.4)
  out <- capture.output(expect_invisible(print(summary(x))))

  expect_match(out[1L], "Stratified allocation summary")
  expect_true(any(grepl("Overall field design", out, fixed = TRUE)))
  expect_true(any(grepl("Continuous optimum", out, fixed = TRUE)))
  expect_true(any(grepl("Allocation by stratum", out, fixed = TRUE)))
  expect_true(any(grepl("Cost and weights by stratum", out, fixed = TRUE)))
  expect_true(any(grepl("Achieved precision by stratum", out, fixed = TRUE)))
  expect_true(any(grepl("Allocation bounds", out, fixed = TRUE)))
  expect_true(any(grepl("Expected respondents", out, fixed = TRUE)))
  expect_true(all(nchar(out) <= 80L))
})

test_that("prec_alloc summary evaluates the supplied allocation exactly", {
  frame <- summary_alloc_frame()
  n <- c(171.5, 717.25, 311.75)
  x <- prec_alloc(frame, n = n, deff = 1.4)
  s <- summary(x)

  expect_identical(s$kind, "assessment")
  expect_identical(s$mode, "assessment")
  expect_match(s$question, "assess a supplied")
  expect_null(s$continuous)
  expect_null(s$bounds)
  expect_equal(s$overall$n, sum(n))
  expect_equal(s$overall$se, x$se)
  expect_equal(s$overall$moe, x$moe)
  expect_equal(s$overall$cv, x$cv)
  expect_equal(s$allocation$n_supplied, n)
  expect_false("n_field" %in% names(s$allocation))
  expect_equal(s$precision$se, x$detail$.se)
  expect_equal(s$precision$variance_share, x$detail$.share)

  compact <- capture.output(print(x))
  expect_true(any(grepl("n = 1200.5 (3 strata)", compact, fixed = TRUE)))
  expect_true(any(grepl("expected respondents = 960.7", compact,
                        fixed = TRUE)))

  out <- capture.output(print(s))
  expect_match(out[1L], "Allocation precision summary")
  expect_false(any(grepl("Continuous optimum", out, fixed = TRUE)))
  expect_false(any(grepl("Allocation bounds", out, fixed = TRUE)))
  expect_true(any(grepl("Supplied", out, fixed = TRUE)))
  expect_true(all(nchar(out) <= 80L))
})

test_that("cluster summaries use whole operational stages and stage costs", {
  frame <- summary_alloc_cluster_frame()
  x <- n_alloc(frame, budget = 100000, deff = 1.2)
  s <- summary(x)

  expect_equal(s$allocation$n_field, x$detail$n_int)
  expect_equal(s$allocation$n_psu_field, x$detail$n_psu_int)
  expect_equal(s$allocation$n_per_psu_field, x$detail$n_per_psu_int)
  expect_equal(
    s$allocation$n_field,
    s$allocation$n_psu_field * s$allocation$n_per_psu_field
  )
  expect_equal(s$overall$se, x$operational$se)
  expect_equal(s$overall$moe, x$operational$moe)
  expect_equal(s$overall$cv, x$operational$cv)
  expect_equal(s$overall$cost, x$operational$cost)
  expect_equal(sum(s$allocation$cost), x$operational$cost)
  expect_equal(s$assumptions$response_rate, c(0.855, 0.72))

  out <- capture.output(print(s))
  expect_match(out[1L], "Stratified cluster allocation summary")
  expect_true(any(grepl("Cluster stages by stratum", out, fixed = TRUE)))
  expect_true(any(grepl("PSUs field", out, fixed = TRUE)))
  expect_true(any(grepl("Take field", out, fixed = TRUE)))
  expect_true(all(nchar(out) <= 80L))
})

test_that("domain summaries are recomputed on the field design", {
  frame <- summary_alloc_cluster_frame(domains = TRUE)
  x <- n_alloc(frame, cv = 0.05, domains = "region", deff = c(1.1, 1.3))
  s <- summary(x)

  expect_equal(s$domains$n, x$detail$n_int)
  expect_equal(s$overall$worst_domain_cv, max(s$domains$cv))
  expect_equal(s$overall$worst_domain_cv, x$operational$cv)
  expect_equal(s$continuous$worst_domain_cv, max(x$domains$.cv))

  out <- capture.output(print(s))
  expect_true(any(grepl("Worst domain CV", out, fixed = TRUE)))
  expect_true(any(grepl("Domain allocation", out, fixed = TRUE)))
  expect_true(any(grepl("not additive", out, fixed = TRUE)))
  expect_true(any(grepl("Domain precision", out, fixed = TRUE)))
  expect_false(any(grepl("1_N", out, fixed = TRUE)))
  expect_true(all(nchar(out) <= 80L))
})

test_that("active allocation bounds are named in the summary", {
  frame <- data.frame(
    stratum = c("Small", "Large"),
    N = c(20, 1000),
    sd = c(5, 10),
    mean = c(20, 30),
    take_all = c(TRUE, FALSE)
  )
  s <- summary(n_alloc(frame, n = 100))

  expect_true(s$bounds$binding[s$bounds$stratum == "Small"])
  expect_identical(s$bounds$source[s$bounds$stratum == "Small"], "take_all")
  out <- capture.output(print(s))
  expect_true(any(grepl("1 allocation bound is active", out, fixed = TRUE)))
  expect_true(any(grepl("take_all", out, fixed = TRUE)))
})

test_that("summary fallback remains unchanged for non-allocation results", {
  ordinary <- n_prop(p = 0.3, moe = 0.05)
  expect_s3_class(summary(ordinary), "summaryDefault")
})
