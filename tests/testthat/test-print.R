test_that("print.svyplan_n outputs correctly", {
  result <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
  out <- capture.output(print(result))
  expect_match(out[1], "Sample size for proportion")
  expect_match(out[2], "n = ")
  expect_match(out[2], "deff = 1.50")
})

test_that("print.svyplan_cluster outputs correctly", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  out <- capture.output(print(result))
  expect_match(out[1], "Optimal 2-stage allocation")
  expect_match(out[2], "n_psu")
  expect_match(out[2], "n_per_psu")
  expect_match(out[3], "cv =")
})

test_that("print.svyplan_varcomp outputs correctly", {
  set.seed(1019)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  result <- varcomp(income ~ district, data = frame)
  out <- capture.output(print(result))
  expect_match(out[1], "Variance components")
  expect_match(out[2], "varb =")
})

test_that("format methods return strings", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_true(is.character(format(result)))
  expect_match(format(result), "svyplan_n")
})

test_that("as.integer.svyplan_n returns ceiling", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(as.integer(result), ceiling(result$n))
})

test_that("as.double.svyplan_n returns raw n", {
  result <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(as.double(result), result$n)
})

test_that("as.integer.svyplan_cluster returns the operational stage vector", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_identical(as.integer(result), as.integer(result$operational$n))
  expect_length(as.integer(result), 2L)
  expect_lte(result$operational$cost, 100000)
})

test_that("print returns invisible(x)", {
  result <- n_prop(p = 0.3, moe = 0.05)
  capture.output(expect_invisible(print(result)))
})

test_that("print shows fixed_cost when > 0", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                       fixed_cost = 5000)
  out <- capture.output(print(result))
  expect_true(any(grepl("fixed: 5000", out)))
})

test_that("print hides fixed_cost when 0", {
  result <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  out <- capture.output(print(result))
  expect_false(any(grepl("fixed", out)))
})

test_that("print shows field design and continuous optimum for cluster", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  out <- capture.output(print(x))
  expect_true(any(grepl("field design", out)))
  expect_true(any(grepl("continuous optimum", out)))
  expect_true(any(grepl(sprintf("total n = %d", x$operational$total_n), out)))
})

test_that("as.double.svyplan_cluster returns the continuous stage vector", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  expect_equal(as.double(x), x$n)
  expect_length(as.double(x), length(as.integer(x)))
})

test_that("format.svyplan_cluster shows unrounded", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  fmt <- format(x)
  expected_unrounded <- format(signif(x$total_n, 8), trim = TRUE, scientific = FALSE)
  expect_match(fmt, "unrounded")
  expect_true(grepl(expected_unrounded, fmt, fixed = TRUE))
})

test_that("cluster print and format handle pathological stage sizes", {
  x <- structure(
    list(
      n = c(n_psu = 3.99303258914344e-15, n_per_psu = 150518950596651648),
      stages = 2L,
      total_n = 601.027075016102,
      se = NA_real_,
      moe = NA_real_,
      cv = 0.2,
      cost = 30051.3537508051,
      params = list(resp_rate = 1),
      targets = NULL,
      detail = NULL,
      binding = NULL,
      domains = NULL
    ),
    class = c("svyplan_cluster", "list")
  )

  expect_no_error(capture.output(print(x)))
  expect_no_error(format(x))
})

test_that("print.svyplan_power shows vector n", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  out <- capture.output(print(res))
  expect_true(any(grepl("n1 = ", out)))
  expect_true(any(grepl("n2 = ", out)))
  expect_true(any(grepl("total = ", out)))
})

test_that("print.svyplan_power shows method for arcsine", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, method = "arcsine")
  out <- capture.output(print(res))
  expect_true(any(grepl("method = arcsine", out)))
})

test_that("print.svyplan_power shows one-sided", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, alternative = "one.sided")
  out <- capture.output(print(res))
  expect_true(any(grepl("one-sided", out)))
})

test_that("format.svyplan_power handles vector n", {
  res <- power_mean(effect = 5, var = 100, ratio = 2)
  fmt <- format(res)
  expect_match(fmt, ",")
})

test_that("as.data.frame.svyplan_n returns the alloc detail table", {
  frame <- data.frame(
    stratum = c("A", "B"),
    N = c(100, 200),
    sd = c(5, 10)
  )
  res <- n_alloc(frame, n = 50, alloc = "neyman")
  df <- as.data.frame(res)
  expect_identical(df, res$detail)
  expect_true(all(c("stratum", "n", "n_int", "weight") %in% names(df)))
})

test_that("as.data.frame.svyplan_n returns domains for n_multi", {
  targets <- data.frame(
    name = "stunting",
    domain = c("urban", "rural"),
    p = c(0.25, 0.35),
    cv = 0.08
  )
  res <- n_multi(targets, domains = "domain")
  expect_identical(as.data.frame(res), res$domains)

  no_dom <- n_multi(data.frame(name = "x", p = 0.3, cv = 0.08))
  expect_identical(as.data.frame(no_dom), no_dom$detail)
})

test_that("as.data.frame.svyplan_n returns one-row summary otherwise", {
  res <- n_prop(p = 0.3, moe = 0.05)
  df <- as.data.frame(res)
  expect_equal(nrow(df), 1L)
  expect_equal(df$n, res$n)
  expect_equal(df$n_int, as.integer(res))
})

test_that("as.data.frame.svyplan_cluster returns the stage table", {
  res <- n_cluster(budget = 100000, icc = 0.05, unit_relvar = 1,
                   stage_cost = c(500, 50))
  df <- as.data.frame(res)
  expect_equal(df$stage, c("n_psu", "n_per_psu"))
  expect_equal(df$n, as.numeric(res$n))
  expect_identical(df$n_int, as.integer(res))
  expect_false(identical(as.integer(ceiling(res$n)), as.integer(res)))
})

test_that("as.data.frame.svyplan_cluster returns domains when present", {
  targets <- data.frame(
    name = "s",
    domain = c("u", "r"),
    p = c(0.25, 0.35),
    cv = 0.08,
    icc_psu = 0.05
  )
  res <- n_multi_cluster(targets, domains = "domain", stage_cost = c(500, 50))
  expect_identical(as.data.frame(res), res$domains)
})

test_that("as.data.frame.svyplan_prec returns a single-indicator summary", {
  res <- prec_prop(p = 0.3, n = 400)
  df <- as.data.frame(res)

  expect_identical(names(df), c("n", "se", "moe", "cv"))
  expect_equal(nrow(df), 1L)
  expect_equal(df$n, 400)
  expect_equal(df$se, res$se)
})

test_that("as.data.frame.svyplan_prec returns cluster stages", {
  res <- prec_cluster(
    n = c(n_psu = 30, n_per_psu = 10),
    icc = 0.05,
    unit_relvar = 1
  )
  df <- as.data.frame(res)

  expect_identical(
    names(df),
    c("n_psu", "n_per_psu", "total_n", "se", "moe", "cv")
  )
  expect_equal(df$total_n, 300)
  expect_equal(df$cv, res$cv)
})

test_that("as.data.frame.svyplan_prec returns detail when available", {
  targets <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.4),
    n = c(400, 500)
  )
  res <- prec_multi(targets)

  expect_identical(as.data.frame(res), res$detail)
})

test_that("as.data.frame.svyplan_power has a stable two-group schema", {
  equal <- power_prop(p1 = 0.30, p2 = 0.35)
  unequal <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  expected_names <- c(
    "n1", "n2", "n1_int", "n2_int", "power", "effect", "type", "solved"
  )

  equal_df <- as.data.frame(equal)
  unequal_df <- as.data.frame(unequal)
  expect_identical(names(equal_df), expected_names)
  expect_identical(names(unequal_df), expected_names)
  expect_equal(equal_df$n1, equal_df$n2)
  expect_equal(c(unequal_df$n1, unequal_df$n2), unname(unequal$n))
})

test_that("a single binding constraint prints on one line with the total", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  expect_equal(sum(fit$constraints$.binding), 1L)
  expect_lt(sum(fit$constraints$.binding), nrow(fit$constraints))

  out <- capture.output(print(fit))
  binding <- fit$constraints$constraint[fit$constraints$.binding]
  expect_true(any(grepl(paste0("binding: ", binding), out, fixed = TRUE)))
  # the constraint total is on the field-design line, so the reader can see
  # that one of several is named without a separate see-also line
  expect_true(any(grepl(
    sprintf("%d targets, all pass", nrow(fit$constraints)), out, fixed = TRUE
  )))
  # confirmations that never vary are gone
  expect_false(any(grepl("status: optimal", out, fixed = TRUE)))
  expect_false(any(grepl("active bounds: 0 lower, 0 upper", out, fixed = TRUE)))

  prec_out <- capture.output(print(prec_alloc(fit)))
  expect_true(any(grepl("binding of", prec_out, fixed = TRUE)))
  expect_true(any(grepl("see $detail", prec_out, fixed = TRUE)))
})

test_that("active bounds are reported only when some are active", {
  z <- .bethel_fixture()
  bounded <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    min_n_stratum = min(z$frame$N)
  )
  expect_true(any(grepl("active bounds:", capture.output(print(bounded)),
                        fixed = TRUE)))
})

test_that("print omits the subset label when every constraint is shown", {
  z <- .bethel_fixture()
  z$targets <- z$targets[1L, , drop = FALSE]
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  expect_identical(nrow(fit$constraints), 1L)

  out <- capture.output(print(fit))
  expect_false(any(grepl("showing", out, fixed = TRUE)))
  expect_false(any(grepl("see $constraints", out, fixed = TRUE)))
})

## Probabilities never print as a bound they are not

test_that("a probability that rounds onto a bound keeps its own digits", {
  # Two decimals wherever they read well, and more where rounding would name
  # a value the validators refuse.
  expect_identical(svyplan:::.fmt_prob(0.95), "0.95")
  expect_identical(svyplan:::.fmt_prob(0.8), "0.80")
  expect_identical(svyplan:::.fmt_prob(0.05), "0.05")
  expect_identical(svyplan:::.fmt_prob(0.999), "0.999")
  expect_identical(svyplan:::.fmt_prob(0.9999999), "0.9999999")
  expect_identical(svyplan:::.fmt_prob(0.001), "0.001")
  # The bounds themselves are legal for some arguments and print plainly.
  expect_identical(svyplan:::.fmt_prob(1), "1.00")
  expect_identical(svyplan:::.fmt_prob(0), "0.00")
})

test_that("a response rate below 1 never prints as 1", {
  # The line appears only because the rate is below 1, so showing 1.00 would
  # contradict the test that produced it.
  out <- capture.output(print(n_prop(0.3, moe = 0.05, resp_rate = 0.999)))
  expect_true(any(grepl("resp_rate = 0.999", out, fixed = TRUE)))
  expect_false(any(grepl("resp_rate = 1.00", out, fixed = TRUE)))
  expect_true(any(grepl(
    "resp_rate = 0.80",
    capture.output(print(n_prop(0.3, moe = 0.05, resp_rate = 0.8))),
    fixed = TRUE
  )))
})

test_that("a small alpha prints as itself", {
  out <- capture.output(print(power_mean(var = 100, effect = 2, alpha = 0.001)))
  expect_true(any(grepl("alpha = 0.001", out, fixed = TRUE)))
  expect_false(any(grepl("alpha = 0.00,", out, fixed = TRUE)))
})
