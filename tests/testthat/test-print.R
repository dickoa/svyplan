test_that("print.svyplan_n outputs correctly", {
  result <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
  out <- capture.output(print(result))
  expect_match(out[1], "Sample size for proportion")
  expect_match(out[2], "n = ")
  expect_match(out[2], "deff = 1.50")
})

test_that("print.svyplan_cluster outputs correctly", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  out <- capture.output(print(result))
  expect_match(out[1], "Optimal 2-stage allocation")
  expect_match(out[2], "n_psu")
  expect_match(out[2], "psu_size")
  expect_match(out[3], "cv =")
})

test_that("print.svyplan_varcomp outputs correctly", {
  set.seed(42)
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

test_that("as.integer.svyplan_cluster returns product of ceiled stages", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_equal(as.integer(result), as.integer(prod(ceiling(result$n))))
})

test_that("print returns invisible(x)", {
  result <- n_prop(p = 0.3, moe = 0.05)
  capture.output(expect_invisible(print(result)))
})

test_that("print shows fixed_cost when > 0", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05,
                       fixed_cost = 5000)
  out <- capture.output(print(result))
  expect_true(any(grepl("fixed: 5000", out)))
})

test_that("print hides fixed_cost when 0", {
  result <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  out <- capture.output(print(result))
  expect_false(any(grepl("fixed", out)))
})

test_that("print shows unrounded total for cluster", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  out <- capture.output(print(x))
  expect_true(any(grepl("unrounded", out)))
  expect_true(any(grepl(as.character(floor(x$total_n)), out)))
})

test_that("as.double.svyplan_cluster returns continuous total_n", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_equal(as.double(x), x$total_n)
  expect_true(as.double(x) < as.integer(x))
})

test_that("format.svyplan_cluster shows unrounded", {
  x <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  fmt <- format(x)
  expect_match(fmt, "unrounded")
  expect_match(fmt, as.character(floor(x$total_n)))
})
