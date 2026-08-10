test_that("svyplan_n constructors return canonical fields", {
  single <- n_prop(p = 0.3, moe = 0.05)
  alloc <- n_alloc(
    data.frame(
      stratum = c("a", "b"),
      N = c(1000, 2000),
      sd = c(5, 8)
    ),
    n = 100
  )

  expected <- c(
    "n", "type", "method", "params", "se", "moe", "cv", "rmoe", "indicators",
    "targets", "detail", "binding", "domains", "operational"
  )
  # a proportion carries the one conditional field the schema adds, its
  # expected count of positive cases; an allocation has no proportion to
  # count and stays on the base schema
  expect_identical(names(single), c(expected, "expected_cases"))
  expect_identical(names(alloc), expected)
  expect_equal(single$expected_cases, single$n * 0.3)
  expect_null(alloc$expected_cases)
  expect_null(single$operational)
  expect_type(alloc$operational, "list")
  expect_true(is.numeric(alloc$se))
  expect_true(is.numeric(alloc$operational$se))
})

test_that("design-df constructor produces a numeric result class", {
  result <- design_df(n_psu = 300, n_strata = 20)

  expect_s3_class(result, "svyplan_df")
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_equal(as.double(result), 280)
  expect_equal(result * 2, 560)
  expect_false(inherits(result * 2, "svyplan_df"))
  expect_equal(sqrt(result), sqrt(280))
  expect_false(inherits(sqrt(result), "svyplan_df"))
  expect_identical(
    names(as.list(result)),
    c("df", "n_units", "n_strata", "stage", "strata", "domains")
  )
})

test_that("design-effect constructor produces a numeric result class", {
  result <- design_effect(icc = 0.05, n_per_psu = 20)

  expect_s3_class(result, "svyplan_deff")
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_equal(as.double(result), 1.95)
  expect_equal(result * 2, 3.9)
  expect_false(inherits(result * 2, "svyplan_deff"))
  expect_equal(sqrt(result), sqrt(1.95))
  expect_false(inherits(sqrt(result), "svyplan_deff"))
})
