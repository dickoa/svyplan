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
    "n", "type", "method", "params", "se", "moe", "cv", "targets",
    "detail", "binding", "domains", "operational"
  )
  expect_identical(names(single), expected)
  expect_identical(names(alloc), expected)
  expect_null(single$operational)
  expect_type(alloc$operational, "list")
  expect_true(is.numeric(alloc$se))
  expect_true(is.numeric(alloc$operational$se))
})

test_that("design-effect constructor produces a numeric result class", {
  result <- design_effect(delta = 0.05, psu_size = 20)

  expect_s3_class(result, "svyplan_design_effect")
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_equal(as.double(result), 1.95)
  expect_equal(result * 2, 3.9)
  expect_false(inherits(result * 2, "svyplan_design_effect"))
  expect_equal(sqrt(result), sqrt(1.95))
  expect_false(inherits(sqrt(result), "svyplan_design_effect"))
})
