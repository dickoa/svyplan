## The rotation is where a pattern is parsed
##
## Every notation, every refusal and every derived property of a unit's life
## belong to design_rotation(). design_overlap() and design_schedule() read
## the object it returns, so a pattern is parsed once and the two cannot
## disagree about what a string means.

## R1. The two notations, and the string that means both

test_that("an even spell count is refused rather than read as a cycle", {
  expect_error(design_rotation("4-8"), "repeating cycle")
  expect_error(design_rotation("4-8-4-8"), "repeating cycle")
  expect_silent(design_rotation("4-8-4"))
  expect_silent(design_rotation("4-8-4-8-4"))
})

test_that("a rotation must start and end in sample", {
  expect_error(design_rotation(c(0, 1, 1)), "start and end in sample")
  expect_error(design_rotation(c(1, 1, 0)), "start and end in sample")
})

test_that("a rotation must be a spec or a finite numeric vector", {
  expect_error(design_rotation(list(1, 2)), "compact string")
  expect_error(design_rotation(c(1, NA, 1)), "compact string")
  expect_error(design_rotation(c(1, Inf, 1)), "compact string")
  expect_error(design_rotation(numeric(0)), "compact string")
  expect_error(design_rotation(c(1, -1, 1)),
               "must not take a negative number")
})

test_that("a malformed compact spec names itself", {
  expect_error(design_rotation("4-x-4"), "is not a rotation spec")
  expect_error(design_rotation("4.5-8-4"), "is not a rotation spec")
  # a 0 makes it a per-occasion pattern, where 8 is not a flag
  expect_error(design_rotation("0-8-4"), "other than 0 or 1")
  expect_error(design_rotation(""), "single non-empty string")
  expect_error(design_rotation(c("4", "5")), "single non-empty string")
})

test_that("a life of one occasion has no lag and is refused", {
  expect_error(design_rotation("1"), "at least two occasions")
  expect_error(design_rotation(1), "at least two occasions")
  expect_error(design_rotation(TRUE), "at least two occasions")
  expect_silent(design_rotation("2"))
})

test_that("no wave may take more units than the cohort recruited", {
  expect_error(design_rotation(c(0.5, 1)), "more units at a later occasion")
  expect_error(design_rotation(c(1, 2, 1)), "more units at a later occasion")
  # An interim wave subsampled and then the whole cohort again is nested and
  # stays accepted: the take never exceeds the recruitment.
  r <- design_rotation(c(1, 0.5, 1))
  expect_equal(attr(r, "n_occasion"), 2.5, tolerance = 1e-12)
  expect_false(attr(r, "equal_take"))
  expect_true(attr(r, "unbroken"))
})

test_that("an all-1s spec is refused, naming both readings", {
  # "1-1-1" is three consecutive occasions in one notation and in-out-in in
  # the other, and the two differ at every lag
  expect_error(design_rotation("1-1-1"), "is ambiguous")
  expect_error(design_rotation("1-1-1"), "consecutive occasions in sample")
  expect_error(design_rotation("1-1-1"), "\"3\"")
  expect_error(design_rotation("1-1-1"), "\"1-0-1\"")
  # the two readings it names are the two designs, and neither is silent
  expect_equal(as.double(design_rotation("3")), rep(1, 3), tolerance = 1e-12)
  expect_equal(as.double(design_rotation("1-0-1")), c(1, 0, 1),
               tolerance = 1e-12)
})

test_that("an even all-1s spec offers the pattern reading only", {
  # Lynn Figures 4 and 5: six consecutive years
  expect_error(design_rotation("1-1-1-1-1-1"), "is ambiguous")
  expect_error(design_rotation("1-1-1-1-1-1"), "\"6\"")
  expect_error(design_rotation("1-1-1-1-1-1"), "repeating cycle")
  expect_equal(as.double(design_rotation("6")), rep(1, 6), tolerance = 1e-12)
})

test_that("a pattern is validated as a life like any other rotation", {
  expect_error(design_rotation("0-1-1"), "must start and end in sample")
  expect_error(design_rotation("1-1-0"), "must start and end in sample")
  expect_error(design_rotation("2-2-0"), "other than 0 or 1")
  expect_error(design_rotation("1-0.5-1"), "is not a rotation spec")
})

## R2. The object surface

test_that("a rotation carries the life, the take and the two shape flags", {
  cps <- design_rotation("4-8-4")
  expect_s3_class(cps, "svyplan_rotation")
  expect_equal(as.double(cps), c(rep(1, 4), rep(0, 8), rep(1, 4)))
  expect_identical(attr(cps, "life"), 16L)
  expect_equal(attr(cps, "n_occasion"), 8)
  expect_false(attr(cps, "unbroken"))
  expect_true(attr(cps, "equal_take"))
  expect_identical(attr(cps, "spec"), "4-8-4")

  panel <- design_rotation("5")
  expect_true(attr(panel, "unbroken"))
  expect_true(attr(panel, "equal_take"))
})

test_that("as.data.frame reports one row per occasion of the life", {
  d <- as.data.frame(design_rotation("1-1-0-0-1-1"))
  expect_named(d, c("occasion", "in_sample", "take"))
  expect_equal(d$occasion, 1:6)
  expect_equal(d$in_sample, c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(d$take, c(1, 1, 0, 0, 1, 1))
})

test_that("print and format name the spells", {
  out <- capture.output(print(design_rotation("4-8-4")))
  expect_true(any(grepl("16-occasion life", out)))
  expect_true(any(grepl("8 in sample each occasion", out)))
  expect_true(any(grepl("4 in, 8 out, 4 in", out)))
  expect_match(format(design_rotation("4-8-4")), "^svyplan_rotation \\[")
})

test_that("the methods reject unused arguments", {
  cps <- design_rotation("4-8-4")
  expect_error(print(cps, nope = 1), "unused argument")
  expect_error(as.double(cps, nope = 1), "unused argument")
  expect_error(as.data.frame(cps, nope = 1), "unused argument")
})

## R3. Construction is idempotent, and an overlap is not a rotation

test_that("a rotation passed back in is returned unchanged", {
  cps <- design_rotation("4-8-4")
  expect_identical(design_rotation(cps), cps)
})

test_that("an overlap handed in as a rotation names the fix", {
  expect_error(design_rotation(design_overlap("4")), "not a rotation")
})

## R4. design_overlap() reads the rotation, and the string is sugar for it

test_that("wrapping the input first leaves the overlap unchanged", {
  grid <- list(
    "4-8-4", "5", "3", "27", "1-1-0-0-1-1", "2-2-2", "1-0-1",
    c(1, 1, 0, 1), c(1, 0.5, 1), c(1, 1, 0.5, 0.5), c(TRUE, TRUE, FALSE, TRUE)
  )
  for (spec in grid) {
    label <- paste(format(spec), collapse = ",")
    expect_identical(
      as.double(design_overlap(design_rotation(spec))),
      as.double(design_overlap(spec)),
      info = label
    )
  }
})

test_that("the overlap carries the rotation it was computed from", {
  cps <- design_overlap("4-8-4")
  expect_s3_class(cps$rotation, "svyplan_rotation")
  expect_identical(cps$rotation, design_rotation("4-8-4"))
  expect_equal(cps$life, 16L)
  expect_equal(cps$n_occasion, 8)
})

test_that("the old argument name is refused rather than matched", {
  expect_error(design_overlap(schedule = "4"), "unused argument")
})
