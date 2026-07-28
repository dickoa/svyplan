test_that("fixed-take continuous and operational precision round-trip", {
  for (stages in 2:3) {
    z <- .bethel_multistage_fixture(stages)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    continuous <- prec_alloc(fit)
    operational <- prec_alloc(fit, n = fit$detail$n_int)

    expect_equal(continuous$detail$.achieved, fit$constraints$.achieved)
    expect_equal(
      operational$detail$.achieved,
      fit$operational$constraints$.achieved
    )
    expect_equal(continuous$params$achieved$n, sum(fit$detail$n))
    expect_equal(operational$params$achieved$n, sum(fit$detail$n_int))
  }
})

test_that("named ultimate-unit allocations convert to whole PSU counts", {
  z <- .bethel_multistage_fixture(3L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  n <- setNames(fit$detail$n_int, fit$detail$stratum)
  direct <- prec_alloc(fit, n = n)
  reordered <- prec_alloc(fit, n = n[rev(names(n))])
  expect_equal(direct$detail$.achieved, reordered$detail$.achieved)

  bad <- n
  bad[1] <- bad[1] + 1
  expect_error(
    prec_alloc(fit, n = bad),
    "whole number of PSUs"
  )
})

test_that("fixed-take precision inverts to the same continuous allocation", {
  for (stages in 2:3) {
    z <- .bethel_multistage_fixture(stages)
    first <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
    second <- n_alloc(prec_alloc(first))
    expect_equal(second$detail$n_psu, first$detail$n_psu, tolerance = 1e-5)
    expect_equal(
      second$params$achieved$cost,
      first$params$achieved$cost,
      tolerance = 1e-6
    )
  }
})
