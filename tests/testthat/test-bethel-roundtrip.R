test_that("continuous and operational precision round-trip exactly", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  continuous <- prec_alloc(fit)
  operational <- prec_alloc(fit, n = fit$detail$n_int)

  expect_s3_class(continuous, "svyplan_prec")
  expect_identical(continuous$method, "bethel")
  expect_equal(continuous$detail$.achieved, fit$constraints$.achieved)
  expect_equal(
    operational$detail$.achieved,
    fit$operational$constraints$.achieved
  )
  expect_equal(as.data.frame(continuous), continuous$detail)
})

test_that("named adopted allocations are matched by stratum", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  n <- setNames(fit$detail$n_int, fit$detail$stratum)
  reordered <- n[rev(names(n))]
  x <- prec_alloc(fit, n = n)
  y <- prec_alloc(fit, n = reordered)

  expect_equal(x$detail$.achieved, y$detail$.achieved)
  expect_error(prec_alloc(fit, n = unname(n[-1])), "one value per stratum")
  bad_names <- n
  names(bad_names)[1] <- "unknown"
  expect_error(prec_alloc(fit, n = bad_names), "match every frame stratum")
})

test_that("generalized precision inverts to the same continuous optimum", {
  z <- .bethel_fixture()
  first <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  precision <- prec_alloc(first)
  second <- n_alloc(precision)

  expect_equal(second$detail$n, first$detail$n, tolerance = 1e-5)
  expect_equal(
    second$params$achieved$cost,
    first$params$achieved$cost,
    tolerance = 1e-6
  )
  expect_identical(
    second$constraints$constraint,
    first$constraints$constraint
  )
})

test_that("explicit default-method generalized assessment works", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  assessed <- prec_alloc(
    z$frame,
    n = fit$detail$n,
    measures = z$measures,
    targets = z$targets
  )
  expect_equal(assessed$detail$.achieved, fit$constraints$.achieved)
})


test_that("a reassessment's alpha and df reach the targets that inherited them", {
  fr <- data.frame(stratum = c("A", "B"), N = c(1000, 2000))
  m <- data.frame(stratum = fr$stratum, name = "y", p = c(0.3, 0.5))
  f <- n_alloc(fr, measures = m, targets = data.frame(name = "y", cv = 0.05),
               df = 30)
  expect_identical(f$params$df, 30)
  expect_true(is.na(f$params$targets$alpha) && is.na(f$params$targets$df))
  e <- prec_alloc(f)
  expect_equal(e$detail$.moe, qt(0.975, 30) * e$detail$.se)
  e <- prec_alloc(f, alpha = 0.1, df = 10)
  expect_equal(e$detail$.moe, qt(0.95, 10) * e$detail$.se)
  # A row's own alpha keeps its precedence over the scalar.
  own <- n_alloc(fr, measures = m,
                 targets = data.frame(name = "y", cv = 0.05, alpha = 0.2))
  e <- prec_alloc(own, alpha = 0.01)
  expect_equal(e$detail$.moe, qnorm(0.9) * e$detail$.se)
  expect_identical(own$params$targets$alpha, 0.2)
})
