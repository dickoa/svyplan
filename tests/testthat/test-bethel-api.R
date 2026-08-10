test_that("joint constrained allocation has stable public result semantics", {
  z <- .bethel_fixture()
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 2
  )

  expect_s3_class(fit, "svyplan_n")
  expect_identical(fit$type, "alloc")
  expect_identical(fit$method, "bethel")
  expect_true(all(is.na(c(fit$se, fit$moe, fit$cv))))
  expect_true(all(c("n", "n_int", ".lower", ".upper") %in%
                    names(fit$detail)))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
  expect_true(all(fit$detail$n_int == floor(fit$detail$n_int)))
  expect_true(all(fit$detail$n_int >= ceiling(fit$detail$.lower - 1e-9)))
  expect_true(all(fit$detail$n_int <= floor(fit$detail$.upper + 1e-9)))
  expect_identical(fit$optimization$classification, "optimal")
  expect_equal(as.integer(fit), fit$operational$n)
  expect_equal(as.double(fit), fit$n)
  expect_equal(as.data.frame(fit), fit$detail)
})

test_that("joint constrained allocation dispatch is strict", {
  z <- .bethel_fixture()
  expect_error(
    n_alloc(z$frame, measures = z$measures),
    "supplied together"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets, n = 100),
    "cannot be combined"
  )
  expect_error(
    n_alloc(
      z$frame, measures = z$measures, targets = z$targets,
      domains = "region"
    ),
    "targets\\$domain"
  )
  expect_error(
    n_alloc(
      z$frame, measures = z$measures, targets = z$targets,
      alloc = "neyman"
    ),
    "not used"
  )
  expect_error(
    n_alloc(
      z$frame, measures = z$measures, targets = z$targets,
      alloc_q = 0.3
    ),
    "alloc_q.*not used"
  )
})

test_that("joint allocation honors plan defaults without exposing a method", {
  z <- .bethel_fixture()
  direct <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    deff = 1.5, resp_rate = 0.8, min_n_stratum = 3
  )
  planned <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets,
    plan = svyplan(deff = 1.5, resp_rate = 0.8, min_n_stratum = 3)
  )
  expect_equal(planned$detail$n, direct$detail$n)
})

test_that("joint allocation print and predict methods are deliberate", {
  z <- .bethel_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  output <- capture.output(print(fit))
  expect_true(any(grepl("Joint constrained allocation", output)))
  expect_true(any(grepl("continuous optimum", output)))
  expect_true(any(grepl("field design", output)))
  expect_error(
    predict(fit, data.frame(cv = 0.1)),
    "modify 'targets'"
  )
})

test_that("incomplete generalized multistage contracts are rejected", {
  z <- .bethel_fixture()
  z$frame$n_per_psu <- 10
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets),
    "fixed-take 2-stage.*must contain"
  )
})

test_that("adopted allocations report planning-bound violations", {
  z <- .bethel_fixture()
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  assessed <- prec_alloc(
    fit,
    n = stats::setNames(rep(5, nrow(z$frame)), z$frame$stratum)
  )

  expect_true(any(!assessed$bounds$.pass))
  expect_true(all(assessed$bounds$.lower_violation))
  expect_true(any(grepl(
    "allocation bounds: 4 violated",
    capture.output(print(assessed)),
    fixed = TRUE
  )))
})

## Degrees of freedom in the joint allocation

test_that("an all-cv Bethel problem is bit-identical with and without df", {
  frame <- data.frame(
    stratum = c("a", "b", "c"),
    N = c(4000, 9000, 2500),
    unit_cost = c(12, 9, 15)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("y1", "y2"), each = 3),
    mean = c(0.30, 0.22, 0.41, 55, 62, 48),
    sd = c(0.46, 0.41, 0.49, 20, 24, 18)
  )
  targets <- data.frame(
    name = c("y1", "y2"),
    domain = ".overall",
    level = NA_character_,
    cv = c(0.05, 0.04)
  )

  plain <- n_alloc(frame, measures = measures, targets = targets)
  with_df <- n_alloc(frame, measures = measures, targets = targets, df = 30)

  # a cv-metric constraint sets its variance ceiling with no quantile in it,
  # so nothing the solver sees can move
  expect_identical(with_df$detail$n, plain$detail$n)
  expect_identical(with_df$detail$n_int, plain$detail$n_int)
  expect_identical(with_df$n, plain$n)
  expect_identical(with_df$constraints$.se, plain$constraints$.se)
  expect_identical(with_df$constraints$.cv, plain$constraints$.cv)
  expect_identical(with_df$constraints$.achieved, plain$constraints$.achieved)
  expect_identical(with_df$constraints$.multiplier,
                   plain$constraints$.multiplier)
  expect_identical(with_df$constraints$.sensitivity,
                   plain$constraints$.sensitivity)

  # only the reported margin of error widens, by exactly the quantile ratio
  ratio <- qt(0.975, 30) / qnorm(0.975)
  expect_equal(with_df$constraints$.moe, plain$constraints$.moe * ratio)
  expect_false(identical(with_df$constraints$.moe, plain$constraints$.moe))
})

test_that("a moe-metric constraint reads df in its ceiling and its sensitivity", {
  frame <- data.frame(
    stratum = c("a", "b"),
    N = c(6000, 9000),
    unit_cost = c(10, 10)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 1),
    name = rep("y1", 2),
    mean = c(0.30, 0.22),
    sd = c(0.46, 0.41)
  )
  targets <- data.frame(
    name = "y1", domain = ".overall", level = NA_character_, moe = 0.03
  )

  plain <- n_alloc(frame, measures = measures, targets = targets)
  with_df <- n_alloc(frame, measures = measures, targets = targets, df = 12)

  # a t quantile widens the interval, so holding the same moe costs sample
  expect_gt(with_df$n, plain$n)
  expect_equal(with_df$constraints$.achieved, 0.03, tolerance = 1e-6)

  # the ceiling and the sensitivity reported for it read one quantile: the
  # multiplier must still predict the variance the ceiling admits
  q <- qt(0.975, 12)
  expect_equal(
    with_df$constraints$.sensitivity,
    -2 * with_df$constraints$.multiplier * sum(frame$N)^2 * 0.03 / q^2
  )
  expect_false(isTRUE(all.equal(
    with_df$constraints$.sensitivity, plain$constraints$.sensitivity
  )))
})

test_that("a per-target df column beats the scalar argument", {
  frame <- data.frame(stratum = c("a", "b"), N = c(6000, 9000),
                      unit_cost = c(10, 10))
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("y1", "y2"), each = 2),
    mean = c(0.30, 0.22, 0.5, 0.4),
    sd = c(0.46, 0.41, 0.5, 0.49)
  )
  targets <- data.frame(
    name = c("y1", "y2"), domain = ".overall", level = NA_character_,
    moe = c(0.03, 0.04), df = c(12, NA)
  )
  res <- n_alloc(frame, measures = measures, targets = targets, df = 40)
  # row 1 keeps its own 12, row 2 takes the scalar 40
  per_row <- n_alloc(
    frame, measures = measures,
    targets = transform(targets, df = c(12, 40))
  )
  expect_equal(res$detail$n, per_row$detail$n)
  expect_equal(res$constraints$.moe, per_row$constraints$.moe)
})
