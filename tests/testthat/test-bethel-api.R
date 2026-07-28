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
