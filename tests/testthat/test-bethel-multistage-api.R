test_that("two-stage generalized allocation exposes both public units", {
  z <- .bethel_multistage_fixture(2L)
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  d <- fit$detail

  expect_s3_class(fit, "svyplan_n")
  expect_identical(fit$method, "bethel")
  expect_identical(fit$params$stages, 2L)
  expect_true(all(c("n_psu", "n_psu_int", "n_per_psu") %in% names(d)))
  expect_equal(d$n, d$n_psu * d$n_per_psu)
  expect_equal(d$n_int, d$n_psu_int * d$n_per_psu)
  expect_type(d$n_int, "double")
  expect_type(d$n_psu_int, "double")
  expect_true(all(d$n_psu_int == floor(d$n_psu_int)))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
  expect_equal(
    fit$params$achieved$cost,
    sum(d$n_psu * (d$cost_psu + d$cost_ssu * d$n_per_psu))
  )
})

test_that("three-stage generalized allocation preserves fixed-take identities", {
  z <- .bethel_multistage_fixture(3L)
  fit <- n_alloc(
    z$frame, measures = z$measures, targets = z$targets, min_n_stratum = 20
  )
  d <- fit$detail

  expect_identical(fit$params$stages, 3L)
  expect_true(all(c("N_ssu", "n_per_ssu", "cost_tsu") %in% names(d)))
  expect_equal(d$n, d$n_psu * d$n_per_psu * d$n_per_ssu)
  expect_equal(d$n_int, d$n_psu_int * d$n_per_psu * d$n_per_ssu)
  expect_type(d$n_int, "double")
  expect_type(d$n_psu_int, "double")
  expect_equal(fit$n, sum(d$n))
  expect_equal(fit$operational$n, sum(d$n_int))
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
})
