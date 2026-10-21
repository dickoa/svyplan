## A fixed n allocated from stratum sizes alone, with no variance input

.size_frame <- function() {
  data.frame(stratum = c("A", "B", "C"), N = c(4e5, 1e5, 2.5e4))
}

test_that("power allocation on a size measure needs no placeholder", {
  fr <- .size_frame()
  fr$alloc_measure <- fr$N
  sized <- n_alloc(fr, n = 60, alloc = "power", alloc_q = 0.25,
                   min_n_stratum = 15)
  placeholder <- n_alloc(transform(fr, sd = 1, mean = 1), n = 60,
                         alloc = "power", alloc_q = 0.25, min_n_stratum = 15)
  expect_identical(sized$detail$n_int, c(26L, 19L, 15L))
  expect_identical(sized$detail$n_int, placeholder$detail$n_int)
  expect_equal(sized$detail$n, placeholder$detail$n)
})

test_that("the size measure defaults to N", {
  fr <- .size_frame()
  expect_equal(
    n_alloc(fr, n = 60, alloc = "power", alloc_q = 0.25)$detail$n,
    n_alloc(transform(fr, alloc_measure = N), n = 60, alloc = "power",
            alloc_q = 0.25)$detail$n
  )
})

test_that("the response and design-effect adjustment stays in the weights", {
  fr <- .size_frame()
  rates <- c(0.9, 0.6, 0.8)
  deffs <- c(1.2, 2, 1.5)
  sized <- n_alloc(fr, n = 300, alloc = "power", alloc_q = 0.5,
                   resp_rate = rates, deff = deffs)
  placeholder <- n_alloc(transform(fr, sd = 1, mean = 1), n = 300,
                         alloc = "power", alloc_q = 0.5,
                         resp_rate = rates, deff = deffs)
  expect_equal(sized$detail$n, placeholder$detail$n)
  expect_false(isTRUE(all.equal(
    sized$detail$n,
    n_alloc(fr, n = 300, alloc = "power", alloc_q = 0.5)$detail$n
  )))
})

test_that("proportional allocation needs no variance", {
  fr <- .size_frame()
  prop <- n_alloc(fr, n = 60, alloc = "proportional")
  expect_equal(prop$detail$n, 60 * fr$N / sum(fr$N))
})

test_that("a size-only result reports no precision", {
  x <- n_alloc(.size_frame(), n = 60, alloc = "power", alloc_q = 0.25,
               min_n_stratum = 15)
  expect_identical(c(x$se, x$moe, x$cv), rep(NA_real_, 3))
  expect_true(all(is.na(x$detail[c(".se", ".moe", ".cv", ".share")])))
  expect_identical(c(x$operational$se, x$operational$cv), rep(NA_real_, 2))
  expect_equal(x$operational$cost, 60)

  printed <- capture.output(print(x))
  expect_true(any(
    printed == "precision: not assessed (no 'sd' or 'var' in frame)"
  ))
  expect_false(any(grepl("cv =|se =", printed)))
  summarized <- capture.output(summary(x))
  expect_true(any(grepl("Precision not assessed", summarized)))
  expect_false(any(grepl("Achieved precision by stratum", summarized)))

  swept <- predict(x, newdata = data.frame(n = c(60, 90)))
  expect_true(all(is.na(swept$cv)))
  expect_identical(swept$cost, c(60, 90))
})

test_that("a size-only cluster allocation reports no precision either", {
  cl <- data.frame(
    stratum = c("A", "B"), N = c(4e5, 1e5), N_psu = c(800, 200),
    n_per_psu = 20, icc_psu = 0.05
  )
  x <- n_alloc(cl, n = 600, alloc = "proportional")
  expect_identical(x$detail$n_psu_int, c(24L, 6L))
  expect_true(is.na(x$cv))
  expect_true(is.na(x$operational$cv))
})

test_that("precision targets and variance rules still need a variance", {
  fr <- .size_frame()
  expect_error(n_alloc(fr, cv = 0.05, alloc = "proportional"),
               "'cv' or 'budget' target needs 'sd' or 'var'")
  expect_error(n_alloc(fr, budget = 100, alloc = "power"),
               "'cv' or 'budget' target needs 'sd' or 'var'")
  expect_error(n_alloc(fr, n = 60),
               "alloc = \"neyman\" needs 'sd' or 'var'")
  expect_error(n_alloc(fr, n = 60, alloc = "optimal"),
               "alloc = \"optimal\" needs 'sd' or 'var'")
  expect_error(
    n_alloc(transform(fr, alloc_measure = c(1, -1, 1)), n = 60,
            alloc = "power"),
    "'alloc_measure' must contain one positive finite value per stratum"
  )
})

test_that("assessing a size-only allocation asks for the variance", {
  x <- n_alloc(.size_frame(), n = 60, alloc = "proportional")
  expect_error(prec_alloc(x), "'sd' or 'var'")
})

## Weighting loss of a size-only allocation

test_that("design_effect() of a size-only allocation is Kish's loss on n_int", {
  sizes <- data.frame(stratum = c("A", "B", "C"), N = c(4e5, 1e5, 2.5e4))
  fit <- n_alloc(sizes, n = 60, alloc = "power", alloc_q = 0.25)
  w <- rep(sizes$N / fit$detail$n_int, fit$detail$n_int)
  d <- design_effect(fit)
  expect_equal(as.double(d), as.double(design_effect(weights = w)))
  expect_equal(d$deff_weighting, as.double(d))
  prop <- n_alloc(sizes, n = 60, alloc = "proportional")
  expect_equal(as.double(design_effect(prop)), 1, tolerance = 0.01)
  expect_equal(
    as.double(design_effect(fit, weights = c(1, 1.2, 1.4))),
    as.double(d) * as.double(design_effect(weights = c(1, 1.2, 1.4)))
  )
})
