## A precision print shows the sizes it evaluated

printed_numbers <- function(line) {
  as.numeric(regmatches(line, gregexpr("[0-9]+(\\.[0-9]+)?", line))[[1]])
}

prec_fits <- function() {
  frame <- data.frame(N = c(4000, 3000, 3000), sd = c(10, 15, 8),
                      mean = c(50, 60, 55))
  list(
    prop = list(fit = n_prop(p = 0.3, moe = 0.15), prec = prec_prop),
    mean = list(fit = n_mean(var = 100, moe = 5), prec = prec_mean),
    ratio = list(fit = n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6,
                               component_cor = 0.7, cv = 0.05),
                 prec = prec_ratio),
    change = list(fit = n_change(var = 100, moe = 2, overlap = 0.5,
                                 overlap_cor = 0.6),
                  prec = prec_change),
    pooled = list(fit = n_pooled(var = 100, moe = 2, occasions = 4,
                                 overlap = 0.75, cor_decay = 0.8),
                  prec = prec_pooled),
    alloc = list(fit = n_alloc(frame, cv = 0.01), prec = prec_alloc)
  )
}

test_that("a fit prints its continuous size, labeled, beside its precision", {
  for (case in prec_fits()) {
    fit <- case$fit
    expect_false(all(fit$n == floor(fit$n)))
    pr <- case$prec(fit)
    if (!is.null(fit$moe)) expect_equal(pr$moe, fit$moe, tolerance = 1e-10)
    out <- capture.output(print(pr))
    expect_match(out[1], "continuous optimum\\)$")
    shown <- printed_numbers(out[2])
    expect_equal(shown[1], sum(fit$n), tolerance = 1e-6)
  }
})

test_that("a size override prints the whole size, unlabeled, at its precision", {
  for (case in prec_fits()) {
    fit <- case$fit
    whole <- if (inherits(fit, "svyplan_n") && identical(fit$type, "alloc")) {
      fit$detail$n_int
    } else {
      ceiling(fit$n)
    }
    pr <- case$prec(fit, n = whole)
    expect_null(pr$params$continuous_optimum)
    out <- capture.output(print(pr))
    expect_false(grepl("continuous optimum", out[1]))
    expect_equal(printed_numbers(out[2])[1], sum(whole))
  }
  fit <- n_prop(p = 0.3, moe = 0.15)
  expect_equal(prec_prop(fit, n = 36)$moe, prec_prop(p = 0.3, n = 36)$moe)
  expect_lt(prec_prop(fit, n = 36)$moe, fit$moe)
})

test_that("the allocation header names the design", {
  frame <- data.frame(N = c(4000, 3000, 3000), sd = c(10, 15, 8),
                      mean = c(50, 60, 55))
  fit <- n_alloc(frame, cv = 0.01)
  expect_identical(
    capture.output(print(prec_alloc(fit)))[1],
    "Sampling precision for stratum allocation (continuous optimum)"
  )
  expect_identical(
    capture.output(print(prec_alloc(frame, n = c(200, 150, 100))))[1],
    "Sampling precision for stratum allocation"
  )
})

test_that("an override of the assumptions keeps the label", {
  fit <- n_prop(p = 0.3, moe = 0.15)
  out <- capture.output(print(prec_prop(fit, method = "wilson")))
  expect_identical(
    out[1], "Sampling precision for proportion (wilson, continuous optimum)"
  )
})

test_that("typed fractional sizes print as given, with no label", {
  out <- capture.output(print(prec_prop(p = 0.3, n = 35.3)))
  expect_identical(out[1], "Sampling precision for proportion (wald)")
  expect_identical(out[2], "n = 35.3")
  out <- capture.output(print(prec_prop(p = 0.3, n = 35.3, resp_rate = 0.8)))
  expect_equal(printed_numbers(out[2]), c(35.3, 35.3 * 0.8))
  out <- capture.output(print(prec_mean(var = 100, n = 15.5)))
  expect_identical(out[2], "n = 15.5")
  out <- capture.output(print(prec_change(var = 100, n = 134.5,
                                          resp_rate = 0.8)))
  expect_equal(printed_numbers(out[2])[1:2], c(134.5, 134.5 * 0.8))
  out <- capture.output(print(prec_pooled(var = 100, n = 61.5,
                                          occasions = 4)))
  expect_equal(printed_numbers(out[2])[1:2], c(61.5, 4))
})

test_that("typed whole sizes print as counts, net counts rounded up", {
  out <- capture.output(print(prec_prop(p = 0.3, n = 36, resp_rate = 0.7)))
  expect_identical(out[2], "n = 36 (net: 26)")
  out <- capture.output(print(prec_change(var = 100, n = c(800, 400),
                                          resp_rate = 0.7)))
  expect_match(out[2], "n = 800 then 400 \\(net: 560 then 280\\)")
})

test_that("the change fit print keeps its whole sample to field", {
  fit <- n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6,
                  resp_rate = 0.8)
  out <- capture.output(print(fit))
  expect_match(out[2], sprintf("n = %d per occasion \\(net: %d\\)",
                               ceiling(fit$n), ceiling(fit$n * 0.8)))
})

test_that("the two-phase print shows both phases as evaluated", {
  frame <- data.frame(stratum = c("A", "B"), N = c(6000, 4000),
                      sd = c(12, 25), mean = c(40, 70), unit_cost = c(2, 5),
                      resp_rate = c(0.8, 0.7))
  fit <- n_twophase(frame, cv = 0.02, phase1_cost = 1, resp_rate = 0.9)
  pr <- prec_twophase(fit)
  expect_equal(pr$cv, fit$cv, tolerance = 1e-10)
  out <- capture.output(print(pr))
  expect_identical(out[1],
                   "Sampling precision for two-phase (continuous optimum)")
  shown <- printed_numbers(out[2])
  expect_equal(shown[c(2, 4)], unname(fit$n), tolerance = 1e-6)
  expect_equal(shown[c(5, 6)],
               c(fit$n[["n_phase1"]] * 0.9, sum(pr$detail$n_resp)),
               tolerance = 1e-6)
  frame$nu <- c(0.4, 0.6)
  frame$resp_rate <- NULL
  out <- capture.output(print(prec_twophase(frame, n_phase1 = 2000,
                                            phase1_cost = 1)))
  expect_identical(out[1], "Sampling precision for two-phase")
  expect_identical(out[2], "n_phase1 = 2000 | n_phase2 = 960")
})
