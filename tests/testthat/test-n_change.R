## N1. The closed-form inversion, checked against the family it inverts

test_that("n_change at overlap 0 is twice the single-occasion size", {
  # A difference of two independent means has twice the variance, so each
  # occasion needs twice what one mean needs at the same margin of error.
  res <- n_change(var = 100, moe = 2)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$type, "change")
  expect_equal(res$n, 2 * n_mean(var = 100, moe = 2)$n, tolerance = 1e-12)
})

test_that("prec_change inverts n_change exactly in every target mode", {
  cases <- list(
    n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6),
    n_change(var = 100, change = 5, cv = 0.1, overlap = 0.5, overlap_cor = 0.6),
    n_change(var = 100, change = 5, rmoe = 0.25, overlap = 0.5,
             overlap_cor = 0.6),
    n_change(p = c(0.30, 0.36), moe = 0.02, overlap = 0.75, overlap_cor = 0.5),
    n_change(p = c(0.30, 0.36), cv = 0.15),
    n_change(var = 100, moe = 2, N = 20000, deff = 1.5, resp_rate = 0.8,
             overlap = 0.5, overlap_cor = 0.6),
    n_change(var = 100, moe = 2, ratio = 2, overlap = 0.4, overlap_cor = 0.9),
    n_change(var = c(100, 144), moe = 2, N = 50000, df = 30,
             overlap = 0.6, overlap_cor = 0.7)
  )
  for (res in cases) {
    back <- prec_change(res)
    expect_equal(back$se, res$se, tolerance = 1e-12)
    expect_equal(back$moe, res$moe, tolerance = 1e-12)
    if (!is.na(res$cv)) expect_equal(back$cv, res$cv, tolerance = 1e-12)
  }
})

test_that("each target mode reaches the target it was given", {
  expect_equal(n_change(var = 100, moe = 2)$moe, 2, tolerance = 1e-12)
  expect_equal(n_change(var = 100, change = 5, cv = 0.1)$cv, 0.1,
               tolerance = 1e-12)
  expect_equal(n_change(var = 100, change = 5, rmoe = 0.25)$rmoe, 0.25,
               tolerance = 1e-12)
  expect_equal(n_change(p = c(0.3, 0.36), moe = 0.02)$moe, 0.02,
               tolerance = 1e-12)
})

test_that("n_change round-trips a prec_change result", {
  res <- prec_change(var = 100, n = 500, overlap = 0.5, overlap_cor = 0.6)
  expect_equal(n_change(res)$n, 500, tolerance = 1e-9)

  pair <- prec_change(var = 100, n = c(1000, 500), overlap = 0.4,
                      overlap_cor = 0.7)
  expect_equal(n_change(pair)$n, c(1000, 500), tolerance = 1e-9)
})

## N2. What overlap buys

test_that("a correlated overlap reduces the size, and nothing else does", {
  flat <- n_change(var = 100, moe = 2)$n
  expect_lt(n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6)$n,
            flat)
  expect_equal(n_change(var = 100, moe = 2, overlap = 1, overlap_cor = 0)$n,
               flat, tolerance = 1e-12)
  expect_equal(n_change(var = 100, moe = 2, overlap = 0, overlap_cor = 1)$n,
               flat, tolerance = 1e-12)
})

test_that("the size falls with the product overlap * overlap_cor", {
  size <- function(ov, rho) {
    n_change(var = 100, moe = 2, overlap = ov, overlap_cor = rho)$n
  }
  expect_equal(size(0.8, 0.5), size(0.5, 0.8), tolerance = 1e-12)
  expect_lt(size(0.9, 0.9), size(0.5, 0.5))
})

test_that("no size is identified when the overlap leaves nothing to sample", {
  expect_error(n_change(var = 100, moe = 2, overlap = 1, overlap_cor = 1),
               "no sampling variance at any size")
})

## N3. Finite populations

test_that("the size rises towards N but never past it at equal occasions", {
  sizes <- vapply(
    c(1, 0.1, 0.01, 0.001),
    function(m) n_change(var = 100, moe = m, N = 5000)$n,
    numeric(1L)
  )
  expect_true(all(diff(sizes) > 0))
  expect_true(all(sizes < 5000))
  expect_gt(sizes[4L], 4999)
})

test_that("unequal occasions or nonresponse make a target unattainable", {
  expect_error(n_change(var = 100, moe = 0.01, N = 5000, resp_rate = 0.5),
               "unattainable")
  expect_error(n_change(var = 100, moe = 0.05, N = 5000, ratio = 3),
               "unattainable")
})

## N4. ratio, and the direction overlap is measured in

test_that("ratio sizes the first occasion against the second", {
  res <- n_change(var = 100, moe = 2, ratio = 2)
  expect_length(res$n, 2L)
  expect_equal(res$n[1L] / res$n[2L], 2, tolerance = 1e-12)
})

test_that("ratio 1 returns a single size", {
  expect_length(n_change(var = 100, moe = 2, ratio = 1)$n, 1L)
})

test_that("overlap cannot exceed 1/ratio once the baseline is larger", {
  expect_error(
    n_change(var = 100, moe = 2, ratio = 2, overlap = 0.9, overlap_cor = 0.5),
    "must be <= 1/ratio"
  )
  expect_silent(
    n_change(var = 100, moe = 2, ratio = 2, overlap = 0.4, overlap_cor = 0.5)
  )
})

## N5. deff and response enter where the rest of the package puts them

test_that("n is gross, carrying deff and the 1/resp_rate inflation", {
  base <- n_change(var = 100, moe = 2)$n
  expect_equal(n_change(var = 100, moe = 2, deff = 2)$n, 2 * base,
               tolerance = 1e-12)
  expect_equal(n_change(var = 100, moe = 2, resp_rate = 0.8)$n, base / 0.8,
               tolerance = 1e-12)
})

test_that("deff inflates the variance the fpc then corrects", {
  # The same ordering n_mean documents: inflating a corrected size instead
  # would give a different answer at a non-negligible sampling fraction.
  res <- n_change(var = 100, moe = 2, N = 2000, deff = 2)
  expect_equal(prec_change(res)$moe, 2, tolerance = 1e-12)
  expect_lt(res$n, 2 * n_change(var = 100, moe = 2, N = 2000)$n)
})

## N6. Rejected inputs

test_that("relative targets need a change to be relative to", {
  expect_error(n_change(var = 100, cv = 0.1), "'change' is required")
  expect_error(n_change(var = 100, rmoe = 0.25), "'change' is required")
  expect_error(n_change(p = c(0.3, 0.3), cv = 0.1), "undefined at 'change' = 0")
  expect_error(n_change(p = c(0.3, 0.3), rmoe = 0.25), "undefined at 'change'")
})

test_that("exactly one precision target is required", {
  expect_error(n_change(var = 100), "exactly one of 'moe'")
  expect_error(n_change(var = 100, moe = 2, cv = 0.1), "exactly one of 'moe'")
})

test_that("n_change rejects unused arguments", {
  expect_error(n_change(var = 100, moe = 2, nope = 1), "unused argument")
})

test_that("n_change refuses a prec object of another type", {
  expect_error(n_change(prec_mean(var = 100, n = 400)),
               "requires a svyplan_prec of type 'change'")
  expect_error(prec_change(n_mean(var = 100, moe = 2)),
               "requires a svyplan_n of type 'change'")
})

## N7. Object surface

test_that("the reported target is kept as supplied", {
  expect_equal(n_change(var = 100, moe = 2)$params$moe, 2)
  expect_equal(n_change(var = 100, change = 5, rmoe = 0.25)$params$rmoe, 0.25)
  expect_equal(n_change(var = 100, change = 5, cv = 0.1)$params$cv, 0.1)
  expect_null(n_change(var = 100, moe = 2)$params$cv)
})

test_that("as.integer rounds the gross size up", {
  res <- n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6)
  expect_equal(as.double(res), res$n)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
})

test_that("confint centres on the change the size was built for", {
  res <- n_change(var = 100, moe = 2, change = 5)
  ci <- confint(res)
  expect_equal(unname(ci[1L]), 3, tolerance = 1e-9)
  expect_equal(unname(ci[2L]), 7, tolerance = 1e-9)
})

test_that("print reports the size per occasion and the overlap", {
  out <- capture.output(print(n_change(var = 100, moe = 2, overlap = 0.5,
                                       overlap_cor = 0.6)))
  expect_match(out[1L], "Sample size for change \\(mean scale\\)")
  expect_match(out[2L], "per occasion")
  expect_true(any(grepl("overlap = 0.5", out)))

  pair <- capture.output(print(n_change(var = 100, moe = 2, ratio = 2)))
  expect_match(pair[2L], "then")

  net <- capture.output(print(n_change(var = 100, moe = 2, resp_rate = 0.8)))
  expect_match(net[2L], "net:")
})

test_that("a plan supplies overlap defaults", {
  pl <- svyplan(overlap = 0.5, overlap_cor = 0.6, alpha = 0.10)
  expect_equal(
    n_change(var = 100, moe = 2, plan = pl)$n,
    n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6,
             alpha = 0.10)$n,
    tolerance = 1e-12
  )
})

## N8. Round-trip overrides reach every stored target

test_that("the prec method accepts an rmoe override", {
  x <- prec_change(var = 100, n = 500, change = 5)
  res <- n_change(x, rmoe = 0.2)
  expect_equal(res$rmoe, 0.2, tolerance = 1e-12)
  expect_equal(res$params$rmoe, 0.2)
  expect_null(res$params$moe)
})

test_that("the prec method still implies moe when no target is named", {
  x <- prec_change(var = 100, n = 500, change = 5)
  expect_equal(n_change(x)$n, 500, tolerance = 1e-9)
  expect_equal(n_change(x, cv = 0.05)$cv, 0.05, tolerance = 1e-12)
})

## N9. The proportion scale agrees with n_prop

test_that("n_change at overlap 0 is twice the n_prop size at a finite N", {
  target <- 0.02
  nc <- n_change(p = c(0.3, 0.3), moe = target * sqrt(2), N = 40000)
  np <- n_prop(p = 0.3, moe = target, N = 40000)
  # Two independent occasions carry twice the variance, so a margin wider
  # by sqrt(2) lands on the same size per occasion.
  expect_equal(nc$n, np$n, tolerance = 1e-6)
})

## N10. predict over the design

test_that("predict sweeps overlap and reports one size per occasion", {
  x <- n_change(var = 100, moe = 2, overlap = 0.5, overlap_cor = 0.6)
  grid <- predict(x, data.frame(overlap = c(0, 0.25, 0.5, 0.75)))
  expect_equal(nrow(grid), 4L)
  expect_true(all(c("n1", "n2", "se", "moe") %in% names(grid)))
  expect_true(all(diff(grid$n2) < 0))
  expect_equal(grid$moe, rep(2, 4L), tolerance = 1e-12)
})

test_that("predict keeps the proportion scale and refuses change there", {
  x <- n_change(p = c(0.3, 0.36), moe = 0.02)
  grid <- predict(x, data.frame(overlap_cor = c(0, 0.3, 0.6)))
  expect_equal(nrow(grid), 3L)
  expect_error(predict(x, data.frame(change = c(0.04, 0.06))),
               "change")
})

test_that("predict varies change on the mean scale", {
  x <- n_change(var = 100, change = 5, rmoe = 0.25)
  grid <- predict(x, data.frame(change = c(2, 5, 10)))
  expect_equal(nrow(grid), 3L)
  expect_true(all(diff(grid$n2) < 0))
})
