## Element allocation without a finite population correction

.none_frame <- data.frame(
  stratum = c("A", "B", "C"), N = c(200, 600, 1500), p = c(0.3, 0.4, 0.5)
)

.none_cv <- function(N, p, n, deff = 1, r = 1, e = 1, fpc = TRUE) {
  Ne <- N * e
  W <- Ne / sum(Ne)
  S2 <- if (fpc) Ne * p * (1 - p) / (Ne - 1) else p * (1 - p)
  f <- if (fpc) 1 - n * r / N else 1
  sqrt(sum(W^2 * S2 * deff * f / (n * r * e))) / sum(W * p)
}

test_that("fpc = none gives the uncorrected stratified variance", {
  n <- c(100, 200, 300)
  fr <- .none_frame
  expect_equal(prec_alloc(fr, n = n, fpc = "none")$cv,
               .none_cv(fr$N, fr$p, n, fpc = FALSE), tolerance = 1e-10)
  expect_equal(prec_alloc(fr, n = n)$cv,
               .none_cv(fr$N, fr$p, n), tolerance = 1e-10)

  deff <- c(1.5, 2, 1.2)
  r <- c(0.9, 0.8, 0.85)
  e <- c(0.95, 0.9, 0.8)
  fr$deff <- deff
  fr$resp_rate <- r
  fr$eligibility_rate <- e
  expect_equal(
    prec_alloc(fr, n = n, fpc = "none")$cv,
    .none_cv(fr$N, fr$p, n, deff = deff, r = r, e = e, fpc = FALSE),
    tolerance = 1e-10
  )

  sd_frame <- data.frame(stratum = c("A", "B"), N = c(300, 900),
                         mean = c(10, 12), sd = c(3, 4))
  W <- sd_frame$N / sum(sd_frame$N)
  v <- sum(W^2 * sd_frame$sd^2 / c(50, 80))
  expect_equal(prec_alloc(sd_frame, n = c(50, 80), fpc = "none")$se,
               sqrt(v), tolerance = 1e-10)
})

test_that("fpc = none never reports less variance than fpc = unit", {
  fr <- .none_frame
  for (n in list(c(10, 20, 30), c(100, 200, 300), c(190, 590, 1490))) {
    expect_gt(prec_alloc(fr, n = n, fpc = "none")$cv,
              prec_alloc(fr, n = n)$cv)
  }
  big <- transform(fr, N = N * 1e4)
  expect_equal(prec_alloc(big, n = c(100, 200, 300), fpc = "none")$cv,
               prec_alloc(big, n = c(100, 200, 300))$cv, tolerance = 1e-3)
})

test_that("classic and joint allocation agree under fpc = none", {
  fr <- .none_frame
  measures <- data.frame(stratum = fr$stratum, name = "x", p = fr$p)
  targets <- data.frame(name = "x", domain = ".overall", cv = 0.03)
  classic <- n_alloc(fr, cv = 0.03, fpc = "none")
  joint <- n_alloc(fr[, c("stratum", "N")], measures = measures,
                   targets = targets, fpc = "none")
  expect_equal(joint$detail$n, classic$detail$n, tolerance = 1e-5)
  expect_gt(classic$n, n_alloc(fr, cv = 0.03)$n)

  assessed <- prec_alloc(fr[, c("stratum", "N")], measures = measures,
                         targets = targets, n = c(100, 200, 300),
                         fpc = "none")
  expect_equal(assessed$detail$.achieved,
               .none_cv(fr$N, fr$p, c(100, 200, 300), fpc = FALSE),
               tolerance = 1e-8)
})

test_that("round trips and predict keep fpc = none", {
  fit <- n_alloc(.none_frame, cv = 0.03, fpc = "none")
  expect_identical(fit$params$fpc, "none")
  back <- prec_alloc(fit)
  expect_identical(back$params$fpc, "none")
  expect_equal(back$cv, 0.03, tolerance = 1e-8)
  again <- n_alloc(back)
  expect_identical(again$params$fpc, "none")
  expect_equal(again$n, fit$n, tolerance = 1e-8)
  expect_identical(n_alloc(back, cv = 0.05)$params$fpc, "none")

  grid <- predict(fit, data.frame(cv = c(0.03, 0.05)))
  expect_equal(grid$n[1], fit$n, tolerance = 1e-8)
  expect_equal(grid$n[2], n_alloc(.none_frame, cv = 0.05, fpc = "none")$n,
               tolerance = 1e-8)

  measures <- data.frame(stratum = .none_frame$stratum, name = "x",
                         p = .none_frame$p)
  targets <- data.frame(name = "x", domain = ".overall", cv = 0.03)
  joint <- n_alloc(.none_frame[, c("stratum", "N")], measures = measures,
                   targets = targets, fpc = "none")
  expect_identical(prec_alloc(joint)$params$fpc, "none")
  expect_identical(n_alloc(prec_alloc(joint))$params$fpc, "none")
})

test_that("fpc = none takes a proportion frame with a one-unit stratum", {
  tiny <- data.frame(stratum = c("A", "B"), N = c(1, 500), p = c(0.3, 0.4))
  expect_error(prec_alloc(tiny, n = c(1, 50)), "greater than 1")
  expect_equal(prec_alloc(tiny, n = c(1, 50), fpc = "none")$cv,
               .none_cv(tiny$N, tiny$p, c(1, 50), fpc = FALSE),
               tolerance = 1e-10)
})

test_that("a take-all stratum keeps its sampling variance under fpc = none", {
  fr <- transform(.none_frame, take_all = c(TRUE, FALSE, FALSE))
  unit <- prec_alloc(fr, n = c(200, 200, 300))
  none <- prec_alloc(fr, n = c(200, 200, 300), fpc = "none")
  expect_equal(unit$detail$.se[1], 0)
  expect_equal(none$detail$.se[1], sqrt(0.3 * 0.7 / 200), tolerance = 1e-10)
})

test_that("summary evaluates the field design without the correction", {
  fit <- n_alloc(.none_frame, cv = 0.03, fpc = "none")
  s <- summary(fit)
  n_int <- fit$detail$n_int
  expect_equal(s$overall$cv,
               prec_alloc(.none_frame, n = n_int, fpc = "none")$cv,
               tolerance = 1e-10)
  p <- .none_frame$p
  expect_equal(s$precision$se, sqrt(p * (1 - p) / n_int), tolerance = 1e-10)
})

test_that("summary names the correction away from the default", {
  fr <- .none_frame
  measures <- data.frame(stratum = fr$stratum, name = "x", p = fr$p)
  targets <- data.frame(name = "x", domain = ".overall", cv = 0.03)
  line <- "finite population correction: none"
  expect_true(any(grepl(line, capture.output(
    summary(n_alloc(fr, cv = 0.03, fpc = "none"))
  ), fixed = TRUE)))
  expect_true(any(grepl(line, capture.output(summary(
    n_alloc(fr[, c("stratum", "N")], measures = measures, targets = targets,
            fpc = "none")
  )), fixed = TRUE)))
  expect_false(any(grepl("finite population correction", capture.output(
    summary(n_alloc(fr, cv = 0.03))
  ), fixed = TRUE)))
})

## Strata larger than the integer range

test_that("a stratum above the integer range allocates", {
  huge <- data.frame(stratum = c("a", "b"), N = c(3e9, 1e6), p = c(0.3, 0.4))
  for (fit in list(n_alloc(huge, n = 600),
                   n_alloc(huge, n = 600, alloc = "proportional"),
                   n_alloc(huge, cv = 0.05))) {
    expect_type(fit$detail$n_int, "integer")
    expect_false(anyNA(fit$detail$n_int))
  }
  expect_equal(sum(n_alloc(huge, n = 600)$detail$n_int), 600L)
  expect_equal(prec_alloc(huge, n = c(300, 300))$cv,
               prec_alloc(huge, n = c(300, 300), fpc = "none")$cv,
               tolerance = 1e-3)
  stepped <- transform(huge, n_step = 4)
  expect_equal(n_alloc(stepped, n = 600)$detail$n_int %% 4, c(0, 0))

  measures <- data.frame(stratum = huge$stratum, name = "x", p = huge$p)
  targets <- data.frame(name = "x", domain = ".overall", cv = 0.05)
  joint <- n_alloc(huge[, c("stratum", "N")], measures = measures,
                   targets = targets)
  expect_false(anyNA(joint$detail$n_int))
  expect_equal(joint$detail$n_int, n_alloc(huge, cv = 0.05)$detail$n_int)
})
