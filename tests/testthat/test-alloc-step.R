## Allocation in steps of n_step

.step_frame <- function(H = 30, seed = 1) {
  set.seed(seed)
  f <- data.frame(
    stratum = sprintf("d%02d", seq_len(H)),
    N = round(runif(H, 2e4, 2e5)),
    n_step = ifelse(seq_len(H) %% 5 == 0, 8, 4)
  )
  f$min_n_stratum <- f$n_step * sample(2:4, H, TRUE)
  f
}

## S1. Bounds in whole steps

test_that("a floor rises to the next step and keeps its source", {
  f <- data.frame(stratum = c("a", "b"), N = c(1000, 2000), n_step = 4,
                  min_n_stratum = c(13, NA))
  d <- n_alloc(f, n = 40, alloc = "proportional")$detail
  expect_equal(d$.lower, c(16, 4))
  expect_identical(d$n_int, c(16L, 24L))
  expect_identical(d$.bound_source[1], "min_n_stratum")
})

test_that("the default floor of one unit becomes one step", {
  f <- data.frame(stratum = c("a", "b"), N = c(100, 1e5), n_step = 4)
  d <- n_alloc(f, n = 40, alloc = "proportional")$detail
  expect_identical(d$n_int, c(4L, 36L))
  expect_identical(d$.bound_source[1], "n_step")
})

test_that("the population caps a stratum at its last whole step", {
  f <- data.frame(stratum = c("a", "b"), N = c(30, 1e5), n_step = 4)
  d <- n_alloc(f, n = 100, alloc = "power", alloc_q = 0)$detail
  expect_equal(d$.upper[1], 28)
  expect_identical(d$n_int, c(28L, 72L))
})

test_that("a floor raised past the cap is refused by name", {
  f <- data.frame(stratum = c("a", "b"), N = c(30, 1e5), n_step = 4,
                  min_n_stratum = c(29, NA))
  expect_error(n_alloc(f, n = 100, alloc = "proportional"),
               "no whole number of steps fits the bounds of stratum a")
})

test_that("a census stratum keeps its whole population", {
  f <- data.frame(stratum = c("a", "b"), N = c(30, 1e5), n_step = 4,
                  take_all = c(TRUE, FALSE))
  d <- n_alloc(f, n = 70, alloc = "proportional")$detail
  expect_identical(d$n_int, c(30L, 40L))
  expect_identical(d$n_step, c(1L, 4L))
})

test_that("a stratum smaller than one step is refused by name", {
  f <- data.frame(stratum = c("a", "b"), N = c(3, 2000), n_step = 4)
  expect_error(n_alloc(f, n = 40, alloc = "proportional"),
               "no whole number of steps fits the bounds of stratum a")
})

test_that("a floor above the population keeps the general message", {
  f <- data.frame(stratum = c("a", "b"), N = c(10, 2000), n_step = 4,
                  min_n_stratum = c(20, NA))
  expect_error(n_alloc(f, n = 40, alloc = "proportional"),
               "lower bounds exceed stratum population")
})

test_that("n_step takes whole numbers of at least one, NA for none", {
  f <- data.frame(stratum = c("a", "b"), N = c(1000, 2000))
  for (bad in list(c(0, 4), c(2.5, 4), c(-4, 4), c(Inf, 4), c("4", "4"))) {
    expect_error(n_alloc(transform(f, n_step = bad), n = 40,
                         alloc = "proportional"),
                 "'n_step' must contain whole numbers >= 1")
  }
  plain <- n_alloc(f, n = 41, alloc = "proportional")
  for (none in list(c(NA, NA), c(1, 1))) {
    x <- n_alloc(transform(f, n_step = none), n = 41, alloc = "proportional")
    expect_identical(x$detail$n_int, plain$detail$n_int)
    expect_null(x$detail$n_step)
  }
  mixed <- n_alloc(transform(f, n_step = c(NA, 4)), n = 41,
                   alloc = "proportional")
  expect_identical(mixed$detail$n_step, c(1L, 4L))
  expect_identical(mixed$detail$n_int %% c(1L, 4L), c(0L, 0L))
})

## S2. The n mode

test_that("steps of 4 and 8 reach the total with every floor met", {
  f <- .step_frame()
  for (alloc in c("proportional", "power")) {
    d <- n_alloc(f, n = 552, alloc = alloc, alloc_q = 0.25)$detail
    expect_identical(sum(d$n_int), 552L)
    expect_true(all(d$n_int %% f$n_step == 0))
    expect_true(all(d$n_int >= f$min_n_stratum))
    expect_true(all(abs(d$n_int - d$n) < d$n_step))
  }
})

test_that("a total the steps cannot reach names the nearest ones", {
  expect_error(
    n_alloc(.step_frame(), n = 550, alloc = "proportional"),
    "reaches the total 550.*nearest totals are 548 and 552"
  )
})

test_that("one step for every stratum equals allocating in sets", {
  set.seed(42)
  for (i in 1:30) {
    H <- sample(3:25, 1)
    k <- sample(c(2, 3, 4, 6, 12), 1)
    N <- round(runif(H, 5e3, 3e5))
    fl <- sample(1:3, H, TRUE)
    K <- sum(fl) + sample(5:60, 1)
    alloc <- sample(c("proportional", "power"), 1)
    stepped <- n_alloc(
      data.frame(stratum = seq_len(H), N = N, n_step = k,
                 min_n_stratum = k * fl),
      n = K * k, alloc = alloc, alloc_q = 0.3
    )$detail$n_int
    sets <- n_alloc(
      data.frame(stratum = seq_len(H), N = N, min_n_stratum = fl),
      n = K, alloc = alloc, alloc_q = 0.3
    )$detail$n_int
    expect_identical(stepped, as.integer(k * sets), info = paste("frame", i))
  }
})

## S3. Rule-rounded cv and budget modes

.step_variance_frame <- function() {
  f <- .step_frame()
  set.seed(7)
  f$p <- runif(nrow(f), 0.1, 0.5)
  f$unit_cost <- runif(nrow(f), 1, 3)
  f$region <- rep(c("A", "B", "C"), length.out = nrow(f))
  f
}

test_that("a rule-rounded cv target is met in whole steps", {
  f <- .step_variance_frame()
  x <- n_alloc(f, cv = 0.05, domains = "region")
  expect_true(all(x$detail$n_int %% f$n_step == 0))
  dom <- prec_alloc(f, n = x$detail$n_int, domains = "region")$domains
  expect_true(all(dom$.cv <= 0.05 + 1e-12))
})

test_that("a rule-rounded budget is never exceeded in whole steps", {
  f <- .step_variance_frame()
  x <- n_alloc(f, budget = 3000, alloc = "neyman")
  expect_true(all(x$detail$n_int %% f$n_step == 0))
  expect_lte(sum(x$detail$n_int * f$unit_cost), 3000)
  spare <- 3000 - sum(x$detail$n_int * f$unit_cost)
  room <- x$detail$n_int + f$n_step <= x$detail$.upper
  expect_true(all(spare < (f$unit_cost * f$n_step)[room]))
})

## S4. Searched modes

test_that("the searched cv and budget modes return whole steps", {
  f <- .step_variance_frame()
  f$region <- NULL
  cv_fit <- n_alloc(f, cv = 0.02, alloc = "neyman")
  expect_true(all(cv_fit$detail$n_int %% f$n_step == 0))
  expect_lte(cv_fit$operational$cv, 0.02)
  budget_fit <- n_alloc(f, budget = 3000, alloc = "optimal")
  expect_true(all(budget_fit$detail$n_int %% f$n_step == 0))
  expect_lte(sum(budget_fit$detail$n_int * f$unit_cost), 3000)
})

test_that("the searched cv mode is the smallest design in steps", {
  f <- data.frame(stratum = c("a", "b", "c"), N = c(4000, 2500, 900),
                  p = c(0.2, 0.45, 0.3), n_step = c(4, 6, 8))
  target <- 0.06
  fit <- n_alloc(f, cv = target, alloc = "neyman")
  grid <- expand.grid(a = seq(4, 400, 4), b = seq(6, 402, 6),
                      c = seq(8, 400, 8))
  grid <- grid[rowSums(grid) <= sum(fit$detail$n_int), ]
  cv <- apply(grid, 1, function(n) prec_alloc(f, n = n)$cv)
  best <- min(rowSums(grid)[cv <= target])
  expect_identical(sum(fit$detail$n_int), as.integer(best))
})

## S5. Refusals, detail, print and round trips

test_that("a two-stage frame and joint allocation refuse n_step", {
  cl <- data.frame(stratum = c("a", "b"), N = c(1000, 2000), sd = 1,
                   icc_psu = 0.1, n_per_psu = 10, n_step = 4)
  expect_error(n_alloc(cl, n = 40), "'n_step' applies to element frames")
  expect_error(prec_alloc(cl, n = c(20, 20)),
               "'n_step' applies to element frames")
  jf <- data.frame(stratum = c("a", "b"), N = c(1000, 2000), n_step = 4)
  m <- data.frame(stratum = c("a", "b"), name = "p", p = 0.3)
  t <- data.frame(name = "p", cv = 0.1)
  expect_error(n_alloc(jf, measures = m, targets = t),
               "'n_step' is not used for joint constrained allocation")
  expect_error(prec_alloc(jf, n = c(20, 20), measures = m, targets = t),
               "'n_step' is not used for joint constrained allocation")
})

test_that("prec_alloc() assesses any n on a frame with steps", {
  f <- transform(.step_variance_frame()[1:3, ], region = NULL)
  expect_equal(
    prec_alloc(f, n = c(13, 21, 30))$cv,
    prec_alloc(transform(f, n_step = NULL), n = c(13, 21, 30))$cv
  )
})

test_that("the steps reach the detail, the print and the round trips", {
  f <- .step_variance_frame()
  f$region <- NULL
  x <- n_alloc(f, n = 552)
  expect_identical(x$detail$n_step, as.integer(f$n_step))
  expect_output(print(x), "field sizes in steps of 4, 8")
  expect_equal(prec_alloc(x, n = x$detail$n_int)$cv, x$operational$cv)
  again <- n_alloc(prec_alloc(x), n = 552)
  expect_identical(again$detail$n_int, x$detail$n_int)
  expect_no_warning(predict(x, data.frame(n = c(552, 600))))
})
