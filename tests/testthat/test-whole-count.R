## A size solved back from a whole size's precision fields that whole size

# Invariant: a precision computed at a whole size k is met by k, so sizing
# the same target again must report k, not the k + 1 a raw ceiling() gives
# when the inverse lands a rounding step above k.

rotation <- design_overlap(c(1, 0, 1), max_lag = 1)
round_trips <- list(
  n_prop = function(k) n_prop(prec_prop(p = 0.3, n = k, N = 5000)),
  n_prop_inf = function(k) n_prop(prec_prop(p = 0.3, n = k)),
  n_prop_rates = function(k) {
    n_prop(prec_prop(p = 0.3, n = k, N = 5000, resp_rate = 0.7, deff = 1.3))
  },
  n_mean = function(k) {
    n_mean(prec_mean(var = 100, mu = 50, n = k, N = 5000))
  },
  n_ratio = function(k) {
    n_ratio(prec_ratio(r = 0.5, cv_num = 0.4, cv_den = 0.3,
                       component_cor = 0.5, n = k, N = 5000))
  },
  n_change = function(k) {
    n_change(prec_change(var = 100, n = k, N = 5000, overlap = 0.5,
                         overlap_cor = 0.8))
  },
  n_change_tight = function(k) {
    n_change(prec_change(var = 100, n = k, N = 1.5 * k, overlap = 0.5,
                         overlap_cor = 0.8))
  },
  n_pooled = function(k) {
    n_pooled(prec_pooled(var = 3.5, n = k, N = 5000, occasions = 3,
                         overlap = 0.5, overlap_cor = 0.6))
  },
  n_pooled_rotation = function(k) {
    n_pooled(prec_pooled(var = 3.5, n = k, N = 2 * k, occasions = 2,
                         overlap = rotation, overlap_cor = 0.8))
  },
  n_multi = function(k) {
    moe <- prec_prop(p = 0.3, n = k, N = 5000)$moe
    n_multi(data.frame(name = c("a", "b"), p = c(0.3, 0.4),
                       moe = c(moe, 0.5), N = 5000))
  },
  power_prop = function(k) {
    target <- power_prop(p1 = 0.3, p2 = 0.4, n = k, N = 5000,
                         power = NULL)$power
    power_prop(p1 = 0.3, p2 = 0.4, N = 5000, power = target, n = NULL)
  },
  power_mean = function(k) {
    target <- power_mean(effect = 1, var = 100, n = k, N = 5000,
                         power = NULL)$power
    power_mean(effect = 1, var = 100, N = 5000, power = target, n = NULL)
  }
)

test_that("the reviewer's tight frames field the size they came from", {
  p <- prec_pooled(var = 3.5, n = 4, N = 8, occasions = 2, overlap = rotation,
                   overlap_cor = 0.8)
  fit <- n_pooled(p)
  expect_identical(as.integer(fit), 4L)
  expect_true("n = 4 per occasion, 2 occasions (var = 3.5, deff = 1)" %in%
                capture.output(print(fit)))
  ch <- prec_change(var = 100, n = 20, N = 30, overlap = 0.5,
                    overlap_cor = 0.8)
  fit <- n_change(ch)
  expect_identical(as.integer(fit), 20L)
  expect_equal(prec_change(fit, n = as.integer(fit))$se, ch$se,
               tolerance = 1e-12)
})

test_that("every sizing route reports the whole size a precision came from", {
  sizes <- c(4:40, seq(100, 1000, 75))
  seen <- character(0)
  for (route in names(round_trips)) {
    for (k in sizes) {
      fit <- round_trips[[route]](k)
      label <- sprintf("%s at %d", route, k)
      expect_identical(as.integer(fit)[1L], as.integer(k), label = label)
    }
    seen <- c(seen, route)
  }
  expect_setequal(seen, names(round_trips))
})

test_that("the summary and print of a round trip read the same whole size", {
  for (k in c(7L, 29L, 361L)) {
    fit <- n_prop(prec_prop(p = 0.3, n = k, N = 5000))
    expect_identical(summary(fit)$sizes[["field"]], as.double(k))
    expect_true(any(startsWith(capture.output(print(fit)),
                               sprintf("n = %d (", k))))
  }
})

test_that("a size a fraction above a whole number still rounds up", {
  expect_identical(.whole_count(c(400.25, 400 + 1e-9, 4 + 4e-16, 400)),
                   c(401, 401, 4, 400))
  fit <- n_prop(p = 0.3, moe = 0.05)
  expect_gt(fit$n, floor(fit$n))
  expect_identical(as.integer(fit), as.integer(ceiling(fit$n)))
})

test_that("whole numbers of any size are kept and a tiny size takes a unit", {
  cases <- c(0, 1e-14, 4.000000000000001, 400 + 1e-9, 400.25, 1e12,
             1e12 + 0.5, 1e12 - 1, 2^53, 3.9999999999999996)
  expect_identical(.whole_count(cases),
                   c(0, 1, 4, 401, 401, 1e12, 1e12 + 1, 1e12 - 1, 2^53, 4))
  expect_identical(.whole_count(c(NA, Inf, 1e-300)), c(NA, Inf, 1))
  x <- n_mean(var = 1e-14, moe = 1)
  expect_lt(x$n, 1e-12)
  expect_identical(as.integer(x), 1L)
  expect_identical(summary(x)$sizes[["field"]], 1)
})

## Counts at the integer range

test_that("a large PSU frame does not overflow a modest multistage design", {
  big <- data.frame(stratum = "A", N = 2e10, N_1 = 3e9, N_2 = 6e9, mean = 10,
                    sd = 2, icc_1 = .1, icc_2 = .1, take_2 = 2, take_3 = 2)
  fit <- expect_no_warning(n_alloc(big, n = 400))
  expect_identical(fit$detail$n_int, 400L)
  expect_no_warning(expect_error(
    n_alloc(transform(big, N = 4e9, N_1 = 1e9, N_2 = 2e9),
            n = .Machine$integer.max - 1),
    "counts at most 2,147,483,647 units at the final stage"
  ))
})

test_that("an integer count past the range is refused, not made NA", {
  j <- n_alloc_joint(
    data.frame(stratum = "A", N = 1e10),
    measures = data.frame(stratum = "A", name = "y", p = 0.3),
    targets = data.frame(name = "y", cv = 0.1), min_n_stratum = 3e9
  )
  expect_error(as.integer(j), class = "svyplan_count_range")
  expect_identical(j$operational$n, 3e9)
  frac <- n_alloc_joint(
    data.frame(stratum = "A", N = 1e10),
    measures = data.frame(stratum = "A", name = "y", p = 0.3),
    targets = data.frame(name = "y", cv = 0.1), min_n_stratum = 3e9 + .25
  )
  expect_identical(frac$operational$n, 3000000001)
  expect_equal(frac$n, 3e9 + .25)
  at_max <- n_alloc(data.frame(stratum = c("A", "B"), N = 3e9, mean = 10,
                               sd = 2), n = .Machine$integer.max)
  expect_identical(as.integer(at_max), .Machine$integer.max)
})
