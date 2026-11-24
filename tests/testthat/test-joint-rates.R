## Response and eligibility rates of joint allocations

rate_frame <- function() data.frame(stratum = c("A", "B"), N = c(4000, 3000))
rate_measures <- function() {
  data.frame(stratum = c("A", "B"), name = "y", p = c(0.5, 0.4))
}
rate_targets <- data.frame(name = "y", cv = 0.1)

## A scalar the columns override

test_that("a joint scalar rate that no row reads is refused", {
  fit <- function(frame, measures, ...) {
    n_alloc_joint(frame, measures, targets = rate_targets, ...)
  }
  fr <- rate_frame()
  m <- rate_measures()
  expect_error(fit(transform(fr, resp_rate = 1), m, resp_rate = 0.65),
               "frame column 'resp_rate' has a value on every row")
  expect_error(fit(fr, transform(m, resp_rate = c(0.8, 0.9)), resp_rate = 0.65),
               "uses takes 'resp_rate'")
  expect_error(
    fit(transform(fr, resp_rate = c(NA, 0.7)),
        transform(m, resp_rate = c(0.8, NA)), resp_rate = 0.65),
    "uses takes 'resp_rate'"
  )
  expect_error(
    fit(transform(fr, eligibility_rate = 0.9), m, eligibility_rate = 0.8),
    "frame column 'eligibility_rate'"
  )
  expect_error(fit(fr, transform(m, deff = 1.2), deff = 1.5),
               "uses takes 'deff'")
  expect_error(
    prec_alloc_joint(transform(fr, resp_rate = 1), m, targets = rate_targets,
                     n = c(100, 100), resp_rate = 0.65),
    "frame column 'resp_rate'"
  )
})

test_that("an indicator that no target uses cannot keep an ignored scalar", {
  fr <- rate_frame()
  m <- transform(rate_measures(), resp_rate = c(0.8, 0.9))
  extra <- rbind(m, transform(m[1L, ], name = "unused", resp_rate = NA_real_))
  expect_error(
    n_alloc_joint(fr, extra, targets = rate_targets, resp_rate = 0.65),
    "uses takes 'resp_rate'"
  )
  expect_error(
    prec_alloc_joint(fr, extra, targets = rate_targets, n = c(100, 100),
                     resp_rate = 0.65),
    "uses takes 'resp_rate'"
  )
  two <- rbind(m, transform(m, name = "z", p = c(0.3, 0.2),
                            resp_rate = c(NA, 0.7)))
  objective <- data.frame(name = "z")
  used <- n_alloc_joint(fr, two, targets = rate_targets, objective = objective,
                        budget = 300, resp_rate = 0.65)
  expect_identical(used$rates$resp_rate[used$rates$name == "z"], c(0.65, 0.7))
})

test_that("a fit refuses an override that its columns would ignore", {
  fr <- rate_frame()
  m <- transform(rate_measures(), resp_rate = c(0.8, 0.9))
  stage <- transform(fr, N_1 = c(40, 30), take_2 = c(10, 8), cost_1 = 100,
                     cost_2 = 10)
  ms <- transform(m, icc_1 = 0.05)
  psu <- data.frame(stratum = rep(c("A", "B"), each = 8),
                    N = c(600, rep(200, 7), 450, 220, 165, rep(153, 5)))
  register <- data.frame(stratum = c("A", "B"), N = c(2000, 1600),
                         take_2 = c(10, 8))
  fits <- list(
    element = n_alloc_joint(fr, m, targets = rate_targets),
    fixed = n_alloc_joint(stage, ms, targets = rate_targets),
    free = n_alloc_joint(transform(stage, take_2 = NULL), ms,
                         targets = rate_targets),
    register = n_alloc_joint(register, ms, targets = rate_targets, psu = psu)
  )
  for (route in names(fits)) {
    expect_error(prec_alloc_joint(fits[[route]], resp_rate = 0.65),
                 "uses takes 'resp_rate'", label = route)
    expect_error(prec_alloc_joint(fits[[route]], n = c(100, 100),
                                  resp_rate = 0.65),
                 "uses takes 'resp_rate'", label = route)
  }
  by_arg <- n_alloc_joint(fr, rate_measures(), targets = rate_targets,
                          resp_rate = 0.65)
  expect_error(prec_alloc_joint(by_arg, measures = m), "uses takes 'resp_rate'")

  partial <- n_alloc_joint(fr, transform(m, resp_rate = c(0.8, NA)),
                           targets = rate_targets)
  assessed <- prec_alloc_joint(partial, resp_rate = 0.7)
  expect_identical(assessed$params$resp_rate, 0.7)
  expect_no_error(n_alloc_joint(assessed))
})

test_that("a joint scalar rate that one row reads is kept", {
  fr <- transform(rate_frame(), resp_rate = c(NA, 0.7))
  m <- rate_measures()
  used <- n_alloc_joint(fr, m, targets = rate_targets, resp_rate = 0.65)
  spelled <- n_alloc_joint(transform(fr, resp_rate = c(0.65, 0.7)), m,
                           targets = rate_targets)
  expect_equal(used$detail$n, spelled$detail$n)
  expect_no_error(
    n_alloc_joint(transform(rate_frame(), resp_rate = 0.7), m,
                  targets = rate_targets, resp_rate = 1)
  )
})

## The rates a fit used

## Oracle: refitting with the table's rates written as plain measure or
## frame columns, every argument left at 1, gives the same allocation, so
## the table holds the rates the variance read.

joint_rate_fixture <- function() {
  list(
    measures = data.frame(stratum = rep(c("A", "B"), 2),
                          name = rep(c("y", "z"), each = 2),
                          p = c(0.5, 0.4, 0.3, 0.2),
                          resp_rate = c(NA, NA, 0.3, 0.3), icc_1 = 0.05),
    frame = data.frame(stratum = c("A", "B"), N = c(4000, 3000),
                       resp_rate = c(0.5, NA), eligibility_rate = c(0.9, NA))
  )
}

refit_from_rates <- function(fit, frame, measures, ...) {
  r <- fit$rates
  key <- paste(measures$stratum, measures$name)
  at <- match(key, paste(r$stratum, r$name))
  rate_cols <- setdiff(names(r), c("stratum", "name"))
  frame <- frame[setdiff(names(frame), rate_cols)]
  for (nm in setdiff(rate_cols, "eligibility_rate")) {
    measures[[nm]] <- r[[nm]][at]
  }
  frame$eligibility_rate <- r$eligibility_rate[match(frame$stratum, r$stratum)]
  used <- !is.na(at)
  n_alloc_joint(frame, measures[used, , drop = FALSE],
                targets = rate_targets, ...)
}

test_that("a joint fit reports the rates it used on every route", {
  fx <- joint_rate_fixture()
  stage <- transform(fx$frame, N_1 = c(40, 30), take_2 = c(10, 8),
                     cost_1 = 100, cost_2 = 10, resp_rate_1 = c(0.9, NA))
  psu <- data.frame(stratum = rep(c("A", "B"), each = 8),
                    N = c(600, rep(200, 7), 450, 220, 165, rep(153, 5)))
  register <- data.frame(stratum = c("A", "B"), N = c(2000, 1600),
                         take_2 = c(10, 8), resp_rate = c(0.5, NA),
                         eligibility_rate = c(0.9, NA))
  fits <- list(
    element = list(fx$frame, fx$measures[-5], list()),
    fixed = list(stage, fx$measures, list()),
    free = list(transform(stage, take_2 = NULL), fx$measures, list()),
    register = list(register, fx$measures, list(psu = psu))
  )
  checked <- character()
  for (route in names(fits)) {
    f <- fits[[route]]
    fit <- do.call(n_alloc_joint, c(list(f[[1]], f[[2]], targets = rate_targets,
                                         resp_rate = 0.65), f[[3]]))
    r <- fit$rates
    expect_identical(r$name, c("y", "y"), label = route)
    expect_identical(r$resp_rate, c(0.5, 0.65), label = route)
    expect_identical(r$eligibility_rate, c(0.9, 1), label = route)
    if (route %in% c("fixed", "free")) {
      expect_identical(r$resp_rate_1, c(0.9, 1), label = route)
    } else {
      expect_null(r$resp_rate_1, label = route)
    }
    again <- do.call(refit_from_rates, c(list(fit, f[[1]], f[[2]]), f[[3]]))
    expect_equal(again$detail$n, fit$detail$n, label = route)
    expect_equal(again$detail$n_int, fit$detail$n_int, label = route)
    checked <- c(checked, route)
  }
  expect_setequal(checked, names(fits))
})

test_that("a classic allocation reports the rates it used per stratum", {
  frame <- data.frame(stratum = c("A", "B"), N = c(4000, 3000),
                      sd = c(10, 12), resp_rate = c(0.5, NA),
                      eligibility_rate = c(0.9, NA))
  scalar <- n_alloc(frame, n = 400, resp_rate = 0.7, eligibility_rate = 0.8)
  expect_identical(scalar$rates$resp_rate, c(0.5, 0.7))
  expect_identical(scalar$rates$eligibility_rate, c(0.9, 0.8))
  vector <- n_alloc(frame, n = 400, resp_rate = c(0.6, 0.65))
  expect_identical(vector$rates$resp_rate, c(0.6, 0.65))

  cluster <- n_alloc(
    transform(frame, icc_1 = 0.05, take_2 = c(10, 8), cost_1 = 500,
              cost_2 = 10, N_1 = c(80, 60), resp_rate_1 = c(0.95, NA)),
    n = 400, resp_rate = 0.7
  )
  expect_identical(cluster$rates$resp_rate_1, c(0.95, 1))
  expect_identical(cluster$rates$resp_rate, c(0.5, 0.7))

  for (fit in list(scalar, vector)) {
    spelled <- frame
    spelled$resp_rate <- fit$rates$resp_rate
    spelled$eligibility_rate <- fit$rates$eligibility_rate
    expect_equal(n_alloc(spelled, n = 400)$detail$n, fit$detail$n)
  }
})
