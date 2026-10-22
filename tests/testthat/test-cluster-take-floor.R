## Whole and continuous takes keep one expected respondent

.floor_cases <- function(r) {
  i2 <- data.frame(name = "p", p = 0.5, cv = 0.05, icc_psu = 0.9,
                   resp_rate = r)
  i3 <- data.frame(name = "p", p = 0.5, cv = 0.05, icc_psu = 0.05,
                   icc_ssu = 0.8, resp_rate = r)
  i3s <- data.frame(name = "p", p = 0.5, cv = 0.05, icc_psu = 0.6,
                    icc_ssu = 0.01, resp_rate_ssu = r)
  list(
    two_cv = list(
      scalar = function() n_cluster(stage_cost = c(1, 100), icc = 0.9,
                                    cv = 0.05, resp_rate = r),
      table = function() n_cluster(indicators = i2, stage_cost = c(1, 100)),
      icc = 0.9, rates = c(resp_rate = r)
    ),
    two_budget = list(
      scalar = function() n_cluster(stage_cost = c(1, 100), icc = 0.9,
                                    budget = 1e5, resp_rate = r),
      table = function() n_cluster(indicators = i2, stage_cost = c(1, 100),
                                   budget = 1e5),
      icc = 0.9, rates = c(resp_rate = r)
    ),
    three_cv = list(
      scalar = function() n_cluster(stage_cost = c(500, 50, 500),
                                    icc = c(0.05, 0.8), cv = 0.05,
                                    resp_rate = r),
      table = function() n_cluster(indicators = i3,
                                   stage_cost = c(500, 50, 500)),
      icc = c(0.05, 0.8), rates = c(resp_rate = r)
    ),
    three_budget = list(
      scalar = function() n_cluster(stage_cost = c(500, 50, 500),
                                    icc = c(0.05, 0.8), budget = 3e5,
                                    resp_rate = r),
      table = function() n_cluster(indicators = i3,
                                   stage_cost = c(500, 50, 500),
                                   budget = 3e5),
      icc = c(0.05, 0.8), rates = c(resp_rate = r)
    ),
    three_ssu = list(
      scalar = function() n_cluster(stage_cost = c(5000, 500, 5),
                                    icc = c(0.6, 0.01), cv = 0.05,
                                    resp_rate_ssu = r),
      table = function() n_cluster(indicators = i3s,
                                   stage_cost = c(5000, 500, 5)),
      icc = c(0.6, 0.01), rates = c(resp_rate_ssu = r)
    )
  )
}

.assess <- function(n, case) {
  args <- c(list(n = n, icc = case$icc), as.list(case$rates))
  do.call(prec_cluster, args)
}

test_that("every whole design is one prec_cluster() accepts, at its cv", {
  for (r in c(0.5, 0.2, 0.1)) {
    for (nm in names(cases <- .floor_cases(r))) {
      case <- cases[[nm]]
      fit <- suppressWarnings(case$scalar())
      info <- paste(nm, r)
      assessed <- expect_no_error(.assess(fit$operational$n, case))
      expect_equal(assessed$cv, fit$operational$cv, tolerance = 1e-10,
                   info = info)
      expect_no_error(.assess(fit$n, case))
    }
  }
})

test_that("one indicator in a table is the scalar design", {
  for (r in c(0.5, 0.2, 0.1)) {
    for (nm in names(cases <- .floor_cases(r))) {
      case <- cases[[nm]]
      info <- paste(nm, r)
      scalar <- suppressWarnings(case$scalar())
      table <- expect_no_error(case$table())
      expect_equal(unname(table$n), unname(scalar$n), tolerance = 1e-5,
                   info = info)
      expect_identical(unname(table$operational$n),
                       unname(scalar$operational$n), info = info)
    }
  }
})

test_that("the review's table examples solve at a take of two", {
  ind <- data.frame(name = "p", p = 0.5, cv = 0.05, icc_psu = 0.9,
                    resp_rate = 0.5)
  free <- n_cluster(indicators = ind, stage_cost = c(1, 100))
  fixed <- n_cluster(indicators = ind, stage_cost = c(1, 100), n_per_psu = 2)
  expect_equal(free$n, fixed$n)
  expect_equal(unname(free$n), c(400, 2))
  budget <- n_cluster(indicators = ind, stage_cost = c(1, 100), budget = 1e5)
  expect_equal(unname(budget$n), c(1e5 / 201, 2))
  expect_identical(unname(budget$operational$n), c(497L, 2L))
})

test_that("a fractional fixed take rounds up to a supported whole take", {
  scalar <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                      resp_rate = 0.42, n_per_psu = 2.4)
  expect_identical(unname(scalar$operational$n[["n_per_psu"]]), 3L)
  ind <- data.frame(name = "p", p = 0.5, cv = 0.05, icc_psu = 0.05,
                    resp_rate = 0.42)
  table <- n_cluster(indicators = ind, stage_cost = c(500, 50),
                     n_per_psu = 2.4)
  expect_identical(table$operational$n, scalar$operational$n)
  expect_no_error(prec_cluster(indicators = transform(
    ind, n = table$operational$n[["n_psu"]],
    n_per_psu = table$operational$n[["n_per_psu"]]
  )))
})

test_that("a whole fixed take below one respondent stays an error", {
  expect_error(
    n_cluster(stage_cost = c(1, 100), icc = 0.9, cv = 0.05, resp_rate = 0.5,
              n_per_psu = 1),
    "unsupported below one"
  )
  ind <- data.frame(name = "p", p = 0.5, cv = 0.05, icc_psu = 0.9,
                    resp_rate = 0.5)
  expect_error(
    n_cluster(indicators = ind, stage_cost = c(1, 100), n_per_psu = 1),
    "unsupported below one"
  )
})

test_that("eligibility enters the floor through the ultimate rate", {
  fit <- suppressWarnings(n_cluster(stage_cost = c(1, 100), icc = 0.9,
                                    cv = 0.05, eligible_rate = 0.25))
  expect_identical(unname(fit$operational$n[["n_per_psu"]]), 4L)
  above <- suppressWarnings(n_cluster(stage_cost = c(1, 100), icc = 0.9,
                                      cv = 0.05, eligible_rate = 2))
  expect_identical(unname(above$operational$n[["n_per_psu"]]), 1L)
})

test_that("the take floors and the bounded optimizer", {
  expect_identical(svyplan:::.take_floor(0.42), 3L)
  expect_identical(svyplan:::.take_floor(1 / 3), 3L)
  expect_identical(svyplan:::.take_floor(c(0.9, 0.2)), 5L)
  expect_identical(svyplan:::.take_floor(2), 1L)
  rising <- svyplan:::.optimize_from(function(x) x^2, 2, 10)
  expect_identical(rising$minimum, 2)
  interior <- svyplan:::.optimize_from(function(x) (x - 5)^2, 2, 10)
  expect_equal(interior$minimum, 5, tolerance = 1e-4)
})
