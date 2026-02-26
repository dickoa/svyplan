test_that("prec_multi simple mode computes per-indicator precision", {
  targets <- data.frame(
    name = c("stunting", "vaccination", "anemia"),
    p    = c(0.30, 0.70, 0.10),
    n    = c(400, 400, 400)
  )
  result <- prec_multi(targets)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "multi")
  expect_equal(nrow(result$detail), 3L)

  for (i in 1:3) {
    p <- targets$p[i]
    se_exp <- sqrt(p * (1 - p) / 400)
    expect_equal(result$detail$.se[i], se_exp, tolerance = 1e-6)
  }
})

test_that("prec_multi simple mode with deff and resp_rate", {
  targets <- data.frame(
    name      = c("stunting", "vaccination"),
    p         = c(0.30, 0.70),
    n         = c(500, 500),
    deff      = c(1.5, 2.0),
    resp_rate = c(0.8, 0.9)
  )
  result <- prec_multi(targets)
  expect_equal(nrow(result$detail), 2L)

  n_eff <- 500 * 0.8 / 1.5
  se_exp <- sqrt(0.3 * 0.7 / n_eff)
  expect_equal(result$detail$.se[1], se_exp, tolerance = 1e-6)
})

test_that("prec_multi multistage mode computes per-indicator CV", {
  targets <- data.frame(
    name   = c("stunting", "anemia"),
    p      = c(0.30, 0.10),
    n      = c(50, 50),
    n2     = c(12, 12),
    delta1 = c(0.02, 0.05)
  )
  result <- prec_multi(targets, cost = c(500, 50))
  expect_equal(result$type, "multi")
  expect_equal(nrow(result$detail), 2L)
  expect_true(all(!is.na(result$detail$.cv)))
})

test_that("prec_multi requires n column", {
  targets <- data.frame(p = 0.3, moe = 0.05)
  expect_error(prec_multi(targets), "must contain an 'n' column")
})

test_that("prec_multi simple rejects invalid p", {
  expect_error(
    prec_multi(data.frame(p = 0, n = 100)),
    "p.*\\(0, 1\\)"
  )
  expect_error(
    prec_multi(data.frame(p = 1, n = 100)),
    "p.*\\(0, 1\\)"
  )
  expect_error(
    prec_multi(data.frame(p = -0.1, n = 100)),
    "p.*\\(0, 1\\)"
  )
})

test_that("prec_multi simple rejects invalid var", {
  expect_error(
    prec_multi(data.frame(var = -10, n = 100)),
    "var.*positive"
  )
  expect_error(
    prec_multi(data.frame(var = 0, n = 100)),
    "var.*positive"
  )
})

test_that("prec_multi simple rejects non-positive mu", {
  expect_error(
    prec_multi(data.frame(var = 10, mu = 0, n = 100)),
    "mu.*positive"
  )
  expect_error(
    prec_multi(data.frame(var = 10, mu = -5, n = 100)),
    "mu.*positive"
  )
})

test_that("prec_multi cluster requires n2 column", {
  targets <- data.frame(
    p = 0.3, n = 50, delta1 = 0.05
  )
  expect_error(
    prec_multi(targets, cost = c(500, 50)),
    "n2.*required"
  )
})

test_that("prec_multi 3-stage requires n3 column", {
  targets <- data.frame(
    p = 0.3, n = 50, n2 = 10, delta1 = 0.05
  )
  expect_error(
    prec_multi(targets, cost = c(500, 100, 50)),
    "n3.*required"
  )
})

test_that("prec_multi cluster rejects NA in n2", {
  targets <- data.frame(
    p = c(0.3, 0.1), n = c(50, 50), n2 = c(10, NA), delta1 = c(0.05, 0.02)
  )
  expect_error(
    prec_multi(targets, cost = c(500, 50)),
    "n2.*required"
  )
})
