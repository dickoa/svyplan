test_that("two-stage fixed-take coefficients match the cluster convention", {
  z <- .bethel_multistage_fixture(2L)
  p <- .build_bethel_problem(
    z$frame, z$measures, z$targets, deff = 1.4, resp_rate = 0.8
  )
  k <- match("vaccination@.overall:cv", p$constraint_ids)
  rows <- seq_len(nrow(z$frame))
  mi <- rows
  variance <- z$measures$p[mi] * (1 - z$measures$p[mi])
  # The clustering penalty is paid on the take that responds, so the bracket
  # reads 'n_per_psu * resp_rate'. Reading the gross take here would inflate
  # the between-PSU component by 1/resp_rate, which is whole-PSU loss and a
  # mechanism 'resp_rate_psu' carries separately.
  inflation <- z$measures$var_ratio_psu[mi] *
    (1 + z$measures$icc_psu[mi] * (z$frame$n_per_psu * 0.8 - 1))

  expect_identical(p$stages, 2L)
  expect_equal(
    p$A[, k],
    z$frame$N^2 * variance * inflation * 1.4 /
      (0.8 * z$frame$n_per_psu)
  )
  expect_equal(p$B[k], -sum(z$frame$N * variance * inflation * 1.4))
  expect_equal(
    p$cost,
    z$frame$cost_psu + z$frame$cost_ssu * z$frame$n_per_psu
  )

  n_psu <- c(15, 20, 18, 12)
  got <- .precision_from_allocation(p, n_psu)
  legacy <- .alloc_metrics(
    z$frame$N,
    sqrt(variance * inflation),
    z$measures$p[mi],
    n_psu * z$frame$n_per_psu,
    alpha = 0.05,
    deff = 1.4,
    resp_rate = 0.8,
    cost_h = rep(1, 4)
  )
  expect_equal(got$cv[k], legacy$cv)
  expect_equal(got$se[k], legacy$se)
})

test_that("three-stage coefficients reproduce the established stage formula", {
  z <- .bethel_multistage_fixture(3L)
  z$measures$resp_rate <- 0.7
  z$measures$resp_rate_ssu <- 0.9
  z$measures$resp_rate_psu <- 0.8
  p <- .build_bethel_problem(z$frame, z$measures, z$targets)
  k <- match("vaccination@.overall:cv", p$constraint_ids)
  rows <- seq_len(nrow(z$frame))
  variance <- z$measures$p[rows] * (1 - z$measures$p[rows])
  m <- z$frame$n_per_psu
  q <- z$frame$n_per_ssu
  # Each stage's realized take, as in n_cluster(): the SSU rate shrinks the
  # PSU's realized SSU count and the ultimate rate shrinks the final take.
  mr <- m * 0.9
  qr <- q * 0.7
  bracket <- z$measures$var_ratio_psu[rows] *
    z$measures$icc_psu[rows] * mr * qr +
    z$measures$var_ratio_ssu[rows] *
      (1 + z$measures$icc_ssu[rows] * (qr - 1))
  resp <- 0.8 * 0.9 * 0.7

  expect_identical(p$stages, 3L)
  expect_equal(p$A[, k], z$frame$N^2 * variance * bracket / (resp * m * q))
  expect_equal(p$B[k], -sum(z$frame$N * variance * bracket))
  expect_equal(
    p$cost,
    z$frame$cost_psu + z$frame$cost_ssu * m +
      z$frame$cost_tsu * m * q
  )

  # Removing the FPC term leaves exactly the n_cluster() CV formula
  # within each atomic stratum.
  n_psu <- c(8, 12, 10, 7)
  no_fpc_var <- p$A[, k] / n_psu
  expected <- z$frame$N^2 * variance /
    (resp * n_psu * m * q) * bracket
  expect_equal(no_fpc_var, expected)

  got <- .precision_from_allocation(p, n_psu)
  legacy <- .alloc_metrics(
    z$frame$N,
    sqrt(variance * bracket),
    z$measures$p[rows],
    n_psu * m * q,
    alpha = 0.05,
    deff = 1,
    resp_rate = resp,
    cost_h = rep(1, 4)
  )
  expect_equal(got$cv[k], legacy$cv)
  expect_equal(got$se[k], legacy$se)
})

test_that("stage parameters use measure, frame, then default precedence", {
  z <- .bethel_multistage_fixture(3L)
  z$frame$icc_psu <- c(0.2, 0.2, 0.2, 0.2)
  z$frame$icc_ssu <- c(0.15, 0.15, 0.15, 0.15)
  z$measures$icc_psu[c(1, 5)] <- NA_real_
  z$measures$icc_ssu[c(1, 5)] <- NA_real_
  z$measures$var_ratio_psu <- NULL
  z$measures$var_ratio_ssu <- NULL
  p <- .build_bethel_problem(z$frame, z$measures, z$targets)
  k <- match("vaccination@.overall:cv", p$constraint_ids)
  variance <- 0.5 * 0.5
  m <- z$frame$n_per_psu[1]
  q <- z$frame$n_per_ssu[1]
  bracket <- 0.2 * m * q + (1 - 0.2) * (1 + 0.15 * (q - 1))
  expect_equal(p$A[1, k], z$frame$N[1]^2 * variance * bracket / (m * q))
})

test_that("multistage population and ultimate-unit bounds are converted", {
  z <- .bethel_multistage_fixture(3L)
  z$frame$max_weight <- c(20, NA, NA, NA)
  p <- .build_bethel_problem(
    z$frame, z$measures, z$targets, min_n_stratum = 40
  )
  take <- z$frame$n_per_psu * z$frame$n_per_ssu
  expect_equal(p$lower, pmax(1, pmax(40, z$frame$N / c(20, Inf, Inf, Inf)) / take))
  expect_equal(
    p$upper,
    pmin(z$frame$N_psu, z$frame$N_ssu / z$frame$n_per_psu,
         z$frame$N / take)
  )
})

test_that("fixed-take stage validation rejects ambiguous and free-stage inputs", {
  z <- .bethel_multistage_fixture(2L)
  z$frame$unit_cost <- 1
  expect_error(
    .build_bethel_problem(z$frame, z$measures, z$targets),
    "uses stage costs"
  )
  z <- .bethel_multistage_fixture(2L)
  z$frame$n_psu <- 10
  expect_error(
    .build_bethel_problem(z$frame, z$measures, z$targets),
    "free decision"
  )
  z <- .bethel_multistage_fixture(2L)
  z$frame$take_all <- TRUE
  expect_error(
    .build_bethel_problem(z$frame, z$measures, z$targets),
    "does not imply an ultimate-unit census"
  )
  z <- .bethel_multistage_fixture(3L)
  z$frame$n_per_ssu[1] <- 2.5
  expect_error(
    .build_bethel_problem(z$frame, z$measures, z$targets),
    "positive whole numbers"
  )
})
