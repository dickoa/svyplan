test_that("power_mean solve-n matches formula", {
  icc <- 5; var <- 100; power <- 0.80; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  V <- 2 * var
  expected <- (z_a + z_b)^2 * V / icc^2

  res <- power_mean(effect = icc, var = var)
  expect_s3_class(res, "svyplan_power")
  expect_equal(res$n, expected, tolerance = 1e-6)
  expect_equal(res$solved, "n")
  expect_equal(res$type, "mean")
  expect_equal(res$effect, 5)
})

test_that("power_mean solve-power matches formula", {
  icc <- 5; var <- 100; n <- 200; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  V <- 2 * var
  se <- sqrt(V / n)
  expected <- pnorm(icc / se - z_a) + pnorm(-icc / se - z_a)

  res <- power_mean(effect = icc, var = var, n = n, power = NULL)
  expect_equal(res$power, expected, tolerance = 1e-6)
  expect_equal(res$solved, "power")
})

test_that("power_mean solve-mde matches analytical formula", {
  var <- 100; n <- 500; power <- 0.80; alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  V <- 2 * var
  expected <- (z_a + z_b) * sqrt(V / n)

  res <- power_mean(var = var, n = n)
  expect_equal(res$effect, expected, tolerance = 1e-6)
  expect_equal(res$solved, "mde")
})

test_that("power_mean deff multiplier", {
  base <- power_mean(effect = 5, var = 100)
  deff2 <- power_mean(effect = 5, var = 100, deff = 2)
  expect_equal(deff2$n, base$n * 2, tolerance = 1e-6)
})

test_that("power_mean FPC reduces n", {
  base <- power_mean(effect = 5, var = 100)
  fpc <- power_mean(effect = 5, var = 100, N = 2000)
  expect_true(fpc$n < base$n)
})

test_that("power_mean panel overlap reduces n", {
  base <- power_mean(effect = 5, var = 100)
  panel <- power_mean(effect = 5, var = 100, overlap = 0.5, overlap_cor = 0.6)
  expect_true(panel$n < base$n)
})

test_that("power_mean one-sided needs smaller n", {
  two <- power_mean(effect = 5, var = 100)
  one <- power_mean(effect = 5, var = 100, alternative = "one.sided")
  expect_true(one$n < two$n)
})

test_that("power_mean round-trip n -> power", {
  res_n <- power_mean(effect = 5, var = 100, power = 0.80)
  res_pw <- power_mean(effect = 5, var = 100, n = res_n$n, power = NULL)
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("power_mean round-trip n -> mde", {
  res_n <- power_mean(effect = 5, var = 100, power = 0.80)
  res_mde <- power_mean(var = 100, n = res_n$n, power = 0.80)
  expect_equal(res_mde$effect, 5, tolerance = 1e-3)
})

test_that("power_mean coercion methods", {
  res <- power_mean(effect = 5, var = 100)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
  expect_equal(as.double(res), res$n)
})

test_that("power_mean validates inputs", {
  expect_error(power_mean(var = 100), "exactly one")
  expect_error(power_mean(effect = 5, var = 100, n = 200), "exactly one")
  expect_error(
    power_mean(effect = 5, var = 100, alternative = "bad"),
    "should be one of"
  )
  expect_error(power_mean(effect = 5, var = 100, overlap = -1), "overlap")
})

test_that("power_mean format and print", {
  res <- power_mean(effect = 5, var = 100)
  expect_match(format(res), "svyplan_power")
  expect_output(print(res), "Power analysis for means")
  expect_output(print(res), "solved for sample size")
})

test_that("unequal var scalar reduction", {
  res1 <- power_mean(effect = 5, var = c(100, 100))
  res2 <- power_mean(effect = 5, var = 100)
  expect_equal(res1$n, res2$n, tolerance = 1e-8)
})

test_that("unequal var gives different power than equal var at same total", {
  pw1 <- power_mean(effect = 5, var = c(100, 200), n = c(80, 120), power = NULL)
  pw2 <- power_mean(effect = 5, var = 150, n = c(80, 120), power = NULL)
  expect_false(isTRUE(all.equal(pw1$power, pw2$power)))
})

test_that("unequal n power computation", {
  res <- power_mean(effect = 5, var = 100, n = c(100, 200), power = NULL)
  res_eq <- power_mean(effect = 5, var = 100, n = 150, power = NULL)
  expect_false(isTRUE(all.equal(res$power, res_eq$power)))
  expect_equal(length(res$n), 2L)
})

test_that("ratio returns length-2 n with correct proportion", {
  res <- power_mean(effect = 5, var = 100, ratio = 2)
  expect_equal(length(res$n), 2L)
  expect_equal(res$n[1] / res$n[2], 2, tolerance = 1e-8)
})

test_that("swap-group invariance for solve-power", {
  pw1 <- power_mean(effect = 5, var = c(100, 200), n = c(50, 100), power = NULL)
  pw2 <- power_mean(effect = 5, var = c(200, 100), n = c(100, 50), power = NULL)
  expect_equal(pw1$power, pw2$power, tolerance = 1e-8)
})

test_that("round-trip with unequal var", {
  res_n <- power_mean(effect = 5, var = c(80, 120), power = 0.80)
  res_pw <- power_mean(effect = 5, var = c(80, 120), n = res_n$n, power = NULL)
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("overlap + unequal n increases power", {
  pw0 <- power_mean(effect = 5, var = c(100, 200), n = c(100, 200), power = NULL)
  pw_ov <- power_mean(effect = 5, var = c(100, 200), n = c(100, 200),
                       power = NULL, overlap = 0.3, overlap_cor = 0.7)
  expect_true(pw_ov$power > pw0$power)
})

test_that("overlap + ratio reduces n", {
  res0 <- power_mean(effect = 5, var = 100, ratio = 2)
  res_ov <- power_mean(effect = 5, var = 100, ratio = 2,
                        overlap = 0.4, overlap_cor = 0.5)
  expect_true(res_ov$n[2] < res0$n[2])
})

test_that("overlap exceeds 1/ratio errors", {
  expect_error(
    power_mean(effect = 5, var = 100, ratio = 2, overlap = 0.6, overlap_cor = 0.5),
    "overlap.*must be.*1/ratio"
  )
})

test_that("overlap exceeds n[2]/n[1] errors", {
  expect_error(
    power_mean(effect = 5, var = 100, n = c(100, 50), power = NULL,
               overlap = 0.6, overlap_cor = 0.5),
    "overlap.*must be.*n\\[2\\]/n\\[1\\]"
  )
})

test_that("Valliant Example 4.6", {
  res <- power_mean(effect = 5, var = 200, power = 0.80,
                     alternative = "one.sided")
  expect_equal(ceiling(res$n), 99L, tolerance = 1)

  res_deff <- power_mean(effect = 5, var = 200, power = 0.80,
                          alternative = "one.sided", deff = 1.6)
  expect_equal(ceiling(res_deff$n), 159L, tolerance = 2)
})

test_that("ratio errors when n provided", {
  expect_error(
    power_mean(effect = 5, var = 100, n = 200, power = NULL, ratio = 2),
    "ratio.*cannot be used"
  )
})

test_that("ratio validation is clean when n provided", {
  expect_error(
    power_mean(effect = 5, var = 100, n = 200, power = NULL, ratio = NA_real_),
    "ratio.*positive finite scalar"
  )
  expect_error(
    power_mean(effect = 5, var = 100, n = 200, power = NULL, ratio = c(1, 2)),
    "ratio.*positive finite scalar"
  )
})

test_that("ratio solve-n respects finite scalar N bounds", {
  res <- power_mean(effect = 8, var = 100, power = 0.80, ratio = 2, N = 100)
  expect_equal(length(res$n), 2L)
  expect_lte(res$n[1], 100)
  expect_lte(res$n[2], 100)
})

test_that("ratio solve-n respects length-2 N bounds", {
  res <- power_mean(effect = 8, var = 100, power = 0.80,
                    ratio = 2, N = c(80, 120))
  expect_equal(length(res$n), 2L)
  expect_lte(res$n[1], 80)
  expect_lte(res$n[2], 120)
})

test_that("power_mean errors when target is unattainable under finite N", {
  expect_error(
    power_mean(effect = 1, var = 100, power = 0.80, ratio = 2, N = 100),
    "unattainable"
  )
})

test_that("finite N ratio solve-n round-trips to target power", {
  res_n <- power_mean(effect = 5, var = 100, power = 0.80,
                      ratio = 2, N = 2000)
  res_pw <- power_mean(effect = 5, var = 100, n = res_n$n,
                       power = NULL, N = 2000)
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("finite vector N solve-n round-trips to target power", {
  res_n <- power_mean(effect = 3, var = c(120, 80), power = 0.80,
                      ratio = 1.5, N = c(3000, 1500))
  res_pw <- power_mean(effect = 3, var = c(120, 80), n = res_n$n,
                       power = NULL, N = c(3000, 1500))
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("power_mean handles full and partial censuses componentwise", {
  res <- power_mean(effect = 5, var = 100, n = c(400, 400),
                    power = NULL, N = 400)
  expect_equal(res$power, 1)

  partial <- power_mean(effect = 0.1, var = 1, n = c(100, 10),
                        power = NULL, N = c(100, 100))
  se <- sqrt((1 - 10 / 100) / 10)
  expected <- pnorm(0.1 / se - qnorm(0.975)) +
    pnorm(-0.1 / se - qnorm(0.975))
  expect_equal(partial$power, expected, tolerance = 1e-12)
  expect_lt(partial$power, 1)

  expect_error(
    power_mean(effect = 5, var = 100, n = c(500, 500), power = NULL, N = 400),
    "cannot draw"
  )
})

test_that("the overlap covariance carries one 1/N, not sample 2's correction", {
  # two SRSWOR samples from one population of size N sharing k = overlap * n1
  # units have Cov = rho * S1 * S2 * (k/(n1 n2) - 1/N)
  za <- qnorm(0.975)
  for (rho in c(0.3, 0.7, 1)) {
    for (ov in c(0.25, 0.6)) {
      v <- c(4, 9); nn <- c(50, 40); NN <- 500
      exact <- v[1] * (1 / nn[1] - 1 / NN) + v[2] * (1 / nn[2] - 1 / NN) -
        2 * rho * sqrt(v[1] * v[2]) * (ov / nn[2] - 1 / NN)
      got <- power_mean(var = v, effect = 1, n = nn, power = NULL, N = NN,
                        overlap = ov, overlap_cor = rho)
      se <- 1 / (qnorm(got$power - pnorm(-1 / sqrt(exact) - za)) + za)
      expect_equal(se^2, exact, tolerance = 1e-7)
    }
  }
})

test_that("a perfectly correlated full-overlap difference is free of N", {
  # the population terms cancel: V = 2 S^2 (1 - overlap) / n at rho = 1
  n <- 40; S2 <- 1; effect <- 0.2; overlap <- 0.5
  za <- qnorm(0.975)
  V <- 2 * S2 * (1 - overlap) / n
  reference <- pnorm(abs(effect) / sqrt(V) - za) +
    pnorm(-abs(effect) / sqrt(V) - za)
  for (N in c(100, 1000, Inf)) {
    got <- power_mean(var = S2, effect = effect, n = n, power = NULL, N = N,
                      overlap = overlap, overlap_cor = 1)
    expect_equal(got$power, reference, tolerance = 1e-9)
  }
})

test_that("overlapping samples must come from one population", {
  expect_error(
    power_mean(var = 1, effect = 1, n = 40, power = NULL, N = c(100, 200),
               overlap = 0.5, overlap_cor = 1),
    "one population"
  )
  # without overlap the two groups may be different populations
  expect_silent(
    power_mean(var = 1, effect = 1, n = 40, power = NULL, N = c(100, 200))
  )
})
