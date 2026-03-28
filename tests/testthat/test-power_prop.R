test_that("power_prop solve-n matches formula", {
  p1 <- 0.30
  p2 <- 0.35
  power <- 0.80
  alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  V <- p1 * (1 - p1) + p2 * (1 - p2)
  delta <- abs(p1 - p2)
  expected <- (z_a + z_b)^2 * V / delta^2

  res <- power_prop(p1 = p1, p2 = p2)
  expect_s3_class(res, "svyplan_power")
  expect_equal(res$n, expected, tolerance = 1e-6)
  expect_equal(res$solved, "n")
  expect_equal(res$type, "proportion")
  expect_equal(res$power, 0.80)
  expect_equal(res$effect, 0.05)
})

test_that("power_prop solve-power matches formula", {
  p1 <- 0.30
  p2 <- 0.35
  n <- 1500
  alpha <- 0.05
  z_a <- qnorm(1 - alpha / 2)
  V <- p1 * (1 - p1) + p2 * (1 - p2)
  delta <- abs(p1 - p2)
  se <- sqrt(V / n)
  expected <- pnorm(delta / se - z_a) + pnorm(-delta / se - z_a)

  res <- power_prop(p1 = p1, p2 = p2, n = n, power = NULL)
  expect_equal(res$power, expected, tolerance = 1e-6)
  expect_equal(res$solved, "power")
  expect_equal(res$n, 1500)
})

test_that("power_prop solve-mde returns valid p2", {
  res <- power_prop(p1 = 0.30, n = 1000)
  expect_equal(res$solved, "mde")
  expect_true(res$params$p2 > 0 && res$params$p2 < 1)
  expect_equal(res$effect, abs(res$params$p2 - 0.30), tolerance = 1e-8)
})

test_that("power_prop deff multiplier", {
  base <- power_prop(p1 = 0.30, p2 = 0.35)
  deff2 <- power_prop(p1 = 0.30, p2 = 0.35, deff = 2)
  expect_equal(deff2$n, base$n * 2, tolerance = 1e-6)
})

test_that("power_prop FPC reduces n", {
  base <- power_prop(p1 = 0.30, p2 = 0.35)
  fpc <- power_prop(p1 = 0.30, p2 = 0.35, N = 5000)
  expect_true(fpc$n < base$n)
})

test_that("power_prop panel overlap reduces n", {
  base <- power_prop(p1 = 0.30, p2 = 0.35)
  panel <- power_prop(p1 = 0.30, p2 = 0.35, overlap = 0.5, rho = 0.6)
  expect_true(panel$n < base$n)
})

test_that("power_prop one-sided needs smaller n", {
  two <- power_prop(p1 = 0.30, p2 = 0.35)
  one <- power_prop(p1 = 0.30, p2 = 0.35, alternative = "one.sided")
  expect_true(one$n < two$n)
})

test_that("power_prop round-trip n -> power", {
  res_n <- power_prop(p1 = 0.30, p2 = 0.35, power = 0.80)
  res_pw <- power_prop(p1 = 0.30, p2 = 0.35, n = res_n$n, power = NULL)
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("power_prop round-trip n -> mde", {
  res_n <- power_prop(p1 = 0.30, p2 = 0.35, power = 0.80)
  res_mde <- power_prop(p1 = 0.30, n = res_n$n, power = 0.80)
  expect_true(res_mde$effect <= abs(0.35 - 0.30) + 1e-3)
  expect_true(res_mde$effect > 0)
})

test_that("power_prop coercion methods", {
  res <- power_prop(p1 = 0.30, p2 = 0.35)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
  expect_equal(as.double(res), res$n)
})

test_that("power_prop validates inputs", {
  expect_error(power_prop(p1 = 0), "must be in \\(0, 1\\)")
  expect_error(power_prop(p1 = 0.3, p2 = 0.4, n = 100), "exactly one")
  expect_error(power_prop(p1 = 0.3, p2 = 0.3), "must differ")
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.4, n = 100, power = NULL,
               alternative = "bad"),
    "should be one of"
  )
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.4, n = 100, power = NULL, overlap = -1),
    "overlap"
  )
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.4, n = 100, power = NULL, rho = 2),
    "rho"
  )
})

test_that("power_prop format and print", {
  res <- power_prop(p1 = 0.30, p2 = 0.35)
  expect_match(format(res), "svyplan_power")
  expect_output(print(res), "Power analysis for proportions")
  expect_output(print(res), "solved for sample size")
})

test_that("power_prop MDE works for extreme high p1 (downward solve)", {
  res <- power_prop(p1 = 0.999, n = 100)
  expect_true(res$params$p2 < 0.999)
  expect_true(res$effect > 0)
})

test_that("power_prop MDE works for extreme low p1 (upward solve)", {
  res <- power_prop(p1 = 0.001, n = 100)
  expect_true(res$params$p2 > 0.001)
  expect_true(res$effect > 0)
})

test_that("power_prop MDE returns closest root", {
  res <- power_prop(p1 = 0.5, n = 500)
  expect_true(res$effect > 0)
  expect_true(res$params$p2 > 0 && res$params$p2 < 1)
})

test_that("power_prop MDE errors when no solution exists", {
  expect_error(
    power_prop(p1 = 0.5, n = 2, power = 0.999),
    "no detectable alternative"
  )
})

test_that("power_prop ratio returns length-2 n", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  expect_equal(length(res$n), 2L)
  expect_equal(res$n[1] / res$n[2], 2, tolerance = 1e-8)
})

test_that("power_prop unequal n power computation", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, n = c(500, 1000), power = NULL)
  res_eq <- power_prop(p1 = 0.30, p2 = 0.35, n = 750, power = NULL)
  expect_false(isTRUE(all.equal(res$power, res_eq$power)))
  expect_equal(length(res$n), 2L)
})

test_that("power_prop overlap + ratio reduces n", {
  res0 <- power_prop(p1 = 0.3, p2 = 0.35, ratio = 2)
  res_ov <- power_prop(p1 = 0.3, p2 = 0.35, ratio = 2,
                        overlap = 0.3, rho = 0.5)
  expect_true(res_ov$n[2] < res0$n[2])
})

test_that("power_prop overlap exceeds 1/ratio errors", {
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.35, ratio = 2, overlap = 0.6, rho = 0.5),
    "overlap.*must be.*1/ratio"
  )
})

test_that("power_prop overlap not supported with arcsine", {
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.35, overlap = 0.5, rho = 0.5,
               method = "arcsine"),
    "overlap is only supported with method = 'wald'"
  )
})

test_that("power_prop overlap not supported with logodds", {
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.35, overlap = 0.5, rho = 0.5,
               method = "logodds"),
    "overlap is only supported with method = 'wald'"
  )
})

test_that("Valliant Example 4.8 (Wald, one-sided)", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, power = 0.80,
                     alternative = "one.sided")
  expect_equal(ceiling(res$n), 1890L, tolerance = 2)
})

test_that("arcsine method matches Valliant Example 4.11", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, power = 0.80,
                     alternative = "one.sided", method = "arcsine")
  expect_equal(ceiling(res$n), 1889L, tolerance = 2)
})

test_that("logodds method matches Valliant Example 4.11", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, power = 0.80,
                     alternative = "one.sided", method = "logodds")
  expect_equal(ceiling(res$n), 1889L, tolerance = 2)
})

test_that("arcsine and logodds near-equivalence", {
  arc <- power_prop(p1 = 0.15, p2 = 0.18, power = 0.80,
                     alternative = "one.sided", method = "arcsine")
  lo <- power_prop(p1 = 0.15, p2 = 0.18, power = 0.80,
                    alternative = "one.sided", method = "logodds")
  expect_equal(ceiling(arc$n), ceiling(lo$n), tolerance = 5)
})

test_that("arcsine solve-power works", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, n = 1889,
                     power = NULL, alternative = "one.sided",
                     method = "arcsine")
  expect_equal(res$power, 0.80, tolerance = 0.01)
})

test_that("logodds solve-power works", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, n = 1889,
                     power = NULL, alternative = "one.sided",
                     method = "logodds")
  expect_equal(res$power, 0.80, tolerance = 0.01)
})

test_that("arcsine MDE returns valid p2", {
  res <- power_prop(p1 = 0.30, n = 1000, method = "arcsine")
  expect_equal(res$solved, "mde")
  expect_true(res$params$p2 > 0 && res$params$p2 < 1)
  expect_true(res$effect > 0)
})

test_that("logodds MDE returns valid p2", {
  res <- power_prop(p1 = 0.30, n = 1000, method = "logodds")
  expect_equal(res$solved, "mde")
  expect_true(res$params$p2 > 0 && res$params$p2 < 1)
  expect_true(res$effect > 0)
})

test_that("arcsine round-trip n -> power", {
  res_n <- power_prop(p1 = 0.20, p2 = 0.25, power = 0.80,
                       method = "arcsine")
  res_pw <- power_prop(p1 = 0.20, p2 = 0.25, n = res_n$n,
                        power = NULL, method = "arcsine")
  expect_equal(res_pw$power, 0.80, tolerance = 1e-3)
})

test_that("logodds round-trip n -> power", {
  res_n <- power_prop(p1 = 0.20, p2 = 0.25, power = 0.80,
                       method = "logodds")
  res_pw <- power_prop(p1 = 0.20, p2 = 0.25, n = res_n$n,
                        power = NULL, method = "logodds")
  expect_equal(res_pw$power, 0.80, tolerance = 1e-3)
})

test_that("arcsine ratio returns correct structure", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2, method = "arcsine")
  expect_equal(length(res$n), 2L)
  expect_equal(res$n[1] / res$n[2], 2, tolerance = 1e-8)
})

test_that("logodds ratio returns correct structure", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2, method = "logodds")
  expect_equal(length(res$n), 2L)
  expect_equal(res$n[1] / res$n[2], 2, tolerance = 1e-8)
})

test_that("swap-group invariance for solve-power (wald)", {
  pw1 <- power_prop(p1 = 0.30, p2 = 0.40, n = c(500, 800), power = NULL)
  pw2 <- power_prop(p1 = 0.40, p2 = 0.30, n = c(800, 500), power = NULL)
  expect_equal(pw1$power, pw2$power, tolerance = 1e-6)
})

test_that("power_prop ratio errors when n provided", {
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.35, n = 500, power = NULL, ratio = 2),
    "ratio.*cannot be used"
  )
})

test_that("power_prop ratio validation is clean when n provided", {
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.35, n = 500, power = NULL, ratio = NA_real_),
    "ratio.*positive finite scalar"
  )
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.35, n = 500, power = NULL, ratio = c(1, 2)),
    "ratio.*positive finite scalar"
  )
})

test_that("power_prop ratio solve-n respects finite scalar N bounds", {
  res <- power_prop(p1 = 0.3, p2 = 0.45, power = 0.80, ratio = 2, N = 100)
  expect_equal(length(res$n), 2L)
  expect_lte(res$n[1], 100)
  expect_lte(res$n[2], 100)
})

test_that("power_prop ratio solve-n respects length-2 N bounds across methods", {
  methods <- c("wald", "arcsine", "logodds")
  for (m in methods) {
    res <- power_prop(p1 = 0.3, p2 = 0.5, power = 0.80,
                      ratio = 2, N = c(80, 120), method = m)
    expect_equal(length(res$n), 2L)
    expect_lte(res$n[1], 80)
    expect_lte(res$n[2], 120)
  }
})

test_that("power_prop errors when target is unattainable under finite N", {
  expect_error(
    power_prop(p1 = 0.3, p2 = 0.31, power = 0.80, ratio = 2, N = 100),
    "unattainable"
  )
})

test_that("power_prop finite N ratio solve-n round-trips (wald)", {
  res_n <- power_prop(p1 = 0.30, p2 = 0.35, power = 0.80,
                      ratio = 2, N = 2000, method = "wald")
  res_pw <- power_prop(p1 = 0.30, p2 = 0.35, n = res_n$n,
                       power = NULL, N = 2000, method = "wald")
  expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
})

test_that("power_prop finite vector N ratio solve-n round-trips across methods", {
  methods <- c("wald", "arcsine", "logodds")
  for (m in methods) {
    res_n <- power_prop(p1 = 0.30, p2 = 0.35, power = 0.80,
                        ratio = 1.5, N = c(2500, 1800), method = m)
    res_pw <- power_prop(p1 = 0.30, p2 = 0.35, n = res_n$n,
                         power = NULL, N = c(2500, 1800), method = m)
    expect_equal(res_pw$power, 0.80, tolerance = 1e-4)
  }
})

test_that("print shows method for non-wald", {
  res <- power_prop(p1 = 0.15, p2 = 0.18, method = "arcsine")
  expect_output(print(res), "method = arcsine")
})

test_that("print shows ratio when != 1", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  expect_output(print(res), "ratio")
})

test_that("print shows vector n", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  expect_output(print(res), "n1 = ")
  expect_output(print(res), "n2 = ")
  expect_output(print(res), "total = ")
})

test_that("format handles vector n", {
  res <- power_prop(p1 = 0.30, p2 = 0.35, ratio = 2)
  fmt <- format(res)
  expect_match(fmt, ",")
})
