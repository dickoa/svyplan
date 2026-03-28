test_that("n_prop resp_rate inflates sample size", {
  base <- n_prop(p = 0.3, moe = 0.05)
  rr <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.8)
  expect_equal(rr$n, base$n / 0.8, tolerance = 1e-6)
  expect_equal(rr$params$resp_rate, 0.8)
})

test_that("n_prop resp_rate = 1 gives same result", {
  base <- n_prop(p = 0.3, moe = 0.05)
  rr <- n_prop(p = 0.3, moe = 0.05, resp_rate = 1)
  expect_equal(rr$n, base$n, tolerance = 1e-6)
})

test_that("n_mean resp_rate inflates sample size", {
  base <- n_mean(var = 100, moe = 2)
  rr <- n_mean(var = 100, moe = 2, resp_rate = 0.8)
  expect_equal(rr$n, base$n / 0.8, tolerance = 1e-6)
})

test_that("n_cluster resp_rate inflates stage-1 in CV mode", {
  base <- n_cluster(stage_cost = c(500, 50), delta = 0.05, cv = 0.05)
  rr <- n_cluster(stage_cost = c(500, 50), delta = 0.05, cv = 0.05,
                   resp_rate = 0.8)
  expect_equal(rr$n[1], base$n[1] / 0.8, tolerance = 1e-6)
  expect_equal(unname(rr$n[2]), unname(base$n[2]), tolerance = 1e-6)
})

test_that("power_prop resp_rate works for solve-n", {
  base <- power_prop(p1 = 0.3, p2 = 0.35, n = NULL, power = 0.8)
  rr <- power_prop(p1 = 0.3, p2 = 0.35, n = NULL, power = 0.8,
                    resp_rate = 0.8)
  expect_equal(rr$n, base$n / 0.8, tolerance = 1e-6)
})

test_that("power_mean resp_rate works for solve-n", {
  base <- power_mean(effect = 5, var = 100, n = NULL, power = 0.8)
  rr <- power_mean(effect = 5, var = 100, n = NULL, power = 0.8,
                    resp_rate = 0.8)
  expect_equal(rr$n, base$n / 0.8, tolerance = 1e-6)
})

test_that("resp_rate validation works", {
  expect_error(n_prop(p = 0.3, moe = 0.05, resp_rate = 0), "resp_rate")
  expect_error(n_prop(p = 0.3, moe = 0.05, resp_rate = 1.5), "resp_rate")
  expect_error(n_prop(p = 0.3, moe = 0.05, resp_rate = -0.1), "resp_rate")
  expect_error(n_prop(p = 0.3, moe = 0.05, resp_rate = NA), "resp_rate")
})

test_that("print shows net when resp_rate < 1", {
  result <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.8)
  out <- capture.output(print(result))
  expect_match(out[2], "net:")
  expect_match(out[2], "resp_rate = 0.80")
})

test_that("print hides net when resp_rate = 1", {
  result <- n_prop(p = 0.3, moe = 0.05)
  out <- capture.output(print(result))
  expect_no_match(out[2], "net:")
})

test_that("resp_rate near boundary (0.01) inflates heavily", {
  base <- n_prop(p = 0.3, moe = 0.05)
  rr <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.01)
  expect_equal(rr$n, base$n / 0.01, tolerance = 1e-6)
})

test_that("resp_rate = 0.99 barely inflates", {
  base <- n_prop(p = 0.3, moe = 0.05)
  rr <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.99)
  expect_equal(rr$n, base$n / 0.99, tolerance = 1e-6)
  expect_true(abs(rr$n - base$n) < 5)
})

test_that("prec_prop resp_rate deflates effective n", {
  base <- prec_prop(p = 0.3, n = 400)
  rr <- prec_prop(p = 0.3, n = 400, resp_rate = 0.8)
  expect_true(rr$se > base$se)
  n_eff_base <- 400
  n_eff_rr <- 400 * 0.8
  expect_equal(rr$se / base$se, sqrt(n_eff_base / n_eff_rr), tolerance = 1e-6)
})

test_that("prec_mean resp_rate deflates effective n", {
  base <- prec_mean(var = 100, n = 400)
  rr <- prec_mean(var = 100, n = 400, resp_rate = 0.8)
  expect_true(rr$se > base$se)
  expect_equal(rr$se / base$se, sqrt(1 / 0.8), tolerance = 1e-6)
})

test_that("n_prop and prec_prop round-trip with resp_rate", {
  s <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.8)
  p <- prec_prop(p = 0.3, n = s$n, resp_rate = 0.8)
  expect_equal(p$moe, 0.05, tolerance = 1e-6)
})

test_that("n_mean and prec_mean round-trip with resp_rate", {
  s <- n_mean(var = 100, moe = 2, resp_rate = 0.8)
  p <- prec_mean(var = 100, n = s$n, resp_rate = 0.8)
  expect_equal(p$moe, 2, tolerance = 1e-6)
})

test_that("resp_rate interacts correctly with deff", {
  base <- n_prop(p = 0.3, moe = 0.05)
  both <- n_prop(p = 0.3, moe = 0.05, deff = 2, resp_rate = 0.8)
  expect_equal(both$n, base$n * 2 / 0.8, tolerance = 1e-6)
})

test_that("resp_rate interacts correctly with FPC", {
  s1 <- n_prop(p = 0.3, moe = 0.05, N = 1000)
  s2 <- n_prop(p = 0.3, moe = 0.05, N = 1000, resp_rate = 0.8)
  expect_equal(s2$n, s1$n / 0.8, tolerance = 1e-6)
})

test_that("power_prop resp_rate works for solve-power", {
  base <- power_prop(p1 = 0.3, p2 = 0.35, n = 500, power = NULL)
  rr <- power_prop(p1 = 0.3, p2 = 0.35, n = 500, power = NULL,
                   resp_rate = 0.8)
  expect_true(rr$power < base$power)
})
