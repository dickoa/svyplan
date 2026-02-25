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
  base <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05)
  rr <- n_cluster(cost = c(500, 50), delta = 0.05, cv = 0.05,
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
  base <- power_mean(delta = 5, var = 100, n = NULL, power = 0.8)
  rr <- power_mean(delta = 5, var = 100, n = NULL, power = 0.8,
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
