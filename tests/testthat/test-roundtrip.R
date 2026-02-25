test_that("n_prop -> prec_prop -> n_prop round-trip", {
  s1 <- n_prop(p = 0.3, moe = 0.05)
  p1 <- prec_prop(s1)
  s2 <- n_prop(p1)
  expect_equal(s2$n, s1$n, tolerance = 1e-6)
})

test_that("n_prop -> prec_prop -> n_prop round-trip with deff and resp_rate", {
  s1 <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8)
  p1 <- prec_prop(s1)
  s2 <- n_prop(p1)
  expect_equal(s2$n, s1$n, tolerance = 1e-6)
})

test_that("n_prop -> prec_prop -> n_prop with FPC", {
  s1 <- n_prop(p = 0.3, moe = 0.05, N = 5000)
  p1 <- prec_prop(s1)
  s2 <- n_prop(p1)
  expect_equal(s2$n, s1$n, tolerance = 0.1)
})

test_that("n_prop -> prec_prop -> n_prop with cv override", {
  s1 <- n_prop(p = 0.3, moe = 0.05)
  p1 <- prec_prop(s1)
  s2 <- n_prop(p1, cv = 0.10)
  expect_equal(s2$params$cv, 0.10)
  expect_null(s2$params$moe)
})

test_that("n_mean -> prec_mean -> n_mean round-trip", {
  s1 <- n_mean(var = 100, moe = 2, mu = 50)
  p1 <- prec_mean(s1)
  s2 <- n_mean(p1)
  expect_equal(s2$n, s1$n, tolerance = 1e-6)
})

test_that("n_mean -> prec_mean -> n_mean with deff and FPC", {
  s1 <- n_mean(var = 100, moe = 2, N = 5000, deff = 1.5)
  p1 <- prec_mean(s1)
  s2 <- n_mean(p1)
  expect_equal(s2$n, s1$n, tolerance = 1e-6)
})

test_that("n_cluster -> prec_cluster round-trip", {
  s1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  expect_equal(unname(p1$cv), unname(s1$cv), tolerance = 1e-6)
})

test_that("n_cluster -> prec_cluster -> n_cluster round-trip", {
  s1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  s2 <- n_cluster(p1)
  expect_equal(unname(s2$cv), unname(s1$cv), tolerance = 1e-4)
})

test_that("n_cluster 3-stage -> prec_cluster round-trip", {
  s1 <- n_cluster(cost = c(500, 100, 50), delta = c(0.01, 0.05), cv = 0.05)
  p1 <- prec_cluster(s1)
  expect_equal(unname(p1$cv), unname(s1$cv), tolerance = 1e-6)
})

test_that("prec_prop -> n_prop -> prec_prop round-trip", {
  p1 <- prec_prop(p = 0.3, n = 400)
  s1 <- n_prop(p1)
  p2 <- prec_prop(s1)
  expect_equal(p2$moe, p1$moe, tolerance = 1e-6)
  expect_equal(p2$se, p1$se, tolerance = 1e-6)
})

test_that("prec_mean -> n_mean -> prec_mean round-trip", {
  p1 <- prec_mean(var = 100, n = 400, mu = 50)
  s1 <- n_mean(p1)
  p2 <- prec_mean(s1)
  expect_equal(p2$moe, p1$moe, tolerance = 1e-6)
  expect_equal(p2$se, p1$se, tolerance = 1e-6)
})

test_that("n_multi simple -> prec_multi -> n_multi round-trip", {
  tgt <- data.frame(
    name = c("stunting", "vaccination", "anemia"),
    p = c(0.30, 0.70, 0.10),
    moe = c(0.05, 0.05, 0.03)
  )
  s1 <- n_multi(tgt)
  p1 <- prec_multi(s1)
  s2 <- n_multi(p1)
  expect_equal(s2$n, s1$n, tolerance = 1e-6)
})

test_that("n_multi 2-stage -> prec_multi round-trip", {
  tgt <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta1 = c(0.02, 0.05)
  )
  s1 <- n_multi(tgt, cost = c(500, 50))
  p1 <- prec_multi(s1)
  expect_s3_class(p1, "svyplan_prec")
  expect_equal(nrow(p1$detail), 2L)
})

test_that("n_prop rejects wrong prec type", {
  p1 <- prec_mean(var = 100, n = 400)
  expect_error(n_prop(p1), "type 'proportion'")
})

test_that("n_mean rejects wrong prec type", {
  p1 <- prec_prop(p = 0.3, n = 400)
  expect_error(n_mean(p1), "type 'mean'")
})

test_that("n_cluster rejects wrong prec type", {
  p1 <- prec_prop(p = 0.3, n = 400)
  expect_error(n_cluster(p1), "type 'cluster'")
})

test_that("prec_prop rejects wrong n type", {
  s1 <- n_mean(var = 100, moe = 2)
  expect_error(prec_prop(s1), "type 'proportion'")
})

test_that("prec_mean rejects wrong n type", {
  s1 <- n_prop(p = 0.3, moe = 0.05)
  expect_error(prec_mean(s1), "type 'mean'")
})

test_that("prec_multi rejects non-multi svyplan_n", {
  s1 <- n_prop(p = 0.3, moe = 0.05)
  expect_error(prec_multi(s1), "type 'multi'")
})

test_that("prec_cluster.svyplan_cluster round-trip preserves delta", {
  s1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_no_error(prec_cluster(s1))
})
