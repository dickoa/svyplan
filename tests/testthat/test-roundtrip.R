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
  s1 <- n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  expect_equal(unname(p1$cv), unname(s1$cv), tolerance = 1e-6)
})

test_that("n_cluster -> prec_cluster -> n_cluster round-trip", {
  s1 <- n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  s2 <- n_cluster(p1)
  expect_equal(unname(s2$cv), unname(s1$cv), tolerance = 1e-4)
})

test_that("n_cluster 3-stage -> prec_cluster round-trip", {
  s1 <- n_cluster(stage_cost = c(500, 100, 50), delta = c(0.01, 0.05), cv = 0.05)
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

test_that("n_multi.svyplan_prec forwards prop_method overrides", {
  p1 <- prec_multi(data.frame(p = 0.05, n = 400), prop_method = "wilson")
  s_wilson <- n_multi(p1)
  s_wald <- n_multi(p1, prop_method = "wald")

  expect_equal(s_wilson$n, 400, tolerance = 1e-6)
  expect_equal(
    s_wald$n,
    n_prop(p = 0.05, moe = p1$detail$.moe[1], method = "wald")$n,
    tolerance = 1e-6
  )
})

test_that("n_multi 2-stage -> prec_multi round-trip", {
  tgt <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  s1 <- n_multi(tgt, stage_cost = c(500, 50))
  p1 <- prec_multi(s1)
  expect_s3_class(p1, "svyplan_prec")
  expect_equal(nrow(p1$detail), 2L)
})

test_that("prec_multi.svyplan_n forwards prop_method overrides", {
  s1 <- n_multi(data.frame(p = 0.05, moe = 0.02), prop_method = "wilson")
  p_wilson <- prec_multi(s1)
  p_wald <- prec_multi(s1, prop_method = "wald")

  expect_equal(p_wilson$detail$.moe[1], 0.02, tolerance = 1e-6)
  expect_equal(
    p_wald$detail$.moe[1],
    prec_prop(p = 0.05, n = s1$n, method = "wald")$moe,
    tolerance = 1e-6
  )
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
  s1 <- n_cluster(stage_cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_no_error(prec_cluster(s1))
})

test_that("n_multi 2-stage -> prec_multi -> n_multi preserves cluster class", {
  tgt <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  s1 <- n_multi(tgt, stage_cost = c(500, 50))
  p1 <- prec_multi(s1)
  s2 <- n_multi(p1)
  expect_s3_class(s2, "svyplan_cluster")
  expect_equal(s2$stages, s1$stages)
})

test_that("n_multi budget round-trip preserves budget param", {
  tgt <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  s1 <- n_multi(tgt, stage_cost = c(500, 50), budget = 100000)
  p1 <- prec_multi(s1)
  expect_equal(p1$params$budget, 100000)
  s2 <- n_multi(p1)
  expect_equal(s2$params$budget, 100000)
})

test_that("n_multi round-trip does not leak n2/n3 as domain columns", {
  tgt <- data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    cv = c(0.10, 0.15),
    delta_psu = c(0.02, 0.05)
  )
  s1 <- n_multi(tgt, stage_cost = c(500, 50))
  p1 <- prec_multi(s1)
  expect_no_message(s2 <- n_multi(p1))
})

test_that("n_prop.svyplan_prec respects explicit method override", {
  p_wilson <- prec_prop(p = 0.3, n = 400, method = "wilson")
  s_same <- n_prop(p_wilson)
  s_wald <- n_prop(p_wilson, method = "wald")
  s_logodds <- n_prop(p_wilson, method = "logodds")

  expect_equal(s_same$method, "wilson")
  expect_equal(s_same$n, 400, tolerance = 1e-6)

  expect_equal(s_wald$method, "wald")
  expect_true(s_wald$n != s_same$n)

  expect_equal(s_logodds$method, "logodds")
  expect_true(s_logodds$n != s_same$n)
})

test_that("prec_prop.svyplan_n respects explicit method override", {
  s_wilson <- n_prop(p = 0.3, moe = 0.05, method = "wilson")
  p_same <- prec_prop(s_wilson)
  p_wald <- prec_prop(s_wilson, method = "wald")
  p_logodds <- prec_prop(s_wilson, method = "logodds")

  expect_equal(p_same$method, "wilson")
  expect_equal(p_same$moe, 0.05, tolerance = 1e-6)

  expect_equal(p_wald$method, "wald")
  expect_equal(p_wald$params$n, s_wilson$n)
  expect_true(abs(p_wald$moe - 0.05) > 1e-4)

  expect_equal(p_logodds$method, "logodds")
  expect_equal(p_logodds$params$n, s_wilson$n)
  expect_true(abs(p_logodds$moe - 0.05) > 1e-4)
})
