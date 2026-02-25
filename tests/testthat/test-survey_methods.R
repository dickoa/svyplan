test_that("SE and cv work via $ access (no survey needed)", {
  s1 <- n_prop(p = 0.3, moe = 0.05)
  expect_false(is.na(s1$se))
  expect_false(is.na(s1$cv))
  expect_false(is.na(s1$moe))

  p1 <- prec_prop(p = 0.3, n = 400)
  expect_false(is.na(p1$se))
  expect_false(is.na(p1$cv))

  c1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_true(is.na(c1$se))
  expect_false(is.na(c1$cv))
})

test_that("svyplan_n stores consistent precision measures", {
  s1 <- n_prop(p = 0.3, moe = 0.05)
  z <- qnorm(0.975)
  expect_equal(s1$moe, z * s1$se, tolerance = 1e-6)
  expect_equal(s1$cv, s1$se / 0.3, tolerance = 1e-6)
})

test_that("svyplan_n for mean stores consistent precision measures", {
  s1 <- n_mean(var = 100, moe = 2, mu = 50)
  z <- qnorm(0.975)
  expect_equal(s1$moe, z * s1$se, tolerance = 1e-6)
  expect_equal(s1$cv, s1$se / 50, tolerance = 1e-6)
})

test_that("svyplan_n for mean without mu gives NA cv", {
  s1 <- n_mean(var = 100, moe = 2)
  expect_false(is.na(s1$se))
  expect_false(is.na(s1$moe))
  expect_true(is.na(s1$cv))
})

test_that("svyplan_n for multi gives NA precision", {
  tgt <- data.frame(p = c(0.3, 0.5), moe = c(0.05, 0.05))
  s1 <- n_multi(tgt)
  expect_true(is.na(s1$se))
  expect_true(is.na(s1$moe))
  expect_true(is.na(s1$cv))
})

test_that("survey::SE and survey::cv work when survey is loaded", {
  skip_if_not_installed("survey")

  s1 <- n_prop(p = 0.3, moe = 0.05)
  expect_equal(survey::SE(s1), s1$se)
  expect_equal(survey::cv(s1), s1$cv)

  p1 <- prec_prop(p = 0.3, n = 400)
  expect_equal(survey::SE(p1), p1$se)
  expect_equal(survey::cv(p1), p1$cv)

  c1 <- n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
  expect_true(is.na(survey::SE(c1)))
  expect_equal(survey::cv(c1), c1$cv)
})
