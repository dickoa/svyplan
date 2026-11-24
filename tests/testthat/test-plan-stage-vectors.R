## A profile's per-stage values reach only the functions that read stages

test_that("a per-stage profile value is refused where stages are not read", {
  p2 <- svyplan(stage_cost = c(500, 50), icc = 0.05, resp_rate = c(0.9, 0.8))
  fc <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12))
  expect_error(n_alloc(fc, n = 400, plan = p2),
               "one value per stage.*n_alloc\\(\\).*resp_rate_<k>")
  expect_error(p2 |> n_alloc(fc, n = 400), "one value per stage")
  expect_error(n_prop(p = 0.3, moe = 0.05, plan = p2),
               "one value per stage.*n_prop\\(\\)")
  expect_error(n_prop(p = 0.3, moe = 0.05, plan = svyplan(N = c(50, 5000))),
               "'N' one value per stage")
})

test_that("cluster functions still read a profile by stage", {
  p2 <- svyplan(stage_cost = c(500, 50), icc = 0.05, resp_rate = c(0.9, 0.8))
  a <- n_cluster(cv = 0.05, plan = p2)
  b <- n_cluster(stage_cost = c(500, 50), icc = 0.05, resp_rate = c(0.9, 0.8),
                 cv = 0.05)
  expect_equal(a$n, b$n)
  expect_equal(prec_cluster(n = b$n, plan = p2)$cv,
               prec_cluster(n = b$n, stage_cost = c(500, 50), icc = 0.05,
                            resp_rate = c(0.9, 0.8))$cv)
  fc <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12))
  expect_equal(n_alloc(fc, n = 400, plan = svyplan(resp_rate = 0.8))$n,
               n_alloc(fc, n = 400, resp_rate = 0.8)$n)
})

test_that("explicit scalar values override unused stage vectors in a profile", {
  profile <- svyplan(resp_rate = c(0.9, 0.8), N = c(100, 1000), deff = 1.5)
  for (N in c(5000, Inf)) {
    from_profile <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.85, N = N,
                           plan = profile)
    direct <- n_prop(p = 0.3, moe = 0.05, resp_rate = 0.85, N = N, deff = 1.5)
    expect_equal(from_profile$n, direct$n)
  }
  expect_error(n_prop(p = 0.3, moe = 0.05, resp_rate = 0.85, plan = profile),
               "'N' one value per stage")
})
