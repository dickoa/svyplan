test_that("effective_n from weights is the Kish effective sample size", {
  set.seed(1009)
  w <- runif(100, 1, 5)
  expect_equal(effective_n(weights = w), sum(w)^2 / sum(w^2))
})

test_that("effective_n equals n when weights are equal", {
  w <- rep(2.5, 60)
  expect_equal(effective_n(weights = w), 60)
})

test_that("effective_n is n divided by design_effect", {
  set.seed(444)
  w <- runif(80, 1, 6)
  expect_equal(
    effective_n(weights = w),
    length(w) / as.double(design_effect(weights = w))
  )
  expect_equal(
    effective_n(n = 1200, icc = 0.05, n_per_psu = 25),
    1200 / as.double(design_effect(icc = 0.05, n_per_psu = 25))
  )
})

test_that("effective_n accepts a design effect object or a plain number", {
  deff <- design_effect(icc = 0.05, n_per_psu = 25)
  expect_equal(effective_n(deff, n = 1200), 1200 / 2.2)
  expect_equal(effective_n(n = 1200, deff = deff), 1200 / 2.2)
  expect_equal(effective_n(n = 1200, deff = 2.2), 1200 / 2.2)
})

test_that("deff and its components are mutually exclusive", {
  expect_error(
    effective_n(n = 100, deff = 2, icc = 0.05, n_per_psu = 10),
    "either 'deff' or the design components"
  )
})

test_that("effective_n derives n where it can and requires it otherwise", {
  strata <- data.frame(N = c(50000, 120000), n = c(600, 400))
  expect_equal(
    effective_n(strata = strata),
    1000 / as.double(design_effect(strata = strata))
  )
  expect_error(effective_n(icc = 0.05, n_per_psu = 25), "'n' is required")
  expect_error(
    effective_n(design_effect(icc = 0.05, n_per_psu = 25)), "'n' is required"
  )
})

test_that("cluster planning bounds behave at the extremes", {
  expect_equal(effective_n(n = 500, icc = 0, n_per_psu = 30), 500)
  expect_equal(effective_n(n = 500, icc = 1, n_per_psu = 20), 25)
  expect_equal(effective_n(n = 500, icc = 0.05, n_per_psu = 1), 500)
})

test_that("effective_n reads plans and allocations", {
  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  expect_equal(effective_n(plan), plan$total_n / as.double(design_effect(plan)))
  expect_equal(effective_n(plan, n = 1000), 1000 / as.double(design_effect(plan)))

  alloc <- n_alloc(
    data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 15),
               mean = c(50, 60)),
    n = 500
  )
  expect_equal(effective_n(alloc), alloc$n / as.double(design_effect(alloc)))
})

test_that("effective_n validates its inputs", {
  expect_error(effective_n(n = -5, icc = 0.05, n_per_psu = 10), "positive")
  expect_error(effective_n(weights = c(1, Inf)), "only finite")
  expect_error(effective_n(weights = c(1, -Inf)), "only finite")
  expect_error(effective_n(n = 100, icc = 0.05), "both 'icc' and 'n_per_psu'")
  expect_error(effective_n(n = 100, icc = 1.5, n_per_psu = 10), "in \\[0, 1\\]")
})

test_that("effective_n rejects unused arguments", {
  expect_error(effective_n(n = 100, icc = 0.05, n_per_psu = 10,
                           methd = "kish"),
               "unused argument.*methd")
})

test_that("effective_n nets down for response, matching the planning identity", {
  frame <- data.frame(
    stratum = c("A", "B"),
    N = c(1e5, 1e5),
    sd = c(1, 1),
    mean = c(1, 1)
  )
  alloc <- n_alloc(frame, n = 200, alloc = "neyman", resp_rate = 0.8)
  deff <- as.double(design_effect(alloc))

  # n * resp_rate / deff, the identity the n_eff column reports
  expect_equal(as.double(effective_n(alloc)), 200 * 0.8 / deff)
  expect_equal(as.double(effective_n(alloc)), sum(alloc$detail$n_eff))

  # an explicit rate overrides the plan's
  expect_equal(as.double(effective_n(alloc, resp_rate = 1)), 200 / deff)

  # a plan without nonresponse is unchanged
  full <- n_alloc(frame, n = 200, alloc = "neyman")
  expect_equal(as.double(effective_n(full)),
               200 / as.double(design_effect(full)))
})

test_that("effective_n takes a response rate for a bare n", {
  expect_equal(effective_n(n = 1000, deff = 2), 500)
  expect_equal(effective_n(n = 1000, deff = 2, resp_rate = 0.8), 400)
  expect_error(effective_n(n = 1000, deff = 2, resp_rate = 0), "resp_rate")
})

test_that("a cluster plan carries its own response rate", {
  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                    resp_rate = 0.8)
  expect_equal(
    as.double(effective_n(plan)),
    plan$total_n * 0.8 / as.double(design_effect(plan))
  )
})
