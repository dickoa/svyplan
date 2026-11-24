## Inputs a configuration cannot use are refused, the usable forms kept

test_that("single-stage indicator tables refuse cluster stage columns", {
  ind <- data.frame(name = "a", p = 0.3, cv = 0.1)
  for (col in c("icc_1", "var_ratio_1", "N_1", "take_2", "n_1")) {
    bad <- ind
    bad[[col]] <- 0.05
    expect_error(n_multi(bad), "describe cluster stages.*n_multi\\(\\)")
  }
  expect_error(prec_multi(transform(ind, cv = NULL, n = 400, icc_1 = 0.05)),
               "describe cluster stages.*prec_multi\\(\\)")
  expect_equal(n_multi(ind)$n, n_prop(p = 0.3, cv = 0.1)$n)
})

test_that("cluster tables refuse columns past the last stage and deff", {
  ind <- data.frame(name = "a", p = 0.3, cv = 0.1, icc_1 = 0.05)
  sc <- c(500, 50)
  expect_error(n_cluster(transform(ind, icc_2 = 0.3), stage_cost = sc),
               "'icc_2'.*past the last stage of this 2-stage")
  expect_error(n_cluster(transform(ind, var_ratio_2 = 2), stage_cost = sc),
               "past the last stage")
  expect_error(n_cluster(transform(ind, resp_rate_2 = 0.9), stage_cost = sc),
               "not applicable for 2-stage")
  expect_error(n_cluster(transform(ind, deff = 3), stage_cost = sc),
               "'deff' column")
  expect_error(prec_cluster(transform(ind, cv = NULL, n_1 = 20, take_2 = 10,
                                      icc_2 = 0.1)),
               "'icc_2'.*past the last stage")
  ok <- n_cluster(transform(ind, icc_2 = 0.1), stage_cost = c(500, 100, 20))
  expect_length(ok$n, 3L)
  expect_s3_class(prec_cluster(transform(ind, cv = NULL, n_1 = 20,
                                         take_2 = 10)), "svyplan_prec")
})

test_that("domain options without domains are refused", {
  ind <- data.frame(name = "a", p = 0.3, cv = 0.1, icc_1 = 0.05)
  sc <- c(500, 50)
  expect_error(n_cluster(ind, stage_cost = sc, allocation = "joint"),
               "needs 'domains' and 'budget'")
  dom <- rbind(transform(ind, region = "N"), transform(ind, region = "S"))
  expect_error(n_cluster(dom, stage_cost = sc, domains = "region",
                         allocation = "joint"),
               "needs 'budget'$")
  expect_error(n_cluster(ind, stage_cost = sc, min_n_domain = 100),
               "'min_n_domain' needs 'domains'")
  expect_error(n_multi(transform(ind, icc_1 = NULL), min_n_domain = 500),
               "'min_n_domain' needs 'domains'")
  j <- n_cluster(dom, stage_cost = sc, domains = "region",
                 allocation = "joint", budget = 60000)
  expect_equal(j$params$allocation, "joint")
})

test_that("an argument beside a frame column of the same name is refused", {
  fc <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12),
                   resp_rate = c(0.5, 0.9), deff = c(2, 1.5))
  expect_error(n_alloc(fc, n = 400, resp_rate = 0.8),
               "column 'resp_rate' has a value on every row")
  expect_error(n_alloc(fc, n = 400, deff = 3), "column 'deff' has a value")
  expect_error(n_alloc(fc, n = 400, plan = svyplan(resp_rate = 0.8)),
               "column 'resp_rate' has a value")
  expect_error(svyplan(deff = 3) |> n_alloc(fc, n = 400),
               "column 'deff' has a value")
  expect_error(prec_alloc(fc, n = c(150, 250), resp_rate = 0.8),
               "column 'resp_rate' has a value")
  expect_equal(n_alloc(fc, n = 400, resp_rate = 1)$n, n_alloc(fc, n = 400)$n)
  part <- transform(fc, resp_rate = c(0.5, NA))
  expect_equal(n_alloc(part, n = 400, resp_rate = 0.8)$n,
               n_alloc(fc[, -4], n = 400, resp_rate = c(0.5, 0.8))$n)
  expect_equal(n_alloc(fc, n = 400, resp_rate = c(0.7, 0.8))$n,
               n_alloc(fc[, -4], n = 400, resp_rate = c(0.7, 0.8))$n)
  f0 <- fc
  f0$resp_rate <- NULL
  expect_s3_class(n_alloc(f0, n = 400, resp_rate = 0.8), "svyplan_alloc")
})

test_that("predict() varies resp_rate and deff over a frame column", {
  fc <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12),
                   mean = c(50, 60), resp_rate = c(0.5, 0.9), deff = c(2, 1.5))
  fit <- n_alloc(fc, cv = 0.02)
  g <- predict(fit, data.frame(resp_rate = 0.7, deff = 1.2))
  f0 <- fc[, c("stratum", "N", "sd", "mean")]
  expect_equal(g$n, n_alloc(f0, cv = 0.02, resp_rate = 0.7, deff = 1.2)$n,
               tolerance = 1e-10)
  pr <- prec_alloc(fit)
  gp <- predict(pr, data.frame(resp_rate = 0.7))
  expect_equal(gp$cv, prec_alloc(transform(fc, resp_rate = NULL),
                                 n = pr$params$n, resp_rate = 0.7)$cv,
               tolerance = 1e-10)
})

test_that("eligibility arguments beside complete frame columns are refused", {
  fc <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12),
                   eligibility_rate = c(0.5, 0.9),
                   eligible_per_unit = c(2, 1.5))
  expect_error(n_alloc(fc, n = 400, eligibility_rate = 0.8),
               "column 'eligibility_rate' has a value")
  expect_error(n_alloc(fc, n = 400, eligible_per_unit = 3),
               "column 'eligible_per_unit' has a value")
  expect_error(prec_alloc(fc, n = c(150, 250), eligibility_rate = 0.8),
               "column 'eligibility_rate' has a value")
})

## An indicator row without a target

test_that("an indicator row without a target is refused", {
  tab <- data.frame(name = c("a", "b"), mean = 1, sd = 1, cv = c(NA, NA),
                    icc_1 = 0.05)
  expect_error(n_cluster(tab, stage_cost = c(100, 10), budget = 1e4),
               "each row must have either 'moe' or 'cv'")
  tab$cv <- c(0.1, NA)
  expect_error(n_cluster(tab, stage_cost = c(100, 10), budget = 1e4),
               "each row must have either 'moe' or 'cv'")
  ele <- data.frame(name = c("a", "b"), p = c(0.3, 0.4), cv = c(0.1, NA))
  expect_error(n_multi(ele), "each row must have either 'moe' or 'cv'")
})

## A classic precision assessment refuses an allocation the frame cannot hold

test_that("prec_alloc() refuses a stratum count above its population", {
  f <- data.frame(stratum = c("A", "B"), N = c(100, 200), sd = 1, mean = 1)
  bound <- "population bound for stratum: A$"
  expect_error(prec_alloc(f, n = c(110, 190)), bound)
  expect_error(prec_alloc(f, n = c(110, 150)), bound)
  expect_error(prec_alloc(f, n = c(110, 200)), bound)
  expect_error(prec_alloc(f, n = c(110, 150), fpc = "none"), bound)
  expect_error(prec_alloc(f, n = c(110, 150), resp_rate = 0.8), bound)
  expect_error(prec_alloc(f, n = c(110, 150), eligibility_rate = 0.5), bound)
  expect_error(prec_alloc(f, n = c(110, 150), eligible_per_unit = 2), bound)
  expect_s3_class(prec_alloc(f, n = c(80, 150), eligibility_rate = 0.5),
                  "svyplan_prec")
  expect_error(prec_alloc(f, n = c(120, 250)),
               "population bound for stratum: A, B$")
  expect_error(prec_alloc(f[, -1], n = c(110, 150)),
               "population bound for stratum: 1$")
  g <- data.frame(stratum = c("A", "B"), N = c(2e9, 100), sd = 1, mean = 1)
  expect_error(prec_alloc(g, n = c(2e9 + 1, 50)), bound)
  expect_s3_class(prec_alloc(g, n = c(2e9, 50)), "svyplan_prec")
})

test_that("prec_alloc() keeps a census and a fit at the bound", {
  f <- data.frame(stratum = c("A", "B"), N = c(100, 200), sd = 1, mean = 1)
  expect_equal(prec_alloc(f, n = c(100, 200))$se, 0)
  near <- prec_alloc(f, n = c(100, 200) * (1 + 1e-12))
  expect_equal(near$se, 0)
  expect_equal(near$detail$n_int, c(100, 200))
  expect_equal(prec_alloc(f, n = c(100, 200), fpc = "none")$se,
               0.0577350269, tolerance = 1e-8)
  expect_equal(prec_alloc(n_alloc(f, n = 300))$se, 0)
})

test_that("prec_alloc() refuses more PSUs than N_1 holds at the takes", {
  fc <- data.frame(stratum = c("A", "B"), N = c(1000, 2000),
                   N_1 = c(50, 100), take_2 = 10, icc_1 = 0.05,
                   sd = 1, mean = 1)
  expect_equal(prec_alloc(fc, n = c(500, 500))$se, 0.0335824, tolerance = 1e-6)
  psu <- "needs more PSUs than 'N_1' holds at its takes for stratum: A$"
  for (fpc in c("unit", "stage", "none")) {
    expect_error(prec_alloc(fc, n = c(1000, 500), fpc = fpc), psu)
  }
  expect_error(prec_alloc(fc, n = c(1100, 500)), "population bound")
  free <- data.frame(stratum = "A", N = 1000, N_1 = 50, icc_1 = 0.05,
                     cost_1 = 1, cost_2 = 1, sd = 1, mean = 1)
  rise <- prec_alloc(free, n = 980)
  expect_equal(rise$detail$take_2, 19.6)
  expect_equal(rise$detail$n_1, 50)
  room <- data.frame(stratum = "A", N = 40, N_1 = 10, take_2 = 4, take_3 = 2,
                     icc_1 = 0.05, icc_2 = 0.05, sd = 1, mean = 1)
  expect_warning(prec_alloc(room, n = 40), "puts more units in its PSUs")
  open <- data.frame(stratum = "A", N = 1000, take_2 = 10, icc_1 = 0.05,
                     sd = 1, mean = 1)
  expect_equal(prec_alloc(open, n = 200, fpc = "none")$se, 0.08514693,
               tolerance = 1e-6)
})

test_that("a precision grid fails only the rows above the population", {
  f <- data.frame(stratum = c("A", "B"), N = c(100, 200), sd = 1, mean = 1)
  pr <- prec_alloc(f, n = c(50, 80))
  expect_warning(
    g <- predict(pr, data.frame(n = I(list(c(50, 80), c(110, 150))))),
    "predict row 2 failed: 'n' exceeds the population bound for stratum: A"
  )
  expect_equal(g$se[1], pr$se)
  expect_true(is.na(g$se[2]))
  expect_warning(g <- predict(pr, data.frame(n = c(130, 400))),
                 "predict row 2 failed: 'n' exceeds the population bound")
  expect_equal(g$se[1], pr$se)
  expect_true(is.na(g$se[2]))
})
