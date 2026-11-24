## Stage columns named by unit are refused at every table entry point

test_that("a frame with unit-named stage columns is refused", {
  fr <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12),
                   N_psu = c(80, 100), cost_psu = c(120, 200),
                   cost_ssu = c(15, 12), n_per_psu = c(10, 12),
                   icc_psu = c(0.05, 0.05))
  expect_error(n_alloc(fr, n = 400),
               "'frame' has .*icc_psu.*numbered by stage.*icc_1.*for.*icc_psu")
  expect_error(prec_alloc(fr, n = c(150, 250)), "numbered by stage")
  ok <- data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12),
                   N_1 = c(80, 100), cost_1 = c(120, 200),
                   cost_2 = c(15, 12), take_2 = c(10, 12),
                   icc_1 = c(0.05, 0.05))
  expect_s3_class(n_alloc(ok, n = 400), "svyplan_alloc")
})

test_that("joint tables with unit-named stage columns are refused", {
  frame <- data.frame(stratum = c("A", "B"), N = c(4000, 6000),
                      N_1 = c(80, 100), cost_1 = c(120, 200),
                      cost_2 = c(15, 12), take_2 = c(10, 12))
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = c(0.4, 0.3),
                         icc_1 = c(0.06, 0.04))
  targets <- data.frame(name = "y", cv = 0.06)
  bad_m <- measures
  names(bad_m)[names(bad_m) == "icc_1"] <- "icc_psu"
  expect_error(n_alloc_joint(frame, bad_m, targets = targets),
               "'measures' has .*icc_psu")
  bad_f <- frame
  names(bad_f)[names(bad_f) == "take_2"] <- "n_per_psu"
  expect_error(n_alloc_joint(bad_f, measures, targets = targets),
               "'frame' has .*n_per_psu.*take_2")
  expect_error(prec_alloc_joint(bad_f, n = c(150, 240), measures = measures,
                                targets = targets),
               "'frame' has .*n_per_psu")
  psu <- data.frame(psu_id = 1:4, stratum = c("A", "A", "B", "B"),
                    N = c(2000, 2000, 3000, 3000), N_psu = 10)
  expect_error(n_alloc_joint(frame, measures, targets = targets, psu = psu),
               "'psu' has .*N_psu")
})

test_that("indicator tables with unit-named stage columns are refused", {
  ind <- data.frame(name = "a", p = 0.3, cv = 0.1, icc_psu = 0.05)
  expect_error(n_cluster(ind, stage_cost = c(500, 50)), "'x' has .*icc_psu")
  expect_error(n_multi(ind), "'indicators' has .*icc_psu")
  expect_error(prec_multi(transform(ind, cv = NULL, n = 400)),
               "'indicators' has .*icc_psu")
  expect_error(
    prec_cluster(data.frame(name = "a", p = 0.3, icc_1 = 0.05,
                            var_ratio_psu = 1, n_1 = 20, take_2 = 10)),
    "'x' has .*var_ratio_psu"
  )
})

## Names on stage vectors place the values

stored_rates <- function(...) {
  c(resp_rate_1 = ..1, resp_rate_2 = ..2, resp_rate = ..3)
}

test_that("named stage vectors are put in stage order", {
  sc <- c(500, 100, 20)
  icc <- c(0.05, 0.1)
  a <- n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                 stage_n = c(stage_2 = 4, stage_1 = NA, stage_3 = NA))
  b <- n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                 stage_n = c(NA, 4, NA))
  expect_equal(a$n, b$n)
  expect_equal(unname(a$n[2]), 4)
  r <- n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                 resp_rate = c(stage_3 = 0.6, stage_2 = 1, stage_1 = 1))
  expect_equal(r$params$resp_rate, stored_rates(1, 1, 0.6))
  r2 <- n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                  resp_rate = c(resp_rate = 0.6, resp_rate_2 = 0.9,
                                resp_rate_1 = 1))
  expect_equal(r2$params$resp_rate, stored_rates(1, 0.9, 0.6))
  named <- c(district = 500, ea = 100, hh = 20)
  r3 <- n_cluster(stage_cost = named, icc = icc, cv = 0.05,
                  resp_rate = c(hh = 0.6, ea = 1, district = 1),
                  N = c(N = 12000, N_1 = 40, N_2 = 400))
  r4 <- n_cluster(stage_cost = named, icc = icc, cv = 0.05,
                  resp_rate = c(1, 1, 0.6), N = c(40, 400, 12000))
  expect_equal(r3$n, r4$n)
  p <- prec_cluster(n = c(20, 4, 10), icc = icc,
                    resp_rate = c(stage_3 = 0.6, stage_1 = 1, stage_2 = 1))
  expect_equal(p$params$resp_rate, stored_rates(1, 1, 0.6))
  expect_equal(n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                         resp_rate = c(final = 0.8))$params$resp_rate,
               stored_rates(1, 1, 0.8))
})

test_that("unknown names on stage vectors are refused", {
  sc <- c(500, 100, 20)
  icc <- c(0.05, 0.1)
  expect_error(n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                         resp_rate = c(x = 0.6, y = 1, z = 1)),
               "names of 'resp_rate' must match the stages")
  expect_error(n_cluster(stage_cost = sc, icc = icc, cv = 0.05,
                         stage_n = c(a = NA, b = 4, c = NA)),
               "names of 'stage_n'")
  expect_error(prec_cluster(n = c(20, 4, 10), icc = icc,
                            N = c(stage_1 = 40, stage_1 = 400, stage_3 = 9e3)),
               "names of 'N'")
})
