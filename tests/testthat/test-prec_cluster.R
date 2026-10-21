test_that("prec_cluster computes CV for 2-stage", {
  result <- prec_cluster(n = c(50, 12), icc = 0.05)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "cluster")

  cv_exp <- sqrt(1 * 1 / (50 * 12) * (1 + 0.05 * (12 - 1)))
  expect_equal(result$cv, cv_exp, tolerance = 1e-6)
  expect_true(is.na(result$se))
  expect_true(is.na(result$moe))
})

test_that("prec_cluster computes CV for 3-stage", {
  result <- prec_cluster(n = c(50, 12, 8), icc = c(0.01, 0.05))
  expect_s3_class(result, "svyplan_prec")

  cv_exp <- sqrt(1 / (50 * 12 * 8) *
    (1 * 0.01 * 12 * 8 + (1 - 0.01) * (1 + 0.05 * (8 - 1))))
  expect_equal(result$cv, cv_exp, tolerance = 1e-6)
})

test_that("prec_cluster with resp_rate_psu deflates stage-1", {
  base <- prec_cluster(n = c(50, 12), icc = 0.05)
  rr <- prec_cluster(n = c(50, 12), icc = 0.05, resp_rate_psu = 0.8)
  cv_exp <- sqrt(1 / (40 * 12) * (1 + 0.05 * 11))
  expect_equal(rr$cv, cv_exp, tolerance = 1e-6)
  expect_true(rr$cv > base$cv)
})

test_that("prec_cluster.svyplan_cluster carries stage_cost metadata", {
  s1 <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
  p1 <- prec_cluster(s1)
  expect_equal(unname(p1$params$stage_cost), c(500, 50))
  expect_equal(names(p1$params$stage_cost), c("cost_psu", "cost_ssu"))
  expect_equal(p1$params$budget, 100000)
})

test_that("prec_cluster validates inputs", {
  expect_error(prec_cluster(n = 10, icc = 0.05), "length >= 2")
  expect_error(prec_cluster(n = c(50, -1), icc = 0.05), "positive")
  expect_error(prec_cluster(n = c(50, 12), icc = c(0.05, 0.1)),
               "must have length")
})

test_that("prec_cluster rejects NA in n", {
  expect_error(prec_cluster(n = c(50, NA), icc = 0.05), "NA")
})

test_that("prec_cluster rejects infinite n", {
  expect_error(
    prec_cluster(n = c(Inf, 10), icc = 0.05),
    "finite"
  )
})

test_that("prec_cluster prints cluster format", {
  result <- prec_cluster(n = c(50, 12), icc = 0.05)
  out <- capture.output(print(result))
  expect_match(out[1], "Sampling precision for 2-stage cluster")
  expect_match(out[2], "n_psu = 50")
  expect_match(out[2], "n_per_psu = 12")
  expect_match(out[3], "cv =")
})

test_that("prec_cluster rejects negative unit_relvar", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, unit_relvar = -1),
    "unit_relvar.*positive"
  )
})

test_that("prec_cluster rejects zero unit_relvar", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, unit_relvar = 0),
    "unit_relvar.*positive"
  )
})

test_that("prec_cluster rejects negative var_ratio", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, var_ratio = -1),
    "var_ratio.*positive"
  )
})

test_that("prec_cluster rejects NA var_ratio", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, var_ratio = NA),
    "var_ratio.*positive"
  )
})

test_that("prec_cluster accepts named n vectors and reorders", {
  ref2 <- prec_cluster(n = c(50, 12), icc = 0.05)
  named2 <- prec_cluster(n = c(n_per_psu = 12, n_psu = 50), icc = 0.05)
  alias2 <- prec_cluster(n = c(n_per_psu = 12, n_psu = 50), icc = 0.05)
  expect_equal(named2$cv, ref2$cv, tolerance = 1e-10)
  expect_equal(alias2$cv, ref2$cv, tolerance = 1e-10)
  expect_equal(names(named2$params$n), c("n_psu", "n_per_psu"))

  ref3 <- prec_cluster(n = c(50, 12, 8), icc = c(0.01, 0.05))
  named3 <- prec_cluster(
    n = c(n_per_ssu = 8, n_psu = 50, n_per_psu = 12),
    icc = c(0.01, 0.05)
  )
  alias3 <- prec_cluster(
    n = c(n_per_ssu = 8, n_psu = 50, n_per_psu = 12),
    icc = c(0.01, 0.05)
  )
  expect_equal(named3$cv, ref3$cv, tolerance = 1e-10)
  expect_equal(alias3$cv, ref3$cv, tolerance = 1e-10)
  expect_equal(names(named3$params$n), c("n_psu", "n_per_psu", "n_per_ssu"))

  expect_error(
    prec_cluster(n = c(n1 = 50, n2 = 12), icc = 0.05),
    "unrecognized names"
  )
})

## The two entry modes refuse each other's arguments

test_that("the several-indicators mode refuses the scalar sizes and moments", {
  ind <- data.frame(name = "y", p = 0.3, n = 60, n_per_psu = 12,
                    icc_psu = 0.05)
  refused <- list(n = c(60, 12), icc = 0.05, unit_relvar = 2, var_ratio = 1.2)
  columns <- c(n = "n", icc = "icc_psu", unit_relvar = "unit_relvar",
               var_ratio = "var_ratio_psu")
  for (nm in names(refused)) {
    expect_error(
      do.call(prec_cluster, c(list(indicators = ind), refused[nm])),
      sprintf("'%s' is carried by the 'indicators' column '%s'",
              nm, columns[[nm]]),
      info = nm
    )
  }
})

test_that("the scalar mode refuses the several-indicators arguments", {
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, domains = "region"),
    "'domains' applies to 'indicators'"
  )
})

test_that("an indicator table handed to the first slot is named", {
  ind <- data.frame(name = "y", p = 0.3, n = 60, n_per_psu = 12,
                    icc_psu = 0.05)
  expect_error(
    prec_cluster(ind, icc = 0.05),
    "table of indicators goes to 'indicators'"
  )
})

test_that("scalar costs are recorded for the round trip", {
  p <- prec_cluster(n = c(50, 12), icc = 0.05, stage_cost = c(500, 50))
  expect_equal(unname(p$params$stage_cost), c(500, 50))
  expect_s3_class(n_cluster(p), "svyplan_cluster")
  expect_error(
    prec_cluster(n = c(50, 12), icc = 0.05, stage_cost = c(500, 50, 20)),
    "must have length 2"
  )
})

## Expected takes and printed counts

test_that("the indicator table refuses a take below one response per row", {
  expect_error(
    prec_cluster(c(100, 2), icc = 0.05, resp_rate = 0.1),
    "expects 0.2 responses per unit"
  )
  expect_error(
    prec_cluster(indicators = data.frame(name = "a", p = 0.5, n = 100,
      n_per_psu = 2, resp_rate = 0.1, icc_psu = 0.05)),
    "'resp_rate' = 0.1 for indicator 'a' expects 0.2"
  )
  expect_error(
    n_cluster(indicators = data.frame(name = "a", p = 0.5, cv = 0.1,
      resp_rate = 0.1, icc_psu = 0.05), stage_cost = c(500, 50),
      n_per_psu = 2),
    "for indicator 'a' expects 0.2"
  )
  expect_error(
    prec_cluster(indicators = data.frame(name = c("a", "b"), p = 0.5,
      n = 100, n_per_psu = 4, n_per_ssu = 5, resp_rate_ssu = c(1, 0.2),
      icc_psu = 0.05, icc_ssu = 0.1)),
    "'n_per_psu' = 4 at 'resp_rate_ssu' = 0.2 for indicator 'b'"
  )
  expect_no_error(
    prec_cluster(indicators = data.frame(name = "a", p = 0.5, n = 100,
      n_per_psu = 10, resp_rate = 0.1, icc_psu = 0.05))
  )
})

test_that("the printed net count carries every stage rate", {
  out <- capture.output(print(prec_cluster(c(100, 20), icc = 0.05,
                                           resp_rate_psu = 0.5,
                                           resp_rate = 0.5)))
  expect_match(out[2], "total n = 2,?000 \\(net: 500\\)")
  out3 <- capture.output(print(prec_cluster(c(10, 5, 4),
                                            icc = c(0.05, 0.1),
                                            resp_rate_ssu = 0.5)))
  expect_match(out3[2], "\\(net: 100\\)")
})

test_that("the table and scalar cluster paths agree where a take falls below one", {
  cases <- list(
    list(stage_cost = c(500, 50), icc = 0.05, ind = list(icc_psu = 0.05),
         cv = 0.05, n_psu = 2000),
    list(stage_cost = c(500, 50, 5), icc = c(0.05, 0.1),
         ind = list(icc_psu = 0.05, icc_ssu = 0.1), cv = 0.1, n_psu = 80,
         n_per_ssu = 5),
    list(stage_cost = c(500, 50, 5), icc = c(0.05, 0.1),
         ind = list(icc_psu = 0.05, icc_ssu = 0.1), cv = 0.1, n_psu = 80,
         n_per_psu = 10),
    list(stage_cost = c(500, 50, 5), icc = c(0.8, 0.1),
         ind = list(icc_psu = 0.8, icc_ssu = 0.1), cv = 0.05)
  )
  for (k in cases) {
    fixed <- k[intersect(c("n_psu", "n_per_psu", "n_per_ssu"), names(k))]
    s <- suppressWarnings(do.call(n_cluster, c(
      list(stage_cost = k$stage_cost, icc = k$icc, cv = k$cv,
           unit_relvar = 7 / 3), fixed)))
    m <- do.call(n_cluster, c(
      list(indicators = do.call(data.frame, c(list(p = 0.3, cv = k$cv),
                                              k$ind)),
           stage_cost = k$stage_cost), fixed))
    expect_equal(unname(m$n), unname(s$n), tolerance = 1e-6)
    expect_equal(m$cv, s$cv, tolerance = 1e-6)
    expect_equal(m$cost, s$cost, tolerance = 1e-6)
  }
})

test_that("a supplied fractional take is refused at every stage", {
  expect_error(
    prec_cluster(indicators = data.frame(name = "p", p = 0.5, n = 100,
                                         n_per_psu = 0.5, icc_psu = 0.05)),
    "'n_per_psu' = 0.5 at 'resp_rate' = 1 for indicator 'p'"
  )
  expect_error(
    prec_cluster(indicators = data.frame(name = "p", p = 0.5, n = 100,
      n_per_psu = 4, n_per_ssu = 0.5, icc_psu = 0.05, icc_ssu = 0.1)),
    "'n_per_ssu' = 0.5"
  )
  expect_error(
    prec_cluster(indicators = data.frame(name = "p", p = 0.5, n = 100,
      n_per_psu = 0.5, n_per_ssu = 4, icc_psu = 0.05, icc_ssu = 0.1)),
    "'n_per_psu' = 0.5 at 'resp_rate_ssu'"
  )
})
