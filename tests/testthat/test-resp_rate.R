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

test_that("n_cluster resp_rate_psu inflates stage-1 in CV mode", {
  base <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  rr <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                   resp_rate_psu = 0.8)
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
  expect_match(out[2], "gross")
  expect_match(out[2], "net:")
  expect_match(out[2], "resp_rate = 0.80")
})

test_that("single-size print labels gross n and handles large finite counts", {
  result <- n_prop(p = 0.5, moe = 0.05, deff = 1.8, resp_rate = 0.1)
  out <- capture.output(print(result))
  expect_match(out[2], "n = 6915 gross \\(net: 692\\)")

  huge <- n_prop(p = 0.5, moe = 0.05, resp_rate = 1e-9)
  expect_no_error(capture.output(print(huge)))
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

test_that("n_multi accepts a scalar response-rate default", {
  indicators <- data.frame(
    name = c("a", "b"),
    p = c(0.3, 0.5),
    moe = c(0.05, 0.04)
  )
  scalar <- n_multi(indicators, resp_rate = 0.8)
  column <- n_multi(transform(indicators, resp_rate = 0.8))
  expect_equal(scalar$n, column$n)
  expect_equal(scalar$detail$.n, column$detail$.n)
  expect_equal(scalar$indicators$resp_rate, c(0.8, 0.8))
})

test_that("multi response-rate columns override the scalar row by row", {
  indicators <- data.frame(
    p = c(0.3, 0.5),
    moe = c(0.05, 0.04),
    resp_rate = c(0.6, NA)
  )
  result <- n_multi(indicators, resp_rate = 0.8)
  expect_equal(result$indicators$resp_rate, c(0.6, 0.8))
  expect_equal(
    result$detail$.n,
    c(
      n_prop(0.3, moe = 0.05, resp_rate = 0.6)$n,
      n_prop(0.5, moe = 0.04, resp_rate = 0.8)$n
    )
  )
})

test_that("prec_multi accepts response rate directly and through a plan", {
  indicators <- data.frame(p = c(0.3, 0.5), n = c(500, 500))
  direct <- prec_multi(indicators, resp_rate = 0.8)
  column <- prec_multi(transform(indicators, resp_rate = 0.8))
  planned <- prec_multi(indicators, plan = svyplan(resp_rate = 0.8))
  expect_equal(direct$detail, column$detail)
  expect_equal(planned$detail, column$detail)

  sized <- n_multi(
    transform(indicators, n = NULL, moe = c(0.05, 0.04)),
    plan = svyplan(resp_rate = 0.8)
  )
  expect_equal(sized$indicators$resp_rate, c(0.8, 0.8))
})

test_that("multi round trips can override the response rate", {
  sized <- n_multi(data.frame(p = 0.3, moe = 0.05), resp_rate = 0.8)
  precision <- prec_multi(sized, resp_rate = 0.5)
  expect_equal(precision$params$indicators$resp_rate, 0.5)
  expect_gt(precision$detail$.moe, 0.05)

  resized <- n_multi(precision, resp_rate = 0.9)
  expect_equal(resized$indicators$resp_rate, 0.9)
})

test_that("multi-cluster scalar stage rates match indicator columns", {
  indicators <- data.frame(p = 0.3, cv = 0.1, icc_psu = 0.05)
  direct <- n_cluster(
    indicators = indicators,
    stage_cost = c(500, 50),
    resp_rate_psu = 0.9,
    resp_rate = 0.8
  )
  column <- n_cluster(
    indicators = transform(indicators, resp_rate_psu = 0.9, resp_rate = 0.8),
    stage_cost = c(500, 50)
  )
  expect_equal(direct$n, column$n)
  expect_equal(direct$total_n, column$total_n)

  achieved <- prec_cluster(
    indicators = data.frame(p = 0.3, n = 40, n_per_psu = 10, icc_psu = 0.05),
    resp_rate_psu = 0.9,
    resp_rate = 0.8
  )
  achieved_column <- prec_cluster(indicators = data.frame(
    p = 0.3, n = 40, n_per_psu = 10, icc_psu = 0.05,
    resp_rate_psu = 0.9, resp_rate = 0.8
  ))
  expect_equal(achieved$detail, achieved_column$detail)
})

test_that("multi-cluster validates every stage response rate", {
  three_stage <- data.frame(
    p = 0.3, cv = 0.1, icc_psu = 0.05, icc_ssu = 0.1
  )
  for (rate in c("resp_rate_psu", "resp_rate_ssu", "resp_rate")) {
    bad <- three_stage
    bad[[rate]] <- 0
    expect_error(
      n_cluster(indicators = bad, stage_cost = c(500, 100, 50)),
      paste0("'", rate, "' values must be in")
    )
  }

  expect_error(
    prec_cluster(indicators = data.frame(
      p = 0.3, n = 20, n_per_psu = 10, icc_psu = 0.05,
      resp_rate_psu = 0
    )),
    "'resp_rate_psu' values must be in"
  )
  expect_error(
    n_cluster(
      indicators = data.frame(p = 0.3, cv = 0.1, icc_psu = 0.05),
      stage_cost = c(500, 50),
      resp_rate_ssu = 0.9
    ),
    "not applicable for 2-stage"
  )
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

## Stage-specific response rates: the cluster family names its stage

test_that("the cluster family's two rates act at different stages", {
  base <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  psu <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                   resp_rate_psu = 0.8)
  unit <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                    resp_rate = 0.8)

  # Losing whole PSUs is a pure sample-size loss: the take is unchanged and
  # stage 1 is inflated to cover it.
  expect_equal(unname(psu$n[["n_per_psu"]]), unname(base$n[["n_per_psu"]]),
               tolerance = 1e-9)
  expect_equal(unname(psu$n[["n_psu"]]), unname(base$n[["n_psu"]]) / 0.8,
               tolerance = 1e-9)

  # Losing ultimate units shrinks the realized cluster and so moves the
  # cost-optimal take itself; the two are not interchangeable.
  expect_false(isTRUE(all.equal(unname(unit$n[["n_per_psu"]]),
                                unname(base$n[["n_per_psu"]]))))
  expect_false(isTRUE(all.equal(unname(unit$n[["n_psu"]]),
                                unname(psu$n[["n_psu"]]))))
})

test_that("stage columns are refused where the design has no such stage", {
  # A single-stage design has no PSU or SSU stage to lose units at.
  tg2 <- data.frame(name = c("a", "b"), p = c(0.30, 0.10),
                    moe = c(0.05, 0.05), resp_rate_psu = c(0.5, 0.9))
  expect_error(n_multi(tg2), "stages this design does not have")
  expect_error(
    prec_multi(transform(tg2, moe = NULL, n = 500)),
    "stages this design does not have"
  )

  # A 2-stage cluster design has no SSU stage; those units are the ultimate
  # ones, so their nonresponse is 'resp_rate'.
  tg3 <- data.frame(name = c("a", "b"), p = c(0.30, 0.10),
                    cv = c(0.10, 0.15), icc_psu = c(0.02, 0.05),
                    resp_rate_ssu = c(0.9, 0.9))
  expect_error(n_cluster(indicators = tg3, stage_cost = c(500, 50)),
               "not applicable for 2-stage")
  expect_error(
    prec_cluster(indicators = transform(tg3, cv = NULL, n = 50, n_per_psu = 10)),
    "not applicable for 2-stage"
  )
})

test_that("a cluster frame's two rates act at different stages", {
  base <- data.frame(name = "a", p = 0.30, cv = 0.10, icc_psu = 0.02)
  psu <- transform(base, resp_rate_psu = 0.8)
  unit <- transform(base, resp_rate = 0.8)

  r0 <- n_cluster(indicators = base, stage_cost = c(500, 50))
  r1 <- n_cluster(indicators = psu, stage_cost = c(500, 50))
  r2 <- n_cluster(indicators = unit, stage_cost = c(500, 50))

  expect_equal(unname(r1$n[["n_per_psu"]]), unname(r0$n[["n_per_psu"]]),
               tolerance = 1e-9)
  expect_false(isTRUE(all.equal(unname(r2$n[["n_per_psu"]]),
                                unname(r0$n[["n_per_psu"]]))))
})

test_that("resp_rate_psu validation names the argument it rejects", {
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
              resp_rate_psu = 0),
    "'resp_rate_psu' must be a number"
  )
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
              resp_rate_psu = 1.5),
    "'resp_rate_psu' must be a number"
  )
})

test_that("a plan carries resp_rate_psu to the cluster family", {
  plan <- svyplan(stage_cost = c(500, 50), icc = 0.05, resp_rate_psu = 0.85)
  expect_equal(plan$defaults$resp_rate_psu, 0.85)
  expect_equal(
    n_cluster(cv = 0.05, plan = plan)$n,
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
              resp_rate_psu = 0.85)$n
  )
  expect_error(svyplan(resp_rate_psu = 0), "'resp_rate_psu' must be a number")
})

test_that("the three stage rates enter the CV where decision 7 derives them", {
  V <- 1.5; k1 <- 1; k2 <- 1; d1 <- 0.03; d2 <- 0.10
  r1 <- 0.95; r2 <- 0.90; r3 <- 0.85
  a <- 200; m <- 25; q <- 2

  got <- prec_cluster(n = c(a, m, q), icc = c(d1, d2), unit_relvar = V,
                      var_ratio = c(k1, k2), resp_rate_psu = r1,
                      resp_rate_ssu = r2, resp_rate = r3)$cv
  mr <- m * r2; qr <- q * r3
  manual <- sqrt(V * k1 * d1 / (a * r1) +
                   V * k2 * (1 + d2 * (qr - 1)) / (a * r1 * mr * qr))
  expect_equal(got, manual, tolerance = 1e-12)

  # Lumping the two later rates into the PSU rate is the error the split
  # exists to prevent: it is conservative but wrong, by about 5 percent here.
  lumped <- prec_cluster(n = c(a, m, q), icc = c(d1, d2), unit_relvar = V,
                         var_ratio = c(k1, k2), resp_rate_psu = r1 * r2 * r3)$cv
  expect_gt(lumped, got)
  expect_gt(lumped / got, 1.03)
})

test_that("the cost-optimal take follows sqrt(C1(1-icc)/(C2 icc r))", {
  C1 <- 500; C2 <- 50; icc <- 0.05
  for (r in c(1, 0.9, 0.8, 0.6)) {
    got <- n_cluster(stage_cost = c(C1, C2), icc = icc, cv = 0.05,
                     resp_rate = r)$n[["n_per_psu"]]
    expect_equal(unname(got), sqrt(C1 * (1 - icc) / (C2 * icc * r)),
                 tolerance = 1e-9)
  }

  # resp_rate_psu is a pure 1/n_psu factor and must not move the take.
  base <- n_cluster(stage_cost = c(C1, C2), icc = icc, cv = 0.05)
  psu <- n_cluster(stage_cost = c(C1, C2), icc = icc, cv = 0.05,
                   resp_rate_psu = 0.7)
  expect_equal(unname(psu$n[["n_per_psu"]]), unname(base$n[["n_per_psu"]]),
               tolerance = 1e-12)
  expect_equal(unname(psu$n[["n_psu"]]), unname(base$n[["n_psu"]]) / 0.7,
               tolerance = 1e-9)
})

test_that("cost is charged on the gross take and n stays gross", {
  C1 <- 500; C2 <- 50
  res <- n_cluster(stage_cost = c(C1, C2), icc = 0.05, cv = 0.05,
                   resp_rate = 0.8)
  a <- res$n[["n_psu"]]; m <- res$n[["n_per_psu"]]
  expect_equal(res$cost, a * (C1 + C2 * m), tolerance = 1e-8)
  expect_equal(res$total_n, a * m, tolerance = 1e-9)

  # and the operational design pays for whole issued units
  op <- res$operational
  expect_equal(op$cost, op$n[["n_psu"]] * (C1 + C2 * op$n[["n_per_psu"]]),
               tolerance = 1e-8)
})

test_that("n_alloc rebuilds its inflation factor from the responding take", {
  fr <- data.frame(stratum = c("A", "B"), N = c(50000, 60000), sd = c(10, 12),
                   icc_psu = c(0.05, 0.05),
                   cost_psu = c(500, 500), cost_ssu = c(10, 10))
  base <- n_alloc(fr, n = 2000)$detail$n_per_psu
  unit <- n_alloc(transform(fr, resp_rate = 0.8), n = 2000)$detail$n_per_psu
  psu <- n_alloc(transform(fr, resp_rate_psu = 0.8), n = 2000)$detail$n_per_psu

  expect_equal(unname(unit[1]), sqrt(500 / 10 * (1 - 0.05) / (0.05 * 0.8)),
               tolerance = 1e-9)
  expect_equal(psu, base, tolerance = 1e-12)
})

test_that("the stratum SD is inflated by the responding take, not the issued one", {
  # Takes are fixed, so the take is not the channel here. The strata differ
  # only in response, so their realized cluster sizes differ and the
  # clustering penalty must differ with them; using the gross take would
  # inflate both alike and misallocate between them.
  fr <- data.frame(stratum = c("A", "B"), N = c(50000, 60000),
                   sd = c(10, 10), icc_psu = c(0.05, 0.05),
                   n_per_psu = c(20, 20), resp_rate = c(0.5, 1.0))
  res <- n_alloc(fr, n = 2000)

  # Stratum A realizes clusters of 10 and B of 20, so A suffers less
  # clustering per issued unit and takes the larger share.
  expect_gt(res$detail$n[1], res$detail$n[2])

  # The decisive check: the cluster frame must match an element-level frame
  # whose SDs are inflated by hand at the *responding* take. Using the issued
  # take instead inflates both strata alike and misallocates between them.
  rates <- c(0.5, 1.0)
  equivalent <- n_alloc(
    data.frame(stratum = c("A", "B"), N = c(50000, 60000),
               sd = 10 * sqrt(1 + 0.05 * (20 * rates - 1)),
               resp_rate = rates),
    n = 2000
  )
  expect_equal(res$detail$n, equivalent$detail$n, tolerance = 1e-10)
  expect_equal(res$se, equivalent$se, tolerance = 1e-12)
})

test_that("an expected take below one unit is refused as unsupported", {
  err <- tryCatch(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
              n_per_psu = 1, resp_rate = 0.4),
    error = conditionMessage
  )
  expect_match(err, "expected-take planning approximation is unsupported")
  # It says the approximation is unsupported, not that the design is invalid.
  expect_false(grepl("invalid|impossible", err))
})

test_that("effective_n nets every stage a cluster design loses units at", {
  plan <- n_cluster(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05),
                    cv = 0.05, resp_rate_psu = 0.95, resp_rate_ssu = 0.9,
                    resp_rate = 0.85)
  expect_equal(
    as.double(effective_n(plan)),
    plan$total_n * (0.95 * 0.9 * 0.85) / as.double(design_effect(plan)),
    tolerance = 1e-8
  )
})

test_that("a plan and predict() carry every stage rate", {
  plan <- svyplan(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05),
                  resp_rate_psu = 0.95, resp_rate_ssu = 0.9, resp_rate = 0.85)
  expect_equal(
    n_cluster(cv = 0.05, plan = plan)$n,
    n_cluster(stage_cost = c(500, 100, 50), icc = c(0.02, 0.05), cv = 0.05,
              resp_rate_psu = 0.95, resp_rate_ssu = 0.9, resp_rate = 0.85)$n
  )

  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  grid <- predict(x, data.frame(resp_rate = c(1, 0.8, 0.6)))
  expect_equal(nrow(grid), 3L)
  expect_true(all(diff(grid$n_per_psu) > 0))
})
