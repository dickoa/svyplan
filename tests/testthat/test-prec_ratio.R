## L_R = 1.21 + 0.36 - 2 * 0.7 * 1.1 * 0.6 = 0.646 exactly.
RELVAR <- 0.646

test_that("prec_ratio reports the closed-form precision", {
  res <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                    component_cor = 0.7)
  expected_se <- 2 * sqrt(RELVAR / 400)

  expect_equal(res$se, expected_se, tolerance = 1e-10)
  expect_equal(res$moe, qnorm(0.975) * expected_se, tolerance = 1e-10)
  expect_equal(res$cv, expected_se / 2, tolerance = 1e-10)
  expect_equal(res$rmoe, res$moe / 2, tolerance = 1e-10)
  expect_s3_class(res, "svyplan_prec")
  expect_equal(res$type, "ratio")
  expect_equal(res$method, "linearization")
})

test_that("prec_ratio carries deff, response, and the FPC", {
  res <- prec_ratio(r = 2, n = 500, cv_num = 1.1, cv_den = 0.6,
                    component_cor = 0.7, N = 5000, deff = 1.6,
                    resp_rate = 0.8)
  n_net <- 500 * 0.8
  n_eff <- n_net / 1.6
  fpc <- 1 - n_net / 5000

  expect_equal(res$se, 2 * sqrt(RELVAR * fpc / n_eff), tolerance = 1e-10)
})

test_that("prec_ratio uses a t quantile when df is supplied", {
  normal <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7)
  student <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                        component_cor = 0.7, df = 9)

  expect_equal(normal$se, student$se)
  expect_gt(student$moe, normal$moe)
  expect_equal(student$moe, qt(0.975, 9) * student$se, tolerance = 1e-10)
})

test_that("a census leaves no sampling error", {
  expect_warning(
    res <- prec_ratio(r = 2, n = 1000, cv_num = 1.1, cv_den = 0.6,
                      component_cor = 0.7, N = 1000),
    "census"
  )
  expect_equal(res$se, 0)
  expect_equal(res$moe, 0)
  expect_equal(res$cv, 0)
})

test_that("prec_ratio has no solved field", {
  res <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                    component_cor = 0.7)
  expect_null(res$solved)
})

test_that("prec_ratio rejects a sample larger than the frame", {
  expect_error(
    prec_ratio(r = 2, n = 900, cv_num = 1.1, cv_den = 0.6,
               component_cor = 0.7, N = 800),
    "exceeds the population"
  )
})

test_that("prec_ratio validates its moments and its n", {
  expect_error(
    prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6),
    "needs 'r'"
  )
  expect_error(
    prec_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7),
    "'n' must be a numeric scalar"
  )
  expect_error(
    prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
               component_cor = 0.7, extra = 1),
    "unused"
  )
})

test_that("prec_ratio equals prec_mean on the equivalent variance", {
  ratio <- prec_ratio(r = 2, n = 500, cv_num = 1.1, cv_den = 0.6,
                      component_cor = 0.7, N = 9000, deff = 1.2,
                      resp_rate = 0.85)
  mean_route <- prec_mean(var = 2^2 * RELVAR, mu = 2, n = 500, N = 9000,
                          deff = 1.2, resp_rate = 0.85)

  expect_equal(ratio$se, mean_route$se, tolerance = 1e-12)
  expect_equal(ratio$cv, mean_route$cv, tolerance = 1e-12)
})

test_that("the round trip is exact in every target mode", {
  for (target in list(list(cv = 0.05), list(moe = 0.2), list(rmoe = 0.1))) {
    args <- c(list(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7),
              target)
    size <- do.call(n_ratio, args)
    prec <- prec_ratio(size)
    back <- n_ratio(prec)

    expect_equal(prec$cv, size$cv, tolerance = 1e-10,
                 info = names(target))
    expect_equal(prec$moe, size$moe, tolerance = 1e-10,
                 info = names(target))
    expect_equal(back$n, size$n, tolerance = 1e-10,
                 info = names(target))
  }
})

test_that("the round trip preserves every design parameter", {
  size <- n_ratio(r = -3, cv_num = 0.9, cv_den = 0.4, component_cor = -0.2,
                  rmoe = 0.08, alpha = 0.10, N = 7000, deff = 1.35,
                  resp_rate = 0.72, df = 18)
  prec <- prec_ratio(size)
  back <- n_ratio(prec)

  for (nm in c("r", "cv_num", "cv_den", "component_cor", "unit_relvar",
               "alpha", "N", "deff", "resp_rate", "df")) {
    expect_equal(back$params[[nm]], size$params[[nm]], info = nm)
  }
  expect_equal(back$n, size$n, tolerance = 1e-10)
})

test_that("a round trip can switch the target through dots", {
  size <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                  cv = 0.05)
  prec <- prec_ratio(size)
  switched <- n_ratio(prec, rmoe = 0.2)

  expect_equal(switched$params$rmoe, 0.2)
  expect_null(switched$params$cv)
  expect_equal(switched$rmoe, 0.2, tolerance = 1e-8)
})

test_that("a round trip can override a design parameter", {
  size <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                  cv = 0.05)
  prec <- prec_ratio(size, deff = 2)

  expect_equal(prec$params$deff, 2)
  expect_equal(prec$se, 2 * sqrt(RELVAR * 2 / size$n), tolerance = 1e-10)
})

test_that("round-trip methods refuse the wrong estimand", {
  mean_size <- n_mean(var = 100, mu = 50, cv = 0.05)
  mean_prec <- prec_mean(var = 100, mu = 50, n = 400)

  expect_error(prec_ratio(mean_size), "requires a svyplan_n of type 'ratio'")
  expect_error(n_ratio(mean_prec), "requires a svyplan_prec of type 'ratio'")
})

test_that("prec_ratio takes design defaults from a plan, named and piped", {
  plan <- svyplan(alpha = 0.10, deff = 1.5)
  direct <- prec_ratio(r = 2, n = 500, cv_num = 1.1, cv_den = 0.6,
                       component_cor = 0.7, alpha = 0.10, deff = 1.5)

  named <- prec_ratio(r = 2, n = 500, cv_num = 1.1, cv_den = 0.6,
                      component_cor = 0.7, plan = plan)
  piped <- plan |> prec_ratio(r = 2, n = 500, cv_num = 1.1, cv_den = 0.6,
                              component_cor = 0.7)

  expect_equal(named$se, direct$se)
  expect_equal(named$moe, direct$moe)
  expect_equal(piped$moe, direct$moe)
})

test_that("every prec_ratio formal is reachable through a plan", {
  plan <- svyplan(deff = 1.2)
  for (drop in c("r", "n", "cv_num", "cv_den", "component_cor")) {
    args <- list(r = 2, n = 500, cv_num = 1.1, cv_den = 0.6,
                 component_cor = 0.7, plan = plan)
    args[[drop]] <- NULL
    err <- tryCatch(do.call(prec_ratio, args),
                    error = function(e) conditionMessage(e))
    expect_false(grepl("is missing, with no default", err),
                 info = paste("dropping", drop))
  }
})

test_that("confint on a precision result is centred on r", {
  res <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                    component_cor = 0.7)
  ci <- confint(res)

  expect_equal(mean(as.numeric(ci)), 2, tolerance = 1e-10)
  expect_equal(as.numeric(ci)[2L] - 2, res$moe, tolerance = 1e-10)
})

test_that("predict on a precision result varies n and the moments", {
  res <- prec_ratio(r = 2, n = 400, cv_num = 1.1, cv_den = 0.6,
                    component_cor = 0.7)

  by_n <- predict(res, data.frame(n = c(400, 1600)))
  expect_equal(by_n$se[1L], res$se, tolerance = 1e-10)
  expect_equal(by_n$se[1L] / by_n$se[2L], 2, tolerance = 1e-8)

  by_cor <- predict(res, data.frame(component_cor = c(0.3, 0.7, 0.95)))
  expect_true(all(diff(by_cor$se) < 0))
})

test_that("existing prec schemas are untouched", {
  mean_prec <- prec_mean(var = 100, mu = 50, n = 400)
  expect_false("unit_relvar" %in% names(mean_prec$params))
  expect_false("component_cor" %in% names(mean_prec$params))
  expect_equal(mean_prec$type, "mean")
})
