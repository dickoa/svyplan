## Expected sizes come from the documented formula, never from the result.
## r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7 gives
## L_R = 1.21 + 0.36 - 2 * 0.7 * 1.1 * 0.6 = 0.646 exactly.
RELVAR <- 0.646

test_that("n_ratio CV mode matches the sizing equation", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)

  expect_equal(res$n, RELVAR / 0.05^2, tolerance = 1e-10)
  expect_equal(res$params$unit_relvar, RELVAR)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$type, "ratio")
  expect_equal(res$method, "linearization")
})

test_that("n_ratio MOE mode converts on the ratio scale", {
  q <- qnorm(0.975)
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 moe = 0.2)

  expect_equal(res$n, RELVAR / (0.2 / (q * 2))^2, tolerance = 1e-10)
  expect_equal(res$moe, 0.2, tolerance = 1e-8)
})

test_that("n_ratio relative MOE mode is moe / abs(r)", {
  by_rmoe <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                     rmoe = 0.1)
  by_moe <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                    moe = 0.2)

  expect_equal(by_rmoe$n, by_moe$n, tolerance = 1e-10)
  expect_equal(by_rmoe$rmoe, 0.1, tolerance = 1e-8)
})

test_that("n_ratio applies the FPC with deff inside it", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05, N = 3000, deff = 1.5)

  expect_equal(res$n, 1.5 * RELVAR / (0.05^2 + 1.5 * RELVAR / 3000),
               tolerance = 1e-10)
})

test_that("n_ratio returns a gross size and reports the net one", {
  net <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)
  gross <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                   cv = 0.05, resp_rate = 0.75)

  expect_equal(gross$n, net$n / 0.75, tolerance = 1e-10)
  expect_equal(gross$cv, 0.05, tolerance = 1e-10)
})

test_that("n_ratio switches to a t quantile with df", {
  q <- qt(0.975, 10)
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 moe = 0.2, df = 10)

  expect_equal(res$n, RELVAR / (0.2 / (q * 2))^2, tolerance = 1e-10)
})

test_that("a negative ratio sizes like its mirror image", {
  positive <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                      cv = 0.05)
  negative <- n_ratio(r = -2, cv_num = 1.1, cv_den = 0.6,
                      component_cor = -0.7, cv = 0.05)

  expect_equal(negative$n, positive$n, tolerance = 1e-10)
  expect_equal(negative$params$unit_relvar, RELVAR)
})

test_that("n_ratio equals n_mean on the equivalent variance", {
  ratio <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                   cv = 0.05, N = 8000, deff = 1.3, resp_rate = 0.8)
  mean_route <- n_mean(var = 2^2 * RELVAR, mu = 2, cv = 0.05, N = 8000,
                       deff = 1.3, resp_rate = 0.8)

  expect_equal(ratio$n, mean_route$n, tolerance = 1e-10)
  expect_equal(ratio$se, mean_route$se, tolerance = 1e-10)
})

test_that("n_ratio rejects an incomplete or invalid moment set", {
  expect_error(n_ratio(r = 2, cv_num = 1.1, cv = 0.05), "needs 'r'")
  expect_error(
    n_ratio(r = 0, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7, cv = 0.05),
    "must not be zero"
  )
  expect_error(
    n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 1.2, cv = 0.05),
    "must be a number in \\[-1, 1\\]"
  )
})

test_that("n_ratio requires exactly one precision target", {
  expect_error(
    n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7),
    "exactly one of 'moe', 'cv', or 'rmoe'"
  )
  expect_error(
    n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
            cv = 0.05, moe = 0.2),
    "exactly one of 'moe', 'cv', or 'rmoe'"
  )
})

test_that("the component moments are named-only", {
  expect_error(n_ratio(2, 1.1, 0.6, 0.7, cv = 0.05), "unused")
})

test_that("n_ratio rejects an unused argument", {
  expect_error(
    n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
            cv = 0.05, mu = 3),
    "unused"
  )
})

test_that("n_ratio refuses a degenerate ratio and warns near it", {
  expect_error(
    n_ratio(r = 2, cv_num = 0.8, cv_den = 0.8, component_cor = 1, cv = 0.05),
    "no sampling variance"
  )
  expect_warning(
    n_ratio(r = 2, cv_num = 1, cv_den = 1, component_cor = 0.9996, cv = 0.05),
    "known to that precision"
  )
})

test_that("n_ratio refuses a target the frame cannot deliver", {
  expect_error(
    n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
            cv = 0.001, N = 400, resp_rate = 0.5),
    "unattainable"
  )
})

test_that("under full response the FPC keeps any CV target attainable", {
  # n_net stays below N for every positive target, so only nonresponse can
  # push the gross size past the frame.
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 1e-8, N = 400)
  expect_lt(res$n, 400)
})

test_that("n_ratio takes design defaults from a plan, named and piped", {
  plan <- svyplan(alpha = 0.10, deff = 1.5, resp_rate = 0.8)
  direct <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                    cv = 0.05, alpha = 0.10, deff = 1.5, resp_rate = 0.8)

  named <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                   cv = 0.05, plan = plan)
  piped <- plan |> n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6,
                           component_cor = 0.7, cv = 0.05)

  expect_equal(named$n, direct$n)
  expect_equal(piped$n, direct$n)
})

test_that("every formal is reachable through a plan without an R-level error", {
  # A defaultless formal raises R's own "argument is missing" from mget()
  # inside .merge_plan_args(), before any validator runs.
  plan <- svyplan(deff = 1.2)
  for (drop in c("r", "cv_num", "cv_den", "component_cor")) {
    args <- list(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05, plan = plan)
    args[[drop]] <- NULL
    err <- tryCatch(do.call(n_ratio, args), error = function(e) conditionMessage(e))
    expect_false(grepl("is missing, with no default", err),
                 info = paste("dropping", drop))
  }
})

test_that("the result carries the moments and the derived coefficient", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)
  p <- res$params

  expect_equal(p$r, 2)
  expect_equal(p$cv_num, 1.1)
  expect_equal(p$cv_den, 0.6)
  expect_equal(p$component_cor, 0.7)
  expect_equal(p$unit_relvar, RELVAR)
  expect_equal(p$cv, 0.05)
  expect_null(p$moe)
  expect_null(p$rmoe)
  # The equivalent variance is a step, not an observed quantity.
  expect_null(p$var)
})

test_that("printing reports r and the derived coefficient, not the moments", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05, deff = 1.4, resp_rate = 0.9)
  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  line <- capture.output(print(res))[2L]

  expect_match(line, "r = 2")
  expect_match(line, "unit_relvar = 0.646")
  expect_false(grepl("cv_num", line))
  expect_false(grepl("component_cor", line))
})

test_that("a ratio print stays close to the mean's at the same size", {
  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05, deff = 1.4, resp_rate = 0.9)
  reference <- n_mean(var = 100, mu = 50, cv = 0.05, deff = 1.4,
                      resp_rate = 0.9)
  reference$n <- res$n

  ratio_line <- nchar(capture.output(print(res))[2L])
  mean_line <- nchar(capture.output(print(reference))[2L])

  expect_lte(ratio_line - mean_line, 14L)
})

test_that("no other estimand's print block gains a ratio field", {
  # 'p$r' partial-matches 'resp_rate', so an unguarded read prints
  # 'r = 0.9' on means and proportions.
  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  mean_line <- capture.output(
    print(n_mean(var = 100, mu = 50, cv = 0.05, resp_rate = 0.9))
  )[2L]
  prop_line <- capture.output(
    print(n_prop(p = 0.3, moe = 0.05, resp_rate = 0.9))
  )[2L]

  # Word-bounded: "var = 100.00" contains "r = ".
  expect_false(grepl("\\br = ", mean_line))
  expect_false(grepl("\\br = ", prop_line))
  expect_false(grepl("unit_relvar", mean_line))
  expect_false(grepl("unit_relvar", prop_line))
})

test_that("existing result schemas are untouched by ratio support", {
  mean_res <- n_mean(var = 100, mu = 50, cv = 0.05)
  prop_res <- n_prop(p = 0.3, moe = 0.05)

  expect_named(
    mean_res$params,
    c("var", "alpha", "N", "deff", "resp_rate", "df", "mu", "cv")
  )
  expect_false("unit_relvar" %in% names(prop_res$params))
  expect_false("r" %in% names(mean_res$params))
})

test_that("confint is symmetric and centred on r", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)
  ci <- confint(res)

  expect_equal(mean(as.numeric(ci)), 2, tolerance = 1e-10)
  expect_equal(as.numeric(ci)[2L] - 2, res$moe, tolerance = 1e-10)
  expect_equal(2 - as.numeric(ci)[1L], res$moe, tolerance = 1e-10)
})

test_that("confint at another level widens the interval", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)
  narrow <- confint(res, level = 0.90)
  wide <- confint(res, level = 0.99)

  expect_lt(diff(as.numeric(narrow)), diff(as.numeric(wide)))
})

test_that("coercion works without a ratio-specific branch", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)

  expect_equal(as.double(res), res$n)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
  expect_s3_class(as.data.frame(res), "data.frame")
})

test_that("predict varies the moments and recomputes the coefficient", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)
  grid <- predict(res, data.frame(component_cor = c(0.5, 0.7, 0.9)))

  expect_equal(nrow(grid), 3L)
  expect_equal(grid$n[2L], res$n, tolerance = 1e-10)
  # A stale unit_relvar would leave the column flat.
  expect_true(all(diff(grid$n) < 0))
})

test_that("predict switches the target and honours df", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)

  switched <- predict(res, data.frame(moe = c(0.2, 0.4)))
  expect_equal(switched$moe, c(0.2, 0.4), tolerance = 1e-8)
})

test_that("df moves the size only when the target is a margin of error", {
  # The quantile enters only where a moe is converted to a relative standard
  # error. In CV mode df still moves the reported moe, an interval width.
  by_cv <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                   cv = 0.05)
  cv_grid <- predict(by_cv, data.frame(df = c(5, 1000)))
  expect_equal(cv_grid$n[1L], cv_grid$n[2L])
  expect_gt(cv_grid$moe[1L], cv_grid$moe[2L])

  by_moe <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                    moe = 0.2)
  moe_grid <- predict(by_moe, data.frame(df = c(5, 1000)))
  expect_gt(moe_grid$n[1L], moe_grid$n[2L])
})

test_that("predict rejects a column the estimand does not have", {
  res <- n_ratio(r = 2, cv_num = 1.1, cv_den = 0.6, component_cor = 0.7,
                 cv = 0.05)
  expect_error(predict(res, data.frame(var = c(1, 2))), "var")
})
