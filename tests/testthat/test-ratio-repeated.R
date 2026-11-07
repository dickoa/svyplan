## A ratio of totals in a panel, a change and a pooled average

# An unemployment rate: unemployed over the labor force, both totals over
# all persons, 4 percent unemployed and 65 percent in the labor force. The
# unit coefficients and their correlation follow from the two shares.
pu <- 0.04
pl <- 0.65
ur <- list(
  r = pu / pl,
  cv_num = sqrt((1 - pu) / pu),
  cv_den = sqrt((1 - pl) / pl),
  component_cor = sqrt(pu * (1 - pl) / (pl * (1 - pu)))
)
ur_L <- ur$cv_num^2 + ur$cv_den^2 -
  2 * ur$component_cor * ur$cv_num * ur$cv_den
ret <- c(0.9, 0.92, 0.94, 0.95)

rate_n <- function(...) {
  do.call(n_ratio, c(ur, list(...)))
}

## R1. A ratio target in a panel

test_that("a panel's head precision is prec_ratio() at the pooled respondents", {
  tgt <- rate_n(cv = 0.05, deff = 1.5, resp_rate = 0.6)
  rot <- n_panel(tgt, retention = ret, resp_rate = 0.6, design = "rotating")
  direct <- do.call(prec_ratio, c(ur, list(n = rot$n_resp, deff = 1.5)))
  expect_equal(rot$cv, 0.05, tolerance = 1e-10)
  expect_equal(rot$se, direct$se, tolerance = 1e-10)
  expect_identical(rot$type, "ratio")

  fx <- n_panel(tgt, retention = ret, resp_rate = 0.6)
  wave5 <- do.call(prec_ratio, c(ur, list(n = fx$n_resp, deff = 1.5)))
  expect_equal(fx$se, wave5$se, tolerance = 1e-10)
  per_wave <- vapply(fx$waves$n_resp, function(m) {
    do.call(prec_ratio, c(ur, list(n = m, deff = 1.5)))$cv
  }, numeric(1))
  expect_equal(fx$waves$cv, per_wave, tolerance = 1e-10)
})

test_that("a subset rate matches its proportion among the labor force", {
  # The rate is a proportion among labor force members, so the two plans
  # agree on entrants and on every wave's precision. Only the count of
  # respondents differs, by the labor force share, since the ratio counts
  # persons and the proportion counts members.
  tr <- rate_n(cv = 0.05, deff = 1.5, resp_rate = 0.6)
  tp <- n_prop(p = pu / pl, cv = 0.05, deff = 1.5, resp_rate = 0.6,
               eligibility_rate = pl)
  expect_equal(tr$n, tp$n, tolerance = 1e-10)
  rot <- n_panel(tr, retention = ret, resp_rate = 0.6, design = "rotating")
  alt <- suppressWarnings(n_panel(tp, retention = ret, resp_rate = 0.6,
                                  eligibility_rate = pl, design = "rotating"))
  expect_equal(rot$n_entrants, alt$n_entrants, tolerance = 1e-10)
  expect_equal(rot$waves$cv, alt$waves$cv, tolerance = 1e-10)
  expect_equal(alt$n_resp / rot$n_resp, pl, tolerance = 1e-10)
})

test_that("a ratio panel round-trips and takes every target form", {
  cont <- n_ratio(r = 1200, cv_num = 0.9, cv_den = 0.45, component_cor = 0.6,
                  cv = 0.02, deff = 2)
  fx <- n_panel(cont, retention = c(0.9, 0.9), resp_rate = 0.7,
                assurance = 0.9)
  back <- prec_panel(fx)
  expect_equal(back$cv, fx$cv, tolerance = 1e-10)
  expect_equal(back$n_issued, fx$n_issued, tolerance = 1e-10)

  t_moe <- rate_n(moe = 0.005, df = 30)
  expect_equal(n_panel(t_moe, retention = 0.9, resp_rate = 0.6)$moe, 0.005,
               tolerance = 1e-10)
  t_rmoe <- rate_n(rmoe = 0.1)
  expect_equal(n_panel(t_rmoe, retention = 0.9, resp_rate = 0.6)$rmoe, 0.1,
               tolerance = 1e-10)
  t_prec <- do.call(prec_ratio, c(ur, list(n = 20000)))
  expect_equal(n_panel(t_prec, retention = 0.9, resp_rate = 0.6)$cv,
               t_prec$cv, tolerance = 1e-10)
})

test_that("a ratio panel keeps the target's population in scope", {
  tgt <- rate_n(cv = 0.05, N = 60000, eligibility_rate = 0.8)
  fx <- suppressWarnings(n_panel(tgt, retention = 0.9, resp_rate = 0.6))
  direct <- do.call(prec_ratio, c(ur, list(n = fx$n_resp, N = 60000 * 0.8)))
  expect_equal(fx$cv, direct$cv, tolerance = 1e-10)
})

test_that("a ratio panel's schedule records the estimand", {
  imm <- n_panel(rate_n(cv = 0.05), retention = rep(0.9, 3), resp_rate = 0.6,
                 design = "rotating", rotation = "2-2-2", start = "immediate")
  sc <- design_schedule(
    imm, horizon = 8, horizon_policy = "continuing",
    refreshment = "entrant_register",
    frame_vintage = c(startup = "2026Q1", period_2 = "2026Q2")
  )
  expect_identical(sc$target_estimand,
                   list(type = "ratio", method = "linearization"))
})

test_that("the panel refusals name the ratio functions", {
  expect_error(n_panel(list(), retention = 0.9), "n_ratio\\(\\)")
  ch <- n_change(var = 100, moe = 2)
  expect_error(n_panel(ch, retention = 0.9),
               "a mean, a proportion or a ratio, not 'change'")
})

## R2. A change in a ratio

manual_change <- function(r, ...) {
  n_change(var = r^2 * ur_L, change = r[2L] - r[1L], ...)
}

test_that("a ratio change is the change engine on the linearized variable", {
  r <- c(ur$r, 1.05 * ur$r)
  for (target in list(list(moe = 0.005), list(cv = 0.4), list(rmoe = 0.8))) {
    args <- c(target, list(overlap = 0.8, overlap_cor = 0.54, deff = 1.4,
                           resp_rate = 0.7, N = 2e6))
    a <- do.call(n_change, c(list(r = r), ur[-1L], args))
    b <- do.call(manual_change, c(list(r = r), args))
    expect_equal(a$n, b$n, tolerance = 1e-10)
    expect_equal(c(a$se, a$moe, a$cv, a$rmoe), c(b$se, b$moe, b$cv, b$rmoe),
                 tolerance = 1e-10)
  }
  pa <- do.call(prec_change, c(list(r = r, n = c(9000, 7000)), ur[-1L],
                               list(overlap = 0.5, overlap_cor = 0.6)))
  pb <- prec_change(var = r^2 * ur_L, change = r[2L] - r[1L],
                    n = c(9000, 7000), overlap = 0.5, overlap_cor = 0.6)
  expect_equal(c(pa$se, pa$moe, pa$cv), c(pb$se, pb$moe, pb$cv),
               tolerance = 1e-10)
})

test_that("each occasion of a ratio change takes its own moments", {
  a <- n_change(r = c(0.06, 0.07), cv_num = c(4.9, 4.6), cv_den = 0.73,
                component_cor = c(0.15, 0.17), moe = 0.005)
  L <- c(4.9^2 + 0.73^2 - 2 * 0.15 * 4.9 * 0.73,
         4.6^2 + 0.73^2 - 2 * 0.17 * 4.6 * 0.73)
  b <- n_change(var = c(0.06, 0.07)^2 * L, change = 0.01, moe = 0.005)
  expect_equal(a$n, b$n, tolerance = 1e-10)
  expect_equal(a$params$unit_relvar, L, tolerance = 1e-12)
  expect_equal(a$params$change, 0.01, tolerance = 1e-12)
})

test_that("a ratio change round-trips with its moments", {
  fit <- do.call(n_change, c(list(r = c(ur$r, 1.1 * ur$r)), ur[-1L],
                             list(moe = 0.005, overlap = 0.8,
                                  overlap_cor = 0.54, ratio = 1.2)))
  back <- prec_change(fit)
  expect_equal(back$moe, 0.005, tolerance = 1e-10)
  expect_equal(back$params[["r"]], fit$params[["r"]])
  expect_equal(back$params$component_cor, ur$component_cor)
  expect_equal(n_change(back)$n, fit$n, tolerance = 1e-10)
  # An override of a moment reaches the rebuild, so the stored one does too
  worse <- prec_change(fit, component_cor = 0)
  expect_gt(worse$se, back$se)
})

test_that("predict() varies a ratio change's moments", {
  fit <- do.call(n_change, c(list(r = c(ur$r, 1.1 * ur$r)), ur[-1L],
                             list(moe = 0.005, overlap = 0.8,
                                  overlap_cor = 0.54)))
  grid <- predict(fit, data.frame(component_cor = c(0, 0.3)))
  direct <- vapply(c(0, 0.3), function(cc) {
    n_change(r = c(ur$r, 1.1 * ur$r), cv_num = ur$cv_num,
             cv_den = ur$cv_den, component_cor = cc, moe = 0.005,
             overlap = 0.8, overlap_cor = 0.54)$n
  }, numeric(1))
  expect_equal(grid$n1, direct, tolerance = 1e-10)
  expect_error(predict(fit, data.frame(change = 0.01)), "change")
  pf <- prec_change(fit)
  pgrid <- predict(pf, data.frame(cv_num = c(4, 5)))
  expect_equal(nrow(pgrid), 2L)
})

test_that("a ratio change prints and centers its interval on r[2] - r[1]", {
  fit <- n_change(r = c(0.06, 0.07), cv_num = 4.9, cv_den = 0.73,
                  component_cor = 0.15, moe = 0.005)
  out <- capture.output(print(fit))
  expect_match(out[1L], "ratio of totals")
  expect_true(any(grepl("r = 0.06 to 0.07", out, fixed = TRUE)))
  ci <- confint(prec_change(fit))
  expect_equal(mean(ci), 0.01, tolerance = 1e-10)
})

test_that("a ratio change is stated by r alone", {
  m <- ur[-1L]
  expect_error(do.call(n_change, c(list(r = c(0.06, 0.07), change = 0.01,
                                        moe = 0.005), m)),
               "r\\[2\\] - r\\[1\\]")
  expect_error(do.call(n_change, c(list(r = c(0.06, 0.07), var = 1,
                                        moe = 0.005), m)),
               "exactly one of 'var' \\(or 'sd'\\), 'p', or 'r'")
  expect_error(n_change(var = 1, cv_num = 2, moe = 0.1), "need 'r'")
  expect_error(do.call(n_change, c(list(r = 0.06, moe = 0.005), m)),
               "two ratios")
  expect_error(n_change(r = c(0.06, 0.07), cv_num = c(1, 2, 3), cv_den = 1,
                        component_cor = 0, moe = 0.005), "one per occasion")
  expect_error(n_change(r = c(0.06, 0.07), cv_den = 1, component_cor = 0,
                        moe = 0.005), "a ratio needs")
  expect_error(do.call(n_change, c(list(r = c(0.06, 0.06), cv = 0.1), m)),
               "r\\[1\\] equals r\\[2\\]")
})

test_that("a ratio change keeps its standard error at an extreme scale", {
  # Forming r^2 * L would overflow at 1e200 and underflow at 1e-200.
  base <- n_change(r = c(1, 1.1), cv_num = 1, cv_den = 0.5,
                   component_cor = 0.3, rmoe = 0.2, overlap = 0.5,
                   overlap_cor = 0.5)
  for (s in c(1e200, 1e-200)) {
    big <- n_change(r = s * c(1, 1.1), cv_num = 1, cv_den = 0.5,
                    component_cor = 0.3, rmoe = 0.2, overlap = 0.5,
                    overlap_cor = 0.5)
    expect_equal(big$n, base$n, tolerance = 1e-10)
    expect_equal(big$se / s, base$se, tolerance = 1e-10)
    m <- n_change(r = s * c(1, 1.1), cv_num = 1, cv_den = 0.5,
                  component_cor = 0.3, moe = s * base$moe, overlap = 0.5,
                  overlap_cor = 0.5)
    expect_equal(m$n, base$n, tolerance = 1e-10)
    p <- prec_change(r = s * c(1, 1.1), cv_num = 1, cv_den = 0.5,
                     component_cor = 0.3, n = base$n, overlap = 0.5,
                     overlap_cor = 0.5)
    expect_equal(p$moe / s, base$moe, tolerance = 1e-10)
  }
})

## R3. A pooled ratio

test_that("a pooled ratio is the pooled engine on the linearized variable", {
  args <- list(occasions = 4, overlap = 0.75, cor_decay = 0.8, deff = 1.3,
               resp_rate = 0.7, N = 1e6)
  for (target in list(list(moe = 0.003), list(cv = 0.03),
                      list(rmoe = 0.05))) {
    a <- do.call(n_pooled, c(ur, target, args))
    b <- do.call(n_pooled, c(list(var = ur$r^2 * ur_L, mu = ur$r), target,
                             args))
    expect_equal(a$n, b$n, tolerance = 1e-10)
    expect_equal(c(a$se, a$moe, a$cv), c(b$se, b$moe, b$cv),
                 tolerance = 1e-10)
  }
  pa <- do.call(prec_pooled, c(ur, list(n = 15000), args))
  pb <- do.call(prec_pooled, c(list(var = ur$r^2 * ur_L, mu = ur$r,
                                    n = 15000), args))
  expect_equal(c(pa$se, pa$cv), c(pb$se, pb$cv), tolerance = 1e-10)
})

test_that("a pooled ratio round-trips, predicts and prints", {
  fit <- do.call(n_pooled, c(ur, list(cv = 0.03, occasions = 4,
                                      overlap = 0.75, cor_decay = 0.8)))
  back <- prec_pooled(fit)
  expect_equal(back$cv, 0.03, tolerance = 1e-10)
  expect_equal(back$params[["r"]], ur$r)
  expect_equal(n_pooled(back, cv = 0.03)$n, fit$n, tolerance = 1e-10)
  grid <- predict(fit, data.frame(r = c(0.05, 0.08)))
  direct <- vapply(c(0.05, 0.08), function(rr) {
    n_pooled(r = rr, cv_num = ur$cv_num, cv_den = ur$cv_den,
             component_cor = ur$component_cor, cv = 0.03, occasions = 4,
             overlap = 0.75, cor_decay = 0.8)$n
  }, numeric(1))
  expect_equal(grid$n, direct, tolerance = 1e-10)
  expect_error(predict(fit, data.frame(mu = 0.05)), "mu")
  pgrid <- predict(back, data.frame(r = c(0.05, 0.08)))
  pdirect <- vapply(c(0.05, 0.08), function(rr) {
    prec_pooled(r = rr, cv_num = ur$cv_num, cv_den = ur$cv_den,
                component_cor = ur$component_cor, n = fit$n, occasions = 4,
                overlap = 0.75, cor_decay = 0.8)$se
  }, numeric(1))
  expect_equal(pgrid$se, pdirect, tolerance = 1e-10)
  expect_equal(mean(confint(fit)), ur$r, tolerance = 1e-10)
  out <- capture.output(print(fit))
  expect_match(out[1L], "ratio of totals")
})

test_that("a pooled ratio is stated by r alone", {
  m <- ur[-1L]
  expect_error(do.call(n_pooled, c(list(r = 0.06, mu = 0.06, cv = 0.03,
                                        occasions = 4), m)),
               "'mu' is 'r'")
  expect_error(do.call(n_pooled, c(list(r = c(0.06, 0.07), cv = 0.03,
                                        occasions = 4), m)),
               "single ratio")
  expect_error(n_pooled(var = 1, component_cor = 0.2, moe = 0.1,
                        occasions = 4), "need 'r'")
})

test_that("variance-scale and proportion-scale results are unchanged", {
  # The ratio scale threads a scale factor through both families, and it is
  # exactly 1 elsewhere. These are the documented example values.
  expect_equal(n_change(var = 100, moe = 2)$n, 192.0729, tolerance = 1e-6)
  expect_equal(n_change(p = c(0.30, 0.36), moe = 0.02, overlap = 0.75,
                        overlap_cor = 0.5)$n,
               n_change(var = c(0.21, 0.36 * 0.64), change = 0.06,
                        moe = 0.02, overlap = 0.75, overlap_cor = 0.5)$n,
               tolerance = 1e-10)
  expect_null(n_change(var = 100, moe = 2)$params[["r"]])
  expect_null(n_pooled(p = 0.3, moe = 0.01, occasions = 4)$params[["r"]])
})
