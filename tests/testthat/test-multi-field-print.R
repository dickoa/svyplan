## Field design of a multi-indicator cluster plan, read on the print

# Closed form (oracle rung 1): two-stage cv of a proportion with var_ratio 1,
# relvariance (1 - p) / p inflated by 1 + (m - 1) * icc over n_psu * m units.
cv_two_stage <- function(p, icc, n_psu, m) {
  sqrt((1 - p) / p * (1 + (m - 1) * icc) / (n_psu * m))
}

coincide <- data.frame(
  name = c("stunting", "anemia"),
  p = c(0.30, 0.10),
  moe = c(0.05, 0.03),
  icc_psu = c(0.02, 0.05)
)
diverge <- data.frame(
  name = c("a", "b"),
  p = c(0.35, 0.40),
  cv = c(0.15, 0.03),
  icc_psu = 0.05
)

table_rows <- function(out, names) {
  rows <- out[seq(which(grepl("^ name", out)) + 1L, length(out))]
  fields <- strsplit(trimws(rows), " +")
  stats::setNames(fields, vapply(fields, `[`, "", 1L))[names]
}

test_that("the indicator table reads every cv at the field design", {
  fit <- n_cluster(stage_cost = c(500, 50), indicators = coincide,
                   budget = 1e5)
  n <- fit$operational$n
  expected <- cv_two_stage(coincide$p, coincide$icc_psu, n[["n_psu"]],
                           n[["n_per_psu"]])
  expect_equal(fit$operational$cv_by_target, expected, tolerance = 1e-10)
  out <- capture.output(print(fit))
  rows <- table_rows(out, coincide$name)
  expect_identical(
    vapply(rows, `[`, "", 4L),
    stats::setNames(sprintf("%.4f", expected), coincide$name)
  )
  expect_false(isTRUE(all.equal(fit$detail$.cv_achieved, expected)))
})

test_that("worst cv and binding indicator are reported apart", {
  fit <- n_cluster(stage_cost = c(500, 50), indicators = diverge,
                   budget = 1e5)
  n <- fit$operational$n
  cv <- cv_two_stage(diverge$p, diverge$icc_psu, n[["n_psu"]],
                     n[["n_per_psu"]])
  expect_identical(which.max(cv), 1L)
  expect_identical(which.max(cv / diverge$cv), 2L)
  expect_identical(fit$operational$binding, "b")

  out <- capture.output(print(fit))
  expect_true(sprintf("worst cv = %.4f (a), cost = 100000", cv[1L]) %in% out)
  expect_true("binding: b" %in% out)
  expect_true("short of target: b" %in% out)
  rows <- table_rows(out, diverge$name)
  expect_identical(rows$b[length(rows$b)], "*")
  expect_identical(length(rows$a), 4L)
})

test_that("a design that meets every target names no shortfall", {
  fit <- n_cluster(stage_cost = c(500, 50), indicators = coincide,
                   budget = 1e5)
  out <- capture.output(print(fit))
  expect_true(any(grepl("^worst cv = [0-9.]+ \\(anemia\\)", out)))
  expect_true("binding: anemia" %in% out)
  expect_false(any(grepl("^short of target", out)))

  cv_mode <- n_cluster(stage_cost = c(500, 50),
                       indicators = transform(coincide, moe = NULL,
                                              cv = c(0.08, 0.15)))
  expect_true(all(cv_mode$operational$cv_by_target <= c(0.08, 0.15)))
  expect_false(any(grepl("^short of target",
                         capture.output(print(cv_mode)))))
})

test_that("every no-domain cluster route stores the field binding", {
  rows <- transform(coincide, moe = NULL, cv = c(0.08, 0.15))
  fits <- list(
    two_stage = n_cluster(stage_cost = c(500, 50), indicators = rows),
    three_stage = n_cluster(
      stage_cost = c(500, 100, 50),
      indicators = transform(rows, icc_ssu = c(0.05, 0.10))
    ),
    stage_fpc = n_cluster(stage_cost = c(500, 50), indicators = coincide,
                          budget = 1e5, fpc = "stage", N = 30000,
                          N_psu = 300)
  )
  seen <- character(0)
  for (route in names(fits)) {
    op <- fits[[route]]$operational
    target <- fits[[route]]$detail$.cv_target
    expect_identical(
      op$binding,
      fits[[route]]$detail$name[which.max(op$cv_by_target / target)],
      label = route
    )
    seen <- c(seen, route)
  }
  expect_setequal(seen, c("two_stage", "three_stage", "stage_fpc"))
})

## Separate domain quotas field the sum of their whole quotas

test_that("the separate-quota total is the sum of the ceiled quotas", {
  dom <- data.frame(
    name = c("a", "a", "b", "b"),
    region = c("N", "S", "N", "S"),
    p = c(0.30, 0.35, 0.10, 0.12),
    moe = c(0.05, 0.05, 0.03, 0.03)
  )
  z2 <- stats::qnorm(0.975)^2
  quota <- tapply(z2 * dom$p * (1 - dom$p) / dom$moe^2, dom$region, max)
  field <- sum(ceiling(quota))
  expect_gt(field, ceiling(sum(quota)))

  fit <- n_multi(dom, domains = "region")
  expect_equal(fit$n, sum(quota), tolerance = 1e-10)
  expect_identical(as.integer(fit), as.integer(field))
  expect_identical(format(fit), sprintf("svyplan_n [n=%d, multi]", field))
  out <- capture.output(print(fit))
  expect_true(sprintf("n = %d (sum of domain quotas)", field) %in% out)

  natural <- n_multi(transform(dom, share = ifelse(region == "N", 0.4, 0.6)),
                     domains = "region", domain_sampling = "natural")
  expect_identical(as.integer(natural), as.integer(ceiling(natural$n)))
})

## Indicator names identify one row per domain

test_that("a repeated indicator name is refused at every entry point", {
  twice <- data.frame(name = c("a", "a"), p = c(0.3, 0.1))
  msg <- "label more than one row"
  expect_error(n_multi(transform(twice, moe = 0.05)), msg)
  expect_error(prec_multi(transform(twice, n = 400)), msg)
  expect_error(
    n_cluster(stage_cost = c(500, 50), budget = 1e5,
              indicators = transform(twice, moe = 0.05, icc_psu = 0.05)),
    msg
  )
  expect_error(
    prec_cluster(indicators = transform(twice, n = 60, n_per_psu = 12,
                                        icc_psu = 0.05)),
    msg
  )
  expect_error(
    n_multi(data.frame(name = c("a", "a"), p = c(0.3, NA), moe = c(0.05, NA),
                       min_respondents = c(NA, 800))),
    msg
  )
  expect_error(
    n_multi(data.frame(name = "a", region = "N", p = c(0.3, 0.2),
                       moe = 0.05), domains = "region"),
    "of the same domain"
  )
  across <- data.frame(name = "a", region = c("N", "S"), p = c(0.3, 0.2),
                       moe = 0.05)
  expect_s3_class(n_multi(across, domains = "region"), "svyplan_n")
  expect_s3_class(
    prec_multi(transform(across, moe = NULL, n = 400), domains = "region"),
    "svyplan_prec"
  )
})

## The stored cv is the worst continuous cv

test_that("cv is the worst continuous cv, not the binding row's", {
  fit <- n_cluster(stage_cost = c(500, 50), indicators = diverge,
                   budget = 1e5)
  n <- fit$n
  cv <- cv_two_stage(diverge$p, diverge$icc_psu, n[["n_psu"]],
                     n[["n_per_psu"]])
  expect_equal(fit$detail$.cv_achieved, cv, tolerance = 1e-10)
  expect_identical(fit$binding, "b")
  expect_equal(fit$cv, max(fit$detail$.cv_achieved))
  expect_equal(fit$cv, fit$detail$.cv_achieved[fit$detail$name == "a"])
  expect_equal(survey::cv(fit), fit$cv)
  expect_equal(with(fit$detail, .cv_achieved[.binding]), cv[2L],
               tolerance = 1e-10)
})

test_that("cv is the largest precision of the design on every route", {
  rows <- transform(coincide, moe = NULL, cv = c(0.08, 0.15))
  dom <- rbind(transform(rows, region = "N"),
               transform(rows, region = "S", p = c(0.32, 0.12)))
  fits <- list(
    two_stage_cv = n_cluster(stage_cost = c(500, 50), indicators = rows),
    two_stage_budget = n_cluster(stage_cost = c(500, 50),
                                 indicators = diverge, budget = 1e5),
    stage_fpc_diverge = n_cluster(stage_cost = c(500, 50),
                                  indicators = diverge, budget = 1e5,
                                  fpc = "stage", N = 30000, N_psu = 300),
    three_stage_diverge = n_cluster(stage_cost = c(500, 100, 50),
                                    indicators = transform(diverge,
                                                           icc_ssu = 0.1),
                                    budget = 2e5),
    three_stage = n_cluster(stage_cost = c(500, 100, 50),
                            indicators = transform(rows,
                                                   icc_ssu = c(0.05, 0.10))),
    stage_fpc = n_cluster(stage_cost = c(500, 50), indicators = coincide,
                          budget = 1e5, fpc = "stage", N = 30000,
                          N_psu = 300),
    domains = n_cluster(stage_cost = c(500, 50), indicators = dom,
                        domains = "region"),
    joint = n_cluster(stage_cost = c(500, 50), indicators = dom,
                      domains = "region", budget = 2e5,
                      allocation = "joint")
  )
  seen <- character(0)
  for (route in names(fits)) {
    fit <- fits[[route]]
    prec <- prec_cluster(fit)
    expect_equal(fit$cv, max(prec$cv), tolerance = 1e-10, label = route)
    if (is.null(fit$domains)) {
      expect_equal(fit$cv, max(fit$detail$.cv_achieved), label = route)
    } else {
      expect_equal(fit$cv, max(fit$domains$.cv), label = route)
      by_domain <- tapply(prec$cv, prec$detail$region, max)
      expect_equal(fit$domains$.cv,
                   as.vector(by_domain[fit$domains$region]),
                   tolerance = 1e-10, label = route)
    }
    expect_equal(fit$operational$cv,
                 max(fit$operational$cv_by_target %||%
                       fit$operational$domains$.cv), label = route)
    seen <- c(seen, route)
  }
  expect_setequal(seen, names(fits))
  for (route in c("two_stage_budget", "stage_fpc_diverge",
                  "three_stage_diverge")) {
    d <- fits[[route]]$detail
    expect_false(d$.binding[which.max(d$.cv_achieved)], label = route)
  }
})

test_that("one indicator row keeps the scalar cv", {
  one <- data.frame(name = "a", p = 0.3, cv = 0.05, icc_psu = 0.05)
  for (budget in list(NULL, 1e5)) {
    ind <- n_cluster(stage_cost = c(500, 50), indicators = one,
                     budget = budget)
    cv_target <- if (is.null(budget)) 0.05
    scalar <- n_cluster(stage_cost = c(500, 50), icc = 0.05,
                        unit_relvar = 0.7 / 0.3, cv = cv_target,
                        budget = budget)
    expect_equal(ind$cv, scalar$cv, tolerance = 1e-6)
  }
})
