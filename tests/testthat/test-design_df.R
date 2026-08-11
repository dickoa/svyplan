## T1. Counting rules, one per design shape

alloc_cluster_frame <- function() {
  data.frame(
    stratum = c("a", "b", "c"),
    N = c(12000, 30000, 8000),
    sd = c(5, 8, 3),
    mean = c(10, 20, 5),
    icc_psu = 0.05,
    n_per_psu = 12
  )
}

alloc_element_frame <- function() {
  data.frame(
    stratum = c("a", "b", "c"),
    N = c(1000, 2000, 500),
    sd = c(5, 8, 3),
    mean = c(10, 20, 5)
  )
}

test_that("an unstratified element design counts its units", {
  res <- n_prop(p = 0.3, moe = 0.05)
  d <- design_df(res)
  expect_s3_class(d, "svyplan_df")
  expect_equal(as.double(d), ceiling(res$n) - 1)
  expect_identical(d$stage, "element")
  expect_identical(d$n_strata, 1L)
  expect_null(d$strata)
  expect_null(d$domains)
})

test_that("an unstratified cluster design counts its whole PSUs", {
  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  d <- design_df(plan)
  expect_equal(as.double(d), plan$operational$n[["n_psu"]] - 1)
  expect_identical(d$stage, "psu")
  # the whole-unit count, never the continuous optimum
  expect_false(isTRUE(all.equal(as.double(d), plan$n[["n_psu"]] - 1)))
})

test_that("a stratified element allocation counts units minus strata", {
  alloc <- n_alloc(alloc_element_frame(), n = 300)
  d <- design_df(alloc)
  expect_equal(as.double(d), sum(alloc$detail$n_int) - 3)
  expect_identical(d$stage, "element")
})

test_that("a stratified cluster allocation counts PSUs minus strata", {
  alloc <- n_alloc(alloc_cluster_frame(), n = 3000)
  d <- design_df(alloc)
  expect_equal(as.double(d), sum(alloc$detail$n_psu_int) - 3)
  expect_identical(d$stage, "psu")
  expect_equal(d$n_units, sum(alloc$detail$n_psu_int))
})

test_that("a stratification counts its allocation", {
  set.seed(9)
  sb <- strata_bound(rlnorm(500, 5, 1), n_strata = 3, n = 200)
  d <- design_df(sb)
  expect_equal(as.double(d), sum(sb$strata$n) - 3)
  expect_identical(d$stage, "element")
})

test_that("a two-phase allocation counts its phase-2 units", {
  frame <- data.frame(
    stratum = c("a", "b"), N = c(1000, 2000), sd = c(5, 8),
    mean = c(10, 20), unit_cost = c(20, 20)
  )
  tp <- n_twophase(frame, phase1_cost = 2, budget = 20000)
  d <- design_df(tp)
  expect_equal(as.double(d), sum(tp$detail$n_int) - 2)
})

test_that("counts can be given directly with no plan", {
  expect_equal(as.double(design_df(n_psu = 300, n_strata = 20)), 280)
  expect_equal(as.double(design_df(n = 1200)), 1199)
  expect_identical(design_df(n_psu = 300)$stage, "psu")
  expect_identical(design_df(n = 300)$stage, "element")

  expect_error(design_df(n_psu = 300, n = 1200), "exactly one of")
  expect_error(design_df(), "exactly one of")
  expect_error(design_df(n_psu = 20, n_strata = 25), "cannot exceed")
  expect_error(design_df(n_psu = 0), "must be an integer >= 1")
  expect_error(design_df(1:3, n_psu = 10), "must be a svyplan result")
})

## T2. The tables and the identities that tie them to the scalar

test_that("the per-stratum table sums to the scalar", {
  alloc <- n_alloc(alloc_cluster_frame(), n = 3000)
  d <- design_df(alloc)
  expect_identical(
    names(d$strata), c("stratum", "n_units", "df", ".status")
  )
  expect_equal(sum(d$strata$df), as.double(d))
  expect_true(all(d$strata$.status == "ok"))
})

test_that("per-domain df is exact and sums to the scalar when domains partition", {
  frame <- alloc_cluster_frame()
  frame$province <- c("N", "N", "S")
  alloc <- n_alloc(frame, n = 3000, domains = "province")
  d <- design_df(alloc)

  expect_identical(
    names(d$domains), c("province", ".domain", ".n_units", ".df", ".status")
  )
  expect_equal(sum(d$domains$.df), as.double(d))

  # a domain's df is its own strata's contribution, computed independently
  north <- d$strata[c(1, 2), ]
  expect_equal(
    d$domains$.df[d$domains$province == "N"],
    sum(north$n_units) - 2
  )
})

test_that("domains covering part of the frame give less than the design df", {
  frame <- alloc_cluster_frame()
  frame$province <- c("N", "N", "S")
  alloc <- n_alloc(frame, n = 3000, domains = "province")
  d <- design_df(alloc)
  expect_true(all(d$domains$.df < as.double(d)))
})

test_that("a take-all stratum is a census in element mode and drops out", {
  frame <- data.frame(
    stratum = c("a", "b", "c"), N = c(100, 2000, 500),
    sd = c(5, 8, 3), mean = c(10, 20, 5),
    take_all = c(TRUE, FALSE, FALSE)
  )
  alloc <- n_alloc(frame, n = 400)
  d <- design_df(alloc)
  expect_identical(d$strata$.status, c("census", "ok", "ok"))
  expect_equal(d$strata$df[1L], 0)
  expect_identical(d$n_strata, 2L)
  expect_equal(as.double(d), sum(alloc$detail$n_int[2:3]) - 2)
  expect_equal(sum(d$strata$df), as.double(d))
})

test_that("take_all is refused in cluster mode rather than half-honoured", {
  # It used to be accepted and treated as an element-level census while the
  # PSU stage kept sampling, which is not a design: the within-PSU take still
  # applies, so no stratum is enumerated.
  frame <- alloc_cluster_frame()
  frame$N[1L] <- 600
  frame$take_all <- c(TRUE, FALSE, FALSE)
  expect_error(n_alloc(frame, n = 3000), "take_all.*not supported")
})

test_that("cluster strata contribute their PSUs minus one to design df", {
  frame <- alloc_cluster_frame()
  frame$N[1L] <- 600
  alloc <- n_alloc(frame, n = 3000)
  d <- design_df(alloc)
  expect_false(any(d$strata$.status == "census"))
  expect_equal(as.double(d), sum(alloc$detail$n_psu_int) - 3)
})

## T3. Singleton strata

test_that("a one-PSU stratum is marked and warns naming the stratum", {
  frame <- data.frame(
    stratum = c("a", "b", "c"), N = c(200, 60000, 8000),
    sd = c(1, 20, 3), mean = c(10, 20, 5),
    icc_psu = 0.05, n_per_psu = 40
  )
  expect_warning(alloc <- n_alloc(frame, n = 3000), "single PSU")
  expect_warning(d <- design_df(alloc), "'a', 'c'")
  expect_identical(d$strata$.status, c("singleton", "ok", "singleton"))
  expect_equal(d$strata$df, c(0, d$strata$n_units[2L] - 1, 0))
  expect_equal(sum(d$strata$df), as.double(d))

  # it fires whether or not the caller ever looks at the table
  expect_warning(design_df(alloc), "no within-stratum variance")
})

test_that("a whole design on a single PSU warns at construction", {
  expect_warning(
    n_cluster(stage_cost = c(5000, 50), icc = 0.05, budget = 9000),
    "single PSU"
  )
  expect_warning(
    n_multi_cluster(
      data.frame(name = "a", p = 0.3, cv = 0.5, icc_psu = 0.05),
      stage_cost = c(5000, 50), budget = 9000
    ),
    "single PSU"
  )
})

test_that("a small but positive df is reported, never warned about", {
  frame <- data.frame(
    stratum = c("a", "b"), N = c(4000, 6000), sd = c(5, 8),
    mean = c(10, 20), icc_psu = 0.05, n_per_psu = 30
  )
  alloc <- n_alloc(frame, n = 300)
  expect_silent(d <- design_df(alloc))
  expect_lt(as.double(d), 20)
  expect_match(paste(capture.output(print(d)), collapse = " "), "df = ")
})

## T4. The object behaves as the number it is

test_that("svyplan_df is drop-in numeric", {
  alloc <- n_alloc(alloc_cluster_frame(), n = 3000)
  d <- design_df(alloc)
  value <- sum(alloc$detail$n_psu_int) - 3

  expect_true(is.numeric(d))
  expect_length(d, 1L)
  expect_equal(as.double(d), value)
  expect_equal(d + 1, value + 1)
  expect_false(inherits(d + 1, "svyplan_df"))
  expect_equal(sqrt(d), sqrt(value))
  expect_true(d > 10)
  expect_equal(qt(0.975, d), qt(0.975, value))
  expect_s3_class(n_prop(p = 0.3, moe = 0.05, df = d), "svyplan_n")
})

test_that("svyplan_df fields are reached by name and unknown ones error", {
  d <- design_df(n_psu = 300, n_strata = 20)
  expect_equal(d$df, 280)
  expect_equal(d[["df"]], 280)
  expect_equal(d[[1L]], 280)
  expect_equal(as.list(d)$n_units, 300)
  expect_error(d$deff, "no field 'deff' in a design df")
  expect_error(d$deff, "available: df, n_units, n_strata, stage")

  df_row <- as.data.frame(d)
  expect_identical(nrow(df_row), 1L)
  expect_identical(names(df_row), c("df", "n_units", "n_strata", "stage"))
  expect_match(format(d), "svyplan_df \\[psu, 280\\]")
})

test_that("summary replaces the meaningless scalar numeric summary", {
  d <- design_df(n_psu = 300, n_strata = 20)
  s <- summary(d)

  expect_s3_class(s, "summary.svyplan_df")
  expect_identical(
    names(s),
    c("df", "n_units", "n_strata", "stage", "strata", "domains", "basis")
  )
  expect_equal(s$df, 280)
  expect_equal(s$n_units, 300)
  expect_identical(s$n_strata, 20L)
  expect_identical(s$stage, "psu")
  expect_null(s$strata)
  expect_null(s$domains)
  expect_match(s$basis, "PSUs minus contributing strata")

  out <- capture.output(expect_invisible(print(s)))
  expect_match(out[1L], "Analysis of design degrees of freedom")
  expect_true(any(grepl("Overall", out, fixed = TRUE)))
  expect_true(any(grepl("Counted stage: PSU", out, fixed = TRUE)))
  expect_true(any(grepl("Per-stratum counts were not supplied", out,
                        fixed = TRUE)))
  expect_false(any(grepl("Min.", out, fixed = TRUE)))
  expect_error(summary(d, digits = 2), "unused argument")
})

test_that("summary prints additive strata and separate domain detail", {
  frame <- alloc_cluster_frame()
  frame$province <- c("N", "N", "S")
  d <- design_df(n_alloc(frame, n = 3000, domains = "province"))
  s <- summary(d)

  expect_identical(s$strata, d$strata)
  expect_identical(s$domains, d$domains)
  expect_equal(sum(s$strata$df), s$df)

  out <- capture.output(print(s))
  expect_true(any(grepl("Stratum", out, fixed = TRUE)))
  expect_true(any(grepl("Constraints", out, fixed = TRUE)))
  expect_true(any(grepl("Design df", out, fixed = TRUE)))
  expect_true(any(grepl("Overall", out, fixed = TRUE)))
  expect_true(any(grepl("Domain degrees of freedom", out, fixed = TRUE)))
  expect_true(any(grepl("not additive to the overall row", out,
                        fixed = TRUE)))
  expect_true(any(grepl("province", out, fixed = TRUE)))
})

test_that("summary distinguishes sampled and counted units for a census", {
  frame <- data.frame(
    stratum = c("a", "b", "c"), N = c(100, 2000, 500),
    sd = c(5, 8, 3), mean = c(10, 20, 5),
    take_all = c(TRUE, FALSE, FALSE)
  )
  s <- summary(design_df(n_alloc(frame, n = 400)))
  out <- capture.output(print(s))

  expect_true(any(grepl("Sampled", out, fixed = TRUE)))
  expect_true(any(grepl("Counted", out, fixed = TRUE)))
  expect_true(any(grepl("census", out, fixed = TRUE)))
})

test_that("a multi-indicator result has no single design to count", {
  res <- n_multi(data.frame(name = "a", p = 0.3, moe = 0.05))
  expect_error(design_df(res), "sizes several indicators against one design")
})

## T5. df in the interval arithmetic

test_that("df widens the interval for every proportion method and the mean", {
  for (m in c("wald", "wilson", "logodds", "beta")) {
    plain <- n_prop(p = 0.2, moe = 0.04, method = m)
    tight <- n_prop(p = 0.2, moe = 0.04, method = m, df = 12)
    expect_gt(tight$n, plain$n)
  }
  expect_gt(n_mean(var = 100, moe = 1, df = 12)$n,
            n_mean(var = 100, moe = 1)$n)
})

test_that("df = Inf reproduces the normal result where the quantile is substituted", {
  for (m in c("wald", "wilson", "logodds")) {
    expect_equal(
      n_prop(p = 0.2, moe = 0.04, method = m, df = Inf)$n,
      n_prop(p = 0.2, moe = 0.04, method = m)$n,
      tolerance = 1e-6
    )
  }
  expect_equal(n_mean(var = 100, moe = 1, df = Inf)$n,
               n_mean(var = 100, moe = 1)$n, tolerance = 1e-8)
  expect_equal(prec_prop(p = 0.2, n = 500, df = Inf)$moe,
               prec_prop(p = 0.2, n = 500)$moe, tolerance = 1e-8)
})

test_that("the beta baseline is the SRS df, not infinity", {
  # Korn-Graubard rescales n_eff by the ratio of two t quantiles, and the
  # reference is the SRS value n - 1. So an unset df matches df = n - 1,
  # and df = Inf is a claim of more degrees of freedom than an SRS has,
  # which narrows the interval rather than reproducing it.
  n <- 500
  expect_equal(
    prec_prop(p = 0.2, n = n, method = "beta", df = n - 1)$moe,
    prec_prop(p = 0.2, n = n, method = "beta")$moe,
    tolerance = 1e-10
  )
  expect_lt(
    prec_prop(p = 0.2, n = n, method = "beta", df = Inf)$moe,
    prec_prop(p = 0.2, n = n, method = "beta")$moe
  )
  expect_gt(
    prec_prop(p = 0.2, n = n, method = "beta", df = 20)$moe,
    prec_prop(p = 0.2, n = n, method = "beta")$moe
  )
})

test_that("moe = q * se holds for wald and the mean engine only", {
  wald <- prec_prop(p = 0.2, n = 500, method = "wald", df = 15)
  expect_equal(wald$moe, qt(0.975, 15) * wald$se, tolerance = 1e-10)

  pm <- prec_mean(var = 100, n = 500, df = 15)
  expect_equal(pm$moe, qt(0.975, 15) * pm$se, tolerance = 1e-10)

  # The other three draw a different interval around the same standard error,
  # so their half-width is not q * se and must not be read back as one.
  for (m in c("wilson", "logodds", "beta")) {
    res <- prec_prop(p = 0.2, n = 500, method = m, df = 15)
    expect_equal(res$se, wald$se, tolerance = 1e-12)
    expect_false(isTRUE(all.equal(res$moe, qt(0.975, 15) * res$se)))
  }
})

test_that("cv is invariant to df and to the interval method", {
  # df switches the interval quantile; it does not enter the sampling variance
  for (m in c("wald", "wilson", "logodds", "beta")) {
    expect_equal(prec_prop(p = 0.2, n = 500, method = m, df = 15)$cv,
                 prec_prop(p = 0.2, n = 500, method = m)$cv, tolerance = 1e-12)
  }
  expect_equal(prec_mean(var = 100, mu = 20, n = 500, df = 15)$cv,
               prec_mean(var = 100, mu = 20, n = 500)$cv, tolerance = 1e-12)

  base <- prec_prop(p = 0.2, n = 500, method = "wald", df = 15)$cv
  for (m in c("wilson", "logodds", "beta")) {
    expect_equal(prec_prop(p = 0.2, n = 500, method = m, df = 15)$cv, base,
                 tolerance = 1e-12)
  }
})

test_that("Korn-Graubard and Wald agree on the direction and rough size of the widening", {
  ratio <- function(m) {
    n_prop(p = 0.05, moe = 0.02, method = m, df = 20)$n /
      n_prop(p = 0.05, moe = 0.02, method = m)$n
  }
  expect_gt(ratio("beta"), 1)
  expect_gt(ratio("wald"), 1)
  expect_equal(ratio("beta"), ratio("wald"), tolerance = 0.1)
})

test_that("a round trip with df set is exact", {
  for (m in c("wald", "wilson", "logodds", "beta")) {
    res <- n_prop(p = 0.2, moe = 0.04, method = m, df = 18)
    back <- prec_prop(res)
    expect_equal(back$params$df, 18)
    expect_equal(back$moe, 0.04, tolerance = 1e-8)
    expect_equal(n_prop(back)$n, res$n, tolerance = 1e-8)
  }
  mres <- n_mean(var = 100, mu = 20, moe = 1, df = 18)
  expect_equal(prec_mean(mres)$params$df, 18)
  expect_equal(prec_mean(mres)$moe, 1, tolerance = 1e-8)
  expect_equal(n_mean(prec_mean(mres))$n, mres$n, tolerance = 1e-8)
})

test_that("confint reads the design's own df", {
  wide <- confint(prec_prop(p = 0.3, n = 200, df = 8))
  narrow <- confint(prec_prop(p = 0.3, n = 200))
  expect_gt(diff(as.numeric(wide)), diff(as.numeric(narrow)))

  wide_mean <- confint(prec_mean(var = 100, mu = 20, n = 200, df = 8))
  narrow_mean <- confint(prec_mean(var = 100, mu = 20, n = 200))
  expect_gt(diff(as.numeric(wide_mean)), diff(as.numeric(narrow_mean)))
})

test_that("df is validated wherever it is accepted", {
  expect_error(n_prop(p = 0.3, moe = 0.05, df = 0.5), "must be a number >= 1")
  expect_error(prec_prop(p = 0.3, n = 100, df = -1), "must be a number >= 1")
  expect_error(n_mean(var = 100, moe = 1, df = 0), "must be a number >= 1")
  expect_error(prec_mean(var = 100, n = 100, df = NA), "must be a number >= 1")
})

test_that("the power functions refuse df and say why", {
  expect_error(power_prop(p1 = 0.3, p2 = 0.4, df = 10),
               "not accepted by the power functions")
  expect_error(power_mean(var = 100, effect = 5, df = 10),
               "a t-based power calculation is a different procedure")
  expect_error(
    power_did(p1_pre = 0.3, p1_post = 0.4, p2_pre = 0.3, p2_post = 0.35,
              df = 10),
    "not accepted by the power functions"
  )
})

## T6. Plans, domains, and the loop back to design_df()

test_that("a plan carries df into the sizing functions", {
  alloc <- n_alloc(alloc_cluster_frame(), n = 3000)
  plan <- svyplan(deff = 1.8, df = design_df(alloc))
  expect_equal(
    n_prop(p = 0.3, moe = 0.05, plan = plan)$n,
    n_prop(p = 0.3, moe = 0.05, deff = 1.8,
           df = as.double(design_df(alloc)))$n
  )
  expect_equal(
    n_mean(var = 100, moe = 1, plan = plan)$n,
    n_mean(var = 100, moe = 1, deff = 1.8,
           df = as.double(design_df(alloc)))$n
  )
  expect_error(svyplan(df = 0.5), "must be a number >= 1")
})

test_that("per-domain df feeds an indicators table row by row", {
  frame <- alloc_cluster_frame()
  frame$province <- c("N", "N", "S")
  alloc <- n_alloc(frame, n = 3000, domains = "province")
  doms <- design_df(alloc)$domains

  ind <- data.frame(
    name = c("a", "b"),
    province = c("N", "S"),
    p = c(0.3, 0.3),
    moe = c(0.05, 0.05)
  )
  ind$df <- doms$.df[match(ind$province, doms$province)]
  res <- n_multi(ind, domains = "province")

  # each row reproduces the single-indicator size at its own domain's df
  expect_equal(
    res$domains$.n[res$domains$province == "N"],
    n_prop(p = 0.3, moe = 0.05, df = doms$.df[doms$province == "N"])$n
  )
  expect_equal(
    res$domains$.n[res$domains$province == "S"],
    n_prop(p = 0.3, moe = 0.05, df = doms$.df[doms$province == "S"])$n
  )
})

test_that("the cluster and allocation prints report their design df", {
  alloc <- n_alloc(alloc_cluster_frame(), n = 3000)
  out <- capture.output(print(alloc))
  expect_true(any(grepl("^design df = ", out)))
  expect_true(any(grepl(sprintf("design df = %g", as.double(design_df(alloc))),
                        out)))

  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  expect_true(any(grepl("^design df = ", capture.output(print(plan)))))
})

test_that("design degrees of freedom cannot be modified in place", {
  f <- design_df(n_psu = 40, n_strata = 4)
  expect_error({f[1] <- 99}, "cannot be modified in place")
  expect_error({f[[1]] <- 99}, "cannot be modified in place")
  expect_equal(as.double(f), 36)
})
