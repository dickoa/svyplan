test_that("clustering component follows the two-stage variance model", {
  expect_equal(as.double(design_effect(icc = 0.05, n_per_psu = 25)), 2.2)
  expect_equal(as.double(design_effect(icc = 0, n_per_psu = 25)), 1)
  expect_equal(as.double(design_effect(icc = 1, n_per_psu = 25)), 25)
  expect_equal(as.double(design_effect(icc = 0.05, n_per_psu = 1)), 1)
  expect_equal(
    as.double(design_effect(icc = 0.05, n_per_psu = 25, var_ratio = 1.5)),
    1.5 * 2.2
  )
})

test_that("clustering component follows the three-stage variance model", {
  d <- c(0.01, 0.05)
  # var_ratio_ssu is the within-PSU counterpart of var_ratio_psu, so the pair is not free.
  var_ratio <- c(1.2, 1.2 * (1 - d[1]))
  expected <- var_ratio[1] * d[1] * 10 * 4 + var_ratio[2] * (1 + d[2] * 3)
  expect_equal(
    as.double(design_effect(icc = d, n_per_psu = 10, n_per_ssu = 4, var_ratio = var_ratio)),
    expected
  )
  # With var_ratio unset, var_ratio_ssu follows from var_ratio_psu * (1 - icc_psu): the
  # decomposition leaves it no freedom.
  expect_equal(
    as.double(design_effect(icc = d, n_per_psu = 10, n_per_ssu = 4)),
    d[1] * 40 + (1 - d[1]) * (1 + d[2] * 3)
  )
})

test_that("clustering component reproduces prec_cluster precision exactly", {
  specs <- list(
    list(cost = c(500, 50), icc = 0.05, var_ratio = 1, unit_relvar = 1),
    list(cost = c(500, 100, 50), icc = c(0.01, 0.05), var_ratio = 1, unit_relvar = 1),
    list(cost = c(500, 100, 50), icc = c(0.01, 0.05),
         var_ratio = c(1.2, 1.2 * 0.99), unit_relvar = 2)
  )
  for (spec in specs) {
    plan <- n_cluster(
      stage_cost = spec$cost, icc = spec$icc, var_ratio = spec$var_ratio,
      unit_relvar = spec$unit_relvar, cv = 0.05
    )
    deff <- as.double(design_effect(plan))
    expect_equal(deff, plan$cv^2 * plan$total_n / spec$unit_relvar)
    expect_equal(effective_n(plan), spec$unit_relvar / plan$cv^2)
  }
})

test_that("named icc and var_ratio vectors are reordered to canonical stage order", {
  straight <- design_effect(
    icc = c(0.01, 0.05), n_per_psu = 10, n_per_ssu = 4,
    var_ratio = c(1.2, 1.2 * 0.99)
  )
  swapped <- design_effect(
    icc = c(icc_ssu = 0.05, icc_psu = 0.01),
    n_per_psu = 10, n_per_ssu = 4,
    var_ratio = c(var_ratio_ssu = 1.2 * 0.99, var_ratio_psu = 1.2)
  )
  expect_equal(as.double(swapped), as.double(straight))
})

test_that("clustering component validates its inputs", {
  expect_error(design_effect(icc = 0.05), "both 'icc' and 'n_per_psu'")
  expect_error(design_effect(n_per_psu = 25), "both 'icc' and 'n_per_psu'")
  expect_error(design_effect(icc = 1.2, n_per_psu = 25), "in \\[0, 1\\]")
  expect_error(design_effect(icc = -0.1, n_per_psu = 25), "in \\[0, 1\\]")
  expect_error(
    design_effect(icc = c(0.1, 0.2, 0.3), n_per_psu = 5, n_per_ssu = 2),
    "length 1 .* or 2"
  )
  expect_error(
    design_effect(icc = c(0.01, 0.05), n_per_psu = 10),
    "also needs 'n_per_ssu'"
  )
  expect_error(
    design_effect(icc = 0.05, n_per_psu = 10, n_per_ssu = 4),
    "three stages only"
  )
  expect_error(
    design_effect(icc = 0.05, n_per_psu = 0.5), "at least 1"
  )
  expect_error(
    design_effect(icc = c(0.01, 0.05), n_per_psu = 10, n_per_ssu = 4,
                  var_ratio = c(1, 2, 3)),
    "'var_ratio' must contain 1 or 2"
  )
  expect_error(
    design_effect(icc = 0.05, n_per_psu = 25, var_ratio = -1), "positive finite"
  )
})

test_that("weighting component is Kish's weighting loss", {
  w <- c(1, 1, 1, 1, 5)
  expect_equal(
    as.double(design_effect(weights = w)),
    length(w) * sum(w^2) / sum(w)^2
  )
  expect_equal(as.double(design_effect(weights = rep(3, 100))), 1)
  expect_equal(
    as.double(design_effect(weights = w)),
    as.double(design_effect(weights = 7 * w))
  )
  set.seed(9)
  expect_gte(as.double(design_effect(weights = runif(50, 1, 9))), 1)
})

test_that("weighting component from an allocation matches expanded weights", {
  N <- c(50000, 120000)
  n <- c(600, 400)
  w <- rep(N / n, times = n)
  expect_equal(
    as.double(design_effect(strata = data.frame(N = N, n = n))),
    as.double(design_effect(weights = w))
  )
})

test_that("proportional allocation costs nothing in weighting", {
  N <- c(4000, 3000, 3000)
  prop <- 600 * N / sum(N)
  expect_equal(
    as.double(design_effect(strata = data.frame(N = N, n = prop))), 1
  )
})

test_that("stratification component is the proportional-allocation gain", {
  frame <- data.frame(N = c(50000, 120000), sd = c(12, 20), mean = c(55, 48))
  share <- frame$N / sum(frame$N)
  within <- sum(share * frame$sd^2)
  mu <- sum(share * frame$mean)
  between <- sum(share * (frame$mean - mu)^2)
  expect_equal(
    as.double(design_effect(strata = frame)), within / (within + between)
  )
  expect_lte(as.double(design_effect(strata = frame)), 1)
})

test_that("stratification gain is 1 when stratum means are equal", {
  frame <- data.frame(N = c(100, 200), sd = c(3, 5), mean = c(7, 7))
  expect_equal(as.double(design_effect(strata = frame)), 1)
})

test_that("var and p are accepted as aliases of sd and mean", {
  a <- design_effect(strata = data.frame(N = c(10, 20), sd = c(2, 3),
                                         mean = c(5, 8)))
  b <- design_effect(strata = data.frame(N = c(10, 20), var = c(4, 9),
                                         p = c(5, 8)))
  expect_equal(as.double(a), as.double(b))
})

test_that("components multiply", {
  frame <- data.frame(N = c(50000, 120000), n = c(600, 400),
                      sd = c(12, 20), mean = c(55, 48))
  full <- design_effect(icc = 0.05, n_per_psu = 25, strata = frame)
  parts <- attr(full, "components")
  expect_named(parts, c("cluster", "weight", "strata"))
  expect_equal(as.double(full), prod(parts))
  expect_equal(
    unname(parts[["cluster"]]),
    as.double(design_effect(icc = 0.05, n_per_psu = 25))
  )
})

test_that("the weighting component may only be supplied once", {
  frame <- data.frame(N = c(10, 20), n = c(5, 5))
  expect_error(
    design_effect(weights = c(1, 2), strata = frame),
    "either 'weights' or 'n' in 'strata'"
  )
})

test_that("at least one component is required", {
  expect_error(design_effect(), "no design components")
})

test_that("strata table is validated", {
  expect_error(design_effect(strata = list(N = 1)), "non-empty data frame")
  expect_error(
    design_effect(strata = data.frame(n = 1:2)), "must have an 'N' column"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, 20))),
    "needs 'n' .* or 'sd' and 'mean'"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, 20), n = c(20, 5))),
    "cannot exceed"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, -1), n = c(5, 1))),
    "must be positive"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, 20), n = c(NA, 5))),
    "finite non-missing"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, 20), sd = c(-1, 2),
                                      mean = c(1, 2))),
    "must be non-negative"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, 20), sd = c(1, 2))),
    "needs 'mean'"
  )
  expect_error(
    design_effect(strata = data.frame(N = c(10, 20), sd = c(0, 0),
                                      mean = c(0, 0))),
    "no variability"
  )
})

test_that("weights are validated", {
  expect_error(design_effect(weights = c(1, 0, 2)), "only positive")
  expect_error(design_effect(weights = c(1, NA)), "must not contain NA")
  expect_error(design_effect(weights = c(1, Inf)), "only finite")
  expect_error(design_effect(weights = c(1, -Inf)), "only finite")
  expect_error(design_effect(weights = character(0)), "non-empty numeric")
})

test_that("a bare numeric first argument is rejected with a pointer", {
  expect_error(design_effect(runif(10, 1, 5)), "components by name")
  expect_error(design_effect(runif(10, 1, 5)), "survey::svymean")
})

test_that("varcomp supplies icc and var_ratio, the design supplies the takes", {
  set.seed(314)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  vc <- varcomp(income ~ district, data = frame)
  expect_equal(
    as.double(design_effect(vc, n_per_psu = 25)),
    unname(vc$var_ratio * (1 + vc$icc * 24))
  )
  expect_error(design_effect(vc), "both 'icc' and 'n_per_psu'")
})

test_that("three-stage varcomp needs both takes", {
  # A frame with genuine nesting, so the estimated components are stable
  # enough to be internally consistent.
  set.seed(1618)
  district <- rep(1:20, each = 60)
  village <- rep(1:200, each = 6)
  frame <- data.frame(
    income = 50000 + rnorm(20, 0, 4000)[district] +
      rnorm(200, 0, 3000)[village] + rnorm(1200, 0, 9000),
    district = district,
    village = village
  )
  vc <- varcomp(income ~ district / village, data = frame)
  expect_error(design_effect(vc, n_per_psu = 5), "also needs 'n_per_ssu'")
  expect_s3_class(design_effect(vc, n_per_psu = 5, n_per_ssu = 4),
                  "svyplan_deff")
})

test_that("stratified varcomp is rejected with a pointer to its table", {
  set.seed(11)
  d <- data.frame(
    y = rnorm(120),
    ea = rep(1:12, each = 10),
    region = rep(c("N", "S"), each = 60)
  )
  vc <- varcomp(y ~ ea, data = d, strata = ~region)
  expect_error(design_effect(vc, n_per_psu = 10), "stratified varcomp")
})

test_that("an allocation returns its own variance ratio, not a decomposition", {
  frame <- data.frame(
    stratum = c("A", "B", "C"),
    N = c(4000, 3000, 3000),
    sd = c(10, 15, 8),
    mean = c(50, 60, 55)
  )
  alloc <- n_alloc(frame, n = 600)
  deff <- design_effect(alloc)
  parts <- attr(deff, "components")
  expect_named(parts, "allocation")

  d <- alloc$detail
  W <- d$N / sum(d$N)
  total <- sum(W * d$sd^2) + sum(W * (d$mean - sum(W * d$mean))^2)
  expect_equal(
    as.double(deff),
    sum(d$n) * sum(W^2 * d$sd^2 / d$n) / total
  )
})

test_that("extra weights are an additional Kish factor on the allocation", {
  frame <- data.frame(
    stratum = c("A", "B"),
    N = c(1e5, 1e5),
    sd = c(1, 2),
    mean = c(1, 1)
  )
  alloc <- n_alloc(frame, n = 200, alloc = "neyman")
  extra <- rep(c(1, 4), c(150, 50))
  deff <- design_effect(alloc, weights = extra)
  parts <- attr(deff, "components")

  expect_named(parts, c("allocation", "weight"))
  expect_equal(unname(parts[["allocation"]]), as.double(design_effect(alloc)))
  expect_equal(unname(parts[["weight"]]),
               as.double(design_effect(weights = extra)))
  expect_equal(as.double(deff), prod(unname(parts)))
})

test_that("a generalized allocation has no single variance ratio", {
  frame <- data.frame(
    stratum = c("North urban", "North rural", "South urban"),
    region = c("North", "North", "South"),
    N = c(1000, 1800, 1200)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("coverage", "income"), each = 3),
    p = c(0.5, 0.4, 0.6, NA, NA, NA),
    mean = c(NA, NA, NA, 50, 55, 60),
    sd = c(NA, NA, NA, 10, 12, 15)
  )
  targets <- data.frame(
    name = c("coverage", "income"),
    domain = c(".overall", "region"),
    level = c(NA, "North"),
    cv = c(0.06, NA),
    moe = c(NA, 3)
  )
  fit <- n_alloc(frame, measures = measures, targets = targets)
  expect_error(design_effect(fit), "do not share one variance ratio")
  expect_error(effective_n(fit), "do not share one variance ratio")
})

test_that("an allocation without stratum means cannot form a ratio", {
  frame <- data.frame(stratum = c("A", "B"), N = c(1e5, 1e5), sd = c(1, 2))
  alloc <- n_alloc(frame, n = 200, alloc = "neyman")
  expect_error(design_effect(alloc), "not identified")
  expect_error(design_effect(alloc), "prec_alloc")
})

test_that("a clustered allocation contributes its stratum cluster factors", {
  frame <- data.frame(
    stratum = c("Urban", "Rural"),
    N = c(50000, 150000),
    sd = c(0.45, 0.48),
    mean = c(0.35, 0.25),
    icc_psu = c(0.03, 0.08),
    cost_psu = c(300, 600),
    cost_ssu = c(40, 60)
  )
  alloc <- n_alloc(frame, cv = 0.05)
  deff <- design_effect(alloc)
  d <- alloc$detail
  # the per-stratum factors enter the ratio stratum by stratum, not averaged
  factor_h <- 1 + frame$icc_psu * (d$n_per_psu - 1)
  W <- d$N / sum(d$N)
  total <- sum(W * d$sd^2) + sum(W * (d$mean - sum(W * d$mean))^2)
  expect_equal(
    as.double(deff),
    sum(d$n) * sum(W^2 * d$sd^2 * factor_h / d$n) / total
  )
})

test_that("a clustered allocation deff is the exact variance ratio when sd agree", {
  # with equal stratum sd the Kish product is exact, whatever the sizes,
  # means, takes and var_ratio are
  base <- data.frame(
    stratum = c("A", "B"),
    N = c(3e5, 1.7e6),
    sd = c(2, 2),
    mean = c(10, 14),
    icc_psu = c(0.01, 0.20),
    var_ratio_psu = c(1.3, 0.8),
    n_per_psu = 10,
    cost_psu = 100,
    cost_ssu = 1
  )
  for (tweak in list(
    base,
    transform(base, mean = c(12, 12)),
    transform(base, N = c(1e6, 1e6)),
    transform(base, var_ratio_psu = c(1, 1))
  )) {
    fit <- n_alloc(tweak, n = 2000, alloc = "neyman")
    d <- fit$detail
    W <- d$N / sum(d$N)
    D_h <- tweak$var_ratio_psu * (1 + tweak$icc_psu * (d$n_per_psu - 1))
    S2 <- sum(W * d$sd^2) + sum(W * (d$mean - sum(W * d$mean))^2)
    exact <- sum(d$n) * sum(W^2 * d$sd^2 * D_h / d$n) / S2
    expect_equal(as.double(design_effect(fit)), exact, tolerance = 1e-8)
  }
})

test_that("non-allocation sizing results are rejected", {
  expect_error(design_effect(n_prop(p = 0.3, moe = 0.05)),
               "reads a stratified allocation")
})

test_that("design effects can be used directly and combined", {
  deff <- design_effect(icc = 0.05, n_per_psu = 20)
  expect_s3_class(deff, "svyplan_deff")
  expect_true(is.numeric(deff))
  expect_length(deff, 1L)
  expect_equal(as.double(deff), 1.95)

  by_deff <- n_prop(p = 0.3, moe = 0.05, deff = deff)
  by_number <- n_prop(p = 0.3, moe = 0.05, deff = 1.95)
  expect_equal(by_deff$n, by_number$n)

  expect_equal(deff * 2, 3.9)
  expect_false(inherits(deff * 2, "svyplan_deff"))
  expect_equal(sqrt(deff), sqrt(1.95))
  expect_false(inherits(sqrt(deff), "svyplan_deff"))
  expect_true(deff > 1)
})

test_that("display and export methods are explicit", {
  deff <- design_effect(
    icc = 0.05, n_per_psu = 20,
    strata = data.frame(N = c(10, 20), n = c(5, 5), sd = c(2, 3),
                        mean = c(5, 8))
  )
  expect_output(print(deff), "Design effect \\(planning\\)")
  expect_output(print(deff), "clustering")
  expect_output(print(deff), "stratification")
  expect_output(print(deff), "overall")
  expect_match(format(deff), "^svyplan_deff \\[")

  df <- as.data.frame(deff)
  expect_equal(nrow(df), 1L)
  expect_named(df, c("deff", "deff_cluster", "deff_weight", "deff_strata"))
  expect_equal(df$deff, as.double(deff))
  expect_equal(df$deff, df$deff_cluster * df$deff_weight * df$deff_strata)
})

test_that("unused arguments are rejected", {
  expect_error(design_effect(icc = 0.05, n_per_psu = 20, methd = "kish"),
               "unused argument.*methd")
  expect_error(
    design_effect(n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05),
                  methd = "kish"),
    "unused argument.*methd"
  )
})

test_that("the three-stage factor reduces to var_ratio_psu without clustering", {
  # D = var_ratio_psu icc_psu m q + var_ratio_ssu (1 + icc_ssu (q - 1)). At m = q = 1 no
  # clustering is left, so D must be var_ratio_psu. That holds only when
  # var_ratio_ssu = var_ratio_psu (1 - icc_psu), which is why var_ratio_ssu is derived rather than
  # defaulted to 1.
  for (d1 in c(0.01, 0.05, 0.2304)) {
    for (kp in c(1, 1.4)) {
      expect_equal(
        as.double(design_effect(icc = c(d1, 0.15), n_per_psu = 1,
                                n_per_ssu = 1, var_ratio = kp)),
        kp
      )
    }
  }
})

test_that("the three-stage factor matches the exact variance decomposition", {
  # With components (sigma1^2, sigma2^2, sigma3^2) summing to S^2, the design
  # effect of an a x m x q design is 1 + d1(mq - 1) + d2(q - 1) on
  # total-referenced deltas. The package's within-PSU-referenced icc_ssu is
  # d2 / (1 - d1).
  s1 <- 1.2; s2 <- 0.9; s3 <- 2.0
  total <- s1^2 + s2^2 + s3^2
  d1 <- s1^2 / total
  d2 <- s2^2 / total
  icc_ssu <- d2 / (1 - d1)

  for (spec in list(c(m = 1, q = 1), c(m = 4, q = 3), c(m = 8, q = 5),
                    c(m = 10, q = 1))) {
    m <- spec[["m"]]; q <- spec[["q"]]
    expect_equal(
      as.double(design_effect(icc = c(d1, icc_ssu), n_per_psu = m,
                              n_per_ssu = q)),
      1 + d1 * (m * q - 1) + d2 * (q - 1),
      tolerance = 1e-10
    )
  }
})

test_that("a icc_psu of 1 leaves no three-stage design to identify", {
  expect_error(
    design_effect(icc = c(1, 0.15), n_per_psu = 4, n_per_ssu = 3),
    "leaves no within-PSU variance"
  )
})

test_that("as.list gives the same fields as as.data.frame", {
  deff <- design_effect(
    icc = 0.05, n_per_psu = 20,
    strata = data.frame(N = c(10, 20), n = c(5, 5), sd = c(2, 3),
                        mean = c(5, 8))
  )

  lst <- as.list(deff)
  expect_type(lst, "list")
  expect_identical(lst, as.list(as.data.frame(deff)))
  expect_named(lst, c("deff", "deff_cluster", "deff_weight", "deff_strata"))
  expect_equal(lst$deff, as.double(deff))
  expect_equal(lst$deff, lst$deff_cluster * lst$deff_weight * lst$deff_strata)
})

test_that("as.list on a one-component design effect repeats the overall value", {
  lst <- as.list(design_effect(icc = 0.05, n_per_psu = 25))
  expect_named(lst, c("deff", "deff_cluster"))
  expect_equal(lst$deff, 2.2)
  expect_equal(lst$deff_cluster, 2.2)
})

test_that("components are reachable by name through $ and [[", {
  deff <- design_effect(
    icc = 0.05, n_per_psu = 20,
    strata = data.frame(N = c(10, 20), n = c(5, 5), sd = c(2, 3),
                        mean = c(5, 8))
  )

  expect_equal(deff$deff, as.double(deff))
  expect_equal(deff$deff_cluster, 1.95)
  expect_equal(deff$deff_strata, as.list(deff)$deff_strata)

  # every access route agrees on names and values
  for (nm in names(as.list(deff))) {
    expect_equal(deff[[nm]], as.list(deff)[[nm]])
    expect_equal(deff[[nm]], as.data.frame(deff)[[nm]])
  }

  # a numeric index still reads the underlying vector
  expect_equal(deff[[1L]], as.double(deff))
})

test_that("naming a field a design effect lacks is an error", {
  deff <- design_effect(icc = 0.05, n_per_psu = 25)
  expect_error(deff$deff_strata, "no field 'deff_strata'")
  expect_error(deff$deff_strata, "available: deff, deff_cluster")
  expect_error(deff[["nope"]], "no field 'nope'")
})

test_that("as.list rejects unused arguments", {
  deff <- design_effect(icc = 0.05, n_per_psu = 25)
  expect_error(as.list(deff, methd = "kish"), "unused argument.*methd")
})

test_that("design_effect and effective_n read a prec_cluster result", {
  p <- prec_cluster(n = c(50, 12), icc = 0.05)

  deff <- design_effect(p)
  expect_s3_class(deff, "svyplan_deff")
  expect_equal(as.double(deff), 1 + 0.05 * (12 - 1))
  expect_equal(effective_n(p), 600 / (1 + 0.05 * 11))

  # the n_cluster route agrees on the same design
  expect_equal(as.double(design_effect(p)), as.double(deff))
})

test_that("other prec types carry the deff you supplied, and are rejected", {
  pp <- prec_prop(p = 0.3, n = 500)
  expect_error(design_effect(pp), "reads a cluster design")
  expect_error(effective_n(pp), "reads a cluster design")
})
