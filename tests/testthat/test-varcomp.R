test_that("varcomp 2-stage SRS formula interface works", {
  set.seed(1031)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  result <- varcomp(income ~ district, data = frame)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)
  expect_true(result$icc >= 0 && result$icc <= 1)
  expect_true(result$var_ratio > 0)
  expect_true(result$unit_relvar > 0)
})

test_that("varcomp 2-stage SRS vector interface matches formula", {
  set.seed(1033)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  result_formula <- varcomp(income ~ district, data = frame)
  result_vector <- varcomp(frame$income, stage_id = list(frame$district))

  expect_equal(
    result_vector$varb,
    result_formula$varb,
    tolerance = 1e-10
  )
  expect_equal(
    result_vector$varw,
    result_formula$varw,
    tolerance = 1e-10
  )
  expect_equal(result_vector$icc, result_formula$icc, tolerance = 1e-10)
})

test_that("unstratified variance components have explicit export schemas", {
  set.seed(1217)
  y <- rnorm(200, 50, 10)
  psu <- rep(seq_len(20), each = 10)
  ssu <- rep(seq_len(100), each = 2)

  two_stage <- varcomp(y, stage_id = list(psu))
  three_stage <- varcomp(y, stage_id = list(psu, ssu))
  two_df <- as.data.frame(two_stage)
  three_df <- as.data.frame(three_stage)

  expect_identical(
    names(two_df),
    c("stages", "varb", "varw", "icc", "var_ratio", "unit_relvar")
  )
  expect_identical(
    names(three_df),
    c(
      "stages", "varb", "varw_psu", "varw_ssu", "icc_psu",
      "icc_ssu", "var_ratio_psu", "var_ratio_ssu", "unit_relvar"
    )
  )
  expect_equal(two_df$icc, two_stage$icc)
  expect_equal(three_df$icc_psu, three_stage$icc[["icc_psu"]])
})

test_that("varcomp 2-stage SRS known values", {
  # Manually construct a known case
  set.seed(123)
  psu_id <- rep(1:10, each = 5)
  y <- rnorm(50, 100, 20)

  result <- varcomp(y, stage_id = list(psu_id))

  # Replicate ANOVA decomposition
  M <- 10L
  Ni <- as.numeric(table(psu_id))
  ti <- as.numeric(by(y, INDICES = psu_id, FUN = sum))
  S2Ui <- as.numeric(by(y, INDICES = psu_id, FUN = var))
  tbarU <- mean(ti)
  tU <- M * tbarU
  S2U1 <- var(ti)
  B2 <- S2U1 / tbarU^2
  W2 <- M * sum(Ni^2 * S2Ui) / tU^2

  expect_equal(result$varb, B2, tolerance = 1e-10)
  expect_equal(result$varw, W2, tolerance = 1e-10)
  expect_equal(result$icc, B2 / (B2 + W2), tolerance = 1e-10)
})

test_that("varcomp 2-stage PPS known values", {
  set.seed(456)
  psu_id <- rep(1:10, each = 5)
  y <- rnorm(50, 100, 20)
  pp <- rep(1 / 10, 10)

  result <- varcomp(y, stage_id = list(psu_id), prob = pp)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)

  # Replicate PPS decomposition
  Ni <- as.numeric(table(psu_id))
  cl_tots <- as.numeric(by(y, INDICES = psu_id, FUN = sum))
  cl_vars <- as.numeric(by(y, INDICES = psu_id, FUN = var))
  tU <- sum(cl_tots)
  S2U1 <- sum(pp * (cl_tots / pp - tU)^2)
  B2 <- S2U1 / tU^2
  W2 <- sum(Ni^2 * cl_vars / pp) / tU^2

  expect_equal(result$varb, B2, tolerance = 1e-10)
  expect_equal(result$varw, W2, tolerance = 1e-10)
})

test_that("varcomp 2-stage PPS with formula prob", {
  set.seed(789)
  frame <- data.frame(
    income = rnorm(50, 50000, 10000),
    district = rep(1:10, each = 5),
    pp = rep(1 / 10, 50)
  )
  result <- varcomp(income ~ district, data = frame, prob = ~pp)
  expect_s3_class(result, "svyplan_varcomp")
})

test_that("varcomp integrates with n_cluster", {
  set.seed(1039)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  vc <- varcomp(income ~ district, data = frame)
  plan <- n_cluster(stage_cost = c(500, 50), icc = vc, budget = 100000)
  expect_s3_class(plan, "svyplan_cluster")
})

test_that("varcomp validates inputs", {
  expect_error(varcomp("not_numeric"), "must be a formula")
  expect_error(varcomp(1:10), "'stage_id' must be a list")
  expect_error(varcomp(~x, data = data.frame(x = 1:10)), "must have a response")
  expect_error(
    varcomp(income ~ district, data = data.frame(x = 1:10)),
    "not found"
  )
})

test_that("varcomp handles lonely SSUs", {
  # One cluster with a single element
  psu_id <- c(1, 2, 2, 3, 3, 3)
  y <- c(10, 20, 25, 30, 35, 40)
  result <- varcomp(y, stage_id = list(psu_id))
  expect_s3_class(result, "svyplan_varcomp")
  expect_false(is.na(result$icc))
})

test_that("varcomp.survey.design matches formula interface", {
  skip_if_not_installed("survey")
  set.seed(1049)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  ref <- varcomp(income ~ district, data = frame)

  dsgn <- survey::svydesign(
    ids = ~district,
    data = frame,
    weights = rep(1, 200)
  )
  result <- varcomp(dsgn, ~income)

  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$varb, ref$varb, tolerance = 1e-10)
  expect_equal(result$varw, ref$varw, tolerance = 1e-10)
  expect_equal(result$icc, ref$icc, tolerance = 1e-10)
})

test_that("varcomp.survey.design works with strata", {
  skip_if_not_installed("survey")
  set.seed(1)
  frame <- data.frame(
    y = rnorm(200, 50, 10),
    cluster = rep(1:20, each = 10),
    stratum = rep(1:4, each = 50)
  )
  dsgn <- survey::svydesign(
    ids = ~cluster,
    strata = ~stratum,
    data = frame,
    weights = rep(1, 200),
    nest = TRUE
  )
  result <- varcomp(dsgn, ~y)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)

  # The design's own strata are used, so the result is per-stratum without
  # having to name them again. Asserting only that an icc is in [0, 1] would
  # pass whether or not they were read.
  expect_false(is.null(result$strata))
  expect_equal(nrow(result$strata), 4L)
  expect_setequal(as.character(result$strata$stratum), as.character(1:4))
  expect_true(all(result$strata$icc_psu >= 0 & result$strata$icc_psu <= 1))

  # and naming them explicitly is the same calculation
  explicit <- varcomp(dsgn, ~y, strata = ~stratum)
  expect_equal(result$strata$icc_psu, explicit$strata$icc_psu)
})

test_that("an explicit strata argument overrides the design's own", {
  skip_if_not_installed("survey")
  set.seed(1)
  frame <- data.frame(
    y = rnorm(200, 50, 10),
    cluster = rep(1:20, each = 10),
    stratum = rep(1:4, each = 50),
    other = rep(c("a", "b"), each = 100)
  )
  dsgn <- survey::svydesign(
    ids = ~cluster, strata = ~stratum, data = frame,
    weights = rep(1, 200), nest = TRUE
  )
  native <- varcomp(dsgn, ~y)
  override <- varcomp(dsgn, ~y, strata = ~other)
  expect_equal(nrow(native$strata), 4L)
  expect_equal(nrow(override$strata), 2L)
  expect_setequal(as.character(override$strata$stratum), c("a", "b"))
})

test_that("varcomp.survey.design feeds into n_cluster", {
  skip_if_not_installed("survey")
  set.seed(2)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  dsgn <- survey::svydesign(
    ids = ~district,
    data = frame,
    weights = rep(1, 200)
  )
  vc <- varcomp(dsgn, ~income)
  plan <- n_cluster(stage_cost = c(500, 50), icc = vc, budget = 100000)
  expect_s3_class(plan, "svyplan_cluster")
})

test_that("3-stage varcomp with non-nested SSU IDs matches nested", {
  set.seed(3)
  psu_id <- rep(1:5, each = 10)
  ssu_id_non_nested <- rep(rep(1:2, each = 5), 5)
  ssu_id_nested <- interaction(psu_id, ssu_id_non_nested, drop = TRUE)
  y <- rnorm(50, 100, 20)
  pp <- rep(0.2, 5)

  res_non <- varcomp(y, stage_id = list(psu_id, ssu_id_non_nested), prob = pp)
  res_nested <- varcomp(y, stage_id = list(psu_id, ssu_id_nested), prob = pp)

  expect_equal(res_non$icc, res_nested$icc, tolerance = 1e-10)
  expect_equal(res_non$varb, res_nested$varb, tolerance = 1e-10)
  expect_equal(res_non$varw, res_nested$varw, tolerance = 1e-10)
})

test_that("all-singleton 2-stage returns icc = 1 with warning", {
  psu_id <- 1:10
  y <- rnorm(10, 100, 20)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "all clusters are singletons"
  )
  expect_equal(res$icc, 1)
  expect_false(is.nan(res$icc))
  expect_true(is.finite(res$var_ratio))
})

test_that("varcomp.survey.design validates inputs", {
  skip_if_not_installed("survey")
  frame <- data.frame(y = rnorm(50), cl = rep(1:10, each = 5))
  dsgn <- survey::svydesign(ids = ~cl, data = frame, weights = rep(1, 50))

  expect_error(varcomp(dsgn), "one-sided formula")
  expect_error(varcomp(dsgn, ~nonexistent), "not found")

  dsgn_nocl <- survey::svydesign(ids = ~1, data = frame, weights = rep(1, 50))
  expect_error(varcomp(dsgn_nocl, ~y), "no clusters")
})

test_that("varcomp rejects mismatched stage_id length", {
  expect_error(
    varcomp(1:10, stage_id = list(1:5)),
    "same as outcome vector"
  )
})

test_that("varcomp rejects empty stage_id", {
  expect_error(
    varcomp(1:10, stage_id = list()),
    "must not be empty"
  )
})

test_that("varcomp 2-stage SRS constant y gives icc = 0 with warning", {
  y <- rep(42, 50)
  psu_id <- rep(1:10, each = 5)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "no variance to split"
  )
  expect_equal(res$icc, 0)
  expect_equal(res$var_ratio, 1)
  expect_true(is.finite(res$icc))
  expect_true(is.finite(res$var_ratio))
})

test_that("varcomp 2-stage PPS constant y gives icc = 0 with warning", {
  y <- rep(42, 50)
  psu_id <- rep(1:10, each = 5)
  pp <- rep(0.1, 10)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id), prob = pp),
    "no variance to split"
  )
  expect_equal(res$icc, 0)
  expect_equal(res$var_ratio, 1)
  expect_true(is.finite(res$icc))
  expect_true(is.finite(res$var_ratio))
})

test_that("varcomp 3-stage PPS constant y gives finite icc with warning", {
  y <- rep(42, 60)
  psu_id <- rep(1:6, each = 10)
  ssu_id <- rep(rep(1:2, each = 5), 6)
  pp <- rep(1 / 6, 6)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id, ssu_id), prob = pp),
    "no variance to split"
  )
  expect_length(res$icc, 2)
  expect_true(all(is.finite(res$icc)))
  expect_true(all(is.finite(res$var_ratio)))
})

test_that("varcomp rejects NA in outcome vector", {
  expect_error(
    varcomp(c(1, 2, NA), stage_id = list(c(1, 1, 2))),
    "must not contain NA"
  )
})

test_that("varcomp rejects empty outcome vector", {
  expect_error(
    varcomp(numeric(0), stage_id = list(integer(0))),
    "non-empty numeric"
  )
})

test_that("varcomp accepts '/' and '%in%' for multi-stage formula", {
  set.seed(1)
  frame <- data.frame(
    y = rnorm(40),
    psu = rep(1:5, each = 8),
    ssu = rep(1:20, each = 2),
    pp = rep(1 / 5, 40)
  )
  expect_error(
    varcomp(y ~ psu + ssu, data = frame, prob = ~pp),
    "must express nesting"
  )
  expect_error(
    varcomp(y ~ psu * ssu, data = frame, prob = ~pp),
    "must express nesting"
  )
  expect_error(
    varcomp(y ~ psu + ssu + psu:ssu, data = frame, prob = ~pp),
    "must express nesting"
  )
  res_slash <- varcomp(y ~ psu/ssu, data = frame, prob = ~pp)
  res_in <- varcomp(y ~ ssu %in% psu, data = frame, prob = ~pp)
  expect_equal(res_slash$icc, res_in$icc, tolerance = 1e-10)
  expect_equal(res_slash$varb, res_in$varb, tolerance = 1e-10)
  expect_equal(res_slash$varw, res_in$varw, tolerance = 1e-10)
})

test_that("3-stage SRS works without prob", {
  set.seed(99)
  frame <- data.frame(
    y = rnorm(400, 50, 10),
    psu = rep(1:20, each = 20),
    ssu = rep(1:100, each = 4)
  )
  vc <- varcomp(y ~ psu/ssu, data = frame)
  expect_s3_class(vc, "svyplan_varcomp")
  expect_equal(vc$stages, 3L)
  expect_length(vc$icc, 2L)
  expect_length(vc$var_ratio, 2L)
  expect_true(all(vc$icc >= 0 & vc$icc <= 1))
})

test_that("3-stage SRS matches PPS with uniform prob", {
  set.seed(99)
  M <- 20
  frame <- data.frame(
    y = rnorm(400, 50, 10),
    psu = rep(1:M, each = 20),
    ssu = rep(1:100, each = 4),
    pp = rep(1 / M, 400)
  )
  vc_srs <- varcomp(y ~ psu/ssu, data = frame)
  vc_pps <- varcomp(y ~ psu/ssu, data = frame, prob = ~pp)
  expect_equal(vc_srs$icc, vc_pps$icc, tolerance = 1e-10)
  expect_equal(vc_srs$var_ratio, vc_pps$var_ratio, tolerance = 1e-10)
  expect_equal(vc_srs$unit_relvar, vc_pps$unit_relvar, tolerance = 1e-10)
  expect_equal(vc_srs$varb, vc_pps$varb, tolerance = 1e-10)
  expect_equal(vc_srs$varw, vc_pps$varw, tolerance = 1e-10)
})

test_that("3-stage SRS works with vector interface", {
  set.seed(99)
  y <- rnorm(400, 50, 10)
  psu <- rep(1:20, each = 20)
  ssu <- rep(1:100, each = 4)
  vc_vec <- varcomp(y, stage_id = list(psu, ssu))
  vc_frm <- varcomp(y ~ psu/ssu, data = data.frame(y, psu, ssu))
  expect_equal(vc_vec$icc, vc_frm$icc, tolerance = 1e-10)
  expect_equal(vc_vec$var_ratio, vc_frm$var_ratio, tolerance = 1e-10)
})

test_that("survey.design method with unit weights matches frame", {
  skip_if_not_installed("survey")
  set.seed(1051)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  ref <- varcomp(income ~ district, data = frame)
  dsgn <- survey::svydesign(
    ids = ~district,
    data = frame,
    weights = rep(1, 200)
  )
  result <- varcomp(dsgn, ~income)
  expect_equal(result$varb, ref$varb, tolerance = 1e-12)
  expect_equal(result$varw, ref$varw, tolerance = 1e-12)
  expect_equal(result$icc, ref$icc, tolerance = 1e-12)
})

test_that("weighted correction shrinks the between component", {
  skip_if_not_installed("survey")
  set.seed(1061)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  ref <- varcomp(income ~ district, data = frame)
  dsgn <- survey::svydesign(
    ids = ~district,
    data = frame,
    weights = rep(7, 200)
  )
  result <- varcomp(dsgn, ~income, weights = rep(7, 200))
  expect_lt(result$varb, ref$varb)
  expect_equal(result$varw, ref$varw, tolerance = 1e-12)
})

test_that("survey.design method applies within-cluster weights", {
  skip_if_not_installed("survey")
  set.seed(7)
  M <- 50
  pop <- data.frame(psu = rep(seq_len(M), each = 80))
  pop$y <- rnorm(nrow(pop), rnorm(M, 100, 15)[pop$psu], 20)
  truth <- varcomp(y ~ psu, data = pop)

  s <- do.call(rbind, lapply(split(pop, pop$psu), function(d) {
    q <- if (d$psu[1] %% 2 == 0) 8L else 24L
    out <- d[sample.int(nrow(d), q), ]
    out$w <- nrow(d) / q
    out
  }))
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = s)
  wtd <- varcomp(dsgn, ~y, weights = ~w)
  unw <- varcomp(s$y, stage_id = list(s$psu))

  expect_lt(abs(wtd$icc - truth$icc), 0.1)
  expect_lt(abs(wtd$icc - truth$icc), abs(unw$icc - truth$icc))
})

test_that("weighted 2-stage PPS recovers the population icc", {
  skip_if_not_installed("survey")
  set.seed(21)
  M <- 40
  sizes <- rep(c(60, 120), length.out = M)
  pop <- data.frame(psu = rep(seq_len(M), times = sizes))
  pop$y <- rnorm(nrow(pop), rnorm(M, 100, 15)[pop$psu], 20)
  pp <- sizes / sum(sizes)
  truth <- varcomp(pop$y, stage_id = list(pop$psu), prob = pp)

  s <- do.call(rbind, lapply(split(pop, pop$psu), function(d) {
    q <- if (nrow(d) > 100) 10L else 30L
    out <- d[sample.int(nrow(d), q), ]
    out$w <- nrow(d) / q
    out
  }))
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = s)
  wtd <- varcomp(dsgn, ~y, prob = pp, weights = ~w)
  unw <- varcomp(s$y, stage_id = list(s$psu), prob = pp)

  expect_lt(abs(wtd$icc - truth$icc), 0.1)
  expect_lt(abs(wtd$icc - truth$icc), abs(unw$icc - truth$icc))
})

test_that("per-stratum weighted results do not depend on other strata", {
  skip_if_not_installed("survey")
  set.seed(1063)
  d <- data.frame(
    y = rnorm(400, rep(c(10, 20), each = 200)),
    psu = rep(1:40, each = 10),
    region = rep(c("N", "S"), each = 200),
    w = rep(c(2, 8), each = 200)
  )
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = d)
  vc <- varcomp(dsgn, ~y, strata = ~region, weights = ~w)

  for (s in c("N", "S")) {
    sub <- d[d$region == s, ]
    dsub <- survey::svydesign(ids = ~psu, weights = ~w, data = sub)
    ref <- varcomp(dsub, ~y, weights = ~w)
    i <- match(s, vc$strata$stratum)
    expect_equal(vc$strata$icc_psu[i], ref$icc)
    expect_equal(vc$strata$var_ratio_psu[i], ref$var_ratio)
  }
})

test_that("3-stage survey.design with unit weights matches frame", {
  skip_if_not_installed("survey")
  set.seed(11)
  M <- 30
  frame <- data.frame(
    psu = rep(1:M, each = 60),
    ssu = rep(rep(1:6, each = 10), times = M)
  )
  frame$y <- rnorm(nrow(frame), rnorm(M, 100, 12)[frame$psu], 15)
  ref <- varcomp(y ~ psu / ssu, data = frame)
  dsgn <- survey::svydesign(
    ids = ~psu + ssu,
    data = frame,
    weights = rep(1, nrow(frame))
  )
  result <- varcomp(dsgn, ~y)
  expect_equal(result$icc, ref$icc, tolerance = 1e-12)
  expect_equal(result$varw, ref$varw, tolerance = 1e-12)
})

test_that("3-stage weighted components stay finite and bounded", {
  skip_if_not_installed("survey")
  set.seed(13)
  M <- 20
  frame <- data.frame(
    psu = rep(1:M, each = 40),
    ssu = rep(rep(1:4, each = 10), times = M)
  )
  frame$y <- rnorm(nrow(frame), rnorm(M, 100, 12)[frame$psu], 15)
  s <- frame[unlist(lapply(split(seq_len(nrow(frame)), frame$psu),
                           function(ix) sample(ix, 24))), ]
  s$w <- 40 / 24
  dsgn <- survey::svydesign(ids = ~psu + ssu, data = s, weights = ~w)
  res <- varcomp(dsgn, ~y, weights = ~w)
  expect_true(all(is.finite(res$icc)))
  expect_true(all(res$icc >= 0 & res$icc <= 1))
  expect_true(all(res$var_ratio > 0))
})

test_that("survey.design method rejects non-finite weights", {
  skip_if_not_installed("survey")
  frame <- data.frame(y = rnorm(50), cl = rep(1:10, each = 5))
  dsgn <- survey::svydesign(ids = ~cl, data = frame, weights = rep(1, 50))
  dsgn$prob <- rep(0, 50)
  expect_error(varcomp(dsgn, ~y), "positive and finite")
})

test_that("strata gives per-stratum components matching manual splits", {
  set.seed(3)
  d <- data.frame(
    region = rep(c("N", "S"), each = 400),
    ea = rep(1:40, each = 20),
    y = rnorm(800, rep(c(50, 70), each = 400), 15)
  )
  vc <- varcomp(y ~ ea, data = d, strata = ~region)
  expect_null(vc$icc)
  expect_equal(nrow(vc$strata), 2L)
  expect_equal(vc$strata$stratum, c("N", "S"))

  for (s in c("N", "S")) {
    sub <- d[d$region == s, ]
    ref <- varcomp(sub$y, stage_id = list(sub$ea))
    i <- match(s, vc$strata$stratum)
    expect_equal(vc$strata$icc_psu[i], ref$icc)
    expect_equal(vc$strata$var_ratio_psu[i], ref$var_ratio)
    expect_equal(vc$strata$unit_relvar[i], ref$unit_relvar)
    expect_equal(vc$strata$sd[i], sd(sub$y))
    expect_equal(vc$strata$mean[i], mean(sub$y))
  }

  expect_identical(as.data.frame(vc), vc$strata)
  expect_output(print(vc), "2-stage, 2 strata")
})

test_that("strata works for 3-stage and the vector interface", {
  set.seed(5)
  d <- data.frame(
    dom = rep(c("A", "B"), each = 240),
    psu = rep(1:24, each = 20),
    ssu = rep(rep(1:4, each = 5), times = 24),
    y = rnorm(480, 100, 20)
  )
  vc <- varcomp(d$y, stage_id = list(d$psu, d$ssu), strata = d$dom)
  expect_equal(vc$stages, 3L)
  expect_true(all(c("icc_psu", "icc_ssu", "var_ratio_psu", "var_ratio_ssu",
                    "varw_psu", "varw_ssu") %in% names(vc$strata)))
})

test_that("strata works on survey.design objects", {
  skip_if_not_installed("survey")
  set.seed(9)
  d <- data.frame(
    region = rep(c("N", "S"), each = 300),
    ea = rep(1:30, each = 20),
    y = rnorm(600, 60, 12)
  )
  dsgn <- survey::svydesign(ids = ~ea, data = d, weights = rep(1, 600))
  vc <- varcomp(dsgn, ~y, strata = ~region)
  ref <- varcomp(y ~ ea, data = d, strata = ~region)
  expect_equal(vc$strata, ref$strata)
})

test_that("stratified varcomp is rejected where a pooled one is needed", {
  set.seed(3)
  d <- data.frame(
    region = rep(c("N", "S"), each = 200),
    ea = rep(1:20, each = 20),
    y = rnorm(400, 50, 10)
  )
  vc <- varcomp(y ~ ea, data = d, strata = ~region)
  expect_error(n_cluster(stage_cost = c(500, 50), icc = vc, cv = 0.05),
               "stratified varcomp")
  expect_error(design_effect(vc, n_per_psu = 10), "stratified varcomp")
  expect_error(prec_cluster(n = c(20, 10), icc = vc),
               "stratified varcomp")
  pooled <- varcomp(y ~ ea, data = d)
  pooled_df <- as.data.frame(pooled)
  expect_equal(nrow(pooled_df), 1L)
  expect_equal(pooled_df$icc, pooled$icc)
})

test_that("strata validates inputs", {
  y <- rnorm(40)
  id <- rep(1:4, each = 10)
  expect_error(varcomp(y, stage_id = list(id), strata = rep("a", 39)),
               "same length")
  expect_error(
    varcomp(y, stage_id = list(rep(1:6, each = 5)),
            strata = rep(c("a", "b"), each = 20)),
    "must have length 40"
  )
  expect_error(varcomp(y, stage_id = list(id), strata = c(NA, rep("a", 39))),
               "NA")
  expect_error(
    varcomp(y, stage_id = list(id), strata = rep(c("a", "b"), each = 20),
            prob = rep(0.25, 4)),
    "one value per observation"
  )
})

test_that("single-PSU data is rejected with a clear message", {
  y <- rnorm(20)
  expect_error(varcomp(y, stage_id = list(rep(1, 20))),
               "at least two PSUs")
  d <- data.frame(
    y = rnorm(60),
    psu = c(rep(1, 20), rep(2:3, each = 20)),
    region = c(rep("Solo", 20), rep("Duo", 40))
  )
  expect_error(varcomp(y ~ psu, data = d, strata = ~region),
               "stratum 'Solo': at least two PSUs")
  expect_error(varcomp(y ~ psu, data = d, strata = ~region),
               "n_alloc\\(psu =\\), leave its icc_psu NA")
  # The register hint is for strata only. An unstratified call has no
  # stratum to leave out.
  err <- tryCatch(varcomp(y, stage_id = list(rep(1, 20))),
                  error = conditionMessage)
  expect_false(grepl("n_alloc", err))
})

test_that("varcomp validates outcomes, stage ids, and probabilities", {
  psu <- rep(1:3, each = 2)
  expect_error(varcomp(c(1, 2, Inf, 4, 5, 6), stage_id = list(psu)),
               "finite")
  expect_error(varcomp(1:6, stage_id = list(c(1, 1, NA, NA, 3, 3))),
               "must not contain NA")
  expect_error(varcomp(1:6, stage_id = list(psu), prob = c(1.2, -0.1, -0.1)),
               "strictly between 0 and 1")
  expect_error(varcomp(1:6, stage_id = list(psu), prob = c(0.5, 0.5, 0)),
               "strictly between 0 and 1")
  expect_error(
    varcomp(1:6, stage_id = list(psu),
            prob = c(0.2, 0.8, 0.3, 0.3, 0.5, 0.5)),
    "constant within each PSU"
  )
})

test_that("named per-PSU probabilities are matched by PSU id", {
  psu <- rep(1:3, each = 2)
  a <- varcomp(1:6, stage_id = list(psu), prob = c("1" = 0.2, "2" = 0.3, "3" = 0.5))
  b <- varcomp(1:6, stage_id = list(psu), prob = c("3" = 0.5, "1" = 0.2, "2" = 0.3))
  expect_equal(a$varb, b$varb)
  expect_equal(a$icc, b$icc)
  expect_error(
    varcomp(1:6, stage_id = list(psu), prob = c("1" = 0.2, "2" = 0.3, "9" = 0.5)),
    "must match the PSU identifiers"
  )
})

test_that("formula and vector 'weights' match the survey.design method", {
  skip_if_not_installed("survey")
  set.seed(11)
  s <- data.frame(
    y = rnorm(120, rep(c(10, 20, 30), each = 40)),
    psu = rep(1:12, each = 10),
    w = rep(c(5, 9), each = 60)
  )
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = s)
  ref <- varcomp(dsgn, ~y, weights = ~w)

  frm <- varcomp(y ~ psu, data = s, weights = ~w)
  expect_equal(frm$varb, ref$varb, tolerance = 1e-12)
  expect_equal(frm$varw, ref$varw, tolerance = 1e-12)
  expect_equal(frm$icc, ref$icc, tolerance = 1e-12)

  vec <- varcomp(s$y, stage_id = list(s$psu), weights = s$w)
  expect_equal(vec$icc, ref$icc, tolerance = 1e-12)
  expect_equal(vec$var_ratio, ref$var_ratio, tolerance = 1e-12)
})

test_that("'weights' works with a PPS first stage and with strata", {
  skip_if_not_installed("survey")
  set.seed(12)
  s <- data.frame(
    y = rnorm(160, rep(c(5, 15, 25, 35), each = 40)),
    psu = rep(1:16, each = 10),
    region = rep(c("N", "S"), each = 80),
    w = rep(c(4, 8), times = 80)
  )
  pp <- rep(1 / 16, 16)
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = s)
  ref <- varcomp(dsgn, ~y, prob = pp, weights = ~w)
  frm <- varcomp(y ~ psu, data = s, weights = ~w, prob = pp)
  expect_equal(frm$varb, ref$varb, tolerance = 1e-12)
  expect_equal(frm$icc, ref$icc, tolerance = 1e-12)

  ref_s <- varcomp(dsgn, ~y, strata = ~region, weights = ~w)
  frm_s <- varcomp(y ~ psu, data = s, weights = ~w, strata = ~region)
  expect_equal(frm_s$strata$icc_psu, ref_s$strata$icc_psu,
               tolerance = 1e-12)
  expect_equal(frm_s$strata$sd, ref_s$strata$sd, tolerance = 1e-12)
})

test_that("unit 'weights' collapse to the exact frame formulas", {
  set.seed(13)
  frame <- data.frame(
    y = rnorm(100, 50, 10),
    psu = rep(1:10, each = 10)
  )
  ref <- varcomp(y ~ psu, data = frame)
  wtd <- varcomp(y ~ psu, data = frame, weights = rep(1, 100))
  expect_identical(wtd$varb, ref$varb)
  expect_identical(wtd$icc, ref$icc)
})

test_that("'weights' argument is validated", {
  s <- data.frame(y = rnorm(20), psu = rep(1:4, each = 5), w = 2)
  expect_error(varcomp(y ~ psu, data = s, weights = rep(-1, 20)),
               "positive and finite")
  expect_error(varcomp(y ~ psu, data = s, weights = rep(NA_real_, 20)),
               "positive and finite")
  expect_error(varcomp(s$y, stage_id = list(s$psu), weights = c(2, 2)),
               "same length as the outcome")
  expect_error(varcomp(y ~ psu, data = s, weights = ~nope),
               "not found in 'data'")
  expect_error(varcomp(y ~ psu, data = s, weights = ~w + y),
               "exactly one variable")
})

test_that("formula interface warns on a samplyr sample without weights", {
  s <- data.frame(y = rnorm(20), psu = rep(1:4, each = 5), w = 2)
  class(s) <- c("tbl_sample", "data.frame")
  expect_warning(varcomp(y ~ psu, data = s), "population frame")
  expect_no_warning(varcomp(y ~ psu, data = s, weights = ~w))
})

test_that("a zero-mean outcome keeps icc and var_ratio, which do not use the mean", {
  # nonconstant, but centred on zero: every relvariance is infinite while the
  # split between the components is unaffected
  y <- rep(c(-2, -1, 1, 2), each = 4) + rep(c(-0.2, -0.1, 0.1, 0.2), 4)
  psu_id <- rep(seq_len(4), each = 4)
  expect_gt(var(y), 1)
  expect_equal(mean(y), 0)

  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "mean is approximately zero"
  )
  expect_true(is.infinite(res$varb))
  expect_true(is.infinite(res$unit_relvar))

  # the same components read off a shifted copy, where nothing is degenerate
  shifted <- varcomp(y + 100, stage_id = list(psu_id))
  expect_equal(res$icc, shifted$icc)
  expect_equal(res$var_ratio, shifted$var_ratio)
  expect_gt(res$icc, 0.9)
})

test_that("a mean left just off zero by rounding is still treated as zero", {
  # whether the terms of a centred outcome cancel to exactly zero depends on
  # the width of the accumulator, so the same data gives a mean of 0 on one
  # platform and of 1e-17 on another. Both are the same zero mean.
  y <- rep(c(-2, -1, 1, 2), each = 4) + rep(c(-0.2, -0.1, 0.1, 0.2), 4)
  y[1] <- y[1] + 1e-15
  psu_id <- rep(seq_len(4), each = 4)
  expect_gt(abs(mean(y)), 0)

  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "mean is approximately zero"
  )
  expect_true(is.infinite(res$varb))
  expect_true(is.infinite(res$varw))
  expect_true(is.infinite(res$unit_relvar))

  shifted <- varcomp(y + 100, stage_id = list(psu_id))
  expect_equal(res$icc, shifted$icc)
  expect_equal(res$var_ratio, shifted$var_ratio)
})

test_that("relvariances do not depend on the outcome's unit of measurement", {
  y <- rep(c(-2, -1, 1, 2), each = 4) + rep(c(-0.2, -0.1, 0.1, 0.2), 4) + 3
  psu_id <- rep(seq_len(4), each = 4)
  ref <- varcomp(y, stage_id = list(psu_id))
  tiny <- varcomp(y * 1e-9, stage_id = list(psu_id))

  expect_equal(tiny$varb, ref$varb)
  expect_equal(tiny$unit_relvar, ref$unit_relvar)
  expect_equal(tiny$icc, ref$icc)
})

test_that("the two degenerate outcomes are told apart", {
  psu_id <- rep(seq_len(4), each = 4)
  collect <- function(expr) {
    msgs <- character()
    withCallingHandlers(
      expr,
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    msgs
  }

  zero_mean <- collect(varcomp(rep(c(-1, 1), 8), stage_id = list(psu_id)))
  constant <- collect(varcomp(rep(7, 16), stage_id = list(psu_id)))

  expect_true(any(grepl("mean is approximately zero", zero_mean)))
  expect_false(any(grepl("no variance to split", zero_mean)))
  expect_true(any(grepl("no variance to split", constant)))
  expect_false(any(grepl("mean is approximately zero", constant)))
})

test_that("a constant outcome is degenerate whether or not weights are given", {
  # centring a constant outcome cancels exactly under equal weights only, so
  # the weighted variance bottoms out at rounding noise rather than at zero
  psu_id <- rep(seq_len(4), each = 3)
  for (seed in c(3L, 4L, 8L, 17L)) {
    set.seed(seed)
    w <- runif(12, 0.3, 9)
    expect_warning(
      res <- varcomp(rep(7, 12), stage_id = list(psu_id), weights = w),
      "no variance to split"
    )
    expect_identical(res$icc, 0)
    expect_identical(res$var_ratio, 1)
  }
})

## Back-out of an icc from a published design effect

test_that("varcomp inverts the two-stage clustering identity", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  expect_s3_class(vc, "svyplan_varcomp")
  expect_identical(vc$source, "deff")
  expect_identical(vc$stages, 2L)
  expect_equal(vc$icc, (1.8 - 1) / 19)
  expect_identical(vc$var_ratio, 1)
  expect_identical(vc$params$n_per_psu, 20)

  vk <- varcomp(deff = 1.8, n_per_psu = 20, var_ratio = 1.2)
  expect_equal(vk$icc, (1.8 / 1.2 - 1) / 19)
  expect_identical(vk$var_ratio, 1.2)
})

test_that("a design effect identifies icc and var_ratio and nothing else", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  expect_identical(vc$varb, NA_real_)
  expect_identical(vc$varw, NA_real_)
  expect_identical(vc$unit_relvar, NA_real_)
  expect_null(vc$strata)
})

test_that("the round trip returns the design effect at the source take", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  expect_equal(as.double(design_effect(vc, n_per_psu = 20)), 1.8)
  vk <- varcomp(deff = 2.4, n_per_psu = 15, var_ratio = 1.3)
  expect_equal(as.double(design_effect(vk, n_per_psu = 15)), 2.4)
})

test_that("re-planning at another take moves the design effect", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  smaller <- as.double(design_effect(vc, n_per_psu = 12))
  expect_lt(smaller, 1.8)
  expect_equal(smaller, 1 + vc$icc * 11)
})

test_that("design_effect does not read the source take back", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  expect_error(design_effect(vc), "supply the take you are planning for")
  expect_error(design_effect(vc), "n_per_psu = 20")
})

test_that("the identifying take is the size-weighted average", {
  takes <- c(18, 22, 20, 25, 15)
  vc <- varcomp(deff = 1.8, n_per_psu = takes)
  b_star <- sum(takes^2) / sum(takes)
  expect_equal(vc$params$n_per_psu, b_star)
  expect_equal(vc$params$n_per_psu_nominal, mean(takes))
  expect_equal(vc$icc, (1.8 - 1) / (b_star - 1))

  # b* = b_bar (1 + cv_b^2) with the population cv, and never below b_bar
  cv_b <- sqrt(mean((takes - mean(takes))^2)) / mean(takes)
  expect_equal(b_star, mean(takes) * (1 + cv_b^2))
  expect_gt(b_star, mean(takes))
})

test_that("a constant take reduces to itself", {
  vc <- varcomp(deff = 1.8, n_per_psu = rep(20, 7))
  expect_equal(vc$params$n_per_psu, 20)
  expect_equal(vc$icc, varcomp(deff = 1.8, n_per_psu = 20)$icc)
  expect_equal(varcomp(deff = 1.8, n_per_psu = 20, cv_take = 0)$icc,
               varcomp(deff = 1.8, n_per_psu = 20)$icc)
})

test_that("cv_take is the summary form of the same reduction", {
  takes <- c(18, 22, 20, 25, 15)
  cv_b <- sqrt(mean((takes - mean(takes))^2)) / mean(takes)
  expect_equal(
    varcomp(deff = 1.8, n_per_psu = mean(takes), cv_take = cv_b)$icc,
    varcomp(deff = 1.8, n_per_psu = takes)$icc
  )
})

test_that("varying takes raise the icc relative to the nominal one", {
  nominal <- varcomp(deff = 1.8, n_per_psu = 20)$icc
  weighted <- varcomp(deff = 1.8, n_per_psu = 20, cv_take = 0.3)$icc
  expect_lt(weighted, nominal)
  # the relative understatement tracks the relative gap in the take
  expect_equal((nominal - weighted) / weighted, (20 * 1.09 - 20) / (20 - 1),
               tolerance = 0.02)
})

test_that("a published standard error forms the same design effect", {
  se <- 0.021
  p <- 0.30
  n <- 1200
  deff <- se^2 / (p * (1 - p) / n)
  expect_equal(
    varcomp(se = se, p = p, n = n, n_per_psu = 20)$icc,
    varcomp(deff = deff, n_per_psu = 20)$icc
  )
  expect_equal(
    varcomp(se = 1.2, var = 900, n = 1200, n_per_psu = 20)$icc,
    varcomp(deff = 1.2^2 / (900 / 1200), n_per_psu = 20)$icc
  )
})

test_that("a design effect below var_ratio gives a negative icc and warns", {
  expect_warning(vc <- varcomp(deff = 0.8, n_per_psu = 20),
                 "beat simple random sampling")
  expect_lt(vc$icc, 0)
  expect_equal(vc$icc, (0.8 - 1) / 19)
  # the warning says where the value will be refused
  expect_error(design_effect(vc, n_per_psu = 20), "must be in \\[0, 1\\]")
  expect_error(
    n_cluster(stage_cost = c(500, 50), icc = vc, budget = 1e5,
              unit_relvar = 0.6),
    "must be in \\[0, 1\\]"
  )
})

test_that("an icc above 1 warns rather than being clamped", {
  expect_warning(vc <- varcomp(deff = 30, n_per_psu = 20),
                 "unequal weighting")
  expect_gt(vc$icc, 1)
})

test_that("the back-out validates its inputs", {
  expect_error(varcomp(deff = 1.8), "'n_per_psu' is required")
  expect_error(varcomp(deff = 1.8, n_per_psu = 1), "must exceed 1")
  expect_error(varcomp(deff = 1.8, n_per_psu = 0.5), "must exceed 1")
  expect_error(varcomp(deff = 1.8, n_per_psu = -2), "positive and finite")
  expect_error(varcomp(deff = 1.8, n_per_psu = NA_real_), "positive and finite")
  expect_error(varcomp(deff = 0, n_per_psu = 20), "'deff' must be positive")
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, cv_take = -0.1),
               "non-negative")
  expect_error(varcomp(deff = 1.8, n_per_psu = c(10, 20), cv_take = 0.2),
               "takes that were not supplied")
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, var_ratio = 0),
               "'var_ratio' must be positive")
})

test_that("a three-stage back-out is refused as under-identified", {
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, n_per_ssu = 4),
               "does not identify a three-stage design")
})

test_that("the two back-out entry points do not mix", {
  expect_error(varcomp(deff = 1.8, se = 0.02, n_per_psu = 20),
               "either 'deff' or 'se'")
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, p = 0.3),
               "'deff' is already the design effect")
  expect_error(varcomp(se = 0.021, n = 1200, n_per_psu = 20),
               "'p' for a proportion or 'var' for a mean")
  expect_error(varcomp(se = 0.021, p = 0.3, var = 4, n = 1200, n_per_psu = 20),
               "'p' for a proportion or 'var' for a mean")
  expect_error(varcomp(se = 0.021, p = 1.3, n = 1200, n_per_psu = 20),
               "'p' must be in \\(0, 1\\)")
})

test_that("back-out arguments and data do not mix", {
  set.seed(88)
  y <- rnorm(40)
  psu <- rep(1:8, each = 5)
  expect_error(varcomp(y, stage_id = list(psu), deff = 1.8, n_per_psu = 20),
               "takes no data")
  expect_error(varcomp(y, stage_id = list(psu), n_per_psu = 20),
               "belongs to the design effect back-out")
  expect_error(varcomp(y, stage_id = list(psu), cv_take = 0.2, var_ratio = 2),
               "belong to the design effect back-out")
})

test_that("an estimated varcomp keeps the data provenance", {
  set.seed(404)
  frame <- data.frame(y = rnorm(200, 50, 10), psu = rep(1:20, each = 10))
  vc <- varcomp(y ~ psu, data = frame)
  expect_identical(vc$source, "data")
  expect_identical(vc$params, list())
  expect_identical(varcomp(y ~ psu, data = frame, strata = NULL)$source, "data")
})

test_that("the unit relvariance a design effect cannot give is demanded", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  expect_error(n_cluster(stage_cost = c(500, 50), icc = vc, budget = 1e5),
               "not identified by a design effect")
  expect_error(prec_cluster(n = c(40, 20), icc = vc),
               "not identified by a design effect")

  res <- n_cluster(stage_cost = c(500, 50), icc = vc, budget = 1e5,
                   unit_relvar = 0.6)
  expect_s3_class(res, "svyplan_cluster")
  expect_equal(res$params$unit_relvar, 0.6)
  expect_equal(unname(res$params$icc), vc$icc)
})

test_that("a plan supplies the unit relvariance a design effect cannot", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20)
  bare <- svyplan(stage_cost = c(500, 50))
  expect_error(n_cluster(icc = vc, budget = 1e5, plan = bare),
               "not identified by a design effect")
  filled <- svyplan(stage_cost = c(500, 50), unit_relvar = 0.6)
  expect_equal(
    n_cluster(icc = vc, budget = 1e5, plan = filled)$cv,
    n_cluster(stage_cost = c(500, 50), icc = vc, budget = 1e5,
              unit_relvar = 0.6)$cv
  )
})

test_that("an estimated varcomp still overrides a supplied unit relvariance", {
  set.seed(77)
  frame <- data.frame(y = rnorm(200, 50, 10), psu = rep(1:20, each = 10))
  vc <- varcomp(y ~ psu, data = frame)
  expect_equal(
    n_cluster(stage_cost = c(500, 50), icc = vc, budget = 1e5)$cv,
    n_cluster(stage_cost = c(500, 50), icc = vc, budget = 1e5,
              unit_relvar = 99)$cv
  )
})

test_that("a backed-out varcomp prints and exports its provenance", {
  vc <- varcomp(deff = 1.8, n_per_psu = 20, cv_take = 0.3)
  out <- capture.output(print(vc))
  expect_match(out[1], "from a design effect")
  expect_match(out[4], "deff = 1.8000 at n_per_psu = 21.8")
  expect_match(out[4], "nominal 20")
  expect_match(out[5], "not identified by a design effect")
  expect_match(format(vc), "from deff")

  df <- as.data.frame(vc)
  expect_identical(nrow(df), 1L)
  expect_identical(df$source, "deff")
  expect_equal(df$n_per_psu, 21.8)
  expect_true(is.na(df$unit_relvar))

  plain <- capture.output(print(varcomp(deff = 1.8, n_per_psu = 20)))
  expect_false(grepl("nominal", plain[4]))
})

test_that("a backed-out icc drives the cluster functions", {
  vc <- varcomp(deff = 2.0, n_per_psu = 25)
  deff_at_10 <- as.double(design_effect(vc, n_per_psu = 10))
  expect_equal(effective_n(vc, n = 1000, n_per_psu = 10), 1000 / deff_at_10)
  prec <- prec_cluster(n = c(50, 10), icc = vc, unit_relvar = 0.5)
  expect_equal(
    prec$cv,
    sqrt(0.5 / 500 * (1 + vc$icc * 9)),
    tolerance = 1e-8
  )
})

test_that("the back-out refuses data arguments alongside a design effect", {
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, weights = rep(1, 4)),
               "without 'weights'")
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, strata = ~region),
               "without 'strata'")
  expect_error(varcomp(deff = 1.8, n_per_psu = 20, prob = rep(0.25, 4)),
               "without 'prob'")
})

test_that("two-stage probs recover the conditional weights exactly", {
  skip_if_not_installed("survey")
  set.seed(4242)
  NP <- 40; nb <- 20
  psu <- rep(seq_len(NP), each = nb)
  y <- rnorm(NP, 0, 1)[psu] + rnorm(NP * nb, 0, 3)
  keep <- sort(sample(NP, 10))
  s <- do.call(rbind, lapply(keep, function(g) {
    ix <- sample(which(psu == g), 5)
    data.frame(psu = g, y = y[ix], p1 = 10 / NP, p2 = 5 / nb,
               w_full = (NP / 10) * (nb / 5), w_within = nb / 5)
  }))

  # The stage probabilities are kept apart, so the method recovers exactly the
  # within-cluster weights the estimator documents that it needs.
  auto <- varcomp(survey::svydesign(ids = ~psu, probs = ~p1 + p2, data = s), ~y)
  ref <- varcomp(s$y, stage_id = list(s$psu), weights = s$w_within)
  expect_equal(auto$icc, ref$icc, tolerance = 1e-12)
  expect_equal(auto$varb, ref$varb, tolerance = 1e-12)

  # The full design weight is a different, wrong answer, which is what makes
  # the recovery worth doing rather than a refactor.
  wrong <- varcomp(s$y, stage_id = list(s$psu), weights = s$w_full)
  expect_false(isTRUE(all.equal(wrong$icc, ref$icc, tolerance = 1e-3)))
})

test_that("three-stage probs recover the conditional weights", {
  skip_if_not_installed("survey")
  set.seed(4243)
  s <- expand.grid(el = 1:3, ssu = 1:4, psu = 1:12)
  s$y <- rnorm(12, 0, 1)[s$psu] + rnorm(48, 0, 0.7)[with(s, (psu - 1) * 4 + ssu)] +
    rnorm(nrow(s), 0, 2)
  s$p1 <- 0.25; s$p2 <- 0.5; s$p3 <- 0.6

  auto <- varcomp(
    survey::svydesign(ids = ~psu + ssu, probs = ~p1 + p2 + p3, data = s), ~y
  )
  ref <- varcomp(s$y, stage_id = list(s$psu, s$ssu),
                 weights = rep(1 / (0.5 * 0.6), nrow(s)))
  expect_equal(auto$icc, ref$icc, tolerance = 1e-12)
  expect_equal(auto$stages, 3L)
})

test_that("a weights-only multistage design is refused, not guessed at", {
  skip_if_not_installed("survey")
  set.seed(4244)
  s <- data.frame(psu = rep(1:10, each = 5), y = rnorm(50), w = 16)
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = s)

  # A constant full weight is equally consistent with a census inside each
  # cluster and with a two-stage design that folded stage 2 into the weight.
  expect_error(varcomp(dsgn, ~y), "cannot be recovered")
  expect_error(varcomp(dsgn, ~y), "weights=")

  # Naming the conditional weights resolves it.
  expect_s3_class(varcomp(dsgn, ~y, weights = rep(4, 50)), "svyplan_varcomp")

  # Unit weights carry no stage to separate, so they still go straight through.
  unit <- survey::svydesign(ids = ~psu, weights = rep(1, 50), data = s)
  expect_equal(varcomp(unit, ~y)$icc,
               varcomp(y ~ psu, data = s)$icc, tolerance = 1e-12)
})

test_that("recovered components are near-unbiased for the frame components", {
  skip_if_not_installed("survey")
  set.seed(4245)
  NP <- 60; nb <- 30

  # The estimand is the frame's own B^2 / W^2 relvariance icc, recomputed per
  # replicate, not the superpopulation ratio sigma_b^2 / (sigma_b^2 +
  # sigma_w^2). Benchmarking against the latter measures the wrong thing and
  # makes the biased estimator look like the better one.
  one <- function() {
    psu <- rep(seq_len(NP), each = nb)
    y <- rnorm(NP, 0, 1)[psu] + rnorm(NP * nb, 0, 3)
    keep <- sort(sample(NP, 15))
    s <- do.call(rbind, lapply(keep, function(g) {
      ix <- sample(which(psu == g), 6)
      data.frame(psu = g, y = y[ix], p1 = 15 / NP, p2 = 6 / nb,
                 w_full = (NP / 15) * (nb / 6))
    }))
    c(
      frame = varcomp(y, stage_id = list(psu))$icc,
      recovered = varcomp(
        survey::svydesign(ids = ~psu, probs = ~p1 + p2, data = s), ~y
      )$icc,
      full = varcomp(s$y, stage_id = list(s$psu), weights = s$w_full)$icc
    )
  }
  est <- vapply(seq_len(60), function(i) one(), numeric(3L))
  bias_recovered <- mean(est["recovered", ] - est["frame", ])
  bias_full <- mean(est["full", ] - est["frame", ])

  # Averaged over replicates, so the test does not turn on one draw.
  expect_lt(abs(bias_recovered), 0.015)
  expect_lt(abs(bias_recovered), abs(bias_full))
})
