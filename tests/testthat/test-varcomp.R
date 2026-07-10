test_that("varcomp 2-stage SRS formula interface works", {
  set.seed(42)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  result <- varcomp(income ~ district, data = frame)
  expect_s3_class(result, "svyplan_varcomp")
  expect_equal(result$stages, 2L)
  expect_true(result$delta >= 0 && result$delta <= 1)
  expect_true(result$k > 0)
  expect_true(result$rel_var > 0)
})

test_that("varcomp 2-stage SRS vector interface matches formula", {
  set.seed(42)
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
  expect_equal(result_vector$delta, result_formula$delta, tolerance = 1e-10)
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
  expect_equal(result$delta, B2 / (B2 + W2), tolerance = 1e-10)
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
  set.seed(42)
  frame <- data.frame(
    income = rnorm(200, 50000, 10000),
    district = rep(1:20, each = 10)
  )
  vc <- varcomp(income ~ district, data = frame)
  plan <- n_cluster(stage_cost = c(500, 50), delta = vc, budget = 100000)
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
  expect_false(is.na(result$delta))
})

test_that("varcomp.survey.design matches formula interface", {
  skip_if_not_installed("survey")
  set.seed(42)
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
  expect_equal(result$delta, ref$delta, tolerance = 1e-10)
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
  expect_true(result$delta >= 0 && result$delta <= 1)
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
  plan <- n_cluster(stage_cost = c(500, 50), delta = vc, budget = 100000)
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

  expect_equal(res_non$delta, res_nested$delta, tolerance = 1e-10)
  expect_equal(res_non$varb, res_nested$varb, tolerance = 1e-10)
  expect_equal(res_non$varw, res_nested$varw, tolerance = 1e-10)
})

test_that("all-singleton 2-stage returns delta = 1 with warning", {
  psu_id <- 1:10
  y <- rnorm(10, 100, 20)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "all clusters are singletons"
  )
  expect_equal(res$delta, 1)
  expect_false(is.nan(res$delta))
  expect_true(is.finite(res$k))
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

test_that("varcomp 2-stage SRS constant y gives delta = 0 with warning", {
  y <- rep(42, 50)
  psu_id <- rep(1:10, each = 5)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id)),
    "approximately zero"
  )
  expect_equal(res$delta, 0)
  expect_equal(res$k, 1)
  expect_true(is.finite(res$delta))
  expect_true(is.finite(res$k))
})

test_that("varcomp 2-stage PPS constant y gives delta = 0 with warning", {
  y <- rep(42, 50)
  psu_id <- rep(1:10, each = 5)
  pp <- rep(0.1, 10)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id), prob = pp),
    "approximately zero"
  )
  expect_equal(res$delta, 0)
  expect_equal(res$k, 1)
  expect_true(is.finite(res$delta))
  expect_true(is.finite(res$k))
})

test_that("varcomp 3-stage PPS constant y gives finite delta with warning", {
  y <- rep(42, 60)
  psu_id <- rep(1:6, each = 10)
  ssu_id <- rep(rep(1:2, each = 5), 6)
  pp <- rep(1 / 6, 6)
  expect_warning(
    res <- varcomp(y, stage_id = list(psu_id, ssu_id), prob = pp),
    "approximately zero"
  )
  expect_length(res$delta, 2)
  expect_true(all(is.finite(res$delta)))
  expect_true(all(is.finite(res$k)))
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
  expect_equal(res_slash$delta, res_in$delta, tolerance = 1e-10)
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
  expect_length(vc$delta, 2L)
  expect_length(vc$k, 2L)
  expect_true(all(vc$delta >= 0 & vc$delta <= 1))
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
  expect_equal(vc_srs$delta, vc_pps$delta, tolerance = 1e-10)
  expect_equal(vc_srs$k, vc_pps$k, tolerance = 1e-10)
  expect_equal(vc_srs$rel_var, vc_pps$rel_var, tolerance = 1e-10)
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
  expect_equal(vc_vec$delta, vc_frm$delta, tolerance = 1e-10)
  expect_equal(vc_vec$k, vc_frm$k, tolerance = 1e-10)
})

test_that("survey.design method with unit weights matches frame", {
  skip_if_not_installed("survey")
  set.seed(42)
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
  expect_equal(result$delta, ref$delta, tolerance = 1e-12)
})

test_that("weighted correction shrinks the between component", {
  skip_if_not_installed("survey")
  set.seed(42)
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
  result <- varcomp(dsgn, ~income)
  expect_lt(result$varb, ref$varb)
  expect_equal(result$varw, ref$varw, tolerance = 1e-12)
})

test_that("survey.design method uses design weights", {
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
  wtd <- varcomp(dsgn, ~y)
  unw <- varcomp(s$y, stage_id = list(s$psu))

  expect_lt(abs(wtd$delta - truth$delta), 0.1)
  expect_lt(abs(wtd$delta - truth$delta), abs(unw$delta - truth$delta))
})

test_that("weighted 2-stage PPS recovers the population delta", {
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
  wtd <- varcomp(dsgn, ~y, prob = pp)
  unw <- varcomp(s$y, stage_id = list(s$psu), prob = pp)

  expect_lt(abs(wtd$delta - truth$delta), 0.1)
  expect_lt(abs(wtd$delta - truth$delta), abs(unw$delta - truth$delta))
})

test_that("per-stratum weighted results do not depend on other strata", {
  skip_if_not_installed("survey")
  set.seed(42)
  d <- data.frame(
    y = rnorm(400, rep(c(10, 20), each = 200)),
    psu = rep(1:40, each = 10),
    region = rep(c("N", "S"), each = 200),
    w = rep(c(2, 8), each = 200)
  )
  dsgn <- survey::svydesign(ids = ~psu, weights = ~w, data = d)
  vc <- varcomp(dsgn, ~y, strata = ~region)

  for (s in c("N", "S")) {
    sub <- d[d$region == s, ]
    dsub <- survey::svydesign(ids = ~psu, weights = ~w, data = sub)
    ref <- varcomp(dsub, ~y)
    i <- match(s, vc$strata$stratum)
    expect_equal(vc$strata$delta_psu[i], ref$delta)
    expect_equal(vc$strata$k_psu[i], ref$k)
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
  expect_equal(result$delta, ref$delta, tolerance = 1e-12)
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
  res <- varcomp(dsgn, ~y)
  expect_true(all(is.finite(res$delta)))
  expect_true(all(res$delta >= 0 & res$delta <= 1))
  expect_true(all(res$k > 0))
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
  expect_null(vc$delta)
  expect_equal(nrow(vc$strata), 2L)
  expect_equal(vc$strata$stratum, c("N", "S"))

  for (s in c("N", "S")) {
    sub <- d[d$region == s, ]
    ref <- varcomp(sub$y, stage_id = list(sub$ea))
    i <- match(s, vc$strata$stratum)
    expect_equal(vc$strata$delta_psu[i], ref$delta)
    expect_equal(vc$strata$k_psu[i], ref$k)
    expect_equal(vc$strata$rel_var[i], ref$rel_var)
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
  expect_true(all(c("delta_psu", "delta_ssu", "k_psu", "k_ssu",
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
  expect_error(n_cluster(stage_cost = c(500, 50), delta = vc, cv = 0.05),
               "stratified varcomp")
  expect_error(design_effect(delta = vc, psu_size = 10, method = "cluster"),
               "stratified varcomp")
  expect_error(prec_cluster(n = c(20, 10), stage_cost = c(500, 50), delta = vc),
               "stratified varcomp")
  pooled <- varcomp(y ~ ea, data = d)
  expect_error(as.data.frame(pooled), "stratified")
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
  expect_equal(a$delta, b$delta)
  expect_error(
    varcomp(1:6, stage_id = list(psu), prob = c("1" = 0.2, "2" = 0.3, "9" = 0.5)),
    "must match the PSU identifiers"
  )
})
