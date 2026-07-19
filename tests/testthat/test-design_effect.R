test_that("design_effect kish formula", {
  set.seed(1069)
  w <- runif(100, 1, 5)
  result <- design_effect(w, method = "kish")

  # Kish: 1 + sum((w - mean(w))^2) / n / mean(w)^2
  n <- length(w)
  expected <- 1 + sum((w - mean(w))^2) / n / mean(w)^2

  expect_s3_class(result, "svyplan_design_effect")
  expect_equal(as.double(result), expected, tolerance = 1e-10)
})

test_that("design_effect cluster planning formula", {
  result <- design_effect(delta = 0.05, psu_size = 25, method = "cluster")
  expect_s3_class(result, "svyplan_design_effect")
  expect_equal(as.double(result), 1 + (25 - 1) * 0.05)
  expect_equal(as.double(result), 2.2)
})

test_that("design_effect has explicit display and export methods", {
  result <- design_effect(delta = 0.05, psu_size = 25, method = "cluster")
  df <- as.data.frame(result)

  expect_identical(names(df), c("method", "design_effect"))
  expect_equal(df$method, "cluster")
  expect_equal(df$design_effect, as.double(result))
  expect_match(format(result), "svyplan_design_effect")
  expect_output(print(result), "overall")
})

test_that("design_effect can be passed directly to sample-size functions", {
  result <- design_effect(delta = 0.05, psu_size = 25)
  direct <- n_prop(p = 0.3, moe = 0.05, deff = result)
  numeric <- n_prop(p = 0.3, moe = 0.05, deff = as.double(result))

  expect_equal(direct$n, numeric$n)
  expect_false(inherits(direct$n, "svyplan_design_effect"))
})

test_that("design_effect cluster with svyplan_varcomp", {
  vc <- .new_svyplan_varcomp(
    varb = 0.01,
    varw = 1.0,
    delta = 0.05,
    k = 1.0,
    rel_var = 1.0,
    stages = 2L
  )
  result <- design_effect(delta = vc, psu_size = 25, method = "cluster")
  expect_equal(as.double(result), 1 + (25 - 1) * 0.05)
})

test_that("design_effect henry method", {
  set.seed(1087)
  n <- 100
  x_cal <- runif(n, 10, 100)
  y <- 5 + 0.3 * x_cal + rnorm(n, 0, 2)
  w <- runif(n, 1, 5)

  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1L)
})

test_that("design_effect spencer method", {
  set.seed(1091)
  n <- 100
  p_sel <- runif(n, 0.01, 0.2)
  y <- 50 + rnorm(n, 0, 10)
  w <- 1 / p_sel

  result <- design_effect(w, y = y, prob = p_sel, method = "spencer")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1L)
})

test_that("design_effect henry returns 1 for equal weights", {
  set.seed(1093)
  n <- 100
  w <- rep(5, n)
  x_cal <- runif(n, 10, 100)
  y <- 5 + 0.3 * x_cal + rnorm(n, 0, 2)
  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_equal(as.double(result), 1.0)
})

test_that("design_effect spencer returns 1 for equal weights", {
  set.seed(1097)
  n <- 100
  w <- rep(5, n)
  p <- rep(0.1, n)
  y <- 50 + rnorm(n, 0, 10)
  result <- design_effect(w, y = y, prob = p, method = "spencer")
  expect_equal(as.double(result), 1.0)
})

test_that("design_effect henry returns 1 for constant y", {
  set.seed(1103)
  n <- 100
  w <- runif(n, 1, 5)
  x_cal <- runif(n, 10, 100)
  y <- rep(50, n)
  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_equal(as.double(result), 1.0)
})

test_that("design_effect spencer returns 1 for constant y", {
  set.seed(1109)
  n <- 100
  w <- runif(n, 1, 5)
  p <- runif(n, 0.01, 0.2)
  y <- rep(50, n)
  result <- design_effect(w, y = y, prob = p, method = "spencer")
  expect_equal(as.double(result), 1.0)
})

test_that("design_effect spencer handles equal probabilities", {
  set.seed(1117)
  n <- 100
  w <- runif(n, 1, 5)
  p <- rep(0.1, n)
  y <- 50 + rnorm(n, 0, 10)
  result <- design_effect(w, y = y, prob = p, method = "spencer")
  expect_true(is.numeric(result))
  expect_false(is.nan(result))
})

test_that("design_effect cluster rejects invalid delta", {
  expect_error(
    design_effect(delta = -0.2, psu_size = 10, method = "cluster"),
    "\\[0, 1\\]"
  )
  expect_error(
    design_effect(delta = 1.5, psu_size = 10, method = "cluster"),
    "\\[0, 1\\]"
  )
})

test_that("design_effect rejects zero weights", {
  expect_error(
    design_effect(c(1, 0, 2), method = "kish"),
    "positive"
  )
})

test_that("design_effect validates inputs", {
  expect_error(
    design_effect(delta = 0.05, method = "cluster"),
    "'delta' and 'psu_size' are required"
  )
  expect_error(
    design_effect(psu_size = 25, method = "cluster"),
    "'delta' and 'psu_size' are required"
  )
  w <- runif(10, 1, 5)
  expect_error(
    design_effect(w, method = "henry"),
    "'y' and 'x_cal' are required"
  )
  expect_error(
    design_effect(w, method = "spencer"),
    "'y' and 'prob' are required"
  )
  expect_error(
    design_effect(w, y = rnorm(10), method = "cr"),
    "'strata_id' or 'cluster_id'"
  )
})

test_that("design_effect cr stratified + clustered", {
  set.seed(1123)
  n <- 200
  strata <- rep(1:5, each = 40)
  cluster <- rep(1:50, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- rep(2L, 5)

  result <- design_effect(
    w,
    y = y,
    strata_id = strata,
    cluster_id = cluster,
    stages = stages,
    method = "cr"
  )
  components <- as.data.frame(result)
  expect_s3_class(result, "svyplan_design_effect")
  expect_true(is.numeric(result))
  expect_equal(nrow(components), 5L)
  expect_true(all(c("deff_w", "deff_c", "deff_s") %in% names(components)))
  expect_equal(
    as.double(result),
    sum(
      components$deff_w *
        components$deff_c *
        components$deff_s
    )
  )
  expect_true(all(components$overall == as.double(result)))
})

test_that("design_effect cr stratified, no clusters", {
  set.seed(1129)
  n <- 100
  strvar <- rep(1:5, each = 20)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, strata_id = strvar, method = "cr")
  components <- as.data.frame(result)
  expect_equal(
    as.double(result),
    sum(components$deff_w * components$deff_s)
  )
})

test_that("design_effect cr unstratified + clustered", {
  set.seed(1151)
  n <- 100
  clvar <- rep(1:25, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, cluster_id = clvar, method = "cr")
  components <- as.data.frame(result)
  expect_true("rho" %in% names(components))
  expect_equal(
    as.double(result),
    components$deff_w * components$deff_c
  )
})

test_that("design_effect cr mixed stages", {
  set.seed(1153)
  n <- 160
  strvar <- rep(1:4, each = 40)
  clvar <- rep(1:40, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- c(1L, 2L, 2L, 1L)

  result <- design_effect(
    w,
    y = y,
    strata_id = strvar,
    cluster_id = clvar,
    stages = stages,
    method = "cr"
  )
  components <- as.data.frame(result)
  expect_equal(components$deff_c[1], 1)
  expect_equal(components$deff_c[4], 1)
})

test_that("design_effect cr equal weights gives deff_w near 1", {
  set.seed(1163)
  n <- 100
  clvar <- rep(1:25, each = 4)
  w <- rep(50, n)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, cluster_id = clvar, method = "cr")
  expect_equal(as.data.frame(result)$deff_w, 1, tolerance = 1e-10)
})

test_that("design_effect cr agrees with survey package", {
  skip_if_not_installed("survey")
  set.seed(102)
  n <- 200
  strvar <- rep(1:4, each = 50)
  clvar <- rep(1:50, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- rep(2L, 4)

  our <- design_effect(
    w,
    y = y,
    strata_id = strvar,
    cluster_id = clvar,
    stages = stages,
    method = "cr"
  )

  dsgn <- survey::svydesign(
    ids = ~clvar,
    strata = ~strvar,
    data = data.frame(y = y),
    weights = w,
    nest = TRUE
  )
  mn <- survey::svymean(~y, design = dsgn, deff = TRUE)
  survey_deff <- as.numeric(survey::deff(mn))

  expect_equal(as.double(our), survey_deff, tolerance = 0.05)
})

test_that("design_effect cluster rejects vector delta", {
  expect_error(
    design_effect(delta = c(0.05, 0.10), psu_size = 25),
    "length 1"
  )
})

test_that("design_effect cluster accepts scalar delta", {
  d <- design_effect(delta = 0.05, psu_size = 25)
  expect_equal(as.double(d), 1 + 24 * 0.05)
})

test_that("design_effect cr constant y stratified + clustered returns 1", {
  set.seed(1)
  n <- 100
  w <- runif(n, 1, 5)
  y <- rep(42, n)
  strvar <- rep(1:2, each = 50)
  clvar <- rep(1:20, each = 5)
  expect_warning(
    res <- design_effect(
      w,
      y = y,
      strata_id = strvar,
      cluster_id = clvar,
      stages = c(2L, 2L),
      method = "cr"
    ),
    "variance is approximately zero"
  )
  components <- as.data.frame(res)
  expect_equal(as.double(res), 1)
  expect_true(all(is.na(components$deff_s)))
  expect_true(all(is.na(components$deff_c)))
  expect_true(all(is.na(components$rho)))
  expect_true(all(!is.na(components$cv2_w)))
})

test_that("design_effect cr constant y unstratified + clustered returns 1", {
  set.seed(2)
  n <- 60
  w <- runif(n, 1, 3)
  y <- rep(10, n)
  clvar <- rep(1:12, each = 5)
  expect_warning(
    res <- design_effect(w, y = y, cluster_id = clvar, method = "cr"),
    "variance is approximately zero"
  )
  components <- as.data.frame(res)
  expect_equal(as.double(res), 1)
  expect_true(is.na(components$rho))
  expect_true(is.na(components$deff_c))
})

test_that("design_effect cr constant y stratified no clusters returns 1", {
  set.seed(3)
  n <- 80
  w <- runif(n, 1, 4)
  y <- rep(7, n)
  strvar <- rep(1:4, each = 20)
  expect_warning(
    res <- design_effect(w, y = y, strata_id = strvar, method = "cr"),
    "variance is approximately zero"
  )
  components <- as.data.frame(res)
  expect_equal(as.double(res), 1)
  expect_true(all(is.na(components$deff_s)))
})

test_that(".stratum_cluster_deff returns NA for constant y in stratum", {
  w <- c(1, 2, 3, 4, 5, 6)
  y <- rep(5, 6)
  clvar <- c(1, 1, 2, 2, 3, 3)
  res <- svyplan:::.stratum_cluster_deff(w, y, clvar, sig2h = 0, nh = 6)
  expect_true(is.na(res))
})

test_that("design_effect CR rejects mismatched strata_id length", {
  w <- runif(10, 1, 2)
  y <- rnorm(10)
  expect_error(
    design_effect(w, y = y, strata_id = 1:9, method = "cr"),
    "strata_id.*same length"
  )
})

test_that("design_effect CR rejects mismatched cluster_id length", {
  w <- runif(10, 1, 2)
  y <- rnorm(10)
  expect_error(
    design_effect(
      w,
      y = y,
      strata_id = rep(1:2, 5),
      cluster_id = 1:9,
      method = "cr"
    ),
    "cluster_id.*same length"
  )
})

test_that("design_effect CR rejects mismatched y length", {
  w <- runif(10, 1, 2)
  expect_error(
    design_effect(w, y = rnorm(9), strata_id = rep(1:2, 5), method = "cr"),
    "y.*same length"
  )
})

test_that("design_effect rejects Inf weights", {
  expect_error(
    design_effect(c(1, Inf, 2), method = "kish"),
    "finite"
  )
})

test_that("design_effect rejects -Inf weights", {
  expect_error(
    design_effect(c(1, -Inf, 2), method = "kish"),
    "finite"
  )
})

test_that("henry rejects y with NA", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, NA, 3), x_cal = c(1, 2, 3), method = "henry"),
    "'y'.*NA"
  )
})

test_that("henry rejects x_cal with Inf", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, 2, 3), x_cal = c(1, Inf, 3), method = "henry"),
    "'x_cal'.*finite"
  )
})

test_that("henry rejects length mismatch", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, 2), x_cal = c(1, 2, 3), method = "henry"),
    "'y'.*same length"
  )
  expect_error(
    design_effect(w, y = c(1, 2, 3), x_cal = c(1, 2), method = "henry"),
    "'x_cal'.*same length"
  )
})

test_that("spencer rejects y with NA", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(
      w,
      y = c(1, NA, 3),
      prob = c(0.1, 0.2, 0.3),
      method = "spencer"
    ),
    "'y'.*NA"
  )
})

test_that("spencer rejects prob with Inf", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(
      w,
      y = c(1, 2, 3),
      prob = c(0.1, Inf, 0.3),
      method = "spencer"
    ),
    "'prob'.*finite"
  )
})

test_that("spencer rejects prob outside (0, 1]", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, 2, 3), prob = c(0.1, 0, 0.3), method = "spencer"),
    "'prob'.*\\(0, 1\\]"
  )
  expect_error(
    design_effect(
      w,
      y = c(1, 2, 3),
      prob = c(0.1, 1.5, 0.3),
      method = "spencer"
    ),
    "'prob'.*\\(0, 1\\]"
  )
})

test_that("henry rejects single-observation input cleanly", {
  expect_error(
    design_effect(c(1), y = c(1), x_cal = c(1), method = "henry"),
    "length >= 2"
  )
})

test_that("spencer rejects single-observation input cleanly", {
  expect_error(
    design_effect(c(1), y = c(1), prob = c(1), method = "spencer"),
    "length >= 2"
  )
})

test_that("spencer rejects length mismatch", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, 2, 3), prob = c(0.1, 0.2), method = "spencer"),
    "'prob'.*same length"
  )
})

test_that("cr rejects y with NA", {
  w <- c(1, 2, 3, 4)
  expect_error(
    design_effect(
      w,
      y = c(1, NA, 3, 4),
      strata_id = c(1, 1, 2, 2),
      method = "cr"
    ),
    "'y'.*NA"
  )
})

test_that("cr rejects y with Inf", {
  w <- c(1, 2, 3, 4)
  expect_error(
    design_effect(
      w,
      y = c(1, Inf, 3, 4),
      strata_id = c(1, 1, 2, 2),
      method = "cr"
    ),
    "'y'.*finite"
  )
})

test_that("CR rejects weights below the population scale", {
  set.seed(1181)
  y <- rnorm(100)
  cl <- rep(1:25, each = 4)
  expect_error(
    design_effect(rep(1 / 100, 100), y = y, cluster_id = cl, method = "cr"),
    "population scale"
  )
})

test_that("CR rejects data with no within-cluster replication", {
  set.seed(1)
  expect_error(
    design_effect(rep(100, 10), y = rnorm(10), cluster_id = 1:10,
                  method = "cr"),
    "single observation"
  )
})

test_that("CR rejects weights at or below the sample-size scale", {
  set.seed(1)
  n <- 20
  y <- rnorm(n)
  cl <- rep(1:5, each = 4)
  expect_error(design_effect(rep(1, n), y = y, cluster_id = cl, method = "cr"),
               "must exceed the sample size")
  expect_error(design_effect(rep(0.9, n), y = y, cluster_id = cl,
                             method = "cr"),
               "must exceed the sample size")
  res <- design_effect(rep(1.5, n), y = y, cluster_id = cl, method = "cr")
  expect_true(is.finite(as.double(res)) && as.double(res) >= 0)
})

test_that("CR names an invalidly scaled stratum", {
  set.seed(2)
  y <- rnorm(20)
  stratum <- rep(c("bad", "good"), each = 10)
  w <- c(rep(0.05, 10), rep(10, 10))
  expect_error(design_effect(w, y = y, strata_id = stratum, method = "cr"),
               "stratum 'bad'")
  w_ok <- c(rep(5, 10), rep(10, 10))
  res <- design_effect(w_ok, y = y, strata_id = stratum, method = "cr")
  components <- as.data.frame(res)
  expect_true(all(is.finite(components$deff_s)))
  expect_true(all(components$deff_s >= 0))
  expect_true(is.finite(as.double(res)) && as.double(res) >= 0)
})
