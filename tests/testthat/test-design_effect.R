test_that("design_effect kish formula", {
  set.seed(42)
  w <- runif(100, 1, 5)
  result <- design_effect(w, method = "kish")

  # Kish: 1 + sum((w - mean(w))^2) / n / mean(w)^2
  n <- length(w)
  expected <- 1 + sum((w - mean(w))^2) / n / mean(w)^2

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("design_effect cluster planning formula", {
  result <- design_effect(delta = 0.05, psu_size =25, method = "cluster")
  expect_equal(result, 1 + (25 - 1) * 0.05)
  expect_equal(result, 2.2)
})

test_that("design_effect cluster with svyplan_varcomp", {
  vc <- .new_svyplan_varcomp(
    varb =0.01, varw =1.0, delta = 0.05, k = 1.0,
    rel_var = 1.0, stages = 2L
  )
  result <- design_effect(delta = vc, psu_size =25, method = "cluster")
  expect_equal(result, 1 + (25 - 1) * 0.05)
})

test_that("design_effect henry method", {
  set.seed(42)
  n <- 100
  x_cal <- runif(n, 10, 100)
  y <- 5 + 0.3 * x_cal + rnorm(n, 0, 2)
  w <- runif(n, 1, 5)

  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1L)
})

test_that("design_effect spencer method", {
  set.seed(42)
  n <- 100
  p_sel <- runif(n, 0.01, 0.2)
  y <- 50 + rnorm(n, 0, 10)
  w <- 1 / p_sel

  result <- design_effect(w, y = y, prob = p_sel, method = "spencer")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1L)
})

test_that("design_effect henry returns 1 for equal weights", {
  set.seed(42)
  n <- 100
  w <- rep(5, n)
  x_cal <- runif(n, 10, 100)
  y <- 5 + 0.3 * x_cal + rnorm(n, 0, 2)
  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_equal(result, 1.0)
})

test_that("design_effect spencer returns 1 for equal weights", {
  set.seed(42)
  n <- 100
  w <- rep(5, n)
  p <- rep(0.1, n)
  y <- 50 + rnorm(n, 0, 10)
  result <- design_effect(w, y = y, prob = p, method = "spencer")
  expect_equal(result, 1.0)
})

test_that("design_effect henry returns 1 for constant y", {
  set.seed(42)
  n <- 100
  w <- runif(n, 1, 5)
  x_cal <- runif(n, 10, 100)
  y <- rep(50, n)
  result <- design_effect(w, y = y, x_cal = x_cal, method = "henry")
  expect_equal(result, 1.0)
})

test_that("design_effect spencer returns 1 for constant y", {
  set.seed(42)
  n <- 100
  w <- runif(n, 1, 5)
  p <- runif(n, 0.01, 0.2)
  y <- rep(50, n)
  result <- design_effect(w, y = y, prob = p, method = "spencer")
  expect_equal(result, 1.0)
})

test_that("design_effect spencer handles equal probabilities", {
  set.seed(42)
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
    design_effect(delta = -0.2, psu_size =10, method = "cluster"),
    "\\[0, 1\\]"
  )
  expect_error(
    design_effect(delta = 1.5, psu_size =10, method = "cluster"),
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
  expect_error(design_effect(delta = 0.05, method = "cluster"),
               "'delta' and 'psu_size' are required")
  expect_error(design_effect(psu_size = 25, method = "cluster"),
               "'delta' and 'psu_size' are required")
  w <- runif(10, 1, 5)
  expect_error(design_effect(w, method = "henry"),
               "'y' and 'x_cal' are required")
  expect_error(design_effect(w, method = "spencer"),
               "'y' and 'prob' are required")
  expect_error(design_effect(w, y = rnorm(10), method = "cr"),
               "'strata_id' or 'cluster_id'")
})

test_that("design_effect cr stratified + clustered", {
  set.seed(42)
  n <- 200
  strvar <- rep(1:5, each = 40)
  clvar <- rep(1:50, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- rep(2L, 5)

  result <- design_effect(w, y = y, strata_id = strvar, cluster_id = clvar,
                          stages = stages, method = "cr")
  expect_true(is.list(result))
  expect_true("strata" %in% names(result))
  expect_true("overall" %in% names(result))
  expect_true(is.numeric(result$overall))
  expect_equal(nrow(result$strata), 5L)
  expect_true(all(c("deff_w", "deff_c", "deff_s") %in% names(result$strata)))
  expect_equal(result$overall, sum(result$strata$deff_w *
                                     result$strata$deff_c *
                                     result$strata$deff_s))
})

test_that("design_effect cr stratified, no clusters", {
  set.seed(42)
  n <- 100
  strvar <- rep(1:5, each = 20)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, strata_id = strvar, method = "cr")
  expect_true(is.list(result))
  expect_equal(result$overall, sum(result$strata$deff_w * result$strata$deff_s))
})

test_that("design_effect cr unstratified + clustered", {
  set.seed(42)
  n <- 100
  clvar <- rep(1:25, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, cluster_id = clvar, method = "cr")
  expect_true(is.list(result))
  expect_true("rho" %in% names(result$strata))
  expect_equal(result$overall, result$strata$deff_w * result$strata$deff_c)
})

test_that("design_effect cr mixed stages", {
  set.seed(42)
  n <- 160
  strvar <- rep(1:4, each = 40)
  clvar <- rep(1:40, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- c(1L, 2L, 2L, 1L)

  result <- design_effect(w, y = y, strata_id = strvar, cluster_id = clvar,
                          stages = stages, method = "cr")
  expect_equal(result$strata$deff_c[1], 1)
  expect_equal(result$strata$deff_c[4], 1)
})

test_that("design_effect cr equal weights gives deff_w near 1", {
  set.seed(42)
  n <- 100
  clvar <- rep(1:25, each = 4)
  w <- rep(50, n)
  y <- rnorm(n, 50, 10)

  result <- design_effect(w, y = y, cluster_id = clvar, method = "cr")
  expect_equal(result$strata$deff_w, 1, tolerance = 1e-10)
})

test_that("design_effect cr agrees with survey package", {
  skip_if_not_installed("survey")
  set.seed(42)
  n <- 200
  strvar <- rep(1:4, each = 50)
  clvar <- rep(1:50, each = 4)
  w <- runif(n, 10, 50)
  y <- rnorm(n, 50, 10)
  stages <- rep(2L, 4)

  our <- design_effect(w, y = y, strata_id = strvar, cluster_id = clvar,
                       stages = stages, method = "cr")

  dsgn <- survey::svydesign(ids = ~clvar, strata = ~strvar,
                            data = data.frame(y = y), weights = w,
                            nest = TRUE)
  mn <- survey::svymean(~y, design = dsgn, deff = TRUE)
  survey_deff <- as.numeric(survey::deff(mn))

  expect_equal(our$overall, survey_deff, tolerance = 0.05)
})

test_that("design_effect cluster rejects vector delta", {
  expect_error(
    design_effect(delta = c(0.05, 0.10), psu_size =25),
    "length 1"
  )
})

test_that("design_effect cluster accepts scalar delta", {
  d <- design_effect(delta = 0.05, psu_size =25)
  expect_equal(d, 1 + 24 * 0.05)
})

test_that("design_effect cr constant y stratified + clustered returns 1", {
  set.seed(1)
  n <- 100
  w <- runif(n, 1, 5)
  y <- rep(42, n)
  strvar <- rep(1:2, each = 50)
  clvar <- rep(1:20, each = 5)
  expect_warning(
    res <- design_effect(w, y = y, strata_id = strvar, cluster_id = clvar,
                         stages = c(2L, 2L), method = "cr"),
    "variance is approximately zero"
  )
  expect_equal(res$overall, 1)
  expect_true(all(is.na(res$strata$deff_s)))
  expect_true(all(is.na(res$strata$deff_c)))
  expect_true(all(is.na(res$strata$rho_h)))
  expect_true(all(!is.na(res$strata$cv2_w)))
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
  expect_equal(res$overall, 1)
  expect_true(is.na(res$strata$rho))
  expect_true(is.na(res$strata$deff_c))
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
  expect_equal(res$overall, 1)
  expect_true(all(is.na(res$strata$deff_s)))
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
    design_effect(w, y = y, strata_id = rep(1:2, 5), cluster_id = 1:9, method = "cr"),
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

# --- Henry input validation ---

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

# --- Spencer input validation ---

test_that("spencer rejects y with NA", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, NA, 3), prob = c(0.1, 0.2, 0.3), method = "spencer"),
    "'y'.*NA"
  )
})

test_that("spencer rejects prob with Inf", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, 2, 3), prob = c(0.1, Inf, 0.3), method = "spencer"),
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
    design_effect(w, y = c(1, 2, 3), prob = c(0.1, 1.5, 0.3), method = "spencer"),
    "'prob'.*\\(0, 1\\]"
  )
})

test_that("spencer rejects length mismatch", {
  w <- c(1, 2, 3)
  expect_error(
    design_effect(w, y = c(1, 2, 3), prob = c(0.1, 0.2), method = "spencer"),
    "'prob'.*same length"
  )
})

# --- CR input validation ---

test_that("cr rejects y with NA", {
  w <- c(1, 2, 3, 4)
  expect_error(
    design_effect(w, y = c(1, NA, 3, 4), strata_id = c(1, 1, 2, 2), method = "cr"),
    "'y'.*NA"
  )
})

test_that("cr rejects y with Inf", {
  w <- c(1, 2, 3, 4)
  expect_error(
    design_effect(w, y = c(1, Inf, 3, 4), strata_id = c(1, 1, 2, 2), method = "cr"),
    "'y'.*finite"
  )
})
