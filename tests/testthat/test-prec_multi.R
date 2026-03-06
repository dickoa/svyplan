test_that("prec_multi simple mode computes per-indicator precision", {
  targets <- data.frame(
    name = c("stunting", "vaccination", "anemia"),
    p    = c(0.30, 0.70, 0.10),
    n    = c(400, 400, 400)
  )
  result <- prec_multi(targets)
  expect_s3_class(result, "svyplan_prec")
  expect_equal(result$type, "multi")
  expect_equal(nrow(result$detail), 3L)

  for (i in 1:3) {
    p <- targets$p[i]
    se_exp <- sqrt(p * (1 - p) / 400)
    expect_equal(result$detail$.se[i], se_exp, tolerance = 1e-6)
  }
})

test_that("prec_multi simple mode with deff and resp_rate", {
  targets <- data.frame(
    name      = c("stunting", "vaccination"),
    p         = c(0.30, 0.70),
    n         = c(500, 500),
    deff      = c(1.5, 2.0),
    resp_rate = c(0.8, 0.9)
  )
  result <- prec_multi(targets)
  expect_equal(nrow(result$detail), 2L)

  n_eff <- 500 * 0.8 / 1.5
  se_exp <- sqrt(0.3 * 0.7 / n_eff)
  expect_equal(result$detail$.se[1], se_exp, tolerance = 1e-6)
})

test_that("prec_multi simple mode honors global prop_method", {
  targets <- data.frame(p = 0.05, n = 400)

  res_wilson <- prec_multi(targets, prop_method = "wilson")
  ref_wilson <- prec_prop(p = 0.05, n = 400, method = "wilson")
  expect_equal(res_wilson$detail$.se[1], ref_wilson$se, tolerance = 1e-6)
  expect_equal(res_wilson$detail$.moe[1], ref_wilson$moe, tolerance = 1e-6)
  expect_equal(res_wilson$detail$.cv[1], ref_wilson$cv, tolerance = 1e-6)

  res_logodds <- prec_multi(targets, prop_method = "logodds")
  ref_logodds <- prec_prop(p = 0.05, n = 400, method = "logodds")
  expect_equal(res_logodds$detail$.se[1], ref_logodds$se, tolerance = 1e-6)
  expect_equal(res_logodds$detail$.moe[1], ref_logodds$moe, tolerance = 1e-6)
  expect_equal(res_logodds$detail$.cv[1], ref_logodds$cv, tolerance = 1e-6)
})

test_that("prec_multi targets prop_method overrides the global default", {
  targets <- data.frame(
    name = c("rare", "common"),
    p = c(0.05, 0.5),
    n = c(400, 400),
    prop_method = c("wilson", NA)
  )

  res <- prec_multi(targets, prop_method = "wald")
  ref_rare <- prec_prop(p = 0.05, n = 400, method = "wilson")
  ref_common <- prec_prop(p = 0.5, n = 400, method = "wald")

  expect_equal(res$detail$.moe[1], ref_rare$moe, tolerance = 1e-6)
  expect_equal(res$detail$.moe[2], ref_common$moe, tolerance = 1e-6)
})

test_that("prec_multi prop_method is ignored for mean rows", {
  targets <- data.frame(
    name = c("prop_ind", "mean_ind"),
    p = c(0.05, NA),
    var = c(NA, 100),
    mu = c(NA, 50),
    n = c(400, 400),
    prop_method = c("wilson", "logodds")
  )

  res <- prec_multi(targets)
  ref_prop <- prec_prop(p = 0.05, n = 400, method = "wilson")
  ref_mean <- prec_mean(var = 100, mu = 50, n = 400)

  expect_equal(res$detail$.moe[1], ref_prop$moe, tolerance = 1e-6)
  expect_equal(res$detail$.moe[2], ref_mean$moe, tolerance = 1e-6)
  expect_equal(res$detail$.cv[2], ref_mean$cv, tolerance = 1e-6)
})

test_that("prec_multi validates prop_method values in simple mode", {
  expect_error(
    prec_multi(data.frame(p = 0.3, n = 100), prop_method = "bad"),
    "'prop_method' must be one of"
  )
  expect_error(
    prec_multi(data.frame(p = 0.3, n = 100, prop_method = "bad")),
    "'prop_method' values must be one of"
  )
})

test_that("prec_multi multistage mode computes per-indicator CV", {
  targets <- data.frame(
    name   = c("stunting", "anemia"),
    p      = c(0.30, 0.10),
    n      = c(50, 50),
    psu_size     = c(12, 12),
    delta_psu = c(0.02, 0.05)
  )
  result <- prec_multi(targets, stage_cost = c(500, 50))
  expect_equal(result$type, "multi")
  expect_equal(nrow(result$detail), 2L)
  expect_true(all(!is.na(result$detail$.cv)))
})

test_that("prec_multi requires n column", {
  targets <- data.frame(p = 0.3, moe = 0.05)
  expect_error(prec_multi(targets), "must contain an 'n' column")
})

test_that("prec_multi.svyplan_n preserves proportion method from n_multi()", {
  targets <- data.frame(p = 0.05, moe = 0.02)
  n_res <- n_multi(targets, prop_method = "wilson")
  p_res <- prec_multi(n_res)
  ref <- prec_prop(p = 0.05, n = n_res$n, method = "wilson")

  expect_equal(p_res$detail$.moe[1], ref$moe, tolerance = 1e-6)
  expect_equal(p_res$detail$.cv[1], ref$cv, tolerance = 1e-6)
})

test_that("prec_multi simple rejects invalid p", {
  expect_error(
    prec_multi(data.frame(p = 0, n = 100)),
    "p.*\\(0, 1\\)"
  )
  expect_error(
    prec_multi(data.frame(p = 1, n = 100)),
    "p.*\\(0, 1\\)"
  )
  expect_error(
    prec_multi(data.frame(p = -0.1, n = 100)),
    "p.*\\(0, 1\\)"
  )
})

test_that("prec_multi simple rejects invalid var", {
  expect_error(
    prec_multi(data.frame(var = -10, n = 100)),
    "var.*positive"
  )
  expect_error(
    prec_multi(data.frame(var = 0, n = 100)),
    "var.*positive"
  )
})

test_that("prec_multi simple rejects non-positive mu", {
  expect_error(
    prec_multi(data.frame(var = 10, mu = 0, n = 100)),
    "mu.*positive"
  )
  expect_error(
    prec_multi(data.frame(var = 10, mu = -5, n = 100)),
    "mu.*positive"
  )
})

test_that("prec_multi cluster requires psu_size column", {
  targets <- data.frame(
    p = 0.3, n = 50, delta_psu = 0.05
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "psu_size.*required"
  )
})

test_that("prec_multi 3-stage requires ssu_size column", {
  targets <- data.frame(
    p = 0.3, n = 50, psu_size = 10, delta_psu = 0.05
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 100, 50)),
    "ssu_size.*required"
  )
})

test_that("prec_multi cluster rejects NA in psu_size", {
  targets <- data.frame(
    p = c(0.3, 0.1), n = c(50, 50), psu_size = c(10, NA), delta_psu = c(0.05, 0.02)
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "psu_size.*required"
  )
})

test_that("prec_multi cluster rejects negative rel_var", {
  targets <- data.frame(
    p = 0.3, n = 50, psu_size = 10, delta_psu = 0.05, rel_var = -1
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "rel_var.*positive"
  )
})

test_that("prec_multi cluster rejects zero rel_var", {
  targets <- data.frame(
    p = 0.3, n = 50, psu_size = 10, delta_psu = 0.05, rel_var = 0
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "rel_var.*positive"
  )
})

test_that("prec_multi cluster rejects negative k_psu", {
  targets <- data.frame(
    p = 0.3, n = 50, psu_size = 10, delta_psu = 0.05, k_psu = -1
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "k_psu.*positive"
  )
})

test_that("prec_multi cluster rejects negative k_ssu", {
  targets <- data.frame(
    p = 0.3, n = 50, psu_size = 10, delta_psu = 0.05, k_ssu = -1
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "k_ssu.*positive"
  )
})

test_that("prec_multi cluster rejects missing delta_psu", {
  targets <- data.frame(p = 0.3, n = 50, psu_size = 10)
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "delta_psu.*required"
  )
})

test_that("prec_multi cluster rejects delta_psu out of range", {
  targets <- data.frame(p = 0.3, n = 50, psu_size = 2, delta_psu = -2)
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "delta_psu.*\\[0, 1\\]"
  )
  targets2 <- data.frame(p = 0.3, n = 50, psu_size = 2, delta_psu = 1.5)
  expect_error(
    prec_multi(targets2, stage_cost = c(500, 50)),
    "delta_psu.*\\[0, 1\\]"
  )
})

test_that("prec_multi cluster rejects NA delta_psu", {
  targets <- data.frame(p = 0.3, n = 50, psu_size = 10, delta_psu = NA_real_)
  expect_error(
    prec_multi(targets, stage_cost = c(500, 50)),
    "delta_psu.*NA"
  )
})

test_that("prec_multi 3-stage rejects missing delta_ssu", {
  targets <- data.frame(p = 0.3, n = 50, psu_size = 10, ssu_size = 5, delta_psu = 0.05)
  expect_error(
    prec_multi(targets, stage_cost = c(500, 100, 50)),
    "delta_ssu.*required"
  )
})

test_that("prec_multi 3-stage rejects NA delta_ssu", {
  targets <- data.frame(
    p = 0.3, n = 50, psu_size = 10, ssu_size = 5, delta_psu = 0.05, delta_ssu = NA_real_
  )
  expect_error(
    prec_multi(targets, stage_cost = c(500, 100, 50)),
    "delta_ssu.*NA"
  )
})
