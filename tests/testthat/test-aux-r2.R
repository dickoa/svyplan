test_that("design_effect() supports optional aux_r2 calibration factor", {
  # Standalone aux_r2
  d0 <- design_effect(aux_r2 = 0)
  expect_equal(as.double(d0), 1)
  expect_equal(d0$deff_calibration, 1)

  d <- design_effect(aux_r2 = 0.25)
  expect_equal(as.double(d), 0.75)
  expect_equal(d$deff_calibration, 0.75)
  expect_no_error(capture.output(print(summary(d))))
  expect_match(summary(d)$basis, "calibration")

  # Combined with clustering and weighting
  w <- rep(c(1, 3), c(300, 100))
  d_base <- design_effect(icc = 0.05, take = 20, weights = w)
  d_cal <- design_effect(icc = 0.05, take = 20, weights = w, aux_r2 = 0.3)
  expect_equal(as.double(d_cal), as.double(d_base) * 0.7, tolerance = 1e-12)
  expect_equal(d_cal$deff_calibration, 0.7)
  expect_equal(as.data.frame(d_cal)$deff_calibration, 0.7)

  # Composes on svyplan_cluster, svyplan_varcomp, svyplan_alloc, svyplan_prec
  cl <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  expect_equal(
    as.double(design_effect(cl, aux_r2 = 0.2)),
    as.double(design_effect(cl)) * 0.8,
    tolerance = 1e-12
  )
  pr_cl <- prec_cluster(cl)
  expect_equal(
    as.double(design_effect(pr_cl, aux_r2 = 0.2)),
    as.double(design_effect(pr_cl)) * 0.8,
    tolerance = 1e-12
  )
  vc <- varcomp(deff = 2.0, take = 21)
  expect_equal(
    as.double(design_effect(vc, take = 15, aux_r2 = 0.4)),
    as.double(design_effect(vc, take = 15)) * 0.6,
    tolerance = 1e-12
  )
  fr <- data.frame(N = c(1000, 2000), sd = c(10, 20), mean = c(50, 70))
  al <- n_alloc(fr, n = 300, alloc = "neyman")
  d_al <- design_effect(al, aux_r2 = 0.25)
  expect_equal(as.double(d_al), as.double(design_effect(al)) * 0.75, tolerance = 1e-12)
  expect_match(summary(d_al)$basis, "calibration")
  d_al_w <- design_effect(al, weights = w, aux_r2 = 0.25)
  expect_equal(
    as.double(d_al_w),
    as.double(design_effect(al, weights = w)) * 0.75,
    tolerance = 1e-12
  )
  expect_match(summary(d_al_w)$basis, "weighting and calibration")
})

test_that("design_effect() validates aux_r2 bounds cleanly", {
  expect_error(design_effect(aux_r2 = -0.1), "aux_r2")
  expect_error(design_effect(aux_r2 = 1), "aux_r2")
  expect_error(design_effect(aux_r2 = 1.2), "aux_r2")
  expect_error(design_effect(aux_r2 = NA_real_), "aux_r2")
  expect_error(design_effect(aux_r2 = c(0.1, 0.2)), "aux_r2")
  expect_error(design_effect(aux_r2 = "0.2"), "aux_r2")
})
