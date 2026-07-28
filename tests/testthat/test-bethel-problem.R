test_that("problem builder creates correct proportion coefficients", {
  z <- .bethel_fixture()
  p <- .build_bethel_problem(z$frame, z$measures, z$targets)
  k <- match("vaccination@.overall:cv", p$constraint_ids)
  variance <- z$measures$p[1:4] * (1 - z$measures$p[1:4])

  expect_equal(p$A[, k], z$frame$N^2 * variance)
  expect_equal(p$B[k], -sum(z$frame$N * variance))
  expect_true(all(p$membership[, k]))
})

test_that("problem builder supports overlapping domain classifications", {
  z <- .bethel_fixture()
  p <- .build_bethel_problem(z$frame, z$measures, z$targets)

  north <- match("vaccination@region=North:cv", p$constraint_ids)
  urban <- match("income@residence=Urban:moe", p$constraint_ids)
  expect_identical(which(p$membership[, north]), c(1L, 2L))
  expect_identical(which(p$membership[, urban]), c(1L, 3L))
})

test_that("proportion shorthand equals explicit mean and variance", {
  z <- .bethel_fixture()
  p1 <- .build_bethel_problem(z$frame, z$measures, z$targets)
  m2 <- z$measures
  prop <- !is.na(m2$p)
  m2$mean[prop] <- m2$p[prop]
  m2$var <- ifelse(prop, m2$p * (1 - m2$p), NA)
  m2$p[prop] <- NA
  p2 <- .build_bethel_problem(z$frame, m2, z$targets)

  expect_equal(p1$A, p2$A)
  expect_equal(p1$B, p2$B)
  expect_equal(p1$total, p2$total)
})

test_that("row-specific deff and response rates override scalar defaults", {
  z <- .bethel_fixture()
  z$measures$deff <- NA_real_
  z$measures$resp_rate <- NA_real_
  z$measures$deff[1] <- 2
  z$measures$resp_rate[1] <- 0.5
  p <- .build_bethel_problem(
    z$frame, z$measures, z$targets,
    deff = 1.5, resp_rate = 0.8
  )
  k <- match("vaccination@.overall:cv", p$constraint_ids)
  variance <- 0.5 * 0.5
  expect_equal(p$A[1, k], 1000^2 * variance * 2 / 0.5)
  expect_equal(p$A[2, k], 2000^2 * 0.4 * 0.6 * 1.5 / 0.8)
})

test_that("builder rejects incomplete and ambiguous public tables", {
  z <- .bethel_fixture()
  expect_error(
    .build_bethel_problem(z$frame[-1, ], z$measures, z$targets),
    "not found in frame"
  )
  expect_error(
    .build_bethel_problem(
      z$frame,
      rbind(z$measures, z$measures[1, ]),
      z$targets
    ),
    "unique stratum x name"
  )
  expect_error(
    .build_bethel_problem(
      z$frame, z$measures,
      rbind(z$targets, z$targets[1, ])
    ),
    "duplicate"
  )
  missing_measure <- z$measures[-1, ]
  expect_error(
    .build_bethel_problem(z$frame, missing_measure, z$targets),
    "missing measures"
  )
})

test_that("unused measure values cannot change the allocation problem", {
  z <- .bethel_fixture()
  baseline <- .build_bethel_problem(z$frame, z$measures, z$targets)
  unused <- z$measures$name == "income" &
    z$measures$stratum %in% c("NR", "SR")
  augmented <- z$measures
  augmented$mean[unused] <- Inf
  augmented$sd[unused] <- NA_real_
  augmented$icc_psu <- NA_real_
  augmented$icc_psu[unused] <- 0.5

  got <- .build_bethel_problem(
    z$frame,
    augmented[nrow(augmented):1, ],
    z$targets
  )

  expect_identical(got$stages, 1L)
  expect_equal(got$A, baseline$A)
  expect_equal(got$B, baseline$B)
  expect_equal(got$bound, baseline$bound)
  expect_equal(got$cost, baseline$cost)
  expect_false(any(got$measures$stratum %in% c("NR", "SR") &
                     got$measures$name == "income"))
})

test_that("unused multistage rows need no stage planning values", {
  for (stages in 2:3) {
    z <- .bethel_multistage_fixture(stages)
    baseline <- .build_bethel_problem(z$frame, z$measures, z$targets)
    unused <- z$measures$name == "income" &
      z$measures$stratum %in% c("NR", "SR")
    augmented <- z$measures
    augmented$mean[unused] <- Inf
    augmented$sd[unused] <- NA_real_
    augmented$icc_psu[unused] <- NA_real_
    if (stages == 3L) augmented$icc_ssu[unused] <- NA_real_

    got <- .build_bethel_problem(z$frame, augmented, z$targets)

    expect_identical(got$stages, stages)
    expect_equal(got$A, baseline$A)
    expect_equal(got$B, baseline$B)
    expect_equal(got$cost, baseline$cost)
  }
})

test_that("shared evaluator reproduces legacy one-indicator precision", {
  frame <- data.frame(
    stratum = c("A", "B", "C"),
    N = c(1000, 1500, 2000)
  )
  measures <- data.frame(
    stratum = frame$stratum,
    name = "income",
    mean = c(50, 60, 55),
    sd = c(10, 15, 8)
  )
  targets <- data.frame(name = "income", cv = 0.1)
  problem <- .build_bethel_problem(frame, measures, targets)
  n <- c(100, 150, 200)
  got <- .precision_from_allocation(problem, n)
  legacy <- .alloc_metrics(
    frame$N, measures$sd, measures$mean, n,
    alpha = 0.05, deff = 1, resp_rate = 1, cost_h = rep(1, 3)
  )

  expect_equal(got$cv, legacy$cv)
  expect_equal(got$moe, legacy$moe)
  expect_equal(got$se, legacy$se)
})

test_that("CV and MOE bounds follow their public definitions", {
  z <- .bethel_fixture()
  p <- .build_bethel_problem(z$frame, z$measures, z$targets)

  cv_k <- match("vaccination@.overall:cv", p$constraint_ids)
  cv_total <- sum(z$frame$N * z$measures$p[1:4])
  expect_equal(p$bound[cv_k], (0.05 * cv_total)^2 - p$B[cv_k])

  moe_k <- match("income@residence=Urban:moe", p$constraint_ids)
  domain_N <- sum(z$frame$N[c(1, 3)])
  z_alpha <- qnorm(0.975)
  expect_equal(p$bound[moe_k], (domain_N * 2 / z_alpha)^2 - p$B[moe_k])
})

test_that("builder applies costs, allocation bounds, and stable identifiers", {
  z <- .bethel_fixture()
  z$frame$max_weight <- c(100, 50, NA, NA)
  z$frame$take_all <- c(FALSE, FALSE, TRUE, FALSE)
  p <- .build_bethel_problem(
    z$frame, z$measures, z$targets,
    unit_cost = c(2, 3, 4, 5), min_n_stratum = 3
  )
  reversed <- .build_bethel_problem(
    z$frame, z$measures, z$targets[nrow(z$targets):1, ],
    unit_cost = c(2, 3, 4, 5), min_n_stratum = 3
  )

  expect_equal(p$cost, c(2, 3, 4, 5))
  expect_equal(p$lower[1:2], c(10, 40))
  expect_equal(p$lower[3], z$frame$N[3])
  expect_equal(p$upper[3], z$frame$N[3])
  expect_setequal(p$constraint_ids, reversed$constraint_ids)
})

test_that("builder rejects invalid 0/1 flags and undefined constraints", {
  z <- .bethel_fixture()
  bad_take <- z$frame
  bad_take$take_all <- c(0, 0, 2, 0)
  expect_error(
    .build_bethel_problem(bad_take, z$measures, z$targets),
    "logical \\(or 0/1\\)"
  )

  zero <- z$measures
  zero$p[1:4] <- 0
  expect_error(
    .build_bethel_problem(z$frame, zero, z$targets),
    "no positive variance"
  )

  negligible <- z$measures
  negligible$p[1:4] <- c(0, 0, 0, 1)
  targets <- z$targets[1, , drop = FALSE]
  targets$name <- "signed"
  signed <- data.frame(
    stratum = z$frame$stratum,
    name = "signed",
    mean = c(1, -0.5, 1, -1.5),
    sd = 1
  )
  expect_error(
    .build_bethel_problem(z$frame, signed, targets),
    "zero or negligible total"
  )
})

test_that("materially negative evaluated variance is an internal error", {
  z <- .bethel_fixture()
  p <- .build_bethel_problem(z$frame, z$measures, z$targets[1, ])
  p$upper[] <- 1e9
  expect_error(
    .precision_from_allocation(p, rep(1e9, nrow(z$frame))),
    "materially negative variance"
  )
})

test_that("a negligible total is judged against the terms that formed it", {
  # rescaling every mean and sd leaves the CV constraint unchanged, so the
  # allocation must not depend on the unit the measure is expressed in
  z <- .bethel_fixture()
  z$targets <- data.frame(
    name = "income", domain = ".overall", level = NA_character_,
    cv = 0.05, moe = NA_real_
  )
  ref <- n_alloc(z$frame, measures = z$measures, targets = z$targets)
  for (scale in c(1e-6, 1e-14, 1e-30)) {
    small <- z
    small$measures$mean <- small$measures$mean * scale
    small$measures$sd <- small$measures$sd * scale
    fit <- n_alloc(small$frame, measures = small$measures,
                   targets = small$targets)
    expect_equal(fit$n, ref$n)
  }

  # a total that is zero because its terms cancel is still rejected
  zero <- z
  zero$measures$mean[5:8] <- c(50, -25, 60, -90)
  expect_error(
    n_alloc(zero$frame, measures = zero$measures, targets = zero$targets),
    "negligible total"
  )
})
