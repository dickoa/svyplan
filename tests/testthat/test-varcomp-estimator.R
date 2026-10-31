## Components of the ratio estimator beside the expansion estimator

.unequal_frame <- function(seed = 11, A = 40) {
  set.seed(seed)
  Ni <- sample(8:79, A, replace = TRUE)
  pop <- lapply(Ni, function(n) rnorm(n, mean = 0.4 * n, sd = 5))
  data.frame(y = unlist(pop), psu = rep(seq_len(A), Ni))
}

## The with-replacement relvariance of each estimator at a common take
.closed_form <- function(fr, a, m, ratio) {
  pop <- split(fr$y, fr$psu)
  A <- length(pop)
  Ni <- lengths(pop)
  ti <- vapply(pop, sum, 0)
  tU <- sum(ti)
  S2 <- vapply(pop, var, 0)
  spread <- if (ratio) sum((ti - Ni * tU / sum(Ni))^2) / (A - 1) else var(ti)
  (A^2 * spread / a + A * sum(Ni^2 * S2) / (a * m)) / tU^2
}

test_that("each estimator's components reproduce its closed form", {
  fr <- .unequal_frame()
  expansion <- varcomp(y ~ psu, data = fr)
  ratio <- varcomp(y ~ psu, data = fr, estimator = "ratio")
  expect_identical(expansion$estimator, "expansion")
  expect_identical(ratio$estimator, "ratio")
  for (design in list(c(8, 5), c(15, 3))) {
    expect_equal(prec_cluster(n = design, icc = expansion)$cv^2,
                 .closed_form(fr, design[1], design[2], ratio = FALSE),
                 tolerance = 1e-10)
    expect_equal(prec_cluster(n = design, icc = ratio)$cv^2,
                 .closed_form(fr, design[1], design[2], ratio = TRUE),
                 tolerance = 1e-10)
  }
  expect_lt(ratio$varb, expansion$varb)
  expect_equal(ratio$varw, expansion$varw)
  expect_equal(ratio$unit_relvar, expansion$unit_relvar)

  vec <- varcomp(fr$y, stage_id = list(fr$psu), estimator = "ratio")
  expect_equal(vec$varb, ratio$varb)
})

test_that("the estimators agree at equal sizes and under size-proportional draws", {
  set.seed(3)
  eq <- data.frame(psu = rep(1:20, each = 10))
  eq$y <- 5 + rep(rnorm(20), each = 10) + rnorm(200)
  a <- varcomp(y ~ psu, data = eq)
  b <- varcomp(y ~ psu, data = eq, estimator = "ratio")
  expect_equal(b$varb, a$varb, tolerance = 1e-12)
  expect_equal(b$icc, a$icc, tolerance = 1e-12)

  fr <- .unequal_frame()
  fr$pp <- ave(fr$y, fr$psu, FUN = length) / nrow(fr)
  a <- varcomp(y ~ psu, data = fr, prob = ~pp)
  b <- varcomp(y ~ psu, data = fr, prob = ~pp, estimator = "ratio")
  expect_equal(b$varb, a$varb, tolerance = 1e-10)
  expect_equal(b$icc, a$icc, tolerance = 1e-10)
})

test_that("three stages, strata and weights carry the ratio estimator", {
  fr <- .unequal_frame()
  fr$ssu <- ave(seq_along(fr$y), fr$psu, FUN = function(i) {
    (seq_along(i) - 1) %/% 4
  })
  two <- varcomp(y ~ psu, data = fr, estimator = "ratio")
  three <- varcomp(y ~ psu / ssu, data = fr, estimator = "ratio")
  expect_equal(three$varb, two$varb, tolerance = 1e-12)
  expect_identical(three$estimator, "ratio")
  expect_lt(three$varw[["varw_psu"]],
            varcomp(y ~ psu / ssu, data = fr)$varw[["varw_psu"]])

  fr$st <- ifelse(fr$psu <= 20, "a", "b")
  by_stratum <- varcomp(y ~ psu, data = fr, strata = ~st, estimator = "ratio")
  expect_identical(by_stratum$estimator, "ratio")
  expect_equal(
    by_stratum$strata$icc_psu[1],
    varcomp(y ~ psu, data = fr[fr$st == "a", ], estimator = "ratio")$icc
  )

  set.seed(5)
  kept <- unlist(lapply(split(seq_len(nrow(fr)), fr$psu), function(i) {
    i[sort(sample.int(length(i), ceiling(length(i) / 2)))]
  }))
  sub <- fr[kept, ]
  sub$w <- ave(fr$y, fr$psu, FUN = length)[kept] /
    ave(sub$y, sub$psu, FUN = length)
  weighted <- varcomp(y ~ psu, data = sub, weights = ~w, estimator = "ratio")
  expect_identical(weighted$input, "weighted")
  sub$z <- sub$y - sum(sub$w * sub$y) / sum(sub$w)
  centered <- suppressWarnings(varcomp(z ~ psu, data = sub, weights = ~w))
  expect_equal(weighted$icc, centered$icc, tolerance = 1e-10)
  expect_equal(weighted$var_ratio, centered$var_ratio, tolerance = 1e-10)
})

test_that("a survey design passes the estimator through", {
  skip_if_not_installed("survey")
  fr <- .unequal_frame()
  fr$w <- 1
  design <- survey::svydesign(ids = ~psu, weights = ~w, data = fr)
  expect_equal(
    varcomp(design, ~y, estimator = "ratio")$icc,
    varcomp(y ~ psu, data = fr, estimator = "ratio")$icc,
    tolerance = 1e-10
  )
})

test_that("the estimator is named in the print and refused on the back-out", {
  fr <- .unequal_frame()
  expect_true(any(grepl("ratio estimator", capture.output(
    print(varcomp(y ~ psu, data = fr, estimator = "ratio"))
  ), fixed = TRUE)))
  expect_false(any(grepl("estimator", capture.output(
    print(varcomp(y ~ psu, data = fr))
  ), fixed = TRUE)))
  expect_error(varcomp(deff = 2, n_per_psu = 10, estimator = "ratio"),
               "without 'estimator'")
  expect_true(is.na(varcomp(deff = 2, n_per_psu = 10)$estimator))
  expect_error(varcomp(y ~ psu, data = fr, estimator = "mean"), "should be one of")
})
