## Whole-cluster designs in n mode

cluster_frame <- function(N_psu = NULL) {
  fr <- data.frame(
    stratum = c("North", "South"), N = c(40000, 60000),
    sd = c(12, 15), mean = c(50, 60),
    icc_psu = c(0.05, 0.08), var_ratio_psu = 1,
    cost_psu = c(400, 550), cost_ssu = c(45, 60)
  )
  if (!is.null(N_psu)) fr$N_psu <- N_psu
  fr
}

test_that("the whole take stays at its optimum whatever the total factors into", {
  fr <- cluster_frame()
  totals <- 1900:2060
  fits <- lapply(totals, function(n) n_alloc(fr, n = n))
  drift <- vapply(fits, function(x) {
    max(abs(x$detail$n_per_psu_int - x$detail$n_per_psu))
  }, numeric(1L))
  expect_true(all(drift < 1))

  cost <- vapply(fits, function(x) x$operational$cost, numeric(1L))
  psu_cost <- max(fr$cost_psu + fr$cost_ssu * 20)
  expect_true(all(abs(diff(cost)) <= psu_cost))
})

test_that("the field total is the closest the whole takes can reach", {
  fr <- cluster_frame()
  for (n in c(1999, 2000, 2001, 2003, 2311)) {
    d <- n_alloc(fr, n = n)$detail
    gap <- abs(n - sum(d$n_int))
    expect_lt(gap, max(d$n_per_psu_int))
    moved <- abs(n - sum(d$n_int) + outer(c(-1, 1), d$n_per_psu_int))
    expect_true(all(moved >= gap))
    expect_equal(d$n_int, d$n_psu_int * d$n_per_psu_int)
  }
})

test_that("a binding N_psu bounds the PSU count and loses no units", {
  fr <- cluster_frame(N_psu = c(40, 900))
  for (n in c(2000, 2500, 3000, 3500)) {
    x <- n_alloc(fr, n = n)
    d <- x$detail
    expect_identical(d$.bound_source[[1L]], "N_psu")
    expect_true(all(d$n_psu_int <= fr$N_psu))
    expect_equal(d$n_psu_int[[1L]], 40L)
    expect_lt(abs(sum(d$n_int) - n), max(d$n_per_psu_int))
    expect_equal(x$operational$n, sum(d$n_int))
  }
})

test_that("moves in several strata at once close the total", {
  fit <- n_alloc(cluster_frame(), n = 1999)
  d <- fit$detail
  expect_equal(d$n_per_psu_int, c(13L, 10L))
  expect_equal(sum(d$n_int), 1999L)
  expect_equal(fit$operational$n, 1999)
})

test_that("one PSU per stratum reaches the closest total those moves allow", {
  set.seed(8)
  for (i in seq_len(300)) {
    H <- sample(2:6, 1)
    b <- sample(1:30, H, replace = TRUE)
    e <- sample(50:2000, H, replace = TRUE)
    a <- pmax(1L, as.integer(round(e / b)))
    lo <- pmax(1L, a - sample(0:2, H, replace = TRUE))
    hi <- a + sample(0:2, H, replace = TRUE)
    got <- svyplan:::.cluster_close_total(a, b, lo, hi, e)
    expect_true(all(got >= lo & got <= hi & abs(got - a) <= 1L))
    grid <- as.matrix(expand.grid(lapply(seq_len(H), function(h) {
      seq.int(max(lo[h], a[h] - 1L), min(hi[h], a[h] + 1L))
    })))
    gap <- abs(sum(e) - sum(got * b))
    expect_equal(gap, min(abs(sum(e) - grid %*% b)))
    if (all(lo < a & hi > a)) expect_lte(gap, max(b) / 2)
  }
})

test_that("a take is chosen only if its PSU count fits the universe", {
  south <- cluster_frame(N_psu = c(40, 40))[2, ]
  fit <- n_alloc(south, n = 405, min_n_stratum = 405)
  expect_lte(fit$detail$n_psu_int, 40L)
  expect_gte(fit$operational$n, 405)
  expect_equal(fit$detail$n_int, fit$detail$n_psu_int * fit$detail$n_per_psu_int)
})

test_that("a capped stratum's shortfall moves to strata with room", {
  fr <- cluster_frame(N_psu = c(5000, 1000))
  fit <- n_alloc(fr, n = 20000)
  d <- fit$detail
  expect_identical(d$.bound_source[[2L]], "N_psu")
  expect_true(all(d$n_psu_int <= fr$N_psu))
  expect_equal(fit$operational$n, 20000)
  expect_gt(d$n_int[[1L]], round(d$n[[1L]]))
})

test_that("a take rises rather than the total falling when every stratum is capped", {
  fr <- cluster_frame(N_psu = c(1000, 1000))
  fit <- n_alloc(fr, n = 23200)
  d <- fit$detail
  expect_equal(fit$operational$n, 23200)
  expect_true(all(d$n_psu_int <= fr$N_psu))
  expect_equal(d$n_int, d$n_psu_int * d$n_per_psu_int)
  for (n in seq(22800, 23250, by = 50)) {
    op <- n_alloc(fr, n = n)$operational$n
    expect_lte(abs(op - n), 7, label = paste("n =", n))
  }
})

test_that("a shortfall within half a take keeps the cost-optimal take", {
  fr <- data.frame(
    stratum = "A", N = 997, sd = 10, mean = 50, icc_psu = 0.05,
    var_ratio_psu = 1, cost_psu = 1500, cost_ssu = 20, N_psu = 50
  )
  fit <- n_alloc(fr, n = 997)
  d <- fit$detail
  expect_lt(abs(d$n_per_psu_int - d$n_per_psu), 1)
  expect_lte(997 - fit$operational$n, d$n_per_psu_int / 2)
  expect_lte(d$n_psu_int, 50L)
})

test_that("an overshoot beyond half a take is an error naming the smallest design", {
  fr <- data.frame(
    stratum = c("a", "b"), N = c(1000, 1000), sd = 10, mean = 50,
    icc_psu = 0.05, n_per_psu = 10, cost_psu = 500, cost_ssu = 50
  )
  expect_error(n_alloc(fr, n = 30, min_n_stratum = 15), "holds 40 units")
  expect_equal(n_alloc(fr, n = 36, min_n_stratum = 15)$operational$n, 40)

  five <- data.frame(
    stratum = letters[1:5], N = 1000, sd = 10, mean = 50, icc_psu = 0.05,
    n_per_psu = 10, cost_psu = 500, cost_ssu = 50
  )
  expect_error(n_alloc(five, n = 30), "holds 50 units")
})

## The field design and prec_alloc() agree

test_that("prec_alloc() at the fielded stage sizes reproduces $operational", {
  fr <- cluster_frame()
  fits <- list(
    n = n_alloc(fr, n = 2000),
    cv = n_alloc(fr, cv = 0.01),
    budget = n_alloc(fr, budget = 3e5),
    stage = n_alloc(transform(fr, N_psu = c(2000, 3000)), n = 2000,
                    fpc = "stage"),
    capped = n_alloc(cluster_frame(N_psu = c(40, 900)), n = 3000)
  )
  for (mode in names(fits)) {
    fit <- fits[[mode]]
    d <- fit$detail
    field <- prec_alloc(fit, n = d$n_int, n_per_psu = d$n_per_psu_int)
    expect_equal(field$cv, fit$operational$cv, tolerance = 1e-12,
                 label = mode)
    expect_equal(field$se, fit$operational$se, tolerance = 1e-12,
                 label = mode)
  }
})

test_that("an element allocation's field design needs no take", {
  fr <- data.frame(stratum = c("a", "b", "c"), N = c(1000, 2000, 500),
                   sd = c(10, 20, 5))
  fit <- n_alloc(fr, n = 301)
  expect_equal(prec_alloc(fit, n = fit$detail$n_int)$cv, fit$operational$cv,
               tolerance = 1e-12)
  expect_error(prec_alloc(fit, n_per_psu = 5), "no PSU stage")
})

test_that("the n_per_psu override is honoured and validated", {
  fit <- n_alloc(cluster_frame(), n = 2000)
  d <- fit$detail
  base <- prec_alloc(fit, n = d$n_int)$cv
  expect_false(isTRUE(all.equal(
    prec_alloc(fit, n = d$n_int, n_per_psu = d$n_per_psu_int)$cv, base
  )))
  expect_equal(
    prec_alloc(fit, n = d$n_int, n_per_psu = 10)$cv,
    prec_alloc(fit, n = d$n_int, n_per_psu = c(10, 10))$cv
  )
  expect_error(prec_alloc(fit, n_per_psu = c(1, 2, 3)), "one per stratum")
  expect_error(prec_alloc(fit, n_per_psu = 0), "positive")

  joint <- n_alloc(
    data.frame(stratum = c("a", "b"), N = c(2e5, 1e5), N_psu = c(5e3, 3e3),
               n_per_psu = 20, cost_psu = 300, cost_ssu = 25),
    measures = data.frame(stratum = c("a", "b"), name = "y", p = c(0.3, 0.4),
                          icc_psu = 0.05, var_ratio_psu = 1),
    targets = data.frame(name = "y", cv = 0.05)
  )
  same <- prec_alloc(joint, n_per_psu = 20)$detail$.achieved
  expect_equal(same, prec_alloc(joint)$detail$.achieved)
  expect_false(isTRUE(all.equal(
    prec_alloc(joint, n_per_psu = 40)$detail$.achieved, same
  )))
})

test_that("a PSU register keeps its own takes", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  expect_error(prec_alloc(fit, n_per_psu = 10), "PSU register")
})
