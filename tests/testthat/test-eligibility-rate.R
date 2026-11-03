## Identity with the folded rate at an infinite population

test_that("eligibility_rate equals folding it into resp_rate when N is infinite", {
  for (e in c(0.05, 0.3, 0.7, 1)) {
    r <- 0.9
    expect_equal(
      n_prop(p = 0.3, moe = 0.03, resp_rate = r, eligibility_rate = e)$n,
      n_prop(p = 0.3, moe = 0.03, resp_rate = r * e)$n
    )
    expect_equal(
      n_mean(var = 4, moe = 0.2, resp_rate = r, eligibility_rate = e)$n,
      n_mean(var = 4, moe = 0.2, resp_rate = r * e)$n
    )
    expect_equal(
      n_ratio(r = 0.5, cv_num = 0.8, cv_den = 0.5, component_cor = 0.3,
              cv = 0.05, resp_rate = r, eligibility_rate = e)$n,
      n_ratio(r = 0.5, cv_num = 0.8, cv_den = 0.5, component_cor = 0.3,
              cv = 0.05, resp_rate = r * e)$n
    )
    expect_equal(
      prec_prop(p = 0.3, n = 500, resp_rate = r, eligibility_rate = e)$se,
      prec_prop(p = 0.3, n = 500, resp_rate = r * e)$se
    )
    expect_equal(
      prec_mean(var = 4, n = 500, resp_rate = r, eligibility_rate = e)$se,
      prec_mean(var = 4, n = 500, resp_rate = r * e)$se
    )
    expect_equal(
      power_prop(p1 = 0.3, p2 = 0.4, resp_rate = r, eligibility_rate = e)$n,
      power_prop(p1 = 0.3, p2 = 0.4, resp_rate = r * e)$n
    )
    expect_equal(
      power_mean(effect = 0.5, var = 4, resp_rate = r,
                 eligibility_rate = e)$n,
      power_mean(effect = 0.5, var = 4, resp_rate = r * e)$n
    )
    expect_equal(
      power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                outcome = "prop", resp_rate = r, eligibility_rate = e)$n,
      power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                outcome = "prop", resp_rate = r * e)$n
    )
    expect_equal(
      n_change(var = 4, change = 0.5, moe = 0.3, overlap = 0.5,
               overlap_cor = 0.4, resp_rate = r, eligibility_rate = e)$n,
      n_change(var = 4, change = 0.5, moe = 0.3, overlap = 0.5,
               overlap_cor = 0.4, resp_rate = r * e)$n
    )
    expect_equal(
      n_pooled(p = 0.3, occasions = 4, moe = 0.02, overlap = 0.5,
               overlap_cor = 0.5, resp_rate = r, eligibility_rate = e)$n,
      n_pooled(p = 0.3, occasions = 4, moe = 0.02, overlap = 0.5,
               overlap_cor = 0.5, resp_rate = r * e)$n
    )
    expect_equal(
      n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                resp_rate = r, eligibility_rate = e)$n,
      n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                resp_rate = r * e)$n
    )
    expect_equal(
      n_cluster(stage_cost = c(500, 100, 20), icc = c(0.05, 0.02), cv = 0.05,
                resp_rate = r, eligibility_rate = e)$n,
      n_cluster(stage_cost = c(500, 100, 20), icc = c(0.05, 0.02), cv = 0.05,
                resp_rate = r * e)$n
    )
    expect_equal(
      effective_n(n = 1000, deff = 2, resp_rate = r, eligibility_rate = e),
      effective_n(n = 1000, deff = 2, resp_rate = r * e)
    )
  }
})

## Exactness at a finite population

test_that("with a finite N the population in scope is N * eligibility_rate", {
  for (e in c(0.3, 0.7, 0.95)) {
    for (N in c(800, 3000, 1e5)) {
      expect_equal(
        n_prop(p = 0.3, moe = 0.05, N = N, resp_rate = 0.9,
               eligibility_rate = e)$n,
        n_prop(p = 0.3, moe = 0.05, N = N * e, resp_rate = 0.9)$n / e
      )
      expect_equal(
        n_mean(var = 4, moe = 0.3, N = N, resp_rate = 0.9,
               eligibility_rate = e)$n,
        n_mean(var = 4, moe = 0.3, N = N * e, resp_rate = 0.9)$n / e
      )
      expect_equal(
        prec_mean(var = 4, n = 200, N = N, resp_rate = 0.9,
                  eligibility_rate = e)$se,
        prec_mean(var = 4, n = 200 * e, N = N * e, resp_rate = 0.9)$se
      )
      if (N * e > 1000) {
        expect_equal(
          power_mean(effect = 0.5, var = 4, N = N, resp_rate = 0.9,
                     eligibility_rate = e)$n,
          power_mean(effect = 0.5, var = 4, N = N * e, resp_rate = 0.9)$n / e
        )
      }
      if (N * e > 1000) {
        expect_equal(
          n_change(p = c(0.3, 0.35), moe = 0.06, N = N, resp_rate = 0.9,
                   eligibility_rate = e, overlap = 0.5, overlap_cor = 0.4)$n,
          n_change(p = c(0.3, 0.35), moe = 0.06, N = N * e, resp_rate = 0.9,
                   overlap = 0.5, overlap_cor = 0.4)$n / e
        )
      }
      expect_equal(
        n_pooled(p = 0.3, occasions = 3, moe = 0.05, N = N, resp_rate = 0.9,
                 eligibility_rate = e, overlap = 0.5, overlap_cor = 0.5)$n,
        n_pooled(p = 0.3, occasions = 3, moe = 0.05, N = N * e,
                 resp_rate = 0.9, overlap = 0.5, overlap_cor = 0.5)$n / e
      )
    }
  }
})

test_that("an infeasible overlap names the population in scope", {
  expect_error(
    n_change(p = c(0.3, 0.35), moe = 0.06, N = 800, resp_rate = 0.9,
             eligibility_rate = 0.3, overlap = 0.5, overlap_cor = 0.4),
    "population in scope .* is 240"
  )
})

test_that("a census of the frame caps the issued size at N", {
  res <- n_prop(p = 0.3, moe = 0.01, N = 900, resp_rate = 1,
                eligibility_rate = 1)
  expect_lte(res$n, 900)
  expect_error(
    n_prop(p = 0.3, moe = 0.01, N = 900, resp_rate = 0.9,
           eligibility_rate = 0.5),
    "eligibility_rate = 0.5"
  )
})

## Equivalence with eligible_per_unit at or below 1

test_that("eligibility_rate combines with eligible_per_unit as one yield", {
  expect_equal(
    n_prop(p = 0.3, moe = 0.03, N = 5000, resp_rate = 0.9,
           eligibility_rate = 0.9, eligible_per_unit = 1.7)$n,
    n_prop(p = 0.3, moe = 0.03, N = 5000, resp_rate = 0.9,
           eligible_per_unit = 0.9 * 1.7)$n
  )
  for (fpc in c("none", "stage")) {
    pop <- if (fpc == "stage") list(N = 40000, N_psu = 400) else list()
    a <- do.call(n_cluster, c(list(stage_cost = c(500, 50), icc = 0.05,
                                   cv = 0.05, resp_rate = 0.9,
                                   eligibility_rate = 0.85, fpc = fpc), pop))
    b <- do.call(n_cluster, c(list(stage_cost = c(500, 50), icc = 0.05,
                                   cv = 0.05, resp_rate = 0.9,
                                   eligible_per_unit = 0.85, fpc = fpc), pop))
    expect_equal(a$n, b$n)
    expect_equal(a$cost, b$cost)
  }
})

## Round trips

test_that("the inverse functions carry eligibility_rate both ways", {
  x <- n_prop(p = 0.3, moe = 0.03, N = 5000, resp_rate = 0.9,
              eligibility_rate = 0.9)
  expect_equal(prec_prop(x)$moe, 0.03)
  expect_equal(n_prop(prec_prop(x))$n, x$n)
  m <- n_mean(var = 4, moe = 0.2, N = 800, resp_rate = 0.9,
              eligibility_rate = 0.8)
  expect_equal(prec_mean(m)$moe, 0.2)
  q <- n_ratio(r = 0.5, cv_num = 0.8, cv_den = 0.5, component_cor = 0.3,
               cv = 0.05, N = 900, resp_rate = 0.9, eligibility_rate = 0.9)
  expect_equal(prec_ratio(q)$cv, 0.05)
  ch <- n_change(p = c(0.3, 0.35), moe = 0.04, N = 3000, resp_rate = 0.9,
                 eligibility_rate = 0.9, overlap = 0.5, overlap_cor = 0.4)
  expect_equal(prec_change(ch)$moe, 0.04)
  po <- n_pooled(p = 0.3, occasions = 4, moe = 0.02, N = 5000,
                 resp_rate = 0.9, eligibility_rate = 0.9, overlap = 0.5,
                 overlap_cor = 0.5)
  expect_equal(prec_pooled(po)$moe, 0.02)
  cl <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                  resp_rate = 0.9, eligibility_rate = 0.9, fpc = "stage",
                  N = 40000, N_psu = 400)
  expect_equal(n_cluster(prec_cluster(cl))$n, cl$n)
  pa <- n_panel(n_prop(p = 0.3, moe = 0.03, N = 20000),
                retention = c(0.9, 0.9), resp_rate = 0.8,
                eligibility_rate = 0.95)
  expect_equal(prec_panel(pa)$n_resp, pa$n_resp)
})

test_that("a panel target sized with the panel's rates raises no warning", {
  tgt <- n_prop(p = 0.3, moe = 0.03, resp_rate = 0.8, eligibility_rate = 0.95)
  expect_no_warning(
    pa <- n_panel(tgt, retention = c(0.9, 0.9), resp_rate = 0.8,
                  eligibility_rate = 0.95)
  )
  expect_equal(pa$n_target, n_prop(p = 0.3, moe = 0.03)$n)
  expect_warning(
    n_panel(tgt, retention = c(0.9, 0.9), resp_rate = 0.8),
    "response 0.8 times eligibility 0.95 \\(0.76\\) is not used"
  )
  expect_warning(
    n_panel(tgt, retention = c(0.9, 0.9), resp_rate = 0.8,
            eligibility_rate = 0.9),
    "'resp_rate' times 'eligibility_rate'"
  )
})

## Several indicators

test_that("n_multi() sizes every row with its own eligibility_rate", {
  ind <- data.frame(
    name = c("a", "b", "c"), p = c(0.3, 0.1, NA), var = c(NA, NA, 4),
    mu = c(NA, NA, 10), moe = c(0.03, 0.02, 0.3), N = 5000,
    eligibility_rate = c(0.8, NA, 1)
  )
  res <- n_multi(ind, resp_rate = 0.9, eligibility_rate = 0.9)
  rows <- c(
    n_prop(p = 0.3, moe = 0.03, N = 5000, resp_rate = 0.9,
           eligibility_rate = 0.8)$n,
    n_prop(p = 0.1, moe = 0.02, N = 5000, resp_rate = 0.9,
           eligibility_rate = 0.9)$n,
    n_mean(var = 4, mu = 10, moe = 0.3, N = 5000, resp_rate = 0.9)$n
  )
  expect_equal(res$detail$.n, rows)
  expect_equal(prec_multi(res)$detail$.moe[1L], 0.03)
  floor <- n_multi(data.frame(p = 0.3, moe = 0.03, min_cases = 200),
                   resp_rate = 0.9, eligibility_rate = 0.8)
  expect_equal(
    floor$n,
    max(200 / (0.3 * 0.9 * 0.8),
        n_prop(p = 0.3, moe = 0.03, resp_rate = 0.9,
               eligibility_rate = 0.8)$n)
  )
  scarce <- n_multi(data.frame(p = 0.3, moe = 0.1, min_cases = 200),
                    resp_rate = 0.9, eligibility_rate = 0.8)
  expect_equal(scarce$n, 200 / (0.3 * 0.9 * 0.8))
})

test_that("cluster indicators take eligibility_rate as one argument", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.1), cv = c(0.05, 0.1),
                    icc_psu = c(0.05, 0.02), resp_rate = c(0.9, 0.8))
  folded <- transform(ind, resp_rate = resp_rate * 0.9)
  a <- n_cluster(stage_cost = c(500, 50), indicators = ind,
                 eligibility_rate = 0.9, fpc = "stage", N = 40000,
                 N_psu = 400)
  b <- n_cluster(stage_cost = c(500, 50), indicators = folded,
                 fpc = "stage", N = 36000, N_psu = 400)
  expect_equal(a$n, b$n)
  expect_equal(a$indicators$resp_rate, ind$resp_rate)
  expect_equal(n_cluster(prec_cluster(a))$n, a$n, tolerance = 1e-6)
  expect_error(
    n_cluster(stage_cost = c(500, 50),
              indicators = transform(ind, eligibility_rate = 0.9)),
    "give the argument 'eligibility_rate'"
  )
})

## Stratified allocation under the exact model

test_that("n_alloc() weights strata by their population in scope", {
  fr <- data.frame(stratum = c("u", "r"), N = c(5000, 5000), sd = c(10, 10),
                   mean = c(50, 50))
  e <- c(0.95, 0.75)
  ney <- n_alloc(fr, n = 800, resp_rate = 0.9, eligibility_rate = e)
  expect_equal(ney$detail$n, 800 * sqrt(e) / sum(sqrt(e)))
  prop <- n_alloc(fr, n = 800, alloc = "proportional", resp_rate = 0.9,
                  eligibility_rate = e)
  resp <- prop$detail$n * 0.9 * e
  expect_equal(resp / sum(resp), e / sum(e))
  W <- e / sum(e)
  n_h <- ney$detail$n
  se <- sqrt(sum(W^2 * 100 / (n_h * 0.9 * e) * (1 - n_h * 0.9 / 5000)))
  expect_equal(prec_alloc(ney)$se, se)
  expect_equal(
    n_alloc(transform(fr, eligibility_rate = e), n = 800, resp_rate = 0.9)$detail$n,
    ney$detail$n
  )
})

test_that("a scalar eligibility_rate in n_alloc() matches a frame of N * e", {
  fr <- data.frame(stratum = c("u", "r"), N = c(5000, 8000), sd = c(10, 14),
                   mean = c(50, 60))
  a <- n_alloc(fr, cv = 0.01, resp_rate = 0.9, eligibility_rate = 0.8)
  b <- n_alloc(transform(fr, N = N * 0.8), cv = 0.01, resp_rate = 0.9)
  expect_equal(a$detail$n, b$detail$n / 0.8)
  cl <- data.frame(stratum = c("u", "r"), N = c(50000, 50000),
                   N_psu = c(400, 500), sd = c(10, 10), mean = c(50, 50),
                   icc_psu = c(0.05, 0.08), n_per_psu = c(20, 20))
  for (fpc in c("unit", "stage")) {
    a <- n_alloc(cl, cv = 0.01, resp_rate = 0.9, eligibility_rate = 0.9,
                 fpc = fpc)
    b <- n_alloc(transform(cl, N = N * 0.9, n_per_psu = n_per_psu * 0.9),
                 cv = 0.01, resp_rate = 0.9, fpc = fpc)
    expect_equal(a$detail$n, b$detail$n / 0.9)
  }
  ms <- data.frame(stratum = rep(c("u", "r"), 2),
                   name = rep(c("y", "z"), each = 2),
                   mean = c(50, 50, 10, 12), sd = c(10, 10, 4, 5))
  tg <- data.frame(name = c("y", "z"), cv = c(0.01, 0.02))
  a <- n_alloc(fr[1:2], measures = ms, targets = tg, resp_rate = 0.9,
               eligibility_rate = 0.9)
  b <- n_alloc(transform(fr[1:2], N = N * 0.9), measures = ms, targets = tg,
               resp_rate = 0.9)
  expect_equal(a$detail$n, b$detail$n / 0.9)
})

test_that("allocation output reads the population in scope", {
  fr <- data.frame(stratum = c("u", "r"), N = c(5000, 5000), sd = c(10, 10),
                   mean = c(50, 50), eligibility_rate = c(0.95, 0.75))
  a <- n_alloc(fr, n = 800, resp_rate = 0.9)
  expect_equal(prec_alloc(a)$cv, a$cv)
  expect_equal(
    effective_n(a),
    sum(a$detail$n * 0.9 * c(0.95, 0.75)) / as.double(design_effect(a))
  )
  g <- predict(n_alloc(fr[1:4], n = 800, resp_rate = 0.9),
               data.frame(eligibility_rate = c(0.8, 1)))
  expect_gt(g$se[1L], g$se[2L])
  out <- capture.output(summary(a))
  expect_true(any(grepl("eligibility rate: 0.75 to 0.95", out)))
})

## Display and storage

test_that("eligibility_rate is printed and stored only when it is set", {
  x <- n_prop(p = 0.3, moe = 0.03, resp_rate = 0.9, eligibility_rate = 0.95)
  expect_equal(x$params[["eligibility_rate"]], 0.95)
  expect_true(any(grepl("\\beligibility_rate = 0.95", capture.output(x))))
  y <- n_prop(p = 0.3, moe = 0.03, resp_rate = 0.9)
  expect_null(y$params[["eligibility_rate"]])
  expect_false(any(grepl("eligibility_rate", capture.output(y))))
  expect_equal(x$expected_cases, x$n * 0.9 * 0.95 * 0.3)
  pw <- power_prop(p1 = 0.3, p2 = 0.4, resp_rate = 0.9, eligibility_rate = 0.9)
  expect_true(any(grepl("\\beligibility_rate = 0.90", capture.output(pw))))
  g <- predict(x, data.frame(eligibility_rate = c(0.9, 1)))
  expect_equal(g$n[2L], n_prop(p = 0.3, moe = 0.03, resp_rate = 0.9)$n)
})

## Refusals

test_that("eligibility_rate is refused outside (0, 1]", {
  for (bad in list(0, 1.2, NA_real_, "a", c(0.9, 0.8))) {
    expect_error(n_prop(p = 0.3, moe = 0.03, eligibility_rate = bad),
                 "eligibility_rate")
  }
  fr <- data.frame(stratum = c("u", "r"), N = c(5000, 5000), sd = c(10, 10))
  expect_error(n_alloc(fr, n = 800, eligibility_rate = c(0.9, 0.8, 0.7)),
               "eligibility_rate")
  expect_error(n_alloc(transform(fr, eligibility_rate = c(0.9, 1.1)), n = 800),
               "eligibility_rate")
  plan <- svyplan(eligibility_rate = 0.9)
  expect_equal(n_prop(p = 0.3, moe = 0.03, plan = plan)$n,
               n_prop(p = 0.3, moe = 0.03, eligibility_rate = 0.9)$n)
  expect_error(svyplan(eligibility_rate = 1.5), "eligibility_rate")
})

## PSU register

test_that("the register design effect reads the responding take in scope", {
  z <- .psu_fixture()
  sizes <- z$psu$N[z$psu$stratum == "A"]
  certain <- sizes >= sort(sizes, decreasing = TRUE)[6]
  per <- pmin(ceiling(0.02 * sizes[certain]), sizes[certain])
  direct <- function(e, k = 20, take = 12, icc = 0.05, r = 0.9) {
    Ne <- sum(sizes) * e
    held <- sizes[certain] * e
    rest <- Ne - sum(held)
    n <- sum(per) + k * take
    v <- sum(held^2 * (1 / (r * e * per) - 1 / held)) +
      rest^2 * (1 + icc * (take * r * e - 1)) *
      (1 / (r * e * k * take) - 1 / rest)
    max(v / (Ne^2 * (1 / (r * e * n) - 1 / Ne)), 1)
  }
  for (e in c(1, 0.9, 0.6, 0.3)) {
    expect_equal(
      svyplan:::.psu_field_deff(sizes, certain, per, 20, 12, 0.05, 0.9, e),
      direct(e)
    )
  }
})

test_that("a register allocation takes eligibility_rate per stratum", {
  z <- .psu_fixture()
  base <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = z$psu, resp_rate = 0.9)
  same <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = z$psu, resp_rate = 0.9, eligibility_rate = 1)
  expect_identical(same$detail$n, base$detail$n)
  fr <- transform(z$frame, eligibility_rate = c(0.95, 0.9, 0.8, 0.7))
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 psu = z$psu, resp_rate = 0.9)
  expect_false(isTRUE(all.equal(fit$detail$n, base$detail$n)))
  assessed <- prec_alloc(fit)$detail
  expect_true(all(assessed$.cv <= z$targets$cv + 1e-9))
  bad <- transform(z$measures, eligibility_rate = 0.9)
  expect_error(
    n_alloc(z$frame, measures = bad, targets = z$targets, psu = z$psu),
    "stratum rate"
  )
})

test_that("a register solve is the joint solve at the register design effect", {
  z <- .psu_fixture()
  e <- c(0.95, 0.9, 0.8, 0.7)
  fr <- transform(z$frame, eligibility_rate = e)
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 psu = z$psu, resp_rate = 0.9)
  certain <- fit$psu$certainty
  h_of <- match(z$measures$stratum, fr$stratum)
  deff_r <- vapply(seq_len(nrow(z$measures)), function(r) {
    h <- h_of[r]
    i <- z$psu$stratum == fr$stratum[h]
    svyplan:::.psu_deff(fr$N[h], sum(z$psu$N[i & certain]),
                        sum(z$psu$N[i & !certain]), z$measures$icc_psu[r],
                        12 * 0.9 * e[h])
  }, numeric(1))
  joint <- n_alloc(
    fr[c("stratum", "N", "eligibility_rate")],
    measures = transform(z$measures[setdiff(names(z$measures), "icc_psu")],
                         deff = deff_r),
    targets = z$targets, resp_rate = 0.9
  )
  expect_equal(fit$detail$n, joint$detail$n, tolerance = 1e-3)
})

test_that("a register state's score is its variance over the population in scope", {
  sizes <- c(4000, 3000, 900, 800, 700, 600)
  held <- sizes[1:2]
  per <- c(40, 30)
  a <- 3
  b <- 12
  r <- 0.9
  icc <- 0.05
  N <- sum(sizes)
  states <- data.frame(Q = sum(per), a = a, held = sum(held),
                       A = sum(held^2 / per), n_c = 2)
  for (e in c(1, 0.8, 0.5)) {
    ev <- list(
      N = N, take = b, unit = 1, cost_psu = NULL, objective = NULL,
      target = list(w = matrix(1), weight = 1, resp = matrix(r),
                    scope = matrix(e), icc = matrix(icc))
    )
    rest <- (N - sum(held)) * e
    direct <- sum((held * e)^2 * (1 / (r * e * per) - 1 / (held * e))) +
      rest^2 * (1 + icc * (b * r * e - 1)) * (1 / (r * e * a * b) - 1 / rest)
    expect_equal(svyplan:::.psu_eval_states(ev, 1L, states)$C[1, 1], direct)
  }
})

test_that("the register search bound is the scorer's stratum-rate floor", {
  N <- 10000
  r <- 0.9
  n <- 400
  for (e in c(1, 0.8, 0.5)) {
    ev <- list(
      N = N, take = 12, unit = 1, cost_psu = NULL, objective = NULL,
      lower = 1,
      target = list(w = matrix(2), weight = 0.5, resp = matrix(r),
                    scope = matrix(e), icc = matrix(0.05))
    )
    bd <- svyplan:::.psu_menu_bound(ev, 1L, list(mu = 1, price = 0, obj = 0))
    expect_equal(bd$W / n - bd$W0, 2 * 0.5 * e * N^2 * (1 / (r * n) - 1 / N))
  }
})

## The design effect and the frame bounds read the rate

test_that("a cluster design effect reads the take in scope", {
  for (cost in list(c(500, 50), c(500, 100, 50))) {
    icc <- if (length(cost) == 2L) 0.05 else c(0.05, 0.1)
    scoped <- n_cluster(stage_cost = cost, icc = icc, cv = 0.05,
                        resp_rate = 0.9, eligibility_rate = 0.8)
    folded <- n_cluster(stage_cost = cost, icc = icc, cv = 0.05,
                        resp_rate = 0.72)
    expect_equal(as.double(design_effect(scoped)),
                 as.double(design_effect(folded)))
    expect_equal(effective_n(scoped), effective_n(folded))
    expect_equal(effective_n(scoped), 1 / 0.05^2)
  }
  scoped <- prec_cluster(icc = 0.05, n = c(n_psu = 50, n_per_psu = 10),
                         unit_relvar = 1, resp_rate = 0.9,
                         eligibility_rate = 0.8)
  folded <- prec_cluster(icc = 0.05, n = c(n_psu = 50, n_per_psu = 10),
                         unit_relvar = 1, resp_rate = 0.72)
  expect_equal(as.double(design_effect(scoped)),
               as.double(design_effect(folded)))
  expect_equal(effective_n(scoped), 1 / scoped$cv^2)
})

test_that("a power solve returns the size whose power it reached", {
  cases <- list(
    function(...) power_prop(p1 = 0.3, p2 = 0.33, N = 1000, ...),
    function(...) power_prop(p1 = 0.3, p2 = 0.33, N = 1000, overlap = 0.5,
                             overlap_cor = 0.5, ...),
    function(...) power_mean(var = 100, effect = 1.5, N = 1000, ...),
    function(...) power_mean(var = 100, effect = 1.5, N = 1000,
                             overlap = 0.5, overlap_cor = 0.5, ...),
    function(...) power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                            outcome = "prop", N = 6000, resp_rate = 0.9,
                            overlap = 0.5, overlap_cor = 0.4, ...)
  )
  for (call in cases) {
    for (e in c(1, 0.5, 0.4)) {
      at <- call(eligibility_rate = e, n = 300, power = NULL)$power
      back <- call(eligibility_rate = e, power = at)
      expect_equal(back$n, 300, tolerance = 1e-6)
      expect_equal(call(eligibility_rate = e, n = back$n, power = NULL)$power,
                   at, tolerance = 1e-8)
    }
  }
})

test_that("a gross sample may exceed the population in scope", {
  at <- power_prop(p1 = 0.3, p2 = 0.33, N = 1000, eligibility_rate = 0.5,
                   n = 800, power = NULL)$power
  fit <- power_prop(p1 = 0.3, p2 = 0.33, N = 1000, eligibility_rate = 0.5,
                    power = at)
  expect_equal(fit$n, 800, tolerance = 1e-6)
  fit <- power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                   outcome = "prop", N = 6000, resp_rate = 0.9,
                   eligibility_rate = 0.4, overlap = 0.5, overlap_cor = 0.4)
  expect_gt(fit$n, 6000 * 0.4 / 2)
  expect_equal(
    power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32), outcome = "prop",
              n = fit$n, power = NULL, N = 6000, resp_rate = 0.9,
              eligibility_rate = 0.4, overlap = 0.5, overlap_cor = 0.4)$power,
    0.8, tolerance = 1e-6
  )
})

test_that("an issued overlap profile is accepted at full response in scope", {
  rot <- design_overlap("4")
  scoped <- prec_pooled(var = 100, n = 500, occasions = 4, resp_rate = 1,
                        eligibility_rate = 0.8, overlap = rot,
                        cor_decay = 0.8)
  eligible <- prec_pooled(var = 100, n = 400, occasions = 4, overlap = rot,
                          cor_decay = 0.8)
  expect_equal(scoped$se, eligible$se)
  expect_error(
    prec_pooled(var = 100, n = 500, occasions = 4, resp_rate = 0.9,
                eligibility_rate = 0.8, overlap = rot, cor_decay = 0.8),
    "coincide only at resp_rate = 1"
  )
  fit <- n_pooled(var = 100, moe = 1, occasions = 4, eligibility_rate = 0.8,
                  overlap = rot, cor_decay = 0.8)
  res <- predict(fit, data.frame(eligibility_rate = c(0.7, 0.8)))
  expect_equal(res$n[2L], fit$n)
  expect_error(predict(fit, data.frame(resp_rate = 0.9)),
               "coincide only at resp_rate = 1")
})

test_that("expected cases and the case floor read every yield factor", {
  x <- n_prop(p = 0.2, moe = 0.05, resp_rate = 0.8, eligibility_rate = 0.9,
              eligible_per_unit = 1.3, N = 5000, min_cases = 200)
  expect_identical(x$binding, "min_cases")
  expect_equal(x$n, 200 / (0.2 * 0.8 * 0.9 * 1.3))
  expect_equal(x$expected_cases, x$n * 0.8 * 0.9 * 1.3 * 0.2)
})

test_that("the proportion variance reads the eligible respondents and population", {
  se <- prec_prop(p = 0.2, n = 100, N = 1000, resp_rate = 0.8,
                  eligibility_rate = 0.9, eligible_per_unit = 1.5)$se
  n_net <- 100 * 0.8 * 0.9 * 1.5
  N_elig <- 1000 * 0.9 * 1.5
  expect_equal(se, sqrt(0.2 * 0.8 * (N_elig - n_net) / (N_elig - 1) / n_net))
})
