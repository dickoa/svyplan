test_that("effective_n from weights is the Kish effective sample size", {
  set.seed(1009)
  w <- runif(100, 1, 5)
  expect_equal(effective_n(weights = w), sum(w)^2 / sum(w^2))
})

test_that("effective_n equals n when weights are equal", {
  w <- rep(2.5, 60)
  expect_equal(effective_n(weights = w), 60)
})

test_that("effective_n is n divided by design_effect", {
  set.seed(444)
  w <- runif(80, 1, 6)
  expect_equal(
    effective_n(weights = w),
    length(w) / as.double(design_effect(weights = w))
  )
  expect_equal(
    effective_n(n = 1200, icc = 0.05, n_per_psu = 25),
    1200 / as.double(design_effect(icc = 0.05, n_per_psu = 25))
  )
})

test_that("effective_n accepts a design effect object or a plain number", {
  deff <- design_effect(icc = 0.05, n_per_psu = 25)
  expect_equal(effective_n(deff, n = 1200), 1200 / 2.2)
  expect_equal(effective_n(n = 1200, deff = deff), 1200 / 2.2)
  expect_equal(effective_n(n = 1200, deff = 2.2), 1200 / 2.2)
})

test_that("deff and its components are mutually exclusive", {
  expect_error(
    effective_n(n = 100, deff = 2, icc = 0.05, n_per_psu = 10),
    "either 'deff' or the design components"
  )
})

test_that("effective_n derives n where it can and requires it otherwise", {
  strata <- data.frame(N = c(50000, 120000), n = c(600, 400))
  expect_equal(
    effective_n(strata = strata),
    1000 / as.double(design_effect(strata = strata))
  )
  expect_error(effective_n(icc = 0.05, n_per_psu = 25), "'n' is required")
  expect_error(
    effective_n(design_effect(icc = 0.05, n_per_psu = 25)), "'n' is required"
  )
})

test_that("cluster planning bounds behave at the extremes", {
  expect_equal(effective_n(n = 500, icc = 0, n_per_psu = 30), 500)
  expect_equal(effective_n(n = 500, icc = 1, n_per_psu = 20), 25)
  expect_equal(effective_n(n = 500, icc = 0.05, n_per_psu = 1), 500)
})

test_that("effective_n reads plans and allocations", {
  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
  expect_equal(effective_n(plan), plan$total_n / as.double(design_effect(plan)))
  expect_equal(effective_n(plan, n = 1000), 1000 / as.double(design_effect(plan)))

  alloc <- n_alloc(
    data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 15),
               mean = c(50, 60)),
    n = 500
  )
  expect_equal(effective_n(alloc), alloc$n / as.double(design_effect(alloc)))
})

test_that("effective_n validates its inputs", {
  expect_error(effective_n(n = -5, icc = 0.05, n_per_psu = 10), "positive")
  expect_error(effective_n(weights = c(1, Inf)), "only finite")
  expect_error(effective_n(weights = c(1, -Inf)), "only finite")
  expect_error(effective_n(n = 100, icc = 0.05), "both 'icc' and 'n_per_psu'")
  expect_error(effective_n(n = 100, icc = 1.5, n_per_psu = 10), "in \\[0, 1\\]")
})

test_that("effective_n rejects unused arguments", {
  expect_error(effective_n(n = 100, icc = 0.05, n_per_psu = 10,
                           methd = "kish"),
               "unused argument.*methd")
})

test_that("effective_n nets down for response, matching the planning identity", {
  frame <- data.frame(
    stratum = c("A", "B"),
    N = c(1e5, 1e5),
    sd = c(1, 1),
    mean = c(1, 1)
  )
  alloc <- n_alloc(frame, n = 200, alloc = "neyman", resp_rate = 0.8)
  deff <- as.double(design_effect(alloc))

  # n * resp_rate / deff, the identity the n_eff column reports
  expect_equal(as.double(effective_n(alloc)), 200 * 0.8 / deff)
  expect_equal(as.double(effective_n(alloc)), sum(alloc$detail$n_eff))

  # an explicit rate overrides the plan's
  expect_equal(as.double(effective_n(alloc, resp_rate = 1)), 200 / deff)

  # a plan without nonresponse is unchanged
  full <- n_alloc(frame, n = 200, alloc = "neyman")
  expect_equal(as.double(effective_n(full)),
               200 / as.double(design_effect(full)))
})

test_that("effective_n takes a response rate for a bare n", {
  expect_equal(effective_n(n = 1000, deff = 2), 500)
  expect_equal(effective_n(n = 1000, deff = 2, resp_rate = 0.8), 400)
  expect_error(effective_n(n = 1000, deff = 2, resp_rate = 0), "resp_rate")
})

test_that("a cluster plan carries its own response rate", {
  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                    resp_rate_psu = 0.8)
  expect_equal(
    as.double(effective_n(plan)),
    plan$total_n * 0.8 / as.double(design_effect(plan))
  )
})

test_that("an allocation's effective size reads frame response and PSU loss", {
  frame <- data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                      mean = c(10, 20), sd = c(5, 10))
  from_frame <- n_alloc(transform(frame, resp_rate = c(0.5, 0.8)), n = 200)
  from_arg <- n_alloc(frame, n = 200, resp_rate = c(0.5, 0.8))
  expect_equal(from_frame$detail$n, from_arg$detail$n)
  expect_equal(effective_n(from_frame), effective_n(from_arg))
  expect_equal(
    effective_n(from_frame),
    sum(from_frame$detail$n * c(0.5, 0.8)) /
      as.double(design_effect(from_frame))
  )
  expect_equal(effective_n(from_frame), 202.0254, tolerance = 1e-6)
  # an explicit rate still replaces the stored one
  expect_equal(effective_n(from_frame, resp_rate = 1),
               200 / as.double(design_effect(from_frame)))

  # PSU nonresponse removes whole clusters from the same total
  cluster <- transform(frame, n_per_psu = 10, icc_psu = c(0.03, 0.06),
                       resp_rate_psu = c(0.5, 0.9))
  fit <- n_alloc(cluster, n = 1000)
  expect_equal(
    effective_n(fit),
    sum(fit$detail$n * c(0.5, 0.9)) / as.double(design_effect(fit))
  )
  expect_equal(effective_n(fit), 751.1398, tolerance = 1e-6)

  # eligibility enters once, beside unit and PSU response
  both <- n_alloc(
    transform(cluster, resp_rate = c(0.8, 0.7), resp_rate_psu = NULL,
              eligibility_rate = c(0.9, 0.6)),
    n = 1000
  )
  expect_equal(
    effective_n(both),
    sum(both$detail$n * c(0.8, 0.7) * c(0.9, 0.6)) /
      as.double(design_effect(both))
  )
})

## Effective size against the plan's own variance

pop_var_alloc <- function(frame) {
  e <- frame$eligibility_rate %||% rep(1, nrow(frame))
  W <- frame$N * e / sum(frame$N * e)
  mu <- sum(W * frame$mean)
  sum(W * frame$sd^2) + sum(W * (frame$mean - mu)^2)
}

test_that("stratum response rates enter the variance stratum by stratum", {
  frame <- data.frame(stratum = c("A", "B"), N = 10000, sd = 1, mean = 1,
                      resp_rate = c(0.1, 1))
  x <- n_alloc(frame, n = 200, min_n_stratum = 100, fpc = "none")
  expect_equal(x$detail$n, c(100, 100))
  V <- 0.5^2 / 10 + 0.5^2 / 100
  expect_equal(x$se^2, V)
  expect_equal(effective_n(x), 1 / V)
  expect_equal(as.double(design_effect(x)), 110 * V)
  expect_equal(prec_mean(var = 1, n = effective_n(x))$se, x$se)
})

test_that("effective_n() equals S^2 / V on constrained frames", {
  base <- data.frame(stratum = c("a", "b", "c"), N = c(800, 2500, 1200),
                     mean = c(10, 18, 25), sd = c(4, 9, 6),
                     resp_rate = c(0.35, 0.9, 0.6))
  frames <- list(
    element = transform(base, unit_cost = c(1, 2, 3)),
    eligibility = transform(base, eligibility_rate = c(0.7, 0.95, 0.8),
                            deff = c(1.2, 2, 1.5)),
    cluster = transform(base, icc_psu = c(0.04, 0.1, 0.07),
                        cost_psu = c(80, 150, 120), cost_ssu = c(5, 10, 8),
                        resp_rate_psu = c(0.85, 1, 0.9),
                        eligibility_rate = c(0.8, 0.9, 1))
  )
  for (nm in names(frames)) {
    fr <- frames[[nm]]
    x <- n_alloc(fr, n = 300, min_n_stratum = 80, fpc = "none")
    expect_equal(effective_n(x), pop_var_alloc(fr) / x$se^2,
                 label = nm)
    expect_equal(effective_n(x), pop_var_alloc(fr) / prec_alloc(x)$se^2,
                 label = nm)
  }
})

test_that("one rate in every stratum leaves the variance ratio unchanged", {
  fr <- data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                   mean = c(10, 20), sd = c(5, 10))
  x <- n_alloc(fr, n = 300, min_n_stratum = 120, fpc = "none")
  xr <- n_alloc(transform(fr, resp_rate = 0.6), n = 300, min_n_stratum = 120,
                fpc = "none")
  W <- fr$N / sum(fr$N)
  gross <- 300 * sum(W^2 * fr$sd^2 / x$detail$n) / pop_var_alloc(fr)
  expect_equal(as.double(design_effect(x)), gross)
  expect_equal(as.double(design_effect(xr)), gross)
  expect_equal(effective_n(xr), 300 * 0.6 / gross)
})

test_that("a rate override replaces each stratum's rate in the ratio too", {
  fr <- data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                   mean = c(10, 20), sd = c(5, 10))
  x <- n_alloc(fr, n = 300, min_n_stratum = 120, fpc = "none")
  n_h <- x$detail$n
  W <- fr$N / sum(fr$N)
  r <- c(0.3, 0.9)
  V <- sum(W^2 * fr$sd^2 / (n_h * r))
  expect_equal(effective_n(x, resp_rate = r), pop_var_alloc(fr) / V)
})

test_that("a size-only allocation charges Kish's loss on the units that respond", {
  fr <- data.frame(stratum = c("a", "b"), N = c(1000, 3000),
                   resp_rate = c(0.4, 0.9))
  x <- n_alloc(fr, n = 200, min_n_stratum = 90, alloc = "proportional")
  m <- x$detail$n_int * c(0.4, 0.9)
  W <- fr$N / sum(fr$N)
  expect_equal(as.double(design_effect(x)), sum(m) * sum(W^2 / m))

  scoped <- n_alloc(transform(fr, eligibility_rate = c(0.5, 1)), n = 200,
                    min_n_stratum = 90, alloc = "proportional")
  e <- c(0.5, 1)
  m <- scoped$detail$n_int * c(0.4, 0.9) * e
  W <- fr$N * e / sum(fr$N * e)
  expect_equal(as.double(design_effect(scoped)), sum(m) * sum(W^2 / m))
})

test_that("a precision allocation evaluates its own whole stratum counts", {
  fr <- data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                   mean = c(10, 20), sd = c(5, 10), deff = c(2, 3))
  x <- n_alloc(fr, cv = 0.3, fpc = "none")
  expect_equal(x$detail$n_int, c(1, 7))
  whole <- prec_alloc(x, n = x$detail$n_int)
  expect_equal(effective_n(whole), pop_var_alloc(fr) / whole$se^2)
  expect_equal(effective_n(whole), 3.951613, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(effective_n(whole),
                                effective_n(x, n = sum(x$detail$n_int)))))
  expect_equal(as.double(design_effect(whole)),
               8 / effective_n(whole))
  expect_error(design_effect(whole, strata = data.frame(N = 1, n = 1)),
               "carries its own strata")
  expect_error(effective_n(whole, strata = data.frame(N = 1, n = 1)),
               "carries its own strata")
})

test_that("an indicators-table cluster plan has no single design effect", {
  tab <- data.frame(name = c("y", "z"), p = c(0.3, 0.2), cv = 0.1,
                    icc_psu = c(0.05, 0.1))
  plain <- n_cluster(indicators = tab, stage_cost = c(100, 10))
  by_dom <- n_cluster(
    indicators = data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
                            region = c("A", "B")),
    stage_cost = c(100, 10), domains = "region"
  )
  for (x in list(plain, by_dom)) {
    expect_error(design_effect(x), "no\\s+single design effect")
    expect_error(effective_n(x), "no\\s+single design effect")
  }
})

## Rate overrides re-evaluate the plan

test_that("a response override moves the responding take and its penalty", {
  fr <- data.frame(stratum = "A", N = 10000, mean = 1, sd = 1,
                   icc_psu = 0.2, n_per_psu = 10)
  x <- n_alloc(fr, n = 200, fpc = "none")
  expect_equal(effective_n(x, resp_rate = 0.5), 100 / 1.8)
  p <- prec_alloc(transform(fr, resp_rate = 0.5), n = x$detail$n,
                  fpc = "none")
  expect_equal(effective_n(x, resp_rate = 0.5), 1 / p$se^2)
})

test_that("an eligibility override moves the weights and the population variance", {
  fr <- data.frame(stratum = c("a", "b"), N = c(1000, 2000),
                   mean = c(10, 20), sd = c(5, 10),
                   eligibility_rate = c(0.2, 1))
  x <- n_alloc(fr, n = 200, min_n_stratum = 100, fpc = "none")
  nf <- transform(fr, eligibility_rate = 1)
  p <- prec_alloc(nf, n = x$detail$n, fpc = "none")
  expect_equal(effective_n(x, eligibility_rate = 1),
               pop_var_alloc(nf) / p$se^2)
  expect_equal(effective_n(x, eligibility_rate = 1), 205.8824,
               tolerance = 1e-6)
})

test_that("overrides agree with prec_alloc() on the re-rated frame", {
  base <- data.frame(stratum = c("a", "b", "c"), N = c(3000, 5000, 4000),
                     mean = c(10, 18, 25), sd = c(4, 9, 6),
                     resp_rate = c(0.35, 0.9, 0.6),
                     eligibility_rate = c(0.5, 0.95, 0.8))
  frames <- list(
    element = base,
    cluster = transform(base, icc_psu = c(0.04, 0.15, 0.08),
                        n_per_psu = c(8, 12, 10),
                        resp_rate_psu = c(0.8, 1, 0.9))
  )
  overrides <- list(
    list(resp_rate = 0.7),
    list(eligibility_rate = 0.9),
    list(resp_rate = c(0.5, 0.6, 0.9), eligibility_rate = c(1, 0.6, 0.7))
  )
  for (nm in names(frames)) {
    fr <- frames[[nm]]
    x <- n_alloc(fr, n = 400, min_n_stratum = 90, fpc = "none")
    for (o in overrides) {
      nf <- fr
      for (rate in names(o)) nf[[rate]] <- rep_len(o[[rate]], nrow(fr))
      p <- prec_alloc(nf, n = x$detail$n, fpc = "none")
      expect_equal(do.call(effective_n, c(list(x), o)),
                   pop_var_alloc(nf) / p$se^2,
                   label = paste(nm, paste(names(o), collapse = "+")))
    }
  }
})

test_that("a cluster plan's override is the unit rate, stage response kept", {
  pc <- prec_cluster(n = c(39, 10), icc = 0.2)
  at_half <- prec_cluster(n = c(39, 10), icc = 0.2, resp_rate = 0.5)
  expect_equal(effective_n(pc, resp_rate = 0.5), 1 / at_half$cv^2)
  expect_equal(effective_n(pc, resp_rate = 0.5), 390 * 0.5 / 1.8)

  stored <- prec_cluster(n = c(20, 4, 8), icc = c(0.05, 0.1),
                         resp_rate = 0.9, resp_rate_psu = 0.85,
                         eligibility_rate = 0.8)
  moved <- prec_cluster(n = c(20, 4, 8), icc = c(0.05, 0.1),
                        resp_rate = 0.6, resp_rate_psu = 0.85,
                        eligibility_rate = 0.7)
  expect_equal(effective_n(stored, resp_rate = 0.6, eligibility_rate = 0.7),
               1 / moved$cv^2)

  plan <- n_cluster(stage_cost = c(500, 50), icc = 0.1, cv = 0.05,
                    resp_rate_psu = 0.9)
  fielded <- prec_cluster(n = plan$n, icc = 0.1, resp_rate = 0.5,
                          resp_rate_psu = 0.9)
  expect_equal(effective_n(plan, resp_rate = 0.5), 1 / fielded$cv^2)
  expect_error(effective_n(plan, resp_rate = 1.5), "resp_rate")
})

test_that("a generalized precision allocation is refused before any rate is read", {
  fr <- data.frame(stratum = c("A", "B"), N = c(1000, 2000))
  x <- prec_alloc(fr, n = c(50, 60),
                  measures = data.frame(stratum = c("A", "B"), name = "y",
                                        p = c(0.3, 0.4)),
                  targets = data.frame(name = "y", cv = 0.1))
  expect_error(effective_n(x), "generalized allocation")
  expect_error(effective_n(x, resp_rate = 0.5), "generalized allocation")
})

test_that("a size-only allocation's effective size reads its whole counts", {
  fr <- data.frame(stratum = c("a", "b"), N = c(1000, 3000),
                   resp_rate = c(0.4, 0.9))
  x <- n_alloc(fr, n = 10, alloc = "proportional")
  expect_equal(x$detail$n_int, c(4, 6))
  m <- c(4, 6) * c(0.4, 0.9)
  W <- c(0.25, 0.75)
  expect_equal(effective_n(x), 1 / sum(W^2 / m))
  expect_equal(effective_n(x), 6.981818, tolerance = 1e-6)
  short <- n_alloc(fr, n = 10.4, alloc = "proportional")
  expect_equal(sum(short$detail$n_int), 10)
  expect_equal(effective_n(short), effective_n(x))
})

test_that("an eligibility override derives a p-based sd from the new scope", {
  fr <- data.frame(stratum = c("a", "b"), N = c(100, 1000),
                   p = c(0.5, 0.2), eligibility_rate = c(0.02, 1))
  x <- n_alloc(fr, n = 20, min_n_stratum = 10)
  expect_equal(x$detail$sd[1], sqrt(0.5))
  fresh <- prec_alloc(transform(fr, eligibility_rate = 1), n = x$detail$n)
  expect_equal(effective_n(x, eligibility_rate = 1), effective_n(fresh))
  expect_equal(effective_n(x, eligibility_rate = 1), 13.08992,
               tolerance = 1e-6)
  expect_error(effective_n(x, eligibility_rate = 0.005),
               "greater than 1 in every stratum")

  by_hand <- function(S2_h, n_h) {
    W <- c(100, 1000) / 1100
    pbar <- sum(W * c(0.5, 0.2))
    (sum(W * S2_h) + sum(W * (c(0.5, 0.2) - pbar)^2)) / sum(W^2 * S2_h / n_h)
  }
  pq <- c(0.5, 0.2) * c(0.5, 0.8)
  expect_equal(effective_n(x, eligibility_rate = 1),
               by_hand(pq * c(100 / 99, 1000 / 999), x$detail$n))

  given <- transform(fr, sd = c(0.4, 0.4))
  xs <- n_alloc(given, n = 20, min_n_stratum = 10)
  expect_equal(xs$detail$sd, c(0.4, 0.4))
  W <- c(100, 1000) / 1100
  expect_equal(
    effective_n(xs, eligibility_rate = 1),
    (0.16 + sum(W * (c(0.5, 0.2) - sum(W * c(0.5, 0.2)))^2)) /
      sum(W^2 * 0.16 / xs$detail$n)
  )
  none <- n_alloc(fr, n = 20, min_n_stratum = 10, fpc = "none")
  expect_equal(effective_n(none, eligibility_rate = 1),
               by_hand(pq, none$detail$n))
})
