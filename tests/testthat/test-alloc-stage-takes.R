## Free takes of a classic cluster allocation

take_product_oracle <- function(cost, icc, r, M, cap, lo = 1) {
  S <- length(cost)
  objective <- function(t) {
    eff <- t * r[-1L]
    Q <- rev(cumprod(rev(c(eff, 1))))
    R <- numeric(S)
    keep <- 1
    for (s in seq_len(S - 1L)) {
      R[s] <- keep * icc[s]
      keep <- keep * (1 - icc[s])
    }
    R[S] <- keep
    f <- c(0, pmin(1, eff / M))
    inflation <- R[1L] * Q[1L] + sum((1 - f[-1L]) * R[-1L] * Q[-1L])
    psu_cost <- sum(cost * c(1, cumprod(t)))
    psu_cost * inflation / Q[1L]
  }
  if (S == 2L) {
    o <- stats::optimize(objective, c(lo, cap), tol = 1e-10)
    return(o$minimum)
  }
  best <- NULL
  for (s1 in exp(seq(log(lo), log(cap[1L]), length.out = 7L))) {
    for (s2 in exp(seq(log(lo), log(cap[2L]), length.out = 7L))) {
      o <- stats::optim(
        log(c(s1, s2)), function(z) objective(exp(z)),
        method = "L-BFGS-B", lower = log(c(lo, lo)), upper = log(cap)
      )
      if (is.null(best) || o$value < best$value) best <- o
    }
  }
  exp(best$par)
}

test_that("the free take under the stage correction counts eligible units", {
  frame <- data.frame(stratum = c("U", "R"), N = c(10000, 20000),
                      N_1 = c(200, 400), sd = c(10, 12), mean = c(50, 45),
                      icc_1 = c(0.05, 0.08), cost_1 = c(100, 200),
                      cost_2 = c(20, 25))
  for (yield in list(c(0.25, 1), c(1, 0.5), c(0.5, 0.5))) {
    fit <- suppressWarnings(n_alloc(
      frame, cv = 0.05, alloc = "optimal", fpc = "stage",
      eligibility_rate = yield[1L], eligible_per_unit = yield[2L]
    ))
    y <- prod(yield)
    for (h in 1:2) {
      oracle <- take_product_oracle(
        cost = c(frame$cost_1[h], frame$cost_2[h]), icc = frame$icc_1[h],
        r = c(1, y), M = frame$N[h] * y / frame$N_1[h],
        cap = frame$N[h] / frame$N_1[h]
      )
      expect_equal(fit$detail$take_2[h], oracle, tolerance = 1e-6,
                   info = paste(yield, collapse = "/"))
    }
  }
})

test_that("three-stage free takes under the stage correction count eligible units", {
  frame <- data.frame(stratum = "A", N = 40000, N_1 = 100, N_2 = 2000,
                      sd = 10, mean = 50, icc_1 = 0.05, icc_2 = 0.1,
                      cost_1 = 300, cost_2 = 40, cost_3 = 6)
  fit <- suppressWarnings(n_alloc(frame, cv = 0.05, alloc = "optimal",
                                  fpc = "stage", eligibility_rate = 0.25))
  oracle <- take_product_oracle(
    cost = c(300, 40, 6), icc = c(0.05, 0.1), r = c(1, 1, 0.25),
    M = c(2000 / 100, 40000 * 0.25 / 2000), cap = c(2000 / 100, 40000 / 2000)
  )
  expect_equal(unname(unlist(fit$detail[c("take_2", "take_3")])), oracle,
               tolerance = 1e-4)
})

test_that("a free two-stage take keeps one expected response", {
  frame <- data.frame(stratum = c("A", "B"), N = c(6000, 9000), sd = c(10, 14),
                      mean = c(40, 50), icc_1 = c(0.3, 0.1),
                      cost_1 = c(20, 600), cost_2 = c(100, 40))
  expect_warning(
    fit <- n_alloc(frame, cv = 0.03, alloc = "optimal", resp_rate = 0.5),
    "below 1 expected response.*'A'"
  )
  expect_equal(fit$detail$take_2[1], 2)
  expect_gte(fit$detail$take_2_int[1], 2)
  expect_gt(fit$detail$take_2[2], 2)
  expect_lte(fit$operational$cv, 0.03 + 1e-8)
})

test_that("free takes at three stages keep one expected response, as in n_cluster()", {
  frame <- data.frame(stratum = "A", N = 1e6, N_1 = 1e4, N_2 = 1e5, mean = 50,
                      sd = 10, icc_1 = 0.4, icc_2 = 0.5, cost_1 = 10,
                      cost_2 = 1, cost_3 = 100)
  expect_warning(
    fit <- n_alloc(frame, cv = 0.05, alloc = "optimal", resp_rate = 0.2),
    "take_3 is below 1 expected response"
  )
  ref <- suppressWarnings(n_cluster(
    stage_cost = c(10, 1, 100), icc = c(0.4, 0.5), cv = 0.05,
    unit_relvar = (10 / 50)^2, resp_rate = 0.2
  ))
  expect_equal(unname(unlist(fit$detail[c("take_2", "take_3")])),
               unname(ref$n[2:3]), tolerance = 1e-6)
  expect_gte(fit$detail$take_3_int, 5)
  expect_gte(fit$detail$take_2_int, 1)
})

test_that("a three-stage take rises past the PSU universe as a two-stage one does", {
  frame <- data.frame(stratum = c("A", "B"), N = c(6000, 9000), sd = c(10, 14),
                      mean = c(40, 50), icc_1 = c(0.05, 0.1), icc_2 = 0.1,
                      cost_1 = c(400, 600), cost_2 = c(50, 40), cost_3 = 5,
                      N_1 = c(3, 8))
  fit <- n_alloc(frame, cv = 0.025, alloc = "optimal")
  expect_lte(fit$operational$cv, 0.025 + 1e-8)
  expect_true(all(fit$detail$n_1 <= frame$N_1 + 1e-9))
  expect_true(all(fit$detail$n_1_int <= frame$N_1))
  risen <- fit$detail$n_1 >= frame$N_1 - 1e-9
  expect_true(any(risen))
  expect_equal(fit$detail$n[risen],
               frame$N_1[risen] * fit$detail$take_2[risen] *
                 fit$detail$take_3[risen], tolerance = 1e-9)

  two <- frame[, setdiff(names(frame), c("icc_2", "cost_3"))]
  fit_2 <- n_alloc(two, cv = 0.025, alloc = "optimal")
  expect_lte(fit_2$operational$cv, 0.025 + 1e-8)
})

## A response floor the stage size cannot hold

test_that("a free take whose floor exceeds its stage size is refused", {
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, N_2 = 20, mean = 10,
                  sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 10, cost_2 = 1,
                  cost_3 = 100)
  for (fpc in c("none", "unit", "stage")) {
    expect_error(
      n_alloc(f, n = 20, alloc = "optimal", fpc = fpc, resp_rate = .2),
      "stage 3 of stratum 'A' holds 2 units per parent unit", info = fpc
    )
  }
  f25 <- transform(f, N = 20, N_1 = 4, N_2 = 8)
  expect_error(
    n_alloc(f25, n = 8, alloc = "optimal", fpc = "stage", resp_rate = .4),
    "holds 2.5 units per parent unit, which cannot deliver one expected response in whole units"
  )
  at_cap <- suppressWarnings(
    n_alloc(f, n = 20, alloc = "optimal", fpc = "none", resp_rate = .5)
  )
  expect_identical(at_cap$detail$take_3_int, 2L)
  fixed <- n_alloc(transform(f, take_2 = 1, take_3 = 2), n = 20, fpc = "none",
                   resp_rate = .2)
  expect_identical(fixed$detail$take_3_int, 2L)
})

test_that("a two-stage whole take stays within the PSU size under the stage correction", {
  f2 <- data.frame(stratum = "A", N = 20, N_1 = 10, mean = 10, sd = 2,
                   icc_1 = .1, cost_1 = 1000, cost_2 = 1)
  fit <- suppressWarnings(n_alloc(f2, n = 10, alloc = "optimal",
                                  fpc = "stage"))
  expect_identical(c(fit$detail$n_1_int, fit$detail$take_2_int), c(5L, 2L))
  expect_equal(fit$operational$cost, 5010)
  for (budget in c(3000, 6000)) {
    fb <- suppressWarnings(n_alloc(f2, budget = budget, alloc = "optimal",
                                   fpc = "stage"))
    expect_lte(fb$detail$take_2_int, 2L)
    expect_lte(fb$operational$cost, budget)
  }
  expect_error(
    n_alloc(f2, n = 10, alloc = "optimal", fpc = "stage", resp_rate = .2),
    "stage 2 of stratum 'A' holds 2 units per parent unit"
  )
})

## A fixed take above its stage size under the stage correction

test_that("a fixed take above its stage size is refused under the stage correction", {
  f2 <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                   icc_1 = .1, cost_1 = 1000, cost_2 = 1)
  f3 <- data.frame(stratum = "A", N = 40, N_1 = 10, N_2 = 20, mean = 10,
                   sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                   cost_3 = 1)
  ## N / N_1 = 4 units per PSU at two stages, N_2 / N_1 = 2 SSUs at three.
  expect_error(n_alloc(transform(f2, take_2 = 8), n = 16, alloc = "optimal",
                       fpc = "stage"),
               "take_2 = 8 in stratum 'A' exceeds the 4 units")
  expect_error(n_alloc(transform(f3, take_2 = 4, take_3 = 2), n = 16,
                       alloc = "optimal", fpc = "stage"),
               "take_2 = 4 in stratum 'A' exceeds the 2 units")
  at_cap <- suppressWarnings(n_alloc(transform(f2, take_2 = 4), n = 16,
                                     alloc = "optimal", fpc = "stage"))
  expect_equal(c(at_cap$detail$take_2, at_cap$detail$take_2_int), c(4, 4))
})

## Stage sizes cap every take under every correction

test_that("stage sizes cap continuous and whole takes under every correction", {
  f2 <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                   icc_1 = .1, cost_1 = 1000, cost_2 = 1)
  f3 <- transform(f2, N_2 = 20, icc_2 = .1, cost_3 = 1)
  ## 4 units per PSU at two stages, 2 SSUs of 2 units at three.
  caps <- list(f2 = 4, f3 = c(2, 2))
  checked <- character(0)
  for (nm in names(caps)) {
    fr <- list(f2 = f2, f3 = f3)[[nm]]
    S <- length(caps[[nm]]) + 1L
    for (fpc in c("none", "unit", "stage")) {
      for (target in list(list(n = 20), list(cv = .09), list(budget = 1020))) {
        id <- paste(nm, fpc, names(target))
        z <- suppressWarnings(do.call(n_alloc, c(list(fr, alloc = "optimal",
                                                      fpc = fpc), target)))
        cont <- unlist(z$detail[paste0("take_", 2:S)])
        whole <- unlist(z$detail[paste0("take_", 2:S, "_int")])
        expect_true(all(cont <= caps[[nm]] + 1e-9), info = id)
        expect_true(all(whole <= caps[[nm]]), info = id)
        expect_lte(z$detail$n_1, 10 + 1e-9)
        expect_gt(z$operational$cv, 0)
        expect_equal(suppressWarnings(prec_alloc(z))$cv, z$cv, info = id)
        checked <- c(checked, id)
      }
    }
  }
  expect_length(checked, 18L)
})

test_that("a total past what the stages hold is refused before allocation", {
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, N_2 = 20, mean = 10,
                  sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1, take_2 = 1)
  ## With one SSU per PSU, 10 PSUs of 1 SSU of 2 units hold 20.
  for (fpc in c("none", "unit", "stage")) {
    full <- suppressWarnings(n_alloc(f, n = 20, alloc = "optimal", fpc = fpc))
    expect_equal(c(full$detail$n_1, full$detail$take_3), c(10, 2), info = fpc)
    expect_equal(full$detail$n_int, 20, info = fpc)
    expect_error(suppressWarnings(n_alloc(f, n = 30, alloc = "optimal",
                                          fpc = fpc)),
                 "ceiling comes from the stage sizes", info = fpc)
  }
  ## At four stages 1.6 SSUs and 2 units fit continuously, 10 * 1.6 * 2 =
  ## 32, and the whole stages hold 10 * 1 * 2 = 20.
  f4 <- transform(f, N_3 = 32, N = 64, icc_3 = .1, cost_4 = 1)
  fit <- suppressWarnings(n_alloc(f4, n = 20, alloc = "optimal",
                                  fpc = "stage"))
  expect_true(all(c(fit$detail$take_3, fit$detail$take_4) <=
                    c(1.6, 2) + 1e-9))
  expect_error(suppressWarnings(n_alloc(f4, n = 32, alloc = "optimal",
                                        fpc = "stage")),
               "feasible totals: 1 to 20")
  expect_error(suppressWarnings(n_alloc(f4, n = 50, alloc = "optimal",
                                        fpc = "stage")),
               "ceiling comes from the stage sizes")
})

test_that("fixed takes, fractional caps and missing sizes meet the capacity rule", {
  base <- data.frame(stratum = "A", N = 40, N_1 = 10, N_2 = 20, mean = 10,
                     sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                     cost_2 = 1, cost_3 = 1)
  for (fpc in c("none", "unit")) {
    expect_error(n_alloc(transform(base, take_2 = 4, take_3 = 2), n = 16,
                         alloc = "optimal", fpc = fpc),
                 "take_2 = 4 in stratum 'A' exceeds the 2 units", info = fpc)
  }
  half <- suppressWarnings(n_alloc(transform(base, take_2 = 1.5), n = 20,
                                   alloc = "optimal", fpc = "none"))
  expect_equal(c(half$detail$take_2, half$detail$take_2_int), c(1.5, 2))
  ## 2.5 units per SSU: one expected response at rate .4 needs 2.5, which
  ## whole units cannot give, and at rate .5 needs 2.
  wide <- transform(base, N = 50)
  expect_error(n_alloc(wide, n = 20, alloc = "optimal", fpc = "none",
                       resp_rate = .4), "in whole units")
  ok <- suppressWarnings(n_alloc(wide, n = 20, alloc = "optimal",
                                 fpc = "none", resp_rate = .5))
  expect_equal(c(ok$detail$take_3, ok$detail$take_3_int), c(2.5, 2))
  ## Without N_2 neither ratio is known, so neither take is capped.
  open <- base
  open$N_2 <- NULL
  wide_open <- suppressWarnings(n_alloc(open, n = 20, alloc = "optimal",
                                        fpc = "unit"))
  expect_gt(wide_open$detail$take_2, 2)
})

test_that("whole capacity clips a fixed fractional take before multiplying", {
  f <- data.frame(stratum = "A", N = 52, N_1 = 10, N_2 = 26, mean = 10,
                  sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                  cost_2 = 1, cost_3 = 1, take_2 = 2.6)
  for (fpc in c("none", "unit", "stage")) {
    for (take_3 in c(NA_real_, 2)) {
      frame <- transform(f, take_3 = take_3)
      fit <- n_alloc(frame, n = 40, alloc = "optimal", fpc = fpc)
      expect_equal(fit$detail$n_int, 40)
      expect_equal(c(fit$detail$take_2_int, fit$detail$take_3_int), c(2, 2))
      expect_error(n_alloc(frame, n = 50, alloc = "optimal", fpc = fpc),
                   "feasible totals: 1 to 40")
    }
  }
})

test_that("every PSU of the frame keeps a unit when a stage size is missing", {
  ## 10 PSUs hold 40 units, so a sample of a PSUs holds at most 30 + a.
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  room <- function(fit) 30 + fit$detail$n_1_int
  checked <- character(0)
  for (fpc in c("unit", "none")) {
    for (target in list(list(n = 20), list(cv = .05), list(budget = 1040))) {
      id <- paste(fpc, names(target))
      fit <- suppressWarnings(do.call(n_alloc, c(list(f, alloc = "optimal",
                                                      fpc = fpc), target)))
      expect_lte(fit$detail$n_int, room(fit))
      expect_lte(fit$detail$n, 30 + fit$detail$n_1 + 1e-9)
      if (fit$detail$n_1_int < 10) expect_gt(fit$operational$cv, 0)
      checked <- c(checked, id)
    }
  }
  expect_length(checked, 6L)

  partial <- transform(f, take_2 = 5, take_3 = 4)
  for (fpc in c("unit", "none")) {
    expect_error(n_alloc(partial, n = 40, fpc = fpc),
                 "every other PSU of the frame holds one", info = fpc)
  }
  expect_equal(suppressWarnings(n_alloc(partial, n = 20, fpc = "unit"))$detail$n_int,
               20)
  expect_error(n_alloc(transform(f, take_2 = 8, take_3 = 4), n = 20),
               "the units one PSU holds while every other PSU holds one")
  msg <- "more units in its PSUs than the frame leaves them"
  for (args in list(list(fpc = "unit"), list(fpc = "none"),
                    list(fpc = "unit", resp_rate = .9))) {
    expect_warning(do.call(prec_alloc, c(list(partial, n = 40), args)), msg)
  }
  every <- transform(partial, take_2 = 2, take_3 = 2)
  expect_no_warning(n_alloc(every, n = 40, fpc = "unit"))
  expect_no_warning(prec_alloc(every, n = 40, fpc = "unit"))
})

test_that("searched designs at a missing stage size match enumeration under the room bound", {
  ## Oracle: every whole design that leaves each other PSU a unit, priced
  ## from the stage costs and assessed by prec_alloc() at its fixed takes.
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  grid <- expand.grid(a = 1:10, b = 1:31, c = 1:31)
  grid <- grid[grid$a * grid$b * grid$c <= 30 + grid$a, ]
  grid$cost <- with(grid, a * 1000 + a * b + a * b * c)
  checked <- character(0)
  for (fpc in c("unit", "none")) {
    cv <- vapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      suppressWarnings(prec_alloc(transform(f, take_2 = g$b, take_3 = g$c),
                                  n = g$a * g$b * g$c, fpc = fpc))$cv
    }, 1)
    fit <- suppressWarnings(n_alloc(f, budget = 1040, alloc = "optimal",
                                    fpc = fpc))
    expect_equal(fit$operational$cv, min(cv[grid$cost <= 1040]),
                 tolerance = 1e-9, info = fpc)
    fit <- suppressWarnings(n_alloc(f, cv = .05, alloc = "optimal", fpc = fpc))
    expect_equal(fit$operational$cost, min(grid$cost[cv <= .05 * (1 + 1e-9)]),
                 info = fpc)
    checked <- c(checked, fpc)
  }
  expect_identical(checked, c("unit", "none"))
})

## A fixed total reached exactly where whole takes allow it

test_that("a fixed total is fielded exactly whenever some whole design reaches it", {
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  ## Oracle: a PSUs of b SSUs of c units, within the 10 PSUs and the room
  ## that leaves each other PSU a unit.
  grid <- expand.grid(a = 1:10, b = 1:31, c = 1:31)
  grid <- grid[grid$a * (grid$b * grid$c - 1) <= 30, ]
  reach <- sort(unique(with(grid, a * b * c)))
  checked <- character(0)
  for (fpc in c("unit", "none")) {
    for (n in 2:40) {
      fit <- suppressWarnings(n_alloc(f, n = n, alloc = "optimal", fpc = fpc))
      d <- fit$detail
      expect_lte(d$n_int, 30 + d$n_1_int)
      if (n %in% reach) expect_equal(d$n_int, n, info = paste(fpc, n))
      checked <- c(checked, paste(fpc, n))
    }
  }
  expect_length(checked, 78L)
  ## 37 is prime and past one PSU's 31, and 34, 38 and 39 need two or
  ## three PSUs that hold more than the room.
  expect_identical(setdiff(2:40, reach), c(34L, 37L, 38L, 39L))
  ## 40 needs every PSU, so the takes fall to 4 units per PSU.
  census <- suppressWarnings(n_alloc(f, n = 40, alloc = "optimal",
                                     fpc = "none"))
  expect_equal(census$detail$n_1, 10)
  expect_equal(census$detail$take_2 * census$detail$take_3, 4)
  ## With N_2 = 20 the takes are 2 and 2, so 30 is out of reach and the
  ## rounding keeps its nearest total.
  full <- suppressWarnings(n_alloc(transform(f, N_2 = 20), n = 30,
                                   alloc = "optimal", fpc = "stage"))
  expect_equal(c(full$detail$take_2_int, full$detail$take_3_int), c(2, 2))
  expect_equal(full$detail$n_int, 32)
})

test_that("an exact fixed total keeps the strata that already met theirs", {
  two <- data.frame(stratum = c("A", "B"), N = c(40, 400), N_1 = c(10, 50),
                    mean = 10, sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                    cost_2 = 1, cost_3 = 1)
  for (n in c(25, 60, 101)) {
    fit <- suppressWarnings(n_alloc(two, n = n, alloc = "optimal",
                                    fpc = "none"))
    expect_equal(sum(fit$detail$n_int), n, info = n)
    d <- fit$detail
    expect_true(all(d$n_int <= d$N - d$N_1 + d$n_1_int), info = n)
  }
})

## Capacity at the take product that holds most, and the fixed-size repair

test_that("a response floor above the mean PSU size prices a reachable census", {
  ## At rate .2 a final take holds at least 5 units, more than the 4 per PSU
  ## of the frame. Oracle: whole designs within the room, priced and
  ## assessed at fixed takes.
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  grid <- expand.grid(a = 1:10, b = 1:31, c = 5:31)
  grid <- grid[grid$a * (grid$b * grid$c - 1) <= 30, ]
  grid$cost <- with(grid, a * 1000 + a * b + a * b * c)
  grid <- grid[grid$cost <= 1040, ]
  for (fpc in c("unit", "none")) {
    cv <- vapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      suppressWarnings(prec_alloc(transform(f, take_2 = g$b, take_3 = g$c),
                                  n = g$a * g$b * g$c, fpc = fpc,
                                  resp_rate = .2))$cv
    }, 1)
    fit <- suppressWarnings(n_alloc(f, budget = 1040, alloc = "optimal",
                                    fpc = fpc, resp_rate = .2))
    expect_equal(fit$operational$cv, min(cv), tolerance = 1e-9, info = fpc)
    expect_lte(fit$operational$cost, 1040)
  }
})

test_that("a smaller take product counts toward whole capacity when it holds more", {
  ## take_2 may be 1 to 3 SSUs. With 4 x 2 fixed below, 1 SSU leaves room for
  ## 4 PSUs of 8 (32), 2 SSUs for 2 PSUs of 16 (32), 3 SSUs for 1 of 24.
  f4 <- data.frame(stratum = "A", N = 40, N_1 = 10, N_2 = 30, mean = 10,
                   sd = 2, icc_1 = .1, icc_2 = .1, icc_3 = .1, cost_1 = 1000,
                   cost_2 = 1, cost_3 = 1, cost_4 = 1, take_3 = 4, take_4 = 2)
  b <- 1:3
  cap <- max(pmin(10, floor(30 / (8 * b - 1))) * 8 * b)
  expect_equal(cap, 32)
  fit <- suppressWarnings(n_alloc(f4, n = 32, alloc = "optimal", fpc = "none"))
  expect_equal(fit$detail$n_int, 32)
  expect_equal(c(fit$detail$take_3_int, fit$detail$take_4_int), c(4, 2))
  expect_error(suppressWarnings(n_alloc(f4, n = 33, alloc = "optimal",
                                        fpc = "none")),
               "feasible totals: 1 to 32|maximum feasible total")
})

test_that("fixed takes beyond one PSU's room are refused, and the field keeps them", {
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  expect_error(n_alloc(transform(f, take_2 = 32), n = 8, alloc = "optimal",
                       fpc = "none"),
               "need 32 units per PSU, more than one PSU holds")
  ## 31 fits with a final take of one, which the clamp may not go below.
  fit <- suppressWarnings(n_alloc(transform(f, take_2 = 31), n = 31,
                                  alloc = "optimal", fpc = "none"))
  expect_equal(c(fit$detail$take_2, fit$detail$take_3), c(31, 1))
  expect_equal(c(fit$detail$take_2_int, fit$detail$take_3_int), c(31, 1))
  expect_equal(fit$detail$n_int, 31)
})

test_that("the fixed-size repair runs when rounding fails and moves units across strata", {
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  floored <- suppressWarnings(n_alloc(f, n = 25, alloc = "optimal",
                                      fpc = "none", min_n_stratum = 25))
  expect_equal(floored$detail$n_int, 25)
  ## Oracle: the whole totals each stratum reaches within the room, and the
  ## pairs of them that make n.
  grid <- expand.grid(a = 1:10, b = 1:31, c = 1:31)
  grid <- grid[grid$a * (grid$b * grid$c - 1) <= 30, ]
  reach <- sort(unique(with(grid, a * b * c)))
  two <- rbind(transform(f, stratum = "A"), transform(f, stratum = "B"))
  fit <- suppressWarnings(n_alloc(two, n = 68, alloc = "optimal",
                                  fpc = "none"))
  expect_equal(sum(fit$detail$n_int), 68)
  expect_true(all(fit$detail$n_int %in% reach))
  expect_false(34 %in% reach)
  ## Far from the continuous design an exact total is not worth taking: a
  ## prime total of 2000003 would need one PSU of two million.
  big <- data.frame(stratum = "A", N = 4e6, N_1 = 1000, mean = 10, sd = 2,
                    icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                    cost_3 = 1)
  prime <- suppressWarnings(n_alloc(big, n = 2000003, alloc = "optimal",
                                    fpc = "none"))
  expect_equal(prime$detail$n_1_int, 1000)
  expect_lte(abs(prime$detail$n_int - 2000003), 50)
})

## Whole capacity and the fixed-size selection across strata and domains

test_that("whole capacity matches every take product the stage bounds allow", {
  ## Oracle: each product of whole takes within the stage intervals, at the
  ## most PSUs that leave every other PSU of the frame a unit.
  brute <- function(N, A, lo, hi) {
    g <- as.matrix(expand.grid(lapply(seq_along(lo), function(k) lo[k]:hi[k])))
    P <- unique(apply(g, 1L, prod))
    a <- vapply(P, function(p) sum(seq_len(A) * (p - 1) <= N - A), 1)
    max(a * P)
  }
  set.seed(20261122)
  multi <- 0L
  interrupted <- 0L
  for (i in 1:300) {
    K <- sample(1:3, 1)
    A <- sample(1:12, 1)
    N <- A + sample(0:60, 1) + if (i %% 5 == 0) .5 else 0
    lo <- sample(1:4, K, TRUE)
    hi <- lo + sample(0:8, K, TRUE)
    if (i %% 4 == 0) hi[1L] <- Inf
    want <- brute(N, A, lo, pmin(hi, N + 1))
    expect_identical(svyplan:::.ms_whole_capacity(N, A, lo, hi), want,
                     info = i)
    multi <- multi + (sum(lo < hi) >= 2L)
    cut <- svyplan:::.ms_whole_capacity(N, A, lo, hi, max_boxes = 0L)
    if (isFALSE(attr(cut, "exact"))) {
      interrupted <- interrupted + 1L
      expect_gte(as.numeric(cut), want)
    } else {
      expect_identical(cut, want, info = i)
    }
  }
  expect_gt(multi, 80L)
  expect_gt(interrupted, 20L)
})

test_that("a large fixed take prices its capacity without listing the products", {
  f <- data.frame(stratum = "A", N = 30e6, N_1 = 2, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1, take_2 = 10e6)
  setTimeLimit(elapsed = 60, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  ## Two PSUs of ten million hold 20 million, and one PSU of 20 million
  ## leaves the other a unit.
  expect_identical(
    svyplan:::.ms_whole_capacity(3e8, 2, c(1e8, 1), c(1e8, Inf)), 2e8
  )
  fit <- suppressWarnings(n_alloc(f, n = 20e6, alloc = "optimal",
                                  fpc = "none"))
  expect_equal(c(fit$detail$n_1_int, fit$detail$take_2_int,
                 fit$detail$take_3_int), c(2, 10e6, 1))
  expect_error(suppressWarnings(n_alloc(f, n = 20e6 + 1, alloc = "optimal",
                                        fpc = "none")),
               "feasible totals: 1 to 20000000")
})

test_that("a stratum with every take fixed offers its own totals to the fixed-size repair", {
  f <- data.frame(stratum = c("A", "B"), N = c(40, 4), N_1 = c(10, 1),
                  mean = 10, sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                  cost_2 = 1, cost_3 = 1, take_2 = c(NA, 2), take_3 = c(NA, 2))
  for (fpc in c("none", "unit")) {
    fit <- suppressWarnings(n_alloc(f, n = 29, alloc = "optimal", fpc = fpc,
                                    min_n_stratum = c(25, 4)))
    d <- fit$detail
    expect_equal(d$n_int, c(25, 4), info = fpc)
    expect_equal(c(d$n_1_int[2], d$take_2_int[2], d$take_3_int[2]),
                 c(1, 2, 2), info = fpc)
  }
})

test_that("domain floors are met while the strata totals are chosen", {
  ## Each stratum alone reaches 33 and 35 but not 34 (three PSUs of 11 and
  ## five of 7 leave every other PSU a unit, 34 needs more).
  f <- data.frame(stratum = c("A", "B"), N = 40, N_1 = 10, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 1,
                  cost_3 = 1)
  f$d <- f$stratum
  grid <- expand.grid(a = 1:10, b = 1:31, c = 1:31)
  reach <- sort(unique(with(grid[grid$a * (grid$b * grid$c - 1) <= 30, ],
                            a * b * c)))
  expect_true(all(c(33, 35) %in% reach))
  expect_false(34 %in% reach)
  checked <- character(0)
  for (lv in c("A", "B")) {
    for (rows in list(1:2, 2:1)) {
      fit <- suppressWarnings(n_alloc(
        f[rows, ], n = 68, domains = "d", alloc = "optimal", fpc = "none",
        min_n_domain = data.frame(domain = "d", level = lv, min_n = 34)
      ))
      d <- fit$detail
      id <- paste(lv, rows[1])
      expect_equal(sum(d$n_int), 68, info = id)
      expect_equal(d$n_int[d$stratum == lv], 35, info = id)
      expect_true(all(d$n_int %in% reach), info = id)
      checked <- c(checked, id)
    }
  }
  expect_setequal(checked, c("A 1", "A 2", "B 1", "B 2"))
  ## 35 needs A at 35 or more and B at one at least, so 36 comes closest. 79
  ## needs a total of 39, 37 or 38 in one stratum, so 80 comes closest.
  floor_a <- data.frame(domain = "d", level = "A", min_n = 34)
  for (n in c(35, 79)) {
    fit <- suppressWarnings(n_alloc(f, n = n, domains = "d", alloc = "optimal",
                                    fpc = "none", min_n_domain = floor_a))
    d <- fit$detail
    expect_equal(sum(d$n_int), n + 1, info = n)
    expect_gte(d$n_int[1], 34)
    expect_true(all(d$n_int %in% reach), info = n)
  }
})

test_that("below one continuous PSU the exact design may take more PSUs", {
  ## N_2 = 20 caps both takes at 2, so 3 units need three PSUs of one.
  f <- data.frame(stratum = "A", N = 40, N_1 = 10, N_2 = 20, mean = 10,
                  sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000, cost_2 = 10,
                  cost_3 = 1)
  for (fpc in c("none", "unit", "stage")) {
    fit <- suppressWarnings(n_alloc(f, n = 3, alloc = "optimal", fpc = fpc))
    d <- fit$detail
    expect_lt(d$n_1, 1)
    expect_equal(c(d$n_int, d$n_1_int, d$take_2_int, d$take_3_int),
                 c(3, 3, 1, 1), info = fpc)
  }
})

test_that("fixed-size designs over random small frames stay within every bound", {
  ## Oracle: every whole design of each stratum within its adjacent stage
  ## caps, its fixed takes, N_1, its floor and the room that leaves every
  ## other PSU a unit. The returned design must be one of them.
  legal <- function(N, A, N_2, t2, t3, lo) {
    c2 <- if (is.na(N_2)) N else floor(N_2 / A)
    c3 <- if (is.na(N_2)) N else floor(N / N_2)
    g <- expand.grid(a = seq_len(A),
                     b = if (is.na(t2)) seq_len(c2) else t2,
                     c = if (is.na(t3)) seq_len(c3) else t3)
    g[g$a * (g$b * g$c - 1) <= N - A & g$a * g$b * g$c >= lo, ]
  }
  set.seed(4656)
  checked <- 0L
  exact <- 0L
  floored <- 0L
  for (i in 1:120) {
    H <- sample(2:3, 1)
    sized <- i %% 2 == 0
    fr <- do.call(rbind, lapply(seq_len(H), function(h) {
      N <- sample(8:40, 1)
      A <- sample(2:min(10, N - 1), 1)
      N_2 <- if (sized) as.numeric(sample(A:N, 1)) else NA_real_
      t2 <- if (runif(1) < .3) 2 else NA_real_
      t3 <- if (runif(1) < .3) 1 else NA_real_
      if (sized && !is.na(t2) && t2 > floor(N_2 / A)) t2 <- NA_real_
      data.frame(stratum = LETTERS[h], N = N, N_1 = A, N_2 = N_2, mean = 10,
                 sd = sample(1:4, 1), icc_1 = .1, icc_2 = .1,
                 cost_1 = sample(c(50, 1000), 1), cost_2 = sample(1:10, 1),
                 cost_3 = 1, take_2 = t2, take_3 = t3,
                 d = if (h == 1L) "x" else "y")
    }))
    if (!sized) fr$N_2 <- NULL
    lo <- sample(1:6, H, TRUE)
    sets <- lapply(seq_len(H), function(h) {
      legal(fr$N[h], fr$N_1[h], if (sized) fr$N_2[h] else NA_real_,
            fr$take_2[h], fr$take_3[h], lo[h])
    })
    if (any(vapply(sets, nrow, 1L) == 0L)) next
    top <- vapply(sets, function(g) max(g$a * g$b * g$c), 1)
    n <- sample(seq(sum(lo), sum(top)), 1)
    args <- list(fr, n = n, alloc = "optimal", fpc = "none",
                 min_n_stratum = lo)
    need <- 0
    if (i %% 3 == 0) {
      need <- sample(seq(lo[1], top[1]), 1)
      args$domains <- "d"
      args$min_n_domain <- data.frame(domain = "d", level = "x", min_n = need)
    }
    fit <- tryCatch(suppressWarnings(do.call(n_alloc, args)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    d <- fit$detail
    for (h in seq_len(H)) {
      g <- sets[[h]]
      expect_true(any(g$a == d$n_1_int[h] & g$b == d$take_2_int[h] &
                        g$c == d$take_3_int[h]), info = paste(i, h))
    }
    expect_equal(d$n_int, d$n_1_int * d$take_2_int * d$take_3_int)
    expect_gte(d$n_int[1], need)
    checked <- checked + 1L
    exact <- exact + (sum(d$n_int) == n)
    floored <- floored + (need > 0)
  }
  expect_gte(checked, 100L)
  expect_gte(exact, 90L)
  expect_gte(floored, 25L)
})

test_that("strata targets follow what whole takes reach, and their shares stay in reach", {
  ## Oracle: the whole totals each stratum reaches within its stage caps and
  ## the room, and the sums of one total per stratum.
  reach <- function(N, A, c2, c3, t2 = NA, lo = 1) {
    g <- expand.grid(a = seq_len(A), b = if (is.na(t2)) seq_len(c2) else t2,
                     c = seq_len(c3))
    g <- g[g$a * (g$b * g$c - 1) <= N - A, ]
    n <- unique(g$a * g$b * g$c)
    sort(n[n >= lo])
  }
  sums <- function(...) sort(unique(rowSums(expand.grid(list(...)))))
  ## B reaches 30 at most, so A takes the rest of 39, and 38 or 40 is closest.
  short <- data.frame(stratum = c("A", "B"), N = c(12, 35), N_1 = 2,
                      N_2 = c(12, 7), mean = 10, sd = 1, icc_1 = .1,
                      icc_2 = c(.15, .1), cost_1 = c(1000, 200),
                      cost_2 = c(10, 4), cost_3 = 1)
  totals <- sums(reach(12, 2, 6, 1), reach(35, 2, 3, 5))
  expect_false(39 %in% totals)
  expect_true(all(c(38, 40) %in% totals))
  fit <- suppressWarnings(n_alloc(short, n = 39, alloc = "optimal",
                                  fpc = "none"))
  expect_equal(abs(sum(fit$detail$n_int) - 39), 1)
  ## C holds 24 whole units, less than its share of 32.7. A then fields 10,
  ## near its own share of 9.6, for an exact 66.
  three <- data.frame(stratum = c("A", "B", "C"), N = c(14, 41, 35),
                      N_1 = c(7, 4, 5), mean = 10, sd = c(4, 1, 4),
                      icc_1 = c(.1, .1, .05), icc_2 = c(.1, .1, .15),
                      cost_1 = c(1000, 1000, 50), cost_2 = c(1, 6, 3),
                      cost_3 = 1, take_2 = c(2, NA, 4), take_3 = c(NA, NA, 3))
  totals <- sums(reach(14, 7, 14, 14, t2 = 2), reach(41, 4, 41, 41, lo = 3),
                 c(12, 24))
  expect_true(66 %in% totals)
  fit <- suppressWarnings(n_alloc(three, n = 66, alloc = "optimal",
                                  fpc = "unit", min_n_stratum = c(1, 3, 3)))
  expect_equal(sum(fit$detail$n_int), 66)
  expect_equal(fit$detail$n_int[3], 24)
})

test_that("a floor the selected totals meet leaves the fixed-size design as chosen", {
  ## A reaches 2 to 4 and B reaches 3, 4, 5, 6, 8 or 10, so 2 + 5 fields 7.
  ## A floor of 2 on A holds already and must not move the design.
  f <- data.frame(stratum = c("A", "B"), N = c(6, 10), N_1 = c(4, 5),
                  mean = 10, sd = c(1, 2), icc_1 = c(.05, .1),
                  icc_2 = c(.1, .15), cost_1 = c(200, 50), cost_2 = c(1, 4),
                  cost_3 = 1, d = c("y", "x"))
  free <- suppressWarnings(n_alloc(f, n = 7, alloc = "optimal", fpc = "none",
                                   min_n_stratum = c(2, 3)))
  held <- suppressWarnings(n_alloc(
    f, n = 7, alloc = "optimal", fpc = "none", min_n_stratum = c(2, 3),
    domains = "d", min_n_domain = data.frame(domain = "d", level = "y",
                                             min_n = 2)
  ))
  expect_equal(free$detail$n_int, c(2, 5))
  expect_equal(held$detail$n_int, c(2, 5))
})

## The whole PSU room at an exact counting boundary

test_that("whole PSU room respects the integer boundary at large takes", {
  ## Oracle: a PSUs of P units fit when a * (P - 1) <= N - N_1, with every
  ## count whole and exact in doubles.
  room <- svyplan:::.psu_room
  P <- 1000000002
  at_two <- 3 + 2 * (P - 1)
  expect_equal(room(at_two + c(-1, 0, 1), 3, rep(P, 3)), c(1, 2, 2))
  expect_equal(room(10.5, 3, c(1, 2, 4, 10, Inf)), c(Inf, 7, 2, 0, Inf))
  expect_equal(room(10, NULL, c(1, 2)), c(Inf, Inf))
  set.seed(2061)
  for (i in 1:200) {
    N_1 <- sample(2:50, 1)
    P <- sample(c(2:40, 1e6 + sample(0:99, 1), 1e9 + sample(0:99, 1)), 1)
    a <- sample(1:N_1, 1)
    N <- N_1 + a * (P - 1) + sample(-1:1, 1)
    want <- floor((N - N_1) / (P - 1))
    if ((want + 1) * (P - 1) <= N - N_1) want <- want + 1
    if (want * (P - 1) > N - N_1) want <- want - 1
    got <- room(N, N_1, P)
    expect_lte(got * (P - 1), N - N_1)
    expect_gt((got + 1) * (P - 1), N - N_1)
    expect_equal(got, want, info = i)
  }
})

test_that("a large whole budget fit leaves a unit for each unsampled PSU", {
  P <- 1000000002
  f <- data.frame(stratum = "A", N = 2 * P, N_1 = 3, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                  cost_2 = 1, cost_3 = 1, take_2 = P, take_3 = 1)
  fit <- suppressWarnings(n_alloc(f, budget = 4000002008, alloc = "optimal",
                                  fpc = "unit"))
  d <- fit$detail
  expect_equal(d$n_1_int, 1)
  expect_equal(d$n_int, P)
  expect_lte(d$n_int, f$N - f$N_1 + d$n_1_int)
  expect_gt(fit$operational$cv, 0)
  ## One PSU of 2e12 leaves the other PSU empty, so two PSUs of 1e12 hold most.
  expect_identical(
    svyplan:::.ms_whole_capacity(3e12, 2, c(1e12, 1), c(1e12, Inf)), 2e12
  )
})

## A nearby window that misses the row bounds

test_that("a whole-size window disjoint from row bounds offers no targets", {
  ## A holds 2,000 whole units (takes of at most 2 and 1) against a
  ## continuous share near 2,200, and B fields even totals only, so the
  ## rounding misses 2401 and the nearby search runs.
  f <- data.frame(stratum = c("A", "B"), N = 2250, N_1 = 1000, N_2 = 2000,
                  mean = 10, sd = c(10, 1), icc_1 = .1, icc_2 = .1,
                  cost_1 = 1000, cost_2 = 1, cost_3 = 1, take_2 = c(NA, 2),
                  take_3 = c(NA, 1))
  original <- svyplan:::.ms_exact_takes
  targets_a <- numeric()
  local_mocked_bindings(
    .ms_exact_takes = function(prep, h, target, b_row, a_cont) {
      if (h == 1L) targets_a <<- c(targets_a, target)
      original(prep, h, target, b_row, a_cont)
    }
  )
  fit <- suppressWarnings(n_alloc(f, n = 2401, alloc = "optimal",
                                  fpc = "none"))
  expect_gt(fit$detail$n[1], 2150)
  expect_equal(fit$detail$n_int, c(2000, 400))
  expect_true(2000 %in% targets_a)
  expect_true(all(targets_a <= 2000))
})

test_that("a domain floor above what its strata hold in whole units is refused first", {
  ## Each stratum holds 1,000 whole units, one per PSU, since N_2 / N_1 and
  ## N / N_2 both round down to 1, while 2,250 fit continuously.
  f <- data.frame(stratum = c("A", "B"), d = c("A", "B"),
                  N = 2250, N_1 = 1000, N_2 = 1500, mean = 10, sd = 2,
                  icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                  cost_2 = 1, cost_3 = 1)
  floor_at <- function(x) data.frame(domain = "d", level = "A", min_n = x)
  for (target in list(list(n = 1900), list(cv = .05))) {
    expect_error(
      suppressWarnings(do.call(n_alloc, c(list(f, alloc = "optimal",
        fpc = "none", domains = "d", min_n_domain = floor_at(1500)), target))),
      "'min_n_domain' exceeds what domain(s) 'A' can hold in whole units (1000)",
      fixed = TRUE
    )
  }
  fit <- suppressWarnings(n_alloc(f, n = 1900, alloc = "optimal",
                                  fpc = "none", domains = "d",
                                  min_n_domain = floor_at(1000)))
  expect_equal(fit$detail$n_int, c(1000, 900))
})

## Stage caps from adjacent sizes at large counts

test_that("whole stage caps follow the size ratio at a billion units per parent", {
  ## Oracle: a take b fits a stage when b * N_{s-1} <= N_s. Here N / N_2 is
  ## 1.999999999, so no take of 2 fits the last stage.
  N2 <- 1e9 + 1
  f <- data.frame(stratum = "A", N = 2 * N2 - 1, N_1 = 2, N_2 = N2,
                  mean = 10, sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                  cost_2 = 1, cost_3 = 1)
  expect_gt(2 * f$N_2, f$N)
  for (fpc in c("none", "stage")) {
    expect_error(n_alloc(transform(f, take_3 = 2), n = 2e9, alloc = "optimal",
                         fpc = fpc),
                 "exceeds the 1.999999999 units per parent unit", info = fpc)
  }
  expect_error(n_alloc(f, n = 1e9, alloc = "optimal", fpc = "none",
                       resp_rate = .5),
               "holds 1.999999999 units per parent unit")
  ## One expected response at rate .6 needs 1.67 units, which fits, but
  ## 2 in whole units, which does not.
  expect_error(n_alloc(f, n = 1e6, alloc = "optimal", fpc = "none",
                       resp_rate = .6),
               "holds 1.999999999 units per parent unit, which cannot deliver one expected response in whole units")
  ## Two PSUs of 5e8 SSUs of one unit hold 1e9 at most.
  expect_error(suppressWarnings(n_alloc(f, n = 1e9 + 1, alloc = "optimal",
                                        fpc = "none")),
               "feasible totals: 1 to 1000000000")
  for (target in list(list(n = 1e6), list(budget = 1e9))) {
    fit <- suppressWarnings(do.call(n_alloc, c(list(transform(f, cost_2 = 1e6),
      alloc = "optimal", fpc = "none"), target)))
    expect_equal(fit$detail$take_3_int, 1, info = names(target))
  }
  two <- data.frame(stratum = "A", N = 2e9 - 1, N_1 = 1e9, mean = 10, sd = 2,
                    icc_1 = .1, cost_1 = 1e6, cost_2 = 1)
  fit <- suppressWarnings(n_alloc(two, n = 1e9, alloc = "optimal",
                                  fpc = "none"))
  expect_equal(fit$detail$take_2_int, 1)
  ## A take at an exact ratio still fits.
  exact <- suppressWarnings(n_alloc(transform(f, N = 2 * N2, take_3 = 2),
                                    n = 2e9, alloc = "optimal", fpc = "none"))
  expect_equal(exact$detail$take_3_int, 2)
  ## A whole take is held to the whole cap, so 1e6 does not fit 999999.9999995
  ## even though the two differ by less than 1e-12 of the cap.
  big <- data.frame(stratum = "A", N = 2e12 - 1, N_1 = 2, N_2 = 2e6,
                    mean = 10, sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                    cost_2 = 1, cost_3 = 1, take_3 = 1e6)
  expect_gt(1e6 * big$N_2, big$N)
  for (fpc in c("none", "unit", "stage")) {
    expect_error(
      n_alloc(big, n = 1e6, alloc = "optimal", fpc = fpc),
      "take_3 = 1000000 in stratum 'A' exceeds the 999999.9999995 units",
      info = fpc
    )
  }
  fit <- suppressWarnings(n_alloc(transform(big, N = 2e12), n = 1e6,
                                  alloc = "optimal", fpc = "none"))
  expect_equal(c(fit$detail$take_3, fit$detail$take_3_int), c(1e6, 1e6))
  ## A fractional take is compared with the ratio itself.
  small <- data.frame(stratum = "A", N = 52, N_1 = 10, N_2 = 26, mean = 10,
                      sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 1000,
                      cost_2 = 1, cost_3 = 1)
  expect_equal(n_alloc(transform(small, take_2 = 2.6), n = 40,
                       alloc = "optimal", fpc = "stage")$detail$take_2, 2.6)
  expect_error(n_alloc(transform(small, take_2 = 2.6000001), n = 40,
                       alloc = "optimal", fpc = "stage"),
               "take_2 = 2.6000001 in stratum 'A' exceeds the 2.6 units")
})
