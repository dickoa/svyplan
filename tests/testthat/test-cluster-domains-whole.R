## Domain floors and whole designs in several-indicator cluster allocation

.dom_tab <- function(N = 1e6, levels = c("A", "B")) {
  data.frame(name = "y", p = 0.3, cv = 0.1, icc_1 = 0.05,
             region = levels, N = N, N_1 = 100)
}

.dom_three <- function() {
  tab <- data.frame(
    name = c("y", "z"), p = c(0.3, 0.2), cv = c(0.1, 0.15),
    icc_1 = 0.05, icc_2 = 0.1, region = rep(c("A", "B"), each = 2)
  )
  tab$icc_1[3:4] <- 0.1
  tab
}

.whole_cost <- function(w, stage_cost) {
  if (length(stage_cost) == 2L) {
    w$stage_1 * (stage_cost[1] + stage_cost[2] * w$stage_2)
  } else {
    w$stage_1 * (stage_cost[1] + stage_cost[2] * w$stage_2 +
                 stage_cost[3] * w$stage_2 * w$stage_3)
  }
}

## The floor in a joint split

test_that("a joint floor holds when a PSU cost makes the total less than proportional", {
  expect_no_warning(fit <- n_cluster(
    .dom_tab(), stage_cost = c(100, 10), domains = "region",
    allocation = "joint", budget = 2450, stage_n = c(2, NA), min_n_domain = 100,
    fpc = "stage"
  ))
  expect_true(all(fit$domains$.total_n >= 100 - 1e-8))
  expect_true(all(fit$operational$domains$.total_n >= 100))
  expect_equal(sum(fit$domains$.cost), 2450, tolerance = 1e-8)
  expect_lte(fit$operational$cost, 2450)
})

test_that("a joint floor is found under population caps and across three domains", {
  capped <- n_cluster(
    .dom_tab(N = 10000), stage_cost = c(100, 10),
    domains = "region", allocation = "joint", budget = 3000, stage_n = c(2, NA),
    min_n_domain = 100, fpc = "stage"
  )
  expect_true(all(capped$domains$.total_n >= 100 - 1e-8))
  three <- n_cluster(
    .dom_tab(N = 10000, levels = c("A", "B", "C")),
    stage_cost = c(100, 10), domains = "region", allocation = "joint",
    budget = 3700, stage_n = c(2, NA), min_n_domain = 100, fpc = "stage"
  )
  expect_true(all(three$domains$.total_n >= 100 - 1e-8))
  expect_true(all(three$operational$domains$.total_n >= 100))
  expect_lte(three$operational$cost, 3700)

  expect_error(
    n_cluster(.dom_tab(N = 10000), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 3000,
              stage_n = c(2, NA), min_n_domain = 250, fpc = "stage"),
    "not achievable for domain 'A'"
  )
  expect_error(
    n_cluster(.dom_tab(), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 2390,
              stage_n = c(2, NA), min_n_domain = 100, fpc = "stage"),
    "not achievable for all domains"
  )
})

test_that("a floor only a continuous design reaches is refused at its whole cost", {
  ## Two PSUs need a take of 51 for 101 units, 1220 per domain, where the
  ## continuous take of 50.5 costs 1210.
  expect_error(
    n_cluster(.dom_tab(), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 2420,
              stage_n = c(2, NA), min_n_domain = 101, fpc = "stage"),
    "need 2440 together"
  )

  ## A take of 10.5 per SSU reaches 42, and the whole take of 10 gives 40.
  three <- data.frame(name = "y", p = 0.3, cv = 0.1, icc_1 = 0.05,
                      icc_2 = 0.1, region = c("A", "B"))
  expect_error(
    n_cluster(three, stage_cost = c(200, 40, 5),
              domains = "region", allocation = "joint", budget = 1540,
              stage_n = c(2, 2, NA), min_n_domain = 42),
    "need 1560 together"
  )
  room <- n_cluster(three, stage_cost = c(200, 40, 5),
                    domains = "region", allocation = "joint", budget = 1560,
                    stage_n = c(2, 2, NA), min_n_domain = 42)
  expect_equal(room$operational$domains$stage_3, c(11L, 11L))
})

## Floors reached by raising the takes

.floor_tab <- function(fpc, cv = 0.1, three = FALSE) {
  tab <- data.frame(name = "y", p = 0.3, cv = cv, icc_1 = if (three) 0.3 else 0.5,
                    region = c("A", "B"))
  if (three) tab$icc_2 <- 0.3
  if (fpc == "stage") {
    tab$N <- 1e6
    tab$N_1 <- if (three) 200 else 100
    if (three) tab$N_2 <- 20000
  }
  tab
}

.floor_fit <- function(tab, fpc, budget = 4800, stage_cost = c(500, 10), ...) {
  n_cluster(tab, stage_cost = stage_cost, domains = "region",
            allocation = "joint", budget = budget, fpc = fpc, ...)
}

.worst_ratio <- function(fit, cv) max(fit$domains$.cv / cv)

## A fixed take whose whole design misses the floor is refused, and it
## scores no better than any design.
.grid_ratio <- function(cv, ...) {
  tryCatch(.worst_ratio(suppressWarnings(.floor_fit(...)), cv),
           error = function(e) Inf)
}

test_that("a joint floor is reached by a larger take than the precision optimum", {
  for (fpc in c("stage", "none")) {
    tab <- .floor_tab(fpc)
    for (fl in c(30, 40, 80)) {
      expect_no_warning(fit <- .floor_fit(tab, fpc, min_n_domain = fl))
      expect_true(all(fit$domains$.total_n >= fl - 1e-8))
      expect_true(all(fit$operational$domains$.total_n >= fl))
      expect_lte(fit$operational$cost, 4800)
      ## The least take reaching the floor at each domain's share.
      expect_equal(fit$domains$stage_2,
                   fl * 500 / (fit$domains$.cost - fl * 10), tolerance = 1e-8)
      fixed <- .floor_fit(tab, fpc, min_n_domain = fl, stage_n = c(NA, fl / 2))
      expect_lte(max(fit$domains$.cv), max(fixed$domains$.cv) + 1e-8)
    }
  }
})

test_that("free takes do no worse than any fixed take under a joint floor", {
  for (fpc in c("stage", "none")) {
    tab <- .floor_tab(fpc, cv = c(0.1, 0.05))
    fit <- .floor_fit(tab, fpc, budget = 9600, min_n_domain = 30)
    expect_equal(min(fit$domains$.total_n), 30, tolerance = 1e-8)
    grid <- vapply(4:30, function(b) {
      .grid_ratio(c(0.1, 0.05), tab, fpc, budget = 9600, min_n_domain = 30,
                  stage_n = c(NA, b))
    }, numeric(1L))
    expect_lte(.worst_ratio(fit, c(0.1, 0.05)), min(grid) + 1e-6)
  }

  for (fpc in c("stage", "none")) {
    tab <- .floor_tab(fpc, cv = c(0.1, 0.08), three = TRUE)
    grid <- expand.grid(b = c(2, 3, 4, 6), cc = c(4, 8, 12, 16, 24))
    for (fl in c(80, 150)) {
      fit <- .floor_fit(tab, fpc, budget = 8000, stage_cost = c(400, 40, 4),
                        min_n_domain = fl)
      expect_true(all(fit$domains$.total_n >= fl - 1e-8))
      expect_true(all(fit$operational$domains$.total_n >= fl))
      ratios <- mapply(function(b, cc) {
        .grid_ratio(c(0.1, 0.08), tab, fpc, budget = 8000,
                    stage_cost = c(400, 40, 4), min_n_domain = fl,
                    stage_n = c(NA, b, cc))
      }, grid$b, grid$cc)
      expect_lte(.worst_ratio(fit, c(0.1, 0.08)), min(ratios) + 1e-6)
    }
  }
})

test_that("a joint floor with a fixed PSU count raises the free takes", {
  tab <- .floor_tab("none", cv = c(0.1, 0.08), three = TRUE)
  fit <- .floor_fit(tab, "none", budget = 8000, stage_cost = c(400, 40, 4),
                    stage_n = c(4, NA, NA), min_n_domain = 150)
  expect_equal(min(fit$domains$.total_n), 150, tolerance = 1e-8)
  grid <- vapply(1:8, function(b) {
    .grid_ratio(c(0.1, 0.08), tab, "none", budget = 8000,
                stage_cost = c(400, 40, 4), stage_n = c(4, b, NA),
                min_n_domain = 150)
  }, numeric(1L))
  expect_lte(.worst_ratio(fit, c(0.1, 0.08)), min(grid) + 1e-6)
})

test_that("a joint floor that does not bind leaves the design as it was", {
  for (fpc in c("stage", "none")) {
    tab <- .floor_tab(fpc)
    free <- .floor_fit(tab, fpc)
    floored <- .floor_fit(tab, fpc, min_n_domain = 25)
    expect_equal(floored$domains$stage_2, free$domains$stage_2,
                 tolerance = 1e-6)
    expect_equal(floored$operational$domains, free$operational$domains)
  }
})

test_that("a joint floor out of reach names the population or the cost", {
  capped <- .floor_tab("stage")
  capped$N <- 2000
  capped$N_1 <- 20
  expect_error(
    .floor_fit(capped, "stage", min_n_domain = 2500),
    "not achievable for domain 'A' within its population"
  )
  expect_error(
    .floor_fit(.floor_tab("none"), "none", min_n_domain = 500),
    "not achievable for domain 'A' within budget: its cheapest design costs 5500"
  )
})

## A budget that pays exactly for the least designs

test_that("an exactly affordable joint design returns its unique split", {
  fit <- n_cluster(
    .dom_tab(N = 10000), stage_cost = c(100, 10),
    domains = "region", allocation = "joint", budget = 440, stage_n = c(2, NA),
    fpc = "stage"
  )
  expect_equal(fit$domains$stage_2, c(1, 1), tolerance = 1e-8)
  expect_equal(fit$operational$domains$stage_1, c(2L, 2L))
  expect_equal(fit$operational$domains$stage_2, c(1L, 1L))
  expect_equal(fit$operational$cost, 440)

  plain <- n_cluster(
    .dom_tab()[, c("name", "p", "cv", "icc_1", "region")],
    stage_cost = c(100, 10), domains = "region", allocation = "joint",
    budget = 440, stage_n = c(2, NA)
  )
  expect_equal(plain$operational$domains$stage_2, c(1L, 1L))

  three <- n_cluster(
    .dom_tab(N = 10000, levels = c("A", "B", "C")),
    stage_cost = c(100, 10), domains = "region", allocation = "joint",
    budget = 660, stage_n = c(2, NA), fpc = "stage"
  )
  expect_equal(three$operational$cost, 660)
})

## Whole designs of domain fits

test_that("a joint fit carries a whole design within the budget", {
  tab <- .dom_tab()[, c("name", "p", "cv", "icc_1", "region")]
  fit <- n_cluster(tab, stage_cost = c(100, 10),
                   domains = "region", allocation = "joint", budget = 450,
                   stage_n = c(2, NA))
  w <- fit$operational$domains
  expect_type(w$stage_1, "integer")
  expect_equal(w$.cost, .whole_cost(w, c(100, 10)))
  expect_lte(fit$operational$cost, 450)
  expect_equal(fit$operational$cost, sum(w$.cost))
  expect_equal(w$.cv, sqrt((0.7 / 0.3) * (1 + 0.05 * (w$stage_2 - 1)) /
                             (w$stage_1 * w$stage_2)), tolerance = 1e-10)
  expect_equal(fit$total_n, sum(fit$domains$.total_n))

  ## Each domain is rounded as the scalar path rounds its share.
  for (d in c("A", "B")) {
    share <- fit$domains$.cost[fit$domains$region == d]
    one <- n_cluster(tab[tab$region == d, -5],
                     stage_cost = c(100, 10), budget = share, stage_n = c(2, NA))
    row <- w[w$region == d, ]
    expect_equal(c(row$stage_1, row$stage_2), unname(one$operational$n))
    expect_equal(row$.cv, one$operational$cv, tolerance = 1e-12)
  }

  fc <- n_cluster(.dom_tab(), stage_cost = c(100, 10),
                  domains = "region", allocation = "joint", budget = 1500,
                  fixed_cost = 50, fpc = "stage")
  expect_equal(fc$operational$cost, 50 + sum(fc$operational$domains$.cost))
  expect_lte(fc$operational$cost, 1500)
})

test_that("a separate fit keeps each domain's scalar whole design", {
  tab <- .dom_tab()
  fit <- n_cluster(tab, stage_cost = c(100, 10),
                   domains = "region", budget = 450, stage_n = c(2, NA),
                   fpc = "stage")
  one <- n_cluster(tab[1, names(tab) != "region"],
                   stage_cost = c(100, 10), budget = 450, stage_n = c(2, NA),
                   fpc = "stage")
  w <- fit$operational$domains
  expect_equal(c(w$stage_1[1], w$stage_2[1]), unname(one$operational$n))
  expect_equal(w$.cost[1], one$operational$cost)
  expect_equal(w$.cv[1], one$operational$cv, tolerance = 1e-12)
  expect_true(all(w$.cost <= 450))

  three <- .dom_three()
  sep <- n_cluster(three, stage_cost = c(200, 40, 5),
                   domains = "region", budget = 3000, fixed_cost = 100)
  for (d in c("A", "B")) {
    alone <- n_cluster(three[three$region == d, -6],
                       stage_cost = c(200, 40, 5), budget = 3000,
                       fixed_cost = 100)
    row <- sep$operational$domains[sep$operational$domains$region == d, ]
    expect_equal(c(row$stage_1, row$stage_2, row$stage_3),
                 unname(alone$operational$n))
    expect_equal(row$.cost, alone$operational$cost)
    ## z has the larger cv and y the larger ratio to its target.
    ratio <- alone$operational$cv_by_target / alone$detail$.cv_target
    expect_identical(alone$detail$name[which.max(ratio)], "y")
    expect_gt(alone$operational$cv_by_target[2], alone$operational$cv_by_target[1])
    expect_identical(row$.binding, "y")
  }
})

test_that("a three-stage joint fit stays within the budget", {
  fit <- n_cluster(.dom_three(), stage_cost = c(200, 40, 5),
                   domains = "region", allocation = "joint", budget = 6000)
  w <- fit$operational$domains
  expect_equal(w$.cost, .whole_cost(w, c(200, 40, 5)))
  expect_lte(sum(w$.cost), 6000)
  expect_equal(w$.total_n, w$stage_1 * w$stage_2 * w$stage_3)
})

test_that("the print shows the whole domain designs and their cost", {
  fit <- n_cluster(.dom_tab(), stage_cost = c(100, 10),
                   domains = "region", allocation = "joint", budget = 450,
                   stage_n = c(2, NA), fpc = "stage")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Total n = 4, cost = ", out, fixed = TRUE)))
  expect_false(any(grepl("unrounded", out, fixed = TRUE)))
  rows <- grep("^ [AB] ", out, value = TRUE)
  expect_length(rows, 2L)
  expect_true(all(grepl("^ [AB] +2 +1 +2 +1\\.0795 +220 ", rows)))
  expect_match(format(fit), "total_n=4 (unrounded: 5)", fixed = TRUE)

  sep <- n_cluster(.dom_tab(), stage_cost = c(100, 10),
                   domains = "region", budget = 450, stage_n = c(2, NA),
                   fpc = "stage")
  expect_true(any(grepl("budget = 450 for each domain, total cost = 880",
                        capture.output(print(sep)), fixed = TRUE)))
})

test_that("joint domain shares reserve the cost of a whole design at the floor", {
  # a floor of 41 on PSUs of 20 units needs 3 PSUs whole (cost 342 at take
  # 14) against 2.05 continuous (246), and a share between the two used to
  # be rejected although the budget paid for both domains
  tab <- data.frame(name = "y", region = c("A", "B"), p = c(0.2, 0.5),
                    cv = c(0.15, 0.2), icc_1 = 0.4, N = 400, N_1 = 20)
  fit_at <- function(budget) {
    n_cluster(tab, stage_cost = c(100, 1), domains = "region",
              allocation = "joint", budget = budget, min_n_domain = 41,
              fpc = "stage")
  }
  fit <- fit_at(690)
  d <- fit$operational$domains
  expect_true(all(d$.total_n >= 41))
  expect_lte(fit$operational$cost, 690)
  expect_equal(d$stage_1, c(3, 3))

  # within 1 percent of the best whole design by enumeration
  cvw <- function(a, b, p) {
    prec_cluster(n = c(a, b), icc = 0.4, unit_relvar = (1 - p) / p,
                 N = c(20, 400), fpc = "stage")$cv
  }
  g <- expand.grid(a = 1:20, b = 1:20)
  g <- g[g$a * g$b >= 41, ]
  g$cost <- g$a * (100 + g$b)
  g$rA <- mapply(cvw, g$a, g$b, 0.2) / 0.15
  g$rB <- mapply(cvw, g$a, g$b, 0.5) / 0.2
  for (budget in c(690, 1000, 2000)) {
    best <- min(vapply(seq_len(nrow(g)), function(i) {
      ok <- g$cost <= budget - g$cost[i]
      if (any(ok)) max(g$rA[i], min(g$rB[ok])) else Inf
    }, numeric(1)))
    d <- fit_at(budget)$operational$domains
    expect_lte(max(d$.cv / c(0.15, 0.2)), best * 1.01)
  }

  expect_error(fit_at(683), "need 684 together")
  expect_equal(
    svyplan:::.stage_min_total_cost(41, c(100, 1), c(1, 1), c(20, 20), whole = TRUE),
    342
  )
})

test_that("a fixed PSU count above what the floor needs still prices the floor", {
  tab <- data.frame(name = "y", region = c("A", "B"), p = c(0.2, 0.5),
                    cv = c(0.15, 0.2), icc_1 = 0.1, N = 400, N_1 = 20)
  fit_at <- function(...) {
    n_cluster(tab, stage_cost = c(100, 1), domains = "region",
              allocation = "joint", budget = 1200, stage_n = c(5, NA), fpc = "stage",
              ...)
  }
  free <- fit_at()
  floored <- fit_at(min_n_domain = 1)
  expect_equal(floored$operational$domains$.total_n,
               free$operational$domains$.total_n)
  expect_true(all(fit_at(min_n_domain = 50)$operational$domains$.total_n >= 50))

  least <- function(m, co, lo, hi) {
    svyplan:::.stage_min_total_cost(m, co, lo, hi, whole = TRUE)
  }
  # five PSUs fixed, a take of one: 5 * (100 + 1)
  expect_equal(least(1, c(100, 1), c(5, 1), c(5, 20)), 505)
  # an SSU take below its box is never priced
  expect_equal(least(1, c(100, 1, 5), c(1, 1, 2), c(20, 20, 10)), 111)
  set.seed(9)
  for (i in 1:60) {
    m <- sample(1:120, 1)
    co <- c(runif(1, 10, 200), runif(1, 1, 10), sample(c(0, 2), 1))
    alo <- sample(1:6, 1)
    ahi <- alo + sample(0:15, 1)
    blo <- sample(1:4, 1)
    bhi <- blo + sample(0:20, 1)
    clo <- if (co[3] > 0) sample(1:3, 1) else 1
    chi <- if (co[3] > 0) clo + sample(0:6, 1) else 1
    g <- expand.grid(a = alo:ahi, b = blo:bhi, c = clo:chi)
    g <- g[g$a * g$b * g$c >= m, ]
    ref <- if (nrow(g)) min(g$a * (co[1] + co[2] * g$b + co[3] * g$b * g$c)) else Inf
    if (co[3] > 0) {
      expect_equal(least(m, co, c(alo, blo, clo), c(ahi, bhi, chi)), ref)
    } else {
      expect_equal(least(m, co[1:2], c(alo, blo), c(ahi, bhi)), ref)
    }
  }
})

## The joint split

.dom_ratio <- function(fit, tab, stage_cost, ...) {
  vapply(seq_len(nrow(fit$domains)), function(d) {
    sub <- tab[tab$region == fit$domains$region[d], names(tab) != "region"]
    one <- suppressWarnings(n_cluster(sub, stage_cost = stage_cost,
                                      budget = fit$domains$.cost[d], ...))
    max(one$detail$.cv_achieved / one$detail$.cv_target)
  }, numeric(1L))
}

.dom_four <- function() {
  data.frame(
    name = paste0(rep(c("y", "z"), 4), rep(1:4, each = 2)),
    p = c(0.3, 0.1, 0.25, 0.4, 0.5, 0.2, 0.15, 0.35),
    cv = c(0.1, 0.15, 0.12, 0.1, 0.08, 0.2, 0.15, 0.1),
    icc_1 = c(0.05, 0.02, 0.1, 0.08, 0.03, 0.06, 0.12, 0.04),
    region = rep(c("A", "B", "C", "D"), each = 2)
  )
}

test_that("a joint split without correction is proportional to each domain's cost", {
  tab <- .dom_four()[1:6, ]
  cost <- c(400, 40)
  fit <- n_cluster(tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 1e5)
  ## The ratio falls as 1 / sqrt(budget), so b * ratio^2 is each domain's
  ## cost at its targets, read here from a fit at the whole budget.
  k <- vapply(c("A", "B", "C"), function(lev) {
    one <- n_cluster(tab[tab$region == lev, -5], stage_cost = cost,
                     budget = 1e5)
    1e5 * max(one$detail$.cv_achieved / one$detail$.cv_target)^2
  }, numeric(1L))
  expect_equal(fit$domains$.cost, unname(1e5 * k / sum(k)), tolerance = 1e-8)
  r <- .dom_ratio(fit, tab, cost)
  expect_equal(r, rep(max(r), 3), tolerance = 1e-8)
})

test_that("two domains split as a search over the share finds", {
  tab <- .dom_four()[1:4, ]
  cost <- c(400, 40)
  stage <- transform(tab, N = rep(c(6000, 30000), each = 2),
                     N_1 = rep(c(50, 250), each = 2))
  cases <- list(
    list(tab = tab, args = list()),
    list(tab = tab, args = list(stage_n = c(20, NA))),
    list(tab = stage, args = list(fpc = "stage"))
  )
  for (cs in cases) {
    fit <- do.call(n_cluster, c(list(x = cs$tab, stage_cost = cost,
                                     domains = "region",
                                     allocation = "joint", budget = 6e4),
                                cs$args))
    worst <- function(w) {
      max(vapply(1:2, function(d) {
        sub <- cs$tab[cs$tab$region == c("A", "B")[d], names(cs$tab) != "region"]
        one <- suppressWarnings(do.call(n_cluster, c(
          list(x = sub, stage_cost = cost,
               budget = c(w, 1 - w)[d] * 6e4), cs$args)))
        max(one$detail$.cv_achieved / one$detail$.cv_target)
      }, numeric(1L)))
    }
    best <- optimize(worst, c(0.05, 0.95), tol = 1e-10)
    expect_equal(worst(fit$domains$.cost[1] / 6e4), best$objective,
                 tolerance = 1e-7)
  }
})

test_that("a joint split under a population correction equalizes the ratios", {
  tab <- transform(.dom_four(), N_1 = rep(c(60, 150, 400, 1000), each = 2))
  tab$N <- tab$N_1 * 120
  cost <- c(400, 40)
  fit <- n_cluster(tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 1.5e5, fpc = "stage")
  r <- .dom_ratio(fit, tab, cost, fpc = "stage")
  expect_equal(r, rep(max(r), 4), tolerance = 1e-7)
  expect_lte(sum(fit$domains$.cost), 1.5e5 * (1 + 1e-12))
})

test_that("a joint split keeps a whole design for a domain that needs almost nothing", {
  tab <- data.frame(name = paste0("y", 1:3), p = c(0.3, 0.3, 0.99999),
                    cv = 0.1, icc_1 = 0.05, region = c("A", "B", "C"))
  expect_warning(
    fit <- n_cluster(tab, stage_cost = c(500, 50),
                     domains = "region", allocation = "joint", budget = 1e5),
    "single PSU"
  )
  expect_gte(fit$domains$.cost[3], 550)
  w <- fit$operational$domains
  expect_gte(w$stage_1[3], 1)
  expect_lte(sum(w$.cost), 1e5)
  r <- .dom_ratio(fit, tab, c(500, 50))
  expect_equal(r[1], r[2], tolerance = 1e-8)
  expect_lt(r[3], r[1])
})

test_that("a joint floor is held and the other domains share the worst ratio", {
  tab <- .dom_four()[1:6, ]
  cost <- c(400, 40)
  fit <- n_cluster(tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 1e5, min_n_domain = 400)
  expect_true(all(fit$domains$.total_n >= 400 * (1 - 1e-8)))
  above <- fit$domains$.total_n > 400 * (1 + 1e-6)
  expect_identical(sum(above), 2L)
  expect_true(any(!above))
  r <- .dom_ratio(fit, tab, cost)
  expect_equal(r[above], rep(max(r), sum(above)), tolerance = 1e-7)
  expect_true(all(r[!above] <= max(r) * (1 + 1e-9)))
})

test_that("a joint split follows its domains and not their order or the currency", {
  tab <- .dom_four()
  cost <- c(400, 40, 4)
  tab$icc_2 <- 0.1
  fit <- n_cluster(tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 2e5)
  rev_tab <- tab[nrow(tab):1, ]
  back <- n_cluster(rev_tab, stage_cost = cost, domains = "region",
                    allocation = "joint", budget = 2e5)
  expect_equal(back$domains$.cost[match(fit$domains$region, back$domains$region)],
               fit$domains$.cost, tolerance = 1e-8)
  ten <- n_cluster(tab, stage_cost = 10 * cost, domains = "region",
                   allocation = "joint", budget = 2e6)
  expect_equal(ten$domains$.cost / 10, fit$domains$.cost, tolerance = 1e-8)
})

test_that("the split reaches the least worst ratio on flat and jumping curves", {
  split <- svyplan:::.joint_split
  curves <- list(
    flat = list(function(b) max(1, 2 / sqrt(b)), function(b) 1 / sqrt(b)),
    step = list(function(b) if (b < 4) 3 else 1, function(b) 2 / sqrt(b)),
    floor = list(function(b) if (b < 3) Inf else 1 / sqrt(b),
                 function(b) 3 / sqrt(b))
  )
  best <- c(flat = 1, step = 1, floor = NA)
  for (nm in names(curves)) {
    f <- curves[[nm]]
    b <- split(function(d, x) f[[d]](x), lower = c(1, 1), total = 10,
               labels = c("a", "b"))
    expect_lte(sum(b), 10 * (1 + 1e-12))
    expect_true(all(b >= 1))
    worst <- max(f[[1]](b[1]), f[[2]](b[2]))
    if (is.na(best[[nm]])) {
      ## Equal ratios at 1 / sqrt(b1) = 3 / sqrt(b2) give b1 = 1.
      ## The floor at 3 holds domain a there instead.
      expect_equal(worst, 3 / sqrt(7), tolerance = 1e-8)
    } else {
      expect_lte(worst, best[[nm]] * (1 + 1e-8))
    }
  }
})

test_that("a split taken as scaling is kept only when its ratios agree", {
  ## The first ratio falls faster than 1 / sqrt(b), so the proportional
  ## split holds it at its lower budget with a ratio above the other's.
  f <- list(function(b) 100 / b^3, function(b) 1 / sqrt(b))
  b <- svyplan:::.joint_split(function(d, x) f[[d]](x), lower = c(4, 1),
                              total = 10, labels = c("a", "b"),
                              scaling = TRUE)
  ## Both ratios are falling, so the optimum is where they cross.
  b1 <- uniroot(function(x) f[[1]](x) - f[[2]](10 - x), c(4, 9),
                tol = 1e-14)$root
  expect_equal(b, c(b1, 10 - b1), tolerance = 1e-8)
  expect_equal(max(f[[1]](b[1]), f[[2]](b[2])), f[[1]](b1), tolerance = 1e-9)
})

## Fractional fixed sizes in a joint split

.dom_frac <- function() {
  data.frame(name = rep(c("a", "b"), 3), p = c(0.3, 0.1, 0.4, 0.2, 0.5, 0.25),
             cv = c(0.1, 0.15, 0.2, 0.2, 0.3, 0.3), icc_1 = 0.05,
             region = rep(c("A", "B", "C"), each = 2))
}

test_that("a fixed PSU count that rounds down floors each share at its continuous cost", {
  ## 2.4 PSUs at a take of 1 cost 264, and the 2 PSUs fielded cost 220.
  tab <- .dom_frac()
  fit <- n_cluster(tab, stage_cost = c(100, 10), budget = 800,
                   stage_n = c(2.4, NA), domains = "region", allocation = "joint")
  dom <- fit$domains
  expect_equal(sum(dom$.cost), 800, tolerance = 1e-10)
  expect_true(all(dom$.cost >= 264 * (1 - 1e-12)))
  for (lev in dom$region) {
    at <- dom$region == lev
    one <- n_cluster(tab[tab$region == lev, 1:4],
                     stage_cost = c(100, 10), budget = dom$.cost[at],
                     stage_n = c(2.4, NA))
    expect_equal(dom$stage_2[at], one$n[["stage_2"]], tolerance = 1e-8)
    expect_equal(dom$.cv[at], max(one$detail$.cv_achieved), tolerance = 1e-8)
  }
  worst <- dom$.cv / vapply(dom$region, function(l)
    tab$cv[tab$region == l & tab$name == "b"], numeric(1L))
  above <- dom$.cost > 264 * (1 + 1e-9)
  expect_true(any(above))
  expect_true(all(worst[!above] <= max(worst[above]) * (1 + 1e-8)))
  expect_lte(fit$operational$cost, 800)
})

test_that("a joint budget at exactly the summed floors fits and one just below is refused", {
  fit <- n_cluster(.dom_frac(), stage_cost = c(100, 10),
                   budget = 792, stage_n = c(2.4, NA), domains = "region",
                   allocation = "joint")
  expect_equal(fit$domains$.cost, rep(264, 3), tolerance = 1e-10)
  expect_equal(fit$domains$stage_2, rep(1, 3), tolerance = 1e-10)
  expect_error(
    n_cluster(.dom_frac(), stage_cost = c(100, 10),
              budget = 792 * (1 - 1e-6), stage_n = c(2.4, NA), domains = "region",
              allocation = "joint"),
    "costs 792 plus fixed_cost = 0"
  )
})

test_that("a fixed PSU count that rounds up floors each share at its whole cost", {
  ## 2.6 PSUs at a take of 1 cost 286, and the 3 PSUs fielded cost 330.
  tab <- .dom_tab(levels = c("A", "B"))[, 1:5]
  expect_error(
    n_cluster(tab, stage_cost = c(100, 10), budget = 659,
              stage_n = c(2.6, NA), domains = "region", allocation = "joint"),
    "costs 660 plus"
  )
  fit <- n_cluster(tab, stage_cost = c(100, 10), budget = 660,
                   stage_n = c(2.6, NA), domains = "region", allocation = "joint")
  expect_identical(fit$operational$domains$stage_1, c(3L, 3L))
  expect_lte(fit$operational$cost, 660)
})

test_that("fractional fixed takes under the stage correction fit at their continuous floor", {
  ## One PSU at a take of 2.4 costs 124 and fields a take of 2 at 120.
  tab <- .dom_tab(levels = c("A", "B"))
  args <- list(x = tab, stage_cost = c(100, 10), stage_n = c(NA, 2.4),
               domains = "region", allocation = "joint", fpc = "stage")
  fits <- function(...) suppressWarnings(n_cluster(...))
  fit <- do.call(fits, c(args, budget = 248))
  expect_equal(fit$domains$.cost, c(124, 124), tolerance = 1e-10)
  expect_lte(fit$operational$cost, 248)
  expect_no_error(do.call(fits, c(args, budget = 260)))
  expect_error(do.call(fits, c(args, budget = 247.99)), "costs 248 plus")

  ## Three stages at takes 2.4 and 3.4 cost 132.16 per PSU, 126 whole. At
  ## the exact share the SSU take bound is met to rounding only.
  tab3 <- cbind(tab, icc_2 = 0.1, N_2 = 1000)
  args3 <- list(x = tab3, stage_cost = c(100, 10, 1),
                stage_n = c(NA, 2.4, 3.4), domains = "region",
                allocation = "joint", fpc = "stage")
  fit3 <- do.call(fits, c(args3, budget = 264.32))
  expect_equal(fit3$domains$.cost, c(132.16, 132.16), tolerance = 1e-10)
  expect_lte(fit3$operational$cost, 264.32)
  expect_no_error(do.call(fits, c(args3, budget = 270)))
  expect_error(do.call(fits, c(args3, budget = 264.3)),
               "costs 264.3 plus")
})

test_that("a three-stage fixed PSU count that rounds down floors at its continuous cost", {
  ## 2.4 PSUs at takes of 1 cost 266.4, and the 2 PSUs fielded cost 222.
  tab3 <- cbind(.dom_tab(levels = c("A", "B"))[, 1:5], icc_2 = 0.1)
  args <- list(x = tab3, stage_cost = c(100, 10, 1), stage_n = c(2.4, NA, NA),
               domains = "region", allocation = "joint")
  fit <- do.call(n_cluster, c(args, budget = 532.8))
  expect_equal(fit$domains$.cost, c(266.4, 266.4), tolerance = 1e-10)
  expect_error(do.call(n_cluster, c(args, budget = 532.7)), "costs 532.8 plus")
})

test_that("a domain floor is reached under a fractional fixed PSU count", {
  tab <- .dom_tab(levels = c("A", "B"))[, 1:5]
  for (np in c(2.4, 2.6)) {
    fit <- n_cluster(tab, stage_cost = c(100, 10), budget = 3000,
                     stage_n = c(np, NA), min_n_domain = 100, domains = "region",
                     allocation = "joint")
    expect_true(all(fit$domains$.total_n >= 100 - 1e-8))
    expect_true(all(fit$operational$domains$.total_n >= 100))
    expect_lte(fit$operational$cost, 3000)
  }
  ## 2.4 PSUs reach 100 at 240 + 10 * 100 = 1240, and the 2 PSUs fielded
  ## at a take of 50 cost 1200.
  args <- list(x = tab, stage_cost = c(100, 10), stage_n = c(2.4, NA),
               min_n_domain = 100, domains = "region", allocation = "joint")
  fit <- do.call(n_cluster, c(args, budget = 2480))
  expect_equal(fit$domains$.cost, c(1240, 1240), tolerance = 1e-10)
  expect_error(do.call(n_cluster, c(args, budget = 2479)),
               "need 2480 together")
  ## A fixed take of 2.4 fields 2, and 10 PSUs reach 20.
  fit <- n_cluster(.dom_tab(levels = c("A", "B")),
                   stage_cost = c(100, 10), budget = 5000, stage_n = c(NA, 2.4),
                   min_n_domain = 20, domains = "region",
                   allocation = "joint", fpc = "stage")
  expect_true(all(fit$operational$domains$.total_n >= 20))
})

test_that("free PSUs keep at least one PSU at the continuous floor", {
  ## One PSU at a continuous take of 2.4 costs 124 (fielded take of 2 costs 120).
  tab <- .dom_tab(levels = c("A", "B"))[, 1:5]
  args <- list(x = tab, stage_cost = c(100, 10), stage_n = c(NA, 2.4),
               domains = "region", allocation = "joint")
  fits <- function(...) suppressWarnings(n_cluster(...))
  fit <- do.call(fits, c(args, budget = 248))
  expect_equal(fit$domains$.cost, c(124, 124), tolerance = 1e-10)
  expect_error(do.call(fits, c(args, budget = 247.99)), "costs 248 plus")
})

## Floors against the solvers

.dom_two <- function(...) {
  data.frame(name = "y", p = 0.3, cv = 0.1, icc_1 = 0.05,
             region = c("A", "B"), ...)
}

## Each domain's whole design equals the scalar fit of that domain at its
## share, so the joint floor let through a share the solver can field.
.expect_domains_as_scalar <- function(fit, tab, ...) {
  for (lev in fit$domains$region) {
    share <- fit$domains$.cost[fit$domains$region == lev]
    one <- suppressWarnings(n_cluster(
      tab[tab$region == lev, names(tab) != "region"],
      budget = share, ...
    ))
    op <- fit$operational$domains[fit$operational$domains$region == lev, ]
    expect_equal(unname(unlist(op[names(one$operational$n)])),
                 unname(one$operational$n))
  }
}

test_that("a domain floor already met keeps the minimum-PSU allowance", {
  tab <- .dom_two()
  args <- list(x = tab, stage_cost = c(100, 10), budget = 248,
               stage_n = c(NA, 2.4), domains = "region", allocation = "joint")
  plain <- suppressWarnings(do.call(n_cluster, args))
  floored <- suppressWarnings(do.call(n_cluster, c(args, min_n_domain = 1)))
  expect_equal(floored$domains$.cost, plain$domains$.cost)
  expect_equal(floored$domains$stage_1, c(1, 1), tolerance = 1e-10)
  .expect_domains_as_scalar(floored, tab, stage_cost = c(100, 10),
                            stage_n = c(NA, 2.4))

  tab3 <- .dom_two(icc_2 = 0.1)
  args3 <- list(x = tab3, stage_cost = c(100, 10, 1), budget = 264.32,
                stage_n = c(NA, 2.4, 3.4), domains = "region",
                allocation = "joint")
  plain3 <- suppressWarnings(do.call(n_cluster, args3))
  floored3 <- suppressWarnings(do.call(n_cluster, c(args3, min_n_domain = 1)))
  expect_equal(floored3$domains$.cost, plain3$domains$.cost)

  ## A floor of 30 at a fixed take of 2.4 and a free SSU take. One PSU at
  ## continuous take 2.4 and SSU take 12.5 costs 154 and reaches 30.
  args30 <- list(x = tab3, stage_cost = c(100, 10, 1),
                 stage_n = c(NA, 2.4, NA), min_n_domain = 30, domains = "region",
                 allocation = "joint")
  fit <- suppressWarnings(do.call(n_cluster, c(args30, budget = 308)))
  expect_true(all(fit$domains$.total_n >= 30 - 1e-8))
  expect_true(all(fit$operational$domains$.total_n >= 30))
  expect_error(suppressWarnings(do.call(n_cluster, c(args30, budget = 307))),
               "need 308 together")
})

test_that("a fixed PSU count below one is refused", {
  tab <- .dom_two()
  expect_error(
    n_cluster(tab, stage_cost = c(100, 10), budget = 220,
              stage_n = c(0.4, NA), domains = "region", allocation = "joint"),
    ">= 1"
  )
})

test_that("a fixed take above the whole population cap is priced at the cap", {
  ## PSUs of 2.8 units field a take of 2 for a fixed 2.6, so each domain's
  ## floor is the continuous 126, not 130 for a take of 3.
  tab <- .dom_two(N = 28, N_1 = 10)
  args <- list(x = tab, stage_cost = c(100, 10), stage_n = c(NA, 2.6),
               domains = "region", allocation = "joint", fpc = "stage")
  fit <- suppressWarnings(do.call(n_cluster, c(args, budget = 252)))
  expect_identical(fit$operational$domains$stage_2, c(2L, 2L))
  .expect_domains_as_scalar(fit, tab, stage_cost = c(100, 10),
                            stage_n = c(NA, 2.6), fpc = "stage")
  expect_error(suppressWarnings(do.call(n_cluster, c(args, budget = 251.9))),
               "costs 252 plus")
})

test_that("free PSUs under the stage correction keep at least one PSU at a floor", {
  ## The continuous design holds one PSU at least, and at a share of 150 it
  ## reaches 26 of the 30. The least share reaching 30 is 154.
  tab <- .dom_two(icc_2 = 0.1, N = 1e6, N_1 = 100, N_2 = 2000)
  args <- list(x = tab, stage_cost = c(100, 10, 1), stage_n = c(NA, 2.4, NA),
               min_n_domain = 30, domains = "region", allocation = "joint",
               fpc = "stage")
  fit <- suppressWarnings(do.call(n_cluster, c(args, budget = 308)))
  expect_true(all(fit$domains$stage_1 >= 1 - 1e-12))
  expect_true(all(fit$domains$.total_n >= 30 - 1e-8))
  expect_error(suppressWarnings(do.call(n_cluster, c(args, budget = 307))),
               "need 308 together")
})

## The floor precheck of a joint budget

test_that("a floor search stopped at its limit does not refuse a funded budget", {
  d <- data.frame(region = c("A", "B"), name = "y", p = .3, cv = .2,
                  icc_1 = .15, icc_2 = .10, N_1 = 8, N_2 = 120, N = 1320)
  joint <- function(budget, max_nodes) {
    n_cluster(d, stage_cost = c(18, 9, 22), domains = "region",
              allocation = "joint", budget = budget, min_n_domain = 451,
              control = list(max_nodes = max_nodes))
  }
  expect_error(joint(21192, 1), class = "svyplan_plan_unknown")
  fit <- joint(21192, 1e5)
  expect_identical(unname(as.matrix(fit$operational$domains[
    , c("stage_1", "stage_2", "stage_3")])), matrix(c(3L, 3L, 14L, 14L, 11L, 11L), 2))
  expect_equal(fit$operational$cost, 21192)
  for (max_nodes in c(1, 1e5)) {
    expect_error(joint(20000, max_nodes), "not achievable for all domains",
                 info = max_nodes)
  }
})

test_that("the relaxed whole floor never exceeds the whole minimum", {
  checked <- 0L
  for (take in c(2.4, 2.6, 3, 3.5)) {
    for (m in c(6, 20, 45)) {
      for (cost in list(c(10, 1), c(3, 5))) {
        lower_int <- c(1, round(take))
        upper_int <- c(Inf, round(take))
        relaxed <- .stage_min_total_cost(m, cost, lower_int, upper_int,
                                         whole = FALSE)
        whole <- .stage_min_total_cost(m, cost, lower_int, upper_int,
                                       whole = TRUE)
        expect_false(is.na(whole))
        expect_lte(relaxed, whole * (1 + 1e-12))
        checked <- checked + 1L
      }
    }
  }
  expect_identical(checked, 24L)
  expect_equal(.stage_min_total_cost(6, c(10, 1), c(1, 2.6), c(Inf, 2.6)),
               29.0769230769, tolerance = 1e-9)
  expect_equal(.stage_min_total_cost(6, c(10, 1), c(1, 3), c(Inf, 3),
                                    whole = TRUE), 26)
  expect_true(is.na(.stage_min_total_cost(451, c(18, 9, 22), c(1, 1, 1),
                                          c(8, 15, 11), whole = TRUE,
                                          max_nodes = 1)))
})

## The whole design of a joint budget

## Oracle (rung 2): every two-stage design of each domain within the
## budget, priced from the CV formula of ?n_cluster.
joint_points <- function(p, icc, cv, c1, c2, budget, r = 1, M1 = Inf,
                         M2 = Inf, m = 0, take = NULL) {
  V <- (1 - p) / p
  out <- list()
  for (a in seq_len(min(floor(M1 + 1e-9), floor(budget / (c1 + c2))))) {
    b <- if (is.null(take)) {
      b_max <- min(floor(M2 + 1e-9), floor((budget / a - c1) / c2))
      if (b_max < ceiling(1 / r - 1e-9)) next
      seq.int(ceiling(1 / r - 1e-9), b_max)
    } else {
      take
    }
    b <- b[a * b >= m & a * (c1 + b * c2) <= budget]
    if (!length(b)) next
    cv2 <- V * icc / a * (1 - a / M1) +
      V * (1 - icc) / (a * b * r) * (1 - b * r / M2)
    out[[length(out) + 1L]] <- data.frame(cost = a * (c1 + b * c2),
                                          ratio = sqrt(cv2) / cv)
  }
  do.call(rbind, out)
}

joint_best <- function(pa, pb, budget) {
  best <- Inf
  for (i in seq_len(nrow(pa))) {
    ok <- pb$cost <= budget - pa$cost[i] + 1e-9
    if (any(ok)) best <- min(best, max(pa$ratio[i], min(pb$ratio[ok])))
  }
  best
}

test_that("a joint budget fields the least worst ratio over whole designs", {
  di <- data.frame(region = c("A", "B"), name = "y", p = c(.2, .4),
                   cv = c(.3, .2), icc_1 = c(.1, .2))
  fit <- suppressWarnings(n_cluster(di, stage_cost = c(4, 1),
                                    domains = "region", allocation = "joint",
                                    budget = 12))
  best <- joint_best(joint_points(.2, .1, .3, 4, 1, 12),
                     joint_points(.4, .2, .2, 4, 1, 12), 12)
  expect_equal(best, 4.944132, tolerance = 1e-6)
  expect_equal(max(fit$operational$domains$.cv / di$cv), best,
               tolerance = 1e-9)
  expect_equal(fit$operational$cost, 12)
  expect_identical(fit$operational$search$status, "optimal")
  expect_lte(fit$operational$search$bound, fit$operational$search$ratio)
  expect_lte(fit$operational$search$gap, 1e-9)
})

test_that("joint whole designs match enumeration over floors, costs and caps", {
  set.seed(7)
  checked <- 0L
  refused <- 0L
  for (i in 1:28) {
    stage <- i %% 3 == 0
    floor_m <- if (i %% 4 == 0) sample(4:10, 1) else 0
    fc <- if (i %% 5 == 0) 20 else 0
    r <- if (i %% 2) 1 else .8
    take <- if (i > 24) 2.6 else NULL
    p <- runif(2, .1, .5)
    icc <- runif(2, .02, .3)
    cv <- runif(2, .1, .4)
    c1 <- sample(5:30, 1)
    c2 <- sample(1:5, 1)
    budget <- sample(80:220, 1)
    M1 <- if (stage) sample(6:12, 2, TRUE) else c(Inf, Inf)
    M2 <- if (stage) sample(3:8, 2, TRUE) else c(Inf, Inf)
    tab <- data.frame(region = c("A", "B"), name = "y", p = p, cv = cv,
                      icc_1 = icc)
    args <- list(tab, stage_cost = c(c1, c2), domains = "region",
                 allocation = "joint", budget = budget, resp_rate = r,
                 fixed_cost = fc)
    if (stage) {
      args[[1]]$N_1 <- M1
      args[[1]]$N <- M1 * M2
      args$fpc <- "stage"
    }
    if (floor_m > 0) args$min_n_domain <- floor_m
    if (!is.null(take)) args$stage_n <- c(NA, take)
    fit <- tryCatch(suppressWarnings(do.call(n_cluster, args)),
                    error = function(e) e)
    whole_take <- if (is.null(take)) NULL else round(take)
    pts <- lapply(1:2, function(d) {
      joint_points(p[d], icc[d], cv[d], c1, c2, budget - fc, r, M1[d], M2[d],
                   floor_m, whole_take)
    })
    best <- if (any(vapply(pts, is.null, TRUE))) Inf else
      joint_best(pts[[1]], pts[[2]], budget - fc)
    info <- paste("fixture", i)
    if (!is.finite(best)) {
      expect_true(inherits(fit, "error"), info = info)
      refused <- refused + 1L
      next
    }
    expect_false(inherits(fit, "error"), info = info)
    if (inherits(fit, "error")) next
    expect_equal(max(fit$operational$domains$.cv / cv), best,
                 tolerance = 1e-9, info = info)
    expect_lte(fit$operational$cost, budget * (1 + 1e-9))
    expect_identical(fit$operational$search$status, "optimal", info = info)
    checked <- checked + 1L
  }
  expect_identical(c(checked, refused), c(25L, 3L))
})

test_that("a joint whole search stopped at its node limit says so", {
  di <- data.frame(region = c("A", "B"), name = "y", p = c(.2, .4),
                   cv = c(.3, .2), icc_1 = c(.1, .2))
  fit <- suppressWarnings(n_cluster(di, stage_cost = c(4, 1),
                                    domains = "region", allocation = "joint",
                                    budget = 12, control = list(max_nodes = 8)))
  expect_identical(fit$operational$search$status, "limit")
  expect_lte(fit$operational$search$bound, fit$operational$search$ratio)
  expect_lte(fit$operational$cost, 12)
})

test_that("the joint minimax improves a poor first design and says when it stops", {
  domain <- function(p, icc, cv) {
    V <- (1 - p) / p
    list(R = matrix(c(V * icc, V * (1 - icc)), 1), rates = matrix(1, 1, 2),
         M = matrix(Inf, 1, 2), cv_t = cv, lower_int = c(1, 1),
         upper_int = c(Inf, Inf))
  }
  dd <- list(domain(.2, .1, .3), domain(.4, .2, .2))
  poor <- list(c(1, 1), c(1, 1))
  mm <- .joint_whole_minimax(dd, c(4, 1), 12, 0, poor, 2e7)
  best <- joint_best(joint_points(.2, .1, .3, 4, 1, 12),
                     joint_points(.4, .2, .2, 4, 1, 12), 12)
  expect_identical(mm$status, "optimal")
  expect_equal(mm$ratio, best, tolerance = 1e-9)
  stopped <- .joint_whole_minimax(dd, c(4, 1), 12, 0, poor, 2)
  expect_identical(stopped$status, "limit")
  expect_lte(stopped$bound, stopped$ratio)
})

## The continuous and the whole floors are checked apart

test_that("a joint budget funds whole floors another split would refuse", {
  d <- data.frame(region = c("A", "B"), name = "y", p = .3, cv = .2,
                  icc_1 = .1, icc_2 = .1, N_1 = 12, N_2 = c(12, 24),
                  N = c(120, 240))
  joint <- function(budget, fixed_cost = 0) {
    n_cluster(d, stage_cost = c(10, 1, 1), stage_n = c(NA, NA, 2.6),
              domains = "region", allocation = "joint", budget = budget,
              min_n_domain = 9, fixed_cost = fixed_cost)
  }
  ## Whole takes: A holds one SSU per PSU and B two, and the take of 2.6
  ## becomes 3, so 9 units cost 14 a in A, and in B 14 a (b = 1) or
  ## 18 a (b = 2, a >= 2). The only pair within 78 is 42 + 36.
  for (fc in c(0, 500)) {
    fit <- joint(78 + fc, fc)
    op <- fit$operational$domains
    expect_identical(unname(as.matrix(op[, c("stage_1", "stage_2", "stage_3")])),
                     matrix(c(3L, 2L, 1L, 2L, 3L, 3L), 2), info = fc)
    expect_equal(fit$operational$cost, 78 + fc)
    expect_lte(fit$cost, 78 + fc)
    expect_true(all(fit$domains$.total_n >= 9 - 1e-8))
  }
  ## The continuous floors cost 76.85 and the whole ones 78.
  expect_error(joint(77.5), "within budget in whole units: they need 78")
})

test_that("the joint record describes the returned design and its work", {
  di <- data.frame(region = c("A", "B"), name = "y", p = c(.2, .4),
                   cv = c(.3, .2), icc_1 = c(.1, .2))
  fit <- suppressWarnings(n_cluster(di, stage_cost = c(4, 1),
                                    domains = "region", allocation = "joint",
                                    budget = 12))
  s <- fit$operational$search
  expect_equal(s$ratio, max(fit$operational$domains$.cv / di$cv),
               tolerance = 1e-12)
  expect_equal(sum(s$nodes_by_domain), s$nodes)
  expect_identical(unname(s$by_domain), rep(s$status, 2))
  stopped <- suppressWarnings(n_cluster(di, stage_cost = c(4, 1),
                                        domains = "region",
                                        allocation = "joint", budget = 12,
                                        control = list(max_nodes = 8)))
  expect_lte(stopped$operational$search$nodes, 8)
  expect_equal(sum(stopped$operational$search$nodes_by_domain),
               stopped$operational$search$nodes)
})

test_that("the minimax reports the ratio its design reaches", {
  domain <- function(p, icc, cv) {
    V <- (1 - p) / p
    list(R = matrix(c(V * icc, V * (1 - icc)), 1), rates = matrix(1, 1, 2),
         M = matrix(Inf, 1, 2), cv_t = cv, lower_int = c(1, 1),
         upper_int = c(Inf, Inf))
  }
  dd <- list(domain(.2, .1, .3), domain(.4, .2, .2))
  mm <- .joint_whole_minimax(dd, c(4, 1), 12, 0, list(c(1, 1), c(1, 1)), 2e7)
  reached <- max(vapply(1:2, function(d) {
    max(sqrt(.stage_cv2(mm$n[[d]], dd[[d]]$R, dd[[d]]$rates, dd[[d]]$M)) /
          dd[[d]]$cv_t)
  }, 1))
  expect_identical(mm$ratio, reached)
  expect_lte(mm$bound, mm$ratio)
  expect_equal(sum(mm$nodes_by_domain), mm$nodes)
})
