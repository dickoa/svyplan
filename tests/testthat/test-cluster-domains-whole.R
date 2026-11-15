## Domain floors and whole designs in several-indicator cluster allocation

.dom_tab <- function(N = 1e6, levels = c("A", "B")) {
  data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
             region = levels, N = N, N_psu = 100)
}

.dom_three <- function() {
  tab <- data.frame(
    name = c("y", "z"), p = c(0.3, 0.2), cv = c(0.1, 0.15),
    icc_psu = 0.05, icc_ssu = 0.1, region = rep(c("A", "B"), each = 2)
  )
  tab$icc_psu[3:4] <- 0.1
  tab
}

.whole_cost <- function(w, stage_cost) {
  if (length(stage_cost) == 2L) {
    w$n_psu * (stage_cost[1] + stage_cost[2] * w$n_per_psu)
  } else {
    w$n_psu * (stage_cost[1] + stage_cost[2] * w$n_per_psu +
                 stage_cost[3] * w$n_per_psu * w$n_per_ssu)
  }
}

## The floor in a joint split

test_that("a joint floor holds when a PSU cost makes the total less than proportional", {
  expect_no_warning(fit <- n_cluster(
    indicators = .dom_tab(), stage_cost = c(100, 10), domains = "region",
    allocation = "joint", budget = 2450, n_psu = 2, min_n_domain = 100,
    fpc = "stage"
  ))
  expect_true(all(fit$domains$.total_n >= 100 - 1e-8))
  expect_true(all(fit$operational$domains$.total_n >= 100))
  expect_equal(sum(fit$domains$.cost), 2450, tolerance = 1e-8)
  expect_lte(fit$operational$cost, 2450)
})

test_that("a joint floor is found under population caps and across three domains", {
  capped <- n_cluster(
    indicators = .dom_tab(N = 10000), stage_cost = c(100, 10),
    domains = "region", allocation = "joint", budget = 3000, n_psu = 2,
    min_n_domain = 100, fpc = "stage"
  )
  expect_true(all(capped$domains$.total_n >= 100 - 1e-8))
  three <- n_cluster(
    indicators = .dom_tab(N = 10000, levels = c("A", "B", "C")),
    stage_cost = c(100, 10), domains = "region", allocation = "joint",
    budget = 3700, n_psu = 2, min_n_domain = 100, fpc = "stage"
  )
  expect_true(all(three$domains$.total_n >= 100 - 1e-8))
  expect_true(all(three$operational$domains$.total_n >= 100))
  expect_lte(three$operational$cost, 3700)

  expect_error(
    n_cluster(indicators = .dom_tab(N = 10000), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 3000,
              n_psu = 2, min_n_domain = 250, fpc = "stage"),
    "not achievable for domain 'A'"
  )
  expect_error(
    n_cluster(indicators = .dom_tab(), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 2390,
              n_psu = 2, min_n_domain = 100, fpc = "stage"),
    "not achievable for all domains"
  )
})

test_that("a floor only a continuous design reaches is refused at its whole cost", {
  ## Two PSUs need a take of 51 for 101 units, 1220 per domain, where the
  ## continuous take of 50.5 costs 1210.
  expect_error(
    n_cluster(indicators = .dom_tab(), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 2420,
              n_psu = 2, min_n_domain = 101, fpc = "stage"),
    "together at a cost of 2440"
  )

  ## A take of 10.5 per SSU reaches 42, and the whole take of 10 gives 40.
  three <- data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
                      icc_ssu = 0.1, region = c("A", "B"))
  expect_error(
    n_cluster(indicators = three, stage_cost = c(200, 40, 5),
              domains = "region", allocation = "joint", budget = 1540,
              n_psu = 2, n_per_psu = 2, min_n_domain = 42),
    "together at a cost of 1560"
  )
  room <- n_cluster(indicators = three, stage_cost = c(200, 40, 5),
                    domains = "region", allocation = "joint", budget = 1560,
                    n_psu = 2, n_per_psu = 2, min_n_domain = 42)
  expect_equal(room$operational$domains$n_per_ssu, c(11L, 11L))
})

## Floors reached by raising the takes

.floor_tab <- function(fpc, cv = 0.1, three = FALSE) {
  tab <- data.frame(name = "y", p = 0.3, cv = cv, icc_psu = if (three) 0.3 else 0.5,
                    region = c("A", "B"))
  if (three) tab$icc_ssu <- 0.3
  if (fpc == "stage") {
    tab$N <- 1e6
    tab$N_psu <- if (three) 200 else 100
    if (three) tab$N_ssu <- 20000
  }
  tab
}

.floor_fit <- function(tab, fpc, budget = 4800, stage_cost = c(500, 10), ...) {
  n_cluster(indicators = tab, stage_cost = stage_cost, domains = "region",
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
      expect_equal(fit$domains$n_per_psu,
                   fl * 500 / (fit$domains$.cost - fl * 10), tolerance = 1e-8)
      fixed <- .floor_fit(tab, fpc, min_n_domain = fl, n_per_psu = fl / 2)
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
                  n_per_psu = b)
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
                    n_per_psu = b, n_per_ssu = cc)
      }, grid$b, grid$cc)
      expect_lte(.worst_ratio(fit, c(0.1, 0.08)), min(ratios) + 1e-6)
    }
  }
})

test_that("a joint floor with a fixed PSU count raises the free takes", {
  tab <- .floor_tab("none", cv = c(0.1, 0.08), three = TRUE)
  fit <- .floor_fit(tab, "none", budget = 8000, stage_cost = c(400, 40, 4),
                    n_psu = 4, min_n_domain = 150)
  expect_equal(min(fit$domains$.total_n), 150, tolerance = 1e-8)
  grid <- vapply(1:8, function(b) {
    .grid_ratio(c(0.1, 0.08), tab, "none", budget = 8000,
                stage_cost = c(400, 40, 4), n_psu = 4, n_per_psu = b,
                min_n_domain = 150)
  }, numeric(1L))
  expect_lte(.worst_ratio(fit, c(0.1, 0.08)), min(grid) + 1e-6)
})

test_that("a joint floor that does not bind leaves the design as it was", {
  for (fpc in c("stage", "none")) {
    tab <- .floor_tab(fpc)
    free <- .floor_fit(tab, fpc)
    floored <- .floor_fit(tab, fpc, min_n_domain = 25)
    expect_equal(floored$domains$n_per_psu, free$domains$n_per_psu,
                 tolerance = 1e-6)
    expect_equal(floored$operational$domains, free$operational$domains)
  }
})

test_that("a joint floor out of reach names the population or the cost", {
  capped <- .floor_tab("stage")
  capped$N <- 2000
  capped$N_psu <- 20
  expect_error(
    .floor_fit(capped, "stage", min_n_domain = 2500),
    "not achievable for domain 'A': its largest design within the population has total_n = 2000"
  )
  expect_error(
    .floor_fit(.floor_tab("none"), "none", min_n_domain = 500),
    "not achievable for domain 'A' within budget: its cheapest design reaching it costs 5500"
  )
})

## A budget that pays exactly for the least designs

test_that("an exactly affordable joint design returns its unique split", {
  fit <- n_cluster(
    indicators = .dom_tab(N = 10000), stage_cost = c(100, 10),
    domains = "region", allocation = "joint", budget = 440, n_psu = 2,
    fpc = "stage"
  )
  expect_equal(fit$domains$n_per_psu, c(1, 1), tolerance = 1e-8)
  expect_equal(fit$operational$domains$n_psu, c(2L, 2L))
  expect_equal(fit$operational$domains$n_per_psu, c(1L, 1L))
  expect_equal(fit$operational$cost, 440)

  plain <- n_cluster(
    indicators = .dom_tab()[, c("name", "p", "cv", "icc_psu", "region")],
    stage_cost = c(100, 10), domains = "region", allocation = "joint",
    budget = 440, n_psu = 2
  )
  expect_equal(plain$operational$domains$n_per_psu, c(1L, 1L))

  three <- n_cluster(
    indicators = .dom_tab(N = 10000, levels = c("A", "B", "C")),
    stage_cost = c(100, 10), domains = "region", allocation = "joint",
    budget = 660, n_psu = 2, fpc = "stage"
  )
  expect_equal(three$operational$cost, 660)
})

## Whole designs of domain fits

test_that("a joint fit carries a whole design within the budget", {
  tab <- .dom_tab()[, c("name", "p", "cv", "icc_psu", "region")]
  fit <- n_cluster(indicators = tab, stage_cost = c(100, 10),
                   domains = "region", allocation = "joint", budget = 450,
                   n_psu = 2)
  w <- fit$operational$domains
  expect_type(w$n_psu, "integer")
  expect_equal(w$.cost, .whole_cost(w, c(100, 10)))
  expect_lte(fit$operational$cost, 450)
  expect_equal(fit$operational$cost, sum(w$.cost))
  expect_equal(w$.cv, sqrt((0.7 / 0.3) * (1 + 0.05 * (w$n_per_psu - 1)) /
                             (w$n_psu * w$n_per_psu)), tolerance = 1e-10)
  expect_equal(fit$total_n, sum(fit$domains$.total_n))

  ## Each domain is rounded as the scalar path rounds its share.
  for (d in c("A", "B")) {
    share <- fit$domains$.cost[fit$domains$region == d]
    one <- n_cluster(indicators = tab[tab$region == d, -5],
                     stage_cost = c(100, 10), budget = share, n_psu = 2)
    row <- w[w$region == d, ]
    expect_equal(c(row$n_psu, row$n_per_psu), unname(one$operational$n))
    expect_equal(row$.cv, one$operational$cv, tolerance = 1e-12)
  }

  fc <- n_cluster(indicators = .dom_tab(), stage_cost = c(100, 10),
                  domains = "region", allocation = "joint", budget = 1500,
                  fixed_cost = 50, fpc = "stage")
  expect_equal(fc$operational$cost, 50 + sum(fc$operational$domains$.cost))
  expect_lte(fc$operational$cost, 1500)
})

test_that("a separate fit keeps each domain's scalar whole design", {
  tab <- .dom_tab()
  fit <- n_cluster(indicators = tab, stage_cost = c(100, 10),
                   domains = "region", budget = 450, n_psu = 2,
                   fpc = "stage")
  one <- n_cluster(indicators = tab[1, names(tab) != "region"],
                   stage_cost = c(100, 10), budget = 450, n_psu = 2,
                   fpc = "stage")
  w <- fit$operational$domains
  expect_equal(c(w$n_psu[1], w$n_per_psu[1]), unname(one$operational$n))
  expect_equal(w$.cost[1], one$operational$cost)
  expect_equal(w$.cv[1], one$operational$cv, tolerance = 1e-12)
  expect_true(all(w$.cost <= 450))

  three <- .dom_three()
  sep <- n_cluster(indicators = three, stage_cost = c(200, 40, 5),
                   domains = "region", budget = 3000, fixed_cost = 100)
  for (d in c("A", "B")) {
    alone <- n_cluster(indicators = three[three$region == d, -6],
                       stage_cost = c(200, 40, 5), budget = 3000,
                       fixed_cost = 100)
    row <- sep$operational$domains[sep$operational$domains$region == d, ]
    expect_equal(c(row$n_psu, row$n_per_psu, row$n_per_ssu),
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
  fit <- n_cluster(indicators = .dom_three(), stage_cost = c(200, 40, 5),
                   domains = "region", allocation = "joint", budget = 6000)
  w <- fit$operational$domains
  expect_equal(w$.cost, .whole_cost(w, c(200, 40, 5)))
  expect_lte(sum(w$.cost), 6000)
  expect_equal(w$.total_n, w$n_psu * w$n_per_psu * w$n_per_ssu)
})

test_that("the print shows the whole domain designs and their cost", {
  fit <- n_cluster(indicators = .dom_tab(), stage_cost = c(100, 10),
                   domains = "region", allocation = "joint", budget = 450,
                   n_psu = 2, fpc = "stage")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Total n = 4, cost = ", out, fixed = TRUE)))
  expect_false(any(grepl("unrounded", out, fixed = TRUE)))
  rows <- grep("^ [AB] ", out, value = TRUE)
  expect_length(rows, 2L)
  expect_true(all(grepl("^ [AB] +2 +1 +2 +1\\.0795 +220 ", rows)))
  expect_match(format(fit), "total_n=4 (unrounded: 5)", fixed = TRUE)

  sep <- n_cluster(indicators = .dom_tab(), stage_cost = c(100, 10),
                   domains = "region", budget = 450, n_psu = 2,
                   fpc = "stage")
  expect_true(any(grepl("budget = 450 for each domain, total cost = 880",
                        capture.output(print(sep)), fixed = TRUE)))
})

test_that("joint domain shares reserve the cost of a whole design at the floor", {
  # a floor of 41 on PSUs of 20 units needs 3 PSUs whole (cost 342 at take
  # 14) against 2.05 continuous (246), and a share between the two used to
  # be rejected although the budget paid for both domains
  tab <- data.frame(name = "y", region = c("A", "B"), p = c(0.2, 0.5),
                    cv = c(0.15, 0.2), icc_psu = 0.4, N = 400, N_psu = 20)
  fit_at <- function(budget) {
    n_cluster(indicators = tab, stage_cost = c(100, 1), domains = "region",
              allocation = "joint", budget = budget, min_n_domain = 41,
              fpc = "stage")
  }
  fit <- fit_at(690)
  d <- fit$operational$domains
  expect_true(all(d$.total_n >= 41))
  expect_lte(fit$operational$cost, 690)
  expect_equal(d$n_psu, c(3, 3))

  # within 1 percent of the best whole design by enumeration
  cvw <- function(a, b, p) {
    prec_cluster(n = c(a, b), icc = 0.4, unit_relvar = (1 - p) / p,
                 N = 400, N_psu = 20, fpc = "stage")$cv
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

  expect_error(fit_at(683), "together at a cost of 684")
  expect_equal(
    svyplan:::.floor_least_whole_cost(41, c(100, 1, 0), 1, 20, c(1, 20),
                                      c(1, 1)),
    342
  )
})

test_that("a fixed PSU count above what the floor needs still prices the floor", {
  tab <- data.frame(name = "y", region = c("A", "B"), p = c(0.2, 0.5),
                    cv = c(0.15, 0.2), icc_psu = 0.1, N = 400, N_psu = 20)
  fit_at <- function(...) {
    n_cluster(indicators = tab, stage_cost = c(100, 1), domains = "region",
              allocation = "joint", budget = 1200, n_psu = 5, fpc = "stage",
              ...)
  }
  free <- fit_at()
  floored <- fit_at(min_n_domain = 1)
  expect_equal(floored$operational$domains$.total_n,
               free$operational$domains$.total_n)
  expect_true(all(fit_at(min_n_domain = 50)$operational$domains$.total_n >= 50))

  least <- svyplan:::.floor_least_whole_cost
  # five PSUs fixed, a take of one: 5 * (100 + 1)
  expect_equal(least(1, c(100, 1, 0), 5, 5, c(1, 20), c(1, 1)), 505)
  # an SSU take below its box is never priced
  expect_equal(least(1, c(100, 1, 5), 1, 20, c(1, 20), c(2, 10)), 111)
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
    expect_equal(least(m, co, alo, ahi, c(blo, bhi), c(clo, chi)), ref)
  }
})

## The joint split

.dom_ratio <- function(fit, tab, stage_cost, ...) {
  vapply(seq_len(nrow(fit$domains)), function(d) {
    sub <- tab[tab$region == fit$domains$region[d], names(tab) != "region"]
    one <- suppressWarnings(n_cluster(stage_cost = stage_cost, indicators = sub,
                                      budget = fit$domains$.cost[d], ...))
    max(one$detail$.cv_achieved / one$detail$.cv_target)
  }, numeric(1L))
}

.dom_four <- function() {
  data.frame(
    name = paste0(rep(c("y", "z"), 4), rep(1:4, each = 2)),
    p = c(0.3, 0.1, 0.25, 0.4, 0.5, 0.2, 0.15, 0.35),
    cv = c(0.1, 0.15, 0.12, 0.1, 0.08, 0.2, 0.15, 0.1),
    icc_psu = c(0.05, 0.02, 0.1, 0.08, 0.03, 0.06, 0.12, 0.04),
    region = rep(c("A", "B", "C", "D"), each = 2)
  )
}

test_that("a joint split without correction is proportional to each domain's cost", {
  tab <- .dom_four()[1:6, ]
  cost <- c(400, 40)
  fit <- n_cluster(indicators = tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 1e5)
  ## The ratio falls as 1 / sqrt(budget), so b * ratio^2 is each domain's
  ## cost at its targets, read here from a fit at the whole budget.
  k <- vapply(c("A", "B", "C"), function(lev) {
    one <- n_cluster(indicators = tab[tab$region == lev, -5], stage_cost = cost,
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
                     N_psu = rep(c(50, 250), each = 2))
  cases <- list(
    list(tab = tab, args = list()),
    list(tab = tab, args = list(n_psu = 20)),
    list(tab = stage, args = list(fpc = "stage"))
  )
  for (cs in cases) {
    fit <- do.call(n_cluster, c(list(indicators = cs$tab, stage_cost = cost,
                                     domains = "region",
                                     allocation = "joint", budget = 6e4),
                                cs$args))
    worst <- function(w) {
      max(vapply(1:2, function(d) {
        sub <- cs$tab[cs$tab$region == c("A", "B")[d], names(cs$tab) != "region"]
        one <- suppressWarnings(do.call(n_cluster, c(
          list(stage_cost = cost, indicators = sub,
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
  tab <- transform(.dom_four(), N_psu = rep(c(60, 150, 400, 1000), each = 2))
  tab$N <- tab$N_psu * 120
  cost <- c(400, 40)
  fit <- n_cluster(indicators = tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 1.5e5, fpc = "stage")
  r <- .dom_ratio(fit, tab, cost, fpc = "stage")
  expect_equal(r, rep(max(r), 4), tolerance = 1e-7)
  expect_lte(sum(fit$domains$.cost), 1.5e5 * (1 + 1e-12))
})

test_that("a joint split keeps a whole design for a domain that needs almost nothing", {
  tab <- data.frame(name = paste0("y", 1:3), p = c(0.3, 0.3, 0.99999),
                    cv = 0.1, icc_psu = 0.05, region = c("A", "B", "C"))
  expect_warning(
    fit <- n_cluster(indicators = tab, stage_cost = c(500, 50),
                     domains = "region", allocation = "joint", budget = 1e5),
    "single PSU"
  )
  expect_gte(fit$domains$.cost[3], 550)
  w <- fit$operational$domains
  expect_gte(w$n_psu[3], 1)
  expect_lte(sum(w$.cost), 1e5)
  r <- .dom_ratio(fit, tab, c(500, 50))
  expect_equal(r[1], r[2], tolerance = 1e-8)
  expect_lt(r[3], r[1])
})

test_that("a joint floor is held and the other domains share the worst ratio", {
  tab <- .dom_four()[1:6, ]
  cost <- c(400, 40)
  fit <- n_cluster(indicators = tab, stage_cost = cost, domains = "region",
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
  tab$icc_ssu <- 0.1
  fit <- n_cluster(indicators = tab, stage_cost = cost, domains = "region",
                   allocation = "joint", budget = 2e5)
  rev_tab <- tab[nrow(tab):1, ]
  back <- n_cluster(indicators = rev_tab, stage_cost = cost, domains = "region",
                    allocation = "joint", budget = 2e5)
  expect_equal(back$domains$.cost[match(fit$domains$region, back$domains$region)],
               fit$domains$.cost, tolerance = 1e-8)
  ten <- n_cluster(indicators = tab, stage_cost = 10 * cost, domains = "region",
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
             cv = c(0.1, 0.15, 0.2, 0.2, 0.3, 0.3), icc_psu = 0.05,
             region = rep(c("A", "B", "C"), each = 2))
}

test_that("a fixed PSU count that rounds down floors each share at its continuous cost", {
  ## 2.4 PSUs at a take of 1 cost 264, and the 2 PSUs fielded cost 220.
  tab <- .dom_frac()
  fit <- n_cluster(indicators = tab, stage_cost = c(100, 10), budget = 800,
                   n_psu = 2.4, domains = "region", allocation = "joint")
  dom <- fit$domains
  expect_equal(sum(dom$.cost), 800, tolerance = 1e-10)
  expect_true(all(dom$.cost >= 264 * (1 - 1e-12)))
  for (lev in dom$region) {
    at <- dom$region == lev
    one <- n_cluster(indicators = tab[tab$region == lev, 1:4],
                     stage_cost = c(100, 10), budget = dom$.cost[at],
                     n_psu = 2.4)
    expect_equal(dom$n_per_psu[at], one$n[["n_per_psu"]], tolerance = 1e-8)
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
  fit <- n_cluster(indicators = .dom_frac(), stage_cost = c(100, 10),
                   budget = 792, n_psu = 2.4, domains = "region",
                   allocation = "joint")
  expect_equal(fit$domains$.cost, rep(264, 3), tolerance = 1e-10)
  expect_equal(fit$domains$n_per_psu, rep(1, 3), tolerance = 1e-10)
  expect_error(
    n_cluster(indicators = .dom_frac(), stage_cost = c(100, 10),
              budget = 792 * (1 - 1e-6), n_psu = 2.4, domains = "region",
              allocation = "joint"),
    "costs 792 plus fixed_cost = 0, counting each domain at the larger cost"
  )
})

test_that("a fixed PSU count that rounds up floors each share at its whole cost", {
  ## 2.6 PSUs at a take of 1 cost 286, and the 3 PSUs fielded cost 330.
  tab <- .dom_tab(levels = c("A", "B"))[, 1:5]
  expect_error(
    n_cluster(indicators = tab, stage_cost = c(100, 10), budget = 659,
              n_psu = 2.6, domains = "region", allocation = "joint"),
    "costs 660 plus"
  )
  fit <- n_cluster(indicators = tab, stage_cost = c(100, 10), budget = 660,
                   n_psu = 2.6, domains = "region", allocation = "joint")
  expect_identical(fit$operational$domains$n_psu, c(3L, 3L))
  expect_lte(fit$operational$cost, 660)
})

test_that("fractional fixed takes under the stage correction fit at their continuous floor", {
  ## One PSU at a take of 2.4 costs 124 and fields a take of 2 at 120.
  tab <- .dom_tab(levels = c("A", "B"))
  args <- list(indicators = tab, stage_cost = c(100, 10), n_per_psu = 2.4,
               domains = "region", allocation = "joint", fpc = "stage")
  fits <- function(...) suppressWarnings(n_cluster(...))
  fit <- do.call(fits, c(args, budget = 248))
  expect_equal(fit$domains$.cost, c(124, 124), tolerance = 1e-10)
  expect_lte(fit$operational$cost, 248)
  expect_no_error(do.call(fits, c(args, budget = 260)))
  expect_error(do.call(fits, c(args, budget = 247.99)), "costs 248 plus")

  ## Three stages at takes 2.4 and 3.4 cost 132.16 per PSU, 126 whole. At
  ## the exact share the SSU take bound is met to rounding only.
  tab3 <- cbind(tab, icc_ssu = 0.1, N_ssu = 1000)
  args3 <- list(indicators = tab3, stage_cost = c(100, 10, 1),
                n_per_psu = 2.4, n_per_ssu = 3.4, domains = "region",
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
  tab3 <- cbind(.dom_tab(levels = c("A", "B"))[, 1:5], icc_ssu = 0.1)
  args <- list(indicators = tab3, stage_cost = c(100, 10, 1), n_psu = 2.4,
               domains = "region", allocation = "joint")
  fit <- do.call(n_cluster, c(args, budget = 532.8))
  expect_equal(fit$domains$.cost, c(266.4, 266.4), tolerance = 1e-10)
  expect_error(do.call(n_cluster, c(args, budget = 532.7)), "costs 532.8 plus")
})

test_that("a domain floor is reached under a fractional fixed PSU count", {
  tab <- .dom_tab(levels = c("A", "B"))[, 1:5]
  for (np in c(2.4, 2.6)) {
    fit <- n_cluster(indicators = tab, stage_cost = c(100, 10), budget = 3000,
                     n_psu = np, min_n_domain = 100, domains = "region",
                     allocation = "joint")
    expect_true(all(fit$domains$.total_n >= 100 - 1e-8))
    expect_true(all(fit$operational$domains$.total_n >= 100))
    expect_lte(fit$operational$cost, 3000)
  }
  ## 2.4 PSUs reach 100 at 240 + 10 * 100 = 1240, and the 2 PSUs fielded
  ## at a take of 50 cost 1200.
  args <- list(indicators = tab, stage_cost = c(100, 10), n_psu = 2.4,
               min_n_domain = 100, domains = "region", allocation = "joint")
  fit <- do.call(n_cluster, c(args, budget = 2480))
  expect_equal(fit$domains$.cost, c(1240, 1240), tolerance = 1e-10)
  expect_error(do.call(n_cluster, c(args, budget = 2479)),
               "together at a cost of 2480")
  ## A fixed take of 2.4 fields 2, and 10 PSUs reach 20.
  fit <- n_cluster(indicators = .dom_tab(levels = c("A", "B")),
                   stage_cost = c(100, 10), budget = 5000, n_per_psu = 2.4,
                   min_n_domain = 20, domains = "region",
                   allocation = "joint", fpc = "stage")
  expect_true(all(fit$operational$domains$.total_n >= 20))
})

test_that("free PSUs without the stage correction keep only the whole floor", {
  ## One PSU at a fielded take of 2 costs 120, below the 124 a continuous
  ## take of 2.4 would need on a whole PSU.
  tab <- .dom_tab(levels = c("A", "B"))[, 1:5]
  args <- list(indicators = tab, stage_cost = c(100, 10), n_per_psu = 2.4,
               domains = "region", allocation = "joint")
  fits <- function(...) suppressWarnings(n_cluster(...))
  fit <- do.call(fits, c(args, budget = 240))
  expect_equal(fit$domains$.cost, c(120, 120), tolerance = 1e-10)
  expect_error(do.call(fits, c(args, budget = 239.99)), "costs 240 plus")
})

## Floors against the solvers

.dom_two <- function(...) {
  data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
             region = c("A", "B"), ...)
}

## Each domain's whole design equals the scalar fit of that domain at its
## share, so the joint floor let through a share the solver can field.
.expect_domains_as_scalar <- function(fit, tab, ...) {
  for (lev in fit$domains$region) {
    share <- fit$domains$.cost[fit$domains$region == lev]
    one <- suppressWarnings(n_cluster(
      indicators = tab[tab$region == lev, names(tab) != "region"],
      budget = share, ...
    ))
    op <- fit$operational$domains[fit$operational$domains$region == lev, ]
    expect_equal(unname(unlist(op[names(one$operational$n)])),
                 unname(one$operational$n))
  }
}

test_that("a domain floor already met keeps the fractional-PSU allowance", {
  ## Without the stage correction free PSUs go below one, so a floor of
  ## one unit must not cost more than the fit without it.
  tab <- .dom_two()
  args <- list(indicators = tab, stage_cost = c(100, 10), budget = 240,
               n_per_psu = 2.4, domains = "region", allocation = "joint")
  plain <- suppressWarnings(do.call(n_cluster, args))
  floored <- suppressWarnings(do.call(n_cluster, c(args, min_n_domain = 1)))
  expect_equal(floored$domains$.cost, plain$domains$.cost)
  expect_lt(max(floored$domains$n_psu), 1)
  .expect_domains_as_scalar(floored, tab, stage_cost = c(100, 10),
                            n_per_psu = 2.4)

  tab3 <- .dom_two(icc_ssu = 0.1)
  args3 <- list(indicators = tab3, stage_cost = c(100, 10, 1), budget = 252,
                n_per_psu = 2.4, n_per_ssu = 3.4, domains = "region",
                allocation = "joint")
  plain3 <- suppressWarnings(do.call(n_cluster, args3))
  floored3 <- suppressWarnings(do.call(n_cluster, c(args3, min_n_domain = 1)))
  expect_equal(floored3$domains$.cost, plain3$domains$.cost)

  ## A floor of 30 at a fixed take of 2.4 and a free SSU take. One PSU at
  ## the whole takes 2 and 15 costs 150 and reaches 30, and the continuous
  ## design reaches it on that share with fewer than one PSU.
  args30 <- list(indicators = tab3, stage_cost = c(100, 10, 1),
                 n_per_psu = 2.4, min_n_domain = 30, domains = "region",
                 allocation = "joint")
  fit <- suppressWarnings(do.call(n_cluster, c(args30, budget = 300)))
  expect_true(all(fit$domains$.total_n >= 30 - 1e-8))
  expect_true(all(fit$operational$domains$.total_n >= 30))
  expect_error(suppressWarnings(do.call(n_cluster, c(args30, budget = 299))),
               "together at a cost of 300")
})

test_that("a fixed PSU count below one fields one PSU with or without a floor", {
  tab <- .dom_two()
  args <- list(indicators = tab, stage_cost = c(100, 10), budget = 220,
               n_psu = 0.4, domains = "region", allocation = "joint")
  plain <- suppressWarnings(do.call(n_cluster, args))
  floored <- suppressWarnings(do.call(n_cluster, c(args, min_n_domain = 1)))
  expect_equal(floored$domains, plain$domains)
  expect_identical(floored$operational$domains$n_psu, c(1L, 1L))
  .expect_domains_as_scalar(floored, tab, stage_cost = c(100, 10),
                            n_psu = 0.4)
})

test_that("a fixed take above the whole population cap is priced at the cap", {
  ## PSUs of 2.8 units field a take of 2 for a fixed 2.6, so each domain's
  ## floor is the continuous 126, not 130 for a take of 3.
  tab <- .dom_two(N = 28, N_psu = 10)
  args <- list(indicators = tab, stage_cost = c(100, 10), n_per_psu = 2.6,
               domains = "region", allocation = "joint", fpc = "stage")
  fit <- suppressWarnings(do.call(n_cluster, c(args, budget = 252)))
  expect_identical(fit$operational$domains$n_per_psu, c(2L, 2L))
  .expect_domains_as_scalar(fit, tab, stage_cost = c(100, 10),
                            n_per_psu = 2.6, fpc = "stage")
  expect_error(suppressWarnings(do.call(n_cluster, c(args, budget = 251.9))),
               "costs 252 plus")
})

test_that("free PSUs under the stage correction keep at least one PSU at a floor", {
  ## The continuous design holds one PSU at least, and at a share of 150 it
  ## reaches 26 of the 30. The least share reaching 30 is 154.
  tab <- .dom_two(icc_ssu = 0.1, N = 1e6, N_psu = 100, N_ssu = 2000)
  args <- list(indicators = tab, stage_cost = c(100, 10, 1), n_per_psu = 2.4,
               min_n_domain = 30, domains = "region", allocation = "joint",
               fpc = "stage")
  fit <- suppressWarnings(do.call(n_cluster, c(args, budget = 308)))
  expect_true(all(fit$domains$n_psu >= 1 - 1e-12))
  expect_true(all(fit$domains$.total_n >= 30 - 1e-8))
  expect_error(suppressWarnings(do.call(n_cluster, c(args, budget = 307))),
               "together at a cost of 308")
})
