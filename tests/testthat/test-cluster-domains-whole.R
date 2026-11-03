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
  expect_true(any(grepl("Total n = 4 (unrounded: 5)", out, fixed = TRUE)))
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
