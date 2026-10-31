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

test_that("a floor no whole design in the share reaches is refused", {
  expect_error(
    n_cluster(indicators = .dom_tab(), stage_cost = c(100, 10),
              domains = "region", allocation = "joint", budget = 2420,
              n_psu = 2, min_n_domain = 101, fpc = "stage"),
    "domain 'A': no whole-unit design within its share"
  )

  ## A take of 10.5 per SSU reaches 42, and the whole take of 10 gives 40.
  three <- data.frame(name = "y", p = 0.3, cv = 0.1, icc_psu = 0.05,
                      icc_ssu = 0.1, region = c("A", "B"))
  expect_error(
    n_cluster(indicators = three, stage_cost = c(200, 40, 5),
              domains = "region", allocation = "joint", budget = 1540,
              n_psu = 2, n_per_psu = 2, min_n_domain = 42),
    "domain 'A': no whole-unit design within its share"
  )
  room <- n_cluster(indicators = three, stage_cost = c(200, 40, 5),
                    domains = "region", allocation = "joint", budget = 1560,
                    n_psu = 2, n_per_psu = 2, min_n_domain = 42)
  expect_equal(room$operational$domains$n_per_ssu, c(11L, 11L))
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
