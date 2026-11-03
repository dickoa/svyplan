## The cv mode with domains, each domain at its own minimum

ten_regions <- function() {
  set.seed(4)
  R <- 10
  fr <- data.frame(region = rep(sprintf("R%02d", 1:R), each = 2),
                   area = rep(c("urban", "rural"), R))
  fr$stratum <- paste(fr$region, fr$area)
  fr$N <- round(rep(exp(seq(log(4e5), log(2e4), length.out = R)), each = 2) *
                  rep(c(0.4, 0.6), R))
  fr$p <- round(runif(2 * R, 0.15, 0.6), 2)
  fr
}

per_domain <- function(fr, dom, ...) {
  parts <- lapply(split(fr, fr[[dom]]), function(f) n_alloc(f, ...))
  list(n = sum(vapply(parts, function(x) x$n, numeric(1))),
       n_int = sum(vapply(parts, function(x) sum(x$detail$n_int), numeric(1))))
}

## C1. The continuous design

test_that("each domain sits at its own smallest total", {
  fr <- ten_regions()
  for (rule in c("neyman", "proportional", "power")) {
    fit <- n_alloc(fr, domains = "region", cv = 0.05, alloc = rule)
    expect_equal(fit$n, per_domain(fr, "region", cv = 0.05, alloc = rule)$n,
                 info = rule)
    expect_equal(fit$domains$.cv, rep(0.05, 10), tolerance = 1e-6,
                 info = rule)
  }
})

test_that("the domain minimum is the joint mode's minimum", {
  fr <- ten_regions()
  fit <- n_alloc(fr, domains = "region", cv = 0.05)
  joint <- n_alloc(fr[, c("stratum", "region", "N")],
                   measures = data.frame(stratum = fr$stratum, name = "y",
                                         p = fr$p),
                   targets = data.frame(name = "y", domain = "region",
                                        cv = 0.05))
  expect_equal(fit$n, joint$n, tolerance = 1e-6)
  expect_equal(sum(fit$detail$n_int), sum(joint$detail$n_int))
})

test_that("the documented example needs 104 units", {
  fd <- data.frame(province = c("North", "North", "South", "South"),
                   stratum = c("Urban", "Rural", "Urban", "Rural"),
                   N = c(2000, 3000, 1800, 3200), sd = c(12, 18, 10, 16),
                   mean = c(55, 48, 58, 50))
  fit <- n_alloc(fd, domains = "province", cv = 0.04, alloc = "power",
                 alloc_q = 0.3)
  expect_identical(sum(fit$detail$n_int), 104L)
})

test_that("a floor raises only the domains the target leaves below it", {
  fr <- data.frame(stratum = paste(rep(c("A", "B", "C"), each = 2),
                                   c("u", "r")),
                   province = rep(c("A", "B", "C"), each = 2),
                   N = c(20000, 30000, 8000, 12000, 3000, 4000),
                   p = c(0.2, 0.35, 0.25, 0.4, 0.3, 0.45))
  free <- n_alloc(fr, domains = "province", cv = 0.09)
  held <- n_alloc(fr, domains = "province", cv = 0.09, min_n_domain = 250)
  low <- free$domains$.n < 250
  expect_equal(held$domains$.n[low], rep(250, sum(low)))
  expect_equal(held$domains$.n[!low], free$domains$.n[!low])
  expect_identical(held$domains$.floor_binding, low)
})

test_that("a domain that cannot reach the target is named", {
  # A census always reaches a target, so the cap here is the PSU universe.
  fr <- data.frame(stratum = c("a", "b", "c"), region = c("N", "N", "S"),
                   N = c(3000, 5000, 2000), N_psu = c(100, 100, 3),
                   p = c(0.2, 0.4, 0.3), icc_psu = 0.1, n_per_psu = 10,
                   cost_psu = 300, cost_ssu = 20)
  expect_error(suppressWarnings(n_alloc(fr, domains = "region", cv = 0.05)),
               "unattainable in domain")
})

test_that("a cluster frame takes each domain's minimum", {
  fr <- data.frame(stratum = paste(rep(c("A", "B", "C"), each = 2),
                                   c("u", "r")),
                   province = rep(c("A", "B", "C"), each = 2),
                   N = c(20000, 30000, 8000, 12000, 3000, 4000),
                   N_psu = c(200, 300, 80, 120, 30, 40),
                   p = c(0.2, 0.35, 0.25, 0.4, 0.3, 0.45),
                   icc_psu = 0.05, n_per_psu = 12, cost_psu = 300,
                   cost_ssu = 20)
  fit <- suppressWarnings(n_alloc(fr, domains = "province", cv = 0.09))
  own <- vapply(split(fr, fr$province), function(f) {
    suppressWarnings(n_alloc(f, cv = 0.09))$n
  }, numeric(1))
  expect_equal(fit$domains$.n, unname(own))
  expect_true(all(as.numeric(tapply(fit$detail$n_int, fr$province, sum)) >=
                    unname(own) - 1e-9))
})

## C2. The whole design

test_that("each domain's whole design is its cheapest", {
  fr <- data.frame(stratum = c("a", "b", "c", "d"),
                   region = c("N", "N", "S", "S"),
                   N = c(400, 600, 300, 500), sd = c(4, 9, 6, 2),
                   mean = c(20, 25, 30, 22), unit_cost = c(1, 2, 1, 3))
  fit <- n_alloc(fr, domains = "region", cv = 0.03, alloc = "optimal")
  # Exhaustive: the cheapest whole pair in each domain meeting its cv.
  for (r in c("N", "S")) {
    idx <- which(fr$region == r)
    W <- fr$N[idx] / sum(fr$N[idx])
    ybar <- sum(W * fr$mean[idx])
    best <- Inf
    for (a in 1:120) for (b in 1:120) {
      n <- c(a, b)
      v <- sum(W^2 * fr$sd[idx]^2 * (1 / n - 1 / fr$N[idx]))
      if (sqrt(v) / ybar <= 0.03 + 1e-12) {
        best <- min(best, sum(fr$unit_cost[idx] * n))
      }
    }
    expect_equal(sum(fr$unit_cost[idx] * fit$detail$n_int[idx]), best,
                 info = r)
  }
})

test_that("a domain its stratum minimums already satisfy stays at them", {
  fr <- data.frame(stratum = c("a", "b", "c", "d"),
                   region = c("N", "N", "S", "S"),
                   N = c(4000, 6000, 3000, 5000), p = c(0.5, 0.5, 0.5, 0.5),
                   min_n_stratum = c(1, 1, 200, 200))
  fit <- n_alloc(fr, domains = "region", cv = 0.08)
  expect_equal(fit$domains$.n[2L], 400)
  expect_lt(fit$domains$.cv[2L], 0.08)
  expect_equal(fit$domains$.cv[1L], 0.08, tolerance = 1e-6)
})

## C3. The whole-cluster design per domain

cluster_regions <- function(free_take = FALSE) {
  set.seed(11)
  R <- 3
  fr <- data.frame(region = rep(sprintf("R%02d", 1:R), each = 2),
                   area = rep(c("urban", "rural"), R))
  fr$stratum <- paste(fr$region, fr$area)
  fr$N <- round(rep(c(1e5, 5e4, 2e4), each = 2) * rep(c(0.4, 0.6), R))
  fr$N_psu <- round(fr$N / 120)
  fr$p <- round(runif(2 * R, 0.15, 0.6), 2)
  fr$icc_psu <- round(runif(2 * R, 0.02, 0.15), 3)
  if (!free_take) fr$n_per_psu <- 15
  fr$cost_psu <- round(runif(2 * R, 200, 600))
  fr$cost_ssu <- round(runif(2 * R, 10, 30))
  fr
}

expect_domains_searched <- function(fr, cv, alloc) {
  fit <- suppressWarnings(n_alloc(fr, domains = "region", cv = cv,
                                  alloc = alloc))
  for (r in unique(fr$region)) {
    idx <- which(fr$region == r)
    own <- suppressWarnings(n_alloc(fr[idx, ], cv = cv, alloc = alloc))
    expect_identical(fit$detail$n_psu_int[idx], own$detail$n_psu_int,
                     info = paste(alloc, r))
    expect_identical(fit$detail$n_per_psu_int[idx],
                     own$detail$n_per_psu_int, info = paste(alloc, r))
  }
  expect_named(fit$optimization$domains, names(fit$params$domain_idx))
  fit
}

test_that("each domain's whole cluster design is the one it gets alone", {
  for (alloc in c("neyman", "optimal")) {
    expect_domains_searched(cluster_regions(), 0.08, alloc)
  }
  # A free take under "optimal" may rise past N_psu * n_per_psu, as the
  # national search lets it.
  expect_domains_searched(cluster_regions(free_take = TRUE), 0.06, "optimal")
})

test_that("a free take rises past N_psu * n_per_psu in a domain search", {
  fr <- cluster_regions(free_take = TRUE)
  fr$N_psu <- c(12, 15, 8, 10, 5, 6)
  fit <- expect_domains_searched(fr, 0.05, "optimal")
  # The PSU universe binds, so the take rises above the rule's take.
  expect_true(any(fit$detail$n_psu_int == fr$N_psu))
})

test_that("a domain mixing fixed and free takes is searched as one", {
  fr <- cluster_regions(free_take = TRUE)
  fr$n_per_psu <- c(15, NA, NA, 20, NA, NA)
  expect_domains_searched(fr, 0.06, "optimal")
})

test_that("a floor still holds the cluster domains it binds", {
  fr <- cluster_regions()
  free <- suppressWarnings(n_alloc(fr, domains = "region", cv = 0.08))
  held <- suppressWarnings(n_alloc(fr, domains = "region", cv = 0.08,
                                   min_n_domain = 700))
  totals <- as.numeric(tapply(held$detail$n_int, fr$region, sum))
  expect_true(all(totals >= 700))
  low <- free$domains$.n < 700
  expect_identical(held$domains$.floor_binding, low)
  # A domain the floor leaves alone keeps its searched design.
  keep <- which(fr$region %in% unique(fr$region)[!low])
  expect_identical(held$detail$n_psu_int[keep], free$detail$n_psu_int[keep])
  # A held domain keeps the rule's shares at the floor, each stratum within
  # one PSU of its continuous size.
  floored <- which(fr$region %in% unique(fr$region)[low])
  d <- held$detail[floored, ]
  expect_true(all(abs(d$n_psu_int * d$n_per_psu_int - d$n) <=
                    d$n_per_psu_int + 1e-9))
})
