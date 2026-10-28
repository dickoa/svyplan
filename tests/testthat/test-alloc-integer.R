## The classic allocation's whole design against every whole allocation

classic_frame <- function(seed) {
  set.seed(seed)
  H <- sample(2:3, 1)
  data.frame(N = sample(20:70, H), sd = stats::runif(H, 1, 9),
             mean = stats::runif(H, 8, 12),
             unit_cost = round(stats::runif(H, 1, 10), 1))
}

# Every whole allocation within the fit's bounds, with its cost, sample
# and variance, sum W^2 S^2 deff (1 / (r n) - 1 / N).
classic_plans <- function(fit, deff, resp) {
  d <- fit$detail
  grid <- as.matrix(expand.grid(lapply(seq_len(nrow(d)), function(h) {
    max(1, ceiling(d$.lower[h] - 1e-9)):floor(d$.upper[h] + 1e-9)
  })))
  W <- d$N / sum(d$N)
  a <- W^2 * d$sd^2 * deff
  V <- drop((1 / grid) %*% (a / resp)) - sum(a / d$N)
  list(cost = drop(grid %*% d$unit_cost), n = rowSums(grid), V = V,
       ybar = sum(W * d$mean))
}

test_that("a classic cv allocation is the cheapest whole allocation", {
  for (seed in 1:15) {
    fr <- classic_frame(seed)
    deff <- stats::runif(nrow(fr), 1, 2)
    resp <- stats::runif(nrow(fr), 0.6, 1)
    for (alloc in c("optimal", "neyman")) {
      fit <- n_alloc(fr, cv = 0.12, alloc = alloc, deff = deff,
                     resp_rate = resp)
      plans <- classic_plans(fit, deff, resp)
      ok <- sqrt(plans$V) / plans$ybar <= 0.12
      label <- paste("seed", seed, alloc)
      if (alloc == "optimal") {
        expect_equal(fit$operational$cost, min(plans$cost[ok]),
                     tolerance = 1e-9, label = label)
      } else {
        expect_equal(fit$operational$n, min(plans$n[ok]), label = label)
      }
      expect_lte(fit$operational$cv, 0.12)
      expect_identical(fit$operational$status, "optimal", label = label)
    }
  }
})

test_that("a classic budget buys the smallest whole variance", {
  for (seed in 1:15) {
    fr <- classic_frame(seed)
    deff <- stats::runif(nrow(fr), 1, 2)
    resp <- stats::runif(nrow(fr), 0.6, 1)
    budget <- 0.5 * sum(fr$N * fr$unit_cost)
    fit <- n_alloc(fr, budget = budget, alloc = "optimal", deff = deff,
                   resp_rate = resp)
    plans <- classic_plans(fit, deff, resp)
    ok <- plans$cost <= budget + 1e-9
    expect_equal(fit$operational$se^2, min(plans$V[ok]), tolerance = 1e-9,
                 label = paste("seed", seed))
    expect_lte(fit$operational$cost, budget)
  }
})

test_that("rules and domains keep their rounding", {
  fr <- data.frame(N = rep(100, 3), sd = 1:3, mean = 10)
  # The search finds 16 units where rounding each stratum up gives 17.
  expect_identical(n_alloc(fr, cv = 0.05)$detail$n_int, c(3L, 5L, 8L))
  prop <- n_alloc(fr, cv = 0.05, alloc = "proportional")
  expect_identical(prop$detail$n_int, as.integer(ceiling(prop$detail$n - 1e-9)))
  expect_null(prop$operational$status)
  fr$region <- c("a", "a", "b")
  dom <- n_alloc(fr, cv = 0.05, domains = "region")
  expect_identical(dom$detail$n_int, as.integer(ceiling(dom$detail$n - 1e-9)))
  expect_null(dom$optimization$search)
  expect_error(n_alloc(fr[, -5], n = 30, control = list(max_nodes = 5)),
               "'control' applies where n_alloc\\(\\) searches")
})

test_that("a classic search at its node limit keeps its limit everywhere", {
  set.seed(3)
  fr <- data.frame(N = sample(200:900, 12), sd = stats::runif(12, 1, 9),
                   mean = 10, unit_cost = round(stats::runif(12, 1, 10), 2))
  fit <- n_alloc(fr, cv = 0.02, alloc = "optimal",
                 control = list(max_nodes = 1))
  expect_identical(fit$operational$status, "limit")
  expect_lte(fit$optimization$search$bound, fit$operational$cost)
  expect_true(any(grepl("^search: stopped at the node limit",
                        utils::capture.output(print(fit)))))
  again <- n_alloc(prec_alloc(fit), cv = 0.02)
  expect_identical(again$optimization$search$max_nodes, 1)
  # A rebuild in n mode does not search, so it carries no limit.
  expect_null(n_alloc(prec_alloc(fit))$optimization$search)
})

## The two-stage field design against every whole design

cluster_frame <- function(seed, fixed = FALSE) {
  set.seed(seed)
  H <- sample(2:3, 1)
  N <- sample(150:400, H)
  fr <- data.frame(N = N, N_psu = sample(8:16, H), icc_psu = stats::runif(H, 0.02, 0.2),
                   sd = stats::runif(H, 1, 9), mean = 10,
                   cost_psu = sample(c(50, 100, 200), H, TRUE),
                   cost_ssu = sample(c(5, 10, 20), H, TRUE))
  if (fixed) fr$n_per_psu <- sample(4:9, H, TRUE)
  fr
}

# Every state (a PSUs at take b) of every stratum, with its cost and its
# variance term W^2 S^2 deff (icc b r + (1 - f) (1 - icc)) (1 / (a b r) -
# 1 / N_fpc), f the within-PSU fraction of fpc = "stage".
cluster_states <- function(fr, fit, deff, resp, fpc, takes = NULL) {
  W <- fr$N / sum(fr$N)
  do.call(rbind, lapply(seq_len(nrow(fr)), function(h) {
    # A free take stays within the average PSU size.
    b_set <- if (is.null(takes)) seq_len(max(1, floor(fr$N[h] / fr$N_psu[h]))) else takes[h]
    g <- expand.grid(a = seq_len(fr$N_psu[h]), b = b_set)
    g <- g[g$a * g$b <= fr$N[h] & g$a * g$b >= ceiling(fit$detail$.lower[h] - 1e-9), ]
    f <- if (fpc == "stage") pmin(1, g$b * resp[h] / (fr$N[h] / fr$N_psu[h])) else 0
    infl <- fr$icc_psu[h] * g$b * resp[h] + (1 - f) * (1 - fr$icc_psu[h])
    K <- W[h]^2 * fr$sd[h]^2 * deff[h] * infl
    D <- if (fpc == "stage") K * fr$icc_psu[h] / (fr$N_psu[h] * infl) else K / fr$N[h]
    data.frame(h = h, a = g$a, b = g$b,
               cost = g$a * (fr$cost_psu[h] + fr$cost_ssu[h] * g$b),
               n = g$a * g$b, V = K / (g$a * g$b * resp[h]) - D)
  }))
}

# The best plan of up to three strata: for each plan of all strata but the
# last, the last stratum's best state within what is left.
best_plan <- function(st, objective, limit, cap) {
  parts <- split(st, st$h)
  head_plans <- parts[[1]][, c(objective, limit)]
  if (length(parts) == 3L) {
    p2 <- parts[[2]]
    idx <- expand.grid(i = seq_len(nrow(head_plans)), j = seq_len(nrow(p2)))
    head_plans <- data.frame(head_plans[idx$i, objective] + p2[idx$j, objective],
                             head_plans[idx$i, limit] + p2[idx$j, limit])
    names(head_plans) <- c(objective, limit)
  }
  last <- parts[[length(parts)]]
  last <- last[order(last[[limit]]), ]
  best_by_limit <- cummin(last[[objective]])
  k <- findInterval(cap - head_plans[[limit]] + 1e-12 * abs(cap),
                    last[[limit]])
  ok <- k > 0
  min(head_plans[[objective]][ok] + best_by_limit[k[ok]])
}

test_that("a two-stage cv design is the cheapest whole design", {
  for (seed in 1:12) {
    fr <- cluster_frame(seed)
    H <- nrow(fr)
    fpc <- if (seed %% 3 == 0) "stage" else "unit"
    deff <- stats::runif(H, 1, 1.5)
    resp <- if (seed %% 2 == 0) stats::runif(H, 0.7, 1) else rep(1, H)
    fit <- suppressWarnings(n_alloc(fr, cv = 0.12, alloc = "optimal", deff = deff,
                                    resp_rate = resp, fpc = fpc))
    st <- cluster_states(fr, fit, deff, resp, fpc)
    target <- (0.12 * 10)^2
    label <- paste("seed", seed, fpc)
    expect_equal(fit$operational$cost, best_plan(st, "cost", "V", target),
                 tolerance = 1e-9, label = label)
    expect_lte(fit$operational$cv, 0.12)
    expect_identical(fit$operational$status, "optimal", label = label)
    expect_identical(fit$optimization$search$family, "whole units")
  }
})

test_that("a two-stage budget buys the smallest whole variance", {
  for (seed in 1:12) {
    fr <- cluster_frame(seed)
    H <- nrow(fr)
    fpc <- if (seed %% 3 == 0) "stage" else "unit"
    budget <- sample(c(800, 1500, 3000), 1)
    fit <- suppressWarnings(n_alloc(fr, budget = budget, alloc = "optimal",
                                    fpc = fpc))
    st <- cluster_states(fr, fit, rep(1, H), rep(1, H), fpc)
    expect_equal(fit$operational$se^2, best_plan(st, "V", "cost", budget),
                 tolerance = 1e-9, label = paste("seed", seed, fpc))
    expect_lte(fit$operational$cost, budget)
  }
})

test_that("a fixed take stays fixed and neyman keeps the rule's takes", {
  for (seed in 1:8) {
    fr <- cluster_frame(seed, fixed = TRUE)
    H <- nrow(fr)
    fit <- suppressWarnings(n_alloc(fr, cv = 0.12, alloc = "optimal"))
    expect_identical(fit$detail$n_per_psu_int, as.integer(fr$n_per_psu))
    st <- cluster_states(fr, fit, rep(1, H), rep(1, H), "unit",
                         takes = fr$n_per_psu)
    expect_equal(fit$operational$cost, best_plan(st, "cost", "V", 1.2^2),
                 tolerance = 1e-9, label = paste("seed", seed))
    fr$n_per_psu <- NULL
    ney <- suppressWarnings(n_alloc(fr, cv = 0.12, alloc = "neyman"))
    take <- ney$detail$n_per_psu
    expect_true(all(ney$detail$n_per_psu_int %in% c(floor(take), ceiling(take))))
    st <- cluster_states(fr, ney, rep(1, H), rep(1, H), "unit",
                         takes = ney$detail$n_per_psu_int)
    expect_equal(ney$operational$n, best_plan(st, "n", "V", 1.2^2),
                 label = paste("seed", seed, "neyman"))
  }
})

test_that("a capped PSU count raises the take where that is cheaper", {
  fr <- data.frame(N = c(4000, 6000), N_psu = c(18, 120), icc_psu = 0.05,
                   sd = c(5, 3), mean = 10, cost_psu = 200, cost_ssu = 10)
  fit <- n_alloc(fr, budget = 20000, alloc = "optimal")
  # The cost-optimal take is 19, and the capped stratum takes 35 at its 18.
  expect_identical(fit$detail$n_psu_int[1], 18L)
  expect_gt(fit$detail$n_per_psu_int[1], round(fit$detail$n_per_psu[1]))
  again <- n_alloc(prec_alloc(fit), budget = 20000)
  expect_identical(again$detail$n_per_psu_int, fit$detail$n_per_psu_int)
  limited <- n_alloc(fr, cv = 0.02, alloc = "optimal",
                     control = list(max_nodes = 1))
  expect_identical(n_alloc(prec_alloc(limited), cv = 0.02)$optimization$search$max_nodes, 1)
  expect_null(n_alloc(fr, n = 300)$optimization$search)
})

## A free take past the PSU universe

test_that("a free take rises past the PSU universe to the cheapest design", {
  for (seed in 1:12) {
    fr <- cluster_frame(seed)
    H <- nrow(fr)
    set.seed(seed + 100)
    fr$N_psu <- sample(2:4, H, TRUE)
    fpc <- if (seed %% 3 == 0) "stage" else "unit"
    label <- paste("seed", seed, fpc)
    fit <- suppressWarnings(n_alloc(fr, cv = 0.08, alloc = "optimal",
                                    fpc = fpc))
    st <- cluster_states(fr, fit, rep(1, H), rep(1, H), fpc)
    expect_equal(fit$operational$cost, best_plan(st, "cost", "V", 0.8^2),
                 tolerance = 1e-9, label = label)
    expect_true(all(fit$detail$n_psu <= fr$N_psu + 1e-9), label = label)
    budget <- if (seed %% 2 == 0) 1500 else 3000
    fb <- suppressWarnings(n_alloc(fr, budget = budget, alloc = "optimal",
                                   fpc = fpc))
    st <- cluster_states(fr, fb, rep(1, H), rep(1, H), fpc)
    expect_equal(fb$operational$se^2, best_plan(st, "V", "cost", budget),
                 tolerance = 1e-9, label = label)
    expect_equal(prec_alloc(fb)$cv, fb$cv, tolerance = 1e-10, label = label)
  }
})

test_that("the rule's continuous design fields every PSU at a risen take", {
  fr <- data.frame(N = 1000, N_psu = 5, icc_psu = 0.05, sd = 5, mean = 10,
                   cost_psu = 100, cost_ssu = 10)
  fit <- n_alloc(fr, cv = 0.05, alloc = "optimal")
  expect_equal(fit$detail$n_psu, 5)
  expect_equal(fit$detail$n_per_psu, fit$detail$n / 5)
  expect_equal(fit$cv, 0.05, tolerance = 1e-8)
  expect_equal(prec_alloc(fit)$cv, fit$cv, tolerance = 1e-10)
  expect_lt(fit$operational$cost,
            n_alloc(transform(fr, n_per_psu = 100), cv = 0.05,
                    alloc = "optimal")$operational$cost)
  fb <- n_alloc(fr, budget = 2000, alloc = "optimal")
  expect_equal(fb$operational$cost, 2000)
  expect_equal(fb$detail$n_per_psu_int, 30L)
  expect_equal(prec_alloc(fr, n = 400)$detail$n_psu, 5)
  expect_error(n_alloc(fr, cv = 0.05, alloc = "neyman"),
               "At the take 'n_per_psu' this allocation holds")
})

## Targets near zero variance

test_that("a target near zero variance gets the census, not a refusal", {
  fr <- data.frame(N = c(40, 20, 40), sd = c(11.7, 5.9, 11.7),
                   mean = c(20.5, 50.5, 80.5))
  for (cv in c(1e-7, 1e-12)) {
    fit <- n_alloc(fr, cv = cv)
    expect_identical(fit$detail$n_int, c(40L, 20L, 40L), label = cv)
    expect_equal(fit$operational$cv, 0)
    sb <- strata_bound(1:100, n_strata = 3, cv = cv, method = "cumrootf")
    expect_equal(sb$n, 100)
    expect_equal(sb$cv, 0)
  }
  expect_error(n_alloc(fr, cv = 1e-7, resp_rate = 0.9), "unattainable")
  cl <- data.frame(N = c(400, 600), N_psu = c(20, 30), icc_psu = 0.05,
                   sd = c(5, 3), mean = 10, cost_psu = 200, cost_ssu = 10)
  fit <- n_alloc(cl, cv = 1e-7, alloc = "optimal")
  expect_identical(fit$detail$n_psu_int, c(20L, 30L))
  expect_identical(fit$detail$n_per_psu_int, c(20L, 20L))
})

test_that("a target just above the smallest attainable cv is met", {
  fr <- data.frame(N = 100, N_psu = 10, n_per_psu = 10, icc_psu = 0.05,
                   sd = 1, mean = 10, cost_psu = 1, cost_ssu = 1)
  # Nonresponse keeps the variance of every PSU at its full take positive.
  for (f in list(fr, fr[names(fr) != "n_per_psu"])) {
    floor_cv <- prec_alloc(f, n = 100, resp_rate = 0.8)$cv
    fit <- n_alloc(f, cv = floor_cv * (1 + 1e-10), alloc = "optimal",
                   resp_rate = 0.8)
    expect_identical(fit$detail$n_psu_int, 10L)
    expect_identical(fit$detail$n_per_psu_int, 10L)
    expect_lte(fit$operational$cv, floor_cv * (1 + 1e-10))
    expect_error(n_alloc(f, cv = floor_cv * (1 - 1e-6), alloc = "optimal",
                         resp_rate = 0.8))
  }
  el <- data.frame(N = c(40, 60), sd = c(3, 5), mean = 10)
  floor_cv <- prec_alloc(el, n = c(40, 60), resp_rate = 0.8)$cv
  fit <- n_alloc(el, cv = floor_cv * (1 + 1e-10), resp_rate = 0.8)
  expect_identical(fit$detail$n_int, c(40L, 60L))
})

test_that("a target a rounding step past the floor is met, a larger one refused", {
  fr <- data.frame(N = 100, N_psu = 10, n_per_psu = 10, icc_psu = 0.05,
                   sd = 1, mean = 1, cost_psu = 1, cost_ssu = 1)
  floor_cv <- prec_alloc(fr, n = 100, resp_rate = 0.8)$cv
  el <- data.frame(N = c(40, 60), sd = c(3, 5), mean = 10)
  el_cv <- prec_alloc(el, n = c(40, 60), resp_rate = 0.8)$cv
  below <- 1 - .Machine$double.eps
  fit <- n_alloc(fr, cv = floor_cv * below, alloc = "optimal",
                 resp_rate = 0.8)
  expect_identical(fit$detail$n_psu_int, 10L)
  expect_identical(n_alloc(el, cv = el_cv * below, resp_rate = 0.8)$detail$n_int,
                   c(40L, 60L))
  # Past the rounding tolerance but within the continuous one, the
  # whole-unit step refuses in its own words.
  expect_error(n_alloc(fr, cv = floor_cv * (1 - 1e-10), alloc = "optimal",
                       resp_rate = 0.8),
               "no whole-cluster design found that attains target CV")
  expect_error(n_alloc(el, cv = el_cv * (1 - 1e-10), resp_rate = 0.8),
               "unattainable with whole units")
  el$unit_cost <- c(1, 2)
  census <- sum(el$N * el$unit_cost)
  expect_identical(n_alloc(el, budget = census + 1e-9, alloc = "optimal")$detail$n_int,
                   c(40L, 60L))
})

test_that("a target at the floor takes each stratum's cheapest state there", {
  fr <- data.frame(N = c(100, 100), N_psu = 10, n_per_psu = 10,
                   icc_psu = 0.05, sd = c(1, 0), mean = 1, cost_psu = 1,
                   cost_ssu = 1)
  floor_cv <- prec_alloc(fr, n = c(100, 100), resp_rate = 0.8)$cv
  # A stratum without variance is at its floor at its fewest PSUs, here 2
  # of 10 under min_n_stratum = 20.
  for (k in c(1, 1 + 1e-10)) {
    fit <- n_alloc(fr, cv = floor_cv * k, alloc = "optimal", resp_rate = 0.8,
                   min_n_stratum = 20)
    expect_identical(fit$detail$n_psu_int, c(10L, 2L), label = k)
    expect_equal(fit$operational$cost, 132)
    expect_equal(fit$optimization$search$bound, 132)
    expect_equal(fit$optimization$search$gap, 0)
    expect_identical(fit$operational$status, "optimal")
    expect_lte(fit$operational$cv, floor_cv * k)
  }
})

test_that("a target assessed at the floor is met despite rounding", {
  # Each frame's floor, read by the search, sits a few rounding steps above
  # the cv prec_alloc() assesses for the census.
  for (seed in 6:10) {
    set.seed(seed)
    H <- sample(2:3, 1)
    b <- sample(3:12, H)
    P <- sample(4:12, H)
    sd <- stats::runif(H, 0.5, 5)
    sd[sample(H, 1)] <- if (seed %% 2) 0 else sd[1]
    fr <- data.frame(N = b * P, N_psu = P, n_per_psu = b,
                     icc_psu = stats::runif(H, 0.02, 0.2), sd = sd,
                     mean = stats::runif(H, 1, 10),
                     cost_psu = round(stats::runif(H, 1, 50)),
                     cost_ssu = round(stats::runif(H, 1, 5)))
    resp <- stats::runif(1, 0.6, 0.95)
    fpc <- sample(c("unit", "stage"), 1)
    mins <- sample(1:3, 1) * max(b)
    cv <- prec_alloc(fr, n = fr$N, resp_rate = resp, fpc = fpc)$cv
    fit <- function(k) {
      n_alloc(fr, cv = cv * k, alloc = "optimal", resp_rate = resp,
              fpc = fpc, min_n_stratum = mins)
    }
    at <- fit(1)
    label <- paste("seed", seed)
    expect_equal(at$operational$cost, fit(1 + 1e-10)$operational$cost,
                 label = label)
    expect_equal(at$optimization$search$bound, at$operational$cost,
                 label = label)
    expect_lte(at$operational$cv, cv * (1 + 1e-12))
  }
})

test_that("a target just above a near-census plan finds that plan", {
  for (seed in 1:6) {
    fr <- classic_frame(seed)
    fit <- n_alloc(fr, cv = 0.12, alloc = "optimal")
    plans <- classic_plans(fit, 1, 1)
    # The smallest variance short of zero, with every plan at or below it.
    v <- min(plans$V[plans$V > 1e-12 * max(plans$V)])
    cv <- sqrt(v) / plans$ybar * (1 + 1e-8)
    fit <- n_alloc(fr, cv = cv, alloc = "optimal")
    ok <- sqrt(pmax(plans$V, 0)) / plans$ybar <= cv
    expect_equal(fit$operational$cost, min(plans$cost[ok]), tolerance = 1e-9,
                 label = paste("seed", seed))
    expect_lte(fit$operational$cv, cv)
  }
})

## The take within the PSU size, and the n mode at a risen take

test_that("a free take stays within the average PSU size", {
  cl <- data.frame(N = c(400, 600), N_psu = c(20, 30), icc_psu = 0.05,
                   sd = c(5, 3), mean = 10, cost_psu = 200, cost_ssu = 10)
  # Under the unit correction one PSU of 400 units was a census.
  fit <- suppressWarnings(n_alloc(cl, cv = 0.01, alloc = "optimal"))
  expect_true(all(fit$detail$n_per_psu_int <= cl$N / cl$N_psu))
  for (seed in 1:10) {
    set.seed(seed)
    H <- sample(2:4, 1)
    N <- sample(200:3000, H)
    fr <- data.frame(N = N, N_psu = pmax(2, round(N / stats::runif(H, 10, 60))),
                     icc_psu = stats::runif(H, 0.02, 0.2), sd = stats::runif(H, 2, 8),
                     mean = 10, cost_psu = stats::runif(H, 50, 500),
                     cost_ssu = stats::runif(H, 2, 20))
    fit <- tryCatch(suppressWarnings(n_alloc(fr, cv = stats::runif(1, 0.02, 0.1),
                                             alloc = "optimal")),
                    error = function(e) NULL)
    if (is.null(fit)) next
    expect_true(all(fit$detail$n_per_psu_int <= pmax(1, floor(fr$N / fr$N_psu))),
                label = paste("seed", seed))
  }
  # Without N_psu no PSU size bounds the take, so it stays at the rule's.
  cl$N_psu <- NULL
  free <- suppressWarnings(n_alloc(cl, cv = 0.01, alloc = "optimal"))
  take <- free$detail$n_per_psu
  expect_true(all(free$detail$n_per_psu_int %in% c(floor(take), ceiling(take))))
})

test_that("the n mode raises a free take, so the round trip holds", {
  fr <- data.frame(N = 1000, N_psu = 5, icc_psu = 0.05, sd = 5, mean = 10,
                   cost_psu = 100, cost_ssu = 10)
  for (fit in list(n_alloc(fr, cv = 0.05, alloc = "optimal"),
                   n_alloc(fr, budget = 2000, alloc = "optimal"))) {
    again <- n_alloc(prec_alloc(fit))
    expect_identical(again$detail$n_psu_int, fit$detail$n_psu_int)
    expect_identical(again$detail$n_per_psu_int, fit$detail$n_per_psu_int)
    expect_equal(again$operational$cost, fit$operational$cost)
  }
  at <- n_alloc(fr, n = 400, alloc = "optimal")
  expect_equal(at$cv, prec_alloc(fr, n = 400)$cv, tolerance = 1e-10)
  expect_identical(at$detail$n_psu_int, 5L)
  expect_identical(at$detail$n_per_psu_int, 80L)
  expect_true(any(grepl("^continuous design: n = ",
                        utils::capture.output(print(at)))))
  out <- utils::capture.output(print(summary(at)))
  expect_true("Continuous design" %in% out)
  expect_false(any(grepl("optimum", out, ignore.case = TRUE)))
})

test_that("a frame without variance takes its lower bounds", {
  z <- data.frame(N = c(100, 200), sd = 0, mean = 10)
  for (fit in list(
    suppressWarnings(n_alloc(z, cv = 0.05)),
    suppressWarnings(n_alloc(z, budget = 50, alloc = "optimal"))
  )) {
    expect_identical(fit$detail$n_int, c(1L, 1L))
    expect_identical(fit$operational$cv, 0)
    expect_identical(fit$optimization$search$status, "optimal")
    expect_identical(fit$optimization$search$nodes, 0)
  }
})

## The search's status within a gap
##
## A design within 1e-6 of its bound at the node limit is reported "gap"
## rather than "limit", and control$gap stops the search there on purpose.

.gap_frame <- function() {
  set.seed(5)
  H <- 150
  data.frame(stratum = paste0("s", seq_len(H)), N = sample(200:60000, H),
             sd = runif(H, 5, 40), mean = runif(H, 20, 80),
             unit_cost = runif(H, 1, 10))
}

test_that("a stop within 1e-6 of the bound reads gap, a wider one limit", {
  near <- n_alloc(.gap_frame(), cv = 0.005, alloc = "optimal",
                  control = list(max_nodes = 1000))
  s <- near$optimization$search
  expect_identical(s$status, "gap")
  expect_lte(s$gap, 1e-6)
  expect_identical(near$operational$status, "gap")
  expect_output(print(near), "search: within")
  far <- n_alloc(.gap_frame(), cv = 0.005, alloc = "optimal",
                 control = list(max_nodes = 50))
  expect_identical(far$optimization$search$status, "limit")
  expect_gt(far$optimization$search$gap, 1e-6)
})

test_that("control$gap stops the search within the gap it sets", {
  fit <- n_alloc(.gap_frame(), cv = 0.005, alloc = "optimal",
                 control = list(gap = 1e-5))
  s <- fit$optimization$search
  expect_identical(s$status, "gap")
  expect_lte(s$gap, 1e-5)
  expect_lt(s$nodes, 5000)
  expect_lte(fit$operational$cv, 0.005)
  expect_error(n_alloc(.gap_frame(), cv = 0.005, alloc = "optimal",
                       control = list(gap = 1)),
               "must be one number in \\[0, 1\\)")
  expect_error(n_alloc(.gap_frame(), cv = 0.005, alloc = "optimal",
                       control = list(tolerance = 1)),
               "the known ones are 'max_nodes' and 'gap'")
})
