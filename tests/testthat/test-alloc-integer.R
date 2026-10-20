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
    b_set <- if (is.null(takes)) seq_len(fr$N[h]) else takes[h]
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
