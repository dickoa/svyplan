## The menu listing against brute force over every part

list_register <- function(seed) {
  set.seed(seed)
  H <- sample(2:3, 1)
  take <- sample(c(5, 10), 1)
  psu <- do.call(rbind, lapply(seq_len(H), function(h) {
    k <- sample(c(1, 4, 6, 8), 1)
    sizes <- pmax(round(stats::rlnorm(k, log(20), 0.9)), take)
    if (k > 4) sizes[1:2] <- sizes[1:2] * 4
    data.frame(stratum = LETTERS[h], N = sizes)
  }))
  frame <- data.frame(
    stratum = LETTERS[seq_len(H)],
    N = as.numeric(tapply(psu$N, psu$stratum, sum)),
    n_per_psu = take
  )
  m <- c(0L, 1L, 2L)[seed %% 3 + 1]
  if (m > 0 && seed %% 4 == 0) psu$zone_order <- sample(nrow(psu))
  if (seed %% 5 == 0) {
    psu$certainty <- seq_len(nrow(psu)) %in% sample(nrow(psu), H)
  }
  if (seed %% 2 == 0) {
    frame$cost_psu <- 50
    frame$cost_ssu <- 5
  }
  if (seed %% 6 == 1) frame$max_weight <- c(4, rep(NA, H - 1))
  z <- list(
    frame = frame, psu = psu,
    measures = data.frame(stratum = frame$stratum, name = "y",
                          p = stats::runif(H, 0.1, 0.6),
                          icc_psu = stats::runif(H, 0.02, 0.2)),
    targets = data.frame(name = "y", cv = 0.1)
  )
  ev <- .psu_setup(
    z$frame, .check_psu_table(z$psu, z$frame, z$measures), z$measures,
    z$targets, min_n_stratum = if (seed %% 3 == 2) min(frame$N),
    certainty_cutoff = if (seed %% 2) 0.8, n_psu_per_zone = if (m > 0) m,
    resp_rate = 0.9
  )
  list(z = z, ev = ev)
}

# Every state of every part from a fraction near 0 to 1, by the R closure
# at every event and midpoint, kept where it reaches the lower bound on
# units. Part t forces the t largest PSUs the caller did not flag, and
# under unit cost only part 0 is planned.
list_oracle <- function(ev, h) {
  M <- ev$sizes[[h]]
  free <- which(!ev$flagged[[h]])
  cand <- free[order(-M[free], free)]
  parts <- if (is.null(ev$cost_psu)) 0L else 0:length(cand)
  keys <- unlist(lapply(parts, function(t) {
    forced <- ev$flagged[[h]]
    forced[cand[seq_len(t)]] <- TRUE
    menu_oracle(M, ev$N[h], ev$take[h], 1e-9, 1, ev$cutoff[h], forced,
                ev$zone_m, ev$zone_key[[h]])
  }))
  keys <- unique(keys)
  parts <- strsplit(keys, "|", fixed = TRUE)
  n <- vapply(parts, function(x) as.numeric(x[2]) + as.numeric(x[3]) *
                ev$take[h], numeric(1))
  keys[n >= ev$lower[h] - 1e-8]
}

reduced <- function(menu, h, weights) {
  mh <- menu[[h]]
  weights$price * mh$states$cost + drop(mh$C %*% weights$mu) +
    weights$obj * mh$states$objective
}

test_that("a listing without a stop gives every state of every part", {
  for (seed in 1:12) {
    r <- list_register(seed)
    weights <- list(price = 1, mu = 1, obj = 0)
    menu <- .psu_menu_round(.psu_menu_new(r$ev), r$ev, weights,
                            rep(Inf, length(r$ev$N)))
    for (h in seq_along(menu)) {
      expect_setequal(menu[[h]]$states$key, list_oracle(r$ev, h))
      expect_false(any(menu[[h]]$open))
    }
    expect_false(.psu_menu_short(menu, r$ev, weights,
                                 rep(Inf, length(r$ev$N))))
  }
})

test_that("a stopped listing leaves out only states above the limit", {
  # Each limit is the reduced cost of a state, some of them in the part
  # where every PSU is held, where the bound is nearly tight.
  cut <- 0
  for (seed in 1:12) {
    r <- list_register(seed)
    H <- length(r$ev$N)
    set.seed(seed)
    weights <- list(price = 1, mu = stats::runif(1, 100, 5000), obj = 0)
    full <- .psu_menu_round(.psu_menu_new(r$ev), r$ev, weights, rep(Inf, H))
    rc <- lapply(seq_len(H), function(h) reduced(full, h, weights))
    pick <- lapply(seq_len(H), function(h) {
      held <- which(full[[h]]$states$n_c == length(r$ev$sizes[[h]]))
      c(sample(length(rc[[h]]), min(3, length(rc[[h]]))),
        held[order(rc[[h]][held])][seq_len(min(2, length(held)))])
    })
    for (j in seq_len(max(lengths(pick)))) {
      limit <- vapply(seq_len(H), function(h) {
        rc[[h]][pick[[h]][min(j, length(pick[[h]]))]]
      }, numeric(1))
      menu <- .psu_menu_round(.psu_menu_new(r$ev), r$ev, weights, limit)
      for (h in seq_len(H)) {
        wanted <- full[[h]]$states$key[rc[[h]] <= limit[h]]
        expect_true(all(wanted %in% menu[[h]]$states$key))
        expect_true(all(menu[[h]]$states$key %in% full[[h]]$states$key))
        cut <- cut + NROW(full[[h]]$states) - NROW(menu[[h]]$states)
      }
      expect_false(.psu_menu_short(menu, r$ev, weights, limit))
    }
  }
  expect_gt(cut, 0)
})

test_that("a resumed listing equals one listed at the higher limit", {
  for (seed in 1:12) {
    r <- list_register(seed)
    H <- length(r$ev$N)
    weights <- list(price = 1, mu = 800, obj = 0)
    full <- .psu_menu_round(.psu_menu_new(r$ev), r$ev, weights, rep(Inf, H))
    rc <- lapply(seq_len(H), function(h) reduced(full, h, weights))
    low <- vapply(rc, min, numeric(1))
    high <- vapply(rc, function(x) stats::quantile(x, 0.4), numeric(1))
    menu <- .psu_menu_round(.psu_menu_new(r$ev), r$ev, weights, low)
    menu <- .psu_menu_round(menu, r$ev, weights, high)
    once <- .psu_menu_round(.psu_menu_new(r$ev), r$ev, weights, high)
    for (h in seq_len(H)) {
      expect_setequal(menu[[h]]$states$key, once[[h]]$states$key)
    }
  }
})

test_that("the stop rule holds for objective rows and a budget price", {
  z <- .psu_fixture(4)
  z$psu <- z$psu[z$psu$stratum %in% c("C", "D"), ]
  z$psu <- do.call(rbind, lapply(split(z$psu, z$psu$stratum),
                                 function(d) d[1:12, ]))
  z$frame <- z$frame[3:4, ]
  z$frame$N <- as.numeric(tapply(z$psu$N, z$psu$stratum, sum))
  z$frame$cost_psu <- 80
  z$frame$cost_ssu <- 8
  z$measures <- z$measures[z$measures$stratum %in% c("C", "D"), ]
  ev <- .psu_setup(z$frame, z$psu, z$measures, z$targets,
                   objective = c("literacy", "income"))
  weights <- list(price = 2e-6, mu = c(0.003, 0.001), obj = 1)
  full <- .psu_menu_round(.psu_menu_new(ev), ev, weights, c(Inf, Inf))
  rc <- lapply(1:2, function(h) reduced(full, h, weights))
  limit <- vapply(rc, function(x) stats::quantile(x, 0.2), numeric(1))
  menu <- .psu_menu_round(.psu_menu_new(ev), ev, weights, limit)
  for (h in 1:2) {
    wanted <- full[[h]]$states$key[rc[[h]] <= limit[h]]
    expect_true(all(wanted %in% menu[[h]]$states$key))
    expect_lt(nrow(menu[[h]]$states), nrow(full[[h]]$states))
  }
})

test_that("a stratum the caller holds whole carries no remainder visit", {
  # One flagged PSU takes j units at f = j / M, so its variance is the
  # stratum-rate one and the bound is exact there.
  frame <- data.frame(stratum = c("A", "B"), N = c(400, 900), n_per_psu = 50,
                      cost_psu = 60, cost_ssu = 1)
  psu <- data.frame(stratum = c("A", "B", "B"), N = c(400, 600, 300),
                    certainty = c(TRUE, FALSE, FALSE))
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = 0.3,
                         icc_psu = 0.1)
  ev <- .psu_setup(frame, psu, measures, data.frame(name = "y", cv = 0.1))
  weights <- list(price = 1, mu = 0.05, obj = 0)
  full <- .psu_menu_round(.psu_menu_new(ev), ev, weights, c(Inf, Inf))
  rc <- reduced(full, 1, weights)
  for (j in 4:49) {
    key <- sprintf("1|%d|0", j)
    menu <- .psu_menu_round(.psu_menu_new(ev), ev, weights,
                            c(rc[full[[1]]$states$key == key], Inf))
    expect_true(key %in% menu[[1]]$states$key)
  }
})
