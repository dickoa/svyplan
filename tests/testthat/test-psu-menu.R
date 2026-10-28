## The menu sweep in C against the closure the R code computes

test_that("each part lists every state the R closure gives, and only those", {
  for (seed in 1:18) {
    r <- menu_register(seed)
    m <- c(0L, 1L, 2L)[seed %% 3 + 1]
    cutoff <- if (seed %% 2) 1 else 0.8
    zone_key <- if (m > 0 && seed %% 4 == 0) sample(length(r$sizes))
    N <- sum(r$sizes)
    o <- order(-r$sizes)
    for (t in c(0L, 1L, 3L)) {
      forced <- logical(length(r$sizes))
      forced[o[seq_len(t)]] <- TRUE
      hi <- if (t > 0L) min(0.5, cutoff * r$take / r$sizes[o[t]]) else 0.5
      got <- .psu_menu_part(r$sizes, N, r$take, 1e-9, hi, cutoff, forced, m,
                            zone_key)$states$key
      want <- menu_oracle(r$sizes, N, r$take, 1e-9, hi, cutoff, forced, m,
                          zone_key)
      expect_setequal(got, want)
    }
  }
})

test_that("the sweep reproduces the R closure on sizes that are not whole", {
  sizes <- c(410.5, 96.25, 61.75, 40.5, 33.25, 22.75)
  N <- sum(sizes)
  for (m in 0:2) {
    got <- .psu_menu_part(sizes, N, 5, 1e-9, 0.6, m = m)$states$key
    expect_setequal(got, menu_oracle(sizes, N, 5, 1e-9, 0.6, m = m))
  }
})

test_that("each state carries the sums its variance is built from", {
  r <- menu_register(3)
  N <- sum(r$sizes)
  part <- .psu_menu_part(r$sizes, N, r$take, 1e-9, 0.5, m = 2L)
  s <- part$states
  for (i in seq_len(nrow(s))) {
    held <- part$sets[, s$set[i]]
    q <- pmin(ceiling(s$f[i] * N / N * r$sizes[held]), r$sizes[held])
    expect_equal(s$Q[i], sum(q))
    expect_equal(s$A[i], sum(r$sizes[held]^2 / q))
    expect_equal(s$held[i], sum(r$sizes[held]))
    expect_equal(s$n_c[i], sum(held))
  }
  # A drawn remainder comes in whole zones of two, unless it is all drawn.
  rest <- vapply(s$set, function(j) sum(!part$sets[, j]), 0)
  expect_true(all(s$a %% 2 == 0 | s$a == rest))
})

test_that("one PSU, and a stratum held whole, give one draw each", {
  one <- .psu_menu_part(500, 500, 10, 1e-9, 1)$states
  expect_setequal(one$key, menu_oracle(500, 500, 10, 1e-9, 1))
  expect_true(all(one$set == one$set[1]))

  sizes <- c(120, 80, 60)
  all_held <- .psu_menu_part(sizes, 260, 10, 1e-9, 1, forced = rep(TRUE, 3))
  expect_true(all(all_held$states$a == 0))
  expect_identical(ncol(all_held$sets), 1L)
})

test_that("the sweep stops at the first event at or above f_stop", {
  r <- menu_register(5)
  N <- sum(r$sizes)
  full <- .psu_menu_part(r$sizes, N, r$take, 1e-9, 0.5)
  cut <- .psu_menu_part(r$sizes, N, r$take, 1e-9, 0.5, f_stop = 0.2)
  expect_false(full$stopped)
  expect_true(cut$stopped)
  expect_gte(cut$f_end, 0.2)
  # Every state met below the stop is kept, none past the stop's event.
  expect_true(all(full$states$key[full$states$f < 0.2] %in% cut$states$key))
  expect_true(all(cut$states$f <= cut$f_end))
})

test_that("the sweep refuses settings it cannot list", {
  expect_error(.psu_menu_part(c(10, Inf), 20, 5, 1e-9, 0.5), "invalid PSU")
  expect_error(.psu_menu_part(c(10, 10), 20, 5, 0.6, 0.5), "invalid part")
  expect_error(.psu_menu_part(c(10, 10), 20, 5, 1e-9, 0.5, cutoff = 1.2),
               "invalid part")
})

test_that("every listed fraction lies in the part's interval", {
  s <- .psu_menu_part(c(100, 80, 60), 240, 10, 0.001, 0.011,
                      forced = c(TRUE, FALSE, FALSE))$states
  expect_true(all(s$f >= 0.001 & s$f <= 0.011))
  set.seed(11)
  for (i in 1:60) {
    M <- pmax(round(stats::rlnorm(sample(3:9, 1), log(60), 0.8)), 10)
    lo <- stats::runif(1, 0.001, 0.05)
    hi <- min(1, lo * stats::runif(1, 1.2, 6))
    forced <- logical(length(M))
    if (i %% 2) forced[which.max(M)] <- TRUE
    s <- .psu_menu_part(M, sum(M), 10, lo, hi, forced = forced,
                        m = c(0L, 1L, 2L)[i %% 3 + 1])$states
    expect_true(all(s$f >= lo & s$f <= hi), label = paste("part", i))
  }
})

test_that("the sweep keeps a state that two close events bound", {
  sizes <- c(999999, rep(1, 160001))
  N <- sum(sizes)
  # A take event at 517240 / 999999 and a draw event at 599999 / N lie
  # 8.6e-13 apart, which a window of 1e-11 steps over.
  f <- 517240 / 999999
  got <- .psu_menu_part(sizes, N, 1, f - 5e-6, f + 5e-6)$states$key
  want <- menu_oracle(sizes, N, 1, f - 5e-6, f + 5e-6)
  expect_setequal(got, want)
  expect_true("1|517241|82758" %in% got)
})

test_that("the sweep does not emit consecutive copies of a state", {
  sizes <- c(30, 30, 40)
  out <- .Call(
    C_psu_menu_part, sizes, 100, 10, 0.001, 1, 1,
    logical(3), 0L, as.integer(order(-sizes) - 1L), Inf, TRUE
  )
  s <- out$states
  same <- s$set[-1] == s$set[-length(s$set)] &
    s$Q[-1] == s$Q[-length(s$Q)] & s$a[-1] == s$a[-length(s$a)]
  expect_false(any(same))
})
