## The searched whole design of a classic allocation with three stages

without_singleton_psu_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("a single PSU, so no within-stratum variance", conditionMessage(w),
              fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

three_stage_points <- function(fr, h, W, r, fpc, a_max = 60, t2_max = 20,
                               t3_max = 25, lower = 1) {
  g <- expand.grid(a = seq_len(a_max), t2 = seq_len(t2_max),
                   t3 = seq.int(ceiling(1 / r - 1e-9), t3_max))
  n <- g$a * g$t2 * g$t3
  g <- g[n <= fr$N[h] & n >= lower, ]
  R1 <- fr$icc_1[h]
  R2 <- (1 - fr$icc_1[h]) * fr$icc_2[h]
  R3 <- (1 - fr$icc_1[h]) * (1 - fr$icc_2[h])
  core <- W[h]^2 * fr$sd[h]^2 *
    (R1 / g$a + R2 / (g$a * g$t2) + R3 / (g$a * g$t2 * g$t3 * r))
  infl <- R1 * g$t2 * g$t3 * r + R2 * g$t3 * r + R3
  g$var <- if (fpc == "unit") core - W[h]^2 * fr$sd[h]^2 * infl / fr$N[h] else core
  g$cost <- g$a * (fr$cost_1[h] + fr$cost_2[h] * g$t2 +
                     fr$cost_3[h] * g$t2 * g$t3)
  g <- g[order(g$cost, g$var), ]
  g[g$var < c(Inf, cummin(g$var)[-nrow(g)]), ]
}

three_stage_best <- function(p1, p2, v_max = NULL, budget = NULL) {
  best <- list(value = Inf)
  for (i in seq_len(nrow(p1))) {
    if (!is.null(v_max)) {
      room <- v_max - p1$var[i]
      j <- which(p2$var <= room * (1 + 1e-12))
      if (!length(j)) next
      j <- j[1L]
      value <- p1$cost[i] + p2$cost[j]
    } else {
      room <- budget - p1$cost[i]
      j <- which(p2$cost <= room + 1e-9)
      if (!length(j)) next
      j <- j[length(j)]
      value <- p1$var[i] + p2$var[j]
    }
    if (value < best$value) best <- list(value = value, i = i, j = j)
  }
  best
}

test_that("one stratum at three stages takes the design n_cluster() finds", {
  frame <- data.frame(stratum = "A", N = 1e6, N_1 = 1e4, N_2 = 1e5, mean = 50,
                      sd = 10, icc_1 = 0.05, icc_2 = 0.10, cost_1 = 500,
                      cost_2 = 100, cost_3 = 20)
  fit <- n_alloc(frame, cv = 0.03, alloc = "optimal", fpc = "none")
  expect_equal(unname(unlist(fit$detail[c("n_1_int", "take_2_int",
                                          "take_3_int")])), c(5, 4, 6))
  expect_equal(fit$operational$cost, 6900)
  expect_identical(fit$operational$status, "optimal")
})

test_that("three strata at three stages meet the cv at the searched cost", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(20000, 30000, 12000),
                      sd = c(10, 14, 8), mean = c(40, 50, 45),
                      icc_1 = c(0.05, 0.1, 0.08), icc_2 = c(0.15, 0.2, 0.1),
                      cost_1 = c(400, 600, 300), cost_2 = c(50, 40, 60),
                      cost_3 = c(10, 12, 8))
  fit <- n_alloc(frame, cv = 0.02, alloc = "optimal")
  expect_lte(fit$operational$cv, 0.02 + 1e-8)
  expect_lte(fit$operational$cost, 22700)
  expect_identical(fit$operational$status, "optimal")
  expect_identical(fit$optimization$search$status, "optimal")
})

test_that("the three-stage whole design matches enumeration over two strata", {
  set.seed(31)
  checked <- 0L
  for (i in seq_len(25)) {
    fpc <- if (i %% 2) "unit" else "none"
    r <- if (i %% 3 == 0) 0.8 else 1
    fr <- data.frame(
      stratum = c("A", "B"), N = round(runif(2, 5000, 50000)),
      sd = runif(2, 5, 20), mean = runif(2, 30, 60),
      icc_1 = runif(2, 0.02, 0.15), icc_2 = runif(2, 0.05, 0.3),
      cost_1 = runif(2, 100, 800), cost_2 = runif(2, 10, 80),
      cost_3 = runif(2, 2, 20)
    )
    W <- fr$N / sum(fr$N)
    ybar <- sum(W * fr$mean)
    cv <- runif(1, 0.03, 0.06)
    fit <- without_singleton_psu_warning(
      n_alloc(fr, cv = cv, alloc = "optimal", fpc = fpc, resp_rate = r)
    )
    info <- sprintf("instance %d", i)
    p1 <- three_stage_points(fr, 1, W, r, fpc, lower = fit$detail$.lower[1])
    p2 <- three_stage_points(fr, 2, W, r, fpc, lower = fit$detail$.lower[2])
    best <- three_stage_best(p1, p2, v_max = (cv * ybar)^2)
    expect_equal(fit$operational$cost, best$value, tolerance = 1e-9,
                 info = info)
    expect_identical(fit$operational$status, "optimal", info = info)

    budget <- fit$operational$cost * runif(1, 0.7, 1.1)
    fb <- without_singleton_psu_warning(
      n_alloc(fr, budget = budget, alloc = "optimal", fpc = fpc, resp_rate = r)
    )
    best_b <- three_stage_best(p1, p2, budget = budget)
    expect_equal((fb$operational$cv * ybar)^2, best_b$value,
                 tolerance = 1e-9, info = info)
    expect_lte(fb$operational$cost, budget + 1e-8)
    checked <- checked + 1L
  }
  expect_identical(checked, 25L)
})

test_that("a node limit at three stages is reported", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(20000, 30000, 12000),
                      sd = c(10, 14, 8), mean = c(40, 50, 45),
                      icc_1 = c(0.05, 0.1, 0.08), icc_2 = c(0.15, 0.2, 0.1),
                      cost_1 = c(400, 600, 300), cost_2 = c(50, 40, 60),
                      cost_3 = c(10, 12, 8))
  fit <- n_alloc(frame, cv = 0.01, alloc = "optimal",
                 control = list(max_nodes = 1))
  expect_identical(fit$operational$status, "limit")
  expect_lte(fit$operational$cv, 0.01 + 1e-8)
})

four_stage_points <- function(fr, h, W, r, fpc, lower) {
  g <- expand.grid(a = 1:40, t3 = 1:15, t4 = seq.int(ceiling(1 / r - 1e-9), 40))
  t2 <- fr$take_2[h]
  n <- g$a * t2 * g$t3 * g$t4
  g <- g[n <= fr$N[h] & n >= lower, ]
  icc <- c(fr$icc_1[h], fr$icc_2[h], fr$icc_3[h])
  R <- c(icc[1], (1 - icc[1]) * icc[2], (1 - icc[1]) * (1 - icc[2]) * icc[3],
         prod(1 - icc))
  core <- R[1] / g$a + R[2] / (g$a * t2) + R[3] / (g$a * t2 * g$t3) +
    R[4] / (g$a * t2 * g$t3 * g$t4 * r)
  infl <- R[1] * t2 * g$t3 * g$t4 * r + R[2] * g$t3 * g$t4 * r +
    R[3] * g$t4 * r + R[4]
  g$var <- W[h]^2 * fr$sd[h]^2 *
    (core - if (fpc == "unit") infl / fr$N[h] else 0)
  g$cost <- g$a * (fr$cost_1[h] + fr$cost_2[h] * t2 +
                     fr$cost_3[h] * t2 * g$t3 + fr$cost_4[h] * t2 * g$t3 * g$t4)
  g <- g[order(g$cost, g$var), ]
  g[g$var < c(Inf, cummin(g$var)[-nrow(g)]), ]
}

test_that("four-stage whole designs with a fixed take match enumeration", {
  set.seed(41)
  for (i in seq_len(6)) {
    fpc <- if (i %% 2) "unit" else "none"
    r <- if (i <= 3) 0.75 else 1
    fr <- data.frame(
      stratum = c("A", "B"), region = c("N", "S"),
      N = round(runif(2, 20000, 60000)), sd = runif(2, 5, 20),
      mean = runif(2, 30, 60), icc_1 = runif(2, 0.02, 0.1),
      icc_2 = runif(2, 0.05, 0.2), icc_3 = runif(2, 0.05, 0.2), take_2 = 3,
      cost_1 = runif(2, 300, 900), cost_2 = runif(2, 40, 120),
      cost_3 = runif(2, 10, 40), cost_4 = runif(2, 2, 10)
    )
    W <- fr$N / sum(fr$N)
    ybar <- sum(W * fr$mean)
    cv <- runif(1, 0.04, 0.07)
    info <- sprintf("instance %d", i)
    fit <- without_singleton_psu_warning(
      n_alloc(fr, cv = cv, alloc = "optimal", fpc = fpc, resp_rate = r)
    )
    p1 <- four_stage_points(fr, 1, W, r, fpc, fit$detail$.lower[1])
    p2 <- four_stage_points(fr, 2, W, r, fpc, fit$detail$.lower[2])
    best <- three_stage_best(p1, p2, v_max = (cv * ybar)^2)
    expect_equal(fit$operational$cost, best$value, tolerance = 1e-9,
                 info = info)
    expect_true(all(fit$detail$take_2_int == 3))

    dom <- without_singleton_psu_warning(
      n_alloc(fr, cv = cv, alloc = "optimal", fpc = fpc, resp_rate = r,
              domains = "region")
    )
    expect_identical(dom$operational$status, "optimal", info = info)
    for (h in 1:2) {
      Wd <- c(0, 0)
      Wd[h] <- 1
      ph <- four_stage_points(fr, h, Wd, r, fpc, dom$detail$.lower[h])
      ok <- ph$var <= (cv * fr$mean[h])^2 * (1 + 1e-12)
      expect_equal(dom$detail$n_1_int[h] *
                     (fr$cost_1[h] + fr$cost_2[h] * 3 +
                        fr$cost_3[h] * 3 * dom$detail$take_3_int[h] +
                        fr$cost_4[h] * 3 * dom$detail$take_3_int[h] *
                          dom$detail$take_4_int[h]),
                   min(ph$cost[ok]), tolerance = 1e-9, info = info)
    }
  }
})

test_that("fixed n under proportional allocation keeps its total at three stages", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(20000, 30000, 12000),
                      sd = c(10, 14, 8), mean = c(40, 50, 45),
                      icc_1 = c(0.05, 0.1, 0.08), icc_2 = c(0.15, 0.2, 0.1),
                      cost_1 = c(400, 600, 300), cost_2 = c(50, 40, 60),
                      cost_3 = c(10, 12, 8))
  fit <- n_alloc(frame, n = 600, alloc = "proportional")
  expect_lte(abs(sum(fit$detail$n_int) - 600),
             max(fit$detail$take_2_int * fit$detail$take_3_int) / 2)
  expect_null(fit$operational$status)
})

test_that("a three-stage design with too many states keeps the greedy one, labelled", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(20000, 30000, 12000),
                      sd = c(10, 14, 8), mean = c(40, 50, 45),
                      icc_1 = c(0.05, 0.1, 0.08), icc_2 = c(0.15, 0.2, 0.1),
                      cost_1 = c(400, 600, 300), cost_2 = c(50, 40, 60),
                      cost_3 = c(10, 12, 8))
  local_mocked_bindings(
    .alloc_whole_cluster_multistage = function(...) {
      list(status = "heuristic", reason = "listing")
    }
  )
  expect_warning(
    fit <- n_alloc(frame, cv = 0.02, alloc = "optimal"),
    "was not searched"
  )
  expect_identical(fit$operational$status, "heuristic")
  expect_lte(fit$operational$cv, 0.02 + 1e-8)
  expect_output(print(fit), "search: not run \\(too many take vectors to list\\)")
})

## Whole designs under binding size floors and stage caps

## Oracle (rung 2): every PSU count and take within the stage caps, priced
## from the definition of the multistage variance and cost.
floor_points <- function(fr, h, S, r, fpc, lower, W) {
  Ns <- c(vapply(seq_len(S - 1L), function(s) fr[[paste0("N_", s)]][h], 1),
          fr$N[h])
  caps <- c(Ns[1L], Ns[-1L] / Ns[-S])
  ranges <- lapply(seq_len(S), function(s) {
    fixed <- if (s > 1L) fr[[paste0("take_", s)]] else NULL
    if (!is.null(fixed) && !is.na(fixed[h])) return(fixed[h])
    from <- if (s == S) ceiling(1 / r - 1e-9) else 1
    if (from > floor(caps[s] + 1e-9)) integer(0) else
      seq.int(from, floor(caps[s] + 1e-9))
  })
  z <- as.matrix(expand.grid(ranges))
  p <- matrix(t(apply(z, 1L, cumprod)), ncol = S)
  ok <- p[, S] >= lower & p[, S] <= fr$N[h]
  z <- z[ok, , drop = FALSE]
  p <- p[ok, , drop = FALSE]
  if (!nrow(z)) return(data.frame(cost = numeric(0), var = numeric(0)))
  icc <- vapply(seq_len(S - 1L), function(s) fr[[paste0("icc_", s)]][h], 1)
  R <- c(icc, 1) * cumprod(c(1, 1 - icc))
  eff <- p
  eff[, S] <- eff[, S] * r
  term <- sweep(1 / eff, 2L, R, "*")
  if (fpc == "stage") {
    rr <- matrix(rep(c(rep(1, S - 1L), r), each = nrow(z)), nrow(z))
    term <- term * (1 - sweep(z * rr, 2L, caps, "/"))
  }
  v <- rowSums(term) * W[h]^2 * fr$sd[h]^2
  if (fpc == "unit") v <- v * (1 - eff[, S] / fr$N[h])
  cost <- drop(p %*% vapply(seq_len(S), function(s)
    fr[[paste0("cost_", s)]][h], 1))
  data.frame(cost = cost, var = v)
}

floor_best <- function(pts, budget = NULL, v_max = NULL) {
  cost <- pts[[1L]]$cost
  var <- pts[[1L]]$var
  if (length(pts) == 2L) {
    cost <- outer(cost, pts[[2L]]$cost, "+")
    var <- outer(var, pts[[2L]]$var, "+")
  }
  if (!is.null(budget)) {
    min(var[cost <= budget * (1 + 1e-9)])
  } else {
    min(cost[var <= v_max * (1 + 1e-12)])
  }
}

test_that("a size floor does not refuse a budget that other takes fund", {
  f <- data.frame(stratum = "A", N = 120, N_1 = 6, N_2 = 24, mean = 10,
                  sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 100, cost_2 = 1,
                  cost_3 = 1)
  for (budget in c(120, 115)) {
    expect_warning(
      fit <- n_alloc(f, budget = budget, min_n_stratum = 10,
                     alloc = "optimal", fpc = "none"),
      "holds a single PSU"
    )
    pts <- list(floor_points(f, 1, 3, 1, "none", 10, 1))
    expect_equal((fit$operational$cv * 10)^2, floor_best(pts, budget),
                 tolerance = 1e-9, info = budget)
    expect_identical(fit$operational$status, "optimal")
    expect_lte(fit$operational$cost, budget)
    expect_gte(fit$operational$n, 10)
  }
})

test_that("a floored four-stage optimum is not certified below the enumeration", {
  f4 <- data.frame(
    stratum = c("A", "B"), N = c(180, 96), N_1 = c(5, 6),
    N_2 = c(15, 12), N_3 = c(60, 24), mean = c(17, 15), sd = c(4, 9),
    icc_1 = c(.243364150484558, .302603007026482),
    icc_2 = c(.289602816966362, .110131041507702),
    icc_3 = c(.269387253315654, .091465327376500),
    cost_1 = c(51, 54), cost_2 = c(17, 1), cost_3 = c(13, 8),
    cost_4 = c(8, 8)
  )
  W <- f4$N / sum(f4$N)
  pts <- lapply(1:2, function(h)
    floor_points(f4, h, 4, .8, "none", c(10, 15)[h], W))
  best <- floor_best(pts, budget = 1186)
  expect_equal(best, 1.886599, tolerance = 1e-6)
  z <- suppressWarnings(n_alloc(f4, budget = 1186, min_n_stratum = c(10, 15),
                                alloc = "optimal", fpc = "none",
                                resp_rate = .8))
  expect_equal(z$operational$se^2, best, tolerance = 1e-9)
  expect_identical(z$optimization$search$status, "optimal")
  expect_lte(z$operational$cost, 1186)
})

test_that("a target at the variance floor is met, not refused", {
  frame <- data.frame(
    stratum = "A", N = 16, N_1 = 4, N_2 = 8, mean = 1, sd = 1,
    icc_1 = .25, icc_2 = 2 / 3, take_2 = 2, take_3 = 2,
    cost_1 = 5, cost_2 = 2, cost_3 = 1
  )
  expect_no_warning(
    fit <- n_alloc(frame, cv = .375, alloc = "optimal", fpc = "none")
  )
  expect_identical(fit$operational$status, "optimal")
  expect_equal(fit$operational$cost, 52)
  expect_equal(fit$detail$n_1_int, 4)
  expect_lte(fit$operational$cv, .375 * (1 + 1e-9))
  expect_error(
    n_alloc(frame, cv = .375 * (1 - 1e-6), alloc = "optimal", fpc = "none"),
    "unattainable|attains the target"
  )
})

floor_fixture <- function(i) {
  S <- if (i %% 2) 3 else 4
  H <- if (i %% 4 < 2) 1 else 2
  fpc <- c("none", "unit", "stage")[(i - 1) %% 3 + 1]
  r <- sample(c(.6, .8, 1), 1)
  cap <- cbind(sample(3:6, H, TRUE),
               matrix(sample(2:4, H * (S - 2), TRUE), H),
               sample(2:5, H, TRUE))
  pop <- t(apply(cap, 1, cumprod))
  fr <- data.frame(stratum = LETTERS[1:H], N = pop[, S],
                   mean = sample(10:25, H, TRUE), sd = sample(2:10, H, TRUE))
  costs <- c(sample(20:100, H, TRUE), sample(1:18, H * (S - 1), TRUE))
  for (j in 1:S) {
    fr[[paste0("cost_", j)]] <- costs[(j - 1) * H + 1:H]
    if (j < S) {
      fr[[paste0("N_", j)]] <- pop[, j]
      fr[[paste0("icc_", j)]] <- runif(H, .05, .5)
    }
  }
  if (i %% 5 == 0) fr$take_2 <- 2
  list(fr = fr, S = S, H = H, fpc = fpc, r = r,
       floor = sample(8:15, H, TRUE))
}

test_that("floored whole designs match enumeration at three and four stages", {
  set.seed(20261122)
  checked <- character(0)
  for (i in 1:40) {
    x <- floor_fixture(i)
    fr <- x$fr
    W <- fr$N / sum(fr$N)
    mu <- sum(W * fr$mean)
    pts <- lapply(seq_len(x$H), function(h)
      floor_points(fr, h, x$S, x$r, x$fpc, x$floor[h], W))
    if (any(vapply(pts, nrow, 1L) == 0L)) next
    cst <- pts[[1]]$cost
    var <- pts[[1]]$var
    if (x$H == 2) {
      cst <- outer(cst, pts[[2]]$cost, "+")
      var <- outer(var, pts[[2]]$var, "+")
    }
    budget <- floor(min(cst) + .21 * (max(cst) - min(cst)))
    v_t <- min(var) + .33 * (max(var) - min(var))
    for (mode in c("budget", "cv")) {
      id <- paste(i, mode)
      args <- c(list(frame = fr, alloc = "optimal", fpc = x$fpc,
                     resp_rate = x$r, min_n_stratum = x$floor),
                if (mode == "budget") list(budget = budget) else
                  list(cv = sqrt(v_t) / mu))
      fit <- tryCatch(suppressWarnings(do.call(n_alloc, args)),
                      error = function(e) e)
      expect_false(inherits(fit, "error"), info = id)
      if (inherits(fit, "error")) next
      if (mode == "budget") {
        expect_equal((fit$operational$cv * mu)^2,
                     floor_best(pts, budget = budget), tolerance = 1e-9,
                     info = id)
        expect_lte(fit$operational$cost, budget)
        expect_lte(fit$params$achieved$cost, budget * (1 + 1e-8))
      } else {
        expect_equal(fit$operational$cost, floor_best(pts, v_max = v_t),
                     tolerance = 1e-9, info = id)
      }
      expect_identical(fit$operational$status, "optimal", info = id)
      expect_true(all(fit$detail$n_int >= x$floor), info = id)
      checked <- c(checked, id)
    }
  }
  expect_length(checked, 78L)
})

test_that("a budget below the cheapest whole design names its cost", {
  f <- data.frame(stratum = "A", N = 120, N_1 = 6, N_2 = 24, mean = 10,
                  sd = 2, icc_1 = .1, icc_2 = .1, cost_1 = 100, cost_2 = 1,
                  cost_3 = 1)
  run <- function(budget) {
    suppressWarnings(n_alloc(f, budget = budget, min_n_stratum = 10,
                             alloc = "optimal", fpc = "none"))
  }
  # Enumeration gives 1 x 2 x 5 at 100 + 2 + 10 = 112 as the cheapest
  # design with 10 units.
  pts <- list(floor_points(f, 1, 3, 1, "none", 10, 1))
  expect_identical(min(pts[[1]]$cost), 112)
  for (budget in c(101, 111)) {
    expect_error(run(budget), "(minimum whole-design cost = 112)",
                 fixed = TRUE)
  }
  expect_identical(run(112)$operational$cost, 112)
})

test_that("every refused budget is one that enumeration cannot fund", {
  set.seed(20261123)
  checked <- integer(0)
  for (i in 1:40) {
    x <- floor_fixture(i)
    fr <- x$fr
    W <- fr$N / sum(fr$N)
    mu <- sum(W * fr$mean)
    pts <- lapply(seq_len(x$H), function(h)
      floor_points(fr, h, x$S, x$r, x$fpc, x$floor[h], W))
    if (any(vapply(pts, nrow, 1L) == 0L)) next
    cst <- pts[[1]]$cost
    if (x$H == 2) cst <- outer(cst, pts[[2]]$cost, "+")
    cheapest <- min(cst)
    run <- function(budget) {
      suppressWarnings(n_alloc(fr, budget = budget, alloc = "optimal",
                               fpc = x$fpc, resp_rate = x$r,
                               min_n_stratum = x$floor))
    }
    id <- paste("fixture", i)
    err <- tryCatch(run(cheapest * (1 - 1e-7)), error = function(e) e)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)
    # The integer pre-check refuses first when its weaker bound already
    # exceeds the budget. Otherwise the search states the cheapest design.
    if (grepl("whole-design cost", msg, fixed = TRUE)) {
      stated <- as.numeric(sub(".*minimum whole-design cost = ([0-9.]+)\\).*",
                               "\\1", msg))
      expect_equal(stated, cheapest, tolerance = 1e-7, info = id)
      checked <- c(checked, i)
    } else {
      expect_match(msg, "cannot fund the integer lower bounds", fixed = TRUE,
                   info = id)
    }
    # A budget at the cheapest design, or inside the comparison margin,
    # funds it.
    for (budget in c(cheapest * (1 - 5e-10), cheapest)) {
      fit <- run(budget)
      expect_equal((fit$operational$cv * mu)^2,
                   floor_best(pts, budget = budget), tolerance = 1e-9,
                   info = id)
    }
  }
  expect_identical(checked, c(1L, 11L, 13L, 14L, 15L, 16L, 19L, 21L, 23L, 33L))
})

test_that("target pruning keeps a feasible design at the decimal cost cap", {
  frame <- data.frame(stratum = c("A", "B"), N = 48, N_1 = 3, N_2 = 12,
                      mean = 10, sd = 2, icc_1 = c(.1, .2), icc_2 = c(.15, .2),
                      cost_1 = c(999, .1), cost_2 = c(.15, .02),
                      cost_3 = c(.01, .004))
  cv <- 0.03925
  points <- lapply(1:2, function(h) {
    floor_points(frame, h, 3, 1, "none", 1, c(0.5, 0.5))
  })
  best <- floor_best(points, v_max = (10 * cv)^2)
  cost_cols <- c("cost_1", "cost_2", "cost_3")
  for (scale in c(1, 1000)) {
    scaled <- frame
    scaled[cost_cols] <- frame[cost_cols] * scale
    fit <- n_alloc(scaled, cv = cv, alloc = "optimal", fpc = "none")
    expect_identical(fit$operational$status, "optimal")
    expect_equal(fit$operational$cost / scale, best, tolerance = 1e-12)
    expect_lte(fit$operational$cv, cv)
  }
})

test_that("a Neyman search counts final-stage units regardless of stage costs", {
  frame <- data.frame(stratum = c("A", "B"), N = c(120, 180),
                      N_1 = c(3, 5), N_2 = c(12, 15),
                      mean = c(10, 12), sd = c(2, 4), icc_1 = c(.1, .2),
                      icc_2 = c(.15, .2), take_2 = c(2, 3), take_3 = c(3, 4))
  for (scale in c(0, 1e-5, 1)) {
    priced <- frame
    if (scale > 0) {
      priced$cost_1 <- 1000 * scale
      priced$cost_2 <- 100 * scale
      priced$cost_3 <- 10 * scale
    }
    fit <- without_singleton_psu_warning(
      n_alloc(priced, cv = 0.08, alloc = "neyman", fpc = "none")
    )
    expect_identical(fit$operational$status, "optimal")
    expect_equal(fit$operational$n, 42)
    expect_lte(fit$operational$cv, 0.08)
  }
})

test_that("a Neyman incumbent holds the rounded takes used in its search", {
  frame <- data.frame(stratum = c("A", "B"), N = 48, N_1 = 3, N_2 = 12,
                      mean = 10, sd = 2, icc_1 = c(.204, .034), icc_2 = .035,
                      cost_1 = c(1.95, 1.38), cost_2 = c(3.93, 1.27),
                      cost_3 = c(1.69, .28))
  expect_warning(
    fit <- n_alloc(frame, cv = 0.0566, alloc = "neyman", fpc = "none"),
    "cost-optimal take_2 is below 1 expected response"
  )
  expect_identical(fit$operational$status, "optimal")
  expect_equal(fit$detail$take_2_int, c(1, 2))
  expect_equal(fit$detail$take_3_int, c(3, 4))
  expect_equal(fit$operational$n, 25)
  expect_lte(fit$operational$cv, 0.0566)
})

test_that("an empty frontier preserves a validated design with an incomplete search", {
  expect_identical(.ms_frontier(list(P = 1), 2, 1), integer(0))
  frame <- data.frame(stratum = c("A", "B"), N = 48, N_1 = 3, N_2 = 12,
                      mean = 10, sd = 2, icc_1 = c(.1, .2), icc_2 = c(.15, .2),
                      take_2 = 4, take_3 = 4)
  local_mocked_bindings(.ms_frontier = function(...) integer(0))
  expect_warning(
    fit <- n_alloc(frame, cv = 0.03925, alloc = "neyman", fpc = "none"),
    "numerical comparisons could not preserve the feasible starting design"
  )
  expect_identical(fit$operational$status, "heuristic")
  expect_equal(fit$operational$n, 96)
  expect_output(print(fit), "search: incomplete")

  local_mocked_bindings(.ms_incumbent = function(...) NULL)
  expect_error(n_alloc(frame, cv = 0.03925, alloc = "neyman", fpc = "none"),
               "no whole-cluster design found")
})

## The fallback fields the best validated starting design, which can
## differ from the greedy rounding when the takes vary. Mocking the listing
## away gives the greedy rounding by an independent route.
fallback_frame <- function() {
  data.frame(stratum = c("A", "B"), N = 48, N_1 = 3, N_2 = 12, mean = 10,
             sd = 2, icc_1 = c(.15, .21), icc_2 = c(.18, .15),
             cost_1 = c(46, 49), cost_2 = c(4, 10), cost_3 = 1)
}

fallback_design <- function(fit) {
  d <- fit$detail
  cbind(d$n_1_int, d$take_2_int, d$take_3_int)
}

expect_fallback_design <- function(fit, frame, full) {
  d <- fit$detail
  expect_identical(fit$operational$status, "heuristic")
  expect_equal(d$n_int, d$n_1_int * d$take_2_int * d$take_3_int)
  expect_true(all(d$n_1_int >= 1 & d$n_1_int <= frame$N_1))
  expect_true(all(d$n_1_int * d$take_2_int <= frame$N_2))
  expect_equal(fallback_design(fit), fallback_design(full))
}

test_that("a pruned search fields the cheapest validated design, not the rounding", {
  frame <- fallback_frame()
  run <- function() n_alloc(frame, cv = 0.052, alloc = "optimal", fpc = "none")
  full <- run()
  greedy <- suppressWarnings(
    with_mocked_bindings(run(), .ms_enumerate = function(...) NULL)
  )
  expect_identical(greedy$operational$status, "heuristic")
  expect_warning(
    fit <- with_mocked_bindings(run(), .ms_frontier = function(...) integer(0)),
    "numerical comparisons could not preserve the feasible starting design"
  )
  expect_fallback_design(fit, frame, full)
  expect_lte(fit$operational$cv, 0.052)
  expect_lt(fit$operational$cost, greedy$operational$cost)
  expect_lte(greedy$operational$cv, 0.052)
})

test_that("a pruned budget search fields the most precise validated design", {
  frame <- fallback_frame()
  run <- function() n_alloc(frame, budget = 357, alloc = "optimal", fpc = "none")
  full <- run()
  greedy <- suppressWarnings(
    with_mocked_bindings(run(), .ms_enumerate = function(...) NULL)
  )
  expect_warning(
    fit <- with_mocked_bindings(run(), .ms_frontier = function(...) integer(0)),
    "numerical comparisons could not preserve the feasible starting design"
  )
  expect_fallback_design(fit, frame, full)
  expect_lte(fit$operational$cost, 357)
  expect_lt(greedy$operational$cost, fit$operational$cost)
  expect_lt(fit$operational$cv, greedy$operational$cv)
})

test_that("every lost search state falls back to the validated design", {
  frame <- fallback_frame()
  run <- function() n_alloc(frame, cv = 0.052, alloc = "optimal", fpc = "none")
  full <- run()
  no_takes <- function(prep, h, b_rule, ...) {
    matrix(numeric(0), 0L, length(b_rule))
  }
  for (lost in c("listing", "infeasible", "limit")) {
    expect_warning(
      fit <- if (lost == "listing") {
        with_mocked_bindings(run(), .ms_enumerate = no_takes)
      } else {
        with_mocked_bindings(
          run(), .ms_search = function(...) list(status = lost, a = NULL)
        )
      },
      "numerical comparisons could not preserve the feasible starting design",
      info = lost
    )
    expect_fallback_design(fit, frame, full)
    expect_lte(fit$operational$cv, 0.052)
  }
})

test_that("a lost design is never refused as unaffordable", {
  # With no validated starting design, losing every state used to read as a
  # proof that the budget funds nothing. The cheapest design is found apart
  # from the pruning, so the fit falls back to the rounding instead.
  frame <- fallback_frame()
  run <- function() n_alloc(frame, budget = 357, alloc = "optimal", fpc = "none")
  lost <- list(
    frontier = function() {
      with_mocked_bindings(run(), .ms_frontier = function(...) integer(0),
                           .ms_incumbent = function(...) NULL)
    },
    search = function() {
      with_mocked_bindings(
        run(), .ms_search = function(...) list(status = "infeasible", a = NULL),
        .ms_incumbent = function(...) NULL
      )
    }
  )
  for (route in names(lost)) {
    expect_warning(
      fit <- lost[[route]](),
      "numerical comparisons could not preserve the feasible starting design",
      info = route
    )
    expect_identical(fit$operational$status, "heuristic", info = route)
    expect_lte(fit$operational$cost, 357)
  }
})

test_that("a starting design outside the PSU bounds is never fielded", {
  frame <- data.frame(stratum = c("A", "B"), N = 48, N_1 = 3, N_2 = 12,
                      mean = 10, sd = 2, icc_1 = c(.1, .2), icc_2 = c(.15, .2),
                      take_2 = 4, take_3 = 4)
  start <- .ms_incumbent
  local_mocked_bindings(
    .ms_frontier = function(...) integer(0),
    .ms_incumbent = function(prep, alloc, b_mat, a_h, ...) {
      start(prep, alloc, b_mat, a_h + 100, ...)
    }
  )
  expect_error(n_alloc(frame, cv = 0.03925, alloc = "neyman", fpc = "none"),
               "no whole-cluster design found")
})

test_that("a fallback warning names why the design was not searched", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(20000, 30000, 12000),
                      sd = c(10, 14, 8), mean = c(40, 50, 45),
                      icc_1 = c(0.05, 0.1, 0.08), icc_2 = c(0.15, 0.2, 0.1),
                      cost_1 = c(400, 600, 300), cost_2 = c(50, 40, 60),
                      cost_3 = c(10, 12, 8))
  local_mocked_bindings(.ms_frontier = function(...) NULL)
  expect_warning(
    fit <- n_alloc(frame, cv = 0.02, alloc = "optimal"),
    "frontier has too many states"
  )
  expect_identical(fit$operational$status, "heuristic")
  expect_lte(fit$operational$cv, 0.02 + 1e-8)
})

test_that("a numerical fallback at the variance floor is named", {
  frame <- data.frame(stratum = c("A", "B", "C"), N = c(20000, 30000, 12000),
                      sd = c(10, 14, 8), mean = c(40, 50, 45),
                      icc_1 = c(0.05, 0.1, 0.08), icc_2 = c(0.15, 0.2, 0.1),
                      cost_1 = c(400, 600, 300), cost_2 = c(50, 40, 60),
                      cost_3 = c(10, 12, 8))
  local_mocked_bindings(.alloc_whole_cluster_multistage = function(...) {
    list(status = "heuristic", reason = "numerical")
  })
  expect_warning(
    fit <- n_alloc(frame, cv = 0.02, alloc = "optimal"),
    "the target sits at the variance floor"
  )
  expect_output(print(fit), "too many designs at the variance floor")
})

## A target at the variance floor

test_that("the floor design meets the target as one sum, not stratum by stratum", {
  prep <- list(cost_mat = matrix(c(1, 10), nrow = 1))
  B2 <- matrix(c(1, 1, 2, 3), nrow = 2)
  ## Stratum 1's cheap branch sits .5e-12 above its least variance, inside
  ## the target's 1e-12 allowance once stratum 2 sits at its least.
  states <- list(
    list(P = c(2, 3), pc = c(10, 1), C = c(.5 - .9e-12, .5 + .5e-12),
         D = c(0, 0), lo = c(1, 1), hi = c(1, 1), B = B2),
    list(P = 2, pc = 1, C = .5, D = 0, lo = 1, hi = 1,
         B = matrix(c(1, 2), nrow = 1))
  )
  z <- .ms_search(states, "cv", 1, NULL, prep, TRUE, 100, 0, NULL)
  expect_identical(z$status, "optimal")
  expect_equal(z$search$bound, 2)
  expect_equal(z$b, rbind(c(1, 3), c(1, 2)))
  ## Two cheap branches .8e-12 above their least each fit one by one but
  ## not together, so one stratum takes its dear branch.
  twin <- list(P = c(2, 3), pc = c(10, 1), C = c(.5, .5 + .8e-12),
               D = c(0, 0), lo = c(1, 1), hi = c(1, 1), B = B2)
  z <- .ms_search(list(twin, twin), "cv", 1, NULL, prep, TRUE, 100, 0, NULL)
  expect_equal(z$search$bound, 11)
  expect_lte(sum(z$b[, 2] == 3), 1)
  z <- .ms_search(list(twin, twin), "cv", sqrt(1 - 1e-11), NULL, prep, TRUE,
                  100, 0, NULL)
  expect_identical(z$status, "infeasible")
})

## A continuous precheck at the cost-optimal takes

test_that("a floored census design is found in every target mode", {
  a5 <- data.frame(stratum = "A", N = 18, N_1 = 3, N_2 = 9, mean = 10, sd = 2,
                   icc_1 = .087580317875836, icc_2 = .337548924167640,
                   cost_1 = 83, cost_2 = 3, cost_3 = 9)
  fit <- function(...) {
    suppressWarnings(n_alloc(a5, alloc = "optimal", fpc = "none",
                             resp_rate = .6, min_n_stratum = 13, ...))
  }
  ## Takes are capped at 3 and 2, so 13 units need the census of 3 PSUs of
  ## 3 SSUs of 2, at 3 * 83 + 9 * 3 + 18 * 9 = 438.
  cv_census <- fit(budget = 438)$cv
  for (target in list(list(budget = 438), list(cv = cv_census),
                      list(moe = cv_census * 10 * qnorm(.975)),
                      list(rmoe = cv_census * qnorm(.975)))) {
    z <- do.call(fit, target)
    info <- names(target)
    expect_equal(c(z$detail$n_1_int, z$detail$take_2_int,
                   z$detail$take_3_int), c(3, 3, 2), info = info)
    expect_equal(z$operational$cost, 438, info = info)
    expect_equal(z$params$achieved$cost, 438, info = info)
    expect_identical(z$operational$status, "optimal", info = info)
    expect_null(z$params$frame$take_2, info = info)
    expect_equal(suppressWarnings(prec_alloc(z))$cv, z$cv, info = info)
  }
  z <- fit(budget = 438)
  p <- suppressWarnings(predict(z, data.frame(budget = 438)))
  expect_equal(p$cv, z$cv)
  expect_equal(p$cost, 438)
  expect_error(fit(budget = 437), "budget")
})

test_that("a takes-first search that stops with no design says so", {
  f <- data.frame(stratum = c("A", "B"), N = c(48, 60), N_1 = 3, N_2 = 12,
                  mean = c(12, 14), sd = c(2, 9),
                  icc_1 = c(.317279242968652, .0821973482845351),
                  icc_2 = c(.174681426165625, .395742602960672),
                  cost_1 = c(100, 37), cost_2 = c(6, 5), cost_3 = c(6, 10))
  local_mocked_bindings(.alloc_whole_cluster_multistage = function(...) {
    list(status = "limit", a = NULL)
  })
  expect_error(
    suppressWarnings(n_alloc(f, budget = 610, alloc = "optimal", fpc = "none",
                             min_n_stratum = c(12, 14))),
    "at the continuous takes", class = "svyplan_plan_unknown"
  )
})

test_that("a takes-first fit resizes over free takes and fields its first search", {
  f <- data.frame(stratum = c("A", "B"), N = c(48, 60), N_1 = 3, N_2 = 12,
                  mean = c(12, 14), sd = c(2, 9),
                  icc_1 = c(.317279242968652, .0821973482845351),
                  icc_2 = c(.174681426165625, .395742602960672),
                  cost_1 = c(100, 37), cost_2 = c(6, 5), cost_3 = c(6, 10))
  fit <- function(...) {
    suppressWarnings(n_alloc(f, alloc = "optimal", fpc = "none",
                             min_n_stratum = c(12, 14), ...))
  }
  z <- fit(budget = 610)
  expect_false(is.null(z$params$searched_takes))
  pr <- suppressWarnings(prec_alloc(z))
  expect_equal(pr$cv, z$cv)
  expect_null(pr$params$frame$take_2)
  scenario <- predict(pr, data.frame(resp_rate = c(1, .8)))
  expect_equal(scenario$n, rep(pr$params$achieved$n, 2))
  expect_equal(scenario$cv[1], pr$cv)
  expect_equal(scenario$cost[1], pr$params$achieved$cost)
  held <- transform(f, take_2 = pr$detail$take_2, take_3 = pr$detail$take_3)
  reduced <- prec_alloc(held, n = pr$params$n, fpc = "none", resp_rate = .8)
  expect_equal(scenario$cv[2], reduced$cv)
  expect_equal(scenario$cost[2], reduced$params$achieved$cost)
  ## A changed target sizes the original request, not the chosen takes.
  expect_equal(suppressWarnings(n_alloc(pr, cv = .13))$operational$cost,
               fit(cv = .13)$operational$cost)
  expect_equal(suppressWarnings(n_alloc(pr, budget = 600))$operational$cv,
               fit(budget = 600)$operational$cv)
  same <- suppressWarnings(n_alloc(pr))
  expect_equal(same$operational$cost, z$operational$cost)
  expect_equal(same$cv, z$cv)
  expect_error(suppressWarnings(prec_alloc(z, take_2 = 2)),
               "needs more PSUs than 'N_1' holds at its takes for stratum: B")
  over <- suppressWarnings(prec_alloc(z, take_2 = c(2, 4)))
  expect_equal(over$detail$take_2, c(2, 4))
  expect_equal(over$params$frame$take_2, c(2, 4))
  expect_identical(colnames(over$params$searched_takes), "take_3")
  expect_equal(predict(over, data.frame(resp_rate = 1))$cv, over$cv)

  searches <- 0L
  real <- .alloc_whole_cluster_multistage
  local_mocked_bindings(.alloc_whole_cluster_multistage = function(...) {
    searches <<- searches + 1L
    real(...)
  })
  fit(budget = 610)
  expect_identical(searches, 1L)
})
