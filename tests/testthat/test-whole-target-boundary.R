## A whole-unit search returns a design whose precision equals the target

# Oracle (rung 2): every whole design of a small problem is enumerated and
# the cheapest one meeting the target is read off. A target counts as met
# within a relative 1e-12 of the variance, the rounding the searches allow.
meets <- function(cv, target) cv^2 <= target^2 * (1 + 1e-12)

## Two-phase designs

twophase_frame <- data.frame(
  stratum = letters[1:3], N = 6,
  sd = c(0, 3, 1.5), mean = c(7, 12, 9),
  unit_cost = c(4, 3, 3), resp_rate = c(1, 1, 0.8)
)

# Every whole two-phase design up to a phase-1 size of n1_max: phase 1 of n1
# units, then n_h of the pool floor(resp_rate * W_h * n1) in each stratum,
# priced at phase1_cost per phase-1 unit and unit_cost per phase-2 unit.
twophase_designs <- function(f, n1_max, resp_rate, N, phase1_cost) {
  W <- f$N / sum(f$N)
  rows <- list()
  for (n1 in seq_len(n1_max)) {
    cap <- floor(resp_rate * W * n1)
    if (any(cap < 1)) next
    grid <- as.matrix(expand.grid(lapply(cap, seq_len)))
    cv <- apply(grid, 1L, function(n_h) {
      prec_twophase(transform(f, nu = n_h / (W * n1)), n_phase1 = n1, N = N,
                    resp_rate = resp_rate, phase1_cost = phase1_cost)$cv
    })
    rows[[length(rows) + 1L]] <- data.frame(
      n1 = n1, cost = phase1_cost * n1 + drop(grid %*% f$unit_cost), cv = cv
    )
  }
  do.call(rbind, rows)
}

test_that("a two-phase target read back from prec_twophase() is met at its cost", {
  known <- prec_twophase(transform(twophase_frame, nu = c(1, 3, 2) / 4),
                         n_phase1 = 12, N = 18, resp_rate = 0.75,
                         phase1_cost = 1)
  expect_equal(known$params$cost, 31)
  fit <- n_twophase(twophase_frame, phase1_cost = 1, N = 18, cv = known$cv,
                    resp_rate = 0.75)
  designs <- twophase_designs(twophase_frame, 18, 0.75, 18, 1)
  best <- min(designs$cost[meets(designs$cv, known$cv)])
  expect_equal(best, 31)
  expect_equal(fit$operational$cost, best)
  expect_identical(fit$operational$status, "optimal")
  expect_true(meets(fit$operational$cv, known$cv))
})

test_that("every attainable two-phase cv is met at the least cost", {
  designs <- twophase_designs(twophase_frame, 18, 0.75, 18, 1)
  picks <- c(5L, 40L, 120L, nrow(designs) - 3L)
  seen <- character(0)
  for (i in picks) {
    target <- designs$cv[i]
    fit <- n_twophase(twophase_frame, phase1_cost = 1, N = 18, cv = target,
                      resp_rate = 0.75)
    label <- sprintf("design %d", i)
    expect_equal(fit$operational$cost,
                 min(designs$cost[meets(designs$cv, target)]), label = label)
    expect_true(meets(fit$operational$cv, target), label = label)
    seen <- c(seen, label)
  }
  expect_length(unique(seen), length(picks))
})

## Classic allocations

alloc_frame <- data.frame(
  stratum = c("a", "b", "c"), N = c(24, 31, 18),
  sd = c(4.5, 9.1, 2.3), mean = c(31, 52, 44),
  unit_cost = c(2, 5, 1)
)

# The stratified SRSWOR variance of the mean, sum_h W_h^2 S_h^2 (1/n_h - 1/N_h).
alloc_cv <- function(f, n) {
  W <- f$N / sum(f$N)
  n <- matrix(n, ncol = nrow(f))
  fpc <- 1 / n - matrix(1 / f$N, nrow(n), nrow(f), byrow = TRUE)
  sqrt(drop(fpc %*% (W^2 * f$sd^2))) / sum(W * f$mean)
}

test_that("the allocation variance matches prec_alloc() on the frame", {
  rows <- rbind(c(3, 5, 2), c(10, 17, 4), c(24, 31, 18), c(1, 1, 1))
  expect_equal(alloc_cv(alloc_frame, rows),
               apply(rows, 1L, function(n) prec_alloc(alloc_frame, n = n)$cv),
               tolerance = 1e-12)
})

test_that("an allocation's own cv as the target returns the least design", {
  grid <- as.matrix(expand.grid(lapply(alloc_frame$N, seq_len)))
  cv <- alloc_cv(alloc_frame, grid)
  seen <- character(0)
  for (alloc in c("neyman", "optimal")) {
    weight <- if (alloc == "neyman") rep(1, 3) else alloc_frame$unit_cost
    objective <- drop(grid %*% weight)
    for (start in c(0.03, 0.06, 0.1)) {
      fit <- n_alloc(alloc_frame, cv = start, alloc = alloc)
      target <- fit$operational$cv
      again <- n_alloc(alloc_frame, cv = target, alloc = alloc)
      label <- sprintf("%s at %.2f", alloc, start)
      best <- min(objective[meets(cv, target)])
      expect_equal(sum(again$detail$n_int * weight), best, label = label)
      expect_equal(sum(fit$detail$n_int * weight), best, label = label)
      seen <- c(seen, label)
    }
  }
  expect_length(unique(seen), 6L)
})

test_that("a two-stage allocation's own cv as the target costs no more", {
  fr <- data.frame(
    stratum = c("North", "South"), N = c(40000, 60000),
    sd = c(12, 15), mean = c(50, 60),
    icc_psu = c(0.05, 0.08), var_ratio_psu = 1,
    cost_psu = c(400, 550), cost_ssu = c(45, 60)
  )
  seen <- character(0)
  for (alloc in c("neyman", "optimal")) {
    for (start in c(0.01, 0.02)) {
      fit <- n_alloc(fr, cv = start, alloc = alloc)
      again <- n_alloc(fr, cv = fit$operational$cv, alloc = alloc)
      label <- sprintf("%s at %.2f", alloc, start)
      # The design that set the target meets it, so the least design costs
      # no more than it does.
      if (alloc == "neyman") {
        expect_lte(sum(again$detail$n_int), sum(fit$detail$n_int),
                   label = label)
      } else {
        expect_lte(again$operational$cost, fit$operational$cost,
                   label = label)
      }
      expect_true(meets(again$operational$cv, fit$operational$cv),
                  label = label)
      seen <- c(seen, label)
    }
  }
  expect_length(unique(seen), 4L)
})
