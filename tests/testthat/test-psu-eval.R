## The state evaluator against the shipped assessment of a register design

# A stratum's state at f under the R closure, with the takes and draw the
# shipped assessment reads and the sums .psu_eval_states() reads.
eval_state_at <- function(f, sizes, N, take, m = 0L) {
  n_h <- f * N
  held <- .psu_classify(sizes, take / f)$certain
  closed <- .psu_draw_closure(n_h, N, take, sizes, held, 1, m)
  counts <- .psu_stratum_counts(n_h, N, take, sizes, closed$certain, m)
  M <- sizes[closed$certain]
  q <- counts$per_psu
  list(
    certain = closed$certain, per_psu = q, draw = counts$n_psu_draw,
    sums = data.frame(Q = sum(q), A = sum(M^2 / q), held = sum(M),
                      n_c = length(M), a = counts$n_psu_draw)
  )
}

# The assessment .n_alloc_psu() applies to a whole-unit design: each row's
# design effect from .psu_field_deff(), then the element assessment.
eval_shipped <- function(z, states, unit_cost = NULL, alpha = 0.05,
                         deff = 1, resp_rate = 1, objective = NULL,
                         df = NULL) {
  stratum <- z$frame$stratum
  idx_of <- split(seq_len(nrow(z$psu)), factor(z$psu$stratum, stratum))
  take <- rep_len(z$frame$n_per_psu, length(stratum))
  icc_row <- .psu_icc(z$measures, stratum, idx_of)
  resp_row <- .joint_row_value("resp_rate", z$measures, z$frame, resp_rate)
  user_deff <- .joint_row_value("deff", z$measures, z$frame, deff)
  m <- z$measures[, setdiff(names(z$measures), "icc_psu"), drop = FALSE]
  m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
    h <- match(m$stratum[r], stratum)
    s <- states[[h]]
    .psu_field_deff(z$psu$N[idx_of[[h]]], s$certain, s$per_psu, s$draw,
                    take[h], icc_row[r], resp_row[r])
  }, numeric(1))
  n <- vapply(seq_along(states), function(h) {
    sum(states[[h]]$per_psu) + states[[h]]$draw * take[h]
  }, numeric(1))
  base_frame <- z$frame[, setdiff(names(z$frame),
                                  c("n_per_psu", "cost_psu", "cost_ssu"))]
  .prec_alloc_bethel(
    frame = base_frame, n = n, measures = m, targets = z$targets,
    unit_cost = unit_cost, alpha = alpha, deff = deff, resp_rate = resp_rate,
    objective = objective, df = df, .allow_fractional_stages = TRUE,
    .finite_prop_var = FALSE
  )
}

eval_sum <- function(ev, sums) {
  out <- lapply(seq_along(sums), function(h) .psu_eval_states(ev, h, sums[[h]]))
  list(
    ratio2 = Reduce(`+`, lapply(out, function(e) e$C[1, ])),
    cost = sum(vapply(out, `[[`, numeric(1), "cost")),
    objective = sum(vapply(out, function(e) e$objective %||% 0, numeric(1)))
  )
}

## Fitted plans

test_that("the evaluator gives a fitted register plan's cost and precision", {
  plan_sums <- function(fit, ev) {
    lapply(seq_along(ev$idx_of), function(h) {
      i <- ev$idx_of[[h]]
      held <- fit$psu$certainty[i]
      q <- fit$psu$n_take[i][held]
      M <- fit$psu$N[i][held]
      data.frame(Q = sum(q), A = sum(M^2 / q), held = sum(M),
                 n_c = sum(held), a = fit$detail$n_psu_draw[h])
    })
  }
  base <- .psu_fixture(4)
  staged <- base
  staged$frame$cost_psu <- c(100, 120, 90, 150)
  staged$frame$cost_ssu <- 10
  flagged <- base
  flagged$psu$certainty <- seq_len(nrow(base$psu)) %in% c(7, 305, 520)
  mixed <- base
  mixed$targets <- data.frame(name = c("literacy", "income"),
                              moe = c(0.02, NA), rmoe = c(NA, 0.09))
  cases <- list(
    list(z = base),
    list(z = base, resp_rate = 0.7, deff = 1.3),
    list(z = base, unit_cost = c(1, 1.5, 2, 2.5), certainty_cutoff = 0.8),
    list(z = staged),
    list(z = staged, n_psu_per_zone = 1L),
    list(z = staged, n_psu_per_zone = 2L, resp_rate = 0.85),
    list(z = flagged),
    list(z = mixed, df = 30, alpha = 0.1),
    list(z = .psu_random_register(3)),
    list(z = .psu_wide_register(2))
  )
  for (cs in cases) {
    z <- cs$z
    opts <- cs[setdiff(names(cs), "z")]
    fit <- suppressWarnings(do.call(n_alloc, c(
      list(z$frame, measures = z$measures, targets = z$targets, psu = z$psu),
      opts
    )))
    ev <- .psu_setup(
      z$frame, .check_psu_table(z$psu, z$frame, z$measures), z$measures,
      z$targets, unit_cost = opts$unit_cost, alpha = opts$alpha %||% 0.05,
      deff = opts$deff %||% 1, resp_rate = opts$resp_rate %||% 1,
      df = opts$df
    )
    got <- eval_sum(ev, plan_sums(fit, ev))
    expect_equal(got$ratio2, fit$operational$constraints$.ratio^2,
                 tolerance = 1e-12)
    expect_equal(got$cost, fit$operational$cost, tolerance = 1e-12)
  }
})

## Any state

test_that("the evaluator matches the shipped assessment on any state", {
  z <- .psu_fixture(7)
  z$frame$region <- c("north", "north", "south", "south")
  z$measures$resp_rate <- c(0.9, 0.8, NA, 0.7, NA, 0.75, 0.95, 0.6)
  z$measures$deff <- c(1.2, NA, 1.1, NA, 1.4, NA, NA, 1.05)
  z$targets <- data.frame(
    name = c("literacy", "income", "literacy", "income"),
    domain = c(".overall", ".overall", "region", "region"),
    level = c(NA, NA, "north", "south"),
    cv = c(0.04, NA, NA, NA), moe = c(NA, NA, 0.03, NA),
    rmoe = c(NA, 0.1, NA, 0.12), alpha = c(NA, 0.1, NA, NA)
  )
  objective <- data.frame(
    name = c("literacy", "income", "income"),
    domain = c("region", ".overall", "region"),
    level = c("south", NA, "north"), priority = c(2, 1, 0.5)
  )
  staged <- z
  staged$frame$cost_psu <- c(100, 120, 90, 150)
  staged$frame$cost_ssu <- c(10, 12, 9, 15)
  idx_of <- split(seq_len(nrow(z$psu)), factor(z$psu$stratum, z$frame$stratum))
  set.seed(11)
  for (draw in 1:40) {
    zz <- if (draw %% 2) z else staged
    unit_cost <- if (draw %% 2) c(1, 1.5, 2, 2.5)
    m <- c(0L, 1L, 2L)[draw %% 3 + 1]
    resp <- if (draw %% 4 < 2) 1 else 0.8
    # Fractions from a sparse remainder to a census, where the stratum-rate
    # variance of a row at full response reaches zero.
    f <- sample(c(0.002, 0.01, 0.04, 0.2, 0.6, 1), 4, replace = TRUE)
    states <- lapply(1:4, function(h) {
      eval_state_at(f[h], zz$psu$N[idx_of[[h]]], zz$frame$N[h],
                    zz$frame$n_per_psu[h], m)
    })
    want <- eval_shipped(zz, states, unit_cost = unit_cost, deff = 1.15,
                         resp_rate = resp, objective = objective, df = 25)
    ev <- .psu_setup(zz$frame, zz$psu, zz$measures, zz$targets,
                          unit_cost = unit_cost, deff = 1.15,
                          resp_rate = resp, objective = objective, df = 25)
    got <- eval_sum(ev, lapply(states, `[[`, "sums"))
    expect_equal(got$ratio2, want$detail$.ratio^2, tolerance = 1e-12)
    expect_equal(got$objective, want$objective_value, tolerance = 1e-12)
    n <- vapply(states, function(s) sum(s$per_psu) + s$draw * 12, numeric(1))
    cost <- if (draw %% 2) {
      sum(n * unit_cost)
    } else {
      sum(vapply(states, function(s) s$sums$n_c + s$draw, numeric(1)) *
            staged$frame$cost_psu + n * staged$frame$cost_ssu)
    }
    expect_equal(got$cost, cost, tolerance = 1e-12)
  }
})

test_that("a field variance below the stratum rate is raised to it", {
  # Only a take times response rate below 1 shrinks the remainder's term
  # below its stratum-rate value, so only there does the floor bind.
  z <- .psu_fixture(4, take = 1)
  idx_of <- split(seq_len(nrow(z$psu)), factor(z$psu$stratum, z$frame$stratum))
  ev <- .psu_setup(z$frame, z$psu, z$measures, z$targets,
                        resp_rate = 0.5)
  floored <- 0
  for (f in c(0.001, 0.004, 0.02)) {
    states <- lapply(1:4, function(h) {
      eval_state_at(f, z$psu$N[idx_of[[h]]], z$frame$N[h], 1)
    })
    want <- eval_shipped(z, states, resp_rate = 0.5)
    got <- eval_sum(ev, lapply(states, `[[`, "sums"))
    expect_equal(got$ratio2, want$detail$.ratio^2, tolerance = 1e-12)
    floored <- floored + sum(want$params$measures$deff == 1)
  }
  expect_gt(floored, 0)
})

test_that("the evaluator scores a listed menu state by state", {
  z <- .psu_fixture(4)
  ev <- .psu_setup(z$frame, z$psu, z$measures, z$targets)
  i <- ev$idx_of[[2]]
  sizes <- z$psu$N[i]
  part <- .psu_menu_part(sizes, z$frame$N[2], 12, 1e-6, 0.05)$states
  got <- .psu_eval_states(ev, 2, part)
  pick <- unique(round(seq(1, nrow(part), length.out = 25)))
  for (j in pick) {
    s <- eval_state_at(part$f[j], sizes, z$frame$N[2], 12)
    one <- .psu_eval_states(ev, 2, s$sums)
    expect_equal(got$C[j, ], one$C[1, ], tolerance = 1e-12)
    expect_equal(got$cost[j], one$cost)
  }
})
