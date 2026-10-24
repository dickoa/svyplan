## Budget mode on the register planner against every plan of the menus

budget_register <- function(seed) {
  set.seed(seed)
  H <- if (seed %% 3 == 0) 3L else 2L
  take <- sample(c(5, 10), 1)
  psu <- do.call(rbind, lapply(seq_len(H), function(h) {
    k <- if (H == 3L) sample(2:3, 1) else sample(c(3, 5, 7), 1)
    sizes <- pmax(round(stats::rlnorm(k, log(if (H == 3L) 12 else 20), 0.9)),
                  take)
    if (k > 4) sizes[1:2] <- sizes[1:2] * 4
    data.frame(stratum = LETTERS[h], N = sizes)
  }))
  frame <- data.frame(
    stratum = LETTERS[seq_len(H)],
    N = as.numeric(tapply(psu$N, psu$stratum, sum)),
    n_per_psu = take
  )
  if (seed %% 2 == 0) {
    frame$cost_psu <- 50
    frame$cost_ssu <- 5
  }
  list(
    frame = frame, psu = psu,
    measures = data.frame(
      stratum = rep(frame$stratum, 2),
      name = rep(c("y", "x"), each = H),
      p = c(stats::runif(H, 0.1, 0.6), rep(NA, H)),
      mean = c(rep(NA, H), stats::runif(H, 50, 100)),
      sd = c(rep(NA, H), stats::runif(H, 20, 60)),
      icc_psu = stats::runif(2 * H, 0.02, 0.2)
    ),
    targets = data.frame(name = c("y", "x"), cv = c(0.12, 0.1)),
    objective = data.frame(name = c("y", "x"), priority = c(1, 3)),
    zone = c(0L, 1L, 2L)[seed %% 4 %% 3 + 1]
  )
}

budget_fit <- function(z, budget, hard) {
  n_alloc(z$frame, measures = z$measures, targets = if (hard) z$targets,
          objective = z$objective, budget = budget, psu = z$psu,
          n_psu_per_zone = if (z$zone > 0) z$zone)
}

# Every plan of the complete menus (test-psu-list.R holds them against brute
# force): each one's objective, field cost and target loads.
budget_plans <- function(z) {
  ev <- .psu_setup(
    z$frame, .check_psu_table(z$psu, z$frame, z$measures), z$measures,
    z$targets, objective = z$objective,
    n_psu_per_zone = if (z$zone > 0) z$zone
  )
  menu <- .psu_menu_round(.psu_menu_new(ev), ev,
                          list(price = 1, mu = c(0, 0), obj = 0),
                          rep(Inf, length(ev$N)))
  idx <- expand.grid(lapply(menu, function(m) seq_len(nrow(m$states))))
  add <- function(f) Reduce(`+`, Map(f, menu, idx))
  list(
    objective = add(function(m, i) m$states$objective[i]),
    cost = add(function(m, i) m$states$cost[i]),
    load = sapply(1:2, function(k) add(function(m, i) m$C[i, k]))
  )
}

best_within <- function(plans, budget, hard) {
  ok <- plans$cost <= budget * (1 + 1e-9)
  if (hard) ok <- ok & rowSums(plans$load <= 1 + 1e-9) == 2
  min(plans$objective[ok])
}

test_that("a budget plan is the best plan of the whole menus", .quiet_single_draw({
  for (seed in 1:12) {
    z <- budget_register(seed)
    plans <- budget_plans(z)
    cheapest <- min(plans$cost)
    meets <- rowSums(plans$load <= 1 + 1e-9) == 2
    cheapest_hard <- min(plans$cost[meets])
    cases <- list(
      list(budget = cheapest, hard = FALSE),
      list(budget = 1.5 * cheapest, hard = FALSE),
      list(budget = cheapest_hard, hard = TRUE),
      list(budget = 1.3 * cheapest_hard, hard = TRUE)
    )
    for (case in cases) {
      fit <- budget_fit(z, case$budget, case$hard)
      op <- fit$operational
      label <- sprintf("seed %d, budget %g, hard %s", seed, case$budget,
                       case$hard)
      expect_equal(op$objective_value,
                   best_within(plans, case$budget, case$hard),
                   tolerance = 1e-9, label = label)
      expect_identical(op$status, "optimal", label = label)
      expect_lte(op$cost, case$budget * (1 + 1e-9), label = label)
      expect_equal(op$budget_residual, case$budget - op$cost)
      expect_true(all(op$constraints$.pass), label = label)
      expect_equal(prec_alloc(fit)$objective_value, op$objective_value,
                   tolerance = 1e-12)
    }
  }
}))

test_that("a budget below the cheapest plan is refused with its cost", {
  for (seed in c(2, 3, 5)) {
    z <- budget_register(seed)
    plans <- budget_plans(z)
    cheapest <- min(plans$cost)
    expect_error(
      budget_fit(z, 0.99 * cheapest, FALSE),
      sprintf("'budget' is below the cheapest field design, which costs %.6g",
              cheapest),
      fixed = TRUE
    )
    meets <- rowSums(plans$load <= 1 + 1e-9) == 2
    cheapest_hard <- min(plans$cost[meets])
    expect_error(
      budget_fit(z, 0.99 * cheapest_hard, TRUE),
      sprintf("the cheapest field design meeting them costs %.6g",
              cheapest_hard),
      fixed = TRUE
    )
  }
})

test_that("the cheapest plan can hold a PSU below the threshold", .quiet_single_draw({
  # With 12 units a stratum, holding the large PSU for 2 and drawing one
  # costs 2 visits and 12 units (100), drawing two costs 2 and 20 (140).
  frame <- data.frame(stratum = c("A", "B"), N = c(1360, 64), n_per_psu = 10,
                      cost_psu = 20, cost_ssu = 5)
  psu <- data.frame(stratum = rep(c("A", "B"), each = 4),
                    N = c(1330, 10, 10, 10, 14, 17, 23, 10))
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = 0.3,
                         icc_psu = 0.1)
  fit <- n_alloc(frame, measures = measures, objective = "y", budget = 200,
                 psu = psu, min_n_stratum = 12)
  expect_equal(fit$operational$cost, 200)
  expect_equal(fit$detail$n_int, c(12, 12))
  expect_identical(fit$psu$.certainty_source[fit$psu$certainty],
                   c("planned", "planned"))
  expect_error(
    n_alloc(frame, measures = measures, objective = "y", budget = 199,
            psu = psu, min_n_stratum = 12),
    "'budget' is below the cheapest field design, which costs 200",
    fixed = TRUE
  )
}))

test_that("a register budget fit reports the field design it plans", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = c("literacy", "income"), priority = c(1, 2))
  budget <- 40000
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 objective = objective, budget = budget, psu = z$psu)
  op <- fit$operational
  expect_identical(fit$params$mode, "budget_objective")
  expect_identical(fit$params$budget, budget)
  expect_identical(fit$optimization$budget, budget)
  expect_false(any(c("start", "spend_iterations", "exchange_iterations",
                     "repair_iterations") %in% names(op)))
  expect_identical(op$status, "optimal")
  expect_identical(fit$optimization$search$family, "voluntary")
  expect_equal(fit$optimization$search$gap, 0)
  expect_equal(op$budget_residual, budget - op$cost)
  expect_gte(op$budget_residual, 0)
  expect_lte(fit$params$achieved$cost, budget * (1 + 1e-9))
  expect_equal(fit$optimization$budget_residual,
               budget - fit$params$achieved$cost)
  expect_equal(op$objective_value, sum(op$objective$.contribution))
  expect_equal(op$cost, sum(
    (fit$detail$n_psu_certain + fit$detail$n_psu_draw) * 400 +
      fit$detail$n_int * 20
  ))

  pa <- prec_alloc(fit)
  expect_equal(pa$objective_value, op$objective_value, tolerance = 1e-12)
  expect_equal(pa$params$budget_residual, op$budget_residual)
  expect_equal(pa$psu, fit$psu)
  again <- n_alloc(pa)
  expect_equal(again$operational$objective_value, op$objective_value,
               tolerance = 1e-12)
  expect_lte(again$operational$cost, budget)

  # The same budget buys no worse a design without the hard targets, and
  # the closure family at unit cost is planned the same way.
  free <- n_alloc(fr, measures = z$measures, objective = objective,
                  budget = budget, psu = z$psu)
  expect_lte(free$operational$objective_value, op$objective_value)
  unit <- n_alloc(z$frame, measures = z$measures, objective = objective,
                  budget = 900, psu = z$psu)
  expect_identical(unit$optimization$search$family, "closure")
  expect_lte(unit$operational$cost, 900)
})

test_that("the relaxation is solved at a larger budget when the targets need it", {
  # At the cheapest field design meeting the targets, the continuous model
  # under the same certainty prices visits by size share and needs more.
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = c("literacy", "income"), priority = c(1, 2))
  target <- n_alloc(fr, measures = z$measures, targets = z$targets,
                    psu = z$psu)
  budget <- target$operational$cost
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 objective = objective, budget = budget, psu = z$psu)
  expect_equal(fit$operational$cost, budget)
  expect_equal(fit$operational$budget_residual, 0)
  expect_identical(fit$optimization$budget, budget)
  expect_lt(fit$optimization$budget_residual, 0)
  expect_equal(fit$optimization$budget_residual,
               budget - fit$params$achieved$cost)
})

test_that("a register assessment prices certainty visits in its residual", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = "literacy", priority = 1)
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets, psu = z$psu)
  pa <- prec_alloc(fr, n = fit$detail$n, measures = z$measures,
                   targets = z$targets, objective = objective,
                   budget = 40000, psu = z$psu)
  expect_gt(sum(pa$params$certainty$n_psu_certain), 0)
  expect_equal(pa$params$budget_residual, 40000 - pa$params$achieved$cost)
})

## Register grids

test_that("a register grid reports each row's search status", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets, psu = z$psu)
  expect_warning(
    g <- predict(fit, data.frame(n_per_psu = c(8, 12, 400))),
    "n_per_psu = 400 is infeasible"
  )
  expect_identical(g$.status, c("optimal", "optimal", "infeasible"))
  expect_identical(g$.feasible, c(TRUE, TRUE, FALSE))
  expect_identical(names(g)[ncol(g) - 1:0], c(".status", ".feasible"))
})

test_that("a register budget grid reports the field design", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = c("literacy", "income"), priority = c(1, 2))
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 objective = objective, budget = 40000, psu = z$psu)
  expect_warning(
    g <- predict(fit, data.frame(budget = c(30000, 40000))),
    "budget = 30000 is infeasible: 'budget' cannot fund the precision targets"
  )
  expect_identical(g$.status, c("infeasible", "optimal"))
  expect_identical(g$.feasible, c(FALSE, TRUE))
  expect_equal(g$objective_value[2], fit$operational$objective_value)
  expect_equal(g$cv[2], sqrt(fit$operational$objective_value))
  expect_equal(g$cost_int[2], fit$operational$cost)
  expect_equal(g$n[2], fit$n)
})

test_that("a node limit without a plan is unknown, not infeasible", {
  # At one node the target search ends at 35,100 with a bound of 34,326,
  # so a budget of 34,800 is neither met nor ruled out.
  z <- .psu_fixture(3)
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  ev <- .psu_setup(fr, .check_psu_table(z$psu, fr, z$measures), z$measures,
                   z$targets, objective = "literacy")
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets, psu = z$psu)
  e <- tryCatch(.psu_plan(ev, fit$detail$n / fr$N, 1, budget = 34800),
                error = function(e) e)
  expect_s3_class(e, "svyplan_plan_unknown")
  expect_match(conditionMessage(e), "found in 1 search node, so whether")
  expect_identical(.plan_error_status(e), "unknown")
  for (msg in c("internal error: x",
                "joint allocation solver failed certification: x")) {
    expect_identical(.plan_error_status(simpleError(msg)), "unknown")
  }
  expect_identical(
    .plan_error_status(simpleError("'budget' is below the cheapest")),
    "infeasible"
  )
  cols <- .bethel_predict_cols(TRUE, 1L, certainty = TRUE)
  row <- .bethel_predict_row(NULL, cols, "unknown")
  expect_identical(row$.status, "unknown")
  expect_identical(row$.feasible, NA)
  expect_identical(.bethel_predict_row(NULL, cols)$.feasible, FALSE)
})

test_that("a budget grid row never reports a worse plan than a lower budget", {
  newdata <- data.frame(budget = c(3, 1, 2, 4, 1, 2), zone = c(1, 1, 1, 1, 2, 2))
  values <- data.frame(
    n = c(30, 10, 20, 40, 11, 21),
    objective_value = c(0.5, 0.9, 0.4, 0.45, 0.8, NA),
    n_int = c(31, 11, 21, 41, 12, NA),
    .status = c("limit", "optimal", "optimal", "optimal", "limit", "unknown"),
    .feasible = c(TRUE, TRUE, TRUE, TRUE, TRUE, NA),
    stringsAsFactors = FALSE
  )
  out <- .predict_budget_monotone(newdata, values)
  # Budget 3 stopped at 0.5 takes budget 2's 0.4. Budget 4 is optimal and
  # keeps its own, worse than nothing below it can be.
  expect_equal(out$objective_value[1], 0.4)
  expect_equal(out$n_int[1], 21)
  expect_equal(out$n[1], 30)
  expect_identical(out$.status[1], "limit")
  expect_identical(out[2:4, ], values[2:4, ])
  # In the second group the unknown row takes the lower budget's plan.
  expect_equal(out$objective_value[6], 0.8)
  expect_identical(out$.status[6], "limit")
  expect_identical(out$.feasible[6], TRUE)
  expect_equal(out$n[6], 21)
  expect_identical(.predict_budget_monotone(newdata["zone"], values), values)
})

test_that("predict() reuses a lower budget's plan in a register grid", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = c("literacy", "income"), priority = c(1, 2))
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 objective = objective, budget = 40000, psu = z$psu)
  # Stand-ins for a search stopped at the node limit, one on a worse plan
  # and one with none.
  worse <- fit
  worse$operational$status <- "limit"
  worse$operational$objective_value <- 2 * fit$operational$objective_value
  local_mocked_bindings(n_alloc.default = function(..., budget) {
    if (budget == 40000) return(fit)
    if (budget == 50000) return(worse)
    stop(structure(list(message = "no plan", call = NULL),
                   class = c("svyplan_plan_unknown", "error", "condition")))
  })
  expect_no_warning(
    g <- predict(fit, data.frame(budget = c(60000, 50000, 40000)))
  )
  expect_equal(g$objective_value, rep(fit$operational$objective_value, 3))
  expect_identical(g$.status, c("limit", "limit", "optimal"))
  expect_identical(g$.feasible, rep(TRUE, 3))
  expect_true(is.na(g$n[1]))
  expect_warning(
    g <- predict(fit, data.frame(budget = 60000)),
    "feasibility of budget = 60000 is unknown: no plan"
  )
  expect_identical(g$.status, "unknown")
  expect_identical(g$.feasible, NA)
})

test_that("a register fit prints its continuous solve as a relaxation", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = c("literacy", "income"), priority = c(1, 2))
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 objective = objective, budget = 33760, psu = z$psu)
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl(
    "^continuous relaxation: n = .*, cost = 33823 \\(field design -0.19%\\)$",
    out
  )))
  expect_false(any(grepl("integerizing|one more unit|not binding", out)))
  expect_true(any(grepl("budget sensitivity",
                        utils::capture.output(summary(fit)))))
})

test_that("a register frontier draws the field designs from the register floor", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  objective <- data.frame(name = c("literacy", "income"), priority = c(1, 2))
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets,
                 objective = objective, budget = 40000, psu = z$psu)
  target <- n_alloc(fr, measures = z$measures, targets = z$targets,
                    psu = z$psu)
  grid <- .budget_grid(fit, 3)
  expect_equal(grid, seq(target$operational$cost, 80000, length.out = 3))
  g <- predict(fit, data.frame(budget = grid))
  xy <- .frontier_xy(fit, g)
  expect_equal(xy$cost, g$cost_int)
  expect_equal(xy$cv, g$cv)
  expect_equal(xy$fitted, c(fit$operational$cost,
                            sqrt(fit$operational$objective_value)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit, npoints = 3), fit)
})

test_that("an assessment override keeps a register fit's field design", {
  z <- .psu_fixture(4)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  a <- prec_alloc(fit)
  same <- prec_alloc(fit, alpha = fit$params$alpha)
  expect_equal(same$params$achieved$n, a$params$achieved$n)
  expect_equal(same$detail$.se, a$detail$.se)
  wide <- prec_alloc(fit, alpha = 0.1)
  expect_equal(wide$detail$.se, a$detail$.se)
  expect_equal(wide$detail$.moe, qnorm(0.95) * a$detail$.se)
  expect_identical(wide$psu, fit$psu)
  heavier <- prec_alloc(fit, deff = 1.5)
  expect_equal(heavier$params$achieved$n, a$params$achieved$n)
  expect_equal(heavier$detail$.se, sqrt(1.5) * a$detail$.se)
  expect_equal(prec_alloc(fit, unit_cost = 2)$params$achieved$cost,
               2 * fit$operational$cost)
  # A new allocation is a new design, built from the totals.
  # An override is checked as a supplied allocation would be.
  m <- z$measures
  m$var_ratio_psu <- 2
  expect_error(prec_alloc(fit, measures = m),
               "'var_ratio_psu' is not available with a PSU register")
  m <- z$measures
  m$resp_rate_psu <- 0.9
  expect_error(prec_alloc(fit, measures = m), "'resp_rate_psu' and")
  m <- z$measures
  m$eligible_rate <- 0.5
  for (call in list(quote(prec_alloc(fit, measures = m)),
                    quote(prec_alloc(fit, n = fit$detail$n, measures = m)))) {
    expect_error(eval(call),
                 "'eligible_rate' is not available with a PSU register")
  }
  expect_error(prec_alloc(fit, eligible_rate = 0.5), "unused argument")
  fr <- z$frame
  fr$cost_psu <- 100
  fr$cost_ssu <- 10
  staged <- n_alloc(fr, measures = z$measures, targets = z$targets,
                    psu = z$psu)
  expect_error(prec_alloc(staged, unit_cost = 2),
               "uses stage costs; do not supply 'unit_cost'")
  rebuilt <- prec_alloc(fit, n = fit$detail$n)
  expect_false(isTRUE(all.equal(rebuilt$params$achieved$n,
                                a$params$achieved$n)))
  expect_error(prec_alloc(fit, nonsense = 1), "unused argument")
})

test_that("the cost grain never claims a step the costs do not have", {
  expect_identical(.menu_grain(c(1e9, 1e9 + 0.1)), 0.1)
  expect_identical(.menu_grain(c(100, 10)), 10)
  expect_identical(.menu_grain(c(0.25, 1.5)), 0.25)
  expect_identical(.menu_grain(1 / 3), 0)
  expect_identical(.menu_grain(2^60), 0)
  expect_identical(.menu_grain(c(1, 1.0000005)), 0)
  expect_true(.on_grain(c(30, 40.5), 0.5))
  expect_false(.on_grain(c(30, 40.5), 1))
  expect_false(.on_grain(1000000.5, 1))
  # A grain read off near-whole prices would certify the dearer plan.
  p <- .menu_problem(c(1.0000005, 1), matrix(c(0.1, 0.2)), c(1, 1))
  r <- .menu_search(p, 0, pick = 1L, grain = .menu_grain(p$cost))
  expect_identical(r$value, 1)
  # Fractional PSU sizes put certainty takes, and costs, off the price
  # lattice. The fit is the cheapest plan of the whole menus.
  psu <- data.frame(stratum = rep(c("A", "B"), c(3, 4)),
                    N = c(10.25, 12.5, 17.75, 9.5, 11, 14.25, 8.5))
  frame <- data.frame(stratum = c("A", "B"),
                      N = as.numeric(tapply(psu$N, psu$stratum, sum)[c("A", "B")]),
                      n_per_psu = 4, cost_psu = 10, cost_ssu = 1)
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = c(0.3, 0.4),
                         icc_psu = 0.1)
  targets <- data.frame(name = "y", cv = 0.15)
  fit <- n_alloc(frame, measures = measures, targets = targets, psu = psu)
  ev <- .psu_setup(frame, .check_psu_table(psu, frame, measures), measures,
                   targets)
  menu <- .psu_menu_round(.psu_menu_new(ev), ev, list(price = 1, mu = 0),
                          c(Inf, Inf))
  grid <- expand.grid(lapply(menu, function(m) seq_len(nrow(m$states))))
  cost <- menu[[1]]$states$cost[grid[[1]]] + menu[[2]]$states$cost[grid[[2]]]
  load <- menu[[1]]$C[grid[[1]], 1] + menu[[2]]$C[grid[[2]], 1]
  expect_equal(fit$operational$cost, min(cost[load <= 1 + 1e-9]))
  expect_identical(fit$operational$status, "optimal")
})
