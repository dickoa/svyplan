#' One reciprocal cap per objective component under the worst rule
#'
#' Component `j` has variance `B_j + sum_h A_hj / n_h` and relative variance
#' that over `total_j^2`. The cap `priority_j * relvar_j <= t` is
#' `sum_h A_hj / n_h <= t * rate_j - B_j` with `rate_j = total_j^2 /
#' priority_j`. A component with no variance term is constant at zero and
#' gets no cap.
#' @keywords internal
#' @noRd
.bethel_worst_caps <- function(problem) {
  obj <- problem$objective
  live <- which(colSums(obj$A) > 0)
  rate <- obj$total^2 / obj$priority
  contribution <- function(allocation) {
    v <- .bethel_objective_variance(obj, allocation)
    v[v < 0] <- 0
    v / rate
  }
  list(
    A = obj$A[, live, drop = FALSE],
    B = obj$B[live],
    rate = rate[live],
    live = live,
    J = length(obj$priority),
    ids = paste0(".worst:", obj$component[live]),
    contribution = contribution,
    level = function(allocation) max(contribution(allocation)),
    bound = function(t) t * rate[live] - obj$B[live]
  )
}

#' The least worst weighted relative variance within a budget, continuous
#'
#' The first stage root-searches the level `t` whose minimum cost, with every
#' component capped at `t`, equals the budget. The minimum cost does not
#' rise as `t` grows. The second stage minimizes the priority-weighted sum
#' with the caps held at `t*`, which spends any budget the first stage leaves
#' and returns a design no other affordable design dominates. The prices
#' come from the first stage. `objective_multiplier` is `-dC/dt`, so the
#' slopes of a budget fit keep their form.
#' @keywords internal
#' @noRd
.bethel_worst_solve <- function(problem, budget, control = .bethel_control()) {
  A <- problem$A
  bound <- problem$bound
  cost <- problem$cost
  lower <- problem$lower
  upper <- problem$upper
  K <- ncol(A)
  caps <- .bethel_worst_caps(problem)
  accept <- .bethel_canonical_tolerance(problem, control$tolerance)
  ids <- c(problem$constraint_ids, caps$ids)
  cap_accept <- function(t) {
    vmax <- t * caps$rate
    pmax(0, vmax * ((1 + control$tolerance)^2 - 1) / caps$bound(t))
  }
  solve_at <- function(t) {
    b <- c(bound, caps$bound(t))
    if (!all(b > 0)) {
      return(list(classification = "failed", message = "a zero cap"))
    }
    .bethel_solve(
      cbind(A, caps$A), b, cost, lower, upper,
      control = control, constraint_ids = ids,
      accept = c(accept, cap_accept(t)),
      floors = problem$floors
    )
  }

  if (K > 0L) {
    base <- .bethel_solve(
      A, bound, cost, lower, upper,
      control = control, constraint_ids = problem$constraint_ids,
      accept = accept,
      floors = problem$floors
    )
    if (identical(base$classification, "infeasible")) {
      return(list(status = "infeasible_targets", infeasible = base$infeasible))
    }
    if (!identical(base$classification, "optimal")) {
      return(list(status = "failed", message = base$message))
    }
    base_allocation <- base$allocation
  } else {
    base_allocation <- .bethel_floor_cheapest(problem)
  }
  cost_min <- sum(cost * base_allocation)
  cost_max <- sum(cost * upper)
  budget_tol <- control$tolerance * abs(budget)
  b_eff <- .budget_room(budget, cost_min)
  if (budget < cost_min - budget_tol) {
    return(list(status = "unaffordable", cost_min = cost_min,
                shortfall = cost_min - budget))
  }
  ## The largest level proven out of reach: the floor, or a certified trial
  ## whose dual bound on the minimum cost exceeds the budget.
  proven <- caps$level(upper)
  finish <- function(first, t_star, solved, binding, at_minimum,
                     iterations = 0L, precision = 0) {
    lam <- solved$lambda %||% rep(NA_real_, K + length(caps$live))
    m <- lam[K + seq_along(caps$live)]
    M <- sum(m * caps$rate)
    weights <- rep(0, caps$J)
    weights[caps$live] <- m * caps$rate * problem$objective$priority[caps$live]
    weights <- if (anyNA(weights) || !(sum(weights) > 0)) {
      rep(NA_real_, caps$J)
    } else {
      weights / sum(weights)
    }
    at_cap <- caps$contribution(first)[caps$live] >=
      t_star * (1 - 1e-5)^2
    ## A solve stops just short of an upper bound, so a census sits a hair
    ## below N. Rounding such counts up lowers every variance and is kept
    ## when it stays within the budget.
    snap <- function(n) {
      near <- abs(n - upper) <= control$tolerance * upper
      snapped <- n
      snapped[near] <- upper[near]
      if (sum(cost * snapped) <= b_eff) snapped else n
    }
    first <- snap(first)
    second <- .bethel_worst_second(problem, caps, t_star, budget, control)
    allocation <- first
    if (!is.null(second$allocation)) {
      candidate <- snap(second$allocation)
      if (sum(cost * candidate) <= b_eff &&
            caps$level(candidate) <= t_star * (1 + 1e-9)) {
        allocation <- candidate
      }
    }
    list(
      status = "optimal",
      allocation = allocation,
      first = first,
      base_allocation = base_allocation,
      solved = solved,
      lambda = lam[seq_len(K)],
      objective_multiplier = if (is.finite(M)) M else NA_real_,
      budget_sensitivity = if (is.finite(M) && M > 0) -1 / M else NA_real_,
      objective_weights = weights,
      objective_cols = caps$A[, at_cap, drop = FALSE],
      objective_rates = caps$rate[at_cap],
      certified_lower = proven,
      objective_bound = t_star,
      root_precision = precision,
      budget_binding = binding,
      at_minimum = at_minimum,
      cost = sum(cost * allocation),
      root_iterations = iterations,
      second_stage = second$status
    )
  }

  if (cost_max <= b_eff) {
    return(finish(upper, caps$level(upper), NULL, binding = FALSE,
                  at_minimum = FALSE))
  }
  t_hi <- caps$level(base_allocation)
  t_lo <- caps$level(upper)
  span <- t_hi - t_lo
  certified <- function(s) identical(s$classification, "optimal")
  ## A budget within the tolerance above the cheapest cost still buys a
  ## lower level when a large fixed cost sets the tolerance, so only a
  ## budget at or below it takes the cheapest design.
  near_min <- budget <= cost_min + budget_tol
  if (budget <= cost_min || !(span > 0)) {
    solved <- solve_at(t_hi)
    return(finish(base_allocation, t_hi, if (certified(solved)) solved,
                  binding = near_min, at_minimum = near_min))
  }
  ## An uncertified trial still has a sign when its dual bound, a lower bound
  ## on the minimum cost for any multipliers, exceeds the budget. Near a zero
  ## floor the solver often stops short of certifying such a trial.
  ## The smallest level whose certified design fits, the answer when the
  ## root's own design overspends.
  affordable <- NULL
  residual <- function(t) {
    s <- solve_at(t)
    ## Any multipliers give a lower bound on the minimum cost, so a trial
    ## proves its level out of reach whether or not it certified.
    if ((s$dual_value %||% -Inf) * (s$cost_scale %||% 1) > budget) {
      proven <<- max(proven, t)
    }
    if (certified(s) && sum(cost * s$allocation) <= b_eff &&
          (is.null(affordable) || t < affordable$t)) {
      affordable <<- list(t = t, s = s)
    }
    if (certified(s)) return(sum(cost * s$allocation) - budget)
    least <- (s$dual_value %||% NA_real_) * (s$cost_scale %||% NA_real_)
    if (is.finite(least) && least > b_eff) return(least - budget)
    ## A design that meets every cap and fits the budget proves the level
    ## affordable, certified or not.
    spent <- if (is.null(s$allocation)) Inf else sum(cost * s$allocation)
    if (isTRUE(s$primal_residual <= control$kkt_tolerance) && spent <= budget) {
      return(spent - budget)
    }
    NA_real_
  }
  ## The floor itself, such as a census under unit FPC, is tried first, so
  ## an affordable floor is returned exactly.
  at_floor <- solve_at(t_lo)
  if (certified(at_floor) && sum(cost * at_floor$allocation) <= b_eff) {
    return(finish(at_floor$allocation, t_lo, at_floor, binding = FALSE,
                  at_minimum = FALSE))
  }
  bracket_lo <- t_lo + 1e-10 * span
  f_lo <- residual(bracket_lo)
  if (is.na(f_lo)) {
    return(list(status = "failed",
                message = "the worst-level bracket has no certified lower end"))
  }
  if (f_lo <= 0) {
    solved <- solve_at(bracket_lo)
    return(finish(solved$allocation, bracket_lo, solved, binding = FALSE,
                  at_minimum = FALSE))
  }
  failed <- FALSE
  guarded <- function(t) {
    f <- residual(t)
    if (is.na(f)) {
      failed <<- TRUE
      stop("an uncertified trial", call. = FALSE)
    }
    f
  }
  root <- tryCatch(
    stats::uniroot(guarded, interval = c(bracket_lo, t_hi), f.lower = f_lo,
                   f.upper = cost_min - budget, tol = span * 1e-12,
                   maxiter = 1000L),
    error = function(e) NULL
  )
  if (is.null(root)) {
    return(list(status = "failed", message = if (failed) {
      "a worst-level trial failed certification"
    } else {
      "the worst-level root search did not converge"
    }))
  }
  t_star <- root$root
  solved <- solve_at(t_star)
  if (certified(solved) && sum(cost * solved$allocation) > b_eff) {
    t_star <- min(t_hi, root$root + root$estim.prec)
    solved <- solve_at(t_star)
  }
  if ((!certified(solved) || sum(cost * solved$allocation) > b_eff) &&
        !is.null(affordable)) {
    t_star <- affordable$t
    solved <- affordable$s
  }
  if (!certified(solved) || sum(cost * solved$allocation) > b_eff) {
    return(list(status = "failed",
                message = solved$message %||% "the worst-level design exceeds the budget"))
  }
  finish(solved$allocation, t_star, solved, binding = TRUE,
         at_minimum = near_min, iterations = root$iter,
         precision = root$estim.prec)
}

#' The second stage, the weighted sum within the budget with caps at `t*`
#'
#' Returns the allocation and a status: `"optimal"`, `"unchanged"` when the
#' caps leave no room (the first-stage design is then the only one), or
#' `"failed"`.
#' @keywords internal
#' @noRd
.bethel_worst_second <- function(problem, caps, t_star, budget, control) {
  held <- problem
  held$A <- cbind(problem$A, caps$A)
  held$B <- c(problem$B, caps$B)
  ## The solve accepts a column at Vmax (1 + tol)^2, so this cap accepts a
  ## component at most at `t_star`, and the second stage never worsens it.
  vmax <- t_star * caps$rate / (1 + control$tolerance)^2
  if (!all(vmax - caps$B > 0)) {
    return(list(allocation = NULL, status = "unchanged"))
  }
  held$bound <- c(problem$bound, vmax - caps$B)
  held$constraint_ids <- c(problem$constraint_ids, caps$ids)
  fit <- tryCatch(.bethel_budget_solve(held, budget, control),
                  error = function(e) NULL)
  if (is.null(fit) || identical(fit$status, "unaffordable")) {
    return(list(allocation = NULL, status = "unchanged"))
  }
  if (!identical(fit$status, "optimal")) {
    return(list(allocation = NULL, status = "failed"))
  }
  list(allocation = fit$allocation, status = "optimal")
}

#' The least worst weighted relative variance over whole designs
#'
#' Bisects on the level `t` between a proven lower bound and an affordable
#' whole design. A trial is the exact minimum-cost search with every
#' component capped at `t`. It is feasible when it returns a design within
#' the budget, infeasible when its search proves the least whole cost above
#' the budget, and unresolved otherwise, which ends the bisection with
#' status `"limit"`. The lower bound is the continuous level at the budget
#' and targets the whole search accepts, or the level of the whole upper
#' bounds. `max_nodes` bounds all trials and the second stage together.
#' @keywords internal
#' @noRd
.bethel_whole_worst <- function(problem, budget, control = .bethel_control(),
                                max_nodes = 5000, max_trials = 200L,
                                continuous = NULL, solved = NULL) {
  tolerance <- control$tolerance
  cost <- problem$cost
  H <- length(cost)
  K <- ncol(problem$A)
  caps <- .bethel_worst_caps(problem)
  opening <- .bethel_whole_start(problem, budget, control, max_nodes, 0)
  limit <- .budget_room(budget, opening$cost)
  used <- opening$nodes
  accept <- 1 + .bethel_canonical_tolerance(problem, tolerance)
  target_load <- sweep(problem$A, 2L, problem$bound * accept, "/") *
    .whole_target_scale()
  ## A cap at `t` accepts a component whose weighted relative variance is
  ## at most `t`, the line the assessment reports.
  cap_load <- function(t) {
    sweep(caps$A, 2L, t * caps$rate - caps$B, "/") * .whole_target_scale()
  }
  evaluate <- function(x) {
    list(x = .whole_storage(x), cost = sum(cost * x), level = caps$level(x))
  }
  best <- evaluate(opening$start)
  trials <- 0L
  trial <- function(t) {
    left <- max_nodes - used
    if (left < 1 || !all(t * caps$rate - caps$B > 0)) {
      return(list(out = "unknown"))
    }
    ## With floors, the continuous solve's multipliers, scaled to this
    ## level's loads, start the search.
    held <- problem
    held$A <- cbind(problem$A, caps$A)
    held$bound <- c(problem$bound, t * caps$rate - caps$B)
    q <- .bethel_whole_problem(problem, cost, 0,
                               matrix(0, H, K + length(caps$live)),
                               cbind(target_load, cap_load(t)))
    lambda <- if (!is.null(problem$floors)) {
      .bethel_whole_lambda(held, q, solved,
                           c(accept, rep(1, length(caps$live))))
    }
    r <- .whole_search(q, grain = .menu_grain(cost), max_nodes = left,
                       lambda = lambda)
    used <<- used + r$nodes
    if (identical(r$status, "infeasible")) return(list(out = "infeasible"))
    if (!is.null(r$x) && sum(cost * r$x) <= limit) {
      return(c(list(out = "feasible"), evaluate(r$x)))
    }
    if (identical(r$status, "optimal") || r$bound > .budget_proof(limit)) {
      return(list(out = "infeasible"))
    }
    list(out = "unknown")
  }

  floor_level <- caps$level(problem$upper_int)
  if (sum(cost * problem$upper_int) <= limit) {
    best <- evaluate(problem$upper_int)
    lo <- floor_level
  } else {
    relaxed <- problem
    relaxed$bound <- problem$bound * accept * (1 + 1e-12)
    cont <- tryCatch(.bethel_worst_solve(relaxed, .budget_proof(limit),
                                         control),
                     error = function(e) NULL)
    ## Only a proven level bounds the whole designs. The solve's own t* can
    ## overstate it, as when a budget within the tolerance of the cheapest
    ## cost still buys cheaper strata.
    lo <- if (!is.null(cont) && identical(cont$status, "optimal")) {
      max(floor_level, cont$certified_lower)
    } else {
      floor_level
    }
  }
  if (!is.null(continuous)) {
    rounded <- .bethel_worst_round(problem, caps, continuous, limit, tolerance)
    if (!is.null(rounded) && rounded$level < best$level) best <- rounded
  }
  status <- "optimal"
  while (best$level > 0 && best$level > lo * (1 + 1e-9)) {
    if (trials >= max_trials) {
      status <- "limit"
      break
    }
    m <- if (lo > 0) sqrt(lo * best$level) else best$level / 2
    trials <- trials + 1L
    tr <- trial(m)
    if (identical(tr$out, "feasible")) {
      if (tr$level < best$level) best <- tr[c("x", "cost", "level")]
    } else if (identical(tr$out, "infeasible")) {
      lo <- m
    } else {
      status <- "limit"
      break
    }
  }
  bound <- min(lo, best$level)
  primary <- best$x
  second <- if (max_nodes - used >= 1) {
    .bethel_whole_worst_second(problem, caps, best$level, limit, control,
                               max_nodes - used, start = best$x)
  } else {
    list(allocation = NULL, status = "skipped", nodes = 0)
  }
  used <- used + second$nodes
  allocation <- primary
  if (!is.null(second$allocation) &&
        sum(cost * second$allocation) <= limit &&
        caps$level(second$allocation) <= best$level) {
    allocation <- second$allocation
  }
  evaluation <- .precision_from_allocation(problem, allocation,
                                           tolerance = tolerance)
  if (!evaluation$all_pass || !evaluation$bounds_pass ||
        !.bethel_floors_met(problem, allocation)) {
    stop("internal error: integer allocation failed final validation",
         call. = FALSE)
  }
  if (evaluation$cost > limit) {
    stop("internal error: integer allocation exceeded the budget",
         call. = FALSE)
  }
  objective <- .bethel_objective_value(problem, allocation)
  list(
    allocation = .whole_storage(allocation),
    cost = evaluation$cost,
    constraints = evaluation$constraints,
    all_pass = evaluation$all_pass,
    objective = objective$components,
    objective_value = objective$value,
    budget_residual = budget - evaluation$cost,
    status = status,
    search = list(
      family = "whole units",
      status = status,
      bound = bound,
      gap = if (objective$value > 0) {
        max(0, (objective$value - bound) / objective$value)
      } else 0,
      nodes = used,
      max_nodes = max_nodes,
      states = NA_real_,
      kept = NA_real_,
      trials = trials,
      second_stage = second$status
    )
  )
}

#' A whole design near the continuous worst-rule design
#'
#' Rounds the continuous design down within the whole bounds, then adds one
#' unit at a time while the budget allows, to the stratum that lowers the
#' worst level most (the weighted sum breaks ties). Returns `NULL` when the
#' rounded design misses a target or the budget.
#' @keywords internal
#' @noRd
.bethel_worst_round <- function(problem, caps, continuous, limit, tolerance) {
  cost <- problem$cost
  x <- pmin(pmax(floor(continuous + 1e-9), problem$lower_int),
            problem$upper_int)
  if (sum(cost * x) > limit) return(NULL)
  sum_at <- function(n) sum(caps$contribution(n))
  repeat {
    room <- limit - sum(cost * x)
    open <- which(x < problem$upper_int & cost <= room)
    if (length(open) == 0L) break
    score <- vapply(open, function(h) {
      y <- x
      y[h] <- y[h] + 1
      c(caps$level(y), sum_at(y))
    }, numeric(2))
    pick <- open[order(score[1, ], score[2, ])[1L]]
    x[pick] <- x[pick] + 1
  }
  if (!.precision_from_allocation(problem, x, tolerance = tolerance)$all_pass ||
        !.bethel_floors_met(problem, x)) {
    return(NULL)
  }
  list(x = .whole_storage(x), cost = sum(cost * x), level = caps$level(x))
}

#' The whole second stage, the weighted sum with every cap at the level found
#'
#' The caps accept a component at most at `level`, so the returned design
#' never worsens the worst level. A budget search that cannot run keeps the
#' first-stage design.
#' @keywords internal
#' @noRd
.bethel_whole_worst_second <- function(problem, caps, level, budget, control,
                                       max_nodes, start) {
  held <- problem
  held$A <- cbind(problem$A, caps$A)
  held$B <- c(problem$B, caps$B)
  ## The search accepts a column at Vmax (1 + tol)^2, so a cap of
  ## level * rate / (1 + tol)^2 accepts exactly the components at `level`.
  vmax <- level * caps$rate / (1 + control$tolerance)^2
  held$bound <- c(problem$bound, vmax - caps$B)
  held$constraint_ids <- c(problem$constraint_ids, caps$ids)
  ## The search starts from the design found, which meets every cap at its
  ## own level, so it never searches for a starting design.
  H <- length(problem$cost)
  accept <- 1 + .bethel_canonical_tolerance(held, control$tolerance)
  load <- sweep(held$A, 2L, held$bound * accept, "/") * .whole_target_scale()
  q <- .bethel_whole_problem(
    problem, 0, problem$objective$q_h,
    cbind(.whole_budget_load(problem$cost, budget),
          matrix(0, H, ncol(held$A))),
    cbind(0, load)
  )
  r <- tryCatch(.whole_search(q, start = start, max_nodes = max_nodes),
                error = function(e) NULL)
  if (is.null(r) || is.null(r$x)) {
    return(list(allocation = NULL, status = "failed", nodes = r$nodes %||% 0))
  }
  list(allocation = .whole_storage(r$x), status = r$status, nodes = r$nodes)
}
