#' The worst-rule caps of a free-take problem, one per component
#'
#' Component `j`'s columns are scaled by `priority_j / total_j^2`, so the cap
#' `priority_j * relvar_j <= t` is the column at the bound `t` for every
#' component. `floor` is the largest of the per-component minima over the
#' stratum boxes, a lower bound on the worst level even though the minima
#' come from different takes. Each stratum's minimum is read from the exact
#' terms at the points where it can lie, as the level of a design is read.
#' @keywords internal
#' @noRd
.freetake_worst_caps <- function(prob) {
  obj <- prob$objective
  w <- obj$priority / obj$total^2
  a <- sweep(obj$a_mat, 2L, w, "*")
  b <- sweep(obj$b_mat, 2L, w, "*")
  d <- sweep(obj$d_mat, 2L, w, "*")
  e <- sweep(obj$e_mat, 2L, w, "*")
  ft <- prob$ft_spec
  floors <- vapply(seq_len(ncol(a)), function(j) {
    pts <- lapply(seq_len(prob$H), function(h) {
      .freetake_stratum_obj_points(
        b[h, j], d[h, j], prob$N_1[h], prob$N[h], prob$L_ult[h],
        ft$take_2_min[h], ft$take_2_max_cont[h]
      )
    })
    at <- vapply(1:3, function(k) {
      x <- vapply(pts, function(p) p$x[min(k, length(p$x))], numeric(1))
      y <- vapply(pts, function(p) p$y[min(k, length(p$y))], numeric(1))
      .freetake_terms(obj$exact, x, y)[, j]
    }, numeric(prob$H))
    w[j] * sum(apply(matrix(at, prob$H), 1L, min))
  }, numeric(1))
  level <- function(x, y) {
    max(colSums(.freetake_terms(obj$exact, x, y)) * w)
  }
  list(a = a, b = b, d = d, e = e, J = ncol(a), floor = max(0, floors),
       level = level, ids = paste0(".worst:", obj$component))
}

#' A free-take problem with every component capped at `vmax`
#' @keywords internal
#' @noRd
.freetake_with_caps <- function(prob, caps, vmax) {
  out <- prob
  out$a <- cbind(prob$a, caps$a)
  out$b <- cbind(prob$b, caps$b)
  out$d <- cbind(prob$d, caps$d)
  out$e <- cbind(prob$e, caps$e)
  out$Vmax <- c(prob$Vmax, rep(vmax, caps$J))
  out$K <- prob$K + caps$J
  out$constraint_ids <- c(prob$constraint_ids, caps$ids)
  out
}

#' The least worst weighted relative variance within a budget, free takes
#'
#' The first stage root-searches the level `t` whose minimum cost, every
#' component capped at `t`, equals the budget. A trial proves its sign by a
#' certified cost, an affordable design or an infeasibility proof. A trial
#' with none of these near the floor raises the lower end. The dual bound of the solve at `t*` gives the lower bound
#' `t* + (dual - max(budget, cost) - slack) / M`, with `M` the sum of the cap
#' multipliers, and the fit is `"optimal"` when the level lies within a
#' relative `kkt_tolerance` of it. The second stage minimizes the weighted
#' sum with every cap at `t*`.
#' @keywords internal
#' @noRd
.freetake_worst_solve <- function(prob, budget, control = .bethel_control()) {
  K <- prob$K
  ft <- prob$ft_spec
  m_min <- ft$take_2_min
  m_max <- ft$take_2_max_cont
  caps <- .freetake_worst_caps(prob)
  J <- caps$J
  opening <- .freetake_budget_open(prob, budget, control)
  if (!is.null(opening$early)) return(opening$early)
  base_x <- opening$x
  base_y <- opening$y
  cost_min <- opening$cost_min
  budget_tol <- opening$budget_tol
  b_eff <- opening$b_eff
  acc <- (1 + control$tolerance)^2 - 1
  accept <- rep(acc, K + J)
  everything <- sum(prob$c1 * prob$N_1 + prob$c2 * prob$N)

  solve_at <- function(t) {
    if (!(t > 0)) {
      return(list(classification = "failed", message = "a zero cap"))
    }
    capped <- .freetake_with_caps(prob, caps, t)
    .freetake_solve_core(
      capped$a, capped$b, capped$d, capped$e, capped$Vmax,
      prob$c1, prob$c2, prob$N_1, prob$N, prob$L_ult, m_min, m_max,
      capped$constraint_ids, control = control, accept = accept
    )
  }
  best <- NULL
  residual <- function(t) {
    s <- solve_at(t)
    if (identical(s$classification, "infeasible")) {
      return(max(1, everything - budget))
    }
    if (!s$classification %in% c("optimal", "feasible")) return(NA_real_)
    if (s$cost <= b_eff) {
      if (is.null(best) || t < best$t) best <<- list(t = t, s = s)
      return(s$cost - budget)
    }
    if (identical(s$classification, "optimal")) return(s$cost - budget)
    NA_real_
  }
  finish <- function(s, t0, binding, at_minimum, iterations = 0L) {
    second <- .freetake_worst_second(prob, caps, t0, budget, control)
    x <- s$x
    y <- s$y
    first_value <- caps$level(x, y)
    if (!is.null(second$x) && sum(prob$c1 * second$x + prob$c2 * second$y) <=
          b_eff && caps$level(second$x, second$y) <= first_value * (1 + 1e-9)) {
      x <- second$x
      y <- second$y
    }
    value <- caps$level(x, y)
    dl <- s$dual_lambda
    M <- if (length(dl) == K + J && all(is.finite(dl))) {
      sum(dl[K + seq_len(J)])
    } else {
      NA_real_
    }
    lower <- caps$floor
    if (isTRUE(M > 0) && is.finite(s$dual_value %||% NA_real_)) {
      slack <- sum(dl[seq_len(K)] * prob$Vmax * acc)
      lower <- max(lower,
                   t0 + (s$dual_value - max(budget, s$cost) - slack) / M)
    }
    lower <- min(lower, value)
    gap <- if (value > 0) (value - lower) / value else 0
    lam <- s$lambda %||% rep(NA_real_, K + J)
    m <- lam[K + seq_len(J)]
    multiplier <- sum(m)
    weights <- m * prob$objective$priority
    weights <- if (length(weights) == J && !anyNA(weights) &&
                     sum(weights) > 0) {
      weights / sum(weights)
    } else {
      rep(NA_real_, J)
    }
    list(
      status = if (gap <= control$kkt_tolerance) "optimal" else "feasible",
      x = x,
      y = y,
      first = list(x = s$x, y = s$y),
      solved = s,
      lambda = lam[seq_len(K)],
      objective_multiplier = if (is.finite(multiplier)) multiplier else NA_real_,
      budget_sensitivity = if (is.finite(multiplier) && multiplier > 0) {
        -1 / multiplier
      } else NA_real_,
      objective_weights = weights,
      objective_bound = t0,
      objective_lower_bound = lower,
      objective_gap = gap,
      budget_binding = binding,
      at_minimum = at_minimum,
      cost = sum(prob$c1 * x + prob$c2 * y),
      root_iterations = iterations,
      second_stage = second$status
    )
  }

  t_hi <- caps$level(base_x, base_y)
  t_lo <- caps$floor
  span <- t_hi - t_lo
  ## Only a budget at or below the cheapest cost takes the cheapest design,
  ## as in the fixed-take search.
  near_min <- budget <= cost_min + budget_tol
  if (budget <= cost_min || !(span > 0)) {
    s <- solve_at(t_hi)
    if (!s$classification %in% c("optimal", "feasible")) {
      s <- list(x = base_x, y = base_y, cost = cost_min,
                lambda = rep(NA_real_, K + J))
    }
    return(finish(s, t_hi, binding = near_min, at_minimum = near_min))
  }
  ## Near a zero floor the core can fail at the tightest level, so a trial
  ## with no sign is raised, never taken as the answer.
  bracket_lo <- t_lo + 1e-10 * span
  f_lo <- residual(bracket_lo)
  while (is.na(f_lo) && bracket_lo < t_lo + 0.1 * span) {
    bracket_lo <- t_lo + 10 * (bracket_lo - t_lo)
    f_lo <- residual(bracket_lo)
  }
  if (is.na(f_lo)) {
    return(list(status = "failed",
                message = "the worst-level bracket has no proven lower end"))
  }
  if (f_lo <= 0) {
    return(finish(best$s, best$t, binding = FALSE, at_minimum = FALSE))
  }
  guarded <- function(t) {
    f <- residual(t)
    if (is.na(f)) stop("a trial with no sign", call. = FALSE)
    f
  }
  root <- tryCatch(
    stats::uniroot(guarded, interval = c(bracket_lo, t_hi), f.lower = f_lo,
                   f.upper = cost_min - budget, tol = span * 1e-12,
                   maxiter = 1000L),
    error = function(e) NULL
  )
  if (!is.null(root)) {
    s <- solve_at(root$root)
    if (s$classification %in% c("optimal", "feasible") && s$cost <= b_eff) {
      return(finish(s, root$root, binding = TRUE, at_minimum = near_min,
                    iterations = root$iter))
    }
  }
  if (is.null(best)) {
    return(list(status = "failed",
                message = "no worst-level trial found an affordable design"))
  }
  finish(best$s, best$t, binding = TRUE, at_minimum = near_min,
         iterations = root$iter %||% 0L)
}

#' The free-take second stage, the weighted sum with every cap at `t_star`
#'
#' The core accepts a column at `Vmax (1 + tol)^2`, so the caps are set at
#' `t_star / (1 + tol)^2` and accept a component at most at `t_star`.
#' @keywords internal
#' @noRd
.freetake_worst_second <- function(prob, caps, t_star, budget, control) {
  if (!(t_star > 0)) return(list(x = NULL, status = "unchanged"))
  held <- .freetake_with_caps(prob, caps,
                              t_star / (1 + control$tolerance)^2)
  fit <- tryCatch(.freetake_budget_solve_continuous(held, budget, control),
                  error = function(e) NULL)
  if (is.null(fit) || !fit$status %in% c("optimal", "feasible")) {
    status <- if (!is.null(fit) && fit$status %in% c("unaffordable", "unknown")) {
      "unchanged"
    } else "failed"
    return(list(x = NULL, status = status))
  }
  list(x = fit$x, y = fit$y, status = fit$status)
}

#' The least worst weighted relative variance over whole free-take designs
#'
#' Bisects on the level as the fixed-take search does. A trial is the exact
#' minimum-cost search over the `(stratum, take_2)` branches with every
#' component capped at the trial level. The lower end is the continuous
#' lower bound at the budget and targets the whole search accepts. The
#' search starts from the cheapest whole design and from the continuous
#' design rounded within the budget, and a second stage takes the least
#' weighted sum at the level found.
#' @keywords internal
#' @noRd
.freetake_whole_worst <- function(prob, budget, control, max_nodes,
                                  continuous = NULL, max_trials = 200L) {
  tolerance <- control$tolerance
  caps <- .freetake_worst_caps(prob)
  bt0 <- .freetake_branch_table(prob, tolerance)
  opening <- .freetake_whole_open(prob, bt0, budget, control, max_nodes, 0)
  limit <- .budget_room(budget, opening$cost)
  used <- opening$nodes
  evaluate <- function(x, take) {
    list(x = as.numeric(x), take = as.numeric(take),
         cost = sum((prob$c1 + prob$c2 * take) * x),
         level = caps$level(x, x * take))
  }
  best <- evaluate(opening$x, bt0$m[opening$branch])
  grain <- .menu_grain(c(prob$c1, prob$c2))
  ## A cap at `t / (1 + tol)^2` accepts a component at most at `t`. A trial
  ## may use half the nodes left, so one hard trial near the level leaves
  ## room for easier ones above it.
  trial <- function(t) {
    left <- max_nodes - used
    if (left < 1 || !(t > 0)) return(list(out = "unknown"))
    capped <- .freetake_with_caps(prob, caps, t / (1 + tolerance)^2)
    bt <- .freetake_branch_table(capped, tolerance)
    if (any(bt$denom <= 0)) return(list(out = "infeasible"))
    B <- length(bt$owner)
    q <- .whole_problem(o1 = bt$cost, o2 = 0, g1 = matrix(0, B, capped$K),
                        g2 = bt$g2, lo = bt$lo, hi = bt$hi,
                        owner = bt$owner, o0 = 0, g0 = bt$g0)
    r <- .whole_search(q, grain = grain, max_nodes = max(1, ceiling(left / 2)))
    used <<- used + r$nodes
    if (identical(r$status, "infeasible")) return(list(out = "infeasible"))
    if (!is.null(r$x)) {
      tr <- evaluate(r$x, bt$m[r$branch])
      if (tr$cost <= limit) return(c(list(out = "feasible"), tr))
    }
    if (identical(r$status, "optimal") || r$bound > .budget_proof(limit)) {
      return(list(out = "infeasible"))
    }
    list(out = "unknown")
  }

  relaxed <- prob
  relaxed$Vmax <- prob$Vmax * (1 + tolerance)^2 * (1 + 1e-12)
  cont <- tryCatch(.freetake_worst_solve(relaxed, .budget_proof(limit),
                                         control),
                   error = function(e) NULL)
  lo <- if (!is.null(cont) && cont$status %in% c("optimal", "feasible")) {
    max(caps$floor, cont$objective_lower_bound)
  } else {
    caps$floor
  }
  if (!is.null(continuous)) {
    rounded <- .freetake_worst_round(prob, caps, bt0, continuous, limit,
                                     tolerance)
    if (!is.null(rounded) && rounded$level < best$level) best <- rounded
  }
  ## An unresolved trial proves nothing. The bracket keeps its lower end,
  ## and the next trial moves between the unresolved level and the best.
  status <- "optimal"
  trials <- 0L
  floor_try <- lo
  while (best$level > 0 && best$level > floor_try * (1 + 1e-9) &&
           max_nodes - used >= 1) {
    if (trials >= max_trials) break
    m <- if (floor_try > 0) sqrt(floor_try * best$level) else best$level / 2
    trials <- trials + 1L
    tr <- trial(m)
    if (identical(tr$out, "feasible")) {
      if (tr$level < best$level) best <- tr[c("x", "take", "cost", "level")]
    } else if (identical(tr$out, "infeasible")) {
      lo <- max(lo, m)
      floor_try <- max(floor_try, m)
    } else {
      floor_try <- m
    }
  }
  if (best$level > lo * (1 + 1e-9)) status <- "limit"
  bound <- min(lo, best$level)
  final <- best
  second_status <- "skipped"
  ## The second stage starts from the design found, which meets every cap
  ## at its own level, so it never searches for a starting design.
  if (max_nodes - used >= 1 && best$level > 0) {
    held <- .freetake_with_caps(prob, caps, best$level / (1 + tolerance)^2)
    bt <- .freetake_branch_table(held, tolerance)
    B <- length(bt$owner)
    q <- .whole_problem(
      o1 = 0, o2 = bt$obj_o2,
      g1 = cbind(.whole_budget_load(bt$cost, limit), matrix(0, B, held$K)),
      g2 = cbind(0, bt$g2), lo = bt$lo, hi = bt$hi, owner = bt$owner,
      o0 = bt$obj_o0, g0 = cbind(0, bt$g0)
    )
    branch <- vapply(seq_len(prob$H), function(h) {
      which(bt$owner == h & abs(bt$m - best$take[h]) < 1e-9)[1L]
    }, integer(1))
    r <- if (anyNA(branch)) NULL else tryCatch(
      .whole_search(q, start = .whole_plan(q, best$x, branch = branch),
                    max_nodes = max_nodes - used),
      error = function(e) NULL
    )
    second_status <- if (is.null(r)) "failed" else r$status
    if (!is.null(r) && !is.null(r$x)) {
      used <- used + r$nodes
      cand <- evaluate(r$x, bt$m[r$branch])
      if (cand$cost <= limit && cand$level <= best$level) final <- cand
    }
  }
  list(
    n_1_int = final$x,
    take_2_int = final$take,
    r = list(status = status),
    search = list(
      family = "whole units (free take)",
      status = status,
      bound = bound,
      gap = if (final$level > 0) max(0, (final$level - bound) / final$level) else 0,
      nodes = used,
      max_nodes = max_nodes,
      states = NA_real_,
      kept = NA_real_,
      trials = trials,
      second_stage = second_status
    )
  )
}

#' A whole free-take design near the continuous worst-rule design
#'
#' Each stratum takes the allowed take closest to its continuous take and
#' the PSU count rounded down within that branch's bounds. One PSU at a time
#' is then added while the budget allows, where it lowers the worst level
#' most. Returns `NULL` when the design misses a target or the budget.
#' @keywords internal
#' @noRd
.freetake_worst_round <- function(prob, caps, bt, continuous, limit,
                                  tolerance) {
  H <- prob$H
  pick <- integer(H)
  for (h in seq_len(H)) {
    idx <- which(bt$owner == h)
    take <- continuous$y[h] / continuous$x[h]
    pick[h] <- idx[which.min(abs(log(bt$m[idx] / take)))]
  }
  take <- bt$m[pick]
  x <- pmin(pmax(floor(continuous$x + 1e-9), bt$lo[pick]), bt$hi[pick])
  unit <- prob$c1 + prob$c2 * take
  if (sum(unit * x) > limit) return(NULL)
  repeat {
    room <- limit - sum(unit * x)
    open <- which(x < bt$hi[pick] & unit <= room)
    if (length(open) == 0L) break
    score <- vapply(open, function(h) {
      z <- x
      z[h] <- z[h] + 1
      caps$level(z, z * take)
    }, numeric(1))
    h <- open[which.min(score)]
    x[h] <- x[h] + 1
  }
  if (prob$K > 0L) {
    v <- .freetake_var(prob$a, prob$b, prob$d, prob$e, x, x * take)
    if (any(v > prob$Vmax * (1 + tolerance)^2)) return(NULL)
  }
  list(x = x, take = take, cost = sum(unit * x),
       level = caps$level(x, x * take))
}

#' Slopes of a free-take fit, from verified multipliers
#'
#' A slope is reported only for a fit classified `"optimal"` whose
#' multipliers are verified against the full free-take problem. A candidate
#' vector, the solve's own or the dual's, must be nonnegative, vanish on
#' slack rows, and give a free-take dual value equal to the primal cost, so
#' that the design minimizes the full Lagrangian and not only a fixed-take
#' restriction of it. `.sensitivity` is recomputed from the verified vector.
#' A row has no slope of its own when its gradient is in a dependency on
#' the face the active bounds leave, take bounds included, or when two
#' verified vectors disagree on it. `cols` holds the objective columns of a
#' budget fit and is `NULL` for a minimum-cost fit.
#' @keywords internal
#' @noRd
.freetake_slopes <- function(prob, fit, constraints, cols, control) {
  K <- prob$K
  n_rows <- nrow(constraints)
  active <- constraints$.binding
  budget_mode <- !is.null(cols)
  ## An uncertified design can leave slack a target that binds at the
  ## optimum, so a withheld fit reports no slope on any row, zero included.
  withheld <- function(status) {
    out <- list(sensitivity = rep(NA_real_, n_rows), status = status)
    if (budget_mode) {
      out$objective <- out$sensitivity
      out$budget <- NA_real_
    }
    out
  }
  if (budget_mode && isTRUE(fit$at_minimum) && isTRUE(fit$budget_binding)) {
    return(withheld("budget_at_minimum"))
  }
  if (!identical(fit$status, "optimal")) return(withheld("not_optimal"))
  s <- fit$solved
  if (is.null(s$x)) return(withheld("multipliers_not_verified"))
  x <- s$x
  y <- s$y
  J <- if (budget_mode) ncol(cols$a) else 0L
  a <- cbind(prob$a, if (budget_mode) cols$a)
  b <- cbind(prob$b, if (budget_mode) cols$b)
  d <- cbind(prob$d, if (budget_mode) cols$d)
  e <- cbind(prob$e, if (budget_mode) cols$e)
  vmax <- c(prob$Vmax, if (budget_mode) cols$vmax)
  if (K + J == 0L) return(list(sensitivity = numeric(0), status = "available"))
  ft <- prob$ft_spec
  cost <- sum(prob$c1 * x + prob$c2 * y)
  resid <- .freetake_var(a, b, d, e, x, y) / vmax - 1
  dual_at <- function(v) {
    A <- drop(a %*% v)
    B <- drop(b %*% v)
    D <- drop(d %*% v)
    value <- 0
    for (h in seq_len(prob$H)) {
      value <- value + .freetake_stratum_min(
        prob$c1[h], prob$c2[h], A[h], B[h], D[h], prob$N_1[h], prob$N[h],
        prob$L_ult[h], ft$take_2_min[h], ft$take_2_max_cont[h]
      )$val
    }
    value - sum(v * (vmax + colSums(e)))
  }
  verified <- function(v) {
    if (length(v) != K + J || anyNA(v) || any(!is.finite(v)) ||
          any(v < -1e-12 * max(1, abs(v)))) {
      return(FALSE)
    }
    v <- pmax(v, 0)
    if (any(v * vmax * abs(resid) > control$kkt_tolerance * cost &
              resid < -1e-6)) {
      return(FALSE)
    }
    abs(cost - dual_at(v)) <= control$kkt_tolerance * cost
  }
  candidates <- Filter(verified, list(s$lambda, s$dual_lambda))
  if (length(candidates) == 0L) return(withheld("multipliers_not_verified"))
  v <- pmax(candidates[[1L]], 0)
  disagree <- logical(K + J)
  if (length(candidates) == 2L) {
    w <- pmax(candidates[[2L]], 0)
    disagree <- abs(v - w) * vmax > 1e-6 * cost
  }
  rows <- c(which(active),
            if (budget_mode) K + which(resid[K + seq_len(J)] >= -1e-6))
  G <- .freetake_face_gradients(prob, a[, rows, drop = FALSE],
                                b[, rows, drop = FALSE],
                                d[, rows, drop = FALSE], x, y)
  shared <- .ambiguous_columns(G) | disagree[rows]
  n_obj <- length(rows) - sum(active)
  target_shared <- shared[seq_len(sum(active))]
  lambda <- v[seq_len(K)][active]
  meta <- constraints[active, , drop = FALSE]
  q <- prob$q_target[active]
  sens <- ifelse(
    meta$.metric == "moe",
    -2 * lambda * prob$domain_N[active]^2 * meta$.target / q^2,
    -2 * lambda * meta$.target * prob$total[active]^2 *
      ifelse(meta$.metric == "rmoe", 1 / q^2, 1)
  )
  sens[target_shared] <- NA_real_
  sensitivity <- numeric(n_rows)
  sensitivity[active] <- sens
  if (!budget_mode) return(list(sensitivity = sensitivity, status = "available"))
  out <- list(sensitivity = sensitivity, objective = numeric(n_rows),
              budget = NA_real_)
  if (!isTRUE(fit$budget_binding)) {
    out$status <- "budget_not_binding"
    return(out)
  }
  mu <- sum(v[K + seq_len(J)])
  out$status <- if (!(mu > 0)) {
    "objective_multiplier_not_positive"
  } else if (any(disagree[K + seq_len(J)]) ||
             .price_undetermined(G, c(rep(0, sum(active)), rep(1, n_obj)))) {
    "not_identifiable"
  } else {
    "available"
  }
  if (identical(out$status, "available")) {
    out$objective[active] <- sens / mu
    out$budget <- -1 / mu
  } else {
    out$objective[active] <- NA_real_
  }
  out
}

#' Gradients of active free-take rows on the face their bounds leave
#'
#' The variables are `(x_h, y_h)`. Each active bound, `x` at 1 or `N_1`, `y`
#' at its floor or `N`, or `y` at `m_min x` or `m_max x`, removes one
#' direction. Returns one column per row, the row's gradient on what
#' remains.
#' @keywords internal
#' @noRd
.freetake_face_gradients <- function(prob, a, b, d, x, y) {
  H <- prob$H
  r <- ncol(a)
  if (r == 0L) return(matrix(0, 2L * H, 0L))
  ft <- prob$ft_spec
  near <- function(u, v) abs(u - v) <= 1e-7 * max(1, abs(v))
  bound_rows <- list()
  bound <- function(g) bound_rows[[length(bound_rows) + 1L]] <<- g
  for (h in seq_len(H)) {
    if (near(x[h], 1) || near(x[h], prob$N_1[h])) {
      g <- numeric(2L * H)
      g[h] <- 1
      bound(g)
    }
    if (near(y[h], prob$L_ult[h]) || near(y[h], prob$N[h])) {
      g <- numeric(2L * H)
      g[H + h] <- 1
      bound(g)
    }
    for (m in unique(c(ft$take_2_min[h], ft$take_2_max_cont[h]))) {
      if (near(y[h], m * x[h])) {
        g <- numeric(2L * H)
        g[h] <- -m
        g[H + h] <- 1
        bound(g)
      }
    }
  }
  ## Row k holds dV_k / dx_h, then dV_k / dy_h.
  G <- cbind(t((-a + d * (y / x)) / x^2), t(-b / y^2 - d / x))
  if (length(bound_rows)) {
    C <- do.call(rbind, bound_rows)
    qc <- qr(t(C))
    if (qc$rank >= 2L * H) return(matrix(0, 0L, r))
    Q <- qr.Q(qc, complete = TRUE)
    G <- G %*% Q[, (qc$rank + 1L):(2L * H), drop = FALSE]
  }
  t(G)
}

#' Which active free-take rows have no multiplier of their own on their face
#' @keywords internal
#' @noRd
.freetake_shared_rows <- function(prob, a, b, d, x, y) {
  .ambiguous_columns(.freetake_face_gradients(prob, a, b, d, x, y))
}
