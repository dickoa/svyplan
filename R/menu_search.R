#' One state per group under column limits
#'
#' The integer problem an exact allocation reduces to. Each group, a
#' stratum, takes one of its states. A state has a cost and a load on each
#' column, and a plan minimizes its total cost with every column's load at
#' most 1. Target mode reads columns as shares of the targets' ceilings and
#' cost as field cost. Budget mode reads cost as the objective and adds a
#' column for field cost over the budget. `owner` gives each state's group,
#' numbered from 1. States in group order carry the ends of their runs,
#' which lets the group minima run in one pass.
#' @keywords internal
#' @noRd
.menu_problem <- function(cost, C, owner) {
  C <- matrix(as.double(C), nrow = length(cost))
  owner <- as.integer(owner)
  list(cost = as.double(cost), C = C, owner = owner,
       groups = split(seq_along(owner), owner), H = max(owner),
       ends = if (!is.unsorted(owner)) as.integer(cumsum(tabulate(owner))))
}

#' Each run's smallest value, its position and the second smallest
#'
#' `x` is a vector or a matrix whose rows fall in consecutive runs ending
#' at `ends`. Each result has one row per run and one column per column of
#' `x`, a vector for a vector.
#' @keywords internal
#' @noRd
.menu_group_min <- function(x, ends) {
  out <- .Call(C_menu_group_min, x, ends)
  if (is.null(dim(x))) out <- lapply(out, as.vector)
  out
}

#' Which states of `ord`, visited in order, no earlier kept state dominates
#'
#' A state is dominated when a state kept before it loads no column more.
#' Visited in cost order, the kept state also costs no more.
#' @keywords internal
#' @noRd
.menu_dominated <- function(C, ord) {
  .Call(C_menu_dominated, C, as.integer(ord))
}

.menu_argmin <- function(r, groups) {
  vapply(groups, function(j) j[which.min(r[j])], 1L, USE.NAMES = FALSE)
}

#' Lagrangian bound, each group's best state and the subgradient
#' @keywords internal
#' @noRd
.menu_lagrangian <- function(p, lambda) {
  r <- p$cost + drop(p$C %*% lambda)
  best <- if (is.null(p$ends)) {
    .menu_argmin(r, p$groups)
  } else {
    .menu_group_min(r, p$ends)$arg
  }
  list(L = sum(r[best]) - sum(lambda), r = r, best = best,
       g = colSums(p$C[best, , drop = FALSE]) - 1)
}

.menu_feasible <- function(p, pick) {
  all(colSums(p$C[pick, , drop = FALSE]) <= 1 + 1e-9)
}

#' Move groups to cheaper states while every column holds
#'
#' The largest saving first, until none is left.
#' @keywords internal
#' @noRd
.menu_step_down <- function(p, pick) {
  load <- colSums(p$C[pick, , drop = FALSE])
  repeat {
    best <- NULL
    for (h in seq_len(p$H)) {
      j <- p$groups[[h]]
      j <- j[p$cost[j] < p$cost[pick[h]]]
      if (!length(j)) next
      new <- sweep(p$C[j, , drop = FALSE], 2, load - p$C[pick[h], ], "+")
      ok <- j[rowSums(new <= 1 + 1e-9) == ncol(p$C)]
      if (!length(ok)) next
      s <- ok[which.min(p$cost[ok])]
      save <- p$cost[pick[h]] - p$cost[s]
      if (is.null(best) || save > best$save) {
        best <- list(h = h, s = s, save = save)
      }
    }
    if (is.null(best)) break
    load <- load - p$C[pick[best$h], ] + p$C[best$s, ]
    pick[best$h] <- best$s
  }
  pick
}

#' A plan from multipliers, or `NULL`
#'
#' The multipliers are doubled until every group's best state together
#' meets every column, and the plan is then stepped down.
#' @keywords internal
#' @noRd
.menu_plan_from <- function(p, lambda) {
  if (all(lambda == 0)) lambda <- rep(1, length(lambda))
  for (i in seq_len(200L)) {
    b <- .menu_lagrangian(p, lambda)$best
    if (.menu_feasible(p, b)) return(.menu_step_down(p, b))
    lambda <- lambda * 2
  }
  NULL
}

#' Multipliers by subgradient, with the best plan met on the way
#'
#' A cold start runs up to 1,500 steps from a uniform guess. Given
#' `lambda`, the multipliers are taken as close and a warm run takes up to
#' `iters` smaller steps. `start`, a plan meeting every column, is kept
#' when no cheaper plan is found. `ub` is `Inf` when no plan is known.
#' @keywords internal
#' @noRd
.menu_subgradient <- function(p, lambda = NULL, start = NULL,
                              iters = if (is.null(lambda)) 1500L else 100L) {
  K <- ncol(p$C)
  warm <- !is.null(lambda)
  if (!warm) lambda <- rep(abs(mean(p$cost)) * p$H / K + 1e-12, K)
  theta <- if (warm) 0.5 else 2
  patience <- if (warm) 10L else 25L
  pick <- .menu_plan_from(p, lambda)
  if (!is.null(start) &&
        (is.null(pick) || sum(p$cost[start]) < sum(p$cost[pick]))) {
    pick <- start
  }
  ub <- if (is.null(pick)) Inf else sum(p$cost[pick])
  # A step needs a finite target. Every plan costs less than this one.
  ceiling_cost <- sum(vapply(p$groups, function(j) max(p$cost[j]), 0)) + 1
  best_L <- -Inf
  best_lambda <- lambda
  stall <- 0L
  for (i in seq_len(iters)) {
    l <- .menu_lagrangian(p, lambda)
    if (!is.finite(best_L) || l$L > best_L + 1e-9 * abs(best_L)) {
      best_L <- l$L
      best_lambda <- lambda
      stall <- 0L
    } else if ((stall <- stall + 1L) >= patience) {
      theta <- theta / 2
      stall <- 0L
    }
    if (all(l$g <= 1e-9) && sum(p$cost[l$best]) < 1.02 * ub) {
      cand <- .menu_step_down(p, l$best)
      if (sum(p$cost[cand]) < ub) {
        ub <- sum(p$cost[cand])
        pick <- cand
      }
    }
    if (theta < 1e-5 || ub - best_L <= 1e-9 * abs(ub)) break
    g <- l$g
    g[lambda <= 0 & g < 0] <- 0
    if (sum(g^2) == 0) break
    lambda <- pmax(0, lambda + theta * (min(ub, ceiling_cost) - l$L) /
                     sum(g^2) * g)
  }
  cand <- .menu_plan_from(p, best_lambda)
  if (!is.null(cand) && sum(p$cost[cand]) < ub) {
    ub <- sum(p$cost[cand])
    pick <- cand
  }
  list(lambda = best_lambda, L = best_L, ub = ub, pick = pick)
}

#' The exact plan, or the best one found within a node budget
#'
#' Groups that share no column with a nonzero load are separate problems,
#' each solved on its own. Inside one, the reduced-cost test at the root
#' multipliers keeps only states that can be in a plan cheaper than the
#' incumbent by more than `grain`, the step between attainable costs (0
#' for a continuous cost), and a state is dropped when another of its
#' group costs no more and loads no column more.
#'
#' The search is best-first branch and bound. Each node fixes some groups
#' to one state, re-solves the Lagrangian bound on the rest by a few
#' subgradient steps from its parent's multipliers, and drops states by
#' the same reduced-cost test. It branches on the group whose two best
#' reduced costs are closest, and solves the last one or two free groups
#' outright. Past `max_open` open nodes it continues depth-first from the
#' node with the smallest bound, so memory stays bounded. It stops after
#' `max_nodes` nodes, and the smallest bound left open is then a lower
#' bound on the optimum.
#'
#' Returns the plan (`pick`, one state per group), its `value`, the lower
#' `bound`, the `status` (`"optimal"`, `"limit"` or `"infeasible"`) and
#' the number of `nodes` and of states `kept` at the root. `pick` is
#' `NULL` when no plan was found.
#' @keywords internal
#' @noRd
.menu_search <- function(p, lambda, pick = NULL, grain = 0,
                         max_nodes = 5000, max_open = 1e4) {
  K <- ncol(p$C)
  if (K == 0L) {
    pick <- .menu_argmin(p$cost, p$groups)
    v <- sum(p$cost[pick])
    return(list(pick = pick, value = v, bound = v, status = "optimal",
                nodes = 0, kept = p$H))
  }
  inc <- matrix(FALSE, p$H, K)
  for (h in seq_len(p$H)) {
    inc[h, ] <- colSums(p$C[p$groups[[h]], , drop = FALSE] != 0) > 0
  }
  comp <- seq_len(p$H)
  repeat {
    old <- comp
    for (k in seq_len(K)) {
      h <- which(inc[, k])
      if (length(h)) comp[comp %in% comp[h]] <- min(comp[h])
    }
    if (identical(old, comp)) break
  }
  out <- list(pick = integer(p$H), value = 0, bound = 0, status = "optimal",
              nodes = 0, kept = 0)
  for (cc in unique(comp)) {
    S <- which(comp == cc)
    Kc <- which(colSums(inc[S, , drop = FALSE]) > 0)
    if (!length(Kc)) {
      s <- .menu_argmin(p$cost, p$groups[S])
      out$pick[S] <- s
      out$value <- out$value + sum(p$cost[s])
      out$bound <- out$bound + sum(p$cost[s])
      out$kept <- out$kept + length(S)
      next
    }
    idx <- unlist(p$groups[S], use.names = FALSE)
    q <- .menu_problem(p$cost[idx], p$C[idx, Kc, drop = FALSE],
                       match(p$owner[idx], S))
    start <- if (!is.null(pick)) match(pick[S], idx)
    r <- .menu_component(q, lambda[Kc], start, grain,
                         max(max_nodes - out$nodes, 1), max_open)
    out$nodes <- out$nodes + r$nodes
    out$kept <- out$kept + r$kept
    if (is.null(r$pick)) {
      return(list(pick = NULL, value = Inf, bound = r$bound,
                  status = r$status, nodes = out$nodes))
    }
    out$pick[S] <- idx[r$pick]
    out$value <- out$value + r$value
    out$bound <- out$bound + r$bound
    if (r$status == "limit") out$status <- "limit"
  }
  out
}

#' `.menu_search()` on one component
#' @keywords internal
#' @noRd
.menu_component <- function(q, lambda, start, grain, max_nodes, max_open) {
  K <- ncol(q$C)
  best <- if (is.null(start)) Inf else sum(q$cost[start])
  plan <- start
  ceiling_cost <- sum(vapply(q$groups, function(j) max(q$cost[j]), 0)) + 1
  # Relative to the incumbent, not the ceiling. A budget-mode objective near
  # 1e-4 sits far below the states at a tiny fraction.
  tol <- 1e-9 * max(abs(if (is.finite(best)) best else ceiling_cost),
                    1e-300)
  # A node is worth opening only if it can hold a plan cheaper than the
  # incumbent by a grain, or by the tolerance where cost is continuous.
  limit <- function() {
    if (is.finite(best)) best - max(grain - tol, tol) else ceiling_cost
  }
  improve <- function(cost, states) {
    if (cost < best - tol) {
      best <<- cost
      plan <<- states
    }
  }
  fits <- function(j, cap) {
    j[rowSums(q$C[j, , drop = FALSE] <= rep(cap, each = length(j)) + 1e-9) ==
        K]
  }

  l <- .menu_lagrangian(q, lambda)
  rmin <- l$r[l$best][q$owner]
  keep <- l$r - rmin <= limit() - l$L
  for (j in q$groups) {
    j <- j[keep[j]]
    if (length(j) < 2L) next
    # In cost order, so a kept state never costs more than the next.
    o <- j[order(q$cost[j], rowSums(q$C[j, , drop = FALSE]))]
    keep[o] <- .menu_dominated(q$C, o)
  }

  # Fixes groups down to one state, solves one or two free groups
  # outright, and otherwise re-solves the bound and narrows each group.
  # Returns NULL when the node holds no cheaper plan.
  evaluate <- function(nd) {
    repeat {
      one <- which(is.na(nd$plan) & lengths(nd$alive) == 1L)
      if (!length(one)) break
      for (h in one) {
        s <- nd$alive[[h]]
        nd$c_fix <- nd$c_fix + q$cost[s]
        nd$cap <- nd$cap - q$C[s, ]
        nd$plan[h] <- s
      }
    }
    if (any(nd$cap < -1e-9)) return(NULL)
    free <- which(is.na(nd$plan))
    if (!length(free)) {
      improve(nd$c_fix, nd$plan)
      return(NULL)
    }
    if (length(free) == 1L) {
      j <- fits(nd$alive[[free]], nd$cap)
      if (length(j)) {
        s <- j[which.min(q$cost[j])]
        nd$plan[free] <- s
        improve(nd$c_fix + q$cost[s], nd$plan)
      }
      return(NULL)
    }
    if (length(free) == 2L) {
      a <- free[order(lengths(nd$alive[free]))]
      j1 <- nd$alive[[a[1]]]
      j2 <- nd$alive[[a[2]]]
      j2 <- j2[order(q$cost[j2])]
      for (s in j1[order(q$cost[j1])]) {
        if (nd$c_fix + q$cost[s] + q$cost[j2[1]] >= best - tol) break
        j <- fits(j2, nd$cap - q$C[s, ])
        if (length(j)) {
          nd$plan[a] <- c(s, j[1])
          improve(nd$c_fix + q$cost[s] + q$cost[j[1]], nd$plan)
        }
      }
      return(NULL)
    }
    sizes <- lengths(nd$alive[free])
    ends <- as.integer(cumsum(sizes))
    grp <- rep.int(seq_along(free), sizes)
    j_all <- unlist(nd$alive[free], use.names = FALSE)
    Cj <- q$C[j_all, , drop = FALSE]
    cj <- q$cost[j_all]
    cmin <- .menu_group_min(Cj, ends)$min
    if (any(colSums(cmin) > nd$cap + 1e-9)) return(NULL)
    mu <- nd$mu
    best_L <- -Inf
    best_mu <- mu
    theta <- 1
    for (it in seq_len(30L)) {
      r <- cj + drop(Cj %*% mu)
      gm <- .menu_group_min(r, ends)
      arg <- gm$arg
      L <- nd$c_fix + sum(gm$min) - sum(mu * nd$cap)
      if (L > best_L) {
        best_L <- L
        best_mu <- mu
      } else {
        theta <- theta * 0.7
      }
      g <- colSums(Cj[arg, , drop = FALSE]) - nd$cap
      if (all(g <= 1e-9)) {
        states <- nd$plan
        states[free] <- j_all[arg]
        improve(nd$c_fix + sum(cj[arg]), states)
      }
      if (best_L > limit()) break
      g[mu <= 0 & g < 0] <- 0
      if (sum(g^2) == 0) break
      mu <- pmax(0, mu + theta * (limit() - L) / sum(g^2) * g)
    }
    slack <- limit() - best_L
    if (slack < 0) return(NULL)
    r <- cj + drop(Cj %*% best_mu)
    gm <- .menu_group_min(r, ends)
    k <- r - gm$min[grp] <= slack
    nd$alive[free] <- split(j_all[k], grp[k])
    gap <- ifelse(gm$second - gm$min <= slack, gm$second - gm$min, Inf)
    nd$L <- best_L
    nd$mu <- best_mu
    if (all(!is.finite(gap))) return(evaluate(nd))
    nd$branch <- free[which.min(gap)]
    nd
  }

  root <- if (all(vapply(q$groups, function(j) any(keep[j]), logical(1)))) {
    evaluate(list(alive = lapply(q$groups, function(j) j[keep[j]]),
                  plan = rep(NA_integer_, q$H), mu = lambda,
                  cap = rep(1, K), c_fix = 0, L = l$L))
  }
  open <- if (is.null(root)) list() else list(root)
  bounds <- if (is.null(root)) numeric(0) else root$L
  stack <- list()
  nodes <- 1
  while (length(open) || length(stack)) {
    if (nodes >= max_nodes) break
    if (length(stack)) {
      nd <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      if (nd$L > limit()) next
      deep <- TRUE
    } else {
      i <- which.min(bounds)
      if (bounds[i] > limit()) {
        open <- list()
        bounds <- numeric(0)
        break
      }
      nd <- open[[i]]
      open[[i]] <- NULL
      bounds <- bounds[-i]
      deep <- length(open) >= max_open
    }
    h <- nd$branch
    kids <- list()
    for (s in nd$alive[[h]]) {
      ch <- nd
      ch$alive[[h]] <- s
      nodes <- nodes + 1
      ch <- evaluate(ch)
      if (!is.null(ch) && ch$L <= limit()) kids[[length(kids) + 1L]] <- ch
    }
    if (!length(kids)) next
    kl <- vapply(kids, `[[`, 0, "L")
    if (deep) {
      stack <- c(stack, kids[order(kl, decreasing = TRUE)])
    } else {
      open <- c(open, kids)
      bounds <- c(bounds, kl)
    }
  }
  left <- c(bounds, vapply(stack, `[[`, 0, "L"))
  left <- left[left <= limit()]
  if (!is.finite(best)) {
    if (!length(left)) {
      return(list(pick = NULL, value = Inf, bound = Inf,
                  status = "infeasible", nodes = nodes, kept = sum(keep)))
    }
    return(list(pick = NULL, value = Inf, bound = min(left), status = "limit",
                nodes = nodes, kept = sum(keep)))
  }
  list(
    pick = plan,
    value = best,
    bound = if (length(left)) min(min(left), best) else best,
    status = if (length(left)) "limit" else "optimal",
    nodes = nodes,
    kept = sum(keep)
  )
}
