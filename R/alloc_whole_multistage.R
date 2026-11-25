#' Variance and cost terms of one stratum's whole take vectors
#'
#' A state's variance at `a` PSUs is `C / a - D`, and its PSU cost is `pc`.
#' @keywords internal
#' @noRd
.ms_state_terms <- function(prep, h, B, lo_i, hi_i, costed, w2s2d, r) {
  P <- apply(B, 1L, prod)
  rows <- .ms_stage_rows(prep, h, B)
  infl <- rows$inflation
  N_psu <- prep$N_psu_h[h] %||% Inf
  fpc <- prep$fpc %||% "unit"
  R1 <- prep$var_ratio_mat[h, 1L] * prep$icc_mat[h, 1L]
  N_fpc <- switch(
    fpc,
    unit = rep(.alloc_scope_N(prep)[h], length(P)),
    none = rep(Inf, length(P)),
    stage = if (R1 > 0) N_psu * infl / R1 else rep(Inf, length(P))
  )
  list(
    B = B,
    P = P,
    lo = pmax(1, ceiling(lo_i / P - 1e-9)),
    hi = pmin(floor(hi_i / P + 1e-9), N_psu,
              .psu_room(prep$N_h[h], prep$N_psu_h[h], P)),
    pc = if (costed) rows$psu_cost else P,
    C = w2s2d * infl / (P * r),
    D = w2s2d * infl / N_fpc
  )
}

#' `.cluster_stage_row()` for every row of a take matrix at once
#' @keywords internal
#' @noRd
.ms_stage_rows <- function(prep, h, B) {
  S <- prep$stages
  n <- nrow(B)
  storage.mode(B) <- "double"
  eff <- B * rep(prep$stage_resp_mat[h, -1L], each = n)
  Q <- matrix(1, n, S)
  for (s in (S - 1L):1L) Q[, s] <- Q[, s + 1L] * eff[, s]
  vr <- prep$var_ratio_mat[h, ]
  icc <- prep$icc_mat[h, ]
  R_h <- c(vr * icc, vr[S - 1L] * (1 - icc[S - 1L]))
  keep <- matrix(1, n, S - 1L)
  if (identical(prep$fpc %||% "unit", "stage")) {
    N_st <- prep$N_stages_mat[h, ]
    M <- c(N_st[-1L] / N_st[-(S - 1L)], .alloc_scope_N(prep)[h] / N_st[S - 1L])
    M <- M[seq_len(S - 1L)]
    keep <- 1 - pmin(1, eff / rep(M, each = n))
  }
  inflation <- R_h[1L] * Q[, 1L] +
    drop((keep * Q[, -1L, drop = FALSE]) %*% R_h[-1L])
  psu_cost <- if (isTRUE(prep$has_stage_costs)) {
    cum <- B
    if (S > 2L) for (s in 2L:(S - 1L)) cum[, s] <- cum[, s - 1L] * B[, s]
    prep$cost_mat[h, 1L] + drop(cum %*% prep$cost_mat[h, -1L])
  } else {
    rep(NA_real_, n)
  }
  list(inflation = inflation, psu_cost = psu_cost)
}

#' A lower bound on the cost of every design completing the takes `t[1..j]`
#' @keywords internal
#' @noRd
.ms_cost_floor <- function(cost, t, j, low, high, lo_i, hi_i) {
  S1 <- length(t)
  rest <- seq_len(S1) > j
  t_low <- ifelse(rest, low, t)
  t_high <- ifelse(rest, high, t)
  P_max <- min(prod(t_high), hi_i)
  a_min <- max(1, ceiling(lo_i / P_max - 1e-9))
  Q_low <- c(1, cumprod(t_low))
  suffix <- c(rev(cumprod(rev(t_high))), 1)
  per_stage <- pmax(a_min * Q_low, ifelse(is.finite(suffix), lo_i / suffix, 0))
  sum(cost * per_stage) * (1 - 1e-12)
}

#' Whole take vectors a stratum can field within a cost slack
#'
#' A take whose completions all cost more is skipped, not a stop, because a
#' larger take can need fewer PSUs to reach a size floor.
#' @keywords internal
#' @noRd
.ms_enumerate <- function(prep, h, base, free, floor_t, caps, lo_i, hi_i,
                          slack, costed, keep, limit) {
  S1 <- length(base)
  cost <- if (costed) prep$cost_mat[h, ] else NULL
  low <- ifelse(free, floor_t, base)
  high <- ifelse(free, caps, base)
  floor_cost <- if (costed) cost else c(rep(0, S1), 1)
  pc_of <- function(t) if (costed) sum(cost * c(1, cumprod(t))) else prod(t)
  last <- max(which(free))
  out <- list()
  n <- 0L
  block <- function(t) {
    t0 <- t
    t0[last] <- 0
    t1 <- t
    t1[last] <- 1
    A <- pc_of(t0)
    B <- pc_of(t1) - A
    P_other <- prod(t[-last])
    v_max <- min(caps[last], floor(hi_i / P_other + 1e-9),
                 floor((slack - A) / B + 1e-9))
    if (!(v_max >= low[last])) return(TRUE)
    from <- low[last]
    while (from <= v_max) {
      to <- min(v_max, from + 1e5 - 1)
      v <- seq.int(from, to)
      M <- matrix(t, length(v), S1, byrow = TRUE)
      M[, last] <- v
      M <- M[keep(M), , drop = FALSE]
      if (nrow(M)) {
        n <<- n + nrow(M)
        if (n > limit) return(FALSE)
        out[[length(out) + 1L]] <<- M
      }
      from <- to + 1
    }
    TRUE
  }
  walk <- function(t, j) {
    if (j == last) return(block(t))
    if (!free[j]) return(walk(t, j + 1L))
    v <- low[j]
    repeat {
      if (v > caps[j]) break
      t[j] <- v
      tt <- t
      if (j < S1) tt[(j + 1L):S1] <- low[(j + 1L):S1]
      P <- prod(tt)
      if (P > hi_i) break
      if (pc_of(tt) > slack) break
      if (.ms_cost_floor(floor_cost, t, j, low, high, lo_i, hi_i) <= slack) {
        if (!walk(t, j + 1L)) return(FALSE)
      }
      v <- v + 1
    }
    TRUE
  }
  if (!walk(low, 1L)) return(NULL)
  if (!length(out)) return(matrix(numeric(0), 0L, S1))
  do.call(rbind, out)
}

#' The fewest PSUs a state needs to leave `room` for its variance
#' @keywords internal
#' @noRd
.ms_need <- function(t, room) {
  need <- ifelse(t$C <= 0, t$lo, ceiling(t$C / pmax(room + t$D, 0) - 1e-9))
  need[!is.finite(need)] <- Inf
  pmax(t$lo, need)
}

#' The take vectors on a stratum's cost-variance frontier
#'
#' A point another matches or beats on cost and variance is never needed.
#' @keywords internal
#' @noRd
.ms_frontier <- function(t, a_lo, a_hi) {
  n <- length(t$P)
  len <- pmax(0, a_hi - a_lo + 1)
  total <- sum(len)
  if (total == 0) return(integer(0))
  if (total > 2e6) return(NULL)
  id <- rep(seq_len(n), len)
  a <- a_lo[id] + sequence(len) - 1
  cost <- a * t$pc[id]
  var <- t$C[id] / a - t$D[id]
  ord <- order(cost, var)
  v <- var[ord]
  before <- c(Inf, cummin(v)[-length(v)])
  sort(unique(id[ord[v < before]]))
}

#' Exact whole design of a stratified allocation with three or more stages
#'
#' Status `"heuristic"` with a `reason` when it cannot search, `"infeasible"` when none fits.
#' A budget within the relative `1e-8` that `n_alloc()` accepts below the
#' cheapest whole design buys that design: the search is repeated at its
#' cost.
#' @keywords internal
#' @noRd
.alloc_whole_cluster_multistage <- function(prep, alloc, mode, cv, budget,
                                            ...) {
  w <- .alloc_whole_cluster_multistage_at(prep, alloc, mode, cv, budget, ...)
  if (identical(mode, "budget") && identical(w$status, "infeasible") &&
        isTRUE(w$min_cost <= budget * (1 + 1e-8))) {
    w <- .alloc_whole_cluster_multistage_at(prep, alloc, mode, cv,
                                            w$min_cost, ...)
  }
  w
}

#' The whole design search at one budget
#' @keywords internal
#' @noRd
.alloc_whole_cluster_multistage_at <- function(prep, alloc, mode, cv, budget,
                                               deff, resp_rate, lo_i, hi_i,
                                               b_rule, max_nodes, gap = 0,
                                               target_se = NULL,
                                               incumbent = NULL,
                                               state_limit = 2000L) {
  H <- length(prep$N_h)
  W <- .alloc_scope_N(prep) / sum(.alloc_scope_N(prep))
  r <- rep_len(resp_rate, H)
  w2s2d <- W^2 * prep$S_raw_h^2 * rep_len(deff, H)
  costed <- alloc == "optimal" && isTRUE(prep$has_stage_costs)
  free_mat <- costed & !prep$takes_fixed_mat
  floor_mat <- .cluster_take_floor(prep)
  se_goal <- if (mode == "cv") {
    if (!is.null(target_se)) {
      target_se
    } else {
      cv * abs(.aggregate_mean(W, prep$mean_h))
    }
  }
  room_total <- if (mode == "cv") se_goal^2 * (1 + 1e-12)
  cap <- if (mode == "cv") (incumbent$cost %||% Inf) * (1 + 1e-9) else
    .budget_room(budget)

  empty <- function() {
    if (is.null(incumbent)) return(refused())
    list(status = "heuristic", reason = "pruning", a = incumbent$a,
         b = incumbent$b)
  }
  ## A budget refusal names the cheapest whole design. Found apart from the
  ## pruning, it also checks the refusal: a cheapest design within the cap
  ## means the pruning lost a design, which is not a proof of infeasibility.
  refused <- function() {
    if (mode != "budget") return(list(status = "infeasible", a = NULL))
    least_cost <- .ms_cheapest_cost(
      prep, least, terms_of, free_mat, floor_mat, b_rule, lo_i, hi_i, costed
    )
    if (!is.na(least_cost) && least_cost <= cap) {
      return(list(status = "heuristic", reason = "pruning", a = NULL))
    }
    list(status = "infeasible", a = NULL, min_cost = least_cost)
  }

  terms_of <- function(h, B) {
    .ms_state_terms(prep, h, B, lo_i[h], hi_i[h], costed, w2s2d[h], r[h])
  }
  least <- vapply(seq_len(H), function(h) {
    low <- ifelse(free_mat[h, ], floor_mat[h, ], b_rule[h, ])
    high <- ifelse(free_mat[h, ],
                   .whole_caps(.cluster_take_caps(prep, h)), b_rule[h, ])
    unit_cost <- if (costed) prep$cost_mat[h, ] else c(rep(0, length(low)), 1)
    .ms_cost_floor(unit_cost, low, 0L, low, high, lo_i[h], hi_i[h])
  }, 1)
  vmin <- rep(0, H)
  states <- vector("list", H)
  for (pass in seq_len(20L)) {
    changed <- FALSE
    for (h in seq_len(H)) {
      slack <- cap - sum(least[-h])
      room <- if (mode == "cv") room_total - sum(vmin[-h])
      fits <- function(B) {
        t <- terms_of(h, B)
        ok <- t$lo <= t$hi
        if (mode == "cv") {
          need <- .ms_need(t, room)
          ok & need <= t$hi & need * t$pc <= slack
        } else {
          ok & t$lo * t$pc <= slack
        }
      }
      if (!any(free_mat[h, ])) {
        B <- matrix(b_rule[h, ], nrow = 1L)
        B <- B[fits(B), , drop = FALSE]
      } else {
        caps <- .whole_caps(.cluster_take_caps(prep, h))
        B <- .ms_enumerate(prep, h, b_rule[h, ], free_mat[h, ], floor_mat[h, ],
                           caps, lo_i[h], hi_i[h], slack, costed, fits, 1e6)
        if (is.null(B)) return(list(status = "heuristic", reason = "listing"))
      }
      t <- terms_of(h, B)
      if (!length(t$P)) return(empty())
      new_least <- if (mode == "cv") {
        min(.ms_need(t, room) * t$pc)
      } else {
        min(t$lo * t$pc)
      }
      new_vmin <- min(t$C / t$hi - t$D)
      if (new_least > least[h] * (1 + 1e-12) ||
          new_vmin > vmin[h] * (1 + 1e-12)) {
        changed <- TRUE
      }
      least[h] <- max(least[h], new_least)
      vmin[h] <- max(vmin[h], new_vmin)
      states[[h]] <- t
    }
    if (!changed) break
  }
  a_range <- function(h, t) {
    slack <- cap - sum(least[-h])
    a_lo <- if (mode == "cv") .ms_need(t, room_total - sum(vmin[-h])) else t$lo
    list(lo = a_lo, hi = pmin(t$hi, floor(slack / t$pc + 1e-9)))
  }
  for (h in seq_len(H)) {
    t <- states[[h]]
    ar <- a_range(h, t)
    keep <- .ms_frontier(t, ar$lo, ar$hi)
    if (is.null(keep) || length(keep) > state_limit) {
      return(list(status = "heuristic", reason = "frontier"))
    }
    if (!length(keep)) return(empty())
    t <- lapply(t, function(v) {
      if (is.matrix(v)) v[keep, , drop = FALSE] else v[keep]
    })
    t$lo <- ar$lo[keep]
    t$hi <- ar$hi[keep]
    states[[h]] <- t
  }
  result <- .ms_search(states, mode, se_goal, budget, prep, costed, max_nodes,
                       gap, incumbent)
  if (!is.null(incumbent) && is.null(result$a) &&
      result$status %in% c("infeasible", "limit")) {
    return(empty())
  }
  if (identical(result$status, "infeasible") && is.null(result$a)) {
    return(refused())
  }
  result
}

#' The cost of the cheapest whole design within the bounds, or NA
#'
#' The strata share only the budget, so the cheapest design takes each
#' stratum's cheapest. A stratum's is listed under a cost cap that doubles
#' from the pruning's lower bound `least[h]`. The listing keeps every take
#' vector whose designs can cost less than the cap, so the first cap that
#' lists one holds the cheapest. NA when a listing is too long or no cap
#' finds a design.
#' @keywords internal
#' @noRd
.ms_cheapest_cost <- function(prep, least, terms_of, free_mat, floor_mat,
                              b_rule, lo_i, hi_i, costed) {
  total <- 0
  for (h in seq_len(length(least))) {
    found <- NA_real_
    cap_h <- max(least[h], 1e-12)
    for (k in seq_len(128L)) {
      fits <- function(B) {
        t <- terms_of(h, B)
        t$lo <= t$hi & t$lo * t$pc <= cap_h
      }
      B <- if (!any(free_mat[h, ])) {
        matrix(b_rule[h, ], nrow = 1L)
      } else {
        .ms_enumerate(prep, h, b_rule[h, ], free_mat[h, ], floor_mat[h, ],
                      .whole_caps(.cluster_take_caps(prep, h)), lo_i[h],
                      hi_i[h], cap_h, costed, fits, 1e6)
      }
      if (is.null(B)) return(NA_real_)
      B <- B[fits(B), , drop = FALSE]
      if (nrow(B)) {
        t <- terms_of(h, B)
        found <- min(t$lo * t$pc)
        break
      }
      cap_h <- cap_h * 2
    }
    if (is.na(found)) return(NA_real_)
    total <- total + found
  }
  total
}

#' The cheapest design for a target at the variance floor
#' @keywords internal
#' @noRd
.ms_floor_design <- function(states, owner, vmin, se_goal, pc, C, D, lo, hi,
                             max_nodes, state_limit = 2000L) {
  H <- length(vmin)
  allow <- se_goal^2 * (1 + 1e-12) - sum(vmin)
  if (allow < 0) return(list(status = "infeasible", a = NULL))
  base <- C / hi - D - vmin[owner]
  front <- list(excess = 0, cost = 0, pick = matrix(0L, 1L, 0L),
                a = matrix(0, 1L, 0L))
  for (h in seq_len(H)) {
    rows <- list()
    for (j in which(owner == h & base <= allow)) {
      a_lo <- if (C[j] > 0) {
        max(lo[j], ceiling(C[j] * hi[j] / (C[j] + (allow - base[j]) * hi[j]) -
                             1e-9))
      } else {
        lo[j]
      }
      a_top <- if (C[j] > 0) hi[j] else lo[j]
      if (a_top - a_lo + 1 > state_limit) {
        return(list(status = "heuristic", reason = "numerical"))
      }
      if (a_lo > a_top) next
      a <- seq(a_lo, a_top)
      rows[[length(rows) + 1L]] <- data.frame(
        j = j, a = a, cost = a * pc[j],
        excess = base[j] + C[j] * (hi[j] - a) / (a * hi[j])
      )
    }
    if (!length(rows)) return(list(status = "infeasible", a = NULL))
    cand <- do.call(rbind, rows)
    cand <- cand[cand$excess <= allow, , drop = FALSE]
    pair <- expand.grid(f = seq_along(front$cost), k = seq_len(nrow(cand)))
    ex <- front$excess[pair$f] + cand$excess[pair$k]
    keep <- ex <= allow
    if (!any(keep)) return(list(status = "infeasible", a = NULL))
    pair <- pair[keep, , drop = FALSE]
    ex <- ex[keep]
    cost <- front$cost[pair$f] + cand$cost[pair$k]
    ord <- order(cost, ex)
    nd <- ord[ex[ord] < c(Inf, cummin(ex[ord])[-length(ord)])]
    if (length(nd) > state_limit) {
      return(list(status = "heuristic", reason = "numerical"))
    }
    front <- list(
      excess = ex[nd], cost = cost[nd],
      pick = cbind(front$pick[pair$f[nd], , drop = FALSE], cand$j[pair$k[nd]]),
      a = cbind(front$a[pair$f[nd], , drop = FALSE], cand$a[pair$k[nd]])
    )
  }
  best <- which.min(front$cost)
  pick <- front$pick[best, ]
  B <- do.call(rbind, lapply(states, `[[`, "B"))
  cost <- front$cost[best]
  list(
    a = front$a[best, ], b = B[pick, , drop = FALSE], status = "optimal",
    search = list(family = "whole units", status = "optimal", bound = cost,
                  gap = 0, nodes = 0, max_nodes = max_nodes,
                  states = length(owner), kept = length(front$cost))
  )
}

#' Search the whole design over the frontier states of every stratum
#'
#' The problem shapes are those of `.alloc_whole_cluster()`.
#' @keywords internal
#' @noRd
.ms_search <- function(states, mode, se_goal, budget, prep, costed,
                       max_nodes, gap, incumbent) {
  H <- length(states)
  owner <- rep(seq_len(H), vapply(states, function(t) length(t$P), 1L))
  col <- function(name) unlist(lapply(states, `[[`, name), use.names = FALSE)
  pc <- col("pc")
  C <- col("C")
  D <- col("D")
  lo <- col("lo")
  hi <- col("hi")
  if (mode == "cv") {
    v_at_hi <- C / hi - D
    lightest <- vapply(seq_len(H), function(h) {
      i <- which(owner == h)
      i[which.min(v_at_hi[i])]
    }, 1L)
    vmin <- v_at_hi[lightest]
    room <- se_goal^2 - sum(vmin)
    if (room <= 1e-12 * se_goal^2) {
      return(.ms_floor_design(states, owner, vmin, se_goal, pc, C, D, lo, hi,
                              max_nodes))
    }
    m <- .whole_target_scale() / room
    g2 <- C * m
    g0 <- (-D - vmin[owner]) * m
    g0[lightest] <- -(g2[lightest] / hi[lightest])
    q <- .whole_problem(pc, 0, 0, g2, lo, hi, owner = owner, o0 = 0, g0 = g0)
    grain <- if (costed) .menu_grain(as.vector(prep$cost_mat)) else 1
    offset <- 0
  } else {
    room <- .budget_room(budget)
    Dmax <- vapply(seq_len(H), function(h) max(D[owner == h]), 1)
    q <- .whole_problem(0, C, .whole_budget_load(pc, room), 0, lo, hi,
                        owner = owner,
                        o0 = Dmax[owner] - D, g0 = 0)
    grain <- 0
    offset <- -sum(Dmax)
  }
  start <- NULL
  if (!is.null(incumbent$a)) {
    pick <- integer(H)
    a0 <- numeric(H)
    for (h in seq_len(H)) {
      i <- which(owner == h)
      cost_h <- incumbent$a[h] * incumbent$pc[h]
      var_h <- incumbent$C[h] / incumbent$a[h] - incumbent$D[h]
      need <- .ms_need(list(C = C[i], D = D[i], lo = lo[i]), var_h)
      fits <- need <= hi[i] & need * pc[i] <= cost_h * (1 + 1e-12) &
        C[i] / need - D[i] <= var_h + 1e-12 * abs(var_h)
      best <- if (any(fits)) {
        k <- which(fits)[which.min((need * pc[i])[fits])]
        c(i[k], need[k])
      }
      if (is.null(best)) {
        start <- NULL
        break
      }
      pick[h] <- best[1L]
      a0[h] <- best[2L]
      start <- TRUE
    }
    if (isTRUE(start)) start <- .whole_plan(q, a0, branch = pick)
  }
  r_s <- .whole_search(q, start = start, grain = grain, max_nodes = max_nodes,
                       gap = gap)
  if (is.null(r_s$x)) return(list(status = r_s$status, a = NULL))
  B <- do.call(rbind, lapply(states, `[[`, "B"))
  list(
    a = as.numeric(r_s$x),
    b = B[r_s$branch, , drop = FALSE],
    status = r_s$status,
    search = .whole_record(r_s, "whole units", max_nodes, offset = offset)
  )
}
