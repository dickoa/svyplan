#' A whole-unit problem with one convex term per branch
#'
#' Each stratum takes one of its branches and a whole `x` in that branch's
#' `lo..hi`, with `lo >= 1`. A branch is a fixed choice that shapes the
#' stratum's terms, such as a cluster take. On branch `i` the cost is
#' `o0 + o1 x + o2 / x` and the load on column `k` is
#' `g0[, k] + g1[, k] x + g2[, k] / x`, with `o1`, `o2`, `g1` and `g2`
#' nonnegative, and a plan keeps every column's load at most 1. Under
#' multipliers `lambda` a state's reduced cost is
#' `gamma + alpha x + beta / x`, convex in `x` on each branch, so the
#' Lagrangian bound, its subgradient and the sizes reduced-cost fixing
#' keeps, one interval per branch, need no listing. With one branch per
#' stratum (`owner = seq_along(lo)`) a stratum's choice is its size alone.
#' Branches are stored sorted by stratum, and `orig` maps them back.
#' @keywords internal
#' @noRd
.whole_problem <- function(o1, o2, g1, g2, lo, hi, owner = seq_along(lo),
                           o0 = 0, g0 = 0) {
  B <- length(lo)
  K <- ncol(matrix(as.numeric(g1), B))
  ord <- order(owner, seq_len(B))
  pick <- function(v) rep_len(as.numeric(v), B)[ord]
  mat <- function(m) matrix(as.numeric(m), B, K)[ord, , drop = FALSE]
  owner <- as.integer(owner)[ord]
  list(
    o0 = pick(o0), o1 = pick(o1), o2 = pick(o2),
    g0 = mat(g0), g1 = mat(g1), g2 = mat(g2),
    lo = as.numeric(lo)[ord], hi = as.numeric(hi)[ord],
    owner = owner, ends = as.integer(cumsum(tabulate(owner))),
    orig = ord, pos = order(ord), H = max(owner), B = B, K = K
  )
}

.whole_abc <- function(q, lambda) {
  list(alpha = q$o1 + drop(q$g1 %*% lambda),
       beta = q$o2 + drop(q$g2 %*% lambda),
       gamma = q$o0 + drop(q$g0 %*% lambda))
}

.whole_cost <- function(q, plan) {
  b <- plan$b
  sum(q$o0[b] + q$o1[b] * plan$x + q$o2[b] / plan$x)
}

.whole_load <- function(q, plan) {
  b <- plan$b
  colSums(q$g0[b, , drop = FALSE] + q$g1[b, , drop = FALSE] * plan$x +
            q$g2[b, , drop = FALSE] / plan$x)
}

#' Each branch's best whole size under `alpha` and `beta`
#' @keywords internal
#' @noRd
.whole_best <- function(alpha, beta, lo, hi) {
  star <- ifelse(alpha > 0, sqrt(beta / alpha), ifelse(beta > 0, hi, lo))
  x1 <- pmin(pmax(floor(star), lo), hi)
  x2 <- pmin(pmax(ceiling(star), lo), hi)
  ifelse(alpha * x1 + beta / x1 <= alpha * x2 + beta / x2, x1, x2)
}

#' Lagrangian bound over every branch and whole size, its plan and
#' subgradient
#' @keywords internal
#' @noRd
.whole_lagrangian <- function(q, lambda) {
  abc <- .whole_abc(q, lambda)
  x <- .whole_best(abc$alpha, abc$beta, q$lo, q$hi)
  r <- abc$gamma + abc$alpha * x + abc$beta / x
  gm <- .menu_group_min(r, q$ends)
  plan <- list(b = gm$arg, x = x[gm$arg])
  list(L = sum(gm$min) - sum(lambda), plan = plan, rmin = gm$min,
       g = .whole_load(q, plan) - 1)
}

.whole_feasible <- function(q, plan) all(.whole_load(q, plan) <= 1 + 1e-9)

#' Factor that puts a target at the search's acceptance line
#'
#' `.whole_feasible()` accepts a load up to `1 + 1e-9`. A load measured as a
#' share of a target and multiplied by this factor is accepted when it is at
#' most `1 + 1e-12` of the target, the rounding at which a target counts as
#' met. A plan whose precision equals the target, such as one whose cv was
#' read back from a `prec_*()` call, is then one the search can return.
#' @keywords internal
#' @noRd
.whole_target_scale <- function() (1 + 1e-9) / (1 + 1e-12)

#' A cost's load against a monetary cap, accepted when the cost is at most the cap
#' @keywords internal
#' @noRd
.whole_budget_load <- function(cost, cap) cost * (1 + 1e-9) / cap

#' A plan from multipliers, doubled until the best states fit, or `NULL`
#' @keywords internal
#' @noRd
.whole_plan_from <- function(q, lambda) {
  if (all(lambda == 0)) lambda <- rep(1, q$K)
  for (i in seq_len(200L)) {
    plan <- .whole_lagrangian(q, lambda)$plan
    if (.whole_feasible(q, plan)) return(plan)
    lambda <- lambda * 2
  }
  NULL
}

#' Multipliers by subgradient, keeping the best plan met on the way
#'
#' `start`, a plan meeting every column, is kept when no cheaper plan is
#' found. With no plan, from `start` or from the multipliers, `pick` is
#' `NULL`. With one column the dual is maximized exactly by `.dual_one()`.
#' @keywords internal
#' @noRd
.whole_subgradient <- function(q, start = NULL, iters = 1500L, lambda = NULL) {
  at_lo <- abs(sum(q$o0 + q$o1 * q$lo + q$o2 / q$lo)) * q$H / q$B
  if (q$K == 1L) {
    one <- .whole_subgradient_one(q, start, at_lo)
    if (!is.null(one)) return(one)
  }
  warm <- !is.null(lambda) && length(lambda) == q$K &&
    all(is.finite(lambda)) && all(lambda >= 0)
  if (!warm) {
    lambda <- rep(at_lo / q$K + 1e-12, q$K)
  }
  pick <- .whole_plan_from(q, lambda)
  if (!is.null(start) &&
        (is.null(pick) || .whole_cost(q, start) < .whole_cost(q, pick))) {
    pick <- start
  }
  # With no plan yet the steps aim at a cost every plan stays below.
  ceiling_cost <- .whole_ceiling_cost(q)
  if (is.null(pick) && !is.finite(ceiling_cost)) {
    return(list(lambda = lambda, L = -Inf, pick = NULL))
  }
  ub <- if (is.null(pick)) Inf else .whole_cost(q, pick)
  best_L <- -Inf
  best_lambda <- lambda
  theta <- if (warm) 0.5 else 2
  patience <- if (warm) 15L else 25L
  stall <- 0L
  for (i in seq_len(iters)) {
    l <- .whole_lagrangian(q, lambda)
    if (!is.finite(best_L) || l$L > best_L + 1e-12 * abs(best_L)) {
      best_L <- l$L
      best_lambda <- lambda
      stall <- 0L
    } else if ((stall <- stall + 1L) >= patience) {
      theta <- theta / 2
      stall <- 0L
    }
    if (all(l$g <= 1e-9) && .whole_cost(q, l$plan) < ub) {
      ub <- .whole_cost(q, l$plan)
      pick <- l$plan
    }
    if (theta < 1e-6 || (is.finite(ub) && ub - best_L <= 1e-12 * abs(ub))) {
      break
    }
    g <- l$g
    g[lambda <= 0 & g < 0] <- 0
    if (sum(g^2) == 0) break
    lambda <- pmax(0, lambda + theta * (min(ub, ceiling_cost) - l$L) /
                     sum(g^2) * g)
  }
  cand <- .whole_plan_from(q, best_lambda)
  if (!is.null(cand) && .whole_cost(q, cand) < ub) {
    ub <- .whole_cost(q, cand)
    pick <- cand
  }
  if (is.null(pick)) return(list(lambda = best_lambda, L = best_L, pick = NULL))
  list(lambda = best_lambda, L = best_L, ub = ub, pick = pick)
}

#' A cost above every plan's, `Inf` when a branch is unbounded
#'
#' Each branch's cost is convex in `x`, so its largest value is at an end.
#' @keywords internal
#' @noRd
.whole_ceiling_cost <- function(q) {
  at <- function(x) q$o0 + q$o1 * x + q$o2 / x
  top <- pmax(at(q$lo), at(q$hi))
  if (any(!is.finite(top))) return(Inf)
  -sum(.menu_group_min(-top, q$ends)$min) + abs(sum(top)) + 1
}

#' The smaller root of `a x^2 - top x + b`, without cancellation
#' @keywords internal
#' @noRd
.whole_small_root <- function(a, b, top, d) {
  ifelse(top + d > 0, 2 * b / (top + d), (top - d) / (2 * a))
}

#' `.whole_subgradient()` with one column, or `NULL`
#' @keywords internal
#' @noRd
.whole_subgradient_one <- function(q, start, at_lo) {
  pick <- start
  ub <- if (is.null(pick)) Inf else .whole_cost(q, pick)
  f <- function(lambda) {
    l <- .whole_lagrangian(q, lambda)
    if (l$g <= 1e-9) {
      v <- .whole_cost(q, l$plan)
      if (v < ub) {
        ub <<- v
        pick <<- l$plan
      }
    }
    list(L = l$L, g = l$g)
  }
  d <- .dual_one(f, at_lo + 1e-12)
  if (is.null(d)) return(NULL)
  cand <- .whole_plan_from(q, d$lambda)
  if (!is.null(cand) && .whole_cost(q, cand) < ub) {
    ub <- .whole_cost(q, cand)
    pick <- cand
  }
  if (is.null(pick)) return(list(lambda = d$lambda, L = -Inf, pick = NULL))
  list(lambda = d$lambda, L = d$L, ub = ub, pick = pick)
}

#' The sizes reduced-cost fixing keeps, one interval per branch
#'
#' A state whose reduced cost exceeds its stratum's minimum by more than
#' `slack` is in no plan cheaper than the bound plus `slack`. On a branch
#' the kept sizes lie between the roots of `alpha x^2 - T x + beta`, with
#' `T` the stratum's minimum plus `slack` less the branch's constant. A
#' branch with no kept size comes back with `lo > hi`.
#' @keywords internal
#' @noRd
.whole_intervals <- function(q, lambda, slack) {
  abc <- .whole_abc(q, lambda)
  rmin <- .whole_lagrangian(q, lambda)$rmin
  top <- rmin[q$owner] + slack - abc$gamma
  a <- abc$alpha
  b <- abc$beta
  lo <- q$lo
  hi <- q$hi
  both <- a > 0 & b > 0
  d <- sqrt(pmax(top^2 - 4 * a * b, 0))
  lo[both] <- pmax(lo, ceiling(.whole_small_root(a, b, top, d) - 1e-9))[both]
  hi[both] <- pmin(hi, floor((top + d) / (2 * a) + 1e-9))[both]
  only_a <- a > 0 & !(b > 0)
  hi[only_a] <- pmin(hi, floor(top / a + 1e-9))[only_a]
  only_b <- !(a > 0) & b > 0
  lo[only_b] <- pmax(lo, ceiling(b / top - 1e-9))[only_b]
  # A branch whose best state already passes the threshold keeps nothing.
  none <- top < 0 | (both & top^2 < 4 * a * b - 1e-12 * abs(top)^2)
  lo[none] <- Inf
  list(lo = lo, hi = hi)
}

#' A plan given by stratum sizes, one branch per stratum, or by
#' `branch`, branches as the caller numbered them
#' @keywords internal
#' @noRd
.whole_plan <- function(q, x, branch = NULL) {
  if (is.null(branch)) branch <- q$orig[match(seq_len(q$H), q$owner)]
  list(b = q$pos[branch], x = as.numeric(x))
}

#' The best whole plan, by multipliers, fixing and the menu search
#'
#' The subgradient gives multipliers and a plan, the plan's cost bounds
#' what fixing keeps, and `.menu_search()` runs on the kept intervals. A
#' state left out is in no plan cheaper than the plan by `grain`, so a
#' proof on the intervals is a proof over every branch and whole size.
#' `start` is a plan from `.whole_plan()`. Returns each stratum's size `x`
#' and `branch` (as the caller numbered them), the plan's `value`, a lower
#' `bound` and the search's `status` (`"optimal"`, `"limit"` or
#' `"infeasible"`), with the nodes and the states listed and kept. With a
#' plan it also returns the multipliers `lambda` and each stratum's `top`,
#' the reduced cost above which a state is in no better plan, so a caller
#' that listed only some branches can prove the others out.
#' @keywords internal
#' @noRd
.whole_search <- function(q, start = NULL, grain = 0, max_nodes = 5000,
                          gap = 0, lambda = NULL) {
  if (is.numeric(start)) start <- .whole_plan(q, start)
  none <- list(x = NULL, branch = NULL, value = Inf, bound = Inf,
               status = "infeasible", nodes = 0, states = 0, kept = 0)
  # When every load falls as x grows, each stratum's smallest load per
  # column, at its largest size, bounds what any plan can reach.
  if (all(q$g1 == 0)) {
    floor_load <- .menu_group_min(q$g0 + q$g2 / q$hi, q$ends)$min
    if (any(colSums(matrix(floor_load, q$H)) > 1 + 1e-9)) return(none)
  }
  sg <- .whole_subgradient(q, start, lambda = lambda)
  if (is.null(sg$pick)) {
    none$status <- "limit"
    none$bound <- sg$L
    return(none)
  }
  slack <- max(sg$ub - grain - sg$L, 0) + 1e-9 * abs(sg$ub)
  iv <- .whole_intervals(q, sg$lambda, slack)
  pb <- sg$pick$b
  iv$lo[pb] <- pmin(iv$lo[pb], sg$pick$x)
  iv$hi[pb] <- pmax(iv$hi[pb], sg$pick$x)
  keep <- which(iv$lo <= iv$hi)
  xs <- lapply(keep, function(i) seq(iv$lo[i], iv$hi[i]))
  branch <- rep(keep, lengths(xs))
  x <- unlist(xs, use.names = FALSE)
  p <- .menu_problem(
    q$o0[branch] + q$o1[branch] * x + q$o2[branch] / x,
    q$g0[branch, , drop = FALSE] + q$g1[branch, , drop = FALSE] * x +
      q$g2[branch, , drop = FALSE] / x,
    q$owner[branch]
  )
  at <- vapply(seq_len(q$H), function(h) {
    j <- p$groups[[h]]
    j[branch[j] == pb[h] & x[j] == sg$pick$x[h]][1]
  }, 1L)
  sg2 <- .menu_subgradient(p, sg$lambda, start = at, exact = TRUE)
  sr <- .menu_search(p, sg2$lambda, sg2$pick, grain = grain,
                     max_nodes = max_nodes, gap = gap, exact = TRUE)
  if (is.null(sr$pick)) {
    stop("internal error: the whole-unit search lost the plan it started from",
         call. = FALSE)
  }
  list(x = x[sr$pick], branch = q$orig[branch[sr$pick]], value = sr$value,
       bound = sr$bound, status = sr$status, nodes = sr$nodes,
       states = length(x), kept = sr$kept, lambda = sg$lambda,
       top = .whole_lagrangian(q, sg$lambda)$rmin + slack)
}

#' The search record a fit carries in `$optimization$search`
#'
#' `offset` is a constant added to every plan's value, such as the part of
#' an objective no size changes.
#' @keywords internal
#' @noRd
.whole_record <- function(r, family, max_nodes, offset = 0) {
  value <- r$value + offset
  bound <- r$bound + offset
  list(
    family = family,
    status = r$status,
    bound = bound,
    gap = if (value > 0) (value - bound) / value else 0,
    nodes = r$nodes,
    max_nodes = max_nodes,
    states = r$states,
    kept = r$kept
  )
}

#' The status and record that speak for a whole design searched by domain
#' @keywords internal
#' @noRd
.whole_domains_summary <- function(records) {
  status <- vapply(records, function(r) r$status, character(1))
  measured <- Filter(function(r) is.numeric(r$gap), records)
  worst <- if (length(measured)) {
    measured[[which.max(vapply(measured, function(r) r$gap, numeric(1)))]]
  } else {
    records[[1L]]
  }
  list(
    status = if (any(status == "heuristic")) {
      "heuristic"
    } else if (all(status == "optimal")) {
      "optimal"
    } else {
      worst$status
    },
    search = worst,
    domains = records
  )
}

#' Mark the searched domains whose design a floor repair changed
#' @keywords internal
#' @noRd
.whole_repaired <- function(records, domain_idx, before, after) {
  for (d in seq_along(domain_idx)) {
    nm <- names(domain_idx)[d]
    idx <- domain_idx[[d]]
    if (!identical(records[[nm]]$status, "heuristic") &&
          any(before[idx] != after[idx])) {
      records[[nm]] <- list(status = "heuristic", reason = "floor repair",
                            search = records[[nm]])
    }
  }
  records
}

#' `as.integer()` of a whole count, refused past the integer range
#' @keywords internal
#' @noRd
.as_count <- function(x, what, where) {
  if (any(x > .Machine$integer.max, na.rm = TRUE)) .stop_count_range(what, where)
  as.integer(x)
}

#' The integer form of a result's count, refused when it has none
#' @keywords internal
#' @noRd
.coerce_count <- function(x) {
  if (any(x > .Machine$integer.max, na.rm = TRUE)) {
    stop(structure(
      list(message = sprintf(
        "the count %s is past the integer range (%s), so it has no integer form. Read it as a double, from $operational$n where the result has a field design",
        format(max(x, na.rm = TRUE), big.mark = ",", scientific = FALSE),
        format(.Machine$integer.max, big.mark = ",")
      ), call = NULL),
      class = c("svyplan_count_range", "error", "condition")
    ))
  }
  as.integer(x)
}

#' The search settings of `control`, with their defaults filled in
#' @keywords internal
#' @noRd
.search_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control) || is.data.frame(control) ||
        (length(control) > 0L &&
           (is.null(names(control)) || any(!nzchar(names(control)))))) {
    stop("'control' must be a named list", call. = FALSE)
  }
  unknown <- setdiff(names(control), c("max_nodes", "gap"))
  if (length(unknown) > 0L) {
    stop(sprintf(
      "unknown 'control' setting%s %s, the known ones are 'max_nodes' and 'gap'",
      if (length(unknown) == 1L) "" else "s",
      paste0("'", unknown, "'", collapse = ", ")
    ), call. = FALSE)
  }
  nodes <- control[["max_nodes"]] %||% 5000
  if (!is.numeric(nodes) || length(nodes) != 1L || !is.finite(nodes) ||
        nodes < 1 || abs(nodes - round(nodes)) > 1e-8) {
    stop("'control$max_nodes' must be one positive whole number",
         call. = FALSE)
  }
  gap <- control[["gap"]] %||% 0
  if (!is.numeric(gap) || length(gap) != 1L || !is.finite(gap) || gap < 0 ||
        gap >= 1) {
    stop("'control$gap' must be one number in [0, 1)", call. = FALSE)
  }
  list(max_nodes = round(nodes), gap = gap)
}

#' Whole counts kept as doubles, exact up to 2^53
#' @keywords internal
#' @noRd
.whole_storage <- function(x) {
  if (any(abs(x) > 2^53)) {
    stop("a whole count exceeds 2^53, where doubles stop being exact",
         call. = FALSE)
  }
  round(as.numeric(x))
}

#' Refuse counts past the integer range of a search that counts in integers
#' @keywords internal
#' @noRd
.stop_count_range <- function(what, where) {
  stop(
    sprintf(
      "%s counts at most %s units %s, and this design needs more. Plan with fewer units there",
      what, format(.Machine$integer.max, big.mark = ","), where
    ),
    call. = FALSE
  )
}
