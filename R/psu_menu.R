#' One part of a stratum's menu of field designs
#'
#' Every distinct state the closure gives for sampling fractions from
#' `f_lo` to `f_hi`, with `forced` PSUs held in certainty: the threshold of
#' `.psu_classify()`, the whole-PSU draw of `.psu_stratum_counts()` and the
#' crossing of `.psu_crossing()`, zone by zone under zones. The sweep runs
#' in C (`src/psu_menu.c`) and evaluates each point with the R code's
#' float expressions and summation, so a listed state is the state these
#' functions give at its `f`.
#'
#' The sweep stops at the first event at or above `f_stop`, for a caller
#' that knows no state from there on is wanted.
#'
#' Each row of `states` is one state: the `f` it was first met at, its
#' certainty `set` (a column of `sets`), the sum `Q` of its certainty
#' takes, `A = sum(M_i^2 / q_i)` over its certainty PSUs, their total size
#' `held`, their number `n_c`, and the remainder draw `a`. `key` names the
#' state.
#' @keywords internal
#' @noRd
.psu_menu_part <- function(sizes, N, take, f_lo, f_hi, cutoff = 1,
                           forced = NULL, m = 0L, zone_key = NULL,
                           f_stop = Inf) {
  n <- length(sizes)
  if (is.null(forced)) forced <- logical(n)
  # The order .psu_zones() cuts in, of every PSU: the remainder's order is
  # this one with the certainty PSUs left out.
  ord <- if (is.null(zone_key)) {
    order(-sizes, seq_len(n))
  } else {
    order(zone_key, seq_len(n))
  }
  out <- .Call(
    C_psu_menu_part, as.double(sizes), as.double(N), as.double(take),
    as.double(f_lo), as.double(f_hi), as.double(cutoff), as.logical(forced),
    as.integer(m), as.integer(ord - 1L), as.double(f_stop),
    isTRUE(capabilities("long.double"))
  )
  states <- as.data.frame(out$states)
  sig <- apply(out$sets, 2L, function(x) paste(which(x), collapse = ","))
  states$key <- paste(sig[states$set], states$Q, states$a, sep = "|")
  states <- states[!duplicated(states$key), , drop = FALSE]
  rownames(states) <- NULL
  list(states = states, sets = out$sets, f_end = out$f_end,
       stopped = out$stopped)
}

#' A register planning problem, from the caller's checked tables
#'
#' The element problem of `.build_bethel_problem()` on the caller's tables,
#' with the take, the stage costs and the ICC set aside, as the planner and
#' the shipped assessment both read it: one working variance, design
#' effect, response rate and membership per stratum and row, each target's
#' ceiling and each objective row's weight `priority / total^2`. With them
#' come each stratum's PSU sizes, supplied flags, zone order and cutoff,
#' and its lower bound on units.
#' @keywords internal
#' @noRd
.psu_setup <- function(frame, psu, measures, targets, unit_cost = NULL,
                       alpha = 0.05, deff = 1, resp_rate = 1,
                       min_n_stratum = NULL, objective = NULL, df = NULL,
                       certainty_cutoff = NULL, n_psu_per_zone = NULL) {
  stratum <- as.character(frame[["stratum"]] %||% seq_len(nrow(frame)))
  H <- nrow(frame)
  idx_of <- split(seq_len(nrow(psu)), factor(psu$stratum, levels = stratum))
  take <- .psu_take(frame, psu, stratum)
  zone_m <- .psu_zone_m(n_psu_per_zone)
  zone_key <- .psu_zone_key(psu, zone_m)
  flagged <- if ("certainty" %in% names(psu)) psu$certainty else
    logical(nrow(psu))
  icc_row <- .psu_icc(measures, stratum, idx_of)
  stage_cost <- c("cost_psu", "cost_ssu") %in% names(frame)
  if (all(stage_cost) && !is.null(unit_cost)) {
    stop(
      "certainty-aware allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }
  cost_psu <- NULL
  cost_ssu <- NULL
  if (all(stage_cost)) {
    cost_psu <- rep_len(as.numeric(frame$cost_psu), H)
    cost_ssu <- rep_len(as.numeric(frame$cost_ssu), H)
    if (any(!is.finite(cost_psu)) || any(cost_psu <= 0) ||
          any(!is.finite(cost_ssu)) || any(cost_ssu <= 0)) {
      stop("'cost_psu' and 'cost_ssu' must be positive and finite",
           call. = FALSE)
    }
  }
  base_frame <- frame[
    , setdiff(names(frame),
              c("n_per_psu", "cost_psu", "cost_ssu", "certainty_cutoff")),
    drop = FALSE
  ]
  base_measures <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  problem <- .build_bethel_problem(
    frame = base_frame, measures = base_measures, targets = targets,
    unit_cost = unit_cost, deff = deff, resp_rate = resp_rate, alpha = alpha,
    df = df, min_n_stratum = min_n_stratum, objective = objective,
    fpc = "unit", .finite_prop_var = FALSE
  )
  key <- .bethel_key(as.character(measures$stratum),
                     as.character(measures$name))
  rows <- function(name, member, var, deff_hk, resp_hk, weight) {
    icc <- matrix(0, H, length(name))
    for (k in seq_along(name)) {
      icc[member[, k], k] <-
        icc_row[match(.bethel_key(stratum, name[k]), key)][member[, k]]
    }
    w <- var * deff_hk
    w[!member] <- 0
    list(w = w, resp = resp_hk, icc = icc, weight = weight)
  }
  obj <- problem$objective
  list(
    problem = problem,
    N = frame$N,
    take = take,
    lower = problem$lower,
    idx_of = idx_of,
    sizes = lapply(idx_of, function(i) psu$N[i]),
    flagged = lapply(idx_of, function(i) flagged[i]),
    zone_key = lapply(idx_of, function(i) zone_key[i]),
    zone_m = zone_m,
    cutoff = .psu_cutoff(certainty_cutoff, frame),
    cost_psu = cost_psu,
    cost_ssu = cost_ssu,
    unit = problem$cost,
    target = rows(problem$constraint_meta$name, problem$membership,
                  problem$var_hk, problem$deff_hk, problem$resp_hk,
                  1 / problem$Vmax),
    objective = if (!is.null(obj)) {
      rows(obj$name, obj$membership, obj$var_hk, obj$deff_hk, obj$resp_hk,
           obj$priority / obj$total^2)
    }
  )
}

#' Cost and variance shares of a stratum's states
#'
#' Stratum `h`'s part of a row's variance is `w * max(v, s)`, where `v` is
#' the variance of its certainty takes and remainder draw and
#' `s = N^2 (1 / (r n) - 1 / N)` the variance at the stratum rate, and
#' `w * s` where `s <= 0`. This is `.psu_field_deff()` times `s`, from the
#' sums `.psu_menu_part()` returns: `sum(M_i^2 (1 / (r q_i) - 1 / M_i))`
#' is `A / r - held`. Target columns are shares of the row's ceiling, so a
#' plan meets every target when each column sums to at most 1 over the
#' strata, and `objective` is the stratum's part of the priority-weighted
#' relvariance.
#' @keywords internal
#' @noRd
.psu_eval_states <- function(ev, h, states) {
  N <- ev$N[h]
  b <- ev$take[h]
  n <- states$Q + states$a * b
  rest <- N - states$held
  drawn <- states$a > 0
  shares <- function(x) {
    out <- matrix(0, length(n), ncol(x$w))
    for (k in which(x$w[h, ] > 0)) {
      r <- x$resp[h, k]
      s <- N^2 * (1 / (r * n) - 1 / N)
      v <- states$A / r - states$held
      v[drawn] <- v[drawn] + rest[drawn]^2 * (1 + x$icc[h, k] * (b * r - 1)) *
        (1 / (r * states$a[drawn] * b) - 1 / rest[drawn])
      out[, k] <- x$w[h, k] * x$weight[k] * ifelse(s <= 0, s, pmax(v, s))
    }
    out
  }
  cost <- if (is.null(ev$cost_psu)) {
    ev$unit[h] * n
  } else {
    ev$cost_psu[h] * (states$n_c + states$a) + ev$cost_ssu[h] * n
  }
  list(
    n = n,
    cost = cost,
    C = shares(ev$target),
    objective = if (!is.null(ev$objective)) rowSums(shares(ev$objective))
  )
}

#' The parts of every stratum's menu, none listed yet
#'
#' Part `t` holds the `t` largest PSUs not flagged by the caller in
#' certainty, with the flagged ones, and closes the rest as the fraction
#' runs. It ends at `cutoff b / M_(t)`, where the fraction alone makes
#' those `t` certain and the part repeats part 0 (the skip rule), and part
#' 0 ends at 1. No part holds every candidate: one candidate left in the
#' remainder is drawn with probability 1, so the part before already
#' holds it. Under unit cost a visit has no price, so voluntary certainty
#' would take nearly every PSU whole at the stratum rate. Only part 0, the
#' closure family, is listed there. A state has at most
#' `f N + K + b max(m, 1)` units, so none below the stratum's lower bound
#' on units is met before `(lower - K - b max(m, 1)) / N`, where every part
#' starts.
#' @keywords internal
#' @noRd
.psu_menu_new <- function(ev) {
  lapply(seq_along(ev$N), function(h) {
    M <- ev$sizes[[h]]
    b <- ev$take[h]
    free <- which(!ev$flagged[[h]])
    cand <- free[order(-M[free], free)]
    f_lo <- max(1e-9, (ev$lower[h] - length(M) - b * max(ev$zone_m, 1)) /
                  ev$N[h])
    parts <- if (is.null(ev$cost_psu)) 1L else max(1L, length(cand))
    hi <- c(1, pmin(1, ev$cutoff[h] * b / M[cand]))[seq_len(parts)]
    n_forced <- sum(ev$flagged[[h]]) + seq_along(hi) - 1L
    list(
      cand = cand,
      hi = hi,
      visits = n_forced + (n_forced < length(M)),
      f_end = rep(f_lo, length(hi)),
      open = hi > f_lo,
      states = NULL,
      C = NULL
    )
  })
}

#' The stop rule's bound for a stratum under given weights
#'
#' `weights` prices a state at `price * cost + sum(mu * C) + obj *
#' objective`, its reduced cost under the multipliers. A state of part `t`
#' met at fraction `f` or above has at least `max(f N, lower)` units and
#' the part's visit floor, and each of its shares is at least the
#' stratum-rate share `w N^2 (1 / (r n) - 1 / N)`, the floor of
#' `.psu_eval_states()`. Its reduced cost is therefore at least `visit
#' visits + P1 n + W / n - W0`, minimized over `n` at `max(n0, min(n*,
#' N))` with `n* = sqrt(W / P1)`.
#' @keywords internal
#' @noRd
.psu_menu_bound <- function(ev, h, weights) {
  N <- ev$N[h]
  term <- function(x, mu) {
    if (is.null(x) || !length(mu)) return(c(0, 0))
    a <- mu * x$w[h, ] * x$weight
    k <- which(a > 0)
    c(sum(a[k] * N^2 / x$resp[h, k]), sum(a[k] * N))
  }
  tw <- term(ev$target, weights$mu)
  if (!is.null(ev$objective)) {
    tw <- tw + term(ev$objective, rep(weights$obj %||% 0,
                                      ncol(ev$objective$w)))
  }
  stage <- !is.null(ev$cost_psu)
  P1 <- weights$price * if (stage) ev$cost_ssu[h] else ev$unit[h]
  list(
    N = N,
    lower = ev$lower[h],
    visit = if (stage) weights$price * ev$cost_psu[h] else 0,
    P1 = P1,
    W = tw[1],
    W0 = tw[2],
    n_star = if (P1 > 0) min(sqrt(tw[1] / P1), N) else N
  )
}

#' The bound at fractions `f` for a part with `visits` visits
#'
#' `f N` carries a relative slack of 1e-9 for the float rounding of a
#' draw read off it.
#' @keywords internal
#' @noRd
.psu_menu_lb <- function(bd, visits, f) {
  n <- pmax(f * bd$N * (1 - 1e-9), bd$lower, bd$n_star)
  bd$visit * visits + bd$P1 * n + bd$W / n - bd$W0
}

#' A limit raised by a relative 1e-9 of the bound's terms
#'
#' A state can sit exactly on the bound, a census at `n* = N` for one, and
#' the two sides of that tie are computed along different float paths.
#' Every comparison with the bound uses this limit, so the stop and the
#' check that a stop is final agree.
#' @keywords internal
#' @noRd
.psu_menu_slack <- function(bd, limit) {
  limit + 1e-9 * (abs(limit) + bd$W0)
}

#' The fraction from which the bound exceeds `limit`, or `Inf`
#'
#' The bound is flat up to `n*` and rises after it, so it passes the limit
#' at the larger root of `P1 n^2 - L n + W`. The caller has checked that
#' it does not at the part's start.
#' @keywords internal
#' @noRd
.psu_menu_f_stop <- function(bd, visits, limit) {
  if (bd$P1 <= 0 || bd$n_star >= bd$N) return(Inf)
  L <- limit - bd$visit * visits + bd$W0
  n <- (L + sqrt(max(L^2 - 4 * bd$P1 * bd$W, 0))) / (2 * bd$P1)
  f <- n / (bd$N * (1 - 1e-9))
  step <- 1e-9
  while (f < 1 && !(.psu_menu_lb(bd, visits, f) > limit)) {
    f <- f * (1 + step)
    step <- 2 * step
  }
  if (f >= 1) Inf else f
}

#' List every open part of every stratum up to its stop
#'
#' A part is listed, or resumed from where it last stopped, unless the
#' bound already exceeds `limit[h]` there, and it stops at the first event
#' past which the bound exceeds it. Every state left out therefore has a
#' reduced cost above `limit[h]` under `weights`. New states are kept once
#' per stratum, those below the lower bound on units dropped, and scored
#' as they are met.
#' @keywords internal
#' @noRd
.psu_menu_round <- function(menu, ev, weights, limit) {
  for (h in seq_along(menu)) {
    mh <- menu[[h]]
    bd <- .psu_menu_bound(ev, h, weights)
    lim <- .psu_menu_slack(bd, limit[h])
    M <- ev$sizes[[h]]
    found <- list()
    for (p in which(mh$open)) {
      if (.psu_menu_lb(bd, mh$visits[p], mh$f_end[p]) > lim) next
      forced <- ev$flagged[[h]]
      forced[mh$cand[seq_len(p - 1L)]] <- TRUE
      r <- .psu_menu_part(
        M, ev$N[h], ev$take[h], mh$f_end[p], mh$hi[p], ev$cutoff[h], forced,
        ev$zone_m, ev$zone_key[[h]],
        f_stop = .psu_menu_f_stop(bd, mh$visits[p], lim)
      )
      mh$f_end[p] <- r$f_end
      mh$open[p] <- r$stopped
      found[[length(found) + 1L]] <- r$states[, c("f", "Q", "A", "held",
                                                  "n_c", "a", "key")]
    }
    if (length(found)) {
      d <- do.call(rbind, found)
      d <- d[!duplicated(d$key) & !d$key %in% mh$states$key, , drop = FALSE]
      d <- d[d$Q + d$a * ev$take[h] >= ev$lower[h] - 1e-8, , drop = FALSE]
      if (nrow(d)) {
        e <- .psu_eval_states(ev, h, d)
        d$n <- e$n
        d$cost <- e$cost
        d$objective <- e$objective %||% 0
        mh$states <- rbind(mh$states, d)
        mh$C <- rbind(mh$C, e$C)
      }
    }
    menu[[h]] <- mh
  }
  menu
}

#' Whether a round under `limit` would list anything more
#' @keywords internal
#' @noRd
.psu_menu_short <- function(menu, ev, weights, limit) {
  any(vapply(seq_along(menu), function(h) {
    mh <- menu[[h]]
    p <- which(mh$open)
    if (!length(p)) return(FALSE)
    bd <- .psu_menu_bound(ev, h, weights)
    lim <- .psu_menu_slack(bd, limit[h])
    any(.psu_menu_lb(bd, mh$visits[p], mh$f_end[p]) <= lim)
  }, logical(1)))
}

#' Each stratum's listed states with the plan's, as one menu problem
#'
#' `extra` holds states that must stay in the menus whatever the stop rule
#' leaves out, the plan found first among them, so it stays a candidate.
#' @keywords internal
#' @noRd
.psu_menu_rows <- function(menu, extra = NULL) {
  lapply(seq_along(menu), function(h) {
    d <- menu[[h]]$states
    C <- menu[[h]]$C
    e <- extra[[h]]
    if (!is.null(e)) {
      new <- !e$states$key %in% d$key
      d <- rbind(d, e$states[new, , drop = FALSE])
      C <- rbind(C, e$C[new, , drop = FALSE])
    }
    list(states = d, C = C)
  })
}

.psu_target_problem <- function(rows) {
  .menu_problem(
    unlist(lapply(rows, function(r) r$states$cost), use.names = FALSE),
    do.call(rbind, lapply(rows, `[[`, "C")),
    rep(seq_along(rows), vapply(rows, function(r) nrow(r$states), 1L))
  )
}

#' The budget problem over the same menus
#'
#' The cost a plan minimizes is the objective, and field cost over the
#' budget is the first column, ahead of the target columns.
#' @keywords internal
#' @noRd
.psu_budget_problem <- function(rows, budget) {
  money <- unlist(lapply(rows, function(r) r$states$cost), use.names = FALSE)
  C <- do.call(rbind, lapply(rows, `[[`, "C"))
  if (!NCOL(C)) C <- matrix(0, length(money), 0L)
  .menu_problem(
    unlist(lapply(rows, function(r) r$states$objective), use.names = FALSE),
    cbind(money / budget, C),
    rep(seq_along(rows), vapply(rows, function(r) nrow(r$states), 1L))
  )
}

#' The step between attainable plan costs, 0 when there is none
#'
#' Every plan cost is a whole combination of the per-unit and per-visit
#' prices, so their greatest common divisor is the step, once the prices
#' are scaled to whole numbers by a power of 10 up to 1e6. A scaled price
#' counts as whole only to double precision, as a looser test rounds away a
#' real decimal part, of a large price or of one near a whole number, and a
#' price scaled past exact double integers gives no step. 0 is always safe,
#' as it only slows the search.
#' @keywords internal
#' @noRd
.menu_grain <- function(x) {
  x <- unique(x[x != 0])
  if (!length(x)) return(0)
  gcd <- function(a, b) if (b == 0) a else gcd(b, a %% b)
  for (d in 0:6) {
    y <- x * 10^d
    if (any(abs(y) >= 2^50)) return(0)
    if (all(abs(y - round(y)) <= 64 * .Machine$double.eps * pmax(1, abs(y)))) {
      return(Reduce(gcd, round(abs(y))) / 10^d)
    }
  }
  0
}

#' Whether every cost is a whole multiple of `grain`
#' @keywords internal
#' @noRd
.on_grain <- function(x, grain) {
  k <- x / grain
  all(abs(k - round(k)) <= 64 * .Machine$double.eps * pmax(1, abs(k)))
}

#' The best field design, by listing and search
#'
#' In target mode the plan is the cheapest design meeting every target. In
#' budget mode it is the design with the smallest priority-weighted
#' relvariance whose field cost stays within `budget`, every target in
#' `ev$target` held as a hard constraint (`.psu_budget_problem()`).
#'
#' A first pass lists each stratum's closure family (part 0) up to four
#' times `f_anchor`, the fraction of a continuous solve, and its
#' multipliers and plan start the listing. In target mode the band widens
#' when it holds no plan, up to a census of every stratum, which gives each
#' row its smallest variance, so a census that misses the targets proves
#' that no design meets them. In budget mode the start comes from
#' `.psu_budget_start()`. The rounds then list every part under the stop
#' rule at `rmin + gamma`, first with `gamma = 0` and then with the `gamma
#' = UB - grain - L` the plan certifies, so a state left out is in no plan
#' better than UB by a grain. The search runs on what was listed. Budget
#' mode prices a state at `objective + nu cost / budget + lambda C`, so
#' the stop rule takes `price = nu / budget` and `obj = 1`, and its
#' objective has no grain.
#' @keywords internal
#' @noRd
.psu_plan <- function(ev, f_anchor, max_nodes, budget = NULL, gap = 0) {
  H <- length(ev$N)
  K <- ncol(ev$target$w)
  f_lo <- vapply(.psu_menu_new(ev), function(mh) mh$f_end[1], 0)
  band <- 4
  repeat {
    first <- .psu_menu_new(ev)
    for (h in seq_len(H)) {
      hi <- min(1, band * max(f_anchor[h], f_lo[h]))
      first[[h]]$hi[1] <- hi
      first[[h]]$open <- seq_along(first[[h]]$hi) == 1L & hi > f_lo[h]
    }
    first <- .psu_menu_round(first, ev, list(price = 1, mu = numeric(K)),
                             rep(Inf, H))
    rows1 <- .psu_menu_rows(first)
    if (!is.null(budget)) break
    p1 <- .psu_target_problem(rows1)
    sg1 <- .menu_subgradient(p1)
    if (!is.null(sg1$pick)) break
    if (all(band * pmax(f_anchor, f_lo) >= 1)) {
      stop(
        "no field design meets every target, not even a census of every PSU, under the response rates and design effects given",
        call. = FALSE
      )
    }
    band <- band * 4
  }
  if (is.null(budget)) {
    problem_of <- .psu_target_problem
    weights_of <- function(lambda) list(price = 1, mu = lambda, obj = 0)
    prices <- if (is.null(ev$cost_psu)) ev$unit else
      c(ev$cost_psu, ev$cost_ssu)
    # A certainty take capped at a fractional PSU size has a fractional
    # cost, off the lattice of the prices.
    sizes <- unlist(ev$sizes, use.names = FALSE)
    whole <- all(abs(sizes - round(sizes)) <= 1e-9) &&
      all(abs(ev$take - round(ev$take)) <= 1e-9)
    grain <- if (whole) .menu_grain(prices) else 0
    if (grain > 0 && !.on_grain(p1$cost, grain)) grain <- 0
  } else {
    problem_of <- function(rows) .psu_budget_problem(rows, budget)
    weights_of <- function(lambda) {
      list(price = lambda[1] / budget, mu = lambda[-1], obj = 1)
    }
    grain <- 0
    start <- .psu_budget_start(ev, rows1, f_anchor, max_nodes, budget, gap)
    rows1 <- .psu_menu_rows(first, start)
    p1 <- problem_of(rows1)
    sg1 <- .menu_subgradient(p1, start = .psu_plan_index(
      p1, rows1, vapply(start, function(e) e$states$key, "")
    ))
  }
  lambda <- sg1$lambda
  weights <- weights_of(lambda)
  extra <- .psu_plan_rows(p1, rows1, sg1$pick)
  plan_keys <- vapply(extra, function(e) e$states$key, "")
  rmin_of <- function(p) {
    l <- .menu_lagrangian(p, lambda)
    l$r[l$best]
  }
  rmin <- rmin_of(p1)
  ub <- sg1$ub
  gamma <- 0
  menu <- .psu_menu_new(ev)
  certified <- FALSE
  rounds <- 0L
  repeat {
    rounds <- rounds + 1L
    if (rounds > 20L) {
      stop("internal error: the menu listing did not settle", call. = FALSE)
    }
    menu <- .psu_menu_round(menu, ev, weights, rmin + gamma)
    rows <- .psu_menu_rows(menu, extra)
    p <- problem_of(rows)
    rmin <- pmin(rmin, rmin_of(p))
    if (.psu_menu_short(menu, ev, weights, rmin + gamma)) next
    if (certified) break
    start <- .psu_plan_index(p, rows, plan_keys)
    sg <- .menu_subgradient(p, lambda, start = start)
    plan_keys <- .psu_plan_keys(p, rows, sg$pick)
    ub <- min(ub, sg$ub)
    gamma <- max(0, ub - grain - (sum(rmin) - sum(lambda)))
    certified <- TRUE
    if (!.psu_menu_short(menu, ev, weights, rmin + gamma)) break
  }
  start <- .psu_plan_index(p, rows, plan_keys)
  sg2 <- .menu_subgradient(p, sg$lambda, start = start)
  sr <- .menu_search(p, sg2$lambda, sg2$pick, grain = grain,
                     max_nodes = max_nodes, gap = gap)
  if (is.null(sr$pick)) {
    stop("internal error: the search lost the plan it started from",
         call. = FALSE)
  }
  list(
    rows = lapply(.psu_plan_rows(p, rows, sr$pick), `[[`, "states"),
    value = sr$value,
    bound = sr$bound,
    status = sr$status,
    nodes = sr$nodes,
    kept = sr$kept,
    states = length(p$cost),
    max_nodes = max_nodes
  )
}

#' A plan inside the budget to start the budget search from
#'
#' The cheapest plan takes each stratum's cheapest state. One round prices
#' states on field cost alone, with each stratum's cheapest state of the
#' first pass as its limit, so it lists every state that costs no more,
#' and a budget below that plan is an error. With hard targets the start
#' is a plan meeting them inside the budget: the first pass's cheapest one
#' when it fits, and otherwise the target-mode plan, whose bound proves a
#' budget too small to fund the targets.
#' @keywords internal
#' @noRd
.psu_budget_start <- function(ev, rows1, f_anchor, max_nodes, budget,
                              gap = 0) {
  H <- length(ev$N)
  K <- ncol(ev$target$w)
  fits <- function(cost) cost <= budget * (1 + 1e-9)
  lim <- vapply(rows1, function(r) {
    if (NROW(r$states)) min(r$states$cost) else Inf
  }, 0)
  cheap <- .psu_menu_rows(.psu_menu_round(
    .psu_menu_new(ev), ev, list(price = 1, mu = numeric(K), obj = 0), lim
  ))
  plan <- lapply(cheap, function(r) {
    j <- order(r$states$cost, r$states$objective)[1]
    list(states = r$states[j, , drop = FALSE], C = r$C[j, , drop = FALSE])
  })
  cheapest <- sum(vapply(plan, function(e) e$states$cost, 0))
  if (!fits(cheapest)) {
    stop(sprintf(
      "'budget' is below the cheapest field design, which costs %.6g",
      cheapest
    ), call. = FALSE)
  }
  if (K == 0L) return(plan)
  p1 <- .psu_target_problem(rows1)
  sg <- .menu_subgradient(p1)
  if (!is.null(sg$pick) && fits(sg$ub)) {
    return(.psu_plan_rows(p1, rows1, sg$pick))
  }
  tp <- .psu_plan(ev, f_anchor, max_nodes, gap = gap)
  if (fits(tp$value)) {
    return(lapply(seq_len(H), function(h) {
      s <- tp$rows[[h]]
      list(states = s, C = .psu_eval_states(ev, h, s)$C)
    }))
  }
  if (!fits(tp$bound)) {
    stop(sprintf(
      "'budget' cannot fund the precision targets: the cheapest field design meeting them costs %.6g, short by %.6g",
      tp$value, tp$value - budget
    ), call. = FALSE)
  }
  # Classed, so a predict() grid reads the row as unknown, not infeasible.
  stop(structure(
    list(message = sprintf(
      "no field design meeting the precision targets within 'budget' was found in %d search node%s, so whether one exists is unknown. The cheapest found costs %.6g and none costs less than %.6g. Raise 'control$max_nodes'",
      as.integer(max_nodes), if (max_nodes == 1) "" else "s", tp$value,
      tp$bound
    ), call = NULL),
    class = c("svyplan_plan_unknown", "error", "condition")
  ))
}

#' Each stratum's state and loads in a plan of a menu problem
#' @keywords internal
#' @noRd
.psu_plan_rows <- function(p, rows, pick) {
  lapply(seq_along(rows), function(h) {
    j <- pick[h] - p$groups[[h]][1] + 1L
    list(states = rows[[h]]$states[j, , drop = FALSE],
         C = rows[[h]]$C[j, , drop = FALSE])
  })
}

.psu_plan_index <- function(p, rows, keys) {
  vapply(seq_along(rows), function(h) {
    p$groups[[h]][match(keys[h], rows[[h]]$states$key)]
  }, 1L)
}

.psu_plan_keys <- function(p, rows, pick) {
  vapply(seq_along(rows), function(h) {
    rows[[h]]$states$key[pick[h] - p$groups[[h]][1] + 1L]
  }, "")
}

#' The field design of a planned state, rebuilt with the R closure
#'
#' The state's certainty set is read from its key and its takes and draw
#' recomputed at its fraction `f`, where the sweep met it, as
#' `.psu_stratum_counts()` gives them. A PSU at or above the threshold at
#' `f` reads `"threshold"`, one the caller flagged `"supplied"`, one the
#' draw pushes to probability one without any planned PSU
#' `"operational"`, and the rest `"planned"`.
#' @keywords internal
#' @noRd
.psu_design <- function(ev, h, state) {
  M <- ev$sizes[[h]]
  N <- ev$N[h]
  b <- ev$take[h]
  f <- state$f
  n_h <- f * N
  set <- strsplit(state$key, "|", fixed = TRUE)[[1]][1]
  certain <- logical(length(M))
  if (nzchar(set)) certain[as.integer(strsplit(set, ",", fixed = TRUE)[[1]])] <- TRUE
  counts <- .psu_stratum_counts(n_h, N, b, M, certain, ev$zone_m)
  if (sum(counts$per_psu) != state$Q || counts$n_psu_draw != state$a) {
    stop("internal error: a planned state does not rebuild", call. = FALSE)
  }
  threshold <- ev$cutoff[h] * b / f
  above <- M >= threshold
  plain <- .psu_draw_closure(n_h, N, b, M, above | ev$flagged[[h]],
                             ev$cutoff[h], ev$zone_m, ev$zone_key[[h]])$certain
  source <- ifelse(above, "threshold",
                   ifelse(ev$flagged[[h]], "supplied",
                          ifelse(plain, "operational", "planned")))
  source[!certain] <- NA_character_
  list(certain = certain, per_psu = counts$per_psu,
       n_psu_draw = counts$n_psu_draw, n_certain_int = sum(counts$per_psu),
       threshold = threshold, source = source)
}
