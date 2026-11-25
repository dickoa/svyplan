#' @rdname n_cluster
#' @name n_cluster
#' @aliases n_cluster.data.frame
#' @export
n_cluster.data.frame <- function(
  x,
  stage_cost = NULL,
  ...,
  domains = NULL,
  budget = NULL,
  stage_n = NULL,
  allocation = c("separate", "joint"),
  domain_sampling = c("separate", "natural"),
  min_n_domain = NULL,
  fixed_cost = 0,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  N = NULL,
  fpc = c("auto", "none", "stage"),
  control = NULL,
  plan = NULL
) {
  dots <- list(...)
  col_carried <- c(
    icc = "icc_1", unit_relvar = "unit_relvar",
    var_ratio = "var_ratio_1", cv = "cv"
  )
  hit <- intersect(names(dots), names(col_carried))
  if (length(hit) > 0L) {
    .stop_carried_by_column(hit, col_carried)
  }
  fpc <- match.arg(fpc)
  sc_len <- if (is.numeric(stage_cost)) length(stage_cost) else NULL
  plan <- .plan_for_cluster_mode(plan, TRUE, fpc, stages = sc_len)
  .plan <- .merge_plan_args(plan, n_cluster.data.frame, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_cluster.data.frame, c(.plan, dots)))
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(x, "x")
  for (nm in c("budget", "fixed_cost", "eligibility_rate", "eligible_per_unit")) {
    v <- get(nm)
    if (is.numeric(v)) assign(nm, unname(v))
  }
  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)

  .n_cluster_indicators(
    as.data.frame(x),
    stage_cost = stage_cost,
    domains = domains,
    budget = budget,
    stage_n = stage_n,
    allocation = allocation,
    domain_sampling = domain_sampling,
    min_n_domain = min_n_domain,
    fixed_cost = fixed_cost,
    resp_rate = resp_rate,
    eligibility_rate = eligibility_rate,
    eligible_per_unit = eligible_per_unit,
    N_arg = N,
    fpc = fpc,
    control = control
  )
}

#' Extract stage matrices (R, rates, M) from a validated indicator data frame
#' @keywords internal
#' @noRd
.cluster_ind_matrices <- function(indicators, stages, pop_collapsed = NULL) {
  J <- nrow(indicators)
  n_b <- stages - 1L
  icc_mat <- matrix(0, nrow = J, ncol = n_b)
  vr_mat <- matrix(0, nrow = J, ncol = n_b)
  for (s in seq_len(n_b)) {
    icc_mat[, s] <- as.numeric(indicators[[paste0("icc_", s)]])
    vr_mat[, s] <- as.numeric(indicators[[paste0("var_ratio_", s)]])
  }
  R <- .stage_R_matrix(icc_mat, indicators$unit_relvar, vr_mat)
  rate_cols <- c(paste0("resp_rate_", seq_len(n_b)), "resp_rate")
  rates <- matrix(0, nrow = J, ncol = stages)
  for (s in seq_len(stages)) {
    rates[, s] <- as.numeric(indicators[[rate_cols[s]]])
  }
  M <- if (!is.null(pop_collapsed)) {
    pop_collapsed$M_mat
  } else {
    matrix(Inf, nrow = J, ncol = stages)
  }
  list(R = R, rates = rates, M = M)
}

#' Check expected takes across all rows of an indicator table for stages 2..S
#' @keywords internal
#' @noRd
.check_ind_expected_takes <- function(indicators, stage_n, stages, stage_nms) {
  if (is.null(stage_n)) return(invisible(TRUE))
  nr <- nrow(indicators)
  labels <- if ("name" %in% names(indicators)) indicators$name else seq_len(nr)
  rate_cols <- c(paste0("resp_rate_", seq_len(stages - 1L)), "resp_rate")
  for (s in 2L:stages) {
    if (is.na(stage_n[s])) next
    r_col <- indicators[[rate_cols[s]]]
    for (j in seq_len(nr)) {
      where <- sprintf(" for indicator %s", sQuote(labels[j]))
      .check_expected_take(stage_n[s], r_col[j], stage_nms[s], rate_cols[s],
                           where)
    }
  }
  invisible(TRUE)
}

#' Minimum variable cost to reach total_n >= min_total in S stages
#'
#' Continuous minimum cost greedily fills innermost stages first (since
#' earlier-stage unit costs are amortized by all deeper takes). Whole-unit
#' minimum cost uses the exact C branch-and-bound solver with zero variance,
#' and is `NA` when that search stopped at its node limit.
#' @keywords internal
#' @noRd
.stage_min_total_cost <- function(min_total, cost_vec, lower, upper,
                                  whole = FALSE, max_nodes = 2e7) {
  stages <- length(cost_vec)
  if (min_total <= 0) {
    return(.stage_cost(lower, cost_vec, 0))
  }
  if (all(is.finite(upper)) && prod(upper) < min_total * (1 - 1e-12)) {
    return(Inf)
  }
  if (!whole) {
    n <- lower
    # Fill from stage S back to stage 1
    for (s in rev(seq_len(stages))) {
      cur <- prod(n)
      if (cur >= min_total * (1 - 1e-12)) break
      need <- n[s] * (min_total / cur)
      n[s] <- min(upper[s], max(lower[s], need))
    }
    if (prod(n) < min_total * (1 - 1e-9)) return(Inf)
    return(.stage_cost(n, cost_vec, 0))
  }
  w <- .stage_min_total_whole(min_total, cost_vec, lower, upper, max_nodes)
  if (!w$exact) return(NA_real_)
  if (is.null(w$n)) return(Inf)
  w$cost
}

#' The cheapest whole design reaching `min_total`, `n = NULL` when none is found
#' @keywords internal
#' @noRd
.stage_min_total_whole <- function(min_total, cost_vec, lower, upper,
                                   max_nodes = 2e7) {
  stages <- length(cost_vec)
  if (min_total <= 0) {
    return(list(exact = TRUE, n = lower, cost = .stage_cost(lower, cost_vec, 0),
                nodes = 0))
  }
  if (all(is.finite(upper)) && prod(upper) < min_total * (1 - 1e-12)) {
    return(list(exact = TRUE, n = NULL, cost = Inf, nodes = 0))
  }
  R0 <- matrix(0, nrow = 1L, ncol = stages)
  rates0 <- matrix(1, nrow = 1L, ncol = stages)
  M0 <- matrix(Inf, nrow = 1L, ncol = stages)
  n_init <- lower
  for (s in rev(seq_len(stages))) {
    cur <- prod(n_init)
    if (cur >= min_total) break
    n_init[s] <- min(upper[s], max(lower[s], ceiling(n_init[s] * min_total / cur)))
  }
  disc <- tryCatch(
    .stage_solve_discrete(
      R0, rates0, M0, 1.0, cost_vec, lower, upper,
      budget_var = NA_real_, min_total = min_total, n_init = n_init,
      max_nodes = max_nodes
    ),
    svyplan_plan_unknown = function(e) list(found = FALSE, exact = FALSE,
                                            nodes = e$nodes %||% max_nodes)
  )
  n <- if (isTRUE(disc$found)) as.numeric(disc$n)
  list(exact = isTRUE(disc$exact), n = n,
       cost = if (is.null(n)) Inf else .stage_cost(n, cost_vec, 0),
       nodes = disc$nodes)
}

#' The whole stage bounds of one domain, from its response floors, caps and fixed stages
#' @keywords internal
#' @noRd
.stage_whole_box <- function(rates, pop_c, stage_n) {
  stages <- ncol(rates)
  lower <- c(1, vapply(2L:stages, function(s) .take_floor(rates[, s]),
                       integer(1L)))
  upper <- if (!is.null(pop_c)) floor(pop_c$M_take + 1e-9) else rep(Inf, stages)
  if (!is.null(stage_n)) {
    for (s in which(!is.na(stage_n))) {
      w <- if (s == 1L) {
        min(upper[1L], .whole_psu(stage_n[1L]))
      } else {
        .whole_take(stage_n[s], lower[s], upper[s])
      }
      lower[s] <- w
      upper[s] <- w
    }
  }
  list(lower = lower, upper = upper)
}

#' The least worst ratio over whole designs within a joint budget
#'
#' Bisects on the ratio `t` between a proven infeasible ratio and a feasible
#' design. At each `t`, every domain's exact cost-mode search gives its least
#' whole cost at targets `t * cv`. `max_nodes` bounds all the trials together.
#' @keywords internal
#' @noRd
.joint_whole_minimax <- function(dom_data, cost_vec, var_budget, min_total,
                                 witness, max_nodes) {
  nd <- length(dom_data)
  limit <- .budget_room(var_budget)
  used <- numeric(nd)
  ratio_of <- function(d, n) {
    dd <- dom_data[[d]]
    max(sqrt(.stage_cv2(n, dd$R, dd$rates, dd$M)) / dd$cv_t)
  }
  trial <- function(t) {
    designs <- vector("list", nd)
    exact <- TRUE
    for (d in seq_len(nd)) {
      dd <- dom_data[[d]]
      left <- max_nodes - sum(used)
      if (left < 1) return(list(out = "unknown"))
      disc <- tryCatch(
        .stage_solve_discrete(
          dd$R, dd$rates, dd$M, (t * dd$cv_t)^2, cost_vec, dd$lower_int,
          dd$upper_int, budget_var = NA_real_, min_total = min_total,
          max_nodes = left
        ),
        svyplan_plan_unknown = function(e) list(nodes = e$nodes %||% left)
      )
      used[d] <<- used[d] + disc$nodes
      if (is.null(disc$found)) return(list(out = "unknown"))
      if (!isTRUE(disc$found)) return(list(out = "infeasible"))
      exact <- exact && isTRUE(disc$exact)
      designs[[d]] <- disc
    }
    n <- lapply(designs, function(x) as.numeric(x$n))
    cost <- sum(vapply(n, .stage_cost, 1, cost_vec, 0))
    if (cost <= limit) {
      return(list(out = "feasible", n = n, cost = cost, discs = designs,
                  ratio = max(vapply(seq_len(nd), function(d) ratio_of(d, n[[d]]),
                                     1))))
    }
    list(out = if (exact) "infeasible" else "unknown")
  }
  best <- list(n = witness, ratio = max(vapply(seq_len(nd), function(d)
    ratio_of(d, witness[[d]]), 1)), discs = NULL)
  keep <- function(tr) if (tr$ratio < best$ratio) tr else best
  ## The native search accepts a cv up to sqrt(1 + 1e-9) of its target, so
  ## best$ratio may sit just above the endpoint t_hi.
  t_hi <- best$ratio
  lo <- NA_real_
  status <- "optimal"
  trials <- 0L
  if (t_hi > 0) {
    t <- t_hi
    for (k in seq_len(60L)) {
      t <- t / 2
      trials <- trials + 1L
      tr <- trial(t)
      if (tr$out == "feasible") {
        best <- keep(tr)
        t_hi <- min(t, tr$ratio)
        t <- t_hi
        if (best$ratio == 0) break
      } else {
        if (tr$out == "infeasible") lo <- t else status <- "limit"
        break
      }
    }
    if (is.na(lo) && best$ratio > 0) status <- "limit"
  }
  while (identical(status, "optimal") && t_hi > 0 && t_hi / lo > 1 + 1e-9) {
    m <- sqrt(lo * t_hi)
    trials <- trials + 1L
    tr <- trial(m)
    if (tr$out == "feasible") {
      best <- keep(tr)
      t_hi <- min(m, tr$ratio)
    } else if (tr$out == "infeasible") {
      lo <- m
    } else {
      status <- "limit"
    }
  }
  bound <- if (best$ratio == 0 || is.na(lo)) 0 else min(lo, best$ratio)
  list(n = best$n, discs = best$discs, ratio = best$ratio, bound = bound,
       gap = if (best$ratio > 0) (best$ratio - bound) / best$ratio else 0,
       status = status, trials = trials, nodes = sum(used),
       nodes_by_domain = used)
}

#' Cluster allocation from a table of indicators (arbitrary S >= 2)
#' @keywords internal
#' @noRd
.n_cluster_indicators <- function(
  indicators,
  stage_cost = NULL,
  domains = NULL,
  budget = NULL,
  stage_n = NULL,
  allocation = c("separate", "joint"),
  domain_sampling = c("separate", "natural"),
  min_n_domain = NULL,
  fixed_cost = 0,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  N_arg = NULL,
  fpc = "auto",
  control = NULL
) {
  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }
  indicators_input <- indicators
  max_nodes <- .stage_max_nodes(control)
  if (is.null(stage_cost)) {
    stop("'stage_cost' is required (directly or via plan)", call. = FALSE)
  }
  if (!is.character(allocation) || anyNA(allocation)) {
    stop("'allocation' must be \"separate\" or \"joint\"", call. = FALSE)
  }
  allocation <- match.arg(allocation)
  domain_sampling <- match.arg(domain_sampling)
  joint <- identical(allocation, "joint")
  if (joint && (is.null(domains) || is.null(budget))) {
    stop(
      sprintf(
        "allocation = \"joint\" splits one budget across domains, so it needs %s",
        paste(sQuote(c("domains", "budget")[c(is.null(domains), is.null(budget))]),
              collapse = " and ")
      ),
      call. = FALSE
    )
  }
  if (identical(domain_sampling, "natural")) {
    stop(
      "domain_sampling = \"natural\" is not supported for multistage cluster designs; use domain_sampling = \"separate\" (or n_multi() for single-stage designs)",
      call. = FALSE
    )
  }
  if (!is.null(min_n_domain) &&
      (!is.numeric(min_n_domain) || length(min_n_domain) != 1L ||
       is.na(min_n_domain) || min_n_domain <= 0)) {
    stop("'min_n_domain' must be a positive numeric scalar", call. = FALSE)
  }
  if (!is.null(min_n_domain) && is.null(domains)) {
    stop("'min_n_domain' needs 'domains'", call. = FALSE)
  }

  check_stage_cost(stage_cost)
  stage_cost <- .reorder_stage_cost(stage_cost)
  stages <- length(stage_cost)
  stage_nms <- .stage_names(stages, names(stage_cost))
  names(stage_cost) <- stage_nms
  resp_vec <- .expand_stage_resp_rate(resp_rate, stages, stage_nms)

  if (!is.null(budget)) check_scalar(budget, "budget")
  check_fixed_cost(fixed_cost, budget)
  sn <- .validate_stage_n(stage_n, stages, stage_nms)

  indicators <- .indicators_var_from_sd(indicators, domains)
  .stop_min_cases_column(
    indicators,
    "a multistage design sizes stages against a cost, not a single total"
  )
  provenance <- .target_provenance(indicators)
  indicators <- .indicators_moe_from_rmoe(indicators, domains)
  info <- .validate_targets(
    indicators,
    TRUE,
    domains = domains,
    stages = stages,
    context = "n_cluster()"
  )
  indicators <- .fill_defaults(
    indicators,
    TRUE,
    stages = stages,
    resp_rate = resp_vec
  )
  indicators <- .fill_eligible_per_unit(indicators, eligible_per_unit)
  indicators <- .convert_moe_to_cv(indicators)

  given_resp <- indicators$resp_rate
  epu_vec <- .indicator_rate(indicators, "eligible_per_unit", nrow(indicators))
  yield_vec <- eligibility_rate * epu_vec
  indicators$resp_rate <- given_resp * yield_vec

  rv_final <- indicators$unit_relvar[!is.na(indicators$unit_relvar)]
  if (length(rv_final) > 0L && (any(rv_final <= 0) || any(!is.finite(rv_final)))) {
    stop("'unit_relvar' values must be positive and finite", call. = FALSE)
  }
  for (s in seq_len(stages - 1L)) {
    col <- paste0("var_ratio_", s)
    if (any(indicators[[col]] <= 0) || any(!is.finite(indicators[[col]]))) {
      stop(sprintf("'%s' values must be positive and finite", col), call. = FALSE)
    }
  }

  domain_cols <- info$domain_cols
  mode <- if (!is.null(budget)) {
    "budget"
  } else if ("moe" %in% names(indicators) && any(!is.na(indicators$moe))) {
    "moe"
  } else {
    "cv"
  }

  row_pops <- .stage_row_pops(indicators, domain_cols, fpc, N_arg, stages,
                              yield_vec)
  res <- if (length(domain_cols) == 0L) {
    pop_c <- if (!is.null(row_pops)) .collapse_row_pops(row_pops) else NULL
    .n_cluster_multi_solve(
      indicators, stage_cost, budget, sn, fixed_cost,
      pop_collapsed = pop_c, domain_cols = domain_cols, mode = mode,
      record_pop = !is.null(N_arg), max_nodes = max_nodes
    )
  } else {
    .n_cluster_domains(
      indicators, stage_cost, budget, sn, domain_cols,
      joint = joint, min_n_domain = min_n_domain, fixed_cost = fixed_cost,
      mode = mode, domain_sampling = domain_sampling, row_pops = row_pops,
      max_nodes = max_nodes
    )
  }

  res$indicators$resp_rate <- given_resp
  res$params$resp_rate <- .name_stage_resp_rate(resp_vec)
  res$params$prop_method <- "wald"
  if (eligibility_rate != 1) res$params$eligibility_rate <- eligibility_rate
  if (eligible_per_unit != 1 || "eligible_per_unit" %in% names(indicators)) {
    res$params$eligible_per_unit <- eligible_per_unit
  }
  res$params$indicators_input <- indicators_input
  if (!is.null(control)) res$params$control <- control
  .attach_target_provenance(res, provenance, indicators, domain_cols)
}

#' Refuse cv targets below their floor, naming the indicator and the remedy
#'
#' The stage-1 size that meets a target is `G / (cv^2 + kappa)`.
#' @keywords internal
#' @noRd
.stage_floor_stop <- function(R, rates, M, cv_t, lower, upper, labels,
                              stage_nms) {
  floor_cv <- sqrt(pmax(0, .stage_cv2(upper, R, rates, M)))
  bad <- which(floor_cv >= cv_t * (1 - 1e-9))
  if (!length(bad)) {
    stop(
      "n_cluster(): precision targets are not achievable under the given stage constraints",
      call. = FALSE
    )
  }
  j <- bad[1L]
  tel <- .stage_telescoping(R, rates, M)
  S <- ncol(R)
  q <- c(1, cumprod(upper[-1L]))
  G <- sum(tel$alpha[j, ] / q)
  need <- ceiling(G / (cv_t[j]^2 + tel$kappa[j]) - 1e-9)
  fixed <- which(lower == upper)
  at <- if (length(fixed)) {
    sprintf(" at %s", paste(stage_nms[fixed], "=", .fmt_count_n(upper[fixed]),
                            collapse = ", "))
  } else {
    ""
  }
  remedy <- if (1L %in% fixed && is.finite(need) && need > upper[1L] &&
                need <= (M[j, 1L] %||% Inf)) {
    sprintf("raise stage_n[1] to at least %s, or relax", .fmt_count_n(need))
  } else {
    "relax"
  }
  stop(
    sprintf(
      "n_cluster(): target CV %.4g for indicator '%s' is below the achievable floor %.4g%s; %s the target above %.4g%s",
      cv_t[j], labels[j], floor_cv[j], at, remedy, floor_cv[j],
      if (length(bad) > 1L) {
        sprintf(" (%d other indicator%s also miss%s)", length(bad) - 1L,
                if (length(bad) > 2L) "s" else "",
                if (length(bad) > 2L) "" else "es")
      } else {
        ""
      }
    ),
    call. = FALSE
  )
}

#' Solve a single-population S-stage multi-indicator cluster problem
#' @keywords internal
#' @noRd
.n_cluster_multi_solve <- function(
  indicators,
  stage_cost,
  budget = NULL,
  stage_n = NULL,
  fixed_cost = 0,
  pop_collapsed = NULL,
  domain_cols = character(0),
  mode = "cv",
  record_pop = FALSE,
  min_total = 0,
  max_nodes = 2e7,
  whole_optional = FALSE
) {
  stages <- length(stage_cost)
  stage_nms <- .stage_names(stages, names(stage_cost))
  cost_vec <- unname(as.numeric(stage_cost))
  names(stage_cost) <- stage_nms

  mats <- .cluster_ind_matrices(indicators, stages, pop_collapsed)
  R <- mats$R
  rates <- mats$rates
  M <- mats$M
  cv_t <- indicators$cv

  # Check fixed stage sizes against population caps & expected takes
  if (!is.null(stage_n)) {
    if (!is.null(pop_collapsed)) {
      for (s in seq_len(stages)) {
        if (!is.na(stage_n[s]) && stage_n[s] > pop_collapsed$M_take[s] + 1e-9) {
          stop(
            sprintf(
              "fixed %s (%g) exceeds population capacity %g at stage %d",
              stage_nms[s], stage_n[s], pop_collapsed$M_take[s], s
            ),
            call. = FALSE
          )
        }
      }
    }
    .check_ind_expected_takes(indicators, stage_n, stages, stage_nms)
  }

  # Continuous bounds
  lower <- c(1, vapply(2L:stages, function(s) .take_bound(rates[, s]),
                       numeric(1L)))
  upper <- if (!is.null(pop_collapsed)) pop_collapsed$M_take else rep(Inf, stages)
  if (!is.null(stage_n)) {
    fixed <- !is.na(stage_n)
    lower[fixed] <- stage_n[fixed]
    upper[fixed] <- stage_n[fixed]
  }
  if (any(lower > upper + 1e-12)) {
    bad <- which(lower > upper + 1e-12)[1L]
    stop(
      sprintf(
        "stage %d (%s) population capacity %g cannot deliver one expected response at these response rates",
        bad, stage_nms[bad], upper[bad]
      ),
      call. = FALSE
    )
  }

  sol <- .stage_solve_continuous_multi(
    R, rates, M, cv_t, cost_vec, lower, upper,
    budget = budget, fixed_cost = fixed_cost, min_total = min_total
  )
  if (!sol$ok) {
    if (identical(sol$reason, "budget")) {
      stop(
        sprintf(
          "'budget' is too small for any realizable design: the cheapest design costs %s for total_n = %s, plus fixed_cost = %.4g",
          format(signif(.stage_cost(lower, cost_vec, 0), 6), scientific = FALSE),
          .fmt_count_n(prod(lower)), fixed_cost
        ),
        call. = FALSE
      )
    }
    if (identical(sol$reason, "min_total")) {
      stop(
        sprintf(
          "min_n_domain = %g exceeds the maximum population total %.0f",
          min_total, prod(upper)
        ),
        call. = FALSE
      )
    }
    .stage_floor_stop(R, rates, M, cv_t, lower, upper,
                      as.character(indicators[["name"]] %||% seq_len(nrow(R))),
                      stage_nms)
  }

  cont_n <- as.numeric(sol$n)
  names(cont_n) <- stage_nms
  cont_cv_by_t <- if (!is.null(sol$cv_by_target)) {
    sol$cv_by_target
  } else {
    sqrt(.stage_cv2(cont_n, R, rates, M))
  }
  cont_ratios <- cont_cv_by_t / cv_t
  cont_binding_idx <- which(cont_ratios >= max(cont_ratios) - 1e-7)[1L]

  # Whole-unit operational design via C branch-and-bound
  box <- .stage_whole_box(rates, pop_collapsed, stage_n)
  lower_int <- box$lower
  upper_int <- box$upper
  if (any(lower_int > upper_int)) {
    bad <- which(lower_int > upper_int)[1L]
    stop(
      sprintf(
        "stage %d (%s) population capacity %g cannot deliver one expected response in whole units",
        bad, stage_nms[bad], upper_int[bad]
      ),
      call. = FALSE
    )
  }

  vb <- if (!is.null(budget)) budget - fixed_cost else NA_real_
  disc <- tryCatch(
    .stage_solve_discrete(
      R, rates, M, cv_t^2, cost_vec, lower_int, upper_int,
      budget_var = vb, min_total = min_total, n_init = cont_n,
      max_nodes = max_nodes
    ),
    svyplan_plan_unknown = function(e) {
      if (!whole_optional) stop(e)
      list(found = FALSE)
    }
  )
  if (!isTRUE(disc$found) && whole_optional) {
    operational <- NULL
  } else if (!isTRUE(disc$found)) {
    stop(
      if (!is.null(budget)) {
        "no whole-unit design fits the budget"
      } else {
        "precision targets are not achievable with whole units under the given stage constraints"
      },
      call. = FALSE
    )
  }

  labels <- if ("name" %in% names(indicators)) {
    as.character(indicators$name)
  } else {
    as.character(seq_len(nrow(indicators)))
  }
  if (isTRUE(disc$found)) {
    op_n <- as.integer(disc$n)
    names(op_n) <- stage_nms
    op_cvs <- sqrt(.stage_cv2(op_n, R, rates, M))
    operational <- list(
      n = op_n,
      total_n = prod(as.numeric(op_n)),
      cost = .stage_cost(op_n, cost_vec, fixed_cost),
      cv = max(op_cvs),
      cv_by_target = op_cvs,
      binding = labels[which.max(op_cvs / cv_t)],
      search = .stage_search_record(disc, max_nodes)
    )
  }

  # Per-row standalone stage-1 requirement at the optimal takes (for $detail$.n)
  tel <- .stage_telescoping(R, rates, M)
  q_opt <- c(1, cumprod(cont_n[2L:stages]))
  G_row <- as.numeric(tel$alpha %*% (1 / q_opt))
  row_n1 <- pmax(1, G_row / (cv_t^2 + tel$kappa))
  if (!is.null(pop_collapsed)) row_n1 <- pmin(pop_collapsed$M_take[1L], row_n1)
  row_total <- row_n1 * prod(cont_n[2L:stages])

  detail <- data.frame(
    name = labels,
    .n = row_total,
    .cv_target = cv_t,
    .cv_achieved = cont_cv_by_t,
    .binding = seq_len(nrow(indicators)) == cont_binding_idx,
    stringsAsFactors = FALSE
  )

  params <- list(
    stage_cost = stage_cost,
    mode = mode,
    domain_cols = domain_cols
  )
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(stage_n)) {
    names(stage_n) <- stage_nms
    params$stage_n <- stage_n
  }
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost
  if (!is.null(pop_collapsed)) {
    params$fpc <- "stage"
    if (record_pop) params$N <- pop_collapsed$N
  }

  res <- .new_svyplan_cluster(
    n = cont_n,
    stages = stages,
    total_n = prod(cont_n),
    cv = max(cont_cv_by_t),
    cost = sol$cost,
    params = params,
    indicators = indicators,
    detail = detail,
    operational = operational
  )
  res$binding <- labels[cont_binding_idx]
  res
}

#' Solve multistage indicator table across domains (separate or joint budget)
#' @keywords internal
#' @noRd
.n_cluster_domains <- function(
  indicators,
  stage_cost,
  budget,
  stage_n,
  domain_cols,
  joint = FALSE,
  min_n_domain = NULL,
  fixed_cost = 0,
  mode = "cv",
  domain_sampling = "separate",
  row_pops = NULL,
  max_nodes = 2e7
) {
  key <- .domain_key(indicators, domain_cols)
  domain_levels <- unique(key)
  domain_keys <- factor(key, levels = domain_levels)
  split_idx <- split(seq_len(nrow(indicators)), domain_keys)

  if (joint && !is.null(budget)) {
    return(.n_cluster_domains_joint(
      indicators, stage_cost, budget, stage_n, domain_cols,
      domain_levels, split_idx, min_n_domain, fixed_cost, mode, row_pops,
      max_nodes = max_nodes
    ))
  }

  results <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- indicators[rows, , drop = FALSE]
    sub_inner <- .label_by_row(
      sub[, !names(sub) %in% domain_cols, drop = FALSE], rows
    )
    pop_c <- if (!is.null(row_pops)) .collapse_row_pops(row_pops[rows]) else NULL
    .n_cluster_multi_solve(
      sub_inner, stage_cost, budget, stage_n, fixed_cost,
      pop_collapsed = pop_c, domain_cols = character(0), mode = mode,
      record_pop = FALSE, max_nodes = max_nodes
    )
  })
  names(results) <- domain_levels

  res <- .aggregate_cluster_domains_s(
    results, indicators, domain_cols, domain_levels, split_idx,
    stage_cost, budget, stage_n, min_n_domain, fixed_cost,
    mode = mode, domain_sampling = domain_sampling, joint = FALSE
  )
  if (!is.null(row_pops)) res$params$fpc <- "stage"
  res
}

#' Joint budget split across domains for arbitrary S >= 2
#' @keywords internal
#' @noRd
.n_cluster_domains_joint <- function(
  indicators,
  stage_cost,
  budget,
  stage_n,
  domain_cols,
  domain_levels,
  split_idx,
  min_n_domain = NULL,
  fixed_cost = 0,
  mode = "budget",
  row_pops = NULL,
  max_nodes = 2e7
) {
  stages <- length(stage_cost)
  stage_nms <- .stage_names(stages, names(stage_cost))
  cost_vec <- unname(as.numeric(stage_cost))
  nd <- length(domain_levels)
  var_budget <- budget - fixed_cost
  min_total <- min_n_domain %||% 0

  domain_label <- function(d) {
    paste(
      unlist(lapply(
        indicators[split_idx[[domain_levels[d]]][1L], domain_cols, drop = FALSE],
        as.character
      )),
      collapse = ":"
    )
  }

  # Precompute per-domain matrices and continuous/whole bounds
  dom_data <- lapply(seq_len(nd), function(d) {
    rows <- split_idx[[domain_levels[d]]]
    sub <- indicators[rows, , drop = FALSE]
    pop_c <- if (!is.null(row_pops)) .collapse_row_pops(row_pops[rows]) else NULL
    mats <- .cluster_ind_matrices(sub, stages, pop_c)
    if (!is.null(stage_n)) {
      if (!is.null(pop_c)) {
        for (s in seq_len(stages)) {
          if (!is.na(stage_n[s]) && stage_n[s] > pop_c$M_take[s] + 1e-9) {
            stop(
              sprintf(
                "domain '%s': fixed %s (%g) exceeds population capacity %g",
                domain_label(d), stage_nms[s], stage_n[s], pop_c$M_take[s]
              ),
              call. = FALSE
            )
          }
        }
      }
      .check_ind_expected_takes(sub, stage_n, stages, stage_nms)
    }
    lower <- c(1, vapply(2L:stages, function(s) .take_bound(mats$rates[, s]),
                         numeric(1L)))
    upper <- if (!is.null(pop_c)) pop_c$M_take else rep(Inf, stages)
    if (!is.null(stage_n)) {
      fixed <- !is.na(stage_n)
      lower[fixed] <- stage_n[fixed]
      upper[fixed] <- stage_n[fixed]
    }
    box <- .stage_whole_box(mats$rates, pop_c, stage_n)
    lower_int <- box$lower
    upper_int <- box$upper
    list(
      rows = rows, sub = sub, pop_c = pop_c,
      R = mats$R, rates = mats$rates, M = mats$M, cv_t = sub$cv,
      lower = lower, upper = upper,
      lower_int = lower_int, upper_int = upper_int
    )
  })

  ## The continuous and the whole design each need their own floors funded.
  cont_floor <- vapply(dom_data, function(dd) .stage_cost(dd$lower, cost_vec, 0),
                       numeric(1L))
  whole_floor <- vapply(dom_data, function(dd) {
    .stage_cost(dd$lower_int, cost_vec, 0)
  }, numeric(1L))
  floor_design <- lapply(dom_data, `[[`, "lower_int")
  refuse_sum <- function(need, whole) {
    if (sum(need) <= var_budget * (1 + 1e-12)) return(invisible(NULL))
    unit <- if (whole) " in whole units" else ""
    stop(
      if (is.null(min_n_domain)) {
        sprintf(
          "'budget' is too small for a joint allocation%s: the cheapest design of every domain together costs %.4g plus fixed_cost = %.4g",
          unit, sum(need), fixed_cost
        )
      } else {
        sprintf(
          "min_n_domain = %g not achievable for all domains within budget%s: they need %.4g together plus fixed_cost = %.4g",
          min_n_domain, unit, sum(need), fixed_cost
        )
      },
      call. = FALSE
    )
  }
  if (!is.null(min_n_domain)) {
    for (d in seq_len(nd)) {
      dd <- dom_data[[d]]
      w <- .stage_min_total_whole(min_n_domain, cost_vec, dd$lower_int,
                                  dd$upper_int, max_nodes)
      whole_c <- if (w$exact) {
        w$cost
      } else {
        .stage_min_total_cost(min_n_domain, cost_vec, dd$lower_int,
                              dd$upper_int, whole = FALSE)
      }
      cont_c <- .stage_min_total_cost(min_n_domain, cost_vec, dd$lower,
                                      dd$upper, whole = FALSE)
      if (!is.finite(cont_c) || !is.finite(whole_c)) {
        stop(
          sprintf(
            "min_n_domain = %g not achievable for domain '%s' within its population",
            min_n_domain, domain_label(d)
          ),
          call. = FALSE
        )
      }
      least <- max(whole_c, cont_c)
      if (least > var_budget * (1 + 1e-12)) {
        stop(
          sprintf(
            "min_n_domain = %g not achievable for domain '%s' within budget%s: its cheapest design costs %.4g plus fixed_cost = %.4g",
            min_n_domain, domain_label(d),
            if (whole_c > cont_c) " in whole units" else "", least, fixed_cost
          ),
          call. = FALSE
        )
      }
      cont_floor[d] <- max(cont_floor[d], cont_c)
      whole_floor[d] <- max(whole_floor[d], whole_c)
      floor_design[d] <- list(w$n)
    }
  }
  refuse_sum(cont_floor, whole = FALSE)
  refuse_sum(whole_floor, whole = TRUE)

  lower_bounds <- pmax(1e-4, cont_floor / var_budget)
  floor_tol <- 1e-8

  ratio_at <- function(d, b) {
    dd <- dom_data[[d]]
    sol <- .stage_solve_continuous_multi(
      dd$R, dd$rates, dd$M, dd$cv_t, cost_vec, dd$lower, dd$upper,
      budget = b, fixed_cost = 0, min_total = min_total
    )
    if (!sol$ok) return(Inf)
    if (!is.null(min_n_domain) && prod(sol$n) < min_n_domain * (1 - floor_tol)) {
      return(Inf)
    }
    if (!is.null(sol$ratio)) sol$ratio else sol$cv / dd$cv_t[1L]
  }

  if (sum(lower_bounds) > 1 + 1e-12) {
    stop(
      "cumulative lower bounds exceed 1; joint allocation is infeasible",
      call. = FALSE
    )
  }
  no_fixed_n1 <- is.null(stage_n) || is.na(stage_n[1L])
  budgets <- .joint_split(
    ratio_at, lower_bounds * var_budget, var_budget,
    vapply(seq_len(nd), domain_label, character(1L)),
    scaling = is.null(row_pops) && no_fixed_n1 && is.null(min_n_domain)
  )

  results <- lapply(seq_len(nd), function(d) {
    dd <- dom_data[[d]]
    sub_inner <- .label_by_row(
      dd$sub[, !names(dd$sub) %in% domain_cols, drop = FALSE], dd$rows
    )
    fit <- tryCatch(
      .n_cluster_multi_solve(
        sub_inner, stage_cost, budgets[d], stage_n, 0,
        pop_collapsed = dd$pop_c, domain_cols = character(0), mode = mode,
        record_pop = FALSE, min_total = min_total, max_nodes = max_nodes,
        whole_optional = TRUE
      ),
      error = function(e) {
        e$message <- sprintf("domain '%s': %s", domain_label(d),
                             conditionMessage(e))
        e$call <- NULL
        stop(e)
      }
    )
    fit
  })
  names(results) <- domain_levels

  shares_n <- lapply(results, function(f) f$operational$n)
  witness <- if (!any(vapply(shares_n, is.null, logical(1L)))) {
    lapply(shares_n, as.numeric)
  } else if (!any(vapply(floor_design, is.null, logical(1L))) &&
               sum(vapply(floor_design, .stage_cost, 1, cost_vec, 0)) <=
                 .budget_room(var_budget)) {
    floor_design
  }
  if (is.null(witness)) {
    stop(structure(
      list(message = sprintf(
        "no whole design within the joint budget was found in %s search nodes, so whether one exists is unknown. Raise 'control$max_nodes'",
        format(max_nodes, big.mark = ",", scientific = FALSE)
      ), call = NULL),
      class = c("svyplan_plan_unknown", "error", "condition")
    ))
  }
  mm <- .joint_whole_minimax(dom_data, cost_vec, var_budget, min_total,
                             witness, max_nodes)
  n_whole <- mm$n
  used <- mm$nodes_by_domain
  ratio_of <- function(d, n) {
    dd <- dom_data[[d]]
    max(sqrt(.stage_cv2(n, dd$R, dd$rates, dd$M)) / dd$cv_t)
  }
  if (!is.null(mm$discs)) {
    spare <- var_budget -
      sum(vapply(n_whole, function(n) .stage_cost(n, cost_vec, 0), 1))
    share <- budgets / sum(budgets)
    for (d in seq_len(nd)) {
      left <- max_nodes - sum(used)
      if (left < 1) break
      dd <- dom_data[[d]]
      b_d <- .stage_cost(n_whole[[d]], cost_vec, 0) + spare * share[d]
      up <- tryCatch(
        .stage_solve_discrete(
          dd$R, dd$rates, dd$M, dd$cv_t^2, cost_vec, dd$lower_int,
          dd$upper_int, budget_var = b_d, min_total = min_total,
          n_init = n_whole[[d]], max_nodes = left
        ),
        svyplan_plan_unknown = function(e) list(nodes = e$nodes %||% left)
      )
      used[d] <- used[d] + up$nodes
      if (isTRUE(up$found) &&
            .stage_cost(up$n, cost_vec, 0) <= b_d * (1 + 1e-12) &&
            ratio_of(d, up$n) <= mm$ratio * (1 + 1e-12)) {
        n_whole[[d]] <- as.numeric(up$n)
      }
    }
  }
  ratio <- max(vapply(seq_len(nd), function(d) ratio_of(d, n_whole[[d]]), 1))
  for (d in seq_len(nd)) {
    dd <- dom_data[[d]]
    op <- results[[d]]$operational %||% list()
    n_d <- as.integer(n_whole[[d]])
    names(n_d) <- names(results[[d]]$n)
    cvs <- sqrt(.stage_cv2(n_d, dd$R, dd$rates, dd$M))
    op$n <- n_d
    op$total_n <- prod(as.numeric(n_d))
    op$cost <- .stage_cost(n_d, cost_vec, 0)
    op$cv <- max(cvs)
    op$cv_by_target <- cvs
    op$binding <- results[[d]]$detail$name[which.max(cvs / dd$cv_t)]
    op$search <- list(status = mm$status, nodes = used[d],
                      max_nodes = max_nodes)
    results[[d]]$operational <- op
  }

  res <- .aggregate_cluster_domains_s(
    results, indicators, domain_cols, domain_levels, split_idx,
    stage_cost, budget, stage_n, min_n_domain, fixed_cost,
    mode = mode, domain_sampling = "separate", joint = TRUE
  )
  per_domain <- res$operational$search
  res$operational$search <- list(
    status = mm$status, nodes = sum(used), max_nodes = max_nodes,
    by_domain = per_domain$by_domain,
    nodes_by_domain = per_domain$nodes_by_domain, bound = mm$bound,
    gap = if (ratio > 0) max(0, (ratio - mm$bound) / ratio) else 0,
    ratio = ratio, trials = mm$trials,
    shares = vapply(n_whole, function(n) .stage_cost(n, cost_vec, 0), 1)
  )
  names(res$operational$search$shares) <- domain_levels
  res$params$joint <- TRUE
  if (!is.null(row_pops)) {
    res$params$fpc <- "stage"
    res$cost <- fixed_cost + sum(res$domains$.cost)
  } else if (fixed_cost > 0) {
    res$cost <- budget
  }
  res$operational$cost <- fixed_cost + sum(res$operational$domains$.cost)
  res
}

#' Aggregate per-domain S-stage cluster results into one svyplan_cluster object
#' @keywords internal
#' @noRd
.aggregate_cluster_domains_s <- function(
  results,
  indicators,
  domain_cols,
  domain_levels,
  split_idx,
  stage_cost,
  budget,
  stage_n,
  min_n_domain,
  fixed_cost,
  mode = "cv",
  domain_sampling = "separate",
  joint = FALSE
) {
  stages <- length(stage_cost)
  stage_nms <- .stage_names(stages, names(stage_cost))
  nd <- length(domain_levels)

  dom_rows <- vector("list", nd)
  op_rows <- vector("list", nd)
  searches <- vector("list", nd)
  search_labels <- character(nd)
  detail_list <- vector("list", nd)
  cont_ratio <- numeric(nd)
  op_ratio <- numeric(nd)

  for (d in seq_len(nd)) {
    lev <- domain_levels[d]
    rows <- split_idx[[lev]]
    fit <- results[[lev]]
    dom_vals <- indicators[rows[1L], domain_cols, drop = FALSE]

    row_df <- dom_vals
    op_df <- dom_vals
    for (s in seq_len(stages)) {
      row_df[[stage_nms[s]]] <- unname(fit$n[s])
      op_df[[stage_nms[s]]] <- unname(fit$operational$n[s])
    }
    row_df$.total_n <- fit$total_n
    row_df$.cv <- fit$cv
    row_df$.cost <- fit$cost
    bind_name <- fit$detail$name[which(fit$detail$.binding)[1L]] %||%
      as.character(which(fit$detail$.binding)[1L])
    row_df$.binding <- bind_name
    dom_rows[[d]] <- row_df

    op_df$.total_n <- fit$operational$total_n
    op_df$.cv <- fit$operational$cv
    op_df$.cost <- fit$operational$cost
    op_df$.binding <- fit$operational$binding
    op_rows[[d]] <- op_df
    searches[[d]] <- fit$operational$search
    search_labels[d] <- paste(unlist(dom_vals), collapse = "/")

    det <- fit$detail
    cont_ratio[d] <- max(det$.cv_achieved / det$.cv_target)
    op_ratio[d] <- max(fit$operational$cv_by_target / det$.cv_target)
    det <- cbind(indicators[rows, domain_cols, drop = FALSE], det)
    detail_list[[d]] <- det
  }

  dom_df <- do.call(rbind, dom_rows)
  rownames(dom_df) <- NULL
  op_dom_df <- do.call(rbind, op_rows)
  rownames(op_dom_df) <- NULL
  detail_df <- do.call(rbind, detail_list)
  rownames(detail_df) <- NULL

  binding_cause <- "precision"
  if (!is.null(min_n_domain)) {
    if (joint) {
      if (any(dom_df$.total_n <= min_n_domain * (1 + 1e-6))) {
        binding_cause <- "min_n_domain"
      }
    } else {
      below <- which(dom_df$.total_n < min_n_domain * (1 - 1e-8))
      if (length(below) > 0L) {
        warning(
          sprintf(
            "%d domain(s) fell below min_n_domain = %g; increase targets or use allocation = \"joint\" in budget mode",
            length(below), min_n_domain
          ),
          call. = FALSE
        )
      }
    }
  }

  na_n <- setNames(rep(NA_real_, stages), stage_nms)
  na_op_n <- setNames(rep(NA_integer_, stages), stage_nms)
  operational <- list(
    n = na_op_n,
    total_n = sum(op_dom_df$.total_n),
    cost = sum(op_dom_df$.cost),
    cv = max(op_dom_df$.cv),
    binding = op_dom_df$.binding[which.max(op_ratio)],
    domains = op_dom_df,
    search = .stage_search_domains(searches, search_labels)
  )

  params <- list(
    stage_cost = stage_cost,
    mode = mode,
    domain_cols = domain_cols,
    domain_sampling = domain_sampling,
    allocation = if (joint) "joint" else "separate",
    joint = joint
  )
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(stage_n)) {
    names(stage_n) <- stage_nms
    params$stage_n <- stage_n
  }
  if (!is.null(min_n_domain)) params$min_n_domain <- min_n_domain
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost

  res <- .new_svyplan_cluster(
    n = na_n,
    stages = stages,
    total_n = sum(dom_df$.total_n),
    cv = max(dom_df$.cv),
    cost = sum(dom_df$.cost),
    params = params,
    indicators = indicators,
    detail = detail_df,
    domains = dom_df,
    operational = operational
  )
  res$binding <- dom_df$.binding[which.max(cont_ratio)]
  res$binding_cause <- binding_cause
  res
}

#' @rdname prec_cluster
#' @name prec_cluster
#' @aliases prec_cluster.data.frame
#' @export
prec_cluster.data.frame <- function(
  x,
  ...,
  domains = NULL,
  stage_cost = NULL,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  N = NULL,
  fpc = c("auto", "none", "stage"),
  plan = NULL
) {
  dots <- list(...)
  .check_prec_cluster_indicator_args(
    n = dots[["n"]],
    icc = dots[["icc"]],
    unit_relvar = dots[["unit_relvar"]],
    var_ratio = dots[["var_ratio"]]
  )
  fpc <- match.arg(fpc)
  sc_len <- if (is.numeric(stage_cost)) length(stage_cost) else NULL
  plan <- .plan_for_cluster_mode(plan, TRUE, fpc, stages = sc_len)
  .plan <- .merge_plan_args(plan, prec_cluster.data.frame, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_cluster.data.frame, c(.plan, dots)))
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(x, "x")

  x <- as.data.frame(x)
  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)

  if (!is.data.frame(x) || nrow(x) == 0L) {
    stop("'x' must be a non-empty data frame", call. = FALSE)
  }

  # Detect stages from take_2..take_S columns (and stage_cost if supplied)
  take_cols <- grep("^take_[0-9]+$", names(x), value = TRUE)
  if (!"n_1" %in% names(x) || length(take_cols) == 0L) {
    stop(
      "prec_cluster() on an indicator table requires 'n_1' and stage take columns 'take_2', ..., 'take_S'",
      call. = FALSE
    )
  }
  take_nums <- sort(as.integer(sub("^take_", "", take_cols)))
  stages <- max(take_nums)
  implied <- max(stages, if (is.numeric(stage_cost)) length(stage_cost))
  if (implied > stages) {
    stop(
      sprintf(
        "a %d-stage table needs take column(s) %s",
        implied, paste(sQuote(paste0("take_", (stages + 1L):implied)),
                       collapse = ", ")
      ),
      call. = FALSE
    )
  }
  expected_takes <- seq.int(2L, stages)
  if (!identical(take_nums, expected_takes)) {
    stop(
      sprintf(
        "take columns must be contiguous from 'take_2' to 'take_%d'",
        stages
      ),
      call. = FALSE
    )
  }

  if (!is.null(stage_cost)) {
    check_stage_cost(stage_cost)
    stage_cost <- .reorder_stage_cost(stage_cost)
    if (length(stage_cost) != stages) {
      stop(sprintf("'stage_cost' must have length %d", stages), call. = FALSE)
    }
    names(stage_cost) <- .stage_names(stages, names(stage_cost))
  }

  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(x))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in indicators: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  domain_cols <- domains %||% character(0)

  .prec_multi_cluster(
    indicators = x,
    stages = stages,
    stage_cost = stage_cost,
    domain_cols = domain_cols,
    fpc = fpc,
    N_arg = N,
    resp_rate = resp_rate,
    eligibility_rate = eligibility_rate,
    eligible_per_unit = eligible_per_unit
  )
}

#' Multi-indicator precision for arbitrary S >= 2 stages
#' @keywords internal
#' @noRd
.prec_multi_cluster <- function(
  indicators,
  stages,
  stage_cost = NULL,
  domain_cols = character(0),
  fpc = "auto",
  N_arg = NULL,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1
) {
  indicators <- .indicators_var_from_sd(indicators, domain_cols)
  .check_indicator_columns(indicators, domain_cols, eligible = FALSE,
                           stages = stages)
  .check_resp_rate_column(indicators, TRUE, stages)
  .stop_min_cases_column(indicators, "is not supported by prec_cluster()")

  # Validate stage sizes
  size_cols <- c("n_1", paste0("take_", 2L:stages))
  for (col in size_cols) {
    vals <- indicators[[col]]
    if (!is.numeric(vals) || anyNA(vals) || any(!is.finite(vals)) || any(vals <= 0)) {
      stop(sprintf("'%s' values must be positive, finite, and non-NA", col),
           call. = FALSE)
    }
  }

  # Validate estimands
  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)
  .check_ratio_rows(indicators, sizing = FALSE)
  if (!has_p && !has_var && !any(is_ratio)) {
    stop("'indicators' must contain a 'p', 'var', or 'r' column", call. = FALSE)
  }
  .check_estimand_markers(indicators)
  if (has_p) {
    p_vals <- indicators$p[!is.na(indicators$p)]
    if (any(p_vals <= 0 | p_vals >= 1)) {
      stop("all 'p' values must be in (0, 1)", call. = FALSE)
    }
  }
  if (has_var) {
    var_vals <- indicators$var[!is.na(indicators$var)]
    if (any(var_vals <= 0) || any(!is.finite(var_vals))) {
      stop("all 'var' values must be positive and finite", call. = FALSE)
    }
    if (!"mean" %in% names(indicators) ||
        any(is.na(indicators$mean[!is.na(indicators$var)]))) {
      stop("'mean' is required for mean indicators in prec_cluster()",
           call. = FALSE)
    }
    if (any(indicators$mean[!is.na(indicators$var)] == 0)) {
      stop("'mean' must be non-zero for mean indicators in prec_cluster()",
           call. = FALSE)
    }
  }

  # Validate icc_1..icc_{S-1}
  n_b <- stages - 1L
  req_icc <- paste0("icc_", seq_len(n_b))
  missing_icc <- setdiff(req_icc, names(indicators))
  if (length(missing_icc) > 0L) {
    stop(
      sprintf(
        "%d-stage prec_cluster() requires column(s) %s",
        stages, paste(sQuote(missing_icc), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  for (col in req_icc) {
    vals <- indicators[[col]]
    if (!is.numeric(vals) || anyNA(vals) || any(!is.finite(vals)) ||
        any(vals < 0 | vals > 1)) {
      stop(sprintf("'%s' values must be finite numbers in [0, 1]", col),
           call. = FALSE)
    }
  }

  .validate_common_columns(indicators)
  resp_vec <- .expand_stage_resp_rate(resp_rate, stages)
  indicators <- .fill_defaults(indicators, TRUE, stages = stages,
                               resp_rate = resp_vec)
  indicators <- .fill_eligible_per_unit(indicators, eligible_per_unit)

  ## A constant ratio has zero relvariance, and prec_ratio() reports its
  ## precision as se = 0. On any other row a zero contradicts the inputs.
  rv_vals <- indicators$unit_relvar
  if (!is.numeric(rv_vals) || anyNA(rv_vals) ||
      any(!is.finite(rv_vals)) || any(rv_vals < 0) ||
      any(rv_vals == 0 & !.is_ratio_row(indicators))) {
    stop("'unit_relvar' values must be positive and finite", call. = FALSE)
  }
  for (s in seq_len(n_b)) {
    col <- paste0("var_ratio_", s)
    vals <- indicators[[col]]
    if (!is.numeric(vals) || anyNA(vals) ||
        any(!is.finite(vals)) || any(vals <= 0)) {
      stop(sprintf("'%s' values must be positive and finite", col),
           call. = FALSE)
    }
  }

  given_resp <- indicators$resp_rate
  epu_vec <- .indicator_rate(indicators, "eligible_per_unit", nrow(indicators))
  yield_vec <- eligibility_rate * epu_vec
  indicators$resp_rate <- given_resp * yield_vec

  row_pops <- .stage_row_pops(indicators, domain_cols, fpc, N_arg, stages,
                              yield_vec)
  nr <- nrow(indicators)
  rate_cols <- c(paste0("resp_rate_", seq_len(n_b)), "resp_rate")
  labels <- if ("name" %in% names(indicators)) indicators$name else seq_len(nr)

  cv_vals <- numeric(nr)
  se_vals <- numeric(nr)
  moe_vals <- numeric(nr)
  rmoe_vals <- numeric(nr)

  for (i in seq_len(nr)) {
    n_i <- vapply(size_cols, function(col) indicators[[col]][i], numeric(1L))
    rates_i <- vapply(rate_cols, function(col) indicators[[col]][i], numeric(1L))
    for (s in 2L:stages) {
      .check_expected_take(
        n_i[s], rates_i[s], size_cols[s], rate_cols[s],
        sprintf(" for indicator %s", sQuote(labels[i]))
      )
    }
    if (!is.null(row_pops)) {
      .check_stage_sizes(n_i, row_pops[[i]], size_cols)
    }
    icc_i <- vapply(req_icc, function(col) indicators[[col]][i], numeric(1L))
    vr_i <- vapply(paste0("var_ratio_", seq_len(n_b)),
                   function(col) indicators[[col]][i], numeric(1L))
    R_i <- .stage_R_matrix(icc_i, indicators$unit_relvar[i], vr_i)
    M_i <- matrix(if (!is.null(row_pops)) row_pops[[i]]$M else rep(Inf, stages),
                  nrow = 1L)
    cv_i <- sqrt(.stage_cv2(n_i, R_i, matrix(rates_i, nrow = 1L), M_i)[1L])
    cv_vals[i] <- cv_i

    alpha_i <- indicators$alpha[i]
    df_i <- .row_df(indicators, i)
    z_i <- .q_alpha(alpha_i, df_i)
    if (is_ratio[i]) {
      scale_i <- abs(indicators$r[i])
      se_vals[i] <- cv_i * scale_i
      moe_vals[i] <- z_i * se_vals[i]
      rmoe_vals[i] <- z_i * cv_i
    } else if (has_p && !is.na(indicators$p[i])) {
      p_i <- indicators$p[i]
      se_vals[i] <- cv_i * p_i
      pm_i <- indicators$prop_method[i]
      if (identical(pm_i, "wald") || cv_i == 0) {
        moe_vals[i] <- z_i * se_vals[i]
      } else {
        n_eff_i <- (1 - p_i) / (p_i * cv_i^2)
        moe_vals[i] <- prec_prop(p = p_i, n = n_eff_i, alpha = alpha_i,
                                 method = pm_i, df = df_i)$moe
      }
      rmoe_vals[i] <- moe_vals[i] / p_i
    } else {
      scale_i <- abs(indicators$mean[i])
      se_vals[i] <- cv_i * scale_i
      moe_vals[i] <- z_i * se_vals[i]
      rmoe_vals[i] <- z_i * cv_i
    }
  }

  indicators$resp_rate <- given_resp
  detail <- indicators
  if (length(domain_cols) > 0L) {
    other_cols <- setdiff(names(detail), domain_cols)
    detail <- detail[, c(domain_cols, other_cols), drop = FALSE]
  }
  detail$.se <- se_vals
  detail$.moe <- moe_vals
  detail$.rmoe <- rmoe_vals
  detail$.cv <- cv_vals

  params <- list(
    indicators = indicators,
    stages = stages,
    domain_cols = domain_cols,
    resp_rate = .name_stage_resp_rate(resp_vec)
  )
  if (!is.null(stage_cost)) params$stage_cost <- stage_cost
  if (eligibility_rate != 1) params$eligibility_rate <- eligibility_rate
  if (eligible_per_unit != 1 || "eligible_per_unit" %in% names(indicators)) {
    params$eligible_per_unit <- eligible_per_unit
  }
  if (!is.null(row_pops)) {
    params$fpc <- "stage"
    if (!is.null(N_arg)) params$N <- N_arg
  }

  .new_svyplan_prec(
    se = se_vals,
    moe = moe_vals,
    cv = cv_vals,
    type = "cluster",
    params = params,
    detail = detail
  )
}
