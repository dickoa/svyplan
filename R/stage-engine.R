#' Default or caller-supplied stage names for an S-stage design
#'
#' Preserves valid unique names supplied on `stage_cost` (or `n`), falling
#' back to `stage_1`, ..., `stage_S` when unnamed.
#' @keywords internal
#' @noRd
.stage_names <- function(stages, nms = NULL) {
  default <- paste0("stage_", seq_len(stages))
  if (is.null(nms) || length(nms) != stages) {
    return(default)
  }
  if (anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    return(default)
  }
  nms
}

#' Expand scalar or length-(S-1) var_ratio to length S - 1
#' @keywords internal
#' @noRd
.expand_var_ratio <- function(var_ratio, icc) {
  n_icc <- length(icc)
  if (length(var_ratio) == 1L) {
    var_ratio <- c(var_ratio, rep(NA_real_, n_icc - 1L))
  }
  if (length(var_ratio) != n_icc) {
    stop(
      sprintf(
        "'var_ratio' must have length 1 or %d (matching 'icc')",
        n_icc
      ),
      call. = FALSE
    )
  }
  out <- unname(as.numeric(var_ratio))
  if (is.na(out[1L])) out[1L] <- 1
  for (s in seq_len(n_icc)[-1L]) {
    if (!is.na(out[s])) next
    out[s] <- out[s - 1L] * (1 - icc[s - 1L])
    if (out[s] <= 0) {
      stop(
        sprintf(
          "'icc[%d]' of 1 leaves no within-stage variance, so deeper stages are not identified; supply 'var_ratio' explicitly or reduce the number of stages",
          s - 1L
        ),
        call. = FALSE
      )
    }
  }
  out
}

#' Put a named stage vector in stage order
#' @keywords internal
#' @noRd
.match_stage_vec <- function(x, name, stages, stage_nms = NULL,
                             keys = list()) {
  nms <- names(x)
  if (is.null(nms) || length(x) == 1L) return(unname(x))
  if (length(x) != stages) return(unname(x))
  stage_nms <- stage_nms %||% paste0("stage_", seq_len(stages))
  for (key in c(list(stage_nms, paste0("stage_", seq_len(stages))), keys)) {
    if (!anyNA(nms) && !anyDuplicated(nms) && setequal(nms, key)) {
      return(unname(x[match(key, nms)]))
    }
  }
  stop(
    sprintf(
      "names of '%s' must match the stages (%s) or be absent",
      name, paste(sQuote(stage_nms), collapse = ", ")
    ),
    call. = FALSE
  )
}

#' Check a user's var_ratio, `NA` allowed past the first entry
#' @keywords internal
#' @noRd
.check_var_ratio <- function(var_ratio) {
  if (!is.numeric(var_ratio) || length(var_ratio) == 0L ||
        is.na(var_ratio[1L]) || any(var_ratio <= 0, na.rm = TRUE) ||
        any(is.infinite(var_ratio))) {
    stop("'var_ratio' must contain positive finite values, with NA only ",
         "past the first, where it is derived", call. = FALSE)
  }
}

#' The var_ratio values a user gave, `NA` where the stage is derived
#' @keywords internal
#' @noRd
.var_ratio_given <- function(var_ratio, n_icc) {
  if (length(var_ratio) == 1L) {
    return(c(unname(var_ratio), rep(NA_real_, n_icc - 1L)))
  }
  unname(as.numeric(var_ratio))
}

#' Expand scalar or length-S stage response rates
#'
#' A scalar applies to the final stage (ultimate units) with earlier stages
#' set to 1. A length-S vector supplies stage-1 through stage-S rates.
#' @keywords internal
#' @noRd
.expand_stage_resp_rate <- function(resp_rate, stages, stage_nms = NULL) {
  resp_rate <- .match_stage_vec(
    resp_rate, "resp_rate", stages, stage_nms,
    list(c(paste0("resp_rate_", seq_len(stages - 1L)), "resp_rate"),
         paste0("resp_rate_", seq_len(stages)))
  )
  if (!is.numeric(resp_rate) || length(resp_rate) == 0L ||
      anyNA(resp_rate) || any(!is.finite(resp_rate)) ||
      any(resp_rate <= 0) || any(resp_rate > 1)) {
    stop("'resp_rate' must contain values in (0, 1]", call. = FALSE)
  }
  if (length(resp_rate) == 1L) {
    out <- rep(1, stages)
    out[stages] <- unname(resp_rate)
    return(out)
  }
  if (length(resp_rate) != stages) {
    stop(
      sprintf(
        "'resp_rate' must be a scalar (final stage) or a vector of length %d (one per stage)",
        stages
      ),
      call. = FALSE
    )
  }
  unname(as.numeric(resp_rate))
}

#' Stage response rates as stored in `params`, named as the allocations name
#' their columns: `resp_rate_1`, ..., `resp_rate_{S-1}` and `resp_rate`
#' @keywords internal
#' @noRd
.name_stage_resp_rate <- function(rate) {
  stages <- length(rate)
  names(rate) <- c(if (stages > 1L) paste0("resp_rate_", seq_len(stages - 1L)),
                   "resp_rate")
  rate
}

#' Resolve fpc = c("auto", "none", "stage") against population vector N
#' @keywords internal
#' @noRd
.resolve_stage_fpc <- function(fpc, N, has_pop_cols = FALSE) {
  fpc <- match.arg(fpc, c("auto", "none", "stage"))
  has_pop <- !is.null(N) || has_pop_cols
  if (identical(fpc, "auto")) {
    return(if (has_pop) "stage" else "none")
  }
  if (identical(fpc, "none") && has_pop) {
    stop(
      "'N' is supplied for finite population correction, but fpc = \"none\" disables it. Use fpc = \"auto\" (the default) or fpc = \"stage\", or drop 'N'",
      call. = FALSE
    )
  }
  if (identical(fpc, "stage") && !has_pop) {
    stop(
      "fpc = \"stage\" requires stage population counts 'N' (length equal to the number of stages)",
      call. = FALSE
    )
  }
  fpc
}

#' Validate and build an S-stage population structure from aggregate counts N
#'
#' `N` is a length-S vector of aggregate population counts ordered from
#' outermost stage 1 to innermost stage S (ultimate units), so
#' `N[1] <= N[2] <= ... <= N[S]`. Per-parent capacities are `M_take[1] = N[1]`
#' and `M_take[s] = N[s] / N[s-1]` for `s >= 2`. For variance evaluation, the
#' final-stage capacity scales by `eligibility_rate * eligible_per_unit`.
#' @keywords internal
#' @noRd
.stage_population <- function(N, stages, yield = 1, stage_nms = NULL) {
  if (is.null(N)) return(NULL)
  N <- .match_stage_vec(
    N, "N", stages, stage_nms,
    list(c(paste0("N_", seq_len(stages - 1L)), "N"),
         paste0("N_", seq_len(stages)))
  )
  if (!is.numeric(N) || length(N) != stages || anyNA(N) ||
      any(!is.finite(N)) || any(N <= 0)) {
    stop(
      sprintf(
        "'N' must be a positive finite numeric vector of length %d (aggregate population counts from stage 1 to stage %d)",
        stages, stages
      ),
      call. = FALSE
    )
  }
  if (any(abs(N - round(N)) > 1e-8)) {
    stop("'N' must contain whole numbers", call. = FALSE)
  }
  N <- unname(round(as.numeric(N)))
  if (any(diff(N) < 0)) {
    stop(
      "'N' must be non-decreasing from stage 1 to stage S (at least one sub-unit per parent at each stage)",
      call. = FALSE
    )
  }
  M_take <- numeric(stages)
  M_take[1L] <- N[1L]
  if (stages > 1L) {
    M_take[2L:stages] <- N[2L:stages] / N[seq_len(stages - 1L)]
  }
  M <- M_take
  M[stages] <- M_take[stages] * yield
  list(N = N, M_take = M_take, M = M)
}

#' Stage variance components R_1, ..., R_S for single or multiple indicators
#'
#' `icc` and `var_ratio` may be numeric vectors of length S - 1 (single
#' indicator) or J x (S - 1) matrices (J indicators). Returns a J x S matrix
#' where `R[j, s] = V_j * k_{j,s} * delta_{j,s}` for `s < S` and
#' `R[j, S] = V_j * k_{j,S-1} * (1 - delta_{j,S-1})`.
#' @keywords internal
#' @noRd
.stage_R_matrix <- function(icc, unit_relvar, var_ratio) {
  if (!is.matrix(icc)) {
    icc <- matrix(as.numeric(icc), nrow = 1L)
  }
  J <- nrow(icc)
  S <- ncol(icc) + 1L
  if (!is.matrix(var_ratio)) {
    if (length(var_ratio) == 1L) {
      vr_mat <- matrix(0, nrow = J, ncol = S - 1L)
      for (j in seq_len(J)) {
        vr_mat[j, ] <- .expand_var_ratio(var_ratio, icc[j, ])
      }
      var_ratio <- vr_mat
    } else if (length(var_ratio) == S - 1L) {
      var_ratio <- matrix(rep(as.numeric(var_ratio), each = J), nrow = J)
    } else {
      stop("incompatible 'var_ratio' dimensions", call. = FALSE)
    }
  }
  V <- rep_len(as.numeric(unit_relvar), J)
  R <- matrix(0, nrow = J, ncol = S)
  for (s in seq_len(S - 1L)) {
    R[, s] <- V * var_ratio[, s] * icc[, s]
  }
  R[, S] <- V * var_ratio[, S - 1L] * (1 - icc[, S - 1L])
  R
}

#' Telescoping stage-FPC coefficients alpha (J x S) and kappa (J)
#'
#' Given J x S matrices `R`, `rates`, and `M` (`Inf` where no FPC applies),
#' computes exact telescoping coefficients satisfying
#' `cv2_j(n) = sum_s alpha[j, s] / prod(n[1:s]) - kappa[j]`.
#' @keywords internal
#' @noRd
.stage_telescoping <- function(R, rates, M) {
  J <- nrow(R)
  S <- ncol(R)
  rbar <- matrix(0, nrow = J, ncol = S)
  rbar[, 1L] <- rates[, 1L]
  for (s in 2L:S) {
    rbar[, s] <- rbar[, s - 1L] * rates[, s]
  }
  alpha <- matrix(0, nrow = J, ncol = S)
  for (s in seq_len(S - 1L)) {
    sub <- ifelse(is.finite(M[, s + 1L]), R[, s + 1L] / M[, s + 1L], 0)
    alpha[, s] <- (R[, s] - sub) / rbar[, s]
  }
  alpha[, S] <- R[, S] / rbar[, S]
  kappa <- ifelse(is.finite(M[, 1L]), R[, 1L] / M[, 1L], 0)
  list(alpha = alpha, kappa = kappa, rbar = rbar)
}

#' Exact squared CV vector for an S-stage allocation n (length S)
#'
#' Uses the positive stage-by-stage FPC form so near-census allocations suffer
#' no floating-point cancellation.
#' @keywords internal
#' @noRd
.stage_cv2 <- function(n, R, rates, M) {
  n <- unname(as.numeric(n))
  J <- nrow(R)
  S <- ncol(R)
  out <- numeric(J)
  n_eff_cum <- rep(1, J)
  for (s in seq_len(S)) {
    n_eff_s <- n[s] * rates[, s]
    n_eff_cum <- n_eff_cum * n_eff_s
    fpc_s <- ifelse(is.finite(M[, s]), pmax(0, 1 - n_eff_s / M[, s]), 1)
    out <- out + R[, s] * fpc_s / n_eff_cum
  }
  pmax(out, 0)
}

#' Total cost of an S-stage design
#' @keywords internal
#' @noRd
.stage_cost <- function(n, stage_cost, fixed_cost = 0) {
  fixed_cost + sum(stage_cost * cumprod(n))
}

#' General S-stage clustering design effect
#'
#' Evaluates `sum_{s=1}^{S-1} k_s * delta_s * prod_{j=s+1}^S take_{j-1} +
#' k_{S-1} * (1 - delta_{S-1})`, where `take` has length `S - 1` (stages 2..S).
#' @keywords internal
#' @noRd
.stage_deff_cluster <- function(icc, take, var_ratio = 1) {
  S_minus_1 <- length(icc)
  k <- .expand_var_ratio(var_ratio, icc)
  suf <- numeric(S_minus_1)
  p <- 1
  for (s in rev(seq_len(S_minus_1))) {
    p <- p * take[s]
    suf[s] <- p
  }
  sum(k * icc * suf) + k[S_minus_1] * (1 - icc[S_minus_1])
}

#' Minimize G(n_2..n_K) * H(n_2..n_K) on box `(lo, hi)` for a single indicator
#'
#' `alpha` and `cost` have length K >= 1. `lo` and `hi` have length K - 1
#' (for takes n_2..n_K). Handles non-positive intermediate alpha_s exactly via
#' active-face elimination, then solves the remaining strictly log-convex
#' problem by closed-form Cauchy-Schwarz coordinate descent.
#' @keywords internal
#' @noRd
.stage_min_gh <- function(alpha, cost, lo, hi) {
  K <- length(alpha)
  if (K == 1L) {
    return(list(takes = numeric(0), gh = alpha[1L] * cost[1L]))
  }
  neg <- which(alpha[seq_len(K - 1L)] <= 0)
  if (length(neg) > 0L) {
    s <- neg[1L]
    cands <- list()
    u <- hi[s]
    if (is.finite(u)) {
      a1 <- alpha[-(s + 1L)]
      c1 <- cost[-(s + 1L)]
      a1[s] <- alpha[s] + alpha[s + 1L] / u
      c1[s] <- cost[s] + cost[s + 1L] * u
      if (s + 1L < K) {
        idx_old <- (s + 2L):K
        idx_new <- (s + 1L):(K - 1L)
        a1[idx_new] <- alpha[idx_old] / u
        c1[idx_new] <- cost[idx_old] * u
      }
      sub1 <- .stage_min_gh(a1, c1, lo[-s], hi[-s])
      t1 <- append(sub1$takes, u, after = s - 1L)
      cands[[length(cands) + 1L]] <- list(takes = t1, gh = sub1$gh)
    }
    if (s >= 2L) {
      l <- lo[s - 1L]
      a2 <- alpha[-s]
      c2 <- cost[-s]
      a2[s - 1L] <- alpha[s - 1L] + alpha[s] / l
      c2[s - 1L] <- cost[s - 1L] + cost[s] * l
      if (s < K) {
        idx_old <- (s + 1L):K
        idx_new <- s:(K - 1L)
        a2[idx_new] <- alpha[idx_old] / l
        c2[idx_new] <- cost[idx_old] * l
      }
      sub2 <- .stage_min_gh(a2, c2, lo[-(s - 1L)], hi[-(s - 1L)])
      t2 <- append(sub2$takes, l, after = s - 2L)
      cands[[length(cands) + 1L]] <- list(takes = t2, gh = sub2$gh)
    }
    best_i <- 1L
    best_val <- Inf
    for (i in seq_along(cands)) {
      q <- c(1, cumprod(cands[[i]]$takes))
      val <- sum(alpha / q) * sum(cost * q)
      if (val < best_val) {
        best_val <- val
        best_i <- i
      }
    }
    return(list(takes = cands[[best_i]]$takes, gh = best_val))
  }

  takes <- sqrt(cost[seq_len(K - 1L)] * alpha[2L:K] /
                  (cost[2L:K] * alpha[seq_len(K - 1L)]))
  takes <- pmin(hi, pmax(lo, takes))

  free <- which(hi > lo * (1 + 1e-14))
  if (length(free) > 0L && (length(free) < K - 1L ||
      any(takes[free] <= lo[free] * (1 + 1e-12) |
          takes[free] >= hi[free] * (1 - 1e-12)))) {
    for (iter in seq_len(200L)) {
      max_rel <- 0
      for (idx in free) {
        k <- idx + 1L # stage index in 2..K
        q <- c(1, cumprod(takes))
        A_k <- sum(alpha[seq_len(k - 1L)] / q[seq_len(k - 1L)])
        B_k <- takes[idx] * sum(alpha[k:K] / q[k:K])
        U_k <- sum(cost[seq_len(k - 1L)] * q[seq_len(k - 1L)])
        V_k <- sum(cost[k:K] * q[k:K]) / takes[idx]
        t_new <- if (A_k > 0 && B_k > 0) {
          min(hi[idx], max(lo[idx], sqrt(B_k * U_k / (A_k * V_k))))
        } else if (A_k <= 0) {
          hi[idx]
        } else {
          lo[idx]
        }
        rel <- abs(t_new - takes[idx]) / max(1, takes[idx])
        if (rel > max_rel) max_rel <- rel
        takes[idx] <- t_new
      }
      if (max_rel < 1e-13) break
    }
  }
  q <- c(1, cumprod(takes))
  list(takes = takes, gh = sum(alpha / q) * sum(cost * q))
}

#' Recursive single-indicator continuous S-stage solver
#'
#' Solves the exact continuous problem for K stages with cost vector `cost`,
#' telescoping coefficients `alpha`, lower/upper bounds `lo`/`hi` (length K),
#' and either a budget limit `vb` (`sum(cost * cumprod(n)) <= vb`) or a
#' variance target `g_target` (`sum(alpha / cumprod(n)) <= g_target`).
#' Returns a list with `ok`, `n` (length K), `reason` (if `!ok`), and
#' `open_takes` (unclamped stage-optimal takes for clamping warnings).
#' @keywords internal
#' @noRd
.stage_solve_one_rec <- function(alpha, cost, lo, hi, g_target = NULL,
                                 vb = NULL) {
  K <- length(alpha)
  eval_G <- function(n) sum(alpha / cumprod(n))
  eval_H <- function(n) sum(cost * cumprod(n))

  if (K == 1L) {
    if (!is.null(vb)) {
      if (cost[1L] * lo[1L] > vb * (1 + 1e-12)) {
        return(list(ok = FALSE, reason = "budget"))
      }
      n1 <- min(hi[1L], max(lo[1L], vb / cost[1L]))
      return(list(ok = TRUE, n = n1, open_takes = numeric(0)))
    } else {
      if (alpha[1L] / hi[1L] > g_target * (1 + 1e-12)) {
        return(list(ok = FALSE, reason = "cv"))
      }
      n1 <- if (alpha[1L] <= 0) lo[1L] else min(hi[1L], max(lo[1L], alpha[1L] / g_target))
      return(list(ok = TRUE, n = n1, open_takes = numeric(0)))
    }
  }

  solve_at_n1 <- function(n1) {
    if (!is.null(vb)) {
      rem_vb <- vb / n1 - cost[1L]
      if (rem_vb < -1e-12 * max(1, vb / n1)) {
        return(list(ok = FALSE, reason = "budget"))
      }
      sub <- .stage_solve_one_rec(
        alpha[-1L], cost[-1L], lo[-1L], hi[-1L],
        vb = max(0, rem_vb)
      )
    } else {
      rem_g <- g_target * n1 - alpha[1L]
      if (rem_g <= 0) {
        at_lo <- c(n1, lo[-1L])
        if (eval_G(at_lo) <= g_target * (1 + 1e-12)) {
          return(list(ok = TRUE, n = at_lo,
                      open_takes = rep(NA_real_, K - 1L)))
        }
        at_hi <- c(n1, hi[-1L])
        if (all(is.finite(hi[-1L])) &&
              eval_G(at_hi) <= g_target * (1 + 1e-12)) {
          return(list(ok = TRUE, n = at_hi,
                      open_takes = rep(NA_real_, K - 1L)))
        }
        return(list(ok = FALSE, reason = "cv"))
      }
      sub <- .stage_solve_one_rec(
        alpha[-1L], cost[-1L], lo[-1L], hi[-1L],
        g_target = rem_g
      )
    }
    if (!sub$ok) return(sub)
    list(ok = TRUE, n = c(n1, sub$n), open_takes = c(NA_real_, sub$open_takes))
  }

  if (abs(hi[1L] - lo[1L]) <= 1e-12 * max(1, lo[1L])) {
    return(solve_at_n1(lo[1L]))
  }

  if (!is.null(vb)) {
    if (eval_H(lo) > vb * (1 + 1e-12)) {
      return(list(ok = FALSE, reason = "budget"))
    }
    if (all(is.finite(hi)) && eval_H(hi) <= vb * (1 + 1e-12)) {
      return(list(ok = TRUE, n = hi, open_takes = rep(NA_real_, K - 1L)))
    }
  } else {
    if (eval_G(hi) > g_target * (1 + 1e-12)) {
      return(list(ok = FALSE, reason = "cv"))
    }
    if (eval_G(lo) <= g_target * (1 + 1e-12)) {
      return(list(ok = TRUE, n = lo, open_takes = rep(NA_real_, K - 1L)))
    }
  }

  open_takes <- rep(NA_real_, K - 1L)
  best_gh <- .stage_min_gh(alpha, cost, lo[-1L], hi[-1L])
  takes <- best_gh$takes
  q <- c(1, cumprod(takes))
  for (idx in seq_len(K - 1L)) {
    k <- idx + 1L
    if (hi[k] > lo[k] * (1 + 1e-14)) {
      A_k <- sum(alpha[seq_len(k - 1L)] / q[seq_len(k - 1L)])
      B_k <- takes[idx] * sum(alpha[k:K] / q[k:K])
      U_k <- sum(cost[seq_len(k - 1L)] * q[seq_len(k - 1L)])
      V_k <- sum(cost[k:K] * q[k:K]) / takes[idx]
      open_takes[idx] <- if (A_k > 0 && B_k > 0) {
        sqrt(B_k * U_k / (A_k * V_k))
      } else {
        Inf
      }
    }
  }

  G_val <- sum(alpha / q)
  H_val <- sum(cost * q)
  n1 <- if (!is.null(vb)) vb / H_val else G_val / g_target

  if (n1 < lo[1L]) {
    res <- solve_at_n1(lo[1L])
    res$open_takes <- open_takes
    return(res)
  }
  if (n1 > hi[1L]) {
    res <- solve_at_n1(hi[1L])
    res$open_takes <- open_takes
    return(res)
  }
  list(ok = TRUE, n = c(n1, takes), open_takes = open_takes)
}

#' Continuous single-indicator S-stage solver wrapper
#' @keywords internal
#' @noRd
.stage_solve_continuous_one <- function(R, rates, M, stage_cost, lower, upper,
                                        cv = NULL, budget = NULL,
                                        fixed_cost = 0) {
  if (!is.matrix(R)) R <- matrix(R, nrow = 1L)
  if (!is.matrix(rates)) rates <- matrix(rates, nrow = 1L)
  if (!is.matrix(M)) M <- matrix(M, nrow = 1L)
  tel <- .stage_telescoping(R, rates, M)
  alpha <- as.numeric(tel$alpha[1L, ])
  kappa <- tel$kappa[1L]
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  g_target <- if (!is.null(cv)) cv^2 + kappa else NULL
  sol <- .stage_solve_one_rec(alpha, stage_cost, lower, upper,
                              g_target = g_target, vb = vb)
  if (!sol$ok) {
    return(sol)
  }
  cv_val <- sqrt(.stage_cv2(sol$n, R, rates, M)[1L])
  cost_val <- .stage_cost(sol$n, stage_cost, fixed_cost)
  list(
    ok = TRUE,
    n = sol$n,
    cv = cv_val,
    cost = cost_val,
    open_takes = sol$open_takes
  )
}

#' Continuous multi-indicator S-stage solver
#'
#' Minimizes cost subject to `cv_j <= cv_target[j]` (CV mode) or minimizes
#' `max_j (cv_j / cv_target[j])^2` subject to `cost <= budget` and
#' `prod(n) >= min_total` (budget mode), for arbitrary S >= 2 and J >= 1.
#' @keywords internal
#' @noRd
.stage_solve_continuous_multi <- function(R, rates, M, cv_target, stage_cost,
                                          lower, upper, budget = NULL,
                                          fixed_cost = 0, min_total = 0) {
  J <- nrow(R)
  S <- ncol(R)
  cv2_target <- cv_target^2
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL

  if (J == 1L && min_total <= 0) {
    return(.stage_solve_continuous_one(
      R, rates, M, stage_cost, lower, upper,
      cv = if (is.null(budget)) cv_target[1L] else NULL,
      budget = budget, fixed_cost = fixed_cost
    ))
  }

  tel <- .stage_telescoping(R, rates, M)
  alpha <- tel$alpha
  kappa <- tel$kappa
  fixed_mask <- abs(upper - lower) <= 1e-12 * pmax(1, lower)
  free_idx <- which(!fixed_mask)
  if (length(free_idx) == 0L) {
    n <- lower
    cvs <- sqrt(.stage_cv2(n, R, rates, M))
    return(list(
      ok = TRUE, n = n, cv_by_target = cvs,
      ratio = max(cvs / cv_target),
      cost = .stage_cost(n, stage_cost, fixed_cost)
    ))
  }

  pivot <- free_idx[1L]
  other_free <- free_idx[-1L]

  solve_pivot <- function(other_vals) {
    n <- lower
    if (length(other_free) > 0L) n[other_free] <- other_vals

    P_pre <- if (pivot == 1L) 1 else prod(n[seq_len(pivot - 1L)])
    C_pre <- if (pivot == 1L) 0 else sum(stage_cost[seq_len(pivot - 1L)] *
                                           cumprod(n[seq_len(pivot - 1L)]))
    q_suf <- c(1, if (pivot < S) cumprod(n[(pivot + 1L):S]) else numeric(0))
    H_suf <- sum(stage_cost[pivot:S] * q_suf) # cost = C_pre + P_pre * n[pivot] * H_suf
    prod_suf <- q_suf[length(q_suf)]          # total_n = P_pre * n[pivot] * prod_suf

    lo_p <- lower[pivot]
    if (min_total > 0) {
      lo_p <- max(lo_p, min_total / (P_pre * prod_suf))
    }
    hi_p <- upper[pivot]

    if (!is.null(vb)) {
      avail <- vb - C_pre
      if (avail < 0) return(NULL)
      n_p <- min(hi_p, avail / (P_pre * H_suf))
      if (n_p < lo_p * (1 - 1e-11)) return(NULL)
      n_p <- max(lo_p, n_p)
      n[pivot] <- n_p
      cv2s <- .stage_cv2(n, R, rates, M)
      return(list(
        n = n,
        obj = max(cv2s / cv2_target),
        cost = C_pre + P_pre * n_p * H_suf,
        cv2 = cv2s
      ))
    } else {
      base_j <- -kappa
      if (pivot > 1L) {
        P_vec <- cumprod(n[seq_len(pivot - 1L)])
        for (s in seq_len(pivot - 1L)) {
          base_j <- base_j + alpha[, s] / P_vec[s]
        }
      }
      G_suf_j <- numeric(J)
      for (k in seq_along(q_suf)) {
        G_suf_j <- G_suf_j + alpha[, pivot + k - 1L] / q_suf[k]
      }
      room_j <- cv2_target - base_j
      n_p <- lo_p
      for (j in seq_len(J)) {
        if (G_suf_j[j] > 0) {
          if (room_j[j] <= 0) return(NULL)
          req <- G_suf_j[j] / (P_pre * room_j[j])
          if (req > n_p) n_p <- req
        } else {
          if (base_j[j] > cv2_target[j] * (1 + 1e-11)) return(NULL)
        }
      }
      if (n_p > hi_p * (1 + 1e-11)) return(NULL)
      n_p <- min(hi_p, max(lo_p, n_p))
      n[pivot] <- n_p
      cv2s <- .stage_cv2(n, R, rates, M)
      if (any(cv2s > cv2_target * (1 + 1e-7) + 1e-15)) return(NULL)
      cost_val <- C_pre + P_pre * n_p * H_suf
      return(list(n = n, obj = cost_val, cost = cost_val, cv2 = cv2s))
    }
  }

  if (length(other_free) == 0L) {
    hit <- solve_pivot(numeric(0))
    if (is.null(hit)) {
      return(list(ok = FALSE, reason = if (!is.null(vb)) "budget" else "cv"))
    }
    cvs <- sqrt(hit$cv2)
    return(list(
      ok = TRUE, n = hit$n, cv_by_target = cvs,
      ratio = max(cvs / cv_target),
      cost = fixed_cost + hit$cost
    ))
  }

  lo_other <- lower[other_free]
  if (!is.null(vb)) {
    if (.stage_cost(lower, stage_cost, 0) > vb * (1 + 1e-12)) {
      return(list(ok = FALSE, reason = "budget"))
    }
    if (min_total > 0 && all(is.finite(upper)) && prod(upper) < min_total * (1 - 1e-12)) {
      return(list(ok = FALSE, reason = "min_total"))
    }
  } else {
    if (any(.stage_cv2(upper, R, rates, M) > cv2_target * (1 + 1e-9) + 1e-15)) {
      return(list(ok = FALSE, reason = "cv"))
    }
  }

  starts <- list()
  for (j in seq_len(J)) {
    s1 <- .stage_solve_continuous_one(
      R[j, , drop = FALSE], rates[j, , drop = FALSE], M[j, , drop = FALSE],
      stage_cost, lower, upper,
      cv = if (is.null(budget)) cv_target[j] else NULL,
      budget = budget, fixed_cost = fixed_cost
    )
    if (s1$ok) {
      starts[[length(starts) + 1L]] <- pmin(upper[other_free],
                                            pmax(lo_other, s1$n[other_free]))
    }
  }
  if (length(starts) == 0L) {
    mid <- ifelse(is.finite(upper[other_free]),
                  sqrt(lo_other * upper[other_free]),
                  pmax(lo_other, 10))
    starts[[1L]] <- mid
  } else if (length(starts) > 1L) {
    mat <- do.call(rbind, starts)
    starts[[length(starts) + 1L]] <- exp(colMeans(log(mat)))
  }
  if (all(is.finite(upper[other_free]))) {
    starts[[length(starts) + 1L]] <- upper[other_free]
  }
  starts[[length(starts) + 1L]] <- lo_other

  ## No design worth having passes these caps (budget or incumbent cost).
  room <- vb
  if (is.null(room)) {
    obj_at <- function(v) {
      hit <- solve_pivot(pmin(upper[other_free], pmax(lo_other, v)))
      if (is.null(hit)) Inf else hit$obj
    }
    room <- min(vapply(starts, obj_at, 1))
    grow <- lo_other
    for (k in seq_len(12L)) {
      if (is.finite(room)) break
      grow <- grow * 10
      room <- obj_at(grow)
    }
    if (!is.finite(room)) return(list(ok = FALSE, reason = "cv"))
  }
  reach <- vapply(other_free, function(s) {
    room / (prod(lower[seq_len(s - 1L)]) * stage_cost[s])
  }, 1)
  upper_cap <- pmin(upper[other_free], pmax(lo_other, reach))
  starts <- lapply(starts, function(v) pmin(upper_cap, v))

  log_lo <- log(lo_other)
  log_hi <- log(upper_cap)
  d <- length(other_free)

  eval_log <- function(v) {
    v <- pmin(log_hi, pmax(log_lo, v))
    hit <- solve_pivot(exp(v))
    if (is.null(hit)) return(1e300)
    hit$obj
  }

  best_hit <- NULL
  update_best <- function(vals) {
    hit <- solve_pivot(pmin(upper_cap, pmax(lo_other, vals)))
    if (!is.null(hit) && (is.null(best_hit) || hit$obj < best_hit$obj)) {
      best_hit <<- hit
    }
    invisible(NULL)
  }

  for (st in starts) update_best(st)

  if (d == 1L) {
    lo_1d <- lo_other[1L]
    hi_1d <- upper_cap[1L]
    if (!is.null(vb) && other_free[1L] > pivot) {
      f_idx <- other_free[1L]
      n_base <- lower
      n_base[f_idx] <- 1
      P_pre <- if (pivot == 1L) 1 else prod(n_base[seq_len(pivot - 1L)])
      C_pre <- if (pivot == 1L) 0 else sum(stage_cost[seq_len(pivot - 1L)] *
                                             cumprod(n_base[seq_len(pivot - 1L)]))
      avail <- vb - C_pre
      q_at1 <- c(1, if (pivot < S) cumprod(n_base[(pivot + 1L):S]) else numeric(0))
      rel_pos <- f_idx - pivot + 1L
      H0 <- sum(stage_cost[pivot:(f_idx - 1L)] * q_at1[seq_len(rel_pos - 1L)])
      H1 <- sum(stage_cost[f_idx:S] * q_at1[rel_pos:length(q_at1)])
      P0 <- q_at1[length(q_at1)]
      if (H1 > 0 && avail > 0) {
        hi_budget <- (avail / (P_pre * lower[pivot]) - H0) / H1
        hi_1d <- min(hi_1d, hi_budget)
      }
      if (min_total > 0) {
        denom <- avail * P0 - min_total * H1
        if (denom > 0) {
          lo_min <- (min_total * H0) / denom
          lo_1d <- max(lo_1d, lo_min)
        } else {
          lo_1d <- Inf
        }
      }
    }
    if (is.finite(lo_1d) && lo_1d <= hi_1d * (1 + 1e-12)) {
      lo_1d <- min(lo_1d, hi_1d)
      update_best(lo_1d)
      update_best(hi_1d)
      l_lo <- log(lo_1d)
      l_hi <- log(hi_1d)
      st_in <- vapply(starts, function(s) log(s[1L]), numeric(1L))
      st_in <- st_in[st_in >= l_lo & st_in <= l_hi]
      grid <- sort(unique(c(seq(l_lo, l_hi, length.out = 65L), st_in)))
      vals <- vapply(grid, eval_log, numeric(1L))
      for (i in which(vals < 1e299)) update_best(exp(grid[i]))
      ord <- order(vals)
      for (idx in ord[seq_len(min(3L, length(ord)))]) {
        if (vals[idx] >= 1e299) break
        a_int <- if (idx > 1L) grid[idx - 1L] else grid[idx]
        b_int <- if (idx < length(grid)) grid[idx + 1L] else grid[idx]
        if (b_int > a_int) {
          opt <- stats::optimize(eval_log, interval = c(a_int, b_int), tol = 1e-12)
          update_best(exp(opt$minimum))
        }
      }
      if (l_hi > l_lo) {
        opt_all <- stats::optimize(eval_log, interval = c(l_lo, l_hi), tol = 1e-12)
        update_best(exp(opt_all$minimum))
      }
    }
  } else {
    smooth_eval <- function(v, tau) {
      v <- pmin(log_hi, pmax(log_lo, v))
      vals_o <- exp(v)
      n <- lower
      n[other_free] <- vals_o
      P_pre <- if (pivot == 1L) 1 else prod(n[seq_len(pivot - 1L)])
      C_pre <- if (pivot == 1L) 0 else sum(stage_cost[seq_len(pivot - 1L)] *
                                             cumprod(n[seq_len(pivot - 1L)]))
      q_suf <- c(1, if (pivot < S) cumprod(n[(pivot + 1L):S]) else numeric(0))
      H_suf <- sum(stage_cost[pivot:S] * q_suf)
      prod_suf <- q_suf[length(q_suf)]
      lo_p <- lower[pivot]
      if (min_total > 0) lo_p <- max(lo_p, min_total / (P_pre * prod_suf))
      hi_p <- upper[pivot]

      if (!is.null(vb)) {
        avail <- vb - C_pre
        if (avail <= 0) return(1e12 + sum((v - log_lo)^2))
        n_p <- min(hi_p, avail / (P_pre * H_suf))
        if (n_p < lo_p) {
          return(1e8 + 1e6 * (lo_p / max(n_p, 1e-12) - 1)^2)
        }
        n[pivot] <- n_p
        rats <- .stage_cv2(n, R, rates, M) / cv2_target
        mx <- max(rats)
        if (mx <= 0) return(0)
        mx * (sum((rats / mx)^tau))^(1 / tau)
      } else {
        base_j <- -kappa
        if (pivot > 1L) {
          P_vec <- cumprod(n[seq_len(pivot - 1L)])
          for (s in seq_len(pivot - 1L)) base_j <- base_j + alpha[, s] / P_vec[s]
        }
        G_suf_j <- numeric(J)
        for (k in seq_along(q_suf)) {
          G_suf_j <- G_suf_j + alpha[, pivot + k - 1L] / q_suf[k]
        }
        room_j <- cv2_target - base_j
        if (any(G_suf_j > 0 & room_j <= 0)) {
          return(1e12 + sum(pmax(0, -room_j) * 1e6))
        }
        reqs <- ifelse(G_suf_j > 0, G_suf_j / (P_pre * room_j), lo_p)
        reqs <- pmax(reqs, lo_p)
        mx <- max(reqs)
        n_p_sm <- mx * (sum((reqs / mx)^tau))^(1 / tau)
        pen <- if (n_p_sm > hi_p) 1e8 * (n_p_sm / hi_p - 1)^2 else 0
        C_pre + P_pre * min(hi_p, n_p_sm) * H_suf + pen
      }
    }

    seed_logs <- list()
    for (st in starts) {
      sv <- pmin(log_hi, pmax(log_lo, log(st)))
      if (length(seed_logs) == 0L ||
          all(vapply(seed_logs, function(z) max(abs(z - sv)) > 1e-3, logical(1L)))) {
        seed_logs[[length(seed_logs) + 1L]] <- sv
      }
    }
    if (length(seed_logs) > 5L) seed_logs <- seed_logs[seq_len(5L)]

    for (v0 in seed_logs) {
      v <- v0
      for (tau in c(20, 200)) {
        opt <- tryCatch(
          stats::optim(
            v, smooth_eval, tau = tau, method = "L-BFGS-B",
            lower = log_lo, upper = log_hi,
            control = list(maxit = 120L, factr = 1e6)
          ),
          error = function(e) NULL
        )
        if (!is.null(opt) && all(is.finite(opt$par))) {
          v <- opt$par
          update_best(exp(v))
        }
      }
    }

    if (!is.null(best_hit)) {
      v_best <- pmin(log_hi, pmax(log_lo, log(best_hit$n[other_free])))
      nm <- tryCatch(
        stats::optim(
          v_best,
          function(v) {
            if (any(v < log_lo - 1e-12) || any(v > log_hi + 1e-12)) return(1e300)
            eval_log(v)
          },
          method = "Nelder-Mead",
          control = list(maxit = 400L, reltol = 1e-11)
        ),
        error = function(e) NULL
      )
      if (!is.null(nm) && all(is.finite(nm$par))) {
        update_best(exp(nm$par))
      }

      v_cur <- pmin(log_hi, pmax(log_lo, log(best_hit$n[other_free])))
      for (pass in seq_len(2L)) {
        for (k in seq_len(d)) {
          if (log_hi[k] <= log_lo[k]) next
          fn1d <- function(x) {
            vv <- v_cur
            vv[k] <- x
            eval_log(vv)
          }
          lo_k <- max(log_lo[k], v_cur[k] - 2)
          hi_k <- min(log_hi[k], v_cur[k] + 2)
          opt1 <- stats::optimize(fn1d, interval = c(lo_k, hi_k), tol = 1e-11)
          if (opt1$objective < best_hit$obj) {
            v_cur[k] <- opt1$minimum
            update_best(exp(v_cur))
          }
        }
      }
    }
  }

  if (is.null(best_hit)) {
    return(list(ok = FALSE, reason = if (!is.null(vb)) "budget" else "cv"))
  }

  if (!is.null(vb) && is.finite(upper[pivot]) &&
      best_hit$n[pivot] >= upper[pivot] * (1 - 1e-10) &&
      best_hit$cost < vb * (1 - 1e-9) && length(other_free) >= 1L) {
    lo2 <- lower
    hi2 <- upper
    lo2[pivot] <- upper[pivot]
    hi2[pivot] <- upper[pivot]
    sub <- .stage_solve_continuous_multi(
      R, rates, M, cv_target, stage_cost, lo2, hi2,
      budget = budget, fixed_cost = fixed_cost, min_total = min_total
    )
    if (sub$ok && sub$ratio < sqrt(best_hit$obj)) {
      return(sub)
    }
  }

  cvs <- sqrt(best_hit$cv2)
  list(
    ok = TRUE,
    n = best_hit$n,
    cv_by_target = cvs,
    ratio = max(cvs / cv_target),
    cost = fixed_cost + best_hit$cost
  )
}

#' The stage search's node limit from `control`
#' @keywords internal
#' @noRd
.stage_max_nodes <- function(control) {
  if (is.list(control) && "gap" %in% names(control)) {
    stop(
      "'control$gap' does not apply to n_cluster(): its search keeps no lower bound to measure a gap against",
      call. = FALSE
    )
  }
  checked <- .search_control(control)
  if (is.null(control[["max_nodes"]])) 2e7 else checked$max_nodes
}

#' What the stage search proved, kept on the result
#' @keywords internal
#' @noRd
.stage_search_record <- function(disc, max_nodes) {
  list(
    status = if (isTRUE(disc$exact)) "optimal" else "limit",
    nodes = disc$nodes,
    max_nodes = max_nodes
  )
}

#' The domain searches as one record, `"limit"` when any stopped early
#' @keywords internal
#' @noRd
.stage_search_domains <- function(searches, domain_levels) {
  status <- vapply(searches, function(s) s$status, character(1))
  names(status) <- domain_levels
  list(
    status = if (any(status == "limit")) "limit" else "optimal",
    nodes = sum(vapply(searches, function(s) s$nodes, numeric(1))),
    max_nodes = searches[[1L]]$max_nodes,
    by_domain = status,
    nodes_by_domain = stats::setNames(
      vapply(searches, function(s) s$nodes, numeric(1)), domain_levels
    )
  )
}

#' Whole-unit S-stage discrete solver via C branch-and-bound
#' @keywords internal
#' @noRd
.stage_solve_discrete <- function(R, rates, M, cv2_target, stage_cost,
                                  lower, upper, budget_var = NA_real_,
                                  min_total = 0, n_init = NULL,
                                  max_nodes = 2e7) {
  if (!is.matrix(R)) R <- matrix(as.numeric(R), nrow = 1L)
  if (!is.matrix(rates)) rates <- matrix(as.numeric(rates), nrow = 1L)
  if (!is.matrix(M)) M <- matrix(as.numeric(M), nrow = 1L)
  storage.mode(R) <- "double"
  storage.mode(rates) <- "double"
  storage.mode(M) <- "double"
  cap <- .Machine$integer.max
  upper <- pmin(as.numeric(upper), cap)
  out <- .Call(
    C_stage_search,
    R,
    rates,
    M,
    as.double(cv2_target),
    as.double(stage_cost),
    as.double(lower),
    as.double(upper),
    as.double(budget_var),
    as.double(min_total),
    if (is.null(n_init)) NULL else as.double(n_init),
    as.double(max_nodes)
  )
  big_init <- !is.null(n_init) && any(n_init > cap, na.rm = TRUE)
  if ((isTRUE(out$found) && any(out$n >= cap)) ||
        (!isTRUE(out$found) && big_init)) {
    s <- if (isTRUE(out$found)) which(out$n >= cap)[1L] else
      which(n_init > cap)[1L]
    .stop_count_range("the whole-unit search", sprintf("at stage %d", s))
  }
  if (!isTRUE(out$found) && !isTRUE(out$exact)) {
    stop(structure(
      list(message = sprintf(
        "no whole-unit design was found within %s search nodes, so whether one exists is unknown. Raise 'control$max_nodes'",
        format(out$nodes, big.mark = ",", scientific = FALSE)
      ), call = NULL, nodes = out$nodes),
      class = c("svyplan_plan_unknown", "error", "condition")
    ))
  }
  out
}
