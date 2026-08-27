#' @keywords internal
#' @noRd
.alloc_weights <- function(alloc, q, N_h, S_h, cost_h,
                           deff_h = 1, resp_rate_h = 1) {
  # Proportional is a count rule rather than an optimum, so it takes 1 / R
  # and no deff. Both collapse to the plain weights under scalar inputs.
  var_adj <- sqrt(deff_h / resp_rate_h)
  switch(
    alloc,
    proportional = N_h / resp_rate_h,
    neyman = N_h * S_h * var_adj,
    optimal = N_h * S_h * var_adj / sqrt(cost_h),
    power = S_h * N_h^q * var_adj
  )
}

#' @keywords internal
#' @noRd
.strata_precompute <- function(x_sort) {
  list(
    cs = c(0, cumsum(x_sort)),
    cs2 = c(0, cumsum(x_sort^2)),
    n = length(x_sort),
    x_sort = x_sort
  )
}

#' @keywords internal
#' @noRd
.bk_to_idx <- function(x_sort, bk) {
  c(0L, findInterval(bk, x_sort), length(x_sort))
}

#' @keywords internal
#' @noRd
.strata_stats_from_prefix <- function(pre, idx) {
  L <- length(idx) - 1L
  N_h <- diff(idx)
  lo <- idx[seq_len(L)] + 1L
  hi <- idx[seq_len(L) + 1L] + 1L
  sum_h <- pre$cs[hi] - pre$cs[lo]
  sum2_h <- pre$cs2[hi] - pre$cs2[lo]
  mean_h <- ifelse(N_h > 0L, sum_h / N_h, 0)
  var_h <- ifelse(N_h > 1L, pmax(0, (sum2_h - N_h * mean_h^2) / (N_h - 1L)), 0)
  list(
    N_h = N_h,
    W_h = N_h / pre$n,
    S_h = sqrt(var_h),
    mean_h = mean_h
  )
}

#' @keywords internal
#' @noRd
.rna_alloc <- function(a_h, n_total, m_h, M_h) {
  if (!is.numeric(a_h) || !is.numeric(m_h) || !is.numeric(M_h) ||
      length(a_h) == 0L || length(m_h) != length(a_h) ||
      length(M_h) != length(a_h) || anyNA(a_h) || anyNA(m_h) ||
      anyNA(M_h) || any(!is.finite(a_h)) || any(!is.finite(m_h)) ||
      any(!is.finite(M_h)) || any(a_h < 0) || any(m_h > M_h)) {
    stop("internal error: invalid bounded-allocation inputs", call. = FALSE)
  }

  if (!is.numeric(n_total) || length(n_total) != 1L || is.na(n_total) ||
      !is.finite(n_total)) {
    stop("internal error: invalid bounded-allocation total", call. = FALSE)
  }
  tol <- max(1e-10, 1e-10 * max(1, abs(n_total)))
  lo_total <- sum(m_h)
  hi_total <- sum(M_h)
  if (n_total < lo_total - tol ||
      n_total > hi_total + tol) {
    stop(structure(
      list(message = sprintf(
        "bounded allocation total %.10g is infeasible; feasible range is %.10g to %.10g",
        n_total, lo_total, hi_total
      )),
      class = c("svyplan_alloc_infeasible", "error", "condition")
    ))
  }

  n_total <- min(max(n_total, lo_total), hi_total)
  if (abs(n_total - lo_total) <= tol) {
    return(.rna_check(as.numeric(m_h), n_total, m_h, M_h, tol))
  }
  if (abs(n_total - hi_total) <= tol) {
    return(.rna_check(as.numeric(M_h), n_total, m_h, M_h, tol))
  }

  # The plain ratio allocation, when it already lies inside every bound, is
  # the bounded solution: no breakpoint can bind, so the search below would
  # return this same vector. Skipping it is exact rather than approximate.
  # A zero weight sits at 0 here, which satisfies its own lower bound only
  # when that bound is 0, so a fixed or take-all component declines this
  # branch on its own and takes the general path. Placed ahead of the
  # positive-weight bookkeeping so an unconstrained call never pays for it.
  sa <- sum(a_h)
  if (sa > 0) {
    raw <- n_total * a_h / sa
    if (all(raw >= m_h) && all(raw <= M_h) &&
        all(is.finite(raw)) && abs(sum(raw) - n_total) <= tol) {
      return(.rna_check(raw, n_total, m_h, M_h, tol))
    }
  }

  # Positive weights follow the bounded ratio rule n_h = clamp(lambda a_h).
  # Zero-weight strata remain at their lower bounds unless all positive-weight
  # capacity is exhausted, after which they share the unavoidable residual.
  pos <- a_h > 0
  reachable <- sum(ifelse(pos, M_h, m_h))
  target_pos <- min(n_total, reachable)
  n_h <- as.numeric(m_h)

  if (any(pos) && target_pos > lo_total + tol) {
    idx_pos <- which(pos)
    budget <- target_pos - sum(m_h[!pos])
    lower_bp <- m_h[idx_pos] / a_h[idx_pos]
    upper_bp <- M_h[idx_pos] / a_h[idx_pos]
    breakpoints <- sort(unique(c(lower_bp, upper_bp)))
    at_lambda <- function(lambda) {
      sum(pmin(pmax(lambda * a_h[idx_pos], m_h[idx_pos]), M_h[idx_pos]))
    }
    values <- vapply(breakpoints, at_lambda, numeric(1L))
    hit <- which(values >= budget - tol)[1L]

    if (abs(values[hit] - budget) <= tol) {
      lambda <- breakpoints[hit]
    } else {
      lambda_lo <- if (hit == 1L) 0 else breakpoints[hit - 1L]
      lambda_hi <- breakpoints[hit]
      lambda_mid <- (lambda_lo + lambda_hi) / 2
      raw_mid <- lambda_mid * a_h[idx_pos]
      active <- raw_mid > m_h[idx_pos] & raw_mid < M_h[idx_pos]
      fixed <- ifelse(raw_mid <= m_h[idx_pos], m_h[idx_pos], M_h[idx_pos])
      fixed[active] <- 0
      lambda <- (budget - sum(fixed)) / sum(a_h[idx_pos][active])
    }
    n_h[idx_pos] <- pmin(
      pmax(lambda * a_h[idx_pos], m_h[idx_pos]),
      M_h[idx_pos]
    )
  }

  residual <- n_total - sum(n_h)
  if (residual > tol) {
    zero <- which(!pos & M_h > n_h + tol)
    capacity <- M_h[zero] - n_h[zero]
    caps <- sort(unique(capacity))
    used <- 0
    remaining <- length(capacity)
    level <- 0
    for (cap in caps) {
      if (used + remaining * (cap - level) >= residual - tol) {
        level <- level + (residual - used) / remaining
        break
      }
      used <- used + remaining * (cap - level)
      remaining <- remaining - sum(capacity == cap)
      level <- cap
    }
    n_h[zero] <- n_h[zero] + pmin(level, capacity)
  }

  total_gap <- n_total - sum(n_h)
  if (abs(total_gap) <= tol && total_gap != 0) {
    candidates <- which(
      n_h + total_gap >= m_h - tol & n_h + total_gap <= M_h + tol
    )
    if (length(candidates) > 0L) {
      idx <- candidates[1L]
      n_h[idx] <- n_h[idx] + total_gap
    }
  }
  .rna_check(n_h, n_total, m_h, M_h, tol)
}

#' Post-condition every bounded allocation must satisfy
#'
#' Both the fast path and the breakpoint search return through here, so no
#' route out of `.rna_alloc()` can skip the total-and-bounds guarantee its
#' callers rely on.
#' @keywords internal
#' @noRd
.rna_check <- function(n_h, n_total, m_h, M_h, tol) {
  if (abs(sum(n_h) - n_total) > tol || any(n_h < m_h - tol) ||
      any(n_h > M_h + tol)) {
    stop("internal error: bounded allocation violated its total or bounds",
         call. = FALSE)
  }
  n_h
}

#' Reject a stratification that cannot support its own allocation
#'
#' A sampled stratum needs two units. `m_h = pmin(2, N_h)` falls to one for a
#' single-unit stratum, so shrinking a stratum to one unit lowers the minimum
#' feasible total and lets a search satisfy a total that should have been
#' refused. Such a stratum also carries `sd = 0` and supports no within-stratum
#' variance estimate. A take-all stratum is enumerated rather than sampled, so
#' one unit is enough there.
#' @keywords internal
#' @noRd
.strata_degenerate <- function(N_h, take_all_idx = NULL) {
  floor_h <- rep(2L, length(N_h))
  if (!is.null(take_all_idx)) {
    floor_h[take_all_idx] <- 1L
  }
  any(N_h < floor_h)
}

#' Stratified variance of the mean under the package's shared convention
#'
#' `V = deff * sum(W_h^2 S_h^2 (1 / n_net_h - 1 / N_h))` with
#' `n_net_h = n_h * resp_rate`, the same expression `.alloc_metrics()` uses,
#' so that `strata_bound()` and `n_alloc()` report the same `cv` for the same
#' design. A scalar `deff` scales every candidate boundary set equally and so
#' leaves the optimal boundaries unchanged. What it changes is the `n` a `cv`
#' target requires and the `cv` a given `n` achieves.
#' @keywords internal
#' @noRd
.strata_variance <- function(W_h, S_h, n_h, N_h, deff = 1, resp_rate = 1) {
  n_net <- n_h * resp_rate
  n_eff <- n_net / deff
  active <- N_h > 0 & n_eff > 0
  fpc <- pmax(0, 1 - n_net[active] / N_h[active])
  sum(W_h[active]^2 * S_h[active]^2 * fpc / n_eff[active])
}

#' Smallest total whose allocation reaches a target variance
#'
#' The one sizing routine behind both the LH/general path and the Kozak inner
#' objective, so the two cannot drift apart on formula, tolerance, or edge
#' behaviour.
#'
#' While no bound binds, the allocation is proportional, \eqn{n_h = n c_h}
#' with \eqn{c_h = a_h / \sum a_h}, and the package's variance is
#' \eqn{V(n) = A/n - B} with
#' \eqn{A = (deff / resp\_rate) \sum W_h^2 S_h^2 / c_h} and
#' \eqn{B = deff \sum W_h^2 S_h^2 / N_h}. That inverts exactly to
#' \eqn{n = A / (V_{target} + B)}, so the ordinary case needs no search at
#' all. The candidate is accepted only when it lands inside the feasible
#' total, keeps every \eqn{n c_h} within its bounds, and is confirmed to meet
#' the target once. A fixed or take-all component breaks the proportional
#' form, so those cases take the bracket below.
#'
#' The fallback keeps a directional bracket: `lo` is known to miss the target
#' and `hi` to meet it, so returning `hi` preserves the contract that the
#' continuous size achieves the requested precision. It stops on interval
#' width rather than a fixed iteration count, with a cap as a backstop
#' against a non-monotone objective.
#' @keywords internal
#' @noRd
.strata_size_for_target <- function(W_h, S_h, N_h, a_h, m_h, M_h, target_V,
                                    deff = 1, resp_rate = 1,
                                    abs_tol = 1e-8, rel_tol = 1e-12,
                                    max_iter = 200L) {
  lo <- sum(m_h)
  hi <- sum(M_h)

  variance_at <- function(n_total) {
    n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
    .strata_variance(W_h, S_h, n_h, N_h, deff, resp_rate)
  }

  if (variance_at(lo) <= target_V) {
    return(lo)
  }
  if (variance_at(hi) > target_V) {
    return(Inf)
  }

  if (all(a_h > 0)) {
    c_h <- a_h / sum(a_h)
    ws <- W_h^2 * S_h^2
    A <- deff / resp_rate * sum(ws / c_h)
    B <- deff * sum(ws / N_h)
    if (is.finite(A) && A > 0 && is.finite(B)) {
      cand <- A / (target_V + B)
      if (is.finite(cand) && cand >= lo && cand <= hi) {
        n_cand <- cand * c_h
        if (all(n_cand >= m_h) && all(n_cand <= M_h) &&
            variance_at(cand) <= target_V) {
          return(cand)
        }
      }
    }
  }

  for (i in seq_len(max_iter)) {
    if (hi - lo <= max(abs_tol, rel_tol * max(1, hi))) {
      break
    }
    mid <- (lo + hi) / 2
    if (variance_at(mid) > target_V) lo <- mid else hi <- mid
  }
  hi
}

#' Evaluate stratified allocation for given boundaries
#' @keywords internal
#' @noRd
.strata_alloc <- function(
  x,
  bk,
  n_total,
  alloc,
  q,
  cost_h,
  take_all_idx = NULL,
  .pre = NULL,
  deff = 1,
  resp_rate = 1
) {
  L <- length(bk) + 1L

  if (!is.null(.pre)) {
    idx <- .bk_to_idx(.pre$x_sort, bk)
    stats <- .strata_stats_from_prefix(.pre, idx)
    N_h <- stats$N_h
    N <- .pre$n
    W_h <- stats$W_h
    S_h <- stats$S_h
    mean_h <- stats$mean_h
    breaks <- c(.pre$x_sort[1L], bk, .pre$x_sort[N])
  } else {
    x_range <- range(x)
    breaks <- c(x_range[1L], bk, x_range[2L])
    bins <- findInterval(x, bk, left.open = TRUE) + 1L
    if (!is.null(take_all_idx)) {
      # Ordinary cutpoints define right-closed strata, but the documented
      # take-all rule is x >= the threshold.  Move equality at that threshold into
      # the take-all stratum explicitly.
      bins[x >= bk[take_all_idx - 1L]] <- take_all_idx
    }
    N_h <- tabulate(bins, nbins = L)
    N <- length(x)
    W_h <- N_h / N

    bins_f <- factor(bins, levels = seq_len(L))
    x_split <- split(x, bins_f)
    S_h <- vapply(
      x_split,
      function(xi) {
        if (length(xi) < 2L) {
          return(0)
        }
        sqrt(var(xi))
      },
      numeric(1L)
    )
    mean_h <- vapply(x_split, mean, numeric(1L))
    mean_h[is.nan(mean_h)] <- 0
  }

  a_h <- .alloc_weights(alloc, q, N_h, S_h, cost_h)
  m_h <- pmin(rep(2, L), N_h)
  M_h <- N_h
  if (!is.null(take_all_idx)) {
    a_h[take_all_idx] <- 0
    m_h[take_all_idx] <- N_h[take_all_idx]
  }

  n_h <- .rna_alloc(a_h, n_total, m_h, M_h)

  V <- .strata_variance(W_h, S_h, n_h, N_h, deff, resp_rate)
  ybar <- .aggregate_mean(W_h, mean_h)
  cv <- if (ybar == 0) Inf else sqrt(V) / abs(ybar)

  list(
    N_h = N_h,
    W_h = W_h,
    S_h = S_h,
    mean_h = mean_h,
    n_h = n_h,
    cv = cv,
    V = V,
    lower = breaks[-length(breaks)],
    upper = breaks[-1L]
  )
}

#' Objective function: CV given boundaries (for minimization)
#' @keywords internal
#' @noRd
.strata_obj <- function(
  x,
  bk,
  n_total,
  alloc,
  q,
  cost_h,
  take_all_idx = NULL,
  .pre = NULL,
  deff = 1,
  resp_rate = 1
) {
  if (is.unsorted(bk)) {
    return(Inf)
  }
  res <- tryCatch(
    .strata_alloc(x, bk, n_total, alloc, q, cost_h, take_all_idx, .pre,
                  deff, resp_rate),
    svyplan_alloc_infeasible = function(e) NULL
  )
  if (is.null(res)) {
    return(Inf)
  }
  if (.strata_degenerate(res$N_h, take_all_idx)) {
    return(Inf)
  }
  res$cv
}

#' Required n to achieve target CV
#' @keywords internal
#' @noRd
.strata_n_for_cv <- function(
  x,
  bk,
  target_cv,
  alloc,
  q,
  cost_h,
  take_all_idx = NULL,
  .pre = NULL,
  deff = 1,
  resp_rate = 1,
  target_V = NULL
) {
  if (is.unsorted(bk)) {
    return(Inf)
  }
  L <- length(bk) + 1L

  if (!is.null(.pre)) {
    idx <- .bk_to_idx(.pre$x_sort, bk)
    stats <- .strata_stats_from_prefix(.pre, idx)
    N_h <- stats$N_h
    if (.strata_degenerate(N_h, take_all_idx)) {
      return(Inf)
    }
    N <- .pre$n
    W_h <- stats$W_h
    S_h <- stats$S_h
    mean_h <- stats$mean_h
  } else {
    bins <- findInterval(x, bk, left.open = TRUE) + 1L
    if (!is.null(take_all_idx)) {
      bins[x >= bk[take_all_idx - 1L]] <- take_all_idx
    }
    N_h <- tabulate(bins, nbins = L)
    if (.strata_degenerate(N_h, take_all_idx)) {
      return(Inf)
    }
    N <- length(x)
    W_h <- N_h / N

    bins_f <- factor(bins, levels = seq_len(L))
    x_split <- split(x, bins_f)
    S_h <- vapply(
      x_split,
      function(xi) {
        if (length(xi) < 2L) {
          return(0)
        }
        sqrt(var(xi))
      },
      numeric(1L)
    )
    mean_h <- vapply(x_split, mean, numeric(1L))
    mean_h[is.nan(mean_h)] <- 0
  }

  ybar <- .aggregate_mean(W_h, mean_h)
  if (is.null(target_V)) {
    target_V <- (target_cv * abs(ybar))^2
  }

  a_h <- .alloc_weights(alloc, q, N_h, S_h, cost_h)
  m_h <- pmin(rep(2, L), N_h)
  M_h <- N_h
  if (!is.null(take_all_idx)) {
    a_h[take_all_idx] <- 0
    m_h[take_all_idx] <- N_h[take_all_idx]
  }

  .strata_size_for_target(W_h, S_h, N_h, a_h, m_h, M_h, target_V,
                          deff, resp_rate)
}

#' Dalenius-Hodges cumulative sqrt(f) rule
#' @keywords internal
#' @noRd
.strata_cumrootf <- function(x, L, n_class = NULL, warn = TRUE) {
  if (is.null(n_class)) {
    n_class <- nclass.FD(x)
  }
  n_class <- max(n_class, L)
  h <- hist(x, breaks = n_class, plot = FALSE)
  csf <- cumsum(sqrt(h$counts))
  total <- csf[length(csf)]
  targets <- total * seq_len(L - 1L) / L
  idx <- vapply(targets, function(t) which.min(abs(csf - t)), integer(1L))
  idx <- pmin(idx, length(h$breaks) - 1L)
  bk <- unique(h$breaks[idx + 1L])
  if (length(bk) < L - 1L) {
    bk <- .strata_cumrootf_discrete(x, L, warn = warn)
  }
  bk
}

#' Fallback boundaries for concentrated or discrete distributions
#'
#' Applies the cumulative-root-frequency rule to the distinct observed
#' values and forces strictly increasing cuts, so any input with at
#' least `L` distinct values yields `L` nonempty strata. Boundaries are
#' placed midway between adjacent distinct values.
#' @keywords internal
#' @noRd
.strata_cumrootf_discrete <- function(x, L, warn = TRUE) {
  ux <- sort(unique(x))
  nu <- length(ux)
  if (nu < L) {
    stop(
      sprintf(
        "'cumrootf' cannot form %d strata from %d distinct values; reduce 'n_strata'",
        L, nu
      ),
      call. = FALSE
    )
  }
  if (warn) {
    warning(
      "cumulative-root-frequency boundaries are degenerate for this distribution; boundaries placed between adjacent distinct values",
      call. = FALSE
    )
  }
  counts <- tabulate(match(x, ux), nbins = nu)
  csf <- cumsum(sqrt(counts))
  targets <- csf[nu] * seq_len(L - 1L) / L
  idx <- vapply(targets, function(t) which.min(abs(csf - t)), integer(1L))
  for (j in seq_along(idx)) {
    lo <- if (j == 1L) 1L else idx[j - 1L] + 1L
    hi <- nu - (L - 1L) + j - 1L
    idx[j] <- min(max(idx[j], lo), hi)
  }
  (ux[idx] + ux[idx + 1L]) / 2
}

#' Geometric progression boundaries
#' @keywords internal
#' @noRd
.strata_geo <- function(x, L) {
  b0 <- min(x)
  bL <- max(x)
  r <- (bL / b0)^(1 / L)
  b0 * r^seq_len(L - 1L)
}

#' LH-inspired coordinate-wise boundary optimization
#' @keywords internal
#' @noRd
.strata_lh <- function(
  x_sort,
  L,
  n_total,
  target_cv,
  alloc,
  q,
  cost_h,
  max_iter,
  take_all_idx = NULL,
  deff = 1,
  resp_rate = 1,
  target_V = NULL
) {
  x_uniq <- sort(unique(x_sort))
  nu <- length(x_uniq)

  if (nu < L) {
    stop("fewer unique values than requested strata", call. = FALSE)
  }

  quant_p <- seq(0, 1, length.out = L + 1L)[2:L]
  bk <- unname(quantile(x_sort, probs = quant_p))
  bk <- pmin(bk, x_uniq[nu - 1L])
  bk <- pmax(bk, x_uniq[2L])

  pre <- .strata_precompute(x_sort)

  use_cv <- !is.null(target_cv)
  raw_obj <- if (use_cv) {
    function(bk_) {
      .strata_n_for_cv(
        x_sort,
        bk_,
        target_cv,
        alloc,
        q,
        cost_h,
        take_all_idx,
        .pre = pre,
        deff = deff,
        resp_rate = resp_rate,
        target_V = target_V
      )
    }
  } else {
    function(bk_) {
      .strata_obj(
        x_sort,
        bk_,
        n_total,
        alloc,
        q,
        cost_h,
        take_all_idx,
        .pre = pre,
        deff = deff,
        resp_rate = resp_rate
      )
    }
  }

  # The objective reads the population partition, so any two boundary vectors
  # falling between the same pair of observations describe the same design and
  # share an answer. `optimize()` probes a whole interval per coordinate and
  # revisits the same partitions across sweeps, so a fit-local cache keyed by
  # the cut indices removes most of the work. It is fit-local because the
  # objective also depends on the data, allocation, target, costs, deff and
  # response rate, none of which vary within one call. Infinite objectives are
  # cached too: a degenerate partition is just as repeatable as a good one.
  memo <- new.env(hash = TRUE, parent = emptyenv())
  obj_fn <- function(bk_) {
    key <- paste0(findInterval(bk_, x_sort), collapse = ",")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    val <- raw_obj(bk_)
    assign(key, val, envir = memo)
    val
  }

  # A boundary anywhere between two adjacent observations describes the same
  # design, so the search is really over cut indices and only their motion is
  # progress. Pinning each boundary to the observation it sits above keeps the
  # next sweep's coordinate intervals stable, which is what makes an unchanged
  # index vector mean the search has actually stopped rather than merely
  # drifted inside one partition. Boundaries that would not round-trip, at a
  # run of tied values, are left where they are.
  canonicalize <- function(bk_) {
    p <- findInterval(bk_, x_sort)
    keep <- p >= 1L
    if (!any(keep)) {
      return(bk_)
    }
    cand <- x_sort[pmax(p, 1L)]
    ok <- keep & findInterval(cand, x_sort) == p
    bk_[ok] <- cand[ok]
    bk_
  }
  key_of <- function(bk_) paste0(findInterval(bk_, x_sort), collapse = ",")

  # Canonicalize before recording anything. A sweep from non-canonical
  # boundaries searches different coordinate intervals than a sweep from the
  # canonical representative of the same partition, so comparing the two
  # would compare sweeps from different starting points and could call a
  # first sweep converged when another sweep from the boundaries actually
  # returned still moves.
  bk <- canonicalize(bk)
  best_obj <- obj_fn(bk)
  if (!is.finite(best_obj)) {
    best_obj <- Inf
  }
  best_bk <- bk
  best_key <- key_of(bk)
  converged <- FALSE
  tol <- diff(range(x_sort)) * 1e-6

  seen <- new.env(hash = TRUE, parent = emptyenv())
  prev_key <- best_key
  assign(prev_key, TRUE, envir = seen)

  for (iter in seq_len(max_iter)) {
    for (h in seq_len(L - 1L)) {
      lo <- if (h == 1L) x_uniq[1L] + tol else bk[h - 1L] + tol
      hi <- if (h == L - 1L) x_uniq[nu] - tol else bk[h + 1L] - tol
      if (lo >= hi) {
        next
      }
      opt <- suppressWarnings(optimize(
        function(b) {
          bk_try <- bk
          bk_try[h] <- b
          obj_fn(bk_try)
        },
        interval = c(lo, hi)
      ))
      # `optimize()` assumes a continuous unimodal objective and this one is a
      # step function, so its proposal can be worse than where the coordinate
      # already sits. Take it only when it does not worsen the objective, or
      # when the current state is infeasible and any move is an escape.
      cand <- bk
      cand[h] <- opt$minimum
      cur_obj <- obj_fn(bk)
      if (!is.finite(cur_obj) || obj_fn(cand) <= cur_obj) {
        bk <- cand
      }
    }
    bk <- canonicalize(bk)
    new_obj <- obj_fn(bk)
    key <- key_of(bk)
    if (is.finite(new_obj) && (is.infinite(best_obj) || new_obj < best_obj)) {
      best_obj <- new_obj
      best_bk <- bk
      best_key <- key
    }

    if (identical(key, prev_key)) {
      # Convergence has to describe the design actually returned, which is the
      # best partition seen rather than necessarily the last one visited.
      converged <- is.finite(new_obj) && identical(key, best_key)
      break
    }
    if (!is.null(seen[[key]])) {
      # The sweep returned to a partition it has already left. Further sweeps
      # retrace the same cycle, so stop with the best partition seen and say
      # plainly that this is not a fixed point.
      break
    }
    assign(key, TRUE, envir = seen)
    prev_key <- key
  }
  list(
    bk = canonicalize(best_bk),
    converged = converged,
    feasible = is.finite(best_obj)
  )
}

#' Kozak-inspired random-restart adjacent-boundary local search
#' @keywords internal
#' @noRd
.strata_kozak <- function(
  x_sort,
  L,
  n_total,
  target_cv,
  alloc,
  q,
  cost_h,
  max_iter,
  n_restart,
  take_all_idx = NULL,
  deff = 1,
  resp_rate = 1,
  target_V = NULL
) {
  x_uniq <- sort(unique(x_sort))
  nu <- length(x_uniq)

  if (nu < L) {
    stop("fewer unique values than requested strata", call. = FALSE)
  }

  pre <- .strata_precompute(x_sort)
  N <- pre$n
  cs <- pre$cs
  cs2 <- pre$cs2
  use_cv <- !is.null(target_cv)

  # Map each unique value to its rightmost position in x_sort (prefix index)
  uniq_pidx <- findInterval(x_uniq, x_sort)

  obj_from_pidx <- function(pidx) {
    lo <- pidx[seq_len(L)] + 1L
    hi <- pidx[seq_len(L) + 1L] + 1L
    N_h <- diff(pidx)
    if (.strata_degenerate(N_h, take_all_idx)) {
      return(Inf)
    }
    sum_h <- cs[hi] - cs[lo]
    sum2_h <- cs2[hi] - cs2[lo]
    mean_h <- ifelse(N_h > 0L, sum_h / N_h, 0)
    var_h <- ifelse(
      N_h > 1L,
      pmax(0, (sum2_h - N_h * mean_h^2) / (N_h - 1L)),
      0
    )
    S_h <- sqrt(var_h)
    W_h <- N_h / N
    a_h <- .alloc_weights(alloc, q, N_h, S_h, cost_h)
    if (!is.null(take_all_idx)) {
      a_h[take_all_idx] <- 0
    }
    m_h <- pmin(rep(2, L), N_h)
    M_h <- N_h
    if (!is.null(take_all_idx)) {
      m_h[take_all_idx] <- N_h[take_all_idx]
    }

    if (use_cv) {
      ybar <- .aggregate_mean(W_h, mean_h)
      tgt_V <- if (is.null(target_V)) {
        (target_cv * abs(ybar))^2
      } else {
        target_V
      }
      .strata_size_for_target(W_h, S_h, N_h, a_h, m_h, M_h, tgt_V,
                              deff, resp_rate)
    } else {
      if (n_total < sum(m_h) || n_total > sum(M_h)) return(Inf)
      n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
      V <- .strata_variance(W_h, S_h, n_h, N_h, deff, resp_rate)
      ybar <- .aggregate_mean(W_h, mean_h)
      if (ybar == 0) Inf else sqrt(V) / abs(ybar)
    }
  }

  init_bk <- .strata_cumrootf(x_sort, L, warn = FALSE)
  if (length(init_bk) < L - 1L) {
    quant_p <- seq(0, 1, length.out = L + 1L)[2:L]
    init_bk <- unname(quantile(x_sort, probs = quant_p))
  }

  idx_of <- function(bk_) {
    k <- findInterval(bk_, x_uniq, all.inside = TRUE)
    k1 <- k + 1L
    ifelse(abs(x_uniq[k1] - bk_) < abs(x_uniq[k] - bk_), k1, k)
  }

  init_idx <- idx_of(init_bk)
  init_pidx <- c(0L, uniq_pidx[init_idx], N)
  best_obj <- obj_from_pidx(init_pidx)
  best_idx <- init_idx

  total_steps <- as.integer(n_restart) * as.integer(max_iter)
  rand_h <- sample.int(L - 1L, total_steps, replace = TRUE)
  rand_dir <- sample(c(-1L, 1L), total_steps, replace = TRUE)
  ri <- 0L

  for (restart in seq_len(n_restart)) {
    if (restart == 1L) {
      cur_idx <- init_idx
    } else {
      cur_idx <- sort(sample.int(nu - 1L, L - 1L) + 0L)
      cur_idx <- pmin(cur_idx, nu - 1L)
      cur_idx <- pmax(cur_idx, 2L)
    }
    cur_pidx <- c(0L, uniq_pidx[cur_idx], N)
    cur_obj <- obj_from_pidx(cur_pidx)
    if (!is.finite(cur_obj)) {
      cur_obj <- Inf
    }

    for (step in seq_len(max_iter)) {
      ri <- ri + 1L
      h <- rand_h[ri]
      new_h <- cur_idx[h] + rand_dir[ri]

      if (new_h < 2L || new_h >= nu) {
        next
      }
      if (h > 1L && new_h <= cur_idx[h - 1L]) {
        next
      }
      if (h < L - 1L && new_h >= cur_idx[h + 1L]) {
        next
      }

      new_pidx <- cur_pidx
      new_pidx[h + 1L] <- uniq_pidx[new_h]
      new_obj <- obj_from_pidx(new_pidx)

      if (is.finite(new_obj) && new_obj < cur_obj) {
        cur_idx[h] <- new_h
        cur_pidx <- new_pidx
        cur_obj <- new_obj
      }
    }

    if (cur_obj < best_obj) {
      best_idx <- cur_idx
      best_obj <- cur_obj
    }
  }
  list(
    bk = x_uniq[best_idx],
    converged = NA,
    feasible = is.finite(best_obj)
  )
}

#' Largest-remainder rounding constrained to integer bounds
#'
#' Preserves the (rounded) total while keeping every element inside
#' `[lo, hi]`. Errors when no integer vector can satisfy both.
#' @keywords internal
#' @noRd
.round_oric_bounded <- function(x, lo, hi, total = NULL) {
  T <- as.integer(round(if (is.null(total)) sum(x) else total))
  lo <- as.integer(lo)
  hi <- as.integer(hi)
  if (sum(lo) > T || sum(hi) < T) {
    stop(
      sprintf(
        "no integer allocation preserves the total (%d) within the integer bounds (feasible totals: %d to %d); relax the constraints or change the target",
        T, sum(lo), sum(hi)
      ),
      call. = FALSE
    )
  }
  m <- pmin(pmax(as.integer(floor(x)), lo), hi)
  d <- T - sum(m)
  while (d != 0L) {
    if (d > 0L) {
      cand <- which(m < hi)
      j <- cand[which.max(x[cand] - m[cand])]
      m[j] <- m[j] + 1L
      d <- d - 1L
    } else {
      cand <- which(m > lo)
      j <- cand[which.max(m[cand] - x[cand])]
      m[j] <- m[j] - 1L
      d <- d + 1L
    }
  }
  m
}

#' @keywords internal
#' @noRd
.round_oric <- function(x) {
  m <- floor(x)
  frac <- x - m
  k <- as.integer(round(sum(frac)))
  if (k == 0L) {
    return(as.integer(m))
  }
  idx <- order(frac, m, decreasing = TRUE, method = "radix")[seq_len(k)]
  m[idx] <- m[idx] + 1L
  as.integer(m)
}
