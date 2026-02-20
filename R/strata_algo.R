.strata_precompute <- function(x_sort) {
  list(
    cs     = c(0, cumsum(x_sort)),
    cs2    = c(0, cumsum(x_sort^2)),
    n      = length(x_sort),
    x_sort = x_sort
  )
}

.bk_to_idx <- function(x_sort, bk) {
  c(0L, findInterval(bk, x_sort), length(x_sort))
}

.strata_stats_from_prefix <- function(pre, idx) {
  L <- length(idx) - 1L
  N_h <- diff(idx)
  lo <- idx[seq_len(L)] + 1L
  hi <- idx[seq_len(L) + 1L] + 1L
  sum_h  <- pre$cs[hi] - pre$cs[lo]
  sum2_h <- pre$cs2[hi] - pre$cs2[lo]
  mean_h <- ifelse(N_h > 0L, sum_h / N_h, 0)
  var_h  <- ifelse(N_h > 1L, pmax(0, (sum2_h - N_h * mean_h^2) / (N_h - 1L)), 0)
  list(
    N_h    = N_h,
    W_h    = N_h / pre$n,
    S_h    = sqrt(var_h),
    mean_h = mean_h
  )
}

.rna_alloc <- function(a_h, n_total, m_h, M_h) {
  n_h <- numeric(length(a_h))
  fixed <- logical(length(a_h))
  budget <- n_total

  for (pass in seq_len(2L * length(a_h))) {
    act <- which(!fixed)
    if (length(act) == 0L) break
    sa <- sum(a_h[act])
    if (sa <= 0) {
      n_h[act] <- budget / length(act)
      break
    }
    n_h[act] <- budget * a_h[act] / sa
    hi <- act[n_h[act] > M_h[act]]
    lo <- act[n_h[act] < m_h[act]]
    if (length(hi) == 0L && length(lo) == 0L) break
    if (length(hi) > 0L) {
      n_h[hi] <- M_h[hi]
      fixed[hi] <- TRUE
      budget <- budget - sum(M_h[hi])
    } else {
      n_h[lo] <- m_h[lo]
      fixed[lo] <- TRUE
      budget <- budget - sum(m_h[lo])
    }
  }
  n_h
}

#' Evaluate stratified allocation for given boundaries
#' @keywords internal
#' @noRd
.strata_alloc <- function(x, bk, n_total, alloc_q, cost_h, certain_idx = NULL,
                           .pre = NULL) {
  L <- length(bk) + 1L

  if (!is.null(.pre)) {
    idx <- .bk_to_idx(.pre$x_sort, bk)
    stats <- .strata_stats_from_prefix(.pre, idx)
    N_h    <- stats$N_h
    N      <- .pre$n
    W_h    <- stats$W_h
    S_h    <- stats$S_h
    mean_h <- stats$mean_h
    breaks <- c(.pre$x_sort[1L], bk, .pre$x_sort[N])
  } else {
    x_range <- range(x)
    breaks <- c(x_range[1L], bk, x_range[2L])
    bins <- findInterval(x, bk, left.open = TRUE) + 1L
    N_h <- tabulate(bins, nbins = L)
    N <- length(x)
    W_h <- N_h / N

    bins_f <- factor(bins, levels = seq_len(L))
    x_split <- split(x, bins_f)
    S_h <- vapply(x_split, function(xi) {
      if (length(xi) < 2L) return(0)
      sqrt(var(xi))
    }, numeric(1L))
    mean_h <- vapply(x_split, mean, numeric(1L))
    mean_h[is.nan(mean_h)] <- 0
  }

  a_h <- N_h^alloc_q[[1L]] * S_h^alloc_q[[2L]] / cost_h^alloc_q[[3L]]
  m_h <- pmin(rep(2, L), N_h)
  M_h <- N_h
  if (!is.null(certain_idx)) {
    a_h[certain_idx] <- 0
    m_h[certain_idx] <- N_h[certain_idx]
  }

  sa <- sum(a_h)
  if (sa <= 0) {
    n_h <- rep(2, L)
    if (!is.null(certain_idx)) n_h[certain_idx] <- N_h[certain_idx]
    n_h <- n_h * (n_total / sum(n_h))
  } else {
    n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
  }

  V <- sum(W_h^2 * S_h^2 / n_h * (1 - n_h / N_h))
  ybar <- sum(W_h * mean_h)
  cv <- if (ybar == 0) Inf else sqrt(V) / abs(ybar)

  list(
    N_h    = N_h,
    W_h    = W_h,
    S_h    = S_h,
    mean_h = mean_h,
    n_h    = n_h,
    cv     = cv,
    V      = V,
    lower  = breaks[-length(breaks)],
    upper  = breaks[-1L]
  )
}

#' Objective function: CV given boundaries (for minimization)
#' @keywords internal
#' @noRd
.strata_obj <- function(x, bk, n_total, alloc_q, cost_h, certain_idx = NULL,
                         .pre = NULL) {
  res <- .strata_alloc(x, bk, n_total, alloc_q, cost_h, certain_idx, .pre)
  res$cv
}

#' Required n to achieve target CV
#' @keywords internal
#' @noRd
.strata_n_for_cv <- function(x, bk, target_cv, alloc_q, cost_h,
                             certain_idx = NULL, .pre = NULL) {
  L <- length(bk) + 1L

  if (!is.null(.pre)) {
    idx <- .bk_to_idx(.pre$x_sort, bk)
    stats <- .strata_stats_from_prefix(.pre, idx)
    N_h    <- stats$N_h
    N      <- .pre$n
    W_h    <- stats$W_h
    S_h    <- stats$S_h
    mean_h <- stats$mean_h
  } else {
    x_range <- range(x)
    bins <- findInterval(x, bk, left.open = TRUE) + 1L
    N_h <- tabulate(bins, nbins = L)
    N <- length(x)
    W_h <- N_h / N

    bins_f <- factor(bins, levels = seq_len(L))
    x_split <- split(x, bins_f)
    S_h <- vapply(x_split, function(xi) {
      if (length(xi) < 2L) return(0)
      sqrt(var(xi))
    }, numeric(1L))
    mean_h <- vapply(x_split, mean, numeric(1L))
    mean_h[is.nan(mean_h)] <- 0
  }

  ybar <- sum(W_h * mean_h)
  target_V <- (target_cv * abs(ybar))^2

  a_h <- N_h^alloc_q[[1L]] * S_h^alloc_q[[2L]] / cost_h^alloc_q[[3L]]
  m_h <- pmin(rep(2, L), N_h)
  M_h <- N_h
  if (!is.null(certain_idx)) {
    a_h[certain_idx] <- 0
    m_h[certain_idx] <- N_h[certain_idx]
  }

  sa <- sum(a_h)
  if (sa <= 0) return(Inf)

  active <- N_h > 0
  lo <- 2 * L
  hi <- N
  for (i in seq_len(60L)) {
    mid <- (lo + hi) / 2
    n_h <- .rna_alloc(a_h, mid, m_h, M_h)
    V <- sum(W_h[active]^2 * S_h[active]^2 / n_h[active] *
               (1 - n_h[active] / N_h[active]))
    if (V > target_V) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

#' Dalenius-Hodges cumulative sqrt(f) rule
#' @keywords internal
#' @noRd
.strata_cumrootf <- function(x, L, nclass = NULL) {
  if (is.null(nclass)) nclass <- nclass.FD(x)
  nclass <- max(nclass, L)
  h <- hist(x, breaks = nclass, plot = FALSE)
  csf <- cumsum(sqrt(h$counts))
  total <- csf[length(csf)]
  targets <- total * seq_len(L - 1L) / L
  idx <- vapply(targets, function(t) which.min(abs(csf - t)), integer(1L))
  idx <- pmin(idx, length(h$breaks) - 1L)
  unique(h$breaks[idx + 1L])
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

#' Lavallée-Hidiroglou iterative boundary optimization
#' @keywords internal
#' @noRd
.strata_lh <- function(x_sort, L, n_total, target_cv, alloc_q, cost_h,
                        maxiter, certain_idx = NULL) {
  N <- length(x_sort)
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
  obj_fn <- if (use_cv) {
    function(bk_) .strata_n_for_cv(x_sort, bk_, target_cv, alloc_q, cost_h,
                                    certain_idx, .pre = pre)
  } else {
    function(bk_) .strata_obj(x_sort, bk_, n_total, alloc_q, cost_h,
                               certain_idx, .pre = pre)
  }

  best_obj <- obj_fn(bk)
  best_bk <- bk
  converged <- FALSE
  tol <- diff(range(x_sort)) * 1e-6

  for (iter in seq_len(maxiter)) {
    bk_old <- bk
    for (h in seq_len(L - 1L)) {
      lo <- if (h == 1L) x_uniq[1L] + tol else bk[h - 1L] + tol
      hi <- if (h == L - 1L) x_uniq[nu] - tol else bk[h + 1L] - tol
      if (lo >= hi) next
      opt <- optimize(function(b) {
        bk_try <- bk
        bk_try[h] <- b
        obj_fn(bk_try)
      }, interval = c(lo, hi))
      bk[h] <- opt$minimum
    }
    new_obj <- obj_fn(bk)
    if (new_obj < best_obj) {
      best_obj <- new_obj
      best_bk <- bk
    }
    rel_change <- abs(new_obj - best_obj) / (abs(best_obj) + 1e-12)
    if (max(abs(bk - bk_old)) < tol || rel_change < 1e-8) {
      converged <- TRUE
      break
    }
  }
  list(bk = best_bk, converged = converged)
}

#' Kozak random search boundary optimization
#' @keywords internal
#' @noRd
.strata_kozak <- function(x_sort, L, n_total, target_cv, alloc_q, cost_h,
                           maxiter, niter, certain_idx = NULL) {
  x_uniq <- sort(unique(x_sort))
  nu <- length(x_uniq)

  if (nu < L) {
    stop("fewer unique values than requested strata", call. = FALSE)
  }

  pre <- .strata_precompute(x_sort)

  use_cv <- !is.null(target_cv)
  obj_fn <- if (use_cv) {
    function(bk_) .strata_n_for_cv(x_sort, bk_, target_cv, alloc_q, cost_h,
                                    certain_idx, .pre = pre)
  } else {
    function(bk_) .strata_obj(x_sort, bk_, n_total, alloc_q, cost_h,
                               certain_idx, .pre = pre)
  }

  init_bk <- .strata_cumrootf(x_sort, L)
  if (length(init_bk) < L - 1L) {
    quant_p <- seq(0, 1, length.out = L + 1L)[2:L]
    init_bk <- unname(quantile(x_sort, probs = quant_p))
  }

  idx_of <- function(bk_) {
    k <- findInterval(bk_, x_uniq, all.inside = TRUE)
    k1 <- k + 1L
    ifelse(abs(x_uniq[k1] - bk_) < abs(x_uniq[k] - bk_), k1, k)
  }

  best_bk <- init_bk
  best_obj <- obj_fn(init_bk)

  for (restart in seq_len(niter)) {
    if (restart == 1L) {
      cur_idx <- idx_of(init_bk)
    } else {
      cur_idx <- sort(sample.int(nu - 1L, L - 1L) + 0L)
      cur_idx <- pmin(cur_idx, nu - 1L)
      cur_idx <- pmax(cur_idx, 2L)
    }
    cur_bk <- x_uniq[cur_idx]
    cur_obj <- obj_fn(cur_bk)

    for (step in seq_len(maxiter)) {
      h <- sample.int(L - 1L, 1L)
      direction <- sample(c(-1L, 1L), 1L)
      new_idx <- cur_idx
      new_idx[h] <- cur_idx[h] + direction

      if (new_idx[h] < 2L || new_idx[h] >= nu) next
      if (h > 1L && new_idx[h] <= new_idx[h - 1L]) next
      if (h < L - 1L && new_idx[h] >= new_idx[h + 1L]) next

      new_bk <- x_uniq[new_idx]
      new_obj <- obj_fn(new_bk)
      if (new_obj < cur_obj) {
        cur_idx <- new_idx
        cur_bk <- new_bk
        cur_obj <- new_obj
      }
    }

    if (cur_obj < best_obj) {
      best_bk <- cur_bk
      best_obj <- cur_obj
    }
  }
  list(bk = best_bk, converged = TRUE)
}

.round_oric <- function(x) {
  m <- floor(x)
  frac <- x - m
  k <- as.integer(round(sum(frac)))
  if (k == 0L) return(as.integer(m))
  idx <- order(frac, decreasing = TRUE)[seq_len(k)]
  m[idx] <- m[idx] + 1L
  as.integer(m)
}
