#' Resolve free-take two-stage specification for joint allocation
#'
#' Detects whether a `n_alloc_joint()` call requests free-take two-stage
#' optimization (`take_2` omitted or containing `NA`, or take bounds supplied),
#' validates stage counts, stage costs, and bounds, and builds the per-stratum
#' candidate integer take grids.
#' @keywords internal
#' @noRd
.resolve_freetake_spec <- function(
  frame,
  measures,
  psu = NULL,
  take_2_min = NULL,
  take_2_max = NULL,
  take_2_step = NULL,
  unit_cost = NULL,
  fpc = "unit"
) {
  if (!is.data.frame(frame) || nrow(frame) == 0L) {
    return(list(active = FALSE))
  }
  H <- nrow(frame)
  frame_names <- names(frame)
  measure_names <- if (is.data.frame(measures)) names(measures) else character(0)
  has_bound_args <- !is.null(take_2_min) || !is.null(take_2_max) || !is.null(take_2_step)
  has_bound_cols <- any(c("take_2_min", "take_2_max", "take_2_step") %in% frame_names)
  has_bounds <- has_bound_args || has_bound_cols

  resolve_bound <- function(arg, col, default_val, allow_na_default = FALSE) {
    val <- if (!is.null(arg)) {
      if (!is.numeric(arg) || !length(arg) %in% c(1L, H) || anyNA(arg) ||
          any(!is.finite(arg))) {
        stop(sprintf("'%s' must be one positive whole number, or one per stratum",
                     col), call. = FALSE)
      }
      rep_len(as.numeric(arg), H)
    } else if (col %in% frame_names) {
      v <- frame[[col]]
      if (is.logical(v) && all(is.na(v))) v <- as.numeric(v)
      if (!is.numeric(v)) {
        stop(sprintf("'frame$%s' must be numeric", col), call. = FALSE)
      }
      as.numeric(v)
    } else {
      rep(default_val, H)
    }
    check_idx <- if (allow_na_default) which(!is.na(val)) else seq_len(H)
    if (length(check_idx) > 0L) {
      v <- val[check_idx]
      if (any(!is.finite(v)) || any(v < 1) || any(abs(v - round(v)) > 1e-8)) {
        stop(sprintf("'%s' must contain positive whole numbers", col),
             call. = FALSE)
      }
      val[check_idx] <- round(v)
    }
    val
  }

  if (!is.null(psu)) {
    deeper <- grep(
      "^(take_([3-9]|[1-9][0-9]+)|(icc|var_ratio|N)_([2-9]|[1-9][0-9]+)|cost_([3-9]|[1-9][0-9]+))$",
      c(frame_names, measure_names),
      value = TRUE
    )
    if (length(deeper) > 0L) {
      stop(
        sprintf("certainty-aware allocation is two-stage; remove: %s",
                paste(sQuote(unique(deeper)), collapse = ", ")),
        call. = FALSE
      )
    }
    has_take_col <- "take_2" %in% frame_names
    take_raw <- if (has_take_col) frame$take_2 else rep(NA_real_, H)
    if (is.logical(take_raw) && all(is.na(take_raw))) {
      take_raw <- as.numeric(take_raw)
    }
    if (!is.numeric(take_raw)) {
      stop("'take_2' must contain positive whole numbers or NA", call. = FALSE)
    }
    free <- is.na(take_raw)
    if (!any(free)) {
      if (has_bounds) {
        stop(
          "'take_2_min', 'take_2_max', and 'take_2_step' require at least one free stratum with 'take_2 = NA' or 'take_2' omitted",
          call. = FALSE
        )
      }
      return(list(active = FALSE))
    }
    has_cost <- c("cost_1", "cost_2") %in% frame_names
    if (!all(has_cost)) {
      if (!has_take_col && !has_bounds) {
        stop("'psu' requires the frame column 'take_2', the within-PSU take",
             call. = FALSE)
      }
      stop(
        "free-take optimization requires stage costs 'cost_1' and 'cost_2' in 'frame'",
        call. = FALSE
      )
    }
    if ("unit_cost" %in% frame_names || !is.null(unit_cost)) {
      stop(
        "certainty-aware allocation uses stage costs; do not supply 'unit_cost'",
        call. = FALSE
      )
    }
    fixed_idx <- which(!free)
    if (length(fixed_idx) > 0L) {
      tf <- take_raw[fixed_idx]
      if (any(!is.finite(tf)) || any(tf < 1) || any(abs(tf - round(tf)) > 1e-8)) {
        stop("'take_2' must contain positive whole numbers or NA", call. = FALSE)
      }
      take_raw[fixed_idx] <- round(tf)
    }
    t_min <- resolve_bound(take_2_min, "take_2_min", 1, allow_na_default = TRUE)
    t_min[is.na(t_min)] <- 1
    t_max_user <- resolve_bound(take_2_max, "take_2_max", NA_real_, allow_na_default = TRUE)
    t_step <- resolve_bound(take_2_step, "take_2_step", 1, allow_na_default = TRUE)
    t_step[is.na(t_step)] <- 1
    if (any(free & !is.na(t_max_user) & t_max_user < t_min)) {
      stop("'take_2_max' must be at least 'take_2_min' in every free stratum",
           call. = FALSE)
    }

    stratum <- as.character(frame[["stratum"]] %||% seq_len(H))
    psu_stratum <- as.character(psu$stratum)
    flagged <- if ("certainty" %in% names(psu)) {
      as.logical(.check_take_all(psu$certainty, nrow(psu)))
    } else {
      logical(nrow(psu))
    }
    idx_of <- split(seq_len(nrow(psu)), factor(psu_stratum, levels = stratum))
    cand_takes <- vector("list", H)
    m_max_int <- numeric(H)
    for (h in seq_len(H)) {
      i <- idx_of[[h]]
      M_h <- psu$N[i]
      fl_h <- flagged[i]
      rem_M <- M_h[!fl_h]
      if (!free[h]) {
        b_h <- take_raw[h]
        if (length(rem_M) > 0L && any(rem_M < b_h)) {
          small <- i[!fl_h][rem_M < b_h]
          label <- if (is.null(psu[["psu_id"]])) sprintf("row %d", small) else
            as.character(psu$psu_id[small])
          shown <- paste(sQuote(utils::head(label, 5L)), collapse = ", ")
          if (length(label) > 5L) {
            shown <- sprintf("%s and %d more", shown, length(label) - 5L)
          }
          stop(
            sprintf(
              "%d PSU%s in 'psu' hold%s fewer units than the take 'take_2' of %s stratum (%s). Merge each with a neighboring PSU, as merge_psus() does, or flag it in 'psu$certainty'",
              length(label), if (length(label) == 1L) "" else "s",
              if (length(label) == 1L) "s" else "",
              if (length(label) == 1L) "its" else "their", shown
            ),
            call. = FALSE
          )
        }
        m_max_int[h] <- b_h
        cand_takes[[h]] <- as.integer(b_h)
      } else {
        if (length(rem_M) == 0L) {
          m_max_int[h] <- t_min[h]
          cand_takes[[h]] <- as.integer(t_min[h])
        } else {
          min_rem <- floor(min(rem_M))
          if (t_min[h] > min_rem) {
            small <- i[!fl_h][rem_M < t_min[h]]
            label <- if (is.null(psu[["psu_id"]])) sprintf("row %d", small) else
              as.character(psu$psu_id[small])
            shown <- paste(sQuote(utils::head(label, 5L)), collapse = ", ")
            if (length(label) > 5L) {
              shown <- sprintf("%s and %d more", shown, length(label) - 5L)
            }
            stop(
              sprintf(
                "%d PSU%s in 'psu' hold%s fewer units than 'take_2_min' (%d) in stratum '%s' (%s). Lower 'take_2_min', merge with merge_psus(), or flag in 'psu$certainty'",
                length(label), if (length(label) == 1L) "" else "s",
                if (length(label) == 1L) "s" else "",
                as.integer(t_min[h]), stratum[h], shown
              ),
              call. = FALSE
            )
          }
          cap_h <- if (is.na(t_max_user[h])) min_rem else min(t_max_user[h], min_rem)
          m_max_int[h] <- cap_h
          cand_takes[[h]] <- seq.int(as.integer(t_min[h]), as.integer(cap_h),
                                     by = as.integer(t_step[h]))
        }
      }
    }
    return(list(
      active = TRUE,
      psu = TRUE,
      free = free,
      take_fixed = take_raw,
      take_2_min = ifelse(free, t_min, take_raw),
      take_2_max_int = m_max_int,
      take_2_max_cont = m_max_int,
      take_2_step = ifelse(free, t_step, 1),
      cand_takes = cand_takes
    ))
  }

  take_cols <- grep("^take_[0-9]+$", frame_names, value = TRUE)
  cost_cols <- grep("^cost_[0-9]+$", frame_names, value = TRUE)
  sub_cols <- unique(c(
    grep("^(N|icc|var_ratio)_[0-9]+$", frame_names, value = TRUE),
    grep("^(icc|var_ratio)_[0-9]+$", measure_names, value = TRUE)
  ))
  stage_hits <- integer(0)
  if (length(take_cols) > 0L) {
    stage_hits <- c(stage_hits, as.integer(sub("^take_", "", take_cols)))
  }
  if (length(cost_cols) > 0L) {
    stage_hits <- c(stage_hits, as.integer(sub("^cost_", "", cost_cols)))
  }
  if (length(sub_cols) > 0L) {
    stage_hits <- c(stage_hits, as.integer(sub("^(N|icc|var_ratio)_", "", sub_cols)) + 1L)
  }

  if (length(stage_hits) == 0L) {
    if (has_bounds) {
      stop(
        "'take_2_min', 'take_2_max', and 'take_2_step' apply only to two-stage designs",
        call. = FALSE
      )
    }
    return(list(active = FALSE))
  }

  stages <- max(2L, max(stage_hits))
  if (stages > 2L) {
    any_free_3p <- has_bounds ||
      any(vapply(paste0("take_", 2L:stages), function(col) {
        !col %in% frame_names || anyNA(frame[[col]])
      }, logical(1L)))
    if (any_free_3p) {
      has_na_take <- any(vapply(intersect(paste0("take_", 2L:stages), frame_names),
                                function(col) anyNA(frame[[col]]), logical(1L)))
      if (has_bounds || has_na_take) {
        stop(
          "free-take optimization in n_alloc_joint() supports two-stage designs; supply fixed takes for three or more stages",
          call. = FALSE
        )
      }
    }
    return(list(active = FALSE))
  }

  has_take_col <- "take_2" %in% frame_names
  take_raw <- if (has_take_col) frame$take_2 else rep(NA_real_, H)
  if (is.logical(take_raw) && all(is.na(take_raw))) {
    take_raw <- as.numeric(take_raw)
  }
  if (!is.numeric(take_raw)) {
    stop("'take_2' must contain positive whole numbers or NA", call. = FALSE)
  }
  free <- is.na(take_raw)
  if (!any(free)) {
    if (has_bounds) {
      stop(
        "'take_2_min', 'take_2_max', and 'take_2_step' require at least one free stratum with 'take_2 = NA' or 'take_2' omitted",
        call. = FALSE
      )
    }
    return(list(active = FALSE))
  }

  has_N1 <- "N_1" %in% frame_names
  has_c1 <- "cost_1" %in% frame_names
  has_c2 <- "cost_2" %in% frame_names
  if (!has_take_col && !has_bounds && !(has_N1 && has_c1 && has_c2)) {
    return(list(active = FALSE))
  }
  if (!has_c1 || !has_c2) {
    stop(
      "free-take two-stage allocation requires stage costs 'cost_1' and 'cost_2' in 'frame'",
      call. = FALSE
    )
  }
  if (!has_N1) {
    stop(
      "free-take two-stage allocation requires 'N_1' in 'frame' (or a 'psu' register)",
      call. = FALSE
    )
  }
  if ("unit_cost" %in% frame_names || !is.null(unit_cost)) {
    stop(
      "two-stage allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }

  fixed_idx <- which(!free)
  if (length(fixed_idx) > 0L) {
    tf <- take_raw[fixed_idx]
    if (any(!is.finite(tf)) || any(tf < 1) || any(abs(tf - round(tf)) > 1e-8)) {
      stop("'take_2' must contain positive whole numbers or NA", call. = FALSE)
    }
    take_raw[fixed_idx] <- round(tf)
  }

  N <- frame$N
  N_1 <- frame$N_1
  if (!is.numeric(N) || anyNA(N) || any(!is.finite(N)) || any(N <= 0)) {
    stop("'N' must contain positive finite values", call. = FALSE)
  }
  if (!is.numeric(N_1) || anyNA(N_1) || any(!is.finite(N_1)) || any(N_1 < 1) ||
      any(abs(N_1 - round(N_1)) > 1e-8)) {
    stop("'N_1' must contain positive whole numbers", call. = FALSE)
  }
  N_1 <- as.numeric(round(N_1))
  if (any(N_1 > N + 1e-8)) {
    stop("'N_1' must not exceed the ultimate-unit population 'N'", call. = FALSE)
  }
  M_bar <- N / N_1

  t_min <- resolve_bound(take_2_min, "take_2_min", 1, allow_na_default = TRUE)
  t_min[is.na(t_min)] <- 1
  t_max_user <- resolve_bound(take_2_max, "take_2_max", NA_real_, allow_na_default = TRUE)
  t_step <- resolve_bound(take_2_step, "take_2_step", 1, allow_na_default = TRUE)
  t_step[is.na(t_step)] <- 1

  if (any(free & !is.na(t_max_user) & t_max_user < t_min)) {
    stop("'take_2_max' must be at least 'take_2_min' in every free stratum",
         call. = FALSE)
  }
  if (identical(fpc, "stage") &&
      any(free & !is.na(t_max_user) & t_max_user > M_bar + 1e-8)) {
    stop(
      "under fpc = \"stage\", 'take_2_max' cannot exceed the average PSU size 'N / N_1'",
      call. = FALSE
    )
  }
  if (any(free & !is.na(t_max_user) & t_max_user > N + 1e-8)) {
    stop("'take_2_max' cannot exceed the stratum population 'N'", call. = FALSE)
  }

  m_max_cont <- ifelse(free, ifelse(is.na(t_max_user), M_bar, t_max_user), take_raw)
  m_max_int <- ifelse(free, ifelse(is.na(t_max_user), floor(M_bar + 1e-9), floor(t_max_user + 1e-9)), take_raw)
  m_min <- ifelse(free, t_min, take_raw)
  m_step <- ifelse(free, t_step, 1)

  if (any(free & m_min > m_max_int)) {
    stop(
      "'take_2_min' exceeds the maximum feasible take in at least one stratum",
      call. = FALSE
    )
  }

  cand_takes <- lapply(seq_len(H), function(h) {
    if (!free[h]) {
      as.integer(take_raw[h])
    } else {
      seq.int(as.integer(m_min[h]), as.integer(m_max_int[h]),
              by = as.integer(m_step[h]))
    }
  })

  list(
    active = TRUE,
    psu = FALSE,
    free = free,
    take_fixed = take_raw,
    take_2_min = m_min,
    take_2_max_cont = m_max_cont,
    take_2_max_int = m_max_int,
    take_2_step = m_step,
    cand_takes = cand_takes
  )
}

#' Exact global minimum of one stratum's 2-stage Lagrangian on its polygon
#'
#' Minimizes `c1 * x + c2 * y + A / x + B / y - D * (y / x)` over
#' `1 <= x <= N1`, `L_ult <= y <= N`, `m_min * x <= y <= m_max * x`,
#' where `x = n_1` and `y = n_1 * take_2`. Handles fixed takes
#' (`m_min == m_max`), `D = 0` (`fpc = "none"` or `"stage"`, including `A <= 0`),
#' and `D > 0` (`fpc = "unit"`).
#' @keywords internal
#' @noRd
.freetake_stratum_min <- function(c1, c2, A, B, D, N1, N, L_ult, m_min, m_max) {
  best_val <- Inf
  best_x <- 1
  best_y <- max(L_ult, m_min)

  eval_pt <- function(x, y) {
    if (!is.finite(x) || !is.finite(y) || x <= 0 || y <= 0) return(invisible(NULL))
    val <- c1 * x + c2 * y + A / x + B / y - D * (y / x)
    if (is.finite(val) && val < best_val) {
      best_val <<- val
      best_x <<- x
      best_y <<- y
    }
  }

  min_1d <- function(p, q, lo, hi) {
    if (lo > hi + 1e-12) return(numeric(0))
    if (lo > hi) lo <- hi <- 0.5 * (lo + hi)
    pts <- c(lo, hi)
    if (p > 0 && q > 0) {
      t0 <- sqrt(q / p)
      if (t0 > lo && t0 < hi) pts <- c(pts, t0)
    }
    pts
  }

  if (abs(m_max - m_min) <= 1e-12) {
    m0 <- 0.5 * (m_min + m_max)
    x_lo <- max(1, L_ult / m0)
    x_hi <- min(N1, N / m0)
    p <- c1 + c2 * m0
    q <- max(0, A + B / m0)
    for (x in min_1d(p, q, x_lo, x_hi)) eval_pt(x, m0 * x)
    return(list(x = best_x, y = best_y, val = best_val))
  }

  {
    x_lo <- max(1, L_ult / m_min)
    x_hi <- min(N1, N / m_min)
    for (x in min_1d(c1 + c2 * m_min, A + B / m_min, x_lo, x_hi)) {
      eval_pt(x, m_min * x)
    }
  }
  {
    x_lo <- max(1, L_ult / m_max)
    x_hi <- min(N1, N / m_max)
    for (x in min_1d(c1 + c2 * m_max, A + B / m_max, x_lo, x_hi)) {
      eval_pt(x, m_max * x)
    }
  }
  {
    y_lo <- max(L_ult, m_min)
    y_hi <- min(N, m_max)
    for (y in min_1d(c2 - D, B, y_lo, y_hi)) eval_pt(1, y)
  }
  {
    y_lo <- max(L_ult, N1 * m_min)
    y_hi <- min(N, N1 * m_max)
    for (y in min_1d(c2 - D / N1, B, y_lo, y_hi)) eval_pt(N1, y)
  }
  {
    x_lo <- max(1, L_ult / m_max)
    x_hi <- min(N1, L_ult / m_min)
    for (x in min_1d(c1, A - D * L_ult, x_lo, x_hi)) eval_pt(x, L_ult)
  }
  {
    x_lo <- max(1, N / m_max)
    x_hi <- min(N1, N / m_min)
    for (x in min_1d(c1, A - D * N, x_lo, x_hi)) eval_pt(x, N)
  }

  if (A > 0 && B > 0) {
    if (D <= 0) {
      x0 <- sqrt(A / c1)
      y0 <- sqrt(B / c2)
      if (x0 >= 1 - 1e-12 && x0 <= N1 + 1e-12 &&
          y0 >= L_ult - 1e-12 && y0 <= N + 1e-12 &&
          y0 >= m_min * x0 - 1e-12 && y0 <= m_max * x0 + 1e-12) {
        eval_pt(
          min(max(x0, 1), N1),
          min(max(y0, max(L_ult, m_min * min(max(x0, 1), N1))),
              min(N, m_max * min(max(x0, 1), N1)))
        )
      }
    } else {
      kappa <- (D^2 * sqrt(c1) / (4 * B))^(1 / 3)
      root <- sqrt(D^2 + 4 * kappa^2 * A)
      ## Each form subtracts nothing close to itself on its own side of 0.
      y_crit <- if (D >= 0) 2 * A / (D + root) else (root - D) / (2 * kappa^2)
      y_lo_root <- sqrt(B / c2)
      if (y_crit > y_lo_root && A - D * y_crit > 0) {
        h_fn <- function(y) B / (y * y) + D * sqrt(c1) / sqrt(A - D * y) - c2
        if (h_fn(y_crit) <= 0) {
          lo_r <- y_lo_root
          hi_r <- y_crit
          y0 <- 0.5 * (lo_r + hi_r)
          for (it in seq_len(40L)) {
            rem <- A - D * y0
            if (rem <= 0) {
              y0 <- 0.5 * (lo_r + hi_r)
              rem <- A - D * y0
            }
            hv <- B / (y0 * y0) + D * sqrt(c1) / sqrt(rem) - c2
            if (hv > 0) lo_r <- y0 else hi_r <- y0
            if (abs(hv) <= 1e-14 * max(1, c2) || hi_r - lo_r <= 1e-13 * max(1, y0)) {
              break
            }
            hp <- -2 * B / (y0^3) + 0.5 * D^2 * sqrt(c1) / (rem * sqrt(rem))
            cand <- if (hp < -1e-30) y0 - hv / hp else 0.5 * (lo_r + hi_r)
            if (cand <= lo_r || cand >= hi_r) cand <- 0.5 * (lo_r + hi_r)
            y0 <- cand
          }
          if (A - D * y0 > 0) {
            x0 <- sqrt((A - D * y0) / c1)
            if (x0 >= 1 - 1e-10 && x0 <= N1 + 1e-10 &&
                y0 >= L_ult - 1e-10 && y0 <= N + 1e-10 &&
                y0 >= m_min * x0 - 1e-10 && y0 <= m_max * x0 + 1e-10) {
              x0 <- min(max(x0, 1), N1)
              y0 <- min(max(y0, max(L_ult, m_min * x0)), min(N, m_max * x0))
              eval_pt(x0, y0)
            }
          }
        }
      }
    }
  }

  list(x = best_x, y = best_y, val = best_val)
}

#' Build a free-take 2-stage joint allocation problem
#' @keywords internal
#' @noRd
.freetake_build_problem <- function(
  frame,
  measures,
  targets,
  ft_spec,
  deff = 1,
  resp_rate = 1,
  eligible_per_unit = 1,
  eligibility_rate = 1,
  alpha = 0.05,
  df = NULL,
  min_n_stratum = NULL,
  objective = NULL,
  fpc = "unit"
) {
  check_deff(deff)
  check_resp_rate(resp_rate)
  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)
  check_alpha(alpha)
  if (!is.null(df)) check_df(df)
  min_n_stratum <- .resolve_min_n_stratum(min_n_stratum, frame)

  required_frame <- c("stratum", "N", "N_1", "cost_1", "cost_2")
  missing_frame <- setdiff(required_frame, names(frame))
  if (length(missing_frame) > 0L) {
    stop(sprintf("free-take two-stage allocation frame must contain: %s",
                 paste(sQuote(missing_frame), collapse = ", ")),
         call. = FALSE)
  }
  stratum <- as.character(frame$stratum)
  if (anyNA(stratum) || any(!nzchar(stratum)) || anyDuplicated(stratum)) {
    stop("'frame$stratum' must contain unique non-empty values", call. = FALSE)
  }
  H <- nrow(frame)
  N <- as.numeric(frame$N)
  N_1 <- as.numeric(round(frame$N_1))
  c1 <- as.numeric(frame$cost_1)
  c2 <- as.numeric(frame$cost_2)
  if (anyNA(c1) || any(!is.finite(c1)) || any(c1 <= 0) ||
      anyNA(c2) || any(!is.finite(c2)) || any(c2 <= 0)) {
    stop("'cost_1' and 'cost_2' must contain positive finite values",
         call. = FALSE)
  }
  if ("cost" %in% names(frame)) {
    stop("the per-stratum cost columns are 'cost_1' and 'cost_2'; 'cost' is the total field cost",
         call. = FALSE)
  }
  max_weight <- if ("max_weight" %in% names(frame)) frame$max_weight else rep(NA_real_, H)
  if (!is.numeric(max_weight)) stop("'max_weight' must be numeric", call. = FALSE)
  if (any(!is.na(max_weight) & (!is.finite(max_weight) | max_weight < 1))) {
    stop("'max_weight' must be >= 1 and finite when provided", call. = FALSE)
  }
  take_all <- .check_take_all(frame[["take_all"]], H)
  if (any(take_all)) {
    stop(
      "'take_all' is not supported for multistage allocation because taking every PSU does not imply an ultimate-unit census",
      call. = FALSE
    )
  }
  indicator_cols <- intersect(c("p", "mean", "sd", "var"), names(frame))
  if (length(indicator_cols) > 0L) {
    stop(
      sprintf(
        "joint allocation reads indicator values from 'measures', one row per stratum and indicator. Move %s from 'frame' to 'measures'",
        paste(sQuote(indicator_cols), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ultimate_bounds <- .alloc_bounds(N, max_weight, rep(FALSE, H), min_n_stratum)
  L_ult <- ultimate_bounds$m_h
  y_max <- pmin(N, N_1 * ft_spec$take_2_max_cont)
  y_min <- pmax(L_ult, ft_spec$take_2_min)
  if (any(y_min > y_max + 1e-8)) {
    bad <- stratum[y_min > y_max + 1e-8]
    stop(
      sprintf("stratum bounds and take bounds are infeasible for stratum: %s",
              paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }

  required <- .bethel_required_measures(
    frame, measures, targets, objective, alpha, df
  )
  measures <- required$measures
  target_spec <- required$target_spec
  objective_spec <- required$objective_spec
  req_indices <- required$req_indices
  req_name <- required$req_name
  m_stratum <- as.character(measures$stratum)
  m_name <- as.character(measures$name)
  measure_key <- .bethel_key(m_stratum, m_name)

  mcol <- function(name) {
    if (name %in% names(measures)) measures[[name]] else rep(NA_real_, nrow(measures))
  }
  p <- mcol("p")
  mean_value <- mcol("mean")
  var_value <- mcol("var")
  sd_value <- mcol("sd")
  for (z in list(p = p, mean = mean_value, var = var_value, sd = sd_value)) {
    if (!is.numeric(z)) stop("measure statistics must be numeric", call. = FALSE)
  }
  has_p <- !is.na(p)
  has_mean <- !is.na(mean_value)
  has_var <- !is.na(var_value)
  has_sd <- !is.na(sd_value)
  if (any(!((has_p & !has_mean & !has_var & !has_sd) |
            (!has_p & has_mean & xor(has_var, has_sd))))) {
    stop("each measure row must specify either 'p', or 'mean' and exactly one of 'var'/'sd'",
         call. = FALSE)
  }
  if (any(has_p & (!is.finite(p) | p < 0 | p > 1))) {
    stop("'p' must contain values in [0, 1]", call. = FALSE)
  }
  if (any(has_mean & !is.finite(mean_value))) {
    stop("'mean' must contain finite values", call. = FALSE)
  }
  if (any(has_var & (!is.finite(var_value) | var_value < 0)) ||
      any(has_sd & (!is.finite(sd_value) | sd_value < 0))) {
    stop("'var' and 'sd' must contain non-negative finite values", call. = FALSE)
  }
  mean_norm <- ifelse(has_p, p, mean_value)
  var_norm <- ifelse(has_p, p * (1 - p), ifelse(has_var, var_value, sd_value^2))
  elig_norm <- .joint_row_value("eligible_per_unit", measures, frame, eligible_per_unit)
  scope_norm <- .joint_row_value("eligibility_rate", measures, frame, eligibility_rate)
  yield_norm <- scope_norm * elig_norm
  deff_norm <- .joint_row_value("deff", measures, frame, deff)
  resp_norm <- .joint_row_value("resp_rate", measures, frame, resp_rate)

  rr_present <- unique(c(
    grep("^resp_rate_[0-9]+$", names(frame), value = TRUE),
    grep("^resp_rate_[0-9]+$", names(measures), value = TRUE)
  ))
  if (length(rr_present) > 0L) {
    rr_idx <- as.integer(sub("^resp_rate_", "", rr_present))
    bad_rr <- rr_present[rr_idx < 1L | rr_idx >= 2L]
    .bethel_reject_stage_rate(bad_rr, measures, frame, 2L)
  }
  icc_1 <- .bethel_stage_measure(
    "icc_1", measures, frame, m_stratum, stratum,
    required = TRUE, interval = c(0, 1)
  )
  vr_1 <- .bethel_stage_measure(
    "var_ratio_1", measures, frame, m_stratum, stratum,
    default = 1, positive = TRUE
  )
  rr_1 <- .bethel_stage_measure(
    "resp_rate_1", measures, frame, m_stratum, stratum,
    default = 1, interval = c(0, 1), positive = TRUE
  )

  K <- length(target_spec$name)
  J <- length(objective_spec$name)
  R <- K + J
  membership <- cbind(target_spec$membership, objective_spec$membership)

  a_all <- matrix(0, H, R)
  b_all <- matrix(0, H, R)
  d_all <- matrix(0, H, R)
  e_all <- matrix(0, H, R)
  total_all <- numeric(R)
  domain_N_all <- numeric(R)
  mean_all <- matrix(NA_real_, H, R)
  var_all <- matrix(NA_real_, H, R)
  deff_all <- matrix(NA_real_, H, R)
  resp_all <- matrix(NA_real_, H, R)
  elig_all <- matrix(NA_real_, H, R)
  scope_all <- matrix(NA_real_, H, R)

  for (k in seq_len(R)) {
    idx <- req_indices[[k]]
    wanted <- .bethel_key(stratum[idx], rep(req_name[k], length(idx)))
    mi <- match(wanted, measure_key)
    mean_all[idx, k] <- mean_norm[mi]
    var_all[idx, k] <- var_norm[mi]
    deff_all[idx, k] <- deff_norm[mi]
    resp_all[idx, k] <- resp_norm[mi]
    elig_all[idx, k] <- elig_norm[mi]
    scope_all[idx, k] <- scope_norm[mi]
    N_k <- N[idx] * yield_norm[mi]
    r1 <- rr_1[mi]
    r2 <- resp_norm[mi] * yield_norm[mi]
    v_k <- var_norm[mi] * deff_norm[mi]
    R1 <- vr_1[mi] * icc_1[mi]
    R2 <- vr_1[mi] * (1 - icc_1[mi])
    M_k <- N_k / N_1[idx]

    if (identical(fpc, "none")) {
      a_all[idx, k] <- N_k^2 * v_k * R1 / r1
      b_all[idx, k] <- N_k^2 * v_k * R2 / (r1 * r2)
      d_all[idx, k] <- 0
      e_all[idx, k] <- 0
    } else if (identical(fpc, "stage")) {
      a_all[idx, k] <- N_k^2 * v_k * (R1 - R2 / M_k) / r1
      b_all[idx, k] <- N_k^2 * v_k * R2 / (r1 * r2)
      d_all[idx, k] <- 0
      e_all[idx, k] <- N_k * v_k * R1 * M_k
    } else {
      a_all[idx, k] <- N_k^2 * v_k * R1 / r1
      b_all[idx, k] <- N_k^2 * v_k * R2 / (r1 * r2)
      d_all[idx, k] <- N_k * v_k * R1 * r2
      e_all[idx, k] <- N_k * v_k * R2
    }
    total_all[k] <- sum(N_k * mean_norm[mi])
    domain_N_all[k] <- sum(N_k)
  }

  m_mid <- 0.5 * (ft_spec$take_2_min + ft_spec$take_2_max_cont)
  pos_var <- colSums(a_all + b_all / m_mid) > 0
  if (K > 0L && any(!pos_var[seq_len(K)])) {
    stop(sprintf("precision constraints have no positive variance: %s",
                 paste(target_spec$constraint[!pos_var[seq_len(K)]], collapse = ", ")),
         call. = FALSE)
  }
  if (J > 0L && any(!pos_var[K + seq_len(J)])) {
    stop(sprintf("objective components have no positive variance: %s",
                 paste(objective_spec$component[!pos_var[K + seq_len(J)]],
                       collapse = ", ")),
         call. = FALSE)
  }

  negligible_total <- function(k) {
    in_k <- membership[, k]
    scale_total <- sum(abs(N[in_k] * scope_all[in_k, k] * elig_all[in_k, k] *
                             mean_all[in_k, k]))
    abs(total_all[k]) <= sqrt(.Machine$double.eps) * scale_total
  }

  Vmax <- numeric(K)
  q_target <- .q_alpha(target_spec$alpha, target_spec$df)
  for (k in seq_len(K)) {
    if (target_spec$metric[k] == "moe") {
      Vmax[k] <- (domain_N_all[k] * target_spec$target[k] / q_target[k])^2
      next
    }
    if (negligible_total(k)) {
      stop(sprintf("%s is undefined for constraint '%s' with zero or negligible total",
                   if (target_spec$metric[k] == "cv") "CV" else "relative margin of error",
                   target_spec$constraint[k]), call. = FALSE)
    }
    scaled <- if (target_spec$metric[k] == "cv") target_spec$target[k] else
      target_spec$target[k] / q_target[k]
    Vmax[k] <- (scaled * abs(total_all[k]))^2
  }

  objective_part <- NULL
  if (J > 0L) {
    for (j in seq_len(J)) {
      if (negligible_total(K + j)) {
        stop(sprintf("the objective is undefined for component '%s' with zero or negligible total",
                     objective_spec$component[j]), call. = FALSE)
      }
    }
    total_obj <- total_all[K + seq_len(J)]
    priority <- objective_spec$priority
    w_obj <- priority / total_obj^2
    a_obj <- a_all[, K + seq_len(J), drop = FALSE]
    b_obj <- b_all[, K + seq_len(J), drop = FALSE]
    d_obj <- d_all[, K + seq_len(J), drop = FALSE]
    e_obj <- e_all[, K + seq_len(J), drop = FALSE]
    objective_part <- list(
      component = objective_spec$component,
      name = objective_spec$name,
      domain = objective_spec$domain,
      level = objective_spec$level,
      priority = priority,
      total = total_obj,
      domain_N = domain_N_all[K + seq_len(J)],
      a_mat = a_obj,
      b_mat = b_obj,
      d_mat = d_obj,
      e_mat = e_obj,
      a_h = as.numeric(a_obj %*% w_obj),
      b_h = as.numeric(b_obj %*% w_obj),
      d_h = as.numeric(d_obj %*% w_obj),
      e_h = as.numeric(e_obj %*% w_obj),
      E0 = sum(e_obj %*% w_obj),
      spec = objective_spec$objective
    )
    if (sum(objective_part$a_h + objective_part$b_h / m_mid) <= 0) {
      stop("the objective has no positive-priority component with positive variance",
           call. = FALSE)
    }
  }

  meta <- data.frame(
    constraint = target_spec$constraint,
    name = target_spec$name,
    domain = target_spec$domain,
    level = target_spec$level,
    .metric = target_spec$metric,
    .target = target_spec$target,
    stringsAsFactors = FALSE
  )

  list(
    H = H,
    K = K,
    J = J,
    stratum_ids = stratum,
    constraint_ids = target_spec$constraint,
    constraint_meta = meta,
    N = N,
    N_1 = N_1,
    c1 = c1,
    c2 = c2,
    L_ult = L_ult,
    ft_spec = ft_spec,
    fpc = fpc,
    a = a_all[, seq_len(K), drop = FALSE],
    b = b_all[, seq_len(K), drop = FALSE],
    d = d_all[, seq_len(K), drop = FALSE],
    e = e_all[, seq_len(K), drop = FALSE],
    Vmax = Vmax,
    total = total_all[seq_len(K)],
    domain_N = domain_N_all[seq_len(K)],
    alpha = target_spec$alpha,
    df = target_spec$df,
    q_target = q_target,
    objective = objective_part,
    membership = membership[, seq_len(K), drop = FALSE],
    mean_hk = mean_all[, seq_len(K), drop = FALSE],
    var_hk = var_all[, seq_len(K), drop = FALSE],
    deff_hk = deff_all[, seq_len(K), drop = FALSE],
    resp_hk = resp_all[, seq_len(K), drop = FALSE],
    elig_hk = elig_all[, seq_len(K), drop = FALSE],
    scope_hk = scope_all[, seq_len(K), drop = FALSE],
    rates = .joint_rates_table(measures, m_stratum, resp_norm, scope_norm,
                               matrix(rr_1, ncol = 1L)),
    frame = frame,
    measures = measures,
    targets = target_spec$targets
  )
}

#' Evaluate variance vector for a free-take 2-stage design `(x, y)`
#' @keywords internal
#' @noRd
.freetake_var <- function(a, b, d, e, x, y) {
  m <- y / x
  v <- drop(crossprod(a, 1 / x)) + drop(crossprod(b, 1 / y)) -
    drop(crossprod(d, m)) - colSums(e)
  scale <- pmax(
    abs(drop(crossprod(a, 1 / x))) + abs(drop(crossprod(b, 1 / y))),
    abs(drop(crossprod(d, m))) + abs(colSums(e))
  )
  tiny <- v < 0 & abs(v) <= 100 * .Machine$double.eps * scale
  v[tiny] <- 0
  v
}

#' Solve a minimum-cost continuous free-take 2-stage problem via dual L-BFGS-B
#' @keywords internal
#' @noRd
.freetake_solve_core <- function(
  a, b, d, e, Vmax, c1, c2, N1, N, L_ult, m_min, m_max,
  constraint_ids, control = .bethel_control(), accept = NULL
) {
  H <- length(c1)
  K <- ncol(a)
  if (is.null(accept)) accept <- rep((1 + control$tolerance)^2 - 1, K)

  eval_at_xy <- function(A_h, B_h, D_h, c1_s, c2_s) {
    x <- numeric(H)
    y <- numeric(H)
    val <- 0
    for (h in seq_len(H)) {
      res <- .freetake_stratum_min(
        c1_s[h], c2_s[h], A_h[h], B_h[h], D_h[h],
        N1[h], N[h], L_ult[h], m_min[h], m_max[h]
      )
      x[h] <- res$x
      y[h] <- res$y
      val <- val + res$val
    }
    list(x = x, y = y, val = val)
  }

  infeasible <- integer(0)
  for (k in seq_len(K)) {
    res_min <- eval_at_xy(a[, k], b[, k], d[, k], rep(1e-30, H), rep(1e-30, H))
    vk_min <- .freetake_var(a[, k, drop = FALSE], b[, k, drop = FALSE],
                            d[, k, drop = FALSE], e[, k, drop = FALSE],
                            res_min$x, res_min$y)
    if (vk_min > Vmax[k] * (1 + accept[k])) {
      infeasible <- c(infeasible, k)
    }
  }
  if (length(infeasible) > 0L) {
    return(list(
      classification = "infeasible",
      converged = FALSE,
      infeasible = constraint_ids[infeasible],
      message = "precision targets are unattainable within the population and take bounds"
    ))
  }

  min_cost_pt <- eval_at_xy(numeric(H), numeric(H), numeric(H), c1, c2)
  v_lo <- .freetake_var(a, b, d, e, min_cost_pt$x, min_cost_pt$y)
  if (all(v_lo <= Vmax * (1 + accept))) {
    return(list(
      classification = "optimal",
      converged = TRUE,
      x = min_cost_pt$x,
      y = min_cost_pt$y,
      cost = sum(c1 * min_cost_pt$x + c2 * min_cost_pt$y),
      lambda = numeric(K),
      iterations = 0L,
      primal_residual = max(0, max(v_lo / Vmax - 1)),
      relative_duality_gap = 0,
      dual_value = sum(c1 * min_cost_pt$x + c2 * min_cost_pt$y),
      message = "lower allocation satisfies every precision target",
      feasibility_tolerance = control$tolerance,
      kkt_tolerance = control$kkt_tolerance
    ))
  }

  m_ref <- sqrt(m_min * m_max)
  cost_scale <- stats::median(c1 + c2 * m_ref)
  if (!is.finite(cost_scale) || cost_scale <= 0) cost_scale <- 1
  c1_s <- c1 / cost_scale
  c2_s <- c2 / cost_scale
  a_s <- sweep(a, 2L, Vmax, "/")
  b_s <- sweep(b, 2L, Vmax, "/")
  d_s <- sweep(d, 2L, Vmax, "/")
  e_s <- sweep(e, 2L, Vmax, "/")
  E_col <- colSums(e_s)

  dual_eval <- function(lam) {
    A_h <- drop(a_s %*% lam)
    B_h <- drop(b_s %*% lam)
    D_h <- drop(d_s %*% lam)
    pt <- eval_at_xy(A_h, B_h, D_h, c1_s, c2_s)
    dual <- pt$val - sum(lam * (1 + E_col))
    resid <- .freetake_var(a_s, b_s, d_s, e_s, pt$x, pt$y) - 1
    list(x = pt$x, y = pt$y, dual = dual, grad = resid)
  }

  lam0 <- numeric(K)
  for (k in seq_len(K)) {
    res_k <- function(lk) {
      pt <- eval_at_xy(a_s[, k] * lk, b_s[, k] * lk, d_s[, k] * lk, c1_s, c2_s)
      .freetake_var(a_s[, k, drop = FALSE], b_s[, k, drop = FALSE],
                    d_s[, k, drop = FALSE], e_s[, k, drop = FALSE],
                    pt$x, pt$y) - 1
    }
    if (res_k(0) <= 1e-8) next
    hi <- 1
    for (i in seq_len(60L)) {
      if (res_k(hi) <= 0) break
      hi <- hi * 10
    }
    if (is.finite(hi) && res_k(hi) <= 0) {
      lam0[k] <- stats::uniroot(res_k, c(0, hi), tol = 1e-8)$root
    } else {
      lam0[k] <- hi
    }
  }
  lam0 <- lam0 / max(1L, K)

  cache_lam <- NULL
  cache_res <- NULL
  get_eval <- function(lam) {
    if (!is.null(cache_lam) && identical(lam, cache_lam)) return(cache_res)
    cache_lam <<- lam
    cache_res <<- dual_eval(lam)
    cache_res
  }
  fn <- function(lam) -get_eval(lam)$dual
  gr <- function(lam) -get_eval(lam)$grad

  opt <- stats::optim(
    par = lam0,
    fn = fn,
    gr = gr,
    method = "L-BFGS-B",
    lower = rep(0, K),
    control = list(maxit = control$max_iterations, factr = 1e2, pgtol = 1e-12)
  )
  ev <- dual_eval(opt$par)

  if (max(ev$grad) > control$tolerance * 0.5) {
    opt2 <- stats::optim(
      par = opt$par,
      fn = fn,
      gr = gr,
      method = "L-BFGS-B",
      lower = rep(0, K),
      control = list(maxit = control$max_iterations, factr = 1, pgtol = 1e-14)
    )
    ev2 <- dual_eval(opt2$par)
    if (ev2$dual >= ev$dual - 1e-12 * max(1, abs(ev$dual))) {
      opt <- opt2
      ev <- ev2
    }
  }

  x <- ev$x
  y <- ev$y
  primal_cost <- sum(c1 * x + c2 * y)
  dual_cost <- ev$dual * cost_scale
  dual_par <- opt$par
  lambda_unscaled <- opt$par * cost_scale / Vmax
  v_now <- .freetake_var(a, b, d, e, x, y)
  resid_rel <- v_now / Vmax - 1
  fixed <- function(m) {
    .freetake_fixed_solve(a, b, d, e, Vmax, c1, c2, N1, N, L_ult, m,
                          constraint_ids, control, accept)
  }
  adopt <- function(s) {
    x <<- s$x
    y <<- s$y
    primal_cost <<- s$cost
    resid_rel <<- .freetake_var(a, b, d, e, x, y) / Vmax - 1
    lambda_unscaled <<- s$lambda
    multiplier <<- s$multiplier
  }
  multiplier <- NULL
  if (any(resid_rel > accept)) {
    s <- fixed(pmin(pmax(y / x, m_min), m_max))
    if (!is.null(s)) adopt(s)
  }
  feasible <- all(resid_rel <= accept + 1e-7)
  gap_of <- function(cost) (cost - dual_cost) / cost
  if (!feasible || !(gap_of(primal_cost) <= control$kkt_tolerance)) {
    start <- if (feasible) pmin(pmax(y / x, m_min), m_max)
    ref <- .freetake_refine(fixed, start, m_min, m_max)
    if (!is.null(ref) && (!feasible || ref$cost < primal_cost)) {
      adopt(ref)
      feasible <- all(resid_rel <= accept + 1e-7)
    }
  }

  if (!feasible) {
    return(list(
      classification = "failed",
      converged = FALSE,
      message = sprintf("continuous free-take solver failed feasibility (max residual %.3e)",
                        max(resid_rel))
    ))
  }

  ## Any lambda >= 0 bounds the cost, so this restart can only tighten it.
  if (!(gap_of(primal_cost) <= control$kkt_tolerance) &&
      length(multiplier) == K && all(is.finite(multiplier))) {
    lam_fix <- pmax(0, multiplier * Vmax / cost_scale)
    opt3 <- stats::optim(
      par = lam_fix, fn = fn, gr = gr, method = "L-BFGS-B",
      lower = rep(0, K),
      control = list(maxit = control$max_iterations, factr = 1, pgtol = 1e-14)
    )
    for (cand in list(lam_fix, opt3$par)) {
      polished <- dual_eval(cand)$dual * cost_scale
      if (is.finite(polished) && polished > dual_cost) {
        dual_cost <- polished
        dual_par <- cand
      }
    }
  }

  rel_gap <- gap_of(primal_cost)
  certified <- is.finite(dual_cost) && rel_gap >= -1e-6 &&
    rel_gap <= control$kkt_tolerance
  list(
    classification = if (certified) "optimal" else "feasible",
    converged = certified,
    x = x,
    y = y,
    cost = primal_cost,
    lambda = lambda_unscaled,
    iterations = as.integer(opt$counts[["function"]]),
    convergence_code = as.integer(opt$convergence),
    primal_residual = max(0, max(resid_rel)),
    relative_duality_gap = rel_gap,
    dual_value = dual_cost,
    dual_lambda = dual_par * cost_scale / Vmax,
    message = if (certified) {
      "optimal free-take two-stage allocation"
    } else {
      sprintf(
        "feasible free-take two-stage allocation within %.3g%% of the dual bound",
        100 * rel_gap
      )
    },
    feasibility_tolerance = control$tolerance,
    kkt_tolerance = control$kkt_tolerance
  )
}

#' Solve the convex free-take problem at fixed per-stratum takes
#' @keywords internal
#' @noRd
.freetake_fixed_solve <- function(a, b, d, e, Vmax, c1, c2, N1, N, L_ult, m,
                                  constraint_ids, control, accept) {
  A_fix <- a + sweep(b, 1L, m, "/")
  bound_fix <- Vmax + colSums(sweep(d, 1L, m, "*") + e)
  c_fix <- c1 + c2 * m
  lo_fix <- pmax(1, L_ult / m)
  hi_fix <- pmin(N1, N / m)
  A_fix[A_fix < 0 & A_fix > -1e-12 * max(abs(A_fix))] <- 0
  if (any(!is.finite(A_fix)) || any(A_fix < 0) ||
      any(!is.finite(bound_fix)) || any(bound_fix <= 0) ||
      any(lo_fix > hi_fix)) {
    return(NULL)
  }
  sol <- .bethel_solve(
    A_fix, bound_fix, c_fix, lo_fix, hi_fix,
    control = control, constraint_ids = constraint_ids,
    accept = pmax(0, Vmax * accept / bound_fix)
  )
  ## An uncertified design still serves, checked against Vmax, not bound_fix.
  x <- sol$allocation
  if (is.null(x)) return(NULL)
  x <- pmin(x * max(1, colSums(A_fix / x) / bound_fix), hi_fix)
  resid <- .freetake_var(a, b, d, e, x, x * m) / Vmax - 1
  if (any(resid > accept + 1e-7)) return(NULL)
  list(m = m, x = x, y = x * m, cost = sum(c_fix * x), lambda = sol$lambda,
       multiplier = sol$cost_scale * sol$lambda_scaled / bound_fix)
}

#' Improve a free-take incumbent by coordinate search over the takes
#'
#' A feasible design, not a proof: only the dual bound certifies.
#' @keywords internal
#' @noRd
.freetake_refine <- function(fixed, m_start, m_min, m_max) {
  cost_of <- function(s) if (is.null(s)) Inf else s$cost
  best <- NULL
  for (m0 in list(m_start, m_max, sqrt(m_min * m_max))) {
    if (is.null(m0)) next
    s <- fixed(m0)
    if (cost_of(s) < cost_of(best)) best <- s
  }
  free <- which(m_max - m_min > 1e-12 * m_max)
  if (is.null(best) || !length(free)) return(best)
  for (pass in seq_len(4L)) {
    before <- best$cost
    for (h in free) {
      at <- function(t) {
        m <- best$m
        m[h] <- t
        fixed(m)
      }
      grid <- exp(seq(log(m_min[h]), log(m_max[h]), length.out = 9L))
      sols <- lapply(grid, at)
      vals <- vapply(sols, cost_of, 1)
      j <- which.min(vals)
      if (vals[j] < best$cost) best <- sols[[j]]
      o <- stats::optimize(
        function(lt) min(cost_of(at(exp(lt))), .Machine$double.xmax),
        log(grid[c(max(1L, j - 1L), min(9L, j + 1L))]),
        tol = 1e-10
      )
      s <- at(exp(o$minimum))
      if (cost_of(s) < best$cost) best <- s
    }
    if (before - best$cost <= 1e-10 * best$cost) break
  }
  best
}

#' Minimize a stratum's budget objective over the continuous take bounds
#'
#' The objective decreases in y when B or D is positive. On each piece of
#' the upper edge it is a constant plus a coefficient times 1/x, so the
#' endpoints and the bend at N / m_max suffice. When B and D are zero, use
#' the cheapest y at each x. Cost breaks equal objective values.
#' @keywords internal
#' @noRd
.freetake_stratum_obj_min <- function(c1, c2, A, B, D, N1, N, L_ult,
                                      m_min, m_max) {
  x_hi <- min(N1, N / m_min)
  x_lo <- min(x_hi, max(1, L_ult / m_max))
  if (B > 0 || D > 0) {
    x <- c(x_lo, x_hi, min(x_hi, max(x_lo, N / m_max)))
    y <- pmin(N, m_max * x)
  } else {
    x <- c(x_lo, x_hi)
    y <- pmin(N, pmax(L_ult, m_min * x))
  }
  value <- A / x + B / y - D * (y / x)
  best <- order(value, c1 * x + c2 * y)[1L]
  list(x = x[best], y = y[best])
}

#' Solve the fixed-budget continuous free-take 2-stage problem by epsilon-constraint
#' @keywords internal
#' @noRd
.freetake_budget_solve_continuous <- function(prob, budget, control = .bethel_control()) {
  H <- prob$H
  K <- prob$K
  obj <- prob$objective
  ft <- prob$ft_spec
  m_min <- ft$take_2_min
  m_max <- ft$take_2_max_cont

  Q_at <- function(x, y) {
    m <- y / x
    sum(obj$a_h / x + obj$b_h / y - obj$d_h * m - obj$e_h)
  }
  cost_at <- function(x, y) sum(prob$c1 * x + prob$c2 * y)

  if (K > 0L) {
    base <- .freetake_solve_core(
      prob$a, prob$b, prob$d, prob$e, prob$Vmax,
      prob$c1, prob$c2, prob$N_1, prob$N, prob$L_ult, m_min, m_max,
      prob$constraint_ids, control = control
    )
    if (identical(base$classification, "infeasible")) {
      return(list(status = "infeasible_targets", infeasible = base$infeasible))
    }
    if (!base$classification %in% c("optimal", "feasible")) {
      return(list(status = "failed", message = base$message))
    }
    base_x <- base$x
    base_y <- base$y
  } else {
    base_x <- numeric(H)
    base_y <- numeric(H)
    for (h in seq_len(H)) {
      pt <- .freetake_stratum_min(
        prob$c1[h], prob$c2[h], 0, 0, 0,
        prob$N_1[h], prob$N[h], prob$L_ult[h], m_min[h], m_max[h]
      )
      base_x[h] <- pt$x
      base_y[h] <- pt$y
    }
  }
  cost_min <- cost_at(base_x, base_y)
  budget_tol <- control$tolerance * abs(budget)
  b_eff <- budget + budget_tol
  if (budget < cost_min - budget_tol) {
    ## The solve accepts each target within acc of its variance, so a bound
    ## from the nominal dual gives back the multipliers' share of that slack.
    acc <- (1 + control$tolerance)^2 - 1
    exact <- K == 0L || is.null(base$dual_lambda)
    cost_bound <- if (exact) {
      cost_min
    } else if (length(base$dual_lambda) == K && all(is.finite(base$dual_lambda))) {
      base$dual_value - sum(pmax(base$dual_lambda, 0) * prob$Vmax * acc)
    } else {
      -Inf
    }
    if (budget < cost_bound - budget_tol) {
      return(list(
        status = "unaffordable",
        cost_min = cost_bound,
        shortfall = cost_bound - budget,
        exact = exact
      ))
    }
    relaxed <- if (!exact) {
      .freetake_solve_core(
        prob$a, prob$b, prob$d, prob$e, prob$Vmax * (1 + acc),
        prob$c1, prob$c2, prob$N_1, prob$N, prob$L_ult, m_min, m_max,
        prob$constraint_ids, control = control, accept = rep(0, K)
      )
    }
    fits <- !is.null(relaxed$x) &&
      cost_at(relaxed$x, relaxed$y) <= b_eff &&
      all(.freetake_var(prob$a, prob$b, prob$d, prob$e, relaxed$x, relaxed$y) <=
            prob$Vmax * (1 + acc))
    if (!fits) {
      return(list(status = "unknown", cost_min = cost_min,
                  cost_bound = cost_bound))
    }
    base_x <- relaxed$x
    base_y <- relaxed$y
    cost_min <- cost_at(base_x, base_y)
  }

  up_x <- numeric(H)
  up_y <- numeric(H)
  for (h in seq_len(H)) {
    pt <- .freetake_stratum_obj_min(
      prob$c1[h], prob$c2[h],
      obj$a_h[h], obj$b_h[h], obj$d_h[h],
      prob$N_1[h], prob$N[h], prob$L_ult[h], m_min[h], m_max[h]
    )
    up_x[h] <- pt$x
    up_y[h] <- pt$y
  }
  cost_max <- cost_at(up_x, up_y)
  q_hi <- Q_at(base_x, base_y)
  q_lo <- Q_at(up_x, up_y)
  up_meets <- K == 0L || all(
    .freetake_var(prob$a, prob$b, prob$d, prob$e, up_x, up_y) <=
      prob$Vmax * (1 + control$tolerance)^2
  )

  if (up_meets && budget >= cost_max - budget_tol) {
    return(list(
      status = "optimal",
      x = up_x,
      y = up_y,
      solved = NULL,
      lambda = rep(NA_real_, K),
      objective_multiplier = NA_real_,
      budget_sensitivity = NA_real_,
      objective_bound = q_lo,
      budget_binding = FALSE,
      cost = cost_max,
      root_iterations = 0L
    ))
  }

  a_ext <- cbind(prob$a, obj$a_h)
  b_ext <- cbind(prob$b, obj$b_h)
  d_ext <- cbind(prob$d, obj$d_h)
  e_ext <- cbind(prob$e, obj$e_h)
  ids_ext <- c(prob$constraint_ids, ".objective")
  accept_ext <- c(rep((1 + control$tolerance)^2 - 1, K), control$tolerance)

  solve_at <- function(qb) {
    .freetake_solve_core(
      a_ext, b_ext, d_ext, e_ext, c(prob$Vmax, qb),
      prob$c1, prob$c2, prob$N_1, prob$N, prob$L_ult, m_min, m_max,
      ids_ext, control = control, accept = accept_ext
    )
  }
  finish <- function(s, qb, lower = NULL, binding = TRUE, iterations = 0L) {
    dl <- s$dual_lambda
    mu <- if (length(dl) == K + 1L && all(is.finite(dl))) dl[K + 1L] else NA_real_
    value <- Q_at(s$x, s$y)
    if (is.null(lower) && isTRUE(mu > 0) && isTRUE(is.finite(s$dual_value))) {
      slack <- sum(dl[seq_len(K)] * prob$Vmax * accept_ext[seq_len(K)])
      lower <- qb + (s$dual_value - max(budget, s$cost) - slack) / mu
    }
    if (!is.null(lower) && is.finite(lower)) {
      lower <- if (lower > value + 1e-9 * abs(value)) NULL else min(lower, value)
    }
    gap <- if (is.null(lower) || !is.finite(lower)) {
      NA_real_
    } else if (value == 0) {
      0
    } else {
      (value - lower) / abs(value)
    }
    obj_mult <- s$lambda[K + 1L]
    list(
      status = if (is.finite(gap) && gap <= control$kkt_tolerance) "optimal"
        else "feasible",
      x = s$x,
      y = s$y,
      solved = s,
      lambda = s$lambda[seq_len(K)],
      objective_multiplier = obj_mult,
      budget_sensitivity = if (is.na(obj_mult) || obj_mult <= 0) NA_real_ else -1 / obj_mult,
      objective_bound = qb,
      objective_lower_bound = lower %||% NA_real_,
      objective_gap = gap,
      budget_binding = binding,
      cost = s$cost,
      root_iterations = iterations
    )
  }

  s_base <- list(x = base_x, y = base_y, cost = cost_min,
                 lambda = rep(NA_real_, K + 1L))
  best <- NULL
  fallback <- function(iterations = 0L) {
    if (is.null(best)) finish(s_base, q_hi, binding = FALSE) else
      finish(best$s, best$q, iterations = iterations)
  }
  span <- q_hi - q_lo
  if (!(span > 0)) return(fallback())

  bracket_lo <- q_lo + 1e-10 * span
  residual <- function(qb) {
    s <- solve_at(qb)
    if (!s$classification %in% c("optimal", "feasible")) return(NA_real_)
    if (s$cost <= b_eff && (is.null(best) || qb < best$q)) {
      best <<- list(q = qb, s = s)
    }
    s$cost - budget
  }
  f_lo <- residual(bracket_lo)
  s_lo <- best
  if (!is.na(f_lo) && f_lo <= 0) {
    return(finish(s_lo$s, bracket_lo, lower = if (up_meets) q_lo,
                  binding = FALSE))
  }
  f_hi <- residual(q_hi)
  if (is.na(f_lo) || is.na(f_hi)) return(fallback())
  if (f_hi >= 0) {
    if (f_hi > budget_tol) return(fallback())
    objective_bound <- q_hi
    root_iterations <- 0L
  } else {
    root <- stats::uniroot(
      residual, interval = c(bracket_lo, q_hi),
      f.lower = f_lo, f.upper = f_hi,
      tol = span * 1e-12, maxiter = 1000L
    )
    objective_bound <- root$root
    root_iterations <- root$iter
  }
  solved <- solve_at(objective_bound)
  if (!solved$classification %in% c("optimal", "feasible") ||
        solved$cost > b_eff) {
    return(fallback(root_iterations))
  }
  finish(solved, objective_bound, iterations = root_iterations)
}

#' Build branch table for exact discrete `(n_1, take_2)` global search
#'
#' Each branch is one `(stratum h, integer take m)` pair with valid PSU bounds
#' `[lo_i, hi_i]`, positive field cost per PSU `c1_h + c2_h * m`, and
#' non-negative variance load terms (`g0 >= 0`, `g2 >= 0`) obtained by shifting
#' each stratum's constant term by its stratum minimum.
#' @keywords internal
#' @noRd
.freetake_branch_table <- function(prob, tolerance) {
  H <- prob$H
  K <- prob$K
  cands <- prob$ft_spec$cand_takes

  owner_list <- vector("list", H)
  take_list <- vector("list", H)
  lo_list <- vector("list", H)
  hi_list <- vector("list", H)
  for (h in seq_len(H)) {
    ms <- as.numeric(cands[[h]])
    lo <- pmax(1, ceiling(prob$L_ult[h] / ms - 1e-9))
    hi <- pmin(prob$N_1[h], floor(prob$N[h] / ms + 1e-9))
    ok <- lo <= hi
    if (!any(ok)) {
      stop(
        sprintf("no integer (n_1, take_2) satisfies the bounds for stratum '%s'",
                prob$stratum_ids[h]),
        call. = FALSE
      )
    }
    owner_list[[h]] <- rep.int(h, sum(ok))
    take_list[[h]] <- ms[ok]
    lo_list[[h]] <- lo[ok]
    hi_list[[h]] <- hi[ok]
  }
  owner <- unlist(owner_list, use.names = FALSE)
  m_vec <- unlist(take_list, use.names = FALSE)
  lo_vec <- unlist(lo_list, use.names = FALSE)
  hi_vec <- unlist(hi_list, use.names = FALSE)
  B <- length(owner)
  cost_vec <- prob$c1[owner] + prob$c2[owner] * m_vec

  g0 <- matrix(0, B, K)
  g2 <- matrix(0, B, K)
  denom <- numeric(K)
  if (K > 0L) {
    accept_mult <- (1 + tolerance)^2
    for (k in seq_len(K)) {
      A_b <- prob$a[owner, k] + prob$b[owner, k] / m_vec
      B_b <- -prob$d[owner, k] * m_vec - prob$e[owner, k]
      B_min_h <- vapply(seq_len(H), function(h) min(B_b[owner == h]), numeric(1L))
      denom_k <- prob$Vmax[k] * accept_mult - sum(B_min_h)
      denom[k] <- denom_k
      ## As for a fixed take: the search's own 1 + 1e-9 must not stack on
      ## the slack already in denom_k.
      g0[, k] <- (B_b - B_min_h[owner]) / denom_k * .whole_target_scale()
      g2[, k] <- pmax(0, A_b) / denom_k * .whole_target_scale()
    }
  }

  obj_o0 <- NULL
  obj_o2 <- NULL
  obj_offset <- 0
  if (!is.null(prob$objective)) {
    obj <- prob$objective
    A_obj <- obj$a_h[owner] + obj$b_h[owner] / m_vec
    B_obj <- -obj$d_h[owner] * m_vec - obj$e_h[owner]
    B_min_h <- vapply(seq_len(H), function(h) min(B_obj[owner == h]), numeric(1L))
    obj_o0 <- B_obj - B_min_h[owner]
    obj_o2 <- pmax(0, A_obj)
    obj_offset <- sum(B_min_h)
  }

  list(
    owner = owner,
    m = m_vec,
    lo = lo_vec,
    hi = hi_vec,
    cost = cost_vec,
    g0 = g0,
    g2 = g2,
    denom = denom,
    obj_o0 = obj_o0,
    obj_o2 = obj_o2,
    obj_offset = obj_offset
  )
}

#' Exact integer `(n_1, take_2)` global search for minimum-cost mode
#' @keywords internal
#' @noRd
.freetake_whole_targets <- function(prob, control, max_nodes, gap,
                                    lambda_cont = NULL) {
  bt <- .freetake_branch_table(prob, control$tolerance)
  B <- length(bt$owner)
  K <- prob$K
  q <- .whole_problem(
    o1 = bt$cost, o2 = 0,
    g1 = matrix(0, B, K), g2 = bt$g2,
    lo = bt$lo, hi = bt$hi,
    owner = bt$owner, o0 = 0, g0 = bt$g0
  )
  lam_hint <- if (!is.null(lambda_cont) && length(lambda_cont) == K) {
    lambda_cont * bt$denom
  }
  grain <- .menu_grain(c(prob$c1, prob$c2))
  r <- .whole_search(q, grain = grain, max_nodes = max_nodes, gap = gap,
                     lambda = lam_hint)
  if (identical(r$status, "infeasible")) {
    bad <- prob$constraint_ids[
      .whole_load(q, .whole_plan(q, q$hi)) > 1 + 1e-9
    ]
    if (!length(bad)) bad <- prob$constraint_ids
    stop(sprintf("no feasible integer allocation meets constraints: %s",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  if (is.null(r$x)) {
    stop(structure(
      list(message = sprintf(
        "no whole allocation meeting every target was found in %d search node%s, so whether one exists is unknown. Raise 'control$max_nodes'",
        as.integer(max_nodes), if (max_nodes == 1) "" else "s"
      ), call = NULL),
      class = c("svyplan_plan_unknown", "error", "condition")
    ))
  }
  list(
    n_1_int = as.numeric(r$x),
    take_2_int = as.numeric(bt$m[r$branch]),
    branch = r$branch,
    r = r,
    bt = bt
  )
}

#' Exact integer `(n_1, take_2)` global search for fixed-budget mode
#' @keywords internal
#' @noRd
.freetake_whole_budget <- function(prob, budget, control, max_nodes, gap) {
  bt <- .freetake_branch_table(prob, control$tolerance)
  H <- prob$H
  B <- length(bt$owner)
  K <- prob$K
  budget_tol <- control$tolerance * abs(budget)

  cheap_b <- integer(H)
  cheap_x <- numeric(H)
  for (h in seq_len(H)) {
    idx <- which(bt$owner == h)
    c_lo <- bt$cost[idx] * bt$lo[idx]
    j <- idx[which.min(c_lo)]
    cheap_b[h] <- j
    cheap_x[h] <- bt$lo[j]
  }
  lowest <- sum(bt$cost[cheap_b] * cheap_x)
  if (lowest > budget + budget_tol) {
    stop(sprintf(
      "no integer allocation fits the budget: the smallest allocation satisfying the stratum bounds costs %.6g against a budget of %.6g",
      lowest, budget
    ), call. = FALSE)
  }

  start_x <- cheap_x
  start_b <- cheap_b
  if (K > 0L) {
    cheapest <- .freetake_whole_targets(prob, control, max_nodes, gap)
    cost_cheap <- sum((prob$c1 + prob$c2 * cheapest$take_2_int) * cheapest$n_1_int)
    if (cost_cheap > budget + budget_tol) {
      if (cheapest$r$bound > budget + budget_tol) {
        stop(sprintf(
          "no integer allocation meets the precision targets within the budget: the cheapest target-feasible integer allocation costs %.6g against a budget of %.6g",
          cost_cheap, budget
        ), call. = FALSE)
      }
      stop(structure(
        list(message = sprintf(
          "no whole allocation meeting the precision targets within 'budget' was found in %d search node%s, so whether one exists is unknown. The cheapest found costs %.6g and none costs less than %.6g. Raise 'control$max_nodes'",
          as.integer(max_nodes), if (max_nodes == 1) "" else "s",
          cost_cheap, cheapest$r$bound
        ), call = NULL),
        class = c("svyplan_plan_unknown", "error", "condition")
      ))
    }
    start_x <- cheapest$n_1_int
    start_b <- cheapest$branch
  }

  q <- .whole_problem(
    o1 = 0, o2 = bt$obj_o2,
    g1 = cbind(.whole_budget_load(bt$cost, budget + budget_tol),
               matrix(0, B, K)),
    g2 = cbind(0, bt$g2),
    lo = bt$lo, hi = bt$hi,
    owner = bt$owner,
    o0 = bt$obj_o0,
    g0 = cbind(0, bt$g0)
  )
  start_plan <- .whole_plan(q, start_x, branch = start_b)
  r <- .whole_search(q, start = start_plan, max_nodes = max_nodes, gap = gap)
  if (is.null(r$x)) {
    stop("internal error: the budget search lost the plan it started from",
         call. = FALSE)
  }
  list(
    n_1_int = as.numeric(r$x),
    take_2_int = as.numeric(bt$m[r$branch]),
    r = r,
    bt = bt
  )
}

#' Evaluate constraints & objective on a 2-stage `(n_1, take_2)` using exact fixed-take path
#' @keywords internal
#' @noRd
.freetake_assess <- function(
  frame, measures, targets, objective,
  n_1, take_2, alpha, deff, resp_rate, eligible_per_unit, eligibility_rate,
  df, min_n_stratum, fpc, lambda = NULL, tolerance = .bethel_control()$tolerance
) {
  H <- nrow(frame)
  prob_eval <- .freetake_build_problem(
    frame = frame, measures = measures, targets = targets,
    ft_spec = list(
      active = TRUE, psu = FALSE, free = rep(FALSE, H),
      take_fixed = take_2, take_2_min = take_2,
      take_2_max_cont = take_2, take_2_max_int = floor(take_2),
      take_2_step = rep(1, H), cand_takes = as.list(round(take_2))
    ),
    deff = deff, resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit, eligibility_rate = eligibility_rate,
    alpha = alpha, df = df, min_n_stratum = min_n_stratum,
    objective = objective, fpc = fpc
  )
  y <- n_1 * take_2
  K <- prob_eval$K
  variance <- if (K > 0L) {
    .freetake_var(prob_eval$a, prob_eval$b, prob_eval$d, prob_eval$e, n_1, y)
  } else numeric(0)
  if (any(variance < 0)) {
    stop("internal error: allocation produced a materially negative variance",
         call. = FALSE)
  }
  se_total <- sqrt(variance)
  se <- se_total / prob_eval$domain_N
  q <- prob_eval$q_target
  moe <- q * se
  cv <- se_total / abs(prob_eval$total)
  metric <- prob_eval$constraint_meta$.metric
  rmoe <- .rmoe_from_moe(moe, prob_eval$total / prob_eval$domain_N)
  achieved <- as.numeric(ifelse(metric == "cv", cv, ifelse(metric == "moe", moe, rmoe)))
  target <- prob_eval$constraint_meta$.target
  ratio <- achieved / target
  out <- prob_eval$constraint_meta
  out$.achieved <- achieved
  out$.ratio <- ratio
  out$.residual <- ratio - 1
  out$.tolerance <- rep(tolerance, nrow(out))
  out$.se <- se
  out$.cv <- cv
  out$.moe <- moe
  out$.rmoe <- rmoe
  out$.pass <- ratio <= 1 + tolerance
  out$.binding <- abs(ratio - 1) <= max(1e-6, 10 * tolerance)
  if (is.null(lambda)) {
    out$.multiplier <- rep(NA_real_, nrow(out))
    out$.sensitivity <- rep(NA_real_, nrow(out))
  } else {
    out$.multiplier <- lambda
    sensitivity <- numeric(nrow(out))
    cv_idx <- metric == "cv"
    sensitivity[cv_idx] <- -2 * lambda[cv_idx] * target[cv_idx] *
      prob_eval$total[cv_idx]^2
    moe_idx <- metric == "moe"
    sensitivity[moe_idx] <- -2 * lambda[moe_idx] *
      prob_eval$domain_N[moe_idx]^2 * target[moe_idx] / q[moe_idx]^2
    rmoe_idx <- metric == "rmoe"
    sensitivity[rmoe_idx] <- -2 * lambda[rmoe_idx] * target[rmoe_idx] *
      prob_eval$total[rmoe_idx]^2 / q[rmoe_idx]^2
    out$.sensitivity <- sensitivity
  }

  obj_res <- NULL
  if (!is.null(prob_eval$objective)) {
    obj <- prob_eval$objective
    v_obj <- .freetake_var(obj$a_mat, obj$b_mat, obj$d_mat, obj$e_mat, n_1, y)
    if (any(v_obj < 0)) {
      stop("internal error: allocation produced a materially negative objective variance",
           call. = FALSE)
    }
    relvar <- v_obj / obj$total^2
    contribution <- obj$priority * relvar
    val <- sum(contribution)
    obj_res <- list(
      value = val,
      components = data.frame(
        component = obj$component,
        name = obj$name,
        domain = obj$domain,
        level = obj$level,
        priority = obj$priority,
        .relvar = relvar,
        .cv = sqrt(relvar),
        .contribution = contribution,
        .share = if (val > 0) contribution / val else rep(NA_real_, length(relvar)),
        stringsAsFactors = FALSE
      )
    )
  }

  list(
    constraints = out,
    cost = sum(prob_eval$c1 * n_1 + prob_eval$c2 * y),
    all_pass = all(out$.pass),
    objective = obj_res
  )
}

#' Public-wrapper engine for free-take 2-stage joint constrained allocation
#' @keywords internal
#' @noRd
.n_alloc_bethel_freetake <- function(
  frame,
  measures,
  targets,
  ft_spec,
  alpha,
  deff,
  resp_rate,
  min_n_stratum,
  objective = NULL,
  budget = NULL,
  df = NULL,
  fpc = "unit",
  eligible_per_unit = 1,
  eligibility_rate = 1,
  take_2_min = NULL,
  take_2_max = NULL,
  take_2_step = NULL,
  max_nodes = 5000,
  gap = 0
) {
  control <- .bethel_control()
  prob <- .freetake_build_problem(
    frame = frame,
    measures = measures,
    targets = targets,
    ft_spec = ft_spec,
    deff = deff,
    resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit,
    eligibility_rate = eligibility_rate,
    alpha = alpha,
    df = df,
    min_n_stratum = min_n_stratum,
    objective = objective,
    fpc = fpc
  )

  mode <- if (is.null(objective)) "targets" else "budget_objective"
  if (mode == "targets") {
    solved <- .freetake_solve_core(
      prob$a, prob$b, prob$d, prob$e, prob$Vmax,
      prob$c1, prob$c2, prob$N_1, prob$N, prob$L_ult,
      ft_spec$take_2_min, ft_spec$take_2_max_cont,
      prob$constraint_ids, control = control
    )
    if (identical(solved$classification, "infeasible")) {
      stop(sprintf("precision targets are infeasible: %s",
                   paste(solved$infeasible, collapse = ", ")), call. = FALSE)
    }
    if (!solved$classification %in% c("optimal", "feasible")) {
      stop(sprintf("joint allocation solver failed certification: %s",
                   solved$message), call. = FALSE)
    }
    x_cont <- solved$x
    y_cont <- solved$y
    take_cont <- y_cont / x_cont
    lambda <- solved$lambda
    optimization <- solved
    optimization$x <- NULL
    optimization$y <- NULL
  } else {
    fit <- .freetake_budget_solve_continuous(prob, budget, control)
    if (identical(fit$status, "infeasible_targets")) {
      stop(sprintf("precision targets are infeasible: %s",
                   paste(fit$infeasible, collapse = ", ")), call. = FALSE)
    }
    if (identical(fit$status, "unaffordable")) {
      stop(sprintf(
        "'budget' cannot fund the precision targets: the cheapest target-feasible design costs %s%.6g, short by %.6g",
        if (isTRUE(fit$exact)) "" else "at least ",
        fit$cost_min, fit$shortfall
      ), call. = FALSE)
    }
    if (identical(fit$status, "unknown")) {
      stop(structure(
        list(message = sprintf(
          "no continuous design meeting the precision targets within 'budget' was found, so whether one exists is unknown. The cheapest found costs %.6g and none costs less than %.6g",
          fit$cost_min, fit$cost_bound
        ), call = NULL),
        class = c("svyplan_plan_unknown", "error", "condition")
      ))
    }
    if (!fit$status %in% c("optimal", "feasible")) {
      stop(sprintf("joint allocation solver failed certification: %s",
                   fit$message %||% "unknown"), call. = FALSE)
    }
    x_cont <- fit$x
    y_cont <- fit$y
    take_cont <- y_cont / x_cont
    lambda <- fit$lambda
    optimization <- if (is.null(fit$solved)) {
      list(
        classification = fit$status,
        converged = identical(fit$status, "optimal"),
        message = "the budget exceeds the cost of the upper allocation bounds",
        feasibility_tolerance = control$tolerance,
        kkt_tolerance = control$kkt_tolerance
      )
    } else {
      opt <- fit$solved
      opt$x <- NULL
      opt$y <- NULL
      opt$classification <- fit$status
      opt$converged <- identical(fit$status, "optimal")
      opt
    }
    optimization$budget <- budget
    optimization$budget_binding <- fit$budget_binding
    optimization$budget_residual <- budget - fit$cost
    optimization$budget_sensitivity <- fit$budget_sensitivity
    optimization$objective_bound <- fit$objective_bound
    optimization$objective_lower_bound <- fit$objective_lower_bound
    optimization$objective_gap <- fit$objective_gap
    optimization$objective_multiplier <- fit$objective_multiplier
    optimization$root_iterations <- fit$root_iterations
  }

  cont_eval <- .freetake_assess(
    frame, measures, targets, objective,
    n_1 = x_cont, take_2 = take_cont,
    alpha = alpha, deff = deff, resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit, eligibility_rate = eligibility_rate,
    df = df, min_n_stratum = min_n_stratum, fpc = fpc,
    lambda = lambda, tolerance = control$tolerance
  )

  if (mode == "targets") {
    disc <- .freetake_whole_targets(prob, control, max_nodes, gap,
                                    lambda_cont = lambda)
    search_rec <- .whole_record(disc$r, "whole units (free take)", max_nodes)
  } else {
    disc <- .freetake_whole_budget(prob, budget, control, max_nodes, gap)
    search_rec <- .whole_record(disc$r, "whole units (free take)", max_nodes,
                                offset = disc$bt$obj_offset)
  }
  optimization$search <- search_rec

  n_1_int <- disc$n_1_int
  take_2_int <- disc$take_2_int
  n_int <- n_1_int * take_2_int
  int_eval <- .freetake_assess(
    frame, measures, targets, objective,
    n_1 = n_1_int, take_2 = take_2_int,
    alpha = alpha, deff = deff, resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit, eligibility_rate = eligibility_rate,
    df = df, min_n_stratum = min_n_stratum, fpc = fpc,
    lambda = NULL, tolerance = control$tolerance
  )
  if (!int_eval$all_pass) {
    stop("internal error: free-take integer allocation failed final validation",
         call. = FALSE)
  }
  if (mode == "budget_objective") {
    budget_tol <- control$tolerance * abs(budget)
    if (int_eval$cost > budget + budget_tol) {
      stop("internal error: free-take integer allocation exceeded the budget",
           call. = FALSE)
    }
  }

  lower_n1 <- pmax(1, prob$L_ult / take_cont)
  upper_n1 <- pmin(prob$N_1, prob$N / take_cont)
  detail <- data.frame(
    stratum = prob$stratum_ids,
    N = prob$N,
    N_1 = prob$N_1,
    cost_1 = prob$c1,
    cost_2 = prob$c2,
    take_2 = take_cont,
    n_1 = x_cont,
    n_1_int = n_1_int,
    take_2_int = take_2_int,
    n = y_cont,
    n_int = n_int,
    weight = prob$N / y_cont,
    .lower = pmax(prob$L_ult, ft_spec$take_2_min),
    .upper = pmin(prob$N, prob$N_1 * ft_spec$take_2_max_cont),
    .lower_n_1 = lower_n1,
    .upper_n_1 = upper_n1,
    .binding = abs(x_cont - lower_n1) <= 1e-6 | abs(x_cont - upper_n1) <= 1e-6 |
      (ft_spec$free & (abs(take_cont - ft_spec$take_2_min) <= 1e-6 |
                       abs(take_cont - ft_spec$take_2_max_cont) <= 1e-6)),
    stringsAsFactors = FALSE
  )

  operational <- list(
    n = sum(n_int),
    cost = int_eval$cost,
    constraints = int_eval$constraints,
    all_pass = int_eval$all_pass,
    status = disc$r$status
  )
  if (mode == "budget_objective") {
    operational$objective <- int_eval$objective$components
    operational$objective_value <- int_eval$objective$value
    operational$budget_residual <- budget - int_eval$cost
  }

  summary_problem <- list(
    stratum_ids = prob$stratum_ids,
    constraint_ids = prob$constraint_ids,
    population_N = prob$N,
    stages = 2L,
    variance_model = "two_stage_free_take_wald",
    cost = prob$c1 + prob$c2 * take_2_int,
    lower = pmax(prob$L_ult, ft_spec$take_2_min) / take_2_int,
    upper = pmin(prob$N, prob$N_1 * ft_spec$take_2_max_cont) / take_2_int,
    df = prob$df,
    membership = prob$membership,
    mean_hk = prob$mean_hk,
    var_hk = prob$var_hk,
    deff_hk = prob$deff_hk,
    resp_hk = prob$resp_hk,
    elig_hk = prob$elig_hk,
    scope_hk = prob$scope_hk,
    stage = list(
      stages = 2L,
      take = take_2_int,
      takes = matrix(take_2_int, ncol = 1L, dimnames = list(NULL, "take_2")),
      N_stages = matrix(prob$N_1, ncol = 1L, dimnames = list(NULL, "N_1")),
      costs = cbind(cost_1 = prob$c1, cost_2 = prob$c2)
    )
  )

  params <- list(
    frame = frame,
    fpc = fpc,
    measures = prob$measures,
    targets = prob$targets,
    alpha = alpha,
    df = df,
    deff = deff,
    resp_rate = resp_rate,
    unit_cost = NULL,
    cost_h = prob$c1 + prob$c2 * take_cont,
    min_n_stratum = min_n_stratum,
    take_2_min = take_2_min,
    take_2_max = take_2_max,
    take_2_step = take_2_step,
    free_take = TRUE,
    mode = mode,
    alloc = "bethel",
    n_h = y_cont,
    n_1_h = x_cont,
    take_2_h = take_cont,
    take_2_int = take_2_int,
    stages = 2L,
    problem = summary_problem,
    achieved = list(n = sum(y_cont), cost = cont_eval$cost)
  )
  if (eligible_per_unit != 1) params$eligible_per_unit <- eligible_per_unit
  if (!identical(eligibility_rate, 1)) params$eligibility_rate <- eligibility_rate
  if (mode == "budget_objective") {
    params$objective <- prob$objective$spec
    params$budget <- budget
    params$achieved$objective_value <- cont_eval$objective$value
  }
  params$feasibility_tolerance <- control$tolerance
  carried <- .alloc_carry_labels(detail, frame)
  detail <- carried$detail
  params$strata_keys <- carried$keys

  .new_svyplan_alloc_joint(
    n = sum(y_cont),
    stages = 2L,
    params = params,
    targets = prob$targets,
    detail = detail,
    rates = prob$rates,
    binding = prob$constraint_ids[cont_eval$constraints$.binding],
    operational = operational,
    constraints = cont_eval$constraints,
    optimization = optimization,
    objective = cont_eval$objective$components,
    objective_value = cont_eval$objective$value
  )
}

#' Precision assessment for a continuous free-take 2-stage joint fit
#' @keywords internal
#' @noRd
.prec_alloc_freetake <- function(
  frame,
  n,
  take_2,
  measures,
  targets,
  alpha,
  deff,
  resp_rate,
  min_n_stratum = NULL,
  objective = NULL,
  budget = NULL,
  df = NULL,
  fpc = "unit",
  eligible_per_unit = 1,
  eligibility_rate = 1,
  take_2_min = NULL,
  take_2_max = NULL,
  take_2_step = NULL
) {
  n_1 <- as.numeric(n) / as.numeric(take_2)
  ev <- .freetake_assess(
    frame = frame, measures = measures, targets = targets,
    objective = objective, n_1 = n_1, take_2 = as.numeric(take_2),
    alpha = alpha, deff = deff, resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit, eligibility_rate = eligibility_rate,
    df = df, min_n_stratum = min_n_stratum, fpc = fpc,
    lambda = NULL, tolerance = .bethel_control()$tolerance
  )
  stratum <- as.character(frame$stratum)
  N <- as.numeric(frame$N)
  N_1 <- as.numeric(round(frame$N_1))
  max_weight <- if ("max_weight" %in% names(frame)) frame$max_weight else rep(NA_real_, nrow(frame))
  ub <- .alloc_bounds(N, max_weight, rep(FALSE, nrow(frame)), min_n_stratum)
  lo <- pmax(ub$m_h, as.numeric(take_2))
  hi <- pmin(N, N_1 * as.numeric(take_2))
  tol <- .bethel_control()$tolerance
  lo_viol <- as.numeric(n) < lo - tol
  hi_viol <- as.numeric(n) > hi + tol
  bounds <- data.frame(
    stratum = stratum,
    n = as.numeric(n),
    .lower = lo,
    .upper = hi,
    .lower_violation = lo_viol,
    .upper_violation = hi_viol,
    .pass = !(lo_viol | hi_viol),
    stringsAsFactors = FALSE
  )
  params <- list(
    frame = frame,
    fpc = fpc,
    measures = measures,
    targets = targets,
    n = as.numeric(n),
    alpha = alpha,
    df = df,
    deff = deff,
    resp_rate = resp_rate,
    unit_cost = NULL,
    min_n_stratum = min_n_stratum,
    take_2_min = take_2_min,
    take_2_max = take_2_max,
    take_2_step = take_2_step,
    take_2_h = as.numeric(take_2),
    free_take = TRUE,
    stages = 2L,
    achieved = list(n = sum(as.numeric(n)), cost = ev$cost),
    feasibility_tolerance = tol
  )
  if (eligible_per_unit != 1) params$eligible_per_unit <- eligible_per_unit
  if (!identical(eligibility_rate, 1)) params$eligibility_rate <- eligibility_rate
  if (!is.null(ev$objective)) {
    params$objective <- objective
    params$budget <- budget
    params$achieved$objective_value <- ev$objective$value
    params$budget_residual <- if (is.null(budget)) NA_real_ else budget - ev$cost
  }
  .new_svyplan_prec(
    se = ev$constraints$.se,
    moe = ev$constraints$.moe,
    cv = ev$constraints$.cv,
    type = "alloc_joint",
    method = "bethel",
    params = params,
    detail = ev$constraints,
    bounds = bounds,
    objective = ev$objective$components,
    objective_value = ev$objective$value
  )
}
