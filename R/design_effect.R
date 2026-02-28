#' Design Effect
#'
#' Estimate the design effect using various methods. This is an S3 generic
#' that dispatches on the class of `x`.
#'
#' @param x A numeric vector of survey weights (for diagnostic methods),
#'   or `NULL` (for the `"cluster"` planning method).
#' @param ... Additional arguments passed to methods.
#'
#' @return For `"kish"`, `"cluster"`, `"henry"`, `"spencer"`: a numeric
#'   scalar. For `"cr"`: a list with `$strata` (data frame) and
#'   `$overall` (numeric scalar).
#'
#' @details
#' The `"kish"` method uses only weights and produces a single survey-wide
#' DEFF. The `"henry"`, `"spencer"`, and `"cr"` methods are
#' **outcome-dependent**: they require `y`, and the resulting DEFF varies
#' by outcome variable.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' Henry, K. A. and Valliant, R. (2015). A design effect measure
#' for calibration weighting in single-stage samples. *Survey
#' Methodology*, 41(2), 315--331.
#'
#' Spencer, B. D. (2000). An approximate design effect for unequal
#' weighting when measurements may correlate with selection
#' probabilities. *Survey Methodology*, 26(2), 137--138.
#'
#' Chen, S. and Rust, K. (2017). An extension of Kish's formula
#' for design effects to two- and three-stage designs with
#' stratification. *Journal of Survey Statistics and Methodology*,
#' 5(2), 111--130.
#'
#' @seealso [effective_n()] for effective sample size, [varcomp()] for
#'   estimating inputs to the `"cluster"` method.
#'
#' @examples
#' # Kish design effect from weights
#' set.seed(2)
#' w <- runif(100, 1, 5)
#' design_effect(w, method = "kish")
#'
#' # Planning: cluster design effect
#' design_effect(delta = 0.05, psu_size = 25, method = "cluster")
#'
#' @name design_effect
#' @export
design_effect <- function(x = NULL, ...) {
  if (is.null(x)) {
    design_effect.default(x, ...)
  } else {
    UseMethod("design_effect")
  }
}

#' @describeIn design_effect Method for numeric weights vector.
#'
#' @param y Outcome variable (required for `"henry"`, `"spencer"`, `"cr"`).
#' @param x_cal Calibration covariate (required for `"henry"`).
#' @param prob 1-draw selection probabilities (required for `"spencer"`).
#' @param strata_id Stratum IDs (required for `"cr"`).
#' @param cluster_id Cluster IDs (required for `"cr"`).
#' @param stages Integer vector of stages per stratum (required for `"cr"`).
#' @param method For numeric weights: one of `"kish"` (default), `"henry"`,
#'   `"spencer"`, or `"cr"`. For planning (no weights): `"cluster"`
#'   (default and only option).
#'
#' @export
design_effect.numeric <- function(
  x,
  ...,
  y = NULL,
  x_cal = NULL,
  prob = NULL,
  strata_id = NULL,
  cluster_id = NULL,
  stages = NULL,
  method = "kish"
) {
  method <- match.arg(method, c("kish", "henry", "spencer", "cr"))
  check_weights(x)

  switch(
    method,
    kish = .deff_kish(x),
    henry = .deff_henry(x, y, x_cal),
    spencer = .deff_spencer(x, y, prob),
    cr = .deff_cr(x, y, strata_id, cluster_id, stages)
  )
}

#' @describeIn design_effect Planning method (no weights needed).
#'
#' @param delta ICC / homogeneity measure, scalar or `svyplan_varcomp`
#'   (extracts `delta[1]`).
#' @param psu_size Mean cluster size (scalar).
#'
#' @export
design_effect.default <- function(
  x = NULL,
  ...,
  delta = NULL,
  psu_size = NULL,
  method = "cluster"
) {
  method <- match.arg(method, "cluster")

  if (is.null(delta) || is.null(psu_size)) {
    stop(
      "'delta' and 'psu_size' are required for the cluster method",
      call. = FALSE
    )
  }

  if (inherits(delta, "svyplan_varcomp")) {
    delta <- delta$delta[1L]
  }

  check_scalar(psu_size, "psu_size")
  check_delta(delta, expected_length = 1L)

  1 + (psu_size - 1) * delta
}

#' Kish (1965) design effect: 1 + CV^2(w)
#' @keywords internal
#' @noRd
.deff_kish <- function(w) {
  n <- length(w)
  1 + sum((w - mean(w))^2) / n / mean(w)^2
}

#' Henry and Valliant (2015) design effect
#' @keywords internal
#' @noRd
.deff_henry <- function(w, y, x_cal) {
  if (is.null(y) || is.null(x_cal)) {
    stop("'y' and 'x_cal' are required for the Henry method", call. = FALSE)
  }
  n <- length(w)
  check_covariate(y, n, "y")
  check_covariate(x_cal, n, "x_cal")
  dK <- .deff_kish(w)

  sw <- sqrt(w)
  fit <- lm.fit(sw * cbind(1, x_cal), sw * y)
  e <- fit$residuals / sw
  A <- fit$coefficients[1L]
  u <- A + e

  Nhat <- sum(w)
  wbar <- weighted.mean(w, w)
  vy <- .wtdvar(y, w)
  if (vy < .Machine$double.eps) {
    return(1.0)
  }
  vw <- .wtdvar(w, w)
  if (vw < .Machine$double.eps) {
    return(1.0)
  }
  vu <- .wtdvar(u, w)
  vu2 <- .wtdvar(u^2, w)
  ubar <- weighted.mean(u, w)
  u2bar <- weighted.mean(u^2, w)

  if (vu < .Machine$double.eps) {
    rho_u2w <- 0
    rho_uw <- 0
  } else {
    rho_u2w <- sum(w * (u^2 - u2bar) * (w - wbar)) / Nhat / sqrt(vu2 * vw)
    rho_uw <- sum(w * (u - ubar) * (w - wbar)) / Nhat / sqrt(vu * vw)
  }

  as.numeric(
    dK *
      vu /
      vy +
      n *
        sqrt(vw) *
        (rho_u2w * sqrt(vu2) - 2 * A * rho_uw * sqrt(vu)) /
        (Nhat * vy)
  )
}

#' Spencer (2000) design effect
#' @keywords internal
#' @noRd
.deff_spencer <- function(w, y, prob) {
  if (is.null(y) || is.null(prob)) {
    stop("'y' and 'prob' are required for the Spencer method", call. = FALSE)
  }
  n <- length(w)
  check_covariate(y, n, "y")
  check_covariate(prob, n, "prob")
  if (any(prob <= 0) || any(prob > 1)) {
    stop("'prob' values must be in (0, 1]", call. = FALSE)
  }
  dK <- .deff_kish(w)

  sw <- sqrt(w)
  fit <- lm.fit(sw * cbind(1, prob), sw * y)
  e <- fit$residuals / sw
  A <- fit$coefficients[1L]

  Nhat <- sum(w)
  ybar <- weighted.mean(y, w)
  pbar <- weighted.mean(prob, w)
  wbar <- weighted.mean(w, w)
  vy <- .wtdvar(y, w)
  if (vy < .Machine$double.eps) {
    return(1.0)
  }
  vw <- .wtdvar(w, w)
  if (vw < .Machine$double.eps) {
    return(1.0)
  }
  vp <- .wtdvar(prob, w)
  ve <- .wtdvar(e, w)
  ve2 <- .wtdvar(e^2, w)
  ebar <- weighted.mean(e, w)
  e2bar <- weighted.mean(e^2, w)

  rho_yp <- if (vp < .Machine$double.eps) {
    0
  } else {
    sum(w * (y - ybar) * (prob - pbar)) / Nhat / sqrt(vy * vp)
  }

  if (ve < .Machine$double.eps) {
    rho_e2w <- 0
    rho_ew <- 0
  } else {
    rho_e2w <- sum(w * (e^2 - e2bar) * (w - wbar)) / Nhat / sqrt(ve2 * vw)
    rho_ew <- sum(w * (e - ebar) * (w - wbar)) / Nhat / sqrt(ve * vw)
  }

  as.numeric(
    A^2 *
      (dK - 1) /
      vy +
      dK * (1 - rho_yp^2) +
      n * rho_e2w * sqrt(ve2 * vw) / (Nhat * vy) +
      2 * n * A * rho_ew * sqrt(ve * vw) / (Nhat * vy)
  )
}

#' Chen and Rust (2017) design effect
#' @keywords internal
#' @noRd
.deff_cr <- function(w, y, strata_id, cluster_id, stages) {
  if (is.null(y)) {
    stop("'y' is required for the CR method", call. = FALSE)
  }
  if (is.null(strata_id) && is.null(cluster_id)) {
    stop(
      "at least one of 'strata_id' or 'cluster_id' is required for the CR method",
      call. = FALSE
    )
  }

  n <- length(w)
  check_covariate(y, n, "y")
  if (!is.null(strata_id) && length(strata_id) != n) {
    stop("'strata_id' must have the same length as weights", call. = FALSE)
  }
  if (!is.null(cluster_id) && length(cluster_id) != n) {
    stop("'cluster_id' must have the same length as weights", call. = FALSE)
  }
  sig2 <- n / (n - 1) * sum(w * (y - sum(w * y) / sum(w))^2) / (sum(w) - 1)

  if (sig2 < .Machine$double.eps) {
    warning(
      "outcome variance is approximately zero; ",
      "CR design effect is undefined, returning 1 by convention",
      call. = FALSE
    )
    cv2h <- if (!is.null(strata_id)) {
      str_idx <- match(strata_id, sort(unique(strata_id)))
      H <- max(str_idx)
      nh <- tabulate(str_idx, nbins = H)
      w_grp <- split(w, str_idx)
      vapply(
        seq_len(H),
        function(i) {
          if (nh[i] <= 1L) {
            return(0)
          }
          wsub <- w_grp[[i]]
          (nh[i] - 1) / nh[i] * var(wsub) / mean(wsub)^2
        },
        numeric(1L)
      )
    } else {
      (n - 1) / n * var(w) / mean(w)^2
    }
    if (!is.null(strata_id)) {
      strat <- sort(unique(strata_id))
      H <- length(strat)
      str_idx <- match(strata_id, strat)
      nh <- tabulate(str_idx, nbins = H)
      na_vec <- rep(NA_real_, H)
      if (is.null(cluster_id)) {
        return(list(
          strata = data.frame(
            stratum = strat,
            n_h = nh,
            cv2_w = cv2h,
            deff_w = 1 + cv2h,
            deff_s = na_vec
          ),
          overall = 1
        ))
      } else {
        return(list(
          strata = data.frame(
            stratum = strat,
            n_h = nh,
            rho_h = na_vec,
            cv2_w = cv2h,
            deff_w = 1 + cv2h,
            deff_c = na_vec,
            deff_s = na_vec
          ),
          overall = 1
        ))
      }
    } else {
      if (is.null(cluster_id)) {
        return(list(
          strata = data.frame(n = n, cv2_w = cv2h, deff_w = 1 + cv2h),
          overall = 1
        ))
      } else {
        return(list(
          strata = data.frame(
            n = n,
            rho = NA_real_,
            cv2_w = cv2h,
            deff_w = 1 + cv2h,
            deff_c = NA_real_
          ),
          overall = 1
        ))
      }
    }
  }

  if (!is.null(strata_id)) {
    .deff_cr_stratified(w, y, strata_id, cluster_id, stages, n, sig2)
  } else {
    .deff_cr_unstratified(w, y, cluster_id, n, sig2)
  }
}

#' Linearization deff for a single stratum with clusters
#' V_COM / V_SRS using the ultimate cluster (Taylor) estimator.
#' @keywords internal
#' @noRd
.stratum_cluster_deff <- function(w, y, cluster_id, sig2h, nh) {
  wsum <- sum(w)
  ybar <- sum(w * y) / wsum

  cl_idx <- split(seq_along(w), cluster_id)
  a <- length(cl_idx)
  if (a <= 1L) {
    return(NA_real_)
  }

  z <- vapply(cl_idx, function(idx) sum(w[idx] * (y[idx] - ybar)), numeric(1L))
  v_com <- a / (a - 1) * sum(z^2) / wsum^2
  v_srs <- sig2h / nh
  if (v_srs < .Machine$double.eps) {
    return(NA_real_)
  }
  v_com / v_srs
}

#' CR: stratified path
#' @keywords internal
#' @noRd
.deff_cr_stratified <- function(w, y, strata_id, cluster_id, stages, n, sig2) {
  strat <- sort(unique(strata_id))
  H <- length(strat)
  str_idx <- match(strata_id, strat)
  nh <- tabulate(str_idx, nbins = H)

  if (!is.null(cluster_id) && (is.null(stages) || length(stages) != H)) {
    stop(
      "'stages' must have length equal to the number of strata",
      call. = FALSE
    )
  }

  w_grp <- split(w, str_idx)
  Wh <- vapply(w_grp, sum, numeric(1L)) / sum(w)

  y_grp <- split(y, str_idx)
  cv2h <- sig2h <- deff_s <- numeric(H)

  for (i in seq_len(H)) {
    if (nh[i] <= 1L) {
      cv2h[i] <- 0
      sig2h[i] <- 0
      deff_s[i] <- 0
      next
    }
    wsub <- w_grp[[i]]
    ysub <- y_grp[[i]]
    mw <- mean(wsub)
    sw <- sum(wsub)
    cv2h[i] <- (nh[i] - 1) / nh[i] * var(wsub) / mw^2
    sig2h[i] <- nh[i] /
      (nh[i] - 1) *
      sum(wsub * (ysub - sum(wsub * ysub) / sw)^2) /
      (sw - 1)
    deff_s[i] <- Wh[i]^2 / nh[i] * n * sig2h[i] / sig2
  }

  singletons <- which(nh <= 1L)
  if (length(singletons) > 0L) {
    warning(
      "singleton stratum/strata (n_h = 1) at position(s) ",
      paste(singletons, collapse = ", "),
      "; contribution set to zero",
      call. = FALSE
    )
  }

  deff_w <- 1 + cv2h

  if (is.null(cluster_id)) {
    out <- list(
      strata = data.frame(
        stratum = strat,
        n_h = nh,
        cv2_w = cv2h,
        deff_w = deff_w,
        deff_s = deff_s
      ),
      overall = sum(deff_w * deff_s)
    )
  } else {
    cl_grp <- split(cluster_id, str_idx)
    nh_star <- rhoh <- numeric(H)
    for (i in seq_len(H)) {
      wsub <- w_grp[[i]]
      clsub <- cl_grp[[i]]
      cl_wt_sums <- vapply(split(wsub, clsub), sum, numeric(1L))
      nh_star[i] <- sum(cl_wt_sums^2) / sum(wsub^2)
    }

    for (i in seq_len(H)) {
      if (stages[i] == 1L) {
        next
      }
      strdeff_i <- .stratum_cluster_deff(
        w_grp[[i]],
        y_grp[[i]],
        cl_grp[[i]],
        sig2h[i],
        nh[i]
      )
      if (is.na(strdeff_i)) {
        next
      }
      rhoh[i] <- (strdeff_i - deff_w[i]) / (deff_w[i] * (nh_star[i] - 1))
    }

    deff_c <- 1 + (nh_star - 1) * rhoh
    deff_c[stages == 1L] <- 1

    out <- list(
      strata = data.frame(
        stratum = strat,
        n_h = nh,
        rho_h = rhoh,
        cv2_w = cv2h,
        deff_w = deff_w,
        deff_c = deff_c,
        deff_s = deff_s
      ),
      overall = sum(deff_w * deff_c * deff_s)
    )
  }

  out
}

#' CR: unstratified path
#' @keywords internal
#' @noRd
.deff_cr_unstratified <- function(w, y, cluster_id, n, sig2) {
  cv2h <- (n - 1) / n * var(w) / mean(w)^2
  deff_w <- 1 + cv2h

  if (is.null(cluster_id)) {
    list(
      strata = data.frame(n = n, cv2_w = cv2h, deff_w = deff_w),
      overall = deff_w
    )
  } else {
    cl_wt_sums <- vapply(split(w, cluster_id), sum, numeric(1L))
    nh_star <- sum(cl_wt_sums^2) / sum(w^2)

    strdeff <- .stratum_cluster_deff(w, y, cluster_id, sig2, n)
    if (is.na(strdeff)) {
      rho <- 0
    } else {
      rho <- (strdeff - deff_w) / (deff_w * (nh_star - 1))
    }
    deff_c <- 1 + (nh_star - 1) * rho

    list(
      strata = data.frame(
        n = n,
        rho = rho,
        cv2_w = cv2h,
        deff_w = deff_w,
        deff_c = deff_c
      ),
      overall = deff_w * deff_c
    )
  }
}
