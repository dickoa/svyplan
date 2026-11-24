#' Validate overlap against ratio or explicit n
#' @keywords internal
#' @noRd
.check_overlap_n <- function(overlap, n = NULL, ratio = NULL) {
  if (overlap == 0) {
    return(invisible(NULL))
  }
  if (!is.null(ratio) && ratio > 1 && overlap > 1 / ratio) {
    stop(
      sprintf("overlap (%.3g) must be <= 1/ratio (%.3g)", overlap, 1 / ratio),
      call. = FALSE
    )
  }
  if (!is.null(n) && length(n) == 2L && overlap > n[2] / n[1]) {
    stop(
      sprintf(
        "overlap (%.3g) must be <= n[2]/n[1] (%.3g)",
        overlap,
        n[2] / n[1]
      ),
      call. = FALSE
    )
  }
}

#' Finite population factors of a change under each overlap draw
#'
#' The population terms of a change variance are
#' \eqn{-(m_1v_1 + m_2v_2 - 2cS_{12})/N}{-(m_1 v_1 + m_2 v_2 - 2 c S_12)/N}.
#' A fixed shared count gives \eqn{m_t = c = 1}. A retained part and fresh
#' parts drawn independently of it from the whole frame, each occasion
#' pooling its units with equal weight, give
#' \eqn{m_t = w_t^2 + (1 - w_t)^2}{m_t = w_t^2 + (1 - w_t)^2} and
#' \eqn{c = w_1w_2}{c = w_1 w_2}, where \eqn{w_t} is the retained share of
#' occasion \eqn{t}: `overlap` on the first and `overlap * size_ratio` on
#' the second, `size_ratio` being \eqn{n_1/n_2}. At zero overlap both reduce
#' to two independent samples.
#' @keywords internal
#' @noRd
.overlap_pop_factors <- function(overlap, size_ratio, overlap_draw) {
  if (overlap_draw == "fixed") {
    return(list(marginal = c(1, 1), cross = 1))
  }
  w <- overlap * c(1, size_ratio)
  list(marginal = w^2 + (1 - w)^2, cross = w[1L] * w[2L])
}

#' The overlap the frame must hold distinct units for
#'
#' Independently drawn parts each fit the frame on their own, so only a
#' fixed shared count makes the two occasions compete for distinct units.
#' @keywords internal
#' @noRd
.frame_overlap <- function(overlap, overlap_draw) {
  if (overlap_draw == "independent") 0 * overlap else overlap
}

#' Positive-overlap samples must fit their shared finite frame
#'
#' Sizes here count respondents. Zero overlap retains the public independent-
#' samples convention. Vector inputs check several pairs (DiD arms or pooled
#' lags), not the union of an entire rotation schedule.
#' @keywords internal
#' @noRd
.check_overlap_frame <- function(n1_net, n2_net, N, overlap,
                                 coupled = FALSE) {
  distinct <- n2_net + (1 - overlap) * n1_net
  bad <- which((overlap > 0 | coupled) & is.finite(N) &
                 distinct > N + 32 * .Machine$double.eps * pmax(1, N))
  if (length(bad)) {
    j <- bad[1L]
    # Recycle scalar arguments explicitly for an informative vector diagnostic.
    counts <- distinct[j]
    population <- rep_len(N, length(distinct))[j]
    stop(sprintf(
      "finite population overlap is infeasible: the two responding samples require %.10g distinct units but the population in scope ('N', times 'eligibility_rate' when given) is %.10g; increase 'overlap', reduce 'n', or relax the sizing target",
      counts, population
    ), call. = FALSE)
  }
  invisible(NULL)
}

#' The whole count a continuous size rounds up to
#'
#' A size solved by inverting a precision lands within rounding of a whole
#' number whenever that number is the size the precision came from, often a
#' step above it, so a raw `ceiling()` reports one unit more than the
#' design that meets the target. A value within a relative 1e-12, at most
#' 1e-6, of a positive whole number is taken as that number before rounding
#' up. Round trips stray by a few rounding steps, under 1e-15 relative, so
#' a real fraction still rounds up, a whole number keeps its value, and a
#' positive size below 1 still takes a unit.
#' @keywords internal
#' @noRd
.whole_count <- function(x) {
  k <- round(x)
  snap <- !is.na(x) & is.finite(x) & k >= 1 &
    abs(x - k) <= pmin(1e-12 * pmax(1, abs(x)), 1e-6)
  x[snap] <- k[snap]
  ceiling(x)
}

#' A solved size must fit the frame once rounded up
#'
#' The frame and covariance checks run on the continuous size, while
#' `print()` and `as.integer()` report `ceiling(n)`, which can need more
#' distinct units than the frame holds. `check` runs on the whole sizes the
#' checks a supplied size meets, and a failure means no whole size reaches
#' the target, since a smaller one misses it.
#' @keywords internal
#' @noRd
.check_whole_size <- function(n, check, unit = "") {
  whole <- .whole_count(n)
  if (all(whole == n)) {
    return(invisible(NULL))
  }
  tryCatch(check(whole), error = function(e) {
    stop(sprintf(
      "the target needs n = %s%s, which the frame holds only before rounding up. At the whole size %s, %s",
      paste(format(signif(n, 7)), collapse = " and "), unit,
      paste(whole, collapse = " and "), conditionMessage(e)
    ), call. = FALSE)
  })
  invisible(NULL)
}

#' Validate and normalize a take_all column
#'
#' Documented as logical or 0/1, so any other number is a mistake rather
#' than a truthy value. `as.logical()` on its own would take 2 for `TRUE`
#' and pin a stratum the caller never meant to pin.
#' @keywords internal
#' @noRd
.check_take_all <- function(take_all, n) {
  if (is.null(take_all)) {
    return(rep(FALSE, n))
  }
  if (is.numeric(take_all)) {
    if (anyNA(take_all) || any(!take_all %in% c(0, 1))) {
      stop("'take_all' must be logical (or 0/1)", call. = FALSE)
    }
    return(take_all != 0)
  }
  if (!is.logical(take_all) || anyNA(take_all)) {
    stop("'take_all' must be logical (or 0/1)", call. = FALSE)
  }
  take_all
}

#' Overlapping samples must share one population
#'
#' Two occasions can only have units in common if they are drawn from the
#' same frame, so a single `N` is the only coherent input once `overlap` is
#' positive. Without overlap the two groups are independent and may
#' legitimately be different populations of different sizes.
#' @keywords internal
#' @noRd
.check_overlap_N <- function(overlap, N_pair) {
  if (overlap > 0 && !isTRUE(all.equal(N_pair[1], N_pair[2]))) {
    stop(
      paste("overlapping samples are drawn from one population, so 'N' must",
            "be a single size or two equal ones"),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Design variance of a difference between two overlapping occasion means
#'
#' The marginal terms carry their own finite population correction, whereas
#' the overlap covariance does not. For two SRSWOR samples drawn from one
#' population of size \eqn{N} and sharing \eqn{k = overlap \cdot n_1}{k = overlap * n_1} units,
#' \deqn{Cov(\bar y_1, \bar y_2) = \rho S_1S_2\{k/(n_1n_2) - 1/N\},}{Cov(ybar_1, ybar_2) = rho S_1S_2\{k/(n_1n_2) - 1/N\},}
#' so the population term enters once as \eqn{1/N} rather than through
#' sample 2's marginal factor. Collecting terms,
#' \deqn{V = \frac{v_1}{n_1} + \frac{v_2}{n_2}
#'       - \frac{2\rho\,overlap\sqrt{v_1v_2}}{n_2}
#'       - \frac{v_1 + v_2 - 2\rho\sqrt{v_1v_2}}{N},}{V = v_1/n_1 + v_2/n_2 - (2 rho overlap sqrt(v_1v_2))/n_2 - (v_1 + v_2 - 2 rho sqrt(v_1v_2))/N,}
#' a per-unit part less a census part. At \eqn{\rho = 1} with equal sizes
#' and variances the census part is zero, the population terms cancel
#' exactly, and the whole reduces to \eqn{2S^2(1 - overlap)/n}.
#'
#' Under `overlap_draw = "independent"` the per-unit part is the same and
#' the census part carries the factors of `.overlap_pop_factors()`. That
#' path is kept apart so the fixed-count result does not move by a bit.
#' @keywords internal
#' @noRd
.diff_var_fpc <- function(n_eff, var_pair, N_pair, deff, overlap, overlap_cor,
                          overlap_draw) {
  if (length(n_eff) == 1L) n_eff <- c(n_eff, n_eff)
  if (overlap_draw == "independent" && overlap > 0) {
    f <- .overlap_pop_factors(overlap, n_eff[1L] / n_eff[2L], overlap_draw)
    pop <- if (is.infinite(N_pair[1L])) 0 else 1 / N_pair[1L]
    cross <- overlap_cor * sqrt(var_pair[1] * var_pair[2])
    V <- sum(var_pair * (1 / n_eff - f$marginal * pop)) -
      2 * cross * (overlap / n_eff[2L] - f$cross * pop)
    return(deff * V)
  }
  .check_overlap_frame(n_eff[1L], n_eff[2L], N_pair[1L], overlap)
  fpc1 <- .fpc_factor(n_eff[1], N_pair[1])
  fpc2 <- .fpc_factor(n_eff[2], N_pair[2])
  V <- var_pair[1] * fpc1 / n_eff[1] + var_pair[2] * fpc2 / n_eff[2]
  if (overlap > 0) {
    cross <- overlap_cor * sqrt(var_pair[1] * var_pair[2])
    V <- V - 2 * overlap * cross / n_eff[2] +
      if (is.infinite(N_pair[1])) 0 else 2 * cross / N_pair[1]
  }
  deff * V
}

#' Finite-population variance of a Bernoulli variable
#'
#' The variance consumed by the package's without-replacement kernels is the
#' population variance on `N - 1` degrees of freedom. For a population
#' proportion `p`, that is `N * p * (1 - p) / (N - 1)`. The infinite-
#' population limit is the usual `p * (1 - p)`.
#' @keywords internal
#' @noRd
.bernoulli_var <- function(p, N) {
  adj <- ifelse(is.infinite(N), 1, N / (N - 1))
  p * (1 - p) * adj
}

#' Resolve the two occasions' variances and the change they bound
#'
#' A change is planned on one of two scales. On the mean scale the
#' dispersion is supplied directly and the expected change is optional,
#' needed only by the relative measures. On the proportion scale the two
#' occasion proportions supply both: the variances are \eqn{p(1-p)} and the
#' change is \eqn{p_2 - p_1}, so supplying `change` as well would let the
#' two disagree and is refused rather than silently preferred.
#'
#' The proportion scale carries the same \eqn{N/(N-1)} adjustment
#' `.prec_engine_prop()` applies, because the variance the change formula
#' wants is the population variance \eqn{S^2} on \eqn{N-1} degrees of
#' freedom and a Bernoulli population's is \eqn{Np(1-p)/(N-1)}, not
#' \eqn{p(1-p)}. Without it a finite `N` would make one occasion of
#' `prec_change(p = )` disagree with `prec_prop()` on the same design. The
#' mean scale needs no adjustment: `var` is already defined on \eqn{N-1}.
#'
#' A ratio of totals is the third scale, through `.ratio_occasions()`. Its
#' variances and change come back divided by the ratio's scale, which is
#' returned as `scale` (1 on the other two) for the caller to multiply the
#' standard error and margin back. `var_stored` and `change_stored` are the
#' unscaled values a result records.
#' @keywords internal
#' @noRd
.change_inputs <- function(var, sd, p, change, N_pair, r = NULL,
                           cv_num = NULL, cv_den = NULL,
                           component_cor = NULL) {
  has_var <- !is.null(var) || !is.null(sd)
  has_p <- !is.null(p)
  has_r <- !is.null(r)
  .check_ratio_trio_needs_r(has_r, cv_num, cv_den, component_cor)
  if (has_var + has_p + has_r != 1L) {
    stop("specify exactly one of 'var' (or 'sd'), 'p', or 'r'", call. = FALSE)
  }
  if (has_r) {
    if (!is.null(change)) {
      stop(
        "'change' is r[2] - r[1] for a ratio; do not supply it",
        call. = FALSE
      )
    }
    rat <- .ratio_occasions(r, cv_num, cv_den, component_cor, 2L,
                            sizing = FALSE)
    return(list(
      var_pair = rat$var,
      change = r[2L] / rat$scale - r[1L] / rat$scale,
      p = NULL,
      ratio = rat,
      scale = rat$scale,
      var_stored = r^2 * rat$unit_relvar,
      change_stored = r[2L] - r[1L]
    ))
  }
  if (has_p) {
    if (!is.null(change)) {
      stop(
        "'change' is p[2] - p[1] on the proportion scale; do not supply it",
        call. = FALSE
      )
    }
    if (!is.numeric(p) || length(p) != 2L) {
      stop("'p' must be two proportions, one per occasion", call. = FALSE)
    }
    check_proportion(p[1L], "p[1]")
    check_proportion(p[2L], "p[2]")
    var_pair <- .bernoulli_var(p, N_pair)
    return(list(
      var_pair = var_pair,
      change = p[2L] - p[1L],
      p = p,
      scale = 1,
      var_stored = var_pair,
      change_stored = p[2L] - p[1L]
    ))
  }
  var <- .resolve_var(var, sd)
  if (!is.null(change)) check_scalar(change, "change", positive = FALSE)
  var_pair <- .as_pair(var, "var")
  list(var_pair = var_pair, change = change, p = NULL, scale = 1,
       var_stored = var_pair, change_stored = change)
}

#' The moments of a ratio of totals on each occasion of a repeated survey
#'
#' A ratio enters the change and pooled engines through its linearized
#' variable \eqn{z = (y - Rx)/\bar X}{z = (y - R x) / Xbar}, whose unit
#' variance is \eqn{R^2 L} with \eqn{L} the unit relative variance of the
#' scalar pair. The engines run on \eqn{z} divided by
#' \eqn{s = \max_t |r_t|}{s = max |r_t|}, so the variance handed over is
#' \eqn{(r_t/s)^2 L_t}{(r_t / s)^2 L_t}, and the caller multiplies the
#' standard error and margin back by \eqn{s}. Forming \eqn{r^2 L} itself
#' overflows for a large finite ratio and underflows for a small one, the
#' trap `.prec_engine_ratio()` avoids for one occasion.
#'
#' `k` is the number of occasions the ratio is stated for: two for a change,
#' whose ratios may differ, and one for a pooled average, whose level is
#' taken as the same on every occasion. The moments take one value or one
#' per occasion. Each occasion's unit relative variance may be zero for a
#' change, since the other can carry the variance, and `sizing` is passed
#' through for a pooled average, where one occasion is all of them.
#' @keywords internal
#' @noRd
.ratio_occasions <- function(r, cv_num, cv_den, component_cor, k, sizing) {
  if (is.null(cv_num) || is.null(cv_den) || is.null(component_cor)) {
    stop(
      "a ratio needs 'r', 'cv_num', 'cv_den', and 'component_cor'",
      call. = FALSE
    )
  }
  if (!is.numeric(r) || length(r) != k) {
    stop(
      if (k == 2L) {
        "'r' must be two ratios, one per occasion"
      } else {
        "'r' must be a single ratio, the level on every occasion"
      },
      call. = FALSE
    )
  }
  spread <- function(x, name) {
    if (!is.numeric(x) || !length(x) %in% unique(c(1L, k))) {
      stop(
        sprintf("'%s' takes %s", name, if (k == 2L) {
          "one value, or one per occasion"
        } else {
          "one value"
        }),
        call. = FALSE
      )
    }
    rep_len(x, k)
  }
  num <- spread(cv_num, "cv_num")
  den <- spread(cv_den, "cv_den")
  rho <- spread(component_cor, "component_cor")
  relvar <- numeric(k)
  for (t in seq_len(k)) {
    .check_ratio_moments(r[t], num[t], den[t], rho[t])
    relvar[t] <- .ratio_unit_relvar(r[t], num[t], den[t], rho[t])
    # Once per distinct set of moments, so moments shared by both occasions
    # warn once rather than twice.
    if (t == 1L || !identical(c(num[t], den[t], rho[t], sign(r[t])),
                              c(num[1L], den[1L], rho[1L], sign(r[1L])))) {
      .check_ratio_relvar(relvar[t], num[t], den[t], rho[t],
                          sizing = sizing)
    }
  }
  s <- max(abs(r))
  list(
    r = r,
    cv_num = cv_num,
    cv_den = cv_den,
    component_cor = component_cor,
    unit_relvar = relvar,
    scale = s,
    var = (r / s)^2 * relvar
  )
}

#' Record a ratio's moments on a change or pooled result
#'
#' The moments as supplied, plus the unit relative variance of each
#' occasion, which a reader would otherwise have to rederive from them.
#' `var` stays the variance of the linearized variable, so the field means
#' the same on every scale.
#' @keywords internal
#' @noRd
.add_ratio_params <- function(params, ratio) {
  if (is.null(ratio)) {
    return(params)
  }
  params$r <- ratio$r
  params$cv_num <- ratio$cv_num
  params$cv_den <- ratio$cv_den
  params$component_cor <- ratio$component_cor
  params$unit_relvar <- ratio$unit_relvar
  params
}

#' Refuse the ratio moments without the ratio they describe
#'
#' Without `r` the three would be ignored in favor of `var` or `p`, and a
#' plan the caller meant for a ratio would silently be sized for something
#' else.
#' @keywords internal
#' @noRd
.check_ratio_trio_needs_r <- function(has_r, cv_num, cv_den, component_cor) {
  if (!has_r &&
        !(is.null(cv_num) && is.null(cv_den) && is.null(component_cor))) {
    stop(
      "'cv_num', 'cv_den' and 'component_cor' describe a ratio and need 'r'",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The arguments that rebuild a change or pooled result on its own scale
#'
#' One of three: the ratio's moments, the proportions, or the variance with
#' its change or level. Read with `[[` throughout, since `$r` would
#' partial-match `resp_rate` on a result that has no ratio.
#' @keywords internal
#' @noRd
.repeated_scale_args <- function(par, level) {
  if (!is.null(par[["r"]])) {
    return(list(
      r = par[["r"]],
      cv_num = par[["cv_num"]],
      cv_den = par[["cv_den"]],
      component_cor = par[["component_cor"]]
    ))
  }
  if (!is.null(par[["p"]])) {
    return(list(p = par[["p"]]))
  }
  out <- list(var = par[["var"]])
  out[[level]] <- par[[level]]
  out
}

#' Is a change or pooled result on the variance scale, where its level or
#' change is a free input rather than determined by `p` or `r`?
#' @keywords internal
#' @noRd
.repeated_free_level <- function(par) {
  is.null(par[["p"]]) && is.null(par[["r"]])
}

#' Refuse a between-occasion correlation two proportions cannot have
#'
#' Two Bernoulli variables with means \eqn{p_1} and \eqn{p_2} admit
#' correlations only up to the Frechet-Hoeffding bound
#' \eqn{(\min(p_1,p_2) - p_1p_2)/\sqrt{p_1q_1p_2q_2}}{(min(p_1,p_2) - p_1p_2)/sqrt(p_1q_1p_2q_2)}. Above it no joint
#' distribution exists, so the covariance the change formula subtracts
#' describes nothing. Checked only where the marginals are known, which is
#' the proportion scale, and only where the correlation is used at all.
#' @keywords internal
#' @noRd
.check_bernoulli_cor <- function(p, overlap_cor, overlap) {
  if (is.null(p) || overlap == 0 || overlap_cor == 0) {
    return(invisible(TRUE))
  }
  bound <- (min(p) - p[1L] * p[2L]) /
    sqrt(p[1L] * (1 - p[1L]) * p[2L] * (1 - p[2L]))
  if (overlap_cor > bound * (1 + 1e-9)) {
    stop(
      sprintf(
        "'overlap_cor' (%.4g) exceeds the largest correlation two proportions of %.4g and %.4g can have (%.4g)",
        overlap_cor, p[1L], p[2L], bound
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Feasible second marginal at a requested Bernoulli correlation
#' @keywords internal
#' @noRd
.bernoulli_cor_interval <- function(p1, overlap_cor, overlap) {
  if (overlap == 0 || overlap_cor == 0) {
    return(c(0, 1))
  }
  rho2 <- overlap_cor^2
  q1 <- 1 - p1
  c(
    rho2 * p1 / (q1 + rho2 * p1),
    p1 / (p1 + rho2 * q1)
  )
}

#' Check a size supplied per occasion
#'
#' One size means both occasions are the same size, which is the ordinary
#' repeated survey. Two mean the second occasion was resized, which is what
#' makes `overlap` directional.
#' @keywords internal
#' @noRd
.check_change_n <- function(n) {
  if (!is.numeric(n) || anyNA(n) || !all(is.finite(n))) {
    stop("'n' must be finite numeric", call. = FALSE)
  }
  if (!length(n) %in% c(1L, 2L)) {
    stop(
      "'n' must be length 1 (equal occasions) or 2 (one size per occasion)",
      call. = FALSE
    )
  }
  if (any(n <= 0)) {
    stop("'n' must be positive", call. = FALSE)
  }
  n
}

#' Sampling variance of a change between two occasions of one population
#'
#' A thin wrapper over `.diff_var_fpc()` that nets the gross sizes down and
#' pairs a scalar `n`, so the size the user states is the size every entry
#' point in the package states: units drawn, not interviews completed.
#' @keywords internal
#' @noRd
.change_var <- function(n, var_pair, N_pair, deff, resp_rate, overlap,
                        overlap_cor, overlap_draw) {
  n_vec <- if (length(n) == 1L) c(n, n) else n
  V <- .diff_var_fpc(
    n_vec * resp_rate, var_pair, N_pair, deff, overlap, overlap_cor,
    overlap_draw
  )
  .safe_variance(V, "change variance")
}

#' Precision of a change at a given size
#' @keywords internal
#' @noRd
.prec_engine_change <- function(n, var_pair, change, N_pair, deff, resp_rate,
                                overlap, overlap_cor, overlap_draw, alpha,
                                df = NULL) {
  se <- sqrt(
    .change_var(n, var_pair, N_pair, deff, resp_rate, overlap, overlap_cor,
                overlap_draw)
  )
  list(
    se = se,
    moe = .q_alpha(alpha, df) * se,
    cv = .cv_of(se, change)
  )
}

#' Refuse to size a change no sample size can measure
#'
#' `A` is the per-unit variance of the difference before any finite
#' population term, and `scale` the sum of the two variances it is built
#' from. At full overlap and correlation the difference has no variance at
#' any size, so a target precision or power identifies no size.
#' @keywords internal
#' @noRd
.check_change_variance <- function(A, scale) {
  if (A <= sqrt(.Machine$double.eps) * scale) {
    stop(
      "the change carries no sampling variance at any size (overlap and overlap_cor leave nothing to sample); reduce either, or size the occasions unequally",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Second-occasion size that reaches a target standard error on the change
#'
#' The change variance is \eqn{A/n_2 + B} in the net second-occasion size,
#' with the finite population terms collected into a constant \eqn{B} that
#' no sample size can move. Both coefficients are exact on \eqn{n \le N},
#' where the correction factor is still \eqn{1 - n/N}, so the inversion is a
#' division rather than a search. \eqn{B \le 0} always, since the population
#' term it carries is \eqn{-(v_1 + v_2 - 2\rho\sqrt{v_1v_2})/N}{-(v_1 + v_2 - 2 rho sqrt(v_1v_2))/N} at the common
#' \eqn{N} that a positive overlap requires, so the divisor cannot vanish.
#' Independent draws keep it so, since \eqn{m_t \ge w_t^2} in the factors
#' of `.overlap_pop_factors()`, which are constant in \eqn{n} at a fixed
#' `ratio`.
#'
#' \eqn{A} can vanish. At equal sizes, equal variances, full overlap and unit
#' correlation the two occasions share every unit and the change is measured
#' without sampling error at any size, so no size is identified and the
#' caller is told which input to relax rather than handed a zero.
#' @keywords internal
#' @noRd
.n_change_from_se <- function(target_se, var_pair, N_pair, deff, resp_rate,
                              ratio, overlap, overlap_cor, overlap_draw) {
  cross <- if (overlap > 0) overlap_cor * sqrt(var_pair[1] * var_pair[2]) else 0
  f <- .overlap_pop_factors(overlap, ratio, overlap_draw)
  A <- var_pair[1] / ratio + var_pair[2] - 2 * overlap * cross
  B <- -f$marginal[1] * var_pair[1] / N_pair[1] -
    f$marginal[2] * var_pair[2] / N_pair[2] +
    if (is.infinite(N_pair[1])) 0 else 2 * f$cross * cross / N_pair[1]
  .check_change_variance(A, sum(var_pair))
  n2_net <- A / (target_se^2 / deff - B)
  n2 <- n2_net / resp_rate
  if (ratio == 1) n2 else c(ratio * n2, n2)
}

#' Number of occasions entering a pooled estimate
#' @keywords internal
#' @noRd
.MAX_OCCASIONS <- 1000L

.check_occasions <- function(occasions) {
  if (
    !is.numeric(occasions) || length(occasions) != 1L || anyNA(occasions) ||
      !is.finite(occasions) || occasions < 2 ||
      occasions != round(occasions)
  ) {
    stop("'occasions' must be a whole number of at least 2", call. = FALSE)
  }
  # Dense assembly is cubic, so a mistyped figure fails rather than allocates.
  if (occasions > .MAX_OCCASIONS) {
    stop(
      sprintf(
        "'occasions' must be at most %d; the covariance is assembled densely, so a larger horizon is a typing error more often than a design",
        .MAX_OCCASIONS
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve the overlap at every lag a pooled estimate spans
#'
#' A bare number is respondent overlap, the reading every other entry point
#' in the package uses. A `svyplan_overlap` from a rotation is issued
#' overlap, which is a different quantity below full response, so it is
#' accepted only where the two coincide and refused with the conversion named
#' otherwise. One from a panel plan is already respondent overlap. The
#' profile is recomputed from the rotation it carries at every lag the
#' horizon spans, since `max_lag` may have cut it short of the life, and
#' beyond the life nothing is shared.
#' @keywords internal
#' @noRd
.resolve_lag_overlap <- function(overlap, occasions, resp_rate) {
  max_lag <- occasions - 1L
  if (inherits(overlap, "svyplan_overlap")) {
    ov_basis <- .overlap_basis(overlap)
    .check_overlap_basis(ov_basis, resp_rate)
    prof <- .overlap_profile(
      as.double(.overlap_rotation(overlap)),
      seq_len(max_lag),
      attr(overlap, "response", exact = TRUE)
    )
    return(list(overlap = prof$shared / prof$in_sample, basis = ov_basis,
                rotation = overlap))
  }
  if (
    !is.numeric(overlap) || anyNA(overlap) || any(overlap < 0) ||
      any(overlap > 1)
  ) {
    stop("'overlap' must be numbers in [0, 1]", call. = FALSE)
  }
  # An issued profile reaches here only as its object, which round trips
  # and grids pass back whole, so a bare number is respondent overlap.
  basis <- "respondent"
  if (length(overlap) == 1L) {
    return(list(overlap = rep(overlap, max_lag), basis = basis,
                rotation = NULL))
  }
  if (length(overlap) != max_lag) {
    stop(
      sprintf(
        "'overlap' must be one number or one per lag (%d for %d occasions)",
        max_lag, occasions
      ),
      call. = FALSE
    )
  }
  list(overlap = overlap, basis = basis, rotation = NULL)
}

#' An issued overlap is only a respondent overlap at full response
#'
#' The refusal names the conversion rather than performing it, since which
#' assumption to make about response persisting across occasions is the
#' planner's to state. Checked wherever a profile is resolved, so a stored
#' one meets it again on a round trip or in a grid.
#' @keywords internal
#' @noRd
.check_overlap_basis <- function(basis, resp_rate) {
  if (identical(basis, "issued") && !isTRUE(all.equal(resp_rate, 1))) {
    stop(
      "design_overlap() of a rotation reports the overlap between issued samples, and this variance is formed on respondents; the two coincide only at resp_rate = 1. Under loss that is permanent, design_overlap() of the rotating n_panel() plan gives the respondent overlap, and recruitment response cancels from it. Under response independent between occasions it is about resp_rate times the issued figure. Supply the one whose assumption holds",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve the correlation at every lag
#'
#' A panel's correlation falls away with distance, so a single number is the
#' wrong shape for anything but a first look and `cor_decay` is the form most
#' panels have. The two are alternatives rather than a value and a modifier,
#' so supplying both is refused.
#' @keywords internal
#' @noRd
.resolve_lag_cor <- function(overlap_cor, cor_decay, max_lag) {
  if (!is.null(cor_decay)) {
    if (!is.null(overlap_cor)) {
      stop("supply exactly one of 'overlap_cor' or 'cor_decay'", call. = FALSE)
    }
    if (
      !is.numeric(cor_decay) || length(cor_decay) != 1L || anyNA(cor_decay) ||
        cor_decay < 0 || cor_decay > 1
    ) {
      stop("'cor_decay' must be a number in [0, 1]", call. = FALSE)
    }
    return(cor_decay^seq_len(max_lag))
  }
  cor <- overlap_cor %||% 0
  if (!is.numeric(cor) || anyNA(cor) || any(cor < 0) || any(cor > 1)) {
    stop("'overlap_cor' must be numbers in [0, 1]", call. = FALSE)
  }
  if (length(cor) == 1L) {
    return(rep(cor, max_lag))
  }
  if (length(cor) != max_lag) {
    stop(
      sprintf(
        "'overlap_cor' must be one number or one per lag (%d for %d occasions)",
        max_lag, max_lag + 1L
      ),
      call. = FALSE
    )
  }
  cor
}

#' Finite population factors of a pooled estimate under each overlap draw
#'
#' The factors of `.overlap_pop_factors()` read lag by lag: `marginal`
#' multiplies \eqn{1/N} in an occasion's variance and `shared[m]` in the
#' covariance at lag `m`. Independent draws are defined at two occasions
#' only, where the lag-1 overlap is the retained share of both, and
#' `.check_pooled_draw()` refuses them beyond.
#' @keywords internal
#' @noRd
.pooled_pop_factors <- function(ov, overlap_draw) {
  if (overlap_draw == "fixed") {
    return(list(marginal = 1, shared = rep(1, length(ov))))
  }
  f <- .overlap_pop_factors(ov[1L], 1, overlap_draw)
  list(marginal = f$marginal[1L], shared = ov^2)
}

#' Independent draws pool over two occasions only
#'
#' Beyond two, an occasion's variance carries the sum of its panels'
#' squared shares and a lag covariance the same sum over the panels it
#' shares, which the lag overlaps do not determine. Without overlap at any
#' lag the two draws are the same design and nothing is refused.
#' @keywords internal
#' @noRd
.check_pooled_draw <- function(overlap_draw, ov, occasions) {
  if (overlap_draw == "independent" && occasions > 2 && any(ov > 0)) {
    stop(
      "overlap_draw = \"independent\" pools two occasions only. Beyond two, each occasion's variance depends on how its sample splits into independently drawn panels, which the lag overlaps do not determine. Use overlap_draw = \"fixed\", or plan two occasions",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The overlap a stored pooled result passes back
#'
#' The rotation object when the result was built from one, so its zero lags
#' stay disjoint cohorts and a new horizon re-derives every lag from the
#' rotation. A caller who replaces `overlap` decides by the replacement's
#' type, and a bare number then clears the rotation. The object carries
#' its own basis, so a round trip that lowers the response rate meets the
#' refusal the first call met.
#' @keywords internal
#' @noRd
.stored_overlap_args <- function(par, replaced = FALSE) {
  if (replaced) {
    return(list(overlap = par$overlap))
  }
  rot <- par[["overlap_rotation"]]
  if (!is.null(rot)) {
    return(list(overlap = rot))
  }
  list(overlap = par[["overlap_input"]] %||% par$overlap)
}

#' The stored lag correlation in the form it was given
#'
#' An override of either form replaces both, since the two cannot be
#' supplied together.
#' @keywords internal
#' @noRd
.stored_lag_cor_args <- function(par, overriding) {
  if (any(c("overlap_cor", "cor_decay") %in% overriding)) {
    return(list())
  }
  if (!is.null(par[["cor_decay"]])) {
    return(list(cor_decay = par[["cor_decay"]]))
  }
  list(overlap_cor = par[["overlap_cor_input"]] %||% par$overlap_cor)
}

#' Whether the occasions are cohorts of one rotation
#'
#' Decided by the input's type, never by its values. A [design_overlap()]
#' object describes a rotation, whose cohorts are drawn without replacement
#' from each other, so under a fixed draw a lag that shares no unit still
#' carries the population term of two disjoint samples, whatever the
#' horizon. A bare number describes overlap only, and zero means independent
#' draws, which keeps `V_level / occasions` exact for fresh samples.
#' Reading the relationship from positive lags instead made a pair's
#' covariance depend on how many occasions were requested.
#' @keywords internal
#' @noRd
.pooled_coupled <- function(rotation, overlap_draw) {
  !is.null(rotation) && overlap_draw == "fixed"
}

#' Covariance between two occasion means at each lag
#'
#' The kernel `.diff_var_fpc()` already commits the package to, read at a
#' general lag. A zero lag contributes exactly zero, population term
#' included, unless the occasions are cohorts of one rotation
#' (`.pooled_coupled()`), where it keeps \eqn{-\rho_m S^2/N}.
#' @keywords internal
#' @noRd
.pooled_lag_cov <- function(var, n_net, N, ov, rho, overlap_draw,
                            rotation = NULL) {
  pop <- if (is.infinite(N)) 0 else 1 / N
  shared <- .pooled_pop_factors(ov, overlap_draw)$shared
  out <- rho * var * (ov / n_net - shared * pop)
  if (!.pooled_coupled(rotation, overlap_draw)) out[ov <= 0] <- 0
  out
}

#' The occasion-by-occasion covariance a rotation implies
#' @keywords internal
#' @noRd
.pooled_kernel <- function(var, n_net, N, occasions, ov, rho, overlap_draw,
                           rotation = NULL) {
  c0 <- if (overlap_draw == "fixed") {
    var * .fpc_factor(n_net, N) / n_net
  } else {
    pop <- if (is.infinite(N)) 0 else 1 / N
    var * (1 / n_net - .pooled_pop_factors(ov, overlap_draw)$marginal * pop)
  }
  lag <- abs(outer(seq_len(occasions), seq_len(occasions), "-"))
  matrix(
    c(c0, .pooled_lag_cov(var, n_net, N, ov, rho, overlap_draw,
                          rotation))[lag + 1L],
    occasions, occasions
  )
}

#' The assembled covariance must be a covariance
#'
#' Checked on the matrix and not on the quantity reported from it. A kernel
#' that is not positive semidefinite can still give a positive pooled
#' variance, so `.safe_variance()` on the result passes cases this rejects.
#' The failure is ordinary rather than exotic: a lag's covariance changes
#' sign once its overlap falls below the sampling fraction, and a Toeplitz
#' kernel whose signs vary across lags need not be positive semidefinite.
#' An AR(1) correlation does not escape it, since the assembled matrix is a
#' difference of two positive semidefinite parts and not the Schur product
#' of the overlap and correlation kernels alone.
#' @keywords internal
#' @noRd
.check_lag_psd <- function(K, n_gross, n_net, N, ov, rho, occasions,
                           overlap_draw, rotation = NULL) {
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  diag_term <- K[1L, 1L]
  # Relative to the spectrum: an absolute floor would tie the verdict to the
  # outcome's units. The zero matrix has scale zero and passes on equality.
  scale <- max(abs(ev))
  if (min(ev) >= -sqrt(.Machine$double.eps) * scale) {
    return(invisible(NULL))
  }
  cov_m <- .pooled_lag_cov(1, n_net, N, ov, rho, overlap_draw, rotation)
  neg <- which(cov_m < 0)
  msg <- sprintf(
    "the assembled between-occasion covariance is not a valid covariance (minimum eigenvalue %.3g against a diagonal of %.3g)",
    min(ev), diag_term
  )
  if (length(neg) && !is.infinite(N)) {
    msg <- paste0(
      msg,
      sprintf(
        "; at a sampling fraction of %.3g the covariance at lag %d is negative, which happens once a lag's overlap falls below n/N",
        n_net / N, neg[1L]
      )
    )
  }
  # Interviews per cohort, not calendar span. Incomplete inside the horizon,
  # where the diagnostic is dropped.
  if (!is.null(rotation)) {
    # The object knows its interviews whatever the horizon.
    complete <- TRUE
    interviews <- sum(as.double(.overlap_rotation(rotation)) > 0)
  } else {
    complete <- length(ov) > 0L && ov[length(ov)] == 0
    interviews <- 1 + 2 * sum(ov)
  }
  if (complete && !is.infinite(N) && interviews > 0 &&
        occasions * n_gross / interviews > N) {
    msg <- paste0(
      msg,
      sprintf(
        ". Over %d occasions a rotation interviewing each cohort %.3g times at this size draws %.2f times the population, so the schedule could not be fielded either",
        occasions, interviews, occasions * n_gross / (interviews * N)
      )
    )
  }
  stop(
    paste0(
      msg,
      ". Reduce the sampling fraction, shorten the horizon, or supply an overlap and correlation that hold together across lags"
    ),
    call. = FALSE
  )
}

#' Variance of the equal-weight mean of the occasion estimates
#' @keywords internal
#' @noRd
.pooled_var <- function(n, var, N, deff, resp_rate, occasions, ov, rho,
                        overlap_draw, rotation = NULL) {
  n_net <- n * resp_rate
  K <- .pooled_kernel(var, n_net, N, occasions, ov, rho, overlap_draw,
                      rotation)
  .check_lag_psd(K, n, n_net, N, ov, rho, occasions, overlap_draw, rotation)
  .check_overlap_frame(n_net, n_net, N, .frame_overlap(ov, overlap_draw),
                       coupled = .pooled_coupled(rotation, overlap_draw))
  .safe_variance(deff * sum(K) / occasions^2, "pooled variance")
}

#' Precision of a pooled estimate at a given size
#' @keywords internal
#' @noRd
.prec_engine_pooled <- function(n, var, mu, N, deff, resp_rate, occasions,
                                ov, rho, overlap_draw, alpha, df = NULL,
                                rotation = NULL) {
  se <- sqrt(.pooled_var(n, var, N, deff, resp_rate, occasions, ov, rho,
                         overlap_draw, rotation))
  list(
    se = se,
    moe = .q_alpha(alpha, df) * se,
    cv = .cv_of(se, mu)
  )
}

#' Gross size per occasion that reaches a target standard error
#'
#' The pooled variance is \eqn{deff\{A/(n r) + B\}} in the gross size, with
#' \eqn{r} the response rate and the finite population terms collected into
#' \eqn{B}. \eqn{B \le 0} always, since its bracket is at least
#' `occasions` under the nonnegative correlation contract, and at least
#' half that under independent draws, so the divisor
#' cannot vanish and no target is unattainable for want of precision. The
#' size can still exceed \eqn{N}, which is the ordinary boundary and is
#' checked where every other size is.
#' @keywords internal
#' @noRd
.n_pooled_from_se <- function(target_se, var, N, deff, resp_rate, occasions,
                              ov, rho, overlap_draw, rotation = NULL) {
  m <- seq_len(occasions - 1L)
  w <- 2 * (occasions - m)
  pos <- ov > 0 | .pooled_coupled(rotation, overlap_draw)
  f <- .pooled_pop_factors(ov, overlap_draw)
  A <- var * (occasions + sum(w[pos] * rho[pos] * ov[pos])) / occasions^2
  B <- if (is.infinite(N)) {
    0
  } else {
    -var * (occasions * f$marginal + sum(w[pos] * rho[pos] * f$shared[pos])) /
      (N * occasions^2)
  }
  deff * A / (resp_rate * (target_se^2 - deff * B))
}

#' Resolve the dispersion and level a pooled estimate is planned on
#'
#' The proportion scale carries the same \eqn{N/(N-1)} adjustment
#' `.change_inputs()` applies and for the same reason: the variance the
#' engine wants is \eqn{S^2} on \eqn{N-1} degrees of freedom.
#'
#' A ratio of totals is the third scale, at one level for every occasion as
#' the other two assume. As in `.change_inputs()`, its variance and level
#' come back divided by the ratio's scale, returned as `scale`.
#' @keywords internal
#' @noRd
.pooled_inputs <- function(var, sd, p, mu, N, r = NULL, cv_num = NULL,
                           cv_den = NULL, component_cor = NULL, sizing = TRUE) {
  has_var <- !is.null(var) || !is.null(sd)
  has_p <- !is.null(p)
  has_r <- !is.null(r)
  .check_ratio_trio_needs_r(has_r, cv_num, cv_den, component_cor)
  if (has_var + has_p + has_r != 1L) {
    stop("specify exactly one of 'var' (or 'sd'), 'p', or 'r'", call. = FALSE)
  }
  if (has_r) {
    if (!is.null(mu)) {
      stop("'mu' is 'r' for a ratio; do not supply it", call. = FALSE)
    }
    rat <- .ratio_occasions(r, cv_num, cv_den, component_cor, 1L,
                            sizing = sizing)
    return(list(
      var = rat$var,
      mu = r / rat$scale,
      p = NULL,
      ratio = rat,
      scale = rat$scale,
      var_stored = r^2 * rat$unit_relvar,
      mu_stored = r
    ))
  }
  if (has_p) {
    if (!is.null(mu)) {
      stop(
        "'mu' is 'p' on the proportion scale; do not supply it",
        call. = FALSE
      )
    }
    check_proportion(p, "p")
    adj <- if (is.infinite(N)) 1 else N / (N - 1)
    v <- p * (1 - p) * adj
    return(list(var = v, mu = p, p = p, scale = 1, var_stored = v,
                mu_stored = p))
  }
  var <- .resolve_var(var, sd)
  check_scalar(var, "var")
  if (!is.null(mu)) check_scalar(mu, "mu", positive = FALSE)
  list(var = var, mu = mu, p = NULL, scale = 1, var_stored = var,
       mu_stored = mu)
}

#' Null-coalescing operator
#' Internal fallback for R versions before base::`%||%` became available.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

