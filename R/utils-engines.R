#' Normalize scalar or length-2 input to always length-2
#' @keywords internal
#' @noRd
.as_pair <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || anyNA(x) || !all(is.finite(x))) {
    stop(sprintf("'%s' must be finite numeric", name), call. = FALSE)
  }
  if (length(x) == 1L) {
    if (positive && x <= 0) {
      stop(sprintf("'%s' must be positive", name), call. = FALSE)
    }
    return(c(x, x))
  }
  if (length(x) == 2L) {
    if (positive && any(x <= 0)) {
      stop(sprintf("all '%s' elements must be positive", name), call. = FALSE)
    }
    return(x)
  }
  stop(sprintf("'%s' must be length 1 or 2", name), call. = FALSE)
}

#' Validate and normalize n for power functions
#' @keywords internal
#' @noRd
.check_power_n <- function(n) {
  if (!is.numeric(n) || anyNA(n) || !all(is.finite(n))) {
    stop("'n' must be finite numeric", call. = FALSE)
  }
  if (!length(n) %in% c(1L, 2L)) {
    stop(
      "'n' must be length 1 (equal groups) or 2 (unequal groups)",
      call. = FALSE
    )
  }
  if (any(n < 2)) {
    stop("'n' must be >= 2", call. = FALSE)
  }
  n
}

#' Resolve n/ratio when solving for n
#' @keywords internal
#' @noRd
.resolve_ratio <- function(n, ratio) {
  if (is.null(ratio)) {
    return(1)
  }
  if (
    !is.numeric(ratio) ||
      length(ratio) != 1L ||
      is.na(ratio) ||
      ratio <= 0 ||
      !is.finite(ratio)
  ) {
    stop("'ratio' must be a positive finite scalar", call. = FALSE)
  }
  if (!is.null(n) && ratio != 1) {
    stop("'ratio' cannot be used when 'n' is provided", call. = FALSE)
  }
  ratio
}

#' Two-sided interval quantile, normal or t
#'
#' Every confidence interval in the package is a half-width `q * se`, and
#' `q` is the normal quantile when the variance is treated as known and the
#' t quantile on `df` when it is estimated from a design with that many
#' degrees of freedom. Routing every interval site through one function is
#' what keeps a target and the sensitivity of that target to it reading the
#' same quantile: they are the same expression differentiated, and a
#' mismatch between them misreports without failing.
#'
#' `.z_alpha()` below is deliberately separate and stays normal-only. It
#' serves the power functions, where the quantile is a normal deviate for
#' an alternative rather than an interval half-width, and a t-based power
#' calculation is a different procedure than a quantile substitution.
#' @keywords internal
#' @noRd
.q_alpha <- function(alpha, df = NULL) {
  q <- qnorm(1 - alpha / 2)
  if (is.null(df) || all(is.na(df))) {
    return(q)
  }
  df <- rep_len(as.double(df), length(alpha))
  use <- !is.na(df)
  check_df(df[use], scalar = FALSE)
  q[use] <- stats::qt(1 - alpha[use] / 2, df[use])
  q
}

#' Check degrees of freedom for an interval quantile
#'
#' Student's t exists for every positive df. The package asks for at least
#' one because a design-based variance estimator with less than one degree
#' of freedom is not a design a planner can field. `Inf` is admitted and
#' gives the normal quantile. Entry points take one value, `scalar = FALSE`
#' serves the per-row columns of a targets table.
#' @keywords internal
#' @noRd
check_df <- function(df, name = "df", scalar = TRUE) {
  if (scalar) {
    if (!is.numeric(df) || length(df) != 1L || is.na(df) || df < 1) {
      stop(sprintf("'%s' must be a single number >= 1", name), call. = FALSE)
    }
  } else if (!is.numeric(df) || anyNA(df) || any(df < 1)) {
    stop(sprintf("'%s' values must be numbers >= 1", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Refuse degrees of freedom where the quantile is not an interval
#'
#' A power calculation's quantile is a normal deviate for an alternative,
#' not the half-width of a confidence interval, so substituting a t
#' quantile for it does not produce a t-based power calculation. That is a
#' different procedure. The argument is refused with the reason rather
#' than reported as an unknown name.
#' @keywords internal
#' @noRd
.stop_power_df <- function(...) {
  if ("df" %in% (...names() %||% character(0))) {
    stop(
      "'df' is not accepted by the power functions: their quantile is a normal deviate for an alternative, not an interval half-width, and a t-based power calculation is a different procedure rather than a substituted quantile. It applies to n_prop(), n_mean(), n_alloc() and their precision counterparts",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Compute z_alpha for power functions
#' @keywords internal
#' @noRd
.z_alpha <- function(alpha, alternative) {
  qnorm(1 - alpha / if (alternative == "two.sided") 2 else 1)
}

#' Check power population size N (scalar or length-2)
#' @keywords internal
#' @noRd
.check_power_N <- function(N) {
  if (!is.numeric(N) || anyNA(N) || !all(is.finite(N) | is.infinite(N))) {
    stop("'N' must be numeric, finite, or Inf", call. = FALSE)
  }
  if (!length(N) %in% c(1L, 2L)) {
    stop(
      "'N' must be length 1 (equal groups) or 2 (group-specific)",
      call. = FALSE
    )
  }
  if (any(N <= 1)) {
    stop("'N' values must be greater than 1 (or Inf)", call. = FALSE)
  }
  if (length(N) == 1L) c(N, N) else N
}

#' Effective N for n2 when ratio constraint links n1 = ratio * n2
#' @keywords internal
#' @noRd
.N2_for_ratio <- function(N_pair, ratio) {
  min(N_pair[2], N_pair[1] / ratio)
}

#' Upper bound for gross n2 from frame and respondent-overlap constraints
#'
#' `N_pair` is the population in scope and `resp_rate` the share of an
#' issued unit that responds in scope, so the overlap bound counts
#' respondents against `N_pair`. An issued sample also holds the units out
#' of scope, so the gross bound is the frame, `N_frame`.
#' @keywords internal
#' @noRd
.n2_upper_bound <- function(N_pair, ratio, resp_rate, overlap = 0,
                            within_arm = FALSE, N_frame = N_pair) {
  ratio_cap <- function(N) {
    min(if (is.infinite(N[1])) Inf else N[1] / ratio, N[2])
  }
  cap <- ratio_cap(N_frame)
  if (overlap > 0) {
    overlap_cap <- if (within_arm) {
      ratio_cap(N_pair) / (resp_rate * (2 - overlap))
    } else {
      N_pair[1L] / (resp_rate * (1 + ratio * (1 - overlap)))
    }
    cap <- min(cap, overlap_cap)
  }
  cap
}

#' Solve n2 by inverting a monotone power function
#' @keywords internal
#' @noRd
.solve_n2_from_power <- function(
  target_power,
  power_fn,
  N_pair,
  ratio,
  resp_rate,
  tol = 1e-8,
  overlap = 0,
  within_arm = FALSE,
  N_frame = N_pair
) {
  lo <- max(2, 2 / ratio)
  hi_cap <- .n2_upper_bound(N_pair, ratio, resp_rate, overlap, within_arm,
                            N_frame)
  if (hi_cap < lo) {
    stop("target power is unattainable under finite population and overlap constraints",
         call. = FALSE)
  }
  p_lo <- suppressWarnings(power_fn(lo))
  if (!is.finite(p_lo)) {
    p_lo <- 0
  }
  if (target_power <= p_lo) {
    return(lo)
  }

  if (is.finite(hi_cap)) {
    hi <- hi_cap
    if (hi <= lo) {
      stop(
        "target power is unattainable under finite population and overlap constraints",
        call. = FALSE
      )
    }
    p_hi <- suppressWarnings(power_fn(hi))
    if (!is.finite(p_hi) || p_hi < target_power - 32 * .Machine$double.eps) {
      stop(
        "target power is unattainable under finite population and overlap constraints",
        call. = FALSE
      )
    }
    if (abs(p_hi - target_power) <= 32 * .Machine$double.eps) return(hi)
  } else {
    hi <- max(2, lo * 2)
    p_hi <- suppressWarnings(power_fn(hi))
    iter <- 0L
    while ((is.na(p_hi) || p_hi < target_power) && iter < 100L) {
      hi <- hi * 2
      p_hi <- suppressWarnings(power_fn(hi))
      iter <- iter + 1L
    }
    if (!is.finite(p_hi) || p_hi < target_power) {
      stop("could not bracket sample size for target power", call. = FALSE)
    }
  }

  root <- uniroot(
    function(n2) power_fn(n2) - target_power,
    interval = c(lo, hi),
    tol = tol
  )$root
  # The root is found to within tol, either side, so a target that a whole
  # size reaches exactly comes back a fraction above or below it. The whole
  # size is returned when it reaches the target to within rounding.
  k <- round(root)
  if (k >= lo && k <= hi && abs(root - k) <= 2 * tol &&
        power_fn(k) >= target_power - 32 * .Machine$double.eps) {
    return(k)
  }
  root
}

#' Validate a power target used for inversion
#' @keywords internal
#' @noRd
.check_power_target <- function(power, alpha) {
  if (power <= alpha) {
    stop(
      sprintf(
        "'power' must be greater than 'alpha' (%.4g) when solving for sample size or minimum detectable effect; alpha is the power at zero effect",
        alpha
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Normal-approximation power for a standardized absolute effect
#' @keywords internal
#' @noRd
.normal_power <- function(effect, se, alpha, alternative) {
  if (se == 0) {
    return(if (effect == 0) alpha else 1)
  }
  z_a <- .z_alpha(alpha, alternative)
  z <- abs(effect) / se
  out <- pnorm(z - z_a)
  if (alternative == "two.sided") {
    out <- out + pnorm(-z - z_a)
  }
  min(out, 1)
}

#' Refuse a minimum detectable effect that a census does not have
#'
#' Enumerating the population leaves no sampling variance, so power is
#' `alpha` at exactly zero effect and 1 at every nonzero effect. No effect
#' has any intermediate power, so the requested target has no solution and
#' reporting zero alongside it would state a power the design never achieves.
#' @keywords internal
#' @noRd
.stop_census_mde <- function() {
  stop(
    "no minimum detectable effect exists: the design enumerates its population, so the sampling variance is zero and power is 'alpha' at zero effect and 1 at any nonzero effect",
    call. = FALSE
  )
}

#' Invert normal-approximation power for an absolute effect
#' @keywords internal
#' @noRd
.solve_normal_mde <- function(power, se, alpha, alternative, tol = 1e-10) {
  if (se == 0) {
    .stop_census_mde()
  }
  if (alternative == "one.sided") {
    return((.z_alpha(alpha, alternative) + qnorm(power)) * se)
  }
  hi <- (.z_alpha(alpha, alternative) + qnorm(power)) * se
  hi <- max(hi, se * sqrt(.Machine$double.eps))
  while (.normal_power(hi, se, alpha, alternative) < power) {
    hi <- hi * 2
  }
  uniroot(
    function(effect) .normal_power(effect, se, alpha, alternative) - power,
    c(0, hi),
    tol = tol
  )$root
}

#' Per-group FPC factor (returns 1 for Inf, 0-clamped)
#' @keywords internal
#' @noRd
.fpc_factor <- function(n, N) {
  if (is.infinite(N)) {
    return(1)
  }
  max(0, 1 - n / N)
}

#' Finite population correction on N - 1 degrees of freedom
#'
#' The transformed proportion scales carry their dispersion through
#' \eqn{Var(\hat p) = S^2(1/n - 1/N)} with the Bernoulli population variance
#' \eqn{S^2 = Np(1-p)/(N-1)}, so the factor that multiplies the per-unit term
#' is \eqn{(N-n)/(N-1)}, not \eqn{1 - n/N}. The delta-method variances on the
#' arcsine and log-odds scales inherit it unchanged, because both divide
#' \eqn{Var(\hat p)} by a function of \eqn{p} alone. Matches
#' `.bernoulli_var()` and `.prec_engine_prop()`, and reduces to
#' `.fpc_factor()` as \eqn{N} grows.
#' @keywords internal
#' @noRd
.fpc_factor_prop <- function(n, N) {
  if (is.infinite(N)) {
    return(1)
  }
  max(0, (N - n) / (N - 1))
}

#' Smallest proportion a design measures to a given CV
#'
#' The relative standard error se(p)/p decreases in p, so the CV target has
#' one root over the range that matters and the equality solution is also the
#' threshold: every larger p meets the target. `se` is the sampling standard
#' error under every interval method, so the root is the same under all four
#' and inverts in closed form from se^2 = p q / n_eff. The interval method
#' governs `moe` rather than `se`, so it does not enter here and the
#' signature does not carry it.
#' @keywords internal
#' @noRd
.prec_solve_prop <- function(cv, n, N, deff, resp_rate) {
  n_net <- n * resp_rate
  if (!is.infinite(N) && n_net >= N) {
    stop(
      "a census carries no sampling variance, so every proportion meets any 'cv'; supply 'p' instead",
      call. = FALSE
    )
  }

  # Every positive cv has a root in (0, 1) here, so there is no unattainable
  # case to report: p -> 1 as cv -> 0 and p -> 0 as cv grows.
  1 / (1 + .effective_from_n(n_net, N, deff) * cv^2)
}

#' Smallest proportion a design measures to a given relative margin of error
#'
#' The counterpart of `.prec_solve_prop()` on the interval row of the
#' package's precision 2x2. `moe(p) / p` reads the interval the chosen
#' method builds rather than the sampling variance, so unlike the `cv`
#' solve this one is genuinely method-specific and has no closed form.
#'
#' Two properties the `cv` solve does not share govern the search:
#'
#' - Wilson and Korn-Graubard half-widths do not vanish as `p` approaches
#'   1, so their relative margin of error has a positive floor and a target
#'   below it is unattainable rather than merely extreme.
#' - The back-transformed log-odds half-width turns upward once the logit
#'   spread outgrows `logit(p)`, near `p = 0.999` at `n = 1500` but as low
#'   as `p = 0.97` at `n = 30`. `moe(p) / p` is therefore decreasing then
#'   increasing, and the bracket ends at that turn: the root below it is
#'   the smallest proportion meeting the target, which is what this solve
#'   reports.
#' @keywords internal
#' @noRd
.prec_solve_prop_rmoe <- function(rmoe, n, alpha, N, deff, resp_rate, method,
                                  df = NULL) {
  n_net <- n * resp_rate
  if (!is.infinite(N) && n_net >= N) {
    stop(
      "a census carries no sampling variance, so every proportion meets any 'rmoe'; supply 'p' instead",
      call. = FALSE
    )
  }
  relative <- function(p) {
    .prec_engine_prop(p, n, alpha, N, deff, resp_rate, method, df)$moe / p
  }
  lower <- 1e-12
  cap <- 1 - 1e-9
  ## The turn sits near p = 1 - O(1 / n), so it is searched on the logit,
  ## which resolves 1 - p there.
  upper <- if (method == "logodds") {
    stats::plogis(stats::optimize(
      function(u) relative(stats::plogis(u)), stats::qlogis(c(lower, cap)),
      tol = 1e-10
    )$minimum)
  } else {
    cap
  }
  attainable <- relative(upper)
  if (rmoe <= attainable) {
    stop(
      sprintf(
        "'rmoe' = %.4g is unattainable: at n = %.4g the '%s' interval's relative margin of error never falls below %.4g",
        rmoe, n, method, attainable
      ),
      call. = FALSE
    )
  }
  if (relative(lower) <= rmoe) {
    return(lower)
  }
  uniroot(
    function(p) relative(p) - rmoe, c(lower, upper),
    tol = .Machine$double.eps^0.75
  )$root
}

#' Shared precision engine for proportions
#'
#' One documented variance equation: n_net = n * resp_rate responding
#' units drive both the leading term and the FPC's sampling fraction, and
#' deff multiplies the SRSWOR variance at n_net, so
#' se^2 = deff * p * q * fpc(n_net) / n_net with fpc(n) = (N - n)/(N - 1).
#' All four methods read that variance through one quantity, the effective
#' size n_eff = n_net / (deff * fpc) at which an infinite-population SRS
#' would reproduce it. A census drives fpc to 0, n_eff to infinity, and
#' every method's margin of error to 0.
#'
#' `se` is that sampling standard error and `cv` is se/p, so both are the
#' same under all four methods: they describe the estimator, not the
#' interval drawn around it. The methods differ only in `moe`, the half-width
#' of the interval they construct at that variance. Deriving `se` as moe/q
#' instead would make the reported standard error move when only the interval
#' construction changed.
#' @keywords internal
#' @noRd
.prec_engine_prop <- function(p, n, alpha, N, deff, resp_rate, method,
                              df = NULL) {
  q <- .q_alpha(alpha, df)
  q_complement <- 1 - p
  n_net <- n * resp_rate
  fpc <- if (is.infinite(N)) 1 else (N - n_net) / (N - 1)
  fpc <- .clamp_fpc(fpc, n_net, N)
  n_eff <- .effective_from_n(n_net, N, deff)
  # A census has no sampling variance under any method. Returning here keeps
  # the single census warning raised by .clamp_fpc() above.
  if (is.infinite(n_eff)) {
    return(list(se = 0, moe = 0, cv = 0))
  }

  se <- sqrt(p * q_complement / n_eff)

  moe <- if (method == "wald") {
    q * se
  } else if (method == "wilson") {
    .wilson_moe(p, n_eff, q)
  } else if (method == "logodds") {
    .logodds_moe(p, n_net, alpha, N, deff, df)
  } else {
    .beta_moe(p, .kg_effective(n_eff, n_net, alpha, df), alpha)
  }

  list(se = se, moe = moe, cv = se / p)
}

#' Convert a net sample size to its infinite-population equivalent
#'
#' `n_eff` is the size at which an infinite-population SRS has the same
#' variance as `n_net` responding units under `deff` and the finite
#' population correction, so `p q / n_eff == deff p q fpc / n_net`. It is
#' infinite for a census, where the design carries no sampling variance.
#' @keywords internal
#' @noRd
.effective_from_n <- function(n_net, N, deff) {
  if (is.infinite(N)) {
    return(n_net / deff)
  }
  if (n_net >= N) {
    return(Inf)
  }
  n_net * (N - 1) / (deff * (N - n_net))
}

#' Invert `.effective_from_n()`
#'
#' Solves `n_eff = n_net (N - 1) / (deff (N - n_net))` for `n_net`. The
#' result approaches but never reaches `N`, so a finite frame can always
#' meet a positive margin of error.
#' @keywords internal
#' @noRd
.n_from_effective <- function(n_eff, N, deff) {
  if (is.infinite(N)) {
    return(n_eff * deff)
  }
  n_eff * deff * N / (N - 1 + n_eff * deff)
}

#' Half-width of the Wilson score interval
#'
#' Centered at `(n p + z^2/2) / (n + z^2)` with half-width
#' `z sqrt(n p q + z^2/4) / (n + z^2)`, evaluated at the effective sample
#' size `n_eff` so that the design effect and the finite population
#' correction enter through the same variance the other methods use. The
#' half-width tends to 1/2 as `n_eff` tends to 0 and to 0 as it grows.
#' @keywords internal
#' @noRd
.wilson_moe <- function(p, n_eff, z) {
  z * sqrt(p * (1 - p) / n_eff + z^2 / (4 * n_eff^2)) / (1 + z^2 / n_eff)
}

#' Degrees-of-freedom adjusted effective sample size (Korn-Graubard 2.2)
#'
#' A variance estimated from `df` degrees of freedom is less stable than one
#' from a simple random sample of `n_net` units, and the Korn-Graubard
#' interval widens for it by scaling the effective size by
#' `(t_{n_net-1}(1 - alpha/2) / t_df(1 - alpha/2))^2`. Leaving `df` `NULL`
#' means "as stable as a simple random sample of the same size" and applies
#' no adjustment.
#' @keywords internal
#' @noRd
.kg_effective <- function(n_eff, n_net, alpha, df = NULL) {
  if (is.null(df) || is.infinite(n_eff)) {
    return(n_eff)
  }
  check_df(df)
  reference <- max(n_net - 1, 1)
  n_eff * (stats::qt(1 - alpha / 2, reference) /
             stats::qt(1 - alpha / 2, df))^2
}

#' Korn-Graubard (1998) confidence limits for a proportion
#'
#' The Clopper-Pearson limits evaluated at the effective sample size, with
#' `x = p * n_eff` playing the role of the observed positive count. Written
#' in the equivalent beta form of Korn and Graubard's equation (1.2).
#' @keywords internal
#' @noRd
.beta_limits <- function(p, n_eff, alpha) {
  if (is.infinite(n_eff)) {
    return(c(p, p))
  }
  x <- p * n_eff
  c(
    stats::qbeta(alpha / 2, x, n_eff - x + 1),
    stats::qbeta(1 - alpha / 2, x + 1, n_eff - x)
  )
}

#' Half-width of the Korn-Graubard interval
#'
#' The interval is asymmetric about `p`, so this half-width is a summary of
#' its length rather than an offset either limit sits at. Use [confint()] for
#' the limits themselves. It decreases monotonically in `n_eff`, which is
#' what lets `n_prop()` invert it.
#' @keywords internal
#' @noRd
.beta_moe <- function(p, n_eff, alpha) {
  limits <- .beta_limits(p, n_eff, alpha)
  (limits[2L] - limits[1L]) / 2
}

#' Half-width of the back-transformed log-odds interval
#'
#' Builds a symmetric interval on the logit scale from the delta-method
#' standard error of `logit(p_hat)`, then maps both endpoints back to the
#' probability scale and halves their distance. The variance of `p_hat`
#' follows the package convention shared with the Wald method,
#' `deff * N / (N - 1) * p q (1/n_net - 1/N)`, so the two agree as the
#' margin of error shrinks.
#' @keywords internal
#' @noRd
.logodds_moe <- function(p, n_net, alpha, N, deff = 1, df = NULL) {
  if (!is.infinite(N) && n_net >= N) {
    warning("net sample size >= population size; moe is 0", call. = FALSE)
    return(0)
  }
  bernoulli <- if (is.infinite(N)) 1 else N / (N - 1)
  fraction <- if (is.infinite(N)) 0 else 1 / N
  var_p <- deff * bernoulli * p * (1 - p) * (1 / n_net - fraction)
  spread <- .q_alpha(alpha, df) * sqrt(var_p) / (p * (1 - p))
  center <- qlogis(p)
  (plogis(center + spread) - plogis(center - spread)) / 2
}

#' Shared precision engine for means
#'
#' Same convention as .prec_engine_prop(), with fpc(n) = 1 - n / N.
#' Smallest mean a design measures to a given CV
#'
#' The standard error of a mean does not involve the mean, so the CV target
#' inverts directly. Only the magnitude is recoverable, `cv` being defined
#' against `abs(mu)`, and the positive root is returned.
#' @keywords internal
#' @noRd
.prec_solve_mean <- function(cv, var, n, alpha, N, deff, resp_rate,
                             df = NULL) {
  se <- .prec_engine_mean(var, NULL, n, alpha, N, deff, resp_rate, df)$se
  if (se == 0) {
    stop(
      "a census carries no sampling variance, so every mean meets any 'cv'; supply 'mu' instead",
      call. = FALSE
    )
  }
  se / cv
}

#' @keywords internal
#' @noRd
.prec_engine_mean <- function(var, mu, n, alpha, N, deff, resp_rate,
                              df = NULL) {
  q <- .q_alpha(alpha, df)
  n_net <- n * resp_rate
  n_eff <- n_net / deff
  fpc <- if (is.infinite(N)) 1 else 1 - n_net / N
  fpc <- .clamp_fpc(fpc, n_net, N)
  se <- sqrt(var * fpc / n_eff)
  list(se = se, moe = q * se, cv = .cv_of(se, mu))
}

#' Collision-free key for domain grouping and matching
#'
#' Length-prefixed encoding of each column value, so pasted keys are
#' injective regardless of separators inside the values. Keys built from
#' different data frames are comparable. Missing domain values error.
#' @keywords internal
#' @noRd
.domain_key <- function(df, cols) {
  vals <- lapply(cols, function(col) {
    v <- df[[col]]
    if (anyNA(v)) {
      stop(
        sprintf("domain column '%s' must not contain missing values", col),
        call. = FALSE
      )
    }
    v <- as.character(v)
    paste0(nchar(v), "_", v)
  })
  do.call(paste, c(vals, list(sep = ":")))
}

#' Error when a supplied gross sample exceeds a finite frame
#'
#' Vectorized over groups and indicators. `label` names the offender.
#' @keywords internal
#' @noRd
.check_gross_n <- function(n, N, label = NULL) {
  len <- max(length(n), length(N))
  n <- rep_len(n, len)
  N <- rep_len(N, len)
  bad <- is.finite(N) & n > N * (1 + 1e-9)
  if (any(bad)) {
    i <- which(bad)[1L]
    where <- if (!is.null(label)) {
      sprintf(" for %s", rep_len(label, len)[i])
    } else {
      ""
    }
    stop(
      sprintf(
        "'n' exceeds the population%s: cannot draw %s units without replacement from N = %s",
        where, round(n[i], 1), N[i]
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' A gross size within rounding of a finite frame is its census
#'
#' A closed-form size reaches `N` only as the target tends to zero, so a
#' tiny target leaves `n` a few rounding steps below `N`. The precision
#' functions then read `N - n` from two nearly equal numbers, lose its
#' digits, and the size misses its own target by up to 13 percent. Within a
#' relative 1e-7 of `N` the gap is a small fraction of a unit, and the
#' census, which meets any target, is the size.
#' @keywords internal
#' @noRd
.census_size <- function(n, N) {
  near <- is.finite(N) & n >= N * (1 - 1e-7) & n <= N * (1 + 1e-9)
  ifelse(near, N, n)
}

#' Error when a gross sample cannot be drawn from a finite frame
#' @keywords internal
#' @noRd
.check_attainable <- function(n, N, resp_rate, eligibility_rate = 1,
                              eligible_per_unit = 1) {
  if (!is.infinite(N) && n > N * (1 + 1e-9)) {
    stop(
      sprintf(
        "required sample (%s units drawn) exceeds the population (N = %s): even a census yields about %s respondents at resp_rate = %s%s%s; the target is unattainable",
        round(n, 1), N,
        round(N * resp_rate * eligibility_rate * eligible_per_unit), resp_rate,
        if (eligibility_rate == 1) "" else
          sprintf(" and eligibility_rate = %s", eligibility_rate),
        if (eligible_per_unit == 1) "" else
          sprintf(" with eligible_per_unit = %s", eligible_per_unit)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Error when a multi-indicator size exceeds a row's finite population
#'
#' The single-indicator engines refuse a target a census cannot meet, but
#' the indicator table raises sizes after them, by a floor or by another
#' row, and a census would then read as cv 0.
#' @keywords internal
#' @noRd
.check_multi_attainable <- function(n, indicators, where = "") {
  N <- indicators$N
  bad <- is.finite(N) & n > N * (1 + 1e-9)
  if (!any(bad)) {
    return(invisible(TRUE))
  }
  i <- which(bad)[1L]
  label <- if ("name" %in% names(indicators)) {
    sQuote(indicators$name[i])
  } else {
    i
  }
  stop(
    sprintf(
      "required sample (%s units drawn%s) exceeds the population of indicator %s (N = %s). A census cannot meet it, so lower the floor or the target",
      round(n, 1), where, label, N[i]
    ),
    call. = FALSE
  )
}

#' Expected number of positive cases a proportion design will yield
#'
#' A count, not a variance: `deff` has no part in it, since a design effect
#' says how precisely the proportion is estimated and nothing about how
#' many cases turn up. The sample is gross, so the response rate nets it
#' down first. `NULL` for anything but a proportion, which is what makes
#' the field conditional on the constructors.
#' @keywords internal
#' @noRd
.expected_cases <- function(n, params, type) {
  if (!identical(type, "proportion") || is.null(params$p) || is.null(n)) {
    return(NULL)
  }
  n * (params$resp_rate %||% 1) * .params_yield(params) * params$p
}

#' Sample size the minimum expected count of positive cases demands
#'
#' Gross units, on the same footing as every other size the package
#' returns: the count is realized among respondents, so the drawn sample
#' carries the `1 / resp_rate` inflation.
#' @keywords internal
#' @noRd
.n_from_min_cases <- function(min_cases, p, resp_rate) {
  check_scalar(min_cases, "min_cases")
  min_cases / (p * resp_rate)
}

#' Smallest issued size whose respondent count clears a target with
#' probability at least `level`
#'
#' Planning at the expected respondent count leaves about half the designs
#' short, so an assurance level converts a required number of respondents
#' into an issued number. Shared by `n_twophase()` and `n_panel()`, whose
#' assurance vocabulary is one thing.
#'
#' The tail is non-decreasing in the issued size, so the answer is bracketed
#' and bisected rather than walked. **Walking up from the expected count
#' returns a size that is not the smallest**, in two directions: below a level
#' of a half the answer lies under the expected count, and at a high response
#' rate the expected count itself overshoots, `need = 2` at a rate of 0.99
#' clearing 0.8 assurance with 2 issued where `ceiling(2 / 0.99)` is 3. The
#' bracket is anchored at `need - 1`, which clears no level at all, fewer
#' issued than needed carrying no chance of reaching it. The expected count is
#' the upper end of the first bracket, which is a starting point rather than a
#' contract: correctness rests on that anchor and on the monotonicity, so any
#' other start converges on the same answer and only the count of
#' `pbinom()` calls changes.
#' @keywords internal
#' @noRd
.assure_size <- function(m, r, level) {
  r <- rep_len(r, length(m))
  vapply(seq_along(m), function(i) {
    need <- .whole_count(m[i])
    if (need <= 0) {
      return(0)
    }
    clears <- function(g) {
      stats::pbinom(need - 1L, g, r[i], lower.tail = FALSE) >= level
    }
    lo <- need - 1L
    # A rate is at most 1, so this is at least `need` and no separate floor is
    # needed.
    hi <- ceiling(need / r[i])
    step <- max(1, hi - lo)
    while (!clears(hi)) {
      lo <- hi
      hi <- hi + step
      step <- step * 2
      if (hi > 1e9) {
        stop("the assurance search did not converge", call. = FALSE)
      }
    }
    while (hi - lo > 1) {
      mid <- lo + (hi - lo) %/% 2
      if (clears(mid)) hi <- mid else lo <- mid
    }
    as.double(hi)
  }, numeric(1))
}

#' Validate computed variance term before sqrt
#' @keywords internal
#' @noRd
.safe_variance <- function(V, what = "variance") {
  if (!is.finite(V)) {
    stop(
      sprintf("%s is not finite; check input parameters", what),
      call. = FALSE
    )
  }
  tol <- sqrt(.Machine$double.eps)
  if (V < -tol) {
    stop(
      sprintf("%s is negative; reduce overlap or overlap_cor, or adjust inputs", what),
      call. = FALSE
    )
  }
  if (V < 0) 0 else V
}


#' A key no two different value tuples share
#' @keywords internal
#' @noRd
.collision_free_key <- function(...) {
  vals <- list(...)
  enc <- lapply(vals, function(x) {
    x <- as.character(x)
    if (anyNA(x)) stop("join keys must not contain missing values", call. = FALSE)
    paste0(nchar(x), "_", x)
  })
  do.call(paste, c(enc, list(sep = ":")))
}
