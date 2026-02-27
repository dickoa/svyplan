#' Plot svyplan objects
#'
#' Visualize sampling fractions per stratum or power curves from svyplan
#' results.
#'
#' @param x A svyplan object.
#' @param npoints Number of points in the power curve grid (default 101).
#' @param ... Additional graphical parameters passed to [barplot()]
#'   (for strata) or [plot()] (for power). These override the defaults,
#'   so you can set `main`, `col`, `ylab`, `xlab`, `ylim`, etc.
#'
#' @return `x`, invisibly.
#'
#' @details
#' `plot.svyplan_strata()` draws a bar chart of per-stratum sampling
#' fractions (`f_h = n_h / N_h`) using [barplot()]. This shows how
#' intensively each stratum is sampled — under Neyman allocation,
#' high-variance strata get higher fractions. A dashed horizontal line
#' marks the overall sampling fraction (`n / N`). Defaults:
#' `col = "grey40"`, `ylab = "Sampling fraction (f_h)"`, `las = 2`.
#'
#' `plot.svyplan_power()` draws the power-vs-sample-size curve using
#' [plot()]. The solved point is shown as a filled dot, with dashed
#' reference lines at the computed power and sample size, and a dotted
#' line at the significance level. Defaults: `ylim = c(0, 1)`,
#' `type = "l"`, `xlab = "Sample size (per group)"`,
#' `ylab = "Power"`.
#'
#' @examples
#' # Sampling fraction per stratum
#' set.seed(1)
#' sb <- strata_bound(rlnorm(2000, 6, 1), n_strata = 4, n = 200,
#'                     method = "cumrootf")
#' plot(sb)
#'
#' # Custom colour
#' plot(sb, col = "steelblue")
#'
#' # Power curve with defaults
#' pw <- power_prop(p1 = 0.30, p2 = 0.40, power = 0.80)
#' plot(pw)
#'
#' # Custom line width and colour
#' plot(pw, lwd = 2, col = "darkred")
#'
#' @name plot.svyplan
NULL

#' @rdname plot.svyplan
#' @export
plot.svyplan_strata <- function(x, ...) {
  strata <- x$strata
  labels <- paste0(
    "[", signif(strata$lower, 3), ", ",
    signif(strata$upper, 3), ")"
  )
  labels[length(labels)] <- sub("\\)$", "]", labels[length(labels)])

  f_h <- strata$n_h / strata$N_h
  names(f_h) <- labels

  method_label <- switch(x$method,
    cumrootf = "Dalenius-Hodges",
    geo      = "Geometric",
    lh       = "Lavallee-Hidiroglou",
    kozak    = "Kozak",
    x$method
  )

  defaults <- list(
    col  = "grey40",
    main = sprintf("%s (%d strata, cv = %.4f)", method_label, x$n_strata, x$cv),
    ylab = "Sampling fraction (f_h)",
    las  = 2
  )
  args <- modifyList(defaults, list(...))
  do.call(barplot, c(list(height = f_h), args))

  f_overall <- sum(strata$n_h) / sum(strata$N_h)
  abline(h = f_overall, lty = 2, col = "grey50")

  invisible(x)
}

#' @rdname plot.svyplan
#' @export
plot.svyplan_power <- function(x, npoints = 101L, ...) {
  n_lo <- max(10, x$n * 0.1)
  n_hi <- x$n * 3
  n_seq <- seq(n_lo, n_hi, length.out = npoints)

  pw_seq <- .power_curve_grid(x, n_seq)

  type_label <- switch(x$type,
    proportion = "proportions",
    mean       = "means",
    x$type
  )

  defaults <- list(
    type = "l",
    ylim = c(0, 1),
    xlab = "Sample size (per group)",
    ylab = "Power",
    main = sprintf("Power curve for %s", type_label)
  )
  args <- modifyList(defaults, list(...))
  do.call(plot, c(list(x = n_seq, y = pw_seq), args))

  abline(h = x$power, lty = 2, col = "grey50")
  abline(v = x$n, lty = 2, col = "grey50")
  abline(h = x$params$alpha, lty = 3, col = "grey70")
  points(x$n, x$power, pch = 19)

  invisible(x)
}

.power_curve_grid <- function(x, n_seq) {
  p <- x$params
  vapply(n_seq, function(ni) {
    tryCatch({
      if (x$type == "proportion") {
        res <- power_prop(
          p1 = p$p1, p2 = p$p2, n = ni, power = NULL,
          alpha = p$alpha, N = p$N, deff = p$deff,
          resp_rate = p$resp_rate,
          sides = p$sides, overlap = p$overlap, rho = p$rho
        )
      } else {
        res <- power_mean(
          delta = x$delta, var = p$var, n = ni, power = NULL,
          alpha = p$alpha, N = p$N, deff = p$deff,
          resp_rate = p$resp_rate,
          sides = p$sides, overlap = p$overlap, rho = p$rho
        )
      }
      res$power
    }, error = function(e) NA_real_)
  }, numeric(1))
}
