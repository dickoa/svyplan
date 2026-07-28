#' Plot svyplan objects
#'
#' Visualize sampling fractions per stratum or power curves from svyplan
#' results.
#'
#' @param x A svyplan object.
#' @param npoints Number of points in the grid: the power curve (default
#'   101) or the budget frontier (default 25, since each point is a solve).
#' @param newdata Optional one-column data frame of `budget` values for
#'   `plot.svyplan_n()`. The default sweeps from the cheapest design that
#'   meets the hard targets up to twice the fitted budget.
#' @param ... Additional graphical parameters passed to [barplot()]
#'   (for strata) or [plot()] (for power and the budget frontier). These
#'   override the defaults, so you can set `main`, `col`, `ylab`, `xlab`,
#'   `ylim`, etc.
#'
#' @return `x`, invisibly.
#'
#' @details
#' `plot.svyplan_strata()` draws a bar chart of per-stratum sampling
#' fractions (`f = n / N`) using [barplot()]. This shows how
#' intensively each stratum is sampled, under Neyman allocation,
#' high-variance strata get higher fractions. A dashed horizontal line
#' marks the overall sampling fraction (`n / N`). Defaults:
#' `col = "grey40"`, `ylab = "Sampling fraction (f)"`, `las = 2`.
#'
#' `plot.svyplan_power()` draws the power-vs-sample-size curve using
#' [plot()]. The solved point is shown as a filled dot, with dashed
#' reference lines at the computed power and sample size, and a dotted
#' line at the significance level. Defaults: `ylim = c(0, 1)`,
#' `type = "l"`, `xlab = "Sample size (per group)"`, `ylab = "Power"`.
#'
#' `plot.svyplan_n()` draws the budget frontier for a fixed-budget joint
#' allocation ([n_alloc()] with `objective` and `budget`): what precision
#' each budget buys on the objective indicator, over the range where the
#' hard targets remain fundable. The fitted design is a filled dot. The
#' curve is the same one [predict()] returns as a table, so read exact
#' numbers there. Its shape is the point: the objective falls as
#' `1 / cost`, so the marginal return on budget flattens, and the plot
#' shows where. Other `svyplan_n` results have no frontier to draw and
#' produce an error naming what is plottable.
#'
#' @examples
#' # Sampling fraction per stratum
#' set.seed(1907)
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
#' # Budget frontier: what each budget buys on the objective indicator
#' frame <- data.frame(
#'   stratum = c("A", "B", "C"),
#'   N = c(4000, 3000, 3000),
#'   unit_cost = c(1, 1.2, 1.5)
#' )
#' measures <- data.frame(
#'   stratum = rep(frame$stratum, 2),
#'   name = rep(c("vaccination", "income"), each = 3),
#'   p = c(0.5, 0.4, 0.6, rep(NA, 3)),
#'   mean = c(rep(NA, 3), 50, 55, 60),
#'   sd = c(rep(NA, 3), 10, 12, 15)
#' )
#' targets <- data.frame(name = "vaccination", cv = 0.05)
#' fit <- n_alloc(frame, measures = measures, targets = targets,
#'                objective = "income", budget = 4000)
#' plot(fit)
#'
#' @seealso [predict.svyplan] for the sensitivity grids these curves are
#'   drawn from, and [strata_bound()], [power_prop()], [n_alloc()] for the
#'   results that are plottable.
#'
#' @name plot.svyplan
NULL

#' @rdname plot.svyplan
#' @export
plot.svyplan_strata <- function(x, ...) {
  strata <- x$strata
  if (!is.null(strata$lower) && !is.null(strata$upper)) {
    labels <- paste0(
      "[",
      signif(strata$lower, 3),
      ", ",
      signif(strata$upper, 3),
      ")"
    )
    labels[length(labels)] <- sub("\\)$", "]", labels[length(labels)])
  } else {
    labels <- paste0("H", strata$stratum)
  }

  f_h <- strata$n / strata$N
  names(f_h) <- labels

  method_label <- switch(
    x$method,
    cumrootf = "Dalenius-Hodges",
    geo = "Geometric",
    lh = "LH-inspired coordinate search",
    kozak = "Kozak-inspired local search",
    x$method
  )

  defaults <- list(
    col = "grey40",
    main = sprintf("%s (%d strata, cv = %.4f)", method_label, x$n_strata, x$cv),
    ylab = "Sampling fraction (f)",
    las = 2
  )
  args <- modifyList(defaults, list(...))
  do.call(barplot, c(list(height = f_h), args))

  f_overall <- sum(strata$n) / sum(strata$N)
  abline(h = f_overall, lty = 2, col = "grey50")

  invisible(x)
}

#' @rdname plot.svyplan
#' @export
plot.svyplan_power <- function(x, npoints = 101L, ...) {
  if (length(x$n) == 2L) {
    stop(
      "plot() does not support power objects with unequal-group n",
      call. = FALSE
    )
  }
  n_lo <- max(10, x$n * 0.1)
  n_hi <- x$n * 3
  n_seq <- seq(n_lo, n_hi, length.out = npoints)

  pw_seq <- .power_curve_grid(x, n_seq)

  type_label <- switch(
    x$type,
    proportion = "proportions",
    mean = "means",
    did_prop = "DiD proportions",
    did_mean = "DiD means",
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

#' @rdname plot.svyplan
#' @export
plot.svyplan_n <- function(x, npoints = 25L, newdata = NULL, ...) {
  if (!identical(x$params$mode, "budget_objective")) {
    stop(
      paste0(
        "plot() is available for fixed-budget joint allocations ",
        "(n_alloc() with 'objective' and 'budget'), which have a budget ",
        "frontier to draw. Use plot() on a strata_bound() or power_*() ",
        "result, or predict() for a sensitivity table on this one."
      ),
      call. = FALSE
    )
  }
  if (!is.numeric(npoints) || length(npoints) != 1L || is.na(npoints) ||
      npoints < 2) {
    stop("'npoints' must be a single number >= 2", call. = FALSE)
  }

  if (is.null(newdata)) newdata <- data.frame(budget = .budget_grid(x, npoints))
  fr <- predict(x, newdata)
  fr <- fr[!is.na(fr$.feasible) & fr$.feasible, , drop = FALSE]
  if (nrow(fr) < 2L) {
    stop("too few feasible budgets to draw a frontier; supply 'newdata'",
         call. = FALSE)
  }

  defaults <- list(
    type = "l",
    xlab = "Cost",
    ylab = "Objective cv",
    main = .frontier_title(x$params$objective)
  )
  args <- modifyList(defaults, list(...))
  do.call(plot, c(list(x = fr$cost, y = fr$cv), args))

  fitted_cv <- sqrt(x$objective_value)
  fitted_cost <- x$params$achieved$cost
  abline(h = fitted_cv, lty = 2, col = "grey50")
  abline(v = fitted_cost, lty = 2, col = "grey50")
  points(fitted_cost, fitted_cv, pch = 19)

  invisible(x)
}

#' Title for the frontier plot
#'
#' `objective` is the normalized component table, not a name, so a single
#' component is named and several are counted rather than concatenated into
#' a title too wide for the device.
#' @keywords internal
#' @noRd
.frontier_title <- function(objective) {
  if (is.null(objective) || nrow(objective) == 0L) return("Budget frontier")
  if (nrow(objective) == 1L) {
    return(sprintf("Budget frontier for %s", objective$component[1L]))
  }
  sprintf("Budget frontier (%d objective components)", nrow(objective))
}

#' Feasible budget sweep for the frontier plot
#'
#' Starts at the cheapest design meeting the hard targets rather than at an
#' arbitrary fraction of the fitted budget, so the grid contains no
#' infeasible points and the curve begins where the frontier actually does.
#' @keywords internal
#' @noRd
.budget_grid <- function(x, npoints) {
  p <- x$params
  targets <- if (is.null(p$targets) || nrow(p$targets) == 0L) NULL else
    p$targets
  lo <- NULL
  if (!is.null(targets)) {
    floor_fit <- tryCatch(
      n_alloc.default(
        frame = p$frame, measures = p$measures, targets = targets,
        unit_cost = p$unit_cost, alpha = p$alpha, deff = p$deff,
        resp_rate = p$resp_rate, min_n_stratum = p$min_n_stratum
      ),
      error = function(e) NULL
    )
    # the integer design is the true floor: the continuous optimum is not
    # itself affordable in whole units
    if (!is.null(floor_fit)) lo <- floor_fit$operational$cost
  }
  hi <- 2 * p$budget
  if (is.null(lo) || !is.finite(lo) || lo >= hi) lo <- 0.5 * p$budget
  seq(lo, hi, length.out = as.integer(npoints))
}

#' @keywords internal
#' @noRd
.power_curve_grid <- function(x, n_seq) {
  p <- x$params
  vapply(
    n_seq,
    function(ni) {
      tryCatch(
        {
          if (x$type == "proportion") {
            res <- power_prop(
              p1 = p$p1,
              p2 = p$p2,
              n = ni,
              power = NULL,
              alpha = p$alpha,
              N = p$N,
              deff = p$deff,
              resp_rate = p$resp_rate,
              alternative = p$alternative,
              overlap = p$overlap,
              overlap_cor = p$overlap_cor,
              method = p$method %||% "wald"
            )
          } else if (x$type %in% c("did_prop", "did_mean")) {
            res <- power_did(
              treat = p$treat,
              control = p$control,
              outcome = p$outcome,
              var = p$var,
              effect = x$effect,
              n = ni,
              power = NULL,
              alpha = p$alpha,
              N = p$N,
              deff = p$deff,
              resp_rate = p$resp_rate,
              alternative = p$alternative,
              overlap = p$overlap,
              overlap_cor = p$overlap_cor
            )
          } else {
            res <- power_mean(
              effect = x$effect,
              var = p$var,
              n = ni,
              power = NULL,
              alpha = p$alpha,
              N = p$N,
              deff = p$deff,
              resp_rate = p$resp_rate,
              alternative = p$alternative,
              overlap = p$overlap,
              overlap_cor = p$overlap_cor
            )
          }
          res$power
        },
        error = function(e) NA_real_
      )
    },
    numeric(1)
  )
}
