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
#' # Custom color
#' plot(sb, col = "steelblue")
#'
#' # Power curve with defaults
#' pw <- power_prop(p1 = 0.30, p2 = 0.40, power = 0.80)
#' plot(pw)
#'
#' # Custom line width and color
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

#' Chart a rotation schedule
#'
#' Draw the rotation chart of a [design_overlap()] schedule, one row per
#' cohort and one column per time period, with a labelled cell wherever
#' that cohort is in sample. This is the figure rotation designs are published as,
#' and it is the fastest way to see that a schedule is the one you meant.
#'
#' @param x A `svyplan_overlap` object from [design_overlap()].
#' @param type `"schedule"` (default) for the rotation chart, `"overlap"` for
#'   a bar chart of the overlap at each lag, which is what [print()] reports
#'   as a table.
#' @param start The launch to draw. `"gradual"` recruits one cohort per
#'   period, so the design fills up over a life. `"immediate"` adds, at the
#'   first period, one launch panel for every possible remaining life length.
#'   All panels begin at their first interview, so the design is full at once.
#'   [design_overlap()] gives the same mature overlap profile under either
#'   launch, although realized overlap is higher early in a gradual launch.
#'   Defaults to the launch `panel` was planned with, or to `"gradual"` when
#'   there is none. Available for a life without a break in it, for the reason
#'   [n_panel()]'s own `start` gives.
#' @param n_period Time periods to draw, with a cohort entering at each one.
#'   Defaults to the life plus four, which reaches the steady state and shows
#'   several periods of it.
#' @param panel Optional [n_panel()] result with `design = "rotating"`, whose
#'   `n_entrants` scales the sample labels and the total row from cohort
#'   shares to units. Its cohort count must match the schedule's life.
#' @param ... Additional graphical parameters. `main` and `col` are honored
#'   by both types, `col` being the cell fill for the chart and the bar fill
#'   for the profile; the profile passes the rest to [barplot()].
#'
#' @return `x`, invisibly.
#'
#' @details
#' Each row is a cohort with an entry period, and its cell at period
#' \eqn{t} is the stage of its life that period reaches, drawn when the
#' schedule puts that stage in sample. Cells are labelled by wave, counting
#' only the occasions in sample, so a schedule with a gap numbers its waves
#' consecutively across the gap, and a cohort launched part-way through a
#' life still starts at W1, its waves being counted from its own first
#' interview.
#'
#' **The chart is one launch, and the overlap is a steady state.** Drawing
#' one cohort entering per period leaves the early periods short of cohorts:
#' the total row climbs until every stage of the life is represented, which
#' happens at the period marked on the axis, and the overlaps
#' [design_overlap()] reports describe the design from that period on.
#'
#' That gradual start is a design decision rather than the only one, and
#' `start` draws either. `"immediate"` splits the first period into panels
#' planned for remaining life lengths from the full life down to one period.
#' Every panel begins at its first interview, so the design holds its whole
#' sample at once and the marked period is the first. [design_overlap()] gives
#' the mature overlap profile under either launch. Before a gradual launch
#' reaches the marked period, its realized overlap is higher because no full
#' set of cohorts has yet rotated through. [n_panel()] reports the response
#' and precision path while either launch settles.
#'
#' An immediate launch is drawn for a life without a break in it, for the
#' reason [n_panel()]'s own `start` gives. When `panel` is supplied and
#' `start` is not, the chart draws the launch that panel was planned with,
#' rather than defaulting past it.
#'
#' Lynn's printed Figure 5 labels Samples 1 through 10 only. Its total row of
#' 1,800 through period 10 nevertheless assumes that a new 300-unit sample
#' continues to enter in periods 6 through 10. This chart draws those implicit
#' Samples 11 through 15 as well, so every displayed total is supported by the
#' cohort rows above it.
#'
#' A take that varies over the life shades its cell in proportion, so a
#' schedule that subsamples later waves is visible as it is drawn.
#'
#' @references
#' Lynn, P. (2012). *Longitudinal Survey Methods for the Household Finance
#' and Consumption Survey*. Report to the European Central Bank. Figures 2 to
#' 5 are charts of this kind.
#'
#' @seealso [design_overlap()] for the schedule and the overlaps it produces,
#'   and [n_panel()] for the recruitment that fills it.
#'
#' @examples
#' # Two occasions in, two out, two in
#' plot(design_overlap("1-1-0-0-1-1"))
#'
#' # The overlap the chart produces, at each lag
#' plot(design_overlap("1-1-0-0-1-1"), type = "overlap")
#'
#' # CPS 4-8-4, one row per monthly cohort
#' plot(design_overlap("4-8-4"))
#'
#' # The same mature design brought up at once instead, full from period 1
#' plot(design_overlap("6"), start = "immediate")
#'
#' # Lynn Figure 5: six 300-unit launch samples, then 300 new units per period
#' lynn_target <- prec_prop(n = 1800, p = 0.5)
#' lynn_panel <- prec_panel(
#'   300, target = lynn_target, retention = rep(1, 5),
#'   design = "rotating", start = "immediate"
#' )
#' plot(design_overlap("6"), panel = lynn_panel, n_period = 10,
#'      main = "1-1-1-1-1-1 rotating panel: immediate start")
#'
#' # Labelled in units rather than cohort shares
#' target <- n_prop(p = 0.5, moe = 0.031)
#' rot <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
#'                resp_rate = 0.728, design = "rotating")
#' plot(design_overlap("5"), panel = rot)
#'
#' @export
plot.svyplan_overlap <- function(x, type = c("schedule", "overlap"),
                                 start = c("gradual", "immediate"),
                                 n_period = NULL, panel = NULL, ...) {
  type <- match.arg(type)
  if (identical(type, "overlap")) {
    return(.plot_overlap_profile(x, ...))
  }
  # NULL here means "not stated", which a panel carrying its own launch is
  # then asked for; stating it explicitly overrides the panel, which is how
  # the two launches of one plan are compared
  .plot_rotation_chart(
    x, n_period = n_period, panel = panel,
    start = if (missing(start)) NULL else match.arg(start), ...
  )
}

#' Bar chart of the overlap at each lag
#' @keywords internal
#' @noRd
.plot_overlap_profile <- function(x, ...) {
  defaults <- list(
    col = "grey40",
    ylim = c(0, 1),
    xlab = "Lag (occasions)",
    ylab = "Issued-sample overlap",
    main = sprintf("Overlap by lag (%s)", .fmt_schedule(attr(x, "schedule")))
  )
  args <- modifyList(defaults, list(...))
  do.call(
    barplot,
    c(list(height = as.double(x), names.arg = names(x)), args)
  )
  invisible(x)
}

#' Rotation chart, cohorts down and time periods across
#'
#' Cohort `c` enters at period `c`, so period `t` shows it at stage
#' `t - c + 1`. Everything drawn follows from that one index; the schedule
#' decides only whether a cell is in sample and how dark it is.
#' @keywords internal
#' @noRd
.plot_rotation_chart <- function(x, n_period, panel, start, ...) {
  w <- attr(x, "schedule")
  life <- attr(x, "life")

  # a cohort enters at every period drawn, which is what holds the sample at
  # a steady state once the life is spanned; drawing a fixed set of cohorts
  # instead would wind the design down again at the right-hand edge
  n_period <- .whole_arg(n_period, "n_period", life + 4L)
  take <- .chart_take(panel, w)
  start <- .check_chart_start(start, w, panel)
  cohorts <- .chart_cohorts(w, start, n_period)
  n_sample <- length(cohorts)

  args <- modifyList(
    list(col = "grey75", main = sprintf("Rotation chart (%s)",
                                        .fmt_schedule(attr(x, "schedule")))),
    list(...)
  )

  op <- par(mar = c(2.6, if (is.null(take)) 6.5 else 9, 4.2, 1.2))
  on.exit(par(op), add = TRUE)

  plot.new()
  plot.window(xlim = c(0.5, n_period + 0.5), ylim = c(-1.2, n_sample + 0.5))

  cell_cex <- min(1, 0.82 / max(strwidth(paste0("W", sum(w > 0))), 1e-8))
  for (i in seq_along(cohorts)) {
    co <- cohorts[[i]]
    # a cohort numbers its own waves from its own first interview, so one
    # launched mid-life starts at W1 like any other
    wave <- cumsum(co$w > 0)
    for (t in seq_len(n_period)) {
      stage <- t - co$entry + 1L
      if (stage < 1L || stage > length(co$w) || co$w[[stage]] <= 0) next
      y <- n_sample - i + 1
      share <- co$w[[stage]] / max(w)
      rect(t - 0.42, y - 0.32, t + 0.42, y + 0.32,
           col = adjustcolor(args$col, alpha.f = 0.35 + 0.65 * share),
           border = "grey35")
      if (cell_cex >= 0.45) {
        text(t, y, paste0("W", wave[[stage]]), cex = cell_cex)
      }
    }
  }

  in_sample <- .chart_totals(cohorts, n_period)
  totals <- if (is.null(take)) {
    .drop_trailing_zeros(in_sample)
  } else {
    format(round(in_sample * take), trim = TRUE)
  }
  text(seq_len(n_period), -0.25, totals, cex = 0.75, col = "grey25")

  labels <- paste0("Sample ", seq_len(n_sample))
  if (!is.null(take)) {
    labels <- paste0(labels, "  (", format(take, trim = TRUE), ")")
  }
  axis(2, at = c(n_sample:1, -0.25), labels = c(labels, "Total"),
       las = 1, tick = FALSE, line = -0.6, cex.axis = 0.8)
  axis(3, at = seq_len(n_period), labels = seq_len(n_period),
       tick = FALSE, line = -0.9, cex.axis = 0.8)
  mtext("Time period", side = 3, line = 0.9, cex = 0.9)
  title(main = args$main, line = 2.4)

  # the overlap describes the design once every stage of the life is
  # represented, which a gradual launch reaches at period `life` and an
  # immediate one holds from the first period
  settled <- if (identical(start, "immediate")) 1L else life
  if (n_period >= settled) {
    abline(v = settled - 0.5, lty = 3, col = "grey55")
    text(max(settled - 0.35, 0.6), -0.95, "steady state from here", adj = 0,
         cex = 0.7, col = "grey40")
  }

  invisible(x)
}

#' The launch a chart draws, inherited from a panel when not stated
#'
#' An immediate launch has to cover every stage of the life, which for a
#' schedule with a gap means cohorts selected before they are first
#' interviewed. That is a longer definition than a launch label carries, and
#' [n_panel()] refuses it for the same reason, so the two surfaces agree on
#' what is expressible.
#' @keywords internal
#' @noRd
.check_chart_start <- function(start, w, panel = NULL) {
  if (is.null(start)) {
    # a panel that was planned with a launch has already answered this, and a
    # chart that ignored it would draw a design its own labels contradict
    start <- panel$start %||% "gradual"
  }
  start <- match.arg(start, c("gradual", "immediate"))
  if (identical(start, "immediate") && any(w <= 0)) {
    stop(
      "an immediate launch of a schedule with a gap needs cohorts selected before their first interview, which is not described yet; chart the gradual launch, or a schedule without a break",
      call. = FALSE
    )
  }
  start
}

#' The cohorts a launch puts on the chart, in drawing order
#'
#' A gradual launch is one cohort an occasion and nothing else. An immediate
#' launch adds, at the first occasion, one component for every possible
#' remaining life length, each at its first displayed stage, so the design
#' holds its whole sample at once. They are ordered like Lynn's Figure 5: the
#' one-occasion tail first, up through the full-life cohort.
#' @keywords internal
#' @noRd
.chart_cohorts <- function(w, start, n_period) {
  life <- length(w)
  entering <- lapply(seq_len(n_period), function(p) list(entry = p, w = w))
  if (identical(start, "gradual")) {
    return(entering)
  }
  launch <- lapply(seq.int(life, 1L), function(k) {
    list(entry = 1L, w = w[k:life])
  })
  c(launch, entering[-1L])
}

#' Units in sample at each drawn period, in cohort shares
#'
#' Summed over the cohorts drawn rather than derived from the schedule, so
#' one reading serves both launches. Under a gradual launch the total climbs
#' until the whole life is spanned and is `sum(w)` from period `life` on;
#' under an immediate one it is `sum(w)` throughout. That figure is the
#' `n_occasion` [design_overlap()] divides by, so where the total reaches it
#' is where the chart's overlaps become the design's.
#' @keywords internal
#' @noRd
.chart_totals <- function(cohorts, n_period) {
  vapply(seq_len(n_period), function(t) {
    sum(vapply(cohorts, function(co) {
      stage <- t - co$entry + 1L
      if (stage >= 1L && stage <= length(co$w)) co$w[[stage]] else 0
    }, numeric(1)))
  }, numeric(1))
}

#' Entrants per cohort for the chart labels, or NULL for cohort shares
#'
#' A panel may only scale a schedule it represents. [n_panel()] models equal
#' cohorts interviewed at every wave of their life, so its `n_in_sample` is
#' `n_cohorts` times the entrants; a schedule that leaves the sample and
#' returns, or that subsamples a later wave, holds fewer than that at an
#' occasion. Scaling one by the other would print two designs on one chart,
#' the label column reading from the panel and the total row from the
#' schedule.
#' @keywords internal
#' @noRd
.chart_take <- function(panel, w) {
  if (is.null(panel)) return(NULL)
  life <- length(w)
  if (!inherits(panel, "svyplan_panel")) {
    stop("'panel' must be an n_panel() or prec_panel() result", call. = FALSE)
  }
  if (!identical(panel$design, "rotating")) {
    stop(
      "'panel' is a fixed panel, which recruits one cohort and has no rotation to chart; pass a result with design = \"rotating\"",
      call. = FALSE
    )
  }
  if (!identical(as.integer(panel$n_cohorts), as.integer(life))) {
    stop(
      sprintf(
        "'panel' has %d live cohorts and the schedule has a life of %d occasions; they must describe the same design",
        as.integer(panel$n_cohorts), as.integer(life)
      ),
      call. = FALSE
    )
  }
  if (!isTRUE(all.equal(as.numeric(w), rep(1, life)))) {
    stop(
      sprintf(
        "'panel' recruits equal cohorts interviewed at every wave, so it holds %d cohorts in sample at an occasion, where this schedule holds %s; %s is not a design n_panel() describes, so its entrants cannot scale this chart",
        as.integer(life), format(sum(w), trim = TRUE),
        if (any(w <= 0)) {
          "a life with a break in it"
        } else {
          "a life that subsamples a later wave"
        }
      ),
      call. = FALSE
    )
  }
  # print.svyplan_panel() issues whole units and reads the recruitment
  # through .panel_recruit(); the chart has to agree with it, or the same
  # object would report two entrant counts
  ceiling(.panel_recruit(panel))
}

#' Whole-number plot argument with a default
#' @keywords internal
#' @noRd
.whole_arg <- function(value, name, default) {
  if (is.null(value)) return(as.integer(default))
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < 1 || value != trunc(value)) {
    stop(sprintf("'%s' must be a whole number >= 1", name), call. = FALSE)
  }
  as.integer(value)
}

#' Cohort-share totals read better without the trailing zeros of a share
#' @keywords internal
#' @noRd
.drop_trailing_zeros <- function(v) {
  out <- format(v, trim = TRUE, drop0trailing = TRUE)
  sub("^0$", "", out)
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
        resp_rate = p$resp_rate, min_n_stratum = p$min_n_stratum,
        fpc = p$fpc %||% "unit"
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
