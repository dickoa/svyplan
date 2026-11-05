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
#' numbers there. Its shape is the point, since the objective falls as
#' `1 / cost`, so the marginal return on budget flattens, and the plot
#' shows where. A fit with a `psu` register draws its field designs, their
#' cost and objective, as its grid reports them. Other `svyplan_n` results
#' have no frontier to draw and produce an error naming what is plottable.
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
    kozak = "Kozak random search",
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
  xy <- .frontier_xy(x, fr)
  do.call(plot, c(list(x = xy$cost, y = xy$cv), args))

  abline(h = xy$fitted[2], lty = 2, col = "grey50")
  abline(v = xy$fitted[1], lty = 2, col = "grey50")
  points(xy$fitted[1], xy$fitted[2], pch = 19)

  invisible(x)
}

#' Rotation chart of a rotation pattern
#'
#' Draw the rotation chart of a [design_rotation()] pattern, one row per
#' panel and one column per period, with a labeled cell wherever the panel
#' is interviewed. This is the figure a rotation design is published as in
#' a survey's methodology report, and the fastest way to see that a pattern
#' is the one you meant. The chart is drawn from `as.data.frame(x)`, the
#' table of cells described on [print.svyplan_rotation()], which can be
#' drawn another way, with ggplot2 for instance.
#'
#' @param x A `svyplan_rotation` object from [design_rotation()].
#' @param start How the first panels enter. `"gradual"` selects one panel
#'   per period, so the sample builds up over a life. `"immediate"` selects,
#'   in the first period, one panel for each possible remaining life, so the
#'   full sample is in the field from period 1. Defaults to the start
#'   `panel` was planned with, or to `"gradual"` when there is none.
#' @param n_period Periods to draw, a new panel entering at each one.
#'   Defaults to the life plus four.
#' @param panel Optional [n_panel()] or [prec_panel()] result with
#'   `design = "rotating"`. Its `n_entrants` labels the rows and the total
#'   in units rather than in panel shares, and the cell table gains a column
#'   `n`. It must have been planned on the rotation drawn.
#' @param label What each cell shows. `"wave"` (default) writes `W1`, `W2`,
#'   and so on, `"none"` leaves the cells blank, any other column of the
#'   cell table writes its values, and a function of the cell table returns
#'   one label per row.
#' @param fill What the cell color shows. `"none"` (default) draws every
#'   cell in one color. A column of the cell table, such as
#'   `"rotation_group"` or `"wave"`, or a function of the cell table
#'   returning one value per row, gives each distinct value its own color.
#' @param col Cell colors. One color when `fill = "none"`, otherwise at
#'   least one color per distinct fill value, taken in sorted order of the
#'   values. Defaults to grey, and to a qualitative palette when `fill` is
#'   set.
#' @param rows `"panel"` (default) draws one row per panel.
#'   `"rotation_group"` draws one row per rotation group, the panels that
#'   take the group over following one another along the row.
#' @param period_labels Optional labels for the periods, one per period
#'   drawn, such as `toupper(month.abb)`.
#' @param marker `TRUE` (default) marks the first period with the full
#'   sample in the field.
#' @param main Chart title.
#' @param ... Not used. Any other argument is an error.
#'
#' @return `x`, invisibly.
#'
#' @details
#' Each row is a panel, the set of units that enter and leave the sample
#' together (Steel and McLaren, 2008). A cell is drawn in each period the
#' panel is interviewed, and by default its label is the wave, the place of
#' that interview in the rotation pattern. That is the meaning of wave in
#' Regulation (EU) 2019/2240, "the sub-sample of observation units due to
#' be interviewed for the same nth time according to the sample rotation
#' scheme", and the convention of the CPS month in sample and of the
#' Statistics Canada rotation chart (Statistics Canada, 2017, Figure 2.1).
#'
#' The wave and the interview number differ only for the panels of an
#' immediate start. Under `"6"`, the panel selected in period 1 with two
#' periods of its life left is first interviewed then, at wave 5.
#' `label = "visit"` shows the interview number instead, which is how Lynn
#' (2012) labels the samples of an immediate start in Figure 5. Attrition
#' follows the interview number, and that is how [n_panel()] plans it.
#'
#' A rotation group is the slot a new panel takes over when the panel
#' before it rotates out (Steel and McLaren, 2008, and Statistics Canada,
#' 2017). Kalton (2009) and Hussmanns, Mehran and Verma (1990) also use the
#' term for a panel itself. Each new panel joins the lowest-numbered group
#' that is free at every one of its interviews. Under patterns such as
#' `"6"`, `"2-2-2"` and `"4-8-4"` every group then holds a panel in sample
#' at every period of a full sample. Under `"3-1-2"` and `"3-2-2"` no
#' assignment can do that, because no other panel can fill the period a
#' panel spends out of sample without overlapping it, so there are more
#' groups than panels in sample and each group has periods with none.
#'
#' The dotted line marks the first period with the full sample in the
#' field, every stage of the life held by one panel. A gradual start
#' reaches it once a whole life has passed, and an immediate start has it
#' from period 1. From that period on, the overlaps [design_overlap()]
#' reports are those of the chart. Under an immediate start the interview
#' mix still changes after that, because every launch panel begins at its
#' first interview. [n_panel()] reports the period at which it reaches
#' equilibrium (Lynn, 2012), the start-up period before it (Kalton, 2009)
#' being when response and precision still move.
#'
#' A launch panel whose remaining life starts out of sample is selected in
#' period 1 and first interviewed when its life returns to sample, as
#' [n_panel()] plans it.
#'
#' Lynn's Figure 5 labels Samples 1 through 10 only. Its total row of 1,800
#' through period 10 nevertheless assumes that a new 300-unit sample
#' continues to enter in periods 6 through 10. This chart draws those
#' implicit panels 11 through 15 as well, so every displayed total is
#' supported by the rows above it.
#'
#' A take that varies over the life shades its cell in proportion, so a
#' pattern that subsamples later waves is visible as it is drawn.
#'
#' @references
#' Commission Implementing Regulation (EU) 2019/2240 of 16 December 2019
#' specifying the technical items of the data set, establishing the
#' technical formats for transmission of information and specifying the
#' detailed arrangements and content of the quality reports on the
#' organisation of a sample survey in the labour force domain in accordance
#' with Regulation (EU) 2019/1700 of the European Parliament and of the
#' Council. Article 2.
#'
#' Hussmanns, R., Mehran, F., and Verma, V. (1990). *Surveys of Economically
#' Active Population, Employment, Unemployment and Underemployment: An ILO
#' Manual on Concepts and Methods*. International Labour Office, Geneva.
#' Chapter 11.
#'
#' Kalton, G. (2009). Designs for surveys over time. In D. Pfeffermann and
#' C. R. Rao (eds.), *Handbook of Statistics 29A, Sample Surveys: Design,
#' Methods and Applications*, 89--108. Elsevier.
#'
#' Lynn, P. (2012). *Longitudinal Survey Methods for the Household Finances
#' and Consumption Survey*. Report prepared for the European Central Bank.
#' Figures 2 to 5 are charts of this kind.
#'
#' Statistics Canada (2017). *Methodology of the Canadian Labour Force
#' Survey*. Catalogue no. 71-526-X.
#'
#' Steel, D. and McLaren, C. (2008). *Design and Analysis of Repeated
#' Surveys*. Working Paper 11-08, Centre for Statistical and Survey
#' Methodology, University of Wollongong.
#'
#' @seealso [print.svyplan_rotation()] for the cell table,
#'   [design_overlap()] for the overlaps the chart produces, and
#'   [n_panel()] for the recruitment that fills it.
#'
#' @examples
#' # Two waves in, two out, two in, gradual start
#' plot(design_rotation("1-1-0-0-1-1"))
#'
#' # Statistics Canada, Figure 2.1: six months in sample from January,
#' # colored by rotation group
#' plot(design_rotation("6"), start = "immediate", n_period = 12,
#'      fill = "rotation_group", period_labels = toupper(month.abb))
#'
#' # Rwanda LFS 2024, Annex 1: panel number and appearance letter
#' plot(design_rotation("2-2-2"), start = "immediate", n_period = 12,
#'      label = function(d) paste0(d$panel, LETTERS[d$wave]),
#'      fill = "wave")
#'
#' # The same design, one row per rotation group
#' plot(design_rotation("2-2-2"), start = "immediate", n_period = 12,
#'      rows = "rotation_group", label = "panel")
#'
#' # CPS 4-8-4, one row per monthly panel
#' plot(design_rotation("4-8-4"))
#'
#' # Lynn Figure 5: six 300-unit launch samples, then 300 new units per
#' # period, labeled by interview number as Lynn labels them
#' lynn_target <- prec_prop(n = 1800, p = 0.5)
#' lynn_panel <- prec_panel(
#'   300, target = lynn_target, retention = rep(1, 5),
#'   design = "rotating", start = "immediate"
#' )
#' plot(design_rotation("6"), panel = lynn_panel, n_period = 10,
#'      label = "visit", main = "1-1-1-1-1-1 rotating panel: immediate start")
#'
#' @export
plot.svyplan_rotation <- function(x, start = c("gradual", "immediate"),
                                  n_period = NULL, panel = NULL,
                                  label = "wave", fill = "none", col = NULL,
                                  rows = c("panel", "rotation_group"),
                                  period_labels = NULL, marker = TRUE,
                                  main = NULL, ...) {
  .check_unused_dots(...)
  rows <- match.arg(rows)
  cells <- as.data.frame(
    x, start = if (missing(start)) NULL else start,
    n_period = n_period, panel = panel
  )
  w <- as.double(x)
  n_period <- attr(cells, "n_period", exact = TRUE)
  entrants <- .chart_take(panel, w)
  lab <- .chart_labels(cells, label)
  fill_col <- .chart_fill(cells, fill, col)
  axis_labels <- .chart_period_labels(period_labels, n_period)
  if (!is.logical(marker) || length(marker) != 1L || is.na(marker)) {
    stop("'marker' must be TRUE or FALSE", call. = FALSE)
  }
  if (is.null(main)) {
    main <- sprintf("Rotation chart (%s)", .fmt_rotation(w))
  }

  row_id <- cells[[rows]]
  n_row <- max(row_id)
  row_labels <- if (identical(rows, "panel")) {
    out <- paste0("Panel ", seq_len(n_row))
    if (!is.null(entrants)) {
      out <- paste0(out, "  (", format(entrants, trim = TRUE), ")")
    }
    out
  } else {
    paste0("Rotation group ", seq_len(n_row))
  }

  op <- par(mar = c(2.6, max(6.5, 0.55 * max(nchar(row_labels)) + 1.5),
                    4.2, 1.2))
  on.exit(par(op), add = TRUE)
  plot.new()
  plot.window(xlim = c(0.5, n_period + 0.5), ylim = c(-1.2, n_row + 0.5))

  y <- n_row - row_id + 1
  share <- cells$take / max(w)
  rect(cells$period - 0.42, y - 0.32, cells$period + 0.42, y + 0.32,
       col = mapply(adjustcolor, fill_col, alpha.f = 0.35 + 0.65 * share,
                    USE.NAMES = FALSE),
       border = "grey35")
  if (!is.null(lab)) {
    cell_cex <- min(1, 0.82 / max(strwidth(lab), 1e-8))
    if (cell_cex >= 0.45) {
      text(cells$period, y, lab, cex = cell_cex)
    }
  }

  in_sample <- vapply(seq_len(n_period), function(t) {
    sum(cells$take[cells$period == t])
  }, numeric(1L))
  totals <- if (is.null(entrants)) {
    .drop_trailing_zeros(in_sample)
  } else {
    format(round(in_sample * entrants), trim = TRUE)
  }
  text(seq_len(n_period), -0.25, totals, cex = 0.75, col = "grey25")

  axis(2, at = c(n_row:1, -0.25), labels = c(row_labels, "Total"),
       las = 1, tick = FALSE, line = -0.6, cex.axis = 0.8)
  axis(3, at = seq_len(n_period), labels = axis_labels,
       tick = FALSE, line = -0.9, cex.axis = 0.8)
  mtext("Period", side = 3, line = 0.9, cex = 0.9)
  title(main = main, line = 2.4)

  full <- attr(cells, "full_sample_from", exact = TRUE)
  if (marker && n_period >= full) {
    abline(v = full - 0.5, lty = 3, col = "grey55")
    text(max(full - 0.35, 0.6), -0.95, "full sample from here", adj = 0,
         cex = 0.7, col = "grey40")
  }

  invisible(x)
}

#' Bar chart of the overlap at each lag
#'
#' Draw the issued-sample overlap of a [design_overlap()] result as a bar
#' per lag, the profile [print()] reports as a table. The rotation chart
#' the overlap comes from is drawn from the rotation, with
#' [plot.svyplan_rotation()].
#'
#' @param x A `svyplan_overlap` object from [design_overlap()].
#' @param ... Graphical parameters passed to [barplot()], which override the
#'   defaults: `col = "grey40"`, `ylim = c(0, 1)`, and axis labels and a
#'   title naming the rotation.
#'
#' @return `x`, invisibly.
#'
#' @seealso [design_overlap()], and [plot.svyplan_rotation()] for the
#'   rotation chart.
#'
#' @examples
#' plot(design_overlap("1-1-0-0-1-1"))
#'
#' # The CPS 4-8-4 profile, with its 50 percent overlap at 12 months
#' plot(design_overlap("4-8-4"), col = "steelblue")
#'
#' @export
plot.svyplan_overlap <- function(x, ...) {
  passed <- names(match.call(expand.dots = FALSE)$...)
  chart <- intersect(passed, c("type", "start", "n_period", "panel", "label",
                               "fill", "rows", "period_labels", "marker"))
  if (length(chart) > 0L) {
    stop(
      sprintf(
        "'%s' is an argument of the rotation chart, which is drawn from the rotation: plot(x$rotation, ...)",
        chart[[1L]]
      ),
      call. = FALSE
    )
  }
  defaults <- list(
    col = "grey40",
    ylim = c(0, 1),
    xlab = "Lag (occasions)",
    ylab = "Issued-sample overlap",
    main = sprintf("Overlap by lag (%s)",
                   .fmt_rotation(as.double(.overlap_rotation(x))))
  )
  args <- modifyList(defaults, list(...))
  do.call(
    barplot,
    c(list(height = as.double(x), names.arg = names(x)), args)
  )
  invisible(x)
}

#' The launch a chart draws, inherited from a panel when not stated
#' @keywords internal
#' @noRd
.check_chart_start <- function(start, panel = NULL) {
  if (is.null(start)) {
    start <- panel$start %||% "gradual"
  }
  match.arg(start, c("gradual", "immediate"))
}

#' The cohorts of a launch, in Lynn's Figure 5 order
#'
#' The one construction the chart, [n_panel()] and [design_schedule()] read.
#' A gradual launch is one cohort an occasion. An immediate one adds, at the
#' first occasion, a panel covering each stage of the life with the rest of
#' the life from it, the one-occasion tail first. A panel whose rest of life
#' starts out of sample is held and first interviewed later.
#' @keywords internal
#' @noRd
.rotation_cohorts <- function(w, start, n_period) {
  life <- length(w)
  entering <- lapply(seq_len(n_period), function(p) {
    list(entry = p, w = w)
  })
  if (identical(start, "gradual")) {
    return(entering)
  }
  launch <- lapply(seq.int(life, 1L), function(k) {
    list(entry = 1L, w = w[k:life])
  })
  c(launch, entering[-1L])
}

#' One row per interview the chart draws
#' @keywords internal
#' @noRd
.rotation_cells <- function(w, start, n_period, entrants = NULL) {
  life <- length(w)
  cohorts <- .rotation_cohorts(w, start, n_period)
  wave_at <- cumsum(w > 0)
  groups <- .rotation_groups(lapply(cohorts, function(co) {
    (co$entry + seq_along(co$w) - 1L)[co$w > 0]
  }))
  rows <- vector("list", length(cohorts))
  for (i in seq_along(cohorts)) {
    co <- cohorts[[i]]
    offset <- life - length(co$w)
    stage <- seq_along(co$w)
    period <- co$entry + stage - 1L
    on <- co$w > 0
    keep <- on & period <= n_period
    rows[[i]] <- data.frame(
      period = as.integer(period[keep]),
      panel = i,
      entry = as.integer(co$entry),
      launch = identical(start, "immediate") && i <= life,
      rotation_group = groups[[i]],
      wave = wave_at[offset + stage[keep]],
      visit = cumsum(on)[keep],
      life_stage = as.integer(offset + stage[keep]),
      take = co$w[keep]
    )
  }
  out <- do.call(rbind, rows)
  if (!is.null(entrants)) {
    out$n <- out$take * entrants
  }
  rownames(out) <- NULL
  out
}

#' The rotation group of each panel, panels taken in order
#'
#' A panel joins the lowest-numbered group free at every period it is
#' interviewed, over its whole life and not only a window, so widening a
#' window never moves a panel. The chart and [design_schedule()] both read
#' this, so their groups agree.
#' @keywords internal
#' @noRd
.rotation_groups <- function(occupancy) {
  busy <- list()
  out <- integer(length(occupancy))
  for (i in seq_along(occupancy)) {
    g <- 1L
    while (g <= length(busy) && any(occupancy[[i]] %in% busy[[g]])) {
      g <- g + 1L
    }
    busy[[g]] <- c(if (g <= length(busy)) busy[[g]], occupancy[[i]])
    out[[i]] <- g
  }
  out
}

#' The text of each cell
#' @keywords internal
#' @noRd
.chart_labels <- function(cells, label) {
  if (is.function(label)) {
    out <- label(cells)
    if (length(out) != nrow(cells)) {
      stop(
        sprintf(
          "'label' returned %d values for %d cells. A label function must return one per row of the cell table",
          length(out), nrow(cells)
        ),
        call. = FALSE
      )
    }
    return(as.character(out))
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop(.chart_arg_message("label", "\"wave\", \"none\""), call. = FALSE)
  }
  if (identical(label, "none")) {
    return(NULL)
  }
  if (identical(label, "wave")) {
    return(paste0("W", cells$wave))
  }
  if (!label %in% names(cells)) {
    stop(.chart_arg_message("label", "\"wave\", \"none\""), call. = FALSE)
  }
  format(cells[[label]], trim = TRUE, drop0trailing = TRUE)
}

#' The color of each cell, before take shading
#' @keywords internal
#' @noRd
.chart_fill <- function(cells, fill, col) {
  if (is.character(fill) && length(fill) == 1L && identical(fill, "none")) {
    col <- col %||% "grey75"
    if (length(col) != 1L) {
      stop("'col' must be one color when 'fill' is \"none\"", call. = FALSE)
    }
    return(rep(col, nrow(cells)))
  }
  key <- if (is.function(fill)) {
    fill(cells)
  } else if (is.character(fill) && length(fill) == 1L && !is.na(fill) &&
             fill %in% names(cells)) {
    cells[[fill]]
  } else {
    stop(.chart_arg_message("fill", "\"none\""), call. = FALSE)
  }
  if (length(key) != nrow(cells)) {
    stop(
      sprintf(
        "'fill' returned %d values for %d cells. A fill function must return one per row of the cell table",
        length(key), nrow(cells)
      ),
      call. = FALSE
    )
  }
  values <- sort(unique(key))
  col <- col %||% hcl.colors(length(values), "Set 2")
  if (length(col) < length(values)) {
    stop(
      sprintf("'col' gives %d colors for %d fill values", length(col),
              length(values)),
      call. = FALSE
    )
  }
  col[match(key, values)]
}

#' @keywords internal
#' @noRd
.chart_arg_message <- function(arg, presets) {
  sprintf(
    "'%s' must be %s, a column of the cell table (period, panel, entry, launch, rotation_group, wave, visit, life_stage, take) or a function of it",
    arg, presets
  )
}

#' @keywords internal
#' @noRd
.chart_period_labels <- function(period_labels, n_period) {
  if (is.null(period_labels)) {
    return(seq_len(n_period))
  }
  if (length(period_labels) != n_period) {
    stop(
      sprintf("'period_labels' has %d labels for %d periods",
              length(period_labels), n_period),
      call. = FALSE
    )
  }
  as.character(period_labels)
}

#' Entrants per cohort for the chart labels, or NULL for cohort shares
#'
#' A panel may only scale the rotation it was planned on, or the label column
#' and the total row would describe two designs.
#' @keywords internal
#' @noRd
.chart_take <- function(panel, w) {
  if (is.null(panel)) return(NULL)
  if (!inherits(panel, "svyplan_panel")) {
    stop("'panel' must be an n_panel() or prec_panel() result", call. = FALSE)
  }
  if (!identical(panel$design, "rotating")) {
    stop(
      "'panel' is a fixed panel, which recruits one cohort and has no rotation to chart; pass a result with design = \"rotating\"",
      call. = FALSE
    )
  }
  rot <- as.double(panel$rotation)
  if (!isTRUE(all.equal(rot, as.numeric(w)))) {
    stop(
      sprintf(
        "'panel' follows the rotation %s and the chart draws %s; they must describe the same design",
        .fmt_rotation(rot), .fmt_rotation(as.numeric(w))
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

#' The frontier's cost and cv, and the fitted design's
#'
#' A register fit plans its field design, which its grid reports, so its
#' frontier and fitted point are that design's.
#' @keywords internal
#' @noRd
.frontier_xy <- function(x, fr) {
  if (is.null(x$params$psu)) {
    return(list(cost = fr$cost, cv = fr$cv,
                fitted = c(x$params$achieved$cost, sqrt(x$objective_value))))
  }
  list(cost = fr$cost_int, cv = fr$cv,
       fitted = c(x$operational$cost,
                  sqrt(x$operational$objective_value)))
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
        psu = p$psu, certainty_cutoff = p$certainty_cutoff,
        n_psu_per_zone = p$n_psu_per_zone, control = p[["control"]],
        unit_cost = p$unit_cost, alpha = p$alpha, deff = p$deff,
        resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p$eligible_per_unit %||% 1,
        df = p$df, min_n_stratum = p$min_n_stratum, fpc = p$fpc %||% "unit"
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
              eligibility_rate = p[["eligibility_rate"]] %||% 1,
              eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
              alternative = p$alternative,
              overlap = p$overlap,
              overlap_cor = p$overlap_cor,
              overlap_draw = p$overlap_draw,
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
              eligibility_rate = p[["eligibility_rate"]] %||% 1,
              eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
              alternative = p$alternative,
              overlap = p$overlap,
              overlap_cor = p$overlap_cor,
              overlap_draw = p$overlap_draw
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
              eligibility_rate = p[["eligibility_rate"]] %||% 1,
              eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
              alternative = p$alternative,
              overlap = p$overlap,
              overlap_cor = p$overlap_cor,
              overlap_draw = p$overlap_draw
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
