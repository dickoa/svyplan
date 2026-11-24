#' Print, format and coerce a rotation pattern
#'
#' A `svyplan_rotation` from [design_rotation()] is a numeric vector of takes,
#' one entry per period of a unit's life. `as.double()` gives that vector
#' bare. `as.data.frame()` lays the pattern out over time as the table of
#' cells [plot.svyplan_rotation()] draws, one row per interview of every
#' panel in the periods drawn.
#'
#' @param x A `svyplan_rotation` object.
#' @param e1,e2 Operands. Arithmetic and comparison return bare numerics,
#'   since moved takes no longer have the life, standing sample and spec the
#'   rotation records.
#' @param i,value Index and replacement value. Replacement is refused: build
#'   the pattern you mean with [design_rotation()].
#' @param name Field to extract: `spec`, `life`, `in_sample`, `unbroken` or
#'   `equal_take`.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   `as.data.frame()` arguments.
#' @param start,horizon,panel_plan How the first panels enter, the
#'   periods covered and an optional [n_panel()] result, as for
#'   [plot.svyplan_rotation()].
#' @param ... Additional arguments are not supported and produce an error.
#' @return The method-specific result:
#' \describe{
#'   \item{`print()`}{Returns `x` invisibly.}
#'   \item{`format()`}{Returns a string.}
#'   \item{`as.double()`}{Returns the take vector.}
#'   \item{`$`}{Returns one field of the pattern.}
#'   \item{`Ops()` and `Math()`}{Return bare numerics.}
#'   \item{`as.data.frame()`}{Returns one row per interview, ordered by
#'     panel and then period, with columns:
#'     \describe{
#'       \item{`period`}{The period of the interview.}
#'       \item{`panel`}{The panel, numbered in the chart's order: under an
#'         immediate start the start-up panels first, the one with the
#'         shortest remaining life first, then one new panel per period.
#'         Each panel here is one cohort, the units at the same wave at
#'         every period, so the chart needs no separate cohort column.}
#'       \item{`entry`}{The period the panel is selected. A start-up panel
#'         whose remaining life starts out of sample is selected in period 1
#'         and first interviewed later.}
#'       \item{`startup`}{`TRUE` for a panel of an immediate start.}
#'       \item{`rotation_group`}{The rotation group the panel takes over,
#'         as defined on [plot.svyplan_rotation()].}
#'       \item{`wave`}{The place of the interview in the rotation pattern,
#'         counting in-sample periods only.}
#'       \item{`visit`}{The interview number of the panel. Equal to `wave`
#'         except for a start-up panel.}
#'       \item{`life_stage`}{The period of the full life, out-of-sample
#'         periods counted.}
#'       \item{`take`}{The take at that period, 1 under an equal take.}
#'       \item{`n`}{Units interviewed, `take` times the panel's entrants.
#'         Present only when `panel_plan` is given.}
#'     }
#'     The attributes `rotation`, `start`, `horizon` and
#'     `full_sample_from`, the first period with the full sample in the
#'     field, describe the table.}
#' }
#'
#' @seealso [design_rotation()], which builds these objects,
#'   [plot.svyplan_rotation()] for the chart of the cell table, and
#'   [design_overlap()] for the lag profile one produces.
#'
#' @examples
#' cps <- design_rotation("4-8-4")
#' cps
#' as.double(cps)
#'
#' # The fields of the pattern
#' cps$life
#' cps$in_sample
#' cps$unbroken
#'
#' # Arithmetic returns plain takes, never a rotation with stale fields
#' 2 * cps
#'
#' # The cells of an immediate start: a start-up panel's first interview can
#' # be at a later wave
#' cells <- as.data.frame(design_rotation("2-2-2"), start = "immediate",
#'                        horizon = 8)
#' head(cells, 8)
#'
#' # Which panels are at their first wave in period 6
#' subset(cells, period == 6 & wave == 1)
#'
#' @name print.svyplan_rotation
NULL

#' @rdname print.svyplan_rotation
#' @export
print.svyplan_rotation <- function(x, ...) {
  .check_unused_dots(...)
  cat("Rotation pattern (planning)\n\n")
  cat(.fmt_rotation_header(x))
  invisible(x)
}

#' @rdname print.svyplan_rotation
#' @export
format.svyplan_rotation <- function(x, ...) {
  .check_unused_dots(...)
  sprintf(
    "svyplan_rotation [%s]",
    .fmt_rotation(as.double(x))
  )
}

#' @rdname print.svyplan_rotation
#' @export
as.double.svyplan_rotation <- function(x, ...) {
  .check_unused_dots(...)
  out <- unclass(x)
  attributes(out) <- NULL
  out
}

#' @rdname print.svyplan_rotation
#' @export
as.data.frame.svyplan_rotation <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  start = c("gradual", "immediate"),
  horizon = NULL,
  panel_plan = NULL,
  ...
) {
  .check_unused_dots(...)
  if (missing(start)) {
    start <- NULL
  }
  start <- .check_chart_start(start, panel_plan)
  w <- as.double(x)
  horizon <- .whole_arg(horizon, "horizon", length(w) + 4L)
  out <- .rotation_cells(w, start, horizon, .chart_take(panel_plan, w))
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  structure(
    out,
    rotation = x,
    start = start,
    horizon = horizon,
    full_sample_from = if (identical(start, "immediate")) 1L else length(w)
  )
}

#' @rdname print.svyplan_rotation
#' @export
`$.svyplan_rotation` <- function(x, name) {
  fields <- c("spec", "life", "in_sample", "unbroken", "equal_take")
  if (length(name) != 1L || !name %in% fields) {
    stop(
      sprintf(
        "no field '%s' in a rotation pattern; available: %s",
        paste(name, collapse = ", "),
        paste(fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  attr(x, name, exact = TRUE)
}

#' @rdname print.svyplan_rotation
#' @export
Ops.svyplan_rotation <- function(e1, e2) {
  e1 <- if (inherits(e1, "svyplan_rotation")) as.double(e1) else e1
  if (missing(e2)) {
    return(do.call(.Generic, list(e1)))
  }
  e2 <- if (inherits(e2, "svyplan_rotation")) as.double(e2) else e2
  do.call(.Generic, list(e1, e2))
}

#' @rdname print.svyplan_rotation
#' @export
Math.svyplan_rotation <- function(x, ...) {
  do.call(.Generic, c(list(as.double(x)), list(...)))
}

#' @rdname print.svyplan_rotation
#' @export
`[<-.svyplan_rotation` <- function(x, i, value) {
  .no_replacement("rotation pattern", "design_rotation()")
}

#' @rdname print.svyplan_rotation
#' @export
`[[<-.svyplan_rotation` <- function(x, i, value) {
  .no_replacement("rotation pattern", "design_rotation()")
}

#' @rdname print.svyplan_rotation
#' @export
`$<-.svyplan_rotation` <- function(x, name, value) {
  .no_replacement("rotation pattern", "design_rotation()")
}

#' Print, format and coerce a rotation overlap
#'
#' A `svyplan_overlap` from [design_overlap()] is a numeric vector of
#' overlap fractions indexed by lag, so `x[1]` and `x[12]` are plain
#' numbers ready for an `overlap` argument. `$` reaches the counts and the
#' rotation behind them.
#'
#' @param x A `svyplan_overlap` object.
#' @param e1,e2 Operands. Arithmetic and comparison return bare numerics,
#'   the counts and the rotation describing the profile as computed and not
#'   whatever it was transformed into.
#' @param i Lag to extract, by position or by its name, so `x[12]` and
#'   `x[["12"]]` are both the twelve-period overlap. The result is a bare
#'   number, carrying neither the class nor the lag as a name.
#' @param value Replacement value. Replacement is refused, an overlap being
#'   computed from a rotation rather than assembled.
#' @param name Field to extract: `overlap`, `shared`, `in_sample`,
#'   `basis`, `rotation`, `lag`, `life`, `window` or `visits`.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   `as.data.frame()` arguments.
#' @param ... Additional arguments are not supported and produce an error.
#' @return The method-specific result:
#' \describe{
#'   \item{`print()`}{Returns `x` invisibly.}
#'   \item{`format()`}{Returns a string.}
#'   \item{`as.double()`}{Returns the overlap vector.}
#'   \item{`as.data.frame()`}{Returns one row per lag.}
#'   \item{`Ops()` and `Math()`}{Return bare numerics.}
#' }
#' Replacement is an error.
#'
#' @details
#' The values, the shared counts and the rotation describe one design, so
#' nothing may move the values while leaving the rest. Subsetting and
#' arithmetic return bare numerics, and replacement is an error. That also
#' settles `pmax()` and `pmin()`, which copy the attributes of their first
#' argument without dispatching to any method and would otherwise return
#' something still labeled an overlap whose counts no longer follow from its
#' values. They assign through `[<-`, so they are refused here too. Convert
#' with `as.double(x)` and they work as usual.
#'
#' @seealso [design_overlap()], which builds these objects,
#'   [design_rotation()] for the pattern behind them,
#'   [plot.svyplan_overlap()] for the bar chart by lag, and
#'   [print.svyplan_schedule()] for the operational schedule a rotation
#'   becomes.
#'
#' @examples
#' cps <- design_overlap("4-8-4")
#' cps
#'
#' # one lag, as a bare number, which is what the overlap arguments take
#' cps[1]
#' cps[["12"]]
#'
#' # the fields behind the values
#' cps$shared
#' cps$rotation
#' as.data.frame(cps)
#'
#' # arithmetic drops the class rather than carrying stale counts
#' as.double(cps)[1:3]
#'
#' @name print.svyplan_overlap
NULL

#' @rdname print.svyplan_overlap
#' @export
print.svyplan_overlap <- function(x, ...) {
  .check_unused_dots(...)
  rot <- .overlap_rotation(x)
  respondent <- identical(.overlap_basis(x), "respondent")
  cat(if (respondent) {
    "Respondent overlap (planning, permanent attrition)\n\n"
  } else {
    "Rotation overlap (planning)\n\n"
  })
  cat(.fmt_rotation_header(rot))
  if (respondent) {
    cat(sprintf(
      "  respondents: %.4g per entrant in sample each period\n",
      .overlap_in_sample(x)
    ))
  }
  cat("\n")
  # in_sample is constant by construction and is already in the header;
  # as.data.frame() keeps it, since a table read into code wants it.
  tab <- as.data.frame(x)[c("lag", "shared", "overlap")]
  shown <- utils::head(tab, 16L)
  if (respondent) {
    shown$shared <- sprintf("%.4g", shown$shared)
  }
  shown$overlap <- sprintf("%.4g", shown$overlap)
  print(shown, row.names = FALSE, right = FALSE)
  if (nrow(tab) > nrow(shown)) {
    cat(sprintf(
      "  %d further lag%s, read them with as.data.frame()\n",
      nrow(tab) - nrow(shown),
      if (nrow(tab) - nrow(shown) > 1L) "s" else ""
    ))
  }
  visits <- attr(x, "visits", exact = TRUE)
  if (!is.null(visits) && !is.na(visits)) {
    cat(sprintf(
      "\n  %.3g visits per unit over a window of %d periods\n",
      visits, attr(x, "window", exact = TRUE)
    ))
  }
  invisible(x)
}

#' The rotation a result was computed from
#'
#' One accessor, so a reader never has to know whether the life and the
#' takes ride on the overlap or on the rotation it came from.
#' @keywords internal
#' @noRd
.overlap_rotation <- function(x) attr(x, "rotation", exact = TRUE)

#' Whether an overlap counts issued units or expected respondents
#' @keywords internal
#' @noRd
.overlap_basis <- function(x) attr(x, "basis", exact = TRUE) %||% "issued"

#' Units, or expected respondents per entrant, in sample at one period
#' @keywords internal
#' @noRd
.overlap_in_sample <- function(x) {
  attr(x, "in_sample", exact = TRUE) %||%
    attr(.overlap_rotation(x), "in_sample", exact = TRUE)
}

#' The two lines that describe a rotation
#'
#' Shared by the rotation's own print method and the overlap's, so the two
#' cannot drift into describing the same design differently.
#' @keywords internal
#' @noRd
.fmt_rotation_header <- function(rot) {
  paste0(
    sprintf(
      "  life: %s periods, %s in sample each period\n",
      .fmt_count_n(attr(rot, "life", exact = TRUE)),
      .fmt_count_n(attr(rot, "in_sample", exact = TRUE))
    ),
    sprintf("  rotation: %s\n", .fmt_rotation(as.double(rot)))
  )
}

#' Render a rotation as the spells a planner declared
#'
#' Run-length form, since that is how a rotation is named and argued about.
#' The per-period vector it expands to is what the arithmetic reads.
#' @keywords internal
#' @noRd
.fmt_rotation <- function(w) {
  r <- rle(w)
  paste(
    vapply(
      seq_along(r$lengths),
      function(i) {
        if (r$values[i] == 0) {
          sprintf("%d out", r$lengths[i])
        } else if (isTRUE(all.equal(r$values[i], 1))) {
          sprintf("%d in", r$lengths[i])
        } else {
          sprintf("%d in at %.4g", r$lengths[i], r$values[i])
        }
      },
      character(1L)
    ),
    collapse = ", "
  )
}

#' @rdname print.svyplan_overlap
#' @export
format.svyplan_overlap <- function(x, ...) {
  .check_unused_dots(...)
  sprintf(
    "svyplan_overlap [life %g, lag 1 = %.4g]",
    attr(.overlap_rotation(x), "life", exact = TRUE),
    unclass(x)[[1L]]
  )
}

#' @rdname print.svyplan_overlap
#' @export
as.double.svyplan_overlap <- function(x, ...) {
  .check_unused_dots(...)
  # unclass() drops only the class; the counts and the schedule ride along
  # as attributes and would surface in anything expecting a bare vector.
  out <- unclass(x)
  attributes(out) <- NULL
  out
}

#' @rdname print.svyplan_overlap
#' @export
as.list.svyplan_overlap <- function(x, ...) {
  .check_unused_dots(...)
  rot <- .overlap_rotation(x)
  list(
    overlap = as.double(x),
    shared = attr(x, "shared", exact = TRUE),
    in_sample = .overlap_in_sample(x),
    basis = .overlap_basis(x),
    rotation = rot,
    lag = as.integer(names(x)),
    life = attr(rot, "life", exact = TRUE),
    window = attr(x, "window", exact = TRUE) %||% NA_integer_,
    visits = attr(x, "visits", exact = TRUE) %||% NA_real_
  )
}

#' @rdname print.svyplan_overlap
#' @export
`[.svyplan_overlap` <- function(x, i) {
  # Names survive arithmetic, and a subset is no longer the profile the
  # counts and schedule describe, so both are dropped.
  v <- unclass(x)
  attributes(v) <- list(names = names(x))
  if (missing(i)) {
    return(unname(v))
  }
  unname(v[i])
}

#' @rdname print.svyplan_overlap
#' @export
Ops.svyplan_overlap <- function(e1, e2) {
  # Bare numerics: keeping the class would leave moved values still claiming
  # the schedule that produced the originals.
  e1 <- if (inherits(e1, "svyplan_overlap")) as.double(e1) else e1
  if (missing(e2)) {
    return(do.call(.Generic, list(e1)))
  }
  e2 <- if (inherits(e2, "svyplan_overlap")) as.double(e2) else e2
  do.call(.Generic, list(e1, e2))
}

#' @rdname print.svyplan_overlap
#' @export
Math.svyplan_overlap <- function(x, ...) {
  do.call(.Generic, c(list(as.double(x)), list(...)))
}

#' @rdname print.svyplan_overlap
#' @export
`[[.svyplan_overlap` <- function(x, i) {
  unname(unclass(x)[[i]])
}

#' @rdname print.svyplan_overlap
#' @export
`$.svyplan_overlap` <- function(x, name) {
  fields <- as.list(x)
  if (length(name) != 1L || !name %in% names(fields)) {
    stop(
      sprintf(
        "no field '%s' in a rotation overlap; available: %s",
        paste(name, collapse = ", "),
        paste(names(fields), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  fields[[name]]
}

#' @rdname print.svyplan_overlap
#' @export
as.data.frame.svyplan_overlap <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  data.frame(
    lag = as.integer(names(x)),
    shared = attr(x, "shared", exact = TRUE),
    in_sample = .overlap_in_sample(x),
    overlap = as.double(x),
    row.names = row.names,
    stringsAsFactors = stringsAsFactors
  )
}

#' Print, summarize, format and coerce a panel recruitment
#'
#' Display and coercion methods for the object [n_panel()] and
#' [prec_panel()] return. `print()` gives the number to recruit, the
#' responding sample it is expected to leave, the precision at the target
#' and the wave-by-wave table. `summary()` adds what the design implies
#' around that answer: the standing sample a rotating design holds, the
#' response and retention it assumes, where the life's loss falls, and the
#' start-up a design passes through on its way to equilibrium.
#' The coercions return the recruitment count, which is a number of units to
#' release and not the analysis sample. Those differ by the whole of the
#' panel's attrition, and it is why `svyplan_panel` is a sibling of
#' `svyplan_n` rather than a subtype.
#'
#' @param x A `svyplan_panel` object, or the `summary.svyplan_panel` object
#'   `summary()` returns.
#' @param object A `svyplan_panel` object.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   `as.data.frame()` arguments.
#' @param ... Additional arguments are not supported and produce an error.
#' @return The method-specific result:
#' \describe{
#'   \item{`print()`}{Returns its argument invisibly.}
#'   \item{`summary()`}{Returns a `summary.svyplan_panel` object carrying the
#'     plan, full wave table, the table by period, and the visits in sample
#'     at each period.}
#'   \item{`format()`}{Returns a string.}
#'   \item{`as.double()`}{Returns the recruitment count.}
#'   \item{`as.integer()`}{Returns the whole units to which the recruitment
#'     count rounds up.}
#'   \item{`as.data.frame()`}{Returns the wave table.}
#' }
#'
#' @examples
#' plan <- n_panel(
#'   n_prop(p = 0.5, moe = 0.031),
#'   retention = c(0.878, 0.963, 0.936, 0.956),
#'   resp_rate = 0.728
#' )
#' plan
#' summary(plan)
#'
#' # A rotating design reaching equilibrium: the table by period shows the
#' # start-up before it gets there
#' rot <- n_panel(
#'   n_prop(p = 0.5, moe = 0.031),
#'   retention = c(0.878, 0.963, 0.936, 0.956),
#'   resp_rate = 0.728,
#'   design = "rotating",
#'   start = "immediate"
#' )
#' summary(rot)$periods
#' summary(rot)$composition
#'
#' as.integer(plan)
#' as.data.frame(plan)
#'
#' @seealso [n_panel()] and [prec_panel()], which build these objects, and
#'   [print.svyplan_schedule()] for the operational schedule a rotating plan
#'   becomes.
#'
#' @name print.svyplan_panel
NULL

#' @rdname print.svyplan_panel
#' @export
print.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  k <- nrow(x$waves)
  cat(sprintf("Panel recruitment (%s, %d-wave life)\n", x$design, k))
  cat(.fmt_panel_rotation(x))
  cat(.fmt_panel_headline(x))
  cat(.fmt_panel_whole_status(x))
  cat(.fmt_panel_precision(x))
  cat(.fmt_panel_assured(x))
  cat("\n")
  print(.fmt_panel_waves(x, brief = TRUE), row.names = FALSE, right = FALSE)
  cat("\n")
  cat("# summary() for the start-up, the loss and per-wave cv\n")
  invisible(x)
}

#' @rdname print.svyplan_panel
#' @export
summary.svyplan_panel <- function(object, ...) {
  .check_unused_dots(...)
  shown <- .panel_shown(object)
  structure(
    list(
      plan = object,
      recruit = shown$recruit,
      n_resp = round(shown$head),
      n_in_sample = shown$in_sample,
      waves = .fmt_panel_waves(object, brief = FALSE),
      periods = .fmt_panel_period_table(object),
      composition = object$period_visits
    ),
    class = "summary.svyplan_panel"
  )
}

#' @rdname print.svyplan_panel
#' @export
print.summary.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  plan <- x$plan
  cat(sprintf(
    "Analysis of a panel recruitment (%s, %d-wave life)\n\n",
    plan$design,
    nrow(plan$waves)
  ))
  cat(.fmt_panel_rotation(plan))
  cat(.fmt_panel_headline(plan))
  cat(.fmt_panel_whole_status(plan))
  cat(.fmt_panel_in_sample(plan))
  cat(.fmt_panel_histories(plan))
  cat(.fmt_panel_rates(plan))
  cat(.fmt_panel_target_rate(plan))
  cat(.fmt_panel_precision(plan))
  cat(.fmt_panel_assured(plan))

  cat(sprintf(
    "\n%s\n",
    if (identical(plan$design, "fixed")) {
      "Waves of the life"
    } else {
      "Waves alive at one period"
    }
  ))
  print(x$waves, row.names = FALSE, right = FALSE)

  if (!is.null(x$periods)) {
    cat(sprintf("\nBy period, %s start\n", plan$start))
    print(x$periods, row.names = FALSE, right = FALSE)
    cat("# $composition for the visits in sample at each period\n")
  }
  invisible(x)
}

#' Counts for the units the printed block describes
#'
#' Every count in the printed block comes from the same recruitment as the
#' precision beside it, the whole units a planner would release, or the
#' continuous requirement when whole units do not fit the frame. The stored
#' fields stay continuous, which is what keeps the round trip exact.
#' @keywords internal
#' @noRd
.panel_shown <- function(x) {
  op <- .panel_operational(x)
  recruit <- .panel_recruit(op)
  wave <- recruit * x$waves$q
  list(
    recruit = recruit,
    whole = attr(op, "whole", exact = TRUE),
    wave = wave,
    head = if (identical(x$design, "fixed")) {
      wave[[x$target_wave]]
    } else {
      sum(wave)
    },
    in_sample = if (identical(x$design, "rotating")) recruit * x$n_cohorts
  )
}

#' The plan evaluated at the whole units it would release
#'
#' Precision shown beside whole-unit counts is the precision of those
#' counts, so the printed block describes one design. A count already whole
#' is the plan itself. Should rounding up leave the plan's frame, the
#' continuous plan is the one shown, and `whole` is `FALSE` so the block
#' says so.
#' @keywords internal
#' @noRd
.panel_operational <- function(x) {
  recruit <- .whole_count(.panel_recruit(x))
  if (recruit == .panel_recruit(x)) {
    return(structure(x, whole = TRUE))
  }
  tryCatch(
    structure(.panel_at(x, recruit, assurance = NULL), whole = TRUE),
    error = function(e) structure(x, whole = FALSE)
  )
}

#' Say why a block shows the continuous requirement
#' @keywords internal
#' @noRd
.fmt_panel_whole_status <- function(x) {
  shown <- .panel_shown(x)
  if (isTRUE(shown$whole)) {
    return("")
  }
  whole <- .whole_count(.panel_recruit(x))
  held <- if (identical(x$design, "fixed")) whole else whole * length(x$rotation)
  sprintf(
    "whole units: %s %s hold %s, above N = %s, so these are continuous\n",
    .fmt_count_n(whole),
    if (identical(x$design, "fixed")) "issued" else "per period",
    .fmt_count_n(held), .fmt_count_n(x$target$params$N)
  )
}

#' State the recruitment and what it leaves, in that order
#'
#' The two designs answer with different quantities, so the line names the
#' quantity rather than printing a number a reader has to attribute. In the
#' reverse direction the recruitment was supplied and may not reach the
#' target, which the line says outright: a plan that misses is a fact about
#' the plan and not a detail, so it stays on the printed line.
#' @keywords internal
#' @noRd
.fmt_panel_headline <- function(x) {
  shown <- .panel_shown(x)
  where <- if (identical(x$design, "fixed")) {
    sprintf(" at wave %d", x$target_wave)
  } else {
    sprintf(", pooled over %d cohorts", x$n_cohorts)
  }
  count <- if (isTRUE(shown$whole)) {
    .fmt_count_n(shown$recruit)
  } else {
    format(signif(shown$recruit, 4))
  }
  lead <- if (identical(x$design, "fixed")) {
    sprintf("issued: %s", count)
  } else {
    sprintf("entrants: %s per period", count)
  }
  short <- shown$head < x$n_target * (1 - 1e-9)
  have <- round(shown$head)
  need <- .whole_count(x$n_target)
  # Integers that collide would print a shortfall between equal numbers, so
  # take the fewest decimals at which the two differ.
  if (short && have >= need) {
    d <- 1L
    while (d < 6L && round(shown$head, d) >= round(x$n_target, d)) {
      d <- d + 1L
    }
    have <- format(round(shown$head, d), nsmall = d)
    need <- format(round(x$n_target, d), nsmall = d)
  } else {
    have <- .fmt_count_n(have)
    need <- .fmt_count_n(need)
  }
  sprintf(
    "%s -> %s responding%s%s%s\n",
    lead,
    have,
    .fmt_panel_eligible(x, shown$head),
    where,
    if (short) sprintf(", short of the %s the target needs", need) else ""
  )
}

#' Name the eligible count a yield per unit implies
#'
#' The headline counts sampled units, and with several eligible units per
#' unit, or fewer than one, the precision beside it is computed on another
#' count. Shown only when the two differ.
#' @keywords internal
#' @noRd
.fmt_panel_eligible <- function(x, head) {
  epu <- x$target$params[["eligible_per_unit"]] %||% 1
  if (isTRUE(all.equal(epu, 1))) {
    return("")
  }
  sprintf(" (%s eligible)", .fmt_count_n(round(head * epu)))
}

#' Name the standing sample a rotating design carries
#'
#' The entrants are the release and the live cohorts are what stands behind
#' them. Only a rotating design has the second quantity.
#' @keywords internal
#' @noRd
.fmt_panel_in_sample <- function(x) {
  if (!identical(x$design, "rotating")) {
    return("")
  }
  out <- sprintf(
    "in sample: %s across %d live cohorts\n",
    .fmt_count_n(round(.panel_shown(x)$in_sample)),
    x$n_cohorts
  )
  life <- length(x$rotation)
  if (life > x$n_cohorts) {
    out <- paste0(out, sprintf(
      "held: %s across %d cohorts, %d of them out of sample\n",
      .fmt_count_n(round(.panel_shown(x)$recruit * life)), life,
      life - x$n_cohorts
    ))
  }
  out
}

#' Complete histories a rotating design delivers
#'
#' The level target pools every live cohort, and only the cohort at its last
#' wave carries a response at every wave of the life. Under permanent loss
#' that is the entrants times the last wave's `q`.
#' @keywords internal
#' @noRd
.fmt_panel_histories <- function(x) {
  if (!identical(x$design, "rotating")) {
    return("")
  }
  k <- nrow(x$waves)
  sprintf(
    "complete histories: %s per period, the cohort at wave %d with all %d responses\n",
    .fmt_count_n(round(.panel_shown(x)$wave[[k]])), k, k
  )
}

#' The rotation of a rotating design whose life has a gap
#' @keywords internal
#' @noRd
.fmt_panel_rotation <- function(x) {
  w <- as.double(x$rotation)
  if (length(w) == 0L || all(w > 0)) {
    return("")
  }
  sprintf(
    "life: %d periods (%s)\n", length(w), .fmt_rotation(w)
  )
}

#' Report an assurance only where one was asked for
#'
#' A level a finite frame cannot supply is named on the line that reports
#' it, the number being a requirement rather than a design.
#' @keywords internal
#' @noRd
.fmt_panel_assured <- function(x) {
  if (is.null(x$n_assured)) {
    return("")
  }
  sprintf(
    "assured (%s): %s %s%s\n",
    .fmt_prob(x$params$assurance),
    .fmt_count_n(.whole_count(x$n_assured)),
    if (identical(x$design, "fixed")) "issued" else "entrants per period",
    if (isFALSE(x$assured_feasible)) {
      sprintf(", beyond the population of %s", .fmt_count_n(x$target$params$N))
    } else {
      ""
    }
  )
}

#' @keywords internal
#' @noRd
.fmt_panel_rates <- function(x) {
  ret <- x$params$retention
  u <- unique(signif(ret, 10))
  ret_txt <- if (length(u) == 1L) {
    sprintf("%.3g", u)
  } else {
    sprintf("%.3g to %.3g", min(ret), max(ret))
  }
  share <- x$waves$loss_share[1L]
  loss_txt <- if (is.na(share)) {
    "no loss over the life"
  } else {
    sprintf("%.0f%% of the loss at wave 1", 100 * share)
  }
  elig <- x$params[["eligibility_rate"]]
  sprintf(
    "rates: response %.3g%s, retention %s (%s)\n",
    x$params$resp_rate,
    if (is.null(elig)) "" else sprintf(", eligibility %.3g", elig),
    ret_txt,
    loss_txt
  )
}

#' Name the response rate the plan did not use
#'
#' The warning at the call is gone by the time anyone reads the object, so
#' the disclosure travels with it. Same condition, from the same helper.
#' @keywords internal
#' @noRd
.fmt_panel_target_rate <- function(x) {
  rr <- .panel_rate_conflict(
    .panel_target_rate(x$target$params),
    x$params$resp_rate * (x$params[["eligibility_rate"]] %||% 1)
  )
  if (is.null(rr)) {
    return("")
  }
  sprintf(
    "target: %s removed, requirement %s responding\n",
    .panel_target_rate_label(x$target$params, lead = "response"),
    .fmt_count_n(round(x$n_target))
  )
}

#' @keywords internal
#' @noRd
.fmt_panel_precision <- function(x) {
  label <- if (is.null(x$method)) {
    x$type
  } else {
    sprintf("%s (%s)", x$type, x$method)
  }
  op <- .panel_operational(x)
  sprintf(
    "%s: se = %.4g, moe = %.4g%s\n",
    label,
    op$se,
    op$moe,
    if (is.na(op$cv)) "" else sprintf(", cv = %.3g", op$cv)
  )
}

#' The wave table at two widths
#'
#' `brief` is what `print()` shows: the retention that produced each wave,
#' the units left and the precision they buy. The cumulative `q`, the cv and
#' the expected cases are all derivable from those, so they belong to
#' `summary()`.
#' @keywords internal
#' @noRd
.fmt_panel_waves <- function(x, brief = FALSE) {
  w <- .panel_operational(x)$waves
  out <- data.frame(
    wave = w$wave,
    retention = ifelse(is.na(w$retention), "", sprintf("%.3g", w$retention)),
    n_resp = format(round(.panel_shown(x)$wave)),
    se = sprintf("%.4g", w$se),
    moe = sprintf("%.4g", w$moe),
    stringsAsFactors = FALSE
  )
  if (brief) {
    return(out)
  }
  out <- cbind(
    out[c("wave", "retention")],
    q = sprintf("%.4g", w$q),
    loss = ifelse(is.na(w$loss_share), "", sprintf("%.3g", w$loss_share)),
    out[c("n_resp", "se", "moe")],
    stringsAsFactors = FALSE
  )
  if (!all(is.na(w$cv))) {
    out$cv <- sprintf("%.3g", w$cv)
  }
  if (!is.null(w$expected_cases)) {
    out$cases <- format(round(w$expected_cases))
  }
  out
}

#' The table by period, on the display path
#'
#' Every count and its precision come from the plan at whole units, as in
#' the wave table, so the period that settles reports the sample the
#' headline promises.
#' @keywords internal
#' @noRd
.fmt_panel_period_table <- function(x) {
  if (is.null(x$periods)) {
    return(NULL)
  }
  l <- .panel_operational(x)$periods
  out <- data.frame(
    period = l$period,
    entrants = format(round(l$n_entrants)),
    in_sample = format(round(l$n_in_sample)),
    n_resp = format(round(l$n_resp)),
    se = sprintf("%.4g", l$se),
    moe = sprintf("%.4g", l$moe),
    stringsAsFactors = FALSE
  )
  if (!all(is.na(l$cv))) {
    out$cv <- sprintf("%.3g", l$cv)
  }
  out$equilibrium <- l$equilibrium
  out
}

#' @rdname print.svyplan_panel
#' @export
format.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  sprintf(
    "svyplan_panel [%s, %d waves, recruit %g]",
    x$design,
    nrow(x$waves),
    .panel_recruit(x)
  )
}

#' @rdname print.svyplan_panel
#' @export
as.double.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  .panel_recruit(x)
}

#' @rdname print.svyplan_panel
#' @export
as.integer.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  # The assured count, when one was asked for, stays in $n_assured: which
  # number to field is the planner's call and must not turn on whether an
  # argument was set.
  .coerce_count(.whole_count(.panel_recruit(x)))
}

#' @rdname print.svyplan_panel
#' @export
as.data.frame.svyplan_panel <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  out <- x$waves
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}

#' Refuse in-place modification of a computed planning quantity
#'
#' These classes carry values and the quantities they were computed from, and
#' a replacement changes the first while leaving the second: a doubled overlap
#' still claiming its schedule, a design effect no longer the product of its
#' components. Neither is a thing the package produces, so the assignment is
#' named rather than performed. Arithmetic returns bare numerics for the same
#' reason, and is the way to work with the values.
#' @keywords internal
#' @noRd
.no_replacement <- function(what, from) {
  stop(
    sprintf(
      "a %s cannot be modified in place: its values and the quantities behind them describe one design, and a replacement would leave them contradicting each other; recompute with %s, or take as.double(x) to work with the numbers",
      what,
      from
    ),
    call. = FALSE
  )
}

#' @rdname print.svyplan_overlap
#' @export
`[<-.svyplan_overlap` <- function(x, i, value) {
  .no_replacement("rotation overlap", "design_overlap()")
}

#' @rdname print.svyplan_overlap
#' @export
`[[<-.svyplan_overlap` <- function(x, i, value) {
  .no_replacement("rotation overlap", "design_overlap()")
}

#' @rdname print.svyplan_overlap
#' @export
`$<-.svyplan_overlap` <- function(x, name, value) {
  .no_replacement("rotation overlap", "design_overlap()")
}

#' @rdname print.svyplan_deff
#' @export
`[<-.svyplan_deff` <- function(x, i, value) {
  .no_replacement("planning design effect", "design_effect()")
}

#' @rdname print.svyplan_deff
#' @export
`[[<-.svyplan_deff` <- function(x, i, value) {
  .no_replacement("planning design effect", "design_effect()")
}

#' @rdname print.svyplan_deff
#' @export
`$<-.svyplan_deff` <- function(x, name, value) {
  .no_replacement("planning design effect", "design_effect()")
}

#' @rdname print.svyplan_df
#' @export
`[<-.svyplan_df` <- function(x, i, value) {
  .no_replacement("design degrees of freedom", "design_df()")
}

#' @rdname print.svyplan_df
#' @export
`[[<-.svyplan_df` <- function(x, i, value) {
  .no_replacement("design degrees of freedom", "design_df()")
}

#' @rdname print.svyplan_df
#' @export
`$<-.svyplan_df` <- function(x, name, value) {
  .no_replacement("design degrees of freedom", "design_df()")
}

