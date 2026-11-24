#' Print and coerce planning design effects
#'
#' Display and coercion methods for the composed design effect that
#' [design_effect()] returns. `print()` gives the overall value, while
#' `summary()` gives an ANOVA-style decomposition and the assumptions behind
#' each component. The coercion and arithmetic methods let the object stand
#' in for the overall value wherever a plain number is expected.
#'
#' @param x A `svyplan_deff` object from [design_effect()], or its summary for
#'   the summary print method.
#' @param object A `svyplan_deff` object to summarize.
#' @param value Replacement value. Replacement is refused, a design effect
#'   being the product of the components it carries.
#' @param row.names,optional Standard [as.data.frame()] arguments.
#' @param stringsAsFactors Logical. Retained for compatibility when a result
#'   is converted through [data.frame()].
#' @param validRN Logical. Accepted for compatibility with [data.frame()] in
#'   R 4.7.0 and later. Svyplan results already have valid row names.
#' @param e1,e2 Objects supplied to an arithmetic or comparison operator.
#' @param name,i A field name, one of those listed under Details. `[[` also
#'   accepts a numeric index, which reads the underlying numeric vector.
#' @param ... For mathematical transformations, additional arguments passed to
#'   the underlying operation. The other methods do not support additional
#'   arguments.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns a
#'   `summary.svyplan_deff` object containing the overall value, named
#'   components, component notes, and the basis of the decomposition.
#'   `format()` returns a character scalar, and `as.double()` returns the
#'   overall design effect.
#'   `as.data.frame()` returns a one-row table with the overall value and one
#'   column per component. `as.list()` returns the same fields as a named
#'   list, and `$` and `[[` return one of them. Arithmetic and mathematical
#'   transformations return ordinary numeric results.
#'
#' @details
#' A `svyplan_deff` behaves as the numeric overall design effect wherever one
#' is expected. It can be passed to any `deff` argument, compared, and
#' arithmetically combined, with the components dropped by any such
#' operation.
#'
#' The decomposition is reached by name, under one set of names shared by
#' every access route: `deff` for the overall value and `deff_<component>`
#' for each part, so `d$deff_cluster`, `d[["deff_cluster"]]`,
#' `as.list(d)$deff_cluster`, and `as.data.frame(d)$deff_cluster` are the
#' same number. Naming a field that this design effect does not have is an
#' error listing the ones it does. The overall value is always present, so a
#' one-component design effect reports it twice, once as `deff` and once as
#' the component it is made of.
#'
#' @seealso [design_effect()], which builds these objects, and
#'   [print.svyplan] for the sample size and precision results.
#'
#' @examples
#' d <- design_effect(icc = 0.03, take = 20,
#'                    weights = rep(c(1, 3), c(400, 100)))
#' d
#' summary(d)
#'
#' # one component, by name
#' d$deff
#' d$deff_cluster
#'
#' # or the whole decomposition, as a list or a one-row table
#' as.list(d)
#' as.data.frame(d)
#'
#' # it is the overall value wherever a number is expected
#' as.double(d)
#' n_prop(p = 0.3, moe = 0.05, deff = d)
#'
#' @name print.svyplan_deff
NULL

#' @rdname print.svyplan_deff
#' @export
print.svyplan_deff <- function(x, ...) {
  .check_unused_dots(...)
  cat(sprintf("Planning design effect: %.4f\n", as.double(x)))
  invisible(x)
}

#' @rdname print.svyplan_deff
#' @export
summary.svyplan_deff <- function(object, ...) {
  .check_unused_dots(...)
  components <- attr(object, "components", exact = TRUE)
  structure(
    list(
      deff = as.double(object),
      components = components,
      notes = attr(object, "notes", exact = TRUE),
      basis = .deff_summary_basis(components)
    ),
    class = "summary.svyplan_deff"
  )
}

#' @rdname print.svyplan_deff
#' @export
print.summary.svyplan_deff <- function(x, ...) {
  .check_unused_dots(...)
  labels <- .deff_component_labels(names(x$components))
  values <- c(unname(x$components), x$deff)
  tab <- data.frame(
    `Design effect` = sprintf("%.4f", values),
    check.names = FALSE,
    row.names = c(labels, "Overall")
  )

  cat("Analysis of design effects\n\n")
  print(tab, quote = FALSE, right = TRUE)
  if (length(x$components) > 1L) {
    cat("\nComponents combine multiplicatively.\n")
  }
  cat(sprintf("Basis: %s\n", x$basis))

  notes <- x$notes
  if (length(notes) > 0L && any(nzchar(notes))) {
    width <- max(nchar(labels))
    cat("\nComponent assumptions:\n")
    for (i in seq_along(notes)) {
      cat(sprintf("  %-*s  %s\n", width, labels[i], notes[[i]]))
    }
  }
  invisible(x)
}

#' Labels shared by the design-effect summary table and its assumptions
#' @keywords internal
#' @noRd
.deff_component_labels <- function(components) {
  labels <- c(
    cluster = "Clustering",
    weight = "Weighting",
    strata = "Stratification",
    calibration = "Calibration",
    allocation = "Allocation"
  )
  out <- unname(labels[components])
  unknown <- is.na(out)
  if (any(unknown)) {
    out[unknown] <- paste0(
      toupper(substr(components[unknown], 1L, 1L)),
      substring(components[unknown], 2L)
    )
  }
  out
}

#' Describe what the rows in a design-effect summary mean
#' @keywords internal
#' @noRd
.deff_summary_basis <- function(components) {
  component_names <- names(components)
  if ("allocation" %in% component_names) {
    if (length(component_names) == 1L) {
      return("direct allocation variance ratio")
    }
    has_w <- "weight" %in% component_names
    has_c <- "calibration" %in% component_names
    if (has_w && has_c) {
      return("direct allocation variance ratio with weighting and calibration adjustments")
    }
    if (has_c) {
      return("direct allocation variance ratio with a calibration adjustment")
    }
    return("direct allocation variance ratio with a Kish weighting adjustment")
  }
  if (length(component_names) > 1L) {
    return("approximate multiplicative decomposition")
  }
  switch(
    component_names,
    cluster = "cluster variance model",
    weight = "Kish weighting factor",
    strata = "proportional-allocation stratification factor",
    calibration = "auxiliary calibration variance reduction",
    "planning component"
  )
}

#' @rdname print.svyplan_deff
#' @export
format.svyplan_deff <- function(x, ...) {
  .check_unused_dots(...)
  components <- attr(x, "components", exact = TRUE)
  sprintf(
    "svyplan_deff [%s, %.4f]",
    paste(names(components), collapse = " x "),
    as.double(x)
  )
}

#' @rdname print.svyplan_deff
#' @export
as.double.svyplan_deff <- function(x, ...) {
  .check_unused_dots(...)
  unclass(x)[[1L]]
}

#' @rdname print.svyplan_deff
#' @export
as.list.svyplan_deff <- function(x, ...) {
  .check_unused_dots(...)
  components <- attr(x, "components", exact = TRUE)
  out <- c(list(deff = as.double(x)), as.list(components))
  names(out) <- c("deff", paste0("deff_", names(components)))
  out
}

#' @rdname print.svyplan_deff
#' @export
`$.svyplan_deff` <- function(x, name) {
  .deff_field(x, name)
}

#' @rdname print.svyplan_deff
#' @export
`[[.svyplan_deff` <- function(x, i, ...) {
  if (is.character(i)) {
    return(.deff_field(x, i))
  }
  unclass(x)[[i, ...]]
}

#' Look a design effect field up by name
#' @keywords internal
#' @noRd
.deff_field <- function(x, name) {
  fields <- as.list(x)
  if (length(name) != 1L || !name %in% names(fields)) {
    stop(
      sprintf(
        "no field '%s' in a design effect; available: %s",
        paste(name, collapse = ", "),
        paste(names(fields), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  fields[[name]]
}

#' @rdname print.svyplan_deff
#' @export
as.data.frame.svyplan_deff <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  as.data.frame(
    as.list(x),
    row.names = row.names,
    optional = optional,
    stringsAsFactors = stringsAsFactors
  )
}

#' @rdname print.svyplan_deff
#' @export
Ops.svyplan_deff <- function(e1, e2) {
  e1 <- if (inherits(e1, "svyplan_deff")) as.double(e1) else e1
  if (missing(e2)) {
    return(do.call(.Generic, list(e1)))
  }
  e2 <- if (inherits(e2, "svyplan_deff")) as.double(e2) else e2
  do.call(.Generic, list(e1, e2))
}

#' @rdname print.svyplan_deff
#' @export
Math.svyplan_deff <- function(x, ...) {
  do.call(.Generic, c(list(as.double(x)), list(...)))
}

#' @rdname print.svyplan
#' @name print.svyplan
#' @aliases print.svyplan_strata
#' @export
print.svyplan_strata <- function(x, ...) {
  .check_unused_dots(...)
  .print_boundary_strata(x)
  invisible(x)
}

#' @keywords internal
#' @noRd
.print_boundary_strata <- function(x) {
  method_label <- switch(
    x$method,
    cumrootf = "Dalenius-Hodges",
    geo = "Geometric",
    lh = "LH-inspired coordinate search",
    kozak = "Kozak random search",
    x$method
  )
  # Part of the method's name, and NA stays silent rather than reading as a
  # failure to converge.
  converged <- if (isTRUE(x$converged)) {
    ", converged"
  } else if (isFALSE(x$converged)) {
    ", not converged"
  } else {
    ""
  }
  cat(sprintf(
    "Strata boundaries (%s, %d strata%s)\n",
    method_label,
    x$n_strata,
    converged
  ))
  # No `Boundaries:` line: the cut points are the lower and upper columns of
  # the table below, and naming them twice at two precisions invites the
  # reader to look for a difference.
  alloc <- if (!is.null(x$alloc) && is.character(x$alloc)) {
    if (x$alloc == "power") {
      sprintf(", allocation: power (alloc_q = %.2f)", x$params$alloc_q)
    } else {
      sprintf(", allocation: %s", x$alloc)
    }
  } else {
    ""
  }
  cat(sprintf("n = %d, cv = %.4f%s\n", .whole_count(x$n), x$cv, alloc))
  .print_search_limit(x$search)
  df <- x$strata
  df$lower <- .fmt_boundary(df$lower)
  df$upper <- .fmt_boundary(df$upper)
  df$share <- sprintf("%.3f", df$share)
  df$sd <- sprintf("%.1f", df$sd)
  df$mean <- sprintf("%.1f", df$mean)
  df$take_all <- NULL
  cat("\n")
  # Every column here is a number, so they line up on the right.
  print(df, row.names = FALSE, right = TRUE)
}

#' A cut point at reading precision
#'
#' Five significant figures, and never in scientific notation: a boundary is
#' a value on the frame's own scale, and `3.018e+04` is not a number anyone
#' will compare against their data.
#' @keywords internal
#' @noRd
.fmt_boundary <- function(v) {
  # Element by element: formatting the vector at once pads every cut point to
  # the widest one's decimals, which puts a 428.4100 next to a 30185.00.
  vapply(
    v,
    function(z) format(signif(z, 5), scientific = FALSE, trim = TRUE),
    character(1L),
    USE.NAMES = FALSE
  )
}

#' @keywords internal
#' @noRd
.print_alloc_n <- function(x) {
  if (identical(x$method, "bethel")) {
    return(.print_bethel_n(x))
  }
  p <- x$params
  alloc <- p$alloc %||% x$method
  detail <- x$detail
  cluster <- !is.null(detail) && "n_1" %in% names(detail)
  H <- if (!is.null(detail)) nrow(detail) else NA
  cat(sprintf("Stratum allocation (%s", alloc))
  if (cluster) {
    S <- x$stages %||% 2L
    cat(", ", switch(as.character(S), "2" = "two", "3" = "three",
                     "4" = "four", as.character(S)), "-stage", sep = "")
  }
  if (!is.na(H)) {
    cat(sprintf(", %d strata", H))
  }
  cat(")\n")
  op <- x$operational
  if (!is.null(op)) {
    cat(sprintf("field design: n = %s", .fmt_count_n(op$n)))
    if (cluster && "n_1_int" %in% names(detail)) {
      cat(sprintf(", n_1 = %s", .fmt_count_n(sum(detail$n_1_int))))
    }
    rates <- .alloc_summary_response_rates(x)$combined
    if (any(rates != 1)) {
      cat(sprintf(", expected respondents = %.1f", sum(detail$n_int * rates)))
    }
    if (!is.na(op$cv)) {
      cat(sprintf(", cv = %.4f", op$cv))
    }
    if (!is.na(op$cost)) {
      cat(sprintf(", cost = %.0f", op$cost))
    }
    cat("\n")
    .print_search_limit(x$optimization$search)
    .print_search_skipped(x$optimization)
    cat(sprintf("continuous design: n = %s", .fmt_continuous_n(x$n)))
    if (!is.na(x$cv)) {
      cat(sprintf(", cv = %.4f", x$cv))
    }
    if (!is.na(x$se)) {
      cat(sprintf(", se = %.4f", x$se))
    }
    cat("\n")
  } else {
    cat(sprintf("n = %d", .whole_count(x$n)))
    if (cluster && "n_1_int" %in% names(detail)) {
      cat(sprintf(", n_1 = %s", .fmt_count_n(sum(detail$n_1_int))))
    }
    rates <- .alloc_summary_response_rates(x)$combined
    if (any(rates != 1)) {
      cat(sprintf(", expected respondents = %.1f", sum(detail$n_int * rates)))
    }
    if (!is.na(x$cv)) {
      cat(sprintf(", cv = %.4f", x$cv))
    }
    if (!is.na(x$se)) {
      cat(sprintf(", se = %.4f", x$se))
    }
    cat("\n")
  }
  parts <- character(0L)
  if (length(p$min_n_stratum) == 1L && isTRUE(p$min_n_stratum > 0)) {
    parts <- c(parts, sprintf("min_n_stratum = %g", p$min_n_stratum))
  } else if (length(p$min_n_stratum) > 1L && !all(is.na(p$min_n_stratum))) {
    parts <- c(parts, "min_n_stratum per stratum")
  }
  resp_rate <- p$resp_rate
  if (!is.null(p[["eligibility_rate"]])) {
    parts <- c(
      parts,
      sprintf("eligibility_rate = %s",
              .fmt_by_stratum(p[["eligibility_rate"]], .fmt_prob))
    )
  }
  epu <- .alloc_summary_response_rates(x)$per_unit
  if (any(epu != 1)) {
    parts <- c(
      parts,
      sprintf("eligible_per_unit = %s",
              .fmt_by_stratum(epu, function(v) sprintf("%.4g", v)))
    )
  }
  if (!is.null(resp_rate) && any(resp_rate < 1)) {
    parts <- c(
      parts,
      sprintf("resp_rate = %s", .fmt_by_stratum(resp_rate, .fmt_prob))
    )
  }
  parts <- c(parts, .fmt_deff(p$deff))
  if (length(parts) > 0L) {
    cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))
  }
  if (!is.null(detail$n_step)) {
    steps <- sort(unique(detail$n_step[detail$n_step > 1L]))
    cat(sprintf("field sizes in steps of %s\n", paste(steps, collapse = ", ")))
  }
  if (.alloc_size_only(x)) {
    cat("precision: not assessed (no 'sd', 'var' or 'p' in frame)\n")
  }
  .print_psu_fraction_note(detail, fpc = p$fpc %||% "unit")
  .print_design_df(x)
  .print_domains_block(x$domains)
}

#' Disclose an appreciable first-stage sampling fraction
#'
#' Under `fpc = "unit"` or `"none"`, `N_1` bounds the allocation but
#' activates no first-stage correction, so a design taking a large share of
#' the available PSUs is planned conservatively. `fpc = "stage"` applies that
#' correction, and the note is not printed.
#' The ultimate-unit `1 - n / N` does attenuate the between-PSU component, but
#' it is smaller than the first-stage fraction this note reports, so the
#' between-PSU term is still overstated. That is a property of the result worth
#' seeing, but not an event to act on, so it is disclosed here rather than
#' warned about. `predict()` re-runs the allocation once per grid row, and a
#' warning would repeat with it.
#' @keywords internal
#' @noRd
.print_psu_fraction_note <- function(detail, fpc = "unit", threshold = 0.1) {
  if (identical(fpc, "stage") || is.null(detail) ||
      !".frac_1" %in% names(detail)) {
    return(invisible(NULL))
  }
  hit <- which(is.finite(detail$.frac_1) & detail$.frac_1 >= threshold)
  if (length(hit) == 0L) {
    return(invisible(NULL))
  }
  bound <- !is.na(detail$.bound_source) & detail$.bound_source == "N_1"
  cat(sprintf(
    paste0(
      "note: samples %s of the available PSUs in %s%s\n",
      "      precision uses a with-replacement first stage and may be ",
      "conservative\n"
    ),
    paste0(
      format(round(100 * detail$.frac_1[hit]), trim = TRUE),
      "%",
      collapse = ", "
    ),
    paste(detail$stratum[hit], collapse = ", "),
    if (any(bound)) {
      sprintf(
        " (bound active in %s)",
        paste(detail$stratum[bound], collapse = ", ")
      )
    } else {
      ""
    }
  ))
  invisible(NULL)
}

#' Print the per-domain precision table an allocation carries
#' @keywords internal
#' @noRd
.print_domains_block <- function(dom) {
  if (is.null(dom)) {
    return(invisible(NULL))
  }
  cat(sprintf("Domains: %d\n", nrow(dom)))
  cat("\n")
  if (".cv" %in% names(dom)) {
    dom$.cv <- sprintf("%.4f", dom$.cv)
  }
  if (".cost" %in% names(dom)) {
    dom$.cost <- sprintf("%.0f", dom$.cost)
  }
  print(dom, row.names = FALSE, right = FALSE)
  invisible(NULL)
}

#' Print a constraint block, naming what was left out
#'
#' `sel` is the subset chosen for display (violations, or binding rows when
#' everything passes) out of `total`. Says so whenever rows are hidden, so a
#' one-row block under a "3 constraints" header does not read as the whole
#' picture.
#' @keywords internal
#' @noRd
.print_constraint_rows <- function(sel, total, label, hint) {
  if (nrow(sel) == 0L) {
    return(invisible(NULL))
  }
  shown <- utils::head(
    sel[c("constraint", ".metric", ".target", ".achieved", ".pass")],
    6L
  )
  if (nrow(sel) < total) {
    cat(sprintf("showing %d %s of %d\n", nrow(shown), label, total))
  }
  cat("\n")
  print(shown, row.names = FALSE, right = FALSE)
  cat("\n")
  if (nrow(shown) < total) {
    cat(hint, "\n", sep = "")
  }
  invisible(NULL)
}

#' The search line of a print, only when the node limit stopped the search
#'
#' A plan the node limit stopped short of proving is news, an optimal one
#' is not.
#' @keywords internal
#' @noRd
.print_search_limit <- function(search) {
  if (identical(search$status, "gap")) {
    cat(sprintf(
      "search: within %.3g%% of a bound of %s\n",
      100 * search$gap,
      format(signif(search$bound, 6), big.mark = ",")
    ))
  }
  if (!identical(search$status, "limit")) return(invisible(NULL))
  cat(sprintf(
    "search: stopped at the node limit, gap %.3g%% to a bound of %s\n",
    100 * search$gap,
    format(signif(search$bound, 6), big.mark = ",")
  ))
}

#' The search line of a whole design some of which was not searched
#' @keywords internal
#' @noRd
.print_search_skipped <- function(optimization, indent = "") {
  records <- optimization$domains %||% list(all = optimization$search)
  skipped <- Filter(function(r) identical(r$status, "heuristic"), records)
  if (!length(skipped)) return(invisible(NULL))
  reason <- vapply(skipped, function(r) r$reason %||% "listing", character(1))
  label <- function(nm) if (identical(nm, "all")) NULL else sQuote(nm)
  out <- character(0)
  pruning <- reason == "pruning"
  repaired <- reason == "floor repair"
  rest <- !(pruning | repaired)
  if (any(rest)) {
    why <- vapply(reason[rest], function(r) {
      switch(r,
             floor = "its size floor set the design",
             listing = "too many take vectors to list",
             frontier = "too many frontier states",
             numerical = "too many designs at the variance floor")
    }, character(1))
    where <- if (identical(names(skipped)[rest], "all")) why else
      paste0(sQuote(names(skipped)[rest]), ": ", why)
    out <- c(out, sprintf(
      "search: not run (%s), so that design is the rounded continuous one",
      paste(where, collapse = "; ")
    ))
  }
  if (any(pruning)) {
    out <- c(out, sprintf(
      "search: incomplete (%s), so the feasible starting design is used",
      paste(vapply(names(skipped)[pruning], function(nm) {
        label(nm) %||% "numerical pruning"
      }, character(1)), collapse = ", ")
    ))
  }
  if (any(repaired)) {
    out <- c(out, sprintf(
      "search: ran, then the floor repair raised %s, so that design is not proven",
      paste(vapply(names(skipped)[repaired], function(nm) {
        label(nm) %||% "the design"
      }, character(1)), collapse = ", ")
    ))
  }
  for (line in out) {
    cat(strwrap(line, width = 80L, indent = nchar(indent),
                exdent = nchar(indent) + 2L), sep = "\n")
  }
  invisible(NULL)
}

#' How far a feasible continuous design may sit above the optimum
#' @keywords internal
#' @noRd
.bethel_gap_note <- function(status, gap) {
  if (!identical(status, "feasible") || !is.numeric(gap) ||
      length(gap) != 1L || !is.finite(gap)) {
    return(NULL)
  }
  sprintf("within %.3g%% of the dual bound", 100 * gap)
}

#' Print a joint constrained allocation
#' @keywords internal
#' @noRd
.print_bethel_n <- function(x) {
  opt <- x$optimization
  op <- x$operational
  budget_mode <- identical(x$params$mode, "budget_objective")
  status <- opt$classification %||% "unknown"
  search <- opt$search
  register <- !is.null(x$params$psu)
  cat(sprintf(
    "Joint constrained allocation (Bethel%s)\n",
    if (register) ", PSU register" else ""
  ))
  # Report exceptions, not confirmations: a line that only ever says
  # "nothing went wrong" buries the lines that do carry news.
  gap_note <- .bethel_gap_note(status, opt$relative_duality_gap)
  if (!is.null(gap_note)) {
    cat(sprintf("continuous: feasible, %s\n", gap_note))
  } else if (!identical(status, "optimal")) {
    cat(sprintf("status: %s\n", status))
  }
  # The two modes return the same class and the same numbers mean different
  # things in each, so the question stays even though it never varies.
  cat(
    if (budget_mode) {
      sprintf(
        "question: best design affordable within a budget of %.6g\n",
        x$params$budget
      )
    } else {
      "question: cheapest design meeting every precision target\n"
    }
  )
  n_targets <- nrow(op$constraints)
  cat(sprintf(
    "field design: n = %s, cost = %.0f%s\n",
    .fmt_count_n(op$n),
    op$cost,
    if (n_targets == 0L) {
      ""
    } else {
      sprintf(
        " (%d target%s, %s)",
        n_targets,
        if (n_targets == 1L) "" else "s",
        if (isTRUE(op$all_pass)) "all pass" else "violations"
      )
    }
  ))
  if (register) {
    # Counts are PSUs available in each part. How many of the remainder are
    # drawn is a selection decision this plan does not make.
    n_certain <- sum(x$psu$certainty)
    n_supplied <- sum(x$psu$.certainty_source %in% "supplied")
    cat(sprintf(
      "PSUs: %d certainty%s, %d to draw from %d in the remainder\n",
      n_certain,
      if (n_supplied > 0L) sprintf(" (%d supplied)", n_supplied) else "",
      op$n_psu_draw,
      nrow(x$psu) - n_certain
    ))
  }
  .print_search_limit(search)
  continuous_cost <- x$params$achieved$cost
  increase <- if (continuous_cost > 0) {
    100 * (op$cost / continuous_cost - 1)
  } else {
    NA_real_
  }
  # A register's field design is planned whole, not rounded from the
  # continuous solve, which is only its relaxation.
  cat(sprintf(
    "continuous %s: n = %s, cost = %.0f%s\n",
    if (register) "relaxation" else if (is.null(gap_note)) "optimum" else "design",
    .fmt_continuous_n(x$n),
    continuous_cost,
    if (!is.na(increase) && abs(increase) >= 0.005) {
      sprintf(
        if (register) " (field design %+.2f%%)" else
          " (integerizing costs %+.2f%%)",
        increase
      )
    } else {
      ""
    }
  ))
  if (budget_mode) {
    cat(sprintf(
      "objective: weighted relative variance %.6g continuous, %.6g operational\n",
      x$objective_value,
      op$objective_value
    ))
    # For a register both describe the relaxation, so summary() keeps them.
    if (!register) {
      if (!isTRUE(opt$budget_binding)) {
        at_upper <- length(opt$active_upper %||% integer(0)) ==
          length(x$detail$n)
        cat(if (at_upper) {
          "the budget is not binding: the allocation sits at its upper bounds\n"
        } else {
          "the budget is not binding: the objective reaches its floor below it\n"
        })
      } else if (is.finite(opt$budget_sensitivity %||% NA_real_)) {
        cat(sprintf(
          "one more unit of budget changes the objective by %.4g\n",
          opt$budget_sensitivity
        ))
      }
    }
    obj <- x$objective
    if (!is.null(obj) && nrow(obj) > 0L) {
      shown <- utils::head(
        obj[c("component", "priority", ".cv", ".share")],
        6L
      )
      cat("\n")
      print(shown, row.names = FALSE, right = FALSE)
      cat("\n")
      if (nrow(obj) > 6L) cat("... see $objective for all components\n")
    }
  }
  violated <- any(!op$constraints$.pass)
  sel <- if (violated) {
    op$constraints[!op$constraints$.pass, , drop = FALSE]
  } else {
    x$constraints[x$constraints$.binding, , drop = FALSE]
  }
  if (violated || nrow(sel) > 1L) {
    # More than one binding constraint, or any failure, is worth a table.
    .print_constraint_rows(
      sel,
      total = nrow(op$constraints),
      label = if (violated) "violated" else "binding",
      hint = "... see $constraints and $operational$constraints for all rows"
    )
  } else if (nrow(sel) == 1L) {
    cat(sprintf(
      "binding: %s (target %.4g, achieved %.4g)\n",
      sel$constraint[1L],
      sel$.target[1L],
      sel$.achieved[1L]
    ))
  }
  n_lower <- length(opt$active_lower %||% integer(0))
  n_upper <- length(opt$active_upper %||% integer(0))
  if (n_lower > 0L || n_upper > 0L) {
    cat(sprintf("active bounds: %d lower, %d upper\n", n_lower, n_upper))
  }
  if (register) {
    cat("# summary() for the certainty split by stratum\n")
  }
  invisible(x)
}

#' Print generalized allocation precision
#' @keywords internal
#' @noRd
.print_bethel_prec <- function(x) {
  d <- x$detail
  cat(sprintf("Joint allocation precision (%d constraints)\n", nrow(d)))
  cat(sprintf(
    "targets: %s\n",
    if (all(d$.pass)) {
      "all pass"
    } else {
      sprintf("%d violated", sum(!d$.pass))
    }
  ))
  violated <- any(!d$.pass)
  sel <- if (violated) {
    d[!d$.pass, , drop = FALSE]
  } else {
    d[d$.binding, , drop = FALSE]
  }
  .print_constraint_rows(
    sel,
    total = nrow(d),
    label = if (violated) "violated" else "binding",
    hint = "... see $detail for all rows"
  )
  if (!is.null(x$objective_value)) {
    cat(sprintf(
      "objective: weighted relative variance %.6g\n",
      x$objective_value
    ))
    if (!is.null(x$params[["budget"]])) {
      cat(sprintf(
        "budget: %.6g, residual %.6g\n",
        x$params[["budget"]],
        x$params[["budget_residual"]]
      ))
    }
  }
  if (!is.null(x$bounds) && any(!x$bounds$.pass)) {
    cat(sprintf(
      "allocation bounds: %d violated (see $bounds)\n",
      sum(!x$bounds$.pass)
    ))
  }
  invisible(x)
}


#' @rdname print.svyplan
#' @export
format.svyplan_strata <- function(x, ...) {
  .check_unused_dots(...)
  paste0(
    "svyplan_strata [",
    x$method,
    ", ",
    x$n_strata,
    " strata, n=",
    .whole_count(x$n),
    "]"
  )
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_strata <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  as.data.frame(x$strata, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_strata <- function(x, ...) {
  .check_unused_dots(...)
  .coerce_count(.whole_count(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_strata <- function(x, ...) {
  .check_unused_dots(...)
  as.double(x$n)
}

#' Assign observations to strata
#'
#' Apply strata boundaries from a [strata_bound()] result to a numeric
#' vector, returning a factor of stratum assignments.
#'
#' @param object A `svyplan_strata` object.
#' @param newdata Numeric vector to stratify.
#' @param labels Labels for the resulting factor levels. Default `NULL`
#'   generates labels from the boundary intervals.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @return A factor of stratum assignments with length equal to
#'   `length(newdata)`. Values beyond the original training range are
#'   assigned to the lowest or highest stratum.
#'
#' @seealso [strata_bound()] to compute the boundaries.
#'
#' @examples
#' set.seed(1103)
#' x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
#' sb <- strata_bound(x, n_strata = 4, n = 100)
#'
#' # Default interval labels
#' head(predict(sb, x))
#'
#' # Custom labels
#' predict(sb, x, labels = c("Low", "Mid-Low", "Mid-High", "High"))
#'
#' @export
predict.svyplan_strata <- function(object, newdata, labels = NULL, ...) {
  .check_unused_dots(...)
  if (!is.numeric(newdata)) {
    stop("'newdata' must be a numeric vector", call. = FALSE)
  }
  breaks <- c(-Inf, object$boundaries, Inf)
  if (!is.null(labels) && length(labels) != object$n_strata) {
    stop(
      sprintf(
        "'labels' must have length %d (number of strata)",
        object$n_strata
      ),
      call. = FALSE
    )
  }
  out <- cut(newdata, breaks = breaks, include.lowest = TRUE, labels = labels)
  take_all_above <- object$params$take_all_above
  if (!is.null(take_all_above)) {
    # cut() uses right-closed intervals, whereas the take-all contract is
    # x >= take_all_above.  Preserve the training assignment at equality.
    out[newdata >= take_all_above] <- levels(out)[object$n_strata]
  }
  out
}

#' @rdname print.svyplan
#' @export
print.svyplan_twophase <- function(x, ...) {
  .check_unused_dots(...)
  cat("Two-phase allocation (", nrow(x$detail), " phase-2 strata)\n", sep = "")
  cat(.fmt_twophase_sizes(x))
  cat(.fmt_twophase_responding(x))
  cat(.fmt_twophase_cost(x))
  cat(.fmt_twophase_deff(x))
  cat("\n")
  print(.fmt_twophase_detail(x, brief = TRUE), row.names = FALSE)
  cat("\n")
  cat(.fmt_twophase_assured(x))
  cat(.fmt_twophase_search(x))
  cat(.fmt_twophase_single(x))
  cat("# summary() for the continuous optimum and the comparator\n")
  invisible(x)
}

#' The design as it would be fielded
#'
#' Every count printed here is a whole unit off `$operational`, because the
#' continuous solution and the fielded one differ by a unit or two and one
#' printed block must not carry both readings of the same design. The
#' continuous values stay in the object, which is what keeps the round trip
#' exact, and `summary()` is where they are read.
#' @keywords internal
#' @noRd
.twophase_shown <- function(x) {
  o <- x$operational
  if (is.null(o)) {
    return(list(
      n = round(x$n),
      n_int = round(x$detail$n_issued),
      cost = x$cost,
      cv = x$cv
    ))
  }
  list(n = o$n, n_int = o$n_int, cost = o$cost, cv = o$cv)
}

#' @keywords internal
#' @noRd
.fmt_twophase_sizes <- function(x) {
  shown <- .twophase_shown(x)
  sprintf(
    "field design: n_phase1 = %s | n_phase2 = %s\n",
    format(shown$n[["n_phase1"]]),
    format(shown$n[["n_phase2"]])
  )
}

#' What the two phases are expected to leave
#'
#' On the fielded sizes, not the continuous ones: a responding count read off
#' a design a planner is not releasing is a number nothing in the block
#' reconciles with.
#' @keywords internal
#' @noRd
.fmt_twophase_responding <- function(x) {
  if (isTRUE(all.equal(unname(x$responding), unname(x$n)))) {
    return("")
  }
  shown <- .twophase_shown(x)
  sprintf(
    "expected responding: n_phase1 = %s | n_phase2 = %s\n",
    format(round(shown$n[["n_phase1"]] *
                   x$responding[["n_phase1"]] / x$n[["n_phase1"]])),
    format(round(sum(shown$n_int * x$detail$resp_rate)))
  )
}

#' @keywords internal
#' @noRd
.fmt_twophase_cost <- function(x) {
  shown <- .twophase_shown(x)
  fixed <- x$params$fixed_cost
  sprintf(
    "cv = %s, cost = %s%s\n",
    if (is.na(shown$cv)) {
      "undefined without 'mu' or a 'mean' column"
    } else {
      formatC(shown$cv, format = "f", digits = 4)
    },
    format(round(shown$cost)),
    if (isTRUE(fixed > 0)) sprintf(" (fixed: %s)", format(fixed)) else ""
  )
}

#' @keywords internal
#' @noRd
.fmt_twophase_deff <- function(x) {
  p <- x$params
  if (
    isTRUE(all.equal(p$phase1_deff, 1)) &&
      isTRUE(all.equal(p$single_deff, 1))
  ) {
    return("")
  }
  sprintf(
    "phase-1 deff = %s, single-phase deff = %s\n",
    formatC(p$phase1_deff, format = "f", digits = 2),
    formatC(p$single_deff, format = "f", digits = 2)
  )
}

#' The stratum table at two widths
#'
#' `brief` is what `print()` shows. The continuous `n_issued` is dropped
#' there: it and `n_int` differ by less than a unit, so side by side they
#' read as two answers to one question.
#' @keywords internal
#' @noRd
.fmt_twophase_detail <- function(x, brief = FALSE) {
  d <- x$detail
  shown <- .twophase_shown(x)
  out <- data.frame(
    stratum = d$stratum,
    share = formatC(d$share, format = "f", digits = 3),
    sd = formatC(d$sd, format = "f", digits = 2),
    unit_cost = formatC(d$unit_cost, format = "f", digits = 2),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(all.equal(d$deff, rep(1, nrow(d))))) {
    out$deff <- formatC(d$deff, format = "f", digits = 2)
  }
  show_resp <- !isTRUE(all.equal(d$resp_rate, rep(1, nrow(d))))
  if (show_resp) {
    out$resp <- formatC(d$resp_rate, format = "f", digits = 2)
  }
  # The print reads the field design, so its fraction is the whole take
  # over the whole pool phase 1 classifies, on nu's scale (a share of the
  # units issued, at most resp_rate). summary() keeps the continuous one
  # beside the continuous issue.
  nu <- if (brief) {
    r1 <- x$params$resp_rate %||% 1
    pool <- floor(r1 * d$share * shown$n[["n_phase1"]] + 1e-9)
    ifelse(pool > 0, r1 * shown$n_int / pool, 0)
  } else {
    d$nu
  }
  out$nu <- formatC(nu, format = "f", digits = 4)
  if (!brief) {
    out$n_issued <- format(round(d$n_issued))
  }
  out$n_int <- format(shown$n_int)
  if (show_resp) {
    out$n_resp <- format(round(shown$n_int * d$resp_rate))
  }
  if (any(d$take_all)) {
    out$take_all <- ifelse(d$take_all, "*", "")
  }
  out
}

#' @keywords internal
#' @noRd
.fmt_twophase_assured <- function(x) {
  o <- x$operational
  if (is.null(o) || is.null(o$assured)) {
    return("")
  }
  sprintf(
    "assured (%s): issue n_phase1 = %s | n_phase2 = %s (cost %s)\n",
    .fmt_prob(x$params$assurance),
    format(o$assured_phase1),
    format(sum(o$assured)),
    format(round(o$assured_cost))
  )
}

#' Analyze a two-phase allocation
#'
#' `print()` on a [n_twophase()] result gives the design as it would be
#' fielded: the two phase sizes in whole units, the precision and cost they
#' buy, the per-stratum subsampling fractions, and which of the two designs
#' to run. Its `nu` is the fielded fraction, each stratum's `n_int` over the
#' whole pool phase 1 classifies, times `resp_rate`, so a stratum that takes
#' its whole pool reads `resp_rate`. `summary()` gives what that answer was
#' chosen against: the continuous optimum beside the fielded one, the full
#' stratum table with the continuous `nu` and issue, the single-phase
#' comparator in whole units and at its continuous optimum, whether it
#' reaches the target, and the planning assumptions.
#'
#' @param object A `svyplan_twophase` object.
#' @param x A `summary.svyplan_twophase` object.
#' @param ... Additional arguments are not supported and produce an error.
#' @return `summary()` returns an object of class
#'   `summary.svyplan_twophase`. Its `print()` method returns it invisibly.
#' @seealso [n_twophase()] for the planner, and [print.svyplan] for the
#'   printed block this expands on.
#'
#' @examples
#' frame <- data.frame(
#'   stratum   = c("A", "B", "C", "D"),
#'   N         = c(3500, 2500, 2500, 1500),
#'   sd        = c(12, 25, 8, 40),
#'   mean      = c(40, 70, 35, 90),
#'   unit_cost = c(2, 5, 1, 9)
#' )
#' plan <- n_twophase(frame, phase1_cost = 1, budget = 50000)
#' plan
#' summary(plan)
#' summary(plan)$detail
#'
#' @name summary.svyplan_twophase
NULL

#' @rdname summary.svyplan_twophase
#' @export
summary.svyplan_twophase <- function(object, ...) {
  .check_unused_dots(...)
  structure(
    list(
      plan = object,
      n = object$n,
      operational = .twophase_shown(object),
      detail = .fmt_twophase_detail(object, brief = FALSE),
      single_phase = object$single_phase,
      params = object$params
    ),
    class = "summary.svyplan_twophase"
  )
}

#' @rdname summary.svyplan_twophase
#' @export
print.summary.svyplan_twophase <- function(x, ...) {
  .check_unused_dots(...)
  plan <- x$plan
  cat(sprintf(
    "Analysis of a two-phase allocation (%d phase-2 strata)\n\n",
    nrow(plan$detail)
  ))
  cat(.fmt_twophase_sizes(plan))
  cat(.fmt_twophase_responding(plan))
  cat(.fmt_twophase_cost(plan))
  cat(sprintf(
    "continuous optimum: n_phase1 = %s | n_phase2 = %s (cv %s, cost %s)\n",
    format(round(x$n[["n_phase1"]])),
    format(round(x$n[["n_phase2"]])),
    if (is.na(plan$cv)) "NA" else formatC(plan$cv, format = "f", digits = 4),
    format(round(plan$cost))
  ))
  cat(.fmt_twophase_deff(plan))
  cat(.fmt_twophase_assured(plan))
  cat(.fmt_twophase_search(plan))

  cat("\nStrata\n")
  print(x$detail, row.names = FALSE)

  s <- x$single_phase
  if (!is.null(s) && is.finite(s$n)) {
    fmt_cv <- function(v) {
      if (is.na(v)) "NA" else formatC(v, format = "f", digits = 4)
    }
    cat(sprintf(
      "\nSingle-phase comparator: n = %s, cv = %s, cost = %s\n",
      format(s$n_int %||% round(s$n)),
      fmt_cv(s$cv_int %||% s$cv),
      format(round(s$cost_int %||% s$cost))
    ))
    cat(sprintf(
      "continuous optimum: n = %s (cv %s, cost %s)\n",
      format(signif(s$n, 6)), fmt_cv(s$cv), format(round(s$cost))
    ))
    cat(sprintf(
      "%s, and it %s the target\n",
      if (isTRUE(s$better)) {
        "It is the better design here"
      } else {
        "Two-phase is the better design here"
      },
      if (isTRUE(s$reaches_target)) "reaches" else "does not reach"
    ))
  }
  invisible(x)
}

#' The search line, only when a limit stopped the whole-unit search
#' @keywords internal
#' @noRd
.fmt_twophase_search <- function(x) {
  search <- x$operational$search
  if (!identical(search$status, "limit")) {
    return("")
  }
  which <- paste(search$limit, collapse = " and ")
  if (is.finite(search$gap)) {
    sprintf("search: stopped at the %s limit, gap %.3g%% to a bound\n",
            which, 100 * search$gap)
  } else {
    sprintf("search: stopped at the %s limit, the best design found\n",
            which)
  }
}

#' The single-phase comparator, as a verdict
#'
#' Which design to field is the whole content, so the line states it and the
#' numbers that decide it. `summary()` carries the comparator's cost and the
#' target it does or does not reach.
#' @keywords internal
#' @noRd
.fmt_twophase_single <- function(x) {
  s <- x$single_phase
  if (is.null(s) || !is.finite(s$n)) {
    return("")
  }
  cv_int <- s$cv_int %||% s$cv
  cv <- if (is.na(cv_int)) "NA" else formatC(cv_int, format = "f", digits = 4)
  n_int <- s$n_int %||% round(s$n)
  if (isTRUE(s$better)) {
    sprintf(
      "single-phase is better here: n = %s at cv %s, so skip phase 1\n",
      format(n_int),
      cv
    )
  } else {
    sprintf(
      "single-phase alternative: n = %s at cv %s, two-phase wins\n",
      format(n_int),
      cv
    )
  }
}

#' Print and coerce design degrees of freedom
#'
#' Display and coercion methods for the count that [design_df()] returns.
#' `print()` shows the headline count and the numbers it was formed from,
#' while `summary()` gives the additive per-stratum decomposition and the
#' separate per-domain counts. The coercion and arithmetic methods let the
#' object stand in for the headline count wherever a plain number is expected.
#'
#' @param x A `svyplan_df` object from [design_df()], or its summary for the
#'   summary print method.
#' @param object A `svyplan_df` object to summarize.
#' @param value Replacement value. Replacement is refused, the count being
#'   derived from the strata and domains it carries.
#' @param row.names,optional Standard [as.data.frame()] arguments.
#' @param stringsAsFactors Logical. Retained for compatibility when a result
#'   is converted through [data.frame()].
#' @param validRN Logical. Accepted for compatibility with [data.frame()] in
#'   R 4.7.0 and later. Svyplan results already have valid row names.
#' @param e1,e2 Objects supplied to an arithmetic or comparison operator.
#' @param name,i A field name, one of those listed under Details. `[[` also
#'   accepts a numeric index, which reads the underlying numeric vector.
#' @param ... For mathematical transformations, additional arguments passed
#'   to the underlying operation. The other methods do not support
#'   additional arguments.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns a
#'   `summary.svyplan_df` object containing the overall count, its inputs,
#'   the per-stratum decomposition when available, the per-domain counts,
#'   and the counting basis. `format()` returns a character scalar, and
#'   `as.double()` returns the degrees of freedom.
#'   `as.data.frame()` returns a one-row table of the scalar fields.
#'   `as.list()` returns every field, including the per-stratum and
#'   per-domain tables, and `$` and `[[` return one of them. Arithmetic and
#'   mathematical transformations return ordinary numeric results.
#'
#' @details
#' A `svyplan_df` behaves as the numeric degrees of freedom wherever one is
#' expected. It can be passed to any `df` argument, compared, and
#' arithmetically combined, with the detail dropped by any such operation.
#'
#' The fields are `df` for the count itself, `n_units` for the units it
#' counts and `stage` for what those units are (`"psu"` or `"element"`),
#' `n_strata` for the constraints subtracted, and the `strata` and
#' `domains` tables, which are `NULL` for an unstratified or domain-free
#' plan. Naming a field the object does not have is an error listing the
#' ones it does.
#'
#' The summary's stratum table is an additive decomposition. A contributing
#' stratum uses one constraint and contributes its sampled units minus one.
#' A census stratum has no units counted and no constraint, while a singleton
#' uses its one constraint and contributes zero. The domain table is printed
#' separately because domain degrees of freedom are alternative subpopulation
#' counts, not additional contributions to the overall design.
#'
#' @seealso [design_df()], which builds these objects, and
#'   [print.svyplan_deff] for the design effect's counterpart.
#'
#' @examples
#' d <- design_df(n_psu = 300, n_strata = 20)
#' d
#' d$df
#' as.double(d) + 1
#' as.data.frame(d)
#' summary(d)
#'
#' @name print.svyplan_df
NULL

#' @rdname print.svyplan_df
#' @export
print.svyplan_df <- function(x, ...) {
  .check_unused_dots(...)
  unit_label <- if (identical(attr(x, "stage", exact = TRUE), "psu")) {
    "PSUs"
  } else {
    "units"
  }
  n_units <- attr(x, "n_units", exact = TRUE)
  n_strata <- attr(x, "n_strata", exact = TRUE)
  strata <- attr(x, "strata", exact = TRUE)
  constraint_label <- if (!is.null(strata$n_groups)) {
    if (n_strata == 1L) "variance group" else "variance groups"
  } else {
    if (n_strata == 1L) "stratum" else "strata"
  }
  cat("Design degrees of freedom (planning)\n\n")
  if (identical(attr(x, "method", exact = TRUE), "satterthwaite")) {
    cat(sprintf(
      "  df = %.4g   (Satterthwaite, nominal %g = %g %s - %d %s)\n",
      as.double(x),
      attr(x, "nominal", exact = TRUE),
      n_units,
      unit_label,
      n_strata,
      constraint_label
    ))
  } else {
    cat(sprintf(
      "  df = %g   (%g %s - %d %s)\n",
      as.double(x),
      n_units,
      unit_label,
      n_strata,
      constraint_label
    ))
  }
  if (!is.null(strata)) {
    flagged <- strata$.status != "ok"
    if (any(flagged)) {
      cat(sprintf(
        "  no df from %s: %s\n",
        if (sum(flagged) > 1L) "these strata" else "this stratum",
        paste(
          sprintf("%s (%s)", strata$stratum[flagged], strata$.status[flagged]),
          collapse = ", "
        )
      ))
    }
  }
  domains <- attr(x, "domains", exact = TRUE)
  if (!is.null(domains)) {
    cat(sprintf(
      "  %d domains, df from %g to %g\n",
      nrow(domains),
      min(domains$.df),
      max(domains$.df)
    ))
  }
  invisible(x)
}

#' @rdname print.svyplan_df
#' @export
summary.svyplan_df <- function(object, ...) {
  .check_unused_dots(...)
  stage <- attr(object, "stage", exact = TRUE)
  unit_label <- if (identical(stage, "psu")) "PSUs" else "units"
  method <- attr(object, "method", exact = TRUE) %||% "nominal"
  nominal <- attr(object, "nominal", exact = TRUE)
  structure(
    list(
      df = as.double(object),
      n_units = attr(object, "n_units", exact = TRUE),
      n_strata = attr(object, "n_strata", exact = TRUE),
      stage = stage,
      strata = attr(object, "strata", exact = TRUE),
      domains = attr(object, "domains", exact = TRUE),
      basis = if (identical(method, "satterthwaite")) {
        sprintf(
          "Satterthwaite, combining the per-stratum df by each stratum's share of the variance. The per-stratum df add to the nominal %s, counted %s minus contributing strata",
          format(nominal), unit_label
        )
      } else if (is.null(attr(object, "strata", exact = TRUE)$n_groups)) {
        sprintf("counted %s minus contributing strata", unit_label)
      } else {
        sprintf("counted %s minus variance groups (zones, or collapsed zones)",
                unit_label)
      },
      method = method,
      nominal = nominal
    ),
    class = "summary.svyplan_df"
  )
}

#' @rdname print.svyplan_df
#' @export
print.summary.svyplan_df <- function(x, ...) {
  .check_unused_dots(...)
  cat("Analysis of design degrees of freedom\n\n")

  if (is.null(x$strata)) {
    tab <- data.frame(
      Units = x$n_units,
      Constraints = x$n_strata,
      `Design df` = x$df,
      check.names = FALSE,
      row.names = "Overall"
    )
    print(tab, right = TRUE)
    cat("\nPer-stratum counts were not supplied.\n")
  } else {
    status <- x$strata$.status
    census <- status == "census"
    counted <- ifelse(census, 0, x$strata$n_units)
    constraints <- x$strata$n_groups %||% as.integer(!census)
    tab <- if (any(census)) {
      data.frame(
        Stratum = x$strata$stratum,
        Sampled = x$strata$n_units,
        Counted = counted,
        Constraints = constraints,
        `Design df` = x$strata$df,
        Status = status,
        check.names = FALSE
      )
    } else {
      data.frame(
        Stratum = x$strata$stratum,
        Units = counted,
        Constraints = constraints,
        `Design df` = x$strata$df,
        Status = status,
        check.names = FALSE
      )
    }
    if (identical(x$method, "satterthwaite") && !is.null(x$strata$.share)) {
      tab[["Variance share"]] <- round(x$strata$.share, 4)
    }
    overall <- tab[1L, , drop = FALSE]
    overall[1L, ] <- NA
    overall[["Stratum"]] <- "Overall"
    if ("Sampled" %in% names(overall)) {
      overall[["Sampled"]] <- sum(x$strata$n_units)
      overall[["Counted"]] <- x$n_units
    } else {
      overall[["Units"]] <- x$n_units
    }
    overall[["Constraints"]] <- x$n_strata
    overall[["Design df"]] <- x$df
    overall[["Status"]] <- ""
    if ("Variance share" %in% names(overall)) {
      overall[["Variance share"]] <- sum(x$strata$.share, na.rm = TRUE)
    }
    print(rbind(tab, overall), row.names = FALSE, quote = FALSE, right = TRUE)
  }

  cat(sprintf(
    "\nCounted stage: %s\n",
    if (identical(x$stage, "psu")) {
      "PSU"
    } else {
      "element"
    }
  ))
  cat(sprintf("Basis: %s\n", x$basis))

  if (!is.null(x$domains)) {
    domains <- x$domains
    id_cols <- setdiff(
      names(domains),
      c(".domain", ".n_units", ".df", ".status")
    )
    if (length(id_cols) > 0L) {
      domains$.domain <- NULL
    }
    names(domains)[names(domains) == ".domain"] <- "Domain"
    names(domains)[names(domains) == ".n_units"] <- "Units"
    names(domains)[names(domains) == ".df"] <- "Design df"
    names(domains)[names(domains) == ".status"] <- "Status"
    cat("\nDomain degrees of freedom (not additive to the overall row)\n\n")
    print(domains, row.names = FALSE, quote = FALSE, right = TRUE)
  }
  invisible(x)
}

#' @rdname print.svyplan_df
#' @export
format.svyplan_df <- function(x, ...) {
  .check_unused_dots(...)
  sprintf("svyplan_df [%s, %g]", attr(x, "stage", exact = TRUE), as.double(x))
}

#' @rdname print.svyplan_df
#' @export
as.double.svyplan_df <- function(x, ...) {
  .check_unused_dots(...)
  unclass(x)[[1L]]
}

#' @rdname print.svyplan_df
#' @export
as.list.svyplan_df <- function(x, ...) {
  .check_unused_dots(...)
  list(
    df = as.double(x),
    n_units = attr(x, "n_units", exact = TRUE),
    n_strata = attr(x, "n_strata", exact = TRUE),
    stage = attr(x, "stage", exact = TRUE),
    strata = attr(x, "strata", exact = TRUE),
    domains = attr(x, "domains", exact = TRUE),
    method = attr(x, "method", exact = TRUE) %||% "nominal",
    nominal = attr(x, "nominal", exact = TRUE)
  )
}

#' @rdname print.svyplan_df
#' @export
`$.svyplan_df` <- function(x, name) {
  .df_field(x, name)
}

#' @rdname print.svyplan_df
#' @export
`[[.svyplan_df` <- function(x, i, ...) {
  if (is.character(i)) {
    return(.df_field(x, i))
  }
  unclass(x)[[i, ...]]
}

#' Look a degrees-of-freedom field up by name
#' @keywords internal
#' @noRd
.df_field <- function(x, name) {
  fields <- as.list(x)
  if (length(name) != 1L || !name %in% names(fields)) {
    stop(
      sprintf(
        "no field '%s' in a design df; available: %s",
        paste(name, collapse = ", "),
        paste(names(fields), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  fields[[name]]
}

#' @rdname print.svyplan_df
#' @export
as.data.frame.svyplan_df <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  # The two tables are data frames in their own right and are read from
  # $strata and $domains; a one-row export carries the scalars only.
  out <- data.frame(
    df = as.double(x),
    n_units = attr(x, "n_units", exact = TRUE),
    n_strata = attr(x, "n_strata", exact = TRUE),
    stage = attr(x, "stage", exact = TRUE),
    method = attr(x, "method", exact = TRUE) %||% "nominal",
    nominal = attr(x, "nominal", exact = TRUE) %||% NA_real_,
    stringsAsFactors = stringsAsFactors
  )
  as.data.frame(out, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan_df
#' @export
Ops.svyplan_df <- function(e1, e2) {
  e1 <- if (inherits(e1, "svyplan_df")) as.double(e1) else e1
  if (missing(e2)) {
    return(do.call(.Generic, list(e1)))
  }
  e2 <- if (inherits(e2, "svyplan_df")) as.double(e2) else e2
  do.call(.Generic, list(e1, e2))
}

#' @rdname print.svyplan_df
#' @export
Math.svyplan_df <- function(x, ...) {
  do.call(.Generic, c(list(as.double(x)), list(...)))
}

