#' Print svyplan objects
#'
#' Display and coercion methods shared by every result the package returns.
#' `print()` gives the human-readable design, while `as.integer()`,
#' `as.double()`, and `as.data.frame()` extract it in the shapes a
#' downstream package can consume.
#'
#' @param x A svyplan object.
#' @param object A svyplan object, for `summary()`.
#' @param row.names,optional Standard [as.data.frame()] arguments.
#' @param stringsAsFactors Logical. Retained for compatibility when a result
#'   is converted through [data.frame()].
#' @param validRN Logical. Accepted for compatibility with [data.frame()] in
#'   R 4.7.0 and later. Svyplan results already have valid row names.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @return `print()` returns `x` invisibly. `format()` returns a character
#'   vector, `as.integer()` and `as.double()` return numeric vectors, and
#'   `as.data.frame()` returns a data frame. The shapes are described under
#'   Details. `summary()` of a size or a cluster design is described in
#'   [summary.svyplan_n]. A result type with no structured summary yet
#'   returns itself from `summary()`, so it prints as `print()` does.
#'
#' @details
#' ## Print and coercion
#'
#' Constrained designs (`n_cluster()`, `n_alloc()`)
#' carry two representations: the continuous design in the top-level
#' fields (`n`, `cv`, `cost`, ...) and the whole-unit field design in
#' `$operational`, whose cost and precision are recomputed from the
#' integer design. `print()` shows the field design. For `n_cluster()` the
#' continuous optimum is read in `summary()` (see [summary.svyplan_n]). For
#' `n_alloc()` the print also shows the allocation rule's continuous design
#' as a diagnostic, which a two-stage frame with a finite population
#' correction can improve on (see [n_alloc()]).
#'
#' `as.integer(x)` returns the operational design in the same shape as
#' `x$n`: the named integer stage vector for `svyplan_cluster` objects,
#' the operational total for allocation results, the sum of the ceiled
#' quotas for `n_multi()` with separate domain quotas, and the ceiled scalar
#' otherwise. A count past R's integer range has no integer form, and
#' `as.integer(x)` then stops with an error of class `svyplan_count_range`.
#' `x$operational$n` keeps such a field count as a double. `as.double(x)`
#' returns the continuous counterpart of the same shape. For `svyplan_strata`, both coercions return the total sample
#' size. Boundary cutpoints remain available in `$boundaries`.
#'
#' `as.data.frame()` returns the tabular form of a result, intended as
#' the stable handoff to downstream packages. For
#' `svyplan_n`: the stratum allocation table (`$detail`) for
#' `n_alloc()` results, the per-domain table (`$domains`, falling back
#' to `$detail`) for `n_multi()` results, and a one-row summary
#' (`n`, `n_int`, `se`, `moe`, `cv`, `rmoe`) otherwise, with an `occasion`
#' column and one row per occasion for a two-occasion `n_change()`. For
#' `svyplan_cluster`:
#' the per-domain table when domains are present, otherwise a stage
#' table with columns `stage`, `n`, and `n_int`, where `n` is the continuous
#' optimum and `n_int` is the constraint-preserving operational design.
#' For `svyplan_prec`, detail-bearing multi-indicator and allocation results
#' return `$detail`. Single-indicator results return their sample size and
#' precision measures in one row. For `svyplan_varcomp`, stratified results
#' return `$strata`, while unstratified results return one row containing the
#' two- or three-stage components. For `svyplan_power`, one row contains the
#' two group sizes, their integer counterparts, power, effect, type, and the
#' quantity that was solved for.
#'
#' @seealso [confint.svyplan] for confidence intervals on a result.
#'
#' @examples
#' # print leads with the operational (whole-unit) design
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
#'
#' # the tabular handoff to downstream packages
#' as.data.frame(n_prop(p = 0.3, moe = 0.05))
#'
#' @name print.svyplan
NULL

#' Name the stage correction a cluster result carries
#' @keywords internal
#' @noRd
.print_stage_fpc <- function(p) {
  if (!identical(p[["fpc"]], "stage")) {
    return(invisible(NULL))
  }
  N_vec <- p[["N"]]
  if (is.null(N_vec)) {
    cat(sprintf(
      "corrected by stage for %s in 'indicators'\n",
      if (length(p[["domain_cols"]]) > 0L) "each domain's population" else
        "the population"
    ))
  } else if (length(N_vec) == length(p[["stage_cost"]])) {
    nms <- .stage_names(length(N_vec), names(p[["stage_cost"]]))
    cat(sprintf(
      "corrected by stage for %s\n",
      paste(nms, "=", vapply(N_vec, .fmt_count_n, character(1L)),
            collapse = " | ")
    ))
  } else {
    cat(sprintf(
      "corrected by stage for N = c(%s)\n",
      paste(vapply(N_vec, .fmt_count_n, character(1L)), collapse = ", ")
    ))
  }
  invisible(NULL)
}

#' Name the stages a cluster design loses units at
#'
#' Three rates acting at three stages are not interchangeable, so the print
#' says which one is set rather than showing a single netted figure.
#' Share of issued units expected to respond, response times eligibility
#' @keywords internal
#' @noRd
.net_rate <- function(p) {
  (p[["resp_rate"]] %||% 1) * (p[["eligibility_rate"]] %||% 1)
}

#' @keywords internal
#' @noRd
.print_stage_rates <- function(p) {
  rr <- p[["resp_rate"]] %||% 1
  parts <- character(0)
  if (length(rr) == 1L) {
    if (rr < 1) {
      parts <- c(parts, sprintf("resp_rate = %s", .fmt_prob(rr)))
    }
  } else if (any(rr < 1)) {
    if (all(rr[-length(rr)] == 1)) {
      parts <- c(parts, sprintf("resp_rate = %s", .fmt_prob(rr[length(rr)])))
    } else {
      parts <- c(parts, sprintf(
        "resp_rate = c(%s)",
        paste(vapply(rr, .fmt_prob, ""), collapse = ", ")
      ))
    }
  }
  if (!is.null(p[["eligibility_rate"]])) {
    parts <- c(parts, sprintf("eligibility_rate = %s",
                              .fmt_prob(p[["eligibility_rate"]])))
  }
  if (!is.null(p[["eligible_per_unit"]])) {
    parts <- c(parts, sprintf("eligible_per_unit = %.4g", p[["eligible_per_unit"]]))
  }
  if (length(parts) == 0L) {
    return(invisible(NULL))
  }
  cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))
  invisible(NULL)
}

#' Format the design effect for printing
#'
#' `deff = 1` is the most consequential silent default in survey planning,
#' so it is always shown rather than suppressed when left at its default.
#' It is printed without decimals so that an assumed value is visually
#' distinct from a supplied one.
#' @keywords internal
#' @noRd
.fmt_deff <- function(deff) {
  if (is.null(deff)) {
    return(NULL)
  }
  if (isTRUE(all.equal(unname(deff), rep(1, length(deff))))) {
    "deff = 1"
  } else {
    sprintf(
      "deff = %s",
      .fmt_by_stratum(deff, function(v) sprintf("%.2f", v))
    )
  }
}

#' Format a probability without letting it round onto a bound
#'
#' Two decimals everywhere it already reads well, and more where rounding
#' would print a value the validators refuse or the surrounding condition
#' denies: an assurance of 0.999 shown as 1.00 names a level `n_panel()` and
#' `n_twophase()` both reject, and a response rate of 0.999 shown as 1.00
#' contradicts the test that made the line appear at all.
#' @keywords internal
#' @noRd
.fmt_prob <- function(p) {
  out <- sprintf("%.2f", p)
  # Enough significant digits that nothing inside (0, 1) can print as a
  # bound: "%g" alone rounds 0.9999999 to 1.
  if (out %in% c("0.00", "1.00") && p > 0 && p < 1) {
    return(sprintf("%.15g", p))
  }
  out
}

#' Format a design parameter that may vary by stratum
#'
#' A single value prints as itself. Values that differ print as their range,
#' so the header stays one line however many strata there are. `fmt`
#' formats one value.
#' @keywords internal
#' @noRd
.fmt_by_stratum <- function(x, fmt) {
  u <- unique(x)
  if (length(u) == 1L) {
    return(fmt(u))
  }
  sprintf("%s to %s by stratum", fmt(min(x)), fmt(max(x)))
}

#' @rdname print.svyplan
#' @export
print.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  if (x$type == "multi") {
    .print_multi_n(x)
  } else if (x$type == "change") {
    .print_change_n(x)
  } else if (x$type == "pooled") {
    .print_pooled_n(x)
  } else {
    .print_single_n(x)
  }
  invisible(x)
}

#' @rdname print.svyplan
#' @export
print.svyplan_alloc <- function(x, ...) {
  .check_unused_dots(...)
  .print_alloc_n(x)
  invisible(x)
}

#' Report a size that buys one estimate averaged over several occasions
#'
#' The size is per occasion. The distinct units the series consumes is
#' smaller than `occasions * n` whenever the occasions overlap, and the
#' number of interviews is not, so the header states which quantity it is.
#' @keywords internal
#' @noRd
.print_pooled_n <- function(x) {
  p <- x$params
  cat(sprintf("Sample size for pooled estimate (%s)\n", .pooled_scale_label(p)))
  cat(sprintf("n = %s per occasion", .fmt_count_n(.whole_count(x$n))))
  if (.net_rate(p) < 1) {
    cat(sprintf(" (net: %s)", .fmt_count_n(.whole_count(x$n * .net_rate(p)))))
  }
  cat(sprintf(", %d occasions", p$occasions))
  parts <- c(.fmt_pooled_estimand(p), .fmt_deff(p$deff))
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_lag_overlap(p))
  cat(sprintf("se = %.4g, moe = %.4g", x$se, x$moe))
  if (!is.na(x$cv)) {
    cat(sprintf(", cv = %.4g", x$cv))
  }
  if (!is.null(x$rmoe) && !is.na(x$rmoe)) {
    cat(sprintf(", rmoe = %.4g", x$rmoe))
  }
  cat("\n")
}

#' The scale a change or pooled estimate is stated on
#'
#' A ratio is checked first, and with `[[`, since `$r` would partial-match
#' `resp_rate` on every result that has no ratio.
#' @keywords internal
#' @noRd
.repeated_scale_label <- function(p) {
  if (!is.null(p[["r"]])) {
    "ratio of totals"
  } else if (!is.null(p[["p"]])) {
    "proportion scale"
  } else {
    "mean scale"
  }
}

#' @keywords internal
#' @noRd
.pooled_scale_label <- function(p) {
  .repeated_scale_label(p)
}

#' @keywords internal
#' @noRd
.fmt_pooled_estimand <- function(p) {
  if (!is.null(p[["r"]])) {
    return(sprintf("r = %.4g", p[["r"]]))
  }
  if (!is.null(p[["p"]])) {
    return(sprintf("p = %.3g", p[["p"]]))
  }
  c(
    sprintf("var = %.4g", p$var),
    if (!is.null(p$mu)) sprintf("mu = %.4g", p$mu)
  )
}

#' Describe a lag profile without printing every entry of it
#'
#' A twelve-occasion average carries eleven overlaps and eleven
#' correlations, which is a table rather than a line. The consecutive figure
#' is the one a planner states, so it leads, and the rest is summarized by
#' how far the shared units reach.
#' @keywords internal
#' @noRd
.fmt_lag_overlap <- function(p) {
  ov <- p$overlap %||% 0
  rho <- p$overlap_cor %||% 0
  rotation_term <- .pooled_coupled(p[["overlap_rotation"]],
                                   p$overlap_draw %||% "fixed") &&
    all(is.finite(p$N %||% Inf)) && any(rho[ov == 0] > 0)
  if (all(ov * rho == 0) && !rotation_term) {
    return(sprintf(
      "No between-occasion covariance (overlap x overlap_cor = 0)%s\n",
      .fmt_overlap_draw(p)
    ))
  }
  if (!any(ov > 0)) {
    return(sprintf(
      "overlap = 0 at every lag, overlap_cor = %.3g at lag 1, disjoint cohorts%s\n",
      rho[1L],
      .fmt_overlap_draw(p)
    ))
  }
  reach <- max(which(ov > 0))
  sprintf(
    "overlap = %.3g, overlap_cor = %.3g at lag 1, shared out to lag %d%s\n",
    ov[1L],
    rho[1L],
    reach,
    .fmt_overlap_draw(p)
  )
}

#' Name the overlap draw only where it moves the result
#'
#' Independent draws change the finite population terms alone, so the
#' default prints nothing and an infinite population or a zero overlap,
#' where the two draws agree, prints nothing either.
#' @keywords internal
#' @noRd
.fmt_overlap_draw <- function(p) {
  if (identical(p$overlap_draw, "independent") && any((p$overlap %||% 0) > 0) &&
        all(is.finite(p$N))) {
    ", independent draws"
  } else {
    ""
  }
}

#' Report a size that buys one change measured on two occasions
#'
#' The size is per occasion, and a positive overlap means the occasions are
#' not disjoint, so the header says which it is rather than leaving a
#' reader to total two numbers that partly count the same units.
#' @keywords internal
#' @noRd
.print_change_n <- function(x) {
  p <- x$params
  cat(sprintf("Sample size for change (%s)\n", .change_scale_label(p)))
  cat(.fmt_change_n(x$n, .net_rate(p)))
  parts <- c(
    .fmt_change_estimand(p),
    if (!is.null(p$moe)) sprintf("moe = %.4g", p$moe),
    if (!is.null(p$rmoe)) sprintf("rmoe = %.3f", p$rmoe),
    if (!is.null(p$cv)) sprintf("cv = %.3f", p$cv),
    .fmt_deff(p$deff),
    if (!is.null(p$resp_rate) && p$resp_rate < 1) {
      sprintf("resp_rate = %s", .fmt_prob(p$resp_rate))
    },
    if (!is.null(p[["eligibility_rate"]])) {
      sprintf("eligibility_rate = %s", .fmt_prob(p[["eligibility_rate"]]))
    }
  )
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_overlap(p))
}

#' @keywords internal
#' @noRd
.change_scale_label <- function(p) {
  .repeated_scale_label(p)
}

#' @keywords internal
#' @noRd
.fmt_change_n <- function(n, resp_rate, exact = FALSE) {
  net <- !is.null(resp_rate) && resp_rate < 1
  fmt_n <- if (exact) {
    .evaluated_n_formatter(n)
  } else {
    function(v) .fmt_count_n(.whole_count(v))
  }
  if (length(n) == 1L) {
    out <- sprintf("n = %s per occasion", fmt_n(n))
    if (net) {
      out <- sprintf("%s (net: %s)", out, fmt_n(n * resp_rate))
    }
    return(out)
  }
  out <- sprintf("n = %s then %s", fmt_n(n[1L]), fmt_n(n[2L]))
  if (net) {
    out <- sprintf(
      "%s (net: %s then %s)",
      out,
      fmt_n(n[1L] * resp_rate),
      fmt_n(n[2L] * resp_rate)
    )
  }
  out
}

#' @keywords internal
#' @noRd
.fmt_change_estimand <- function(p) {
  if (!is.null(p[["r"]])) {
    return(sprintf("r = %.4g to %.4g", p[["r"]][1L], p[["r"]][2L]))
  }
  if (!is.null(p[["p"]])) {
    return(sprintf("p = %.3g to %.3g", p[["p"]][1L], p[["p"]][2L]))
  }
  c(
    sprintf("var = %s", paste(sprintf("%.4g", p$var), collapse = ", ")),
    if (!is.null(p$change)) sprintf("change = %.4g", p$change)
  )
}

#' Report the overlap and the correlation as the one product they are
#'
#' Either at zero leaves the change with no covariance term, so the line
#' names that state rather than printing two numbers whose joint meaning a
#' reader has to reconstruct. It says what the model does, not that the
#' occasions are disjoint: a full overlap at zero correlation shares every
#' unit and still contributes nothing.
#' @keywords internal
#' @noRd
.fmt_overlap <- function(p) {
  overlap <- p$overlap %||% 0
  overlap_cor <- p$overlap_cor %||% 0
  if (overlap == 0 || overlap_cor == 0) {
    return(sprintf(
      "No between-occasion covariance (overlap x overlap_cor = 0)%s\n",
      .fmt_overlap_draw(p)
    ))
  }
  sprintf(
    "overlap = %.3g, overlap_cor = %.3g (%.1f%% of the independent variance)%s\n",
    overlap,
    overlap_cor,
    100 * .change_var_share(p),
    .fmt_overlap_draw(p)
  )
}

#' Share of the two-independent-samples variance the overlap leaves
#'
#' Reported without the finite population terms, which cancel at full
#' overlap and would otherwise make the percentage depend on `N` as well as
#' on the schedule the planner controls.
#' @keywords internal
#' @noRd
.change_var_share <- function(p) {
  v <- p$var
  ratio <- p$ratio %||%
    (if (length(p$n) == 2L) p$n[1L] / p$n[2L] else 1)
  base <- v[1L] / ratio + v[2L]
  cross <- 2 * p$overlap * p$overlap_cor * sqrt(v[1L] * v[2L])
  if (base <= 0) {
    return(NA_real_)
  }
  max(0, base - cross) / base
}

#' @keywords internal
#' @noRd
.print_single_n <- function(x) {
  type_label <- switch(x$type, proportion = "proportion", mean = "mean", x$type)
  method_label <- if (!is.null(x$method)) paste0(" (", x$method, ")") else ""
  cat(sprintf("Sample size for %s%s\n", type_label, method_label))
  p <- x$params
  resp_rate <- p$resp_rate
  eligibility_rate <- p[["eligibility_rate"]]
  eligible_per_unit <- p[["eligible_per_unit"]]
  field <- .whole_count(x$n)
  if (!is.null(eligible_per_unit)) {
    cat(sprintf(
      "n = %s gross (eligible respondents: %s)",
      .fmt_count_n(field),
      .fmt_count_n(.whole_count(field * (resp_rate %||% 1) * .params_yield(p)))
    ))
  } else if (!is.null(eligibility_rate) ||
             (!is.null(resp_rate) && resp_rate < 1)) {
    net_n <- .whole_count(field * (resp_rate %||% 1) * (eligibility_rate %||% 1))
    cat(sprintf(
      "n = %s gross (net: %s)",
      .fmt_count_n(field),
      .fmt_count_n(net_n)
    ))
  } else {
    cat(sprintf("n = %s", .fmt_count_n(field)))
  }

  parts <- character(0L)
  if (!is.null(p$p)) {
    parts <- c(parts, sprintf("p = %s", .fmt_prob(p$p)))
  }
  if (!is.null(p$var)) {
    parts <- c(parts, sprintf("var = %.2f", p$var))
  }
  # A ratio reports its derived coefficient rather than the three moments it
  # was built from, which would not fit the line. Exact extraction and a type
  # gate, because `p$r` partial-matches `resp_rate` on every other estimand.
  if (identical(x$type, "ratio")) {
    parts <- c(parts, sprintf("r = %.4g", p[["r"]]))
    parts <- c(parts, sprintf("unit_relvar = %.3g", p[["unit_relvar"]]))
  }
  if (!is.null(p$moe)) {
    parts <- c(parts, sprintf("moe = %.3f", p$moe))
  }
  if (!is.null(p$rmoe)) {
    parts <- c(parts, sprintf("rmoe = %.3f", p$rmoe))
  }
  if (!is.null(p$cv)) {
    parts <- c(parts, sprintf("cv = %.3f", p$cv))
  }
  parts <- c(parts, .fmt_deff(p$deff))
  if (!is.null(resp_rate) && resp_rate < 1) {
    parts <- c(parts, sprintf("resp_rate = %s", .fmt_prob(resp_rate)))
  }
  if (!is.null(eligibility_rate)) {
    parts <- c(parts, sprintf("eligibility_rate = %s",
                              .fmt_prob(eligibility_rate)))
  }
  if (!is.null(eligible_per_unit)) {
    parts <- c(parts, sprintf("eligible_per_unit = %.4g", eligible_per_unit))
  }
  if (length(parts) > 0L) {
    cat(sprintf(" (%s)", paste(parts, collapse = ", ")))
  }
  cat("\n")

  cases <- if (!is.null(x$expected_cases)) x$expected_cases * field / x$n
  .print_expected_cases(cases, x$binding, p$min_cases)
  if (x$type %in% c("proportion", "mean", "ratio")) {
    cat("# summary() for the precision at the field size and the derived sizes\n")
  }
}

#' Report the expected count of positive cases, and what it bound against
#'
#' Printed for every proportion result, since it is the number the choice
#' between the interval methods turns on. The binding constraint is named
#' only when there were two of them to choose between.
#' Report the design degrees of freedom a plan can derive
#'
#' Counted from the whole-unit design rather than the continuous optimum,
#' so it is the number the fielded design would have. Silent when the plan
#' does not know its PSU or stratum counts.
#' @keywords internal
#' @noRd
.print_design_df <- function(x) {
  value <- tryCatch(suppressWarnings(design_df(x)), error = function(e) NULL)
  if (is.null(value)) {
    return(invisible(NULL))
  }
  cat(sprintf("design df = %g\n", as.double(value)))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.print_expected_cases <- function(cases, binding = NULL, min_cases = NULL) {
  if (is.null(cases)) {
    return(invisible(NULL))
  }
  cat(sprintf("expected cases = %.1f", cases))
  if (!is.null(binding)) {
    cat(sprintf(" (min_cases = %.4g, binding: %s)", min_cases, binding))
  }
  cat("\n")
  invisible(NULL)
}

#' Format continuous total sample size for display
#' @keywords internal
#' @noRd
.fmt_continuous_n <- function(x) {
  format(signif(x, 8), trim = TRUE, scientific = FALSE)
}

#' Formatter for the sizes a precision result evaluated
#'
#' Whole sizes print as counts, and the counts derived from them round up.
#' A fractional size prints as given, and so does every count derived from
#' it, so the printed design is the one the precision belongs to.
#' @keywords internal
#' @noRd
.evaluated_n_formatter <- function(n) {
  if (all(n == floor(n))) {
    function(v) .fmt_count_n(.whole_count(v))
  } else {
    .fmt_continuous_n
  }
}

#' Header note for a fit read back at its continuous optimum
#' @keywords internal
#' @noRd
.continuous_optimum_note <- function(p) {
  if (isTRUE(p[["continuous_optimum"]])) "continuous optimum"
}

#' Format a count-like numeric value defensively
#' @keywords internal
#' @noRd
.fmt_count_n <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    return("NA")
  }
  if (!is.finite(x)) {
    return(format(x, trim = TRUE))
  }
  format(x, trim = TRUE, scientific = FALSE)
}

#' Field size of a multi-indicator result
#'
#' Separate quotas are fielded one by one, so the total is the sum of the
#' whole quotas, which can exceed the rounded continuous total.
#' @keywords internal
#' @noRd
.multi_field_n <- function(x) {
  dom <- x$domains
  if (!is.null(dom) && !identical(x$params$domain_sampling, "natural")) {
    return(sum(.whole_count(dom$.n)))
  }
  .whole_count(x$n)
}

#' @keywords internal
#' @noRd
.print_multi_n <- function(x) {
  if (!is.null(x$domains)) {
    min_n_label <- if (!is.null(x$params$min_n_domain)) {
      sprintf(", min_n_domain = %g", x$params$min_n_domain)
    } else {
      ""
    }
    natural <- identical(x$params$domain_sampling, "natural")
    cat(sprintf(
      "Multi-indicator sample size (%d domains, %s%s)\n",
      nrow(x$domains),
      if (natural) "natural incidence" else "separate quotas",
      min_n_label
    ))
    # Both numbers are labeled: the total is what you field, the largest
    # single domain is a different quantity and was once printed as 'n'.
    cat(sprintf(
      "n = %s (%s)\n",
      .fmt_count_n(.multi_field_n(x)),
      if (natural) "largest quota / share" else "sum of domain quotas"
    ))
    if (!is.null(x$n_domain_max)) {
      cat(sprintf(
        "Largest single domain = %s (binding: %s)\n",
        .fmt_count_n(.whole_count(x$n_domain_max)),
        .summary_label(x$binding)
      ))
    }
    cat("\n")
    dom <- .truncate_labels(x$domains, c(x$params$domain_cols, ".binding"))
    dom$.n <- .whole_count(dom$.n)
    print(dom, row.names = FALSE, right = FALSE)
  } else {
    cat(sprintf("Multi-indicator sample size\n"))
    cat(sprintf(
      "n = %s (binding: %s)\n",
      .fmt_count_n(.multi_field_n(x)),
      .summary_label(x$binding)
    ))
    cat("\n")
    d <- .truncate_labels(x$detail[, c("name", ".n", ".cv_target"),
                                   drop = FALSE], "name")
    d$.n <- .whole_count(d$.n)
    d$.cv_target <- .fmt_cv_column(d$.cv_target)
    d$.cv_achieved <- .fmt_cv_column(x$detail$.cv_field %||%
                                       x$detail$.cv_achieved)
    d$.binding <- ifelse(x$detail$.binding, "*", "")
    print(d, row.names = FALSE, right = FALSE)
  }
  cat("\n# summary() for each indicator at the field size\n")
}

#' Shorten the label columns of a printed table
#' @keywords internal
#' @noRd
.truncate_labels <- function(df, cols) {
  for (col in intersect(cols, names(df))) {
    df[[col]] <- .summary_label(df[[col]])
  }
  df
}

#' Four-decimal cv column, blank where a row has no precision
#' @keywords internal
#' @noRd
.fmt_cv_column <- function(v) {
  ifelse(is.na(v), "", sprintf("%.4f", v))
}

#' @rdname print.svyplan
#' @export
print.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  if (!is.null(x$indicators)) {
    .print_multi_cluster(x)
  } else {
    .print_single_cluster(x)
  }
  invisible(x)
}

#' The search line of a cluster design, only when its node limit stopped it
#' @keywords internal
#' @noRd
.print_stage_search <- function(search, lead = "", mode = "cv") {
  if (!identical(search$status, "limit")) return(invisible(NULL))
  where <- if (!is.null(search$by_domain) && any(search$by_domain == "limit")) {
    sprintf(" in %s", paste(names(search$by_domain)[search$by_domain == "limit"],
                            collapse = ", "))
  } else {
    ""
  }
  cat(lead, paste(strwrap(sprintf(
    "search: stopped at the node limit of %s%s, the best design found, not proven %s",
    format(search$max_nodes, big.mark = ",", scientific = FALSE), where,
    if (identical(mode, "budget")) "most precise" else "cheapest"
  ), width = 80L, exdent = 2L), collapse = "\n"), "\n", sep = "")
}

#' A cluster heading claims optimality only for a proven whole design
#' @keywords internal
#' @noRd
.cluster_heading <- function(search, multi = FALSE) {
  proven <- is.null(search$status) || identical(search$status, "optimal")
  if (multi) {
    if (proven) "Multi-indicator optimal allocation" else
      "Multi-indicator allocation"
  } else {
    if (proven) "Optimal allocation" else "Allocation"
  }
}

#' Whether a cluster result was sized to a budget or to a precision target
#' @keywords internal
#' @noRd
.cluster_mode <- function(params) {
  params[["mode"]] %||% if (!is.null(params[["budget"]])) "budget" else "cv"
}

#' @keywords internal
#' @noRd
.print_single_cluster <- function(x) {
  cat(trimws(sprintf("%s %d-stage allocation",
                     sub(" ?[Aa]llocation$", "",
                         .cluster_heading(x$operational$search)),
                     x$stages)), "\n", sep = "")
  op <- x$operational
  n_display <- if (!is.null(op)) op$n else .whole_count(x$n)
  total_display <- if (!is.null(op)) op$total_n else prod(n_display)
  stage_labels <- names(x$n)
  if (is.null(stage_labels)) {
    stage_labels <- paste0("stage", seq_along(n_display))
  }
  stage_parts <- vapply(
    seq_along(n_display),
    function(i) {
      sprintf("%s = %s", stage_labels[i], .fmt_count_n(n_display[i]))
    },
    character(1L)
  )
  cat("field design: ")
  cat(paste(stage_parts, collapse = " | "))
  # Every stage's loss removes observations from the same total, so the net
  # figure reads their product.
  resp_rate <- prod(x$params[["resp_rate"]] %||% 1) *
    (x$params[["eligibility_rate"]] %||% 1)
  eligible_per_unit <- x$params[["eligible_per_unit"]]
  if (!is.null(eligible_per_unit)) {
    cat(sprintf(
      " -> total n = %s (eligible respondents: %s)\n",
      .fmt_count_n(total_display),
      .fmt_count_n(.whole_count(total_display * resp_rate * eligible_per_unit))
    ))
  } else if (resp_rate < 1) {
    net_total <- .whole_count(total_display * resp_rate)
    cat(sprintf(
      " -> total n = %s (net: %s)\n",
      .fmt_count_n(total_display),
      .fmt_count_n(net_total)
    ))
  } else {
    cat(sprintf(" -> total n = %s\n", .fmt_count_n(total_display)))
  }
  .print_stage_rates(x$params)
  .print_stage_fpc(x$params)
  fc <- x$params$fixed_cost
  op_cv <- if (!is.null(op)) op$cv else x$cv
  op_cost <- if (!is.null(op)) op$cost else x$cost
  if (!is.null(fc) && fc > 0) {
    cat(sprintf("cv = %.4f, cost = %.0f (fixed: %.0f)\n", op_cv, op_cost, fc))
  } else {
    cat(sprintf("cv = %.4f, cost = %.0f\n", op_cv, op_cost))
  }
  .print_stage_search(op$search, mode = .cluster_mode(x$params))
  .print_design_df(x)
  cat("# summary() for the continuous optimum, the design effect and costs\n")
}

#' @keywords internal
#' @noRd
.print_multi_cluster <- function(x) {
  if (!is.null(x$domains)) {
    joint_label <- if (isTRUE(x$params$joint)) ", joint" else ""
    min_n_label <- if (!is.null(x$params$min_n_domain)) {
      sprintf(", min_n_domain = %g", x$params$min_n_domain)
    } else {
      ""
    }
    cat(sprintf(
      "%s (%d-stage, %d domains%s%s)\n",
      .cluster_heading(x$operational$search, multi = TRUE), x$stages,
      nrow(x$domains),
      joint_label,
      min_n_label
    ))
    op <- x$operational
    if (!is.null(op$domains)) {
      dom <- op$domains
      cost_display <- op$cost
    } else {
      dom <- x$domains
      stage_cols <- names(x$n)
      for (col in stage_cols) {
        dom[[col]] <- .whole_count(dom[[col]])
      }
      dom$.total_n <- apply(
        dom[, stage_cols, drop = FALSE],
        1,
        prod
      )
      cost_display <- x$cost
    }
    fc <- x$params$fixed_cost
    if (!is.null(fc) && fc > 0) {
      cat(sprintf(
        "Total n = %s, cost = %.0f (fixed: %.0f)\n",
        .fmt_count_n(sum(dom$.total_n)),
        cost_display,
        fc
      ))
    } else {
      cat(sprintf(
        "Total n = %s, cost = %.0f\n",
        .fmt_count_n(sum(dom$.total_n)),
        cost_display
      ))
    }
    if (!isTRUE(x$params$joint) && !is.null(x$params$budget)) {
      cat(sprintf(
        "budget = %.0f for each domain, total cost = %.0f\n",
        x$params$budget, cost_display
      ))
    }
    .print_stage_fpc(x$params)
    dom$.cv <- sprintf("%.4f", dom$.cv)
    dom$.cost <- sprintf("%.0f", dom$.cost)
    dom <- .truncate_labels(dom, c(x$params$domain_cols, ".binding"))
    cat("\n")
    print(dom, row.names = FALSE, right = FALSE)
    .print_stage_search(op$search, lead = "\n", mode = .cluster_mode(x$params))
    cat("\n# summary() for each indicator's precision at the field design\n")
  } else {
    cat(sprintf("%s (%d-stage)\n",
                .cluster_heading(x$operational$search, multi = TRUE),
                x$stages))
    op <- x$operational
    n_display <- if (!is.null(op)) op$n else .whole_count(x$n)
    total_display <- if (!is.null(op)) op$total_n else prod(n_display)
    stage_labels <- names(x$n)
    if (is.null(stage_labels)) {
      stage_labels <- paste0("stage", seq_along(n_display))
    }
    stage_parts <- vapply(
      seq_along(n_display),
      function(i) {
        sprintf("%s = %s", stage_labels[i], .fmt_count_n(n_display[i]))
      },
      character(1L)
    )
    cat("field design: ")
    cat(paste(stage_parts, collapse = " | "))
    cat(sprintf(" -> total n = %s\n", .fmt_count_n(total_display)))
    .print_stage_fpc(x$params)
    d <- x$detail
    cv_field <- op$cv_by_target %||% d$.cv_achieved
    op_cost <- if (!is.null(op)) op$cost else x$cost
    fc <- x$params$fixed_cost
    cost_label <- if (!is.null(fc) && fc > 0) {
      sprintf("%.0f (fixed: %.0f)", op_cost, fc)
    } else {
      sprintf("%.0f", op_cost)
    }
    worst <- which.max(cv_field)
    bind <- if (!is.null(op$binding)) {
      match(op$binding, d$name)
    } else {
      which.max(cv_field / d$.cv_target)
    }
    cat(sprintf(
      "worst cv = %.4f (%s), cost = %s\n",
      cv_field[worst], .summary_label(d$name[worst]), cost_label
    ))
    cat(sprintf("binding: %s\n", .summary_label(d$name[bind])))
    short <- cv_field > d$.cv_target * (1 + 1e-9)
    if (any(short)) {
      cat(sprintf(
        "short of target: %s\n",
        paste(.summary_label(d$name[short]), collapse = ", ")
      ))
    }
    cat("\n")
    d <- .truncate_labels(d[, c("name", ".n", ".cv_target"), drop = FALSE],
                          "name")
    d$.n <- .whole_count(d$.n)
    d$.cv_target <- sprintf("%.4f", d$.cv_target)
    d$.cv_achieved <- sprintf("%.4f", cv_field)
    d$.binding <- ifelse(seq_len(nrow(d)) == bind, "*", "")
    print(d, row.names = FALSE, right = FALSE)
    .print_stage_search(x$operational$search, lead = "\n",
                        mode = .cluster_mode(x$params))
    cat("\n# summary() for the continuous optimum and each indicator's precision\n")
  }
}

#' @rdname print.svyplan
#' @export
print.svyplan_prec <- function(x, ...) {
  .check_unused_dots(...)
  if (identical(x$type, "alloc_joint") ||
      (x$type == "alloc" && identical(x$method, "bethel"))) {
    .print_bethel_prec(x)
  } else if (x$type == "multi" ||
             (identical(x$type, "cluster") && !is.null(x$detail))) {
    .print_multi_prec(x)
  } else if (x$type == "change") {
    .print_change_prec(x)
  } else if (x$type == "pooled") {
    .print_pooled_prec(x)
  } else {
    .print_single_prec(x)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.print_pooled_prec <- function(x) {
  p <- x$params
  cat(sprintf(
    "Sampling precision for pooled estimate (%s)\n",
    paste(c(.pooled_scale_label(p), .continuous_optimum_note(p)),
          collapse = ", ")
  ))
  fmt_n <- .evaluated_n_formatter(p$n)
  cat(sprintf("n = %s per occasion", fmt_n(p$n)))
  if (.net_rate(p) < 1) {
    cat(sprintf(" (net: %s)", fmt_n(p$n * .net_rate(p))))
  }
  cat(sprintf(", %d occasions", p$occasions))
  parts <- c(.fmt_pooled_estimand(p), .fmt_deff(p$deff))
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_lag_overlap(p))
  cat(sprintf("se = %.4g, moe = %.4g", x$se, x$moe))
  if (!is.na(x$cv)) {
    cat(sprintf(", cv = %.4g", x$cv))
  }
  if (!is.null(x$rmoe) && !is.na(x$rmoe)) {
    cat(sprintf(", rmoe = %.4g", x$rmoe))
  }
  cat("\n")
}

#' @keywords internal
#' @noRd
.print_change_prec <- function(x) {
  p <- x$params
  cat(sprintf(
    "Sampling precision for change (%s)\n",
    paste(c(.change_scale_label(p), .continuous_optimum_note(p)),
          collapse = ", ")
  ))
  cat(.fmt_change_n(p$n, .net_rate(p), exact = TRUE))
  parts <- c(.fmt_change_estimand(p), .fmt_deff(p$deff))
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_overlap(p))
  cat(sprintf("se = %.4g, moe = %.4g", x$se, x$moe))
  if (!is.na(x$cv)) {
    cat(sprintf(", cv = %.4g", x$cv))
  }
  if (!is.null(x$rmoe) && !is.na(x$rmoe)) {
    cat(sprintf(", rmoe = %.4g", x$rmoe))
  }
  cat("\n")
}

#' @keywords internal
#' @noRd
.print_single_prec <- function(x) {
  type_label <- switch(
    x$type,
    proportion = "proportion",
    mean = "mean",
    cluster = "cluster",
    twophase = "two-phase",
    alloc = "stratum allocation",
    x$type
  )
  p <- x$params

  if (x$type == "cluster") {
    stages <- p$stages
    note <- .continuous_optimum_note(p)
    cat(sprintf(
      "Sampling precision for %d-stage cluster%s\n",
      stages,
      if (is.null(note)) "" else sprintf(" (%s)", note)
    ))
    # The cv is that of p$n. Rounding a fractional stage up would print
    # another design beside it.
    fmt_n <- .evaluated_n_formatter(p$n)
    total_display <- prod(p$n)
    stage_labels <- names(p$n)
    if (is.null(stage_labels)) {
      stage_labels <- paste0("stage", seq_along(p$n))
    }
    stage_parts <- vapply(
      seq_along(p$n),
      function(i) {
        sprintf("%s = %s", stage_labels[i], fmt_n(p$n[i]))
      },
      character(1L)
    )
    cat(paste(stage_parts, collapse = " | "))
    cat(sprintf(" -> total n = %s", fmt_n(total_display)))
    resp_all <- prod(p[["resp_rate"]] %||% 1) * (p[["eligibility_rate"]] %||% 1)
    if (!is.null(p[["eligible_per_unit"]])) {
      cat(sprintf(
        " (eligible respondents: %s, eligible_per_unit = %.4g)",
        fmt_n(total_display * resp_all * p[["eligible_per_unit"]]),
        p[["eligible_per_unit"]]
      ))
    } else if (resp_all < 1) {
      cat(sprintf(" (net: %s)", fmt_n(total_display * resp_all)))
    }
    cat("\n")
    .print_stage_fpc(p)
  } else {
    notes <- c(
      x$method,
      if (!is.null(x$solved)) paste("solved for", x$solved),
      .continuous_optimum_note(p)
    )
    note_label <- if (length(notes) > 0L) {
      sprintf(" (%s)", paste(notes, collapse = ", "))
    } else {
      ""
    }
    cat(sprintf("Sampling precision for %s%s\n", type_label, note_label))
    if (identical(x$type, "twophase")) {
      .print_twophase_prec_sizes(x)
    } else {
      .print_single_prec_sizes(x)
    }
  }

  if (!is.null(x$solved)) {
    cat(sprintf("%s = %.4g\n", x$solved, p[[x$solved]]))
  }

  if (!is.na(x$se[1L])) {
    cat(sprintf("se = %.4f, moe = %.4f", x$se[1L], x$moe[1L]))
  }
  if (!is.na(x$cv[1L])) {
    if (!is.na(x$se[1L])) {
      cat(", ")
    }
    cat(sprintf("cv = %.4f", x$cv[1L]))
  }
  if (!is.null(x$rmoe) && !is.na(x$rmoe[1L])) {
    cat(sprintf(", rmoe = %.4f", x$rmoe[1L]))
  }
  cat("\n")

  .print_expected_cases(x$expected_cases)
  .print_domains_block(x$domains)
}

#' Size line of a one-stage or allocation precision result
#'
#' The precision is that of `p$n`. Rounding a fractional size up would print
#' another design beside it.
#' @keywords internal
#' @noRd
.print_single_prec_sizes <- function(x) {
  p <- x$params
  fmt_n <- .evaluated_n_formatter(p$n)
  cat(sprintf("n = %s", fmt_n(sum(p$n))))
  if (length(p$n) > 1L) {
    cat(sprintf(" (%d strata)", length(p$n)))
  }
  if (identical(x$type, "alloc")) {
    resp_rate <- .alloc_summary_response_rates(x)$combined
    if (any(resp_rate != 1)) {
      cat(sprintf(", expected respondents = %.1f", sum(p$n * resp_rate)))
    }
    cat("\n")
    return(invisible(NULL))
  }
  resp_rate <- if (is.null(p[["eligibility_rate"]])) {
    p$resp_rate
  } else {
    (p$resp_rate %||% 1) * p[["eligibility_rate"]]
  }
  eligible_per_unit <- p[["eligible_per_unit"]]
  if (!is.null(eligible_per_unit)) {
    cat(sprintf(
      " (eligible respondents: %s, eligible_per_unit = %.4g)",
      fmt_n(sum(p$n * (resp_rate %||% 1) * eligible_per_unit)),
      eligible_per_unit
    ))
  } else if (!is.null(resp_rate) && all(resp_rate < 1)) {
    cat(sprintf(" (net: %s)", fmt_n(sum(p$n * resp_rate))))
  }
  cat("\n")
  invisible(NULL)
}

#' Size line of a two-phase precision result
#'
#' Both phases as evaluated, the phase-2 sample summed over its strata as
#' the `n_twophase()` field line sums it.
#' @keywords internal
#' @noRd
.print_twophase_prec_sizes <- function(x) {
  p <- x$params
  d <- x$detail
  n <- c(p[["n_phase1"]], sum(d$n_issued))
  fmt_n <- .evaluated_n_formatter(n)
  cat(sprintf("n_phase1 = %s | n_phase2 = %s", fmt_n(n[1L]), fmt_n(n[2L])))
  net <- c(n[1L] * (p[["resp_rate"]] %||% 1), sum(d$n_resp))
  if (any(net < n)) {
    cat(sprintf(" (net: %s | %s)", fmt_n(net[1L]), fmt_n(net[2L])))
  }
  cat("\n")
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.print_multi_prec <- function(x) {
  cat("Multi-indicator sampling precision\n")
  .print_stage_fpc(x$params)
  if (!is.null(x$detail)) {
    tab <- x$detail
    if (!is.null(x$params$domain_cols)) {
      keep <- c(x$params$domain_cols, "name", ".se", ".moe", ".rmoe", ".cv")
      tab <- tab[intersect(keep, names(tab))]
    }
    cat("\n")
    print(tab, row.names = FALSE, right = FALSE)
  }
}

#' @rdname print.svyplan
#' @export
print.svyplan_varcomp <- function(x, ...) {
  .check_unused_dots(...)
  if (!is.null(x$strata)) {
    cat(sprintf(
      "Variance components (%d-stage, %d strata%s)\n",
      x$stages,
      nrow(x$strata),
      .varcomp_design_label(x)
    ))
    tab <- x$strata
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(v) sprintf("%.4f", v))
    cat("\n")
    print(tab, row.names = FALSE, right = FALSE)
    return(invisible(x))
  }
  if (identical(x$source, "deff")) {
    return(.print_varcomp_deff(x))
  }
  cat(sprintf(
    "Variance components (%d-stage%s)\n",
    x$stages,
    .varcomp_design_label(x)
  ))
  cat(sprintf("varb = %.4f", x$varb))
  varw_names <- names(x$varw)
  for (i in seq_along(x$varw)) {
    label <- if (!is.null(varw_names) && nzchar(varw_names[i])) {
      varw_names[i]
    } else {
      "varw"
    }
    cat(sprintf(", %s = %.4f", label, x$varw[i]))
  }
  cat("\n")
  cat(sprintf("icc = %s\n", paste(sprintf("%.4f", x$icc), collapse = ", ")))
  cat(sprintf(
    "var_ratio = %s\n",
    paste(sprintf("%.4f", x$var_ratio), collapse = ", ")
  ))
  cat(sprintf("Unit relvariance = %.4f\n", x$unit_relvar))

  invisible(x)
}

#' How the components were obtained, for the print header
#'
#' The stage-wise selection matters downstream: components from a
#' size-proportional stage do not describe a design that selects that stage
#' by simple random sampling, and the reverse.
#' @keywords internal
#' @noRd
.varcomp_design_label <- function(x) {
  sel <- x$selection
  if (is.null(sel) || anyNA(sel)) {
    return("")
  }
  stage <- sub("^stage_([0-9]+)$", "stage \\1",
               names(sel) %||% paste0("stage_", seq_along(sel)))
  paste0(
    ", ",
    paste(toupper(sel), stage, collapse = ", "),
    if (identical(x$input, "weighted")) ", weighted sample" else ", frame",
    if (identical(x$estimator, "ratio")) ", ratio estimator" else ""
  )
}

#' Report an icc backed out of a published design effect
#' @keywords internal
#' @noRd
.print_varcomp_deff <- function(x) {
  take <- x$params$take %||% x$params$n_per_psu
  nominal <- x$params$take_nominal %||% x$params$n_per_psu_nominal
  cat("Variance components (2-stage, from a design effect)\n")
  cat(sprintf("icc = %.4f\n", x$icc))
  cat(sprintf("var_ratio = %.4f\n", x$var_ratio))
  cat(sprintf(
    "deff = %.4f at take = %.4g%s\n",
    x$var_ratio * (1 + x$icc * (take - 1)),
    take,
    if (isTRUE(abs(take - nominal) > 1e-8)) {
      sprintf(" (nominal %.4g, size-weighted)", nominal)
    } else {
      ""
    }
  ))
  cat("varb, varw and unit_relvar are not identified by a design effect\n")
  invisible(x)
}

#' @rdname print.svyplan
#' @export
print.svyplan_power <- function(x, ...) {
  .check_unused_dots(...)
  type_label <- switch(
    x$type,
    proportion = "proportions",
    mean = "means",
    did_prop = "DiD proportions",
    did_mean = "DiD means",
    x$type
  )
  solved_label <- switch(
    x$solved,
    n = "sample size",
    power = "power",
    mde = "minimum detectable effect"
  )
  cat(sprintf(
    "Power analysis for %s (solved for %s)\n",
    type_label,
    solved_label
  ))

  p <- x$params
  # Respondents in scope per issued unit. With eligible_per_unit they are
  # eligible persons, which is what the power is computed on.
  rate_all <- (p$resp_rate %||% 1) * .params_yield(p)
  show_net <- !is.null(p[["eligible_per_unit"]]) || rate_all < 1
  net_label <- if (is.null(p[["eligible_per_unit"]])) "net" else
    "eligible respondents"

  is_did <- startsWith(x$type, "did")
  if (length(x$n) == 2L) {
    n1 <- .whole_count(x$n[1])
    n2 <- .whole_count(x$n[2])
    lab1 <- if (is_did) "n_treat" else "n1"
    lab2 <- if (is_did) "n_control" else "n2"
    if (show_net) {
      cat(sprintf(
        "%s = %s, %s = %s (total = %s, %s: %s), power = %.3f, effect = %.4f\n",
        lab1,
        .fmt_count_n(n1),
        lab2,
        .fmt_count_n(n2),
        .fmt_count_n(n1 + n2),
        net_label,
        .fmt_count_n(.whole_count((n1 + n2) * rate_all)),
        x$power,
        x$effect
      ))
    } else {
      cat(sprintf(
        "%s = %s, %s = %s (total = %s), power = %.3f, effect = %.4f\n",
        lab1,
        .fmt_count_n(n1),
        lab2,
        .fmt_count_n(n2),
        .fmt_count_n(n1 + n2),
        x$power,
        x$effect
      ))
    }
  } else if (show_net) {
    net_n <- .whole_count(x$n * rate_all)
    cat(sprintf(
      "n = %s (%s: %s, per group), power = %.3f, effect = %.4f\n",
      .fmt_count_n(.whole_count(x$n)),
      net_label,
      .fmt_count_n(net_n),
      x$power,
      x$effect
    ))
  } else {
    cat(sprintf(
      "n = %s (per group), power = %.3f, effect = %.4f\n",
      .fmt_count_n(.whole_count(x$n)),
      x$power,
      x$effect
    ))
  }

  parts <- character(0L)
  if (!is.null(p$treat)) {
    parts <- c(
      parts,
      sprintf(
        "treat = (%.3f, %.3f)",
        p$treat[1],
        p$treat[2]
      )
    )
  }
  if (!is.null(p$control)) {
    parts <- c(
      parts,
      sprintf(
        "control = (%.3f, %.3f)",
        p$control[1],
        p$control[2]
      )
    )
  }
  if (!is.null(p$p1)) {
    parts <- c(parts, sprintf("p1 = %.3f", p$p1))
  }
  if (!is.null(p$p2)) {
    parts <- c(parts, sprintf("p2 = %.3f", p$p2))
  }
  parts <- c(parts, sprintf("alpha = %s", .fmt_prob(p$alpha)))
  parts <- c(parts, .fmt_deff(p$deff))
  if (!is.null(p$resp_rate) && p$resp_rate < 1) {
    parts <- c(parts, sprintf("resp_rate = %s", .fmt_prob(p$resp_rate)))
  }
  if (!is.null(p[["eligibility_rate"]])) {
    parts <- c(parts, sprintf("eligibility_rate = %s",
                              .fmt_prob(p[["eligibility_rate"]])))
  }
  if (!is.null(p[["eligible_per_unit"]])) {
    parts <- c(parts, sprintf("eligible_per_unit = %.4g",
                              p[["eligible_per_unit"]]))
  }
  if (!is.null(p$var) && length(p$var) == 4L) {
    parts <- c(
      parts,
      sprintf(
        "var = (%.2f, %.2f, %.2f, %.2f)",
        p$var[1],
        p$var[2],
        p$var[3],
        p$var[4]
      )
    )
  }
  if (!is.null(p$overlap) && p$overlap > 0) {
    parts <- c(parts, sprintf("overlap = %s", .fmt_prob(p$overlap)))
    parts <- c(parts, sprintf("overlap_cor = %s", .fmt_prob(p$overlap_cor)))
    if (nzchar(.fmt_overlap_draw(p))) parts <- c(parts, "independent draws")
  }
  if (!is.null(p$alternative) && p$alternative == "one.sided") {
    parts <- c(parts, "one-sided")
  }
  if (!is.null(p$method) && p$method != "wald") {
    parts <- c(parts, sprintf("method = %s", p$method))
  }
  if (!is.null(p$ratio) && p$ratio != 1) {
    parts <- c(parts, sprintf("ratio = %.2g", p$ratio))
  }
  cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))

  invisible(x)
}

#' @rdname print.svyplan
#' @export
format.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  n_field <- if (identical(x$type, "multi")) .multi_field_n(x) else .whole_count(x$n)
  paste0("svyplan_n [n=", paste(n_field, collapse = " then "), ", ",
         x$type, "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_alloc <- function(x, ...) {
  .check_unused_dots(...)
  n_field <- if (!is.null(x$operational)) x$operational$n else .whole_count(x$n)
  cls <- if (inherits(x, "svyplan_alloc_joint")) "svyplan_alloc_joint" else "svyplan_alloc"
  paste0(cls, " [n=", n_field, ", ", x$method, "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  total <- if (!is.null(x$operational)) {
    x$operational$total_n
  } else {
    prod(.whole_count(x$n))
  }
  paste0(
    "svyplan_cluster [",
    x$stages,
    "-stage, total_n=",
    .fmt_count_n(total),
    " (unrounded: ",
    .fmt_continuous_n(x$total_n),
    ")]"
  )
}

#' @rdname print.svyplan
#' @export
format.svyplan_prec <- function(x, ...) {
  .check_unused_dots(...)
  value <- if (is.na(x$se[1L]) && !is.na(x$cv[1L])) {
    sprintf("cv=%.4f", x$cv[1L])
  } else {
    sprintf("se=%.4f", x$se[1L])
  }
  paste0("svyplan_prec [", x$type, ", ", value, "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_twophase <- function(x, ...) {
  .check_unused_dots(...)
  o <- .twophase_shown(x)
  paste0(
    "svyplan_twophase [n_phase1=", .fmt_count_n(o$n[["n_phase1"]]),
    ", n_phase2=", .fmt_count_n(o$n[["n_phase2"]]),
    if (!is.na(o$cv)) paste0(", cv=", sprintf("%.4f", o$cv)), "]"
  )
}

#' @rdname print.svyplan
#' @export
format.svyplan_varcomp <- function(x, ...) {
  .check_unused_dots(...)
  paste0(
    "svyplan_varcomp [",
    x$stages,
    "-stage",
    if (identical(x$source, "deff")) ", from deff" else "",
    "]"
  )
}

#' @rdname print.svyplan
#' @export
format.svyplan_power <- function(x, ...) {
  .check_unused_dots(...)
  n_str <- if (length(x$n) == 2L) {
    paste0(.whole_count(x$n[1]), ",", .whole_count(x$n[2]))
  } else {
    as.character(.whole_count(x$n))
  }
  paste0("svyplan_power [", x$type, ", ", x$solved, ", n=", n_str, "]")
}

#' Compute method-specific confidence interval limits for a proportion
#' @keywords internal
#' @noRd
.prop_ci_limits <- function(
  p,
  n,
  alpha,
  N = Inf,
  deff = 1,
  resp_rate = 1,
  method = "wald",
  df = NULL
) {
  # One variance through n_eff, as .prec_engine_prop() does, so the interval
  # and the reported margin of error agree.
  n_net <- n * resp_rate
  n_eff <- .effective_from_n(n_net, N, deff)
  z <- .q_alpha(alpha, df)
  method <- match.arg(method, c("wald", "wilson", "logodds", "beta"))

  if (is.infinite(n_eff)) {
    return(c(p, p))
  }

  if (method == "wald") {
    se <- sqrt(p * (1 - p) / n_eff)
    lo <- p - z * se
    hi <- p + z * se
  } else if (method == "wilson") {
    den <- 1 + z^2 / n_eff
    center <- (p + z^2 / (2 * n_eff)) / den
    half <- z * sqrt(p * (1 - p) / n_eff + z^2 / (4 * n_eff^2)) / den
    lo <- center - half
    hi <- center + half
  } else if (method == "logodds") {
    eta <- qlogis(p)
    se_eta <- sqrt(1 / (n_eff * p * (1 - p)))
    lo <- plogis(eta - z * se_eta)
    hi <- plogis(eta + z * se_eta)
  } else {
    limits <- .beta_limits(p, .kg_effective(n_eff, n_net, alpha, df), alpha)
    lo <- limits[1L]
    hi <- limits[2L]
  }

  c(max(lo, 0), min(hi, 1))
}

#' Build CI matrix with standard labels
#' @keywords internal
#' @noRd
.ci_matrix <- function(lo, hi, alpha) {
  matrix(
    c(lo, hi),
    nrow = 1L,
    dimnames = list(
      "",
      c(
        sprintf("%.1f %%", 100 * alpha / 2),
        sprintf("%.1f %%", 100 * (1 - alpha / 2))
      )
    )
  )
}

#' Confidence intervals for svyplan results
#'
#' Compute a confidence interval for the parameter a sizing or precision
#' result was built around, at the planned sample size.
#'
#' @param object A sizing or precision result from [n_prop()]/[prec_prop()],
#'   [n_mean()]/[prec_mean()], [n_ratio()]/[prec_ratio()],
#'   [n_change()]/[prec_change()], or [n_pooled()]/[prec_pooled()].
#' @param parm Ignored (included for S3 consistency with [confint()]).
#' @param level Confidence level (default 0.95). This is independent of the
#'   `alpha` used to size the design, so a plan built at `alpha = 0.05` can
#'   be reported at any level.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @return A one-row, two-column matrix with the lower and upper confidence
#'   limits, named for the percentiles they correspond to.
#'
#' @details
#' For proportions, the interval type matches the `method` the result was
#' computed with (`"wald"`, `"wilson"`, `"logodds"`, or `"beta"`). All four
#' methods read a stored `df`. Beta applies it through effective-size scaling,
#' while the other three substitute the corresponding t quantile. Only the Wald interval is
#' symmetric about `p`, so for the other three the limits are not
#' `p` plus or minus the reported `moe`. `$moe` remains half the interval
#' width, and `confint()` is the way to read where the interval actually
#' sits. All four apply `deff`, `resp_rate`, and the finite population
#' correction through the same effective size the sizing functions use.
#'
#' For means, ratios, changes and pooled estimates the interval is symmetric
#' about the estimand, at the normal quantile or, when the result carries a
#' `df`, the t quantile it was built with. A mean requires `mu` in the
#' original call, since a size targeted on `moe` alone knows the spread but
#' not the level. A ratio always has `r`, so it needs nothing extra.
#'
#' A ratio interval is the first-order symmetric one, `r` plus or minus
#' `q * se`. It is not a Fieller interval and does not correct the ratio
#' estimator's bias.
#'
#' Multi-indicator results (`n_multi()`, `prec_multi()`) and allocation
#' results have no single parameter to bound, and error rather than
#' returning an interval for an arbitrary component.
#'
#' @seealso [n_prop()] and [prec_prop()] for the methods themselves,
#'   [print.svyplan] for printing and coercion.
#'
#' @examples
#' # confint on a proportion sample size
#' res <- n_prop(p = 0.3, moe = 0.05)
#' confint(res)
#'
#' # confint at 90% level
#' confint(res, level = 0.90)
#'
#' # confint on a mean (requires mu)
#' res_mean <- n_mean(var = 100, mu = 50, moe = 2)
#' confint(res_mean)
#'
#' # confint on a precision result
#' prec <- prec_prop(p = 0.3, n = 400)
#' confint(prec)
#'
#' # The Korn-Graubard interval is asymmetric for a rare outcome
#' confint(prec_prop(p = 0.02, n = 150, method = "beta"))
#'
#' # confint on a ratio, which always carries its center r
#' confint(n_ratio(r = 2, cv_num = 0.5, cv_den = 0.3,
#'                 component_cor = 0.4, cv = 0.05))
#'
#' # confint on a repeated-survey result, the pooled level here
#' confint(prec_pooled(var = 100, n = 300, occasions = 2, mu = 50))
#'
#' @name confint.svyplan
NULL

#' @rdname confint.svyplan
#' @export
confint.svyplan_alloc <- function(object, parm, level = 0.95, ...) {
  .refuse_alloc_object(object, "confint")
}

#' Refuse an allocation where a single-estimate result is expected
#' @keywords internal
#' @noRd
.refuse_alloc_object <- function(x, fn, type = NULL) {
  if (!inherits(x, "svyplan_alloc")) return(invisible(NULL))
  stop(
    sprintf(
      "%s. Read an allocation back with prec_alloc()",
      if (is.null(type)) {
        sprintf("%s not supported for this type", fn)
      } else {
        sprintf("%s requires a svyplan_n of type '%s'", fn, type)
      }
    ),
    call. = FALSE
  )
}

#' @rdname confint.svyplan
#' @export
confint.svyplan_n <- function(object, parm, level = 0.95, ...) {
  .check_unused_dots(...)
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop("'level' must be a number in (0, 1)", call. = FALSE)
  }
  if (object$type == "multi" || is.na(object$se)) {
    stop(
      "confint requires a single-indicator svyplan_n object with computed SE",
      call. = FALSE
    )
  }

  p <- object$params
  if (object$type == "proportion") {
    alpha <- 1 - level
    ci <- .prop_ci_limits(
      p = p$p,
      n = object$n,
      alpha = alpha,
      N = (p$N %||% Inf) * .params_yield(p),
      deff = p$deff %||% 1,
      resp_rate = (p$resp_rate %||% 1) * .params_yield(p),
      method = object$method %||% "wald",
      df = p$df
    )
    return(.ci_matrix(ci[1L], ci[2L], alpha))
  } else if (object$type == "mean") {
    est <- p$mu
    if (is.null(est)) {
      stop(
        "'mu' is required to compute a confidence interval for a mean",
        call. = FALSE
      )
    }
  } else if (object$type == "ratio") {
    # Always present: a ratio cannot be built without it, unlike a mean's 'mu'.
    est <- p$r
  } else if (object$type == "change") {
    est <- .change_ci_estimand(p)
  } else if (object$type == "pooled") {
    est <- .pooled_ci_estimand(p)
  } else {
    stop("confint not supported for this type", call. = FALSE)
  }

  alpha <- 1 - level
  z <- .q_alpha(alpha, p$df)
  moe <- z * object$se

  lo <- est - moe
  hi <- est + moe
  .ci_matrix(lo, hi, alpha)
}

#' The change an interval on a change is centered on
#'
#' A change sized from `moe` alone has no known level, and the interval
#' would then be centered on nothing. The message names `change` rather than
#' the interval so the fix is the argument to add.
#' @keywords internal
#' @noRd
.change_ci_estimand <- function(p) {
  if (is.null(p$change)) {
    stop(
      "'change' (or 'p') is required to compute a confidence interval for a change",
      call. = FALSE
    )
  }
  p$change
}

#' The level an interval on a pooled estimate is centered on
#'
#' A size solved from `moe` alone knows the spread but not the level, and
#' the interval would then be centered on nothing. The message names `mu`
#' rather than the interval so the fix is the argument to add.
#' @keywords internal
#' @noRd
.pooled_ci_estimand <- function(p) {
  if (is.null(p$mu)) {
    stop(
      "'mu' (or 'p') is required to compute a confidence interval for a pooled estimate",
      call. = FALSE
    )
  }
  p$mu
}

#' @rdname confint.svyplan
#' @export
confint.svyplan_prec <- function(object, parm, level = 0.95, ...) {
  .check_unused_dots(...)
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop("'level' must be a number in (0, 1)", call. = FALSE)
  }
  p <- object$params
  if (object$type == "proportion") {
    alpha <- 1 - level
    ci <- .prop_ci_limits(
      p = p$p,
      n = p$n,
      alpha = alpha,
      N = (p$N %||% Inf) * .params_yield(p),
      deff = p$deff %||% 1,
      resp_rate = (p$resp_rate %||% 1) * .params_yield(p),
      method = object$method %||% "wald",
      df = p$df
    )
    return(.ci_matrix(ci[1L], ci[2L], alpha))
  } else if (object$type == "mean") {
    est <- p$mu
    if (is.null(est)) {
      stop(
        "'mu' is required to compute a confidence interval for a mean",
        call. = FALSE
      )
    }
  } else if (object$type == "ratio") {
    # Always present: a ratio cannot be built without it, unlike a mean's 'mu'.
    est <- p$r
  } else if (object$type == "change") {
    est <- .change_ci_estimand(p)
  } else if (object$type == "pooled") {
    est <- .pooled_ci_estimand(p)
  } else {
    stop("confint not supported for this precision type", call. = FALSE)
  }

  alpha <- 1 - level
  z <- .q_alpha(alpha, p$df)
  moe <- z * object$se

  lo <- est - moe
  hi <- est + moe
  .ci_matrix(lo, hi, alpha)
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  if (!is.null(x$operational)) {
    return(.coerce_count(x$operational$n))
  }
  if (identical(x$type, "multi")) {
    return(.coerce_count(.multi_field_n(x)))
  }
  .coerce_count(.whole_count(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  x$n
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_alloc <- function(x, ...) {
  .check_unused_dots(...)
  if (!is.null(x$operational)) {
    return(.coerce_count(x$operational$n))
  }
  .coerce_count(.whole_count(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_alloc <- function(x, ...) {
  .check_unused_dots(...)
  x$n
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_alloc <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  as.data.frame(x$detail, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  n <- if (!is.null(x$operational)) x$operational$n else .whole_count(x$n)
  setNames(.coerce_count(n), names(n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  x$n
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_n <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  out <- if (identical(x$type, "multi")) {
    x$domains %||% x$detail
  } else {
    out <- data.frame(
      n = x$n,
      n_int = .count_column(.whole_count(x$n)),
      se = x$se,
      moe = x$moe,
      cv = x$cv,
      rmoe = x[["rmoe"]] %||% NA_real_,
      stringsAsFactors = stringsAsFactors
    )
    if (nrow(out) == 2L) cbind(occasion = 1:2, out) else out
  }
  as.data.frame(out, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_twophase <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  as.data.frame(x$detail, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_prec <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  if (!is.null(x$detail)) {
    return(as.data.frame(x$detail, row.names = row.names, optional = optional))
  }

  if (identical(x$type, "cluster")) {
    stage_n <- x$params$n
    stage_names <- names(stage_n)
    if (is.null(stage_names) || any(!nzchar(stage_names))) {
      stage_names <- paste0("stage", seq_along(stage_n))
    }
    out <- as.data.frame(
      as.list(stats::setNames(as.numeric(stage_n), stage_names)),
      optional = optional,
      stringsAsFactors = stringsAsFactors
    )
    out$total_n <- prod(stage_n)
    out$se <- x$se
    out$moe <- x$moe
    out$cv <- x$cv
  } else {
    out <- data.frame(
      n = x$params$n,
      se = x$se,
      moe = x$moe,
      cv = x$cv,
      rmoe = x[["rmoe"]] %||% NA_real_,
      stringsAsFactors = stringsAsFactors
    )
    if (nrow(out) == 2L) out <- cbind(occasion = 1:2, out)
  }

  as.data.frame(out, row.names = row.names, optional = optional)
}

#' A whole count as an integer column, or as a double past the integer range
#' @keywords internal
#' @noRd
.count_column <- function(x) {
  if (all(is.na(x) | abs(x) <= .Machine$integer.max)) as.integer(x) else x
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_varcomp <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  out <- if (!is.null(x$strata)) {
    x$strata
  } else if (identical(x$source, "deff")) {
    data.frame(
      stages = x$stages,
      varb = x$varb,
      varw = x$varw,
      icc = x$icc,
      var_ratio = x$var_ratio,
      unit_relvar = x$unit_relvar,
      take = x$params$take %||% x$params$n_per_psu,
      source = x$source,
      stringsAsFactors = stringsAsFactors
    )
  } else if (x$stages == 2L) {
    data.frame(
      stages = x$stages,
      varb = x$varb,
      varw = x$varw,
      icc = x$icc,
      var_ratio = x$var_ratio,
      unit_relvar = x$unit_relvar,
      stringsAsFactors = stringsAsFactors
    )
  } else {
    n_b <- x$stages - 1L
    cols <- list(stages = x$stages, varb = x$varb)
    for (k in seq_len(n_b)) cols[[paste0("varw_", k)]] <- unname(x$varw[k])
    for (k in seq_len(n_b)) cols[[paste0("icc_", k)]] <- unname(x$icc[k])
    for (k in seq_len(n_b)) cols[[paste0("var_ratio_", k)]] <- unname(x$var_ratio[k])
    cols$unit_relvar <- x$unit_relvar
    as.data.frame(cols, stringsAsFactors = stringsAsFactors)
  }
  as.data.frame(out, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_cluster <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  out <- if (!is.null(x$domains)) {
    x$domains
  } else {
    data.frame(
      stage = names(x$n),
      n = as.numeric(x$n),
      n_int = as.integer(x),
      stringsAsFactors = stringsAsFactors
    )
  }
  as.data.frame(out, row.names = row.names, optional = optional)
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_power <- function(x, ...) {
  .check_unused_dots(...)
  .coerce_count(.whole_count(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_power <- function(x, ...) {
  .check_unused_dots(...)
  x$n
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_power <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  stringsAsFactors = FALSE,
  validRN = TRUE,
  ...
) {
  .check_unused_dots(...)
  n_pair <- if (length(x$n) == 1L) rep(x$n, 2L) else x$n
  out <- data.frame(
    n1 = unname(n_pair[1L]),
    n2 = unname(n_pair[2L]),
    n1_int = .count_column(.whole_count(n_pair[1L])),
    n2_int = .count_column(.whole_count(n_pair[2L])),
    power = x$power,
    effect = x$effect,
    type = x$type,
    solved = x$solved,
    stringsAsFactors = stringsAsFactors
  )
  as.data.frame(out, row.names = row.names, optional = optional)
}


#' The result itself, for a type with no structured summary yet
#'
#' A list inventory from `summary.default()` says nothing about a design,
#' so these types print as `print()` does until a summary is written.
#' @keywords internal
#' @noRd
.summary_as_print <- function(object, ...) {
  .check_unused_dots(...)
  object
}

#' @rdname print.svyplan
#' @export
summary.svyplan <- function(object, ...) .summary_as_print(object, ...)

#' @rdname print.svyplan
#' @export
summary.svyplan_power <- function(object, ...) .summary_as_print(object, ...)

#' @rdname print.svyplan
#' @export
summary.svyplan_strata <- function(object, ...) .summary_as_print(object, ...)

#' @rdname print.svyplan
#' @export
summary.svyplan_varcomp <- function(object, ...) .summary_as_print(object, ...)

#' @rdname print.svyplan
#' @export
summary.svyplan_overlap <- function(object, ...) .summary_as_print(object, ...)

#' @rdname print.svyplan
#' @export
summary.svyplan_rotation <- function(object, ...) .summary_as_print(object, ...)
