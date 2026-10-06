#' Print svyplan objects
#'
#' Display and coercion methods shared by every result the package returns.
#' `print()` gives the human-readable design, while `as.integer()`,
#' `as.double()`, and `as.data.frame()` extract it in the shapes a
#' downstream package can consume.
#'
#' @param x A svyplan object.
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
#'   Details.
#'
#' @details
#' ## Print and coercion
#'
#' Constrained designs (`n_cluster()`, `n_alloc()`)
#' carry two representations: the continuous mathematical
#' optimum in the top-level fields (`n`, `cv`, `cost`, ...) and the
#' whole-unit field design in `$operational`, whose cost and precision
#' are recomputed from the integer design. `print()` leads with the
#' field design and shows the continuous optimum as a diagnostic.
#'
#' `as.integer(x)` returns the operational design in the same shape as
#' `x$n`: the named integer stage vector for `svyplan_cluster` objects,
#' the operational total for allocation results, and the ceiled scalar
#' otherwise. `as.double(x)` returns the continuous counterpart of the
#' same shape. For `svyplan_strata`, both coercions return the total sample
#' size. Boundary cutpoints remain available in `$boundaries`.
#'
#' `as.data.frame()` returns the tabular form of a result, intended as
#' the stable handoff to downstream packages. For
#' `svyplan_n`: the stratum allocation table (`$detail`) for
#' `n_alloc()` results, the per-domain table (`$domains`, falling back
#' to `$detail`) for `n_multi()` results, and a one-row summary
#' (`n`, `n_int`, `se`, `moe`, `cv`) otherwise. For `svyplan_cluster`:
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

#' Format the design effect for printing
#'
#' `deff = 1` is the most consequential silent default in survey planning,
#' so it is always shown rather than suppressed when left at its default.
#' It is printed without decimals so that an assumed value is visually
#' distinct from a supplied one.
#'
#' Name the stages a cluster design loses units at
#'
#' Three rates acting at three stages are not interchangeable, so the print
#' says which one is set rather than showing a single netted figure.
#' @keywords internal
#' @noRd
.print_stage_rates <- function(p) {
  rates <- c(
    resp_rate_psu = p$resp_rate_psu %||% 1,
    resp_rate_ssu = p$resp_rate_ssu %||% 1,
    resp_rate = p$resp_rate %||% 1
  )
  set <- rates[rates < 1]
  if (length(set) == 0L) {
    return(invisible(NULL))
  }
  cat(sprintf(
    "(%s)\n",
    paste(sprintf("%s = %.2f", names(set), set), collapse = ", ")
  ))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.fmt_deff <- function(deff) {
  if (is.null(deff)) {
    return(NULL)
  }
  if (isTRUE(all.equal(unname(deff), rep(1, length(deff))))) {
    "deff = 1"
  } else {
    sprintf("deff = %s", .fmt_by_stratum(deff, "%.2f"))
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
#' so the header stays one line however many strata there are.
#' @keywords internal
#' @noRd
.fmt_by_stratum <- function(x, fmt) {
  u <- unique(x)
  if (length(u) == 1L) {
    return(sprintf(fmt, u))
  }
  sprintf(paste0(fmt, " to ", fmt, " by stratum"), min(x), max(x))
}

#' @rdname print.svyplan
#' @export
print.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  if (x$type == "alloc") {
    .print_alloc_n(x)
  } else if (x$type == "multi") {
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
  cat(sprintf("n = %s per occasion", .fmt_count_n(ceiling(x$n))))
  if (!is.null(p$resp_rate) && p$resp_rate < 1) {
    cat(sprintf(" (net: %s)", .fmt_count_n(ceiling(x$n * p$resp_rate))))
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
.pooled_scale_label <- function(p) {
  if (is.null(p$p)) "mean scale" else "proportion scale"
}

#' @keywords internal
#' @noRd
.fmt_pooled_estimand <- function(p) {
  if (!is.null(p$p)) {
    return(sprintf("p = %.3g", p$p))
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
  if (all(ov * rho == 0)) {
    return("No between-occasion covariance (overlap x overlap_cor = 0)\n")
  }
  reach <- max(which(ov > 0))
  sprintf(
    "overlap = %.3g, overlap_cor = %.3g at lag 1, shared out to lag %d\n",
    ov[1L],
    rho[1L],
    reach
  )
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
  cat(.fmt_change_n(x$n, p$resp_rate))
  parts <- c(
    .fmt_change_estimand(p),
    if (!is.null(p$moe)) sprintf("moe = %.4g", p$moe),
    if (!is.null(p$rmoe)) sprintf("rmoe = %.3f", p$rmoe),
    if (!is.null(p$cv)) sprintf("cv = %.3f", p$cv),
    .fmt_deff(p$deff),
    if (!is.null(p$resp_rate) && p$resp_rate < 1) {
      sprintf("resp_rate = %s", .fmt_prob(p$resp_rate))
    }
  )
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_overlap(p))
}

#' @keywords internal
#' @noRd
.change_scale_label <- function(p) {
  if (is.null(p$p)) "mean scale" else "proportion scale"
}

#' @keywords internal
#' @noRd
.fmt_change_n <- function(n, resp_rate) {
  net <- !is.null(resp_rate) && resp_rate < 1
  if (length(n) == 1L) {
    out <- sprintf("n = %s per occasion", .fmt_count_n(ceiling(n)))
    if (net) {
      out <- sprintf(
        "%s (net: %s)",
        out,
        .fmt_count_n(ceiling(n * resp_rate))
      )
    }
    return(out)
  }
  out <- sprintf(
    "n = %s then %s",
    .fmt_count_n(ceiling(n[1L])),
    .fmt_count_n(ceiling(n[2L]))
  )
  if (net) {
    out <- sprintf(
      "%s (net: %s then %s)",
      out,
      .fmt_count_n(ceiling(n[1L] * resp_rate)),
      .fmt_count_n(ceiling(n[2L] * resp_rate))
    )
  }
  out
}

#' @keywords internal
#' @noRd
.fmt_change_estimand <- function(p) {
  if (!is.null(p$p)) {
    return(sprintf("p = %.3g to %.3g", p$p[1L], p$p[2L]))
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
    return("No between-occasion covariance (overlap x overlap_cor = 0)\n")
  }
  sprintf(
    "overlap = %.3g, overlap_cor = %.3g (%.1f%% of the independent variance)\n",
    overlap,
    overlap_cor,
    100 * .change_var_share(p)
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
  if (!is.null(resp_rate) && resp_rate < 1) {
    net_n <- ceiling(x$n * resp_rate)
    cat(sprintf(
      "n = %s gross (net: %s)",
      .fmt_count_n(ceiling(x$n)),
      .fmt_count_n(net_n)
    ))
  } else {
    cat(sprintf("n = %s", .fmt_count_n(ceiling(x$n))))
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
  if (length(parts) > 0L) {
    cat(sprintf(" (%s)", paste(parts, collapse = ", ")))
  }
  cat("\n")

  .print_expected_cases(x$expected_cases, x$binding, p$min_cases)

  if (!is.null(x$domains)) {
    cat(sprintf("Domains: %d\n", nrow(x$domains)))
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
    # Both numbers are labelled: the total is what you field, the largest
    # single domain is a different quantity and was once printed as 'n'.
    cat(sprintf(
      "n = %s (%s)\n",
      .fmt_count_n(ceiling(x$n)),
      if (natural) "largest quota / share" else "sum of domain quotas"
    ))
    if (!is.null(x$n_domain_max)) {
      cat(sprintf(
        "Largest single domain = %s (binding: %s)\n",
        .fmt_count_n(ceiling(x$n_domain_max)),
        x$binding
      ))
    }
    cat("\n")
    dom <- x$domains
    dom$.n <- ceiling(dom$.n)
    print(dom, row.names = FALSE, right = FALSE)
  } else {
    cat(sprintf("Multi-indicator sample size\n"))
    cat(sprintf(
      "n = %s (binding: %s)\n",
      .fmt_count_n(ceiling(x$n)),
      x$binding
    ))
    cat("\n")
    d <- x$detail
    d$.n <- ceiling(d$.n)
    d$.binding <- ifelse(d$.binding, "*", "")
    print(d, row.names = FALSE, right = FALSE)
  }
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

#' @keywords internal
#' @noRd
.print_single_cluster <- function(x) {
  cat(sprintf("Optimal %d-stage allocation\n", x$stages))
  op <- x$operational
  n_display <- if (!is.null(op)) op$n else ceiling(x$n)
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
  resp_rate <- (x$params$resp_rate_psu %||% 1) *
    (x$params$resp_rate_ssu %||% 1) *
    (x$params$resp_rate %||% 1)
  if (resp_rate < 1) {
    net_total <- ceiling(total_display * resp_rate)
    cat(sprintf(
      " -> total n = %s (net: %s)\n",
      .fmt_count_n(total_display),
      .fmt_count_n(net_total)
    ))
  } else {
    cat(sprintf(" -> total n = %s\n", .fmt_count_n(total_display)))
  }
  .print_stage_rates(x$params)
  fc <- x$params$fixed_cost
  op_cv <- if (!is.null(op)) op$cv else x$cv
  op_cost <- if (!is.null(op)) op$cost else x$cost
  if (!is.null(fc) && fc > 0) {
    cat(sprintf("cv = %.4f, cost = %.0f (fixed: %.0f)\n", op_cv, op_cost, fc))
  } else {
    cat(sprintf("cv = %.4f, cost = %.0f\n", op_cv, op_cost))
  }
  cont_parts <- vapply(
    seq_along(x$n),
    function(i) {
      sprintf("%s = %s", stage_labels[i], .fmt_continuous_n(x$n[i]))
    },
    character(1L)
  )
  cat(sprintf(
    "continuous optimum: %s (cv = %.4f, cost = %.0f)\n",
    paste(cont_parts, collapse = " | "),
    x$cv,
    x$cost
  ))
  .print_design_df(x)

  if (!is.null(x$domains)) {
    cat(sprintf("Domains: %d\n", nrow(x$domains)))
  }
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
      "Multi-indicator optimal allocation (%d-stage, %d domains%s%s)\n",
      x$stages,
      nrow(x$domains),
      joint_label,
      min_n_label
    ))
    dom <- x$domains
    stage_cols <- names(x$n)
    for (col in stage_cols) {
      dom[[col]] <- ceiling(dom[[col]])
    }
    dom$.total_n <- apply(
      dom[, stage_cols, drop = FALSE],
      1,
      prod
    )
    unrounded <- .fmt_continuous_n(x$total_n)
    fc <- x$params$fixed_cost
    if (!is.null(fc) && fc > 0) {
      cat(sprintf(
        "Total n = %s (unrounded: %s), cost = %.0f (fixed: %.0f)\n",
        .fmt_count_n(sum(dom$.total_n)),
        unrounded,
        x$cost,
        fc
      ))
    } else {
      cat(sprintf(
        "Total n = %s (unrounded: %s)\n",
        .fmt_count_n(sum(dom$.total_n)),
        unrounded
      ))
    }
    dom$.cv <- sprintf("%.4f", dom$.cv)
    dom$.cost <- sprintf("%.0f", dom$.cost)
    cat("\n")
    print(dom, row.names = FALSE, right = FALSE)
  } else {
    cat(sprintf("Multi-indicator optimal allocation (%d-stage)\n", x$stages))
    op <- x$operational
    n_display <- if (!is.null(op)) op$n else ceiling(x$n)
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
    op_cv <- if (!is.null(op)) op$cv else x$cv
    op_cost <- if (!is.null(op)) op$cost else x$cost
    fc <- x$params$fixed_cost
    if (!is.null(fc) && fc > 0) {
      cat(sprintf(
        "worst cv = %.4f, cost = %.0f (fixed: %.0f, binding: %s)\n",
        op_cv,
        op_cost,
        fc,
        x$binding
      ))
    } else {
      cat(sprintf(
        "worst cv = %.4f, cost = %.0f (binding: %s)\n",
        op_cv,
        op_cost,
        x$binding
      ))
    }
    cont_parts <- vapply(
      seq_along(x$n),
      function(i) {
        sprintf("%s = %s", stage_labels[i], .fmt_continuous_n(x$n[i]))
      },
      character(1L)
    )
    cat(sprintf(
      "continuous optimum: %s (cv = %.4f, cost = %.0f)\n",
      paste(cont_parts, collapse = " | "),
      x$cv,
      x$cost
    ))
    cat("\n")
    d <- x$detail
    d$.cv_achieved <- sprintf("%.4f", d$.cv_achieved)
    d$.binding <- ifelse(d$.binding, "*", "")
    print(d, row.names = FALSE, right = FALSE)
  }
}

#' @rdname print.svyplan
#' @export
print.svyplan_prec <- function(x, ...) {
  .check_unused_dots(...)
  if (x$type == "alloc" && identical(x$method, "bethel")) {
    .print_bethel_prec(x)
  } else if (x$type == "multi") {
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
    .pooled_scale_label(p)
  ))
  cat(sprintf("n = %s per occasion", .fmt_count_n(ceiling(p$n))))
  if (!is.null(p$resp_rate) && p$resp_rate < 1) {
    cat(sprintf(" (net: %s)", .fmt_count_n(ceiling(p$n * p$resp_rate))))
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
  cat(sprintf("Sampling precision for change (%s)\n", .change_scale_label(p)))
  cat(.fmt_change_n(p$n, p$resp_rate))
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
    x$type
  )
  p <- x$params

  if (x$type == "cluster") {
    stages <- p$stages
    cat(sprintf("Sampling precision for %d-stage cluster\n", stages))
    n_display <- ceiling(p$n)
    total_display <- prod(n_display)
    stage_labels <- names(p$n)
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
    cat(paste(stage_parts, collapse = " | "))
    cat(sprintf(" -> total n = %s", .fmt_count_n(total_display)))
    resp_rate <- p$resp_rate_psu
    if (!is.null(resp_rate) && resp_rate < 1) {
      cat(sprintf(
        " (net: %s)",
        .fmt_count_n(ceiling(total_display * resp_rate))
      ))
    }
    cat("\n")
  } else {
    notes <- c(x$method, if (!is.null(x$solved)) paste("solved for", x$solved))
    note_label <- if (length(notes) > 0L) {
      sprintf(" (%s)", paste(notes, collapse = ", "))
    } else {
      ""
    }
    cat(sprintf("Sampling precision for %s%s\n", type_label, note_label))
    # An allocation carries one supplied n per stratum. It may be fractional,
    # so report that allocation exactly rather than ceiling each row into a
    # different design.
    n_display <- if (identical(x$type, "alloc")) {
      sum(p$n)
    } else {
      sum(ceiling(p$n))
    }
    cat(sprintf(
      "n = %s",
      if (identical(x$type, "alloc")) {
        .fmt_continuous_n(n_display)
      } else {
        .fmt_count_n(n_display)
      }
    ))
    if (length(p$n) > 1L) {
      cat(sprintf(" (%d strata)", length(p$n)))
    }
    if (identical(x$type, "alloc")) {
      resp_rate <- .alloc_summary_response_rates(x)$combined
      if (any(resp_rate < 1)) {
        cat(sprintf(", expected respondents = %.1f", sum(p$n * resp_rate)))
      }
    } else {
      resp_rate <- p$resp_rate
    }
    if (
      !identical(x$type, "alloc") && !is.null(resp_rate) && all(resp_rate < 1)
    ) {
      cat(sprintf(
        " (net: %s)",
        .fmt_count_n(ceiling(sum(p$n * resp_rate)))
      ))
    }
    cat("\n")
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

#' @keywords internal
#' @noRd
.print_multi_prec <- function(x) {
  cat("Multi-indicator sampling precision\n")
  if (!is.null(x$detail)) {
    cat("\n")
    print(x$detail, row.names = FALSE, right = FALSE)
  }
}

#' @rdname print.svyplan
#' @export
print.svyplan_varcomp <- function(x, ...) {
  .check_unused_dots(...)
  if (!is.null(x$strata)) {
    cat(sprintf(
      "Variance components (%d-stage, %d strata)\n",
      x$stages,
      nrow(x$strata)
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
  cat(sprintf("Variance components (%d-stage)\n", x$stages))
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

#' Report an icc backed out of a published design effect
#'
#' The design effect is not stored, being recoverable from the identity, so
#' it is reformed here at the take that identifies the icc. That take is
#' shown alongside the nominal one whenever the two differ, which is the
#' only visible sign that realized takes varied.
#' @keywords internal
#' @noRd
.print_varcomp_deff <- function(x) {
  take <- x$params$n_per_psu
  nominal <- x$params$n_per_psu_nominal
  cat("Variance components (2-stage, from a design effect)\n")
  cat(sprintf("icc = %.4f\n", x$icc))
  cat(sprintf("var_ratio = %.4f\n", x$var_ratio))
  cat(sprintf(
    "deff = %.4f at n_per_psu = %.4g%s\n",
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
  resp_rate <- p$resp_rate

  is_did <- startsWith(x$type, "did")
  if (length(x$n) == 2L) {
    n1 <- ceiling(x$n[1])
    n2 <- ceiling(x$n[2])
    lab1 <- if (is_did) "n_treat" else "n1"
    lab2 <- if (is_did) "n_control" else "n2"
    if (!is.null(resp_rate) && resp_rate < 1) {
      cat(sprintf(
        "%s = %s, %s = %s (total = %s, net: %s), power = %.3f, effect = %.4f\n",
        lab1,
        .fmt_count_n(n1),
        lab2,
        .fmt_count_n(n2),
        .fmt_count_n(n1 + n2),
        .fmt_count_n(ceiling((n1 + n2) * resp_rate)),
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
  } else if (!is.null(resp_rate) && resp_rate < 1) {
    net_n <- ceiling(x$n * resp_rate)
    cat(sprintf(
      "n = %s (net: %s, per group), power = %.3f, effect = %.4f\n",
      .fmt_count_n(ceiling(x$n)),
      .fmt_count_n(net_n),
      x$power,
      x$effect
    ))
  } else {
    cat(sprintf(
      "n = %s (per group), power = %.3f, effect = %.4f\n",
      .fmt_count_n(ceiling(x$n)),
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
  if (!is.null(resp_rate) && resp_rate < 1) {
    parts <- c(parts, sprintf("resp_rate = %s", .fmt_prob(resp_rate)))
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
  paste0("svyplan_n [n=", ceiling(x$n), ", ", x$type, "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  total <- if (!is.null(x$operational)) {
    x$operational$total_n
  } else {
    prod(ceiling(x$n))
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
  paste0("svyplan_prec [", x$type, ", se=", sprintf("%.4f", x$se[1L]), "]")
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
    paste0(ceiling(x$n[1]), ",", ceiling(x$n[2]))
  } else {
    as.character(ceiling(x$n))
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
      N = p$N %||% Inf,
      deff = p$deff %||% 1,
      resp_rate = p$resp_rate %||% 1,
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

#' The change an interval on a change is centred on
#'
#' A change sized from `moe` alone has no known level, and the interval
#' would then be centred on nothing. The message names `change` rather than
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

#' The level an interval on a pooled estimate is centred on
#'
#' A size solved from `moe` alone knows the spread but not the level, and
#' the interval would then be centred on nothing. The message names `mu`
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
      N = p$N %||% Inf,
      deff = p$deff %||% 1,
      resp_rate = p$resp_rate %||% 1,
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
    return(as.integer(x$operational$n))
  }
  as.integer(ceiling(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  x$n
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  n <- if (!is.null(x$operational)) x$operational$n else ceiling(x$n)
  setNames(as.integer(n), names(n))
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
  out <- if (identical(x$type, "alloc")) {
    x$detail
  } else if (identical(x$type, "multi")) {
    x$domains %||% x$detail
  } else {
    data.frame(
      n = x$n,
      n_int = as.integer(ceiling(x$n)),
      se = x$se,
      moe = x$moe,
      cv = x$cv,
      stringsAsFactors = stringsAsFactors
    )
  }
  as.data.frame(out, row.names = row.names, optional = optional)
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
      stringsAsFactors = stringsAsFactors
    )
  }

  as.data.frame(out, row.names = row.names, optional = optional)
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
      n_per_psu = x$params$n_per_psu,
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
    data.frame(
      stages = x$stages,
      varb = x$varb,
      varw_psu = unname(x$varw[1L]),
      varw_ssu = unname(x$varw[2L]),
      icc_psu = unname(x$icc[1L]),
      icc_ssu = unname(x$icc[2L]),
      var_ratio_psu = unname(x$var_ratio[1L]),
      var_ratio_ssu = unname(x$var_ratio[2L]),
      unit_relvar = x$unit_relvar,
      stringsAsFactors = stringsAsFactors
    )
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
  as.integer(ceiling(x$n))
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
    n1_int = as.integer(ceiling(n_pair[1L])),
    n2_int = as.integer(ceiling(n_pair[2L])),
    power = x$power,
    effect = x$effect,
    type = x$type,
    solved = x$solved,
    stringsAsFactors = stringsAsFactors
  )
  as.data.frame(out, row.names = row.names, optional = optional)
}

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
#' d <- design_effect(icc = 0.03, n_per_psu = 20,
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
  cat(sprintf("n = %d, cv = %.4f%s\n", ceiling(x$n), x$cv, alloc))
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
  cluster <- !is.null(detail) && "n_psu" %in% names(detail)
  H <- if (!is.null(detail)) nrow(detail) else NA
  cat(sprintf("Stratum allocation (%s", alloc))
  if (cluster) {
    cat(", two-stage")
  }
  if (!is.na(H)) {
    cat(sprintf(", %d strata", H))
  }
  cat(")\n")
  op <- x$operational
  if (!is.null(op)) {
    cat(sprintf("field design: n = %d", op$n))
    if (cluster) {
      cat(sprintf(", n_psu = %d", sum(detail$n_psu_int)))
    }
    rates <- .alloc_summary_response_rates(x)$combined
    if (any(rates < 1)) {
      cat(sprintf(", expected respondents = %.1f", sum(detail$n_int * rates)))
    }
    if (!is.na(op$cv)) {
      cat(sprintf(", cv = %.4f", op$cv))
    }
    if (!is.na(op$cost)) {
      cat(sprintf(", cost = %.0f", op$cost))
    }
    cat("\n")
    cat(sprintf("continuous optimum: n = %s", .fmt_continuous_n(x$n)))
    if (!is.na(x$cv)) {
      cat(sprintf(", cv = %.4f", x$cv))
    }
    if (!is.na(x$se)) {
      cat(sprintf(", se = %.4f", x$se))
    }
    cat("\n")
  } else {
    cat(sprintf("n = %d", ceiling(x$n)))
    if (cluster) {
      cat(sprintf(", n_psu = %d", sum(detail$n_psu_int)))
    }
    rates <- .alloc_summary_response_rates(x)$combined
    if (any(rates < 1)) {
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
  if (!is.null(p$min_n_stratum) && p$min_n_stratum > 0) {
    parts <- c(parts, sprintf("min_n_stratum = %g", p$min_n_stratum))
  }
  resp_rate <- p$resp_rate
  if (!is.null(resp_rate) && any(resp_rate < 1)) {
    parts <- c(
      parts,
      sprintf("resp_rate = %s", .fmt_by_stratum(resp_rate, "%.2f"))
    )
  }
  parts <- c(parts, .fmt_deff(p$deff))
  if (length(parts) > 0L) {
    cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))
  }
  .print_psu_fraction_note(detail)
  .print_design_df(x)
  .print_domains_block(x$domains)
}

#' Disclose an appreciable first-stage sampling fraction
#'
#' `N_psu` bounds the allocation but activates no first-stage correction, so a
#' design taking a large share of the available PSUs is planned conservatively.
#' The ultimate-unit `1 - n / N` does attenuate the between-PSU component, but
#' it is smaller than the first-stage fraction this note reports, so the
#' between-PSU term is still overstated. That is a property of the result worth
#' seeing, but not an event to act on, so it is disclosed here rather than
#' warned about. `predict()` re-runs the allocation once per grid row, and a
#' warning would repeat with it.
#' @keywords internal
#' @noRd
.print_psu_fraction_note <- function(detail, threshold = 0.1) {
  if (is.null(detail) || !".psu_frac" %in% names(detail)) {
    return(invisible(NULL))
  }
  hit <- which(is.finite(detail$.psu_frac) & detail$.psu_frac >= threshold)
  if (length(hit) == 0L) {
    return(invisible(NULL))
  }
  bound <- !is.na(detail$.bound_source) & detail$.bound_source == "N_psu"
  cat(sprintf(
    paste0(
      "note: samples %s of the available PSUs in %s%s\n",
      "      precision uses a with-replacement first stage and may be ",
      "conservative\n"
    ),
    paste0(
      format(round(100 * detail$.psu_frac[hit]), trim = TRUE),
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

#' Print a joint constrained allocation
#' @keywords internal
#' @noRd
.print_bethel_n <- function(x) {
  opt <- x$optimization
  op <- x$operational
  budget_mode <- identical(x$params$mode, "budget_objective")
  status <- opt$classification %||% "unknown"
  certainty <- opt$certainty
  cat(sprintf(
    "Joint constrained allocation (Bethel%s)\n",
    if (is.null(certainty)) "" else ", PSU register"
  ))
  # Report exceptions, not confirmations: a line that only ever says
  # "nothing went wrong" buries the lines that do carry news.
  if (!identical(status, "optimal")) {
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
    "field design: n = %d, cost = %.0f%s\n",
    op$n,
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
  if (!is.null(certainty)) {
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
    # A register with no self-consistent classification is a fact about the
    # plan rather than detail, so it stays on the block, as a panel short of
    # its target does.
    if (!isTRUE(certainty$fixed_point)) {
      cat("no self-consistent classification: the plan meets every target\n")
    }
  }
  continuous_cost <- x$params$achieved$cost
  increase <- if (continuous_cost > 0) {
    100 * (op$cost / continuous_cost - 1)
  } else {
    NA_real_
  }
  cat(sprintf(
    "continuous optimum: n = %s, cost = %.0f%s\n",
    .fmt_continuous_n(x$n),
    continuous_cost,
    if (!is.na(increase) && abs(increase) >= 0.005) {
      sprintf(" (integerizing costs %+.2f%%)", increase)
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
    if (!isTRUE(opt$budget_binding)) {
      cat(
        "the budget is not binding: the allocation sits at its upper bounds\n"
      )
    } else if (is.finite(opt$budget_sensitivity %||% NA_real_)) {
      cat(sprintf(
        "one more unit of budget changes the objective by %.4g\n",
        opt$budget_sensitivity
      ))
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
  if (!is.null(certainty)) {
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

#' Summarize stratified allocation results
#'
#' `summary()` separates a classic [n_alloc()] or [prec_alloc()] result into
#' the overall answer, the per-stratum allocation, achieved precision,
#' allocation bounds, domains, and assumptions. The design summary evaluates
#' the whole-unit field allocation. The precision summary evaluates the
#' supplied allocation exactly. Generalized Bethel allocations return a
#' separate diagnostic summary containing their constraint, objective,
#' optimization, and bound tables.
#'
#' @param object An allocation result from [n_alloc()] or [prec_alloc()].
#' @param x A `summary.svyplan_alloc` or `summary.svyplan_bethel` object.
#' @param ... Additional arguments are not supported for allocation summaries.
#'   For other `svyplan_n` and `svyplan_prec` results they are passed to the
#'   default summary method.
#'
#' @return For a classic allocation, `summary()` returns a
#'   `summary.svyplan_alloc` object with fields `kind`, `question`, `method`,
#'   `mode`, `overall`, `continuous`, `allocation`, `precision`, `bounds`,
#'   `domains`, and `assumptions`. For a generalized allocation it returns the
#'   `summary.svyplan_bethel` fields described in Details. Numeric quantities
#'   remain unformatted. Formatting is applied only by the summary print
#'   methods. `print()` returns the summary invisibly.
#'
#' @details
#' For `n_alloc()`, `overall`, `allocation`, `precision`, and `domains`
#' describe `$detail$n_int` and the operational cluster takes when present.
#' The mathematical optimum remains in `continuous`. For `prec_alloc()`, the
#' same fields describe `$detail$n`, the allocation actually supplied to the
#' precision calculation. `continuous` and `bounds` are `NULL`.
#'
#' Stratum standard errors, margins of error, and CVs do not add to their
#' overall counterparts. `variance_share` is the stratum's contribution to
#' the overall variance and sums to one when that variance is positive.
#' Domain rows are separate subpopulation assessments, not contributions to
#' add to the overall row.
#'
#' A generalized Bethel summary has fields `kind`, `question`, `mode`,
#' `stages`, `status`, `overall`, `continuous`, `allocation`, `constraints`,
#' `operational_constraints`, `objective`, `operational_objective`, `bounds`,
#' `optimization`, and `assumptions`. For a fitted design the headline and
#' allocation describe the integer field recommendation, while `continuous`
#' and `constraints` retain the mathematical optimum. For a precision
#' assessment they describe the supplied allocation exactly.
#'
#' Constraint rows are simultaneous requirements, not an additive
#' decomposition. `.ratio` is achieved divided by target, `.residual` is that
#' ratio minus one, and `.sensitivity` is the local change in minimum cost per
#' unit relaxation of the stated target. Objective `.contribution` and
#' `.share` columns are genuinely additive across components.
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("Urban", "Rural", "Remote"),
#'   N = c(12000, 30000, 8000),
#'   sd = c(5, 8, 12),
#'   mean = c(20, 18, 15),
#'   unit_cost = c(1, 1.4, 3)
#' )
#' allocation <- n_alloc(frame, n = 1200)
#' summary(allocation)
#' summary(prec_alloc(allocation, n = allocation$detail$n_int))
#'
#' @seealso [n_alloc()] and [prec_alloc()], which build these objects,
#'   [strata_bound()] for constructing the strata to allocate over, and
#'   [design_df()] for the degrees of freedom the allocation leaves.
#'
#' @name summary.svyplan_alloc
NULL

#' @rdname summary.svyplan_alloc
#' @export
summary.svyplan_n <- function(object, ...) {
  if (identical(object$type, "alloc")) {
    .check_unused_dots(...)
    if (identical(object$method, "bethel")) {
      return(.summary_bethel_alloc(object, kind = "design"))
    }
    return(.summary_classic_alloc(object, kind = "design"))
  }
  base::summary.default(object, ...)
}

#' @rdname summary.svyplan_alloc
#' @export
summary.svyplan_prec <- function(object, ...) {
  if (identical(object$type, "alloc")) {
    .check_unused_dots(...)
    if (identical(object$method, "bethel")) {
      return(.summary_bethel_alloc(object, kind = "assessment"))
    }
    return(.summary_classic_alloc(object, kind = "assessment"))
  }
  base::summary.default(object, ...)
}

#' Stages of the planned design, not of the solve
#'
#' A register fit is solved as one stage, its clustering handed over as a
#' design effect, but the plan it describes has two.
#' @keywords internal
#' @noRd
.bethel_plan_stages <- function(object, problem) {
  if (is.null(object$params$psu)) problem$stages else 2L
}

#' Build a structured diagnostic summary for a generalized allocation
#' @keywords internal
#' @noRd
.summary_bethel_alloc <- function(object, kind) {
  p <- object$params
  problem <- p$problem
  budget <- p[["budget"]]
  budget_residual <- p[["budget_residual"]]
  if (is.null(problem) || is.null(problem$stratum_ids)) {
    stop(
      "generalized allocation result is missing its planning problem",
      call. = FALSE
    )
  }
  design <- identical(kind, "design")
  constraints <- if (design) object$constraints else object$detail
  operational_constraints <- if (design) {
    object$operational$constraints
  } else {
    NULL
  }
  mode <- if (design) {
    p$mode
  } else if (!is.null(object$objective) && !is.null(budget)) {
    "budget_objective"
  } else if (!is.null(object$objective)) {
    "objective"
  } else {
    "targets"
  }
  allocation <- .bethel_summary_allocation(object, design)
  bounds <- .bethel_summary_bounds(object, allocation, design)
  objective <- object$objective
  operational_objective <- if (design) object$operational$objective else NULL

  primary_constraints <- if (design) operational_constraints else constraints
  primary_objective <- if (design) operational_objective else objective
  overall <- .bethel_summary_overall(
    n = if (design) object$operational$n else p$achieved$n,
    cost = if (design) object$operational$cost else p$achieved$cost,
    constraints = primary_constraints,
    bounds = bounds,
    objective = primary_objective,
    objective_value = if (design) {
      object$operational$objective_value
    } else {
      object$objective_value
    },
    budget = budget,
    budget_residual = if (design) {
      object$operational$budget_residual
    } else {
      budget_residual
    }
  )
  continuous <- if (design) {
    .bethel_summary_overall(
      n = object$n,
      cost = p$achieved$cost,
      constraints = constraints,
      bounds = bounds,
      objective = objective,
      objective_value = object$objective_value,
      budget = budget,
      budget_residual = object$optimization$budget_residual
    )
  } else {
    NULL
  }

  # `certainty` is attached only when there is one, so a plan with no register
  # keeps the schema it has always had.
  out <- list(
    kind = kind,
    question = .bethel_summary_question(mode, kind, budget),
    mode = mode,
    stages = .bethel_plan_stages(object, problem),
    status = if (design) object$optimization$classification else "assessed",
    overall = overall,
    continuous = continuous,
    allocation = allocation,
    constraints = constraints,
    operational_constraints = operational_constraints,
    objective = objective,
    operational_objective = operational_objective,
    bounds = bounds,
    optimization = if (design) {
      .bethel_summary_optimization(object)
    } else {
      NULL
    },
    assumptions = .bethel_summary_assumptions(object)
  )
  cert <- .bethel_summary_certainty(object)
  if (!is.null(cert)) {
    out$certainty <- cert
    sources <- c("threshold", "operational", "supplied", "orbit")
    out$certainty_sources <- vapply(
      sources,
      function(src) sum(object$psu$.certainty_source %in% src),
      integer(1)
    )
  }
  structure(out, class = c("summary.svyplan_bethel", "list"))
}

#' The certainty split by stratum, which the printed line points at
#'
#' Counts are PSUs available in each part. The threshold and the design effect
#' are read from the plan as returned, so the table describes the answer rather
#' than any iterate on the way to it.
#' @keywords internal
#' @noRd
.bethel_summary_certainty <- function(object) {
  if (is.null(object$optimization$certainty) || is.null(object$psu)) {
    return(NULL)
  }
  detail <- object$detail
  out <- data.frame(
    stratum = detail$stratum,
    certainty = detail$n_psu_certain,
    remainder = detail$n_psu_rest,
    threshold = detail$threshold,
    stringsAsFactors = FALSE
  )
  # The threshold is the cutoff's, so a relaxed rule is shown beside it.
  cutoff <- .psu_cutoff(object$params$certainty_cutoff, object$params$frame)
  if (any(cutoff < 1)) out$cutoff <- cutoff
  if (!is.null(detail$n_zone)) out$zones <- detail$n_zone
  # The design effects the split produced already have a home, one row per
  # stratum and constraint in the resolved planning inputs below, so they are
  # not repeated here.
  out
}

#' Headline values for one numerical version of a generalized allocation
#' @keywords internal
#' @noRd
.bethel_summary_overall <- function(
  n,
  cost,
  constraints,
  bounds,
  objective,
  objective_value,
  budget = NULL,
  budget_residual = NULL
) {
  n_constraints <- if (is.null(constraints)) 0L else nrow(constraints)
  n_pass <- if (n_constraints == 0L) 0L else sum(constraints$.pass)
  n_bounds_bad <- if (is.null(bounds)) {
    0L
  } else {
    status <- if ("field_status" %in% names(bounds)) {
      bounds$field_status
    } else {
      bounds$status
    }
    sum(grepl("violation", status), na.rm = TRUE)
  }
  list(
    n = n,
    cost = cost,
    n_constraints = n_constraints,
    n_pass = n_pass,
    n_violated = n_constraints - n_pass,
    all_pass = n_constraints == 0L || n_pass == n_constraints,
    n_bound_violations = n_bounds_bad,
    objective_value = objective_value,
    objective_components = if (is.null(objective)) 0L else nrow(objective),
    budget = budget,
    budget_residual = budget_residual
  )
}

#' State which generalized allocation question was answered
#' @keywords internal
#' @noRd
.bethel_summary_question <- function(mode, kind, budget) {
  if (identical(kind, "assessment")) {
    if (identical(mode, "budget_objective")) {
      return(
        "assess a supplied joint allocation against its targets, objective, and budget"
      )
    }
    if (identical(mode, "objective")) {
      return(
        "assess a supplied joint allocation's precision and weighted objective"
      )
    }
    return("assess a supplied joint allocation against every precision target")
  }
  if (identical(mode, "budget_objective")) {
    return(sprintf(
      "find the best joint allocation affordable within budget = %.6g",
      budget
    ))
  }
  "find the cheapest joint allocation meeting every precision target"
}

#' Per-stratum decisions on their public ultimate-unit scale
#' @keywords internal
#' @noRd
.bethel_summary_allocation <- function(object, design) {
  p <- object$params
  problem <- p$problem
  take <- problem$stage$take
  if (design) {
    d <- object$detail
    n_continuous <- d$n
    n_primary <- d$n_int
  } else {
    n_continuous <- NULL
    n_primary <- as.numeric(p$n)
  }
  decision <- n_primary / take
  out <- data.frame(
    stratum = problem$stratum_ids,
    population = problem$population_N,
    stringsAsFactors = FALSE
  )
  if (design) {
    out$n_continuous <- n_continuous
  }
  out[[if (design) "n_field" else "n_supplied"]] <- n_primary
  out$weight <- problem$population_N / n_primary
  out$cost <- problem$cost * decision
  # A register's solver prices interviews at the marginal cost of its held
  # classification, so its certainty visits are added here, and the field
  # design is priced as it is drawn.
  if (!is.null(p$psu) && all(c("cost_psu", "cost_ssu") %in% names(p$frame))) {
    cost_psu <- rep_len(as.numeric(p$frame$cost_psu), nrow(out))
    cost_ssu <- rep_len(as.numeric(p$frame$cost_ssu), nrow(out))
    out$cost <- if (design) {
      d <- object$detail
      (d$n_psu_certain + d$n_psu_draw) * cost_psu + d$n_int * cost_ssu
    } else {
      p$certainty$n_psu_certain * cost_psu + out$cost
    }
  }
  if (problem$stages > 1L) {
    out$population_psu <- problem$stage$N_psu
    if (problem$stages == 3L) {
      out$population_ssu <- problem$stage$N_ssu
    }
    if (design) {
      out$psu_continuous <- n_continuous / take
    }
    out[[if (design) "psu_field" else "psu_supplied"]] <- decision
    out$take_per_psu <- problem$stage$n_per_psu
    if (problem$stages == 3L) {
      out$take_per_ssu <- problem$stage$n_per_ssu
    }
  }
  out
}

#' Bound positions and violations on the ultimate-unit scale
#' @keywords internal
#' @noRd
.bethel_summary_bounds <- function(object, allocation, design) {
  p <- object$params
  problem <- p$problem
  lower <- problem$lower * problem$stage$take
  upper <- problem$upper * problem$stage$take
  if (design) {
    continuous <- object$detail$n
    field <- object$detail$n_int
    data.frame(
      stratum = problem$stratum_ids,
      lower = lower,
      continuous = continuous,
      field = field,
      upper = upper,
      continuous_status = .bethel_bound_status(continuous, lower, upper),
      field_status = .bethel_bound_status(field, lower, upper),
      stringsAsFactors = FALSE
    )
  } else {
    supplied <- allocation$n_supplied
    data.frame(
      stratum = problem$stratum_ids,
      lower = lower,
      supplied = supplied,
      upper = upper,
      lower_violation = supplied < lower - p$feasibility_tolerance,
      upper_violation = supplied > upper + p$feasibility_tolerance,
      status = .bethel_bound_status(supplied, lower, upper),
      pass = !(supplied < lower - p$feasibility_tolerance |
        supplied > upper + p$feasibility_tolerance),
      stringsAsFactors = FALSE
    )
  }
}

#' Label one allocation's relation to each pair of bounds
#' @keywords internal
#' @noRd
.bethel_bound_status <- function(value, lower, upper, tolerance = 1e-6) {
  tol <- tolerance * pmax(1, abs(value), abs(lower), abs(upper))
  out <- rep("interior", length(value))
  out[value < lower - tol] <- "lower violation"
  out[value > upper + tol] <- "upper violation"
  fixed <- abs(upper - lower) <= tol
  out[fixed & abs(value - lower) <= tol] <- "fixed"
  out[!fixed & abs(value - lower) <= tol] <- "lower"
  out[!fixed & abs(value - upper) <= tol] <- "upper"
  out
}

#' Make solver diagnostics readable without discarding their raw values
#' @keywords internal
#' @noRd
.bethel_summary_optimization <- function(object) {
  opt <- object$optimization
  ids <- object$params$problem$constraint_ids
  solver_ids <- c(
    ids,
    if (length(opt$constraint_scale %||% numeric(0)) > length(ids)) {
      "objective"
    }
  )
  active <- opt$active_precision %||% integer(0)
  active_precision <- ifelse(
    active <= length(ids),
    ids[active],
    "objective"
  )
  diagnostic_names <- c(
    "primal_residual",
    "projected_dual_residual",
    "stationarity_residual",
    "complementarity_residual",
    "relative_duality_gap"
  )
  present <- diagnostic_names[diagnostic_names %in% names(opt)]
  diagnostics <- data.frame(
    diagnostic = present,
    value = unname(unlist(opt[present])),
    tolerance = rep(opt$kkt_tolerance %||% NA_real_, length(present)),
    stringsAsFactors = FALSE
  )
  diagnostics$pass <- diagnostics$value <= diagnostics$tolerance
  multiplier_ids <- c(
    ids,
    if (
      length(
        opt$multiplier_identifiable %||%
          logical(0)
      ) >
        length(ids)
    ) {
      "objective"
    }
  )
  identifiable <- opt$multiplier_identifiable %||% logical(0)
  if (length(identifiable) > 0L) {
    names(identifiable) <- multiplier_ids
  }
  constraint_scale <- opt$constraint_scale %||% numeric(0)
  constraint_residual <- opt$constraint_residual %||% numeric(0)
  if (length(constraint_scale) > 0L) {
    names(constraint_scale) <- solver_ids
  }
  if (length(constraint_residual) > 0L) {
    names(constraint_residual) <- solver_ids
  }
  list(
    classification = opt$classification,
    converged = opt$converged,
    certainty = opt$certainty,
    message = opt$message,
    iterations = opt$iterations,
    convergence_code = opt$convergence_code,
    feasibility_tolerance = opt$feasibility_tolerance,
    kkt_tolerance = opt$kkt_tolerance,
    diagnostics = diagnostics,
    active_precision = unname(active_precision),
    active_lower = object$params$problem$stratum_ids[
      opt$active_lower %||% integer(0)
    ],
    active_upper = object$params$problem$stratum_ids[
      opt$active_upper %||% integer(0)
    ],
    multiplier_identifiable = identifiable,
    continuous_cost = opt$cost,
    cost_scale = opt$cost_scale,
    constraint_scale = constraint_scale,
    constraint_residual = constraint_residual,
    dual_value = opt$dual_value,
    budget = opt$budget,
    budget_binding = opt$budget_binding,
    budget_residual = opt$budget_residual,
    budget_sensitivity = opt$budget_sensitivity,
    objective_bound = opt$objective_bound,
    objective_multiplier = opt$objective_multiplier,
    root_iterations = opt$root_iterations
  )
}

#' Resolved numerical assumptions behind each target-stratum contribution
#' @keywords internal
#' @noRd
.bethel_summary_assumptions <- function(object) {
  p <- object$params
  problem <- p$problem
  K <- length(problem$constraint_ids)
  H <- length(problem$stratum_ids)
  model <- NULL
  if (K > 0L) {
    grid <- expand.grid(
      stratum_index = seq_len(H),
      constraint_index = seq_len(K),
      KEEP.OUT.ATTRS = FALSE
    )
    keep <- problem$membership[cbind(grid$stratum_index, grid$constraint_index)]
    grid <- grid[keep, , drop = FALSE]
    at <- cbind(grid$stratum_index, grid$constraint_index)
    model <- data.frame(
      stratum = problem$stratum_ids[grid$stratum_index],
      constraint = problem$constraint_ids[grid$constraint_index],
      mean = problem$mean_hk[at],
      variance = problem$var_hk[at],
      deff = problem$deff_hk[at],
      response_rate = problem$resp_hk[at],
      stringsAsFactors = FALSE
    )
  }
  list(
    alpha = p$alpha,
    df = problem$df,
    stages = .bethel_plan_stages(object, problem),
    # The solver saw a one-stage problem because the certainty loop hands it
    # the clustering as a design effect, so the model it recorded describes
    # the solve rather than the plan.
    variance_model = if (is.null(object$optimization$certainty)) {
      problem$variance_model
    } else {
      "two_part_certainty_wald"
    },
    feasibility_tolerance = p$feasibility_tolerance,
    min_n_stratum = p$min_n_stratum,
    fixed_takes = if (problem$stages > 1L) {
      data.frame(
        stratum = problem$stratum_ids,
        n_per_psu = problem$stage$n_per_psu,
        n_per_ssu = if (problem$stages == 3L) {
          problem$stage$n_per_ssu
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    },
    model = model
  )
}

#' Build one structured summary for design and assessment directions
#' @keywords internal
#' @noRd
.summary_classic_alloc <- function(object, kind) {
  detail <- object$detail
  if (is.null(detail) || !all(c("stratum", "N", "n") %in% names(detail))) {
    stop(
      "classic allocation result is missing its per-stratum detail",
      call. = FALSE
    )
  }

  design <- identical(kind, "design")
  n_used <- if (design) detail$n_int else detail$n
  state <- .alloc_summary_state(object, n_used, operational = design)
  overall <- .alloc_summary_overall(state, n_used)
  if (design) {
    df <- tryCatch(
      suppressWarnings(as.double(design_df(object))),
      error = function(e) NA_real_
    )
    overall$design_df <- df
  }

  continuous <- NULL
  if (design) {
    continuous <- list(
      n = object$n,
      expected_respondents = sum(detail$n * state$response_rate),
      cost = if (all(is.na(state$cost_h))) {
        NA_real_
      } else {
        object$params$achieved$cost
      },
      se = object$se,
      moe = object$moe,
      rmoe = object$rmoe,
      cv = object$cv
    )
  }

  domains <- .alloc_summary_domains(object, state, n_used)
  if (!is.null(domains)) {
    overall$worst_domain_cv <- max(domains$cv)
    if (design && !is.null(object$domains)) {
      continuous$worst_domain_cv <- max(object$domains$.cv)
    }
  }

  mode <- if (design) object$params$mode else "assessment"
  method <- if (design) object$method else NULL
  structure(
    list(
      kind = kind,
      question = .alloc_summary_question(object, kind),
      method = method,
      mode = mode,
      overall = overall,
      continuous = continuous,
      allocation = .alloc_summary_allocation(
        object,
        state,
        n_used,
        kind = kind
      ),
      precision = .alloc_summary_precision(state, n_used),
      bounds = if (design) .alloc_summary_bounds(detail, n_used) else NULL,
      domains = domains,
      assumptions = list(
        alpha = object$params$alpha,
        deff = state$deff,
        response_rate = state$response_rate,
        unit_response_rate = state$unit_response_rate,
        psu_response_rate = state$psu_response_rate,
        df = object$params$df,
        min_n_stratum = object$params$min_n_stratum,
        alloc_q = object$params$alloc_q
      )
    ),
    class = c("summary.svyplan_alloc", "list")
  )
}

#' Resolve the variance, response, and cost model at the displayed allocation
#' @keywords internal
#' @noRd
.alloc_summary_state <- function(object, n_h, operational) {
  p <- object$params
  d <- object$detail
  frame <- p$frame
  H <- nrow(d)
  deff <- .alloc_resolve_h(p$deff, frame, "deff", H, .check_deff_h)
  rates <- .alloc_summary_response_rates(object)
  unit_response <- rates$unit
  cluster <- .alloc_is_cluster(frame)
  psu_response <- rates$psu
  response <- rates$combined
  S_h <- d$sd
  cost_h <- d$unit_cost
  N_fpc <- d$N

  if (cluster) {
    take <- if (operational && "n_per_psu_int" %in% names(d)) {
      d$n_per_psu_int
    } else {
      d$n_per_psu
    }
    var_ratio <- frame[["var_ratio_psu"]] %||% rep(1, H)
    icc <- frame$icc_psu
    fpc_mode <- p$fpc %||% "unit"
    responding_take <- take * unit_response
    # Rebuilt on the take being displayed, so the operational block reads the
    # correction the operational design gets. See .alloc_cluster_prep().
    within_frac <- if (identical(fpc_mode, "stage")) {
      pmin(1, responding_take / (d$N / frame$N_psu))
    } else {
      rep(0, H)
    }
    inflation <- var_ratio *
      (icc * responding_take + (1 - within_frac) * (1 - icc))
    S_h <- d$sd * sqrt(inflation)
    N_fpc <- switch(
      fpc_mode,
      unit = d$N,
      none = rep(Inf, H),
      stage = ifelse(icc > 0, frame$N_psu * inflation / (var_ratio * icc), Inf)
    )
    if (all(c("cost_psu", "cost_ssu") %in% names(frame))) {
      cost_h <- frame$cost_psu / take + frame$cost_ssu
    } else {
      cost_h <- rep(NA_real_, H)
    }
  }

  mean_h <- if ("mean" %in% names(d)) d$mean else rep(NA_real_, H)
  metrics <- .alloc_metrics(
    N_h = d$N,
    S_h = S_h,
    mean_h = mean_h,
    n_h = n_h,
    alpha = p$alpha,
    deff = deff,
    resp_rate = response,
    cost_h = cost_h,
    df = p$df,
    N_fpc = N_fpc
  )
  list(
    metrics = metrics,
    stratum = d$stratum,
    population = d$N,
    N_fpc = N_fpc,
    sd = S_h,
    mean = mean_h,
    deff = deff,
    unit_response_rate = unit_response,
    psu_response_rate = psu_response,
    response_rate = response,
    cost_h = cost_h,
    cluster = cluster
  )
}

#' Resolve the stage rates whose product leaves responding ultimate units
#' @keywords internal
#' @noRd
.alloc_summary_response_rates <- function(object) {
  p <- object$params
  frame <- p$frame
  H <- nrow(object$detail)
  unit <- .alloc_resolve_h(
    p$resp_rate,
    frame,
    "resp_rate",
    H,
    .check_resp_rate_h
  )
  psu <- if (.alloc_is_cluster(frame)) {
    .alloc_resolve_h(1, frame, "resp_rate_psu", H, .check_resp_rate_h)
  } else {
    rep(1, H)
  }
  list(unit = unit, psu = psu, combined = unit * psu)
}

#' Overall fields on the same numerical basis as the stratum tables
#' @keywords internal
#' @noRd
.alloc_summary_overall <- function(state, n_h) {
  m <- state$metrics
  list(
    n = sum(n_h),
    expected_respondents = sum(n_h * state$response_rate),
    cost = m$cost,
    se = m$se,
    moe = m$moe,
    rmoe = m$rmoe,
    cv = m$cv
  )
}

#' State the allocation question rather than only its algorithm
#' @keywords internal
#' @noRd
.alloc_summary_question <- function(object, kind) {
  if (identical(kind, "assessment")) {
    return("assess a supplied stratified allocation")
  }
  p <- object$params
  method <- paste0(
    toupper(substr(object$method, 1L, 1L)),
    substring(object$method, 2L)
  )
  switch(
    p$mode,
    n = sprintf(
      "distribute a fixed sample of %g using %s allocation",
      p$n,
      method
    ),
    cv = sprintf(
      "find the smallest %s allocation attaining cv = %.4g",
      method,
      p$cv
    ),
    budget = sprintf(
      "find the best %s allocation within budget = %.6g",
      method,
      p$budget
    ),
    sprintf("construct a %s allocation", method)
  )
}

#' Field or supplied allocation table
#' @keywords internal
#' @noRd
.alloc_summary_allocation <- function(object, state, n_h, kind) {
  d <- object$detail
  out <- data.frame(
    stratum = d$stratum,
    population = d$N,
    response_rate = state$response_rate,
    unit_cost = state$cost_h,
    stringsAsFactors = FALSE
  )
  if (identical(kind, "design")) {
    out$n_continuous <- d$n
    out$n_field <- n_h
  } else {
    out$n_supplied <- n_h
  }
  if (state$cluster) {
    if (identical(kind, "design")) {
      out$n_psu_continuous <- d$n_psu
      out$n_psu_field <- d$n_psu_int
      out$n_per_psu_continuous <- d$n_per_psu
      out$n_per_psu_field <- d$n_per_psu_int
    } else {
      out$n_psu <- d$n_psu
      out$n_per_psu <- d$n_per_psu
    }
  }
  out$expected_respondents <- n_h * state$response_rate
  out$weight <- d$N / n_h
  out$cost <- n_h * state$cost_h
  if (all(is.na(out$unit_cost))) {
    out$unit_cost <- NULL
  }
  if (all(is.na(out$cost))) {
    out$cost <- NULL
  }
  out
}

#' Per-stratum precision at the field or supplied allocation
#' @keywords internal
#' @noRd
.alloc_summary_precision <- function(state, n_h) {
  m <- state$metrics
  data.frame(
    stratum = state$stratum,
    effective_n = n_h * state$response_rate / state$deff,
    se = m$se_h,
    moe = m$moe_h,
    rmoe = m$rmoe_h,
    cv = m$cv_h,
    variance_share = m$share_h,
    stringsAsFactors = FALSE
  )
}

#' Bound table for a fitted classic allocation
#' @keywords internal
#' @noRd
.alloc_summary_bounds <- function(detail, n_h) {
  data.frame(
    stratum = detail$stratum,
    lower = detail$.lower,
    field = n_h,
    upper = detail$.upper,
    binding = detail$.binding,
    source = detail$.bound_source,
    stringsAsFactors = FALSE
  )
}

#' Recompute domains on the same basis as the headline
#' @keywords internal
#' @noRd
.alloc_summary_domains <- function(object, state, n_h) {
  p <- object$params
  idx_list <- p$domain_idx
  if (is.null(idx_list) || length(idx_list) == 0L) {
    return(NULL)
  }
  frame <- p$frame
  ids <- p$domain_cols %||% character(0)
  rows <- lapply(seq_along(idx_list), function(i) {
    idx <- idx_list[[i]]
    m <- .alloc_metrics(
      N_h = state$population[idx],
      S_h = state$sd[idx],
      mean_h = state$mean[idx],
      n_h = n_h[idx],
      alpha = p$alpha,
      deff = state$deff[idx],
      resp_rate = state$response_rate[idx],
      cost_h = state$cost_h[idx],
      df = p$df,
      N_fpc = state$N_fpc[idx]
    )
    id <- if (length(ids) > 0L) {
      frame[idx[1L], ids, drop = FALSE]
    } else {
      data.frame()
    }
    id$.domain <- names(idx_list)[i]
    id$population <- sum(state$population[idx])
    id$n <- sum(n_h[idx])
    id$expected_respondents <- sum(n_h[idx] * state$response_rate[idx])
    id$se <- m$se
    id$moe <- m$moe
    id$rmoe <- m$rmoe
    id$cv <- m$cv
    id$cost <- m$cost
    id
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (all(is.na(out$cost))) {
    out$cost <- NULL
  }
  out
}

#' @rdname summary.svyplan_alloc
#' @export
print.summary.svyplan_bethel <- function(x, ...) {
  .check_unused_dots(...)
  design <- identical(x$kind, "design")
  cat(
    if (design) {
      "Generalized Bethel allocation summary\n\n"
    } else {
      "Generalized Bethel precision summary\n\n"
    }
  )
  cat(sprintf("Question: %s\n", x$question))
  cat(sprintf("Status: %s\n", x$status))
  cat(sprintf("Stages: %d\n\n", x$stages))

  .print_bethel_overall(
    x$overall,
    if (design) "Operational field design" else "Supplied allocation"
  )
  if (!is.null(x$continuous)) {
    cat("\nContinuous optimum\n")
    .print_bethel_overall(x$continuous, title = NULL)
  }

  cat("\nAllocation by stratum\n\n")
  main <- intersect(
    c(
      "stratum",
      "population",
      "n_continuous",
      "n_field",
      "n_supplied",
      "weight",
      "cost"
    ),
    names(x$allocation)
  )
  .print_bethel_table(x$allocation[main])
  stages <- intersect(
    c(
      "stratum",
      "population_psu",
      "population_ssu",
      "psu_continuous",
      "psu_field",
      "psu_supplied",
      "take_per_psu",
      "take_per_ssu"
    ),
    names(x$allocation)
  )
  if (length(stages) > 1L) {
    cat("\nStage decisions by stratum\n\n")
    .print_bethel_table(x$allocation[stages])
  }

  .print_bethel_constraints(
    x$constraints,
    if (design) {
      "Continuous precision constraints"
    } else {
      "Assessed precision constraints"
    }
  )
  if (design) {
    .print_bethel_constraints(
      x$operational_constraints,
      "Operational precision constraints"
    )
  }

  .print_bethel_objective(
    x$objective,
    if (design) {
      "Continuous objective components"
    } else {
      "Assessed objective components"
    }
  )
  if (design) {
    .print_bethel_objective(
      x$operational_objective,
      "Operational objective components"
    )
  }

  cat("\nAllocation bounds\n\n")
  .print_bethel_table(x$bounds)
  if (x$overall$n_bound_violations == 0L) {
    cat("\nNo allocation bounds are violated.\n")
  }

  if (!is.null(x$optimization)) {
    .print_bethel_optimization(x$optimization)
  }

  cat("\nAssumptions\n")
  cat(sprintf("  variance model: %s\n", x$assumptions$variance_model))
  cat(sprintf("  alpha: %s\n", .fmt_prob(x$assumptions$alpha)))
  cat(sprintf(
    "  feasibility tolerance: %.3g\n",
    x$assumptions$feasibility_tolerance
  ))
  if (!is.null(x$assumptions$df) && any(is.finite(x$assumptions$df))) {
    cat(sprintf(
      "  interval df: %s\n",
      paste(
        unique(x$assumptions$df[is.finite(x$assumptions$df)]),
        collapse = ", "
      )
    ))
  }
  if (!is.null(x$assumptions$fixed_takes)) {
    cat("\nFixed later-stage takes\n\n")
    takes <- x$assumptions$fixed_takes
    if (all(is.na(takes$n_per_ssu))) {
      takes$n_per_ssu <- NULL
    }
    .print_bethel_table(takes)
  }
  if (!is.null(x$certainty)) {
    cat("\nCertainty split by stratum\n")
    by_source <- x$certainty_sources[x$certainty_sources > 0L]
    if (length(by_source)) {
      cat(sprintf(
        "certainty PSUs by source: %s\n",
        paste(names(by_source), by_source, collapse = ", ")
      ))
    }
    cat("\n")
    tab <- x$certainty
    tab$threshold <- round(tab$threshold)
    .print_bethel_table(tab)
  }
  if (!is.null(x$assumptions$model) && nrow(x$assumptions$model) > 0L) {
    cat("\nResolved target-stratum planning inputs\n\n")
    .print_bethel_table(x$assumptions$model)
  }
  invisible(x)
}

#' Print generalized allocation headline values
#' @keywords internal
#' @noRd
.print_bethel_overall <- function(values, title) {
  if (!is.null(title)) {
    cat(title, "\n", sep = "")
  }
  rows <- list(
    Sample = .fmt_continuous_n(values$n),
    Cost = sprintf("%.2f", values$cost)
  )
  rows[["Hard targets"]] <- if (values$n_constraints == 0L) {
    "none"
  } else {
    sprintf(
      "%d of %d pass%s",
      values$n_pass,
      values$n_constraints,
      if (values$n_violated == 0L) {
        ""
      } else {
        sprintf(", %d violated", values$n_violated)
      }
    )
  }
  if (!is.null(values$objective_value)) {
    rows[["Weighted rel. variance"]] <- sprintf(
      "%.6g (%d component%s)",
      values$objective_value,
      values$objective_components,
      if (values$objective_components == 1L) "" else "s"
    )
  }
  if (!is.null(values$budget)) {
    rows[["Budget"]] <- sprintf("%.6g", values$budget)
    if (!is.null(values$budget_residual)) {
      rows[["Residual"]] <- sprintf("%.6g", values$budget_residual)
    }
  }
  width <- max(nchar(names(rows)))
  for (name in names(rows)) {
    cat(sprintf("  %-*s  %s\n", width, name, rows[[name]]))
  }
  invisible(NULL)
}

#' Print both the result and numerical diagnostics for every hard target
#' @keywords internal
#' @noRd
.print_bethel_constraints <- function(tab, title) {
  cat("\n", title, "\n\n", sep = "")
  if (is.null(tab) || nrow(tab) == 0L) {
    cat("No hard precision constraints.\n")
    return(invisible(NULL))
  }
  result <- intersect(
    c(
      "constraint",
      ".metric",
      ".target",
      ".achieved",
      ".ratio",
      ".pass",
      ".binding"
    ),
    names(tab)
  )
  .print_bethel_table(tab[result])
  diagnostics <- intersect(
    c("constraint", ".residual", ".tolerance", ".multiplier", ".sensitivity"),
    names(tab)
  )
  if (length(diagnostics) > 1L) {
    cat("\nConstraint diagnostics\n\n")
    .print_bethel_table(tab[diagnostics])
  }
  invisible(NULL)
}

#' Print a genuinely additive generalized-allocation objective
#' @keywords internal
#' @noRd
.print_bethel_objective <- function(tab, title) {
  if (is.null(tab) || nrow(tab) == 0L) {
    return(invisible(NULL))
  }
  cat("\n", title, " (contributions are additive)\n\n", sep = "")
  cols <- intersect(
    c("component", "priority", ".relvar", ".cv", ".contribution", ".share"),
    names(tab)
  )
  .print_bethel_table(tab[cols])
  invisible(NULL)
}

#' Print solver certification and local sensitivity
#' @keywords internal
#' @noRd
.print_bethel_optimization <- function(opt) {
  cat("\nOptimization diagnostics\n")
  cat(sprintf("  classification: %s\n", opt$classification))
  cat(sprintf("  converged: %s\n", if (isTRUE(opt$converged)) "yes" else "no"))
  if (!is.null(opt$iterations)) {
    cat(sprintf("  solver iterations: %d\n", opt$iterations))
  }
  if (length(opt$active_precision) > 0L) {
    cat(sprintf(
      "  active precision: %s\n",
      paste(opt$active_precision, collapse = ", ")
    ))
  }
  cat(sprintf(
    "  active bounds: %d lower, %d upper\n",
    length(opt$active_lower),
    length(opt$active_upper)
  ))
  if (!is.null(opt$certainty)) {
    cc <- opt$certainty
    cat(sprintf("  certainty verdict: %s\n", cc$verdict))
    cat(sprintf(
      "  classification orbit: %d, reached in %d iteration%s\n",
      cc$orbit,
      cc$iterations,
      if (cc$iterations == 1L) "" else "s"
    ))
    cat(sprintf(
      "  allocation settled in %d iteration%s%s\n",
      cc$settle_iterations,
      if (cc$settle_iterations == 1L) "" else "s",
      if (cc$absorbed > 0L) {
        sprintf(
          ", absorbing %d PSU(s) the settle put above a threshold",
          cc$absorbed
        )
      } else {
        ""
      }
    ))
    cat(sprintf(
      "  returned plan is a fixed point: %s\n",
      if (isTRUE(cc$fixed_point)) "yes" else "no"
    ))
  }
  if (!is.null(opt$budget_binding)) {
    cat(sprintf(
      "  budget binding: %s\n",
      if (isTRUE(opt$budget_binding)) "yes" else "no"
    ))
  }
  if (!is.null(opt$budget_sensitivity) && is.finite(opt$budget_sensitivity)) {
    cat(sprintf("  budget sensitivity: %.6g\n", opt$budget_sensitivity))
  }
  if (nrow(opt$diagnostics) > 0L) {
    cat("\nKKT certification\n\n")
    .print_bethel_table(opt$diagnostics)
  }
  invisible(NULL)
}

#' Format and print one bounded-width diagnostic table
#' @keywords internal
#' @noRd
.print_bethel_table <- function(tab, max_rows = 20L) {
  total <- nrow(tab)
  shown <- utils::head(tab, max_rows)
  print(
    .fmt_bethel_summary_table(shown),
    row.names = FALSE,
    quote = FALSE,
    right = TRUE
  )
  if (total > nrow(shown)) {
    cat(sprintf(
      "... %d of %d rows shown (all rows retained in the object)\n",
      nrow(shown),
      total
    ))
  }
  invisible(NULL)
}

#' Format generalized diagnostic tables without changing stored values
#' @keywords internal
#' @noRd
.fmt_bethel_summary_table <- function(tab) {
  out <- tab
  id_cols <- intersect(c("constraint", "component"), names(out))
  for (name in id_cols) {
    out[[name]] <- vapply(out[[name]], .abbrev_bethel_id, character(1L))
  }
  count_cols <- intersect(
    c("population", "population_psu", "population_ssu", "n_field", "field"),
    names(out)
  )
  for (name in count_cols) {
    out[[name]] <- vapply(out[[name]], .fmt_count_n, character(1L))
  }
  two_cols <- intersect(
    c(
      "n_continuous",
      "n_supplied",
      "continuous",
      "supplied",
      "lower",
      "upper",
      "psu_continuous",
      "psu_supplied",
      "cost"
    ),
    names(out)
  )
  for (name in two_cols) {
    out[[name]] <- sprintf("%.2f", out[[name]])
  }
  if ("weight" %in% names(out)) {
    out$weight <- sprintf("%.3f", out$weight)
  }
  number_cols <- intersect(
    c(
      "mean",
      "variance",
      "deff",
      "response_rate",
      "priority",
      ".metric",
      ".target",
      ".achieved",
      ".ratio",
      ".residual",
      ".tolerance",
      ".multiplier",
      ".sensitivity",
      ".relvar",
      ".cv",
      ".contribution",
      ".share",
      "value",
      "tolerance"
    ),
    names(out)
  )
  number_cols <- setdiff(number_cols, ".metric")
  for (name in number_cols) {
    out[[name]] <- ifelse(is.na(out[[name]]), "", sprintf("%.4g", out[[name]]))
  }
  logical_cols <- names(out)[vapply(out, is.logical, logical(1L))]
  for (name in logical_cols) {
    out[[name]] <- ifelse(
      is.na(out[[name]]),
      "",
      ifelse(out[[name]], "yes", "no")
    )
  }
  labels <- c(
    stratum = "Stratum",
    population = "Pop.",
    n_continuous = "Cont.",
    n_field = "Field",
    n_supplied = "Supplied",
    weight = "Weight",
    cost = "Cost",
    population_psu = "PSUs pop.",
    population_ssu = "SSUs pop.",
    psu_continuous = "PSUs cont.",
    psu_field = "PSUs field",
    psu_supplied = "PSUs supplied",
    take_per_psu = "Take/PSU",
    take_per_ssu = "Take/SSU",
    constraint = "Constraint",
    component = "Component",
    .metric = "Metric",
    .target = "Target",
    .achieved = "Ach.",
    .ratio = "Ratio",
    .pass = "Pass",
    .binding = "Bind",
    .residual = "Residual",
    .tolerance = "Tol.",
    .multiplier = "Mult.",
    .sensitivity = "Sens.",
    priority = "Priority",
    .relvar = "Rel. var.",
    .cv = "CV",
    .contribution = "Contribution",
    .share = "Share",
    lower = "Lower",
    continuous = "Cont.",
    field = "Field",
    supplied = "Supplied",
    upper = "Upper",
    continuous_status = "Cont. status",
    field_status = "Field status",
    lower_violation = "Lower bad",
    upper_violation = "Upper bad",
    status = "Status",
    pass = "Pass",
    mean = "Mean",
    variance = "Variance",
    deff = "Deff",
    response_rate = "Resp.",
    n_per_psu = "Take/PSU",
    n_per_ssu = "Take/SSU",
    diagnostic = "Diagnostic",
    value = "Value",
    tolerance = "Tolerance"
  )
  hit <- names(out) %in% names(labels)
  names(out)[hit] <- unname(labels[names(out)[hit]])
  out
}

#' Keep identifier columns from determining console width
#' @keywords internal
#' @noRd
.abbrev_bethel_id <- function(x, width = 20L) {
  x <- as.character(x)
  if (is.na(x) || nchar(x) <= width) {
    return(x)
  }
  paste0(substr(x, 1L, width - 3L), "...")
}

#' @rdname summary.svyplan_alloc
#' @export
print.summary.svyplan_alloc <- function(x, ...) {
  .check_unused_dots(...)
  cluster <- any(grepl("^n_psu", names(x$allocation)))
  cat(
    if (identical(x$kind, "design")) {
      if (cluster) {
        "Stratified cluster allocation summary\n\n"
      } else {
        "Stratified allocation summary\n\n"
      }
    } else {
      if (cluster) {
        "Cluster allocation precision summary\n\n"
      } else {
        "Allocation precision summary\n\n"
      }
    }
  )
  cat(sprintf("Question: %s\n\n", x$question))

  .print_alloc_summary_overall(
    x$overall,
    if (identical(x$kind, "design")) {
      "Overall field design"
    } else {
      "Overall precision"
    },
    show_response = any(x$assumptions$response_rate < 1)
  )
  if (!is.null(x$continuous)) {
    cat("\nContinuous optimum\n")
    .print_alloc_summary_values(
      x$continuous,
      show_response = any(x$assumptions$response_rate < 1),
      indent = "  "
    )
  }

  cat("\nAllocation by stratum\n\n")
  alloc_main <- intersect(
    c(
      "stratum",
      "population",
      "response_rate",
      "n_continuous",
      "n_field",
      "n_supplied",
      "expected_respondents"
    ),
    names(x$allocation)
  )
  print(
    .fmt_alloc_summary_table(x$allocation[alloc_main]),
    row.names = FALSE,
    quote = FALSE,
    right = TRUE
  )

  stage_cols <- intersect(
    c(
      "stratum",
      "n_psu_continuous",
      "n_psu_field",
      "n_per_psu_continuous",
      "n_per_psu_field",
      "n_psu",
      "n_per_psu"
    ),
    names(x$allocation)
  )
  if (length(stage_cols) > 1L) {
    cat("\nCluster stages by stratum\n\n")
    print(
      .fmt_alloc_summary_table(x$allocation[stage_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  fieldwork_cols <- intersect(
    c("stratum", "unit_cost", "weight", "cost"),
    names(x$allocation)
  )
  if (length(fieldwork_cols) > 1L) {
    cat("\nCost and weights by stratum\n\n")
    print(
      .fmt_alloc_summary_table(x$allocation[fieldwork_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  cat("\nAchieved precision by stratum\n\n")
  print(
    .fmt_alloc_summary_table(x$precision),
    row.names = FALSE,
    quote = FALSE,
    right = TRUE
  )

  if (!is.null(x$bounds)) {
    cat("\nAllocation bounds\n\n")
    print(
      .fmt_alloc_summary_table(x$bounds),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
    active <- sum(x$bounds$binding, na.rm = TRUE)
    cat(sprintf(
      "\n%s\n",
      if (active == 0L) {
        "No allocation bounds are active."
      } else {
        sprintf(
          "%d allocation bound%s active.",
          active,
          if (active == 1L) " is" else "s are"
        )
      }
    ))
  }

  if (!is.null(x$domains)) {
    measures <- c(
      ".domain",
      "population",
      "n",
      "expected_respondents",
      "se",
      "moe",
      "rmoe",
      "cv",
      "cost"
    )
    ids <- setdiff(names(x$domains), measures)
    key <- if (length(ids) > 0L) ids else ".domain"
    allocation_cols <- intersect(
      c(key, "population", "n", "expected_respondents", "cost"),
      names(x$domains)
    )
    precision_cols <- intersect(
      c(key, "se", "moe", "rmoe", "cv"),
      names(x$domains)
    )
    cat("\nDomain allocation (not additive to the overall row)\n\n")
    print(
      .fmt_alloc_summary_table(x$domains[allocation_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
    cat("\nDomain precision\n\n")
    print(
      .fmt_alloc_summary_table(x$domains[precision_cols]),
      row.names = FALSE,
      quote = FALSE,
      right = TRUE
    )
  }

  cat("\nAssumptions\n")
  cat(sprintf("  alpha: %s\n", .fmt_prob(x$assumptions$alpha)))
  cat(sprintf(
    "  design effect: %s\n",
    .fmt_by_stratum(x$assumptions$deff, "%.2f")
  ))
  cat(sprintf(
    "  response rate: %s\n",
    .fmt_by_stratum(x$assumptions$response_rate, "%.2f")
  ))
  if (!is.null(x$assumptions$df)) {
    cat(sprintf("  interval df: %g\n", x$assumptions$df))
  }
  invisible(x)
}

#' Print named overall allocation values
#' @keywords internal
#' @noRd
.print_alloc_summary_overall <- function(values, title, show_response) {
  cat(title, "\n", sep = "")
  .print_alloc_summary_values(values, show_response, indent = "  ")
}

#' @keywords internal
#' @noRd
.print_alloc_summary_values <- function(values, show_response, indent) {
  labels <- c(
    n = "Sample",
    expected_respondents = "Expected respondents",
    cost = "Cost",
    se = "SE",
    moe = "MOE",
    rmoe = "Relative MOE",
    cv = "CV",
    worst_domain_cv = "Worst domain CV",
    design_df = "Design df"
  )
  keep <- names(values)
  keep <- keep[keep %in% names(labels)]
  if (!show_response) {
    keep <- setdiff(keep, "expected_respondents")
  }
  keep <- keep[vapply(
    values[keep],
    function(v) {
      length(v) == 1L && !is.null(v) && !is.na(v)
    },
    logical(1L)
  )]
  width <- max(nchar(labels[keep]))
  for (name in keep) {
    value <- values[[name]]
    shown <- switch(
      name,
      n = .fmt_continuous_n(value),
      expected_respondents = sprintf("%.1f", value),
      cost = sprintf("%.1f", value),
      design_df = sprintf("%g", value),
      sprintf("%.4f", value)
    )
    cat(sprintf("%s%-*s  %s\n", indent, width, labels[[name]], shown))
  }
  invisible(NULL)
}

#' Format allocation summary tables without changing stored values
#' @keywords internal
#' @noRd
.fmt_alloc_summary_table <- function(tab) {
  out <- tab
  count_cols <- intersect(
    c("population", "n_field", "n_psu_field", "field"),
    names(out)
  )
  for (name in count_cols) {
    out[[name]] <- vapply(out[[name]], .fmt_count_n, character(1L))
  }
  two_cols <- intersect(
    c(
      "unit_cost",
      "n_continuous",
      "n_supplied",
      "n_psu_continuous",
      "n_per_psu_continuous",
      "n_per_psu",
      "expected_respondents",
      "weight",
      "cost",
      "n"
    ),
    names(out)
  )
  for (name in two_cols) {
    out[[name]] <- ifelse(is.na(out[[name]]), "", sprintf("%.2f", out[[name]]))
  }
  if ("n_per_psu_field" %in% names(out)) {
    out$n_per_psu_field <- vapply(
      out$n_per_psu_field,
      .fmt_count_n,
      character(1L)
    )
  }
  for (name in intersect(c("lower", "upper"), names(out))) {
    out[[name]] <- ifelse(
      is.na(out[[name]]),
      "",
      vapply(out[[name]], .fmt_continuous_n, character(1L))
    )
  }
  if ("response_rate" %in% names(out)) {
    out$response_rate <- sprintf("%.3g", out$response_rate)
  }
  four_cols <- intersect(
    c("effective_n", "se", "moe", "rmoe", "cv"),
    names(out)
  )
  for (name in four_cols) {
    out[[name]] <- ifelse(is.na(out[[name]]), "", sprintf("%.4f", out[[name]]))
  }
  if ("variance_share" %in% names(out)) {
    out$variance_share <- ifelse(
      is.na(out$variance_share),
      "",
      sprintf("%.3f", out$variance_share)
    )
  }
  if ("binding" %in% names(out)) {
    out$binding <- ifelse(is.na(out$binding), "", ifelse(out$binding, "*", ""))
  }
  if ("source" %in% names(out)) {
    out$source[is.na(out$source)] <- ""
  }

  labels <- c(
    stratum = "Stratum",
    population = "Pop.",
    response_rate = "Resp.",
    unit_cost = "Cost/unit",
    n_continuous = "Cont.",
    n_field = "Field",
    n_supplied = "Supplied",
    n_psu_continuous = "PSUs cont.",
    n_psu_field = "PSUs field",
    n_per_psu_continuous = "Take cont.",
    n_per_psu_field = "Take field",
    n_psu = "PSUs",
    n_per_psu = "Take",
    expected_respondents = "Exp. resp.",
    weight = "Weight",
    cost = "Cost",
    effective_n = "Eff. n",
    se = "SE",
    moe = "MOE",
    rmoe = "Rel. MOE",
    cv = "CV",
    variance_share = "Var. share",
    lower = "Lower",
    field = "Field",
    upper = "Upper",
    binding = "Binding",
    source = "Source",
    .domain = "Domain",
    n = "Sample"
  )
  hit <- names(out) %in% names(labels)
  names(out)[hit] <- unname(labels[names(out)[hit]])
  out
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
    ceiling(x$n),
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
  as.integer(ceiling(x$n))
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
  take_all_above <- object$params$take_all_above %||% NULL
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
  cat("Two-phase allocation (", nrow(x$detail), " phase-2 strata)\n", sep = "")
  cat(.fmt_twophase_sizes(x))
  cat(.fmt_twophase_responding(x))
  cat(.fmt_twophase_cost(x))
  cat(.fmt_twophase_deff(x))
  cat("\n")
  print(.fmt_twophase_detail(x, brief = TRUE), row.names = FALSE)
  cat("\n")
  cat(.fmt_twophase_assured(x))
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
    format(round(shown$n[["n_phase1"]] * (x$params$resp_rate %||% 1))),
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
    if (is.na(shown$cv)) "NA" else formatC(shown$cv, format = "f", digits = 4),
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
  out$nu <- formatC(d$nu, format = "f", digits = 4)
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

#' Analyse a two-phase allocation
#'
#' `print()` on a [n_twophase()] result gives the design as it would be
#' fielded: the two phase sizes in whole units, the precision and cost they
#' buy, the per-stratum subsampling fractions, and which of the two designs
#' to run. `summary()` gives what that answer was chosen against: the
#' continuous optimum beside the fielded one, the full stratum table
#' including the continuous issue, the single-phase comparator with its cost
#' and whether it reaches the target, and the planning assumptions.
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

  cat("\nStrata\n")
  print(x$detail, row.names = FALSE)

  s <- x$single_phase
  if (!is.null(s) && is.finite(s$n)) {
    cat(sprintf(
      "\nSingle-phase comparator: n = %s, cv = %s, cost = %s\n",
      format(round(s$n)),
      if (is.na(s$cv)) "NA" else formatC(s$cv, format = "f", digits = 4),
      format(round(s$cost))
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
  cv <- if (is.na(s$cv)) "NA" else formatC(s$cv, format = "f", digits = 4)
  if (isTRUE(s$better)) {
    sprintf(
      "single-phase is better here: n = %s at cv %s, so skip phase 1\n",
      format(round(s$n)),
      cv
    )
  } else {
    sprintf(
      "single-phase alternative: n = %s at cv %s, two-phase wins\n",
      format(round(s$n)),
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
  cat(sprintf(
    "  df = %g   (%g %s - %d %s)\n",
    as.double(x),
    n_units,
    unit_label,
    n_strata,
    constraint_label
  ))
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
  structure(
    list(
      df = as.double(object),
      n_units = attr(object, "n_units", exact = TRUE),
      n_strata = attr(object, "n_strata", exact = TRUE),
      stage = stage,
      strata = attr(object, "strata", exact = TRUE),
      domains = attr(object, "domains", exact = TRUE),
      basis = if (is.null(attr(object, "strata", exact = TRUE)$n_groups)) {
        sprintf("counted %s minus contributing strata", unit_label)
      } else {
        sprintf("counted %s minus variance groups (zones, or collapsed zones)",
                unit_label)
      }
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
    domains = attr(x, "domains", exact = TRUE)
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

#' Print, format and coerce a rotation pattern
#'
#' A `svyplan_rotation` from [design_rotation()] is a numeric vector of takes,
#' one entry per occasion of a unit's life. `as.double()` gives that vector
#' bare, and `as.data.frame()` gives it one row per occasion.
#'
#' @param x A `svyplan_rotation` object.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   `as.data.frame()` arguments.
#' @param ... Additional arguments are not supported and produce an error.
#' @return The method-specific result:
#' \describe{
#'   \item{`print()`}{Returns `x` invisibly.}
#'   \item{`format()`}{Returns a string.}
#'   \item{`as.double()`}{Returns the take vector.}
#'   \item{`as.data.frame()`}{Returns one row per occasion, with `occasion`,
#'     `in_sample` and `take`.}
#' }
#'
#' @seealso [design_rotation()], which builds these objects, and
#'   [design_overlap()] for the lag profile one produces.
#'
#' @examples
#' cps <- design_rotation("4-8-4")
#' cps
#' as.double(cps)
#' as.data.frame(design_rotation("1-1-0-0-1-1"))
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
  ...
) {
  .check_unused_dots(...)
  take <- as.double(x)
  data.frame(
    occasion = seq_along(take),
    in_sample = take > 0,
    take = take,
    row.names = row.names,
    stringsAsFactors = stringsAsFactors
  )
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
#'   `x[["12"]]` are both the twelve-occasion overlap. The result is a bare
#'   number, carrying neither the class nor the lag as a name.
#' @param value Replacement value. Replacement is refused, an overlap being
#'   computed from a rotation rather than assembled.
#' @param name Field to extract: `overlap`, `shared`, `n_occasion`,
#'   `rotation`, `lag` or `life`.
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
#' nothing may move the values while leaving the rest: subsetting and
#' arithmetic return bare numerics, and replacement is an error. That also
#' settles `pmax()` and `pmin()`, which copy the attributes of their first
#' argument without dispatching to any method and would otherwise return
#' something still labelled an overlap whose counts no longer follow from its
#' values. They assign through `[<-`, so they are refused here too. Convert
#' with `as.double(x)` and they work as usual.
#'
#' @seealso [design_overlap()], which builds these objects,
#'   [design_rotation()] for the pattern behind them,
#'   [plot.svyplan_overlap()] for the rotation chart, and
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
  cat("Rotation overlap (planning)\n\n")
  cat(.fmt_rotation_header(rot))
  cat("\n")
  # n_occasion is constant by construction and is already in the header;
  # as.data.frame() keeps it, since a table read into code wants it.
  tab <- as.data.frame(x)[c("lag", "shared", "overlap")]
  shown <- utils::head(tab, 16L)
  shown$overlap <- sprintf("%.4g", shown$overlap)
  print(shown, row.names = FALSE, right = FALSE)
  if (nrow(tab) > nrow(shown)) {
    cat(sprintf(
      "  %d further lag%s, read them with as.data.frame()\n",
      nrow(tab) - nrow(shown),
      if (nrow(tab) - nrow(shown) > 1L) "s" else ""
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

#' The two lines that describe a rotation
#'
#' Shared by the rotation's own print method and the overlap's, so the two
#' cannot drift into describing the same design differently.
#' @keywords internal
#' @noRd
.fmt_rotation_header <- function(rot) {
  paste0(
    sprintf(
      "  %s-occasion life, %s in sample each occasion\n",
      .fmt_count_n(attr(rot, "life", exact = TRUE)),
      .fmt_count_n(attr(rot, "n_occasion", exact = TRUE))
    ),
    sprintf("  rotation: %s\n", .fmt_rotation(as.double(rot)))
  )
}

#' Render a rotation as the spells a planner declared
#'
#' Run-length form, since that is how a rotation is named and argued about.
#' The per-occasion vector it expands to is what the arithmetic reads.
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
    n_occasion = attr(rot, "n_occasion", exact = TRUE),
    rotation = rot,
    lag = as.integer(names(x)),
    life = attr(rot, "life", exact = TRUE)
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
    n_occasion = attr(.overlap_rotation(x), "n_occasion", exact = TRUE),
    overlap = as.double(x),
    row.names = row.names,
    stringsAsFactors = stringsAsFactors
  )
}

#' Print, summarise, format and coerce a panel recruitment
#'
#' Display and coercion methods for the object [n_panel()] and
#' [prec_panel()] return. `print()` gives the number to recruit, the
#' responding sample it is expected to leave, the precision at the target
#' and the wave-by-wave table. `summary()` adds what the design implies
#' around that answer: the standing sample a rotating design holds, the
#' response and retention it assumes, where the life's loss falls, and the
#' launch a design reaching its steady state passes through.
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
#'     plan, full wave table, launch path, and cohort composition.}
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
#' # A rotating design reaching its steady state: the launch table is the
#' # occasions before it gets there
#' rot <- n_panel(
#'   n_prop(p = 0.5, moe = 0.031),
#'   retention = c(0.878, 0.963, 0.936, 0.956),
#'   resp_rate = 0.728,
#'   design = "rotating",
#'   start = "immediate"
#' )
#' summary(rot)$launch
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
  cat(.fmt_panel_headline(x))
  cat(.fmt_panel_precision(x))
  cat(.fmt_panel_assured(x))
  cat("\n")
  print(.fmt_panel_waves(x, brief = TRUE), row.names = FALSE, right = FALSE)
  cat("\n")
  cat("# summary() for the launch, the loss and per-wave cv\n")
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
      launch = .fmt_panel_launch_table(object),
      composition = object$launch_waves
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
  cat(.fmt_panel_headline(plan))
  cat(.fmt_panel_in_sample(plan))
  cat(.fmt_panel_rates(plan))
  cat(.fmt_panel_target_rate(plan))
  cat(.fmt_panel_precision(plan))
  cat(.fmt_panel_assured(plan))

  cat(sprintf(
    "\n%s\n",
    if (identical(plan$design, "fixed")) {
      "Waves of the life"
    } else {
      "Waves alive at one occasion"
    }
  ))
  print(x$waves, row.names = FALSE, right = FALSE)

  if (!is.null(x$launch)) {
    cat(sprintf("\nLaunch (%s)\n", plan$start))
    print(x$launch, row.names = FALSE, right = FALSE)
    cat("# $composition for the cohort mix at each launch occasion\n")
  }
  invisible(x)
}

#' Counts for the whole units a planner would actually release
#'
#' Every count in the printed block is derived from the rounded-up
#' recruitment, so the headline and the wave it names agree. The stored
#' fields stay continuous, which is what keeps the round trip exact.
#' @keywords internal
#' @noRd
.panel_shown <- function(x) {
  recruit <- ceiling(.panel_recruit(x))
  wave <- recruit * x$waves$q
  list(
    recruit = recruit,
    wave = wave,
    head = if (identical(x$design, "fixed")) {
      wave[[x$target_wave]]
    } else {
      sum(wave)
    },
    in_sample = if (identical(x$design, "rotating")) recruit * x$n_cohorts
  )
}

#' Name the launch and the two counts that bracket it
#'
#' The occasion-1 figure against the steady-state one is the whole content:
#' a gradual launch opens below the design and climbs, an immediate one opens
#' at or above it, every unit there being at wave 1 and no later wave holding
#' more than wave 1 does. The two coincide where nothing is lost after
#' recruitment. The table carries the rest.
#' @keywords internal
#' @noRd
.fmt_panel_launch <- function(x) {
  if (is.null(x$launch)) {
    return("")
  }
  # whole units on the display path, as every other count printed here is
  recruit <- ceiling(.panel_recruit(x))
  scale <- recruit / .panel_recruit(x)
  settled <- which(x$launch$steady_state)[[1L]]
  sprintf(
    "launch (%s): %s responding at occasion 1, %s from occasion %d\n",
    x$start,
    .fmt_count_n(round(x$launch$n_resp[[1L]] * scale)),
    .fmt_count_n(round(x$launch$n_resp[[settled]] * scale)),
    x$launch$period[[settled]]
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
  lead <- if (identical(x$design, "fixed")) {
    sprintf("issued: %s", .fmt_count_n(shown$recruit))
  } else {
    sprintf("entrants: %s per occasion", .fmt_count_n(shown$recruit))
  }
  need <- ceiling(x$n_target)
  sprintf(
    "%s -> %s responding%s%s\n",
    lead,
    .fmt_count_n(round(shown$head)),
    where,
    if (round(shown$head) < need) {
      sprintf(", short of the %s the target needs", .fmt_count_n(need))
    } else {
      ""
    }
  )
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
  sprintf(
    "in sample: %s across %d live cohorts\n",
    .fmt_count_n(.panel_shown(x)$in_sample),
    x$n_cohorts
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
    .fmt_count_n(ceiling(x$n_assured)),
    if (identical(x$design, "fixed")) "issued" else "entrants per occasion",
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
  sprintf(
    "rates: response %.3g, retention %s (%s)\n",
    x$params$resp_rate,
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
  rr <- .panel_rate_conflict(x$target$params$resp_rate, x$params$resp_rate)
  if (is.null(rr)) {
    return("")
  }
  sprintf(
    "target: response %.3g removed, requirement %s responding\n",
    rr,
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
  sprintf(
    "%s: se = %.4g, moe = %.4g%s\n",
    label,
    x$se,
    x$moe,
    if (is.na(x$cv)) "" else sprintf(", cv = %.3g", x$cv)
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
  w <- x$waves
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
    scale <- .panel_shown(x)$recruit / .panel_recruit(x)
    out$cases <- format(round(w$expected_cases * scale))
  }
  out
}

#' The launch path as a table, on the display path
#'
#' Every count is scaled by the same whole-unit factor the wave table uses,
#' so the occasion that settles reports the sample the headline promises.
#' @keywords internal
#' @noRd
.fmt_panel_launch_table <- function(x) {
  if (is.null(x$launch)) {
    return(NULL)
  }
  scale <- .panel_shown(x)$recruit / .panel_recruit(x)
  l <- x$launch
  out <- data.frame(
    occasion = l$period,
    entrants = format(round(l$n_entrants * scale)),
    in_sample = format(round(l$n_in_sample * scale)),
    n_resp = format(round(l$n_resp * scale)),
    se = sprintf("%.4g", l$se),
    moe = sprintf("%.4g", l$moe),
    stringsAsFactors = FALSE
  )
  if (!all(is.na(l$cv))) {
    out$cv <- sprintf("%.3g", l$cv)
  }
  out$steady <- l$steady_state
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
  as.integer(ceiling(.panel_recruit(x)))
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
