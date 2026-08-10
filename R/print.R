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
#'   `as.data.frame()` returns a data frame; the shapes are described under
#'   Details.
#'
#' @details
#' ## Print and coercion
#'
#' Constrained designs (`n_cluster()`, `n_alloc()`, `n_multi_cluster()`)
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
#' the stable handoff to downstream packages (e.g. samplyr). For
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
  rates <- c(resp_rate_psu = p$resp_rate_psu %||% 1,
             resp_rate_ssu = p$resp_rate_ssu %||% 1,
             resp_rate = p$resp_rate %||% 1)
  set <- rates[rates < 1]
  if (length(set) == 0L) {
    return(invisible(NULL))
  }
  cat(sprintf("(%s)\n", paste(sprintf("%s = %.2f", names(set), set),
                              collapse = ", ")))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.fmt_deff <- function(deff) {
  if (is.null(deff)) return(NULL)
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
#' A single value prints as itself; values that differ print as their range,
#' so the header stays one line however many strata there are.
#' @keywords internal
#' @noRd
.fmt_by_stratum <- function(x, fmt) {
  u <- unique(x)
  if (length(u) == 1L) return(sprintf(fmt, u))
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
  cat(sprintf("n = %d per occasion", ceiling(x$n)))
  if (!is.null(p$resp_rate) && p$resp_rate < 1) {
    cat(sprintf(" (net: %d)", ceiling(x$n * p$resp_rate)))
  }
  cat(sprintf(", %d occasions", p$occasions))
  parts <- c(.fmt_pooled_estimand(p), .fmt_deff(p$deff))
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_lag_overlap(p))
  cat(sprintf("se = %.4g, moe = %.4g", x$se, x$moe))
  if (!is.na(x$cv)) cat(sprintf(", cv = %.4g", x$cv))
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
    ov[1L], rho[1L], reach
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
    out <- sprintf("n = %d per occasion", ceiling(n))
    if (net) {
      out <- sprintf("%s (net: %d)", out, ceiling(n * resp_rate))
    }
    return(out)
  }
  out <- sprintf("n = %d then %d", ceiling(n[1L]), ceiling(n[2L]))
  if (net) {
    out <- sprintf(
      "%s (net: %d then %d)", out,
      ceiling(n[1L] * resp_rate), ceiling(n[2L] * resp_rate)
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
    overlap, overlap_cor, 100 * .change_var_share(p)
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
  if (base <= 0) return(NA_real_)
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
    cat(sprintf("n = %d (net: %d)", ceiling(x$n), net_n))
  } else {
    cat(sprintf("n = %d", ceiling(x$n)))
  }

  parts <- character(0L)
  if (!is.null(p$p)) {
    parts <- c(parts, sprintf("p = %s", .fmt_prob(p$p)))
  }
  if (!is.null(p$var)) {
    parts <- c(parts, sprintf("var = %.2f", p$var))
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
      "n = %d (%s)\n",
      ceiling(x$n),
      if (natural) "largest quota / share" else "sum of domain quotas"
    ))
    if (!is.null(x$n_domain_max)) {
      cat(sprintf(
        "Largest single domain = %d (binding: %s)\n",
        ceiling(x$n_domain_max), x$binding
      ))
    }
    cat("---\n")
    dom <- x$domains
    dom$.n <- ceiling(dom$.n)
    print(dom, row.names = FALSE, right = FALSE)
  } else {
    cat(sprintf("Multi-indicator sample size\n"))
    cat(sprintf("n = %d (binding: %s)\n", ceiling(x$n), x$binding))
    cat("---\n")
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
    (x$params$resp_rate_ssu %||% 1) * (x$params$resp_rate %||% 1)
  if (resp_rate < 1) {
    net_total <- ceiling(total_display * resp_rate)
    cat(sprintf(
      " -> total n = %s (net: %s)\n",
      .fmt_count_n(total_display), .fmt_count_n(net_total)
    ))
  } else {
    cat(sprintf(" -> total n = %s\n", .fmt_count_n(total_display)))
  }
  .print_stage_rates(x$params)
  fc <- x$params$fixed_cost
  op_cv <- if (!is.null(op)) op$cv else x$cv
  op_cost <- if (!is.null(op)) op$cost else x$cost
  if (!is.null(fc) && fc > 0) {
    cat(sprintf("cv = %.4f, cost = %.0f (fixed: %.0f)\n",
                op_cv, op_cost, fc))
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
    paste(cont_parts, collapse = " | "), x$cv, x$cost
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
    cat("---\n")
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
        "Total n = %d (unrounded: %s), cost = %.0f (fixed: %.0f)\n",
        sum(dom$.total_n),
        unrounded,
        x$cost,
        fc
      ))
    } else {
      cat(sprintf(
        "Total n = %d (unrounded: %s)\n",
        sum(dom$.total_n), unrounded
      ))
    }
    dom$.cv <- sprintf("%.4f", dom$.cv)
    dom$.cost <- sprintf("%.0f", dom$.cost)
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
        op_cv, op_cost, fc, x$binding
      ))
    } else {
      cat(sprintf(
        "worst cv = %.4f, cost = %.0f (binding: %s)\n",
        op_cv, op_cost, x$binding
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
      paste(cont_parts, collapse = " | "), x$cv, x$cost
    ))
    cat("---\n")
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
    "Sampling precision for pooled estimate (%s)\n", .pooled_scale_label(p)
  ))
  cat(sprintf("n = %d per occasion", ceiling(p$n)))
  if (!is.null(p$resp_rate) && p$resp_rate < 1) {
    cat(sprintf(" (net: %d)", ceiling(p$n * p$resp_rate)))
  }
  cat(sprintf(", %d occasions", p$occasions))
  parts <- c(.fmt_pooled_estimand(p), .fmt_deff(p$deff))
  cat(sprintf(" (%s)\n", paste(parts, collapse = ", ")))
  cat(.fmt_lag_overlap(p))
  cat(sprintf("se = %.4g, moe = %.4g", x$se, x$moe))
  if (!is.na(x$cv)) cat(sprintf(", cv = %.4g", x$cv))
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
  if (!is.na(x$cv)) cat(sprintf(", cv = %.4g", x$cv))
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
      cat(sprintf(" (net: %s)", .fmt_count_n(ceiling(total_display * resp_rate))))
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
    # An allocation carries one n per stratum; the header reports the total.
    n_display <- sum(ceiling(p$n))
    cat(sprintf("n = %d", n_display))
    if (length(p$n) > 1L) cat(sprintf(" (%d strata)", length(p$n)))
    resp_rate <- p$resp_rate
    if (!is.null(resp_rate) && all(resp_rate < 1)) {
      cat(sprintf(" (net: %d)", ceiling(sum(p$n * resp_rate))))
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
      x$stages, nrow(x$strata)
    ))
    tab <- x$strata
    num <- vapply(tab, is.numeric, logical(1L))
    tab[num] <- lapply(tab[num], function(v) sprintf("%.4f", v))
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
  cat(sprintf("var_ratio = %s\n", paste(sprintf("%.4f", x$var_ratio), collapse = ", ")))
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
      sprintf(" (size-weighted; nominal %.4g)", nominal)
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
        "%s = %d, %s = %d (total = %d, net: %d), power = %.3f, effect = %.4f\n",
        lab1, n1, lab2, n2, n1 + n2, ceiling((n1 + n2) * resp_rate),
        x$power, x$effect
      ))
    } else {
      cat(sprintf(
        "%s = %d, %s = %d (total = %d), power = %.3f, effect = %.4f\n",
        lab1, n1, lab2, n2, n1 + n2, x$power, x$effect
      ))
    }
  } else if (!is.null(resp_rate) && resp_rate < 1) {
    net_n <- ceiling(x$n * resp_rate)
    cat(sprintf(
      "n = %d (net: %d, per group), power = %.3f, effect = %.4f\n",
      ceiling(x$n), net_n, x$power, x$effect
    ))
  } else {
    cat(sprintf(
      "n = %d (per group), power = %.3f, effect = %.4f\n",
      ceiling(x$n), x$power, x$effect
    ))
  }

  parts <- character(0L)
  if (!is.null(p$treat)) {
    parts <- c(parts, sprintf(
      "treat = (%.3f, %.3f)", p$treat[1], p$treat[2]
    ))
  }
  if (!is.null(p$control)) {
    parts <- c(parts, sprintf(
      "control = (%.3f, %.3f)", p$control[1], p$control[2]
    ))
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
    parts <- c(parts, sprintf(
      "var = (%.2f, %.2f, %.2f, %.2f)", p$var[1], p$var[2], p$var[3], p$var[4]
    ))
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
    "svyplan_varcomp [", x$stages, "-stage",
    if (identical(x$source, "deff")) ", from deff" else "", "]"
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
  # Every method reads one variance through the effective size
  # n_eff = n_net / (deff * fpc), exactly as .prec_engine_prop() does, so a
  # confidence interval and the margin of error reported for the same design
  # always agree.
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
#' @param object A [n_prop()], [n_mean()], [prec_prop()], or [prec_mean()]
#'   result.
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
#' computed with (`"wald"`, `"wilson"`, `"logodds"`, or `"beta"`), including
#' its `df` when the beta method carries one. Only the Wald interval is
#' symmetric about `p`, so for the other three the limits are not
#' `p` plus or minus the reported `moe`; `$moe` remains half the interval
#' width, and `confint()` is the way to read where the interval actually
#' sits. All four apply `deff`, `resp_rate`, and the finite population
#' correction through the same effective size the sizing functions use.
#'
#' For means, a symmetric z-interval is used, which requires `mu` in the
#' original call.
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
#' @name confint.svyplan
NULL

#' @rdname confint.svyplan
#' @export
confint.svyplan_n <- function(object, parm, level = 0.95, ...) {
  .check_unused_dots(...)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
    level <= 0 || level >= 1) {
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
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
    level <= 0 || level >= 1) {
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
  if (!is.null(x$operational)) {
    return(as.integer(x$operational$n))
  }
  as.integer(ceiling(x$n))
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
#' [design_effect()] returns. `print()` itemizes the components and the
#' overall value; the coercion and arithmetic methods let the object stand
#' in for that overall value wherever a plain number is expected.
#'
#' @param x A `svyplan_deff` object from [design_effect()].
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
#' @return `print()` returns `x` invisibly, `format()` returns a character
#'   scalar, and `as.double()` returns the overall design effect.
#'   `as.data.frame()` returns a one-row table with the overall value and one
#'   column per component; `as.list()` returns the same fields as a named
#'   list, and `$` and `[[` return one of them. Arithmetic and mathematical
#'   transformations return ordinary numeric results.
#'
#' @details
#' A `svyplan_deff` behaves as the numeric overall design effect wherever one
#' is expected: it can be passed to any `deff` argument, compared, and
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
  components <- attr(x, "components", exact = TRUE)
  notes <- attr(x, "notes", exact = TRUE)
  labels <- c(cluster = "clustering", weight = "weighting",
              strata = "stratification", allocation = "allocation")
  shown <- unname(labels[names(components)])
  width <- max(nchar(c(shown, "overall")))
  cat("Design effect (planning)\n\n")
  for (i in seq_along(components)) {
    cat(sprintf(
      "  %-*s  %7.4f   %s\n", width, shown[i], components[[i]], notes[[i]]
    ))
  }
  cat(sprintf("  %s\n", strrep("-", width + 11L)))
  # one component is the result itself; multiplying several is Kish's
  # composition, which is exact only when the stratum sd agree. The
  # components carry no sd, so the marker is conservative by design.
  overall_note <- if (length(components) > 1L) "   approx. (Kish)" else ""
  cat(sprintf("  %-*s  %7.4f%s\n", width, "overall", as.double(x), overall_note))
  invisible(x)
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
        paste(name, collapse = ", "), paste(names(fields), collapse = ", ")
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
    as.list(x), row.names = row.names, optional = optional,
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
    kozak = "Kozak-inspired local search",
    x$method
  )
  cat(sprintf("Strata boundaries (%s, %d strata)\n", method_label, x$n_strata))
  cat(sprintf(
    "Boundaries: %s\n",
    paste(sprintf("%.1f", x$boundaries), collapse = ", ")
  ))
  cat(sprintf("n = %d, cv = %.4f\n", ceiling(x$n), x$cv))
  if (!is.null(x$alloc) && is.character(x$alloc)) {
    if (x$alloc == "power") {
      cat(sprintf("Allocation: power (alloc_q = %.2f)\n", x$params$alloc_q))
    } else {
      cat(sprintf("Allocation: %s\n", x$alloc))
    }
  }
  if (isTRUE(x$converged)) {
    cat("Converged: yes\n")
  } else if (isFALSE(x$converged)) {
    cat("Converged: no\n")
  }
  cat("---\n")
  df <- x$strata
  df$share <- sprintf("%.3f", df$share)
  df$sd <- sprintf("%.1f", df$sd)
  df$take_all <- NULL
  print(df, row.names = FALSE, right = FALSE)
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
  if (cluster) cat(", two-stage")
  if (!is.na(H)) cat(sprintf(", %d strata", H))
  cat(")\n")
  op <- x$operational
  if (!is.null(op)) {
    cat(sprintf("field design: n = %d", op$n))
    if (cluster) cat(sprintf(", n_psu = %d", sum(detail$n_psu_int)))
    if (!is.na(op$cv)) cat(sprintf(", cv = %.4f", op$cv))
    if (!is.na(op$cost)) cat(sprintf(", cost = %.0f", op$cost))
    cat("\n")
    cat(sprintf("continuous optimum: n = %s", .fmt_continuous_n(x$n)))
    if (!is.na(x$cv)) cat(sprintf(", cv = %.4f", x$cv))
    if (!is.na(x$se)) cat(sprintf(", se = %.4f", x$se))
    cat("\n")
  } else {
    cat(sprintf("n = %d", ceiling(x$n)))
    if (cluster) cat(sprintf(", n_psu = %d", sum(detail$n_psu_int)))
    if (!is.na(x$cv)) cat(sprintf(", cv = %.4f", x$cv))
    if (!is.na(x$se)) cat(sprintf(", se = %.4f", x$se))
    cat("\n")
  }
  parts <- character(0L)
  if (!is.null(p$min_n_stratum) && p$min_n_stratum > 0)
    parts <- c(parts, sprintf("min_n_stratum = %g", p$min_n_stratum))
  resp_rate <- p$resp_rate
  if (!is.null(resp_rate) && any(resp_rate < 1))
    parts <- c(parts, sprintf("resp_rate = %s", .fmt_by_stratum(resp_rate, "%.2f")))
  parts <- c(parts, .fmt_deff(p$deff))
  if (length(parts) > 0L)
    cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))
  .print_psu_fraction_note(detail)
  .print_design_df(x)
  .print_domains_block(x$domains)
}

#' Disclose an appreciable first-stage sampling fraction
#'
#' `N_psu` bounds the allocation but activates no first-stage correction, so a
#' design taking a large share of the available PSUs is planned conservatively:
#' the between-PSU term keeps its full with-replacement size. That is a
#' property of the result worth seeing, but not an event to act on, so it is
#' disclosed here rather than warned about. `predict()` re-runs the allocation
#' once per grid row, and a warning would repeat with it.
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
    "note: samples %s of the available PSUs in %s%s; precision uses a with-replacement first stage and may be conservative\n",
    paste0(format(round(100 * detail$.psu_frac[hit]), trim = TRUE), "%",
           collapse = ", "),
    paste(detail$stratum[hit], collapse = ", "),
    if (any(bound)) sprintf(" (bound active in %s)",
                            paste(detail$stratum[bound], collapse = ", ")) else ""
  ))
  invisible(NULL)
}

#' Print the per-domain precision table an allocation carries
#' @keywords internal
#' @noRd
.print_domains_block <- function(dom) {
  if (is.null(dom)) return(invisible(NULL))
  cat(sprintf("Domains: %d\n", nrow(dom)))
  cat("---\n")
  if (".cv" %in% names(dom)) dom$.cv <- sprintf("%.4f", dom$.cv)
  if (".cost" %in% names(dom)) dom$.cost <- sprintf("%.0f", dom$.cost)
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
  if (nrow(sel) == 0L) return(invisible(NULL))
  shown <- utils::head(
    sel[c("constraint", ".metric", ".target", ".achieved", ".pass")],
    6L
  )
  if (nrow(sel) < total) {
    cat(sprintf("showing %d %s of %d\n", nrow(shown), label, total))
  }
  print(shown, row.names = FALSE, right = FALSE)
  if (nrow(shown) < total) cat(hint, "\n", sep = "")
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
  cat("Joint constrained allocation (Bethel)\n")
  # Report exceptions, not confirmations: a line that only ever says
  # "nothing went wrong" buries the lines that do carry news.
  if (!identical(status, "optimal")) {
    cat(sprintf("status: %s\n", status))
  }
  # The two modes return the same class and the same numbers mean different
  # things in each, so the question stays even though it never varies.
  cat(if (budget_mode) {
    sprintf("question: best design affordable within a budget of %.6g\n",
            x$params$budget)
  } else {
    "question: cheapest design meeting every precision target\n"
  })
  n_targets <- nrow(op$constraints)
  cat(sprintf(
    "field design: n = %d, cost = %.0f%s\n",
    op$n, op$cost,
    if (n_targets == 0L) {
      ""
    } else {
      sprintf(" (%d target%s, %s)", n_targets,
              if (n_targets == 1L) "" else "s",
              if (isTRUE(op$all_pass)) "all pass" else "violations")
    }
  ))
  continuous_cost <- x$params$achieved$cost
  increase <- if (continuous_cost > 0) {
    100 * (op$cost / continuous_cost - 1)
  } else {
    NA_real_
  }
  cat(sprintf(
    "continuous optimum: n = %s, cost = %.0f%s\n",
    .fmt_continuous_n(x$n), continuous_cost,
    if (!is.na(increase) && abs(increase) >= 0.005) {
      sprintf(" (integerizing costs %+.2f%%)", increase)
    } else {
      ""
    }
  ))
  if (budget_mode) {
    cat(sprintf(
      "objective: weighted relative variance %.6g continuous, %.6g operational\n",
      x$objective_value, op$objective_value
    ))
    if (!isTRUE(opt$budget_binding)) {
      cat("the budget is not binding: the allocation sits at its upper bounds\n")
    } else if (is.finite(opt$budget_sensitivity %||% NA_real_)) {
      cat(sprintf("one more unit of budget changes the objective by %.4g\n",
                  opt$budget_sensitivity))
    }
    obj <- x$objective
    if (!is.null(obj) && nrow(obj) > 0L) {
      shown <- utils::head(
        obj[c("component", "priority", ".cv", ".share")], 6L
      )
      print(shown, row.names = FALSE, right = FALSE)
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
      sel$constraint[1L], sel$.target[1L], sel$.achieved[1L]
    ))
  }
  n_lower <- length(opt$active_lower %||% integer(0))
  n_upper <- length(opt$active_upper %||% integer(0))
  if (n_lower > 0L || n_upper > 0L) {
    cat(sprintf("active bounds: %d lower, %d upper\n", n_lower, n_upper))
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
    if (all(d$.pass)) "all pass" else
      sprintf("%d violated", sum(!d$.pass))
  ))
  violated <- any(!d$.pass)
  sel <- if (violated) d[!d$.pass, , drop = FALSE] else
    d[d$.binding, , drop = FALSE]
  .print_constraint_rows(
    sel,
    total = nrow(d),
    label = if (violated) "violated" else "binding",
    hint = "... see $detail for all rows"
  )
  if (!is.null(x$objective_value)) {
    cat(sprintf("objective: weighted relative variance %.6g\n",
                x$objective_value))
    if (!is.null(x$params$budget)) {
      cat(sprintf("budget: %.6g, residual %.6g\n",
                  x$params$budget, x$params$budget_residual))
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
  p <- x$params
  cat("Two-phase allocation (", nrow(x$detail), " phase-2 strata)\n", sep = "")
  cat(sprintf(
    "issued: n_phase1 = %s | n_phase2 = %s\n",
    format(round(x$n[["n_phase1"]])), format(round(x$n[["n_phase2"]]))
  ))
  if (!isTRUE(all.equal(unname(x$responding), unname(x$n)))) {
    cat(sprintf(
      "expected responding: n_phase1 = %s | n_phase2 = %s\n",
      format(round(x$responding[["n_phase1"]])),
      format(round(x$responding[["n_phase2"]]))
    ))
  }
  cat(sprintf(
    "cv = %s, cost = %s%s\n",
    if (is.na(x$cv)) "NA" else formatC(x$cv, format = "f", digits = 4),
    format(round(x$cost)),
    if (isTRUE(p$fixed_cost > 0)) sprintf(" (fixed: %s)", format(p$fixed_cost)) else ""
  ))
  if (!isTRUE(all.equal(p$phase1_deff, 1)) ||
      !isTRUE(all.equal(p$single_deff, 1))) {
    cat(sprintf("phase-1 deff = %s, single-phase deff = %s\n",
                formatC(p$phase1_deff, format = "f", digits = 2),
                formatC(p$single_deff, format = "f", digits = 2)))
  }
  d <- x$detail
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
  out$n_issued <- format(round(d$n_issued))
  out$n_int <- format(d$n_int)
  if (show_resp) {
    out$n_resp <- format(round(d$n_resp))
  }
  if (any(d$take_all)) {
    out$take_all <- ifelse(d$take_all, "*", "")
  }
  cat("---\n")
  print(out, row.names = FALSE)
  o <- x$operational
  if (!is.null(o)) {
    cat(sprintf("field design: n_phase1 = %s | n_phase2 = %s (cost %s, cv %s)\n",
                format(o$n[["n_phase1"]]), format(o$n[["n_phase2"]]),
                format(round(o$cost)),
                if (is.na(o$cv)) "NA" else formatC(o$cv, format = "f", digits = 4)))
    if (!is.null(o$assured)) {
      cat(sprintf(
        "assured (%s): issue n_phase1 = %s | n_phase2 = %s (cost %s)\n",
        .fmt_prob(x$params$assurance),
        format(o$assured_phase1), format(sum(o$assured)),
        format(round(o$assured_cost))))
    }
  }
  s <- x$single_phase
  if (isTRUE(s$better)) {
    cat(sprintf(
      "\nSingle-phase is better here: n = %s, cv = %s, cost = %s\n",
      format(round(s$n)),
      if (is.na(s$cv)) "NA" else formatC(s$cv, format = "f", digits = 4),
      format(round(s$cost))
    ))
    cat("Skip phase 1 and measure directly.\n")
  } else if (is.finite(s$n)) {
    cat(sprintf(
      "\nSingle-phase alternative: n = %s, cv = %s, cost = %s (two-phase wins)\n",
      format(round(s$n)),
      if (is.na(s$cv)) "NA" else formatC(s$cv, format = "f", digits = 4),
      format(round(s$cost))
    ))
  }
  invisible(x)
}

#' Print and coerce design degrees of freedom
#'
#' Display and coercion methods for the count that [design_df()] returns.
#' `print()` shows the count and the numbers it was formed from; the
#' coercion and arithmetic methods let the object stand in for that count
#' wherever a plain number is expected.
#'
#' @param x A `svyplan_df` object from [design_df()].
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
#' @return `print()` returns `x` invisibly, `format()` returns a character
#'   scalar, and `as.double()` returns the degrees of freedom.
#'   `as.data.frame()` returns a one-row table of the scalar fields;
#'   `as.list()` returns every field, including the per-stratum and
#'   per-domain tables, and `$` and `[[` return one of them. Arithmetic and
#'   mathematical transformations return ordinary numeric results.
#'
#' @details
#' A `svyplan_df` behaves as the numeric degrees of freedom wherever one is
#' expected: it can be passed to any `df` argument, compared, and
#' arithmetically combined, with the detail dropped by any such operation.
#'
#' The fields are `df` for the count itself, `n_units` for the units it
#' counts and `stage` for what those units are (`"psu"` or `"element"`),
#' `n_strata` for the constraints subtracted, and the `strata` and
#' `domains` tables, which are `NULL` for an unstratified or domain-free
#' plan. Naming a field the object does not have is an error listing the
#' ones it does.
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
  cat("Design degrees of freedom (planning)\n\n")
  cat(sprintf("  df = %g   (%g %s - %d strat%s)\n",
              as.double(x), n_units, unit_label, n_strata,
              if (n_strata == 1L) "um" else "a"))
  strata <- attr(x, "strata", exact = TRUE)
  if (!is.null(strata)) {
    flagged <- strata$.status != "ok"
    if (any(flagged)) {
      cat(sprintf("  no df from %s: %s\n",
                  if (sum(flagged) > 1L) "these strata" else "this stratum",
                  paste(sprintf("%s (%s)", strata$stratum[flagged],
                                strata$.status[flagged]), collapse = ", ")))
    }
  }
  domains <- attr(x, "domains", exact = TRUE)
  if (!is.null(domains)) {
    cat(sprintf("  %d domains, df from %g to %g\n", nrow(domains),
                min(domains$.df), max(domains$.df)))
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
        paste(name, collapse = ", "), paste(names(fields), collapse = ", ")
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

#' Print, format and coerce a rotation overlap
#'
#' A `svyplan_overlap` from [design_overlap()] is a numeric vector of
#' overlap fractions indexed by lag, so `x[1]` and `x[12]` are plain
#' numbers ready for an `overlap` argument. `$` reaches the counts and the
#' schedule behind them.
#'
#' @param x A `svyplan_overlap` object.
#' @param e1,e2 Operands. Arithmetic and comparison return bare numerics,
#'   the counts and the schedule describing the profile as computed and not
#'   whatever it was transformed into.
#' @param i Lag to extract, by position or by its name, so `x[12]` and
#'   `x[["12"]]` are both the twelve-occasion overlap. The result is a bare
#'   number, carrying neither the class nor the lag as a name.
#' @param value Replacement value. Replacement is refused, an overlap being
#'   computed from a schedule rather than assembled.
#' @param name Field to extract: `overlap`, `shared`, `n_occasion`,
#'   `schedule`, `lag` or `life`.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   `as.data.frame()` arguments.
#' @param ... Additional arguments are not supported and produce an error.
#' @return `print()` returns `x` invisibly; `format()` a string;
#'   `as.double()` the overlap vector; `as.data.frame()` one row per lag;
#'   `Ops()` and `Math()` bare numerics. Replacement is an error.
#'
#' @details
#' The values, the shared counts and the schedule describe one design, so
#' nothing may move the values while leaving the rest: subsetting and
#' arithmetic return bare numerics, and replacement is an error. That also
#' settles `pmax()` and `pmin()`, which copy the attributes of their first
#' argument without dispatching to any method and would otherwise return
#' something still labelled an overlap whose counts no longer follow from its
#' values. They assign through `[<-`, so they are refused here too. Convert
#' with `as.double(x)` and they work as usual.
#' @name print.svyplan_overlap
NULL

#' @rdname print.svyplan_overlap
#' @export
print.svyplan_overlap <- function(x, ...) {
  .check_unused_dots(...)
  w <- attr(x, "schedule", exact = TRUE)
  cat("Rotation overlap (planning)\n\n")
  cat(sprintf(
    "  %s-occasion life, %s in sample each occasion\n",
    .fmt_count_n(attr(x, "life", exact = TRUE)),
    .fmt_count_n(attr(x, "n_occasion", exact = TRUE))
  ))
  cat(sprintf("  schedule: %s\n\n", .fmt_schedule(w)))
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

#' Render a schedule as the spells a planner declared
#'
#' Run-length form, since that is how a rotation is named and argued about;
#' the per-occasion vector it expands to is what the arithmetic reads.
#' @keywords internal
#' @noRd
.fmt_schedule <- function(w) {
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
    attr(x, "life", exact = TRUE), unclass(x)[[1L]]
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
  list(
    overlap = as.double(x),
    shared = attr(x, "shared", exact = TRUE),
    n_occasion = attr(x, "n_occasion", exact = TRUE),
    schedule = attr(x, "schedule", exact = TRUE),
    lag = as.integer(names(x)),
    life = attr(x, "life", exact = TRUE)
  )
}

#' @rdname print.svyplan_overlap
#' @export
`[.svyplan_overlap` <- function(x, i) {
  # A lag is picked out to be passed as `overlap`, so it has to come back
  # bare. The default would carry the lag along as a name, and a name on a
  # numeric survives arithmetic: it would reappear on the `$effect` or
  # `$moe` of whatever the number was handed to. The counts and the schedule
  # go too, and for a stronger reason: a subset is no longer the overlap
  # profile they describe, so an `x[]` still carrying them would report a
  # schedule against values that no longer come from it.
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
  # Arithmetic returns bare numerics, as it does for the other classed
  # numerics in the package. Keeping the class would leave a vector whose
  # values had moved while `shared`, `n_occasion` and `schedule` had not, so
  # a doubled overlap would still claim to have come from the schedule that
  # produced the original, and could sit above 1.
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
        paste(name, collapse = ", "), paste(names(fields), collapse = ", ")
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
    n_occasion = attr(x, "n_occasion", exact = TRUE),
    overlap = as.double(x),
    row.names = row.names,
    stringsAsFactors = stringsAsFactors
  )
}

#' Print, format and coerce a panel recruitment
#'
#' Display and coercion methods for the object [n_panel()] and
#' [prec_panel()] return. `print()` leads with the number to recruit and the
#' responding sample it is expected to leave, then the wave-by-wave table.
#' The coercions return the recruitment count, which is a number of units to
#' release and not the analysis sample: those differ by the whole of the
#' panel's attrition, and it is why `svyplan_panel` is a sibling of
#' `svyplan_n` rather than a subtype.
#'
#' @param x A `svyplan_panel` object.
#' @param row.names,optional,stringsAsFactors,validRN Standard
#'   `as.data.frame()` arguments.
#' @param ... Additional arguments are not supported and produce an error.
#' @return `print()` returns `x` invisibly; `format()` a string;
#'   `as.double()` the recruitment count and `as.integer()` the whole units
#'   that count rounds up to; `as.data.frame()` the wave table.
#' @name print.svyplan_panel
NULL

#' @rdname print.svyplan_panel
#' @export
print.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  k <- nrow(x$waves)
  cat(sprintf("Panel recruitment (%s, %d-wave life)\n", x$design, k))
  cat(.fmt_panel_headline(x))
  cat(.fmt_panel_launch(x))
  cat(.fmt_panel_rates(x))
  cat(.fmt_panel_precision(x))
  if (!is.null(x$n_assured)) {
    cat(sprintf(
      "assured (%s): %s %s%s\n",
      .fmt_prob(x$params$assurance), .fmt_count_n(ceiling(x$n_assured)),
      if (identical(x$design, "fixed")) "issued" else "entrants per occasion",
      # A level a finite frame cannot supply is named on the line that
      # reports it, the number being a requirement rather than a design.
      if (isFALSE(x$assured_feasible)) {
        sprintf(", beyond the population of %s", .fmt_count_n(x$target$params$N))
      } else {
        ""
      }
    ))
  }
  cat(if (identical(x$design, "fixed")) {
    "---\n"
  } else {
    "--- (the cohorts alive at one occasion)\n"
  })
  print(.fmt_panel_waves(x), row.names = FALSE, right = FALSE)
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
#' target, which the line says outright.
#' @keywords internal
#' @noRd
.fmt_panel_headline <- function(x) {
  shown <- .panel_shown(x)
  where <- if (identical(x$design, "fixed")) {
    sprintf("at wave %d", x$target_wave)
  } else {
    sprintf("per occasion, pooled over %d cohorts", x$n_cohorts)
  }
  lead <- if (identical(x$design, "fixed")) {
    sprintf("issue %s to hold", .fmt_count_n(shown$recruit))
  } else {
    sprintf("%s entrants per occasion to hold", .fmt_count_n(shown$recruit))
  }
  need <- ceiling(x$n_target)
  out <- sprintf(
    "%s %s responding %s%s\n", lead, .fmt_count_n(round(shown$head)), where,
    if (round(shown$head) < need) {
      sprintf(", short of the %s the target needs", .fmt_count_n(need))
    } else {
      ""
    }
  )
  if (identical(x$design, "rotating")) {
    out <- paste0(out, sprintf(
      "%s in sample across %d live cohorts\n",
      .fmt_count_n(shown$in_sample), x$n_cohorts
    ))
  }
  out
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
    sprintf("%.0f%% of the life's loss at wave 1", 100 * share)
  }
  sprintf(
    "recruitment response %.3g, retention %s (%s)\n",
    x$params$resp_rate, ret_txt, loss_txt
  )
}

#' @keywords internal
#' @noRd
.fmt_panel_precision <- function(x) {
  label <- if (is.null(x$method)) x$type else sprintf("%s (%s)", x$type, x$method)
  sprintf(
    "%s: se = %.4g, moe = %.4g%s\n", label, x$se, x$moe,
    if (is.na(x$cv)) "" else sprintf(", cv = %.3g", x$cv)
  )
}

#' @keywords internal
#' @noRd
.fmt_panel_waves <- function(x) {
  w <- x$waves
  out <- data.frame(
    wave = w$wave,
    retention = ifelse(is.na(w$retention), "", sprintf("%.3g", w$retention)),
    q = sprintf("%.4g", w$q),
    n_resp = format(round(.panel_shown(x)$wave)),
    se = sprintf("%.4g", w$se),
    moe = sprintf("%.4g", w$moe),
    stringsAsFactors = FALSE
  )
  if (!all(is.na(w$cv))) {
    out$cv <- sprintf("%.3g", w$cv)
  }
  out
}

#' @rdname print.svyplan_panel
#' @export
format.svyplan_panel <- function(x, ...) {
  .check_unused_dots(...)
  sprintf(
    "svyplan_panel [%s, %d waves, recruit %g]",
    x$design, nrow(x$waves), .panel_recruit(x)
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
      what, from
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
