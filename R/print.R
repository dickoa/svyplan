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
#' @keywords internal
#' @noRd
.fmt_deff <- function(deff) {
  if (is.null(deff)) return(NULL)
  if (isTRUE(all.equal(unname(deff), 1))) "deff = 1" else sprintf("deff = %.2f", deff)
}

#' @rdname print.svyplan
#' @export
print.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  if (x$type == "alloc") {
    .print_alloc_n(x)
  } else if (x$type == "multi") {
    .print_multi_n(x)
  } else {
    .print_single_n(x)
  }
  invisible(x)
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
    parts <- c(parts, sprintf("p = %.2f", p$p))
  }
  if (!is.null(p$var)) {
    parts <- c(parts, sprintf("var = %.2f", p$var))
  }
  if (!is.null(p$moe)) {
    parts <- c(parts, sprintf("moe = %.3f", p$moe))
  }
  if (!is.null(p$cv)) {
    parts <- c(parts, sprintf("cv = %.3f", p$cv))
  }
  parts <- c(parts, .fmt_deff(p$deff))
  if (!is.null(resp_rate) && resp_rate < 1) {
    parts <- c(parts, sprintf("resp_rate = %.2f", resp_rate))
  }
  if (length(parts) > 0L) {
    cat(sprintf(" (%s)", paste(parts, collapse = ", ")))
  }
  cat("\n")

  if (!is.null(x$domains)) {
    cat(sprintf("Domains: %d\n", nrow(x$domains)))
  }
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
    cat(sprintf(
      "Multi-indicator sample size (%d domains%s)\n",
      nrow(x$domains),
      min_n_label
    ))
    cat(sprintf("n = %d (binding: %s)\n", ceiling(x$n), x$binding))
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
  resp_rate <- x$params$resp_rate
  if (!is.null(resp_rate) && resp_rate < 1) {
    net_total <- ceiling(total_display * resp_rate)
    cat(sprintf(
      " -> total n = %s (net: %s)\n",
      .fmt_count_n(total_display), .fmt_count_n(net_total)
    ))
  } else {
    cat(sprintf(" -> total n = %s\n", .fmt_count_n(total_display)))
  }
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
  } else {
    .print_single_prec(x)
  }
  invisible(x)
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
    resp_rate <- p$resp_rate
    if (!is.null(resp_rate) && resp_rate < 1) {
      cat(sprintf(" (net: %s)", .fmt_count_n(ceiling(total_display * resp_rate))))
    }
    cat("\n")
  } else {
    method_label <- if (!is.null(x$method)) paste0(" (", x$method, ")") else ""
    cat(sprintf("Sampling precision for %s%s\n", type_label, method_label))
    cat(sprintf("n = %d", ceiling(p$n)))
    resp_rate <- p$resp_rate
    if (!is.null(resp_rate) && resp_rate < 1) {
      cat(sprintf(" (net: %d)", ceiling(p$n * resp_rate)))
    }
    cat("\n")
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
  cat("\n")
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
  parts <- c(parts, sprintf("alpha = %.2f", p$alpha))
  parts <- c(parts, .fmt_deff(p$deff))
  if (!is.null(resp_rate) && resp_rate < 1) {
    parts <- c(parts, sprintf("resp_rate = %.2f", resp_rate))
  }
  if (!is.null(p$var) && length(p$var) == 4L) {
    parts <- c(parts, sprintf(
      "var = (%.2f, %.2f, %.2f, %.2f)", p$var[1], p$var[2], p$var[3], p$var[4]
    ))
  }
  if (!is.null(p$overlap) && p$overlap > 0) {
    parts <- c(parts, sprintf("overlap = %.2f", p$overlap))
    parts <- c(parts, sprintf("overlap_cor = %.2f", p$overlap_cor))
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
  paste0("svyplan_varcomp [", x$stages, "-stage]")
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
  z <- qnorm(1 - alpha / 2)
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
  } else {
    stop("confint not supported for this type", call. = FALSE)
  }

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  moe <- z * object$se

  lo <- est - moe
  hi <- est + moe
  .ci_matrix(lo, hi, alpha)
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
  } else {
    stop("confint not supported for this precision type", call. = FALSE)
  }

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
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
  if (!is.null(resp_rate) && resp_rate < 1)
    parts <- c(parts, sprintf("resp_rate = %.2f", resp_rate))
  parts <- c(parts, .fmt_deff(p$deff))
  if (length(parts) > 0L)
    cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))
  if (!is.null(x$domains)) {
    cat(sprintf("Domains: %d\n", nrow(x$domains)))
    cat("---\n")
    dom <- x$domains
    if (".cv" %in% names(dom)) dom$.cv <- sprintf("%.4f", dom$.cv)
    if (".cost" %in% names(dom)) dom$.cost <- sprintf("%.0f", dom$.cost)
    print(dom, row.names = FALSE, right = FALSE)
  }
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

#' Assign Observations to Strata
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
        formatC(x$params$assurance, format = "f", digits = 2),
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
