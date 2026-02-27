#' Print svyplan objects
#'
#' @param x A svyplan object.
#' @param object A svyplan object (for `summary` and `confint` methods).
#' @param parm Ignored (included for S3 consistency with [confint()]).
#' @param level Confidence level (default 0.95).
#' @param ... Additional arguments (currently unused).
#'
#' @return `x` (or `object`), invisibly.
#'
#' @name print.svyplan
NULL

#' @rdname print.svyplan
#' @export
print.svyplan_n <- function(x, ...) {
  if (x$type == "multi") {
    .print_multi_n(x)
  } else {
    .print_single_n(x)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.print_single_n <- function(x) {
  type_label <- switch(x$type,
    proportion = "proportion",
    mean = "mean",
    x$type
  )
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
  if (!is.null(p$p)) parts <- c(parts, sprintf("p = %.2f", p$p))
  if (!is.null(p$var)) parts <- c(parts, sprintf("var = %.2f", p$var))
  if (!is.null(p$moe)) parts <- c(parts, sprintf("moe = %.3f", p$moe))
  if (!is.null(p$cv)) parts <- c(parts, sprintf("cv = %.3f", p$cv))
  if (!is.null(p$deff) && p$deff != 1) {
    parts <- c(parts, sprintf("deff = %.2f", p$deff))
  }
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

#' @keywords internal
#' @noRd
.print_multi_n <- function(x) {
  if (!is.null(x$domains)) {
    min_n_label <- if (!is.null(x$params$min_n)) {
      sprintf(", min_n = %g", x$params$min_n)
    } else {
      ""
    }
    cat(sprintf("Multi-indicator sample size (%d domains%s)\n",
                nrow(x$domains), min_n_label))
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
  if (!is.null(x$targets)) {
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
  n_display <- ceiling(x$n)
  total_display <- prod(n_display)
  stage_parts <- vapply(seq_along(n_display), function(i) {
    sprintf("Stage %d: n%d = %d", i, i, n_display[i])
  }, character(1L))
  cat(paste(stage_parts, collapse = " | "))
  resp_rate <- x$params$resp_rate
  if (!is.null(resp_rate) && resp_rate < 1) {
    net_total <- ceiling(total_display * resp_rate)
    cat(sprintf(" -> total n = %d (net: %d)\n", total_display, net_total))
  } else {
    cat(sprintf(" -> total n = %d\n", total_display))
  }
  fc <- x$params$fixed_cost
  if (!is.null(fc) && fc > 0) {
    cat(sprintf("cv = %.4f, cost = %.0f (fixed: %.0f)\n", x$cv, x$cost, fc))
  } else {
    cat(sprintf("cv = %.4f, cost = %.0f\n", x$cv, x$cost))
  }

  if (!is.null(x$domains)) {
    cat(sprintf("Domains: %d\n", nrow(x$domains)))
  }
}

#' @keywords internal
#' @noRd
.print_multi_cluster <- function(x) {
  if (!is.null(x$domains)) {
    joint_label <- if (isTRUE(x$params$joint)) ", joint" else ""
    min_n_label <- if (!is.null(x$params$min_n)) {
      sprintf(", min_n = %g", x$params$min_n)
    } else {
      ""
    }
    cat(sprintf("Multi-indicator optimal allocation (%d-stage, %d domains%s%s)\n",
                x$stages, nrow(x$domains), joint_label, min_n_label))
    cat("---\n")
    dom <- x$domains
    for (s in seq_len(x$stages)) {
      col <- paste0("n", s)
      dom[[col]] <- ceiling(dom[[col]])
    }
    dom$.total_n <- apply(dom[, paste0("n", seq_len(x$stages)), drop = FALSE], 1, prod)
    fc <- x$params$fixed_cost
    if (!is.null(fc) && fc > 0) {
      cat(sprintf("Total n = %d, cost = %.0f (fixed: %.0f)\n",
                  sum(dom$.total_n), x$cost, fc))
    } else {
      cat(sprintf("Total n = %d\n", sum(dom$.total_n)))
    }
    dom$.cv <- sprintf("%.4f", dom$.cv)
    dom$.cost <- sprintf("%.0f", dom$.cost)
    print(dom, row.names = FALSE, right = FALSE)
  } else {
    cat(sprintf("Multi-indicator optimal allocation (%d-stage)\n", x$stages))
    n_display <- ceiling(x$n)
    total_display <- prod(n_display)
    stage_parts <- vapply(seq_along(n_display), function(i) {
      sprintf("n%d = %d", i, n_display[i])
    }, character(1L))
    cat(paste(stage_parts, collapse = " | "))
    cat(sprintf(" -> total n = %d\n", total_display))
    fc <- x$params$fixed_cost
    if (!is.null(fc) && fc > 0) {
      cat(sprintf("cv = %.4f, cost = %.0f (fixed: %.0f, binding: %s)\n",
                  x$cv, x$cost, fc, x$binding))
    } else {
      cat(sprintf("cv = %.4f, cost = %.0f (binding: %s)\n",
                  x$cv, x$cost, x$binding))
    }
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
  if (x$type == "multi") {
    .print_multi_prec(x)
  } else {
    .print_single_prec(x)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.print_single_prec <- function(x) {
  type_label <- switch(x$type,
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
    stage_parts <- vapply(seq_along(n_display), function(i) {
      sprintf("n%d = %d", i, n_display[i])
    }, character(1L))
    cat(paste(stage_parts, collapse = " | "))
    cat(sprintf(" -> total n = %d", total_display))
    resp_rate <- p$resp_rate
    if (!is.null(resp_rate) && resp_rate < 1) {
      cat(sprintf(" (net: %d)", ceiling(total_display * resp_rate)))
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
    if (!is.na(x$se[1L])) cat(", ") else cat("")
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
  cat(sprintf("Variance components (%d-stage)\n", x$stages))
  cat(sprintf("var_between = %.4f", x$var_between))
  for (i in seq_along(x$var_within)) {
    cat(sprintf(", var_within%s = %.4f",
                if (length(x$var_within) > 1L) paste0("[", i, "]") else "",
                x$var_within[i]))
  }
  cat("\n")
  cat(sprintf("delta = %s\n", paste(sprintf("%.4f", x$delta), collapse = ", ")))
  cat(sprintf("k = %s\n", paste(sprintf("%.4f", x$k), collapse = ", ")))
  cat(sprintf("Unit relvariance = %.4f\n", x$rel_var))

  invisible(x)
}

#' @rdname print.svyplan
#' @export
print.svyplan_power <- function(x, ...) {
  type_label <- switch(x$type,
    proportion = "proportions",
    mean = "means",
    x$type
  )
  solved_label <- switch(x$solved,
    n = "sample size",
    power = "power",
    mde = "minimum detectable effect"
  )
  cat(sprintf("Power analysis for %s (solved for %s)\n", type_label, solved_label))

  p <- x$params
  resp_rate <- p$resp_rate
  if (!is.null(resp_rate) && resp_rate < 1) {
    net_n <- ceiling(x$n * resp_rate)
    cat(sprintf("n = %d (net: %d, per group), power = %.3f, delta = %.4f\n",
                ceiling(x$n), net_n, x$power, x$delta))
  } else {
    cat(sprintf("n = %d (per group), power = %.3f, delta = %.4f\n",
                ceiling(x$n), x$power, x$delta))
  }

  parts <- character(0L)
  if (!is.null(p$p1)) parts <- c(parts, sprintf("p1 = %.3f", p$p1))
  if (!is.null(p$p2)) parts <- c(parts, sprintf("p2 = %.3f", p$p2))
  parts <- c(parts, sprintf("alpha = %.2f", p$alpha))
  if (!is.null(p$deff) && p$deff != 1) {
    parts <- c(parts, sprintf("deff = %.2f", p$deff))
  }
  if (!is.null(resp_rate) && resp_rate < 1) {
    parts <- c(parts, sprintf("resp_rate = %.2f", resp_rate))
  }
  if (!is.null(p$overlap) && p$overlap > 0) {
    parts <- c(parts, sprintf("overlap = %.2f", p$overlap))
    parts <- c(parts, sprintf("rho = %.2f", p$rho))
  }
  if (!is.null(p$sides) && p$sides == 1) {
    parts <- c(parts, "one-sided")
  }
  cat(sprintf("(%s)\n", paste(parts, collapse = ", ")))

  invisible(x)
}

#' @rdname print.svyplan
#' @export
format.svyplan_n <- function(x, ...) {
  paste0("svyplan_n [n=", ceiling(x$n), ", ", x$type, "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_cluster <- function(x, ...) {
  paste0("svyplan_cluster [", x$stages, "-stage, total_n=",
         prod(ceiling(x$n)), "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_prec <- function(x, ...) {
  paste0("svyplan_prec [", x$type, ", se=", sprintf("%.4f", x$se[1L]), "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_varcomp <- function(x, ...) {
  paste0("svyplan_varcomp [", x$stages, "-stage]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_power <- function(x, ...) {
  paste0("svyplan_power [", x$type, ", ", x$solved,
         ", n=", ceiling(x$n), "]")
}

#' @rdname print.svyplan
#' @export
summary.svyplan_n <- function(object, ...) {
  print(object, ...)
}

#' @rdname print.svyplan
#' @export
summary.svyplan_cluster <- function(object, ...) {
  print(object, ...)
}

#' @rdname print.svyplan
#' @export
summary.svyplan_prec <- function(object, ...) {
  print(object, ...)
}

#' @rdname print.svyplan
#' @export
summary.svyplan_varcomp <- function(object, ...) {
  print(object, ...)
}

#' @rdname print.svyplan
#' @export
summary.svyplan_power <- function(object, ...) {
  print(object, ...)
}

#' @rdname print.svyplan
#' @export
confint.svyplan_n <- function(object, parm, level = 0.95, ...) {
  if (object$type == "multi" || is.na(object$se)) {
    stop("confint requires a single-indicator svyplan_n object with computed SE",
         call. = FALSE)
  }

  p <- object$params
  if (object$type == "proportion") {
    est <- p$p
  } else if (object$type == "mean") {
    est <- p$mu
    if (is.null(est)) {
      stop("'mu' is required to compute a confidence interval for a mean",
           call. = FALSE)
    }
  } else {
    stop("confint not supported for this type", call. = FALSE)
  }

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  moe <- z * object$se

  lo <- est - moe
  hi <- est + moe
  if (object$type == "proportion") {
    lo <- max(lo, 0)
    hi <- min(hi, 1)
  }

  ci <- matrix(c(lo, hi), nrow = 1L,
               dimnames = list("", c(
                 sprintf("%.1f %%", 100 * alpha / 2),
                 sprintf("%.1f %%", 100 * (1 - alpha / 2))
               )))
  ci
}

#' @rdname print.svyplan
#' @export
confint.svyplan_prec <- function(object, parm, level = 0.95, ...) {
  p <- object$params
  if (object$type == "proportion") {
    est <- p$p
  } else if (object$type == "mean") {
    est <- p$mu
    if (is.null(est)) {
      stop("'mu' is required to compute a confidence interval for a mean",
           call. = FALSE)
    }
  } else {
    stop("confint not supported for this precision type", call. = FALSE)
  }

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  moe <- z * object$se

  lo <- est - moe
  hi <- est + moe
  if (object$type == "proportion") {
    lo <- max(lo, 0)
    hi <- min(hi, 1)
  }

  ci <- matrix(c(lo, hi), nrow = 1L,
               dimnames = list("", c(
                 sprintf("%.1f %%", 100 * alpha / 2),
                 sprintf("%.1f %%", 100 * (1 - alpha / 2))
               )))
  ci
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_n <- function(x, ...) {
  as.integer(ceiling(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_n <- function(x, ...) {
  x$n
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_cluster <- function(x, ...) {
  as.integer(ceiling(x$total_n))
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_power <- function(x, ...) {
  as.integer(ceiling(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_power <- function(x, ...) {
  x$n
}

#' @rdname print.svyplan
#' @export
print.svyplan_strata <- function(x, ...) {
  method_label <- switch(x$method,
    cumrootf = "Dalenius-Hodges",
    geo      = "Geometric",
    lh       = "Lavallee-Hidiroglou",
    kozak    = "Kozak",
    x$method
  )
  cat(sprintf("Strata boundaries (%s, %d strata)\n", method_label, x$n_strata))
  cat(sprintf("Boundaries: %s\n",
              paste(sprintf("%.1f", x$boundaries), collapse = ", ")))
  cat(sprintf("n = %d, cv = %.4f\n", ceiling(x$n), x$cv))
  if (!is.null(x$alloc) && is.character(x$alloc)) {
    if (x$alloc == "power") {
      cat(sprintf("Allocation: power (q = %.2f)\n", x$params$q))
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
  df$W_h <- sprintf("%.3f", df$W_h)
  df$S_h <- sprintf("%.1f", df$S_h)
  df$certain <- NULL
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @rdname print.svyplan
#' @export
format.svyplan_strata <- function(x, ...) {
  paste0("svyplan_strata [", x$method, ", ", x$n_strata, " strata, n=",
         ceiling(x$n), "]")
}

#' @rdname print.svyplan
#' @export
summary.svyplan_strata <- function(object, ...) {
  print(object, ...)
}

#' @rdname print.svyplan
#' @export
as.data.frame.svyplan_strata <- function(x, ...) {
  x$strata
}

#' @rdname print.svyplan
#' @export
as.integer.svyplan_strata <- function(x, ...) {
  as.integer(ceiling(x$n))
}

#' @rdname print.svyplan
#' @export
as.double.svyplan_strata <- function(x, ...) {
  x$boundaries
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
#' @param ... Additional arguments (currently unused).
#'
#' @return A factor of stratum assignments with length equal to
#'   `length(newdata)`. Values beyond the original training range are
#'   assigned to the lowest or highest stratum.
#'
#' @seealso [strata_bound()] to compute the boundaries.
#'
#' @examples
#' set.seed(42)
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
  if (!is.numeric(newdata)) {
    stop("'newdata' must be a numeric vector", call. = FALSE)
  }
  breaks <- c(-Inf, object$boundaries, Inf)
  if (!is.null(labels) && length(labels) != object$n_strata) {
    stop(sprintf("'labels' must have length %d (number of strata)", object$n_strata),
         call. = FALSE)
  }
  cut(newdata, breaks = breaks, include.lowest = TRUE, labels = labels)
}
