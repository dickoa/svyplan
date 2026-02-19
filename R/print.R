#' Print svyplan objects
#'
#' @param x A svyplan object.
#' @param object A svyplan object (for `summary` methods).
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
  cat(sprintf("n = %d", ceiling(x$n)))

  p <- x$params
  parts <- character(0L)
  if (!is.null(p$p)) parts <- c(parts, sprintf("p = %.2f", p$p))
  if (!is.null(p$var)) parts <- c(parts, sprintf("var = %.2f", p$var))
  if (!is.null(p$moe)) parts <- c(parts, sprintf("moe = %.3f", p$moe))
  if (!is.null(p$cv)) parts <- c(parts, sprintf("cv = %.3f", p$cv))
  if (!is.null(p$deff) && p$deff != 1) {
    parts <- c(parts, sprintf("deff = %.2f", p$deff))
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
    cat(sprintf("Multi-indicator sample size (%d domains)\n", nrow(x$domains)))
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
  stage_parts <- vapply(seq_along(x$n), function(i) {
    sprintf("Stage %d: n%d = %d", i, i, ceiling(x$n[i]))
  }, character(1L))
  cat(paste(stage_parts, collapse = " | "))
  cat(sprintf(" -> total n = %d\n", ceiling(x$total_n)))
  cat(sprintf("CV = %.4f, cost = %.0f\n", x$cv, x$cost))

  if (!is.null(x$domains)) {
    cat(sprintf("Domains: %d\n", nrow(x$domains)))
  }
}

#' @keywords internal
#' @noRd
.print_multi_cluster <- function(x) {
  if (!is.null(x$domains)) {
    joint_label <- if (isTRUE(x$params$joint)) ", joint" else ""
    cat(sprintf("Multi-indicator optimal allocation (%d-stage, %d domains%s)\n",
                x$stages, nrow(x$domains), joint_label))
    cat(sprintf("Total n = %d\n", ceiling(x$total_n)))
    cat("---\n")
    dom <- x$domains
    dom$.total_n <- ceiling(dom$.total_n)
    for (s in seq_len(x$stages)) {
      col <- paste0("n", s)
      dom[[col]] <- ceiling(dom[[col]])
    }
    dom$.cv <- sprintf("%.4f", dom$.cv)
    dom$.cost <- sprintf("%.0f", dom$.cost)
    print(dom, row.names = FALSE, right = FALSE)
  } else {
    cat(sprintf("Multi-indicator optimal allocation (%d-stage)\n", x$stages))
    stage_parts <- vapply(seq_along(x$n), function(i) {
      sprintf("n%d = %d", i, ceiling(x$n[i]))
    }, character(1L))
    cat(paste(stage_parts, collapse = " | "))
    cat(sprintf(" -> total n = %d\n", ceiling(x$total_n)))
    cat(sprintf("CV = %.4f, cost = %.0f (binding: %s)\n",
                x$cv, x$cost, x$binding))
    cat("---\n")
    d <- x$detail
    d$.cv_achieved <- sprintf("%.4f", d$.cv_achieved)
    d$.binding <- ifelse(d$.binding, "*", "")
    print(d, row.names = FALSE, right = FALSE)
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
format.svyplan_n <- function(x, ...) {
  paste0("svyplan_n [n=", ceiling(x$n), ", ", x$type, "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_cluster <- function(x, ...) {
  paste0("svyplan_cluster [", x$stages, "-stage, total_n=",
         ceiling(x$total_n), "]")
}

#' @rdname print.svyplan
#' @export
format.svyplan_varcomp <- function(x, ...) {
  paste0("svyplan_varcomp [", x$stages, "-stage]")
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
summary.svyplan_varcomp <- function(object, ...) {
  print(object, ...)
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

  cat(sprintf("n = %d (per group), power = %.3f, delta = %.4f\n",
              ceiling(x$n), x$power, x$delta))

  p <- x$params
  parts <- character(0L)
  if (!is.null(p$p1)) parts <- c(parts, sprintf("p1 = %.3f", p$p1))
  if (!is.null(p$p2)) parts <- c(parts, sprintf("p2 = %.3f", p$p2))
  parts <- c(parts, sprintf("alpha = %.2f", p$alpha))
  if (!is.null(p$deff) && p$deff != 1) {
    parts <- c(parts, sprintf("deff = %.2f", p$deff))
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
format.svyplan_power <- function(x, ...) {
  paste0("svyplan_power [", x$type, ", ", x$solved,
         ", n=", ceiling(x$n), "]")
}

#' @rdname print.svyplan
#' @export
summary.svyplan_power <- function(object, ...) {
  print(object, ...)
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
  cat(sprintf("n = %d, CV = %.4f\n", ceiling(x$n), x$cv))
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
#'   `length(newdata)`. Observations outside the original range receive
#'   `NA`.
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
  s <- object$strata
  breaks <- c(s$lower[1L], object$boundaries, s$upper[object$n_strata])
  if (!is.null(labels) && length(labels) != object$n_strata) {
    stop(sprintf("'labels' must have length %d (number of strata)", object$n_strata),
         call. = FALSE)
  }
  cut(newdata, breaks = breaks, include.lowest = TRUE, labels = labels)
}
