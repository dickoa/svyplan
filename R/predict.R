#' Grid Exploration for svyplan Objects
#'
#' Evaluate a svyplan result at new parameter combinations.
#' Returns a data frame with the varied parameters and resulting
#' quantities, suitable for sensitivity analysis or plotting.
#'
#' @param object A svyplan object (`svyplan_n`, `svyplan_cluster`,
#'   `svyplan_power`, or `svyplan_prec`).
#' @param newdata A data frame of parameter combinations to evaluate.
#'   Column names must be valid parameters for the object type (see
#'   Details). Parameters not in `newdata` stay at their original
#'   values from the object.
#' @param ... Ignored.
#'
#' @return A data frame with `newdata` columns followed by result
#'   columns. The result columns depend on the object type:
#'
#'   - `svyplan_n`: `n`, `se`, `moe`, `cv`
#'   - `svyplan_cluster`: `n1`, `n2`, (opt. `n3`), `total_n`, `cv`, `cost`
#'   - `svyplan_power`: `n`, `power`, `delta`
#'   - `svyplan_prec`: `se`, `moe`, `cv`
#'
#' @details
#' Valid parameters for `newdata` by object type:
#'
#' - **`n_prop`**: `p`, `moe`, `cv`, `alpha`, `N`, `deff`, `resp_rate`
#' - **`n_mean`**: `var`, `mu`, `moe`, `cv`, `alpha`, `N`, `deff`,
#'   `resp_rate`
#' - **`n_cluster`**: `cv`, `budget`, `resp_rate`, `fixed_cost`
#' - **`power_prop`**: `p1`, `p2`, `n`, `power`, `alpha`, `N`, `deff`,
#'   `sides`, `overlap`, `rho`, `resp_rate` (excluding the solved-for
#'   parameter)
#' - **`power_mean`**: `delta`, `var`, `n`, `power`, `alpha`, `N`,
#'   `deff`, `sides`, `overlap`, `rho`, `resp_rate` (excluding the
#'   solved-for parameter)
#' - **`prec_prop`**: `p`, `n`, `alpha`, `N`, `deff`, `resp_rate`
#' - **`prec_mean`**: `var`, `n`, `mu`, `alpha`, `N`, `deff`, `resp_rate`
#'
#' For `svyplan_n` objects, `moe` and `cv` are mutually exclusive in
#' `newdata`. If one appears, that mode is used; if neither, the
#' original mode is preserved.
#'
#' Similarly, for `svyplan_cluster` objects, `cv` and `budget` are
#' mutually exclusive.
#'
#' Multi-indicator (`n_multi`) and multi-indicator cluster results are
#' not supported; use the underlying single-indicator functions instead.
#'
#' If evaluation fails for a particular row (e.g. invalid parameter
#' combinations), that row's result columns are `NA` and a warning is
#' issued.
#'
#' @examples
#' # Sensitivity of sample size to deff and response rate
#' x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
#' predict(x, expand.grid(
#'   deff = seq(1, 3, 0.5),
#'   resp_rate = c(0.7, 0.8, 0.9)
#' ))
#'
#' # Power curve: how does power vary with sample size?
#' pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
#' predict(pw, data.frame(n = seq(100, 1000, 100)))
#'
#' @name predict.svyplan
NULL

#' @rdname predict.svyplan
#' @export
predict.svyplan_n <- function(object, newdata, ...) {
  if (!is.null(object$targets)) {
    stop(
      "predict() is not supported for multi-indicator results; ",
      "use the underlying single-indicator functions instead",
      call. = FALSE
    )
  }

  if (object$type == "proportion") {
    allowed <- c("p", "moe", "cv", "alpha", "N", "deff", "resp_rate")
    base <- object$params
    method <- object$method %||% "wald"

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv")

    .predict_grid(newdata, base, function(p) {
      res <- n_prop.default(
        p = p$p, moe = p$moe, cv = p$cv,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        method = method
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv)
    })

  } else if (object$type == "mean") {
    allowed <- c("var", "mu", "moe", "cv", "alpha", "N", "deff", "resp_rate")
    base <- object$params

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv")

    .predict_grid(newdata, base, function(p) {
      res <- n_mean.default(
        var = p$var, mu = p$mu, moe = p$moe, cv = p$cv,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv)
    })

  } else {
    stop(
      sprintf(
        "predict() is not supported for svyplan_n of type '%s'",
        object$type
      ),
      call. = FALSE
    )
  }
}

#' @rdname predict.svyplan
#' @export
predict.svyplan_cluster <- function(object, newdata, ...) {
  if (!is.null(object$targets)) {
    stop(
      "predict() is not supported for multi-indicator cluster results; ",
      "use n_cluster() directly",
      call. = FALSE
    )
  }

  allowed <- c("cv", "budget", "resp_rate", "fixed_cost")
  .validate_newdata(newdata, allowed)

  p <- object$params
  has_cv_col <- "cv" %in% names(newdata)
  has_budget_col <- "budget" %in% names(newdata)
  if (has_cv_col && has_budget_col) {
    stop("newdata cannot contain both 'cv' and 'budget'", call. = FALSE)
  }

  base <- list(
    cost = p$cost, delta = p$delta, rel_var = p$rel_var,
    k = p$k, resp_rate = p$resp_rate %||% 1, m = p$m,
    fixed_cost = p$fixed_cost %||% 0
  )

  if (has_cv_col) {
    # cv mode from newdata
  } else if (has_budget_col) {
    # budget mode from newdata
  } else if (!is.null(p$cv)) {
    base$cv <- p$cv
  } else if (!is.null(p$budget)) {
    base$budget <- p$budget
  } else {
    base$cv <- object$cv
  }

  .predict_grid(newdata, base, function(params) {
    res <- n_cluster.default(
      cost = params$cost, delta = params$delta,
      rel_var = params$rel_var, k = params$k,
      cv = params$cv, budget = params$budget,
      m = params$m, resp_rate = params$resp_rate,
      fixed_cost = params$fixed_cost
    )
    out <- as.list(res$n)
    out$total_n <- res$total_n
    out$cv <- res$cv
    out$cost <- res$cost
    as.data.frame(out)
  })
}

#' @rdname predict.svyplan
#' @export
predict.svyplan_power <- function(object, newdata, ...) {
  solved <- object$solved

  if (object$type == "proportion") {
    all_params <- c("p1", "p2", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "sides", "overlap", "rho")
    excluded <- switch(solved, n = "n", power = "power", mde = "p2")
    allowed <- setdiff(all_params, excluded)

    .validate_newdata(newdata, allowed)
    base <- object$params

    .predict_grid(newdata, base, function(p) {
      args <- list(
        p1 = p$p1, p2 = p$p2, n = p$n, power = p$power,
        alpha = p$alpha, N = p$N, deff = p$deff,
        resp_rate = p$resp_rate,
        sides = p$sides, overlap = p$overlap, rho = p$rho
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_prop, args)
      data.frame(n = res$n, power = res$power, delta = res$delta)
    })

  } else if (object$type == "mean") {
    all_params <- c("delta", "var", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "sides", "overlap", "rho")
    excluded <- switch(solved, n = "n", power = "power", mde = "delta")
    allowed <- setdiff(all_params, excluded)

    .validate_newdata(newdata, allowed)
    base <- object$params

    .predict_grid(newdata, base, function(p) {
      args <- list(
        delta = p$delta, var = p$var, n = p$n, power = p$power,
        alpha = p$alpha, N = p$N, deff = p$deff,
        resp_rate = p$resp_rate,
        sides = p$sides, overlap = p$overlap, rho = p$rho
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_mean, args)
      data.frame(n = res$n, power = res$power, delta = res$delta)
    })

  } else {
    stop(
      sprintf(
        "predict() is not supported for svyplan_power of type '%s'",
        object$type
      ),
      call. = FALSE
    )
  }
}

#' @rdname predict.svyplan
#' @export
predict.svyplan_prec <- function(object, newdata, ...) {
  if (object$type %in% c("multi", "cluster")) {
    stop(
      sprintf(
        "predict() is not supported for svyplan_prec of type '%s'",
        object$type
      ),
      call. = FALSE
    )
  }

  if (object$type == "proportion") {
    allowed <- c("p", "n", "alpha", "N", "deff", "resp_rate")
    base <- object$params
    method <- object$method %||% "wald"

    .validate_newdata(newdata, allowed)

    .predict_grid(newdata, base, function(p) {
      res <- prec_prop.default(
        p = p$p, n = p$n, alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        method = method
      )
      data.frame(se = res$se, moe = res$moe, cv = res$cv)
    })

  } else if (object$type == "mean") {
    allowed <- c("var", "n", "mu", "alpha", "N", "deff", "resp_rate")
    base <- object$params

    .validate_newdata(newdata, allowed)

    .predict_grid(newdata, base, function(p) {
      res <- prec_mean.default(
        var = p$var, n = p$n, mu = p$mu, alpha = p$alpha,
        N = p$N, deff = p$deff, resp_rate = p$resp_rate
      )
      data.frame(se = res$se, moe = res$moe, cv = res$cv)
    })

  } else {
    stop(
      sprintf(
        "predict() is not supported for svyplan_prec of type '%s'",
        object$type
      ),
      call. = FALSE
    )
  }
}

#' @keywords internal
#' @noRd
.validate_newdata <- function(newdata, allowed) {
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame", call. = FALSE)
  }
  if (nrow(newdata) == 0L) {
    stop("'newdata' must have at least one row", call. = FALSE)
  }
  if (ncol(newdata) == 0L) {
    stop("'newdata' must have at least one column", call. = FALSE)
  }

  unknown <- setdiff(names(newdata), allowed)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "unknown parameter(s) in newdata: %s\nvalid parameters: %s",
        paste(sQuote(unknown), collapse = ", "),
        paste(sQuote(allowed), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  non_numeric <- names(newdata)[!vapply(newdata, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "non-numeric columns in newdata: %s",
        paste(sQuote(non_numeric), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.resolve_exclusive <- function(newdata, base, name1, name2) {
  has1 <- name1 %in% names(newdata)
  has2 <- name2 %in% names(newdata)
  if (has1 && has2) {
    stop(
      sprintf("newdata cannot contain both '%s' and '%s'", name1, name2),
      call. = FALSE
    )
  }
  if (has1) {
    base[[name2]] <- NULL
  } else if (has2) {
    base[[name1]] <- NULL
  }
  base
}

#' @keywords internal
#' @noRd
.predict_grid <- function(newdata, base, eval_fn) {
  nrows <- nrow(newdata)
  nd_names <- names(newdata)
  results <- vector("list", nrows)

  for (i in seq_len(nrows)) {
    params <- base
    for (nm in nd_names) {
      params[[nm]] <- newdata[[nm]][i]
    }

    results[[i]] <- tryCatch(
      eval_fn(params),
      error = function(e) {
        warning(
          sprintf("predict row %d failed: %s", i, conditionMessage(e)),
          call. = FALSE
        )
        NULL
      }
    )
  }

  ok <- !vapply(results, is.null, logical(1))
  if (!any(ok)) {
    stop("all rows failed evaluation", call. = FALSE)
  }

  template <- results[[which(ok)[1L]]]
  result_names <- names(template)
  na_vals <- as.list(rep(NA_real_, length(result_names)))
  names(na_vals) <- result_names
  na_row <- as.data.frame(na_vals)

  result_df <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r)) na_row else r
  }))
  rownames(result_df) <- NULL

  dup_cols <- intersect(names(result_df), nd_names)
  if (length(dup_cols) > 0L) {
    result_df <- result_df[, !names(result_df) %in% dup_cols, drop = FALSE]
  }

  out <- cbind(newdata, result_df)
  rownames(out) <- NULL
  out
}
