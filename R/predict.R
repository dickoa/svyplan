#' Grid exploration for svyplan objects
#'
#' Evaluate a svyplan result at new parameter combinations.
#' Returns a data frame with the varied parameters and resulting
#' quantities, suitable for sensitivity analysis or plotting.
#'
#' @param object A svyplan object (`svyplan_n`, `svyplan_cluster`,
#'   `svyplan_power`, or `svyplan_prec`). For `svyplan_prec`, only
#'   types `"proportion"` and `"mean"` are supported (not `"cluster"` or
#'   `"multi"`).
#' @param newdata A data frame of parameter combinations to evaluate.
#'   Column names must be valid parameters for the object type (see
#'   Details). Parameters not in `newdata` stay at their original
#'   values from the object.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @return A data frame with `newdata` columns followed by result
#'   columns. The result columns depend on the object type:
#'
#'   - `svyplan_n`: `n`, `se`, `moe`, `cv`, `rmoe`, except for a
#'     `type = "change"` result, which reports one size per occasion as
#'     `n1` and `n2` in place of `n`
#'   - `svyplan_cluster`: `n_psu`, `n_per_psu`, (opt. `n_per_ssu`), `total_n`, `cv`, `cost`
#'   - `svyplan_power`: `n`, `power`, `effect`
#'   - `svyplan_prec`: `se`, `moe`, `cv`, `rmoe`
#'
#' @details
#' Valid parameters for `newdata` by object type:
#'
#' - **`n_prop`**: `p`, `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`,
#'   `resp_rate`, `df` (`method = "beta"` only)
#' - **`n_mean`**: `var`, `mu`, `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`,
#'   `resp_rate`
#' - **`n_change`**: `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `ratio`, `overlap`, `overlap_cor`, and `change` on the mean scale only,
#'   since two proportions determine it. The occasion variances are held at
#'   the values the result was built from, so vary the design rather than the
#'   estimand: sweeping `overlap` and `overlap_cor` is what prices a rotation
#'   against a fresh sample each round.
#' - **`n_pooled`**: `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `mu` on the mean scale, and `occasions` only where the overlap and
#'   correlation are flat across lags, a shaped profile being a statement
#'   about a horizon that changing the horizon would contradict. The lag
#'   profiles themselves cannot be varied, a grid row having nowhere to hold
#'   a vector.
#' - **`n_cluster`**: `cv`, `budget`, `unit_relvar`, `resp_rate`, `fixed_cost`,
#'   stage deltas (`icc` or `icc_psu`, plus `icc_ssu` for 3-stage),
#'   stage ratios (`var_ratio` or `var_ratio_psu`, plus `var_ratio_ssu` for 3-stage),
#'   and stage costs (`cost_psu`, `cost_ssu`, `cost_tsu`). For 2-stage
#'   designs, `cost_tsu` aliases `cost_ssu`.
#' - **`power_prop`**: `p1`, `p2`, `n`, `power`, `alpha`, `N`, `deff`,
#'   `alternative`, `overlap`, `overlap_cor`, `resp_rate` (excluding the solved-for
#'   parameter). Not supported for objects with vector `n`.
#' - **`power_mean`**: `effect`, `var`, `n`, `power`, `alpha`, `N`,
#'   `deff`, `alternative`, `overlap`, `overlap_cor`, `resp_rate` (excluding the
#'   solved-for parameter). Not supported for objects with vector `n`.
#' - **`prec_prop`**: `p`, `n`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `df` (`method = "beta"` only). A result that solved for `p` varies the
#'   target it solved from, `cv` or `rmoe`, in place of `p`.
#' - **`prec_mean`**: `var`, `n`, `mu`, `alpha`, `N`, `deff`, `resp_rate`
#'
#' For `svyplan_n` objects, `moe`, `rmoe`, and `cv` are mutually exclusive
#' in `newdata`. If one appears, that mode is used. If none appears, the
#' original mode is preserved.
#'
#' Similarly, for `svyplan_cluster` objects, `cv` and `budget` are
#' mutually exclusive.
#'
#' Multi-indicator (`n_multi`) and multi-indicator cluster results are
#' not supported. Use the underlying single-indicator functions instead.
#'
#' A joint constrained allocation ([n_alloc()] with `measures` and `targets`)
#' is supported only in fixed-budget objective mode, where `newdata` varies
#' `budget` alone and the result is the cost-versus-objective frontier: `n`,
#' continuous `cost`, `objective_value` and its equivalent `cv`, the
#' operational `n_int` and `cost_int`, whether the budget binds, and
#' `.feasible`. Budgets that cannot fund the hard targets give an all-`NA` row
#' with `.feasible = FALSE` and a warning, so one infeasible point does not
#' discard the rest of the frontier. In minimum-cost mode there is no scalar to
#' vary; modify `targets` and call [n_alloc()] again.
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
#' # Cluster design: sensitivity to icc (homogeneity)
#' cl <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
#' predict(cl, data.frame(icc = c(0.01, 0.03, 0.05, 0.10, 0.15)))
#'
#' # Allocation: how does the CV change with sample size?
#' frame <- data.frame(
#'   N    = c(4000, 3000, 3000),
#'   sd   = c(10, 15, 8),
#'   mean = c(50, 60, 55)
#' )
#' alloc <- n_alloc(frame, n = 600)
#' predict(alloc, data.frame(n = seq(200, 1000, 200)))
#'
#' @seealso [plot.svyplan] to draw a one-parameter grid, and
#'   [confint.svyplan] for the interval implied by a single result.
#'
#' @name predict.svyplan
NULL

#' @rdname predict.svyplan
#' @export
predict.svyplan_n <- function(object, newdata, ...) {
  .check_unused_dots(...)
  if (identical(object$method, "bethel")) {
    if (!identical(object$params$mode, "budget_objective")) {
      stop(
        "predict() is not supported for joint constrained allocations; modify 'targets' and rerun n_alloc()",
        call. = FALSE
      )
    }
    return(.predict_bethel_budget(object, newdata))
  }
  if (!is.null(object$indicators)) {
    stop(
      "predict() is not supported for multi-indicator results; ",
      "use the underlying single-indicator functions instead",
      call. = FALSE
    )
  }

  if (object$type == "proportion") {
    allowed <- c("p", "moe", "cv", "rmoe", "alpha", "N", "deff", "resp_rate",
                 "df", "min_cases")
    base <- object$params
    method <- object$method %||% "wald"

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- n_prop.default(
        p = p$p, moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        method = method, df = p$df, min_cases = p$min_cases
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
    })

  } else if (object$type == "mean") {
    allowed <- c("var", "mu", "moe", "cv", "rmoe", "alpha", "N", "deff",
                 "resp_rate")
    base <- object$params

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- n_mean.default(
        var = p$var, mu = p$mu, moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
    })

  } else if (object$type == "change") {
    base <- object$params
    # 'change' is p[2] - p[1] on the proportion scale, so varying it there
    # would ask for a change the two proportions contradict.
    allowed <- c("moe", "cv", "rmoe", "alpha", "N", "deff", "resp_rate",
                 "ratio", "overlap", "overlap_cor",
                 if (is.null(base$p)) "change")

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- n_change.default(
        var = if (is.null(p$p)) p$var else NULL,
        p = p$p,
        change = if (is.null(p$p)) p$change else NULL,
        moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N, deff = p$deff, resp_rate = p$resp_rate,
        ratio = p$ratio, overlap = p$overlap, overlap_cor = p$overlap_cor,
        df = p$df
      )
      data.frame(n1 = res$n[1L], n2 = res$n[length(res$n)],
                 se = res$se, moe = res$moe, cv = res$cv, rmoe = res$rmoe)
    })

  } else if (object$type == "pooled") {
    base <- object$params
    # The lag profiles are stored one entry per lag, so a grid can only vary
    # 'occasions' where they are flat and a shorter or longer profile is the
    # same statement. Varying the profiles themselves would need a column
    # holding a vector, which a data frame row cannot.
    # An issued profile is only a respondent one at full response, so a
    # grid cannot move the response rate under it any more than a direct
    # call could.
    if (identical(base$overlap_basis, "issued") &&
          "resp_rate" %in% names(newdata)) {
      .check_overlap_basis("issued", min(newdata$resp_rate))
    }
    flat <- length(unique(base$overlap)) == 1L &&
      length(unique(base$overlap_cor)) == 1L
    allowed <- c("moe", "cv", "rmoe", "alpha", "N", "deff", "resp_rate",
                 if (flat) "occasions",
                 if (is.null(base$p)) "mu")

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")
    if (flat) {
      base$overlap <- base$overlap[1L]
      base$overlap_cor <- base$overlap_cor[1L]
    }

    .predict_grid(newdata, base, function(p) {
      res <- n_pooled.default(
        var = if (is.null(p$p)) p$var else NULL,
        p = p$p,
        mu = if (is.null(p$p)) p$mu else NULL,
        occasions = p$occasions,
        moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N, deff = p$deff, resp_rate = p$resp_rate,
        overlap = p$overlap, overlap_cor = p$overlap_cor, df = p$df,
        .overlap_basis = p$overlap_basis
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
    })

  } else if (object$type == "alloc") {
    allowed <- c("n", "cv", "budget", "alpha", "deff", "resp_rate", "min_n_stratum", "alloc_q")
    base <- object$params
    base_small <- list(
      n = base$n, cv = base$cv, budget = base$budget,
      alpha = base$alpha, deff = base$deff, resp_rate = base$resp_rate,
      min_n_stratum = base$min_n_stratum, alloc_q = base$alloc_q %||% 0.5
    )

    .validate_newdata(newdata, allowed)

    has_n <- "n" %in% names(newdata)
    has_cv <- "cv" %in% names(newdata)
    has_budget <- "budget" %in% names(newdata)
    if ((has_n + has_cv + has_budget) > 1L) {
      stop("newdata can vary at most one of 'n', 'cv', or 'budget'",
           call. = FALSE)
    }
    if (has_n) {
      base_small$cv <- NULL
      base_small$budget <- NULL
    } else if (has_cv) {
      base_small$n <- NULL
      base_small$budget <- NULL
    } else if (has_budget) {
      base_small$n <- NULL
      base_small$cv <- NULL
    }

    .predict_grid(newdata, base_small, function(p) {
      count_mode <- (!is.null(p$n)) + (!is.null(p$cv)) + (!is.null(p$budget))
      if (count_mode != 1L) {
        stop("each predict row must define exactly one of n/cv/budget",
             call. = FALSE)
      }
      alloc_args <- list(
        frame = base$frame,
        n = p$n, cv = p$cv, budget = p$budget,
        alloc = base$alloc %||% "neyman",
        alpha = p$alpha,
        deff = p$deff,
        resp_rate = p$resp_rate,
        min_n_stratum = p$min_n_stratum,
        alloc_q = p$alloc_q
      )
      if (!.alloc_is_cluster(base$frame)) {
        alloc_args$unit_cost <- base$cost_h
      }
      res <- do.call(n_alloc.default, alloc_args)
      data.frame(
        n = res$n, se = res$se, moe = res$moe, cv = res$cv,
        cost = res$params$achieved$cost
      )
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
  .check_unused_dots(...)
  if (!is.null(object$indicators)) {
    stop(
      "predict() is not supported for multi-indicator cluster results; ",
      "use n_cluster() directly",
      call. = FALSE
    )
  }

  stages <- length(object$params$stage_cost)
  if (stages == 3L && "icc" %in% names(newdata)) {
    stop(
      "for 3-stage cluster predict, vary 'icc_psu' and 'icc_ssu' instead of 'icc'",
      call. = FALSE
    )
  }
  if (stages == 3L && "var_ratio" %in% names(newdata)) {
    stop(
      "for 3-stage cluster predict, vary 'var_ratio_psu' and 'var_ratio_ssu' instead of 'var_ratio'",
      call. = FALSE
    )
  }

  cost_meta <- .cluster_cost_col_map(names(newdata), stages)
  icc_meta <- .cluster_stage_col_map(
    names(newdata),
    "icc",
    stage_count = stages - 1L,
    allow_scalar_alias = stages == 2L
  )
  k_meta <- .cluster_stage_col_map(
    names(newdata),
    "var_ratio",
    stage_count = stages - 1L,
    allow_scalar_alias = stages == 2L
  )

  allowed <- unique(c(
    "cv", "budget", "unit_relvar", "resp_rate_psu", "resp_rate_ssu",
    "resp_rate", "fixed_cost",
    "n_psu", "n_per_psu", "n_per_ssu",
    icc_meta$allowed, k_meta$allowed, cost_meta$allowed
  ))
  .validate_newdata(newdata, allowed)

  p <- object$params
  has_cv_col <- "cv" %in% names(newdata)
  has_budget_col <- "budget" %in% names(newdata)
  if (has_cv_col && has_budget_col) {
    stop("newdata cannot contain both 'cv' and 'budget'", call. = FALSE)
  }

  base <- list(
    stage_cost = p$stage_cost, icc = p$icc, unit_relvar = p$unit_relvar,
    var_ratio = p$var_ratio, resp_rate_psu = p$resp_rate_psu %||% 1,
    resp_rate_ssu = p$resp_rate_ssu %||% 1, resp_rate = p$resp_rate %||% 1,
    n_psu = p$n_psu,
    n_per_psu = p$n_per_psu, n_per_ssu = p$n_per_ssu,
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
    row_cost <- .apply_cluster_cost_cols(base$stage_cost, params, cost_meta$map)
    row_icc <- .apply_cluster_stage_cols(
      base$icc,
      params,
      icc_meta$map,
      icc_meta$canonical
    )
    row_k <- .apply_cluster_stage_cols(
      base$var_ratio,
      params,
      k_meta$map,
      k_meta$canonical
    )
    res <- n_cluster.default(
      stage_cost = row_cost, icc = row_icc,
      unit_relvar = params$unit_relvar, var_ratio = row_k,
      cv = params$cv, budget = params$budget,
      n_psu = params$n_psu, n_per_psu = params$n_per_psu,
      n_per_ssu = params$n_per_ssu, resp_rate_psu = params$resp_rate_psu,
      resp_rate_ssu = params$resp_rate_ssu %||% 1,
      resp_rate = params$resp_rate %||% 1,
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
  .check_unused_dots(...)
  if (length(object$n) == 2L)
    stop("predict() does not support power objects with unequal-group n; call power_*() directly",
         call. = FALSE)

  solved <- object$solved

  if (object$type == "proportion") {
    all_params <- c("p1", "p2", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "alternative", "overlap", "overlap_cor")
    excluded <- switch(solved, n = "n", power = "power", mde = "p2")
    allowed <- setdiff(all_params, excluded)

    .validate_newdata(newdata, allowed, allowed_character = "alternative")
    base <- object$params
    method <- base$method %||% "wald"

    .predict_grid(newdata, base, function(p) {
      args <- list(
        p1 = p$p1, p2 = p$p2, n = p$n, power = p$power,
        alpha = p$alpha, N = p$N, deff = p$deff,
        resp_rate = p$resp_rate,
        alternative = p$alternative, overlap = p$overlap, overlap_cor = p$overlap_cor,
        method = method
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_prop, args)
      data.frame(n = res$n, power = res$power, effect = res$effect)
    })

  } else if (object$type == "mean") {
    all_params <- c("effect", "var", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "alternative", "overlap", "overlap_cor")
    excluded <- switch(solved, n = "n", power = "power", mde = "effect")
    allowed <- setdiff(all_params, excluded)

    .validate_newdata(newdata, allowed, allowed_character = "alternative")
    base <- object$params

    .predict_grid(newdata, base, function(p) {
      args <- list(
        effect = p$effect, var = p$var, n = p$n, power = p$power,
        alpha = p$alpha, N = p$N, deff = p$deff,
        resp_rate = p$resp_rate,
        alternative = p$alternative, overlap = p$overlap, overlap_cor = p$overlap_cor
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_mean, args)
      data.frame(n = res$n, power = res$power, effect = res$effect)
    })

  } else if (object$type %in% c("did_prop", "did_mean")) {
    all_params <- c("effect", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "alternative", "overlap", "overlap_cor", "ratio")
    excluded <- switch(solved, n = "n", power = "power", mde = "effect")
    allowed <- setdiff(all_params, excluded)

    .validate_newdata(newdata, allowed, allowed_character = "alternative")
    base <- object$params

    .predict_grid(newdata, base, function(p) {
      args <- list(
        treat = p$treat, control = p$control,
        outcome = p$outcome, var = p$var,
        effect = p$effect, n = p$n, power = p$power,
        alpha = p$alpha, N = p$N, deff = p$deff,
        resp_rate = p$resp_rate,
        alternative = p$alternative, ratio = p$ratio,
        overlap = p$overlap, overlap_cor = p$overlap_cor
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_did, args)
      data.frame(n = res$n, power = res$power, effect = res$effect)
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
  .check_unused_dots(...)
  if (object$type %in% c("multi", "cluster")) {
    stop(
      sprintf(
        "predict() is not supported for svyplan_prec of type '%s'",
        object$type
      ),
      call. = FALSE
    )
  }

  # A result that solved for the level varies 'cv' instead of the level, and
  # reports the level back, mirroring how power objects treat solved-for
  # parameters.
  solved <- object$solved

  if (object$type == "proportion") {
    allowed <- c("p", "n", "alpha", "N", "deff", "resp_rate", "df")
    base <- object$params
    method <- object$method %||% "wald"
    if (identical(solved, "p")) {
      # The grid varies whichever target the level was solved from, so an
      # 'rmoe' solve stays an 'rmoe' solve: re-solving it through 'cv'
      # would drop the interval method the target was stated against.
      metric <- if (is.null(base$rmoe)) "cv" else "rmoe"
      allowed <- c(setdiff(allowed, "p"), metric)
      base$p <- NULL
      base[[metric]] <- object[[metric]]
    }

    .validate_newdata(newdata, allowed)

    .predict_grid(newdata, base, function(p) {
      res <- prec_prop.default(
        p = p$p, n = p$n, cv = p$cv, rmoe = p$rmoe, alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        method = method, df = p$df
      )
      out <- data.frame(se = res$se, moe = res$moe, cv = res$cv,
                        rmoe = res$rmoe)
      if (identical(solved, "p")) out <- cbind(data.frame(p = res$params$p), out)
      out
    })

  } else if (object$type == "mean") {
    allowed <- c("var", "n", "mu", "alpha", "N", "deff", "resp_rate")
    base <- object$params
    if (identical(solved, "mu")) {
      allowed <- c(setdiff(allowed, "mu"), "cv")
      base$mu <- NULL
      base$cv <- object$cv
    }

    .validate_newdata(newdata, allowed)

    .predict_grid(newdata, base, function(p) {
      res <- prec_mean.default(
        var = p$var, n = p$n, mu = p$mu, cv = p$cv, alpha = p$alpha,
        N = p$N, deff = p$deff, resp_rate = p$resp_rate
      )
      out <- data.frame(se = res$se, moe = res$moe, cv = res$cv,
                        rmoe = res$rmoe)
      if (identical(solved, "mu")) {
        out <- cbind(data.frame(mu = res$params$mu), out)
      }
      out
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
.validate_newdata <- function(newdata, allowed, allowed_character = character(0)) {
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
    char_ok <- non_numeric[
      non_numeric %in% allowed_character &
      vapply(newdata[non_numeric], function(x) is.character(x) || is.factor(x), logical(1))
    ]
    bad <- setdiff(non_numeric, char_ok)
    if (length(bad) == 0L) return(invisible(TRUE))
    stop(
      sprintf(
        "non-numeric columns in newdata: %s",
        paste(sQuote(bad), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.resolve_exclusive <- function(newdata, base, ...) {
  names <- c(...)
  present <- names[names %in% names(newdata)]
  if (length(present) > 1L) {
    stop(
      sprintf("newdata cannot contain more than one of %s",
              paste(sQuote(names), collapse = ", ")),
      call. = FALSE
    )
  }
  if (length(present) == 1L) {
    for (dropped in setdiff(names, present)) base[[dropped]] <- NULL
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
      val <- newdata[[nm]][i]
      if (is.factor(val)) val <- as.character(val)
      params[[nm]] <- val
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

#' Budget frontier for a joint budget-objective allocation
#'
#' Re-solves the stored problem at each requested budget. The root search
#' already sweeps the minimum-cost problem across objective bounds, so the
#' cost-versus-objective frontier costs little beyond the solves themselves.
#' Budgets that cannot fund the hard targets yield an all-`NA` row with
#' `.feasible = FALSE` rather than aborting the grid.
#' @keywords internal
#' @noRd
.predict_bethel_budget <- function(object, newdata) {
  .validate_newdata(newdata, "budget")
  if (!"budget" %in% names(newdata)) {
    stop("newdata must vary 'budget' for a joint budget-objective allocation",
         call. = FALSE)
  }
  p <- object$params
  targets <- if (is.null(p$targets) || nrow(p$targets) == 0L) NULL else
    p$targets
  rows <- lapply(newdata$budget, function(b) {
    fit <- tryCatch(
      n_alloc.default(
        frame = p$frame,
        measures = p$measures,
        targets = targets,
        objective = p$objective,
        budget = b,
        unit_cost = p$unit_cost,
        alpha = p$alpha,
        deff = p$deff,
        resp_rate = p$resp_rate,
        min_n_stratum = p$min_n_stratum
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      warning(sprintf("budget %s is infeasible: %s", format(b),
                      conditionMessage(fit)), call. = FALSE)
      return(data.frame(
        n = NA_real_, cost = NA_real_, objective_value = NA_real_,
        cv = NA_real_, n_int = NA_real_, cost_int = NA_real_,
        .binding = NA, .feasible = FALSE
      ))
    }
    data.frame(
      n = fit$n,
      cost = fit$params$achieved$cost,
      objective_value = fit$objective_value,
      cv = sqrt(fit$objective_value),
      n_int = as.numeric(fit$operational$n),
      cost_int = fit$operational$cost,
      .binding = isTRUE(fit$optimization$budget_binding),
      .feasible = TRUE
    )
  })
  out <- cbind(newdata, do.call(rbind, rows))
  rownames(out) <- NULL
  out
}
