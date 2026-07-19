#' Multi-Indicator Sample Size
#'
#' Compute the sample size that satisfies precision requirements for
#' multiple survey indicators simultaneously under a simple sampling design.
#' Optional domain columns support separate requirements by subpopulation.
#'
#' @param targets For the default method: a data frame where **each row
#'   is one survey indicator** you want to measure. For example,
#'   a prevalence (proportion) or a population mean. Surveys typically
#'   track several indicators simultaneously and the sample must be large
#'   enough for the most demanding one. `n_multi()` finds that size.
#'
#'   See the Details section for the full column reference. At minimum,
#'   each row needs:
#'   \itemize{
#'     \item **What to measure**: `p` for a proportion (e.g. 0.30 for
#'       30\% stunting) **or** `var` for a continuous variable's
#'       population variance. Each row must use exactly one.
#'     \item **How precise**: `moe` (margin of error) **or**
#'       `cv` (coefficient of variation). Each row must
#'       specify exactly one.
#'   }
#'
#'   For `svyplan_prec` objects: a precision result from [prec_multi()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param domains Character vector of column names in `targets` to treat
#'   as domain variables, or `NULL` (default) for no domains. All names
#'   must exist in `targets`. When specified, sizing runs
#'   independently for each domain combination.
#' @param min_n Numeric scalar or `NULL` (default). Minimum total sample
#'   size per domain. It applies only when domains are present. Per-domain
#'   sample sizes are floored to `min_n`.
#' @param prop_method Proportion CI method, one of `"wald"`
#'   (default), `"wilson"`, or `"logodds"`. This is passed to [n_prop()]
#'   for proportion rows and ignored for mean rows.
#'   An optional `prop_method` column in `targets` overrides this default
#'   on a per-row basis.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object. The output class is the same with or without
#'   domains.
#'
#'   **Without domains**, the object contains:
#'   \describe{
#'     \item{`n`}{The sample size required by the binding indicator.}
#'     \item{`detail`}{Per-indicator sample-size results.}
#'     \item{`binding`}{Name or index of the binding (most demanding) indicator.}
#'     \item{`targets`}{The input targets data frame.}
#'   }
#'
#'   **With domains**, the object additionally contains:
#'   \describe{
#'     \item{`n`}{The largest sample size required across domains.}
#'     \item{`domains`}{Data frame with one row per domain, including
#'       domain variables, `.n`, and `.binding`.}
#'   }
#'
#' @details
#' ## Building the targets data frame
#'
#' Each row of `targets` represents one survey indicator. The two key
#' decisions per row are:
#'
#' 1. **Type of indicator**: is it a proportion (binary variable like
#'    "stunted yes/no") or a mean (continuous variable like "household
#'    expenditure")? This determines whether you fill the `p` or `var`
#'    column.
#' 2. **Precision target**: do you want an absolute margin of error
#'    (`moe`, e.g. +/- 5 percentage points) or a relative
#'    coefficient of variation (`cv`, e.g. 10 percent relative error)?
#'
#' A minimal example for three health indicators:
#'
#' ```
#' targets <- data.frame(
#'   name = c("stunting", "vaccination", "expenditure"),
#'   p    = c(0.30, 0.70, NA),
#'   var  = c(NA, NA, 2500),
#'   moe  = c(0.05, 0.05, 10)
#' )
#' ```
#'
#' Rows with `p` are treated as proportions, whereas rows with `var` (and
#' `p = NA`) are treated as means. You cannot have both `p` and `var` non-`NA`
#' in the same row.
#'
#' ## Column reference
#'
#' \describe{
#'   \item{`name`}{Indicator label (optional). If omitted, row numbers
#'     are used in output.}
#'   \item{`p`}{Expected proportion, in (0, 1). Use this for binary
#'     indicators such as prevalences or coverage rates. The value is
#'     your best prior guess (e.g. from a previous survey or literature).
#'     One of `p` or `var` per row.}
#'   \item{`var`}{Population variance of a continuous indicator. Use
#'     this for means (e.g. income, expenditure, weight). One of `p`
#'     or `var` per row.}
#'   \item{`mu`}{Population mean (positive). It is required when `var` is
#'     used with `cv` because CV = SE / mean.}
#'   \item{`moe`}{Margin of error, the half-width of the confidence
#'     interval you want. For proportions, this is on the probability
#'     scale (e.g. 0.05 for +/- 5 percentage points). For means,
#'     it is in the same units as the variable (e.g. 10 dollars).}
#'   \item{`cv`}{Target coefficient of variation (relative standard
#'     error). For example, 0.10 means the SE should be at most 10\%
#'     of the estimate.}
#'   \item{`alpha`}{Significance level for the confidence interval
#'     (default 0.05, giving a 95 percent CI).}
#'   \item{`deff`}{Design effect multiplier (default 1). Set > 1 to inflate
#'     the sample size for complex
#'     designs (e.g. 1.5 for a cluster design).}
#'   \item{`N`}{Population size (default `Inf`).
#'     A finite value applies a finite population correction, reducing
#'     the required sample size.}
#'   \item{`prop_method`}{Proportion CI method:
#'     `"wald"` (default), `"wilson"`, or `"logodds"`. `"wilson"` is
#'     recommended for rare proportions (below 0.1 or above 0.9).
#'     It is used only for rows with `p`.}
#'   \item{`rel_var`}{Unit relvariance. If omitted, derived
#'     automatically from `p` (as `(1 - p) / p`) or from
#'     `var` / `mu^2`.}
#'   \item{`resp_rate`}{Expected response rate, in (0, 1\]. Default 1
#'     (no adjustment). A value of 0.90 inflates the sample size by
#'     `1 / 0.90` to compensate for 10 percent non-response.}
#' }
#'
#' Domain columns are specified via the `domains` parameter. When domains
#' are present, sizing runs independently for each domain combination.
#'
#' `n_multi()` computes sample size per indicator by delegating proportion
#' rows to [n_prop()] and mean rows to [n_mean()],
#' then takes the maximum per domain. Use `prop_method` or a
#' `targets$prop_method` column to choose `"wald"`, `"wilson"`, or
#' `"logodds"` for proportion rows.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @seealso [n_prop()] and [n_mean()] for single-indicator sizing,
#'   [n_multi_cluster()] for multistage cluster designs, and [prec_multi()]
#'   for the inverse.
#'
#' @examples
#' # Simple mode: three indicators, take the max
#' targets <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   moe  = c(0.05, 0.05, 0.03)
#' )
#' n_multi(targets)
#'
#' # MICS/DHS-style: specify precision as a relative margin of error (RME).
#' # RME = moe / p, so convert with moe = RME * p before calling n_multi().
#' rme <- 0.12
#' targets_rme <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   deff = c(2.0, 1.5, 2.5)
#' )
#' targets_rme$moe <- rme * targets_rme$p
#' n_multi(targets_rme)
#'
#' # Rare proportion: use Wilson globally in simple mode
#' n_multi(targets[3, , drop = FALSE], prop_method = "wilson")
#'
#' # Per-row proportion methods in a mixed target table
#' targets_mixed <- data.frame(
#'   name = c("rare_prop", "mean_ind"),
#'   p = c(0.05, NA),
#'   var = c(NA, 100),
#'   moe = c(0.02, 2),
#'   prop_method = c("wilson", NA)
#' )
#' n_multi(targets_mixed)
#'
#' # Simple mode with domains
#' targets_dom <- data.frame(
#'   name   = rep(c("stunting", "anemia"), each = 2),
#'   p      = c(0.30, 0.25, 0.10, 0.15),
#'   moe    = c(0.05, 0.05, 0.03, 0.03),
#'   region = rep(c("North", "South"), 2)
#' )
#' n_multi(targets_dom, domains = "region")
#'
#' # Two-stage CV mode
#' targets_cl <- data.frame(
#'   name   = c("stunting", "anemia"),
#'   p      = c(0.30, 0.10),
#'   cv     = c(0.10, 0.15),
#'   delta_psu = c(0.02, 0.05)
#' )
#' n_multi_cluster(targets_cl, stage_cost = c(500, 50))
#'
#' # Two-stage with MOE (converted to CV internally)
#' targets_moe <- data.frame(
#'   name   = c("stunting", "anemia"),
#'   p      = c(0.30, 0.10),
#'   moe    = c(0.05, 0.03),
#'   delta_psu = c(0.02, 0.05)
#' )
#' n_multi_cluster(targets_moe, stage_cost = c(500, 50))
#'
#' # Joint budget allocation across domains
#' targets_jnt <- data.frame(
#'   name   = rep(c("stunting", "anemia"), each = 2),
#'   p      = c(0.30, 0.25, 0.10, 0.15),
#'   cv     = c(0.10, 0.10, 0.15, 0.15),
#'   delta_psu = c(0.02, 0.03, 0.05, 0.04),
#'   region = rep(c("Urban", "Rural"), 2)
#' )
#' n_multi_cluster(
#'   targets_jnt,
#'   stage_cost = c(500, 50),
#'   domains = "region",
#'   budget = 100000,
#'   joint = TRUE
#' )
#'
#' @export
n_multi <- function(targets, ...) {
  if (!missing(targets)) {
    .res <- .dispatch_plan(targets, "targets", n_multi.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_multi")
}

#' @rdname n_multi
#' @export
n_multi.default <- function(
  targets,
  ...,
  domains = NULL,
  min_n = NULL,
  prop_method = c("wald", "wilson", "logodds"),
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_multi.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_multi.default, c(.plan, list(...))))
  .check_multi_split_args(list(...), "n_multi_cluster()")
  .check_unused_dots(...)
  if (!is.data.frame(targets) || nrow(targets) == 0L) {
    stop("'targets' must be a non-empty data frame", call. = FALSE)
  }
  if (!is.null(min_n)) {
    if (
      !is.numeric(min_n) || length(min_n) != 1L || is.na(min_n) || min_n <= 0
    ) {
      stop("'min_n' must be a positive numeric scalar", call. = FALSE)
    }
  }
  if (missing(prop_method)) {
    prop_method <- prop_method[[1L]]
  }
  if (
    !is.character(prop_method) ||
      length(prop_method) != 1L ||
      is.na(prop_method) ||
      !prop_method %in% c("wald", "wilson", "logodds")
  ) {
    stop(
      "'prop_method' must be one of 'wald', 'wilson', or 'logodds'",
      call. = FALSE
    )
  }

  info <- .validate_targets(targets, FALSE, domains = domains)
  targets <- .fill_defaults(targets, FALSE, prop_method = prop_method)

  rv_final <- targets$rel_var[!is.na(targets$rel_var)]
  if (
    length(rv_final) > 0L &&
      (any(rv_final <= 0) || any(!is.finite(rv_final)))
  ) {
    stop("'rel_var' values must be positive and finite", call. = FALSE)
  }
  domain_cols <- info$domain_cols
  mode <- if ("moe" %in% names(targets) && any(!is.na(targets$moe))) {
    "moe"
  } else {
    "cv"
  }

  if (length(domain_cols) == 0L) {
    .n_multi_simple(targets, domain_cols = domain_cols, mode = mode,
                    prop_method = prop_method)
  } else {
    .n_multi_domains(
      targets,
      stage_cost = NULL,
      budget = NULL,
      n_psu = NULL,
      psu_size = NULL,
      ssu_size = NULL,
      domain_cols,
      multistage = FALSE,
      joint = FALSE,
      min_n,
      fixed_cost = 0,
      mode = mode,
      prop_method = prop_method
    )
  }
}

#' Multi-Indicator Sample Size for Cluster Designs
#'
#' Compute a two- or three-stage cluster allocation that satisfies precision
#' requirements for several survey indicators. Domain-level planning and a
#' shared budget across domains are supported.
#'
#' @param targets For the default method, a non-empty data frame with one row
#'   per indicator. Each row requires `p` or `var`, a `cv` or `moe` target,
#'   and `delta_psu`. Three-stage designs also require `delta_ssu`. For the
#'   `svyplan_prec` method, a result from [prec_multi_cluster()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param stage_cost Numeric vector of per-stage costs with length 2 or 3.
#' @param domains Optional character vector naming domain columns in
#'   `targets`. The function solves each domain independently unless `joint`
#'   is `TRUE` in budget mode.
#' @param budget Optional total budget. Supply precision targets or a budget,
#'   according to the target schema described in Details.
#' @param n_psu Optional fixed stage-1 sample size.
#' @param psu_size Optional fixed stage-2 sample size per PSU.
#' @param ssu_size Optional fixed stage-3 sample size per SSU. This is valid
#'   only for three-stage designs.
#' @param joint If `TRUE`, split one budget across domains to minimize the
#'   worst precision ratio. This applies only when domains and `budget` are
#'   supplied.
#' @param min_n Optional positive minimum total sample size per domain. In
#'   joint budget mode it is a constraint. In independent domain mode,
#'   domains below the floor produce a warning.
#' @param fixed_cost Non-negative fixed overhead cost. The default is 0.
#' @param plan Optional [svyplan()] profile providing `stage_cost` and other
#'   applicable defaults.
#'
#' @return A `svyplan_cluster` object. The output class does not depend on
#'   which optional arguments are supplied.
#'
#' @details
#' Margin-of-error targets are converted to CV before optimization. For each
#' candidate allocation, the required stage-1 size is the maximum across all
#' indicators. The solver minimizes total cost for precision targets or the
#' worst precision ratio under a fixed budget.
#'
#' Homogeneity values numerically close to 0 or 1 are rejected because they
#' make the analytical cluster optimum degenerate. The result includes an
#' integer `$operational` allocation that preserves the applicable precision
#' or budget constraint. See [n_multi()] for shared indicator columns and
#' [n_cluster()] for the cluster cost model.
#'
#' @seealso [n_multi()] for simple designs and [prec_multi_cluster()] for the
#'   inverse calculation.
#'
#' @examples
#' targets <- data.frame(
#'   name = c("stunting", "anemia"),
#'   p = c(0.30, 0.10),
#'   cv = c(0.10, 0.15),
#'   delta_psu = c(0.02, 0.05)
#' )
#' n_multi_cluster(targets, stage_cost = c(500, 50))
#'
#' @export
n_multi_cluster <- function(targets, ...) {
  if (!missing(targets)) {
    .res <- .dispatch_plan(targets, "targets", n_multi_cluster.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_multi_cluster")
}

#' @rdname n_multi_cluster
#' @export
n_multi_cluster.default <- function(
  targets,
  ...,
  stage_cost = NULL,
  domains = NULL,
  budget = NULL,
  n_psu = NULL,
  psu_size = NULL,
  ssu_size = NULL,
  joint = FALSE,
  min_n = NULL,
  fixed_cost = 0,
  plan = NULL
) {
  .plan <- .merge_plan_args(
    plan,
    n_multi_cluster.default,
    match.call(),
    environment()
  )
  if (!is.null(.plan)) {
    return(do.call(n_multi_cluster.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)

  if (!is.data.frame(targets) || nrow(targets) == 0L) {
    stop("'targets' must be a non-empty data frame", call. = FALSE)
  }
  if (is.null(stage_cost)) {
    stop("'stage_cost' is required (directly or via plan)", call. = FALSE)
  }
  if (!is.logical(joint) || length(joint) != 1L || is.na(joint)) {
    stop("'joint' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(min_n) &&
      (!is.numeric(min_n) || length(min_n) != 1L || is.na(min_n) ||
       min_n <= 0)) {
    stop("'min_n' must be a positive numeric scalar", call. = FALSE)
  }

  check_stage_cost(stage_cost)
  stage_cost <- .reorder_stage_cost(stage_cost)
  stages <- length(stage_cost)
  if (!is.null(budget)) check_scalar(budget, "budget")
  if (!is.null(n_psu)) check_scalar(n_psu, "n_psu")
  if (!is.null(psu_size)) check_scalar(psu_size, "psu_size")
  if (!is.null(ssu_size)) check_scalar(ssu_size, "ssu_size")
  if (stages == 2L && !is.null(ssu_size)) {
    stop("'ssu_size' is not applicable for 2-stage designs", call. = FALSE)
  }
  n_fixed <- sum(!is.null(n_psu), !is.null(psu_size), !is.null(ssu_size))
  if (n_fixed >= stages) {
    stop("cannot fix all stages; use prec_multi_cluster() instead",
         call. = FALSE)
  }
  check_fixed_cost(fixed_cost, budget)

  info <- .validate_targets(
    targets,
    TRUE,
    domains = domains,
    stages = stages,
    context = "n_multi_cluster()"
  )
  targets <- .fill_defaults(targets, TRUE)
  targets <- .convert_moe_to_cv(targets)

  rv_final <- targets$rel_var[!is.na(targets$rel_var)]
  if (length(rv_final) > 0L &&
      (any(rv_final <= 0) || any(!is.finite(rv_final)))) {
    stop("'rel_var' values must be positive and finite", call. = FALSE)
  }
  if (any(targets$k_psu <= 0) || any(!is.finite(targets$k_psu))) {
    stop("'k_psu' values must be positive and finite", call. = FALSE)
  }
  if (any(targets$k_ssu <= 0) || any(!is.finite(targets$k_ssu))) {
    stop("'k_ssu' values must be positive and finite", call. = FALSE)
  }

  domain_cols <- info$domain_cols
  mode <- if (!is.null(budget)) {
    "budget"
  } else if ("moe" %in% names(targets) && any(!is.na(targets$moe))) {
    "moe"
  } else {
    "cv"
  }

  if (length(domain_cols) == 0L) {
    .n_multi_cluster(
      targets,
      stage_cost,
      budget,
      n_psu,
      psu_size,
      ssu_size,
      fixed_cost,
      domain_cols = domain_cols,
      mode = mode
    )
  } else {
    .n_multi_domains(
      targets,
      stage_cost,
      budget,
      n_psu,
      psu_size,
      ssu_size,
      domain_cols,
      multistage = TRUE,
      joint,
      min_n,
      fixed_cost,
      mode = mode
    )
  }
}

#' Report arguments that moved to the cluster-specific API
#' @keywords internal
#' @noRd
.check_multi_split_args <- function(dots, replacement) {
  moved <- intersect(
    names(dots) %||% character(0),
    c("stage_cost", "budget", "n_psu", "psu_size", "ssu_size", "joint",
      "fixed_cost")
  )
  if (length(moved) > 0L) {
    stop(
      sprintf(
        "cluster argument%s %s moved to %s",
        if (length(moved) > 1L) "s" else "",
        paste(sQuote(moved), collapse = ", "),
        replacement
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Validate targets data frame
#' @return List with `indicator_type` (per-row "p" or "var") and `domain_cols`.
#' @keywords internal
#' @noRd
.validate_targets <- function(targets, multistage, domains = NULL,
                              stages = NULL,
                              context = "n_multi_cluster()") {
  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(targets))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in targets: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)
  if (!has_p && !has_var) {
    stop("'targets' must contain 'p' or 'var' column", call. = FALSE)
  }

  has_moe <- "moe" %in% names(targets)
  has_cv <- "cv" %in% names(targets)
  if (!has_moe && !has_cv) {
    stop("'targets' must contain 'moe' or 'cv' column", call. = FALSE)
  }

  if (has_p) {
    p_vals <- targets$p[!is.na(targets$p)]
    if (any(p_vals <= 0 | p_vals >= 1)) {
      stop("all 'p' values must be in (0, 1)", call. = FALSE)
    }
  }

  if (has_var) {
    var_vals <- targets$var[!is.na(targets$var)]
    if (any(var_vals <= 0) || any(!is.finite(var_vals))) {
      stop("all 'var' values must be positive and finite", call. = FALSE)
    }
  }

  # Each row needs at least one of p or var (non-NA)
  has_indicator <- rep(FALSE, nrow(targets))
  if (has_p) {
    has_indicator <- has_indicator | !is.na(targets$p)
  }
  if (has_var) {
    has_indicator <- has_indicator | !is.na(targets$var)
  }
  if (any(!has_indicator)) {
    stop(
      sprintf(
        "row(s) %s must have a non-NA 'p' or 'var' value",
        paste(which(!has_indicator), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (has_p && has_var) {
    both_set <- !is.na(targets$p) & !is.na(targets$var)
    if (any(both_set)) {
      stop("each row must have only one of 'p' or 'var'", call. = FALSE)
    }
  }

  if (has_moe && has_cv) {
    both_na <- is.na(targets$moe) & is.na(targets$cv)
    both_set <- !is.na(targets$moe) & !is.na(targets$cv)
    if (any(both_na)) {
      stop("each row must have either 'moe' or 'cv' specified", call. = FALSE)
    }
    if (any(both_set)) {
      stop("each row must have only one of 'moe' or 'cv'", call. = FALSE)
    }
  }

  if (has_moe) {
    moe_vals <- targets$moe[!is.na(targets$moe)]
    if (any(moe_vals <= 0) || any(!is.finite(moe_vals))) {
      stop("'moe' values must be positive and finite", call. = FALSE)
    }
  }

  if ("prop_method" %in% names(targets)) {
    method_vals <- targets$prop_method[!is.na(targets$prop_method)]
    bad_methods <- !method_vals %in% c("wald", "wilson", "logodds")
    if (any(bad_methods)) {
      stop(
        "'prop_method' values must be one of 'wald', 'wilson', or 'logodds'",
        call. = FALSE
      )
    }
  }

  if (multistage) {
    if (!has_cv && !has_moe) {
      stop("multistage mode requires 'cv' or 'moe' column in targets",
           call. = FALSE)
    }
    if (has_moe && any(!is.na(targets$moe))) {
      moe_rows <- !is.na(targets$moe)
      var_moe <- has_var & moe_rows & !is.na(targets$var)
      if (any(var_moe)) {
        if (!"mu" %in% names(targets) || any(is.na(targets$mu[var_moe]))) {
          stop(
            "'mu' is required to convert 'moe' to 'cv' for mean indicators ",
            "in multistage mode",
            call. = FALSE
          )
        }
      }
    }
    if (!"delta_psu" %in% names(targets)) {
      stop(
        "multistage mode requires 'delta_psu' column in targets",
        call. = FALSE
      )
    }
    if (isTRUE(stages == 3L)) {
      if (!"delta_ssu" %in% names(targets)) {
        stop(
          "3-stage mode requires a 'delta_ssu' column in targets",
          call. = FALSE
        )
      }
      if (anyNA(targets$delta_ssu) || any(!is.finite(targets$delta_ssu))) {
        stop("'delta_ssu' must contain finite non-missing values",
             call. = FALSE)
      }
    }
    if (has_cv) {
      cv_vals <- targets$cv[!is.na(targets$cv)]
      if (length(cv_vals) > 0L && (any(cv_vals <= 0) || any(!is.finite(cv_vals)))) {
        stop("'cv' values must be positive and finite", call. = FALSE)
      }
    }
    d1_vals <- targets$delta_psu[!is.na(targets$delta_psu)]
    if (any(d1_vals < 0 | d1_vals > 1)) {
      stop("'delta_psu' values must be in [0, 1]", call. = FALSE)
    }
    .check_cluster_delta_open(d1_vals, context = context)
    if ("delta_ssu" %in% names(targets)) {
      d2_vals <- targets$delta_ssu[!is.na(targets$delta_ssu)]
      if (any(d2_vals < 0 | d2_vals > 1)) {
        stop("'delta_ssu' values must be in [0, 1]", call. = FALSE)
      }
      .check_cluster_delta_open(d2_vals, context = context)
    }
  }

  if (has_var && has_cv) {
    needs_mu <- !is.na(targets$var) & !is.na(targets$cv)
    if (any(needs_mu)) {
      if (!"mu" %in% names(targets) || any(is.na(targets$mu[needs_mu]))) {
        stop(
          "'mu' is required when 'var' and 'cv' are specified",
          call. = FALSE
        )
      }
    }
  }

  if ("mu" %in% names(targets)) {
    mu_vals <- targets$mu[!is.na(targets$mu)]
    if (any(mu_vals <= 0) || any(!is.finite(mu_vals))) {
      stop("'mu' values must be positive and finite", call. = FALSE)
    }
  }

  domain_cols <- domains %||% character(0)

  if ("rel_var" %in% names(targets)) {
    rv_vals <- targets$rel_var[!is.na(targets$rel_var)]
    if (any(rv_vals <= 0) || any(!is.finite(rv_vals))) {
      stop("'rel_var' values must be positive and finite", call. = FALSE)
    }
  }
  if (multistage) {
    if ("k_psu" %in% names(targets)) {
      k_psu_vals <- targets$k_psu[!is.na(targets$k_psu)]
      if (any(k_psu_vals <= 0) || any(!is.finite(k_psu_vals))) {
        stop("'k_psu' values must be positive and finite", call. = FALSE)
      }
    }
    if ("k_ssu" %in% names(targets)) {
      k_ssu_vals <- targets$k_ssu[!is.na(targets$k_ssu)]
      if (any(k_ssu_vals <= 0) || any(!is.finite(k_ssu_vals))) {
        stop("'k_ssu' values must be positive and finite", call. = FALSE)
      }
    }
  }

  .validate_common_columns(targets)

  list(domain_cols = domain_cols)
}

#' Fill default values in targets
#' @keywords internal
#' @noRd
.fill_defaults <- function(targets, multistage, prop_method = "wald") {
  if (!"alpha" %in% names(targets)) {
    targets$alpha <- 0.05
  } else {
    targets$alpha[is.na(targets$alpha)] <- 0.05
  }

  if (!multistage) {
    if (!"deff" %in% names(targets)) {
      targets$deff <- 1
    } else {
      targets$deff[is.na(targets$deff)] <- 1
    }

    if (!"N" %in% names(targets)) {
      targets$N <- Inf
    } else {
      targets$N[is.na(targets$N)] <- Inf
    }
  }

  if (multistage) {
    if (!"k_psu" %in% names(targets)) {
      targets$k_psu <- 1
    } else {
      targets$k_psu[is.na(targets$k_psu)] <- 1
    }
    if (!"k_ssu" %in% names(targets)) {
      targets$k_ssu <- 1
    } else {
      targets$k_ssu[is.na(targets$k_ssu)] <- 1
    }
  }

  if (!"resp_rate" %in% names(targets)) {
    targets$resp_rate <- 1
  } else {
    targets$resp_rate[is.na(targets$resp_rate)] <- 1
  }

  if (!"prop_method" %in% names(targets)) {
    targets$prop_method <- prop_method
  } else {
    targets$prop_method[is.na(targets$prop_method)] <- prop_method
  }

  if (!"rel_var" %in% names(targets)) {
    targets$rel_var <- NA_real_
  }
  targets$rel_var <- .derive_rel_var(targets, require_all = multistage)

  targets
}


#' Derive unit relvariance from p or var/mu
#' @keywords internal
#' @noRd
.derive_rel_var <- function(targets, require_all = FALSE) {
  rv <- targets$rel_var
  needs <- is.na(rv)

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)
  has_mu <- "mu" %in% names(targets)

  for (i in which(needs)) {
    if (has_p && !is.na(targets$p[i])) {
      rv[i] <- (1 - targets$p[i]) / targets$p[i]
    } else if (has_var && !is.na(targets$var[i])) {
      if (has_mu && !is.na(targets$mu[i])) {
        rv[i] <- targets$var[i] / targets$mu[i]^2
      } else if (require_all) {
        stop(
          sprintf("row %d: 'mu' is required to derive 'rel_var' from 'var'", i),
          call. = FALSE
        )
      }
    }
  }

  if (require_all && anyNA(rv)) {
    stop(
      sprintf(
        "row(s) %s: could not derive 'rel_var' - provide 'p', or 'var'+'mu', or 'rel_var' directly",
        paste(which(is.na(rv)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  rv
}

#' Convert moe to cv in multistage targets
#'
#' For proportion rows: cv = moe / (z * p).
#' For mean rows: cv = moe / (z * mu).
#' Rows that already have cv are left unchanged.
#' @keywords internal
#' @noRd
.convert_moe_to_cv <- function(targets) {
  has_moe <- "moe" %in% names(targets)
  if (!has_moe) return(targets)

  moe_rows <- !is.na(targets$moe)
  if (!any(moe_rows)) return(targets)

  if (!"cv" %in% names(targets)) {
    targets$cv <- NA_real_
  }

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)

  for (i in which(moe_rows)) {
    z <- qnorm(1 - targets$alpha[i] / 2)
    if (has_p && !is.na(targets$p[i])) {
      targets$cv[i] <- targets$moe[i] / (z * targets$p[i])
    } else if (has_var && !is.na(targets$var[i])) {
      targets$cv[i] <- targets$moe[i] / (z * targets$mu[i])
    }
  }

  targets
}

#' Simple mode: compute n per indicator, take max
#' @keywords internal
#' @noRd
.n_multi_simple <- function(targets, domain_cols = character(0), mode = "moe",
                           prop_method = "wald") {
  simple <- .compute_simple_n(targets)
  n_vec <- simple$n
  cv_target_vec <- simple$cv_target

  idx <- which.max(n_vec)
  n_max <- n_vec[idx]

  labels <- if ("name" %in% names(targets)) {
    targets$name
  } else {
    seq_len(nrow(targets))
  }
  binding_name <- labels[idx]

  has_p <- "p" %in% names(targets)
  has_mu <- "mu" %in% names(targets)
  cv_achieved_vec <- vapply(seq_len(nrow(targets)), function(i) {
    prec <- suppressWarnings(
      if (has_p && !is.na(targets$p[i])) {
        .prec_engine_prop(targets$p[i], n_max, targets$alpha[i],
                          targets$N[i], targets$deff[i],
                          targets$resp_rate[i], targets$prop_method[i])
      } else {
        .prec_engine_mean(targets$var[i],
                          if (has_mu) targets$mu[i] else NULL,
                          n_max, targets$alpha[i], targets$N[i],
                          targets$deff[i], targets$resp_rate[i])
      }
    )
    prec$cv
  }, numeric(1L))

  detail <- data.frame(
    name = labels,
    .n = n_vec,
    .cv_target = cv_target_vec,
    .cv_achieved = cv_achieved_vec,
    .binding = seq_len(nrow(targets)) == idx
  )

  .new_svyplan_n(
    n = n_max,
    type = "multi",
    params = list(domain_cols = domain_cols, mode = mode,
                  prop_method = prop_method),
    targets = targets,
    detail = detail,
    binding = binding_name
  )
}

#' Compute per-indicator n by delegating to the single-indicator engines
#'
#' Returns list with `n` (per-indicator sample size) and `cv_target` (the CV
#' each indicator achieves at its own n; equals the target CV for cv-mode and
#' the CV implied by the target MOE for moe-mode).
#' @keywords internal
#' @noRd
.compute_simple_n <- function(targets) {
  nr <- nrow(targets)
  n_vec <- numeric(nr)
  cv_vec <- numeric(nr)

  has_p <- "p" %in% names(targets)
  has_moe <- "moe" %in% names(targets)
  has_mu <- "mu" %in% names(targets)

  for (i in seq_len(nr)) {
    is_prop <- has_p && !is.na(targets$p[i])
    use_moe <- has_moe && !is.na(targets$moe[i])

    if (is_prop) {
      res_i <- n_prop.default(
        p = targets$p[i],
        moe = if (use_moe) targets$moe[i] else NULL,
        cv = if (use_moe) NULL else targets$cv[i],
        alpha = targets$alpha[i],
        N = targets$N[i],
        deff = targets$deff[i],
        resp_rate = targets$resp_rate[i],
        method = targets$prop_method[i]
      )
    } else {
      res_i <- n_mean.default(
        var = targets$var[i],
        mu = if (has_mu && !is.na(targets$mu[i])) targets$mu[i] else NULL,
        moe = if (use_moe) targets$moe[i] else NULL,
        cv = if (use_moe) NULL else targets$cv[i],
        alpha = targets$alpha[i],
        N = targets$N[i],
        deff = targets$deff[i],
        resp_rate = targets$resp_rate[i]
      )
    }
    n_vec[i] <- res_i$n
    cv_vec[i] <- if (!is.null(res_i$cv)) res_i$cv else NA_real_
  }

  list(n = n_vec, cv_target = cv_vec)
}

#' Multistage cluster mode dispatcher
#' @keywords internal
#' @noRd
.n_multi_cluster <- function(targets, stage_cost, budget, n_psu,
                            psu_size = NULL, ssu_size = NULL, fixed_cost = 0,
                            domain_cols = character(0), mode = "cv",
                            prop_method = "wald") {
  stages <- length(stage_cost)
  if (stages == 2L) {
    .n_multi_2stage(targets, stage_cost, budget, n_psu, psu_size, fixed_cost,
                    domain_cols = domain_cols, mode = mode,
                    prop_method = prop_method)
  } else {
    .n_multi_3stage(targets, stage_cost, budget, n_psu, psu_size, ssu_size,
                    fixed_cost, domain_cols = domain_cols, mode = mode,
                    prop_method = prop_method)
  }
}

#' Candidate whole stage sizes around a continuous optimum
#' @keywords internal
#' @noRd
.multi_stage_candidates <- function(fixed, continuous, upper = Inf,
                                    base_limit = 400L) {
  if (!is.null(fixed)) {
    return(max(1L, as.integer(round(fixed))))
  }
  upper_i <- if (is.finite(upper)) {
    max(1L, as.integer(floor(upper + 1e-9)))
  } else {
    .Machine$integer.max
  }
  base <- seq_len(min(base_limit, upper_i))
  centre <- max(1L, as.integer(round(continuous)))
  local <- seq.int(max(1L, centre - 20L), min(upper_i, centre + 20L))
  unique(c(base, local))
}

#' Constraint-preserving whole-unit design for two-stage n_multi
#' @keywords internal
#' @noRd
.op_multi_2stage <- function(n1_required, cv_fn, cv_t, stage_cost, budget,
                             n_psu, psu_size, fixed_cost, cont_m) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  variable_budget <- if (is.null(budget)) NULL else budget - fixed_cost

  upper_m <- if (is.null(variable_budget)) {
    Inf
  } else if (!is.null(n_psu)) {
    (variable_budget / max(1L, as.integer(round(n_psu))) - C1) / C2
  } else {
    (variable_budget - C1) / C2
  }
  m_cand <- .multi_stage_candidates(
    psu_size, cont_m, upper = upper_m, base_limit = 100000L
  )

  if (!is.null(budget)) {
    a <- if (!is.null(n_psu)) {
      rep(max(1L, as.integer(round(n_psu))), length(m_cand))
    } else {
      as.integer(floor(variable_budget / (C1 + C2 * m_cand)))
    }
    cost <- a * (C1 + C2 * m_cand)
    keep <- a >= 1L & cost <= variable_budget + 1e-8
    if (!any(keep)) {
      stop("no whole-unit n_multi design fits the budget", call. = FALSE)
    }
    a <- a[keep]
    m_cand <- m_cand[keep]
    cost <- cost[keep]
    ratios <- vapply(seq_along(a), function(i) {
      max(cv_fn(a[i], m_cand[i]) / cv_t)
    }, numeric(1L))
    j <- which.min(ratios + 1e-12 * m_cand)
  } else {
    need <- vapply(m_cand, function(m) max(n1_required(m)), numeric(1L))
    if (!is.null(n_psu)) {
      a <- rep(max(1L, as.integer(round(n_psu))), length(m_cand))
      keep <- a + 1e-9 >= need
      if (!any(keep)) {
        stop("target CV is not achievable with the fixed whole PSU count",
             call. = FALSE)
      }
      a <- a[keep]
      m_cand <- m_cand[keep]
    } else {
      a <- pmax(1L, as.integer(ceiling(need - 1e-9)))
    }
    cost <- a * (C1 + C2 * m_cand)
    j <- which.min(cost + 1e-9 * a * m_cand + 1e-12 * m_cand)
  }

  a_best <- a[j]
  m_best <- m_cand[j]
  cvs <- cv_fn(a_best, m_best)
  list(
    n = c(n_psu = a_best, psu_size = m_best),
    total_n = a_best * m_best,
    cost = fixed_cost + a_best * (C1 + C2 * m_best),
    cv = max(cvs),
    cv_by_target = cvs
  )
}

#' Constraint-preserving whole-unit design for three-stage n_multi
#' @keywords internal
#' @noRd
.op_multi_3stage <- function(n1_required, cv_fn, cv_t, stage_cost, budget,
                             n_psu, psu_size, ssu_size, fixed_cost,
                             cont_m, cont_q) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  variable_budget <- if (is.null(budget)) NULL else budget - fixed_cost
  a_fixed <- if (is.null(n_psu)) NULL else max(1L, as.integer(round(n_psu)))
  a_min <- a_fixed %||% 1L

  if (is.null(variable_budget)) {
    upper_m <- upper_q <- Inf
  } else {
    per <- variable_budget / a_min
    upper_m <- (per - C1) / (C2 + C3)
    upper_q <- (per - C1 - C2) / C3
  }
  m_cand <- .multi_stage_candidates(psu_size, cont_m, upper_m)
  q_cand <- .multi_stage_candidates(ssu_size, cont_q, upper_q)
  grid <- expand.grid(m = m_cand, q = q_cand)
  per_psu <- C1 + C2 * grid$m + C3 * grid$m * grid$q

  if (!is.null(budget)) {
    a <- if (is.null(a_fixed)) {
      as.integer(floor(variable_budget / per_psu))
    } else {
      rep(a_fixed, nrow(grid))
    }
    cost <- a * per_psu
    keep <- a >= 1L & cost <= variable_budget + 1e-8
    if (!any(keep)) {
      stop("no whole-unit n_multi design fits the budget", call. = FALSE)
    }
    grid <- grid[keep, , drop = FALSE]
    a <- a[keep]
    cost <- cost[keep]
    ratios <- vapply(seq_along(a), function(i) {
      max(cv_fn(a[i], grid$m[i], grid$q[i]) / cv_t)
    }, numeric(1L))
    j <- which.min(ratios + 1e-12 * (grid$m + grid$q))
  } else {
    need <- vapply(seq_len(nrow(grid)), function(i) {
      max(n1_required(grid$m[i], grid$q[i]))
    }, numeric(1L))
    if (is.null(a_fixed)) {
      a <- pmax(1L, as.integer(ceiling(need - 1e-9)))
    } else {
      a <- rep(a_fixed, nrow(grid))
      keep <- a + 1e-9 >= need
      if (!any(keep)) {
        stop("target CV is not achievable with the fixed whole PSU count",
             call. = FALSE)
      }
      grid <- grid[keep, , drop = FALSE]
      per_psu <- per_psu[keep]
      a <- a[keep]
    }
    cost <- a * (C1 + C2 * grid$m + C3 * grid$m * grid$q)
    j <- which.min(cost + 1e-9 * a * grid$m * grid$q +
                     1e-12 * (grid$m + grid$q))
  }

  a_best <- a[j]
  m_best <- as.integer(grid$m[j])
  q_best <- as.integer(grid$q[j])
  cvs <- cv_fn(a_best, m_best, q_best)
  list(
    n = c(n_psu = a_best, psu_size = m_best, ssu_size = q_best),
    total_n = a_best * m_best * q_best,
    cost = fixed_cost + a_best * (C1 + C2 * m_best + C3 * m_best * q_best),
    cv = max(cvs),
    cv_by_target = cvs
  )
}

#' 2-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_2stage <- function(targets, stage_cost, budget, n_psu,
                           psu_size = NULL, fixed_cost = 0,
                           domain_cols = character(0), mode = "cv",
                           prop_method = "wald") {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  nr <- nrow(targets)

  cv_t <- targets$cv
  delta <- targets$delta_psu
  rel_var <- targets$rel_var
  k <- targets$k_psu
  rr <- targets$resp_rate
  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  # n1_actual_j = rel_var_j * k_j * (1 + delta_j*(psu_size-1)) / (psu_size * cv_j^2 * rr_j)
  n1_required <- function(psu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        rel_var[j] *
          k[j] *
          (1 + delta[j] * (psu_size - 1)) /
          (psu_size * cv_t[j]^2 * rr[j])
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, psu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] *
            k[j] /
            (n1 * rr[j] * psu_size) *
            (1 + delta[j] * (psu_size - 1))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(budget)) {
    if (!is.null(psu_size)) {
      psu_size_opt <- psu_size
      n1_vals <- n1_required(psu_size_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost + n1_opt * (C1 + C2 * psu_size_opt)
      binding_idx <- which.max(n1_vals)
    } else {
      cost_fn <- function(psu_size) {
        n1 <- max(n1_required(psu_size))
        n1 * (C1 + C2 * psu_size)
      }

      upper <- max(10, sqrt(C1 / C2 * (1 - delta) / delta))
      opt <- optimize(cost_fn, interval = c(1, upper))
      psu_size_opt <- opt$minimum

      n1_vals <- n1_required(psu_size_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost + n1_opt * (C1 + C2 * psu_size_opt)
      binding_idx <- which.max(n1_vals)
    }

    cv_achieved <- cv_achieved_fn(n1_opt, psu_size_opt)
  } else {
    var_budget <- budget - fixed_cost
    bres <- .eval_2stage_budget(
      cv_t,
      delta,
      rel_var,
      k,
      rr,
      C1,
      C2,
      var_budget,
      n_psu,
      psu_size
    )
    n1_opt <- bres$n1
    psu_size_opt <- bres$psu_size
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- budget
  }

  n_vec <- c(n_psu = n1_opt, psu_size = psu_size_opt)
  total_n <- prod(n_vec)
  operational <- .op_multi_2stage(
    n1_required, cv_achieved_fn, cv_t, stage_cost, budget,
    n_psu, psu_size, fixed_cost, cont_m = psu_size_opt
  )

  n1_per <- n1_required(psu_size_opt)
  n_per <- n1_per * psu_size_opt

  detail <- data.frame(
    name = labels,
    .n = n_per,
    .cv_target = cv_t,
    .cv_achieved = cv_achieved,
    .binding = seq_len(nr) == binding_idx
  )

  params <- list(stage_cost = stage_cost, domain_cols = domain_cols,
                  mode = mode, prop_method = prop_method)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) {
    params$n_psu <- n_psu
  }
  if (!is.null(psu_size)) {
    params$psu_size <- psu_size
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 2L,
    total_n = total_n,
    cv = cv_achieved[binding_idx],
    cost = total_cost,
    params = params,
    targets = targets,
    detail = detail,
    binding = labels[binding_idx],
    operational = operational
  )
}

#' 3-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_3stage <- function(targets, stage_cost, budget, n_psu,
                           psu_size = NULL, ssu_size = NULL, fixed_cost = 0,
                           domain_cols = character(0), mode = "cv",
                           prop_method = "wald") {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  nr <- nrow(targets)

  cv_t <- targets$cv
  delta_psu <- targets$delta_psu
  delta_ssu <- if ("delta_ssu" %in% names(targets)) {
    targets$delta_ssu
  } else {
    rep(0, nr)
  }
  rel_var <- targets$rel_var
  k_psu <- targets$k_psu
  k_ssu <- targets$k_ssu
  rr <- targets$resp_rate
  labels <- if ("name" %in% names(targets)) {
    targets$name
  } else {
    seq_len(nr)
  }

  n1_required <- function(psu_size, ssu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        rel_var[j] /
          (cv_t[j]^2 * psu_size * ssu_size * rr[j]) *
          (k_psu[j] *
            delta_psu[j] *
            psu_size *
            ssu_size +
            k_ssu[j] * (1 + delta_ssu[j] * (ssu_size - 1)))
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, psu_size, ssu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] /
            (n1 * rr[j] * psu_size * ssu_size) *
            (k_psu[j] *
              delta_psu[j] *
              psu_size *
              ssu_size +
              k_ssu[j] * (1 + delta_ssu[j] * (ssu_size - 1)))
        )
      },
      numeric(1L)
    )
  }

  solve_for <- if (!is.null(n_psu) && !is.null(psu_size)) {
    "n3"
  } else if (!is.null(n_psu)) {
    "n2"
  } else {
    "n1"
  }

  n_free <- 3L - sum(!is.null(n_psu), !is.null(psu_size), !is.null(ssu_size))

  if (is.null(budget)) {
    if (n_free == 3L) {
      cost_fn <- function(par) {
        ps <- par[1L]
        ss <- par[2L]
        n1 <- max(n1_required(ps, ss))
        n1 * (C1 + C2 * ps + C3 * ps * ss)
      }

      init_ssu_size <- max(2, sqrt(C2 / C3))
      init_psu_size <- max(2, sqrt(C1 / C2))
      upper_ps <- max(1000, 10 * init_psu_size)
      upper_ss <- max(1000, 10 * init_ssu_size)
      opt <- optim(
        par = c(init_psu_size, init_ssu_size),
        fn = cost_fn,
        method = "L-BFGS-B",
        lower = c(1, 1),
        upper = c(upper_ps, upper_ss)
      )
      if (opt$convergence != 0L) {
        warning(
          "L-BFGS-B did not converge (code ",
          opt$convergence,
          "): ",
          "allocation may be approximate",
          call. = FALSE
        )
      }
      if (opt$par[1L] >= upper_ps * (1 - 1e-6) ||
        opt$par[2L] >= upper_ss * (1 - 1e-6)) {
        warning(
          sprintf(
            "optimal stage size reached the search upper bound (psu_size <= %.0f, ssu_size <= %.0f); result may be unreliable -- review stage costs and target CVs",
            upper_ps,
            upper_ss
          ),
          call. = FALSE
        )
      }

      psu_size_opt <- opt$par[1L]
      ssu_size_opt <- opt$par[2L]
      n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost +
        n1_opt * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
      binding_idx <- which.max(n1_vals)
    } else if (n_free == 2L) {
      if (solve_for == "n2" && is.null(ssu_size)) {
        n1_opt <- n_psu

        cv_floor <- .multistage_cv_floor(
          rel_var, k_psu, delta_psu,
          rr = rr, n_psu = n_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          rel_var, k_psu, delta_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )

        psu_size_required_fn <- function(ss) {
          psu_size_per <- vapply(
            seq_len(nr),
            function(j) {
              denom <- cv_t[j]^2 *
                n_psu *
                rr[j] /
                (rel_var[j] * k_ssu[j]) -
                k_psu[j] * delta_psu[j] / k_ssu[j]
              if (denom <= 0) Inf
              else (1 + delta_ssu[j] * (ss - 1)) / (ss * denom)
            },
            numeric(1L)
          )
          max(psu_size_per)
        }

        cost_fn_fixed <- function(ss) {
          ps <- psu_size_required_fn(ss)
          n_psu * (C1 + C2 * ps + C3 * ps * ss)
        }

        ssu_size_analytic <- vapply(
          seq_len(nr),
          function(j) {
            if (delta_ssu[j] <= 0) 1
            else sqrt((1 - delta_ssu[j]) / delta_ssu[j] * C2 / C3)
          },
          numeric(1L)
        )
        upper_ssu_size <- max(10, 3 * max(ssu_size_analytic))

        opt <- optimize(cost_fn_fixed, interval = c(1, upper_ssu_size))
        ssu_size_opt <- opt$minimum
        psu_size_opt <- psu_size_required_fn(ssu_size_opt)

        if (!is.finite(psu_size_opt) || psu_size_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }

        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
      } else if (!is.null(ssu_size) && is.null(n_psu) && is.null(psu_size)) {
        ssu_size_opt <- ssu_size
        cost_fn_ss <- function(ps) {
          n1 <- max(n1_required(ps, ssu_size_opt))
          n1 * (C1 + C2 * ps + C3 * ps * ssu_size_opt)
        }
        psu_size_analytic <- vapply(
          seq_len(nr),
          function(j) {
            sqrt(
              C1 * k_ssu[j] * (1 + delta_ssu[j] * (ssu_size_opt - 1)) /
                (k_psu[j] * delta_psu[j] * ssu_size_opt *
                   (C2 + C3 * ssu_size_opt))
            )
          },
          numeric(1L)
        )
        upper <- max(10, psu_size_analytic)
        opt <- optimize(cost_fn_ss, interval = c(1, upper))
        psu_size_opt <- opt$minimum
        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        n1_opt <- max(n1_vals)
        total_cost <- fixed_cost +
          n1_opt * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
        binding_idx <- which.max(n1_vals)
      } else if (!is.null(psu_size) && is.null(n_psu) && is.null(ssu_size)) {
        psu_size_opt <- psu_size
        cost_fn_ps <- function(ss) {
          n1 <- max(n1_required(psu_size_opt, ss))
          n1 * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ss)
        }
        ssu_size_analytic <- vapply(
          seq_len(nr),
          function(j) {
            if (delta_ssu[j] <= 0) 1
            else sqrt((1 - delta_ssu[j]) / delta_ssu[j] * C2 / C3)
          },
          numeric(1L)
        )
        upper <- max(10, 3 * max(ssu_size_analytic))
        opt <- optimize(cost_fn_ps, interval = c(1, upper))
        ssu_size_opt <- opt$minimum
        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        n1_opt <- max(n1_vals)
        total_cost <- fixed_cost +
          n1_opt * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
        binding_idx <- which.max(n1_vals)
      } else {
        n1_opt <- n_psu
        ssu_size_opt <- ssu_size
        cv_floor <- .multistage_cv_floor(
          rel_var, k_psu, delta_psu,
          rr = rr, n_psu = n_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          rel_var, k_psu, delta_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )
        psu_size_required_fn2 <- function(j) {
          denom <- cv_t[j]^2 * n_psu * rr[j] / (rel_var[j] * k_ssu[j]) -
            k_psu[j] * delta_psu[j] / k_ssu[j]
          if (denom <= 0) Inf
          else (1 + delta_ssu[j] * (ssu_size - 1)) / (ssu_size * denom)
        }
        psu_per <- vapply(seq_len(nr), psu_size_required_fn2, numeric(1L))
        psu_size_opt <- max(psu_per)
        if (!is.finite(psu_size_opt) || psu_size_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }
        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
      }
    } else {
      if (solve_for == "n1") {
        psu_size_opt <- psu_size
        ssu_size_opt <- ssu_size
        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        n1_opt <- max(n1_vals)
        total_cost <- fixed_cost +
          n1_opt * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
        binding_idx <- which.max(n1_vals)
      } else if (solve_for == "n2") {
        n1_opt <- n_psu
        ssu_size_opt <- ssu_size
        cv_floor <- .multistage_cv_floor(
          rel_var, k_psu, delta_psu,
          rr = rr, n_psu = n_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          rel_var, k_psu, delta_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )
        psu_per <- vapply(
          seq_len(nr),
          function(j) {
            denom <- cv_t[j]^2 * n_psu * rr[j] / (rel_var[j] * k_ssu[j]) -
              k_psu[j] * delta_psu[j] / k_ssu[j]
            if (denom <= 0) Inf
            else (1 + delta_ssu[j] * (ssu_size - 1)) / (ssu_size * denom)
          },
          numeric(1L)
        )
        psu_size_opt <- max(psu_per)
        if (!is.finite(psu_size_opt) || psu_size_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }
        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
      } else {
        n1_opt <- n_psu
        psu_size_opt <- psu_size
        cv_floor <- .multistage_cv_floor(
          rel_var, k_psu, delta_psu,
          k_ssu = k_ssu, delta_ssu = delta_ssu,
          rr = rr, n_psu = n_psu, psu_size = psu_size
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          rel_var, k_psu, delta_psu,
          k_ssu = k_ssu, delta_ssu = delta_ssu,
          psu_size = psu_size,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )
        ssu_per <- vapply(
          seq_len(nr),
          function(j) {
            a <- k_psu[j] * delta_psu[j]
            b <- k_ssu[j]
            d2 <- delta_ssu[j]
            denom <- psu_size * (cv_t[j]^2 * n_psu * rr[j] / rel_var[j] - a) -
              b * d2
            if (denom <= 0) Inf
            else b * (1 - d2) / denom
          },
          numeric(1L)
        )
        ssu_size_opt <- max(ssu_per)
        if (!is.finite(ssu_size_opt) || ssu_size_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }
        n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
      }
    }

    cv_achieved <- cv_achieved_fn(n1_opt, psu_size_opt, ssu_size_opt)
  } else {
    var_budget <- budget - fixed_cost
    bres <- .eval_3stage_budget(
      cv_t,
      delta_psu,
      delta_ssu,
      rel_var,
      k_psu,
      k_ssu,
      rr,
      C1,
      C2,
      C3,
      var_budget,
      n_psu,
      psu_size,
      ssu_size
    )
    n1_opt <- bres$n1
    psu_size_opt <- bres$psu_size
    ssu_size_opt <- bres$ssu_size
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- fixed_cost + bres$cost
  }

  n_vec <- c(n_psu = n1_opt, psu_size = psu_size_opt, ssu_size = ssu_size_opt)
  total_n <- prod(n_vec)
  operational <- .op_multi_3stage(
    n1_required, cv_achieved_fn, cv_t, stage_cost, budget,
    n_psu, psu_size, ssu_size, fixed_cost,
    cont_m = psu_size_opt, cont_q = ssu_size_opt
  )

  n1_per <- n1_required(psu_size_opt, ssu_size_opt)
  n_per <- n1_per * psu_size_opt * ssu_size_opt

  detail <- data.frame(
    name = labels,
    .n = n_per,
    .cv_target = cv_t,
    .cv_achieved = cv_achieved,
    .binding = seq_len(nr) == binding_idx
  )

  params <- list(stage_cost = stage_cost, domain_cols = domain_cols,
                  mode = mode, prop_method = prop_method)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(psu_size)) params$psu_size <- psu_size
  if (!is.null(ssu_size)) params$ssu_size <- ssu_size
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 3L,
    total_n = total_n,
    cv = cv_achieved[binding_idx],
    cost = total_cost,
    params = params,
    targets = targets,
    detail = detail,
    binding = labels[binding_idx],
    operational = operational
  )
}

#' Solve n_multi independently per domain, then aggregate
#' @keywords internal
#' @noRd
.n_multi_domains <- function(
  targets,
  stage_cost,
  budget,
  n_psu,
  psu_size = NULL,
  ssu_size = NULL,
  domain_cols,
  multistage,
  joint = FALSE,
  min_n = NULL,
  fixed_cost = 0,
  mode = "moe",
  prop_method = "wald"
) {
  key <- .domain_key(targets, domain_cols)
  domain_levels <- unique(key)
  domain_keys <- factor(key, levels = domain_levels)
  split_idx <- split(seq_len(nrow(targets)), domain_keys)

  if (multistage && joint && !is.null(budget) && length(domain_levels) > 1L) {
    return(.n_multi_domains_joint(
      targets,
      stage_cost,
      budget,
      n_psu,
      psu_size,
      ssu_size,
      domain_cols,
      domain_levels,
      split_idx,
      min_n,
      fixed_cost,
      mode = mode,
      prop_method = prop_method
    ))
  }

  results <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- targets[rows, , drop = FALSE]
    # Drop domain columns for inner solve
    sub_inner <- sub[, !names(sub) %in% domain_cols, drop = FALSE]
    if (multistage) {
      .n_multi_cluster(sub_inner, stage_cost, budget, n_psu, psu_size,
                       ssu_size, fixed_cost)
    } else {
      .n_multi_simple(sub_inner)
    }
  })
  names(results) <- domain_levels

  if (multistage) {
    .aggregate_cluster_domains(
      results,
      targets,
      domain_cols,
      domain_levels,
      split_idx,
      stage_cost,
      budget,
      n_psu,
      psu_size,
      ssu_size,
      min_n,
      fixed_cost,
      mode = mode,
      prop_method = prop_method
    )
  } else {
    .aggregate_simple_domains(
      results,
      targets,
      domain_cols,
      domain_levels,
      split_idx,
      min_n,
      mode = mode,
      prop_method = prop_method
    )
  }
}

#' Aggregate simple-mode domain results
#' @keywords internal
#' @noRd
.aggregate_simple_domains <- function(
  results,
  targets,
  domain_cols,
  domain_levels,
  split_idx,
  min_n = NULL,
  mode = "moe",
  prop_method = "wald"
) {
  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- targets[rows[1L], domain_cols, drop = FALSE]
    dom_vals$.n <- res$n
    dom_vals$.binding <- res$binding
    dom_vals
  })
  domains <- do.call(rbind, domain_rows)
  rownames(domains) <- NULL

  if (!is.null(min_n)) {
    floored <- domains$.n < min_n
    domains$.n[floored] <- min_n
    domains$.binding[floored] <- "(min_n)"
  }

  n_max <- max(domains$.n)
  overall_binding_idx <- which.max(domains$.n)
  binding_label <- domains$.binding[overall_binding_idx]

  res <- .new_svyplan_n(
    n = n_max,
    type = "multi",
    params = list(domain_cols = domain_cols, mode = mode,
                  prop_method = prop_method),
    targets = targets,
    detail = NULL,
    binding = binding_label,
    domains = domains
  )
  if (!is.null(min_n)) {
    res$params$min_n <- min_n
  }
  res
}

#' Aggregate multistage domain results
#' @keywords internal
#' @noRd
.aggregate_cluster_domains <- function(
  results,
  targets,
  domain_cols,
  domain_levels,
  split_idx,
  stage_cost,
  budget = NULL,
  n_psu = NULL,
  psu_size = NULL,
  ssu_size = NULL,
  min_n = NULL,
  fixed_cost = 0,
  mode = "cv",
  prop_method = "wald"
) {
  stages <- results[[1L]]$stages
  stage_names <- names(results[[1L]]$n)

  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- targets[rows[1L], domain_cols, drop = FALSE]
    for (s in seq_len(stages)) {
      dom_vals[[stage_names[s]]] <- res$n[s]
    }
    dom_vals$.total_n <- res$total_n
    dom_vals$.cv <- res$cv
    dom_vals$.cost <- res$cost
    dom_vals$.binding <- res$binding
    dom_vals
  })
  domains <- do.call(rbind, domain_rows)
  rownames(domains) <- NULL

  if (!is.null(min_n)) {
    below <- which(domains$.total_n < min_n)
    if (length(below) > 0L) {
      dom_labels <- vapply(
        below,
        function(i) {
          paste(domains[i, domain_cols, drop = TRUE], collapse = ":")
        },
        character(1L)
      )
      warning(
        sprintf(
          "domain(s) %s have total_n below min_n = %g",
          paste(sQuote(dom_labels), collapse = ", "),
          min_n
        ),
        call. = FALSE
      )
    }
  }

  total_n <- sum(ceiling(domains$.total_n))
  total_cost <- sum(domains$.cost)
  worst_cv_idx <- which.max(domains$.cv)

  n_vec <- vapply(
    seq_len(stages),
    function(s) {
      max(domains[[stage_names[s]]])
    },
    numeric(1L)
  )
  names(n_vec) <- stage_names

  params <- list(stage_cost = stage_cost, domain_cols = domain_cols,
                  mode = mode, prop_method = prop_method)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(psu_size)) params$psu_size <- psu_size
  if (!is.null(ssu_size)) params$ssu_size <- ssu_size
  if (!is.null(min_n)) {
    params$min_n <- min_n
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = stages,
    total_n = total_n,
    cv = domains$.cv[worst_cv_idx],
    cost = total_cost,
    params = params,
    targets = targets,
    detail = NULL,
    binding = domains$.binding[worst_cv_idx],
    domains = domains
  )
}

#' @keywords internal
#' @noRd
.n_multi_domains_joint <- function(
  targets,
  stage_cost,
  budget,
  n_psu,
  psu_size = NULL,
  ssu_size = NULL,
  domain_cols,
  domain_levels,
  split_idx,
  min_n = NULL,
  fixed_cost = 0,
  mode = "cv",
  prop_method = "wald"
) {
  stages <- length(stage_cost)
  nd <- length(domain_levels)
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- if (stages == 3L) stage_cost[3L] else NULL

  domain_params <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- targets[rows, , drop = FALSE]
    labels <- if ("name" %in% names(sub)) {
      sub$name
    } else {
      seq_along(rows)
    }
    list(
      cv_t = sub$cv,
      delta_psu = sub$delta_psu,
      delta_ssu = if ("delta_ssu" %in% names(sub)) {
        sub$delta_ssu
      } else {
        rep(0, length(rows))
      },
      rel_var = sub$rel_var,
      k_psu = sub$k_psu,
      k_ssu = sub$k_ssu,
      resp_rate = sub$resp_rate,
      labels = labels
    )
  })

  var_budget <- budget - fixed_cost

  eval_domain <- function(d, budget_d) {
    p <- domain_params[[d]]
    if (stages == 2L) {
      .eval_2stage_budget(
        p$cv_t,
        p$delta_psu,
        p$rel_var,
        p$k_psu,
        p$resp_rate,
        C1,
        C2,
        budget_d,
        n_psu,
        psu_size
      )
    } else {
      .eval_3stage_budget(
        p$cv_t,
        p$delta_psu,
        p$delta_ssu,
        p$rel_var,
        p$k_psu,
        p$k_ssu,
        p$resp_rate,
        C1,
        C2,
        C3,
        budget_d,
        n_psu,
        psu_size,
        ssu_size
      )
    }
  }

  total_n_for <- function(bres) {
    if (stages == 2L) {
      bres$n1 * bres$psu_size
    } else {
      bres$n1 * bres$psu_size * bres$ssu_size
    }
  }

  if (!is.null(min_n)) {
    full_total <- vapply(
      seq_len(nd),
      function(d) {
        tryCatch(total_n_for(eval_domain(d, var_budget)), error = function(e) 0)
      },
      numeric(1L)
    )
    for (d in seq_len(nd)) {
      if (full_total[d] < min_n) {
        lab <- paste(
          unlist(lapply(targets[split_idx[[domain_levels[d]]][1L],
                                domain_cols, drop = FALSE], as.character)),
          collapse = ":"
        )
        stop(
          sprintf(
            "min_n = %g not achievable for domain '%s' (max total_n = %.0f at full budget)",
            min_n,
            lab,
            full_total[d]
          ),
          call. = FALSE
        )
      }
    }
    min_fracs <- min_n / full_total
    if (sum(min_fracs) > 1) {
      stop(
        sprintf(
          "min_n = %g not achievable for all domains within budget",
          min_n
        ),
        call. = FALSE
      )
    }
  }

  lower_bounds <- rep(1e-4, nd)
  if (!is.null(min_n)) {
    lower_bounds <- pmax(lower_bounds, min_fracs)
  }

  outer_obj <- function(w) {
    w_last <- 1 - sum(w)
    if (w_last < lower_bounds[nd]) {
      return(1e12)
    }
    fracs <- c(w, w_last)
    tryCatch(
      {
        vals <- vapply(
          seq_len(nd),
          function(d) {
            bres <- eval_domain(d, fracs[d] * var_budget)
            if (!is.null(min_n) && total_n_for(bres) < min_n) {
              return(1e12)
            }
            bres$ratio
          },
          numeric(1L)
        )
        max(vals)
      },
      error = function(e) 1e12
    )
  }

  if (nd == 2L) {
    opt <- optimize(
      function(w1) outer_obj(w1),
      interval = c(lower_bounds[1], 1 - lower_bounds[2])
    )
    w_opt <- c(opt$minimum, 1 - opt$minimum)
  } else {
    slack <- 1 - sum(lower_bounds)
    if (slack < 0) {
      stop(
        "cumulative lower bounds exceed 1; joint allocation is infeasible",
        call. = FALSE
      )
    }
    init_w <- lower_bounds[-nd] + slack / nd
    opt <- optim(
      par = init_w,
      fn = outer_obj,
      method = "L-BFGS-B",
      lower = lower_bounds[-nd],
      upper = rep(1 - lower_bounds[nd], nd - 1L)
    )
    if (opt$convergence != 0L) {
      warning(
        "L-BFGS-B did not converge (code ",
        opt$convergence,
        "): ",
        "allocation may be approximate",
        call. = FALSE
      )
    }
    w_opt <- c(opt$par, 1 - sum(opt$par))
  }

  budgets <- w_opt * var_budget
  results <- lapply(seq_len(nd), function(d) {
    bres <- eval_domain(d, budgets[d])
    p <- domain_params[[d]]
    if (stages == 2L) {
      n_vec <- c(n_psu = bres$n1, psu_size = bres$psu_size)
    } else {
      n_vec <- c(
        n_psu = bres$n1,
        psu_size = bres$psu_size,
        ssu_size = bres$ssu_size
      )
    }
    .new_svyplan_cluster(
      n = n_vec,
      stages = stages,
      total_n = prod(n_vec),
      cv = bres$cv_achieved[bres$binding_idx],
      cost = bres$cost,
      params = list(stage_cost = stage_cost),
      targets = NULL,
      detail = NULL,
      binding = p$labels[bres$binding_idx]
    )
  })
  names(results) <- domain_levels

  res <- .aggregate_cluster_domains(
    results,
    targets,
    domain_cols,
    domain_levels,
    split_idx,
    stage_cost,
    budget,
    n_psu,
    psu_size,
    ssu_size,
    min_n,
    fixed_cost,
    mode = mode,
    prop_method = prop_method
  )
  res$params$joint <- TRUE
  if (fixed_cost > 0) {
    res$cost <- budget
  }
  res
}

#' @keywords internal
#' @noRd
.eval_2stage_budget <- function(
  cv_t,
  delta,
  rel_var,
  k,
  resp_rate,
  C1,
  C2,
  budget,
  n_psu,
  psu_size = NULL
) {
  nr <- length(cv_t)

  if (!is.null(psu_size)) {
    psu_size_opt <- psu_size
    n1_opt <- budget / (C1 + C2 * psu_size)
    if (n1_opt <= 0) {
      stop(
        "budget is too small for the given fixed cluster size",
        call. = FALSE
      )
    }
  } else if (is.null(n_psu)) {
    obj_fn <- function(psu_size) {
      n1 <- budget / (C1 + C2 * psu_size)
      if (n1 <= 0) {
        return(1e12)
      }
      cv_ratios <- vapply(
        seq_len(nr),
        function(j) {
          sqrt(
            rel_var[j] *
              k[j] /
              (n1 * resp_rate[j] * psu_size) *
              (1 + delta[j] * (psu_size - 1))
          ) /
            cv_t[j]
        },
        numeric(1L)
      )
      max(cv_ratios)
    }

    upper <- max(10, budget / (C1 + C2))
    opt <- optimize(obj_fn, interval = c(1, upper))
    psu_size_opt <- opt$minimum
    n1_opt <- budget / (C1 + C2 * psu_size_opt)
  } else {
    n1_opt <- n_psu
    psu_size_opt <- (budget - C1 * n_psu) / (C2 * n_psu)
    if (psu_size_opt <= 0) {
      stop(
        "budget is too small for the given fixed stage-1 size",
        call. = FALSE
      )
    }
  }

  cv_achieved <- vapply(
    seq_len(nr),
    function(j) {
      sqrt(
        rel_var[j] *
          k[j] /
          (n1_opt * resp_rate[j] * psu_size_opt) *
          (1 + delta[j] * (psu_size_opt - 1))
      )
    },
    numeric(1L)
  )
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    psu_size = psu_size_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = budget
  )
}

#' @keywords internal
#' @noRd
.eval_3stage_budget <- function(
  cv_t,
  delta_psu,
  delta_ssu,
  rel_var,
  k_psu,
  k_ssu,
  resp_rate,
  C1,
  C2,
  C3,
  budget,
  n_psu,
  psu_size = NULL,
  ssu_size = NULL
) {
  nr <- length(cv_t)

  cv_fn <- function(n1, psu_size, ssu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] /
            (n1 * resp_rate[j] * psu_size * ssu_size) *
            (k_psu[j] *
              delta_psu[j] *
              psu_size *
              ssu_size +
              k_ssu[j] * (1 + delta_ssu[j] * (ssu_size - 1)))
        )
      },
      numeric(1L)
    )
  }

  n_free <- 3L - sum(!is.null(n_psu), !is.null(psu_size), !is.null(ssu_size))

  if (n_free == 3L) {
    obj_fn_2d <- function(par) {
      ps <- par[1L]
      ss <- par[2L]
      n1 <- budget / (C1 + C2 * ps + C3 * ps * ss)
      if (n1 <= 0) return(1e12)
      max(cv_fn(n1, ps, ss) / cv_t)
    }
    init_ps <- max(2, sqrt(C1 / C2))
    init_ss <- max(2, sqrt(C2 / C3))
    upper_ps <- max(1000, 10 * init_ps)
    upper_ss <- max(1000, 10 * init_ss)
    opt <- optim(
      par = c(init_ps, init_ss),
      fn = obj_fn_2d,
      method = "L-BFGS-B",
      lower = c(1, 1),
      upper = c(upper_ps, upper_ss)
    )
    if (opt$convergence != 0L) {
      warning(
        "L-BFGS-B did not converge (code ",
        opt$convergence,
        "): allocation may be approximate",
        call. = FALSE
      )
    }
    if (opt$par[1L] >= upper_ps * (1 - 1e-6) ||
      opt$par[2L] >= upper_ss * (1 - 1e-6)) {
      warning(
        sprintf(
          "optimal stage size reached the search upper bound (psu_size <= %.0f, ssu_size <= %.0f); result may be unreliable -- review stage costs and target CVs",
          upper_ps,
          upper_ss
        ),
        call. = FALSE
      )
    }
    psu_size_opt <- opt$par[1L]
    ssu_size_opt <- opt$par[2L]
    n1_opt <- budget / (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
    total_cost <- budget
  } else if (n_free == 2L) {
    if (!is.null(n_psu)) {
      n1_opt <- n_psu
      avail <- budget / n_psu - C1
      if (avail < C2 + C3) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      max_ss <- (avail - C2) / C3
      if (!is.finite(max_ss) || max_ss < 1) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      ps_from_ss <- function(ss) avail / (C2 + C3 * ss)
      obj_fn_ss <- function(ss) {
        ps <- ps_from_ss(ss)
        if (!is.finite(ps) || ps < 1) return(Inf)
        max(cv_fn(n_psu, ps, ss) / cv_t)
      }
      if (max_ss <= 1 + 1e-10) {
        ssu_size_opt <- 1
      } else {
        opt <- optimize(obj_fn_ss, interval = c(1, max_ss))
        ssu_size_opt <- opt$minimum
      }
      psu_size_opt <- ps_from_ss(ssu_size_opt)
      if (!is.finite(psu_size_opt) || psu_size_opt < 1) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
    } else if (!is.null(psu_size)) {
      obj_fn_ps <- function(ss) {
        n1 <- budget / (C1 + C2 * psu_size + C3 * psu_size * ss)
        if (n1 <= 0) return(1e12)
        max(cv_fn(n1, psu_size, ss) / cv_t)
      }
      upper <- max(10, budget / (C1 + C2 * psu_size + C3 * psu_size))
      opt <- optimize(obj_fn_ps, interval = c(1, upper))
      ssu_size_opt <- opt$minimum
      psu_size_opt <- psu_size
      n1_opt <- budget / (C1 + C2 * psu_size + C3 * psu_size * ssu_size_opt)
    } else {
      obj_fn_ss2 <- function(ps) {
        n1 <- budget / (C1 + C2 * ps + C3 * ps * ssu_size)
        if (n1 <= 0) return(1e12)
        max(cv_fn(n1, ps, ssu_size) / cv_t)
      }
      upper <- max(10, budget / (C1 + C2 + C3 * ssu_size))
      opt <- optimize(obj_fn_ss2, interval = c(1, upper))
      psu_size_opt <- opt$minimum
      ssu_size_opt <- ssu_size
      n1_opt <- budget / (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size)
    }
    total_cost <- C1 * n1_opt +
      C2 * n1_opt * psu_size_opt +
      C3 * n1_opt * psu_size_opt * ssu_size_opt
  } else {
    if (!is.null(n_psu) && !is.null(psu_size)) {
      ssu_size_opt <- (budget - C1 * n_psu - C2 * n_psu * psu_size) /
        (C3 * n_psu * psu_size)
      if (ssu_size_opt <= 0) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      n1_opt <- n_psu
      psu_size_opt <- psu_size
    } else if (!is.null(n_psu) && !is.null(ssu_size)) {
      psu_size_opt <- (budget / n_psu - C1) / (C2 + C3 * ssu_size)
      if (psu_size_opt <= 0) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      n1_opt <- n_psu
      ssu_size_opt <- ssu_size
    } else {
      n1_opt <- budget / (C1 + C2 * psu_size + C3 * psu_size * ssu_size)
      if (n1_opt <= 0) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      psu_size_opt <- psu_size
      ssu_size_opt <- ssu_size
    }
    total_cost <- C1 * n1_opt +
      C2 * n1_opt * psu_size_opt +
      C3 * n1_opt * psu_size_opt * ssu_size_opt
  }

  cv_achieved <- cv_fn(n1_opt, psu_size_opt, ssu_size_opt)
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    psu_size = psu_size_opt,
    ssu_size = ssu_size_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = total_cost
  )
}

#' @rdname n_multi
#' @export
n_multi.svyplan_prec <- function(targets, ...) {
  x <- targets
  dots <- list(...)
  if (x$type != "multi") {
    stop("n_multi requires a svyplan_prec of type 'multi'", call. = FALSE)
  }
  if (identical(x$params$design, "cluster") ||
      !is.null(x$params$stage_cost)) {
    stop(
      "cluster precision must be passed to n_multi_cluster()",
      call. = FALSE
    )
  }
  tgt <- x$params$targets
  if ("prop_method" %in% names(dots)) {
    tgt$prop_method <- NA_character_
  }
  tgt$n <- NULL
  tgt$psu_size <- NULL
  tgt$ssu_size <- NULL

  stored_mode <- x$params$mode
  if (!is.null(stored_mode) && stored_mode == "moe") {
    tgt$moe <- x$detail$.moe
    tgt$cv <- NULL
  } else if (!is.null(stored_mode) && stored_mode %in% c("cv", "budget")) {
    tgt$cv <- x$detail$.cv
    tgt$moe <- NULL
  } else if (!is.null(x$detail)) {
    if (".moe" %in% names(x$detail) && !all(is.na(x$detail$.moe))) {
      tgt$moe <- x$detail$.moe
      tgt$cv <- NULL
    } else if (".cv" %in% names(x$detail) && !all(is.na(x$detail$.cv))) {
      tgt$cv <- x$detail$.cv
      tgt$moe <- NULL
    }
  }

  args <- list(
    targets = tgt,
    domains = x$params$domain_cols,
    min_n = x$params$min_n,
    prop_method = x$params$prop_method %||% "wald"
  )
  do.call(n_multi.default, .roundtrip_args(args, dots, n_multi.default))
}

#' @rdname n_multi_cluster
#' @export
n_multi_cluster.svyplan_prec <- function(targets, ...) {
  x <- targets
  dots <- list(...)
  if (x$type != "multi" || !identical(x$params$design, "cluster")) {
    stop(
      "n_multi_cluster requires cluster precision from prec_multi_cluster()",
      call. = FALSE
    )
  }

  tgt <- x$params$targets
  tgt$n <- NULL
  tgt$psu_size <- NULL
  tgt$ssu_size <- NULL

  stored_mode <- x$params$mode
  if (identical(stored_mode, "moe")) {
    tgt$moe <- x$detail$.moe
    tgt$cv <- NULL
  } else if (!is.null(x$detail) && ".cv" %in% names(x$detail)) {
    tgt$cv <- x$detail$.cv
    tgt$moe <- NULL
  }

  args <- list(
    targets = tgt,
    stage_cost = x$params$stage_cost,
    domains = x$params$domain_cols,
    budget = x$params$budget,
    n_psu = x$params$n_psu,
    psu_size = x$params$psu_size,
    ssu_size = x$params$ssu_size,
    joint = x$params$joint %||% FALSE,
    min_n = x$params$min_n,
    fixed_cost = x$params$fixed_cost %||% 0
  )
  do.call(
    n_multi_cluster.default,
    .roundtrip_args(args, dots, n_multi_cluster.default)
  )
}
