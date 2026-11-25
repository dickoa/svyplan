#' Grid exploration for svyplan objects
#'
#' Evaluate a svyplan result at new parameter combinations.
#' Returns a data frame with the varied parameters and resulting
#' quantities, suitable for sensitivity analysis or plotting.
#'
#' @param object A svyplan object (`svyplan_n`, `svyplan_multi`,
#'   `svyplan_alloc`, `svyplan_alloc_joint`, `svyplan_cluster`,
#'   `svyplan_power`, or `svyplan_prec`). For `svyplan_prec`, types
#'   `"proportion"`, `"mean"`, `"ratio"`, `"change"`, `"pooled"`, `"alloc"`,
#'   and `"alloc_joint"` (including PSU register allocations) are supported.
#'   Types `"cluster"` and `"multi"` are not.
#' @param newdata A data frame of parameter combinations to evaluate.
#'   Column names must be valid parameters for the object type (see
#'   Details). Parameters not in `newdata` stay at their original
#'   values from the object.
#' @param ... Additional arguments are not supported and produce an error.
#' @param control For a [n_cluster()] result only, the whole-unit search
#'   settings of each refit, as in [n_cluster()]. Defaults to the fit's own.
#'
#' @return A data frame with `newdata` columns followed by result
#'   columns. The result columns depend on the object type:
#'
#'   - `svyplan_n`: `n`, `se`, `moe`, `cv`, `rmoe`, except for a
#'     `type = "change"` result, which reports one size per occasion as
#'     `n1` and `n2` in place of `n`
#'   - `svyplan_alloc` from [n_alloc()]: `n` (the total), `se`, `moe`,
#'     `cv`, `rmoe` and `cost`
#'   - `svyplan_alloc_joint` from [n_alloc_joint()]: the design and cost
#'     columns described under Details
#'   - `svyplan_multi`: `n` and `binding`, the name of the binding indicator
#'   - `svyplan_cluster`: stage counts (`stage_1`, `stage_2`, ..., or custom stage names), `total_n`, `cv`, `cost`
#'   - `svyplan_power`: `n`, `power`, `effect`
#'   - `svyplan_prec`: `se`, `moe`, `cv`, `rmoe`, except for an allocation.
#'     A classic allocation reports `n` (the total), `se`, `moe`, `cv`,
#'     `rmoe` and `cost`, with `.max_domain_cv` when it has domains. A joint
#'     allocation reports `n`, `cost`, `objective_value` and `cv` when it has
#'     an objective, `.max_ratio` (the largest achieved over target ratio),
#'     `.n_fail` and `.feasible`. Given a list column `n`, the total is
#'     `total_n`.
#'
#'   When `newdata` sets a target, `cv`, `moe`, `rmoe` or `power`, the
#'   achieved value is reported beside it as `.cv_achieved`, `.moe_achieved`,
#'   `.rmoe_achieved` or `.power_achieved`. It equals the target unless a
#'   floor such as `min_cases` or `min_n_stratum` binds and carries the
#'   precision past it. A result column that repeats an input set in
#'   `newdata`, such as `n` or `take_2`, is not shown twice.
#'
#' @details
#' Valid parameters for `newdata` by object type:
#'
#' - **`n_prop`**: `p`, `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`,
#'   `resp_rate`, `eligibility_rate`, `eligible_per_unit`, `df`, `min_cases`
#' - **`n_mean`**: `var`, `mu`, `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`,
#'   `resp_rate`, `eligibility_rate`, `eligible_per_unit`, `df`
#' - **`n_ratio`**: `r`, `cv_num`, `cv_den`, `component_cor`, `moe`, `rmoe`,
#'   `cv`, `alpha`, `N`, `deff`, `resp_rate`, `eligibility_rate`, `df`
#' - **`n_change`**: `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `df`, `ratio`, `overlap`, `overlap_cor`, and `change` on
#'   the mean scale only, since two proportions or two ratios determine it.
#'   A ratio change also takes `cv_num`, `cv_den` and `component_cor`, one
#'   value per row applied to both occasions, but not `r`, a pair that a
#'   grid row cannot hold. Otherwise the occasion variances are held at the
#'   values the result was built from, so vary the design rather than the
#'   estimand. Sweeping `overlap` and `overlap_cor`
#'   is what prices a rotation against a fresh sample each round.
#' - **`n_pooled`**: `moe`, `rmoe`, `cv`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `df`, `mu` on the mean scale, `r`, `cv_num`,
#'   `cv_den` and `component_cor` on the ratio scale, and `occasions` only where
#'   the overlap and correlation are flat across lags, a shaped profile being
#'   a statement about a horizon that changing the horizon would contradict.
#'   The lag profiles themselves cannot be varied, a grid row having nowhere
#'   to hold a vector. `cor_decay` can, being one number that replaces the
#'   correlation profile, whatever the overlap's shape.
#' - **`n_multi`**: `target_mult` (multiplier on all `cv`, `moe`, `rmoe`
#'   targets), `deff`, `resp_rate`, `eligibility_rate`, `eligible_per_unit`,
#'   `alpha`, `df`.
#' - **`n_cluster`** (single-indicator): `cv`, `budget`, `unit_relvar`,
#'   `resp_rate`, stage response rates (`resp_rate_1`, ..., `resp_rate_{S-1}`),
#'   `eligibility_rate`, `eligible_per_unit`, `fixed_cost`, fixed stage
#'   sizes (`n_1`, `take_2`, ..., `take_S`), stage deltas (`icc` for 2-stage
#'   or `icc_1`, ..., `icc_{S-1}`), stage ratios (`var_ratio` for 2-stage or
#'   `var_ratio_1`, ..., `var_ratio_{S-1}`), and stage costs (`cost_1`, ...,
#'   `cost_S`). Fits with stage-by-stage finite population correction also
#'   accept `N_1`, ..., `N_{S-1}`, and `N`.
#' - **`n_cluster`** (multi-indicator `indicators` table): `target_mult`,
#'   `cost_mult`, `icc_mult`, `budget`, `resp_rate`, `eligibility_rate`,
#'   `eligible_per_unit`, `fixed_cost`. `cost_mult` scales the stage costs
#'   and leaves `fixed_cost` as it is, which its own column varies.
#' - **`power_prop`**: `p1`, `p2`, `n`, `power`, `alpha`, `N`, `deff`,
#'   `alternative`, `overlap`, `overlap_cor`, `resp_rate`, `eligibility_rate`,
#'   `eligible_per_unit` (excluding the solved-for parameter). Not supported
#'   for objects with vector `n`.
#' - **`power_mean`**: `effect`, `var`, `n`, `power`, `alpha`, `N`,
#'   `deff`, `alternative`, `overlap`, `overlap_cor`, `resp_rate`,
#'   `eligibility_rate`, `eligible_per_unit` (excluding the solved-for
#'   parameter). Not supported for objects with vector `n`.
#' - **`n_alloc`** (classic modes): `n`, `cv`, `moe`, `rmoe`, `budget`,
#'   `alpha`, `deff`, `resp_rate`, `eligibility_rate`, `eligible_per_unit`,
#'   `df`, `min_n_stratum`, `alloc_q`. A joint
#'   constrained allocation varies other controls, described below.
#' - **`power_did`**: `effect`, `n`, `power`, `alpha`, `N`, `deff`,
#'   `resp_rate`, `eligibility_rate`, `eligible_per_unit`, `alternative`,
#'   `overlap`, `overlap_cor`, `ratio` (excluding the solved-for parameter).
#' - **`prec_prop`**: `p`, `n`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `eligible_per_unit`, `df`. A result that solved for
#'   `p` varies the target it solved from, `cv` or `rmoe`, in place of `p`.
#' - **`prec_mean`**: `var`, `n`, `mu`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `eligible_per_unit`, `df`. A result that solved for `mu` varies `cv` in
#'   place of `mu`.
#' - **`prec_ratio`**: `r`, `n`, `cv_num`, `cv_den`, `component_cor`,
#'   `alpha`, `N`, `deff`, `resp_rate`, `eligibility_rate`, `df`
#' - **`prec_change`**: `n` for both occasions, or `n1` and `n2` for each
#'   (not both forms), `alpha`, `N`, `deff`, `resp_rate`, `eligibility_rate`,
#'   `df`, `overlap`, `overlap_cor`, `change` on the mean scale, and
#'   `cv_num`, `cv_den` and `component_cor` on the ratio scale, each applied
#'   to both occasions.
#' - **`prec_pooled`**: `n`, `alpha`, `N`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `df`, `mu` on the mean scale, and `r`, `cv_num`,
#'   `cv_den` and `component_cor` on the ratio scale. Where the
#'   overlap is the same at every lag, also `overlap`. Where the correlation
#'   is, also `overlap_cor`. Where both are, also `occasions`. `cor_decay`
#'   whatever the overlap's shape, since it replaces the correlation
#'   profile.
#' - **`prec_alloc`**: `n`, `alpha`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `eligible_per_unit`, `df`. A number in `n` rescales
#'   the stratum sizes to that total and keeps their shares. A list column `n`
#'   gives each row its own stratum sizes, one per row of the frame, to
#'   compare allocations. Takes held by the assessment after a whole-design
#'   search also apply to each scenario.
#' - **`prec_alloc_joint`**: fixed-design robustness sweeps over `alpha`,
#'   `deff`, `resp_rate`, `eligibility_rate`, `df`, `var_mult`, `cost_mult`,
#'   plus `icc_mult` (for multistage and PSU register fits) and `n` and
#'   `eligible_per_unit` (when `psu` is not used). The sample allocation is
#'   held fixed while assessment inputs vary.
#'
#' For `svyplan_n` objects, `moe`, `rmoe`, and `cv` are mutually exclusive
#' in `newdata`. If one appears, that mode is used. If none appears, the
#' original mode is preserved.
#'
#' Similarly, for `svyplan_cluster` objects, `cv` and `budget` are
#' mutually exclusive.
#'
#' A joint constrained allocation ([n_alloc_joint()]) re-solves the stored
#' problem at each scenario row. Valid columns include `target_mult` (when
#' targets are present), `cost_mult`, `deff`, `resp_rate`, `eligibility_rate`,
#' `alpha`, and `df`, `eligible_per_unit` without a PSU register, plus
#' `budget` in budget mode, `take_2`, ..., `take_S`
#' and `icc_mult` for multistage and register fits, `take_2_min` and `take_2_max`
#' for free-take fits, and `certainty_cutoff` and `n_psu_per_zone` for PSU
#' register fits.
#'
#' A fixed-budget objective fit varies `budget`, tracing the
#' cost-versus-objective frontier. A fixed-take multistage fit varies the
#' takes (`take_2`, ..., `take_S`), tracing what the fieldwork decision
#' costs. A scalar take applies to every stratum, so the per-stratum takes the
#' fit was built from are replaced rather than scaled. Sweep them
#' asymmetrically by editing `frame` and calling [n_alloc_joint()] per point. Both
#' controls may appear together on a budget-objective fit, one row per
#' combination.
#'
#' The result reports the continuous `n` and `cost`, the operational `n_int`
#' and `cost_int`, and `.feasible`. A multistage fit adds the PSU counts
#' `n_1` and `n_1_int`, and budget mode adds `objective_value`, its
#' equivalent `cv`, and whether the budget binds. On grids without a PSU
#' register, a rejected scenario or a failed solve returns missing results,
#' `.feasible = FALSE`, and a warning. Read the warning to distinguish a
#' failed solve from proven infeasibility. The grid checks before solving
#' that each requested take is a positive whole number.
#'
#' A certainty-aware fit ([n_alloc_joint()] with `psu`) varies `take_2` too,
#' and splits the PSU count `n_1_int` into `n_psu_certain` and `n_psu_draw`.
#' This is how the take is priced. It also sets the certainty threshold, so
#' a larger take leaves fewer certainty PSUs and costs more sample, and the
#' resulting cost curve has a minimum
#' worth finding rather than assuming. Such a fit also varies
#' `certainty_cutoff`, one value per row applied to every stratum, to price
#' a relaxed certainty rule the same way, and `n_psu_per_zone`, with `NA` for
#' no zones, to compare zone settings. A grid that varies `n_psu_per_zone`,
#' or a fit with zones, adds `n_zone`. A register fit sized to an
#' `objective` within a `budget` varies `budget` as well.
#'
#' A register grid adds `.status`. It reads `"optimal"`, `"gap"` or
#' `"limit"` for a row with a plan, as `$operational$status` does, `"infeasible"` where no
#' design exists, and `"unknown"` where the search stopped before finding
#' one or the solve failed, a row whose `.feasible` is `NA`. In budget mode
#' `objective_value` and `cv` are those of the field design. A row the node
#' limit stopped with a worse plan than a lower budget's, or with none,
#' takes that plan, which fits its budget, so the objective never worsens
#' as the budget grows. Its status reads `"limit"`.
#'
#' A cluster grid from [n_cluster()] adds `.status` and `.nodes`, the
#' status and node count of each row's whole-unit search. A row whose
#' search stopped before finding a design reads `"unknown"`, with its
#' sizes missing.
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
#' # A floor of 400 expected cases holds n at 1334 for any cv target above
#' # 0.042, so the achieved cv stays at 0.042 there
#' rare <- n_prop(p = 0.3, cv = 0.05, min_cases = 400)
#' predict(rare, data.frame(cv = c(0.03, 0.05, 0.08)))
#'
#' # Power curve: how does power vary with sample size?
#' pw <- power_prop(p1 = 0.30, p2 = 0.35, n = 500, power = NULL)
#' predict(pw, data.frame(n = seq(100, 1000, 100)))
#'
#' # Cluster design: sensitivity to icc (homogeneity)
#' cl <- n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
#' predict(cl, data.frame(icc = c(0.01, 0.03, 0.05, 0.10, 0.15)))
#'
#' # Indicator names repeated across domains: each domain is sized on its
#' # own, and the identity row target_mult = 1 returns the fitted n
#' dtab <- data.frame(region = rep(c("N", "S"), each = 2),
#'                    name = rep(c("y1", "y2"), 2),
#'                    p = c(0.3, 0.5, 0.35, 0.45), cv = 0.1)
#' dm <- n_multi(dtab, domains = "region")
#' predict(dm, data.frame(target_mult = c(1, 1.5)))
#'
#' # Stage costs 20% higher, the fixed cost of 10000 unchanged
#' tab <- data.frame(name = c("a", "b"), p = c(0.3, 0.4), cv = 0.05,
#'                   icc_1 = c(0.05, 0.1))
#' fit <- n_cluster(tab, stage_cost = c(500, 50), budget = 1e5,
#'                  fixed_cost = 1e4)
#' predict(fit, data.frame(cost_mult = c(1, 1.2)))
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
#' # The precision of that allocation if the budget is cut or raised, the
#' # stratum shares held, then three candidate allocations side by side
#' fixed <- prec_alloc(alloc, n = alloc$detail$n_int)
#' predict(fixed, data.frame(n = c(450, 600, 750)))
#' predict(fixed, data.frame(n = I(list(c(200, 200, 200),
#'                                      alloc$detail$n_int,
#'                                      c(100, 200, 300)))))
#'
#' # The precision of a change between two rounds as the panel overlap
#' # and its correlation vary
#' ch <- prec_change(p = c(0.30, 0.36), n = 1000, overlap = 0.5,
#'                   overlap_cor = 0.6)
#' predict(ch, expand.grid(n = c(800, 1200), overlap_cor = c(0.4, 0.7)))
#'
#' # Two-stage joint allocation: what does the per-PSU take cost?
#' cluster_frame <- data.frame(
#'   stratum = c("Urban", "Rural"),
#'   N       = c(60000, 40000),
#'   N_1     = c(600, 500),
#'   take_2  = c(10, 10),
#'   cost_1  = c(400, 550),
#'   cost_2  = c(20, 25)
#' )
#' cluster_measures <- data.frame(
#'   stratum = rep(cluster_frame$stratum, 2),
#'   name    = rep(c("literacy", "income"), each = 2),
#'   p       = c(0.62, 0.48, NA, NA),
#'   mean    = c(NA, NA, 520, 380),
#'   sd      = c(NA, NA, 210, 160),
#'   icc_1   = rep(c(0.04, 0.07), 2)
#' )
#' cluster_targets <- data.frame(
#'   name = c("literacy", "income"),
#'   cv   = c(0.03, 0.04)
#' )
#' joint <- n_alloc_joint(cluster_frame, measures = cluster_measures,
#'                        targets = cluster_targets)
#' predict(joint, data.frame(take_2 = c(6, 10, 14, 20)))
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
  if (!is.null(object$indicators)) {
    stop(
      "predict() is not supported for multi-indicator results; ",
      "use the underlying single-indicator functions instead",
      call. = FALSE
    )
  }

  if (object$type == "proportion") {
    allowed <- c("p", "moe", "cv", "rmoe", "alpha", "N", "deff", "resp_rate",
                 "eligibility_rate", "eligible_per_unit", "df", "min_cases")
    base <- object$params
    method <- object$method %||% "wald"

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- n_prop.default(
        p = p$p, moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p$eligible_per_unit %||% 1,
        method = method, df = p$df, min_cases = p$min_cases
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
    })

  } else if (object$type == "mean") {
    allowed <- c("var", "mu", "moe", "cv", "rmoe", "alpha", "N", "deff",
                 "resp_rate",
                 "eligibility_rate", "eligible_per_unit", "df")
    base <- object$params

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- n_mean.default(
        var = p$var, mu = p$mu, moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p[["eligible_per_unit"]] %||% 1, df = p$df
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
    })

  } else if (object$type == "ratio") {
    allowed <- c("r", "cv_num", "cv_den", "component_cor", "moe", "cv",
                 "rmoe", "alpha", "N", "deff", "resp_rate",
                 "eligibility_rate", "df")
    base <- object$params

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- n_ratio.default(
        r = p$r, cv_num = p$cv_num, cv_den = p$cv_den,
        component_cor = p$component_cor,
        moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N,
        deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1, df = p$df
      )
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
    })

  } else if (object$type == "change") {
    base <- object$params
    # 'change' is p[2] - p[1] on the proportion scale, so varying it there
    # would ask for a change the two proportions contradict.
    # A ratio's moments may be varied, its two ratios not, since a grid row
    # holds one number and the change needs a pair.
    allowed <- c("moe", "cv", "rmoe", "alpha", "N", "deff", "resp_rate",
                 "eligibility_rate", "df", "ratio", "overlap", "overlap_cor",
                 if (.repeated_free_level(base)) "change",
                 if (!is.null(base[["r"]])) c("cv_num", "cv_den",
                                              "component_cor"))

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")

    .predict_grid(newdata, base, function(p) {
      res <- do.call(n_change.default, c(.repeated_scale_args(p, "change"), list(
        moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N, deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        ratio = p[["ratio"]], overlap = p$overlap,
        overlap_cor = p$overlap_cor,
        overlap_draw = p$overlap_draw, df = p$df
      )))
      data.frame(n1 = res$n[1L], n2 = res$n[length(res$n)],
                 se = res$se, moe = res$moe, cv = res$cv, rmoe = res$rmoe)
    })

  } else if (object$type == "pooled") {
    base <- object$params
    # An issued profile is a respondent one only at full response. The
    # stored rotation carries that basis, so each row is checked on its own.
    flat <- length(unique(base$overlap)) == 1L &&
      length(unique(base$overlap_cor)) == 1L
    allowed <- c("moe", "cv", "rmoe", "alpha", "N", "deff", "resp_rate",
                 "eligibility_rate", "df", "cor_decay",
                 if (flat) "occasions", if (.repeated_free_level(base)) "mu",
                 if (!is.null(base[["r"]])) c("r", "cv_num", "cv_den",
                                              "component_cor"))

    .validate_newdata(newdata, allowed)
    base <- .resolve_exclusive(newdata, base, "moe", "cv", "rmoe")
    if (flat) {
      base$overlap <- base$overlap[1L]
      base$overlap_cor <- base$overlap_cor[1L]
    }
    ## The grid varies the expanded lags, not the inputs a round trip reads.
    base$overlap_input <- NULL
    base$overlap_cor_input <- NULL
    if (!"cor_decay" %in% names(newdata)) base$cor_decay <- NULL
    if ("cor_decay" %in% names(newdata)) {
      base$overlap_cor <- NULL
    }

    .predict_grid(newdata, base, function(p) {
      res <- do.call(n_pooled.default, c(.repeated_scale_args(p, "mu"), list(
        occasions = p$occasions,
        moe = p$moe, cv = p$cv, rmoe = p$rmoe,
        alpha = p$alpha, N = p$N, deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        overlap_cor = p$overlap_cor,
        cor_decay = p$cor_decay, overlap_draw = p$overlap_draw, df = p$df
      ), .stored_overlap_args(p)))
      data.frame(n = res$n, se = res$se, moe = res$moe, cv = res$cv,
                 rmoe = res$rmoe)
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
predict.svyplan_alloc <- function(object, newdata, ...) {
  .check_unused_dots(...)
  allowed <- c("n", "cv", "moe", "rmoe", "budget", "alpha", "deff", "resp_rate",
               "eligibility_rate", "eligible_per_unit", "df", "min_n_stratum",
               "alloc_q")
  base <- object$params
  base_small <- list(
    n = base[["n"]], cv = base[["cv"]], moe = base[["moe"]],
    rmoe = base[["rmoe"]], budget = base[["budget"]],
    alpha = base$alpha, deff = base$deff, resp_rate = base$resp_rate,
    eligibility_rate = base[["eligibility_rate"]] %||% 1,
    eligible_per_unit = base[["eligible_per_unit"]] %||% 1,
    df = base$df,
    min_n_stratum = base$min_n_stratum, alloc_q = base$alloc_q %||% 0.5
  )

  .validate_newdata(newdata, allowed)

  mode_cols <- c("n", "cv", "moe", "rmoe", "budget")
  present <- mode_cols[mode_cols %in% names(newdata)]
  if (length(present) > 1L) {
    stop(
      "newdata can vary at most one of 'n', 'cv', 'moe', 'rmoe', or 'budget'",
      call. = FALSE
    )
  }
  if (length(present) == 1L) {
    for (col in setdiff(mode_cols, present)) base_small[[col]] <- NULL
  }

  .predict_grid(newdata, base_small, function(p) {
    count_mode <- (!is.null(p$n)) + (!is.null(p$cv)) + (!is.null(p$moe)) +
      (!is.null(p$rmoe)) + (!is.null(p$budget))
    if (count_mode != 1L) {
      stop("each predict row must define exactly one of n/cv/moe/rmoe/budget",
           call. = FALSE)
    }
    frame_i <- base$frame
    for (col in intersect(
      c("resp_rate", "deff", "eligibility_rate", "eligible_per_unit"),
      names(newdata)
    )) {
      frame_i[[col]] <- NULL
    }
    alloc_args <- list(
      frame = frame_i,
      domains = base$domain_cols,
      n = p$n, cv = p$cv, moe = p$moe, rmoe = p$rmoe, budget = p$budget,
      alloc = base$alloc %||% "neyman",
      alpha = p$alpha,
      deff = p$deff,
      resp_rate = p$resp_rate,
      eligibility_rate = p[["eligibility_rate"]] %||% 1,
      eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
      df = p$df,
      min_n_stratum = p$min_n_stratum,
      min_n_domain = base[["min_n_domain"]],
      alloc_q = p$alloc_q
    )
    row_mode <- if (!is.null(p$n)) {
      "n"
    } else if (!is.null(p$cv) || !is.null(p$moe) || !is.null(p$rmoe)) {
      "cv"
    } else {
      "budget"
    }
    if (.classic_searched(alloc_args$alloc, row_mode,
                          length(base$domain_cols) > 0L)) {
      alloc_args$control <- base[["control"]]
    }
    if (!.alloc_is_cluster(base$frame)) {
      alloc_args$unit_cost <- base$cost_h
    }
    alloc_args$fpc <- base$fpc %||% "unit"
    res <- do.call(n_alloc.default, alloc_args)
    data.frame(
      n = res$n, se = res$se, moe = res$moe, cv = res$cv, rmoe = res$rmoe,
      cost = res$params$achieved$cost
    )
  })
}

#' @rdname predict.svyplan
#' @export
predict.svyplan_multi <- function(object, newdata, ...) {
  .check_unused_dots(...)
  allowed <- c("target_mult", "deff", "resp_rate", "eligibility_rate",
               "eligible_per_unit", "alpha", "df")
  .validate_newdata(newdata, allowed)
  p <- object$params
  ind_base <- object$indicators
  base <- list(
    target_mult = 1,
    deff = 1,
    resp_rate = p$resp_rate %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    alpha = 0.05,
    df = NA_real_
  )
  varied <- names(newdata)
  dom_cols <- .predict_domain_cols(object)
  .predict_grid(newdata, base, function(r) {
    tbl <- .predict_scale_targets(ind_base, r$target_mult)
    if ("deff" %in% varied) tbl$deff <- r$deff
    if ("alpha" %in% varied) tbl$alpha <- r$alpha
    if ("df" %in% varied) tbl$df <- r$df
    if ("resp_rate" %in% varied && "resp_rate" %in% names(tbl)) tbl$resp_rate <- NULL
    if ("eligibility_rate" %in% varied && "eligibility_rate" %in% names(tbl)) {
      tbl$eligibility_rate <- NULL
    }
    if ("eligible_per_unit" %in% varied && "eligible_per_unit" %in% names(tbl)) {
      tbl$eligible_per_unit <- NULL
    }
    res <- n_multi.default(
      indicators = tbl,
      domains = dom_cols,
      domain_sampling = p[["domain_sampling"]] %||% "separate",
      min_n_domain = p[["min_n_domain"]],
      prop_method = p[["prop_method"]] %||% "wald",
      resp_rate = r$resp_rate,
      eligibility_rate = r$eligibility_rate,
      eligible_per_unit = r$eligible_per_unit
    )
    data.frame(
      n = res$n,
      binding = paste(res$binding, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
}

#' @rdname predict.svyplan
#' @export
predict.svyplan_alloc_joint <- function(object, newdata, ...) {
  .check_unused_dots(...)
  .predict_bethel(object, newdata)
}

#' @rdname predict.svyplan
#' @export
predict.svyplan_cluster <- function(object, newdata, ..., control = NULL) {
  .check_unused_dots(...)
  control <- control %||% object$params[["control"]]
  .stage_max_nodes(control)
  if (!is.null(object$indicators)) {
    return(.predict_cluster_multi(object, newdata, control))
  }

  stages <- length(object$params$stage_cost)
  if (stages > 2L && "icc" %in% names(newdata)) {
    stop(
      sprintf(
        "for %d-stage cluster predict, vary %s instead of 'icc'",
        stages,
        paste(sQuote(paste0("icc_", seq_len(stages - 1L))), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (stages > 2L && "var_ratio" %in% names(newdata)) {
    stop(
      sprintf(
        "for %d-stage cluster predict, vary %s instead of 'var_ratio'",
        stages,
        paste(sQuote(paste0("var_ratio_", seq_len(stages - 1L))), collapse = ", ")
      ),
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
  rr_stage_cols <- paste0("resp_rate_", seq_len(stages - 1L))
  sn_cols <- c("n_1", paste0("take_", 2L:stages))
  N_stage_cols <- c(paste0("N_", seq_len(stages - 1L)), "N")

  allowed <- unique(c(
    "cv", "budget", "unit_relvar",
    rr_stage_cols, "resp_rate", "eligibility_rate", "eligible_per_unit",
    "fixed_cost", sn_cols,
    if (identical(object$params[["fpc"]], "stage")) N_stage_cols,
    icc_meta$allowed, k_meta$allowed, cost_meta$allowed
  ))
  .validate_newdata(newdata, allowed)

  p <- object$params
  has_cv_col <- "cv" %in% names(newdata)
  has_budget_col <- "budget" %in% names(newdata)
  if (has_cv_col && has_budget_col) {
    stop("newdata cannot contain both 'cv' and 'budget'", call. = FALSE)
  }

  base_rr <- .expand_stage_resp_rate(p[["resp_rate"]] %||% 1, stages)
  base_sn <- p[["stage_n"]] %||% rep(NA_real_, stages)
  base_N <- p[["N"]]

  base <- list(
    stage_cost = p$stage_cost, icc = p$icc, unit_relvar = p$unit_relvar,
    var_ratio = p[["var_ratio_input"]] %||% p$var_ratio,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
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
      base$icc, params, icc_meta$map, icc_meta$canonical
    )
    row_k <- .apply_cluster_stage_cols(
      base$var_ratio, params, k_meta$map, k_meta$canonical
    )
    row_rr <- base_rr
    for (s in seq_len(stages - 1L)) {
      if (rr_stage_cols[s] %in% names(newdata)) {
        row_rr[s] <- params[[rr_stage_cols[s]]]
      }
    }
    if ("resp_rate" %in% names(newdata)) {
      row_rr[stages] <- params[["resp_rate"]]
    }
    row_sn <- base_sn
    for (s in seq_len(stages)) {
      if (sn_cols[s] %in% names(newdata)) {
        row_sn[s] <- params[[sn_cols[s]]]
      }
    }
    if (all(is.na(row_sn))) row_sn <- NULL
    row_N <- base_N
    if (!is.null(row_N)) {
      for (s in seq_len(stages)) {
        if (N_stage_cols[s] %in% names(newdata)) {
          row_N[s] <- params[[N_stage_cols[s]]]
        }
      }
    }
    res <- tryCatch(n_cluster.default(
      stage_cost = row_cost, icc = row_icc,
      unit_relvar = params$unit_relvar, var_ratio = row_k,
      cv = params$cv, budget = params$budget,
      stage_n = row_sn,
      resp_rate = row_rr,
      eligibility_rate = params[["eligibility_rate"]] %||% 1,
      eligible_per_unit = params[["eligible_per_unit"]] %||% 1,
      fixed_cost = params$fixed_cost,
      fpc = p[["fpc"]] %||% "auto",
      N = row_N,
      control = control
    ), svyplan_plan_unknown = function(e) NULL)
    if (is.null(res)) {
      return(.predict_cluster_unknown(names(object$n), c("total_n", "cv", "cost")))
    }
    out <- as.list(res$n)
    out$total_n <- res$total_n
    out$cv <- res$cv
    out$cost <- res$cost
    out$.status <- res$operational$search$status %||% NA_character_
    out$.nodes <- res$operational$search$nodes %||% NA_real_
    as.data.frame(out)
  })
}

#' A scenario row whose whole search found no design within its node limit
#' @keywords internal
#' @noRd
.predict_cluster_unknown <- function(stage_names, cols, text = character(0)) {
  out <- as.list(rep(NA_real_, length(stage_names) + length(cols)))
  names(out) <- c(stage_names, cols)
  for (nm in text) out[[nm]] <- NA_character_
  out$.status <- "unknown"
  out$.nodes <- NA_real_
  as.data.frame(out, stringsAsFactors = FALSE)
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
                    "resp_rate", "eligibility_rate", "eligible_per_unit", "alternative",
                    "overlap", "overlap_cor")
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
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
        alternative = p$alternative, overlap = p$overlap, overlap_cor = p$overlap_cor,
        overlap_draw = p$overlap_draw, method = method
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_prop, args)
      data.frame(n = res$n, power = res$power, effect = res$effect)
    })

  } else if (object$type == "mean") {
    all_params <- c("effect", "var", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "eligibility_rate", "eligible_per_unit", "alternative",
                    "overlap", "overlap_cor")
    excluded <- switch(solved, n = "n", power = "power", mde = "effect")
    allowed <- setdiff(all_params, excluded)

    .validate_newdata(newdata, allowed, allowed_character = "alternative")
    base <- object$params

    .predict_grid(newdata, base, function(p) {
      args <- list(
        effect = p$effect, var = p$var, n = p$n, power = p$power,
        alpha = p$alpha, N = p$N, deff = p$deff,
        resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
        alternative = p$alternative, overlap = p$overlap, overlap_cor = p$overlap_cor,
        overlap_draw = p$overlap_draw
      )
      args[excluded] <- list(NULL)
      res <- do.call(power_mean, args)
      data.frame(n = res$n, power = res$power, effect = res$effect)
    })

  } else if (object$type %in% c("did_prop", "did_mean")) {
    all_params <- c("effect", "n", "power", "alpha", "N", "deff",
                    "resp_rate", "eligibility_rate", "eligible_per_unit", "alternative",
                    "overlap", "overlap_cor", "ratio")
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
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
        alternative = p$alternative, ratio = p$ratio,
        overlap = p$overlap, overlap_cor = p$overlap_cor,
        overlap_draw = p$overlap_draw
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
    allowed <- c("p", "n", "alpha", "N", "deff", "resp_rate",
                 "eligibility_rate", "eligible_per_unit", "df")
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
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p$eligible_per_unit %||% 1,
        method = method, df = p$df
      )
      out <- data.frame(se = res$se, moe = res$moe, cv = res$cv,
                        rmoe = res$rmoe)
      if (identical(solved, "p")) out <- cbind(data.frame(p = res$params$p), out)
      out
    })

  } else if (object$type == "mean") {
    allowed <- c("var", "n", "mu", "alpha", "N", "deff", "resp_rate",
                 "eligibility_rate", "eligible_per_unit", "df")
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
        N = p$N, deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        eligible_per_unit = p[["eligible_per_unit"]] %||% 1, df = p$df
      )
      out <- data.frame(se = res$se, moe = res$moe, cv = res$cv,
                        rmoe = res$rmoe)
      if (identical(solved, "mu")) {
        out <- cbind(data.frame(mu = res$params$mu), out)
      }
      out
    })

  } else if (object$type == "ratio") {
    allowed <- c("r", "n", "cv_num", "cv_den", "component_cor", "alpha", "N",
                 "deff", "resp_rate", "eligibility_rate", "df")
    base <- object$params

    .validate_newdata(newdata, allowed)

    .predict_grid(newdata, base, function(p) {
      res <- prec_ratio.default(
        r = p$r, n = p$n, cv_num = p$cv_num, cv_den = p$cv_den,
        component_cor = p$component_cor, alpha = p$alpha,
        N = p$N, deff = p$deff, resp_rate = p$resp_rate,
        eligibility_rate = p[["eligibility_rate"]] %||% 1, df = p$df
      )
      data.frame(se = res$se, moe = res$moe, cv = res$cv, rmoe = res$rmoe)
    })

  } else if (object$type == "pooled") {
    .predict_prec_pooled(object, newdata)
  } else if (object$type == "change") {
    .predict_prec_change(object, newdata)
  } else if (object$type %in% c("alloc", "alloc_joint")) {
    .predict_prec_alloc(object, newdata)
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

#' A varied response rate replaces the final stage's rate only
#'
#' A stored per-stage vector keeps its earlier stages, as the
#' single-indicator method does.
#' @keywords internal
#' @noRd
.predict_final_rate <- function(stored, value, varied) {
  if (!varied || length(stored) < 2L || length(value) != 1L) {
    return(value)
  }
  stored[length(stored)] <- value
  stored
}

#' Grid exploration for a pooled precision
#'
#' A flat lag profile is one number, so it can be varied, along with the
#' number of occasions it spans. A shaped profile is a statement about a
#' horizon and stays as it is. `cor_decay` replaces the correlation profile
#' at any overlap shape.
#' @keywords internal
#' @noRd
.predict_prec_pooled <- function(object, newdata) {
  base <- object$params
  flat_ov <- length(unique(base$overlap)) == 1L
  flat_cor <- length(unique(base$overlap_cor)) == 1L
  allowed <- c("n", "alpha", "N", "deff", "resp_rate", "eligibility_rate",
               "df", if (.repeated_free_level(base)) "mu",
               if (!is.null(base[["r"]])) c("r", "cv_num", "cv_den",
                                            "component_cor"),
               "cor_decay",
               if (flat_ov && flat_cor) "occasions",
               if (flat_ov) "overlap",
               if (flat_cor) "overlap_cor")
  .validate_newdata(newdata, allowed)
  if (flat_ov) base$overlap <- base$overlap[1L]
  if (flat_cor) base$overlap_cor <- base$overlap_cor[1L]
  ## The grid varies the expanded lags, not the inputs a round trip reads.
  base$overlap_input <- NULL
  base$overlap_cor_input <- NULL
  if (!"cor_decay" %in% names(newdata)) base$cor_decay <- NULL
  if ("cor_decay" %in% names(newdata)) {
    if ("overlap_cor" %in% names(newdata)) {
      stop("newdata cannot contain both 'overlap_cor' and 'cor_decay'",
           call. = FALSE)
    }
    base$overlap_cor <- NULL
  }

  .predict_grid(newdata, base, function(p) {
    res <- do.call(prec_pooled.default, c(.repeated_scale_args(p, "mu"), list(
      n = p$n, occasions = p$occasions,
      alpha = p$alpha, N = p$N, deff = p$deff, resp_rate = p$resp_rate,
      eligibility_rate = p[["eligibility_rate"]] %||% 1,
      overlap_cor = p$overlap_cor,
      cor_decay = p$cor_decay, overlap_draw = p$overlap_draw, df = p$df
    ), .stored_overlap_args(p, "overlap" %in% names(newdata))))
    data.frame(se = res$se, moe = res$moe, cv = res$cv, rmoe = res$rmoe)
  })
}

#' Grid exploration for the precision of a change
#'
#' The two occasions can be sized alike through `n`, or apart through `n1`
#' and `n2`, which start from the object's own sizes.
#' @keywords internal
#' @noRd
.predict_prec_change <- function(object, newdata) {
  base <- object$params
  allowed <- c("n", "n1", "n2", "alpha", "N", "deff", "resp_rate",
               "eligibility_rate", "df", "overlap", "overlap_cor",
               if (.repeated_free_level(base)) "change",
               if (!is.null(base[["r"]])) c("cv_num", "cv_den",
                                            "component_cor"))
  .validate_newdata(newdata, allowed)
  pair <- any(c("n1", "n2") %in% names(newdata))
  if (pair && "n" %in% names(newdata)) {
    stop("newdata cannot contain 'n' together with 'n1' or 'n2'",
         call. = FALSE)
  }
  n_obj <- base$n
  base$n1 <- n_obj[1L]
  base$n2 <- n_obj[length(n_obj)]

  .predict_grid(newdata, base, function(p) {
    n <- if (pair) c(p$n1, p$n2) else p$n
    res <- do.call(prec_change.default, c(.repeated_scale_args(p, "change"), list(
      n = n,
      alpha = p$alpha, N = p$N, deff = p$deff, resp_rate = p$resp_rate,
      eligibility_rate = p[["eligibility_rate"]] %||% 1,
      overlap = p$overlap, overlap_cor = p$overlap_cor,
      overlap_draw = p$overlap_draw, df = p$df
    )))
    data.frame(se = res$se, moe = res$moe, cv = res$cv, rmoe = res$rmoe)
  })
}

#' Grid exploration for the precision of a stratum allocation
#'
#' A grid row cannot hold one size per stratum, so `n` is read two ways. A
#' number rescales the object's stratum sizes to that total, keeping their
#' shares, which is the allocation that was decided at another budget. A
#' list column gives each row its own stratum vector, to compare candidate
#' allocations side by side.
#' @keywords internal
#' @noRd
.predict_prec_alloc <- function(object, newdata) {
  p <- object$params
  joint <- identical(object$type, "alloc_joint") || identical(object$method, "bethel")
  is_psu <- !is.null(p[["psu"]])
  stages <- p$stages %||% 1L
  if (is_psu && "n" %in% names(newdata)) {
    stop(
      "cannot vary 'n' on a PSU register precision result; vary assessment inputs such as 'icc_mult', 'var_mult', 'cost_mult', 'deff', 'resp_rate', or 'alpha'",
      call. = FALSE
    )
  }
  allowed <- if (joint) {
    c(
      if (!is_psu) "n",
      "alpha", "deff", "resp_rate", "eligibility_rate",
      if (!is_psu) "eligible_per_unit",
      "df", "var_mult", "cost_mult",
      if (stages >= 2L || is_psu) "icc_mult"
    )
  } else {
    c("n", "alpha", "deff", "resp_rate", "eligibility_rate",
      "eligible_per_unit", "df")
  }
  .validate_newdata(newdata, allowed, allowed_list = "n")
  H <- length(p$n)
  varied <- names(newdata)
  if ("icc_mult" %in% varied) {
    for (m in newdata$icc_mult) {
      .predict_scale_icc(p$measures, m)
      .predict_scale_icc(p$frame, m)
    }
  }
  given_n <- "n" %in% varied
  list_n <- given_n && is.list(newdata[["n"]])
  if (list_n) {
    ok <- vapply(newdata[["n"]], function(v) {
      is.numeric(v) && length(v) == H
    }, logical(1))
    if (!all(ok)) {
      stop(
        sprintf(
          "each element of the list column 'n' must give %d stratum sizes, one per row of the frame",
          H
        ),
        call. = FALSE
      )
    }
  }
  base <- list(
    n = p$n, alpha = p$alpha, deff = p$deff, resp_rate = p$resp_rate,
    eligibility_rate = p[["eligibility_rate"]] %||% 1, df = p$df,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    var_mult = 1, cost_mult = 1, icc_mult = 1
  )
  if (!joint) assessment_frame <- .alloc_frame_at_takes(p$frame, p$searched_takes)

  .predict_grid(newdata, base, function(r) {
    n_h <- if (list_n) {
      as.numeric(r$n)
    } else if (given_n) {
      p$n * r$n / sum(p$n)
    } else {
      p$n
    }
    if (!joint) {
      frame_i <- assessment_frame
      for (col in intersect(
        c("resp_rate", "deff", "eligibility_rate", "eligible_per_unit"), varied
      )) {
        frame_i[[col]] <- NULL
      }
      args <- list(
        frame = frame_i, n = n_h, alpha = r$alpha, deff = r$deff,
        resp_rate = r$resp_rate, eligibility_rate = r$eligibility_rate,
        df = r$df, fpc = p$fpc %||% "unit",
        domains = p$domain_cols,
        eligible_per_unit = r$eligible_per_unit
      )
      if (!.alloc_is_cluster(p$frame)) args$unit_cost <- p$cost_h
      res <- do.call(prec_alloc.default, args)
      total <- res$params$achieved$n
      out <- data.frame(n = total, se = res$se, moe = res$moe, cv = res$cv,
                        rmoe = res$rmoe, cost = res$params$achieved$cost)
      if (!is.null(res$domains)) {
        out$.max_domain_cv <- max(res$domains$.cv, na.rm = TRUE)
      }
      if (list_n) names(out)[1L] <- "total_n"
      return(out)
    }

    frame_i <- p$frame
    measures_i <- p$measures
    unit_cost_i <- p$unit_cost
    deff_i <- r$deff
    if ("deff" %in% varied) {
      if ("deff" %in% names(frame_i)) frame_i$deff <- NULL
      if ("deff" %in% names(measures_i)) measures_i$deff <- NULL
    }
    if ("resp_rate" %in% varied) {
      if ("resp_rate" %in% names(frame_i)) frame_i$resp_rate <- NULL
      if ("resp_rate" %in% names(measures_i)) measures_i$resp_rate <- NULL
    }
    if ("eligibility_rate" %in% varied) {
      if ("eligibility_rate" %in% names(frame_i)) frame_i$eligibility_rate <- NULL
      if ("eligibility_rate" %in% names(measures_i)) measures_i$eligibility_rate <- NULL
    }
    if ("eligible_per_unit" %in% varied) {
      if ("eligible_per_unit" %in% names(frame_i)) frame_i$eligible_per_unit <- NULL
      if ("eligible_per_unit" %in% names(measures_i)) measures_i$eligible_per_unit <- NULL
    }
    if ("var_mult" %in% varied) {
      scaled_v <- .predict_scale_deff(deff_i, frame_i, measures_i, r$var_mult)
      deff_i <- scaled_v$deff
      frame_i <- scaled_v$frame
      measures_i <- scaled_v$measures
    }
    if ("cost_mult" %in% varied) {
      scaled_c <- .predict_scale_costs(frame_i, unit_cost_i, r$cost_mult)
      frame_i <- scaled_c$frame
      unit_cost_i <- scaled_c$unit_cost
    }
    if ("icc_mult" %in% varied) {
      measures_i <- .predict_scale_icc(measures_i, r$icc_mult)
      frame_i <- .predict_scale_icc(frame_i, r$icc_mult)
    }

    if (is_psu) {
      if (!is.null(object$psu) && !is.null(p$certainty$n_psu_draw)) {
        stratum <- as.character(frame_i[["stratum"]] %||% seq_len(nrow(frame_i)))
        idx_of <- split(seq_len(nrow(p$psu)),
                        factor(as.character(p$psu$stratum), levels = stratum))
        certain <- object$psu$certainty
        per_psu <- lapply(idx_of, function(i) object$psu[["take_2"]][i][certain[i]])
        n_psu_draw <- p$certainty$n_psu_draw
        take_ov <- if (isTRUE(p$free_take)) p[["take_2_int"]] else NULL
        res <- .psu_field_assess(
          frame = frame_i, psu = p$psu, measures = measures_i,
          targets = p$targets, unit_cost = unit_cost_i, alpha = r$alpha,
          deff = deff_i, resp_rate = r$resp_rate,
          min_n_stratum = p$min_n_stratum, objective = p$objective,
          budget = p$budget, df = r$df, certain = certain, per_psu = per_psu,
          n_psu_draw = n_psu_draw, eligibility_rate = r$eligibility_rate,
          take_override = take_ov
        )
        if (all(c("cost_1", "cost_2") %in% names(frame_i))) {
          n_cert <- vapply(idx_of, function(i) sum(certain[i]), numeric(1))
          take_use <- take_ov %||% .psu_take(frame_i, p$psu, stratum)
          n_int <- vapply(per_psu, sum, numeric(1)) + n_psu_draw * take_use
          res$params$achieved$cost <- sum(
            (n_cert + n_psu_draw) * as.numeric(frame_i$cost_1) +
              n_int * as.numeric(frame_i$cost_2)
          )
        }
      } else {
        held <- p$psu
        if (!is.null(object$psu)) held$certainty <- object$psu$certainty
        res <- .prec_alloc_psu(
          frame = frame_i, psu = held, n = p$n, measures = measures_i,
          targets = p$targets, objective = p$objective, budget = p$budget,
          unit_cost = unit_cost_i, alpha = r$alpha, deff = deff_i,
          resp_rate = r$resp_rate, min_n_stratum = p$min_n_stratum,
          df = r$df, certainty_cutoff = p$certainty_cutoff,
          n_psu_per_zone = p$n_psu_per_zone,
          eligibility_rate = r$eligibility_rate
        )
      }
    } else if (isTRUE(p$free_take)) {
      res <- .prec_alloc_freetake(
        frame = frame_i, n = n_h, take_2 = p$take_2_h,
        measures = measures_i, targets = p$targets,
        alpha = r$alpha, deff = deff_i, resp_rate = r$resp_rate,
        min_n_stratum = p$min_n_stratum, objective = p$objective,
        budget = p$budget, df = r$df, fpc = p$fpc %||% "unit",
        eligible_per_unit = r$eligible_per_unit,
        eligibility_rate = r$eligibility_rate,
        take_2_min = p[["take_2_min"]], take_2_max = p[["take_2_max"]], take_2_step = p[["take_2_step"]],
        objective_rule = p$objective_rule %||% "sum"
      )
    } else {
      res <- .prec_alloc_bethel(
        frame = frame_i, n = n_h, measures = measures_i,
        targets = p$targets, objective = p$objective, budget = p$budget,
        unit_cost = unit_cost_i, alpha = r$alpha, deff = deff_i,
        resp_rate = r$resp_rate, eligible_per_unit = r$eligible_per_unit,
        eligibility_rate = r$eligibility_rate,
        min_n_stratum = .resolve_min_n_stratum(p$min_n_stratum, frame_i),
        df = r$df, fpc = p$fpc %||% "unit",
        .allow_fractional_stages = !list_n,
        .finite_prop_var = TRUE,
        objective_rule = p$objective_rule %||% "sum"
      )
    }

    total <- res$params$achieved$n
    out <- data.frame(n = total, cost = res$params$achieved$cost)
    if (!is.null(res$objective_value)) {
      out$objective_value <- res$objective_value
      out$cv <- sqrt(res$objective_value)
    }
    if (nrow(res$detail) > 0L) {
      out$.max_ratio <- max(res$detail$.ratio, na.rm = TRUE)
      out$.n_fail <- sum(!res$detail$.pass, na.rm = TRUE)
      out$.feasible <- all(res$detail$.pass, na.rm = TRUE) &&
        all(res$bounds$.pass, na.rm = TRUE)
    } else {
      out$.max_ratio <- NA_real_
      out$.n_fail <- 0L
      out$.feasible <- all(res$bounds$.pass, na.rm = TRUE)
    }
    if (list_n) names(out)[1L] <- "total_n"
    out
  })
}

#' @keywords internal
#' @noRd
.validate_newdata <- function(newdata, allowed, allowed_character = character(0),
                              allowed_list = character(0)) {
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
  non_numeric <- setdiff(
    non_numeric,
    intersect(allowed_list, names(newdata)[vapply(newdata, is.list, logical(1))])
  )
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
      if (is.list(val)) val <- val[[1L]]
      params[[nm]] <- val
    }

    results[i] <- list(tryCatch(
      eval_fn(params),
      error = function(e) {
        warning(
          sprintf("predict row %d failed: %s", i, conditionMessage(e)),
          call. = FALSE
        )
        NULL
      }
    ))
  }

  ok <- !vapply(results, is.null, logical(1))
  if (!any(ok)) {
    stop("all rows failed evaluation", call. = FALSE)
  }

  template <- results[[which(ok)[1L]]]
  ## Indexing by NA keeps each column's type, logical and integer included.
  na_row <- template[NA_integer_, , drop = FALSE]

  result_df <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r)) na_row else r
  }))
  rownames(result_df) <- NULL

  dup_cols <- intersect(names(result_df), nd_names)
  achieved <- intersect(dup_cols, c("cv", "moe", "rmoe", "power"))
  result_df <- result_df[, !names(result_df) %in% setdiff(dup_cols, achieved),
                         drop = FALSE]
  hit <- names(result_df) %in% achieved
  names(result_df)[hit] <- paste0(".", names(result_df)[hit], "_achieved")

  out <- cbind(newdata, result_df)
  rownames(out) <- NULL
  out
}

#' Scale target columns (`cv`, `moe`, `rmoe`) by a positive multiplier
#' @keywords internal
#' @noRd
.predict_scale_targets <- function(tbl, mult) {
  if (is.null(tbl) || nrow(tbl) == 0L) return(tbl)
  if (!is.numeric(mult) || length(mult) != 1L || is.na(mult) ||
      !is.finite(mult) || mult <= 0) {
    stop("'target_mult' must be positive and finite", call. = FALSE)
  }
  if (mult == 1) return(tbl)
  for (col in intersect(c("cv", "moe", "rmoe"), names(tbl))) {
    idx <- !is.na(tbl[[col]])
    tbl[[col]][idx] <- tbl[[col]][idx] * mult
  }
  tbl
}

#' Scale ICC columns in a table by a non-negative multiplier
#' @keywords internal
#' @noRd
.predict_scale_icc <- function(tbl, mult) {
  if (is.null(tbl) || !is.data.frame(tbl)) return(tbl)
  if (!is.numeric(mult) || length(mult) != 1L || is.na(mult) ||
      !is.finite(mult) || mult < 0) {
    stop("'icc_mult' must be non-negative and finite", call. = FALSE)
  }
  if (mult == 1) return(tbl)
  cols <- grep("^icc(_[0-9]+)?$", names(tbl), value = TRUE)
  for (col in cols) {
    idx <- !is.na(tbl[[col]])
    scaled <- tbl[[col]][idx] * mult
    if (any(scaled > 1)) {
      stop(sprintf("scaled '%s' exceeds 1 under 'icc_mult = %g'", col, mult),
           call. = FALSE)
    }
    tbl[[col]][idx] <- scaled
  }
  tbl
}

#' Scale stratum/stage costs by a positive multiplier
#' @keywords internal
#' @noRd
.predict_scale_costs <- function(frame, unit_cost, mult) {
  if (!is.numeric(mult) || length(mult) != 1L || is.na(mult) ||
      !is.finite(mult) || mult <= 0) {
    stop("'cost_mult' must be positive and finite", call. = FALSE)
  }
  if (mult == 1) return(list(frame = frame, unit_cost = unit_cost))
  cost_cols <- grep("^cost_[0-9]+$", names(frame), value = TRUE)
  if (length(cost_cols) > 0L) {
    for (col in cost_cols) frame[[col]] <- frame[[col]] * mult
  } else if (!is.null(unit_cost)) {
    unit_cost <- unit_cost * mult
  } else if ("unit_cost" %in% names(frame)) {
    frame$unit_cost <- frame$unit_cost * mult
  } else {
    unit_cost <- mult
  }
  list(frame = frame, unit_cost = unit_cost)
}

#' Scale variance via design effect multiplier across scalar/frame/measures
#' @keywords internal
#' @noRd
.predict_scale_deff <- function(deff, frame, measures, mult) {
  if (!is.numeric(mult) || length(mult) != 1L || is.na(mult) ||
      !is.finite(mult) || mult <= 0) {
    stop("'var_mult' must be positive and finite", call. = FALSE)
  }
  if (mult == 1) return(list(deff = deff, frame = frame, measures = measures))
  deff <- deff * mult
  if ("deff" %in% names(frame)) {
    idx <- !is.na(frame$deff)
    frame$deff[idx] <- frame$deff[idx] * mult
  }
  if ("deff" %in% names(measures)) {
    idx <- !is.na(measures$deff)
    measures$deff[idx] <- measures$deff[idx] * mult
  }
  list(deff = deff, frame = frame, measures = measures)
}

#' Domain columns of a fit, refused when its domains cannot be rebuilt
#' @keywords internal
#' @noRd
.predict_domain_cols <- function(object) {
  cols <- object$params[["domain_cols"]]
  if (!is.null(object$domains) && length(cols) == 0L) {
    stop("this result has domains but records no domain columns, so ",
         "predict() cannot refit it", call. = FALSE)
  }
  cols
}

#' Scenario grid for multi-indicator cluster results (`n_cluster(indicators, ...)`)
#' @keywords internal
#' @noRd
.predict_cluster_multi <- function(object, newdata, control = NULL) {
  p <- object$params
  ind_base <- p[["indicators_input"]] %||% object$indicators
  allowed <- c("target_mult", "cost_mult", "icc_mult", "budget",
               "resp_rate", "eligibility_rate", "eligible_per_unit",
               "fixed_cost")
  .validate_newdata(newdata, allowed)
  varied <- names(newdata)
  dom_cols <- .predict_domain_cols(object)
  if ("icc_mult" %in% varied) {
    for (m in newdata$icc_mult) .predict_scale_icc(ind_base, m)
  }
  base <- list(
    target_mult = 1,
    cost_mult = 1,
    icc_mult = 1,
    budget = p$budget,
    resp_rate = p$resp_rate %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    fixed_cost = p$fixed_cost %||% 0
  )
  .predict_grid(newdata, base, function(r) {
    ind <- .predict_scale_targets(ind_base, r$target_mult)
    ind <- .predict_scale_icc(ind, r$icc_mult)
    if (!is.numeric(r$cost_mult) || is.na(r$cost_mult) ||
        !is.finite(r$cost_mult) || r$cost_mult <= 0) {
      stop("'cost_mult' must be positive and finite", call. = FALSE)
    }
    sc <- p$stage_cost * r$cost_mult
    if ("resp_rate" %in% varied && "resp_rate" %in% names(ind)) {
      ind$resp_rate <- NULL
    }
    if ("eligibility_rate" %in% varied && "eligibility_rate" %in% names(ind)) {
      ind$eligibility_rate <- NULL
    }
    if ("eligible_per_unit" %in% varied && "eligible_per_unit" %in% names(ind)) {
      ind$eligible_per_unit <- NULL
    }
    res <- tryCatch(n_cluster.data.frame(
      x = ind,
      stage_cost = sc,
      domains = dom_cols,
      budget = if ("budget" %in% varied) r$budget else p[["budget"]],
      stage_n = p[["stage_n"]],
      allocation = p[["allocation"]] %||% "separate",
      domain_sampling = p[["domain_sampling"]] %||% "separate",
      min_n_domain = p[["min_n_domain"]],
      fixed_cost = r$fixed_cost,
      resp_rate = .predict_final_rate(p$resp_rate %||% 1, r$resp_rate,
                                      "resp_rate" %in% varied),
      eligibility_rate = r$eligibility_rate,
      eligible_per_unit = r$eligible_per_unit,
      N = p$N,
      fpc = p$fpc %||% "auto",
      control = control
    ), svyplan_plan_unknown = function(e) NULL)
    if (is.null(res)) {
      return(.predict_cluster_unknown(names(object$n), c("total_n", "cost"),
                                      text = "binding"))
    }
    out <- if (!is.null(res$n) && !is.null(names(res$n))) {
      as.list(res$n)
    } else {
      list()
    }
    out$total_n <- res$total_n
    out$cost <- res$cost
    out$binding <- paste(res$binding, collapse = ", ")
    out$.status <- res$operational$search$status %||% NA_character_
    out$.nodes <- res$operational$search$nodes %||% NA_real_
    as.data.frame(out, stringsAsFactors = FALSE)
  })
}

#' Grid exploration for joint constrained allocations
#' @keywords internal
#' @noRd
.predict_bethel <- function(object, newdata) {
  p <- object$params
  budget_mode <- identical(p$mode, "budget_objective")
  stages <- p$stages %||% 1L
  # A certainty fit is one-stage to the solver, so `stages` cannot find the
  # take, but it is still a frame column and still worth sweeping.
  certainty <- !is.null(object$params$psu)
  free_take <- isTRUE(p$free_take)
  take_cols <- if (certainty) {
    "take_2"
  } else if (stages >= 2L) {
    paste0("take_", 2L:stages)
  } else {
    character(0)
  }
  targets <- if (is.null(p$targets) || nrow(p$targets) == 0L) NULL else
    p$targets
  allowed <- c(
    if (budget_mode) "budget",
    take_cols,
    if (free_take) c("take_2_min", "take_2_max"),
    if (certainty) c("certainty_cutoff", "n_psu_per_zone"),
    if (!is.null(targets)) "target_mult",
    "cost_mult",
    if (stages >= 2L || certainty) "icc_mult",
    "deff", "resp_rate", "eligibility_rate",
    if (!certainty) "eligible_per_unit", "alpha", "df"
  )
  .validate_newdata(newdata, allowed)
  varied <- names(newdata)

  varied_takes <- intersect(take_cols, varied)
  for (nm in varied_takes) {
    value <- newdata[[nm]]
    if (anyNA(value) || any(!is.finite(value)) || any(value < 1) ||
          any(abs(value - round(value)) > 1e-8)) {
      stop(sprintf("newdata '%s' must contain positive whole numbers", nm),
           call. = FALSE)
    }
  }
  for (nm in intersect(c("take_2_min", "take_2_max"), varied)) {
    value <- newdata[[nm]]
    if (anyNA(value) || any(!is.finite(value)) || any(value < 1) ||
          any(abs(value - round(value)) > 1e-8)) {
      stop(sprintf("newdata '%s' must contain positive whole numbers", nm),
           call. = FALSE)
    }
  }
  if ("target_mult" %in% varied) {
    value <- newdata$target_mult
    if (anyNA(value) || any(!is.finite(value)) || any(value <= 0)) {
      stop("newdata 'target_mult' must contain positive finite values",
           call. = FALSE)
    }
  }
  if ("cost_mult" %in% varied) {
    value <- newdata$cost_mult
    if (anyNA(value) || any(!is.finite(value)) || any(value <= 0)) {
      stop("newdata 'cost_mult' must contain positive finite values",
           call. = FALSE)
    }
  }
  if ("icc_mult" %in% varied) {
    value <- newdata$icc_mult
    if (anyNA(value) || any(!is.finite(value)) || any(value < 0)) {
      stop("newdata 'icc_mult' must contain non-negative finite values",
           call. = FALSE)
    }
    for (m in value) {
      .predict_scale_icc(p$measures, m)
      .predict_scale_icc(p$frame, m)
    }
  }
  # One value per row applies to every stratum. A per-stratum cutoff is a
  # frame column, as for any other stratum input.
  if ("certainty_cutoff" %in% varied) {
    value <- newdata$certainty_cutoff
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
          any(value <= 0 | value > 1)) {
      stop("newdata 'certainty_cutoff' must contain values in (0, 1]",
           call. = FALSE)
    }
  }

  # NA is the design without zones, which NULL means to n_alloc().
  if ("n_psu_per_zone" %in% varied) {
    value <- newdata$n_psu_per_zone
    if (any(!is.na(value) & !value %in% c(1, 2))) {
      stop("newdata 'n_psu_per_zone' must contain NA, 1 or 2", call. = FALSE)
    }
  }

  zoned <- certainty &&
    ("n_psu_per_zone" %in% varied || !is.null(p$n_psu_per_zone))
  cols <- .bethel_predict_cols(budget_mode, stages, certainty, zoned)

  rows <- lapply(seq_len(nrow(newdata)), function(i) {
    frame <- p$frame
    measures <- p$measures
    targets_i <- targets
    unit_cost_i <- p$unit_cost
    deff_i <- if ("deff" %in% varied) newdata$deff[i] else p$deff
    resp_i <- if ("resp_rate" %in% varied) newdata$resp_rate[i] else p$resp_rate
    elig_r_i <- if ("eligibility_rate" %in% varied) {
      newdata$eligibility_rate[i]
    } else {
      p[["eligibility_rate"]] %||% 1
    }
    epu_i <- if ("eligible_per_unit" %in% varied) {
      newdata$eligible_per_unit[i]
    } else {
      p$eligible_per_unit %||% 1
    }
    alpha_i <- if ("alpha" %in% varied) newdata$alpha[i] else p$alpha
    df_i <- if ("df" %in% varied) newdata$df[i] else p$df

    for (nm in varied_takes) frame[[nm]] <- round(newdata[[nm]][i])
    if ("target_mult" %in% varied) {
      targets_i <- .predict_scale_targets(targets_i, newdata$target_mult[i])
    }
    if ("cost_mult" %in% varied) {
      sc <- .predict_scale_costs(frame, unit_cost_i, newdata$cost_mult[i])
      frame <- sc$frame
      unit_cost_i <- sc$unit_cost
    }
    if ("icc_mult" %in% varied) {
      measures <- .predict_scale_icc(measures, newdata$icc_mult[i])
      frame <- .predict_scale_icc(frame, newdata$icc_mult[i])
    }
    if ("deff" %in% varied) {
      if ("deff" %in% names(frame)) frame$deff <- NULL
      if ("deff" %in% names(measures)) measures$deff <- NULL
    }
    if ("resp_rate" %in% varied) {
      if ("resp_rate" %in% names(frame)) frame$resp_rate <- NULL
      if ("resp_rate" %in% names(measures)) measures$resp_rate <- NULL
    }
    if ("eligibility_rate" %in% varied) {
      if ("eligibility_rate" %in% names(frame)) frame$eligibility_rate <- NULL
      if ("eligibility_rate" %in% names(measures)) measures$eligibility_rate <- NULL
    }
    if ("eligible_per_unit" %in% varied) {
      if ("eligible_per_unit" %in% names(frame)) frame$eligible_per_unit <- NULL
      if ("eligible_per_unit" %in% names(measures)) {
        measures$eligible_per_unit <- NULL
      }
    }
    if ("alpha" %in% varied && !is.null(targets_i) && "alpha" %in% names(targets_i)) {
      targets_i$alpha <- NA_real_
    }
    if ("df" %in% varied && !is.null(targets_i) && "df" %in% names(targets_i)) {
      targets_i$df <- NA_real_
    }

    # If fixed take_2 is varied on a free-take fit, drop take bounds so fixed-take runs
    use_ft <- free_take && !("take_2" %in% varied_takes)
    t_min_i <- if (!use_ft) NULL else if ("take_2_min" %in% varied) {
      round(newdata[["take_2_min"]][i])
    } else {
      p[["take_2_min"]]
    }
    t_max_i <- if (!use_ft) NULL else if ("take_2_max" %in% varied) {
      round(newdata[["take_2_max"]][i])
    } else {
      p[["take_2_max"]]
    }
    t_step_i <- if (!use_ft) NULL else p[["take_2_step"]]
    if (!use_ft) {
      frame <- frame[, setdiff(names(frame),
                               c("take_2_min", "take_2_max", "take_2_step")),
                     drop = FALSE]
    }

    budget <- if (!budget_mode) {
      NULL
    } else if ("budget" %in% varied) {
      newdata$budget[i]
    } else {
      p$budget
    }
    cutoff <- if ("certainty_cutoff" %in% varied) {
      newdata$certainty_cutoff[i]
    } else {
      p$certainty_cutoff
    }
    zone_m <- if ("n_psu_per_zone" %in% varied) {
      if (is.na(newdata$n_psu_per_zone[i])) NULL else newdata$n_psu_per_zone[i]
    } else {
      p$n_psu_per_zone
    }
    fit <- tryCatch(
      n_alloc_joint.default(
        frame = frame,
        measures = measures,
        targets = targets_i,
        psu = p$psu,
        certainty_cutoff = cutoff,
        n_psu_per_zone = zone_m,
        control = p[["control"]],
        objective = p$objective,
        budget = budget,
        objective_rule = p$objective_rule %||% "sum",
        unit_cost = unit_cost_i,
        alpha = alpha_i,
        deff = deff_i,
        resp_rate = resp_i,
        eligibility_rate = elig_r_i,
        eligible_per_unit = epu_i,
        df = df_i,
        min_n_stratum = p$min_n_stratum,
        take_2_min = t_min_i,
        take_2_max = t_max_i,
        take_2_step = t_step_i,
        fpc = p$fpc %||% "unit"
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      status <- if (certainty) .plan_error_status(fit) else "infeasible"
      return(list(values = .bethel_predict_row(NULL, cols, status),
                  status = status, message = conditionMessage(fit)))
    }
    list(values = .bethel_predict_row(fit, cols))
  })
  values <- do.call(rbind, lapply(rows, `[[`, "values"))
  if (certainty && budget_mode) {
    values <- .predict_budget_monotone(newdata, values)
  }
  for (i in seq_along(rows)) {
    if (is.null(rows[[i]]$message) || isTRUE(values$.feasible[i])) next
    label <- .predict_row_label(newdata[i, , drop = FALSE])
    warning(
      if (identical(rows[[i]]$status, "unknown")) {
        sprintf("feasibility of %s is unknown: %s", label, rows[[i]]$message)
      } else {
        sprintf("%s is infeasible: %s", label, rows[[i]]$message)
      },
      call. = FALSE
    )
  }
  out <- cbind(newdata, values)
  rownames(out) <- NULL
  out
}

#' Whether a failed grid row is proved infeasible or left unknown
#'
#' A node limit reached with no plan, an internal error and a solver that
#' fails its certificate prove nothing about the design point.
#' @keywords internal
#' @noRd
.plan_error_status <- function(e) {
  unknown <- inherits(e, "svyplan_plan_unknown") ||
    grepl("^internal error|failed certification", conditionMessage(e))
  if (unknown) "unknown" else "infeasible"
}

#' A register budget grid whose objective never worsens with the budget
#'
#' Rows equal but for the budget are read in budget order. A plan fits
#' every larger budget, so a row stopped at the node limit with a worse
#' plan, or left without one, takes the best plan of a lower budget. Its
#' status reads `"limit"`, as no search proved it at this budget. The
#' continuous columns stay the row's own.
#' @keywords internal
#' @noRd
.predict_budget_monotone <- function(newdata, values) {
  if (!"budget" %in% names(newdata)) return(values)
  other <- setdiff(names(newdata), "budget")
  key <- if (length(other)) {
    do.call(paste, c(lapply(newdata[other], format), sep = "\r"))
  } else {
    rep("", nrow(newdata))
  }
  field <- intersect(
    c("n_psu_certain", "n_psu_draw", "n_zone", "objective_value", "cv",
      "n_1_int", "n_int", "cost_int", ".status", ".feasible"),
    names(values)
  )
  planned <- function(r) values$.status[r] %in% c("optimal", "gap", "limit")
  for (k in unique(key)) {
    rows <- which(key == k)
    best <- NA_integer_
    for (r in rows[order(newdata$budget[rows])]) {
      if (!is.na(best) && values$.status[r] %in% c("gap", "limit", "unknown") &&
            (!planned(r) ||
               values$objective_value[best] < values$objective_value[r])) {
        values[r, field] <- values[best, field]
        values$.status[r] <- "limit"
      }
      if (planned(r) && (is.na(best) ||
                           values$objective_value[r] <
                             values$objective_value[best])) {
        best <- r
      }
    }
  }
  values
}

#' Result columns of a joint-allocation grid, in report order
#' @keywords internal
#' @noRd
.bethel_predict_cols <- function(budget_mode, stages, certainty = FALSE,
                                 zoned = FALSE) {
  if (certainty) {
    return(c(
      "n_psu_certain", "n_psu_draw", if (zoned) "n_zone", "n", "cost",
      if (budget_mode) c("objective_value", "cv"),
      "n_1_int", "n_int", "cost_int", if (budget_mode) ".binding", ".status",
      ".feasible"
    ))
  }
  c(
    if (stages > 1L) "n_1",
    "n", "cost",
    if (budget_mode) c("objective_value", "cv"),
    if (stages > 1L) "n_1_int",
    "n_int", "cost_int",
    if (budget_mode) ".binding",
    ".feasible"
  )
}

#' One grid row, or its all-NA counterpart when the solve failed
#' @keywords internal
#' @noRd
.bethel_predict_row <- function(fit, cols, status = "infeasible") {
  if (is.null(fit)) {
    values <- as.list(rep(NA_real_, length(cols)))
    names(values) <- cols
    if (".binding" %in% cols) values$.binding <- NA
    if (".status" %in% cols) values$.status <- status
    values$.feasible <- if (identical(status, "unknown")) NA else FALSE
    return(as.data.frame(values, stringsAsFactors = FALSE))
  }
  # A register fit plans its field design, so its objective is that one's.
  objective <- if (is.null(fit$params$psu)) {
    fit$objective_value
  } else {
    fit$operational$objective_value
  }
  value_of <- function(name) {
    switch(
      name,
      n_1 = sum(fit$detail$n_1),
      n_psu_certain = sum(fit$detail$n_psu_certain),
      n_psu_draw = sum(fit$detail$n_psu_draw),
      n_zone = sum(fit$detail$n_zone %||% 0),
      n = fit$n,
      cost = fit$params$achieved$cost,
      objective_value = objective,
      cv = sqrt(objective),
      n_1_int = sum(fit$detail$n_1_int),
      n_int = as.numeric(fit$operational$n),
      cost_int = fit$operational$cost,
      .binding = isTRUE(fit$optimization$budget_binding),
      .status = fit$operational$status,
      .feasible = isTRUE(fit$operational$all_pass)
    )
  }
  values <- lapply(cols, value_of)
  names(values) <- cols
  as.data.frame(values, stringsAsFactors = FALSE)
}

#' Name the grid row a warning is about
#' @keywords internal
#' @noRd
.predict_row_label <- function(row) {
  paste(
    sprintf("%s = %s", names(row), vapply(row, format, character(1))),
    collapse = ", "
  )
}
