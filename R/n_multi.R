#' Multi-indicator sample size
#'
#' Compute the sample size that satisfies precision requirements for
#' multiple survey indicators simultaneously under a simple sampling design.
#' Optional domain columns support separate requirements by subpopulation.
#'
#' @param indicators For the default method: a data frame where **each row
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
#'     \item **How precise**: `moe` (margin of error), `rmoe` (margin of
#'       error relative to the estimand) **or** `cv` (coefficient of
#'       variation). Each row must specify exactly one.
#'   }
#'
#'   For `svyplan_prec` objects: a precision result from [prec_multi()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param domains Character vector of column names in `indicators` to treat
#'   as domain variables, or `NULL` (default) for no domains. All names
#'   must exist in `indicators`. When specified, sizing runs
#'   independently for each domain combination.
#' @param domain_sampling How the per-domain requirements combine into one
#'   overall size, either `"separate"` (default) or `"natural"`. It applies
#'   only when domains are present.
#'
#'   `"separate"` treats the domains as disjoint quotas drawn independently,
#'   the usual case for regional or urban/rural domains in one national
#'   survey, and reports `n` as their sum. Every quota has to be met, so no
#'   smaller total delivers the design.
#'
#'   `"natural"` treats them as analytic domains that arise at their own rate
#'   inside one sample. It requires a `share` column giving each domain's
#'   expected share of the population, constant within a domain, and reports
#'   `n` as `max(.n / share)`. That is an *expected yield*: a sample of that
#'   size delivers each domain's quota on average, not with certainty in any
#'   one realized sample.
#' @param min_n_domain Numeric scalar or `NULL` (default). Minimum total sample
#'   size per domain. It applies only when domains are present. Per-domain
#'   sample sizes are floored to `min_n_domain`.
#' @param prop_method Proportion CI method, one of `"wald"`
#'   (default), `"wilson"`, `"logodds"`, or `"beta"`. This is passed to
#'   [n_prop()] for proportion rows and ignored for mean rows.
#'   An optional `prop_method` column in `indicators` overrides this default
#'   on a per-row basis. `"wilson"`, `"logodds"` and `"beta"` size from an
#'   interval half-width, so those rows need `moe` rather than `cv`.
#' @param resp_rate Default expected response rate at the ultimate unit, in
#'   (0, 1\]. Used for rows whose `resp_rate` column is absent or `NA`, and
#'   a non-missing row value overrides it.
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
#'     \item{`indicators`}{The input indicators data frame.}
#'   }
#'
#'   **With domains**, the object additionally contains:
#'   \describe{
#'     \item{`n`}{The overall size, read according to `domain_sampling`:
#'       the sum of the domain quotas under `"separate"`, or the size whose
#'       expected yield meets every quota under `"natural"`. It is not the
#'       largest domain requirement, which is a different and smaller
#'       number.}
#'     \item{`n_domain_max`}{The largest single domain requirement, and the
#'       one `binding` refers to.}
#'     \item{`domains`}{Data frame with one row per domain, including
#'       domain variables, `.n`, `.binding`, and `.share` under
#'       `"natural"`.}
#'   }
#'
#' @details
#' ## Building the indicators data frame
#'
#' Each row of `indicators` represents one survey indicator. The two key
#' decisions per row are:
#'
#' 1. **Type of indicator**: is it a proportion (binary variable like
#'    "stunted yes/no") or a mean (continuous variable like "household
#'    expenditure")? This determines whether you fill the `p` or `var`
#'    column.
#' 2. **Precision target**: do you want an absolute margin of error
#'    (`moe`, e.g. +/- 5 percentage points), the same margin stated
#'    relative to the estimand (`rmoe`, e.g. 12 percent of the
#'    proportion, as MICS and DHS state it), or a relative coefficient of
#'    variation (`cv`, e.g. 10 percent relative error)?
#'
#' A minimal example for three health indicators:
#'
#' ```
#' indicators <- data.frame(
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
#'   \item{`mu`}{Population mean, finite and non-zero. It may be
#'     negative. It is required when `var` is used with `cv` because
#'     CV = SE / abs(mean).}
#'   \item{`moe`}{Margin of error, the half-width of the confidence
#'     interval you want. For proportions, this is on the probability
#'     scale (e.g. 0.05 for +/- 5 percentage points). For means,
#'     it is in the same units as the variable (e.g. 10 dollars).}
#'   \item{`rmoe`}{Margin of error relative to the row's estimand, so
#'     0.12 asks for a half-width of 12 percent of `p` or of `abs(mu)`.
#'     It is converted to `moe` on ingestion, so the row needs `p` or
#'     `mu`, and it is read under the row's own `prop_method`.}
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
#'     `"wald"` (default), `"wilson"`, `"logodds"`, or `"beta"`.
#'     `"wilson"` is recommended for rare proportions (below 0.1 or above
#'     0.9), and `"beta"` (Korn-Graubard) when the expected number of
#'     positive counts is small enough that the interval should stay
#'     inside `[0, 1]` by construction. It is used only for rows with
#'     `p`. See [n_prop()] for how to choose.}
#'   \item{`df`}{Degrees of freedom of the planned variance estimator,
#'     typically sampled PSUs minus strata, and available from
#'     [design_df()]. It switches that row's interval quantile from normal
#'     to t, under every method. `NA` (the default) applies no adjustment.
#'     A df is a property of the design rather than of an indicator, so
#'     every row of a single-domain table shares one value. The column
#'     earns its place when the rows are *domains*, each covering its own
#'     set of strata: `design_df(alloc)$domains$.df` gives one value per
#'     domain to match in.}
#'   \item{`min_cases`}{Minimum expected number of positive cases the row
#'     must yield, a floor on that row's own size, as in
#'     [n_prop()]`(min_cases = )`. Proportion rows only, and `NA` (the
#'     default) sizes the row on precision alone. The row that ends up
#'     largest is still the binding one, whichever constraint raised it.
#'     Single-stage `n_multi()` only: a multistage design sizes stages
#'     against a cost and has no single total for a count to raise.}
#'   \item{`unit_relvar`}{Unit relvariance. If omitted, derived
#'     automatically from `p` (as `(1 - p) / p`) or from
#'     `var` / `mu^2`.}
#'   \item{`resp_rate`}{Expected response rate at the ultimate unit, in
#'     (0, 1\]. Default 1 (no adjustment). A value of 0.90 inflates the
#'     sample size by `1 / 0.90` to compensate for 10 percent
#'     non-response. It means the same thing in [n_multi_cluster()],
#'     which additionally takes `resp_rate_psu` for whole clusters that
#'     cannot be worked and `resp_rate_ssu` for second-stage units in a
#'     three-stage design. Naming a stage the design does not have is an
#'     error rather than a column carried along and ignored.}
#' }
#'
#' Domain columns are specified via the `domains` parameter. When domains
#' are present, sizing runs independently for each domain combination, and
#' `domain_sampling` decides how those requirements combine into the one
#' number `n` reports. The two readings answer different questions and give
#' different sizes, so the choice belongs to the design rather than to a
#' default: separate quotas need every requirement met and therefore their
#' sum, while natural domains need a sample large enough that the rarest
#' demanding domain turns up often enough. Neither is the largest single
#' requirement, which `n_domain_max` reports separately.
#'
#' Columns beyond those listed are carried along and ignored, so the table
#' can keep questionnaire modules, sources, or other bookkeeping. The
#' exception is a name the rest of the package would lead you to expect
#' here: `mean` (this table takes `mu`, while the [n_alloc()] frame takes
#' `mean`), `indicator` or `label` (it takes `name`), `method` (it takes
#' `prop_method`), and `icc` or `var_ratio` (which are per stage here, so
#' `icc_psu` and `var_ratio_psu`). Those are rejected with the name they
#' should have carried rather than silently ignored. A column named
#' through `domains` is exempt, so a domain may be called `label` or
#' `method`.
#'
#' The dispersion of a continuous indicator may be given as `var` or as
#' `sd`, whichever the source reports, exactly as in [n_mean()]. They are
#' different quantities rather than two names for one, so supplying both
#' in the same table is an error rather than a preference. `sd` is squared
#' on the way in and everything downstream reads `var`.
#'
#' `n_multi()` computes sample size per indicator by delegating proportion
#' rows to [n_prop()] and mean rows to [n_mean()],
#' then takes the maximum per domain. Use `prop_method` or a
#' `indicators$prop_method` column to choose `"wald"`, `"wilson"`,
#' `"logodds"` or `"beta"` for proportion rows.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @family sample size functions
#' @seealso [n_prop()] and [n_mean()] for single-indicator sizing,
#'   [n_alloc()] to split a multi-indicator size across strata or domains,
#'   [n_multi_cluster()] for multistage cluster designs, and [prec_multi()]
#'   for the inverse.
#'
#' @examples
#' # Simple mode: three indicators, take the max
#' indicators <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   moe  = c(0.05, 0.05, 0.03)
#' )
#' n_multi(indicators)
#'
#' # MICS/DHS-style: state precision as a relative margin of error
#' targets_rmoe <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   rmoe = 0.12,
#'   deff = c(2.0, 1.5, 2.5)
#' )
#' n_multi(targets_rmoe)
#'
#' # Continuous indicators: 'var' for the dispersion, 'mu' for the mean.
#' # A CV target needs 'mu', because a relative standard error is relative
#' # to something; a MOE target does not.
#' targets_mean <- data.frame(
#'   name = c("expenditure", "hh_size"),
#'   var  = c(250000, 4.0),
#'   mu   = c(1200, 5.4),
#'   cv   = c(0.05, 0.03)
#' )
#' n_multi(targets_mean)
#'
#' # Proportions and means in one table, sized to a common CV
#' targets_both <- data.frame(
#'   name = c("stunting", "expenditure"),
#'   p    = c(0.30, NA),
#'   var  = c(NA, 250000),
#'   mu   = c(NA, 1200),
#'   cv   = c(0.08, 0.05)
#' )
#' n_multi(targets_both)
#'
#' # Rare proportion: use Wilson globally in simple mode
#' n_multi(indicators[3, , drop = FALSE], prop_method = "wilson")
#'
#' # Korn-Graubard for an indicator whose expected count is small, with the
#' # degrees of freedom of the planned variance estimator. 'df' is read on
#' # beta rows only, so the other rows carry NA.
#' targets_rare <- data.frame(
#'   name = c("stunting", "cocaine_use"),
#'   p    = c(0.30, 0.02),
#'   moe  = c(0.05, 0.01),
#'   prop_method = c("wald", "beta"),
#'   df   = c(NA, 25)
#' )
#' n_multi(targets_rare)
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
#'   icc_psu = c(0.02, 0.05)
#' )
#' n_multi_cluster(targets_cl, stage_cost = c(500, 50))
#'
#' # Two-stage with MOE (converted to CV internally)
#' targets_moe <- data.frame(
#'   name   = c("stunting", "anemia"),
#'   p      = c(0.30, 0.10),
#'   moe    = c(0.05, 0.03),
#'   icc_psu = c(0.02, 0.05)
#' )
#' n_multi_cluster(targets_moe, stage_cost = c(500, 50))
#'
#' # Joint budget allocation across domains
#' targets_jnt <- data.frame(
#'   name   = rep(c("stunting", "anemia"), each = 2),
#'   p      = c(0.30, 0.25, 0.10, 0.15),
#'   cv     = c(0.10, 0.10, 0.15, 0.15),
#'   icc_psu = c(0.02, 0.03, 0.05, 0.04),
#'   region = rep(c("Urban", "Rural"), 2)
#' )
#' n_multi_cluster(
#'   targets_jnt,
#'   stage_cost = c(500, 50),
#'   domains = "region",
#'   budget = 100000,
#'   allocation = "joint"
#' )
#'
#' @export
n_multi <- function(indicators, ...) {
  if (!missing(indicators)) {
    .res <- .dispatch_plan(indicators, "indicators", n_multi.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_multi")
}

#' @rdname n_multi
#' @export
n_multi.default <- function(
  indicators,
  ...,
  domains = NULL,
  domain_sampling = c("separate", "natural"),
  min_n_domain = NULL,
  prop_method = c("wald", "wilson", "logodds", "beta"),
  resp_rate = 1,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_multi.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_multi.default, c(.plan, list(...))))
  .check_multi_split_args(list(...), "n_multi_cluster()")
  .check_unused_dots(...)
  domain_sampling <- match.arg(domain_sampling)
  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }
  if (!is.null(min_n_domain)) {
    if (
      !is.numeric(min_n_domain) || length(min_n_domain) != 1L || is.na(min_n_domain) || min_n_domain <= 0
    ) {
      stop("'min_n_domain' must be a positive numeric scalar", call. = FALSE)
    }
  }
  check_resp_rate(resp_rate)
  # An unresolved default arrives either as a missing argument or, through
  # the plan-merge path, as the full choice vector itself.
  if (identical(prop_method, c("wald", "wilson", "logodds", "beta"))) {
    prop_method <- "wald"
  }
  if (
    !is.character(prop_method) ||
      length(prop_method) != 1L ||
      is.na(prop_method) ||
      !prop_method %in% c("wald", "wilson", "logodds", "beta")
  ) {
    stop(
      "'prop_method' must be one of 'wald', 'wilson', 'logodds', or 'beta'",
      call. = FALSE
    )
  }

  indicators <- .indicators_var_from_sd(indicators, domains)
  indicators <- .indicators_moe_from_rmoe(indicators, domains)
  info <- .validate_targets(indicators, FALSE, domains = domains)
  indicators <- .fill_defaults(
    indicators,
    FALSE,
    prop_method = prop_method,
    resp_rate = resp_rate
  )

  rv_final <- indicators$unit_relvar[!is.na(indicators$unit_relvar)]
  if (
    length(rv_final) > 0L &&
      (any(rv_final <= 0) || any(!is.finite(rv_final)))
  ) {
    stop("'unit_relvar' values must be positive and finite", call. = FALSE)
  }
  domain_cols <- info$domain_cols
  mode <- if ("moe" %in% names(indicators) && any(!is.na(indicators$moe))) {
    "moe"
  } else {
    "cv"
  }

  if (length(domain_cols) == 0L) {
    .n_multi_simple(indicators, domain_cols = domain_cols, mode = mode,
                    prop_method = prop_method)
  } else {
    .n_multi_domains(
      indicators,
      stage_cost = NULL,
      budget = NULL,
      n_psu = NULL,
      n_per_psu = NULL,
      n_per_ssu = NULL,
      domain_cols,
      multistage = FALSE,
      joint = FALSE,
      min_n_domain,
      fixed_cost = 0,
      mode = mode,
      prop_method = prop_method,
      domain_sampling = domain_sampling
    )
  }
}

#' Multi-indicator sample size for cluster designs
#'
#' Compute a two- or three-stage cluster allocation that satisfies precision
#' requirements for several survey indicators. Domain-level planning and a
#' shared budget across domains are supported.
#'
#' @param indicators For the default method, a non-empty data frame with one row
#'   per indicator. Each row requires `p` or `var`, a `cv`, `moe`, or
#'   `rmoe` target, and `icc_psu`. Three-stage designs also require
#'   `icc_ssu`. Optional
#'   `var_ratio_psu` defaults to 1; three-stage `var_ratio_ssu` is derived as
#'   `var_ratio_psu * (1 - icc_psu)` when absent, the value the variance
#'   decomposition implies (see [design_effect()]). For the `svyplan_prec`
#'   method, a result from [prec_multi_cluster()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param stage_cost Numeric vector of per-stage costs with length 2 or 3.
#' @param domains Optional character vector naming domain columns in
#'   `indicators`. The function solves each domain independently unless
#'   `allocation = "joint"` in budget mode. Domains are sized as separate
#'   quotas, so `$total_n` is their sum and `$n` carries no aggregate stage
#'   vector; the fieldable per-domain stage sizes are in `$domains`.
#'   `domain_sampling = "natural"` is refused here, because a domain's
#'   expected yield in a multistage design depends on how its members sit
#'   inside PSUs and SSUs rather than on its population share alone.
#' @param budget Optional total budget. Supply precision indicators or a budget,
#'   according to the target schema described in Details.
#' @param n_psu Optional fixed stage-1 sample size.
#' @param n_per_psu Optional fixed stage-2 sample size per PSU.
#' @param n_per_ssu Optional fixed stage-3 sample size per SSU. This is valid
#'   only for three-stage designs.
#' @param allocation How a budget is split across domains, either
#'   `"separate"` (default, each domain sized on its own) or `"joint"`
#'   (one budget split across domains to minimize the worst precision
#'   ratio). `"joint"` applies only when `domains` and `budget` are
#'   supplied.
#' @param domain_sampling How the per-domain requirements combine into one
#'   overall size. Only `"separate"` (the default) is available here, and
#'   `$total_n` is then the sum of the per-domain totals. `"natural"` is
#'   refused: a domain's expected yield in a multistage design depends on
#'   how its members sit inside PSUs and SSUs rather than on its share of
#'   the population alone, and that model is not in the package.
#' @param min_n_domain Optional positive minimum total sample size per domain. In
#'   joint budget mode it is a constraint. In independent domain mode,
#'   domains below the floor produce a warning.
#' @param fixed_cost Non-negative fixed overhead cost. The default is 0.
#' @param resp_rate_psu Default expected PSU response rate, in (0, 1\].
#'   Used where the indicator column is absent or `NA`.
#' @param resp_rate_ssu Default expected SSU response rate for a three-stage
#'   design, in (0, 1\]. It is not applicable to a two-stage design.
#' @param resp_rate Default expected ultimate-unit response rate, in (0, 1\].
#'   Non-missing indicator columns override these three defaults row by row.
#' @param plan Optional [svyplan()] profile providing `stage_cost` and other
#'   applicable defaults.
#'
#' @return A `svyplan_cluster` object. The output class does not depend on
#'   which optional arguments are supplied.
#'
#' @details
#' The indicator columns follow [n_multi()], with one difference that matters.
#' Nonresponse is named for the stage it acts on: `resp_rate_psu` for
#' clusters that cannot be worked at all, `resp_rate_ssu` for second-stage
#' units in a three-stage design, and `resp_rate` for the ultimate units.
#' They are not interchangeable, and [n_cluster()] sets out why. A column
#' naming a stage the design does not have is an error rather than a column
#' carried along and ignored, since a silently dropped response rate plans a
#' design with none.
#'
#' Margin-of-error indicators are converted to CV before optimization. For each
#' candidate allocation, the required stage-1 size is the maximum across all
#' indicators. The solver minimizes total cost for precision indicators or the
#' worst precision ratio under a fixed budget.
#'
#' That conversion respects the row's `prop_method`. The multistage model is
#' driven by a relative standard error, and only the Wald interval has
#' half-width `z * se`, so a proportion row's `moe` is restated as the
#' sampling CV at the effective sample size its own method needs to close the
#' interval to that margin. A stricter interval therefore asks for a larger
#' design, in the same order it does in [n_prop()]. Under `"wald"` the
#' restatement is `moe / (z * p)` exactly. A mean row converts as
#' `moe / (z * |mu|)`, on the magnitude so that a negative mean yields a
#' positive target. [prec_multi_cluster()] inverts the same way, so a design
#' sized from a `moe` target reports that `moe` back.
#'
#' Homogeneity values numerically close to 0 or 1 are rejected because they
#' make the analytical cluster optimum degenerate. The result includes an
#' integer `$operational` allocation that preserves the applicable precision
#' or budget constraint. See [n_multi()] for shared indicator columns and
#' [n_cluster()] for the cluster cost model.
#'
#' ## How strong the optimum is
#'
#' With a stage size fixed, the remaining problem is solved from the
#' closed-form cluster optimum. With all stage sizes free, the objective is
#' a maximum over indicator requirements, which is not smooth, and it is
#' minimized by a bounded quasi-Newton search. That search warns when it
#' fails to converge or lands on a bound, but it carries no
#' global-optimality or KKT certificate: a successful return means the best
#' design this search found, not a proven minimum-cost one. The
#' `$operational` allocation can always be checked against its own
#' precision or budget constraint, which is a separate and exact statement.
#' [n_alloc()] gives the stronger guarantee where it applies, returning
#' feasibility and KKT diagnostics for its convex continuous problem.
#'
#' @family sample size functions
#' @seealso [n_multi()] for simple designs, [n_cluster()] for a single
#'   indicator, [n_alloc()] for stratified multistage allocation, and
#'   [prec_multi_cluster()] for the inverse calculation.
#'
#' @examples
#' indicators <- data.frame(
#'   name = c("stunting", "anemia"),
#'   p = c(0.30, 0.10),
#'   cv = c(0.10, 0.15),
#'   icc_psu = c(0.02, 0.05)
#' )
#' n_multi_cluster(indicators, stage_cost = c(500, 50))
#'
#' @export
n_multi_cluster <- function(indicators, ...) {
  if (!missing(indicators)) {
    .res <- .dispatch_plan(indicators, "indicators", n_multi_cluster.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_multi_cluster")
}

#' @rdname n_multi_cluster
#' @export
n_multi_cluster.default <- function(
  indicators,
  ...,
  stage_cost = NULL,
  domains = NULL,
  budget = NULL,
  n_psu = NULL,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  allocation = c("separate", "joint"),
  domain_sampling = c("separate", "natural"),
  min_n_domain = NULL,
  fixed_cost = 0,
  resp_rate_psu = 1,
  resp_rate_ssu = 1,
  resp_rate = 1,
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

  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }
  if (is.null(stage_cost)) {
    stop("'stage_cost' is required (directly or via plan)", call. = FALSE)
  }
  if (!is.character(allocation) || anyNA(allocation)) {
    stop("'allocation' must be \"separate\" or \"joint\"", call. = FALSE)
  }
  allocation <- match.arg(allocation)
  domain_sampling <- match.arg(domain_sampling)
  joint <- identical(allocation, "joint")
  if (!is.null(min_n_domain) &&
      (!is.numeric(min_n_domain) || length(min_n_domain) != 1L || is.na(min_n_domain) ||
       min_n_domain <= 0)) {
    stop("'min_n_domain' must be a positive numeric scalar", call. = FALSE)
  }

  check_stage_cost(stage_cost)
  stage_cost <- .reorder_stage_cost(stage_cost)
  stages <- length(stage_cost)
  check_resp_rate(resp_rate_psu, "resp_rate_psu")
  check_resp_rate(resp_rate_ssu, "resp_rate_ssu")
  check_resp_rate(resp_rate)
  if (stages == 2L && resp_rate_ssu != 1) {
    stop("'resp_rate_ssu' is not applicable for 2-stage designs", call. = FALSE)
  }
  if (!is.null(budget)) check_scalar(budget, "budget")
  if (!is.null(n_psu)) check_scalar(n_psu, "n_psu")
  if (!is.null(n_per_psu)) check_scalar(n_per_psu, "n_per_psu")
  if (!is.null(n_per_ssu)) check_scalar(n_per_ssu, "n_per_ssu")
  if (stages == 2L && !is.null(n_per_ssu)) {
    stop("'n_per_ssu' is not applicable for 2-stage designs", call. = FALSE)
  }
  n_fixed <- sum(!is.null(n_psu), !is.null(n_per_psu), !is.null(n_per_ssu))
  if (n_fixed >= stages) {
    stop("cannot fix all stages; use prec_multi_cluster() instead",
         call. = FALSE)
  }
  check_fixed_cost(fixed_cost, budget)

  indicators <- .indicators_var_from_sd(indicators, domains)
  indicators <- .indicators_moe_from_rmoe(indicators, domains)
  info <- .validate_targets(
    indicators,
    TRUE,
    domains = domains,
    stages = stages,
    context = "n_multi_cluster()"
  )
  indicators <- .fill_defaults(
    indicators,
    TRUE,
    stages = stages,
    resp_rate_psu = resp_rate_psu,
    resp_rate_ssu = resp_rate_ssu,
    resp_rate = resp_rate
  )
  indicators <- .convert_moe_to_cv(indicators)

  rv_final <- indicators$unit_relvar[!is.na(indicators$unit_relvar)]
  if (length(rv_final) > 0L &&
      (any(rv_final <= 0) || any(!is.finite(rv_final)))) {
    stop("'unit_relvar' values must be positive and finite", call. = FALSE)
  }
  if (any(indicators$var_ratio_psu <= 0) || any(!is.finite(indicators$var_ratio_psu))) {
    stop("'var_ratio_psu' values must be positive and finite", call. = FALSE)
  }
  if (any(indicators$var_ratio_ssu <= 0) || any(!is.finite(indicators$var_ratio_ssu))) {
    stop("'var_ratio_ssu' values must be positive and finite", call. = FALSE)
  }

  domain_cols <- info$domain_cols
  mode <- if (!is.null(budget)) {
    "budget"
  } else if ("moe" %in% names(indicators) && any(!is.na(indicators$moe))) {
    "moe"
  } else {
    "cv"
  }

  if (length(domain_cols) == 0L) {
    .n_multi_cluster(
      indicators,
      stage_cost,
      budget,
      n_psu,
      n_per_psu,
      n_per_ssu,
      fixed_cost,
      domain_cols = domain_cols,
      mode = mode
    )
  } else {
    .n_multi_domains(
      indicators,
      stage_cost,
      budget,
      n_psu,
      n_per_psu,
      n_per_ssu,
      domain_cols,
      multistage = TRUE,
      joint,
      min_n_domain,
      fixed_cost,
      mode = mode,
      domain_sampling = domain_sampling
    )
  }
}

#' A per-row rate column, defaulting to full response when absent
#'
#' Read rather than filled, so a stage that a design does not have never
#' materializes a column that a round trip would then have to explain.
#' @keywords internal
#' @noRd
.indicator_rate <- function(indicators, name, nr) {
  if (!name %in% names(indicators)) {
    return(rep(1, nr))
  }
  v <- indicators[[name]]
  v[is.na(v)] <- 1
  v
}

#' Refuse the response-rate column that belongs to the other interface
#'
#' The two interfaces spend a response rate at different stages, so each
#' names its own column. An unrecognized column would otherwise be carried
#' along and ignored, which for a response rate means silently planning a
#' design with none.
#' @keywords internal
#' @noRd
.check_resp_rate_column <- function(indicators, multistage, stages = 2L) {
  if (!multistage) {
    wrong <- intersect(c("resp_rate_psu", "resp_rate_ssu"), names(indicators))
    if (length(wrong) == 0L) {
      return(invisible(TRUE))
    }
    stop(
      sprintf(
        "column(s) %s name sampling stages this design does not have; a single-stage design loses ultimate units, so its rate is the plain 'resp_rate'",
        paste(sQuote(wrong), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (stages == 2L && "resp_rate_ssu" %in% names(indicators)) {
    stop(
      "column 'resp_rate_ssu' is not applicable for 2-stage designs: the units inside a PSU are the ultimate ones, so their nonresponse is 'resp_rate'",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Report arguments that moved to the cluster-specific API
#' @keywords internal
#' @noRd
.check_multi_split_args <- function(dots, replacement) {
  moved <- intersect(
    names(dots) %||% character(0),
    c("stage_cost", "budget", "n_psu", "n_per_psu", "n_per_ssu", "allocation",
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

#' Validate indicators data frame
#' @return List with `indicator_type` (per-row "p" or "var") and `domain_cols`.
#' @keywords internal
#' @noRd
.validate_targets <- function(indicators, multistage, domains = NULL,
                              stages = NULL,
                              context = "n_multi_cluster()") {
  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(indicators))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in indicators: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  .check_indicator_columns(indicators, domains)
  .check_resp_rate_column(indicators, multistage, stages %||% 2L)

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  if (!has_p && !has_var) {
    stop("'indicators' must contain 'p' or 'var' column", call. = FALSE)
  }

  has_moe <- "moe" %in% names(indicators)
  has_cv <- "cv" %in% names(indicators)
  if (!has_moe && !has_cv) {
    stop("'indicators' must contain a 'moe', 'rmoe', or 'cv' column",
         call. = FALSE)
  }

  if (has_p) {
    p_vals <- indicators$p[!is.na(indicators$p)]
    if (any(p_vals <= 0 | p_vals >= 1)) {
      stop("all 'p' values must be in (0, 1)", call. = FALSE)
    }
  }

  if (has_var) {
    var_vals <- indicators$var[!is.na(indicators$var)]
    if (any(var_vals <= 0) || any(!is.finite(var_vals))) {
      stop("all 'var' values must be positive and finite", call. = FALSE)
    }
  }

  # Each row needs at least one of p or var (non-NA)
  has_indicator <- rep(FALSE, nrow(indicators))
  if (has_p) {
    has_indicator <- has_indicator | !is.na(indicators$p)
  }
  if (has_var) {
    has_indicator <- has_indicator | !is.na(indicators$var)
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
    both_set <- !is.na(indicators$p) & !is.na(indicators$var)
    if (any(both_set)) {
      stop("each row must have only one of 'p' or 'var'", call. = FALSE)
    }
  }

  if (has_moe && has_cv) {
    both_na <- is.na(indicators$moe) & is.na(indicators$cv)
    both_set <- !is.na(indicators$moe) & !is.na(indicators$cv)
    if (any(both_na)) {
      stop("each row must have either 'moe' or 'cv' specified", call. = FALSE)
    }
    if (any(both_set)) {
      stop("each row must have only one of 'moe' or 'cv'", call. = FALSE)
    }
  }

  if (has_moe) {
    moe_vals <- indicators$moe[!is.na(indicators$moe)]
    if (any(moe_vals <= 0) || any(!is.finite(moe_vals))) {
      stop("'moe' values must be positive and finite", call. = FALSE)
    }
  }

  if ("prop_method" %in% names(indicators)) {
    method_vals <- indicators$prop_method[!is.na(indicators$prop_method)]
    bad_methods <- !method_vals %in% c("wald", "wilson", "logodds", "beta")
    if (any(bad_methods)) {
      stop(
        "'prop_method' values must be one of 'wald', 'wilson', 'logodds', or 'beta'",
        call. = FALSE
      )
    }
  }

  if (multistage) {
    if (!has_cv && !has_moe) {
      stop("multistage mode requires a 'cv', 'moe', or 'rmoe' column in indicators",
           call. = FALSE)
    }
    if (has_moe && any(!is.na(indicators$moe))) {
      moe_rows <- !is.na(indicators$moe)
      var_moe <- has_var & moe_rows & !is.na(indicators$var)
      if (any(var_moe)) {
        if (!"mu" %in% names(indicators) || any(is.na(indicators$mu[var_moe]))) {
          stop(
            "'mu' is required to convert 'moe' to 'cv' for mean indicators ",
            "in multistage mode",
            call. = FALSE
          )
        }
      }
    }
    if (!"icc_psu" %in% names(indicators)) {
      stop(
        "multistage mode requires 'icc_psu' column in indicators",
        call. = FALSE
      )
    }
    if (isTRUE(stages == 3L)) {
      if (!"icc_ssu" %in% names(indicators)) {
        stop(
          "3-stage mode requires a 'icc_ssu' column in indicators",
          call. = FALSE
        )
      }
      if (anyNA(indicators$icc_ssu) || any(!is.finite(indicators$icc_ssu))) {
        stop("'icc_ssu' must contain finite non-missing values",
             call. = FALSE)
      }
    }
    if (has_cv) {
      cv_vals <- indicators$cv[!is.na(indicators$cv)]
      if (length(cv_vals) > 0L && (any(cv_vals <= 0) || any(!is.finite(cv_vals)))) {
        stop("'cv' values must be positive and finite", call. = FALSE)
      }
    }
    d1_vals <- indicators$icc_psu[!is.na(indicators$icc_psu)]
    if (any(d1_vals < 0 | d1_vals > 1)) {
      stop("'icc_psu' values must be in [0, 1]", call. = FALSE)
    }
    .check_cluster_icc_open(d1_vals, context = context)
    if ("icc_ssu" %in% names(indicators)) {
      d2_vals <- indicators$icc_ssu[!is.na(indicators$icc_ssu)]
      if (any(d2_vals < 0 | d2_vals > 1)) {
        stop("'icc_ssu' values must be in [0, 1]", call. = FALSE)
      }
      .check_cluster_icc_open(d2_vals, context = context)
    }
  }

  if (has_var && has_cv) {
    needs_mu <- !is.na(indicators$var) & !is.na(indicators$cv)
    if (any(needs_mu)) {
      if (!"mu" %in% names(indicators) || any(is.na(indicators$mu[needs_mu]))) {
        stop(
          "'mu' is required when 'var' and 'cv' are specified",
          call. = FALSE
        )
      }
    }
  }

  if ("mu" %in% names(indicators)) {
    mu_vals <- indicators$mu[!is.na(indicators$mu)]
    if (any(mu_vals == 0) || any(!is.finite(mu_vals))) {
      stop("'mu' values must be finite and non-zero", call. = FALSE)
    }
  }

  domain_cols <- domains %||% character(0)

  if ("unit_relvar" %in% names(indicators)) {
    rv_vals <- indicators$unit_relvar[!is.na(indicators$unit_relvar)]
    if (any(rv_vals <= 0) || any(!is.finite(rv_vals))) {
      stop("'unit_relvar' values must be positive and finite", call. = FALSE)
    }
  }
  if (multistage) {
    if ("var_ratio_psu" %in% names(indicators)) {
      var_ratio_psu_vals <- indicators$var_ratio_psu[!is.na(indicators$var_ratio_psu)]
      if (any(var_ratio_psu_vals <= 0) || any(!is.finite(var_ratio_psu_vals))) {
        stop("'var_ratio_psu' values must be positive and finite", call. = FALSE)
      }
    }
    if ("var_ratio_ssu" %in% names(indicators)) {
      var_ratio_ssu_vals <- indicators$var_ratio_ssu[!is.na(indicators$var_ratio_ssu)]
      if (any(var_ratio_ssu_vals <= 0) || any(!is.finite(var_ratio_ssu_vals))) {
        stop("'var_ratio_ssu' values must be positive and finite", call. = FALSE)
      }
    }
  }

  .validate_common_columns(indicators)

  list(domain_cols = domain_cols)
}

#' Fill default values in indicators
#' @keywords internal
#' @noRd
.fill_defaults <- function(
  indicators,
  multistage,
  prop_method = "wald",
  stages = if (multistage) 2L else 1L,
  resp_rate_psu = 1,
  resp_rate_ssu = 1,
  resp_rate = 1
) {
  if (!"alpha" %in% names(indicators)) {
    indicators$alpha <- 0.05
  } else {
    indicators$alpha[is.na(indicators$alpha)] <- 0.05
  }

  if (!multistage) {
    if (!"deff" %in% names(indicators)) {
      indicators$deff <- 1
    } else {
      indicators$deff[is.na(indicators$deff)] <- 1
    }

    if (!"N" %in% names(indicators)) {
      indicators$N <- Inf
    } else {
      indicators$N[is.na(indicators$N)] <- Inf
    }
  }

  if (multistage) {
    if (!"var_ratio_psu" %in% names(indicators)) {
      indicators$var_ratio_psu <- 1
    } else {
      indicators$var_ratio_psu[is.na(indicators$var_ratio_psu)] <- 1
    }
    # var_ratio_ssu is the within-PSU share of unit variance, not a free parameter:
    # var_ratio_psu * icc_psu + var_ratio_ssu must equal 1. See .var_ratio_ssu_default().
    if ("icc_ssu" %in% names(indicators) && "icc_psu" %in% names(indicators)) {
      implied <- .var_ratio_ssu_default(indicators$var_ratio_psu, indicators$icc_psu)
      if (!"var_ratio_ssu" %in% names(indicators)) {
        indicators$var_ratio_ssu <- implied
      } else {
        indicators$var_ratio_ssu[is.na(indicators$var_ratio_ssu)] <- implied[is.na(indicators$var_ratio_ssu)]
      }
    } else if (!"var_ratio_ssu" %in% names(indicators)) {
      indicators$var_ratio_ssu <- 1
    } else {
      indicators$var_ratio_ssu[is.na(indicators$var_ratio_ssu)] <- 1
    }
  }

  # A simple design has only ultimate-unit nonresponse. A multistage design
  # keeps a separate rate for every stage it actually has.
  rate_cols <- if (multistage) {
    c(
      resp_rate_psu = resp_rate_psu,
      if (stages == 3L) c(resp_rate_ssu = resp_rate_ssu),
      resp_rate = resp_rate
    )
  } else {
    c(resp_rate = resp_rate)
  }
  for (rate_col in names(rate_cols)) {
    if (!rate_col %in% names(indicators)) {
      indicators[[rate_col]] <- unname(rate_cols[[rate_col]])
    } else {
      indicators[[rate_col]][is.na(indicators[[rate_col]])] <-
        unname(rate_cols[[rate_col]])
    }
  }

  if (!"prop_method" %in% names(indicators)) {
    indicators$prop_method <- prop_method
  } else {
    indicators$prop_method[is.na(indicators$prop_method)] <- prop_method
  }

  if (!"df" %in% names(indicators)) {
    indicators$df <- NA_real_
  }

  if (!"unit_relvar" %in% names(indicators)) {
    indicators$unit_relvar <- NA_real_
  }
  indicators$unit_relvar <- .derive_unit_relvar(indicators, require_all = multistage)

  indicators
}


#' Derive unit relvariance from p or var/mu
#' @keywords internal
#' @noRd
.derive_unit_relvar <- function(indicators, require_all = FALSE) {
  rv <- indicators$unit_relvar
  needs <- is.na(rv)

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  has_mu <- "mu" %in% names(indicators)

  for (i in which(needs)) {
    if (has_p && !is.na(indicators$p[i])) {
      rv[i] <- (1 - indicators$p[i]) / indicators$p[i]
    } else if (has_var && !is.na(indicators$var[i])) {
      if (has_mu && !is.na(indicators$mu[i])) {
        rv[i] <- indicators$var[i] / indicators$mu[i]^2
      } else if (require_all) {
        stop(
          sprintf("row %d: 'mu' is required to derive 'unit_relvar' from 'var'", i),
          call. = FALSE
        )
      }
    }
  }

  if (require_all && anyNA(rv)) {
    stop(
      sprintf(
        "row(s) %s: could not derive 'unit_relvar' - provide 'p', or 'var'+'mu', or 'unit_relvar' directly",
        paste(which(is.na(rv)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  rv
}

#' Convert moe to cv in multistage indicators
#'
#' The multistage model is driven by a relative standard error, so a margin
#' of error target has to be restated as one. For a mean row that is exact:
#' `cv = moe / (z |mu|)`, on the magnitude because a mean may be negative
#' while a CV may not.
#'
#' For a proportion row it is exact only under `"wald"`, the one method whose
#' half-width is `z se`. The others close their interval at a different
#' standard error, so the restatement goes through the effective sample size
#' the row's own method needs to reach that margin and reports the sampling
#' CV there. Under `"wald"` this reproduces `moe / (z p)`, so no Wald result
#' moves; under the other three it is what makes `prop_method` reach the
#' multistage path at all.
#'
#' Rows that already have cv are left unchanged.
#' @keywords internal
#' @noRd
.convert_moe_to_cv <- function(indicators) {
  has_moe <- "moe" %in% names(indicators)
  if (!has_moe) return(indicators)

  moe_rows <- !is.na(indicators$moe)
  if (!any(moe_rows)) return(indicators)

  if (!"cv" %in% names(indicators)) {
    indicators$cv <- NA_real_
  }

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  has_method <- "prop_method" %in% names(indicators)

  for (i in which(moe_rows)) {
    df_i <- .row_df(indicators, i)
    if (has_p && !is.na(indicators$p[i])) {
      p_i <- indicators$p[i]
      method_i <- if (has_method) indicators$prop_method[i] else "wald"
      n_eff <- .n_prop_effective(p_i, indicators$moe[i], indicators$alpha[i],
                                 method_i, df_i)
      indicators$cv[i] <- sqrt(p_i * (1 - p_i) / n_eff) / p_i
    } else if (has_var && !is.na(indicators$var[i])) {
      z <- .q_alpha(indicators$alpha[i], df_i)
      indicators$cv[i] <- indicators$moe[i] / (z * abs(indicators$mu[i]))
    }
  }

  indicators
}

#' Simple mode: compute n per indicator, take max
#' @keywords internal
#' @noRd
.n_multi_simple <- function(indicators, domain_cols = character(0), mode = "moe",
                           prop_method = "wald") {
  simple <- .compute_simple_n(indicators)
  n_vec <- simple$n
  cv_target_vec <- simple$cv_target

  # A case floor raises the row's own size, so the binding row is chosen
  # among sizes that already satisfy both constraints. `.cv_target` keeps
  # its precision meaning; `.cv_achieved` is read at the size that wins.
  floor_n <- .multi_min_cases_n(indicators)
  if (!is.null(floor_n)) {
    n_vec <- pmax(n_vec, floor_n, na.rm = TRUE)
  }

  idx <- which.max(n_vec)
  n_max <- n_vec[idx]

  labels <- if ("name" %in% names(indicators)) {
    indicators$name
  } else {
    seq_len(nrow(indicators))
  }
  binding_name <- labels[idx]

  has_p <- "p" %in% names(indicators)
  has_mu <- "mu" %in% names(indicators)
  cv_achieved_vec <- vapply(seq_len(nrow(indicators)), function(i) {
    prec <- suppressWarnings(
      if (has_p && !is.na(indicators$p[i])) {
        .prec_engine_prop(indicators$p[i], n_max, indicators$alpha[i],
                          indicators$N[i], indicators$deff[i],
                          indicators$resp_rate[i], indicators$prop_method[i],
                          .row_df(indicators, i))
      } else {
        .prec_engine_mean(indicators$var[i],
                          if (has_mu) indicators$mu[i] else NULL,
                          n_max, indicators$alpha[i], indicators$N[i],
                          indicators$deff[i], indicators$resp_rate[i],
                          .row_df(indicators, i))
      }
    )
    prec$cv
  }, numeric(1L))

  detail <- data.frame(
    name = labels,
    .n = n_vec,
    .cv_target = cv_target_vec,
    .cv_achieved = cv_achieved_vec,
    .binding = seq_len(nrow(indicators)) == idx
  )

  .new_svyplan_n(
    n = n_max,
    type = "multi",
    params = list(domain_cols = domain_cols, mode = mode,
                  prop_method = prop_method),
    indicators = indicators,
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
.compute_simple_n <- function(indicators) {
  nr <- nrow(indicators)
  n_vec <- numeric(nr)
  cv_vec <- numeric(nr)

  has_p <- "p" %in% names(indicators)
  has_moe <- "moe" %in% names(indicators)
  has_mu <- "mu" %in% names(indicators)

  for (i in seq_len(nr)) {
    is_prop <- has_p && !is.na(indicators$p[i])
    use_moe <- has_moe && !is.na(indicators$moe[i])

    if (is_prop) {
      res_i <- n_prop.default(
        p = indicators$p[i],
        moe = if (use_moe) indicators$moe[i] else NULL,
        cv = if (use_moe) NULL else indicators$cv[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        method = indicators$prop_method[i],
        df = .row_df(indicators, i)
      )
    } else {
      res_i <- n_mean.default(
        var = indicators$var[i],
        mu = if (has_mu && !is.na(indicators$mu[i])) indicators$mu[i] else NULL,
        moe = if (use_moe) indicators$moe[i] else NULL,
        cv = if (use_moe) NULL else indicators$cv[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        df = .row_df(indicators, i)
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
.n_multi_cluster <- function(indicators, stage_cost, budget, n_psu,
                            n_per_psu = NULL, n_per_ssu = NULL, fixed_cost = 0,
                            domain_cols = character(0), mode = "cv",
                            prop_method = "wald") {
  .stop_min_cases_column(
    indicators,
    "a multistage design sizes stages against a cost, with no single total for a count to raise. Size the count-driven indicator with n_prop(min_cases = ) and set the stage takes around it"
  )
  stages <- length(stage_cost)
  if (stages == 2L) {
    .n_multi_2stage(indicators, stage_cost, budget, n_psu, n_per_psu, fixed_cost,
                    domain_cols = domain_cols, mode = mode,
                    prop_method = prop_method)
  } else {
    .n_multi_3stage(indicators, stage_cost, budget, n_psu, n_per_psu, n_per_ssu,
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
  center <- max(1L, as.integer(round(continuous)))
  local <- seq.int(max(1L, center - 20L), min(upper_i, center + 20L))
  unique(c(base, local))
}

#' Constraint-preserving whole-unit design for two-stage n_multi
#' @keywords internal
#' @noRd
.op_multi_2stage <- function(n1_required, cv_fn, cv_t, stage_cost, budget,
                             n_psu, n_per_psu, fixed_cost, cont_m) {
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
    n_per_psu, cont_m, upper = upper_m, base_limit = 100000L
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
    n = c(n_psu = a_best, n_per_psu = m_best),
    total_n = a_best * m_best,
    cost = fixed_cost + a_best * (C1 + C2 * m_best),
    cv = max(cvs),
    cv_by_target = cvs
  )
}

#' Constraint-preserving whole-unit design for three-stage n_multi
#' @keywords internal
#' @noRd
.op_multi_3stage <- function(coef, cv_fn, cv_t, stage_cost, budget,
                             n_psu, n_per_psu, n_per_ssu, fixed_cost,
                             cont_m) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  variable_budget <- if (is.null(budget)) NULL else budget - fixed_cost
  a_fixed <- if (is.null(n_psu)) NULL else max(1L, as.integer(round(n_psu)))
  a_min <- a_fixed %||% 1L

  upper_m <- if (is.null(variable_budget)) {
    Inf
  } else {
    (variable_budget / a_min - C1) / (C2 + C3)
  }
  m_cand <- .multi_stage_candidates(n_per_psu, cont_m, upper_m)

  hit <- .op_cluster3_search(
    coef$alpha, coef$beta, coef$gamma, stage_cost, variable_budget,
    n_psu, n_per_ssu, m_cand,
    msg_budget = "no whole-unit n_multi design fits the budget",
    msg_cv = "target CV is not achievable with the fixed whole PSU count"
  )
  a_best <- hit[["n_psu"]]
  m_best <- hit[["n_per_psu"]]
  q_best <- hit[["n_per_ssu"]]
  cvs <- cv_fn(a_best, m_best, q_best)
  list(
    n = c(n_psu = a_best, n_per_psu = m_best, n_per_ssu = q_best),
    total_n = a_best * m_best * q_best,
    cost = fixed_cost + a_best * (C1 + C2 * m_best + C3 * m_best * q_best),
    cv = max(cvs),
    cv_by_target = cvs
  )
}

#' 2-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_2stage <- function(indicators, stage_cost, budget, n_psu,
                           n_per_psu = NULL, fixed_cost = 0,
                           domain_cols = character(0), mode = "cv",
                           prop_method = "wald") {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  nr <- nrow(indicators)

  cv_t <- indicators$cv
  icc <- indicators$icc_psu
  unit_relvar <- indicators$unit_relvar
  var_ratio <- indicators$var_ratio_psu
  rr <- indicators$resp_rate_psu
  ru <- indicators$resp_rate
  labels <- if ("name" %in% names(indicators)) indicators$name else seq_len(nr)

  # n_per_psu is the gross take; each row reads the take it realizes.
  n1_required <- function(n_per_psu) {
    vapply(
      seq_len(nr),
      function(j) {
        mr <- n_per_psu * ru[j]
        unit_relvar[j] *
          var_ratio[j] *
          (1 + icc[j] * (mr - 1)) /
          (mr * cv_t[j]^2 * rr[j])
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, n_per_psu) {
    vapply(
      seq_len(nr),
      function(j) {
        mr <- n_per_psu * ru[j]
        sqrt(
          unit_relvar[j] *
            var_ratio[j] /
            (n1 * rr[j] * mr) *
            (1 + icc[j] * (mr - 1))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(budget)) {
    if (!is.null(n_per_psu)) {
      n_per_psu_opt <- n_per_psu
      n1_vals <- n1_required(n_per_psu_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost + n1_opt * (C1 + C2 * n_per_psu_opt)
      binding_idx <- which.max(n1_vals)
    } else if (!is.null(n_psu)) {
      # A fixed PSU count leaves the take as the only free stage, so it is
      # solved rather than cost-optimized. Feasible only above the
      # between-PSU floor A * icc.
      A <- unit_relvar * var_ratio / (cv_t^2 * rr)
      floor_psu <- A * icc
      if (any(n_psu <= floor_psu + 1e-12)) {
        stop(
          sprintf(
            "a fixed 'n_psu' of %g cannot reach the target for %s: the between-PSU term alone needs more than %g clusters however large the take",
            n_psu,
            paste(sQuote(labels[n_psu <= floor_psu + 1e-12]), collapse = ", "),
            max(floor_psu[n_psu <= floor_psu + 1e-12])
          ),
          call. = FALSE
        )
      }
      takes <- A * (1 - icc) / (ru * (n_psu - floor_psu))
      binding_idx <- which.max(takes)
      n_per_psu_opt <- max(takes)
      n1_opt <- n_psu
      total_cost <- fixed_cost + n1_opt * (C1 + C2 * n_per_psu_opt)
    } else {
      cost_fn <- function(n_per_psu) {
        n1 <- max(n1_required(n_per_psu))
        n1 * (C1 + C2 * n_per_psu)
      }

      # The take grows as r falls, so the bracket is drawn at r and doubled to
      # keep the optimum interior across indicators.
      upper <- max(10, 2 * max(sqrt(C1 / C2 * (1 - icc) / (icc * ru))))
      opt <- optimize(cost_fn, interval = c(1, upper),
                      tol = .Machine$double.eps^0.5)
      n_per_psu_opt <- opt$minimum

      n1_vals <- n1_required(n_per_psu_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost + n1_opt * (C1 + C2 * n_per_psu_opt)
      binding_idx <- which.max(n1_vals)
    }

    cv_achieved <- cv_achieved_fn(n1_opt, n_per_psu_opt)
  } else {
    var_budget <- budget - fixed_cost
    bres <- .eval_2stage_budget(
      cv_t,
      icc,
      unit_relvar,
      var_ratio,
      rr,
      ru,
      C1,
      C2,
      var_budget,
      n_psu,
      n_per_psu
    )
    n1_opt <- bres$n1
    n_per_psu_opt <- bres$n_per_psu
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- budget
  }

  n_vec <- c(n_psu = n1_opt, n_per_psu = n_per_psu_opt)
  total_n <- prod(n_vec)
  operational <- .op_multi_2stage(
    n1_required, cv_achieved_fn, cv_t, stage_cost, budget,
    n_psu, n_per_psu, fixed_cost, cont_m = n_per_psu_opt
  )

  n1_per <- n1_required(n_per_psu_opt)
  n_per <- n1_per * n_per_psu_opt

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
  if (!is.null(n_per_psu)) {
    params$n_per_psu <- n_per_psu
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
    indicators = indicators,
    detail = detail,
    binding = labels[binding_idx],
    operational = operational
  )
}

#' 3-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_3stage <- function(indicators, stage_cost, budget, n_psu,
                           n_per_psu = NULL, n_per_ssu = NULL, fixed_cost = 0,
                           domain_cols = character(0), mode = "cv",
                           prop_method = "wald") {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  nr <- nrow(indicators)

  cv_t <- indicators$cv
  icc_psu <- indicators$icc_psu
  icc_ssu <- if ("icc_ssu" %in% names(indicators)) {
    indicators$icc_ssu
  } else {
    rep(0, nr)
  }
  unit_relvar <- indicators$unit_relvar
  var_ratio_psu <- indicators$var_ratio_psu
  var_ratio_ssu <- indicators$var_ratio_ssu
  rr <- indicators$resp_rate_psu
  rs <- .indicator_rate(indicators, "resp_rate_ssu", nr)
  ru <- indicators$resp_rate
  labels <- if ("name" %in% names(indicators)) {
    indicators$name
  } else {
    seq_len(nr)
  }

  # Both stage takes are gross; each row reads the takes it realizes.
  n1_required <- function(n_per_psu, n_per_ssu) {
    vapply(
      seq_len(nr),
      function(j) {
        mr <- n_per_psu * rs[j]
        qr <- n_per_ssu * ru[j]
        unit_relvar[j] /
          (cv_t[j]^2 * mr * qr * rr[j]) *
          (var_ratio_psu[j] *
            icc_psu[j] *
            mr *
            qr +
            var_ratio_ssu[j] * (1 + icc_ssu[j] * (qr - 1)))
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, n_per_psu, n_per_ssu) {
    vapply(
      seq_len(nr),
      function(j) {
        mr <- n_per_psu * rs[j]
        qr <- n_per_ssu * ru[j]
        sqrt(
          unit_relvar[j] /
            (n1 * rr[j] * mr * qr) *
            (var_ratio_psu[j] *
              icc_psu[j] *
              mr *
              qr +
              var_ratio_ssu[j] * (1 + icc_ssu[j] * (qr - 1)))
        )
      },
      numeric(1L)
    )
  }

  # n1_required() in the gross takes, so the whole-unit search can invert it.
  search_coef <- list(
    alpha = unit_relvar * var_ratio_psu * icc_psu / (rr * cv_t^2),
    beta = unit_relvar * var_ratio_ssu * icc_ssu / (rr * rs * cv_t^2),
    gamma = unit_relvar * var_ratio_ssu * (1 - icc_ssu) / (rr * rs * ru * cv_t^2)
  )

  # One representation for every fixed-stage branch, rather than re-deriving
  # the algebra per branch and dropping a stage rate.
  ps_required <- function(ss, n1) {
    room <- n1 - search_coef$alpha
    per <- (search_coef$beta + search_coef$gamma / ss) / room
    per[room <= 0] <- Inf
    max(per)
  }

  ss_required <- function(ps, n1) {
    room <- n1 - search_coef$alpha - search_coef$beta / ps
    per <- search_coef$gamma / (ps * room)
    per[room <= 0] <- Inf
    max(per)
  }

  solve_for <- if (!is.null(n_psu) && !is.null(n_per_psu)) {
    "n3"
  } else if (!is.null(n_psu)) {
    "n2"
  } else {
    "n1"
  }

  n_free <- 3L - sum(!is.null(n_psu), !is.null(n_per_psu), !is.null(n_per_ssu))

  if (is.null(budget)) {
    if (n_free == 3L) {
      cost_fn <- function(par) {
        ps <- par[1L]
        ss <- par[2L]
        n1 <- max(n1_required(ps, ss))
        n1 * (C1 + C2 * ps + C3 * ps * ss)
      }

      init_n_per_ssu <- max(2, sqrt(C2 / C3))
      init_n_per_psu <- max(2, sqrt(C1 / C2))
      upper_ps <- max(1000, 10 * init_n_per_psu)
      upper_ss <- max(1000, 10 * init_n_per_ssu)
      opt <- optim(
        par = c(init_n_per_psu, init_n_per_ssu),
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
            "optimal stage size reached the search upper bound (n_per_psu <= %.0f, n_per_ssu <= %.0f); result may be unreliable -- review stage costs and target CVs",
            upper_ps,
            upper_ss
          ),
          call. = FALSE
        )
      }

      n_per_psu_opt <- opt$par[1L]
      n_per_ssu_opt <- opt$par[2L]
      n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost +
        n1_opt * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
      binding_idx <- which.max(n1_vals)
    } else if (n_free == 2L) {
      if (solve_for == "n2" && is.null(n_per_ssu)) {
        n1_opt <- n_psu

        cv_floor <- .multistage_cv_floor(
          unit_relvar, var_ratio_psu, icc_psu,
          rr = rr, n_psu = n_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          unit_relvar, var_ratio_psu, icc_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )

        n_per_psu_required_fn <- function(ss) ps_required(ss, n_psu)

        cost_fn_fixed <- function(ss) {
          ps <- n_per_psu_required_fn(ss)
          n_psu * (C1 + C2 * ps + C3 * ps * ss)
        }

        n_per_ssu_analytic <- vapply(
          seq_len(nr),
          function(j) {
            if (icc_ssu[j] <= 0) 1
            else sqrt((1 - icc_ssu[j]) / icc_ssu[j] * C2 / C3)
          },
          numeric(1L)
        )
        upper_n_per_ssu <- max(10, 3 * max(n_per_ssu_analytic))

        opt <- optimize(cost_fn_fixed, interval = c(1, upper_n_per_ssu))
        n_per_ssu_opt <- opt$minimum
        n_per_psu_opt <- n_per_psu_required_fn(n_per_ssu_opt)

        if (!is.finite(n_per_psu_opt) || n_per_psu_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }

        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
      } else if (!is.null(n_per_ssu) && is.null(n_psu) && is.null(n_per_psu)) {
        n_per_ssu_opt <- n_per_ssu
        cost_fn_ss <- function(ps) {
          n1 <- max(n1_required(ps, n_per_ssu_opt))
          n1 * (C1 + C2 * ps + C3 * ps * n_per_ssu_opt)
        }
        n_per_psu_analytic <- vapply(
          seq_len(nr),
          function(j) {
            sqrt(
              C1 * var_ratio_ssu[j] * (1 + icc_ssu[j] * (n_per_ssu_opt - 1)) /
                (var_ratio_psu[j] * icc_psu[j] * n_per_ssu_opt *
                   (C2 + C3 * n_per_ssu_opt))
            )
          },
          numeric(1L)
        )
        upper <- max(10, n_per_psu_analytic)
        opt <- optimize(cost_fn_ss, interval = c(1, upper))
        n_per_psu_opt <- opt$minimum
        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        n1_opt <- max(n1_vals)
        total_cost <- fixed_cost +
          n1_opt * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
      } else if (!is.null(n_per_psu) && is.null(n_psu) && is.null(n_per_ssu)) {
        n_per_psu_opt <- n_per_psu
        cost_fn_ps <- function(ss) {
          n1 <- max(n1_required(n_per_psu_opt, ss))
          n1 * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * ss)
        }
        n_per_ssu_analytic <- vapply(
          seq_len(nr),
          function(j) {
            if (icc_ssu[j] <= 0) 1
            else sqrt((1 - icc_ssu[j]) / icc_ssu[j] * C2 / C3)
          },
          numeric(1L)
        )
        upper <- max(10, 3 * max(n_per_ssu_analytic))
        opt <- optimize(cost_fn_ps, interval = c(1, upper))
        n_per_ssu_opt <- opt$minimum
        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        n1_opt <- max(n1_vals)
        total_cost <- fixed_cost +
          n1_opt * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
      } else {
        n1_opt <- n_psu
        n_per_ssu_opt <- n_per_ssu
        cv_floor <- .multistage_cv_floor(
          unit_relvar, var_ratio_psu, icc_psu,
          rr = rr, n_psu = n_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          unit_relvar, var_ratio_psu, icc_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )
        n_per_psu_required_fn2 <- function(j) {
          denom <- cv_t[j]^2 * n_psu * rr[j] / (unit_relvar[j] * var_ratio_ssu[j]) -
            var_ratio_psu[j] * icc_psu[j] / var_ratio_ssu[j]
          if (denom <= 0) Inf
          else (1 + icc_ssu[j] * (n_per_ssu - 1)) / (n_per_ssu * denom)
        }
        psu_per <- vapply(seq_len(nr), n_per_psu_required_fn2, numeric(1L))
        n_per_psu_opt <- max(psu_per)
        if (!is.finite(n_per_psu_opt) || n_per_psu_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }
        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
      }
    } else {
      if (solve_for == "n1") {
        n_per_psu_opt <- n_per_psu
        n_per_ssu_opt <- n_per_ssu
        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        n1_opt <- max(n1_vals)
        total_cost <- fixed_cost +
          n1_opt * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
      } else if (solve_for == "n2") {
        n1_opt <- n_psu
        n_per_ssu_opt <- n_per_ssu
        cv_floor <- .multistage_cv_floor(
          unit_relvar, var_ratio_psu, icc_psu,
          rr = rr, n_psu = n_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          unit_relvar, var_ratio_psu, icc_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )
        n_per_psu_opt <- ps_required(n_per_ssu, n_psu)
        if (!is.finite(n_per_psu_opt) || n_per_psu_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }
        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
      } else {
        n1_opt <- n_psu
        n_per_psu_opt <- n_per_psu
        cv_floor <- .multistage_cv_floor(
          unit_relvar, var_ratio_psu, icc_psu,
          var_ratio_ssu = var_ratio_ssu, icc_ssu = icc_ssu,
          rr = rr, n_psu = n_psu, n_per_psu = n_per_psu
        )
        .check_multistage_feasibility(
          cv_t, cv_floor, n_psu,
          unit_relvar, var_ratio_psu, icc_psu,
          var_ratio_ssu = var_ratio_ssu, icc_ssu = icc_ssu,
          n_per_psu = n_per_psu,
          rr = rr, labels = labels, context = "n_multi_cluster()"
        )
        n_per_ssu_opt <- ss_required(n_per_psu, n_psu)
        if (!is.finite(n_per_ssu_opt) || n_per_ssu_opt <= 0) {
          stop(
            "target CV is too small for the given fixed stage sizes and parameters",
            call. = FALSE
          )
        }
        n1_vals <- n1_required(n_per_psu_opt, n_per_ssu_opt)
        binding_idx <- which.max(n1_vals)
        total_cost <- fixed_cost +
          n_psu * (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
      }
    }

    cv_achieved <- cv_achieved_fn(n1_opt, n_per_psu_opt, n_per_ssu_opt)
  } else {
    var_budget <- budget - fixed_cost
    bres <- .eval_3stage_budget(
      cv_t,
      icc_psu,
      icc_ssu,
      unit_relvar,
      var_ratio_psu,
      var_ratio_ssu,
      rr,
      rs,
      ru,
      C1,
      C2,
      C3,
      var_budget,
      n_psu,
      n_per_psu,
      n_per_ssu
    )
    n1_opt <- bres$n1
    n_per_psu_opt <- bres$n_per_psu
    n_per_ssu_opt <- bres$n_per_ssu
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- fixed_cost + bres$cost
  }

  n_vec <- c(n_psu = n1_opt, n_per_psu = n_per_psu_opt, n_per_ssu = n_per_ssu_opt)
  total_n <- prod(n_vec)
  operational <- .op_multi_3stage(
    search_coef, cv_achieved_fn, cv_t, stage_cost, budget,
    n_psu, n_per_psu, n_per_ssu, fixed_cost,
    cont_m = n_per_psu_opt
  )

  n1_per <- n1_required(n_per_psu_opt, n_per_ssu_opt)
  n_per <- n1_per * n_per_psu_opt * n_per_ssu_opt

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
  if (!is.null(n_per_psu)) params$n_per_psu <- n_per_psu
  if (!is.null(n_per_ssu)) params$n_per_ssu <- n_per_ssu
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
    indicators = indicators,
    detail = detail,
    binding = labels[binding_idx],
    operational = operational
  )
}

#' Solve n_multi independently per domain, then aggregate
#' @keywords internal
#' @noRd
.n_multi_domains <- function(
  indicators,
  stage_cost,
  budget,
  n_psu,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  domain_cols,
  multistage,
  joint = FALSE,
  min_n_domain = NULL,
  fixed_cost = 0,
  mode = "moe",
  prop_method = "wald",
  domain_sampling = "separate"
) {
  key <- .domain_key(indicators, domain_cols)
  domain_levels <- unique(key)
  domain_keys <- factor(key, levels = domain_levels)
  split_idx <- split(seq_len(nrow(indicators)), domain_keys)

  if (multistage && joint && !is.null(budget) && length(domain_levels) > 1L) {
    return(.n_multi_domains_joint(
      indicators,
      stage_cost,
      budget,
      n_psu,
      n_per_psu,
      n_per_ssu,
      domain_cols,
      domain_levels,
      split_idx,
      min_n_domain,
      fixed_cost,
      mode = mode,
      prop_method = prop_method
    ))
  }

  results <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- indicators[rows, , drop = FALSE]
    # Drop domain columns for inner solve
    sub_inner <- sub[, !names(sub) %in% domain_cols, drop = FALSE]
    if (multistage) {
      .n_multi_cluster(sub_inner, stage_cost, budget, n_psu, n_per_psu,
                       n_per_ssu, fixed_cost)
    } else {
      .n_multi_simple(sub_inner)
    }
  })
  names(results) <- domain_levels

  if (multistage) {
    .aggregate_cluster_domains(
      results,
      indicators,
      domain_cols,
      domain_levels,
      split_idx,
      stage_cost,
      budget,
      n_psu,
      n_per_psu,
      n_per_ssu,
      min_n_domain,
      fixed_cost,
      mode = mode,
      prop_method = prop_method,
      domain_sampling = domain_sampling
    )
  } else {
    .aggregate_simple_domains(
      results,
      indicators,
      domain_cols,
      domain_levels,
      split_idx,
      min_n_domain,
      mode = mode,
      prop_method = prop_method,
      domain_sampling = domain_sampling
    )
  }
}

#' Domain shares for a natural-incidence total
#'
#' Under natural incidence each domain turns up in one population sample at
#' its own rate, so the size that yields every domain quota is
#' `max(.n / share)`: the binding domain is the one whose requirement is
#' largest relative to how often it appears. The shares have to describe
#' parts of one population, so they may sum to less than 1 when the domains
#' cover only part of it, but not to more.
#' @keywords internal
#' @noRd
.domain_shares <- function(indicators, domain_levels, split_idx) {
  if (!"share" %in% names(indicators)) {
    stop(
      "domain_sampling = \"natural\" needs a 'share' column giving each domain's expected share of the population. Without it there is no one overall size to report, only the per-domain requirements in $domains",
      call. = FALSE
    )
  }
  shares <- vapply(domain_levels, function(lev) {
    v <- indicators$share[split_idx[[lev]]]
    if (anyNA(v) || any(!is.finite(v)) || any(v <= 0) || any(v > 1)) {
      stop("'share' must be in (0, 1] for every indicator row", call. = FALSE)
    }
    if (diff(range(v)) > 1e-8) {
      stop(
        sprintf("'share' identifies a domain, so it must be constant within one; it varies within %s", lev),
        call. = FALSE
      )
    }
    v[1L]
  }, numeric(1L))

  total <- sum(shares)
  if (total > 1 + 1e-6) {
    stop(
      sprintf("domain 'share' values sum to %.4f; they are shares of one population and cannot exceed 1", total),
      call. = FALSE
    )
  }
  as.numeric(shares)
}

#' Aggregate simple-mode domain results
#' @keywords internal
#' @noRd
.aggregate_simple_domains <- function(
  results,
  indicators,
  domain_cols,
  domain_levels,
  split_idx,
  min_n_domain = NULL,
  mode = "moe",
  prop_method = "wald",
  domain_sampling = "separate"
) {
  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- indicators[rows[1L], domain_cols, drop = FALSE]
    dom_vals$.n <- res$n
    dom_vals$.binding <- res$binding
    dom_vals
  })
  domains <- do.call(rbind, domain_rows)
  rownames(domains) <- NULL

  if (!is.null(min_n_domain)) {
    floored <- domains$.n < min_n_domain
    domains$.n[floored] <- min_n_domain
    domains$.binding[floored] <- "(min_n_domain)"
  }

  n_domain_max <- max(domains$.n)
  binding_label <- domains$.binding[which.max(domains$.n)]

  if (identical(domain_sampling, "natural")) {
    domains$.share <- .domain_shares(indicators, domain_levels, split_idx)
    yields <- domains$.n / domains$.share
    n_overall <- max(yields)
    binding_label <- domains$.binding[which.max(yields)]
  } else {
    n_overall <- sum(domains$.n)
  }

  res <- .new_svyplan_n(
    n = n_overall,
    type = "multi",
    params = list(domain_cols = domain_cols, mode = mode,
                  prop_method = prop_method,
                  domain_sampling = domain_sampling),
    indicators = indicators,
    detail = NULL,
    binding = binding_label,
    domains = domains
  )
  res$n_domain_max <- n_domain_max
  if (!is.null(min_n_domain)) {
    res$params$min_n_domain <- min_n_domain
  }
  res
}

#' Aggregate multistage domain results
#' @keywords internal
#' @noRd
.aggregate_cluster_domains <- function(
  results,
  indicators,
  domain_cols,
  domain_levels,
  split_idx,
  stage_cost,
  budget = NULL,
  n_psu = NULL,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  min_n_domain = NULL,
  fixed_cost = 0,
  mode = "cv",
  prop_method = "wald",
  domain_sampling = "separate"
) {
  if (identical(domain_sampling, "natural")) {
    stop(
      "domain_sampling = \"natural\" is not available for multistage designs: the expected yield of a domain then depends on how its members sit inside PSUs and SSUs, not on its share of the population alone, and that model is not in the package. Size the domains separately and read $domains",
      call. = FALSE
    )
  }
  stages <- results[[1L]]$stages
  stage_names <- names(results[[1L]]$n)

  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- indicators[rows[1L], domain_cols, drop = FALSE]
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

  if (!is.null(min_n_domain)) {
    below <- which(domains$.total_n < min_n_domain)
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
          "domain(s) %s have total_n below min_n_domain = %g",
          paste(sQuote(dom_labels), collapse = ", "),
          min_n_domain
        ),
        call. = FALSE
      )
    }
  }

  total_n <- sum(ceiling(domains$.total_n))
  total_cost <- sum(domains$.cost)
  worst_cv_idx <- which.max(domains$.cv)

  # No single stage vector: $domains holds the fieldable per-domain sizes.
  n_vec <- rep(NA_real_, stages)
  names(n_vec) <- stage_names

  params <- list(stage_cost = stage_cost, domain_cols = domain_cols,
                  mode = mode, prop_method = prop_method,
                  domain_sampling = domain_sampling)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(n_per_psu)) params$n_per_psu <- n_per_psu
  if (!is.null(n_per_ssu)) params$n_per_ssu <- n_per_ssu
  if (!is.null(min_n_domain)) {
    params$min_n_domain <- min_n_domain
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
    indicators = indicators,
    detail = NULL,
    binding = domains$.binding[worst_cv_idx],
    domains = domains
  )
}

#' @keywords internal
#' @noRd
.n_multi_domains_joint <- function(
  indicators,
  stage_cost,
  budget,
  n_psu,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  domain_cols,
  domain_levels,
  split_idx,
  min_n_domain = NULL,
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
    sub <- indicators[rows, , drop = FALSE]
    labels <- if ("name" %in% names(sub)) {
      sub$name
    } else {
      seq_along(rows)
    }
    list(
      cv_t = sub$cv,
      icc_psu = sub$icc_psu,
      icc_ssu = if ("icc_ssu" %in% names(sub)) {
        sub$icc_ssu
      } else {
        rep(0, length(rows))
      },
      unit_relvar = sub$unit_relvar,
      var_ratio_psu = sub$var_ratio_psu,
      var_ratio_ssu = sub$var_ratio_ssu,
      resp_rate_psu = sub$resp_rate_psu,
      resp_rate_ssu = if ("resp_rate_ssu" %in% names(sub)) {
        sub$resp_rate_ssu
      } else {
        rep(1, length(rows))
      },
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
        p$icc_psu,
        p$unit_relvar,
        p$var_ratio_psu,
        p$resp_rate_psu,
        p$resp_rate,
        C1,
        C2,
        budget_d,
        n_psu,
        n_per_psu
      )
    } else {
      .eval_3stage_budget(
        p$cv_t,
        p$icc_psu,
        p$icc_ssu,
        p$unit_relvar,
        p$var_ratio_psu,
        p$var_ratio_ssu,
        p$resp_rate_psu,
        p$resp_rate_ssu,
        p$resp_rate,
        C1,
        C2,
        C3,
        budget_d,
        n_psu,
        n_per_psu,
        n_per_ssu
      )
    }
  }

  total_n_for <- function(bres) {
    if (stages == 2L) {
      bres$n1 * bres$n_per_psu
    } else {
      bres$n1 * bres$n_per_psu * bres$n_per_ssu
    }
  }

  if (!is.null(min_n_domain)) {
    full_total <- vapply(
      seq_len(nd),
      function(d) {
        tryCatch(total_n_for(eval_domain(d, var_budget)), error = function(e) 0)
      },
      numeric(1L)
    )
    for (d in seq_len(nd)) {
      if (full_total[d] < min_n_domain) {
        lab <- paste(
          unlist(lapply(indicators[split_idx[[domain_levels[d]]][1L],
                                domain_cols, drop = FALSE], as.character)),
          collapse = ":"
        )
        stop(
          sprintf(
            "min_n_domain = %g not achievable for domain '%s' (max total_n = %.0f at full budget)",
            min_n_domain,
            lab,
            full_total[d]
          ),
          call. = FALSE
        )
      }
    }
    min_fracs <- min_n_domain / full_total
    if (sum(min_fracs) > 1) {
      stop(
        sprintf(
          "min_n_domain = %g not achievable for all domains within budget",
          min_n_domain
        ),
        call. = FALSE
      )
    }
  }

  lower_bounds <- rep(1e-4, nd)
  if (!is.null(min_n_domain)) {
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
            if (!is.null(min_n_domain) && total_n_for(bres) < min_n_domain) {
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
      n_vec <- c(n_psu = bres$n1, n_per_psu = bres$n_per_psu)
    } else {
      n_vec <- c(
        n_psu = bres$n1,
        n_per_psu = bres$n_per_psu,
        n_per_ssu = bres$n_per_ssu
      )
    }
    .new_svyplan_cluster(
      n = n_vec,
      stages = stages,
      total_n = prod(n_vec),
      cv = bres$cv_achieved[bres$binding_idx],
      cost = bres$cost,
      params = list(stage_cost = stage_cost),
      indicators = NULL,
      detail = NULL,
      binding = p$labels[bres$binding_idx]
    )
  })
  names(results) <- domain_levels

  res <- .aggregate_cluster_domains(
    results,
    indicators,
    domain_cols,
    domain_levels,
    split_idx,
    stage_cost,
    budget,
    n_psu,
    n_per_psu,
    n_per_ssu,
    min_n_domain,
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
  icc,
  unit_relvar,
  var_ratio,
  resp_rate_psu,
  resp_rate,
  C1,
  C2,
  budget,
  n_psu,
  n_per_psu = NULL
) {
  nr <- length(cv_t)

  # Gross take in, realized take read, as in the CV-mode solver.
  cv_fn <- function(n1, take) {
    vapply(
      seq_len(nr),
      function(j) {
        mr <- take * resp_rate[j]
        sqrt(
          unit_relvar[j] *
            var_ratio[j] /
            (n1 * resp_rate_psu[j] * mr) *
            (1 + icc[j] * (mr - 1))
        )
      },
      numeric(1L)
    )
  }

  if (!is.null(n_per_psu)) {
    n_per_psu_opt <- n_per_psu
    n1_opt <- budget / (C1 + C2 * n_per_psu)
    if (n1_opt <= 0) {
      stop(
        "budget is too small for the given fixed cluster size",
        call. = FALSE
      )
    }
  } else if (is.null(n_psu)) {
    obj_fn <- function(n_per_psu) {
      n1 <- budget / (C1 + C2 * n_per_psu)
      if (n1 <= 0) {
        return(1e12)
      }
      max(cv_fn(n1, n_per_psu) / cv_t)
    }

    upper <- max(10, budget / (C1 + C2))
    opt <- optimize(obj_fn, interval = c(1, upper))
    n_per_psu_opt <- opt$minimum
    n1_opt <- budget / (C1 + C2 * n_per_psu_opt)
  } else {
    n1_opt <- n_psu
    n_per_psu_opt <- (budget - C1 * n_psu) / (C2 * n_psu)
    if (n_per_psu_opt <= 0) {
      stop(
        "budget is too small for the given fixed stage-1 size",
        call. = FALSE
      )
    }
  }

  cv_achieved <- cv_fn(n1_opt, n_per_psu_opt)
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    n_per_psu = n_per_psu_opt,
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
  icc_psu,
  icc_ssu,
  unit_relvar,
  var_ratio_psu,
  var_ratio_ssu,
  resp_rate_psu,
  resp_rate_ssu,
  resp_rate,
  C1,
  C2,
  C3,
  budget,
  n_psu,
  n_per_psu = NULL,
  n_per_ssu = NULL
) {
  nr <- length(cv_t)

  # Mirrors cv_achieved_fn() in the CV-mode solver: both stage takes are
  # gross and each row reads the takes it realizes.
  cv_fn <- function(n1, n_per_psu, n_per_ssu) {
    vapply(
      seq_len(nr),
      function(j) {
        mr <- n_per_psu * resp_rate_ssu[j]
        qr <- n_per_ssu * resp_rate[j]
        sqrt(
          unit_relvar[j] /
            (n1 * resp_rate_psu[j] * mr * qr) *
            (var_ratio_psu[j] *
              icc_psu[j] *
              mr *
              qr +
              var_ratio_ssu[j] * (1 + icc_ssu[j] * (qr - 1)))
        )
      },
      numeric(1L)
    )
  }

  n_free <- 3L - sum(!is.null(n_psu), !is.null(n_per_psu), !is.null(n_per_ssu))

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
          "optimal stage size reached the search upper bound (n_per_psu <= %.0f, n_per_ssu <= %.0f); result may be unreliable -- review stage costs and target CVs",
          upper_ps,
          upper_ss
        ),
        call. = FALSE
      )
    }
    n_per_psu_opt <- opt$par[1L]
    n_per_ssu_opt <- opt$par[2L]
    n1_opt <- budget / (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu_opt)
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
        n_per_ssu_opt <- 1
      } else {
        opt <- optimize(obj_fn_ss, interval = c(1, max_ss))
        n_per_ssu_opt <- opt$minimum
      }
      n_per_psu_opt <- ps_from_ss(n_per_ssu_opt)
      if (!is.finite(n_per_psu_opt) || n_per_psu_opt < 1) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
    } else if (!is.null(n_per_psu)) {
      obj_fn_ps <- function(ss) {
        n1 <- budget / (C1 + C2 * n_per_psu + C3 * n_per_psu * ss)
        if (n1 <= 0) return(1e12)
        max(cv_fn(n1, n_per_psu, ss) / cv_t)
      }
      upper <- max(10, budget / (C1 + C2 * n_per_psu + C3 * n_per_psu))
      opt <- optimize(obj_fn_ps, interval = c(1, upper))
      n_per_ssu_opt <- opt$minimum
      n_per_psu_opt <- n_per_psu
      n1_opt <- budget / (C1 + C2 * n_per_psu + C3 * n_per_psu * n_per_ssu_opt)
    } else {
      obj_fn_ss2 <- function(ps) {
        n1 <- budget / (C1 + C2 * ps + C3 * ps * n_per_ssu)
        if (n1 <= 0) return(1e12)
        max(cv_fn(n1, ps, n_per_ssu) / cv_t)
      }
      upper <- max(10, budget / (C1 + C2 + C3 * n_per_ssu))
      opt <- optimize(obj_fn_ss2, interval = c(1, upper))
      n_per_psu_opt <- opt$minimum
      n_per_ssu_opt <- n_per_ssu
      n1_opt <- budget / (C1 + C2 * n_per_psu_opt + C3 * n_per_psu_opt * n_per_ssu)
    }
    total_cost <- C1 * n1_opt +
      C2 * n1_opt * n_per_psu_opt +
      C3 * n1_opt * n_per_psu_opt * n_per_ssu_opt
  } else {
    if (!is.null(n_psu) && !is.null(n_per_psu)) {
      n_per_ssu_opt <- (budget - C1 * n_psu - C2 * n_psu * n_per_psu) /
        (C3 * n_psu * n_per_psu)
      if (n_per_ssu_opt <= 0) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      n1_opt <- n_psu
      n_per_psu_opt <- n_per_psu
    } else if (!is.null(n_psu) && !is.null(n_per_ssu)) {
      n_per_psu_opt <- (budget / n_psu - C1) / (C2 + C3 * n_per_ssu)
      if (n_per_psu_opt <= 0) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      n1_opt <- n_psu
      n_per_ssu_opt <- n_per_ssu
    } else {
      n1_opt <- budget / (C1 + C2 * n_per_psu + C3 * n_per_psu * n_per_ssu)
      if (n1_opt <= 0) {
        stop("budget is too small for the given fixed stage sizes", call. = FALSE)
      }
      n_per_psu_opt <- n_per_psu
      n_per_ssu_opt <- n_per_ssu
    }
    total_cost <- C1 * n1_opt +
      C2 * n1_opt * n_per_psu_opt +
      C3 * n1_opt * n_per_psu_opt * n_per_ssu_opt
  }

  cv_achieved <- cv_fn(n1_opt, n_per_psu_opt, n_per_ssu_opt)
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    n_per_psu = n_per_psu_opt,
    n_per_ssu = n_per_ssu_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = total_cost
  )
}

#' @rdname n_multi
#' @export
n_multi.svyplan_prec <- function(indicators, ...) {
  x <- indicators
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
  tgt <- x$params$indicators
  if ("prop_method" %in% names(dots)) {
    tgt$prop_method <- NA_character_
  }
  if ("resp_rate" %in% names(dots)) {
    tgt$resp_rate <- NA_real_
  }
  tgt$n <- NULL
  tgt$n_per_psu <- NULL
  tgt$n_per_ssu <- NULL

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
    indicators = tgt,
    domains = x$params$domain_cols,
    min_n_domain = x$params$min_n_domain,
    domain_sampling = x$params$domain_sampling %||% "separate",
    prop_method = x$params$prop_method %||% "wald"
  )
  do.call(n_multi.default, .roundtrip_args(args, dots, n_multi.default))
}

#' @rdname n_multi_cluster
#' @export
n_multi_cluster.svyplan_prec <- function(indicators, ...) {
  x <- indicators
  dots <- list(...)
  if (x$type != "multi" || !identical(x$params$design, "cluster")) {
    stop(
      "n_multi_cluster requires cluster precision from prec_multi_cluster()",
      call. = FALSE
    )
  }

  tgt <- x$params$indicators
  for (rate in intersect(
    c("resp_rate_psu", "resp_rate_ssu", "resp_rate"),
    names(dots)
  )) {
    tgt[[rate]] <- NA_real_
  }
  tgt$n <- NULL
  tgt$n_per_psu <- NULL
  tgt$n_per_ssu <- NULL

  stored_mode <- x$params$mode
  if (identical(stored_mode, "moe")) {
    tgt$moe <- x$detail$.moe
    tgt$cv <- NULL
  } else if (!is.null(x$detail) && ".cv" %in% names(x$detail)) {
    tgt$cv <- x$detail$.cv
    tgt$moe <- NULL
  }

  args <- list(
    indicators = tgt,
    stage_cost = x$params$stage_cost,
    domains = x$params$domain_cols,
    budget = x$params$budget,
    n_psu = x$params$n_psu,
    n_per_psu = x$params$n_per_psu,
    n_per_ssu = x$params$n_per_ssu,
    allocation = if (isTRUE(x$params$joint)) "joint" else "separate",
    min_n_domain = x$params$min_n_domain,
    domain_sampling = x$params$domain_sampling %||% "separate",
    fixed_cost = x$params$fixed_cost %||% 0
  )
  do.call(
    n_multi_cluster.default,
    .roundtrip_args(args, dots, n_multi_cluster.default)
  )
}
