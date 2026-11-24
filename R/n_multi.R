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
#'       30% stunting), `var` for a continuous variable's population
#'       variance, **or** `r` for a ratio of two means, with its moments
#'       `cv_num`, `cv_den` and `component_cor`. Each row must use exactly
#'       one.
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
#'   independently for each domain combination. A domain column cannot be
#'   called `name` or start with a dot, the names of result columns.
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
#'   `n` as `max(.n / share)`. That is an *expected yield*. A sample of that
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
#'   interval half-width, so those rows need `moe` or `rmoe` rather than
#'   `cv`.
#' @param resp_rate Default expected response rate among eligible ultimate
#'   units, in (0, 1\]. Used for rows whose `resp_rate` column is absent or
#'   `NA`, and a non-missing row value overrides it. See the nonresponse
#'   section of [svyplan-package] for what the rate counts and what this
#'   adjustment does and does not claim.
#' @param eligibility_rate Default share of issued units in scope, in (0, 1\].
#'   Used for rows whose `eligibility_rate` column is absent or `NA`, and a
#'   non-missing row value overrides it. With a finite `N`, a row's population
#'   in scope is `N * eligibility_rate`. See the nonresponse section of
#'   [svyplan-package].
#' @param eligible_per_unit Default eligible units per sampled unit, a positive
#'   number such as married women per household, which can exceed 1. Used
#'   for rows whose `eligible_per_unit` column is absent or `NA`, and a
#'   non-missing row value overrides it. Proportion and mean rows. See the
#'   `eligible_per_unit` column under Details.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object. The output class is the same with or without
#'   domains.
#'
#'   **Without domains**, the object contains:
#'   \describe{
#'     \item{`n`}{The sample size required by the binding indicator.}
#'     \item{`detail`}{One row per indicator: `name`, its requirement `.n`,
#'       `.cv_target`, `.cv_achieved` at `n`, `.cv_field` at `ceiling(n)`,
#'       `.binding`, `.constraint` (what sets the row's requirement:
#'       `"precision"`, `"min_cases"` or `"min_respondents"`), and the
#'       target as stated, `.target_scale` (`"moe"`, `"rmoe"` or `"cv"`)
#'       and `.target`. A row with only a respondent floor has `NA`
#'       precision.}
#'     \item{`binding`}{Name or index of the binding (most demanding) indicator.}
#'     \item{`binding_cause`}{What made it bind, one of the `.constraint`
#'       values.}
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
#'     \item{`n_domain_max`}{The largest single domain requirement. Under
#'       `domain_sampling = "natural"`, this need not be the domain that
#'       determines `n`, which is selected after division by `.share`.}
#'     \item{`binding`}{The binding indicator within the domain that
#'       determines `n`: the largest domain quota under `"separate"`, or the
#'       largest share-adjusted requirement under `"natural"`. Use
#'       `$domains` to identify the domain itself.}
#'     \item{`domains`}{Data frame with one row per domain, including
#'       domain variables, `.n`, `.binding`, and `.share` under
#'       `"natural"`.}
#'     \item{`detail`}{One row per indicator and domain, the domain
#'       variables first, then `name`, `.n`, `.cv_target`, `.binding`,
#'       `.constraint`, `.target_scale` and `.target`. These are each
#'       domain's own requirements, before `min_n_domain` or natural
#'       incidence raises its size. [summary()] reads what each row achieves
#'       at the fielded size.}
#'     \item{`binding_cause`}{As without domains, or `"min_n_domain"` when
#'       that floor sets the binding domain's quota.}
#'   }
#'
#' @details
#' ## Building the indicators data frame
#'
#' Each row of `indicators` represents one survey indicator. The two key
#' decisions per row are:
#'
#' 1. **Type of indicator**: is it a proportion (binary variable like
#'    "stunted yes/no"), a mean (continuous variable like "household
#'    expenditure"), or a ratio of two means (such as expenditure per
#'    adult equivalent)? This determines whether you fill the `p`, `var`
#'    or `r` column.
#' 2. **Precision target**: do you want an absolute margin of error
#'    (`moe`, e.g. +/- 5 percentage points), the same margin stated
#'    relative to the estimand (`rmoe`, e.g. 12 percent of the
#'    proportion), or a relative coefficient of
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
#' Rows with `p` are treated as proportions, rows with `var` (and
#' `p = NA`) as means, and rows with `r` as ratios. A row carries only one
#' of `p`, `var` and `r`.
#'
#' A row may instead carry only a floor on the number of respondents,
#' `min_respondents`, with no `p`, `var` or target. It stands for an
#' estimate whose precision is not modeled, such as the fertility rates a
#' DHS sizes by a count of interviews with women. Its `resp_rate` and
#' `eligible_per_unit` turn the count into sampled units, and it raises the
#' size the other rows set. Every domain needs at least one row with an
#' estimate and a target.
#'
#' ## Column reference
#'
#' \describe{
#'   \item{`name`}{Indicator label (optional). If omitted, each row is
#'     labeled by its row number in `indicators`, with domains as without.
#'     When given, every row needs a name, unique within a domain.}
#'   \item{`p`}{Expected proportion, in (0, 1). Use this for binary
#'     indicators such as prevalences or coverage rates. The value is
#'     your best prior guess (e.g. from a previous survey or literature).
#'     One of `p`, `var` or `r` per row.}
#'   \item{`var`}{Population variance of a continuous indicator. Use
#'     this for means (e.g. income, expenditure, weight). One of `p`,
#'     `var` or `r` per row.}
#'   \item{`r`}{Anticipated ratio of two means, marking a ratio row, sized
#'     as [n_ratio()] sizes it. It requires `cv_num`, `cv_den` and
#'     `component_cor` alongside it. One of `p`, `var` or `r` per row.}
#'   \item{`cv_num`, `cv_den`}{Coefficients of variation of a ratio row's
#'     numerator and denominator, strictly positive.}
#'   \item{`component_cor`}{Correlation between a ratio row's numerator and
#'     denominator across units, in \eqn{[-1, 1]}.}
#'   \item{`mean`}{Population mean, finite. It may be negative. It is
#'     required when `var` is used with `cv` because CV = SE / abs(mean),
#'     and it must then be non-zero, as it must in multistage mode, which
#'     works on the cv scale. With a `moe` target a zero mean is accepted
#'     and the row's `cv` and `rmoe` are `NA`.}
#'   \item{`moe`}{Margin of error, the half-width of the confidence
#'     interval you want. For proportions, this is on the probability
#'     scale (e.g. 0.05 for +/- 5 percentage points). For means,
#'     it is in the same units as the variable (e.g. 10 dollars).}
#'   \item{`rmoe`}{Margin of error relative to the row's estimand, so
#'     0.12 asks for a half-width of 12 percent of `p` or of `abs(mean)`.
#'     It is converted to `moe` on ingestion, so the row needs `p` or
#'     `mean`, and it is read under the row's own `prop_method`.}
#'   \item{`cv`}{Target coefficient of variation (relative standard
#'     error). For example, 0.10 means the SE should be at most 10%
#'     of the estimate.}
#'   \item{`alpha`}{Significance level for the confidence interval
#'     (default 0.05, giving a 95 percent CI). `2 * pnorm(-2)` gives the
#'     multiplier of 2 that the DHS and MICS manuals use in place of 1.96.}
#'   \item{`deff`}{Design effect multiplier (default 1). Set > 1 to inflate
#'     the sample size for complex
#'     designs (e.g. 1.5 for a cluster design). A published DEFT is
#'     squared first. See [varcomp()] for that and for removing the
#'     published survey's weighting loss.}
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
#'     Single-stage `n_multi()` only, because a multistage design sizes
#'     stages against a cost and has no single total for a count to raise.}
#'   \item{`min_respondents`}{Minimum expected number of eligible
#'     respondents, `n * resp_rate * eligibility_rate * eligible_per_unit`,
#'     a floor on the row's size. Unlike `min_cases`, it counts every respondent of the row's
#'     group, not those with the attribute, and it applies to any row. A row
#'     with this floor and no estimand is the floor-only row described
#'     above. `NA` (the default) sets no floor. Single-stage `n_multi()`
#'     only, as for `min_cases`.}
#'   \item{`unit_relvar`}{Unit relvariance. If omitted, derived
#'     automatically from `p` (as `(1 - p) / p`) or from
#'     `var` / `mean^2`. A ratio row derives it from its moments, and a value
#'     given there must agree with them.}
#'   \item{`resp_rate`}{Expected response rate at the ultimate unit, in
#'     (0, 1\]. Default 1 (no adjustment). A value of 0.90 inflates the
#'     sample size by `1 / 0.90` to compensate for 10 percent
#'     nonresponse. It means the same thing in a cluster design, where
#'     [n_cluster()] additionally takes `resp_rate_1` for whole clusters
#'     that cannot be worked and `resp_rate_2` for second-stage units in a
#'     three-stage design. Naming a stage the design does not have is an
#'     error rather than a column carried along and ignored.}
#'   \item{`eligibility_rate`}{Share of the row's issued units in scope, in
#'     (0, 1\]. It inflates the row's size by `1 / eligibility_rate`, and
#'     with a finite `N` the population in scope is `N * eligibility_rate`.
#'     `NA` takes the `eligibility_rate` argument.}
#'   \item{`eligible_per_unit`}{Eligible units per sampled unit for the row's
#'     estimate, such as married women or children aged 12 to 23 months per
#'     household. The row is sized as [n_prop()] or [n_mean()] size it with
#'     `eligible_per_unit`. Sizes stay in sampled units, and the variance reads
#'     `N * eligibility_rate * eligible_per_unit` eligible units and
#'     `resp_rate * eligibility_rate * eligible_per_unit` respondents per
#'     sampled unit. Each row
#'     can carry its own group.
#'     Proportion and mean rows, so a ratio row takes 1. The
#'     several-indicators mode of [n_cluster()] reads the same column.}
#' }
#'
#' Domain columns are specified via the `domains` parameter. When domains
#' are present, sizing runs independently for each domain combination, and
#' `domain_sampling` decides how those requirements combine into the one
#' number `n` reports. The two readings answer different questions and give
#' different sizes, so the choice belongs to the design rather than to a
#' default. Separate quotas need every requirement met and therefore their
#' sum, while natural domains need a sample large enough that the rarest
#' demanding domain turns up often enough. Neither is the largest single
#' requirement, which `n_domain_max` reports separately.
#'
#' Columns beyond those listed are carried along and ignored, so the table
#' can keep questionnaire modules, sources, or other bookkeeping. The
#' exception is a name the rest of the package would lead you to expect
#' here: `mu` (tables use `mean` across the package, while scalar functions
#' such as [n_mean()] take `mu`), `indicator` or `label` (it takes `name`),
#' `method` (it takes `prop_method`), and `icc` or `var_ratio` (which are
#' per stage here, so `icc_1` and `var_ratio_1`). Those are rejected
#' with the name they should have carried rather than silently ignored. A
#' column named through `domains` is exempt, so a domain may be called
#' `label` or `method`.
#'
#' The dispersion of a continuous indicator may be given as `var` or as
#' `sd`, whichever the source reports, exactly as in [n_mean()]. They are
#' different quantities rather than two names for one, so supplying both
#' in the same table is an error rather than a preference. `sd` is squared
#' on the way in and everything downstream reads `var`.
#'
#' `n_multi()` computes sample size per indicator by delegating proportion
#' rows to [n_prop()], mean rows to [n_mean()] and ratio rows to
#' [n_ratio()], then takes the maximum per domain. Use `prop_method` or a
#' `indicators$prop_method` column to choose `"wald"`, `"wilson"`,
#' `"logodds"` or `"beta"` for proportion rows.
#'
#' ## A domain that takes most of its PSUs
#'
#' A `deff` column carries the clustering of a design whose first stage is
#' treated as drawn with replacement, and `N` corrects only for the units.
#' Nothing here knows how many PSUs a domain holds. A small domain, such as
#' an island district of 16 enumeration areas, can then be given more PSUs
#' than it has, and a sample that takes most of them is assessed as if
#' they were a few among many. [n_alloc()] and [prec_alloc()] read `N_1`
#' in the frame and, with `fpc = "stage"`, apply the stage-by-stage
#' correction. They size and assess such a domain, as the Details of
#' [n_alloc()] describe.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer.
#'
#' @family multi-indicator functions
#' @seealso [n_prop()] and [n_mean()] for single-indicator sizing,
#'   [n_alloc()] to split a multi-indicator size across strata or domains,
#'   [n_cluster()] for a multistage cluster design over the same indicator
#'   table, and [prec_multi()] for the inverse.
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
#' # MICS-style: state precision as a relative margin of error
#' targets_rmoe <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   rmoe = 0.12,
#'   deff = c(2.0, 1.5, 2.5)
#' )
#' n_multi(targets_rmoe)
#'
#' # Continuous indicators: 'var' for the dispersion, 'mean' for the mean.
#' # A CV target needs 'mean', because a relative standard error is relative
#' # to something. A MOE target does not.
#' targets_mean <- data.frame(
#'   name = c("expenditure", "hh_size"),
#'   var  = c(250000, 4.0),
#'   mean = c(1200, 5.4),
#'   cv   = c(0.05, 0.03)
#' )
#' n_multi(targets_mean)
#'
#' # Proportions and means in one table, sized to a common CV
#' targets_both <- data.frame(
#'   name = c("stunting", "expenditure"),
#'   p    = c(0.30, NA),
#'   var  = c(NA, 250000),
#'   mean = c(NA, 1200),
#'   cv   = c(0.08, 0.05)
#' )
#' n_multi(targets_both)
#'
#' # Rare proportion: use Wilson globally in simple mode
#' n_multi(indicators[3, , drop = FALSE], prop_method = "wilson")
#'
#' # Korn-Graubard for an indicator whose expected count is small, with the
#' # degrees of freedom of the planned variance estimator. 'df' is read on
#' # every row that sets it, and NA keeps the normal quantile on the other.
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
#' # One indicator in several domains, each with its own design effect. The
#' # 'name' column is optional with a single indicator.
#' poverty <- data.frame(
#'   district = c("North", "Centre", "South"),
#'   p = c(0.31, 0.45, 0.52),
#'   moe = 0.10,
#'   deff = c(4.4, 3.1, 2.8)
#' )
#' n_multi(poverty, domains = "district")
#'
#' # A household survey whose estimates are over persons. Each row counts
#' # its own group per household, and a row with only 'min_respondents'
#' # asks for 800 interviewed women whatever the targets need.
#' dhs <- data.frame(
#'   region = rep(c("Coast", "Inland"), 3),
#'   name = rep(c("mcpr", "vaccination", "women"), each = 2),
#'   p = c(0.22, 0.08, 0.77, 0.85, NA, NA),
#'   cv = c(0.10, 0.10, 0.10, 0.10, NA, NA),
#'   deff = c(1.03, 2.71, 1.32, 2.37, NA, NA),
#'   resp_rate = 0.88,
#'   eligible_per_unit = c(0.59, 2.30, 0.09, 0.47, 1.13, 2.89),
#'   min_respondents = c(NA, NA, NA, NA, 800, 800)
#' )
#' n_multi(dhs, domains = "region")
#'
#' @export
n_multi <- function(indicators, ...) {
  if (!missing(indicators)) {
    .refuse_result_as_input(indicators, "indicators", "n_multi")
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
  eligibility_rate = 1,
  eligible_per_unit = 1,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_multi.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_multi.default, c(.plan, list(...))))
  }
  .check_multi_split_args(list(...), "n_cluster(<table>, ...)")
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(indicators, "indicators")
  domain_sampling <- match.arg(domain_sampling)
  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }
  if (!is.null(min_n_domain)) {
    if (
      !is.numeric(min_n_domain) ||
        length(min_n_domain) != 1L ||
        is.na(min_n_domain) ||
        min_n_domain <= 0
    ) {
      stop("'min_n_domain' must be a positive numeric scalar", call. = FALSE)
    }
  }
  if (!is.null(min_n_domain) && is.null(domains)) {
    stop("'min_n_domain' needs 'domains'", call. = FALSE)
  }
  check_resp_rate(resp_rate)
  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)
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
  provenance <- .target_provenance(indicators)
  indicators <- .indicators_moe_from_rmoe(indicators, domains)
  .check_unique_indicators(indicators, domains)
  floor_only <- .floor_only_rows(indicators)
  .check_respondent_floors(indicators, floor_only, domains)
  .validate_common_columns(indicators)
  info <- .validate_targets(
    indicators[!floor_only, , drop = FALSE],
    FALSE,
    domains = domains
  )
  indicators <- .fill_defaults(
    indicators,
    FALSE,
    prop_method = prop_method,
    resp_rate = resp_rate
  )
  indicators <- .fill_eligible_per_unit(indicators, eligible_per_unit)
  indicators <- .fill_eligibility_rate(indicators, eligibility_rate)

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

  res <- if (length(domain_cols) == 0L) {
    .n_multi_simple(
      indicators,
      domain_cols = domain_cols,
      mode = mode,
      prop_method = prop_method
    )
  } else {
    .n_multi_domains(
      indicators,
      domain_cols,
      min_n_domain,
      mode = mode,
      prop_method = prop_method,
      domain_sampling = domain_sampling
    )
  }
  .attach_target_provenance(res, provenance, indicators, domain_cols)
}

#' The target each row was stated with, before normalization
#'
#' Normalization rewrites an `rmoe` target as a `moe`, so the scale the user
#' wrote is recorded first for the summary to report it.
#' @keywords internal
#' @noRd
.target_provenance <- function(indicators) {
  nr <- nrow(indicators)
  scale <- rep(NA_character_, nr)
  value <- rep(NA_real_, nr)
  for (col in intersect(c("cv", "moe", "rmoe"), names(indicators))) {
    set <- !is.na(indicators[[col]])
    scale[set] <- col
    value[set] <- indicators[[col]][set]
  }
  list(scale = scale, value = value)
}

#' Write the stated targets onto the per-row detail
#'
#' With domains the detail stacks the domains in the order the solver split
#' them, so the rows are reordered the same way.
#' @keywords internal
#' @noRd
.attach_target_provenance <- function(res, provenance, indicators,
                                      domain_cols) {
  if (is.null(res$detail)) {
    return(res)
  }
  ord <- seq_len(nrow(indicators))
  if (length(domain_cols) > 0L) {
    key <- .domain_key(indicators, domain_cols)
    ord <- unlist(split(ord, factor(key, levels = unique(key))),
                  use.names = FALSE)
  }
  res$detail$.target_scale <- provenance$scale[ord]
  res$detail$.target <- provenance$value[ord]
  res
}

#' Per-row eligibility rate, materialized only when one is set
#' @keywords internal
#' @noRd
.fill_eligibility_rate <- function(indicators, eligibility_rate) {
  has_col <- "eligibility_rate" %in% names(indicators)
  if (!has_col && eligibility_rate == 1) {
    return(indicators)
  }
  e <- if (has_col) indicators$eligibility_rate else
    rep(NA_real_, nrow(indicators))
  e[is.na(e)] <- eligibility_rate
  indicators$eligibility_rate <- e
  indicators
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

#' Per-row eligible rate, materialized only when one is set
#'
#' A table without the column under the default of 1 is left as it is, so a
#' default result keeps its columns. Otherwise every row carries its rate,
#' which a round trip through [prec_multi()] reads back.
#' @keywords internal
#' @noRd
.fill_eligible_per_unit <- function(indicators, eligible_per_unit) {
  has_col <- "eligible_per_unit" %in% names(indicators)
  if (!has_col && eligible_per_unit == 1) {
    return(indicators)
  }
  e <- if (has_col) indicators$eligible_per_unit else rep(NA_real_, nrow(indicators))
  if (!is.numeric(e) && !all(is.na(e))) {
    stop("'eligible_per_unit' values must be positive and finite", call. = FALSE)
  }
  e <- as.numeric(e)
  e[is.na(e)] <- eligible_per_unit
  if (any(!is.finite(e)) || any(e <= 0)) {
    stop("'eligible_per_unit' values must be positive and finite", call. = FALSE)
  }
  is_ratio <- .is_ratio_row(indicators)
  bad <- which(e != 1 & is_ratio)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "'eligible_per_unit' applies to proportion and mean rows, as in n_prop() and n_mean(), and to rows with only 'min_respondents'. Row(s) %s carry a ratio, so set their 'eligible_per_unit' to 1",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  indicators$eligible_per_unit <- e
  indicators
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
  stage_rr <- grep("^resp_rate_[0-9]+$", names(indicators), value = TRUE)
  if (!multistage) {
    if (length(stage_rr) == 0L) {
      return(invisible(TRUE))
    }
    stop(
      sprintf(
        "column(s) %s name sampling stages this design does not have; a single-stage design loses ultimate units, so its rate is the plain 'resp_rate'",
        paste(sQuote(stage_rr), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  allowed <- paste0("resp_rate_", seq_len(stages - 1L))
  bad <- setdiff(stage_rr, allowed)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "column(s) %s are not applicable for %d-stage designs (allowed stage rate columns: %s, plus 'resp_rate' for stage %d)",
        paste(sQuote(bad), collapse = ", "),
        stages,
        paste(sQuote(allowed), collapse = ", "),
        stages
      ),
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
    c(
      "stage_cost",
      "budget",
      "stage_n",
      "allocation",
      "fixed_cost"
    )
  )
  if (length(moved) > 0L) {
    stop(
      sprintf(
        "cluster argument%s %s belong%s to %s",
        if (length(moved) > 1L) "s" else "",
        paste(sQuote(moved), collapse = ", "),
        if (length(moved) > 1L) "" else "s",
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
.validate_targets <- function(
  indicators,
  multistage,
  domains = NULL,
  stages = NULL,
  context = "n_cluster(<table>, ...)"
) {
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

  .check_indicator_columns(indicators, domains, eligible = !multistage,
                           stages = if (multistage) stages %||% 2L,
                           context = if (multistage) context else "n_multi()")
  .check_resp_rate_column(indicators, multistage, stages %||% 2L)

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)
  .check_ratio_rows(indicators)
  if (!has_p && !has_var && !any(is_ratio)) {
    stop("'indicators' must contain a 'p', 'var', or 'r' column", call. = FALSE)
  }

  has_moe <- "moe" %in% names(indicators)
  has_cv <- "cv" %in% names(indicators)
  if (!has_moe && !has_cv) {
    stop(
      "'indicators' must contain a 'moe', 'rmoe', or 'cv' column",
      call. = FALSE
    )
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

  # Each row is exactly one estimand: a proportion, a mean, or a ratio.
  .check_estimand_markers(indicators)

  given <- (has_moe & !is.na(indicators[["moe"]] %||% NA)) |
    (has_cv & !is.na(indicators[["cv"]] %||% NA))
  if (any(!given)) {
    stop("each row must have either 'moe' or 'cv' specified", call. = FALSE)
  }
  if (has_moe && has_cv) {
    both_set <- !is.na(indicators$moe) & !is.na(indicators$cv)
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
      stop(
        "multistage mode requires a 'cv', 'moe', or 'rmoe' column in indicators",
        call. = FALSE
      )
    }
    if (has_moe && any(!is.na(indicators$moe))) {
      moe_rows <- !is.na(indicators$moe)
      var_moe <- has_var & moe_rows & !is.na(indicators[["var"]])
      if (any(var_moe)) {
        if (
          !"mean" %in% names(indicators) || any(is.na(indicators$mean[var_moe]))
        ) {
          stop(
            "'mean' is required to convert 'moe' to 'cv' for mean indicators ",
            "in multistage mode",
            call. = FALSE
          )
        }
      }
    }
    n_b <- (stages %||% 2L) - 1L
    req_icc <- paste0("icc_", seq_len(n_b))
    missing_icc <- setdiff(req_icc, names(indicators))
    if (length(missing_icc) > 0L) {
      stop(
        sprintf(
          "%d-stage mode requires column(s) %s in indicators",
          stages %||% 2L, paste(sQuote(missing_icc), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (has_cv) {
      cv_vals <- indicators$cv[!is.na(indicators$cv)]
      if (
        length(cv_vals) > 0L && (any(cv_vals <= 0) || any(!is.finite(cv_vals)))
      ) {
        stop("'cv' values must be positive and finite", call. = FALSE)
      }
    }
    for (col in req_icc) {
      vals <- indicators[[col]]
      if (anyNA(vals) || any(!is.finite(vals))) {
        stop(sprintf("'%s' must contain finite non-missing values", col),
             call. = FALSE)
      }
      if (any(vals < 0 | vals > 1)) {
        stop(sprintf("'%s' values must be in [0, 1]", col), call. = FALSE)
      }
      .check_cluster_icc_open(vals, context = context)
    }
  }

  if (has_var && has_cv) {
    needs_mean <- !is.na(indicators$var) & !is.na(indicators$cv)
    if (any(needs_mean)) {
      if (!"mean" %in% names(indicators) || any(is.na(indicators$mean[needs_mean]))) {
        stop(
          "'mean' is required when 'var' and 'cv' are specified",
          call. = FALSE
        )
      }
    }
  }

  if ("mean" %in% names(indicators)) {
    mean_vals <- indicators$mean[!is.na(indicators$mean)]
    if (any(!is.finite(mean_vals))) {
      stop("'mean' values must be finite", call. = FALSE)
    }
    # A zero level serves an absolute target. A cv target, and the cv a
    # multistage solve converts every target to, have no scale there.
    zero <- !is.na(indicators$mean) & indicators$mean == 0
    relative <- if (multistage) {
      rep(TRUE, nrow(indicators))
    } else if (has_cv) {
      !is.na(indicators$cv)
    } else {
      rep(FALSE, nrow(indicators))
    }
    if (any(zero & relative)) {
      stop(
        if (multistage) {
          "'mean' must be non-zero for a mean indicator in multistage mode, which works on the cv scale"
        } else {
          "'mean' must be non-zero on a row with a 'cv' target"
        },
        call. = FALSE
      )
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
    n_b <- (stages %||% 2L) - 1L
    for (s in seq_len(n_b)) {
      col <- paste0("var_ratio_", s)
      if (col %in% names(indicators)) {
        vals <- indicators[[col]][!is.na(indicators[[col]])]
        if (any(vals <= 0) || any(!is.finite(vals))) {
          stop(sprintf("'%s' values must be positive and finite", col),
               call. = FALSE)
        }
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
    n_b <- stages - 1L
    if (!"var_ratio_1" %in% names(indicators)) {
      indicators$var_ratio_1 <- 1
    } else {
      indicators$var_ratio_1[is.na(indicators$var_ratio_1)] <- 1
    }
    if (n_b >= 2L) {
      for (s in 2L:n_b) {
        col <- paste0("var_ratio_", s)
        prev_vr <- indicators[[paste0("var_ratio_", s - 1L)]]
        prev_icc <- indicators[[paste0("icc_", s - 1L)]]
        implied <- prev_vr * (1 - prev_icc)
        if (!col %in% names(indicators)) {
          indicators[[col]] <- implied
        } else {
          miss <- is.na(indicators[[col]])
          indicators[[col]][miss] <- implied[miss]
        }
      }
    }
  }

  # Stage response rates: resp_rate_1..resp_rate_{S-1} plus resp_rate (stage S)
  if (multistage) {
    rr_vec <- .expand_stage_resp_rate(resp_rate, stages)
    rate_names <- c(paste0("resp_rate_", seq_len(stages - 1L)), "resp_rate")
    for (s in seq_len(stages)) {
      col <- rate_names[s]
      if (!col %in% names(indicators)) {
        indicators[[col]] <- rr_vec[s]
      } else {
        indicators[[col]][is.na(indicators[[col]])] <- rr_vec[s]
      }
    }
  } else {
    if (!"resp_rate" %in% names(indicators)) {
      indicators$resp_rate <- resp_rate[1L]
    } else {
      indicators$resp_rate[is.na(indicators$resp_rate)] <- resp_rate[1L]
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
  indicators$unit_relvar <- .derive_unit_relvar(
    indicators,
    require_all = multistage
  )

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
  has_mean <- "mean" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)

  for (i in which(needs)) {
    if (is_ratio[i]) {
      rv[i] <- .ratio_unit_relvar(
        indicators$r[i],
        indicators$cv_num[i],
        indicators$cv_den[i],
        indicators$component_cor[i]
      )
    } else if (has_p && !is.na(indicators$p[i])) {
      rv[i] <- (1 - indicators$p[i]) / indicators$p[i]
    } else if (has_var && !is.na(indicators$var[i])) {
      # A zero mean has no relvariance. Only an absolute target reaches
      # here with one, and that target reads 'var' directly.
      if (has_mean && !is.na(indicators$mean[i]) && indicators$mean[i] != 0) {
        rv[i] <- indicators$var[i] / indicators$mean[i]^2
      } else if (require_all) {
        stop(
          sprintf(
            "row %d: 'mean' is required to derive 'unit_relvar' from 'var'",
            i
          ),
          call. = FALSE
        )
      }
    }
  }

  if (require_all && anyNA(rv)) {
    stop(
      sprintf(
        "row(s) %s: could not derive 'unit_relvar' - provide 'p', or 'var' and 'mean', or 'unit_relvar' directly",
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
#' moves. Under the other three it is what makes `prop_method` reach the
#' multistage path at all.
#'
#' Rows that already have cv are left unchanged.
#' @keywords internal
#' @noRd
.convert_moe_to_cv <- function(indicators) {
  has_moe <- "moe" %in% names(indicators)
  if (!has_moe) {
    return(indicators)
  }

  moe_rows <- !is.na(indicators$moe)
  if (!any(moe_rows)) {
    return(indicators)
  }

  if (!"cv" %in% names(indicators)) {
    indicators$cv <- NA_real_
  }

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  has_method <- "prop_method" %in% names(indicators)

  is_ratio <- .is_ratio_row(indicators)

  for (i in which(moe_rows)) {
    df_i <- .row_df(indicators, i)
    if (is_ratio[i]) {
      # Exact under the first-order model, the same conversion the mean row
      # gets, on the ratio's own scale.
      z <- .q_alpha(indicators$alpha[i], df_i)
      indicators$cv[i] <- indicators$moe[i] / (z * abs(indicators$r[i]))
    } else if (has_p && !is.na(indicators$p[i])) {
      p_i <- indicators$p[i]
      method_i <- if (has_method) indicators$prop_method[i] else "wald"
      n_eff <- .n_prop_effective(
        p_i,
        indicators$moe[i],
        indicators$alpha[i],
        method_i,
        df_i
      )
      indicators$cv[i] <- sqrt(p_i * (1 - p_i) / n_eff) / p_i
    } else if (has_var && !is.na(indicators$var[i])) {
      z <- .q_alpha(indicators$alpha[i], df_i)
      indicators$cv[i] <- indicators$moe[i] / (z * abs(indicators$mean[i]))
    }
  }

  indicators
}

#' Simple mode: compute n per indicator, take max
#' @keywords internal
#' @noRd
.n_multi_simple <- function(
  indicators,
  domain_cols = character(0),
  mode = "moe",
  prop_method = "wald"
) {
  simple <- .compute_simple_n(indicators)
  n_vec <- simple$n
  cv_target_vec <- simple$cv_target

  # A case floor raises the row's own size, so the binding row is chosen
  # among sizes that already satisfy both constraints. `.cv_target` keeps
  # its precision meaning; `.cv_achieved` is read at the size that wins.
  constraint <- rep("precision", nrow(indicators))
  floor_n <- .multi_min_cases_n(indicators)
  if (!is.null(floor_n)) {
    constraint[!is.na(floor_n) & (is.na(n_vec) | floor_n > n_vec)] <-
      "min_cases"
    n_vec <- pmax(n_vec, floor_n, na.rm = TRUE)
  }
  resp_floor_n <- .multi_min_respondents_n(indicators)
  if (!is.null(resp_floor_n)) {
    constraint[!is.na(resp_floor_n) &
                 (is.na(n_vec) | resp_floor_n > n_vec)] <- "min_respondents"
    n_vec <- pmax(n_vec, resp_floor_n, na.rm = TRUE)
  }

  idx <- which.max(n_vec)
  n_max <- n_vec[idx]
  .check_multi_attainable(n_max, indicators)

  labels <- if ("name" %in% names(indicators)) {
    indicators$name
  } else {
    seq_len(nrow(indicators))
  }
  binding_name <- labels[idx]

  has_p <- "p" %in% names(indicators)
  has_mean <- "mean" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)
  elig <- .indicator_rate(indicators, "eligible_per_unit", nrow(indicators))
  scope <- .indicator_rate(indicators, "eligibility_rate", nrow(indicators))
  floor_only <- .floor_only_rows(indicators)
  cv_at <- function(n_eval) vapply(
    seq_len(nrow(indicators)),
    function(i) {
      if (floor_only[i]) {
        return(NA_real_)
      }
      prec <- if (is_ratio[i]) {
        .prec_engine_ratio(
          indicators$r[i],
          indicators$unit_relvar[i],
          n_eval,
          indicators$alpha[i],
          indicators$N[i] * scope[i],
          indicators$deff[i],
          indicators$resp_rate[i] * scope[i],
          .row_df(indicators, i)
        )
      } else if (has_p && !is.na(indicators$p[i])) {
        .prec_engine_prop(
          indicators$p[i],
          n_eval,
          indicators$alpha[i],
          indicators$N[i] * scope[i] * elig[i],
          indicators$deff[i],
          indicators$resp_rate[i] * scope[i] * elig[i],
          indicators$prop_method[i],
          .row_df(indicators, i)
        )
      } else {
        .prec_engine_mean(
          indicators$var[i],
          if (has_mean) indicators$mean[i] else NULL,
          n_eval,
          indicators$alpha[i],
          indicators$N[i] * scope[i] * elig[i],
          indicators$deff[i],
          indicators$resp_rate[i] * scope[i] * elig[i],
          .row_df(indicators, i)
        )
      }
      prec$cv
    },
    numeric(1L)
  )

  cv_achieved <- cv_at(n_max)
  detail <- data.frame(
    name = labels,
    .n = n_vec,
    .cv_target = cv_target_vec,
    .cv_achieved = cv_achieved,
    .cv_field = if (.whole_count(n_max) == n_max) cv_achieved else
      .muffle_census(cv_at(.whole_count(n_max))),
    .binding = seq_len(nrow(indicators)) == idx,
    .constraint = constraint
  )

  res <- .new_svyplan_multi(
    n = n_max,
    params = list(
      domain_cols = domain_cols,
      mode = mode,
      prop_method = prop_method
    ),
    indicators = indicators,
    detail = detail,
    binding = binding_name
  )
  res$binding_cause <- constraint[idx]
  res
}

#' Compute per-indicator n by delegating to the single-indicator engines
#'
#' Returns list with `n` (per-indicator sample size) and `cv_target` (the CV
#' each indicator achieves at its own n. It equals the target CV for cv-mode and
#' the CV implied by the target MOE for moe-mode).
#' @keywords internal
#' @noRd
.compute_simple_n <- function(indicators) {
  nr <- nrow(indicators)
  n_vec <- numeric(nr)
  cv_vec <- numeric(nr)

  has_p <- "p" %in% names(indicators)
  has_moe <- "moe" %in% names(indicators)
  has_mean <- "mean" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)
  elig <- .indicator_rate(indicators, "eligible_per_unit", nr)
  scope <- .indicator_rate(indicators, "eligibility_rate", nr)
  floor_only <- .floor_only_rows(indicators)

  for (i in seq_len(nr)) {
    if (floor_only[i]) {
      n_vec[i] <- NA_real_
      cv_vec[i] <- NA_real_
      next
    }
    is_prop <- has_p && !is.na(indicators$p[i])
    use_moe <- has_moe && !is.na(indicators$moe[i])

    if (is_ratio[i]) {
      res_i <- n_ratio.default(
        r = indicators$r[i],
        cv_num = indicators$cv_num[i],
        cv_den = indicators$cv_den[i],
        component_cor = indicators$component_cor[i],
        moe = if (use_moe) indicators$moe[i] else NULL,
        cv = if (use_moe) NULL else indicators$cv[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        eligibility_rate = scope[i],
        df = .row_df(indicators, i)
      )
    } else if (is_prop) {
      res_i <- n_prop.default(
        p = indicators$p[i],
        moe = if (use_moe) indicators$moe[i] else NULL,
        cv = if (use_moe) NULL else indicators$cv[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        eligibility_rate = scope[i],
        eligible_per_unit = elig[i],
        method = indicators$prop_method[i],
        df = .row_df(indicators, i)
      )
    } else {
      res_i <- n_mean.default(
        var = indicators$var[i],
        mu = if (has_mean && !is.na(indicators$mean[i])) indicators$mean[i] else NULL,
        moe = if (use_moe) indicators$moe[i] else NULL,
        cv = if (use_moe) NULL else indicators$cv[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        eligibility_rate = scope[i],
        eligible_per_unit = elig[i],
        df = .row_df(indicators, i)
      )
    }
    n_vec[i] <- res_i$n
    cv_vec[i] <- if (!is.null(res_i$cv)) res_i$cv else NA_real_
  }

  list(n = n_vec, cv_target = cv_vec)
}

#' Solve n_multi independently per domain, then aggregate
#' @keywords internal
#' @noRd
.n_multi_domains <- function(
  indicators,
  domain_cols,
  min_n_domain = NULL,
  mode = "moe",
  prop_method = "wald",
  domain_sampling = "separate"
) {
  key <- .domain_key(indicators, domain_cols)
  domain_levels <- unique(key)
  domain_keys <- factor(key, levels = domain_levels)
  split_idx <- split(seq_len(nrow(indicators)), domain_keys)

  results <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- indicators[rows, , drop = FALSE]
    sub_inner <- .label_by_row(sub[, !names(sub) %in% domain_cols,
                                   drop = FALSE], rows)
    .n_multi_simple(sub_inner)
  })
  names(results) <- domain_levels

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
  shares <- vapply(
    domain_levels,
    function(lev) {
      v <- indicators$share[split_idx[[lev]]]
      if (anyNA(v) || any(!is.finite(v)) || any(v <= 0) || any(v > 1)) {
        stop("'share' must be in (0, 1] for every indicator row", call. = FALSE)
      }
      if (diff(range(v)) > 1e-8) {
        stop(
          sprintf(
            "'share' identifies a domain, so it must be constant within one; it varies within %s",
            lev
          ),
          call. = FALSE
        )
      }
      v[1L]
    },
    numeric(1L)
  )

  total <- sum(shares)
  if (total > 1 + 1e-6) {
    stop(
      sprintf(
        "domain 'share' values sum to %.4f; they are shares of one population and cannot exceed 1",
        total
      ),
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
    for (d in which(floored)) {
      .check_multi_attainable(
        domains$.n[d],
        indicators[split_idx[[domain_levels[d]]], , drop = FALSE],
        sprintf(
          " in domain %s",
          sQuote(paste(unlist(domains[d, domain_cols]), collapse = "/"))
        )
      )
    }
  }

  n_domain_max <- max(domains$.n)
  binding_domain <- which.max(domains$.n)

  if (identical(domain_sampling, "natural")) {
    domains$.share <- .domain_shares(indicators, domain_levels, split_idx)
    yields <- domains$.n / domains$.share
    n_overall <- max(yields)
    binding_domain <- which.max(yields)
    for (d in seq_along(domain_levels)) {
      .check_multi_attainable(
        n_overall * domains$.share[d],
        indicators[split_idx[[domain_levels[d]]], , drop = FALSE],
        sprintf(
          " expected in domain %s",
          sQuote(paste(unlist(domains[d, domain_cols]), collapse = "/"))
        )
      )
    }
  } else {
    n_overall <- sum(domains$.n)
  }
  binding_label <- domains$.binding[binding_domain]
  binding_cause <- if (identical(binding_label, "(min_n_domain)")) {
    "min_n_domain"
  } else {
    results[[domain_levels[binding_domain]]]$binding_cause
  }

  res <- .new_svyplan_multi(
    n = n_overall,
    params = list(
      domain_cols = domain_cols,
      mode = mode,
      prop_method = prop_method,
      domain_sampling = domain_sampling
    ),
    indicators = indicators,
    detail = .stack_domain_detail(results, indicators, domain_cols,
                                  domain_levels, split_idx),
    binding = binding_label,
    domains = domains
  )
  res$binding_cause <- binding_cause
  res$n_domain_max <- n_domain_max
  if (!is.null(min_n_domain)) {
    res$params$min_n_domain <- min_n_domain
  }
  res
}

#' Per-row requirements of every domain, stacked
#'
#' Each domain's own solve, before `min_n_domain` or natural incidence
#' raises its size. What the rows achieve at the fielded size is a
#' precision result and is not stored here.
#' @keywords internal
#' @noRd
.stack_domain_detail <- function(results, indicators, domain_cols,
                                 domain_levels, split_idx) {
  keep <- c("name", ".n", ".cv_target", ".binding", ".constraint")
  rows <- lapply(domain_levels, function(lev) {
    d <- results[[lev]]$detail
    d <- d[, intersect(keep, names(d)), drop = FALSE]
    cbind(indicators[split_idx[[lev]], domain_cols, drop = FALSE], d)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Label unnamed indicator rows by their row in the input
#'
#' A domain is solved on its own rows, which would otherwise be numbered
#' from 1 within the domain, so two domains would both report a row "1" and
#' neither label would say which input row it is.
#' @keywords internal
#' @noRd
.label_by_row <- function(sub, rows) {
  if (!"name" %in% names(sub)) sub$name <- rows
  sub
}

#' Least root of a falling function inside a bracket
#'
#' `f(lo) > 0 >= f(hi)`. The Illinois form of regula falsi, with a
#' bisection whenever three steps leave the bracket more than half as wide,
#' which also covers a flat or jumping `f`. Stops when the bracket is
#' narrower than `tol`, or when `f(hi)` is within `ftol` below zero, and
#' returns the side where `f <= 0`. `NULL` when `max_iter` steps do not.
#' @keywords internal
#' @noRd
.falling_root <- function(f, lo, hi, f_lo, f_hi, tol, ftol = 0,
                          max_iter = 200L) {
  side <- 0L
  width <- hi - lo
  for (i in seq_len(max_iter)) {
    if (hi - lo <= tol || -f_hi <= ftol) {
      return(list(x = hi, f = f_hi))
    }
    halve <- i %% 3L == 0L && hi - lo > width / 2
    if (i %% 3L == 0L) width <- hi - lo
    x <- if (!halve && is.finite(f_lo) && f_lo > f_hi) {
      lo + f_lo * (hi - lo) / (f_lo - f_hi)
    } else {
      NA_real_
    }
    if (is.na(x) || x <= lo || x >= hi) x <- (lo + hi) / 2
    fx <- f(x)
    if (fx > 0) {
      lo <- x
      f_lo <- fx
      if (side < 0L) f_hi <- f_hi / 2
      side <- -1L
    } else {
      hi <- x
      f_hi <- fx
      if (side > 0L) f_lo <- f_lo / 2
      side <- 1L
    }
  }
  NULL
}

#' Split a budget to minimize the worst of several falling ratios
#'
#' `ratio(d, b)` is domain `d`'s worst achieved-to-target ratio at budget
#' `b`, never rising with `b`, and `Inf` where `b` reaches no design.
#' Every domain gets at least `lower[d]`. A worst ratio `t` is reachable
#' when the least budgets `h_d(t)` at which each domain reaches it fit in
#' `total`, so the optimum is the least such `t`, found by a bracketed root
#' search on `sum(h_d(t)) - total` with each `h_d(t)` found the same way
#' from the evaluations already made. No root of `ratio(d, b) = t` is
#' assumed unique, since a ratio is flat once a domain takes its whole
#' population.
#'
#' With `scaling`, each ratio is taken to fall as `1 / sqrt(b)`, which holds
#' without a stage correction, a fixed PSU count or a size floor. The split
#' is then proportional to `b r_d(b)^2` above the lower budgets, and it is
#' kept only when the ratios it gives agree.
#'
#' Money left once every domain reaches `t` goes to the domains above their
#' lower budget, since more never raises a ratio.
#' @keywords internal
#' @noRd
.joint_split <- function(ratio, lower, total, labels, scaling = FALSE) {
  nd <- length(lower)
  if (total - sum(lower) <= 1e-12 * total) {
    return(lower * total / sum(lower))
  }
  upper <- total - (sum(lower) - lower)
  seen <- lapply(seq_len(nd), function(d) list(b = numeric(0), r = numeric(0)))
  r_at <- function(d, b) {
    k <- match(b, seen[[d]]$b)
    if (!is.na(k)) {
      return(seen[[d]]$r[k])
    }
    v <- ratio(d, b)
    seen[[d]]$b <<- c(seen[[d]]$b, b)
    seen[[d]]$r <<- c(seen[[d]]$r, v)
    v
  }
  agree <- function(b, tol) {
    r <- vapply(seq_len(nd), function(d) r_at(d, b[d]), numeric(1L))
    free <- b > lower * (1 + 1e-12)
    if (!any(free) || !all(is.finite(r))) {
      return(FALSE)
    }
    t <- max(r[free])
    all(r[free] >= t * (1 - tol)) && all(r[!free] <= t * (1 + tol))
  }

  ## The split proportional to b r_d(b)^2 is exact when each ratio falls as
  ## 1 / sqrt(b). Otherwise its worst ratio still bounds the optimum above.
  propose <- NULL
  K <- vapply(seq_len(nd), function(d) upper[d] * r_at(d, upper[d])^2,
              numeric(1L))
  if (all(is.finite(K)) && all(K > 0)) {
    held <- rep(FALSE, nd)
    repeat {
      b <- lower
      rest <- total - sum(lower[held])
      b[!held] <- rest * K[!held] / sum(K[!held])
      low <- !held & b < lower
      if (!any(low)) break
      held <- held | low
    }
    b[held] <- lower[held]
    if (scaling && agree(b, 1e-9)) {
      return(b)
    }
    propose <- b
  }

  ## The least budget in [lower, upper] reaching t, NA when upper does not.
  need <- function(d, t) {
    if (r_at(d, lower[d]) <= t) {
      return(lower[d])
    }
    if (r_at(d, upper[d]) > t) {
      return(NA_real_)
    }
    b <- seen[[d]]$b
    r <- seen[[d]]$r
    lo <- max(b[r > t])
    hi <- min(b[r <= t])
    hit <- .falling_root(function(x) r_at(d, x) - t, lo, hi,
                         r_at(d, lo) - t, r_at(d, hi) - t,
                         tol = 1e-10 * upper[d], ftol = 1e-12 * t)
    if (is.null(hit)) {
      stop(sprintf("the joint split found no budget for domain '%s' at a worst ratio of %.6g",
                   labels[d], t), call. = FALSE)
    }
    hit$x
  }
  over <- function(t) {
    h <- vapply(seq_len(nd), function(d) need(d, t), numeric(1L))
    if (anyNA(h)) Inf else sum(h) - total
  }

  start <- lower + (total - sum(lower)) / nd
  t_hi <- max(vapply(seq_len(nd), function(d) r_at(d, start[d]), numeric(1L)))
  if (!is.null(propose)) {
    t_hi <- min(t_hi, max(vapply(seq_len(nd), function(d) r_at(d, propose[d]),
                                 numeric(1L))))
  }
  t_lo <- max(vapply(seq_len(nd), function(d) r_at(d, upper[d]), numeric(1L)))
  if (!is.finite(t_hi)) {
    stop("the joint split found no design for every domain at an equal share of the budget",
         call. = FALSE)
  }
  g_lo <- over(t_lo)
  t <- if (g_lo <= 0) {
    t_lo
  } else {
    hit <- .falling_root(over, t_lo, t_hi, g_lo, over(t_hi),
                         tol = 1e-9 * t_hi, ftol = 1e-10 * total)
    if (is.null(hit)) {
      stop(sprintf("the joint split did not settle between worst ratios %.6g and %.6g",
                   t_lo, t_hi), call. = FALSE)
    }
    hit$x
  }
  b <- vapply(seq_len(nd), function(d) need(d, t), numeric(1L))
  left <- total - sum(b)
  if (left > 0) {
    up <- b > lower * (1 + 1e-12)
    if (!any(up)) up <- rep(TRUE, nd)
    b[up] <- b[up] + left * b[up] / sum(b[up])
  }
  b
}

#' @rdname n_multi
#' @export
n_multi.svyplan_prec <- function(indicators, ...) {
  x <- indicators
  dots <- list(...)
  if (
    identical(x$type, "cluster") ||
      identical(x$params$design, "cluster") ||
      !is.null(x$params$stage_cost)
  ) {
    stop(
      "cluster precision must be passed to n_cluster()",
      call. = FALSE
    )
  }
  if (x$type != "multi") {
    stop("n_multi requires a svyplan_prec of type 'multi'", call. = FALSE)
  }
  tgt <- x$params$indicators
  if ("prop_method" %in% names(dots)) {
    tgt$prop_method <- NA_character_
  }
  if ("resp_rate" %in% names(dots)) {
    tgt$resp_rate <- NA_real_
  }
  if ("eligible_per_unit" %in% names(dots) && "eligible_per_unit" %in% names(tgt)) {
    tgt$eligible_per_unit <- NA_real_
  }
  if ("eligibility_rate" %in% names(dots) && "eligibility_rate" %in% names(tgt)) {
    tgt$eligibility_rate <- NA_real_
  }
  tgt$n <- NULL
  tgt$n_per_psu <- NULL
  tgt$n_per_ssu <- NULL

  stored_mode <- x$params$mode
  ## A non-Wald interval takes a margin target only, so its achieved margin
  ## is the target whatever mode the first fit used.
  method <- c(dots[["prop_method"]] %||% x$params$prop_method,
              if (!"prop_method" %in% names(dots)) tgt$prop_method)
  non_wald <- any(!is.na(method) & method != "wald")
  if (non_wald && !is.null(x$detail) && ".moe" %in% names(x$detail)) {
    stored_mode <- "moe"
  }
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

#' Evaluate a field size that rounding brought up to the census
#'
#' The census is then the rounded size, not a request beyond the frame, so
#' its notice is not raised again.
#' @keywords internal
#' @noRd
.muffle_census <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("FPC set to 0 (census)", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}
