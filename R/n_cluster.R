#' Optimal multistage cluster allocation
#'
#' Compute optimal per-stage sample sizes for a multistage cluster
#' design, minimizing cost for a given precision or minimizing
#' variance for a given budget.
#'
#' @param stage_cost For the default method: numeric vector of per-stage
#'   costs. Length determines the number of stages (2 or 3). Named vectors
#'   are accepted with stage names `cost_psu`, `cost_ssu`, `cost_tsu`
#'   (`cost_tsu` aliases `cost_ssu` in 2-stage).
#'   For `svyplan_prec` objects: a precision result from [prec_cluster()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param icc Numeric vector of homogeneity measures (length = stages - 1),
#'   or a `svyplan_varcomp` object. The ICC quantifies how similar units
#'   within the same cluster are: 0 means no similarity (clusters are as
#'   variable as the whole population), 1 means perfect similarity (all
#'   units in a cluster are identical). Typical values in household
#'   surveys range from 0.01 to 0.10. Higher icc means more clusters
#'   are needed for the same precision. Use [varcomp()] to estimate
#'   icc from a previous survey or pilot data.
#' @param unit_relvar Unit relvariance (default 1). For most applications,
#'   the default of 1 is appropriate. Non-unit values arise when working
#'   with variance components from [varcomp()] that separate the total
#'   variance into stage-specific pieces.
#' @param var_ratio Ratio of the stage components' unit variance to the analysis
#'   variable's, default 1. A scalar names `var_ratio_psu` and, for a three-stage
#'   design, derives `var_ratio_ssu = var_ratio_psu * (1 - icc_psu)`, the identity the
#'   variance decomposition imposes. Supply a length-2 vector only to
#'   override it, which is meaningful when the two stages' ratios come from
#'   different decompositions. See [design_effect()] for the identity and
#'   why the design effect would not reduce to `var_ratio_psu` without it. The
#'   default of 1 is appropriate for most designs.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.05` means the standard error of the estimate
#'   should be at most 5 percent of the estimate itself. Specify exactly
#'   one of `cv` or `budget`.
#' @param budget Total budget. Specify exactly one of `cv` or `budget`.
#' @param n_psu Fixed number of PSUs (stage-1 sample size). `NULL` (default)
#'   means optimize. For 2-stage, at most one of `n_psu` or `n_per_psu` may be
#'   specified. For 3-stage, up to two of `n_psu`, `n_per_psu`, `n_per_ssu` may
#'   be fixed.
#' @param n_per_psu Fixed cluster size (stage-2 sample size per PSU). `NULL`
#'   (default) means optimize. This is the typical MICS/DHS parameterization
#'   where the number of households per cluster is fixed.
#' @param n_per_ssu Fixed SSU take size (stage-3 sample size per SSU). `NULL`
#'   (default) means optimize. Only valid for 3-stage designs.
#' @param resp_rate_psu Expected **PSU-level** response rate, in (0, 1\].
#'   Default 1 (no adjustment). It is the share of selected clusters that can
#'   be worked at all, and the stage-1 sample size is inflated by
#'   `1 / resp_rate_psu` to cover the loss. Nonresponse among the ultimate
#'   units inside a cluster is a different quantity that acts at a different
#'   stage, and the plain `resp_rate` carries that meaning elsewhere in the
#'   package.
#' @param resp_rate_ssu Expected **SSU-level** response rate, in (0, 1\].
#'   Three-stage designs only, default 1. It is the share of selected
#'   second-stage units, typically households, that can be interviewed at
#'   all. A 2-stage design has no such stage, since the units inside a PSU
#'   are the ultimate ones, and supplying it there is an error.
#' @param resp_rate Expected **ultimate-unit** response rate, in (0, 1\].
#'   Default 1. It is the share of selected final-stage units, typically
#'   persons, that respond. This is the plain `resp_rate` the rest of the
#'   package uses. See the nonresponse section of [svyplan-package] for what
#'   this adjustment does and does not claim.
#' @param fixed_cost Fixed overhead cost (C0). Default 0.
#'   The total cost model becomes
#'   `C = C0 + c1*n_psu + c2*n_psu*n_per_psu [+ c3*n_psu*n_per_psu*n_per_ssu]`.
#'   In budget mode, only `budget - fixed_cost` is available for variable
#'   costs. In CV mode, `fixed_cost` is added to the variable cost.
#' @param indicators Optional data frame with one row per indicator, which
#'   switches the function to the several-indicators mode of Details. Each
#'   row carries its own homogeneity and its own precision target, and the
#'   allocation returned satisfies every row. Leave it `NULL` for the single
#'   indicator that `icc`, `unit_relvar`, `var_ratio` and `cv` describe.
#' @param domains Optional character vector naming domain columns in
#'   `indicators`. Each domain is solved on its own unless
#'   `allocation = "joint"` in budget mode. Domains are sized as separate
#'   quotas, so `$total_n` is their sum and `$n` carries no aggregate stage
#'   vector. The fieldable per-domain stage sizes are in `$domains`, and
#'   `domain_sampling` covers why the quota reading is the only one offered
#'   for a multistage design.
#' @param allocation How a budget is split across domains, either
#'   `"separate"` (default, each domain sized on its own) or `"joint"` (one
#'   budget split across domains to minimize the worst precision ratio).
#'   `"joint"` applies only when `indicators`, `domains` and `budget` are all
#'   supplied.
#' @param domain_sampling How the per-domain requirements combine into the
#'   one overall size. Only `"separate"` (the default) is available for a
#'   cluster design, so the domains are quotas fielded in their own right and
#'   `$total_n` is the sum of the per-domain totals. `"natural"`, which
#'   [n_multi()] offers for a single-stage design, is refused, because a
#'   domain's expected yield then depends on how its members sit inside PSUs
#'   and SSUs rather than on its share of the population alone, and that
#'   model is not in the package.
#' @param min_n_domain Optional positive minimum total sample size per
#'   domain. In joint budget mode it is a constraint. With domains solved
#'   separately, domains below the floor produce a warning.
#' @param plan Optional [svyplan()] object providing design defaults
#'   (including `stage_cost`, `icc`, `unit_relvar`, `var_ratio`, `resp_rate_psu`, `fixed_cost`).
#'   In the several-indicators mode a profile supplies only the defaults that
#'   mode accepts, so `icc`, `unit_relvar` and `var_ratio` are left to the
#'   indicator columns.
#'
#' @return A `svyplan_cluster` object with components:
#' \describe{
#'   \item{`n`}{Named numeric vector of continuous per-stage sample sizes
#'     (e.g. `c(n_psu = 84.1, n_per_psu = 13.8)`), the mathematical
#'     optimum.}
#'   \item{`stages`}{Number of stages (2 or 3).}
#'   \item{`total_n`}{Continuous total sample size (`prod(n)`).}
#'   \item{`cv`}{Coefficient of variation of the continuous optimum.}
#'   \item{`cost`}{Cost of the continuous optimum.}
#'   \item{`operational`}{The whole-unit field design, found by a
#'     discrete search: `n` (named integer stage sizes), `total_n`,
#'     `cost`, and `cv`, all recomputed from the integer design. In
#'     `budget` mode its cost never exceeds the budget, whereas in `cv` mode it
#'     meets the target at the lowest cost among the designs searched.
#'     `as.integer()` returns `operational$n`, and `as.double()` returns the
#'     continuous `n`.}
#'   \item{`params`}{List of input parameters.}
#' }
#'   With `indicators`, the same class carries `$indicators`, the per-row
#'   `$detail` and, where domains were named, `$domains`. The class does not
#'   depend on which optional arguments were supplied.
#'
#' @details
#' ## Getting started
#'
#' A typical 2-stage household survey workflow:
#'
#' 1. **Decide what you know.** You need the cost per cluster visit
#'    (`stage_cost[1]`, e.g. travel + logistics) and the cost per
#'    interview (`stage_cost[2]`), plus an estimate of within-cluster
#'    homogeneity (`icc`). Estimate icc from a pilot or previous
#'    survey with [varcomp()], or use a plausible range (0.01--0.10
#'    for most household indicators).
#' 2. **Choose a mode.** If you have a target precision, set `cv`. If you
#'    have a fixed budget, set `budget`. Never set both.
#' 3. **Fix stages or let the optimizer decide.** In MICS/DHS-style
#'    designs, the number of households per cluster is fixed by
#'    fieldwork logistics (e.g. `n_per_psu = 20`). The optimizer then
#'    solves for how many clusters (`n_psu`) to visit. If no stage is
#'    fixed, both are optimized jointly.
#'
#' ## How it works
#'
#' Stage count is determined by `length(stage_cost)`:
#' - **2-stage** (e.g. clusters then households): `stage_cost` has 2
#'   elements, `icc` is a scalar.
#' - **3-stage** (e.g. districts, clusters, households): `stage_cost`
#'   has 3 elements, `icc` is length 2.
#'
#' Two solving modes:
#' - **CV mode**: minimize total cost subject to achieving the target CV.
#' - **Budget mode**: minimize the CV (maximize precision) within the
#'   available budget.
#'
#' ## Fixing stage sizes
#'
#' One or more stage sizes can be fixed, leaving the remaining stage(s) to be
#' optimized or derived from the constraint. For 2-stage designs, at most one
#' stage may be fixed. For 3-stage designs, up to two stages may be fixed.
#' The remaining free stage is derived from the budget or CV constraint.
#'
#' If `icc` is a `svyplan_varcomp` object, `icc`, `unit_relvar`, and `var_ratio`
#' are extracted automatically.
#'
#' ## Several indicators at once
#'
#' `indicators` sizes one design against several requirements at once. It is a
#' data frame with one row per indicator, and the design returned satisfies
#' every row. The scalar arguments `icc`, `unit_relvar`, `var_ratio` and `cv`
#' are refused in this mode, because the frame carries each of them per row.
#'
#' Each row needs `p`, `var`, or the ratio quartet (`r`, `cv_num`, `cv_den`,
#' `component_cor`), a `cv`, `moe` or `rmoe` target, and `icc_psu`. A
#' three-stage design also needs `icc_ssu`. Optional `var_ratio_psu` defaults
#' to 1, and three-stage `var_ratio_ssu` is derived as
#' `var_ratio_psu * (1 - icc_psu)` when absent, the value the variance
#' decomposition implies (see [design_effect()]). On a ratio row `icc_psu`,
#' `icc_ssu` and the stage variance ratios describe the linearized variable
#' `e = y - r * x`, not either component. The remaining columns follow
#' [n_multi()].
#'
#' Nonresponse is still named for the stage it acts on. `resp_rate_psu`,
#' `resp_rate_ssu` and `resp_rate` are the defaults used where a row's own
#' column is absent or `NA`. A column naming a stage the design does not have
#' is an error rather than a column carried along and ignored, since a
#' silently dropped response rate plans a design with none.
#'
#' Margin-of-error rows are converted to CV before optimization, respecting
#' the row's `prop_method`. The multistage model is driven by a relative
#' standard error, and only the Wald interval has half-width `z * se`, so a
#' proportion row's `moe` is restated as the sampling CV at the effective
#' sample size its own method needs to close the interval to that margin. A
#' stricter interval therefore asks for a larger design, in the same order it
#' does in [n_prop()]. Under `"wald"` the restatement is `moe / (z * p)`
#' exactly. A mean row converts as `moe / (z * |mu|)`, on the magnitude so
#' that a negative mean yields a positive target. [prec_cluster()] inverts the
#' same way, so a design sized from a `moe` target reports that `moe` back.
#'
#' For each candidate allocation the required stage-1 size is the maximum
#' across rows. The solver minimizes total cost against the precision targets,
#' or the worst precision ratio under a fixed budget. `domains` sizes each
#' domain on its own, `allocation = "joint"` splits one budget across them to
#' minimize the worst precision ratio, and `min_n_domain` sets a floor.
#'
#' ## How strong the several-indicators optimum is
#'
#' With a stage size fixed, the remaining problem is solved from the
#' closed-form cluster optimum. With all stage sizes free, the objective is a
#' maximum over indicator requirements, which is not smooth, and it is
#' minimized by a bounded quasi-Newton search. That search warns when it fails
#' to converge or lands on a bound, but it carries no global-optimality or KKT
#' certificate, so a successful return means the best design this search
#' found, not a proven minimum-cost one. The `$operational` allocation can
#' always be checked against its own precision or budget constraint, which is
#' a separate and exact statement. [n_alloc()] gives the stronger guarantee
#' where it applies, returning feasibility and KKT diagnostics for its convex
#' continuous problem.
#'
#' ## Boundary icc values
#'
#' Boundary and near-boundary homogeneity values are not supported by the
#' analytical optimum used here. When `icc` is near 0, most variability is
#' within PSUs, so the closed-form optimum collapses toward taking many units
#' in very few PSUs. When `icc` is near 1, most variability is between PSUs,
#' so the optimum collapses toward taking very few units in many PSUs. In both
#' cases the analytical allocation becomes degenerate, so `n_cluster()`
#' rejects values numerically too close to 0 or 1.
#'
#' ## Nonresponse acts at whichever stage it happens
#'
#' A cluster design can lose units at more than one stage, and the losses are
#' not interchangeable, so each has its own argument named for the stage it
#' acts on. Writing the three-stage CV out with realized takes shows why:
#'
#' \deqn{cv^2 = \frac{V k_1 \delta_1}{a\,r_{psu}}
#'   + \frac{V k_2 (1 + \delta_2 (q r - 1))}{a\,r_{psu}\, m\,r_{ssu}\, q\,r}}{cv^2 = (V k_1 delta_1)/(a r_psu) + (V k_2 (1 + delta_2 (q r - 1)))/(a r_psu m r_ssu q r)}
#'
#' `resp_rate_psu` divides both terms, so losing whole PSUs is a pure
#' sample-size loss. `resp_rate_ssu` divides only the second. `resp_rate`
#' divides the second *and* enters the \eqn{\delta_2} bracket, because it
#' shrinks the realized final-stage take and so changes the clustering
#' penalty itself. None is recoverable from the others.
#'
#' It follows that the later-stage rates move the cost-optimal design while
#' the PSU rate does not. For two stages the optimal take becomes
#' \deqn{b^* = \sqrt{\frac{C_1 (1 - \mathrm{icc})}{C_2\,\mathrm{icc}\,r}},}{b^* = sqrt((C_1 (1 - icc))/(C_2 icc r)),}
#' with `r` the ultimate-unit rate. `resp_rate_psu` is absent because it
#' scales cost without moving that trade-off. Sizes and costs stay **gross**:
#' `$n` counts units to issue and `$cost` pays for them, while the variance
#' reads what they realize.
#'
#' This is a deterministic expected-take approximation. It substitutes the
#' expected realized stage counts into the planning variance and assumes
#' response is noninformative for it. Below one expected respondent per unit
#' the substitution stops describing the design, and that is an error naming
#' the approximation rather than the design.
#'
#' These functions assume sampling fractions are negligible at each stage
#' (equivalent to sampling with replacement). No finite population correction
#' is applied. This is standard for multistage planning when cluster
#' populations are large relative to the sample.
#'
#' [n_alloc()] differs here. Its cluster mode carries the ultimate-unit
#' correction `1 - n / N`, so the two agree only where the ultimate-unit
#' sampling fraction is negligible. They are the same variance model
#' otherwise, and the clustering bracket, the stage response rates and the
#' cost-optimal take are shared exactly. Size a design through both at an
#' appreciable sampling fraction and `n_alloc()` returns the smaller answer,
#' by that correction and nothing else.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 9.
#'
#' @family cluster design functions
#' @seealso [prec_cluster()] for the inverse, [varcomp()] for estimating
#'   variance components, [n_multi()] for several indicators without
#'   clustering, and [n_alloc()] for a stratified multistage allocation.
#'
#' @examples
#' # 2-stage, budget mode
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
#'
#' # 2-stage, CV mode
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
#'
#' # 2-stage, fixed n_psu
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000, n_psu = 40)
#'
#' # 2-stage, fixed n_per_psu (MICS/DHS style: 20 households per cluster)
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000, n_per_psu = 20)
#'
#' # 3-stage
#' n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05), cv = 0.05)
#'
#' # 3-stage, fixed n_psu + n_per_ssu (solve for n_per_psu)
#' n_cluster(
#'   stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
#'   budget = 500000, n_psu = 50, n_per_ssu = 8
#' )
#'
#' # With fixed overhead cost
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000, fixed_cost = 5000)
#'
#' # Several indicators, CV mode
#' indicators <- data.frame(
#'   name    = c("stunting", "anemia"),
#'   p       = c(0.30, 0.10),
#'   cv      = c(0.10, 0.15),
#'   icc_psu = c(0.02, 0.05)
#' )
#' n_cluster(indicators = indicators, stage_cost = c(500, 50))
#'
#' # Several indicators with MOE targets, converted to CV internally
#' by_moe <- data.frame(
#'   name    = c("stunting", "anemia"),
#'   p       = c(0.30, 0.10),
#'   moe     = c(0.05, 0.03),
#'   icc_psu = c(0.02, 0.05)
#' )
#' n_cluster(indicators = by_moe, stage_cost = c(500, 50))
#'
#' # One budget split across domains
#' by_domain <- data.frame(
#'   name    = rep(c("stunting", "anemia"), each = 2),
#'   p       = c(0.30, 0.25, 0.10, 0.15),
#'   cv      = c(0.10, 0.10, 0.15, 0.15),
#'   icc_psu = c(0.02, 0.03, 0.05, 0.04),
#'   region  = rep(c("Urban", "Rural"), 2)
#' )
#' n_cluster(
#'   indicators = by_domain,
#'   stage_cost = c(500, 50),
#'   domains = "region",
#'   budget = 100000,
#'   allocation = "joint"
#' )
#'
#' @export
n_cluster <- function(stage_cost = NULL, ...) {
  if (!missing(stage_cost)) {
    .res <- .dispatch_plan(stage_cost, "stage_cost", n_cluster.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_cluster")
}

#' @rdname n_cluster
#' @export
n_cluster.default <- function(
  stage_cost = NULL,
  ...,
  icc = NULL,
  unit_relvar = NULL,
  var_ratio = NULL,
  cv = NULL,
  budget = NULL,
  n_psu = NULL,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  resp_rate_psu = 1,
  resp_rate_ssu = 1,
  resp_rate = 1,
  fixed_cost = 0,
  indicators = NULL,
  domains = NULL,
  allocation = c("separate", "joint"),
  domain_sampling = c("separate", "natural"),
  min_n_domain = NULL,
  plan = NULL
) {
  plan <- .plan_for_cluster_mode(plan, !is.null(indicators))
  # Before the merge, which is the only point at which the profile's own
  # defaults are still distinguishable from the formals.
  .check_relvar_identified(
    icc,
    !is.null(unit_relvar) || "unit_relvar" %in% names(plan$defaults),
    "n_cluster()"
  )
  .plan <- .merge_plan_args(plan, n_cluster.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_cluster.default, c(.plan, list(...))))
  .check_unused_dots(...)
  .check_cluster_mode_args(
    indicators,
    icc = icc, unit_relvar = unit_relvar, var_ratio = var_ratio, cv = cv,
    domains = domains, allocation = allocation,
    domain_sampling = domain_sampling, min_n_domain = min_n_domain
  )
  if (!is.null(indicators)) {
    return(.n_cluster_indicators(
      indicators,
      stage_cost = stage_cost,
      domains = domains,
      budget = budget,
      n_psu = n_psu,
      n_per_psu = n_per_psu,
      n_per_ssu = n_per_ssu,
      allocation = allocation,
      domain_sampling = domain_sampling,
      min_n_domain = min_n_domain,
      fixed_cost = fixed_cost,
      resp_rate_psu = resp_rate_psu,
      resp_rate_ssu = resp_rate_ssu,
      resp_rate = resp_rate
    ))
  }
  if (is.data.frame(stage_cost)) {
    .stop_frame_in_first_slot("stage_cost", "per-stage costs")
  }
  unit_relvar <- unit_relvar %||% 1
  var_ratio <- var_ratio %||% 1
  if (is.null(stage_cost))
    stop("'stage_cost' is required (directly or via plan)", call. = FALSE)
  if (is.null(icc))
    stop("'icc' is required (directly or via plan)", call. = FALSE)
  check_stage_cost(stage_cost)
  stage_cost <- .reorder_stage_cost(stage_cost)
  check_resp_rate(resp_rate_psu, "resp_rate_psu")
  check_resp_rate(resp_rate_ssu, "resp_rate_ssu")
  check_resp_rate(resp_rate, "resp_rate")

  if (inherits(icc, "svyplan_varcomp")) {
    vc <- icc
    if (!is.null(vc$strata)) {
      stop(
        "stratified varcomp: merge its $strata columns into an n_alloc() frame, or pass one stratum's values",
        call. = FALSE
      )
    }
    icc <- vc$icc
    var_ratio <- vc$var_ratio
    if (!is.na(vc$unit_relvar)) unit_relvar <- vc$unit_relvar
  }

  check_scalar(unit_relvar, "unit_relvar")
  if (
    !is.numeric(var_ratio) ||
      length(var_ratio) == 0L ||
      anyNA(var_ratio) ||
      any(var_ratio <= 0) ||
      any(!is.finite(var_ratio))
  ) {
    stop("'var_ratio' must contain positive finite values", call. = FALSE)
  }

  stages <- length(stage_cost)
  icc <- .reorder_stage_vec(icc, "icc")
  var_ratio <- .reorder_stage_vec(var_ratio, "var_ratio")
  check_icc(icc, expected_length = stages - 1L)
  .check_cluster_icc_open(icc, context = "n_cluster()")

  has_cv <- !is.null(cv)
  has_budget <- !is.null(budget)
  if (has_cv == has_budget) {
    stop("specify exactly one of 'cv' or 'budget'", call. = FALSE)
  }
  if (has_cv) {
    check_scalar(cv, "cv")
  }
  if (has_budget) {
    check_scalar(budget, "budget")
  }
  if (!is.null(n_psu)) check_scalar(n_psu, "n_psu")
  if (!is.null(n_per_psu)) check_scalar(n_per_psu, "n_per_psu")
  if (!is.null(n_per_ssu)) check_scalar(n_per_ssu, "n_per_ssu")
  for (nm in c("n_psu", "n_per_psu", "n_per_ssu")) {
    v <- get(nm)
    if (!is.null(v) && v < 1) {
      stop(sprintf("'%s' must be at least 1", nm), call. = FALSE)
    }
  }
  if (stages == 2L && !is.null(n_per_ssu)) {
    stop("'n_per_ssu' is not applicable for 2-stage designs", call. = FALSE)
  }
  if (stages == 2L && !isTRUE(all.equal(resp_rate_ssu, 1))) {
    stop(
      "'resp_rate_ssu' is not applicable for 2-stage designs: the units inside a PSU are the ultimate ones, so their nonresponse is 'resp_rate'",
      call. = FALSE
    )
  }
  n_fixed <- sum(!is.null(n_psu), !is.null(n_per_psu), !is.null(n_per_ssu))
  if (n_fixed >= stages) {
    stop("cannot fix all stages; use prec_cluster() instead", call. = FALSE)
  }
  check_fixed_cost(fixed_cost, budget)

  res <- if (stages == 2L) {
    .n_cluster_2stage(
      stage_cost,
      icc,
      unit_relvar,
      var_ratio,
      cv,
      budget,
      n_psu,
      n_per_psu,
      resp_rate_psu,
      resp_rate,
      fixed_cost
    )
  } else {
    var_ratio <- .stage_k_pair(var_ratio, icc)
    .n_cluster_3stage(
      stage_cost,
      icc,
      unit_relvar,
      var_ratio,
      cv,
      budget,
      n_psu,
      n_per_psu,
      n_per_ssu,
      resp_rate_psu,
      resp_rate_ssu,
      resp_rate,
      fixed_cost
    )
  }
  res
}

#' @rdname n_cluster
#' @export
n_cluster.svyplan_prec <- function(stage_cost, ..., cv = NULL, budget = NULL) {
  x <- stage_cost
  if (
    identical(x[["type"]], "multi") &&
      identical(x[["params"]][["design"]], "cluster")
  ) {
    return(.n_cluster_from_indicator_prec(x, list(...), cv, budget))
  }
  if (x$type != "cluster") {
    stop(
      "n_cluster requires a svyplan_prec of type 'cluster', from prec_cluster()",
      call. = FALSE
    )
  }
  p <- x$params
  if (is.null(cv) && is.null(budget)) {
    cv <- x$cv
  }
  args <- list(
    stage_cost = p$stage_cost,
    icc = p$icc,
    unit_relvar = p$unit_relvar,
    var_ratio = p$var_ratio,
    cv = cv,
    budget = budget,
    n_psu = p$n_psu,
    n_per_psu = p$n_per_psu,
    n_per_ssu = p$n_per_ssu,
    resp_rate_psu = p$resp_rate_psu %||% 1,
    resp_rate_ssu = p$resp_rate_ssu %||% 1,
    resp_rate = p$resp_rate %||% 1,
    fixed_cost = p$fixed_cost %||% 0
  )
  do.call(n_cluster.default, .roundtrip_args(args, list(...), n_cluster.default))
}

#' Rebuild a several-indicators cluster allocation from its precision
#'
#' The target the design was sized against is read back from `$detail`, in
#' whichever of `cv` and `moe` the fit recorded, so the round trip returns
#' the allocation rather than a design sized against a different quantity.
#' @keywords internal
#' @noRd
.n_cluster_from_indicator_prec <- function(x, dots, cv, budget) {
  if (!is.null(cv)) {
    stop(
      "'cv' is carried by the 'indicators' column 'cv', one value per row",
      call. = FALSE
    )
  }
  p <- x[["params"]]
  tgt <- p[["indicators"]]
  for (rate in intersect(
    c("resp_rate_psu", "resp_rate_ssu", "resp_rate"),
    names(dots)
  )) {
    tgt[[rate]] <- NA_real_
  }
  tgt$n <- NULL
  tgt$n_per_psu <- NULL
  tgt$n_per_ssu <- NULL

  stored_mode <- p[["mode"]]
  if (identical(stored_mode, "moe")) {
    tgt$moe <- x$detail$.moe
    tgt$cv <- NULL
  } else if (!is.null(x$detail) && ".cv" %in% names(x$detail)) {
    tgt$cv <- x$detail$.cv
    tgt$moe <- NULL
  }

  args <- list(
    indicators = tgt,
    stage_cost = p[["stage_cost"]],
    domains = p[["domain_cols"]],
    budget = budget %||% p[["budget"]],
    n_psu = p[["n_psu"]],
    n_per_psu = p[["n_per_psu"]],
    n_per_ssu = p[["n_per_ssu"]],
    allocation = p[["allocation"]] %||%
      if (isTRUE(p[["joint"]])) "joint" else "separate",
    min_n_domain = p[["min_n_domain"]],
    domain_sampling = p[["domain_sampling"]] %||% "separate",
    fixed_cost = p[["fixed_cost"]] %||% 0
  )
  do.call(
    n_cluster.default,
    .roundtrip_args(args, dots, n_cluster.default)
  )
}

#' Discrete (integer) 2-stage design
#'
#' Enumerates whole n_per_psu values (or uses rounded fixed sizes),
#' computes the largest affordable whole n_psu in budget mode or the
#' smallest sufficient whole n_psu in cv mode, and returns the best
#' integer design with its own recomputed cost and cv. Budget-mode
#' designs never exceed the budget. CV-mode designs meet the target.
#' @keywords internal
#' @noRd
.op_cluster_2stage <- function(stage_cost, icc, unit_relvar, var_ratio, cv, budget,
                               n_psu, n_per_psu, resp_rate_psu, resp_rate,
                               fixed_cost, cont_m) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  # m is enumerated gross, because that is what a field team can be sent to
  # do; the variance reads the realized take m * resp_rate, and the cost the
  # gross one.
  cv_fn <- function(a, m) {
    mr <- m * resp_rate
    sqrt(unit_relvar * var_ratio * (1 + icc * (mr - 1)) / (a * resp_rate_psu * mr))
  }

  m_cand <- if (!is.null(n_per_psu)) {
    max(1L, as.integer(round(n_per_psu)))
  } else {
    seq_len(max(100L, min(10L * as.integer(ceiling(cont_m)) + 1L, 100000L)))
  }

  if (!is.null(budget)) {
    if (!is.null(n_psu)) {
      a <- max(1L, as.integer(round(n_psu)))
      m <- as.integer(floor((vb / a - C1) / C2))
      if (m < 1L) {
        stop("no whole-unit design fits the budget for the given fixed stage-1 size",
             call. = FALSE)
      }
      a_cand <- a
      m_best <- m
    } else {
      a_cand <- as.integer(floor(vb / (C1 + C2 * m_cand)))
      keep <- a_cand >= 1L
      m_cand <- m_cand[keep]
      a_cand <- a_cand[keep]
      if (length(m_cand) == 0L) {
        stop("no whole-unit design fits the budget", call. = FALSE)
      }
      cvs <- cv_fn(a_cand, m_cand)
      j <- which.min(cvs + 1e-12 * m_cand)
      m_best <- m_cand[j]
      a_cand <- a_cand[j]
    }
    a_best <- a_cand
  } else {
    if (!is.null(n_psu)) {
      a_best <- max(1L, as.integer(round(n_psu)))
      need <- unit_relvar * var_ratio * (1 - icc) /
        (cv^2 * a_best * resp_rate_psu - unit_relvar * var_ratio * icc)
      if (!is.finite(need) || need <= 0) {
        stop("target CV is not achievable with whole units at the given fixed stage-1 size",
             call. = FALSE)
      }
      # 'need' counts responding units; the take that delivers them is gross.
      m_best <- max(1L, as.integer(ceiling(need / resp_rate - 1e-9)))
    } else {
      a_cand <- vapply(m_cand, function(m) {
        mr <- m * resp_rate
        as.integer(ceiling(
          unit_relvar * var_ratio * (1 + icc * (mr - 1)) /
            (mr * cv^2 * resp_rate_psu) - 1e-9
        ))
      }, integer(1L))
      a_cand <- pmax(a_cand, 1L)
      costs <- a_cand * (C1 + C2 * m_cand)
      j <- which.min(costs + 1e-9 * (a_cand * m_cand) + 1e-12 * m_cand)
      m_best <- m_cand[j]
      a_best <- a_cand[j]
    }
  }

  list(
    n = c(n_psu = a_best, n_per_psu = m_best),
    total_n = a_best * m_best,
    cost = fixed_cost + a_best * (C1 + C2 * m_best),
    cv = cv_fn(a_best, m_best)
  )
}

#' Best whole-unit 3-stage design, searched over PSU count and SSU take
#'
#' Every three-stage planning variance is
#' `cv_j^2 = (alpha_j + beta_j / m + gamma_j / (m q)) / a` in the gross takes
#' `a` (PSUs), `m` (SSUs per PSU) and `q` (units per SSU), so the ultimate
#' take that reaches a given PSU count is available in closed form. Searching
#' `(a, m)` and solving for `q` covers the whole `q` axis, which a window on
#' the continuous optimum cannot: in cv mode the cost is a step function of
#' `q` whose minimum sits where `ceiling()` of the required PSU count is
#' tight, and that can be an order of magnitude away from the continuous take.
#'
#' The coefficients arrive already divided by the indicator's squared target
#' cv, so `max_j(alpha_j + beta_j / m + gamma_j / (m q))` is the PSU count the
#' design requires, and dividing it by `a` gives the squared ratio of achieved
#' to target cv. A single-indicator budget-mode caller, which has no target,
#' passes a scale of 1 and reads that ratio as the squared cv itself.
#'
#' `m_cand` is supplied by the caller and holds the fixed take when
#' `n_per_psu` is pinned, so the search is exact in `a` and `q` and bounded in
#' `m`.
#'
#' @return Named integer vector `c(n_psu, n_per_psu, n_per_ssu)`.
#' @keywords internal
#' @noRd
.op_cluster3_search <- function(alpha, beta, gamma, stage_cost,
                                variable_budget, n_psu, n_per_ssu, m_cand,
                                msg_budget, msg_cv) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  tol <- 1e-9
  m_cand <- as.numeric(m_cand)
  nm <- length(m_cand)
  nj <- length(alpha)
  budget_mode <- !is.null(variable_budget)
  a_fixed <- if (is.null(n_psu)) NULL else max(1, round(n_psu))
  q_fixed <- if (is.null(n_per_ssu)) NULL else max(1, round(n_per_ssu))

  need_at <- function(m, q) {
    out <- alpha[1L] + beta[1L] / m + gamma[1L] / (m * q)
    for (j in seq_len(nj)[-1L]) {
      out <- pmax(out, alpha[j] + beta[j] / m + gamma[j] / (m * q))
    }
    out
  }

  # Smallest whole ultimate take reaching a PSU count of 'a'; Inf where the
  # PSU count is out of reach at that SSU take whatever 'q' does.
  q_needed <- function(a, m) {
    out <- rep(1, length(m))
    for (j in seq_len(nj)) {
      d <- a - alpha[j] - beta[j] / m
      out <- pmax(out, ifelse(d > 0, gamma[j] / (m * d), Inf))
    }
    ceiling(out - tol)
  }

  best <- list(key = Inf)
  consider <- function(a, m, q) {
    valid <- a >= 1 & is.finite(q) & q >= 1 & q <= .Machine$integer.max
    q[!valid] <- 1
    cost <- a * (C1 + C2 * m + C3 * m * q)
    if (budget_mode) {
      valid <- valid & cost <= variable_budget + 1e-8
      key <- need_at(m, q) / a + 1e-12 * (m + q)
    } else {
      valid <- valid & need_at(m, q) <= a + tol + 1e-12 * a
      key <- cost + 1e-9 * (a * m * q) + 1e-12 * (m + q)
    }
    key[!valid] <- Inf
    i <- which.min(key)
    if (length(i) == 1L && is.finite(key[i]) && key[i] < best$key) {
      best <<- list(key = key[i], cost = cost[i], a = a[i], m = m[i], q = q[i])
    }
    invisible(NULL)
  }

  # Evaluate a block of PSU counts against every candidate SSU take at once.
  scan_a <- function(a_seq) {
    block <- max(1L, as.integer(2e5 %/% nm))
    i <- 1L
    while (i <= length(a_seq)) {
      A <- a_seq[seq.int(i, min(length(a_seq), i + block - 1L))]
      av <- rep(A, times = nm)
      mv <- rep(m_cand, each = length(A))
      qv <- if (budget_mode) {
        floor((variable_budget / av - C1 - C2 * mv) / (C3 * mv) + tol)
      } else {
        q_needed(av, mv)
      }
      consider(av, mv, qv)
      i <- i + block
    }
  }

  if (!is.null(q_fixed)) {
    q <- rep(q_fixed, nm)
    a <- if (!is.null(a_fixed)) {
      rep(a_fixed, nm)
    } else if (budget_mode) {
      floor((variable_budget + 1e-8) / (C1 + C2 * m_cand + C3 * m_cand * q))
    } else {
      pmax(1, ceiling(need_at(m_cand, q) - tol))
    }
    consider(a, m_cand, q)
  } else if (!is.null(a_fixed)) {
    a <- rep(a_fixed, nm)
    q <- if (budget_mode) {
      floor((variable_budget / a_fixed - C1 - C2 * m_cand) / (C3 * m_cand) + tol)
    } else {
      q_needed(a_fixed, m_cand)
    }
    consider(a, m_cand, q)
  } else {
    # Cheapest a PSU can be, used to bound the PSU count from the cost side.
    m_lo <- min(m_cand)
    psu_floor <- C1 + m_lo * (C2 + C3)
    if (budget_mode) {
      a_hi <- floor((variable_budget + 1e-8) / psu_floor)
      if (a_hi < 1) stop(msg_budget, call. = FALSE)
      # Thin the sweep only if the grid would be enormous, then refine.
      stride <- max(1, ceiling(a_hi * nm / 5e6))
      scan_a(seq.int(1, a_hi, by = stride))
      if (stride > 1 && is.finite(best$key)) {
        scan_a(seq.int(max(1, best$a - stride), min(a_hi, best$a + stride)))
      }
    } else {
      # Above max(alpha + beta / m_lo) some take is feasible, so this seeds a
      # finite cost and with it the cost bound the sweep prunes on.
      a_lo <- max(1, ceiling(max(alpha) - tol))
      a_seed <- max(a_lo, floor(max(alpha + beta / m_lo) + tol) + 1)
      consider(rep(a_seed, nm), m_cand, q_needed(a_seed, m_cand))
      block <- max(1L, as.integer(2e5 %/% nm))
      a_cur <- a_lo
      a_hi <- floor((best$cost + 1e-8) / psu_floor)
      while (a_cur <= a_hi) {
        scan_a(seq.int(a_cur, min(a_hi, a_cur + block - 1L)))
        a_cur <- a_cur + block
        a_hi <- min(a_hi, floor((best$cost + 1e-8) / psu_floor))
      }
    }
  }

  if (!is.finite(best$key)) {
    stop(if (budget_mode) msg_budget else msg_cv, call. = FALSE)
  }
  c(
    n_psu = as.integer(best$a),
    n_per_psu = as.integer(best$m),
    n_per_ssu = as.integer(best$q)
  )
}

#' Discrete (integer) 3-stage design
#' @keywords internal
#' @noRd
.op_cluster_3stage <- function(stage_cost, icc, unit_relvar, var_ratio, cv, budget,
                               n_psu, n_per_psu, n_per_ssu, resp_rate_psu,
                               resp_rate_ssu, resp_rate,
                               fixed_cost, cont_m) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  delta1 <- icc[1L]
  delta2 <- icc[2L]
  k1 <- var_ratio[1L]
  k2 <- var_ratio[2L]
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  # m and q are enumerated gross; the variance reads the realized takes and
  # the cost the gross ones.
  cv_fn <- function(a, m, q) {
    mr <- m * resp_rate_ssu
    qr <- q * resp_rate
    sqrt(
      unit_relvar / (a * resp_rate_psu * mr * qr) *
        (k1 * delta1 * mr * qr + k2 * (1 + delta2 * (qr - 1)))
    )
  }

  # Same variance, written in the gross takes as (alpha + beta/m + gamma/(mq))/a.
  # Budget mode has no target to divide by and reads the numerator as cv^2.
  scale <- if (is.null(cv)) 1 else cv^2
  alpha <- unit_relvar * k1 * delta1 / (resp_rate_psu * scale)
  beta <- unit_relvar * k2 * delta2 / (resp_rate_psu * resp_rate_ssu * scale)
  gamma <- unit_relvar * k2 * (1 - delta2) /
    (resp_rate_psu * resp_rate_ssu * resp_rate * scale)

  m_cand <- if (!is.null(n_per_psu)) {
    max(1L, as.integer(round(n_per_psu)))
  } else {
    seq_len(max(60L, min(6L * as.integer(ceiling(cont_m)) + 1L, 400L)))
  }

  hit <- .op_cluster3_search(
    alpha, beta, gamma, stage_cost, vb, n_psu, n_per_ssu, m_cand,
    msg_budget = "no whole-unit design fits the budget",
    msg_cv = "target CV is not achievable with whole units at the given fixed stage-1 size"
  )
  a_best <- hit[["n_psu"]]
  m_best <- hit[["n_per_psu"]]
  q_best <- hit[["n_per_ssu"]]
  list(
    n = c(n_psu = a_best, n_per_psu = m_best, n_per_ssu = q_best),
    total_n = a_best * m_best * q_best,
    cost = fixed_cost +
      a_best * (C1 + C2 * m_best + C3 * m_best * q_best),
    cv = cv_fn(a_best, m_best, q_best)
  )
}

#' @keywords internal
#' @noRd
.n_cluster_2stage <- function(
  stage_cost,
  icc,
  unit_relvar,
  var_ratio,
  cv,
  budget,
  n_psu,
  n_per_psu,
  resp_rate_psu,
  resp_rate = 1,
  fixed_cost = 0
) {
  C1 <- stage_cost[1L]
  # Solved in realized units per PSU, so the closed forms need no separate
  # derivation. Only the take is converted back at the end.
  C2 <- stage_cost[2L] / resp_rate
  n_per_psu_gross <- n_per_psu
  if (!is.null(n_per_psu)) n_per_psu <- n_per_psu * resp_rate
  var_budget <- if (!is.null(budget)) budget - fixed_cost else NULL

  if (is.null(n_psu) && is.null(n_per_psu)) {
    n2_opt <- sqrt(C1 / C2 * (1 - icc) / icc)
    if (n2_opt < 1) {
      warning(
        "cost-optimal n_per_psu is below 1 (high 'icc' relative to stage costs); clamped to 1",
        call. = FALSE
      )
      n2_opt <- 1
    }

    if (!is.null(budget)) {
      n1_opt <- var_budget / (C1 + C2 * n2_opt)
      n1_eff <- n1_opt * resp_rate_psu
      cv_achieved <- sqrt(
        unit_relvar / (n1_eff * n2_opt) * var_ratio * (1 + icc * (n2_opt - 1))
      )
      total_cost <- budget
    } else {
      n1_eff_needed <- unit_relvar *
        var_ratio *
        (1 + icc * (n2_opt - 1)) /
        (n2_opt * cv^2)
      n1_opt <- n1_eff_needed / resp_rate_psu
      total_cost <- fixed_cost + C1 * n1_opt + C2 * n1_opt * n2_opt
      cv_achieved <- cv
    }
  } else if (!is.null(n_psu)) {
    n1_opt <- n_psu
    n1_eff <- n_psu * resp_rate_psu
    if (!is.null(budget)) {
      n2_opt <- (var_budget - C1 * n_psu) / (C2 * n_psu)
      if (n2_opt < 1) {
        stop(
          "budget affords fewer than one unit per PSU for the given fixed stage-1 size",
          call. = FALSE
        )
      }
      cv_achieved <- sqrt(
        unit_relvar * var_ratio / (n1_eff * n2_opt) * (1 + icc * (n2_opt - 1))
      )
      total_cost <- budget
    } else {
      cv_floor <- sqrt(unit_relvar * var_ratio * icc / (n_psu * resp_rate_psu))
      if (cv_floor >= cv) {
        required_n_psu <- ceiling(unit_relvar * var_ratio * icc / (cv^2 * resp_rate_psu))
        stop(
          sprintf(
            "n_cluster(): target CV %.4g is below the achievable floor %.4g at n_psu = %d; increase n_psu to at least %d, or relax target CV above %.4g",
            cv, cv_floor, as.integer(n_psu), required_n_psu, cv_floor
          ),
          call. = FALSE
        )
      }
      n2_opt <- (1 - icc) / (cv^2 * n1_eff / (unit_relvar * var_ratio) - icc)
      if (n2_opt <= 0) {
        stop(
          "target CV is too small for the given fixed stage-1 size and parameters",
          call. = FALSE
        )
      }
      if (n2_opt < 1) {
        n2_opt <- 1
        cv_achieved <- sqrt(unit_relvar * var_ratio / n1_eff)
      } else {
        cv_achieved <- cv
      }
      total_cost <- fixed_cost + C1 * n_psu + C2 * n_psu * n2_opt
    }
  } else {
    n2_opt <- n_per_psu
    if (!is.null(budget)) {
      n1_opt <- var_budget / (C1 + C2 * n_per_psu)
      n1_eff <- n1_opt * resp_rate_psu
      cv_achieved <- sqrt(
        unit_relvar * var_ratio / (n1_eff * n2_opt) * (1 + icc * (n2_opt - 1))
      )
      total_cost <- budget
    } else {
      n1_eff_needed <- unit_relvar *
        var_ratio *
        (1 + icc * (n2_opt - 1)) /
        (n2_opt * cv^2)
      n1_opt <- n1_eff_needed / resp_rate_psu
      total_cost <- fixed_cost + C1 * n1_opt + C2 * n1_opt * n2_opt
      cv_achieved <- cv
    }
  }

  if (!is.null(budget) && n1_opt < 1) {
    stop(
      sprintf(
        "'budget' is too small for any realizable design: one PSU of size %.3g costs %.4g plus fixed_cost = %.4g",
        n2_opt, C1 + C2 * n2_opt, fixed_cost
      ),
      call. = FALSE
    )
  }

  # Back to gross units: the field team visits n2_opt / resp_rate units to
  # obtain n2_opt responses.
  n2_opt <- n2_opt / resp_rate
  .check_expected_take(n2_opt, resp_rate, "n_per_psu", "resp_rate")

  operational <- .op_cluster_2stage(
    stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu,
    n_per_psu_gross, resp_rate_psu, resp_rate, fixed_cost, cont_m = n2_opt
  )

  n_vec <- c(n_psu = n1_opt, n_per_psu = n2_opt)
  total_n <- prod(n_vec)

  params <- list(
    stage_cost = c(cost_psu = stage_cost[1L], cost_ssu = stage_cost[2L]),
    icc = c(icc_psu = icc[1L]),
    unit_relvar = unit_relvar,
    var_ratio = c(var_ratio_psu = var_ratio[1L]),
    resp_rate_psu = resp_rate_psu,
    resp_rate = resp_rate
  )
  if (!is.null(cv)) {
    params$cv <- cv
  }
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) {
    params$n_psu <- n_psu
  }
  if (!is.null(n_per_psu_gross)) {
    params$n_per_psu <- n_per_psu_gross
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 2L,
    total_n = total_n,
    cv = cv_achieved,
    cost = total_cost,
    params = params,
    operational = operational
  )
}

#' @keywords internal
#' @noRd
.n_cluster_3stage <- function(
  stage_cost,
  icc,
  unit_relvar,
  var_ratio,
  cv,
  budget,
  n_psu,
  n_per_psu,
  n_per_ssu,
  resp_rate_psu,
  resp_rate_ssu = 1,
  resp_rate = 1,
  fixed_cost = 0
) {
  C1 <- stage_cost[1L]
  # As two-stage, solved in realized units. Total cost is unchanged.
  C2 <- stage_cost[2L] / resp_rate_ssu
  C3 <- stage_cost[3L] / (resp_rate_ssu * resp_rate)
  n_per_psu_gross <- n_per_psu
  n_per_ssu_gross <- n_per_ssu
  if (!is.null(n_per_psu)) n_per_psu <- n_per_psu * resp_rate_ssu
  if (!is.null(n_per_ssu)) n_per_ssu <- n_per_ssu * resp_rate
  delta1 <- icc[1L]
  delta2 <- icc[2L]
  k1 <- var_ratio[1L]
  k2 <- var_ratio[2L]
  var_budget <- if (!is.null(budget)) budget - fixed_cost else NULL

  .cv3 <- function(n1e, n2v, n3v) {
    sqrt(
      unit_relvar / (n1e * n2v * n3v) *
        (k1 * delta1 * n2v * n3v + k2 * (1 + delta2 * (n3v - 1)))
    )
  }

  solve_for <- if (!is.null(n_psu) && !is.null(n_per_psu)) {
    "n3"
  } else if (!is.null(n_psu)) {
    "n2"
  } else {
    "n1"
  }

  n3 <- if (!is.null(n_per_ssu)) {
    n_per_ssu
  } else if (solve_for != "n3") {
    n3_free <- if (!is.null(n_per_psu)) {
      # Restricted optimum: a fixed n2 changes the final-stage take.
      a_fixed <- k1 * delta1 + k2 * delta2 / n_per_psu
      g_fixed <- k2 * (1 - delta2) / n_per_psu
      sqrt(g_fixed * (C1 + C2 * n_per_psu) /
             (a_fixed * C3 * n_per_psu))
    } else {
      sqrt((1 - delta2) / delta2 * C2 / C3)
    }
    if (n3_free < 1) {
      warning(
        "cost-optimal n_per_ssu is below 1 (high 'icc_ssu' relative to stage costs); clamped to 1",
        call. = FALSE
      )
      n3_free <- 1
    }
    n3_free
  }

  n2 <- if (!is.null(n_per_psu)) {
    n_per_psu
  } else if (solve_for != "n2") {
    n2_free <- if (!is.null(n_per_ssu)) {
      sqrt(
        k2 * (1 + delta2 * (n3 - 1)) * C1 /
          (n3 * k1 * delta1 * (C2 + C3 * n3))
      )
    } else {
      1 / n3 * sqrt((1 - delta2) / delta1 * C1 / C3 * k2 / k1)
    }
    if (n2_free < 1) {
      warning(
        "cost-optimal n_per_psu is below 1 (high 'icc_psu' relative to stage costs); clamped to 1",
        call. = FALSE
      )
      n2_free <- 1
    }
    n2_free
  }

  n1 <- n_psu

  if (solve_for == "n1") {
    if (!is.null(budget)) {
      n1 <- var_budget / (C1 + C2 * n2 + C3 * n2 * n3)
      cv_achieved <- .cv3(n1 * resp_rate_psu, n2, n3)
      total_cost <- budget
    } else {
      n1_eff <- unit_relvar / (cv^2 * n2 * n3) *
        (k1 * delta1 * n2 * n3 + k2 * (1 + delta2 * (n3 - 1)))
      n1 <- n1_eff / resp_rate_psu
      total_cost <- fixed_cost +
        C1 * n1 + C2 * n1 * n2 + C3 * n1 * n2 * n3
      cv_achieved <- cv
    }
  } else if (solve_for == "n2") {
    n1_eff <- n_psu * resp_rate_psu
    if (!is.null(budget)) {
      n2 <- (var_budget / n_psu - C1) / (C2 + C3 * n3)
      if (n2 < 1) {
        stop(
          "budget affords fewer than one SSU per PSU for the given fixed stage sizes",
          call. = FALSE
        )
      }
      cv_achieved <- .cv3(n1_eff, n2, n3)
      total_cost <- budget
    } else {
      cv_floor <- sqrt(unit_relvar * k1 * delta1 / n1_eff)
      if (cv_floor >= cv) {
        required_n_psu <- ceiling(unit_relvar * k1 * delta1 / (cv^2 * resp_rate_psu))
        stop(
          sprintf(
            "n_cluster(): target CV %.4g is below the achievable floor %.4g at n_psu = %d; increase n_psu to at least %d, or relax target CV above %.4g",
            cv, cv_floor, as.integer(n_psu), required_n_psu, cv_floor
          ),
          call. = FALSE
        )
      }
      n2 <- k2 * (1 + delta2 * (n3 - 1)) /
        (n3 * (cv^2 * n1_eff / unit_relvar - k1 * delta1))
      if (n2 <= 0) {
        stop(
          "target CV is too small for the given fixed stage sizes and parameters",
          call. = FALSE
        )
      }
      if (n2 < 1) {
        n2 <- 1
        cv_achieved <- .cv3(n1_eff, 1, n3)
      } else {
        cv_achieved <- cv
      }
      total_cost <- fixed_cost +
        C1 * n_psu + C2 * n_psu * n2 + C3 * n_psu * n2 * n3
    }
  } else {
    n1_eff <- n_psu * resp_rate_psu
    if (!is.null(budget)) {
      n3 <- (var_budget - C1 * n_psu - C2 * n_psu * n2) / (C3 * n_psu * n2)
      if (n3 < 1) {
        stop(
          "budget affords fewer than one unit per SSU for the given fixed stage sizes",
          call. = FALSE
        )
      }
      cv_achieved <- .cv3(n1_eff, n2, n3)
      total_cost <- budget
    } else {
      cv_floor <- sqrt(unit_relvar * (k1 * delta1 + k2 * delta2 / n2) / n1_eff)
      if (cv_floor >= cv) {
        required_n_psu <- ceiling(
          unit_relvar * (k1 * delta1 + k2 * delta2 / n2) / (cv^2 * resp_rate_psu)
        )
        stop(
          sprintf(
            "n_cluster(): target CV %.4g is below the achievable floor %.4g at n_psu = %d, n_per_psu = %.0f; increase n_psu to at least %d, or relax target CV above %.4g",
            cv, cv_floor, as.integer(n_psu), n2, required_n_psu, cv_floor
          ),
          call. = FALSE
        )
      }
      denom <- n2 * (cv^2 * n1_eff / unit_relvar - k1 * delta1) - k2 * delta2
      n3 <- k2 * (1 - delta2) / denom
      if (n3 <= 0) {
        stop(
          "target CV is too small for the given fixed stage sizes and parameters",
          call. = FALSE
        )
      }
      if (n3 < 1) {
        n3 <- 1
        cv_achieved <- .cv3(n1_eff, n2, 1)
      } else {
        cv_achieved <- cv
      }
      total_cost <- fixed_cost +
        C1 * n_psu + C2 * n_psu * n2 + C3 * n_psu * n2 * n3
    }
  }

  if (!is.null(budget) && n1 < 1) {
    stop(
      sprintf(
        "'budget' is too small for any realizable design: one PSU with n_per_psu = %.3g and n_per_ssu = %.3g costs %.4g plus fixed_cost = %.4g",
        n2, n3, C1 + C2 * n2 + C3 * n2 * n3, fixed_cost
      ),
      call. = FALSE
    )
  }

  n2 <- n2 / resp_rate_ssu
  n3 <- n3 / resp_rate
  .check_expected_take(n2, resp_rate_ssu, "n_per_psu", "resp_rate_ssu")
  .check_expected_take(n3, resp_rate, "n_per_ssu", "resp_rate")

  operational <- .op_cluster_3stage(
    stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu,
    n_per_psu_gross, n_per_ssu_gross,
    resp_rate_psu, resp_rate_ssu, resp_rate, fixed_cost,
    cont_m = n2
  )

  n_vec <- c(n_psu = n1, n_per_psu = n2, n_per_ssu = n3)
  total_n <- prod(n_vec)

  params <- list(
    stage_cost = c(cost_psu = stage_cost[1L], cost_ssu = stage_cost[2L],
                   cost_tsu = stage_cost[3L]),
    icc = c(icc_psu = delta1, icc_ssu = delta2),
    unit_relvar = unit_relvar,
    var_ratio = c(var_ratio_psu = k1, var_ratio_ssu = k2),
    resp_rate_psu = resp_rate_psu,
    resp_rate_ssu = resp_rate_ssu,
    resp_rate = resp_rate
  )
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(n_per_psu_gross)) params$n_per_psu <- n_per_psu_gross
  if (!is.null(n_per_ssu_gross)) params$n_per_ssu <- n_per_ssu_gross
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost

  .new_svyplan_cluster(
    n = n_vec,
    stages = 3L,
    total_n = total_n,
    cv = cv_achieved,
    cost = total_cost,
    params = params,
    operational = operational
  )
}
