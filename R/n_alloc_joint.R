#' Joint constrained stratified allocation
#'
#' Allocate a sample across atomic strata to meet precision requirements for
#' multiple indicators and domains at once, or to minimize a priority-weighted
#' relative variance, summed or at its largest, within a fixed budget. Supports element allocation,
#' fixed-take multistage allocation, and two-stage allocation with a free
#' within-PSU take. A PSU register adds certainty-aware two-stage planning.
#'
#' @param frame For the default method, a stratum-level data frame describing
#'   atomic allocation strata. Required columns are unique, non-missing
#'   `stratum` and positive finite `N` (ultimate units). Optional columns
#'   include `max_weight`, `min_n_stratum`, `deff`, `resp_rate`,
#'   `eligibility_rate`, `eligible_per_unit`, and complete domain columns named
#'   by `targets$domain` or `objective$domain`. Indicator values (`p`, `mean`,
#'   `sd`, `var`) belong in `measures`, not `frame`.
#'
#'   * For element designs, use the stratum and ultimate-unit totals without
#'     numbered stage columns. Optional `unit_cost` prices each unit and
#'     `take_all` marks strata to enumerate.
#'   * For fixed takes, supply whole stage populations `N_1`, ...,
#'     `N_{S-1}`, whole takes `take_2`, ..., `take_S`, and positive stage
#'     costs `cost_1`, ..., `cost_S`. Stage costs replace `unit_cost`.
#'   * For a free two-stage take, supply `N_1` or `psu` and positive
#'     `cost_1` and `cost_2`. Omit `take_2` or set it to `NA` in the
#'     strata whose takes should be optimized. Bound the take with
#'     `take_2_min`, `take_2_max`, and `take_2_step`, as arguments or as
#'     frame columns of the same names.
#'   * For a PSU register, supply `psu` in place of `N_1`. A fixed take uses
#'     `take_2` and can use element costs or the pair `cost_1`, `cost_2`.
#'     A free take uses the stage costs and bounds above.
#'
#'   For `svyplan_prec` objects, a precision result from [prec_alloc_joint()].
#' @param measures Long data frame keyed by `stratum` and indicator `name`.
#'   Supply one row for each indicator and stratum covered by its targets or
#'   objective. Each row specifies either `p` in `[0, 1]`, or `mean` and one of
#'   non-negative `sd` or `var`. Optional columns `deff`, `resp_rate`,
#'   `eligibility_rate`, `eligible_per_unit`, `icc_1`, ..., `icc_{S-1}`,
#'   `var_ratio_1`, ..., `var_ratio_{S-1}`, and `resp_rate_1`, ...,
#'   `resp_rate_{S-1}` override stratum defaults in `frame` and scalar
#'   arguments where supported. A scalar `deff`, `resp_rate`,
#'   `eligibility_rate` or `eligible_per_unit` other than 1 is refused when
#'   every row that a target or the objective uses takes the value from its
#'   own column or the `frame` column. Register fits use `icc_1` and ultimate-unit
#'   `resp_rate`. They require eligibility as an argument or frame column,
#'   a variance ratio of 1, and no upper-stage response-rate columns.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param targets Optional long data frame with one precision target per row.
#'   Requires `name` and exactly one positive `cv`, `moe`, or `rmoe`. Optional
#'   columns are `domain` (default `".overall"`), `level` (`NA` for
#'   `".overall"`), `alpha`, `df`, and `constraint`. Leaving out `level` applies
#'   a domain row to every level of that domain in `frame`. Required unless both
#'   `objective` and `budget` are supplied.
#' @param objective Optional estimates whose priority-weighted relative variance
#'   is minimized within `budget`. Either a character vector of indicator names
#'   (overall domain, equal priority) or a data frame with `name` and optional
#'   `domain`, `level`, `priority`, and `component` columns. Requires `budget`.
#' @param budget Optional positive field budget. Must be supplied with
#'   `objective`.
#' @param objective_rule How `objective` combines its components. `"sum"`
#'   (default) minimizes the sum of `priority * relvar` over the components.
#'   `"worst"` minimizes the largest `priority * relvar`, so with equal
#'   priorities it gives the smallest CV that every component is held to.
#'   Under `"worst"` every `priority` must be positive, `control$gap` is
#'   refused, and a PSU register is not available.
#' @param psu Optional PSU register for certainty-aware two-stage allocation,
#'   one row per PSU, with `stratum` and `N` (ultimate units in the PSU). PSU
#'   sizes must sum to `frame$N` within each stratum. Supply unique
#'   `psu_id` values to join the plan back to the register. Optional
#'   `certainty` flags PSUs that must be included. `zone_order` requires
#'   zoning and orders PSUs within each stratum without missing values or ties.
#'   The default zone order is decreasing PSU size.
#' @param certainty_cutoff Optional inclusion probability threshold in `(0, 1]`
#'   at which a register PSU is taken with certainty. Scalar, length-`nrow(frame)`
#'   vector, or `frame$certainty_cutoff`. Requires `psu`.
#' @param n_psu_per_zone Optional `1` or `2`, the number of PSUs drawn per zone
#'   in each stratum remainder. Requires `psu`. A selection zone is a group
#'   of noncertainty PSUs, formed by ordering the remainder of a stratum
#'   (`psu$zone_order`, by default decreasing size) and cutting it by
#'   cumulative size. No geographic variable is needed. With two PSUs per
#'   zone the zones are the variance groups. With one PSU per zone the zones
#'   are collapsed into variance groups of two or three, which can span
#'   input strata, and `$psu$.variance_group` gives each PSU's group.
#' @param control Optional list with `max_nodes` (default `5000L`) and `gap`
#'   (default `0`) controlling the whole-unit search.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of per-unit
#'   stratum costs, overriding `frame$unit_cost`. Not used with stage costs.
#' @param alpha Significance level, default `0.05`.
#' @param deff Default design effect multiplier (`> 0`), scalar or length
#'   `nrow(frame)`.
#' @param resp_rate Default expected response rate among eligible units in
#'   `(0, 1]`.
#' @param eligibility_rate Default share of issued units in scope in `(0, 1]`.
#'   For a register, supply the rate by argument or in `frame`.
#' @param eligible_per_unit Default eligible units per sampled unit (`> 0`).
#'   Not available with `psu`.
#' @param df Optional degrees of freedom for turning `moe` and `rmoe` targets
#'   into variances with a `t` quantile.
#' @param min_n_stratum Optional minimum sample size per stratum. Scalar,
#'   length-`nrow(frame)` vector, or `frame$min_n_stratum`.
#' @param min_n_domain Optional sample floors on the levels of one domain
#'   column of `frame`, as a data frame with columns `domain` (that
#'   column's name), `level` and `min_n`. A floor counts issued units
#'   before response. The continuous design meets `min_n`, and the whole
#'   design meets it rounded up. Element frames only, under both objective
#'   rules and in cost mode.
#' @param take_2_min Optional minimum within-PSU take for free-take two-stage
#'   optimization. Positive whole number, scalar or length-`nrow(frame)` vector,
#'   overriding `frame$take_2_min`. Default is `1`.
#' @param take_2_max Optional maximum within-PSU take for free-take two-stage
#'   optimization. Positive whole number, scalar or length-`nrow(frame)` vector,
#'   overriding `frame$take_2_max`. Defaults to `N / N_1` (or the smallest
#'   noncertainty PSU size when `psu` is supplied).
#' @param take_2_step Optional integer step between candidate within-PSU
#'   takes for free-take two-stage optimization. Positive whole number, scalar
#'   or length-`nrow(frame)` vector, overriding `frame$take_2_step`. Default
#'   is `1`.
#' @param fpc Finite population correction. `"unit"` (default) uses
#'   the ultimate-unit correction. `"stage"` uses stage-by-stage corrections
#'   with `N_1`, ..., `N_{S-1}`. `"none"` applies no correction.
#'   Leave this at its default for a PSU register, which uses its own
#'   certainty and remainder variance model.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_alloc_joint` object inheriting from `svyplan_alloc`.
#'   * `stages` counts field stages: 1 for elements, the stage count of a
#'     cluster frame, and 2 for a PSU register. In a register fit,
#'     `params$stages` is 1 because clustering enters the precision model
#'     through a design effect.
#'   * `n` is the continuous total. In `detail`, `n` and `n_int` give
#'     continuous and whole ultimate units by stratum. Cluster frames also
#'     report PSU counts `n_1` and `n_1_int`, and stage takes `take_2`, ...,
#'     `take_S` with their whole values `take_2_int`, ..., `take_S_int`.
#'     The whole field design is `n_1_int` PSUs with `take_s_int` units at
#'     each later stage.
#'   * `constraints` reports each target, achieved precision, and pass
#'     indicator at the continuous design. `.multiplier` prices each target
#'     in cost, and `.sensitivity` is the change in minimum cost per unit
#'     relaxation of the target. In budget mode, `.sensitivity` holds the
#'     objective fixed, and `.objective_sensitivity` gives the change in the
#'     objective per unit relaxation with the budget fixed, so
#'     `.objective_sensitivity * 0.1 * .target` is the first-order effect of
#'     relaxing a target by 10 percent. It is 0 for a target with slack and
#'     for every target when the budget does not bind. An active target
#'     whose price it shares with another active target has no slope of its
#'     own, so both of its slopes are `NA`.
#'     `$optimization$slopes` reads `"available"`, or names why slopes are
#'     withheld. `"budget_at_minimum"` (the budget equals the cheapest cost
#'     meeting the targets), `"register"`, and for free takes
#'     `"not_optimal"` (the continuous design is only `"feasible"`) and
#'     `"multipliers_not_verified"` (no multipliers price the full free-take
#'     problem at the design) withhold every slope.
#'     `"objective_multiplier_not_positive"` and `"not_identifiable"` (the
#'     objective shares its price with the targets) withhold
#'     `.objective_sensitivity` and `$optimization$budget_sensitivity`, and
#'     `"budget_not_binding"` withholds `budget_sensitivity`. A refit with
#'     the changed target or budget then measures its effect. `operational` gives the whole
#'     total, cost, constraint results, `all_pass`, and search status.
#'     Whole counts are stored as doubles and can represent integers
#'     exactly up to 2^53.
#'   * `optimization` describes the continuous solve and holds the
#'     whole-unit search record in `search`. Fixed-budget runs also report
#'     `objective` and `objective_value`.
#'   * With `min_n_domain`, `domain_floors` gives each floor's `required`
#'     and `achieved` count and `.pass` at the continuous design, and
#'     `operational$domain_floors` the same at the whole design.
#'     `.sensitivity` is the change in minimum cost per unit of floor, and
#'     in budget mode `.objective_sensitivity` the change in the objective
#'     per unit of floor at the fixed budget. Both are 0 for a floor with
#'     slack, and follow `$optimization$slopes` as the target slopes do.
#'   * For registers, `detail$n_psu_certain` and `detail$n_psu_draw` count
#'     the certainty PSUs and the remainder draws, and `detail$n_1_int` is
#'     their sum. `detail$take_2` and `take_2_int` give the remainder take.
#'     Zoned fits add `detail$n_zone`. The `psu` table gives `certainty` and
#'     `take_2` per PSU, plus `psu_id` if supplied, `.zone` for zoned fits,
#'     and `.variance_group` when drawing one PSU per zone. On certainty
#'     rows `take_2` is the PSU's own take. On remainder rows it is the take
#'     if the PSU is drawn.
#'   * `rates` gives the rates the fit used, one row per stratum and
#'     indicator that a target or the objective reads, with `stratum`,
#'     `name`, `resp_rate_1`, ..., `resp_rate_{S-1}` for the stages of the
#'     precision model, `resp_rate` and `eligibility_rate`. Each rate comes
#'     from the measure row, then the stratum's `frame` column, then the
#'     argument.
#'   * `targets`, `binding`, and `params` retain the requirements,
#'     binding constraints, and inputs.
#'
#' @details
#' ## Minimum-cost and fixed-budget problems
#'
#' Without `objective`, `n_alloc_joint()` finds the cheapest allocation meeting
#' every precision requirement in `targets`. With `objective` and `budget`, it
#' minimizes the priority-weighted sum of squared CVs across the objective
#' components subject to `budget` and any hard requirements in `targets`.
#' Replanning from [prec_alloc_joint()] uses the achieved metrics as
#' targets in minimum-cost mode and the assessed cost as the budget in
#' objective mode. It optimizes the resulting problem again.
#'
#' ## Hold every component to one level
#'
#' Under the sum rule one estimate can lose precision while the others gain.
#' With `objective_rule = "worst"`, the budget buys the smallest level `t`
#' with `priority * relvar <= t` for every component, the weighted
#' Chebyshev allocation of Willems and Münnich (2026). `$objective_value`
#' is that level, `.contribution` in `$objective` is each component's
#' `priority * relvar`, and `.binding` marks the components at the level.
#' A second stage then minimizes the priority-weighted sum within the
#' budget, with every component held at the level. It spends any budget the
#' level leaves, and it returns a design that no other affordable design
#' improves in every component (Steuer and Choo, 1983, Theorem 4.6).
#'
#' `$optimization$objective_weights` gives sum-rule priorities, summing to 1,
#' under which the sum rule returns the same continuous design. A component
#' with weight 0 does not drive the design at this budget. The prices and
#' slopes come from the first stage, so `budget_sensitivity` and
#' `.objective_sensitivity` are changes in the level. The whole design comes
#' from a bisection on the level, each trial an exact minimum-cost search
#' with every component capped. `$optimization$search` records the proven
#' `bound`, the `gap`, the number of `trials`, and the status of the whole
#' second stage in `second_stage`. `$optimization$second_stage` gives that of
#' the continuous one.
#'
#' With free takes, the continuous level carries the same kind of
#' certificate as the free-take sum rule. Weak duality gives
#' `$optimization$objective_lower_bound`, and the design is `"optimal"` when
#' its level lies within a relative `1e-6` of that bound and `"feasible"`
#' otherwise. A free-take budget fit reports `.objective_sensitivity` and
#' `budget_sensitivity`, under either rule, only when it is `"optimal"` and
#' its multipliers are verified against the full free-take problem, so that
#' they do not price only one choice of takes.
#'
#' ## Multistage and register designs
#'
#' Fixed-take multistage designs keep later-stage takes fixed and optimize the
#' first-stage PSU counts (`n_1`) per stratum. Population columns (`N`, `N_1`,
#' ..., `N_{S-1}`) are stratum totals, while takes (`take_2`, ..., `take_S`)
#' count units per parent unit. Stage response rates (`resp_rate_1`, ...,
#' `resp_rate_{S-1}`, `resp_rate`) act at the stage they name.
#'
#' With a `psu` register, `n_alloc_joint()` separates certainty PSUs from the
#' remainder, prices certainty visits alongside drawn PSUs, and chooses the
#' whole-unit field design by branch-and-bound search over per-stratum designs.
#' A cluster frame or register cannot use `take_all`.
#'
#' ## Read the field design and its certificate
#'
#' Use `$operational` for fieldwork. Its whole counts may differ from the
#' continuous design. Read `$operational$status` and `$optimization$search`
#' to see whether the search proved optimality, reached the requested gap,
#' or stopped at its node limit.
#'
#' For a non-register minimum-cost fit with fixed takes, the continuous
#' solver checks KKT conditions. With a free take, the continuous
#' minimum-cost certificate is a duality gap. `$optimization$classification` is
#' `"optimal"` when the design costs within a relative `1e-6` of the Lagrangian
#' dual bound in `$optimization$dual_value`, and `"feasible"` otherwise, with
#' the gap in `$optimization$relative_duality_gap`. That gap is
#' `(cost - dual_value) / cost`, so a feasible design meets every target and
#' its cost exceeds the optimum by at most the gap times its own cost.
#'
#' In non-register free-take budget mode, a budget below the dual cost
#' bound is refused. A budget between the bound and the cheapest design
#' found raises an error of
#' class `svyplan_plan_unknown`. Otherwise the design is the most precise
#' the solver found within the budget, and its certificate is a bound on
#' precision. Weak duality gives `$optimization$objective_lower_bound`, below
#' which no design can go that costs at most the budget, or at most the
#' returned design when that costs more. The design is `"optimal"` when its
#' objective lies within a relative `1e-6` of that bound
#' (`$optimization$objective_gap`). It is `"feasible"`, with both fields
#' missing, when the solver has no multiplier for the objective, as at a
#' budget that only the cheapest design fits.
#'
#' Every budget fit, continuous or whole and with fixed or free takes,
#' buys its design within the budget itself, up to floating-point rounding
#' of the budget. The cheapest design that meets the targets is accepted
#' when it costs at most a relative `1e-6` more than the budget, and the
#' search then stays at that cost.
#'
#' See `vignette("joint-allocation")` for the three input tables, assessment
#' of whole designs, and an example of a fixed-budget grid.
#'
#' @references
#' Bethel, J. (1989). Sample allocation in multivariate surveys. *Survey
#' Methodology*, 15(1), 47--57.
#'
#' Steuer, R. E. and Choo, E.-U. (1983). An interactive weighted Tchebycheff
#' procedure for multiple objective programming. *Mathematical Programming*,
#' 26(3), 326--344.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer. Chapter 5.
#'
#' Willems, F. and Münnich, R. (2026). Sensitivity analysis for
#' multivariate optimal sample allocation. *arXiv preprint*
#' arXiv:2610.26093.
#'
#' @family stratified design functions
#' @seealso [prec_alloc_joint()] for the inverse assessment, [n_alloc()] for
#'   classic single-indicator stratum allocation, [merge_psus()] for merging
#'   small PSUs before register planning, [predict.svyplan()] for budget
#'   and take frontiers, and [tradeoff_alloc()] for trading the precision of
#'   one objective component for another's.
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("North urban", "North rural", "South urban"),
#'   region = c("North", "North", "South"),
#'   N = c(1000, 1800, 1200)
#' )
#' measures <- data.frame(
#'   stratum = rep(frame$stratum, 2),
#'   name = rep(c("coverage", "income"), each = 3),
#'   p = c(0.5, 0.4, 0.6, NA, NA, NA),
#'   mean = c(NA, NA, NA, 50, 55, 60),
#'   sd = c(NA, NA, NA, 10, 12, 15)
#' )
#' targets <- data.frame(
#'   name = c("coverage", "income"),
#'   domain = c(".overall", "region"),
#'   level = c(NA, "North"),
#'   cv = c(0.06, NA),
#'   moe = c(NA, 3)
#' )
#' fit <- n_alloc_joint(frame, measures, targets = targets)
#' fit
#' prec_alloc_joint(fit, n = fit$detail$n_int)
#'
#' # Fixed-budget weighted objective subject to the same targets
#' bought <- n_alloc_joint(
#'   frame, measures, targets = targets,
#'   objective = "income", budget = 400
#' )
#' bought$objective
#'
#' # The same budget held to the smallest CV both estimates meet
#' level <- n_alloc_joint(
#'   frame, measures, targets = targets,
#'   objective = c("coverage", "income"), budget = 400,
#'   objective_rule = "worst"
#' )
#' level$objective
#'
#' # At least 150 interviews in the South, which the budget alone would not give
#' floored <- n_alloc_joint(
#'   frame, measures, targets = targets,
#'   objective = c("coverage", "income"), budget = 400,
#'   min_n_domain = data.frame(domain = "region", level = "South", min_n = 150)
#' )
#' floored$operational$domain_floors
#'
#' # Fixed and free takes for a two-stage coverage plan
#' cluster <- data.frame(
#'   stratum = c("A", "B"), N = c(20000, 8000), N_1 = c(400, 160),
#'   cost_1 = 400, cost_2 = 40, take_2 = c(12, 8)
#' )
#' coverage <- data.frame(
#'   stratum = c("A", "B"), name = "coverage", p = c(0.5, 0.4), icc_1 = 0.05
#' )
#' goal <- data.frame(name = "coverage", cv = 0.10)
#' held <- n_alloc_joint(cluster, coverage, targets = goal)
#' open <- n_alloc_joint(transform(cluster, take_2 = NULL), coverage,
#'                       targets = goal, take_2_min = 2, take_2_max = 20)
#' held$detail[c("n_1_int", "take_2")]
#' open$detail[c("n_1_int", "take_2_int")]
#'
#' # A register with two remainder PSUs drawn per zone
#' register <- data.frame(stratum = "A", psu_id = seq_len(100), N = 1000)
#' zoned <- n_alloc_joint(
#'   data.frame(stratum = "A", N = 100000, take_2 = 10),
#'   coverage[1, ], targets = transform(goal, cv = 0.08),
#'   psu = register, n_psu_per_zone = 2
#' )
#' zoned$detail[c("n_psu_certain", "n_psu_draw", "n_zone", "n_int")]
#' prec_alloc_joint(zoned)
#'
#' @export
n_alloc_joint <- function(frame, ...) {
  if (!missing(frame)) {
    .refuse_result_as_input(frame, "frame", "n_alloc_joint")
    .res <- .dispatch_plan(frame, "frame", n_alloc_joint.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_alloc_joint")
}

#' @rdname n_alloc_joint
#' @export
n_alloc_joint.default <- function(
  frame,
  measures = NULL,
  ...,
  targets = NULL,
  objective = NULL,
  budget = NULL,
  objective_rule = c("sum", "worst"),
  psu = NULL,
  certainty_cutoff = NULL,
  n_psu_per_zone = NULL,
  control = NULL,
  unit_cost = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  df = NULL,
  min_n_stratum = NULL,
  min_n_domain = NULL,
  take_2_min = NULL,
  take_2_max = NULL,
  take_2_step = NULL,
  fpc = c("unit", "stage", "none"),
  plan = NULL
) {
  .plan <- .merge_plan_args(
    plan,
    n_alloc_joint.default,
    match.call(),
    environment()
  )
  if (!is.null(.plan)) {
    return(do.call(n_alloc_joint.default, c(.plan, list(...))))
  }
  .check_alloc_classic_args(list(...), "n_alloc()")
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(frame, "frame")
  .refuse_unnumbered_stage_columns(measures, "measures")
  .refuse_unnumbered_stage_columns(psu, "psu")
  .check_register_args(certainty_cutoff, n_psu_per_zone, frame, psu)
  .check_eligible_scope(eligible_per_unit, frame, measures, psu)
  .refuse_joint_argument_beside_columns(
    frame, measures, targets, objective, alpha, df,
    resp_rate = resp_rate, deff = deff,
    eligibility_rate = eligibility_rate, eligible_per_unit = eligible_per_unit
  )
  min_n_stratum <- .resolve_min_n_stratum(min_n_stratum, frame)
  fpc <- match.arg(fpc)
  objective_rule <- match.arg(objective_rule)
  .check_objective_rule(objective_rule, objective, control)

  if (is.null(measures)) {
    stop("joint constrained allocation requires 'measures'", call. = FALSE)
  }
  if (is.null(objective)) {
    if (is.null(targets)) {
      stop(
        "'measures' and 'targets' must be supplied together",
        call. = FALSE
      )
    }
    if (!is.null(budget)) {
      stop(
        "'budget' requires 'objective': name the estimates whose relative variance the budget should buy down, or drop 'budget' for the minimum-cost design meeting 'targets'",
        call. = FALSE
      )
    }
  } else {
    if (is.null(budget)) {
      stop(
        "'objective' requires 'budget': an objective always improves with more sample, so it needs a budget to trade against",
        call. = FALSE
      )
    }
    check_scalar(budget, "budget")
  }
  if ("alloc_measure" %in% names(frame)) {
    stop(
      "'alloc_measure' is not used for joint constrained allocation, it belongs to n_alloc(alloc = \"power\")",
      call. = FALSE
    )
  }
  .alloc_refuse_step(frame, "joint constrained allocation")
  ft_spec <- .resolve_freetake_spec(
    frame = frame, measures = measures, psu = psu,
    take_2_min = take_2_min, take_2_max = take_2_max, take_2_step = take_2_step,
    unit_cost = unit_cost, fpc = fpc
  )

  if (!is.null(min_n_domain) && (!is.null(psu) || isTRUE(ft_spec$active))) {
    stop("joint 'min_n_domain' is available for element frames",
         call. = FALSE)
  }
  if (!is.null(psu)) {
    if (!identical(fpc, "unit")) {
      stop(
        "'fpc' is not available with a PSU register: the certainty split and its size-proportional remainder need their own variance decomposition, which is not derived here. Drop 'fpc', or drop 'psu' and give 'N_1' instead",
        call. = FALSE
      )
    }
    .refuse_objective_rule_route(objective_rule, "register")
    return(.psu_result(.n_alloc_psu(
      frame = frame,
      psu = psu,
      measures = measures,
      targets = targets,
      unit_cost = unit_cost,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      min_n_stratum = min_n_stratum,
      objective = objective,
      budget = budget,
      df = df,
      certainty_cutoff = certainty_cutoff,
      n_psu_per_zone = n_psu_per_zone,
      control = control,
      eligibility_rate = eligibility_rate,
      ft_spec = ft_spec,
      take_2_min = take_2_min,
      take_2_max = take_2_max,
      take_2_step = take_2_step
    )))
  }

  sc <- .search_control(control)
  if (isTRUE(ft_spec$active)) {
    fit <- .n_alloc_bethel_freetake(
      frame = frame,
      measures = measures,
      targets = targets,
      ft_spec = ft_spec,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      min_n_stratum = min_n_stratum,
      objective = objective,
      budget = budget,
      df = df,
      fpc = fpc,
      eligible_per_unit = eligible_per_unit,
      eligibility_rate = eligibility_rate,
      take_2_min = take_2_min,
      take_2_max = take_2_max,
      take_2_step = take_2_step,
      max_nodes = sc$max_nodes,
      gap = sc$gap,
      objective_rule = objective_rule
    )
    fit$params$control <- control
    return(fit)
  }
  fit <- .n_alloc_bethel(
    frame = frame,
    measures = measures,
    targets = targets,
    unit_cost = unit_cost,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit,
    eligibility_rate = eligibility_rate,
    min_n_stratum = min_n_stratum,
    objective = objective,
    budget = budget,
    df = df,
    fpc = fpc,
    .finite_prop_var = TRUE,
    max_nodes = sc$max_nodes,
    gap = sc$gap,
    objective_rule = objective_rule,
    min_n_domain = min_n_domain
  )
  fit$params$control <- control
  fit
}

#' Refuse the worst rule where it is not defined
#' @keywords internal
#' @noRd
.check_objective_rule <- function(objective_rule, objective, control) {
  if (!identical(objective_rule, "worst")) {
    return(invisible(NULL))
  }
  if (is.null(objective)) {
    stop(
      "objective_rule = \"worst\" needs 'objective': name the components whose largest weighted relative variance the budget should minimize",
      call. = FALSE
    )
  }
  priority <- if (is.data.frame(objective)) objective[["priority"]]
  if (is.numeric(priority) && any(priority == 0, na.rm = TRUE)) {
    stop(
      "under objective_rule = \"worst\" every objective 'priority' must be positive. Leave out a component that is only reported, and read its precision from prec_alloc_joint()",
      call. = FALSE
    )
  }
  gap <- if (is.list(control)) control[["gap"]]
  if (is.numeric(gap) && length(gap) == 1L && !is.na(gap) && gap != 0) {
    stop(
      "'control$gap' is refused under objective_rule = \"worst\": the whole search bisects on the worst level and stops only when its bracket closes. Its measured gap is in $optimization$search$gap",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Refuse the worst rule on a route that has no worst-rule solver
#' @keywords internal
#' @noRd
.refuse_objective_rule_route <- function(objective_rule, route) {
  if (!identical(objective_rule, "worst")) {
    return(invisible(NULL))
  }
  stop(
    switch(
      route,
      register = "objective_rule = \"worst\" is not available with a PSU register, whose design is planned for one certainty configuration at a time. Use objective_rule = \"sum\""
    ),
    call. = FALSE
  )
}

#' @rdname n_alloc_joint
#' @export
n_alloc_joint.svyplan_prec <- function(frame, ...) {
  x <- frame
  if (!identical(x$type, "alloc_joint")) {
    stop(
      "n_alloc_joint requires a svyplan_prec of type 'alloc_joint', from prec_alloc_joint()",
      call. = FALSE
    )
  }
  p <- x$params
  targets <- p$targets
  if (is.null(p$objective)) {
    for (column in c("cv", "moe", "rmoe")) {
      if (!column %in% names(targets)) {
        targets[[column]] <- NA_real_
      }
      targets[[column]][] <- NA_real_
    }
    matched <- match(targets$constraint, x$detail$constraint)
    if (anyNA(matched)) {
      stop(
        "precision result is missing stored constraint identifiers",
        call. = FALSE
      )
    }
    for (column in c("cv", "moe", "rmoe")) {
      rows <- x$detail$.metric[matched] == column
      targets[[column]][rows] <- x$detail$.achieved[matched][rows]
    }
    budget_arg <- NULL
  } else {
    budget_arg <- p$achieved$cost
    if (nrow(targets) == 0L) targets <- NULL
  }
  args <- list(
    frame = p$frame,
    measures = p$measures,
    targets = targets,
    objective = p$objective,
    budget = budget_arg,
    objective_rule = p$objective_rule %||% "sum",
    unit_cost = p$unit_cost,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    eligible_per_unit = p$eligible_per_unit %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    df = p$df,
    min_n_stratum = p$min_n_stratum,
    min_n_domain = p$min_n_domain,
    fpc = p$fpc %||% "unit"
  )
  if (!is.null(p$psu)) {
    args$psu <- p$psu
    args$certainty_cutoff <- p$certainty_cutoff
    args$n_psu_per_zone <- p$n_psu_per_zone
  }
  if (isTRUE(p$free_take)) {
    args$take_2_min <- p[["take_2_min"]]
    args$take_2_max <- p[["take_2_max"]]
    args$take_2_step <- p[["take_2_step"]]
  }
  args$control <- p[["control"]]
  do.call(
    n_alloc_joint.default,
    .roundtrip_args(args, list(...), n_alloc_joint.default)
  )
}

#' Report classic allocation arguments passed to joint allocation
#' @keywords internal
#' @noRd
.check_alloc_classic_args <- function(dots, replacement) {
  classic <- intersect(
    names(dots) %||% character(0),
    c("n", "cv", "moe", "rmoe", "domains", "alloc", "alloc_q", "min_n_domain")
  )
  if (length(classic) > 0L) {
    stop(
      sprintf(
        "classic allocation argument%s %s belong%s to %s",
        if (length(classic) > 1L) "s" else "",
        paste(sQuote(classic), collapse = ", "),
        if (length(classic) > 1L) "" else "s",
        replacement
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Report joint allocation arguments passed to classic allocation
#' @keywords internal
#' @noRd
.check_alloc_joint_args <- function(dots, replacement) {
  joint <- intersect(
    names(dots) %||% character(0),
    c("measures", "targets", "objective", "psu", "certainty_cutoff",
      "n_psu_per_zone")
  )
  if (length(joint) > 0L) {
    stop(
      sprintf(
        "joint allocation argument%s %s belong%s to %s",
        if (length(joint) > 1L) "s" else "",
        paste(sQuote(joint), collapse = ", "),
        if (length(joint) > 1L) "" else "s",
        replacement
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}
