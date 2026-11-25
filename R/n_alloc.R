#' Stratified sample allocation
#'
#' Distribute a sample across strata for a single indicator under a fixed total
#' `n`, a precision target (`cv`, `moe`, or `rmoe`), or a `budget`, using
#' Neyman, optimal, proportional, or Bankier power allocation. For multiple
#' indicators, overlapping domain targets, or a PSU register with certainty
#' selection, use [n_alloc_joint()].
#'
#' @param frame For the default method, a stratum-level data frame describing
#'   the population to sample. Each row represents one stratum. When a design
#'   stratifies by several variables (such as region and urbanicity), cross
#'   them into a single stratification column first so that each row maps to
#'   one population cell.
#'
#'   **Required columns**
#'   \describe{
#'     \item{`N`}{Population count of units in each stratum. Must be positive
#'       and finite.}
#'     \item{`sd` or `var`}{Within-stratum spread of the key variable, required
#'       except for a fixed `n` under `alloc = "proportional"` or `"power"`.
#'       Supply `sd` or `var`, not both. A frame with `p` and neither column
#'       derives its variance from `p`, as `N * p * (1 - p) / (N - 1)` in an
#'       element design under `fpc = "unit"` and `p * (1 - p)` in a cluster
#'       design or under `fpc = "none"`. Without `sd`, `var`, or `p`, fixed-`n`
#'       proportional and power rules allocate from stratum sizes alone and
#'       report `NA` precision.}
#'   }
#'
#'   **Optional columns**
#'   \describe{
#'     \item{`stratum`}{Stratum label. Defaults to row numbers. Must be unique,
#'       or unique within each domain when `domains` is set.}
#'     \item{Label columns}{Other character, factor, or logical columns are
#'       carried into `$detail` after `stratum`.}
#'     \item{`mean` or `p`}{Stratum population mean or proportion (in `[0, 1]`).
#'       Required when solving for `cv` or `rmoe`, and by `alloc = "power"`
#'       when `sd` or `var` is given.}
#'     \item{`alloc_measure`}{Positive stratum size measure used only by
#'       `alloc = "power"`. Defaults to `N * abs(mean)` (or `N * p`) when a
#'       variance is given, or `N` otherwise.}
#'     \item{`unit_cost`}{Positive per-unit interviewing cost per stratum.
#'       Defaults to 1. Not allowed in cluster mode, which uses stage costs.}
#'     \item{`max_weight`}{Maximum allowed gross sampling weight `N / n` per
#'       stratum, `NA` for none.}
#'     \item{`min_n_stratum`}{Minimum sample size per stratum, `NA` for none.}
#'     \item{`n_step`}{Whole step size per stratum (`NA` or 1 for none) so
#'       `n_int` is a multiple of `n_step`. Element frames only.}
#'     \item{`take_all`}{Logical (or 0/1). If `TRUE`, every unit in the stratum
#'       is selected (census stratum). Element frames only.}
#'     \item{`deff`, `resp_rate`, `eligibility_rate`, `eligible_per_unit`}{
#'       Per-stratum design effect, response rate, eligibility rate, and
#'       eligible units per sampled unit. `NA` entries fall back to the scalar
#'       argument, so a scalar other than 1 beside a column with no `NA` is
#'       refused, since it would have no effect. Vector arguments override
#'       these columns.}
#'     \item{`icc_1`, ..., `icc_{S-1}`, `var_ratio_1`, ..., `var_ratio_{S-1}`,
#'       `take_2`, ..., `take_S`, `cost_1`, ..., `cost_S`, `N_1`, ...,
#'       `N_{S-1}`, `resp_rate_1`, ..., `resp_rate_{S-1}`}{Multistage cluster
#'       columns for arbitrary stage depth S >= 2. Adding `icc_1` (and
#'       `icc_2`, ..., `icc_{S-1}` for deeper designs) activates stratified
#'       multistage allocation.}
#'   }
#'
#'   For `svyplan_prec` objects, a precision result from [prec_alloc()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param domains Character vector of column names in `frame` that partition
#'   strata into reporting domains, or `NULL` (default). When solving for `cv`,
#'   `moe`, or `rmoe`, the target is enforced inside every domain. In `n` and
#'   `budget` modes, domains are reported in `$domains` and constrained only if
#'   `min_n_domain` is supplied.
#' @param n Total sample size. Specify exactly one of `n`, `cv`, `moe`, `rmoe`,
#'   or `budget`.
#' @param cv Target coefficient of variation (relative standard error).
#'   Requires `mean` or `p` in `frame`. Specify exactly one of `n`, `cv`,
#'   `moe`, `rmoe`, or `budget`.
#' @param moe Target absolute margin of error (half-width of the confidence
#'   interval at `1 - alpha`). Specify exactly one of `n`, `cv`, `moe`, `rmoe`,
#'   or `budget`.
#' @param rmoe Target relative margin of error (`moe / abs(mean)`). Requires
#'   `mean` or `p` in `frame`. Specify exactly one of `n`, `cv`, `moe`, `rmoe`,
#'   or `budget`.
#' @param budget Total field budget. Caps the affordable allocation. Specify
#'   exactly one of `n`, `cv`, `moe`, `rmoe`, or `budget`.
#' @param alloc Allocation rule, one of `"neyman"` (default), `"optimal"`,
#'   `"proportional"`, or `"power"`.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of
#'   per-stratum unit costs overriding `frame$unit_cost`. Element frames only.
#' @param alpha Significance level for `moe` and `rmoe`, default 0.05.
#' @param deff Design effect multiplier (> 0). A scalar fills the `NA` rows
#'   of `frame$deff`, and a length-`nrow(frame)` vector overrides it.
#' @param resp_rate Expected response rate among eligible units in `(0, 1]`.
#'   Default 1. A scalar fills the `NA` rows of `frame$resp_rate`, and a
#'   length-`nrow(frame)` vector overrides it.
#' @param eligibility_rate Share of issued units in scope in `(0, 1]`.
#'   Default 1. A scalar fills the `NA` rows of `frame$eligibility_rate`, and
#'   a length-`nrow(frame)` vector overrides it.
#' @param eligible_per_unit Expected eligible units per sampled unit (> 0).
#'   Default 1. A scalar fills the `NA` rows of `frame$eligible_per_unit`,
#'   and a length-`nrow(frame)` vector overrides it. May exceed 1 (for example, adults per
#'   household). Sample counts, bounds, weights, and costs remain in sampled
#'   units while variance terms use `N * eligibility_rate * eligible_per_unit`
#'   and `resp_rate * eligibility_rate * eligible_per_unit`.
#' @param df Degrees of freedom of the variance estimator for computing `moe`
#'   and `rmoe` via the Student t quantile. `NULL` (default) uses the normal
#'   quantile.
#' @param min_n_stratum Optional minimum sample size per stratum. Scalar or
#'   length-`nrow(frame)` vector overriding `frame$min_n_stratum`.
#' @param min_n_domain Optional floor on each domain's total sample. Either a
#'   single positive number applied to every domain, or a data frame with
#'   columns `domain`, `level`, and `min_n`. Requires a single `domains`
#'   column. A floor above what the domain's strata can hold, continuously or
#'   in whole units, is refused.
#' @param alloc_q Bankier power allocation exponent in `[0, 1]`, used only when
#'   `alloc = "power"`. Default 0.5.
#' @param fpc Finite population correction mode, one of `"unit"` (default),
#'   `"stage"`, or `"none"`. `"unit"` applies `1 - n / N` to the stratum
#'   variance. `"stage"` applies the stage-by-stage correction in cluster mode
#'   and requires `N_1`, ..., `N_{S-1}`. `"none"` drops the finite population
#'   correction.
#' @param control `NULL` (default) or a list of branch-and-bound settings
#'   (`max_nodes`, `gap`) for searched whole-unit modes (`cv`/`moe`/`rmoe`
#'   under `"optimal"` or `"neyman"`, and `"optimal"` `budget` mode without
#'   domains). `max_nodes` bounds the combination search. With three or more
#'   stages, the listing of take vectors before it has its own limit.
#'   An unfinished listing can return a `"heuristic"` design.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_alloc` object with `type = "alloc"`, `stages` (1 for
#'   element frames, S >= 2 for cluster frames), top-level continuous totals
#'   and precision (`n`, `se`, `moe`, `cv`, `rmoe`), a stratum-level data
#'   frame in `$detail` (also returned by `as.data.frame()`), an `$operational`
#'   list for the integer field design, and `$domains` when `domains` is set.
#'
#'   Key columns of `$detail` include `stratum`, `N`, `sd`, `unit_cost`,
#'   continuous `n`, operational `n_int`, base `weight`, effective size
#'   `n_eff`, bounds `.lower` and `.upper`, `.binding`, `.bound_source`,
#'   stratum precision (`.se`, `.moe`, `.rmoe`, `.cv`), and variance share
#'   `.share`. Cluster mode adds `take_2`, ..., `take_S`, `n_1`, `n_1_int`,
#'   and `take_2_int`, ..., `take_S_int`. When stage populations are
#'   supplied it also adds `N_1`, ..., `N_{S-1}` and `.frac_1`, the continuous
#'   first-stage sampling fraction `n_1 / N_1`.
#'
#'   Top-level fields report the continuous allocation. For fieldwork, use
#'   `$detail$n_int` (and `$detail$n_1_int`, `$detail$take_2_int`, ... in
#'   cluster mode). `as.integer()` returns the operational integer total.
#'
#'   `$rates` gives the rates the allocation used, one row per stratum, after
#'   the argument and `frame` columns are resolved: `resp_rate_1`, ...,
#'   `resp_rate_{S-1}` in cluster mode, then `resp_rate` and
#'   `eligibility_rate`. `eligible_per_unit` is a yield and is not listed.
#'
#' @details
#' ## Allocation rules
#'
#' Allocation is controlled by `alloc` (matching [strata_bound()]).
#' \itemize{
#'   \item **proportional**: \eqn{n_h \propto N_h / r_h}{n_h ~ N_h / r_h}
#'   \item **neyman**: \eqn{n_h \propto N_h S_h \sqrt{d_h / r_h}}{n_h ~ N_h S_h sqrt(d_h / r_h)}
#'   \item **optimal**: \eqn{n_h \propto N_h S_h \sqrt{d_h / r_h} / \sqrt{c_h}}{n_h ~ N_h S_h sqrt(d_h / r_h) / sqrt(c_h)}
#'   \item **power**: Bankier (1988) allocation,
#'     \eqn{n_h \propto (S_h / |\mu_h|) X_h^q \sqrt{d_h / r_h}}{n_h ~ (S_h / |mean_h|) * alloc_measure_h^alloc_q * sqrt(deff_h / resp_rate_h)},
#'     where \eqn{X_h} is `alloc_measure` (default `N * abs(mean)`).
#' }
#'
#' Bounds (`min_n_stratum`, `max_weight`, `take_all`, `N_1`) are enforced by
#' recursive Neyman allocation (Wesolowski et al., 2022). Fixed-`n` integer
#' rounding uses bounded largest-remainder / ORIC rounding (Cont and Heidari,
#' 2014).
#'
#' ## Cluster designs within strata
#'
#' Adding `icc_1` (and `icc_2`, ..., `icc_{S-1}` for deeper designs) to `frame`
#' switches `n_alloc()` to a stratified S-stage design. Each stratum's
#' effective SD is inflated by its S-stage clustering factor and its effective
#' per-element cost is `(cost_1 + sum_{s=2}^S cost_s * prod_{j=2}^s take_j) /
#' prod_{j=2}^S take_j`. Missing `take_s` entries use the cost-optimal takes
#' from the stage costs and homogeneities, never below one expected response
#' per unit (a warning names the strata where that floor binds), and the
#' whole takes never fall below it either. A free take whose stage size per
#' parent unit cannot hold that floor is refused, in whole units as well as
#' continuously. `N_1`, ..., `N_{S-1}` cap the stage counts per stratum and
#' activate stage-level FPC when `fpc = "stage"`. Under every `fpc`, a take
#' is at most the ratio of adjacent stage sizes (`N_s / N_{s-1}`, with `N`
#' for the last stage), continuously and, rounded down, in whole units. A
#' fixed take above it is refused, and so is a total that every PSU at those
#' takes cannot hold. The ratios describe balanced stages, so a design that
#' takes more from its larger parent units needs capacity data the frame does
#' not carry. A ratio whose sizes are not both given leaves its take uncapped.
#' Every PSU counted in `N_1` is taken to contain at least one of the units
#' counted in `N`. The sampled PSUs of a stratum then hold at most `N - N_1`
#' units more than their count, under every `fpc`. A design within known
#' stage sizes meets this already. Where a stage size is missing, the bound
#' limits the takes, the totals and the whole PSU counts, free takes fall
#' when a total needs more PSUs than the room at the cost-optimal takes
#' leaves, and [prec_alloc()] warns about a supplied design that breaks it.
#' In fixed-`n` mode, when the rounded takes miss the whole total, each
#' stratum offers the whole totals near its share of `n` that its stage
#' sizes, PSU counts and bounds allow. These are other PSU counts at its
#' rounded takes, also for a stratum whose takes are all fixed, and other
#' free takes that field a nearby total exactly. Those keep the PSU count
#' and free takes within a factor of two of the continuous design, and take
#' any PSU count and takes the stage sizes allow when that design holds less
#' than one PSU. The strata totals are chosen to meet every `min_n_domain`
#' floor and to add up to `n`, or to come as close as they can. Totals far
#' from a stratum's share are not considered, so a design that reaches `n`
#' only by moving a stratum far from its share is not found, and a floor
#' that only such designs meet is refused.
#' Under `"optimal"`, a stratum whose allocation needs more than `N_1` PSUs
#' fields every PSU and raises its free takes instead.
#'
#' In the searched modes, the whole design is found by branch-and-bound over
#' whole PSU counts and whole takes at every stage count, with `control`
#' setting its limits. Its status is in `$operational$status` (`"optimal"`,
#' or `"limit"` when the node limit stopped it) and its record in
#' `$optimization$search`. With three or more stages, a stratum whose take
#' combinations are too many to list keeps a greedy rounding of the
#' continuous design. If numerical pruning loses a feasible starting design,
#' that design is retained. Both cases have status `"heuristic"` and a warning.
#' A target or budget that the cost-optimal takes cannot meet is not refused
#' when the search can still vary a take. The search then runs over every
#' take first, and the continuous design is the one at the takes of the whole
#' design, which `$params$searched_takes` records. [prec_alloc()] assesses
#' the design at those takes, and a new target given to [n_alloc()] on that
#' assessment sizes the original request again. The refusal stands when
#' the search proves that no whole design exists, and a search that stops
#' with no design raises an error of class `svyplan_plan_unknown`. With
#' `domains`, a domain whose `min_n_domain` floor sets its design, or whose
#' searched design the floor repair changes, is not proven either. Its
#' record in `$optimization$domains` reads `"heuristic"` with a `reason`,
#' and so does `$operational$status`.
#'
#' ## Domains and steps
#'
#' When `domains` is set and a precision target (`cv`, `moe`, or `rmoe`) is
#' given, every domain meets the target on its own strata. National precision
#' is reported in `$se`, `$moe`, and `$cv` without acting as a constraint. To
#' enforce national and domain targets simultaneously, or multiple indicators,
#' use [n_alloc_joint()].
#'
#' When `frame$n_step` is supplied on an element frame, every stratum's `n_int`
#' is a whole multiple of `n_step`.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer. Chapters 5
#' and 9.
#'
#' Bankier, M. D. (1988). Power allocations: determining sample sizes for
#' subnational areas. *The American Statistician*, 42(3), 174--177.
#'
#' Cont, R. and Heidari, M. (2014). Optimal rounding under integer
#' constraints. *arXiv preprint* arXiv:1501.00014.
#'
#' Weso\enc{ł}{l}owski, J., Wieczorkowski, R., and W\enc{ó}{o}jciak, W.
#' (2022). Optimality of the recursive Neyman allocation. *Journal of Survey
#' Statistics and Methodology*, 10(5), 1263--1275.
#'
#' @family stratified design functions
#' @seealso [n_alloc_joint()] for joint multi-indicator, overlapping-domain,
#'   and PSU-register allocation, [prec_alloc()] for the inverse precision
#'   assessment, [strata_bound()] for building strata, and [n_multi()] for
#'   unstratified multi-indicator sizing.
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("A", "B", "C"),
#'   N    = c(4000, 3000, 3000),
#'   sd   = c(10, 15, 8),
#'   mean = c(50, 60, 55),
#'   unit_cost = c(1, 1.5, 1)
#' )
#'
#' n_alloc(frame, n = 600)
#' n_alloc(frame, cv = 0.03)
#' n_alloc(frame, moe = 1.5)
#'
#' frame_constraints <- transform(
#'   frame,
#'   max_weight = c(25, 20, NA),
#'   take_all = c(FALSE, FALSE, TRUE)
#' )
#' n_alloc(frame_constraints, budget = 3500, alloc = "optimal", min_n_stratum = 40)
#'
#' frame_domains <- data.frame(
#'   province = c("North", "North", "South", "South"),
#'   stratum = c("Urban", "Rural", "Urban", "Rural"),
#'   N    = c(2000, 3000, 1800, 3200),
#'   sd   = c(12, 18, 10, 16),
#'   mean = c(55, 48, 58, 50)
#' )
#' n_alloc(frame_domains, domains = "province",
#'         cv = 0.04, alloc = "power", alloc_q = 0.3)
#'
#' # Power allocation on a size measure alone
#' sizes <- data.frame(stratum = c("A", "B", "C"), N = c(4e5, 1e5, 2.5e4))
#' n_alloc(sizes, n = 60, alloc = "power", alloc_q = 0.25, min_n_stratum = 15)
#'
#' # Stratified two-stage design (EAs then households per stratum)
#' frame_cluster <- data.frame(
#'   stratum = c("Urban", "Rural"),
#'   N       = c(50000, 150000),
#'   sd      = c(0.45, 0.48),
#'   mean    = c(0.35, 0.25),
#'   icc_1   = c(0.03, 0.08),
#'   cost_1  = c(300, 600),
#'   cost_2  = c(40, 60)
#' )
#' n_alloc(frame_cluster, cv = 0.05)
#'
#' @export
n_alloc <- function(frame, ...) {
  if (!missing(frame)) {
    .refuse_result_as_input(frame, "frame", "n_alloc")
    .res <- .dispatch_plan(frame, "frame", n_alloc.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_alloc")
}

#' @rdname n_alloc
#' @export
n_alloc.default <- function(
  frame,
  ...,
  domains = NULL,
  n = NULL,
  cv = NULL,
  moe = NULL,
  rmoe = NULL,
  budget = NULL,
  alloc = c("neyman", "optimal", "proportional", "power"),
  unit_cost = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  df = NULL,
  min_n_stratum = NULL,
  min_n_domain = NULL,
  alloc_q = 0.5,
  fpc = c("unit", "stage", "none"),
  control = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_alloc.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_alloc.default, c(.plan, list(...))))
  }
  .check_alloc_joint_args(list(...), "n_alloc_joint()")
  .check_unused_dots(...)
  call_args <- mget(setdiff(names(formals(n_alloc.default)),
                            c("frame", "...", "plan")), environment())
  .refuse_unnumbered_stage_columns(frame, "frame")
  .refuse_argument_beside_column(
    frame, resp_rate = resp_rate, deff = deff,
    eligibility_rate = eligibility_rate, eligible_per_unit = eligible_per_unit
  )
  min_n_stratum <- .resolve_min_n_stratum(min_n_stratum, frame)
  fpc <- match.arg(fpc)
  alloc <- match.arg(alloc)
  check_alpha(alpha)
  alloc_measure <- if ("alloc_measure" %in% names(frame)) {
    frame[["alloc_measure"]]
  }
  if (!is.null(alloc_measure) && alloc != "power") {
    stop(
      sprintf(
        "'alloc_measure' is used only by alloc = \"power\", not alloc = \"%s\". Drop the column or set alloc = \"power\"",
        alloc
      ),
      call. = FALSE
    )
  }

  if (alloc == "power") {
    if (
      !is.numeric(alloc_q) ||
        length(alloc_q) != 1L ||
        is.na(alloc_q) ||
        alloc_q < 0 ||
        alloc_q > 1
    ) {
      stop("'alloc_q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }

  mode_count <- (!is.null(n)) + (!is.null(cv)) + (!is.null(moe)) +
    (!is.null(rmoe)) + (!is.null(budget))
  if (mode_count != 1L) {
    stop(
      "specify exactly one of 'n', 'cv', 'moe', 'rmoe', or 'budget'",
      call. = FALSE
    )
  }

  if (!any(c("sd", "var", "p") %in% names(frame))) {
    if (!alloc %in% c("proportional", "power")) {
      stop(
        sprintf(
          "alloc = \"%s\" needs 'sd', 'var' or 'p' in 'frame'. Without them, only a fixed 'n' can be allocated, by alloc = \"proportional\" or \"power\"",
          alloc
        ),
        call. = FALSE
      )
    }
    if (is.null(n)) {
      stop(
        "a precision or 'budget' target needs 'sd', 'var' or 'p' in 'frame'. Without them, only a fixed 'n' can be allocated",
        call. = FALSE
      )
    }
  }
  prep <- .alloc_prepare_frame(
    frame,
    domains = domains,
    unit_cost = unit_cost,
    deff = deff,
    resp_rate = resp_rate,
    fpc = fpc,
    size_only = alloc %in% c("proportional", "power") && !is.null(n),
    eligibility_rate = eligibility_rate,
    eligible_per_unit = eligible_per_unit
  )
  N_h <- prep$N_h
  S_h <- prep$S_h
  mean_h <- prep$mean_h
  cost_h <- prep$cost_h
  deff_arg <- deff
  resp_rate_arg <- resp_rate
  deff <- prep$deff_h
  resp_rate <- prep$resp_rate_h

  target_mode <- if (!is.null(n)) {
    "n"
  } else if (!is.null(cv)) {
    "cv"
  } else if (!is.null(moe)) {
    "moe"
  } else if (!is.null(rmoe)) {
    "rmoe"
  } else {
    "budget"
  }
  mode <- if (target_mode %in% c("cv", "moe", "rmoe")) "cv" else target_mode
  target_se <- NULL
  if (target_mode == "rmoe") {
    check_scalar(rmoe, "rmoe")
    cv <- rmoe / .q_alpha(alpha, df)
  } else if (target_mode == "moe") {
    check_scalar(moe, "moe")
    target_se <- moe / .q_alpha(alpha, df)
  }
  searched <- .classic_searched(alloc, mode, length(prep$domain_idx) > 0L)
  per_domain <- .domain_searched(alloc, mode, length(prep$domain_idx) > 0L)
  ## A take-rule precheck does not stop a search over free takes, which
  ## may find a design at other takes.
  at_whole_takes <- !is.null(attr(frame, "svyplan_whole_takes"))
  takes_first <- searched && alloc == "optimal" && isTRUE(prep$cluster) &&
    (prep$stages %||% 2L) >= 3L && isTRUE(prep$has_stage_costs) &&
    any(!prep$takes_fixed_mat)
  take_rule_refusal <- NULL
  refuse <- function(msg) {
    if (!takes_first) stop(msg, call. = FALSE)
    take_rule_refusal <<- take_rule_refusal %||% msg
  }
  psu_cap <- NULL
  whole_cap <- rep(Inf, length(prep$N_h))
  if (!is.null(prep$N_psu_h)) {
    P0 <- prep$take_prod_h
    room <- ifelse(P0 > 1, (prep$N_h - prep$N_psu_h) / (P0 - 1), Inf)
    psu_cap <- pmin(prep$N_psu_h, room) * P0
    prep$psu_room_h <- room < prep$N_psu_h
    # A free take rises past N_1 under "optimal", up to its stage capacity.
    rise <- alloc == "optimal" && (searched || per_domain || mode == "n")
    rise <- rise & rowSums(!prep$takes_fixed_mat) > 0L
    whole_rows <- rise | prep$stages >= 3L
    if (any(whole_rows)) {
      cap <- .cluster_capacity(prep)
      psu_cap[rise] <- cap$continuous[rise]
      whole_cap[whole_rows] <- cap$whole[whole_rows]
      prep$takes_rise_h <- rise
    }
  }
  prep$whole_cap_h <- whole_cap
  bounds <- .alloc_bounds(
    N_h = N_h,
    max_weight = prep$max_weight,
    take_all = prep$take_all,
    min_n_stratum = min_n_stratum,
    psu_cap = psu_cap,
    step = prep$step_h,
    stratum = prep$stratum
  )
  m_h <- bounds$m_h
  M_h <- bounds$M_h
  prep$m_src <- bounds$m_src
  prep$M_src <- bounds$M_src
  lo_i <- ceiling(m_h - 1e-9)
  hi_i <- pmin(floor(M_h + 1e-9), whole_cap)
  if (any(lo_i > hi_i)) {
    bad <- prep$stratum[lo_i > hi_i]
    stop(
      sprintf(
        "no integer sample size satisfies the bounds for stratum: %s",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  lo <- sum(m_h)
  hi <- sum(M_h)
  tol <- 1e-8

  floors <- .alloc_domain_floors(min_n_domain, prep)
  alloc_fn <- function(total) .rna_alloc(a_h, total, m_h, M_h)
  if (!is.null(floors)) {
    dom <- floors$dom
    D <- length(floors$floor)
    lo_d <- vapply(seq_len(D), function(d) sum(m_h[dom == d]), numeric(1))
    hi_d <- vapply(seq_len(D), function(d) sum(M_h[dom == d]), numeric(1))
    over <- floors$floor > hi_d + tol
    if (any(over)) {
      stop(
        sprintf(
          "'min_n_domain' exceeds what domain(s) %s can hold under the stratum bounds (%s)",
          paste(sQuote(floors$levels[over]), collapse = ", "),
          paste(format(hi_d[over]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    ## A continuous cap from N_1 at the rule's takes does not bound whole designs.
    whole_h <- pmin(
      ifelse(prep$M_src %in% "N_1", floor(N_h + 1e-9), floor(M_h + 1e-9)),
      whole_cap
    )
    whole_d <- vapply(seq_len(D), function(d) sum(whole_h[dom == d]),
                      numeric(1))
    over <- floors$floor > whole_d + tol
    if (any(over)) {
      stop(
        sprintf(
          "'min_n_domain' exceeds what domain(s) %s can hold in whole units (%s)",
          paste(sQuote(floors$levels[over]), collapse = ", "),
          paste(format(whole_d[over]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    idle <- floors$floor > 0 & floors$floor <= lo_d + tol
    if (any(idle)) {
      warning(
        sprintf(
          "'min_n_domain' is at or below the stratum minimums of domain(s) %s, so it changes nothing there",
          paste(sQuote(floors$levels[idle]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    floor_d <- pmax(floors$floor, lo_d)
    lo <- sum(floor_d)
    alloc_fn <- function(total) {
      .rna_alloc_floored(a_h, total, m_h, M_h, dom, floor_d)
    }
  }

  a_h <- .alloc_weights(alloc, alloc_q, .alloc_scope_N(prep), S_h, cost_h,
                         deff, resp_rate, mean_h = mean_h,
                         alloc_measure = alloc_measure)
  if (!is.finite(sum(a_h)) || sum(a_h) <= 0) {
    a_h <- .alloc_scope_N(prep)
  }

  if (
    mode == "budget" && isTRUE(prep$cluster) && !isTRUE(prep$has_stage_costs)
  ) {
    stop(
      sprintf(
        "'budget' mode for a %d-stage cluster frame requires stage cost columns %s",
        prep$stages,
        paste(sQuote(paste0("cost_", seq_len(prep$stages))), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.null(control) && !searched && !per_domain) {
    stop(
      "'control' applies where n_alloc() searches the whole design: joint allocation, the 'cv' mode of alloc = \"optimal\" or \"neyman\" (with domains, on an element frame), and without domains the 'budget' mode of alloc = \"optimal\"",
      call. = FALSE
    )
  }
  target_total <- NA_real_
  domain_n_h <- NULL

  if (mode == "n") {
    check_scalar(n, "n")
    if (n < lo - tol) {
      stop(
        if (is.null(floors)) {
          "'n' is below the minimum feasible total under constraints"
        } else {
          sprintf("'n' is below the %s the domain floors need", format(lo))
        },
        call. = FALSE
      )
    }
    if (n > hi + tol) {
      stop(
        sprintf(
          "'n' exceeds the maximum feasible total (census bound)%s",
          .alloc_psu_bound_note(prep)
        ),
        call. = FALSE
      )
    }
    if (round(n) < sum(lo_i) || round(n) > sum(hi_i)) {
      stop(
        sprintf(
          "no integer allocation reaches total %.0f within the integer bounds (feasible totals: %.0f to %.0f)",
          round(n),
          sum(lo_i),
          sum(hi_i)
        ),
        call. = FALSE
      )
    }
    target_total <- n
  } else if (mode == "cv" && length(prep$domain_idx) > 0L) {
    if (is.null(target_se)) {
      check_scalar(cv, "cv")
      if (anyNA(mean_h)) {
        stop(
          sprintf(
            "'mean' (or 'p') is required in frame when solving for '%s'",
            target_mode
          ),
          call. = FALSE
        )
      }
      if (cv <= 0) {
        stop("'cv' must be positive", call. = FALSE)
      }
    }
    domain_n_h <- .alloc_domain_minimum(
      prep, a_h, m_h, M_h, cv, alpha, deff, resp_rate,
      floor_d = if (!is.null(floors)) floor_d, tol = tol,
      target_se = target_se, target_mode = target_mode
    )
    target_total <- sum(domain_n_h)
  } else if (mode == "cv") {
    if (is.null(target_se)) {
      check_scalar(cv, "cv")
      if (anyNA(mean_h)) {
        stop(
          sprintf(
            "'mean' (or 'p') is required in frame when solving for '%s'",
            target_mode
          ),
          call. = FALSE
        )
      }
      if (cv <= 0) {
        stop("'cv' must be positive", call. = FALSE)
      }
    }
    goal <- if (!is.null(target_se)) target_se else cv
    metric_for_total <- function(n_total) {
      n_h <- alloc_fn(n_total)
      at <- .cluster_take_rise(prep, n_h)
      m <- .alloc_metrics(
        N_h = N_h,
        N_pop = .alloc_scope_N(prep),
        S_h = at$S_h,
        mean_h = mean_h,
        n_h = n_h,
        alpha = alpha,
        deff = deff,
        resp_rate = resp_rate,
        cost_h = at$cost_h,
        N_fpc = at$N_fpc_h
      )
      if (!is.null(target_se)) m$se else m$cv
    }

    val_lo <- metric_for_total(lo)
    if (val_lo <= goal + tol) {
      target_total <- lo
    } else {
      val_hi <- metric_for_total(hi)
      if (!is.finite(val_hi) || val_hi > goal + tol) {
        refuse(sprintf(
          "target '%s' is unattainable under current constraints%s",
          target_mode,
          .alloc_psu_bound_note(prep)
        ))
      }
      target_total <- if (!is.null(take_rule_refusal) || val_hi >= goal) {
        hi
      } else {
        uniroot(
          function(x) metric_for_total(x) - goal,
          interval = c(lo, hi),
          tol = 1e-8
        )$root
      }
    }
  } else {
    check_scalar(budget, "budget")
    btol <- 1e-8 * budget
    cost_at <- function(n_h) {
      sum(n_h * .cluster_take_rise(prep, n_h)$cost_h)
    }
    cost_lo <- cost_at(if (is.null(floors)) m_h else alloc_fn(lo))
    cost_hi <- cost_at(M_h)
    if (budget < cost_at(lo_i) - btol) {
      refuse(sprintf(
        "'budget' cannot fund the integer lower bounds (minimum integer cost = %.4g)",
        cost_at(lo_i)
      ))
    }
    if (budget < cost_lo - btol) {
      refuse("'budget' is below the minimum feasible cost under constraints")
    }
    ## At the takes of a searched whole design, a budget above the census
    ## leaves money unspent.
    if (budget > cost_hi + btol && !at_whole_takes) {
      refuse(sprintf(
        "'budget' exceeds the maximum feasible cost (census bound)%s",
        .alloc_psu_bound_note(prep)
      ))
    }

    if (!is.null(take_rule_refusal) || abs(budget - cost_lo) <= btol) {
      target_total <- lo
    } else if (budget >= cost_hi) {
      target_total <- hi
    } else {
      cost_for_total <- function(n_total) {
        cost_at(alloc_fn(n_total))
      }
      target_total <- uniroot(
        function(x) cost_for_total(x) - budget,
        interval = c(lo, hi),
        tol = 1e-8
      )$root
    }
  }

  if (!is.null(take_rule_refusal)) {
    return(.alloc_takes_first(
      frame, call_args, prep, alloc, mode, cv, budget, deff, resp_rate,
      lo_i, hi_i, control, target_se, take_rule_refusal
    ))
  }
  n_h <- if (!is.null(domain_n_h)) {
    as.numeric(domain_n_h)
  } else {
    alloc_fn(target_total)
  }
  prep <- .cluster_take_rise(prep, n_h)
  metrics <- .alloc_metrics(
    N_h = N_h,
    N_pop = .alloc_scope_N(prep),
    S_h = prep$S_h,
    mean_h = mean_h,
    n_h = n_h,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = prep$cost_h,
    df = df,
    N_fpc = prep$N_fpc_h
  )

  detail <- .alloc_detail(
    prep = prep,
    n_h = n_h,
    m_h = m_h,
    M_h = M_h,
    resp_rate = resp_rate,
    deff = deff,
    mode = mode,
    budget = budget,
    lo_i = lo_i,
    hi_i = hi_i,
    metrics = metrics
  )
  domain_summary <- .alloc_domain_summary(
    prep = prep,
    n_h = n_h,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  params <- list(
    frame = frame,
    alloc = alloc,
    alpha = alpha,
    fpc = fpc,
    deff = deff_arg,
    resp_rate = resp_rate_arg,
    cost_h = cost_h,
    min_n_stratum = min_n_stratum,
    df = df,
    domain_cols = prep$domain_cols,
    domain_idx = prep$domain_idx,
    mode = target_mode,
    alloc_q = if (alloc == "power") alloc_q else NULL,
    n_h = n_h,
    achieved = list(n = sum(n_h), cv = metrics$cv, cost = metrics$cost)
  )
  if (!identical(eligibility_rate, 1)) params$eligibility_rate <- eligibility_rate
  if (!identical(eligible_per_unit, 1)) params$eligible_per_unit <- eligible_per_unit
  params[[target_mode]] <- switch(
    target_mode,
    n = n,
    cv = cv,
    moe = moe,
    rmoe = rmoe,
    budget = budget
  )
  if (!is.null(floors)) {
    params$min_n_domain <- min_n_domain
    total_d <- vapply(seq_along(floor_d), function(d) sum(n_h[dom == d]),
                      numeric(1))
    binds <- floors$floor > lo_d + tol & abs(total_d - floor_d) <= 1e-6 *
      pmax(1, floor_d)
    raised <- binds[dom] & !detail$.binding
    detail$.binding[raised] <- TRUE
    detail$.bound_source[raised] <- "min_n_domain"
    if (!is.null(domain_summary)) {
      domain_summary$.min_n <- ifelse(floors$floor > 0, floors$floor,
                                      NA_real_)
      domain_summary$.floor_binding <- binds
    }
  }

  whole <- NULL
  if (isTRUE(prep$cluster)) {
    prep$given_whole <- attr(frame, "svyplan_whole_takes")
    opc <- .alloc_operational_cluster(
      prep,
      detail,
      n_h,
      mode,
      budget,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      target_cv = if (mode == "cv" && is.null(target_se)) cv else NULL,
      target_se = if (mode == "cv") target_se else NULL,
      df = df,
      alloc_weights = a_h,
      alloc = alloc,
      max_nodes = if (searched || per_domain) {
        .search_control(control)$max_nodes
      },
      gap = .search_control(control)$gap,
      floors = if (!is.null(floors)) {
        list(dom = dom, need = floor_d, levels = floors$levels)
      },
      domain_search = if (per_domain) {
        list(floored = attr(domain_n_h, "floored"))
      }
    )
    detail <- opc$detail
    operational <- opc$operational
    if (searched || per_domain) {
      whole <- list(search = opc$search, domains = opc$domains)
      params$control <- control
    }
  } else {
    if (searched) {
      whole <- .alloc_whole_element(
        prep, alloc, mode, cv, budget, deff, resp_rate, lo_i, hi_i,
        max_nodes = .search_control(control)$max_nodes,
        gap = .search_control(control)$gap,
        target_se = target_se
      )
      detail$n_int <- whole$n_int
      params$control <- control
    } else if (per_domain) {
      # A domain its floor holds above the target keeps the rounded design,
      # which the floor repair below raises to the floor.
      floored <- attr(domain_n_h, "floored")
      records <- list()
      for (d in seq_along(prep$domain_idx)) {
        nm <- names(prep$domain_idx)[d]
        if (floored[d]) {
          records[[nm]] <- list(status = "heuristic", reason = "floor")
          next
        }
        idx <- prep$domain_idx[[d]]
        w <- .alloc_whole_element(
          .alloc_prep_slice(prep, idx), alloc, mode, cv, budget,
          .subset_h(deff, idx), .subset_h(resp_rate, idx), lo_i[idx],
          hi_i[idx],
          max_nodes = .search_control(control)$max_nodes,
          gap = .search_control(control)$gap,
          target_se = target_se
        )
        detail$n_int[idx] <- w$n_int
        records[[nm]] <- w$search %||% list(status = w$status)
      }
      whole <- list(domains = records)
      params$control <- control
    }
    if (!is.null(floors)) {
      # in whole steps, so the repair keeps every n_int on its step
      step <- prep$step_h %||% rep(1L, length(n_h))
      total <- sum(detail$n_int)
      a_lo <- ceiling(lo_i / step - 1e-9)
      a_hi <- floor(hi_i / step + 1e-9)
      need <- ceiling(floor_d - 1e-9)
      a <- .alloc_floor_repair(
        round(detail$n_int / step), cont = n_h / step, unit = step,
        lo = a_lo, hi = a_hi, dom = dom, need = need,
        keep = switch(mode, n = "total", budget = "budget", "none"),
        total = total, cost_step = prep$cost_h * step, budget = budget
      )
      if (identical(mode, "n") && sum(a * step) != total) {
        # whole steps larger than the excess cannot give it back one at a
        # time, so the total is reached by moving several strata together
        closed <- .floor_close_total(a, step, a_lo, a_hi, dom, need, n_h,
                                     total)
        if (length(closed) == 0L) {
          reached <- attr(closed, "totals", exact = TRUE)
          stop(
            sprintf(
              "no allocation in steps of 'n_step' meets every 'min_n_domain' floor at the total %d. Within one step of each stratum's floored design, the nearest totals are %s",
              as.integer(total),
              paste(c(
                if (any(reached < total)) max(reached[reached < total]),
                if (any(reached > total)) min(reached[reached > total])
              ), collapse = " and ")
            ),
            call. = FALSE
          )
        }
        a <- closed
      }
      repaired <- as.integer(a * step)
      if (!is.null(whole$domains)) {
        whole$domains <- .whole_repaired(whole$domains, prep$domain_idx,
                                         detail$n_int, repaired)
      }
      detail$n_int <- repaired
    }
    if (per_domain && !searched) {
      whole <- .whole_domains_summary(whole$domains)
    }
    operational <- .alloc_operational_element(
      prep,
      detail,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      df = df
    )
    operational$status <- whole$status
  }

  carried <- .alloc_carry_labels(detail, frame)
  detail <- carried$detail
  params$strata_keys <- carried$keys

  obj <- .new_svyplan_alloc(
    n = sum(n_h),
    method = alloc,
    stages = prep$stages %||% 1L,
    params = params,
    se = metrics$se,
    moe = metrics$moe,
    cv = metrics$cv,
    rmoe = metrics$rmoe,
    detail = detail,
    rates = .alloc_rates_table(prep),
    binding = .alloc_binding_label(n_h, m_h, M_h),
    domains = domain_summary,
    operational = operational,
    optimization = if (!is.null(whole)) {
      list(search = whole$search, domains = whole$domains)
    }
  )

  obj
}

#' @rdname n_alloc
#' @export
n_alloc.svyplan_prec <- function(
  frame,
  ...,
  n = NULL,
  cv = NULL,
  moe = NULL,
  rmoe = NULL,
  budget = NULL
) {
  x <- frame
  if (identical(x$type, "alloc_joint")) {
    stop(
      "use n_alloc_joint() for a svyplan_prec of type 'alloc_joint', from prec_alloc_joint()",
      call. = FALSE
    )
  }
  if (!identical(x$type, "alloc")) {
    stop("n_alloc requires a svyplan_prec of type 'alloc'", call. = FALSE)
  }
  p <- x$params

  frame <- p$frame
  if (is.null(n) && is.null(cv) && is.null(moe) && is.null(rmoe) &&
      is.null(budget)) {
    n <- p$achieved$n
    ## Resizing the same design keeps the takes its search chose.
    frame <- .alloc_frame_at_takes(frame, p$searched_takes)
  }

  args <- list(
    frame = frame,
    domains = p$domain_cols,
    n = n,
    cv = cv,
    moe = moe,
    rmoe = rmoe,
    budget = budget,
    alloc = p$alloc %||% "neyman",
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    df = p$df,
    min_n_stratum = p$min_n_stratum,
    min_n_domain = p[["min_n_domain"]],
    alloc_q = p$alloc_q %||% 0.5
  )
  if (!.alloc_is_cluster(p$frame)) {
    args$unit_cost <- p$cost_h
  }
  args$fpc <- p$fpc %||% "unit"
  rebuilt <- if (!is.null(n)) {
    "n"
  } else if (!is.null(cv) || !is.null(moe) || !is.null(rmoe)) {
    "cv"
  } else {
    "budget"
  }
  if (.classic_searched(args$alloc, rebuilt, length(p$domain_cols) > 0L) ||
      .domain_searched(args$alloc, rebuilt, length(p$domain_cols) > 0L)) {
    args$control <- p[["control"]]
  }
  do.call(n_alloc.default, .roundtrip_args(args, list(...), n_alloc.default))
}

#' @keywords internal
#' @noRd
.alloc_prepare_frame <- function(
  frame,
  domains = NULL,
  unit_cost = NULL,
  deff = 1,
  resp_rate = 1,
  fpc = "unit",
  size_only = FALSE,
  eligibility_rate = 1,
  eligible_per_unit = 1
) {
  if (!is.data.frame(frame) || nrow(frame) == 0L) {
    stop("'frame' must be a non-empty data frame", call. = FALSE)
  }

  if (!"N" %in% names(frame)) {
    stop("'frame' must contain an 'N' column", call. = FALSE)
  }
  N_h <- frame$N
  if (!is.numeric(N_h) || anyNA(N_h) || any(!is.finite(N_h)) || any(N_h <= 0)) {
    stop("'N' must contain positive finite values", call. = FALSE)
  }
  e_h <- .alloc_resolve_h(
    eligibility_rate,
    frame,
    "eligibility_rate",
    nrow(frame),
    .check_resp_rate_h
  )
  elig_h <- .alloc_resolve_h(
    eligible_per_unit,
    frame,
    "eligible_per_unit",
    nrow(frame),
    .check_deff_h
  )
  yield_h <- e_h * elig_h

  if ("sd" %in% names(frame) && "var" %in% names(frame)) {
    stop("'frame' must contain either 'sd' or 'var', not both", call. = FALSE)
  }
  if ("mean" %in% names(frame) && "p" %in% names(frame)) {
    stop("'frame' must contain either 'mean' or 'p', not both", call. = FALSE)
  }

  if ("sd" %in% names(frame)) {
    S_h <- frame$sd
    if (
      !is.numeric(S_h) || anyNA(S_h) || any(!is.finite(S_h)) || any(S_h < 0)
    ) {
      stop("'sd' must contain non-negative finite values", call. = FALSE)
    }
  } else if ("var" %in% names(frame)) {
    vv <- frame$var
    if (!is.numeric(vv) || anyNA(vv) || any(!is.finite(vv)) || any(vv < 0)) {
      stop("'var' must contain non-negative finite values", call. = FALSE)
    }
    S_h <- sqrt(vv)
  } else if ("p" %in% names(frame)) {
    S_h <- .alloc_sd_from_p(.alloc_check_p(frame$p), N_h, yield_h, frame, fpc)
  } else if (size_only) {
    S_h <- rep(NA_real_, nrow(frame))
  } else {
    stop("'frame' must contain 'sd', 'var' or 'p'", call. = FALSE)
  }
  size_only <- all(is.na(S_h))
  if (!size_only && all(S_h == 0)) {
    warning(
      "all 'sd' values are zero, allocation has no variability to distribute",
      call. = FALSE
    )
  }

  mean_h <- rep(NA_real_, nrow(frame))
  if ("mean" %in% names(frame)) {
    mean_h <- frame$mean
  } else if ("p" %in% names(frame)) {
    mean_h <- .alloc_check_p(frame$p)
  }
  if (any(!is.na(mean_h) & !is.finite(mean_h))) {
    stop("'mean' (or 'p') must contain finite values", call. = FALSE)
  }
  if (!all(is.na(mean_h)) && all(mean_h == 0, na.rm = TRUE)) {
    warning(
      "all 'mean' (or 'p') values are zero, CV will be Inf",
      call. = FALSE
    )
  }

  if ("cost" %in% names(frame)) {
    stop(
      "the per-stratum cost column is 'unit_cost', 'cost' is the total field cost",
      call. = FALSE
    )
  }

  if (!is.null(unit_cost)) {
    if (
      !is.numeric(unit_cost) ||
        anyNA(unit_cost) ||
        any(!is.finite(unit_cost)) ||
        any(unit_cost <= 0)
    ) {
      stop("'unit_cost' must contain positive finite values", call. = FALSE)
    }
    if (length(unit_cost) == 1L) {
      cost_h <- rep(unit_cost, nrow(frame))
    } else if (length(unit_cost) == nrow(frame)) {
      cost_h <- unit_cost
    } else {
      stop("'unit_cost' must have length 1 or nrow(frame)", call. = FALSE)
    }
  } else if ("unit_cost" %in% names(frame)) {
    cost_h <- frame$unit_cost
    if (
      !is.numeric(cost_h) ||
        anyNA(cost_h) ||
        any(!is.finite(cost_h)) ||
        any(cost_h <= 0)
    ) {
      stop("'unit_cost' must contain positive finite values", call. = FALSE)
    }
  } else {
    cost_h <- rep(1, nrow(frame))
  }

  deff_h <- .alloc_resolve_h(deff, frame, "deff", nrow(frame), .check_deff_h)
  resp_given_h <- .alloc_resolve_h(
    resp_rate,
    frame,
    "resp_rate",
    nrow(frame),
    .check_resp_rate_h
  )
  resp_rate_h <- resp_given_h * yield_h

  max_weight <- rep(NA_real_, nrow(frame))
  if ("max_weight" %in% names(frame)) {
    max_weight <- frame$max_weight
    bad <- !is.na(max_weight) & (!is.finite(max_weight) | max_weight < 1)
    if (any(bad)) {
      stop("'max_weight' must be >= 1 and finite when provided", call. = FALSE)
    }
  }

  take_all <- .check_take_all(frame[["take_all"]], nrow(frame))

  stratum <- if ("stratum" %in% names(frame)) {
    as.character(frame$stratum)
  } else {
    as.character(seq_len(nrow(frame)))
  }

  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(frame))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in frame: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  domain_cols <- domains %||% character(0)
  domain_idx <- list()
  domain_values <- NULL
  if (length(domain_cols) > 0L) {
    key <- .domain_key(frame, domain_cols)
    lev <- unique(key)
    domain_idx <- setNames(lapply(lev, function(l) which(key == l)), lev)
    domain_values <- frame[match(lev, key), domain_cols, drop = FALSE]
    rownames(domain_values) <- NULL

    for (i in seq_along(domain_idx)) {
      s_d <- stratum[domain_idx[[i]]]
      if (anyDuplicated(s_d)) {
        lab <- paste(
          unlist(lapply(domain_values[i, , drop = FALSE], as.character)),
          collapse = ":"
        )
        stop(
          sprintf(
            "duplicate stratum labels in domain '%s': %s",
            lab,
            paste(s_d[duplicated(s_d)], collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  } else {
    if (anyDuplicated(stratum)) {
      dups <- stratum[duplicated(stratum)]
      stop(
        sprintf(
          "duplicate stratum labels: %s",
          paste(unique(dups), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  step_h <- .alloc_resolve_step(frame, take_all)

  prep <- list(
    frame = frame,
    size_only = size_only,
    step_h = step_h,
    N_h = as.numeric(N_h),
    Ne_h = as.numeric(N_h * yield_h),
    eligibility_h = e_h,
    eligible_per_unit_h = elig_h,
    S_h = as.numeric(S_h),
    mean_h = as.numeric(mean_h),
    cost_h = as.numeric(cost_h),
    deff_h = deff_h,
    resp_rate_h = resp_rate_h,
    resp_given_h = resp_given_h,
    max_weight = as.numeric(max_weight),
    take_all = as.logical(take_all),
    stratum = stratum,
    domain_cols = domain_cols,
    domain_idx = domain_idx,
    domain_values = domain_values
  )
  .alloc_cluster_prep(prep, frame, unit_cost, fpc)
}

#' Resolve a frame's per-stratum step
#'
#' `NA` and 1 mean no step. A census stratum keeps its whole population, so
#' its step is 1 whatever the column says.
#' @keywords internal
#' @noRd
.alloc_resolve_step <- function(frame, take_all) {
  H <- nrow(frame)
  if (!"n_step" %in% names(frame)) {
    return(rep(1L, H))
  }
  if (.alloc_is_cluster(frame)) {
    stop(
      "'n_step' applies to element frames, where it counts the units 'n' counts. A cluster frame does not take it. To allocate PSU counts in steps, allocate them on a frame without 'icc_1', with 'n' counting PSUs",
      call. = FALSE
    )
  }
  k <- frame[["n_step"]]
  if (is.logical(k) && all(is.na(k))) {
    k <- rep(NA_real_, H)
  }
  if (!is.numeric(k) || any(!is.na(k) & (!is.finite(k) | k < 1 | k != trunc(k)))) {
    stop("'n_step' must contain whole numbers >= 1, or NA for no step",
         call. = FALSE)
  }
  k[is.na(k) | take_all] <- 1
  as.integer(k)
}

#' Refuse a step where the allocation rounds on its own
#' @keywords internal
#' @noRd
.alloc_refuse_step <- function(frame, where) {
  if ("n_step" %in% names(frame)) {
    stop(
      sprintf(
        "'n_step' is not used for %s, which builds its own whole design. It applies to the classic modes of n_alloc() on an element frame",
        where
      ),
      call. = FALSE
    )
  }
}

#' Validate a frame's stratum proportions
#' @keywords internal
#' @noRd
.alloc_check_p <- function(ph) {
  if (!is.numeric(ph) || anyNA(ph) || any(!is.finite(ph))) {
    stop("'p' must contain finite numeric values", call. = FALSE)
  }
  if (any(ph < 0 | ph > 1)) {
    stop("'p' must contain values in [0, 1]", call. = FALSE)
  }
  ph
}

#' Carry a frame's label columns into an allocation's stratum detail
#' @keywords internal
#' @noRd
.alloc_carry_labels <- function(detail, frame) {
  is_label <- vapply(
    frame,
    function(v) is.character(v) || is.factor(v) || is.logical(v),
    logical(1)
  )
  labels <- setdiff(names(frame)[is_label], c("stratum", "take_all", "n_step"))
  if (length(labels) == 0L) return(list(detail = detail, keys = NULL))
  clash <- intersect(labels, names(detail))
  if (length(clash) > 0L) {
    stop(
      sprintf(
        "'frame' label column%s %s share%s a name with the allocation's own columns. Rename %s",
        if (length(clash) == 1L) "" else "s",
        paste(sQuote(clash), collapse = ", "),
        if (length(clash) == 1L) "s" else "",
        if (length(clash) == 1L) "it" else "them"
      ),
      call. = FALSE
    )
  }
  carried <- as.data.frame(frame)[labels]
  rownames(carried) <- NULL
  at <- match("stratum", names(detail))
  detail <- cbind(
    detail[seq_len(at)], carried, detail[-seq_len(at)]
  )
  keys <- if (anyDuplicated(carried) == 0L) labels else NULL
  list(detail = detail, keys = keys)
}

#' Whether a classic allocation was made from sizes alone, with no variance
#' @keywords internal
#' @noRd
.alloc_size_only <- function(x) {
  !identical(x$method, "bethel") &&
    !any(c("sd", "var", "p") %in% names(x$params$frame))
}

#' Detect the number of sampling stages from a frame's numbered `icc_s` columns
#' @keywords internal
#' @noRd
.alloc_detect_stages <- function(frame) {
  nms <- names(frame)
  icc_cols <- grep("^icc_[0-9]+$", nms, value = TRUE)
  if (length(icc_cols) == 0L) {
    return(1L)
  }
  idx <- sort(as.integer(sub("^icc_", "", icc_cols)))
  if (!identical(idx, seq_along(idx))) {
    stop(
      sprintf(
        "stage ICC columns must be numbered consecutively from 1 ('icc_1' .. 'icc_%d')",
        max(idx)
      ),
      call. = FALSE
    )
  }
  length(idx) + 1L
}

#' Multi-stage-within-strata transformation (arbitrary S >= 2 stages)
#'
#' When the frame carries numbered `icc_1` (and `icc_2`, ..., `icc_{S-1}`)
#' columns, the stratified S-stage variance model reduces to the element model
#' with an inflated stratum SD, S_h * sqrt(inflation_h), an effective FPC
#' population N_fpc_h, and (with stage costs) an effective per-element cost.
#' Downstream allocation, constraints, and metrics then apply unchanged.
#' @keywords internal
#' @noRd
.alloc_cluster_prep <- function(prep, frame, unit_cost, fpc = "unit") {
  nms <- names(frame)
  S <- .alloc_detect_stages(frame)
  if (S == 1L) {
    if (identical(fpc, "stage")) {
      stop(
        "'fpc = \"stage\"' applies to cluster allocation, where the correction has stages to choose between. An element allocation takes 'fpc = \"unit\"' for the correction 1 - n/N or 'fpc = \"none\"' for no correction. Add 'icc_1', or choose one of those",
        call. = FALSE
      )
    }
    orphan <- grep(
      "^(var_ratio|take|cost|N|resp_rate)_[0-9]+$",
      nms,
      value = TRUE
    )
    if (length(orphan) > 0L) {
      stop(
        sprintf(
          "cluster column(s) %s require an 'icc_1' column",
          paste(sQuote(orphan), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    prep$fpc <- fpc
    if (identical(fpc, "none")) {
      prep$N_fpc_h <- rep(Inf, nrow(frame))
    }
    prep$cluster <- FALSE
    prep$stages <- 1L
    return(prep)
  }

  # Reject numbered stage columns outside the detected S-stage range
  bad_vr <- setdiff(
    grep("^var_ratio_[0-9]+$", nms, value = TRUE),
    paste0("var_ratio_", seq_len(S - 1L))
  )
  bad_take <- setdiff(
    grep("^take_[0-9]+$", nms, value = TRUE),
    paste0("take_", 2L:S)
  )
  bad_cost <- setdiff(
    grep("^cost_[0-9]+$", nms, value = TRUE),
    paste0("cost_", seq_len(S))
  )
  bad_N <- setdiff(
    grep("^N_[0-9]+$", nms, value = TRUE),
    paste0("N_", seq_len(S - 1L))
  )
  bad_rr <- setdiff(
    grep("^resp_rate_[0-9]+$", nms, value = TRUE),
    paste0("resp_rate_", seq_len(S - 1L))
  )
  bad_cols <- c(bad_vr, bad_take, bad_cost, bad_N, bad_rr)
  if (length(bad_cols) > 0L) {
    stop(
      sprintf(
        "stage column(s) %s do not match a %d-stage design (detected from 'icc_1'%s)",
        paste(sQuote(bad_cols), collapse = ", "),
        S,
        if (S > 2L) sprintf(" .. 'icc_%d'", S - 1L) else ""
      ),
      call. = FALSE
    )
  }

  if ("unit_cost" %in% nms || !is.null(unit_cost)) {
    stop(
      sprintf(
        "with 'icc_1' in the frame, use 'cost_1' .. 'cost_%d' columns instead of 'unit_cost'",
        S
      ),
      call. = FALSE
    )
  }

  H <- nrow(frame)
  icc_mat <- matrix(0, nrow = H, ncol = S - 1L)
  for (s in seq_len(S - 1L)) {
    col <- paste0("icc_", s)
    v <- frame[[col]]
    if (!is.numeric(v) || anyNA(v) || any(!is.finite(v))) {
      stop(sprintf("'%s' must contain finite numeric values", col), call. = FALSE)
    }
    .check_cluster_icc_open(v, context = "n_alloc()")
    icc_mat[, s] <- as.numeric(v)
  }

  vr_mat <- matrix(0, nrow = H, ncol = S - 1L)
  for (s in seq_len(S - 1L)) {
    col <- paste0("var_ratio_", s)
    if (col %in% nms) {
      v <- frame[[col]]
      if (!is.numeric(v) || anyNA(v) || any(!is.finite(v)) || any(v <= 0)) {
        stop(
          sprintf("'%s' must contain positive finite values", col),
          call. = FALSE
        )
      }
      vr_mat[, s] <- as.numeric(v)
    } else if (s == 1L) {
      vr_mat[, 1L] <- 1
    } else {
      vr_mat[, s] <- vr_mat[, s - 1L] * (1 - icc_mat[, s - 1L])
    }
  }

  cost_cols <- paste0("cost_", seq_len(S))
  has_costs <- any(cost_cols %in% nms)
  cost_mat <- NULL
  if (has_costs) {
    if (!all(cost_cols %in% nms)) {
      stop(
        sprintf(
          "%s must be supplied together",
          paste(sQuote(cost_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    cost_mat <- matrix(0, nrow = H, ncol = S)
    for (s in seq_len(S)) {
      col <- cost_cols[s]
      cc <- frame[[col]]
      if (!is.numeric(cc) || anyNA(cc) || any(!is.finite(cc)) || any(cc <= 0)) {
        stop(
          sprintf("'%s' must contain positive finite values", col),
          call. = FALSE
        )
      }
      cost_mat[, s] <- as.numeric(cc)
    }
  }

  takes_mat <- matrix(NA_real_, nrow = H, ncol = S - 1L)
  for (s in 2L:S) {
    col <- paste0("take_", s)
    if (col %in% nms) {
      ps <- frame[[col]]
      if (is.logical(ps) && all(is.na(ps))) {
        ps <- rep(NA_real_, H)
      }
      if (!is.numeric(ps) || any(!is.na(ps) & (!is.finite(ps) | ps < 1))) {
        stop(
          sprintf("'%s' must contain values >= 1 (NA for cost-optimal)", col),
          call. = FALSE
        )
      }
      takes_mat[, s - 1L] <- as.numeric(ps)
    }
  }
  takes_fixed_mat <- !is.na(takes_mat)

  N_cols <- paste0("N_", seq_len(S - 1L))
  present_N <- N_cols[N_cols %in% nms]
  N_stages_mat <- NULL
  N_1_h <- NULL
  if (length(present_N) > 0L) {
    if (!"N_1" %in% present_N) {
      stop(
        sprintf(
          "stage population columns %s require 'N_1'",
          paste(sQuote(present_N), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    N_stages_mat <- matrix(NA_real_, nrow = H, ncol = S - 1L)
    for (s in seq_len(S - 1L)) {
      col <- N_cols[s]
      if (!col %in% nms) next
      nv <- frame[[col]]
      if (
        !is.numeric(nv) ||
          anyNA(nv) ||
          any(!is.finite(nv)) ||
          any(nv < 1) ||
          any(abs(nv - round(nv)) > 1e-8)
      ) {
        stop(
          sprintf("'%s' must contain positive whole numbers, one per stratum", col),
          call. = FALSE
        )
      }
      N_stages_mat[, s] <- as.numeric(round(nv))
    }
    N_1_h <- N_stages_mat[, 1L]
    # Validate monotone hierarchy N_1 <= N_2 <= ... <= N_{S-1} <= N
    prev <- N_stages_mat[, 1L]
    prev_name <- "N_1"
    for (s in seq_len(S - 1L)) {
      if (s > 1L && !anyNA(N_stages_mat[, s])) {
        if (any(prev > N_stages_mat[, s] + 1e-8)) {
          stop(
            sprintf("'%s' must not exceed 'N_%d'", prev_name, s),
            call. = FALSE
          )
        }
        prev <- N_stages_mat[, s]
        prev_name <- paste0("N_", s)
      }
    }
    if (any(prev > prep$N_h + 1e-8)) {
      stop(
        sprintf("'%s' must not exceed the ultimate-unit population 'N'", prev_name),
        call. = FALSE
      )
    }
  }
  prep$N_stages_mat <- N_stages_mat
  prep$N_psu_h <- N_1_h

  if (identical(fpc, "stage")) {
    missing_N <- setdiff(N_cols, nms)
    if (length(missing_N) > 0L) {
      stop(
        sprintf(
          "fpc = \"stage\" needs stage population column(s) %s, or use the default fpc = \"unit\"",
          paste(sQuote(missing_N), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  stage_resp_mat <- matrix(1, nrow = H, ncol = S)
  for (s in seq_len(S - 1L)) {
    stage_resp_mat[, s] <- .alloc_resolve_h(
      1,
      frame,
      paste0("resp_rate_", s),
      H,
      .check_resp_rate_h
    )
  }
  stage_resp_mat[, S] <- prep$resp_rate_h

  if (!is.null(N_stages_mat)) {
    for (h in which(apply(takes_fixed_mat, 1L, any))) {
      caps <- .stage_caps_of(prep$N_h[h], N_stages_mat, h, S)
      t <- takes_mat[h, ]
      over <- which(takes_fixed_mat[h, ] & ifelse(
        t == round(t), t > .whole_caps(caps), t > caps * (1 + 1e-12)
      ))
      if (length(over)) {
        s <- over[1L]
        stop(
          sprintf(
            "take_%d = %s in stratum '%s' exceeds the %s units per parent unit at that stage",
            s + 1L, format(takes_mat[h, s], digits = 15, scientific = FALSE),
            prep$stratum[h], format(caps[s], digits = 15, scientific = FALSE)
          ),
          call. = FALSE
        )
      }
    }
  }

  # Solve optimal free takes where NA
  need_opt_row <- apply(!takes_fixed_mat, 1L, any)
  if (any(need_opt_row)) {
    if (!has_costs) {
      stop(
        sprintf(
          "%s are required when %s %s not fixed for every stratum",
          paste(sQuote(cost_cols), collapse = ", "),
          paste(sQuote(paste0("take_", 2L:S)), collapse = ", "),
          if (S == 2L) "is" else "are"
        ),
        call. = FALSE
      )
    }
    .cluster_take_room(prep$stratum, prep$N_h, N_stages_mat,
                       stage_resp_mat[, -1L, drop = FALSE], takes_fixed_mat,
                       need_opt_row, fpc)
    clamped <- matrix(FALSE, nrow = H, ncol = S - 1L)
    if (S == 2L) {
      open_take <- .cluster_free_take(
        cost_mat[need_opt_row, 1L],
        cost_mat[need_opt_row, 2L],
        icc_mat[need_opt_row, 1L],
        stage_resp_mat[need_opt_row, 2L],
        M = if (identical(fpc, "stage")) {
          .alloc_scope_N(prep)[need_opt_row] / N_1_h[need_opt_row]
        },
        cap = if (!is.null(N_1_h)) prep$N_h[need_opt_row] / N_1_h[need_opt_row]
      )
      floor_2 <- vapply(stage_resp_mat[need_opt_row, 2L], .take_bound, 1)
      clamped[need_opt_row, 1L] <- open_take < floor_2 - 1e-9
      takes_mat[need_opt_row, 1L] <- pmax(open_take, floor_2)
    } else {
      for (h in which(need_opt_row)) {
        R_h <- matrix(0, nrow = 1L, ncol = S)
        for (s in seq_len(S - 1L)) {
          R_h[1L, s] <- vr_mat[h, s] * icc_mat[h, s]
        }
        R_h[1L, S] <- vr_mat[h, S - 1L] * (1 - icc_mat[h, S - 1L])
        rates_h <- matrix(stage_resp_mat[h, ], nrow = 1L)
        M_h <- matrix(Inf, nrow = 1L, ncol = S)
        if (identical(fpc, "stage")) {
          M_h[1L, 1L] <- N_stages_mat[h, 1L]
          for (s in 2L:(S - 1L)) {
            M_h[1L, s] <- N_stages_mat[h, s] / N_stages_mat[h, s - 1L]
          }
          M_h[1L, S] <- .alloc_scope_N(prep)[h] / N_stages_mat[h, S - 1L]
        }
        alpha_h <- .stage_telescoping(R_h, rates_h, M_h)$alpha[1L, ]
        fixed_h <- takes_fixed_mat[h, ]
        floor_h <- vapply(stage_resp_mat[h, -1L], .take_bound, 1)
        lo_h <- ifelse(fixed_h, takes_mat[h, ], floor_h)
        cap_h <- .stage_caps_of(prep$N_h[h], N_stages_mat, h, S)
        hi_h <- ifelse(
          fixed_h,
          takes_mat[h, ],
          ifelse(is.finite(cap_h), cap_h, Inf)
        )
        takes_mat[h, ] <- .stage_min_gh(alpha_h, cost_mat[h, ], lo_h, hi_h)$takes
        open_h <- .stage_min_gh(
          alpha_h, cost_mat[h, ], ifelse(fixed_h, takes_mat[h, ], 1e-9), hi_h
        )$takes
        clamped[h, ] <- !fixed_h & open_h < floor_h - 1e-9
      }
    }
    for (s in which(colSums(clamped) > 0L)) {
      warning(
        sprintf(
          "cost-optimal take_%d is below 1 expected response (high 'icc' relative to stage costs) in %s; clamped to one expected response",
          s + 1L,
          paste(sQuote(prep$stratum[clamped[, s]]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  take_prod <- apply(takes_mat, 1L, prod)
  ## Every other PSU of the frame holds a unit, so one holds N - N_1 + 1.
  unit_max <- if (is.null(N_1_h)) prep$N_h else prep$N_h - N_1_h + 1
  too_big <- take_prod > unit_max * (1 + 1e-12)
  all_fixed_row <- apply(takes_fixed_mat, 1L, all)
  if (any(too_big & all_fixed_row)) {
    stop(
      sprintf(
        "%s exceeds %s for: %s",
        if (S == 2L) "'take_2'" else "product of stage takes",
        if (is.null(N_1_h)) "the stratum population" else
          "the units one PSU holds while every other PSU holds one",
        paste(prep$stratum[too_big & all_fixed_row], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(too_big & !all_fixed_row)) {
    for (h in which(too_big & !all_fixed_row)) {
      free_s <- which(!takes_fixed_mat[h, ])
      fixed_prod <- prod(takes_mat[h, takes_fixed_mat[h, ]])
      floor_s <- vapply(stage_resp_mat[h, free_s + 1L], .take_bound, 1)
      if (fixed_prod * prod(floor_s) > unit_max[h] * (1 + 1e-12)) {
        stop(
          sprintf(
            "the fixed takes of stratum '%s' with one expected response at its free stages need %g units per PSU, more than %s (%g)",
            prep$stratum[h], fixed_prod * prod(floor_s),
            if (is.null(N_1_h)) "the stratum population" else
              "one PSU holds while every other PSU holds one",
            unit_max[h]
          ),
          call. = FALSE
        )
      }
      takes_mat[h, free_s] <- .clamp_takes(takes_mat[h, free_s], floor_s,
                                           unit_max[h] / fixed_prod)
    }
    warning(
      sprintf(
        "cost-optimal %s %s",
        if (S == 2L) "'take_2' exceeds" else "stage takes exceed",
        if (is.null(N_1_h)) "the stratum population, clamped to 'N'" else
          "the units one PSU holds while every other PSU holds one, clamped to N - N_1 + 1"
      ),
      call. = FALSE
    )
  }

  if (any(prep$take_all)) {
    stop(
      "'take_all' is not supported for cluster allocation because taking every PSU does not imply an ultimate-unit census: within-PSU takes still apply. Drop 'take_all', or drop 'icc_1' to allocate at the element level",
      call. = FALSE
    )
  }

  prep$S_raw_h <- prep$S_h
  prep$fpc <- fpc
  prep$stages <- S
  prep$unit_resp_rate_h <- stage_resp_mat[, S]
  prep$stage_resp_mat <- stage_resp_mat
  prep$resp_rate_h <- apply(stage_resp_mat, 1L, prod)
  prep$icc_mat <- icc_mat
  prep$var_ratio_mat <- vr_mat
  prep$takes_fixed_mat <- takes_fixed_mat
  prep$cost_mat <- cost_mat
  prep$has_stage_costs <- has_costs
  prep$cluster <- TRUE
  # 2-stage convenience aliases used by 2-stage whole-cluster menu search
  prep$icc_psu_h <- icc_mat[, 1L]
  prep$var_ratio_psu_h <- vr_mat[, 1L]
  prep$n_per_psu_fixed_h <- takes_fixed_mat[, 1L]
  if (has_costs && S == 2L) {
    prep$cost_psu_h <- cost_mat[, 1L]
    prep$cost_ssu_h <- cost_mat[, 2L]
  }
  .cluster_at_takes(prep, takes_mat)
}

#' The cost-optimal free take of a two-stage stratum
#'
#' Minimizing cost times variance gives
#' `sqrt(cost_1 (1 - icc_1) / (cost_2 r D))`, with `D = icc_1` when no stage
#' correction applies. Under `fpc = "stage"` the within-PSU fraction
#' `take * r / M` lowers the between-PSU term to `D = icc_1 - (1 - icc_1) / M`.
#' When that is not positive, a further unit inside a PSU always buys more
#' precision than a further PSU, so the take is the whole PSU, `cap`, which
#' is `M` unless the take counts other units than `M` does.
#' @keywords internal
#' @noRd
.cluster_free_take <- function(cost_1, cost_2, icc, resp_rate, M = NULL,
                               cap = M) {
  if (is.null(M)) {
    take <- sqrt(cost_1 / cost_2 * (1 - icc) / (icc * resp_rate))
    return(if (is.null(cap)) take else pmin(take, cap))
  }
  D <- icc - (1 - icc) / M
  take <- rep_len(cap, length(D))
  inner <- D > 0
  take[inner] <- pmin(
    take[inner],
    sqrt(
      cost_1[inner] / cost_2[inner] * (1 - icc[inner]) /
        (D[inner] * resp_rate[inner])
    )
  )
  take
}

#' One stratum's S-stage inflation and PSU cost at a take vector
#' @keywords internal
#' @noRd
.cluster_stage_row <- function(prep, h, t_h,
                               scope_N_h = .alloc_scope_N(prep)[h]) {
  S <- prep$stages
  r_h <- prep$stage_resp_mat[h, ]
  eff_t <- t_h * r_h[2L:S]
  Q <- numeric(S)
  Q[S] <- 1
  for (s in (S - 1L):1L) {
    Q[s] <- Q[s + 1L] * eff_t[s]
  }
  vr <- prep$var_ratio_mat[h, ]
  icc <- prep$icc_mat[h, ]
  R_h <- c(vr * icc, vr[S - 1L] * (1 - icc[S - 1L]))
  f_stage <- numeric(S)
  if (identical(prep$fpc %||% "unit", "stage")) {
    N_st <- prep$N_stages_mat[h, ]
    for (s in 2L:S) {
      M_s <- if (s < S) N_st[s] / N_st[s - 1L] else scope_N_h / N_st[S - 1L]
      f_stage[s] <- min(1, eff_t[s - 1L] / M_s)
    }
  }
  list(
    inflation = R_h[1L] * Q[1L] +
      sum((1 - f_stage[2L:S]) * R_h[2L:S] * Q[2L:S]),
    psu_cost = if (isTRUE(prep$has_stage_costs)) {
      sum(prep$cost_mat[h, ] * c(1, cumprod(t_h)))
    } else {
      NA_real_
    }
  )
}

#' Compute exact S-stage stratum inflation, FPC population, and costs
#' @keywords internal
#' @noRd
.cluster_stage_terms <- function(prep, takes_mat) {
  H <- nrow(takes_mat)
  S <- prep$stages
  icc_mat <- prep$icc_mat
  vr_mat <- prep$var_ratio_mat
  rr_mat <- prep$stage_resp_mat
  scope_N <- .alloc_scope_N(prep)
  use_stage_fpc <- identical(prep$fpc %||% "unit", "stage")

  inflation <- numeric(H)
  R1_vec <- vr_mat[, 1L] * icc_mat[, 1L]
  take_prod <- numeric(H)
  psu_cost <- if (isTRUE(prep$has_stage_costs)) numeric(H) else rep(NA_real_, H)
  elem_cost <- if (isTRUE(prep$has_stage_costs)) numeric(H) else NULL

  for (h in seq_len(H)) {
    row <- .cluster_stage_row(prep, h, takes_mat[h, ], scope_N[h])
    inflation[h] <- row$inflation
    take_prod[h] <- prod(takes_mat[h, ])
    if (isTRUE(prep$has_stage_costs)) {
      psu_cost[h] <- row$psu_cost
      elem_cost[h] <- psu_cost[h] / take_prod[h]
    }
  }

  N_fpc_h <- switch(
    prep$fpc %||% "unit",
    unit = scope_N,
    none = rep(Inf, H),
    stage = ifelse(R1_vec > 0, prep$N_psu_h * inflation / R1_vec, Inf)
  )

  list(
    inflation = inflation,
    N_fpc_h = N_fpc_h,
    take_prod = take_prod,
    psu_cost = psu_cost,
    elem_cost = elem_cost
  )
}

#' Set an S-stage stratum's element terms at given stage takes (H x (S-1))
#' @keywords internal
#' @noRd
.cluster_at_takes <- function(prep, takes_mat) {
  terms <- .cluster_stage_terms(prep, takes_mat)
  prep$S_h <- prep$S_raw_h * sqrt(terms$inflation)
  prep$N_fpc_h <- terms$N_fpc_h
  prep$takes_mat <- takes_mat
  prep$take_prod_h <- terms$take_prod
  prep$n_per_psu_h <- takes_mat[, 1L]
  prep$responding_take_h <- takes_mat[, 1L] * prep$stage_resp_mat[, 2L]
  if (isTRUE(prep$has_stage_costs)) {
    prep$cost_h <- terms$elem_cost
  }
  prep
}

#' A two-stage stratum's element terms at a given take vector
#' @keywords internal
#' @noRd
.cluster_at_take <- function(prep, take) {
  takes_mat <- prep$takes_mat
  takes_mat[, 1L] <- take
  .cluster_at_takes(prep, takes_mat)
}

#' Free takes raised where the rule's allocation passes the PSU universe
#'
#' Past `N_1` PSUs the free takes carry `n_h / N_1` units per PSU.
#' @keywords internal
#' @noRd
.cluster_take_rise <- function(prep, n_h) {
  if (!isTRUE(prep$cluster) || is.null(prep$N_psu_h)) {
    return(prep)
  }
  free <- rowSums(!prep$takes_fixed_mat) > 0L
  P <- prep$take_prod_h
  spare <- prep$N_h - prep$N_psu_h
  rise <- free & n_h > prep$N_psu_h * P
  ## Past the room every other PSU leaves, the takes fall until n_h fits.
  fall <- free & !rise & P > 1 & n_h * (P - 1) / P > spare * (1 + 1e-12)
  if (!any(rise | fall)) {
    return(prep)
  }
  takes <- prep$takes_mat
  for (h in which(rise)) {
    takes[h, ] <- .cluster_split_takes(prep, h, n_h[h] / prep$N_psu_h[h])
  }
  for (h in which(fall)) {
    takes[h, ] <- .cluster_split_takes(prep, h, n_h[h] / (n_h[h] - spare[h]))
  }
  .cluster_at_takes(prep, takes)
}

#' Split a stratum's units per PSU across its free takes
#'
#' Cost times variance is log-convex, so pairwise moves along the product are exact.
#' @keywords internal
#' @noRd
.cluster_split_takes <- function(prep, h, per_psu) {
  t <- prep$takes_mat[h, ]
  free <- which(!prep$takes_fixed_mat[h, ])
  goal <- log(per_psu) - sum(log(t[-free]))
  lo <- log(vapply(prep$stage_resp_mat[h, free + 1L], .take_bound, 1))
  hi <- log(.cluster_take_caps(prep, h)[free])
  if (goal > sum(hi) + 1e-9 || goal < sum(lo) - 1e-9) {
    stop(
      sprintf(
        "stratum '%s' needs %g units per PSU, outside what its stage sizes and response floors allow",
        prep$stratum[h], per_psu
      ),
      call. = FALSE
    )
  }
  if (length(free) == 1L) {
    t[free] <- exp(min(max(goal, lo), hi))
    return(t)
  }
  u0 <- log(t[free])
  at <- function(c) pmin(pmax(u0 + c, lo), hi)
  u <- if (goal <= sum(lo)) {
    lo
  } else if (goal >= sum(hi)) {
    hi
  } else {
    span <- c(min(lo - u0) - 1, max(hi - u0) + 1)
    span[!is.finite(span)] <- sign(span[!is.finite(span)]) * 1e3
    at(stats::uniroot(function(c) sum(at(c)) - goal, span, tol = 1e-12)$root)
  }
  value <- function(u) {
    tt <- t
    tt[free] <- exp(u)
    row <- .cluster_stage_row(prep, h, tt)
    log(row$psu_cost) + log(row$inflation)
  }
  for (pass in seq_len(50L)) {
    before <- value(u)
    for (i in seq_len(length(u) - 1L)) {
      for (j in (i + 1L):length(u)) {
        pair <- u[i] + u[j]
        a <- max(lo[i], pair - hi[j])
        b <- min(hi[i], pair - lo[j])
        if (!(a < b)) next
        o <- stats::optimize(function(x) {
          v <- u
          v[i] <- x
          v[j] <- pair - x
          value(v)
        }, c(a, b), tol = 1e-10)
        u[i] <- o$minimum
        u[j] <- pair - o$minimum
      }
    }
    if (before - value(u) <= 1e-12 * max(1, abs(before))) break
  }
  t[free] <- exp(u)
  t
}

#' Physical caps on a stratum's takes, `Inf` where no stage size is known
#' @keywords internal
#' @noRd
.cluster_take_caps <- function(prep, h) {
  .stage_caps_of(prep$N_h[h], prep$N_stages_mat, h, prep$stages)
}

#' A stratum's units per parent unit at each stage past the first, `Inf` where unknown
#' @keywords internal
#' @noRd
.stage_caps_of <- function(N_h, N_st, h, S) {
  caps <- rep(Inf, S - 1L)
  if (is.null(N_st)) return(caps)
  for (s in 2L:S) {
    upper <- if (s < S) N_st[h, s] else N_h
    if (!is.na(upper) && !is.na(N_st[h, s - 1L])) {
      caps[s - 1L] <- upper / N_st[h, s - 1L]
    }
  }
  caps
}

#' The whole take a ratio of stage sizes admits, exact for whole sizes below 2^53
#' @keywords internal
#' @noRd
.whole_caps <- function(caps) {
  floor(caps)
}

#' Most PSUs of `P` units each that leave every other PSU of the frame one unit
#' @keywords internal
#' @noRd
.psu_room <- function(N, N_1, P) {
  if (is.null(N_1)) return(rep(Inf, length(P)))
  uncapped <- !(P > 1 + 1e-12) | !is.finite(P) | !is.finite(N_1)
  spare <- floor(N) - N_1
  step <- P - 1
  out <- floor(spare / step)
  over <- !uncapped & out * step > spare
  out[over] <- out[over] - 1
  under <- !uncapped & (out + 1) * step <= spare
  out[under] <- out[under] + 1
  out[uncapped] <- Inf
  out
}

#' The whole PSU counts each stratum may field at take products `P`
#' @keywords internal
#' @noRd
.ms_psu_range <- function(prep, P, lo_i, hi_i) {
  hi <- floor(hi_i / P)
  if (!is.null(prep$N_psu_h)) {
    hi <- pmin(hi, prep$N_psu_h, .psu_room(prep$N_h, prep$N_psu_h, P))
  }
  list(lo = ceiling(lo_i / P), hi = hi)
}

#' Free takes scaled down to a product, none below its floor
#' @keywords internal
#' @noRd
.clamp_takes <- function(t, floor_s, goal) {
  low <- rep(FALSE, length(t))
  repeat {
    rest <- goal / prod(floor_s[low])
    t[!low] <- t[!low] * (rest / prod(t[!low]))^(1 / sum(!low))
    under <- !low & t < floor_s
    if (!any(under)) break
    low <- low | under
    t[low] <- floor_s[low]
    if (all(low)) break
  }
  t
}

#' Whole takes within stage capacities
#' @keywords internal
#' @noRd
.cluster_round_takes <- function(takes, caps, lower = 1) {
  pmin(pmax(lower, round(takes)), pmax(1, .whole_caps(caps)))
}

#' The most units a stratum's stages hold, continuously and in whole units
#' @keywords internal
#' @noRd
.cluster_capacity <- function(prep) {
  H <- length(prep$N_h)
  cont <- whole <- numeric(H)
  floor_int <- .cluster_take_floor(prep)
  for (h in seq_len(H)) {
    caps <- .cluster_take_caps(prep, h)
    free <- !prep$takes_fixed_mat[h, ]
    t <- prep$takes_mat[h, ]
    N_1 <- prep$N_psu_h[h]
    ## Capacity min(N_1, room(P)) * P rises to P = N / N_1 and falls after it.
    floor_c <- vapply(prep$stage_resp_mat[h, -1L], .take_bound, 1)
    P_lo <- prod(ifelse(free, floor_c, t))
    P_hi <- prod(ifelse(free, caps, t))
    P <- min(max(prep$N_h[h] / N_1, P_lo), P_hi)
    a <- N_1
    if (is.finite(P) && P > 1) a <- min(a, (prep$N_h[h] - N_1) / (P - 1))
    cont[h] <- if (is.finite(P)) a * P else prep$N_h[h]
    fixed_w <- .cluster_round_takes(t, caps)
    whole[h] <- .ms_whole_capacity(
      prep$N_h[h], N_1,
      ifelse(free, floor_int[h, ], fixed_w),
      ifelse(free, .whole_caps(caps), fixed_w)
    )
  }
  list(continuous = cont, whole = whole)
}

#' The most whole units a stratum fields over every take product its stage bounds allow
#' @keywords internal
#' @noRd
.ms_whole_capacity <- function(N, N_1, lo_s, hi_s, max_boxes = 5000L) {
  A <- floor(N_1 + 1e-9)
  D <- N - A
  fit <- function(P) min(A, .psu_room(N, N_1, P))
  if (any(lo_s > hi_s) || fit(prod(lo_s)) < 1) return(0)
  for (k in seq_along(lo_s)) {
    hi_s[k] <- min(hi_s[k], floor((D + 2) / prod(lo_s[-k])))
  }
  P_top <- prod(hi_s)
  if (fit(P_top) >= A) return(A * P_top)
  ## Largest a * F * b over b in [L, U]: the best b falls as a grows, so
  ## each step jumps to the last PSU count that can use it.
  line <- function(F, L, U) {
    best <- 0
    a <- 1
    while (a <= A) {
      b <- min(U, floor((floor(D / a + 1e-9) + 1) / F))
      while (b >= L && fit(F * b) < a) b <- b - 1
      while (b + 1 <= U && fit(F * (b + 1)) >= a) b <- b + 1
      if (b < L) break
      last <- fit(F * b)
      best <- max(best, last * F * b)
      a <- last + 1
    }
    best
  }
  leaf <- function(lo, hi) {
    v <- which(lo < hi)
    if (length(v) == 0L) return(fit(prod(lo)) * prod(lo))
    line(prod(lo[-v]), lo[v], hi[v])
  }
  bound <- function(lo, hi) line(1, prod(lo), min(prod(hi), D + 2))
  best <- 0
  for (k in seq_along(lo_s)) {
    for (rest in list(lo_s, hi_s)) {
      lo <- rest
      hi <- rest
      lo[k] <- lo_s[k]
      hi[k] <- hi_s[k]
      best <- max(best, leaf(lo, hi))
    }
  }
  open <- list(list(lo = lo_s, hi = hi_s))
  boxes <- 0L
  while (length(open)) {
    box <- open[[length(open)]]
    open[[length(open)]] <- NULL
    if (bound(box$lo, box$hi) <= best) next
    v <- which(box$lo < box$hi)
    if (length(v) <= 1L) {
      best <- max(best, leaf(box$lo, box$hi))
      next
    }
    boxes <- boxes + 1L
    if (boxes > max_boxes) {
      open[[length(open) + 1L]] <- box
      rest <- vapply(open, function(b) bound(b$lo, b$hi), 1)
      return(structure(max(best, rest), exact = FALSE))
    }
    k <- v[which.min(box$hi[v] - box$lo[v])]
    mid <- floor((box$lo[k] + box$hi[k]) / 2)
    low <- high <- box
    low$hi[k] <- mid
    high$lo[k] <- mid + 1
    turn <- N / N_1
    far <- function(b) {
      max(prod(b$lo) - turn, turn - min(prod(b$hi), D + 2), 0)
    }
    kids <- if (far(low) <= far(high)) list(high, low) else list(low, high)
    open <- c(open, kids)
  }
  best
}

#' The smallest whole take with one expected response, per stratum and stage
#'
#' A fixed take keeps its own value, so its floor is one.
#' @keywords internal
#' @noRd
.cluster_take_floor <- function(prep) {
  rates <- prep$stage_resp_mat[, -1L, drop = FALSE]
  out <- matrix(vapply(rates, .take_floor, 1L), nrow = nrow(rates))
  out[prep$takes_fixed_mat] <- 1L
  out
}

#' @keywords internal
#' @noRd
.alloc_is_cluster <- function(frame) {
  .alloc_detect_stages(frame) >= 2L
}

#' Restore the takes held by an assessed or searched design
#' @keywords internal
#' @noRd
.alloc_frame_at_takes <- function(frame, takes) {
  for (col in colnames(takes)) frame[[col]] <- takes[, col]
  frame
}

#' Apply named `take_2`, ..., `take_S` overrides to a re-evaluated cluster allocation
#' @keywords internal
#' @noRd
.alloc_take_override <- function(frame, dots, cluster) {
  take_args <- grep("^take_[0-9]+$", names(dots), value = TRUE)
  if (length(take_args) == 0L) {
    return(list(frame = frame, dots = dots))
  }
  if (!cluster) {
    stop(
      sprintf(
        "'%s' applies to a cluster allocation, and this allocation has no cluster stages",
        take_args[1L]
      ),
      call. = FALSE
    )
  }
  H <- nrow(frame)
  for (nm in take_args) {
    take <- dots[[nm]]
    if (
      !is.numeric(take) ||
        !length(take) %in% c(1L, H) ||
        anyNA(take) ||
        any(!is.finite(take)) ||
        any(take <= 0)
    ) {
      stop(
        sprintf(
          "'%s' must be positive and finite, one value or one per stratum (%d)",
          nm,
          H
        ),
        call. = FALSE
      )
    }
    frame[[nm]] <- rep_len(as.numeric(take), H)
    dots[[nm]] <- NULL
  }
  list(frame = frame, dots = dots)
}

#' Explain an upper-bound failure in the language of what caused it
#'
#' A total, budget, or CV target the allocation cannot reach may be blocked by
#' the population, or by the PSU universe well below it. The second is not a
#' claim that the precision is impossible: this path keeps a between-PSU term
#' even when every PSU is taken, so a finite-population first stage could
#' reach targets the with-replacement approximation reports as out of range.
#' @keywords internal
#' @noRd
.alloc_psu_bound_note <- function(prep) {
  if (is.null(prep$M_src) || !any(prep$M_src == "N_1")) {
    return("")
  }
  room <- prep$M_src == "N_1" & (prep$psu_room_h %||% FALSE)
  if (any(room)) {
    return(sprintf(
      ". The ceiling comes from %s, whose sampled PSUs hold at most N - N_1 units more than their count, since every other PSU of the frame holds one",
      paste(sQuote(prep$stratum[room]), collapse = ", ")
    ))
  }
  rise <- prep$M_src == "N_1" & (prep$takes_rise_h %||% FALSE)
  if (any(rise)) {
    return(sprintf(
      ". The ceiling comes from the stage sizes in %s, where every PSU is taken and each free take is at its stage size per parent unit",
      paste(sQuote(prep$stratum[rise]), collapse = ", ")
    ))
  }
  if (identical(prep$fpc, "stage")) {
    return(sprintf(
      ". The ceiling comes from the PSU universe in %s, where every PSU is taken",
      paste(sQuote(prep$stratum[prep$M_src == "N_1"]), collapse = ", ")
    ))
  }
  sprintf(
    ". The ceiling comes from the PSU universe in %s rather than from the population. At the stage takes this allocation holds, only 'N_1' times the take product ultimate units are available there. That is what this path's with-replacement first-stage approximation and the supplied PSU availability allow, not a statement that the target is unattainable under the real design",
    paste(sQuote(prep$stratum[prep$M_src == "N_1"]), collapse = ", ")
  )
}

.alloc_bounds <- function(
  N_h,
  max_weight,
  take_all,
  min_n_stratum = NULL,
  psu_cap = NULL,
  step = NULL,
  stratum = NULL
) {
  H <- length(N_h)
  m_h <- pmin(rep(1, H), N_h)
  M_h <- as.numeric(N_h)
  # Which constraint produced each bound, so that .binding says not only that
  # the allocation sits on a bound but which one put it there.
  m_src <- rep(NA_character_, H)
  M_src <- rep("N", H)

  if (!is.null(min_n_stratum)) {
    floor_h <- rep_len(min_n_stratum, H)
    floor_h[is.na(floor_h)] <- 0
    hit <- floor_h > m_h
    m_h <- pmax(m_h, floor_h)
    m_src[hit] <- "min_n_stratum"
  }

  has_wmax <- !is.na(max_weight)
  if (any(has_wmax)) {
    idx <- which(has_wmax)
    wbound <- N_h[idx] / max_weight[idx]
    m_src[idx[wbound > m_h[idx]]] <- "max_weight"
    m_h[idx] <- pmax(m_h[idx], wbound)
  }

  if (!is.null(psu_cap)) {
    M_src[psu_cap < M_h] <- "N_1"
    M_h <- pmin(M_h, psu_cap)
  }

  if (any(take_all)) {
    m_h[take_all] <- N_h[take_all]
    M_h[take_all] <- N_h[take_all]
    m_src[take_all] <- "take_all"
    M_src[take_all] <- "take_all"
  }

  if (!is.null(step) && any(step > 1L)) {
    fits <- m_h <= M_h + 1e-8
    raised <- step * ceiling(m_h / step - 1e-9)
    m_src[raised > m_h + 1e-9 & is.na(m_src)] <- "n_step"
    m_h <- raised
    M_h <- step * floor(M_h / step + 1e-9)
    short <- fits & m_h > M_h + 1e-8
    if (any(short)) {
      stop(
        sprintf(
          "no whole number of steps fits the bounds of stratum %s. Its 'n_step' exceeds the room between its lower and upper bounds, so lower 'n_step' or relax the bounds",
          paste(stratum[short] %||% which(short), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  infeasible <- m_h > M_h + 1e-8
  if (any(infeasible)) {
    if (any(infeasible & M_src == "N_1")) {
      stop(
        "constraints are infeasible: the PSU universe 'N_1' caps the stratum below its own lower bound. Raise 'N_1' or stage takes, or relax 'min_n_stratum' or 'max_weight'",
        call. = FALSE
      )
    }
    stop(
      "constraints are infeasible: lower bounds exceed stratum population",
      call. = FALSE
    )
  }

  list(m_h = m_h, M_h = M_h, m_src = m_src, M_src = M_src)
}

#' Stratum standard deviation of a binary variable given as `p`
#'
#' An element frame under the unit FPC reads the finite-population variance
#' of the eligible population `N_h * e_h`, so the value moves with the
#' eligibility rate. A cluster frame and `fpc = "none"` read `p (1 - p)`.
#' @keywords internal
#' @noRd
.alloc_sd_from_p <- function(ph, N_h, e_h, frame, fpc) {
  if (.alloc_is_cluster(frame) || identical(fpc, "none")) {
    return(sqrt(ph * (1 - ph)))
  }
  if (any(N_h * e_h <= 1)) {
    stop(
      sprintf(
        "'%s' must be greater than 1 in every stratum when the variance is derived from 'p'",
        if (all(e_h == 1)) "N" else "N * eligibility_rate"
      ),
      call. = FALSE
    )
  }
  sqrt(.bernoulli_var(ph, N_h * e_h))
}

#' Stratum population in scope, N_h * e_h, or N_h for a frame prepared
#' without an eligibility rate
#' @keywords internal
#' @noRd
.alloc_scope_N <- function(prep) {
  prep$Ne_h %||% prep$N_h
}

#' Flag an assessed design whose PSUs hold more than the frame leaves them
#' @keywords internal
#' @noRd
.alloc_warn_psu_room <- function(prep, detail) {
  if (!anyNA(prep$N_stages_mat)) return(invisible(NULL))
  room <- prep$N_h - prep$N_psu_h + detail$n_1
  affected <- detail$n_1 < prep$N_psu_h - 1e-9 & detail$n > room * (1 + 1e-9)
  if (any(affected)) {
    warning(
      sprintf(
        "the sample in %s puts more units in its PSUs than the frame leaves them while every other PSU holds one (N - N_1 + n_1). Supply the missing intermediate stage populations to verify capacity",
        paste(sQuote(prep$stratum[affected]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.alloc_metrics <- function(
  N_h,
  S_h,
  mean_h,
  n_h,
  alpha,
  deff,
  resp_rate,
  cost_h,
  df = NULL,
  N_fpc = NULL,
  N_pop = NULL
) {
  N_pop <- N_pop %||% N_h
  n_net <- n_h * resp_rate
  n_eff <- n_net / deff
  W_h <- N_pop / sum(N_pop)
  if (all(is.na(S_h))) {
    na_h <- rep(NA_real_, length(n_h))
    return(list(
      se = NA_real_, moe = NA_real_, cv = NA_real_,
      cost = sum(n_h * cost_h), rmoe = NA_real_, term = na_h,
      se_h = na_h, moe_h = na_h, cv_h = na_h, rmoe_h = na_h, share_h = na_h
    ))
  }

  term <- numeric(length(n_h))
  good <- n_eff > 0
  # N_fpc is the population the correction divides by, which is the stratum
  # population unless a cluster allocation has folded a first-stage
  # correction into it. N_pop is the population in scope, N_h * e_h, which
  # the weights read. Bounds read the frame count N_h.
  fpc <- pmax(0, 1 - n_net / (N_fpc %||% N_pop))
  term[good] <- W_h[good]^2 * S_h[good]^2 * fpc[good] / n_eff[good]

  zero_ok <- !good & S_h == 0
  bad <- !good & !zero_ok
  if (any(bad)) {
    stop(
      "allocation implies zero effective sample in a stratum with positive variability",
      call. = FALSE
    )
  }

  V <- sum(term)
  se <- sqrt(max(V, 0))
  moe <- .q_alpha(alpha, df) * se

  cv <- NA_real_
  if (!all(is.na(mean_h))) {
    if (anyNA(mean_h)) {
      stop(
        "'mean' (or 'p') must be complete to compute aggregate CV",
        call. = FALSE
      )
    }
    ybar <- .aggregate_mean(W_h, mean_h)
    cv <- if (ybar == 0) Inf else se / abs(ybar)
  }

  # Stratum precision is the stratum's own mean, so the W_h^2 weight that
  # makes 'term' a contribution to the overall variance is divided back out.
  se_h <- numeric(length(n_h))
  se_h[good] <- sqrt(S_h[good]^2 * fpc[good] / n_eff[good])
  cv_h <- rep(NA_real_, length(n_h))
  if (!all(is.na(mean_h))) {
    cv_h <- ifelse(mean_h == 0, Inf, se_h / abs(mean_h))
  }

  moe_h <- .q_alpha(alpha, df) * se_h
  list(
    se = se,
    moe = moe,
    cv = cv,
    cost = sum(n_h * cost_h),
    rmoe = .rmoe_from_moe(
      moe,
      if (all(is.na(mean_h))) {
        NULL
      } else {
        .aggregate_mean(W_h, mean_h)
      }
    ),
    term = term,
    se_h = se_h,
    moe_h = moe_h,
    cv_h = cv_h,
    rmoe_h = .rmoe_from_moe(moe_h, if (all(is.na(mean_h))) NULL else mean_h),
    share_h = if (V > 0) term / V else rep(NA_real_, length(n_h))
  )
}

#' @keywords internal
#' @noRd
.alloc_detail <- function(
  prep,
  n_h,
  m_h,
  M_h,
  resp_rate,
  deff,
  mode = "n",
  budget = NULL,
  lo_i = NULL,
  hi_i = NULL,
  metrics = NULL
) {
  if (sum(n_h) > .Machine$integer.max) {
    .stop_count_range("the whole-unit rounding of n_alloc()", "in total")
  }
  if (is.null(lo_i)) {
    lo_i <- ceiling(pmax(ifelse(is.na(m_h), 0, m_h), 0) - 1e-9)
  }
  if (is.null(hi_i)) {
    hi_i <- floor(M_h + 1e-9)
  }
  step <- prep$step_h %||% rep(1L, length(n_h))
  stepped <- any(step > 1L)
  n_int <- switch(
    mode,
    cv = as.integer(pmin(
      pmax(step * ceiling(n_h / step - 1e-9), lo_i),
      hi_i
    )),
    budget = .round_within_budget(
      n_h,
      prep$cost_h,
      budget,
      lo_i,
      hi_i,
      (.alloc_scope_N(prep) / sum(.alloc_scope_N(prep)))^2 * prep$S_h^2 *
        rep_len(deff, length(n_h)) / rep_len(resp_rate, length(n_h)),
      step = step
    ),
    if (stepped) {
      .round_steps(n_h, step, lo_i, hi_i)
    } else {
      .round_oric_bounded(n_h, lo_i, hi_i)
    }
  )
  out <- data.frame(
    stratum = prep$stratum,
    N = prep$N_h,
    sd = prep$S_raw_h %||% prep$S_h,
    unit_cost = prep$cost_h,
    n = n_h,
    n_int = n_int,
    weight = prep$N_h / n_h,
    n_eff = n_h * resp_rate / deff,
    .lower = m_h,
    .upper = M_h,
    .binding = abs(n_h - m_h) < 1e-6 | abs(n_h - M_h) < 1e-6
  )
  on_lower <- abs(n_h - m_h) < 1e-6
  out$.bound_source <- ifelse(
    !out$.binding,
    NA_character_,
    ifelse(
      on_lower,
      prep$m_src %||% NA_character_,
      prep$M_src %||% NA_character_
    )
  )
  if (!is.null(metrics)) {
    out$.se <- metrics$se_h
    out$.moe <- metrics$moe_h
    out$.rmoe <- metrics$rmoe_h
    out$.cv <- metrics$cv_h
    out$.share <- metrics$share_h
  }
  if (isTRUE(prep$cluster)) {
    stages <- prep$stages %||% 2L
    for (s in 2L:stages) {
      out[[paste0("take_", s)]] <- prep$takes_mat[, s - 1L]
    }
    out$n_1 <- n_h / prep$take_prod_h
    n_1_int <- ceiling(out$n_1)
    if (!is.null(prep$N_stages_mat)) {
      n_1_int <- pmin(n_1_int, prep$N_stages_mat[, 1L])
      for (s in seq_len(stages - 1L)) {
        out[[paste0("N_", s)]] <- prep$N_stages_mat[, s]
      }
      out$.frac_1 <- out$n_1 / prep$N_stages_mat[, 1L]
    }
    out$n_1_int <- as.integer(n_1_int)
  }
  if (any(prep$take_all)) {
    out$take_all <- prep$take_all
  }
  if (stepped) {
    out$n_step <- step
  }
  if (!all(is.na(prep$mean_h))) {
    out$mean <- prep$mean_h
  }
  out
}

#' Integerize a budget-mode allocation without exceeding the budget
#'
#' Floors the continuous allocation (respecting integer lower bounds),
#' then greedily adds units where the variance reduction per unit cost
#' is largest while the budget allows. `a_h` is each stratum's variance
#' coefficient, `W_h^2 S_h^2 deff_h / resp_rate_h`, the term the reported
#' variance divides by `n_h`.
#' @keywords internal
#' @noRd
.round_within_budget <- function(n_h, cost_h, budget, lower, upper, a_h,
                                step = 1L) {
  step <- rep_len(as.integer(step), length(n_h))
  n_int <- pmin(pmax(as.integer(step * floor(n_h / step + 1e-9)), lower),
                upper)
  # The variance change of moving stratum h from n to n + its step. A
  # stratum without variance gains and loses nothing, and one left empty
  # with variance is infinitely worse, so no 0 * Inf reaches the comparison.
  change <- function(h, n) {
    out <- a_h[h] * (1 / n - 1 / (n + step[h]))
    out[a_h[h] == 0] <- 0
    out[a_h[h] > 0 & n == 0] <- Inf
    out
  }
  step_cost <- cost_h * step
  room <- .budget_room(budget, sum(lower * cost_h))
  # Lifting a stratum to its whole lower bound can spend more than the
  # floors elsewhere save. The caller checked that the bounds alone fit, so
  # trimming above them always comes back inside the budget.
  while (sum(n_int * cost_h) > room) {
    can <- which(n_int > lower)
    loss <- change(can, n_int[can] - step[can]) / step_cost[can]
    j <- can[which.min(loss)]
    n_int[j] <- n_int[j] - step[j]
  }
  repeat {
    can <- which(n_int < upper & step_cost <= room - sum(n_int * cost_h))
    if (length(can) == 0L) {
      break
    }
    gain <- change(can, n_int[can]) / step_cost[can]
    j <- can[which.max(gain)]
    n_int[j] <- n_int[j] + step[j]
  }
  as.integer(n_int)
}

#' Whole sizes in steps, each within one step of its share, at the total
#'
#' `.cluster_close_total()` rounds counts of blocks, here steps of `step`
#' units. A total the steps cannot reach within one step of each share is
#' an error naming the nearest totals they do reach.
#' @keywords internal
#' @noRd
.round_steps <- function(x, step, lo, hi) {
  total <- as.integer(round(sum(x)))
  a_min <- as.integer(lo %/% step)
  a_max <- as.integer(pmin(hi, total) %/% step)
  a0 <- pmin(pmax(as.integer(round(x / step)), a_min), a_max)
  a <- .cluster_close_total(a0, step, a_min, a_max, x)
  n_int <- as.integer(a * step)
  if (sum(n_int) != total) {
    reach <- attr(a, "totals", exact = TRUE)
    below <- reach[reach < total]
    above <- reach[reach > total]
    stop(
      sprintf(
        "no allocation in steps of 'n_step' reaches the total %d. Within one step of each stratum's share, the nearest totals are %s",
        total,
        paste(c(
          if (length(below)) max(below),
          if (length(above)) min(above)
        ), collapse = " and ")
      ),
      call. = FALSE
    )
  }
  n_int
}

#' The prepared frame of some strata, as a frame without domains
#' @keywords internal
#' @noRd
.alloc_prep_slice <- function(prep, idx) {
  H <- length(prep$N_h)
  out <- lapply(prep, function(v) {
    if (is.matrix(v) && nrow(v) == H) {
      v[idx, , drop = FALSE]
    } else if (is.atomic(v) && !is.null(v) && length(v) == H && is.null(dim(v))) {
      v[idx]
    } else {
      v
    }
  })
  out$frame <- prep$frame[idx, , drop = FALSE]
  out$domain_cols <- character(0)
  out$domain_idx <- list()
  out$domain_values <- NULL
  out
}

#' Whether the cv mode with domains searches each domain's whole design
#'
#' The domain cv mode is each domain at its own minimum, so the exact
#' search of the mode without domains applies to each domain's strata, on
#' element and cluster frames alike.
#' @keywords internal
#' @noRd
.domain_searched <- function(alloc, mode, domains) {
  domains && mode == "cv" && alloc %in% c("optimal", "neyman")
}

#' Whether a classic allocation searches its whole design
#'
#' The cv mode of `"optimal"` and `"neyman"`, whose continuous allocation is
#' the minimum-cost one, and the budget mode of `"optimal"`, without
#' domains, for element and two-stage frames alike. The other modes and
#' rules are rounded as rules.
#' @keywords internal
#' @noRd
.classic_searched <- function(alloc, mode, domains) {
  !domains &&
    ((mode == "cv" && alloc %in% c("optimal", "neyman")) ||
       (mode == "budget" && alloc == "optimal"))
}

#' The whole allocation of a classic cv or budget mode, by exact search
#'
#' A stratum's variance term is `a_h / n_h - c_h`, with
#' `a_h = W_h^2 S_h^2 deff_h / resp_h` and `c_h = W_h^2 S_h^2 deff_h / N_h`,
#' so a cv target is one column. It is measured from the upper bounds `hi`,
#' `sum_h a_h (1 / n_h - 1 / hi_h) <= (cv ybar)^2 - V(hi)`, so a census
#' loads exactly zero however small the target. A plan meets the target
#' when its variance is within a relative 1e-12 of it, a matter of rounding,
#' so the design whose cv was read back as the target is attainable.
#' `"optimal"` minimizes cost and `"neyman"` the sample. Under a budget `"optimal"`
#' minimizes the variance, with cost over the budget as the column. A
#' stratum without variance takes its lower bound. Returns the whole sizes,
#' the search's status and its record.
#' @keywords internal
#' @noRd
.alloc_whole_element <- function(prep, alloc, mode, cv, budget, deff,
                                 resp_rate, lo_i, hi_i, max_nodes, gap = 0,
                                 target_se = NULL) {
  H <- length(prep$N_h)
  W <- .alloc_scope_N(prep) / sum(.alloc_scope_N(prep))
  S2 <- W^2 * prep$S_h^2 * rep_len(deff, H)
  a <- S2 / rep_len(resp_rate, H)
  c0 <- S2 / (prep$N_fpc_h %||% .alloc_scope_N(prep))
  cost <- if (identical(alloc, "neyman")) rep(1, H) else prep$cost_h
  # The search runs on whole steps, n = step * x. A stratum's cost per step
  # is step * cost and its variance coefficient per step a / step.
  step <- as.numeric(prep$step_h %||% rep(1L, H))
  n <- as.numeric(lo_i)
  on <- a > 0
  k <- step[on]
  lo <- pmax(lo_i[on] / k, 1)
  hi <- hi_i[on] / k
  if (!any(on)) {
    value <- if (mode == "cv") sum(cost * n) else 0
    r <- list(value = value, bound = value, status = "optimal", nodes = 0,
              states = 0, kept = 0)
    return(list(n_int = as.integer(n), status = r$status,
                search = .whole_record(r, "whole units", max_nodes)))
  }
  if (mode == "cv") {
    se_goal <- if (!is.null(target_se)) {
      target_se
    } else {
      cv * abs(.aggregate_mean(W, prep$mean_h))
    }
    room <- se_goal^2 - sum(a[on] / hi_i[on] - c0[on])
    # Within rounding of the target, the upper bounds meet it.
    if (room < -1e-12 * se_goal^2) {
      stop(sprintf(
        "target precision is unattainable with whole units under current constraints%s",
        .alloc_psu_bound_note(prep)
      ), call. = FALSE)
    }
    if (room <= 0) {
      n[on] <- hi_i[on]
      r <- list(value = sum(cost * n), bound = sum(cost * n),
                status = "optimal", nodes = 0, states = 0, kept = 0)
      return(list(n_int = as.integer(n), status = r$status,
                  search = .whole_record(r, "whole units", max_nodes)))
    }
    # Measured from the upper bounds, so a census loads exactly zero however
    # small the target.
    g2 <- a[on] / k / room * .whole_target_scale()
    q <- .whole_problem(cost[on] * k, 0, 0, g2, lo, hi, g0 = -(g2 / hi))
    r <- .whole_search(q, grain = .menu_grain(cost[on] * k),
                       max_nodes = max_nodes, gap = gap)
    offset <- sum(cost[!on] * n[!on])
  } else {
    fixed <- sum(prep$cost_h[!on] * n[!on])
    room <- .budget_room(budget, fixed + sum(prep$cost_h[on] * k * lo)) - fixed
    q <- .whole_problem(0, a[on] / k,
                        .whole_budget_load(prep$cost_h[on] * k, room), 0,
                        lo, hi)
    r <- .whole_search(q, start = lo, max_nodes = max_nodes, gap = gap)
    offset <- -sum(c0)
  }
  if (identical(r$status, "infeasible")) {
    stop(sprintf(
      "target 'cv' is unattainable with whole units under current constraints%s",
      .alloc_psu_bound_note(prep)
    ), call. = FALSE)
  }
  if (is.null(r$x)) {
    stop(structure(
      list(message = sprintf(
        "no whole allocation was found in %d search node%s, so whether one exists is unknown. Raise 'control$max_nodes'",
        as.integer(max_nodes), if (max_nodes == 1) "" else "s"
      ), call = NULL),
      class = c("svyplan_plan_unknown", "error", "condition")
    ))
  }
  n[on] <- r$x * k
  list(n_int = as.integer(round(n)), status = r$status,
       search = .whole_record(r, "whole units", max_nodes, offset = offset))
}

#' The whole field design of a two-stage cv or budget mode, by exact search
#'
#' A stratum's state is `a` whole PSUs at a whole take `b`. At take `b` its
#' variance term is `C_b / a - D_b`, with `D_b` the finite population part,
#' and its cost is `a` times the PSU's cost at `b`. Each take is a branch of
#' `.whole_search()`. In cv mode a stratum's load is measured from its
#' lightest listed state, `C_b / a - D_b - vmin`, against the target less
#' the sum of `vmin`, so no listed load is negative and a target near the
#' smallest reachable variance stays resolvable. In budget mode `D_b`
#' enters the objective as `Dmax - D_b`, `Dmax` the largest over the
#' listed takes. Under `"optimal"` with stage costs and `N_1`, a take the
#' frame leaves free ranges over `1..b_max`, `b_max` the whole part of the
#' average PSU size `N / N_1`. Every other take stays at
#' `b_rule`, the whole take of the rule, since without a cost per PSU the
#' smallest sample would always take one unit a PSU, and without `N_1`
#' no PSU size bounds the take.
#'
#' Free takes are listed up to a window, first the rule's take and 50
#' more, doubled while the lightest listed states miss a cv target. The
#' takes above it are proven out. On each piece above the window `alpha` and `gamma` of the reduced
#' cost `gamma + alpha a + beta / a` are linear in the take and `beta`
#' falls, so `gamma + alpha + beta_lin / hi` and
#' `gamma + 2 sqrt(alpha beta_end)` are concave lower bounds, with
#' `beta_lin` the part of `beta b` linear in the take and `beta_end` the
#' piece's last `beta`. A piece whose bound passes the stratum's `top` at
#' both ends holds no state of a better plan. Otherwise the window grows
#' over that piece and the search runs again.
#' @keywords internal
#' @noRd
.alloc_whole_cluster <- function(prep, alloc, mode, cv, budget, deff,
                                 resp_rate, lo_i, hi_i, b_rule, max_nodes,
                                 gap = 0, target_se = NULL) {
  H <- length(prep$N_h)
  W <- .alloc_scope_N(prep) / sum(.alloc_scope_N(prep))
  r_u <- rep_len(prep$unit_resp_rate_h %||% resp_rate, H)
  r <- rep_len(resp_rate, H)
  K <- W^2 * prep$S_raw_h^2 * rep_len(deff, H) * prep$var_ratio_psu_h
  icc <- prep$icc_psu_h
  fpc <- prep$fpc %||% "unit"
  N_psu <- prep$N_psu_h %||% rep(Inf, H)
  costed <- alloc == "optimal" && isTRUE(prep$has_stage_costs)
  # A take moves only where N_psu gives a PSU size to bound it. Under the
  # unit correction an unbounded take lets one PSU draw a census.
  free <- costed & !prep$n_per_psu_fixed_h & is.finite(N_psu)
  b_max <- ifelse(
    free,
    pmin(hi_i, pmax(1, floor(prep$N_h / N_psu + 1e-9))),
    b_rule
  )
  b_min <- pmin(.cluster_take_floor(prep)[, 1L], b_max)

  terms <- function(h, b) {
    wf <- if (fpc == "stage") {
      pmin(1, b * r_u[h] / (.alloc_scope_N(prep)[h] / N_psu[h]))
    } else {
      0
    }
    g <- icc[h] * b * r_u[h] + (1 - wf) * (1 - icc[h])
    list(
      b = b,
      lo = pmax(1, ceiling(lo_i[h] / b)),
      hi = pmin(floor(hi_i[h] / b), N_psu[h]),
      pc = if (costed) prep$cost_psu_h[h] + prep$cost_ssu_h[h] * b else b,
      C = K[h] * g / (b * r[h]),
      D = switch(fpc,
                 unit = K[h] * g / .alloc_scope_N(prep)[h],
                 none = 0 * b,
                 stage = rep(K[h] * icc[h] / N_psu[h], length(b))),
      lin = K[h] * icc[h] * r_u[h] * b / (r[h] * hi_i[h])
    )
  }
  takes <- function(h, w) {
    if (free[h]) seq.int(b_min[h], max(b_min[h], w[h])) else b_rule[h]
  }
  # D_b grows with the take, or stays constant, so its largest listed value
  # is at the window's end.
  d_max <- function(w) {
    vapply(seq_len(H), function(h) terms(h, w[h])$D, 1)
  }
  # Each stratum's state with the smallest load or cost, among the listed
  # takes, and that value.
  lightest <- function(w, value) {
    lapply(seq_len(H), function(h) {
      t <- terms(h, takes(h, w))
      ok <- t$lo <= t$hi
      v <- value(h, t)
      i <- which(ok)[which.min(v[ok])]
      list(b = t$b[i], v = v[i], a = if (mode == "cv") t$hi[i] else t$lo[i])
    })
  }
  w <- ifelse(free, pmin(b_max, b_rule + 50), b_rule)

  if (mode == "cv") {
    se_goal <- if (!is.null(target_se)) {
      target_se
    } else {
      cv * abs(.aggregate_mean(W, prep$mean_h))
    }
    # The largest sample a take allows has the smallest variance. When the
    # window's floor misses the target, a larger take may meet it.
    repeat {
      low <- lightest(w, function(h, t) t$C / t$hi - t$D)
      # Within rounding of the target, the floor meets it.
      if (sum(vapply(low, `[[`, 1, "v")) <= se_goal^2 * (1 + 1e-12)) {
        break
      }
      if (all(w >= b_max)) {
        stop(sprintf(
          "no whole-cluster design found that attains the target precision under the integer bounds%s",
          .alloc_psu_bound_note(prep)
        ), call. = FALSE)
      }
      w <- pmin(b_max, ifelse(free, 2 * w, w))
    }
    grain <- if (costed) {
      .menu_grain(c(prep$cost_psu_h, prep$cost_ssu_h))
    } else {
      1
    }
  } else {
    # Past the lower bound one PSU suffices and cost grows with the take.
    w <- pmax(w, ifelse(free, pmin(b_max, pmax(1, lo_i)), w))
    low <- lightest(w, function(h, t) t$lo * t$pc)
    min_cost <- sum(vapply(low, `[[`, 1, "v"))
    if (min_cost > budget * (1 + 1e-8)) {
      stop(sprintf(
        "'budget' cannot fund one whole PSU per stratum while satisfying the lower bounds (minimum field cost = %.4g)",
        min_cost
      ), call. = FALSE)
    }
    room <- .budget_room(budget, min_cost)
    grain <- 0
  }
  shape <- function(h, t, Dmax) {
    if (mode == "cv") {
      list(o0 = 0, o1 = t$pc, o2 = 0, g0 = (-t$D - vmin[h]) * m, g1 = 0,
           g2 = t$C * m)
    } else {
      list(o0 = Dmax[h] - t$D, o1 = 0, o2 = t$C, g0 = 0,
           g1 = .whole_budget_load(t$pc, room),
           g2 = 0)
    }
  }

  repeat {
    Dmax <- d_max(w)
    low <- lightest(w, if (mode == "cv") {
      function(h, t) t$C / t$hi - t$D
    } else {
      function(h, t) t$lo * t$pc
    })
    if (mode == "cv") {
      vmin <- vapply(low, `[[`, 1, "v")
      room <- se_goal^2 - sum(vmin)
      if (room <= 0) {
        if (any(w < b_max)) {
          w <- b_max
          next
        }
        # A target at the smallest variance, to within rounding, holds every
        # stratum at its own smallest, so each takes its cheapest state there.
        # A stratum without variance has it at its fewest PSUs.
        best <- lapply(seq_len(H), function(h) {
          t <- terms(h, takes(h, w))
          ok <- t$lo <= t$hi
          t <- lapply(t, function(v) v[ok])
          a <- vapply(seq_along(t$b), function(i) {
            fits <- function(a) t$C[i] / a - t$D[i] <= vmin[h]
            if (!fits(t$hi[i])) return(NA_real_)
            if (t$C[i] == 0) return(t$lo[i])
            a <- min(t$hi[i],
                     max(t$lo[i], ceiling(t$C[i] / (vmin[h] + t$D[i]))))
            while (a > t$lo[i] && fits(a - 1)) a <- a - 1
            while (!fits(a)) a <- a + 1
            a
          }, 1)
          i <- which.min(a * t$pc)
          list(a = a[i], b = t$b[i], cost = a[i] * t$pc[i])
        })
        value <- sum(vapply(best, `[[`, 1, "cost"))
        r_s <- list(value = value, bound = value, status = "optimal",
                    nodes = 0, states = 0, kept = 0)
        return(list(a = as.integer(vapply(best, `[[`, 1, "a")),
                    b = as.integer(vapply(best, `[[`, 1, "b")),
                    status = "optimal",
                    search = .whole_record(r_s, "whole units", max_nodes)))
      }
      m <- .whole_target_scale() / room
    }
    rows <- lapply(seq_len(H), function(h) {
      t <- terms(h, takes(h, w))
      ok <- t$lo <= t$hi
      t <- lapply(t, function(v) v[ok])
      z <- shape(h, t, Dmax)
      if (mode == "cv") {
        # The lightest state loads exactly zero, so a target near the
        # smallest reachable variance stays within the search's tolerance.
        j <- match(low[[h]]$b, t$b)
        z$g0 <- rep_len(z$g0, length(t$b))
        z$g0[j] <- -(z$g2[j] / t$hi[j])
      }
      c(list(h = rep(h, sum(ok))), t, z)
    })
    col <- function(name) {
      unlist(lapply(rows, function(z) rep_len(z[[name]], length(z$h))),
             use.names = FALSE)
    }
    owner <- col("h")
    q <- .whole_problem(col("o1"), col("o2"), col("g1"), col("g2"),
                        col("lo"), col("hi"), owner = owner,
                        o0 = col("o0"), g0 = col("g0"))
    first <- c(0, cumsum(lengths(lapply(rows, `[[`, "h"))))[seq_len(H)]
    at <- first + vapply(seq_len(H), function(h) {
      match(low[[h]]$b, rows[[h]]$b)
    }, 1L)
    start_a <- vapply(low, `[[`, 1, "a")
    r_s <- .whole_search(q, start = .whole_plan(q, start_a, branch = at),
                         grain = grain, max_nodes = max_nodes, gap = gap)
    if (is.null(r_s$x)) {
      stop(structure(
        list(message = sprintf(
          "no whole-cluster design was found in %d search node%s, so whether one exists is unknown. Raise 'control$max_nodes'",
          as.integer(max_nodes), if (max_nodes == 1) "" else "s"
        ), call = NULL),
        class = c("svyplan_plan_unknown", "error", "condition")
      ))
    }
    grown <- w
    for (h in which(free & w < b_max)) {
      lam <- r_s$lambda
      reduced <- function(b) {
        t <- terms(h, b)
        z <- shape(h, t, Dmax)
        list(alpha = z$o1 + z$g1 * lam, beta = z$o2 + z$g2 * lam,
             gamma = z$o0 + z$g0 * lam,
             lin = t$lin * (if (mode == "cv") m * lam else 1))
      }
      from <- w[h] + 1
      while (from <= b_max[h]) {
        to <- min(b_max[h], 2 * from)
        e <- reduced(c(from, to))
        bound_a <- min(e$gamma + e$alpha + e$lin)
        bound_b <- min(e$gamma + 2 * sqrt(e$alpha * e$beta[2]))
        if (max(bound_a, bound_b) <= r_s$top[h]) {
          grown[h] <- to
          break
        }
        from <- to + 1
      }
    }
    if (all(grown == w)) break
    w <- grown
  }
  list(
    a = as.integer(r_s$x),
    b = as.integer(col("b")[r_s$branch]),
    status = r_s$status,
    search = .whole_record(r_s, "whole units", max_nodes,
                           offset = if (mode == "cv") 0 else -sum(Dmax))
  )
}

#' Operational (integer) design for an element allocation
#'
#' Recomputes cost and precision from the integer allocation so the
#' operational metrics describe the fieldable design, not the
#' continuous optimum.
#' @keywords internal
#' @noRd
.alloc_operational_element <- function(
  prep,
  detail,
  alpha,
  deff,
  resp_rate,
  df = NULL
) {
  n_int <- detail$n_int
  m <- .alloc_metrics(
    N_h = prep$N_h,
    N_pop = .alloc_scope_N(prep),
    S_h = prep$S_h,
    mean_h = prep$mean_h,
    n_h = n_int,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = prep$cost_h,
    df = df,
    N_fpc = prep$N_fpc_h
  )
  list(
    n = sum(n_int),
    cost = sum(n_int * prep$cost_h),
    se = m$se,
    moe = m$moe,
    cv = m$cv,
    rmoe = m$rmoe
  )
}

#' Domain floors resolved against the prepared domains
#'
#' `min_n_domain` is one number for every domain or a table of `domain`,
#' `level` and `min_n`, where `domain` names the domain column. Levels the
#' table leaves out have no floor. Returns the floor of each domain in the
#' order of `prep$domain_idx` and the domain index of each stratum, or
#' `NULL` without floors.
#' @keywords internal
#' @noRd
.alloc_domain_floors <- function(min_n_domain, prep) {
  if (is.null(min_n_domain)) {
    return(NULL)
  }
  cols <- prep$domain_cols
  if (length(cols) == 0L) {
    stop("'min_n_domain' needs 'domains', the column whose levels it floors",
         call. = FALSE)
  }
  if (length(cols) > 1L) {
    stop(
      "'min_n_domain' applies with one 'domains' column, whose domains partition the strata. Floors on several columns overlap",
      call. = FALSE
    )
  }
  levels_d <- as.character(prep$domain_values[[cols]])
  D <- length(levels_d)
  if (is.data.frame(min_n_domain)) {
    need <- c("domain", "level", "min_n")
    if (!all(need %in% names(min_n_domain))) {
      stop("a 'min_n_domain' table needs columns 'domain', 'level' and 'min_n'",
           call. = FALSE)
    }
    if (any(as.character(min_n_domain[["domain"]]) != cols)) {
      stop(sprintf("'min_n_domain' column 'domain' must name '%s'", cols),
           call. = FALSE)
    }
    lev <- as.character(min_n_domain[["level"]])
    if (anyDuplicated(lev)) {
      stop("'min_n_domain' names a level more than once", call. = FALSE)
    }
    unknown <- setdiff(lev, levels_d)
    if (length(unknown) > 0L) {
      stop(
        sprintf("'min_n_domain' names level(s) not in '%s': %s", cols,
                paste(sQuote(unknown), collapse = ", ")),
        call. = FALSE
      )
    }
    value <- min_n_domain[["min_n"]]
    floor_d <- rep(0, D)
    floor_d[match(lev, levels_d)] <- value
  } else {
    value <- min_n_domain
    if (!is.numeric(value) || length(value) != 1L) {
      stop("'min_n_domain' must be one number or a table of domain floors",
           call. = FALSE)
    }
    floor_d <- rep(value, D)
  }
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
      any(value < 0)) {
    stop("'min_n_domain' floors must be non-negative finite numbers",
         call. = FALSE)
  }
  dom <- integer(length(prep$N_h))
  for (d in seq_len(D)) dom[prep$domain_idx[[d]]] <- d
  list(floor = as.numeric(floor_d), dom = dom, levels = levels_d)
}

#' Restore domain floors a rounding lost
#'
#' Works in steps: an element, or a whole PSU of `unit` elements in a
#' cluster stratum. A short domain gains a step where the rounding fell
#' furthest below the continuous design. To keep a total (`n` mode) or a
#' budget, a domain with room above its floor gives a step back where the
#' rounding rose furthest above it.
#' @keywords internal
#' @noRd
.alloc_floor_repair <- function(x, cont, unit, lo, hi, dom, need,
                                keep = "none", total = NULL,
                                cost_step = NULL, budget = NULL) {
  D <- length(need)
  dom_total <- function(x) {
    vapply(seq_len(D), function(d) sum((x * unit)[dom == d]), numeric(1))
  }
  repeat {
    short <- which(dom_total(x) < need - 1e-9)
    if (length(short) == 0L) break
    d <- short[1L]
    cand <- which(dom == d & x < hi)
    if (length(cand) == 0L) {
      stop(
        sprintf(
          "no whole design meets the floor of domain %d within its stratum bounds",
          d
        ),
        call. = FALSE
      )
    }
    h <- cand[which.max(cont[cand] - x[cand])]
    x[h] <- x[h] + 1
  }
  over <- function(x) {
    if (keep == "total") {
      sum(x * unit) - total
    } else if (keep == "budget") {
      sum(x * cost_step) - budget
    } else {
      -Inf
    }
  }
  slack <- if (keep == "budget") .budget_slack(budget) else 1e-9
  while (over(x) > slack) {
    room <- dom_total(x)[dom] - unit >= need[dom] - 1e-9 & x > lo
    if (keep == "total") room <- room & unit <= over(x) + 1e-9
    cand <- which(room)
    if (length(cand) == 0L) {
      if (keep == "budget") {
        stop(
          "no whole design meets every 'min_n_domain' floor within the budget",
          call. = FALSE
        )
      }
      break
    }
    h <- cand[which.min(cont[cand] - x[cand])]
    x[h] <- x[h] - 1
  }
  x
}

#' Whole steps that meet every domain floor and the total exactly
#'
#' The floors are met first, by `.alloc_floor_repair()`. Each stratum then
#' moves at most `width` steps from that design, and the move closest to the
#' continuous design that keeps each domain at or above its floor and the
#' sum at `total` is kept. A dynamic program per domain lists the domain
#' totals its strata reach, with the least squared distance for each, and a
#' second one over domains reaches the total. `NULL` when none does, with
#' the totals reached as attribute `totals`.
#' @keywords internal
#' @noRd
.floor_close_total <- function(a0, step, a_min, a_max, dom, need, target,
                               total, width = 1L) {
  reach <- function(states, values, costs) {
    out <- list()
    for (key in names(states)) {
      s <- states[[key]]
      for (k in seq_along(values)) {
        t <- s$tot + values[k]
        c <- s$cost + costs[k]
        id <- as.character(t)
        if (is.null(out[[id]]) || c < out[[id]]$cost - 1e-12) {
          out[[id]] <- list(tot = t, cost = c, pick = c(s$pick, k))
        }
      }
    }
    out
  }
  D <- length(need)
  moves <- vector("list", length(a0))
  per_domain <- vector("list", D)
  for (d in seq_len(D)) {
    states <- list("0" = list(tot = 0, cost = 0, pick = integer(0)))
    for (h in which(dom == d)) {
      m <- seq.int(max(a_min[h] - a0[h], -width), min(a_max[h] - a0[h], width))
      moves[[h]] <- m
      e <- step[h] * (a0[h] + m)
      states <- reach(states, e, (e - target[h])^2)
    }
    states <- states[vapply(states, function(s) s$tot >= need[d] - 1e-9,
                            logical(1))]
    if (length(states) == 0L) {
      return(NULL)
    }
    per_domain[[d]] <- states
  }
  combined <- list("0" = list(tot = 0, cost = 0, pick = integer(0)))
  for (d in seq_len(D)) {
    opts <- per_domain[[d]]
    combined <- reach(combined, vapply(opts, `[[`, numeric(1), "tot"),
                      vapply(opts, `[[`, numeric(1), "cost"))
  }
  hit <- combined[[as.character(total)]]
  if (is.null(hit)) {
    return(structure(list(),
                     totals = sort(vapply(combined, `[[`, numeric(1), "tot"))))
  }
  a <- a0
  for (d in seq_len(D)) {
    chosen <- per_domain[[d]][[hit$pick[d]]]
    idx <- which(dom == d)
    a[idx] <- a0[idx] + mapply(function(h, k) moves[[h]][k], idx, chosen$pick)
  }
  a
}

#' Whole PSU counts, each within one PSU of its rounding, nearest the total
#' @keywords internal
#' @noRd
.cluster_close_total <- function(a_h, b_h, a_min, a_max, e_target) {
  H <- length(a_h)
  target <- sum(e_target)
  base <- sum(a_h * b_h)
  if (base == target) {
    return(structure(a_h, totals = base))
  }
  lo_d <- as.integer(pmax(a_min - a_h, -1L))
  hi_d <- as.integer(pmin(a_max - a_h, 1L))
  offset <- sum(-lo_d * b_h)
  size <- offset + sum(hi_d * b_h) + 1L
  cost <- rep(Inf, size)
  cost[offset + 1L] <- 0
  moves <- rep(0L, size)
  back <- matrix(NA_integer_, size, H)
  for (h in seq_len(H)) {
    next_cost <- rep(Inf, size)
    next_moves <- rep(0L, size)
    arg <- rep(NA_integer_, size)
    d_all <- seq.int(lo_d[h], hi_d[h])
    for (d in d_all[order(abs(d_all))]) {
      shift <- d * b_h[h]
      from <- max(1L, 1L - shift):min(size, size - shift)
      to <- from + shift
      cand <- cost[from] + (b_h[h] * (a_h[h] + d) - e_target[h])^2
      cand_moves <- moves[from] + abs(d)
      better <- cand < next_cost[to] |
        (cand == next_cost[to] & cand_moves < next_moves[to])
      next_cost[to[better]] <- cand[better]
      next_moves[to[better]] <- cand_moves[better]
      arg[to[better]] <- d
    }
    cost <- next_cost
    moves <- next_moves
    back[, h] <- arg
  }
  dev <- seq_len(size) - 1L - offset
  gap <- ifelse(is.finite(cost), abs(target - base - dev), Inf)
  t <- order(gap, cost, moves)[1L]
  d <- integer(H)
  for (h in rev(seq_len(H))) {
    d[h] <- back[t, h]
    t <- t - d[h] * b_h[h]
  }
  structure(as.integer(a_h + d), totals = base + dev[is.finite(cost)])
}

#' Whole-unit field design of a stratified cluster allocation
#' @keywords internal
#' @noRd
.alloc_operational_cluster <- function(
  prep,
  detail,
  n_h,
  mode,
  budget,
  alpha,
  deff,
  resp_rate,
  target_cv = NULL,
  target_se = NULL,
  df = NULL,
  alloc_weights = NULL,
  alloc = NULL,
  max_nodes = NULL,
  gap = 0,
  floors = NULL,
  domain_search = NULL
) {
  if ((prep$stages %||% 2L) >= 3L) {
    return(.alloc_operational_cluster_multistage(
      prep, detail, n_h, mode, target_cv = target_cv, target_se = target_se,
      budget = budget, alpha = alpha, deff = deff, resp_rate = resp_rate,
      df = df, alloc_weights = alloc_weights, floors = floors, alloc = alloc,
      max_nodes = max_nodes, gap = gap, domain_search = domain_search
    ))
  }
  H <- length(n_h)
  ps <- prep$takes_mat[, 1L]
  # Ultimate-unit rate, not the combined one: PSU loss stays out of the bracket.
  resp_h <- rep_len(prep$unit_resp_rate_h %||% resp_rate, H)
  icc <- prep$icc_psu_h
  var_ratio <- prep$var_ratio_psu_h
  S_raw <- prep$S_raw_h
  W <- .alloc_scope_N(prep) / sum(.alloc_scope_N(prep))
  lo_i <- ceiling(detail$.lower - 1e-9)
  hi_i <- floor(detail$.upper + 1e-9)
  # N_1 bounds the PSU count through a_max, not through this element cap.
  psu_capped <- prep$M_src %in% "N_1"
  hi_i[psu_capped] <- floor(prep$N_h[psu_capped] + 1e-9)
  psu_max <- if (is.null(prep$N_psu_h)) rep(Inf, H) else prep$N_psu_h
  take_floor <- .cluster_take_floor(prep)[, 1L]
  take_cap <- ifelse(is.finite(psu_max), .whole_caps(prep$N_h / psu_max), Inf)
  take_options <- lapply(
    seq_len(H),
    function(h) {
      center <- max(1L, as.integer(round(ps[h])))
      cand <- unique(c(
        max(1L, as.integer(floor(ps[h]))),
        max(1L, as.integer(ceiling(ps[h]))),
        seq.int(max(1L, center - 50L), min(hi_i[h], center + 50L)),
        1L
      ))
      # A take-all or other exact bound must be representable as a product
      # of whole PSUs and a whole, common take.  Including the exact bound
      # guarantees a feasible divisor when it is an integer.
      if (lo_i[h] == hi_i[h]) {
        cand <- unique(c(cand, hi_i[h]))
      }
      cand <- cand[cand >= 1L & cand <= min(hi_i[h], take_cap[h])]
      if (any(cand >= take_floor[h])) cand <- cand[cand >= take_floor[h]]
      feasible <- vapply(
        cand,
        function(b) {
          ceiling(lo_i[h] / b) <= min(floor(hi_i[h] / b), psu_max[h])
        },
        logical(1L)
      )
      cand <- cand[feasible]
      if (length(cand) == 0L) {
        stop(
          sprintf(
            "no whole-cluster design satisfies the bounds for stratum '%s'",
            prep$stratum[h]
          ),
          call. = FALSE
        )
      }

      # Fixed cluster sizes are rounded to the nearest whole size whenever
      # that is compatible with the hard element bounds.  Only exact bounds
      # (notably take-all) may require a divisor farther away.
      if (isTRUE(prep$n_per_psu_fixed_h[h])) {
        near <- cand[cand %in% unique(c(floor(ps[h]), ceiling(ps[h])))]
        if (length(near) > 0L) cand <- near
      }
      score <- if (isTRUE(prep$has_stage_costs)) {
        (var_ratio[h] * (1 + icc[h] * (cand * resp_h[h] - 1))) *
          (prep$cost_psu_h[h] / cand + prep$cost_ssu_h[h])
      } else {
        abs(cand - ps[h])
      }
      list(take = as.integer(cand), score = score)
    }
  )
  b_h <- vapply(
    take_options,
    function(o) o$take[which.min(o$score)],
    integer(1L)
  )

  if (mode == "n") {
    # A take risen past the PSU universe fields every PSU, so it starts from
    # its nearest whole value rather than the cost-optimal one.
    risen <- !prep$n_per_psu_fixed_h & n_h >= psu_max * ps - 1e-9
    for (h in which(risen)) {
      o <- take_options[[h]]$take
      b_h[h] <- o[which.min(abs(o - ps[h]))]
    }
    # Raise a take only when capacity falls over half a take short of n.
    capacity <- function(h, b) min(floor(hi_i[h] / b), psu_max[h]) * b
    need <- round(sum(n_h))
    repeat {
      cap <- vapply(seq_len(H), function(h) capacity(h, b_h[h]), numeric(1L))
      if (need - sum(cap) <= max(b_h) / 2) break
      best <- NULL
      for (h in seq_len(H)) {
        o <- take_options[[h]]
        gain <- vapply(o$take, function(b) capacity(h, b), numeric(1L)) > cap[h]
        if (!any(gain)) next
        cur <- o$score[match(b_h[h], o$take)]
        pen <- o$score[gain] - cur
        if (isTRUE(prep$has_stage_costs)) pen <- pen / cur
        k <- which.min(pen)
        if (is.null(best) || pen[k] < best$pen) {
          best <- list(h = h, take = o$take[gain][k], pen = pen[k])
        }
      }
      if (is.null(best)) {
        stop(
          sprintf(
            "no whole-cluster design fields n = %d within half a take: whole takes within 'N' and 'N_1' hold at most %d units. Raise 'N_1' or lower 'n'",
            as.integer(need),
            as.integer(sum(cap))
          ),
          call. = FALSE
        )
      }
      b_h[best$h] <- best$take
    }
  }

  whole <- NULL
  searched_a <- NULL
  if (!is.null(max_nodes) && is.null(domain_search)) {
    whole <- .alloc_whole_cluster(prep, alloc, mode, target_cv, budget, deff,
                                  resp_rate, lo_i, hi_i, b_h, max_nodes, gap,
                                  target_se = target_se)
    b_h <- whole$b
  } else if (!is.null(max_nodes)) {
    # Each domain is searched as its own frame. A domain its floor holds
    # keeps the rounding below, which the floor repair raises.
    searched_a <- rep(NA_integer_, H)
    records <- list()
    for (d in seq_along(prep$domain_idx)) {
      nm <- names(prep$domain_idx)[d]
      if (domain_search$floored[d]) {
        records[[nm]] <- list(status = "heuristic", reason = "floor")
        next
      }
      idx <- prep$domain_idx[[d]]
      w <- .alloc_whole_cluster(
        .alloc_prep_slice(prep, idx), alloc, mode, target_cv, budget,
        .subset_h(deff, idx), .subset_h(resp_rate, idx), lo_i[idx],
        hi_i[idx], b_h[idx], max_nodes, gap, target_se = target_se
      )
      b_h[idx] <- w$b
      searched_a[idx] <- w$a
      records[[nm]] <- w$search %||% list(status = w$status)
    }
    domain_whole <- list(domains = records)
  }

  a_min <- as.integer(ceiling(lo_i / b_h))
  a_max <- floor(hi_i / b_h)
  # hi_i counts ultimate units, and the integer take b_h chosen here need not
  # be the one that set it, so the PSU universe has to bound a_max directly.
  if (!is.null(prep$N_psu_h)) {
    a_max <- pmin(a_max, prep$N_psu_h)
  }
  if (any(a_min > a_max)) {
    stop(
      "no whole-cluster design satisfies the integer allocation bounds",
      call. = FALSE
    )
  }

  e_target <- NULL
  if (mode == "n") {
    # Re-solve within what each stratum can field at its whole take.
    cap_lo <- a_min * b_h
    cap_hi <- a_max * b_h
    e_cont <- n_h
    outside <- n_h < cap_lo - 1e-9 | n_h > cap_hi + 1e-9
    if (any(outside) && !is.null(alloc_weights)) {
      total <- min(max(sum(n_h), sum(cap_lo)), sum(cap_hi))
      e_cont <- .rna_alloc(alloc_weights, total, cap_lo, cap_hi)
    }
    e_target <- .round_oric_bounded(
      e_cont, pmax(lo_i, cap_lo), pmin(hi_i, cap_hi),
      total = round(sum(e_cont))
    )
  }

  psu_cost <- if (isTRUE(prep$has_stage_costs)) {
    prep$cost_psu_h + prep$cost_ssu_h * b_h
  } else {
    rep(NA_real_, H)
  }
  # The whole-unit design fields its own take, so both halves of the
  # correction have to be rebuilt on 'b_h'. Reading the continuous take here
  # would score the operational design under a correction it does not get.
  op_take <- b_h * resp_h
  op_within_frac <- if (identical(prep$fpc %||% "unit", "stage")) {
    pmin(1, op_take / (.alloc_scope_N(prep) / prep$N_psu_h))
  } else {
    rep(0, H)
  }
  op_inflation <- var_ratio *
    (icc * op_take + (1 - op_within_frac) * (1 - icc))
  S_op <- S_raw * sqrt(op_inflation)
  N_fpc_op <- switch(
    prep$fpc %||% "unit",
    unit = .alloc_scope_N(prep),
    none = rep(Inf, H),
    stage = ifelse(
      icc > 0,
      prep$N_psu_h * op_inflation / (var_ratio * icc),
      Inf
    )
  )
  Cj <- W^2 * S_op^2 * deff / (b_h * resp_rate)

  if (!is.null(whole)) {
    a_h <- whole$a
  } else if (mode == "cv") {
    a_h <- as.integer(ceiling(
      n_h *
        (1 + icc * (b_h * resp_h - 1)) /
        (1 + icc * (ps * resp_h - 1)) /
        b_h -
        1e-9
    ))
    a_h <- pmin(pmax(a_h, a_min), a_max)
    # The rescaling reads the uncorrected clustering factor, and rounding the
    # take moves both stage fractions under fpc = "stage", so it can stop
    # short. Add the PSU that lowers the cv most per unit of cost until the
    # target is met or every stratum is at its bound.
    goal_op <- if (!is.null(target_se)) target_se else target_cv
    metric_at <- function(a) {
      .alloc_op_metric(prep, a * b_h, S_op, N_fpc_op, alpha, deff, resp_rate,
                       df, metric = if (!is.null(target_se)) "se" else "cv")
    }
    per_psu <- if (all(is.finite(psu_cost))) psu_cost else rep(1, H)
    repeat {
      if (is.null(goal_op)) break
      now <- metric_at(a_h)
      if (is.finite(now) && now <= goal_op + 1e-8) break
      cand <- which(a_h < a_max)
      if (length(cand) == 0L) break
      after <- vapply(cand, function(h) {
        a_try <- a_h
        a_try[h] <- a_try[h] + 1L
        metric_at(a_try)
      }, numeric(1L))
      j <- cand[which.max((now - after) / per_psu[cand])]
      a_h[j] <- a_h[j] + 1L
    }
  } else if (mode == "budget") {
    min_cost <- sum(a_min * psu_cost)
    if (budget < min_cost * (1 - 1e-8)) {
      stop(
        sprintf(
          "'budget' cannot fund one whole PSU per stratum while satisfying the lower bounds (minimum field cost = %.4g)",
          min_cost
        ),
        call. = FALSE
      )
    }
    room <- .budget_room(budget, min_cost)
    a_h <- pmin(pmax(as.integer(floor(n_h / b_h + 1e-9)), a_min), a_max)
    repeat {
      if (sum(a_h * psu_cost) <= room) {
        break
      }
      cand <- which(a_h > a_min)
      if (length(cand) == 0L) {
        stop(
          "no whole-cluster design satisfies both the budget and lower bounds",
          call. = FALSE
        )
      }
      loss <- Cj[cand] * (1 / (a_h[cand] - 1L) - 1 / a_h[cand]) / psu_cost[cand]
      j <- cand[which.min(loss)]
      a_h[j] <- a_h[j] - 1L
    }
    repeat {
      cand <- which(psu_cost <= room - sum(a_h * psu_cost) & a_h < a_max)
      if (length(cand) == 0L) {
        break
      }
      gain <- Cj[cand] * (1 / a_h[cand] - 1 / (a_h[cand] + 1L)) / psu_cost[cand]
      j <- cand[which.max(gain)]
      a_h[j] <- a_h[j] + 1L
    }
  } else {
    a_h <- as.integer(round(e_target / b_h))
    a_h <- pmin(pmax(a_h, a_min), a_max)
    # A take chosen to divide the total fields a prime total as 1-unit PSUs.
    a_h <- as.integer(.cluster_close_total(a_h, b_h, a_min, a_max, e_target))
    field <- sum(a_h * b_h)
    if (field - need > max(b_h) / 2) {
      smallest <- sum(a_min * b_h)
      stop(
        sprintf(
          "no whole-cluster design fields n = %d within half a take: at whole takes, with at least one PSU per stratum and the strata's lower bounds, the smallest design holds %d units. Raise 'n' to %d, or relax 'min_n_stratum', 'max_weight' or the take",
          as.integer(need),
          as.integer(smallest),
          as.integer(smallest)
        ),
        call. = FALSE
      )
    }
    if (need - field > max(b_h) / 2) {
      stop(
        sprintf(
          "no whole-cluster design fields n = %d within half a take: whole takes within 'N' and 'N_1' hold at most %d units. Raise 'N_1' or lower 'n'",
          as.integer(need),
          as.integer(sum(a_max * b_h))
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(searched_a)) {
    held <- !is.na(searched_a)
    a_h[held] <- searched_a[held]
    a_h <- as.integer(a_h)
  }
  if (!is.null(floors)) {
    before <- a_h
    a_h <- .cluster_floor_repair(a_h, b_h, n_h, a_min, a_max, floors, mode,
                                 isTRUE(prep$has_stage_costs), psu_cost,
                                 budget, e_target)
    if (!is.null(searched_a)) {
      domain_whole$domains <- .whole_repaired(domain_whole$domains,
                                              prep$domain_idx, before, a_h)
    }
  }
  if (!is.null(searched_a)) {
    domain_whole <- .whole_domains_summary(domain_whole$domains)
  }

  e_h <- a_h * b_h
  m <- .alloc_metrics(
    N_h = prep$N_h,
    N_pop = .alloc_scope_N(prep),
    S_h = S_op,
    mean_h = prep$mean_h,
    n_h = e_h,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = prep$cost_h,
    df = df,
    N_fpc = N_fpc_op
  )
  detail$n_int <- as.integer(e_h)
  detail$n_1_int <- a_h
  detail$take_2_int <- b_h

  # A one-PSU stratum leaves the design with no within-stratum variance
  # estimate at all, which the caller should hear about here rather than
  # discovering it from design_df() later.
  .warn_singleton_strata(prep$stratum[a_h == 1L], "psu")

  if (any(e_h < lo_i | e_h > hi_i)) {
    stop(
      "internal error: operational cluster allocation violates its bounds",
      call. = FALSE
    )
  }
  if (!is.null(target_se)) {
    op_se <- .alloc_op_metric(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate,
                              df, metric = "se")
    if (!is.finite(op_se) || op_se > target_se + 1e-8) {
      stop(
        sprintf(
          "no whole-cluster design found that attains target precision under the integer bounds%s",
          .alloc_psu_bound_note(prep)
        ),
        call. = FALSE
      )
    }
    op_cv <- m$cv
  } else if (!is.null(target_cv)) {
    op_cv <- .alloc_op_cv(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate,
                          df)
    if (!is.finite(op_cv) || op_cv > target_cv + 1e-8) {
      stop(
        sprintf(
          "no whole-cluster design found that attains target CV %.4g under the integer bounds%s",
          target_cv,
          .alloc_psu_bound_note(prep)
        ),
        call. = FALSE
      )
    }
  } else {
    op_cv <- m$cv
  }
  operational <- list(
    n = sum(e_h),
    cost = if (isTRUE(prep$has_stage_costs)) {
      sum(a_h * psu_cost)
    } else {
      NA_real_
    },
    se = m$se,
    moe = m$moe,
    cv = m$cv,
    rmoe = m$rmoe
  )
  if (!is.null(searched_a)) whole <- domain_whole
  operational$status <- whole$status
  list(operational = operational, detail = detail, search = whole$search,
       domains = whole$domains)
}

#' Whole-unit field design for a stratified S-stage allocation (S >= 3)
#'
#' Searched modes run the exact search, the others round greedily.
#' @keywords internal
#' @noRd
.alloc_operational_cluster_multistage <- function(
  prep,
  detail,
  n_h,
  mode,
  target_cv = NULL,
  target_se = NULL,
  budget = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  df = NULL,
  alloc_weights = NULL,
  floors = NULL,
  alloc = NULL,
  max_nodes = NULL,
  gap = 0,
  domain_search = NULL
) {
  H <- length(n_h)
  S <- prep$stages
  lo_i <- ceiling(detail$.lower - 1e-9)
  hi_i <- floor(detail$.upper + 1e-9)
  psu_capped <- prep$M_src %in% "N_1"
  hi_i[psu_capped] <- floor(prep$N_h[psu_capped] + 1e-9)
  hi_i <- pmin(hi_i, prep$whole_cap_h %||% Inf)

  b_mat <- .ms_round_takes(prep, hi_i)
  searched <- !is.null(max_nodes) && mode %in% c("cv", "budget")
  ## A searched budget may fund smaller takes than the rounded ones.
  g <- tryCatch(
    .ms_greedy(prep, b_mat, n_h, mode, lo_i, hi_i, target_cv, target_se,
               budget, alpha, deff, resp_rate, df, alloc_weights),
    error = function(e) {
      if (searched || identical(mode, "n")) return(NULL)
      stop(e)
    }
  )
  a_h <- g$a %||% rep(NA_integer_, H)
  selected <- FALSE
  shift <- NULL
  if (identical(mode, "n")) {
    shift <- .ms_exact_total(prep, b_mat, a_h, n_h, lo_i, hi_i, detail$n_1,
                             floors)
    if (is.null(shift$a)) {
      if (anyNA(a_h)) stop(shift$message, call. = FALSE)
    } else {
      a_h <- .as_count(shift$a, "the whole-unit rounding of n_alloc()",
                       "at the PSU stage")
      b_mat <- shift$b
      selected <- TRUE
    }
  }
  status <- NULL
  search <- NULL
  domains <- NULL
  if (searched) {
    groups <- if (is.null(domain_search)) {
      list(seq_len(H))
    } else {
      prep$domain_idx[!domain_search$floored]
    }
    records <- list()
    if (!is.null(domain_search)) {
      for (nm in names(prep$domain_idx)[domain_search$floored]) {
        records[[nm]] <- list(status = "heuristic", reason = "floor")
      }
    }
    heuristic <- character(0)
    for (k in seq_along(groups)) {
      idx <- groups[[k]]
      sub <- if (is.null(domain_search)) prep else .alloc_prep_slice(prep, idx)
      ## The takes-first route fields the design its first search proved.
      w <- prep$given_whole
      if (is.null(w)) {
        inc <- if (!anyNA(a_h[idx])) {
          .ms_incumbent(sub, alloc, b_mat[idx, , drop = FALSE], a_h[idx], mode,
                        lo_i[idx], hi_i[idx], target_cv, target_se, budget,
                        alpha, .subset_h(deff, idx), .subset_h(resp_rate, idx),
                        df, alloc_weights[idx], n_h[idx])
        }
        w <- .alloc_whole_cluster_multistage(
          sub, alloc, mode, target_cv, budget, .subset_h(deff, idx),
          .subset_h(resp_rate, idx), lo_i[idx], hi_i[idx],
          b_mat[idx, , drop = FALSE], max_nodes, gap, target_se = target_se,
          incumbent = inc
        )
      } else {
        range <- .ms_psu_range(prep, apply(matrix(w$b, nrow = H), 1L, prod),
                               lo_i, hi_i)
        if (any(w$a < range$lo | w$a > range$hi)) {
          stop("internal error: the searched design left its bounds at its own takes",
               call. = FALSE)
        }
      }
      if (identical(w$status, "heuristic")) {
        if (!is.null(w$a)) {
          a_h[idx] <- w$a
          b_mat[idx, ] <- w$b
        }
        if (anyNA(a_h[idx])) {
          stop(
            sprintf(
              "the whole design was not searched, because %s, and the greedy rounding does not meet the bounds",
              .ms_reason_text(w$reason)
            ),
            call. = FALSE
          )
        }
        heuristic <- c(heuristic, stats::setNames(
          w$reason, paste(sQuote(prep$stratum[idx]), collapse = ", ")
        ))
        records[[if (is.null(domain_search)) "all" else names(groups)[k]]] <-
          list(status = "heuristic", reason = w$reason)
        next
      }
      if (identical(w$status, "infeasible")) {
        stop(
          if (mode == "budget") {
            sprintf(
              "'budget' cannot fund a whole-cluster design within the lower bounds at any whole takes%s",
              if (is.null(w$min_cost) || is.na(w$min_cost)) "" else {
                sprintf(" (minimum whole-design cost = %s)",
                        format(w$min_cost, digits = 7, scientific = FALSE))
              }
            )
          } else {
            sprintf(
              "no whole-cluster design found that attains the target precision under the integer bounds%s",
              .alloc_psu_bound_note(prep)
            )
          },
          call. = FALSE
        )
      }
      if (is.null(w$a)) {
        stop(structure(
          list(message = sprintf(
            "no whole-cluster design was found in %d search node%s, so whether one exists is unknown. Raise 'control$max_nodes'",
            as.integer(max_nodes), if (max_nodes == 1) "" else "s"
          ), call = NULL),
          class = c("svyplan_plan_unknown", "error", "condition")
        ))
      }
      a_h[idx] <- w$a
      b_mat[idx, ] <- w$b
      records[[if (is.null(domain_search)) "all" else names(groups)[k]]] <-
        w$search %||% list(status = w$status)
    }
    for (why in unique(heuristic)) {
      warning(
        sprintf(
          "the whole design of %s was not searched to completion: %s, so the design found before the search is used",
          paste(names(heuristic)[heuristic == why], collapse = ", "),
          .ms_reason_text(why)
        ),
        call. = FALSE
      )
    }
  }
  P_h <- apply(b_mat, 1L, prod)
  what <- "the whole-unit rounding of n_alloc()"
  range <- .ms_psu_range(prep, P_h, lo_i, hi_i)
  a_min <- .as_count(range$lo, what, "at the PSU stage")
  a_max <- range$hi
  op_terms <- .cluster_stage_terms(prep, b_mat)
  S_op <- prep$S_raw_h * sqrt(op_terms$inflation)
  N_fpc_op <- op_terms$N_fpc_h
  psu_cost <- op_terms$psu_cost
  elem_cost <- op_terms$elem_cost %||% prep$cost_h

  if (!is.null(floors) && !selected) {
    before <- a_h
    a_h <- tryCatch(
      .cluster_floor_repair(a_h, P_h, n_h, a_min, a_max, floors, mode,
                            isTRUE(prep$has_stage_costs), psu_cost, budget,
                            g$e_target),
      error = function(e) {
        if (is.null(shift$message)) stop(e)
        stop(shift$message, call. = FALSE)
      }
    )
    if (searched && !is.null(domain_search)) {
      records <- .whole_repaired(records, prep$domain_idx, before, a_h)
    }
  }
  short <- !is.null(floors) && any(
    vapply(seq_along(floors$need), function(d) sum((a_h * P_h)[floors$dom == d]),
           1) < ceiling(floors$need - 1e-9)
  )
  if (anyNA(a_h) || any(a_h < a_min | a_h > a_max) || short) {
    stop("internal error: the whole design left its bounds or domain floors",
         call. = FALSE)
  }
  if (searched && length(records)) {
    summ <- .whole_domains_summary(records)
    status <- summ$status
    search <- summ$search
    if (!is.null(domain_search)) domains <- records
  }

  e_h <- .as_count(a_h * P_h, what, "at the final stage")
  m <- .alloc_metrics(
    N_h = prep$N_h,
    N_pop = .alloc_scope_N(prep),
    S_h = S_op,
    mean_h = prep$mean_h,
    n_h = e_h,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = elem_cost,
    df = df,
    N_fpc = N_fpc_op
  )
  detail$n_int <- e_h
  detail$n_1_int <- a_h
  for (s in 2L:S) {
    detail[[paste0("take_", s, "_int")]] <- b_mat[, s - 1L]
  }
  .warn_singleton_strata(prep$stratum[a_h == 1L], "psu")

  if (!is.null(target_se)) {
    op_se <- .alloc_op_metric(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate,
                              df, metric = "se")
    if (!is.finite(op_se) || op_se > target_se + 1e-8) {
      stop(
        sprintf(
          "no whole-cluster design found that attains target precision under the integer bounds%s",
          .alloc_psu_bound_note(prep)
        ),
        call. = FALSE
      )
    }
    op_cv <- m$cv
  } else if (!is.null(target_cv)) {
    op_cv <- .alloc_op_cv(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate, df)
    if (!is.finite(op_cv) || op_cv > target_cv + 1e-8) {
      stop(
        sprintf(
          "no whole-cluster design found that attains target CV %.4g under the integer bounds%s",
          target_cv,
          .alloc_psu_bound_note(prep)
        ),
        call. = FALSE
      )
    }
  } else {
    op_cv <- m$cv
  }

  operational <- list(
    n = sum(e_h),
    cost = if (isTRUE(prep$has_stage_costs)) sum(a_h * psu_cost) else NA_real_,
    se = m$se,
    moe = m$moe,
    cv = m$cv,
    rmoe = m$rmoe,
    status = status
  )
  list(operational = operational, detail = detail, search = search,
       domains = domains)
}

#' Refuse a free take whose response floor its stage size cannot hold
#'
#' The continuous take is capped under the stage correction only, the whole
#' take at three or more stages under every correction.
#' @keywords internal
#' @noRd
.cluster_take_room <- function(stratum, N_h, N_st, rates, fixed, rows, fpc) {
  S <- ncol(rates) + 1L
  for (h in which(rows)) {
    caps <- .stage_caps_of(N_h[h], N_st, h, S)
    cont_cap <- caps
    whole_cap <- caps
    for (s in which(!fixed[h, ])) {
      in_whole <- .take_floor(rates[h, s]) > .whole_caps(whole_cap[s])
      if (.take_bound(rates[h, s]) > cont_cap[s] + 1e-12 || in_whole) {
        stop(
          sprintf(
            "stage %d of stratum '%s' holds %s units per parent unit, which cannot deliver one expected response %s",
            s + 1L, stratum[h],
            format(if (in_whole) whole_cap[s] else cont_cap[s], digits = 15),
            if (.take_bound(rates[h, s]) > cont_cap[s] + 1e-12) {
              "at this response rate"
            } else {
              "in whole units"
            }
          ),
          call. = FALSE
        )
      }
    }
  }
  invisible(NULL)
}

#' Why a classic multistage whole design was not searched
#' @keywords internal
#' @noRd
.ms_reason_text <- function(reason) {
  switch(
    reason,
    listing = "its take vectors are too many to list",
    frontier = "its cost and variance frontier has too many states",
    numerical = "the target sits at the variance floor and its designs there are too many to list",
    pruning = "numerical comparisons could not preserve the feasible starting design"
  )
}

#' Search the whole design first when its continuous takes are refused
#' @keywords internal
#' @noRd
.alloc_takes_first <- function(frame, call_args, prep, alloc, mode, cv, budget,
                               deff, resp_rate, lo_i, hi_i, control, target_se,
                               refusal) {
  psu_capped <- prep$M_src %in% "N_1"
  hi_i[psu_capped] <- floor(prep$N_h[psu_capped] + 1e-9)
  sc <- .search_control(control)
  w <- .alloc_whole_cluster_multistage(
    prep, alloc, mode, cv, budget, deff, resp_rate, lo_i, hi_i,
    .ms_round_takes(prep, hi_i), sc$max_nodes, sc$gap, target_se = target_se
  )
  if (identical(w$status, "infeasible")) stop(refusal, call. = FALSE)
  if (is.null(w$a)) {
    stop(structure(
      list(message = sprintf(
        "%s at the continuous takes, and the whole search over the other takes found no design%s, so whether one exists is unknown. Raise 'control$max_nodes'",
        refusal,
        if (identical(w$status, "heuristic")) {
          sprintf(" because %s", .ms_reason_text(w$reason))
        } else {
          sprintf(" in %s nodes", format(sc$max_nodes, big.mark = ",",
                                          scientific = FALSE))
        }
      ), call = NULL),
      class = c("svyplan_plan_unknown", "error", "condition")
    ))
  }
  takes <- matrix(w$b, nrow = nrow(frame))
  colnames(takes) <- paste0("take_", seq_len(ncol(takes)) + 1L)
  fixed <- .alloc_frame_at_takes(frame, takes)
  attr(fixed, "svyplan_whole_takes") <- w
  fit <- do.call(n_alloc.default, c(list(frame = fixed), call_args))
  fit$params$frame <- frame
  fit$params$searched_takes <- takes
  fit
}

#' Whole designs near each stratum's share, chosen to add up to the total within the domain floors
#'
#' `NULL` designs and a `message` when no choice meets the floors.
#' @keywords internal
#' @noRd
.ms_exact_total <- function(prep, b_mat, a_h, n_h, lo_i, hi_i, a_cont,
                            floors = NULL) {
  H <- length(n_h)
  goal <- round(sum(n_h))
  P_h <- apply(b_mat, 1L, prod)
  dom <- floors$dom %||% rep(1L, H)
  need <- if (is.null(floors)) 0 else ceiling(floors$need - 1e-9)
  groups <- split(seq_len(H), factor(dom, levels = seq_along(need)))
  meets <- function(e) all(vapply(groups, function(i) sum(e[i]), 1) >= need)
  now <- if (!anyNA(a_h) && meets(a_h * P_h)) sum(a_h * P_h) else NA
  if (!is.na(now) && now == goal) return(list(a = a_h, b = b_mat))
  t_h <- .round_oric_bounded(n_h, lo_i, hi_i, total = goal)
  free <- !prep$takes_fixed_mat
  take_floor <- .cluster_take_floor(prep)
  range <- .ms_psu_range(prep, P_h, lo_i, hi_i)
  opts <- lapply(seq_len(H), function(h) {
    b <- b_mat[h, ]
    ## Near the rounded target and the continuous share, which differ where
    ## the rounding moved a shortfall in whole capacity.
    mid <- unique(c(t_h[h], round(n_h[h])))
    w <- max(3, min(20, ceiling(.05 * t_h[h])))
    totals <- unique(unlist(lapply(mid, function(m) {
      from <- max(lo_i[h], m - w)
      to <- min(hi_i[h], m + w)
      if (from <= to) seq(from, to)
    })))
    a <- if (!is.na(a_h[h])) a_h[h]
    legal <- all(b >= ifelse(free[h, ], take_floor[h, ], 1)) &&
      all(b <= .whole_caps(.cluster_take_caps(prep, h)))
    if (legal && range$lo[h] <= range$hi[h]) {
      P <- P_h[h]
      fit <- totals[totals %% P == 0] / P
      near <- pmin(pmax(c(floor(mid / P), ceiling(mid / P)), range$lo[h]),
                   range$hi[h])
      a <- unique(c(a, fit[fit >= range$lo[h] & fit <= range$hi[h]], near))
    }
    cand <- lapply(a, function(x) list(a = x, b = b, rank = 0))
    if (any(free[h, ])) {
      for (T in sort(totals)) {
        x <- .ms_exact_takes(prep, h, T, b, a_cont[h])
        if (!is.null(x)) {
          cand[[length(cand) + 1L]] <- list(a = x$a, b = x$b, rank = 1)
        }
      }
    }
    if (!length(cand)) return(NULL)
    total <- vapply(cand, function(o) o$a * prod(o$b), 1)
    rank <- vapply(cand, `[[`, 1, "rank")
    score <- vapply(cand, function(o) {
      log(o$a / a_cont[h])^2 + sum(log(o$b / prep$takes_mat[h, ])^2)
    }, 1)
    o <- order(total, rank, score)
    o <- o[!duplicated(total[o])]
    list(total = total[o], cost = (total[o] - n_h[h])^2, design = cand[o])
  })
  empty <- vapply(opts, is.null, TRUE)
  if (any(empty)) {
    return(list(message = sprintf(
      "no whole design within the stage sizes and bounds of stratum %s fields a total near its share of 'n'",
      paste(sQuote(prep$stratum[empty]), collapse = ", ")
    )))
  }
  ## Each total keeps its closest choice; ties go to the earlier option.
  reach <- function(items) {
    tot <- 0
    cost <- 0
    back <- vector("list", length(items))
    for (j in seq_along(items)) {
      S <- length(tot)
      K <- length(items[[j]]$total)
      tt <- rep(tot, times = K) + rep(items[[j]]$total, each = S)
      cc <- rep(cost, times = K) + rep(items[[j]]$cost, each = S)
      from <- rep(seq_len(S), times = K)
      k <- rep(seq_len(K), each = S)
      o <- order(tt, cc, from, k)
      o <- o[!duplicated(tt[o])]
      back[[j]] <- list(from = from[o], k = k[o])
      tot <- tt[o]
      cost <- cc[o]
    }
    list(total = tot, cost = cost, back = back)
  }
  trace <- function(res, i) {
    picks <- integer(length(res$back))
    for (j in rev(seq_along(res$back))) {
      picks[j] <- res$back[[j]]$k[i]
      i <- res$back[[j]]$from[i]
    }
    picks
  }
  inner <- lapply(seq_along(groups), function(d) {
    res <- reach(opts[groups[[d]]])
    res$ok <- which(res$total >= need[d])
    res
  })
  short <- vapply(inner, function(r) length(r$ok) == 0L, TRUE)
  if (any(short)) {
    return(list(message = sprintf(
      "no whole design among the totals near each stratum's share of 'n' meets the 'min_n_domain' floor of %s",
      paste(sQuote(floors$levels[short]), collapse = ", ")
    )))
  }
  outer <- reach(lapply(inner, function(r) {
    list(total = r$total[r$ok], cost = r$cost[r$ok])
  }))
  i <- order(abs(outer$total - goal), outer$cost)[1L]
  if (!is.na(now) && abs(outer$total[i] - goal) >= abs(now - goal)) {
    return(list(a = a_h, b = b_mat))
  }
  pick_d <- trace(outer, i)
  for (d in seq_along(groups)) {
    rows <- groups[[d]]
    pk <- trace(inner[[d]], inner[[d]]$ok[pick_d[d]])
    for (j in seq_along(rows)) {
      o <- opts[[rows[j]]]$design[[pk[j]]]
      a_h[rows[j]] <- o$a
      b_mat[rows[j], ] <- o$b
    }
  }
  list(a = a_h, b = b_mat)
}

#' The divisors of a positive whole number, in increasing order
#' @keywords internal
#' @noRd
.divisors <- function(x) {
  d <- seq_len(floor(sqrt(x) + 1e-9))
  d <- d[x %% d == 0]
  sort(unique(c(d, x / d)))
}

#' Every split of `P` into whole takes within `[lo_s, hi_s]` (the first with `first = TRUE`), or `NULL`
#' @keywords internal
#' @noRd
.ms_take_split <- function(P, lo_s, hi_s, first = TRUE) {
  K <- length(lo_s)
  out <- list()
  walk <- function(P, k, pre) {
    if (k == K) {
      if (P >= lo_s[k] && P <= hi_s[k]) out[[length(out) + 1L]] <<- c(pre, P)
      return(invisible())
    }
    for (d in .divisors(P)) {
      if (d < lo_s[k]) next
      if (d > hi_s[k]) break
      walk(P / d, k + 1L, c(pre, d))
      if (first && length(out)) return(invisible())
    }
  }
  walk(P, 1L, numeric(0))
  if (!length(out)) return(NULL)
  if (first) out[[1L]] else out
}

#' The whole PSU count and takes of one stratum that field `target` exactly, near its continuous design
#' @keywords internal
#' @noRd
.ms_exact_takes <- function(prep, h, target, b_row, a_cont) {
  if (target < 1 || target != round(target)) return(NULL)
  free <- !prep$takes_fixed_mat[h, ]
  lo_s <- ifelse(free, .cluster_take_floor(prep)[h, ], b_row)
  hi_s <- ifelse(free, .whole_caps(.cluster_take_caps(prep, h)), b_row)
  N_1 <- prep$N_psu_h[h] %||% Inf
  most <- prep$N_h[h] - if (is.finite(N_1)) N_1 - 1 else 0
  t_cont <- prep$takes_mat[h, ]
  near <- a_cont >= 1
  a_lo <- if (near) a_cont / 2 else 1
  a_hi <- if (near) 2 * a_cont else Inf
  if (near) {
    lo_s <- pmax(lo_s, ifelse(free, ceiling(t_cont / 2 - 1e-9), lo_s))
    hi_s <- pmin(hi_s, ifelse(free, floor(2 * t_cont + 1e-9), hi_s))
    if (any(lo_s > hi_s)) return(NULL)
  }
  best <- NULL
  for (P in .divisors(target)) {
    a <- target / P
    if (P > most || a > N_1 || a < a_lo - 1e-9 || a > a_hi + 1e-9 ||
          a > .psu_room(prep$N_h[h], prep$N_psu_h[h], P)) {
      next
    }
    for (b in .ms_take_split(P, lo_s, hi_s, first = FALSE)) {
      score <- log(a / a_cont)^2 + sum(log(b / t_cont)^2)
      if (is.null(best) || score < best$score) {
        best <- list(a = a, b = b, score = score)
      }
    }
  }
  best
}

#' Round a stratum's continuous takes to whole ones within its stage sizes
#' @keywords internal
#' @noRd
.ms_round_takes <- function(prep, hi_i) {
  H <- nrow(prep$takes_mat)
  S <- prep$stages
  b_mat <- matrix(1L, nrow = H, ncol = S - 1L)
  take_floor <- .cluster_take_floor(prep)
  for (h in seq_len(H)) {
    rounded <- .cluster_round_takes(prep$takes_mat[h, ],
                                    .cluster_take_caps(prep, h), take_floor[h, ])
    for (s in 2L:S) {
      b_mat[h, s - 1L] <- .as_count(rounded[s - 1L],
                                   "the whole-unit rounding of n_alloc()",
                                   sprintf("in take_%d", s))
    }
    most <- hi_i[h]
    if (!is.null(prep$N_psu_h)) most <- min(most, prep$N_h[h] - prep$N_psu_h[h] + 1)
    while (prod(b_mat[h, ]) > most) {
      idx <- which(b_mat[h, ] > take_floor[h, ] & !prep$takes_fixed_mat[h, ])
      if (length(idx) == 0L) break
      k <- idx[length(idx)]
      b_mat[h, k] <- b_mat[h, k] - 1L
    }
  }
  b_mat
}

#' PSU counts of the greedy whole design at given whole takes
#' @keywords internal
#' @noRd
.ms_greedy <- function(prep, b_mat, n_h, mode, lo_i, hi_i, target_cv,
                       target_se, budget, alpha, deff, resp_rate, df,
                       alloc_weights) {
  H <- length(n_h)
  P_h <- apply(b_mat, 1L, prod)
  what <- "the whole-unit rounding of n_alloc()"
  range <- .ms_psu_range(prep, P_h, lo_i, hi_i)
  a_min <- .as_count(range$lo, what, "at the PSU stage")
  a_max <- range$hi
  if (any(a_min > a_max)) {
    stop(
      "no whole-cluster design satisfies the integer allocation bounds",
      call. = FALSE
    )
  }
  op_terms <- .cluster_stage_terms(prep, b_mat)
  S_op <- prep$S_raw_h * sqrt(op_terms$inflation)
  N_fpc_op <- op_terms$N_fpc_h
  psu_cost <- op_terms$psu_cost
  W <- .alloc_scope_N(prep) / sum(.alloc_scope_N(prep))
  Cj <- W^2 * S_op^2 * deff / (P_h * resp_rate)
  ok <- TRUE
  if (mode == "cv") {
    a_h <- .as_count(ceiling(n_h / P_h - 1e-9), what, "at the PSU stage")
    a_h <- pmin(pmax(a_h, a_min), a_max)
    goal_op <- if (!is.null(target_se)) target_se else target_cv
    metric_at <- function(a) {
      .alloc_op_metric(
        prep, a * P_h, S_op, N_fpc_op, alpha, deff, resp_rate,
        df, metric = if (!is.null(target_se)) "se" else "cv"
      )
    }
    per_psu <- if (all(is.finite(psu_cost))) psu_cost else rep(1, H)
    repeat {
      if (is.null(goal_op)) break
      now <- metric_at(a_h)
      if (is.finite(now) && now <= goal_op + 1e-8) break
      cand <- which(a_h < a_max)
      if (length(cand) == 0L) {
        ok <- FALSE
        break
      }
      after <- vapply(cand, function(h) {
        a_try <- a_h
        a_try[h] <- a_try[h] + 1L
        metric_at(a_try)
      }, numeric(1L))
      j <- cand[which.max((now - after) / per_psu[cand])]
      a_h[j] <- a_h[j] + 1L
    }
  } else if (mode == "budget") {
    min_cost <- sum(a_min * psu_cost)
    if (budget < min_cost * (1 - 1e-8)) {
      stop(
        sprintf(
          "'budget' cannot fund one whole PSU per stratum while satisfying the lower bounds (minimum field cost = %.4g)",
          min_cost
        ),
        call. = FALSE
      )
    }
    room <- .budget_room(budget, min_cost)
    a_h <- .as_count(floor(n_h / P_h + 1e-9), what, "at the PSU stage")
    a_h <- pmin(pmax(a_h, a_min), a_max)
    repeat {
      if (sum(a_h * psu_cost) <= room) break
      cand <- which(a_h > a_min)
      if (length(cand) == 0L) {
        stop(
          "no whole-cluster design satisfies both the budget and lower bounds",
          call. = FALSE
        )
      }
      loss <- Cj[cand] * (1 / (a_h[cand] - 1L) - 1 / a_h[cand]) / psu_cost[cand]
      j <- cand[which.min(loss)]
      a_h[j] <- a_h[j] - 1L
    }
    repeat {
      cand <- which(psu_cost <= room - sum(a_h * psu_cost) & a_h < a_max)
      if (length(cand) == 0L) break
      gain <- Cj[cand] * (1 / a_h[cand] - 1 / (a_h[cand] + 1L)) / psu_cost[cand]
      j <- cand[which.max(gain)]
      a_h[j] <- a_h[j] + 1L
    }
  } else {
    cap_lo <- a_min * P_h
    cap_hi <- a_max * P_h
    e_cont <- n_h
    outside <- n_h < cap_lo - 1e-9 | n_h > cap_hi + 1e-9
    if (any(outside) && !is.null(alloc_weights)) {
      total <- min(max(sum(n_h), sum(cap_lo)), sum(cap_hi))
      e_cont <- .rna_alloc(alloc_weights, total, cap_lo, cap_hi)
    }
    e_target <- .round_oric_bounded(
      e_cont, pmax(lo_i, cap_lo), pmin(hi_i, cap_hi),
      total = round(sum(e_cont))
    )
    a_h <- .as_count(round(e_target / P_h), what, "at the PSU stage")
    a_h <- pmin(pmax(a_h, a_min), a_max)
    a_h <- .as_count(.cluster_close_total(a_h, P_h, a_min, a_max, e_target),
                     what, "at the PSU stage")
  }
  list(a = a_h, ok = ok, cost = sum(a_h * psu_cost),
       e_target = if (mode == "n") e_target)
}

#' Raise a cluster design to its domain floors, keeping its total or budget
#' @keywords internal
#' @noRd
.cluster_floor_repair <- function(a_h, unit, n_h, a_min, a_max, floors, mode,
                                  costed, psu_cost, budget, e_target = NULL) {
  need <- ceiling(floors$need - 1e-9)
  total <- if (mode == "n" && !is.null(e_target)) sum(e_target)
  keep <- if (mode == "budget" && costed) {
    "budget"
  } else if (!is.null(total)) {
    "total"
  } else {
    "none"
  }
  a <- .alloc_floor_repair(
    a_h, cont = n_h / unit, unit = unit, lo = a_min, hi = a_max,
    dom = floors$dom, need = need, keep = keep, total = total,
    cost_step = psu_cost, budget = budget
  )
  if (identical(keep, "total") && sum(a * unit) != total) {
    closed <- .floor_close_total(a, unit, a_min, a_max, floors$dom, need,
                                 e_target, total)
    if (length(closed) == 0L) {
      reached <- attr(closed, "totals", exact = TRUE)
      best <- reached[which.min(abs(reached - total))]
      closed <- if (length(best) &&
                      abs(best - total) < abs(sum(a * unit) - total)) {
        .floor_close_total(a, unit, a_min, a_max, floors$dom, need, e_target,
                           best)
      } else {
        a
      }
    }
    a <- closed
  }
  as.integer(a)
}

#' A feasible whole design to bound and start the multistage search
#'
#' `NULL` when no candidate meets the target.
#' @keywords internal
#' @noRd
.ms_incumbent <- function(prep, alloc, b_mat, a_h, mode, lo_i, hi_i, target_cv,
                          target_se, budget, alpha, deff, resp_rate, df,
                          alloc_weights, n_h) {
  H <- nrow(b_mat)
  costed <- alloc == "optimal" && isTRUE(prep$has_stage_costs)
  free <- costed & !prep$takes_fixed_mat
  tries <- list(list(a = a_h, b = b_mat))
  if (any(free)) {
    caps <- t(vapply(seq_len(H), function(h) .cluster_take_caps(prep, h),
                     numeric(ncol(b_mat))))
    caps <- .whole_caps(caps)
    up <- b_mat
    up[free] <- pmin(caps, pmax(b_mat, ceiling(prep$takes_mat - 1e-9)))[free]
    big <- up
    at_cap <- free & is.finite(caps)
    big[at_cap] <- caps[at_cap]
    for (b in list(up, big)) {
      if (any(apply(b, 1L, prod) > hi_i)) next
      g <- tryCatch(
        .ms_greedy(prep, b, n_h, mode, lo_i, hi_i, target_cv, target_se,
                   budget, alpha, deff, resp_rate, df, alloc_weights),
        error = function(e) NULL
      )
      if (!is.null(g)) tries[[length(tries) + 1L]] <- list(a = g$a, b = b)
    }
  }
  W <- .alloc_scope_N(prep) / sum(.alloc_scope_N(prep))
  r <- rep_len(resp_rate, H)
  w2s2d <- W^2 * prep$S_raw_h^2 * rep_len(deff, H)
  goal <- if (mode == "cv") {
    se <- target_se %||% (target_cv * abs(.aggregate_mean(W, prep$mean_h)))
    se^2 * (1 + 1e-12)
  }
  best <- NULL
  best_value <- Inf
  for (tr in tries) {
    terms <- lapply(seq_len(H), function(h) {
      .ms_state_terms(prep, h, tr$b[h, , drop = FALSE], lo_i[h], hi_i[h],
                      costed, w2s2d[h], r[h])
    })
    col <- function(name) vapply(terms, `[[`, 1, name)
    if (any(!is.finite(tr$a) | tr$a < col("lo") | tr$a > col("hi"))) next
    pc <- col("pc")
    C <- col("C")
    D <- col("D")
    cost <- sum(tr$a * pc)
    variance <- sum(C / tr$a - D)
    if (mode == "cv") {
      if (!isTRUE(variance <= goal)) next
      value <- cost
    } else {
      if (!isTRUE(cost <= .budget_room(budget))) next
      value <- variance
    }
    if (value < best_value) {
      best_value <- value
      best <- list(a = tr$a, b = tr$b, cost = cost, pc = pc, C = C, D = D)
    }
  }
  best
}

#' The cv of a whole-cluster design, over its domains where it has them
#' @keywords internal
#' @noRd
.alloc_op_cv <- function(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate,
                         df) {
  if (length(prep$domain_idx) == 0L) {
    return(.alloc_metrics(
      N_h = prep$N_h,
      N_pop = .alloc_scope_N(prep), S_h = S_op, mean_h = prep$mean_h,
      n_h = e_h,
      alpha = alpha, deff = deff, resp_rate = resp_rate,
      cost_h = prep$cost_h, df = df, N_fpc = N_fpc_op
    )$cv)
  }
  prep_op <- prep
  prep_op$S_h <- S_op
  prep_op$N_fpc_h <- N_fpc_op
  .alloc_domain_cv_max(prep_op, e_h, alpha, deff, resp_rate)
}

#' @keywords internal
#' @noRd
.alloc_op_metric <- function(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate,
                             df, metric = "cv") {
  if (identical(metric, "cv")) {
    return(.alloc_op_cv(prep, e_h, S_op, N_fpc_op, alpha, deff, resp_rate, df))
  }
  if (length(prep$domain_idx) == 0L) {
    return(.alloc_metrics(
      N_h = prep$N_h,
      N_pop = .alloc_scope_N(prep), S_h = S_op, mean_h = prep$mean_h,
      n_h = e_h,
      alpha = alpha, deff = deff, resp_rate = resp_rate,
      cost_h = prep$cost_h, df = df, N_fpc = N_fpc_op
    )$se)
  }
  prep_op <- prep
  prep_op$S_h <- S_op
  prep_op$N_fpc_h <- N_fpc_op
  max(vapply(
    prep_op$domain_idx,
    function(idx) .alloc_domain_se(prep_op, e_h, alpha, deff, resp_rate, idx),
    numeric(1)
  ))
}

#' The se of one domain's estimate under an allocation
#' @keywords internal
#' @noRd
.alloc_domain_se <- function(prep, n_h, alpha, deff, resp_rate, idx) {
  .alloc_metrics(
    N_h = prep$N_h[idx],
    N_pop = .alloc_scope_N(prep)[idx],
    S_h = prep$S_h[idx],
    mean_h = prep$mean_h[idx],
    n_h = n_h[idx],
    alpha = alpha,
    deff = .subset_h(deff, idx),
    resp_rate = .subset_h(resp_rate, idx),
    cost_h = prep$cost_h[idx],
    N_fpc = prep$N_fpc_h[idx]
  )$se
}

#' @keywords internal
#' @noRd
.alloc_binding_label <- function(n_h, m_h, M_h) {
  on_lower <- abs(n_h - m_h) < 1e-6
  on_upper <- abs(n_h - M_h) < 1e-6
  if (any(on_lower)) {
    return("lower_bound")
  }
  if (any(on_upper)) {
    return("upper_bound")
  }
  "none"
}

#' @keywords internal
#' @noRd
.alloc_domain_cv_max <- function(prep, n_h, alpha, deff, resp_rate) {
  cvs <- vapply(
    prep$domain_idx,
    function(idx) .alloc_domain_cv(prep, n_h, alpha, deff, resp_rate, idx),
    numeric(1)
  )
  max(cvs)
}

#' The cv of one domain's estimate under an allocation
#' @keywords internal
#' @noRd
.alloc_domain_cv <- function(prep, n_h, alpha, deff, resp_rate, idx) {
  .alloc_metrics(
    N_h = prep$N_h[idx],
    N_pop = .alloc_scope_N(prep)[idx],
    S_h = prep$S_h[idx],
    mean_h = prep$mean_h[idx],
    n_h = n_h[idx],
    alpha = alpha,
    deff = .subset_h(deff, idx),
    resp_rate = .subset_h(resp_rate, idx),
    cost_h = prep$cost_h[idx],
    N_fpc = prep$N_fpc_h[idx]
  )$cv
}

#' The domain cv mode, each domain at the smallest total meeting the target
#'
#' The domains partition the strata, and a domain's cv reads only its own
#' strata, so the smallest design meeting every domain's target is each
#' domain at its own smallest total. Inside a domain the shares follow the
#' rule's weights `a_h`, which are scale free, so this is the domain
#' allocated as its own frame. A floor from `min_n_domain` raises a domain's
#' total to the floor where the target alone asks for less.
#' @return The stratum sizes, with the domains whose floor set their total
#'   as attribute `floored`.
#' @keywords internal
#' @noRd
.alloc_domain_minimum <- function(prep, a_h, m_h, M_h, cv, alpha, deff,
                                  resp_rate, floor_d = NULL, tol = 1e-8,
                                  target_se = NULL, target_mode = "cv") {
  n_h <- m_h
  floored <- rep(FALSE, length(prep$domain_idx))
  goal <- if (!is.null(target_se)) target_se else cv
  for (d in seq_along(prep$domain_idx)) {
    idx <- prep$domain_idx[[d]]
    val_at <- function(t) {
      n_try <- n_h
      n_try[idx] <- .rna_alloc(a_h[idx], t, m_h[idx], M_h[idx])
      at <- .cluster_take_rise(prep, n_try)
      if (!is.null(target_se)) {
        .alloc_domain_se(at, n_try, alpha, deff, resp_rate, idx)
      } else {
        .alloc_domain_cv(at, n_try, alpha, deff, resp_rate, idx)
      }
    }
    lo <- sum(m_h[idx])
    hi <- sum(M_h[idx])
    t <- if (val_at(lo) <= goal + tol) {
      lo
    } else {
      val_hi <- val_at(hi)
      if (!is.finite(val_hi) || val_hi > goal + tol) {
        stop(
          sprintf(
            "target '%s' is unattainable in domain '%s' under current constraints%s",
            target_mode,
            names(prep$domain_idx)[d] %||% d,
            .alloc_psu_bound_note(prep)
          ),
          call. = FALSE
        )
      }
      if (val_hi >= goal) {
        hi
      } else {
        uniroot(function(x) val_at(x) - goal, interval = c(lo, hi),
                tol = 1e-8)$root
      }
    }
    if (!is.null(floor_d) && floor_d[d] > t) {
      t <- floor_d[d]
      floored[d] <- TRUE
    }
    n_h[idx] <- .rna_alloc(a_h[idx], t, m_h[idx], M_h[idx])
  }
  structure(n_h, floored = floored)
}

#' @keywords internal
#' @noRd
.alloc_domain_summary <- function(
  prep,
  n_h,
  alpha,
  deff,
  resp_rate,
  df = NULL
) {
  if (length(prep$domain_idx) == 0L) {
    return(NULL)
  }

  out <- vector("list", length(prep$domain_idx))
  dom_names <- names(prep$domain_idx)
  for (i in seq_along(prep$domain_idx)) {
    idx <- prep$domain_idx[[i]]
    met <- .alloc_metrics(
      N_h = prep$N_h[idx],
      N_pop = .alloc_scope_N(prep)[idx],
      S_h = prep$S_h[idx],
      mean_h = prep$mean_h[idx],
      n_h = n_h[idx],
      alpha = alpha,
      deff = .subset_h(deff, idx),
      resp_rate = .subset_h(resp_rate, idx),
      cost_h = prep$cost_h[idx],
      df = df,
      N_fpc = prep$N_fpc_h[idx]
    )
    row <- if (!is.null(prep$domain_values)) {
      prep$domain_values[i, , drop = FALSE]
    } else {
      data.frame()
    }
    row$.domain <- dom_names[i]
    row$.n <- sum(n_h[idx])
    row$.se <- met$se
    row$.moe <- met$moe
    row$.rmoe <- met$rmoe
    row$.cv <- met$cv
    row$.cost <- met$cost
    out[[i]] <- row
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
