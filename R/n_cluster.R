#' Optimal multistage cluster allocation
#'
#' Compute optimal per-stage sample sizes for an S-stage cluster design
#' (S >= 2), minimizing cost for a given precision or minimizing variance
#' for a given budget.
#'
#' @param x First argument used for S3 dispatch. Supply a numeric vector of
#'   per-stage costs (same as `stage_cost`), a data frame of indicators (one
#'   row per indicator with stage homogeneity and precision target columns),
#'   a `svyplan_prec` result from [prec_cluster()], or a [svyplan()] profile.
#'   Leave `NULL` when supplying `stage_cost` by name.
#' @param stage_cost Numeric vector of per-stage costs (length S >= 2, ordered
#'   from outermost stage 1 to innermost stage S). Optional unique names on
#'   `stage_cost` (such as `c(district = 500, ea = 100, hh = 50)`) flow into
#'   output stage names. Unnamed vectors use `stage_1`, ..., `stage_S`.
#'   Optional when passed as `x` in the default method, required (directly or
#'   via `plan`) for the `data.frame` method.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param icc Numeric vector of stage homogeneities of length `stages - 1`
#'   (canonical conditional measures from [varcomp()]), or a `svyplan_varcomp`
#'   object. Typical values in household surveys range from 0.01 to 0.10.
#' @param unit_relvar Unit relvariance, the unit variance of the analysis
#'   variable over its squared mean, default 1. A proportion `p` has
#'   `(1 - p) / p`, so the default plans a proportion as if `p = 0.5`. A
#'   [varcomp()] result estimated from data sets it, and an explicit value is
#'   then refused. One backed out of a design effect sets none, so
#'   `unit_relvar` must be supplied beside it.
#' @param var_ratio Ratio of stage component unit variance to the analysis
#'   variable's unit variance, default 1. A scalar supplies stage 1 and
#'   derives deeper ratios via `k_s = k_{s-1} * (1 - icc_{s-1})`, exact for
#'   equal-sized or PPS sub-units. Supply a length `stages - 1` vector when
#'   unequal-sized sub-units are drawn by SRS. In that vector, an `NA` after
#'   the first entry is derived the same way. Refused when `icc` is a
#'   [varcomp()] result, which sets it.
#' @param cv Target coefficient of variation (relative standard error).
#'   Specify exactly one of `cv` or `budget` in the single-indicator method.
#' @param budget Total budget. Specify exactly one of `cv` or `budget` in the
#'   single-indicator method, or supply alongside an indicator table to
#'   minimize the worst precision ratio within budget.
#' @param stage_n Optional numeric vector of length `stages` fixing selected
#'   stage sizes. Use `NA` for free stages to optimize and positive numbers
#'   (`>= 1`) for fixed stages. At least one stage must remain `NA` (use
#'   [prec_cluster()] when all stages are fixed). Names, if given, place each
#'   value: the stage names, `stage_<k>`, or `n_1`, `take_2`, and so on.
#'   Other names are refused.
#' @param resp_rate Expected response rate(s) in `(0, 1]`, default 1. A scalar
#'   applies to the final stage (ultimate units) with earlier stages set to 1.
#'   A vector of length `stages` supplies stage-1 through stage-S rates. Its
#'   names, if given, place each rate: the stage names, `stage_<k>`, or
#'   `resp_rate_1`, ..., `resp_rate_{S-1}` with `resp_rate` for the last.
#' @param eligibility_rate Share of issued final-stage units in scope, such as
#'   occupied dwellings, in `(0, 1]`. Default 1.
#' @param eligible_per_unit Eligible units per final-stage unit, a positive
#'   number, default 1, such as married women per household. Can exceed 1.
#'   Multiplies the final-stage response rate in the realized take while sizes
#'   stay in final-stage units. With an indicator table, an
#'   `eligible_per_unit` column gives each row its own value, and the
#'   argument fills rows where it is absent or `NA`.
#' @param fixed_cost Fixed overhead cost (C0). Default 0. Total cost is
#'   `C0 + sum_{s=1}^S c_s * prod_{j=1}^s n_j`. With `domains`, the separate
#'   allocation sizes each domain as its own design and charges `C0` once per
#'   domain, while `allocation = "joint"` charges it once for the whole
#'   budget.
#' @param N Optional numeric vector of length `stages` giving aggregate
#'   population counts from stage 1 to stage S (outermost to innermost, so
#'   `N[1] <= N[2] <= ... <= N[S]`). For example, 40 districts, 400 EAs, and
#'   12000 households is `N = c(40, 400, 12000)`. Names, if given, place each
#'   count: the stage names, `stage_<k>`, or `N_1`, ..., `N_{S-1}` with `N`
#'   for the last.
#' @param fpc Which finite population correction the variance carries.
#'   `"auto"` (default) applies stage-by-stage correction (`"stage"`) whenever
#'   `N` (or population columns) is supplied, and `"none"` otherwise.
#'   `"stage"` requires `N` (or population columns). Supply population
#'   counts only with `"auto"` or `"stage"`.
#' @param domains Optional character vector naming domain columns in `x`
#'   (data frame method). Each domain is solved on its own unless
#'   `allocation = "joint"` in budget mode.
#' @param allocation How a budget is split across domains, either
#'   `"separate"` (default) or `"joint"` (one budget split across domains to
#'   minimize the worst precision ratio). Applies when `x` is a data frame
#'   with `domains` and `budget`. Under `"joint"`, every domain keeps at least
#'   one PSU in the continuous design as in the whole one, so the least a
#'   domain receives is the cost of its one-PSU design. The continuous split
#'   gives `$n` and `$cost`. The whole design in `$operational` is searched
#'   over the whole budget, not domain by domain within that split. It is
#'   the least worst ratio over whole designs that fit the budget, found by
#'   bisection between a feasible design and a ratio that no whole design
#'   reaches within the budget. The floors of each design (one PSU, or
#'   `min_n_domain`, in every domain) are checked against the budget on
#'   their own, and a refusal names the whole design when only its floors
#'   are out of reach.
#' @param domain_sampling How per-domain requirements combine. Only
#'   `"separate"` (default) is supported for multistage cluster designs.
#' @param min_n_domain Optional positive minimum total sample size per domain.
#' @param control `NULL` (default) or a list with `max_nodes`, the node limit
#'   of the whole-unit search (default `2e7`). Under `allocation = "joint"`,
#'   it bounds the joint search as a whole (the bisection and the spending of
#'   the leftover budget), and each domain's fit at its continuous share
#'   separately. `gap` is refused. A single search keeps no lower bound to
#'   stop at, and the joint search stops only when its bracket closes. Its
#'   measured gap is reported in `$operational$search$gap`.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_cluster` object with components:
#' \describe{
#'   \item{`n`}{Named numeric vector of continuous per-stage sample sizes
#'     (length `stages`). With several indicators the continuous solve
#'     returns the best design it found.}
#'   \item{`stages`}{Integer number of stages (S >= 2).}
#'   \item{`total_n`}{Continuous total sample size (`prod(n)`).}
#'   \item{`cv`}{Coefficient of variation of the continuous optimum.}
#'   \item{`cost`}{Cost of the continuous optimum.}
#'   \item{`operational`}{Whole-unit field design: `n` (named integer stage
#'     sizes), `total_n`, `cost`, and `cv`, all recomputed from the integer
#'     design, and `search`, with `status` (`"optimal"` when the search
#'     proved the design cheapest for a precision target, or most precise
#'     for a budget, `"limit"` when it stopped at its node limit with the
#'     best design found), `nodes` and `max_nodes`. With domains,
#'     `search$by_domain` gives each domain's status and
#'     `search$nodes_by_domain` its node count. Under
#'     `allocation = "joint"`, a domain's design is proven only jointly, so
#'     every domain carries the joint status, and `nodes_by_domain` splits the
#'     joint search's `nodes` among the domains. `search` also holds the worst
#'     ratio `ratio` of the returned design, the proven lower bound `bound` on
#'     it, their relative `gap`, the number of bisection `trials`, and each
#'     domain's whole cost in `shares`.}
#'   \item{`params`}{List of input parameters. Stage response rates are
#'     stored as `resp_rate_1`, ..., `resp_rate_{S-1}` and `resp_rate`, the
#'     names the allocation frames use.}
#' }
#'   With an indicator data frame, the same class also carries `$indicators`,
#'   `$detail`, and (when `domains` is set) `$domains` and
#'   `$operational$domains`.
#'
#' @details
#' ## General S-stage model
#'
#' For an \eqn{S}-stage design with stage variance components
#' \eqn{R_s = V k_s \delta_s} for \eqn{s = 1, \dots, S - 1} and
#' \eqn{R_S = V k_{S-1}(1 - \delta_{S-1})}{R_S = V k_(S-1) (1 - delta_(S-1))},
#' gross stage sizes \eqn{n_1, \dots, n_S}, stage response multipliers
#' \eqn{r_1, \dots, r_S} (with \eqn{r_S} including
#' `eligibility_rate * eligible_per_unit`), and per-parent capacities
#' \eqn{M_1 = N_1}, \eqn{M_s = N_s / N_{s-1}}{M_s = N_s / N_(s-1)} for
#' \eqn{1 < s < S}{1 < s < S} and
#' \eqn{M_S = N e / N_{S-1}}{M_S = N e / N_(S-1)}, where \eqn{e} is
#' `eligibility_rate * eligible_per_unit` (all infinite under
#' `fpc = "none"`), the squared CV is
#' \deqn{CV^2 = \sum_{s=1}^{S} \frac{R_s}{\prod_{j=1}^{s} n_j r_j}
#'   \left(1 - \frac{n_s r_s}{M_s}\right).}{CV^2 = sum_{s=1}^S (R_s / prod_{j=1}^s (n_j r_j)) (1 - n_s r_s / M_s).}
#' Regrouping terms by \eqn{1 / \prod_{j=1}^s n_j}{1 / prod_(j=1)^s n_j}
#' yields the telescoping form
#' \deqn{CV^2 = \sum_{s=1}^{S} \frac{\alpha_s}{\prod_{j=1}^{s} n_j} - \kappa,}{CV^2 = sum_{s=1}^S alpha_s / prod_{j=1}^s n_j - kappa,}
#' where \eqn{\kappa = R_1 / M_1},
#' \eqn{\alpha_s = (R_s - R_{s+1}/M_{s+1}) / \prod_{j=1}^s r_j}{alpha_s =
#' (R_s - R_(s+1) / M_(s+1)) / prod_(j=1)^s r_j} for \eqn{s < S}, and
#' \eqn{\alpha_S = R_S / \prod_{j=1}^S r_j}{alpha_S = R_S / prod_(j=1)^S r_j}.
#'
#' When all \eqn{\alpha_s > 0} and takes are unconstrained, Cauchy-Schwarz
#' gives the exact cost-optimal per-stage takes
#' \deqn{n_s^* = \sqrt{\frac{c_{s-1} \alpha_s}{c_s \alpha_{s-1}}}, \quad
#'   s = 2, \dots, S.}{n_s* = sqrt((c_{s-1} alpha_s) / (c_s alpha_{s-1})), s = 2, ..., S.}
#' Box constraints, fixed stages (`stage_n`), and non-positive \eqn{\alpha_s}
#' (when a small parent cluster makes full enumeration optimal) are solved
#' exactly via active-face elimination and log-convex coordinate optimization.
#' Whole-unit operational designs in `$operational` come from a
#' branch-and-bound search that proves the design cheapest for a precision
#' target, or most precise for a budget, or reports
#' `status = "limit"` when it stops at its node limit. A search that stops
#' before finding any design raises an error of class
#' `svyplan_plan_unknown`, since whether one exists is then unknown.
#'
#' ## Several indicators at once
#'
#' Pass an indicator data frame as the first argument to size one S-stage
#' design against multiple precision targets simultaneously. Required stage
#' homogeneity columns are `icc_1`, ..., `icc_{S-1}`. Optional stage columns
#' are `var_ratio_1`, ..., `var_ratio_{S-1}`, stage response rates
#' `resp_rate_1`, ..., `resp_rate_{S-1}` and `resp_rate` (final stage), and
#' population counts `N_1`, ..., `N_{S-1}` and `N` (final stage).
#'
#' The exact continuous solution above holds for one indicator. With
#' several, the continuous optimum comes from a local search in log takes
#' from several starting points, and it is not certified. The whole design
#' in `$operational` keeps its own proof, since its branch-and-bound search
#' does not depend on that starting point.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer. Chapter 9.
#'
#' @family cluster design functions
#' @seealso [prec_cluster()] for the inverse, [varcomp()] for estimating
#'   variance components, [design_effect()] for planning design effects,
#'   [n_multi()] for single-stage multi-indicator sizing, and [n_alloc()] for
#'   stratified multistage allocation.
#'
#' @examples
#' # 2-stage, budget mode
#' n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
#'
#' # 2-stage, CV mode with named stages
#' n_cluster(c(ea = 500, hh = 50), icc = 0.05, cv = 0.05)
#'
#' # 2-stage, fixed cluster take (MICS/DHS style: 20 households per EA)
#' n_cluster(c(500, 50), icc = 0.05, budget = 100000, stage_n = c(NA, 20))
#'
#' # 3-stage with finite population correction (fpc = "auto" enables stage FPC)
#' n_cluster(c(300, 50, 10), icc = c(0.04, 0.10), unit_relvar = 0.7 / 0.3,
#'           cv = 0.06, N = c(40, 400, 12000))
#'
#' # 4-stage design with custom stage names
#' n_cluster(
#'   c(district = 1000, ea = 200, segment = 60, hh = 15),
#'   icc = c(0.01, 0.03, 0.06),
#'   cv = 0.05
#' )
#'
#' # Several indicators, 2-stage CV mode
#' indicators <- data.frame(
#'   name  = c("stunting", "anemia"),
#'   p     = c(0.30, 0.10),
#'   cv    = c(0.10, 0.15),
#'   icc_1 = c(0.02, 0.05)
#' )
#' n_cluster(indicators, stage_cost = c(500, 50))
#'
#' # Design separately for two regions
#' by_region <- data.frame(region = c("North", "South"), name = "coverage",
#'                         p = c(0.5, 0.4), cv = 0.08, icc_1 = 0.05)
#' n_cluster(by_region, stage_cost = c(50, 10), domains = "region")$operational$domains
#'
#' @export
n_cluster <- function(x = NULL, ...) {
  if ("indicators" %in% ...names()) {
    stop(
      "pass the indicator table as the first argument: n_cluster(indicators, stage_cost = ...)",
      call. = FALSE
    )
  }
  if (!missing(x)) {
    .refuse_result_as_input(x, "x", "n_cluster")
    .res <- .dispatch_plan(x, "x", n_cluster, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_cluster")
}

#' @rdname n_cluster
#' @export
n_cluster.default <- function(
  x = NULL,
  ...,
  stage_cost = NULL,
  icc = NULL,
  unit_relvar = NULL,
  var_ratio = NULL,
  cv = NULL,
  budget = NULL,
  stage_n = NULL,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  fixed_cost = 0,
  N = NULL,
  fpc = c("auto", "none", "stage"),
  control = NULL,
  plan = NULL
) {
  dots <- list(...)
  if ("indicators" %in% names(dots)) {
    stop(
      "pass the indicator table as the first argument: n_cluster(indicators, stage_cost = ...)",
      call. = FALSE
    )
  }
  df_only <- intersect(
    names(dots),
    c("domains", "allocation", "domain_sampling", "min_n_domain")
  )
  if (length(df_only) > 0L) {
    stop(
      sprintf(
        "%s %s only when the first argument to n_cluster() is an indicator data frame",
        .quote_names(df_only),
        if (length(df_only) == 1L) "applies" else "apply"
      ),
      call. = FALSE
    )
  }
  if (!is.null(x)) {
    if (!is.null(stage_cost)) {
      stop(
        "supply stage costs either as the first argument or via 'stage_cost', not both",
        call. = FALSE
      )
    }
    stage_cost <- x
    x <- NULL
    if (inherits(plan, "svyplan") && !is.null(plan$defaults)) {
      plan$defaults$stage_cost <- NULL
    }
  }
  fpc <- match.arg(fpc)
  sc_len <- if (is.numeric(stage_cost)) length(stage_cost) else NULL
  plan <- .plan_for_cluster_mode(plan, FALSE, fpc, stages = sc_len)
  plan <- .check_varcomp_components(icc, unit_relvar, var_ratio, plan,
                                    "n_cluster()")
  .check_relvar_identified(
    icc,
    !is.null(unit_relvar) || "unit_relvar" %in% names(plan$defaults),
    "n_cluster()"
  )
  .plan <- .merge_plan_args(plan, n_cluster.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_cluster.default, c(.plan, dots)))
  .check_unused_dots(...)
  max_nodes <- .stage_max_nodes(control)

  for (nm in c("cv", "budget", "fixed_cost", "unit_relvar",
               "eligibility_rate", "eligible_per_unit")) {
    v <- get(nm)
    if (is.numeric(v)) assign(nm, unname(v))
  }

  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)
  unit_relvar <- unit_relvar %||% 1
  var_ratio <- var_ratio %||% 1
  if (is.null(stage_cost)) {
    stop("'stage_cost' is required (directly or via plan)", call. = FALSE)
  }
  if (is.null(icc)) {
    stop("'icc' is required (directly or via plan)", call. = FALSE)
  }
  check_stage_cost(stage_cost)
  stage_cost <- .reorder_stage_cost(stage_cost)
  stages <- length(stage_cost)
  stage_nms <- .stage_names(stages, names(stage_cost))
  cost_vec <- unname(as.numeric(stage_cost))
  names(stage_cost) <- stage_nms

  fpc_resolved <- .resolve_stage_fpc(fpc, N)
  yield <- eligibility_rate * eligible_per_unit
  pop <- if (identical(fpc_resolved, "stage")) {
    .stage_population(N, stages, yield = yield, stage_nms = stage_nms)
  } else {
    NULL
  }

  resp_vec <- .expand_stage_resp_rate(resp_rate, stages, stage_nms)

  if (inherits(icc, "svyplan_varcomp")) {
    vc <- icc
    if (!is.null(vc$strata)) {
      stop(
        "stratified varcomp: merge its $strata columns into an n_alloc() frame, or pass one stratum's values",
        call. = FALSE
      )
    }
    if (!is.null(pop)) .check_stage_fpc_varcomp(vc)
    icc <- vc$icc
    var_ratio <- vc$var_ratio
    if (!is.na(vc$unit_relvar)) unit_relvar <- vc$unit_relvar
  }

  check_scalar(unit_relvar, "unit_relvar")
  .check_var_ratio(var_ratio)

  icc <- .reorder_stage_vec(icc, "icc")
  var_ratio <- .reorder_stage_vec(var_ratio, "var_ratio")
  check_icc(icc, expected_length = stages - 1L)
  .check_cluster_icc_open(icc, context = "n_cluster()")
  var_ratio_input <- .var_ratio_given(var_ratio, length(icc))
  var_ratio <- .expand_var_ratio(var_ratio, icc)

  has_cv <- !is.null(cv)
  has_budget <- !is.null(budget)
  if (has_cv == has_budget) {
    stop("specify exactly one of 'cv' or 'budget'", call. = FALSE)
  }
  if (has_cv) check_scalar(cv, "cv")
  if (has_budget) check_scalar(budget, "budget")
  check_fixed_cost(fixed_cost, budget)

  sn <- .validate_stage_n(stage_n, stages, stage_nms)

  rates <- resp_vec
  rates[stages] <- rates[stages] * yield
  rate_label <- paste(
    c("resp_rate", if (eligibility_rate != 1) "eligibility_rate",
      if (eligible_per_unit != 1) "eligible_per_unit"),
    collapse = " * "
  )

  # Check fixed stage sizes against population caps and expected takes
  if (!is.null(sn)) {
    for (s in seq_len(stages)) {
      if (is.na(sn[s])) next
      if (!is.null(pop) && sn[s] > pop$M_take[s] + 1e-9) {
        stop(
          sprintf(
            "fixed %s (%g) exceeds population capacity %g at stage %d",
            stage_nms[s], sn[s], pop$M_take[s], s
          ),
          call. = FALSE
        )
      }
      if (s >= 2L) {
        .check_expected_take(
          sn[s], rates[s], stage_nms[s],
          if (s == stages) rate_label else sprintf("resp_rate[%d]", s)
        )
      }
    }
  }

  # Continuous bounds [lower, upper]
  lower <- c(1, vapply(rates[2L:stages], .take_bound, numeric(1L)))
  upper <- if (!is.null(pop)) pop$M_take else rep(Inf, stages)
  if (!is.null(sn)) {
    fixed <- !is.na(sn)
    lower[fixed] <- sn[fixed]
    upper[fixed] <- sn[fixed]
  }
  if (any(lower > upper + 1e-12)) {
    bad <- which(lower > upper + 1e-12)[1L]
    stop(
      sprintf(
        "stage %d (%s) population capacity %g cannot deliver one expected response at this response rate",
        bad, stage_nms[bad], upper[bad]
      ),
      call. = FALSE
    )
  }

  R_mat <- .stage_R_matrix(icc, unit_relvar, var_ratio)
  rates_mat <- matrix(rates, nrow = 1L)
  M_mat <- matrix(if (!is.null(pop)) pop$M else rep(Inf, stages), nrow = 1L)

  sol <- .stage_solve_continuous_one(
    R_mat, rates_mat, M_mat, cost_vec, lower, upper,
    cv = cv, budget = budget, fixed_cost = fixed_cost
  )
  if (!sol$ok) {
    if (identical(sol$reason, "budget")) {
      cheapest <- .stage_cost(lower, cost_vec, 0)
      stop(
        sprintf(
          "'budget' is too small for any realizable design: the cheapest design costs %.4g plus fixed_cost = %.4g",
          cheapest, fixed_cost
        ),
        call. = FALSE
      )
    }
    floor_cv <- sqrt(.stage_cv2(upper, R_mat, rates_mat, M_mat)[1L])
    stop(
      sprintf(
        "n_cluster(): target CV %.4g is below the achievable floor %.4g under the given stage constraints",
        cv, floor_cv
      ),
      call. = FALSE
    )
  }

  # Warn when unconstrained cost-optimal take fell below one expected response
  if (!is.null(sol$open_takes) && length(sol$open_takes) == stages - 1L) {
    for (idx in seq_len(stages - 1L)) {
      s <- idx + 1L
      ot <- sol$open_takes[idx]
      if (!is.na(ot) && is.finite(ot) && (is.null(sn) || is.na(sn[s])) &&
          ot * rates[s] < 1 - 1e-9) {
        warning(
          sprintf(
            "cost-optimal %s is below 1 expected response (high 'icc' relative to stage costs); clamped to %.4g",
            stage_nms[s], lower[s]
          ),
          call. = FALSE
        )
      }
    }
  }

  # Discrete whole-unit operational design via C branch-and-bound
  lower_int <- c(1, vapply(rates[2L:stages], .take_floor, integer(1L)))
  upper_int <- if (!is.null(pop)) floor(pop$M_take + 1e-9) else rep(Inf, stages)
  if (!is.null(sn)) {
    for (s in seq_len(stages)) {
      if (is.na(sn[s])) next
      w <- if (s == 1L) {
        min(upper_int[1L], .whole_psu(sn[1L]))
      } else {
        .whole_take(sn[s], lower_int[s], upper_int[s])
      }
      lower_int[s] <- w
      upper_int[s] <- w
    }
  }
  if (any(lower_int > upper_int)) {
    bad <- which(lower_int > upper_int)[1L]
    stop(
      sprintf(
        "stage %d (%s) population capacity %g cannot deliver one expected response in whole units",
        bad, stage_nms[bad], upper_int[bad]
      ),
      call. = FALSE
    )
  }

  vb <- if (!is.null(budget)) budget - fixed_cost else NA_real_
  cv2_tgt <- if (!is.null(cv)) cv^2 else 1.0
  disc <- .stage_solve_discrete(
    R_mat, rates_mat, M_mat, cv2_tgt, cost_vec,
    lower_int, upper_int, budget_var = vb, min_total = 0, n_init = sol$n,
    max_nodes = max_nodes
  )
  if (!isTRUE(disc$found)) {
    stop(
      if (!is.null(budget)) {
        "no whole-unit design fits the budget"
      } else {
        "target CV is not achievable with whole units under the given stage constraints"
      },
      call. = FALSE
    )
  }

  op_n <- as.integer(disc$n)
  names(op_n) <- stage_nms
  operational <- list(
    n = op_n,
    total_n = prod(as.numeric(op_n)),
    cost = .stage_cost(op_n, cost_vec, fixed_cost),
    cv = sqrt(.stage_cv2(op_n, R_mat, rates_mat, M_mat)[1L]),
    search = .stage_search_record(disc, max_nodes)
  )

  n_vec <- as.numeric(sol$n)
  names(n_vec) <- stage_nms
  names(icc) <- paste0("icc_", seq_len(stages - 1L))
  names(var_ratio) <- paste0("var_ratio_", seq_len(stages - 1L))

  params <- list(
    stage_cost = stage_cost,
    icc = icc,
    unit_relvar = unit_relvar,
    var_ratio = var_ratio,
    resp_rate = .name_stage_resp_rate(resp_vec)
  )
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(sn)) {
    names(sn) <- stage_nms
    params$stage_n <- sn
  }
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost
  if (eligibility_rate != 1) params$eligibility_rate <- eligibility_rate
  if (eligible_per_unit != 1) params$eligible_per_unit <- eligible_per_unit
  if (!is.null(pop)) {
    params$fpc <- "stage"
    params$N <- pop$N
  }
  if (anyNA(var_ratio_input)) params$var_ratio_input <- var_ratio_input
  if (!is.null(control)) params$control <- control

  .new_svyplan_cluster(
    n = n_vec,
    stages = stages,
    total_n = prod(n_vec),
    cv = sol$cv,
    cost = sol$cost,
    params = params,
    operational = operational
  )
}

#' @rdname n_cluster
#' @export
n_cluster.svyplan_prec <- function(x, ..., cv = NULL, budget = NULL) {
  if (identical(x[["type"]], "cluster") &&
      !is.null(x[["params"]][["indicators"]])) {
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
    if (x$cv > 0) {
      cv <- x$cv
    } else if (!is.null(p$budget)) {
      budget <- p$budget
    } else {
      .stop_census_target()
    }
  }
  args <- list(
    stage_cost = p$stage_cost,
    icc = p$icc,
    unit_relvar = p$unit_relvar,
    var_ratio = p[["var_ratio_input"]] %||% p$var_ratio,
    cv = cv,
    budget = budget,
    stage_n = p[["stage_n"]],
    resp_rate = p[["resp_rate"]] %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    fixed_cost = p[["fixed_cost"]] %||% 0,
    N = p[["N"]],
    fpc = p[["fpc"]] %||% "auto",
    control = p[["control"]]
  )
  do.call(n_cluster.default, .roundtrip_args(args, list(...), n_cluster.default))
}

#' Rebuild a several-indicators cluster allocation from its precision
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
  stages <- p[["stages"]]
  if ("resp_rate" %in% names(dots)) {
    for (s in seq_len(stages - 1L)) {
      col <- paste0("resp_rate_", s)
      if (col %in% names(tgt)) tgt[[col]] <- NA_real_
    }
    if ("resp_rate" %in% names(tgt)) tgt$resp_rate <- NA_real_
  }
  if ("eligible_per_unit" %in% names(dots) &&
      "eligible_per_unit" %in% names(tgt)) {
    tgt$eligible_per_unit <- NA_real_
  }
  tgt$n_1 <- NULL
  for (s in 2L:stages) {
    tgt[[paste0("take_", s)]] <- NULL
  }

  stored_mode <- p[["mode"]]
  census <- !is.null(x$detail) && any(x$detail$.cv <= 0)
  if (census) {
    budget <- budget %||% p[["budget"]]
    if (is.null(budget) || is.null(p[["cv_target"]])) .stop_census_target()
    tgt$cv <- p[["cv_target"]]
    tgt$moe <- NULL
  } else if (identical(stored_mode, "moe")) {
    tgt$moe <- x$detail$.moe
    tgt$cv <- NULL
  } else if (!is.null(x$detail) && ".cv" %in% names(x$detail)) {
    tgt$cv <- x$detail$.cv
    tgt$moe <- NULL
  }

  args <- list(
    x = tgt,
    stage_cost = p[["stage_cost"]],
    domains = p[["domain_cols"]],
    budget = budget %||% p[["budget"]],
    stage_n = p[["stage_n"]],
    allocation = p[["allocation"]] %||%
      if (isTRUE(p[["joint"]])) "joint" else "separate",
    min_n_domain = p[["min_n_domain"]],
    domain_sampling = p[["domain_sampling"]] %||% "separate",
    fixed_cost = p[["fixed_cost"]] %||% 0,
    resp_rate = p[["resp_rate"]] %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    eligible_per_unit = p[["eligible_per_unit"]] %||% 1,
    N = p[["N"]],
    fpc = p[["fpc"]] %||% "auto",
    control = p[["control"]]
  )
  do.call(
    n_cluster.data.frame,
    .roundtrip_args(args, dots, n_cluster.data.frame)
  )
}

#' Validate stage_n vector for n_cluster()
#' @keywords internal
#' @noRd
.validate_stage_n <- function(stage_n, stages, stage_nms) {
  if (is.null(stage_n)) return(NULL)
  if (!is.numeric(stage_n) || length(stage_n) != stages) {
    stop(
      sprintf(
        "'stage_n' must be a numeric vector of length %d (NA for free stages)",
        stages
      ),
      call. = FALSE
    )
  }
  sn <- as.numeric(.match_stage_vec(
    stage_n, "stage_n", stages, stage_nms,
    list(c("n_1", paste0("take_", 2L:stages)))
  ))
  if (all(is.na(sn))) return(NULL)
  if (all(!is.na(sn))) {
    stop("cannot fix all stages; use prec_cluster() instead", call. = FALSE)
  }
  fixed <- !is.na(sn)
  if (any(!is.finite(sn[fixed]) | sn[fixed] < 1)) {
    stop("fixed entries of 'stage_n' must be finite numbers >= 1",
         call. = FALSE)
  }
  sn
}

#' Check that every stage size stays within its population capacity
#' @keywords internal
#' @noRd
.check_stage_sizes <- function(n, pop, stage_names = NULL) {
  stages <- length(n)
  nms <- stage_names %||% names(n) %||% paste0("stage_", seq_len(stages))
  for (s in seq_len(stages)) {
    if (n[s] > pop$M_take[s] + 1e-9) {
      stop(
        sprintf(
          "'%s' (%g) exceeds population capacity %g at stage %d",
          nms[s], n[s], pop$M_take[s], s
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Refuse PPS varcomp components under stage-by-stage FPC
#' @keywords internal
#' @noRd
.check_stage_fpc_varcomp <- function(vc) {
  sel <- vc$selection
  if (!is.null(sel) && any(sel == "pps", na.rm = TRUE)) {
    bad <- names(sel)[which(sel == "pps")[1L]]
    stop(
      sprintf(
        "fpc = \"stage\" assumes equal-probability selection at every stage, and these components describe a size-proportional %s. Estimate them without 'prob', or plan with fpc = \"none\"",
        bad
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' A census reports a cv of zero, which no target can be rebuilt from
#' @keywords internal
#' @noRd
.stop_census_target <- function() {
  stop(
    "the design takes every unit at every stage, so its cv is 0 and states no target to size against. Pass 'cv' or 'budget'",
    call. = FALSE
  )
}

#' Resolve per-row S-stage population structures for an indicator table
#'
#' Reads aggregate population columns `N_1, ..., N_{S-1}, N` from `indicators`
#' (constant within each domain) or the vector `N` argument (when no domains).
#' Returns a list of length `nrow(indicators)` of `.stage_population()` objects,
#' or `NULL` when no stage FPC applies.
#' @keywords internal
#' @noRd
.stage_row_pops <- function(indicators, domain_cols, fpc, N_arg, stages,
                            yield_vec = 1) {
  pop_cols <- c(paste0("N_", seq_len(stages - 1L)), "N")
  present_cols <- intersect(pop_cols, names(indicators))
  # Also detect any out-of-range N_k column
  extra_n <- grep("^N_[0-9]+$", names(indicators), value = TRUE)
  bad_n <- setdiff(extra_n, pop_cols)
  if (length(bad_n) > 0L) {
    stop(
      sprintf(
        "unexpected population column(s) %s for a %d-stage design (expected %s)",
        .quote_names(bad_n), stages, .quote_names(pop_cols)
      ),
      call. = FALSE
    )
  }

  fpc_resolved <- .resolve_stage_fpc(
    fpc, N_arg, has_pop_cols = length(present_cols) > 0L
  )
  if (!identical(fpc_resolved, "stage")) {
    return(NULL)
  }

  nr <- nrow(indicators)
  if (length(yield_vec) == 1L) yield_vec <- rep(yield_vec, nr)

  if (length(present_cols) == 0L) {
    if (length(domain_cols) > 0L) {
      stop(
        sprintf(
          "with 'domains', fpc = \"stage\" reads each domain's population from 'indicators' columns %s, constant within a domain",
          .quote_names(pop_cols)
        ),
        call. = FALSE
      )
    }
    one <- .stage_population(N_arg, stages, yield = 1)
    return(lapply(yield_vec, function(y) {
      p <- one
      p$M[stages] <- p$M_take[stages] * y
      p
    }))
  }

  if (!is.null(N_arg)) {
    stop(
      "the population is given twice, as argument 'N' and as 'indicators' columns. Keep the columns or the argument",
      call. = FALSE
    )
  }
  missing_cols <- setdiff(pop_cols, present_cols)
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "fpc = \"stage\" for %d stages needs all population columns %s; missing %s",
        stages, .quote_names(pop_cols), .quote_names(missing_cols)
      ),
      call. = FALSE
    )
  }

  key <- if (length(domain_cols) > 0L) {
    .domain_key(indicators, domain_cols)
  } else {
    rep("all", nr)
  }
  out <- vector("list", nr)
  for (k in unique(key)) {
    rows <- which(key == k)
    N_dom <- numeric(stages)
    for (s in seq_len(stages)) {
      vals <- unique(indicators[[pop_cols[s]]][rows])
      if (length(vals) != 1L) {
        stop(
          sprintf(
            "'%s' must be constant within %s",
            pop_cols[s],
            if (length(domain_cols) > 0L) {
              sprintf("each domain, and varies in domain %s", sQuote(k))
            } else {
              "the table, which has no domains"
            }
          ),
          call. = FALSE
        )
      }
      N_dom[s] <- vals
    }
    one <- .stage_population(N_dom, stages, yield = 1)
    for (r in rows) {
      p <- one
      p$M[stages] <- p$M_take[stages] * yield_vec[r]
      out[[r]] <- p
    }
  }
  out
}

#' Collapse per-row populations of one domain into shared M_take and J x S M matrix
#' @keywords internal
#' @noRd
.collapse_row_pops <- function(pops) {
  base <- pops[[1L]]
  J <- length(pops)
  S <- length(base$N)
  M_mat <- matrix(0, nrow = J, ncol = S)
  for (j in seq_len(J)) {
    M_mat[j, ] <- pops[[j]]$M
  }
  list(N = base$N, M_take = base$M_take, M_mat = M_mat)
}
