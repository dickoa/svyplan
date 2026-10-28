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
#' @param unit_relvar Unit relvariance, the unit variance of the analysis
#'   variable over its squared mean, default 1. The cv a design reaches
#'   scales with it. A proportion `p` has `(1 - p) / p`, so the default plans
#'   a proportion as if `p = 0.5`. A [varcomp()] result estimated from data
#'   sets it, and an explicit value is then refused. One backed out of a
#'   design effect sets none, and `unit_relvar` must then be supplied beside
#'   it.
#' @param var_ratio Ratio of the stage components' unit variance to the analysis
#'   variable's, default 1. A scalar names `var_ratio_psu` and, for a three-stage
#'   design, derives `var_ratio_ssu = var_ratio_psu * (1 - icc_psu)`, the
#'   value for SSUs of equal size or drawn proportional to size. Supply a
#'   length-2 vector when SSUs of unequal size are drawn by simple random
#'   sampling, where it is larger. See [design_effect()]. The
#'   default of 1 is appropriate for most designs. A [varcomp()] result passed
#'   as `icc` sets it, and an explicit value is then refused.
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
#' @param eligible_rate Eligible units per final-stage unit, a positive
#'   number, default 1, such as married women per household when households
#'   are the final stage. It can exceed 1. It multiplies `resp_rate` in the
#'   realized take and in its cost, and sizes stay in final-stage units. Not
#'   available with `indicators`. See the nonresponse section of Details.
#' @param fixed_cost Fixed overhead cost (C0). Default 0.
#'   The total cost model becomes
#'   `C = C0 + c1*n_psu + c2*n_psu*n_per_psu [+ c3*n_psu*n_per_psu*n_per_ssu]`.
#'   In budget mode, only `budget - fixed_cost` is available for variable
#'   costs. In CV mode, `fixed_cost` is added to the variable cost.
#' @param fpc Which finite population correction the variance carries.
#'   `"none"` (default) drops every sampling fraction, the usual model when
#'   each stage samples a small share of its units. `"stage"` corrects each
#'   component by its own stage's sampling fraction and needs `N` and
#'   `N_psu`, and `N_ssu` at three stages. See "Finite populations".
#' @param N Number of final-stage units in the population, the units the
#'   takes count. With `eligible_rate` the correction reads
#'   `N * eligible_rate` eligible units. Read only with `fpc = "stage"` and
#'   refused otherwise.
#' @param N_psu Number of PSUs in the population. Read only with
#'   `fpc = "stage"` and refused otherwise.
#' @param N_ssu Number of SSUs in the population, across all PSUs. Read only
#'   with `fpc = "stage"` at three stages and refused otherwise.
#' @param indicators Optional data frame with one row per indicator, which
#'   switches the function to the several-indicators mode of Details. Each
#'   row carries its own homogeneity and its own precision target, and the
#'   allocation returned satisfies every row. Leave it `NULL` for the single
#'   indicator that `icc`, `unit_relvar`, `var_ratio` and `cv` describe.
#' @param domains Optional character vector naming domain columns in
#'   `indicators`. Each domain is solved on its own unless
#'   `allocation = "joint"` in budget mode, so a `budget` then applies to
#'   each domain and `$cost` is their sum. Domains are sized as separate
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
#' If `icc` is a `svyplan_varcomp` object, `icc` and `var_ratio` are taken
#' from it, and so is `unit_relvar` when the components were estimated from
#' data. A design effect identifies no unit relvariance, so a result of
#' `varcomp(deff = )` needs an explicit `unit_relvar`, such as `(1 - p) / p`
#' for a proportion. That pair is how a published design effect sizes a
#' design for a `cv` target.
#'
#' A value the result sets cannot also be passed, since the result would
#' override it. To plan with other values, pass `icc` as numbers, or give
#' `var_ratio` to `varcomp(deff = , var_ratio = )` for a back-out. A
#' [svyplan()] profile's `unit_relvar` and `var_ratio` give way to the
#' result, as any explicit argument overrides a profile.
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
#' `var_ratio_psu * (1 - icc_psu)` when absent, the value for SSUs of equal
#' size or drawn proportional to size (see [design_effect()]). On a ratio row `icc_psu`,
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
#' `eligible_rate` enters where `resp_rate` does at the final stage. The
#' realized take is the eligible respondents per PSU,
#' `n_per_psu * resp_rate * eligible_rate`, and `r` above is that product.
#' `icc` is then the homogeneity of eligible persons within a PSU, which a
#' design effect published at the person level identifies through
#' [varcomp()]. Above 1, the persons of one household resemble each other
#' more than persons of different households, and one `icc` averages over
#' that. For an outcome the whole household shares, such as poverty counted
#' per person, the clustering penalty is then understated.
#'
#' With households as the final stage, a survey that reports a household and
#' an individual response rate passes their product as `resp_rate`, as in
#' [n_prop()].
#'
#' This is a deterministic expected-take approximation. It substitutes the
#' expected realized stage counts into the planning variance and assumes
#' response is noninformative for it. Below one expected respondent per unit
#' the substitution stops describing the design, and that is an error naming
#' the approximation rather than the design. A take the optimum would put
#' below one expected respondent is raised to one, and the continuous and
#' whole-unit designs both keep it there. The smallest whole take is
#' `ceiling(1 / rate)`, so 2 units per PSU at `resp_rate = 0.5`, and a fixed
#' fractional take whose rounding would fall short rounds up to it. With
#' `indicators`, the rows share the takes, so the row with the lowest rate
#' sets the bound.
#'
#' ## Finite populations
#'
#' By default no finite population correction is applied at any stage,
#' which amounts to drawing every stage with replacement. That is the usual
#' planning model when each stage samples a small share of its units, as
#' with a few households from enumeration areas of several hundred.
#'
#' With `fpc = "stage"` each component is corrected by its own stage. For
#' PSUs drawn with equal probability and a common take,
#' \deqn{CV^2 = \frac{Vk}{a}\Big[\delta\Big(1 - \frac{a}{A}\Big)
#'   + \frac{1-\delta}{m}\Big(1 - \frac{m}{M}\Big)\Big],}{CV^2 = (V k / a) [delta (1 - a/A) + (1 - delta)/m (1 - m/M)],}
#' with \eqn{A} = `N_psu`, \eqn{M} = `N / N_psu` units per PSU, and
#' \eqn{a} and \eqn{m} the responding PSU count and take. A three-stage
#' design also needs `N_ssu`. Its SSU component is corrected by the share
#' taken of the `N_ssu / N_psu` SSUs per PSU, and its unit component by the
#' share taken of the `N / N_ssu` units per SSU. This is exact for PSUs and
#' SSUs of equal size drawn with equal probability, when `icc` and
#' `var_ratio` are the frame components [varcomp()] returns without `prob`.
#' Components from a size-proportional stage describe another design and
#' are refused.
#'
#' The correction moves the cost-optimal take to
#' \deqn{m = \sqrt{\frac{c_1(1-\delta)}{c_2\,r\,\{\delta - (1-\delta)/M\}}},}{m = sqrt(c_1 (1 - delta) / (c_2 r (delta - (1 - delta)/M))),}
#' and to the whole PSU when \eqn{\delta \le 1/(M+1)}{delta <= 1/(M + 1)}. A
#' further unit inside a PSU then always buys more precision than a further
#' PSU, so the efficient design takes whole clusters. The move matters for
#' small PSUs, such as plots or clusters of a few dozen units, and is
#' negligible at a few hundred units per PSU. At three stages the same
#' happens one level down. The unit fraction lowers the SSU component and
#' the SSU fraction lowers the PSU component, and when either reaches zero
#' the efficient design takes whole SSUs, or every SSU of a PSU.
#'
#' A design never takes more PSUs, SSUs or units than the population holds.
#' Once a budget or a target calls for more PSUs than `N_psu`, every PSU is
#' taken and the rest goes to the takes. A fixed take that would need more
#' PSUs than exist is refused, since a larger take may still reach the
#' target. `N` counts final-stage units, as the takes do, and with
#' `eligible_rate` the correction reads `N * eligible_rate` eligible units.
#'
#' With `indicators`, the arguments `N`, `N_psu` and `N_ssu` describe the
#' one population every row is measured on, and the design meets every
#' row's corrected target. With `domains`, each domain is a population of
#' its own, and its counts are `indicators` columns of the same names,
#' constant within the domain. That needs domains made of whole PSUs, such
#' as regions or urban and rural strata. A domain that cuts across PSUs,
#' such as an age group, has no PSU count of its own. The sizes in
#' `$domains` are continuous and the print rounds them up, so a printed take
#' can pass a mean PSU size that is not a whole number by less than one
#' unit.
#'
#' [n_alloc()] and [prec_alloc()] apply the same correction with
#' `fpc = "stage"`, and a one-row allocation returns the same continuous
#' design. Their default `fpc = "unit"` applies `1 - n / N` to the whole
#' variance instead, so without the argument the two functions agree only
#' where that fraction is negligible.
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
#' # Households per cluster for married women, 0.41 per household, for a
#' # proportion near 0.28, whose unit relvariance is (1 - p) / p
#' n_cluster(stage_cost = c(500, 50), icc = 0.07, unit_relvar = 0.72 / 0.28,
#'           cv = 0.05, resp_rate = 0.9, eligible_rate = 0.41)
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
#' # 100 plots of 49 trees: the stage correction moves the cost-optimal
#' # take from about 6.3 to 7
#' n_cluster(stage_cost = c(12, 3), icc = 5 / 54, unit_relvar = 0.54,
#'           budget = 1000)
#' n_cluster(stage_cost = c(12, 3), icc = 5 / 54, unit_relvar = 0.54,
#'           budget = 1000, fpc = "stage", N = 4900, N_psu = 100)
#'
#' # A capital of 30 EAs and a rural domain of 900, each its own population
#' rows <- data.frame(name = c("stunting", "anemia"), p = c(0.30, 0.10),
#'                    cv = c(0.08, 0.12), icc_psu = c(0.03, 0.12))
#' by_region <- rbind(
#'   transform(rows, region = "capital", N = 4500, N_psu = 30),
#'   transform(rows, region = "rural", N = 135000, N_psu = 900)
#' )
#' n_cluster(indicators = by_region, stage_cost = c(120, 10),
#'           domains = "region", fpc = "stage")
#'
#' # 40 villages of 10 hamlets of 30 households, three stages
#' n_cluster(stage_cost = c(300, 50, 10), icc = c(0.04, 0.10),
#'           unit_relvar = 0.7 / 0.3, cv = 0.06, fpc = "stage",
#'           N = 12000, N_psu = 40, N_ssu = 400)
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
  eligible_rate = 1,
  fixed_cost = 0,
  fpc = c("none", "stage"),
  N = NULL,
  N_psu = NULL,
  N_ssu = NULL,
  indicators = NULL,
  domains = NULL,
  allocation = c("separate", "joint"),
  domain_sampling = c("separate", "natural"),
  min_n_domain = NULL,
  plan = NULL
) {
  fpc <- match.arg(fpc)
  plan <- .plan_for_cluster_mode(plan, !is.null(indicators), fpc)
  plan <- .check_varcomp_components(icc, unit_relvar, var_ratio, plan,
                                    "n_cluster()")
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
  # A name on a scalar would ride into the stage names the result is read
  # by. Stage vectors keep theirs, which say which stage a value is for.
  for (nm in c("cv", "budget", "n_psu", "n_per_psu", "n_per_ssu",
               "fixed_cost", "unit_relvar", "resp_rate_psu", "resp_rate_ssu",
               "resp_rate", "eligible_rate", "N", "N_psu", "N_ssu")) {
    v <- get(nm)
    if (is.numeric(v)) assign(nm, unname(v))
  }
  # With indicators the population can come from columns instead, which
  # .stage_row_pops() reads.
  pop <- if (!is.null(indicators) &&
             all(vapply(list(N, N_psu, N_ssu), is.null, logical(1L)))) {
    .cluster_stage_population("none", NULL, NULL, "")
  } else {
    .cluster_stage_population(fpc, N, N_psu, "the default fpc = \"none\"",
                              N_ssu = N_ssu)
  }
  .check_cluster_mode_args(
    indicators,
    icc = icc, unit_relvar = unit_relvar, var_ratio = var_ratio, cv = cv,
    domains = domains, allocation = allocation,
    domain_sampling = domain_sampling, min_n_domain = min_n_domain
  )
  if (!is.null(indicators) && !identical(eligible_rate, 1)) {
    stop(
      "'eligible_rate' applies to the single-indicator call. With 'indicators', fold a rate at most 1 into each row's 'resp_rate'",
      call. = FALSE
    )
  }
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
      resp_rate = resp_rate,
      fpc = fpc,
      pop = pop
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
  check_eligible_rate(eligible_rate)
  pop <- .stage_eligible(pop, eligible_rate)

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
  .check_stage_population_depth(pop, stages)
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

  # The engines solve in eligible respondents at the ultimate stage, so the
  # take and its cost per respondent read the product of the two rates.
  ultimate_rate <- resp_rate * eligible_rate
  rate_label <- if (eligible_rate == 1) "resp_rate" else
    "resp_rate * eligible_rate"
  res <- if (stages == 2L && !is.null(pop)) {
    .n_cluster_2stage_stage(
      stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu, n_per_psu,
      resp_rate_psu, ultimate_rate, fixed_cost, pop, rate_label = rate_label
    )
  } else if (stages == 2L) {
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
      ultimate_rate,
      fixed_cost,
      rate_label = rate_label
    )
  } else if (!is.null(pop)) {
    .n_cluster_3stage_stage(
      stage_cost, icc, unit_relvar, .stage_k_pair(var_ratio, icc), cv, budget,
      n_psu, n_per_psu, n_per_ssu, resp_rate_psu, resp_rate_ssu,
      ultimate_rate, fixed_cost, pop, rate_label = rate_label
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
      ultimate_rate,
      fixed_cost,
      rate_label = rate_label
    )
  }
  res$params$resp_rate <- resp_rate
  if (eligible_rate != 1) {
    res$params$eligible_rate <- eligible_rate
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
    var_ratio = p$var_ratio,
    cv = cv,
    budget = budget,
    n_psu = p$n_psu,
    n_per_psu = p$n_per_psu,
    n_per_ssu = p$n_per_ssu,
    resp_rate_psu = p$resp_rate_psu %||% 1,
    resp_rate_ssu = p$resp_rate_ssu %||% 1,
    resp_rate = p$resp_rate %||% 1,
    eligible_rate = p$eligible_rate %||% 1,
    fixed_cost = p$fixed_cost %||% 0,
    fpc = p[["fpc"]] %||% "none",
    N = p[["N"]],
    N_psu = p[["N_psu"]],
    N_ssu = p[["N_ssu"]]
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
  # A census row has no variance left to state a target with. With the
  # budget at hand the targets the fit was sized against still rank rows.
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
    fixed_cost = p[["fixed_cost"]] %||% 0,
    fpc = p[["fpc"]] %||% "none",
    N = p[["N"]],
    N_psu = p[["N_psu"]],
    N_ssu = p[["N_ssu"]]
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

  m_floor <- .take_floor(resp_rate)
  m_cand <- if (!is.null(n_per_psu)) {
    max(m_floor, as.integer(round(n_per_psu)))
  } else {
    seq.int(
      m_floor,
      max(m_floor + 99L, min(10L * as.integer(ceiling(cont_m)) + 1L, 100000L))
    )
  }

  if (!is.null(budget)) {
    if (!is.null(n_psu)) {
      a <- max(1L, as.integer(round(n_psu)))
      m <- as.integer(floor((vb / a - C1) / C2))
      if (m < m_floor) {
        stop(
          sprintf(
            "no whole-unit design fits the budget for the given fixed stage-1 size with at least one expected response per PSU, which takes %d units",
            m_floor
          ),
          call. = FALSE
        )
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
      m_best <- max(m_floor, as.integer(ceiling(need / resp_rate - 1e-9)))
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
#' `m`. With the PSU count fixed and the take free (`m_free`), the search
#' ignores `m_cand` and covers every take that can matter, with the
#' ultimate take free or fixed at `q`. In cv mode a take reaches the target
#' only above `beta_j / (a - alpha_j)`, or `(beta_j + gamma_j / q) /
#' (a - alpha_j)` at a fixed `q`, for every `j`, and a design with take `m`
#' costs at least `a (C1 + C2 m)`, so the scan runs from the smallest
#' feasible take until that floor passes the best cost found. In budget mode
#' it runs to the largest take the budget pays with one unit per SSU, or
#' `q` at a fixed `q`. A target out of reach is then proved so.
#'
#' `a_max`, `m_max` and `q_max` bound the PSU count and the two takes by the
#' population under the stage correction. A budget-mode count or take above
#' its bound is lowered to it, and a cv-mode one is out of reach. The
#' correction also subtracts a constant `kappa_j` from each squared cv,
#' which moves the budget-mode ranking across several targets, so the key
#' is then `max_j (need_j / a - kappa_j)`.
#'
#' @return Named integer vector `c(n_psu, n_per_psu, n_per_ssu)`.
#' @keywords internal
#' @noRd
.op_cluster3_search <- function(alpha, beta, gamma, stage_cost,
                                variable_budget, n_psu, n_per_ssu, m_cand,
                                msg_budget, msg_cv, m_free = FALSE,
                                m_floor = 1L, q_floor = 1L,
                                a_max = Inf, m_max = Inf, q_max = Inf,
                                kappa = NULL) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  tol <- 1e-9
  m_cand <- as.numeric(m_cand)
  m_cand <- m_cand[m_cand >= m_floor]
  nm <- length(m_cand)
  nj <- length(alpha)
  budget_mode <- !is.null(variable_budget)
  a_fixed <- if (is.null(n_psu)) NULL else max(1, round(n_psu))
  q_fixed <- if (is.null(n_per_ssu)) NULL else max(q_floor, round(n_per_ssu))

  need_at <- function(m, q) {
    out <- alpha[1L] + beta[1L] / m + gamma[1L] / (m * q)
    for (j in seq_len(nj)[-1L]) {
      out <- pmax(out, alpha[j] + beta[j] / m + gamma[j] / (m * q))
    }
    out
  }
  ratio_at <- if (is.null(kappa)) {
    function(a, m, q) need_at(m, q) / a
  } else {
    function(a, m, q) {
      out <- (alpha[1L] + beta[1L] / m + gamma[1L] / (m * q)) / a - kappa[1L]
      for (j in seq_len(nj)[-1L]) {
        out <- pmax(out, (alpha[j] + beta[j] / m + gamma[j] / (m * q)) / a -
                      kappa[j])
      }
      out
    }
  }

  # Smallest whole ultimate take reaching a PSU count of 'a'; Inf where the
  # PSU count is out of reach at that SSU take whatever 'q' does.
  q_needed <- function(a, m) {
    out <- rep(q_floor, length(m))
    for (j in seq_len(nj)) {
      d <- a - alpha[j] - beta[j] / m
      out <- pmax(out, ifelse(d > 0, gamma[j] / (m * d), Inf))
    }
    ceiling(out - tol)
  }

  best <- list(key = Inf)
  consider <- function(a, m, q) {
    valid <- a >= 1 & is.finite(q) & q >= q_floor & q <= .Machine$integer.max &
      a <= a_max & q <= q_max
    q[!valid] <- q_floor
    cost <- a * (C1 + C2 * m + C3 * m * q)
    if (budget_mode) {
      valid <- valid & cost <= variable_budget + 1e-8
      key <- ratio_at(a, m, q) + 1e-12 * (m + q)
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
        pmin(q_max, floor((variable_budget / av - C1 - C2 * mv) / (C3 * mv) + tol))
      } else {
        q_needed(av, mv)
      }
      consider(av, mv, qv)
      i <- i + block
    }
  }

  if (!is.null(a_fixed) && m_free) {
    # The SSU take is the free stage, whether or not the ultimate take is
    # fixed with the PSU count.
    q_min <- q_fixed %||% q_floor
    take_at <- function(m) {
      if (!is.null(q_fixed)) {
        rep(q_fixed, length(m))
      } else if (budget_mode) {
        pmin(q_max, floor((variable_budget / a_fixed - C1 - C2 * m) / (C3 * m) + tol))
      } else {
        q_needed(a_fixed, m)
      }
    }
    scan_m <- function(m) {
      consider(rep(a_fixed, length(m)), m, take_at(m))
    }
    if (budget_mode) {
      m_hi <- min(m_max,
                  floor((variable_budget / a_fixed - C1) / (C2 + C3 * q_min) + tol))
      m_cur <- m_floor
      while (m_cur <= m_hi) {
        scan_m(seq.int(m_cur, min(m_hi, m_cur + 1e5 - 1)))
        m_cur <- m_cur + 1e5
      }
    } else if (all(a_fixed - alpha > tol)) {
      # The smallest take with every indicator in reach. A fixed ultimate
      # take adds its own term to the bound.
      m_cur <- if (is.null(q_fixed)) {
        max(m_floor, floor(max(beta / (a_fixed - alpha))) + 1)
      } else {
        max(m_floor, ceiling(max((beta + gamma / q_fixed) / (a_fixed - alpha)) - tol))
      }
      m_hi <- m_max
      while (m_cur <= m_hi) {
        m <- seq.int(m_cur, m_cur + 1e5 - 1)
        if (is.finite(m_hi)) m <- m[m <= m_hi]
        scan_m(m)
        m_cur <- m_cur + 1e5
        if (is.finite(best$key)) {
          m_hi <- min(m_max, floor((best$cost / a_fixed - C1) / C2 + 1e-9))
        }
      }
    }
  } else if (!is.null(q_fixed)) {
    q <- rep(q_fixed, nm)
    a <- if (!is.null(a_fixed)) {
      rep(a_fixed, nm)
    } else if (budget_mode) {
      pmin(a_max, floor((variable_budget + 1e-8) / (C1 + C2 * m_cand + C3 * m_cand * q)))
    } else {
      pmax(1, ceiling(need_at(m_cand, q) - tol))
    }
    consider(a, m_cand, q)
  } else if (!is.null(a_fixed)) {
    a <- rep(a_fixed, nm)
    q <- if (budget_mode) {
      pmin(q_max, floor((variable_budget / a_fixed - C1 - C2 * m_cand) / (C3 * m_cand) + tol))
    } else {
      q_needed(a_fixed, m_cand)
    }
    consider(a, m_cand, q)
  } else {
    # Cheapest a PSU can be, used to bound the PSU count from the cost side.
    m_lo <- min(m_cand)
    psu_floor <- C1 + m_lo * (C2 + C3 * q_floor)
    if (budget_mode) {
      a_hi <- min(a_max, floor((variable_budget + 1e-8) / psu_floor))
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
      a_seed <- min(a_max, max(a_lo, floor(max(alpha + beta / m_lo) + tol) + 1))
      consider(rep(a_seed, nm), m_cand, q_needed(a_seed, m_cand))
      block <- max(1L, as.integer(2e5 %/% nm))
      a_cur <- a_lo
      cost_bound <- function() {
        if (is.finite(best$key)) floor((best$cost + 1e-8) / psu_floor) else Inf
      }
      a_hi <- min(a_max, cost_bound())
      while (a_cur <= a_hi) {
        scan_a(seq.int(a_cur, min(a_hi, a_cur + block - 1L)))
        a_cur <- a_cur + block
        a_hi <- min(a_hi, cost_bound())
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

  m_floor <- .take_floor(resp_rate_ssu)
  m_cand <- if (!is.null(n_per_psu)) {
    max(m_floor, as.integer(round(n_per_psu)))
  } else {
    seq.int(
      m_floor,
      max(m_floor + 59L, min(6L * as.integer(ceiling(cont_m)) + 1L, 400L))
    )
  }

  hit <- .op_cluster3_search(
    alpha, beta, gamma, stage_cost, vb, n_psu, n_per_ssu, m_cand,
    msg_budget = "no whole-unit design fits the budget",
    msg_cv = "target CV is not achievable with whole units at the given fixed stage-1 size",
    m_free = is.null(n_per_psu),
    m_floor = m_floor,
    q_floor = .take_floor(resp_rate)
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
  fixed_cost = 0,
  rate_label = "resp_rate"
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
    # Less than one PSU at the cost-optimal take does not mean no design:
    # one PSU at a smaller take may fit, so the cheapest one is priced first.
    take_min <- n_per_psu_gross %||% 1
    cheapest <- stage_cost[1L] + stage_cost[2L] * take_min
    if (is.null(n_psu) && var_budget >= cheapest - 1e-8) {
      return(.n_cluster_2stage(
        stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu = 1,
        n_per_psu = n_per_psu_gross, resp_rate_psu = resp_rate_psu,
        resp_rate = resp_rate, fixed_cost = fixed_cost, rate_label = rate_label
      ))
    }
    stop(
      sprintf(
        "'budget' is too small for any realizable design: the cheapest design, one PSU with %s, costs %.4g plus fixed_cost = %.4g",
        if (is.null(n_per_psu_gross)) "one unit" else sprintf("n_per_psu = %.3g", take_min),
        cheapest, fixed_cost
      ),
      call. = FALSE
    )
  }

  # Back to gross units: the field team visits n2_opt / resp_rate units to
  # obtain n2_opt responses.
  n2_opt <- n2_opt / resp_rate
  .check_expected_take(n2_opt, resp_rate, "n_per_psu", rate_label)

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

#' Validate the populations a stage correction divides by
#'
#' `N` counts final-stage units, `N_psu` the PSUs and, at three stages,
#' `N_ssu` the SSUs. They are needed under `fpc = "stage"` and refused
#' otherwise, since a population that is read nowhere would suggest a
#' correction that was not applied. `M` is the mean PSU size in final-stage
#' units, and with `N_ssu`, `Q` and `C` are the SSUs per PSU and the units
#' per SSU. `.stage_eligible()` turns `M` and `C` into the eligible counts
#' the variance reads.
#' @keywords internal
#' @noRd
.cluster_stage_population <- function(fpc, N, N_psu, context, N_ssu = NULL) {
  if (!identical(fpc, "stage")) {
    given <- c("N", "N_psu", "N_ssu")[
      !vapply(list(N, N_psu, N_ssu), is.null, logical(1L))
    ]
    if (length(given) > 0L) {
      stop(
        sprintf(
          "%s %s the population%s that fpc = \"stage\" corrects for, and %s applies no correction. Add fpc = \"stage\", or drop %s",
          .quote_names(given),
          if (length(given) > 1L) "are" else "is",
          if (length(given) > 1L) "s" else "",
          context,
          .quote_names(given)
        ),
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (is.null(N) || is.null(N_psu)) {
    stop(
      "fpc = \"stage\" needs the population: 'N' ultimate units in 'N_psu' PSUs",
      call. = FALSE
    )
  }
  check_scalar(N_psu, "N_psu")
  check_scalar(N, "N")
  if (abs(N_psu - round(N_psu)) > 1e-8 || abs(N - round(N)) > 1e-8) {
    stop("'N' and 'N_psu' must be whole numbers", call. = FALSE)
  }
  if (N < N_psu) {
    stop("'N' must be at least 'N_psu', one ultimate unit per PSU", call. = FALSE)
  }
  pop <- list(N = round(N), A = round(N_psu), M = round(N) / round(N_psu))
  if (!is.null(N_ssu)) {
    check_scalar(N_ssu, "N_ssu")
    if (abs(N_ssu - round(N_ssu)) > 1e-8) {
      stop("'N_ssu' must be a whole number", call. = FALSE)
    }
    if (N_ssu < N_psu || N_ssu > N) {
      stop("'N_ssu' must lie between 'N_psu' and 'N'", call. = FALSE)
    }
    pop$N_ssu <- round(N_ssu)
    pop$Q <- pop$N_ssu / pop$A
    pop$C <- round(N) / pop$N_ssu
  }
  pop
}

#' The stage populations a design of this depth reads
#'
#' Three stages correct the SSU stage too, so they need `N_ssu`, and two
#' stages have no SSU stage for it to describe.
#' @keywords internal
#' @noRd
.check_stage_population_depth <- function(pop, stages) {
  if (is.null(pop)) return(invisible(TRUE))
  if (stages == 3L && is.null(pop$N_ssu)) {
    stop(
      "fpc = \"stage\" at three stages also needs 'N_ssu', the number of SSUs in the population",
      call. = FALSE
    )
  }
  if (stages == 2L && !is.null(pop$N_ssu)) {
    stop(
      "'N_ssu' applies to three-stage designs. A two-stage design reads 'N' and 'N_psu'",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The counts a stage population gives the takes and the variance
#'
#' A take counts final-stage units and is bounded by the `M` a PSU holds,
#' or the `C` an SSU holds. The variance reads eligible units, which number
#' `eligible_rate` per final-stage unit, so `M` and `C` become the eligible
#' counts and the bounds move to `M_take` and `C_take`.
#' @keywords internal
#' @noRd
.stage_eligible <- function(pop, eligible_rate = 1) {
  if (is.null(pop)) return(NULL)
  pop$M_take <- pop$M
  pop$M <- pop$M * eligible_rate
  if (!is.null(pop$C)) {
    pop$C_take <- pop$C
    pop$C <- pop$C * eligible_rate
  }
  pop
}

#' The population each indicator row is corrected for
#'
#' `N`, `N_psu` and `N_ssu` columns give each domain its own population,
#' constant within the domain. The scalar arguments give the whole table
#' one population and are refused with `domains`, where each domain needs
#' its own, and beside the columns. Columns without `fpc = "stage"` are
#' refused like the scalar arguments. Returns one population per row, with
#' the eligible counts of `.stage_eligible()`, or `NULL` without the
#' correction.
#' @keywords internal
#' @noRd
.stage_row_pops <- function(indicators, domain_cols, fpc, pop, stages) {
  cols <- intersect(c("N", "N_psu", "N_ssu"), names(indicators))
  if (!identical(fpc, "stage")) {
    if (length(cols) > 0L) {
      stop(
        sprintf(
          "the 'indicators' column%s %s %s the population that fpc = \"stage\" corrects for, and the default fpc = \"none\" applies no correction. Add fpc = \"stage\", or drop %s",
          if (length(cols) > 1L) "s" else "",
          .quote_names(cols),
          if (length(cols) > 1L) "are" else "is",
          if (length(cols) > 1L) "them" else "it"
        ),
        call. = FALSE
      )
    }
    return(NULL)
  }
  nr <- nrow(indicators)
  if (length(cols) == 0L) {
    if (length(domain_cols) > 0L) {
      stop(
        "with 'domains', fpc = \"stage\" reads each domain's population from 'indicators' columns 'N' and 'N_psu', and 'N_ssu' at three stages, constant within a domain",
        call. = FALSE
      )
    }
    if (is.null(pop)) {
      stop(
        "fpc = \"stage\" needs the population: 'N' ultimate units in 'N_psu' PSUs, as arguments or as 'indicators' columns",
        call. = FALSE
      )
    }
    .check_stage_population_depth(pop, stages)
    return(rep(list(pop), nr))
  }
  if (!is.null(pop)) {
    stop(
      "the population is given twice, as arguments and as 'indicators' columns. Keep the columns, which can differ by domain, or the arguments",
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
    vals <- lapply(cols, function(cl) unique(indicators[[cl]][rows]))
    names(vals) <- cols
    if (any(lengths(vals) != 1L)) {
      stop(
        sprintf(
          "%s must be constant within %s",
          .quote_names(cols),
          if (length(domain_cols) > 0L) {
            sprintf("each domain, and vary in domain %s", sQuote(k))
          } else {
            "the table, which has no domains"
          }
        ),
        call. = FALSE
      )
    }
    one <- .cluster_stage_population("stage", vals[["N"]], vals[["N_psu"]],
                                     "", N_ssu = vals[["N_ssu"]])
    .check_stage_population_depth(one, stages)
    out[rows] <- list(.stage_eligible(one))
  }
  out
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

#' Stage sizes no larger than the population they are drawn from
#'
#' `n` is the gross stage vector. The takes are bounded in final-stage
#' units, see `.stage_eligible()`.
#' @keywords internal
#' @noRd
.check_stage_sizes <- function(n, pop) {
  if (n[[1L]] > pop$A) {
    stop(sprintf("'n_psu' (%g) exceeds 'N_psu' (%g)", n[[1L]], pop$A),
         call. = FALSE)
  }
  if (length(n) == 2L && n[[2L]] > pop$M_take) {
    stop(
      sprintf("'n_per_psu' (%g) exceeds the units per PSU, %g", n[[2L]],
              pop$M_take),
      call. = FALSE
    )
  }
  if (length(n) == 3L && n[[2L]] > pop$Q) {
    stop(
      sprintf("'n_per_psu' (%g) exceeds the SSUs per PSU, N_ssu / N_psu = %g",
              n[[2L]], pop$Q),
      call. = FALSE
    )
  }
  if (length(n) == 3L && n[[3L]] > pop$C_take) {
    stop(
      sprintf("'n_per_ssu' (%g) exceeds the units per SSU, %g", n[[3L]],
              pop$C_take),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Components from a size-proportional stage do not describe the
#' equal-probability selection the stage correction assumes
#' @keywords internal
#' @noRd
.check_stage_fpc_varcomp <- function(vc) {
  if (identical(unname(vc$selection["psu"]), "pps")) {
    stop(
      "fpc = \"stage\" assumes PSUs drawn with equal probability, and these components describe a size-proportional first stage. Estimate them without 'prob', or plan without the stage correction",
      call. = FALSE
    )
  }
  if (identical(unname(vc$selection["ssu"]), "pps")) {
    stop(
      "fpc = \"stage\" assumes SSUs drawn with equal probability, and these components describe a size-proportional second stage. Estimate them without a second-stage 'prob', or plan without the stage correction",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Squared cv of a two-stage design with each stage's own correction
#'
#' `a_eff` and `b_eff` are the responding PSU count and take. Each component
#' is corrected by its own stage's responding fraction.
#' @keywords internal
#' @noRd
.cv2_cluster_stage <- function(a_eff, b_eff, icc, unit_relvar, var_ratio, A, M) {
  uk <- unit_relvar * var_ratio
  uk * icc * (1 - a_eff / A) / a_eff +
    uk * (1 - icc) * (1 - b_eff / M) / (a_eff * b_eff)
}

#' Two-stage design under the stage-by-stage correction
#'
#' Gross PSU count `a` and take `b` give
#' `cv^2 = [B - W/M + W/(b r)] / (a r1) - B/A`, with `B = U k icc` and
#' `W = U k (1 - icc)`, so the take that minimizes cost times variance is
#' `.cluster_free_take()` with the PSU size `M`. A design never takes more
#' than `N_psu` PSUs or more than `M` units per PSU. Past `N_psu` PSUs the
#' rest of the budget, or of the precision, goes to the take.
#' @keywords internal
#' @noRd
.n_cluster_2stage_stage <- function(
  stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu, n_per_psu,
  resp_rate_psu, resp_rate, fixed_cost, pop, rate_label = "resp_rate"
) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  A <- pop$A
  M <- pop$M
  M_take <- pop$M_take %||% M
  r1 <- resp_rate_psu
  r <- resp_rate
  B <- unit_relvar * var_ratio * icc
  W <- unit_relvar * var_ratio * (1 - icc)
  D <- B - W / M
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  target <- if (!is.null(cv)) cv^2 + B / A else NULL
  cv2 <- function(a, b) .cv2_cluster_stage(a * r1, b * r, icc, unit_relvar,
                                           var_ratio, A, M)
  a_for <- function(b) (D + W / (b * r)) / (r1 * target)
  b_for <- function(a) W / (r * (target * a * r1 - D))

  if (!is.null(n_psu) && n_psu > A) {
    stop(sprintf("'n_psu' (%g) exceeds 'N_psu' (%g)", n_psu, A), call. = FALSE)
  }
  if (!is.null(n_per_psu) && n_per_psu > M_take) {
    stop(
      sprintf("'n_per_psu' (%g) exceeds the units per PSU, %g", n_per_psu,
              M_take),
      call. = FALSE
    )
  }

  if (is.null(n_psu) && is.null(n_per_psu)) {
    b <- .cluster_free_take(C1, C2, icc, r, M, cap = M_take)
    if (b * r < 1) {
      warning(
        "cost-optimal n_per_psu is below 1 (high 'icc' relative to stage costs); clamped to 1",
        call. = FALSE
      )
      b <- min(M_take, 1 / r)
    }
    if (!is.null(budget)) {
      a <- vb / (C1 + C2 * b)
      if (a > A) {
        a <- A
        b <- min(M_take, (vb / A - C1) / C2)
      }
    } else {
      a <- a_for(b)
      if (a > A) {
        a <- A
        b <- b_for(A)
        if (!is.finite(b) || b <= 0 || b > M_take + 1e-9) {
          stop(
            sprintf(
              "target CV %.4g is out of reach: all %g PSUs at the full take of %g units reach %.4g",
              cv, A, M_take, sqrt(max(cv2(A, M_take), 0))
            ),
            call. = FALSE
          )
        }
      }
    }
  } else if (!is.null(n_psu)) {
    a <- n_psu
    if (!is.null(budget)) {
      b <- (vb - C1 * a) / (C2 * a)
      if (b * r < 1) {
        stop(
          "budget affords fewer than one unit per PSU for the given fixed stage-1 size",
          call. = FALSE
        )
      }
      b <- min(b, M_take)
    } else {
      floor_cv <- sqrt(max(cv2(a, M_take), 0))
      if (cv < floor_cv) {
        stop(
          sprintf(
            "n_cluster(): target CV %.4g is below the floor %.4g that n_psu = %g reaches at the full take of %g units",
            cv, floor_cv, a, M_take
          ),
          call. = FALSE
        )
      }
      b <- max(b_for(a), 1 / r)
    }
  } else {
    b <- n_per_psu
    if (!is.null(budget)) {
      a <- min(A, vb / (C1 + C2 * b))
    } else {
      a <- a_for(b)
      if (a > A) {
        stop(
          sprintf(
            "target CV %.4g needs %.4g PSUs at n_per_psu = %g, more than N_psu = %g. A larger take may reach it",
            cv, a, b, A
          ),
          call. = FALSE
        )
      }
    }
  }

  # One PSU at a smaller take may still fit a budget that the cost-optimal
  # take cannot spread over a whole PSU.
  smallest <- n_per_psu %||% .take_bound(r)
  if (!is.null(budget) && a < 1 && is.null(n_psu) && is.null(n_per_psu) &&
      vb >= C1 + C2 * smallest - 1e-8) {
    a <- 1
    b <- min(M_take, (vb - C1) / C2)
  }
  if (!is.null(budget) && a < 1) {
    stop(
      sprintf(
        "'budget' is too small for any realizable design: the cheapest design, one PSU with %s, costs %.4g plus fixed_cost = %.4g",
        if (is.null(n_per_psu)) sprintf("%.3g units", smallest) else
          sprintf("n_per_psu = %.3g", n_per_psu),
        C1 + C2 * smallest, fixed_cost
      ),
      call. = FALSE
    )
  }
  .check_expected_take(b, r, "n_per_psu", rate_label)

  operational <- .op_cluster_2stage_stage(
    stage_cost, cv2, a_for, b_for, cv, budget, n_psu, n_per_psu, r, A,
    M_take, fixed_cost, cont_m = b
  )

  params <- list(
    stage_cost = c(cost_psu = C1, cost_ssu = C2),
    icc = c(icc_psu = icc[1L]),
    unit_relvar = unit_relvar,
    var_ratio = c(var_ratio_psu = var_ratio[1L]),
    resp_rate_psu = resp_rate_psu,
    resp_rate = resp_rate,
    fpc = "stage",
    N = pop$N,
    N_psu = A
  )
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(n_per_psu)) params$n_per_psu <- n_per_psu
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost

  .new_svyplan_cluster(
    n = c(n_psu = a, n_per_psu = b),
    stages = 2L,
    total_n = a * b,
    cv = sqrt(max(cv2(a, b), 0)),
    cost = fixed_cost + a * (C1 + C2 * b),
    params = params,
    operational = operational
  )
}

#' Whole-unit two-stage design under the stage-by-stage correction
#'
#' The search of `.op_cluster_2stage()`, with the variance and the required
#' sizes of the corrected model and the population bounds `a <= N_psu` and
#' `m <= M`, where `M` is the PSU size in final-stage units.
#' @keywords internal
#' @noRd
.op_cluster_2stage_stage <- function(stage_cost, cv2, a_for, b_for, cv, budget,
                                     n_psu, n_per_psu, resp_rate, A, M,
                                     fixed_cost, cont_m) {
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  m_floor <- .take_floor(resp_rate)
  m_max <- floor(M + 1e-9)
  if (m_floor > m_max) {
    stop(
      sprintf(
        "a PSU of %g units cannot deliver one expected response at this response rate",
        M
      ),
      call. = FALSE
    )
  }
  m_cand <- if (!is.null(n_per_psu)) {
    min(m_max, max(m_floor, as.integer(round(n_per_psu))))
  } else {
    seq.int(
      m_floor,
      min(m_max, max(m_floor + 99L, min(10L * as.integer(ceiling(cont_m)) + 1L,
                                        100000L)))
    )
  }

  if (!is.null(budget)) {
    if (!is.null(n_psu)) {
      a_best <- min(A, max(1L, as.integer(round(n_psu))))
      m_best <- min(m_max, as.integer(floor((vb / a_best - C1) / C2)))
      if (m_best < m_floor) {
        stop(
          sprintf(
            "no whole-unit design fits the budget for the given fixed stage-1 size with at least one expected response per PSU, which takes %d units",
            m_floor
          ),
          call. = FALSE
        )
      }
    } else {
      a_cand <- pmin(A, floor(vb / (C1 + C2 * m_cand)))
      keep <- a_cand >= 1
      m_cand <- m_cand[keep]
      a_cand <- a_cand[keep]
      if (length(m_cand) == 0L) {
        stop("no whole-unit design fits the budget", call. = FALSE)
      }
      j <- which.min(cv2(a_cand, m_cand) + 1e-12 * m_cand)
      m_best <- m_cand[j]
      a_best <- a_cand[j]
    }
  } else {
    if (!is.null(n_psu)) {
      a_best <- min(A, max(1L, as.integer(round(n_psu))))
      need <- b_for(a_best)
      if (!is.finite(need) || need <= 0 || ceiling(need - 1e-9) > m_max) {
        stop(
          "target CV is not achievable with whole units at the given fixed stage-1 size",
          call. = FALSE
        )
      }
      m_best <- max(m_floor, as.integer(ceiling(need - 1e-9)))
    } else {
      a_cand <- pmax(1, ceiling(a_for(m_cand) - 1e-9))
      keep <- a_cand <= A
      if (!any(keep)) {
        stop(
          "target CV is not achievable with whole units within the PSU population",
          call. = FALSE
        )
      }
      m_cand <- m_cand[keep]
      a_cand <- a_cand[keep]
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
    cv = sqrt(max(cv2(a_best, m_best), 0))
  )
}

#' Three-stage coefficients under the stage-by-stage correction
#'
#' With gross sizes `a`, `b`, `c` and stage rates `r1`, `r2`, `r`, the
#' corrected variance is `cv^2 = G(b, c) / a - kappa`, where
#' `G = alpha + beta / b + gamma / (b c)`. Each fraction lowers the
#' component of the stage above it, `R1 - R2 / Q` and `R2 - R3 / C`, so
#' `alpha` and `beta` can be zero or negative. Vectorized over indicators.
#' @keywords internal
#' @noRd
.stage3_coef <- function(icc_psu, icc_ssu, unit_relvar, var_ratio_psu,
                         var_ratio_ssu, r1, r2, r, pop) {
  R1 <- unit_relvar * var_ratio_psu * icc_psu
  R2 <- unit_relvar * var_ratio_ssu * icc_ssu
  R3 <- unit_relvar * var_ratio_ssu * (1 - icc_ssu)
  list(
    alpha = (R1 - R2 / pop$Q) / r1,
    beta = (R2 - R3 / pop$C) / (r1 * r2),
    gamma = R3 / (r1 * r2 * r),
    kappa = R1 / pop$A
  )
}

#' Squared cv of a three-stage design with each stage's own correction
#'
#' `a_eff`, `b_eff` and `c_eff` are the responding sizes, and each
#' component is corrected by its own stage's responding fraction.
#' @keywords internal
#' @noRd
.cv2_cluster3_stage <- function(a_eff, b_eff, c_eff, icc, unit_relvar,
                                var_ratio, pop) {
  R1 <- unit_relvar * var_ratio[1L] * icc[1L]
  R2 <- unit_relvar * var_ratio[2L] * icc[2L]
  R3 <- unit_relvar * var_ratio[2L] * (1 - icc[2L])
  R1 * (1 - a_eff / pop$A) / a_eff +
    R2 * (1 - b_eff / pop$Q) / (a_eff * b_eff) +
    R3 * (1 - c_eff / pop$C) / (a_eff * b_eff * c_eff)
}

#' Cost-optimal takes of a three-stage design within population bounds
#'
#' Minimizes `G(b, c) H(b, c)`, with `H = c1 + c2 b + c3 b c` the cost per
#' PSU, over the box `box$b` by `box$c`. A fixed take has a box of one
#' point. The product has one stationary point, the no-FPC optimum on the
#' corrected components, when `alpha` and `beta` are positive. Otherwise,
#' or when that point lies outside the box, the minimum is on an edge, and
#' along each edge the product has the form `(K + u / t)(L + v t)`, whose
#' minimum is closed. A nonpositive `beta` sends `c` to its upper bound and
#' a nonpositive `alpha` sends `b` to its upper bound. `b_open` and
#' `c_open` are the edge minima before the box was applied.
#' @keywords internal
#' @noRd
.stage3_best_takes <- function(co, stage_cost, box) {
  al <- co$alpha
  be <- co$beta
  ga <- co$gamma
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  clamp <- function(x, lim) min(lim[2L], max(lim[1L], x))
  c_open <- function(b) {
    K <- al + be / b
    if (K > 0) sqrt(ga * (C1 + C2 * b) / (K * C3)) / b else Inf
  }
  b_open <- function(cc) {
    if (al > 0) sqrt((be + ga / cc) * C1 / (al * (C2 + C3 * cc))) else Inf
  }
  cand <- rbind(
    c(box$b[1L], clamp(c_open(box$b[1L]), box$c)),
    c(box$b[2L], clamp(c_open(box$b[2L]), box$c)),
    c(clamp(b_open(box$c[1L]), box$b), box$c[1L]),
    c(clamp(b_open(box$c[2L]), box$b), box$c[2L])
  )
  if (al > 0 && be > 0) {
    b0 <- sqrt(be * C1 / (al * C2))
    c0 <- sqrt(ga * C2 / (be * C3))
    if (b0 >= box$b[1L] && b0 <= box$b[2L] &&
        c0 >= box$c[1L] && c0 <= box$c[2L]) {
      cand <- rbind(c(b0, c0), cand)
    }
  }
  prod_gh <- (al + be / cand[, 1L] + ga / (cand[, 1L] * cand[, 2L])) *
    (C1 + C2 * cand[, 1L] + C3 * cand[, 1L] * cand[, 2L])
  j <- which.min(prod_gh)
  b <- cand[j, 1L]
  cc <- cand[j, 2L]
  list(b = b, c = cc, b_open = b_open(cc), c_open = c_open(b))
}

#' Three-stage takes at a set PSU count
#'
#' With `a` PSUs, a cv target `cv2_target = cv^2 + kappa` asks for the
#' cheapest takes with `G <= a cv2_target`, and a budget asks for the
#' smallest `G` with `H <= per_psu`. Either way the constraint binds unless
#' the box corner settles it, and along it both problems minimize
#' `(beta + gamma / c)(c2 + c3 c)`, so `c` is the free-take ratio
#' `sqrt(gamma c2 / (beta c3))` within the range that keeps `b` in its box.
#' `NULL` when no takes in the box meet the constraint.
#' @keywords internal
#' @noRd
.stage3_takes_at <- function(co, stage_cost, box, a, cv2_target = NULL,
                             per_psu = NULL) {
  al <- co$alpha
  be <- co$beta
  ga <- co$gamma
  C1 <- stage_cost[1L]
  C2 <- stage_cost[2L]
  C3 <- stage_cost[3L]
  G <- function(b, cc) al + be / b + ga / (b * cc)
  H <- function(b, cc) C1 + C2 * b + C3 * b * cc
  lo <- c(box$b[1L], box$c[1L])
  hi <- c(box$b[2L], box$c[2L])
  if (!is.null(cv2_target)) {
    g0 <- a * cv2_target
    if (G(hi[1L], hi[2L]) > g0 * (1 + 1e-12)) return(NULL)
    if (G(lo[1L], lo[2L]) <= g0) return(list(b = lo[1L], c = lo[2L]))
    room <- g0 - al
    b_of <- function(cc) (be + ga / cc) / room
    c_range <- c(
      if (hi[1L] * room - be > 0) ga / (hi[1L] * room - be) else Inf,
      if (lo[1L] * room - be > 0) ga / (lo[1L] * room - be) else Inf
    )
  } else {
    if (H(lo[1L], lo[2L]) > per_psu * (1 + 1e-12)) return(NULL)
    if (H(hi[1L], hi[2L]) <= per_psu) return(list(b = hi[1L], c = hi[2L]))
    spend <- per_psu - C1
    b_of <- function(cc) spend / (C2 + C3 * cc)
    c_range <- c((spend / hi[1L] - C2) / C3, (spend / lo[1L] - C2) / C3)
  }
  c_range <- c(max(lo[2L], c_range[1L]), min(hi[2L], c_range[2L]))
  c_free <- if (be > 0) sqrt(ga * C2 / (be * C3)) else Inf
  cc <- min(c_range[2L], max(c_range[1L], c_free))
  list(b = min(hi[1L], max(lo[1L], b_of(cc))), c = cc)
}

#' Three-stage design under the stage-by-stage correction
#'
#' The cost-optimal takes come from `.stage3_best_takes()`. When they call
#' for more PSUs than `N_psu`, or `n_psu` is fixed, the PSU count is set and
#' `.stage3_takes_at()` spends the rest of the budget, or of the precision,
#' on the takes. A design never takes more than `N_psu` PSUs, `Q = N_ssu /
#' N_psu` SSUs per PSU or `C = N / N_ssu` units per SSU.
#' @keywords internal
#' @noRd
.n_cluster_3stage_stage <- function(
  stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu, n_per_psu,
  n_per_ssu, resp_rate_psu, resp_rate_ssu, resp_rate, fixed_cost, pop,
  rate_label = "resp_rate"
) {
  A <- pop$A
  Q <- pop$Q
  C <- pop$C_take %||% pop$C
  r1 <- resp_rate_psu
  r2 <- resp_rate_ssu
  r <- resp_rate
  co <- .stage3_coef(icc[1L], icc[2L], unit_relvar, var_ratio[1L],
                     var_ratio[2L], r1, r2, r, pop)
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  cv2_target <- if (!is.null(cv)) cv^2 + co$kappa else NULL
  G <- function(b, cc) co$alpha + co$beta / b + co$gamma / (b * cc)
  H <- function(b, cc) sum(stage_cost * c(1, b, b * cc))
  cv_at <- function(a, b, cc) sqrt(max(G(b, cc) / a - co$kappa, 0))

  if (!is.null(n_psu) && n_psu > A) {
    stop(sprintf("'n_psu' (%g) exceeds 'N_psu' (%g)", n_psu, A), call. = FALSE)
  }
  if (!is.null(n_per_psu) && n_per_psu > Q) {
    stop(
      sprintf("'n_per_psu' (%g) exceeds the SSUs per PSU, N_ssu / N_psu = %g",
              n_per_psu, Q),
      call. = FALSE
    )
  }
  if (!is.null(n_per_ssu) && n_per_ssu > C) {
    stop(
      sprintf("'n_per_ssu' (%g) exceeds the units per SSU, %g", n_per_ssu, C),
      call. = FALSE
    )
  }
  box <- list(
    b = if (!is.null(n_per_psu)) rep(n_per_psu, 2L) else c(.take_bound(r2), Q),
    c = if (!is.null(n_per_ssu)) rep(n_per_ssu, 2L) else c(.take_bound(r), C)
  )
  if (box$b[1L] > box$b[2L] || box$c[1L] > box$c[2L]) {
    stop(
      sprintf(
        "a PSU of %g SSUs of %g units cannot deliver one expected response per stage at these response rates",
        Q, C
      ),
      call. = FALSE
    )
  }

  if (is.null(n_psu)) {
    best <- .stage3_best_takes(co, stage_cost, box)
    b <- best$b
    cc <- best$c
    if (is.null(n_per_ssu) && best$c_open * r < 1 - 1e-9) {
      warning(
        "cost-optimal n_per_ssu is below 1 (high 'icc_ssu' relative to stage costs); clamped to 1",
        call. = FALSE
      )
    }
    if (is.null(n_per_psu) && best$b_open * r2 < 1 - 1e-9) {
      warning(
        "cost-optimal n_per_psu is below 1 (high 'icc_psu' relative to stage costs); clamped to 1",
        call. = FALSE
      )
    }
    a <- if (!is.null(budget)) vb / H(b, cc) else G(b, cc) / cv2_target
    if (!is.null(budget) && a < 1) {
      at_one <- .stage3_takes_at(co, stage_cost, box, 1, per_psu = vb)
      if (is.null(at_one)) {
        stop(
          sprintf(
            "'budget' is too small for any realizable design: the cheapest design, one PSU at the smallest takes, costs %.4g plus fixed_cost = %.4g",
            H(box$b[1L], box$c[1L]), fixed_cost
          ),
          call. = FALSE
        )
      }
      a <- 1
      b <- at_one$b
      cc <- at_one$c
    }
    if (a > A) {
      a <- A
      at_all <- .stage3_takes_at(co, stage_cost, box, A,
                                 cv2_target = cv2_target,
                                 per_psu = if (!is.null(budget)) vb / A)
      if (is.null(at_all)) {
        stop(
          sprintf(
            "target CV %.4g is out of reach: all %g PSUs at the largest takes allowed reach %.4g%s",
            cv, A, cv_at(A, box$b[2L], box$c[2L]),
            if (!is.null(n_per_psu) || !is.null(n_per_ssu)) {
              ". A larger fixed take may reach it"
            } else {
              ""
            }
          ),
          call. = FALSE
        )
      }
      b <- at_all$b
      cc <- at_all$c
    }
  } else {
    a <- n_psu
    at_n <- .stage3_takes_at(co, stage_cost, box, a,
                             cv2_target = cv2_target,
                             per_psu = if (!is.null(budget)) vb / a)
    if (is.null(at_n)) {
      if (!is.null(budget)) {
        stop(
          "budget affords fewer than one unit per SSU for the given fixed stage sizes",
          call. = FALSE
        )
      }
      stop(
        sprintf(
          "n_cluster(): target CV %.4g is below the floor %.4g that n_psu = %g reaches at the largest takes allowed",
          cv, cv_at(a, box$b[2L], box$c[2L]), a
        ),
        call. = FALSE
      )
    }
    b <- at_n$b
    cc <- at_n$c
  }
  .check_expected_take(b, r2, "n_per_psu", "resp_rate_ssu")
  .check_expected_take(cc, r, "n_per_ssu", rate_label)

  operational <- .op_cluster_3stage_stage(
    stage_cost, co, cv, budget, n_psu, n_per_psu, n_per_ssu, r2, r, pop,
    fixed_cost, cont_m = b
  )

  params <- list(
    stage_cost = c(cost_psu = stage_cost[1L], cost_ssu = stage_cost[2L],
                   cost_tsu = stage_cost[3L]),
    icc = c(icc_psu = icc[1L], icc_ssu = icc[2L]),
    unit_relvar = unit_relvar,
    var_ratio = c(var_ratio_psu = var_ratio[1L], var_ratio_ssu = var_ratio[2L]),
    resp_rate_psu = resp_rate_psu,
    resp_rate_ssu = resp_rate_ssu,
    resp_rate = resp_rate,
    fpc = "stage",
    N = pop$N,
    N_psu = A,
    N_ssu = pop$N_ssu
  )
  if (!is.null(cv)) params$cv <- cv
  if (!is.null(budget)) params$budget <- budget
  if (!is.null(n_psu)) params$n_psu <- n_psu
  if (!is.null(n_per_psu)) params$n_per_psu <- n_per_psu
  if (!is.null(n_per_ssu)) params$n_per_ssu <- n_per_ssu
  if (fixed_cost > 0) params$fixed_cost <- fixed_cost

  .new_svyplan_cluster(
    n = c(n_psu = a, n_per_psu = b, n_per_ssu = cc),
    stages = 3L,
    total_n = a * b * cc,
    cv = cv_at(a, b, cc),
    cost = fixed_cost + a * H(b, cc),
    params = params,
    operational = operational
  )
}

#' Whole-unit three-stage design under the stage-by-stage correction
#'
#' `.op_cluster3_search()` on the corrected coefficients, bounded by the
#' population at every stage.
#' @keywords internal
#' @noRd
.op_cluster_3stage_stage <- function(stage_cost, co, cv, budget, n_psu,
                                     n_per_psu, n_per_ssu, resp_rate_ssu,
                                     resp_rate, pop, fixed_cost, cont_m) {
  vb <- if (!is.null(budget)) budget - fixed_cost else NULL
  scale <- if (is.null(cv)) 1 else cv^2 + co$kappa
  m_floor <- .take_floor(resp_rate_ssu)
  q_floor <- .take_floor(resp_rate)
  m_max <- floor(pop$Q + 1e-9)
  q_max <- floor((pop$C_take %||% pop$C) + 1e-9)
  if (m_floor > m_max || q_floor > q_max) {
    stop(
      sprintf(
        "a PSU of %g SSUs of %g units cannot deliver one expected response per stage at these response rates",
        pop$Q, pop$C_take %||% pop$C
      ),
      call. = FALSE
    )
  }
  m_cand <- if (!is.null(n_per_psu)) {
    min(m_max, max(m_floor, as.integer(round(n_per_psu))))
  } else {
    .multi_stage_candidates(NULL, cont_m, upper = m_max, base_limit = 400L,
                            floor = m_floor)
  }
  hit <- .op_cluster3_search(
    co$alpha / scale, co$beta / scale, co$gamma / scale, stage_cost, vb,
    n_psu, if (!is.null(n_per_ssu)) min(q_max, n_per_ssu), m_cand,
    msg_budget = "no whole-unit design fits the budget",
    msg_cv = "target CV is not achievable with whole units within the population",
    m_free = is.null(n_per_psu),
    m_floor = m_floor,
    q_floor = q_floor,
    a_max = pop$A,
    m_max = m_max,
    q_max = q_max
  )
  a <- hit[["n_psu"]]
  m <- hit[["n_per_psu"]]
  q <- hit[["n_per_ssu"]]
  G <- co$alpha + co$beta / m + co$gamma / (m * q)
  list(
    n = c(n_psu = a, n_per_psu = m, n_per_ssu = q),
    total_n = a * m * q,
    cost = fixed_cost + a * sum(stage_cost * c(1, m, m * q)),
    cv = sqrt(max(G / a - co$kappa, 0))
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
  fixed_cost = 0,
  rate_label = "resp_rate"
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

  n3_clamped <- FALSE
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
      n3_clamped <- TRUE
    }
    n3_free
  }

  n2 <- if (!is.null(n_per_psu)) {
    n_per_psu
  } else if (solve_for != "n2") {
    # A clamped n3 is a fixed one for n2, so n2 takes the restricted optimum.
    n2_free <- if (!is.null(n_per_ssu) || n3_clamped) {
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
      if (is.null(n_per_ssu)) {
        # The free n3 assumed the unclamped n2, so it is re-solved at n2 = 1
        # as the restricted optimum above solves it for a fixed n2.
        n3 <- max(1, sqrt(k2 * (1 - delta2) * (C1 + C2) /
                            ((k1 * delta1 + k2 * delta2) * C3)))
      }
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
      # One SSU at a smaller free take may still fit.
      if (n2 < 1 && is.null(n_per_ssu) &&
          (var_budget / n_psu - C1 - C2) / C3 >= 1 - 1e-9) {
        n2 <- 1
        n3 <- (var_budget / n_psu - C1 - C2) / C3
      }
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
    # As in two stages: one PSU at smaller takes may fit.
    m_min <- n_per_psu_gross %||% 1
    q_min <- n_per_ssu_gross %||% 1
    cheapest <- stage_cost[1L] + stage_cost[2L] * m_min + stage_cost[3L] * m_min * q_min
    if (is.null(n_psu) && var_budget >= cheapest - 1e-8) {
      return(.n_cluster_3stage(
        stage_cost, icc, unit_relvar, var_ratio, cv, budget, n_psu = 1,
        n_per_psu = n_per_psu_gross, n_per_ssu = n_per_ssu_gross,
        resp_rate_psu = resp_rate_psu, resp_rate_ssu = resp_rate_ssu,
        resp_rate = resp_rate, fixed_cost = fixed_cost, rate_label = rate_label
      ))
    }
    stop(
      sprintf(
        "'budget' is too small for any realizable design: the cheapest design, one PSU with %s, costs %.4g plus fixed_cost = %.4g",
        if (is.null(n_per_psu_gross) && is.null(n_per_ssu_gross)) {
          "one SSU and one unit"
        } else {
          sprintf("n_per_psu = %.3g and n_per_ssu = %.3g", m_min, q_min)
        },
        cheapest, fixed_cost
      ),
      call. = FALSE
    )
  }

  n2 <- n2 / resp_rate_ssu
  n3 <- n3 / resp_rate
  .check_expected_take(n2, resp_rate_ssu, "n_per_psu", "resp_rate_ssu")
  .check_expected_take(n3, resp_rate, "n_per_ssu", rate_label)

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
