#' Constrained stratified allocation
#'
#' Distribute a total sample size across strata defined by a single
#' stratification variable, under a fixed total `n`, target CV, or budget.
#' When the design uses multiple stratification variables (e.g. region and
#' urbanicity), cross them into a single variable beforehand so that each row
#' of `frame` represents one unique stratum.
#'
#' @param frame For the default method: a **stratum-level** data frame
#'   describing the population you want to sample. Each row represents one
#'   stratum, a subgroup of the population defined by a stratification
#'   variable such as region, age group, or urbanicity. The values in this
#'   frame typically come from a census, a population register, or a
#'   previous survey. Any stratum table with the columns below works,
#'   including a pool summary carried back from an executed sample once the
#'   measure columns are added.
#'
#'   When a design stratifies by several variables at once (e.g. region
#'   \eqn{\times}{*} urbanicity), cross them into a single variable before
#'   calling `n_alloc` (e.g. with [interaction()]) so that each row maps
#'   to exactly one population cell.
#'
#'   When `measures` and `targets` are supplied, `frame` instead describes
#'   atomic allocation strata for joint constrained allocation. It requires
#'   unique, non-missing `stratum` and positive finite `N` columns. It may
#'   contain `unit_cost`, `max_weight`, `take_all`, and complete domain
#'   classification columns named by `targets$domain`. In that mode `sd`,
#'   `var`, `mean` and `p` belong in `measures` instead.
#'   Each row must be homogeneous for every
#'   domain classification used by a target. For a fixed-take two-stage
#'   design, also supply whole `N_psu` and `n_per_psu`, plus positive
#'   `cost_psu` and `cost_ssu`. A three-stage design additionally requires
#'   whole `N_ssu` and `n_per_ssu`, and positive `cost_tsu`. `N_psu` is the
#'   number of PSUs in the stratum. In three-stage mode, `N_ssu` is the
#'   aggregate number of SSUs across all PSUs in the stratum (not the number
#'   per PSU). Stage costs replace `unit_cost`.
#'
#'   **Required columns:**
#'   \describe{
#'     \item{`N`}{Number of units (e.g. households, individuals) in
#'       each stratum. These are population counts, not sample sizes.
#'       Must be positive and finite.}
#'     \item{`sd` **or** `var` (legacy mode)}{A measure of how spread out the variable
#'       of interest is within each stratum. Provide **exactly one**:
#'       \itemize{
#'         \item `sd`: the stratum standard deviation
#'           (\eqn{\sqrt{\text{variance}}}{sqrt(variance)}), or
#'         \item `var`: the stratum variance.
#'       }
#'       Both must be non-negative and finite. When all strata have
#'       equal variability (or variability is unknown), a constant
#'       column (e.g. `sd = 1`) yields proportional-to-size allocation.}
#'   }
#'
#'   **Optional columns:**
#'   \describe{
#'     \item{`stratum`}{A label identifying each stratum (e.g.
#'       `"Urban"`, `"Rural"`). If omitted, row numbers are used. Must
#'       be unique, or unique within each domain when `domains` is set.}
#'     \item{`mean` **or** `p`}{The stratum population mean or
#'       proportion of the variable of interest. **Required when solving
#'       for `cv`**, because the coefficient of variation is defined
#'       relative to the mean. Use `mean` for continuous variables
#'       and `p` (in \eqn{[0, 1]}) for binary (yes/no) variables.}
#'     \item{`unit_cost`}{Per-unit interviewing cost in each stratum
#'       (positive, finite). Set higher values for strata that are more
#'       expensive to reach. Defaults to 1 everywhere (equal cost).}
#'     \item{`max_weight`}{Maximum allowed sampling weight
#'       \eqn{N_h / n_h}. Caps how under-represented a stratum can be.
#'       Use `NA` for strata without a cap. It bounds the gross weight,
#'       the units drawn, not the response-adjusted \eqn{N_h/(n_hr_h)},
#'       which differs once `resp_rate` varies by stratum.}
#'     \item{`take_all`}{Logical (or 0/1). If `TRUE`, every unit in the
#'       stratum is included, a census stratum. Useful for small strata
#'       whose total population is tiny enough to enumerate.}
#'     \item{`deff`, `resp_rate`}{Per-stratum design effect and expected
#'       response rate, for designs whose fieldwork differs across
#'       strata. `NA` in either column falls back to the scalar argument.
#'       The `deff` and `resp_rate` arguments override these columns when
#'       given as vectors.}
#'   }
#'
#'   For `svyplan_prec` objects: a precision result from [prec_alloc()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param domains Character vector of column names in `frame` to treat as
#'   domain identifiers, or `NULL` (default) for no domains. All names
#'   must exist in `frame`. Domains define sub-populations that each
#'   contain one or more strata. When `cv` is the target, precision is
#'   enforced *within every domain* (see Details).
#' @param n Total sample size. Specify exactly one of `n`, `cv`, or `budget`.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.05` means the standard error of the estimated
#'   population mean or total should be at most 5 percent of the estimate.
#'   Requires `mean` or `p` in `frame`. When domain columns are
#'   present, this target is enforced in each domain. Specify exactly one
#'   of `n`, `cv`, or `budget`.
#' @param budget Total field budget. In the classic modes, specify exactly one
#'   of `n`, `cv`, or `budget`, and `budget` caps the affordable allocation.
#'   In joint constrained allocation it must be accompanied by `objective`:
#'   the allocation minimizes the priority-weighted
#'   relative variance of the objective components subject to the budget, every
#'   hard target in `targets`, and the stratum bounds. Supplying `budget`
#'   without `objective` in joint mode is an error, because the two readings
#'   are incompatible.
#' @param measures Optional long data frame for joint constrained allocation,
#'   with one row per required `stratum` and indicator `name`. Each row must
#'   contain either `p` in `[0, 1]`, or `mean` and exactly one of non-negative
#'   `sd` or `var`. Optional row-specific `deff` and `resp_rate` values override
#'   the scalar arguments, and `NA` uses the scalar default. In fixed-take
#'   multistage mode, `icc_psu` is required and `var_ratio_psu` defaults to 1.
#'   Three-stage mode also requires `icc_ssu` and derives `var_ratio_ssu` as
#'   `var_ratio_psu * (1 - icc_psu)` when it is absent, which is the value the
#'   variance decomposition implies.
#'   Stage parameters may instead be supplied as stratum defaults in `frame`,
#'   in which case non-missing measure-row values take precedence.
#'   Only rows selected by a
#'   target and its domain are value-validated. Unused rows are ignored after
#'   their non-empty, unique keys and frame-stratum membership are checked.
#'   Must be supplied with `targets`.
#' @param targets Optional long data frame with one precision requirement per
#'   row. It requires indicator `name` and exactly one positive `cv`, `moe`,
#'   or `rmoe`.
#'   Optional `domain` and `level` select a frame domain (`.overall` and
#'   `level = NA` denote the whole population), `alpha` overrides the scalar
#'   value for MOE, and `constraint` supplies a stable unique identifier.
#'   Must be supplied with `measures`. Optional when `objective` and `budget`
#'   are given, which requests the best design a budget can buy with no hard
#'   precision requirement.
#' @param psu Optional PSU register for certainty-aware allocation, one row
#'   per PSU, with `stratum` and `N`, the count of ultimate units in that PSU.
#'   The PSU sizes must sum to `frame$N` within each stratum.
#'   Optional `psu_id` is carried through and optional `certainty` forces a
#'   PSU into the certainty part whatever its size. It is an alternative to
#'   the frame column `N_psu`, which says only how many PSUs a stratum has,
#'   and requires the frame column `n_per_psu` and `measures$icc_psu`. Joint
#'   constrained allocation only, and two stages only. See the certainty
#'   section under Details.
#' @param objective Optional estimates whose priority-weighted relative
#'   variance is minimized among the allocations that meet `targets` and the
#'   budget. Either a character vector of indicator names (overall domain,
#'   equal priority) or a data frame with `name` plus optional `domain`,
#'   `level`, non-negative `priority`, and `component` columns, using the same
#'   indicator-domain identifiers as `targets`. At least one priority must be
#'   positive. Priorities need not be normalized. Requires `budget` and must be
#'   supplied with `measures`.
#' @param alloc Allocation rule: `"neyman"` (default), `"optimal"`,
#'   `"proportional"`, or `"power"`.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of
#'   per-stratum unit costs, overriding `frame$unit_cost`. Not used with
#'   fixed-take multistage joint allocation, which requires stage costs.
#' @param alpha Significance level, default 0.05.
#' @param deff Design effect multiplier (> 0). A scalar applies to every
#'   stratum. A length-`nrow(frame)` vector gives one per stratum and
#'   overrides a `deff` column in `frame`.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1. A scalar
#'   applies to every stratum. A length-`nrow(frame)` vector gives one per
#'   stratum and overrides a `resp_rate` column in `frame`. See the
#'   nonresponse section of [svyplan-package] for what this adjustment does
#'   and does not claim.
#' @param df Degrees of freedom of the variance estimator the allocation
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()]. It switches the quantile used to translate a `moe`
#'   target in and to report `moe` out. A `cv` target carries no quantile
#'   and is unaffected. `NULL` (default) applies no adjustment.
#' @param min_n_stratum Optional minimum sample size per stratum.
#' @param alloc_q Bankier power parameter, used only when `alloc = "power"`.
#'   Numeric scalar in \eqn{[0, 1]}. At `alloc_q = 1` the allocation equals
#'   Neyman. At `alloc_q = 0` it yields near-equal subnational CVs.
#'   Default 0.5.
#' @param fpc Which finite population correction the variance carries in
#'   cluster allocation: `"unit"` (default) for the ultimate-unit
#'   `1 - n / N` on the whole inflated variance, `"stage"` for the
#'   stage-by-stage correction, exact under the equal-size, equal-take
#'   variance-component model this path assumes, or `"none"` for
#'   [n_cluster()]'s with-replacement model. `"stage"` requires `N_psu`, and
#'   `N_ssu` as well at three stages. Refused outside cluster allocation and
#'   with a `psu` register. See the correction section in Details.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_n` object with `type = "alloc"` and a stratum-level
#'   allocation table in `$detail` (also available via `as.data.frame()`),
#'   with columns:
#'   \describe{
#'     \item{`stratum`, `N`, `sd`, `unit_cost`}{Stratum identifiers and inputs
#'       carried over from the frame. In cluster mode `unit_cost` is the
#'       derived effective per-element cost
#'       `cost_psu / n_per_psu + cost_ssu` (1 when no stage costs were
#'       given), while `sd` stays the input stratum SD.}
#'     \item{`n`}{Allocated sample size (continuous).}
#'     \item{`n_int`}{Integer allocation. In `n` mode the requested
#'       total is preserved (bounded largest-remainder rounding), except
#'       under a fixed cluster take, where the take is held and the total
#'       moves to the nearest whole multiple of it instead. In
#'       `cv` mode each stratum is rounded up so the integer design
#'       meets the target. In `budget` mode units are added by variance
#'       reduction per unit cost so the integer design stays within
#'       budget. Always inside the integerized bounds
#'       (`ceiling(.lower)`, `floor(.upper)`). An error is raised when
#'       no integer allocation can satisfy them.}
#'     \item{`weight`}{Design weight `N / n`.}
#'     \item{`n_eff`}{Effective sample size `n * resp_rate / deff`.}
#'     \item{`.lower`, `.upper`}{Bounds applied to the stratum
#'       (from `min_n_stratum`, `max_weight`, `take_all`, `N_psu`, or
#'       `N`).}
#'     \item{`.binding`}{Whether the allocation sits on one of its
#'       bounds.}
#'     \item{`.bound_source`}{Which constraint produced the bound the
#'       allocation sits on, and `NA` when it sits on none.}
#'     \item{`.psu_frac`}{Cluster mode with `N_psu` only: the share of
#'       the stratum's PSU universe the design selects, `n_psu / N_psu`.}
#'     \item{`.se`, `.moe`, `.rmoe`, `.cv`}{Precision of the stratum's own
#'       mean under the allocation, on the same scale as `mean` (or `p`).
#'       Multiply `.se` by `N` for the stratum total. `.rmoe` is `.moe`
#'       as a fraction of the stratum mean. `.cv` and `.rmoe` are `NA`
#'       unless the frame supplies `mean` or `p`.}
#'     \item{`.share`}{The stratum's share of the design variance of the
#'       overall mean. Sums to 1 and identifies where precision is
#'       actually being bought. It is `NA` when the design variance is
#'       zero, which a census is, since no stratum then holds a share of
#'       nothing.}
#'     \item{`take_all`, `mean`}{Present when take-all strata or stratum
#'       means were supplied.}
#'     \item{`n_per_psu`, `n_psu`}{Cluster mode only: the continuous
#'       per-stratum take and implied number of PSUs (`n / n_per_psu`).}
#'     \item{`n_psu_int`, `n_per_psu_int`}{Cluster mode only: the
#'       whole-unit field design (whole PSUs and whole takes), chosen so
#'       that `budget` designs stay within budget and `cv` designs meet
#'       the target. Here `n_int = n_psu_int * n_per_psu_int`.}
#'     \item{`n_psu`, `n_psu_int`, `n_per_psu`, `n_per_ssu`}{Joint fixed-take
#'       multistage mode: continuous and whole first-stage decisions and the
#'       fixed later-stage takes. `n_per_ssu` is present only for three stages.
#'       Public `n` remains in ultimate-unit units, so
#'       `n = n_psu * n_per_psu` for two stages and
#'       `n = n_psu * n_per_psu * n_per_ssu` for three stages. The same exact
#'       identities hold for `n_int` and `n_psu_int`.}
#'   }
#'
#'   The result also carries an `$operational` list describing the integer
#'   field design. In joint mode it contains total `n`, variable field `cost`,
#'   the full target table in `constraints`, `all_pass`, and the number of
#'   integer-repair additions. In legacy modes it also reports the applicable
#'   aggregate `se`, `moe`, and `cv`. Top-level quantities describe the
#'   continuous design, while `as.integer()` returns the operational total.
#'
#'   In fixed-budget objective mode the result adds `$objective`, one row per
#'   component with its `priority`, relative variance `.relvar`, equivalent
#'   `.cv`, weighted `.contribution`, and `.share` of the total, plus the scalar
#'   `$objective_value`. `$operational` carries the same two fields for the
#'   integer design together with `budget_residual`. `$params$mode` is
#'   `"budget_objective"` rather than `"targets"`, and `$optimization` adds the
#'   budget, whether it binds, its residual, the local `budget_sensitivity`
#'   (change in objective value per extra unit of budget), and the
#'   epsilon-constraint `objective_bound` and its multiplier.
#'
#' @details
#' ## Joint constrained allocation
#'
#' Supplying both `measures` and `targets` requests a minimum-cost allocation
#' meeting every indicator-domain precision requirement at once, and adding
#' `objective` with `budget` asks instead for the best design that budget can
#' buy. This generalized mode is mutually exclusive with scalar `n`, `cv`,
#' `budget`, `domains` and an explicit `alloc` rule, because the target table
#' defines the constraints and the allocation is solved jointly. Its input
#' contract, its multistage variance model, its infeasibility messages and the
#' `predict()` frontier are on [n_alloc-generalized].
#'
#' For fieldwork, use `$detail$n_int`, not the analytical continuous `n`.
#' Fixed-take count columns are numeric vectors whose values are exact whole
#' numbers. The `_int` suffix denotes operational integer semantics rather
#' than R's 32-bit `integer` storage type. Result lists contain diagnostics and
#' are not promised as a long-term serialization format.
#'
#' ## Building the frame
#'
#' The `frame` is a data frame where **each row is one stratum** of
#' your target population. It summarizes what you know about each
#' subgroup *before* sampling. A typical workflow:
#'
#' 1. **Identify strata** from a census or register (e.g. provinces,
#'    urban/rural areas, age groups).
#' 2. **Look up `N`**: the population count per stratum.
#' 3. **Estimate `sd`**: the standard deviation of your key variable
#'    within each stratum (from a pilot survey, a previous census, or
#'    expert judgement). If unknown, set `sd = 1` everywhere for
#'    proportional allocation.
#' 4. **Add `mean` or `p`** if you want to solve for a target CV.
#'
#' A minimal frame:
#'
#' ```
#' frame <- data.frame(
#'   stratum = c("Urban", "Rural"),
#'   N       = c(50000, 120000),
#'   sd      = c(12, 20)
#' )
#' ```
#'
#' When a design stratifies by several variables (e.g. region
#' \eqn{\times}{*} urbanicity), cross them into one variable first:
#'
#' ```
#' frame$stratum <- interaction(frame$region, frame$urban, drop = TRUE)
#' ```
#'
#' This ensures that each row maps to exactly one population cell and that
#' the allocation formulas apply to the correct per-stratum `N` and `sd` pairs.
#'
#' ## Cluster designs within strata
#'
#' Adding a `icc_psu` column to the frame turns the allocation into a
#' stratified **two-stage** design (e.g. enumeration areas then
#' households within each stratum). Under the cluster variance model
#' the problem reduces to the element allocation above with the stratum
#' SD inflated to
#' `sd * sqrt(var_ratio_psu * (1 + icc_psu * (n_per_psu * resp_rate - 1)))`
#' and, when stage costs are given, a per-element cost of
#' `cost_psu / n_per_psu + cost_ssu`. The clustering penalty is paid on the
#' take that responds, `n_per_psu * resp_rate`, since a unit that does not
#' respond contributes no within-cluster observation. The cost is paid on
#' the gross take, since it is issued either way. The two coincide at
#' `resp_rate = 1`. The whole-cluster operational search reads the same
#' responding take, so a design it accepts is one the continuous reduction
#' also accepts. All solve modes, allocation
#' methods, and constraints work unchanged. The `n`, `cv`, and `budget` modes
#' keep their meanings.
#'
#' Cluster-mode columns:
#' - `icc_psu` (required): within-PSU homogeneity per stratum,
#'   e.g. from [varcomp()] with `strata`.
#' - `var_ratio_psu` (optional, default 1): variance ratio per stratum.
#' - `n_per_psu` (optional): fixes the per-stratum take, in the
#'   operational whole-unit design as much as in the continuous one. Any
#'   `NA` entries are replaced by the cost-optimal take
#'   `sqrt(cost_psu / cost_ssu * (1 - icc_psu) / icc_psu)`. Holding the
#'   take is what `n` mode gives up its exact total for. The stratum is
#'   fielded as a whole number of clusters of that size, so its element
#'   count lands on the nearest multiple.
#' - `cost_psu`, `cost_ssu` (together): per-PSU and per-element costs.
#'   These are required for `budget` mode or when `n_per_psu` is not fixed. They
#'   replace `unit_cost`, which is not allowed in this mode.
#' - `N_psu` (optional): the number of PSUs available in the stratum.
#'
#' ## What `N_psu` does, and what it does not
#'
#' `N_psu` is a feasibility constraint and nothing more. It caps the
#' allocation at `N_psu * n_per_psu` ultimate units and caps the whole-unit
#' design at `N_psu` clusters, matching the bound the fixed-take path
#' already applies. Leaving it out preserves the unbounded behaviour, in
#' which the allocation may ask for more PSUs than a stratum contains.
#'
#' On its own it does **not** activate a first-stage finite population
#' correction. Under the default `fpc = "unit"` no factor keyed to the PSU
#' sampling fraction enters the variance, however large a share of the PSU
#' universe the design takes, and `print()` says so once that share is
#' appreciable. Supplying `N_psu` bounds the design. Asking for the
#' correction as well is what `fpc = "stage"` does.
#'
#' So at `n_psu == N_psu` the default still carries between-PSU variance even
#' though every PSU has been selected. That is deliberately conservative
#' and is no longer a literal variance representation. For the same
#' reason, a target this path reports as unreachable at the PSU bound may
#' be reachable under a finite-population first stage, and the error says so
#' rather than claiming the precision is impossible. Under `fpc = "stage"`
#' the between-PSU component does reach zero there, and the whole variance
#' reaches zero only when every stage is enumerated.
#'
#' Because taking every PSU leaves the within-PSU take in force, it does not
#' enumerate a stratum, and `take_all` is refused in cluster mode for the
#' same reason the fixed-take path refuses it.
#'
#' ## Which correction the variance carries
#'
#' `fpc` chooses, and the three choices are nested rather than arbitrary.
#' Write \eqn{f_1} for the PSU sampling fraction and \eqn{f_2} for the
#' within-PSU one, both measured on what responds. Under the equal-size,
#' equal-take variance-component model this path assumes, the two-stage
#' without-replacement variance corrects each component by its own stage:
#'
#' \deqn{V = (1-f_1)\frac{S_1^2}{a} + (1-f_2)\frac{S_2^2}{am}.}{V = (1-f_1) S_1^2/a + (1-f_2) S_2^2/(am).}
#'
#' - `"unit"` (default) applies the ultimate-unit `1 - n / N` to the whole
#'   clustering-inflated variance, so the between-PSU component is attenuated
#'   by it too. Since \eqn{n/N = f_1f_2}, that factor is at least as large as
#'   either exact one, so this **overstates both components** rather than
#'   either being understated. It is conservative against the exact
#'   without-replacement variance at every sampling fraction.
#' - `"stage"` applies the expression above, and its three-stage analogue.
#'   It needs `N_psu`, and `N_ssu` as well at three stages. It is exact for
#'   that model and not beyond it: unequal PSU sizes, a size-proportional
#'   first stage, and stochastic rather than expected nonresponse each need
#'   a more general variance expression, so it is refused with a `psu`
#'   register.
#' - `"none"` applies no correction at any stage, which is [n_cluster()]'s
#'   model exactly. Use it to compare the two interfaces at a sampling
#'   fraction where the correction would otherwise separate them.
#'
#' So `none` \eqn{\ge} `unit` \eqn{\ge} `stage` at every sampling fraction,
#' and the three coincide as the fractions vanish.
#'
#' The default stays `"unit"` deliberately. `"stage"` is exact under that
#' model, but it buys the exactness by depending on the mean PSU size
#' \eqn{M = N/N_{psu}}{M = N/N_psu}, which the with-replacement form does not
#' use at all. Under unequal PSU sizes that mean is a stand-in, and planners
#' usually know it less well than they know `icc_psu` and the take. Choose
#' `"stage"` when the PSU sizes really are close to equal and the fractions
#' are large enough to matter.
#'
#' Note that the two fractions are controlled by different quantities.
#' \eqn{f_1} falls as the design takes a smaller share of the PSU universe,
#' but \eqn{f_2} is the realized take over \eqn{M} and does not move with `N`
#' at all. A large population reached through a few large PSUs leaves
#' \eqn{f_2} wherever it was, so `"stage"` is not in general a small
#' correction to `"none"`.
#'
#' `fpc` is refused outside cluster allocation, where there are no stages to
#' choose between and the correction is always `1 - n / N`.
#'
#' Because `icc_psu` already accounts for the clustering, leave
#' `deff` at 1 unless it captures a *different* source of design
#' effect (e.g. weighting loss). A clustering `deff` on top of
#' `icc_psu` would double-count. The constraints `min_n_stratum` and
#' `max_weight` stay in element units. For fielding,
#' use the whole-unit design in `n_psu_int` and `n_per_psu_int`
#' (`n_int = n_psu_int * n_per_psu_int`). Its actual field cost and
#' precision are reported in `$operational` and, in `budget` mode,
#' never exceed the budget.
#'
#' ## Certainty PSUs from a register
#'
#' `N_psu` gives only the number of PSUs in a stratum. Supply `psu` to give
#' their sizes instead. It has one row per PSU with `stratum` and `N`, and the
#' sizes must sum to `frame$N` in every stratum. This two-stage joint-allocation
#' mode requires `n_per_psu` in `frame` and `icc_psu` in `measures`.
#'
#' With take \eqn{b_h} and sampling fraction \eqn{f_h = n_h / N_h}, a PSU is
#' certain when its size reaches
#'
#' \deqn{N_{hi} \ge b_h / f_h.}{N_hi >= b_h / f_h.}
#'
#' The certainty part has no first-stage sampling variance. The remainder has
#' the usual clustering component, and the two are combined into the
#' anticipated design effect. `psu$certainty` can add a PSU to the certainty
#' part. It cannot remove a PSU above the threshold.
#'
#' The threshold and allocation determine each other, so the solver iterates.
#' `$optimization$certainty` records whether it converged, cycled, or reached
#' its iteration limit. `$detail` gives the split and threshold by stratum,
#' while `$psu` gives each PSU's classification, its source, and `n_take`,
#' the take the operational design fields in it. A certainty PSU carries its
#' own whole take at the stratum rate, capped at its size, and a remainder
#' PSU carries `n_per_psu`. The takes are the operational design's own
#' numbers, so `sum(n_take[certainty]) = n_certain_int` in every stratum.
#'
#' A fit solved with a register keeps it in `$params$psu`. Its presence is
#' the stable test for a certainty-aware fit, and it's possible to key on it
#' rather than on `$detail` column names.
#'
#' The operational design is fieldable: certainty PSUs use their whole take
#' at the stratum rate and the remainder uses whole PSUs at `n_per_psu`. Thus
#' `n_int = n_certain_int + n_psu_draw * n_per_psu` in every stratum.
#' `$operational` reports the count, cost, and precision of that same design.
#' Use `cost_psu` and `cost_ssu` together to price PSU visits and interviews.
#' [predict.svyplan] can compare fixed values of `n_per_psu`.
#'
#' ## Domains vs. strata
#'
#' Domains are specified via the `domains` parameter. Domain columns
#' partition strata into sub-populations. Each domain groups
#' one or more strata. When `cv` is specified, the algorithm finds the
#' minimum total \eqn{n} such that the *worst-case* domain CV meets the
#' target, i.e. every domain achieves the required precision.
#'
#' In `n` or `budget` mode, domains affect reporting only: per-domain
#' precision metrics appear in `$domains` but the allocation itself treats
#' all strata globally.
#'
#' ## Allocation methods
#'
#' Allocation is controlled by the `alloc` parameter (same methods as
#' [strata_bound()]):
#' - **proportional**: \eqn{n_h \propto N_h / r_h}{n_h proportional to N_h / r_h}
#' - **neyman**: \eqn{n_h \propto N_h S_h \sqrt{d_h / r_h}}{n_h proportional to N_h S_h sqrt(d_h / r_h)}
#' - **optimal**: \eqn{n_h \propto N_h S_h \sqrt{d_h / r_h} / \sqrt{c_h}}{n_h proportional to N_h S_h sqrt(d_h / r_h) / sqrt(c_h)}
#' - **power**: Bankier (1988), \eqn{n_h \propto S_h N_h^{q}\sqrt{d_h/r_h}}{n_h ~ S_h * N_h^alloc_q * sqrt(deff_h / resp_rate_h)},
#'   with exponent `alloc_q`
#'
#' The design effect \eqn{d_h} and response rate \eqn{r_h} enter only when
#' they vary by stratum. A value shared by every stratum is a constant
#' factor and cancels out of a proportional weighting, leaving the classical
#' rules above it. The three variance-based rules carry
#' \eqn{\sqrt{d_h/r_h}}{sqrt(d_h/r_h)} because they minimize
#' \eqn{\sum W_h^2S_h^2d_h/(r_hn_h)} against a constraint on the units
#' *drawn*, which is what a budget pays for. Proportional carries
#' \eqn{1/r_h} instead, and no \eqn{d_h} at all, because it is a count rule, not a
#' variance optimum, and its purpose is a self-weighting sample, so it is
#' the responding sample it holds proportional to \eqn{N_h}.
#'
#' Stratum allocations are rounded to integers using the ORIC method
#' (Cont and Heidari, 2015). Constraints (`min_n_stratum`, `max_weight`, `take_all`)
#' are enforced via recursive Neyman allocation (RNA, Wesolowski et al., 2021).
#'
#' When `budget` is specified, the algorithm finds the maximum affordable
#' allocation under unit costs.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Chapters 5
#'   and 9.
#'
#' Bankier, M. D. (1988). Power allocations: determining sample sizes for
#'   subnational areas. *The American Statistician*, 42(3), 174--177.
#'
#' Bethel, J. (1989). Sample allocation in multivariate surveys.
#'   *Survey Methodology*, 15(1), 47--57.
#'
#' Cont, R. and Heidari, M. (2015). Optimal rounding under integer
#'   constraints. *arXiv preprint* arXiv:1501.00014.
#'
#' Wesolowski, J., Wieczorkowski, R. and Wojciak, W. (2021). Optimality of
#'   the recursive Neyman allocation. *Journal of Survey Statistics and
#'   Methodology*, 10(5), 1263--1275.
#'
#' @family stratified design functions
#' @seealso [n_alloc-generalized] for the joint mode in full,
#'   [prec_alloc()] for the inverse, [strata_bound()] for constructing the
#'   strata to allocate over, and [n_multi()] for the unstratified
#'   multi-indicator size that joint allocation refines.
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
#'
#' frame_constraints <- transform(
#'   frame,
#'   max_weight = c(25, 20, NA),
#'   take_all = c(FALSE, FALSE, TRUE)
#' )
#'
#' n_alloc(frame_constraints, budget = 3500, alloc = "optimal", min_n_stratum = 40)
#'
#' frame_domains <- data.frame(
#'   province = c("North", "North", "South", "South"),
#'   stratum = c("Urban", "Rural", "Urban", "Rural"),
#'   N    = c(2000, 3000, 1800, 3200),
#'   sd   = c(12, 18, 10, 16),
#'   mean = c(55, 48, 58, 50)
#' )
#'
#' n_alloc(frame_domains, domains = "province",
#'        cv = 0.04, alloc = "power", alloc_q = 0.3)
#'
#' # Joint targets for two indicators and an overlapping domain
#' joint_frame <- data.frame(
#'   stratum = c("North urban", "North rural", "South urban"),
#'   region = c("North", "North", "South"),
#'   N = c(1000, 1800, 1200)
#' )
#' joint_measures <- data.frame(
#'   stratum = rep(joint_frame$stratum, 2),
#'   name = rep(c("coverage", "income"), each = 3),
#'   p = c(0.5, 0.4, 0.6, NA, NA, NA),
#'   mean = c(NA, NA, NA, 50, 55, 60),
#'   sd = c(NA, NA, NA, 10, 12, 15)
#' )
#' joint_targets <- data.frame(
#'   name = c("coverage", "income"),
#'   domain = c(".overall", "region"),
#'   level = c(NA, "North"),
#'   cv = c(0.06, NA),
#'   moe = c(NA, 3)
#' )
#' joint <- n_alloc(
#'   joint_frame, measures = joint_measures, targets = joint_targets
#' )
#' prec_alloc(joint, n = joint$detail$n_int)
#'
#' # The best design a fixed budget can buy, subject to the same targets
#' bought <- n_alloc(
#'   joint_frame, measures = joint_measures, targets = joint_targets,
#'   objective = "income", budget = 4000
#' )
#' bought$objective
#'
#' # Unequal priorities, and an objective on a specific domain
#' priorities <- data.frame(
#'   name     = c("income", "coverage"),
#'   domain   = c(".overall", "region"),
#'   level    = c(NA, "North"),
#'   priority = c(2, 1)
#' )
#' n_alloc(
#'   joint_frame, measures = joint_measures, targets = joint_targets,
#'   objective = priorities, budget = 4000
#' )
#'
#' # What would other budgets buy? (cost-versus-objective frontier)
#' predict(bought, data.frame(budget = c(3000, 4000, 5000)))
#'
#' # The same target tables with fixed two-stage takes
#' joint_frame$N_psu <- c(100, 150, 100)
#' joint_frame$n_per_psu <- c(8, 10, 12)
#' joint_frame$cost_psu <- c(300, 400, 450)
#' joint_frame$cost_ssu <- c(25, 30, 35)
#' joint_measures$icc_psu <- rep(c(0.03, 0.05, 0.08), 2)
#' joint_measures$var_ratio_psu <- 1
#' joint_cluster <- n_alloc(
#'   joint_frame, measures = joint_measures, targets = joint_targets
#' )
#' joint_cluster$detail[, c("stratum", "n_psu_int", "n_per_psu", "n_int")]
#'
#' # Fixed three-stage takes: PSUs, SSUs within PSU, ultimate units within SSU
#' joint_frame$N_ssu <- c(600, 1000, 750) # aggregate SSUs by stratum
#' joint_frame$n_per_ssu <- c(3, 4, 3)
#' joint_frame$cost_tsu <- c(5, 6, 7)
#' joint_measures$icc_ssu <- rep(c(0.10, 0.08, 0.12), 2)
#' joint_measures$var_ratio_ssu <- 1
#' joint_three <- n_alloc(
#'   joint_frame, measures = joint_measures, targets = joint_targets,
#'   min_n_stratum = 20
#' )
#' joint_three$detail[, c(
#'   "stratum", "n_psu_int", "n_per_psu", "n_per_ssu", "n_int"
#' )]
#' prec_alloc(joint_three, n = joint_three$detail$n_int)
#'
#' # Stratified two-stage design (EAs then households per stratum)
#' frame_cluster <- data.frame(
#'   stratum   = c("Urban", "Rural"),
#'   N         = c(50000, 150000),
#'   sd        = c(0.45, 0.48),
#'   mean      = c(0.35, 0.25),
#'   icc_psu = c(0.03, 0.08),
#'   cost_psu  = c(300, 600),
#'   cost_ssu  = c(40, 60)
#' )
#'
#' n_alloc(frame_cluster, cv = 0.05)
#'
#' @export
n_alloc <- function(frame, ...) {
  if (!missing(frame)) {
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
  budget = NULL,
  measures = NULL,
  targets = NULL,
  psu = NULL,
  objective = NULL,
  alloc = c("neyman", "optimal", "proportional", "power"),
  unit_cost = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  df = NULL,
  min_n_stratum = NULL,
  alloc_q = 0.5,
  fpc = c("unit", "stage", "none"),
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_alloc.default, match.call(), environment())
  if (!is.null(.plan)) {
    return(do.call(n_alloc.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  alloc_default <- c("neyman", "optimal", "proportional", "power")
  alloc_explicit <- !missing(alloc) && !identical(alloc, alloc_default)
  alloc_q_explicit <- !missing(alloc_q) && !identical(alloc_q, 0.5)
  fpc <- match.arg(fpc)
  joint_any <- !is.null(measures) || !is.null(targets) || !is.null(objective)
  if (joint_any) {
    if (is.null(measures)) {
      stop("joint constrained allocation requires 'measures'", call. = FALSE)
    }
    if (!is.null(n) || !is.null(cv)) {
      stop(
        "joint constrained allocation cannot be combined with 'n' or 'cv'",
        call. = FALSE
      )
    }
    if (is.null(objective)) {
      if (is.null(targets)) {
        stop(
          "'measures' and 'targets' must be supplied together",
          call. = FALSE
        )
      }
      if (!is.null(budget)) {
        # 'budget' means "ceiling on the cheapest feasible design" in classic
        # allocation and "spend this" in objective mode. A silent third
        # meaning would be worse than an error.
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
    if (!is.null(domains)) {
      stop(
        "use 'targets$domain' for joint constrained allocation, not 'domains'",
        call. = FALSE
      )
    }
    if (alloc_explicit) {
      stop(
        "'alloc' is not used for joint constrained allocation, it is determined by 'targets'",
        call. = FALSE
      )
    }
    if (alloc_q_explicit) {
      stop(
        "'alloc_q' is not used for joint constrained allocation",
        call. = FALSE
      )
    }
    if (!is.null(psu)) {
      # A register is not the equal-size, equal-take design the stage
      # correction is derived for.
      if (!identical(fpc, "unit")) {
        stop(
          "'fpc' is not available with a PSU register: the certainty split and its size-proportional remainder need their own variance decomposition, which is not derived here. Drop 'fpc', or drop 'psu' and give 'N_psu' instead",
          call. = FALSE
        )
      }
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
        df = df
      )))
    }
    return(.n_alloc_bethel(
      frame = frame,
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
      fpc = fpc
    ))
  }
  if (!is.null(psu)) {
    stop(
      "'psu' requires a joint constrained allocation: supply 'measures' and 'targets'",
      call. = FALSE
    )
  }
  alloc <- match.arg(alloc)
  check_alpha(alpha)

  if (!is.null(min_n_stratum)) {
    check_scalar(min_n_stratum, "min_n_stratum")
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

  mode_count <- (!is.null(n)) + (!is.null(cv)) + (!is.null(budget))
  if (mode_count != 1L) {
    stop("specify exactly one of 'n', 'cv', or 'budget'", call. = FALSE)
  }

  prep <- .alloc_prepare_frame(
    frame,
    domains = domains,
    unit_cost = unit_cost,
    deff = deff,
    resp_rate = resp_rate,
    fpc = fpc
  )
  N_h <- prep$N_h
  S_h <- prep$S_h
  mean_h <- prep$mean_h
  cost_h <- prep$cost_h
  deff_arg <- deff
  resp_rate_arg <- resp_rate
  deff <- prep$deff_h
  resp_rate <- prep$resp_rate_h

  bounds <- .alloc_bounds(
    N_h = N_h,
    max_weight = prep$max_weight,
    take_all = prep$take_all,
    min_n_stratum = min_n_stratum,
    psu_cap = if (!is.null(prep$N_psu_h)) prep$N_psu_h * prep$n_per_psu_h
  )
  m_h <- bounds$m_h
  M_h <- bounds$M_h
  prep$m_src <- bounds$m_src
  prep$M_src <- bounds$M_src
  lo_i <- as.integer(ceiling(m_h - 1e-9))
  hi_i <- as.integer(floor(M_h + 1e-9))
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

  a_h <- .alloc_weights(alloc, alloc_q, N_h, S_h, cost_h, deff, resp_rate)
  if (!is.finite(sum(a_h)) || sum(a_h) <= 0) {
    a_h <- N_h
  }

  mode <- if (!is.null(n)) {
    "n"
  } else if (!is.null(cv)) {
    "cv"
  } else {
    "budget"
  }
  if (
    mode == "budget" && isTRUE(prep$cluster) && !isTRUE(prep$has_stage_costs)
  ) {
    stop(
      "'budget' mode with 'icc_psu' requires 'cost_psu' and 'cost_ssu' columns",
      call. = FALSE
    )
  }
  target_total <- NA_real_

  if (mode == "n") {
    check_scalar(n, "n")
    if (n < lo - tol) {
      stop(
        "'n' is below the minimum feasible total under constraints",
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
          "no integer allocation reaches total %d within the integer bounds (feasible totals: %d to %d)",
          as.integer(round(n)),
          sum(lo_i),
          sum(hi_i)
        ),
        call. = FALSE
      )
    }
    target_total <- n
  } else if (mode == "cv") {
    check_scalar(cv, "cv")
    if (anyNA(mean_h)) {
      stop(
        "'mean' (or 'p') is required in frame when solving for 'cv'",
        call. = FALSE
      )
    }

    cv_for_total <- function(n_total) {
      n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
      if (length(prep$domain_idx) == 0L) {
        .alloc_metrics(
          N_h = N_h,
          S_h = S_h,
          mean_h = mean_h,
          n_h = n_h,
          alpha = alpha,
          deff = deff,
          resp_rate = resp_rate,
          cost_h = cost_h,
          N_fpc = prep$N_fpc_h
        )$cv
      } else {
        .alloc_domain_cv_max(
          prep = prep,
          n_h = n_h,
          alpha = alpha,
          deff = deff,
          resp_rate = resp_rate
        )
      }
    }

    cv_lo <- cv_for_total(lo)
    if (cv <= 0) {
      stop("'cv' must be positive", call. = FALSE)
    }
    if (cv_lo <= cv + tol) {
      target_total <- lo
    } else {
      cv_hi <- cv_for_total(hi)
      if (!is.finite(cv_hi) || cv_hi > cv + tol) {
        stop(
          sprintf(
            "target 'cv' is unattainable under current constraints%s",
            .alloc_psu_bound_note(prep)
          ),
          call. = FALSE
        )
      }
      target_total <- uniroot(
        function(x) cv_for_total(x) - cv,
        interval = c(lo, hi),
        tol = 1e-8
      )$root
    }
  } else {
    check_scalar(budget, "budget")
    cost_lo <- sum(m_h * cost_h)
    cost_hi <- sum(M_h * cost_h)
    if (budget < sum(lo_i * cost_h) - tol) {
      stop(
        sprintf(
          "'budget' cannot fund the integer lower bounds (minimum integer cost = %.4g)",
          sum(lo_i * cost_h)
        ),
        call. = FALSE
      )
    }

    if (budget < cost_lo - tol) {
      stop(
        "'budget' is below the minimum feasible cost under constraints",
        call. = FALSE
      )
    }
    if (budget > cost_hi + tol) {
      stop(
        sprintf(
          "'budget' exceeds the maximum feasible cost (census bound)%s",
          .alloc_psu_bound_note(prep)
        ),
        call. = FALSE
      )
    }

    if (abs(budget - cost_lo) <= tol) {
      target_total <- lo
    } else {
      cost_for_total <- function(n_total) {
        n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
        sum(n_h * cost_h)
      }
      target_total <- uniroot(
        function(x) cost_for_total(x) - budget,
        interval = c(lo, hi),
        tol = 1e-8
      )$root
    }
  }

  n_h <- .rna_alloc(a_h, target_total, m_h, M_h)
  metrics <- .alloc_metrics(
    N_h = N_h,
    S_h = S_h,
    mean_h = mean_h,
    n_h = n_h,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = cost_h,
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
    mode = mode,
    alloc_q = if (alloc == "power") alloc_q else NULL,
    n_h = n_h,
    achieved = list(n = sum(n_h), cv = metrics$cv, cost = metrics$cost)
  )
  params[[mode]] <- switch(mode, n = n, cv = cv, budget = budget)

  if (isTRUE(prep$cluster)) {
    opc <- .alloc_operational_cluster(
      prep,
      detail,
      n_h,
      mode,
      budget,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      target_cv = if (mode == "cv") cv else NULL,
      df = df
    )
    detail <- opc$detail
    operational <- opc$operational
  } else {
    operational <- .alloc_operational_element(
      prep,
      detail,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      df = df
    )
  }

  obj <- .new_svyplan_n(
    n = sum(n_h),
    type = "alloc",
    method = alloc,
    params = params,
    se = metrics$se,
    moe = metrics$moe,
    cv = metrics$cv,
    detail = detail,
    binding = .alloc_binding_label(n_h, m_h, M_h),
    domains = domain_summary,
    operational = operational
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
  budget = NULL
) {
  x <- frame
  if (x$type != "alloc") {
    stop("n_alloc requires a svyplan_prec of type 'alloc'", call. = FALSE)
  }
  p <- x$params

  if (identical(x$method, "bethel")) {
    if (!is.null(n) || !is.null(cv) || !is.null(budget)) {
      stop(
        "generalized precision results cannot be combined with scalar 'n', 'cv', or 'budget'",
        call. = FALSE
      )
    }
    targets <- p$targets
    if (is.null(p$objective)) {
      # Each row is re-pinned in the units it was stated in.
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
      # Inverts through cost: pinning precision would over-constrain a design
      # already spending its budget.
      budget_arg <- p$achieved$cost
      if (nrow(targets) == 0L) targets <- NULL
    }
    args <- list(
      frame = p$frame,
      measures = p$measures,
      targets = targets,
      objective = p$objective,
      budget = budget_arg,
      unit_cost = p$unit_cost,
      alpha = p$alpha,
      deff = p$deff,
      resp_rate = p$resp_rate,
      df = p$df,
      min_n_stratum = p$min_n_stratum,
      fpc = p$fpc %||% "unit"
    )
    return(do.call(
      n_alloc.default,
      .roundtrip_args(args, list(...), n_alloc.default)
    ))
  }

  if (is.null(n) && is.null(cv) && is.null(budget)) {
    n <- p$achieved$n
  }

  args <- list(
    frame = p$frame,
    domains = p$domain_cols,
    n = n,
    cv = cv,
    budget = budget,
    alloc = p$alloc %||% "neyman",
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    df = p$df,
    min_n_stratum = p$min_n_stratum,
    alloc_q = p$alloc_q %||% 0.5
  )
  if (!.alloc_is_cluster(p$frame)) {
    args$unit_cost <- p$cost_h
  } else {
    args$fpc <- p$fpc %||% "unit"
  }
  do.call(n_alloc.default, .roundtrip_args(args, list(...), n_alloc.default))
}

#' Precision for a constrained allocation
#'
#' Compute aggregate precision for a stratum allocation. For a joint
#' constrained allocation, return one achieved-precision row per target and
#' separate stratum-bound diagnostics.
#'
#' @param frame For the default method: a stratum-level data frame in the same
#'   format as the `frame` argument to [n_alloc()] (one row per stratum,
#'   with at least `N` and `sd` or `var` columns). See [n_alloc()] for
#'   the full column reference).
#'   For `svyplan_n` objects: an allocation result from [n_alloc()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param n Stratum sample sizes, length `nrow(frame)`. For a fitted joint
#'   allocation, omission uses its continuous allocation. Pass `$detail$n_int`
#'   to assess the operational recommendation. A named vector is matched to
#'   `frame$stratum`, an unnamed vector is positional. In fixed-take
#'   multistage mode these are ultimate-unit sizes. Explicit adopted sizes must
#'   correspond to whole PSU counts. Omission from a fitted result retains its
#'   continuous first-stage allocation.
#' @param measures Optional long indicator table for joint assessment. See the
#'   `measures` argument to [n_alloc()]. It must be supplied with `targets` in
#'   the default method and is recovered automatically from a fitted result.
#' @param targets Optional long precision-requirement table for joint
#'   assessment. See the `targets` argument to [n_alloc()]. It must be supplied
#'   with `measures` in the default method and is recovered automatically from
#'   a fitted result.
#' @param psu Optional PSU register for certainty-aware assessment. See the
#'   `psu` argument to [n_alloc()]. No loop is needed here, since the
#'   allocation is supplied, so the threshold it implies is supplied with it and the
#'   classification is read off the design being assessed. The `$psu` table
#'   carries the same columns as a fitted plan's, with `n_take` read off the
#'   supplied allocation. A fitted result
#'   carries its own register and its held classification, so
#'   `prec_alloc(fit)` reproduces the plan's precision exactly.
#' @param objective Optional objective components to report alongside the
#'   targets. See the `objective` argument to [n_alloc()]. Recovered
#'   automatically from a fitted budget-objective result.
#' @param budget Optional budget, used only to report the residual against the
#'   assessed allocation's cost. Requires `objective`.
#' @param domains Character vector of column names in `frame` to treat as
#'   domain identifiers, or `NULL` (default) for no domains.
#' @param alpha Significance level, default 0.05.
#' @param deff Design effect multiplier (> 0). A scalar applies to every
#'   stratum. A length-`nrow(frame)` vector gives one per stratum and
#'   overrides a `deff` column in `frame`.
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1. A scalar
#'   applies to every stratum. A length-`nrow(frame)` vector gives one per
#'   stratum and overrides a `resp_rate` column in `frame`.
#' @param df Degrees of freedom of the variance estimator the allocation
#'   will have, typically sampled PSUs minus strata, and available from
#'   [design_df()]. It switches the quantile used to translate a `moe`
#'   target in and to report `moe` out. A `cv` target carries no quantile
#'   and is unaffected. `NULL` (default) applies no adjustment.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of
#'   per-stratum unit costs, overriding `frame$unit_cost`. Fixed-take
#'   multistage joint assessment instead uses the stage costs stored in
#'   `frame`.
#' @param fpc Which finite population correction the variance carries in
#'   cluster allocation: `"unit"` (default) for the ultimate-unit
#'   `1 - n / N` on the whole inflated variance, `"stage"` for the
#'   stage-by-stage correction, exact under the equal-size, equal-take
#'   variance-component model this path assumes, or `"none"` for
#'   [n_cluster()]'s with-replacement model. `"stage"` requires `N_psu`, and
#'   `N_ssu` as well at three stages. Refused outside cluster allocation and
#'   with a `psu` register. See the correction section in Details.
#' @param min_n_stratum Optional minimum sample size per stratum, applied as the
#'   lower bound the assessment reports against in `$bounds`. It is the same
#'   argument [n_alloc()] takes, so a design and its assessment can be held
#'   to one floor. Joint assessment only: supplying it without `measures` and
#'   `targets` is an error.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "alloc"`. Top-level `se`,
#'   `moe`, and `cv` describe the whole population. `$detail` carries the
#'   stratum table documented in [n_alloc()], including per-stratum `.se`,
#'   `.moe`, `.rmoe`, `.cv`, and the variance `.share`. When `domains` is
#'   given, `$domains` reports `.n`, `.se`, `.moe`, `.rmoe`, `.cv`, and
#'   `.cost` per domain,
#'   the same table [n_alloc()] returns, so a design and its assessment can
#'   be compared row for row.
#'
#'   For joint
#'   assessment, `$detail` is the constraint dictionary described in
#'   [n_alloc()]: target and achieved precision, ratio/residual/tolerance,
#'   pass/binding flags, and (when available) multiplier/sensitivity columns.
#'   `$bounds` reports ultimate-unit `n`, lower/upper limits, each violation
#'   flag, and `.pass` by stratum. `$params$achieved$cost` is variable field
#'   cost under the supplied allocation. When the fitted object carries an
#'   objective, `$objective` and `$objective_value` report its components and
#'   weighted value under the assessed allocation.
#'
#'   Passing a joint precision result back to [n_alloc()] round trips the
#'   design. Minimum-cost results invert through precision, pinning the
#'   achieved values as the requirement. Budget-objective results invert
#'   through cost instead. The targets stay as specified and the assessed
#'   allocation's cost becomes the budget, because pinning achieved precision
#'   as hard targets would over-constrain a design that already spends its
#'   whole budget.
#'
#' @family stratified design functions
#' @seealso [n_alloc()], and [n_alloc-generalized] for the joint mode.
#'
#' @examples
#' frame <- data.frame(
#'   N    = c(4000, 3000, 3000),
#'   sd   = c(10, 15, 8),
#'   mean = c(50, 60, 55)
#' )
#' res <- n_alloc(frame, n = 600)
#' prec_alloc(res)
#'
#' # Direct joint assessment (the allocation is in ultimate-unit units)
#' jf <- data.frame(
#'   stratum = c("A", "B"), region = c("North", "South"),
#'   N = c(1000, 1500), unit_cost = c(1, 1.2)
#' )
#' jm <- data.frame(
#'   stratum = rep(jf$stratum, 2),
#'   name = rep(c("coverage", "income"), each = 2),
#'   p = c(0.4, 0.6, NA, NA),
#'   mean = c(NA, NA, 50, 55),
#'   sd = c(NA, NA, 10, 12)
#' )
#' jt <- data.frame(
#'   name = c("coverage", "income"),
#'   cv = c(0.10, NA), moe = c(NA, 2.5)
#' )
#' prec_alloc(jf, n = c(120, 150), measures = jm, targets = jt)
#'
#' @export
prec_alloc <- function(frame, ...) {
  if (!missing(frame)) {
    .res <- .dispatch_plan(frame, "frame", prec_alloc.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_alloc")
}

#' @rdname prec_alloc
#' @export
prec_alloc.default <- function(
  frame,
  n = NULL,
  ...,
  measures = NULL,
  targets = NULL,
  psu = NULL,
  objective = NULL,
  budget = NULL,
  domains = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  df = NULL,
  unit_cost = NULL,
  min_n_stratum = NULL,
  fpc = c("unit", "stage", "none"),
  plan = NULL
) {
  .plan <- .merge_plan_args(
    plan,
    prec_alloc.default,
    match.call(),
    environment()
  )
  if (!is.null(.plan)) {
    return(do.call(prec_alloc.default, c(.plan, list(...))))
  }
  .check_unused_dots(...)
  fpc <- match.arg(fpc)
  joint_any <- !is.null(measures) || !is.null(targets) || !is.null(objective)
  if (joint_any) {
    if (is.null(measures)) {
      stop("joint precision assessment requires 'measures'", call. = FALSE)
    }
    if (is.null(targets) && is.null(objective)) {
      stop("'measures' and 'targets' must be supplied together", call. = FALSE)
    }
    if (!is.null(budget) && is.null(objective)) {
      stop(
        "'budget' requires 'objective' for joint precision assessment",
        call. = FALSE
      )
    }
    if (!is.null(budget)) {
      check_scalar(budget, "budget")
    }
    if (!is.null(min_n_stratum)) {
      check_scalar(min_n_stratum, "min_n_stratum")
    }
    if (!is.null(domains)) {
      stop(
        "use 'targets$domain' for joint constrained allocation, not 'domains'",
        call. = FALSE
      )
    }
    if (!is.null(psu)) {
      # Refused for the same reason n_alloc() refuses it. See there.
      if (!identical(fpc, "unit")) {
        stop(
          "'fpc' is not available with a PSU register: the certainty split and its size-proportional remainder need their own variance decomposition, which is not derived here. Drop 'fpc', or drop 'psu' and give 'N_psu' instead",
          call. = FALSE
        )
      }
      return(.prec_alloc_psu(
        frame = frame,
        psu = psu,
        n = n,
        measures = measures,
        targets = targets,
        objective = objective,
        budget = budget,
        unit_cost = unit_cost,
        alpha = alpha,
        deff = deff,
        resp_rate = resp_rate,
        min_n_stratum = min_n_stratum,
        df = df
      ))
    }
    return(.prec_alloc_bethel(
      frame = frame,
      n = n,
      measures = measures,
      targets = targets,
      objective = objective,
      budget = budget,
      unit_cost = unit_cost,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      min_n_stratum = min_n_stratum,
      df = df,
      fpc = fpc
    ))
  }
  if (!is.null(psu)) {
    stop(
      "'psu' requires a joint constrained allocation: supply 'measures' and 'targets'",
      call. = FALSE
    )
  }
  if (!is.null(min_n_stratum)) {
    stop(
      "'min_n_stratum' applies to joint precision assessment. Supply 'measures' and 'targets', or drop it",
      call. = FALSE
    )
  }
  check_alpha(alpha)

  prep <- .alloc_prepare_frame(
    frame,
    domains = domains,
    unit_cost = unit_cost,
    deff = deff,
    resp_rate = resp_rate,
    fpc = fpc
  )
  deff_arg <- deff
  resp_rate_arg <- resp_rate
  deff <- prep$deff_h
  resp_rate <- prep$resp_rate_h
  H <- length(prep$N_h)

  if (!is.numeric(n) || anyNA(n) || any(!is.finite(n)) || length(n) != H) {
    stop(
      "'n' must be a finite numeric vector with length nrow(frame)",
      call. = FALSE
    )
  }
  if (any(n <= 0)) {
    stop("all 'n' elements must be positive", call. = FALSE)
  }

  metrics <- .alloc_metrics(
    N_h = prep$N_h,
    S_h = prep$S_h,
    mean_h = prep$mean_h,
    n_h = n,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = prep$cost_h,
    df = df,
    N_fpc = prep$N_fpc_h
  )

  detail <- .alloc_detail(
    prep = prep,
    n_h = n,
    m_h = rep(NA_real_, H),
    M_h = prep$N_h,
    resp_rate = resp_rate,
    deff = deff,
    metrics = metrics
  )
  domain_summary <- .alloc_domain_summary(
    prep = prep,
    n_h = n,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    df = df
  )

  .new_svyplan_prec(
    se = metrics$se,
    moe = metrics$moe,
    cv = metrics$cv,
    type = "alloc",
    params = list(
      frame = frame,
      n = n,
      alpha = alpha,
      fpc = fpc,
      deff = deff_arg,
      resp_rate = resp_rate_arg,
      df = df,
      cost_h = prep$cost_h,
      domain_cols = prep$domain_cols,
      domain_idx = prep$domain_idx,
      achieved = list(n = sum(n), cv = metrics$cv, cost = metrics$cost)
    ),
    detail = detail,
    domains = domain_summary
  )
}

#' @rdname prec_alloc
#' @export
prec_alloc.svyplan_n <- function(frame, ...) {
  obj <- frame
  if (obj$type != "alloc") {
    stop("prec_alloc requires a svyplan_n of type 'alloc'", call. = FALSE)
  }
  p <- obj$params
  if (identical(obj$method, "bethel")) {
    dots <- list(...)
    dot_names <- names(dots)
    n_explicit <- length(dots) > 0L &&
      (is.null(dot_names) || any(!nzchar(dot_names)) || "n" %in% dot_names)
    if (!is.null(p$psu)) {
      # The fit's classification travels: a fresh derivation would drop the
      # PSUs the absorb loop added.
      held <- p$psu
      held$certainty <- obj$psu$certainty
      args <- list(
        frame = p$frame,
        psu = held,
        n = obj$detail$n,
        measures = p$measures,
        targets = p$targets,
        objective = p$objective,
        budget = p$budget,
        alpha = p$alpha,
        deff = p$deff,
        resp_rate = p$resp_rate,
        df = p$df,
        unit_cost = p$unit_cost,
        min_n_stratum = p$min_n_stratum,
        .allow_fractional_stages = !n_explicit
      )
      return(do.call(
        .prec_alloc_psu,
        .roundtrip_args(args, dots, .prec_alloc_psu)
      ))
    }
    args <- list(
      frame = p$frame,
      n = obj$detail$n,
      measures = p$measures,
      targets = p$targets,
      objective = p$objective,
      budget = p$budget,
      alpha = p$alpha,
      deff = p$deff,
      resp_rate = p$resp_rate,
      df = p$df,
      unit_cost = p$unit_cost,
      min_n_stratum = p$min_n_stratum,
      fpc = p$fpc %||% "unit",
      .allow_fractional_stages = !n_explicit
    )
    return(do.call(
      .prec_alloc_bethel,
      .roundtrip_args(args, dots, .prec_alloc_bethel)
    ))
  }
  n_h <- p$n_h
  if (is.null(n_h) && !is.null(obj$detail) && "n" %in% names(obj$detail)) {
    n_h <- obj$detail$n
  }
  if (is.null(n_h)) {
    stop(
      "allocation detail does not contain stratum sample sizes",
      call. = FALSE
    )
  }

  args <- list(
    frame = p$frame,
    n = n_h,
    domains = p$domain_cols,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    df = p$df
  )
  if (!.alloc_is_cluster(p$frame)) {
    args$unit_cost <- p$cost_h
  } else {
    args$fpc <- p$fpc %||% "unit"
  }
  do.call(
    prec_alloc.default,
    .roundtrip_args(args, list(...), prec_alloc.default)
  )
}

#' @keywords internal
#' @noRd
.alloc_prepare_frame <- function(
  frame,
  domains = NULL,
  unit_cost = NULL,
  deff = 1,
  resp_rate = 1,
  fpc = "unit"
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
  } else {
    stop("'frame' must contain either 'sd' or 'var'", call. = FALSE)
  }
  if (all(S_h == 0)) {
    warning(
      "all 'sd' values are zero, allocation has no variability to distribute",
      call. = FALSE
    )
  }

  mean_h <- rep(NA_real_, nrow(frame))
  if ("mean" %in% names(frame)) {
    mean_h <- frame$mean
  } else if ("p" %in% names(frame)) {
    ph <- frame$p
    if (!is.numeric(ph) || anyNA(ph) || any(!is.finite(ph))) {
      stop("'p' must contain finite numeric values", call. = FALSE)
    }
    if (any(ph < 0 | ph > 1)) {
      stop("'p' must contain values in [0, 1]", call. = FALSE)
    }
    mean_h <- ph
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
  resp_rate_h <- .alloc_resolve_h(
    resp_rate,
    frame,
    "resp_rate",
    nrow(frame),
    .check_resp_rate_h
  )

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

  prep <- list(
    frame = frame,
    N_h = as.numeric(N_h),
    S_h = as.numeric(S_h),
    mean_h = as.numeric(mean_h),
    cost_h = as.numeric(cost_h),
    deff_h = deff_h,
    resp_rate_h = resp_rate_h,
    max_weight = as.numeric(max_weight),
    take_all = as.logical(take_all),
    stratum = stratum,
    domain_cols = domain_cols,
    domain_idx = domain_idx,
    domain_values = domain_values
  )
  .alloc_cluster_prep(prep, frame, unit_cost, fpc)
}

#' Two-stage-within-strata transformation
#'
#' When the frame carries a 'icc_psu' column, the stratified two-stage
#' variance model reduces to the element model with an inflated stratum
#' SD, S_h * sqrt(k_h (1 + icc_h (n_per_psu_h - 1))), and an effective
#' per-element cost, cost_psu / n_per_psu + cost_ssu. Downstream
#' allocation, constraints, and metrics then apply unchanged.
#' @keywords internal
#' @noRd
.alloc_cluster_prep <- function(prep, frame, unit_cost, fpc = "unit") {
  if (!.alloc_is_cluster(frame)) {
    if (!identical(fpc, "unit")) {
      stop(
        "'fpc' applies to cluster allocation only, where the correction has stages to choose between. An element allocation always carries 1 - n/N. Add 'icc_psu', or drop 'fpc'",
        call. = FALSE
      )
    }
    orphan <- intersect(
      c(
        "var_ratio_psu",
        "n_per_psu",
        "cost_psu",
        "cost_ssu",
        "N_psu",
        "resp_rate_psu"
      ),
      names(frame)
    )
    if (length(orphan) > 0L) {
      stop(
        sprintf(
          "cluster column(s) %s require a 'icc_psu' column",
          paste(sQuote(orphan), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    prep$cluster <- FALSE
    return(prep)
  }

  if ("unit_cost" %in% names(frame) || !is.null(unit_cost)) {
    stop(
      "with 'icc_psu' in the frame, use 'cost_psu' and 'cost_ssu' columns instead of 'unit_cost'",
      call. = FALSE
    )
  }

  H <- nrow(frame)
  icc <- frame$icc_psu
  if (!is.numeric(icc) || anyNA(icc) || any(!is.finite(icc))) {
    stop("'icc_psu' must contain finite numeric values", call. = FALSE)
  }
  .check_cluster_icc_open(icc, context = "n_alloc()")

  var_ratio <- if ("var_ratio_psu" %in% names(frame)) {
    frame[["var_ratio_psu"]]
  } else {
    rep(1, H)
  }
  if (
    !is.numeric(var_ratio) ||
      anyNA(var_ratio) ||
      any(!is.finite(var_ratio)) ||
      any(var_ratio <= 0)
  ) {
    stop("'var_ratio_psu' must contain positive finite values", call. = FALSE)
  }

  has_costs <- any(c("cost_psu", "cost_ssu") %in% names(frame))
  if (has_costs) {
    if (!all(c("cost_psu", "cost_ssu") %in% names(frame))) {
      stop("'cost_psu' and 'cost_ssu' must be supplied together", call. = FALSE)
    }
    for (col in c("cost_psu", "cost_ssu")) {
      cc <- frame[[col]]
      if (!is.numeric(cc) || anyNA(cc) || any(!is.finite(cc)) || any(cc <= 0)) {
        stop(
          sprintf("'%s' must contain positive finite values", col),
          call. = FALSE
        )
      }
    }
  }

  n_per_psu_h <- rep(NA_real_, H)
  if ("n_per_psu" %in% names(frame)) {
    ps <- frame$n_per_psu
    if (!is.numeric(ps) || any(!is.na(ps) & (!is.finite(ps) | ps < 1))) {
      stop(
        "'n_per_psu' must contain values >= 1 (NA for cost-optimal)",
        call. = FALSE
      )
    }
    n_per_psu_h <- as.numeric(ps)
  }

  need_opt <- is.na(n_per_psu_h)
  if (any(need_opt)) {
    if (!has_costs) {
      stop(
        "'cost_psu' and 'cost_ssu' are required when 'n_per_psu' is not fixed for every stratum",
        call. = FALSE
      )
    }
    # b* = sqrt(C1 (1 - icc) / (C2 icc r)). resp_rate_psu scales cost without
    # moving the trade-off, so it is absent.
    n_per_psu_h[need_opt] <- sqrt(
      frame$cost_psu[need_opt] /
        frame$cost_ssu[need_opt] *
        (1 - icc[need_opt]) /
        (icc[need_opt] * prep$resp_rate_h[need_opt])
    )
  }

  too_big <- n_per_psu_h > prep$N_h
  if (any(too_big & !need_opt)) {
    stop(
      sprintf(
        "'n_per_psu' exceeds the stratum population for: %s",
        paste(prep$stratum[too_big & !need_opt], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(too_big & need_opt)) {
    n_per_psu_h[too_big & need_opt] <- prep$N_h[too_big & need_opt]
    warning(
      "cost-optimal 'n_per_psu' exceeds the stratum population, clamped to 'N'",
      call. = FALSE
    )
  }

  if (any(prep$take_all)) {
    stop(
      "'take_all' is not supported for cluster allocation because taking every PSU does not imply an ultimate-unit census: the within-PSU take 'n_per_psu' still applies. Drop 'take_all', or drop 'icc_psu' to allocate at the element level",
      call. = FALSE
    )
  }

  resp_rate_psu_h <- .alloc_resolve_h(
    1,
    frame,
    "resp_rate_psu",
    H,
    .check_resp_rate_h
  )

  N_psu_h <- NULL
  if ("N_psu" %in% names(frame)) {
    N_psu <- frame$N_psu
    if (
      !is.numeric(N_psu) ||
        anyNA(N_psu) ||
        any(!is.finite(N_psu)) ||
        any(N_psu < 1) ||
        any(abs(N_psu - round(N_psu)) > 1e-8)
    ) {
      stop(
        "'N_psu' must contain positive whole numbers, one per stratum",
        call. = FALSE
      )
    }
    N_psu_h <- as.numeric(round(N_psu))
    if (any(N_psu_h > prep$N_h + 1e-8)) {
      stop(
        "'N_psu' must not exceed the ultimate-unit population 'N'",
        call. = FALSE
      )
    }
  }
  prep$N_psu_h <- N_psu_h

  # The clustering penalty is set by the take that is realized, not the one
  # that is issued, so the inflation reads 'n_per_psu * resp_rate'.
  responding_take <- n_per_psu_h * prep$resp_rate_h
  prep$S_raw_h <- prep$S_h
  if (identical(fpc, "stage") && is.null(N_psu_h)) {
    stop(
      "fpc = \"stage\" needs the PSU population: add an 'N_psu' column, or use the default fpc = \"unit\"",
      call. = FALSE
    )
  }
  # Two components, since each stage's correction multiplies its own.
  within_frac <- if (identical(fpc, "stage")) {
    pmin(1, responding_take / (prep$N_h / N_psu_h))
  } else {
    rep(0, length(responding_take))
  }
  inflation <- var_ratio *
    (icc * responding_take + (1 - within_frac) * (1 - icc))
  prep$S_h <- prep$S_h * sqrt(inflation)
  # A constant in n_h, so it travels as an effective population rather than
  # as a separate stage.
  prep$N_fpc_h <- switch(
    fpc,
    unit = prep$N_h,
    none = rep(Inf, length(prep$N_h)),
    stage = ifelse(
      icc > 0,
      N_psu_h * inflation / (var_ratio * icc),
      Inf
    )
  )
  prep$fpc <- fpc
  prep$unit_resp_rate_h <- prep$resp_rate_h
  prep$resp_rate_h <- prep$resp_rate_h * resp_rate_psu_h
  prep$responding_take_h <- responding_take
  if (has_costs) {
    prep$cost_h <- frame$cost_psu / n_per_psu_h + frame$cost_ssu
    prep$cost_psu_h <- frame$cost_psu
    prep$cost_ssu_h <- frame$cost_ssu
  }
  prep$cluster <- TRUE
  prep$has_stage_costs <- has_costs
  prep$n_per_psu_fixed_h <- !need_opt
  prep$n_per_psu_h <- n_per_psu_h
  prep$icc_psu_h <- icc
  prep$var_ratio_psu_h <- var_ratio
  prep
}

#' @keywords internal
#' @noRd
.alloc_is_cluster <- function(frame) {
  "icc_psu" %in% names(frame)
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
  if (is.null(prep$M_src) || !any(prep$M_src == "N_psu")) {
    return("")
  }
  sprintf(
    ". The ceiling comes from the PSU universe in %s rather than from the population: with 'n_per_psu' fixed, only 'N_psu * n_per_psu' ultimate units are available there. That is what this path's with-replacement first-stage approximation and the supplied PSU availability allow, not a statement that the target is unattainable under the real design",
    paste(sQuote(prep$stratum[prep$M_src == "N_psu"]), collapse = ", ")
  )
}

.alloc_bounds <- function(
  N_h,
  max_weight,
  take_all,
  min_n_stratum = NULL,
  psu_cap = NULL
) {
  H <- length(N_h)
  m_h <- pmin(rep(1, H), N_h)
  M_h <- as.numeric(N_h)
  # Which constraint produced each bound, so that .binding says not only that
  # the allocation sits on a bound but which one put it there.
  m_src <- rep(NA_character_, H)
  M_src <- rep("N", H)

  if (!is.null(min_n_stratum)) {
    hit <- min_n_stratum > m_h
    m_h <- pmax(m_h, min_n_stratum)
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
    M_src[psu_cap < M_h] <- "N_psu"
    M_h <- pmin(M_h, psu_cap)
  }

  if (any(take_all)) {
    m_h[take_all] <- N_h[take_all]
    M_h[take_all] <- N_h[take_all]
    m_src[take_all] <- "take_all"
    M_src[take_all] <- "take_all"
  }

  infeasible <- m_h > M_h + 1e-8
  if (any(infeasible)) {
    if (any(infeasible & M_src == "N_psu")) {
      stop(
        "constraints are infeasible: the PSU universe 'N_psu' caps the stratum below its own lower bound. Raise 'N_psu' or 'n_per_psu', or relax 'min_n_stratum' or 'max_weight'",
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
  N_fpc = NULL
) {
  n_net <- n_h * resp_rate
  n_eff <- n_net / deff
  W_h <- N_h / sum(N_h)

  term <- numeric(length(n_h))
  good <- n_eff > 0
  # N_fpc is the population the correction divides by, which is the stratum
  # population unless a cluster allocation has folded a first-stage
  # correction into it. Weights and bounds always read the real N_h.
  fpc <- pmax(0, 1 - n_net / (N_fpc %||% N_h))
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
  if (is.null(lo_i)) {
    lo_i <- as.integer(ceiling(pmax(ifelse(is.na(m_h), 0, m_h), 0) - 1e-9))
  }
  if (is.null(hi_i)) {
    hi_i <- as.integer(floor(M_h + 1e-9))
  }
  n_int <- switch(
    mode,
    cv = pmin(pmax(as.integer(ceiling(n_h - 1e-9)), lo_i), hi_i),
    budget = .round_within_budget(
      n_h,
      prep$cost_h,
      budget,
      lo_i,
      hi_i,
      prep$N_h,
      prep$S_h
    ),
    .round_oric_bounded(n_h, lo_i, hi_i)
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
    out$n_per_psu <- prep$n_per_psu_h
    out$n_psu <- n_h / prep$n_per_psu_h
    n_psu_int <- ceiling(out$n_psu)
    if (!is.null(prep$N_psu_h)) {
      # The element bound is expressed in ultimate units, so rounding the PSU
      # count up can still land past the universe when the take is fractional.
      n_psu_int <- pmin(n_psu_int, prep$N_psu_h)
      out$N_psu <- prep$N_psu_h
      out$.psu_frac <- out$n_psu / prep$N_psu_h
    }
    out$n_psu_int <- as.integer(n_psu_int)
  }
  if (any(prep$take_all)) {
    out$take_all <- prep$take_all
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
#' is largest while the budget allows.
#' @keywords internal
#' @noRd
.round_within_budget <- function(n_h, cost_h, budget, lower, upper, N_h, S_h) {
  n_int <- pmin(pmax(as.integer(floor(n_h + 1e-9)), lower), upper)
  W2S2 <- (N_h / sum(N_h))^2 * S_h^2
  repeat {
    spare <- budget - sum(n_int * cost_h)
    can <- which(n_int < upper & cost_h <= spare + 1e-9)
    if (length(can) == 0L) {
      break
    }
    gain <- W2S2[can] * (1 / n_int[can] - 1 / (n_int[can] + 1)) / cost_h[can]
    j <- can[which.max(gain)]
    n_int[j] <- n_int[j] + 1L
  }
  as.integer(n_int)
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
    cv = m$cv
  )
}

#' Operational (integer) design for a stratified two-stage allocation
#'
#' Chooses a whole per-stratum take b_h (floor/ceiling candidate with the
#' lowest variance-cost product when stage costs are known) and a whole
#' PSU count a_h per stratum. In cv mode a_h matches or beats each
#' stratum's continuous variance contribution. In budget mode PSUs are
#' removed/added greedily so the field cost sum(a_h * (cost_psu +
#' cost_ssu * b_h)) never exceeds the budget. In n mode a_h approximates
#' the continuous element total. Returns the design plus recomputed
#' metrics and updates the detail integers.
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
  df = NULL
) {
  H <- length(n_h)
  ps <- prep$n_per_psu_h
  # Ultimate-unit rate, not the combined one: PSU loss stays out of the bracket.
  resp_h <- rep_len(prep$unit_resp_rate_h %||% resp_rate, H)
  icc <- prep$icc_psu_h
  var_ratio <- prep$var_ratio_psu_h
  S_raw <- prep$S_raw_h
  W <- prep$N_h / sum(prep$N_h)
  lo_i <- as.integer(ceiling(detail$.lower - 1e-9))
  hi_i <- as.integer(floor(detail$.upper + 1e-9))
  b_h <- vapply(
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
      cand <- cand[cand >= 1L & cand <= hi_i[h]]
      feasible <- vapply(
        cand,
        function(b) {
          ceiling(lo_i[h] / b) <= floor(hi_i[h] / b)
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
      if (isTRUE(prep$has_stage_costs)) {
        score <- (var_ratio[h] * (1 + icc[h] * (cand * resp_h[h] - 1))) *
          (prep$cost_psu_h[h] / cand + prep$cost_ssu_h[h])
        cand[which.min(score)]
      } else {
        cand[which.min(abs(cand - ps[h]))]
      }
    },
    integer(1L)
  )

  e_target <- NULL
  if (mode == "n") {
    # A fixed n_per_psu keeps its take and moves its element count instead.
    b_free <- b_h
    e_target <- .round_oric_bounded(n_h, lo_i, hi_i)
    b_h <- vapply(
      seq_len(H),
      function(h) {
        if (isTRUE(prep$n_per_psu_fixed_h[h])) {
          return(b_free[h])
        }
        e <- e_target[h]
        root <- seq_len(max(1L, as.integer(floor(sqrt(e)))))
        small <- root[e %% root == 0L]
        divisors <- sort(unique(c(small, e %/% small)))
        if (isTRUE(prep$has_stage_costs)) {
          score <- (var_ratio[h] * (1 + icc[h] * (divisors * resp_h[h] - 1))) *
            (prep$cost_psu_h[h] / divisors + prep$cost_ssu_h[h])
          return(divisors[which.min(score)])
        }
        divisors[which.min(abs(divisors - ps[h]))]
      },
      integer(1L)
    )
  }

  a_min <- as.integer(ceiling(lo_i / b_h))
  a_max <- as.integer(floor(hi_i / b_h))
  # hi_i counts ultimate units, and the integer take b_h chosen here need not
  # be the one that set it, so the PSU universe has to bound a_max directly.
  if (!is.null(prep$N_psu_h)) {
    a_max <- pmin(a_max, as.integer(prep$N_psu_h))
  }
  if (any(a_min > a_max)) {
    stop(
      "no whole-cluster design satisfies the integer allocation bounds",
      call. = FALSE
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
    pmin(1, op_take / (prep$N_h / prep$N_psu_h))
  } else {
    rep(0, H)
  }
  op_inflation <- var_ratio *
    (icc * op_take + (1 - op_within_frac) * (1 - icc))
  S_op <- S_raw * sqrt(op_inflation)
  N_fpc_op <- switch(
    prep$fpc %||% "unit",
    unit = prep$N_h,
    none = rep(Inf, H),
    stage = ifelse(
      icc > 0,
      prep$N_psu_h * op_inflation / (var_ratio * icc),
      Inf
    )
  )
  Cj <- W^2 * S_op^2 * deff / (b_h * resp_rate)

  if (mode == "cv") {
    a_h <- as.integer(ceiling(
      n_h *
        (1 + icc * (b_h * resp_h - 1)) /
        (1 + icc * (ps * resp_h - 1)) /
        b_h -
        1e-9
    ))
    a_h <- pmin(pmax(a_h, a_min), a_max)
  } else if (mode == "budget") {
    min_cost <- sum(a_min * psu_cost)
    if (budget < min_cost - 1e-8) {
      stop(
        sprintf(
          "'budget' cannot fund one whole PSU per stratum while satisfying the lower bounds (minimum field cost = %.4g)",
          min_cost
        ),
        call. = FALSE
      )
    }
    a_h <- pmin(pmax(as.integer(floor(n_h / b_h + 1e-9)), a_min), a_max)
    repeat {
      if (sum(a_h * psu_cost) - budget <= 1e-8) {
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
      spare <- budget - sum(a_h * psu_cost)
      cand <- which(psu_cost <= spare + 1e-8 & a_h < a_max)
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
  }

  e_h <- a_h * b_h
  m <- .alloc_metrics(
    N_h = prep$N_h,
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
  detail$n_psu_int <- a_h
  detail$n_per_psu_int <- b_h

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
  if (!is.null(target_cv)) {
    op_cv <- if (length(prep$domain_idx) == 0L) {
      m$cv
    } else {
      prep_op <- prep
      prep_op$S_h <- S_op
      prep_op$N_fpc_h <- N_fpc_op
      .alloc_domain_cv_max(prep_op, e_h, alpha, deff, resp_rate)
    }
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
  list(
    operational = list(
      n = sum(e_h),
      cost = if (isTRUE(prep$has_stage_costs)) {
        sum(a_h * psu_cost)
      } else {
        NA_real_
      },
      se = m$se,
      moe = m$moe,
      cv = op_cv
    ),
    detail = detail
  )
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
    function(idx) {
      .alloc_metrics(
        N_h = prep$N_h[idx],
        S_h = prep$S_h[idx],
        mean_h = prep$mean_h[idx],
        n_h = n_h[idx],
        alpha = alpha,
        deff = .subset_h(deff, idx),
        resp_rate = .subset_h(resp_rate, idx),
        cost_h = prep$cost_h[idx],
        N_fpc = prep$N_fpc_h[idx]
      )$cv
    },
    numeric(1)
  )
  max(cvs)
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
