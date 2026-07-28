#' Constrained Stratified Allocation
#'
#' Distribute a total sample size across strata defined by a single
#' stratification variable, under a fixed total \eqn{n}, target CV, or budget.
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
#'   for example the pool summary of an executed samplyr sample
#'   (`samplyr::frame_summary()`), once the measure columns are added.
#'
#'   When a design stratifies by several variables at once (e.g. region
#'   \eqn{\times} urbanicity), cross them into a single variable before
#'   calling `n_alloc` (e.g. with [interaction()]) so that each row maps
#'   to exactly one population cell.
#'
#'   When `measures` and `targets` are supplied, `frame` instead describes
#'   atomic allocation strata for joint constrained allocation. It requires
#'   unique, non-missing `stratum` and positive finite `N` columns. It may
#'   contain `unit_cost`, `max_weight`, `take_all`, and complete domain
#'   classification columns named by `targets$domain`; `sd`, `var`, `mean`,
#'   and `p` then belong in `measures`. Each row must be homogeneous for every
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
#'           (\eqn{\sqrt{\text{variance}}}), or
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
#'       Use `NA` for strata without a cap.}
#'     \item{`take_all`}{Logical (or 0/1). If `TRUE`, every unit in the
#'       stratum is included, a census stratum. Useful for small strata
#'       whose total population is tiny enough to enumerate.}
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
#'   the scalar arguments; `NA` uses the scalar default. In fixed-take
#'   multistage mode, `icc_psu` is required and `var_ratio_psu` defaults to 1;
#'   three-stage mode also requires `icc_ssu` and derives `var_ratio_ssu` as
#'   `var_ratio_psu * (1 - icc_psu)` when it is absent, which is the value the
#'   variance decomposition implies.
#'   Stage parameters may instead be supplied as stratum defaults in `frame`;
#'   non-missing measure-row values take precedence. Only rows selected by a
#'   target and its domain are value-validated. Unused rows are ignored after
#'   their non-empty, unique keys and frame-stratum membership are checked.
#'   Must be supplied with `targets`.
#' @param targets Optional long data frame with one precision requirement per
#'   row. It requires indicator `name` and exactly one positive `cv` or `moe`.
#'   Optional `domain` and `level` select a frame domain (`.overall` and
#'   `level = NA` denote the whole population), `alpha` overrides the scalar
#'   value for MOE, and `constraint` supplies a stable unique identifier.
#'   Must be supplied with `measures`. Optional when `objective` and `budget`
#'   are given, which requests the best design a budget can buy with no hard
#'   precision requirement.
#' @param objective Optional estimates whose priority-weighted relative
#'   variance is minimized among the allocations that meet `targets` and the
#'   budget. Either a character vector of indicator names (overall domain,
#'   equal priority) or a data frame with `name` plus optional `domain`,
#'   `level`, non-negative `priority`, and `component` columns, using the same
#'   indicator-domain identifiers as `targets`. At least one priority must be
#'   positive; priorities need not be normalized. Requires `budget` and must be
#'   supplied with `measures`.
#' @param alloc Allocation rule: `"neyman"` (default), `"optimal"`,
#'   `"proportional"`, or `"power"`.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of
#'   per-stratum unit costs, overriding `frame$unit_cost`. Not used with
#'   fixed-take multistage joint allocation, which requires stage costs.
#' @param alpha Significance level, default 0.05.
#' @param deff Design effect multiplier (> 0).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1.
#' @param min_n_stratum Optional minimum sample size per stratum.
#' @param alloc_q Bankier power parameter from 0 to 1, used when
#'   `alloc = "power"`.
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
#'       moves to the nearest whole multiple of it instead; in
#'       `cv` mode each stratum is rounded up so the integer design
#'       meets the target; in `budget` mode units are added by variance
#'       reduction per unit cost so the integer design stays within
#'       budget. Always inside the integerized bounds
#'       (`ceiling(.lower)`, `floor(.upper)`); an error is raised when
#'       no integer allocation can satisfy them.}
#'     \item{`weight`}{Design weight `N / n`.}
#'     \item{`n_eff`}{Effective sample size `n * resp_rate / deff`.}
#'     \item{`.lower`, `.upper`}{Bounds applied to the stratum
#'       (from `min_n_stratum`, `max_weight`, `take_all`, or `N`).}
#'     \item{`.binding`}{Whether the allocation sits on one of its
#'       bounds.}
#'     \item{`take_all`, `mean`}{Present when take-all strata or stratum
#'       means were supplied.}
#'     \item{`n_per_psu`, `n_psu`}{Cluster mode only: the continuous
#'       per-stratum take and implied number of PSUs (`n / n_per_psu`).}
#'     \item{`n_psu_int`, `n_per_psu_int`}{Cluster mode only: the
#'       whole-unit field design (whole PSUs and whole takes), chosen so
#'       that `budget` designs stay within budget and `cv` designs meet
#'       the target; `n_int = n_psu_int * n_per_psu_int`.}
#'     \item{`n_psu`, `n_psu_int`, `n_per_psu`, `n_per_ssu`}{Joint fixed-take
#'       multistage mode: continuous and whole first-stage decisions and the
#'       fixed later-stage takes. `n_per_ssu` is present only for three stages.
#'       Public `n` remains in ultimate-unit units, so
#'       `n = n_psu * n_per_psu` for two stages and
#'       `n = n_psu * n_per_psu * n_per_ssu` for three stages; the same exact
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
#' that meets all indicator-domain precision requirements simultaneously.
#' This mode is mutually exclusive with scalar `n`, `cv`, `budget`, `domains`,
#' and an explicit `alloc` rule: the target table defines the constraints and
#' the allocation is solved jointly.
#'
#' CV targets apply to domain means or totals (the CV is the same on either
#' scale). MOE targets are absolute margins of error for domain means. The
#' variance model uses Wald/linearized variances, the package's gross-sample
#' response-rate convention, design effects, and finite population correction.
#' Domain classifications may overlap; for example, separate `region` and
#' `residence` targets can use the same atomic `region x residence` frame rows.
#'
#' Fixed-take two- and three-stage joint designs keep only the first-stage PSU
#' count free. Later-stage takes must be positive whole numbers. The canonical
#' optimizer works in PSU units, while public `n`, `min_n_stratum`, `max_weight`, and
#' precision assessments remain in ultimate-unit units. PSU upper bounds use
#' `min(N_psu, N / n_per_psu)` for two stages and
#' `min(N_psu, N_ssu / n_per_psu,
#' N / (n_per_psu * n_per_ssu))` for three stages. Because sampling
#' every PSU does not establish an ultimate-unit census, `take_all` is not
#' supported in this mode. The response rate applies once to the gross
#' ultimate-unit take. Jointly optimized later-stage takes are not supported.
#'
#' What the optimum then means is worth stating plainly. The convex problem
#' is solved to optimality for the planning coefficients supplied, and the
#' feasibility and KKT diagnostics report on that problem. It is a planning
#' model rather than an exact design-based variance for every multistage
#' design: the later-stage takes are fixed inputs rather than decision
#' variables, and the stage variances enter through the design-effect and
#' FPC adapter described above. A design that is optimal here is optimal
#' within that family, not across all multistage designs.
#'
#' The three population columns are all **stratum totals**, not per-parent
#' counts: `N` is the number of ultimate units in the stratum, `N_psu` the
#' number of PSUs in it, and `N_ssu` the number of SSUs across *all* of its
#' PSUs. The take columns are the mirror image, counting per parent unit:
#' `n_per_psu` units (or SSUs) inside one selected PSU and `n_per_ssu` units
#' inside one selected SSU. The `N_ssu / n_per_psu` bound above reads
#' correctly only under this convention, since it divides a stratum total by
#' a per-PSU take to bound the number of PSUs. A frame that stores SSUs per
#' PSU instead must be converted first, typically
#' `N_ssu = N_psu * ssu_per_psu`.
#'
#' In two-stage mode `n_per_psu` is the number of ultimate units sampled in a
#' selected PSU, `cost_psu` is the cost per selected PSU, and `cost_ssu` is the
#' cost per sampled ultimate unit. The cost of one selected PSU is therefore
#' `cost_psu + cost_ssu * n_per_psu`. In three-stage mode `n_per_psu` is the
#' number of SSUs sampled per selected PSU, `n_per_ssu` is the number of
#' ultimate units sampled per selected SSU, `cost_ssu` is per sampled SSU,
#' and `cost_tsu` is per sampled ultimate unit. The corresponding cost is
#' `cost_psu + cost_ssu * n_per_psu +
#' cost_tsu * n_per_psu * n_per_ssu`.
#'
#' For indicator \eqn{k} and stratum \eqn{h}, write \eqn{m_h} for
#' `n_per_psu`, \eqn{q_h} for `n_per_ssu`, and \eqn{S^2_{hk}} for `var` (or
#' `sd^2`, or `p * (1 - p)`). The fixed variance multiplier is
#' \eqn{D_{hk}=k_{1,hk}(1+\delta_{1,hk}(m_h-1))} at two stages and
#' \eqn{D_{hk}=k_{1,hk}\delta_{1,hk}m_hq_h+
#' k_{2,hk}(1+\delta_{2,hk}(q_h-1))} at three stages. At three stages
#' `var_ratio_ssu` is not a free parameter: `var_ratio_psu` rescales the components' unit
#' variance to the analysis variable and `var_ratio_ssu` does the same for the
#' within-PSU part, so
#' \eqn{k_{2}=k_{1}(1-\delta_{1})}. That identity is what makes
#' \eqn{D_{hk}} collapse to \eqn{k_{1}} at \eqn{m_h=q_h=1}, where no
#' clustering is left to inflate anything. Leaving `var_ratio_ssu` out of `measures`
#' applies it; supplying a value overrides it, which is only meaningful when
#' the two stages' ratios come from different decompositions. With ultimate-unit take
#' \eqn{t_h=m_h} or \eqn{m_hq_h}, first-stage decision \eqn{a_h}, response
#' rate \eqn{r_{hk}}, and extra design effect \eqn{d_{hk}}, the implemented
#' total-variance contribution is
#' \deqn{N_h^2 S_{hk}^2 D_{hk}d_{hk}/(r_{hk}a_ht_h)
#'       - N_h S_{hk}^2 D_{hk}d_{hk}.}
#' Thus response inflation is applied once. Leave `deff = 1` unless it
#' represents a source not already captured by the stage deltas and var_ratio factors.
#' Deltas may include 0 and 1; var_ratio values must be positive finite. All stage
#' populations and fixed takes must be positive whole numbers.
#'
#' The continuous result is a KKT-certified global optimum of the convex,
#' fixed-coefficient problem. `$detail$n_int` is a deterministic, feasible,
#' locally cleaned operational recommendation; it is not claimed to be the
#' globally optimal integer allocation. Full achieved precision is in
#' `$constraints` and `$operational$constraints`. Use [prec_alloc()] to assess
#' the fitted, operational, or a modified allocation.
#'
#' A constraint table records `.metric`, `.target`, `.achieved`, `.ratio`, and
#' `.residual = .ratio - 1`; `.pass` means `.ratio <= 1 + .tolerance`.
#' `.binding` identifies a numerically active target. `.multiplier` is the
#' continuous Lagrange multiplier and `.sensitivity` is the local derivative
#' of minimum variable cost with respect to the target. Both are `NA` for a
#' pure precision assessment; multipliers are also `NA` when duplicate
#' normalized constraints make an individual split unidentified. The stored
#' `feasibility_tolerance` and per-row `.tolerance` are the public numerical
#' acceptance contract.
#'
#' ## Fixed budget and a weighted objective
#'
#' Joint allocation answers two different planning questions, distinguished by
#' whether `objective` is supplied. Without it, `n_alloc()` returns the
#' cheapest design meeting every requirement in `targets`. With `objective`
#' and `budget`, it returns the best design that budget can buy: among the
#' allocations that meet every hard target and cost no more than `budget`, the
#' one minimizing
#' \deqn{Q(n) = \sum_j w_j \, \mathrm{cv}_j^2(n),}
#' the priority-weighted sum of the objective components' relative variances.
#' The three roles stay distinct: `targets` are pass/fail requirements,
#' `objective` expresses preference among the allocations that pass, and
#' `budget` is a hard ceiling on variable cost.
#'
#' | Inputs | Meaning |
#' | --- | --- |
#' | `measures` + `targets` | cheapest design meeting every target |
#' | `measures` + `objective` + `budget` | best design the budget can buy |
#' | `measures` + `targets` + `objective` + `budget` | best design the budget can buy, subject to the targets |
#' | `measures` + `objective` without `budget` | error |
#' | `measures` + `targets` + `budget` without `objective` | error |
#'
#' Two consequences of the objective being in squared-CV units are worth
#' stating. Objective rows are always relative-variance, even when the same
#' indicator carries an MOE hard target, so they take no `cv` or `moe` column
#' and they need a non-negligible domain total. And `priority = c(2, 1)`
#' weights *variances* 2 to 1, not CVs. Scaling every priority by a constant
#' leaves the allocation unchanged, and a zero-priority component does not
#' affect it at all.
#'
#' The solver appends the objective as one more reciprocal constraint column
#' and root searches the objective bound whose minimum cost equals the budget.
#' The model is convex with fixed coefficients, so this recovers the global
#' continuous optimum with the same KKT certification as the minimum-cost mode.
#' Solving the minimum-cost problem with the added requirement
#' \eqn{Q \le Q^{*}} returns the same allocation at a cost equal to the budget.
#'
#' Infeasibility takes three distinct forms, each with its own message: the
#' targets may be unattainable even at the stratum upper bounds; they may be
#' attainable but unaffordable, in which case the cheapest target-feasible cost
#' and the shortfall are reported; or the continuous problem may be feasible
#' while no whole-unit allocation meets every target inside the budget.
#'
#' `ceiling()` is not a safe integer start under a budget, so the operational
#' design is built from several roundings of the continuous optimum, repaired
#' and trimmed to be both target-feasible and affordable, then improved by
#' spending residual budget and by pairwise exchanges. It satisfies
#' `$operational$cost <= budget` and `all($operational$constraints$.pass)`, and
#' it is a feasible, locally improved recommendation rather than a globally
#' optimal integer allocation. The greedy construction is incomplete: it can
#' report failure on a problem that does have a feasible integer point.
#'
#' [predict()] on a fixed-budget result varies `budget` in `newdata` and
#' returns the cost-versus-objective frontier, reusing the same solves.
#'
#' For fieldwork, use `$detail$n_int`, not the analytical continuous `n`.
#' Fixed-take count columns are numeric vectors whose values are exact whole
#' numbers; the `_int` suffix denotes operational integer semantics rather
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
#' \eqn{\times} urbanicity), cross them into one variable first:
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
#' SD inflated to `sd * sqrt(var_ratio_psu * (1 + icc_psu * (n_per_psu - 1)))`
#' and, when stage costs are given, a per-element cost of
#' `cost_psu / n_per_psu + cost_ssu`. All solve modes, allocation
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
#'   take is what `n` mode gives up its exact total for: the stratum is
#'   fielded as a whole number of clusters of that size, so its element
#'   count lands on the nearest multiple.
#' - `cost_psu`, `cost_ssu` (together): per-PSU and per-element costs.
#'   These are required for `budget` mode or when `n_per_psu` is not fixed. They
#'   replace `unit_cost`, which is not allowed in this mode.
#'
#' The finite population correction stays at the element level, an
#' approximation consistent with [n_cluster()]'s variance model.
#'
#' Because `icc_psu` already accounts for the clustering, leave
#' `deff` at 1 unless it captures a *different* source of design
#' effect (e.g. weighting loss). A clustering `deff` on top of
#' `icc_psu` would double-count. The constraints `min_n_stratum`,
#' `max_weight`, and `take_all` stay in element units. For fielding,
#' use the whole-unit design in `n_psu_int` and `n_per_psu_int`
#' (`n_int = n_psu_int * n_per_psu_int`). Its actual field cost and
#' precision are reported in `$operational` and, in `budget` mode,
#' never exceed the budget.
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
#' - **proportional**: \eqn{n_h \propto N_h}
#' - **neyman**: \eqn{n_h \propto N_h S_h}
#' - **optimal**: \eqn{n_h \propto N_h S_h / \sqrt{c_h}}
#' - **power**: Bankier (1988), \eqn{n_h \propto S_h N_h^{q}}{n_h ~ S_h * N_h^alloc_q},
#'   with exponent `alloc_q`
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
#' @seealso [prec_alloc()] for the inverse, [strata_bound()] for
#'   constructing the strata to allocate over, and [n_multi()] for the
#'   unstratified multi-indicator size that joint allocation refines.
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
  objective = NULL,
  alloc = c("neyman", "optimal", "proportional", "power"),
  unit_cost = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  min_n_stratum = NULL,
  alloc_q = 0.5,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, n_alloc.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(n_alloc.default, c(.plan, list(...))))
  .check_unused_dots(...)
  alloc_default <- c("neyman", "optimal", "proportional", "power")
  alloc_explicit <- !missing(alloc) && !identical(alloc, alloc_default)
  alloc_q_explicit <- !missing(alloc_q) && !identical(alloc_q, 0.5)
  joint_any <- !is.null(measures) || !is.null(targets) || !is.null(objective)
  if (joint_any) {
    if (is.null(measures)) {
      stop("joint constrained allocation requires 'measures'", call. = FALSE)
    }
    if (!is.null(n) || !is.null(cv)) {
      stop("joint constrained allocation cannot be combined with 'n' or 'cv'",
           call. = FALSE)
    }
    if (is.null(objective)) {
      if (is.null(targets)) {
        stop("'measures' and 'targets' must be supplied together",
             call. = FALSE)
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
      stop("use 'targets$domain' for joint constrained allocation, not 'domains'",
           call. = FALSE)
    }
    if (alloc_explicit) {
      stop("'alloc' is not used for joint constrained allocation; allocation is determined by 'targets'",
           call. = FALSE)
    }
    if (alloc_q_explicit) {
      stop("'alloc_q' is not used for joint constrained allocation",
           call. = FALSE)
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
      budget = budget
    ))
  }
  alloc <- match.arg(alloc)
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)

  if (!is.null(min_n_stratum)) check_scalar(min_n_stratum, "min_n_stratum")
  if (alloc == "power") {
    if (!is.numeric(alloc_q) || length(alloc_q) != 1L || is.na(alloc_q) ||
        alloc_q < 0 || alloc_q > 1) {
      stop("'alloc_q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }

  mode_count <- (!is.null(n)) + (!is.null(cv)) + (!is.null(budget))
  if (mode_count != 1L) {
    stop("specify exactly one of 'n', 'cv', or 'budget'", call. = FALSE)
  }

  prep <- .alloc_prepare_frame(frame, domains = domains, unit_cost = unit_cost)
  N_h <- prep$N_h
  S_h <- prep$S_h
  mean_h <- prep$mean_h
  cost_h <- prep$cost_h

  bounds <- .alloc_bounds(
    N_h = N_h,
    max_weight = prep$max_weight,
    take_all = prep$take_all,
    min_n_stratum = min_n_stratum
  )
  m_h <- bounds$m_h
  M_h <- bounds$M_h
  lo_i <- as.integer(ceiling(m_h - 1e-9))
  hi_i <- as.integer(floor(M_h + 1e-9))
  if (any(lo_i > hi_i)) {
    bad <- prep$stratum[lo_i > hi_i]
    stop(
      sprintf("no integer sample size satisfies the bounds for stratum: %s",
              paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }

  lo <- sum(m_h)
  hi <- sum(M_h)
  tol <- 1e-8

  a_h <- .alloc_weights(alloc, alloc_q, N_h, S_h, cost_h)
  if (!is.finite(sum(a_h)) || sum(a_h) <= 0) a_h <- N_h

  mode <- if (!is.null(n)) "n" else if (!is.null(cv)) "cv" else "budget"
  if (mode == "budget" && isTRUE(prep$cluster) &&
      !isTRUE(prep$has_stage_costs)) {
    stop(
      "'budget' mode with 'icc_psu' requires 'cost_psu' and 'cost_ssu' columns",
      call. = FALSE
    )
  }
  target_total <- NA_real_

  if (mode == "n") {
    check_scalar(n, "n")
    if (n < lo - tol) {
      stop("'n' is below the minimum feasible total under constraints",
           call. = FALSE)
    }
    if (n > hi + tol) {
      stop("'n' exceeds the maximum feasible total (census bound)",
           call. = FALSE)
    }
    if (round(n) < sum(lo_i) || round(n) > sum(hi_i)) {
      stop(
        sprintf(
          "no integer allocation reaches total %d within the integer bounds (feasible totals: %d to %d)",
          as.integer(round(n)), sum(lo_i), sum(hi_i)
        ),
        call. = FALSE
      )
    }
    target_total <- n
  } else if (mode == "cv") {
    check_scalar(cv, "cv")
    if (anyNA(mean_h)) {
      stop("'mean' (or 'p') is required in frame when solving for 'cv'",
           call. = FALSE)
    }

    cv_for_total <- function(n_total) {
      n_h <- .rna_alloc(a_h, n_total, m_h, M_h)
      if (length(prep$domain_idx) == 0L) {
        .alloc_metrics(
          N_h = N_h, S_h = S_h, mean_h = mean_h, n_h = n_h,
          alpha = alpha, deff = deff, resp_rate = resp_rate,
          cost_h = cost_h
        )$cv
      } else {
        .alloc_domain_cv_max(
          prep = prep, n_h = n_h,
          alpha = alpha, deff = deff, resp_rate = resp_rate
        )
      }
    }

    cv_lo <- cv_for_total(lo)
    if (cv <= 0) stop("'cv' must be positive", call. = FALSE)
    if (cv_lo <= cv + tol) {
      target_total <- lo
    } else {
      cv_hi <- cv_for_total(hi)
      if (!is.finite(cv_hi) || cv_hi > cv + tol) {
        stop("target 'cv' is unattainable under current constraints",
             call. = FALSE)
      }
      target_total <- uniroot(
        function(x) cv_for_total(x) - cv,
        interval = c(lo, hi), tol = 1e-8
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
      stop("'budget' is below the minimum feasible cost under constraints",
           call. = FALSE)
    }
    if (budget > cost_hi + tol) {
      stop("'budget' exceeds the maximum feasible cost (census bound)",
           call. = FALSE)
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
        interval = c(lo, hi), tol = 1e-8
      )$root
    }
  }

  n_h <- .rna_alloc(a_h, target_total, m_h, M_h)
  metrics <- .alloc_metrics(
    N_h = N_h, S_h = S_h, mean_h = mean_h, n_h = n_h,
    alpha = alpha, deff = deff, resp_rate = resp_rate, cost_h = cost_h
  )

  detail <- .alloc_detail(
    prep = prep, n_h = n_h, m_h = m_h, M_h = M_h,
    resp_rate = resp_rate, deff = deff,
    mode = mode, budget = budget, lo_i = lo_i, hi_i = hi_i
  )
  domain_summary <- .alloc_domain_summary(
    prep = prep, n_h = n_h,
    alpha = alpha, deff = deff, resp_rate = resp_rate
  )

  params <- list(
    frame = frame,
    alloc = alloc,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    cost_h = cost_h,
    min_n_stratum = min_n_stratum,
    domain_cols = prep$domain_cols,
    mode = mode,
    alloc_q = if (alloc == "power") alloc_q else NULL,
    n_h = n_h,
    achieved = list(n = sum(n_h), cv = metrics$cv, cost = metrics$cost)
  )
  params[[mode]] <- switch(mode, n = n, cv = cv, budget = budget)

  if (isTRUE(prep$cluster)) {
    opc <- .alloc_operational_cluster(
      prep, detail, n_h, mode, budget,
      alpha = alpha, deff = deff, resp_rate = resp_rate,
      target_cv = if (mode == "cv") cv else NULL
    )
    detail <- opc$detail
    operational <- opc$operational
  } else {
    operational <- .alloc_operational_element(
      prep, detail, alpha = alpha, deff = deff, resp_rate = resp_rate
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
      stop("generalized precision results cannot be combined with scalar 'n', 'cv', or 'budget'",
           call. = FALSE)
    }
    targets <- p$targets
    if (is.null(p$objective)) {
      # Minimum-cost mode inverts by pinning the achieved precision as the
      # requirement.
      if (!"cv" %in% names(targets)) targets$cv <- NA_real_
      if (!"moe" %in% names(targets)) targets$moe <- NA_real_
      targets$cv[] <- NA_real_
      targets$moe[] <- NA_real_
      matched <- match(targets$constraint, x$detail$constraint)
      if (anyNA(matched)) {
        stop("precision result is missing stored constraint identifiers",
             call. = FALSE)
      }
      is_cv <- x$detail$.metric[matched] == "cv"
      targets$cv[is_cv] <- x$detail$.achieved[matched][is_cv]
      targets$moe[!is_cv] <- x$detail$.achieved[matched][!is_cv]
      budget_arg <- NULL
    } else {
      # Budget-objective mode inverts through cost, not precision: pinning
      # achieved precision as hard targets would over-constrain a design that
      # already spends its whole budget. The targets stay as specified and the
      # assessed allocation's cost becomes the budget.
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
      min_n_stratum = p$min_n_stratum
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
    n = n, cv = cv, budget = budget,
    alloc = p$alloc %||% "neyman",
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    min_n_stratum = p$min_n_stratum,
    alloc_q = p$alloc_q %||% 0.5
  )
  if (!.alloc_is_cluster(p$frame)) {
    args$unit_cost <- p$cost_h
  }
  do.call(n_alloc.default, .roundtrip_args(args, list(...), n_alloc.default))
}

#' Precision for a Constrained Allocation
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
#'   allocation, omission uses its continuous allocation; pass `$detail$n_int`
#'   to assess the operational recommendation. A named vector is matched to
#'   `frame$stratum`; an unnamed vector is positional. In fixed-take
#'   multistage mode these are ultimate-unit sizes. Explicit adopted sizes must
#'   correspond to whole PSU counts; omission from a fitted result retains its
#'   continuous first-stage allocation.
#' @param measures Optional long indicator table for joint assessment. See the
#'   `measures` argument to [n_alloc()]. It must be supplied with `targets` in
#'   the default method and is recovered automatically from a fitted result.
#' @param targets Optional long precision-requirement table for joint
#'   assessment. See the `targets` argument to [n_alloc()]. It must be supplied
#'   with `measures` in the default method and is recovered automatically from
#'   a fitted result.
#' @param objective Optional objective components to report alongside the
#'   targets. See the `objective` argument to [n_alloc()]. Recovered
#'   automatically from a fitted budget-objective result.
#' @param budget Optional budget, used only to report the residual against the
#'   assessed allocation's cost. Requires `objective`.
#' @param domains Character vector of column names in `frame` to treat as
#'   domain identifiers, or `NULL` (default) for no domains.
#' @param alpha Significance level, default 0.05.
#' @param deff Design effect multiplier (> 0).
#' @param resp_rate Expected response rate, in (0, 1\]. Default 1.
#' @param unit_cost Optional scalar or length-`nrow(frame)` vector of
#'   per-stratum unit costs, overriding `frame$unit_cost`. Fixed-take
#'   multistage joint assessment instead uses the stage costs stored in
#'   `frame`.
#' @param min_n_stratum Optional minimum sample size per stratum, applied as the
#'   lower bound the assessment reports against in `$bounds`. It is the same
#'   argument [n_alloc()] takes, so a design and its assessment can be held
#'   to one floor. Joint assessment only; supplying it without `measures` and
#'   `targets` is an error.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "alloc"`. For joint
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
#'   through cost instead: the targets stay as specified and the assessed
#'   allocation's cost becomes the budget, because pinning achieved precision
#'   as hard targets would over-constrain a design that already spends its
#'   whole budget.
#'
#' @seealso [n_alloc()].
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
  objective = NULL,
  budget = NULL,
  domains = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  unit_cost = NULL,
  min_n_stratum = NULL,
  plan = NULL
) {
  .plan <- .merge_plan_args(plan, prec_alloc.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_alloc.default, c(.plan, list(...))))
  .check_unused_dots(...)
  joint_any <- !is.null(measures) || !is.null(targets) || !is.null(objective)
  if (joint_any) {
    if (is.null(measures)) {
      stop("joint precision assessment requires 'measures'", call. = FALSE)
    }
    if (is.null(targets) && is.null(objective)) {
      stop("'measures' and 'targets' must be supplied together", call. = FALSE)
    }
    if (!is.null(budget) && is.null(objective)) {
      stop("'budget' requires 'objective' for joint precision assessment",
           call. = FALSE)
    }
    if (!is.null(budget)) check_scalar(budget, "budget")
    if (!is.null(min_n_stratum)) check_scalar(min_n_stratum, "min_n_stratum")
    if (!is.null(domains)) {
      stop("use 'targets$domain' for joint constrained allocation, not 'domains'",
           call. = FALSE)
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
      min_n_stratum = min_n_stratum
    ))
  }
  if (!is.null(min_n_stratum)) {
    # Aggregate precision does not depend on a floor that the supplied
    # allocation either already meets or does not; only the joint mode
    # reports per-stratum bound diagnostics.
    stop(
      "'min_n_stratum' applies to joint precision assessment; supply 'measures' and 'targets', or drop it",
      call. = FALSE
    )
  }
  check_alpha(alpha)
  check_deff(deff)
  check_resp_rate(resp_rate)

  prep <- .alloc_prepare_frame(frame, domains = domains, unit_cost = unit_cost)
  H <- length(prep$N_h)

  if (!is.numeric(n) || anyNA(n) || any(!is.finite(n)) || length(n) != H) {
    stop("'n' must be a finite numeric vector with length nrow(frame)",
         call. = FALSE)
  }
  if (any(n <= 0)) {
    stop("all 'n' elements must be positive", call. = FALSE)
  }

  metrics <- .alloc_metrics(
    N_h = prep$N_h, S_h = prep$S_h, mean_h = prep$mean_h, n_h = n,
    alpha = alpha, deff = deff, resp_rate = resp_rate,
    cost_h = prep$cost_h
  )

  detail <- .alloc_detail(
    prep = prep, n_h = n,
    m_h = rep(NA_real_, H), M_h = prep$N_h,
    resp_rate = resp_rate, deff = deff
  )

  .new_svyplan_prec(
    se = metrics$se,
    moe = metrics$moe,
    cv = metrics$cv,
    type = "alloc",
    params = list(
      frame = frame, n = n, alpha = alpha, deff = deff,
      resp_rate = resp_rate, cost_h = prep$cost_h,
      domain_cols = prep$domain_cols,
      achieved = list(n = sum(n), cv = metrics$cv, cost = metrics$cost)
    ),
    detail = detail
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
      unit_cost = p$unit_cost,
      min_n_stratum = p$min_n_stratum,
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
    stop("allocation detail does not contain stratum sample sizes",
         call. = FALSE)
  }

  args <- list(
    frame = p$frame,
    n = n_h,
    domains = p$domain_cols,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate
  )
  if (!.alloc_is_cluster(p$frame)) {
    args$unit_cost <- p$cost_h
  }
  do.call(prec_alloc.default, .roundtrip_args(args, list(...), prec_alloc.default))
}

#' @keywords internal
#' @noRd
.alloc_prepare_frame <- function(frame, domains = NULL, unit_cost = NULL) {
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
    if (!is.numeric(S_h) || anyNA(S_h) || any(!is.finite(S_h)) || any(S_h < 0)) {
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
    warning("all 'sd' values are zero; allocation has no variability to distribute",
            call. = FALSE)
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
    warning("all 'mean' (or 'p') values are zero; CV will be Inf",
            call. = FALSE)
  }

  if ("cost" %in% names(frame)) {
    stop(
      "the per-stratum cost column is 'unit_cost'; 'cost' is the total field cost",
      call. = FALSE
    )
  }

  if (!is.null(unit_cost)) {
    if (!is.numeric(unit_cost) || anyNA(unit_cost) || any(!is.finite(unit_cost)) || any(unit_cost <= 0)) {
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
    if (!is.numeric(cost_h) || anyNA(cost_h) || any(!is.finite(cost_h)) ||
        any(cost_h <= 0)) {
      stop("'unit_cost' must contain positive finite values", call. = FALSE)
    }
  } else {
    cost_h <- rep(1, nrow(frame))
  }

  max_weight <- rep(NA_real_, nrow(frame))
  if ("max_weight" %in% names(frame)) {
    max_weight <- frame$max_weight
    bad <- !is.na(max_weight) & (!is.finite(max_weight) | max_weight < 1)
    if (any(bad)) {
      stop("'max_weight' must be >= 1 and finite when provided",
           call. = FALSE)
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
        lab <- paste(unlist(lapply(domain_values[i, , drop = FALSE],
                                   as.character)), collapse = ":")
        stop(
          sprintf("duplicate stratum labels in domain '%s': %s", lab,
                  paste(s_d[duplicated(s_d)], collapse = ", ")),
          call. = FALSE
        )
      }
    }
  } else {
    if (anyDuplicated(stratum)) {
      dups <- stratum[duplicated(stratum)]
      stop(
        sprintf("duplicate stratum labels: %s", paste(unique(dups), collapse = ", ")),
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
    max_weight = as.numeric(max_weight),
    take_all = as.logical(take_all),
    stratum = stratum,
    domain_cols = domain_cols,
    domain_idx = domain_idx,
    domain_values = domain_values
  )
  .alloc_cluster_prep(prep, frame, unit_cost)
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
.alloc_cluster_prep <- function(prep, frame, unit_cost) {
  if (!.alloc_is_cluster(frame)) {
    orphan <- intersect(c("var_ratio_psu", "n_per_psu", "cost_psu", "cost_ssu"),
                        names(frame))
    if (length(orphan) > 0L) {
      stop(
        sprintf("cluster column(s) %s require a 'icc_psu' column",
                paste(sQuote(orphan), collapse = ", ")),
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

  var_ratio <- if ("var_ratio_psu" %in% names(frame)) frame[["var_ratio_psu"]] else rep(1, H)
  if (!is.numeric(var_ratio) || anyNA(var_ratio) || any(!is.finite(var_ratio)) || any(var_ratio <= 0)) {
    stop("'var_ratio_psu' must contain positive finite values", call. = FALSE)
  }

  has_costs <- any(c("cost_psu", "cost_ssu") %in% names(frame))
  if (has_costs) {
    if (!all(c("cost_psu", "cost_ssu") %in% names(frame))) {
      stop("'cost_psu' and 'cost_ssu' must be supplied together",
           call. = FALSE)
    }
    for (col in c("cost_psu", "cost_ssu")) {
      cc <- frame[[col]]
      if (!is.numeric(cc) || anyNA(cc) || any(!is.finite(cc)) ||
          any(cc <= 0)) {
        stop(sprintf("'%s' must contain positive finite values", col),
             call. = FALSE)
      }
    }
  }

  n_per_psu_h <- rep(NA_real_, H)
  if ("n_per_psu" %in% names(frame)) {
    ps <- frame$n_per_psu
    if (!is.numeric(ps) || any(!is.na(ps) & (!is.finite(ps) | ps < 1))) {
      stop("'n_per_psu' must contain values >= 1 (NA for cost-optimal)",
           call. = FALSE)
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
    n_per_psu_h[need_opt] <- sqrt(
      frame$cost_psu[need_opt] / frame$cost_ssu[need_opt] *
        (1 - icc[need_opt]) / icc[need_opt]
    )
  }

  too_big <- n_per_psu_h > prep$N_h
  if (any(too_big & !need_opt)) {
    stop(
      sprintf("'n_per_psu' exceeds the stratum population for: %s",
              paste(prep$stratum[too_big & !need_opt], collapse = ", ")),
      call. = FALSE
    )
  }
  if (any(too_big & need_opt)) {
    n_per_psu_h[too_big & need_opt] <- prep$N_h[too_big & need_opt]
    warning(
      "cost-optimal 'n_per_psu' exceeds the stratum population; clamped to 'N'",
      call. = FALSE
    )
  }

  prep$S_raw_h <- prep$S_h
  prep$S_h <- prep$S_h * sqrt(var_ratio * (1 + icc * (n_per_psu_h - 1)))
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

#' @keywords internal
#' @noRd
.alloc_bounds <- function(N_h, max_weight, take_all, min_n_stratum = NULL) {
  H <- length(N_h)
  m_h <- pmin(rep(1, H), N_h)
  M_h <- as.numeric(N_h)

  if (!is.null(min_n_stratum)) m_h <- pmax(m_h, min_n_stratum)

  has_wmax <- !is.na(max_weight)
  if (any(has_wmax)) {
    m_h[has_wmax] <- pmax(m_h[has_wmax], N_h[has_wmax] / max_weight[has_wmax])
  }

  if (any(take_all)) {
    m_h[take_all] <- N_h[take_all]
    M_h[take_all] <- N_h[take_all]
  }

  if (any(m_h > M_h + 1e-8)) {
    stop("constraints are infeasible: lower bounds exceed stratum population",
         call. = FALSE)
  }

  list(m_h = m_h, M_h = M_h)
}

#' @keywords internal
#' @noRd
.alloc_metrics <- function(N_h, S_h, mean_h, n_h, alpha, deff,
                           resp_rate, cost_h) {
  n_net <- n_h * resp_rate
  n_eff <- n_net / deff
  W_h <- N_h / sum(N_h)

  term <- numeric(length(n_h))
  good <- n_eff > 0
  fpc <- pmax(0, 1 - n_net / N_h)
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
  moe <- qnorm(1 - alpha / 2) * se

  cv <- NA_real_
  if (!all(is.na(mean_h))) {
    if (anyNA(mean_h)) {
      stop("'mean' (or 'p') must be complete to compute aggregate CV",
           call. = FALSE)
    }
    ybar <- .aggregate_mean(W_h, mean_h)
    cv <- if (ybar == 0) Inf else se / abs(ybar)
  }

  list(se = se, moe = moe, cv = cv, cost = sum(n_h * cost_h))
}

#' @keywords internal
#' @noRd
.alloc_detail <- function(prep, n_h, m_h, M_h, resp_rate, deff,
                          mode = "n", budget = NULL,
                          lo_i = NULL, hi_i = NULL) {
  if (is.null(lo_i)) {
    lo_i <- as.integer(ceiling(pmax(ifelse(is.na(m_h), 0, m_h), 0) - 1e-9))
  }
  if (is.null(hi_i)) {
    hi_i <- as.integer(floor(M_h + 1e-9))
  }
  n_int <- switch(
    mode,
    cv = pmin(pmax(as.integer(ceiling(n_h - 1e-9)), lo_i), hi_i),
    budget = .round_within_budget(n_h, prep$cost_h, budget, lo_i, hi_i,
                                  prep$N_h, prep$S_h),
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
  if (isTRUE(prep$cluster)) {
    out$n_per_psu <- prep$n_per_psu_h
    out$n_psu <- n_h / prep$n_per_psu_h
    out$n_psu_int <- as.integer(ceiling(out$n_psu))
  }
  if (any(prep$take_all)) out$take_all <- prep$take_all
  if (!all(is.na(prep$mean_h))) out$mean <- prep$mean_h
  out
}

#' Integerize a budget-mode allocation without exceeding the budget
#'
#' Floors the continuous allocation (respecting integer lower bounds),
#' then greedily adds units where the variance reduction per unit cost
#' is largest while the budget allows.
#' @keywords internal
#' @noRd
.round_within_budget <- function(n_h, cost_h, budget, lower, upper,
                                 N_h, S_h) {
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
.alloc_operational_element <- function(prep, detail, alpha, deff, resp_rate) {
  n_int <- detail$n_int
  m <- .alloc_metrics(
    N_h = prep$N_h, S_h = prep$S_h, mean_h = prep$mean_h, n_h = n_int,
    alpha = alpha, deff = deff, resp_rate = resp_rate, cost_h = prep$cost_h
  )
  list(
    n = sum(n_int),
    cost = sum(n_int * prep$cost_h),
    se = m$se, moe = m$moe, cv = m$cv
  )
}

#' Operational (integer) design for a stratified two-stage allocation
#'
#' Chooses a whole per-stratum take b_h (floor/ceiling candidate with the
#' lowest variance-cost product when stage costs are known) and a whole
#' PSU count a_h per stratum. In cv mode a_h matches or beats each
#' stratum's continuous variance contribution; in budget mode PSUs are
#' removed/added greedily so the field cost sum(a_h * (cost_psu +
#' cost_ssu * b_h)) never exceeds the budget; in n mode a_h approximates
#' the continuous element total. Returns the design plus recomputed
#' metrics and updates the detail integers.
#' @keywords internal
#' @noRd
.alloc_operational_cluster <- function(prep, detail, n_h, mode, budget,
                                       alpha, deff, resp_rate,
                                       target_cv = NULL) {
  H <- length(n_h)
  ps <- prep$n_per_psu_h
  icc <- prep$icc_psu_h
  var_ratio <- prep$var_ratio_psu_h
  S_raw <- prep$S_raw_h
  W <- prep$N_h / sum(prep$N_h)
  lo_i <- as.integer(ceiling(detail$.lower - 1e-9))
  hi_i <- as.integer(floor(detail$.upper + 1e-9))
  b_h <- vapply(seq_len(H), function(h) {
    centre <- max(1L, as.integer(round(ps[h])))
    cand <- unique(c(
      max(1L, as.integer(floor(ps[h]))),
      max(1L, as.integer(ceiling(ps[h]))),
      seq.int(max(1L, centre - 50L), min(hi_i[h], centre + 50L)),
      1L
    ))
    # A take-all or other exact bound must be representable as a product
    # of whole PSUs and a whole, common take.  Including the exact bound
    # guarantees a feasible divisor when it is an integer.
    if (lo_i[h] == hi_i[h]) {
      cand <- unique(c(cand, hi_i[h]))
    }
    cand <- cand[cand >= 1L & cand <= hi_i[h]]
    feasible <- vapply(cand, function(b) {
      ceiling(lo_i[h] / b) <= floor(hi_i[h] / b)
    }, logical(1L))
    cand <- cand[feasible]
    if (length(cand) == 0L) {
      stop(sprintf("no whole-cluster design satisfies the bounds for stratum '%s'",
                   prep$stratum[h]), call. = FALSE)
    }

    # Fixed cluster sizes are rounded to the nearest whole size whenever
    # that is compatible with the hard element bounds.  Only exact bounds
    # (notably take-all) may require a divisor farther away.
    if (isTRUE(prep$n_per_psu_fixed_h[h])) {
      near <- cand[cand %in% unique(c(floor(ps[h]), ceiling(ps[h])))]
      if (length(near) > 0L) cand <- near
    }
    if (isTRUE(prep$has_stage_costs)) {
      score <- (var_ratio[h] * (1 + icc[h] * (cand - 1))) *
        (prep$cost_psu_h[h] / cand + prep$cost_ssu_h[h])
      cand[which.min(score)]
    } else {
      cand[which.min(abs(cand - ps[h]))]
    }
  }, integer(1L))

  e_target <- NULL
  if (mode == "n") {
    # Round the requested element total first, then factor each stratum's
    # integer element count into whole PSUs and a common whole take. Every
    # positive integer has divisor 1, so this cannot violate bounds merely
    # because a preferred cluster size does not divide the allocated elements.
    #
    # A fixed n_per_psu is exempt: it keeps the take chosen above and lets its
    # element count move to the nearest whole multiple instead. Factoring would
    # otherwise hand a fixed take of 10 the nearest divisor of the rounded
    # total, which is 1 whenever that total is prime.
    b_free <- b_h
    e_target <- .round_oric_bounded(n_h, lo_i, hi_i)
    b_h <- vapply(seq_len(H), function(h) {
      if (isTRUE(prep$n_per_psu_fixed_h[h])) {
        return(b_free[h])
      }
      e <- e_target[h]
      root <- seq_len(max(1L, as.integer(floor(sqrt(e)))))
      small <- root[e %% root == 0L]
      divisors <- sort(unique(c(small, e %/% small)))
      if (isTRUE(prep$has_stage_costs)) {
        score <- (var_ratio[h] * (1 + icc[h] * (divisors - 1))) *
          (prep$cost_psu_h[h] / divisors + prep$cost_ssu_h[h])
        return(divisors[which.min(score)])
      }
      divisors[which.min(abs(divisors - ps[h]))]
    }, integer(1L))
  }

  a_min <- as.integer(ceiling(lo_i / b_h))
  a_max <- as.integer(floor(hi_i / b_h))
  if (any(a_min > a_max)) {
    stop("no whole-cluster design satisfies the integer allocation bounds",
         call. = FALSE)
  }

  psu_cost <- if (isTRUE(prep$has_stage_costs)) {
    prep$cost_psu_h + prep$cost_ssu_h * b_h
  } else {
    rep(NA_real_, H)
  }
  S_op <- S_raw * sqrt(var_ratio * (1 + icc * (b_h - 1)))
  Cj <- W^2 * S_op^2 * deff / (b_h * resp_rate)

  if (mode == "cv") {
    a_h <- as.integer(ceiling(
      n_h * (1 + icc * (b_h - 1)) / (1 + icc * (ps - 1)) / b_h - 1e-9
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
        stop("no whole-cluster design satisfies both the budget and lower bounds",
             call. = FALSE)
      }
      loss <- Cj[cand] * (1 / (a_h[cand] - 1L) - 1 / a_h[cand]) /
        psu_cost[cand]
      j <- cand[which.min(loss)]
      a_h[j] <- a_h[j] - 1L
    }
    repeat {
      spare <- budget - sum(a_h * psu_cost)
      cand <- which(psu_cost <= spare + 1e-8 & a_h < a_max)
      if (length(cand) == 0L) {
        break
      }
      gain <- Cj[cand] * (1 / a_h[cand] - 1 / (a_h[cand] + 1L)) /
        psu_cost[cand]
      j <- cand[which.max(gain)]
      a_h[j] <- a_h[j] + 1L
    }
  } else {
    # Exact where the take divides the target; the nearest whole number of
    # PSUs where a fixed take does not, which is what moves the total.
    a_h <- as.integer(round(e_target / b_h))
    a_h <- pmin(pmax(a_h, a_min), a_max)
  }

  e_h <- a_h * b_h
  m <- .alloc_metrics(
    N_h = prep$N_h, S_h = S_op, mean_h = prep$mean_h, n_h = e_h,
    alpha = alpha, deff = deff, resp_rate = resp_rate, cost_h = prep$cost_h
  )
  detail$n_int <- as.integer(e_h)
  detail$n_psu_int <- a_h
  detail$n_per_psu_int <- b_h

  if (any(e_h < lo_i | e_h > hi_i)) {
    stop("internal error: operational cluster allocation violates its bounds",
         call. = FALSE)
  }
  if (!is.null(target_cv)) {
    op_cv <- if (length(prep$domain_idx) == 0L) {
      m$cv
    } else {
      prep_op <- prep
      prep_op$S_h <- S_op
      .alloc_domain_cv_max(prep_op, e_h, alpha, deff, resp_rate)
    }
    if (!is.finite(op_cv) || op_cv > target_cv + 1e-8) {
      stop(
        sprintf("no whole-cluster design found that attains target CV %.4g under the integer bounds",
                target_cv),
        call. = FALSE
      )
    }
  } else {
    op_cv <- m$cv
  }
  list(
    operational = list(
      n = sum(e_h),
      cost = if (isTRUE(prep$has_stage_costs)) sum(a_h * psu_cost)
             else NA_real_,
      se = m$se, moe = m$moe, cv = op_cv
    ),
    detail = detail
  )
}

#' @keywords internal
#' @noRd
.alloc_binding_label <- function(n_h, m_h, M_h) {
  on_lower <- abs(n_h - m_h) < 1e-6
  on_upper <- abs(n_h - M_h) < 1e-6
  if (any(on_lower)) return("lower_bound")
  if (any(on_upper)) return("upper_bound")
  "none"
}

#' @keywords internal
#' @noRd
.alloc_domain_cv_max <- function(prep, n_h, alpha, deff, resp_rate) {
  cvs <- vapply(
    prep$domain_idx,
    function(idx) {
      .alloc_metrics(
        N_h = prep$N_h[idx], S_h = prep$S_h[idx],
        mean_h = prep$mean_h[idx], n_h = n_h[idx],
        alpha = alpha, deff = deff, resp_rate = resp_rate,
        cost_h = prep$cost_h[idx]
      )$cv
    },
    numeric(1)
  )
  max(cvs)
}

#' @keywords internal
#' @noRd
.alloc_domain_summary <- function(prep, n_h, alpha, deff, resp_rate) {
  if (length(prep$domain_idx) == 0L) return(NULL)

  out <- vector("list", length(prep$domain_idx))
  dom_names <- names(prep$domain_idx)
  for (i in seq_along(prep$domain_idx)) {
    idx <- prep$domain_idx[[i]]
    met <- .alloc_metrics(
      N_h = prep$N_h[idx], S_h = prep$S_h[idx],
      mean_h = prep$mean_h[idx], n_h = n_h[idx],
      alpha = alpha, deff = deff, resp_rate = resp_rate,
      cost_h = prep$cost_h[idx]
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
    row$.cv <- met$cv
    row$.cost <- met$cost
    out[[i]] <- row
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
