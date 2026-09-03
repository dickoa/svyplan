#' Joint constrained allocation, and the design a budget buys
#'
#' The generalized mode of [n_alloc()] and [prec_alloc()], reached by
#' supplying `measures` and `targets` rather than a scalar `n`, `cv` or
#' `budget`. It sizes one stratified design against several indicators and
#' domains at once, either as the cheapest design meeting every requirement
#' or as the best design a fixed budget can buy. This page holds the input
#' contract, the variance model, the infeasibility messages and the
#' `predict()` frontier. The arguments themselves are documented on
#' [n_alloc()] and [prec_alloc()].
#'
#' @details
#' ## The target table
#'
#' Supplying both `measures` and `targets` requests a minimum-cost allocation
#' that meets all indicator-domain precision requirements simultaneously.
#' This mode is mutually exclusive with scalar `n`, `cv`, `budget`, `domains`,
#' and an explicit `alloc` rule, because the target table defines the
#' constraints and the allocation is solved jointly.
#'
#' CV targets apply to domain means or totals (the CV is the same on either
#' scale). MOE targets are absolute margins of error for domain means. An
#' `rmoe` target is a margin of error for the domain mean stated as a
#' fraction of it, which is `q` times the CV of the same
#' estimate. The solver takes it through the CV branch and reports
#' `.achieved` and `.sensitivity` back in `rmoe` units. The
#' quantile `q` is `qnorm(1 - alpha / 2)` by default and
#' `qt(1 - alpha / 2, df)` when `df` is supplied. The
#' variance model uses Wald/linearized variances, the package's gross-sample
#' response-rate convention, design effects, and finite population correction.
#' Domain classifications may overlap. For example, separate `region` and
#' `residence` targets can use the same atomic `region x residence` frame rows.
#'
#' ## Fixed-take multistage joint designs
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
#' design. The later-stage takes are fixed inputs rather than decision
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
#' `n_per_psu`, \eqn{q_h} for `n_per_ssu`, and \eqn{S^2_{hk}}{S^2_hk} for `var` (or
#' `sd^2`, or `p * (1 - p)`). Each stage's *realized* take is what inflates the
#' variance, so with ultimate-unit response \eqn{r_{hk}}{r_hk} and SSU response
#' \eqn{r^{ssu}_{hk}}{r^ssu_hk} write
#' \eqn{\tilde m_h=m_hr^{ssu}_{hk}}{m~_h=m_h r^ssu_hk} and
#' \eqn{\tilde q_h=q_hr_{hk}}{q~_h=q_h r_hk}. At two stages the ultimate-unit
#' rate acts on \eqn{m_h} directly. The fixed variance multiplier is
#' \eqn{D_{hk}=k_{1,hk}(1+\delta_{1,hk}(m_hr_{hk}-1))}{D_hk=k_(1,hk)(1+delta_(1,hk)(m_h r_hk-1))} at two stages and
#' \eqn{D_{hk}=k_{1,hk}\delta_{1,hk}\tilde m_h\tilde q_h+
#' k_{2,hk}(1+\delta_{2,hk}(\tilde q_h-1))}{D_hk=k_(1,hk) delta_(1,hk)m~_hq~_h+ k_(2,hk)(1+delta_(2,hk)(q~_h-1))} at three stages. At three stages
#' `var_ratio_ssu` is not a free parameter: `var_ratio_psu` rescales the components' unit
#' variance to the analysis variable and `var_ratio_ssu` does the same for the
#' within-PSU part, so
#' \eqn{k_{2}=k_{1}(1-\delta_{1})}{k_2=k_1(1-delta_1)}. That identity is what makes
#' \eqn{D_{hk}}{D_hk} collapse to \eqn{k_{1}}{k_1} at \eqn{m_h=q_h=1}, where no
#' clustering is left to inflate anything. Leaving `var_ratio_ssu` out of `measures`
#' applies it. Supplying a value overrides it, which is only meaningful when
#' the two stages' ratios come from different decompositions. With **gross**
#' ultimate-unit take \eqn{t_h=m_h} or \eqn{m_hq_h}, first-stage decision
#' \eqn{a_h}, combined response
#' \eqn{R_{hk}=r^{psu}_{hk}r^{ssu}_{hk}r_{hk}}{R_hk=r^psu_hk r^ssu_hk r_hk},
#' and extra design effect \eqn{d_{hk}}{d_hk}, the implemented
#' total-variance contribution is
#' \deqn{N_h^2 S_{hk}^2 D_{hk}d_{hk}/(R_{hk}a_ht_h)
#'       - N_h S_{hk}^2 D_{hk}d_{hk}.}{N_h^2 S_hk^2 D_hkd_hk/(R_hka_ht_h) - N_h S_hk^2 D_hkd_hk.}
#' Thus response inflation is applied once. Takes and costs stay gross, since
#' a unit is issued and paid for either way, while the variance reads what
#' the design realizes.
#'
#' ## Response acts at whichever stage it happens
#'
#' The three rates are not interchangeable and none is recoverable from the
#' others, so each names the stage it acts on. `resp_rate_psu` divides the
#' whole requirement, because losing a PSU is a pure sample-size loss.
#' `resp_rate_ssu` divides it and shrinks the realized SSU count per PSU.
#' `resp_rate` divides it and enters the \eqn{\delta}{delta} bracket, because
#' it shrinks the realized final-stage take and so changes the clustering
#' penalty itself. Reading a gross take in that bracket would charge the
#' between-PSU component a penalty only whole-cluster loss produces.
#'
#' All three go in `measures`, or in `frame` as a stratum default. A rate
#' naming a stage the design does not have is an error rather than a
#' silently ignored column. This is the same decomposition [n_cluster()] and
#' cluster-mode `n_alloc()` use, so a design described identically to any of
#' them sizes identically in all of them.
#'
#' Leave `deff = 1` unless it
#' represents a source not already captured by the stage deltas and var_ratio factors.
#' Deltas may include 0 and 1. The var_ratio values must be positive finite.
#' All stage
#' populations and fixed takes must be positive whole numbers.
#'
#' The continuous result is a KKT-certified global optimum of the convex,
#' fixed-coefficient problem. `$detail$n_int` is a deterministic, feasible,
#' locally cleaned operational recommendation. It is not claimed to be the
#' globally optimal integer allocation. Full achieved precision is in
#' `$constraints` and `$operational$constraints`. Use [prec_alloc()] to assess
#' the fitted, operational, or a modified allocation.
#'
#' A constraint table records `.metric`, `.target`, `.achieved`, `.ratio`, and
#' `.residual = .ratio - 1`. A `.pass` means `.ratio <= 1 + .tolerance`.
#' `.binding` identifies a numerically active target. `.multiplier` is the
#' continuous Lagrange multiplier and `.sensitivity` is the local derivative
#' of minimum variable cost with respect to the target. Both are `NA` for a
#' pure precision assessment. Multipliers are also `NA` when duplicate
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
#' \deqn{Q(n) = \sum_j w_j \, \mathrm{cv}_j^2(n),}{Q(n) = sum_j w_j cv_j^2(n),}
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
#' \eqn{Q \le Q^{*}}{Q <= Q^*} returns the same allocation at a cost equal to the budget.
#'
#' Infeasibility takes three distinct forms, each with its own message: the
#' targets may be unattainable even at the stratum upper bounds. They may be
#' attainable but unaffordable, in which case the cheapest target-feasible cost
#' and the shortfall are reported. Or the continuous problem may be feasible
#' while no whole-unit allocation meets every target inside the budget.
#'
#' `ceiling()` is not a safe integer start under a budget, so the operational
#' design is built from several roundings of the continuous optimum, repaired
#' and trimmed to be both target-feasible and affordable, then improved by
#' spending residual budget and by pairwise exchanges. It satisfies
#' `$operational$cost <= budget` and `all($operational$constraints$.pass)`, and
#' it is a feasible, locally improved recommendation rather than a globally
#' optimal integer allocation. The greedy construction is incomplete and can
#' report failure on a problem that does have a feasible integer point.
#'
#' [predict()] on a fixed-budget result varies `budget` in `newdata` and
#' returns the cost-versus-objective frontier, reusing the same solves.
#'
#'
#' @seealso [n_alloc()] for the allocator and its arguments, [prec_alloc()]
#'   for the assessment direction, [predict.svyplan()] for the frontier, and
#'   [n_multi()] for the unstratified multi-indicator size this refines.
#'
#' @family stratified design functions
#'
#' @name n_alloc-generalized
NULL
