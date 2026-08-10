# svyplan 0.11.0

Initial CRAN release.

## Sample size determination

* `n_prop()`: sample size for a proportion (Wald, Wilson, log-odds, and
  Korn-Graubard beta methods).
  The Wilson, log-odds, and beta sizes are obtained by inverting the interval
  half-width that `prec_prop()` reports, so each `n_prop()` method is an
  exact inverse of its `prec_prop()` counterpart. All four methods share
  one variance, `deff * N / (N - 1) * p (1 - p) (1 / n - 1 / N)`, reading it
  through the effective size `n_eff = n_net / (deff * fpc)`, so the four
  sizes converge on each other as the margin of error shrinks and differ
  materially only when the expected number of positive counts is small.
  `?n_prop` explains when that difference is worth acting on: for sizing it
  rarely is, for assessing an achieved design it often is. All four methods
  apply the finite population correction, and the three interval-based
  (Wilson, log-odds, beta) require a margin of error below 0.5,
  the widest half-width any sample size can produce.
* `method = "beta"` is the Korn and Graubard (1998) interval: the
  Clopper-Pearson limits evaluated at the effective sample size, for a
  proportion whose expected number of positive counts is small. It matches
  `survey::svyciprop(method = "beta")` exactly on a design where the two see
  the same effective size. A `df` argument supplies the degrees of freedom
  of the planned variance estimator (typically sampled PSUs minus strata)
  and widens the interval for a variance estimated from few clusters, which
  is equation (2.2) of the paper; omitting it applies no adjustment. The
  interval is asymmetric and stays inside `[0, 1]` by construction rather
  than by clamping, so `confint()` is the way to read it: `$moe` is half the
  interval's length and `p` plus or minus it is not the pair of limits.
  Reach for it when the expected count of positive cases is small in
  absolute terms, which is the regime where normality fails however large
  the sample is.
* `min_cases` sizes a proportion against a minimum expected number of
  positive cases, the constraint that actually binds for a rare outcome in
  a small domain, where the requirement is a count rather than a margin of
  error. The size it demands is `min_cases / (p * resp_rate)` and the
  result is the larger of that and the size precision asks for, with
  `binding` naming the one that decided. `deff` has no part in it: a design
  effect says how precisely the proportion is estimated, while the number
  of cases that turn up depends on the sample size and the prevalence
  alone. The response rate does enter, since the cases are counted among
  respondents and the returned size is gross. `n_multi()` takes the same
  constraint as an optional `min_cases` column, a floor on that row's own
  size, so the row that ends up largest is still the binding one whichever
  constraint raised it.
* Every proportion result reports `expected_cases`, the count
  `n * resp_rate * p` the design expects to yield. This is the number the
  choice between the interval methods turns on, `"beta"` earning its place
  when the expected count is small in absolute terms, and it is reported
  by `n_prop()` and `prec_prop()` alike.
* The multi-indicator functions take the same four proportion methods as
  `n_prop()`, including `"beta"`, so a table can mix an ordinary Wald
  indicator with a rare one sized from a Korn-Graubard interval. A `df`
  column sets each row's interval quantile. A df is a property of the
  design rather than of an indicator, so the column earns its place when
  the rows are domains: `design_df(alloc)$domains$.df` supplies one value
  per domain.
* `design_df()`: the degrees of freedom a planned design's variance
  estimator will have, the planning analogue of `survey::degf()`. Counted
  from the whole-unit design a plan reports rather than from its continuous
  optimum: sampled PSUs minus strata for a clustered design, units minus
  strata for an element one, phase-2 units minus strata for a two-phase
  one. A `take_all` stratum is a census in element mode and drops out of
  both terms, while in cluster mode it marks an element-level census that
  leaves the PSU stage sampling, so the stratum still contributes. The
  result is a `svyplan_df`, a numeric scalar usable as any `df` argument
  and carrying `$strata` and `$domains` tables whose entries sum to it.
  Per-domain df is exact rather than approximate, the allocation API
  expressing a domain as a union of whole strata; an analytic domain
  cutting across strata is not derivable from a plan and is documented as
  such. A stratum holding a single PSU supports no within-stratum variance
  estimate and warns, naming the stratum, both here and at construction in
  `n_alloc()`, as does a whole design that rests on one PSU. A small but
  positive df is reported and never warned about.
* `df` sets the interval quantile in every sizing and precision function
  that builds an interval, `n_prop()`, `prec_prop()`, `n_mean()`,
  `prec_mean()`, `n_alloc()`, `prec_alloc()` and the multi-indicator
  family, switching it from the normal to `t` on that many degrees of
  freedom. It is opt-in rather than derived: in `n_cluster()` the df
  depends on the PSU count being solved for, so an automatic version would
  need a fixed point. `design_df()` is where the number comes from, and
  `svyplan(df = )` carries it as a design default. The Korn-Graubard
  `"beta"` interval reaches the same widening by scaling the effective
  sample size rather than the quantile; the two agree in direction and
  magnitude, and both hold `moe = q * se`. The coefficient of variation is
  not uniformly invariant to `df`: only `"wald"` and the mean engine build
  `se` without a quantile in it. A joint allocation reads a scalar `df` or
  a per-target `df` column, with one quantile shared by a constraint's
  variance ceiling and the sensitivity reported for it; a `cv`-metric
  constraint carries no quantile, so an all-`cv` problem is unchanged by
  `df` in everything but its reported `moe`. The power functions do not
  accept `df` and say why: their quantile is a normal deviate for an
  alternative, not an interval half-width.
* `n_twophase()` and `prec_twophase()`: allocation and precision for a
  two-phase (double) sample. One allocator covers two designs that are
  usually treated separately. Double sampling for stratification
  subsamples every phase-2 stratum; nonresponse follow-up carries the
  phase-1 respondents through untouched and subsamples only the
  nonrespondents. They differ only in which strata are marked
  `take_all`, and the follow-up case reduces to the standard
  `sqrt(c1 / (c2 * theta))`. The optimum is closed form, Neyman-shaped in
  `sd / sqrt(unit_cost)` but scaled by the between-stratum variance, so
  weak stratification pushes every fraction toward 1. Fractions above 1
  are truncated one stratum at a time, since pinning one changes the
  scale for the rest. Every result reports the single-phase design that
  skips phase 1 entirely, because two-phase sampling is not always an
  improvement and the comparison is the only thing that says which way it
  falls. Neither phase need be a simple random sample: `phase1_deff`
  inflates the between-stratum component, a `deff` frame column inflates
  each stratum's within-stratum residual, and `single_deff` describes the
  comparator. They are design effects for different variables, so a
  single one on the combined variance would be a different and wrong
  model. Response works the same way: `resp_rate` is the phase-1 rate, a
  `resp_rate` frame column the phase-2 rate per stratum, and each divides
  its own component. Phase 2 cannot outdraw the pool phase 1 classified,
  so the subsampling fractions are capped at the phase-1 rate. Issued and
  expected-responding counts are reported separately, since neither is an
  effective sample size and a two-phase design has no single one. The
  division by a response rate is an expected-information calculation that
  assumes ignorable response within the supplied strata, not a
  nonresponse-bias adjustment. `n_phase1` fixes the phase-1 size when the
  screener has already run or phase 1 is an existing survey, leaving only
  the subsampling depth to choose; the relative allocation is unchanged
  and only its scale moves. Fixing it brings its own two infeasibilities,
  both refused with the binding quantity named: a budget that the
  screening and the take-all strata exhaust before phase 2 begins, and a
  target below the variance floor left once every classified unit is
  carried into phase 2. With `N` set, the phase-1 finite population
  correction applies component by component, so a phase-1 census removes
  the between-stratum term but leaves the residual phase 2 still has to
  subsample. The census term it subtracts carries the design effects, as
  one multiplies an SRSWOR variance everywhere else here, but not the
  response rates: only a design that actually measures every unit reaches
  zero variance, while issuing to every unit and losing some of them to
  nonresponse does not. The single-phase comparator follows the same
  convention, is capped at `N` when a target sits below its own census
  floor, and does not win a `cv` comparison on cost unless it reaches the
  target. `nu` is a share of the phase-1 units issued, so `prec_twophase()`
  refuses one above the classification rate, and `take_all` carries every
  *classified* unit into phase 2, which is `nu = resp_rate` rather than 1.
  Stratum labels, when supplied, must be unique and non-missing, since they
  identify strata in warnings and round trips. `$operational` carries the whole-unit field design,
  bounded by the pool phase 1 supplies, so a stratum whose expected
  phase-1 yield rounds below one unit is named in a warning rather than
  quietly making the reported `cv` infinite. `assurance` converts
  required respondents into an issued count that clears a stated
  probability from the binomial distribution rather than planning at the
  expectation, where about half of designs fall short. That level is
  marginal, holding stratum by stratum, so the chance every stratum clears
  at once is the product over strata and is lower; it is also conditional
  on the phase-1 pool, and a stratum whose assured issue exceeds what
  phase 1 supplies is named in a warning. `prec_twophase()` carries the
  design decisions a stratum table cannot express, so
  `n_twophase(prec_twophase(fit))` rebuilds the same problem rather than a
  new unconstrained one: the `take_all` set as supplied, a fixed
  `n_phase1`, `assurance`, and the comparator settings all survive the
  round trip. A fully pinned frame can be priced by `prec_twophase()` and
  is refused by `n_twophase()`, which would have nothing left to choose.
* `n_mean()`: sample size for a mean (`moe`, `rmoe`, and `cv` modes).
* `rmoe` states a precision target as a margin of error relative to the
  estimand, `moe / p` for a proportion and `moe / abs(mu)` for a mean,
  which is how MICS and DHS state one. It is an argument on `n_prop()`,
  `n_mean()` and `prec_prop()`, an indicator column on `n_multi()` and
  `n_multi_cluster()`, and a target column on the generalized `n_alloc()`.
  Supplying it replaces the arithmetic a plan used to carry: `rmoe = 0.12`
  says what `moe = 0.12 * p` said. It also completes the package's four
  precision quantities on two axes: `se` and `cv` describe the sampling
  variance and are the same under every interval method, while `moe` and
  `rmoe` describe the interval the chosen `prop_method` builds and are
  specific to it. `rmoe = q * cv`, with `q` the interval quantile the design
  uses, holds under `"wald"` alone; the other three miss it by 4 to 10
  percent at `p = 0.02`, `n = 900`, `deff = 2`, which is why a relative
  margin of error is a target in its own right rather than a rescaled `cv`. Every result reports `rmoe` wherever it reports
  `moe`, including the `.rmoe` stratum, domain, indicator, and constraint
  columns, so a target stated relatively reads back relatively. It is `NA`
  where the estimand has no known scale, as `cv` already is. `rme` and
  `RMoE` are recognized as near misses and point at `rmoe`.
* A coefficient of variation is a magnitude, so a mean may be negative
  wherever one is planned on: `n_mean()`, `prec_mean()`, the `mu` column of
  an indicator table, and the `mean` column of an allocation or two-phase
  frame all take the sign as immaterial and use `abs(mu)`. Zero is the
  value with no relative scale, and it is the one rejected.
* Indicator tables take the dispersion of a continuous indicator as `var`
  or as `sd`, matching the scalar methods, so a table written for one can
  be passed to the other. They are different quantities rather than
  synonyms, so supplying both is an error.
* `n_cluster()`: optimal multistage cluster allocation (2- and 3-stage,
  budget and cv modes).
* A cluster design can lose units at more than one stage, and the losses are
  not interchangeable, so each is named for the stage it acts on:
  `resp_rate_psu` for whole clusters that cannot be worked at all,
  `resp_rate_ssu` for second-stage units in a three-stage design, and bare
  `resp_rate` for the ultimate units. Losing a PSU removes its whole
  contribution and so deflates every variance term alike, which leaves the
  cost-optimal take untouched; losing ultimate units shrinks the realized
  cluster and so changes the clustering penalty too, which moves it. The
  optimal take is `sqrt(cost_psu * (1 - icc) / (cost_ssu * icc * resp_rate))`,
  and lumping the later rates into `resp_rate_psu` is conservative but wrong,
  by a margin that grows with the between-PSU homogeneity. Sizes and costs
  stay gross throughout: `$n` counts the units to issue and `$cost` pays for
  them, while the variance reads what they realize. Naming a stage the design
  does not have is an error rather than a silently ignored argument, because
  a dropped response rate plans a design with none.
* Every branch that chooses a clustered design reads the responding take, not
  only the branch that reports its variance. That covers the two-stage
  cost optimum in both cv and budget modes, the budget evaluators at two and
  three stages, the per-domain payload of a joint budget allocation, and the
  whole-cluster operational search in the classic stratified allocator. The
  invariant that pins them together is that one indicator is one design:
  `n_multi_cluster()` returns what `n_cluster()` returns for the same
  problem, at every response rate and in every mode. Before, a low ultimate
  response rate left the optimal take at its full-response value, a
  budget-mode design reported a CV its own precision function did not
  reproduce, a joint-domain budget allocation ignored the rate entirely, and
  the operational search could reject whole-cluster designs that do reach the
  target.
* A fixed stage size changes what the remaining stages should be, and every
  planner now solves the restricted problem rather than reusing the
  unrestricted answer. With `n_per_psu` fixed, `n_cluster()`'s cost-optimal
  final take is `sqrt(g (cost_psu + cost_ssu n2) / (a cost_tsu n2))` with
  `a = var_ratio_psu icc_psu + var_ratio_ssu icc_ssu / n2` and
  `g = var_ratio_ssu (1 - icc_ssu) / n2`, not the unrestricted
  `sqrt((1 - icc_ssu) / icc_ssu * cost_ssu / cost_tsu)`, which ignores both
  the fixed take and the per-PSU cost. The target CV cancels out of that
  ratio, so the same take minimizes cost at a fixed CV and CV at a fixed
  budget. The three-stage `n_multi_cluster()` fixed-stage inversions are read
  off the same response-aware representation the whole-unit search uses, so
  every allowed combination of fixed stages meets an attainable target
  exactly instead of approaching it.
* A fixed `n_psu` is solved for rather than ignored in the two-stage
  `n_multi_cluster()` continuous solution: the take becomes the only free
  stage and is inverted out of the precision requirement. A PSU count below
  the between-PSU floor is refused with the floor named, since no take
  however large gets under it.
* `n_multi()`: multi-indicator sample size under a simple design, with
  optional per-domain sizing via `domains` and a `min_n_domain` floor.
* `n_multi_cluster()`: explicit two- or three-stage multi-indicator cluster
  allocation. Multistage designs have their own function and result class,
  so cluster arguments never change what `n_multi()` returns.
* `n_alloc()`: stratified sample allocation given a frame with stratum
  sizes and variabilities. Three solve modes: fixed total `n`, target `cv`,
  or `budget` constraint. Four allocation methods: proportional, Neyman,
  optimal (cost-weighted), and Bankier power allocation. An `icc_psu`
  frame column switches to a stratified two-stage design (PSUs then
  elements per stratum) with cost-optimal or fixed per-stratum takes.
  All solve modes, methods, and constraints apply unchanged.
  `deff` and `resp_rate` take one value per stratum, as vectors or as frame
  columns, for designs whose fieldwork is not uniform. The three
  variance-based methods weight by `sqrt(deff_h / resp_rate_h)`, which is
  the optimum when the constraint is on the units drawn rather than the
  units that respond; proportional allocation weights by `1 / resp_rate_h`
  and takes no design effect, being a count rule whose purpose is a
  self-weighting responding sample. A value shared by every stratum is a
  constant factor that cancels out of a proportional weighting, so scalar
  arguments give exactly the classical rules. `max_weight` continues to
  bound the gross weight, which is no longer the response-adjusted one.
  Supplying long `measures` and `targets` tables requests joint constrained
  allocation across multiple indicators and overlapping domains, with CV and
  MOE requirements, costs, response/design effects, and stratum bounds. This
  includes one-stage designs and fixed-take two- or three-stage designs in
  which only PSU counts are optimized. Multistage output retains ultimate-unit
  `n` while exposing continuous and whole `n_psu` plus fixed later-stage takes.
  The result separates a KKT-certified continuous minimum-cost optimum from a
  deterministic feasible, locally cleaned integer operational recommendation;
  no claim of global integer optimality is made.
  Joint allocation answers either of the two questions planners ask. Without
  `objective` it returns the cheapest design meeting every requirement in
  `targets`. With `objective` and `budget` it returns the best design that
  budget can buy: among the allocations meeting every hard target and costing
  no more than `budget`, the one minimizing the priority-weighted sum of the
  named estimates' relative variances. Objective components are given either
  as indicator names or as a table of `name`, `domain`, `level`, and
  `priority`, so priorities weight variances rather than CVs and are scale
  free. The solver appends the objective as one more reciprocal constraint
  and root searches the objective bound whose minimum cost equals the budget,
  which keeps the same convexity and KKT certification as the minimum-cost
  mode. Results carry `$objective`, `$objective_value`, the budget residual,
  whether the budget binds, and the local change in the objective per extra
  unit of budget. Infeasibility is reported as one of three distinct cases:
  targets unattainable at the stratum upper bounds, targets attainable but
  unaffordable (with the cheapest feasible cost and the shortfall), or no
  whole-unit allocation meeting every target inside the budget.

* `n_pooled()`: size per occasion for an estimate pooled across occasions,
  the equal-weight mean of the occasion estimates, which is what an annual
  average built from quarterly rounds is. Three things share the word and
  this is none of the other two: it is not a pooled variance, and it is not
  the pooling of the rotating cohorts that make up one occasion, which
  `n_panel()` describes. `overlap` and `overlap_cor` are stated per lag
  rather than once, since a panel's correlation falls away with distance and
  a single number is rarely right over a long horizon; `cor_decay` gives the
  AR(1) shape `rho^m` in one argument. Targets are `moe`, `cv`, or `rmoe`,
  and the inversion is closed form, the finite population terms collecting
  into a constant no size can move. That constant is never positive under
  the nonnegative correlation contract, so no pooled target is out of reach
  for want of precision however far the occasions overlap: overlap raises
  the size a target needs and does not put a floor under the precision.

* `n_change()`: size per occasion for a change in a mean or a proportion
  between two occasions of one population, given how far the two samples
  overlap. `overlap` is the share of the first occasion carried into the
  second and `overlap_cor` the correlation among those shared units; only
  their product buys precision, so a full panel of uncorrelated measurements
  saves nothing. Supply `var`/`sd` for the mean scale or `p = c(p1, p2)` for
  the proportion scale, where the occasion variances and the change are both
  determined. The proportion scale carries the same `N / (N - 1)` adjustment
  `n_prop()` applies, because the variance the change formula wants is the
  population variance on `N - 1` degrees of freedom while a Bernoulli
  population's is `N p (1 - p) / (N - 1)`; without it one occasion of a
  change would disagree with `n_prop()` on the same finite population.
  Targets are `moe`, `cv`, or `rmoe`, and `ratio` sizes a large
  baseline against a smaller follow-up, which is where `overlap` becomes
  directional: it is measured against the first occasion, so it cannot
  exceed `1 / ratio`. The change variance is `A/n + B` in the net
  second-occasion size, with the finite population terms collected into a
  constant no sample size can move, so the target inverts in closed form
  rather than by search. At `overlap = 0` the size is exactly twice the
  single-occasion requirement, the two-independent-samples result. Full
  overlap with unit correlation at equal sizes leaves the change with no
  sampling variance at any size, and is reported as such rather than
  returning zero. On the proportion scale `overlap_cor` is checked against
  the largest correlation two Bernoulli marginals admit, above which no
  joint distribution exists. `predict()` sweeps the design: `overlap`,
  `overlap_cor`, `ratio`, `deff`, `resp_rate` and the target, reporting one
  size per occasion as `n1` and `n2`. At a finite `N` with equal occasions and full response the
  size rises towards `N` but never past it, because two censuses of one
  population measure the change exactly.
* `n_panel()`: how many units to recruit so that a panel still delivers a
  required responding sample after several waves of attrition. It embeds an
  existing `n_mean()`, `n_prop()`, `prec_mean()` or `prec_prop()` result
  rather than restating its arguments, reading it as both the responding
  count to reach and the estimand whose precision is reported at every
  wave; the target's own `resp_rate` is removed first, exactly, since the
  `n_*` functions apply it last. `resp_rate` is response at recruitment and
  `retention[j]` is conditional retention from wave `j` to `j + 1`, kept
  separate because a panel's loss concentrates at the first wave, 61
  percent of it on the UK LFS rates the examples use, and an average rate
  spread over the waves under-issues. **The two designs return different
  quantities and the object never writes one over the other**: a fixed
  panel reports `n_issued`, the whole issue to one cohort, and a rotating
  one reports `n_entrants` per occasion together with the separate
  `n_in_sample` its live cohorts hold between them, which differ by a factor
  of five on those rates. `n_in_sample` is also the cumulative recruitment
  reaching a steady state takes, whether the cohorts are taken on at once or
  phased in. A rotating occasion pools every cohort alive, so its
  responding sample is `e * sum(q)` and not one cohort's count at its final
  wave. `target_wave` sizes for a precision reached at an earlier wave,
  with the later waves still reported. `assurance` converts the
  expected-value recruitment, which leaves about half of all panels short,
  into one that clears the target with a stated probability: one binomial
  tail for a fixed panel, and for a rotating one the exact
  Poisson-binomial over cohorts at different cumulative probabilities,
  computed by convolution rather than approximated by a binomial on their
  mean. Both searches return the **smallest** sufficient recruitment,
  bracketing and bisecting rather than walking up from the expected count:
  that starting point is above the answer whenever the level is below a half,
  and it overshoots at a high response rate, `need = 2` at `r = 0.99`
  clearing 0.8 assurance with 2 issued where the expected count rounds to 3.
  The fixed search is the one `n_twophase()` uses, so its assured issue is
  now minimal in those two regimes as well. Where the assured recruitment
  exceeds a finite `N` the level is unattainable, not merely large: the
  recruitment it would take is still reported, `assured_feasible` records the
  state, and the expected design is untouched. Feasibility is reported
  against the design a planner would field, not only the continuous one: a
  whole-unit design that leaves a finite population is named even where the
  continuous plan fits inside it, which
  equal cohorts rounded up can do by as much as one unit per cohort, and the
  assured recruitment is checked in its own right, since it can exceed the
  expected one and above a level of a half usually does. Waves whose
  expected respondents fall below one are named too, their standard errors
  otherwise reading like any other wave's. A proportion target adds
  `expected_cases` to the wave table, which is where a `min_cases` floor can
  be watched stopping short after the target wave.
* `start` reports what a rotating design delivers while it is being brought
  up to its steady state, in `$launch` by occasion and `$launch_waves` by
  occasion and wave. `"gradual"` recruits one cohort an occasion, so the
  sample climbs over a full life; `"immediate"` splits the first occasion
  into one cohort per stage of the life and holds the whole sample at once.
  The membership overlap `design_overlap()` reports is the same either way,
  being a sum over life stages that entry dates do not enter, so this is a
  question about response and precision rather than about who is interviewed
  with whom. That is the reason it is worth reporting. An immediate launch is
  in membership equilibrium at its first occasion and **not** in response
  equilibrium, every unit there being at wave 1, so the occasion holds at
  least as many respondents as the design ever holds again, and strictly more
  as soon as any wave retains less than all of the one before, which makes it
  temporarily more precise. The two coincide exactly when every retention is
  1, whatever the response rate is. On the LFS rates its first
  occasion holds 1172 respondents against the design's 1001, `moe` 0.0287
  against 0.0310, converging down as the age mix matures; a gradual launch
  approaches the same figure from below, at 234. Either way the early
  occasions rest on a different response composition from the rest of the
  series, which is what `launch_waves` exposes and what nonresponse weighting
  has to carry. The default plans no launch, leaving every other result
  exactly as it was, and a fixed panel refuses the argument, recruiting one
  cohort that is its whole design from wave 1. The tables are continuous, as
  the wave table is, so the last occasion lands on `n_target` rather than on
  a rounded count; whole units are printed. Described for a life without a
  break in it: a schedule that leaves the sample and returns needs launch
  cohorts selected before they are first interviewed, which is a longer
  definition than this argument carries.

## Precision analysis

* `prec_prop()`: sampling precision (se, moe, cv) for a proportion given a
  sample size. Precision counterpart to `n_prop()`.
* `prec_mean()`: sampling precision for a mean given a sample size. Precision
  counterpart to `n_mean()`.
* `prec_prop()` and `prec_mean()` also solve in the remaining direction.
  The proportion, the sample size, and the precision are three readings of
  one equation: `n_prop()` fixes the first two and returns the third,
  `prec_prop()` fixes the first and the last, and supplying `cv` in place of
  `p` fixes the last two and returns the smallest proportion the design
  measures that precisely. Because `se(p) / p` falls as `p` rises, that
  answer is a floor, and every larger proportion meets the target, which is
  the planning form of a publication threshold. `expected_cases` is reported
  beside it, so the precision floor and the `min_cases` floor sit on one
  footing. Because `cv` reads the sampling variance rather than an interval,
  that solution is closed form and the same under all four methods.
  `solved` records the direction and `predict()` varies `cv` in place of the
  level it solved for. `prec_mean()` takes `cv` the same way, returning the
  magnitude, the standard error of a mean carrying no mean of its own.
* `prec_prop()` asks the same question of the interval when given `rmoe`
  instead of `cv`: the smallest proportion the design can report at that
  relative half-width under the chosen `prop_method`. This one is
  method-specific and solved numerically. Two properties the `cv` solve
  does not share bound it, and both are reported rather than silently
  approximated: the Wilson and Korn-Graubard half-widths do not vanish as
  `p` approaches 1, so their relative margin of error has a positive floor
  and a target below it is refused with the floor named; and the
  back-transformed log-odds half-width turns upward once the logit spread
  outgrows `logit(p)`, near `p = 0.999` at `n = 1500` but as low as
  `p = 0.97` at `n = 30`, so the search brackets at that turn and returns
  its lower root.
* `prec_cluster()`: sampling precision (cv) for a multistage cluster
  allocation. Precision counterpart to `n_cluster()`.
* `prec_multi()`: per-indicator sampling precision for multi-indicator
  survey designs. Precision counterpart to `n_multi()`.
* `prec_multi_cluster()`: per-indicator precision for multistage
  multi-indicator designs. Stage costs are optional round-trip metadata and
  do not enter the precision calculation.
* `prec_alloc()`: sampling precision for a stratified allocation. Precision
  counterpart to `n_alloc()`. It reports at three levels: the population in
  `se`, `moe` and `cv`; each stratum in `$detail`, as `.se`, `.moe`,
  `.rmoe`, `.cv` and the variance `.share`; and each domain in `$domains`, the same table
  `n_alloc()` returns, so a design and its assessment compare row for row.
  `.se` is the stratum's own mean, on the scale of `mean` or `p`, and
  `.share` is the stratum's contribution to the design variance of the
  overall mean, summing to 1, which locates where precision is being bought.
  `n_alloc()` carries the same stratum columns.
  Joint allocation results retain one row per
  target and report stratum-bound violations when assessing an adopted sample.
  Joint assessment takes the same `min_n_stratum` floor as `n_alloc()`, so a design
  and a direct assessment of it report against identical bounds.
  When the fitted design carries an objective, the assessment also reports its
  components and weighted value. Minimum-cost designs invert through
  precision, pinning the achieved values as the requirement; budget-objective
  designs invert through cost, because pinning achieved precision as hard
  targets would over-constrain a design that already spends its whole budget.

* `prec_pooled()`: precision of an estimate pooled across occasions given
  the size of each, the counterpart to `n_pooled()`. It reads the same
  covariance kernel `prec_change()` does, so the two are contrasts on one
  object rather than two models: a change is a difference of two occasion
  means and a pooled estimate is their average. Two boundaries fix it. A
  fresh sample each occasion pools to a single occasion's variance divided
  by the number of them, independent averaging; a full panel measured with
  correlation 1 pools to a single occasion's variance unchanged, the same
  units every time and averaging buying nothing. A lag whose overlap is zero
  contributes exactly zero covariance, population term included, which is
  `prec_change()`'s rule at `overlap = 0` read lag by lag and the only
  reading under which the first boundary is exact. A covariance is added in
  a sum and subtracted in a difference, so whichever sign it carries it
  moves a pooled estimate and a change in opposite directions; when it is
  positive, the ordinary case, overlap works against the pooled estimate and
  for the change. It is positive exactly when the overlap exceeds the
  sampling fraction `n / N`, which every positive overlap does when there is
  no finite population correction. Below that fraction the occasions are
  negatively coordinated, sharing fewer units than chance would give them,
  and the directions reverse. The level at a single occasion is unaffected
  either way, so what trades off is the change against the pooled estimate
  and the level is a reference line rather than a third position.
  A `design_overlap()` result is accepted directly at `resp_rate = 1` and
  refused below it, since that function reports the overlap between issued
  samples while this variance is formed on respondents. The assembled
  covariance is checked for validity on the matrix rather than on the
  variance read off it, because a lag's covariance changes sign once its
  overlap falls below the sampling fraction and a kernel that fails the
  check can still return a positive pooled variance; the verdict is
  invariant to the units the outcome is measured in, and the refusal names
  the fraction, the lag, and whether the rotation would also exhaust its own
  population, counting interviews per cohort rather than calendar span so
  that a schedule with a rest period in it is read correctly. `confint()` centres the interval on the level.

* `prec_change()`: precision of a change between two occasions given the
  size of each, the counterpart to `n_change()`. The marginal terms carry
  their own finite population correction while the overlap covariance
  carries a single `1/N`, so at full overlap and unit correlation with equal
  sizes and variances the population terms cancel exactly and the whole
  reduces to `2 S^2 (1 - overlap) / n`. `deff` multiplies the assembled
  variance and is the design effect of the change, which is not in general
  either occasion's: a clustered design revisiting the same clusters has a
  smaller one than a design that reclusters. The covariance is a model form,
  assuming a correlation among the shared units and simple random sampling
  otherwise; it is not the design-based covariance of two waves of a complex
  design, and `overlap` is asserted rather than derived from a rotation
  pattern. The sizes enter the covariance after `resp_rate` has netted them
  down, so `overlap` is the overlap between the two responding samples
  rather than between the issued ones; the two coincide at `resp_rate = 1`,
  and converting between them otherwise needs an assumption about how
  response at the second occasion depends on response at the first, which
  the function does not make. The variance is piecewise in `overlap`: at
  zero the covariance is dropped entirely, population term included, which
  is what lets the two occasions come from different populations.
  `confint()` centres the interval on the change.
* `prec_panel()`: the responding sample and the precision a recruitment
  already budgeted leaves at each wave, the counterpart to `n_panel()`. It
  takes the recruitment directly, so it is an inverse rather than a
  formatter, and reports the requirement alongside it: a recruitment below
  what the target asks for is named as short rather than left to be
  compared. Each wave is the embedded estimand evaluated at that wave's
  expected respondents with `resp_rate = 1`, the panel's losses having
  already been applied, which is what makes
  `prec_panel(n_panel(target, retention, target_wave = w))` reproduce the
  target's own precision at wave `w`. Both directions return
  `svyplan_panel`, unlike the rest of the package, because a panel's
  recruitment and its precision move together across the waves and are only
  readable side by side; `$solved` records which direction produced the
  object. The `svyplan_panel` method re-reads a stored plan under named
  overrides, so retention that turned out worse than planned needs no
  second call to `n_panel()`. `design` is the one stored argument it will
  not take: the recruitment it holds is issued to one cohort under `"fixed"`
  and taken on every occasion under `"rotating"`, so reading it under the
  other design changes the unit rather than the assumption. Restate the
  count through the direct form to do that deliberately.

All `n_*` and `prec_*` functions are S3 generics with bidirectional
round-trip: passing a precision object to the corresponding `n_*` function
recovers the continuous sample-size target under the same method and design
assumptions, and vice versa. Method changes and operational integer rounding
can change the result. Round-trip methods accept named
`...` overrides of any stored argument (e.g. `prec_prop(x, deff = 2)`).
A `NULL` value unsets a stored argument, and unknown names raise an error
instead of being silently ignored. The three targets `moe`, `cv` and `rmoe`
are all reachable as overrides: the achieved margin of error is restored as
the implied target only when the caller names none of them, so
`n_mean(x, rmoe = 0.05)` solves for the relative margin instead of colliding
with the stored absolute one.

All public calculation, prediction, coercion, and display methods reject
unused arguments passed through `...`. This catches misspelled names instead
of silently computing a result with an unintended default. Legitimate plan,
round-trip, plotting, and standard data-frame arguments remain supported.

Optional arguments in the public calculation methods follow `...`, so
they must be fully named. This makes calls resilient to later additions and
prevents positional or partial matching of design controls.

`power_mean()` takes the always-required `var` as its primary argument and
keeps `effect` optional for minimum-detectable-effect calculations.
`n_cluster()` displays `stage_cost = NULL` because that value may come from a
plan. `strata_bound()` requires `n_strata`, which determines the shape of
its result and, like every other function, accepts a `plan`. Short method and allocation choices are displayed in function
signatures where applicable.

## Power analysis

* `power_prop()`: power analysis for two-sample proportion tests. Solves
  for sample size, power, or minimum detectable effect (MDE). Supports
  panel overlap for repeated surveys. MDE mode searches both directions
  (`p2 > p1` and `p2 < p1`) and returns the closest detectable alternative.
  The two-sided/one-sided switch is `alternative`, matching base R. Supports unequal group
  sizes (`n = c(n1, n2)`), allocation ratio (`ratio`), and arcsine and
  log-odds transform methods (Valliant, 2018, Sections 4.3.4--4.3.5).
* `power_mean()`: power analysis for two-sample mean tests. Same solve
  modes and features as `power_prop()`, including `alternative`.
  Supports unequal group variances (`var = c(v1, v2)`), unequal group
  sizes (`n = c(n1, n2)`), and allocation ratio (`ratio`). Cohen's d
  conversion documented.
* Panel overlap under a finite `N` corrects the marginal terms but not the
  overlap covariance, which carries a single `1/N`: two samples sharing
  `k = overlap * n1` units have
  `Cov = overlap_cor * S1 * S2 * (k/(n1 n2) - 1/N)`. At `overlap_cor = 1`
  with equal sizes and variances the population terms cancel and the
  difference variance is `2 * var * (1 - overlap) / n`, free of `N`. A
  positive `overlap` describes a coordinated design, so it requires both
  occasions to sample one population; the default 0 is the ordinary
  two-group comparison, where the groups are independent and may be
  different populations of different sizes. In `power_did()` the overlap
  is within each arm, so each arm's covariance uses that arm's own `N`.
* `power_did()`: power analysis for difference-in-differences designs.
  Parametrized via `treat = c(baseline, endline)` and
  `control = c(baseline, endline)` vectors. Supports both proportion and
  mean outcomes, cell-specific variances, panel overlap, and all common
  design parameters. Those two paths already determine the contrast, so
  `effect` is optional and defaults to the difference-in-differences they
  imply; supplying a value that disagrees with them warns rather than
  silently planning for an effect the displayed inputs do not produce.
  Leaving both `n` and `power` supplied while `effect` is `NULL` still
  solves for the minimum detectable effect.

## Stratification

* `strata_bound()`: candidate strata boundary construction for a continuous
  stratification variable. Four methods: Dalenius-Hodges cumulative root
  frequency (`"cumrootf"`), geometric progression (`"geo"`),
  LH-inspired coordinate optimization (`"lh"`), and Kozak-inspired
  random-restart local search
  (`"kozak"`). Four allocation methods: proportional, Neyman, optimal
  (cost-weighted), and Bankier (1988) power allocation (`"power"`) with
  parameter `q` controlling the national/subnational precision trade-off.
  Take-all (certainty) strata via the `take_all` argument.
  `deff` and `resp_rate` place the reported `n` and `cv` on the same scale
  as every other function, so `cv = 0.05` means one thing across the package
  and `n` is always a fielded sample. Both default to 1, the identity. A
  scalar `deff` does not move the boundaries, since it scales the variance of
  every candidate set equally; it changes the `n` a `cv` target requires and
  the `cv` a given `n` achieves.
  `$strata` carries a `mean` column, so the table has the `N`, `sd`, and
  `mean` columns [`n_alloc()`] expects and can be handed straight to it. With
  matching `cv`, `deff`, and `resp_rate` the two functions agree on the
  continuous total.
  `strata_bound()` accepts `plan` like every other function. A profile can
  supply `alloc`, `alloc_q`, `deff`, `resp_rate`, and a scalar `unit_cost`;
  a profile carrying a vector `unit_cost` is rejected, because
  `strata_bound()` orders costs from the lowest to the highest stratum while
  `n_alloc()` orders them by frame row, so a vector written for one would be
  silently misapplied by the other.
* `predict.svyplan_strata()`: apply strata boundaries to new data,
  returning a factor.
* Strata results retain full-precision `share` and `sd` values, validate
  count inputs as whole numbers, and reject nonfinite stratification values,
  thresholds, and costs. `as.double()` returns the total sample size,
  consistently with `as.integer()`. Cutpoints remain in `$boundaries`.

Cluster stage tables returned by `as.data.frame()` take `n_int` from the
constraint-preserving operational design rather than rounding each
continuous stage size upward independently.

## Design components

* `design_overlap()`: the sample overlap a rotation schedule produces, at
  every lag it reaches. This is the overlap of the **issued** sample, a
  property of the schedule alone, and it is the `overlap` that `n_change()`,
  `prec_change()` and the power family take **at full response**, replacing
  a number the planner would otherwise have to already know. Those
  functions read `overlap` among respondents, netting each occasion by
  `resp_rate` before the covariance forms, so below full response the two
  differ and converting between them needs an assumption the package does
  not make: under response independent between occasions at rate `r`, an
  issued overlap `f` leaves about `f * r` among respondents, rising towards
  `f` as response persists. A schedule is given as a
  compact spec in either notation the literature uses, told apart by whether
  it contains a `0`. Without one it names the spells in order from in
  sample, `"4-8-4"` for the CPS pattern or `"5"` for an unbroken five-wave
  panel; with one it is a flag per occasion, `"1-1-0-0-1-1"` for the
  two-in-two-out-two-in design, which is the same life as `"2-2-2"`. It may
  also be given as one entry per occasion of a unit's life, where unequal
  positive entries describe a design that subsamples later waves. A string
  of all 1s is the one form both notations claim, and it is refused rather
  than resolved by precedence, since `"1-1-1"` is three consecutive occasions
  under one reading and in-out-in under the other, whose consecutive
  overlaps are 2/3 and 0, and the message names both readings and the
  unambiguous spelling of each, `"3"` and `"1-0-1"`. One cohort enters each occasion, so
  at a steady state the sample is `sum(w)` and the units in sample at both
  `t` and `t + m` are `sum(min(w[s], w[s + m]))`, the `min` assuming the
  smaller take at the two occasions is nested in the larger, whichever
  occasion that is. The result is a numeric vector indexed by
  lag, so `design_overlap("4-8-4")[1]` is the consecutive overlap and
  `[12]` the overlap a year apart on monthly occasions, both plain numbers
  ready to pass on. **A schedule is a finite life, not a repeating cycle**,
  and the difference is not visible in the answer: read as a cycle, only
  the cohort finishing its stint leaves each occasion and the CPS
  consecutive overlap comes out at 87.5%, while the finite life loses a
  cohort at the end of each of its two spells and gives the published 75%.
  A compact spec with an even number of spells could only mean a cycle and
  is refused, and so is the profile passed whole where a single `overlap` is
  wanted: it is numeric and inside `[0, 1]`, so the generic message would be
  true and useless, and the refusal names the lag to extract instead. So are a life shorter than two occasions, which has no lag at
  which any sample is shared, and a take larger than the first, a wave being
  able only to re-interview part of the cohort it recruited. Where the takes
  at two occasions are not nested, `min` is the largest overlap they admit
  rather than the overlap the design has, nesting being an assumption about
  the design and not something the takes identify. Overlap is one number here rather than two directional ones
  because a stationary rotation holds the same mix of life stages at every
  occasion, so the occasion sizes are equal and the two directions
  coincide. The correlation between occasions is not a property of the
  schedule and still has to come from a previous round.
* `plot()` on a `design_overlap()` result draws the rotation chart the
  design would be published as: one row per cohort, one column per time
  period, and a cell labelled by wave wherever that cohort is in sample. It
  is the fastest way to see that a schedule is the one that was meant, which
  matters for a notation whose readings differ silently. A cohort enters at
  every period drawn, so the total row climbs through the launch and settles
  at the `n_occasion` the overlaps divide by, with the period it settles at
  marked. The chart's right-hand side is the steady state the overlaps
  describe, and its left-hand side is the gradual start a rotating design
  has unless the first period's sample is split into cohorts of unequal
  life, which one schedule does not express. A take that varies over the
  life shades its cell in proportion. An `n_panel()` result with
  `design = "rotating"` labels the rows and the total row in units instead
  of cohort shares, issuing whole units as `print()` does so that one object
  reports one recruitment. `type = "overlap"` draws the lag profile instead,
  which is the table `print()` reports.
* Three-stage designs derive `var_ratio_ssu` instead of defaulting it to 1. In the
  multiplier
  `D = var_ratio_psu icc_psu m q + var_ratio_ssu (1 + icc_ssu (q - 1))`, `var_ratio_psu`
  rescales the components' unit variance to the analysis variable and
  `var_ratio_ssu` does the same for the within-PSU part, which is `1 - icc_psu`
  of it, so the two are linked by `var_ratio_ssu = var_ratio_psu * (1 - icc_psu)`. That
  identity is what makes `D` collapse to `var_ratio_psu` at `m = q = 1`, where one
  unit is taken per SSU and one SSU per PSU and no clustering remains.
  Defaulting `var_ratio_ssu` to 1 asserted instead that the within-PSU variance was
  the whole variance, contradicting any positive `icc_psu` in the same
  expression and inflating the design effect by
  `var_ratio_psu icc_psu (1 + icc_ssu (q - 1))`. This is a breaking change:
  three-stage sample sizes fall, by about 1.4 to 4.4 percent for
  `icc_psu` between 0.02 and 0.10, and by more when the takes are small.
  It affects `design_effect()`, `n_cluster()`, `prec_cluster()`,
  `n_multi_cluster()`, `prec_multi_cluster()`, and the fixed-take
  three-stage mode of `n_alloc()` and `prec_alloc()`. A scalar `var_ratio` now names
  `var_ratio_psu` and derives `var_ratio_ssu`; supplying both explicitly still overrides the
  identity. Two-stage designs are unaffected, since `var_ratio_psu` there is already
  1 by the same argument. `varcomp()` estimates the two `var_ratio` values from
  their own stage decompositions rather than imposing the identity, so its
  pair can differ by several percent on small clusters; passing a
  `svyplan_varcomp` uses the estimated values as given.
  Written on components referenced to the total unit variance the same
  quantity is the familiar `1 + icc_1 (m q - 1) + icc_2 (q - 1)`; the
  package's `icc_ssu` is referenced to the within-PSU variance instead,
  so `icc_ssu = icc_2 / (1 - icc_1)`.

* `varcomp()`: variance component estimation from frame data via nested
  ANOVA (SRS and PPS). S3 generic with methods for formulas, numeric
  vectors, and `survey::svydesign` objects. The `svydesign` method
  treats design weights as inverse inclusion probabilities: cluster
  sizes and totals are estimated by summed weights and the estimation
  variance of the weighted totals is subtracted from the between-stage
  terms, so unequal-probability samples from a previous round give
  approximately design-unbiased components (unit weights recover the
  frame formulas exactly). The formula and vector interfaces accept the
  same weights directly via a `weights` argument, so sample-based
  estimation does not require constructing a `svydesign` object. The
  documentation specifies the within-cluster weight scale and the
  renormalization of `prob` over sampled PSUs. The formula interface
  warns when `data` is a samplyr sample and no weights are supplied.
  A `strata` argument estimates components per stratum, returning a
  table whose columns match the `n_alloc()` frame contract.
  `icc` and `var_ratio` are ratios of components sharing a denominator,
  so neither depends on the outcome mean, and both are reported for an
  outcome centered on zero even though every relvariance there is infinite.
  That case carries its own warning and is kept apart from a constant
  outcome, which is the one condition that genuinely leaves the components
  unidentified and `icc` at 0 by convention.
* `varcomp(deff = , n_per_psu = )` backs an `icc` out of a previous round's
  published design effect, the most widely available source of one, by
  inverting the same two-stage identity `design_effect()` applies forward.
  `se` with `p` or `var` and `n` starts from a published standard error
  instead. A scalar design effect fixes one product, so `varb`, `varw` and
  `unit_relvar` come back `NA_real_` and `$source` records where the `icc`
  came from; `n_cluster()` and `prec_cluster()` ask for `unit_relvar`
  directly rather than pricing a CV off the default. The take that
  identifies the `icc` is the size-weighted average
  `sum(b_i^2) / sum(b_i) = b_bar (1 + cv_b^2)`, which a planned constant
  take collapses to but a realized one does not, so `n_per_psu` also
  accepts the realized takes or a nominal take with `cv_take`. It is stored
  as provenance and `design_effect()` does not read it back: re-planning at
  a different take is the reason to back an `icc` out, and defaulting to
  the source take would return the previous design's number. A design
  effect below `var_ratio` yields a negative `icc`, which is returned with
  a warning naming where it will be refused rather than clamped.
* `design_effect()`: S3 generic building the design effect a planned design
  is expected to produce. Components are selected by the arguments supplied
  and multiply together: clustering from `icc` with `n_per_psu` (and
  `n_per_ssu`), unequal weighting from planned `weights` or from a `strata`
  table's `N` and `n`, and the stratification gain from that table's `sd`
  and `mean`. The clustering component uses the same variance model as
  `n_cluster()` and `prec_cluster()`, including three stages and non-unit
  `var_ratio`, so `design_effect()` and the cluster functions always agree.
  Methods read the design features directly off a `n_cluster()` or
  `prec_cluster()` allocation, a `varcomp()` estimate, or an `n_alloc()`
  allocation. The result is a `svyplan_deff` object, usable as a plain
  number wherever `deff` is expected, whose print and `as.data.frame()`
  show the component decomposition.

  Given an `n_alloc()` result, `design_effect()` returns the allocation's
  own variance ratio rather than a product of approximations: the
  per-stratum cluster factors enter stratum by stratum, and the scalar
  `deff` the allocation was built under is counted once. Multiplying
  components supplied by hand is still Kish's approximation, exact when the
  stratum standard deviations agree and not otherwise, and the gap runs in
  *either* direction, so it is not a bound: a Neyman allocation over very
  unequal `sd` reads far above its true ratio, while a constraint forcing
  units into a low-variance stratum reads below it. A printed result with
  more than one component is marked `approx. (Kish)` for that reason. The
  ratio is a without-FPC planning quantity, so it need not match one
  derived from `prec_alloc()`, which applies each stratum's FPC. A
  generalized allocation optimizing several measures has no single ratio
  and is refused by both `design_effect()` and `effective_n()`. Supplying
  `weights` for an allocation charges an *additional* anticipated
  adjustment on top of it.

  `design_effect()` is a planning tool only: it anticipates a design effect
  from design parameters, and does not estimate one from collected data.
  Estimating a realized design effect is
  `survey::svymean(..., deff = TRUE)` and `survey::deff()`, which use the
  actual weights, strata, and clusters. Passing a bare numeric vector
  raises an error pointing at both.
* `effective_n()`: mirror of `design_effect()`, returning
  `n * resp_rate / deff`. It accepts the same components, a ready-made
  `deff`, or a plan, and derives `n` from `weights`, `strata$n`, or the
  plan when it can. Sizes count units issued, so nonresponse comes off
  before the design effect is applied; this is the identity the allocation
  and precision functions plan on and the one reported in the `n_eff`
  column of an `n_alloc()` table, so `effective_n(plan)` and that column
  agree. A plan supplies its own `resp_rate`, and `resp_rate` can be given
  directly alongside a bare `n`.

## Survey plan profiles

* `svyplan()`: create a reusable profile capturing shared design defaults
  (`deff`, `N`, `resp_rate`, `alpha`, `stage_cost`, `unit_cost`, etc.). Pass
  to any function via `plan = plan` or pipe with `plan |> n_prop(...)`.
  Piping works with both positional and named arguments
  (e.g. `plan |> n_prop(p = 0.3, moe = 0.05)`).
  Explicit arguments always override plan defaults.
* `svyplan()` and `update.svyplan()` validate every supported default at
  construction time. Related cluster defaults such as `stage_cost`, `icc`,
  and `var_ratio` are also checked for compatible stage counts.

## Naming

* Argument names say what a quantity is rather than reproducing one
  textbook's symbol. `icc` is the design-based measure of homogeneity
  within clusters (written the Greek delta in Valliant, Dever and Kreuter),
  `var_ratio` the ratio of the stage components' unit variance to the
  analysis variable's, and `unit_relvar` the unit relvariance. See the
  Notation section of `?svyplan-package` for the full symbol map. `icc` is
  constrained to `[0, 1]` and is not interchangeable with a mixed-model
  ICC, which can be negative.
* Stage sample sizes distinguish counts from takes. `n_psu` is the number
  of PSUs selected, while `n_per_psu` and `n_per_ssu` are the units selected
  per PSU and per SSU. These are sample takes, not the population size of a
  cluster, which is what a name like `psu_size` would suggest to readers who
  size clusters from a listing.
* `overlap_cor` is the correlation between occasions in a repeated survey,
  named as the partner of `overlap`. It is not a measure of homogeneity.
* `alloc_q` is the Bankier power-allocation exponent, used only when
  `alloc = "power"`. It is unrelated to the statistical power computed by
  the `power_*()` functions.
* The multi-indicator functions take `indicators`, a wide table with one row
  per indicator fusing what to measure with how precisely. `n_alloc()` and
  `prec_alloc()` take `measures` plus `targets`, a long table splitting the
  same information. The two spellings name two different shapes, and the
  fitted objects carry them in separate fields.
* Wherever a mean is planned, dispersion may be given as either `var` or
  `sd`; supply exactly one. Stratum frames and published survey reports
  usually quote standard deviations.
* `n_multi_cluster(allocation = c("separate", "joint"))` names the two
  strategies in the signature rather than hiding them behind a flag.
* Cluster/multistage functions use `stage_cost` for per-stage cost vectors
  (`n_cluster()`, `n_multi()`, `prec_multi()`).
* Stratified allocation functions use `unit_cost` for per-stratum unit costs
  (`n_alloc()`, `prec_alloc()`, `strata_bound()`).
  The optional frame column and the per-stratum detail column carry that
  same name; a `cost` frame column is rejected with a message naming the
  replacement, while `$operational$cost` is the total field cost.
* The first argument of `prec_alloc()` is `frame`, matching `n_alloc()`.
* `strata_bound()` uses `n_class` and `max_iter` for its public controls.
* Take-all strata are named for what they hold. `strata_bound(take_all_above =)`
  is the numeric threshold, so the name says the argument is a cutoff and not
  a flag; the logical output column and the `n_alloc()` frame column stay
  `take_all`. Every path that reads a `take_all` column validates it the
  same way: logical, or numeric restricted to 0 and 1. Any other number is
  a mistake rather than a truthy value, so it is rejected instead of
  silently pinning a stratum.
* Sample-size floors are named for what they count. `min_n_domain`
  (`n_multi()`, `n_multi_cluster()`) is a per-domain floor and `min_n_stratum`
  (`n_alloc()`, `prec_alloc()`) a per-stratum one. `svyplan()` accepts both
  as defaults.

## Domain handling

* `n_multi()`, `prec_multi()`, `n_alloc()`, and `prec_alloc()` require an
  explicit `domains` parameter naming the domain columns. Columns not
  listed there are ignored, so no unrecognized column is ever treated as a
  domain variable by accident.
* `n_multi()` and `prec_multi()` results store `params$domain_cols`,
  `params$mode` (`"moe"`, `"cv"`, or `"budget"`), and `params$prop_method`
  for clean round-trip conversion. The round-trip methods
  (`prec_multi.svyplan_n`, `prec_multi.svyplan_cluster`,
  `n_multi.svyplan_prec`) read these fields directly instead of
  reverse-engineering domain columns from output tables.
* `svyplan()` does not accept `method` as a plan default (its meaning is
  ambiguous across function families). The `prop_method` default
  (validated at construction: `"wald"`, `"wilson"`, or `"logodds"`) covers
  `n_multi()`/`prec_multi()` and also fills the `method` argument of
  `n_prop()`, `prec_prop()`, and `power_prop()` when the value is valid
  for that function.

## Common features

* Result constructors return complete objects with canonical fields.
  Allocation results pass their continuous precision and operational design
  into the constructor rather than mutating the result afterward.
* `as.data.frame()` has explicit schemas for `svyplan_prec` and
  `svyplan_power` objects. Unstratified `svyplan_varcomp` results also export
  as one-row two- or three-stage component tables.
* All `as.data.frame()` methods accept the `validRN` argument forwarded by
  `data.frame()` in R 4.7.0 and later.
* `design_effect()` returns a numeric `svyplan_deff` object. `as.double()`
  extracts the overall design effect. The component decomposition is
  reached by name under one set of names shared by every access route,
  `deff` for the overall value and `deff_<component>` for each part, so
  `d$deff_cluster`, `d[["deff_cluster"]]`, `as.list(d)$deff_cluster`, and
  `as.data.frame(d)$deff_cluster` are the same number; `as.list()` and
  `as.data.frame()` return the whole decomposition as a named list and a
  one-row table. Naming a component the design effect does not have is an
  error listing the ones it does. The object is directly usable as the
  `deff` argument in sample-size, precision, and power calculations.
* All sample size, precision, and power functions accept `deff` (design
  effect), `N` (finite population correction), and `resp_rate` (response
  rate adjustment). One shared variance equation is used throughout:
  with `n_net = n * resp_rate` responding units, `deff` multiplies the
  SRS variance at `n_net` and the finite population correction uses the
  actual sampling fraction `n_net / N` (so a census has zero sampling
  variance under any design effect). Every proportion method reads that
  variance through the effective size
  `n_eff = n_net / (deff * fpc(n_net))`, so the Wald, Wilson, log-odds, and
  beta intervals all respond to `deff` and `N`, and `confint()` and the
  reported margin of error always describe the same interval. Targets that would require
  drawing more than `N` units from a finite frame raise an error instead
  of returning an impossible sample size, and the precision and power
  evaluators (`prec_prop()`, `prec_mean()`, `prec_multi()`, supplied-`n`
  `power_*()`) likewise reject a supplied gross `n` above the frame size
  of any group or indicator.
* All functions accept `plan`, a `svyplan()` profile providing shared
  design defaults.
* `predict()` methods for sensitivity analysis: evaluate any result at new
  parameter combinations. For a fixed-budget joint allocation this returns the
  cost-versus-objective frontier, answering what other budgets would buy;
  budgets that cannot fund the hard targets give an `NA` row flagged
  `.feasible = FALSE` rather than discarding the rest of the frontier.
* `plot()` methods for strata boundaries (per-stratum sampling fractions),
  power results (the power curve), and fixed-budget joint allocations (the
  budget frontier). The frontier grid starts at the cheapest design that
  meets the hard targets, so the drawn curve spans only fundable budgets.
* `confint()` methods for `svyplan_n` and `svyplan_prec` objects, documented
  at `?confint.svyplan`. The interval type follows the `method` the result
  was computed with, and `level` is independent of the `alpha` the design was
  sized at.
* When the survey package is installed, `survey::SE()` and `survey::cv()`
  methods are registered automatically.

## S3 classes

* `svyplan`: survey plan profile (reusable design defaults).
* `svyplan_n`: sample size results (with se, moe, cv fields).
* `svyplan_cluster`: multistage allocation results.
* `svyplan_prec`: precision results.
* `svyplan_varcomp`: variance component estimates.
* `svyplan_strata`: strata boundary results.
* `svyplan_power`: power analysis results.
* `svyplan_deff`: planning design effects, carrying their component
  decomposition.
* `svyplan_overlap`: rotation overlaps indexed by lag, carrying the shared
  counts and the schedule they came from. Subsetting a lag returns a bare
  number, so the lag never travels on as a name into the result of whatever
  it was passed to, and arithmetic, comparison and the `Math` group return
  bare numerics: a transformed overlap is no longer the profile the counts
  and the schedule describe, and keeping the class would leave the two
  contradicting each other. Replacement is refused for the same reason,
  which also settles `pmax()` and `pmin()`: they copy their first argument's
  attributes without dispatching, and assign through `[<-`, so neither can
  return a mislabelled overlap. `svyplan_deff` and `svyplan_df` refuse
  replacement too, a design effect being the product of the components it
  carries and a df the count derived from its strata.
* `svyplan_panel`: panel recruitment and the wave-by-wave sample it leaves,
  returned by both `n_panel()` and `prec_panel()`. It is a sibling of
  `svyplan_n` and not a subtype: its headline number is units to release,
  and the methods registered for an analysis sample would each answer a
  different question about it. `as.integer()` returns the whole units to
  recruit, never the assured count, so which number to field stays the
  planner's call.

All classes have print and format methods.

## Input validation

* `n_multi()` and `prec_multi()` reject non-positive `unit_relvar`,
  `var_ratio_psu`, and `var_ratio_ssu` values in multistage mode.
* `prec_multi()` multistage mode validates `icc_psu` (and `icc_ssu` for
  3-stage) presence, type, NA, and range.
* Arguments a configuration cannot use are rejected rather than silently
  discarded, extending to named arguments the contract the package already
  applies to `...`. `strata_bound()` errors when a control belongs to a
  different method (`n_class` is `"cumrootf"` only, `max_iter` is `"lh"` and
  `"kozak"`, `n_restart` is `"kozak"` only), and `power_did()` errors on
  `var` or `sd` under `outcome = "prop"`, where the cell variances follow
  from `treat` and `control`. `max_iter` now defaults to `NULL` and resolves
  to 200 for the methods that use it, so that supplying it is detectable.
* `prec_cluster()` validates that `unit_relvar` and `var_ratio` are positive and finite.
* `varcomp()` rejects NA and empty outcome vectors, and data with a
  single PSU (overall or within a stratum) with a clear error naming
  the stratum.
* `confint()` methods for `svyplan_n` and `svyplan_prec` validate that
  `level` is in (0, 1).
* `design_effect()` and `effective_n()` reject weights and stratum columns
  holding non-finite, missing, or non-positive values, and reject a stratum
  table that supplies the weighting component twice.

## Feasibility and integer designs

* Constrained designs separate the continuous mathematical optimum
  (top-level fields) from the whole-unit field design (`$operational`,
  with cost, cv, and se recomputed from the integer design). `print()`
  leads with the field design, and `as.integer()` returns the operational
  design in the same shape as `n` (stage vector for cluster plans,
  total for allocations) and `as.double()` its continuous counterpart.
  Joint allocations print the violated constraints, or the binding ones
  when everything passes, and say which subset is on screen out of how many
  so a short block does not read as the whole set.
* `n_cluster()` finds the operational design by discrete search
  (enumerated whole stage sizes): budget-mode field designs never
  exceed the budget and cv-mode field designs meet the target with
  whole units. The three-stage search enumerates PSU counts and SSU takes
  and solves the ultimate take in closed form, so the whole final-stage
  axis is covered. That axis needs it: the cheapest whole-unit design sits
  where the required PSU count is a tight ceiling, which for a cheap final
  stage lands several times below the continuous take and past the reach of
  any window drawn around it.
* Fixed-budget joint allocations cannot start integerization from
  `ceiling()`, which may bust the budget. Several roundings of the continuous
  optimum are repaired and trimmed to be both target-feasible and affordable,
  then improved by spending residual budget and by pairwise exchanges, and the
  best is kept. The field design satisfies `$operational$cost <= budget` with
  every target passing. The construction is greedy and can report failure on a
  problem that does have a feasible integer point.
* Cluster-mode `n_alloc()` integerizes at the PSU level: whole PSUs
  (`n_psu_int`) and whole takes (`n_per_psu_int`) per stratum, with the
  field cost `n_psu_int * (cost_psu + cost_ssu * n_per_psu_int)` kept
  within budget-mode budgets, and a targeted error when the budget
  cannot fund one PSU per stratum.
* Element allocations are rounded to match the solve mode: a target-`cv`
  solve rounds each stratum up (the integer design meets the target), a
  `budget` solve floors and then adds units by variance reduction per
  unit cost (the integer design stays within budget), and a fixed-`n`
  solve preserves the total with bounded largest-remainder rounding.
  Bounds are integerized first (`ceiling` of lower bounds, `floor` of
  upper bounds). Infeasible integer designs raise clear errors instead
  of silently violating `min_n_stratum`, `max_weight`, or the budget.
* `strata_bound()` reports the cv achieved by its integer allocation
  (the continuous optimum's cv is kept in `params$cv_continuous`).
* `n_cluster()` enforces realizable designs: fixed stage sizes must be
  at least 1, cost-optimal stage sizes below 1 are clamped to 1 with a
  warning, solved stage sizes below 1 are clamped with the achieved CV
  reported, and budgets too small for a single PSU raise an error.
* Domain identifiers are collision-free (values containing the display
  separator cannot merge distinct domains) and missing domain values
  raise an error instead of being silently dropped.
* `design_effect(method = "cr")` requires weights on the population
  scale (inverse inclusion probabilities): the sum of weights must
  exceed the sample size, overall and within every stratum (the
  offending stratum is named). Data where no cluster has within-cluster
  replication are rejected, and the returned components are validated
  as finite and non-negative.
* `varcomp()` rejects non-finite outcomes, missing stage identifiers,
  and out-of-range PPS probabilities. Per-observation probabilities must
  be constant within PSU, and named per-PSU probabilities are matched by
  PSU identifier.
* Three-stage `n_multi()` requires `icc_ssu` (matching `prec_multi()`),
  and simple-mode achieved precision is recomputed per indicator with the
  indicator's own method instead of sqrt-n rescaling.
* `strata_bound()` validates `unit_cost` length (1 or `n_strata`).
  When the cumulative-root-frequency boundaries degenerate on
  concentrated or discrete data, `cumrootf` falls back (with a warning)
  to boundaries between adjacent distinct values, so any input with at
  least `n_strata` distinct values yields nonempty strata. Fewer
  distinct values than strata is a targeted error.

## Display

* `print()` for a joint constrained allocation reports exceptions rather than
  confirmations. `status:` appears only when the solve is not optimal,
  `active bounds:` only when some are active, and a single binding constraint
  prints on one line rather than as a one-row table with a see-also pointer;
  the constraint total sits on the `field design:` line. A violated or
  multiply-binding solve still prints the full table. The common case drops
  from nine lines to five, and cost is formatted like every other print
  method rather than in scientific notation.
* `print()` always reports `deff`, including when it is left at its default
  of 1. A design effect of 1 asserts simple random sampling, which is the
  most consequential assumption a planning call can make silently, so the
  value is stated alongside `p`, `moe`, and the rest rather than suppressed.
* `svyplan_cluster` total sample size is the product of ceiled per-stage
  sizes in `print()`, `format()`, and `as.integer()`, so displayed totals
  are consistent with displayed stage sizes.
* `print()` and `format()` for `svyplan_cluster` show the unrounded
  continuous optimum alongside the operational total
  (e.g. `total n = 1190 (unrounded: 1159)`).
* New `as.double.svyplan_cluster()` method returns the continuous total
  (`x$total_n`). Use `as.integer()` for the operational total (fieldwork)
  and `as.double()` for the continuous optimum (mathematical solution).
* New `as.data.frame()` methods for `svyplan_n` and `svyplan_cluster`
  return the tabular form of a result (allocation detail, per-domain
  table, or stage table), the supported handoff to sampling packages
  such as `samplyr`.
