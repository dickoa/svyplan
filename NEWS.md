# svyplan 0.9.0

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
* The multi-indicator functions take the same four proportion methods as
  `n_prop()`, including `"beta"`. A `df` column carries the degrees of
  freedom of the planned variance estimator for the beta rows, so a table
  can mix an ordinary Wald indicator with a rare one sized from a
  Korn-Graubard interval. `df` is read on beta rows only and is rejected
  elsewhere rather than silently ignored.
* `n_twophase()` and `prec_twophase()`: allocation and precision for a
  two-phase (double) sample. One allocator covers two designs that are
  usually treated separately. Double sampling for stratification
  subsamples every phase-2 stratum; nonresponse follow-up carries the
  phase-1 respondents through untouched and subsamples only the
  nonrespondents. They differ only in which strata are marked
  `take_all`, and the follow-up case reduces to the standard
  `sqrt(c1 / (c2 * theta))`. The optimum is closed form, Neyman-shaped in
  `sd / sqrt(unit_cost)` but scaled by the between-stratum variance, so
  weak stratification pushes every fraction towards 1. Fractions above 1
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
* `n_mean()`: sample size for a mean (moe and cv modes).
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

## Precision analysis

* `prec_prop()`: sampling precision (se, moe, cv) for a proportion given a
  sample size. Precision counterpart to `n_prop()`.
* `prec_mean()`: sampling precision for a mean given a sample size. Precision
  counterpart to `n_mean()`.
* `prec_cluster()`: sampling precision (cv) for a multistage cluster
  allocation. Precision counterpart to `n_cluster()`.
* `prec_multi()`: per-indicator sampling precision for multi-indicator
  survey designs. Precision counterpart to `n_multi()`.
* `prec_multi_cluster()`: per-indicator precision for multistage
  multi-indicator designs. Stage costs are optional round-trip metadata and
  do not enter the precision calculation.
* `prec_alloc()`: sampling precision for a stratified allocation. Precision
  counterpart to `n_alloc()`. Joint allocation results retain one row per
  target and report stratum-bound violations when assessing an adopted sample.
  Joint assessment takes the same `min_n_stratum` floor as `n_alloc()`, so a design
  and a direct assessment of it report against identical bounds.
  When the fitted design carries an objective, the assessment also reports its
  components and weighted value. Minimum-cost designs invert through
  precision, pinning the achieved values as the requirement; budget-objective
  designs invert through cost, because pinning achieved precision as hard
  targets would over-constrain a design that already spends its whole budget.

All `n_*` and `prec_*` functions are S3 generics with bidirectional
round-trip: passing a precision object to the corresponding `n_*` function
recovers the continuous sample-size target under the same method and design
assumptions, and vice versa. Method changes and operational integer rounding
can change the result. Round-trip methods accept named
`...` overrides of any stored argument (e.g. `prec_prop(x, deff = 2)`).
A `NULL` value unsets a stored argument, and unknown names raise an error
instead of being silently ignored.

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
  outcome centred on zero even though every relvariance there is infinite.
  That case carries its own warning and is kept apart from a constant
  outcome, which is the one condition that genuinely leaves the components
  unidentified and `icc` at 0 by convention.
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
  listed there are ignored, so no unrecognised column is ever treated as a
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
  whole units.
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
