# svyplan 0.8.7

Initial CRAN release.

## Sample size determination

* `n_prop()`: sample size for a proportion (Wald, Wilson, log-odds methods).
* `n_mean()`: sample size for a mean (moe and cv modes).
* `n_cluster()`: optimal multistage cluster allocation (2- and 3-stage,
  budget and cv modes).
* `n_multi()`: multi-indicator sample size for surveys with several
  precision targets. Supports simple (single-stage) and multistage modes,
  with explicit per-domain optimization via `domains` and `min_n` floor.
* `n_alloc()`: stratified sample allocation given a frame with stratum
  sizes and variabilities. Three solve modes: fixed total `n`, target `cv`,
  or `budget` constraint. Four allocation methods: proportional, Neyman,
  optimal (cost-weighted), and Bankier power allocation. A `delta_psu`
  frame column switches to a stratified two-stage design (PSUs then
  elements per stratum) with cost-optimal or fixed per-stratum takes;
  all solve modes, methods, and constraints apply unchanged.

## Precision analysis

* `prec_prop()`: sampling precision (se, moe, cv) for a proportion given a
  sample size. Inverse of `n_prop()`.
* `prec_mean()`: sampling precision for a mean given a sample size. Inverse
  of `n_mean()`.
* `prec_cluster()`: sampling precision (cv) for a multistage cluster
  allocation. Inverse of `n_cluster()`.
* `prec_multi()`: per-indicator sampling precision for multi-indicator
  survey designs. Inverse of `n_multi()`.
* `prec_alloc()`: sampling precision for a stratified allocation. Inverse
  of `n_alloc()`.

All `n_*` and `prec_*` functions are S3 generics with bidirectional
round-trip: passing a precision object to the corresponding `n_*` function
recovers the sample size, and vice versa. Round-trip methods accept named
`...` overrides of any stored argument (e.g. `prec_prop(x, deff = 2)`);
a `NULL` value unsets a stored argument, and unknown names raise an error
instead of being silently ignored.

## Power analysis

* `power_prop()`: power analysis for two-sample proportion tests. Solves
  for sample size, power, or minimum detectable effect (MDE). Supports
  panel overlap for repeated surveys. MDE mode searches both directions
  (`p2 > p1` and `p2 < p1`) and returns the closest detectable alternative.
  `alternative` replaces `sides` (R standard naming). Supports unequal group
  sizes (`n = c(n1, n2)`), allocation ratio (`ratio`), and arcsine and
  log-odds transform methods (Valliant, 2018, Sections 4.3.4--4.3.5).
* `power_mean()`: power analysis for two-sample mean tests. Same solve
  modes and features as `power_prop()`. `alternative` replaces `sides`.
  Supports unequal group variances (`var = c(v1, v2)`), unequal group
  sizes (`n = c(n1, n2)`), and allocation ratio (`ratio`). Cohen's d
  conversion documented.
* `power_did()`: power analysis for difference-in-differences designs.
  Parametrized via `treat = c(baseline, endline)` and
  `control = c(baseline, endline)` vectors. Supports both proportion and
  mean outcomes, cell-specific variances, panel overlap, and all common
  design parameters.

## Stratification

* `strata_bound()`: optimal strata boundary determination for a continuous
  stratification variable. Four methods: Dalenius-Hodges cumulative root
  frequency (`"cumrootf"`), geometric progression (`"geo"`),
  Lavallee-Hidiroglou iterative (`"lh"`), and Kozak random search
  (`"kozak"`). Four allocation methods: proportional, Neyman, optimal
  (cost-weighted), and Bankier (1988) power allocation (`"power"`) with
  parameter `q` controlling the national/subnational precision trade-off.
  Take-all (certainty) strata via the `take_all` argument.
* `predict.svyplan_strata()`: apply strata boundaries to new data,
  returning a factor.

## Design components

* `varcomp()`: variance component estimation from frame data via nested
  ANOVA (SRS and PPS). S3 generic with methods for formulas, numeric
  vectors, and `survey::svydesign` objects. The `svydesign` method
  treats design weights as inverse inclusion probabilities: cluster
  sizes and totals are estimated by summed weights and the estimation
  variance of the weighted totals is subtracted from the between-stage
  terms, so unequal-probability samples from a previous round give
  approximately design-unbiased components (unit weights recover the
  frame formulas exactly). A `strata` argument estimates components per
  stratum, returning a table whose columns match the `n_alloc()` frame
  contract.
* `design_effect()`: S3 generic for design effect estimation (Kish, Henry,
  Spencer, Chen-Rust decomposition, and cluster planning mode).
* `effective_n()`: S3 generic for effective sample size.

## Survey plan profiles

* `svyplan()`: create a reusable profile capturing shared design defaults
  (`deff`, `N`, `resp_rate`, `alpha`, `stage_cost`, `unit_cost`, etc.). Pass
  to any function via `plan = plan` or pipe with `plan |> n_prop(...)`.
  Piping works with both positional and named arguments
  (e.g. `plan |> n_prop(p = 0.3, moe = 0.05)`).
  Explicit arguments always override plan defaults.

## Naming

* Cluster/multistage functions use `stage_cost` for per-stage cost vectors
  (`n_cluster()`, `n_multi()`, `prec_multi()`).
* Stratified allocation functions use `unit_cost` for per-stratum unit costs
  (`n_alloc()`, `prec_alloc()`, `strata_bound()`).

## Domain handling

* `n_multi()`, `prec_multi()`, `n_alloc()`, and `prec_alloc()` now require
  an explicit `domains` parameter to specify domain columns. Columns not
  listed in `domains` are silently ignored, eliminating the previous
  behaviour where any unrecognised column was automatically treated as a
  domain variable.
* `n_multi()` and `prec_multi()` results now store `params$domain_cols`,
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

* All sample size, precision, and power functions accept `deff` (design
  effect), `N` (finite population correction), and `resp_rate` (response
  rate adjustment). One shared variance equation is used throughout:
  with `n_net = n * resp_rate` responding units, `deff` multiplies the
  SRS variance at `n_net` and the finite population correction uses the
  actual sampling fraction `n_net / N` (so a census has zero sampling
  variance under any design effect). The Wilson proportion method has no
  finite-population form and ignores `N`. Targets that would require
  drawing more than `N` units from a finite frame raise an error instead
  of returning an impossible sample size, and the precision and power
  evaluators (`prec_prop()`, `prec_mean()`, `prec_multi()`, supplied-`n`
  `power_*()`) likewise reject a supplied gross `n` above the frame size
  of any group or indicator.
* All functions accept `plan`, a `svyplan()` profile providing shared
  design defaults.
* `predict()` methods for sensitivity analysis: evaluate any result at new
  parameter combinations.
* `confint()` methods for `svyplan_n` and `svyplan_prec` objects.
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

All classes have print and format methods.

## Input validation

* `n_multi()` and `prec_multi()` now reject non-positive `rel_var`, `k1`,
  and `k2` values in multistage mode.
* `prec_multi()` multistage mode now validates `delta1` (and `delta2` for
  3-stage) presence, type, NA, and range.
* `prec_cluster()` now validates that `rel_var` and `k` are positive and
  finite.
* `varcomp()` now rejects NA and empty outcome vectors, and data with a
  single PSU (overall or within a stratum) with a clear error naming
  the stratum.
* `confint()` methods for `svyplan_n` and `svyplan_prec` now validate that
  `level` is in (0, 1).
* `design_effect()` CR method now rejects mismatched vector lengths for
  `y`, `strvar`, and `clvar`.
* Weight validators (`design_effect()`, `effective_n()`) now reject
  non-finite values (`Inf`, `-Inf`).

## Feasibility and integer designs

* Constrained designs separate the continuous mathematical optimum
  (top-level fields) from the whole-unit field design (`$operational`,
  with cost, cv, and se recomputed from the integer design). `print()`
  leads with the field design; `as.integer()` returns the operational
  design in the same shape as `n` (stage vector for cluster plans,
  total for allocations) and `as.double()` its continuous counterpart.
* `n_cluster()` finds the operational design by discrete search
  (enumerated whole stage sizes): budget-mode field designs never
  exceed the budget and cv-mode field designs meet the target with
  whole units.
* Cluster-mode `n_alloc()` integerizes at the PSU level: whole PSUs
  (`n_psu_int`) and whole takes (`psu_size_int`) per stratum, with the
  field cost `n_psu_int * (cost_psu + cost_ssu * psu_size_int)` kept
  within budget-mode budgets, and a targeted error when the budget
  cannot fund one PSU per stratum.
* Element allocations are rounded to match the solve mode: a target-`cv`
  solve rounds each stratum up (the integer design meets the target), a
  `budget` solve floors and then adds units by variance reduction per
  unit cost (the integer design stays within budget), and a fixed-`n`
  solve preserves the total with bounded largest-remainder rounding.
  Bounds are integerized first (`ceiling` of lower bounds, `floor` of
  upper bounds); infeasible integer designs raise clear errors instead
  of silently violating `min_n`, `max_weight`, or the budget.
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
  and out-of-range PPS probabilities; per-observation probabilities must
  be constant within PSU, and named per-PSU probabilities are matched by
  PSU identifier.
* Three-stage `n_multi()` requires `delta_ssu` (matching `prec_multi()`),
  and simple-mode achieved precision is recomputed per indicator with the
  indicator's own method instead of sqrt-n rescaling.
* `strata_bound()` validates `unit_cost` length (1 or `n_strata`).
  When the cumulative-root-frequency boundaries degenerate on
  concentrated or discrete data, `cumrootf` falls back (with a warning)
  to boundaries between adjacent distinct values, so any input with at
  least `n_strata` distinct values yields nonempty strata; fewer
  distinct values than strata is a targeted error.

## Display

* `svyplan_cluster` total sample size is now computed as the product of
  ceiled per-stage sizes in `print()`, `format()`, and `as.integer()`,
  ensuring displayed totals are consistent with displayed stage sizes.
* `print()` and `format()` for `svyplan_cluster` now show the unrounded
  continuous optimum alongside the operational total
  (e.g. `total n = 1190 (unrounded: 1159)`).
* New `as.double.svyplan_cluster()` method returns the continuous total
  (`x$total_n`). Use `as.integer()` for the operational total (fieldwork)
  and `as.double()` for the continuous optimum (mathematical solution).
* New `as.data.frame()` methods for `svyplan_n` and `svyplan_cluster`
  return the tabular form of a result (allocation detail, per-domain
  table, or stage table), the supported handoff to sampling packages
  such as `samplyr`.
