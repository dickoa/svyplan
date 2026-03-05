# svyplan 0.8.0

Initial CRAN release.

## Sample size determination

* `n_prop()` — sample size for a proportion (Wald, Wilson, log-odds methods).
* `n_mean()` — sample size for a mean (moe and cv modes).
* `n_cluster()` — optimal multistage cluster allocation (2- and 3-stage,
  budget and cv modes).
* `n_multi()` — multi-indicator sample size for surveys with several
  precision targets. Supports simple (single-stage) and multistage modes,
  with automatic per-domain optimization and `min_n` floor.
* `n_alloc()` — stratified sample allocation given a frame with stratum
  sizes and variabilities. Three solve modes: fixed total `n`, target `cv`,
  or `budget` constraint. Four allocation methods: proportional, Neyman,
  optimal (cost-weighted), and Bankier power allocation.

## Precision analysis

* `prec_prop()` — sampling precision (se, moe, cv) for a proportion given a
  sample size. Inverse of `n_prop()`.
* `prec_mean()` — sampling precision for a mean given a sample size. Inverse
  of `n_mean()`.
* `prec_cluster()` — sampling precision (cv) for a multistage cluster
  allocation. Inverse of `n_cluster()`.
* `prec_multi()` — per-indicator sampling precision for multi-indicator
  survey designs. Inverse of `n_multi()`.
* `prec_alloc()` — sampling precision for a stratified allocation. Inverse
  of `n_alloc()`.

All `n_*` and `prec_*` functions are S3 generics with bidirectional
round-trip: passing a precision object to the corresponding `n_*` function
recovers the sample size, and vice versa.

## Power analysis

* `power_prop()` — power analysis for two-sample proportion tests. Solves
  for sample size, power, or minimum detectable effect (MDE). Supports
  panel overlap for repeated surveys. MDE mode searches both directions
  (`p2 > p1` and `p2 < p1`) and returns the closest detectable alternative.
  `alternative` replaces `sides` (R standard naming). Supports unequal group
  sizes (`n = c(n1, n2)`), allocation ratio (`ratio`), and `method`
  parameter for arcsine and log-odds transforms (Valliant, 2018,
  Sections 4.3.4--4.3.5).
* `power_mean()` — power analysis for two-sample mean tests. Same solve
  modes and features as `power_prop()`. `alternative` replaces `sides`.
  Supports unequal group variances (`var = c(v1, v2)`), unequal group
  sizes (`n = c(n1, n2)`), and allocation ratio (`ratio`). Cohen's d
  conversion documented.
* `power_did()` — power analysis for difference-in-differences designs.
  Parametrized via `treat = c(baseline, endline)` and
  `control = c(baseline, endline)` vectors. Supports both proportion and
  mean outcomes, cell-specific variances, panel overlap, and all common
  design parameters.

## Stratification

* `strata_bound()` — optimal strata boundary determination for a continuous
  stratification variable. Four methods: Dalenius-Hodges cumulative root
  frequency (`"cumrootf"`), geometric progression (`"geo"`),
  Lavallee-Hidiroglou iterative (`"lh"`), and Kozak random search
  (`"kozak"`). Four allocation methods: proportional, Neyman, optimal
  (cost-weighted), and Bankier (1988) power allocation (`"power"`) with
  parameter `q` controlling the national/subnational precision trade-off.
  Take-all (certainty) strata via the `take_all` argument.
* `predict.svyplan_strata()` — apply strata boundaries to new data,
  returning a factor.

## Design components

* `varcomp()` — variance component estimation from frame data via nested
  ANOVA (SRS and PPS). S3 generic with methods for formulas, numeric
  vectors, and `survey::svydesign` objects.
* `design_effect()` — S3 generic for design effect estimation (Kish, Henry,
  Spencer, Chen-Rust decomposition, and cluster planning mode).
* `effective_n()` — S3 generic for effective sample size.

## Survey plan profiles

* `svyplan()` — create a reusable profile capturing shared design defaults
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

## Common features

* All sample size, precision, and power functions accept `deff` (design
  effect), `N` (finite population correction), and `resp_rate` (response
  rate adjustment).
* All functions accept `plan` — a `svyplan()` profile providing shared
  design defaults.
* `predict()` methods for sensitivity analysis: evaluate any result at new
  parameter combinations.
* `confint()` methods for `svyplan_n` and `svyplan_prec` objects.
* When the survey package is installed, `survey::SE()` and `survey::cv()`
  methods are registered automatically.

## S3 classes

* `svyplan` — survey plan profile (reusable design defaults).
* `svyplan_n` — sample size results (with se, moe, cv fields).
* `svyplan_cluster` — multistage allocation results.
* `svyplan_prec` — precision results.
* `svyplan_varcomp` — variance component estimates.
* `svyplan_strata` — strata boundary results.
* `svyplan_power` — power analysis results.

All classes have print and format methods.

## Input validation

* `n_multi()` and `prec_multi()` now reject non-positive `rel_var`, `k1`,
  and `k2` values in multistage mode.
* `prec_multi()` multistage mode now validates `delta1` (and `delta2` for
  3-stage) presence, type, NA, and range.
* `prec_cluster()` now validates that `rel_var` and `k` are positive and
  finite.
* `varcomp()` now rejects NA and empty outcome vectors.
* `confint()` methods for `svyplan_n` and `svyplan_prec` now validate that
  `level` is in (0, 1).
* `design_effect()` CR method now rejects mismatched vector lengths for
  `y`, `strvar`, and `clvar`.
* Weight validators (`design_effect()`, `effective_n()`) now reject
  non-finite values (`Inf`, `-Inf`).

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
