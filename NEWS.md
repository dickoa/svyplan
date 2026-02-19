# svyplan 0.3.9999

## Power analysis

* `power_prop()` — power analysis for two-sample proportion tests. Solves for
  sample size, power, or minimum detectable effect (MDE). Supports design
  effect adjustment, finite population correction, one-sided and two-sided
  tests, and panel overlap for repeated surveys.

* `power_mean()` — power analysis for two-sample mean tests. Same solve modes
  and features as `power_prop()`, with analytical MDE solution.

* New S3 class `svyplan_power` with print, format, summary, as.integer, and
  as.double methods.

## Stratification

* `predict.svyplan_strata()` — apply strata boundaries to a numeric vector,
  returning a factor with customizable labels. Enables integration with
  samplyr's `stratify_by()` pipeline.

## Variance components

* `varcomp()` is now an S3 generic with methods for formulas, numeric vectors,
  and survey design objects.

* `varcomp.survey.design()` — extract cluster structure and outcome from a
  `survey::svydesign()` object. Requires the survey package.

## Design effects

* Chen-Rust design effect decomposition no longer requires the survey package.
  Uses a direct linearization (ultimate cluster) variance estimator.

* `design_effect()` with Henry or Spencer method now returns 1.0 for equal
  weights instead of NaN.

# svyplan 0.2.0

## Stratification

* `strata_bound()` — optimal strata boundary determination for a continuous
  stratification variable. Four methods: Dalenius-Hodges cumulative root
  frequency (`"cumrootf"`), geometric progression (`"geo"`),
  Lavallée-Hidiroglou iterative (`"lh"`), and Kozak random search (`"kozak"`).
  Supports proportional, Neyman, optimal, and custom power allocation.
  Take-all (certainty) strata via the `certain` argument.

* New S3 class `svyplan_strata` with print, format, summary, as.data.frame,
  as.integer, and as.double methods.

# svyplan 0.1.0

Initial release.

## Sample size determination

* `n_prop()` — sample size for a proportion (Wald, Wilson, log-odds methods).
* `n_mean()` — sample size for a mean (MOE and CV modes).
* `n_cluster()` — optimal multistage cluster allocation (2- and 3-stage,
  budget and CV modes).
* `n_multi()` — multi-indicator sample size for surveys with several
  precision targets. Supports simple (single-stage) and multistage modes,
  with automatic per-domain optimization.

## Design analysis

* `cv_cluster()` — coefficient of variation for a given multistage
  allocation (inverse of `n_cluster()`).
* `varcomp()` — variance component estimation from frame data via nested
  ANOVA (SRS and PPS, formula and vector interfaces).

## Design effects

* `design_effect()` — S3 generic for design effect estimation (Kish, Henry,
  Spencer, Chen-Rust methods, plus cluster planning mode).
* `effective_n()` — S3 generic for effective sample size.

## S3 classes

* `svyplan_n` — single-stage sample size results.
* `svyplan_cluster` — multistage allocation results.
* `svyplan_varcomp` — variance component estimates with fields
  `var_between`, `var_within`, `delta`, `k`, and `rel_var`.
