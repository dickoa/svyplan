# svyplan 0.4.9999

## Precision analysis

* `prec_prop()` — sampling precision (se, moe, cv) for a proportion given a
  sample size. Inverse of `n_prop()`. Uses the Cochran finite-population
  correction for Bernoulli proportions.

* `prec_mean()` — sampling precision for a mean given a sample size. Inverse
  of `n_mean()`.

* `prec_cluster()` — sampling precision (cv) for a multistage cluster
  allocation. Inverse of `n_cluster()`. Accepts `svyplan_cluster` objects
  directly and carries cost metadata for round-trip.

* `prec_multi()` — per-indicator sampling precision for multi-indicator
  survey designs. Inverse of `n_multi()`. Supports simple and multistage
  modes.

* New S3 class `svyplan_prec` with print, format, summary, and confint
  methods.

## Bidirectional S3 dispatch

* `n_prop()`, `n_mean()`, `n_cluster()`, and `n_multi()` are now S3 generics.
  Passing a `svyplan_prec` object recovers the sample size that produces the
  given precision. Passing a `svyplan_n` or `svyplan_cluster` object to the
  corresponding `prec_*()` function computes the achieved precision.

## Response rate adjustment

* All sample size (`n_prop`, `n_mean`, `n_cluster`, `n_multi`), precision
  (`prec_prop`, `prec_mean`, `prec_cluster`, `prec_multi`), and power
  (`power_prop`, `power_mean`) functions accept a `resp_rate` parameter.
  Sample sizes are inflated by `1 / resp_rate`; precision is computed on the
  effective sample `n * resp_rate`.

## Sensitivity analysis

* `predict()` methods for `svyplan_n`, `svyplan_cluster`, `svyplan_power`,
  and `svyplan_prec` objects. Evaluate the result at new parameter
  combinations (e.g. varying `deff`, `resp_rate`, `budget`) and return a
  data frame suitable for plotting.

## Confidence intervals

* `confint()` methods for `svyplan_n` and `svyplan_prec` objects. Returns
  a matrix with the expected confidence interval at the design's alpha level
  or a user-specified level.

## survey package integration

* When the survey package is installed, `survey::SE()` and `survey::cv()`
  methods are registered for `svyplan_n`, `svyplan_prec`, and
  `svyplan_cluster` objects.

## Precision fields on sample size objects

* `svyplan_n` objects now include `$se`, `$moe`, and `$cv` fields, computed
  automatically from the sample size and design parameters. Multi-indicator
  results have `NA` for these fields.

* `svyplan_cluster` objects include `$se` (NA) and `$cv` fields.

## Removed

* `cv_cluster()` has been removed. Use `prec_cluster(...)$cv` instead.

# svyplan 0.3.0

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
  Lavallee-Hidiroglou iterative (`"lh"`), and Kozak random search (`"kozak"`).
  Supports proportional, Neyman, optimal, and custom power allocation.
  Take-all (certainty) strata via the `certain` argument.

* Kozak random search optimized: inlined prefix-sum evaluation, pre-generated
  random numbers, and reduced bisection precision during search (full
  precision for final result). ~3x faster in cv-mode.

* New S3 class `svyplan_strata` with print, format, summary, as.data.frame,
  as.integer, and as.double methods.

# svyplan 0.1.0

Initial release.

## Sample size determination

* `n_prop()` — sample size for a proportion (Wald, Wilson, log-odds methods).
* `n_mean()` — sample size for a mean (moe and cv modes).
* `n_cluster()` — optimal multistage cluster allocation (2- and 3-stage,
  budget and cv modes).
* `n_multi()` — multi-indicator sample size for surveys with several
  precision targets. Supports simple (single-stage) and multistage modes,
  with automatic per-domain optimization.

## Design analysis

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
