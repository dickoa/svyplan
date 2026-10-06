# svyplan 0.15.0

Initial CRAN release.

svyplan sizes survey samples and reports the precision a size buys. It
reports stratum-level design quantities such as sizes, sampling fractions and
base weights. It never draws a sample and never attaches selection
probabilities or weights to units, which is the sampling stage's work.

Single-stage size, precision and power functions accept `deff`, `N` and
`resp_rate` and take a `svyplan()` profile as `plan`. Cluster functions take
`icc` and stage costs in place of `deff`. Every `n_` function has a `prec_`
counterpart, so passing a sized result to it reports the precision the design
was sized for, and named `...` overrides evaluate a variation of it.
Results carry `print()`, `summary()`, `format()`
and `as.data.frame()` methods, and `confint()`, `plot()` or `predict()` where
those apply. Optimization results with a whole-unit recommendation keep it in
`$operational`, separate from the continuous solution. Integer coercion uses
that recommendation when available, or rounds the continuous size up.

## Sample size determination

* `n_prop()`, `n_mean()` and `n_ratio()` size a proportion, a mean, or a ratio
  of two totals observed on the same units. Targets are `moe`, `cv`, or `rmoe`,
  a margin of error relative to the estimand. MICS uses relative margins of
  error, while DHS also describes precision using relative standard errors
  (`cv`). These targets differ by the confidence-level multiplier.
* Proportions offer four intervals, Wald, Wilson, log-odds, and the Korn and
  Graubard (1998) beta. Each size inverts the half-width its `prec_prop()`
  counterpart reports, so the four converge as the margin of error shrinks and
  separate only where the expected number of positive cases is small.
* `min_cases` sizes a proportion against a minimum expected number of positive
  cases, the constraint that binds for a rare outcome in a small domain. Every
  proportion result reports the `expected_cases` a design should yield.
* `n_ratio()` takes the quantities a pilot yields, the anticipated ratio with
  each component's coefficient of variation and their correlation across units.
  A ratio's design effect and homogeneity belong to the linearized variable
  `y - r * x` rather than to either component.
* `n_cluster()` allocates two- and three-stage cluster designs in budget and cv
  modes. Losses are named for the stage they act on, `resp_rate_psu`,
  `resp_rate_ssu`, and `resp_rate`, because a lost cluster and a lost ultimate
  unit move the cost-optimal take differently. Any stage size may be fixed, and
  the remaining stages solve that restricted problem.
* `n_alloc()` allocates a stratified sample from a frame of stratum sizes and
  variabilities, in fixed-`n`, target-`cv`, or `budget` mode, by proportional,
  Neyman, optimal, or Bankier power allocation. An `icc_psu` frame column makes
  it a stratified two-stage design, where `fpc` picks between the ultimate-unit
  correction, the exact stage-by-stage one, and none.
  Bankier allocation weights each stratum's population CV by a size measure
  raised to `alloc_q`. The measure is the stratum total `N * abs(mean)` unless
  a positive `alloc_measure` frame column supplies another. With that default
  measure, `alloc_q = 1` gives Neyman allocation. In cluster mode each stratum
  fields whole PSUs at a whole take near its cost-optimal value.
* Long `measures` and `targets` tables request joint constrained allocation
  across indicators and overlapping domains. Without `objective` the result is
  the cheapest design meeting every requirement, and with `objective` and
  `budget` it is the best design that budget buys. A KKT-certified continuous
  optimum is reported separately from the feasible whole-unit recommendation.
  An element-level proportion measure uses the population variance
  `N * p * (1 - p) / (N - 1)`, so a one-stratum joint allocation equals
  `n_prop()`. Multistage designs use the working variance `p * (1 - p)`.
  A `frame` column `deff` or `resp_rate` is the stratum default under the
  `measures` rows, as for the stage parameters, and indicator values in
  `frame` are refused.
* `psu` gives `n_alloc()` a PSU register, which separates the PSUs a PPS design
  would select with certainty and reports the take each fielded PSU carries.
  A stratum with a single PSU carries no between-PSU variance, so its
  `icc_psu` may be `NA`. `certainty_cutoff` takes PSUs with certainty from an
  inclusion probability below one, for every stratum or per stratum.
  `n_psu_per_zone` cuts each stratum's remainder into zones of about equal
  size and draws one or two PSUs from each. With two, every stratum's
  variance can be estimated. With one, `$psu$.pair` groups the zones for a
  collapsed variance. The whole-unit design is assessed on its own rounded
  takes and repaired where a target fails, or the fit stops naming it.
  Certainty visits are priced as a fixed cost, and the clustering reads the
  responding take. A `budget`, a `var_ratio_psu` other than 1 and
  a PSU-level response rate are refused with a register. A register PSU
  smaller than its take is refused, and `merge_psus()` merges small PSUs with
  their neighbours in row order to a minimum size.
* `n_multi()` sizes a table of indicators, with optional per-domain sizing and
  a `min_n_domain` floor. `n_cluster(indicators = )` and
  `prec_cluster(indicators = )` take the same table for a multistage design, so
  a cluster argument never changes what `n_multi()` returns.
* `n_twophase()` allocates a two-phase sample. One allocator covers double
  sampling for stratification and nonresponse follow-up, which differ only in
  which strata are marked `take_all`. Every result prices the single-phase
  design that skips phase 1, since two-phase sampling is not always the better
  buy. `assurance` converts required respondents into an issue clearing a
  stated binomial probability.
* `n_change()` and `n_pooled()` size an occasion of a repeated survey, for a
  change between two occasions and for an estimate pooled across several.
  Precision comes from the product of `overlap` and `overlap_cor`, so a full
  panel of uncorrelated measurements saves nothing.
* `n_panel()` sizes the recruitment a panel needs to still deliver a required
  responding sample after several waves of attrition. It embeds an `n_mean()`,
  `n_prop()`, `prec_mean()`, or `prec_prop()` result rather than restating its
  arguments. A fixed panel reports `n_issued`, and a rotating one `n_entrants`
  per occasion beside the `n_in_sample` its live cohorts hold between them.
  `start` reports what a rotating design delivers while it climbs to its steady
  state, and `assurance` again converts an expected-value recruitment into an
  assured one.
* `design_schedule()` turns a rotating launch into a versioned, package-neutral
  operational manifest, validated field by field, for a sampling package to
  field.

## Precision analysis

* Each `n_*` function has a `prec_*` counterpart reporting `se`, `moe`, `cv`,
  and `rmoe` for a size already chosen, namely `prec_prop()`, `prec_mean()`,
  `prec_ratio()`, `prec_cluster()`, `prec_alloc()`, `prec_multi()`,
  `prec_twophase()`, `prec_change()`, `prec_pooled()`, and `prec_panel()`.
* `prec_prop()` and `prec_mean()` also solve the remaining direction. Supplying
  `cv` or `rmoe` in place of the level returns the smallest estimand the design
  measures that precisely, which is the planning form of a publication
  threshold.
* `prec_alloc()` reports at three levels, the population overall, each stratum
  in `$detail` with its share of the design variance, and each domain in
  `$domains`, using the same table `n_alloc()` returns so a design and its
  assessment compare row for row. For a cluster allocation, `n_per_psu` sets
  the take, so the whole-unit counts and whole take of `$operational`
  reproduce its precision exactly.

## Power analysis

* `power_prop()` and `power_mean()` solve a two-sample test for sample size,
  power, or minimum detectable effect, with unequal group sizes and variances,
  an allocation ratio, a two-sided or one-sided `alternative`, and arcsine and
  log-odds transforms.
* `power_did()` covers difference-in-differences, parametrized by `treat` and
  `control` baseline and endline pairs, for proportion and mean outcomes.
* All three accept panel overlap between the two occasions of a repeated
  survey. Bernoulli outcomes use the same finite-population variance as the
  precision family on every scale, Wald, arcsine, and log-odds alike, and
  correlations between repeated Bernoulli outcomes are restricted to the
  attainable Frechet-Hoeffding range.
* Solved sample sizes are at least two units per group and report the power
  achieved at that size. A requested power used to solve a size or minimum
  detectable effect must exceed `alpha`, the power at zero effect, and
  two-sided solutions include both rejection tails. A design that enumerates
  its population has no minimum detectable effect and says so.

## Stratification

* `strata_bound()` builds candidate boundaries for a continuous stratification
  variable by cumulative root frequency, geometric progression, LH coordinate
  optimization, or Kozak's random search, under the same four allocation
  methods `n_alloc()` offers. `$strata` matches the `n_alloc()` frame contract,
  so the table can be handed straight on.
* Bounded allocation satisfies the requested total together with minimum,
  maximum, and take-all constraints. Boundary searches reserve the census
  take correctly, including the full-population scaling of a target CV.
* `predict()` applies fitted boundaries to new data, returning a factor.

## Design components

* `varcomp()` estimates variance components by nested ANOVA, with methods for
  formulas, numeric vectors, and `survey::svydesign` objects. Weights are read
  as inverse inclusion probabilities, so a previous round's unequal-probability
  sample gives approximately design-unbiased components. `varcomp(deff = )`
  backs an `icc` out of a published design effect instead.
* `design_effect()` builds the design effect a planned design should produce
  from clustering, unequal weighting, and stratification components, reading
  them off a plan where one is supplied. It anticipates a design effect from
  design parameters and does not estimate a realized one from collected data,
  which is `survey::svymean(deff = TRUE)`.
* `effective_n()` mirrors it, returning `n * resp_rate / deff`.
* `design_df()` counts the degrees of freedom a planned design's variance
  estimator will have, the planning analogue of `survey::degf()`, with additive
  per-stratum and non-additive per-domain detail. A zoned PSU register counts
  its zones, or its collapsed zone groups at one PSU per zone, whose
  per-stratum counts need not add up when a group joins strata. `df` switches
  the interval quantile from the normal to `t` in every function that builds
  an interval.
* `design_rotation()` declares the occasions a unit spends in and out of
  sample over its whole life, from a compact spec in either notation the
  literature uses or from an explicit per-occasion take. It is the object
  `design_overlap()` and `design_schedule()` both read, so a pattern is parsed
  in one place.
* `design_overlap()` reports the issued-sample overlap that rotation produces
  at every lag. `plot()` draws the rotation chart the design would be
  published as.

## Survey plan profiles

* `svyplan()` captures shared design defaults such as `deff`, `N`, `resp_rate`,
  stage-specific response rates, `df`, `alpha`, `stage_cost`, and `unit_cost`,
  validated at construction. Pass it as `plan = plan` or pipe with
  `plan |> n_prop(...)`. Explicit arguments override a profile, and `update()`
  revises one.

## Naming

* Argument names say what a quantity is rather than reproducing one textbook's
  symbol. `icc` is the design-based measure of homogeneity within clusters,
  `var_ratio` the ratio of a stage component's unit variance to the analysis
  variable's, and `unit_relvar` the unit relvariance. `?svyplan-package` carries
  the full symbol map.
* Stage sizes distinguish counts from takes, `n_psu` against `n_per_psu` and
  `n_per_ssu`. Sample-size floors are named for what they count, `min_n_domain`
  against `min_n_stratum`. Cluster functions cost stages through `stage_cost`
  and stratified ones cost strata through `unit_cost`.
* Wherever a mean is planned, dispersion is given as `var` or `sd`, and exactly
  one of the two.

## Validation

* Arguments a configuration cannot use are rejected rather than silently
  discarded, whether passed by name or through `...`, so a misspelling raises an
  error instead of computing under an unintended default. Optional arguments sit
  after `...` and must be fully named.
* Targets requiring more than `N` units, sizes above a frame, out-of-range rates
  and correlations, non-whole counts, and stage parameters a design does not
  have are all refused where they are supplied, rather than surfacing later as
  an optimizer failure or an infinite precision.
* Infeasibility is reported with the binding quantity named. A joint allocation
  distinguishes targets unattainable at the stratum bounds, targets attainable
  but unaffordable, and the absence of a whole-unit design meeting every target
  inside the budget.

## S3 classes

* `svyplan` for a design profile, `svyplan_n` and `svyplan_prec` for sizes and
  precision, `svyplan_cluster` for multistage allocations, `svyplan_twophase`,
  `svyplan_panel` and `svyplan_schedule` for the longitudinal family,
  `svyplan_power`, `svyplan_strata`, `svyplan_varcomp`, `svyplan_deff`,
  `svyplan_df`, and `svyplan_overlap`.
* All classes have `print()` and `format()` methods. `print()` reports the
  answer and the design behind it, and `summary()` carries the derived detail,
  the assumptions, and the diagnostics.
* `svyplan_deff`, `svyplan_df`, and `svyplan_overlap` are usable as the plain
  numbers they stand for while carrying the decomposition each was built from.
