---
output: github_document
---



# svyplan

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/svyplan)](https://CRAN.R-project.org/package=svyplan)
[![R-CMD-check](https://gitlab.com/dickoa/svyplan/badges/main/pipeline.svg)](https://gitlab.com/dickoa/svyplan/-/pipelines)
[![Codecov test coverage](https://codecov.io/gl/dickoa/svyplan/branch/main/graph/badge.svg)](https://app.codecov.io/gl/dickoa/svyplan?branch=main)
<!-- badges: end -->

Survey sample size determination, precision analysis, optimal and joint
multivariate/multidomain allocation, stratification, and power analysis for R.

## Installation

```r
# From GitLab
pak::pkg_install("gitlab::dickoa/svyplan")
```

## Sample sizes


``` r
library(svyplan)

# Proportion with margin of error
n_prop(p = 0.3, moe = 0.05)
#> Sample size for proportion (wald)
#> n = 323 (p = 0.30, moe = 0.050, deff = 1)

# Mean with finite population and design effect
n_mean(var = 100, moe = 2, N = 5000, deff = 1.5)
#> Sample size for mean
#> n = 141 (var = 100.00, moe = 2.000, deff = 1.50)
```

### Response rate adjustment

Most sizing and precision functions accept `resp_rate`. In sample-size
mode, the required sample is inflated by `1 / resp_rate` to account for
expected non-response:


``` r
n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8)
#> Sample size for proportion (wald)
#> n = 606 (net: 485) (p = 0.30, moe = 0.050, deff = 1.50, resp_rate = 0.80)
```

### Survey plan profiles

When the same design parameters apply across many calls, bundle them into
a `svyplan()` profile:


``` r
plan <- svyplan(deff = 1.5, resp_rate = 0.85, N = 50000)

# Pass as argument
n_prop(p = 0.3, moe = 0.05, plan = plan)
#> Sample size for proportion (wald)
#> n = 564 (net: 480) (p = 0.30, moe = 0.050, deff = 1.50, resp_rate = 0.85)

# Or pipe (positional or named args)
plan |> n_mean(100, moe = 2)
#> Sample size for mean
#> n = 169 (net: 144) (var = 100.00, moe = 2.000, deff = 1.50, resp_rate = 0.85)
plan |> n_mean(var = 100, moe = 2)
#> Sample size for mean
#> n = 169 (net: 144) (var = 100.00, moe = 2.000, deff = 1.50, resp_rate = 0.85)

# Explicit args always override plan defaults
n_prop(p = 0.3, moe = 0.05, plan = plan, deff = 2.0)
#> Sample size for proportion (wald)
#> n = 750 (net: 638) (p = 0.30, moe = 0.050, deff = 2.00, resp_rate = 0.85)
```

## Precision analysis

Given a sample size, how precise will your estimates be? The `prec_*()`
functions are the inverse of `n_*()`:


``` r
prec_prop(p = 0.3, n = 400)
#> Sampling precision for proportion (wald)
#> n = 400
#> se = 0.0229, moe = 0.0449, cv = 0.0764

prec_mean(var = 100, n = 400, mu = 50)
#> Sampling precision for mean
#> n = 400
#> se = 0.5000, moe = 0.9800, cv = 0.0100
```

### Round-trip between size and precision

All `n_*()` and `prec_*()` functions are S3 generics. Pass a precision
result to `n_*()` to recover the sample size, or pass a sample size result
to `prec_*()` to compute the achieved precision:


``` r
# Start with a precision target
s <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)

# What precision does this n achieve?
p <- prec_prop(s)
p
#> Sampling precision for proportion (wald)
#> n = 485
#> se = 0.0255, moe = 0.0500, cv = 0.0850

# Recover the original n
n_prop(p)
#> Sample size for proportion (wald)
#> n = 485 (p = 0.30, moe = 0.050, deff = 1.50)
```

## Multi-indicator surveys

Household surveys track many indicators at once. `n_multi()` finds the
sample size that satisfies all precision targets simultaneously.


``` r
targets <- data.frame(
  name = c("stunting", "vaccination", "anemia"),
  p    = c(0.25, 0.70, 0.12),
  moe  = c(0.05, 0.05, 0.03),
  deff = c(2.0, 1.5, 2.5)
)

n_multi(targets)
#> Multi-indicator sample size
#> n = 1127 (binding: anemia)
#> ---
#>  name        .n   .cv_target .cv_achieved .binding
#>  stunting     577 0.10204269 0.07297042           
#>  vaccination  485 0.03644382 0.02388518           
#>  anemia      1127 0.12755336 0.12755336   *
```

Per-domain optimization works by specifying domain columns via the
`domains` parameter.

### MICS/DHS-style relative margin of error

Programmes like UNICEF MICS and DHS express precision as a **relative
margin of error** (RME = MOE / p). To use this with svyplan, convert to
an absolute margin of error: `moe = RME * p`.


``` r
# RME = 12% for each indicator
rme <- 0.12
targets_rme <- data.frame(
  name = c("stunting", "vaccination", "anemia"),
  p    = c(0.25, 0.70, 0.12),
  deff = c(2.0, 1.5, 2.5)
)
targets_rme$moe <- rme * targets_rme$p

n_multi(targets_rme)
#> Multi-indicator sample size
#> n = 4891 (binding: anemia)
#> ---
#>  name        .n   .cv_target .cv_achieved .binding
#>  stunting    1601 0.06122561 0.03502580           
#>  vaccination  172 0.06122561 0.01146488           
#>  anemia      4891 0.06122561 0.06122561   *
```

The MICS template reports sample size in **households**. svyplan returns
the number of **individuals** in the target population. To convert,
divide by the expected number of eligible individuals per household:
`n_hh = ceiling(n / (pb * hh_size))`, where `pb` is the share of the
target population and `hh_size` is the average household size.

## Multistage cluster designs


``` r
# Optimal 2-stage allocation within a budget
n_cluster(stage_cost = c(500, 50), icc = 0.05, budget = 100000)
#> Optimal 2-stage allocation
#> field design: n_psu = 80 | n_per_psu = 15 -> total n = 1200
#> cv = 0.0376, cost = 100000
#> continuous optimum: n_psu = 84.08997 | n_per_psu = 13.78405 (cv = 0.0376, cost = 100000)

# Precision for a given allocation
prec_cluster(n = c(50, 12), icc = 0.05)
#> Sampling precision for 2-stage cluster
#> n_psu = 50 | n_per_psu = 12 -> total n = 600
#> cv = 0.0508
```

Variance components can be estimated from frame data and passed
directly to `n_cluster()`:


``` r
set.seed(104)
frame <- data.frame(
  district = rep(1:40, each = 20),
  income = rep(rnorm(40, 500, 100), each = 20) + rnorm(800, 0, 50)
)

vc <- varcomp(income ~ district, data = frame)
vc
#> Variance components (2-stage)
#> varb = 0.0255, varw = 0.0099
#> icc = 0.7210
#> var_ratio = 1.0317
#> Unit relvariance = 0.0343
as.data.frame(vc)
#>   stages       varb        varw       icc var_ratio unit_relvar
#> 1      2 0.02554834 0.009886668 0.7209915   1.03174   0.0343449

n_cluster(stage_cost = c(500, 50), icc = vc, cv = 0.05)
#> Optimal 2-stage allocation
#> field design: n_psu = 13 | n_per_psu = 2 -> total n = 26
#> cv = 0.0484, cost = 7800
#> continuous optimum: n_psu = 12.22966 | n_per_psu = 1.967178 (cv = 0.0500, cost = 7318)
```

`icc` is the survey-planning measure of homogeneity used by
`varcomp()`, `n_cluster()`, and `design_effect()`. It is not the same as
a generic mixed-model ICC, and values near 0 or 1 correspond to
degenerate boundary cases for the closed-form cluster optimizer.

## Sensitivity analysis

`predict()` evaluates a result at new parameter combinations, returning
a data frame suitable for plotting:


``` r
x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5)
predict(x, expand.grid(
  deff = c(1, 1.5, 2, 2.5),
  resp_rate = c(0.7, 0.8, 0.9, 1.0)
))
#>    deff resp_rate         n         se  moe         cv
#> 1   1.0       0.7  460.9751 0.02551067 0.05 0.08503558
#> 2   1.5       0.7  691.4626 0.02551067 0.05 0.08503558
#> 3   2.0       0.7  921.9501 0.02551067 0.05 0.08503558
#> 4   2.5       0.7 1152.4376 0.02551067 0.05 0.08503558
#> 5   1.0       0.8  403.3532 0.02551067 0.05 0.08503558
#> 6   1.5       0.8  605.0298 0.02551067 0.05 0.08503558
#> 7   2.0       0.8  806.7064 0.02551067 0.05 0.08503558
#> 8   2.5       0.8 1008.3829 0.02551067 0.05 0.08503558
#> 9   1.0       0.9  358.5362 0.02551067 0.05 0.08503558
#> 10  1.5       0.9  537.8042 0.02551067 0.05 0.08503558
#> 11  2.0       0.9  717.0723 0.02551067 0.05 0.08503558
#> 12  2.5       0.9  896.3404 0.02551067 0.05 0.08503558
#> 13  1.0       1.0  322.6825 0.02551067 0.05 0.08503558
#> 14  1.5       1.0  484.0238 0.02551067 0.05 0.08503558
#> 15  2.0       1.0  645.3651 0.02551067 0.05 0.08503558
#> 16  2.5       1.0  806.7064 0.02551067 0.05 0.08503558
```

Sensitivity analysis is available for single-indicator sample-size and
precision results, cluster designs, power analyses, and strata boundaries.
Multi-indicator results are not currently supported by `predict()`.

## Strata boundaries

`strata_bound()` constructs candidate boundaries for a continuous
stratification variable.


``` r
set.seed(905)
x <- rlnorm(5000, meanlog = 6, sdlog = 1.2)

strata_bound(x, n_strata = 4, n = 300, method = "cumrootf")
#> Strata boundaries (Dalenius-Hodges, 4 strata)
#> Boundaries: 400.0, 1300.0, 3200.0
#> n = 300, cv = 0.0205
#> Allocation: neyman
#> ---
#>  stratum lower       upper    N    share sd     mean      n  
#>  1          8.380438   400.00 2492 0.498 104.2   186.9325  46
#>  2        400.000000  1300.00 1647 0.329 246.5   724.3771  71
#>  3       1300.000000  3200.00  639 0.128 500.9  1946.0399  56
#>  4       3200.000000 28909.53  222 0.044 3264.3 5471.5402 127
```

Four methods are available: Dalenius-Hodges (`"cumrootf"`), geometric
(`"geo"`), LH-inspired coordinate optimization (`"lh"`), and
Kozak-inspired random-restart local search (`"kozak"`). The latter two are
heuristics and do not claim global optimality or exact implementation of the
published algorithms.

## Two-phase designs

A large cheap phase 1, then a subsample measured on the expensive
variable. `n_twophase()` allocates both at once. The frame is one row per
phase-2 stratum, in the `n_alloc()` column vocabulary but on a narrower
contract: `sd` is required where `n_alloc()` also accepts `var`.


``` r
frame <- data.frame(
  stratum   = c("A", "B", "C", "D"),
  N         = c(3500, 2500, 2500, 1500),
  sd        = c(12, 25, 8, 40),
  mean      = c(40, 70, 35, 90),
  unit_cost = c(2, 5, 1, 9)
)

n_twophase(frame, phase1_cost = 1, budget = 50000)
#> Two-phase allocation (4 phase-2 strata)
#> issued: n_phase1 = 16925 | n_phase2 = 8092
#> cv = 0.0050, cost = 50000
#> ---
#>  stratum share    sd unit_cost     nu n_issued n_int
#>        A 0.350 12.00      2.00 0.4154     2461  2462
#>        B 0.250 25.00      5.00 0.5474     2316  2316
#>        C 0.250  8.00      1.00 0.3917     1657  1659
#>        D 0.150 40.00      9.00 0.6528     1657  1657
#> field design: n_phase1 = 16924 | n_phase2 = 8094 (cost 50000, cv 0.0050)
#> 
#> Single-phase is better here: n = 14085, cv = 0.0046, cost = 50000
#> Skip phase 1 and measure directly.
```

Nonresponse follow-up is the same problem with two strata: respondents
are already measured, so they cost nothing more and are all kept, and
only the nonrespondents are subsampled.


``` r
theta <- 0.5

nrfu <- data.frame(
  stratum   = c("respondents", "nonrespondents"),
  N         = c(theta, 1 - theta),
  sd        = c(1, 1),
  unit_cost = c(0, 200),
  take_all  = c(TRUE, FALSE)
)

n_twophase(nrfu, phase1_cost = 50, budget = 100000,
           mu = 1, single_cost = 50 / theta)
#> Two-phase allocation (2 phase-2 strata)
#> issued: n_phase1 = 828 | n_phase2 = 707
#> cv = 0.0382, cost = 1e+05
#> ---
#>         stratum share   sd unit_cost     nu n_issued n_int
#>     respondents 0.500 1.00      0.00 1.0000      414   414
#>  nonrespondents 0.500 1.00    200.00 0.7071      293   293
#>  take_all
#>         *
#>          
#> field design: n_phase1 = 828 | n_phase2 = 707 (cost 1e+05, cv 0.0382)
#> 
#> Single-phase is better here: n = 1000, cv = 0.0316, cost = 1e+05
#> Skip phase 1 and measure directly.
```

Neither phase has to be a simple random sample, and both can lose sample
to nonresponse. `phase1_deff` and `resp_rate` act on the between-stratum
component; `deff` and `resp_rate` frame columns act on each stratum's
residual; `single_deff` and `single_resp_rate` describe the comparator.
Each divides its own component, so a single factor on the combined
variance would be the wrong model. Issued and expected-responding counts
are reported separately.

Every result carries the single-phase comparison, because two-phase
sampling is not always an improvement.

## Power analysis

Solve for sample size, power, or minimum detectable effect. Supports
design effects, finite population correction, response rate adjustment,
panel overlap, unequal groups, and allocation ratios. Arcsine and
log-odds methods available for rare proportions.


``` r
# Sample size to detect a 5pp change from 70% with deff = 2
power_prop(p1 = 0.70, p2 = 0.75, deff = 2.0)
#> Power analysis for proportions (solved for sample size)
#> n = 2496 (per group), power = 0.800, effect = 0.0500
#> (p1 = 0.700, p2 = 0.750, alpha = 0.05, deff = 2.00)

# MDE with n = 1500 per group
power_prop(p1 = 0.70, n = 1500, deff = 2.0)
#> Power analysis for proportions (solved for minimum detectable effect)
#> n = 1500 (per group), power = 0.800, effect = 0.0639
#> (p1 = 0.700, p2 = 0.764, alpha = 0.05, deff = 2.00)

# Means
power_mean(200, effect = 5)
#> Power analysis for means (solved for sample size)
#> n = 126 (per group), power = 0.800, effect = 5.0000
#> (alpha = 0.05, deff = 1)

# Arcsine method for rare proportions
power_prop(p1 = 0.15, p2 = 0.18, alternative = "one.sided",
           method = "arcsine")
#> Power analysis for proportions (solved for sample size)
#> n = 1890 (per group), power = 0.800, effect = 0.0300
#> (p1 = 0.150, p2 = 0.180, alpha = 0.05, deff = 1, one-sided, method = arcsine)

# Difference-in-differences
power_did(treat = c(0.50, 0.55), control = c(0.50, 0.48),
          outcome = "prop", effect = 0.07)
#> Power analysis for DiD proportions (solved for sample size)
#> n = 1598 (per group), power = 0.800, effect = 0.0700
#> (treat = (0.500, 0.550), control = (0.500, 0.480), alpha = 0.05, deff = 1)
```

`plot()` draws the power-vs-sample-size curve with reference lines at the solved point:


``` r
pw <- power_prop(p1 = 0.70, p2 = 0.75, power = 0.80, deff = 2.0)
plot(pw)
```

<div class="figure">
<img src="man/figures/README-power-plot-1.png" alt="Power increases with total sample size. Dashed reference lines mark 80 percent power at the required sample size for detecting a change from 70 to 75 percent with a design effect of 2."  />
<p class="caption">plot of chunk power-plot</p>
</div>


## Stratified allocation

Given a sampling frame with stratum sizes and variabilities, `n_alloc()`
distributes the total sample across strata:


``` r
frame <- data.frame(
  N = c(4000, 3000, 3000),
  sd = c(10, 15, 8),
  mean = c(50, 60, 55)
)

n_alloc(frame, n = 600, alloc = "neyman")
#> Stratum allocation (neyman, 3 strata)
#> field design: n = 600, cv = 0.0079, cost = 600
#> continuous optimum: n = 600, cv = 0.0079, se = 0.4305
#> (deff = 1)
```

Constraints and alternative solve modes are also supported:


``` r
frame_constraints <- transform(
  frame,
  unit_cost = c(1, 1.5, 1),
  max_weight = c(25, 20, NA),
  take_all = c(FALSE, FALSE, TRUE)
)

# Budget-constrained allocation with weight and take-all constraints
n_alloc(frame_constraints, budget = 3500, alloc = "optimal", min_n_stratum = 40)
#> Stratum allocation (optimal, 3 strata)
#> field design: n = 3403, cv = 0.0076, cost = 3500
#> continuous optimum: n = 3403.425, cv = 0.0076, se = 0.4125
#> (min_n_stratum = 40, deff = 1)
```

Domain-level CV targets can be enforced via the `domains` parameter:


``` r
frame_domains <- data.frame(
  province = c("North", "North", "South", "South"),
  stratum = c("Urban", "Rural", "Urban", "Rural"),
  N = c(2000, 3000, 1800, 3200),
  sd = c(12, 18, 10, 16),
  mean = c(55, 48, 58, 50)
)

# Minimum total n such that each province meets the CV target
n_alloc(frame_domains, domains = "province",
        cv = 0.04, alloc = "power", alloc_q = 0.3)
#> Stratum allocation (power, 4 strata)
#> field design: n = 112, cv = 0.0270, cost = 112
#> continuous optimum: n = 110.7422, cv = 0.0272, se = 1.4076
#> (deff = 1)
#> Domains: 2
#> ---
#>  province .domain .n       .se      .moe     .cv    .cost
#>  North    5_North 59.23404 2.032000 3.982647 0.0400 59   
#>  South    5_South 51.50815 1.948447 3.818886 0.0368 52
```

For several indicators and overlapping domains, pass long `measures` and
`targets` tables. The allocation rule is then determined jointly by the
precision requirements:


``` r
joint_frame <- data.frame(
  stratum = c("North urban", "North rural", "South urban", "South rural"),
  region = rep(c("North", "South"), each = 2),
  residence = rep(c("Urban", "Rural"), 2),
  N = c(1000, 1800, 1200, 1600),
  unit_cost = c(1, 1.2, 1.4, 1.1)
)
joint_measures <- data.frame(
  stratum = rep(joint_frame$stratum, 2),
  name = rep(c("coverage", "income"), each = 4),
  p = c(0.50, 0.40, 0.60, 0.35, rep(NA, 4)),
  mean = c(rep(NA, 4), 50, 55, 60, 48),
  sd = c(rep(NA, 4), 10, 12, 15, 9)
)
joint_targets <- data.frame(
  name = c("coverage", "coverage", "income"),
  domain = c(".overall", "region", "residence"),
  level = c(NA, "North", "Urban"),
  cv = c(0.05, 0.08, NA),
  moe = c(NA, NA, 2)
)

joint_fit <- n_alloc(
  joint_frame, measures = joint_measures, targets = joint_targets
)
joint_fit
#> Joint constrained allocation (Bethel)
#> question: cheapest design meeting every precision target
#> field design: n = 442, cost = 518 (3 targets, all pass)
#> continuous optimum: n = 441.7832, cost = 517 (integerizing costs +0.05%)
#> binding: coverage@.overall:cv (target 0.05, achieved 0.05)
```

The fitted object distinguishes the certified continuous optimum from a
feasible integer operational recommendation. Use `prec_alloc(joint_fit)` or
`prec_alloc(joint_fit, n = joint_fit$detail$n_int)` to inspect either design.

The same API handles fixed-take multistage designs. Stage populations, fixed
takes, and costs belong in `frame`; indicator-specific homogeneity parameters
belong in `measures`. Here only the PSU counts are optimized:


``` r
cluster_frame <- within(joint_frame, {
  unit_cost <- NULL
  N_psu <- c(100, 150, 100, 120)
  n_per_psu <- c(8, 10, 12, 8)
  cost_psu <- c(300, 400, 450, 350)
  cost_ssu <- c(25, 30, 35, 28)
})
cluster_measures <- transform(
  joint_measures,
  icc_psu = rep(c(0.03, 0.05, 0.08, 0.04), 2),
  var_ratio_psu = 1
)
cluster_fit <- n_alloc(
  cluster_frame, measures = cluster_measures, targets = joint_targets
)
cluster_fit$detail[, c("stratum", "n_psu_int", "n_per_psu", "n_int")]
#>       stratum n_psu_int n_per_psu n_int
#> 1 North urban        14         8   112
#> 2 North rural        20        10   200
#> 3 South urban        13        12   156
#> 4 South rural        19         8   152
```

Public `n` remains the ultimate-unit sample size. Thus the field design obeys
`n_int = n_psu_int * n_per_psu`; three-stage designs additionally multiply by
the fixed `n_per_ssu`.

`prec_alloc()` computes the precision for a given allocation (inverse of
`n_alloc()`).

Adding a `icc_psu` column to the frame (for example from `varcomp()`
with its `strata` argument) turns the allocation into a stratified
two-stage design, with cost-optimal or fixed cluster takes per stratum.
See `?n_alloc` and the vignette for the full workflow.

## Design effects

`design_effect()` anticipates the design effect of a plan by combining the
design features you are choosing. Components multiply, and each is selected
by the arguments you supply.


``` r
# Clustering alone: 20 households per cluster
design_effect(icc = 0.05, n_per_psu = 20)
#> Design effect (planning)
#> 
#>   clustering   1.9500   icc = 0.05, n_per_psu = 20
#>   ---------------------
#>   overall      1.9500

# Clustering, unequal weighting, and the stratification gain together
frame <- data.frame(
  N = c(50000, 120000), n = c(600, 400), sd = c(12, 20), mean = c(55, 48)
)
deff <- design_effect(icc = 0.05, n_per_psu = 20, strata = frame)
deff
#> Design effect (planning)
#> 
#>   clustering       1.9500   icc = 0.05, n_per_psu = 20
#>   weighting        1.3899   2 strata, n = 1000
#>   stratification   0.9696   2 strata, between-stratum share 0.03038
#>   -------------------------
#>   overall          2.6279   approx. (Kish)

# Use it wherever a deff is expected
n_prop(p = 0.3, moe = 0.05, deff = deff)
#> Sample size for proportion (wald)
#> n = 848 (p = 0.30, moe = 0.050, deff = 2.63)
effective_n(deff, n = 1000)
#> [1] 380.5354

# Or read the features off a plan you already built
design_effect(n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05))
#> Design effect (planning)
#> 
#>   clustering   1.6392   icc = 0.05, n_per_psu = 13.78
#>   ---------------------
#>   overall      1.6392
```

This is a planning tool. To measure the design effect a *collected* sample
actually achieved, use `survey::svymean(..., deff = TRUE)`, which computes
it from the realized weights, strata, and clusters.

## References

Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.

Kish, L. (1965). *Survey Sampling*. Wiley.

Valliant, R., Dever, J. A., and Kreuter, F. (2018).
*Practical Tools for Designing and Weighting Survey Samples*
(2nd ed.). Springer.
