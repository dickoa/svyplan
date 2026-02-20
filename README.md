
# svyplan

Survey sample size determination, optimal allocation, stratification,
and power analysis for R.

## Installation

``` r
# From GitLab
pak::pkg_install("gitlab::dickoa/svyplan")
```

## Sample sizes

``` r
library(svyplan)

# Proportion with margin of error
n_prop(p = 0.3, moe = 0.05)
#> Sample size for proportion (wald)
#> n = 323 (p = 0.30, moe = 0.050)

# Mean with finite population and design effect
n_mean(var = 100, moe = 2, N = 5000, deff = 1.5)
#> Sample size for mean
#> n = 142 (var = 100.00, moe = 2.000, deff = 1.50)
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
#>  name        .n   .binding
#>  stunting     577         
#>  vaccination  485         
#>  anemia      1127 *
```

Per-domain optimization works by adding domain columns to the targets
data frame.

## Multistage cluster designs

``` r
# Optimal 2-stage allocation within a budget
n_cluster(cost = c(500, 50), delta = 0.05, budget = 100000)
#> Optimal 2-stage allocation
#> Stage 1: n1 = 85 | Stage 2: n2 = 14 -> total n = 1160
#> CV = 0.0376, cost = 100000

# CV for a given allocation
cv_cluster(n = c(50, 12), delta = 0.05)
#> [1] 0.0508265
```

Variance components can be estimated from frame data and passed directly
to `n_cluster()`:

``` r
set.seed(1)
frame <- data.frame(
  district = rep(1:40, each = 20),
  income = rep(rnorm(40, 500, 100), each = 20) + rnorm(800, 0, 50)
)

vc <- varcomp(income ~ district, data = frame)
vc
#> Variance components (2-stage)
#> var_between = 0.0307, var_within = 0.0105
#> delta = 0.7448
#> k = 1.0311
#> Unit relvariance = 0.0400

n_cluster(cost = c(500, 50), delta = vc, cv = 0.05)
#> Optimal 2-stage allocation
#> Stage 1: n1 = 15 | Stage 2: n2 = 2 -> total n = 27
#> CV = 0.0500, cost = 8635
```

## Strata boundaries

`strata_bound()` finds optimal boundaries for a continuous
stratification variable.

``` r
set.seed(1)
x <- rlnorm(5000, meanlog = 6, sdlog = 1.2)

strata_bound(x, n_strata = 4, n = 300, method = "cumrootf")
#> Strata boundaries (Dalenius-Hodges, 4 strata)
#> Boundaries: 400.0, 1300.0, 3300.0
#> n = 300, CV = 0.0211
#> ---
#>  stratum lower      upper    N_h  W_h   S_h    n_h
#>  1          4.92557   400.00 2523 0.505 107.1   44
#>  2        400.00000  1300.00 1610 0.322 250.3   66
#>  3       1300.00000  3300.00  647 0.129 544.0   58
#>  4       3300.00000 39039.61  220 0.044 3675.1 132
```

Four methods are available: Dalenius-Hodges (`"cumrootf"`), geometric
(`"geo"`), Lavallee-Hidiroglou (`"lh"`), and Kozak (`"kozak"`).

## Power analysis

Solve for sample size, power, or minimum detectable effect. Supports
design effects, finite population correction, and panel overlap.

``` r
# Sample size to detect a 5pp change from 70% with deff = 2
power_prop(p1 = 0.70, p2 = 0.75, deff = 2.0)
#> Power analysis for proportions (solved for sample size)
#> n = 2496 (per group), power = 0.800, delta = 0.0500
#> (p1 = 0.700, p2 = 0.750, alpha = 0.05, deff = 2.00)

# MDE with n = 1500 per group
power_prop(p1 = 0.70, n = 1500, deff = 2.0)
#> Power analysis for proportions (solved for minimum detectable effect)
#> n = 1500 (per group), power = 0.800, delta = 0.0639
#> (p1 = 0.700, p2 = 0.764, alpha = 0.05, deff = 2.00)

# Means
power_mean(delta = 5, var = 200)
#> Power analysis for means (solved for sample size)
#> n = 126 (per group), power = 0.800, delta = 5.0000
#> (alpha = 0.05)
```

## Design effects

``` r
# Planning: expected cluster design effect
design_effect(delta = 0.05, m = 20, method = "cluster")
#> [1] 1.95

# Diagnostic: Kish design effect from weights
set.seed(1)
w <- runif(500, 0.5, 4)
design_effect(w, method = "kish")
#> [1] 1.196494
effective_n(w, method = "kish")
#> [1] 417.8876
```

## References

Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.

Kish, L. (1965). *Survey Sampling*. Wiley.

Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
Designing and Weighting Survey Samples* (2nd ed.). Springer.
