## One indicator on one design is one answer, whichever interface states it.
##
## svyplan reaches the same two-stage and three-stage cluster design from four
## directions: n_cluster(), n_multi_cluster(), n_alloc() cluster mode, and
## n_alloc() generalized fixed-take allocation. They share a variance model, so
## a design described identically to all four must size identically in all
## four. Nothing else in the suite compares them, and a coefficient that drifts
## in one of them stays green everywhere else.
##
## The comparisons run at a negligible sampling fraction. That is the regime
## every one of these interfaces claims: n_cluster() applies no finite
## population correction at any stage, while n_alloc() applies the
## ultimate-unit one. Away from negligible f they are answering different
## questions and are not expected to agree. The FPC block at the end pins that
## difference rather than hiding from it.

.eq_cv <- 0.10
.eq_p <- 0.30
.eq_relvar <- (1 - .eq_p) / .eq_p

## A population large enough that 1 - n / N is inert to the tolerances used.
.eq_N <- 2e8
.eq_N_psu <- 1e7

.eq_alloc_frame <- function(icc, m) {
  data.frame(
    stratum = "S", N = .eq_N, mean = .eq_p, sd = sqrt(.eq_p * (1 - .eq_p)),
    icc_psu = icc, n_per_psu = m, N_psu = .eq_N_psu,
    cost_psu = 300, cost_ssu = 25
  )
}

.eq_targets <- data.frame(
  name = "y", domain = ".overall", level = NA, cv = .eq_cv
)

## T1. Two-stage: all four interfaces return the same PSU count

test_that("two-stage cluster interfaces agree at full response", {
  icc <- 0.05
  m <- 20

  by_cluster <- n_cluster(
    icc = icc, n_per_psu = m, cv = .eq_cv,
    unit_relvar = .eq_relvar, stage_cost = c(300, 25)
  )$n[["n_psu"]]

  by_multi <- n_multi_cluster(
    data.frame(name = "y", p = .eq_p, cv = .eq_cv, icc_psu = icc),
    n_per_psu = m, stage_cost = c(300, 25)
  )$n[["n_psu"]]

  by_alloc <- n_alloc(.eq_alloc_frame(icc, m), cv = .eq_cv)$detail$n / m

  by_generalized <- n_alloc(
    data.frame(
      stratum = "S", N = .eq_N, N_psu = .eq_N_psu, n_per_psu = m,
      cost_psu = 300, cost_ssu = 25
    ),
    measures = data.frame(
      stratum = "S", name = "y", p = .eq_p, icc_psu = icc, var_ratio_psu = 1
    ),
    targets = .eq_targets
  )$detail$n / m

  expect_equal(by_multi, by_cluster, tolerance = 1e-6)
  expect_equal(by_alloc, by_cluster, tolerance = 1e-4)
  expect_equal(by_generalized, by_cluster, tolerance = 1e-4)
})

test_that("two-stage cluster interfaces agree below full response", {
  icc <- 0.05
  m <- 20
  rr <- 0.5

  by_cluster <- n_cluster(
    icc = icc, n_per_psu = m, cv = .eq_cv, resp_rate = rr,
    unit_relvar = .eq_relvar, stage_cost = c(300, 25)
  )$n[["n_psu"]]

  by_multi <- n_multi_cluster(
    data.frame(name = "y", p = .eq_p, cv = .eq_cv, icc_psu = icc),
    n_per_psu = m, resp_rate = rr, stage_cost = c(300, 25)
  )$n[["n_psu"]]

  by_alloc <- n_alloc(
    .eq_alloc_frame(icc, m), cv = .eq_cv, resp_rate = rr
  )$detail$n / m

  by_generalized <- n_alloc(
    data.frame(
      stratum = "S", N = .eq_N, N_psu = .eq_N_psu, n_per_psu = m,
      cost_psu = 300, cost_ssu = 25
    ),
    measures = data.frame(
      stratum = "S", name = "y", p = .eq_p, icc_psu = icc,
      var_ratio_psu = 1, resp_rate = rr
    ),
    targets = .eq_targets
  )$detail$n / m

  expect_equal(by_multi, by_cluster, tolerance = 1e-6)
  expect_equal(by_alloc, by_cluster, tolerance = 1e-4)
  expect_equal(by_generalized, by_cluster, tolerance = 1e-4)
})

## T2. Three-stage, where the correction is not a single take substitution

test_that("three-stage cluster interfaces agree below full response", {
  d1 <- 0.05
  d2 <- 0.10
  m <- 5
  q <- 4
  rr <- 0.5

  by_cluster <- n_cluster(
    icc = c(d1, d2), n_per_psu = m, n_per_ssu = q, cv = .eq_cv,
    resp_rate = rr, unit_relvar = .eq_relvar, stage_cost = c(300, 50, 25)
  )$n[["n_psu"]]

  by_generalized <- n_alloc(
    data.frame(
      stratum = "S", N = .eq_N, N_psu = .eq_N_psu, N_ssu = 4e7,
      n_per_psu = m, n_per_ssu = q,
      cost_psu = 300, cost_ssu = 50, cost_tsu = 25
    ),
    measures = data.frame(
      stratum = "S", name = "y", p = .eq_p, icc_psu = d1, icc_ssu = d2,
      var_ratio_psu = 1, resp_rate = rr
    ),
    targets = .eq_targets
  )$detail$n / (m * q)

  expect_equal(by_generalized, by_cluster, tolerance = 1e-4)
})

## T3. Each rate acts at the stage it names, and none stands in for another

test_that("ultimate-unit response enters the clustering bracket", {
  icc <- 0.05
  m <- 20
  rr <- 0.5
  gross <- 1 + icc * (m - 1)
  responding <- 1 + icc * (m * rr - 1)

  n_psu <- n_alloc(
    data.frame(
      stratum = "S", N = .eq_N, N_psu = .eq_N_psu, n_per_psu = m,
      cost_psu = 300, cost_ssu = 25
    ),
    measures = data.frame(
      stratum = "S", name = "y", p = .eq_p, icc_psu = icc,
      var_ratio_psu = 1, resp_rate = rr
    ),
    targets = .eq_targets
  )$detail$n / m

  # a = V * bracket / (cv^2 * m * r), so the bracket is recoverable exactly.
  recovered <- n_psu * .eq_cv^2 * m * rr / .eq_relvar
  expect_equal(recovered, responding, tolerance = 1e-4)
  expect_false(isTRUE(all.equal(recovered, gross, tolerance = 1e-3)))
})

test_that("PSU response scales the design without touching the bracket", {
  icc <- 0.05
  m <- 20
  rr <- 0.5

  size <- function(psu_rate) {
    n_alloc(
      data.frame(
        stratum = "S", N = .eq_N, N_psu = .eq_N_psu, n_per_psu = m,
        cost_psu = 300, cost_ssu = 25
      ),
      measures = data.frame(
        stratum = "S", name = "y", p = .eq_p, icc_psu = icc,
        var_ratio_psu = 1, resp_rate = rr, resp_rate_psu = psu_rate
      ),
      targets = .eq_targets
    )$detail$n / m
  }

  # Losing whole PSUs is a pure sample-size loss, so it divides the whole
  # requirement and leaves the clustering penalty where it was.
  expect_equal(size(0.8), size(1) / 0.8, tolerance = 1e-4)
  expect_equal(
    size(0.8),
    n_cluster(
      icc = icc, n_per_psu = m, cv = .eq_cv, resp_rate = rr,
      resp_rate_psu = 0.8, unit_relvar = .eq_relvar,
      stage_cost = c(300, 25)
    )$n[["n_psu"]],
    tolerance = 1e-4
  )
})

test_that("the two stage rates are not recoverable from one another", {
  icc <- 0.05
  m <- 20

  size <- function(...) {
    n_alloc(
      data.frame(
        stratum = "S", N = .eq_N, N_psu = .eq_N_psu, n_per_psu = m,
        cost_psu = 300, cost_ssu = 25
      ),
      measures = data.frame(
        stratum = "S", name = "y", p = .eq_p, icc_psu = icc,
        var_ratio_psu = 1, ...
      ),
      targets = .eq_targets
    )$detail$n / m
  }

  # Same product of rates, different designs, because only the ultimate-unit
  # rate shrinks the realized cluster.
  expect_false(
    isTRUE(all.equal(
      size(resp_rate = 0.5, resp_rate_psu = 0.8),
      size(resp_rate = 0.8, resp_rate_psu = 0.5),
      tolerance = 1e-3
    ))
  )
})

test_that("a stage rate naming an absent stage is refused", {
  frame_1 <- data.frame(stratum = "S", N = 1e5)
  measures_1 <- data.frame(
    stratum = "S", name = "y", p = .eq_p, resp_rate_psu = 0.8
  )
  expect_error(
    n_alloc(frame_1, measures = measures_1, targets = .eq_targets),
    "does not have"
  )

  frame_2 <- data.frame(
    stratum = "S", N = .eq_N, N_psu = .eq_N_psu, n_per_psu = 20,
    cost_psu = 300, cost_ssu = 25
  )
  measures_2 <- data.frame(
    stratum = "S", name = "y", p = .eq_p, icc_psu = 0.05,
    var_ratio_psu = 1, resp_rate_ssu = 0.9
  )
  expect_error(
    n_alloc(frame_2, measures = measures_2, targets = .eq_targets),
    "does not have"
  )
})

## T4. The reported design effect describes the variance model that was used

test_that("design_effect() of a cluster allocation reads the responding take", {
  frame <- data.frame(
    stratum = "S", N = 2e5, mean = .eq_p, sd = sqrt(.eq_p * (1 - .eq_p)),
    icc_psu = 0.05, n_per_psu = 20, N_psu = 1e4,
    cost_psu = 300, cost_ssu = 25
  )
  expect_equal(
    design_effect(n_alloc(frame, cv = .eq_cv))$deff_allocation,
    1 + 0.05 * (20 - 1)
  )
  expect_equal(
    design_effect(n_alloc(frame, cv = .eq_cv, resp_rate = 0.5))$deff_allocation,
    1 + 0.05 * (20 * 0.5 - 1)
  )
})

## T5. Where the interfaces part company, and by how much
##
## n_alloc() carries the ultimate-unit correction 1 - n / N and n_cluster()
## carries none, so the two separate as the sampling fraction grows. The
## correction is applied to the whole ICC-inflated variance rather than stage
## by stage, which makes it conservative against the exact two-stage
## without-replacement variance on both components. These bounds are the
## claim, and they hold at every fraction rather than only at negligible ones.

test_that("the allocation FPC sits between no correction and the exact one", {
  icc <- 0.05
  m <- 20
  for (N in c(2e5, 5e4, 2e4)) {
    frame <- data.frame(
      stratum = "S", N = N, mean = .eq_p, sd = sqrt(.eq_p * (1 - .eq_p)),
      icc_psu = icc, n_per_psu = m, N_psu = N / m,
      cost_psu = 300, cost_ssu = 25
    )
    got <- n_alloc(frame, cv = .eq_cv)$detail$n

    no_fpc <- .eq_relvar * (1 + icc * (m - 1)) / .eq_cv^2
    expect_lte(got, no_fpc + 1e-6)

    # A lower reference, not the model under test. The frame sets
    # N_psu = N / m, so every PSU holds exactly the take and this is a
    # one-stage sample of whole PSUs. It is also written on the exact
    # finite-population between-PSU variance S1^2 = S^2(1+icc(M-1))/M rather
    # than on the icc decomposition the package plans with, so the two agree
    # only as M grows. What is asserted is the bound, not equality: the
    # default correction must not fall below it. The component-wise oracle
    # for the implemented model is in test-cluster-fpc.R.
    M <- m
    A_pop <- N / M
    a <- got / m
    S2 <- .eq_p * (1 - .eq_p)
    between <- (1 - a / A_pop) * S2 * (1 + icc * (M - 1)) / (a * M)
    within <- (1 - m / M) * S2 * (1 - icc) / (a * m)
    exact_var <- between + within
    got_var <- S2 * (1 + icc * (m - 1)) * (1 / got - 1 / N)
    expect_gte(got_var, exact_var - 1e-12)
  }
})
