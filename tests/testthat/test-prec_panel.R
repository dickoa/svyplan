## Q1. The round trip against the target's own precision

lfs_ret <- c(0.878, 0.963, 0.936, 0.956)
lfs_rr <- 0.728

test_that("a panel sized at a wave reports the target's precision there", {
  for (target in list(
    n_mean(var = 100, moe = 2),
    n_prop(p = 0.3, moe = 0.03),
    n_prop(p = 0.2, moe = 0.03, N = 20000, deff = 1.5, method = "wilson"),
    n_prop(p = 0.05, moe = 0.02, df = 12, method = "beta"),
    n_mean(var = 100, cv = 0.05, mu = 40),
    prec_mean(var = 100, n = 1000)
  )) {
    for (w in seq_len(5L)) {
      plan <- n_panel(
        target, retention = lfs_ret, resp_rate = lfs_rr, target_wave = w
      )
      back <- prec_panel(plan)
      expect_equal(back$moe, target$moe, tolerance = 1e-9)
      expect_equal(back$se, target$se, tolerance = 1e-9)
      expect_equal(back$cv, target$cv, tolerance = 1e-9)
      expect_equal(back$n_resp, plan$n_resp, tolerance = 1e-12)
      expect_identical(back$target_wave, plan$target_wave)
    }
  }
})

test_that("a rotating occasion round-trips on the pooled sample", {
  target <- n_prop(p = 0.3, moe = 0.03)
  plan <- n_panel(
    target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  back <- prec_panel(plan)
  expect_equal(back$moe, target$moe, tolerance = 1e-9)
  expect_equal(back$n_entrants, plan$n_entrants, tolerance = 1e-12)
  expect_equal(back$n_in_sample, plan$n_in_sample, tolerance = 1e-12)
})

test_that("only the forward direction claims to have solved for a size", {
  target <- n_mean(var = 100, moe = 2)
  plan <- n_panel(target, retention = 0.9, resp_rate = 0.8)
  expect_identical(plan$solved, "n_recruit")
  expect_null(prec_panel(plan)$solved)
  expect_null(
    prec_panel(500, target, retention = 0.9, resp_rate = 0.8)$solved
  )
  expect_identical(class(prec_panel(plan)), c("svyplan_panel", "list"))
})

## Q2. The direct form reads a recruitment nobody solved for

test_that("the recruitment the target asks for reproduces the target", {
  target <- n_prop(p = 0.3, moe = 0.03)
  plan <- n_panel(target, retention = lfs_ret, resp_rate = lfs_rr)
  direct <- prec_panel(
    plan$n_issued, target, retention = lfs_ret, resp_rate = lfs_rr
  )
  expect_equal(direct$moe, target$moe, tolerance = 1e-9)
  expect_equal(direct$waves$se, plan$waves$se, tolerance = 1e-12)
})

test_that("a smaller recruitment is a wider interval and a shortfall", {
  target <- prec_mean(var = 100, n = 1000)
  short <- prec_panel(1500, target, retention = lfs_ret, resp_rate = lfs_rr)
  expect_lt(short$n_resp, short$n_target)
  expect_gt(short$moe, target$moe)
  expect_equal(short$n_target, 1000, tolerance = 1e-12)
  # The requirement is still reported, so the gap is readable.
  expect_equal(short$n_resp, 1500 * short$waves$q[5L], tolerance = 1e-12)
})

test_that("assurance in the reverse direction is the same requirement", {
  target <- prec_mean(var = 100, n = 1000)
  plan <- n_panel(target, retention = lfs_ret, resp_rate = lfs_rr,
                  assurance = 0.9)
  direct <- prec_panel(
    1500, target, retention = lfs_ret, resp_rate = lfs_rr, assurance = 0.9
  )
  expect_equal(direct$n_assured, plan$n_assured, tolerance = 1e-12)
  expect_gt(direct$n_assured, 1500)
})

test_that("the direct form validates what n_panel validates", {
  target <- prec_mean(var = 100, n = 100)
  expect_error(prec_panel(0, target, retention = 0.9), "single positive size")
  expect_error(prec_panel(-5, target, retention = 0.9), "single positive size")
  expect_error(
    prec_panel(c(100, 200), target, retention = 0.9), "single positive size"
  )
  expect_error(prec_panel(500, target, retention = 2), "conditional retention")
  expect_error(
    prec_panel(500, target, retention = 0.9, design = "rotating",
               target_wave = 2),
    "does not apply to a rotating panel"
  )
  expect_error(
    prec_panel(500, 1000, retention = 0.9), "must be a result from"
  )
  expect_error(
    prec_panel(500, target, retention = 0.9, resp_rte = 0.8),
    "unused argument.*resp_rte"
  )
})

## Q3. Re-reading a stored plan under other assumptions

test_that("named overrides replace the stored assumptions", {
  target <- prec_mean(var = 100, n = 1000)
  plan <- n_panel(target, retention = lfs_ret, resp_rate = lfs_rr)
  worse <- prec_panel(plan, retention = c(0.80, 0.90, 0.90, 0.92))
  expect_lt(worse$n_resp, plan$n_resp)
  expect_gt(worse$moe, plan$moe)
  # The recruitment is what was stored: it is the plan being re-read, not a
  # new one being solved.
  expect_equal(worse$n_issued, plan$n_issued, tolerance = 1e-12)
  expect_equal(
    prec_panel(plan, resp_rate = 0.6)$waves$q[1L], 0.6, tolerance = 1e-12
  )
})

test_that("design cannot be overridden on a stored panel", {
  # The stored count is the whole issue to one cohort. Reading it as entrants
  # per occasion changes the unit the number is in, not an assumption about
  # it, so 1816 issued would silently become 1816 an occasion.
  target <- prec_mean(var = 100, n = 1000)
  fixed <- n_panel(target, retention = lfs_ret, resp_rate = lfs_rr)
  expect_error(
    prec_panel(fixed, design = "rotating"),
    "'design' cannot be overridden"
  )
  expect_error(
    prec_panel(fixed, design = "rotating"), "whole issue to one cohort"
  )
  rot <- n_panel(target, retention = lfs_ret, resp_rate = lfs_rr,
                 design = "rotating")
  expect_error(prec_panel(rot, design = "fixed"), "entrants per occasion")
  # Restating the same design is refused too: the argument is the stored
  # one's, and an override that happens to agree still reads as permission.
  expect_error(prec_panel(fixed, design = "fixed"), "cannot be overridden")
  # The reinterpretation stays available where the count has to be restated.
  direct <- prec_panel(
    fixed$n_issued, target, retention = lfs_ret, resp_rate = lfs_rr,
    design = "rotating"
  )
  expect_identical(direct$design, "rotating")
  expect_equal(direct$n_entrants, fixed$n_issued, tolerance = 1e-12)
})

test_that("an unknown override is an error, not a silent drop", {
  plan <- n_panel(prec_mean(var = 100, n = 100), retention = 0.9)
  expect_error(prec_panel(plan, retenton = 0.9), "unused argument.*retenton")
  expect_error(prec_panel(plan, 0.9), "must be named")
})

## Q4. What a wave row means under each design

test_that("a fixed wave row is that wave's estimate", {
  target <- n_prop(p = 0.3, moe = 0.03, N = 50000, deff = 1.2)
  plan <- n_panel(target, retention = lfs_ret, resp_rate = lfs_rr)
  for (w in seq_len(5L)) {
    expect_equal(
      plan$waves$moe[w],
      prec_prop(p = 0.3, n = plan$waves$n_resp[w], N = 50000,
                deff = 1.2)$moe,
      tolerance = 1e-12
    )
  }
})

test_that("a rotating wave row is one cohort, and the headline is the pool", {
  target <- n_prop(p = 0.3, moe = 0.03)
  plan <- n_panel(
    target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  expect_equal(
    plan$se, prec_prop(p = 0.3, n = plan$n_resp)$se, tolerance = 1e-12
  )
  for (w in seq_len(5L)) {
    expect_equal(
      plan$waves$se[w],
      prec_prop(p = 0.3, n = plan$waves$n_resp[w])$se,
      tolerance = 1e-12
    )
  }
  # A single cohort at a single stage is far wider than the occasion.
  expect_gt(min(plan$waves$se), 2 * plan$se)
})

test_that("the panel applies no response rate of its own on top of attrition", {
  # The losses are already in the wave counts, so a second inflation would
  # count response twice.
  target <- n_mean(var = 100, moe = 2, resp_rate = 0.5)
  plan <- n_panel(target, retention = 0.9, resp_rate = 0.8)
  expect_equal(
    plan$se,
    prec_mean(var = 100, n = plan$n_resp)$se,
    tolerance = 1e-12
  )
  expect_equal(plan$n_resp, target$n * 0.5, tolerance = 1e-9)
})
