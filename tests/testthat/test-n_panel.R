## P1. The UK LFS five-wave validation

# Smith, Lynn and Elliot (2009), chapter 2 section 2.3.3: quarterly rotating
# address panel, wave-1 response 0.728 then conditional retention 0.878,
# 0.963, 0.936, 0.956.
lfs_ret <- c(0.878, 0.963, 0.936, 0.956)
lfs_rr <- 0.728
# A precision result carries an exact size, so the published figures are not
# read against a size that had to be solved first.
lfs_target <- prec_mean(var = 100, n = 1000)

test_that("the cumulative response probabilities are the published ones", {
  q <- .panel_q(lfs_rr, lfs_ret)
  expect_equal(
    q,
    c(0.728000, 0.639184, 0.615534, 0.576140, 0.550790),
    tolerance = 1e-6
  )
  expect_equal(sum(q), 3.109648, tolerance = 1e-6)
})

test_that("a fixed panel issues 1816 to hold 1000 at wave 5", {
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  expect_s3_class(plan, "svyplan_panel")
  expect_equal(1 / plan$waves$q[5L], 1.8156, tolerance = 1e-4)
  expect_equal(plan$n_issued, 1815.58, tolerance = 1e-2)
  expect_equal(plan$n_resp, 1000, tolerance = 1e-9)
  expect_identical(plan$target_wave, 5L)
})

test_that("a rotating panel takes 322 entrants per period and holds 1608", {
  plan <- n_panel(
    lfs_target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  expect_equal(plan$n_entrants, 321.58, tolerance = 1e-2)
  expect_equal(plan$n_in_sample, 1607.90, tolerance = 1e-2)
  expect_identical(plan$n_cohorts, 5L)
  expect_equal(sum(plan$waves$n_resp), 1000, tolerance = 1e-9)
})

test_that("wave 1 carries 61 percent of the whole life's loss", {
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  expect_equal(plan$waves$loss_share[1L], 0.6055, tolerance = 1e-4)
  expect_equal(sum(plan$waves$loss_share), 1, tolerance = 1e-12)
})

## P2. The two recruitment numbers are different quantities

test_that("each design reports its own recruitment under its own name", {
  fixed <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  rot <- n_panel(
    lfs_target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  # A single cohort's issue and the entrants per period differ by a factor
  # of five here, so neither name may hold the other's number.
  expect_null(fixed$n_entrants)
  expect_null(fixed$n_in_sample)
  expect_null(rot$n_issued)
  expect_null(rot$target_wave)
  expect_gt(fixed$n_issued / rot$n_entrants, 5)
})

test_that("the in-sample total is the cohorts a running design already holds", {
  rot <- n_panel(
    lfs_target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  expect_equal(rot$n_in_sample, rot$n_cohorts * rot$n_entrants,
               tolerance = 1e-12)
})

test_that("a rotating period pools its cohorts and a fixed wave does not", {
  rot <- n_panel(
    lfs_target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  expect_equal(rot$n_resp, sum(rot$waves$n_resp), tolerance = 1e-12)
  # The pooled period is the target; one cohort at one stage is a fifth of
  # it, which is the reading the plan warns against.
  expect_lt(max(rot$waves$n_resp), rot$n_resp / 4)
})

## P3. A panel of one loss is the response inflation the package already has

test_that("full retention reduces to the resp_rate inflation of n_prop", {
  # q = (r, r), so the requirement at wave 2 is n_target / r, which is what
  # resp_rate does on a single period.
  plan <- n_panel(
    n_prop(p = 0.3, moe = 0.03), retention = 1, resp_rate = 0.6
  )
  expect_equal(
    plan$n_issued,
    n_prop(p = 0.3, moe = 0.03, resp_rate = 0.6)$n,
    tolerance = 1e-12
  )
})

test_that("no loss anywhere leaves the target untouched", {
  plan <- n_panel(n_mean(var = 100, moe = 2), retention = c(1, 1))
  expect_equal(plan$n_issued, n_mean(var = 100, moe = 2)$n, tolerance = 1e-12)
  expect_true(all(is.na(plan$waves$loss_share)))
})

## P4. Normalizing the target's own response rate is exact

test_that("a target's resp_rate is removed exactly, not approximately", {
  ret <- c(0.9, 0.85)
  cases <- list(
    plain = list(n_prop(p = 0.3, moe = 0.03),
                 n_prop(p = 0.3, moe = 0.03, resp_rate = 0.7)),
    finite_N = list(n_prop(p = 0.3, moe = 0.03, N = 5000),
                    n_prop(p = 0.3, moe = 0.03, N = 5000, resp_rate = 0.7)),
    wilson = list(n_prop(p = 0.3, moe = 0.03, N = 5000, method = "wilson"),
                  n_prop(p = 0.3, moe = 0.03, N = 5000, method = "wilson",
                         resp_rate = 0.7)),
    min_cases = list(n_prop(p = 0.02, moe = 0.05, min_cases = 50),
                     n_prop(p = 0.02, moe = 0.05, min_cases = 50,
                            resp_rate = 0.7)),
    mean = list(n_mean(var = 100, moe = 2),
                n_mean(var = 100, moe = 2, resp_rate = 0.7))
  )
  for (nm in names(cases)) {
    full <- n_panel(cases[[nm]][[1L]], retention = ret, resp_rate = 0.8)
    # 0.7 against the panel's 0.8 is the disagreement that gets reported.
    expect_warning(
      netted <- n_panel(cases[[nm]][[2L]], retention = ret, resp_rate = 0.8),
      "response rate 0.7 is not used"
    )
    expect_equal(netted$n_issued, full$n_issued, tolerance = 1e-12,
                 info = nm)
    expect_equal(netted$moe, full$moe, tolerance = 1e-12, info = nm)
  }
})

test_that("a response rate the panel does not use is reported, once named", {
  ret <- c(0.9, 0.85)
  # Two different recruitment responses named, only one used.
  expect_warning(
    plan <- n_panel(n_prop(p = 0.3, moe = 0.03, resp_rate = 0.5),
                    retention = ret, resp_rate = 0.8),
    "the target's own response rate 0.5 is not used"
  )
  # The panel's own rate is named too, so the message says which one won.
  expect_warning(
    n_panel(n_prop(p = 0.3, moe = 0.03, resp_rate = 0.5),
            retention = ret, resp_rate = 0.8),
    "'resp_rate' \\(0.8\\) at wave 1"
  )
  # A default resp_rate of 1 is still a rate the target's does not survive.
  expect_warning(
    n_panel(n_prop(p = 0.3, moe = 0.03, resp_rate = 0.7), retention = ret),
    "is not used"
  )
  # The other direction shares the assembly, so it shares the report.
  expect_warning(
    prec_panel(1000, n_prop(p = 0.3, moe = 0.03, resp_rate = 0.5),
               retention = ret, resp_rate = 0.8),
    "is not used"
  )
  # The warning is gone by the time the object is read, so summary() carries
  # it. print() does not: the plan itself is correct.
  sm <- capture.output(print(summary(plan)))
  expect_true(any(grepl("^target: response 0.5 removed, requirement [0-9]+ responding$",
                        sm)))
  expect_false(any(grepl("^target:", capture.output(print(plan)))))
})

test_that("a target's eligibility rate is named apart from its response", {
  target <- n_prop(p = 0.3, moe = 0.03, N = 2000, resp_rate = 0.8,
                   eligibility_rate = 0.5)
  expect_warning(
    plan <- n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8),
    "own response 0.8 times eligibility 0.5 \\(0.4\\) is not used"
  )
  sm <- capture.output(print(summary(plan)))
  expect_true(any(grepl(
    "^target: response 0.8 times eligibility 0.5 \\(0.4\\) removed", sm
  )))
  expect_no_warning(n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8,
                            eligibility_rate = 0.5))
})

test_that("a response rate the panel agrees with is not reported", {
  ret <- c(0.9, 0.85)
  # The panel removes 0.8 and re-applies 0.8 at wave 1: nothing to say.
  expect_silent(
    plan <- n_panel(n_prop(p = 0.3, moe = 0.03, resp_rate = 0.8),
                    retention = ret, resp_rate = 0.8)
  )
  expect_false(any(grepl("^target:", capture.output(print(summary(plan))))))
  # A target with no rate of its own, and one that states the null rate.
  expect_silent(
    n_panel(n_prop(p = 0.3, moe = 0.03), retention = ret, resp_rate = 0.8)
  )
  expect_silent(
    n_panel(n_prop(p = 0.3, moe = 0.03, resp_rate = 1), retention = ret,
            resp_rate = 0.8)
  )
})

test_that("the warning and the disclosure cannot drift apart", {
  # One helper decides both, so every pair of rates agrees on whether it is
  # a conflict, including the floating-point near-misses.
  expect_null(.panel_rate_conflict(NULL, 0.8))
  expect_null(.panel_rate_conflict(1, 0.8))
  expect_null(.panel_rate_conflict(0.8, 0.8))
  expect_null(.panel_rate_conflict(0.1 + 0.7, 0.8))
  expect_equal(.panel_rate_conflict(0.5, 0.8), 0.5)
  expect_equal(.panel_rate_conflict(0.7, 1), 0.7)
})

test_that("a min_cases target keeps its floor at the target wave", {
  target <- n_prop(p = 0.02, moe = 0.05, min_cases = 50)
  expect_identical(target$binding, "min_cases")
  plan <- n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8)
  expect_equal(plan$n_resp * 0.02, 50, tolerance = 1e-9)
})

test_that("expected cases are reported by wave, for a proportion only", {
  target <- n_prop(p = 0.02, moe = 0.05, min_cases = 50)
  plan <- n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8,
                  target_wave = 1)
  expect_equal(plan$waves$expected_cases, plan$waves$n_resp * 0.02,
               tolerance = 1e-12)
  # Sized to wave 1, the floor holds there and not after it, which is what
  # the column makes readable.
  expect_equal(plan$waves$expected_cases[1L], 50, tolerance = 1e-9)
  expect_true(all(plan$waves$expected_cases[-1L] < 50))
  expect_null(
    n_panel(n_mean(var = 100, moe = 2), retention = 0.9)$waves$expected_cases
  )
})

## P5. What a target may be

test_that("a target must be a single-indicator mean or proportion result", {
  expect_error(n_panel(1000, retention = 0.9), "n_mean\\(\\)")
  expect_error(
    n_panel(n_change(var = 100, moe = 2), retention = 0.9),
    "two occasions"
  )
  expect_error(
    n_panel(n_cluster(c(500, 50), icc = 0.05, cv = 0.05), retention = 0.9),
    "must be a result from"
  )
  expect_error(
    n_panel(
      n_multi(data.frame(p = c(0.3, 0.5), moe = c(0.05, 0.05))),
      retention = 0.9
    ),
    "not 'multi'"
  )
  expect_error(
    n_panel(
      n_alloc(data.frame(N = c(100, 200), sd = c(10, 15)), n = 50),
      retention = 0.9
    ),
    "not 'alloc'"
  )
})

test_that("both directions of the single-period families are accepted", {
  for (target in list(
    n_mean(var = 100, moe = 2), n_prop(p = 0.3, moe = 0.03),
    prec_mean(var = 100, n = 400), prec_prop(p = 0.3, n = 400)
  )) {
    plan <- n_panel(target, retention = 0.9, resp_rate = 0.8)
    expect_s3_class(plan, "svyplan_panel")
    expect_equal(plan$n_issued * 0.8 * 0.9, plan$n_resp, tolerance = 1e-12)
  }
})

test_that("a domain allocation is refused before its size is read", {
  target <- n_prop(p = 0.3, moe = 0.03)
  target$domains <- data.frame(domain = "a", n = 100)
  expect_error(n_panel(target, retention = 0.9), "single indicator")
})

## P6. Waves, the target wave, and the estimand carried into each

test_that("the wave table describes a unit's whole life", {
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  expect_identical(nrow(plan$waves), 5L)
  expect_identical(plan$waves$wave, 1:5)
  expect_true(is.na(plan$waves$retention[1L]))
  expect_equal(plan$waves$retention[-1L], lfs_ret, tolerance = 1e-12)
  expect_true(all(diff(plan$waves$q) < 0))
  expect_true(all(diff(plan$waves$n_resp) < 0))
  # Fewer respondents is a wider interval, at every wave.
  expect_true(all(diff(plan$waves$moe) > 0))
})

test_that("target_wave selects the wave the target is met at", {
  plan <- n_panel(
    lfs_target, retention = lfs_ret, resp_rate = lfs_rr, target_wave = 3
  )
  expect_identical(plan$target_wave, 3L)
  expect_equal(plan$waves$n_resp[3L], 1000, tolerance = 1e-9)
  # The waves after it are still reported, and they carry fewer units.
  expect_lt(plan$waves$n_resp[5L], 1000)
  expect_gt(plan$waves$n_resp[1L], 1000)
  expect_lt(
    plan$n_issued,
    n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)$n_issued
  )
})

test_that("a rotating design has no wave to select", {
  expect_error(
    n_panel(lfs_target, retention = lfs_ret, design = "rotating",
            target_wave = 3),
    "does not apply to a rotating panel"
  )
  expect_error(
    n_panel(lfs_target, retention = lfs_ret, target_wave = 6),
    "whole number in 1:5"
  )
  expect_error(
    n_panel(lfs_target, retention = lfs_ret, target_wave = 2.5),
    "whole number in 1:5"
  )
})

test_that("the embedded estimand and its design reach every wave", {
  target <- n_prop(p = 0.2, moe = 0.03, N = 5000, deff = 1.5, df = 12,
                   method = "wilson")
  plan <- n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8)
  expect_identical(plan$type, "proportion")
  expect_identical(plan$method, "wilson")
  for (w in seq_len(3L)) {
    direct <- prec_prop(
      p = 0.2, n = plan$waves$n_resp[w], N = 5000, deff = 1.5, df = 12,
      method = "wilson"
    )
    expect_equal(plan$waves$se[w], direct$se, tolerance = 1e-12)
    expect_equal(plan$waves$moe[w], direct$moe, tolerance = 1e-12)
    expect_equal(plan$waves$cv[w], direct$cv, tolerance = 1e-12)
  }
})

test_that("the target's own precision statement never reverses the evaluator", {
  # 'cv' means a target here and an estimand to solve for in prec_mean(), so
  # a stored target arriving in the evaluator would return a solved mean.
  target <- n_mean(var = 100, cv = 0.05, mu = 40)
  plan <- n_panel(target, retention = 0.9, resp_rate = 0.8)
  expect_equal(
    plan$se,
    prec_mean(var = 100, n = plan$n_resp, mu = 40)$se,
    tolerance = 1e-12
  )
  expect_equal(plan$cv, 0.05, tolerance = 1e-9)
})

## P7. Assurance

test_that("a fixed panel inverts one binomial tail", {
  target <- prec_mean(var = 100, n = 100)
  plan <- n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8,
                  assurance = 0.95)
  q_last <- plan$waves$q[3L]
  expect_gt(plan$n_assured, plan$n_issued)
  expect_gte(stats::pbinom(99, plan$n_assured, q_last, lower.tail = FALSE), 0.95)
  expect_lt(
    stats::pbinom(99, plan$n_assured - 1, q_last, lower.tail = FALSE), 0.95
  )
})

test_that("assurance rises with the level and vanishes without one", {
  target <- prec_mean(var = 100, n = 100)
  args <- list(target, retention = c(0.9, 0.9), resp_rate = 0.8)
  assured <- vapply(
    c(0.5, 0.8, 0.95, 0.99),
    function(l) do.call(n_panel, c(args, assurance = l))$n_assured,
    numeric(1L)
  )
  expect_true(all(diff(assured) >= 0))
  expect_null(do.call(n_panel, args)$n_assured)
})

test_that("the assured period clears its target and the count below does not", {
  q <- .panel_q(lfs_rr, lfs_ret)
  got <- .panel_assure_rotating(1000, q, 0.95)
  # Exact by convolution, computed here independently of the search.
  tail_at <- function(e) {
    pmf <- 1
    for (qi in q) pmf <- .conv_pmf(pmf, stats::dbinom(0:e, e, qi))
    sum(pmf[1001:length(pmf)])
  }
  expect_gte(tail_at(got), 0.95)
  expect_lt(tail_at(got - 1), 0.95)
})

test_that("a rotating period is Poisson-binomial, not one binomial", {
  # The LFS rates cannot show this: their cumulative probabilities sit
  # between 0.55 and 0.73, close enough that a binomial matched on their
  # mean returns the same integer at every level. Severe attrition
  # separates the two, the sum of binomials having the smaller variance by
  # Jensen, so the discriminating case is the one tested.
  q <- .panel_q(0.8, c(0.5, 0.25))
  expect_equal(q, c(0.8, 0.4, 0.1), tolerance = 1e-12)
  mean_matched <- function(need, level) {
    e <- 0L
    while (stats::pbinom(need - 1L, length(q) * e, mean(q),
                         lower.tail = FALSE) < level) {
      e <- e + 1L
    }
    as.double(e)
  }
  expect_equal(.panel_assure_rotating(1000, q, 0.95), 794)
  expect_equal(mean_matched(1000, 0.95), 800)
  expect_equal(.panel_assure_rotating(1000, q, 0.99), 805)
  expect_equal(mean_matched(1000, 0.99), 813)
  # It reaches the reported field through the argument, not only the helper.
  plan <- n_panel(
    prec_mean(var = 100, n = 1000), retention = c(0.5, 0.25), resp_rate = 0.8,
    design = "rotating", assurance = 0.95
  )
  expect_equal(plan$n_assured, 794)
  expect_equal(plan$n_entrants, 1000 / sum(q), tolerance = 1e-12)
})

test_that("equal cumulative probabilities collapse onto the binomial", {
  # Retention of 1 makes every live cohort identical, so the period is
  # Binomial(k * e, r) and the exact answer is known independently.
  q <- .panel_q(0.6, c(1, 1, 1))
  got <- .panel_assure_rotating(97, q, 0.95)
  expect_gte(stats::pbinom(96, 4 * got, 0.6, lower.tail = FALSE), 0.95)
  expect_lt(stats::pbinom(96, 4 * (got - 1), 0.6, lower.tail = FALSE), 0.95)
})

test_that("an assurance level must be a probability", {
  target <- prec_mean(var = 100, n = 100)
  for (bad in list(0, 1, -0.1, 1.5, c(0.9, 0.95), NA_real_, "0.95")) {
    expect_error(
      n_panel(target, retention = 0.9, assurance = bad),
      "probability in \\(0, 1\\)"
    )
  }
})

## P8. Input validation

test_that("retention is one conditional rate per transition", {
  target <- prec_mean(var = 100, n = 100)
  for (bad in list(numeric(0), 0, -0.5, 1.2, c(0.9, 0), NA_real_, "0.9",
                   c(0.9, Inf))) {
    expect_error(
      n_panel(target, retention = bad),
      "one conditional retention per wave transition"
    )
  }
  expect_silent(n_panel(target, retention = c(1, 0.5, 1)))
})

test_that("the recruitment response is one rate in (0, 1]", {
  target <- prec_mean(var = 100, n = 100)
  expect_error(n_panel(target, retention = 0.9, resp_rate = 0), "resp_rate")
  expect_error(n_panel(target, retention = 0.9, resp_rate = 1.2), "resp_rate")
  expect_error(
    n_panel(target, retention = 0.9, resp_rate = c(0.8, 0.9)), "resp_rate"
  )
})

test_that("design names the panel type and refuses a design object", {
  target <- prec_mean(var = 100, n = 100)
  expect_error(
    n_panel(target, retention = 0.9, design = "rotational"), "'arg'"
  )
  fake <- structure(list(), class = "survey.design")
  expect_error(
    n_panel(target, retention = 0.9, design = fake),
    "does not take a design object"
  )
  expect_error(
    n_panel(target, retention = 0.9, design = svyplan(deff = 1.5)),
    "does not take a design object"
  )
})

test_that("unused arguments are rejected", {
  target <- prec_mean(var = 100, n = 100)
  expect_error(
    n_panel(target, retention = 0.9, resp_rte = 0.8),
    "unused argument.*resp_rte"
  )
})

## P9. A panel the population cannot supply

test_that("a cohort larger than its population is refused", {
  # 1000 respondents at wave 3 need 1736 issued, more than the frame.
  expect_error(
    n_panel(
      prec_mean(var = 100, n = 1000, N = 1500),
      retention = c(0.9, 0.8), resp_rate = 0.8
    ),
    "more than the population"
  )
})

test_that("cohorts that cannot be in sample at once are refused", {
  # A rotating design needs fewer entrants but holds three cohorts at a
  # time, so the frame it outgrows is the one it fills at one period.
  expect_error(
    n_panel(
      prec_mean(var = 100, n = 1000, N = 1200),
      retention = c(0.9, 0.8), resp_rate = 0.8, design = "rotating"
    ),
    "cohorts alive at one period"
  )
  expect_s3_class(
    n_panel(
      prec_mean(var = 100, n = 1000, N = 1500),
      retention = c(0.9, 0.8), resp_rate = 0.8, design = "rotating"
    ),
    "svyplan_panel"
  )
})

test_that("a whole-unit design the continuous one fits is reported", {
  # 100 respondents over three waves at full response and retention: 33.33
  # entrants per period, which three equal cohorts of 34 cannot field
  # inside a population of 101.
  target <- prec_mean(var = 100, n = 100, N = 101)
  expect_warning(
    plan <- n_panel(target, retention = c(1, 1), design = "rotating"),
    "102 units in sample against a population of 101"
  )
  expect_equal(plan$n_entrants, 100 / 3, tolerance = 1e-12)
  expect_equal(plan$n_cohorts * ceiling(plan$n_entrants), 102)
  # The continuous plan fits exactly, so it is reported and not refused.
  expect_equal(plan$n_cohorts * plan$n_entrants, 100, tolerance = 1e-9)
  expect_lt(plan$n_cohorts * plan$n_entrants, 101)
})

test_that("an assured recruitment beyond the frame is reported in its own right", {
  # The expected design fits and the assured one does not, so the assurance
  # level is the thing that does not hold.
  target <- prec_mean(var = 100, n = 540, N = 1000)
  expect_warning(
    plan <- n_panel(target, retention = c(0.9, 0.8), resp_rate = 0.8,
                    assurance = 0.999),
    "the level is unattainable"
  )
  expect_lt(plan$n_issued, 1000)
  expect_gt(plan$n_assured, 1000)
  # Its own message, not the rounding one: there the continuous design fits
  # and only its rounding does not, and saying so of a figure past N would
  # contradict itself.
  w <- tryCatch(
    n_panel(target, retention = c(0.9, 0.8), resp_rate = 0.8,
            assurance = 0.999),
    warning = conditionMessage
  )
  expect_match(w, "more than the population of 1000")
  expect_false(grepl("fits", w, fixed = TRUE))
  # A warning does not survive into the object, so the state is a field: the
  # level is unattainable, and the recruitment it would take is still what a
  # planner needs in order to argue for a bigger frame.
  expect_false(plan$assured_feasible)
  out <- capture.output(print(plan))
  expect_true(any(grepl("beyond the population of 1000", out, fixed = TRUE)))
  # The plan itself stands: assurance is an addition to it.
  expect_equal(plan$n_resp, 540, tolerance = 1e-9)
  expect_equal(plan$moe, target$moe, tolerance = 1e-9)
})

test_that("an attainable level is marked attainable, and only when asked for", {
  plan <- n_panel(prec_mean(var = 100, n = 100, N = 100000),
                  retention = 0.9, resp_rate = 0.8, assurance = 0.95)
  expect_true(plan$assured_feasible)
  expect_null(
    n_panel(prec_mean(var = 100, n = 100), retention = 0.9)$assured_feasible
  )
  # No finite frame to exceed is attainable by construction.
  expect_true(
    n_panel(prec_mean(var = 100, n = 100), retention = 0.9,
            assurance = 0.95)$assured_feasible
  )
})

test_that("the assured recruitment is the smallest that clears the level", {
  # The fixed search is the shared one, and it is minimal from both sides:
  # a level below a half puts the answer under the expected recruitment.
  target <- prec_mean(var = 100, n = 100)
  for (lvl in c(0.1, 0.3, 0.5, 0.8, 0.95)) {
    plan <- n_panel(target, retention = 0.8, resp_rate = 0.8, assurance = lvl)
    g <- plan$n_assured
    q_w <- plan$waves$q[2L]
    expect_gte(stats::pbinom(99, g, q_w, lower.tail = FALSE), lvl)
    expect_lt(stats::pbinom(99, g - 1, q_w, lower.tail = FALSE), lvl)
  }
  expect_equal(
    n_panel(target, retention = 0.8, resp_rate = 0.8,
            assurance = 0.1)$n_assured,
    144
  )
  # Below a half the assured recruitment sits under the expected one.
  expect_lt(
    n_panel(target, retention = 0.8, resp_rate = 0.8,
            assurance = 0.1)$n_assured,
    ceiling(n_panel(target, retention = 0.8, resp_rate = 0.8)$n_issued)
  )
})

test_that("the rounding gap keeps its own message", {
  # The continuous design does fit here, which is the whole point of the
  # distinction, so this message says so.
  w <- tryCatch(
    n_panel(prec_mean(var = 100, n = 100, N = 101), retention = c(1, 1),
            design = "rotating"),
    warning = conditionMessage
  )
  expect_match(w, "though the continuous design of 33.3 fits")
  expect_false(grepl("unattainable", w, fixed = TRUE))
})

test_that("a level that rounds towards 1 prints as itself", {
  # The API refuses an assurance of exactly 1, so printing 0.999 as 1.00
  # would name a level it rejects.
  plan <- n_panel(prec_mean(var = 100, n = 540, N = 1e6),
                  retention = c(0.9, 0.8), resp_rate = 0.8, assurance = 0.999)
  out <- capture.output(print(plan))
  expect_true(any(grepl("assured (0.999)", out, fixed = TRUE)))
  expect_false(any(grepl("assured (1.00)", out, fixed = TRUE)))
  # An ordinary level is unchanged.
  expect_true(any(grepl(
    "assured (0.95)",
    capture.output(print(n_panel(prec_mean(var = 100, n = 100),
                                 retention = 0.9, assurance = 0.95))),
    fixed = TRUE
  )))
})

test_that("a whole-unit design inside the frame reports nothing", {
  expect_silent(
    n_panel(prec_mean(var = 100, n = 90, N = 200), retention = c(1, 1),
            design = "rotating")
  )
  expect_silent(
    n_panel(prec_mean(var = 100, n = 100), retention = c(1, 1),
            design = "rotating", assurance = 0.95)
  )
})

## P9b. Waves that field nobody

test_that("waves expecting less than one respondent are named", {
  expect_warning(
    plan <- n_panel(
      prec_mean(var = 100, n = 10), retention = c(0.01, 0.01),
      resp_rate = 0.5, target_wave = 1
    ),
    "waves 2, 3 expect less than one respondent"
  )
  expect_true(all(plan$waves$n_resp[2:3] < 1))
  # The calculation is left alone: planning stays continuous throughout.
  expect_true(all(is.finite(plan$waves$se)))
  expect_warning(
    n_panel(prec_mean(var = 100, n = 10), retention = 0.05,
            resp_rate = 0.5, target_wave = 1),
    "wave 2 expects less than one respondent"
  )
})

test_that("a finite population the panel fits inside is planned normally", {
  target <- prec_mean(var = 100, n = 400, N = 20000)
  plan <- n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8)
  expect_lt(plan$n_issued, 20000)
  expect_equal(
    plan$se, prec_mean(var = 100, n = 400, N = 20000)$se, tolerance = 1e-12
  )
})

## P10. Print and coercion

test_that("print leads with the recruitment and what it holds", {
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                  assurance = 0.95)
  out <- capture.output(print(plan))
  expect_match(out[1L], "fixed, 5-wave life")
  expect_match(out[2L], "^issued: 1816 -> 1000 responding at wave 5$")
  expect_true(any(grepl("assured \\(0.95\\)", out)))
  expect_length(grep("^ *[1-5] ", out), 5L)
})

test_that("print carries the answer and nothing that is derived from it", {
  # The rates, the standing sample and the start-up are all implied by the
  # recruitment and the table, so they belong to summary() and the printed
  # block stays the size of the answer.
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                  design = "rotating", start = "immediate")
  out <- capture.output(print(plan))
  # Two blank lines set the wave table off from the header and the pointer.
  expect_length(out, 12L)
  expect_false(any(grepl("loss at wave 1", out, fixed = TRUE)))
  expect_false(any(grepl("in sample:", out, fixed = TRUE)))
  expect_false(any(grepl("start-up (immediate)", out, fixed = TRUE)))
  # The brief table drops the cumulative rate and the cv it implies.
  expect_match(out[5L], "^ wave retention n_resp se +moe *$")
  expect_identical(out[length(out)],
                   "# summary() for the start-up, the loss and per-wave cv")
})

test_that("print names the entrants and summary the standing sample", {
  rot <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                 design = "rotating")
  out <- capture.output(print(rot))
  expect_match(out[2L],
               "^entrants: 322 per period -> .* pooled over 5 cohorts$")
  sm <- capture.output(print(summary(rot)))
  expect_true(any(grepl("^in sample: 1610 across 5 live cohorts$", sm)))
  expect_true(any(grepl("Waves alive at one period", sm, fixed = TRUE)))
  expect_true(any(grepl("61% of the loss at wave 1", sm, fixed = TRUE)))
})

test_that("summary carries the start-up and the cohort composition", {
  rot <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                 design = "rotating", start = "immediate")
  sm <- summary(rot)
  expect_s3_class(sm, "summary.svyplan_panel")
  expect_identical(sm$plan, rot)
  expect_identical(nrow(sm$periods), nrow(rot$periods))
  expect_identical(sm$composition, rot$period_visits)
  expect_true(any(grepl("By period, immediate start",
                        capture.output(print(sm)))))
  # A fixed design has no equilibrium to reach, so it has no period table.
  fixed <- summary(n_panel(lfs_target, retention = lfs_ret,
                           resp_rate = lfs_rr))
  expect_null(fixed$periods)
  expect_null(fixed$composition)
  expect_false(any(grepl("By period", capture.output(print(fixed)))))
})

test_that("summary counts agree with the counts print shows", {
  # Both blocks are on the whole-unit display path, so a reader cannot see
  # two different recruitments for one plan.
  for (design in c("fixed", "rotating")) {
    plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                    design = design)
    sm <- summary(plan)
    expect_identical(sm$recruit, ceiling(as.double(plan)))
    brief <- svyplan:::.fmt_panel_waves(plan, brief = TRUE)
    expect_identical(sm$waves[names(brief)], brief)
  }
  # The expected cases sit on that same path, not on the continuous one, and
  # only a proportion has them at all.
  prop <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = lfs_ret,
                  resp_rate = lfs_rr)
  sm <- summary(prop)
  expect_identical(sm$waves$cases,
                   format(round(prop$waves$expected_cases *
                                  sm$recruit / as.double(prop))))
  expect_null(summary(
    n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  )$waves$cases)
})

test_that("a recruitment short of the target says so", {
  out <- capture.output(print(
    prec_panel(1500, lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  ))
  expect_match(out[2L], "short of the 1000 the target needs")
})

test_that("the coercions return the recruitment, not the analysis sample", {
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                  assurance = 0.95)
  expect_equal(as.double(plan), plan$n_issued, tolerance = 1e-12)
  expect_identical(as.integer(plan), 1816L)
  # An assurance level was set, and the coercion still reports the expected
  # design: which number to field stays the planner's call.
  expect_gt(plan$n_assured, as.double(plan))
  expect_match(format(plan), "^svyplan_panel \\[fixed, 5 waves")
  expect_identical(as.data.frame(plan), plan$waves)
  expect_error(as.integer(plan, digits = 2), "unused argument")
})

test_that("a rotating result coerces to its entrants", {
  rot <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                 design = "rotating")
  expect_equal(as.double(rot), rot$n_entrants, tolerance = 1e-12)
  expect_identical(as.integer(rot), 322L)
})

test_that("a named retention leaves no name on any reported number", {
  # A name on a numeric survives arithmetic, and one riding a wave label into
  # $se or $moe would reappear wherever that number was passed on.
  plan <- n_panel(
    n_mean(var = 100, moe = 2), retention = c(w2 = 0.9, w3 = 0.8),
    resp_rate = 0.8
  )
  for (field in c("n_issued", "n_target", "n_resp", "se", "moe", "cv")) {
    expect_null(names(plan[[field]]), info = field)
  }
  for (col in c("q", "n_resp", "se", "moe")) {
    expect_null(names(plan$waves[[col]]), info = col)
  }
  expect_null(names(as.double(plan)))
})

test_that("svyplan_panel is a sibling of svyplan_n, not a subtype", {
  plan <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr)
  expect_identical(class(plan), c("svyplan_panel", "list"))
  expect_false(inherits(plan, "svyplan_n"))
  expect_false(inherits(plan, "svyplan_prec"))
  # Nothing registered for an analysis sample answers for a recruitment: each
  # would read the issued count as the units carrying the precision.
  expect_error(prec_mean(plan))
  expect_error(prec_change(plan))
  expect_error(effective_n(plan))
  expect_error(confint(plan))
  expect_error(predict(plan))
})

## N-startup. What a rotating design delivers before it reaches equilibrium

.lfs <- function(...) {
  n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(0.878, 0.963, 0.936, 0.956),
          resp_rate = 0.728, design = "rotating", ...)
}

test_that("start = NULL is inert", {
  # the default must add nothing: no field, no print line, no changed number
  a <- .lfs()
  b <- .lfs(start = "immediate")
  expect_null(a$start)
  expect_null(a$periods)
  expect_null(a$period_visits)
  expect_null(a$params$start)
  expect_equal(a$n_entrants, b$n_entrants, tolerance = 1e-12)
  expect_equal(a$n_resp, b$n_resp, tolerance = 1e-12)
  expect_equal(a$moe, b$moe, tolerance = 1e-12)
  expect_equal(a$waves, b$waves)
  # The start-up reads in summary(), so that is where its absence shows too.
  expect_false(any(grepl("By period", capture.output(print(summary(a))))))
  expect_true(any(grepl("By period, immediate start",
                        capture.output(print(summary(b))))))
  # Neither prints a start-up line: print() carries the plan, not its approach.
  expect_identical(capture.output(print(a)), capture.output(print(b)))
})

test_that("the immediate start follows the closed form", {
  x <- .lfs(start = "immediate")
  e <- x$n_entrants
  q <- x$waves$q
  k <- length(q)
  expected <- vapply(seq_len(k), function(t) {
    e * ((k - t + 1) * q[[t]] + sum(q[seq_len(t - 1L)]))
  }, numeric(1L))
  expect_equal(x$periods$n_resp[seq_len(k)], expected, tolerance = 1e-12)
  # period 1 is one binomial: every unit is at wave 1
  expect_equal(x$periods$n_resp[[1L]], k * e * q[[1L]], tolerance = 1e-12)
})

test_that("the gradual start is the cumulative recruitment", {
  x <- .lfs(start = "gradual")
  e <- x$n_entrants
  q <- x$waves$q
  expect_equal(
    x$periods$n_resp[seq_along(q)],
    vapply(seq_along(q), function(t) e * sum(q[seq_len(t)]), numeric(1L)),
    tolerance = 1e-12
  )
  expect_equal(x$periods$n_entrants, rep(e, length(q) + 1L), tolerance = 1e-12)
})

test_that("the start-up is stored continuous and lands on n_target", {
  # the stored table must not carry whole units: they belong to print()
  for (s in c("gradual", "immediate")) {
    x <- .lfs(start = s)
    last <- nrow(x$periods)
    expect_equal(x$periods$n_resp[[last]], x$n_target, tolerance = 1e-9)
    expect_equal(x$periods$n_resp[[last]], x$n_entrants * sum(x$waves$q),
                 tolerance = 1e-9)
    expect_equal(x$periods$moe[[last]], x$moe, tolerance = 1e-12)
    expect_false(isTRUE(all.equal(x$periods$n_entrants[[last]],
                                  ceiling(x$n_entrants))))
  }
})

test_that("the two starts bracket the equilibrium in opposite directions", {
  imm <- .lfs(start = "immediate")$periods
  grad <- .lfs(start = "gradual")$periods
  at_eq <- tail(imm$n_resp, 1L)
  # an immediate start opens above its own equilibrium, every unit being at
  # wave 1, and settles down onto it; a gradual one opens below and climbs
  expect_gt(imm$n_resp[[1L]], at_eq)
  expect_lt(grad$n_resp[[1L]], at_eq)
  expect_false(is.unsorted(rev(imm$n_resp)))
  expect_false(is.unsorted(grad$n_resp))
  expect_lt(imm$moe[[1L]], tail(imm$moe, 1L))
  # and an immediate start holds the whole design from its first period
  expect_equal(imm$n_in_sample, rep(imm$n_in_sample[[1L]], nrow(imm)),
               tolerance = 1e-12)
  expect_lt(grad$n_in_sample[[1L]], tail(grad$n_in_sample, 1L))
})

test_that("equilibrium is flagged where the two starts meet", {
  imm <- .lfs(start = "immediate")
  grad <- .lfs(start = "gradual")
  k <- nrow(imm$waves)
  expect_equal(which(imm$periods$equilibrium)[[1L]], k)
  expect_equal(which(grad$periods$equilibrium)[[1L]], k)
  expect_equal(imm$periods$n_resp[imm$periods$equilibrium],
               grad$periods$n_resp[grad$periods$equilibrium], tolerance = 1e-12)
})

test_that("period_visits decomposes each period", {
  for (s in c("gradual", "immediate")) {
    x <- .lfs(start = s)
    by_period <- tapply(x$period_visits$n_resp, x$period_visits$period, sum)
    expect_equal(as.numeric(by_period), x$periods$n_resp, tolerance = 1e-12)
    issued <- tapply(x$period_visits$n_issued, x$period_visits$period, sum)
    expect_equal(as.numeric(issued), x$periods$n_in_sample, tolerance = 1e-12)
    # the visit column is what makes the early moe readable: an immediate
    # start's first period is entirely at visit 1
    first <- x$period_visits[x$period_visits$period == 1L, ]
    if (identical(s, "immediate")) {
      expect_equal(first$visit, 1L)
      expect_equal(first$n_issued, x$n_in_sample, tolerance = 1e-12)
    }
    expect_true(all(x$period_visits$n_resp <= x$period_visits$n_issued))
  }
})

test_that("a start-up holds for lives other than five, and for a mean", {
  for (k in c(2L, 3L, 8L)) {
    x <- n_panel(n_prop(p = 0.4, moe = 0.03), retention = rep(0.9, k - 1L),
                 resp_rate = 0.7, design = "rotating", start = "immediate")
    e <- x$n_entrants
    q <- x$waves$q
    expect_equal(x$periods$n_resp[[1L]], k * e * q[[1L]], tolerance = 1e-12)
    expect_equal(tail(x$periods$n_resp, 1L), x$n_target, tolerance = 1e-9)
  }
  m <- n_panel(n_mean(var = 100, moe = 2), retention = c(0.9, 0.85),
               resp_rate = 0.75, design = "rotating", start = "gradual")
  expect_equal(tail(m$periods$n_resp, 1L), m$n_target, tolerance = 1e-9)
  expect_true("rmoe" %in% names(m$periods))
})

test_that("start is validated and refused for a fixed panel", {
  expect_error(.lfs(start = "later"), "\"gradual\", \"immediate\" or NULL")
  expect_error(.lfs(start = c("gradual", "immediate")), "must be")
  expect_error(.lfs(start = NA_character_), "must be")
  expect_error(
    n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(0.9, 0.9),
            resp_rate = 0.8, start = "gradual"),
    "fixed panel recruits one cohort"
  )
})

test_that("prec_panel carries the start-up through the round trip", {
  plan <- .lfs(start = "immediate")
  back <- prec_panel(plan)
  expect_identical(back$start, "immediate")
  expect_equal(back$periods, plan$periods, tolerance = 1e-12)

  # start moves no stored quantity, so unlike design it may be overridden
  swapped <- prec_panel(plan, start = "gradual")
  expect_identical(swapped$start, "gradual")
  expect_equal(swapped$n_entrants, plan$n_entrants, tolerance = 1e-12)
  expect_equal(swapped$n_resp, plan$n_resp, tolerance = 1e-12)
  expect_lt(swapped$periods$n_resp[[1L]], plan$periods$n_resp[[1L]])

  # and a plan stored without one stays without one
  expect_null(prec_panel(.lfs())$periods)
  expect_identical(prec_panel(.lfs(), start = "gradual")$start, "gradual")
})

test_that("a start-up given directly to prec_panel reads the same", {
  plan <- .lfs(start = "gradual")
  direct <- prec_panel(plan$n_entrants, n_prop(p = 0.5, moe = 0.031),
                       retention = c(0.878, 0.963, 0.936, 0.956),
                       resp_rate = 0.728, design = "rotating", start = "gradual")
  expect_equal(direct$periods, plan$periods, tolerance = 1e-12)
})

test_that("start = NULL leaves no trace in params", {
  # `list(start = NULL)` carries the name, so a result planned without a
  # start-up would gain a parameter it never had
  a <- .lfs()
  expect_false("start" %in% names(a$params))
  expect_length(a$params, 5L)
  expect_named(a$params,
               c("resp_rate", "retention", "design", "target_wave", "assurance"))
  expect_true("start" %in% names(.lfs(start = "gradual")$params))
})

test_that("an immediate start opens at equilibrium when nothing is lost", {
  # R_1 - equilibrium = e * sum(q_1 - q_s), which is 0 exactly when q is flat, so
  # the claim is "at least as many", not "more"
  flat <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(1, 1, 1),
                  resp_rate = 0.8, design = "rotating", start = "immediate")
  expect_equal(flat$periods$n_resp[[1L]], tail(flat$periods$n_resp, 1L),
               tolerance = 1e-12)
  expect_equal(flat$periods$moe[[1L]], flat$moe, tolerance = 1e-12)
  # every period delivers the equilibrium count here, and the flag still
  # marks the period the wave composition matures at: with nothing lost the
  # counts coincide earlier than the mix does, and it is the mix that decides
  # whether a period is comparable with the rest of the series
  expect_equal(flat$periods$n_resp, rep(flat$n_target, nrow(flat$periods)),
               tolerance = 1e-9)
  expect_equal(which(flat$periods$equilibrium)[[1L]], nrow(flat$waves))
  expect_equal(flat$period_visits$visit[flat$period_visits$period == 1L], 1L)

  # resp_rate cancels from both sides, so it is retention alone that decides
  flat_full <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(1, 1, 1),
                       resp_rate = 1, design = "rotating", start = "immediate")
  expect_equal(flat_full$periods$n_resp[[1L]],
               tail(flat_full$periods$n_resp, 1L), tolerance = 1e-12)

  # and one wave of attrition is enough to separate them
  one <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(1, 0.99, 1),
                 resp_rate = 0.8, design = "rotating", start = "immediate")
  expect_gt(one$periods$n_resp[[1L]], tail(one$periods$n_resp, 1L))
  q <- one$waves$q
  expect_equal(one$periods$n_resp[[1L]] - tail(one$periods$n_resp, 1L),
               one$n_entrants * sum(q[[1L]] - q), tolerance = 1e-12)
})

## Rotations with a gap

.gap_target <- function() n_prop(p = 0.5, moe = 0.031)

test_that("a gap leaves the equilibrium of the same visits unchanged", {
  for (spec in c("2-2-2", "4-8-4", "2-1-3")) {
    v <- sum(as.double(design_rotation(spec)) > 0)
    ret <- seq(0.95, 0.85, length.out = v - 1L)
    plain <- n_panel(.gap_target(), retention = ret, resp_rate = 0.75,
                     design = "rotating", assurance = 0.9)
    gap <- n_panel(.gap_target(), retention = ret, resp_rate = 0.75,
                   design = "rotating", assurance = 0.9, rotation = spec)
    for (f in c("n_entrants", "n_in_sample", "n_cohorts", "n_resp",
                "n_assured", "se", "moe", "cv", "waves")) {
      expect_identical(gap[[f]], plain[[f]], label = paste(spec, f))
    }
    expect_equal(gap$n_held, length(gap$rotation) * gap$n_entrants)
    expect_gt(gap$n_held, gap$n_in_sample)
  }
})

test_that("a rotation is checked against the visits and the design", {
  expect_error(
    n_panel(.gap_target(), retention = c(0.9, 0.9), design = "rotating",
            rotation = "2-2-2"),
    "interviews a unit 4 times over its 6-period life"
  )
  expect_error(
    n_panel(.gap_target(), retention = c(0.9, 0.9, 0.9), design = "rotating",
            rotation = c(1, 1, 0.5, 0.5)),
    "subsamples a later wave"
  )
  expect_error(
    n_panel(.gap_target(), retention = c(0.9, 0.9), rotation = "3"),
    "a fixed panel follows one cohort"
  )
  expect_error(
    n_panel(.gap_target(), retention = c(0.9, 0.9, 0.9), design = "rotating",
            rotation = design_overlap("2-2-2")),
    "an overlap is what a rotation produces"
  )
  plain <- n_panel(.gap_target(), retention = c(0.9, 0.9), design = "rotating")
  expect_identical(as.double(plain$rotation), c(1, 1, 1))
  expect_identical(plain$n_held, plain$n_in_sample)
  expect_false("rotation" %in% names(plain$params))
})

test_that("the start-up of any rotation matches an enumeration of its cohorts", {
  specs <- c("3", "5", "2-2-2", "4-8-4", "2-1-3", "1-0-1", "1-1-0-1")
  for (spec in specs) for (start in c("gradual", "immediate")) {
    w <- as.double(design_rotation(spec))
    L <- length(w)
    v <- sum(w > 0)
    for (ret in list(rep(1, v - 1L), seq(0.9, 0.7, length.out = v - 1L))) {
      x <- n_panel(.gap_target(), retention = ret, resp_rate = 0.8,
                   design = "rotating", start = start, rotation = spec)
      lab <- paste(spec, start)
      o <- startup_oracle(w, start, L + 1L)
      tab <- stats::aggregate(id ~ period + visit, data = o, FUN = length)
      tab <- tab[order(tab$period, tab$visit), ]
      got <- x$period_visits
      expect_identical(got$period, tab$period, label = lab)
      expect_identical(got$visit, tab$visit, label = lab)
      expect_equal(got$n_issued, x$n_entrants * tab$id, tolerance = 1e-12,
                   label = lab)
      expect_equal(got$n_resp, x$n_entrants * tab$id * x$waves$q[tab$visit],
                   tolerance = 1e-12, label = lab)
      in_sample <- tabulate(o$period, L + 1L)
      expect_equal(x$periods$n_in_sample, x$n_entrants * in_sample,
                   tolerance = 1e-12, label = lab)
      expect_identical(x$periods$equilibrium, x$periods$period >= L,
                       label = lab)
      entering <- c(if (start == "immediate") L else 1, rep(1, L))
      expect_equal(x$periods$n_entrants, x$n_entrants * entering,
                   tolerance = 1e-12, label = lab)
    }
  }
})

test_that("an immediate 2-2-2 start-up is full from the first period", {
  x <- n_panel(.gap_target(), retention = c(0.9, 0.85, 0.9), resp_rate = 0.75,
               design = "rotating", start = "immediate", rotation = "2-2-2")
  expect_equal(x$periods$n_in_sample, rep(4 * x$n_entrants, 7L),
               tolerance = 1e-12)
  expect_identical(which(x$periods$equilibrium)[[1L]], 6L)
  grad <- n_panel(.gap_target(), retention = c(0.9, 0.85, 0.9),
                  resp_rate = 0.75, design = "rotating", start = "gradual",
                  rotation = "2-2-2")
  expect_equal(grad$periods$n_in_sample / grad$n_entrants,
               c(1, 2, 2, 2, 3, 4, 4), tolerance = 1e-12)

  # the issued overlap from the first period is the full-sample profile
  o <- startup_oracle(as.double(x$rotation), "immediate", 12L)
  want <- as.double(design_overlap("2-2-2"))
  for (lag in seq_along(want)) {
    a <- o$id[o$period == 1L]
    b <- o$id[o$period == 1L + lag]
    expect_equal(length(intersect(a, b)) / length(a), want[[lag]],
                 tolerance = 1e-12)
  }
})

test_that("a finite population has to supply the cohorts held in a gap", {
  # 4 cohorts in sample fit 1500 units, the 6 held do not
  tgt <- n_prop(p = 0.5, moe = 0.03, N = 1500)
  expect_error(
    n_panel(tgt, retention = c(0.9, 0.8, 0.9), resp_rate = 0.75,
            design = "rotating", rotation = "2-2-2"),
    "6 cohorts held at one period, 4 of them in sample"
  )
  ok <- n_panel(tgt, retention = c(0.9, 0.8, 0.9), resp_rate = 0.75,
                design = "rotating")
  expect_lt(ok$n_in_sample, 1500)
  expect_gt(6 * ok$n_entrants, 1500)
})

test_that("print names a rotation with a gap, and summary the units held", {
  gap <- n_panel(.gap_target(), retention = c(0.9, 0.85, 0.9),
                 resp_rate = 0.75, design = "rotating", rotation = "2-2-2")
  plain <- n_panel(.gap_target(), retention = c(0.9, 0.85, 0.9),
                   resp_rate = 0.75, design = "rotating")
  out <- capture.output(print(gap))
  expect_true(any(grepl("life: 6 periods (2 in, 2 out, 2 in)",
                        out, fixed = TRUE)))
  expect_false(any(grepl("life:", capture.output(print(plain)))))
  s <- capture.output(print(summary(gap)))
  expect_true(any(grepl("held: 2388 across 6 cohorts, 2 of them out of sample",
                        s, fixed = TRUE)))
  expect_false(any(grepl("held:", capture.output(print(summary(plain))))))
})

test_that("the start-up tables are named by what indexes their rows", {
  x <- n_panel(n_prop(p = 0.15, cv = 0.12), retention = c(0.9, 0.85, 0.9),
               resp_rate = 0.75, design = "rotating", start = "immediate",
               rotation = "2-2-2")
  expect_false(any(c("launch", "launch_waves") %in% names(x)))
  expect_named(x$period_visits, c("period", "visit", "n_issued", "n_resp"))
  expect_identical(unique(x$period_visits$visit[x$period_visits$period == 1L]),
                   1L)
  expect_identical(x$periods$period, seq_len(length(x$rotation) + 1L))
  out <- capture.output(print(summary(x)))
  expect_false(any(grepl("occasion|launch|Launch", out)))
  expect_true(max(nchar(out)) <= 80L)
})

## The target's population stays behind its respondent requirement

test_that("a panel keeps the population the target was sized against", {
  for (type in c("prop", "mean")) {
    target <- if (type == "prop") {
      n_prop(p = 0.5, moe = 0.05, N = 100, eligibility_rate = 0.5)
    } else {
      n_mean(var = 1, moe = 0.2, N = 100, eligibility_rate = 0.5)
    }
    expect_warning(
      plan <- n_panel(target, retention = 0.9),
      "sized against 50 units in scope"
    )
    expect_equal(plan$moe, target$moe)
    expect_warning(back <- prec_panel(plan), "keeps that population")
    expect_equal(back$moe, target$moe)
    # matching eligibility: no warning, same precision
    expect_no_warning(
      matched <- n_panel(target, retention = 0.9, eligibility_rate = 0.5)
    )
    expect_equal(matched$moe, target$moe)
    expect_equal(matched$n_issued, 2 * plan$n_issued)
    # rotating and the reverse direction
    expect_warning(
      rot <- n_panel(target, retention = 0.9, design = "rotating"),
      "keeps that population"
    )
    expect_equal(rot$moe, target$moe)
    expect_equal(
      suppressWarnings(prec_panel(rot$n_entrants, target, retention = 0.9,
                                  design = "rotating"))$moe,
      target$moe
    )
  }
  # the panel's eligibility alone does not shrink the target's population
  t1 <- prec_prop(p = 0.5, n = 40, N = 100)
  expect_warning(p1 <- n_panel(t1, retention = 0.9, eligibility_rate = 0.5),
                 "sized against 100 units in scope")
  expect_equal(p1$moe, t1$moe)
  # infinite N has no population to disagree about
  expect_silent(n_panel(prec_prop(p = 0.5, n = 40), retention = 0.9,
                        eligibility_rate = 0.5))
})

## Printed precision describes the whole units printed

test_that("print() shows the precision of the whole-unit recruitment", {
  plan <- n_panel(prec_mean(var = 1, n = 2.1), retention = 1,
                  design = "rotating", rotation = "2")
  expect_equal(plan$se, sqrt(1 / 2.1))
  out <- capture.output(print(plan))
  expect_true(any(grepl("^entrants: 2 per period -> 4 responding", out)))
  expect_true(any(grepl("^mean: se = 0.5, moe = 0.98", out)))
  expect_true(any(grepl("^ 1 +2 +0.7071", out)))
  # a whole recruitment prints the plan itself
  whole <- prec_panel(3, prec_mean(var = 1, n = 6), retention = 1,
                      design = "rotating", rotation = "2")
  expect_true(any(grepl(sprintf("se = %.4g", whole$se),
                        capture.output(print(whole)))))
})

test_that("summary() reports the complete histories of a rotating plan", {
  plan <- n_panel(prec_mean(var = 1, n = 1000), retention = rep(0.9, 3),
                  resp_rate = 0.8, design = "rotating", rotation = "4")
  expect_equal(plan$n_entrants * plan$waves$q[4], 211.9802, tolerance = 1e-6)
  out <- capture.output(print(summary(plan)))
  expect_true(any(grepl(
    "^complete histories: 212 per period, the cohort at wave 4 with all 4 responses$",
    out
  )))
  fixed <- n_panel(prec_mean(var = 1, n = 100), retention = 0.9)
  expect_false(any(grepl("complete histories",
                         capture.output(print(summary(fixed))))))
})

test_that("a whole design beyond the frame prints the continuous one throughout", {
  expect_warning(
    plan <- n_panel(prec_mean(var = 100, n = 100, N = 101),
                    retention = c(1, 1), design = "rotating",
                    start = "immediate"),
    "102 units in sample against a population of 101"
  )
  out <- capture.output(print(summary(plan)))
  expect_true(any(grepl("^entrants: 33.33 per period -> 100 responding", out)))
  expect_true(any(grepl(
    "^whole units: 34 per period hold 102, above N = 101, so these are continuous$",
    out
  )))
  expect_true(any(grepl("^in sample: 100 across 3 live cohorts$", out)))
  expect_true(any(grepl("^ 1 +1 +33 +1.418", out)))
  expect_true(any(grepl("^ 3 +33 +100 +100 +0.0995", out)))
  expect_false(any(grepl("102 responding|^ 1 +34 ", out)))
  # a plan whose whole units fit carries no status line
  fits <- n_panel(prec_mean(var = 100, n = 90, N = 101), retention = c(1, 1),
                  design = "rotating")
  expect_false(any(grepl("^whole units:", capture.output(print(fits)))))
})
