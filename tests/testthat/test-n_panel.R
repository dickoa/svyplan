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

test_that("a rotating panel takes 322 entrants per occasion and holds 1608", {
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
  # A single cohort's issue and the entrants per occasion differ by a factor
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

test_that("a rotating occasion pools its cohorts and a fixed wave does not", {
  rot <- n_panel(
    lfs_target, retention = lfs_ret, resp_rate = lfs_rr, design = "rotating"
  )
  expect_equal(rot$n_resp, sum(rot$waves$n_resp), tolerance = 1e-12)
  # The pooled occasion is the target; one cohort at one stage is a fifth of
  # it, which is the reading the plan warns against.
  expect_lt(max(rot$waves$n_resp), rot$n_resp / 4)
})

## P3. A panel of one loss is the response inflation the package already has

test_that("full retention reduces to the resp_rate inflation of n_prop", {
  # q = (r, r), so the requirement at wave 2 is n_target / r, which is what
  # resp_rate does on a single occasion.
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
    netted <- n_panel(cases[[nm]][[2L]], retention = ret, resp_rate = 0.8)
    expect_equal(netted$n_issued, full$n_issued, tolerance = 1e-12,
                 info = nm)
    expect_equal(netted$moe, full$moe, tolerance = 1e-12, info = nm)
  }
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

test_that("both directions of the single-occasion families are accepted", {
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

test_that("the assured occasion clears its target and the count below does not", {
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

test_that("a rotating occasion is Poisson-binomial, not one binomial", {
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
  # Retention of 1 makes every live cohort identical, so the occasion is
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
  # time, so the frame it outgrows is the one it fills at one occasion.
  expect_error(
    n_panel(
      prec_mean(var = 100, n = 1000, N = 1200),
      retention = c(0.9, 0.8), resp_rate = 0.8, design = "rotating"
    ),
    "cohorts alive at one occasion"
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
  # entrants per occasion, which three equal cohorts of 34 cannot field
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
  expect_match(out[2L], "^issue 1816 to hold 1000 responding at wave 5$")
  expect_true(any(grepl("61% of the life's loss at wave 1", out, fixed = TRUE)))
  expect_true(any(grepl("assured \\(0.95\\)", out)))
  expect_length(grep("^ *[1-5] ", out), 5L)
})

test_that("print names the entrants and the launch separately", {
  rot <- n_panel(lfs_target, retention = lfs_ret, resp_rate = lfs_rr,
                 design = "rotating")
  out <- capture.output(print(rot))
  expect_match(out[2L], "^322 entrants per occasion")
  expect_match(out[3L], "^1610 in sample across 5 live cohorts$")
  expect_true(any(grepl("cohorts alive at one occasion", out)))
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

## N-launch. What a rotating design delivers before it reaches its steady state

.lfs <- function(...) {
  n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(0.878, 0.963, 0.936, 0.956),
          resp_rate = 0.728, design = "rotating", ...)
}

test_that("start = NULL is inert", {
  # the default must add nothing: no field, no print line, no changed number
  a <- .lfs()
  b <- .lfs(start = "immediate")
  expect_null(a$start)
  expect_null(a$launch)
  expect_null(a$launch_waves)
  expect_null(a$params$start)
  expect_equal(a$n_entrants, b$n_entrants, tolerance = 1e-12)
  expect_equal(a$n_resp, b$n_resp, tolerance = 1e-12)
  expect_equal(a$moe, b$moe, tolerance = 1e-12)
  expect_equal(a$waves, b$waves)
  expect_false(any(grepl("launch", capture.output(print(a)))))
  expect_true(any(grepl("launch \\(immediate\\)", capture.output(print(b)))))
})

test_that("the immediate launch follows the closed form", {
  x <- .lfs(start = "immediate")
  e <- x$n_entrants
  q <- x$waves$q
  k <- length(q)
  expected <- vapply(seq_len(k), function(t) {
    e * ((k - t + 1) * q[[t]] + sum(q[seq_len(t - 1L)]))
  }, numeric(1L))
  expect_equal(x$launch$n_resp[seq_len(k)], expected, tolerance = 1e-12)
  # period 1 is one binomial: every unit is at wave 1
  expect_equal(x$launch$n_resp[[1L]], k * e * q[[1L]], tolerance = 1e-12)
})

test_that("the gradual launch is the cumulative recruitment", {
  x <- .lfs(start = "gradual")
  e <- x$n_entrants
  q <- x$waves$q
  expect_equal(
    x$launch$n_resp[seq_along(q)],
    vapply(seq_along(q), function(t) e * sum(q[seq_len(t)]), numeric(1L)),
    tolerance = 1e-12
  )
  expect_equal(x$launch$n_entrants, rep(e, length(q) + 1L), tolerance = 1e-12)
})

test_that("the launch is stored continuous and lands on n_target", {
  # the stored table must not carry whole units: they belong to print()
  for (s in c("gradual", "immediate")) {
    x <- .lfs(start = s)
    last <- nrow(x$launch)
    expect_equal(x$launch$n_resp[[last]], x$n_target, tolerance = 1e-9)
    expect_equal(x$launch$n_resp[[last]], x$n_entrants * sum(x$waves$q),
                 tolerance = 1e-9)
    expect_equal(x$launch$moe[[last]], x$moe, tolerance = 1e-12)
    expect_false(isTRUE(all.equal(x$launch$n_entrants[[last]],
                                  ceiling(x$n_entrants))))
  }
})

test_that("the two launches bracket the steady state in opposite directions", {
  imm <- .lfs(start = "immediate")$launch
  grad <- .lfs(start = "gradual")$launch
  steady <- tail(imm$n_resp, 1L)
  # an immediate start opens above its own steady state, every unit being at
  # wave 1, and settles down onto it; a gradual one opens below and climbs
  expect_gt(imm$n_resp[[1L]], steady)
  expect_lt(grad$n_resp[[1L]], steady)
  expect_false(is.unsorted(rev(imm$n_resp)))
  expect_false(is.unsorted(grad$n_resp))
  expect_lt(imm$moe[[1L]], tail(imm$moe, 1L))
  # and an immediate start holds the whole design from its first occasion
  expect_equal(imm$n_in_sample, rep(imm$n_in_sample[[1L]], nrow(imm)),
               tolerance = 1e-12)
  expect_lt(grad$n_in_sample[[1L]], tail(grad$n_in_sample, 1L))
})

test_that("the steady state is flagged where the two launches meet", {
  imm <- .lfs(start = "immediate")
  grad <- .lfs(start = "gradual")
  k <- nrow(imm$waves)
  expect_equal(which(imm$launch$steady_state)[[1L]], k)
  expect_equal(which(grad$launch$steady_state)[[1L]], k)
  expect_equal(imm$launch$n_resp[imm$launch$steady_state],
               grad$launch$n_resp[grad$launch$steady_state], tolerance = 1e-12)
})

test_that("launch_waves decomposes each occasion", {
  for (s in c("gradual", "immediate")) {
    x <- .lfs(start = s)
    by_period <- tapply(x$launch_waves$n_resp, x$launch_waves$period, sum)
    expect_equal(as.numeric(by_period), x$launch$n_resp, tolerance = 1e-12)
    issued <- tapply(x$launch_waves$n_issued, x$launch_waves$period, sum)
    expect_equal(as.numeric(issued), x$launch$n_in_sample, tolerance = 1e-12)
    # the wave column is what makes the early moe readable: an immediate
    # start's first occasion is entirely wave 1
    first <- x$launch_waves[x$launch_waves$period == 1L, ]
    if (identical(s, "immediate")) {
      expect_equal(first$wave, 1L)
      expect_equal(first$n_issued, x$n_in_sample, tolerance = 1e-12)
    }
    expect_true(all(x$launch_waves$n_resp <= x$launch_waves$n_issued))
  }
})

test_that("a launch holds for lives other than five, and for a mean", {
  for (k in c(2L, 3L, 8L)) {
    x <- n_panel(n_prop(p = 0.4, moe = 0.03), retention = rep(0.9, k - 1L),
                 resp_rate = 0.7, design = "rotating", start = "immediate")
    e <- x$n_entrants
    q <- x$waves$q
    expect_equal(x$launch$n_resp[[1L]], k * e * q[[1L]], tolerance = 1e-12)
    expect_equal(tail(x$launch$n_resp, 1L), x$n_target, tolerance = 1e-9)
  }
  m <- n_panel(n_mean(var = 100, moe = 2), retention = c(0.9, 0.85),
               resp_rate = 0.75, design = "rotating", start = "gradual")
  expect_equal(tail(m$launch$n_resp, 1L), m$n_target, tolerance = 1e-9)
  expect_true("rmoe" %in% names(m$launch))
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

test_that("prec_panel carries the launch through the round trip", {
  plan <- .lfs(start = "immediate")
  back <- prec_panel(plan)
  expect_identical(back$start, "immediate")
  expect_equal(back$launch, plan$launch, tolerance = 1e-12)

  # start moves no stored quantity, so unlike design it may be overridden
  swapped <- prec_panel(plan, start = "gradual")
  expect_identical(swapped$start, "gradual")
  expect_equal(swapped$n_entrants, plan$n_entrants, tolerance = 1e-12)
  expect_equal(swapped$n_resp, plan$n_resp, tolerance = 1e-12)
  expect_lt(swapped$launch$n_resp[[1L]], plan$launch$n_resp[[1L]])

  # and a plan stored without one stays without one
  expect_null(prec_panel(.lfs())$launch)
  expect_identical(prec_panel(.lfs(), start = "gradual")$start, "gradual")
})

test_that("a launch given directly to prec_panel reads the same", {
  plan <- .lfs(start = "gradual")
  direct <- prec_panel(plan$n_entrants, n_prop(p = 0.5, moe = 0.031),
                       retention = c(0.878, 0.963, 0.936, 0.956),
                       resp_rate = 0.728, design = "rotating", start = "gradual")
  expect_equal(direct$launch, plan$launch, tolerance = 1e-12)
})

test_that("start = NULL leaves no trace in params", {
  # `list(start = NULL)` carries the name, so a result planned without a
  # launch would gain a parameter it never had
  a <- .lfs()
  expect_false("start" %in% names(a$params))
  expect_length(a$params, 5L)
  expect_named(a$params,
               c("resp_rate", "retention", "design", "target_wave", "assurance"))
  expect_true("start" %in% names(.lfs(start = "gradual")$params))
})

test_that("an immediate launch opens at the steady state when nothing is lost", {
  # R_1 - steady = e * sum(q_1 - q_s), which is 0 exactly when q is flat, so
  # the claim is "at least as many", not "more"
  flat <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(1, 1, 1),
                  resp_rate = 0.8, design = "rotating", start = "immediate")
  expect_equal(flat$launch$n_resp[[1L]], tail(flat$launch$n_resp, 1L),
               tolerance = 1e-12)
  expect_equal(flat$launch$moe[[1L]], flat$moe, tolerance = 1e-12)
  # every occasion delivers the steady-state count here, and the flag still
  # marks the occasion the wave composition matures at: with nothing lost the
  # counts coincide earlier than the mix does, and it is the mix that decides
  # whether an occasion is comparable with the rest of the series
  expect_equal(flat$launch$n_resp, rep(flat$n_target, nrow(flat$launch)),
               tolerance = 1e-9)
  expect_equal(which(flat$launch$steady_state)[[1L]], nrow(flat$waves))
  expect_equal(flat$launch_waves$wave[flat$launch_waves$period == 1L], 1L)

  # resp_rate cancels from both sides, so it is retention alone that decides
  flat_full <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(1, 1, 1),
                       resp_rate = 1, design = "rotating", start = "immediate")
  expect_equal(flat_full$launch$n_resp[[1L]],
               tail(flat_full$launch$n_resp, 1L), tolerance = 1e-12)

  # and one wave of attrition is enough to separate them
  one <- n_panel(n_prop(p = 0.5, moe = 0.031), retention = c(1, 0.99, 1),
                 resp_rate = 0.8, design = "rotating", start = "immediate")
  expect_gt(one$launch$n_resp[[1L]], tail(one$launch$n_resp, 1L))
  q <- one$waves$q
  expect_equal(one$launch$n_resp[[1L]] - tail(one$launch$n_resp, 1L),
               one$n_entrants * sum(q[[1L]] - q), tolerance = 1e-12)
})
