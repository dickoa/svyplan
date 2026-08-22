## The 'psu' argument: certainty-aware allocation from a PSU register

test_that("a PSU register refines the allocation without changing its class", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  expect_s3_class(fit, "svyplan_n")
  expect_identical(fit$type, "alloc")
  expect_identical(fit$method, "bethel")
  expect_true(all(fit$constraints$.pass))
  # The plan meets its targets exactly, not merely within tolerance: the
  # settle step exists to make that true.
  expect_equal(max(fit$constraints$.achieved / fit$constraints$.target), 1,
               tolerance = 1e-6)
})

test_that("the register in params is the stable mark of a certainty fit", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  # Consumers key on the register's presence, never on detail column names.
  expect_false(is.null(fit$params$psu))
  expect_setequal(fit$params$psu$psu_id, z$psu$psu_id)

  plain <- n_alloc(
    z$frame[, c("stratum", "N")],
    measures = z$measures[, setdiff(names(z$measures), "icc_psu")],
    targets = z$targets
  )
  expect_null(plain$params$psu)
})

test_that("the classification is reported per stratum and per PSU", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  expect_true(all(c("n_psu_certain", "n_psu_rest", "threshold") %in%
                    names(fit$detail)))
  # Counts are PSUs available in each part, so together they are the register.
  expect_equal(
    fit$detail$n_psu_certain + fit$detail$n_psu_rest,
    as.numeric(table(factor(z$psu$stratum, levels = z$frame$stratum))),
    ignore_attr = TRUE
  )

  expect_equal(nrow(fit$psu), nrow(z$psu))
  expect_true(all(c("psu_id", "stratum", "N", "certainty",
                    ".certainty_source", ".threshold", ".distance") %in%
                    names(fit$psu)))
  expect_type(fit$psu$certainty, "logical")
  expect_true(all(is.na(fit$psu$.certainty_source[!fit$psu$certainty])))
  expect_true(all(!is.na(fit$psu$.certainty_source[fit$psu$certainty])))
})

test_that("the threshold is the one the returned allocation implies", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  # threshold = take / sampling fraction, read off the answer itself. R2BEAT
  # returns a threshold from an earlier iterate than its own allocation.
  expect_equal(fit$detail$threshold,
               z$frame$n_per_psu / (fit$detail$n / fit$detail$N))
  expect_equal(fit$psu$.distance,
               fit$psu$N / fit$psu$.threshold - 1)
})

test_that("every PSU above the threshold is certainty", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  # A PSU at or above the threshold has inclusion probability at least one
  # and cannot be sampled less often, so leaving one out is unexecutable.
  above <- fit$psu$N >= fit$psu$.threshold
  expect_true(all(fit$psu$certainty[above]))
  expect_true(all(fit$psu$.certainty_source[above] == "threshold"))
})

test_that("n_take is the take the operational design fields in each PSU", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  expect_true("n_take" %in% names(fit$psu))
  cert <- fit$psu$certainty
  take_of <- z$frame$n_per_psu[match(fit$psu$stratum, z$frame$stratum)]

  # The remainder is fielded at the stated per-PSU take.
  expect_equal(fit$psu$n_take[!cert], take_of[!cert], ignore_attr = TRUE)

  # Every certainty PSU here reached the threshold on its own, so its whole
  # take at the stratum rate is at least the stated take, and never more
  # than the PSU holds.
  expect_true(all(fit$psu$.certainty_source[cert] == "threshold"))
  expect_true(all(fit$psu$n_take[cert] >= take_of[cert]))
  expect_true(all(fit$psu$n_take[cert] <= fit$psu$N[cert]))

  # The takes are the operational design's own numbers, not a re-derivation:
  # they sum to n_certain_int exactly, and the whole-unit identity holds.
  agg <- tapply(fit$psu$n_take * cert, fit$psu$stratum, sum)
  expect_equal(as.numeric(agg[fit$detail$stratum]), fit$detail$n_certain_int,
               ignore_attr = TRUE)
  expect_equal(
    fit$detail$n_int,
    fit$detail$n_certain_int + fit$detail$n_psu_draw * z$frame$n_per_psu
  )
})

test_that("the convergence verdict is reported and never silent", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  cc <- fit$optimization$certainty

  expect_true(cc$verdict %in% c("converged", "cycle", "limit_reached"))
  expect_type(cc$fixed_point, "logical")
  expect_gte(cc$orbit, 1L)
  expect_gte(cc$settle_iterations, 1L)

  # A cycling register still returns a plan, and says it is not a fixed point.
  zc <- .psu_fixture(seed = 3L)
  cyc <- n_alloc(zc$frame, measures = zc$measures, targets = zc$targets,
                 psu = zc$psu)
  expect_identical(cyc$optimization$certainty$verdict, "cycle")
  expect_false(cyc$optimization$certainty$fixed_point)
  expect_true(all(cyc$constraints$.pass))
})

test_that("the certainty flag adds units and never removes them", {
  z <- .psu_fixture()
  bare <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = z$psu)

  # Flag PSUs that sit below the threshold, so the flag has work to do. The
  # largest PSUs are certainty on their own and would test nothing.
  flagged <- z$psu
  ratio <- bare$psu$N / bare$psu$.threshold
  flagged$certainty <- ratio > 0.5 & ratio < 0.9
  expect_gt(sum(flagged$certainty), 0)

  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = flagged)

  expect_gte(sum(fit$psu$certainty), sum(bare$psu$certainty))
  expect_true(all(flagged$certainty <= fit$psu$certainty))
  expect_true("supplied" %in% fit$psu$.certainty_source)

  # A flag cannot take a PSU out of the certainty part: above the threshold
  # its inclusion probability is one whatever the caller asks.
  denied <- z$psu
  denied$certainty <- FALSE
  kept <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = denied)
  expect_true(any(kept$psu$certainty))
})

test_that("more certainty means less clustering and a smaller sample", {
  z <- .psu_fixture()
  bare <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = z$psu)
  flagged <- z$psu
  flagged$certainty <- flagged$N >= stats::quantile(flagged$N, 0.90)
  more <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = flagged)

  expect_gt(sum(more$psu$certainty), sum(bare$psu$certainty))
  expect_lt(more$n, bare$n)
})

test_that("a register and a bare PSU count are alternatives", {
  z <- .psu_fixture()
  with_count <- z$frame
  with_count$N_psu <- c(300, 200, 140, 90)
  expect_error(
    n_alloc(with_count, measures = z$measures, targets = z$targets,
            psu = z$psu),
    "alternatives"
  )
})

test_that("the register is refused when it cannot describe the frame", {
  z <- .psu_fixture()

  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            psu = z$psu[, c("psu_id", "stratum")]),
    "must contain"
  )
  bad <- z$psu
  bad$N[1] <- 0
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = bad),
    "positive and finite"
  )
  stray <- z$psu
  stray$stratum[1] <- "Z"
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = stray),
    "not in 'frame'"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            psu = z$psu[z$psu$stratum != "D", ]),
    "no PSU"
  )
  no_take <- z$frame
  no_take$n_per_psu <- NULL
  expect_error(
    n_alloc(no_take, measures = z$measures, targets = z$targets, psu = z$psu),
    "n_per_psu"
  )
  no_icc <- z$measures
  no_icc$icc_psu <- NULL
  expect_error(
    n_alloc(z$frame, measures = no_icc, targets = z$targets, psu = z$psu),
    "icc_psu"
  )
  wrong_total <- z$psu
  wrong_total$N[wrong_total$stratum == "A"] <-
    wrong_total$N[wrong_total$stratum == "A"] * 2
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            psu = wrong_total),
    "sum of 'psu\\$N'"
  )
})

test_that("the register refuses designs this model does not describe", {
  z <- .psu_fixture()

  three <- z$frame
  three$n_per_ssu <- 4
  expect_error(
    n_alloc(three, measures = z$measures, targets = z$targets, psu = z$psu),
    "two-stage"
  )
  census <- z$frame
  census$take_all <- c(TRUE, FALSE, FALSE, FALSE)
  expect_error(
    n_alloc(census, measures = z$measures, targets = z$targets, psu = z$psu),
    "ultimate-unit census"
  )
  # Outside a joint allocation there are no targets for the loop to solve.
  expect_error(
    n_alloc(z$frame, n = 600, psu = z$psu),
    "requires a joint constrained allocation"
  )
})

test_that("a supplied design effect multiplies the clustering it does not model", {
  z <- .psu_fixture()
  bare <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = z$psu)
  doubled <- z$measures
  doubled$deff <- 2
  more <- n_alloc(z$frame, measures = doubled, targets = z$targets,
                  psu = z$psu)

  # It compounds with the certainty design effect rather than replacing it,
  # so the sample rises.
  expect_gt(more$n, bare$n)
})

## prec_alloc(psu =): assessing a supplied allocation

test_that("a supplied allocation is assessed against its own classification", {
  z <- .psu_fixture()
  n_h <- c(300, 200, 130, 80)
  ev <- prec_alloc(z$frame, n = n_h, measures = z$measures,
                   targets = z$targets, psu = z$psu)

  expect_s3_class(ev, "svyplan_prec")
  expect_identical(ev$type, "alloc")
  expect_equal(nrow(ev$psu), nrow(z$psu))

  # No loop is needed here: the allocation is given, so the threshold it
  # implies is given with it and the split is self-consistent by construction.
  expect_equal(ev$params$certainty$threshold,
               z$frame$n_per_psu / (n_h / z$frame$N))
  above <- ev$psu$N >= ev$psu$.threshold
  expect_identical(ev$psu$certainty, above)
})

test_that("the round trip from n_alloc closes exactly", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  ev <- prec_alloc(fit)

  expect_equal(ev$detail$.achieved, fit$constraints$.achieved,
               tolerance = 1e-10)
  # The held classification travels rather than being derived again, so a
  # PSU the loop absorbed is not silently dropped on the way back.
  expect_identical(ev$psu$certainty, fit$psu$certainty)
})

test_that("the assessed allocation reports its own takes", {
  z <- .psu_fixture()
  n_h <- c(300, 200, 130, 80)
  ev <- prec_alloc(z$frame, n = n_h, measures = z$measures,
                   targets = z$targets, psu = z$psu)

  expect_true("n_take" %in% names(ev$psu))
  cert <- ev$psu$certainty
  take_of <- z$frame$n_per_psu[match(ev$psu$stratum, z$frame$stratum)]
  f_of <- (n_h / z$frame$N)[match(ev$psu$stratum, z$frame$stratum)]

  expect_equal(ev$psu$n_take[!cert], take_of[!cert], ignore_attr = TRUE)
  # The certainty take is read off the supplied allocation, whole and capped
  # at the PSU's size.
  expect_equal(
    ev$psu$n_take[cert],
    pmin(ceiling(f_of[cert] * ev$psu$N[cert]), ev$psu$N[cert]),
    ignore_attr = TRUE
  )
})

test_that("the round trip reproduces the plan's takes exactly", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  ev <- prec_alloc(fit)

  # prec_alloc(fit) assesses the settled continuous allocation under the held
  # classification, which is the pair the operational takes were built from.
  expect_equal(ev$psu$n_take, fit$psu$n_take)
})

test_that("a larger allocation lowers the threshold and buys certainty", {
  z <- .psu_fixture()
  small <- prec_alloc(z$frame, n = c(200, 140, 90, 55), measures = z$measures,
                      targets = z$targets, psu = z$psu)
  large <- prec_alloc(z$frame, n = c(800, 560, 360, 220), measures = z$measures,
                      targets = z$targets, psu = z$psu)

  expect_true(all(large$params$certainty$threshold <
                    small$params$certainty$threshold))
  expect_gt(sum(large$psu$certainty), sum(small$psu$certainty))
})

test_that("prec_alloc refuses the register outside a joint allocation", {
  z <- .psu_fixture()
  expect_error(
    prec_alloc(z$frame, n = c(300, 200, 130, 80), psu = z$psu),
    "requires a joint constrained allocation"
  )
  expect_error(
    prec_alloc(z$frame, n = c(300, 200, 130, 80), measures = z$measures,
               targets = z$targets, psu = z$psu[, c("psu_id", "stratum")]),
    "must contain"
  )
  expect_error(
    prec_alloc(z$frame, n = 300, measures = z$measures,
               targets = z$targets, psu = z$psu),
    "length nrow\\(frame\\)"
  )
})

## Display: print carries the answer, summary carries what it rests on

test_that("print reports the certainty split and stays terse", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  out <- capture.output(print(fit))

  expect_true(any(grepl("^Joint constrained allocation \\(Bethel, PSU register\\)$",
                        out)))
  expect_true(any(grepl(
    "^PSUs: \\d+ certainty, \\d+ to draw from \\d+ in the remainder$", out
  )))
  expect_true(any(grepl("^# summary\\(\\) for the certainty split", out)))
  expect_lte(max(nchar(out)), 80L)
  expect_lte(length(out), 8L)

  # The counts on the block are the plan's own, and together they are the
  # register: the split is the design decision, not a derived figure.
  n_certain <- sum(fit$psu$certainty)
  expect_true(any(grepl(
    sprintf("^PSUs: %d certainty, %d to draw from %d in the remainder$",
            n_certain, fit$operational$n_psu_draw, nrow(z$psu) - n_certain),
    out
  )))
})

test_that("a plan with no fixed point says so on the block", {
  z <- .psu_fixture(seed = 2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  out <- capture.output(print(fit))

  expect_false(fit$optimization$certainty$fixed_point)
  expect_true(any(grepl("^no self-consistent classification", out)))
  expect_lte(max(nchar(out)), 80L)

  # A converged plan does not carry the line.
  ok <- n_alloc(.psu_fixture(seed = 4L)$frame,
                measures = z$measures, targets = z$targets,
                psu = .psu_fixture(seed = 4L)$psu)
  expect_true(ok$optimization$certainty$fixed_point)
  expect_false(any(grepl("no self-consistent", capture.output(print(ok)))))
})

test_that("an allocation with no register prints exactly as before", {
  z <- .bethel_fixture()
  out <- capture.output(print(
    n_alloc(z$frame, measures = z$measures, targets = z$targets)
  ))
  expect_true(any(grepl("^Joint constrained allocation \\(Bethel\\)$", out)))
  expect_false(any(grepl("PSU register|PSUs:|certainty", out)))
})

test_that("summary carries the split by stratum and how the loop resolved", {
  z <- .psu_fixture(seed = 2L)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  s <- summary(fit)

  expect_equal(nrow(s$certainty), nrow(z$frame))
  expect_identical(names(s$certainty),
                   c("stratum", "certainty", "remainder", "threshold"))
  expect_equal(s$certainty$certainty, fit$detail$n_psu_certain)
  expect_equal(s$certainty$certainty + s$certainty$remainder,
               as.numeric(table(factor(z$psu$stratum,
                                       levels = z$frame$stratum))),
               ignore_attr = TRUE)

  out <- capture.output(print(s))
  expect_true(any(grepl("^Certainty split by stratum$", out)))
  expect_true(any(grepl("certainty verdict: cycle", out)))
  expect_true(any(grepl("returned plan is a fixed point: no", out)))

  # The solver saw a one-stage problem, but the plan is not one.
  expect_true(any(grepl("variance model: two_part_certainty_wald", out)))
})


## The fielded design is a whole-PSU design

test_that("the whole-unit design is one a field team can draw", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  d <- fit$detail

  # An element-level rounding of this allocation is not a design: the
  # remainder has to be a whole number of PSUs at the stated take.
  expect_equal(d$n_int, d$n_certain_int + d$n_psu_draw * z$frame$n_per_psu)
  expect_equal(fit$operational$n, sum(d$n_int))
  expect_true(all(d$n_psu_draw <= d$n_psu_rest))
  expect_true(all(d$n_psu_draw == round(d$n_psu_draw)))

  # Rounding both parts up can only add sample, so the fielded design is at
  # least the continuous optimum and still meets every target.
  expect_gte(fit$operational$n, fit$n)
  expect_true(fit$operational$all_pass)
  expect_true(all(fit$operational$constraints$.pass))
})

test_that("every figure on the operational block comes from one design", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  # Cost, count and constraints all describe the whole-PSU design. The
  # solver's element-level integerization describes a total this plan does
  # not field, and mixing the two is what makes a printed block fail to add up.
  expect_equal(fit$operational$cost,
               sum(fit$detail$n_int * (fit$params$cost_h %||% 1)))
  expect_equal(fit$operational$n_certain, sum(fit$detail$n_certain_int))
  expect_equal(fit$operational$n_psu_draw, sum(fit$detail$n_psu_draw))

  out <- capture.output(print(fit))
  expect_true(any(grepl(
    sprintf("field design: n = %d, cost = %d", fit$operational$n,
            round(fit$operational$cost)),
    out
  )))
})

## Stage costs, and pricing the take as a grid

test_that("stage costs price the design and move the allocation", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets, psu = z$psu)
  d <- fit$detail

  # The marginal cost of an ultimate unit is the remainder's, a unit added
  # there bringing a share of a PSU visit with it.
  expect_equal(fit$params$cost_h, rep(400 / 12 + 20, nrow(fr)))

  # The design is priced as it is fielded: a visit to every PSU entered,
  # certainty or drawn, plus every interview.
  expect_equal(
    fit$operational$cost,
    sum((d$n_psu_certain + d$n_psu_draw) * 400 + d$n_int * 20)
  )
  # Both readings on the printed block are on the same basis, so they differ
  # by integerization rather than by what they count.
  expect_gt(fit$operational$cost, fit$params$achieved$cost)

  # A stratum that is dear to enter draws less of the sample.
  dear <- fr
  dear$cost_psu <- c(400, 400, 400, 4000)
  moved <- n_alloc(dear, measures = z$measures, targets = z$targets,
                   psu = z$psu)
  expect_lt(moved$detail$n[4], fit$detail$n[4])
})

test_that("stage costs are refused when they are not a usable pair", {
  z <- .psu_fixture()
  half <- z$frame
  half$cost_psu <- 400
  expect_error(
    n_alloc(half, measures = z$measures, targets = z$targets, psu = z$psu),
    "come as a pair"
  )
  both <- z$frame
  both$cost_psu <- 400
  both$cost_ssu <- 20
  expect_error(
    n_alloc(both, measures = z$measures, targets = z$targets, psu = z$psu,
            unit_cost = rep(1, 4)),
    "do not supply 'unit_cost'"
  )
  bad <- both
  bad$cost_ssu <- 0
  expect_error(
    n_alloc(bad, measures = z$measures, targets = z$targets, psu = z$psu),
    "positive and finite"
  )
})

test_that("predict() sweeps the take on a certainty fit", {
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  fit <- n_alloc(fr, measures = z$measures, targets = z$targets, psu = z$psu)

  # The solver saw a one-stage problem, so the sweep cannot find the take
  # through `stages`; it is reachable because the fit carries a register.
  expect_identical(fit$params$stages, 1L)
  g <- predict(fit, data.frame(n_per_psu = c(6, 12, 24, 40)))

  expect_equal(nrow(g), 4L)
  expect_true(all(g$.feasible))
  expect_true(all(c("n_psu_certain", "n_psu_draw", "n_int", "cost_int") %in%
                    names(g)))

  # A larger take raises the threshold, so fewer PSUs are certainty and the
  # clustering loss costs more sample.
  expect_true(all(diff(g$n_psu_certain) <= 0))
  expect_true(all(diff(g$n) > 0))

  # The row at the fit's own take reproduces the fit.
  at_fit <- predict(fit, data.frame(n_per_psu = 12))
  expect_equal(at_fit$n, fit$n, tolerance = 1e-8)
  expect_equal(at_fit$n_int, fit$operational$n)

  # The cost curve over the take is U-shaped, which is the reason to sweep it.
  fine <- predict(fit, data.frame(n_per_psu = 10:22))
  best <- which.min(fine$cost_int)
  expect_gt(best, 1L)
  expect_lt(best, nrow(fine))
})

## The register describes the frame it refines

test_that("the register total must equal the frame's in every stratum", {
  z <- .psu_fixture()
  # The variance model splits the stratum into a certainty and a remainder
  # size and scales by the frame's N, so a register that does not add up
  # makes the design effect internally inconsistent.
  over <- z$psu
  over$N[1] <- over$N[1] + 1
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = over),
    "must equal 'frame\\$N' in stratum"
  )
  short <- z$psu[-1, ]
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = short),
    "must equal 'frame\\$N'"
  )
})

test_that("the register total is compared at the scale of its terms", {
  # A derived sum against a supplied total: macOS arm64 accumulates sum()
  # without the extended precision x86-64 Linux uses, so an exact comparison
  # is a test of the platform's accumulator. The tolerance is relative, so it
  # holds at any scale and still catches a real mismatch.
  z <- .psu_fixture()
  small <- z$frame
  small$N <- small$N * 1e-9
  psu_small <- z$psu
  psu_small$N <- psu_small$N * 1e-9
  expect_silent(svyplan:::.check_psu_table(psu_small, small, z$measures))

  wrong <- psu_small
  wrong$N[1] <- wrong$N[1] * 1.001
  expect_error(
    svyplan:::.check_psu_table(wrong, small, z$measures),
    "must equal 'frame\\$N'"
  )
})

test_that("prec_alloc requires one allocation per stratum", {
  z <- .psu_fixture()
  # Recycling a scalar would silently read "this many in every stratum" as
  # though it were a total, so the length is required rather than repaired.
  expect_error(
    prec_alloc(z$frame, n = 300, measures = z$measures, targets = z$targets,
               psu = z$psu),
    "length nrow\\(frame\\)"
  )
  expect_error(
    prec_alloc(z$frame, n = c(300, 200), measures = z$measures,
               targets = z$targets, psu = z$psu),
    "length nrow\\(frame\\)"
  )
  expect_s3_class(
    prec_alloc(z$frame, n = c(300, 200, 130, 80), measures = z$measures,
               targets = z$targets, psu = z$psu),
    "svyplan_prec"
  )
})
