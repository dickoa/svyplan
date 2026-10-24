## The 'psu' argument: certainty-aware allocation from a PSU register

test_that("a PSU register refines the allocation without changing its class", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  expect_s3_class(fit, "svyplan_n")
  expect_identical(fit$type, "alloc")
  expect_identical(fit$method, "bethel")
  expect_true(all(fit$constraints$.pass))
  expect_true(all(fit$operational$constraints$.pass))
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

test_that("the threshold is read at the planned design's own fraction", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)

  # threshold = take / f, where f is the fraction the planned design is the
  # closure at. Under unit cost that design is the closure of the threshold
  # classification itself, and it fields at least f N units.
  f <- z$frame$n_per_psu / fit$detail$threshold
  expect_true(all(fit$detail$n_int >= f * fit$detail$N * (1 - 1e-9)))
  for (h in seq_len(nrow(z$frame))) {
    i <- which(z$psu$stratum == z$frame$stratum[h])
    held <- z$psu$N[i] >= fit$detail$threshold[h]
    closed <- .psu_draw_closure(f[h] * fit$detail$N[h], fit$detail$N[h],
                                z$frame$n_per_psu[h], z$psu$N[i], held)
    expect_identical(closed$certain, fit$psu$certainty[i])
  }
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

  # A PSU that reached the threshold takes at least the stated take at the
  # stratum rate. One made certain by the whole-PSU draw sits below the
  # threshold, so its stratum-rate take can be smaller. No take exceeds
  # the PSU.
  thr <- cert & fit$psu$.certainty_source == "threshold"
  expect_true(any(fit$psu$.certainty_source[cert] == "operational"))
  expect_true(all(fit$psu$n_take[thr] >= take_of[thr]))
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

test_that("the search is reported with its bound and never silent", {
  z <- .psu_fixture()
  for (stage in c(FALSE, TRUE)) {
    if (stage) {
      z$frame$cost_psu <- 100
      z$frame$cost_ssu <- 10
    }
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu)
    sr <- fit$optimization$search

    expect_null(fit$optimization$certainty)
    expect_identical(sr$family, if (stage) "voluntary" else "closure")
    expect_identical(sr$status, "optimal")
    expect_identical(fit$operational$status, sr$status)
    expect_equal(sr$bound, fit$operational$cost)
    expect_identical(sr$gap, 0)
    expect_gte(sr$nodes, 1)
    expect_gte(sr$states, sr$kept)
    expect_identical(sr$max_nodes, 5000)
  }
})

test_that("'control' is checked and kept by every rebuild of a register plan", {
  z <- .psu_fixture(3)
  z$frame$cost_psu <- 400
  z$frame$cost_ssu <- 20
  fit_with <- function(control) {
    n_alloc(z$frame, measures = z$measures, targets = z$targets,
            psu = z$psu, control = control)
  }
  frame <- data.frame(stratum = c("A", "B"), N = c(100, 200), sd = 1)
  expect_error(n_alloc(frame, n = 50, control = list(max_nodes = 10)),
               "'control' applies where n_alloc\\(\\) searches")
  expect_error(fit_with("x"), "'control' must be a named list")
  expect_error(fit_with(list(10)), "'control' must be a named list")
  expect_error(fit_with(list(max_nodes = 10, nodes = 1)),
               "unknown 'control' setting 'nodes'")
  for (bad in list(0, 1.5, NA_real_, Inf, c(1, 2), "10")) {
    expect_error(fit_with(list(max_nodes = bad)),
                 "'control$max_nodes' must be one positive whole number",
                 fixed = TRUE)
  }

  # One node leaves this register at the limit and the default proves it,
  # so a rebuild that dropped 'control' would read "optimal".
  fit <- fit_with(list(max_nodes = 1))
  expect_identical(fit$params$control, list(max_nodes = 1))
  expect_identical(fit$operational$status, "limit")
  expect_identical(fit_with(NULL)$operational$status, "optimal")
  again <- n_alloc(prec_alloc(fit))
  expect_identical(again$optimization$search$max_nodes, 1)
  expect_identical(again$operational$status, "limit")
  over <- prec_alloc(fit, n = fit$detail$n_int)
  expect_identical(over$params$control, list(max_nodes = 1))
  expect_identical(n_alloc(over)$optimization$search$max_nodes, 1)
  expect_identical(predict(fit, data.frame(n_per_psu = 12))$.status, "limit")
})

test_that("a search stopped at its node limit reports the gap", {
  z <- .psu_fixture()
  z$frame$cost_psu <- 100
  z$frame$cost_ssu <- 10
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu, control = list(max_nodes = 1))
  sr <- fit$optimization$search
  expect_identical(sr$max_nodes, 1)
  expect_identical(sr$status, "limit")
  expect_identical(fit$operational$status, "limit")
  expect_lt(sr$bound, fit$operational$cost)
  expect_equal(sr$gap, 1 - sr$bound / fit$operational$cost)
  expect_true(all(fit$operational$constraints$.pass))

  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  out <- capture.output(print(fit))
  expect_true(any(grepl("^search: stopped at the node limit, gap ", out)))
  expect_lte(max(nchar(out)), 80L)
  optimal <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     psu = z$psu)
  expect_false(any(grepl("^search:", capture.output(print(optimal)))))
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

test_that("a PSU smaller than its take is refused unless flagged certainty", {
  psu <- data.frame(
    psu_id = c(sprintf("A%02d", 1:20), sprintf("B%02d", 1:10)),
    stratum = rep(c("A", "B"), c(20, 10)),
    N = c(rep(10, 8), rep(300, 12), rep(200, 10))
  )
  frame <- data.frame(
    stratum = c("A", "B"),
    N = as.numeric(tapply(psu$N, psu$stratum, sum)),
    n_per_psu = 25
  )
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = 0.3,
                         icc_psu = 0.05)
  targets <- data.frame(name = "y", cv = 0.05)

  expect_error(
    n_alloc(frame, measures = measures, targets = targets, psu = psu),
    "8 PSUs in 'psu' hold fewer units than the take 'n_per_psu' of their stratum \\(.A01., .A02., .A03., .A04., .A05. and 3 more\\)"
  )
  expect_error(
    prec_alloc(frame, n = c(300, 200), measures = measures,
               targets = targets, psu = psu),
    "8 PSUs in 'psu' hold fewer units"
  )
  one <- psu[-(2:8), ]
  one$psu_id <- NULL
  one_frame <- frame
  one_frame$N[1] <- sum(one$N[one$stratum == "A"])
  expect_error(
    n_alloc(one_frame, measures = measures, targets = targets, psu = one),
    "1 PSU in 'psu' holds fewer units than the take 'n_per_psu' of its stratum \\(.row 1.\\)"
  )

  # The take is per stratum, and a PSU holding exactly the take can field it.
  per_stratum <- frame
  per_stratum$n_per_psu <- c(10, 25)
  expect_s3_class(
    n_alloc(per_stratum, measures = measures, targets = targets, psu = psu),
    "svyplan_n"
  )

  flagged <- psu
  flagged$certainty <- flagged$N < 25
  fit <- n_alloc(frame, measures = measures, targets = targets, psu = flagged)
  small <- fit$psu$N < 25
  expect_true(all(fit$psu$certainty[small]))
  expect_true(all(fit$psu$n_take[small] <= fit$psu$N[small]))
  expect_true(all(fit$psu$N[!fit$psu$certainty] >= 25))
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
  for (stage in c(FALSE, TRUE)) {
    if (stage) {
      z$frame$cost_psu <- 100
      z$frame$cost_ssu <- 10
    }
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu)
    ev <- prec_alloc(fit)

    # prec_alloc(fit) evaluates the stored field design, planned certainty
    # and takes included, so it reports what the fit's own assessment does.
    expect_equal(ev$detail$.achieved, fit$operational$constraints$.achieved,
                 tolerance = 1e-12)
    expect_identical(ev$psu, fit$psu)
    expect_equal(ev$params$achieved$cost, fit$operational$cost)
  }
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

test_that("supplied certainty is reported as supplied, never as planned", {
  z <- .psu_fixture()
  bare <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                  psu = z$psu)

  flagged <- z$psu
  ratio <- bare$psu$N / bare$psu$.threshold
  flagged$certainty <- ratio > 0.5 & ratio < 0.9
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = flagged)
  n_supplied <- sum(fit$psu$.certainty_source %in% "supplied")
  expect_gt(n_supplied, 0)
  expect_identical(fit$psu$.certainty_source[flagged$certainty &
                                               fit$psu$N < fit$psu$.threshold],
                   rep("supplied", n_supplied))
  # Under unit cost the planner chooses no certainty of its own.
  expect_false("planned" %in% fit$psu$.certainty_source)

  old <- options(width = 80)
  on.exit(options(old), add = TRUE)
  out <- capture.output(print(fit))
  expect_true(any(grepl(sprintf("certainty \\(%d supplied\\),", n_supplied),
                        out)))
  expect_false(any(grepl("supplied", capture.output(print(bare)))))
})

test_that("an allocation with no register prints exactly as before", {
  z <- .bethel_fixture()
  out <- capture.output(print(
    n_alloc(z$frame, measures = z$measures, targets = z$targets)
  ))
  expect_true(any(grepl("^Joint constrained allocation \\(Bethel\\)$", out)))
  expect_false(any(grepl("PSU register|PSUs:|certainty", out)))
})

test_that("summary carries the split by stratum and how the search ended", {
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
  expect_true(any(grepl(
    "^  field design search: optimal after \\d+ nodes? of 5,000, bound ", out
  )))
  expect_true(any(grepl("^  states: [0-9,]+ listed, [0-9,]+ kept by the bound$",
                        out)))

  # The solver saw a one-stage problem, but the plan is not one.
  expect_true(any(grepl("variance model: two_part_certainty_wald", out)))
})


test_that("summary reports the two stages of the plan, not the solve", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  expect_identical(summary(fit)$stages, 2L)
  expect_identical(summary(fit)$assumptions$stages, 2L)
  expect_identical(summary(prec_alloc(fit))$stages, 2L)
  expect_true(any(grepl("^Stages: 2$", capture.output(print(summary(fit))))))

  plain <- n_alloc(
    z$frame[, c("stratum", "N")],
    measures = z$measures[, setdiff(names(z$measures), "icc_psu")],
    targets = z$targets
  )
  expect_identical(summary(plain)$stages, 1L)
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

  # A certainty visit is paid whatever its take, so only the remainder's
  # share of a stratum's interviews brings visits with it.
  rest <- tapply(fit$psu$N * !fit$psu$certainty, fit$psu$stratum, sum)
  share <- as.numeric(rest[fr$stratum]) / fr$N
  expect_equal(fit$params$cost_h, 20 + share * 400 / 12)
  expect_true(all(fit$params$cost_h <= 400 / 12 + 20))
  expect_true(any(fit$params$cost_h < 400 / 12 + 20))
  expect_equal(
    fit$params$achieved$cost,
    sum(d$n_psu_certain * 400 + fit$params$cost_h * d$n)
  )
  expect_equal(fit$optimization$cost, fit$params$achieved$cost)

  # The design is priced as it is fielded: a visit to every PSU entered,
  # certainty or drawn, plus every interview.
  expect_equal(
    fit$operational$cost,
    sum((d$n_psu_certain + d$n_psu_draw) * 400 + d$n_int * 20)
  )

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
  expect_warning(
    g <- predict(fit, data.frame(n_per_psu = c(6, 12, 24, 40))),
    "n_per_psu = 40 is infeasible: 5 PSUs in 'psu' hold fewer units"
  )

  # The smallest PSU of stratum D holds 28, so a take of 40 cannot be fielded.
  expect_equal(nrow(g), 4L)
  expect_identical(g$.feasible, c(TRUE, TRUE, TRUE, FALSE))
  g <- g[g$.feasible, ]
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

test_that("the precision a register fit claims is the one it achieves", .quiet_single_draw({
  # Registers on which an earlier planner reported precision computed at
  # another allocation's design effect.
  for (k in c(4, 25, 53, 60, 74, 89, 99, 115, 119)) {
    z <- .psu_random_register(k)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu)
    expect_equal(prec_alloc(fit)$detail$.achieved,
                 fit$operational$constraints$.achieved,
                 tolerance = 1e-12, label = sprintf("register %d", k))
    expect_true(all(prec_alloc(fit)$detail$.pass),
                label = sprintf("register %d passes", k))
  }
}))

## The ICC the register model expects

.psu_population <- function(seed = 11L) {
  set.seed(seed)
  do.call(rbind, lapply(c("A", "B"), function(h) {
    size <- pmax(round(rlnorm(30, log(60), 1)), 2)
    psu <- rep(sprintf("%s%02d", h, seq_along(size)), size)
    effect <- rep(rnorm(length(size), 0, 4), size)
    data.frame(
      stratum = h, psu = psu,
      share = rep(size / sum(size), size),
      y = 100 + effect + rnorm(sum(size), 0, 20)
    )
  }))
}

.between_share <- function(y, psu) {
  ss_between <- sum((ave(y, psu) - mean(y))^2)
  # var() divides by N_i - 1, so the within sum of squares carries
  # N_i / (N_i - 1). That factor is the whole gap to the plain
  # between-to-total share.
  size <- tapply(y, psu, length)
  ss_within <- tapply(y, psu, function(v) sum((v - mean(v))^2))
  ss_between / (ss_between + sum(size / (size - 1) * ss_within))
}

test_that("a PPS varcomp() returns the size-weighted between-PSU share", {
  pop <- .psu_population()
  one <- pop[pop$stratum == "A", ]
  expect_equal(
    varcomp(y ~ psu, data = one, prob = ~share)$icc,
    .between_share(one$y, one$psu),
    tolerance = 1e-12
  )

  by_stratum <- varcomp(y ~ psu, strata = ~stratum, data = pop, prob = ~share)
  expected <- vapply(split(pop, pop$stratum), function(d) {
    .between_share(d$y, d$psu)
  }, numeric(1))
  expect_equal(
    by_stratum$strata$icc_psu,
    unname(expected[by_stratum$strata$stratum]),
    tolerance = 1e-12
  )
})

test_that("an SRS varcomp() overstates it when PSU sizes are unequal", {
  pop <- .psu_population()
  srs <- varcomp(y ~ psu, strata = ~stratum, data = pop)$strata$icc_psu
  pps <- varcomp(y ~ psu, strata = ~stratum, data = pop,
                 prob = ~share)$strata$icc_psu
  expect_true(all(srs > 2 * pps))
})

.single_psu_case <- function(icc_single) {
  list(
    frame = data.frame(stratum = c("A", "B"), N = c(50000, 2000),
                       n_per_psu = 25),
    psu = rbind(
      data.frame(stratum = "A", N = rep(250, 200)),
      data.frame(stratum = "B", N = 2000)
    ),
    measures = data.frame(stratum = c("A", "B"), name = "y", p = 0.3,
                          icc_psu = c(0.05, icc_single)),
    targets = data.frame(name = "y", cv = 0.12)
  )
}

test_that("a single-PSU stratum carries no between-PSU variance", {
  # Stratum B's allocation is below the take, so its PSU is not certainty by
  # the threshold. Drawing one PSU from one gives it probability one, which
  # makes it certainty by the field design. Its ICC is never charged.
  fits <- lapply(c(NA, 0, 0.2, 0.8), function(icc) {
    z <- .single_psu_case(icc)
    n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = z$psu)
  })
  b <- fits[[2]]$psu$stratum == "B"
  expect_true(fits[[2]]$psu$certainty[b])
  expect_identical(fits[[2]]$psu$.certainty_source[b], "operational")
  for (fit in fits[-1]) {
    expect_equal(fit$detail$n, fits[[1]]$detail$n, tolerance = 1e-12)
    expect_identical(fit$detail$n_int, fits[[1]]$detail$n_int)
  }

  cvs <- vapply(c(NA, 0, 0.8), function(icc) {
    z <- .single_psu_case(icc)
    prec_alloc(z$frame, n = c(330, 15), measures = z$measures,
               targets = z$targets, psu = z$psu)$detail$.achieved
  }, numeric(1))
  expect_equal(cvs, rep(cvs[1], 3), tolerance = 1e-12)
})

test_that("a missing ICC is refused where it would be read", {
  z <- .single_psu_case(0.1)
  z$measures$icc_psu[1] <- NA
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = z$psu),
    "missing in stratum .A."
  )
  expect_error(
    prec_alloc(z$frame, n = c(330, 15), measures = z$measures,
               targets = z$targets, psu = z$psu),
    "missing in stratum .A."
  )
})

test_that("a fit with a missing single-PSU ICC re-solves and re-assesses", {
  z <- .single_psu_case(NA)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  expect_true(is.na(fit$params$measures$icc_psu[2]))
  expect_equal(prec_alloc(fit)$detail$.achieved,
               fit$operational$constraints$.achieved, tolerance = 1e-12)
  swept <- predict(fit, data.frame(n_per_psu = c(15, 25)))
  expect_equal(nrow(swept), 2L)
})

## The fielded draw: no remainder PSU reaches probability one

test_that("closing one PSU under the draw can uncover the next", {
  # Stratum C of a random register: drawing 3 of these PSUs gives the
  # largest a probability of 1.33. Once it is certainty the draw falls to 2
  # and the next reaches 1.07. With both certain the draw is 1 and the
  # largest left is at 0.56.
  sizes <- c(1352, 911, 444, 294, 54)
  out <- .psu_draw_closure(67.27, 3055, 30, sizes, rep(FALSE, 5))
  expect_identical(out$certain, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(out$operational, out$certain)
  # Both sit below the size threshold of 30 / (67.27 / 3055) = 1362.
  expect_true(all(sizes[out$certain] < 30 / (67.27 / 3055)))

  lone <- .psu_draw_closure(20, 2000, 25, 2000, FALSE)
  expect_true(lone$certain)

  held <- .psu_draw_closure(50, 600, 10, c(400, 200), c(TRUE, TRUE))
  expect_identical(held$operational, c(FALSE, FALSE))

  calm <- .psu_draw_closure(40, 4000, 10, rep(400, 10), rep(FALSE, 10))
  expect_false(any(calm$certain))
})

test_that("no fitted plan leaves a remainder PSU at probability one", .quiet_single_draw({
  # samplyr refuses such a plan, so this is the contract that lets a fit go
  # straight to draw().
  cases <- c(
    lapply(1:40, .psu_random_register),
    lapply(c(1, 4, 9, 17), function(s) .psu_fixture(seed = s))
  )
  for (z in cases) {
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu)
    expect_false(any(.psu_remainder_crossing(fit)))
    above <- fit$psu$N >= fit$psu$.threshold
    expect_true(all(fit$psu$certainty[above]))
    op <- fit$psu$.certainty_source %in% "operational"
    expect_true(all(!above[op]))
    expect_true(fit$operational$all_pass)
  }
}))

test_that("the classification closes under the planned draw", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  expect_gt(sum(fit$psu$.certainty_source %in% "operational"), 0L)
  expect_false(any(.psu_remainder_crossing(fit)))
})

test_that("wide registers plan no remainder PSU at probability one", .quiet_single_draw({
  # Registers where an earlier planner's repair added a remainder PSU after
  # the classification was closed and the larger draw reached one.
  for (seed in c(95, 106)) {
    z <- .psu_wide_register(seed)
    for (stage in c(FALSE, TRUE)) {
      if (stage) {
        z$frame$cost_psu <- 150
        z$frame$cost_ssu <- 10
      }
      fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     psu = z$psu)
      expect_false(any(.psu_remainder_crossing(fit)))
      expect_true(all(fit$operational$constraints$.pass))
    }
  }
}))

test_that("prec_alloc classifies a supplied allocation under its draw", {
  frame <- data.frame(stratum = c("C", "D"), N = c(3055, 8000),
                      n_per_psu = 30)
  psu <- data.frame(
    psu_id = 1:25,
    stratum = c(rep("C", 5), rep("D", 20)),
    N = c(1352, 911, 444, 294, 54, rep(400, 20))
  )
  measures <- data.frame(stratum = c("C", "D"), name = "y", p = 0.4,
                         icc_psu = 0.1)
  ev <- prec_alloc(frame, n = c(67.27, 120), measures = measures,
                   targets = data.frame(name = "y", cv = 0.1), psu = psu)
  expect_identical(ev$psu$certainty[1:5], c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(ev$psu$.certainty_source[1:2],
                   c("operational", "operational"))
})

test_that("summary counts the certainty PSUs by source", {
  z <- .psu_fixture()
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  s <- summary(fit)
  expect_named(s$certainty_sources,
               c("threshold", "operational", "supplied", "planned"))
  expect_equal(sum(s$certainty_sources), sum(fit$psu$certainty))
  expect_gt(s$certainty_sources[["operational"]], 0L)
  out <- capture.output(print(s))
  expect_true(any(grepl("^certainty PSUs by source: threshold \\d+, operational \\d+$", out)))
  # The print line stays as it was: the breakdown is summary detail.
  expect_false(any(grepl("operational", capture.output(print(fit)))))
})

## A certainty cutoff below one

.fit_psu <- function(z, ...) {
  n_alloc(z$frame, measures = z$measures, targets = z$targets, psu = z$psu,
          ...)
}

test_that("a cutoff of one is the rule without a cutoff", .quiet_single_draw({
  for (z in c(list(.psu_fixture()), lapply(1:8, .psu_random_register))) {
    plain <- .fit_psu(z)
    one <- .fit_psu(z, certainty_cutoff = 1)
    expect_identical(one$detail, plain$detail)
    expect_identical(one$psu, plain$psu)
  }
}))

test_that("a cutoff relaxes both the threshold and the draw", .quiet_single_draw({
  cases <- c(list(.psu_fixture()), lapply(1:15, .psu_random_register))
  for (cutoff in c(0.9, 0.8, 0.7)) {
    for (z in cases) {
      fit <- .fit_psu(z, certainty_cutoff = cutoff)
      expect_false(any(.psu_remainder_crossing(fit, cutoff)))
      above <- fit$psu$N >= fit$psu$.threshold
      expect_true(all(fit$psu$certainty[above]))
      # The threshold is the cutoff's, at the planned design's fraction,
      # and the design fields at least that fraction.
      f <- cutoff * rep_len(z$frame$n_per_psu, nrow(z$frame)) /
        fit$detail$threshold
      expect_true(all(fit$detail$n_int >= f * fit$detail$N * (1 - 1e-9)))
      expect_true(fit$operational$all_pass)
    }
  }
  fit <- .fit_psu(.psu_fixture(seed = 5L), certainty_cutoff = 0.8)
  expect_gt(sum(fit$psu$.certainty_source %in% "operational"), 0L)
}))

test_that("at a held allocation a lower cutoff makes a superset certain", {
  # n_alloc() re-solves under each cutoff and can land on another cycle
  # resolution, so the superset property is a property of a fixed design.
  z <- .psu_fixture()
  n <- .fit_psu(z)$detail$n
  certain <- lapply(c(1, 0.9, 0.8, 0.7), function(cutoff) {
    prec_alloc(z$frame, n = n, measures = z$measures, targets = z$targets,
               psu = z$psu, certainty_cutoff = cutoff)$psu$certainty
  })
  for (j in 2:4) {
    expect_true(all(certain[[j - 1]] <= certain[[j]]))
  }
  expect_gt(sum(certain[[4]]), sum(certain[[1]]))

  # Each rule applies the cutoff on its own: the threshold is scaled, and no
  # remainder PSU reaches the cutoff under the draw the design implies.
  ev <- prec_alloc(z$frame, n = n, measures = z$measures, targets = z$targets,
                   psu = z$psu, certainty_cutoff = 0.7)
  f <- n / z$frame$N
  expect_equal(unique(ev$psu$.threshold[ev$psu$stratum == "A"]),
               0.7 * z$frame$n_per_psu[1] / f[1])
  for (h in seq_len(nrow(z$frame))) {
    i <- ev$psu$stratum == z$frame$stratum[h]
    counts <- .psu_stratum_counts(n[h], z$frame$N[h], z$frame$n_per_psu[h],
                                  ev$psu$N[i], ev$psu$certainty[i])
    expect_false(any(.psu_crossing(counts$n_psu_draw, ev$psu$N[i],
                                   ev$psu$certainty[i], 0.7)))
  }
})

test_that("the cutoff travels with the fit", {
  z <- .psu_fixture()
  fit <- .fit_psu(z, certainty_cutoff = 0.8)
  ev <- prec_alloc(fit)
  expect_equal(ev$detail$.achieved, fit$operational$constraints$.achieved,
               tolerance = 1e-12)
  expect_equal(ev$psu$.threshold, fit$psu$.threshold)
  expect_identical(ev$params$certainty_cutoff, 0.8)
  take <- z$frame$n_per_psu[1]
  expect_equal(predict(fit, data.frame(n_per_psu = take))$n, fit$n)

  swept <- predict(.fit_psu(z), data.frame(certainty_cutoff = c(1, 0.8)))
  expect_equal(swept$n, c(.fit_psu(z)$n, fit$n))

  in_frame <- z
  in_frame$frame$certainty_cutoff <- 0.8
  expect_equal(.fit_psu(in_frame)$n, fit$n)
  # The argument overrides the column, as deff and resp_rate do.
  expect_equal(.fit_psu(in_frame, certainty_cutoff = 1)$n, .fit_psu(z)$n)
  per_stratum <- .fit_psu(z, certainty_cutoff = rep(0.8, nrow(z$frame)))
  expect_equal(per_stratum$n, fit$n)
})

test_that("summary shows the cutoff only where it relaxes the rule", {
  z <- .psu_fixture()
  cutoff <- c(0.8, 0.8, 1, 0.7)
  s <- summary(.fit_psu(z, certainty_cutoff = cutoff))
  expect_equal(s$certainty$cutoff, cutoff)
  expect_false("cutoff" %in% names(summary(.fit_psu(z))$certainty))
})

test_that("a cutoff is refused outside (0, 1] and without a register", {
  z <- .psu_fixture()
  for (bad in list(0, 1.2, -0.5, NA_real_, c(0.8, 0.9), "0.8")) {
    expect_error(.fit_psu(z, certainty_cutoff = bad), "certainty_cutoff")
  }
  plain <- z$measures[, setdiff(names(z$measures), "icc_psu")]
  expect_error(
    n_alloc(z$frame[, c("stratum", "N")], measures = plain,
            targets = z$targets, certainty_cutoff = 0.8),
    "applies to a PSU register"
  )
  in_frame <- z$frame[, c("stratum", "N")]
  in_frame$certainty_cutoff <- 0.8
  expect_error(
    n_alloc(in_frame, measures = plain, targets = z$targets),
    "applies to a PSU register"
  )
  expect_error(
    prec_alloc(z$frame[, c("stratum", "N")], n = c(300, 200, 130, 80),
               measures = plain, targets = z$targets, certainty_cutoff = 0.8),
    "applies to a PSU register"
  )
  expect_error(
    predict(.fit_psu(z), data.frame(certainty_cutoff = 1.5)),
    "certainty_cutoff"
  )
})

## A register precision result plans its way back

test_that("n_alloc() of a register fit's precision re-plans the register", .quiet_single_draw({
  cases <- c(lapply(1:15, .psu_random_register),
             lapply(c(1, 4, 9, 17, 23), function(s) .psu_fixture(seed = s)))
  for (z in cases) {
    fit <- .fit_psu(z)
    back <- n_alloc(prec_alloc(fit))
    expect_false(is.null(back$params$psu))
    expect_true(all(c("n_psu_certain", "n_psu_draw") %in% names(back$detail)))
    # The re-plan solves the same problem, so it proves the same optimum,
    # though it may pick another plan of that cost.
    expect_identical(back$operational$status, "optimal")
    expect_equal(back$operational$cost, fit$operational$cost)
    # It derives its own classification from the register, so it does not
    # report the fit's certainty PSUs as supplied.
    expect_false("supplied" %in% back$psu$.certainty_source)
  }
}))

test_that("the way back keeps the cutoff", {
  fit <- .fit_psu(.psu_fixture(), certainty_cutoff = 0.8)
  back <- n_alloc(prec_alloc(fit))
  expect_identical(back$params$certainty_cutoff, 0.8)
  expect_equal(back$operational$cost, fit$operational$cost)
})

test_that("a supplied allocation's precision re-plans as a register design", {
  z <- .psu_fixture()
  ev <- prec_alloc(z$frame, n = c(300, 200, 130, 80), measures = z$measures,
                   targets = z$targets, psu = z$psu)
  back <- n_alloc(ev)
  expect_false(is.null(back$params$psu))
  expect_true(all(back$constraints$.pass))
  expect_true(back$operational$all_pass)
})

test_that("pinned targets re-plan to a plan of the same cost", {
  # Every constraint of a round trip binds at its achieved value, which an
  # earlier planner's settle step could not close on this register.
  fit <- .fit_psu(.psu_fixture(seed = 19L))
  back <- n_alloc(prec_alloc(fit))
  expect_true(back$operational$all_pass)
  expect_equal(back$operational$cost, fit$operational$cost)
})

## Zones of the remainder

test_that("zones are cut on cumulative size in the key's order", {
  sizes <- c(50, 40, 30, 20, 10)
  # Largest first: cumulative 50, 90, 120, 140, 150 against a width of 75.
  expect_identical(.psu_zones(4, sizes, rep(FALSE, 5), 2L),
                   c(1L, 2L, 2L, 2L, 2L))
  # A key reverses the order: cumulative 10, 30, 60, 100, 150.
  expect_identical(.psu_zones(4, sizes, rep(FALSE, 5), 2L, key = 5:1),
                   c(2L, 2L, 1L, 1L, 1L))
  # Certainty PSUs have no zone, and no draw means no zones.
  expect_identical(.psu_zones(2, sizes, c(TRUE, rep(FALSE, 4)), 2L),
                   c(NA, 1L, 1L, 1L, 1L))
  expect_true(all(is.na(.psu_zones(0, sizes, rep(FALSE, 5), 2L))))
})

test_that("certainty within zones follows the census and probability rules", {
  # The stratum-wide draw already takes the 50 and the 40, at 4 * 50 / 150
  # and 4 * 40 / 150, and the zones wait until they are certain.
  expect_identical(.psu_crossing(4, c(50, 40, 30, 20, 10), rep(FALSE, 5),
                                 m = 2L),
                   c(TRUE, TRUE, FALSE, FALSE, FALSE))
  # A draw that reaches the whole remainder is a census.
  expect_true(all(.psu_crossing(4, c(100, 90, 80), rep(FALSE, 3), m = 2L)))
  # A zone holding no more PSUs than it draws is a census. Zone 1 holds 49
  # and 48, 97 of its width of 100, though neither reaches one stratum-wide
  # (4 * 49 / 200 = 0.98).
  sizes <- c(49, 48, 45, 30, 28)
  expect_identical(.psu_zones(4, sizes, rep(FALSE, 5), 2L),
                   c(1L, 1L, 2L, 2L, 2L))
  expect_identical(.psu_crossing(4, sizes, rep(FALSE, 5), m = 2L),
                   c(TRUE, TRUE, FALSE, FALSE, FALSE))
  # Without zones the same draw takes neither.
  expect_false(any(.psu_crossing(4, sizes, rep(FALSE, 5))))
})

test_that("no zones is the design without the argument", {
  for (z in c(list(.psu_fixture()), lapply(1:6, .psu_random_register))) {
    plain <- .fit_psu(z)
    none <- .fit_psu(z, n_psu_per_zone = NULL)
    expect_identical(none$detail, plain$detail)
    expect_identical(none$psu, plain$psu)
    expect_false(any(c("n_zone", ".zone") %in%
                       c(names(plain$detail), names(plain$psu))))
  }
})

test_that("a zoned plan draws two PSUs from every zone", {
  cases <- c(lapply(1:30, .psu_random_register),
             lapply(c(1, 4, 9, 17, 23), function(s) .psu_fixture(seed = s)))
  for (z in cases) {
    fit <- .fit_psu(z, n_psu_per_zone = 2)
    d <- fit$detail
    expect_true(all(d$n_psu_draw %% 2 == 0))
    expect_identical(d$n_psu_draw, 2 * d$n_zone)
    rest <- !fit$psu$certainty
    drawn <- fit$psu$stratum %in% d$stratum[d$n_psu_draw > 0]
    expect_identical(is.na(fit$psu$.zone), !(rest & drawn))
    cell <- paste(fit$psu$stratum, fit$psu$.zone)[rest & drawn]
    size <- fit$psu$N[rest & drawn]
    expect_true(all(table(cell) >= 3))
    # Zones run from 1 to n_zone in every stratum with none skipped, which
    # samplyr checks when it fields them.
    for (h in d$stratum[d$n_psu_draw > 0]) {
      zone <- fit$psu$.zone[fit$psu$stratum == h & rest]
      expect_identical(sort(unique(zone)),
                       seq_len(d$n_zone[d$stratum == h]))
    }
    expect_true(all(2 * size / ave(size, cell, FUN = sum) < 1))
    expect_true(fit$operational$all_pass)
  }
})

test_that("the zones add certainty of their own", {
  z <- .psu_random_register(1)
  zoned <- .fit_psu(z, n_psu_per_zone = 2)
  plain <- .fit_psu(z)
  expect_gt(sum(zoned$psu$.certainty_source %in% "operational"),
            sum(plain$psu$.certainty_source %in% "operational"))
})

test_that("zones carry through every route back to a plan", {
  z <- .psu_fixture()
  fit <- .fit_psu(z, n_psu_per_zone = 2)
  expect_identical(prec_alloc(fit)$params$n_psu_per_zone, 2)
  expect_equal(prec_alloc(fit)$detail$.achieved,
               fit$operational$constraints$.achieved, tolerance = 1e-12)
  take <- z$frame$n_per_psu[1]
  expect_equal(predict(fit, data.frame(n_per_psu = take))$n, fit$n)
  back <- n_alloc(prec_alloc(fit))
  expect_identical(back$params$n_psu_per_zone, 2)
  expect_equal(back$operational$cost, fit$operational$cost)
  expect_identical(summary(fit)$certainty$zones, fit$detail$n_zone)
})

test_that("one PSU per zone draws one PSU from every zone", {
  cases <- c(lapply(1:30, .psu_random_register),
             lapply(c(1, 4, 9, 17, 23), function(s) .psu_fixture(seed = s)))
  for (z in cases) {
    fit <- .fit_psu(z, n_psu_per_zone = 1)
    d <- fit$detail
    expect_identical(d$n_psu_draw, d$n_zone)
    rest <- !fit$psu$certainty
    drawn <- fit$psu$stratum %in% d$stratum[d$n_psu_draw > 0]
    expect_identical(is.na(fit$psu$.zone), !(rest & drawn))
    for (h in d$stratum[d$n_psu_draw > 0]) {
      zone <- fit$psu$.zone[fit$psu$stratum == h & rest]
      expect_identical(sort(unique(zone)),
                       seq_len(d$n_zone[d$stratum == h]))
    }
    # A zone of one PSU is a census, so every zone left holds two or more
    # and none of them reaches one within its zone.
    cell <- paste(fit$psu$stratum, fit$psu$.zone)[rest & drawn]
    size <- fit$psu$N[rest & drawn]
    expect_true(all(table(cell) >= 2))
    expect_true(all(size / ave(size, cell, FUN = sum) < 1))
    expect_true(fit$operational$all_pass)
  }
})

test_that("a zone of one PSU is a census at one PSU per zone", {
  # A width of 100 / 3: zones 1 and 2 hold one PSU each, though neither PSU
  # reaches one stratum-wide (3 * 30 / 100 = 0.9).
  sizes <- c(30, 30, 20, 20)
  expect_identical(.psu_zones(3, sizes, rep(FALSE, 4), 1L),
                   c(1L, 2L, 3L, 3L))
  expect_identical(.psu_crossing(3, sizes, rep(FALSE, 4), m = 1L),
                   c(TRUE, TRUE, FALSE, FALSE))
  expect_false(any(.psu_crossing(3, sizes, rep(FALSE, 4))))
  # One PSU per zone does not round the draw.
  expect_identical(.psu_round_draw(3, 1L), 3)
})

test_that("each stratum meets its own target on its field takes", {
  # Stratum B needs exactly 12,000 units at its stratum rate and has no
  # clustering. C's rounded takes leave it short at the continuous
  # allocation, so its field design differs from that allocation.
  set.seed(36)
  size <- round(rlnorm(120, log(120), 1.1)) + 15
  frame <- data.frame(stratum = c("B", "C"), N = c(40000, sum(size)),
                      n_per_psu = 10)
  psu <- data.frame(stratum = c(rep("B", 2000), rep("C", length(size))),
                    N = c(rep(20, 2000), size))
  measures <- data.frame(stratum = c("B", "C"), name = "y", p = c(0.5, 0.3),
                         icc_psu = c(0, 0.1))
  targets <- data.frame(name = "y", domain = "stratum", level = c("B", "C"),
                        cv = c(sqrt(1 / 12000 - 1 / 40000), 0.02))
  fit <- expect_no_warning(
    n_alloc(frame, measures = measures, targets = targets, psu = psu)
  )
  expect_identical(fit$detail$n_psu_draw[1], 1200)
  expect_true(all(fit$operational$constraints$.pass))
  expect_true(predict(fit, data.frame(n_per_psu = 10))$.feasible)
})

test_that("a grid row reads feasibility from the field design", {
  fit <- list(operational = list(all_pass = FALSE, n = 10, cost = 10),
              detail = data.frame(n_psu_draw = 1), params = list())
  expect_false(.bethel_predict_row(fit, c("n_int", ".feasible"))$.feasible)
  fit$operational$all_pass <- TRUE
  expect_true(.bethel_predict_row(fit, c("n_int", ".feasible"))$.feasible)
})

test_that("zones are grouped in twos and threes for the variance", {
  zone <- c(1L, 2L, 3L, 4L, 5L, 1L, 1L, 2L, 1L, 1L, NA)
  idx <- list(1:5, 6L, 7:8, 9L, 10L, 11L)
  # Zones pair in order and an odd count ends in a triple. Strata with a
  # single zone are grouped together in frame order.
  expect_identical(.psu_pairs(zone, idx),
                   c(1L, 1L, 2L, 2L, 2L, 3L, 4L, 4L, 3L, 3L, NA))
  # A lone single-zone stratum joins the next stratum with zones, or the
  # previous one when it comes last, and no group exceeds three.
  expect_identical(.psu_pairs(c(1L, 1:4), list(1L, 2:5)),
                   c(1L, 1L, 2L, 2L, 2L))
  expect_identical(.psu_pairs(c(1:3, 1L, NA), list(1:3, 4L, 5L)),
                   c(1L, 1L, 2L, 2L, NA))
  # A design with a single zone has no partner for it.
  expect_identical(.psu_pairs(c(1L, NA), list(1L, 2L)), c(1L, NA))
})

test_that("every zone of a one-per-zone plan has a variance partner", {
  cases <- c(lapply(1:30, .psu_random_register),
             lapply(c(1, 4, 9, 17, 23), function(s) .psu_fixture(seed = s)))
  for (z in cases) {
    fit <- .fit_psu(z, n_psu_per_zone = 1)
    p <- fit$psu
    expect_identical(is.na(p$.pair), is.na(p$.zone))
    zones <- unique(p[!is.na(p$.zone), c("stratum", ".zone", ".pair")])
    if (nrow(zones) == 0L) next
    expect_identical(sort(unique(zones$.pair)),
                     seq_len(max(zones$.pair)))
    size <- table(zones$.pair)
    if (nrow(zones) == 1L) {
      expect_identical(as.integer(size), 1L)
    } else {
      expect_true(all(size %in% 2:3))
    }
  }
  # Two PSUs per zone need no collapsing.
  expect_null(.fit_psu(.psu_fixture(), n_psu_per_zone = 2)$psu$.pair)
})

test_that("one PSU per zone carries through every route back to a plan", {
  z <- .psu_fixture()
  fit <- .fit_psu(z, n_psu_per_zone = 1)
  assessed <- prec_alloc(fit)
  expect_identical(assessed$params$n_psu_per_zone, 1)
  expect_identical(assessed$psu$.zone, fit$psu$.zone)
  expect_identical(assessed$psu$.pair, fit$psu$.pair)
  expect_equal(assessed$detail$.achieved,
               fit$operational$constraints$.achieved, tolerance = 1e-12)
  back <- n_alloc(assessed)
  expect_identical(back$params$n_psu_per_zone, 1)
  expect_equal(back$operational$cost, fit$operational$cost)
})

test_that("predict() sweeps the zone setting, NA meaning no zones", {
  z <- .psu_fixture()
  fit <- .fit_psu(z)
  g <- predict(fit, data.frame(n_psu_per_zone = c(NA, 1, 2)))
  expect_identical(g$.feasible, rep(TRUE, 3))
  for (r in seq_len(nrow(g))) {
    m <- g$n_psu_per_zone[r]
    direct <- .fit_psu(z, n_psu_per_zone = if (is.na(m)) NULL else m)
    expect_equal(g$n[r], direct$n, tolerance = 1e-8)
    expect_identical(g$n_int[r], direct$operational$n)
    expect_equal(g$n_psu_certain[r], sum(direct$psu$certainty))
    expect_equal(g$n_psu_draw[r], sum(direct$detail$n_psu_draw))
    expect_equal(g$n_zone[r], sum(direct$detail$n_zone %||% 0))
  }
  # A zoned fit keeps its zones when the column is absent.
  zoned <- .fit_psu(z, n_psu_per_zone = 1)
  kept <- predict(zoned, data.frame(n_per_psu = z$frame$n_per_psu[1]))
  expect_equal(kept$n_zone, sum(zoned$detail$n_zone))
  expect_null(predict(fit, data.frame(n_per_psu = 12))$n_zone)
  for (bad in list(0, 3, 1.5)) {
    expect_error(predict(fit, data.frame(n_psu_per_zone = bad)),
                 "must contain NA, 1 or 2")
  }
})

test_that("zone_order sets the cut and is refused when it cannot", {
  z <- .psu_fixture()
  by_size <- .fit_psu(z, n_psu_per_zone = 2)
  ordered <- z
  ordered$psu$zone_order <- seq_len(nrow(z$psu))
  by_key <- .fit_psu(ordered, n_psu_per_zone = 2)
  expect_false(identical(by_key$psu$.zone, by_size$psu$.zone))

  expect_error(.fit_psu(ordered), "zone_order")
  missing_key <- ordered
  missing_key$psu$zone_order[3] <- NA
  expect_error(.fit_psu(missing_key, n_psu_per_zone = 2), "missing")
  tied <- ordered
  tied$psu$zone_order[2] <- tied$psu$zone_order[1]
  expect_error(.fit_psu(tied, n_psu_per_zone = 2), "without ties")
  for (bad in list(0, 1.5, 3, NA_real_, "2", c(2, 2))) {
    expect_error(.fit_psu(z, n_psu_per_zone = bad), "must be NULL, 1 or 2")
  }
  plain <- z$measures[, setdiff(names(z$measures), "icc_psu")]
  expect_error(
    n_alloc(z$frame[, c("stratum", "N")], measures = plain,
            targets = z$targets, n_psu_per_zone = 2),
    "applies to a PSU register"
  )
})

## What the register model does not carry

test_that("a register fit sized to targets has no budget to vary", {
  z <- .psu_fixture()
  fit <- .fit_psu(z)
  expect_error(predict(fit, data.frame(budget = 1e5)),
               "unknown parameter(s) in newdata", fixed = TRUE)
})

test_that("a variance ratio other than one is refused with a register", {
  z <- .psu_fixture()
  base <- .fit_psu(z)
  for (one in list(1, NA_real_)) {
    m <- z$measures
    m$var_ratio_psu <- one
    expect_equal(
      n_alloc(z$frame, measures = m, targets = z$targets, psu = z$psu)$n,
      base$n
    )
  }
  for (bad in c(2, -1, 0.5)) {
    m <- z$measures
    m$var_ratio_psu <- bad
    expect_error(
      n_alloc(z$frame, measures = m, targets = z$targets, psu = z$psu),
      "'var_ratio_psu' is not available with a PSU register"
    )
    fr <- z$frame
    fr$var_ratio_psu <- bad
    expect_error(
      n_alloc(fr, measures = z$measures, targets = z$targets, psu = z$psu),
      "'var_ratio_psu' is not available with a PSU register"
    )
    expect_error(
      prec_alloc(z$frame, n = base$detail$n, measures = m,
                 targets = z$targets, psu = z$psu),
      "'var_ratio_psu' is not available with a PSU register"
    )
  }
})

test_that("PSU-level response is refused with a register by its own name", {
  z <- .psu_fixture()
  m <- z$measures
  m$resp_rate_psu <- 0.9
  expect_error(
    n_alloc(z$frame, measures = m, targets = z$targets, psu = z$psu),
    "'resp_rate_psu' and 'resp_rate_ssu' are not available with a PSU register"
  )
})

## Costs, response and the field design

test_that("the continuous optimum prices certainty visits as fixed", {
  # Stratum A is all certainty, so its interviews cost cost_ssu alone, and
  # B's bring a tenth of a visit each. The optimum of sum(a_h / n_h) under
  # those costs is derived here independently of the solver.
  frame <- data.frame(stratum = c("A", "B"), N = 10000, n_per_psu = 10,
                      cost_psu = 100, cost_ssu = 1)
  psu <- data.frame(stratum = rep(c("A", "B"), each = 100), N = 100,
                    certainty = rep(c(TRUE, FALSE), each = 100))
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = 0.5,
                         icc_psu = 0.05)
  targets <- data.frame(name = "y", cv = 0.1)
  fit <- n_alloc(frame, measures = measures, targets = targets, psu = psu)

  a <- c(0.0625, 0.090625)
  cost <- c(1, 11)
  limit <- (0.5 * 0.1)^2 + sum(a / 10000)
  optimum <- sqrt(a / cost) * sum(sqrt(a * cost)) / limit
  expect_equal(fit$detail$n, optimum, tolerance = 1e-6)
  expect_equal(fit$params$cost_h, cost)
  expect_equal(fit$params$achieved$cost, 10000 + sum(cost * optimum),
               tolerance = 1e-8)
  # prec_alloc(fit) prices the field design, as the fit's own cost does.
  expect_equal(prec_alloc(fit)$params$achieved$cost, fit$operational$cost)
})

test_that("the register clusters on the responding take", {
  # At 300 units no PSU reaches the threshold, so the register design is
  # the ordinary two-stage one and both cluster on take times response.
  frame <- data.frame(stratum = "A", N = 100000, n_per_psu = 10,
                      cost_psu = 100, cost_ssu = 10)
  psu <- data.frame(stratum = "A", N = rep(1000, 100))
  measures <- data.frame(stratum = "A", name = "y", p = 0.5, icc_psu = 0.05)
  targets <- data.frame(name = "y", cv = 0.08)
  for (rate in c(1, 0.5)) {
    register <- prec_alloc(frame, n = 300, measures = measures,
                           targets = targets, psu = psu, resp_rate = rate)
    ordinary <- prec_alloc(transform(frame, N_psu = 100), n = 300,
                           measures = measures, targets = targets,
                           resp_rate = rate)
    expect_false(any(register$psu$certainty))
    expect_equal(register$detail$.achieved, ordinary$detail$.achieved,
                 tolerance = 1e-12)
    by_row <- prec_alloc(frame, n = 300,
                         measures = transform(measures, resp_rate = rate),
                         targets = targets, psu = psu)
    expect_equal(by_row$detail$.achieved, register$detail$.achieved,
                 tolerance = 1e-12)
  }
})

test_that("the field design effect is the stratum one at the stratum rate", {
  sizes <- c(400, 300, 50, 60, 70, 80, 40)
  certain <- c(TRUE, TRUE, rep(FALSE, 5))
  f <- 0.1
  for (resp in c(1, 0.6)) {
    field <- .psu_field_deff(sizes, certain, f * sizes[certain],
                             f * sum(sizes[!certain]) / 5, 5, 0.1, resp)
    expect_equal(
      field,
      .psu_deff(sum(sizes), 700, 300, 0.1, 5 * resp)
    )
  }
  # A census of every certainty PSU at full response carries no variance
  # from that part.
  census <- .psu_field_deff(sizes, certain, sizes[certain], 6, 5, 0.1, 1)
  rest <- 300^2 * (1 + 0.1 * 4) * (1 / 30 - 1 / 300)
  expect_equal(census, max(rest / (1000^2 * (1 / 730 - 1 / 1000)), 1))
})

test_that("the field design is assessed on its own takes", .quiet_single_draw({
  # The certainty takes are rounded up one by one and the remainder gets
  # what is left, far below the stratum rate on this register.
  set.seed(36)
  size <- round(rlnorm(120, log(120), 1.1)) + 15
  psu <- data.frame(stratum = "A", N = size)
  frame <- data.frame(stratum = "A", N = sum(size), n_per_psu = 10,
                      cost_psu = 200, cost_ssu = 20)
  measures <- data.frame(stratum = "A", name = "y", p = 0.3, icc_psu = 0.1)
  fit <- n_alloc(frame, measures = measures,
                 targets = data.frame(name = "y", cv = 0.02), psu = psu)

  held <- fit$psu$certainty
  size_c <- fit$psu$N[held]
  size_r <- sum(fit$psu$N[!held])
  draws <- fit$detail$n_psu_draw
  v <- sum(size_c^2 * (1 / fit$psu$n_take[held] - 1 / size_c)) +
    size_r^2 * (1 + 0.1 * 9) * (1 / (10 * draws) - 1 / size_r)
  cv <- sqrt(0.3 * 0.7 * v) / (0.3 * sum(size))
  expect_equal(fit$operational$constraints$.cv, cv, tolerance = 1e-10)
  expect_lte(cv, 0.02)
  expect_true(fit$operational$all_pass)
}))

test_that("no remainder is left without a draw", .quiet_single_draw({
  # On this register the rounded-up certainty takes exceed the allocation
  # in two strata, each with a remainder of one PSU.
  for (seed in c(10, 52, 58)) {
    z <- .psu_wide_register(seed)
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu)
    d <- fit$detail
    expect_true(all(d$n_psu_draw >= 1 | d$n_psu_rest == 0),
                label = sprintf("register %d", seed))
  }
  counts <- .psu_stratum_counts(100, 1000, 10, c(600, 300, 100),
                                c(TRUE, TRUE, FALSE))
  expect_identical(counts$n_rest, 10)
  expect_identical(counts$n_psu_draw, 1)
  over <- .psu_stratum_counts(99, 1000, 10, c(rep(95, 10), 50),
                              c(rep(TRUE, 10), FALSE))
  expect_identical(over$n_certain_int, 100)
  expect_identical(over$n_rest, 0)
  expect_identical(over$n_psu_draw, 1)
}))

test_that("summary stratum costs add up to the plan's cost", {
  frame <- data.frame(stratum = c("A", "B"), N = 10000, n_per_psu = 10,
                      cost_psu = 100, cost_ssu = 1)
  psu <- data.frame(stratum = rep(c("A", "B"), each = 100), N = 100,
                    certainty = rep(c(TRUE, FALSE), each = 100))
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = 0.5,
                         icc_psu = 0.05)
  two <- n_alloc(frame, measures = measures,
                 targets = data.frame(name = "y", cv = 0.1), psu = psu)
  d <- two$detail
  expect_equal(summary(two)$allocation$cost,
               c(100 * 100 + d$n_int[1], d$n_psu_draw[2] * 100 + d$n_int[2]))

  # Strata mixing certainty and remainder PSUs, the takes rounded apart.
  z <- .psu_fixture()
  fr <- z$frame
  fr$cost_psu <- 400
  fr$cost_ssu <- 20
  mixed <- n_alloc(fr, measures = z$measures, targets = z$targets,
                   psu = z$psu)
  expect_true(any(mixed$detail$n_psu_certain > 0 &
                    mixed$detail$n_psu_draw > 0))
  for (fit in list(two, mixed)) {
    s <- summary(fit)
    expect_equal(sum(s$allocation$cost), fit$operational$cost)
    expect_equal(sum(s$allocation$cost), s$overall$cost)
    assessed <- summary(prec_alloc(fit))
    expect_equal(sum(assessed$allocation$cost), fit$operational$cost)
    expect_equal(sum(assessed$allocation$cost), assessed$overall$cost)
  }
})

## A remainder that draws one PSU

test_that("a stratum drawing one PSU from its remainder warns by name", {
  register <- data.frame(
    psu_id = 1:24,
    stratum = rep(c("Large", "Small"), c(16, 8)),
    N = c(rep(c(3000, 4000, 5000, 6000), 4), rep(c(400, 500, 600, 700), 2))
  )
  frame <- data.frame(
    stratum = c("Large", "Small"),
    N = c(sum(register$N[1:16]), sum(register$N[17:24])),
    n_per_psu = 10
  )
  measures <- data.frame(stratum = c("Large", "Small"), name = "y",
                         p = c(0.3, 0.32), icc_psu = 0.05)
  targets <- data.frame(name = "y", cv = 0.12)
  plan <- function(...) {
    n_alloc(frame, measures = measures, targets = targets, psu = register,
            ...)
  }

  expect_warning(fit <- plan(), class = "svyplan_single_draw")
  expect_warning(plan(), "stratum 'Small' draws a single PSU")
  expect_warning(plan(), "n_psu_per_zone = 2")
  expect_identical(fit$detail$n_psu_draw[fit$detail$stratum == "Small"], 1)
  expect_gt(fit$detail$n_psu_draw[fit$detail$stratum == "Large"], 1)

  # Zones of two draw at least two, and zones of one pair the single draws
  # for variance estimation, so neither warns.
  expect_no_warning(two <- plan(n_psu_per_zone = 2))
  expect_true(all(two$detail$n_psu_draw != 1))
  expect_no_warning(one <- plan(n_psu_per_zone = 1))
  expect_true(any(one$detail$n_psu_draw == 1))
})
