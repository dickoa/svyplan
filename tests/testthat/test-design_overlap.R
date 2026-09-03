## O1. The six published overlap figures

test_that("CPS 4-8-4 reproduces 75 percent consecutive and 50 percent annual", {
  cps <- design_overlap("4-8-4")
  expect_s3_class(cps, "svyplan_overlap")
  expect_equal(cps[1], 0.75, tolerance = 1e-12)
  expect_equal(cps[12], 0.50, tolerance = 1e-12)
})

test_that("ONS in-for-15 gives 93 percent and 20 percent", {
  o <- design_overlap("15")
  expect_equal(o[1], 14 / 15, tolerance = 1e-12)
  expect_equal(o[12], 3 / 15, tolerance = 1e-12)
})

test_that("ONS in-for-27 gives 96 percent and 56 percent", {
  o <- design_overlap("27")
  expect_equal(o[1], 26 / 27, tolerance = 1e-12)
  expect_equal(o[12], 15 / 27, tolerance = 1e-12)
})

test_that("an unbroken panel drops one cohort per occasion", {
  expect_equal(design_overlap("5")[1], 0.8, tolerance = 1e-12)
  for (k in c(2L, 4L, 5L, 9L, 20L)) {
    o <- design_overlap(as.character(k))
    expect_equal(o[1], (k - 1) / k, tolerance = 1e-12)
  }
})

## O2. The cycling trap

test_that("a two-spell life is not a repeating cycle", {
  # Read as a cycle, only the cohort finishing its stint leaves and 7 of 8
  # are retained. The finite life loses two cohorts, one at the end of each
  # spell, which is the published 6 of 8.
  cps <- design_overlap("4-8-4")
  expect_equal(cps[1], 6 / 8, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(cps[1], 7 / 8)))
  expect_equal(cps$shared[1L], 6)
  expect_equal(cps$n_occasion, 8)
})

## O3. The arithmetic behind the fractions

test_that("shared counts the cohorts in sample at both stages", {
  # Schedule 1,1,0,1 is in sample at stages 1, 2 and 4, so each lag has
  # exactly one qualifying pair: (1,2) at lag 1, (2,4) at lag 2, (1,4) at
  # lag 3. The pair (1,3) does not count, stage 3 being out of sample.
  o <- design_overlap(c(1, 1, 0, 1))
  expect_equal(o$n_occasion, 3)
  expect_equal(o$shared, c(1, 1, 1))
  expect_equal(as.double(o), rep(1 / 3, 3L), tolerance = 1e-12)
})

test_that("the lag profile is symmetric about the gap for two equal spells", {
  cps <- design_overlap("4-8-4")
  expect_equal(as.double(cps),
               c(6, 4, 2, 0, 0, 0, 0, 0, 1, 2, 3, 4, 3, 2, 1) / 8,
               tolerance = 1e-12)
})

test_that("overlap is zero past the reach of the rotation", {
  o <- design_overlap("5", max_lag = 8)
  expect_equal(as.double(o)[5:8], rep(0, 4L))
})

test_that("a subsampled later wave is nested inside the earlier one", {
  o <- design_overlap(c(1, 1, 0.5, 0.5))
  expect_equal(o$n_occasion, 3)
  expect_equal(o$shared, c(2, 1, 0.5), tolerance = 1e-12)
  expect_equal(o[1], 2 / 3, tolerance = 1e-12)
})

test_that("logical and numeric take vectors agree", {
  a <- design_overlap(c(TRUE, TRUE, FALSE, TRUE))
  b <- design_overlap(c(1, 1, 0, 1))
  expect_equal(as.double(a), as.double(b), tolerance = 1e-12)
})

test_that("the compact string expands to the explicit vector", {
  expect_equal(as.double(design_overlap("4-8-4")$rotation),
               c(rep(1, 4), rep(0, 8), rep(1, 4)))
  expect_equal(as.double(design_overlap("3")$rotation), rep(1, 3))
})

## O4. max_lag

test_that("max_lag defaults to one less than the life", {
  expect_length(design_overlap("4-8-4"), 15L)
  expect_length(design_overlap("5"), 4L)
  expect_length(design_overlap("15"), 14L)
})

test_that("max_lag can be set explicitly", {
  expect_length(design_overlap("4-8-4", max_lag = 12), 12L)
  expect_equal(design_overlap("4-8-4", max_lag = 12)[12], 0.5,
               tolerance = 1e-12)
  expect_error(design_overlap("5", max_lag = 0), "whole number >= 1")
  expect_error(design_overlap("5", max_lag = 2.5), "whole number >= 1")
})

## O5. Feeding the overlap arguments

test_that("subsetting gives a plain number the overlap arguments take", {
  cps <- design_overlap("4-8-4")
  expect_type(cps[1], "double")
  expect_false(inherits(cps[1], "svyplan_overlap"))
  expect_length(cps[1], 1L)
  expect_length(cps[[1]], 1L)
  expect_null(names(cps[[12]]))
})

test_that("a lag drives n_change and prec_change", {
  cps <- design_overlap("4-8-4")
  a <- n_change(p = c(0.30, 0.36), moe = 0.02, overlap = cps[1],
                overlap_cor = 0.5)
  b <- n_change(p = c(0.30, 0.36), moe = 0.02, overlap = 0.75,
                overlap_cor = 0.5)
  expect_equal(a$n, b$n, tolerance = 1e-12)
  expect_equal(a$params$overlap, 0.75, tolerance = 1e-12)

  # The annual lag is a different design question and a different size.
  annual <- n_change(p = c(0.30, 0.36), moe = 0.02, overlap = cps[12],
                     overlap_cor = 0.5)
  expect_gt(annual$n, a$n)
})

test_that("a lag drives the power family unchanged", {
  cps <- design_overlap("4-8-4")
  a <- power_mean(var = 100, n = 900, power = 0.5, overlap = cps[1],
                  overlap_cor = 0.6)
  b <- power_mean(var = 100, n = 900, power = 0.5, overlap = 0.75,
                  overlap_cor = 0.6)
  expect_equal(a$effect, b$effect, tolerance = 1e-12)
})

## O6. Object surface

test_that("as.double strips the counts it carries as attributes", {
  cps <- design_overlap("4-8-4")
  expect_null(attributes(as.double(cps)))
  expect_length(as.double(cps), 15L)
})

test_that("as.data.frame reports one row per lag", {
  tab <- as.data.frame(design_overlap("4-8-4"))
  expect_equal(nrow(tab), 15L)
  expect_named(tab, c("lag", "shared", "n_occasion", "overlap"))
  expect_equal(tab$lag, 1:15)
  expect_equal(tab$overlap, tab$shared / tab$n_occasion, tolerance = 1e-12)
})

test_that("$ reaches the fields and rejects unknown names", {
  cps <- design_overlap("4-8-4")
  expect_equal(cps$life, 16L)
  expect_equal(cps$n_occasion, 8)
  expect_equal(cps$lag, 1:15)
  expect_length(cps$rotation, 16L)
  expect_error(cps$nope, "no field 'nope'")
})

test_that("print reports the life, the take and the spells", {
  out <- capture.output(print(design_overlap("4-8-4")))
  expect_true(any(grepl("16-occasion life", out)))
  expect_true(any(grepl("8 in sample each occasion", out)))
  expect_true(any(grepl("4 in, 8 out, 4 in", out)))
  expect_true(any(grepl("0.75", out)))
})

test_that("print truncates a long lag profile and says so", {
  out <- capture.output(print(design_overlap("27")))
  expect_true(any(grepl("further lag", out)))
})

test_that("format is a one-line summary", {
  expect_match(format(design_overlap("4-8-4")), "^svyplan_overlap \\[life 16")
})

test_that("the methods reject unused arguments", {
  cps <- design_overlap("4-8-4")
  expect_error(print(cps, nope = 1), "unused argument")
  expect_error(as.double(cps, nope = 1), "unused argument")
  expect_error(as.data.frame(cps, nope = 1), "unused argument")
})

## O8. Transformations do not leave an invalid overlap behind

test_that("arithmetic returns bare numerics, not a contradicted object", {
  cps <- design_overlap("4-8-4")
  y <- cps * 2
  expect_false(inherits(y, "svyplan_overlap"))
  expect_null(attributes(y))
  expect_equal(y[1L], 1.5, tolerance = 1e-12)
  # Keeping the class would leave an overlap of 1.5 still claiming the
  # rotation and the shared counts it no longer follows from.
  expect_false(inherits(cps / 2, "svyplan_overlap"))
  expect_false(inherits(1 - cps, "svyplan_overlap"))
  expect_false(inherits(-cps, "svyplan_overlap"))
  expect_type(cps > 0.5, "logical")
  expect_null(attributes(cps > 0.5))
  expect_equal(sum(cps > 0.5), 1L)
})

test_that("Math returns bare numerics", {
  cps <- design_overlap("4-8-4")
  expect_false(inherits(sqrt(cps), "svyplan_overlap"))
  expect_null(attributes(round(cps, 2)))
  expect_equal(sqrt(cps)[1L], sqrt(0.75), tolerance = 1e-12)
})

test_that("an empty subscript strips the structure it no longer describes", {
  cps <- design_overlap("4-8-4")
  expect_null(attributes(cps[]))
  expect_equal(cps[], as.double(cps), tolerance = 1e-12)
  expect_null(attributes(cps[1:3]))
  # Indexing by the lag label still works, and still comes back bare.
  expect_equal(cps[["12"]], 0.5, tolerance = 1e-12)
  expect_equal(cps["12"], 0.5, tolerance = 1e-12)
  expect_null(names(cps["12"]))
})

## O10. In-place modification is refused

test_that("replacement is an error rather than a silent contradiction", {
  x <- design_overlap("4-8-4")
  expect_error({x[1] <- 0.2}, "cannot be modified in place")
  expect_error({x[[1]] <- 0.2}, "cannot be modified in place")
  expect_error({x[1:2] <- 0.2}, "cannot be modified in place")
  expect_error({x[["12"]] <- 0.2}, "cannot be modified in place")
  # The message names the way back to the numbers.
  expect_error({x[1] <- 0.2}, "as.double")
  # x itself is untouched.
  expect_equal(x[1], 0.75, tolerance = 1e-12)
  expect_equal(attr(x, "shared")[1L], 6)
})

test_that("pmax and pmin cannot mislabel a transformed overlap", {
  # They copy the attributes of their first argument without dispatching, so
  # the only thing standing between them and a contradicted object is that
  # they assign through `[<-`.
  x <- design_overlap("4-8-4")
  expect_error(pmax(x, 0.8), "cannot be modified in place")
  expect_error(pmin(x, 0.5), "cannot be modified in place")
  expect_equal(pmax(as.double(x), 0.8)[1L], 0.8, tolerance = 1e-12)
  expect_null(attributes(pmax(as.double(x), 0.8)))
})

test_that("the other reshaping routes already return bare numerics", {
  x <- design_overlap("4-8-4")
  for (v in list(rev(x), c(x, 1), max(x), sum(x), range(x), sort(as.double(x)))) {
    expect_false(inherits(v, "svyplan_overlap"))
  }
})

## O11. A profile is not a lag

test_that("passing the whole profile as overlap names the fix", {
  cps <- design_overlap("4-8-4")
  for (call in list(
    function() n_change(p = c(0.3, 0.36), moe = 0.02, overlap = cps),
    function() prec_change(var = 100, n = 500, overlap = cps),
    function() power_mean(var = 100, n = 900, power = 0.5, overlap = cps),
    function() power_prop(p1 = 0.3, p2 = 0.36, power = 0.8, overlap = cps),
    function() power_did(c(1, 2), c(1, 1), outcome = "mean", var = 1,
                         effect = 0.1, overlap = cps)
  )) {
    expect_error(call(), "'overlap' is one lag")
    expect_error(call(), "design_overlap")
  }
  # A profile passed whole carries no lag, and defaulting to the consecutive
  # figure would answer an annual-lag question with the monthly number.
  expect_error(
    n_change(p = c(0.3, 0.36), moe = 0.02, overlap = design_overlap("2")),
    "'overlap' is one lag"
  )
  expect_silent(
    n_change(p = c(0.3, 0.36), moe = 0.02, overlap = cps[1], overlap_cor = 0.5)
  )
})

## O12. The per-occasion notation, and the string both notations claim

test_that("a spec containing a 0 is one flag per occasion", {
  # Lynn (2012) Figure 3: in for two, out for two, in for two
  lynn <- design_overlap("1-1-0-0-1-1")
  expect_equal(as.double(attr(lynn, "rotation")), c(1, 1, 0, 0, 1, 1))
  expect_identical(as.double(lynn), as.double(design_overlap("2-2-2")))
  expect_identical(as.double(lynn), as.double(design_overlap(c(1, 1, 0, 0, 1, 1))))
})

test_that("the Lynn 1-1-0-0-1-1 lag profile matches the paper's own reading", {
  # Lynn (2012), page 7, on the design of Figure 3: change is estimable
  # between periods 3 and 4, 3 and 6, 3 and 7 and 3 and 8, but not between
  # 3 and 5 or 3 and 9, and half the units of period 3 are in period 7.
  x <- design_overlap("1-1-0-0-1-1")
  expect_equal(as.double(x), c(0.5, 0, 0.25, 0.5, 0.25), tolerance = 1e-12)
  expect_equal(x[4], 0.5, tolerance = 1e-12)   # periods 3 and 7
  expect_equal(x[2], 0, tolerance = 1e-12)     # periods 3 and 5
  expect_true(is.na(x[6]))                     # periods 3 and 9, past the life
})

test_that("Lynn Figure 2, three consecutive waves, is the spell spec 3", {
  # one third of the sample replaced each period
  x <- design_overlap("3")
  expect_equal(as.double(x), c(2 / 3, 1 / 3), tolerance = 1e-12)
})

test_that("spell notation is untouched by the pattern reading", {
  expect_equal(design_overlap("4-8-4")[1], 0.75, tolerance = 1e-12)
  expect_error(design_overlap("4-8"), "repeating cycle")
  expect_equal(as.double(design_overlap("5")), c(0.8, 0.6, 0.4, 0.2),
               tolerance = 1e-12)
})

## O13. A launch does not move the overlap, and this is why start lives on
## n_panel() rather than here

# Explicit cohort bookkeeping, deliberately not calling design_overlap(): a
# cohort is an entry period and a per-stage weight vector, and the overlap
# between two occasions is counted from the cohorts alive at each.
.cohort_overlap <- function(cohorts, t, m) {
  at <- function(p) {
    vapply(cohorts, function(co) {
      s <- p - co$entry + 1L
      if (s >= 1L && s <= length(co$w)) co$w[s] else 0
    }, numeric(1L))
  }
  a <- at(t)
  b <- at(t + m)
  sum(pmin(a, b)) / sum(a)
}

.launch <- function(w, start, n_period = 20L) {
  life <- length(w)
  entering <- lapply(2:n_period, function(p) list(entry = p, w = w))
  first <- switch(
    start,
    # one cohort, the design fills up over its life
    gradual = list(list(entry = 1L, w = w)),
    # one launch component for every remaining life length; for a gapped life,
    # some components are selected before their first in-sample interview
    immediate = lapply(seq_len(life), function(k) list(entry = 1L, w = w[k:life])),
    # the naive reading of an immediate start: split by remaining length
    truncated = lapply(seq_len(life), function(k) list(entry = 1L, w = w[seq_len(k)]))
  )
  c(first, entering)
}

test_that("an immediate start reproduces the steady-state overlap from occasion 1", {
  for (spec in c("6", "5", "1-1-0-0-1-1", "4-8-4")) {
    x <- design_overlap(spec)
    w <- as.double(x$rotation)
    imm <- .launch(w, "immediate", n_period = 3L * length(w))
    for (t in 1:4) {
      expect_equal(
        vapply(seq_along(x), function(m) .cohort_overlap(imm, t, m), numeric(1L)),
        as.double(x),
        tolerance = 1e-12,
        info = sprintf("%s at t = %d", spec, t)
      )
    }
  }
})

test_that("a gradual start reaches the steady-state overlap at the life", {
  x <- design_overlap("6")
  grad <- .launch(as.double(x$rotation), "gradual")
  # nothing has rotated out yet, so early overlaps are higher than the design's
  expect_equal(vapply(1:5, function(m) .cohort_overlap(grad, 1L, m), numeric(1L)),
               rep(1, 5L), tolerance = 1e-12)
  expect_gt(.cohort_overlap(grad, 3L, 1L), x[1])
  # and from the occasion that spans the life they are the design's
  for (t in 6:9) {
    expect_equal(vapply(1:5, function(m) .cohort_overlap(grad, t, m), numeric(1L)),
                 as.double(x), tolerance = 1e-12)
  }
})

test_that("splitting the first occasion by remaining length is not an immediate start", {
  # The trap: under a gapped life, covering every stage requires cohorts that
  # start out of sample. Truncating instead holds six cohorts where the design
  # holds four, and reads 0.83 consecutive against the rotation's 0.50.
  x <- design_overlap("1-1-0-0-1-1")
  w <- as.double(x$rotation)
  naive <- .launch(w, "truncated")
  in_sample <- function(co, t) {
    sum(vapply(co, function(c1) {
      s <- t - c1$entry + 1L
      if (s >= 1L && s <= length(c1$w)) c1$w[s] else 0
    }, numeric(1L)))
  }
  expect_equal(in_sample(naive, 1L), 6)
  expect_equal(in_sample(.launch(w, "immediate"), 1L), x$n_occasion)
  expect_equal(.cohort_overlap(naive, 1L, 1L), 5 / 6, tolerance = 1e-12)
  expect_equal(x[1], 0.5, tolerance = 1e-12)
  # the totals do not even hold constant under the naive launch
  expect_false(isTRUE(all.equal(
    vapply(1:6, function(t) in_sample(naive, t), numeric(1L)),
    rep(x$n_occasion, 6L)
  )))
})

test_that("design_overlap takes no launch argument, and does not need one", {
  # Launch changes the route to the mature life-stage mix, not the overlap that
  # mix produces. Both launch policies therefore reach the same lag profile,
  # so a start argument here would be a second name for one answer.
  expect_error(design_overlap("6", start = "immediate"), "unused argument")
})
