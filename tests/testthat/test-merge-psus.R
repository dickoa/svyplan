## merge_psus(): PSUs below a minimum size merged with their neighbours in
## row order

merged_sizes <- function(id, merged, size = rep(1, length(id))) {
  as.vector(tapply(size, factor(merged, levels = unique(merged)), sum))
}

test_that("small PSUs merge in order, and large ones stay as they are", {
  # Sizes 1, 2, 1, 4, 1 at a floor of 3.
  id <- c("N1", "N2", "N2", "N3", "N4", "N4", "N4", "N4", "N5")
  expect_identical(
    merge_psus(id, 3),
    c("N1", "N1", "N1", "N1", "N4", "N4", "N4", "N4", "N4")
  )
  # A run of small PSUs closes a merged PSU as soon as it reaches the floor.
  sizes <- c(5, 2, 2, 2, 2, 1)
  id <- rep(sprintf("P%d", seq_along(sizes)), sizes)
  merged <- merge_psus(id, 3)
  expect_identical(unique(merged), c("P1", "P2", "P4"))
  expect_identical(merged_sizes(id, merged), c(5, 4, 5))
})

test_that("PSUs are taken in order of first appearance, not of their ids", {
  # Z (1 unit), A (2), M (1): Z opens a short run and joins A, M joins them.
  expect_identical(merge_psus(c("Z", "A", "A", "M"), 2), rep("Z", 4))
})

test_that("a short run joins the PSU before it, or after it when it comes first", {
  expect_identical(.merge_psu_groups(c(5, 1, 5), 3), c(1L, 1L, 3L))
  expect_identical(.merge_psu_groups(c(1, 1, 5, 1, 1, 1, 5), 3),
                   c(1L, 1L, 1L, 4L, 4L, 4L, 7L))
  expect_identical(.merge_psu_groups(c(4, 4), 3), c(1L, 2L))
  # A stratum below the floor in total is one PSU.
  expect_identical(.merge_psu_groups(c(1, 1, 1), 5), c(1L, 1L, 1L))
  # With no PSU at the floor, the whole stratum is one run.
  expect_identical(.merge_psu_groups(c(1, 1, 1, 1, 1), 2),
                   c(1L, 1L, 3L, 3L, 3L))
})

test_that("every merged PSU reaches the floor and is contiguous in order", {
  set.seed(3)
  for (rep in 1:200) {
    k <- sample(1:15, 1)
    psu_size <- sample(1:12, k, replace = TRUE)
    floor <- sample(2:20, 1)
    group <- .merge_psu_groups(psu_size, floor)
    merged <- as.vector(
      tapply(psu_size, factor(group, levels = unique(group)), sum)
    )
    if (sum(psu_size) >= floor) {
      expect_true(all(merged >= floor))
    } else {
      expect_identical(group, rep(1L, k))
    }
    # No merge skips over a PSU: groups are runs, labelled by their first.
    expect_identical(group, cummax(group))
    expect_true(all(group <= seq_len(k)))
    expect_identical(unique(group), which(group == seq_len(k)))
    # Two PSUs at or above the floor never share a group.
    big <- psu_size >= floor
    expect_false(anyDuplicated(group[big]) > 0)
    # Small PSUs close a merged PSU once it reaches the floor, so one made
    # of small PSUs alone stays below three floors.
    small_only <- tapply(big, group, function(b) !any(b))
    if (sum(psu_size) >= floor) {
      expect_true(all(merged[small_only] < 3 * floor))
    }
  }
})

test_that("PSUs never merge across strata when applied per stratum", {
  frame <- data.frame(
    stratum = rep(c("a", "b"), c(4, 4)),
    ea = c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")
  )
  psu <- unsplit(
    lapply(split(frame$ea, frame$stratum), merge_psus, min_size = 2),
    frame$stratum
  )
  expect_identical(psu, c("A1", "A1", "A3", "A3", "B1", "B1", "B3", "B3"))
  # Applied to the whole frame, the merge crosses the stratum boundary.
  expect_identical(merge_psus(c("A1", "B1", "B2"), 2), rep("A1", 3))
})

test_that("size= weighs the rows of a PSU-level frame", {
  psu <- c("E1", "E2", "E3", "E4")
  households <- c(30, 10, 16, 40)
  expect_identical(merge_psus(psu, 25, size = households),
                   c("E1", "E2", "E2", "E4"))
  expect_identical(merge_psus(psu, 25), rep("E1", 4))
})

test_that("the merged id keeps the type and class of the input", {
  f <- factor(c("a", "b", "b", "c"), levels = c("c", "b", "a"))
  out <- merge_psus(f, 2)
  expect_s3_class(out, "factor")
  expect_identical(levels(out), levels(f))
  expect_identical(as.character(out), c("a", "a", "a", "a"))
  expect_identical(merge_psus(c(10L, 11L, 11L, 12L, 12L), 2),
                   c(10L, 10L, 10L, 12L, 12L))
  expect_identical(merge_psus(c(1.5, 2.5, 2.5), 2), c(1.5, 1.5, 1.5))
  expect_identical(merge_psus(character(0), 2), character(0))
})

test_that("merge_psus() refuses inputs it cannot read", {
  expect_error(merge_psus(c("a", NA), 2), "missing")
  expect_error(merge_psus(list("a", "b"), 2), "atomic")
  for (bad in list(0, c(2, 3), NA_real_, "2", Inf)) {
    expect_error(merge_psus(c("a", "b"), bad), "'min_size'")
  }
  for (bad in list(1, c(1, -1), c(1, NA), c("1", "2"))) {
    expect_error(merge_psus(c("a", "b"), 2, size = bad), "'size'")
  }
})

test_that("a register refused for small PSUs plans after merging", {
  # Ten PSUs per stratum, three of them smaller than the take of 8.
  sizes <- c(60, 5, 4, 70, 55, 3, 65, 80, 50, 45)
  register <- data.frame(
    psu_id = c(sprintf("A%02d", 1:10), sprintf("B%02d", 1:10)),
    stratum = rep(c("A", "B"), each = 10),
    N = c(sizes, rev(sizes))
  )
  frame <- data.frame(
    stratum = rep(register$stratum, register$N),
    psu_id = rep(register$psu_id, register$N)
  )
  plan_for <- function(psu) {
    n_alloc(
      data.frame(stratum = c("A", "B"),
                 N = as.numeric(tapply(psu$N, psu$stratum, sum)),
                 n_per_psu = 8),
      measures = data.frame(stratum = c("A", "B"), name = "y", p = 0.5,
                            icc_psu = 0.05),
      targets = data.frame(name = "y", cv = 0.15),
      psu = psu
    )
  }
  expect_error(plan_for(register), "fewer units than the take")

  frame$merged <- unsplit(
    lapply(split(frame$psu_id, frame$stratum), merge_psus, min_size = 8),
    frame$stratum
  )
  merged <- aggregate(N ~ stratum + psu_id,
                      transform(frame, psu_id = merged, N = 1), sum)
  expect_true(all(merged$N >= 8))
  expect_identical(sum(merged$N), sum(register$N))
  plan <- plan_for(merged)
  expect_true(all(plan$psu$n_take <= plan$psu$N))
})

test_that("a register naming a PSU twice is refused", {
  psu <- data.frame(
    psu_id = c("P1", "P2", "P2", "P3"),
    stratum = c("A", "A", "B", "B"),
    N = c(40, 30, 20, 50)
  )
  frame <- data.frame(stratum = c("A", "B"), N = c(70, 70), n_per_psu = 10)
  measures <- data.frame(stratum = c("A", "B"), name = "y", p = 0.5,
                         icc_psu = 0.05)
  targets <- data.frame(name = "y", cv = 0.1)
  expect_error(
    n_alloc(frame, measures = measures, targets = targets, psu = psu),
    "'psu\\$psu_id' must name each PSU once, and .P2. appears more than once"
  )
  expect_error(
    prec_alloc(frame, n = c(30, 30), measures = measures, targets = targets,
               psu = psu),
    "must name each PSU once"
  )
})

test_that("a zero-size PSU joins a group like any small one", {
  expect_identical(merge_psus(c("a", "b"), 10, size = c(10, 0)), c("a", "a"))
  expect_identical(merge_psus(c("a", "b"), 10, size = c(0, 10)), c("a", "a"))
  # A run closes on its third PSU and leaves a zero-size one after it.
  expect_identical(
    merge_psus(c("a", "b", "c", "d"), 10, size = c(4, 3, 3, 0)),
    rep("a", 4)
  )
  expect_identical(
    merge_psus(c("a", "b", "c", "d"), 10, size = c(20, 4, 6, 0)),
    c("a", "b", "b", "b")
  )
})
