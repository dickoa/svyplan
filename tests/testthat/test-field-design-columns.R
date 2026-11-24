## Whole field design columns of every allocation
##
## Oracle: the field design is n_1_int PSUs with take_s_int units at each
## later stage, so n_int is their product within a stratum. A register fields
## its certainty PSUs at their own takes and draws n_psu_draw more at
## take_2_int, so n_int is rebuilt from the PSU table.

field_measures <- function() {
  data.frame(stratum = rep(c("U", "R"), 2),
             name = rep(c("p1", "y1"), each = 2),
             p = c(0.35, 0.45, NA, NA), mean = c(NA, NA, 100, 120),
             sd = c(NA, NA, 25, 30),
             icc_1 = c(0.07, 0.04, 0.03, 0.08))
}
field_targets <- data.frame(name = c("p1", "y1"), cv = c(0.08, 0.04))
field_frame <- data.frame(stratum = c("U", "R"), N = c(1600, 2000),
                          N_1 = c(40, 50), cost_1 = c(110, 190),
                          cost_2 = c(14, 10))

expect_whole_design <- function(d, stages, label) {
  take_int <- paste0("take_", seq_len(stages)[-1L], "_int")
  expect_true(all(take_int %in% names(d)), label = label)
  units <- d$n_1_int * Reduce(`*`, d[take_int])
  expect_equal(units, d$n_int, label = label)
}

test_that("every non-register allocation reports n_1_int and take_s_int", {
  classic <- n_alloc(
    data.frame(stratum = c("A", "B"), N = c(4000, 6000), sd = c(10, 12),
               icc_1 = 0.05, take_2 = c(10, 8), cost_1 = 500, cost_2 = 10,
               N_1 = c(80, 120)),
    n = 600
  )
  fixed_2 <- n_alloc_joint(transform(field_frame, take_2 = c(12, 8)),
                           field_measures(), targets = field_targets)
  fixed_3 <- n_alloc_joint(
    transform(field_frame, N = c(16000, 20000), N_2 = c(400, 500),
              cost_2 = c(30, 20), cost_3 = c(5, 4),
              take_2 = c(4, 3), take_3 = c(6, 8)),
    transform(field_measures(), icc_2 = c(0.1, 0.1, 0.05, 0.05)),
    targets = field_targets
  )
  free <- n_alloc_joint(field_frame, field_measures(), targets = field_targets,
                        take_2_min = 6, take_2_max = 18, take_2_step = 4)

  expect_whole_design(classic$detail, 2L, "classic")
  expect_whole_design(fixed_2$detail, 2L, "fixed, 2 stages")
  expect_whole_design(fixed_3$detail, 3L, "fixed, 3 stages")
  expect_whole_design(free$detail, 2L, "free")
  expect_identical(fixed_3$detail$take_2_int, c(4, 3))
  expect_identical(fixed_3$detail$take_3_int, c(6, 8))
  expect_identical(
    names(free$detail)[names(free$detail) %in%
                         c("take_2", "n_1", "n_1_int", "take_2_int")],
    c("take_2", "n_1", "n_1_int", "take_2_int")
  )
})

register_fits <- function() {
  frame <- data.frame(stratum = c("A", "B"), N = c(2000, 1600),
                      cost_1 = 200, cost_2 = 5)
  psu <- data.frame(stratum = rep(frame$stratum, each = 8),
                    N = c(600, rep(200, 7), 450, 220, 165, rep(153, 5)))
  measures <- data.frame(stratum = frame$stratum, name = "y",
                         p = c(0.5, 0.4), icc_1 = 0.05)
  targets <- data.frame(name = "y", cv = 0.1)
  list(
    fixed = n_alloc_joint(transform(frame, take_2 = c(10, 8)), measures,
                          targets = targets, psu = psu),
    free = n_alloc_joint(frame, measures, targets = targets, psu = psu,
                         take_2_min = 2, take_2_max = 20),
    budget = n_alloc_joint(transform(frame, take_2 = c(10, 8)), measures,
                           objective = "y", budget = 3000, psu = psu)
  )
}

test_that("a register field design is rebuilt from its PSU table", {
  fits <- register_fits()
  checked <- character()
  for (mode in names(fits)) {
    fit <- fits[[mode]]
    d <- fit$detail
    expect_true(all(d$n_psu_certain > 0), label = mode)
    expect_equal(d$n_1_int, d$n_psu_certain + d$n_psu_draw, label = mode)
    expect_equal(d$take_2_int, d$take_2, label = mode)
    held <- fit$psu$certainty
    by_stratum <- factor(fit$psu$stratum, levels = d$stratum)
    certain_take <- as.numeric(tapply(fit$psu$take_2 * held, by_stratum, sum))
    expect_equal(certain_take, d$n_certain_int, label = mode)
    expect_equal(certain_take + d$n_psu_draw * d$take_2_int, d$n_int,
                 label = mode)
    rest <- !held
    expect_equal(fit$psu$take_2[rest],
                 d$take_2_int[match(fit$psu$stratum[rest], d$stratum)],
                 label = mode)
    expect_false("n_take" %in% names(fit$psu))
    checked <- c(checked, mode)
  }
  expect_setequal(checked, c("fixed", "free", "budget"))
})

test_that("register grid rows report n_1_int, missing when the row fails", {
  fit <- register_fits()$budget
  grid <- suppressWarnings(predict(fit, data.frame(budget = c(50, 3000))))
  expect_identical(grid$.status, c("infeasible", "optimal"))
  expect_true(is.na(grid$n_1_int[1]))
  expect_identical(grid$n_1_int[2], grid$n_psu_certain[2] + grid$n_psu_draw[2])
  expect_identical(grid$n_1_int[2], sum(fit$detail$n_1_int))
  cols <- names(grid)
  expect_lt(match("n_1_int", cols), match("n_int", cols))
})
