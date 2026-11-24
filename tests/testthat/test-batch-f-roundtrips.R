## Allocation round trips and totals

test_that("a floored cluster design in n mode keeps its total", {
  f3 <- data.frame(stratum = c("A", "B", "C", "D"),
                   N = c(20000, 30000, 12000, 8000), sd = c(10, 14, 8, 9),
                   mean = c(40, 50, 45, 30), icc_1 = c(0.05, 0.1, 0.08, 0.06),
                   icc_2 = c(0.15, 0.2, 0.1, 0.12),
                   cost_1 = c(400, 600, 300, 500), cost_2 = c(50, 40, 60, 45),
                   cost_3 = c(10, 12, 8, 9), reg = c("x", "x", "y", "y"))
  f2 <- f3[, setdiff(names(f3), c("icc_2", "cost_3"))]
  for (frame in list(f2, f3)) {
    for (a in c("optimal", "neyman", "proportional")) {
      fit <- suppressWarnings(n_alloc(frame, n = 600, domains = "reg",
                                      min_n_domain = 260, alloc = a))
      take <- if ("take_3_int" %in% names(fit$detail)) {
        fit$detail$take_2_int * fit$detail$take_3_int
      } else {
        fit$detail$take_2_int
      }
      info <- paste(ncol(frame), a)
      expect_true(abs(sum(fit$detail$n_int) - 600) <= max(take) / 2,
                  info = info)
      expect_true(all(tapply(fit$detail$n_int, frame$reg, sum) >= 260),
                  info = info)
    }
  }
})

test_that("a domain search keeps its control through a precision round trip", {
  frame <- data.frame(stratum = c("a", "b", "c", "d"),
                      reg = c("N", "N", "S", "S"),
                      N = c(2000, 3000, 2500, 1500), sd = c(10, 14, 12, 9),
                      mean = c(40, 50, 45, 30))
  fit <- n_alloc(frame, cv = 0.05, alloc = "optimal", domains = "reg",
                 control = list(max_nodes = 3))
  back <- n_alloc(prec_alloc(fit), cv = 0.05)
  expect_identical(back$params$control, list(max_nodes = 3))
  expect_equal(back$detail$n_int, fit$detail$n_int)
  expect_null(n_alloc(prec_alloc(fit))$params$control)
})
