test_that("n_alloc: neyman allocation correctness", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$type, "alloc")
  expect_equal(res$method, "neyman")
  expect_true(!is.null(res$detail))
  expect_true(abs(sum(res$detail$n_h) - 600) < 1)

  a_h <- frame$N_h * frame$S_h
  expected_frac <- a_h / sum(a_h)
  actual_frac <- res$detail$n_h / sum(res$detail$n_h)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: proportional allocation", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600, alloc = "proportional")
  expected_frac <- frame$N_h / sum(frame$N_h)
  actual_frac <- res$detail$n_h / sum(res$detail$n_h)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: optimal allocation with cost", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    cost_h = c(1, 4, 1)
  )
  res <- n_alloc(frame, n = 600, alloc = "optimal")
  a_h <- frame$N_h * frame$S_h / sqrt(c(1, 4, 1))
  expected_frac <- a_h / sum(a_h)
  actual_frac <- res$detail$n_h / sum(res$detail$n_h)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: power allocation", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600, alloc = "power", power_q = 0.5)
  a_h <- frame$S_h * frame$N_h^0.5
  expected_frac <- a_h / sum(a_h)
  actual_frac <- res$detail$n_h / sum(res$detail$n_h)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: ORIC rounding preserves sum", {
  frame <- data.frame(
    N_h = c(100, 200, 300, 400),
    S_h = c(5, 10, 15, 8)
  )
  res <- n_alloc(frame, n = 50)
  n_h_int <- res$detail$n_h_int
  expect_equal(sum(n_h_int), 50L)
  expect_true(all(n_h_int == floor(n_h_int)))
})

test_that("n_alloc: min_n floor respected", {
  frame <- data.frame(
    N_h = c(10000, 10, 10),
    S_h = c(100, 1, 1)
  )
  res <- n_alloc(frame, n = 20, min_n = 5)
  expect_true(all(res$detail$n_h >= 5 - 1e-6))
})

test_that("n_alloc: max_weight column constraint", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    max_weight = c(10, 10, 10)
  )
  res <- n_alloc(frame, n = 600)
  weights <- frame$N_h / res$detail$n_h
  expect_true(all(weights <= 10 + 0.5))
})

test_that("n_alloc: take_all column", {
  frame <- data.frame(
    N_h = c(1000, 2000, 50),
    S_h = c(10, 20, 15),
    take_all = c(FALSE, FALSE, TRUE)
  )
  res <- n_alloc(frame, n = 400)
  expect_equal(res$detail$n_h[3], 50, tolerance = 1e-6)
})

test_that("n_alloc: cv target without domains", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, cv = 0.05)
  expect_true(res$cv <= 0.05 + 1e-3)
})

test_that("n_alloc: cv target requires mean_h", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  expect_error(n_alloc(frame, cv = 0.05), "mean_h")
})

test_that("n_alloc: budget mode", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    cost_h = c(1, 2, 1.5)
  )
  res <- n_alloc(frame, budget = 500)
  expect_s3_class(res, "svyplan_n")
  total_cost <- sum(res$detail$n_h * frame$cost_h)
  expect_true(abs(total_cost - 500) < 1)
})

test_that("n_alloc: domain auto-detection", {
  frame <- data.frame(
    region = c("N", "N", "S", "S"),
    N_h = c(1000, 2000, 1500, 500),
    S_h = c(10, 20, 15, 8),
    mean_h = c(50, 70, 60, 45)
  )
  res <- n_alloc(frame, cv = 0.05)
  expect_false(is.null(res$domains))
  expect_equal(nrow(res$domains), 2L)
  expect_true(all(res$domains$.cv <= 0.05 + 1e-3))
})

test_that("n_alloc: domain n target", {
  frame <- data.frame(
    region = c("N", "N", "S", "S"),
    N_h = c(1000, 2000, 1500, 500),
    S_h = c(10, 20, 15, 8),
    mean_h = c(50, 70, 60, 45)
  )
  res <- n_alloc(frame, n = 500)
  expect_false(is.null(res$domains))
})

test_that("n_alloc: detail columns", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  expect_true(all(c("stratum", "N_h", "S_h", "n_h", "n_h_int", "weight")
                  %in% names(res$detail)))
})

test_that("n_alloc: se/moe/cv stored on object", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  expect_true(!is.na(res$se))
  expect_true(!is.na(res$moe))
  expect_true(!is.na(res$cv))
})

test_that("n_alloc: print output", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  out <- capture.output(print(res))
  expect_true(any(grepl("Stratum allocation", out)))
  expect_true(any(grepl("neyman", out)))
})

test_that("n_alloc: format output", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  fmt <- format(res)
  expect_true(grepl("alloc", fmt))
})

test_that("n_alloc: as.integer returns total n", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
})

test_that("n_alloc: validation errors", {
  expect_error(n_alloc(data.frame(N_h = 1:3), n = 10), "S_h.*var")
  expect_error(
    n_alloc(data.frame(N_h = -1, S_h = 1), n = 10),
    "positive"
  )
  expect_error(
    n_alloc(data.frame(N_h = 100, S_h = -1), n = 10),
    "non-negative"
  )
  expect_error(
    n_alloc(data.frame(N_h = c(10, 20), S_h = c(1, 2)),
            n = 5, cv = 0.1),
    "exactly one"
  )
  expect_error(
    n_alloc(data.frame(N_h = c(10, 20), S_h = c(1, 2))),
    "exactly one"
  )
  expect_error(
    n_alloc(data.frame(N_h = c(10, 20), S_h = c(1, 2)), n = 100),
    "maximum feasible"
  )
})

test_that("n_alloc: single stratum", {
  frame <- data.frame(N_h = 1000, S_h = 10)
  res <- n_alloc(frame, n = 100)
  expect_true(abs(res$detail$n_h - 100) < 1)
})

test_that("n_alloc: all take_all", {
  frame <- data.frame(
    N_h = c(50, 30),
    S_h = c(10, 5),
    take_all = c(TRUE, TRUE)
  )
  res <- n_alloc(frame, n = 80)
  expect_equal(res$detail$n_h[1], 50, tolerance = 1e-6)
  expect_equal(res$detail$n_h[2], 30, tolerance = 1e-6)
})

test_that("n_alloc: infeasible n (too small)", {
  frame <- data.frame(
    N_h = c(100, 200, 300),
    S_h = c(10, 20, 15)
  )
  expect_error(n_alloc(frame, n = 0.5), "minimum feasible")
})

test_that("n_alloc: domain print output", {
  frame <- data.frame(
    region = c("N", "N", "S", "S"),
    N_h = c(1000, 2000, 1500, 500),
    S_h = c(10, 20, 15, 8),
    mean_h = c(50, 70, 60, 45)
  )
  res <- n_alloc(frame, cv = 0.05)
  out <- capture.output(print(res))
  expect_true(any(grepl("Domains", out)))
})

test_that("n_alloc: resp_rate and deff params stored", {
  frame <- data.frame(
    N_h = c(1000, 2000),
    S_h = c(10, 20)
  )
  res <- n_alloc(frame, n = 300, resp_rate = 0.8, deff = 1.5)
  expect_equal(res$params$resp_rate, 0.8)
  expect_equal(res$params$deff, 1.5)
})

test_that("n_alloc: cost column in data", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    cost_h = c(1, 4, 1)
  )
  res <- n_alloc(frame, n = 600, alloc = "optimal")
  expect_equal(res$method, "optimal")
})

test_that("n_alloc: cv domain convergence", {
  frame <- data.frame(
    region = c("A", "A", "B", "B", "B"),
    N_h = c(1000, 500, 2000, 1500, 800),
    S_h = c(10, 20, 15, 8, 12),
    mean_h = c(50, 70, 60, 45, 55)
  )
  res <- n_alloc(frame, cv = 0.04)
  expect_true(all(res$domains$.cv <= 0.04 + 1e-3))
})

test_that("n_alloc: var column alternative to S_h", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    var = c(100, 400, 225),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  expect_s3_class(res, "svyplan_n")
})

test_that("n_alloc: p_h column alternative to mean_h", {
  frame <- data.frame(
    N_h = c(1000, 2000),
    S_h = c(0.49, 0.48),
    p_h = c(0.3, 0.5)
  )
  res <- n_alloc(frame, cv = 0.10)
  expect_true(!is.na(res$cv))
})

test_that("prec_alloc: basic usage", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- prec_alloc(frame, n = c(100, 200, 300))
  expect_s3_class(res, "svyplan_prec")
  expect_equal(res$type, "alloc")
  expect_true(!is.na(res$se))
  expect_true(!is.na(res$cv))
})

test_that("prec_alloc: round-trip from n_alloc", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res_n <- n_alloc(frame, n = 600)
  res_p <- prec_alloc(res_n)
  expect_s3_class(res_p, "svyplan_prec")
  expect_equal(res_p$cv, res_n$cv, tolerance = 1e-4)
})

test_that("prec_alloc: validation", {
  frame <- data.frame(
    N_h = c(1000, 2000),
    S_h = c(10, 20)
  )
  expect_error(prec_alloc(frame, n = c(100)), "length")
  expect_error(prec_alloc(frame, n = c(-1, 100)), "positive")
})

test_that("n_alloc: predict grid", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(deff = c(1, 1.5, 2)))
  expect_equal(nrow(grid), 3L)
  expect_true("n" %in% names(grid))
})

test_that("n_alloc: predict with resp_rate", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(resp_rate = c(0.8, 0.9, 1.0)))
  expect_equal(nrow(grid), 3L)
})

test_that("n_alloc: predict returns cost column", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60),
    cost_h = c(1, 2, 1.5)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(deff = c(1, 1.5)))
  expect_true("cost" %in% names(grid))
})

test_that("n_alloc: predict rejects multiple mode switches", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  expect_error(
    predict(res, data.frame(n = 500, cv = 0.05)),
    "at most one"
  )
})

test_that("n_alloc: predict switches from n to cv mode", {
  frame <- data.frame(
    N_h = c(1000, 2000, 3000),
    S_h = c(10, 20, 15),
    mean_h = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(cv = c(0.05, 0.10)))
  expect_equal(nrow(grid), 2L)
  expect_true(grid$n[1] > grid$n[2])
})
