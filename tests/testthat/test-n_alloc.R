test_that("n_alloc: neyman allocation correctness", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  expect_s3_class(res, "svyplan_n")
  expect_equal(res$type, "alloc")
  expect_equal(res$method, "neyman")
  expect_true(!is.null(res$detail))
  expect_true(abs(sum(res$detail$n) - 600) < 1)

  a_h <- frame$N * frame$sd
  expected_frac <- a_h / sum(a_h)
  actual_frac <- res$detail$n / sum(res$detail$n)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: proportional allocation", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600, alloc = "proportional")
  expected_frac <- frame$N / sum(frame$N)
  actual_frac <- res$detail$n / sum(res$detail$n)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: optimal allocation with cost", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    cost = c(1, 4, 1)
  )
  res <- n_alloc(frame, n = 600, alloc = "optimal")
  a_h <- frame$N * frame$sd / sqrt(c(1, 4, 1))
  expected_frac <- a_h / sum(a_h)
  actual_frac <- res$detail$n / sum(res$detail$n)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: power allocation", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600, alloc = "power", power_q = 0.5)
  a_h <- frame$sd * frame$N^0.5
  expected_frac <- a_h / sum(a_h)
  actual_frac <- res$detail$n / sum(res$detail$n)
  expect_equal(actual_frac, expected_frac, tolerance = 0.02)
})

test_that("n_alloc: ORIC rounding preserves sum", {
  frame <- data.frame(
    N = c(100, 200, 300, 400),
    sd = c(5, 10, 15, 8)
  )
  res <- n_alloc(frame, n = 50)
  n_int <- res$detail$n_int
  expect_equal(sum(n_int), 50L)
  expect_true(all(n_int == floor(n_int)))
})

test_that("n_alloc: min_n floor respected", {
  frame <- data.frame(
    N = c(10000, 10, 10),
    sd = c(100, 1, 1)
  )
  res <- n_alloc(frame, n = 20, min_n = 5)
  expect_true(all(res$detail$n >= 5 - 1e-6))
})

test_that("n_alloc: max_weight column constraint", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    max_weight = c(10, 10, 10)
  )
  res <- n_alloc(frame, n = 600)
  weights <- frame$N / res$detail$n
  expect_true(all(weights <= 10 + 0.5))
})

test_that("n_alloc: take_all column", {
  frame <- data.frame(
    N = c(1000, 2000, 50),
    sd = c(10, 20, 15),
    take_all = c(FALSE, FALSE, TRUE)
  )
  res <- n_alloc(frame, n = 400)
  expect_equal(res$detail$n[3], 50, tolerance = 1e-6)
})

test_that("n_alloc: cv target without domains", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, cv = 0.05)
  expect_true(res$cv <= 0.05 + 1e-3)
})

test_that("n_alloc: cv target requires mean", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  expect_error(n_alloc(frame, cv = 0.05), "mean")
})

test_that("n_alloc: budget mode", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    cost = c(1, 2, 1.5)
  )
  res <- n_alloc(frame, budget = 500)
  expect_s3_class(res, "svyplan_n")
  total_cost <- sum(res$detail$n * frame$cost)
  expect_true(abs(total_cost - 500) < 1)
})

test_that("n_alloc: explicit domains parameter", {
  frame <- data.frame(
    region = c("N", "N", "S", "S"),
    N = c(1000, 2000, 1500, 500),
    sd = c(10, 20, 15, 8),
    mean = c(50, 70, 60, 45)
  )
  res <- n_alloc(frame, domains = "region", cv = 0.05)
  expect_false(is.null(res$domains))
  expect_equal(nrow(res$domains), 2L)
  expect_true(all(res$domains$.cv <= 0.05 + 1e-3))
})

test_that("n_alloc: domain n target", {
  frame <- data.frame(
    region = c("N", "N", "S", "S"),
    N = c(1000, 2000, 1500, 500),
    sd = c(10, 20, 15, 8),
    mean = c(50, 70, 60, 45)
  )
  res <- n_alloc(frame, domains = "region", n = 500)
  expect_false(is.null(res$domains))
})

test_that("n_alloc: detail columns", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  expect_true(all(c("stratum", "N", "sd", "n", "n_int", "weight")
                  %in% names(res$detail)))
})

test_that("n_alloc: se/moe/cv stored on object", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  expect_true(!is.na(res$se))
  expect_true(!is.na(res$moe))
  expect_true(!is.na(res$cv))
})

test_that("n_alloc: print output", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  out <- capture.output(print(res))
  expect_true(any(grepl("Stratum allocation", out)))
  expect_true(any(grepl("neyman", out)))
})

test_that("n_alloc: format output", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  fmt <- format(res)
  expect_true(grepl("alloc", fmt))
})

test_that("n_alloc: as.integer returns total n", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15)
  )
  res <- n_alloc(frame, n = 600)
  expect_equal(as.integer(res), as.integer(ceiling(res$n)))
})

test_that("n_alloc: validation errors", {
  expect_error(n_alloc(data.frame(N = 1:3), n = 10), "sd.*var")
  expect_error(
    n_alloc(data.frame(N = -1, sd = 1), n = 10),
    "positive"
  )
  expect_error(
    n_alloc(data.frame(N = 100, sd = -1), n = 10),
    "non-negative"
  )
  expect_error(
    n_alloc(data.frame(N = c(10, 20), sd = c(1, 2)),
            n = 5, cv = 0.1),
    "exactly one"
  )
  expect_error(
    n_alloc(data.frame(N = c(10, 20), sd = c(1, 2))),
    "exactly one"
  )
  expect_error(
    n_alloc(data.frame(N = c(10, 20), sd = c(1, 2)), n = 100),
    "maximum feasible"
  )
})

test_that("n_alloc: rejects conflicting frame columns", {
  f_both_sv <- data.frame(N = c(100, 200), sd = c(1, 2), var = c(1, 4))
  expect_error(n_alloc(f_both_sv, n = 50), "not both")
  f_both_mp <- data.frame(N = c(100, 200), sd = c(1, 2),
                          mean = c(10, 20), p = c(0.1, 0.2))
  expect_error(n_alloc(f_both_mp, cv = 0.05), "not both")
})

test_that("n_alloc: p must be in [0, 1]", {
  expect_error(
    n_alloc(data.frame(N = c(100, 200), sd = c(1, 2), p = c(0.3, 1.5)),
            cv = 0.05),
    "\\[0, 1\\]"
  )
  expect_error(
    n_alloc(data.frame(N = c(100, 200), sd = c(1, 2), p = c(-0.1, 0.5)),
            cv = 0.05),
    "\\[0, 1\\]"
  )
})

test_that("n_alloc: duplicate stratum labels rejected", {
  expect_error(
    n_alloc(
      data.frame(stratum = c("A", "A", "B"), N = c(100, 200, 300), sd = c(1, 2, 3)),
      n = 50
    ),
    "duplicate stratum"
  )
})

test_that("n_alloc: duplicate stratum within domain rejected", {
  expect_error(
    n_alloc(
      data.frame(
        region = c("N", "N", "S"),
        stratum = c("U", "U", "U"),
        N = c(100, 200, 300),
        sd = c(1, 2, 3),
        mean = c(50, 60, 55)
      ),
      domains = "region",
      cv = 0.05
    ),
    "duplicate stratum.*domain"
  )
})

test_that("n_alloc: same stratum label across domains is allowed", {
  frame <- data.frame(
    region = c("N", "S"),
    stratum = c("Urban", "Urban"),
    N = c(1000, 2000),
    sd = c(10, 20),
    mean = c(50, 60)
  )
  res <- n_alloc(frame, domains = "region", cv = 0.05)
  expect_s3_class(res, "svyplan_n")
})

test_that("n_alloc: max_weight < 1 rejected", {
  expect_error(
    n_alloc(
      data.frame(N = c(100, 200), sd = c(1, 2), max_weight = c(0.5, 10)),
      n = 50
    ),
    ">= 1"
  )
})

test_that("n_alloc: all-zero sd warns", {
  expect_warning(
    n_alloc(data.frame(N = c(100, 200), sd = c(0, 0)), n = 50),
    "no variability"
  )
})

test_that("n_alloc: all-zero mean warns about Inf CV", {
  expect_warning(
    n_alloc(
      data.frame(N = c(100, 200), sd = c(1, 2), mean = c(0, 0)),
      n = 50
    ),
    "CV will be Inf"
  )
})

test_that("n_alloc: single stratum", {
  frame <- data.frame(N = 1000, sd = 10)
  res <- n_alloc(frame, n = 100)
  expect_true(abs(res$detail$n - 100) < 1)
})

test_that("n_alloc: all take_all", {
  frame <- data.frame(
    N = c(50, 30),
    sd = c(10, 5),
    take_all = c(TRUE, TRUE)
  )
  res <- n_alloc(frame, n = 80)
  expect_equal(res$detail$n[1], 50, tolerance = 1e-6)
  expect_equal(res$detail$n[2], 30, tolerance = 1e-6)
})

test_that("n_alloc: infeasible n (too small)", {
  frame <- data.frame(
    N = c(100, 200, 300),
    sd = c(10, 20, 15)
  )
  expect_error(n_alloc(frame, n = 0.5), "minimum feasible")
})

test_that("n_alloc: domain print output", {
  frame <- data.frame(
    region = c("N", "N", "S", "S"),
    N = c(1000, 2000, 1500, 500),
    sd = c(10, 20, 15, 8),
    mean = c(50, 70, 60, 45)
  )
  res <- n_alloc(frame, domains = "region", cv = 0.05)
  out <- capture.output(print(res))
  expect_true(any(grepl("Domains", out)))
})

test_that("n_alloc: resp_rate and deff params stored", {
  frame <- data.frame(
    N = c(1000, 2000),
    sd = c(10, 20)
  )
  res <- n_alloc(frame, n = 300, resp_rate = 0.8, deff = 1.5)
  expect_equal(res$params$resp_rate, 0.8)
  expect_equal(res$params$deff, 1.5)
})

test_that("n_alloc: cost column in data", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    cost = c(1, 4, 1)
  )
  res <- n_alloc(frame, n = 600, alloc = "optimal")
  expect_equal(res$method, "optimal")
})

test_that("n_alloc: cv domain convergence", {
  frame <- data.frame(
    region = c("A", "A", "B", "B", "B"),
    N = c(1000, 500, 2000, 1500, 800),
    sd = c(10, 20, 15, 8, 12),
    mean = c(50, 70, 60, 45, 55)
  )
  res <- n_alloc(frame, domains = "region", cv = 0.04)
  expect_true(all(res$domains$.cv <= 0.04 + 1e-3))
})

test_that("n_alloc: var column alternative to sd", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    var = c(100, 400, 225),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  expect_s3_class(res, "svyplan_n")
})

test_that("n_alloc: p column alternative to mean", {
  frame <- data.frame(
    N = c(1000, 2000),
    sd = c(0.49, 0.48),
    p = c(0.3, 0.5)
  )
  res <- n_alloc(frame, cv = 0.10)
  expect_true(!is.na(res$cv))
})

test_that("prec_alloc: basic usage", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- prec_alloc(frame, n = c(100, 200, 300))
  expect_s3_class(res, "svyplan_prec")
  expect_equal(res$type, "alloc")
  expect_true(!is.na(res$se))
  expect_true(!is.na(res$cv))
})

test_that("prec_alloc: round-trip from n_alloc", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res_n <- n_alloc(frame, n = 600)
  res_p <- prec_alloc(res_n)
  expect_s3_class(res_p, "svyplan_prec")
  expect_equal(res_p$cv, res_n$cv, tolerance = 1e-4)
})

test_that("prec_alloc: validation", {
  frame <- data.frame(
    N = c(1000, 2000),
    sd = c(10, 20)
  )
  expect_error(prec_alloc(frame, n = c(100)), "length")
  expect_error(prec_alloc(frame, n = c(-1, 100)), "positive")
})

test_that("n_alloc: predict grid", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(deff = c(1, 1.5, 2)))
  expect_equal(nrow(grid), 3L)
  expect_true("n" %in% names(grid))
})

test_that("n_alloc: predict with resp_rate", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(resp_rate = c(0.8, 0.9, 1.0)))
  expect_equal(nrow(grid), 3L)
})

test_that("n_alloc: predict returns cost column", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60),
    cost = c(1, 2, 1.5)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(deff = c(1, 1.5)))
  expect_true("cost" %in% names(grid))
})

test_that("n_alloc: predict rejects multiple mode switches", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  expect_error(
    predict(res, data.frame(n = 500, cv = 0.05)),
    "at most one"
  )
})

test_that("n_alloc: predict switches from n to cv mode", {
  frame <- data.frame(
    N = c(1000, 2000, 3000),
    sd = c(10, 20, 15),
    mean = c(50, 70, 60)
  )
  res <- n_alloc(frame, n = 600)
  grid <- predict(res, data.frame(cv = c(0.05, 0.10)))
  expect_equal(nrow(grid), 2L)
  expect_true(grid$n[1] > grid$n[2])
})

test_that("n_alloc: domains must be character", {
  frame <- data.frame(N = 100, sd = 10, region = "A")
  expect_error(n_alloc(frame, domains = 1, n = 50), "character")
})

test_that("n_alloc: domains columns must exist in frame", {
  frame <- data.frame(N = 100, sd = 10)
  expect_error(n_alloc(frame, domains = "region", n = 50), "not found")
})

test_that("n_alloc: extra columns ignored when domains is NULL", {
  frame <- data.frame(
    region = c("N", "S"),
    N = c(1000, 2000),
    sd = c(10, 15),
    mean = c(50, 60)
  )
  res <- n_alloc(frame, n = 100)
  expect_null(res$domains)
})
