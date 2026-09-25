test_that("joint proportion precision matches exact finite-population sampling", {
  population <- c(rep(1, 6), rep(0, 14))
  samples <- combn(population, 5, FUN = mean)
  exact <- mean((samples - mean(population))^2)
  frame <- data.frame(stratum = "S", N = length(population))
  measures <- data.frame(stratum = "S", name = "y", p = mean(population))
  targets <- data.frame(name = "y", cv = .5)
  got <- prec_alloc(frame, n = 5, measures = measures, targets = targets)
  expect_equal(got$se^2, exact, tolerance = 1e-12)
  expect_equal(got$se, prec_prop(.3, n = 5, N = 20)$se)
  expect_equal(prec_alloc(frame, n = 20, measures = measures,
                          targets = targets)$se, 0)
})

test_that("joint and scalar proportion planning share response and FPC conventions", {
  for (N in c(20, 2000)) for (response in c(.8, 1)) for (deff in c(1, 1.7)) {
    frame <- data.frame(stratum = "S", N = N)
    measures <- data.frame(stratum = "S", name = "y", p = .3)
    targets <- data.frame(name = "y", cv = .5)
    scalar <- prec_prop(.3, n = 5, N = N, resp_rate = response,
                         deff = deff, df = 9, alpha = .1)
    joint <- prec_alloc(frame, n = 5, measures = measures, targets = targets,
                        resp_rate = response, deff = deff, df = 9, alpha = .1)
    expect_equal(joint$se, scalar$se)
    expect_equal(joint$cv, scalar$cv)
    expect_equal(joint$moe, scalar$moe)

    sized <- n_alloc(frame, measures = measures, targets = targets,
                      resp_rate = response, deff = deff)
    expect_equal(sized$detail$n,
                 n_prop(.3, cv = .5, N = N, resp_rate = response, deff = deff)$n,
                 tolerance = 1e-6)
    expect_equal(prec_alloc(sized)$cv, .5, tolerance = 1e-6)
    expect_true(all(sized$operational$constraints$.pass))
  }
})

test_that("stratum-specific proportion corrections match stratified enumeration", {
  a <- c(1, 0, 0, 0)
  b <- c(1, 1, 0, 0, 0, 0)
  estimates <- outer(combn(a, 2, FUN = mean) * .4,
                     combn(b, 3, FUN = mean) * .6, "+")
  exact <- mean((estimates - mean(c(a, b)))^2)
  frame <- data.frame(stratum = c("A", "B"), N = c(4, 6))
  # Deliberately reverse measure rows to exercise matching by stratum.
  measures <- data.frame(stratum = c("B", "A"), name = "y",
                         p = c(mean(b), mean(a)))
  targets <- data.frame(name = "y", cv = .5)
  got <- prec_alloc(frame, n = c(A = 2, B = 3), measures = measures,
                    targets = targets)
  expect_equal(got$se^2, exact, tolerance = 1e-12)
  targets$domain <- "stratum"
  targets$level <- "A"
  domain <- prec_alloc(frame, n = c(A = 2, B = 3), measures = measures,
                       targets = targets)
  expect_equal(domain$se, prec_prop(mean(a), n = 2, N = 4)$se)
})

test_that("budget objectives use the corrected proportion variance", {
  frame <- data.frame(stratum = "S", N = 20)
  measures <- data.frame(stratum = "S", name = "y", p = .3)
  fit <- n_alloc(frame, measures = measures, objective = "y", budget = 5)
  expect_equal(fit$detail$n, 5)
  expect_equal(fit$objective_value, prec_prop(.3, n = 5, N = 20)$cv^2)
})

test_that("element proportion measures reject undefined finite-population variances", {
  frame <- data.frame(stratum = "S", N = 1)
  measures <- data.frame(stratum = "S", name = "y", p = .3)
  expect_error(n_alloc(frame, measures = measures,
                       targets = data.frame(name = "y", cv = .5)),
               "N.*must be greater than 1.*proportion")
})

test_that("PSU-register designs retain their working proportion variance", {
  frame <- data.frame(stratum = "S", N = 1000, n_per_psu = 10)
  psu <- data.frame(stratum = "S", psu_id = 1:10, N = c(550, rep(50, 9)))
  measures <- data.frame(stratum = "S", name = "y", p = .3, icc_psu = .08)
  targets <- data.frame(name = "y", cv = .2)
  explicit <- measures
  explicit$mean <- explicit$p
  explicit$var <- explicit$p * (1 - explicit$p)
  explicit$p <- NULL
  a <- n_alloc(frame, measures = measures, targets = targets, psu = psu)
  b <- n_alloc(frame, measures = explicit, targets = targets, psu = psu)
  expect_equal(a$detail$n, b$detail$n)
  expect_equal(a$detail$n_int, b$detail$n_int)
  expect_equal(prec_alloc(a)$se, prec_alloc(b)$se)
  expect_equal(prec_alloc(frame, n = 80, measures = measures,
                          targets = targets, psu = psu)$se,
               prec_alloc(frame, n = 80, measures = explicit,
                          targets = targets, psu = psu)$se)
})
