## Floors on domain totals in the classic modes

floor_frame <- function() {
  data.frame(
    stratum = paste(rep(c("A", "B", "C"), each = 2), c("u", "r")),
    province = rep(c("A", "B", "C"), each = 2),
    N = c(20000, 30000, 8000, 12000, 3000, 4000),
    p = c(0.2, 0.35, 0.25, 0.4, 0.3, 0.45)
  )
}

domain_totals <- function(fit, column = "n") {
  as.numeric(tapply(fit$detail[[column]], fit$detail$province, sum))
}

## D1. The floored allocation

test_that("held domains keep the rule inside, free domains a common scale", {
  a_h <- c(4, 6, 1, 2, 0.5, 0.8)
  dom <- c(1, 1, 2, 2, 3, 3)
  m_h <- rep(1, 6)
  M_h <- rep(1e4, 6)
  floor_d <- c(0, 300, 300)
  n_h <- svyplan:::.rna_alloc_floored(a_h, 1500, m_h, M_h, dom, floor_d)
  expect_equal(sum(n_h), 1500)
  totals <- tapply(n_h, dom, sum)
  expect_equal(as.numeric(totals[2:3]), c(300, 300))
  ratio <- n_h / a_h
  expect_equal(ratio[1], ratio[2])
  expect_equal(ratio[3], ratio[4])
  expect_equal(ratio[5], ratio[6])
  # A held domain gets more per unit of a_h than the free scale gives.
  expect_gt(ratio[3], ratio[1])
  expect_gt(ratio[5], ratio[1])
  # No feasible move between the free domain's strata and a held domain
  # lowers sum(a_h^2 / n_h).
  obj <- function(n) sum(a_h^2 / n)
  shift <- n_h
  shift[1] <- shift[1] + 1
  shift[2] <- shift[2] - 1
  expect_gte(obj(shift), obj(n_h))
  shift <- n_h
  shift[3] <- shift[3] + 1
  shift[1] <- shift[1] - 1
  expect_gte(obj(shift), obj(n_h))
  # Without floors it is the plain bounded allocation.
  expect_equal(
    svyplan:::.rna_alloc_floored(a_h, 1500, m_h, M_h, dom, rep(0, 3)),
    svyplan:::.rna_alloc(a_h, 1500, m_h, M_h)
  )
})

## D2. Through n_alloc()

test_that("a floor that does not bind changes nothing", {
  fr <- floor_frame()
  base <- n_alloc(fr, n = 1500, domains = "province")
  slack <- n_alloc(fr, n = 1500, domains = "province", min_n_domain = 100)
  expect_equal(slack$detail$n, base$detail$n)
  expect_equal(slack$detail$n_int, base$detail$n_int)
  expect_false(any(slack$domains$.floor_binding))
})

test_that("binding floors hold in the design and in whole units", {
  fr <- floor_frame()
  fit <- n_alloc(fr, n = 1500, domains = "province", min_n_domain = 400)
  expect_equal(domain_totals(fit), c(700, 400, 400))
  expect_true(all(domain_totals(fit, "n_int") >= 400))
  expect_identical(sum(fit$detail$n_int), 1500L)
  expect_identical(fit$domains$.floor_binding, c(FALSE, TRUE, TRUE))
  expect_identical(fit$detail$.bound_source[3:6], rep("min_n_domain", 4))
  expect_equal(fit$params$min_n_domain, 400)

  tab <- data.frame(domain = "province", level = "C", min_n = 500)
  fit <- n_alloc(fr, n = 1500, domains = "province", min_n_domain = tab)
  expect_equal(domain_totals(fit)[3], 500)
  expect_equal(fit$domains$.min_n, c(NA, NA, 500))
})

test_that("floors in cv and budget modes keep the targets", {
  fr <- floor_frame()
  fit <- n_alloc(fr, cv = 0.09, domains = "province", min_n_domain = 400)
  expect_true(all(fit$domains$.cv <= 0.09 + 1e-8))
  expect_true(all(domain_totals(fit, "n_int") >= 400))
  # The floors alone meet the target, so the smallest design is the floors.
  expect_equal(sum(fit$detail$n), 1200)
  # Here they do not, and the worst domain sits on the target.
  fit <- n_alloc(fr, cv = 0.07, domains = "province", min_n_domain = 400)
  expect_equal(max(fit$domains$.cv), 0.07, tolerance = 1e-6)
  expect_true(all(domain_totals(fit) >= 400 - 1e-8))

  fr$unit_cost <- c(10, 20, 10, 20, 15, 30)
  fit <- n_alloc(fr, budget = 30000, domains = "province",
                 min_n_domain = 300, alloc = "optimal")
  expect_true(all(domain_totals(fit, "n_int") >= 300))
  expect_lte(sum(fit$detail$n_int * fr$unit_cost), 30000)
})

test_that("a cluster frame keeps the floors in whole PSUs", {
  fr <- transform(floor_frame(), N_psu = c(200, 300, 80, 120, 30, 40),
                  icc_psu = 0.05, n_per_psu = 12, cost_psu = 300,
                  cost_ssu = 20)
  for (mode in list(list(n = 1500), list(cv = 0.09), list(budget = 60000))) {
    fit <- suppressWarnings(do.call(
      n_alloc, c(list(fr, domains = "province", min_n_domain = 400), mode)
    ))
    expect_true(all(domain_totals(fit, "n_int") >= 400), info = names(mode))
    if (!is.null(mode$budget)) {
      expect_lte(fit$operational$cost, 60000)
    }
  }
})

test_that("predict() and a summary keep the floors", {
  fr <- floor_frame()
  fit <- n_alloc(fr, n = 1500, domains = "province", min_n_domain = 400)
  res <- predict(fit, data.frame(n = 1500))
  expect_equal(res$cv, fit$cv)
  expect_no_error(capture.output(print(fit), summary(fit)))
})

test_that("floors are refused where they cannot apply", {
  fr <- floor_frame()
  expect_error(n_alloc(fr, n = 1000, domains = "province", min_n_domain = 400),
               "below the 1200 the domain floors need")
  expect_error(n_alloc(fr, n = 1500, min_n_domain = 400),
               "needs 'domains'")
  expect_error(n_alloc(transform(fr, zone = "z"), n = 1500,
                       domains = c("province", "zone"), min_n_domain = 400),
               "one 'domains' column")
  expect_error(n_alloc(fr, n = 1500, domains = "province",
                       min_n_domain = 10000),
               "exceeds what domain")
  expect_error(n_alloc(fr, n = 1500, domains = "province",
                       min_n_domain = data.frame(domain = "province",
                                                 level = "Z", min_n = 10)),
               "not in 'province'")
  expect_error(n_alloc(fr, n = 1500, domains = "province", min_n_domain = -1),
               "non-negative")
  expect_warning(n_alloc(fr, n = 1500, domains = "province",
                         min_n_domain = 1),
                 "changes nothing")
  measures <- data.frame(stratum = fr$stratum, name = "y", p = fr$p)
  targets <- data.frame(name = "y", cv = 0.05)
  expect_error(n_alloc(fr, measures = measures, targets = targets,
                       domains = "province", min_n_domain = 400),
               "classic modes")
})

test_that("floors survive the round trip through prec_alloc()", {
  fr <- floor_frame()
  fit <- n_alloc(fr, n = 1500, domains = "province", min_n_domain = 400)
  fixed <- prec_alloc(fit)
  expect_equal(fixed$params$min_n_domain, 400)
  back <- n_alloc(fixed, n = 1500)
  expect_equal(domain_totals(back), domain_totals(fit))
  back <- n_alloc(fixed, cv = 0.09)
  expect_true(all(domain_totals(back) >= 400 - 1e-8))
  expect_equal(back$params$min_n_domain, 400)
  # A precision result built from sizes carries no floors.
  expect_null(prec_alloc(fr, n = fit$detail$n_int,
                         domains = "province")$params$min_n_domain)
})

test_that("the floor repair keeps every stratum on its n_step", {
  frame <- data.frame(stratum = c("a", "b"), region = c("A", "B"),
                      N = c(2000, 200), p = 0.3, n_step = 12)
  fit <- n_alloc(frame, n = 72, domains = "region", min_n_domain = 25)
  expect_equal(fit$detail$n_int, c(36, 36))

  # unequal steps in each mode
  frame <- data.frame(stratum = letters[1:4], region = rep(c("A", "B"), each = 2),
                      N = c(2000, 400, 150, 2000), p = c(0.2, 0.4, 0.3, 0.1),
                      unit_cost = c(1, 3, 2, 4), n_step = c(4, 6, 12, 1))
  for (args in list(list(n = 240), list(cv = 0.08), list(budget = 900))) {
    fit <- suppressWarnings(do.call(
      n_alloc, c(list(frame, domains = "region", min_n_domain = 80), args)
    ))
    n_int <- fit$detail$n_int
    expect_true(all(n_int %% frame$n_step == 0))
    expect_true(all(tapply(n_int, frame$region, sum) >= 80))
    expect_true(all(n_int <= frame$N))
    if (!is.null(args$n)) expect_equal(sum(n_int), 240)
    if (!is.null(args$budget)) expect_lte(sum(n_int * frame$unit_cost), 900)
  }
})

test_that("a total the floored steps reach only by moving strata together", {
  # domain B moves in steps of 6, so its floor of 40 costs 42, and A must
  # give back 2 with steps of 4 and 6: one more step of 4, one fewer of 6
  frame <- data.frame(stratum = letters[1:4], region = rep(c("A", "B"), each = 2),
                      N = c(2000, 2000, 150, 150),
                      p = c(0.2221891, 0.2567148, 0.1899314, 0.3787320),
                      unit_cost = c(4, 3, 3, 2), n_step = c(4, 6, 6, 6))
  fit <- suppressWarnings(
    n_alloc(frame, n = 276, domains = "region", min_n_domain = 40)
  )
  n_int <- fit$detail$n_int
  expect_equal(sum(n_int), 276)
  expect_true(all(n_int %% frame$n_step == 0))
  expect_true(all(tapply(n_int, frame$region, sum) >= 40))
  expect_true(all(abs(n_int - fit$detail$n) <= 2 * frame$n_step))
})
