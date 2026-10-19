## A frame's label columns travel into the allocation's stratum table

.label_frame <- function() {
  data.frame(
    region = c("R1", "R1", "R2", "R2"),
    residence = c("u", "r", "u", "r"),
    N = c(1e5, 5e4, 3e4, 6e4),
    sd = 1, mean = 1
  )
}

test_that("label columns follow stratum in detail and as.data.frame()", {
  fr <- .label_frame()
  x <- n_alloc(fr, n = 400)
  expect_identical(names(x$detail)[1:3], c("stratum", "region", "residence"))
  expect_identical(x$detail$region, fr$region)
  expect_identical(x$detail$residence, fr$residence)
  expect_identical(names(as.data.frame(x))[1:3],
                   c("stratum", "region", "residence"))
  expect_identical(x$params$strata_keys, c("region", "residence"))
  expect_identical(x$detail$n_int,
                   n_alloc(fr[c("N", "sd", "mean")], n = 400)$detail$n_int)
})

test_that("keys that do not identify the strata are not recorded", {
  fr <- .label_frame()[c("region", "N", "sd", "mean")]
  x <- n_alloc(fr, n = 400)
  expect_identical(x$detail$region, fr$region)
  expect_null(x$params$strata_keys)
})

test_that("numeric columns stay behind and raise nothing", {
  fr <- transform(.label_frame(), weight = 3, n = 9, v024 = c(1L, 1L, 2L, 2L))
  expect_no_error(x <- n_alloc(fr, n = 400))
  expect_false("v024" %in% names(x$detail))
  expect_identical(names(x$detail)[1:3], c("stratum", "region", "residence"))

  coded <- transform(.label_frame(), v024 = factor(c(1, 1, 2, 2)))
  expect_identical(n_alloc(coded, n = 400)$detail$v024, coded$v024)
})

test_that("a logical flag travels and take_all stays an input", {
  fr <- transform(.label_frame(), urban = residence == "u", take_all = FALSE)
  x <- n_alloc(fr, n = 400)
  expect_identical(x$detail$urban, fr$urban)
  expect_false("take_all" %in% names(x$detail))
})

test_that("a label column named like a detail column is refused", {
  expect_error(
    n_alloc(transform(.label_frame(), n_int = "a"), n = 400),
    "label column .n_int. shares a name"
  )
})

test_that("prec_alloc, rebuilds and domains carry the labels", {
  fr <- .label_frame()
  x <- n_alloc(fr, n = 400)
  ev <- prec_alloc(x)
  expect_identical(ev$detail$region, fr$region)
  expect_identical(ev$params$strata_keys, c("region", "residence"))
  expect_identical(n_alloc(ev)$detail$residence, fr$residence)
  expect_identical(n_alloc(fr, n = 400, domains = "region")$detail$region,
                   fr$region)
})

test_that("cluster, joint and register allocations carry the labels", {
  fr <- .label_frame()
  cl <- transform(fr, N_psu = c(200, 100, 60, 120), n_per_psu = 10,
                  icc_psu = 0.05)
  expect_identical(n_alloc(cl, n = 400)$detail$region, fr$region)

  jf <- data.frame(stratum = paste(fr$region, fr$residence),
                   region = fr$region, residence = fr$residence, N = fr$N)
  jm <- data.frame(stratum = jf$stratum, name = "y", p = 0.3)
  tg <- data.frame(name = "y", cv = 0.05)
  joint <- n_alloc(jf, measures = jm, targets = tg)
  expect_identical(joint$detail$residence, fr$residence)
  expect_identical(joint$params$strata_keys, c("region", "residence"))

  jf$n_per_psu <- 10
  jf$cost_psu <- 20
  jf$cost_ssu <- 1
  reg <- data.frame(stratum = rep(jf$stratum, each = 20),
                    N = rep(jf$N / 20, each = 20))
  register <- n_alloc(jf, measures = transform(jm, icc_psu = 0.05),
                      targets = tg, psu = reg)
  expect_identical(register$detail$region, fr$region)
  expect_true("n_psu_certain" %in% names(register$detail))
})

test_that("a strict frame carries its labels without an absent-column read", {
  expect_no_warning(x <- n_alloc(.strict_df(.label_frame()), n = 400))
  expect_identical(x$detail$region, .label_frame()$region)
})
