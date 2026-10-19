## A user table read with `$` only for the columns it must hold

.strict_register <- function() {
  frame <- data.frame(
    stratum = c("A", "B"), N = c(4000, 8000),
    n_per_psu = 10, cost_psu = 20, cost_ssu = 1
  )
  psu <- data.frame(
    stratum = rep(c("A", "B"), each = 40),
    N = rep(c(100, 200), each = 40)
  )
  measures <- data.frame(
    stratum = c("A", "B"), name = "x", p = 0.3, icc_psu = 0.05
  )
  targets <- data.frame(name = "x", cv = 0.1)
  lapply(list(frame = frame, psu = psu, measures = measures,
              targets = targets), .strict_df)
}

test_that("the strict fixture warns on an absent column", {
  z <- .strict_register()
  expect_warning(z$psu$psu_id, "psu_id")
  expect_no_warning(z$psu$N)
})

test_that("a register without psu_id or certainty plans without a warning", {
  z <- .strict_register()
  expect_no_warning(
    fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu)
  )
  expect_false(any(fit$psu$certainty))
  expect_null(fit$psu$psu_id)
  expect_no_warning(prec_alloc(fit))
  expect_no_warning(
    predict(fit, newdata = .strict_df(data.frame(n_per_psu = c(8, 12))))
  )
})

test_that("a zoned register and a certainty cutoff plan without a warning", {
  z <- .strict_register()
  expect_no_warning(
    zoned <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                     psu = z$psu, n_psu_per_zone = 2)
  )
  expect_no_warning(prec_alloc(zoned))
  expect_no_warning(
    cut <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                   psu = z$psu, certainty_cutoff = 0.5)
  )
  expect_no_warning(
    predict(cut, newdata = .strict_df(
      data.frame(certainty_cutoff = c(0.5, 0.8))
    ))
  )
})

test_that("a supplied allocation is assessed from a register without psu_id", {
  z <- .strict_register()
  expect_no_warning(
    ev <- prec_alloc(z$frame, n = c(200, 300), measures = z$measures,
                     targets = z$targets, psu = z$psu)
  )
  expect_null(ev$psu$psu_id)
})

test_that("a PSU below the take is named by its row when psu_id is absent", {
  z <- .strict_register()
  psu <- z$psu
  psu$N[1:2] <- c(5, 195)
  expect_no_warning(
    expect_error(
      n_alloc(z$frame, measures = z$measures, targets = z$targets,
              psu = .strict_df(psu)),
      "row 1"
    )
  )
})

test_that("a register refuses measures without stratum before reading it", {
  z <- .strict_register()
  m <- .strict_df(z$measures[, setdiff(names(z$measures), "stratum")])
  expect_no_warning(
    expect_error(
      n_alloc(z$frame, measures = m, targets = z$targets, psu = z$psu),
      "'measures' must contain 'stratum' and 'name' columns"
    )
  )
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = z$psu)
  expect_no_warning(
    expect_error(
      prec_alloc(fit, measures = m),
      "'measures' must contain 'stratum' and 'name' columns"
    )
  )
})

test_that("a register with psu_id and certainty still carries both", {
  z <- .strict_register()
  psu <- z$psu
  psu$psu_id <- sprintf("p%02d", seq_len(nrow(psu)))
  psu$certainty <- seq_len(nrow(psu)) == 1L
  psu <- .strict_df(psu)
  fit <- n_alloc(z$frame, measures = z$measures, targets = z$targets,
                 psu = psu)
  expect_identical(fit$psu$psu_id, psu$psu_id)
  expect_true(fit$psu$certainty[1L])
  expect_true(fit$psu$.certainty_source[1L] %in% c("supplied", "threshold"))
})

test_that("indicators with moe targets and no var plan without a warning", {
  by_moe <- .strict_df(data.frame(
    name = c("stunting", "anemia"),
    p = c(0.30, 0.10),
    moe = c(0.05, 0.03),
    icc_psu = c(0.02, 0.05)
  ))
  expect_no_warning(n_cluster(indicators = by_moe, stage_cost = c(500, 50)))
})
