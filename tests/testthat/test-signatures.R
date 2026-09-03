test_that("calculation methods place dots after required arguments", {
  expected <- c(
    n_prop.default = 2L,
    prec_prop.default = 3L,
    n_mean.default = 2L,
    prec_mean.default = 3L,
    n_cluster.default = 2L,
    prec_cluster.default = 2L,
    n_multi.default = 2L,
    prec_multi.default = 2L,
    n_alloc.default = 2L,
    prec_alloc.default = 3L,
    n_panel = 3L,
    prec_panel.default = 4L,
    power_prop.default = 2L,
    power_mean.default = 2L,
    power_did.default = 3L
  )

  for (nm in names(expected)) {
    fn <- get(nm, mode = "function")
    expect_identical(match("...", names(formals(fn))), unname(expected[[nm]]))
  }
})

test_that("generic signatures expose required and optional inputs", {
  expect_identical(names(formals(power_mean))[1L], "var")
  expect_identical(names(formals(power_mean.default))[1L], "var")
  expect_null(formals(power_mean.default)$effect)
  expect_null(formals(n_cluster)$stage_cost)
  expect_null(formals(n_cluster.default)$stage_cost)
  # A panel is entered from the object it serves, and read back from the
  # recruitment it asked for.
  expect_identical(names(formals(n_panel))[1L], "target")
  expect_identical(names(formals(prec_panel))[1L], "n_recruit")
  expect_identical(names(formals(prec_panel.default))[2L], "target")
  expect_identical(
    formals(n_panel)$design, quote(c("fixed", "rotating"))
  )
  expect_true(identical(formals(strata_bound)$n_strata, quote(expr = )))
})

test_that("short enumerations advertise their choices", {
  prop_choices <- quote(c("wald", "wilson", "logodds", "beta"))
  expect_identical(formals(n_multi.default)$prop_method, prop_choices)
  expect_identical(formals(prec_multi.default)$prop_method, prop_choices)
  # the multi family offers exactly what the single-indicator engines do
  expect_identical(formals(n_prop.default)$method, prop_choices)
  expect_identical(formals(prec_prop.default)$method, prop_choices)
  expect_identical(
    formals(strata_bound)$method,
    quote(c("lh", "cumrootf", "geo", "kozak"))
  )
  expect_identical(
    formals(strata_bound)$alloc,
    quote(c("neyman", "optimal", "proportional", "power"))
  )
})
