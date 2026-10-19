## A domain without a level column stands for every level of that domain

.every_level_fixture <- function() {
  frame <- data.frame(
    stratum = c("N-u", "N-r", "S-u", "S-r", "E-u", "E-r"),
    region = rep(c("North", "South", "East"), each = 2),
    residence = rep(c("urban", "rural"), 3),
    N = c(1000, 1800, 1200, 900, 700, 1500)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("mcpr", "stunting"), each = 6),
    p = c(0.25, 0.20, 0.30, 0.22, 0.28, 0.18,
          0.20, 0.30, 0.25, 0.35, 0.22, 0.33)
  )
  every <- data.frame(
    name = c("mcpr", "mcpr", "mcpr", "stunting"),
    domain = c(".overall", "region", "residence", "region"),
    cv = c(0.05, 0.10, 0.08, 0.12)
  )
  full <- data.frame(
    name = c(rep("mcpr", 6), rep("stunting", 3)),
    domain = c(".overall", rep("region", 3), rep("residence", 2),
               rep("region", 3)),
    level = c(NA, "North", "South", "East", "urban", "rural",
              "North", "South", "East"),
    cv = c(0.05, rep(0.10, 3), rep(0.08, 2), rep(0.12, 3))
  )
  list(frame = frame, measures = measures, every = every, full = full)
}

test_that("targets without a level column expand to the hand-built table", {
  z <- .every_level_fixture()
  every <- n_alloc(z$frame, measures = z$measures, targets = z$every)
  full <- n_alloc(z$frame, measures = z$measures, targets = z$full)

  expect_identical(every$constraints$constraint, full$constraints$constraint)
  expect_identical(every$detail$n_int, full$detail$n_int)
  expect_identical(every$params$targets, full$params$targets)
  expect_identical(
    every$constraints$constraint[2:4],
    c("mcpr@region=North:cv", "mcpr@region=South:cv", "mcpr@region=East:cv")
  )
})

test_that("an every-level fit round-trips as the hand-built one does", {
  z <- .every_level_fixture()
  every <- n_alloc(z$frame, measures = z$measures, targets = z$every)
  full <- n_alloc(z$frame, measures = z$measures, targets = z$full)

  expect_identical(prec_alloc(every)$constraints$.achieved,
                   prec_alloc(full)$constraints$.achieved)
  expect_identical(n_alloc(prec_alloc(every))$detail$n,
                   n_alloc(prec_alloc(full))$detail$n)
})

test_that("an objective without a level column expands the same way", {
  z <- .every_level_fixture()
  short <- n_alloc(
    z$frame, measures = z$measures, targets = z$every[1, ],
    objective = data.frame(name = "stunting", domain = "region"),
    budget = 4000
  )
  long <- n_alloc(
    z$frame, measures = z$measures, targets = z$every[1, ],
    objective = data.frame(name = "stunting", domain = "region",
                           level = c("North", "South", "East")),
    budget = 4000
  )
  expect_identical(short$params$objective$level, c("North", "South", "East"))
  expect_identical(short$detail$n_int, long$detail$n_int)
  expect_identical(short$params$objective, long$params$objective)
})

test_that("a level column with a missing level on a domain row is refused", {
  z <- .every_level_fixture()
  mixed <- rbind(
    cbind(z$every[1:2, ], level = NA_character_),
    data.frame(name = "stunting", domain = "region", cv = 0.12,
               level = "North")
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = mixed),
    "missing for domain .region.*Leave out the 'level' column.*cannot mix"
  )
  expect_error(
    n_alloc(
      z$frame, measures = z$measures, targets = z$every[1, ],
      objective = data.frame(name = "stunting", domain = "region",
                             level = NA_character_),
      budget = 4000
    ),
    "'objective\\$level' is missing for domain"
  )
})

test_that("an expanded row cannot carry a single identifier", {
  z <- .every_level_fixture()
  labelled <- z$every
  labelled$constraint <- paste0("c", 1:4)
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = labelled),
    "'targets\\$constraint' labels one row"
  )
  expect_error(
    n_alloc(z$frame, measures = z$measures,
            targets = data.frame(name = "mcpr", level = "North", cv = 0.1)),
    "'targets\\$level' requires a 'domain' column"
  )
})

test_that("an every-level row checks its domain column in the frame", {
  z <- .every_level_fixture()
  absent <- data.frame(name = "mcpr", domain = "district", cv = 0.1)
  expect_error(
    n_alloc(z$frame, measures = z$measures, targets = absent),
    "domain column 'district' not found in frame"
  )
  holed <- z$frame
  holed$region[2] <- NA
  expect_error(
    n_alloc(holed, measures = z$measures, targets = z$every),
    "domain column 'region' must not contain missing values"
  )
})

test_that("every-level targets read no absent column", {
  z <- .every_level_fixture()
  expect_no_warning(
    n_alloc(.strict_df(z$frame), measures = .strict_df(z$measures),
            targets = .strict_df(z$every))
  )
})
