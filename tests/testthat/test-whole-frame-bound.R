## A solved size is refused when its whole size breaks the frame bound

# Closed form (oracle rung 1): two responding samples of n that share a
# fraction o of their units need (2 - o) * n distinct units, so a whole size
# fits a frame of N exactly when (2 - o) * ceiling(n) <= N.
fits_whole <- function(n, overlap, N) (2 - overlap) * ceiling(n) <= N

# Each route: the size that a precision call at n0 reports as its target,
# the sizing call that solves it back, and the evaluator at a supplied n.
pc <- function(n) prec_change(p = c(0.3, 0.35), n = n, N = 1001, overlap = 0.5,
                              overlap_cor = 0.5)
pp <- function(n, overlap = 0.5) {
  prec_pooled(var = 1, n = n, N = 1001, occasions = 2, overlap = overlap,
              overlap_cor = 0.5)
}
pwp <- function(n) power_prop(p1 = 0.3, p2 = 0.36, n = n, N = 1001,
                              overlap = 0.5, overlap_cor = 0.5, power = NULL)
pwm <- function(n) power_mean(effect = 2, var = 100, n = n, N = 1001,
                              overlap = 0.5, overlap_cor = 0.5, power = NULL)
pwd <- function(n) power_did(treat = c(0.3, 0.4), control = c(0.3, 0.32),
                             outcome = "prop", n = n, N = 1001, overlap = 0.5,
                             overlap_cor = 0.4, power = NULL)
rotation <- design_overlap(c(1, 0, 1), max_lag = 1)
routes <- list(
  n_change = list(
    eval = pc, overlap = 0.5,
    solve = function(n0) n_change(p = c(0.3, 0.35), moe = pc(n0)$moe,
                                  N = 1001, overlap = 0.5, overlap_cor = 0.5)
  ),
  n_pooled = list(
    eval = pp, overlap = 0.5,
    solve = function(n0) n_pooled(var = 1, moe = pp(n0)$moe, N = 1001,
                                  occasions = 2, overlap = 0.5,
                                  overlap_cor = 0.5)
  ),
  n_pooled_rotation = list(
    eval = function(n) pp(n, rotation), overlap = 0,
    solve = function(n0) n_pooled(var = 1, moe = pp(n0, rotation)$moe,
                                  N = 1001, occasions = 2, overlap = rotation,
                                  overlap_cor = 0.5)
  ),
  power_prop = list(
    eval = pwp, overlap = 0.5,
    solve = function(n0) power_prop(p1 = 0.3, p2 = 0.36, N = 1001,
                                    overlap = 0.5, overlap_cor = 0.5,
                                    power = pwp(n0)$power, n = NULL)
  ),
  power_mean = list(
    eval = pwm, overlap = 0.5,
    solve = function(n0) power_mean(effect = 2, var = 100, N = 1001,
                                    overlap = 0.5, overlap_cor = 0.5,
                                    power = pwm(n0)$power, n = NULL)
  ),
  power_did = list(
    eval = pwd, overlap = 0.5,
    solve = function(n0) power_did(treat = c(0.3, 0.4),
                                   control = c(0.3, 0.32), outcome = "prop",
                                   N = 1001, overlap = 0.5, overlap_cor = 0.4,
                                   power = pwd(n0)$power, n = NULL)
  )
)

test_that("the review's disjoint pair is refused at its whole size", {
  moe <- stats::qnorm(0.975) * sqrt(1 / (2 * 500.25) - 1 / 1001)
  expect_false(fits_whole(500.25, 0, 1001))
  expect_error(
    n_pooled(var = 1, N = 1001, moe = moe, occasions = 2, overlap = rotation,
             overlap_cor = 1),
    "n = 500.25 per occasion, which the frame holds only before rounding up"
  )
})

test_that("a solved size is refused exactly when its whole size cannot fit", {
  seen <- character(0)
  for (route in names(routes)) {
    r <- routes[[route]]
    edge <- 1001 / (2 - r$overlap)
    for (n0 in edge - c(2.6, 1.3, 0.9, 0.5, 0.2)) {
      label <- sprintf("%s at %.2f", route, n0)
      if (fits_whole(n0, r$overlap, 1001)) {
        fit <- r$solve(n0)
        expect_equal(as.double(fit)[1L], n0, tolerance = 1e-6, label = label)
        expect_no_error(r$eval(as.integer(fit)))
        seen <- c(seen, paste(route, "fits"))
      } else {
        expect_error(r$solve(n0), "only before rounding up", label = label)
        expect_error(r$eval(ceiling(n0)), "infeasible|not a valid covariance",
                     label = label)
        seen <- c(seen, paste(route, "refused"))
      }
    }
  }
  expect_setequal(seen, paste(rep(names(routes), each = 2),
                              c("fits", "refused")))
})

test_that("independent draws keep sizes the frame bound does not reach", {
  n0 <- 667.2
  expect_false(fits_whole(n0, 0.5, 1001))
  moe <- pc(n0)$moe
  fit <- n_change(p = c(0.3, 0.35), moe = moe, N = 1001, overlap = 0.5,
                  overlap_cor = 0.5, overlap_draw = "independent")
  expect_no_error(prec_change(fit, n = as.integer(fit)))
})
