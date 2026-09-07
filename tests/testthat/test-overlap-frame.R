test_that("positive overlap must fit the finite responding frame", {
  for (rho in c(0, 0.5, 1)) {
    expect_error(
      prec_change(n=100, var=1, N=100, overlap=.5, overlap_cor=rho),
      "require 150 distinct units"
    )
  }
  for (rho in c(0, 0.1)) {
    expect_error(
      prec_pooled(n=90, var=1, N=100, occasions=2,
                  overlap=.5, overlap_cor=rho),
      "require 135 distinct units"
    )
  }
  expect_error(n_change(var=1, moe=.22, N=100, overlap=.5, overlap_cor=1),
               "overlap is infeasible")
  expect_error(n_pooled(var=1, moe=.01, N=100, occasions=2,
                        overlap=.5, overlap_cor=0), "overlap is infeasible")

  # 80 + 60 - .5*80 = 100, exactly the finite frame, with unequal waves.
  fit <- prec_change(var=1, n=c(80,60), N=100, overlap=.5, overlap_cor=.3)
  expect_gt(fit$se, 0)
  expect_equal(n_change(fit)$n, c(80,60), tolerance=1e-8)
  expect_error(prec_change(var=1, n=c(80,60.001), N=100, overlap=.5),
               "overlap is infeasible")
  # Responding counts, not issued counts, determine the union.
  net <- prec_change(var=1, n=60, N=100, overlap=.5, overlap_cor=.4)
  gross <- prec_change(var=1, n=100, resp_rate=.6, N=100,
                       overlap=.5, overlap_cor=.4)
  expect_equal(gross$se, net$se)
  expect_equal(prec_change(var=1,n=100,N=100,overlap=1)$se, 0)
  expect_equal(prec_change(var=1,n=100,N=100,overlap=0)$se, 0)
  expect_s3_class(prec_change(var=1,n=1000,overlap=.5), "svyplan_prec")
})

test_that("power evaluation and MDE reject contradictory frame overlaps", {
  for (solve_power in c(TRUE, FALSE)) {
    expect_error(power_mean(var=1, n=100, N=100, overlap=.5,
                            effect=if(solve_power) .2 else NULL,
                            power=if(solve_power) NULL else .8),
                 "overlap is infeasible")
    expect_error(power_prop(.2, n=100, N=100, overlap=.5,
                            p2=if(solve_power) .3 else NULL,
                            power=if(solve_power) NULL else .8),
                 "overlap is infeasible")
    expect_error(power_did(c(0,.2), control=c(0,0), var=1,
                           n=100, N=100, overlap=.5,
                           effect=if(solve_power) .2 else NULL,
                           power=if(solve_power) NULL else .8),
                 "overlap is infeasible")
  }
})

test_that("power sizing brackets only feasible sizes including the endpoint", {
  # Maximum size under half overlap is 100/1.5, not the census of 100.
  n_max <- 100/1.5
  for (alternative in c("one.sided", "two.sided")) {
    at_bound <- power_mean(var=1, effect=.2, n=n_max, power=NULL,
                           N=100, overlap=.5, overlap_cor=1,
                           alternative=alternative)
    fit <- power_mean(var=1, effect=.2, power=at_bound$power,
                      N=100, overlap=.5, overlap_cor=1,
                      alternative=alternative)
    expect_equal(fit$n, n_max)
    expect_error(power_mean(var=1,effect=.2,power=.8,N=100,
                            overlap=.5,overlap_cor=1,alternative=alternative),
                 "unattainable.*overlap")
    prop <- power_prop(.3,p2=.4,n=n_max,power=NULL,N=100,
                       overlap=.5,overlap_cor=.3,alternative=alternative)
    back <- power_prop(.3,p2=.4,power=prop$power,N=100,
                       overlap=.5,overlap_cor=.3,alternative=alternative)
    expect_equal(back$n,n_max)
  }
  # Even a huge effect cannot justify the n=2 floor if the union cannot fit.
  expect_error(power_mean(var=1,effect=100,power=.8,N=2,overlap=.5),
               "unattainable.*overlap")
  expect_error(power_prop(.01,p2=.99,power=.8,N=2,overlap=.5),
               "unattainable.*overlap")
  expect_equal(power_mean(var=1,effect=100,power=.8,N=2,overlap=1)$n,2)
})

test_that("DiD overlap is within each arm and supports unequal arm sizes", {
  # Full within-arm panels are valid even with n_treat/n_control > 1.
  fit <- power_did(c(0,.2),control=c(0,0),var=1,n=c(80,40),
                    power=NULL,N=c(100,60),overlap=1,overlap_cor=.5)
  expect_gt(fit$power, .05)
  back <- power_did(c(0,.2),control=c(0,0),var=1,power=fit$power,
                    ratio=2,N=c(100,60),overlap=1,overlap_cor=.5)
  expect_equal(unname(back$n),c(80,40),tolerance=1e-7)
  expect_error(power_did(c(0,.2),control=c(0,0),var=1,n=c(80,40),
                         power=NULL,N=c(100,60),overlap=.5),
               "require 120 distinct units")
  # The control arm binds even though the treated arm has more respondents.
  expect_error(power_did(c(0,.2),control=c(0,0),var=1,n=c(60,50),
                         power=NULL,N=c(100,60),overlap=.5),
               "require 75 distinct units")
  at_bound <- power_did(c(0,.2),control=c(0,0),var=1,n=c(60,30),
                        power=NULL,N=c(90,60),overlap=.5,overlap_cor=.5)
  back <- power_did(c(0,.2),control=c(0,0),var=1,power=at_bound$power,
                    ratio=2,N=c(90,60),overlap=.5,overlap_cor=.5)
  expect_equal(unname(back$n),c(60,30),tolerance=1e-7)
})

test_that("pooled feasibility checks later lags even at zero correlation", {
  expect_error(prec_pooled(var=1,n=80,N=100,occasions=3,
                           overlap=c(.9,.5),overlap_cor=0),
               "require 120 distinct units")
  expect_s3_class(prec_pooled(var=1,n=80,N=100,occasions=3,
                              resp_rate=.5,overlap=c(.9,.5),overlap_cor=0),
                  "svyplan_prec")
})

test_that("feasible change variance agrees with exhaustive coordinated SRS", {
  y1 <- c(1,2,4,5,8,10)
  y2 <- c(2,4,3,8,7,12)
  samples <- combn(6,4,simplify=FALSE)
  changes <- numeric()
  for (s1 in samples) for (s2 in samples) {
    if (length(intersect(s1,s2)) == 2) {
      changes <- c(changes,mean(y2[s2])-mean(y1[s1]))
    }
  }
  fit <- prec_change(n=4,var=c(var(y1),var(y2)),N=6,
                     overlap=.5,overlap_cor=cor(y1,y2))
  expect_equal(fit$se^2,mean((changes-mean(changes))^2),tolerance=1e-10)
})
