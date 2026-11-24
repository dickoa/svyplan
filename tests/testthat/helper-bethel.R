.bethel_fixture <- function() {
  frame <- data.frame(
    stratum = c("NU", "NR", "SU", "SR"),
    region = c("North", "North", "South", "South"),
    residence = c("Urban", "Rural", "Urban", "Rural"),
    N = c(1000, 2000, 1500, 1000),
    unit_cost = c(1, 1.2, 1.5, 1)
  )
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("vaccination", "income"), each = 4),
    p = c(0.5, 0.4, 0.6, 0.3, rep(NA, 4)),
    mean = c(rep(NA, 4), 50, 55, 60, 45),
    sd = c(rep(NA, 4), 10, 12, 15, 9)
  )
  targets <- data.frame(
    name = c("vaccination", "vaccination", "income"),
    domain = c(".overall", "region", "residence"),
    level = c(NA, "North", "Urban"),
    cv = c(0.05, 0.08, NA),
    moe = c(NA, NA, 2)
  )
  list(frame = frame, measures = measures, targets = targets)
}

.bethel_multistage_fixture <- function(stages = 2L) {
  z <- .bethel_fixture()
  z$frame$unit_cost <- NULL
  z$frame$N_1 <- c(100, 160, 120, 90)
  z$frame$take_2 <- c(8, 10, 12, 7)
  z$frame$cost_1 <- c(300, 400, 450, 350)
  z$frame$cost_2 <- c(25, 30, 35, 28)
  z$measures$icc_1 <- rep(c(0.03, 0.05, 0.08, 0.04), 2)
  z$measures$var_ratio_1 <- rep(c(1, 1.1, 0.9, 1.2), 2)
  if (stages == 3L) {
    z$frame$N_2 <- c(600, 1200, 900, 500)
    z$frame$take_3 <- c(4, 5, 3, 4)
    z$frame$cost_3 <- c(5, 6, 7, 5)
    z$measures$icc_2 <- rep(c(0.10, 0.08, 0.12, 0.06), 2)
    # var_ratio_2 is fixed by the decomposition: var_ratio_1 * (1 - icc_1).
    z$measures$var_ratio_2 <- z$measures$var_ratio_1 * (1 - z$measures$icc_1)
  }
  z
}
