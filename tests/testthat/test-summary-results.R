## Structured summaries of sizes and cluster designs

z975 <- stats::qnorm(0.975)

# Closed form (oracle rung 1): Wald margin of a proportion on n_resp
# respondents, with the correction (N - n) / (N - 1) on respondents.
wald_moe <- function(p, n_resp, deff = 1, N = Inf, z = z975) {
  fpc <- if (is.finite(N)) (N - n_resp) / (N - 1) else 1
  z * sqrt(deff * p * (1 - p) / n_resp * fpc)
}

# Closed form (oracle rung 1): two-stage cv with var_ratio 1 on a realized
# take m_r = n_per_psu times the unit rates, without stage corrections.
cluster_cv <- function(relvar, icc, n_psu, m_r, resp_psu = 1) {
  sqrt(relvar / (n_psu * resp_psu * m_r) * (1 + icc * (m_r - 1)))
}

## S1 Scalar sizes

test_that("a proportion summary reads sizes and precision at the field size", {
  x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8, N = 10000)
  s <- summary(x)
  expect_s3_class(s, c("summary.svyplan_n", "list"), exact = TRUE)
  expect_identical(s$kind, "scalar")
  expect_identical(s$sizes[["field"]], 578)
  expect_equal(s$sizes[["continuous"]], x$n)
  expect_equal(s$sizes[["respondents"]], 578 * 0.8)
  expect_equal(s$sizes[["effective"]], 578 * 0.8 / 1.5)
  moe <- wald_moe(0.3, 578 * 0.8, deff = 1.5, N = 10000)
  expect_equal(s$precision[["moe"]], moe, tolerance = 1e-10)
  expect_equal(s$precision[["se"]], moe / z975, tolerance = 1e-10)
  expect_equal(s$precision[["cv"]], moe / z975 / 0.3, tolerance = 1e-10)
  expect_lt(s$precision[["moe"]], 0.05)
  expect_equal(s$expected_cases, 578 * 0.8 * 0.3)
  expect_identical(s$target, list(scale = "moe", value = 0.05))
  expect_null(s$interval)

  out <- capture.output(print(s))
  expect_true("Question: the n that gives moe = 0.05 for p = 0.30" %in% out)
  expect_true(any(grepl("^  Field +578$", out)))
  expect_true(any(grepl("^  Expected respondents +462\\.4$", out)))
  expect_true("  population size: 10000" %in% out)
})

test_that("a field moe never exceeds the target, over methods and rates", {
  seen <- character(0)
  for (method in c("wald", "wilson", "logodds", "beta")) {
    for (p in c(0.05, 0.3, 0.5)) {
      for (resp in c(1, 0.7)) {
        for (N in c(Inf, 5000)) {
          x <- n_prop(p = p, moe = 0.03, method = method, resp_rate = resp,
                      N = N)
          s <- summary(x)
          label <- paste(method, p, resp, N)
          expect_lte(s$precision[["moe"]], 0.03 + 1e-12, label = label)
          expect_identical(s$sizes[["field"]], ceiling(x$n), label = label)
          seen <- c(seen, method)
        }
      }
    }
  }
  expect_setequal(unique(seen), c("wald", "wilson", "logodds", "beta"))
  expect_length(seen, 48L)
})

test_that("an interval is read at the level the size was planned for", {
  x <- n_prop(p = 0.1, moe = 0.03, method = "wilson", alpha = 0.10)
  s <- summary(x)
  n <- ceiling(x$n)
  z <- stats::qnorm(0.95)
  # Closed form (oracle rung 1): Wilson score interval.
  centre <- (0.1 + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z / (1 + z^2 / n) * sqrt(0.1 * 0.9 / n + z^2 / (4 * n^2))
  expect_equal(s$level, 0.90)
  expect_equal(s$interval, c(centre - half, centre + half), tolerance = 1e-8)
  expect_true(any(grepl("\\(90%\\)$", capture.output(print(s)))))

  logodds <- summary(n_prop(p = 0.1, moe = 0.03, method = "logodds",
                            alpha = 0.10))
  expect_length(logodds$interval, 2L)
  expect_gt(0.1 - logodds$interval[1L], 0)
  expect_false(isTRUE(all.equal(0.1 - logodds$interval[1L],
                                logodds$interval[2L] - 0.1)))
})

test_that("a scalar summary states the case floor and the absent scales", {
  rare <- summary(n_prop(p = 0.02, moe = 0.02, min_cases = 30))
  expect_identical(rare$binding, "min_cases")
  expect_true("binding: min_cases (min_cases = 30)" %in%
                capture.output(print(rare)))

  mean_moe <- summary(n_mean(var = 100, moe = 2))
  expect_true(is.na(mean_moe$precision[["cv"]]))
  expect_false(any(grepl("^  CV", capture.output(print(mean_moe)))))

  ratio <- summary(n_ratio(r = 420, cv_num = 1.2, cv_den = 0.45,
                           component_cor = 0.65, rmoe = 0.05))
  expect_identical(ratio$target$scale, "rmoe")
  expect_lte(ratio$precision[["rmoe"]], 0.05 + 1e-12)
  expect_true("  component_cor: 0.65" %in% capture.output(print(ratio)))

  epu <- summary(n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9,
                        eligible_per_unit = 0.41))
  expect_equal(epu$sizes[["eligible"]], 6969 * 0.9 * 0.41)
  expect_equal(epu$sizes[["effective"]], 6969 * 0.9 * 0.41 / 2.5)
})

## S2 Single-indicator cluster

test_that("a cluster summary reads the field design against closed forms", {
  x <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
                 eligibility_rate = 0.95, resp_rate = 0.9)
  s <- summary(x)
  expect_s3_class(s, c("summary.svyplan_cluster", "list"), exact = TRUE)
  a <- x$operational$n[["stage_1"]]
  m <- x$operational$n[["stage_2"]]
  m_r <- m * 0.9 * 0.95
  d <- s$design
  expect_equal(d$field[d$row == "cv"], cluster_cv(1, 0.05, a, m_r),
               tolerance = 1e-10)
  expect_equal(d$continuous[d$row == "cv"], x$cv)
  expect_equal(d$field[d$row == "total_n"], a * m)
  expect_equal(s$field[["respondents"]], a * m_r)
  expect_equal(s$field[["deff"]], 1 + 0.05 * (m_r - 1), tolerance = 1e-10)
  expect_equal(s$field[["effective"]],
               a * m_r / (1 + 0.05 * (m_r - 1)), tolerance = 1e-10)
  expect_equal(s$field[["df"]], a - 1)
  expect_equal(s$cost, c(stage_1 = a * 500, stage_2 = a * m * 50))
  expect_equal(sum(s$cost), x$operational$cost)
})

test_that("the cost by stage adds up on every cluster route", {
  fits <- list(
    two_stage_fixed = n_cluster(stage_cost = c(500, 50), icc = 0.05,
                                budget = 1e5, fixed_cost = 5000),
    three_stage = n_cluster(stage_cost = c(500, 100, 50),
                            icc = c(0.01, 0.05), budget = 5e5,
                            stage_n = c(50, NA, NA)),
    stage_fpc = n_cluster(stage_cost = c(12, 3), icc = 5 / 54,
                          unit_relvar = 0.54, budget = 1000, fpc = "stage",
                          N = c(100, 4900))
  )
  seen <- character(0)
  for (route in names(fits)) {
    s <- summary(fits[[route]])
    expect_equal(sum(s$cost), fits[[route]]$operational$cost, label = route)
    seen <- c(seen, route)
  }
  expect_setequal(seen, names(fits))
  expect_named(summary(fits$three_stage)$cost, c("stage_1", "stage_2", "stage_3"))
})

test_that("the design effect under stage correction is the planning one", {
  x <- n_cluster(stage_cost = c(12, 3), icc = 5 / 54, unit_relvar = 0.54,
                 budget = 1000, fpc = "stage", N = c(100, 4900))
  s <- summary(x)
  m <- x$operational$n[["stage_2"]]
  expect_equal(s$field[["deff"]], 1 + 5 / 54 * (m - 1), tolerance = 1e-10)
  srs_cv2 <- 0.54 / x$operational$total_n
  expect_false(isTRUE(all.equal(s$field[["deff"]],
                                x$operational$cv^2 / srs_cv2)))
  expect_true(any(grepl("Design effect (no FPC)",
                        capture.output(print(s)), fixed = TRUE)))
})

## S3 Multi-indicator sizes

test_that("each indicator is read at the field size on its stated scale", {
  ind <- data.frame(name = c("stunting", "vaccination", "anemia"),
                    p = c(0.30, 0.70, 0.10), moe = c(0.05, 0.05, 0.03),
                    deff = c(2, 1.5, 2.5))
  x <- n_multi(ind, resp_rate = 0.9)
  s <- summary(x)
  expect_identical(s$kind, "multi")
  expect_identical(s$sizes[["field"]], ceiling(x$n))
  expect_equal(s$sizes[["respondents"]], ceiling(x$n) * 0.9)
  expected <- wald_moe(ind$p, ceiling(x$n) * 0.9, deff = ind$deff)
  expect_equal(s$indicators$achieved, expected, tolerance = 1e-10)
  expect_equal(x$detail$.cv_field, expected / z975 / ind$p,
               tolerance = 1e-10)
  expect_identical(s$binding, "anemia")
  expect_identical(s$binding_cause, "precision")

  rmoe <- summary(n_multi(data.frame(name = c("a", "b"), p = c(0.3, 0.1),
                                     rmoe = 0.12)))
  expect_identical(rmoe$indicators$target_scale, c("rmoe", "rmoe"))
  expect_equal(rmoe$indicators$target, c(0.12, 0.12))
  n <- rmoe$sizes[["field"]]
  expect_equal(rmoe$indicators$achieved,
               wald_moe(c(0.3, 0.1), n) / c(0.3, 0.1), tolerance = 1e-10)
  expect_true(any(grepl("rmoe 0.12", capture.output(print(rmoe)))))
})

test_that("the binding cause comes from the solver", {
  floors <- n_multi(data.frame(
    name = c("rare", "common", "floor_row"),
    p = c(0.02, 0.30, NA), moe = c(0.02, 0.05, NA),
    min_cases = c(30, NA, NA), min_respondents = c(NA, NA, 2000)
  ))
  expect_identical(floors$detail$.constraint,
                   c("min_cases", "precision", "min_respondents"))
  expect_identical(floors$binding_cause, "min_respondents")
  s <- summary(floors)
  expect_identical(nrow(s$indicators), 3L)
  expect_true(is.na(s$indicators$achieved[3L]))

  cases <- n_multi(data.frame(name = c("rare", "common"), p = c(0.02, 0.30),
                              moe = c(0.02, 0.05), min_cases = c(30, NA)))
  expect_identical(cases$binding_cause, "min_cases")
  expect_gt(cases$detail$.cv_target[1L], cases$detail$.cv_achieved[1L])

  dom <- data.frame(name = c("a", "a"), region = c("N", "S"),
                    p = c(0.3, 0.35), moe = c(0.05, 0.05))
  floored <- n_multi(dom, domains = "region", min_n_domain = 600)
  expect_identical(floored$binding_cause, "min_n_domain")
})

test_that("domain rows are read at the field size of their domain", {
  dom <- data.frame(name = c("a", "b", "a", "b"),
                    region = c("N", "N", "S", "S"),
                    p = c(0.30, 0.10, 0.35, 0.12),
                    moe = c(0.05, 0.03, 0.05, 0.03))
  sep <- n_multi(dom, domains = "region")
  s <- summary(sep)
  quota <- ceiling(sep$domains$.n)
  expect_equal(s$domains$field, quota)
  expect_identical(s$sizes[["field"]], sum(quota))
  row_field <- quota[match(dom$region, sep$domains$region)]
  expect_equal(s$indicators$achieved, wald_moe(dom$p, row_field),
               tolerance = 1e-10)
  expect_identical(nrow(sep$detail), 4L)
  expect_identical(names(sep$detail)[1:2], c("region", "name"))

  nat <- n_multi(transform(dom, share = ifelse(region == "N", 0.4, 0.6)),
                 domains = "region", domain_sampling = "natural")
  sn <- summary(nat)
  expect_equal(sn$domains$field, ceiling(nat$n) * c(0.4, 0.6))
  expect_equal(sn$indicators$achieved,
               wald_moe(dom$p, ceiling(nat$n) * rep(c(0.4, 0.6), each = 2)),
               tolerance = 1e-10)
})

test_that("row-specific rates are reported per row, common ones once", {
  varied <- summary(n_multi(data.frame(
    name = c("a", "b"), p = c(0.3, 0.1), moe = c(0.05, 0.03),
    resp_rate = c(0.9, 0.7), deff = c(1.5, 2)
  )))
  expect_true(is.na(varied$sizes["respondents"]))
  expect_equal(varied$indicators$respondents,
               varied$sizes[["field"]] * c(0.9, 0.7))
  out <- capture.output(print(varied))
  expect_true(any(grepl("Resp\\.", out)))
  expect_true(any(grepl("Deff", out)))
})

## S4 Multi-indicator cluster

test_that("a multi-indicator cluster summary reports worst cv and binding apart", {
  ab <- data.frame(name = c("a", "b"), p = c(0.35, 0.40), cv = c(0.15, 0.03),
                   icc_1 = 0.05)
  x <- n_cluster(ab, stage_cost = c(500, 50), budget = 1e5)
  s <- summary(x)
  a <- x$operational$n[["stage_1"]]
  m <- x$operational$n[["stage_2"]]
  cv <- cluster_cv((1 - ab$p) / ab$p, 0.05, a, m)
  expect_equal(s$indicators$achieved, cv, tolerance = 1e-10)
  expect_identical(s$worst$name, "a")
  expect_identical(s$binding, "b")
  expect_equal(s$design$continuous[s$design$row == "cv"],
               max(x$detail$.cv_achieved))
  expect_null(s$field)
  out <- capture.output(print(s))
  expect_true(sprintf("worst cv: %.4f (a)", cv[1L]) %in% out)
  expect_true("binding: b (precision)" %in% out)
  expect_true("short of target: b" %in% out)
  expect_false(any(grepl("Design effect", out)))
})

test_that("moe targets read moe at the field design", {
  ci <- data.frame(name = c("stunting", "anemia"), p = c(0.30, 0.10),
                   moe = c(0.05, 0.03), icc_1 = c(0.02, 0.05))
  x <- n_cluster(ci, stage_cost = c(500, 50), budget = 1e5)
  s <- summary(x)
  a <- x$operational$n[["stage_1"]]
  m <- x$operational$n[["stage_2"]]
  cv <- cluster_cv((1 - ci$p) / ci$p, ci$icc_1, a, m)
  expect_identical(s$indicators$target_scale, c("moe", "moe"))
  expect_equal(s$indicators$achieved, z975 * ci$p * cv, tolerance = 1e-10)
  expect_equal(s$indicators$cv, x$operational$cv_by_target,
               tolerance = 1e-12)
})

test_that("cluster domains are read at each domain's whole design", {
  cd <- data.frame(name = c("st", "an", "st", "an"),
                   region = c("N", "N", "S", "S"),
                   p = c(0.30, 0.10, 0.32, 0.12), cv = c(0.08, 0.15),
                   icc_1 = 0.05)
  for (fit in list(
    n_cluster(cd, stage_cost = c(500, 50), domains = "region"),
    n_cluster(cd, stage_cost = c(500, 50), domains = "region",
              budget = 2e5, allocation = "joint")
  )) {
    s <- summary(fit)
    od <- fit$operational$domains
    row <- match(cd$region, od$region)
    cv <- cluster_cv((1 - cd$p) / cd$p, 0.05, od$stage_1[row],
                     od$stage_2[row])
    expect_equal(s$indicators$achieved, cv, tolerance = 1e-10)
    expect_identical(s$domains$binding, od$.binding)
    expect_identical(nrow(fit$detail), 4L)
    for (r in unique(cd$region)) {
      in_r <- s$indicators$domain == r
      ratio <- s$indicators$cv[in_r] / s$indicators$cv_target[in_r]
      expect_identical(which(s$indicators$binding[in_r]), which.max(ratio))
    }
  }
})

## Width over a pinned coverage set

test_that("every new print and summary fits 80 columns", {
  long <- paste(rep("indicator", 5), collapse = "_")
  region <- paste(rep("region", 4), collapse = "_")
  cases <- list(
    scalar = n_prop(p = 0.28, cv = 0.05, deff = 2.5, resp_rate = 0.9,
                    eligible_per_unit = 0.41),
    ratio = n_ratio(r = 420, cv_num = 1.2, cv_den = 0.45,
                    component_cor = 0.65, rmoe = 0.05, deff = 1.3,
                    resp_rate = 0.85),
    mixed = n_multi(data.frame(
      name = c(long, "spend", "floor"), p = c(0.3, NA, NA),
      var = c(NA, 250000, NA), mean = c(NA, 1200, NA),
      moe = c(0.05, NA, NA), cv = c(NA, 0.05, NA),
      min_cases = c(30, NA, NA), min_respondents = c(NA, NA, 900),
      resp_rate = c(0.9, 0.8, 0.7), deff = c(1.5, 2, 1)
    )),
    tie = n_multi(data.frame(name = c("x", "y"), p = 0.3, moe = 0.05)),
    domains = n_multi(data.frame(name = long, region = c(region, "b"),
                                 p = c(0.3, 0.2), moe = 0.05),
                      domains = "region"),
    cluster3 = n_cluster(stage_cost = c(500, 100, 50), icc = c(0.01, 0.05),
                         budget = 5e5, stage_n = c(50, NA, NA),
                         fixed_cost = 5000),
    cluster_ind3 = n_cluster(
      data.frame(name = c(long, "anemia"), p = c(0.3, 0.1),
                              moe = c(0.05, 0.03), icc_1 = c(0.02, 0.05),
                              icc_2 = c(0.05, 0.1)), stage_cost = c(500, 100, 50), budget = 2e5, fpc = "stage",
      N = c(300, 3000, 60000)),
    cluster_joint = n_cluster(
      data.frame(name = c(long, "an"), region = rep(c(region, "b"), each = 2),
                              p = c(0.3, 0.1), cv = c(0.08, 0.15),
                              icc_1 = 0.05), stage_cost = c(500, 50), domains = "region", budget = 2e5,
      allocation = "joint")
  )
  seen <- character(0)
  for (nm in names(cases)) {
    outputs <- list(capture.output(print(summary(cases[[nm]]))))
    if (!nm %in% c("scalar", "ratio")) {
      outputs <- c(outputs, list(capture.output(print(cases[[nm]]))))
    }
    for (out in outputs) {
      expect_lte(max(nchar(out, type = "width")), 80L, label = nm)
    }
    seen <- c(seen, nm)
  }
  expect_setequal(seen, names(cases))
  tie <- summary(cases$tie)
  expect_identical(sum(tie$indicators$binding), 1L)
})

## Prints that gained a summary

test_that("prints carry field counts, field cvs and the summary pointer", {
  x <- n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8, N = 10000)
  out <- capture.output(print(x))
  expect_true(any(grepl("n = 578 gross (net: 463)", out, fixed = TRUE)))
  expect_true("expected cases = 138.7" %in% out)
  expect_match(out[length(out)], "^# summary\\(\\) for ")

  m <- n_multi(data.frame(name = c("stunting", "vaccination", "anemia"),
                          p = c(0.30, 0.70, 0.10), moe = c(0.05, 0.05, 0.03),
                          deff = c(2, 1.5, 2.5)), resp_rate = 0.9)
  out <- capture.output(print(m))
  anemia <- strsplit(trimws(grep("^ anemia", out, value = TRUE)), " +")[[1]]
  expect_identical(anemia[4L], sprintf("%.4f", m$detail$.cv_field[3L]))
  expect_identical(anemia[4L], "0.1530")
  expect_identical(out[length(out) - 1L], "")

  for (fit in list(
    n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05),
    n_cluster(data.frame(name = c("a", "b"), p = c(0.3, 0.1),
                                      moe = c(0.05, 0.03), icc_1 = 0.05), stage_cost = c(500, 50), budget = 1e5)
  )) {
    out <- capture.output(print(fit))
    expect_false(any(grepl("^continuous optimum:", out)))
    expect_match(out[length(out)], "^# summary\\(\\) for ")
  }
})

## Field binding where rounding moves it

flip <- data.frame(name = c("a", "b"), p = c(0.33, 0.20),
                   cv = c(0.110, 0.111), icc_1 = c(0.099, 0.019))

test_that("the summary names the field binding, not the continuous one", {
  x <- n_cluster(flip, stage_cost = c(500, 50))
  expect_identical(x$binding, "a")
  expect_identical(x$operational$binding, "b")
  s <- summary(x)
  expect_identical(s$binding, "b")
  expect_identical(s$indicators$binding, c(FALSE, TRUE))
  expect_true("binding: b (precision)" %in% capture.output(print(s)))
})

test_that("each cluster domain keeps its continuous and its field binding", {
  both <- rbind(transform(flip, region = "N"), transform(flip, region = "S"))
  x <- n_cluster(both, stage_cost = c(500, 50),
                 domains = "region")
  expect_identical(x$detail$name[x$detail$.binding], c("a", "a"))
  expect_identical(x$domains$.binding, c("a", "a"))
  expect_identical(x$operational$domains$.binding, c("b", "b"))
})

test_that("stated targets follow their rows when domains interleave", {
  rows <- data.frame(
    name = c("a", "a", "b", "b"),
    region = c("N", "S", "N", "S"),
    p = c(0.30, 0.35, 0.10, 0.12),
    rmoe = c(0.20, NA, NA, 0.25),
    moe = c(NA, 0.05, 0.03, NA)
  )
  x <- n_multi(rows, domains = "region")
  key <- paste(x$detail$region, x$detail$name)
  expect_identical(key, c("N a", "N b", "S a", "S b"))
  input <- match(key, paste(rows$region, rows$name))
  expect_identical(x$detail$.target_scale,
                   ifelse(is.na(rows$rmoe), "moe", "rmoe")[input])
  expect_equal(x$detail$.target,
               ifelse(is.na(rows$rmoe), rows$moe, rows$rmoe)[input])
})

test_that("the field precision carries the eligibility rate of the fit", {
  x <- n_prop(p = 0.22, moe = 0.04, deff = 1.5, resp_rate = 0.87,
              eligibility_rate = 0.96)
  s <- summary(x)
  n_resp <- ceiling(x$n) * 0.87 * 0.96
  expect_equal(s$sizes[["respondents"]], n_resp)
  expect_equal(s$precision[["moe"]], wald_moe(0.22, n_resp, deff = 1.5),
               tolerance = 1e-10)
  expect_true("  eligibility rate: 0.96" %in% capture.output(print(s)))
})

## Unnamed indicators across interleaved domains

unnamed <- data.frame(
  p = c(0.2, 0.5, 0.3, 0.6), moe = c(0.03, 0.06, 0.04, 0.05),
  region = c("a", "b", "a", "b"),
  resp_rate = c(0.9, 0.8, 0.7, 0.6), deff = c(1, 1.2, 1.4, 1.6)
)

test_that("an unnamed indicator is labeled by its input row in every table", {
  x <- n_multi(unnamed, domains = "region")
  expect_identical(x$detail$region, c("a", "a", "b", "b"))
  expect_identical(as.integer(x$detail$name), c(1L, 3L, 2L, 4L))
  prec <- prec_multi(x)
  expect_setequal(paste(prec$detail$region, prec$detail$name),
                  paste(x$detail$region, x$detail$name))
  expect_identical(x$domains$.binding, x$detail$name[x$detail$.binding])
  expect_null(x$indicators$name)
})

test_that("an unnamed indicator's precision and rates are its own row's", {
  x <- n_multi(unnamed, domains = "region")
  s <- summary(x)$indicators
  row <- as.integer(s$name)
  quota <- z975^2 * unnamed$p * (1 - unnamed$p) * unnamed$deff /
    (unnamed$moe^2 * unnamed$resp_rate)
  field <- as.vector(ceiling(tapply(quota, unnamed$region, max))[
    unnamed$region[row]
  ])
  n_resp <- field * unnamed$resp_rate[row]
  expect_equal(s$field, field)
  expect_equal(s$respondents, n_resp)
  expect_equal(s$deff, unnamed$deff[row])
  expect_equal(s$achieved,
               wald_moe(unnamed$p[row], n_resp, deff = unnamed$deff[row]),
               tolerance = 1e-10)
  expect_false(anyNA(s$achieved))
})

test_that("unnamed cluster indicators are read at their own domain's design", {
  cd <- data.frame(p = c(0.2, 0.5, 0.3, 0.6), cv = c(0.08, 0.1, 0.09, 0.07),
                   region = c("a", "b", "a", "b"),
                   icc_1 = c(0.02, 0.05, 0.08, 0.03))
  seen <- character(0)
  for (route in c("separate", "joint")) {
    fit <- if (route == "joint") {
      n_cluster(cd, stage_cost = c(500, 50), domains = "region",
                budget = 2e5, allocation = "joint")
    } else {
      n_cluster(cd, stage_cost = c(500, 50), domains = "region")
    }
    s <- summary(fit)$indicators
    row <- as.integer(s$name)
    expect_identical(cd$region[row], s$domain, label = route)
    od <- fit$operational$domains
    at <- match(cd$region[row], od$region)
    cv <- cluster_cv((1 - cd$p[row]) / cd$p[row], cd$icc_1[row],
                     od$stage_1[at], od$stage_2[at])
    expect_equal(s$achieved, cv, tolerance = 1e-10, label = route)
    seen <- c(seen, route)
  }
  expect_setequal(seen, c("separate", "joint"))
})

test_that("domain values holding the label separator stay apart", {
  rows <- data.frame(name = "a", r1 = c("x/y", "x"), r2 = c("z", "y/z"),
                     p = c(0.3, 0.2), moe = 0.05)
  x <- n_multi(rows, domains = c("r1", "r2"))
  s <- summary(x)$indicators
  quota <- z975^2 * rows$p * (1 - rows$p) / rows$moe^2
  expect_equal(s$field, ceiling(quota))
  expect_equal(s$achieved, wald_moe(rows$p, ceiling(quota)),
               tolerance = 1e-10)
})

## Missing indicator names

test_that("a name column with a missing or empty name is refused", {
  base <- data.frame(p = c(0.3, 0.1), moe = c(0.05, 0.03))
  cl <- data.frame(p = c(0.33, 0.20), cv = c(0.110, 0.111),
                   icc_1 = c(0.099, 0.019))
  seen <- character(0)
  for (nm in list(c(NA, NA), c("a", NA), c("a", ""), c("a", " "))) {
    bad <- transform(base, name = nm)
    expect_error(n_multi(bad), "indicator name missing on row\\(s\\) ")
    expect_error(prec_multi(transform(bad, n = 400)),
                 "indicator name missing on row\\(s\\) ")
    expect_error(n_multi(rbind(transform(bad, region = "N"),
                               transform(bad, region = "S")),
                         domains = "region"),
                 "indicator name missing on row\\(s\\) ")
    expect_error(n_cluster(transform(cl, name = nm), stage_cost = c(500, 50)),
                 "indicator name missing on row\\(s\\) ")
    seen <- c(seen, paste(nm, collapse = "|"))
  }
  expect_length(seen, 4L)
  expect_error(n_multi(transform(base, name = c("a", NA))),
               "row\\(s\\) 2\\. Give each row a 'name', or drop the column")
})

test_that("the precision join refuses rows it cannot pair one to one", {
  detail <- data.frame(name = c("a", "b"), .target_scale = "moe")
  prec <- data.frame(name = c("a", "a"), .cv = c(0.1, 0.2))
  expect_error(.summary_join_prec(detail, prec, character(0)),
               "internal error")
  expect_error(.summary_join_prec(detail, prec[1L, ], character(0)),
               "internal error")
  floor_only <- data.frame(name = c("a", "b"), .target_scale = c("moe", NA))
  expect_true(is.na(.summary_join_prec(floor_only, prec[1L, ],
                                       character(0))$.cv[2L]))
})

## Headline binding of cluster domains

test_that("the headline binding of cluster domains is the field one", {
  both <- rbind(transform(flip, region = "N"), transform(flip, region = "S"))
  seen <- character(0)
  for (mode in c("cv", "budget")) {
    x <- if (mode == "budget") {
      n_cluster(both, stage_cost = c(500, 50),
                domains = "region", budget = 2e5)
    } else {
      n_cluster(both, stage_cost = c(500, 50),
                domains = "region")
    }
    expect_identical(x$binding, "a", label = mode)
    od <- x$operational$domains
    s <- summary(x)
    expect_identical(s$binding, od$.binding[which.max(od$.cv)], label = mode)
    expect_identical(s$binding, "b", label = mode)
    expect_identical(s$binding_cause, "precision", label = mode)
    expect_true("binding: b (precision)" %in% capture.output(print(s)))
    seen <- c(seen, mode)
  }
  expect_setequal(seen, c("cv", "budget"))
})

test_that("only the joint route reports a domain set by its floor", {
  loose <- rbind(transform(flip, region = "N"),
                 transform(flip, cv = c(0.3, 0.3), region = "S"))
  expect_warning(
    sep <- n_cluster(loose, stage_cost = c(500, 50),
                     domains = "region", min_n_domain = 300),
    "below min_n_domain = 300"
  )
  expect_lt(sep$domains$.total_n[2L], 300)
  expect_identical(sep$binding_cause, "precision")
  expect_identical(summary(sep)$binding_cause, "precision")

  joint <- n_cluster(loose, stage_cost = c(500, 50),
                     domains = "region", min_n_domain = 600, budget = 2e5,
                     allocation = "joint")
  expect_equal(joint$domains$.total_n[2L], 600, tolerance = 1e-8)
  expect_identical(joint$binding_cause, "min_n_domain")
  expect_identical(summary(joint)$binding_cause, "min_n_domain")

  free <- n_cluster(loose, stage_cost = c(500, 50),
                    domains = "region", budget = 2e5, allocation = "joint")
  expect_identical(free$binding_cause, "precision")
})

## Domain labels

test_that("domains whose pasted values collide print apart", {
  rows <- data.frame(name = "a", r1 = c("x/y", "x"), r2 = c("z", "y/z"),
                     p = c(0.3, 0.2), moe = 0.05, cv = 0.2, icc_1 = 0.05)
  fits <- list(
    multi = n_multi(rows[c("name", "r1", "r2", "p", "moe")],
                    domains = c("r1", "r2")),
    cluster = n_cluster(rows[c("name", "r1", "r2", "p", "cv",
                                            "icc_1")], stage_cost = c(500, 50),
                        domains = c("r1", "r2"))
  )
  seen <- character(0)
  for (nm in names(fits)) {
    s <- summary(fits[[nm]])
    expect_identical(s$domains$domain, c("r1=x/y, r2=z", "r1=x, r2=y/z"),
                     label = nm)
    expect_identical(s$domain_values,
                     data.frame(r1 = c("x/y", "x"), r2 = c("z", "y/z")),
                     label = nm)
    out <- capture.output(print(s))
    heads <- grep("^Indicators in ", out, value = TRUE)
    expect_identical(heads, c("Indicators in r1=x/y, r2=z",
                              "Indicators in r1=x, r2=y/z"), label = nm)
    expect_identical(sum(grepl("^ a ", out)), 2L, label = nm)
    expect_true(any(grepl("^ r1 +r2 ", out)), label = nm)
    expect_lte(max(nchar(out, type = "width")), 80L)
    seen <- c(seen, nm)
  }
  expect_setequal(seen, c("multi", "cluster"))
})

test_that("a domain label is distinct for distinct domains", {
  naive <- data.frame(r1 = c("a, r2=b", "a"), r2 = c("c", "b, r2=c"))
  expect_identical(paste0("r1=", naive$r1, ", r2=", naive$r2)[1L],
                   paste0("r1=", naive$r1, ", r2=", naive$r2)[2L])
  lab <- .summary_domain_label(naive, c("r1", "r2"))
  expect_identical(lab, c("r1=\"a, r2=b\", r2=c", "r1=a, r2=\"b, r2=c\""))
  quoted <- data.frame(r1 = c("a\"b", "a\\b"), r2 = "c")
  expect_identical(.summary_domain_label(quoted, c("r1", "r2")),
                   c("r1=\"a\\\"b\", r2=c", "r1=\"a\\\\b\", r2=c"))
  expect_identical(.summary_domain_label(data.frame(region = c("N", "a/b")),
                                         "region"),
                   c("N", "a/b"))
})

test_that("long domain values keep their sections apart", {
  rows <- data.frame(name = "a", p = c(0.3, 0.2), moe = 0.05,
                     region = c("Northern highlands A", "Northern highlands B"),
                     area = "urban")
  out <- capture.output(print(summary(n_multi(rows,
                                              domains = c("region", "area")))))
  heads <- grep("^Indicators in ", out, value = TRUE)
  expect_identical(heads, paste0("Indicators in region=", rows$region,
                                 ", area=urban"))
  expect_lte(max(nchar(out, type = "width")), 80L)
})

## Design quantities of a multi-indicator cluster summary

test_that("per-row design quantities follow the inputs", {
  x <- n_cluster(flip, stage_cost = c(500, 50),
                 resp_rate = 0.9, eligibility_rate = 0.8)
  s <- summary(x)
  ind <- s$indicators
  expect_equal(ind$icc_1, flip$icc_1)
  expect_equal(ind$unit_relvar, (1 - flip$p) / flip$p, tolerance = 1e-12)
  expect_equal(ind$resp_rate, c(0.9, 0.9))
  expect_false(any(c("icc_2", "var_ratio_2", "resp_rate_2",
                     "alpha") %in% names(ind)))
  a <- s$assumptions
  expect_null(a$icc_1)
  expect_null(a$unit_relvar)
  expect_null(a$resp_rate_1)
  expect_identical(a$var_ratio_1, 1)
  expect_identical(a$resp_rate, 0.9)
  expect_identical(a$eligibility_rate, 0.8)
  expect_identical(a$prop_method, "wald")
  out <- capture.output(print(s))
  at <- match("Assumptions by indicator", out)
  expect_false(is.na(at))
  expect_true(any(grepl("^ Indicator +icc_1 +unit_relvar$", out)))
  expect_true(any(grepl("^ a +0\\.099 +2\\.03$", out)))
  expect_true("  var_ratio_1: 1" %in% out)
  expect_true("  eligibility rate: 0.8" %in% out)
})

test_that("a three-stage table of varying quantities wraps within 80 columns", {
  rows <- data.frame(
    name = c("stunting", "anemia", "vacc"), p = c(0.3, 0.1, 0.7),
    moe = c(0.05, 0.03, 0.05), icc_1 = c(0.05, 0.02, 0.1),
    icc_2 = c(0.1, 0.2, 0.05), var_ratio_1 = c(1, 1.1, 0.9),
    resp_rate_1 = c(0.95, 0.9, 1), resp_rate_2 = c(0.9, 0.95, 1),
    resp_rate = c(0.9, 0.8, 0.85), alpha = c(0.05, 0.1, 0.05)
  )
  s <- summary(n_cluster(rows, stage_cost = c(400, 80, 10)))
  ind <- s$indicators
  for (col in c("icc_1", "icc_2", "var_ratio_1", "resp_rate_1",
                "resp_rate_2", "resp_rate", "alpha")) {
    expect_equal(ind[[col]], rows[[col]], label = col)
  }
  expect_equal(ind$var_ratio_2, rows$var_ratio_1 * (1 - rows$icc_1),
               tolerance = 1e-12)
  out <- capture.output(print(s))
  expect_lte(max(nchar(out, type = "width")), 80L)
  tail_out <- out[seq(match("Assumptions by indicator", out), length(out))]
  expect_identical(sum(grepl("^ Indicator ", tail_out)), 2L)
  shown <- unlist(strsplit(trimws(grep("^ Indicator ", tail_out,
                                       value = TRUE)), " +"))
  expect_setequal(setdiff(shown, "Indicator"),
                  c("icc_1", "icc_2", "unit_relvar", "var_ratio_1",
                    "var_ratio_2", "resp_rate_1", "resp_rate_2",
                    "resp_rate", "alpha"))
})

test_that("design quantities follow their rows across domains", {
  rows <- data.frame(name = c("a", "b", "a", "b"),
                     region = c("N", "N", "S", "S"),
                     p = c(0.3, 0.2, 0.35, 0.25), cv = 0.12,
                     icc_1 = c(0.05, 0.02, 0.08, 0.03))
  s <- summary(n_cluster(rows, stage_cost = c(500, 50),
                         domains = "region"))
  ind <- s$indicators
  at <- match(paste(ind$domain, ind$name), paste(rows$region, rows$name))
  expect_equal(ind$icc_1, rows$icc_1[at])
  expect_equal(ind$unit_relvar, (1 - rows$p[at]) / rows$p[at],
               tolerance = 1e-12)
  out <- capture.output(print(s))
  expect_true(any(grepl("^ Domain +Indicator +icc_1 +unit_relvar$", out)))
  expect_true(any(grepl("^ S +a +0\\.08 ", out)))
})

test_that("design quantities of unnamed rows are read from their input rows", {
  rows <- data.frame(p = c(0.3, 0.1), cv = c(0.1, 0.1),
                     icc_1 = c(0.02, 0.05), resp_rate = c(0.8, 0.9))
  ind <- summary(n_cluster(rows, stage_cost = c(500, 50)))$indicators
  expect_identical(ind$name, c("1", "2"))
  expect_equal(ind$icc_1, rows$icc_1)
  expect_equal(ind$unit_relvar, (1 - rows$p) / rows$p, tolerance = 1e-12)
  expect_equal(ind$resp_rate, rows$resp_rate)
  expect_equal(ind$var_ratio_1, c(1, 1))
})

test_that("design quantities of unnamed rows follow them across interleaved domains", {
  rows <- data.frame(p = c(0.3, 0.2, 0.1, 0.4), cv = 0.1,
                     icc_1 = c(0.02, 0.03, 0.04, 0.05),
                     resp_rate = c(0.7, 0.8, 0.9, 0.6),
                     region = c("A", "B", "A", "B"))
  ind <- summary(n_cluster(rows, stage_cost = c(500, 50),
                           domains = "region"))$indicators
  at <- as.integer(ind$name)
  expect_identical(ind$domain, rows$region[at])
  expect_identical(at, c(1L, 3L, 2L, 4L))
  expect_equal(ind$icc_1, rows$icc_1[at])
  expect_equal(ind$unit_relvar, (1 - rows$p[at]) / rows$p[at],
               tolerance = 1e-12)
  expect_equal(ind$resp_rate, rows$resp_rate[at])
})

test_that("third-stage quantities of unnamed rows are read from their input rows", {
  rows <- data.frame(p = c(0.3, 0.1), cv = 0.1, icc_1 = c(0.02, 0.05),
                     icc_2 = c(0.1, 0.2), resp_rate_2 = c(0.9, 0.8))
  ind <- summary(n_cluster(rows, stage_cost = c(500, 100, 20)))$indicators
  expect_equal(ind$icc_1, rows$icc_1)
  expect_equal(ind$icc_2, rows$icc_2)
  expect_equal(ind$resp_rate_2, rows$resp_rate_2)
  expect_equal(ind$resp_rate, c(1, 1))
})

test_that("summary() produces structured summaries for n_change() and n_pooled()", {
  ch <- n_change(p = c(0.30, 0.35), moe = 0.03, overlap = 0.5,
                 overlap_cor = 0.6, resp_rate = 0.85)
  sch <- summary(ch)
  expect_s3_class(sch, "summary.svyplan_n")
  expect_identical(sch$kind, "change")
  expect_equal(sch$n$field, ceiling(ch$n))
  expect_lte(sch$precision[["moe"]], 0.03 + 1e-9)
  out_ch <- capture.output(print(sch))
  expect_true(any(grepl("change between two occasions", out_ch)))
  expect_true(any(grepl("Precision at the field size", out_ch)))

  pl <- n_pooled(p = 0.30, moe = 0.02, occasions = 4, overlap = 0.75,
                 cor_decay = 0.6, resp_rate = 0.9)
  spl <- summary(pl)
  expect_s3_class(spl, "summary.svyplan_n")
  expect_identical(spl$kind, "pooled")
  expect_identical(spl$occasions, 4L)
  expect_equal(spl$sizes[["total_field"]], ceiling(pl$n) * 4)
  expect_lte(spl$precision[["moe"]], 0.02 + 1e-9)
  out_pl <- capture.output(print(spl))
  expect_true(any(grepl("pooled estimate \\(4 occasions\\)", out_pl)))
  expect_true(any(grepl("Total field", out_pl)))
})

## Multi-indicator summaries carry the yields

test_that("multi-indicator summaries show eligibility and eligible units", {
  ind <- data.frame(name = c("a", "b"), p = c(0.3, 0.4), cv = 0.1)
  s <- capture.output(summary(n_multi(ind, eligibility_rate = 0.8,
                                      eligible_per_unit = 2)))
  expect_true(any(grepl("eligibility rate: 0.8", s)))
  expect_true(any(grepl("eligible per unit: 2", s)))
  cl <- n_cluster(transform(ind, icc_1 = 0.05), stage_cost = c(100, 10),
                  eligible_per_unit = 2)
  expect_true(any(grepl("eligible per unit: 2", capture.output(summary(cl)))))
})
