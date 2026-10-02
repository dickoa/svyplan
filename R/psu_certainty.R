#' Validate a PSU register against the frame it refines
#'
#' The register states the size of every PSU in a stratum, where `N_psu`
#' states only how many there are. The two are alternatives, and supplying
#' both leaves it ambiguous which the variance model should read.
#' @keywords internal
#' @noRd
.check_psu_table <- function(psu, frame, measures) {
  if (!is.data.frame(psu) || nrow(psu) == 0L) {
    stop("'psu' must be a data frame with at least one row", call. = FALSE)
  }
  missing_cols <- setdiff(c("stratum", "N"), names(psu))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("'psu' must contain: %s",
              paste(sQuote(missing_cols), collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(psu$N) || anyNA(psu$N) || any(!is.finite(psu$N)) ||
        any(psu$N <= 0)) {
    stop("'psu$N' must be positive and finite for every PSU", call. = FALSE)
  }
  if ("N_psu" %in% names(frame)) {
    stop(
      "'psu' and the frame column 'N_psu' are alternatives: 'psu' gives the size of every PSU and 'N_psu' only how many there are",
      call. = FALSE
    )
  }
  if (!"n_per_psu" %in% names(frame)) {
    stop("'psu' requires the frame column 'n_per_psu', the within-PSU take",
         call. = FALSE)
  }
  if (!"icc_psu" %in% names(measures)) {
    stop(
      "'psu' requires 'measures$icc_psu', the within-PSU homogeneity of the noncertainty part",
      call. = FALSE
    )
  }
  has_cost <- c("cost_psu", "cost_ssu") %in% names(frame)
  if (any(has_cost) && !all(has_cost)) {
    stop("stage costs come as a pair: supply both 'cost_psu' and 'cost_ssu'",
         call. = FALSE)
  }
  if (all(has_cost) && "unit_cost" %in% names(frame)) {
    stop(
      "certainty-aware allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }
  three_stage <- intersect(
    c("n_per_ssu", "icc_ssu", "var_ratio_ssu", "cost_tsu", "N_ssu"),
    c(names(frame), names(measures))
  )
  if (length(three_stage) > 0L) {
    stop(
      sprintf("certainty-aware allocation is two-stage; remove: %s",
              paste(sQuote(three_stage), collapse = ", ")),
      call. = FALSE
    )
  }
  if ("take_all" %in% names(frame) && any(.check_take_all(frame$take_all, nrow(frame)))) {
    stop(
      "'take_all' is not supported with 'psu' because taking every PSU does not imply an ultimate-unit census: the within-PSU take still applies",
      call. = FALSE
    )
  }
  stratum <- as.character(frame$stratum %||% seq_len(nrow(frame)))
  psu_stratum <- as.character(psu$stratum)
  unknown <- setdiff(unique(psu_stratum), stratum)
  if (length(unknown) > 0L) {
    stop(
      sprintf("'psu' stratum not in 'frame': %s",
              paste(sQuote(unknown), collapse = ", ")),
      call. = FALSE
    )
  }
  empty <- setdiff(stratum, unique(psu_stratum))
  if (length(empty) > 0L) {
    stop(
      sprintf("'frame' stratum with no PSU in 'psu': %s",
              paste(sQuote(empty), collapse = ", ")),
      call. = FALSE
    )
  }
  psu_total <- vapply(stratum, function(h) sum(psu$N[psu_stratum == h]),
                      numeric(1))
  frame_total <- as.numeric(frame$N)
  # Relative, never floored: a derived sum and a supplied total agree only to
  # the scale of the terms, and that scale is platform-dependent.
  off <- abs(psu_total - frame_total) > 1e-8 * abs(frame_total)
  if (any(off)) {
    bad <- stratum[off]
    stop(
      sprintf(
        "the sum of 'psu$N' must equal 'frame$N' in stratum: %s",
        paste(sQuote(bad), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if ("certainty" %in% names(psu)) {
    psu$certainty <- as.logical(.check_take_all(psu$certainty, nrow(psu)))
  }
  psu$stratum <- psu_stratum
  psu
}

#' Certainty classification of one stratum's PSUs
#'
#' A PSU is certainty when a proportional-to-size selection would give it an
#' inclusion probability of at least one, which reduces to a size at or above
#' the take over the sampling fraction. A supplied flag adds to that and never
#' removes from it: a PSU above the threshold cannot be sampled less often
#' than certainly, whatever the caller asks for.
#' @keywords internal
#' @noRd
.psu_classify <- function(sizes, threshold, forced = NULL) {
  certain <- sizes >= threshold
  source <- ifelse(certain, "threshold", NA_character_)
  if (!is.null(forced)) {
    added <- forced & !certain
    certain <- certain | forced
    source[added] <- "supplied"
  }
  list(certain = certain, source = source)
}

#' Anticipated design effect of a stratum split into certainty and remainder
#'
#' The variance of a design treating the two parts as strata, over the
#' variance of a simple random sample of the same size. The certainty part has
#' no first-stage sampling. The form is deliberately with-replacement: the
#' finite population correction is applied in the solver's own `B` term, so a
#' design effect carrying one would count it twice.
#'
#' The certainty part gets its proportional share of the allocation, which
#' reduces the design effect to
#' `(size_certain + size_rest * (1 + icc * (take - 1))) / N_h`, free of the
#' allocation. Rounding the share up made the settle step a discontinuous
#' map that could alternate without converging.
#' @keywords internal
#' @noRd
.psu_deff <- function(n_h, N_h, size_certain, size_rest, icc, take) {
  f <- n_h / N_h
  n_certain <- f * size_certain
  n_rest <- n_h - n_certain
  term_certain <- if (size_certain > 0 && n_certain > 0) {
    size_certain^2 / n_certain
  } else {
    0
  }
  term_rest <- if (size_rest > 0 && n_rest > 0) {
    size_rest^2 / n_rest * (1 + icc * (take - 1))
  } else {
    0
  }
  max((n_h / N_h^2) * (term_certain + term_rest), 1)
}

#' Noncertainty ICC per measures row, zero for a single-PSU stratum
#'
#' A stratum with one PSU has no between-PSU variance, so its ICC is zero by
#' definition. The supplied value is replaced rather than honoured and may be
#' missing.
#' @keywords internal
#' @noRd
.psu_icc <- function(measures, stratum, idx_of) {
  icc <- measures$icc_psu
  if (is.logical(icc) && all(is.na(icc))) icc <- as.numeric(icc)
  if (!is.numeric(icc)) {
    stop("'measures$icc_psu' must contain values in [0, 1]", call. = FALSE)
  }
  row_h <- match(as.character(measures$stratum), stratum)
  single <- !is.na(row_h) & lengths(idx_of)[row_h] == 1L
  icc[single] <- 0
  missing <- is.na(icc)
  if (any(missing)) {
    stop(
      sprintf(
        "'measures$icc_psu' is missing in stratum %s. Only a stratum with a single PSU may leave it missing",
        paste(sQuote(unique(as.character(measures$stratum[missing]))),
              collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(icc < 0 | icc > 1)) {
    stop("'measures$icc_psu' must contain values in [0, 1]", call. = FALSE)
  }
  icc
}

#' Certainty-aware joint allocation from a PSU register
#'
#' The classification and the allocation each determine the other, so the
#' solve is a loop. It iterates until the classification orbit repeats, takes
#' the feasible member of that orbit, then holds the classification and settles
#' the allocation under it. The design effect depends on the classification
#' only, so settling takes one confirming solve unless the settled allocation
#' puts a PSU above the threshold and the classification grows.
#' @keywords internal
#' @noRd
.n_alloc_psu <- function(
  frame,
  psu,
  measures,
  targets,
  unit_cost,
  alpha,
  deff,
  resp_rate,
  min_n_stratum,
  objective,
  budget,
  df,
  max_iter = 30L,
  settle_iter = 200L,
  tolerance = 1e-9
) {
  psu <- .check_psu_table(psu, frame, measures)

  stratum <- as.character(frame$stratum %||% seq_len(nrow(frame)))
  H <- nrow(frame)
  N_h <- frame$N
  take <- frame$n_per_psu
  if (!is.numeric(take) || anyNA(take) || any(take < 1) ||
        any(abs(take - round(take)) > 1e-8)) {
    stop("'n_per_psu' must contain positive whole numbers", call. = FALSE)
  }
  take <- rep_len(as.numeric(round(take)), H)
  idx_of <- split(seq_len(nrow(psu)), factor(psu$stratum, levels = stratum))
  forced <- if ("certainty" %in% names(psu)) psu$certainty else NULL

  # Kept as supplied, so the object records the caller's input.
  unit_cost_in <- unit_cost
  stage_cost <- c("cost_psu", "cost_ssu") %in% names(frame)
  if (all(stage_cost) && !is.null(unit_cost)) {
    stop(
      "certainty-aware allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }
  if (all(stage_cost)) {
    cost_psu <- rep_len(as.numeric(frame$cost_psu), H)
    cost_ssu <- rep_len(as.numeric(frame$cost_ssu), H)
    if (any(!is.finite(cost_psu)) || any(cost_psu <= 0) ||
          any(!is.finite(cost_ssu)) || any(cost_ssu <= 0)) {
      stop("'cost_psu' and 'cost_ssu' must be positive and finite",
           call. = FALSE)
    }
    unit_cost <- cost_psu / take + cost_ssu
  } else {
    cost_psu <- NULL
    cost_ssu <- NULL
  }
  base_frame <- frame[
    , setdiff(names(frame), c("n_per_psu", "cost_psu", "cost_ssu")),
    drop = FALSE
  ]
  icc_row <- .psu_icc(measures, stratum, idx_of)
  base_measures <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  # A design effect the caller supplies describes a source this model does not
  # carry, so it multiplies the clustering rather than replacing it.
  user_deff <- if ("deff" %in% names(measures)) {
    ifelse(is.na(measures$deff), deff, measures$deff)
  } else {
    rep(deff, nrow(measures))
  }
  row_h <- match(as.character(measures$stratum), stratum)
  if (anyNA(row_h)) {
    stop("'measures$stratum' must match 'frame$stratum'", call. = FALSE)
  }

  deff_measures <- function(cls) {
    m <- base_measures
    m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
      h <- row_h[r]
      i <- idx_of[[h]]
      .psu_deff(
        cls$n_h[h], N_h[h], sum(psu$N[i][cls$certain[i]]),
        sum(psu$N[i][!cls$certain[i]]), icc_row[r], take[h]
      )
    }, numeric(1))
    m
  }

  solve_at <- function(cls) {
    m <- deff_measures(cls)
    # Although the inner solver receives an element-shaped frame, this is
    # a multistage model. Retain its p(1-p) working variance on every call.
    .n_alloc_bethel(
      frame = base_frame, measures = m, targets = targets,
      unit_cost = unit_cost, alpha = alpha, deff = deff,
      resp_rate = resp_rate, min_n_stratum = min_n_stratum,
      objective = objective, budget = budget, df = df,
      .finite_prop_var = FALSE
    )
  }

  classify_at <- function(n_h) {
    threshold <- take / (n_h / N_h)
    certain <- logical(nrow(psu))
    source <- rep(NA_character_, nrow(psu))
    for (h in seq_len(H)) {
      i <- idx_of[[h]]
      one <- .psu_classify(psu$N[i], threshold[h], forced[i])
      certain[i] <- one$certain
      source[i] <- one$source
    }
    list(n_h = n_h, certain = certain, source = source, threshold = threshold)
  }

  # A first pass with every PSU in the remainder, so full clustering is
  # charged, supplies the allocation the first threshold is read from.
  flat <- list(n_h = N_h, certain = rep(FALSE, nrow(psu)))
  state <- classify_at(solve_at(flat)$detail$n)

  seen <- character(0)
  orbit <- list()
  verdict <- "limit_reached"
  for (it in seq_len(max_iter)) {
    key <- paste(which(state$certain), collapse = ",")
    if (key %in% seen) {
      # Sliced from 'first', not first + 1: seq() counts down on a length-one
      # orbit rather than returning nothing.
      first <- match(key, seen)
      orbit <- orbit[seq.int(first, length(orbit))]
      verdict <- if (length(orbit) == 1L) "converged" else "cycle"
      break
    }
    seen <- c(seen, key)
    fit <- solve_at(state)
    state <- classify_at(fit$detail$n)
    orbit[[length(orbit) + 1L]] <- state
  }
  if (identical(verdict, "limit_reached")) orbit <- list(state)

  # A cycle always holds a feasible member, the most certain having been
  # sized under the largest design effect in the orbit.
  pick <- .psu_pick_feasible(
    orbit, base_frame, base_measures, targets, user_deff, row_h, idx_of,
    psu, N_h, icc_row, take, unit_cost, alpha, deff, resp_rate,
    min_n_stratum, objective, budget, df
  )
  chosen <- orbit[[pick$index]]

  # Absorbing and settling again terminates: the held set only grows.
  settled <- chosen
  settle_used <- 0L
  absorbed <- 0L
  repeat {
    settled_ok <- FALSE
    for (i in seq_len(settle_iter)) {
      fit <- solve_at(settled)
      settle_used <- settle_used + 1L
      if (max(abs(fit$detail$n - settled$n_h)) < tolerance) {
        settled_ok <- TRUE
        break
      }
      settled$n_h <- fit$detail$n
    }
    # The returned constraints would describe the previous allocation's
    # design effect, so an unsettled plan is refused rather than returned.
    if (!settled_ok) {
      stop(
        "internal error: the allocation did not settle under the held classification",
        call. = FALSE
      )
    }
    settled$n_h <- fit$detail$n
    implied <- classify_at(fit$detail$n)
    extra <- implied$certain & !settled$certain
    if (!any(extra)) break
    settled$certain <- settled$certain | implied$certain
    absorbed <- absorbed + sum(extra)
  }
  settled$threshold <- implied$threshold
  # Attribution reads the bare threshold test, not `implied$certain`, which
  # has the caller's flag already folded into it and would report every
  # flagged PSU as reaching the threshold on its own.
  above <- psu$N >= settled$threshold[match(psu$stratum, stratum)]
  settled$source <- ifelse(
    above, "threshold",
    ifelse(if (is.null(forced)) FALSE else forced, "supplied", "orbit")
  )
  settled$source[!settled$certain] <- NA_character_
  settled$above <- above

  # The whole-unit design. An element-level rounding of this allocation is
  # not fieldable: the remainder has to be a whole number of PSUs at the
  # stated take, and the certainty PSUs carry their own whole takes.
  assess <- function(n_try) {
    ev <- tryCatch(
      .prec_alloc_bethel(
        frame = base_frame, n = n_try,
        measures = deff_measures(list(n_h = n_try, certain = settled$certain)),
        targets = targets, unit_cost = unit_cost, alpha = alpha, deff = deff,
        resp_rate = resp_rate, min_n_stratum = min_n_stratum,
        objective = objective, budget = budget, df = df,
        .allow_fractional_stages = TRUE, .finite_prop_var = FALSE
      ),
      error = function(e) NULL
    )
    list(
      pass = !is.null(ev) &&
        !is.null(ev$detail$.pass) && all(ev$detail$.pass),
      constraints = if (is.null(ev)) NULL else ev$detail
    )
  }
  operational <- .psu_operational(
    settled$n_h, N_h, take, psu, idx_of, settled$certain, assess
  )
  operational$constraints <- assess(operational$n_int)$constraints

  list(
    fit = fit,
    frame = frame,
    measures = measures,
    operational = operational,
    take = take,
    unit_cost = unit_cost_in,
    cost_psu = cost_psu,
    cost_ssu = cost_ssu,
    psu = psu,
    certainty = settled$certain,
    source = settled$source,
    threshold = settled$threshold,
    implied = settled$above,
    verdict = verdict,
    orbit = length(orbit),
    iterations = length(seen),
    settle_iterations = settle_used,
    absorbed = absorbed,
    feasible_member = pick$feasible,
    idx_of = idx_of,
    stratum = stratum
  )
}

#' The whole-unit design a certainty plan is fielded as
#'
#' The continuous allocation is a number of ultimate units, but the design
#' draws them two ways: every certainty PSU contributes its own whole take at
#' the stratum rate, and the remainder is a whole number of PSUs at the stated
#' `n_per_psu`. A whole-unit total that is not the sum of those two is not a
#' design anyone can field, which is what an element-level rounding of this
#' allocation would report.
#'
#' Rounding up in both parts can only add sample, and adding sample under a
#' held classification can still leave a target short, because the design
#' effect moves with the allocation. The repair adds whole PSUs to the
#' stratum where one buys the most, in the same spirit as the solver's own
#' integer repair.
#' @keywords internal
#' @noRd
.psu_operational <- function(n_h, N_h, take, psu, idx_of, certain,
                             assess, max_repair = 200L) {
  H <- length(n_h)
  f <- n_h / N_h
  # A certainty PSU is above the threshold, so its own take at the stratum
  # rate is at least `take`; it is not capped at it.
  per_psu <- lapply(seq_len(H), function(h) {
    i <- idx_of[[h]][certain[idx_of[[h]]]]
    if (!length(i)) return(numeric(0))
    pmin(ceiling(f[h] * psu$N[i]), psu$N[i])
  })
  n_certain_int <- vapply(per_psu, sum, numeric(1))
  n_rest <- pmax(n_h - n_certain_int, 0)
  available <- vapply(seq_len(H), function(h) {
    sum(!certain[idx_of[[h]]])
  }, numeric(1))
  n_psu_draw <- pmin(ceiling(n_rest / take), available)

  short <- n_rest > n_psu_draw * take + 1e-8
  if (any(short)) {
    stop(
      sprintf(
        "the noncertainty part of stratum %s cannot supply its allocation: %d PSU(s) at a take of %d fall short",
        paste(sQuote(names(idx_of)[short]), collapse = ", "),
        max(available[short]), max(take[short])
      ),
      call. = FALSE
    )
  }

  repair <- 0L
  for (i in seq_len(max_repair)) {
    n_int <- n_certain_int + n_psu_draw * take
    ok <- assess(n_int)
    if (isTRUE(ok$pass)) break
    room <- n_psu_draw < available
    if (!any(room)) break
    # Add the PSU that buys the most on the worst constraint, which is the
    # stratum with the largest shortfall per unit added.
    gain <- ifelse(room, n_int / N_h, -Inf)
    n_psu_draw[which.max(gain)] <- n_psu_draw[which.max(gain)] + 1L
    repair <- repair + 1L
  }
  n_int <- n_certain_int + n_psu_draw * take

  list(
    n_int = n_int,
    n_certain_int = n_certain_int,
    n_psu_draw = n_psu_draw,
    per_psu = per_psu,
    repair = repair,
    pass = isTRUE(assess(n_int)$pass)
  )
}

#' Fold a certainty solve into the allocation object it refines
#'
#' The class does not change, so `predict()`, `summary()` and `design_df()`
#' keep working. samplyr refuses these fits as stage sizes, because no
#' per-stratum total carries a design of whole certainty takes plus a PPS
#' remainder. It fields the plan from `$psu` instead, which is why that table
#' carries the exact per-PSU takes. `$detail` gains the per-stratum counts and
#' the threshold, `$psu` is the per-PSU classification, and the convergence
#' verdict joins `$optimization`, which already reports how the solve went.
#' @keywords internal
#' @noRd
.psu_result <- function(x) {
  fit <- x$fit
  H <- nrow(fit$detail)
  n_certain <- vapply(seq_len(H), function(h) {
    sum(x$certainty[x$idx_of[[h]]])
  }, numeric(1))
  n_all <- vapply(x$idx_of, length, numeric(1))
  fit$detail$n_psu_certain <- n_certain
  fit$detail$n_psu_rest <- n_all - n_certain
  fit$detail$threshold <- x$threshold

  # Two-stage, not element-level: the solver's integerizer cannot see the
  # take, so its total is not a drawable design.
  op <- x$operational
  fit$detail$n_int <- op$n_int
  fit$detail$n_certain_int <- op$n_certain_int
  fit$detail$n_psu_draw <- op$n_psu_draw
  fit$operational$n <- sum(op$n_int)
  fit$operational$n_certain <- sum(op$n_certain_int)
  fit$operational$n_psu_draw <- sum(op$n_psu_draw)
  # Cost and the constraint table follow the same design as the count above
  # them; the solver's element-level ones describe a total this plan does not
  # field.
  if (is.null(x$cost_psu)) {
    fit$operational$cost <- sum(op$n_int * (fit$params$cost_h %||% 1))
  } else {
    # Under stage costs the design is priced as it is fielded: a visit to
    # every PSU that is entered, certainty or drawn, plus every interview.
    fit$operational$cost <- sum(
      (n_certain + op$n_psu_draw) * x$cost_psu + op$n_int * x$cost_ssu
    )
    # The continuous optimum is priced the same way, so the two readings on a
    # printed block differ by integerization and not by basis.
    rest_cont <- pmax(fit$detail$n - op$n_certain_int, 0)
    fit$params$achieved$cost <- sum(
      (n_certain + rest_cont / x$take) * x$cost_psu +
        fit$detail$n * x$cost_ssu
    )
  }
  if (!is.null(op$constraints)) fit$operational$constraints <- op$constraints
  fit$operational$all_pass <- op$pass
  fit$operational$repair_iterations <- op$repair
  if (!isTRUE(op$pass)) {
    warning(
      "the whole-unit design does not meet every target and the remainder has no PSU left to add",
      call. = FALSE
    )
  }

  # The take each PSU is fielded at, exposed from the operational design
  # rather than recomputed, so sum(n_take[certainty]) equals n_certain_int
  # exactly in every stratum and a consumer never re-derives the rate.
  n_take <- numeric(nrow(x$psu))
  for (h in seq_along(x$idx_of)) {
    i <- x$idx_of[[h]]
    n_take[i] <- x$take[h]
    n_take[i[x$certainty[i]]] <- x$operational$per_psu[[h]]
  }

  stratum_of <- x$psu$stratum
  fit$psu <- data.frame(
    stratum = stratum_of,
    N = x$psu$N,
    certainty = x$certainty,
    n_take = n_take,
    .certainty_source = x$source,
    .threshold = x$threshold[match(stratum_of, x$stratum)],
    stringsAsFactors = FALSE
  )
  if (!is.null(x$psu$psu_id)) {
    fit$psu <- cbind(
      data.frame(psu_id = x$psu$psu_id, stringsAsFactors = FALSE), fit$psu
    )
  }
  fit$psu$.distance <- fit$psu$N / fit$psu$.threshold - 1
  # One direction is unexecutable, the other a design a caller could have
  # written by hand with 'certainty'.
  if (any(x$implied & !x$certainty)) {
    stop(
      "internal error: the settled allocation implies a certainty PSU the held classification excludes",
      call. = FALSE
    )
  }

  # A supplied flag is the caller's design, so only a PSU held by the cycle
  # rule keeps the plan from being a fixed point.
  wanted <- x$implied
  if (!is.null(x$psu$certainty)) wanted <- wanted | x$psu$certainty
  fit$optimization$certainty <- list(
    verdict = x$verdict,
    orbit = x$orbit,
    iterations = x$iterations,
    settle_iterations = x$settle_iterations,
    absorbed = x$absorbed,
    feasible_member = x$feasible_member,
    fixed_point = identical(as.logical(wanted), as.logical(x$certainty))
  )
  # The caller's frame and measures, not the solver's stripped copies, which
  # would lose the take and the ICC.
  fit$params$frame <- x$frame
  fit$params$measures <- x$measures
  fit$params$unit_cost <- x$unit_cost
  fit$params$psu <- x$psu
  fit
}

#' Certainty-aware assessment of a supplied allocation
#'
#' No loop. The allocation is given, so the threshold it implies is given with
#' it, and the classification is a function of the answer being assessed
#' rather than something the answer has to agree with. Reading the
#' classification off the very allocation under assessment is what makes this
#' self-consistent by construction.
#'
#' A `certainty` column reproduces a held classification exactly, which is how
#' the round trip from [n_alloc()] closes: a fitted plan holds a classification
#' that contains every PSU above its threshold, so forcing that set and
#' re-deriving the rest returns the same split.
#' @keywords internal
#' @noRd
.prec_alloc_psu <- function(
  frame,
  psu,
  n,
  measures,
  targets,
  unit_cost,
  alpha,
  deff,
  resp_rate,
  min_n_stratum,
  objective,
  budget,
  df,
  .allow_fractional_stages = FALSE
) {
  if (is.null(n)) stop("'n' is required", call. = FALSE)
  psu <- .check_psu_table(psu, frame, measures)

  stratum <- as.character(frame$stratum %||% seq_len(nrow(frame)))
  H <- nrow(frame)
  N_h <- frame$N
  take <- frame$n_per_psu
  if (!is.numeric(take) || anyNA(take) || any(take < 1) ||
        any(abs(take - round(take)) > 1e-8)) {
    stop("'n_per_psu' must contain positive whole numbers", call. = FALSE)
  }
  take <- rep_len(as.numeric(round(take)), H)
  if (!is.numeric(n) || length(n) != H || anyNA(n) ||
      any(!is.finite(n)) || any(n <= 0)) {
    stop("'n' must be a positive finite numeric vector with length nrow(frame)",
         call. = FALSE)
  }
  n_h <- as.numeric(n)
  idx_of <- split(seq_len(nrow(psu)), factor(psu$stratum, levels = stratum))
  forced <- if ("certainty" %in% names(psu)) psu$certainty else NULL

  icc_row <- .psu_icc(measures, stratum, idx_of)
  # Kept as supplied, so the object records the caller's input.
  unit_cost_in <- unit_cost
  stage_cost <- c("cost_psu", "cost_ssu") %in% names(frame)
  if (all(stage_cost) && !is.null(unit_cost)) {
    stop(
      "certainty-aware allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }
  if (all(stage_cost)) {
    cost_psu <- rep_len(as.numeric(frame$cost_psu), H)
    cost_ssu <- rep_len(as.numeric(frame$cost_ssu), H)
    if (any(!is.finite(cost_psu)) || any(cost_psu <= 0) ||
          any(!is.finite(cost_ssu)) || any(cost_ssu <= 0)) {
      stop("'cost_psu' and 'cost_ssu' must be positive and finite",
           call. = FALSE)
    }
    unit_cost <- cost_psu / take + cost_ssu
  } else {
    cost_psu <- NULL
    cost_ssu <- NULL
  }
  base_frame <- frame[
    , setdiff(names(frame), c("n_per_psu", "cost_psu", "cost_ssu")),
    drop = FALSE
  ]
  base_measures <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  user_deff <- if ("deff" %in% names(measures)) {
    ifelse(is.na(measures$deff), deff, measures$deff)
  } else {
    rep(deff, nrow(measures))
  }
  row_h <- match(as.character(measures$stratum), stratum)
  if (anyNA(row_h)) {
    stop("'measures$stratum' must match 'frame$stratum'", call. = FALSE)
  }

  threshold <- take / (n_h / N_h)
  above <- psu$N >= threshold[match(psu$stratum, stratum)]
  certain <- if (is.null(forced)) above else above | forced
  source <- ifelse(
    above, "threshold",
    ifelse(if (is.null(forced)) FALSE else forced, "supplied", NA_character_)
  )
  source[!certain] <- NA_character_

  m <- base_measures
  m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
    h <- row_h[r]
    i <- idx_of[[h]]
    .psu_deff(
      n_h[h], N_h[h], sum(psu$N[i][certain[i]]), sum(psu$N[i][!certain[i]]),
      icc_row[r], take[h]
    )
  }, numeric(1))

  out <- .prec_alloc_bethel(
    frame = base_frame, n = n_h, measures = m, targets = targets,
    unit_cost = unit_cost, alpha = alpha, deff = deff, resp_rate = resp_rate,
    min_n_stratum = min_n_stratum, objective = objective, budget = budget,
    df = df, .allow_fractional_stages = .allow_fractional_stages,
    .finite_prop_var = FALSE
  )

  n_certain <- vapply(seq_len(H), function(h) sum(certain[idx_of[[h]]]),
                      numeric(1))
  n_all <- vapply(idx_of, length, numeric(1))
  # The take the assessed design fields in each PSU, by the same rule the
  # solver's operational design uses: a certainty PSU carries its own whole
  # take at the stratum rate, the remainder carries the stated n_per_psu.
  f_of <- (n_h / N_h)[match(psu$stratum, stratum)]
  n_take <- ifelse(
    certain,
    pmin(ceiling(f_of * psu$N), psu$N),
    take[match(psu$stratum, stratum)]
  )

  out$psu <- data.frame(
    stratum = psu$stratum,
    N = psu$N,
    certainty = certain,
    n_take = n_take,
    .certainty_source = source,
    .threshold = threshold[match(psu$stratum, stratum)],
    stringsAsFactors = FALSE
  )
  if (!is.null(psu$psu_id)) {
    out$psu <- cbind(
      data.frame(psu_id = psu$psu_id, stringsAsFactors = FALSE), out$psu
    )
  }
  out$psu$.distance <- out$psu$N / out$psu$.threshold - 1
  out$params$psu <- psu
  out$params$certainty <- list(
    n_psu_certain = n_certain,
    n_psu_rest = n_all - n_certain,
    threshold = threshold
  )
  out
}

#' The orbit member to hold, and whether any of them was feasible
#' @keywords internal
#' @noRd
.psu_pick_feasible <- function(
  orbit, base_frame, base_measures, targets, user_deff, row_h, idx_of,
  psu, N_h, icc_row, take, unit_cost, alpha, deff, resp_rate,
  min_n_stratum, objective, budget, df
) {
  ok <- logical(length(orbit))
  size <- vapply(orbit, function(s) sum(s$n_h), numeric(1))
  for (j in seq_along(orbit)) {
    st <- orbit[[j]]
    m <- base_measures
    m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
      h <- row_h[r]
      i <- idx_of[[h]]
      .psu_deff(
        st$n_h[h], N_h[h], sum(psu$N[i][st$certain[i]]),
        sum(psu$N[i][!st$certain[i]]), icc_row[r], take[h]
      )
    }, numeric(1))
    ev <- tryCatch(
      .prec_alloc_bethel(
        frame = base_frame, n = st$n_h, measures = m, targets = targets,
        unit_cost = unit_cost, alpha = alpha, deff = deff,
        resp_rate = resp_rate, min_n_stratum = min_n_stratum,
        objective = objective, budget = budget, df = df,
        .finite_prop_var = FALSE
      ),
      error = function(e) NULL
    )
    ok[j] <- !is.null(ev) && !is.null(ev$detail$.pass) && all(ev$detail$.pass)
  }
  if (any(ok)) {
    list(index = which(ok)[which.min(size[ok])], feasible = TRUE)
  } else {
    # No member met its own requirement, which a converged orbit can do
    # because the design effect moves with the allocation. Hold the most
    # certain one and let the settle step close the gap.
    certain_n <- vapply(orbit, function(s) sum(s$certain), numeric(1))
    list(index = which.max(certain_n), feasible = FALSE)
  }
}
