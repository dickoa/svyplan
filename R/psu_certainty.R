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
  if (!is.null(psu$psu_id) && anyDuplicated(psu$psu_id)) {
    dup <- unique(psu$psu_id[duplicated(psu$psu_id)])
    shown <- paste(sQuote(utils::head(as.character(dup), 5L)), collapse = ", ")
    if (length(dup) > 5L) {
      shown <- sprintf("%s and %d more", shown, length(dup) - 5L)
    }
    stop(
      sprintf(
        "'psu$psu_id' must name each PSU once, and %s %s more than once. A merge_psus() call that crossed strata gives one id to PSUs of two strata",
        shown, if (length(dup) == 1L) "appears" else "appear"
      ),
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
#' map that could alternate without converging. The closed form is computed
#' directly, not through the allocation: a value that moves in its last bits
#' as the allocation moves lets the settle step wander on a flat optimum,
#' which a precision result's pinned targets produce.
#' @keywords internal
#' @noRd
.psu_deff <- function(N_h, size_certain, size_rest, icc, take) {
  max((size_certain + size_rest * (1 + icc * (take - 1))) / N_h, 1)
}

#' Whole-unit counts of one stratum's field design, before any repair
#'
#' Shared by the classification closure and the operational design, so the
#' draw a PSU is classified against is the draw that is fielded. A remainder
#' always draws at least one PSU: the certainty takes are rounded up one by
#' one and can use up the whole allocation, and a remainder left undrawn
#' would be a part of the stratum no sample can reach.
#' @keywords internal
#' @noRd
.psu_stratum_counts <- function(n_h, N_h, take, sizes, certain, m = 0L) {
  held <- sizes[certain]
  per_psu <- pmin(ceiling(n_h / N_h * held), held)
  n_certain_int <- sum(per_psu)
  n_rest <- max(n_h - n_certain_int, 0)
  available <- sum(!certain)
  k <- max(ceiling(n_rest / take), as.numeric(available > 0))
  list(
    per_psu = per_psu,
    n_certain_int = n_certain_int,
    n_rest = n_rest,
    available = available,
    n_psu_draw = min(.psu_round_draw(k, m), available)
  )
}

#' Design effect of a stratum's field design, from its own takes
#'
#' The field design rounds each certainty take up and gives the remainder
#' what is left, so its two parts are not at the stratum rate that
#' `.psu_deff()` assumes. Each part is a stratum of its own, the certainty
#' PSUs sampled within at their takes and the remainder at its draw, and the
#' result is the design effect that puts their variance on the solver's
#' `deff * (1 / (r n) - 1 / N)` scale. At the stratum rate it is
#' `.psu_deff()`.
#' @keywords internal
#' @noRd
.psu_field_deff <- function(sizes, certain, per_psu, n_psu_draw, take, icc,
                            resp) {
  N_h <- sum(sizes)
  held <- sizes[certain]
  n_h <- sum(per_psu) + n_psu_draw * take
  scale <- N_h^2 * (1 / (resp * n_h) - 1 / N_h)
  if (scale <= 0) return(1)
  v_certain <- sum(held^2 * (1 / (resp * per_psu) - 1 / held))
  size_rest <- N_h - sum(held)
  v_rest <- if (n_psu_draw > 0) {
    size_rest^2 * (1 + icc * (take * resp - 1)) *
      (1 / (resp * n_psu_draw * take) - 1 / size_rest)
  } else {
    0
  }
  max((v_certain + v_rest) / scale, 1)
}

#' Refuse what the register model does not carry
#'
#' The register design effect has no slot for a variance ratio or for losing
#' a whole PSU, and its budget would have to price certainty visits as fixed
#' costs and round the remainder within the budget, which the field design
#' does not do.
#' @keywords internal
#' @noRd
.check_register_scope <- function(frame, measures, budget) {
  if (!is.null(budget)) {
    stop(
      "'budget' is not available with a PSU register: the field design rounds its draws up and does not hold certainty visits within a budget. Size the register design to 'targets', or drop 'psu' and give 'N_psu' instead",
      call. = FALSE
    )
  }
  ratio <- c(frame[["var_ratio_psu"]], measures[["var_ratio_psu"]])
  ratio <- ratio[!is.na(ratio)]
  if (length(ratio) > 0L && (!is.numeric(ratio) || any(ratio != 1))) {
    stop(
      "'var_ratio_psu' is not available with a PSU register: its design effect is built from 'icc_psu' and the take alone, so a variance ratio other than 1 would have no effect",
      call. = FALSE
    )
  }
  if (any(c("resp_rate_psu", "resp_rate_ssu") %in%
            c(names(frame), names(measures)))) {
    stop(
      "'resp_rate_psu' and 'resp_rate_ssu' are not available with a PSU register: losing a whole certainty PSU is not part of its variance model. Give the ultimate-unit response in 'resp_rate'",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' A remainder draw in whole zones of m PSUs
#' @keywords internal
#' @noRd
.psu_round_draw <- function(k, m) {
  if (m > 1L && k > 0) m * ceiling(k / m) else k
}

#' Zone of each remainder PSU under a draw of whole zones
#'
#' The remainder is sorted by `key` (size, largest first, when `NULL`) and
#' cut on cumulative size into `n_psu_draw / m` zones of about equal total,
#' a PSU going to the zone its cumulative size ends in (Valliant, Dever and
#' Kreuter 2018, Example 3.13). Certainty PSUs have no zone.
#' @keywords internal
#' @noRd
.psu_zones <- function(n_psu_draw, sizes, certain, m, key = NULL) {
  zone <- rep(NA_integer_, length(sizes))
  rest <- which(!certain)
  n_zone <- n_psu_draw %/% m
  if (!length(rest) || n_zone < 1L) return(zone)
  ord <- if (is.null(key)) order(-sizes[rest], rest) else order(key[rest], rest)
  idx <- rest[ord]
  cum <- cumsum(sizes[idx])
  width <- cum[length(cum)] / n_zone
  zone[idx] <- as.integer(pmin(n_zone, pmax(1, ceiling(cum / width - 1e-9))))
  zone
}

#' Remainder PSUs whose inclusion probability reaches the cutoff under a draw
#'
#' The tolerance is samplyr's, so a plan svyplan returns is one samplyr can
#' field without refusing it.
#' @keywords internal
#' @noRd
.psu_crossing <- function(n_psu_draw, sizes, certain, cutoff = 1, m = 0L,
                          key = NULL) {
  out <- logical(length(sizes))
  rest <- !certain
  if (n_psu_draw <= 0 || !any(rest)) return(out)
  level <- cutoff - 100 * .Machine$double.eps
  out[rest] <- n_psu_draw * sizes[rest] / sum(sizes[rest]) >= level
  if (m < 1L || any(out)) return(out)
  # Zones are cut only once the stratum-wide test passes. That test has then
  # removed every PSU of 1/m of a zone or more, so no PSU spans a zone
  # boundary and every zone is populated. It also covers a draw that would take the
  # whole remainder: probabilities summing to at least the number of PSUs
  # put one of them at one, so the closure works through such a remainder
  # as a census.
  zone <- .psu_zones(n_psu_draw, sizes, certain, m, key)
  for (z in unique(zone[rest])) {
    j <- which(zone == z)
    out[j] <- if (length(j) <= m) TRUE else m * sizes[j] / sum(sizes[j]) >= level
  }
  out
}

#' Close one stratum's classification under its whole-PSU remainder draw
#'
#' Rounding the remainder up to whole PSUs raises every remainder PSU's
#' inclusion probability, so a PSU below the size threshold can reach one in
#' the design that is fielded. It is certainty, and the draw is recomputed
#' without it until no PSU reaches one, the iterative rule for selection with
#' probability proportional to size.
#' @keywords internal
#' @noRd
.psu_draw_closure <- function(n_h, N_h, take, sizes, certain, cutoff = 1,
                              m = 0L, key = NULL) {
  operational <- logical(length(sizes))
  repeat {
    k <- .psu_stratum_counts(n_h, N_h, take, sizes, certain, m)$n_psu_draw
    hit <- .psu_crossing(k, sizes, certain, cutoff, m, key)
    if (!any(hit)) break
    certain <- certain | hit
    operational <- operational | hit
  }
  list(certain = certain, operational = operational)
}

#' Refuse register settings that have no register to apply to
#'
#' Nothing but a PSU register classifies PSUs or forms zones, so these
#' settings anywhere else would be ignored without a word.
#' @keywords internal
#' @noRd
.check_register_args <- function(certainty_cutoff, n_psu_per_zone, frame,
                                 psu) {
  if (!is.null(psu)) return(invisible(NULL))
  if (!is.null(certainty_cutoff) ||
        (is.data.frame(frame) && "certainty_cutoff" %in% names(frame))) {
    stop(
      "'certainty_cutoff' applies to a PSU register: supply 'psu', or drop it",
      call. = FALSE
    )
  }
  if (!is.null(n_psu_per_zone)) {
    stop(
      "'n_psu_per_zone' applies to a PSU register: supply 'psu', or drop it",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' PSUs drawn per zone, 0 meaning no zones
#' @keywords internal
#' @noRd
.psu_zone_m <- function(n_psu_per_zone) {
  if (is.null(n_psu_per_zone)) return(0L)
  if (!is.numeric(n_psu_per_zone) || length(n_psu_per_zone) != 1L ||
        is.na(n_psu_per_zone) || !n_psu_per_zone %in% c(1, 2)) {
    stop("'n_psu_per_zone' must be NULL, 1 or 2", call. = FALSE)
  }
  as.integer(n_psu_per_zone)
}

#' The order zones are cut in, from `psu$zone_order`
#'
#' `NULL` when the column is absent, which sorts by size, largest first.
#' @keywords internal
#' @noRd
.psu_zone_key <- function(psu, m) {
  key <- psu[["zone_order"]]
  if (is.null(key)) return(NULL)
  if (m < 1L) {
    stop(
      "'psu$zone_order' orders zones: supply 'n_psu_per_zone', or drop it",
      call. = FALSE
    )
  }
  if (anyNA(key)) {
    stop("'psu$zone_order' must have no missing values", call. = FALSE)
  }
  tied <- unique(psu$stratum[duplicated(data.frame(psu$stratum, key))])
  if (length(tied)) {
    stop(
      sprintf(
        "'psu$zone_order' must order the PSUs of a stratum without ties, as it does not in %s",
        paste(sQuote(tied), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  key
}

#' Per-stratum certainty cutoff, from the argument or the frame column
#'
#' `NULL` reads `frame$certainty_cutoff` where it exists, with `NA` and an
#' absent column meaning 1. An argument overrides the column, as `deff` and
#' `resp_rate` do.
#' @keywords internal
#' @noRd
.psu_cutoff <- function(cutoff, frame) {
  H <- nrow(frame)
  if (is.null(cutoff)) {
    column <- frame[["certainty_cutoff"]]
    cutoff <- if (is.null(column)) 1 else ifelse(is.na(column), 1, column)
  }
  if (!is.numeric(cutoff) || !length(cutoff) %in% c(1L, H) ||
        anyNA(cutoff) || any(!is.finite(cutoff)) ||
        any(cutoff <= 0 | cutoff > 1)) {
    stop(
      "'certainty_cutoff' must be one value or one per stratum, each in (0, 1]",
      call. = FALSE
    )
  }
  rep_len(as.numeric(cutoff), H)
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

#' The within-PSU take of every stratum, refused where a PSU cannot supply it
#'
#' A remainder PSU is drawn for a fixed take, so one smaller than the take
#' cannot field it. Such a PSU never reaches the size threshold, which is at
#' least the take, so it is refused before classification unless the caller
#' flags it certainty, where it is subsampled at the stratum's fraction.
#' @keywords internal
#' @noRd
.psu_take <- function(frame, psu, stratum) {
  take <- frame$n_per_psu
  if (!is.numeric(take) || anyNA(take) || any(take < 1) ||
        any(abs(take - round(take)) > 1e-8)) {
    stop("'n_per_psu' must contain positive whole numbers", call. = FALSE)
  }
  take <- rep_len(as.numeric(round(take)), length(stratum))
  small <- psu$N < take[match(psu$stratum, stratum)]
  if ("certainty" %in% names(psu)) small <- small & !psu$certainty
  if (any(small)) {
    label <- if (is.null(psu$psu_id)) {
      sprintf("row %d", which(small))
    } else {
      as.character(psu$psu_id[small])
    }
    shown <- paste(sQuote(utils::head(label, 5L)), collapse = ", ")
    if (length(label) > 5L) {
      shown <- sprintf("%s and %d more", shown, length(label) - 5L)
    }
    stop(
      sprintf(
        "%d PSU%s in 'psu' hold%s fewer units than the take 'n_per_psu' of %s stratum (%s). Merge each with a neighbouring PSU, as merge_psus() does, or flag it in 'psu$certainty'",
        length(label), if (length(label) == 1L) "" else "s",
        if (length(label) == 1L) "s" else "",
        if (length(label) == 1L) "its" else "their", shown
      ),
      call. = FALSE
    )
  }
  take
}

#' Variance groups of a one-PSU-per-zone design, fixed before selection
#'
#' Each zone yields one sampled PSU, so its variance is estimated by
#' collapsing zones in groups of two, or three where the count is odd. Zones
#' are grouped in zone order within a stratum. A stratum whose remainder is a
#' single zone has no partner inside it, so such strata are grouped with each
#' other in frame order, and a lone one joins the nearest stratum with zones.
#' The groups depend on the design only, never on sample values (Valliant,
#' Dever and Kreuter 2018, sec. 15.5.3). Certainty PSUs and PSUs outside any
#' zone get `NA`.
#' @keywords internal
#' @noRd
.psu_pairs <- function(zone, idx_of) {
  n_zone <- vapply(idx_of, function(i) {
    z <- zone[i]
    if (all(is.na(z))) 0L else as.integer(max(z, na.rm = TRUE))
  }, integer(1))
  group_of <- function(z, n) pmin(ceiling(z / 2), max(n %/% 2, 1))
  lone <- which(n_zone == 1L)
  zoned <- which(n_zone >= 2L)
  label <- vector("list", length(idx_of))
  for (h in zoned) {
    label[[h]] <- sprintf("%d.%d", h, group_of(seq_len(n_zone[h]), n_zone[h]))
  }
  if (length(lone) >= 2L) {
    g <- group_of(seq_along(lone), length(lone))
    for (s in seq_along(lone)) label[[lone[s]]] <- sprintf("lone.%d", g[s])
  } else if (length(lone) == 1L) {
    # The lone zone is the partner's first or last zone, and the partner's
    # zones are grouped again, so no group exceeds three.
    h <- lone
    after <- zoned[zoned > h]
    before <- zoned[zoned < h]
    if (length(after) || length(before)) {
      partner <- if (length(after)) after[1L] else before[length(before)]
      n <- n_zone[partner] + 1L
      g <- sprintf("%d.%d", partner, group_of(seq_len(n), n))
      if (length(after)) {
        label[[h]] <- g[1L]
        label[[partner]] <- g[-1L]
      } else {
        label[[h]] <- g[n]
        label[[partner]] <- g[-n]
      }
    } else {
      label[[h]] <- "lone.1"
    }
  }
  levels <- unique(unlist(label))
  pair <- rep(NA_integer_, length(zone))
  for (h in which(n_zone > 0L)) {
    i <- idx_of[[h]]
    z <- zone[i]
    pair[i] <- match(label[[h]][z], levels)
  }
  pair
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
  certainty_cutoff = NULL,
  n_psu_per_zone = NULL,
  max_iter = 30L,
  settle_iter = 200L,
  tolerance = 1e-9
) {
  psu <- .check_psu_table(psu, frame, measures)
  .check_register_scope(frame, measures, budget)

  stratum <- as.character(frame$stratum %||% seq_len(nrow(frame)))
  H <- nrow(frame)
  N_h <- frame$N
  take <- .psu_take(frame, psu, stratum)
  cutoff <- .psu_cutoff(certainty_cutoff, frame)
  zone_m <- .psu_zone_m(n_psu_per_zone)
  zone_key <- .psu_zone_key(psu, zone_m)
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
  } else {
    cost_psu <- NULL
    cost_ssu <- NULL
  }
  # Under stage costs a certainty PSU's visit is paid whatever its take, so
  # an interview there costs cost_ssu alone. Only the remainder's share of a
  # stratum's interviews brings visits with it.
  cost_at <- function(cls) {
    if (is.null(cost_psu)) return(unit_cost)
    rest <- vapply(seq_len(H), function(h) {
      i <- idx_of[[h]]
      sum(psu$N[i][!cls$certain[i]])
    }, numeric(1))
    cost_ssu + rest / N_h * cost_psu / take
  }
  base_frame <- frame[
    , setdiff(names(frame),
              c("n_per_psu", "cost_psu", "cost_ssu", "certainty_cutoff")),
    drop = FALSE
  ]
  icc_row <- .psu_icc(measures, stratum, idx_of)
  resp_row <- .joint_row_value("resp_rate", measures, frame, resp_rate)
  base_measures <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  # A design effect the caller supplies describes a source this model does not
  # carry, so it multiplies the clustering rather than replacing it.
  user_deff <- .joint_row_value("deff", measures, frame, deff)
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
        N_h[h], sum(psu$N[i][cls$certain[i]]),
        sum(psu$N[i][!cls$certain[i]]), icc_row[r], take[h] * resp_row[r]
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
      unit_cost = cost_at(cls), alpha = alpha, deff = deff,
      resp_rate = resp_rate, min_n_stratum = min_n_stratum,
      objective = objective, budget = budget, df = df,
      .finite_prop_var = FALSE
    )
  }

  classify_at <- function(n_h) {
    threshold <- cutoff * take / (n_h / N_h)
    certain <- logical(nrow(psu))
    source <- rep(NA_character_, nrow(psu))
    for (h in seq_len(H)) {
      i <- idx_of[[h]]
      one <- .psu_classify(psu$N[i], threshold[h], forced[i])
      closed <- .psu_draw_closure(
        n_h[h], N_h[h], take[h], psu$N[i], one$certain, cutoff[h],
        zone_m, zone_key[i]
      )
      certain[i] <- closed$certain
      source[i] <- ifelse(closed$operational, "operational", one$source)
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
    psu, N_h, icc_row, take, resp_row, alpha, deff, resp_rate,
    min_n_stratum, df
  )
  chosen <- orbit[[pick$index]]

  # Absorbing and settling again terminates: the held set only grows.
  settled <- chosen
  settle_used <- 0L
  absorbed <- 0L
  # PSUs the repaired field design pushed to probability one after the
  # classification was closed, held like the absorbed ones.
  repaired <- logical(nrow(psu))
  supplied <- if (is.null(forced)) logical(nrow(psu)) else forced

  # The whole-unit design. An element-level rounding of this allocation is
  # not fieldable: the remainder has to be a whole number of PSUs at the
  # stated take, and the certainty PSUs carry their own whole takes. It is
  # assessed on those takes, not at the stratum rate.
  assess <- function(n_try, per_psu, n_psu_draw) {
    m <- base_measures
    m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
      h <- row_h[r]
      i <- idx_of[[h]]
      .psu_field_deff(
        psu$N[i], settled$certain[i], per_psu[[h]], n_psu_draw[h], take[h],
        icc_row[r], resp_row[r]
      )
    }, numeric(1))
    ev <- .prec_alloc_bethel(
      frame = base_frame, n = n_try, measures = m,
      targets = targets, unit_cost = cost_at(settled), alpha = alpha,
      deff = deff,
      resp_rate = resp_rate, min_n_stratum = min_n_stratum,
      objective = objective, budget = budget, df = df,
      .allow_fractional_stages = TRUE, .finite_prop_var = FALSE
    )
    list(
      pass = all(ev$detail$.pass),
      constraints = ev$detail,
      feeds = ev$params$problem$A > 0,
      # The price of one more remainder PSU in each stratum.
      added = if (is.null(cost_psu)) {
        take * ev$params$cost_h
      } else {
        cost_psu + take * cost_ssu
      }
    )
  }

  repeat {
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
    by_field <- (implied$source %in% "operational" & implied$certain) |
      repaired
    settled$source <- ifelse(
      above, "threshold",
      ifelse(supplied, "supplied", ifelse(by_field, "operational", "orbit"))
    )
    settled$source[!settled$certain] <- NA_character_
    settled$above <- above | by_field

    operational <- .psu_operational(
      settled$n_h, N_h, take, psu, idx_of, settled$certain, assess, zone_m
    )
    # The closure classified against the draw before repair. A repair that
    # adds a PSU can push another to probability one, which is held and
    # settled like any absorbed PSU.
    crossing <- logical(nrow(psu))
    for (h in seq_len(H)) {
      i <- idx_of[[h]]
      crossing[i] <- .psu_crossing(
        operational$n_psu_draw[h], psu$N[i], settled$certain[i], cutoff[h],
        zone_m, zone_key[i]
      )
    }
    if (!any(crossing)) break
    repaired <- repaired | crossing
    settled$certain <- settled$certain | crossing
    absorbed <- absorbed + sum(crossing)
  }
  operational$constraints <- assess(
    operational$n_int, operational$per_psu, operational$n_psu_draw
  )$constraints

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
    stratum = stratum,
    certainty_cutoff = certainty_cutoff,
    n_psu_per_zone = n_psu_per_zone,
    zone_m = zone_m,
    zone_key = zone_key
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
#' Where the field design misses a target, the repair adds remainder PSUs by
#' the rule of the solver's own integer repair: each stratum that feeds a
#' failing target is tried with one more PSU, or one more zone, and the one
#' that most reduces the failing targets' excess variance per unit of cost is
#' kept. Every step adds a PSU and the remainder is finite, so the repair
#' ends. It stops with an error, naming the failing targets, when the
#' remainder feeding them is used up or no addition reduces them. That is a
#' limit of this repair under the held classification and certainty takes,
#' not a proof that no field design meets the targets.
#' @keywords internal
#' @noRd
.psu_operational <- function(n_h, N_h, take, psu, idx_of, certain,
                             assess, m = 0L) {
  H <- length(n_h)
  counts <- lapply(seq_len(H), function(h) {
    i <- idx_of[[h]]
    .psu_stratum_counts(n_h[h], N_h[h], take[h], psu$N[i], certain[i], m)
  })
  per_psu <- lapply(counts, `[[`, "per_psu")
  n_certain_int <- vapply(counts, `[[`, numeric(1), "n_certain_int")
  n_rest <- vapply(counts, `[[`, numeric(1), "n_rest")
  available <- vapply(counts, `[[`, numeric(1), "available")
  n_psu_draw <- vapply(counts, `[[`, numeric(1), "n_psu_draw")

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

  step <- max(m, 1L)
  grow <- function(draw, h) {
    draw[h] <- min(draw[h] + step, available[h])
    draw
  }
  excess <- function(ev) {
    ifelse(ev$constraints$.pass, 0, ev$constraints$.ratio^2 - 1)
  }
  repair <- 0L
  repeat {
    ev <- assess(n_certain_int + n_psu_draw * take, per_psu, n_psu_draw)
    if (isTRUE(ev$pass)) break
    over <- excess(ev)
    failing <- over > 0
    feeds <- rowSums(ev$feeds[, failing, drop = FALSE]) > 0
    # A stratum outside every failing target's domain cannot move it, since
    # each target's variance is a sum over the strata it covers.
    room <- which(feeds & n_psu_draw < available)
    if (length(room) == 0L) {
      .psu_repair_failure(ev$constraints, "used_up")
    }
    score <- vapply(room, function(h) {
      trial <- grow(n_psu_draw, h)
      tried <- assess(n_certain_int + trial * take, per_psu, trial)
      gain <- pmin(pmax(over - excess(tried), 0), over)
      sum(gain) / ((trial[h] - n_psu_draw[h]) * ev$added[h])
    }, numeric(1))
    if (!any(score > 0)) {
      .psu_repair_failure(ev$constraints, "no_gain")
    }
    n_psu_draw <- grow(n_psu_draw, room[which.max(score)])
    repair <- repair + 1L
  }

  list(
    n_int = n_certain_int + n_psu_draw * take,
    n_certain_int = n_certain_int,
    n_psu_draw = n_psu_draw,
    per_psu = per_psu,
    repair = repair,
    pass = TRUE
  )
}

#' Stop a repair that cannot meet the targets, naming them and why
#' @keywords internal
#' @noRd
.psu_repair_failure <- function(constraints, reason) {
  bad <- constraints[!constraints$.pass, , drop = FALSE]
  where <- ifelse(
    bad$domain == ".overall", "whole population",
    sprintf("%s = %s", bad$domain, bad$level)
  )
  lines <- sprintf(
    "  %s (%s): %s %s, required at most %s",
    bad$constraint, where, bad$.metric, signif(bad$.achieved, 4),
    signif(bad$.target, 4)
  )
  why <- switch(
    reason,
    used_up = "every stratum feeding these targets already draws its whole remainder",
    no_gain = "no remainder PSU or zone that can still be added reduces them"
  )
  stop(
    paste0(
      "no whole-unit repair meets every target under the current certainty classification and fixed takes: ",
      why, ". Failing targets:\n", paste(lines, collapse = "\n")
    ),
    call. = FALSE
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
  if (x$zone_m >= 1L) fit$detail$n_zone <- op$n_psu_draw %/% x$zone_m
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
    # The continuous optimum is priced on the cost the solver minimized, the
    # certainty visits fixed and every interview at its stratum's marginal
    # cost, so the two readings on a printed block differ by integerization
    # and not by basis.
    fit$params$achieved$cost <- sum(
      n_certain * x$cost_psu + fit$params$cost_h * fit$detail$n
    )
    fit$optimization$cost <- fit$params$achieved$cost
  }
  if (!is.null(op$constraints)) fit$operational$constraints <- op$constraints
  fit$operational$all_pass <- op$pass
  fit$operational$repair_iterations <- op$repair

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
  if (x$zone_m >= 1L) {
    # The zones the field design draws from, cut at its final draw.
    zone <- rep(NA_integer_, nrow(x$psu))
    for (h in seq_along(x$idx_of)) {
      i <- x$idx_of[[h]]
      zone[i] <- .psu_zones(
        op$n_psu_draw[h], x$psu$N[i], x$certainty[i], x$zone_m, x$zone_key[i]
      )
    }
    fit$psu$.zone <- zone
    if (x$zone_m == 1L) fit$psu$.pair <- .psu_pairs(zone, x$idx_of)
  }
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
  # As supplied, so a rebuild resolves it against the frame the same way.
  fit$params$certainty_cutoff <- x$certainty_cutoff
  fit$params$n_psu_per_zone <- x$n_psu_per_zone
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
  certainty_cutoff = NULL,
  n_psu_per_zone = NULL,
  .allow_fractional_stages = FALSE
) {
  if (is.null(n)) stop("'n' is required", call. = FALSE)
  psu <- .check_psu_table(psu, frame, measures)
  .check_register_scope(frame, measures, budget)

  stratum <- as.character(frame$stratum %||% seq_len(nrow(frame)))
  H <- nrow(frame)
  N_h <- frame$N
  take <- .psu_take(frame, psu, stratum)
  if (!is.numeric(n) || length(n) != H || anyNA(n) ||
      any(!is.finite(n)) || any(n <= 0)) {
    stop("'n' must be a positive finite numeric vector with length nrow(frame)",
         call. = FALSE)
  }
  n_h <- as.numeric(n)
  cutoff <- .psu_cutoff(certainty_cutoff, frame)
  zone_m <- .psu_zone_m(n_psu_per_zone)
  zone_key <- .psu_zone_key(psu, zone_m)
  idx_of <- split(seq_len(nrow(psu)), factor(psu$stratum, levels = stratum))
  forced <- if ("certainty" %in% names(psu)) psu$certainty else NULL

  icc_row <- .psu_icc(measures, stratum, idx_of)
  resp_row <- .joint_row_value("resp_rate", measures, frame, resp_rate)
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
  } else {
    cost_psu <- NULL
    cost_ssu <- NULL
  }
  base_frame <- frame[
    , setdiff(names(frame),
              c("n_per_psu", "cost_psu", "cost_ssu", "certainty_cutoff")),
    drop = FALSE
  ]
  base_measures <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  user_deff <- .joint_row_value("deff", measures, frame, deff)
  row_h <- match(as.character(measures$stratum), stratum)
  if (anyNA(row_h)) {
    stop("'measures$stratum' must match 'frame$stratum'", call. = FALSE)
  }

  threshold <- cutoff * take / (n_h / N_h)
  above <- psu$N >= threshold[match(psu$stratum, stratum)]
  certain <- if (is.null(forced)) above else above | forced
  source <- ifelse(
    above, "threshold",
    ifelse(if (is.null(forced)) FALSE else forced, "supplied", NA_character_)
  )
  for (h in seq_len(H)) {
    i <- idx_of[[h]]
    closed <- .psu_draw_closure(
      n_h[h], N_h[h], take[h], psu$N[i], certain[i], cutoff[h],
      zone_m, zone_key[i]
    )
    source[i[closed$operational]] <- "operational"
    certain[i] <- closed$certain
  }
  source[!certain] <- NA_character_

  m <- base_measures
  m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
    h <- row_h[r]
    i <- idx_of[[h]]
    .psu_deff(
      N_h[h], sum(psu$N[i][certain[i]]), sum(psu$N[i][!certain[i]]),
      icc_row[r], take[h] * resp_row[r]
    )
  }, numeric(1))

  n_certain <- vapply(seq_len(H), function(h) sum(certain[idx_of[[h]]]),
                      numeric(1))
  # A certainty PSU's visit is paid whatever its take, as in n_alloc().
  if (!is.null(cost_psu)) {
    rest <- vapply(seq_len(H), function(h) {
      i <- idx_of[[h]]
      sum(psu$N[i][!certain[i]])
    }, numeric(1))
    unit_cost <- cost_ssu + rest / N_h * cost_psu / take
  }

  out <- .prec_alloc_bethel(
    frame = base_frame, n = n_h, measures = m, targets = targets,
    unit_cost = unit_cost, alpha = alpha, deff = deff, resp_rate = resp_rate,
    min_n_stratum = min_n_stratum, objective = objective, budget = budget,
    df = df, .allow_fractional_stages = .allow_fractional_stages,
    .finite_prop_var = FALSE
  )

  if (!is.null(cost_psu)) {
    out$params$achieved$cost <- out$params$achieved$cost +
      sum(n_certain * cost_psu)
  }
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
  if (zone_m >= 1L) {
    zone <- rep(NA_integer_, nrow(psu))
    for (h in seq_len(H)) {
      i <- idx_of[[h]]
      k <- .psu_stratum_counts(
        n_h[h], N_h[h], take[h], psu$N[i], certain[i], zone_m
      )$n_psu_draw
      zone[i] <- .psu_zones(k, psu$N[i], certain[i], zone_m, zone_key[i])
    }
    out$psu$.zone <- zone
    if (zone_m == 1L) out$psu$.pair <- .psu_pairs(zone, idx_of)
  }
  out$params$psu <- psu
  out$params$certainty_cutoff <- certainty_cutoff
  out$params$n_psu_per_zone <- n_psu_per_zone
  # The caller's tables, not the solver's stripped copies, which drop the take
  # and the ICC and bake the design effect in, so n_alloc() can rebuild a
  # register plan from the result.
  out$params$frame <- frame
  out$params$measures <- measures
  out$params$unit_cost <- unit_cost_in
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
  psu, N_h, icc_row, take, resp_row, alpha, deff, resp_rate,
  min_n_stratum, df
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
        N_h[h], sum(psu$N[i][st$certain[i]]),
        sum(psu$N[i][!st$certain[i]]), icc_row[r], take[h] * resp_row[r]
      )
    }, numeric(1))
    ev <- tryCatch(
      .prec_alloc_bethel(
        frame = base_frame, n = st$n_h, measures = m, targets = targets,
        alpha = alpha, deff = deff,
        resp_rate = resp_rate, min_n_stratum = min_n_stratum, df = df,
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
