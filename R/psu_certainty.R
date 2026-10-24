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
  if (!is.null(psu[["psu_id"]]) && anyDuplicated(psu$psu_id)) {
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
  stratum <- as.character(frame[["stratum"]] %||% seq_len(nrow(frame)))
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
#' allocation. The closed form is computed directly, not through the
#' allocation, so the continuous solve sees one value per classification and
#' not one that moves in its last bits with the allocation.
#' @keywords internal
#' @noRd
.psu_deff <- function(N_h, size_certain, size_rest, icc, take) {
  max((size_certain + size_rest * (1 + icc * (take - 1))) / N_h, 1)
}

#' Whole-unit counts of one stratum's field design
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
#' The register design effect has no slot for a variance ratio, for losing
#' a whole PSU, or for an eligible-units rate. Every register route reads
#' its tables through here, overrides of a fitted plan included.
#' @keywords internal
#' @noRd
.check_register_scope <- function(frame, measures) {
  if ("eligible_rate" %in% c(names(frame), names(measures))) {
    .check_eligible_scope(1, frame, measures, psu = TRUE)
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
    label <- if (is.null(psu[["psu_id"]])) {
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
        "%d PSU%s in 'psu' hold%s fewer units than the take 'n_per_psu' of %s stratum (%s). Merge each with a neighboring PSU, as merge_psus() does, or flag it in 'psu$certainty'",
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
#' The field design is chosen whole. Each stratum's menu lists the designs
#' the closure gives with its `t` largest PSUs held in certainty, and the
#' search picks one design per stratum, the cheapest meeting every target
#' (`.psu_plan()`). A continuous solve with every PSU in the remainder
#' anchors the listing's first pass. The continuous solve under the planned
#' classification is kept as the relaxation `$detail$n` and
#' `$optimization` describe.
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
  control = NULL
) {
  .check_measures_keys(measures)
  psu <- .check_psu_table(psu, frame, measures)
  .check_register_scope(frame, measures)
  ev <- .psu_setup(
    frame, psu, measures, targets, unit_cost = unit_cost, alpha = alpha,
    deff = deff, resp_rate = resp_rate, min_n_stratum = min_n_stratum,
    objective = objective, df = df, certainty_cutoff = certainty_cutoff,
    n_psu_per_zone = n_psu_per_zone
  )
  max_nodes <- .search_control(control)$max_nodes
  stratum <- as.character(frame[["stratum"]] %||% seq_len(nrow(frame)))
  H <- nrow(frame)
  N_h <- frame$N
  take <- ev$take
  idx_of <- ev$idx_of
  icc_row <- .psu_icc(measures, stratum, idx_of)
  resp_row <- .joint_row_value("resp_rate", measures, frame, resp_rate)
  user_deff <- .joint_row_value("deff", measures, frame, deff)
  row_h <- match(as.character(measures$stratum), stratum)
  if (anyNA(row_h)) {
    stop("'measures$stratum' must match 'frame$stratum'", call. = FALSE)
  }
  base_frame <- frame[
    , setdiff(names(frame),
              c("n_per_psu", "cost_psu", "cost_ssu", "certainty_cutoff")),
    drop = FALSE
  ]
  base_measures <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  cost_psu <- ev$cost_psu
  cost_ssu <- ev$cost_ssu

  # Under stage costs a certainty PSU's visit is paid whatever its take, so
  # an interview there costs cost_ssu alone. Only the remainder's share of a
  # stratum's interviews brings visits with it, and the certainty visits
  # come out of a budget first.
  visits_of <- function(certain) {
    if (is.null(cost_psu)) return(0)
    sum(certain * cost_psu[match(psu$stratum, stratum)])
  }
  solve_at <- function(certain, budget, objective, whole = FALSE) {
    m <- base_measures
    m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
      h <- row_h[r]
      i <- idx_of[[h]]
      .psu_deff(
        N_h[h], sum(psu$N[i][certain[i]]), sum(psu$N[i][!certain[i]]),
        icc_row[r], take[h] * resp_row[r]
      )
    }, numeric(1))
    unit <- if (is.null(cost_psu)) {
      unit_cost
    } else {
      rest <- vapply(seq_len(H), function(h) {
        i <- idx_of[[h]]
        sum(psu$N[i][!certain[i]])
      }, numeric(1))
      cost_ssu + rest / N_h * cost_psu / take
    }
    if (!is.null(budget)) budget <- budget - visits_of(certain)
    # Although the inner solver receives an element-shaped frame, this is
    # a multistage model. Retain its p(1-p) working variance on every call.
    .n_alloc_bethel(
      frame = base_frame, measures = m, targets = targets,
      unit_cost = unit, alpha = alpha, deff = deff,
      resp_rate = resp_rate, min_n_stratum = min_n_stratum,
      objective = objective, budget = budget, df = df,
      .finite_prop_var = FALSE, .whole = whole
    )
  }

  # In budget mode the continuous solve with no certainty can miss targets
  # that certainty makes affordable, so it falls back to target mode, and
  # without targets to the stratum floors.
  flat_certain <- rep(FALSE, nrow(psu))
  flat <- if (is.null(budget)) {
    solve_at(flat_certain, NULL, objective)
  } else {
    tryCatch(solve_at(flat_certain, budget, objective),
             error = function(e) NULL)
  }
  f_anchor <- if (!is.null(flat)) {
    flat$detail$n / N_h
  } else if (!is.null(targets)) {
    solve_at(flat_certain, NULL, NULL)$detail$n / N_h
  } else {
    ev$lower / N_h
  }
  plan <- .psu_plan(ev, f_anchor, max_nodes, budget)
  design <- lapply(seq_len(H), function(h) .psu_design(ev, h, plan$rows[[h]]))
  certain <- logical(nrow(psu))
  source <- rep(NA_character_, nrow(psu))
  for (h in seq_len(H)) {
    certain[idx_of[[h]]] <- design[[h]]$certain
    source[idx_of[[h]]] <- design[[h]]$source
  }
  fit <- tryCatch(solve_at(certain, budget, objective),
                  error = function(e) NULL)
  if (is.null(fit) && !is.null(budget) && !is.null(targets)) {
    # The field design prices its visits exactly, so the continuous model
    # under the same certainty can need more than the budget. It is then
    # solved at the cost of its cheapest whole-unit design meeting the
    # targets, which its own rounding needs.
    fit <- tryCatch({
      least <- solve_at(certain, NULL, NULL, whole = TRUE)$operational$cost +
        visits_of(certain)
      solve_at(certain, least * (1 + 1e-9), objective)
    }, error = function(e) NULL)
  }
  fit <- fit %||% flat
  if (is.null(fit)) {
    stop("internal error: the continuous relaxation has no solution under the budget",
         call. = FALSE)
  }

  per_psu <- lapply(design, `[[`, "per_psu")
  n_psu_draw <- vapply(design, `[[`, numeric(1), "n_psu_draw")
  n_certain_int <- vapply(design, `[[`, numeric(1), "n_certain_int")
  n_int <- n_certain_int + n_psu_draw * take
  assessed <- .psu_field_assess(
    frame, psu, measures, targets, unit_cost, alpha, deff, resp_rate,
    min_n_stratum, objective, budget, df, certain, per_psu, n_psu_draw
  )
  if (!all(assessed$detail$.pass)) {
    stop("internal error: the planned field design misses a target",
         call. = FALSE)
  }
  if (!is.null(budget) &&
        abs(assessed$objective_value - plan$value) >
          1e-9 * max(abs(plan$value), 1e-300)) {
    stop("internal error: the planned field design's objective does not reassess",
         call. = FALSE)
  }

  list(
    fit = fit,
    frame = frame,
    measures = measures,
    operational = list(
      n_int = n_int,
      n_certain_int = n_certain_int,
      n_psu_draw = n_psu_draw,
      per_psu = per_psu,
      constraints = assessed$detail,
      pass = TRUE,
      objective = assessed$objective,
      objective_value = assessed$objective_value
    ),
    budget = budget,
    plan = plan,
    take = take,
    unit_cost = unit_cost,
    cost_psu = cost_psu,
    cost_ssu = cost_ssu,
    psu = psu,
    certainty = certain,
    source = source,
    threshold = vapply(design, `[[`, numeric(1), "threshold"),
    idx_of = idx_of,
    stratum = stratum,
    certainty_cutoff = certainty_cutoff,
    n_psu_per_zone = n_psu_per_zone,
    control = control,
    zone_m = ev$zone_m,
    zone_key = .psu_zone_key(psu, ev$zone_m),
    max_nodes = max_nodes
  )
}

#' Warn once, naming the strata whose remainder draw is a single PSU
#'
#' The register counterpart of `.warn_singleton_strata()`. Certainty PSUs
#' carry their own within-PSU variance, so only a remainder of one drawn PSU
#' leaves the stratum without an estimate.
#' @keywords internal
#' @noRd
.warn_single_draw <- function(labels) {
  if (length(labels) == 0L) {
    return(invisible(NULL))
  }
  many <- length(labels) > 1L
  warning(warningCondition(
    sprintf(
      "%s %s %s a single PSU outside certainty, so no within-stratum variance can be estimated there. Set n_psu_per_zone = 2 to draw at least two, or collapse %s with a neighbor",
      if (many) "strata" else "stratum",
      paste0("'", labels, "'", collapse = ", "),
      if (many) "each draw" else "draws",
      if (many) "them" else "it"
    ),
    class = "svyplan_single_draw"
  ))
  invisible(NULL)
}

#' Fold a planned field design into the allocation object it refines
#'
#' The class does not change, so `predict()`, `summary()` and `design_df()`
#' keep working. samplyr refuses these fits as stage sizes, because no
#' per-stratum total carries a design of whole certainty takes plus a PPS
#' remainder. It fields the plan from `$psu` instead, which is why that table
#' carries the exact per-PSU takes. `$detail` gains the per-stratum counts and
#' the threshold, `$psu` is the per-PSU classification, and the search joins
#' `$optimization`, which already reports how the continuous solve went.
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

  op <- x$operational
  fit$detail$n_int <- op$n_int
  fit$detail$n_certain_int <- op$n_certain_int
  fit$detail$n_psu_draw <- op$n_psu_draw
  if (x$zone_m >= 1L) fit$detail$n_zone <- op$n_psu_draw %/% x$zone_m
  fit$operational$n <- sum(op$n_int)
  fit$operational$n_certain <- sum(op$n_certain_int)
  fit$operational$n_psu_draw <- sum(op$n_psu_draw)
  # Under n_psu_per_zone = 1 a single draw is the design, and $psu$.pair
  # groups it for variance estimation.
  if (x$zone_m != 1L) {
    .warn_single_draw(fit$detail$stratum[op$n_psu_draw == 1])
  }
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
  fit$operational$constraints <- op$constraints
  fit$operational$all_pass <- op$pass
  fit$operational$status <- x$plan$status
  if (!is.null(x$budget)) {
    fit$operational$objective <- op$objective
    fit$operational$objective_value <- op$objective_value
    fit$operational$budget_residual <- x$budget - fit$operational$cost
    if (fit$operational$budget_residual < -1e-9 * x$budget) {
      stop("internal error: the planned field design exceeds the budget",
           call. = FALSE)
    }
    # The continuous solve saw the budget less the certainty visits, or a
    # larger one when the targets needed it.
    fit$params$budget <- x$budget
    fit$optimization$budget <- x$budget
    fit$optimization$budget_residual <- x$budget - fit$params$achieved$cost
  }

  n_take <- numeric(nrow(x$psu))
  for (h in seq_along(x$idx_of)) {
    i <- x$idx_of[[h]]
    n_take[i] <- x$take[h]
    n_take[i[x$certainty[i]]] <- op$per_psu[[h]]
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
  if (!is.null(x$psu[["psu_id"]])) {
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

  plan <- x$plan
  fit$optimization$certainty <- NULL
  fit$optimization$search <- list(
    family = if (is.null(x$cost_psu)) "closure" else "voluntary",
    status = plan$status,
    bound = plan$bound,
    gap = if (plan$value > 0) (plan$value - plan$bound) / plan$value else 0,
    nodes = plan$nodes,
    max_nodes = plan$max_nodes,
    states = plan$states,
    kept = plan$kept
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
  fit$params$control <- x$control
  fit
}

#' Assess a whole-unit field design on its own takes
#'
#' Each row's design effect is `.psu_field_deff()` on the design's
#' certainty takes and remainder draw, times the caller's, and the
#' element assessment does the rest. This is the precision the planner
#' sizes and the one `prec_alloc()` of a register fit reports.
#' @keywords internal
#' @noRd
.psu_field_assess <- function(frame, psu, measures, targets, unit_cost,
                              alpha, deff, resp_rate, min_n_stratum,
                              objective, budget, df, certain, per_psu,
                              n_psu_draw) {
  stratum <- as.character(frame[["stratum"]] %||% seq_len(nrow(frame)))
  idx_of <- split(seq_len(nrow(psu)), factor(psu$stratum, levels = stratum))
  take <- .psu_take(frame, psu, stratum)
  icc_row <- .psu_icc(measures, stratum, idx_of)
  resp_row <- .joint_row_value("resp_rate", measures, frame, resp_rate)
  user_deff <- .joint_row_value("deff", measures, frame, deff)
  row_h <- match(as.character(measures$stratum), stratum)
  m <- measures[
    , setdiff(names(measures), c("icc_psu", "var_ratio_psu")), drop = FALSE
  ]
  m$deff <- user_deff * vapply(seq_len(nrow(m)), function(r) {
    h <- row_h[r]
    i <- idx_of[[h]]
    .psu_field_deff(
      psu$N[i], certain[i], per_psu[[h]], n_psu_draw[h], take[h],
      icc_row[r], resp_row[r]
    )
  }, numeric(1))
  n_int <- vapply(per_psu, sum, numeric(1)) + n_psu_draw * take
  .prec_alloc_bethel(
    frame = frame[
      , setdiff(names(frame),
                c("n_per_psu", "cost_psu", "cost_ssu", "certainty_cutoff")),
      drop = FALSE
    ],
    n = n_int, measures = m, targets = targets, unit_cost = unit_cost,
    alpha = alpha, deff = deff, resp_rate = resp_rate,
    min_n_stratum = min_n_stratum, objective = objective, budget = budget,
    df = df, .allow_fractional_stages = TRUE, .finite_prop_var = FALSE
  )
}

#' The precision of a register fit's own field design
#'
#' `prec_alloc()` of a register fit with nothing that changes the design
#' overridden. The design is the one in `$psu` and `$detail`, so with no
#' override the result reports what `$operational$constraints` reports,
#' with the fit's `$psu` table. `overrides` replace the fit's own
#' assessment inputs, such as `alpha`, `df`, `deff` or `targets`. Under a
#' unit cost the design's takes are priced at the assessed `unit_cost`,
#' and under stage costs at the fit's field cost.
#' @keywords internal
#' @noRd
.prec_alloc_psu_design <- function(fit, overrides = list()) {
  p <- fit$params
  p[names(overrides)] <- overrides
  # An override passes the checks the construction from totals applies.
  .check_measures_keys(p$measures)
  .check_register_scope(p$frame, p$measures)
  if (all(c("cost_psu", "cost_ssu") %in% names(p$frame)) &&
        !is.null(p$unit_cost)) {
    stop(
      "certainty-aware allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }
  psu <- .check_psu_table(p$psu, p$frame, p$measures)
  stratum <- as.character(p$frame[["stratum"]] %||% seq_len(nrow(p$frame)))
  idx_of <- split(seq_len(nrow(psu)), factor(psu$stratum, levels = stratum))
  certain <- fit$psu$certainty
  per_psu <- lapply(idx_of, function(i) fit$psu$n_take[i][certain[i]])
  out <- .psu_field_assess(
    p$frame, psu, p$measures, p$targets, p$unit_cost, p$alpha, p$deff,
    p$resp_rate, p$min_n_stratum, p$objective, p$budget, p$df, certain,
    per_psu, fit$detail$n_psu_draw
  )
  n_certain <- vapply(idx_of, function(i) sum(certain[i]), numeric(1))
  stage_cost <- all(c("cost_psu", "cost_ssu") %in% names(p$frame))
  cost <- if (stage_cost) fit$operational$cost else out$params$achieved$cost
  out$params$achieved$cost <- cost
  if (!is.null(p$budget)) {
    out$params$budget_residual <- p$budget - cost
  }
  out$psu <- fit$psu
  out$params$psu <- p$psu
  out$params$certainty_cutoff <- p$certainty_cutoff
  out$params$n_psu_per_zone <- p$n_psu_per_zone
  out$params$control <- p[["control"]]
  out$params$frame <- p$frame
  out$params$measures <- p$measures
  out$params$unit_cost <- p$unit_cost
  out$params$certainty <- list(
    n_psu_certain = n_certain,
    n_psu_rest = lengths(idx_of) - n_certain,
    threshold = fit$detail$threshold,
    n_psu_draw = fit$detail$n_psu_draw
  )
  out
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
  .check_measures_keys(measures)
  psu <- .check_psu_table(psu, frame, measures)
  .check_register_scope(frame, measures)

  stratum <- as.character(frame[["stratum"]] %||% seq_len(nrow(frame)))
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
    if (!is.null(budget)) {
      out$params$budget_residual <- budget - out$params$achieved$cost
    }
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
  if (!is.null(psu[["psu_id"]])) {
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
