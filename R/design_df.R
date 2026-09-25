#' Design degrees of freedom
#'
#' Count nominal degrees of freedom from a planned design's issued units.
#' This is a planning counterpart to `survey::degf()`. Realized degrees of
#' freedom depend on the responding sample and the analysis variance estimator.
#'
#' @param x A `svyplan` result to count from: a [n_cluster()] or
#'   [prec_cluster()] allocation, a [n_alloc()] allocation, a
#'   [strata_bound()] stratification, a [n_twophase()] allocation, or a
#'   single-indicator [n_prop()] or [n_mean()] size. `NULL` (default)
#'   counts from the components below.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#'
#' @return A `svyplan_df` object: a numeric scalar carrying the counts it
#'   was formed from. Use it directly wherever a `df` argument is expected,
#'   [as.double()] to strip it to a plain number, and `$strata` and
#'   `$domains` for the per-stratum and per-domain tables.
#'
#' @details
#' The degrees of freedom of a design-based variance estimator are the
#' number of independent units contributing to it, less the number of
#' constraints the design imposes: for a stratified design, sampled units
#' minus strata. In a clustered design the units are the PSUs, not the
#' ultimate units, which is why a survey of 12,000 households in 300
#' clusters across 20 strata has 280 degrees of freedom rather than 11,999.
#'
#' Pass the count as `df` to use a `t` quantile in place of the normal one,
#' and to the Korn-Graubard interval, which widens as the count falls. It is
#' the nominal first-stage count. The degrees of freedom of a particular
#' estimator, such as a Satterthwaite approximation combining several
#' variance components, differ from it.
#'
#' ## What is counted, per design shape
#'
#' | Design | Counted | df |
#' | --- | --- | --- |
#' | Unstratified element | units | \eqn{n - 1} |
#' | Unstratified cluster | PSUs | \eqn{m - 1} |
#' | Stratified element | units per stratum | \eqn{\sum_h n_h - H} |
#' | Stratified cluster | PSUs per stratum | \eqn{\sum_h m_h - H} |
#' | Two-phase | phase-2 units per stratum | \eqn{\sum_h n_{2h} - H}{sum_h n_2h - H} |
#'
#' The counts are read from the whole-unit columns a plan reports
#' (`n_int`, `n_psu_int`, `n_psu_draw`, `operational$n`), never from the continuous
#' optimum, since a degree of freedom is a unit that will be fielded. They
#' are issued counts and are not adjusted for nonresponse. For an element
#' design, fewer respondents generally mean fewer degrees of freedom. For
#' a cluster design, element nonresponse need not change the PSU count,
#' but lost PSUs or PSUs with no responding domain members can reduce it.
#' The two-phase count is a phase-2 planning convention, not a derivation of
#' the degrees of freedom of the combined two-phase variance estimator.
#'
#' A stratum that carries no sampling variance drops out of both terms. In
#' element mode a `take_all` stratum is a census and contributes nothing.
#' In cluster mode `take_all` marks an element-level census *inside* the
#' stratum, which leaves the PSU stage sampling as before, so such a
#' stratum still contributes \eqn{m_h - 1}. A phase-2 `take_all` stratum
#' likewise still contributes, since following up every phase-1 unit in a
#' stratum does not remove the phase-1 sampling variance.
#'
#' ## Per-domain degrees of freedom
#'
#' `$domains` applies the same nominal counting rule by domain. The allocation API
#' expresses a domain as a set of whole strata, so a domain's df is that
#' union's own contribution, \eqn{\sum_{h \in d} m_h - H_d}{sum_(h in d) m_h - H_d}, and the
#' per-domain values sum to the overall df whenever the domains partition
#' the frame.
#'
#' The case that is *not* derivable is an analytic domain cutting across
#' strata, such as an age group or sex. Its df is bounded above by the
#' design df, and the Korn-Graubard counting rule for it (PSUs containing
#' domain members, minus strata containing them) needs frame information
#' that a plan does not carry. It is a limit of what is knowable before
#' fielding, not an omission.
#'
#' ## What the count cannot see
#'
#' Certainty, or self-representing, PSUs contribute no between-PSU variance
#' and should not count toward first-stage df. For a PSU register supplied
#' through [n_alloc()] with no certainty PSUs, this helper counts the
#' operational PSUs (`n_psu_draw`) minus strata. Allocations containing
#' certainty PSUs remain unsupported: combining their within-PSU variance
#' contributions with other strata requires an explicit degrees-of-freedom
#' convention that this helper does not yet implement. In other inputs,
#' unrecorded certainty PSUs can make the nominal first-stage count too large.
#'
#' A stratum holding a single PSU supports no within-stratum variance
#' estimate at all, and warns, naming the stratum. Its own contribution is
#' zero and `$strata` marks it `"singleton"`. Whether a small but positive
#' df is a problem is a judgment call, so it is reported and never warned
#' about.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 3.
#'
#' Korn, E. L. and Graubard, B. I. (1998). Confidence intervals for
#' proportions with small expected number of positive counts estimated from
#' survey data. *Survey Methodology*, 24(2), 193--201.
#'
#' @family design input functions
#' @seealso [design_effect()] for the other property of a design read off a
#'   plan, [n_prop()] and [n_mean()], whose `df` argument this is the input
#'   to, and [varcomp()] for the components a clustered plan is built from.
#'
#' @examples
#' # A stratified cluster allocation: PSUs minus strata
#' frame <- data.frame(
#'   stratum = c("a", "b", "c"),
#'   N       = c(12000, 30000, 8000),
#'   sd      = c(5, 8, 3),
#'   mean    = c(10, 20, 5),
#'   icc_psu = 0.05,
#'   n_per_psu = 12
#' )
#' alloc <- n_alloc(frame, n = 3000)
#' design_df(alloc)
#' design_df(alloc)$strata
#'
#' # Drop-in wherever a df is expected
#' n_prop(p = 0.02, moe = 0.01, method = "beta", df = design_df(alloc))
#'
#' # From counts, with no plan in hand
#' design_df(n_psu = 300, n_strata = 20)
#'
#' @export
design_df <- function(x = NULL, ...) {
  # With no object to dispatch on, UseMethod() would dispatch on the first
  # element of ..., which here is a count, not an object.
  if (missing(x) || is.null(x)) {
    return(design_df.default(NULL, ...))
  }
  UseMethod("design_df")
}

#' @describeIn design_df Count from the design's own numbers. Supply
#'   `n_psu` for a clustered design or `n` for an element design, with
#'   `n_strata` where the design is stratified.
#'
#' @param n_psu Number of PSUs the design will select, across all strata.
#' @param n_strata Number of strata, default 1 (unstratified).
#' @param n Number of ultimate units, for an element design. Supply this or
#'   `n_psu`, not both.
#'
#' @export
design_df.default <- function(x = NULL, ..., n_psu = NULL, n_strata = 1,
                              n = NULL) {
  .check_unused_dots(...)
  if (!is.null(x)) {
    stop(
      "'x' must be a svyplan result or NULL. Pass the counts by name (n_psu or n, with n_strata)",
      call. = FALSE
    )
  }
  if (is.null(n_psu) == is.null(n)) {
    stop(
      "supply exactly one of 'n_psu' (a clustered design) or 'n' (an element design)",
      call. = FALSE
    )
  }
  stage <- if (is.null(n_psu)) "element" else "psu"
  units <- n_psu %||% n
  check_count(units, if (stage == "psu") "n_psu" else "n")
  check_count(n_strata, "n_strata")
  if (n_strata > units) {
    stop(
      sprintf(
        "'n_strata' (%g) cannot exceed the %d units counted", n_strata,
        as.integer(units)
      ),
      call. = FALSE
    )
  }
  .df_result(units - n_strata, units, n_strata, stage)
}

#' @describeIn design_df Count the PSUs a [n_cluster()] or [prec_cluster()]
#'   allocation selects. Unstratified, so the count is one design.
#' @export
design_df.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  n_psu <- (x$operational$n %||% x$n)[[1L]]
  .df_result(n_psu - 1, n_psu, 1L, "psu")
}

#' @describeIn design_df Count a [n_alloc()] allocation, at the PSU stage
#'   when the frame carried cluster columns and at the element stage
#'   otherwise, or the units of a single [n_prop()] or [n_mean()] size.
#' @export
design_df.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  if (identical(x$type, "alloc") && !is.null(x$params$psu) &&
      any(x$psu$certainty)) {
    stop(
      paste("design_df() does not yet support allocations containing certainty PSUs:",
            "their within-PSU variance contributions need an explicit",
            "degrees-of-freedom convention"),
      call. = FALSE
    )
  }
  if (identical(x$type, "multi")) {
    stop(
      "a multi-indicator result sizes several indicators against one design; count that design with design_df() on the allocation or pass the counts by name",
      call. = FALSE
    )
  }
  if (!identical(x$type, "alloc")) {
    units <- ceiling(x$n)
    return(.df_result(units - 1, units, 1L, "element"))
  }
  parts <- .df_alloc_parts(x$detail)
  .df_stratified(parts, x, .df_domain_idx(x))
}

#' @describeIn design_df Count a [strata_bound()] stratification, whose
#'   strata are element-stage and whose take-all strata are censuses.
#' @export
design_df.svyplan_strata <- function(x, ...) {
  .check_unused_dots(...)
  tab <- x$strata
  if (is.null(tab) || is.null(tab[["n"]])) {
    stop(
      "this stratification carries no allocation to count; supply 'n' or 'cv' to strata_bound()",
      call. = FALSE
    )
  }
  parts <- list(
    stratum = as.character(tab$stratum),
    units = as.numeric(tab[["n_int"]] %||% tab$n),
    census = .df_take_all(tab),
    stage = "element"
  )
  .df_stratified(parts, x, NULL)
}

#' @describeIn design_df Count the phase-2 units a [n_twophase()]
#'   allocation measures. A take-all phase-2 stratum still carries the
#'   phase-1 sampling variance, so it contributes as any other does.
#' @export
design_df.svyplan_twophase <- function(x, ...) {
  .check_unused_dots(...)
  detail <- x$detail
  parts <- list(
    stratum = as.character(detail$stratum),
    units = as.numeric(detail[["n_int"]] %||% detail$n_issued),
    census = rep(FALSE, nrow(detail)),
    stage = "element"
  )
  .df_stratified(parts, x, NULL)
}

#' Stratum labels, whole-unit counts, and censuses of an allocation
#'
#' The stage a df is counted at is decided by the allocation itself: a
#' frame carrying cluster columns is allocated in PSUs and its df counts
#' those, while an element allocation counts units. `take_all` means
#' different things in the two cases, which is why it is resolved here
#' rather than by the caller.
#' @keywords internal
#' @noRd
.df_alloc_parts <- function(detail) {
  cluster <- any(c("n_psu_draw", "n_psu_int", "n_psu") %in% names(detail))
  if (cluster) {
    units <- as.numeric(detail[["n_psu_draw"]] %||%
                          detail[["n_psu_int"]] %||% ceiling(detail$n_psu))
    # take_all is an element-level census inside the stratum and leaves the
    # PSU stage sampling, so no stratum drops out at this stage.
    census <- rep(FALSE, nrow(detail))
  } else {
    units <- as.numeric(detail[["n_int"]] %||% ceiling(detail$n))
    census <- .df_take_all(detail)
  }
  list(
    stratum = as.character(detail$stratum),
    units = units,
    census = census,
    stage = if (cluster) "psu" else "element"
  )
}

#' Take-all flags of a stratum table, absent meaning none
#' @keywords internal
#' @noRd
.df_take_all <- function(tab) {
  flag <- tab[["take_all"]]
  if (is.null(flag)) rep(FALSE, nrow(tab)) else as.logical(flag)
}

#' The domain-to-stratum mapping an allocation carries, if any
#'
#' Stored by [n_alloc()] as row indices into its own detail table, which is
#' the only exact statement of which strata a domain is made of. The domain
#' identifier columns do not survive into the detail table on their own.
#' @keywords internal
#' @noRd
.df_domain_idx <- function(x) {
  idx <- x$params$domain_idx
  if (is.null(idx) || length(idx) == 0L) NULL else idx
}

#' Assemble a stratified count into its result
#'
#' The scalar, the per-stratum table and the per-domain table are all the
#' same arithmetic over different sets of strata, so they are formed once
#' here: a contributing stratum gives up one degree of freedom to its own
#' mean, and a stratum carrying no sampling variance gives nothing and
#' takes nothing.
#' @keywords internal
#' @noRd
.df_stratified <- function(parts, x, domain_idx) {
  units <- parts$units
  census <- parts$census
  contrib <- !census
  status <- ifelse(census, "census",
                   ifelse(units <= 1, "singleton", "ok"))
  per_stratum <- ifelse(contrib, pmax(units - 1, 0), 0)

  .warn_singleton_strata(parts$stratum[status == "singleton"], parts$stage)

  strata_tab <- data.frame(
    stratum = parts$stratum,
    n_units = units,
    df = per_stratum,
    .status = status,
    stringsAsFactors = FALSE
  )

  domains_tab <- NULL
  if (!is.null(domain_idx)) {
    rows <- lapply(seq_along(domain_idx), function(i) {
      idx <- domain_idx[[i]]
      keep <- idx[contrib[idx]]
      row <- .df_domain_values(x, i)
      row$.domain <- names(domain_idx)[i]
      row$.n_units <- sum(units[keep])
      row$.df <- max(sum(units[keep]) - length(keep), 0)
      row$.status <- if (length(keep) == 0L) {
        "census"
      } else if (row$.df == 0) {
        "singleton"
      } else {
        "ok"
      }
      row
    })
    domains_tab <- do.call(rbind, rows)
    rownames(domains_tab) <- NULL
  }

  .df_result(
    sum(per_stratum), sum(units[contrib]), sum(contrib), parts$stage,
    strata = strata_tab, domains = domains_tab
  )
}

#' Identifier columns of one domain, as a one-row frame
#' @keywords internal
#' @noRd
.df_domain_values <- function(x, i) {
  dom <- x$domains
  cols <- x$params$domain_cols
  if (is.null(dom) || length(cols) == 0L || !all(cols %in% names(dom))) {
    return(data.frame(row.names = 1L)[, 0L, drop = FALSE])
  }
  dom[i, cols, drop = FALSE]
}

#' Warn once, naming the strata that hold a single primary unit
#'
#' A stratum with one PSU supports no within-stratum variance estimate at
#' all, which is a defect in the design rather than a small number to weigh
#' up, and it is the only case here that warrants a condition.
#' @keywords internal
#' @noRd
.warn_singleton_strata <- function(labels, stage = "psu") {
  if (length(labels) == 0L) {
    return(invisible(NULL))
  }
  many <- length(labels) > 1L
  warning(
    sprintf(
      "%s %s %s a single %s, so no within-stratum variance can be estimated there. Collapse %s with a neighbor",
      if (many) "strata" else "stratum",
      paste0("'", labels, "'", collapse = ", "),
      if (many) "each hold" else "holds",
      if (identical(stage, "psu")) "PSU" else "unit",
      if (many) "them" else "it"
    ),
    call. = FALSE
  )
  invisible(NULL)
}

#' Warn when a whole design rests on a single primary unit
#'
#' The unstratified counterpart of a singleton stratum: one PSU leaves no
#' between-PSU variance to estimate anywhere in the design, so the plan is
#' not one a variance estimator can be run on.
#' @keywords internal
#' @noRd
.warn_single_psu <- function(n_psu) {
  if (length(n_psu) != 1L || !is.finite(n_psu) || n_psu > 1) {
    return(invisible(NULL))
  }
  warning(
    "the design selects a single PSU, so no between-PSU variance can be estimated. Raise the budget, the PSU count, or the number of stages",
    call. = FALSE
  )
  invisible(NULL)
}

#' Finish a count into a svyplan_df, refusing an exhausted design
#' @keywords internal
#' @noRd
.df_result <- function(df, n_units, n_strata, stage, strata = NULL,
                       domains = NULL) {
  if (df < 0) {
    stop(
      "the design has fewer sampled units than strata, so its variance estimator has no degrees of freedom",
      call. = FALSE
    )
  }
  .new_svyplan_df(df, n_units, n_strata, stage, strata = strata,
                  domains = domains)
}
