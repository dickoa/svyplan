#' Sample overlap a rotation produces
#'
#' Compute the fraction of the units a rotation puts in sample that it
#' carries from one period to another, at every lag the rotation reaches.
#' This is the overlap of the **issued** sample, a property of the rotation
#' alone. It is the `overlap` that [n_change()], [prec_change()],
#' [power_mean()], [power_prop()] and [power_did()] take when response is
#' complete, replacing a number the planner would otherwise have to know.
#' See Details for what a response rate below 1 does to that reading.
#'
#' @param rotation A `svyplan_rotation` from [design_rotation()], or anything
#'   that function accepts, which is wrapped for you. The pattern vocabulary,
#'   the two notations and every refusal live on [design_rotation()]. A
#'   rotating [n_panel()] plan gives the overlap among respondents instead.
#'   See the section on respondent overlap in Details.
#' @param max_lag Largest lag to report. Defaults to one less than the life,
#'   which is the largest lag at which any overlap is possible.
#' @param window Optional number of consecutive periods, such as 4 for a
#'   year of quarterly periods. When given, the result also carries the
#'   mean number of periods a unit is in sample within such a window. See
#'   the section on visits in Details. Not available for a panel plan, the
#'   visits being those of issued units.
#'
#' @return A `svyplan_overlap` object: a numeric vector of issued-sample
#'   overlap fractions named by lag, so `x[1]` is the consecutive-period
#'   overlap and `x[12]` the overlap a year apart on monthly periods.
#'   Subsetting returns a plain number, which is what the `overlap`
#'   arguments take, and arithmetic returns bare numerics. It also carries:
#' \describe{
#'   \item{`$shared`}{Units in sample in both periods, at each lag. For a
#'     panel plan, expected respondents per entrant.}
#'   \item{`$in_sample`}{Units in sample at one period, or for a panel plan
#'     expected respondents per entrant.}
#'   \item{`$basis`}{`"issued"` for a rotation, `"respondent"` for a panel
#'     plan.}
#'   \item{`$rotation`}{The [design_rotation()] it was computed from, which
#'     carries the life, the per-period takes and the units in sample at
#'     any one period.}
#'   \item{`$window`, `$visits`}{The window length and the mean number of
#'     periods a unit is in sample within it, `NA` when `window` is not
#'     given.}
#' }
#'
#' @details
#' One cohort enters at every period, so once the full sample is in the
#' field the cohorts alive span every stage of the life and the sample at period \eqn{t} is
#' \eqn{\sum_s w_s}, writing \eqn{w} for the rotation. Units in sample at
#' both \eqn{t} and \eqn{t + m} are those whose cohort is in sample at
#' stages \eqn{s} and \eqn{s + m}, so
#'
#' \deqn{\mathrm{shared}(m) = \sum_{s = 1}^{L - m} \min(w_s, w_{s + m}),
#'       \qquad \mathrm{overlap}(m) = \mathrm{shared}(m) / \sum_s w_s.}{
#'       shared(m) = sum_(s = 1)^(L - m) min(w_s, w_(s + m)),
#'       and overlap(m) = shared(m) / sum_s w_s.}
#'
#' The `min` is the count when the smaller take at the two periods is a
#' subset of the larger, which is what attrition and planned subsampling both
#' give, including a design that subsamples an interim wave and returns to
#' the whole cohort afterward. Nesting is an assumption about the design and
#' is not identifiable from the takes. Where the two periods measure
#' partly disjoint subsamples of the cohort, `min` is the largest overlap
#' the takes admit rather than the overlap the design has.
#'
#' These are the overlaps of the design with the full sample in the field,
#' which a rotation reaches once every stage of the life is represented at
#' one period. Which period that is depends on how the design is started,
#' and the overlaps themselves do not: `shared(m)` sums over stages, and
#' entry dates do not enter it. A design that recruits one cohort a period
#' reaches the full sample at the period that spans the life, the overlaps
#' before it being those of a design no cohort has yet aged out of. A design
#' whose first period selects a start-up cohort for every possible remaining
#' life length has the full sample from the first period, and these figures
#' hold throughout. For an unbroken equal-take life, those are equal panels
#' planned for lives from the full life down to one period, all at their
#' first interview.
#'
#' Covering every stage is the exact requirement, and for a life with a gap
#' it is more than one cohort per remaining length. Under `"1-1-0-0-1-1"` the
#' stages a full sample covers include two that are out of sample, so
#' the cohorts starting there are selected but not interviewed until they
#' reach their first in-sample stage. Splitting the first period's
#' interviewed sample alone puts six cohorts in sample where the design holds
#' four, and gives a consecutive overlap of 0.83 against the rotation's 0.50.
#' [plot.svyplan_rotation()] draws either start and marks the first
#' period with the full sample.
#'
#' **The overlap has one direction here, and that is a property of the
#' full sample rather than a simplification.** Elsewhere in the package
#' `overlap` is `n12 / n1` and is not symmetric when the periods differ in
#' size. Under a stationary rotation every period holds the same mix of
#' life stages, so \eqn{n_t} is the same at every period and the two
#' directions coincide. A design whose cohorts differ in size by entry date
#' never holds the same mix at every period and is outside what this
#' function describes.
#'
#' ## Issued overlap, and the overlap a variance formula reads
#'
#' What a rotation fixes is which units are **issued** in both periods.
#' The `overlap` argument of [n_change()], [prec_change()] and the power
#' family is the fraction of the first period's **responding** sample
#' measured again, those functions netting each period down by `resp_rate`
#' before the covariance forms. The two are the same number at
#' `resp_rate = 1`, and this function's result can be passed straight
#' through there. Those functions call the two times compared occasions,
#' the textbook word for any time a survey is taken, whatever its level
#' (Cochran, 1977, sections 12.10 to 12.12). For an overlap read from a
#' rotation, an occasion is a period of that rotation and the lag counts
#' periods.
#'
#' Below full response they are not, and converting one into the other needs
#' an assumption about how response persists across periods, which neither
#' function makes. Under response independent between periods at rate
#' \eqn{r}, an issued overlap \eqn{f} leaves a responding overlap of about
#' \eqn{fr}: a shared unit responds in both periods with probability
#' \eqn{r^2}, against the \eqn{r} of issued units responding at the first. A
#' panel whose wave-1 respondents are much likelier to respond again sits
#' above that, reaching \eqn{f} itself when response persists perfectly.
#' Pass the figure you expect among respondents, and say which assumption
#' produced it.
#'
#' ## Respondent overlap from a panel plan
#'
#' A rotating [n_panel()] plan states the second assumption, loss that is
#' permanent: a unit that stops responding does not return. Then a unit
#' responding at the later of two stages of its life responded at the
#' earlier one, and with \eqn{q_s} the probability of still responding at
#' stage \eqn{s} (zero out of sample),
#'
#' \deqn{\mathrm{overlap}(m) = \sum_{s:\, w_s > 0,\, w_{s + m} > 0}
#'       q_{s + m} \Big/ \sum_s q_s.}{
#'       overlap(m) = sum(q_(s + m), over s with w_s > 0 and w_(s + m) > 0)
#'       / sum_s q_s.}
#'
#' Recruitment response and eligibility multiply every \eqn{q_s} and
#' cancel, so at retention 1 this is the issued overlap whatever `resp_rate`
#' is. Retention below 1 lowers it, and more at long lags. The result can be
#' passed as the `overlap` of [n_pooled()] and [prec_pooled()] at any
#' `resp_rate`, and one lag of it to [n_change()] and [prec_change()]. It
#' describes equal cohorts at equilibrium and does not model response that
#' recovers after a loss, nor units joining a sampled dwelling.
#'
#' ## Visits per unit in a window
#'
#' A window of \eqn{W} consecutive periods holds \eqn{W \sum_s w_s}
#' interviews with the full sample in the field. The units behind them are
#' fewer, since a
#' unit in sample at several periods of the window is counted once. A
#' cohort sees a run of at most \eqn{W} consecutive stages of its life
#' inside the window, and with its units nested across stages, as for the
#' overlap, it contributes the largest take in that run. Summing over the
#' cohorts the window touches gives the distinct units, and
#'
#' \deqn{\mathrm{visits}(W) = W \sum_s w_s \Big/ \sum_k
#'       \max_{s \in [k, k + W - 1] \cap [1, L]} w_s,}{
#'       visits(W) = W sum_s w_s / sum_k max(w_s, s in [k, k + W - 1] and
#'       in [1, L]),}
#'
#' with \eqn{k} running over the \eqn{L + W - 1} cohorts whose life meets
#' the window. Under `"2-2-2"` on quarterly periods, a year holds 16
#' interviews of 9 distinct panels, so 16/9 = 1.78.
#'
#' The figure is a property of the rotation. Some labor force survey
#' methodologies use it to scale the cluster size when a yearly estimate
#' pools quarterly interviews, which folds the correlation between a unit's
#' repeated interviews into the intraclass correlation. [prec_pooled()]
#' treats the average over periods with that correlation stated
#' separately, as `overlap_cor`.
#'
#' ## What it does not give you
#'
#' The correlation between periods. `overlap` is a property of the
#' rotation and is fixed once the design is declared, whereas `overlap_cor`
#' is a
#' property of the variable being measured and has to come from a previous
#' round of the same survey. Both enter the variance of a change, and only
#' their product buys precision, so a rotation alone does not say what a
#' rotation is worth.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley. Sections
#' 12.10 to 12.12, sampling on successive occasions.
#'
#' U.S. Census Bureau and U.S. Bureau of Labor Statistics (2019). *Current
#' Population Survey: Design and Methodology*. Technical Paper 77. The 4-8-4
#' rotation and its 75% and 50% overlaps.
#'
#' Lynn, P. (2012). *Longitudinal Survey Methods for the Household Finances
#' and Consumption Survey*. Report prepared for the European Central Bank.
#' The one-flag-per-period notation, and the 1-1-0-0-1-1 design whose lag
#' profile the examples reproduce.
#'
#' @seealso [design_rotation()] for the pattern this reads,
#'   [plot.svyplan_overlap()] for the bar chart by lag,
#'   [plot.svyplan_rotation()] for the rotation chart, [n_panel()] for the
#'   plan a respondent overlap comes from, [n_change()] and
#'   [prec_change()], which take the result as their `overlap`, and
#'   [design_effect()] and [design_df()] for the other quantities a planned
#'   design determines.
#'
#' @family repeated survey functions
#'
#' @examples
#' # CPS 4-8-4: 75 percent month to month, 50 percent a year apart
#' cps <- design_overlap(design_rotation("4-8-4"))
#' cps[1]
#' cps[12]
#'
#' # A string is wrapped for you, which is the one-liner form
#' design_overlap("5")[1]
#'
#' # Feed it straight into a change requirement, response being complete here
#' n_change(p = c(0.30, 0.36), moe = 0.02,
#'          overlap = cps[1], overlap_cor = 0.5)
#'
#' # The annual lag is the one a year-on-year change uses
#' n_change(p = c(0.30, 0.36), moe = 0.02,
#'          overlap = cps[12], overlap_cor = 0.5)
#'
#' # Under response independent between periods at 80 percent, the overlap
#' # among respondents is the issued one scaled by the response rate
#' n_change(p = c(0.30, 0.36), moe = 0.02, resp_rate = 0.8,
#'          overlap = cps[1] * 0.8, overlap_cor = 0.5)
#'
#' # Under permanent loss, the panel plan gives it: recruitment response
#' # cancels and retention lowers the overlap, more at the annual lag
#' plan <- n_panel(n_prop(p = 0.3, moe = 0.02), resp_rate = 0.8,
#'                 retention = c(0.97, 0.97, 0.97, 0.85, 0.97, 0.97, 0.97),
#'                 design = "rotating", rotation = "4-8-4")
#' resp <- design_overlap(plan)
#' c(resp[1], resp[12])
#' n_change(p = c(0.30, 0.36), moe = 0.02, resp_rate = 0.8,
#'          overlap = resp[1], overlap_cor = 0.5)
#'
#' # One flag per period, the other notation in use: in for two, out for
#' # two, in for two. Change is estimable at lags 1, 3, 4 and 5 but not 2
#' design_overlap("1-1-0-0-1-1")
#'
#' # The same design as spell lengths
#' identical(as.double(design_overlap("2-2-2")),
#'           as.double(design_overlap("1-1-0-0-1-1")))
#'
#' # An explicit take vector, halving the take at later waves
#' design_overlap(c(1, 1, 0.5, 0.5))
#'
#' # Interviews per unit over a year of quarterly periods
#' lfs <- design_overlap("2-2-2", window = 4)
#' lfs$visits
#'
#' # An annual estimate that counts each unit once rests on the distinct
#' # units of the year: four quarters of 1600 interviews, 16/9 per unit
#' distinct <- 4 * 1600 / lfs$visits
#' prec_prop(p = 0.3, n = distinct)
#'
#' # The overlap at each lag, and the rotation chart it comes from
#' ov <- design_overlap("1-1-0-0-1-1")
#' plot(ov)
#' plot(ov$rotation)
#'
#' @export
design_overlap <- function(rotation, max_lag = NULL, window = NULL) {
  q <- NULL
  if (inherits(rotation, "svyplan_panel")) {
    q <- .overlap_panel_response(rotation, window)
    rotation <- rotation$rotation
  }
  rotation <- design_rotation(rotation)
  w <- as.double(rotation)
  life <- length(w)

  if (is.null(max_lag)) {
    max_lag <- max(life - 1L, 1L)
  }
  if (!is.numeric(max_lag) || length(max_lag) != 1L || anyNA(max_lag) ||
      max_lag < 1 || max_lag != trunc(max_lag)) {
    stop("'max_lag' must be a whole number >= 1", call. = FALSE)
  }
  max_lag <- as.integer(max_lag)

  lags <- seq_len(max_lag)
  prof <- .overlap_profile(w, lags, q)

  visits <- NA_real_
  if (!is.null(window)) {
    if (!is.numeric(window) || length(window) != 1L || anyNA(window) ||
        !is.finite(window) || window < 1 || window != trunc(window)) {
      stop("'window' must be a whole number >= 1", call. = FALSE)
    }
    window <- as.integer(window)
    visits <- .overlap_visits(w, window)
  } else {
    window <- NA_integer_
  }

  structure(
    prof$shared / prof$in_sample,
    names = as.character(lags),
    shared = prof$shared,
    rotation = rotation,
    window = window,
    visits = visits,
    basis = if (is.null(q)) "issued" else "respondent",
    response = q,
    in_sample = prof$in_sample,
    class = c("svyplan_overlap", "numeric")
  )
}

#' Shared and in-sample counts at each lag
#'
#' Without `q` these are issued units, with nested takes. With `q`, the
#' probability of still responding at each wave, they are expected
#' respondents per entrant. Loss being permanent, a unit responding at the
#' later of two stages responded at the earlier one, so the shared
#' respondents are those of the later stage.
#' @keywords internal
#' @noRd
.overlap_profile <- function(w, lags, q = NULL) {
  life <- length(w)
  if (is.null(q)) {
    shared <- vapply(lags, function(m) {
      if (m >= life) {
        return(0)
      }
      s <- seq_len(life - m)
      sum(pmin(w[s], w[s + m]))
    }, numeric(1L))
    return(list(shared = shared, in_sample = sum(w)))
  }
  stage_q <- numeric(life)
  stage_q[w > 0] <- q
  shared <- vapply(lags, function(m) {
    if (m >= life) {
      return(0)
    }
    s <- seq_len(life - m)
    s <- s[w[s] > 0 & w[s + m] > 0]
    sum(stage_q[s + m])
  }, numeric(1L))
  list(shared = shared, in_sample = sum(stage_q))
}

#' Response by wave of a rotating panel plan, for its respondent overlap
#' @keywords internal
#' @noRd
.overlap_panel_response <- function(panel, window) {
  if (!identical(panel$design, "rotating")) {
    stop(
      "a fixed panel follows one cohort, so its waves differ in size and the overlap between two of them is the later wave's responding count over the earlier one's, read from its $waves; design_overlap() takes a rotating panel",
      call. = FALSE
    )
  }
  if (!is.null(window)) {
    stop(
      "'window' counts visits per issued unit, a property of the rotation; pass the rotation, design_overlap(panel$rotation, window = ...)",
      call. = FALSE
    )
  }
  panel$waves$q
}

#' Mean periods in sample per distinct unit within a window
#'
#' Each cohort meeting the window sees a run of consecutive stages, and its
#' units being nested across stages, it contributes the largest take in the
#' run.
#' @keywords internal
#' @noRd
.overlap_visits <- function(w, window) {
  life <- length(w)
  distinct <- vapply(seq(2L - window, life), function(k) {
    s <- seq(max(k, 1L), min(k + window - 1L, life))
    max(w[s])
  }, numeric(1L))
  window * sum(w) / sum(distinct)
}
