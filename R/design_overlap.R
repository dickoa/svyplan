#' Sample overlap a rotation produces
#'
#' Compute the fraction of the units a rotation puts in sample that it
#' carries from one occasion to another, at every lag the rotation reaches.
#' This is the overlap of the **issued** sample, a property of the rotation
#' alone. It is the `overlap` that [n_change()], [prec_change()],
#' [power_mean()], [power_prop()] and [power_did()] take when response is
#' complete, replacing a number the planner would otherwise have to know.
#' See Details for what a response rate below 1 does to that reading.
#'
#' @param rotation A `svyplan_rotation` from [design_rotation()], or anything
#'   that function accepts, which is wrapped for you. The pattern vocabulary,
#'   the two notations and every refusal live on [design_rotation()].
#' @param max_lag Largest lag to report. Defaults to one less than the life,
#'   which is the largest lag at which any overlap is possible.
#'
#' @return A `svyplan_overlap` object: a numeric vector of issued-sample
#'   overlap fractions named by lag, so `x[1]` is the consecutive-occasion
#'   overlap and `x[12]` the overlap a year apart on monthly occasions.
#'   Subsetting returns a plain number, which is what the `overlap`
#'   arguments take, and arithmetic returns bare numerics. It also carries:
#' \describe{
#'   \item{`$shared`}{Units in sample on both occasions, at each lag.}
#'   \item{`$rotation`}{The [design_rotation()] it was computed from, which
#'     carries the life, the per-occasion takes and the units in sample at
#'     any one occasion.}
#' }
#'
#' @details
#' One cohort enters at every occasion, so at a steady state the cohorts
#' alive span every stage of the life and the sample at occasion \eqn{t} is
#' \eqn{\sum_s w_s}, writing \eqn{w} for the rotation. Units in sample at
#' both \eqn{t} and \eqn{t + m} are those whose cohort is in sample at
#' stages \eqn{s} and \eqn{s + m}, so
#'
#' \deqn{\mathrm{shared}(m) = \sum_{s = 1}^{L - m} \min(w_s, w_{s + m}),
#'       \qquad \mathrm{overlap}(m) = \mathrm{shared}(m) / \sum_s w_s.}{
#'       shared(m) = sum_(s = 1)^(L - m) min(w_s, w_(s + m)),
#'       and overlap(m) = shared(m) / sum_s w_s.}
#'
#' The `min` is the count when the smaller take at the two occasions is a
#' subset of the larger, which is what attrition and planned subsampling both
#' give, including a design that subsamples an interim wave and returns to
#' the whole cohort afterwards. Nesting is an assumption about the design and
#' is not identifiable from the takes. Where the two occasions measure
#' partly disjoint subsamples of the cohort, `min` is the largest overlap
#' the takes admit rather than the overlap the design has.
#'
#' These are the overlaps of the design at a steady state, which a rotation
#' reaches once every stage of the life is represented at one occasion. Which
#' occasion that is depends on how the design is launched, and the overlaps
#' themselves do not: `shared(m)` sums over stages, and entry dates do not
#' enter it. A design that recruits one cohort an occasion reaches the steady
#' state at the occasion that spans the life, the overlaps before it being
#' those of a design no cohort has yet aged out of. A design whose first
#' occasion contains launch components for every possible remaining life
#' length is at the steady state from the first occasion, and these figures
#' hold throughout. For an unbroken equal-take life, those are equal panels
#' planned for lives from the full life down to one occasion, all at their
#' first interview.
#'
#' Covering every stage is the exact requirement, and for a life with a gap
#' it is more than one cohort per remaining length. Under `"1-1-0-0-1-1"` the
#' stages that hold the steady state include two that are out of sample, so
#' the cohorts starting there are selected but not interviewed until they
#' reach their first in-sample stage. Splitting the first occasion's
#' interviewed sample alone puts six cohorts in sample where the design holds
#' four, and gives a consecutive overlap of 0.83 against the rotation's 0.50.
#' [plot.svyplan_overlap()] draws the one-cohort-an-occasion launch and marks
#' the occasion the steady state begins.
#'
#' **The overlap has one direction here, and that is a property of the
#' steady state rather than a simplification.** Elsewhere in the package
#' `overlap` is `n12 / n1` and is not symmetric when the occasions differ in
#' size. Under a stationary rotation every occasion holds the same mix of
#' life stages, so \eqn{n_t} is the same at every occasion and the two
#' directions coincide. A design whose cohorts differ in size by entry date
#' has no steady state and is outside what this function describes.
#'
#' ## Issued overlap, and the overlap a variance formula reads
#'
#' What a rotation fixes is which units are **issued** at both occasions.
#' The `overlap` argument of [n_change()], [prec_change()] and the power
#' family is the fraction of the first occasion's **responding** sample
#' measured again, those functions netting each occasion down by `resp_rate`
#' before the covariance forms. The two are the same number at
#' `resp_rate = 1`, and this function's result can be passed straight
#' through there.
#'
#' Below full response they are not, and converting one into the other needs
#' an assumption about how response persists across occasions, which neither
#' function makes. Under response independent between occasions at rate
#' \eqn{r}, an issued overlap \eqn{f} leaves a responding overlap of about
#' \eqn{fr}: a shared unit responds at both occasions with probability
#' \eqn{r^2}, against the \eqn{r} of issued units responding at the first. A
#' panel whose wave-1 respondents are much likelier to respond again sits
#' above that, reaching \eqn{f} itself when response persists perfectly.
#' Pass the figure you expect among respondents, and say which assumption
#' produced it.
#'
#' ## What it does not give you
#'
#' The correlation between occasions. `overlap` is a property of the
#' rotation and is fixed once the design is declared, whereas `overlap_cor`
#' is a
#' property of the variable being measured and has to come from a previous
#' round of the same survey. Both enter the variance of a change, and only
#' their product buys precision, so a rotation alone does not say what a
#' rotation is worth.
#'
#' @references
#' U.S. Census Bureau. *Current Population Survey: Design and Methodology*,
#' Technical Paper 77. The 4-8-4 rotation and its 75% and 50% overlaps.
#'
#' Lynn, P. (2012). *Longitudinal Survey Methods for the Household Finance
#' and Consumption Survey*. Report to the European Central Bank. The
#' one-flag-per-occasion notation, and the 1-1-0-0-1-1 design whose lag
#' profile the examples reproduce.
#'
#' @seealso [design_rotation()] for the pattern this reads,
#'   [plot.svyplan_overlap()] for the rotation chart, [n_change()] and
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
#' # Under response independent between occasions at 80 percent, the overlap
#' # among respondents is the issued one scaled by the response rate
#' n_change(p = c(0.30, 0.36), moe = 0.02, resp_rate = 0.8,
#'          overlap = cps[1] * 0.8, overlap_cor = 0.5)
#'
#' # One flag per occasion, the other notation in use: in for two, out for
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
#' # The chart of the rotation, one row per cohort
#' plot(design_overlap("1-1-0-0-1-1"))
#'
#' @export
design_overlap <- function(rotation, max_lag = NULL) {
  rotation <- design_rotation(rotation)
  w <- as.double(rotation)
  life <- length(w)
  n_occasion <- sum(w)

  if (is.null(max_lag)) {
    max_lag <- max(life - 1L, 1L)
  }
  if (!is.numeric(max_lag) || length(max_lag) != 1L || anyNA(max_lag) ||
      max_lag < 1 || max_lag != trunc(max_lag)) {
    stop("'max_lag' must be a whole number >= 1", call. = FALSE)
  }
  max_lag <- as.integer(max_lag)

  lags <- seq_len(max_lag)
  shared <- vapply(lags, function(m) {
    if (m >= life) {
      return(0)
    }
    s <- seq_len(life - m)
    sum(pmin(w[s], w[s + m]))
  }, numeric(1L))

  structure(
    shared / n_occasion,
    names = as.character(lags),
    shared = shared,
    rotation = rotation,
    class = c("svyplan_overlap", "numeric")
  )
}
