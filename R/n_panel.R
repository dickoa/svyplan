#' Recruitment size for a panel that loses units between waves
#'
#' Compute how many units to recruit so that a panel still delivers a
#' required responding sample after several waves of attrition, either for a
#' fixed panel followed to a chosen wave or for a rotating panel at steady
#' state. The requirement and the estimand both come from an existing sizing
#' or precision result, so the panel arithmetic is stated once and the
#' design it serves is stated where it always was.
#'
#' @param target A `svyplan_n` or `svyplan_prec` result for a mean or a
#'   proportion, from [n_mean()], [n_prop()], [prec_mean()] or
#'   [prec_prop()]. It supplies two things: the responding sample the panel
#'   must deliver, and the estimand whose precision is reported at every
#'   wave. Its own `resp_rate` is removed first, so the requirement is a
#'   count of respondents and no response is counted twice. The panel's
#'   `resp_rate` is the only recruitment response that reaches the answer.
#'   Where the two differ the target's is reported as unused, at the call
#'   and again in `summary()`. Clustered, allocation, multi-indicator,
#'   multi-domain, change and two-phase results are refused, their
#'   stage-specific and occasion-specific sizes not being one responding
#'   count.
#' @param retention Conditional retention, one value per wave transition,
#'   each in (0, 1\]. `retention[j]` is the share of wave `j` respondents who
#'   respond again at wave `j + 1`. Its length plus one is the number of
#'   waves in a unit's life, so a five-wave panel takes four values.
#' @param ... Additional arguments are not supported and produce an error.
#' @param resp_rate Response rate at recruitment, wave 1, in (0, 1\], `1` by
#'   default. It is separate from `retention` because the first wave is where
#'   most of a panel's loss happens, and an average rate spread over the waves
#'   would under-issue. See the nonresponse section of [svyplan-package] for
#'   what this adjustment does and does not claim.
#' @param design `"fixed"` for one cohort followed across its waves, or
#'   `"rotating"` for equal cohorts entering every occasion and leaving at
#'   the end of their life. The two return different quantities, described
#'   under Value. This argument names the panel type. It never accepts a
#'   `survey.design` or a [svyplan()] object.
#' @param target_wave The wave at which the target must be met, defaulting
#'   to the last. Sizing to reach a precision at wave 3 of a five-wave panel
#'   is a legitimate request, and waves 4 and 5 are still reported. It must
#'   be absent for a rotating design, whose target is an occasion rather
#'   than a wave.
#' @param assurance Probability in (0, 1), or `NULL` (default). The
#'   recruitment above is an expected-value calculation, which does not
#'   guarantee that the target will be reached. The shortfall probability
#'   depends on the response distribution and integer rounding. Supplying a level reports, next
#'   to it, the recruitment for which the required respondents arrive with
#'   at least that probability. See Details.
#'
#' @param start How a rotating design's cohorts are brought in, or `NULL`
#'   (default). `"gradual"` recruits one cohort an occasion, so the design
#'   fills up over a life. `"immediate"` splits the first occasion into equal
#'   panels planned for life lengths from the full life down to one occasion,
#'   all beginning at wave 1, so it is full at once. Supplying either reports
#'   what the design delivers at each occasion until it settles. The default
#'   plans no launch, which is what the recruitment above describes either
#'   way. Not available for a fixed panel, which recruits one cohort.
#'
#' @return A `svyplan_panel` object. [prec_panel()] returns the same class,
#'   `$solved` naming the direction it was computed in. Fields:
#' \describe{
#'   \item{`n_issued`}{Fixed panels only. Units to issue to the one cohort.}
#'   \item{`n_entrants`, `n_in_sample`, `n_cohorts`}{Rotating panels only.
#'     Entrants per occasion once the design is running, the units the
#'     design holds across every live cohort, and how many cohorts that is.
#'     The first two are different budget lines and are named separately for
#'     that reason: `n_in_sample` is also the cumulative recruitment that
#'     reaching a steady state takes, whether the cohorts are taken on at
#'     once or phased in over the first `n_cohorts` occasions.}
#'   \item{`n_target`}{The responding sample the target requires.}
#'   \item{`n_resp`}{The responding sample the design delivers where the
#'     target is stated: at `target_wave` for a fixed panel, pooled across
#'     the live cohorts for a rotating one. It is what `se`, `moe` and `cv`
#'     are computed on.}
#'   \item{`n_assured`, `assured_feasible`}{Recruitment meeting the target
#'     with probability at least `assurance`, when a level was given, and
#'     whether a finite population can supply it. See Details.}
#'   \item{`se`, `moe`, `cv`, `rmoe`}{Precision of the embedded estimand at
#'     `n_resp`. For a fixed panel sized in the usual direction these
#'     reproduce the target's own precision exactly.}
#'   \item{`waves`}{One row per wave of a unit's life: the conditional
#'     `retention` into it, the cumulative response probability `q`, the
#'     share `loss_share` of the whole life's loss that happens at it,
#'     the expected respondents `n_resp`, and the precision an estimate
#'     using only those respondents would have. For a rotating design the
#'     rows are the cohorts alive at one occasion and `n_resp` sums to the
#'     occasion's sample.}
#'   \item{`start`, `launch`, `launch_waves`}{Present when `start` was
#'     given. `launch` has one row per occasion up to the steady state and
#'     one past it: the entrants taken on, the units in sample, the expected
#'     respondents, the precision they buy, and `steady_state` marking the
#'     occasion from which the design holds one cohort at every wave. That is
#'     the composition rather than the count. Where no wave loses anyone, the
#'     counts coincide from the first occasion and the mix still does not. `launch_waves` decomposes each occasion into the
#'     waves in sample at it, which is where the reason for the precision at
#'     an occasion can be read. Both are continuous, as `waves` is.}
#'   \item{`target`}{The embedded object, unchanged.}
#' }
#'
#' @details
#' Writing \eqn{r} for `resp_rate` and \eqn{c_j} for `retention[j]`, the
#' probability that a unit approached at recruitment is still responding at
#' wave \eqn{w} is
#'
#' \deqn{q_w = r \prod_{j < w} c_j,}{q_w = r * prod_(j < w) c_j,}
#'
#' and the two designs invert it differently:
#'
#' \deqn{\text{fixed: } n_{\text{issued}} = n_{\text{target}} / q_w,
#'       \qquad
#'       \text{rotating: } n_{\text{entrants}} =
#'       n_{\text{target}} / \sum_{s} q_s.}{fixed: n_issued = n_target / q_w,
#'       and rotating: n_entrants = n_target / sum_s q_s.}
#'
#' Both are exact expected values under the stated rates. They are not the
#' same quantity and the object never writes one over the other. The fixed
#' figure is the whole issue to one cohort, while the rotating figure is
#' what enters at each occasion, with `n_in_sample` the separate line for
#' the units every live cohort holds at once.
#'
#' A rotating panel pools every live cohort into one occasion, so its
#' responding sample is \eqn{e \sum_s q_s}{e * sum_s q_s} rather than one
#' cohort's count at
#' its final wave. On the example below the two differ by a factor of five.
#'
#' ## Where a panel loses its sample
#'
#' `waves$loss_share` divides the whole life's loss across the waves. It is
#' usually concentrated at recruitment, 61 percent of it in the example
#' below, which is the argument for `resp_rate` and `retention` being
#' separate arguments rather than one average rate.
#'
#' ## Bringing a rotating design up to its steady state
#'
#' A rotating design's recruitment describes it at a steady state, which it
#' reaches once every stage of the life is represented at one occasion.
#' Getting there is a design decision, and `start` reports what each choice
#' delivers on the way.
#'
#' A **gradual** launch recruits one cohort an occasion, so the sample climbs
#' over a full life before it is the design's own. An **immediate** launch
#' splits the first occasion into equal panels with planned life lengths from
#' the full life down to one occasion. Every panel begins at wave 1, and
#' together they hold the whole sample from the first occasion.
#' [design_overlap()] gives the same mature membership-overlap profile under
#' either launch. During a gradual launch, realized overlap is higher until
#' every life stage is represented because no full set of cohorts has yet
#' rotated through. The launch tables therefore describe both the early
#' membership mix and its response and precision path.
#'
#' The point of reporting it is that an immediate launch is **not** in
#' response equilibrium at its first occasion even though it is in membership
#' equilibrium. Every unit there is at wave 1, so that occasion holds
#'
#' \deqn{R_1 - e \sum_s q_s = e \sum_s (q_1 - q_s) \ge 0,}{
#'       R_1 - e sum_s q_s = e sum_s (q_1 - q_s) >= 0,}
#'
#' at least as many respondents as the design ever holds again, and strictly
#' more as soon as any wave retains less than all of the one before. The two
#' coincide exactly when every `retention` is 1, whatever `resp_rate` is,
#' since it cancels from both sides. Where they differ the early precision is
#' temporarily better. On the example below the first occasion holds 1172
#' respondents against the design's 1001, `moe` 0.0287 against 0.0310,
#' converging down as the interview mix matures. A gradual launch approaches the
#' same figure from below. Either way the early occasions rest on a different
#' response composition from the rest of the series, which `launch_waves` is
#' there to expose and which is what nonresponse weighting has to carry.
#'
#' Both are described for a life without a break in it. A schedule that
#' leaves the sample and returns needs launch cohorts that are selected
#' before they are first interviewed, which is a longer definition than this
#' argument carries.
#'
#' ## Assurance
#'
#' For a fixed panel the assured recruitment is the smallest \eqn{g} with
#' \eqn{P(\mathrm{Binomial}(g, q_w) \ge n_{\text{target}}) \ge}{P(Binomial(g,
#' q_w) >= n_target) >=} `assurance`.
#' For a rotating panel the respondents at one occasion are a sum of
#' binomials at **different** cumulative probabilities, one per live cohort,
#' so the distribution is Poisson-binomial rather than binomial. It is
#' assembled exactly, by convolving one distribution per cohort, rather than
#' approximated by a binomial on their mean rate.
#'
#' The level is marginal. For a rotating design it holds at one occasion,
#' and consecutive occasions share cohorts, so the chance that every
#' occasion of a run clears its target is lower and is not computed here.
#'
#' A requested level can be out of reach. Where the assured recruitment
#' exceeds a finite `N`, not even a census of the frame delivers the target
#' that often. The recruitment it would take is still reported, that being
#' the number a planner needs in order to argue for a larger frame or a
#' smaller target, and `assured_feasible` is `FALSE` alongside a warning.
#' Neither the expected design nor its precision is affected, which is why
#' this is not an error: `assurance` is an addition to a plan that stands
#' without it. A warning alone would not do, though, since it does not
#' survive into the object a script reads.
#'
#' ## What the rates are, and what they are not
#'
#' This is arithmetic on declared rates. Where the rates come from is the
#' planner's problem: attrition modelling, attrition weighting and any
#' adjustment for informative dropout are outside what the package does.
#'
#' Two assumptions are worth stating because they push in opposite
#' directions. A unit lost at wave \eqn{j} is treated as lost for good, so a
#' panel whose nonrespondents return at a later wave will do better than
#' planned here and this over-issues. And the rates are one set for the
#' whole sample. Where response differs by domain, recruitment should be
#' solved separately within groups of similar attrition, which
#' over-represents the hard-to-retain groups at wave 1. Call `n_panel()`
#' once per group and hand the results on as a named vector.
#'
#' A target's constraints hold where the target is stated. Sizing to a
#' `min_cases` floor at wave 3 of a five-wave panel leaves waves 4 and 5
#' with fewer expected cases than the floor, which is a consequence of
#' `target_wave` and not a failure of the constraint.
#'
#' The `deff` of the embedded target is applied unchanged at every wave. A
#' panel adds a within-unit correlation over time that a cross-sectional
#' design effect does not describe, and that correlation works for a change
#' and against a pooled average, so a design serving either should carry a
#' design effect chosen for that estimand.
#'
#' ## One class for both directions
#'
#' `n_panel()` and [prec_panel()] both return `svyplan_panel`, unlike the
#' rest of the package where the two directions return `svyplan_n` and
#' `svyplan_prec`. A panel's recruitment and its precision move together
#' across the waves and are only readable side by side, so the object holds
#' both whichever direction produced it, and `$solved` records which. The
#' class is a sibling of `svyplan_n`, not a subtype: `$n_issued` is a count
#' of units to release, and the methods registered for an analysis sample
#' would each answer a different question about it.
#'
#' @references
#' Smith, P., Lynn, P., and Elliot, D. (2009). Sample design for
#' longitudinal surveys. In P. Lynn (ed.), *Methodology of Longitudinal
#' Surveys*, 21-33. Wiley.
#'
#' @family repeated survey functions
#' @seealso [prec_panel()] for the same design from a recruitment you
#'   already have, [design_overlap()] for the overlap a rotation produces,
#'   [design_schedule()] for turning a rotating plan into an operational
#'   schedule, [n_change()] for sizing the change between two occasions.
#'
#' @examples
#' # UK LFS: five quarterly waves, 73 percent at recruitment then high
#' # conditional retention. Issue 1815 addresses to hold 1000 at wave 5.
#' target <- n_prop(p = 0.5, moe = 0.031)
#' lfs <- n_panel(
#'   target,
#'   retention = c(0.878, 0.963, 0.936, 0.956),
#'   resp_rate = 0.728
#' )
#' lfs
#'
#' # Precision at every wave, and where the loss happens
#' lfs$waves
#'
#' # The same rates run as a rotating panel: entrants per occasion, and the
#' # larger number the five live cohorts hold between them
#' rot <- n_panel(
#'   target,
#'   retention = c(0.878, 0.963, 0.936, 0.956),
#'   resp_rate = 0.728,
#'   design = "rotating"
#' )
#' c(entrants = rot$n_entrants, in_sample = rot$n_in_sample)
#'
#' # Sizing to reach the target earlier: waves past it carry fewer units
#' n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
#'         resp_rate = 0.728, target_wave = 3)$waves$n_resp
#'
#' # What the design delivers while it is being brought up
#' n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
#'         resp_rate = 0.728, design = "rotating",
#'         start = "immediate")$launch
#'
#' # The same recruitment phased in one cohort an occasion instead
#' n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
#'         resp_rate = 0.728, design = "rotating",
#'         start = "gradual")$launch$n_resp
#'
#' # Why the first occasion of an immediate launch is the most precise one:
#' # every unit in it is at wave 1
#' imm <- n_panel(target, retention = c(0.878, 0.963, 0.936, 0.956),
#'                resp_rate = 0.728, design = "rotating", start = "immediate")
#' subset(imm$launch_waves, period <= 2)
#'
#' # Request at least a 95% probability of reaching the respondent target
#' n_panel(target, retention = c(0.9, 0.9), resp_rate = 0.8,
#'         assurance = 0.95)$n_assured
#'
#' @export
n_panel <- function(
  target,
  retention,
  ...,
  resp_rate = 1,
  design = c("fixed", "rotating"),
  target_wave = NULL,
  assurance = NULL,
  start = NULL
) {
  .check_unused_dots(...)
  design <- .check_panel_design(design)
  tgt <- .panel_target(target)
  q <- .panel_q(resp_rate, retention)
  target_wave <- .check_target_wave(target_wave, length(q), design)
  assurance <- .check_assurance(assurance)
  start <- .check_panel_start(start, design)

  n_recruit <- if (design == "fixed") {
    tgt$n_target / q[[target_wave]]
  } else {
    tgt$n_target / sum(q)
  }

  .panel_result(
    tgt, n_recruit, resp_rate, retention, q, design, target_wave,
    assurance, start, solved = "n_recruit"
  )
}

#' The responding requirement and the estimand a panel target carries
#'
#' The target is an evaluator as much as a number: `$params` and `$type`
#' are what every wave's precision is computed from. `resp_rate` is removed
#' here rather than approximately undone later, which is exact because the
#' `n_*` functions apply it last, after the finite population correction
#' and after any `min_cases` floor.
#' @keywords internal
#' @noRd
.panel_target <- function(target) {
  if (!inherits(target, c("svyplan_n", "svyplan_prec"))) {
    stop(
      "'target' must be a result from n_mean(), n_prop(), prec_mean() or prec_prop()",
      call. = FALSE
    )
  }
  if (!target$type %in% c("mean", "proportion")) {
    hint <- if (identical(target$type, "change")) {
      "; a change is measured across two occasions whose sizes the retention chain already fixes, so it is not one responding count"
    } else {
      ""
    }
    stop(
      sprintf(
        "a panel target must estimate a mean or a proportion, not '%s'%s",
        target$type, hint
      ),
      call. = FALSE
    )
  }
  if (!is.null(target$domains) || !is.null(target$indicators) ||
      !is.null(target$operational)) {
    stop(
      "a panel target must be a single indicator over one population; size the domains or the indicators separately and call n_panel() on each",
      call. = FALSE
    )
  }
  n <- if (inherits(target, "svyplan_n")) target$n else target$params$n
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("'target' must carry a single positive sample size", call. = FALSE)
  }
  list(
    object = target,
    type = target$type,
    method = target$method,
    params = target$params,
    N = target$params$N %||% Inf,
    n_target = n * (target$params$resp_rate %||% 1)
  )
}

#' A recruitment response the panel names twice
#'
#' A target sized with its own `resp_rate` carries a gross issued count, and
#' the panel reads it back to the responding requirement. Where the two
#' rates agree that is invisible and unsurprising, the panel re-applying at
#' wave 1 exactly what it removed. Where they disagree the user has named
#' two recruitment responses and only one of them is used, which is the
#' whole of the confusion and the whole of what is worth reporting.
#'
#' The condition lives here rather than at each of its two call sites so
#' that the warning and the line `summary()` prints cannot drift apart.
#' @keywords internal
#' @noRd
.panel_rate_conflict <- function(target_rr, resp_rate) {
  if (is.null(target_rr) || isTRUE(all.equal(target_rr, 1))) {
    return(NULL)
  }
  if (isTRUE(all.equal(target_rr, resp_rate))) {
    return(NULL)
  }
  target_rr
}

#' @keywords internal
#' @noRd
.warn_panel_target_rate <- function(tgt, resp_rate) {
  rr <- .panel_rate_conflict(tgt$params$resp_rate, resp_rate)
  if (is.null(rr)) {
    return(invisible(NULL))
  }
  warning(
    sprintf(
      "the target's own response rate %.3g is not used: a panel applies its own 'resp_rate' (%.3g) at wave 1, and the target is read as %s responding",
      rr, resp_rate, .fmt_count_n(round(tgt$n_target))
    ),
    call. = FALSE
  )
  invisible(NULL)
}

#' Cumulative probability of still responding, by wave
#' @keywords internal
#' @noRd
.panel_q <- function(resp_rate, retention) {
  check_resp_rate(resp_rate)
  if (!is.numeric(retention) || length(retention) < 1L || anyNA(retention) ||
      !all(is.finite(retention)) || any(retention <= 0) || any(retention > 1)) {
    stop(
      "'retention' must give one conditional retention per wave transition, each in (0, 1]",
      call. = FALSE
    )
  }
  resp_rate * cumprod(c(1, as.numeric(retention)))
}

#' Refuse a survey design where the panel type belongs
#'
#' `design` reads like the argument that takes a design object elsewhere in
#' the ecosystem, so the wrong thing arriving here is named rather than
#' being reported as an unmatched string.
#' @keywords internal
#' @noRd
.check_panel_design <- function(design) {
  if (inherits(design, c("survey.design", "survey.design2", "svyplan"))) {
    stop(
      "'design' names the panel type, \"fixed\" or \"rotating\"; it does not take a design object",
      call. = FALSE
    )
  }
  match.arg(design, c("fixed", "rotating"))
}

#' @keywords internal
#' @noRd
.check_target_wave <- function(target_wave, k, design) {
  if (identical(design, "rotating")) {
    if (!is.null(target_wave)) {
      stop(
        "'target_wave' does not apply to a rotating panel: the target is the responding sample at one occasion, pooled over every cohort alive, so there is no wave to select",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (is.null(target_wave)) {
    return(k)
  }
  if (!is.numeric(target_wave) || length(target_wave) != 1L ||
      is.na(target_wave) || target_wave != trunc(target_wave) ||
      target_wave < 1 || target_wave > k) {
    stop(
      sprintf("'target_wave' must be a whole number in 1:%d", k),
      call. = FALSE
    )
  }
  as.integer(target_wave)
}

#' The launch a rotating design is brought up on
#'
#' `NULL` is not "no launch": every rotating design has one. It means the
#' launch is not being planned here, which is the contract every result had
#' before this argument existed, so the default leaves them untouched.
#'
#' A fixed panel recruits one cohort and reaches its own size at its first
#' wave, so there is nothing for this to select between.
#' @keywords internal
#' @noRd
.check_panel_start <- function(start, design) {
  if (is.null(start)) {
    return(NULL)
  }
  if (!is.character(start) || length(start) != 1L || is.na(start) ||
      !start %in% c("gradual", "immediate")) {
    stop(
      "'start' must be \"gradual\", \"immediate\" or NULL",
      call. = FALSE
    )
  }
  if (identical(design, "fixed")) {
    stop(
      "'start' describes how a rotating design's cohorts are brought in; a fixed panel recruits one cohort, which is its whole design from wave 1",
      call. = FALSE
    )
  }
  start
}

#' @keywords internal
#' @noRd
.check_assurance <- function(assurance) {
  if (is.null(assurance)) {
    return(NULL)
  }
  if (!is.numeric(assurance) || length(assurance) != 1L || is.na(assurance) ||
      assurance <= 0 || assurance >= 1) {
    stop("'assurance' must be a probability in (0, 1)", call. = FALSE)
  }
  assurance
}

#' Refuse a panel the population cannot supply
#'
#' A fixed panel issues its whole cohort at once. A rotating one holds every
#' live cohort in sample at once. Either can outgrow a finite frame, and the
#' remedy differs from the one an oversized cross-section has, so the two
#' cases are named.
#' @keywords internal
#' @noRd
.check_panel_frame <- function(n_units, N, design, k) {
  if (!is.finite(N) || n_units <= N * (1 + 1e-9)) {
    return(invisible(TRUE))
  }
  msg <- if (identical(design, "fixed")) {
    sprintf(
      "the cohort would have to issue %s units, more than the population of %s; the target is out of reach at that wave under this retention",
      round(n_units, 1), N
    )
  } else {
    sprintf(
      "the %d cohorts alive at one occasion would hold %s units, more than the population of %s",
      k, round(n_units, 1), N
    )
  }
  stop(msg, call. = FALSE)
}

#' Units the whole-unit design puts in sample, and whether the frame holds them
#'
#' A fixed panel issues its one cohort. A rotating one holds every live cohort
#' at once, so rounding up costs a unit per cohort there.
#' @keywords internal
#' @noRd
.panel_units <- function(n_recruit, design, k) {
  if (identical(design, "fixed")) ceiling(n_recruit) else k * ceiling(n_recruit)
}

#' Report a whole-unit design a continuous one fits inside
#'
#' The continuous plan is checked against the frame before this, and clears
#' it. Rounding up to whole units can still not fit, by as much as one unit
#' per live cohort, and the number a planner would field is the rounded one.
#' Reported rather than refused: the plan the object stores is feasible, and
#' the remedy is either a design decision (cohorts of unequal size, which
#' this function does not describe) or a smaller target.
#' @keywords internal
#' @noRd
.warn_panel_whole_units <- function(n_recruit, N, design, k) {
  if (!is.finite(N)) {
    return(invisible(TRUE))
  }
  units <- .panel_units(n_recruit, design, k)
  if (units <= N) {
    return(invisible(TRUE))
  }
  warning(
    sprintf(
      "the whole-unit design puts %s units in sample against a population of %s, though the continuous design of %s fits; equal cohorts cannot be fielded at this size",
      .fmt_count_n(units), .fmt_count_n(N), format(round(n_recruit, 1))
    ),
    call. = FALSE
  )
  invisible(FALSE)
}

#' Report an assurance level a finite population cannot deliver
#'
#' A separate message from the one above, and a different situation: there the
#' continuous plan fits and only its rounding does not, while here no design
#' inside the frame reaches the target that often. Saying "though the
#' continuous assured design fits" of a figure that exceeds `N` would
#' contradict itself.
#' @keywords internal
#' @noRd
.warn_panel_assurance <- function(n_assured, N, design, k) {
  if (!is.finite(N)) {
    return(invisible(TRUE))
  }
  units <- .panel_units(n_assured, design, k)
  if (units <= N) {
    return(invisible(TRUE))
  }
  warning(
    sprintf(
      "the requested assurance would take %s units in sample, more than the population of %s; the level is unattainable, and no size inside the frame reaches the target that often",
      .fmt_count_n(units), .fmt_count_n(N)
    ),
    call. = FALSE
  )
  invisible(FALSE)
}

#' Name the waves that expect less than one respondent
#'
#' Their standard errors are computed on the same continuous scale as every
#' other wave and look no different, so a design whose later waves field
#' nobody reads as a design. The calculation is left alone, planning being
#' continuous throughout the package, and the waves are named.
#' @keywords internal
#' @noRd
.warn_panel_empty_waves <- function(n_wave) {
  thin <- which(n_wave < 1)
  if (length(thin) == 0L) {
    return(invisible(TRUE))
  }
  warning(
    sprintf(
      "wave%s %s expect%s less than one respondent (%s); the precision reported there describes a sample the design does not field",
      if (length(thin) > 1L) "s" else "", paste(thin, collapse = ", "),
      if (length(thin) > 1L) "" else "s",
      paste(sprintf("%.3g", n_wave[thin]), collapse = ", ")
    ),
    call. = FALSE
  )
  invisible(FALSE)
}

#' Assemble a panel from a recruitment size
#'
#' Shared by both directions: `n_panel()` solves the size from the target
#' and `prec_panel()` is handed one, and everything after that point is the
#' same arithmetic on the same rates.
#' @keywords internal
#' @noRd
.panel_result <- function(tgt, n_recruit, resp_rate, retention, q, design,
                          target_wave, assurance, start = NULL,
                          solved = NULL) {
  k <- length(q)
  if (!is.numeric(n_recruit) || length(n_recruit) != 1L ||
      !is.finite(n_recruit) || n_recruit <= 0) {
    stop("'n_recruit' must be a single positive size", call. = FALSE)
  }
  in_sample <- if (identical(design, "fixed")) n_recruit else k * n_recruit
  .check_panel_frame(in_sample, tgt$N, design, k)
  .warn_panel_target_rate(tgt, resp_rate)
  .warn_panel_whole_units(n_recruit, tgt$N, design, k)

  n_wave <- n_recruit * q
  .warn_panel_empty_waves(n_wave)
  prec <- lapply(n_wave, function(m) .panel_prec_at(tgt, m))
  loss <- c(1, q[-k]) - q
  total_loss <- 1 - q[[k]]

  waves <- data.frame(
    wave = seq_len(k),
    retention = c(NA_real_, as.numeric(retention)),
    q = q,
    loss_share = if (total_loss > 0) loss / total_loss else rep(NA_real_, k),
    n_resp = n_wave,
    se = vapply(prec, function(z) z$se, numeric(1L)),
    moe = vapply(prec, function(z) z$moe, numeric(1L)),
    cv = vapply(prec, function(z) z$cv, numeric(1L))
  )
  # Present for a proportion only, as `expected_cases` is on every other
  # result: it is where a min_cases floor can be read wave by wave, and the
  # waves past the target one are the place it stops holding.
  if (identical(tgt$type, "proportion")) {
    waves$expected_cases <- vapply(
      prec, function(z) z$expected_cases %||% NA_real_, numeric(1L)
    )
  }

  if (identical(design, "fixed")) {
    n_resp <- n_wave[[target_wave]]
    head <- prec[[target_wave]]
  } else {
    n_resp <- sum(n_wave)
    head <- .panel_prec_at(tgt, n_resp)
  }

  assured_feasible <- NULL
  n_assured <- if (!is.null(assurance)) {
    assured <- if (identical(design, "fixed")) {
      .assure_size(tgt$n_target, q[[target_wave]], assurance)
    } else {
      .panel_assure_rotating(tgt$n_target, q, assurance)
    }
    # Recorded, not warned: a warning does not survive into the object.
    assured_feasible <- .warn_panel_assurance(assured, tgt$N, design, k)
    assured
  }

  launch <- if (!is.null(start)) .panel_launch(tgt, n_recruit, q, start)

  params <- list(
    resp_rate = resp_rate,
    retention = as.numeric(retention),
    design = design,
    target_wave = target_wave,
    assurance = assurance
  )
  # Assigning NULL removes the element, so this one line is both branches.
  params$start <- start

  .new_svyplan_panel(
    design = design,
    n_recruit = n_recruit,
    n_target = tgt$n_target,
    n_resp = n_resp,
    n_assured = n_assured,
    assured_feasible = assured_feasible,
    target_wave = target_wave,
    waves = waves,
    prec = head,
    target = tgt$object,
    solved = solved,
    start = start,
    launch = launch$launch,
    launch_waves = launch$launch_waves,
    params = params
  )
}

#' Cohorts alive at each occasion of a launch, by wave
#'
#' One row per occasion, one column per wave, counting cohorts rather than
#' units, so the entrant count multiplies it and the continuous recruitment
#' is never rounded on the way in.
#'
#' A gradual launch recruits one cohort an occasion, so occasion \eqn{t}
#' holds waves 1 to \eqn{t}. An immediate launch splits its first occasion
#' into equal panels planned for life lengths from the full life down to one
#' occasion, all beginning at wave 1. Occasion \eqn{t} therefore holds the
#' \eqn{k - t + 1} launch panels still alive, all at wave \eqn{t}, plus one
#' intake cohort at each earlier wave. Both hold one
#' cohort per wave from occasion \eqn{k}, which is the steady state and is
#' where the two agree.
#' @keywords internal
#' @noRd
.panel_launch_cohorts <- function(k, start, t) {
  out <- numeric(k)
  if (t >= k) {
    out[] <- 1
    return(out)
  }
  if (identical(start, "gradual")) {
    out[seq_len(t)] <- 1
    return(out)
  }
  if (t > 1L) out[seq_len(t - 1L)] <- 1
  out[[t]] <- out[[t]] + (k - t + 1L)
  out
}

#' What a rotating design delivers at each occasion of its launch
#'
#' Reported through the occasion after the steady state is reached, so the
#' table shows it holding rather than leaving the reader to infer it. Counts
#' are continuous, from the continuous entrant size, whole units belonging on
#' the display path, which is where `print()` puts them.
#' @keywords internal
#' @noRd
.panel_launch <- function(tgt, n_recruit, q, start) {
  k <- length(q)
  periods <- seq_len(k + 1L)
  cohorts <- lapply(periods, function(t) .panel_launch_cohorts(k, start, t))

  n_in_sample <- vapply(cohorts, function(m) n_recruit * sum(m), numeric(1L))
  n_resp <- vapply(cohorts, function(m) n_recruit * sum(m * q), numeric(1L))
  prec <- lapply(n_resp, function(m) .panel_prec_at(tgt, m))
  # an immediate launch selects the whole design at its first occasion, and
  # one cohort an occasion after that
  entering <- rep(1, k + 1L)
  if (identical(start, "immediate")) entering[[1L]] <- k

  launch <- data.frame(
    period = periods,
    n_entrants = n_recruit * entering,
    n_in_sample = n_in_sample,
    n_resp = n_resp,
    se = vapply(prec, function(z) z$se, numeric(1L)),
    moe = vapply(prec, function(z) z$moe, numeric(1L)),
    cv = vapply(prec, function(z) z$cv, numeric(1L)),
    rmoe = vapply(prec, function(z) z$rmoe %||% NA_real_, numeric(1L)),
    steady_state = vapply(cohorts, function(m) all(m == 1), logical(1L))
  )
  if (identical(tgt$type, "proportion")) {
    launch$expected_cases <- vapply(
      prec, function(z) z$expected_cases %||% NA_real_, numeric(1L)
    )
  }

  # Long, not one column per wave: the waves present vary by occasion, so
  # columns would be mostly empty and their count would follow the life.
  wide <- do.call(rbind, lapply(periods, function(t) {
    m <- cohorts[[t]]
    w <- which(m > 0)
    data.frame(
      period = t,
      wave = w,
      n_issued = n_recruit * m[w],
      n_resp = n_recruit * m[w] * q[w]
    )
  }))
  rownames(wide) <- NULL

  list(launch = launch, launch_waves = wide)
}

#' Precision of the embedded estimand at a responding count
#'
#' The target's own precision statement is not an input to its evaluator.
#' `cv` and `rmoe` are formals of [prec_mean()] and [prec_prop()] with a
#' different meaning there, solving for the estimand rather than stating a
#' target, so passing a stored target through would silently reverse the
#' direction. Those, and the two arguments the panel sets itself, are the
#' only ones dropped: everything else the target validated is carried.
#' @keywords internal
#' @noRd
.panel_prec_at <- function(tgt, n_resp) {
  fn <- switch(tgt$type, mean = prec_mean.default, prec_prop.default)
  keep <- setdiff(
    intersect(names(tgt$params), names(formals(fn))),
    c("n", "resp_rate", "moe", "cv", "rmoe")
  )
  args <- tgt$params[keep]
  args$n <- n_resp
  args$resp_rate <- 1
  if (identical(tgt$type, "proportion")) {
    args$method <- tgt$method %||% "wald"
  }
  do.call(fn, args)
}

#' Entrants per occasion whose pooled respondents clear a target
#'
#' At steady state the respondents at one occasion are
#' \eqn{\sum_s \mathrm{Binomial}(e, q_s)}{sum_s Binomial(e, q_s)} over the live cohorts, whose
#' cumulative response probabilities differ, so this is Poisson-binomial and
#' not \eqn{\mathrm{Binomial}(ke, \bar{q})}{Binomial(ke, qbar)}. The distribution is assembled
#' exactly, one convolution per cohort, and its tail is monotone in \eqn{e},
#' so the smallest sufficient \eqn{e} is bracketed by doubling and then
#' bisected.
#' @keywords internal
#' @noRd
.panel_assure_rotating <- function(need, q, level) {
  need <- ceiling(need)
  if (need <= 0) {
    return(0)
  }
  clears <- function(e) {
    pmf <- 1
    for (qi in q) {
      pmf <- .conv_pmf(pmf, stats::dbinom(0:e, e, qi))
    }
    if (length(pmf) <= need) {
      return(FALSE)
    }
    sum(pmf[(need + 1L):length(pmf)]) >= level
  }
  lo <- 0L
  hi <- 1L
  while (!clears(hi)) {
    lo <- hi
    hi <- hi * 2L
    if (hi > 1e8) {
      stop("the rotating assurance search did not converge", call. = FALSE)
    }
  }
  while (hi - lo > 1L) {
    mid <- lo + (hi - lo) %/% 2L
    if (clears(mid)) hi <- mid else lo <- mid
  }
  as.double(hi)
}

#' Distribution of a sum of two independent counts
#' @keywords internal
#' @noRd
.conv_pmf <- function(a, b) {
  if (length(a) == 1L) {
    return(a * b)
  }
  pmax(stats::convolve(a, rev(b), type = "open"), 0)
}
