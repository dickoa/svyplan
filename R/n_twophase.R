#' Two-phase sample allocation
#'
#' Allocate a two-phase (double) sample. Phase 1 draws a large sample and
#' measures a variable the frame does not carry. Phase 2 subsamples those
#' units and measures the variable of interest on them. Solves for the
#' phase sizes and the per-stratum subsampling fractions, either to reach
#' a target coefficient of variation or to minimize it under a fixed
#' budget.
#'
#' What phase 1 measures decides which design this is. A stratification
#' variable observed on the phase-1 sample lets phase 2 be stratified on
#' something the frame could not supply. Response status observed in phase
#' 1 makes the follow-up of nonrespondents a second phase. One allocator
#' covers both. Double sampling for stratification subsamples every
#' phase-2 stratum, nonresponse follow-up carries the phase-1 respondents
#' through untouched and subsamples only the nonrespondents, and the two
#' differ only in which strata are marked `take_all`.
#'
#' @param frame For the default method: a data frame with one row per
#'   phase-2 stratum, in the [n_alloc()] column vocabulary (`N`, `sd`,
#'   `mean`, `unit_cost`, `take_all`) but on a narrower contract. `sd` is
#'   required where [n_alloc()] also accepts `var`, and the cost and
#'   `p`/`mean` rules differ, so a table built for one may need adjusting
#'   before it is passed to the other. Columns:
#'   \describe{
#'     \item{`N`}{Stratum size (**required**). Only relative size matters,
#'       since the stratum weights are `N / sum(N)`.}
#'     \item{`sd`}{Within-stratum standard deviation of the variable of
#'       interest (**required**).}
#'     \item{`mean`}{Stratum mean. Supply it to let the between-stratum
#'       component be derived, which is what makes stratifying worthwhile.
#'       Omit it, or pass `between = 0`, when the strata are not expected
#'       to differ in level, as in nonresponse follow-up.}
#'     \item{`unit_cost`}{Incremental phase-2 cost per unit. Default 1.
#'       Zero means the stratum is already measured, which with
#'       `take_all = TRUE` is exactly a phase-1 respondent group.}
#'     \item{`deff`}{Design effect of the phase-2 subsample in this
#'       stratum, default 1. This one applies to the *residual* variation
#'       that phase 2 has to measure, under the phase-2 field design. See
#'       Details for why it is not the same design effect as
#'       `phase1_deff`.}
#'     \item{`resp_rate`}{Expected phase-2 completion rate in this
#'       stratum, in (0, 1], default 1. It divides the stratum's
#'       contribution to the variance, so a stratum that responds poorly
#'       is treated as carrying less information per issued unit.}
#'     \item{`take_all`}{Logical. `TRUE` carries every unit phase 1
#'       successfully classified into phase 2. On the issued-unit scale
#'       `nu` is measured on, that is `nu = resp_rate`, the phase-1
#'       classification rate, and not 1. Among classified units the
#'       fraction is one. The two coincide only when phase 1 classifies
#'       everyone. Default `FALSE`.}
#'     \item{`stratum`}{Optional label.}
#'   }
#'   For `svyplan_prec` objects: a precision result from [prec_twophase()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param phase1_cost Cost per phase-1 unit (> 0). This is paid on every
#'   unit screened, whether or not it reaches phase 2.
#' @param cv Target coefficient of variation. Specify exactly one of `cv`
#'   or `budget`.
#' @param budget Total budget. Specify exactly one of `cv` or `budget`.
#' @param between Between-stratum variance component. `NULL` (default)
#'   derives it from the `mean` column as `sum(W * (mean - ybar)^2)`, or
#'   uses 0 when no `mean` is given. Supply it directly when you know the
#'   stratification gain but not the stratum means.
#' @param mu Population mean, used to turn the variance into a coefficient
#'   of variation. `NULL` (default) derives it from `mean` when present.
#'   Required in `cv` mode when no `mean` column is supplied.
#' @param N Population size for the phase-1 finite population correction.
#'   `Inf` (default) applies none.
#' @param phase1_deff Design effect of the phase-1 sample, default 1. It
#'   multiplies the between-stratum component, which is the part of the
#'   variance phase 1 is responsible for. Set it above 1 when phase 1 is
#'   clustered, as an area screener is.
#' @param resp_rate Expected phase-1 response, screening or successful
#'   classification rate, in (0, 1], default 1. It divides the
#'   between-stratum component and, because phase 2 can only draw from
#'   the units phase 1 actually classified, it also caps every
#'   subsampling fraction at `resp_rate`. In a nonresponse follow-up
#'   frame leave it at 1, because there the strata *are* response status, so
#'   classification succeeds for every unit and setting it again would
#'   count the same loss twice.
#' @param n_phase1 Fixed phase-1 sample size, or `NULL` (default) to solve
#'   for it. Supply it when phase 1 has already run, when it is an existing
#'   survey or panel, or when its size is set by field capacity rather than
#'   by this design. The relative allocation across strata is
#'   `S_h sqrt(d_2h / c_h)` in every mode and does not move. Only the
#'   overall scale does, pinned by the budget left after phase 1 or by what
#'   reaching `cv` requires at that size, so fixing `n_phase1` can only
#'   match or lose to leaving it free. It becomes infeasible in two ways.
#'   The budget may not reach phase 2 once phase 1 and any `take_all`
#'   strata are paid for. Or the target `cv` may sit below the floor left
#'   when phase 2 carries every classified unit through, the
#'   between-stratum component `d_1 A / r_1` plus the phase-2 residual at
#'   `nu_h = resp_rate`, which no amount of subsampling can beat.
#' @param assurance Probability in (0, 1), or `NULL` (default). Planning at
#'   the expected respondent count leaves roughly half of all designs
#'   short. Supplying a level reports, alongside the expected design, the
#'   issued sizes for which the required respondents arrive with at least
#'   that probability, from the binomial distribution of respondents.
#'
#'   The level is **marginal**, holding stratum by stratum. Every stratum
#'   clearing its target at once has the product probability and so is
#'   lower, materially so with many strata, since three strata at 0.95 give
#'   about 0.86 together. Raise the level for a familywise guarantee. It is
#'   also conditional on the phase-1 pool. The assured phase-2 issue is
#'   compared against what phase 1 supplies, and a stratum needing more
#'   than its pool is reported in a warning, since the answer there is to
#'   enlarge phase 1 rather than to over-issue. Phase-1 composition is
#'   itself random unless phase 1 is a census, so the expected stratum
#'   shares behind that pool are not simultaneous lower bounds either.
#' @param single_deff Design effect of the single-phase comparator,
#'   default 1. It is a separate number from the stratum `deff` column
#'   because the comparator need not be fielded the same way as phase 2.
#' @param single_resp_rate Expected response rate of the single-phase
#'   comparator, in (0, 1], default 1.
#' @param single_cost Cost per unit of the single-phase design that skips
#'   phase 1 and measures the variable of interest directly. `NULL` (default)
#'   uses `sum(share * unit_cost)`, which is the right baseline when
#'   `unit_cost` is what measuring one unit costs. In a nonresponse
#'   follow-up design it is not, since a `unit_cost` of 0 there means "already
#'   measured", not "free", and the right baseline is the cost of one
#'   completed interview without any follow-up, `phase1_cost / resp_rate`.
#'   Supply it in that case.
#' @param fixed_cost Fixed overhead, not depending on sample size.
#'   Default 0. In budget mode only `budget - fixed_cost` is allocatable.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_twophase` object with components:
#' \describe{
#'   \item{`n`}{Named numeric vector `c(n_phase1 = , n_phase2 = )`, both
#'     continuous, counting units **issued**. `n_phase2` is the expected
#'     issue, `sum(W_h * nu_h) * n_phase1`.}
#'   \item{`responding`}{The same two quantities in expected
#'     **respondents**, `r_1 n_a` and `sum(r_{2h} W_h nu_h) n_a`. Issued
#'     and responding are different quantities and neither is an
#'     effective sample size. There is no single effective size for a
#'     two-phase design, because the two variance components carry
#'     different design effects.}
#'   \item{`cv`}{Coefficient of variation achieved.}
#'   \item{`cost`}{Total cost, including `fixed_cost`.}
#'   \item{`detail`}{Per-stratum table: `stratum`, `N`, `share`, `sd`,
#'     `unit_cost`, `deff`, `resp_rate`, `nu` (the subsampling fraction,
#'     as a share of the phase-1 units *issued*, and with `resp_rate` below 1
#'     it is capped there, since phase 2 can only draw from the units
#'     phase 1 classified, and the fraction taken among those classified
#'     units is `nu / resp_rate`),
#'     `n_issued`, `n_resp` (its expected responding part), and
#'     `take_all`, which is `TRUE`
#'     for strata that were supplied pinned **and** for those the
#'     allocator truncated at the classification rate.}
#'   \item{`single_phase`}{The design that skips phase 1 and measures the
#'     variable of interest directly, as `c(n = , cv = , cost = )`, plus
#'     `reaches_target`, and `better`, `TRUE` when that plainer design
#'     wins. Two-phase sampling is not always an improvement and this
#'     comparison is the check that says so. In `cv` mode a comparator
#'     that cannot reach the target from a frame of size `N` is capped
#'     there, reports the coefficient of variation a census of it would
#'     achieve, and never wins on cost alone.}
#'   \item{`operational`}{The whole-unit field design: `n` (both phases),
#'     the per-stratum `n_int`, and the `cost` and `cv` it actually
#'     achieves. Budget mode floors and then buys back whole units in order
#'     of variance reduction per unit cost, so the integer design stays
#'     inside the budget, while cv mode rounds up. Both are bounded by the
#'     phase-2 pool, so a stratum whose expected phase-1 yield rounds below
#'     one unit is left unsampled with a warning, and the operational `cv`
#'     is then `Inf`. With `assurance` set it also carries `assured`,
#'     `assured_phase1`, `pool` and `assured_cost`.}
#'   \item{`params`}{Validated inputs, read back by [prec_twophase()]. They
#'     include the design decisions a stratum table cannot express, so that
#'     `n_twophase(prec_twophase(fit))` rebuilds this problem rather than a
#'     new unconstrained one, namely the `take_all` set as *supplied*, a fixed
#'     `n_phase1`, `assurance`, and the comparator settings. There is no
#'     `predict()` method for this class.}
#' }
#'
#' @details
#' Write \eqn{W_h} for the stratum weight, \eqn{S_h} for the within-stratum
#' standard deviation, \eqn{c_h} for the incremental phase-2 unit cost,
#' \eqn{c_a} for the phase-1 unit cost, and \eqn{\nu_h} for the fraction of
#' the phase-1 stratum carried into phase 2. Ignoring the finite population
#' correction,
#'
#' \deqn{V = \frac{1}{n_a}\left[A + \sum_h \frac{W_h S_h^2}{\nu_h}\right],
#'       \qquad C = n_a\left[c_a + \sum_h c_h W_h \nu_h\right],}{V = 1/n_a [A + sum_h (W_h S_h^2)/nu_h ], C = n_a [c_a + sum_h c_h W_h nu_h ],}
#'
#' where \eqn{A} is the between-stratum component. Minimizing \eqn{VC},
#' which is free of \eqn{n_a}, gives
#'
#' \deqn{\nu_h = \frac{S_h}{\sqrt{c_h}}
#'       \sqrt{\frac{\tilde c_a}{\tilde A}}, \qquad
#'       \tilde A = A + \sum_{h \in P} W_h S_h^2, \qquad
#'       \tilde c_a = c_a + \sum_{h \in P} c_h W_h,}{nu_h = S_h/sqrt(c_h) sqrt(ctilde_a/Atilde), Atilde = A + sum_(h in P) W_h S_h^2, ctilde_a = c_a + sum_(h in P) c_h W_h,}
#'
#' where \eqn{P} is the set of pinned strata. The shape is Neyman-like,
#' \eqn{\nu_h \propto S_h/\sqrt{c_h}}{nu_h proportional to S_h/sqrt(c_h)}, but the overall scale is set by the
#' *between*-stratum variance, so weak stratification pushes every
#' \eqn{\nu_h} up, toward keeping everything phase 1 found.
#'
#' Any \eqn{\nu_h} above the cap is truncated there and the stratum joins
#' \eqn{P}, which changes \eqn{\tilde A}{Atilde} and \eqn{\tilde c_a}{ctilde_a} and so the
#' remaining strata are re-solved. The cap is `resp_rate` rather than 1.
#' Phase 2 can only subsample the units phase 1 succeeded in classifying, so a
#' stratum is "take-all" once it keeps all of those, not all of \eqn{N_h}.
#' Strata are pinned one at a time in decreasing \eqn{S_h/\sqrt{c_h}}{S_h/sqrt(c_h)},
#' because pinning lowers the multiplier and can bring others back below
#' the cap.
#'
#' ## The finite population correction
#'
#' With `N` set, the correction applies component by component and not to
#' the variance as a whole:
#'
#' \deqn{V = A\left(\frac{1}{n_a}-\frac{1}{N}\right)
#'       + \sum_h W_h S_h^2
#'         \left\{\frac{1}{\nu_h n_a}-\frac{1}{N}\right\}.}{V = A (1/n_a-1/N ) + sum_h W_h S_h^2 \{1/(nu_h n_a)-1/N \}.}
#'
#' Phase 1 estimates the stratum weights, so a phase-1 census makes the
#' between-stratum term vanish. The phase-2 terms vanish only when phase 2
#' also measures every unit it found. Scaling the whole variance by
#' \eqn{1 - n_a/N} would subtract \eqn{1/(\nu_h N)} in place of
#' \eqn{1/N} and so report a phase-1 census as a zero-variance design.
#' Note that with `resp_rate` below 1 a phase-1 census is not a complete
#' classification of the frame, and the residual it leaves is real rather
#' than an artifact.
#'
#' ## Two design effects, not one
#'
#' The variance has two components and each carries its own design effect,
#'
#' \deqn{V \approx \frac{1}{n_a}\left[d_1 A
#'       + \sum_h \frac{d_{2h} W_h S_h^2}{\nu_h}\right],}{V approx 1/n_a [d_1 A + sum_h (d_2h W_h S_h^2)/nu_h ],}
#'
#' so the optimum becomes \eqn{\nu_h \propto S_h\sqrt{d_{2h}/c_h}}{nu_h proportional to S_h sqrt(d_2h/c_h)} and
#' \eqn{\tilde A}{Atilde} is built from \eqn{d_1 A} and the pinned strata's
#' \eqn{d_{2h} W_h S_h^2}{d_2h W_h S_h^2}. Inflating the combined variance by a single
#' design effect instead is a different and wrong model.
#'
#' The two are design effects **for different variables**. `phase1_deff`
#' applies to what phase 1 explains, the between-stratum contrast.
#' The `deff` column applies to what is left for phase 2 to measure, the
#' within-stratum residual. A clustered phase 1 raises the first. A
#' stratifier that absorbs geographic variation lowers the second.
#'
#' They are not independent, and the trap runs in the direction that looks
#' attractive. A stratifier built purely from between-cluster structure
#' shrinks the residual and so lowers `deff`, but it drives the fitted
#' values toward constancy within clusters and so raises `phase1_deff`
#' toward the cluster size. Under a clustered phase 1 the first term can
#' then dominate and the design loses to a plain single-phase sample even
#' though the residual design effect looks favorable. When the
#' stratification explains little, the phase-2 residual is close to the
#' whole variance and the stratum `deff` should be close to `single_deff`.
#' Setting it well below with a weak stratifier assumes away the cost of
#' the design.
#'
#' ## Response rates
#'
#' All costs and sample sizes count units **issued**, so a cost quoted
#' per completed interview has to be converted before it is passed in:
#' with a contact cost and a completion cost, the issued-basis figure is
#' `c_contact + resp_rate * c_complete`.
#'
#' Response enters each component separately, `d_1` becoming `d_1/r_1`
#' and `d_{2h}` becoming `d_{2h}/r_{2h}`, so a phase-2 stratum that
#' responds poorly is subsampled differently rather than the whole design
#' being inflated by one factor. A response rate common to everything
#' scales the sample needed for a precision target but cancels from the
#' relative allocation. Phase- or stratum-specific rates change it.
#'
#' This is an expected-information calculation and **not** a bias
#' adjustment. Dividing by a response rate assumes response is ignorable
#' within the strata you supplied, which is a substantive assumption
#' about the strata, not a property of the arithmetic. No sample size
#' removes nonresponse bias. Where that assumption is uncomfortable,
#' modeling the nonresponse as a follow-up phase, as below, is the
#' alternative the design itself offers.
#'
#' Nonresponse follow-up is the two-stratum case. Respondents take
#' `unit_cost = 0` and `take_all = TRUE`, the nonrespondent stratum is left
#' free to be subsampled, and there is no between-stratum component. The
#' optimum reduces to
#' \eqn{\nu = \sqrt{c_1/(c_2\theta)}}{nu = sqrt(c_1/(c_2 theta))} for a phase-1 response rate
#' \eqn{\theta}, which is the standard result.
#'
#' @references
#' Fuller, W. A. (2009). \emph{Sampling Statistics}, Sect. 3.3. Wiley.
#'
#' Neyman, J. (1938). Contribution to the theory of sampling human
#' populations. \emph{Journal of the American Statistical Association},
#' 33(201), 101--116.
#'
#' Saerndal, C.-E., Swensson, B., and Wretman, J. (1992). \emph{Model
#' Assisted Survey Sampling}, Sect. 15.4. Springer.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). \emph{Practical
#' Tools for Designing and Weighting Survey Samples}, 2nd edition,
#' Sect. 17.5.2. Springer.
#'
#' @family two-phase design functions
#' @seealso [prec_twophase()] for the inverse, [n_alloc()] for
#'   single-phase stratified allocation, [n_cluster()] for the multistage
#'   analogue.
#'
#' @examples
#' # Double sampling for stratification. A screener splits the frame into
#' # four groups, and the variable of interest is measured on a subsample.
#' frame <- data.frame(
#'   stratum   = c("A", "B", "C", "D"),
#'   N         = c(3500, 2500, 2500, 1500),
#'   sd        = c(12, 25, 8, 40),
#'   mean      = c(40, 70, 35, 90),
#'   unit_cost = c(2, 5, 1, 9)
#' )
#' n_twophase(frame, phase1_cost = 1, budget = 50000)
#'
#' # Weak stratification pushes the fractions up, and strata that reach 1
#' # are reported as take-all.
#' flat <- transform(frame, mean = c(55, 56, 55, 56))
#' n_twophase(flat, phase1_cost = 1, budget = 50000)
#'
#' # Nonresponse follow-up: respondents are already measured, so they cost
#' # nothing more and are all kept. Only the nonrespondents are subsampled.
#' theta <- 0.5
#' nrfu <- data.frame(
#'   stratum   = c("respondents", "nonrespondents"),
#'   N         = c(theta, 1 - theta),
#'   sd        = c(1, 1),
#'   unit_cost = c(0, 200),
#'   take_all  = c(TRUE, FALSE)
#' )
#' # a unit_cost of 0 means "already measured", not "free", so the
#' # single-phase baseline has to be named: one completed interview
#' n_twophase(nrfu, phase1_cost = 50, budget = 100000, mu = 1,
#'            single_cost = 50 / theta)
#'
#' @export
n_twophase <- function(frame, ...) {
  if (!missing(frame)) {
    .res <- .dispatch_plan(frame, "frame", n_twophase.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("n_twophase")
}

#' @rdname n_twophase
#' @export
n_twophase.default <- function(
  frame,
  ...,
  phase1_cost,
  cv = NULL,
  budget = NULL,
  between = NULL,
  mu = NULL,
  N = Inf,
  n_phase1 = NULL,
  assurance = NULL,
  phase1_deff = 1,
  resp_rate = 1,
  single_deff = 1,
  single_resp_rate = 1,
  single_cost = NULL,
  fixed_cost = 0,
  plan = NULL
) {
  merged <- .merge_plan_args(plan, n_twophase.default, match.call(), environment())
  if (!is.null(merged)) {
    return(do.call(n_twophase.default, c(merged, list(...))))
  }
  .check_unused_dots(...)

  spec <- .twophase_frame(frame, between, mu)
  # an evaluation can price a fully pinned design; an allocation has nothing
  # left to choose
  if (all(spec$pinned)) {
    stop("every stratum is 'take_all'; nothing is left to subsample",
         call. = FALSE)
  }
  if (missing(phase1_cost)) {
    stop("'phase1_cost' is required", call. = FALSE)
  }
  check_scalar(phase1_cost, "phase1_cost")
  .check_twophase_deff(phase1_deff, "phase1_deff")
  .check_twophase_deff(single_deff, "single_deff")
  check_resp_rate(resp_rate)
  check_resp_rate(single_resp_rate)
  check_population_size(N)
  check_fixed_cost(fixed_cost, budget)
  if (is.null(cv) == is.null(budget)) {
    stop("supply exactly one of 'cv' or 'budget'", call. = FALSE)
  }
  if (!is.null(n_phase1)) {
    check_scalar(n_phase1, "n_phase1")
    if (n_phase1 > N) {
      stop("'n_phase1' exceeds the population size 'N'", call. = FALSE)
    }
  }
  if (!is.null(assurance)) {
    if (!is.numeric(assurance) || length(assurance) != 1L ||
        is.na(assurance) || assurance <= 0 || assurance >= 1) {
      stop("'assurance' must be a probability in (0, 1)", call. = FALSE)
    }
  }

  d2_eff <- spec$d2 / spec$r2
  a_eff <- phase1_deff * spec$A / resp_rate
  if (!is.null(budget)) check_scalar(budget, "budget")
  if (!is.null(cv)) {
    check_scalar(cv, "cv")
    if (!is.finite(spec$mu) || spec$mu == 0) {
      stop("'mu' is required in cv mode; supply it or a 'mean' column",
           call. = FALSE)
    }
  }

  # Deff-inflated, response-unadjusted: a census removes all the variance,
  # failing to reach some units does not.
  var_pop <- phase1_deff * spec$A + sum(spec$d2 * spec$W * spec$S^2)
  # at nu_max rather than 1, the floor no amount of subsampling can beat
  var_floor <- a_eff + sum(d2_eff * spec$W * spec$S^2 / resp_rate)

  lambda_fn <- NULL
  if (!is.null(n_phase1)) {
    if (!is.null(budget)) {
      allocatable <- (budget - fixed_cost) / n_phase1
      # take-all strata must be measured whatever is left, so their cost is
      # as mandatory as the screening itself
      mandatory <- phase1_cost +
        sum(spec$cost[spec$pinned] * spec$W[spec$pinned] * resp_rate)
      if (allocatable <= mandatory) {
        stop(sprintf(
          paste("the budget does not reach phase 2 at this 'n_phase1':",
                "screening %s units%s costs %s of the %s available"),
          format(n_phase1),
          if (any(spec$pinned)) " and measuring the take-all strata" else " alone",
          format(round(mandatory * n_phase1)),
          format(round(budget - fixed_cost))), call. = FALSE)
      }
      lambda_fn <- .twophase_fixed_lambda(spec$W, spec$S, spec$cost, d2_eff,
                                          allocatable, "budget")
    } else {
      target_var <- n_phase1 * ((cv * spec$mu)^2 + var_pop / N)
      if (target_var <= var_floor) {
        stop(sprintf(
          paste("the target 'cv' is out of reach at this 'n_phase1':",
                "carrying every classified unit into phase 2 still leaves",
                "cv = %s, because the between-stratum component and the",
                "phase-2 residual are both already paid for"),
          formatC(.twophase_cv(n_phase1, var_floor, spec$mu, N, var_pop),
                  format = "f", digits = 5)), call. = FALSE)
      }
      lambda_fn <- .twophase_fixed_lambda(spec$W, spec$S, spec$cost, d2_eff,
                                          target_var, "cv")
    }
  }

  fit <- .twophase_nu(spec$W, spec$S, spec$cost, d2_eff, phase1_cost,
                      a_eff, spec$pinned, nu_max = resp_rate,
                      lambda_fn = lambda_fn)
  nu <- fit$nu
  # a stratum with no variance to measure is legitimately left at zero; any
  # other non-positive fraction means the constraints could not be met
  if (any(!is.finite(nu)) || any(nu < 0) || any(nu[spec$S > 0] <= 0)) {
    stop(paste("the constraints imply a non-positive subsampling fraction, so",
               "there is no feasible allocation; check 'budget' or 'n_phase1'",
               "against the cost of the take-all strata"), call. = FALSE)
  }
  var_unit <- a_eff + sum(ifelse(spec$S > 0, d2_eff * spec$W * spec$S^2 / nu, 0))
  cost_unit <- phase1_cost + sum(spec$cost * spec$W * nu)

  if (!is.null(n_phase1)) {
    n1 <- n_phase1
    cost <- fixed_cost + n1 * cost_unit
  } else if (!is.null(budget)) {
    n1 <- (budget - fixed_cost) / cost_unit
    if (!is.finite(n1) || n1 <= 0) {
      stop("'budget' does not cover a single phase-1 unit", call. = FALSE)
    }
    n1 <- min(n1, N)
    cost <- fixed_cost + n1 * cost_unit
  } else {
    n1 <- var_unit / ((cv * spec$mu)^2 + var_pop / N)
    if (n1 > N) {
      stop("the target 'cv' is not reachable from a frame of size 'N'",
           call. = FALSE)
    }
    cost <- fixed_cost + n1 * cost_unit
  }

  cv_out <- .twophase_cv(n1, var_unit, spec$mu, N, var_pop)
  if (!is.null(cv) && is.finite(cv_out) && cv_out > cv * (1 + 1e-6)) {
    stop(sprintf(
      "the allocation reaches cv = %s, short of the target %s",
      formatC(cv_out, format = "f", digits = 5), format(cv)), call. = FALSE)
  }
  single <- .twophase_single(spec, cv, budget, N, fixed_cost, single_cost,
                             single_deff, single_resp_rate)
  # Budget mode compares precision at equal spend, cv mode compares spend at
  # equal precision. In cv mode a cheaper comparator that never reaches the
  # target is not an alternative, so it does not win on cost alone.
  single$better <- if (!is.null(budget)) {
    isTRUE(is.finite(single$cv) && is.finite(cv_out) && single$cv < cv_out)
  } else {
    isTRUE(single$reaches_target && is.finite(single$cost) && single$cost < cost)
  }

  mode <- if (is.null(budget)) "cv" else "budget"
  operational <- .twophase_integerize(
    n1, nu, spec, phase1_cost, a_eff, d2_eff, resp_rate, mode,
    budget %||% Inf, fixed_cost, N, spec$mu, var_pop
  )
  if (!is.null(assurance)) {
    # the design plans on r * issued respondents; assurance asks what to
    # issue so that many actually arrive with probability `assurance`
    operational$assured <- .assure_size(
      operational$n_int * spec$r2, spec$r2, assurance)
    operational$assured_phase1 <- .assure_size(
      operational$n[["n_phase1"]] * resp_rate, resp_rate, assurance)
    operational$pool <- floor(resp_rate * spec$W * operational$assured_phase1)
    operational$assured_cost <- fixed_cost +
      phase1_cost * operational$assured_phase1 +
      sum(spec$cost * operational$assured)
    short <- operational$assured > operational$pool
    if (any(short)) {
      warning(sprintf(
        paste("the assured phase-2 issue exceeds the pool phase 1 supplies in",
              "%s (needs %s from %s); the assurance level does not hold there",
              "unless phase 1 is enlarged"),
        paste(sprintf("'%s'", spec$stratum[short]), collapse = ", "),
        paste(format(operational$assured[short]), collapse = ", "),
        paste(format(operational$pool[short]), collapse = ", ")),
        call. = FALSE)
    }
  }

  detail <- data.frame(
    stratum = spec$stratum,
    N = spec$N,
    share = spec$W,
    sd = spec$S,
    unit_cost = spec$cost,
    deff = spec$d2,
    resp_rate = spec$r2,
    nu = nu,
    n_issued = nu * spec$W * n1,
    n_int = operational$n_int,
    n_resp = spec$r2 * nu * spec$W * n1,
    take_all = fit$pinned,
    stringsAsFactors = FALSE
  )

  .new_svyplan_twophase(
    n = c(n_phase1 = n1, n_phase2 = sum(spec$W * nu) * n1),
    responding = c(n_phase1 = resp_rate * n1,
                   n_phase2 = sum(spec$r2 * spec$W * nu) * n1),
    cv = cv_out,
    cost = cost,
    detail = detail,
    single_phase = single,
    operational = operational,
    params = list(
      phase1_cost = phase1_cost, between = spec$A, mu = spec$mu, N = N,
      phase1_deff = phase1_deff, single_deff = single_deff,
      resp_rate = resp_rate, single_resp_rate = single_resp_rate,
      n_phase1_fixed = n_phase1, assurance = assurance,
      fixed_cost = fixed_cost, mode = mode,
      cv_target = cv, budget = budget, single_cost = single_cost,
      # the strata the caller pinned, which is not the same set as the
      # `take_all` column: that also flags the ones the allocator truncated
      take_all = spec$pinned
    )
  )
}

#' Validate and normalize the phase-2 stratum frame
#' @keywords internal
#' @noRd
.twophase_frame <- function(frame, between, mu) {
  if (!is.data.frame(frame) || nrow(frame) == 0L) {
    stop("'frame' must be a non-empty data frame", call. = FALSE)
  }
  if ("cost" %in% names(frame)) {
    stop(
      "the per-stratum cost column is 'unit_cost'; 'cost' is the total field cost",
      call. = FALSE
    )
  }
  for (needed in c("N", "sd")) {
    if (!needed %in% names(frame)) {
      stop(sprintf("'frame' must contain a '%s' column", needed), call. = FALSE)
    }
  }
  N_h <- frame$N
  if (!is.numeric(N_h) || anyNA(N_h) || any(!is.finite(N_h)) || any(N_h <= 0)) {
    stop("'N' must contain positive finite values", call. = FALSE)
  }
  S <- frame$sd
  if (!is.numeric(S) || anyNA(S) || any(!is.finite(S)) || any(S < 0)) {
    stop("'sd' must contain non-negative finite values", call. = FALSE)
  }
  cost <- if ("unit_cost" %in% names(frame)) frame$unit_cost else rep(1, nrow(frame))
  if (!is.numeric(cost) || anyNA(cost) || any(!is.finite(cost)) || any(cost < 0)) {
    stop("'unit_cost' must contain non-negative finite values", call. = FALSE)
  }
  d2 <- if ("deff" %in% names(frame)) frame$deff else rep(1, nrow(frame))
  if (!is.numeric(d2) || anyNA(d2) || any(!is.finite(d2)) || any(d2 <= 0)) {
    stop("the 'deff' column must contain positive finite values", call. = FALSE)
  }
  r2 <- if ("resp_rate" %in% names(frame)) frame$resp_rate else rep(1, nrow(frame))
  if (!is.numeric(r2) || anyNA(r2) || any(!is.finite(r2)) ||
      any(r2 <= 0) || any(r2 > 1)) {
    stop("the 'resp_rate' column must contain values in (0, 1]", call. = FALSE)
  }
  pinned <- .check_take_all(frame[["take_all"]], nrow(frame))
  W <- N_h / sum(N_h)

  if (!is.null(between)) {
    check_scalar(between, "between", positive = FALSE)
    if (between < 0) stop("'between' must be non-negative", call. = FALSE)
    A <- between
  } else if ("mean" %in% names(frame)) {
    m <- frame$mean
    if (!is.numeric(m) || anyNA(m) || any(!is.finite(m))) {
      stop("'mean' must contain finite values", call. = FALSE)
    }
    A <- sum(W * (m - sum(W * m))^2)
  } else {
    A <- 0
  }
  if (is.null(mu) && "mean" %in% names(frame)) {
    # the cv paths test this against zero, so a set of stratum means that
    # cancels has to reach them as an exact zero rather than as whatever
    # the summation left behind
    mu <- .aggregate_mean(W, frame$mean)
  }
  list(
    stratum = .twophase_labels(frame),
    N = N_h, W = W, S = S, cost = cost, d2 = d2, r2 = r2, pinned = pinned, A = A,
    mu = if (is.null(mu)) NA_real_ else mu
  )
}

#' Stratum labels for a two-phase frame
#'
#' Labels identify strata in warnings, assurance diagnostics and the round
#' trip, so a duplicated or missing one silently attaches a message to the
#' wrong row. Generated row numbers are used when the column is absent, and
#' those are unambiguous by construction.
#' @keywords internal
#' @noRd
.twophase_labels <- function(frame) {
  if (!"stratum" %in% names(frame)) {
    return(as.character(seq_len(nrow(frame))))
  }
  labels <- as.character(frame$stratum)
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
    stop("'stratum' must not contain missing or empty labels", call. = FALSE)
  }
  if (anyDuplicated(labels) > 0L) {
    stop(
      sprintf("'stratum' must be unique; duplicated: %s",
              paste(sQuote(unique(labels[duplicated(labels)])), collapse = ", ")),
      call. = FALSE
    )
  }
  labels
}

#' Optimal subsampling fractions with truncation at 1
#'
#' Pinning a stratum raises both the effective between component and the
#' effective phase-1 cost, which lowers the multiplier and can bring other
#' strata back under 1, so strata are pinned one at a time. The ordering is
#' taken from `S / sqrt(cost)` rather than from `nu`: the two agree
#' whenever the multiplier is finite, and only the former is defined when
#' it is not.
#' @keywords internal
#' @noRd
.twophase_nu <- function(W, S, cost, d2, phase1_cost, A, pinned, nu_max = 1,
                         lambda_fn = NULL) {
  key <- ifelse(cost > 0, S * sqrt(d2 / cost), Inf)
  if (is.null(lambda_fn)) {
    lambda_fn <- function(a_eff, c_eff, free) {
      if (a_eff > 0) sqrt(c_eff / a_eff) else Inf
    }
  }
  repeat {
    a_eff <- A + sum(d2[pinned] * W[pinned] * S[pinned]^2 / nu_max)
    c_eff <- phase1_cost + sum(cost[pinned] * W[pinned] * nu_max)
    free <- !pinned
    lambda <- lambda_fn(a_eff, c_eff, free)
    nu <- rep(nu_max, length(W))
    nu[free] <- if (is.finite(lambda)) key[free] * lambda else Inf
    over <- free & nu > nu_max
    if (!any(over)) {
      return(list(nu = nu, pinned = pinned))
    }
    pinned[which.max(ifelse(over, key, -Inf))] <- TRUE
  }
}

#' Scale rule when the phase-1 size is fixed
#'
#' With `n_phase1` given, the phase-1 and phase-2 sizes can no longer be
#' traded against each other, so the multiplier is pinned by whichever
#' constraint is active rather than by the variance/cost product. The
#' relative allocation across strata is unaffected: it is
#' `S_h sqrt(d_2h / c_h)` in every mode.
#' @keywords internal
#' @noRd
.twophase_fixed_lambda <- function(W, S, cost, d2, target, what) {
  key <- ifelse(cost > 0, S * sqrt(d2 / cost), Inf)
  function(a_eff, c_eff, free) {
    if (identical(what, "budget")) {
      denom <- sum(cost[free] * W[free] * key[free])
      if (denom <= 0) return(Inf)
      (target - c_eff) / denom
    } else {
      denom <- target - a_eff
      if (denom <= 0) return(Inf)
      sum(W[free] * S[free] * sqrt(d2[free] * cost[free])) / denom
    }
  }
}

#' Whole-unit field design
#'
#' Budget mode floors and then buys back whole phase-2 units in order of
#' variance reduction per unit cost, so the integer design stays inside
#' the budget. CV mode rounds up. Both respect the phase-2 pool, which is
#' what stops either from reaching a stratum whose expected phase-1 yield
#' rounds to zero.
#' @keywords internal
#' @noRd
.twophase_integerize <- function(n1, nu, spec, phase1_cost, a_eff, d2_eff,
                                 nu_max, mode, budget, fixed_cost, N, mu,
                                 var_pop) {
  W <- spec$W
  measured <- spec$S > 0
  if (identical(mode, "budget")) {
    n1_int <- floor(min(n1, N))
    cap <- floor(nu_max * W * n1_int)
    n_int <- pmin(floor(nu * W * n1_int), cap)
    spent <- fixed_cost + phase1_cost * n1_int + sum(spec$cost * n_int)
    repeat {
      room <- n_int < cap & spec$cost > 0 & measured &
        spent + spec$cost <= budget
      if (!any(room)) break
      # a stratum left at zero contributes an infinite term, so its first
      # unit is worth more than any finite marginal gain elsewhere
      first <- room & n_int == 0
      score <- if (any(first)) {
        ifelse(first, d2_eff * W^2 * spec$S^2 / spec$cost, -Inf)
      } else {
        ifelse(room,
               d2_eff * W^2 * spec$S^2 *
                 (1 / n_int - 1 / (n_int + 1)) / spec$cost,
               -Inf)
      }
      h <- which.max(score)
      n_int[h] <- n_int[h] + 1
      spent <- spent + spec$cost[h]
    }
  } else {
    n1_int <- ceiling(min(n1, N))
    cap <- floor(nu_max * W * n1_int)
    n_int <- pmin(ceiling(nu * W * n1_int), cap)
    spent <- fixed_cost + phase1_cost * n1_int + sum(spec$cost * n_int)
  }
  starved <- measured & nu > 0 & n_int == 0
  if (any(starved)) {
    warning(sprintf(
      paste("the phase-1 pool rounds to zero units in %s, so the operational",
            "design cannot estimate %s and its 'cv' is Inf; enlarge phase 1,",
            "or collapse the stratum into a neighbour"),
      paste(sprintf("'%s'", spec$stratum[starved]), collapse = ", "),
      if (sum(starved) > 1L) "those strata" else "that stratum"),
      call. = FALSE)
  }
  # a stratum with no variance needs no sample; one that has variance and no
  # sample leaves the mean unestimable
  contrib <- ifelse(n_int > 0, d2_eff * W^2 * spec$S^2 * n1_int / n_int,
                    ifelse(measured, Inf, 0))
  var_unit <- a_eff + sum(contrib)
  list(
    n = c(n_phase1 = n1_int, n_phase2 = sum(n_int)),
    n_int = n_int,
    cost = spent,
    cv = .twophase_cv(n1_int, var_unit, mu, N, var_pop)
  )
}

#' Two-phase variance of the mean, with a component-wise correction
#'
#' `var_pop` is the census term: what the design would still carry if every
#' unit were measured. It holds the design effects, because one multiplies
#' an SRSWOR variance here as it does throughout the package, but not the
#' response divisors. Measuring every unit removes all the variance.
#' *issuing* to every unit and losing some of them to nonresponse does not,
#' and reusing `var_unit` as the census term would report that case as
#' exact.
#'
#' The phase-1 finite-population correction applies to what phase 1
#' estimates, not to the phase-2 residual. Writing \eqn{S^2} for the unit
#' population variance,
#' \deqn{V = A(1/n_a - 1/N)
#'       + \sum_h W_h S_h^2\{1/(\nu_h n_a) - 1/N\},}{V = A(1/n_a - 1/N) + sum_h W_h S_h^2\{1/(nu_h n_a) - 1/N\},}
#' which is `var_unit / n1 - var_pop / N`. Multiplying the whole variance
#' by \eqn{1 - n_a/N} instead subtracts \eqn{1/(\nu_h N)} where it should
#' subtract \eqn{1/N}, which understates the variance at every sampling
#' fraction and drives it to zero at a phase-1 census even though phase 2
#' still subsamples.
#' @keywords internal
#' @noRd
.twophase_var <- function(n1, var_unit, N, var_pop = var_unit) {
  if (is.infinite(N)) return(var_unit / n1)
  max(var_unit / n1 - var_pop / N, 0)
}

#' Coefficient of variation from the per-phase-1-unit variance
#' @keywords internal
#' @noRd
.twophase_cv <- function(n1, var_unit, mu, N, var_pop = var_unit) {
  if (!is.finite(mu) || mu == 0) return(NA_real_)
  sqrt(.twophase_var(n1, var_unit, N, var_pop)) / abs(mu)
}

#' The single-phase design the two-phase one has to beat
#'
#' Phase 1 is skipped, so the variable of interest is measured directly and
#' no screening cost is paid. Fuller (2009, Sect. 3.3.1) requires this
#' comparison because the interior optimum is not always an improvement.
#' @keywords internal
#' @noRd
.twophase_single <- function(spec, cv, budget, N, fixed_cost, single_cost,
                             single_deff, single_resp_rate) {
  # the same split as the two-phase design: response divides the sampling
  # term only, so issuing to the whole frame and losing units to nonresponse
  # still leaves variance behind
  total <- spec$A + sum(spec$W * spec$S^2)
  var_unit <- single_deff * total / single_resp_rate
  var_pop <- single_deff * total
  unit <- single_cost %||% sum(spec$W * spec$cost)
  if (!is.numeric(unit) || length(unit) != 1L || !is.finite(unit) || unit < 0) {
    stop("'single_cost' must be a non-negative finite scalar", call. = FALSE)
  }
  if (unit <= 0) {
    return(list(n = NA_real_, cv = NA_real_, cost = NA_real_, better = FALSE))
  }
  reaches <- TRUE
  if (!is.null(budget)) {
    n <- min((budget - fixed_cost) / unit, N)
  } else {
    if (!is.finite(spec$mu) || spec$mu == 0) {
      return(list(n = NA_real_, cv = NA_real_, cost = NA_real_, better = FALSE))
    }
    n <- var_unit / ((cv * spec$mu)^2 + var_pop / N)
    # a target under the comparator's own census floor asks for more units
    # than the frame holds; report what a census of it would achieve
    if (n > N) {
      n <- N
      reaches <- FALSE
    }
  }
  cost <- fixed_cost + n * unit
  list(n = n, cv = .twophase_cv(n, var_unit, spec$mu, N, var_pop), cost = cost,
       reaches_target = reaches, better = FALSE)
}

#' Precision of a two-phase allocation
#'
#' Achieved precision for a two-phase design you already have, the inverse
#' of [n_twophase()].
#'
#' @param frame For the default method: the phase-2 stratum frame described
#'   in [n_twophase()], with an extra `nu` column giving the subsampling
#'   fraction actually used in each stratum. For `svyplan_twophase`
#'   objects: a result from [n_twophase()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param n_phase1 Phase-1 sample size, counted as units issued. `nu` is a
#'   share of that same base, so with phase-1 nonresponse the largest
#'   usable `nu` is `resp_rate` rather than 1.
#' @param phase1_cost,between,mu,N,phase1_deff,resp_rate,fixed_cost,plan As
#'   in [n_twophase()]. The per-stratum phase-2 design effect and response
#'   rate travel in the frame's `deff` and `resp_rate` columns, as they do
#'   there.
#'
#' @param alpha Significance level for the reported margin of error. The
#'   default is 0.05. A two-phase design is sized against a cv, but the
#'   half-width it achieves follows from the standard error either way, so
#'   `$moe` is reported alongside `$cv`.
#'
#' @return A `svyplan_prec` object with `type = "twophase"`, carrying
#'   `$se`, `$moe`, `$rmoe`, `$cv`, and the per-stratum table in `$detail`.
#'   `$moe` is the `alpha`-level half-width `z * $se`, and `$rmoe` states it
#'   as a fraction of the population mean the design estimates.
#'
#' @family two-phase design functions
#' @seealso [n_twophase()] for the inverse.
#'
#' @examples
#' frame <- data.frame(
#'   stratum   = c("A", "B"),
#'   N         = c(6000, 4000),
#'   sd        = c(12, 25),
#'   mean      = c(40, 70),
#'   unit_cost = c(2, 5),
#'   nu        = c(0.4, 0.6)
#' )
#' prec_twophase(frame, n_phase1 = 2000, phase1_cost = 1)
#'
#' @export
prec_twophase <- function(frame, ...) {
  if (!missing(frame)) {
    .res <- .dispatch_plan(frame, "frame", prec_twophase.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_twophase")
}

#' @rdname prec_twophase
#' @export
prec_twophase.default <- function(
  frame,
  ...,
  n_phase1,
  phase1_cost = 1,
  between = NULL,
  mu = NULL,
  N = Inf,
  alpha = 0.05,
  phase1_deff = 1,
  resp_rate = 1,
  fixed_cost = 0,
  plan = NULL
) {
  merged <- .merge_plan_args(plan, prec_twophase.default, match.call(), environment())
  if (!is.null(merged)) {
    return(do.call(prec_twophase.default, c(merged, list(...))))
  }
  .check_unused_dots(...)
  if (!is.data.frame(frame) || !"nu" %in% names(frame)) {
    stop("'frame' must contain a 'nu' column giving the subsampling fraction",
         call. = FALSE)
  }
  nu <- frame$nu
  if (!is.numeric(nu) || anyNA(nu) || any(!is.finite(nu)) || any(nu <= 0) || any(nu > 1)) {
    stop("'nu' must contain values in (0, 1]", call. = FALSE)
  }
  spec <- .twophase_frame(frame, between, mu)
  if (missing(n_phase1)) {
    stop("'n_phase1' is required", call. = FALSE)
  }
  check_scalar(n_phase1, "n_phase1")
  check_alpha(alpha)
  .check_twophase_deff(phase1_deff, "phase1_deff")
  check_resp_rate(resp_rate)
  if (any(nu > resp_rate + sqrt(.Machine$double.eps))) {
    stop(sprintf(
      paste("'nu' is a share of the phase-1 units issued, so it cannot exceed",
            "the classified fraction 'resp_rate' (%s); phase 2 has no one else",
            "to draw from"),
      format(resp_rate)), call. = FALSE)
  }
  check_population_size(N)
  if (n_phase1 > N) {
    stop("'n_phase1' exceeds the population size 'N'", call. = FALSE)
  }

  a_eff <- phase1_deff * spec$A / resp_rate
  d2_eff <- spec$d2 / spec$r2
  var_unit <- a_eff + sum(d2_eff * spec$W * spec$S^2 / nu)
  # design effects, but not the response divisors: see .twophase_var()
  var_pop <- phase1_deff * spec$A + sum(spec$d2 * spec$W * spec$S^2)
  se <- sqrt(.twophase_var(n_phase1, var_unit, N, var_pop))
  detail <- data.frame(
    stratum = spec$stratum, N = spec$N, share = spec$W, sd = spec$S,
    unit_cost = spec$cost, deff = spec$d2, resp_rate = spec$r2, nu = nu,
    n_issued = nu * spec$W * n_phase1,
    n_resp = spec$r2 * nu * spec$W * n_phase1,
    take_all = spec$pinned,
    stringsAsFactors = FALSE
  )
  .new_svyplan_prec(
    se = se,
    moe = .q_alpha(alpha) * se,
    cv = .twophase_cv(n_phase1, var_unit, spec$mu, N, var_pop),
    type = "twophase",
    params = list(
      n_phase1 = n_phase1, phase1_cost = phase1_cost, between = spec$A,
      mu = spec$mu, N = N, alpha = alpha,
      phase1_deff = phase1_deff, resp_rate = resp_rate,
      take_all = spec$pinned,
      fixed_cost = fixed_cost,
      cost = fixed_cost + n_phase1 * (phase1_cost + sum(spec$cost * spec$W * nu))
    ),
    detail = detail
  )
}

#' @rdname n_twophase
#' @export
n_twophase.svyplan_prec <- function(frame, ..., cv = NULL, budget = NULL) {
  p <- frame$params
  if (!identical(frame$type, "twophase")) {
    stop("this is a '", frame$type, "' precision result, not a two-phase one",
         call. = FALSE)
  }
  base <- frame$detail
  if (is.null(base[["N"]])) base$N <- base$share
  # the pinned set the design was built under, not the strata the allocator
  # ended up truncating: feeding the latter back would over-constrain
  base$take_all <- p$take_all %||% base[["take_all"]]
  # the achieved precision is the default target, but only when the caller
  # has not asked for the other mode
  if (is.null(budget)) cv <- cv %||% frame$cv
  n_twophase(
    base[, setdiff(names(base), c("nu", "n", "n_issued", "n_resp")),
         drop = FALSE],
    phase1_cost = p$phase1_cost,
    cv = cv, budget = budget,
    between = p$between, mu = p$mu, N = p$N, phase1_deff = p$phase1_deff,
    resp_rate = p$resp_rate, fixed_cost = p$fixed_cost,
    n_phase1 = p$n_phase1_fixed, assurance = p$assurance,
    single_deff = p$single_deff %||% 1,
    single_resp_rate = p$single_resp_rate %||% 1,
    single_cost = p$single_cost, ...
  )
}

#' @rdname prec_twophase
#' @export
prec_twophase.svyplan_twophase <- function(frame, ...) {
  .check_unused_dots(...)
  p <- frame$params
  d <- frame$detail
  out <- prec_twophase(
    d[, c("stratum", "N", "sd", "unit_cost", "deff", "resp_rate", "nu")],
    n_phase1 = frame$n[["n_phase1"]], phase1_cost = p$phase1_cost,
    between = p$between, mu = p$mu, N = p$N, alpha = p$alpha %||% 0.05,
    phase1_deff = p$phase1_deff,
    resp_rate = p$resp_rate, fixed_cost = p$fixed_cost
  )
  # the design decisions a stratum table cannot express, so that
  # n_twophase() rebuilds this problem rather than an unconstrained one
  out$params$take_all <- p$take_all
  out$params$n_phase1_fixed <- p$n_phase1_fixed
  out$params$assurance <- p$assurance
  out$params$single_deff <- p$single_deff
  out$params$single_resp_rate <- p$single_resp_rate
  out$params$single_cost <- p$single_cost
  out$detail$take_all <- p$take_all %||% d$take_all
  out
}

#' Validate a scalar two-phase design effect
#'
#' Mirrors `check_deff()`, which the rest of the package uses for the scalar
#' `deff` argument. Per-stratum variation travels in the frame instead, as
#' it does for `unit_cost`.
#' @keywords internal
#' @noRd
.check_twophase_deff <- function(deff, name) {
  if (!is.numeric(deff) || length(deff) != 1L || is.na(deff) ||
      !is.finite(deff) || deff <= 0) {
    stop(sprintf("'%s' must be a positive finite scalar", name), call. = FALSE)
  }
  invisible(TRUE)
}
