#' Planning design effect
#'
#' Build the design effect you expect a complex design to produce, before
#' any data are collected, by combining the design features you are
#' planning: clustering, unequal weighting, stratification, and optional
#' auxiliary calibration. The result
#' is a multiplier you pass as `deff` to [n_prop()], [n_mean()],
#' [n_alloc()], or any other sizing or precision function.
#'
#' @param x A `svyplan` result to read design features from: a
#'   [n_cluster()] or [prec_cluster()] allocation, a [varcomp()] estimate,
#'   or an [n_alloc()] or [prec_alloc()] allocation. `NULL` (default) builds
#'   the design effect from the component arguments below.
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#'
#' @return A `svyplan_deff` object: a numeric scalar carrying the component
#'   decomposition. Use it directly wherever a `deff` argument is expected,
#'   [as.double()] to strip it to a plain number, and [as.data.frame()] to
#'   export the components.
#'
#' @details
#' The design effect (DEFF) is the ratio of the variance under the planned
#' design to the variance of a simple random sample of the same size.
#' DEFF = 2 means you need twice the sample for the same precision, whereas
#' DEFF below 1 means the design is more efficient than simple random
#' sampling.
#'
#' `design_effect()` is a **planning** tool. It anticipates a design effect
#' from design parameters you choose or estimate beforehand. It does not
#' estimate a design effect from collected data. Once you have a realized
#' sample, use `survey::svymean(..., deff = TRUE)` and `survey::deff()`,
#' which compute it from the actual weights, strata, and clusters.
#'
#' ## Components
#'
#' Components are selected by which arguments you supply, and multiply
#' together:
#'
#' \deqn{DEFF = DEFF_{cluster} \times DEFF_{weight} \times DEFF_{strata} \times DEFF_{calibration}.}{DEFF = DEFF_cluster * DEFF_weight * DEFF_strata * DEFF_calibration.}
#'
#' **Clustering** (`icc`, `take`, `var_ratio`) uses the same variance model as
#' [n_cluster()] and [prec_cluster()], so the two always agree. For an
#' \eqn{S}-stage design with stage homogeneities \eqn{\delta_1, \dots, \delta_{S-1}}{delta_1, ..., delta_{S-1}},
#' variance ratios \eqn{k_1, \dots, k_{S-1}}{k_1, ..., k_{S-1}}, and per-parent
#' takes \eqn{t_2, \dots, t_S}{t_2, ..., t_S} (supplied in `take`, length
#' \eqn{S - 1}),
#' \deqn{DEFF_{cluster} = \sum_{s=1}^{S-1} k_s \delta_s \prod_{j=s+1}^{S} t_j + k_{S-1}(1 - \delta_{S-1}).}{DEFF_cluster = sum_{s=1}^{S-1} k_s delta_s prod_{j=s+1}^S t_j + k_{S-1}(1 - delta_{S-1}).}
#' For two stages (\eqn{t_2 = m}), this is \eqn{k_1(1 + \delta_1(m - 1))}. For
#' three stages (\eqn{t_2 = m, t_3 = q}), it is
#' \eqn{k_1 \delta_1 m q + k_2(1 + \delta_2(q - 1))}. Estimate `icc` (and
#' `var_ratio`) from a previous round or pilot with [varcomp()].
#'
#' In the multistage form \eqn{k_s} rescales the stage-\eqn{s} components to
#' the unit variance of the analysis variable. When sub-units have equal size
#' or are drawn proportional to size, each stage ratio follows recursively from
#' the stage above:
#' \deqn{k_s = k_{s-1}(1 - \delta_{s-1}), \quad s = 2, \dots, S - 1.}{k_s = k_{s-1}(1 - delta_{s-1}), s = 2, ..., S - 1.}
#' This makes the design effect collapse to \eqn{k_1} when every take equals 1,
#' where no clustering remains. A scalar `var_ratio` supplies \eqn{k_1} and
#' derives deeper ratios this way. When unequal-sized sub-units are drawn by
#' simple random sampling, supply the full length-\eqn{S - 1} vector (for
#' instance from [varcomp()]).
#'
#' Written on components referenced to the total unit variance,
#' \eqn{\delta_1=\sigma_1^2/S^2} and \eqn{\delta_2^{tot}=\sigma_2^2/S^2}{delta_2^tot=sigma_2^2/S^2},
#' the same quantity is the familiar
#' \eqn{1+\delta_1(mq-1)+\delta_2^{tot}(q-1)}{1+delta_1(mq-1)+delta_2^tot(q-1)}. The package's `icc_2` is
#' referenced to the within-PSU variance instead, so
#' \eqn{\delta_2=\delta_2^{tot}/(1-\delta_1)}{delta_2=delta_2^tot/(1-delta_1)}.
#'
#' **Unequal weighting** (`weights`, or `N` and `n` in `strata`) is Kish's
#' weighting loss. From a vector of planned weights,
#' \deqn{DEFF_{weight} = n\sum_i w_i^2 / (\sum_i w_i)^2,}{DEFF_weight = n sum_i w_i^2 / (sum_i w_i)^2,}
#' and from a planned stratified allocation with stratum sizes \eqn{N_h}
#' and takes \eqn{n_h} the same quantity is
#' \deqn{DEFF_{weight} = n\sum_h N_h^2/n_h / (\sum_h N_h)^2, \quad
#'       n = \sum_h n_h.}{DEFF_weight = n sum_h N_h^2/n_h / (sum_h N_h)^2, n = sum_h n_h.}
#' Equal weights give exactly 1. This is the cost of a disproportionate
#' allocation, of weighting classes, or of an anticipated nonresponse
#' adjustment, and it is never below 1. The same factor, computed from a
#' past survey's final weights, divides that survey's weighting loss out of
#' its published design effect. See [varcomp()].
#'
#' **Stratification** (`N`, `sd`, and `mean` in `strata`) is the gain from
#' stratifying, the only design-stage component that can fall below 1. With
#' \eqn{W_h = N_h/N},
#' \deqn{DEFF_{strata} = \sum_h W_hS_h^2 /
#'       \left(\sum_h W_hS_h^2 + \sum_h W_h(\bar y_h - \bar y)^2\right).}{DEFF_strata = sum_h W_hS_h^2 / (sum_h W_hS_h^2 + sum_h W_h(ybar_h - ybar)^2 ).}
#' It measures the gain under *proportional* allocation. Any departure from
#' proportional is already charged to the weighting component, so the two
#' compose without double counting.
#'
#' **Auxiliary calibration** (`aux_r2`) approximates the variance reduction
#' from GREG or linear calibration on auxiliary variables with anticipated
#' squared correlation \eqn{R^2 \in [0, 1)}{R^2 in [0, 1)}, multiplying the
#' design effect by \eqn{1 - R^2}. Do not supply `aux_r2` when the `sd`, `var`,
#' or `icc` inputs were already computed from regression residuals rather than
#' the raw survey variable, as that would double-count the auxiliary gain.
#'
#' Because the clustering component already prices the clustering, do not
#' additionally inflate `deff` for it elsewhere. Components you leave
#' unspecified are simply absent, which is equivalent to setting them to 1.
#'
#' ## What a multiplier can and cannot say
#'
#' Multiplying the components you supply by hand is Kish's approximation. It
#' is exact when the stratum standard deviations are equal, whatever the
#' stratum sizes, means and takes, and not otherwise. The gap runs in
#' either direction, so it is not a bound. The weighting component charges
#' a disproportionate allocation as a loss without knowing which strata
#' were favored, so a Neyman allocation over strata that differ sharply in
#' `sd` is reported well above its true variance ratio, while a constraint
#' that forces units into a low-variance stratum is reported well below it.
#' The printed result marks a multi-component product as approximate for
#' that reason.
#'
#' Given an [n_alloc()] result there is no need to approximate at all.
#' Every quantity the ratio needs is already there, so
#' `design_effect(alloc)` returns the allocation's own variance ratio,
#' \deqn{D = \sum_h n_hc_h \sum_h W_h^2S_h^2d_hk_h(1+\delta_h(m_hr_h-1))/(n_hc_h) \Big/
#'        \left(\sum_h W_hS_h^2+\sum_h W_h(\bar y_h-\bar y)^2\right),}{D = sum_h n_h c_h sum_h W_h^2S_h^2d_hk_h(1+delta_h(m_h r_h-1))/(n_h c_h) / (sum_h W_hS_h^2+sum_h W_h(ybar_h-ybar)^2 ),}
#' with the per-stratum cluster factors entering stratum by stratum rather
#' than averaged, and the `deff` the allocation was built under counted
#' once. \eqn{n_hc_h}{n_h c_h} counts the units of stratum \eqn{h} that
#' respond and are in scope, \eqn{c_h}{c_h} being the product of its unit
#' response, PSU response and eligibility rates, since those are the units
#' the plan's variance reads. Where every stratum has the same rate it
#' cancels. Where the rates differ, a stratum that responds poorly carries
#' more of the variance than its issued share suggests, and the ratio
#' includes that loss. A per-stratum `deff` weights its own stratum's contribution and a
#' scalar factors straight back out. Joint designs from [n_alloc_joint()]
#' and their assessments can have different variance ratios by measure.
#' Use [prec_alloc_joint()] for measure-specific precision.
#'
#' The clustering bracket reads \eqn{m_hr_h}, the take that responds, because
#' that is the take [n_alloc()] sized the design on. So the ratio reported
#' here falls as the ultimate-unit response rate falls, and describes the
#' variance model the allocation actually used rather than the one its gross
#' take would suggest. PSU-level response is absent from the bracket,
#' because it removes whole clusters instead of shrinking the ones that
#' are worked.
#'
#' This ratio is computed without finite population corrections, on the
#' same planning scale as the component arguments. [prec_alloc()] applies
#' each stratum's FPC, so a ratio derived from its variance need not equal
#' `design_effect(alloc)` once sampling fractions are material. Use
#' [prec_alloc()] when the design variance itself is what you want.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer. Chapter 9.
#'
#' @family design input functions
#' @seealso [effective_n()] for the sample size this design effect costs
#'   you, [varcomp()] for estimating `icc` and `var_ratio`, [n_cluster()] for
#'   optimizing a multistage design directly.
#'
#' @examples
#' # Clustering: 25 households per cluster, homogeneity 0.05
#' design_effect(icc = 0.05, take = 25)
#'
#' # The same cluster design with an anticipated calibration gain
#' design_effect(icc = 0.05, take = 25, aux_r2 = 0.25)
#'
#' # One design effect per stratum, from a common icc and the takes. The
#' # result of each call is a scalar, so sapply() returns a plain vector.
#' sapply(c(20, 25, 30), function(m) design_effect(icc = 0.03, take = m))
#'
#' # Three stages: 10 SSUs per PSU, 4 units per SSU
#' design_effect(icc = c(0.01, 0.05), take = c(10, 4))
#'
#' # Four stages: 5 EAs per district, 4 segments per EA, 8 households per segment
#' design_effect(icc = c(0.01, 0.03, 0.06), take = c(5, 4, 8))
#'
#' # Weighting loss from a planned set of weights
#' design_effect(weights = rep(c(1, 4), c(300, 100)))
#'
#' # A planned stratified allocation: weighting loss and stratification gain
#' frame <- data.frame(
#'   N    = c(50000, 120000),
#'   n    = c(600, 400),
#'   sd   = c(12, 20),
#'   mean = c(55, 48)
#' )
#' design_effect(strata = frame)
#'
#' # Clustering and weighting together
#' deff <- design_effect(icc = 0.05, take = 25, strata = frame)
#' deff
#' n_prop(p = 0.3, moe = 0.05, deff = deff)
#'
#' # Read the design features off an existing plan
#' plan <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
#' design_effect(plan)
#'
#' @name design_effect
#' @export
design_effect <- function(x = NULL, ...) {
  # With no object to dispatch on, UseMethod() would dispatch on the first
  # element of ..., which here is a planning component, not an object.
  if (missing(x) || is.null(x)) {
    return(design_effect.default(NULL, ...))
  }
  UseMethod("design_effect")
}

#' @describeIn design_effect Build the design effect from planning
#'   components. Supply the components you are designing for and leave the
#'   rest `NULL`.
#'
#' @param icc Homogeneity of units within a cluster, in \eqn{[0, 1]}.
#'   A vector of length `stages - 1` (scalar for two stages, length 2 for
#'   three stages, and so on). This is the survey-planning measure returned by
#'   [varcomp()], not a generic mixed-model ICC.
#' @param take Per-parent take at stages `2..S` (length `stages - 1`, matching
#'   `icc`). Scalar for two stages (units per PSU), length 2 for three stages
#'   (SSUs per PSU, units per SSU), and so on. Required with `icc`.
#' @param var_ratio Ratio of the components' unit variance to the analysis
#'   variable's, default 1. A scalar supplies stage 1 and derives deeper ratios
#'   via `k_s = k_{s-1} * (1 - icc_{s-1})`. Supply a vector matching `icc`
#'   when unequal-sized sub-units are drawn by SRS. See Details.
#' @param weights Numeric vector of planned sampling weights, one per
#'   sampled unit. Only their relative variability matters. On an
#'   [n_alloc()] result they are an *additional* anticipated weighting
#'   adjustment charged on top of the allocation, such as a nonresponse
#'   correction or a calibration step. Do not pass the allocation's own
#'   weights there. The stratum takes already carry those, and supplying
#'   them again counts the same disproportionality twice. The `allocation`
#'   component is exact for the plan, while the weighting one is Kish's
#'   approximation, so with `weights` supplied the overall result is
#'   approximate too and prints as such.
#' @param strata Stratum-level data frame describing a planned stratified
#'   allocation, using the [n_alloc()] column names. `N` with `n` gives
#'   the weighting component. `N` with `sd` and `mean` gives the
#'   stratification component. Supply all four for both.
#' @param aux_r2 Optional anticipated squared correlation \eqn{R^2} in
#'   \eqn{[0, 1)} between the survey variable and auxiliary calibration
#'   variables. Multiplies the design effect by \eqn{1 - R^2}. Avoid supplying
#'   this when variance or ICC inputs are already computed on residuals.
#'
#' @export
design_effect.default <- function(
  x = NULL,
  ...,
  icc = NULL,
  take = NULL,
  var_ratio = 1,
  weights = NULL,
  strata = NULL,
  aux_r2 = NULL
) {
  .check_unused_dots(...)
  if (!is.null(x)) {
    stop(
      "'x' must be a svyplan result or NULL. Pass planning components by name (icc, take, weights, strata, aux_r2)",
      call. = FALSE
    )
  }
  .deff_compose(
    icc = icc,
    take = take,
    var_ratio = var_ratio,
    weights = weights,
    strata = strata,
    aux_r2 = aux_r2
  )
}

#' @describeIn design_effect Reject a bare numeric vector, which says
#'   nothing about which design component it describes.
#' @export
design_effect.numeric <- function(x, ...) {
  stop(
    "design_effect() takes its planning components by name. Use design_effect(weights = w) for the weighting loss of a planned weight vector, or survey::svymean(deff = TRUE) to estimate a design effect from collected data",
    call. = FALSE
  )
}

#' @describeIn design_effect Read `icc`, `var_ratio`, and the stage takes from a
#'   [n_cluster()] or [prec_cluster()] allocation. The takes are those that
#'   respond (including `eligibility_rate` and `eligible_per_unit` at the final
#'   stage), so the design effect and the `cv` describe one design. A plan
#'   sized from an `indicators` table has no single design effect and is
#'   refused.
#' @export
design_effect.svyplan_cluster <- function(x, ..., weights = NULL,
                                          strata = NULL, aux_r2 = NULL) {
  .check_unused_dots(...)
  if (!is.null(x$indicators)) {
    stop(
      paste("this plan sizes an 'indicators' table, where each row carries",
            "its own icc and each domain its own takes, so it has no",
            "single design effect; build one row's with design_effect(icc =",
            "..., take = ...) from its icc and the takes in $operational"),
      call. = FALSE
    )
  }
  take <- .deff_cluster_takes(x$n, x$stages, x$params)
  .deff_compose(
    icc = x$params$icc,
    take = take,
    var_ratio = x$params$var_ratio %||% 1,
    weights = weights,
    strata = strata,
    aux_r2 = aux_r2
  )
}

#' @describeIn design_effect Read the design features off a
#'   [prec_cluster()] result, the same way as for [n_cluster()], or the
#'   stratum counts of a [prec_alloc()] result, the same way as for
#'   [n_alloc()]. Other `svyplan_prec` types carry the `deff` you supplied
#'   rather than one to be derived, and are rejected.
#' @export
design_effect.svyplan_prec <- function(x, ..., weights = NULL,
                                       strata = NULL, aux_r2 = NULL) {
  .check_unused_dots(...)
  if (x$type %in% c("alloc", "alloc_joint")) {
    .deff_no_strata(strata)
    return(.deff_alloc(x, weights, aux_r2 = aux_r2))
  }
  if (!identical(x$type, "cluster")) {
    stop(
      sprintf(
        "design_effect() reads a cluster design or an allocation; this is a '%s' result. Its design effect is the 'deff' you supplied",
        x$type
      ),
      call. = FALSE
    )
  }
  if (!is.null(x$params$indicators)) {
    stop(
      "design_effect() requires a single-indicator cluster result; this result has multiple indicators with different ICCs",
      call. = FALSE
    )
  }
  take <- .deff_cluster_takes(x$params$n, x$params$stages, x$params)
  .deff_compose(
    icc = x$params$icc,
    take = take,
    var_ratio = x$params$var_ratio %||% 1,
    weights = weights,
    strata = strata,
    aux_r2 = aux_r2
  )
}

#' @describeIn design_effect Read `icc` and `var_ratio` from a [varcomp()]
#'   estimate. The stage takes are your design choice, so `take` (length
#'   `stages - 1`) must be supplied.
#' @export
design_effect.svyplan_varcomp <- function(x, ..., take = NULL, weights = NULL,
                                          strata = NULL, aux_r2 = NULL) {
  .check_unused_dots(...)
  if (!is.null(x$strata)) {
    stop(
      "stratified varcomp: pass one stratum's icc and var_ratio (see the $strata table)",
      call. = FALSE
    )
  }
  if (identical(x$source, "deff") && is.null(take)) {
    # The source take is provenance, not a default: re-planning at a
    # different take is the reason to back an icc out at all, so falling
    # back to it would quietly return the previous design's number.
    stop(
      sprintf(
        "'take' is required. This 'icc' was backed out of a design effect at take = %.4g; supply the take you are planning for",
        x$params$take %||% x$params$n_per_psu
      ),
      call. = FALSE
    )
  }
  .deff_compose(
    icc = x$icc,
    take = take,
    var_ratio = x$var_ratio,
    weights = weights,
    strata = strata,
    aux_r2 = aux_r2
  )
}

#' @describeIn design_effect Read the planned allocation from an
#'   [n_alloc()] result. This returns the allocation's own variance ratio,
#'   computed from the stratum sizes, takes, standard deviations, means and
#'   any per-stratum cluster factors, rather than a product of separate
#'   weighting and stratification approximations. An allocation from stratum
#'   sizes alone has no variance to compare, and returns Kish's weighting
#'   loss of its whole takes that respond,
#'   \eqn{\sum_h n_hc_h \sum_h W_h^2 / (n_hc_h)}{sum(n_h c_h) sum(W_h^2 / (n_h c_h))},
#'   the ratio when every stratum has the same variance. Supply `weights` only for
#'   an *additional* anticipated adjustment. See the `weights` argument.
#' @export
design_effect.svyplan_alloc <- function(x, ..., weights = NULL, aux_r2 = NULL) {
  .check_unused_dots(...)
  .deff_alloc(x, weights, aux_r2 = aux_r2)
}

#' @describeIn design_effect Single-indicator sizing results carry the `deff`
#'   supplied at creation rather than deriving one, and are rejected.
#' @export
design_effect.svyplan_n <- function(x, ..., weights = NULL) {
  .check_unused_dots(...)
  stop(
    sprintf(
      "design_effect() reads a stratified allocation; this is a '%s' result. Its design effect is the 'deff' you supplied",
      x$type
    ),
    call. = FALSE
  )
}

#' Design effect of a one-measure allocation at its own stratum rates
#' @keywords internal
#' @noRd
.deff_alloc <- function(x, weights = NULL, aux_r2 = NULL) {
  if (!x$type %in% c("alloc", "alloc_joint")) {
    stop(
      sprintf(
        "design_effect() reads a stratified allocation; this is a '%s' result. Its design effect is the 'deff' you supplied",
        x$type
      ),
      call. = FALSE
    )
  }
  if (identical(x$type, "alloc_joint") || !is.null(x$params$measures)) {
    stop(
      paste("a joint allocation optimizes several measures at once, and",
            "they do not share one variance ratio; take the measure-specific",
            "precision from prec_alloc_joint() instead"),
      call. = FALSE
    )
  }
  if (all(is.na(x$detail[["sd"]]))) {
    # A size-only allocation identifies no variance ratio. Its weights are
    # all it has, and under equal stratum variances their loss is the ratio.
    rates <- .alloc_summary_response_rates(x)
    loss <- .deff_from_allocation(
      x$detail$N * rates$eligibility,
      (x$detail[["n_int"]] %||% x$detail$n) * rates$combined
    )
    components <- c(weighting = loss$value)
    notes <- c(weighting = paste0(
      loss$note, ", whole takes, equal stratum variances assumed (Kish)"
    ))
  } else {
    ratio <- .deff_alloc_ratio(x)
    components <- c(allocation = ratio$value)
    notes <- c(allocation = ratio$note)
  }
  if (!is.null(weights)) {
    extra <- .deff_from_weights(weights)
    components <- c(components, weight = extra$value)
    notes <- c(notes, weight = paste0(extra$note, " (Kish)"))
  }
  if (!is.null(aux_r2)) {
    cal <- .deff_from_calibration(aux_r2)
    components <- c(components, calibration = cal$value)
    notes <- c(notes, calibration = cal$note)
  }
  .new_svyplan_deff(prod(components), components = components, notes = notes)
}

#' Refuse a strata table on an allocation, which carries its own strata
#' @keywords internal
#' @noRd
.deff_no_strata <- function(strata) {
  if (!is.null(strata)) {
    stop(
      "an allocation carries its own strata; drop 'strata'",
      call. = FALSE
    )
  }
}

#' Variance ratio of a one-measure stratified allocation
#'
#' The design variance of the mean under this allocation, over the SRS
#' variance at the same total size. Both sizes count the units that respond
#' and are in scope, `n_h` times the stratum's combined rate, as the plan's
#' variance does.
#' Everything it needs is already in the result, so there is no reason to
#' approximate it by multiplying a Kish weighting factor by a
#' proportional-allocation stratification factor: that product is exact only
#' when the stratum standard deviations agree, and it is not even an upper
#' bound when a constraint forces the allocation to oversample a
#' low-variance stratum.
#'
#' Computed without finite population corrections, on the same scale as the
#' component API. `prec_alloc()` applies each stratum's FPC, so a ratio
#' derived from its variance need not agree when sampling fractions are
#' material.
#' @keywords internal
#' @noRd
.deff_alloc_ratio <- function(x) {
  detail <- x$detail
  missing_cols <- setdiff(c("N", "n", "sd", "mean"), names(detail))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        paste("this allocation has no %s column, so the population variance",
              "the design effect divides by is not identified; add %s to the",
              "n_alloc() frame, or use prec_alloc() for the variance itself"),
        paste(sQuote(missing_cols), collapse = " or "),
        paste(sQuote(missing_cols), collapse = " and ")
      ),
      call. = FALSE
    )
  }
  rates <- .alloc_summary_response_rates(x)
  W <- detail$N * rates$eligibility / sum(detail$N * rates$eligibility)
  n <- detail$n
  factor_h <- .deff_alloc_factor(x$params$frame, detail, rates$stage_mat)
  parts <- .deff_variance_parts(W, detail$sd, detail$mean)
  if (parts$within <= 0) {
    stop(paste("the allocation has no within-stratum variability to summarize:",
               "every 'sd' is zero, so the ratio is 0 or undefined"),
         call. = FALSE)
  }
  # the deff the allocation was built under is part of its variance model, so
  # it belongs in the ratio exactly once. Per-stratum values weight their own
  # stratum's contribution; a scalar factors straight back out.
  deff_h <- x$params$deff %||% 1
  if (length(deff_h) > 1L && length(deff_h) != nrow(detail)) {
    stop("'deff' length does not match the allocation's strata", call. = FALSE)
  }
  deff_h <- .alloc_resolve_h(deff_h, x$params$frame, "deff", nrow(detail),
                             .check_deff_h)
  m <- n * rates$combined
  value <- sum(m) * sum(W^2 * detail$sd^2 * deff_h * factor_h / m) /
    parts$total
  clustered <- is.data.frame(x$params$frame) &&
    .alloc_is_cluster(x$params$frame) &&
    !is.null(detail[["take_2"]])
  qualifications <- c(
    if (clustered) "per-stratum clustering",
    if (!isTRUE(all.equal(unname(deff_h), rep(1, length(deff_h))))) {
      sprintf("input deff = %s", paste(signif(deff_h, 4), collapse = ", "))
    }
  )
  list(
    value = value,
    note = sprintf(
      "%d strata, n = %.4g%s (direct variance ratio)",
      nrow(detail), sum(n),
      if (length(qualifications) > 0L) {
        paste0(", ", paste(qualifications, collapse = ", "))
      } else {
        ""
      }
    )
  )
}

#' Per-stratum clustering factor of an allocation, or one where there is none
#'
#' The factors enter the variance ratio stratum by stratum, so unlike the
#' component API there is nothing to average here.
#'
#' `unit_resp_rate` is the ultimate-unit response rate, and the bracket reads
#' the take that responds, `take_2 * unit_resp_rate`. This is the same
#' quantity the allocation itself was built on, so the ratio reported here
#' describes the variance model the design actually used. Reading the gross
#' take instead would report a clustering penalty the design never pays. PSU
#' response is deliberately absent: it removes whole clusters rather than
#' shrinking the ones that are worked, so it does not enter the bracket.
#' @keywords internal
#' @noRd
.deff_alloc_factor <- function(frame, detail, stage_resp_mat = NULL) {
  H <- nrow(detail)
  if (!is.data.frame(frame) || !.alloc_is_cluster(frame) ||
      is.null(detail[["take_2"]]) || nrow(frame) != H) {
    return(rep(1, H))
  }
  S <- .alloc_detect_stages(frame)
  if (is.null(stage_resp_mat)) {
    stage_resp_mat <- matrix(1, nrow = H, ncol = S)
  }
  out <- numeric(H)
  for (h in seq_len(H)) {
    icc_h <- vapply(seq_len(S - 1L), function(s) {
      as.numeric(frame[[paste0("icc_", s)]][h])
    }, numeric(1L))
    vr_h <- numeric(S - 1L)
    for (s in seq_len(S - 1L)) {
      col <- paste0("var_ratio_", s)
      if (col %in% names(frame)) {
        vr_h[s] <- as.numeric(frame[[col]][h])
      } else if (s == 1L) {
        vr_h[1L] <- 1
      } else {
        vr_h[s] <- vr_h[s - 1L] * (1 - icc_h[s - 1L])
      }
    }
    take_h <- vapply(2L:S, function(s) {
      as.numeric(detail[[paste0("take_", s)]][h]) * stage_resp_mat[h, s]
    }, numeric(1L))
    out[h] <- .stage_deff_cluster(icc_h, take_h, vr_h)
  }
  out
}

#' Combine the requested planning components into a svyplan_deff
#' @keywords internal
#' @noRd
.deff_compose <- function(icc, take, var_ratio, weights, strata, aux_r2 = NULL) {
  parts <- numeric(0)
  notes <- character(0)

  if (!is.null(icc) || !is.null(take)) {
    cl <- .deff_from_cluster(icc, take, var_ratio)
    parts[["cluster"]] <- cl$value
    notes[["cluster"]] <- cl$note
  }

  strata <- .deff_check_strata(strata)
  if (!is.null(weights) && !is.null(strata) && !is.null(strata$n)) {
    stop(
      "supply the weighting component once: either 'weights' or 'n' in 'strata'",
      call. = FALSE
    )
  }

  if (!is.null(weights)) {
    wt <- .deff_from_weights(weights)
    parts[["weight"]] <- wt$value
    notes[["weight"]] <- wt$note
  } else if (!is.null(strata) && !is.null(strata$n)) {
    wt <- .deff_from_allocation(strata$N, strata$n)
    parts[["weight"]] <- wt$value
    notes[["weight"]] <- wt$note
  }

  if (!is.null(strata) && !is.null(strata$sd)) {
    st <- .deff_from_strata(strata$N, strata$sd, strata$mean)
    parts[["strata"]] <- st$value
    notes[["strata"]] <- st$note
  }

  if (!is.null(aux_r2)) {
    cal <- .deff_from_calibration(aux_r2)
    parts[["calibration"]] <- cal$value
    notes[["calibration"]] <- cal$note
  }

  if (length(parts) == 0L) {
    stop(
      "no design components supplied. Give at least one of 'icc' with 'take', 'weights', 'strata', or 'aux_r2'",
      call. = FALSE
    )
  }

  .new_svyplan_deff(prod(parts), components = parts, notes = notes)
}

#' Calibration variance reduction factor from auxiliary R^2
#' @keywords internal
#' @noRd
.deff_from_calibration <- function(aux_r2) {
  if (!is.numeric(aux_r2) || length(aux_r2) != 1L || is.na(aux_r2) ||
      !is.finite(aux_r2) || aux_r2 < 0 || aux_r2 >= 1) {
    stop("'aux_r2' must be a single finite number in [0, 1)", call. = FALSE)
  }
  val <- 1 - as.numeric(aux_r2)
  list(
    value = val,
    note = sprintf("aux_r2 = %.4g (factor 1 - R^2 = %.4g)", as.numeric(aux_r2), val)
  )
}

#' Clustering component for arbitrary S >= 2, matching prec_cluster()
#' @keywords internal
#' @noRd
.deff_from_cluster <- function(icc, take, var_ratio) {
  if (is.null(icc) || is.null(take)) {
    stop(
      "the clustering component needs both 'icc' and 'take'",
      call. = FALSE
    )
  }
  icc <- .reorder_stage_vec(icc, "icc")
  var_ratio <- .reorder_stage_vec(var_ratio, "var_ratio")
  check_icc(icc)
  k_len <- length(icc)
  if (!is.numeric(take) || length(take) != k_len || anyNA(take) ||
      any(!is.finite(take)) || any(take < 1)) {
    if (!is.numeric(take) || anyNA(take) || any(!is.finite(take)) ||
        any(take < 1)) {
      stop("'take' must contain finite numbers >= 1", call. = FALSE)
    }
    stop(
      sprintf("'take' must have length %d (matching 'icc')", k_len),
      call. = FALSE
    )
  }
  if (!is.numeric(var_ratio) || !length(var_ratio) %in% c(1L, k_len) ||
      anyNA(var_ratio) || any(!is.finite(var_ratio)) || any(var_ratio <= 0)) {
    stop(
      sprintf(
        "'var_ratio' must contain 1 or %d positive finite value(s)", k_len
      ),
      call. = FALSE
    )
  }
  var_ratio <- .expand_var_ratio(var_ratio, icc)
  value <- .stage_deff_cluster(icc, take, var_ratio)
  note <- if (k_len == 1L) {
    sprintf("icc = %.4g, take = %.4g, var_ratio = %.4g", icc, take, var_ratio)
  } else {
    sprintf(
      "%d-stage, icc = (%s), take = (%s), var_ratio = (%s)",
      k_len + 1L,
      paste(sprintf("%.4g", icc), collapse = ", "),
      paste(sprintf("%.4g", take), collapse = ", "),
      paste(sprintf("%.4g", var_ratio), collapse = ", ")
    )
  }
  list(value = as.numeric(value), note = note)
}

#' Weighting component from a vector of planned weights
#' @keywords internal
#' @noRd
.deff_from_weights <- function(w) {
  check_weights(w, "weights")
  n <- length(w)
  ## Scale-free, and dividing by the largest keeps w^2 inside the doubles.
  w <- w / max(w)
  value <- n * sum(w^2) / sum(w)^2
  list(
    value = value,
    note = sprintf("%d planned weights, cv(w) = %.4g", n,
                   sqrt(max(value - 1, 0)))
  )
}

#' Weighting component implied by a stratified allocation
#'
#' Equivalent to .deff_from_weights() on the weights N_h / n_h replicated
#' n_h times, without materializing them.
#' @keywords internal
#' @noRd
.deff_from_allocation <- function(N, n) {
  total_n <- sum(n)
  value <- total_n * sum(N^2 / n) / sum(N)^2
  list(
    value = value,
    note = sprintf("%d strata, n = %.4g", length(N), total_n)
  )
}

#' Split a stratified population variance into its two parts
#'
#' The within term is a sum of squares, so it is zero only when every `sd`
#' is. The between term subtracts the overall mean from each stratum mean,
#' and equal means cancel to rounding noise rather than to zero, so the
#' term is compared against the size of that noise instead of against an
#' absolute threshold. A variance expressed in small units is still a
#' variance: nothing here may depend on the unit `sd` and `mean` are
#' measured in.
#' @keywords internal
#' @noRd
.deff_variance_parts <- function(share, sd, mean) {
  within <- sum(share * sd^2)
  deviation <- mean - sum(share * mean)
  between <- sum(share * deviation^2)
  noise <- (.Machine$double.eps * max(abs(mean)))^2
  if (between <= noise) {
    between <- 0
  }
  list(within = within, between = between, total = within + between)
}

#' Stratification component under proportional allocation
#' @keywords internal
#' @noRd
.deff_from_strata <- function(N, sd, mean) {
  if (is.null(mean)) {
    stop(
      "the stratification component needs 'mean' alongside 'N' and 'sd' in 'strata'",
      call. = FALSE
    )
  }
  share <- N / sum(N)
  parts <- .deff_variance_parts(share, sd, mean)
  if (parts$within <= 0) {
    stop(
      paste("'strata' has no within-stratum variability: every 'sd' is zero,",
            "so the stratification component is 0 or undefined"),
      call. = FALSE
    )
  }
  list(
    value = parts$within / parts$total,
    note = sprintf("%d strata, between-stratum share %.4g", length(N),
                   parts$between / parts$total)
  )
}

#' Validate a stage take
#' @keywords internal
#' @noRd
.deff_check_take <- function(size, name) {
  check_scalar(size, name)
  if (size < 1) {
    stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate the stratum table and return only the usable columns
#' @keywords internal
#' @noRd
.deff_check_strata <- function(strata) {
  if (is.null(strata)) {
    return(NULL)
  }
  if (!is.data.frame(strata) || nrow(strata) == 0L) {
    stop("'strata' must be a non-empty data frame", call. = FALSE)
  }
  if (is.null(strata[["N"]])) {
    stop("'strata' must have an 'N' column of stratum population sizes",
         call. = FALSE)
  }
  out <- list(N = .deff_check_column(strata$N, "N", positive = TRUE))
  if (!is.null(strata[["n"]])) {
    out$n <- .deff_check_column(strata$n, "n", positive = TRUE)
    if (any(out$n > out$N)) {
      stop("'strata$n' cannot exceed 'strata$N'", call. = FALSE)
    }
  }
  if (!is.null(strata[["sd"]]) && !is.null(strata[["var"]])) {
    stop("'strata' gives both 'sd' and 'var'. Supply one of them", call. = FALSE)
  }
  if (!is.null(strata[["mean"]]) && !is.null(strata[["p"]])) {
    stop("'strata' gives both 'mean' and 'p'. Supply one of them", call. = FALSE)
  }
  if (!is.null(strata[["p"]])) {
    p <- .deff_check_column(strata$p, "p", positive = FALSE)
    if (any(p < 0 | p > 1)) {
      stop("'strata$p' must be in [0, 1]", call. = FALSE)
    }
  }
  if (!is.null(strata[["sd"]])) {
    out$sd <- .deff_check_column(strata$sd, "sd", positive = FALSE)
    if (any(out$sd < 0)) {
      stop("'strata$sd' must be non-negative", call. = FALSE)
    }
  } else if (!is.null(strata[["var"]])) {
    v <- .deff_check_column(strata$var, "var", positive = FALSE)
    if (any(v < 0)) {
      stop("'strata$var' must be non-negative", call. = FALSE)
    }
    out$sd <- sqrt(v)
  } else if (!is.null(strata[["p"]])) {
    ## The within and between parts then add to pbar * (1 - pbar) exactly,
    ## the split the two-phase frame uses.
    out$sd <- sqrt(strata$p * (1 - strata$p))
  }
  if (!is.null(strata[["mean"]])) {
    out$mean <- .deff_check_column(strata$mean, "mean", positive = FALSE)
  } else if (!is.null(strata[["p"]])) {
    out$mean <- .deff_check_column(strata$p, "p", positive = FALSE)
  }
  if (is.null(out$n) && is.null(out$sd)) {
    stop(
      "'strata' needs 'n' for the weighting component, or 'sd' and 'mean' for the stratification component",
      call. = FALSE
    )
  }
  out
}

#' Validate one stratum column
#' @keywords internal
#' @noRd
.deff_check_column <- function(value, name, positive) {
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
    stop(
      sprintf("'strata$%s' must contain finite non-missing numbers", name),
      call. = FALSE
    )
  }
  if (positive && any(value <= 0)) {
    stop(sprintf("'strata$%s' must be positive", name), call. = FALSE)
  }
  as.numeric(value)
}

#' The responding takes a cluster result's precision reads
#'
#' Each stage keeps the share of its units that respond, and the last stage
#' the ultimate-unit rate times `eligible_per_unit`, as in the precision kernel.
#' PSU response removes whole clusters, so it does not enter a take.
#' @keywords internal
#' @noRd
.deff_cluster_takes <- function(n, stages, params) {
  rr <- params[["resp_rate"]] %||% 1
  rates <- if (length(rr) == stages) {
    as.numeric(rr)
  } else {
    r <- rep(1, stages)
    r[stages] <- as.numeric(rr[1L])
    r
  }
  rates[stages] <- rates[stages] *
    (params[["eligibility_rate"]] %||% 1) *
    (params[["eligible_per_unit"]] %||% 1)
  take <- unname(as.numeric(n[2L:stages])) * rates[2L:stages]
  ## The fit clamps a take to one expected response, and n * r can land a
  ## bit below it, as .check_expected_take() allows.
  take[take < 1 & take >= 1 - 1e-9] <- 1
  take
}
