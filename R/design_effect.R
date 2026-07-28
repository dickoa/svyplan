#' Planning Design Effect
#'
#' Build the design effect you expect a complex design to produce, before
#' any data are collected, by combining the design features you are
#' planning: clustering, unequal weighting, and stratification. The result
#' is a multiplier you pass as `deff` to [n_prop()], [n_mean()],
#' [n_alloc()], or any other sizing or precision function.
#'
#' @param x A `svyplan` result to read design features from: a
#'   [n_cluster()] or [prec_cluster()] allocation, a [varcomp()] estimate,
#'   or an [n_alloc()] allocation. `NULL` (default) builds the design
#'   effect from the component arguments below.
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
#' estimate a design effect from collected data: once you have a realized
#' sample, use `survey::svymean(..., deff = TRUE)` and `survey::deff()`,
#' which compute it from the actual weights, strata, and clusters.
#'
#' ## Components
#'
#' Components are selected by which arguments you supply, and multiply
#' together:
#'
#' \deqn{DEFF = DEFF_{cluster} \times DEFF_{weight} \times DEFF_{strata}.}
#'
#' **Clustering** (`icc`, `n_per_psu`, `n_per_ssu`, `var_ratio`) uses the same
#' variance model as [n_cluster()] and [prec_cluster()], so the two always
#' agree. For a two-stage design with `n_per_psu` units taken per PSU,
#' \deqn{DEFF_{cluster} = k(1 + \delta(m - 1)),}
#' and for a three-stage design taking `n_per_psu` SSUs per PSU and
#' `n_per_ssu` units per SSU,
#' \deqn{DEFF_{cluster} = k_1\delta_1 mq + k_2(1 + \delta_2(q - 1)).}
#' Estimate `icc` (and `var_ratio`) from a previous round or pilot with
#' [varcomp()].
#'
#' In the three-stage form \eqn{k_1} rescales the components' unit variance
#' to the analysis variable and \eqn{k_2} does the same for the within-PSU
#' part, which is \eqn{1-\delta_1} of it. The two are therefore linked,
#' \deqn{k_2 = k_1(1 - \delta_1),}
#' and that identity is what makes the design effect collapse to \eqn{k_1}
#' at \eqn{m=q=1}, where one unit is taken per SSU and one SSU per PSU so no
#' clustering remains. A scalar `var_ratio` supplies \eqn{k_1} and derives
#' \eqn{k_2} from it. Supplying both explicitly overrides the identity,
#' which is meaningful only when the two ratios come from different
#' decompositions; note that `varcomp()` estimates \eqn{k_2} from its own
#' SSU-level decomposition rather than imposing the identity, so its value
#' can differ by several percent on small clusters.
#'
#' Written on components referenced to the total unit variance,
#' \eqn{\delta_1=\sigma_1^2/S^2} and \eqn{\delta_2^{tot}=\sigma_2^2/S^2},
#' the same quantity is the familiar
#' \eqn{1+\delta_1(mq-1)+\delta_2^{tot}(q-1)}. The package's `icc_ssu` is
#' referenced to the within-PSU variance instead, so
#' \eqn{\delta_2=\delta_2^{tot}/(1-\delta_1)}.
#'
#' **Unequal weighting** (`weights`, or `N` and `n` in `strata`) is Kish's
#' weighting loss. From a vector of planned weights,
#' \deqn{DEFF_{weight} = n\sum_i w_i^2 / (\sum_i w_i)^2,}
#' and from a planned stratified allocation with stratum sizes \eqn{N_h}
#' and takes \eqn{n_h} the same quantity is
#' \deqn{DEFF_{weight} = n\sum_h N_h^2/n_h / (\sum_h N_h)^2, \quad
#'       n = \sum_h n_h.}
#' Equal weights give exactly 1. This is the cost of a disproportionate
#' allocation, of weighting classes, or of an anticipated nonresponse
#' adjustment, and it is never below 1.
#'
#' **Stratification** (`N`, `sd`, and `mean` in `strata`) is the gain from
#' stratifying, the only component that can fall below 1. With
#' \eqn{W_h = N_h/N},
#' \deqn{DEFF_{strata} = \sum_h W_hS_h^2 /
#'       \left(\sum_h W_hS_h^2 + \sum_h W_h(\bar y_h - \bar y)^2\right).}
#' It measures the gain under *proportional* allocation; any departure from
#' proportional is already charged to the weighting component, so the two
#' compose without double counting.
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
#' either direction, so it is not a bound: the weighting component charges
#' a disproportionate allocation as a loss without knowing which strata
#' were favoured, so a Neyman allocation over strata that differ sharply in
#' `sd` is reported well above its true variance ratio, while a constraint
#' that forces units into a low-variance stratum is reported well below it.
#' The printed result marks a multi-component product as approximate for
#' that reason.
#'
#' Given an [n_alloc()] result there is no need to approximate at all.
#' Every quantity the ratio needs is already there, so
#' `design_effect(alloc)` returns the allocation's own variance ratio,
#' \deqn{D = d_0\,n \sum_h W_h^2S_h^2k_h(1+\delta_h(m_h-1))/n_h \Big/
#'        \left(\sum_h W_hS_h^2+\sum_h W_h(\bar y_h-\bar y)^2\right),}
#' with the per-stratum cluster factors entering stratum by stratum rather
#' than averaged, and the scalar `deff` the allocation was built under
#' counted once. A generalized allocation optimizing several measures has
#' no single such ratio and is refused; take the measure-specific numbers
#' from [prec_alloc()].
#'
#' This ratio is computed without finite population corrections, on the
#' same planning scale as the component arguments. [prec_alloc()] applies
#' each stratum's FPC, so a ratio derived from its variance need not equal
#' `design_effect(alloc)` once sampling fractions are material; use
#' [prec_alloc()] when the design variance itself is what you want.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 9.
#'
#' @seealso [effective_n()] for the sample size this design effect costs
#'   you, [varcomp()] for estimating `icc` and `var_ratio`, [n_cluster()] for
#'   optimizing a multistage design directly.
#'
#' @examples
#' # Clustering: 25 households per cluster, homogeneity 0.05
#' design_effect(icc = 0.05, n_per_psu = 25)
#'
#' # Three stages: 10 SSUs per PSU, 4 units per SSU
#' design_effect(icc = c(0.01, 0.05), n_per_psu = 10, n_per_ssu = 4)
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
#' deff <- design_effect(icc = 0.05, n_per_psu = 25, strata = frame)
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
#'   A scalar for a two-stage design, or `c(icc_psu, icc_ssu)` for a
#'   three-stage design. This is the survey-planning measure returned by
#'   [varcomp()], not a generic mixed-model ICC.
#' @param n_per_psu Units taken per PSU in a two-stage design, or SSUs
#'   taken per PSU in a three-stage design. Required with `icc`.
#' @param n_per_ssu Units taken per SSU. Supply only for a three-stage
#'   design, where `icc` has length 2.
#' @param var_ratio Ratio of the components' unit variance to the analysis
#'   variable's, default 1. A scalar names `var_ratio_psu` and, for a three-stage
#'   `icc`, derives `var_ratio_ssu = var_ratio_psu * (1 - icc_psu)`. Supply a length-2
#'   vector only to override that identity; see Details.
#' @param weights Numeric vector of planned sampling weights, one per
#'   sampled unit. Only their relative variability matters.
#' @param strata Stratum-level data frame describing a planned stratified
#'   allocation, using the [n_alloc()] column names. `N` with `n` gives
#'   the weighting component; `N` with `sd` and `mean` gives the
#'   stratification component. Supply all four for both.
#'
#' @export
design_effect.default <- function(
  x = NULL,
  ...,
  icc = NULL,
  n_per_psu = NULL,
  n_per_ssu = NULL,
  var_ratio = 1,
  weights = NULL,
  strata = NULL
) {
  .check_unused_dots(...)
  if (!is.null(x)) {
    stop(
      "'x' must be a svyplan result or NULL. Pass planning components by name (icc, n_per_psu, weights, strata)",
      call. = FALSE
    )
  }
  .deff_compose(
    icc = icc,
    n_per_psu = n_per_psu,
    n_per_ssu = n_per_ssu,
    var_ratio = var_ratio,
    weights = weights,
    strata = strata
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
#'   [n_cluster()] or [prec_cluster()] allocation.
#' @export
design_effect.svyplan_cluster <- function(x, ..., weights = NULL,
                                          strata = NULL) {
  .check_unused_dots(...)
  takes <- .deff_cluster_takes(x$n, x$stages)
  .deff_compose(
    icc = x$params$icc,
    n_per_psu = takes$n_per_psu,
    n_per_ssu = takes$n_per_ssu,
    var_ratio = x$params$var_ratio %||% 1,
    weights = weights,
    strata = strata
  )
}

#' @describeIn design_effect Read the design features off a
#'   [prec_cluster()] result, the same way as for [n_cluster()]. Other
#'   `svyplan_prec` types carry the `deff` you supplied rather than one to
#'   be derived, and are rejected.
#' @export
design_effect.svyplan_prec <- function(x, ..., weights = NULL,
                                       strata = NULL) {
  .check_unused_dots(...)
  if (!identical(x$type, "cluster")) {
    stop(
      sprintf(
        "design_effect() reads a cluster design; this is a '%s' result. Its design effect is the 'deff' you supplied",
        x$type
      ),
      call. = FALSE
    )
  }
  takes <- .deff_cluster_takes(x$params$n, x$params$stages)
  .deff_compose(
    icc = x$params$icc,
    n_per_psu = takes$n_per_psu,
    n_per_ssu = takes$n_per_ssu,
    var_ratio = x$params$var_ratio %||% 1,
    weights = weights,
    strata = strata
  )
}

#' @describeIn design_effect Read `icc` and `var_ratio` from a [varcomp()]
#'   estimate. The stage takes are your design choice, so `n_per_psu` (and
#'   `n_per_ssu` for three stages) must be supplied.
#' @export
design_effect.svyplan_varcomp <- function(x, ..., n_per_psu = NULL,
                                          n_per_ssu = NULL, weights = NULL,
                                          strata = NULL) {
  .check_unused_dots(...)
  if (!is.null(x$strata)) {
    stop(
      "stratified varcomp: pass one stratum's icc and var_ratio (see the $strata table)",
      call. = FALSE
    )
  }
  .deff_compose(
    icc = x$icc,
    n_per_psu = n_per_psu,
    n_per_ssu = n_per_ssu,
    var_ratio = x$var_ratio,
    weights = weights,
    strata = strata
  )
}

#' @describeIn design_effect Read the planned allocation from an
#'   [n_alloc()] result. This returns the allocation's own variance ratio,
#'   computed from the stratum sizes, takes, standard deviations, means and
#'   any per-stratum cluster factors, rather than a product of separate
#'   weighting and stratification approximations. Supply `weights` only for
#'   an *additional* anticipated adjustment; see the `weights` argument.
#' @param weights For an [n_alloc()] result only: an *additional*
#'   anticipated weighting adjustment charged on top of the allocation,
#'   such as a nonresponse correction or a calibration step. Do not pass
#'   the allocation's own weights. The stratum takes already carry those,
#'   and supplying them again counts the same disproportionality twice.
#'   The `allocation` component is exact for the plan; this one is Kish's
#'   approximation, so with `weights` supplied the overall result is
#'   approximate too and prints as such.
#' @export
design_effect.svyplan_n <- function(x, ..., weights = NULL) {
  .check_unused_dots(...)
  if (!identical(x$type, "alloc")) {
    stop(
      sprintf(
        "design_effect() reads a stratified allocation; this is a '%s' result. Its design effect is the 'deff' you supplied",
        x$type
      ),
      call. = FALSE
    )
  }
  if (!is.null(x$params$measures)) {
    stop(
      paste("a generalized allocation optimizes several measures at once, and",
            "they do not share one variance ratio; take the measure-specific",
            "precision from prec_alloc() instead"),
      call. = FALSE
    )
  }
  ratio <- .deff_alloc_ratio(x)
  components <- c(allocation = ratio$value)
  notes <- c(allocation = ratio$note)
  if (!is.null(weights)) {
    extra <- .deff_from_weights(weights)
    components <- c(components, weight = extra$value)
    notes <- c(notes, weight = paste0(extra$note, " (Kish)"))
  }
  .new_svyplan_deff(prod(components), components = components, notes = notes)
}

#' Variance ratio of a one-measure stratified allocation
#'
#' The design variance of the mean under this allocation, over the SRS
#' variance at the same total size. Everything it needs is already in the
#' result, so there is no reason to approximate it by multiplying a Kish
#' weighting factor by a proportional-allocation stratification factor: that
#' product is exact only when the stratum standard deviations agree, and it
#' is not even an upper bound when a constraint forces the allocation to
#' oversample a low-variance stratum.
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
  W <- detail$N / sum(detail$N)
  n <- detail$n
  factor_h <- .deff_alloc_factor(x$params$frame, detail)
  within <- sum(W * detail$sd^2)
  between <- sum(W * (detail$mean - sum(W * detail$mean))^2)
  total <- within + between
  if (total < .Machine$double.eps) {
    stop("the allocation has no variability to summarize: 'sd' and 'mean' are constant and zero",
         call. = FALSE)
  }
  # the scalar deff the allocation was built under is part of its variance
  # model, so it belongs in the ratio exactly once
  scalar <- x$params$deff %||% 1
  value <- scalar * sum(n) * sum(W^2 * detail$sd^2 * factor_h / n) / total
  list(
    value = value,
    note = sprintf("%d strata, n = %.4g (direct variance ratio)",
                   nrow(detail), sum(n))
  )
}

#' Clustering component of a stratified two-stage n_alloc() plan
#'
#' Averages the per-stratum clustering factors rather than averaging icc
#' and n_per_psu separately, which would not reproduce any stratum's
#' inflation. The weights are each stratum's share of the design variance,
#' \eqn{W_h^2S_h^2/n_h}, because that is what the factors multiply in
#' \eqn{\sum_h W_h^2S_h^2k_h(1 + \delta_h(m_h-1))/n_h}. Weighting by the
#' allocation instead would answer a different question and does not
#' reproduce the variance ratio even when every stratum is alike in `sd`.
#' @keywords internal
#' @noRd
.deff_alloc_cluster <- function(frame, detail) {
  if (!is.data.frame(frame) || is.null(frame[["icc_psu"]]) ||
      is.null(detail[["n_per_psu"]]) || nrow(frame) != nrow(detail)) {
    return(NULL)
  }
  icc <- as.numeric(frame$icc_psu)
  var_ratio <- if (is.null(frame[["var_ratio_psu"]])) rep(1, nrow(frame)) else
    as.numeric(frame$var_ratio_psu)
  factor_h <- var_ratio * (1 + icc * (detail$n_per_psu - 1))
  share <- .deff_variance_share(detail)
  list(
    value = sum(share * factor_h),
    note = sprintf(
      "%d strata, icc in [%.4g, %.4g]", length(icc),
      min(icc), max(icc)
    )
  )
}

#' Per-stratum clustering factor of an allocation, or one where there is none
#'
#' The factors enter the variance ratio stratum by stratum, so unlike the
#' component API there is nothing to average here.
#' @keywords internal
#' @noRd
.deff_alloc_factor <- function(frame, detail) {
  if (!is.data.frame(frame) || is.null(frame[["icc_psu"]]) ||
      is.null(detail[["n_per_psu"]]) || nrow(frame) != nrow(detail)) {
    return(rep(1, nrow(detail)))
  }
  var_ratio <- if (is.null(frame[["var_ratio_psu"]])) rep(1, nrow(frame)) else
    as.numeric(frame$var_ratio_psu)
  var_ratio * (1 + as.numeric(frame$icc_psu) * (detail$n_per_psu - 1))
}

#' Each stratum's share of the design variance of the overall mean
#'
#' Falls back to the allocation when the frame does not carry the stratum
#' sizes and variabilities the exact weights need.
#' @keywords internal
#' @noRd
.deff_variance_share <- function(detail) {
  n <- detail[["n"]]
  N <- detail[["N"]]
  sd <- detail[["sd"]]
  share <- if (!is.null(N) && !is.null(sd) && !is.null(n)) {
    ifelse(n > 0, N^2 * sd^2 / n, 0)
  } else {
    n
  }
  total <- sum(share)
  if (!is.finite(total) || total <= 0) {
    share <- n
    total <- sum(share)
  }
  share / total
}

#' Combine the requested planning components into a svyplan_deff
#' @keywords internal
#' @noRd
.deff_compose <- function(icc, n_per_psu, n_per_ssu, var_ratio, weights, strata) {
  parts <- numeric(0)
  notes <- character(0)

  if (!is.null(icc) || !is.null(n_per_psu) || !is.null(n_per_ssu)) {
    cl <- .deff_from_cluster(icc, n_per_psu, n_per_ssu, var_ratio)
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

  if (length(parts) == 0L) {
    stop(
      "no design components supplied. Give at least one of 'icc' with 'n_per_psu', 'weights', or 'strata'",
      call. = FALSE
    )
  }

  .new_svyplan_deff(prod(parts), components = parts, notes = notes)
}

#' Clustering component, matching the prec_cluster() variance model
#' @keywords internal
#' @noRd
.deff_from_cluster <- function(icc, n_per_psu, n_per_ssu, var_ratio) {
  if (is.null(icc) || is.null(n_per_psu)) {
    stop(
      "the clustering component needs both 'icc' and 'n_per_psu'",
      call. = FALSE
    )
  }
  icc <- .reorder_stage_vec(icc, "icc")
  var_ratio <- .reorder_stage_vec(var_ratio, "var_ratio")
  check_icc(icc)
  if (length(icc) > 2L) {
    stop(
      "'icc' must have length 1 (two stages) or 2 (three stages)",
      call. = FALSE
    )
  }
  three <- length(icc) == 2L
  if (three && is.null(n_per_ssu)) {
    stop(
      "a three-stage 'icc' also needs 'n_per_ssu' (units taken per SSU)",
      call. = FALSE
    )
  }
  if (!three && !is.null(n_per_ssu)) {
    stop(
      "'n_per_ssu' applies to three stages only; 'icc' must then have length 2",
      call. = FALSE
    )
  }
  if (!is.numeric(var_ratio) || !length(var_ratio) %in% c(1L, length(icc)) || anyNA(var_ratio) ||
      any(!is.finite(var_ratio)) || any(var_ratio <= 0)) {
    stop(
      sprintf(
        "'var_ratio' must contain 1 or %d positive finite value(s)", length(icc)
      ),
      call. = FALSE
    )
  }
  var_ratio <- if (three) .stage_k_pair(var_ratio, icc) else rep_len(var_ratio, length(icc))
  .deff_check_take(n_per_psu, "n_per_psu")
  if (three) .deff_check_take(n_per_ssu, "n_per_ssu")

  if (three) {
    value <- var_ratio[1L] * icc[1L] * n_per_psu * n_per_ssu +
      var_ratio[2L] * (1 + icc[2L] * (n_per_ssu - 1))
    note <- sprintf(
      "icc = (%.4g, %.4g), n_per_psu = %.4g, n_per_ssu = %.4g",
      icc[1L], icc[2L], n_per_psu, n_per_ssu
    )
  } else {
    value <- var_ratio * (1 + icc * (n_per_psu - 1))
    note <- sprintf("icc = %.4g, n_per_psu = %.4g", icc, n_per_psu)
  }
  list(value = as.numeric(value), note = note)
}

#' Weighting component from a vector of planned weights
#' @keywords internal
#' @noRd
.deff_from_weights <- function(w) {
  check_weights(w, "weights")
  n <- length(w)
  value <- n * sum(w^2) / sum(w)^2
  list(value = value, note = sprintf("cv(w) = %.4g", sqrt(max(value - 1, 0))))
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
  within <- sum(share * sd^2)
  between <- sum(share * (mean - sum(share * mean))^2)
  total <- within + between
  if (total < .Machine$double.eps) {
    stop(
      "'strata' has no variability: 'sd' and 'mean' are constant and zero",
      call. = FALSE
    )
  }
  list(
    value = within / total,
    note = sprintf("%d strata, between-stratum share %.4g", length(N),
                   between / total)
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

#' Extract per-stage takes from a svyplan_cluster $n vector
#' @keywords internal
#' @noRd
.deff_cluster_takes <- function(n, stages) {
  if (stages == 2L) {
    list(n_per_psu = unname(n[["n_per_psu"]]), n_per_ssu = NULL)
  } else {
    list(n_per_psu = unname(n[["n_per_psu"]]), n_per_ssu = unname(n[["n_per_ssu"]]))
  }
}
