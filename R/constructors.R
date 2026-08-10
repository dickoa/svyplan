#' Construct a svyplan_n object
#' @keywords internal
#' @noRd
.new_svyplan_n <- function(
  n,
  type,
  method = NULL,
  params = list(),
  se = NULL,
  moe = NULL,
  cv = NULL,
  indicators = NULL,
  targets = NULL,
  detail = NULL,
  binding = NULL,
  domains = NULL,
  operational = NULL,
  constraints = NULL,
  optimization = NULL,
  objective = NULL,
  objective_value = NULL
) {
  stopifnot(
    is.numeric(n),
    is.character(type),
    length(type) == 1L,
    !is.na(type),
    is.list(params),
    is.null(detail) || is.data.frame(detail),
    is.null(domains) || is.data.frame(domains),
    is.null(operational) || is.list(operational),
    is.null(constraints) || is.data.frame(constraints),
    is.null(optimization) || is.list(optimization),
    is.null(objective) || is.data.frame(objective),
    is.null(objective_value) ||
      (is.numeric(objective_value) && length(objective_value) == 1L)
  )

  if (is.null(se) || is.null(moe) || is.null(cv)) {
    prec <- .compute_prec_for_n(n, type, params, method = method)
    se <- se %||% prec$se
    moe <- moe %||% prec$moe
    cv <- cv %||% prec$cv
  }
  stopifnot(is.numeric(se), is.numeric(moe), is.numeric(cv))

  out <- list(
    n = n,
    type = type,
    method = method,
    params = params,
    se = se,
    moe = moe,
    cv = cv,
    rmoe = .rmoe_for_result(moe, type, params),
    indicators = indicators,
    targets = targets,
    detail = detail,
    binding = binding,
    domains = domains,
    operational = operational
  )
  cases <- .expected_cases(n, params, type)
  if (!is.null(cases)) out$expected_cases <- cases
  # Preserve the exact legacy schema unless the optional generalized
  # allocation diagnostics are actually present.
  if (!is.null(constraints)) out$constraints <- constraints
  if (!is.null(optimization)) out$optimization <- optimization
  if (!is.null(objective)) out$objective <- objective
  if (!is.null(objective_value)) out$objective_value <- objective_value

  structure(out, class = c("svyplan_n", "list"))
}

#' The achieved margin of error on the estimand's own scale
#'
#' Reported wherever `moe` is, so a target stated as `rmoe` reads back in
#' the units it was stated in. The scale is the estimand the interval
#' bounds, `p` for a proportion and `mu` for a mean; a result whose
#' estimand has no known scale reports `NA`, which is what `$cv` already
#' does in an allocation over a frame carrying neither.
#' @keywords internal
#' @noRd
.rmoe_for_result <- function(moe, type, params) {
  estimand <- switch(
    type,
    proportion = params$p,
    mean = params$mu,
    change = params$change,
    pooled = params$mu,
    multi = .indicator_estimand(params$indicators, length(moe)),
    alloc = .alloc_estimand(params$frame),
    NULL
  )
  .rmoe_from_moe(moe, estimand)
}

#' The population mean an allocation's overall margin of error refers to
#'
#' The size-weighted mean of the stratum means, which is the estimand the
#' top-level `se` and `moe` describe. `NULL` when the frame carries neither
#' `mean` nor `p`, the same state that leaves `$cv` `NA`.
#' @keywords internal
#' @noRd
.alloc_estimand <- function(frame) {
  if (!is.data.frame(frame) || !"N" %in% names(frame)) {
    return(NULL)
  }
  mean_h <- if ("mean" %in% names(frame)) {
    frame$mean
  } else if ("p" %in% names(frame)) {
    frame$p
  } else {
    return(NULL)
  }
  if (!is.numeric(mean_h) || anyNA(mean_h)) {
    return(NULL)
  }
  .aggregate_mean(frame$N / sum(frame$N), mean_h)
}

#' The per-row estimand of an indicator table, or NULL when it has none
#'
#' A table may hold proportion rows, mean rows, or both, so the scale is
#' taken column by column rather than from the table's type.
#' @keywords internal
#' @noRd
.indicator_estimand <- function(indicators, n_rows) {
  if (!is.data.frame(indicators) || nrow(indicators) != n_rows) {
    return(NULL)
  }
  out <- rep(NA_real_, n_rows)
  if ("p" %in% names(indicators)) out <- as.numeric(indicators$p)
  if ("mu" %in% names(indicators)) {
    out <- ifelse(is.na(out), as.numeric(indicators$mu), out)
  }
  out
}

#' Compute precision measures from n and params
#' @keywords internal
#' @noRd
.compute_prec_for_n <- function(n, type, params, method = NULL) {
  if (type == "multi") {
    return(list(se = NA_real_, moe = NA_real_, cv = NA_real_))
  }

  deff <- if (!is.null(params$deff)) params$deff else 1
  resp_rate <- if (!is.null(params$resp_rate)) params$resp_rate else 1
  N <- if (!is.null(params$N)) params$N else Inf

  if (type == "proportion") {
    .prec_engine_prop(params$p, n, params$alpha, N, deff, resp_rate,
                      method %||% "wald", params$df)
  } else if (type == "mean") {
    .prec_engine_mean(params$var, params$mu, n, params$alpha, N, deff,
                      resp_rate, params$df)
  } else {
    list(se = NA_real_, moe = NA_real_, cv = NA_real_)
  }
}

#' Construct a svyplan_cluster object
#' @keywords internal
#' @noRd
.new_svyplan_cluster <- function(
  n,
  stages,
  total_n,
  cv,
  cost,
  params = list(),
  indicators = NULL,
  detail = NULL,
  binding = NULL,
  domains = NULL,
  operational = NULL
) {
  stopifnot(
    is.numeric(n),
    is.numeric(stages),
    length(stages) == 1L,
    is.numeric(total_n),
    length(total_n) == 1L,
    is.numeric(cv),
    length(cv) == 1L,
    is.numeric(cost),
    length(cost) == 1L,
    is.list(params),
    is.null(detail) || is.data.frame(detail),
    is.null(domains) || is.data.frame(domains),
    is.null(operational) || is.list(operational)
  )

  .warn_single_psu((operational$n %||% n)[[1L]])

  structure(
    list(
      n = n,
      stages = stages,
      total_n = total_n,
      se = NA_real_,
      moe = NA_real_,
      cv = cv,
      cost = cost,
      operational = operational,
      params = params,
      indicators = indicators,
      detail = detail,
      binding = binding,
      domains = domains
    ),
    class = c("svyplan_cluster", "list")
  )
}

#' Construct a svyplan_prec object
#' @keywords internal
#' @noRd
.new_svyplan_prec <- function(
  se,
  moe,
  cv,
  type,
  method = NULL,
  params = list(),
  detail = NULL,
  bounds = NULL,
  domains = NULL,
  solved = NULL,
  objective = NULL,
  objective_value = NULL
) {
  stopifnot(
    is.numeric(se),
    is.numeric(moe),
    is.numeric(cv),
    length(se) == length(moe),
    length(se) == length(cv),
    is.character(type),
    length(type) == 1L,
    !is.na(type),
    is.list(params),
    is.null(detail) || is.data.frame(detail),
    is.null(bounds) || is.data.frame(bounds),
    is.null(domains) || is.data.frame(domains),
    is.null(objective) || is.data.frame(objective),
    is.null(objective_value) ||
      (is.numeric(objective_value) && length(objective_value) == 1L)
  )

  out <- list(
    se = se,
    moe = moe,
    cv = cv,
    rmoe = .rmoe_for_result(moe, type, params),
    type = type,
    method = method,
    params = params,
    detail = detail
  )
  cases <- .expected_cases(params$n, params, type)
  if (!is.null(cases)) out$expected_cases <- cases
  if (!is.null(bounds)) out$bounds <- bounds
  if (!is.null(domains)) out$domains <- domains
  # Recorded only when the level was solved for, so a result computed in
  # the usual direction keeps the schema it has always had.
  if (!is.null(solved)) out$solved <- solved
  if (!is.null(objective)) out$objective <- objective
  if (!is.null(objective_value)) out$objective_value <- objective_value

  structure(out, class = c("svyplan_prec", "list"))
}

#' Construct a svyplan_varcomp object
#' @keywords internal
#' @noRd
.new_svyplan_varcomp <- function(
  varb,
  varw,
  icc,
  var_ratio,
  unit_relvar,
  stages,
  strata = NULL,
  source = "data",
  params = list()
) {
  stopifnot(
    is.numeric(stages),
    length(stages) == 1L,
    is.null(strata) || is.data.frame(strata),
    is.character(source),
    length(source) == 1L,
    source %in% c("data", "deff"),
    is.list(params)
  )
  if (is.null(strata)) {
    stopifnot(
      is.numeric(varb),
      is.numeric(varw),
      is.numeric(icc),
      is.numeric(var_ratio),
      is.numeric(unit_relvar)
    )
  } else {
    stopifnot(
      is.null(varb),
      is.null(varw),
      is.null(icc),
      is.null(var_ratio),
      is.null(unit_relvar)
    )
  }

  structure(
    list(
      varb = varb,
      varw = varw,
      icc = icc,
      var_ratio = var_ratio,
      unit_relvar = unit_relvar,
      stages = stages,
      strata = strata,
      source = source,
      params = params
    ),
    class = c("svyplan_varcomp", "list")
  )
}

#' Construct a svyplan_power object
#' @keywords internal
#' @noRd
.new_svyplan_power <- function(n, power, effect, type, solved, params = list()) {
  stopifnot(
    is.numeric(n),
    length(n) %in% 1:2,
    is.numeric(power),
    length(power) == 1L,
    is.numeric(effect),
    length(effect) == 1L,
    is.character(type),
    length(type) == 1L,
    is.character(solved),
    length(solved) == 1L,
    is.list(params)
  )

  structure(
    list(
      n = n,
      power = power,
      effect = effect,
      type = type,
      solved = solved,
      params = params
    ),
    class = c("svyplan_power", "list")
  )
}

#' Construct a svyplan_strata object
#' @keywords internal
#' @noRd
.new_svyplan_strata <- function(
  boundaries,
  n_strata,
  n,
  cv,
  strata,
  method,
  alloc,
  params,
  converged = NA
) {
  stopifnot(
    is.numeric(boundaries),
    is.numeric(n_strata),
    length(n_strata) == 1L,
    is.numeric(n),
    length(n) == 1L,
    is.numeric(cv),
    length(cv) == 1L,
    is.data.frame(strata),
    is.character(method),
    length(method) == 1L,
    is.character(alloc),
    length(alloc) == 1L,
    is.list(params)
  )

  structure(
    list(
      boundaries = boundaries,
      n_strata = n_strata,
      n = n,
      cv = cv,
      strata = strata,
      method = method,
      alloc = alloc,
      params = params,
      converged = converged
    ),
    class = c("svyplan_strata", "list")
  )
}

#' Construct a svyplan_deff object
#'
#' `components` is a named numeric vector of the multiplicative parts and
#' `notes` a parallel character vector describing each part's inputs.
#' @keywords internal
#' @noRd
.new_svyplan_deff <- function(deff, components, notes) {
  stopifnot(
    is.numeric(deff),
    length(deff) == 1L,
    !is.na(deff),
    is.finite(deff),
    deff > 0,
    is.numeric(components),
    length(components) > 0L,
    !is.null(names(components)),
    all(is.finite(components)),
    is.character(notes),
    identical(names(notes), names(components))
  )

  structure(
    as.double(deff),
    components = components,
    notes = notes,
    class = c("svyplan_deff", "numeric")
  )
}

#' Construct a svyplan_df object
#'
#' A classed numeric scalar, following `svyplan_deff`, so a planned design's
#' degrees of freedom stay usable as a plain number wherever one is
#' expected while still carrying the per-stratum and per-domain detail they
#' were counted from.
#' @keywords internal
#' @noRd
.new_svyplan_df <- function(df, n_units, n_strata, stage, strata = NULL,
                            domains = NULL) {
  stopifnot(
    is.numeric(df),
    length(df) == 1L,
    !is.na(df),
    is.finite(df),
    is.numeric(n_units),
    length(n_units) == 1L,
    is.numeric(n_strata),
    length(n_strata) == 1L,
    is.character(stage),
    length(stage) == 1L,
    stage %in% c("psu", "element"),
    is.null(strata) || is.data.frame(strata),
    is.null(domains) || is.data.frame(domains)
  )

  structure(
    as.double(df),
    n_units = as.double(n_units),
    n_strata = as.integer(n_strata),
    stage = stage,
    strata = strata,
    domains = domains,
    class = c("svyplan_df", "numeric")
  )
}

#' Construct a svyplan_panel object
#'
#' A sibling of `svyplan_n` rather than a subtype: the headline number is a
#' recruitment count, and the methods registered for an analysis sample
#' would each answer a different question about it. The recruitment is
#' stored under a design-specific name, `n_issued` for the one cohort a
#' fixed panel releases and `n_entrants` for what a rotating design takes on
#' each occasion, because those are different quantities and one name over
#' both is how the wrong one gets read.
#' @keywords internal
#' @noRd
.new_svyplan_panel <- function(design, n_recruit, n_target, n_resp, n_assured,
                               assured_feasible, target_wave, waves, prec,
                               target, solved, params, start = NULL,
                               launch = NULL, launch_waves = NULL) {
  stopifnot(
    is.character(design),
    length(design) == 1L,
    design %in% c("fixed", "rotating"),
    is.numeric(n_recruit),
    length(n_recruit) == 1L,
    is.numeric(n_target),
    length(n_target) == 1L,
    is.numeric(n_resp),
    length(n_resp) == 1L,
    is.null(n_assured) || (is.numeric(n_assured) && length(n_assured) == 1L),
    is.null(assured_feasible) ||
      (is.logical(assured_feasible) && length(assured_feasible) == 1L),
    is.data.frame(waves),
    is.list(prec),
    is.list(params),
    is.null(start) ||
      (is.character(start) && length(start) == 1L &&
         start %in% c("gradual", "immediate")),
    is.null(launch) || is.data.frame(launch),
    is.null(launch_waves) || is.data.frame(launch_waves),
    # a launch is a property of a rotating design, and the three fields are
    # one fact, so they arrive together or not at all
    is.null(start) == is.null(launch),
    is.null(start) == is.null(launch_waves),
    is.null(start) || !identical(design, "fixed")
  )
  k <- nrow(waves)

  out <- list(design = design)
  if (!is.null(start)) out$start <- start
  if (identical(design, "fixed")) {
    out$n_issued <- n_recruit
    out$target_wave <- target_wave
  } else {
    out$n_entrants <- n_recruit
    out$n_in_sample <- k * n_recruit
    out$n_cohorts <- k
  }
  out$n_target <- n_target
  out$n_resp <- n_resp
  if (!is.null(n_assured)) {
    out$n_assured <- n_assured
    out$assured_feasible <- isTRUE(assured_feasible)
  }
  out$se <- prec$se
  out$moe <- prec$moe
  out$cv <- prec$cv
  out$rmoe <- prec$rmoe
  out$type <- target$type
  out$method <- target$method
  out$waves <- waves
  if (!is.null(launch)) {
    out$launch <- launch
    out$launch_waves <- launch_waves
  }
  out$target <- target
  # Recorded only in the forward direction, matching svyplan_prec, so a
  # result read from a recruitment carries no claim to have solved for it.
  if (!is.null(solved)) out$solved <- solved
  out$params <- params

  structure(out, class = c("svyplan_panel", "list"))
}

#' Construct a svyplan_twophase object
#' @keywords internal
#' @noRd
.new_svyplan_twophase <- function(n, responding, cv, cost, detail, single_phase,
                                  operational = NULL, params = list()) {
  stopifnot(
    is.numeric(n), length(n) == 2L,
    is.numeric(responding), length(responding) == 2L,
    is.numeric(cost), length(cost) == 1L,
    is.data.frame(detail),
    is.list(single_phase)
  )
  structure(
    list(
      n = n,
      responding = responding,
      cv = cv,
      cost = cost,
      detail = detail,
      single_phase = single_phase,
      operational = operational,
      params = params
    ),
    class = "svyplan_twophase"
  )
}
