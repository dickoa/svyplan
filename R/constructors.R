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
  targets = NULL,
  detail = NULL,
  binding = NULL,
  domains = NULL,
  operational = NULL
) {
  stopifnot(
    is.numeric(n),
    is.character(type),
    length(type) == 1L,
    !is.na(type),
    is.list(params),
    is.null(detail) || is.data.frame(detail),
    is.null(domains) || is.data.frame(domains),
    is.null(operational) || is.list(operational)
  )

  if (is.null(se) || is.null(moe) || is.null(cv)) {
    prec <- .compute_prec_for_n(n, type, params, method = method)
    se <- se %||% prec$se
    moe <- moe %||% prec$moe
    cv <- cv %||% prec$cv
  }
  stopifnot(is.numeric(se), is.numeric(moe), is.numeric(cv))

  structure(
    list(
      n = n,
      type = type,
      method = method,
      params = params,
      se = se,
      moe = moe,
      cv = cv,
      targets = targets,
      detail = detail,
      binding = binding,
      domains = domains,
      operational = operational
    ),
    class = c("svyplan_n", "list")
  )
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
                      method %||% "wald")
  } else if (type == "mean") {
    .prec_engine_mean(params$var, params$mu, n, params$alpha, N, deff,
                      resp_rate)
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
  targets = NULL,
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
      targets = targets,
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
  detail = NULL
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
    is.null(detail) || is.data.frame(detail)
  )

  structure(
    list(
      se = se,
      moe = moe,
      cv = cv,
      type = type,
      method = method,
      params = params,
      detail = detail
    ),
    class = c("svyplan_prec", "list")
  )
}

#' Construct a svyplan_varcomp object
#' @keywords internal
#' @noRd
.new_svyplan_varcomp <- function(
  varb,
  varw,
  delta,
  k,
  rel_var,
  stages,
  strata = NULL
) {
  stopifnot(
    is.numeric(stages),
    length(stages) == 1L,
    is.null(strata) || is.data.frame(strata)
  )
  if (is.null(strata)) {
    stopifnot(
      is.numeric(varb),
      is.numeric(varw),
      is.numeric(delta),
      is.numeric(k),
      is.numeric(rel_var)
    )
  } else {
    stopifnot(
      is.null(varb),
      is.null(varw),
      is.null(delta),
      is.null(k),
      is.null(rel_var)
    )
  }

  structure(
    list(
      varb = varb,
      varw = varw,
      delta = delta,
      k = k,
      rel_var = rel_var,
      stages = stages,
      strata = strata
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

#' Construct a svyplan_design_effect object
#' @keywords internal
#' @noRd
.new_svyplan_design_effect <- function(
  design_effect,
  method,
  components = NULL
) {
  stopifnot(
    is.numeric(design_effect),
    length(design_effect) == 1L,
    !is.na(design_effect),
    is.finite(design_effect),
    design_effect >= 0,
    is.character(method),
    length(method) == 1L,
    !is.na(method),
    is.null(components) || is.data.frame(components)
  )

  structure(
    as.double(design_effect),
    method = method,
    components = components,
    class = c("svyplan_design_effect", "numeric")
  )
}
