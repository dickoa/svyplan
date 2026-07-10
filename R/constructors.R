#' Construct a svyplan_n object
#' @keywords internal
#' @noRd
.new_svyplan_n <- function(
  n,
  type,
  method = NULL,
  params = list(),
  targets = NULL,
  detail = NULL,
  binding = NULL,
  domains = NULL
) {
  prec <- .compute_prec_for_n(n, type, params, method = method)

  structure(
    list(
      n = n,
      type = type,
      method = method,
      params = params,
      se = prec$se,
      moe = prec$moe,
      cv = prec$cv,
      targets = targets,
      detail = detail,
      binding = binding,
      domains = domains
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
  domains = NULL
) {
  structure(
    list(
      n = n,
      stages = stages,
      total_n = total_n,
      se = NA_real_,
      moe = NA_real_,
      cv = cv,
      cost = cost,
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
