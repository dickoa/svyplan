#' Construct a svyplan_n object
#' @keywords internal
#' @noRd
.new_svyplan_n <- function(n, type, method = NULL, params = list(),
                           targets = NULL, detail = NULL,
                           binding = NULL, domains = NULL) {
  structure(
    list(
      n       = n,
      type    = type,
      method  = method,
      params  = params,
      targets = targets,
      detail  = detail,
      binding = binding,
      domains = domains
    ),
    class = c("svyplan_n", "list")
  )
}

#' Construct a svyplan_cluster object
#' @keywords internal
#' @noRd
.new_svyplan_cluster <- function(n, stages, total_n, cv, cost, params = list(),
                                 targets = NULL, detail = NULL,
                                 binding = NULL, domains = NULL) {
  structure(
    list(
      n       = n,
      stages  = stages,
      total_n = total_n,
      cv      = cv,
      cost    = cost,
      params  = params,
      targets = targets,
      detail  = detail,
      binding = binding,
      domains = domains
    ),
    class = c("svyplan_cluster", "list")
  )
}

#' Construct a svyplan_varcomp object
#' @keywords internal
#' @noRd
.new_svyplan_varcomp <- function(var_between, var_within, delta, k,
                                 rel_var, stages) {
  structure(
    list(
      var_between = var_between,
      var_within  = var_within,
      delta       = delta,
      k           = k,
      rel_var     = rel_var,
      stages      = stages
    ),
    class = c("svyplan_varcomp", "list")
  )
}

#' Construct a svyplan_power object
#' @keywords internal
#' @noRd
.new_svyplan_power <- function(n, power, delta, type, solved, params = list()) {
  structure(
    list(
      n      = n,
      power  = power,
      delta  = delta,
      type   = type,
      solved = solved,
      params = params
    ),
    class = c("svyplan_power", "list")
  )
}

#' Construct a svyplan_strata object
#' @keywords internal
#' @noRd
.new_svyplan_strata <- function(boundaries, n_strata, n, cv, strata,
                                method, alloc, params,
                                converged = NA) {
  structure(
    list(
      boundaries = boundaries,
      n_strata   = n_strata,
      n          = n,
      cv         = cv,
      strata     = strata,
      method     = method,
      alloc      = alloc,
      params     = params,
      converged  = converged
    ),
    class = c("svyplan_strata", "list")
  )
}
