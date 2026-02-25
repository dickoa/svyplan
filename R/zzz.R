.onLoad <- function(libname, pkgname) {
  if (requireNamespace("survey", quietly = TRUE)) {
    ns <- asNamespace("survey")
    registerS3method("SE", "svyplan_n", SE.svyplan_n, envir = ns)
    registerS3method("SE", "svyplan_prec", SE.svyplan_prec, envir = ns)
    registerS3method("SE", "svyplan_cluster", SE.svyplan_cluster, envir = ns)
    registerS3method("cv", "svyplan_n", cv.svyplan_n, envir = ns)
    registerS3method("cv", "svyplan_prec", cv.svyplan_prec, envir = ns)
    registerS3method("cv", "svyplan_cluster", cv.svyplan_cluster, envir = ns)
  }
}

#' @keywords internal
#' @noRd
SE.svyplan_n <- function(object, ...) {
  object$se
}

#' @keywords internal
#' @noRd
SE.svyplan_prec <- function(object, ...) {
  object$se
}

#' @keywords internal
#' @noRd
SE.svyplan_cluster <- function(object, ...) {
  object$se
}

#' @keywords internal
#' @noRd
cv.svyplan_n <- function(object, ...) {
  object$cv
}

#' @keywords internal
#' @noRd
cv.svyplan_prec <- function(object, ...) {
  object$cv
}

#' @keywords internal
#' @noRd
cv.svyplan_cluster <- function(object, ...) {
  object$cv
}
