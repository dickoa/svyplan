.onLoad <- function(libname, pkgname) {
  if (requireNamespace("survey", quietly = TRUE)) {
    ns <- asNamespace("survey")
    registerS3method("SE", "svyplan_n", SE.svyplan_n, envir = ns)
    registerS3method("SE", "svyplan_prec", SE.svyplan_prec, envir = ns)
    registerS3method("SE", "svyplan_cluster", SE.svyplan_cluster, envir = ns)
    registerS3method("cv", "svyplan_n", cv.svyplan_n, envir = ns)
    registerS3method("cv", "svyplan_prec", cv.svyplan_prec, envir = ns)
    registerS3method("cv", "svyplan_cluster", cv.svyplan_cluster, envir = ns)
    for (cls in c("svyplan_twophase", "svyplan_panel", "svyplan_strata")) {
      registerS3method("SE", cls, .svyplan_se, envir = ns)
      registerS3method("cv", cls, .svyplan_cv, envir = ns)
    }
  }
}

#' @keywords internal
#' @noRd
SE.svyplan_n <- function(object, ...) {
  .check_unused_dots(...)
  object$se
}

#' @keywords internal
#' @noRd
SE.svyplan_prec <- function(object, ...) {
  .check_unused_dots(...)
  object$se
}

#' @keywords internal
#' @noRd
SE.svyplan_cluster <- function(object, ...) {
  .check_unused_dots(...)
  object$se
}

#' @keywords internal
#' @noRd
cv.svyplan_n <- function(object, ...) {
  .check_unused_dots(...)
  object$cv
}

#' @keywords internal
#' @noRd
cv.svyplan_prec <- function(object, ...) {
  .check_unused_dots(...)
  object$cv
}

#' @keywords internal
#' @noRd
cv.svyplan_cluster <- function(object, ...) {
  .check_unused_dots(...)
  object$cv
}

#' `SE()` and `cv()` for results that may not store a standard error
#'
#' Read with `[[` so that a class without `se` gives `NA` rather than a
#' partial match on another field (`$se` on a stratification result finds
#' `search`).
#' @keywords internal
#' @noRd
.svyplan_se <- function(object, ...) {
  .check_unused_dots(...)
  object[["se"]] %||% NA_real_
}

#' @keywords internal
#' @noRd
.svyplan_cv <- function(object, ...) {
  .check_unused_dots(...)
  object[["cv"]] %||% NA_real_
}
