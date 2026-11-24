## The survey generics are registered when survey's namespace loads, so
## loading svyplan never loads survey.

#' @exportS3Method survey::SE
SE.svyplan_n <- function(object, ...) {
  .check_unused_dots(...)
  object$se
}

#' @exportS3Method survey::SE
SE.svyplan_prec <- function(object, ...) {
  .check_unused_dots(...)
  object$se
}

#' @exportS3Method survey::SE
SE.svyplan_cluster <- function(object, ...) {
  .check_unused_dots(...)
  object$se
}

#' @exportS3Method survey::cv
cv.svyplan_n <- function(object, ...) {
  .check_unused_dots(...)
  object$cv
}

#' @exportS3Method survey::cv
cv.svyplan_prec <- function(object, ...) {
  .check_unused_dots(...)
  object$cv
}

#' @exportS3Method survey::cv
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

#' @exportS3Method survey::SE
SE.svyplan_alloc <- function(object, ...) .svyplan_se(object, ...)

#' @exportS3Method survey::SE
SE.svyplan_twophase <- function(object, ...) .svyplan_se(object, ...)

#' @exportS3Method survey::SE
SE.svyplan_panel <- function(object, ...) .svyplan_se(object, ...)

#' @exportS3Method survey::SE
SE.svyplan_strata <- function(object, ...) .svyplan_se(object, ...)

#' @exportS3Method survey::cv
cv.svyplan_alloc <- function(object, ...) .svyplan_cv(object, ...)

#' @exportS3Method survey::cv
cv.svyplan_twophase <- function(object, ...) .svyplan_cv(object, ...)

#' @exportS3Method survey::cv
cv.svyplan_panel <- function(object, ...) .svyplan_cv(object, ...)

#' @exportS3Method survey::cv
cv.svyplan_strata <- function(object, ...) .svyplan_cv(object, ...)
