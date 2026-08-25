#' Precision a panel recruitment delivers, wave by wave
#'
#' Take a recruitment size a panel already has, or has been budgeted, and
#' report the responding sample and the precision it leaves at each wave of
#' a unit's life. This is the inverse of [n_panel()], running the same rates
#' and the same embedded estimand from a size instead of solving for one.
#'
#' @param n_recruit For the default method: units recruited, meaning the
#'   whole issue to one cohort for a fixed panel and the entrants per
#'   occasion for a rotating one. For `svyplan_panel` objects: a result from
#'   [n_panel()] or from this function.
#' @param target A `svyplan_n` or `svyplan_prec` result for a mean or a
#'   proportion. Required in the default method, and the same object
#'   [n_panel()] takes. It supplies the estimand each wave's precision is
#'   computed for, and the responding sample the design is compared against.
#' @param retention Conditional retention, one value per wave transition,
#'   each in (0, 1\]. See [n_panel()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected. In the `svyplan_panel` method, named arguments override the
#'   stored ones, so a stored plan can be re-read under worse retention.
#' @param resp_rate Response rate at recruitment, wave 1, in (0, 1\].
#'   Default 1.
#' @param design `"fixed"` or `"rotating"`. See [n_panel()].
#' @param target_wave The wave the headline precision is reported at,
#'   defaulting to the last. It must be absent for a rotating design, whose
#'   precision belongs to the pooled occasion.
#' @param start `"gradual"`, `"immediate"` or `NULL`. See [n_panel()]. On a
#'   stored plan this may be overridden, unlike `design`. It selects which
#'   launch is described and moves no stored quantity, the recruitment being
#'   fixed before any of it is computed.
#' @param assurance Probability in (0, 1), or `NULL` (default). Reports the
#'   recruitment the target would need at that level, which is what the
#'   supplied `n_recruit` can then be read against.
#'
#' @return A `svyplan_panel` object, the class [n_panel()] returns, with
#'   `$solved` absent because nothing was solved for. `$n_resp` and the
#'   headline `se`, `moe` and `cv` describe the supplied recruitment, and
#'   `$n_target` the requirement it is being compared against, so a
#'   recruitment below the requirement reports a `moe` above the target's.
#'
#' @details
#' The precision at a wave is the embedded estimand evaluated at that wave's
#' expected respondents with `resp_rate = 1`, the panel's own losses having
#' already been applied. The design effect, population size, interval method
#' and degrees of freedom all come from the target, which is what makes the
#' round trip exact:
#' `prec_panel(n_panel(target, retention, target_wave = w))` reproduces the
#' target's own precision at wave `w`.
#'
#' A rotating panel has one precision rather than one per wave, its estimate
#' pooling every cohort alive at the occasion. The per-wave rows of `$waves`
#' still report what an estimate from a single cohort would carry, which is
#' the wave-1-only estimate a rotating design sometimes publishes, and they
#' are not the occasion's precision.
#'
#' @family precision functions
#' @seealso [n_panel()] for the inverse (solve the recruitment from a
#'   target), [design_schedule()] for turning a rotating plan into an
#'   operational schedule, [prec_mean()] and [prec_prop()] for the
#'   single-occasion precision the waves are evaluated with.
#'
#' @examples
#' target <- n_prop(p = 0.5, moe = 0.031)
#' ret <- c(0.878, 0.963, 0.936, 0.956)
#'
#' # A budget of 1500 addresses rather than the 1816 the target asks for
#' short <- prec_panel(1500, target, retention = ret, resp_rate = 0.728)
#' short
#'
#' # The round trip: precision at the wave the panel was sized for
#' plan <- n_panel(target, retention = ret, resp_rate = 0.728)
#' prec_panel(plan)$moe
#' target$moe
#'
#' # Re-read a stored plan under retention that turned out worse
#' prec_panel(plan, retention = c(0.80, 0.90, 0.90, 0.92))$waves
#'
#' @export
prec_panel <- function(n_recruit, ...) {
  UseMethod("prec_panel")
}

#' @rdname prec_panel
#' @export
prec_panel.default <- function(
  n_recruit,
  target,
  retention,
  ...,
  resp_rate = 1,
  design = c("fixed", "rotating"),
  target_wave = NULL,
  assurance = NULL,
  start = NULL
) {
  .check_unused_dots(...)
  design <- .check_panel_design(design)
  tgt <- .panel_target(target)
  q <- .panel_q(resp_rate, retention)
  target_wave <- .check_target_wave(target_wave, length(q), design)
  assurance <- .check_assurance(assurance)
  start <- .check_panel_start(start, design)

  .panel_result(
    tgt, n_recruit, resp_rate, retention, q, design, target_wave, assurance,
    start
  )
}

#' @rdname prec_panel
#' @export
prec_panel.svyplan_panel <- function(n_recruit, ...) {
  x <- n_recruit
  dots <- list(...)
  p <- x$params
  # The stored count means different units under each design, so overriding
  # 'design' would reinterpret it rather than re-assume it.
  if ("design" %in% names(dots)) {
    stop(
      sprintf(
        "'design' cannot be overridden on a stored panel: the recruitment it holds is %s, which is not what a %s design's recruitment counts; state the count you mean with prec_panel(n_recruit, target, ..., design = \"%s\")",
        if (identical(p$design, "fixed")) {
          "the whole issue to one cohort"
        } else {
          "the entrants per occasion"
        },
        if (identical(p$design, "fixed")) "rotating" else "fixed",
        if (identical(p$design, "fixed")) "rotating" else "fixed"
      ),
      call. = FALSE
    )
  }
  args <- list(
    n_recruit = .panel_recruit(x),
    target = x$target,
    retention = p$retention,
    resp_rate = p$resp_rate,
    design = p$design,
    target_wave = p$target_wave,
    assurance = p$assurance,
    # unlike `design`, this may be overridden: it selects which launch is
    # described and moves no stored quantity, the recruitment being fixed
    # before any of it is computed
    start = p$start
  )
  do.call(
    prec_panel.default,
    .roundtrip_args(args, dots, prec_panel.default)
  )
}

#' The recruitment count, whichever of the two names holds it
#'
#' The fixed and the rotating figures are different quantities and are
#' deliberately not stored under one name, so anything reading the count
#' back has to say which design it means.
#' @keywords internal
#' @noRd
.panel_recruit <- function(x) {
  if (identical(x$design, "fixed")) x$n_issued else x$n_entrants
}
