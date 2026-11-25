#' Precision of a joint constrained stratified allocation
#'
#' Evaluate achieved precision for every indicator and domain requirement under
#' a stratum allocation, or check a fitted [n_alloc_joint()] plan.
#'
#' @param frame For the default method, a stratum-level data frame as in
#'   [n_alloc_joint()]. For `svyplan_alloc_joint` objects, an allocation from
#'   [n_alloc_joint()].
#' @param n Numeric vector of stratum sample sizes (ultimate units) of length
#'   `nrow(frame)`. Required with a raw frame. For a cluster frame, the
#'   supplied sizes must imply whole PSU counts at the specified takes.
#'   See Details for the defaults when passing an allocation.
#' @param measures Long data frame of per-stratum indicator statistics, as in
#'   [n_alloc_joint()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param targets Optional precision targets data frame, as in
#'   [n_alloc_joint()].
#' @param objective Optional objective specification, as in [n_alloc_joint()].
#' @param budget Optional budget used with `objective` to compute budget
#'   residual.
#' @param objective_rule How `objective` combines its components, `"sum"`
#'   (default) or `"worst"`, as in [n_alloc_joint()].
#' @param psu Optional PSU register, as in [n_alloc_joint()].
#' @param certainty_cutoff Optional certainty probability cutoff, as in
#'   [n_alloc_joint()].
#' @param n_psu_per_zone Optional `1` or `2` PSUs per zone, as in
#'   [n_alloc_joint()].
#' @param unit_cost Optional scalar or length-`nrow(frame)` per-unit cost.
#' @param alpha Significance level, default `0.05`.
#' @param deff Default design effect multiplier (`> 0`).
#' @param resp_rate Default expected response rate in `(0, 1]`.
#' @param eligibility_rate Default share of issued units in scope in `(0, 1]`.
#' @param eligible_per_unit Default eligible units per sampled unit (`> 0`).
#'   Available for designs without a PSU register. As in [n_alloc_joint()],
#'   any of these four other than 1 is refused when every row a target or
#'   the objective uses takes the value from a `measures` or `frame`
#'   column, since the argument would be ignored.
#' @param df Optional degrees of freedom for MOE and RMOE targets.
#' @param min_n_stratum Optional minimum sample size per stratum.
#' @param fpc Finite population correction (`"unit"`, `"stage"`, or `"none"`).
#'   Leave the default for a PSU register.
#' @param plan Optional [svyplan()] object providing design defaults.
#'
#' @return A `svyplan_prec` object with `type = "alloc_joint"`, per-constraint
#'   results in `$detail`, stratum bound checks in `$bounds`, and optional
#'   `$objective`, `$objective_value`, and `$psu`.
#'
#' @details
#' * For an element or cluster fit, `prec_alloc_joint(fit)` assesses the
#'   continuous design. Supply `n = fit$detail$n_int` to assess its whole
#'   counts. A free-take fit uses its stored whole takes with explicit `n`.
#' * For a register fit, `prec_alloc_joint(fit)` assesses the stored whole
#'   field design, including its certainty split, remainder draw, and
#'   per-PSU takes. Supplying `n`, `frame`, `psu`, `certainty_cutoff`, or
#'   `n_psu_per_zone` reconstructs the assessment. A raw `frame` with `psu`
#'   needs explicit totals and a fixed `take_2` in the frame.
#'
#' Pass an assessment to [predict.svyplan()] to vary assumptions such as
#' response, variance, or cost. Leave `n` out of the scenario grid to hold the
#' allocation fixed. See `vignette("joint-allocation")` for worked designs.
#'
#' @family stratified design functions
#' @seealso [n_alloc_joint()] to reoptimize an assessment, [design_df()] for
#'   the field design's degrees of freedom, and [prec_alloc()] for
#'   classic single-indicator allocation precision.
#'
#' @examples
#' frame <- data.frame(
#'   stratum = c("North urban", "North rural", "South urban"),
#'   region = c("North", "North", "South"),
#'   N = c(1000, 1800, 1200)
#' )
#' measures <- data.frame(
#'   stratum = rep(frame$stratum, 2),
#'   name = rep(c("coverage", "income"), each = 3),
#'   p = c(0.5, 0.4, 0.6, NA, NA, NA),
#'   mean = c(NA, NA, NA, 50, 55, 60),
#'   sd = c(NA, NA, NA, 10, 12, 15)
#' )
#' targets <- data.frame(
#'   name = c("coverage", "income"),
#'   domain = c(".overall", "region"),
#'   level = c(NA, "North"),
#'   cv = c(0.06, NA),
#'   moe = c(NA, 3)
#' )
#' fit <- n_alloc_joint(frame, measures, targets = targets)
#' prec_alloc_joint(fit)
#' prec_alloc_joint(fit, n = fit$detail$n_int)
#'
#' # Assess the fixed-take field design at two response rates
#' cluster <- data.frame(
#'   stratum = "A", N = 100000, N_1 = 100, take_2 = 10, cost_1 = 50, cost_2 = 1
#' )
#' measure <- data.frame(stratum = "A", name = "coverage", p = 0.5, icc_1 = 0.05)
#' goal <- data.frame(name = "coverage", cv = 0.08)
#' cluster_fit <- n_alloc_joint(cluster, measure, targets = goal)
#' assessment <- prec_alloc_joint(cluster_fit, n = cluster_fit$detail$n_int)
#' predict(assessment, data.frame(resp_rate = c(1, 0.8)))
#'
#' @export
prec_alloc_joint <- function(frame, ...) {
  if (!missing(frame)) {
    .res <- .dispatch_plan(frame, "frame", prec_alloc_joint.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_alloc_joint")
}

#' @rdname prec_alloc_joint
#' @export
prec_alloc_joint.default <- function(
  frame,
  n = NULL,
  measures = NULL,
  ...,
  targets = NULL,
  objective = NULL,
  budget = NULL,
  objective_rule = c("sum", "worst"),
  psu = NULL,
  certainty_cutoff = NULL,
  n_psu_per_zone = NULL,
  unit_cost = NULL,
  alpha = 0.05,
  deff = 1,
  resp_rate = 1,
  eligibility_rate = 1,
  eligible_per_unit = 1,
  df = NULL,
  min_n_stratum = NULL,
  fpc = c("unit", "stage", "none"),
  plan = NULL
) {
  .plan <- .merge_plan_args(
    plan,
    prec_alloc_joint.default,
    match.call(),
    environment()
  )
  if (!is.null(.plan)) {
    return(do.call(prec_alloc_joint.default, c(.plan, list(...))))
  }
  .check_alloc_classic_args(list(...), "prec_alloc()")
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(frame, "frame")
  .refuse_unnumbered_stage_columns(measures, "measures")
  .refuse_unnumbered_stage_columns(psu, "psu")
  .check_register_args(certainty_cutoff, n_psu_per_zone, frame, psu)
  .check_eligible_scope(eligible_per_unit, frame, measures, psu)
  .refuse_joint_argument_beside_columns(
    frame, measures, targets, objective, alpha, df,
    resp_rate = resp_rate, deff = deff,
    eligibility_rate = eligibility_rate, eligible_per_unit = eligible_per_unit
  )
  fpc <- match.arg(fpc)
  objective_rule <- match.arg(objective_rule)
  .check_objective_rule(objective_rule, objective, NULL)

  if (is.null(measures)) {
    stop("joint precision assessment requires 'measures'", call. = FALSE)
  }
  if (is.null(targets) && is.null(objective)) {
    stop("'measures' and 'targets' must be supplied together", call. = FALSE)
  }
  if (!is.null(budget) && is.null(objective)) {
    stop(
      "'budget' requires 'objective' for joint precision assessment",
      call. = FALSE
    )
  }
  if (!is.null(budget)) {
    check_scalar(budget, "budget")
  }
  min_n_stratum <- .resolve_min_n_stratum(min_n_stratum, frame)
  .alloc_refuse_step(frame, "joint constrained allocation")

  if (!is.null(psu)) {
    if (!identical(fpc, "unit")) {
      stop(
        "'fpc' is not available with a PSU register: the certainty split and its size-proportional remainder need their own variance decomposition, which is not derived here. Drop 'fpc', or drop 'psu' and give 'N_1' instead",
        call. = FALSE
      )
    }
    .refuse_objective_rule_route(objective_rule, "register")
    return(.prec_alloc_psu(
      frame = frame,
      psu = psu,
      n = n,
      measures = measures,
      targets = targets,
      objective = objective,
      budget = budget,
      unit_cost = unit_cost,
      alpha = alpha,
      deff = deff,
      resp_rate = resp_rate,
      min_n_stratum = min_n_stratum,
      df = df,
      certainty_cutoff = certainty_cutoff,
      n_psu_per_zone = n_psu_per_zone,
      eligibility_rate = eligibility_rate
    ))
  }

  .prec_alloc_bethel(
    frame = frame,
    n = n,
    measures = measures,
    targets = targets,
    objective = objective,
    budget = budget,
    unit_cost = unit_cost,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    eligible_per_unit = eligible_per_unit,
    eligibility_rate = eligibility_rate,
    min_n_stratum = min_n_stratum,
    df = df,
    fpc = fpc,
    .finite_prop_var = TRUE,
    objective_rule = objective_rule
  )
}

#' @rdname prec_alloc_joint
#' @export
prec_alloc_joint.svyplan_alloc_joint <- function(frame, ...) {
  obj <- frame
  p <- obj$params
  dots <- list(...)
  if (!is.null(p$psu) && any(grepl("^take_[0-9]+$", names(dots)))) {
    stop(
      "stage take overrides cannot override a PSU register, whose takes come from the register itself",
      call. = FALSE
    )
  }
  take <- .alloc_take_override(p$frame, dots, isTRUE(p$stages > 1L))
  given <- function(nm) if (nm %in% names(dots)) dots[[nm]] else p[[nm]]
  .refuse_joint_argument_beside_columns(
    given("frame"), given("measures"), given("targets"), given("objective"),
    given("alpha"), given("df"),
    resp_rate = given("resp_rate"), deff = given("deff"),
    eligibility_rate = given("eligibility_rate"),
    eligible_per_unit = given("eligible_per_unit")
  )
  take_overridden <- !identical(take$frame, p$frame)
  p$frame <- take$frame
  dots <- take$dots
  dot_names <- names(dots)
  n_explicit <- length(dots) > 0L &&
    (is.null(dot_names) || any(!nzchar(dot_names)) || "n" %in% dot_names)
  if (!is.null(p$psu)) {
    design <- c("n", "frame", "psu", "certainty_cutoff", "n_psu_per_zone")
    if (!n_explicit && !any(dot_names %in% design)) {
      .roundtrip_args(list(), dots, .prec_alloc_psu)
      return(.prec_alloc_psu_design(obj, dots))
    }
    if (isTRUE(p$free_take) && !("frame" %in% dot_names)) {
      p$frame$take_2 <- obj$detail$take_2_int
      p$frame <- p$frame[, setdiff(names(p$frame),
                                   c("take_2_min", "take_2_max", "take_2_step")),
                         drop = FALSE]
    }
    held <- p$psu
    held$certainty <- obj$psu$certainty
    args <- list(
      frame = p$frame,
      psu = held,
      n = obj$detail$n,
      measures = p$measures,
      targets = p$targets,
      objective = p$objective,
      budget = p$budget,
      alpha = p$alpha,
      deff = p$deff,
      resp_rate = p$resp_rate,
      eligibility_rate = p[["eligibility_rate"]] %||% 1,
      df = p$df,
      unit_cost = p$unit_cost,
      min_n_stratum = p$min_n_stratum,
      certainty_cutoff = p$certainty_cutoff,
      n_psu_per_zone = p$n_psu_per_zone,
      .allow_fractional_stages = !n_explicit
    )
    out <- do.call(
      .prec_alloc_psu,
      .roundtrip_args(args, dots, .prec_alloc_psu)
    )
    if (!("psu" %in% dot_names)) out$params$psu <- p$psu
    out$params$control <- p[["control"]]
    return(out)
  }
  if (isTRUE(p$free_take) && !take_overridden && !("frame" %in% dot_names)) {
    if (!n_explicit) {
      args_ft <- list(
        frame = p$frame,
        n = obj$detail$n,
        take_2 = obj$detail$take_2,
        measures = p$measures,
        targets = p$targets,
        objective = p$objective,
        budget = p$budget,
        alpha = p$alpha,
        deff = p$deff,
        resp_rate = p$resp_rate,
        eligible_per_unit = p$eligible_per_unit %||% 1,
        eligibility_rate = p[["eligibility_rate"]] %||% 1,
        df = p$df,
        min_n_stratum = p$min_n_stratum,
        fpc = p$fpc %||% "unit",
        take_2_min = p[["take_2_min"]],
        take_2_max = p[["take_2_max"]],
        take_2_step = p[["take_2_step"]],
        objective_rule = p$objective_rule %||% "sum"
      )
      out <- do.call(
        .prec_alloc_freetake,
        .roundtrip_args(args_ft, dots, .prec_alloc_freetake)
      )
      out$params$control <- p[["control"]]
      return(out)
    }
    p$frame$take_2 <- obj$detail$take_2_int
    p$frame <- p$frame[, setdiff(names(p$frame),
                                 c("take_2_min", "take_2_max", "take_2_step")),
                       drop = FALSE]
  }
  args <- list(
    frame = p$frame,
    n = obj$detail$n,
    measures = p$measures,
    targets = p$targets,
    objective = p$objective,
    budget = p$budget,
    alpha = p$alpha,
    deff = p$deff,
    resp_rate = p$resp_rate,
    eligible_per_unit = p$eligible_per_unit %||% 1,
    eligibility_rate = p[["eligibility_rate"]] %||% 1,
    df = p$df,
    unit_cost = p$unit_cost,
    min_n_stratum = p$min_n_stratum,
    fpc = p$fpc %||% "unit",
    .allow_fractional_stages = !n_explicit,
    .finite_prop_var = TRUE,
    objective_rule = p$objective_rule %||% "sum"
  )
  out <- do.call(
    .prec_alloc_bethel,
    .roundtrip_args(args, dots, .prec_alloc_bethel)
  )
  out$params$control <- p[["control"]]
  out
}
