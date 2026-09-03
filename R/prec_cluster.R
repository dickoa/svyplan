#' Sampling precision for a multistage cluster allocation
#'
#' Compute the sampling error (SE, MOE, CV) for a given multistage sample
#' allocation. This is the inverse of [n_cluster()].
#'
#' @param n For the default method: numeric vector of per-stage sample
#'   sizes (`c(n_psu, n_per_psu)` for 2-stage or
#'   `c(n_psu, n_per_psu, n_per_ssu)` for 3-stage). Named vectors are accepted
#'   with stage names `n_psu`, `n_per_psu`, `n_per_ssu`.
#'   For `svyplan_cluster` objects: a cluster allocation from [n_cluster()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param icc Numeric vector of homogeneity measures (length = stages - 1),
#'   or a `svyplan_varcomp` object.
#' @param unit_relvar Unit relvariance (default 1).
#' @param var_ratio Ratio of the stage components' unit variance to the
#'   analysis variable's, default 1. A scalar names `var_ratio_psu`, and for
#'   three stages `var_ratio_ssu = var_ratio_psu * (1 - icc_psu)` follows from
#'   the decomposition. See [design_effect()].
#' @param resp_rate_psu Expected **PSU-level** response rate, in (0, 1\].
#'   Default 1 (no adjustment). The effective stage-1 size is
#'   `n * resp_rate_psu`. It describes clusters that cannot be worked, not
#'   nonresponse among the ultimate units inside a cluster.
#' @param resp_rate_ssu Expected SSU-level response rate, in (0, 1\].
#'   Three-stage designs only, default 1. The effective stage-2 size is
#'   `n[2] * resp_rate_ssu`.
#' @param resp_rate Expected ultimate-unit response rate, in (0, 1\].
#'   Default 1. It scales the final stage, so it also shrinks the realized
#'   cluster and therefore the clustering penalty. See [n_cluster()] for the
#'   decomposition.
#' @param indicators Optional data frame with one row per indicator, which
#'   switches the function to the several-indicators mode of Details. It
#'   carries the stage sizes per row in `n` and `n_per_psu`, plus
#'   `n_per_ssu` for a three-stage design, so `n` is left `NULL`.
#' @param domains Optional character vector naming domain columns in
#'   `indicators`. Precision is computed row by row either way, so naming
#'   domains leaves every `$detail` value unchanged. What it does is record
#'   the domain structure on the result, so that a round trip back to
#'   [n_cluster()] rebuilds the same design. There is no `domain_sampling`
#'   argument here, because that choice governs how per-domain requirements
#'   combine into one size, and this direction reads the sizes you already
#'   have.
#' @param stage_cost Optional per-stage costs, recorded on the result so that
#'   a later round trip to [n_cluster()] can re-solve the design. Costs do not
#'   enter the precision calculation.
#' @param plan Optional [svyplan()] object providing design defaults. In the
#'   several-indicators mode a profile supplies only the defaults that mode
#'   accepts, so `icc`, `unit_relvar` and `var_ratio` are left to the
#'   indicator columns.
#'
#' @return A `svyplan_prec` object with components `$se`, `$moe`, and `$cv`.
#'   Because the cluster model is parameterized with unit relvariance
#'   (`unit_relvar = S^2 / Y_bar^2`), only `$cv` is computable. The `$se` and
#'   `$moe` components are `NA`.
#'
#'   With `indicators`, per-indicator precision is in `$detail`: `.se`,
#'   `.moe`, `.rmoe` and `.cv`. All four are reported whatever target the
#'   design was sized against, since the achieved precision does not depend on
#'   how the requirement was written. `.rmoe` is measured against the row's
#'   own estimand, `p` for a proportion, `abs(mu)` for a mean and `abs(r)` for
#'   a ratio.
#'
#' @details
#' `prec_cluster()` is the inverse of [n_cluster()]: given per-stage
#' sample sizes, it computes the achieved precision. You can pass the
#' result of `n_cluster()` directly: `prec_cluster(n_cluster(...))`.
#'
#' Stage count is determined by `length(n)`.
#'
#' **2-stage** (Valliant et al., 2018, Eq. 9.2.23):
#' \deqn{CV = \sqrt{\frac{V \cdot k}{n_1 \cdot n_2} (1 + \delta (n_2 - 1))}}{CV = sqrt((V * k)/(n_1 * n_2) (1 + delta (n_2 - 1)))}
#'
#' **3-stage**:
#' \deqn{CV = \sqrt{\frac{V}{n_1 \cdot n_2 \cdot n_3} (k_1 \delta_1 n_2 n_3
#'   + k_2 (1 + \delta_2 (n_3 - 1)))}}{CV = sqrt(V/(n_1 * n_2 * n_3) (k_1 delta_1 n_2 n_3 + k_2 (1 + delta_2 (n_3 - 1))))}
#'
#' ## Several indicators at once
#'
#' `indicators` reads the precision of several indicators under one
#' allocation. It is a data frame with one row per indicator, carrying `n` and
#' `n_per_psu`, plus `n_per_ssu` for a three-stage design. The homogeneity and
#' indicator columns follow the schema [n_cluster()] documents, including the
#' derived three-stage `var_ratio_ssu`. The scalar arguments `n`, `icc`,
#' `unit_relvar` and `var_ratio` are refused in this mode, because the frame
#' carries each of them per row.
#'
#' `resp_rate_psu`, `resp_rate_ssu` and `resp_rate` are the defaults used
#' where a row's own column is absent or `NA`, as they are in [n_cluster()].
#' Precision is computed row by row, so an allocation returned by
#' `n_cluster(indicators = )` can be passed straight in.
#'
#' @references
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer. Ch. 9.
#'
#' @family cluster design functions
#' @seealso [n_cluster()] for the inverse operation, [varcomp()] for
#'   estimating variance components, and [prec_multi()] for several
#'   indicators without clustering.
#'
#' @examples
#' # Direct usage
#' prec_cluster(n = c(50, 12), icc = 0.05)
#' prec_cluster(n = c(50, 12, 8), icc = c(0.01, 0.05))
#'
#' # Round-trip from n_cluster
#' res <- n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05)
#' prec_cluster(res)
#'
#' # Several indicators under one allocation
#' indicators <- data.frame(
#'   name      = c("stunting", "anemia"),
#'   p         = c(0.30, 0.10),
#'   n         = c(60, 60),
#'   n_per_psu = c(12, 12),
#'   icc_psu   = c(0.02, 0.05)
#' )
#' prec_cluster(indicators = indicators)
#'
#' @export
prec_cluster <- function(n = NULL, ...) {
  # With 'n' absent, UseMethod() would dispatch on whatever argument came
  # first, so an indicator table would decide the method. The mode belongs to
  # 'indicators', which the default method reads.
  if (missing(n)) return(prec_cluster.default(n, ...))
  .res <- .dispatch_plan(n, "n", prec_cluster.default, ...)
  if (!is.null(.res)) return(.res)
  UseMethod("prec_cluster")
}

#' @rdname prec_cluster
#' @export
prec_cluster.default <- function(
  n = NULL,
  ...,
  icc = NULL,
  unit_relvar = NULL,
  var_ratio = NULL,
  resp_rate_psu = 1,
  resp_rate_ssu = 1,
  resp_rate = 1,
  indicators = NULL,
  domains = NULL,
  stage_cost = NULL,
  plan = NULL
) {
  plan <- .plan_for_cluster_mode(plan, !is.null(indicators))
  # Before the merge, which is the only point at which the profile's own
  # defaults are still distinguishable from the formals.
  .check_relvar_identified(
    icc,
    !is.null(unit_relvar) || "unit_relvar" %in% names(plan$defaults),
    "prec_cluster()"
  )
  .plan <- .merge_plan_args(plan, prec_cluster.default, match.call(), environment())
  if (!is.null(.plan)) return(do.call(prec_cluster.default, c(.plan, list(...))))
  .check_unused_dots(...)
  if (!is.null(indicators)) {
    .check_prec_cluster_indicator_args(n, icc, unit_relvar, var_ratio)
    return(.prec_cluster_indicators(
      indicators,
      domains = domains,
      stage_cost = stage_cost,
      resp_rate_psu = resp_rate_psu,
      resp_rate_ssu = resp_rate_ssu,
      resp_rate = resp_rate
    ))
  }
  if (!is.null(domains)) .stop_applies_to_indicators("domains")
  if (is.data.frame(n)) .stop_frame_in_first_slot("n", "per-stage sample sizes")
  unit_relvar <- unit_relvar %||% 1
  var_ratio <- var_ratio %||% 1
  if (is.null(icc))
    stop("'icc' is required (directly or via plan)", call. = FALSE)
  if (inherits(icc, "svyplan_varcomp")) {
    vc <- icc
    if (!is.null(vc$strata)) {
      stop(
        "stratified varcomp: merge its $strata columns into an n_alloc() frame, or pass one stratum's values",
        call. = FALSE
      )
    }
    icc <- vc$icc
    var_ratio <- vc$var_ratio
    if (!is.na(vc$unit_relvar)) unit_relvar <- vc$unit_relvar
  }

  if (!is.numeric(n) || length(n) < 2L) {
    stop("'n' must be a numeric vector of length >= 2", call. = FALSE)
  }
  if (length(n) > 3L) {
    stop(
      "4+ stage CV calculation is not supported; evaluate the first three stages and fold deeper stages (e.g. persons within households) into the 'deff' passed to prec_prop() or prec_mean()",
      call. = FALSE
    )
  }
  n <- .reorder_n_vec(n)
  if (anyNA(n)) {
    stop("'n' must not contain NA values", call. = FALSE)
  }
  if (any(!is.finite(n))) {
    stop("'n' must contain only finite values", call. = FALSE)
  }

  stages <- length(n)
  icc <- .reorder_stage_vec(icc, "icc")
  var_ratio <- .reorder_stage_vec(var_ratio, "var_ratio")
  check_icc(icc, expected_length = stages - 1L)
  check_resp_rate(resp_rate_psu, "resp_rate_psu")
  check_resp_rate(resp_rate_ssu, "resp_rate_ssu")
  check_resp_rate(resp_rate, "resp_rate")
  if (stages == 2L && !isTRUE(all.equal(resp_rate_ssu, 1))) {
    stop(
      "'resp_rate_ssu' is not applicable for 2-stage designs: the units inside a PSU are the ultimate ones, so their nonresponse is 'resp_rate'",
      call. = FALSE
    )
  }
  check_scalar(unit_relvar, "unit_relvar")
  if (
    !is.numeric(var_ratio) ||
      length(var_ratio) == 0L ||
      anyNA(var_ratio) ||
      any(var_ratio <= 0) ||
      any(!is.finite(var_ratio))
  ) {
    stop("'var_ratio' must contain positive finite values", call. = FALSE)
  }

  if (any(n <= 0)) {
    stop("all elements of 'n' must be positive", call. = FALSE)
  }
  if (!is.null(stage_cost)) {
    check_stage_cost(stage_cost)
    stage_cost <- .reorder_stage_cost(stage_cost)
    if (length(stage_cost) != stages) {
      stop(
        sprintf("'stage_cost' must have length %d", stages),
        call. = FALSE
      )
    }
    names(stage_cost) <- c("cost_psu", "cost_ssu", "cost_tsu")[seq_len(stages)]
  }

  # Each stage keeps the share of its own units that respond. The last stage
  # always carries the ultimate-unit rate, whether the design has two stages
  # or three.
  rates <- if (stages == 2L) {
    c(resp_rate_psu, resp_rate)
  } else {
    c(resp_rate_psu, resp_rate_ssu, resp_rate)
  }
  n_eff <- n * rates
  .check_expected_take(n[[stages]], resp_rate,
                       if (stages == 2L) "n_per_psu" else "n_per_ssu",
                       "resp_rate")
  if (stages == 3L) {
    .check_expected_take(n[[2L]], resp_rate_ssu, "n_per_psu", "resp_rate_ssu")
  }

  if (stages == 2L) {
    var_ratio <- rep_len(var_ratio, 1L)
    cv_val <- unname(.cv_cluster_2stage(n_eff, icc, unit_relvar, var_ratio))
  } else {
    var_ratio <- .stage_k_pair(var_ratio, icc)
    cv_val <- unname(.cv_cluster_3stage(n_eff, icc, unit_relvar, var_ratio))
  }

  params <- list(
    n = n,
    icc = icc,
    unit_relvar = unit_relvar,
    var_ratio = var_ratio,
    resp_rate_psu = resp_rate_psu,
    resp_rate_ssu = resp_rate_ssu,
    resp_rate = resp_rate,
    stages = stages
  )
  if (!is.null(stage_cost)) params$stage_cost <- stage_cost
  .new_svyplan_prec(
    se = NA_real_,
    moe = NA_real_,
    cv = cv_val,
    type = "cluster",
    params = params
  )
}

#' @rdname prec_cluster
#' @export
prec_cluster.svyplan_cluster <- function(n, ...) {
  x <- n
  if (!is.null(x[["indicators"]])) {
    return(.prec_cluster_from_indicator_fit(x, list(...)))
  }
  p <- x$params
  args <- list(
    n = x$n,
    icc = p$icc,
    unit_relvar = p$unit_relvar,
    var_ratio = p$var_ratio,
    resp_rate_psu = p$resp_rate_psu %||% 1,
    resp_rate_ssu = p$resp_rate_ssu %||% 1,
    resp_rate = p$resp_rate %||% 1,
    stage_cost = p$stage_cost
  )
  out <- do.call(
    prec_cluster.default,
    .roundtrip_args(args, list(...), prec_cluster.default)
  )
  if (!is.null(p$budget)) {
    out$params$budget <- p$budget
  }
  if (!is.null(p$n_psu)) {
    out$params$n_psu <- p$n_psu
  }
  if (!is.null(p$n_per_psu)) {
    out$params$n_per_psu <- p$n_per_psu
  }
  if (!is.null(p$n_per_ssu)) {
    out$params$n_per_ssu <- p$n_per_ssu
  }
  if (!is.null(p$fixed_cost)) {
    out$params$fixed_cost <- p$fixed_cost
  }
  out
}

#' Refuse the scalar inputs the indicator frame carries per row
#' @keywords internal
#' @noRd
.check_prec_cluster_indicator_args <- function(n, icc, unit_relvar, var_ratio) {
  cols <- c(
    n = "n", icc = "icc_psu", unit_relvar = "unit_relvar",
    var_ratio = "var_ratio_psu"
  )
  supplied <- !vapply(
    list(n = n, icc = icc, unit_relvar = unit_relvar, var_ratio = var_ratio),
    is.null,
    logical(1L)
  )
  if (any(supplied)) .stop_carried_by_column(names(cols)[supplied], cols)
  invisible(TRUE)
}

#' Read a several-indicators cluster allocation back as precision
#'
#' The fit carries its own stage sizes, so the frame handed to the engine is
#' the one the allocation was solved from with the realized sizes written
#' back onto it, per domain where the design has domains.
#' @keywords internal
#' @noRd
.prec_cluster_from_indicator_fit <- function(x, dots) {
  tgt <- x$indicators
  for (rate in intersect(
    c("resp_rate_psu", "resp_rate_ssu", "resp_rate"),
    names(dots)
  )) {
    tgt[[rate]] <- NA_real_
  }
  tgt$cv <- NULL
  tgt$moe <- NULL
  tgt$n <- x$n[1L]
  if (x$stages >= 2L) {
    tgt$n_per_psu <- x$n[2L]
  }
  if (x$stages >= 3L) {
    tgt$n_per_ssu <- x$n[3L]
  }

  dom_cols <- x$params$domain_cols %||% character(0)
  if (!is.null(x$domains) && length(dom_cols) > 0L) {
    dom <- x$domains
    tgt_key <- .domain_key(tgt, dom_cols)
    dom_key <- .domain_key(dom, dom_cols)
    dom_idx <- match(tgt_key, dom_key)
    tgt$n <- dom$n_psu[dom_idx]
    if (x$stages >= 2L) {
      tgt$n_per_psu <- dom$n_per_psu[dom_idx]
    }
    if (x$stages >= 3L) tgt$n_per_ssu <- dom$n_per_ssu[dom_idx]
  }

  stage_cost <- x$params$stage_cost
  args <- list(
    indicators = tgt,
    stages = x$stages,
    stage_cost = stage_cost,
    domain_cols = dom_cols
  )
  allowed <- c(
    "stage_cost", "domains",
    "resp_rate_psu", "resp_rate_ssu", "resp_rate"
  )
  dot_names <- names(dots) %||% rep("", length(dots))
  unknown <- setdiff(dot_names, allowed)
  if (length(unknown) > 0L || any(!nzchar(dot_names))) {
    do.call(.check_unused_dots, dots)
  }
  if ("stage_cost" %in% names(dots)) {
    check_stage_cost(dots$stage_cost)
    override_cost <- .reorder_stage_cost(dots$stage_cost)
    if (length(override_cost) != x$stages) {
      stop(
        sprintf("'stage_cost' must have length %d", x$stages),
        call. = FALSE
      )
    }
    args$stage_cost <- override_cost
  }
  if ("domains" %in% names(dots)) {
    if (!is.null(dots$domains) &&
        (!is.character(dots$domains) || anyNA(dots$domains) ||
         any(!dots$domains %in% names(tgt)))) {
      stop("'domains' must name columns in indicators", call. = FALSE)
    }
    args$domain_cols <- dots$domains %||% character(0)
  }
  for (rate in intersect(
    c("resp_rate_psu", "resp_rate_ssu", "resp_rate"),
    names(dots)
  )) {
    check_resp_rate(dots[[rate]], rate)
    if (rate == "resp_rate_ssu" && x$stages == 2L && dots[[rate]] != 1) {
      stop("'resp_rate_ssu' is not applicable for 2-stage designs", call. = FALSE)
    }
    args[[rate]] <- dots[[rate]]
  }
  res <- do.call(.prec_multi_cluster, args)
  for (p in c("budget", "n_psu", "n_per_psu", "n_per_ssu", "joint",
              "allocation", "fixed_cost", "min_n_domain", "mode")) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res
}

#' @keywords internal
#' @noRd
.cv_cluster_2stage <- function(n, icc, unit_relvar, var_ratio) {
  n1 <- n[1L]
  n2 <- n[2L]
  sqrt(unit_relvar / (n1 * n2) * var_ratio * (1 + icc * (n2 - 1)))
}

#' @keywords internal
#' @noRd
.cv_cluster_3stage <- function(n, icc, unit_relvar, var_ratio) {
  n1 <- n[1L]
  n2 <- n[2L]
  n3 <- n[3L]
  delta1 <- icc[1L]
  delta2 <- icc[2L]
  k1 <- var_ratio[1L]
  k2 <- var_ratio[2L]
  sqrt(
    unit_relvar /
      (n1 * n2 * n3) *
      (k1 * delta1 * n2 * n3 + k2 * (1 + delta2 * (n3 - 1)))
  )
}
