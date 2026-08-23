#' Multi-indicator sampling precision
#'
#' Compute the sampling error (SE, MOE, CV) for multiple survey indicators
#' given a sample size. This is the inverse of [n_multi()].
#'
#' @param indicators For the default method: a data frame where **each row
#'   is one survey indicator**, in the same format as [n_multi()] but
#'   with an additional `n` column giving the sample size. This lets you
#'   answer: "given this sample size, what precision do I get for each
#'   indicator?"
#'
#'   At minimum, each row needs:
#'   \itemize{
#'     \item `p`, `var`, **or** the ratio quartet `r`, `cv_num`, `cv_den`
#'       and `component_cor`: what you are measuring (see [n_multi()]).
#'       Exactly one estimand per row.
#'     \item `n`: the sample size to evaluate.
#'   }
#'
#'   See the Details section for the full column reference.
#'
#'   For `svyplan_n` objects: a result from [n_multi()].
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param domains Character vector of column names in `indicators` to treat
#'   as domain variables, or `NULL` (default) for no domains. All names
#'   must exist in `indicators`. Domain columns are preserved in the result
#'   for round-trip conversion back to [n_multi()].
#' @param prop_method Proportion CI method, one of `"wald"`
#'   (default), `"wilson"`, `"logodds"`, or `"beta"`. This is passed to
#'   [prec_prop()] for proportion rows and ignored for mean rows.
#'   An optional `prop_method` column in `indicators` overrides this default
#'   on a per-row basis.
#' @param resp_rate Default expected response rate at the ultimate unit, in
#'   (0, 1\]. Used for rows whose `resp_rate` column is absent or `NA`, and
#'   a non-missing row value overrides it.
#' @param plan A [svyplan()] profile providing default design parameters.
#'
#' @return A `svyplan_prec` object with a `$detail` data frame containing
#'   per-indicator precision: `.se`, `.moe`, `.rmoe`, and `.cv`. `.rmoe`
#'   is `.moe` as a fraction of the row's own estimand, `p` for a
#'   proportion, `abs(mu)` for a mean and `abs(r)` for a ratio, and `NA` for
#'   a row carrying none of them.
#'
#' @details
#' ## Building the indicators data frame
#'
#' The `indicators` data frame uses the same structure as [n_multi()],
#' with the addition of a required `n` column specifying the sample
#' size to evaluate. A minimal example:
#'
#' ```
#' indicators <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   n    = c(400, 400, 400)
#' )
#' ```
#'
#' See [n_multi()] for a detailed guide on constructing indicator rows
#' (choosing between `p` and `var`, setting per-row design effects, etc.).
#'
#' ## Column reference
#'
#' \describe{
#'   \item{`name`}{Indicator label (optional).}
#'   \item{`p`}{Expected proportion, in (0, 1). One of `p` or `var`
#'     per row (see [n_multi()]).}
#'   \item{`var`}{Population variance. One of `p`, `var` or `r` per row.}
#'   \item{`mu`}{Population mean. Required for CV output when `var`
#'     is specified, because CV = SE / mean.}
#'   \item{`r`}{Anticipated ratio of two totals, marking a ratio row. It
#'     requires `cv_num`, `cv_den` and `component_cor` alongside it, and no
#'     `unit_relvar`, which the moments derive. See [prec_ratio()].}
#'   \item{`cv_num`, `cv_den`}{Coefficients of variation of a ratio row's
#'     numerator and denominator, strictly positive.}
#'   \item{`component_cor`}{Correlation between a ratio row's numerator and
#'     denominator across units, in \[-1, 1\].}
#'   \item{`n`}{Sample size to evaluate (**required**), with one value per
#'     indicator row.}
#'   \item{`alpha`}{Significance level (default 0.05).}
#'   \item{`deff`}{Design effect multiplier (default 1).}
#'   \item{`N`}{Population size (default `Inf`).}
#'   \item{`prop_method`}{Proportion CI method: `"wald"` (default),
#'     `"wilson"`, `"logodds"`, or `"beta"`. Only for rows with `p`.}
#'   \item{`df`}{Degrees of freedom of the variance estimator, typically
#'     sampled PSUs minus strata, and available from [design_df()]. It
#'     switches that row's interval quantile from normal to t, under every
#'     proportion method and on mean rows alike. `NA` (the default) applies
#'     no adjustment.}
#'   \item{`resp_rate`}{Expected response rate at the ultimate unit
#'     (default 1). [prec_multi_cluster()] spends its rate at stage 1 and
#'     names the column `resp_rate_psu`.}
#' }
#'
#' Domain columns are specified via the `domains` parameter.
#'
#' `prec_multi()` delegates proportion rows to [prec_prop()], mean rows to
#' [prec_mean()] and ratio rows to [prec_ratio()]. Use `prop_method` or a
#' `indicators$prop_method` column to choose `"wald"`, `"wilson"`,
#' `"logodds"` or `"beta"` for proportion rows. It is ignored elsewhere.
#'
#' @family precision functions
#' @seealso [n_multi()] for the inverse, [prec_multi_cluster()] for
#'   multistage cluster designs, and [prec_prop()] and [prec_mean()] for
#'   single-indicator precision.
#'
#' @examples
#' # Simple mode: precision for three indicators at n = 400
#' indicators <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   n    = c(400, 400, 400)
#' )
#' prec_multi(indicators)
#'
#' # Wilson precision for a rare proportion
#' prec_multi(data.frame(p = 0.05, n = 400), prop_method = "wilson")
#'
#' @export
prec_multi <- function(indicators, ...) {
  if (!missing(indicators)) {
    .res <- .dispatch_plan(indicators, "indicators", prec_multi.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_multi")
}

#' @rdname prec_multi
#' @export
prec_multi.default <- function(
  indicators,
  ...,
  domains = NULL,
  prop_method = c("wald", "wilson", "logodds", "beta"),
  resp_rate = 1,
  plan = NULL
) {
  merged <- .merge_plan_args(
    plan,
    prec_multi.default,
    match.call(),
    environment()
  )
  if (!is.null(merged)) {
    return(do.call(prec_multi.default, c(merged, list(...))))
  }
  .check_multi_split_args(list(...), "prec_multi_cluster()")
  .check_unused_dots(...)
  # An unresolved default arrives either as a missing argument or, through
  # the plan-merge path, as the full choice vector itself.
  if (identical(prop_method, c("wald", "wilson", "logodds", "beta"))) {
    prop_method <- "wald"
  }
  if (
    !is.character(prop_method) ||
      length(prop_method) != 1L ||
      is.na(prop_method) ||
      !prop_method %in% c("wald", "wilson", "logodds", "beta")
  ) {
    stop(
      "'prop_method' must be one of 'wald', 'wilson', 'logodds', or 'beta'",
      call. = FALSE
    )
  }
  check_resp_rate(resp_rate)

  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }

  indicators <- .indicators_var_from_sd(indicators, domains)
  .check_indicator_columns(indicators, domains)
  .check_resp_rate_column(indicators, FALSE)
  .stop_min_cases_column(
    indicators,
    "prec_multi() reads the sizes you already have; a case floor is a sizing constraint, so it belongs to n_multi()"
  )

  if (!"n" %in% names(indicators)) {
    stop("'indicators' must contain an 'n' column for prec_multi", call. = FALSE)
  }
  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(indicators))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in indicators: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  domain_cols <- domains %||% character(0)
  .prec_multi_simple(indicators, prop_method = prop_method,
                     resp_rate = resp_rate,
                     domain_cols = domain_cols)
}

#' Multi-indicator precision for cluster designs
#'
#' Compute achieved precision for several indicators under a two- or
#' three-stage cluster allocation. This is the inverse of
#' [n_multi_cluster()].
#'
#' @param indicators For the default method, a non-empty data frame with one row
#'   per indicator. It must contain `n` and `n_per_psu`. Include `n_per_ssu` for
#'   a three-stage design. Cluster homogeneity and indicator columns follow
#'   the schema used by [n_multi_cluster()], including the derived
#'   three-stage `var_ratio_ssu`. For `svyplan_cluster` methods, an allocation
#'   returned by [n_multi_cluster()].
#' @param ... Additional arguments passed to methods. Unused arguments are
#'   rejected.
#' @param domains Optional character vector naming domain columns in
#'   `indicators`.
#' @param stage_cost Optional per-stage costs to retain for a later round trip
#'   to [n_multi_cluster()]. Costs do not enter the precision calculation.
#' @param resp_rate_psu Default expected PSU response rate, in (0, 1\].
#'   Used where the indicator column is absent or `NA`.
#' @param resp_rate_ssu Default expected SSU response rate for a three-stage
#'   design, in (0, 1\]. It is not applicable to a two-stage design.
#' @param resp_rate Default expected ultimate-unit response rate, in (0, 1\].
#'   Non-missing indicator columns override these three defaults row by row.
#' @param plan Optional [svyplan()] profile providing design metadata.
#'
#' @return A `svyplan_prec` object with per-indicator cluster precision in
#'   `$detail`: `.se`, `.moe`, `.rmoe`, and `.cv`, with `.rmoe` measured
#'   against the row's `p` or `abs(mu)`.
#'
#' @examples
#' indicators <- data.frame(
#'   name = c("stunting", "anemia"),
#'   p = c(0.30, 0.10),
#'   n = c(60, 60),
#'   n_per_psu = c(12, 12),
#'   icc_psu = c(0.02, 0.05)
#' )
#' prec_multi_cluster(indicators)
#'
#' @family precision functions
#' @seealso [n_multi_cluster()] for the inverse, [prec_multi()] for the
#'   single-stage counterpart, and [prec_cluster()] for one indicator.
#'
#' @export
prec_multi_cluster <- function(indicators, ...) {
  if (!missing(indicators)) {
    .res <- .dispatch_plan(
      indicators,
      "indicators",
      prec_multi_cluster.default,
      ...
    )
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_multi_cluster")
}

#' @rdname prec_multi_cluster
#' @export
prec_multi_cluster.default <- function(
  indicators,
  ...,
  domains = NULL,
  stage_cost = NULL,
  resp_rate_psu = 1,
  resp_rate_ssu = 1,
  resp_rate = 1,
  plan = NULL
) {
  merged <- .merge_plan_args(
    plan,
    prec_multi_cluster.default,
    match.call(),
    environment()
  )
  if (!is.null(merged)) {
    return(do.call(prec_multi_cluster.default, c(merged, list(...))))
  }
  .check_unused_dots(...)

  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }
  indicators <- .indicators_var_from_sd(indicators, domains)
  .check_indicator_columns(indicators, domains)
  .check_ratio_rows(indicators)
  .check_estimand_markers(indicators)
  .stop_min_cases_column(
    indicators,
    "prec_multi_cluster() reads the sizes you already have; a case floor is a sizing constraint, so it belongs to n_multi()"
  )

  if (!"n" %in% names(indicators)) {
    stop("'indicators' must contain an 'n' column", call. = FALSE)
  }
  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(indicators))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in indicators: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(stage_cost)) {
    check_stage_cost(stage_cost)
    stage_cost <- .reorder_stage_cost(stage_cost)
  }
  stages <- if (!is.null(stage_cost)) {
    length(stage_cost)
  } else if ("n_per_ssu" %in% names(indicators)) {
    3L
  } else {
    2L
  }
  .check_resp_rate_column(indicators, TRUE, stages)
  check_resp_rate(resp_rate_psu, "resp_rate_psu")
  check_resp_rate(resp_rate_ssu, "resp_rate_ssu")
  check_resp_rate(resp_rate)
  if (stages == 2L && resp_rate_ssu != 1) {
    stop("'resp_rate_ssu' is not applicable for 2-stage designs", call. = FALSE)
  }
  if (!is.null(stage_cost)) {
    target_stages <- if ("n_per_ssu" %in% names(indicators)) 3L else 2L
    if (target_stages == 3L && length(stage_cost) != 3L) {
      stop(
        sprintf(
          "'stage_cost' has length %d but target columns describe a %d-stage design",
          length(stage_cost),
          target_stages
        ),
        call. = FALSE
      )
    }
  }

  .prec_multi_cluster(
    indicators,
    stages = stages,
    stage_cost = stage_cost,
    domain_cols = domains %||% character(0),
    resp_rate_psu = resp_rate_psu,
    resp_rate_ssu = resp_rate_ssu,
    resp_rate = resp_rate
  )
}

#' @keywords internal
#' @noRd
.prec_multi_simple <- function(indicators, prop_method = "wald",
                              resp_rate = 1,
                              domain_cols = character(0)) {
  if (!"alpha" %in% names(indicators)) {
    indicators$alpha <- 0.05
  }
  if (!"deff" %in% names(indicators)) {
    indicators$deff <- 1
  }
  if (!"N" %in% names(indicators)) {
    indicators$N <- Inf
  }
  if (!"resp_rate" %in% names(indicators)) {
    indicators$resp_rate <- resp_rate
  } else {
    indicators$resp_rate[is.na(indicators$resp_rate)] <- resp_rate
  }
  if (!"prop_method" %in% names(indicators)) {
    indicators$prop_method <- prop_method
  } else {
    indicators$prop_method[is.na(indicators$prop_method)] <- prop_method
  }
  if (!"df" %in% names(indicators)) {
    indicators$df <- NA_real_
  }

  .validate_common_columns(indicators)

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)
  .check_ratio_rows(indicators)
  if (!has_p && !has_var && !any(is_ratio)) {
    stop("'indicators' must contain a 'p', 'var', or 'r' column",
         call. = FALSE)
  }
  # Exactly one estimand per row, the same rule n_multi() applies, so a table
  # the sizing path refuses cannot be evaluated by the precision path.
  .check_estimand_markers(indicators)

  if (has_p) {
    p_vals <- indicators$p[!is.na(indicators$p)]
    if (any(p_vals <= 0 | p_vals >= 1)) {
      stop("all 'p' values must be in (0, 1)", call. = FALSE)
    }
  }
  if (has_var) {
    var_vals <- indicators$var[!is.na(indicators$var)]
    if (any(var_vals <= 0) || any(!is.finite(var_vals))) {
      stop("all 'var' values must be positive and finite", call. = FALSE)
    }
  }
  if ("mu" %in% names(indicators)) {
    mu_vals <- indicators$mu[!is.na(indicators$mu)]
    if (any(mu_vals == 0) || any(!is.finite(mu_vals))) {
      stop("'mu' values must be finite and non-zero", call. = FALSE)
    }
  }
  method_vals <- indicators$prop_method[!is.na(indicators$prop_method)]
  bad_methods <- !method_vals %in% c("wald", "wilson", "logodds", "beta")
  if (any(bad_methods)) {
    stop(
      "'prop_method' values must be one of 'wald', 'wilson', 'logodds', or 'beta'",
      call. = FALSE
    )
  }

  nr <- nrow(indicators)
  has_mu <- "mu" %in% names(indicators)

  row_labels <- if ("name" %in% names(indicators)) {
    indicators$name
  } else {
    paste("indicator", seq_len(nr))
  }
  .check_gross_n(indicators$n, indicators$N, label = row_labels)

  se_vec <- numeric(nr)
  moe_vec <- numeric(nr)
  cv_vec <- numeric(nr)

  for (i in seq_len(nr)) {
    is_prop <- has_p && !is.na(indicators$p[i])

    if (is_ratio[i]) {
      res_i <- prec_ratio.default(
        r = indicators$r[i],
        n = indicators$n[i],
        cv_num = indicators$cv_num[i],
        cv_den = indicators$cv_den[i],
        component_cor = indicators$component_cor[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        df = .row_df(indicators, i)
      )
    } else if (is_prop) {
      res_i <- prec_prop.default(
        p = indicators$p[i],
        n = indicators$n[i],
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        method = indicators$prop_method[i],
        df = .row_df(indicators, i)
      )
    } else {
      res_i <- prec_mean.default(
        var = indicators$var[i],
        n = indicators$n[i],
        mu = if (has_mu && !is.na(indicators$mu[i])) indicators$mu[i] else NULL,
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        df = .row_df(indicators, i)
      )
    }
    se_vec[i] <- res_i$se
    moe_vec[i] <- res_i$moe
    cv_vec[i] <- res_i$cv
  }

  labels <- if ("name" %in% names(indicators)) indicators$name else seq_len(nr)

  detail <- data.frame(
    name = labels,
    .se = se_vec,
    .moe = moe_vec,
    .rmoe = .rmoe_from_moe(moe_vec, .indicator_estimand(indicators, nr)),
    .cv = cv_vec
  )

  .new_svyplan_prec(
    se = se_vec,
    moe = moe_vec,
    cv = cv_vec,
    type = "multi",
    params = list(
      indicators = indicators,
      domain_cols = domain_cols,
      design = "simple"
    ),
    detail = detail
  )
}

#' @keywords internal
#' @noRd
.prec_multi_cluster <- function(indicators, stages, stage_cost = NULL,
                                domain_cols = character(0),
                                mode = "cv",
                                resp_rate_psu = 1,
                                resp_rate_ssu = 1,
                                resp_rate = 1) {
  if (!"alpha" %in% names(indicators)) {
    indicators$alpha <- 0.05
  }
  if (!"resp_rate_psu" %in% names(indicators)) {
    indicators$resp_rate_psu <- resp_rate_psu
  } else {
    indicators$resp_rate_psu[is.na(indicators$resp_rate_psu)] <- resp_rate_psu
  }
  if (!"resp_rate" %in% names(indicators)) {
    indicators$resp_rate <- resp_rate
  } else {
    indicators$resp_rate[is.na(indicators$resp_rate)] <- resp_rate
  }
  if (stages == 3L) {
    if (!"resp_rate_ssu" %in% names(indicators)) {
      indicators$resp_rate_ssu <- resp_rate_ssu
    } else {
      indicators$resp_rate_ssu[is.na(indicators$resp_rate_ssu)] <-
        resp_rate_ssu
    }
  }
  if (!"var_ratio_psu" %in% names(indicators)) {
    indicators$var_ratio_psu <- 1
  }
  # var_ratio_ssu is derived below, once var_ratio_psu and icc_psu have been validated.
  var_ratio_ssu_supplied <- "var_ratio_ssu" %in% names(indicators)
  if (!var_ratio_ssu_supplied) {
    indicators$var_ratio_ssu <- 1
  }

  .validate_common_columns(indicators)

  if (!"unit_relvar" %in% names(indicators)) {
    indicators$unit_relvar <- NA_real_
  }
  indicators$unit_relvar <- .derive_unit_relvar(indicators, require_all = TRUE)

  rv_check <- indicators$unit_relvar[!is.na(indicators$unit_relvar)]
  if (
    length(rv_check) > 0L &&
      (any(rv_check <= 0) || any(!is.finite(rv_check)))
  ) {
    stop("'unit_relvar' values must be positive and finite", call. = FALSE)
  }
  if (any(indicators$var_ratio_psu <= 0) || any(!is.finite(indicators$var_ratio_psu))) {
    stop("'var_ratio_psu' values must be positive and finite", call. = FALSE)
  }
  if (any(indicators$var_ratio_ssu <= 0) || any(!is.finite(indicators$var_ratio_ssu))) {
    stop("'var_ratio_ssu' values must be positive and finite", call. = FALSE)
  }

  nr <- nrow(indicators)
  labels <- if ("name" %in% names(indicators)) indicators$name else seq_len(nr)

  n1 <- indicators$n
  if (
    stages >= 2L &&
      (!"n_per_psu" %in% names(indicators) ||
        anyNA(indicators$n_per_psu) ||
        any(indicators$n_per_psu <= 0))
  ) {
    stop(
      "'n_per_psu' column is required for multistage precision (positive, no NA)",
      call. = FALSE
    )
  }
  if (
    stages == 3L &&
      (!"n_per_ssu" %in% names(indicators) ||
        anyNA(indicators$n_per_ssu) ||
        any(indicators$n_per_ssu <= 0))
  ) {
    stop(
      "'n_per_ssu' column is required for 3-stage precision (positive, no NA)",
      call. = FALSE
    )
  }
  if (!"icc_psu" %in% names(indicators)) {
    stop(
      "'icc_psu' column is required for multistage precision",
      call. = FALSE
    )
  }
  if (!is.numeric(indicators$icc_psu)) {
    stop("'icc_psu' must be numeric", call. = FALSE)
  }
  if (anyNA(indicators$icc_psu)) {
    stop("'icc_psu' must not contain NA values", call. = FALSE)
  }
  if (any(indicators$icc_psu < 0 | indicators$icc_psu > 1)) {
    stop("'icc_psu' values must be in [0, 1]", call. = FALSE)
  }
  if (stages == 3L) {
    if (!"icc_ssu" %in% names(indicators)) {
      stop(
        "'icc_ssu' column is required for 3-stage precision",
        call. = FALSE
      )
    }
    if (!is.numeric(indicators$icc_ssu)) {
      stop("'icc_ssu' must be numeric", call. = FALSE)
    }
    if (anyNA(indicators$icc_ssu)) {
      stop("'icc_ssu' must not contain NA values", call. = FALSE)
    }
    if (any(indicators$icc_ssu < 0 | indicators$icc_ssu > 1)) {
      stop("'icc_ssu' values must be in [0, 1]", call. = FALSE)
    }
  }

  n2 <- indicators$n_per_psu
  n3 <- if (stages == 3L) indicators$n_per_ssu else rep(NA_real_, nr)

  rr <- indicators$resp_rate_psu
  rs <- .indicator_rate(indicators, "resp_rate_ssu", nr)
  ru <- .indicator_rate(indicators, "resp_rate", nr)
  n1_eff <- n1 * rr

  cv_vec <- numeric(nr)
  delta1 <- indicators$icc_psu
  delta2 <- if ("icc_ssu" %in% names(indicators)) {
    indicators$icc_ssu
  } else {
    rep(0, nr)
  }
  unit_relvar <- indicators$unit_relvar
  k1 <- indicators$var_ratio_psu
  # var_ratio_ssu is the within-PSU counterpart of var_ratio_psu and is fixed by the
  # decomposition: var_ratio_ssu = var_ratio_psu * (1 - icc_psu). See .var_ratio_ssu_default().
  if (stages == 3L) {
    implied <- .var_ratio_ssu_default(k1, delta1)
    if (var_ratio_ssu_supplied) {
      indicators$var_ratio_ssu[is.na(indicators$var_ratio_ssu)] <- implied[is.na(indicators$var_ratio_ssu)]
    } else {
      indicators$var_ratio_ssu <- implied
    }
  }
  k2 <- indicators$var_ratio_ssu

  # Stage sizes arrive gross; the variance reads what each stage realizes.
  n2r <- if (stages == 2L) n2 * ru else n2 * rs
  n3r <- if (stages == 3L) n3 * ru else n3
  for (i in seq_len(nr)) {
    if (stages == 2L) {
      cv_vec[i] <- sqrt(
        unit_relvar[i] * k1[i] / (n1_eff[i] * n2r[i]) * (1 + delta1[i] * (n2r[i] - 1))
      )
    } else {
      cv_vec[i] <- sqrt(
        unit_relvar[i] /
          (n1_eff[i] * n2r[i] * n3r[i]) *
          (k1[i] *
            delta1[i] *
            n2r[i] *
            n3r[i] +
            k2[i] * (1 + delta2[i] * (n3r[i] - 1)))
      )
    }
  }

  se_vec <- rep(NA_real_, nr)
  moe_vec <- rep(NA_real_, nr)

  # The inverse of .convert_moe_to_cv(): the cluster model delivers a sampling
  # CV, and each row's own interval method turns that back into a margin of
  # error. Reading moe as z * se would assume the Wald half-width for all four.
  if (identical(mode, "moe")) {
    has_p <- "p" %in% names(indicators)
    has_mu <- "mu" %in% names(indicators)
    has_method <- "prop_method" %in% names(indicators)
    is_ratio <- .is_ratio_row(indicators)
    for (i in seq_len(nr)) {
      df_i <- .row_df(indicators, i)
      if (is_ratio[i]) {
        se_vec[i] <- cv_vec[i] * abs(indicators$r[i])
        moe_vec[i] <- .q_alpha(indicators$alpha[i], df_i) * se_vec[i]
      } else if (has_p && !is.na(indicators$p[i])) {
        p_i <- indicators$p[i]
        se_vec[i] <- cv_vec[i] * p_i
        method_i <- if (has_method) indicators$prop_method[i] else "wald"
        moe_vec[i] <- .prec_engine_prop(
          p_i, (1 - p_i) / (p_i * cv_vec[i]^2), indicators$alpha[i],
          Inf, 1, 1, method_i, df_i
        )$moe
      } else if (has_mu && !is.na(indicators$mu[i])) {
        se_vec[i] <- cv_vec[i] * abs(indicators$mu[i])
        moe_vec[i] <- .q_alpha(indicators$alpha[i], df_i) * se_vec[i]
      }
    }
  }

  detail <- data.frame(
    name = labels,
    .se = se_vec,
    .moe = moe_vec,
    .rmoe = .rmoe_from_moe(moe_vec, .indicator_estimand(indicators, nr)),
    .cv = cv_vec
  )

  .new_svyplan_prec(
    se = se_vec,
    moe = moe_vec,
    cv = cv_vec,
    type = "multi",
    params = list(
      indicators = indicators,
      stage_cost = stage_cost,
      domain_cols = domain_cols,
      design = "cluster",
      stages = stages
    ),
    detail = detail
  )
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_n <- function(indicators, ...) {
  x <- indicators
  dots <- list(...)
  if (x$type != "multi") {
    stop("prec_multi requires a svyplan_n of type 'multi'", call. = FALSE)
  }
  tgt <- x$indicators
  if ("prop_method" %in% names(dots)) {
    tgt$prop_method <- NA_character_
  }
  if ("resp_rate" %in% names(dots)) {
    tgt$resp_rate <- NA_real_
  }
  tgt$n <- x$n
  tgt$moe <- NULL
  tgt$cv <- NULL
  if (!is.null(x$detail) && ".n" %in% names(x$detail)) {
    tgt$n <- x$detail$.n
  }
  dom_cols <- x$params$domain_cols %||% character(0)
  if (!is.null(x$domains) && length(dom_cols) > 0L) {
    dom <- x$domains
    tgt_key <- .domain_key(tgt, dom_cols)
    dom_key <- .domain_key(dom, dom_cols)
    dom_idx <- match(tgt_key, dom_key)
    tgt$n <- dom$.n[dom_idx]
  }
  res <- do.call(prec_multi.default, c(
    list(indicators = tgt, domains = dom_cols), dots
  ))
  # domain_sampling belongs here for the same reason as the rest: n_multi()
  # cannot rebuild the design it describes without it, and reverts to
  # "separate" rather than failing.
  for (p in c("mode", "prop_method", "min_n_domain", "domain_sampling")) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_cluster <- function(indicators, ...) {
  stop(
    "cluster allocations must be passed to prec_multi_cluster()",
    call. = FALSE
  )
}

#' @rdname prec_multi_cluster
#' @export
prec_multi_cluster.svyplan_cluster <- function(indicators, ...) {
  x <- indicators
  dots <- list(...)
  if (is.null(x$indicators)) {
    stop(
      "prec_multi_cluster requires a svyplan_cluster from n_multi_cluster()",
      call. = FALSE
    )
  }
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
  mode <- x$params$mode %||% "cv"
  args <- list(
    indicators = tgt,
    stages = x$stages,
    stage_cost = stage_cost,
    domain_cols = dom_cols,
    mode = mode
  )
  allowed <- c(
    "stage_cost", "domains",
    "resp_rate_psu", "resp_rate_ssu", "resp_rate"
  )
  dot_names <- names(dots) %||% rep("", length(dots))
  unknown <- setdiff(dot_names, allowed)
  if (length(unknown) > 0L || any(!nzchar(dot_names))) {
    .check_unused_dots(...)
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
  for (p in c("budget", "n_psu", "n_per_psu", "n_per_ssu", "joint", "fixed_cost",
              "min_n_domain", "mode")) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res
}
