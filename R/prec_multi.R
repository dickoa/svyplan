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
#'   For `svyplan_n` objects: a result from [n_multi()]. Every row is
#'   assessed at the size the result chose, the common `n`, a domain's
#'   quota, or under natural domain sampling its expected yield
#'   `n * share`. Rows with only a floor are dropped, since they have no
#'   precision to report.
#' @param ... Additional arguments passed to methods. Unused arguments are rejected.
#' @param domains Character vector of column names in `indicators` to treat
#'   as domain variables, or `NULL` (default) for no domains. All names
#'   must exist in `indicators`. Domain columns are preserved in the result
#'   for round-trip conversion back to [n_multi()], so none can be called
#'   `name` or start with a dot.
#' @param prop_method Proportion CI method, one of `"wald"`
#'   (default), `"wilson"`, `"logodds"`, or `"beta"`. This is passed to
#'   [prec_prop()] for proportion rows and ignored for mean rows.
#'   An optional `prop_method` column in `indicators` overrides this default
#'   on a per-row basis.
#' @param resp_rate Default expected response rate among eligible ultimate
#'   units, in (0, 1\]. Used for rows whose `resp_rate` column is absent or
#'   `NA`, and a non-missing row value overrides it. See the nonresponse
#'   section of [svyplan-package] for what the rate counts.
#' @param eligibility_rate Default share of issued units in scope, in (0, 1\].
#'   Used for rows whose `eligibility_rate` column is absent or `NA`, and a
#'   non-missing row value overrides it. See the nonresponse section of
#'   [svyplan-package].
#' @param eligible_per_unit Default eligible units per sampled unit, a positive
#'   number that can exceed 1. Used for rows whose `eligible_per_unit` column is
#'   absent or `NA`. Proportion and mean rows. See [n_multi()].
#' @param plan A [svyplan()] profile providing default design parameters.
#'
#' @return A `svyplan_prec` object with a `$detail` data frame containing
#'   per-indicator precision: `.se`, `.moe`, `.rmoe`, and `.cv`. `.rmoe`
#'   is `.moe` as a fraction of the row's own estimand, `p` for a
#'   proportion, `abs(mean)` for a mean and `abs(r)` for a ratio, and `NA` for
#'   a row carrying none of them. The `domains` columns come first, so each
#'   row can be joined back to its domain.
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
#'   \item{`p`}{Expected proportion, in (0, 1). One of `p`, `var` or `r`
#'     per row (see [n_multi()]).}
#'   \item{`var`}{Population variance. One of `p`, `var` or `r` per row.}
#'   \item{`mean`}{Population mean. Required for CV output when `var`
#'     is specified, because CV = SE / abs(mean).}
#'   \item{`r`}{Anticipated ratio of two totals, marking a ratio row. It
#'     requires `cv_num`, `cv_den` and `component_cor` alongside it. Its
#'     `unit_relvar` comes from those moments, and a value given in that
#'     column must agree with them. See [prec_ratio()].}
#'   \item{`cv_num`, `cv_den`}{Coefficients of variation of a ratio row's
#'     numerator and denominator, strictly positive.}
#'   \item{`component_cor`}{Correlation between a ratio row's numerator and
#'     denominator across units, in \eqn{[-1, 1]}.}
#'   \item{`n`}{Sample size to evaluate (**required**), with one value per
#'     indicator row.}
#'   \item{`alpha`}{Significance level (default 0.05). `2 * pnorm(-2)`
#'     gives the multiplier of 2 that the DHS and MICS manuals use in place
#'     of 1.96.}
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
#'     (default 1). A cluster design spends a rate at stage 1 too, and
#'     [prec_cluster()] names stage response columns `resp_rate_1`, ..., `resp_rate_{S-1}`.}
#'   \item{`eligibility_rate`}{Share of the row's issued units in scope, in
#'     (0, 1\]. `NA` takes the `eligibility_rate` argument.}
#'   \item{`eligible_per_unit`}{Eligible units per sampled unit for the row's
#'     estimate, as in [prec_prop()] and [prec_mean()]. Proportion and mean
#'     rows. The
#'     several-indicators mode of [prec_cluster()] reads the same column.}
#' }
#'
#' Domain columns are specified via the `domains` parameter.
#'
#' `prec_multi()` delegates proportion rows to [prec_prop()], mean rows to
#' [prec_mean()] and ratio rows to [prec_ratio()]. Use `prop_method` or a
#' `indicators$prop_method` column to choose `"wald"`, `"wilson"`,
#' `"logodds"` or `"beta"` for proportion rows. It is ignored elsewhere.
#'
#' ## A domain that takes most of its PSUs
#'
#' A `deff` column carries the clustering of a design whose first stage is
#' treated as drawn with replacement, and `N` corrects only for the units.
#' Nothing here knows how many PSUs a domain holds. A small domain, such as
#' an island district of 16 enumeration areas, can then be given more PSUs
#' than it has, and a sample that takes most of them is assessed as if
#' they were a few among many. [n_alloc()] and [prec_alloc()] read `N_1`
#' in the frame and, with `fpc = "stage"`, apply the stage-by-stage
#' correction. They size and assess such a domain, as the Details of
#' [n_alloc()] describe.
#'
#' @family multi-indicator functions
#' @seealso [n_multi()] for the inverse, [prec_cluster()] for a multistage
#'   cluster design over the same indicator table, [prec_alloc()] for a
#'   stratified allocation over the whole population and its domains in one
#'   call, and [prec_prop()] and [prec_mean()] for single-indicator
#'   precision.
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
  eligibility_rate = 1,
  eligible_per_unit = 1,
  plan = NULL
) {
  .refuse_alloc_object(indicators, "prec_multi", "multi")
  merged <- .merge_plan_args(
    plan,
    prec_multi.default,
    match.call(),
    environment()
  )
  if (!is.null(merged)) {
    return(do.call(prec_multi.default, c(merged, list(...))))
  }
  .check_multi_split_args(list(...), "prec_cluster(<table>, ...)")
  .check_unused_dots(...)
  .refuse_unnumbered_stage_columns(indicators, "indicators")
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
  check_resp_rate(eligibility_rate, "eligibility_rate")
  check_eligible_per_unit(eligible_per_unit)

  if (!is.data.frame(indicators) || nrow(indicators) == 0L) {
    stop("'indicators' must be a non-empty data frame", call. = FALSE)
  }

  indicators <- .indicators_var_from_sd(indicators, domains)
  .check_indicator_columns(indicators, domains, eligible = TRUE,
                           context = "prec_multi()")
  .check_resp_rate_column(indicators, FALSE)
  .stop_min_cases_column(
    indicators,
    "prec_multi() reads the sizes you already have. A floor is a sizing constraint, so it belongs to n_multi()"
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
                     eligibility_rate = eligibility_rate,
                     eligible_per_unit = eligible_per_unit,
                     domain_cols = domain_cols)
}

#' Put the domain columns of an indicator table before its precision rows
#'
#' Without them a row of `$detail` can be matched back to its domain only by
#' position, which a reordered table breaks.
#' @keywords internal
#' @noRd
.with_domain_keys <- function(detail, indicators, domain_cols) {
  if (length(domain_cols) == 0L) {
    return(detail)
  }
  keys <- as.data.frame(indicators)[, domain_cols, drop = FALSE]
  rownames(keys) <- NULL
  cbind(keys, detail)
}

#' @keywords internal
#' @noRd
.prec_multi_simple <- function(indicators, prop_method = "wald",
                              resp_rate = 1,
                              eligibility_rate = 1,
                              eligible_per_unit = 1,
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
  indicators <- .fill_eligible_per_unit(indicators, eligible_per_unit)
  indicators <- .fill_eligibility_rate(indicators, eligibility_rate)

  .validate_common_columns(indicators)

  has_p <- "p" %in% names(indicators)
  has_var <- "var" %in% names(indicators)
  is_ratio <- .is_ratio_row(indicators)
  .check_ratio_rows(indicators, sizing = FALSE)
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
  if ("mean" %in% names(indicators)) {
    mean_vals <- indicators$mean[!is.na(indicators$mean)]
    if (any(!is.finite(mean_vals))) {
      stop("'mean' values must be finite", call. = FALSE)
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
  has_mean <- "mean" %in% names(indicators)

  row_labels <- if ("name" %in% names(indicators)) {
    indicators$name
  } else {
    paste("indicator", seq_len(nr))
  }
  .check_gross_n(indicators$n, indicators$N, label = row_labels)
  elig <- .indicator_rate(indicators, "eligible_per_unit", nr)
  scope <- .indicator_rate(indicators, "eligibility_rate", nr)

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
        eligibility_rate = scope[i],
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
        eligibility_rate = scope[i],
        eligible_per_unit = elig[i],
        method = indicators$prop_method[i],
        df = .row_df(indicators, i)
      )
    } else {
      res_i <- prec_mean.default(
        var = indicators$var[i],
        n = indicators$n[i],
        mu = if (has_mean && !is.na(indicators$mean[i])) indicators$mean[i] else NULL,
        alpha = indicators$alpha[i],
        N = indicators$N[i],
        deff = indicators$deff[i],
        resp_rate = indicators$resp_rate[i],
        eligibility_rate = scope[i],
        eligible_per_unit = elig[i],
        df = .row_df(indicators, i)
      )
    }
    se_vec[i] <- res_i$se
    moe_vec[i] <- res_i$moe
    cv_vec[i] <- res_i$cv
  }

  labels <- if ("name" %in% names(indicators)) indicators$name else seq_len(nr)

  detail <- .with_domain_keys(
    data.frame(
      name = labels,
      .se = se_vec,
      .moe = moe_vec,
      .rmoe = .rmoe_from_moe(moe_vec, .indicator_estimand(indicators, nr),
                             at_zero = NA_real_),
      .cv = cv_vec
    ),
    indicators,
    domain_cols
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

#' @rdname prec_multi
#' @export
prec_multi.svyplan_multi <- function(indicators, ...) {
  .prec_multi_from_fit(indicators, list(...))
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_n <- function(indicators, ...) {
  .prec_multi_from_fit(indicators, list(...))
}

#' Read a multi-indicator size back as precision
#'
#' At the size the fit chose, or with `field = TRUE` at the size it fields:
#' the whole total, each separate quota ceiled, or the whole total times
#' each domain's share under natural incidence.
#' @keywords internal
#' @noRd
.prec_multi_from_fit <- function(x, dots, field = FALSE) {
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
  if ("eligible_per_unit" %in% names(dots) && "eligible_per_unit" %in% names(tgt)) {
    tgt$eligible_per_unit <- NA_real_
  }
  if ("eligibility_rate" %in% names(dots) && "eligibility_rate" %in% names(tgt)) {
    tgt$eligibility_rate <- NA_real_
  }
  # The design is assessed as sized: every row at the size the fit chose,
  # not at its own requirement, which is the target it was sized to. Floors
  # shaped that size and have no precision of their own to report.
  tgt <- tgt[!.floor_only_rows(tgt), , drop = FALSE]
  tgt$min_cases <- NULL
  tgt$min_respondents <- NULL
  total <- if (field) .whole_count(x$n) else x$n
  tgt$n <- total
  tgt$moe <- NULL
  tgt$cv <- NULL
  dom_cols <- x$params$domain_cols %||% character(0)
  if (!is.null(x$domains) && length(dom_cols) > 0L) {
    dom <- x$domains
    tgt_key <- .domain_key(tgt, dom_cols)
    dom_key <- .domain_key(dom, dom_cols)
    dom_idx <- match(tgt_key, dom_key)
    tgt$n <- if (".share" %in% names(dom)) {
      total * dom$.share[dom_idx]
    } else if (field) {
      .whole_count(dom$.n[dom_idx])
    } else {
      dom$.n[dom_idx]
    }
  }
  res <- do.call(prec_multi.default, c(
    list(indicators = tgt, domains = dom_cols), dots
  ))
  # domain_sampling belongs here for the same reason as the rest: n_multi()
  # cannot rebuild the design it describes without it, and reverts to
  # "separate" rather than failing.
  for (p in setdiff(c("mode", "prop_method", "min_n_domain", "domain_sampling"),
                    names(dots))) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_cluster <- function(indicators, ...) {
  stop(
    "cluster allocations must be passed to prec_cluster()",
    call. = FALSE
  )
}
