#' Multi-Indicator Sampling Precision
#'
#' Compute the sampling error (SE, MOE, CV) for multiple survey indicators
#' given a sample size. This is the inverse of [n_multi()].
#'
#' @param targets For the default method: a data frame where **each row
#'   is one survey indicator**, in the same format as [n_multi()] but
#'   with an additional `n` column giving the sample size. This lets you
#'   answer: "given this sample size, what precision do I get for each
#'   indicator?"
#'
#'   At minimum, each row needs:
#'   \itemize{
#'     \item `p` **or** `var`: what you are measuring (see [n_multi()]).
#'     \item `n`: the sample size to evaluate.
#'   }
#'
#'   See the Details section for the full column reference.
#'
#'   For `svyplan_n` or `svyplan_cluster` objects: a result from
#'   [n_multi()].
#' @param ... Additional arguments passed to methods.
#' @param domains Character vector of column names in `targets` to treat
#'   as domain variables, or `NULL` (default) for no domains. All names
#'   must exist in `targets`. Domain columns are preserved in the result
#'   for round-trip conversion back to [n_multi()].
#' @param stage_cost Numeric vector of per-stage costs. `NULL` (default) for
#'   simple mode; length 2 or 3 for multistage mode.
#' @param budget Total budget. Does not affect precision calculations;
#'   preserved for round-trip conversion back to [n_multi()].
#' @param n_psu Fixed stage-1 sample size. Does not affect precision
#'   calculations; preserved for round-trip conversion back to [n_multi()].
#' @param joint Logical. Does not affect precision calculations;
#'   preserved for round-trip conversion back to [n_multi()].
#' @param prop_method Proportion CI method for simple mode, one of `"wald"`
#'   (default), `"wilson"`, or `"logodds"`. This is passed to [prec_prop()]
#'   for proportion rows and ignored for mean rows and multistage mode.
#'   An optional `prop_method` column in `targets` overrides this default
#'   on a per-row basis.
#' @param plan A [svyplan()] profile providing default design parameters.
#'
#' @return A `svyplan_prec` object with a `$detail` data frame containing
#'   per-indicator precision.
#'
#' @details
#' ## Building the targets data frame
#'
#' The `targets` data frame uses the same structure as [n_multi()],
#' with the addition of a required `n` column specifying the sample
#' size to evaluate. A minimal example:
#'
#' ```
#' targets <- data.frame(
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
#'   \item{`var`}{Population variance. One of `p` or `var` per row.}
#'   \item{`mu`}{Population mean. Required for CV output when `var`
#'     is specified, because CV = SE / mean.}
#'   \item{`n`}{Sample size to evaluate (**required**). In simple mode,
#'     one number per row. In multistage mode, this is the stage-1
#'     (PSU) sample size; add `psu_size` and optionally `ssu_size`
#'     columns for the per-stage cluster sizes.}
#'   \item{`psu_size`}{Stage-2 sample size per PSU (multistage only).
#'     Required for 2+ stage designs in multistage mode.}
#'   \item{`ssu_size`}{Stage-3 sample size per SSU (3-stage only).}
#'   \item{`alpha`}{Significance level (default 0.05).}
#'   \item{`deff`}{Design effect multiplier (simple mode only,
#'     default 1).}
#'   \item{`N`}{Population size (simple mode only, default `Inf`).}
#'   \item{`prop_method`}{Proportion CI method: `"wald"` (default),
#'     `"wilson"`, or `"logodds"`. Only for rows with `p` in
#'     simple mode.}
#'   \item{`resp_rate`}{Expected response rate (default 1).}
#'   \item{`delta_psu`, `delta_ssu`}{Homogeneity measures
#'     (multistage).}
#'   \item{`rel_var`}{Unit relvariance. If omitted, derived from `p`
#'     or `var`/`mu`.}
#'   \item{`k_psu`, `k_ssu`}{Ratio parameters (multistage,
#'     default 1).}
#' }
#'
#' Domain columns are specified via the `domains` parameter.
#'
#' In simple mode, `prec_multi()` delegates proportion rows to [prec_prop()]
#' and mean rows to [prec_mean()]. Use `prop_method` or a
#' `targets$prop_method` column to choose `"wald"`, `"wilson"`, or
#' `"logodds"` for proportion rows.
#'
#' @seealso [n_multi()] for the inverse (compute n from precision targets),
#'   [prec_prop()], [prec_mean()] for single-indicator precision.
#'
#' @examples
#' # Simple mode: precision for three indicators at n = 400
#' targets <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   n    = c(400, 400, 400)
#' )
#' prec_multi(targets)
#'
#' # Wilson precision for a rare proportion
#' prec_multi(data.frame(p = 0.05, n = 400), prop_method = "wilson")
#'
#' @export
prec_multi <- function(targets, ...) {
  if (!missing(targets)) {
    .res <- .dispatch_plan(targets, "targets", prec_multi.default, ...)
    if (!is.null(.res)) return(.res)
  }
  UseMethod("prec_multi")
}

#' @rdname prec_multi
#' @export
prec_multi.default <- function(
  targets,
  domains = NULL,
  stage_cost = NULL,
  budget = NULL,
  n_psu = NULL,
  joint = FALSE,
  prop_method = "wald",
  plan = NULL,
  ...
) {
  merged <- .merge_plan_args(
    plan,
    prec_multi.default,
    match.call(),
    environment()
  )
  if (!is.null(merged)) {
    return(do.call(prec_multi.default, merged))
  }

  if (!is.data.frame(targets) || nrow(targets) == 0L) {
    stop("'targets' must be a non-empty data frame", call. = FALSE)
  }

  if (!"n" %in% names(targets)) {
    stop("'targets' must contain an 'n' column for prec_multi", call. = FALSE)
  }
  if (!is.null(domains)) {
    if (!is.character(domains) || anyNA(domains)) {
      stop("'domains' must be a character vector without NAs", call. = FALSE)
    }
    missing_cols <- setdiff(domains, names(targets))
    if (length(missing_cols) > 0L) {
      stop(
        sprintf(
          "domain column(s) not found in targets: %s",
          paste(sQuote(missing_cols), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  if (
    !is.character(prop_method) ||
      length(prop_method) != 1L ||
      is.na(prop_method) ||
      !prop_method %in% c("wald", "wilson", "logodds")
  ) {
    stop(
      "'prop_method' must be one of 'wald', 'wilson', or 'logodds'",
      call. = FALSE
    )
  }

  domain_cols <- domains %||% character(0)
  multistage <- !is.null(stage_cost)

  if (multistage) {
    check_stage_cost(stage_cost)
    stage_cost <- .reorder_stage_cost(stage_cost)
    .prec_multi_cluster(targets, stage_cost, domain_cols = domain_cols)
  } else {
    .prec_multi_simple(targets, prop_method = prop_method,
                       domain_cols = domain_cols)
  }
}

#' @keywords internal
#' @noRd
.prec_multi_simple <- function(targets, prop_method = "wald",
                              domain_cols = character(0)) {
  if (!"alpha" %in% names(targets)) {
    targets$alpha <- 0.05
  }
  if (!"deff" %in% names(targets)) {
    targets$deff <- 1
  }
  if (!"N" %in% names(targets)) {
    targets$N <- Inf
  }
  if (!"resp_rate" %in% names(targets)) {
    targets$resp_rate <- 1
  }
  if (!"prop_method" %in% names(targets)) {
    targets$prop_method <- prop_method
  } else {
    targets$prop_method[is.na(targets$prop_method)] <- prop_method
  }

  .validate_common_columns(targets)

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)
  if (!has_p && !has_var) {
    stop("'targets' must contain 'p' or 'var' column", call. = FALSE)
  }
  has_indicator <- rep(FALSE, nrow(targets))
  if (has_p) {
    has_indicator <- has_indicator | !is.na(targets$p)
  }
  if (has_var) {
    has_indicator <- has_indicator | !is.na(targets$var)
  }
  if (any(!has_indicator)) {
    stop(
      sprintf(
        "row(s) %s must have a non-NA 'p' or 'var' value",
        paste(which(!has_indicator), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (has_p) {
    p_vals <- targets$p[!is.na(targets$p)]
    if (any(p_vals <= 0 | p_vals >= 1)) {
      stop("all 'p' values must be in (0, 1)", call. = FALSE)
    }
  }
  if (has_var) {
    var_vals <- targets$var[!is.na(targets$var)]
    if (any(var_vals <= 0) || any(!is.finite(var_vals))) {
      stop("all 'var' values must be positive and finite", call. = FALSE)
    }
  }
  if ("mu" %in% names(targets)) {
    mu_vals <- targets$mu[!is.na(targets$mu)]
    if (any(mu_vals <= 0) || any(!is.finite(mu_vals))) {
      stop("'mu' values must be positive and finite", call. = FALSE)
    }
  }
  method_vals <- targets$prop_method[!is.na(targets$prop_method)]
  bad_methods <- !method_vals %in% c("wald", "wilson", "logodds")
  if (any(bad_methods)) {
    stop(
      "'prop_method' values must be one of 'wald', 'wilson', or 'logodds'",
      call. = FALSE
    )
  }

  nr <- nrow(targets)
  has_mu <- "mu" %in% names(targets)

  row_labels <- if ("name" %in% names(targets)) {
    targets$name
  } else {
    paste("indicator", seq_len(nr))
  }
  .check_gross_n(targets$n, targets$N, label = row_labels)

  se_vec <- numeric(nr)
  moe_vec <- numeric(nr)
  cv_vec <- numeric(nr)

  for (i in seq_len(nr)) {
    is_prop <- has_p && !is.na(targets$p[i])

    if (is_prop) {
      res_i <- prec_prop.default(
        p = targets$p[i],
        n = targets$n[i],
        alpha = targets$alpha[i],
        N = targets$N[i],
        deff = targets$deff[i],
        resp_rate = targets$resp_rate[i],
        method = targets$prop_method[i]
      )
    } else {
      res_i <- prec_mean.default(
        var = targets$var[i],
        n = targets$n[i],
        mu = if (has_mu && !is.na(targets$mu[i])) targets$mu[i] else NULL,
        alpha = targets$alpha[i],
        N = targets$N[i],
        deff = targets$deff[i],
        resp_rate = targets$resp_rate[i]
      )
    }
    se_vec[i] <- res_i$se
    moe_vec[i] <- res_i$moe
    cv_vec[i] <- res_i$cv
  }

  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  detail <- data.frame(
    name = labels,
    .se = se_vec,
    .moe = moe_vec,
    .cv = cv_vec
  )

  .new_svyplan_prec(
    se = se_vec,
    moe = moe_vec,
    cv = cv_vec,
    type = "multi",
    params = list(targets = targets, domain_cols = domain_cols),
    detail = detail
  )
}

#' @keywords internal
#' @noRd
.prec_multi_cluster <- function(targets, stage_cost,
                                domain_cols = character(0),
                                mode = "cv") {
  stages <- length(stage_cost)

  if (!"alpha" %in% names(targets)) {
    targets$alpha <- 0.05
  }
  if (!"resp_rate" %in% names(targets)) {
    targets$resp_rate <- 1
  }
  if (!"k_psu" %in% names(targets)) {
    targets$k_psu <- 1
  }
  if (!"k_ssu" %in% names(targets)) {
    targets$k_ssu <- 1
  }

  .validate_common_columns(targets)

  if (!"rel_var" %in% names(targets)) {
    targets$rel_var <- NA_real_
  }
  targets$rel_var <- .derive_rel_var(targets, require_all = TRUE)

  rv_check <- targets$rel_var[!is.na(targets$rel_var)]
  if (
    length(rv_check) > 0L &&
      (any(rv_check <= 0) || any(!is.finite(rv_check)))
  ) {
    stop("'rel_var' values must be positive and finite", call. = FALSE)
  }
  if (any(targets$k_psu <= 0) || any(!is.finite(targets$k_psu))) {
    stop("'k_psu' values must be positive and finite", call. = FALSE)
  }
  if (any(targets$k_ssu <= 0) || any(!is.finite(targets$k_ssu))) {
    stop("'k_ssu' values must be positive and finite", call. = FALSE)
  }

  nr <- nrow(targets)
  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  n1 <- targets$n
  if (
    stages >= 2L &&
      (!"psu_size" %in% names(targets) ||
        anyNA(targets$psu_size) ||
        any(targets$psu_size <= 0))
  ) {
    stop(
      "'psu_size' column is required for multistage precision (positive, no NA)",
      call. = FALSE
    )
  }
  if (
    stages == 3L &&
      (!"ssu_size" %in% names(targets) ||
        anyNA(targets$ssu_size) ||
        any(targets$ssu_size <= 0))
  ) {
    stop(
      "'ssu_size' column is required for 3-stage precision (positive, no NA)",
      call. = FALSE
    )
  }
  if (!"delta_psu" %in% names(targets)) {
    stop(
      "'delta_psu' column is required for multistage precision",
      call. = FALSE
    )
  }
  if (!is.numeric(targets$delta_psu)) {
    stop("'delta_psu' must be numeric", call. = FALSE)
  }
  if (anyNA(targets$delta_psu)) {
    stop("'delta_psu' must not contain NA values", call. = FALSE)
  }
  if (any(targets$delta_psu < 0 | targets$delta_psu > 1)) {
    stop("'delta_psu' values must be in [0, 1]", call. = FALSE)
  }
  if (stages == 3L) {
    if (!"delta_ssu" %in% names(targets)) {
      stop(
        "'delta_ssu' column is required for 3-stage precision",
        call. = FALSE
      )
    }
    if (!is.numeric(targets$delta_ssu)) {
      stop("'delta_ssu' must be numeric", call. = FALSE)
    }
    if (anyNA(targets$delta_ssu)) {
      stop("'delta_ssu' must not contain NA values", call. = FALSE)
    }
    if (any(targets$delta_ssu < 0 | targets$delta_ssu > 1)) {
      stop("'delta_ssu' values must be in [0, 1]", call. = FALSE)
    }
  }

  n2 <- targets$psu_size
  n3 <- if (stages == 3L) targets$ssu_size else rep(NA_real_, nr)

  rr <- targets$resp_rate
  n1_eff <- n1 * rr

  cv_vec <- numeric(nr)
  delta1 <- targets$delta_psu
  delta2 <- if ("delta_ssu" %in% names(targets)) {
    targets$delta_ssu
  } else {
    rep(0, nr)
  }
  rel_var <- targets$rel_var
  k1 <- targets$k_psu
  k2 <- targets$k_ssu

  for (i in seq_len(nr)) {
    if (stages == 2L) {
      cv_vec[i] <- sqrt(
        rel_var[i] * k1[i] / (n1_eff[i] * n2[i]) * (1 + delta1[i] * (n2[i] - 1))
      )
    } else {
      cv_vec[i] <- sqrt(
        rel_var[i] /
          (n1_eff[i] * n2[i] * n3[i]) *
          (k1[i] *
            delta1[i] *
            n2[i] *
            n3[i] +
            k2[i] * (1 + delta2[i] * (n3[i] - 1)))
      )
    }
  }

  se_vec <- rep(NA_real_, nr)
  moe_vec <- rep(NA_real_, nr)

  if (identical(mode, "moe")) {
    has_p <- "p" %in% names(targets)
    has_mu <- "mu" %in% names(targets)
    for (i in seq_len(nr)) {
      z <- qnorm(1 - targets$alpha[i] / 2)
      if (has_p && !is.na(targets$p[i])) {
        moe_vec[i] <- cv_vec[i] * z * targets$p[i]
      } else if (has_mu && !is.na(targets$mu[i])) {
        moe_vec[i] <- cv_vec[i] * z * targets$mu[i]
      }
      se_vec[i] <- moe_vec[i] / z
    }
  }

  detail <- data.frame(
    name = labels,
    .se = se_vec,
    .moe = moe_vec,
    .cv = cv_vec
  )

  .new_svyplan_prec(
    se = se_vec,
    moe = moe_vec,
    cv = cv_vec,
    type = "multi",
    params = list(targets = targets, stage_cost = stage_cost,
                  domain_cols = domain_cols),
    detail = detail
  )
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_n <- function(targets, ...) {
  x <- targets
  dots <- list(...)
  if (x$type != "multi") {
    stop("prec_multi requires a svyplan_n of type 'multi'", call. = FALSE)
  }
  tgt <- x$targets
  if ("prop_method" %in% names(dots)) {
    tgt$prop_method <- NA_character_
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
    list(targets = tgt, domains = dom_cols), dots
  ))
  for (p in c("mode", "prop_method", "min_n")) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_cluster <- function(targets, ...) {
  x <- targets
  if (is.null(x$targets)) {
    stop("prec_multi requires a svyplan_cluster from n_multi()", call. = FALSE)
  }
  tgt <- x$targets
  tgt$cv <- NULL
  tgt$moe <- NULL
  tgt$n <- x$n[1L]
  if (x$stages >= 2L) {
    tgt$psu_size <- x$n[2L]
  }
  if (x$stages >= 3L) {
    tgt$ssu_size <- x$n[3L]
  }

  dom_cols <- x$params$domain_cols %||% character(0)
  if (!is.null(x$domains) && length(dom_cols) > 0L) {
    dom <- x$domains
    tgt_key <- .domain_key(tgt, dom_cols)
    dom_key <- .domain_key(dom, dom_cols)
    dom_idx <- match(tgt_key, dom_key)
    tgt$n <- dom$n_psu[dom_idx]
    if (x$stages >= 2L) {
      tgt$psu_size <- dom$psu_size[dom_idx]
    }
    if (x$stages >= 3L) tgt$ssu_size <- dom$ssu_size[dom_idx]
  }

  stage_cost <- x$params$stage_cost
  mode <- x$params$mode %||% "cv"
  res <- .prec_multi_cluster(tgt, stage_cost, domain_cols = dom_cols,
                             mode = mode)
  for (p in c("budget", "n_psu", "psu_size", "ssu_size", "joint", "fixed_cost",
              "min_n", "mode", "prop_method")) {
    if (!is.null(x$params[[p]])) res$params[[p]] <- x$params[[p]]
  }
  res
}
