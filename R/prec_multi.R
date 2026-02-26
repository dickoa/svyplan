#' Multi-Indicator Sampling Precision
#'
#' Compute the sampling error (SE, MOE, CV) for multiple survey indicators
#' given a sample size. This is the inverse of [n_multi()].
#'
#' @param targets For the default method: data frame with one row per
#'   indicator (must contain an `n` column; see Details). For `svyplan_n`
#'   or `svyplan_cluster` objects: a result from [n_multi()].
#' @param ... Additional arguments passed to methods.
#' @param cost Numeric vector of per-stage costs. `NULL` (default) for
#'   simple mode; length 2 or 3 for multistage mode.
#' @param budget Total budget (currently unused in precision mode).
#' @param m Fixed stage-1 sample size (currently unused in precision mode).
#' @param joint Logical (currently unused in precision mode).
#'
#' @return A `svyplan_prec` object with a `$detail` data frame containing
#'   per-indicator precision.
#'
#' @details
#' The `targets` data frame supports the following columns:
#'
#' \describe{
#'   \item{`name`}{Indicator label (optional).}
#'   \item{`p`}{Expected proportion, in (0, 1). One of `p` or `var` per row.}
#'   \item{`var`}{Population variance. One of `p` or `var` per row.}
#'   \item{`mu`}{Population mean. Required for CV when `var` is specified.}
#'   \item{`n`}{Sample size (required). For simple mode, a scalar per
#'     indicator. For multistage, per-stage sizes can be provided as
#'     `n`, `n2`, `n3` columns.}
#'   \item{`alpha`}{Significance level (default 0.05).}
#'   \item{`deff`}{Design effect multiplier (simple mode only, default 1).}
#'   \item{`N`}{Population size (simple mode only, default Inf).}
#'   \item{`resp_rate`}{Expected response rate (default 1).}
#'   \item{`delta1`, `delta2`}{Homogeneity measures (multistage).}
#'   \item{`rel_var`}{Unit relvariance. If omitted, derived from `p` or
#'     `var`/`mu`.}
#'   \item{`k1`, `k2`}{Ratio parameters (multistage, default 1).}
#' }
#'
#' Any column not in the recognized set is treated as a **domain variable**.
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
#' @export
prec_multi <- function(targets, ...) UseMethod("prec_multi")

#' @rdname prec_multi
#' @export
prec_multi.default <- function(targets, cost = NULL, budget = NULL, m = NULL,
                               joint = FALSE, ...) {
  if (!is.data.frame(targets) || nrow(targets) == 0L) {
    stop("'targets' must be a non-empty data frame", call. = FALSE)
  }

  if (!"n" %in% names(targets)) {
    stop("'targets' must contain an 'n' column for prec_multi", call. = FALSE)
  }

  multistage <- !is.null(cost)

  if (multistage) {
    check_cost(cost)
    .prec_multi_cluster(targets, cost)
  } else {
    .prec_multi_simple(targets)
  }
}

#' @keywords internal
#' @noRd
.prec_multi_simple <- function(targets) {
  if (!"alpha" %in% names(targets)) targets$alpha <- 0.05
  if (!"deff" %in% names(targets)) targets$deff <- 1
  if (!"N" %in% names(targets)) targets$N <- Inf
  if (!"resp_rate" %in% names(targets)) targets$resp_rate <- 1

  .validate_common_columns(targets)

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)
  if (!has_p && !has_var)
    stop("'targets' must contain 'p' or 'var' column", call. = FALSE)
  has_indicator <- rep(FALSE, nrow(targets))
  if (has_p) has_indicator <- has_indicator | !is.na(targets$p)
  if (has_var) has_indicator <- has_indicator | !is.na(targets$var)
  if (any(!has_indicator))
    stop(sprintf("row(s) %s must have a non-NA 'p' or 'var' value",
         paste(which(!has_indicator), collapse = ", ")), call. = FALSE)

  if (has_p) {
    p_vals <- targets$p[!is.na(targets$p)]
    if (any(p_vals <= 0 | p_vals >= 1))
      stop("all 'p' values must be in (0, 1)", call. = FALSE)
  }
  if (has_var) {
    var_vals <- targets$var[!is.na(targets$var)]
    if (any(var_vals <= 0))
      stop("all 'var' values must be positive", call. = FALSE)
  }
  if ("mu" %in% names(targets)) {
    mu_vals <- targets$mu[!is.na(targets$mu)]
    if (any(mu_vals <= 0))
      stop("'mu' values must be positive", call. = FALSE)
  }

  nr <- nrow(targets)
  has_mu <- "mu" %in% names(targets)

  se_vec <- numeric(nr)
  moe_vec <- numeric(nr)
  cv_vec <- numeric(nr)

  for (i in seq_len(nr)) {
    z <- qnorm(1 - targets$alpha[i] / 2)
    n_eff <- targets$n[i] * targets$resp_rate[i] / targets$deff[i]
    N <- targets$N[i]

    is_prop <- has_p && !is.na(targets$p[i])

    if (is_prop) {
      p <- targets$p[i]
      q <- 1 - p
      fpc <- if (is.infinite(N)) 1 else (N - n_eff) / (N - 1)
      fpc <- .clamp_fpc(fpc, n_eff, N)
      se_vec[i] <- sqrt(p * q * fpc / n_eff)
      moe_vec[i] <- z * se_vec[i]
      cv_vec[i] <- se_vec[i] / p
    } else {
      v <- targets$var[i]
      fpc <- if (is.infinite(N)) 1 else 1 - n_eff / N
      fpc <- .clamp_fpc(fpc, n_eff, N)
      se_vec[i] <- sqrt(v * fpc / n_eff)
      moe_vec[i] <- z * se_vec[i]
      if (has_mu && !is.na(targets$mu[i])) {
        cv_vec[i] <- se_vec[i] / targets$mu[i]
      } else {
        cv_vec[i] <- NA_real_
      }
    }
  }

  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  detail <- data.frame(
    name = labels,
    .se = se_vec,
    .moe = moe_vec,
    .cv = cv_vec
  )

  .new_svyplan_prec(
    se     = se_vec,
    moe    = moe_vec,
    cv     = cv_vec,
    type   = "multi",
    params = list(targets = targets),
    detail = detail
  )
}

#' @keywords internal
#' @noRd
.prec_multi_cluster <- function(targets, cost) {
  stages <- length(cost)

  if (!"alpha" %in% names(targets)) targets$alpha <- 0.05
  if (!"resp_rate" %in% names(targets)) targets$resp_rate <- 1
  if (!"k1" %in% names(targets)) targets$k1 <- 1
  if (!"k2" %in% names(targets)) targets$k2 <- 1

  .validate_common_columns(targets)

  if (!"rel_var" %in% names(targets)) targets$rel_var <- NA_real_
  targets$rel_var <- .derive_rel_var(targets, require_all = TRUE)

  nr <- nrow(targets)
  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  n1 <- targets$n
  if (stages >= 2L && (!"n2" %in% names(targets) || anyNA(targets$n2) || any(targets$n2 <= 0)))
    stop("'n2' column is required for multistage precision (positive, no NA)", call. = FALSE)
  if (stages == 3L && (!"n3" %in% names(targets) || anyNA(targets$n3) || any(targets$n3 <= 0)))
    stop("'n3' column is required for 3-stage precision (positive, no NA)", call. = FALSE)
  n2 <- targets$n2
  n3 <- if (stages == 3L) targets$n3 else rep(NA_real_, nr)

  rr <- targets$resp_rate
  n1_eff <- n1 * rr

  cv_vec <- numeric(nr)
  delta1 <- targets$delta1
  delta2 <- if ("delta2" %in% names(targets)) targets$delta2 else rep(0, nr)
  rel_var <- targets$rel_var
  k1 <- targets$k1
  k2 <- targets$k2

  for (i in seq_len(nr)) {
    if (stages == 2L) {
      cv_vec[i] <- sqrt(
        rel_var[i] * k1[i] / (n1_eff[i] * n2[i]) *
          (1 + delta1[i] * (n2[i] - 1))
      )
    } else {
      cv_vec[i] <- sqrt(
        rel_var[i] / (n1_eff[i] * n2[i] * n3[i]) *
          (k1[i] * delta1[i] * n2[i] * n3[i] +
             k2[i] * (1 + delta2[i] * (n3[i] - 1)))
      )
    }
  }

  detail <- data.frame(
    name = labels,
    .cv = cv_vec
  )

  .new_svyplan_prec(
    se     = rep(NA_real_, nr),
    moe    = rep(NA_real_, nr),
    cv     = cv_vec,
    type   = "multi",
    params = list(targets = targets, cost = cost),
    detail = detail
  )
}

#' @rdname prec_multi
#' @export
prec_multi.svyplan_n <- function(targets, ...) {
  x <- targets
  if (x$type != "multi") {
    stop("prec_multi requires a svyplan_n of type 'multi'", call. = FALSE)
  }
  tgt <- x$targets
  tgt$n <- x$n
  tgt$moe <- NULL
  tgt$cv <- NULL
  if (!is.null(x$detail) && ".n" %in% names(x$detail)) {
    tgt$n <- x$detail$.n
  }
  if (!is.null(x$domains)) {
    dom <- x$domains
    known_dom <- c(".n", ".binding")
    dom_cols <- setdiff(names(dom), known_dom)
    tgt_key <- interaction(tgt[dom_cols], drop = FALSE, sep = ":")
    dom_key <- interaction(dom[dom_cols], drop = FALSE, sep = ":")
    dom_idx <- match(tgt_key, dom_key)
    tgt$n <- dom$.n[dom_idx]
  }
  prec_multi.default(targets = tgt)
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
  if (x$stages >= 2L) tgt$n2 <- x$n[2L]
  if (x$stages >= 3L) tgt$n3 <- x$n[3L]

  if (!is.null(x$domains)) {
    dom <- x$domains
    stage_cols <- paste0("n", seq_len(x$stages))
    known_dom <- c(stage_cols, ".total_n", ".cv", ".cost", ".binding")
    dom_cols <- setdiff(names(dom), known_dom)
    tgt_key <- interaction(tgt[dom_cols], drop = FALSE, sep = ":")
    dom_key <- interaction(dom[dom_cols], drop = FALSE, sep = ":")
    dom_idx <- match(tgt_key, dom_key)
    tgt$n <- dom$n1[dom_idx]
    if (x$stages >= 2L) tgt$n2 <- dom$n2[dom_idx]
    if (x$stages >= 3L) tgt$n3 <- dom$n3[dom_idx]
  }

  cost <- x$params$cost
  prec_multi.default(targets = tgt, cost = cost)
}
