#' Multi-Indicator Sample Size
#'
#' Compute the sample size that satisfies precision requirements for
#' multiple survey indicators simultaneously. Supports simple (single-stage)
#' and multistage cluster designs, with optional domain-level planning.
#'
#' @param targets For the default method: data frame with one row per
#'   indicator (see Details). For `svyplan_prec` objects: a precision
#'   result from [prec_multi()].
#' @param ... Additional arguments passed to methods.
#' @param cost Numeric vector of per-stage costs. `NULL` (default) for
#'   simple mode; length 2 or 3 for multistage mode.
#' @param budget Total budget (multistage only). Provide either `cv` values
#'   in the `targets` data frame or a `budget` here, not both.
#' @param n_psu Fixed stage-1 sample size (multistage only).
#' @param joint Logical. If `TRUE`, optimally split a single `budget`
#'   across domains to minimize the worst-case CV ratio. Only applies
#'   to multistage budget mode with multiple domains; ignored otherwise.
#' @param min_n Numeric scalar or `NULL` (default). Minimum total sample
#'   size per domain. Only active when domains are present; silently
#'   ignored otherwise. In simple mode, per-domain sample sizes are
#'   floored to `min_n`. In joint multistage mode, domains that would
#'   receive fewer than `min_n` observations are penalized during
#'   optimization, with an upfront feasibility check. In non-joint
#'   multistage mode, a warning is issued for any domain below the floor.
#' @param fixed_cost Fixed overhead cost (C0). Default 0. Only applies
#'   to multistage mode. See [n_cluster()] for details.
#'
#' @return A `svyplan_n` object (simple mode) or `svyplan_cluster` object
#'   (multistage mode).
#'
#'   **Without domains**, the object contains:
#'   \describe{
#'     \item{`n`}{Sample size (simple) or named per-stage allocation vector
#'       (multistage, e.g. `c(n_psu = 80, psu_size = 12)`).}
#'     \item{`detail`}{Per-indicator results (sample sizes or achieved CVs).}
#'     \item{`binding`}{Name or index of the binding (most demanding) indicator.}
#'     \item{`targets`}{The input targets data frame.}
#'   }
#'
#'   **With domains**, the object additionally contains:
#'   \describe{
#'     \item{`n`}{Maximum per-stage sample size across domains. In simple
#'       mode, a single number; in multistage mode, a named vector
#'       (e.g. `c(n_psu = 120, psu_size = 15)`) giving the conservative
#'       allocation that satisfies all domains.}
#'     \item{`domains`}{Data frame with one row per domain, including
#'       domain variable columns, per-stage allocations (`n_psu`, `psu_size`, ...),
#'       and summary columns (`.total_n`, `.cv`, `.cost`, `.binding`
#'       for multistage; `.n`, `.binding` for simple mode).
#'       Use this for stratum-specific allocations.}
#'     \item{`total_n`}{Total sample size summed across all domains
#'       (multistage only).}
#'     \item{`cost`}{Total cost summed across all domains
#'       (multistage only).}
#'   }
#'
#' @details
#' The `targets` data frame supports the following columns:
#'
#' \describe{
#'   \item{`name`}{Indicator label (optional).}
#'   \item{`p`}{Expected proportion, in (0, 1). One of `p` or `var` per row.}
#'   \item{`var`}{Population variance. One of `p` or `var` per row.}
#'   \item{`mu`}{Population mean magnitude (positive). Required when `var` is specified with `cv`.}
#'   \item{`moe`}{Margin of error (simple mode).}
#'   \item{`cv`}{Target coefficient of variation (either mode).}
#'   \item{`alpha`}{Significance level (default 0.05).}
#'   \item{`deff`}{Design effect multiplier (simple mode only, default 1).}
#'   \item{`N`}{Population size (simple mode only, default Inf).}
#'   \item{`delta_psu`, `delta_ssu`}{Homogeneity measures (multistage).}
#'   \item{`rel_var`}{Unit relvariance. If omitted, derived from `p` or
#'     `var`/`mu`.}
#'   \item{`k_psu`, `k_ssu`}{Ratio parameters (multistage, default 1).}
#'   \item{`resp_rate`}{Expected response rate, in (0, 1\]. Default 1 (no
#'     adjustment). Inflates the required sample size to account for
#'     non-response.}
#' }
#'
#' Any column not in the recognized set is treated as a **domain variable**.
#' When domain columns are present, optimization runs independently per
#' domain combination (default), or jointly when `joint = TRUE`.
#'
#' **Simple mode** (`cost = NULL`): computes sample size per indicator using
#' Wald-type formulas, then takes the maximum per domain.
#'
#' **Multistage mode** (`cost` provided): uses analytical reduction.
#' For each candidate sub-stage allocation, the required stage-1 size is
#' the maximum across all indicators. The total cost is then minimized
#' (CV mode) or the worst-case CV ratio is minimized (budget mode) using
#' numerical optimization.
#'
#' **Joint budget allocation** (`joint = TRUE`): when domains and a budget
#' are specified, the default (`joint = FALSE`) gives each domain the full
#' budget independently. With `joint = TRUE`, a single budget is split
#' optimally across domains using L-BFGS-B optimization of budget fractions,
#' minimizing the worst-case CV ratio across all domains.
#'
#' These functions assume sampling fractions are negligible at each stage
#' (equivalent to sampling with replacement). No finite population correction
#' is applied. This is standard for multistage planning when cluster
#' populations are large relative to the sample.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @seealso [n_prop()], [n_mean()] for single-indicator sizing;
#'   [n_cluster()] for single-indicator multistage allocation;
#'   [prec_multi()] for the inverse.
#'
#' @examples
#' # Simple mode: three indicators, take the max
#' targets <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   moe  = c(0.05, 0.05, 0.03)
#' )
#' n_multi(targets)
#'
#' # Simple mode with domains
#' targets_dom <- data.frame(
#'   name   = rep(c("stunting", "anemia"), each = 2),
#'   p      = c(0.30, 0.25, 0.10, 0.15),
#'   moe    = c(0.05, 0.05, 0.03, 0.03),
#'   region = rep(c("North", "South"), 2)
#' )
#' n_multi(targets_dom)
#'
#' # Two-stage CV mode
#' targets_cl <- data.frame(
#'   name   = c("stunting", "anemia"),
#'   p      = c(0.30, 0.10),
#'   cv     = c(0.10, 0.15),
#'   delta_psu = c(0.02, 0.05)
#' )
#' n_multi(targets_cl, cost = c(500, 50))
#'
#' # Joint budget allocation across domains
#' targets_jnt <- data.frame(
#'   name   = rep(c("stunting", "anemia"), each = 2),
#'   p      = c(0.30, 0.25, 0.10, 0.15),
#'   cv     = c(0.10, 0.10, 0.15, 0.15),
#'   delta_psu = c(0.02, 0.03, 0.05, 0.04),
#'   region = rep(c("Urban", "Rural"), 2)
#' )
#' n_multi(targets_jnt, cost = c(500, 50), budget = 100000, joint = TRUE)
#'
#' @export
n_multi <- function(targets, ...) {
  UseMethod("n_multi")
}

#' @rdname n_multi
#' @export
n_multi.default <- function(
  targets,
  cost = NULL,
  budget = NULL,
  n_psu = NULL,
  joint = FALSE,
  min_n = NULL,
  fixed_cost = 0,
  ...
) {
  if (!is.data.frame(targets) || nrow(targets) == 0L) {
    stop("'targets' must be a non-empty data frame", call. = FALSE)
  }
  if (!is.logical(joint) || length(joint) != 1L || is.na(joint)) {
    stop("'joint' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(min_n)) {
    if (
      !is.numeric(min_n) || length(min_n) != 1L || is.na(min_n) || min_n <= 0
    ) {
      stop("'min_n' must be a positive numeric scalar", call. = FALSE)
    }
  }

  multistage <- !is.null(cost)

  if (multistage) {
    check_cost(cost)
    if (!is.null(budget)) {
      check_scalar(budget, "budget")
    }
    if (!is.null(n_psu)) {
      check_scalar(n_psu, "n_psu")
    }
    check_fixed_cost(fixed_cost, budget)
  } else {
    if (!is.null(budget)) {
      stop("'budget' requires 'cost' to be specified", call. = FALSE)
    }
    if (!is.null(n_psu)) {
      stop("'n_psu' requires 'cost' to be specified", call. = FALSE)
    }
  }

  info <- .validate_targets(targets, multistage)
  targets <- .fill_defaults(targets, multistage)

  rv_final <- targets$rel_var[!is.na(targets$rel_var)]
  if (length(rv_final) > 0L &&
    (any(rv_final <= 0) || any(!is.finite(rv_final)))) {
    stop("'rel_var' values must be positive and finite", call. = FALSE)
  }
  if (multistage) {
    if (any(targets$k_psu <= 0) || any(!is.finite(targets$k_psu))) {
      stop("'k_psu' values must be positive and finite", call. = FALSE)
    }
    if (any(targets$k_ssu <= 0) || any(!is.finite(targets$k_ssu))) {
      stop("'k_ssu' values must be positive and finite", call. = FALSE)
    }
  }

  domain_cols <- info$domain_cols

  if (length(domain_cols) == 0L) {
    if (multistage) {
      .n_multi_cluster(targets, cost, budget, n_psu, fixed_cost)
    } else {
      .n_multi_simple(targets)
    }
  } else {
    .n_multi_domains(
      targets,
      cost,
      budget,
      n_psu,
      domain_cols,
      multistage,
      joint,
      min_n,
      fixed_cost
    )
  }
}

#' Validate targets data frame
#' @return List with `indicator_type` (per-row "p" or "var") and `domain_cols`.
#' @keywords internal
#' @noRd
.validate_targets <- function(targets, multistage) {
  known <- c(
    "name",
    "p",
    "var",
    "mu",
    "moe",
    "cv",
    "alpha",
    "deff",
    "N",
    "delta_psu",
    "delta_ssu",
    "rel_var",
    "k_psu",
    "k_ssu",
    "resp_rate"
  )

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)
  if (!has_p && !has_var) {
    stop("'targets' must contain 'p' or 'var' column", call. = FALSE)
  }

  has_moe <- "moe" %in% names(targets)
  has_cv <- "cv" %in% names(targets)
  if (!has_moe && !has_cv) {
    stop("'targets' must contain 'moe' or 'cv' column", call. = FALSE)
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

  # Each row needs at least one of p or var (non-NA)
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

  if (has_p && has_var) {
    both_set <- !is.na(targets$p) & !is.na(targets$var)
    if (any(both_set)) {
      stop("each row must have only one of 'p' or 'var'", call. = FALSE)
    }
  }

  if (has_moe && has_cv) {
    both_na <- is.na(targets$moe) & is.na(targets$cv)
    both_set <- !is.na(targets$moe) & !is.na(targets$cv)
    if (any(both_na)) {
      stop("each row must have either 'moe' or 'cv' specified", call. = FALSE)
    }
    if (any(both_set)) {
      stop("each row must have only one of 'moe' or 'cv'", call. = FALSE)
    }
  }

  if (has_moe) {
    moe_vals <- targets$moe[!is.na(targets$moe)]
    if (any(moe_vals <= 0) || any(!is.finite(moe_vals))) {
      stop("'moe' values must be positive and finite", call. = FALSE)
    }
  }

  if (multistage) {
    if (!has_cv) {
      stop("multistage mode requires 'cv' column in targets", call. = FALSE)
    }
    if (has_moe) {
      if (any(!is.na(targets$moe))) {
        stop("multistage mode requires 'cv' (not 'moe')", call. = FALSE)
      }
    }
    if (!"delta_psu" %in% names(targets)) {
      stop("multistage mode requires 'delta_psu' column in targets", call. = FALSE)
    }
    cv_vals <- targets$cv[!is.na(targets$cv)]
    if (any(cv_vals <= 0) || any(!is.finite(cv_vals))) {
      stop("'cv' values must be positive and finite", call. = FALSE)
    }
    d1_vals <- targets$delta_psu[!is.na(targets$delta_psu)]
    if (any(d1_vals < 0 | d1_vals > 1)) {
      stop("'delta_psu' values must be in [0, 1]", call. = FALSE)
    }
    if ("delta_ssu" %in% names(targets)) {
      d2_vals <- targets$delta_ssu[!is.na(targets$delta_ssu)]
      if (any(d2_vals < 0 | d2_vals > 1)) {
        stop("'delta_ssu' values must be in [0, 1]", call. = FALSE)
      }
    }
  }

  if (has_var && has_cv) {
    needs_mu <- !is.na(targets$var) & !is.na(targets$cv)
    if (any(needs_mu)) {
      if (!"mu" %in% names(targets) || any(is.na(targets$mu[needs_mu]))) {
        stop(
          "'mu' is required when 'var' and 'cv' are specified",
          call. = FALSE
        )
      }
    }
  }

  if ("mu" %in% names(targets)) {
    mu_vals <- targets$mu[!is.na(targets$mu)]
    if (any(mu_vals <= 0) || any(!is.finite(mu_vals))) {
      stop("'mu' values must be positive and finite", call. = FALSE)
    }
  }

  domain_cols <- setdiff(names(targets), known)

  if (length(domain_cols) > 0L) {
    message(
      "Treating column(s) ",
      paste(sQuote(domain_cols), collapse = ", "),
      " as domain variable(s)"
    )
  }

  if ("rel_var" %in% names(targets)) {
    rv_vals <- targets$rel_var[!is.na(targets$rel_var)]
    if (any(rv_vals <= 0) || any(!is.finite(rv_vals))) {
      stop("'rel_var' values must be positive and finite", call. = FALSE)
    }
  }
  if (multistage) {
    if ("k_psu" %in% names(targets)) {
      k_psu_vals <- targets$k_psu[!is.na(targets$k_psu)]
      if (any(k_psu_vals <= 0) || any(!is.finite(k_psu_vals))) {
        stop("'k_psu' values must be positive and finite", call. = FALSE)
      }
    }
    if ("k_ssu" %in% names(targets)) {
      k_ssu_vals <- targets$k_ssu[!is.na(targets$k_ssu)]
      if (any(k_ssu_vals <= 0) || any(!is.finite(k_ssu_vals))) {
        stop("'k_ssu' values must be positive and finite", call. = FALSE)
      }
    }
  }

  .validate_common_columns(targets)

  list(domain_cols = domain_cols)
}

#' Fill default values in targets
#' @keywords internal
#' @noRd
.fill_defaults <- function(targets, multistage) {
  if (!"alpha" %in% names(targets)) {
    targets$alpha <- 0.05
  } else {
    targets$alpha[is.na(targets$alpha)] <- 0.05
  }

  if (!multistage) {
    if (!"deff" %in% names(targets)) {
      targets$deff <- 1
    } else {
      targets$deff[is.na(targets$deff)] <- 1
    }

    if (!"N" %in% names(targets)) {
      targets$N <- Inf
    } else {
      targets$N[is.na(targets$N)] <- Inf
    }
  }

  if (multistage) {
    if (!"k_psu" %in% names(targets)) {
      targets$k_psu <- 1
    } else {
      targets$k_psu[is.na(targets$k_psu)] <- 1
    }
    if (!"k_ssu" %in% names(targets)) {
      targets$k_ssu <- 1
    } else {
      targets$k_ssu[is.na(targets$k_ssu)] <- 1
    }
  }

  if (!"resp_rate" %in% names(targets)) {
    targets$resp_rate <- 1
  } else {
    targets$resp_rate[is.na(targets$resp_rate)] <- 1
  }

  if (!"rel_var" %in% names(targets)) {
    targets$rel_var <- NA_real_
  }
  targets$rel_var <- .derive_rel_var(targets, require_all = multistage)

  targets
}


#' Derive unit relvariance from p or var/mu
#' @keywords internal
#' @noRd
.derive_rel_var <- function(targets, require_all = FALSE) {
  rv <- targets$rel_var
  needs <- is.na(rv)

  has_p <- "p" %in% names(targets)
  has_var <- "var" %in% names(targets)
  has_mu <- "mu" %in% names(targets)

  for (i in which(needs)) {
    if (has_p && !is.na(targets$p[i])) {
      rv[i] <- (1 - targets$p[i]) / targets$p[i]
    } else if (has_var && !is.na(targets$var[i])) {
      if (has_mu && !is.na(targets$mu[i])) {
        rv[i] <- targets$var[i] / targets$mu[i]^2
      } else if (require_all) {
        stop(
          sprintf("row %d: 'mu' is required to derive 'rel_var' from 'var'", i),
          call. = FALSE
        )
      }
    }
  }

  if (require_all && anyNA(rv)) {
    stop(
      sprintf(
        "row(s) %s: could not derive 'rel_var' - provide 'p', or 'var'+'mu', or 'rel_var' directly",
        paste(which(is.na(rv)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  rv
}

#' Simple mode: compute n per indicator, take max
#' @keywords internal
#' @noRd
.n_multi_simple <- function(targets) {
  n_vec <- .compute_simple_n(targets)

  idx <- which.max(n_vec)
  n_max <- n_vec[idx]

  labels <- if ("name" %in% names(targets)) {
    targets$name
  } else {
    seq_len(nrow(targets))
  }
  binding_name <- labels[idx]

  detail <- data.frame(
    name = labels,
    .n = n_vec,
    .binding = seq_len(nrow(targets)) == idx
  )

  .new_svyplan_n(
    n = n_max,
    type = "multi",
    targets = targets,
    detail = detail,
    binding = binding_name
  )
}

#' Compute per-indicator n using Wald-type formulas (inline)
#' @keywords internal
#' @noRd
.compute_simple_n <- function(targets) {
  nr <- nrow(targets)
  n_vec <- numeric(nr)

  has_p <- "p" %in% names(targets)
  has_moe <- "moe" %in% names(targets)

  for (i in seq_len(nr)) {
    z <- qnorm(1 - targets$alpha[i] / 2)
    deff <- targets$deff[i]
    N <- targets$N[i]

    is_prop <- has_p && !is.na(targets$p[i])
    use_moe <- has_moe && !is.na(targets$moe[i])

    if (is_prop) {
      p <- targets$p[i]
      q <- 1 - p
      if (use_moe) {
        a <- if (is.infinite(N)) 1 else N / (N - 1)
        n_vec[i] <- a * z^2 * p * q / (targets$moe[i]^2 + z^2 * p * q / (N - 1))
      } else {
        cv_t <- targets$cv[i]
        a <- if (is.infinite(N)) 1 else N / (N - 1)
        n_vec[i] <- a * q / p / (cv_t^2 + q / p / (N - 1))
      }
    } else {
      v <- targets$var[i]
      if (use_moe) {
        n_vec[i] <- z^2 * v / (targets$moe[i]^2 + z^2 * v / N)
      } else {
        cv_t <- targets$cv[i]
        mu <- targets$mu[i]
        CVpop <- sqrt(v) / mu
        n_vec[i] <- CVpop^2 / (cv_t^2 + CVpop^2 / N)
      }
    }

    n_vec[i] <- n_vec[i] * deff
    n_vec[i] <- .apply_resp_rate(n_vec[i], targets$resp_rate[i])
  }

  n_vec
}

#' Multistage cluster mode dispatcher
#' @keywords internal
#' @noRd
.n_multi_cluster <- function(targets, cost, budget, n_psu, fixed_cost = 0) {
  stages <- length(cost)
  if (stages == 2L) {
    .n_multi_2stage(targets, cost, budget, n_psu, fixed_cost)
  } else {
    .n_multi_3stage(targets, cost, budget, n_psu, fixed_cost)
  }
}

#' 2-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_2stage <- function(targets, cost, budget, n_psu, fixed_cost = 0) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  nr <- nrow(targets)

  cv_t <- targets$cv
  delta <- targets$delta_psu
  rel_var <- targets$rel_var
  k <- targets$k_psu
  rr <- targets$resp_rate
  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  # n1_actual_j = rel_var_j * k_j * (1 + delta_j*(psu_size-1)) / (psu_size * cv_j^2 * rr_j)
  n1_required <- function(psu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        rel_var[j] * k[j] * (1 + delta[j] * (psu_size - 1)) / (psu_size * cv_t[j]^2 * rr[j])
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, psu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] * k[j] / (n1 * rr[j] * psu_size) * (1 + delta[j] * (psu_size - 1))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(budget)) {
    cost_fn <- function(psu_size) {
      n1 <- max(n1_required(psu_size))
      n1 * (C1 + C2 * psu_size)
    }

    upper <- max(10, 10 * sqrt(C1 / C2))
    opt <- optimize(cost_fn, interval = c(1, upper))
    psu_size_opt <- opt$minimum

    n1_vals <- n1_required(psu_size_opt)
    n1_opt <- max(n1_vals)
    total_cost <- fixed_cost + n1_opt * (C1 + C2 * psu_size_opt)
    binding_idx <- which.max(n1_vals)

    cv_achieved <- cv_achieved_fn(n1_opt, psu_size_opt)
  } else {
    var_budget <- budget - fixed_cost
    bres <- .eval_2stage_budget(
      cv_t,
      delta,
      rel_var,
      k,
      rr,
      C1,
      C2,
      var_budget,
      n_psu
    )
    n1_opt <- bres$n1
    psu_size_opt <- bres$psu_size
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- budget
  }

  n_vec <- c(n_psu = n1_opt, psu_size = psu_size_opt)
  total_n <- prod(n_vec)

  detail <- data.frame(
    name = labels,
    .cv_target = cv_t,
    .cv_achieved = cv_achieved,
    .binding = seq_len(nr) == binding_idx
  )

  params <- list(cost = cost)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) {
    params$n_psu <- n_psu
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 2L,
    total_n = total_n,
    cv = cv_achieved[binding_idx],
    cost = total_cost,
    params = params,
    targets = targets,
    detail = detail,
    binding = labels[binding_idx]
  )
}

#' 3-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_3stage <- function(targets, cost, budget, n_psu, fixed_cost = 0) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  C3 <- cost[3L]
  nr <- nrow(targets)

  cv_t <- targets$cv
  delta_psu <- targets$delta_psu
  delta_ssu <- if ("delta_ssu" %in% names(targets)) {
    targets$delta_ssu
  } else {
    rep(0, nr)
  }
  rel_var <- targets$rel_var
  k_psu <- targets$k_psu
  k_ssu <- targets$k_ssu
  rr <- targets$resp_rate
  labels <- if ("name" %in% names(targets)) {
    targets$name
  } else {
    seq_len(nr)
  }

  n1_required <- function(psu_size, ssu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        rel_var[j] /
          (cv_t[j]^2 * psu_size * ssu_size * rr[j]) *
          (k_psu[j] * delta_psu[j] * psu_size * ssu_size + k_ssu[j] * (1 + delta_ssu[j] * (ssu_size - 1)))
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, psu_size, ssu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] /
            (n1 * rr[j] * psu_size * ssu_size) *
            (k_psu[j] * delta_psu[j] * psu_size * ssu_size + k_ssu[j] * (1 + delta_ssu[j] * (ssu_size - 1)))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(budget)) {
    if (is.null(n_psu)) {
      cost_fn <- function(par) {
        psu_size <- par[1L]
        ssu_size <- par[2L]
        n1 <- max(n1_required(psu_size, ssu_size))
        n1 * (C1 + C2 * psu_size + C3 * psu_size * ssu_size)
      }

      init_ssu_size <- max(2, sqrt(C2 / C3))
      init_psu_size <- max(2, sqrt(C1 / C2))
      opt <- optim(
        par = c(init_psu_size, init_ssu_size),
        fn = cost_fn,
        method = "L-BFGS-B",
        lower = c(1, 1),
        upper = c(Inf, Inf)
      )
      if (opt$convergence != 0L) {
        warning(
          "L-BFGS-B did not converge (code ",
          opt$convergence,
          "): ",
          "allocation may be approximate",
          call. = FALSE
        )
      }

      psu_size_opt <- opt$par[1L]
      ssu_size_opt <- opt$par[2L]
      n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
      n1_opt <- max(n1_vals)
      total_cost <- fixed_cost +
        n1_opt * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
      binding_idx <- which.max(n1_vals)
    } else {
      n1_opt <- n_psu

      psu_size_required_fn <- function(ssu_size) {
        psu_size_per <- vapply(
          seq_len(nr),
          function(j) {
            denom <- cv_t[j]^2 *
              n_psu *
              rr[j] /
              (rel_var[j] * k_ssu[j]) -
              k_psu[j] * delta_psu[j] / k_ssu[j]
            if (denom <= 0) {
              Inf
            } else {
              (1 + delta_ssu[j] * (ssu_size - 1)) / (ssu_size * denom)
            }
          },
          numeric(1L)
        )
        max(psu_size_per)
      }

      cost_fn_fixed <- function(ssu_size) {
        psu_size <- psu_size_required_fn(ssu_size)
        n_psu * (C1 + C2 * psu_size + C3 * psu_size * ssu_size)
      }

      ssu_size_analytic <- vapply(
        seq_len(nr),
        function(j) {
          if (delta_ssu[j] <= 0) 1 else sqrt((1 - delta_ssu[j]) / delta_ssu[j] * C2 / C3)
        },
        numeric(1L)
      )
      upper_ssu_size <- max(10, 3 * max(ssu_size_analytic))

      opt <- optimize(cost_fn_fixed, interval = c(1, upper_ssu_size))
      ssu_size_opt <- opt$minimum
      psu_size_opt <- psu_size_required_fn(ssu_size_opt)

      if (!is.finite(psu_size_opt) || psu_size_opt <= 0) {
        stop(
          "target CV is too small for the given fixed stage-1 size",
          call. = FALSE
        )
      }

      n1_vals <- n1_required(psu_size_opt, ssu_size_opt)
      binding_idx <- which.max(n1_vals)
      total_cost <- fixed_cost +
        n_psu * (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
    }

    cv_achieved <- cv_achieved_fn(n1_opt, psu_size_opt, ssu_size_opt)
  } else {
    var_budget <- budget - fixed_cost
    bres <- .eval_3stage_budget(
      cv_t,
      delta_psu,
      delta_ssu,
      rel_var,
      k_psu,
      k_ssu,
      rr,
      C1,
      C2,
      C3,
      var_budget,
      n_psu
    )
    n1_opt <- bres$n1
    psu_size_opt <- bres$psu_size
    ssu_size_opt <- bres$ssu_size
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- fixed_cost + bres$cost
  }

  n_vec <- c(n_psu = n1_opt, psu_size = psu_size_opt, ssu_size = ssu_size_opt)
  total_n <- prod(n_vec)

  detail <- data.frame(
    name = labels,
    .cv_target = cv_t,
    .cv_achieved = cv_achieved,
    .binding = seq_len(nr) == binding_idx
  )

  params <- list(cost = cost)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) {
    params$n_psu <- n_psu
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = 3L,
    total_n = total_n,
    cv = cv_achieved[binding_idx],
    cost = total_cost,
    params = params,
    targets = targets,
    detail = detail,
    binding = labels[binding_idx]
  )
}

#' Solve n_multi independently per domain, then aggregate
#' @keywords internal
#' @noRd
.n_multi_domains <- function(
  targets,
  cost,
  budget,
  n_psu,
  domain_cols,
  multistage,
  joint = FALSE,
  min_n = NULL,
  fixed_cost = 0
) {
  domain_keys <- interaction(targets[domain_cols], drop = TRUE, sep = ":")
  domain_levels <- levels(domain_keys)
  split_idx <- split(seq_len(nrow(targets)), domain_keys)

  if (multistage && joint && !is.null(budget) && length(domain_levels) > 1L) {
    return(.n_multi_domains_joint(
      targets,
      cost,
      budget,
      n_psu,
      domain_cols,
      domain_levels,
      split_idx,
      min_n,
      fixed_cost
    ))
  }

  results <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- targets[rows, , drop = FALSE]
    # Drop domain columns for inner solve
    sub_inner <- sub[, !names(sub) %in% domain_cols, drop = FALSE]
    if (multistage) {
      .n_multi_cluster(sub_inner, cost, budget, n_psu, fixed_cost)
    } else {
      .n_multi_simple(sub_inner)
    }
  })
  names(results) <- domain_levels

  if (multistage) {
    .aggregate_cluster_domains(
      results,
      targets,
      domain_cols,
      domain_levels,
      split_idx,
      cost,
      budget,
      n_psu,
      min_n,
      fixed_cost
    )
  } else {
    .aggregate_simple_domains(
      results,
      targets,
      domain_cols,
      domain_levels,
      split_idx,
      min_n
    )
  }
}

#' Aggregate simple-mode domain results
#' @keywords internal
#' @noRd
.aggregate_simple_domains <- function(
  results,
  targets,
  domain_cols,
  domain_levels,
  split_idx,
  min_n = NULL
) {
  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- targets[rows[1L], domain_cols, drop = FALSE]
    dom_vals$.n <- res$n
    dom_vals$.binding <- res$binding
    dom_vals
  })
  domains <- do.call(rbind, domain_rows)
  rownames(domains) <- NULL

  if (!is.null(min_n)) {
    floored <- domains$.n < min_n
    domains$.n[floored] <- min_n
    domains$.binding[floored] <- "(min_n)"
  }

  n_max <- max(domains$.n)
  overall_binding_idx <- which.max(domains$.n)
  binding_label <- domains$.binding[overall_binding_idx]

  res <- .new_svyplan_n(
    n = n_max,
    type = "multi",
    targets = targets,
    detail = NULL,
    binding = binding_label,
    domains = domains
  )
  if (!is.null(min_n)) {
    res$params$min_n <- min_n
  }
  res
}

#' Aggregate multistage domain results
#' @keywords internal
#' @noRd
.aggregate_cluster_domains <- function(
  results,
  targets,
  domain_cols,
  domain_levels,
  split_idx,
  cost,
  budget = NULL,
  n_psu = NULL,
  min_n = NULL,
  fixed_cost = 0
) {
  stages <- results[[1L]]$stages
  stage_names <- names(results[[1L]]$n)

  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- targets[rows[1L], domain_cols, drop = FALSE]
    for (s in seq_len(stages)) {
      dom_vals[[stage_names[s]]] <- res$n[s]
    }
    dom_vals$.total_n <- res$total_n
    dom_vals$.cv <- res$cv
    dom_vals$.cost <- res$cost
    dom_vals$.binding <- res$binding
    dom_vals
  })
  domains <- do.call(rbind, domain_rows)
  rownames(domains) <- NULL

  if (!is.null(min_n)) {
    below <- which(domains$.total_n < min_n)
    if (length(below) > 0L) {
      dom_labels <- vapply(
        below,
        function(i) {
          paste(domains[i, domain_cols, drop = TRUE], collapse = ":")
        },
        character(1L)
      )
      warning(
        sprintf(
          "domain(s) %s have total_n below min_n = %g",
          paste(sQuote(dom_labels), collapse = ", "),
          min_n
        ),
        call. = FALSE
      )
    }
  }

  total_n <- sum(ceiling(domains$.total_n))
  total_cost <- sum(domains$.cost)
  worst_cv_idx <- which.max(domains$.cv)

  n_vec <- vapply(
    seq_len(stages),
    function(s) {
      max(domains[[stage_names[s]]])
    },
    numeric(1L)
  )
  names(n_vec) <- stage_names

  params <- list(cost = cost)
  if (!is.null(budget)) {
    params$budget <- budget
  }
  if (!is.null(n_psu)) {
    params$n_psu <- n_psu
  }
  if (!is.null(min_n)) {
    params$min_n <- min_n
  }
  if (fixed_cost > 0) {
    params$fixed_cost <- fixed_cost
  }

  .new_svyplan_cluster(
    n = n_vec,
    stages = stages,
    total_n = total_n,
    cv = domains$.cv[worst_cv_idx],
    cost = total_cost,
    params = params,
    targets = targets,
    detail = NULL,
    binding = domains$.binding[worst_cv_idx],
    domains = domains
  )
}

#' @keywords internal
#' @noRd
.n_multi_domains_joint <- function(
  targets,
  cost,
  budget,
  n_psu,
  domain_cols,
  domain_levels,
  split_idx,
  min_n = NULL,
  fixed_cost = 0
) {
  stages <- length(cost)
  nd <- length(domain_levels)
  C1 <- cost[1L]
  C2 <- cost[2L]
  C3 <- if (stages == 3L) cost[3L] else NULL

  domain_params <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- targets[rows, , drop = FALSE]
    labels <- if ("name" %in% names(sub)) {
      sub$name
    } else {
      seq_along(rows)
    }
    list(
      cv_t = sub$cv,
      delta_psu = sub$delta_psu,
      delta_ssu = if ("delta_ssu" %in% names(sub)) {
        sub$delta_ssu
      } else {
        rep(0, length(rows))
      },
      rel_var = sub$rel_var,
      k_psu = sub$k_psu,
      k_ssu = sub$k_ssu,
      resp_rate = sub$resp_rate,
      labels = labels
    )
  })

  var_budget <- budget - fixed_cost

  eval_domain <- function(d, budget_d) {
    p <- domain_params[[d]]
    if (stages == 2L) {
      .eval_2stage_budget(
        p$cv_t,
        p$delta_psu,
        p$rel_var,
        p$k_psu,
        p$resp_rate,
        C1,
        C2,
        budget_d,
        n_psu
      )
    } else {
      .eval_3stage_budget(
        p$cv_t,
        p$delta_psu,
        p$delta_ssu,
        p$rel_var,
        p$k_psu,
        p$k_ssu,
        p$resp_rate,
        C1,
        C2,
        C3,
        budget_d,
        n_psu
      )
    }
  }

  total_n_for <- function(bres) {
    if (stages == 2L) {
      bres$n1 * bres$psu_size
    } else {
      bres$n1 * bres$psu_size * bres$ssu_size
    }
  }

  if (!is.null(min_n)) {
    full_total <- vapply(
      seq_len(nd),
      function(d) {
        tryCatch(total_n_for(eval_domain(d, var_budget)), error = function(e) 0)
      },
      numeric(1L)
    )
    for (d in seq_len(nd)) {
      if (full_total[d] < min_n) {
        stop(
          sprintf(
            "min_n = %g not achievable for domain '%s' (max total_n = %.0f at full budget)",
            min_n,
            domain_levels[d],
            full_total[d]
          ),
          call. = FALSE
        )
      }
    }
    min_fracs <- min_n / full_total
    if (sum(min_fracs) > 1) {
      stop(
        sprintf(
          "min_n = %g not achievable for all domains within budget",
          min_n
        ),
        call. = FALSE
      )
    }
  }

  lower_bounds <- rep(1e-4, nd)
  if (!is.null(min_n)) {
    lower_bounds <- pmax(lower_bounds, min_fracs)
  }

  outer_obj <- function(w) {
    w_last <- 1 - sum(w)
    if (w_last < lower_bounds[nd]) {
      return(1e12)
    }
    fracs <- c(w, w_last)
    tryCatch(
      {
        vals <- vapply(
          seq_len(nd),
          function(d) {
            bres <- eval_domain(d, fracs[d] * var_budget)
            if (!is.null(min_n) && total_n_for(bres) < min_n) {
              return(1e12)
            }
            bres$ratio
          },
          numeric(1L)
        )
        max(vals)
      },
      error = function(e) 1e12
    )
  }

  if (nd == 2L) {
    opt <- optimize(
      function(w1) outer_obj(w1),
      interval = c(lower_bounds[1], 1 - lower_bounds[2])
    )
    w_opt <- c(opt$minimum, 1 - opt$minimum)
  } else {
    init_w <- rep(1 / nd, nd - 1L)
    init_w <- pmax(init_w, lower_bounds[-nd])
    if (sum(init_w) >= 1 - lower_bounds[nd]) {
      init_w <- lower_bounds[-nd] / sum(lower_bounds) * (1 - lower_bounds[nd])
    }
    opt <- optim(
      par = init_w,
      fn = outer_obj,
      method = "L-BFGS-B",
      lower = lower_bounds[-nd],
      upper = rep(1 - lower_bounds[nd], nd - 1L)
    )
    if (opt$convergence != 0L) {
      warning(
        "L-BFGS-B did not converge (code ",
        opt$convergence,
        "): ",
        "allocation may be approximate",
        call. = FALSE
      )
    }
    w_opt <- c(opt$par, 1 - sum(opt$par))
  }

  budgets <- w_opt * var_budget
  results <- lapply(seq_len(nd), function(d) {
    bres <- eval_domain(d, budgets[d])
    p <- domain_params[[d]]
    if (stages == 2L) {
      n_vec <- c(n_psu = bres$n1, psu_size = bres$psu_size)
    } else {
      n_vec <- c(n_psu = bres$n1, psu_size = bres$psu_size, ssu_size = bres$ssu_size)
    }
    .new_svyplan_cluster(
      n = n_vec,
      stages = stages,
      total_n = prod(n_vec),
      cv = bres$cv_achieved[bres$binding_idx],
      cost = bres$cost,
      params = list(cost = cost),
      targets = NULL,
      detail = NULL,
      binding = p$labels[bres$binding_idx]
    )
  })
  names(results) <- domain_levels

  res <- .aggregate_cluster_domains(
    results,
    targets,
    domain_cols,
    domain_levels,
    split_idx,
    cost,
    budget,
    n_psu,
    min_n,
    fixed_cost
  )
  res$params$joint <- TRUE
  if (fixed_cost > 0) {
    res$cost <- budget
  }
  res
}

#' @keywords internal
#' @noRd
.eval_2stage_budget <- function(
  cv_t,
  delta,
  rel_var,
  k,
  resp_rate,
  C1,
  C2,
  budget,
  n_psu
) {
  nr <- length(cv_t)

  if (is.null(n_psu)) {
    obj_fn <- function(psu_size) {
      n1 <- budget / (C1 + C2 * psu_size)
      if (n1 <= 0) {
        return(1e12)
      }
      cv_ratios <- vapply(
        seq_len(nr),
        function(j) {
          sqrt(
            rel_var[j] *
              k[j] /
              (n1 * resp_rate[j] * psu_size) *
              (1 + delta[j] * (psu_size - 1))
          ) /
            cv_t[j]
        },
        numeric(1L)
      )
      max(cv_ratios)
    }

    upper <- max(10, budget / (C1 + C2))
    opt <- optimize(obj_fn, interval = c(1, upper))
    psu_size_opt <- opt$minimum
    n1_opt <- budget / (C1 + C2 * psu_size_opt)
  } else {
    n1_opt <- n_psu
    psu_size_opt <- (budget - C1 * n_psu) / (C2 * n_psu)
    if (psu_size_opt <= 0) {
      stop(
        "budget is too small for the given fixed stage-1 size",
        call. = FALSE
      )
    }
  }

  cv_achieved <- vapply(
    seq_len(nr),
    function(j) {
      sqrt(
        rel_var[j] *
          k[j] /
          (n1_opt * resp_rate[j] * psu_size_opt) *
          (1 + delta[j] * (psu_size_opt - 1))
      )
    },
    numeric(1L)
  )
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    psu_size = psu_size_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = budget
  )
}

#' @keywords internal
#' @noRd
.eval_3stage_budget <- function(
  cv_t,
  delta_psu,
  delta_ssu,
  rel_var,
  k_psu,
  k_ssu,
  resp_rate,
  C1,
  C2,
  C3,
  budget,
  n_psu
) {
  nr <- length(cv_t)

  cv_fn <- function(n1, psu_size, ssu_size) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] /
            (n1 * resp_rate[j] * psu_size * ssu_size) *
            (k_psu[j] * delta_psu[j] * psu_size * ssu_size + k_ssu[j] * (1 + delta_ssu[j] * (ssu_size - 1)))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(n_psu)) {
    obj_fn_free <- function(par) {
      psu_size <- par[1L]
      ssu_size <- par[2L]
      n1 <- budget / (C1 + C2 * psu_size + C3 * psu_size * ssu_size)
      if (n1 <= 0) {
        return(1e12)
      }
      max(cv_fn(n1, psu_size, ssu_size) / cv_t)
    }

    init_ssu_size <- max(2, sqrt(C2 / C3))
    init_psu_size <- max(2, sqrt(C1 / C2))
    opt <- optim(
      par = c(init_psu_size, init_ssu_size),
      fn = obj_fn_free,
      method = "L-BFGS-B",
      lower = c(1, 1),
      upper = c(Inf, Inf)
    )
    if (opt$convergence != 0L) {
      warning(
        "L-BFGS-B did not converge (code ",
        opt$convergence,
        "): ",
        "allocation may be approximate",
        call. = FALSE
      )
    }

    psu_size_opt <- opt$par[1L]
    ssu_size_opt <- opt$par[2L]
    n1_opt <- budget / (C1 + C2 * psu_size_opt + C3 * psu_size_opt * ssu_size_opt)
    total_cost <- budget
  } else {
    n1_opt <- n_psu
    min_cost <- n_psu * (C1 + C2 + C3)
    if (min_cost > budget) {
      stop(
        "budget is too small for the given fixed stage-1 size",
        call. = FALSE
      )
    }

    max_ssu_size <- (budget - C1 * n_psu - C2 * n_psu) / (C3 * n_psu)
    if (!is.finite(max_ssu_size) || max_ssu_size < 1) {
      stop(
        "budget is too small for the given fixed stage-1 size",
        call. = FALSE
      )
    }

    psu_size_from_ssu <- function(ssu_size) {
      (budget / n_psu - C1) / (C2 + C3 * ssu_size)
    }

    obj_fn_fixed <- function(ssu_size) {
      psu_size <- psu_size_from_ssu(ssu_size)
      if (!is.finite(psu_size) || psu_size < 1) {
        return(Inf)
      }
      max(cv_fn(n_psu, psu_size, ssu_size) / cv_t)
    }

    if (max_ssu_size <= 1 + 1e-10) {
      ssu_size_opt <- 1
    } else {
      opt <- optimize(obj_fn_fixed, interval = c(1, max_ssu_size))
      ssu_size_opt <- opt$minimum
    }
    psu_size_opt <- psu_size_from_ssu(ssu_size_opt)

    if (!is.finite(psu_size_opt) || psu_size_opt < 1) {
      stop(
        "failed to find a feasible allocation under the given budget and fixed stage-1 size",
        call. = FALSE
      )
    }
    total_cost <- C1 * n_psu + C2 * n_psu * psu_size_opt + C3 * n_psu * psu_size_opt * ssu_size_opt

    tol <- max(1e-8, 1e-6 * max(1, budget))
    if (total_cost - budget > tol) {
      stop(
        "failed to find a feasible allocation under the given budget and fixed stage-1 size",
        call. = FALSE
      )
    }
  }

  cv_achieved <- cv_fn(n1_opt, psu_size_opt, ssu_size_opt)
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    psu_size = psu_size_opt,
    ssu_size = ssu_size_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = total_cost
  )
}

#' @rdname n_multi
#' @export
n_multi.svyplan_prec <- function(targets, cost = NULL, ...) {
  x <- targets
  if (x$type != "multi") {
    stop("n_multi requires a svyplan_prec of type 'multi'", call. = FALSE)
  }
  tgt <- x$params$targets
  tgt$n <- NULL
  tgt$psu_size <- NULL
  tgt$ssu_size <- NULL
  if (!is.null(x$detail)) {
    if (".moe" %in% names(x$detail) && !all(is.na(x$detail$.moe))) {
      tgt$moe <- x$detail$.moe
      tgt$cv <- NULL
    } else if (".cv" %in% names(x$detail) && !all(is.na(x$detail$.cv))) {
      tgt$cv <- x$detail$.cv
      tgt$moe <- NULL
    }
  }
  cost <- cost %||% x$params$cost
  budget <- x$params$budget
  n_psu <- x$params$n_psu
  joint <- if (!is.null(x$params$joint)) x$params$joint else FALSE
  min_n <- x$params$min_n
  fixed_cost <- if (!is.null(x$params$fixed_cost)) x$params$fixed_cost else 0
  n_multi.default(
    targets = tgt,
    cost = cost,
    budget = budget,
    n_psu = n_psu,
    joint = joint,
    min_n = min_n,
    fixed_cost = fixed_cost
  )
}
