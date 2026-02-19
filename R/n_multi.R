#' Multi-Indicator Sample Size
#'
#' Compute the sample size that satisfies precision requirements for
#' multiple survey indicators simultaneously. Supports simple (single-stage)
#' and multistage cluster designs, with optional domain-level planning.
#'
#' @param targets Data frame with one row per indicator. See Details for
#'   column specifications.
#' @param cost Numeric vector of per-stage costs. `NULL` (default) for
#'   simple mode; length 2 or 3 for multistage mode.
#' @param budget Total budget (multistage only). Specify exactly one of
#'   `cv` targets or `budget`.
#' @param m Fixed stage-1 sample size (multistage only).
#' @param joint Logical. If `TRUE`, optimally split a single `budget`
#'   across domains to minimize the worst-case CV ratio. Only applies
#'   to multistage budget mode with multiple domains; ignored otherwise.
#'
#' @return A `svyplan_n` object (simple mode) or `svyplan_cluster` object
#'   (multistage mode).
#'
#'   **Without domains**, the object contains:
#'   \describe{
#'     \item{`n`}{Sample size (simple) or named per-stage allocation vector
#'       (multistage, e.g. `c(n1 = 80, n2 = 12)`).}
#'     \item{`detail`}{Per-indicator results (sample sizes or achieved CVs).}
#'     \item{`binding`}{Name or index of the binding (most demanding) indicator.}
#'     \item{`targets`}{The input targets data frame.}
#'   }
#'
#'   **With domains**, the object additionally contains:
#'   \describe{
#'     \item{`n`}{Maximum per-stage sample size across domains. In simple
#'       mode, a single number; in multistage mode, a named vector
#'       (e.g. `c(n1 = 120, n2 = 15)`) giving the conservative
#'       allocation that satisfies all domains.}
#'     \item{`domains`}{Data frame with one row per domain, including
#'       domain variable columns, per-stage allocations (`n1`, `n2`, ...),
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
#'   \item{`mu`}{Population mean. Required when `var` is specified with `cv`.}
#'   \item{`moe`}{Margin of error (simple mode).}
#'   \item{`cv`}{Target coefficient of variation (either mode).}
#'   \item{`alpha`}{Significance level (default 0.05).}
#'   \item{`deff`}{Design effect multiplier (simple mode only, default 1).}
#'   \item{`N`}{Population size (simple mode only, default Inf).}
#'   \item{`delta1`, `delta2`}{Homogeneity measures (multistage).}
#'   \item{`rel_var`}{Unit relvariance. If omitted, derived from `p` or
#'     `var`/`mu`.}
#'   \item{`k1`, `k2`}{Ratio parameters (multistage, default 1).}
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
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018).
#' *Practical Tools for Designing and Weighting Survey Samples*
#' (2nd ed.). Springer.
#'
#' @seealso [n_prop()], [n_mean()] for single-indicator sizing;
#'   [n_cluster()] for single-indicator multistage allocation.
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
#'   delta1 = c(0.02, 0.05)
#' )
#' n_multi(targets_cl, cost = c(500, 50))
#'
#' # Joint budget allocation across domains
#' targets_jnt <- data.frame(
#'   name   = rep(c("stunting", "anemia"), each = 2),
#'   p      = c(0.30, 0.25, 0.10, 0.15),
#'   cv     = c(0.10, 0.10, 0.15, 0.15),
#'   delta1 = c(0.02, 0.03, 0.05, 0.04),
#'   region = rep(c("Urban", "Rural"), 2)
#' )
#' n_multi(targets_jnt, cost = c(500, 50), budget = 100000, joint = TRUE)
#'
#' @export
n_multi <- function(targets, cost = NULL, budget = NULL, m = NULL,
                    joint = FALSE) {
  if (!is.data.frame(targets) || nrow(targets) == 0L) {
    stop("'targets' must be a non-empty data frame", call. = FALSE)
  }
  if (!is.logical(joint) || length(joint) != 1L || is.na(joint)) {
    stop("'joint' must be TRUE or FALSE", call. = FALSE)
  }

  multistage <- !is.null(cost)

  if (multistage) {
    check_cost(cost)
    if (!is.null(budget)) {
      check_scalar(budget, "budget")
    }
    if (!is.null(m)) check_scalar(m, "m")
  } else {
    if (!is.null(budget)) {
      stop("'budget' requires 'cost' to be specified", call. = FALSE)
    }
    if (!is.null(m)) {
      stop("'m' requires 'cost' to be specified", call. = FALSE)
    }
  }

  info <- .validate_targets(targets, multistage)
  targets <- .fill_defaults(targets, multistage)

  domain_cols <- info$domain_cols

  if (length(domain_cols) == 0L) {
    if (multistage) {
      .n_multi_cluster(targets, cost, budget, m)
    } else {
      .n_multi_simple(targets)
    }
  } else {
    .n_multi_domains(targets, cost, budget, m, domain_cols, multistage, joint)
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
    "delta1",
    "delta2",
    "rel_var",
    "k1",
    "k2"
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

  # Each row needs either p or var
  if (has_p && has_var) {
    both_na <- is.na(targets$p) & is.na(targets$var)
    both_set <- !is.na(targets$p) & !is.na(targets$var)
    if (any(both_na)) {
      stop("each row must have either 'p' or 'var' specified", call. = FALSE)
    }
    if (any(both_set)) {
      stop("each row must have only one of 'p' or 'var'", call. = FALSE)
    }
  }

  # Each row needs moe or cv (not both)
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

  # Multistage requires cv and delta1
  if (multistage) {
    if (!has_cv) {
      stop("multistage mode requires 'cv' column in targets", call. = FALSE)
    }
    if (has_moe) {
      if (any(!is.na(targets$moe))) {
        stop("multistage mode requires 'cv' (not 'moe')", call. = FALSE)
      }
    }
    if (!"delta1" %in% names(targets)) {
      stop("multistage mode requires 'delta1' column in targets", call. = FALSE)
    }
  }

  # var + cv requires mu
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

  domain_cols <- setdiff(names(targets), known)

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
    if (!"k1" %in% names(targets)) {
      targets$k1 <- 1
    } else {
      targets$k1[is.na(targets$k1)] <- 1
    }
    if (!"k2" %in% names(targets)) {
      targets$k2 <- 1
    } else {
      targets$k2[is.na(targets$k2)] <- 1
    }
  }

  # Derive rel_var (only needed for multistage, but compute when possible)
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
    .binding = seq_len(nrow(targets)) == idx,
    stringsAsFactors = FALSE
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
  has_var <- "var" %in% names(targets)
  has_moe <- "moe" %in% names(targets)
  has_cv <- "cv" %in% names(targets)
  has_mu <- "mu" %in% names(targets)

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
  }

  n_vec
}

#' Multistage cluster mode dispatcher
#' @keywords internal
#' @noRd
.n_multi_cluster <- function(targets, cost, budget, m) {
  stages <- length(cost)
  if (stages == 2L) {
    .n_multi_2stage(targets, cost, budget, m)
  } else {
    .n_multi_3stage(targets, cost, budget, m)
  }
}

#' 2-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_2stage <- function(targets, cost, budget, m) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  nr <- nrow(targets)

  cv_t <- targets$cv
  delta <- targets$delta1
  rel_var <- targets$rel_var
  k <- targets$k1
  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  # n1 required for given n2, for indicator j
  # n1_j(n2) = rel_var_j * k_j * (1 + delta_j*(n2-1)) / (n2 * cv_j^2)
  n1_required <- function(n2) {
    vapply(
      seq_len(nr),
      function(j) {
        rel_var[j] * k[j] * (1 + delta[j] * (n2 - 1)) / (n2 * cv_t[j]^2)
      },
      numeric(1L)
    )
  }

  if (is.null(budget)) {
    # CV mode: minimize cost
    cost_fn <- function(n2) {
      n1 <- max(n1_required(n2))
      n1 * (C1 + C2 * n2)
    }

    upper <- max(10, 10 * sqrt(C1 / C2))
    opt <- optimize(cost_fn, interval = c(1, upper))
    n2_opt <- opt$minimum

    n1_vals <- n1_required(n2_opt)
    n1_opt <- max(n1_vals)
    total_cost <- n1_opt * (C1 + C2 * n2_opt)
    binding_idx <- which.max(n1_vals)

    cv_achieved <- vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] * k[j] / (n1_opt * n2_opt) * (1 + delta[j] * (n2_opt - 1))
        )
      },
      numeric(1L)
    )
  } else {
    bres <- .eval_2stage_budget(cv_t, delta, rel_var, k, C1, C2, budget, m)
    n1_opt <- bres$n1
    n2_opt <- bres$n2
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- bres$cost
  }

  n_vec <- c(n1 = n1_opt, n2 = n2_opt)
  total_n <- prod(n_vec)

  detail <- data.frame(
    name = labels,
    .cv_target = cv_t,
    .cv_achieved = cv_achieved,
    .binding = seq_len(nr) == binding_idx,
    stringsAsFactors = FALSE
  )

  .new_svyplan_cluster(
    n = n_vec,
    stages = 2L,
    total_n = total_n,
    cv = cv_achieved[binding_idx],
    cost = total_cost,
    params = list(cost = cost),
    targets = targets,
    detail = detail,
    binding = labels[binding_idx]
  )
}


#' 3-stage multi-indicator optimization
#' @keywords internal
#' @noRd
.n_multi_3stage <- function(targets, cost, budget, m) {
  C1 <- cost[1L]
  C2 <- cost[2L]
  C3 <- cost[3L]
  nr <- nrow(targets)

  cv_t <- targets$cv
  delta1 <- targets$delta1
  delta2 <- if ("delta2" %in% names(targets)) targets$delta2 else rep(0, nr)
  rel_var <- targets$rel_var
  k1 <- targets$k1
  k2 <- targets$k2
  labels <- if ("name" %in% names(targets)) targets$name else seq_len(nr)

  # n1 required for given (n2, n3), for indicator j
  n1_required <- function(n2, n3) {
    vapply(
      seq_len(nr),
      function(j) {
        rel_var[j] /
          (cv_t[j]^2 * n2 * n3) *
          (k1[j] * delta1[j] * n2 * n3 + k2[j] * (1 + delta2[j] * (n3 - 1)))
      },
      numeric(1L)
    )
  }

  cv_achieved_fn <- function(n1, n2, n3) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] /
            (n1 * n2 * n3) *
            (k1[j] * delta1[j] * n2 * n3 + k2[j] * (1 + delta2[j] * (n3 - 1)))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(budget)) {
    # CV mode: minimize cost
    cost_fn <- function(par) {
      n2 <- par[1L]
      n3 <- par[2L]
      n1 <- max(n1_required(n2, n3))
      n1 * (C1 + C2 * n2 + C3 * n2 * n3)
    }

    # Starting values from single-indicator heuristic
    init_n3 <- max(2, sqrt(C2 / C3))
    init_n2 <- max(2, sqrt(C1 / C2))
    opt <- optim(
      par = c(init_n2, init_n3),
      fn = cost_fn,
      method = "L-BFGS-B",
      lower = c(1, 1),
      upper = c(1000, 1000)
    )

    n2_opt <- opt$par[1L]
    n3_opt <- opt$par[2L]
    n1_vals <- n1_required(n2_opt, n3_opt)
    n1_opt <- max(n1_vals)
    total_cost <- n1_opt * (C1 + C2 * n2_opt + C3 * n2_opt * n3_opt)
    binding_idx <- which.max(n1_vals)

    cv_achieved <- cv_achieved_fn(n1_opt, n2_opt, n3_opt)
  } else {
    bres <- .eval_3stage_budget(
      cv_t, delta1, delta2, rel_var, k1, k2, C1, C2, C3, budget, m
    )
    n1_opt <- bres$n1
    n2_opt <- bres$n2
    n3_opt <- bres$n3
    cv_achieved <- bres$cv_achieved
    binding_idx <- bres$binding_idx
    total_cost <- bres$cost
  }

  n_vec <- c(n1 = n1_opt, n2 = n2_opt, n3 = n3_opt)
  total_n <- prod(n_vec)

  detail <- data.frame(
    name = labels,
    .cv_target = cv_t,
    .cv_achieved = cv_achieved,
    .binding = seq_len(nr) == binding_idx,
    stringsAsFactors = FALSE
  )

  .new_svyplan_cluster(
    n = n_vec,
    stages = 3L,
    total_n = total_n,
    cv = cv_achieved[binding_idx],
    cost = total_cost,
    params = list(cost = cost),
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
  m,
  domain_cols,
  multistage,
  joint = FALSE
) {
  domain_keys <- interaction(targets[domain_cols], drop = TRUE, sep = ":")
  domain_levels <- levels(domain_keys)
  split_idx <- split(seq_len(nrow(targets)), domain_keys)

  if (multistage && joint && !is.null(budget) && length(domain_levels) > 1L) {
    return(.n_multi_domains_joint(
      targets, cost, budget, m, domain_cols, domain_levels, split_idx
    ))
  }

  results <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- targets[rows, , drop = FALSE]
    # Drop domain columns for inner solve
    sub_inner <- sub[, !names(sub) %in% domain_cols, drop = FALSE]
    if (multistage) {
      .n_multi_cluster(sub_inner, cost, budget, m)
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
      cost
    )
  } else {
    .aggregate_simple_domains(
      results,
      targets,
      domain_cols,
      domain_levels,
      split_idx
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
  split_idx
) {
  # Build domains data frame
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

  # Overall n = max across domains
  n_max <- max(domains$.n)
  overall_binding_idx <- which.max(domains$.n)
  binding_label <- domains$.binding[overall_binding_idx]

  .new_svyplan_n(
    n = n_max,
    type = "multi",
    targets = targets,
    detail = NULL,
    binding = binding_label,
    domains = domains
  )
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
  cost
) {
  stages <- results[[1L]]$stages

  domain_rows <- lapply(domain_levels, function(lev) {
    res <- results[[lev]]
    rows <- split_idx[[lev]]
    dom_vals <- targets[rows[1L], domain_cols, drop = FALSE]
    for (s in seq_len(stages)) {
      dom_vals[[paste0("n", s)]] <- res$n[s]
    }
    dom_vals$.total_n <- res$total_n
    dom_vals$.cv <- res$cv
    dom_vals$.cost <- res$cost
    dom_vals$.binding <- res$binding
    dom_vals
  })
  domains <- do.call(rbind, domain_rows)
  rownames(domains) <- NULL

  # Overall: aggregate totals
  total_n <- sum(ceiling(domains$.total_n))
  total_cost <- sum(domains$.cost)
  worst_cv_idx <- which.max(domains$.cv)

  n_vec <- vapply(seq_len(stages), function(s) {
    max(domains[[paste0("n", s)]])
  }, numeric(1L))
  names(n_vec) <- paste0("n", seq_len(stages))

  .new_svyplan_cluster(
    n = n_vec,
    stages = stages,
    total_n = total_n,
    cv = domains$.cv[worst_cv_idx],
    cost = total_cost,
    params = list(cost = cost),
    targets = targets,
    detail = NULL,
    binding = domains$.binding[worst_cv_idx],
    domains = domains
  )
}

.n_multi_domains_joint <- function(targets, cost, budget, m, domain_cols,
                                   domain_levels, split_idx) {
  stages <- length(cost)
  nd <- length(domain_levels)
  C1 <- cost[1L]
  C2 <- cost[2L]
  C3 <- if (stages == 3L) cost[3L] else NULL

  domain_params <- lapply(domain_levels, function(lev) {
    rows <- split_idx[[lev]]
    sub <- targets[rows, , drop = FALSE]
    labels <- if ("name" %in% names(sub)) sub$name else seq_along(rows)
    list(
      cv_t = sub$cv,
      delta1 = sub$delta1,
      delta2 = if ("delta2" %in% names(sub)) sub$delta2 else rep(0, length(rows)),
      rel_var = sub$rel_var,
      k1 = sub$k1,
      k2 = sub$k2,
      labels = labels
    )
  })

  eval_domain <- function(d, budget_d) {
    p <- domain_params[[d]]
    if (stages == 2L) {
      .eval_2stage_budget(p$cv_t, p$delta1, p$rel_var, p$k1,
                          C1, C2, budget_d, m)
    } else {
      .eval_3stage_budget(p$cv_t, p$delta1, p$delta2, p$rel_var, p$k1, p$k2,
                          C1, C2, C3, budget_d, m)
    }
  }

  outer_obj <- function(w) {
    w_last <- 1 - sum(w)
    if (w_last < 1e-4) return(1e12)
    fracs <- c(w, w_last)
    tryCatch(
      max(vapply(seq_len(nd), function(d) {
        eval_domain(d, fracs[d] * budget)$ratio
      }, numeric(1L))),
      error = function(e) 1e12
    )
  }

  if (nd == 2L) {
    opt <- optimize(
      function(w1) outer_obj(w1),
      interval = c(1e-4, 1 - 1e-4)
    )
    w_opt <- c(opt$minimum, 1 - opt$minimum)
  } else {
    init_w <- rep(1 / nd, nd - 1L)
    opt <- optim(
      par = init_w,
      fn = outer_obj,
      method = "L-BFGS-B",
      lower = rep(1e-4, nd - 1L),
      upper = rep(1 - 1e-4, nd - 1L)
    )
    w_opt <- c(opt$par, 1 - sum(opt$par))
  }

  budgets <- w_opt * budget
  results <- lapply(seq_len(nd), function(d) {
    bres <- eval_domain(d, budgets[d])
    p <- domain_params[[d]]
    if (stages == 2L) {
      n_vec <- c(n1 = bres$n1, n2 = bres$n2)
    } else {
      n_vec <- c(n1 = bres$n1, n2 = bres$n2, n3 = bres$n3)
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
    results, targets, domain_cols, domain_levels, split_idx, cost
  )
  res$params$joint <- TRUE
  res
}

.eval_2stage_budget <- function(cv_t, delta, rel_var, k, C1, C2, budget, m) {
  nr <- length(cv_t)

  if (is.null(m)) {
    obj_fn <- function(n2) {
      n1 <- budget / (C1 + C2 * n2)
      if (n1 <= 0) return(1e12)
      cv_ratios <- vapply(
        seq_len(nr),
        function(j) {
          sqrt(
            rel_var[j] * k[j] / (n1 * n2) * (1 + delta[j] * (n2 - 1))
          ) / cv_t[j]
        },
        numeric(1L)
      )
      max(cv_ratios)
    }

    upper <- max(10, budget / (C1 + C2))
    opt <- optimize(obj_fn, interval = c(1, upper))
    n2_opt <- opt$minimum
    n1_opt <- budget / (C1 + C2 * n2_opt)
  } else {
    n1_opt <- m
    n2_opt <- (budget - C1 * m) / (C2 * m)
    if (n2_opt <= 0) {
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
        rel_var[j] * k[j] / (n1_opt * n2_opt) * (1 + delta[j] * (n2_opt - 1))
      )
    },
    numeric(1L)
  )
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    n2 = n2_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = budget
  )
}

.eval_3stage_budget <- function(cv_t, delta1, delta2, rel_var, k1, k2,
                                C1, C2, C3, budget, m) {
  nr <- length(cv_t)

  cv_fn <- function(n1, n2, n3) {
    vapply(
      seq_len(nr),
      function(j) {
        sqrt(
          rel_var[j] / (n1 * n2 * n3) *
            (k1[j] * delta1[j] * n2 * n3 +
             k2[j] * (1 + delta2[j] * (n3 - 1)))
        )
      },
      numeric(1L)
    )
  }

  if (is.null(m)) {
    obj_fn <- function(par) {
      n2 <- par[1L]
      n3 <- par[2L]
      n1 <- budget / (C1 + C2 * n2 + C3 * n2 * n3)
      if (n1 <= 0) return(1e12)
      max(cv_fn(n1, n2, n3) / cv_t)
    }

    init_n3 <- max(2, sqrt(C2 / C3))
    init_n2 <- max(2, sqrt(C1 / C2))
    opt <- optim(
      par = c(init_n2, init_n3),
      fn = obj_fn,
      method = "L-BFGS-B",
      lower = c(1, 1),
      upper = c(1000, 1000)
    )

    n2_opt <- opt$par[1L]
    n3_opt <- opt$par[2L]
    n1_opt <- budget / (C1 + C2 * n2_opt + C3 * n2_opt * n3_opt)
    total_cost <- budget
  } else {
    n1_opt <- m
    obj_fn <- function(par) {
      n2 <- par[1L]
      n3 <- par[2L]
      cost_check <- C1 * m + C2 * m * n2 + C3 * m * n2 * n3
      if (cost_check > budget) return(1e12)
      max(cv_fn(m, n2, n3) / cv_t)
    }

    init_n3 <- max(2, sqrt(C2 / C3))
    init_n2 <- max(2, (budget / m - C1) / (C2 + C3 * init_n3))
    init_n2 <- max(2, init_n2)
    opt <- optim(
      par = c(init_n2, init_n3),
      fn = obj_fn,
      method = "L-BFGS-B",
      lower = c(1, 1),
      upper = c(1000, 1000)
    )

    n2_opt <- opt$par[1L]
    n3_opt <- opt$par[2L]
    total_cost <- C1 * m + C2 * m * n2_opt + C3 * m * n2_opt * n3_opt
  }

  cv_achieved <- cv_fn(n1_opt, n2_opt, n3_opt)
  ratios <- cv_achieved / cv_t
  binding_idx <- which.max(ratios)

  list(
    n1 = n1_opt,
    n2 = n2_opt,
    n3 = n3_opt,
    cv_achieved = cv_achieved,
    ratio = max(ratios),
    binding_idx = binding_idx,
    cost = total_cost
  )
}
