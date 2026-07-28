#' Internal controls for joint constrained allocation
#' @keywords internal
#' @noRd
.bethel_control <- function(
  tolerance = 1e-6,
  kkt_tolerance = 1e-6,
  optimizer_tolerance = 1e-8,
  max_iterations = 1000L,
  integer_repair_iterations = 100000L
) {
  check_scalar(tolerance, "tolerance")
  check_scalar(kkt_tolerance, "kkt_tolerance")
  check_scalar(optimizer_tolerance, "optimizer_tolerance")
  if (tolerance <= 0 || kkt_tolerance <= 0 || optimizer_tolerance <= 0) {
    stop("Bethel tolerances must be positive", call. = FALSE)
  }
  if (kkt_tolerance > tolerance) {
    stop(
      "'kkt_tolerance' must not exceed public feasibility 'tolerance'",
      call. = FALSE
    )
  }
  max_iterations <- check_count(max_iterations, "max_iterations")
  integer_repair_iterations <- check_count(
    integer_repair_iterations, "integer_repair_iterations"
  )
  list(
    tolerance = tolerance,
    kkt_tolerance = kkt_tolerance,
    optimizer_tolerance = optimizer_tolerance,
    max_iterations = max_iterations,
    integer_repair_iterations = integer_repair_iterations
  )
}

#' Validate a canonical Bethel problem
#' @keywords internal
#' @noRd
.new_bethel_problem <- function(
  A,
  B,
  bound,
  cost,
  lower,
  upper,
  stratum_ids,
  constraint_ids,
  constraint_meta = NULL,
  total = NULL,
  domain_N = NULL,
  frame = NULL,
  measures = NULL,
  targets = NULL,
  alpha = NULL,
  variance_model = "one_stage_wald"
) {
  # Zero columns are allowed: a budget-objective problem may carry no hard
  # precision targets at all.
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) == 0L ||
      anyNA(A) || any(!is.finite(A)) || any(A < 0)) {
    stop("'A' must be a finite non-negative numeric matrix with one row per stratum",
         call. = FALSE)
  }
  H <- nrow(A)
  K <- ncol(A)
  numeric_K <- function(x) {
    is.numeric(x) && length(x) == K && !anyNA(x) && all(is.finite(x))
  }
  numeric_H <- function(x) {
    is.numeric(x) && length(x) == H && !anyNA(x) && all(is.finite(x))
  }
  if (!numeric_K(B)) stop("'B' must be a finite numeric vector", call. = FALSE)
  if (!numeric_K(bound) || any(bound <= 0)) {
    stop("'bound' must contain positive finite values", call. = FALSE)
  }
  if (!numeric_H(cost) || any(cost <= 0)) {
    stop("'cost' must contain positive finite values", call. = FALSE)
  }
  if (!numeric_H(lower) || !numeric_H(upper) || any(lower <= 0) ||
      any(upper < lower)) {
    stop("allocation bounds must be finite, positive, and ordered",
         call. = FALSE)
  }
  if (!is.character(stratum_ids) || length(stratum_ids) != H ||
      anyNA(stratum_ids) || any(!nzchar(stratum_ids)) ||
      anyDuplicated(stratum_ids)) {
    stop("'stratum_ids' must contain unique non-empty values", call. = FALSE)
  }
  if (!is.character(constraint_ids) || length(constraint_ids) != K ||
      anyNA(constraint_ids) || any(!nzchar(constraint_ids)) ||
      anyDuplicated(constraint_ids)) {
    stop("'constraint_ids' must contain unique non-empty values", call. = FALSE)
  }
  if (any(colSums(A) <= 0)) {
    stop("every precision constraint must have a positive coefficient",
         call. = FALSE)
  }
  lower_int <- as.integer(ceiling(lower - 1e-9))
  upper_int <- as.integer(floor(upper + 1e-9))
  if (any(lower_int > upper_int)) {
    bad <- stratum_ids[lower_int > upper_int]
    stop(
      sprintf("no integer allocation satisfies the bounds for stratum: %s",
              paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }
  structure(
    list(
      A = A,
      B = as.numeric(B),
      bound = as.numeric(bound),
      cost = as.numeric(cost),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      lower_int = lower_int,
      upper_int = upper_int,
      stratum_ids = stratum_ids,
      constraint_ids = constraint_ids,
      constraint_meta = constraint_meta,
      total = total,
      domain_N = domain_N,
      frame = frame,
      measures = measures,
      targets = targets,
      alpha = alpha,
      variance_model = variance_model
    ),
    class = c("svyplan_bethel_problem", "list")
  )
}

#' Evaluate canonical reciprocal constraints
#' @keywords internal
#' @noRd
.bethel_constraint_value <- function(A, allocation) {
  drop(crossprod(A, 1 / allocation))
}

#' Initialize non-negative dual multipliers from single constraints
#' @keywords internal
#' @noRd
.bethel_initial_lambda <- function(A, cost, lower, upper, tolerance = 1e-8) {
  K <- ncol(A)
  out <- numeric(K)
  for (k in seq_len(K)) {
    ak <- A[, k]
    residual <- function(lambda) {
      w <- ak * lambda
      n <- pmin(pmax(sqrt(w / cost), lower), upper)
      sum(ak / n) - 1
    }
    if (residual(0) <= tolerance) next
    hi <- 1
    for (i in seq_len(100L)) {
      if (residual(hi) <= 0) break
      hi <- hi * 10
    }
    if (!is.finite(hi) || residual(hi) > tolerance) {
      out[k] <- hi
    } else {
      out[k] <- uniroot(residual, c(0, hi), tol = tolerance)$root
    }
  }
  out / max(1L, K)
}

#' Compute normalized KKT diagnostics for a Bethel solution
#' @keywords internal
#' @noRd
.bethel_kkt <- function(A, cost, lower, upper, allocation, lambda) {
  constraint_residual <- .bethel_constraint_value(A, allocation) - 1
  w <- drop(A %*% lambda)
  derivative <- cost - w / allocation^2
  scale <- max(1, max(abs(cost)))
  tol_bound <- 1e-7 * pmax(1, abs(allocation))
  fixed <- abs(upper - lower) <= tol_bound
  on_lower <- !fixed & abs(allocation - lower) <= tol_bound
  on_upper <- !fixed & abs(allocation - upper) <= tol_bound
  interior <- !(fixed | on_lower | on_upper)
  stationarity <- numeric(length(allocation))
  stationarity[on_lower] <- pmin(derivative[on_lower], 0)
  stationarity[on_upper] <- pmax(derivative[on_upper], 0)
  stationarity[interior] <- derivative[interior]
  projected <- lambda - pmax(0, lambda + constraint_residual)
  dual_value <- sum(cost * allocation + w / allocation) - sum(lambda)
  primal_value <- sum(cost * allocation)
  list(
    primal_residual = max(c(0, constraint_residual)),
    projected_dual_residual = max(abs(projected)),
    stationarity_residual = max(abs(stationarity)) / scale,
    complementarity_residual = max(
      abs(lambda * constraint_residual) / pmax(1, abs(lambda))
    ),
    relative_duality_gap = abs(primal_value - dual_value) /
      max(1, abs(primal_value)),
    constraint_residual = constraint_residual,
    active_precision = which(abs(constraint_residual) <= 1e-6),
    active_lower = which(abs(allocation - lower) <= tol_bound),
    active_upper = which(abs(allocation - upper) <= tol_bound),
    dual_value = dual_value
  )
}

#' Solve a bounded fixed-coefficient Bethel allocation problem
#' @keywords internal
#' @noRd
.bethel_solve <- function(
  A,
  bound,
  cost,
  lower,
  upper,
  control = .bethel_control(),
  constraint_ids = colnames(A) %||% as.character(seq_len(ncol(A)))
) {
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) == 0L || ncol(A) == 0L ||
      anyNA(A) || any(!is.finite(A)) || any(A < 0)) {
    stop("'A' must be a non-empty finite non-negative numeric matrix",
         call. = FALSE)
  }
  H <- nrow(A)
  K <- ncol(A)
  if (!is.numeric(bound) || length(bound) != K || anyNA(bound) ||
      any(!is.finite(bound)) || any(bound <= 0)) {
    stop("'bound' must contain one positive finite value per constraint",
         call. = FALSE)
  }
  for (z in list(cost = cost, lower = lower, upper = upper)) {
    if (!is.numeric(z) || length(z) != H || anyNA(z) || any(!is.finite(z))) {
      stop("cost and bounds must be finite vectors with nrow(A) values",
           call. = FALSE)
    }
  }
  if (any(cost <= 0) || any(lower <= 0) || any(upper < lower)) {
    stop("cost and bounds must be positive and ordered", call. = FALSE)
  }
  if (length(constraint_ids) != K) {
    stop("'constraint_ids' must have length ncol(A)", call. = FALSE)
  }
  tolerance <- control$tolerance
  at_upper <- .bethel_constraint_value(A, upper) / bound
  infeasible <- which(at_upper > 1 + tolerance)
  if (length(infeasible) > 0L) {
    return(list(
      allocation = NULL,
      classification = "infeasible",
      converged = FALSE,
      infeasible = constraint_ids[infeasible],
      message = "precision targets are unattainable at the upper allocation",
      feasibility_tolerance = control$tolerance,
      kkt_tolerance = control$kkt_tolerance
    ))
  }

  cost_scale <- stats::median(cost)
  As <- sweep(A, 2L, bound, "/")
  cs <- cost / cost_scale
  lower_value <- .bethel_constraint_value(As, lower)
  if (all(lower_value <= 1 + tolerance)) {
    lambda <- numeric(K)
    kkt <- .bethel_kkt(As, cs, lower, upper, lower, lambda)
    return(c(
      list(
        allocation = lower,
        cost = sum(cost * lower),
        lambda = lambda,
        iterations = 0L,
        convergence_code = 0L,
        converged = TRUE,
        classification = "optimal",
        message = "lower allocation satisfies every precision target",
        cost_scale = cost_scale,
        constraint_scale = bound,
        feasibility_tolerance = control$tolerance,
        kkt_tolerance = control$kkt_tolerance
      ),
      kkt
    ))
  }

  allocation_for <- function(lambda) {
    w <- drop(As %*% lambda)
    pmin(pmax(sqrt(w / cs), lower), upper)
  }
  objective <- function(lambda) {
    n <- allocation_for(lambda)
    w <- drop(As %*% lambda)
    -(sum(cs * n + w / n) - sum(lambda))
  }
  gradient <- function(lambda) {
    n <- allocation_for(lambda)
    -(drop(crossprod(As, 1 / n)) - 1)
  }
  lambda0 <- .bethel_initial_lambda(
    As, cs, lower, upper, tolerance = control$optimizer_tolerance
  )
  opt <- stats::optim(
    lambda0,
    objective,
    gradient,
    method = "L-BFGS-B",
    lower = rep(0, K),
    control = list(
      maxit = control$max_iterations,
      factr = 1e2,
      pgtol = min(control$optimizer_tolerance,
                  control$kkt_tolerance / 10)
    )
  )
  allocation <- allocation_for(opt$par)
  kkt <- .bethel_kkt(As, cs, lower, upper, allocation, opt$par)
  checks <- unlist(kkt[c(
    "primal_residual", "projected_dual_residual",
    "stationarity_residual", "complementarity_residual",
    "relative_duality_gap"
  )])
  # KKT certification is authoritative. L-BFGS-B can report an abnormal
  # line-search termination on the flat dual face even when the recovered
  # primal-dual pair satisfies every optimality condition.
  converged <- all(checks <= control$kkt_tolerance)
  lambda_original <- cost_scale * opt$par / bound
  normalized_columns <- as.data.frame(t(As), optional = TRUE)
  duplicate_multiplier <- duplicated(normalized_columns) |
    duplicated(normalized_columns, fromLast = TRUE)
  # Exact duplicate normalized constraints identify only their combined
  # multiplier. Do not expose an arbitrary per-row split as sensitivity.
  lambda_original[duplicate_multiplier] <- NA_real_
  c(
    list(
      allocation = allocation,
      cost = sum(cost * allocation),
      lambda = lambda_original,
      multiplier_identifiable = !duplicate_multiplier,
      lambda_scaled = opt$par,
      iterations = unname(opt$counts[["function"]]),
      convergence_code = opt$convergence,
      converged = converged,
      classification = if (converged) "optimal" else "failed",
      message = opt$message %||% if (converged) "converged" else
        "KKT certification failed",
      cost_scale = cost_scale,
      constraint_scale = bound,
      feasibility_tolerance = control$tolerance,
      kkt_tolerance = control$kkt_tolerance
    ),
    kkt
  )
}

#' Collision-free key for generalized allocation joins
#' @keywords internal
#' @noRd
.bethel_key <- function(...) {
  vals <- list(...)
  enc <- lapply(vals, function(x) {
    x <- as.character(x)
    if (anyNA(x)) stop("join keys must not contain missing values", call. = FALSE)
    paste0(nchar(x), "_", x)
  })
  do.call(paste, c(enc, list(sep = ":")))
}

#' Validate the fixed-take stage contract and create decision-unit inputs
#' @keywords internal
#' @noRd
.bethel_stage_spec <- function(
  frame,
  measures,
  unit_cost,
  N,
  max_weight,
  take_all,
  min_n_stratum
) {
  frame_names <- names(frame)
  measure_names <- if (is.data.frame(measures)) {
    names(measures)[vapply(
      measures,
      function(x) any(!is.na(x)),
      logical(1L)
    )]
  } else character(0)
  all_names <- union(frame_names, measure_names)
  stage_names <- c(
    "n_psu", "n_per_psu", "n_per_ssu", "icc_psu", "icc_ssu",
    "var_ratio_psu", "var_ratio_ssu", "cost_psu", "cost_ssu", "cost_tsu",
    "N_psu", "N_ssu"
  )
  present <- intersect(stage_names, all_names)
  if (length(present) == 0L) {
    bounds <- .alloc_bounds(N, max_weight, take_all, min_n_stratum)
    return(list(
      stages = 1L,
      take = rep(1, length(N)),
      cost = NULL,
      lower = bounds$m_h,
      upper = bounds$M_h,
      variance_model = "one_stage_wald"
    ))
  }
  if ("n_psu" %in% present) {
    stop(
      "'n_psu' is the free decision in generalized fixed-take allocation and must not be supplied",
      call. = FALSE
    )
  }
  misplaced <- intersect(
    c("n_per_psu", "n_per_ssu", "cost_psu", "cost_ssu", "cost_tsu",
      "N_psu", "N_ssu"),
    measure_names
  )
  if (length(misplaced) > 0L) {
    stop(
      sprintf("stage frame column(s) must be supplied in 'frame', not 'measures': %s",
              paste(sQuote(misplaced), collapse = ", ")),
      call. = FALSE
    )
  }
  stages <- if (length(intersect(
    present, c("n_per_ssu", "icc_ssu", "var_ratio_ssu", "cost_tsu", "N_ssu")
  )) > 0L) 3L else 2L
  required <- if (stages == 2L) {
    c("N_psu", "n_per_psu", "cost_psu", "cost_ssu")
  } else {
    c("N_psu", "N_ssu", "n_per_psu", "n_per_ssu",
      "cost_psu", "cost_ssu", "cost_tsu")
  }
  missing <- setdiff(required, frame_names)
  if (length(missing) > 0L) {
    stop(
      sprintf("fixed-take %d-stage allocation frame must contain: %s",
              stages, paste(sQuote(missing), collapse = ", ")),
      call. = FALSE
    )
  }
  if ("unit_cost" %in% frame_names || !is.null(unit_cost)) {
    stop(
      "fixed-take multistage allocation uses stage costs; do not supply 'unit_cost'",
      call. = FALSE
    )
  }
  if (any(take_all)) {
    stop(
      "'take_all' is not supported for fixed-take multistage allocation because taking every PSU does not imply an ultimate-unit census",
      call. = FALSE
    )
  }

  whole_positive <- function(x, name) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 1) ||
        any(abs(x - round(x)) > 1e-8)) {
      stop(sprintf("'%s' must contain positive whole numbers", name),
           call. = FALSE)
    }
    as.numeric(round(x))
  }
  positive <- function(x, name) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop(sprintf("'%s' must contain positive finite values", name),
           call. = FALSE)
    }
    as.numeric(x)
  }
  N_psu <- whole_positive(frame$N_psu, "N_psu")
  n_per_psu <- whole_positive(frame$n_per_psu, "n_per_psu")
  if (any(N_psu > N + 1e-8)) {
    stop("'N_psu' must not exceed the ultimate-unit population 'N'",
         call. = FALSE)
  }
  cost_psu <- positive(frame$cost_psu, "cost_psu")
  cost_ssu <- positive(frame$cost_ssu, "cost_ssu")

  if (stages == 2L) {
    take <- n_per_psu
    upper <- pmin(N_psu, N / take)
    cost <- cost_psu + cost_ssu * n_per_psu
    out <- list(
      stages = stages,
      take = take,
      N_psu = N_psu,
      n_per_psu = n_per_psu,
      cost_psu = cost_psu,
      cost_ssu = cost_ssu,
      cost = cost,
      upper = upper,
      variance_model = "two_stage_fixed_take_wald"
    )
  } else {
    N_ssu <- whole_positive(frame$N_ssu, "N_ssu")
    n_per_ssu <- whole_positive(frame$n_per_ssu, "n_per_ssu")
    if (any(N_ssu < N_psu) || any(N_ssu > N + 1e-8)) {
      stop("'N_ssu' must lie between 'N_psu' and 'N'", call. = FALSE)
    }
    take <- n_per_psu * n_per_ssu
    upper <- pmin(N_psu, N_ssu / n_per_psu, N / take)
    cost_tsu <- positive(frame$cost_tsu, "cost_tsu")
    cost <- cost_psu + cost_ssu * n_per_psu + cost_tsu * take
    out <- list(
      stages = stages,
      take = take,
      N_psu = N_psu,
      N_ssu = N_ssu,
      n_per_psu = n_per_psu,
      n_per_ssu = n_per_ssu,
      cost_psu = cost_psu,
      cost_ssu = cost_ssu,
      cost_tsu = cost_tsu,
      cost = cost,
      upper = upper,
      variance_model = "three_stage_fixed_take_wald"
    )
  }
  if (any(floor(upper + 1e-9) < 1)) {
    stop("fixed takes exceed the available stage populations", call. = FALSE)
  }
  ultimate_bounds <- .alloc_bounds(
    N, max_weight, rep(FALSE, length(N)), min_n_stratum
  )
  out$lower <- pmax(1, ultimate_bounds$m_h / take)
  if (any(out$lower > out$upper + 1e-8)) {
    stop(
      "constraints are infeasible after converting ultimate-unit bounds to whole-PSU bounds",
      call. = FALSE
    )
  }
  out
}

#' Normalize an indicator-specific stage parameter
#' @keywords internal
#' @noRd
.bethel_stage_measure <- function(
  name,
  measures,
  frame,
  measure_stratum,
  frame_stratum,
  required = FALSE,
  default = NA_real_,
  interval = NULL,
  positive = FALSE
) {
  nr <- nrow(measures)
  value <- if (name %in% names(measures)) measures[[name]] else
    rep(NA_real_, nr)
  if (!is.numeric(value)) {
    stop(sprintf("'measures$%s' must be numeric", name), call. = FALSE)
  }
  if (name %in% names(frame)) {
    fallback <- frame[[name]]
    if (!is.numeric(fallback)) {
      stop(sprintf("'frame$%s' must be numeric", name), call. = FALSE)
    }
    fallback <- fallback[match(measure_stratum, frame_stratum)]
    value[is.na(value)] <- fallback[is.na(value)]
  }
  value[is.na(value)] <- default
  if (required && anyNA(value)) {
    stop(sprintf("'%s' is required in 'measures' or as a stratum default in 'frame'",
                 name), call. = FALSE)
  }
  bad <- !is.na(value) & !is.finite(value)
  if (!is.null(interval)) {
    bad <- bad | (!is.na(value) &
      (value < interval[1L] | value > interval[2L]))
  }
  if (positive) bad <- bad | (!is.na(value) & value <= 0)
  if (any(bad)) {
    rule <- if (!is.null(interval)) {
      sprintf("values in [%s, %s]", interval[1L], interval[2L])
    } else if (positive) "positive finite values" else "finite values"
    stop(sprintf("'%s' must contain %s", name, rule), call. = FALSE)
  }
  value
}

#' Resolve frame rows for each indicator-domain requirement
#'
#' Shared by hard precision targets and objective components so that both
#' select frame strata by exactly the same rule.
#' @keywords internal
#' @noRd
.bethel_domain_indices <- function(frame, domain, level, what = "target") {
  H <- nrow(frame)
  K <- length(domain)
  G <- matrix(FALSE, H, K)
  indices <- vector("list", K)
  for (k in seq_len(K)) {
    if (domain[k] == ".overall") {
      idx <- seq_len(H)
    } else {
      if (!domain[k] %in% names(frame)) {
        stop(sprintf("domain column '%s' not found in frame", domain[k]),
             call. = FALSE)
      }
      frame_level <- frame[[domain[k]]]
      if (anyNA(frame_level)) {
        stop(sprintf("domain column '%s' must not contain missing values",
                     domain[k]), call. = FALSE)
      }
      idx <- which(as.character(frame_level) == level[k])
      if (length(idx) == 0L) {
        stop(sprintf("%s domain '%s=%s' is empty", what, domain[k], level[k]),
             call. = FALSE)
      }
    }
    indices[[k]] <- idx
    G[idx, k] <- TRUE
  }
  list(membership = G, indices = indices)
}

#' Normalize objective components and their frame-domain membership
#'
#' Accepts a character vector of indicator names (overall domain, equal
#' priority) or a long data frame using the same name/domain/level identifiers
#' as `targets` plus a non-negative `priority`.
#' @keywords internal
#' @noRd
.bethel_objective_spec <- function(frame, objective) {
  if (is.character(objective) || is.factor(objective)) {
    objective <- data.frame(
      name = as.character(objective),
      domain = ".overall",
      level = NA_character_,
      priority = 1,
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(objective) || nrow(objective) == 0L ||
      !"name" %in% names(objective)) {
    stop(
      "'objective' must be a character vector of indicator names or a non-empty data frame with a 'name' column",
      call. = FALSE
    )
  }
  o_name <- as.character(objective$name)
  if (anyNA(o_name) || any(!nzchar(o_name))) {
    stop("'objective$name' must contain non-empty values", call. = FALSE)
  }
  if ("weight" %in% names(objective) && !"priority" %in% names(objective)) {
    stop(
      "objective components are prioritized with 'priority'; 'weight' denotes the sampling weight N / n",
      call. = FALSE
    )
  }
  if (any(c("cv", "moe") %in% names(objective))) {
    stop(
      "objective components are always minimized in relative-variance units; put 'cv' or 'moe' requirements in 'targets'",
      call. = FALSE
    )
  }
  has_domain <- "domain" %in% names(objective)
  has_level <- "level" %in% names(objective)
  if (xor(has_domain, has_level)) {
    stop("'objective$domain' and 'objective$level' must be supplied together",
         call. = FALSE)
  }
  if (!has_domain) {
    domain <- rep(".overall", nrow(objective))
    level <- rep(NA_character_, nrow(objective))
  } else {
    domain <- as.character(objective$domain)
    level <- as.character(objective$level)
    if (anyNA(domain) || any(!nzchar(domain))) {
      stop("'objective$domain' must contain non-empty values", call. = FALSE)
    }
    overall <- domain == ".overall"
    if (any(overall & !is.na(level)) || any(!overall & is.na(level))) {
      stop("overall objective components require level = NA; other domains require a level",
           call. = FALSE)
    }
  }
  priority <- if ("priority" %in% names(objective)) objective$priority else
    rep(1, nrow(objective))
  if (!is.numeric(priority) || anyNA(priority) ||
      any(!is.finite(priority)) || any(priority < 0)) {
    stop("'objective$priority' must contain non-negative finite values",
         call. = FALSE)
  }
  if (sum(priority) <= 0) {
    stop("at least one objective 'priority' must be positive", call. = FALSE)
  }
  component_key <- .bethel_key(
    o_name, domain, ifelse(is.na(level), "<overall>", level)
  )
  if (anyDuplicated(component_key)) {
    stop("duplicate indicator-domain objective components are not allowed",
         call. = FALSE)
  }
  if ("component" %in% names(objective)) {
    component <- as.character(objective$component)
    if (anyNA(component) || any(!nzchar(component)) ||
        anyDuplicated(component)) {
      stop("'objective$component' must contain unique non-empty values",
           call. = FALSE)
    }
  } else {
    component <- paste0(
      o_name, "@", ifelse(domain == ".overall", ".overall",
                          paste0(domain, "=", level))
    )
    if (anyDuplicated(component)) {
      stop("generated objective component identifiers collide; supply 'component'",
           call. = FALSE)
    }
  }
  geometry <- .bethel_domain_indices(
    frame, domain, level, what = "objective component"
  )
  objective$component <- component
  objective$domain <- domain
  objective$level <- level
  objective$priority <- as.numeric(priority)
  list(
    name = o_name,
    domain = domain,
    level = level,
    priority = as.numeric(priority),
    component = component,
    membership = geometry$membership,
    indices = geometry$indices,
    objective = objective
  )
}

#' Normalize target rows and their frame-domain membership
#' @keywords internal
#' @noRd
.bethel_target_spec <- function(frame, targets, alpha, allow_empty = FALSE) {
  if (allow_empty &&
      (is.null(targets) || (is.data.frame(targets) && nrow(targets) == 0L))) {
    return(list(
      name = character(0),
      metric = character(0),
      target = numeric(0),
      domain = character(0),
      level = character(0),
      alpha = numeric(0),
      constraint = character(0),
      membership = matrix(FALSE, nrow(frame), 0L),
      indices = list(),
      targets = data.frame(
        name = character(0), domain = character(0), level = character(0),
        cv = numeric(0), moe = numeric(0), alpha = numeric(0),
        constraint = character(0), stringsAsFactors = FALSE
      )
    ))
  }
  if (!is.data.frame(targets) || nrow(targets) == 0L ||
      !"name" %in% names(targets)) {
    stop("'targets' must be a non-empty data frame with a 'name' column",
         call. = FALSE)
  }
  t_name <- as.character(targets$name)
  if (anyNA(t_name) || any(!nzchar(t_name))) {
    stop("'targets$name' must contain non-empty values", call. = FALSE)
  }
  t_cv <- if ("cv" %in% names(targets)) targets$cv else
    rep(NA_real_, nrow(targets))
  t_moe <- if ("moe" %in% names(targets)) targets$moe else
    rep(NA_real_, nrow(targets))
  if (!is.numeric(t_cv) || !is.numeric(t_moe)) {
    stop("target 'cv' and 'moe' columns must be numeric", call. = FALSE)
  }
  has_cv <- !is.na(t_cv)
  has_moe <- !is.na(t_moe)
  if (any(has_cv == has_moe) ||
      any(has_cv & (!is.finite(t_cv) | t_cv <= 0)) ||
      any(has_moe & (!is.finite(t_moe) | t_moe <= 0))) {
    stop("each target row must specify exactly one positive finite 'cv' or 'moe'",
         call. = FALSE)
  }
  metric <- ifelse(has_cv, "cv", "moe")
  target_value <- ifelse(has_cv, t_cv, t_moe)
  has_domain <- "domain" %in% names(targets)
  has_level <- "level" %in% names(targets)
  if (xor(has_domain, has_level)) {
    stop("'targets$domain' and 'targets$level' must be supplied together",
         call. = FALSE)
  }
  if (!has_domain) {
    domain <- rep(".overall", nrow(targets))
    level <- rep(NA_character_, nrow(targets))
  } else {
    domain <- as.character(targets$domain)
    level <- as.character(targets$level)
    if (anyNA(domain) || any(!nzchar(domain))) {
      stop("'targets$domain' must contain non-empty values", call. = FALSE)
    }
    overall <- domain == ".overall"
    if (any(overall & !is.na(level)) || any(!overall & is.na(level))) {
      stop("overall targets require level = NA; other domains require a level",
           call. = FALSE)
    }
  }
  t_alpha <- if ("alpha" %in% names(targets)) targets$alpha else
    rep(NA_real_, nrow(targets))
  if (!is.numeric(t_alpha)) {
    stop("'targets$alpha' must be numeric", call. = FALSE)
  }
  t_alpha[is.na(t_alpha)] <- alpha
  if (any(!is.finite(t_alpha) | t_alpha <= 0 | t_alpha >= 1)) {
    stop("'targets$alpha' must contain values in (0, 1)", call. = FALSE)
  }
  requirement_key <- .bethel_key(
    t_name, domain, ifelse(is.na(level), "<overall>", level), metric
  )
  if (anyDuplicated(requirement_key)) {
    stop("duplicate indicator-domain precision requirements are not allowed",
         call. = FALSE)
  }
  if ("constraint" %in% names(targets)) {
    constraint <- as.character(targets$constraint)
    if (anyNA(constraint) || any(!nzchar(constraint)) ||
        anyDuplicated(constraint)) {
      stop("'targets$constraint' must contain unique non-empty values",
           call. = FALSE)
    }
  } else {
    constraint <- paste0(
      t_name, "@", ifelse(domain == ".overall", ".overall",
                           paste0(domain, "=", level)), ":", metric
    )
    if (anyDuplicated(constraint)) {
      stop("generated constraint identifiers collide; supply 'constraint'",
           call. = FALSE)
    }
  }

  geometry <- .bethel_domain_indices(frame, domain, level, what = "target")
  targets$constraint <- constraint
  targets$domain <- domain
  targets$level <- level
  targets$alpha <- t_alpha
  list(
    name = t_name,
    metric = metric,
    target = target_value,
    domain = domain,
    level = level,
    alpha = t_alpha,
    constraint = constraint,
    membership = geometry$membership,
    indices = geometry$indices,
    targets = targets
  )
}

#' Build a generalized allocation problem
#' @keywords internal
#' @noRd
.build_bethel_problem <- function(
  frame,
  measures,
  targets,
  unit_cost = NULL,
  deff = 1,
  resp_rate = 1,
  alpha = 0.05,
  min_n_stratum = NULL,
  objective = NULL
) {
  check_deff(deff)
  check_resp_rate(resp_rate)
  check_alpha(alpha)
  if (!is.null(min_n_stratum)) check_scalar(min_n_stratum, "min_n_stratum")
  if (!is.data.frame(frame) || nrow(frame) == 0L) {
    stop("'frame' must be a non-empty data frame", call. = FALSE)
  }
  required_frame <- c("stratum", "N")
  missing_frame <- setdiff(required_frame, names(frame))
  if (length(missing_frame) > 0L) {
    stop(sprintf("generalized allocation frame must contain: %s",
                 paste(sQuote(missing_frame), collapse = ", ")),
         call. = FALSE)
  }
  stratum <- as.character(frame$stratum)
  if (anyNA(stratum) || any(!nzchar(stratum)) || anyDuplicated(stratum)) {
    stop("'frame$stratum' must contain unique non-empty values", call. = FALSE)
  }
  N <- frame$N
  if (!is.numeric(N) || anyNA(N) || any(!is.finite(N)) || any(N <= 0)) {
    stop("'N' must contain positive finite values", call. = FALSE)
  }
  if ("cost" %in% names(frame)) {
    stop(
      "the per-stratum cost column is 'unit_cost'; 'cost' is the total field cost",
      call. = FALSE
    )
  }
  max_weight <- if ("max_weight" %in% names(frame)) {
    frame$max_weight
  } else rep(NA_real_, nrow(frame))
  if (!is.numeric(max_weight)) {
    stop("'max_weight' must be numeric", call. = FALSE)
  }
  bad_weight <- !is.na(max_weight) &
    (!is.finite(max_weight) | max_weight < 1)
  if (any(bad_weight)) {
    stop("'max_weight' must be >= 1 and finite when provided", call. = FALSE)
  }
  take_all <- .check_take_all(frame[["take_all"]], nrow(frame))

  # Determine target relevance before validating measure values. Unused
  # indicator-stratum rows may carry incomplete planning information and must
  # not change the mathematical problem.
  if (!is.data.frame(measures) || nrow(measures) == 0L) {
    stop("'measures' must be a non-empty data frame", call. = FALSE)
  }
  if (!all(c("stratum", "name") %in% names(measures))) {
    stop("'measures' must contain 'stratum' and 'name' columns", call. = FALSE)
  }
  m_stratum_all <- as.character(measures$stratum)
  m_name_all <- as.character(measures$name)
  if (anyNA(m_stratum_all) || any(!nzchar(m_stratum_all)) ||
      anyNA(m_name_all) || any(!nzchar(m_name_all))) {
    stop("measure keys must contain non-missing, non-empty values",
         call. = FALSE)
  }
  unknown_strata <- setdiff(unique(m_stratum_all), stratum)
  if (length(unknown_strata) > 0L) {
    stop(sprintf("measure strata not found in frame: %s",
                 paste(unknown_strata, collapse = ", ")), call. = FALSE)
  }
  measure_key_all <- .bethel_key(m_stratum_all, m_name_all)
  if (anyDuplicated(measure_key_all)) {
    stop("'measures' must have unique stratum x name rows", call. = FALSE)
  }
  target_spec <- .bethel_target_spec(
    frame, targets, alpha, allow_empty = !is.null(objective)
  )
  objective_spec <- if (is.null(objective)) NULL else
    .bethel_objective_spec(frame, objective)
  # Targets and objective components share one coefficient-construction path;
  # the columns are split apart again once the requirement loop has run.
  req_indices <- c(target_spec$indices, objective_spec$indices)
  req_name <- c(target_spec$name, objective_spec$name)
  required_keys <- unique(unlist(Map(
    function(idx, name) {
      .bethel_key(stratum[idx], rep(name, length(idx)))
    },
    req_indices,
    req_name
  ), use.names = FALSE))
  missing_keys <- setdiff(required_keys, measure_key_all)
  if (length(missing_keys) > 0L) {
    for (k in seq_along(req_indices)) {
      idx <- req_indices[[k]]
      wanted <- .bethel_key(stratum[idx], rep(req_name[k], length(idx)))
      missing <- stratum[idx][!wanted %in% measure_key_all]
      if (length(missing) > 0L) {
        stop(
          sprintf("missing measures for indicator '%s' in strata: %s",
                  req_name[k], paste(missing, collapse = ", ")),
          call. = FALSE
        )
      }
    }
  }
  measures <- measures[measure_key_all %in% required_keys, , drop = FALSE]
  m_stratum <- as.character(measures$stratum)
  m_name <- as.character(measures$name)
  measure_key <- .bethel_key(m_stratum, m_name)

  stage <- .bethel_stage_spec(
    frame, measures, unit_cost, N, max_weight, take_all, min_n_stratum
  )

  if (stage$stages > 1L) {
    cost <- stage$cost
  } else if (!is.null(unit_cost)) {
    if (!is.numeric(unit_cost) || anyNA(unit_cost) ||
        any(!is.finite(unit_cost)) || any(unit_cost <= 0) ||
        !length(unit_cost) %in% c(1L, nrow(frame))) {
      stop("'unit_cost' must contain one or nrow(frame) positive finite values",
           call. = FALSE)
    }
    cost <- rep_len(unit_cost, nrow(frame))
  } else if ("unit_cost" %in% names(frame)) {
    cost <- frame$unit_cost
    if (!is.numeric(cost) || anyNA(cost) || any(!is.finite(cost)) ||
        any(cost <= 0)) {
      stop("'unit_cost' must contain positive finite values", call. = FALSE)
    }
  } else {
    cost <- rep(1, nrow(frame))
  }
  mcol <- function(name) {
    if (name %in% names(measures)) measures[[name]] else
      rep(NA_real_, nrow(measures))
  }
  p <- mcol("p")
  mean_value <- mcol("mean")
  var_value <- mcol("var")
  sd_value <- mcol("sd")
  for (x in list(p = p, mean = mean_value, var = var_value, sd = sd_value)) {
    if (!is.numeric(x)) {
      stop("measure statistics must be numeric", call. = FALSE)
    }
  }
  has_p <- !is.na(p)
  has_mean <- !is.na(mean_value)
  has_var <- !is.na(var_value)
  has_sd <- !is.na(sd_value)
  valid_prop <- has_p & !has_mean & !has_var & !has_sd
  valid_mean <- !has_p & has_mean & xor(has_var, has_sd)
  if (any(!(valid_prop | valid_mean))) {
    stop(
      "each measure row must specify either 'p', or 'mean' and exactly one of 'var'/'sd'",
      call. = FALSE
    )
  }
  if (any(has_p & (!is.finite(p) | p < 0 | p > 1))) {
    stop("'p' must contain values in [0, 1]", call. = FALSE)
  }
  if (any(has_mean & !is.finite(mean_value))) {
    stop("'mean' must contain finite values", call. = FALSE)
  }
  if (any(has_var & (!is.finite(var_value) | var_value < 0)) ||
      any(has_sd & (!is.finite(sd_value) | sd_value < 0))) {
    stop("'var' and 'sd' must contain non-negative finite values",
         call. = FALSE)
  }
  mean_norm <- ifelse(has_p, p, mean_value)
  var_norm <- ifelse(has_p, p * (1 - p),
                     ifelse(has_var, var_value, sd_value^2))
  row_deff <- if ("deff" %in% names(measures)) measures$deff else
    rep(NA_real_, nrow(measures))
  row_resp <- if ("resp_rate" %in% names(measures)) measures$resp_rate else
    rep(NA_real_, nrow(measures))
  if (!is.numeric(row_deff) ||
      any(!is.na(row_deff) & (!is.finite(row_deff) | row_deff <= 0))) {
    stop("'measures$deff' must contain positive finite values or NA",
         call. = FALSE)
  }
  if (!is.numeric(row_resp) ||
      any(!is.na(row_resp) &
          (!is.finite(row_resp) | row_resp <= 0 | row_resp > 1))) {
    stop("'measures$resp_rate' must contain values in (0, 1] or NA",
         call. = FALSE)
  }
  deff_norm <- ifelse(is.na(row_deff), deff, row_deff)
  resp_norm <- ifelse(is.na(row_resp), resp_rate, row_resp)
  if (stage$stages >= 2L) {
    icc_psu_norm <- .bethel_stage_measure(
      "icc_psu", measures, frame, m_stratum, stratum,
      required = TRUE, interval = c(0, 1)
    )
    var_ratio_psu_norm <- .bethel_stage_measure(
      "var_ratio_psu", measures, frame, m_stratum, stratum,
      default = 1, positive = TRUE
    )
  }
  if (stage$stages == 3L) {
    icc_ssu_norm <- .bethel_stage_measure(
      "icc_ssu", measures, frame, m_stratum, stratum,
      required = TRUE, interval = c(0, 1)
    )
    # var_ratio_ssu is fixed by var_ratio_psu * icc_psu + var_ratio_ssu = 1; see .var_ratio_ssu_default().
    var_ratio_ssu_norm <- .bethel_stage_measure(
      "var_ratio_ssu", measures, frame, m_stratum, stratum,
      default = NA_real_, positive = TRUE
    )
    implied_var_ratio_ssu <- .var_ratio_ssu_default(var_ratio_psu_norm, icc_psu_norm)
    supplied <- !is.na(var_ratio_ssu_norm)
    var_ratio_ssu_norm[!supplied] <- implied_var_ratio_ssu[!supplied]
  }

  t_name <- target_spec$name
  metric <- target_spec$metric
  target_value <- target_spec$target
  domain <- target_spec$domain
  level <- target_spec$level
  t_alpha <- target_spec$alpha
  constraint <- target_spec$constraint
  H <- nrow(frame)
  K <- length(t_name)
  J <- length(objective_spec$name)
  R <- K + J
  membership <- cbind(target_spec$membership, objective_spec$membership)
  A_all <- matrix(0, H, R)
  B_all <- numeric(R)
  total_all <- numeric(R)
  domain_N_all <- numeric(R)
  mean_all <- matrix(NA_real_, H, R)
  var_all <- matrix(NA_real_, H, R)
  deff_all <- matrix(NA_real_, H, R)
  resp_all <- matrix(NA_real_, H, R)
  for (k in seq_len(R)) {
    idx <- req_indices[[k]]
    wanted <- .bethel_key(stratum[idx], rep(req_name[k], length(idx)))
    mi <- match(wanted, measure_key)
    if (anyNA(mi)) {
      stop("internal error matching required measure rows", call. = FALSE)
    }
    mean_all[idx, k] <- mean_norm[mi]
    var_all[idx, k] <- var_norm[mi]
    deff_all[idx, k] <- deff_norm[mi]
    resp_all[idx, k] <- resp_norm[mi]
    variance_factor <- rep(1, length(idx))
    decision_take <- rep(1, length(idx))
    if (stage$stages == 2L) {
      m <- stage$n_per_psu[idx]
      variance_factor <- var_ratio_psu_norm[mi] *
        (1 + icc_psu_norm[mi] * (m - 1))
      decision_take <- m
    } else if (stage$stages == 3L) {
      m <- stage$n_per_psu[idx]
      q <- stage$n_per_ssu[idx]
      variance_factor <-
        var_ratio_psu_norm[mi] * icc_psu_norm[mi] * m * q +
        var_ratio_ssu_norm[mi] * (1 + icc_ssu_norm[mi] * (q - 1))
      decision_take <- m * q
    }
    variance_inflated <- var_norm[mi] * variance_factor
    A_all[idx, k] <- N[idx]^2 * variance_inflated * deff_norm[mi] /
      (resp_norm[mi] * decision_take)
    B_all[k] <- -sum(N[idx] * variance_inflated * deff_norm[mi])
    total_all[k] <- sum(N[idx] * mean_norm[mi])
    domain_N_all[k] <- sum(N[idx])
  }
  no_variance <- colSums(A_all) <= 0
  if (any(no_variance[seq_len(K)])) {
    stop(sprintf("precision constraints have no positive variance: %s",
                 paste(constraint[no_variance[seq_len(K)]], collapse = ", ")),
         call. = FALSE)
  }
  if (J > 0L && any(no_variance[K + seq_len(J)])) {
    stop(sprintf("objective components have no positive variance: %s",
                 paste(objective_spec$component[no_variance[K + seq_len(J)]],
                       collapse = ", ")),
         call. = FALSE)
  }
  # Relative-variance requirements are undefined against a negligible total.
  # Both CV constraints and objective components use this guard.
  negligible_total <- function(k) {
    scale_total <- sum(abs(N[membership[, k]] * mean_all[membership[, k], k]))
    abs(total_all[k]) <= sqrt(.Machine$double.eps) * max(1, scale_total)
  }
  A <- A_all[, seq_len(K), drop = FALSE]
  B <- B_all[seq_len(K)]
  total <- total_all[seq_len(K)]
  domain_N <- domain_N_all[seq_len(K)]
  mean_hk <- mean_all[, seq_len(K), drop = FALSE]
  var_hk <- var_all[, seq_len(K), drop = FALSE]
  deff_hk <- deff_all[, seq_len(K), drop = FALSE]
  resp_hk <- resp_all[, seq_len(K), drop = FALSE]
  G <- membership[, seq_len(K), drop = FALSE]
  Vmax <- numeric(K)
  for (k in seq_len(K)) {
    if (metric[k] == "cv") {
      if (negligible_total(k)) {
        stop(sprintf("CV is undefined for constraint '%s' with zero or negligible total",
                     constraint[k]), call. = FALSE)
      }
      Vmax[k] <- (target_value[k] * abs(total[k]))^2
    } else {
      z <- stats::qnorm(1 - t_alpha[k] / 2)
      Vmax[k] <- (domain_N[k] * target_value[k] / z)^2
    }
  }
  bound <- Vmax - B
  objective_part <- NULL
  if (J > 0L) {
    for (j in seq_len(J)) {
      if (negligible_total(K + j)) {
        stop(sprintf("the objective is undefined for component '%s' with zero or negligible total",
                     objective_spec$component[j]), call. = FALSE)
      }
    }
    A_obj <- A_all[, K + seq_len(J), drop = FALSE]
    B_obj <- B_all[K + seq_len(J)]
    total_obj <- total_all[K + seq_len(J)]
    priority <- objective_spec$priority
    scaled <- priority / total_obj^2
    objective_part <- list(
      component = objective_spec$component,
      name = objective_spec$name,
      domain = objective_spec$domain,
      level = objective_spec$level,
      priority = priority,
      A = A_obj,
      B = B_obj,
      total = total_obj,
      domain_N = domain_N_all[K + seq_len(J)],
      q_h = as.numeric(A_obj %*% scaled),
      Q0 = sum(scaled * B_obj),
      spec = objective_spec$objective
    )
    if (sum(objective_part$q_h) <= 0) {
      stop("the objective has no positive-priority component with positive variance",
           call. = FALSE)
    }
  }
  meta <- data.frame(
    constraint = constraint,
    name = t_name,
    domain = domain,
    level = level,
    .metric = metric,
    .target = target_value,
    stringsAsFactors = FALSE
  )
  problem <- .new_bethel_problem(
    A = A,
    B = B,
    bound = bound,
    cost = cost,
    lower = stage$lower,
    upper = stage$upper,
    stratum_ids = stratum,
    constraint_ids = constraint,
    constraint_meta = meta,
    total = total,
    domain_N = domain_N,
    frame = frame,
    measures = measures,
    targets = target_spec$targets,
    alpha = t_alpha,
    variance_model = stage$variance_model
  )
  problem$objective <- objective_part
  problem$membership <- G
  problem$mean_hk <- mean_hk
  problem$var_hk <- var_hk
  problem$deff_hk <- deff_hk
  problem$resp_hk <- resp_hk
  problem$take_all <- take_all
  problem$population_N <- N
  problem$stage <- stage
  problem$stages <- stage$stages
  problem
}

#' Normalize and validate an allocation vector for a Bethel problem
#' @keywords internal
#' @noRd
.bethel_match_allocation <- function(
  problem,
  allocation,
  tolerance = .bethel_control()$tolerance
) {
  if (!is.numeric(allocation) || length(allocation) != length(problem$cost) ||
      anyNA(allocation) || any(!is.finite(allocation))) {
    stop("'n' must be a finite numeric vector with one value per stratum",
         call. = FALSE)
  }
  if (!is.null(names(allocation))) {
    nm <- names(allocation)
    if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
        !setequal(nm, problem$stratum_ids)) {
      stop("named 'n' must match every frame stratum exactly", call. = FALSE)
    }
    allocation <- allocation[match(problem$stratum_ids, nm)]
  }
  if (any(allocation <= 0)) stop("all 'n' elements must be positive", call. = FALSE)
  if (any(allocation > problem$upper + tolerance)) {
    bad <- problem$stratum_ids[allocation > problem$upper + tolerance]
    stop(sprintf("'n' exceeds the population bound for stratum: %s",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  as.numeric(allocation)
}

#' Convert a public ultimate-unit allocation to the canonical decision unit
#' @keywords internal
#' @noRd
.bethel_to_decision <- function(
  problem,
  allocation,
  require_whole_stages = FALSE,
  tolerance = 1e-8
) {
  if (!is.numeric(allocation) || length(allocation) != length(problem$cost) ||
      anyNA(allocation) || any(!is.finite(allocation))) {
    stop("'n' must be a finite numeric vector with one value per stratum",
         call. = FALSE)
  }
  if (!is.null(names(allocation))) {
    nm <- names(allocation)
    if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
        !setequal(nm, problem$stratum_ids)) {
      stop("named 'n' must match every frame stratum exactly", call. = FALSE)
    }
    allocation <- allocation[match(problem$stratum_ids, nm)]
  }
  if (any(allocation <= 0)) stop("all 'n' elements must be positive", call. = FALSE)
  public_upper <- problem$upper * problem$stage$take
  if (any(allocation > public_upper + tolerance)) {
    bad <- problem$stratum_ids[allocation > public_upper + tolerance]
    stop(sprintf("'n' exceeds the population bound for stratum: %s",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  decision <- as.numeric(allocation) / problem$stage$take
  nonwhole <- abs(decision - round(decision)) > tolerance
  if (problem$stages > 1L && require_whole_stages && any(nonwhole)) {
    stop(
      sprintf(
        "operational 'n' must correspond to a whole number of PSUs for stratum: %s",
        paste(problem$stratum_ids[nonwhole], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  decision
}

#' Convert a canonical decision allocation to ultimate-unit sample sizes
#' @keywords internal
#' @noRd
.bethel_to_ultimate <- function(problem, allocation) {
  as.numeric(allocation) * problem$stage$take
}

#' Evaluate achieved precision for an allocation
#' @keywords internal
#' @noRd
.precision_from_allocation <- function(
  problem,
  allocation,
  lambda = NULL,
  tolerance = .bethel_control()$tolerance
) {
  allocation <- .bethel_match_allocation(
    problem, allocation, tolerance = tolerance
  )
  reciprocal <- .bethel_constraint_value(problem$A, allocation)
  variance <- problem$B + reciprocal
  cancellation_scale <- pmax(1, abs(problem$B), abs(reciprocal))
  tiny_negative <- variance < 0 &
    abs(variance) <= 100 * .Machine$double.eps * cancellation_scale
  variance[tiny_negative] <- 0
  if (any(variance < 0)) {
    stop("internal error: allocation produced a materially negative variance",
         call. = FALSE)
  }
  se_total <- sqrt(variance)
  se <- se_total / problem$domain_N
  moe <- stats::qnorm(1 - problem$alpha / 2) * se
  cv <- se_total / abs(problem$total)
  metric <- problem$constraint_meta$.metric
  achieved <- as.numeric(ifelse(metric == "cv", cv, moe))
  target <- problem$constraint_meta$.target
  ratio <- achieved / target
  out <- problem$constraint_meta
  out$.achieved <- achieved
  out$.ratio <- ratio
  out$.residual <- ratio - 1
  out$.tolerance <- rep(tolerance, nrow(out))
  out$.se <- se
  out$.cv <- cv
  out$.moe <- moe
  out$.pass <- ratio <= 1 + tolerance
  out$.binding <- abs(ratio - 1) <= max(1e-6, 10 * tolerance)
  if (is.null(lambda)) {
    out$.multiplier <- rep(NA_real_, nrow(out))
    out$.sensitivity <- rep(NA_real_, nrow(out))
  } else {
    if (length(lambda) != nrow(out)) {
      stop("internal error: multiplier length does not match constraints",
           call. = FALSE)
    }
    out$.multiplier <- lambda
    sensitivity <- numeric(nrow(out))
    cv_idx <- metric == "cv"
    sensitivity[cv_idx] <- -2 * lambda[cv_idx] * target[cv_idx] *
      problem$total[cv_idx]^2
    moe_idx <- !cv_idx
    z <- stats::qnorm(1 - problem$alpha[moe_idx] / 2)
    sensitivity[moe_idx] <- -2 * lambda[moe_idx] *
      problem$domain_N[moe_idx]^2 * target[moe_idx] / z^2
    out$.sensitivity <- sensitivity
  }
  lower_violation <- allocation < problem$lower - tolerance
  upper_violation <- allocation > problem$upper + tolerance
  list(
    allocation = allocation,
    constraints = out,
    variance = variance,
    se = se,
    moe = moe,
    cv = cv,
    cost = sum(problem$cost * allocation),
    all_pass = all(out$.pass),
    lower_violation = lower_violation,
    upper_violation = upper_violation,
    bounds_pass = !any(lower_violation | upper_violation),
    feasibility_tolerance = tolerance
  )
}

#' Canonical acceptance tolerance implied by the public feasibility tolerance
#'
#' Integer repair works on the canonical scale `sum_h A_hk / n_h <= bound_k`,
#' while `.precision_from_allocation()` judges the public metric ratio
#' `sqrt(V_k / Vmax_k) <= 1 + tolerance`. Since `bound_k = Vmax_k - B_k` with
#' `B_k <= 0`, a slack of `tolerance` on the canonical scale can exceed the
#' public allowance. Convert it exactly instead, so an accepted integer
#' allocation always passes final validation.
#' @keywords internal
#' @noRd
.bethel_canonical_tolerance <- function(problem, tolerance) {
  Vmax <- problem$bound + problem$B
  pmax(0, Vmax * ((1 + tolerance)^2 - 1) / problem$bound)
}

#' Construct a deterministic feasible integer allocation
#' @keywords internal
#' @noRd
.integerize_bethel <- function(
  continuous,
  problem,
  control = .bethel_control()
) {
  tolerance <- control$tolerance
  n <- pmin(
    pmax(as.integer(ceiling(continuous - tolerance)), problem$lower_int),
    problem$upper_int
  )
  normalized_A <- sweep(problem$A, 2L, problem$bound, "/")
  accept <- 1 + .bethel_canonical_tolerance(problem, tolerance)
  iterations <- 0L
  repeat {
    value <- .bethel_constraint_value(normalized_A, n)
    violation <- pmax(value - 1, 0)
    if (all(value <= accept)) break
    can <- which(n < problem$upper_int)
    if (length(can) == 0L) {
      bad <- problem$constraint_ids[value > accept]
      stop(sprintf("no feasible integer allocation meets constraints: %s",
                   paste(bad, collapse = ", ")), call. = FALSE)
    }
    score <- vapply(can, function(h) {
      gain <- normalized_A[h, ] * (1 / n[h] - 1 / (n[h] + 1))
      sum(pmin(gain, violation)) / problem$cost[h]
    }, numeric(1L))
    h <- can[which.max(score)]
    if (!is.finite(score[which.max(score)]) || score[which.max(score)] <= 0) {
      stop("integer repair could not reduce the violated constraints",
           call. = FALSE)
    }
    n[h] <- n[h] + 1L
    iterations <- iterations + 1L
    if (iterations > control$integer_repair_iterations) {
      stop("integer allocation repair exceeded its iteration limit",
           call. = FALSE)
    }
  }

  repeat {
    can <- which(n > problem$lower_int)
    if (length(can) == 0L) break
    loss <- vapply(can, function(h) {
      sum(normalized_A[h, ] * (1 / (n[h] - 1) - 1 / n[h])) /
        problem$cost[h]
    }, numeric(1L))
    accepted <- FALSE
    for (h in can[order(loss, problem$stratum_ids[can])]) {
      candidate <- n
      candidate[h] <- candidate[h] - 1L
      if (all(.bethel_constraint_value(normalized_A, candidate) <= accept)) {
        n <- candidate
        accepted <- TRUE
        break
      }
    }
    if (!accepted) break
  }
  evaluation <- .precision_from_allocation(problem, n, tolerance = tolerance)
  if (!evaluation$all_pass || !evaluation$bounds_pass) {
    stop("internal error: integer allocation failed final validation",
         call. = FALSE)
  }
  list(
    allocation = as.integer(n),
    cost = evaluation$cost,
    constraints = evaluation$constraints,
    all_pass = evaluation$all_pass,
    repair_iterations = iterations
  )
}

#' Evaluate the priority-weighted objective for an allocation
#'
#' Returns the weighted relative variance and its per-component decomposition.
#' The value equals `Q0 + sum(q_h / n)` by construction.
#' @keywords internal
#' @noRd
.bethel_objective_value <- function(problem, allocation) {
  obj <- problem$objective
  if (is.null(obj)) return(NULL)
  reciprocal <- drop(crossprod(obj$A, 1 / allocation))
  variance <- obj$B + reciprocal
  cancellation_scale <- pmax(1, abs(obj$B), abs(reciprocal))
  tiny_negative <- variance < 0 &
    abs(variance) <= 100 * .Machine$double.eps * cancellation_scale
  variance[tiny_negative] <- 0
  if (any(variance < 0)) {
    stop("internal error: allocation produced a materially negative objective variance",
         call. = FALSE)
  }
  relvar <- variance / obj$total^2
  contribution <- obj$priority * relvar
  value <- sum(contribution)
  components <- data.frame(
    component = obj$component,
    name = obj$name,
    domain = obj$domain,
    level = obj$level,
    priority = obj$priority,
    .relvar = relvar,
    .cv = sqrt(relvar),
    .contribution = contribution,
    .share = if (value > 0) contribution / value else
      rep(NA_real_, length(relvar)),
    stringsAsFactors = FALSE
  )
  list(value = value, components = components)
}

#' Solve the fixed-budget objective problem by epsilon-constraint
#'
#' Appends the objective as one more reciprocal constraint column and root
#' searches the objective bound whose minimum cost equals the budget. Because
#' the model is convex with fixed coefficients this recovers the global
#' continuous optimum and keeps the existing KKT certification.
#' @keywords internal
#' @noRd
.bethel_budget_solve <- function(problem, budget, control = .bethel_control()) {
  A <- problem$A
  bound <- problem$bound
  cost <- problem$cost
  lower <- problem$lower
  upper <- problem$upper
  q_h <- problem$objective$q_h
  Q0 <- problem$objective$Q0
  K <- ncol(A)
  ids <- c(problem$constraint_ids, ".objective")
  Q_at <- function(n) Q0 + sum(q_h / n)
  solve_at <- function(qb) {
    .bethel_solve(
      cbind(A, q_h), c(bound, qb - Q0), cost, lower, upper,
      control = control, constraint_ids = ids
    )
  }

  if (K > 0L) {
    base <- .bethel_solve(
      A, bound, cost, lower, upper,
      control = control, constraint_ids = problem$constraint_ids
    )
    if (identical(base$classification, "infeasible")) {
      return(list(status = "infeasible_targets", infeasible = base$infeasible))
    }
    if (!identical(base$classification, "optimal")) {
      return(list(status = "failed", message = base$message))
    }
    base_allocation <- base$allocation
  } else {
    base_allocation <- lower
  }
  cost_min <- sum(cost * base_allocation)
  cost_max <- sum(cost * upper)
  budget_tol <- control$tolerance * max(1, abs(budget))

  if (budget < cost_min - budget_tol) {
    return(list(
      status = "unaffordable",
      cost_min = cost_min,
      shortfall = cost_min - budget
    ))
  }
  if (budget >= cost_max - budget_tol) {
    return(list(
      status = "optimal",
      allocation = upper,
      base_allocation = base_allocation,
      solved = NULL,
      lambda = rep(NA_real_, K),
      objective_multiplier = NA_real_,
      budget_sensitivity = NA_real_,
      objective_bound = Q_at(upper),
      budget_binding = FALSE,
      cost = cost_max,
      root_iterations = 0L
    ))
  }

  q_hi <- Q_at(base_allocation)
  q_lo <- Q_at(upper)
  span <- q_hi - q_lo
  if (!(span > 0)) {
    return(list(status = "failed",
                message = "objective bracket collapsed to a single point"))
  }
  bracket_lo <- q_lo + 1e-10 * span
  residual <- function(qb) {
    s <- solve_at(qb)
    if (!identical(s$classification, "optimal")) return(NA_real_)
    sum(cost * s$allocation) - budget
  }
  f_lo <- residual(bracket_lo)
  f_hi <- residual(q_hi)
  if (is.na(f_lo) || is.na(f_hi)) {
    return(list(status = "failed",
                message = "the objective bracket endpoints failed certification"))
  }
  if (f_lo <= 0) {
    return(list(
      status = "optimal",
      allocation = upper,
      base_allocation = base_allocation,
      solved = NULL,
      lambda = rep(NA_real_, K),
      objective_multiplier = NA_real_,
      budget_sensitivity = NA_real_,
      objective_bound = q_lo,
      budget_binding = FALSE,
      cost = cost_max,
      root_iterations = 0L
    ))
  }
  if (f_hi >= 0) {
    objective_bound <- q_hi
    root_iterations <- 0L
  } else {
    root <- stats::uniroot(
      residual, interval = c(bracket_lo, q_hi),
      f.lower = f_lo, f.upper = f_hi,
      tol = span * 1e-12, maxiter = 1000L
    )
    objective_bound <- root$root
    root_iterations <- root$iter
  }
  solved <- solve_at(objective_bound)
  if (!identical(solved$classification, "optimal")) {
    return(list(status = "failed", message = solved$message))
  }
  objective_multiplier <- solved$lambda[K + 1L]
  list(
    status = "optimal",
    allocation = solved$allocation,
    base_allocation = base_allocation,
    solved = solved,
    lambda = solved$lambda[seq_len(K)],
    objective_multiplier = objective_multiplier,
    # Local change in the optimal objective per unit of extra budget.
    budget_sensitivity = if (is.na(objective_multiplier) ||
                             objective_multiplier <= 0) NA_real_ else
      -1 / objective_multiplier,
    objective_bound = objective_bound,
    budget_binding = TRUE,
    cost = sum(cost * solved$allocation),
    root_iterations = root_iterations
  )
}

#' Construct a budget-feasible integer allocation
#'
#' `ceiling()` is not a safe start under a budget. This starts from the floor
#' of the continuous solution and repairs hard precision violations by greatest
#' constraint reduction per unit cost. When the budget-optimal point saturates
#' both the budget and a target, that repair cannot fit, so the fallback start
#' is the cheapest target-feasible integer allocation. From whichever start is
#' affordable, residual budget is spent by greatest objective reduction per unit
#' cost and improving pairwise exchanges are applied. The result is a feasible,
#' locally improved recommendation, not a globally optimal integer allocation.
#' @keywords internal
#' @noRd
.integerize_bethel_budget <- function(
  continuous,
  problem,
  budget,
  control = .bethel_control(),
  base_continuous = NULL,
  resolve = NULL
) {
  tolerance <- control$tolerance
  cost <- problem$cost
  lower_int <- problem$lower_int
  upper_int <- problem$upper_int
  q_h <- problem$objective$q_h
  normalized_A <- sweep(problem$A, 2L, problem$bound, "/")
  has_targets <- ncol(normalized_A) > 0L
  accept <- 1 + .bethel_canonical_tolerance(problem, tolerance)
  budget_tol <- tolerance * max(1, abs(budget))
  max_scan <- 200L

  cost_of <- function(x) sum(cost * x)
  feasible <- function(x) {
    !has_targets ||
      all(.bethel_constraint_value(normalized_A, x) <= accept)
  }
  gain_up <- function(h, x) q_h[h] * (1 / x[h] - 1 / (x[h] + 1))
  loss_down <- function(h, x) q_h[h] * (1 / (x[h] - 1) - 1 / x[h])

  if (cost_of(lower_int) > budget + budget_tol) {
    stop(sprintf(
      "no integer allocation fits the budget: the smallest allocation satisfying the stratum bounds costs %.6g against a budget of %.6g",
      cost_of(lower_int), budget
    ), call. = FALSE)
  }

  # Repair upward by greatest constraint reduction per unit cost. Returns NULL
  # when the stratum upper bounds cannot absorb the remaining violation.
  repairs <- 0L
  repair_targets <- function(x) {
    repeat {
      value <- .bethel_constraint_value(normalized_A, x)
      if (all(value <= accept)) return(x)
      violation <- pmax(value - 1, 0)
      can <- which(x < upper_int)
      if (length(can) == 0L) return(NULL)
      score <- vapply(can, function(h) {
        sum(pmin(normalized_A[h, ] * (1 / x[h] - 1 / (x[h] + 1)), violation)) /
          cost[h]
      }, numeric(1L))
      best <- which.max(score)
      if (!is.finite(score[best]) || score[best] <= 0) return(NULL)
      x[can[best]] <- x[can[best]] + 1L
      repairs <<- repairs + 1L
      if (repairs > control$integer_repair_iterations) {
        stop("integer allocation repair exceeded its iteration limit",
             call. = FALSE)
      }
    }
  }

  # Shed the cheapest objective loss per unit of cost saved, never breaking a
  # target. Returns NULL when nothing further can be removed.
  trim_to_budget <- function(x) {
    if (is.null(x)) return(NULL)
    guard <- 0L
    while (cost_of(x) > budget + budget_tol) {
      can <- which(x > lower_int)
      if (length(can) == 0L) return(NULL)
      penalty <- vapply(can, function(h) loss_down(h, x) / cost[h], numeric(1L))
      accepted <- FALSE
      for (h in can[order(penalty, problem$stratum_ids[can])]) {
        candidate <- x
        candidate[h] <- candidate[h] - 1L
        if (feasible(candidate)) {
          x <- candidate
          accepted <- TRUE
          break
        }
      }
      if (!accepted) return(NULL)
      guard <- guard + 1L
      if (guard > control$integer_repair_iterations) return(NULL)
    }
    x
  }

  clip <- function(x) pmin(pmax(as.integer(x), lower_int), upper_int)
  affordable <- function(x) {
    if (is.null(x) || cost_of(x) > budget + budget_tol) NULL else x
  }
  # No single rounding of the continuous optimum is reliable. Rounding down can
  # leave too little slack to repair a binding target; rounding up can be
  # untrimmable back inside the budget. When both fail, re-solve at a budget
  # reduced by the rounding excess so that the rounding does fit, which stays
  # far closer to the optimum than falling back to the cheapest design.
  starts <- list()
  candidate_budget <- budget
  point <- continuous
  for (attempt in seq_len(4L)) {
    up <- repair_targets(clip(ceiling(point - tolerance)))
    down <- affordable(repair_targets(clip(floor(point + tolerance))))
    trimmed <- trim_to_budget(up)
    if (!is.null(down)) starts[[paste0("floor", attempt)]] <- down
    if (!is.null(trimmed)) starts[[paste0("ceiling", attempt)]] <- trimmed
    if (length(starts) > 0L || is.null(resolve) || is.null(up)) break
    excess <- cost_of(up) - budget
    if (!(excess > 0)) break
    candidate_budget <- candidate_budget - excess
    point <- resolve(candidate_budget)
    if (is.null(point)) break
  }
  starts$cheapest <- affordable(if (!has_targets) lower_int else
    .integerize_bethel(base_continuous %||% continuous, problem, control)$allocation)
  starts <- starts[!vapply(starts, is.null, logical(1L))]
  if (length(starts) == 0L) {
    cheapest_cost <- if (has_targets) {
      cost_of(.integerize_bethel(
        base_continuous %||% continuous, problem, control
      )$allocation)
    } else cost_of(lower_int)
    # Round the suggestion up so the printed number is itself affordable,
    # rather than a truncation of a cost it does not actually cover.
    suggest <- ceiling(cheapest_cost * 100) / 100
    stop(sprintf(
      paste0(
        "no integer allocation meets the precision targets within the ",
        "budget: the cheapest target-feasible integer allocation costs ",
        "%.6g against a budget of %.6g.\n",
        "  Raise 'budget' to at least %.2f, or relax the targets.\n",
        "  A budget carried over from a minimum-cost fit falls just short ",
        "like this by construction: that cost is the continuous optimum, ",
        "and whole units cost slightly more."
      ),
      cheapest_cost, budget, suggest
    ), call. = FALSE)
  }

  exchange_once <- function(x) {
    remaining <- budget - cost_of(x)
    donors <- which(x > lower_int)
    receivers <- which(x < upper_int)
    if (length(donors) == 0L || length(receivers) == 0L) return(NULL)
    give <- vapply(donors, function(h) loss_down(h, x), numeric(1L))
    take <- vapply(receivers, function(h) gain_up(h, x), numeric(1L))
    gain <- outer(take, give, "-")
    affordable <- outer(cost[receivers], cost[donors], "-") <=
      remaining + budget_tol
    gain[!(affordable & outer(receivers, donors, "!="))] <- -Inf
    threshold <- 1e-10 * max(sum(q_h / x), .Machine$double.eps)
    if (max(gain) <= threshold) return(NULL)
    ranked <- utils::head(order(gain, decreasing = TRUE), max_scan)
    for (idx in ranked) {
      if (gain[idx] <= threshold) break
      i <- (idx - 1L) %% length(receivers) + 1L
      j <- (idx - 1L) %/% length(receivers) + 1L
      candidate <- x
      candidate[donors[j]] <- candidate[donors[j]] - 1L
      candidate[receivers[i]] <- candidate[receivers[i]] + 1L
      if (feasible(candidate)) return(candidate)
    }
    NULL
  }

  # Spend affordable residual budget by greatest objective reduction per unit
  # cost (adding units can never break a target), then apply improving
  # exchanges, until neither moves.
  improve <- function(x) {
    spends <- 0L
    exchanges <- 0L
    repeat {
      changed <- FALSE
      repeat {
        remaining <- budget - cost_of(x)
        can <- which(x < upper_int & cost <= remaining)
        if (length(can) == 0L) break
        score <- vapply(can, function(h) gain_up(h, x) / cost[h], numeric(1L))
        best <- which.max(score)
        if (!is.finite(score[best]) || score[best] <= 0) break
        x[can[best]] <- x[can[best]] + 1L
        spends <- spends + 1L
        changed <- TRUE
        if (spends > control$integer_repair_iterations) break
      }
      swapped <- exchange_once(x)
      if (!is.null(swapped)) {
        x <- swapped
        exchanges <- exchanges + 1L
        changed <- TRUE
      }
      if (!changed || exchanges > control$max_iterations) break
    }
    list(allocation = x, spends = spends, exchanges = exchanges,
         objective = sum(q_h / x))
  }

  improved <- lapply(starts, improve)
  chosen <- improved[[which.min(
    vapply(improved, function(z) z$objective, numeric(1L))
  )]]
  n <- chosen$allocation
  spends <- chosen$spends
  exchanges <- chosen$exchanges
  start_used <- names(starts)[which.min(
    vapply(improved, function(z) z$objective, numeric(1L))
  )]

  evaluation <- .precision_from_allocation(problem, n, tolerance = tolerance)
  if (!evaluation$all_pass || !evaluation$bounds_pass) {
    stop("internal error: integer allocation failed final validation",
         call. = FALSE)
  }
  if (evaluation$cost > budget + budget_tol) {
    stop("internal error: integer allocation exceeded the budget",
         call. = FALSE)
  }
  objective <- .bethel_objective_value(problem, n)
  list(
    allocation = as.integer(n),
    cost = evaluation$cost,
    constraints = evaluation$constraints,
    all_pass = evaluation$all_pass,
    objective = objective$components,
    objective_value = objective$value,
    budget_residual = budget - evaluation$cost,
    repair_iterations = repairs,
    start = start_used,
    spend_iterations = spends,
    exchange_iterations = exchanges
  )
}

#' Public-wrapper engine for joint constrained allocation
#' @keywords internal
#' @noRd
.n_alloc_bethel <- function(
  frame,
  measures,
  targets,
  unit_cost,
  alpha,
  deff,
  resp_rate,
  min_n_stratum,
  objective = NULL,
  budget = NULL
) {
  control <- .bethel_control()
  problem <- .build_bethel_problem(
    frame = frame,
    measures = measures,
    targets = targets,
    unit_cost = unit_cost,
    deff = deff,
    resp_rate = resp_rate,
    alpha = alpha,
    min_n_stratum = min_n_stratum,
    objective = objective
  )
  mode <- if (is.null(objective)) "targets" else "budget_objective"
  if (mode == "targets") {
    solved <- .bethel_solve(
      problem$A, problem$bound, problem$cost,
      problem$lower, problem$upper,
      control = control,
      constraint_ids = problem$constraint_ids
    )
    if (solved$classification == "infeasible") {
      stop(sprintf("precision targets are infeasible: %s",
                   paste(solved$infeasible, collapse = ", ")), call. = FALSE)
    }
    if (solved$classification != "optimal") {
      stop(sprintf("joint allocation solver failed certification: %s",
                   solved$message), call. = FALSE)
    }
    allocation <- solved$allocation
    lambda <- solved$lambda
    optimization <- solved
    optimization$allocation <- NULL
    optimization$lambda_scaled <- NULL
  } else {
    fit <- .bethel_budget_solve(problem, budget, control)
    if (identical(fit$status, "infeasible_targets")) {
      stop(sprintf("precision targets are infeasible: %s",
                   paste(fit$infeasible, collapse = ", ")), call. = FALSE)
    }
    if (identical(fit$status, "unaffordable")) {
      stop(sprintf(
        "'budget' cannot fund the precision targets: the cheapest target-feasible design costs %.6g, short by %.6g",
        fit$cost_min, fit$shortfall
      ), call. = FALSE)
    }
    if (!identical(fit$status, "optimal")) {
      stop(sprintf("joint allocation solver failed certification: %s",
                   fit$message %||% "unknown"), call. = FALSE)
    }
    allocation <- fit$allocation
    lambda <- fit$lambda
    solved <- fit$solved
    optimization <- if (is.null(solved)) {
      list(classification = "optimal", converged = TRUE,
           message = "the budget exceeds the cost of the upper allocation bounds",
           feasibility_tolerance = control$tolerance,
           kkt_tolerance = control$kkt_tolerance,
           active_lower = integer(0),
           active_upper = seq_along(allocation))
    } else {
      optimization <- solved
      optimization$allocation <- NULL
      optimization$lambda_scaled <- NULL
      optimization
    }
    optimization$budget <- budget
    optimization$budget_binding <- fit$budget_binding
    optimization$budget_residual <- budget - sum(problem$cost * allocation)
    optimization$budget_sensitivity <- fit$budget_sensitivity
    optimization$objective_bound <- fit$objective_bound
    optimization$objective_multiplier <- fit$objective_multiplier
    optimization$root_iterations <- fit$root_iterations
  }
  continuous <- .precision_from_allocation(
    problem, allocation, lambda = lambda,
    tolerance = control$tolerance
  )
  continuous_objective <- .bethel_objective_value(problem, allocation)
  integer <- if (mode == "targets") {
    .integerize_bethel(allocation, problem, control)
  } else {
    .integerize_bethel_budget(
      allocation, problem, budget, control,
      base_continuous = fit$base_allocation,
      resolve = function(b) {
        retry <- .bethel_budget_solve(problem, b, control)
        if (identical(retry$status, "optimal")) retry$allocation else NULL
      }
    )
  }
  decision <- allocation
  decision_int <- integer$allocation
  n <- .bethel_to_ultimate(problem, decision)
  n_int <- .bethel_to_ultimate(problem, decision_int)
  detail <- data.frame(
    stratum = problem$stratum_ids,
    N = problem$population_N,
    stringsAsFactors = FALSE
  )
  if (problem$stages == 1L) {
    detail$unit_cost <- problem$cost
  } else {
    detail$N_psu <- problem$stage$N_psu
    if (problem$stages == 3L) detail$N_ssu <- problem$stage$N_ssu
    detail$cost_psu <- problem$stage$cost_psu
    detail$cost_ssu <- problem$stage$cost_ssu
    if (problem$stages == 3L) detail$cost_tsu <- problem$stage$cost_tsu
    detail$n_per_psu <- problem$stage$n_per_psu
    if (problem$stages == 3L) detail$n_per_ssu <- problem$stage$n_per_ssu
    detail$n_psu <- decision
    # Count columns use R's standard double representation consistently.
    # The "_int" suffix records whole-number semantics, not storage type.
    detail$n_psu_int <- as.numeric(decision_int)
  }
  detail$n <- n
  detail$n_int <- n_int
  detail$weight <- problem$population_N / n
  detail$.lower <- problem$lower * problem$stage$take
  detail$.upper <- problem$upper * problem$stage$take
  if (problem$stages > 1L) {
    detail$.lower_psu <- problem$lower
    detail$.upper_psu <- problem$upper
  }
  detail$.binding <- abs(decision - problem$lower) <= 1e-6 |
    abs(decision - problem$upper) <= 1e-6
  if (any(problem$take_all)) detail$take_all <- problem$take_all
  operational <- list(
    n = sum(n_int),
    cost = integer$cost,
    constraints = integer$constraints,
    all_pass = integer$all_pass,
    repair_iterations = integer$repair_iterations
  )
  if (mode == "budget_objective") {
    operational$objective <- integer$objective
    operational$objective_value <- integer$objective_value
    operational$budget_residual <- integer$budget_residual
    operational$start <- integer$start
    operational$spend_iterations <- integer$spend_iterations
    operational$exchange_iterations <- integer$exchange_iterations
  }
  params <- list(
    frame = frame,
    measures = problem$measures,
    targets = problem$targets,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    unit_cost = unit_cost,
    cost_h = problem$cost,
    min_n_stratum = min_n_stratum,
    mode = mode,
    alloc = "bethel",
    n_h = n,
    n_psu_h = if (problem$stages > 1L) decision else NULL,
    stages = problem$stages,
    problem = problem,
    achieved = list(n = sum(n), cost = continuous$cost)
  )
  if (mode == "budget_objective") {
    params$objective <- problem$objective$spec
    params$budget <- budget
    params$achieved$objective_value <- continuous_objective$value
  }
  params$feasibility_tolerance <- control$tolerance
  .new_svyplan_n(
    n = sum(n),
    type = "alloc",
    method = "bethel",
    params = params,
    se = NA_real_,
    moe = NA_real_,
    cv = NA_real_,
    targets = problem$targets,
    detail = detail,
    binding = problem$constraint_ids[continuous$constraints$.binding],
    operational = operational,
    constraints = continuous$constraints,
    optimization = optimization,
    objective = continuous_objective$components,
    objective_value = continuous_objective$value
  )
}

#' Precision-wrapper engine for joint constrained allocation
#' @keywords internal
#' @noRd
.prec_alloc_bethel <- function(
  frame,
  n,
  measures,
  targets,
  unit_cost,
  alpha,
  deff,
  resp_rate,
  min_n_stratum = NULL,
  objective = NULL,
  budget = NULL,
  .allow_fractional_stages = FALSE
) {
  if (is.null(n)) stop("'n' is required", call. = FALSE)
  problem <- .build_bethel_problem(
    frame = frame,
    measures = measures,
    targets = targets,
    unit_cost = unit_cost,
    deff = deff,
    resp_rate = resp_rate,
    alpha = alpha,
    min_n_stratum = min_n_stratum,
    objective = objective
  )
  decision <- .bethel_to_decision(
    problem,
    n,
    require_whole_stages = !.allow_fractional_stages
  )
  evaluated <- .precision_from_allocation(problem, decision)
  ultimate <- .bethel_to_ultimate(problem, evaluated$allocation)
  bounds <- data.frame(
    stratum = problem$stratum_ids,
    n = ultimate,
    .lower = problem$lower * problem$stage$take,
    .upper = problem$upper * problem$stage$take,
    .lower_violation = evaluated$lower_violation,
    .upper_violation = evaluated$upper_violation,
    .pass = !(evaluated$lower_violation | evaluated$upper_violation),
    stringsAsFactors = FALSE
  )
  assessed <- .bethel_objective_value(problem, evaluated$allocation)
  params <- list(
    frame = frame,
    measures = problem$measures,
    targets = problem$targets,
    n = ultimate,
    alpha = alpha,
    deff = deff,
    resp_rate = resp_rate,
    unit_cost = unit_cost,
    cost_h = problem$cost,
    min_n_stratum = min_n_stratum,
    stages = problem$stages,
    problem = problem,
    achieved = list(n = sum(ultimate), cost = evaluated$cost),
    feasibility_tolerance = .bethel_control()$tolerance
  )
  if (!is.null(assessed)) {
    params$objective <- problem$objective$spec
    params$budget <- budget
    params$achieved$objective_value <- assessed$value
    params$budget_residual <- if (is.null(budget)) NA_real_ else
      budget - evaluated$cost
  }
  .new_svyplan_prec(
    se = evaluated$se,
    moe = evaluated$moe,
    cv = evaluated$cv,
    type = "alloc",
    method = "bethel",
    params = params,
    detail = evaluated$constraints,
    bounds = bounds,
    objective = assessed$components,
    objective_value = assessed$value
  )
}
