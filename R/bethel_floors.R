#' Regional sample floors of a joint element allocation
#'
#' `min_n_domain` is a table with `domain`, `level` and `min_n`, all on one
#' partition column of `frame`. Returns `NULL` without floors, otherwise
#' the table, the stratum membership of each floored level and its indices.
#' @keywords internal
#' @noRd
.bethel_domain_floors <- function(frame, min_n_domain) {
  if (is.null(min_n_domain)) return(NULL)
  columns <- c("domain", "level", "min_n")
  if (!is.data.frame(min_n_domain) ||
        !all(columns %in% names(min_n_domain))) {
    stop("joint 'min_n_domain' must be a table with columns 'domain', 'level' and 'min_n'",
         call. = FALSE)
  }
  spec <- data.frame(domain = as.character(min_n_domain$domain),
                     level = as.character(min_n_domain$level),
                     min_n = min_n_domain$min_n, stringsAsFactors = FALSE)
  if (nrow(spec) == 0L) return(NULL)
  if (anyNA(spec$domain) || any(!nzchar(spec$domain)) ||
        length(unique(spec$domain)) != 1L || any(spec$domain == ".overall")) {
    stop("joint 'min_n_domain' must name one domain column of 'frame'",
         call. = FALSE)
  }
  if (anyNA(spec$level) || any(!nzchar(spec$level)) ||
        anyDuplicated(spec$level)) {
    stop("joint 'min_n_domain' must give each level once",
         call. = FALSE)
  }
  if (!is.numeric(spec$min_n) ||
        any(!is.finite(spec$min_n) | spec$min_n < 0)) {
    stop("'min_n_domain$min_n' must hold nonnegative finite counts",
         call. = FALSE)
  }
  domain <- .bethel_domain_indices(frame, spec$domain, spec$level,
                                   "sample floor")
  list(spec = spec, membership = domain$membership * 1,
       indices = domain$indices)
}

#' Required and achieved counts of each floored level
#'
#' A whole design meets `ceiling(min_n)`. A continuous design meets
#' `min_n` up to rounding of the count.
#' @keywords internal
#' @noRd
.bethel_floor_table <- function(floors, x, whole = FALSE) {
  if (is.null(floors)) return(NULL)
  out <- floors$spec
  out$required <- if (whole) ceiling(out$min_n) else out$min_n
  out$achieved <- drop(crossprod(floors$membership, as.numeric(x)))
  slack <- if (whole) 0 else 16 * .Machine$double.eps * pmax(1, out$required)
  out$.pass <- out$achieved >= out$required - slack
  out
}

#' Refuse floors the strata cannot hold, continuously or in whole units
#' @keywords internal
#' @noRd
.bethel_floor_capacity <- function(floors, upper, upper_int) {
  if (is.null(floors)) return(invisible(NULL))
  for (whole in c(FALSE, TRUE)) {
    held <- .bethel_floor_table(floors, if (whole) upper_int else upper,
                                whole)
    short <- !held$.pass
    if (any(short)) {
      stop(sprintf(
        "joint 'min_n_domain' exceeds the %scapacity of level%s %s",
        if (whole) "whole " else "", if (sum(short) > 1L) "s" else "",
        paste(sQuote(held$level[short]), collapse = ", ")
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Each stratum's least `cost x + weight / x` on its box
#'
#' A nonpositive cost with a positive weight goes to the upper bound, and a
#' zero weight with a positive cost to the lower one.
#' @keywords internal
#' @noRd
.bethel_box_min <- function(cost, weight, lower, upper) {
  weight <- pmax(weight, 0)
  x <- ifelse(cost > 0, sqrt(weight / ifelse(cost > 0, cost, 1)),
              ifelse(weight > 0 | cost < 0, upper, lower))
  pmin(pmax(x, lower), upper)
}

#' The least `sum(cost x + weight / x)` on the box with every floor met
#'
#' The levels are disjoint, so each is solved alone. A binding floor's
#' multiplier `mu` lowers every linear cost of its level by the same amount,
#' and the level's count rises with `mu`, so `mu` is found by bisection.
#' A stratum whose cost the multiplier reaches jumps from its lower toward
#' its upper bound there. Every split of the floor among such strata costs
#' the same, so a floor that stops at one fills them with the least
#' `sum(weight / x)`, or with no weight the least `sum(tie / x)`, the limit
#' of the weighted split as the weights go to zero when `tie` holds the
#' target loads. Returns the sizes `x` and each floor's `mu`, in the
#' units of `cost`.
#' @keywords internal
#' @noRd
.bethel_floor_min <- function(cost, weight, lower, upper, floors,
                              tie = rep(1, length(cost))) {
  x <- .bethel_box_min(cost, weight, lower, upper)
  if (is.null(floors)) return(list(x = x, mu = numeric(0)))
  mu <- numeric(nrow(floors$spec))
  for (k in seq_along(mu)) {
    j <- floors$indices[[k]]
    need <- floors$spec$min_n[k]
    if (sum(x[j]) >= need) next
    cj <- cost[j]
    wj <- weight[j]
    at <- function(m) .bethel_box_min(cj - m, wj, lower[j], upper[j])
    below <- 0
    above <- max(cj)
    for (i in seq_len(200L)) {
      mid <- below + (above - below) / 2
      if (mid <= below || mid >= above) break
      if (sum(at(mid)) >= need) above <- mid else below <- mid
    }
    z <- at(above)
    flat <- which(cj >= below & cj <= above)
    if (length(flat) > 0L) {
      ## A stratum whose cost sits at the break jumps in size there, since
      ## `cost - mu` vanishes, and every split of the rest of the floor among
      ## such strata costs `mu` a unit. The split that least raises their
      ## reciprocal terms follows `sqrt(weight)`, or the target loads when
      ## they carry no weight.
      z[flat] <- lower[j][flat]
      share <- if (any(wj[flat] > 0)) wj[flat] else tie[j][flat]
      z[flat] <- .bethel_fill(need - sum(z[-flat]), lower[j][flat],
                              upper[j][flat], share)
      above <- min(cj[flat])
    } else {
      ## Near `mu = cost` a small weight makes `cost - mu` cancel, so the
      ## interior strata are scaled to meet the floor exactly.
      inside <- wj > 0 & z > lower[j] & z < upper[j]
      want <- need - sum(z[!inside])
      if (any(inside) && want > 0) {
        z[inside] <- pmin(pmax(z[inside] * want / sum(z[inside]),
                               lower[j][inside]), upper[j][inside])
      }
    }
    mu[k] <- above
    x[j] <- z
  }
  list(x = x, mu = mu)
}

#' Sizes with a given sum on their box, least `sum(w / x)`
#'
#' The sizes are `sqrt(w)` scaled and clamped to the box, the scale found
#' by bisection. Sizes with no weight stay at their lower bound until the
#' others reach their upper ones, then fill in order.
#' @keywords internal
#' @noRd
.bethel_fill <- function(want, lower, upper, w) {
  if (want <= sum(lower)) return(lower)
  if (want >= sum(upper)) return(upper)
  r <- sqrt(pmax(w, 0))
  at <- function(t) pmin(pmax(t * r, lower), upper)
  if (any(r > 0) && sum(ifelse(r > 0, upper, lower)) >= want) {
    hi <- 1
    while (sum(at(hi)) < want) hi <- hi * 2
    lo <- 0
    for (i in seq_len(200L)) {
      mid <- lo + (hi - lo) / 2
      if (mid <= lo || mid >= hi) break
      if (sum(at(mid)) >= want) hi <- mid else lo <- mid
    }
    return(at(hi))
  }
  x <- ifelse(r > 0, upper, lower)
  room <- ifelse(r > 0, 0, upper - lower)
  short <- want - sum(x)
  x + pmin(room, pmax(0, short - c(0, utils::head(cumsum(room), -1L))))
}

#' Each floor as a whole-search load on the shortfall from capacity
#'
#' For level `k` with capacity `C` and whole floor `F`, the load is
#' `(C - sum x + 1) / (C - F + 1)`, at most 1 exactly when `sum x >= F`.
#' Each stratum contributes `(U_h - x_h)` scaled, plus 1 on the level's
#' first stratum, so every state's load is nonnegative and the column falls
#' as `x` grows. The scale leaves room for the search's `1 + 1e-9`
#' acceptance below one unit of shortfall.
#' @keywords internal
#' @noRd
.bethel_floor_loads <- function(floors, upper_int) {
  if (is.null(floors)) return(NULL)
  M <- floors$membership
  capacity <- drop(crossprod(M, upper_int))
  denom <- capacity - ceiling(floors$spec$min_n) + 1
  scale <- (1 + 1e-9) / (denom + pmin(1e-9 * denom, 0.25))
  g0 <- M * upper_int
  first <- cbind(vapply(floors$indices, function(i) i[1L], 1L),
                 seq_len(ncol(M)))
  g0[first] <- g0[first] + 1
  list(g0 = sweep(g0, 2L, scale, "*"), g1 = -sweep(M, 2L, scale, "*"))
}

#' The cheapest continuous design meeting every floor
#' @keywords internal
#' @noRd
.bethel_floor_cheapest <- function(problem) {
  .bethel_floor_min(problem$cost, numeric(length(problem$cost)),
                    problem$lower, problem$upper, problem$floors)$x
}

#' A fixed-take whole problem with every floor as one more column
#'
#' The floor columns follow the caller's columns.
#' @keywords internal
#' @noRd
.bethel_whole_problem <- function(problem, o1, o2, g1, g2, g0 = 0) {
  H <- length(problem$cost)
  g1 <- matrix(as.numeric(g1), H)
  g2 <- matrix(as.numeric(g2), H, ncol(g1))
  g0 <- matrix(as.numeric(g0), H, ncol(g1))
  fl <- .bethel_floor_loads(problem$floors, problem$upper_int)
  if (!is.null(fl)) {
    g0 <- cbind(g0, fl$g0)
    g1 <- cbind(g1, fl$g1)
    g2 <- cbind(g2, matrix(0, H, ncol(fl$g1)))
  }
  .whole_problem(o1, o2, g1, g2, problem$lower_int, problem$upper_int,
                 g0 = g0)
}

#' Whether a whole design meets every floor
#' @keywords internal
#' @noRd
.bethel_floors_met <- function(problem, x) {
  is.null(problem$floors) ||
    all(.bethel_floor_table(problem$floors, x, whole = TRUE)$.pass)
}

#' Starting multipliers of a whole search, from a continuous solve
#'
#' `cont` is a cost solve, its multipliers in cost per unit of each row.
#' `lead` holds multipliers of leading columns, such as a budget's, already
#' in load units. A row's load is its variance over `bound * accept`, times
#' the target scale, and a floor's load falls by its scale per unit, so the
#' multipliers are converted by those factors. Missing ones are 0, and
#' `NULL` comes back when the solve has none.
#' @keywords internal
#' @noRd
.bethel_whole_lambda <- function(problem, q, cont, accept,
                                 lead = numeric(0), per = 1) {
  if (is.null(cont) || is.null(cont$lambda)) return(NULL)
  K <- ncol(problem$A)
  lam <- cont$lambda[seq_len(K)] / per
  lam[!is.finite(lam)] <- 0
  targets <- lam * problem$bound * accept / .whole_target_scale()
  floors <- numeric(0)
  if (!is.null(problem$floors)) {
    D <- nrow(problem$floors$spec)
    mu <- (cont$floor_multiplier %||% rep(0, D)) / per
    mu[!is.finite(mu)] <- 0
    scale <- -colSums(q$g1[, ncol(q$g1) - rev(seq_len(D)) + 1L,
                          drop = FALSE]) /
      colSums(problem$floors$membership)
    floors <- mu / scale
  }
  out <- c(lead, targets, floors)
  if (length(out) != q$K || any(!is.finite(out)) || any(out < 0)) return(NULL)
  out
}
