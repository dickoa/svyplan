#' Filter a plan for boundary determination
#'
#' `strata_bound()` orders `unit_cost` by ascending `x`, from the lowest
#' stratum to the highest, while `n_alloc()` orders it by frame row. A vector
#' carried in a profile is written for the latter, so applying it here would
#' attach costs to the wrong strata whenever the two lengths happen to match.
#' A scalar carries no ordering and is a genuine design constant, so it
#' passes through.
#' @keywords internal
#' @noRd
.strata_plan_defaults <- function(plan) {
  if (is.null(plan)) return(NULL)
  if (!inherits(plan, "svyplan")) {
    stop("'plan' must be a svyplan object", call. = FALSE)
  }
  cost <- plan$defaults$unit_cost
  if (!is.null(cost) && length(cost) > 1L) {
    stop(
      "'plan' carries a vector 'unit_cost'; strata_bound() orders costs from the lowest to the highest stratum, not by allocation-frame row, so pass 'unit_cost' explicitly here",
      call. = FALSE
    )
  }
  plan
}

#' Strata boundaries for survey design
#'
#' Determine where to cut a continuous stratification variable to form
#' useful strata. Supports four methods: cumulative root frequency
#' (Dalenius-Hodges), geometric progression, coordinate optimization inspired
#' by Lavall\enc{é}{e}e-Hidiroglou, and random-restart local search inspired
#' by Kozak.
#'
#' @param x Numeric vector: finite stratification variable values. Must not
#'   contain missing values.
#' @param n_strata Required integer: number of strata (including take-all if
#'   `take_all_above` is specified). Must be >= 2. Every sampled stratum must
#'   hold at least two population units, so a request that only sparse or tied
#'   data could satisfy is refused rather than met with a stratum that carries
#'   no within-stratum variance. A take-all stratum is enumerated instead of
#'   sampled, so one unit is enough there.
#' @param ... Unused. Present so that every optional argument must be named.
#'   Unused arguments are rejected.
#' @param n Target total sample size, supplied as a whole number. Specify at
#'   most one of `n` or `cv`. Required for methods `"lh"` and `"kozak"`.
#' @param cv Target coefficient of variation (relative standard error).
#'   For example, `cv = 0.05` means the standard error of the estimated
#'   population total or mean should be at most 5 percent of the estimate.
#'   Specify at most one of `n` or `cv`. Required for methods `"lh"` and
#'   `"kozak"`.
#' @param method Stratification method: `"cumrootf"` (Dalenius-Hodges),
#'   `"geo"` (geometric), `"lh"` (LH-inspired coordinate optimization), or
#'   `"kozak"` (Kozak-inspired random-restart local search). The short method
#'   names are retained for compatibility and do not claim exact
#'   implementations of the published algorithms. Default `"lh"`.
#' @param alloc Allocation rule: `"proportional"`, `"neyman"`,
#'   `"optimal"`, or `"power"` (Bankier compromise). Default `"neyman"`.
#'   See Details.
#' @param alloc_q Bankier power parameter, used only when `alloc = "power"`.
#'   Numeric scalar in \eqn{[0, 1]}. At `alloc_q = 1` the allocation equals
#'   Neyman. At `alloc_q = 0` it yields near-equal subnational CVs.
#'   Default 0.5.
#' @param unit_cost Per-stratum unit costs, ordered from lowest to highest
#'   stratum. Scalar (equal costs) or vector of length `n_strata`.
#'   Default `NULL` (equal unit costs, \eqn{c_h = 1} for all strata),
#'   in which case `"optimal"` and `"neyman"` coincide.
#' @param take_all_above Take-all threshold. Units with `x >= take_all_above` form a
#'   census stratum.
#' @param deff Design effect multiplier (> 0), default 1. Inflates the
#'   variance of the stratified mean exactly as it does in [n_alloc()], so
#'   that `cv` and `n` mean the same thing in both. A scalar value leaves the
#'   boundaries unchanged. See Details.
#' @param resp_rate Expected response rate, in (0, 1\], default 1. With a
#'   response rate below 1, `n` is the sample fielded and `n * resp_rate` the
#'   number expected to respond, matching [n_alloc()].
#' @param n_class Positive whole number of histogram bins. Applies to
#'   `method = "cumrootf"` only. Supplying it with another method is an
#'   error rather than a silent no-op.
#'   Default `NULL` (Freedman-Diaconis rule).
#' @param max_iter Positive whole number of maximum iterations. Applies to
#'   `method = "lh"` and `"kozak"` only. Supplying it with another method is
#'   an error rather than a silent no-op. Default `NULL` (= 200).
#' @param n_restart Positive whole number of random restarts. Applies to
#'   `method = "kozak"` only. Supplying it with another method is an error
#'   rather than a silent no-op. Default `NULL` (= 10 * `n_strata`).
#' @param plan Optional [svyplan()] object providing design defaults. It can
#'   supply `alloc`, `alloc_q`, `deff`, `resp_rate`, and a scalar
#'   `unit_cost`. A profile carrying a *vector* `unit_cost` is rejected,
#'   because this function orders costs from the lowest to the highest
#'   stratum while [n_alloc()] orders them by frame row, so a vector written
#'   for one would be silently misapplied by the other.
#'
#' @return A `svyplan_strata` object with components:
#' \describe{
#'   \item{boundaries}{Numeric vector of cutpoints (length `n_strata - 1`).}
#'   \item{n_strata}{Number of strata.}
#'   \item{n}{Total sample size.}
#'   \item{cv}{Coefficient of variation achieved by the integer
#'     allocation in `$strata$n` (the continuous optimum's cv is kept in
#'     `params$cv_continuous`, and the requested target, if any, in
#'     `params$cv_target`).}
#'   \item{strata}{Data frame with per-stratum summaries: `stratum`,
#'     `lower`, `upper`, `N`, `share`, `sd`, `mean`, `n`, and `take_all`.
#'     The `N`, `sd`, and `mean` columns are exactly what [n_alloc()]
#'     expects, so the table can be passed straight to it.}
#'   \item{method}{Algorithm used.}
#'   \item{alloc}{Allocation method name (character).}
#'   \item{params}{List of additional parameters.}
#'   \item{converged}{Logical (for iterative methods).}
#' }
#'
#' @details
#' ## Choosing a method
#'
#' The four methods differ in approach and use case:
#'
#' - **cumrootf**: Dalenius-Hodges (1959) cumulative root frequency rule.
#'   Non-iterative, does not require `n` or `cv`. Best for a quick
#'   first look at reasonable boundary positions.
#' - **geo**: Gunning-Horgan (2004) geometric progression. Non-iterative,
#'   requires `x > 0`. Works well for right-skewed positive variables
#'   (e.g. income, revenue) where a log-scale spacing is natural.
#' - **lh** (default): coordinate-wise optimization inspired by
#'   Lavall\enc{é}{e}e-Hidiroglou (1988). It starts from quantile boundaries
#'   and repeatedly performs bounded one-dimensional minimizations. Requires
#'   `n` or `cv`. This is a local heuristic, not the published LH recurrence
#'   and not a globally optimal algorithm.
#' - **kozak**: random-restart adjacent-boundary local search inspired by
#'   Kozak (2004). Requires `n` or `cv`. It can explore more starting points
#'   than `"lh"`, but it can still finish at a local optimum and provides no
#'   global-optimality guarantee.
#'
#' In summary, `"lh"` is the faster iterative heuristic. `"kozak"` spends more
#' computation on random restarts and may find a better local solution. Use
#' `"cumrootf"` or `"geo"` when you do not yet have `n` or `cv`.
#'
#' Allocation is controlled by the `alloc` parameter. Four methods are
#' available:
#' - **proportional**: \eqn{n_h \propto N_h}{n_h ~ N_h}.
#' - **neyman**: \eqn{n_h \propto N_h S_h}{n_h ~ N_h * S_h}.
#'   Minimizes the national CV when unit costs are equal.
#' - **optimal**: \eqn{n_h \propto N_h S_h / \sqrt{c_h}}{n_h ~ N_h * S_h / sqrt(c_h)}.
#'   Accounts for differential unit costs.
#' - **power**: Bankier (1988) compromise,
#'   \eqn{n_h \propto S_h N_h^{q}}{n_h ~ S_h * N_h^alloc_q}.
#'   The parameter `alloc_q` controls the trade-off between national precision
#'   (`alloc_q = 1`, equivalent to Neyman) and near-equal subnational CVs
#'   (`alloc_q = 0`).
#'
#' Stratum allocations are rounded to integers using the ORIC method
#' (Cont and Heidari, 2015), which preserves `sum(n) = n` while minimizing
#' rounding distortion.
#'
#' ## What the boundaries are optimal for
#'
#' Every quantity the search reads comes from `x`. The stratum standard
#' deviations \eqn{S_h} are the standard deviations of `x` inside each
#' candidate stratum, and the variance being minimized is that of the
#' estimated total or mean **of `x`**. The boundaries returned are therefore
#' optimal for estimating `x` itself.
#'
#' Surveys rarely measure `x`. It is an auxiliary variable already on the
#' frame, and the study variable `y` is what will be collected. Nothing here
#' models \eqn{E(y \mid x)}{E(y | x)} or \eqn{\mathrm{Var}(y \mid x)}{Var(y | x)}, so for any `y`
#' other than `x` these boundaries are a proxy, good in proportion to how
#' closely `y` tracks `x`. A frame with household expenditure stratified for
#' a poverty rate is the usual case, and it is a reasonable one. A frame
#' stratified on establishment size for a variable unrelated to size is not.
#'
#' Two consequences worth planning around. The reported `$cv` is the CV for
#' `x`, not for the survey's own indicators, so it bounds what the design
#' achieves only to the extent that `y` and `x` share a stratum structure.
#' And with several study variables no single `x` is optimal for all of
#' them: choose the `x` closest to the indicator that matters most, or
#' compare the boundary sets a few candidate `x` produce before committing.
#'
#' ## Design effect, response rate, and the boundaries
#'
#' The variance evaluated here is the one the rest of the package uses,
#' \eqn{\mathrm{deff}\sum_h W_h^2S_h^2(1/(n_h r)-1/N_h)}{deff sum_h W_h^2S_h^2(1/(n_h r)-1/N_h)} with response rate
#' \eqn{r}, so `cv = 0.05` means the same thing to `strata_bound()` and
#' [n_alloc()], and `$n` is a fielded sample in both. Handing `$strata` to
#' [n_alloc()] with the same `cv`, `deff`, and `resp_rate` reproduces this
#' function's continuous total.
#'
#' A *scalar* `deff` scales the variance of every candidate boundary set
#' equally, so it does not move the boundaries. It changes the `n` a `cv`
#' target needs and the `cv` a given `n` achieves. Expect the cutpoints to
#' shift only by the local search's own tolerance. Boundaries would respond
#' to a design effect that varied across strata, which a scalar argument
#' cannot express.
#'
#' The reported `$cv` is that of the integer allocation in `$strata$n`.
#' Because `strata_bound()` rounds up in `cv` mode and [n_alloc()] reports a
#' continuous optimum, the two totals differ by the rounding, not by the
#' model.
#'
#' @seealso [predict.svyplan_strata] to assign new data to strata,
#'   [n_alloc()] to distribute a sample across an existing set of strata,
#'   which accepts `$strata` directly, and [svyplan()] for reusable design
#'   defaults.
#'
#' @references
#' Dalenius, T. and Hodges, J. L. (1959). Minimum variance stratification.
#' \emph{Journal of the American Statistical Association}, 54(285), 88--101.
#'
#' Lavall\enc{é}{e}e, P. and Hidiroglou, M. (1988). On the stratification of skewed
#' populations. \emph{Survey Methodology}, 14(1), 33--43.
#'
#' Kozak, M. (2004). Optimal stratification using random search method in
#' agricultural surveys. \emph{Statistics in Transition}, 6(5), 797--806.
#'
#' Gunning, P. and Horgan, J. M. (2004). A new algorithm for the
#' construction of stratum boundaries in skewed populations.
#' \emph{Survey Methodology}, 30(2), 159--166.
#'
#' Wesolowski, J., Wieczorkowski, R. and Wojciak, W. (2021). Optimality of
#' the recursive Neyman allocation.
#' \emph{Journal of Survey Statistics and Methodology}, 10(5), 1263--1275.
#'
#' Bankier, M. D. (1988). Power allocations: determining sample sizes for
#' subnational areas. \emph{The American Statistician}, 42(3), 174--177.
#'
#' Cont, R. and Heidari, M. (2015). Optimal rounding under integer
#' constraints. \emph{arXiv preprint} arXiv:1501.00014.
#'
#' @examples
#' set.seed(867)
#' x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
#'
#' # Dalenius-Hodges (non-iterative)
#' strata_bound(x, n_strata = 4, method = "cumrootf", n = 100)
#'
#' # LH (default, iterative)
#' strata_bound(x, n_strata = 4, n = 100)
#'
#' # Bankier power allocation (compromise between national and subnational CVs)
#' strata_bound(x, n_strata = 4, n = 100, alloc = "power", alloc_q = 0.5)
#'
#' # With take-all stratum
#' strata_bound(x, n_strata = 3, n = 80, take_all_above = quantile(x, 0.95))
#'
#' # Under a clustered design and imperfect response, on the same scale
#' # n_alloc() uses
#' sb <- strata_bound(x, n_strata = 4, cv = 0.08, deff = 1.8,
#'                    resp_rate = 0.85)
#' sb
#'
#' # The strata table hands straight to n_alloc()
#' n_alloc(sb$strata, cv = 0.08, deff = 1.8, resp_rate = 0.85)
#'
#' # Or carry the design assumptions in a profile
#' plan <- svyplan(deff = 1.8, resp_rate = 0.85)
#' strata_bound(x, n_strata = 4, cv = 0.08, plan = plan)
#'
#' @encoding UTF-8
#' @export
strata_bound <- function(x, n_strata, ..., n = NULL, cv = NULL,
                         method = c("lh", "cumrootf", "geo", "kozak"),
                         alloc = c("neyman", "optimal", "proportional", "power"),
                         alloc_q = 0.5,
                         unit_cost = NULL,
                         take_all_above = NULL,
                         deff = 1,
                         resp_rate = 1,
                         n_class = NULL,
                         max_iter = NULL,
                         n_restart = NULL,
                         plan = NULL) {
  .plan <- .merge_plan_args(
    .strata_plan_defaults(plan), strata_bound, match.call(), environment()
  )
  if (!is.null(.plan)) return(do.call(strata_bound, .plan))
  .check_unused_dots(...)
  if (!is.numeric(x) || length(x) < 2L) {
    stop("'x' must be a numeric vector with at least 2 elements", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("'x' must not contain NA values", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite values", call. = FALSE)
  }

  if (missing(n_strata)) {
    stop("'n_strata' is required", call. = FALSE)
  }
  n_strata <- check_count(n_strata, "n_strata", minimum = 2L)

  method <- match.arg(method)

  has_n <- !is.null(n)
  has_cv <- !is.null(cv)
  if (has_n && has_cv) {
    stop("specify at most one of 'n' or 'cv'", call. = FALSE)
  }
  if (method %in% c("lh", "kozak") && !has_n && !has_cv) {
    stop("method '", method, "' requires 'n' or 'cv'", call. = FALSE)
  }
  if (has_n) n <- check_count(n, "n")
  if (has_cv) check_scalar(cv, "cv")

  if (!is.character(alloc)) {
    stop("'alloc' must be one of \"proportional\", \"neyman\", \"optimal\", \"power\"",
         call. = FALSE)
  }
  alloc <- match.arg(alloc)
  check_deff(deff)
  check_resp_rate(resp_rate)
  if (alloc == "power") {
    if (!is.numeric(alloc_q) || length(alloc_q) != 1L || is.na(alloc_q) || alloc_q < 0 || alloc_q > 1) {
      stop("'alloc_q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }

  x_work <- x
  L_work <- n_strata

  if (!is.null(take_all_above)) {
    if (!is.numeric(take_all_above) || length(take_all_above) != 1L || is.na(take_all_above) ||
        !is.finite(take_all_above)) {
      stop("'take_all_above' must be a finite numeric scalar", call. = FALSE)
    }
    n_take_all <- sum(x >= take_all_above)
    if (n_take_all == 0L) {
      stop("no units meet the 'take_all_above' threshold", call. = FALSE)
    }
    if (n_take_all == length(x)) {
      stop("all units meet the 'take_all_above' threshold", call. = FALSE)
    }
    L_work <- n_strata - 1L
    if (L_work < 1L) {
      stop("'n_strata' must be > 1 when 'take_all_above' is specified", call. = FALSE)
    }
    x_work <- x[x < take_all_above]
  }

  n_uniq <- length(unique(x_work))
  if (n_uniq < L_work) {
    stop("fewer unique values than requested strata", call. = FALSE)
  }

  if (method == "geo" && any(x_work <= 0)) {
    stop("method 'geo' requires all values in 'x' to be positive", call. = FALSE)
  }

  if (is.null(unit_cost)) {
    cost_h <- rep(1, n_strata)
  } else {
    if (!is.numeric(unit_cost) || anyNA(unit_cost) ||
        any(!is.finite(unit_cost)) || any(unit_cost <= 0)) {
      stop("'unit_cost' must be positive and finite", call. = FALSE)
    }
    if (!length(unit_cost) %in% c(1L, n_strata)) {
      stop(
        sprintf("'unit_cost' must have length 1 or %d (one per stratum)",
                n_strata),
        call. = FALSE
      )
    }
    cost_h <- rep_len(unit_cost, n_strata)
  }

  .check_method_controls(method, n_class, max_iter, n_restart)
  if (!is.null(n_class)) n_class <- check_count(n_class, "n_class")
  if (is.null(max_iter)) max_iter <- 200L
  if (is.null(n_restart)) n_restart <- 10 * n_strata
  max_iter <- check_count(max_iter, "max_iter")
  n_restart <- check_count(n_restart, "n_restart")

  x_sort <- sort(x_work)

  n_total <- if (has_n) n else if (has_cv) NULL else length(x) / 2
  search_target_V <- NULL
  if (has_cv && !is.null(take_all_above)) {
    N_full <- length(x)
    N_work <- length(x_work)
    search_target_V <- (cv * abs(mean(x)) * N_full / N_work)^2
  }

  # Every sampled stratum needs two units, and a take-all stratum is enumerated
  # in full, so the total is bounded below before any boundary is chosen.
  min_required <- 2L * L_work +
    if (is.null(take_all_above)) 0L else n_take_all
  if (has_n && n < min_required) {
    stop(
      "'n' is infeasible: minimum feasible is ", min_required,
      call. = FALSE
    )
  }

  search_feasible <- TRUE

  if (method == "cumrootf") {
    bk <- .strata_cumrootf(x_sort, L_work, n_class)
  } else if (method == "geo") {
    bk <- .strata_geo(x_sort, L_work)
  } else if (method == "lh") {
    target_cv <- if (has_cv) cv else NULL
    n_opt <- if (has_n) {
      if (is.null(take_all_above)) n else n - n_take_all
    } else NULL
    res <- .strata_lh(x_sort, L_work, n_opt, target_cv, alloc, alloc_q,
                       cost_h[seq_len(L_work)], max_iter,
                       deff = deff, resp_rate = resp_rate,
                       target_V = search_target_V)
    bk <- res$bk
    converged <- res$converged
    search_feasible <- isTRUE(res$feasible)
  } else {
    target_cv <- if (has_cv) cv else NULL
    n_opt <- if (has_n) {
      if (is.null(take_all_above)) n else n - n_take_all
    } else NULL
    res <- .strata_kozak(x_sort, L_work, n_opt, target_cv, alloc, alloc_q,
                          cost_h[seq_len(L_work)], max_iter, n_restart,
                          deff = deff, resp_rate = resp_rate,
                          target_V = search_target_V)
    bk <- res$bk
    converged <- res$converged
    search_feasible <- isTRUE(res$feasible)
  }

  if (!is.null(take_all_above)) {
    bk <- c(bk, take_all_above)
    bk <- sort(unique(bk))
  }

  if (length(bk) != n_strata - 1L) {
    stop(
      sprintf(
        "internal error: method '%s' returned %d boundaries; expected %d",
        method,
        length(bk),
        n_strata - 1L
      ),
      call. = FALSE
    )
  }

  take_all_strata <- if (!is.null(take_all_above)) n_strata else NULL

  bins <- findInterval(x, bk, left.open = TRUE) + 1L
  if (!is.null(take_all_strata)) {
    bins[x >= bk[take_all_strata - 1L]] <- take_all_strata
  }
  N_h <- tabulate(bins, nbins = n_strata)

  # Reached when no boundary set can satisfy the request, so the search returns
  # its starting point. Every method lands here, including the ones that do not
  # search at all.
  if (.strata_degenerate(N_h, take_all_strata)) {
    stop(
      "computed strata must contain at least two population units, ",
      "except for an explicitly take-all stratum; ",
      "reduce 'n_strata' or set 'take_all_above' to enumerate the ",
      "sparse upper values",
      call. = FALSE
    )
  }

  if (!search_feasible) {
    stop(
      sprintf(
        "no stratification into %d strata satisfies the request; method '%s' found no feasible boundary set",
        n_strata, method
      ),
      call. = FALSE
    )
  }

  if (has_n) {
    m_h <- pmin(rep(2, n_strata), N_h)
    if (!is.null(take_all_strata)) {
      m_h[take_all_strata] <- N_h[take_all_strata]
    }
    min_feasible <- sum(m_h)
    max_feasible <- sum(N_h)
    tol <- max(1e-8, 1e-8 * max(1, n))
    if (n < min_feasible - tol) {
      stop(
        "'n' is infeasible for the computed strata: minimum feasible is ",
        min_feasible,
        call. = FALSE
      )
    }
    if (n > max_feasible + tol) {
      stop(
        "'n' is infeasible for the computed strata: maximum feasible is ",
        max_feasible,
        call. = FALSE
      )
    }
  }

  if (has_cv && !has_n) {
    n_total <- .strata_n_for_cv(
      x, bk, cv, alloc, alloc_q, cost_h, take_all_idx = take_all_strata,
      deff = deff, resp_rate = resp_rate
    )
    if (!is.finite(n_total)) {
      stop("target 'cv' is unattainable for the computed strata",
           call. = FALSE)
    }
  } else if (!has_n) {
    n_total <- length(x) / 2
  }

  alloc_res <- .strata_alloc(x, bk, n_total, alloc, alloc_q, cost_h,
                             take_all_strata, deff = deff,
                             resp_rate = resp_rate)

  strata_df <- data.frame(
    stratum = seq_len(n_strata),
    lower   = alloc_res$lower,
    upper   = alloc_res$upper,
    N       = alloc_res$N_h,
    share   = alloc_res$W_h,
    sd      = alloc_res$S_h,
    mean    = alloc_res$mean_h,
    n       = if (has_cv) {
      as.integer(pmin(ceiling(alloc_res$n_h - 1e-9), alloc_res$N_h))
    } else {
      lo_i <- as.integer(pmin(2, alloc_res$N_h))
      hi_i <- as.integer(alloc_res$N_h)
      if (!is.null(take_all_strata)) {
        lo_i[take_all_strata] <- hi_i[take_all_strata]
      }
      .round_oric_bounded(alloc_res$n_h, lo_i, hi_i)
    },
    take_all = if (!is.null(take_all_above)) {
      seq_len(n_strata) %in% take_all_strata
    } else {
      rep(FALSE, n_strata)
    }
  )

  conv <- if (method %in% c("lh", "kozak")) converged else NA

  n_int <- strata_df$n
  V_int <- .strata_variance(
    alloc_res$W_h, alloc_res$S_h, n_int, alloc_res$N_h, deff, resp_rate
  )
  ybar <- .aggregate_mean(alloc_res$W_h, alloc_res$mean_h)
  cv_int <- if (ybar == 0) Inf else sqrt(V_int) / abs(ybar)

  .new_svyplan_strata(
    boundaries = bk,
    n_strata   = n_strata,
    n          = sum(strata_df$n),
    cv         = cv_int,
    strata     = strata_df,
    method     = method,
    alloc      = alloc,
    params     = list(
      N         = length(x),
      deff      = deff,
      resp_rate = resp_rate,
      unit_cost = unit_cost,
      alloc_q   = if (alloc == "power") alloc_q else NULL,
      max_iter  = if (method %in% c("lh", "kozak")) max_iter else NULL,
      n_restart = if (method == "kozak") n_restart else NULL,
      cv_continuous = alloc_res$cv,
      cv_target = if (has_cv) cv else NULL,
      take_all_above = take_all_above
    ),
    converged  = conv
  )
}

#' Reject stratification controls the chosen method cannot use
#'
#' Each control is consumed by a subset of the methods. Supplying one the
#' chosen method ignores is an error rather than a silent discard, so a
#' user cannot sit tuning a knob that does nothing.
#'
#' @keywords internal
#' @noRd
.check_method_controls <- function(method, n_class, max_iter, n_restart) {
  users <- list(
    n_class = "cumrootf",
    max_iter = c("lh", "kozak"),
    n_restart = "kozak"
  )
  supplied <- list(n_class = n_class, max_iter = max_iter, n_restart = n_restart)
  for (nm in names(users)) {
    if (!is.null(supplied[[nm]]) && !method %in% users[[nm]]) {
      stop(
        sprintf(
          "'%s' applies only to method = %s, not \"%s\"",
          nm,
          paste(sprintf('"%s"', users[[nm]]), collapse = " or "),
          method
        ),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
