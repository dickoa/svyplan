#' Optimal Strata Boundaries
#'
#' Determine where to cut a continuous stratification variable to form
#' optimal strata. Supports four methods: cumulative root frequency
#' (Dalenius-Hodges), geometric progression, Lavallée-Hidiroglou
#' iterative, and Kozak random search.
#'
#' @param x Numeric vector: stratification variable values. Must not
#'   contain `NA`.
#' @param n_strata Integer: number of strata (including take-all if
#'   `certain` is specified). Must be >= 2.
#' @param n Target total sample size. Specify at most one of `n` or `cv`.
#'   Required for methods `"lh"` and `"kozak"`.
#' @param cv Target coefficient of variation. Specify at most one of `n`
#'   or `cv`. Required for methods `"lh"` and `"kozak"`.
#' @param method Stratification method: `"cumrootf"` (Dalenius-Hodges),
#'   `"geo"` (geometric), `"lh"` (Lavallée-Hidiroglou), or `"kozak"`
#'   (random search). Default `"lh"`.
#' @param alloc Allocation rule: `"proportional"`, `"neyman"`,
#'   `"optimal"`, or `"power"` (Bankier compromise). Default `"neyman"`.
#'   See Details.
#' @param q Bankier power parameter, used only when `alloc = "power"`.
#'   Numeric scalar in \eqn{[0, 1]}. At `q = 1` the allocation equals
#'   Neyman; at `q = 0` it yields near-equal subnational CVs.
#'   Default 0.5.
#' @param cost Per-stratum unit costs, ordered from lowest to highest
#'   stratum. Scalar (equal costs) or vector of length `n_strata`.
#'   Default `NULL` (equal unit costs, \eqn{c_h = 1} for all strata),
#'   in which case `"optimal"` and `"neyman"` coincide.
#' @param certain Take-all threshold. Units with `x >= certain` form a
#'   census stratum.
#' @param nclass Number of histogram bins for `"cumrootf"`. Default
#'   `NULL` (Freedman-Diaconis rule).
#' @param maxiter Maximum iterations for `"lh"` and `"kozak"`. Default
#'   200.
#' @param niter Random restarts for `"kozak"`. Default `NULL` (= 10 *
#'   `n_strata`).
#'
#' @return A `svyplan_strata` object with components:
#' \describe{
#'   \item{boundaries}{Numeric vector of cutpoints (length `n_strata - 1`).}
#'   \item{n_strata}{Number of strata.}
#'   \item{n}{Total sample size.}
#'   \item{cv}{Achieved coefficient of variation.}
#'   \item{strata}{Data frame with per-stratum summaries.}
#'   \item{method}{Algorithm used.}
#'   \item{alloc}{Allocation method name (character).}
#'   \item{params}{List of additional parameters.}
#'   \item{converged}{Logical (for iterative methods).}
#' }
#'
#' @details
#' The four methods differ in approach:
#'
#' - **cumrootf**: Dalenius-Hodges (1959) cumulative root frequency rule.
#'   Non-iterative, does not require `n` or `cv`.
#' - **geo**: Gunning-Horgan (2004) geometric progression. Non-iterative,
#'   requires `x > 0`.
#' - **lh**: Lavallée-Hidiroglou (1988) iterative coordinate-wise
#'   optimization. Requires `n` or `cv`.
#' - **kozak**: Kozak (2004) random search, escapes local minima. Requires
#'   `n` or `cv`.
#'
#' When `n` or `cv` is available, `"lh"` (the default) directly
#' minimizes the objective and is fast even on large frames. Use `"kozak"`
#' for highly skewed or multimodal data where coordinate-wise optimization
#' may get trapped in local minima. The non-iterative methods (`"cumrootf"`,
#' `"geo"`) are useful for quick exploratory stratification or when neither
#' `n` nor `cv` is known yet.
#'
#' Allocation is controlled by the `alloc` parameter. Four methods are
#' available:
#' - **proportional**: \eqn{n_h \propto N_h}{n_h ~ N_h}.
#' - **neyman**: \eqn{n_h \propto N_h S_h}{n_h ~ N_h * S_h}.
#'   Minimizes the national CV when unit costs are equal.
#' - **optimal**: \eqn{n_h \propto N_h S_h / \sqrt{c_h}}{n_h ~ N_h * S_h / sqrt(c_h)}.
#'   Accounts for differential unit costs.
#' - **power**: Bankier (1988) compromise,
#'   \eqn{n_h \propto S_h N_h^q}{n_h ~ S_h * N_h^q}.
#'   The parameter `q` controls the trade-off between national precision
#'   (`q = 1`, equivalent to Neyman) and near-equal subnational CVs
#'   (`q = 0`).
#'
#' Stratum allocations `n_h` are rounded to integers using the ORIC method
#' (Cont and Heidari, 2014), which preserves `sum(n_h) = n` while minimizing
#' rounding distortion.
#'
#' @seealso [predict.svyplan_strata] to assign new data to strata.
#'
#' @references
#' Dalenius, T. and Hodges, J. L. (1959). Minimum variance stratification.
#' \emph{Journal of the American Statistical Association}, 54(285), 88--101.
#'
#' Lavallee, P. and Hidiroglou, M. (1988). On the stratification of skewed
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
#' Cont, R. and Heidari, M. (2014). Optimal rounding under integer
#' constraints. \emph{arXiv preprint} arXiv:1501.00014.
#'
#' @examples
#' set.seed(42)
#' x <- rlnorm(500, meanlog = 6, sdlog = 1.5)
#'
#' # Dalenius-Hodges (non-iterative)
#' strata_bound(x, n_strata = 4, method = "cumrootf", n = 100)
#'
#' # LH (default, iterative)
#' strata_bound(x, n_strata = 4, n = 100)
#'
#' # Bankier power allocation (compromise between national and subnational CVs)
#' strata_bound(x, n_strata = 4, n = 100, alloc = "power", q = 0.5)
#'
#' # With take-all stratum
#' strata_bound(x, n_strata = 3, n = 80, certain = quantile(x, 0.95))
#'
#' @export
strata_bound <- function(x, n_strata = 3L, n = NULL, cv = NULL,
                         method = "lh",
                         alloc = "neyman",
                         q = 0.5,
                         cost = NULL,
                         certain = NULL,
                         nclass = NULL,
                         maxiter = 200L,
                         niter = NULL) {
  if (!is.numeric(x) || length(x) < 2L) {
    stop("'x' must be a numeric vector with at least 2 elements", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("'x' must not contain NA values", call. = FALSE)
  }

  n_strata <- as.integer(n_strata)
  if (length(n_strata) != 1L || is.na(n_strata) || n_strata < 2L) {
    stop("'n_strata' must be an integer >= 2", call. = FALSE)
  }

  method <- match.arg(method, c("cumrootf", "geo", "lh", "kozak"))

  has_n <- !is.null(n)
  has_cv <- !is.null(cv)
  if (has_n && has_cv) {
    stop("specify at most one of 'n' or 'cv'", call. = FALSE)
  }
  if (method %in% c("lh", "kozak") && !has_n && !has_cv) {
    stop("method '", method, "' requires 'n' or 'cv'", call. = FALSE)
  }
  if (has_n) check_scalar(n, "n")
  if (has_cv) check_scalar(cv, "cv")

  if (!is.character(alloc)) {
    stop("'alloc' must be one of \"proportional\", \"neyman\", \"optimal\", \"power\"",
         call. = FALSE)
  }
  alloc <- match.arg(alloc, c("proportional", "neyman", "optimal", "power"))
  if (alloc == "power") {
    if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < 0 || q > 1) {
      stop("'q' must be a numeric scalar in [0, 1]", call. = FALSE)
    }
  }

  certain_idx <- NULL
  x_work <- x
  L_work <- n_strata

  if (!is.null(certain)) {
    if (!is.numeric(certain) || length(certain) != 1L || is.na(certain)) {
      stop("'certain' must be a numeric scalar", call. = FALSE)
    }
    n_certain <- sum(x >= certain)
    if (n_certain == 0L) {
      stop("no units meet the 'certain' threshold", call. = FALSE)
    }
    if (n_certain == length(x)) {
      stop("all units meet the 'certain' threshold", call. = FALSE)
    }
    L_work <- n_strata - 1L
    if (L_work < 1L) {
      stop("'n_strata' must be > 1 when 'certain' is specified", call. = FALSE)
    }
    x_work <- x[x < certain]
  }

  n_uniq <- length(unique(x_work))
  if (n_uniq < L_work) {
    stop("fewer unique values than requested strata", call. = FALSE)
  }

  if (method == "geo" && any(x_work <= 0)) {
    stop("method 'geo' requires all values in 'x' to be positive", call. = FALSE)
  }

  if (is.null(cost)) {
    cost_h <- rep(1, n_strata)
  } else {
    if (!is.numeric(cost) || anyNA(cost) || any(cost <= 0)) {
      stop("'cost' must be positive numeric", call. = FALSE)
    }
    cost_h <- rep_len(cost, n_strata)
  }

  if (is.null(niter)) niter <- 10L * n_strata
  maxiter <- as.integer(maxiter)
  niter <- as.integer(niter)

  x_sort <- sort(x_work)

  n_total <- if (has_n) n else if (has_cv) NULL else length(x) / 2

  if (method == "cumrootf") {
    bk <- .strata_cumrootf(x_sort, L_work, nclass)
  } else if (method == "geo") {
    bk <- .strata_geo(x_sort, L_work)
  } else if (method == "lh") {
    target_cv <- if (has_cv) cv else NULL
    n_opt <- if (has_n) n else NULL
    res <- .strata_lh(x_sort, L_work, n_opt, target_cv, alloc, q,
                       cost_h[seq_len(L_work)], maxiter)
    bk <- res$bk
    converged <- res$converged
  } else {
    target_cv <- if (has_cv) cv else NULL
    n_opt <- if (has_n) n else NULL
    res <- .strata_kozak(x_sort, L_work, n_opt, target_cv, alloc, q,
                          cost_h[seq_len(L_work)], maxiter, niter)
    bk <- res$bk
    converged <- res$converged
  }

  if (!is.null(certain)) {
    bk <- c(bk, certain)
    bk <- sort(unique(bk))
  }

  if (length(bk) != n_strata - 1L) {
    bk <- bk[seq_len(n_strata - 1L)]
  }

  if (has_cv && !has_n) {
    n_total <- .strata_n_for_cv(x, bk, cv, alloc, q, cost_h)
  } else if (!has_n) {
    n_total <- length(x) / 2
  }

  certain_strata <- if (!is.null(certain)) n_strata else NULL

  alloc_res <- .strata_alloc(x, bk, n_total, alloc, q, cost_h, certain_strata)

  strata_df <- data.frame(
    stratum = seq_len(n_strata),
    lower   = alloc_res$lower,
    upper   = alloc_res$upper,
    N_h     = alloc_res$N_h,
    W_h     = round(alloc_res$W_h, 6L),
    S_h     = round(alloc_res$S_h, 4L),
    n_h     = .round_oric(alloc_res$n_h),
    certain = if (!is.null(certain)) {
      seq_len(n_strata) %in% certain_strata
    } else {
      rep(FALSE, n_strata)
    }
  )

  conv <- if (method %in% c("lh", "kozak")) converged else NA

  .new_svyplan_strata(
    boundaries = bk,
    n_strata   = n_strata,
    n          = sum(strata_df$n_h),
    cv         = alloc_res$cv,
    strata     = strata_df,
    method     = method,
    alloc      = alloc,
    params     = list(
      N       = length(x),
      q       = if (alloc == "power") q else NULL,
      maxiter = if (method %in% c("lh", "kozak")) maxiter else NULL,
      niter   = if (method == "kozak") niter else NULL
    ),
    converged  = conv
  )
}
