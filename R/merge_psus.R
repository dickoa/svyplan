#' Merge small PSUs up to a minimum size
#'
#' @description
#' `merge_psus()` recodes primary sampling unit (PSU) ids so that every PSU
#' holds at least `min_size` units, merging small PSUs with their neighbours
#' in the current row order. A PSU register for [n_alloc()] must not hold a
#' PSU smaller than its stratum's take, and merging small enumeration areas
#' before the register is built is the usual remedy.
#'
#' It is a vector function. Apply it within each stratum, after sorting the
#' frame so that neighbours in row order are neighbours on the ground, and
#' build the register from the merged ids, so the register and the frame
#' agree.
#'
#' @param id The PSU id of each row. Any atomic type.
#' @param min_size The smallest size a PSU may have, a positive number.
#'   Usually the largest take `n_per_psu`.
#' @param size The size each row contributes, for a frame with one row per
#'   PSU or with weighted rows. `NULL` (default) counts rows.
#'
#' @return A vector like `id`, each row's merged PSU id. A merged PSU takes
#'   the id of its first PSU in row order, so its type and class are those of
#'   `id` and every merged PSU traces to one source PSU.
#'
#' @details
#' PSUs are taken in order of first appearance. A PSU of at least
#' `min_size` stays as it is. A run of smaller PSUs between two of them is
#' merged in order, closing a merged PSU as soon as its size reaches
#' `min_size`, and a remainder below `min_size` at the end of the run joins
#' the merged PSU before it. A run that is below `min_size` in total joins
#' the PSU before it, or the one after it when it opens the stratum. A
#' stratum below `min_size` in total becomes one PSU, which [n_alloc()]
#' refuses as smaller than its take unless the register flags it as
#' certainty.
#'
#' Merging is by order only, so neighbours in row order are neighbours on
#' the ground only as far as the sort makes them so. Valliant, Dever and
#' Kreuter (2018, sec. 10.3) form PSUs to a minimum size from contiguous
#' areas.
#'
#' Called on a whole frame, `merge_psus()` sees one stratum and can merge
#' PSUs of different strata. Apply it per stratum, with `split()` or a
#' grouped `dplyr::mutate()`. A merge that crosses strata gives one id to
#' PSUs of two strata, which [n_alloc()] refuses as a duplicated `psu_id`.
#'
#' @references
#' Valliant, R., Dever, J. A., & Kreuter, F. (2018). *Practical Tools for
#'   Designing and Weighting Survey Samples* (2nd ed.). Springer. Section
#'   10.3.
#'
#' @seealso [n_alloc()] for the register the merged ids feed.
#'
#' @examples
#' listing <- data.frame(
#'   stratum = rep(c("urban", "rural"), c(7, 6)),
#'   ea = c(sprintf("U%d", 1:7), sprintf("R%d", 1:6)),
#'   households = c(120, 6, 4, 95, 110, 3, 80, 70, 8, 60, 5, 2, 90)
#' )
#'
#' # One row per enumeration area, already in geographic order
#' by_stratum <- split(listing, listing$stratum)
#' listing$psu_id <- unsplit(lapply(by_stratum, function(d) {
#'   merge_psus(d$ea, min_size = 10, size = d$households)
#' }), listing$stratum)
#' listing
#'
#' # The register n_alloc() reads
#' aggregate(households ~ stratum + psu_id, listing, sum)
#'
#' @export
merge_psus <- function(id, min_size, size = NULL) {
  if (!is.atomic(id) || is.null(id)) {
    stop("'id' must be an atomic vector of PSU ids", call. = FALSE)
  }
  if (anyNA(id)) {
    stop("'id' must not have missing values", call. = FALSE)
  }
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) ||
        !is.finite(min_size) || min_size <= 0) {
    stop("'min_size' must be a single positive number", call. = FALSE)
  }
  n <- length(id)
  if (is.null(size)) {
    size <- rep(1, n)
  } else if (!is.numeric(size) || length(size) != n || anyNA(size) ||
               any(!is.finite(size)) || any(size < 0)) {
    stop(
      sprintf(
        "'size' must be a nonnegative number for every element of 'id' (it has length %d, 'id' has %d)",
        length(size), n
      ),
      call. = FALSE
    )
  }
  if (n == 0L) {
    return(id)
  }

  key <- match(id, unique(id))
  psu_size <- as.vector(tapply(size, key, sum))
  group <- .merge_psu_groups(psu_size, min_size)
  first <- match(seq_along(psu_size), key)
  id[first[group][key]]
}

#' Group PSUs in order so every group reaches a minimum size
#'
#' Returns, for each PSU in order, the position of the first PSU of its
#' group. PSUs at or above `min_size` stay alone. Each run of smaller PSUs
#' is cut greedily, its short remainder joining the group before it, and a
#' run short in total joins the PSU before it (after it when it is first).
#' @keywords internal
#' @noRd
.merge_psu_groups <- function(psu_size, min_size) {
  k <- length(psu_size)
  group <- seq_len(k)
  big <- psu_size >= min_size
  runs <- rle(big)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  for (r in which(!runs$values)) {
    run <- starts[r]:ends[r]
    start <- run[1L]
    total <- 0
    for (i in run) {
      group[i] <- start
      total <- total + psu_size[i]
      if (total >= min_size) {
        start <- i + 1L
        total <- 0
      }
    }
    # A remainder exists whenever the run did not close on its last PSU,
    # whatever its size: a zero-size PSU still has to join a group.
    if (start <= ends[r]) {
      short <- run[run >= start]
      if (start > run[1L]) {
        group[short] <- group[start - 1L]
      } else if (run[1L] > 1L) {
        group[short] <- group[run[1L] - 1L]
      } else if (ends[r] < k) {
        # The run opens the stratum: the PSU after it joins the run.
        group[ends[r] + 1L] <- run[1L]
      }
    }
  }
  group
}
