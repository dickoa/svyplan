#' Rotation pattern of a repeated survey
#'
#' Declare the occasions a unit spends in and out of sample over its whole
#' life. The rotation is the design input that [design_overlap()] turns into
#' a lag profile and [design_schedule()] turns into a field manifest, and it
#' is the only place a pattern is parsed.
#'
#' @param x The rotation, in one of three forms.
#'
#'   A **compact string** in either of the two notations the literature
#'   uses, told apart by whether it contains a `0`. Without one it names the
#'   **spells** in order, starting in sample: `"4-8-4"` is the CPS pattern,
#'   four occasions in, eight out, four in, and `"5"` is a five-occasion
#'   panel with no break. The count of spells must be odd, so the life starts
#'   and ends in sample. With a `0` it is one **flag per occasion**:
#'   `"1-1-0-0-1-1"` is in for two, out for two, in for two, the same life as
#'   `"2-2-2"`. A string of all 1s is the one form both notations claim and
#'   is refused, naming the two readings and the unambiguous spelling of
#'   each.
#'
#'   A **numeric or logical vector** gives one entry per occasion of a
#'   unit's life, `0`/`FALSE` for out of sample. Positive values may differ,
#'   in which case they are the number of units still measured at that
#'   occasion of the life, which is how a design that subsamples later waves
#'   is described. No entry may exceed the first, a wave being able only to
#'   re-interview part of the cohort it recruited. The first and last entries
#'   must be positive, since an out-of-sample spell before the first or after
#'   the last interview is not part of a unit's life. A life shorter than two
#'   occasions is refused, having no lag at which any sample is shared.
#'
#'   An existing `svyplan_rotation`, which is returned unchanged.
#'
#' @return A `svyplan_rotation` object: a numeric vector of takes, one entry
#'   per occasion of a unit's life. It carries:
#' \describe{
#'   \item{`$spec`}{The input as given, so the printed form can name the
#'     pattern the way it was declared.}
#'   \item{`$life`}{Occasions from a unit's first interview to its last.}
#'   \item{`$n_occasion`}{Units in sample at any one occasion, at a steady
#'     state, which is the sum of the takes.}
#'   \item{`$unbroken`}{`TRUE` when no occasion of the life is out of
#'     sample.}
#'   \item{`$equal_take`}{`TRUE` when every in-sample occasion takes the same
#'     number of units as the first.}
#' }
#'
#' @details
#' A rotation is a finite life, not a cycle. Reading `"4-8-4"` as "four in,
#' eight out, repeat forever" gives the wrong answer, and gives it quietly.
#' Under the cycle only the cohort finishing its four-occasion stint leaves
#' each occasion, so seven of eight are retained and the consecutive overlap
#' comes out at 87.5 percent. Under the finite life two cohorts leave, one
#' reaching the end of its first spell and one reaching the end of its
#' second, so six of eight are retained and the answer is the published 75
#' percent. A compact string with an even number of spells would have to be
#' read as a cycle and is rejected for that reason.
#'
#' ## Two notations, and the one string that means both
#'
#' Rotation designs are named two ways in print. `"4-8-4"` counts occasions
#' per spell, and `"1-1-0-0-1-1"` carries one flag per occasion. Both are
#' accepted, told apart by the `0`, which no spell can be. The exception is a
#' string of all 1s, which is a valid sentence in both notations and a
#' different design in each: `"1-1-1"` is three occasions in sample as a
#' pattern and one in, one out, one in as spells, whose consecutive overlaps
#' are 2/3 and 0. Neither reading is given precedence, because the wrong one
#' is invisible in the answer. Write `"3"` or `"1-0-1"`, which each say one
#' thing only.
#'
#' @references
#' U.S. Census Bureau. *Current Population Survey: Design and Methodology*,
#' Technical Paper 77. The 4-8-4 rotation.
#'
#' Lynn, P. (2012). *Longitudinal Survey Methods for the Household Finance
#' and Consumption Survey*. Report to the European Central Bank. The
#' one-flag-per-occasion notation.
#'
#' @seealso [design_overlap()] for the overlap the rotation produces,
#'   [design_schedule()] for the field manifest it produces, and
#'   [plot.svyplan_overlap()] for the rotation chart.
#'
#' @family repeated survey functions
#'
#' @examples
#' # The CPS rotation, declared once
#' cps <- design_rotation("4-8-4")
#' cps
#'
#' # The overlap and the manifest both read the same object
#' design_overlap(cps)[1]
#' design_overlap(cps, max_lag = 12)[12]
#'
#' # One flag per occasion, the other notation in use
#' design_rotation("1-1-0-0-1-1")
#'
#' # A take that halves at later waves
#' as.data.frame(design_rotation(c(1, 1, 0.5, 0.5)))
#'
#' @export
design_rotation <- function(x) {
  if (inherits(x, "svyplan_rotation")) {
    return(x)
  }
  if (inherits(x, "svyplan_overlap")) {
    stop(
      "an overlap is what a rotation produces, not a rotation; pass the rotation itself, or design_rotation(attr(x, \"rotation\"))",
      call. = FALSE
    )
  }
  w <- .resolve_rotation(x)
  positive <- w[w > 0]
  structure(
    w,
    spec = x,
    life = length(w),
    n_occasion = sum(w),
    unbroken = !any(w == 0),
    equal_take = length(unique(signif(positive, 14L))) == 1L,
    class = c("svyplan_rotation", "numeric")
  )
}

#' Resolve a rotation to one take per occasion of a unit's life
#'
#' The compact string is a constructor, never the representation. It is
#' expanded here and every downstream calculation reads the vector. That is
#' what keeps a two-spell design from being confused with a cycle, since a
#' vector has an end and a cycle does not.
#' @keywords internal
#' @noRd
.resolve_rotation <- function(rotation) {
  # The compact spec is expanded first and then validated with everything
  # else, so a spec cannot reach the arithmetic through a shorter path than a
  # vector does.
  if (is.character(rotation)) {
    rotation <- .expand_rotation_spec(rotation)
  }
  if (is.logical(rotation)) {
    rotation <- as.numeric(rotation)
  }
  if (!is.numeric(rotation) || length(rotation) == 0L || anyNA(rotation) ||
      !all(is.finite(rotation))) {
    stop(
      "'x' must be a compact string like \"4-8-4\" or a finite numeric vector, one entry per occasion of a unit's life",
      call. = FALSE
    )
  }
  if (any(rotation < 0)) {
    stop("a rotation must not take a negative number of units", call. = FALSE)
  }
  if (rotation[1L] <= 0 || rotation[length(rotation)] <= 0) {
    stop(
      "a rotation must start and end in sample. An out-of-sample spell before the first or after the last interview is not part of a unit's life",
      call. = FALSE
    )
  }
  if (length(rotation) < 2L) {
    stop(
      "a rotation must cover at least two occasions. A life of one occasion has no lag at which any sample can be shared",
      call. = FALSE
    )
  }
  if (any(rotation > rotation[1L] * (1 + 1e-9))) {
    stop(
      "a rotation takes more units at a later occasion than it recruits at the first. A wave can only re-interview part of the cohort it started with",
      call. = FALSE
    )
  }
  if (sum(rotation) <= 0) {
    stop("a rotation must place units in sample", call. = FALSE)
  }
  as.numeric(rotation)
}

#' Expand a compact rotation spec such as "4-8-4" or "1-1-0-0-1-1"
#'
#' Two notations are in published use and they look alike. Spell lengths
#' ("4-8-4", the CPS convention) count occasions per spell, alternating from
#' in sample. A per-occasion pattern ("1-1-0-0-1-1", Lynn 2012) carries one
#' flag per occasion. A `0` can only be a pattern, since a spell of no
#' occasions is not a spell, so that case is read as one. A string of all 1s
#' is the one form both notations claim, and it is refused rather than
#' resolved by precedence: "1-1-1" is three occasions in sample under one
#' reading and in-out-in under the other, and each has an unambiguous
#' spelling in the other notation.
#'
#' Under spell lengths an even count would end on an out-of-sample spell and
#' could only mean a repeating cycle, which is refused for the reason the
#' Details of [design_rotation()] give.
#' @keywords internal
#' @noRd
.expand_rotation_spec <- function(spec) {
  if (length(spec) != 1L || is.na(spec) || !nzchar(spec)) {
    stop("'x' must be a single non-empty string", call. = FALSE)
  }
  parts <- strsplit(spec, "-", fixed = TRUE)[[1L]]
  n <- suppressWarnings(as.numeric(parts))
  if (anyNA(n) || any(n < 0) || any(n != trunc(n))) {
    stop(
      sprintf(
        "'%s' is not a rotation spec. Use whole spell lengths separated by '-', as in \"4-8-4\", or one 0/1 flag per occasion, as in \"1-1-0-0-1-1\"",
        spec
      ),
      call. = FALSE
    )
  }
  if (any(n == 0)) {
    # a zero-length spell is not a spell, so this can only be the
    # per-occasion pattern, and every entry has to be a flag
    if (!all(n %in% c(0, 1))) {
      stop(
        sprintf(
          "'%s' reads as one flag per occasion, because of the 0, but takes a value other than 0 or 1; give spell lengths without any 0, as in \"4-8-4\", or a vector for a take that varies over the life",
          spec
        ),
        call. = FALSE
      )
    }
    return(n)
  }
  if (length(n) >= 2L && all(n == 1)) {
    stop(.spec_ambiguous_message(spec, length(n)), call. = FALSE)
  }
  if (length(n) %% 2L == 0L) {
    stop(
      sprintf(
        "'%s' ends on an out-of-sample spell, which can only mean a repeating cycle. A rotation is a finite life, so give the closing in-sample spell, as in \"%s-%g\"",
        spec, spec, n[1L]
      ),
      call. = FALSE
    )
  }
  rep(rep(c(1, 0), length.out = length(n)), times = n)
}

#' Name both readings of an all-1s spec, and the unambiguous form of each
#'
#' Each reading is unambiguous in the other notation, so the remedy is a
#' translation rather than a new argument: "1-1-1" as a pattern is the spell
#' spec "3", and as spell lengths it is the pattern "1-0-1". An even count
#' has no spell reading at all, so only the pattern remedy is offered.
#' @keywords internal
#' @noRd
.spec_ambiguous_message <- function(spec, k) {
  as_pattern <- sprintf(
    "as one flag per occasion it is %d consecutive occasions in sample, which is \"%d\"",
    k, k
  )
  if (k %% 2L == 0L) {
    return(sprintf(
      "'%s' is ambiguous: %s; as spell lengths it ends on an out-of-sample spell, which can only mean a repeating cycle",
      spec, as_pattern
    ))
  }
  sprintf(
    "'%s' is ambiguous: %s; as spell lengths it is one occasion in, one out, and so on, which is \"%s\"",
    spec, as_pattern,
    paste(rep(c("1", "0"), length.out = k), collapse = "-")
  )
}
