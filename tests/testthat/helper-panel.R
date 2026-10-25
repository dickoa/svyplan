## Brute-force launch of a rotating panel, shared by the panel and plot tests

# One row per cohort in sample at each period, with the visit it is at,
# enumerated from entry occasions and stages without the package's builder.
launch_oracle <- function(w, start, n_period) {
  L <- length(w)
  rows <- list()
  add <- function(id, t, visit) {
    rows[[length(rows) + 1L]] <<- data.frame(id = id, period = t, visit = visit)
  }
  first <- if (identical(start, "immediate")) 2L else 1L
  for (e in seq.int(first, n_period)) {
    for (t in seq.int(e, n_period)) {
      stage <- t - e + 1L
      if (stage <= L && w[[stage]] > 0) {
        add(paste0("c", e), t, sum(w[seq_len(stage)] > 0))
      }
    }
  }
  if (identical(start, "immediate")) {
    for (s in seq_len(L)) {
      own <- w[s:L]
      for (t in seq_along(own)) {
        if (t <= n_period && own[[t]] > 0) {
          add(paste0("s", s), t, sum(own[seq_len(t)] > 0))
        }
      }
    }
  }
  do.call(rbind, rows)
}
