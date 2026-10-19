## Brute-force menus of a register stratum, shared by the menu tests

# The state the R code gives at f: the threshold classification, then the
# closure under the whole-PSU draw, named as .psu_menu_part() names it.
menu_state_at <- function(f, sizes, N, take, cutoff, forced, m, zone_key) {
  n_h <- f * N
  held <- .psu_classify(sizes, cutoff * take / (n_h / N), forced)$certain
  closed <- .psu_draw_closure(n_h, N, take, sizes, held, cutoff, m, zone_key)
  counts <- .psu_stratum_counts(n_h, N, take, sizes, closed$certain, m)
  paste(paste(which(closed$certain), collapse = ","),
        sum(counts$per_psu), counts$n_psu_draw, sep = "|")
}

# Every state of the range: the R closure at every event and every midpoint
# between events. A draw changes where f N - Q crosses a multiple of the
# take, on the grid j / N for a whole take, a take q_i at j / M_i, and the
# threshold at cutoff take / M_i.
menu_oracle <- function(sizes, N, take, f_lo, f_hi, cutoff = 1, forced = NULL,
                        m = 0L, zone_key = NULL) {
  ev <- c(seq(ceiling(f_lo * N), floor(f_hi * N)) / N,
          unlist(lapply(sizes, function(x) {
            seq(ceiling(f_lo * x), floor(f_hi * x)) / x
          })),
          cutoff * take / sizes)
  ev <- sort(unique(c(f_lo, f_hi, ev[ev > f_lo & ev < f_hi])))
  pts <- c(ev, (utils::head(ev, -1) + utils::tail(ev, -1)) / 2)
  unique(vapply(pts, menu_state_at, "", sizes = sizes, N = N, take = take,
                cutoff = cutoff, forced = forced, m = m, zone_key = zone_key))
}

menu_register <- function(seed) {
  set.seed(seed)
  k <- sample(c(4, 8, 15), 1)
  take <- sample(c(5, 10), 1)
  sizes <- pmax(round(stats::rlnorm(k, log(40), 0.9)), take)
  if (k > 4) sizes[1:2] <- sizes[1:2] * 6
  list(sizes = sizes, take = take)
}
