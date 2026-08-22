.psu_fixture <- function(seed = 4L, sdlog = 0.35, large = 6L, take = 12) {
  frame <- data.frame(
    stratum = c("A", "B", "C", "D"),
    N = c(60000, 40000, 25000, 15000),
    n_per_psu = take,
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  psu <- do.call(rbind, lapply(seq_len(nrow(frame)), function(h) {
    k <- c(300, 200, 140, 90)[h]
    size <- round(rlnorm(k, log(frame$N[h] / k), sdlog))
    if (large > 0L) size[seq_len(large)] <- size[seq_len(large)] * 25
    size <- pmax(round(size * frame$N[h] / sum(size)), 1)
    size[length(size)] <- size[length(size)] + frame$N[h] - sum(size)
    data.frame(
      psu_id = sprintf("%s%03d", frame$stratum[h], seq_len(k)),
      stratum = frame$stratum[h], N = size, stringsAsFactors = FALSE
    )
  }))
  measures <- data.frame(
    stratum = rep(frame$stratum, 2),
    name = rep(c("literacy", "income"), each = nrow(frame)),
    p = c(0.62, 0.55, 0.48, 0.41, rep(NA_real_, 4)),
    mean = c(rep(NA_real_, 4), 520, 430, 380, 310),
    sd = c(rep(NA_real_, 4), 210, 180, 160, 140),
    icc_psu = rep(c(0.05, 0.08), each = nrow(frame)),
    stringsAsFactors = FALSE
  )
  targets <- data.frame(name = c("literacy", "income"), cv = c(0.04, 0.05))
  list(frame = frame, psu = psu, measures = measures, targets = targets)
}
