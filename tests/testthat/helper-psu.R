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
    size <- pmax(round(size * frame$N[h] / sum(size)), take)
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

.psu_random_register <- function(k) {
  set.seed(99)
  for (i in seq_len(k)) {
    H <- sample(2:6, 1)
    take <- sample(c(5, 10, 20, 30), 1)
    psu <- do.call(rbind, lapply(seq_len(H), function(h) {
      m <- sample(c(1, 2, 5, 20, 60), 1, prob = c(0.1, 0.1, 0.2, 0.3, 0.3))
      data.frame(
        stratum = LETTERS[h],
        N = pmax(round(rlnorm(m, log(400), runif(1, 0.2, 1.5))), take)
      )
    }))
    frame <- data.frame(
      stratum = LETTERS[seq_len(H)],
      N = as.numeric(tapply(psu$N, psu$stratum, sum)[LETTERS[seq_len(H)]]),
      n_per_psu = take
    )
    measures <- data.frame(
      stratum = frame$stratum, name = "y", p = runif(H, 0.1, 0.6),
      icc_psu = runif(H, 0, 0.15)
    )
    targets <- data.frame(name = "y", cv = runif(1, 0.02, 0.12))
  }
  list(frame = frame, psu = psu, measures = measures, targets = targets)
}

.psu_wide_register <- function(seed) {
  set.seed(seed)
  H <- sample(2:6, 1)
  take <- sample(c(5, 10, 20, 30), 1)
  psu <- do.call(rbind, lapply(seq_len(H), function(h) {
    m <- sample(c(1, 2, 3, 4, 5, 8, 20, 60), 1)
    data.frame(
      stratum = LETTERS[h],
      N = pmax(round(rlnorm(m, log(400), runif(1, 0.2, 1.8))), take)
    )
  }))
  frame <- data.frame(
    stratum = LETTERS[seq_len(H)],
    N = as.numeric(tapply(psu$N, psu$stratum, sum)[LETTERS[seq_len(H)]]),
    n_per_psu = take
  )
  measures <- data.frame(
    stratum = frame$stratum, name = "y", p = runif(H, 0.1, 0.6),
    icc_psu = runif(H, 0, 0.3)
  )
  targets <- data.frame(name = "y", cv = runif(1, 0.02, 0.15))
  list(frame = frame, psu = psu, measures = measures, targets = targets)
}

.psu_remainder_crossing <- function(fit, cutoff = 1) {
  cutoff <- rep_len(cutoff, nrow(fit$detail))
  out <- logical(nrow(fit$psu))
  for (j in seq_len(nrow(fit$detail))) {
    h <- fit$detail$stratum[j]
    i <- which(fit$psu$stratum == h & !fit$psu$certainty)
    k <- fit$detail$n_psu_draw[j]
    if (k > 0 && length(i)) {
      out[i] <- k * fit$psu$N[i] / sum(fit$psu$N[i]) >=
        cutoff[j] - 100 * .Machine$double.eps
    }
  }
  out
}
