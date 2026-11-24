## Summaries of single-indicator, multi-indicator and cluster sizes

#' Summarize sample size and cluster results
#'
#' `print()` gives the answer of [n_prop()], [n_mean()], [n_ratio()],
#' [n_multi()], [n_cluster()], [n_change()], and [n_pooled()]. `summary()` gives
#' what follows from it,
#' read at the design that is fielded: the whole sample size or the whole
#' stage sizes in `$operational`, beside the continuous solution they round.
#'
#' @param object A result of [n_prop()], [n_mean()], [n_ratio()],
#'   [n_multi()], [n_cluster()], [n_change()], or [n_pooled()].
#' @param x A `summary.svyplan_n` or `summary.svyplan_cluster` object.
#' @param ... Additional arguments are not supported and produce an error.
#'
#' @return `summary()` returns a list of class `summary.svyplan_n` or
#'   `summary.svyplan_cluster` whose numbers are unformatted. Formatting is
#'   applied only by the print methods, which return the summary invisibly.
#'   A size for a proportion, mean or ratio has `kind = "scalar"` and the
#'   fields `type`, `method`, `target` (its `scale` and `value`),
#'   `estimand`, `binding`, `min_cases`, `sizes`, `precision`,
#'   `expected_cases`, `interval`, `level` and `assumptions`. A
#'   multi-indicator size has `kind = "multi"` and the fields `sizes`,
#'   `binding`, `binding_cause`, `n_domain_max`, `domains`, `domain_values`,
#'   `indicators` and `assumptions`. A single-indicator cluster size has `kind = "single"` and
#'   the fields `stages`, `question`, `design`, `field`, `eligible`, `cost`
#'   and `assumptions`. A multi-indicator cluster size has `kind = "multi"`
#'   and the fields `stages`, `question`, `design`, `domains`,
#'   `domain_values`, `total`, `worst`, `binding`, `binding_cause`,
#'   `indicators`, `cost` and `assumptions`. A change summary has
#'   `kind = "change"` and sizes in `n`. A pooled summary has
#'   `kind = "pooled"` and sizes in `sizes`. Both report the
#'   target, precision at the whole size, and assumptions.
#'
#' @details
#' ## Field size
#'
#' A one-stage size is fielded as `ceiling(n)`. Separate domain quotas are
#' fielded one by one, so their total is the sum of the whole quotas. Under
#' natural incidence the whole total is fielded and each domain is expected
#' at that total times its share, which need not be whole. A cluster design
#' is fielded as `$operational`, and a cluster design with domains as each
#' row of `$operational$domains`.
#'
#' ## Precision at the field size
#'
#' The precision is that of the matching `prec_*()` function evaluated at
#' the field size, so it can be slightly better than the target the
#' continuous size meets exactly. Each indicator's target and achieved
#' precision are reported on the scale the target was stated in (`moe`,
#' `rmoe` or `cv`). A row with only a respondent floor has no precision and
#' reads `NA`. An interval method other than Wald adds its limits, at the
#' confidence level the size was planned for (`1 - alpha`).
#'
#' ## Derived sizes
#'
#' Expected respondents are the field size times the response and
#' eligibility rates, and eligible respondents multiply that by
#' `eligible_per_unit`. The effective size divides the units analyzed by the
#' design effect. For a cluster design the design effect and the effective
#' size are the planning quantities of [design_effect()] and
#' [effective_n()], which leave out finite population corrections. A
#' multi-indicator cluster design has no single design effect and the
#' summary reports none.
#'
#' ## Binding indicator
#'
#' The binding indicator of `n_multi()` is the row whose requirement sets
#' the size, and `binding_cause` says what set it: `"precision"`,
#' `"min_cases"`, `"min_respondents"` or `"min_n_domain"`. For a
#' multi-indicator cluster design it is the row with the largest ratio of
#' achieved to target cv at the field design, which need not be the row with
#' the worst cv. With domains it is that row in the domain with the worst cv
#' at the field design.
#'
#' ## Domains
#'
#' A domain is labeled by its value when one column defines domains, and by
#' `column=value` pairs joined by commas when several do. A value holding a
#' comma, an equals sign, a quote or a backslash is quoted, so two domains
#' never share a label. `domain_values` holds the domain columns, one row
#' per row of `domains`.
#'
#' ## Assumptions of a multi-indicator cluster design
#'
#' The `indicators` table carries each row's design quantities:
#' `unit_relvar`, `icc_1` through `icc_{S-1}`, `var_ratio_1` through
#' `var_ratio_{S-1}`, and the response rates, with `alpha` when a row
#' states a `moe` or `rmoe`. A quantity equal on every row is also
#' listed in `assumptions`, and the print shows the ones that vary in a
#' table of their own.
#'
#' Results of [n_alloc()] have their own summary, described in
#' [summary.svyplan_alloc].
#'
#' @examples
#' summary(n_prop(p = 0.3, moe = 0.05, deff = 1.5, resp_rate = 0.8,
#'                N = 10000))
#'
#' # A Wilson interval at the 90 percent level the size was planned for
#' summary(n_prop(p = 0.1, moe = 0.03, method = "wilson", alpha = 0.10))
#'
#' indicators <- data.frame(
#'   name = c("stunting", "vaccination", "anemia"),
#'   p    = c(0.30, 0.70, 0.10),
#'   moe  = c(0.05, 0.05, 0.03),
#'   deff = c(2.0, 1.5, 2.5)
#' )
#' summary(n_multi(indicators, resp_rate = 0.9))
#'
#' summary(n_cluster(stage_cost = c(500, 50), icc = 0.05, cv = 0.05,
#'                   eligibility_rate = 0.95, resp_rate = 0.9))
#'
#' # Several indicators under one budget, one of them short of its target
#' cluster_rows <- data.frame(
#'   name    = c("a", "b"),
#'   p       = c(0.35, 0.40),
#'   cv      = c(0.15, 0.03),
#'   icc_1 = 0.05
#' )
#' summary(n_cluster(cluster_rows, stage_cost = c(500, 50), budget = 1e5))
#'
#' @seealso [print.svyplan] for the compact display, and
#'   [summary.svyplan_alloc] for stratified allocations.
#' @export
summary.svyplan_n <- function(object, ...) {
  .check_unused_dots(...)
  if (object$type %in% c("proportion", "mean", "ratio")) {
    return(.summary_scalar_n(object))
  }
  if (identical(object$type, "multi")) {
    return(.summary_multi_n(object))
  }
  if (identical(object$type, "change")) {
    return(.summary_change_n(object))
  }
  if (identical(object$type, "pooled")) {
    return(.summary_pooled_n(object))
  }
  object
}

#' @rdname summary.svyplan_alloc
#' @export
summary.svyplan_alloc <- function(object, ...) {
  .check_unused_dots(...)
  .summary_classic_alloc(object, kind = "design")
}

#' @rdname summary.svyplan_alloc
#' @export
summary.svyplan_alloc_joint <- function(object, ...) {
  .check_unused_dots(...)
  .summary_bethel_alloc(object, kind = "design")
}

#' @rdname summary.svyplan_n
#' @export
summary.svyplan_cluster <- function(object, ...) {
  .check_unused_dots(...)
  if (!is.null(object$indicators)) {
    return(.summary_multi_cluster(object))
  }
  .summary_single_cluster(object)
}

#' Target a result was sized to, on the scale it was stated
#' @keywords internal
#' @noRd
.summary_target <- function(p) {
  for (scale in c("moe", "rmoe", "cv")) {
    if (!is.null(p[[scale]])) {
      return(list(scale = scale, value = p[[scale]]))
    }
  }
  NULL
}

#' @keywords internal
#' @noRd
.summary_estimand <- function(x) {
  p <- x$params
  switch(
    x$type,
    proportion = sprintf("p = %s", .fmt_prob(p[["p"]])),
    mean = if (!is.null(p[["mu"]])) {
      sprintf("a mean of %.4g", p[["mu"]])
    } else {
      sprintf("var = %.4g", p[["var"]])
    },
    ratio = sprintf("r = %.4g", p[["r"]])
  )
}

#' Precision of a single-indicator size at the size it fields
#' @keywords internal
#' @noRd
.field_prec_scalar <- function(x, n) {
  switch(
    x$type,
    proportion = prec_prop(x, n = n),
    mean = prec_mean(x, n = n),
    ratio = prec_ratio(x, n = n)
  )
}

#' Structured summary of a proportion, mean or ratio size
#' @keywords internal
#' @noRd
.summary_scalar_n <- function(object) {
  p <- object$params
  field <- .whole_count(object$n)
  fp <- .field_prec_scalar(object, field)
  resp <- (p[["resp_rate"]] %||% 1) * (p[["eligibility_rate"]] %||% 1)
  epu <- p[["eligible_per_unit"]]

  sizes <- c(continuous = object$n, field = field)
  if (resp < 1) {
    sizes[["respondents"]] <- field * resp
  }
  if (!is.null(epu)) {
    sizes[["eligible"]] <- field * resp * epu
  }
  sizes[["effective"]] <- field * resp * (epu %||% 1) / p[["deff"]]

  method <- object$method %||% "wald"
  interval <- NULL
  level <- 1 - p[["alpha"]]
  if (identical(object$type, "proportion") && method != "wald") {
    interval <- as.numeric(stats::confint(fp, level = level))
  }

  assumptions <- list(alpha = p[["alpha"]], deff = p[["deff"]],
                      resp_rate = p[["resp_rate"]] %||% 1)
  if (identical(object$type, "ratio")) {
    for (nm in c("cv_num", "cv_den", "component_cor", "unit_relvar")) {
      assumptions[[nm]] <- p[[nm]]
    }
  }
  assumptions$eligibility_rate <- p[["eligibility_rate"]]
  assumptions$eligible_per_unit <- epu
  if (is.finite(p[["N"]] %||% Inf)) {
    assumptions$N <- p[["N"]]
  }
  assumptions$df <- p[["df"]]

  structure(
    list(
      kind = "scalar",
      type = object$type,
      method = object$method,
      target = .summary_target(p),
      estimand = .summary_estimand(object),
      binding = object$binding,
      min_cases = p[["min_cases"]],
      sizes = sizes,
      precision = c(se = fp$se, moe = fp$moe, rmoe = fp$rmoe %||% NA_real_,
                    cv = fp$cv),
      expected_cases = fp$expected_cases,
      interval = interval,
      level = level,
      assumptions = assumptions
    ),
    class = c("summary.svyplan_n", "list")
  )
}

#' Structured summary of a change sample size
#' @keywords internal
#' @noRd
.summary_change_n <- function(object) {
  p <- object$params
  field <- .whole_count(object$n)
  fp <- prec_change(object, n = field)
  resp <- (p[["resp_rate"]] %||% 1) * (p[["eligibility_rate"]] %||% 1)
  deff <- p[["deff"]] %||% 1

  assumptions <- list(
    alpha = p[["alpha"]],
    deff = deff,
    resp_rate = p[["resp_rate"]] %||% 1,
    eligibility_rate = p[["eligibility_rate"]],
    overlap = p[["overlap"]],
    overlap_cor = p[["overlap_cor"]],
    N = if (is.finite(p[["N"]] %||% Inf)) p[["N"]],
    df = p[["df"]]
  )

  structure(
    list(
      kind = "change",
      type = "change",
      method = object$method,
      target = .summary_target(p),
      n = list(
        continuous = object$n,
        field = field,
        respondents = if (any(resp < 1)) field * resp else NULL,
        effective = field * resp / deff
      ),
      precision = c(se = fp$se, moe = fp$moe, rmoe = fp$rmoe %||% NA_real_,
                    cv = fp$cv),
      assumptions = assumptions
    ),
    class = c("summary.svyplan_n", "list")
  )
}

#' Structured summary of a pooled-occasions sample size
#' @keywords internal
#' @noRd
.summary_pooled_n <- function(object) {
  p <- object$params
  field <- .whole_count(object$n)
  fp <- prec_pooled(object, n = field)
  resp <- (p[["resp_rate"]] %||% 1) * (p[["eligibility_rate"]] %||% 1)
  deff <- p[["deff"]] %||% 1
  occ <- as.integer(p[["occasions"]])

  sizes <- c(continuous = object$n, field = field, total_field = field * occ)
  if (resp < 1) {
    sizes[["respondents"]] <- field * resp
  }
  sizes[["effective"]] <- field * resp / deff

  assumptions <- list(
    occasions = occ,
    alpha = p[["alpha"]],
    deff = deff,
    resp_rate = p[["resp_rate"]] %||% 1,
    eligibility_rate = p[["eligibility_rate"]],
    overlap = p[["overlap"]],
    overlap_cor = p[["overlap_cor"]],
    cor_decay = p[["cor_decay"]],
    N = if (is.finite(p[["N"]] %||% Inf)) p[["N"]],
    df = p[["df"]]
  )

  structure(
    list(
      kind = "pooled",
      type = "pooled",
      method = object$method,
      occasions = occ,
      target = .summary_target(p),
      sizes = sizes,
      overlap_line = .fmt_lag_overlap(p),
      precision = c(se = fp$se, moe = fp$moe, rmoe = fp$rmoe %||% NA_real_,
                    cv = fp$cv),
      assumptions = assumptions
    ),
    class = c("summary.svyplan_n", "list")
  )
}

#' Domain label of each row, from its domain columns
#'
#' One domain column gives its bare value. Several give `col=value` pairs
#' joined by `", "`, a value holding `,`, `=`, `"` or a backslash being
#' quoted with its quotes and backslashes escaped, so two domains never share
#' a label.
#' @keywords internal
#' @noRd
.summary_domain_label <- function(df, domain_cols) {
  if (length(domain_cols) == 0L) {
    return(rep(NA_character_, nrow(df)))
  }
  if (length(domain_cols) == 1L) {
    return(as.character(df[[domain_cols]]))
  }
  parts <- lapply(domain_cols, function(col) {
    v <- as.character(df[[col]])
    quote <- grepl("[,=\"\\\\]", v)
    v[quote] <- paste0("\"", gsub("([\"\\\\])", "\\\\\\1", v[quote]), "\"")
    paste0(col, "=", v)
  })
  do.call(paste, c(parts, sep = ", "))
}

#' Domain columns of a summary table, one row per domain
#' @keywords internal
#' @noRd
.summary_domain_values <- function(df, domain_cols) {
  out <- df[domain_cols]
  rownames(out) <- NULL
  out
}

#' Join a field precision table onto the per-row detail
#'
#' Rows are matched by name within domain, which the input checks make
#' present and unique. A row with only a respondent floor has no precision
#' row and joins to `NA`. Any other key that matches no row, or a repeated
#' precision key, is refused rather than joined to the wrong row.
#' @keywords internal
#' @noRd
.summary_join_prec <- function(detail, prec_detail, domain_cols) {
  key <- .summary_row_key(detail, domain_cols)
  prec_key <- .summary_row_key(prec_detail, domain_cols)
  idx <- match(key, prec_key)
  if (any(is.na(idx) & !is.na(detail$.target_scale)) ||
      anyDuplicated(prec_key) > 0L) {
    stop("internal error: the field precision rows do not match the fit's rows one to one",
         call. = FALSE)
  }
  prec_detail[idx, , drop = FALSE]
}

#' Key of a row's domain alone, or of its domain and name
#'
#' Built on `.domain_key()`, so domain values holding the separator of the
#' printed label cannot collide.
#' @keywords internal
#' @noRd
.summary_domain_key <- function(df, domain_cols) {
  if (length(domain_cols) == 0L) {
    return(rep("", nrow(df)))
  }
  .domain_key(df, domain_cols)
}

#' Key of a row's domain and name
#' @keywords internal
#' @noRd
.summary_row_key <- function(df, domain_cols) {
  paste(.summary_domain_key(df, domain_cols), as.character(df$name),
        sep = "\r")
}

#' Achieved precision on the scale each row's target was stated
#' @keywords internal
#' @noRd
.summary_achieved <- function(scale, prec) {
  out <- rep(NA_real_, length(scale))
  for (s in c("moe", "rmoe", "cv")) {
    hit <- !is.na(scale) & scale == s
    out[hit] <- prec[[paste0(".", s)]][hit]
  }
  out
}

#' Input rows labeled as the fit labels them
#'
#' A fit of unnamed indicator rows names its detail rows by their input
#' positions. The input table gets the same labels, so the summary can
#' match each detail row to its input row, also under interleaved domains.
#' @keywords internal
#' @noRd
.summary_label_rows <- function(object) {
  if (!"name" %in% names(object$indicators)) {
    object$indicators$name <- as.character(seq_len(nrow(object$indicators)))
    if ("name" %in% names(object$detail)) {
      object$detail$name <- as.character(object$detail$name)
    } else {
      object$detail$name <- as.character(seq_len(nrow(object$detail)))
    }
  }
  object
}

#' Structured summary of a multi-indicator size
#' @keywords internal
#' @noRd
.summary_multi_n <- function(object) {
  object <- .summary_label_rows(object)
  ind <- object$indicators
  detail <- object$detail
  dom_cols <- object$params$domain_cols %||% character(0)
  fp <- .prec_multi_from_fit(object, list(), field = TRUE)
  prec <- .summary_join_prec(detail, fp$detail, dom_cols)
  field_total <- .multi_field_n(object)

  dom <- object$domains
  row_field <- rep(field_total, nrow(detail))
  domains <- NULL
  domain_values <- NULL
  if (!is.null(dom)) {
    dom_field <- if (".share" %in% names(dom)) {
      field_total * dom$.share
    } else {
      .whole_count(dom$.n)
    }
    dom_label <- .summary_domain_label(dom, dom_cols)
    row_field <- dom_field[match(.summary_domain_key(detail, dom_cols),
                                 .summary_domain_key(dom, dom_cols))]
    domains <- data.frame(domain = dom_label, field = dom_field,
                          binding = dom$.binding)
    domain_values <- .summary_domain_values(dom, dom_cols)
  }

  ind_row <- match(.summary_row_key(detail, dom_cols),
                   .summary_row_key(ind, dom_cols))
  rate <- function(col) {
    v <- if (col %in% names(ind)) ind[[col]][ind_row] else NULL
    if (is.null(v)) rep(1, nrow(detail)) else ifelse(is.na(v), 1, v)
  }
  respondents <- row_field * rate("resp_rate") * rate("eligibility_rate") *
    rate("eligible_per_unit")
  deff <- rate("deff")

  indicators <- data.frame(
    domain = .summary_domain_label(detail, dom_cols),
    name = as.character(detail$name),
    target_scale = detail$.target_scale,
    target = detail$.target,
    achieved = .summary_achieved(detail$.target_scale, prec),
    cv = prec$.cv,
    deff = deff,
    field = row_field,
    respondents = respondents,
    needs = detail$.n,
    binding = detail$.binding,
    constraint = detail$.constraint %||% rep("precision", nrow(detail))
  )

  sizes <- c(continuous = object$n, field = field_total)
  if (length(unique(signif(respondents, 10))) == 1L &&
      !isTRUE(all.equal(respondents[1L], field_total))) {
    sizes[["respondents"]] <- respondents[1L]
  }

  alpha <- unique(ind$alpha)
  assumptions <- list(
    alpha = if (length(alpha) == 1L) alpha else NULL,
    deff = if (length(unique(deff)) == 1L) deff[1L] else NULL,
    resp_rate = if (length(unique(rate("resp_rate"))) == 1L) {
      rate("resp_rate")[1L]
    },
    eligibility_rate = .summary_constant_rate(rate("eligibility_rate")),
    eligible_per_unit = .summary_constant_rate(rate("eligible_per_unit")),
    prop_method = if ("p" %in% names(ind) && any(!is.na(ind$p))) {
      object$params$prop_method
    },
    domain_sampling = if (!is.null(dom)) {
      object$params$domain_sampling %||% "separate"
    },
    min_n_domain = object$params$min_n_domain
  )

  structure(
    list(
      kind = "multi",
      sizes = sizes,
      binding = object$binding,
      binding_cause = object$binding_cause,
      n_domain_max = object$n_domain_max,
      domains = domains,
      domain_values = domain_values,
      indicators = indicators,
      assumptions = assumptions
    ),
    class = c("summary.svyplan_n", "list")
  )
}

#' Cost of each stage of a whole cluster design
#' @keywords internal
#' @noRd
.summary_stage_cost <- function(n, stage_cost, fixed_cost = 0) {
  units <- cumprod(n)
  out <- c(units * stage_cost)
  names(out) <- names(n) %||% paste0("stage_", seq_along(n))
  if (!is.null(fixed_cost) && fixed_cost > 0) {
    out[["Fixed"]] <- fixed_cost
  }
  out
}

#' Question a cluster size answers
#' @keywords internal
#' @noRd
.summary_cluster_question <- function(p, several = FALSE) {
  target <- .summary_target(p)
  base <- if (!is.null(p$budget)) {
    sprintf("the most precise design for a budget of %s",
            format(p$budget, big.mark = "", scientific = FALSE))
  } else if (several) {
    "the cheapest design that meets every indicator's target"
  } else {
    sprintf("the cheapest design that gives %s = %s", target$scale,
            format(target$value))
  }
  sn <- p[["stage_n"]]
  if (!is.null(sn)) {
    fixed <- sn[!is.na(sn)]
    if (length(fixed) > 0L) {
      fnms <- names(fixed) %||% paste0("stage_", which(!is.na(sn)))
      base <- paste0(base, ", with ", paste(sprintf("%s fixed at %s",
                                                    fnms, fixed),
                                            collapse = " and "))
    }
  }
  base
}

#' Field and continuous designs side by side
#' @keywords internal
#' @noRd
.summary_cluster_design <- function(object) {
  op <- object$operational
  data.frame(
    row = c(names(object$n), "total_n", "cv", "cost"),
    field = c(unname(op$n), op$total_n, op$cv, op$cost),
    continuous = c(unname(object$n), object$total_n, object$cv, object$cost)
  )
}

#' Structured summary of a single-indicator cluster size
#' @keywords internal
#' @noRd
.summary_single_cluster <- function(object) {
  p <- object$params
  op <- object$operational
  fp <- .prec_cluster_from_fit(object, list(), field = TRUE)
  rr <- unname(p[["resp_rate"]] %||% rep(1, object$stages))
  resp_1 <- if (length(rr) == object$stages) rr[1L] else 1
  resp_all <- prod(rr) * (p[["eligibility_rate"]] %||% 1)
  epu <- p[["eligible_per_unit"]]
  respondents <- op$total_n * resp_all * (epu %||% 1)
  df <- tryCatch(as.double(suppressWarnings(design_df(object))),
                 error = function(e) NULL)

  field <- c(respondents = respondents,
             net_take = respondents / (op$n[[1L]] * resp_1),
             deff = as.double(design_effect(fp)),
             effective = as.double(effective_n(fp)))
  if (!is.null(df)) {
    field[["df"]] <- df
  }

  assumptions <- c(as.list(p$icc), list(unit_relvar = p$unit_relvar),
                   as.list(p$var_ratio))
  assumptions$resp_rate <- if (length(rr) == object$stages &&
                               all(rr[-object$stages] == 1)) {
    rr[object$stages]
  } else {
    rr
  }
  assumptions$eligibility_rate <- p[["eligibility_rate"]]
  assumptions$eligible_per_unit <- epu
  if (identical(p[["fpc"]], "stage")) {
    N <- p[["N"]]
    if (length(N) == object$stages) {
      names(N) <- .stage_names(object$stages, names(p$stage_cost))
    }
    assumptions$N <- N
  }

  structure(
    list(
      kind = "single",
      stages = object$stages,
      question = .summary_cluster_question(p),
      design = .summary_cluster_design(object),
      field = field,
      eligible = !is.null(epu),
      cost = .summary_stage_cost(op$n, p$stage_cost, p$fixed_cost),
      assumptions = assumptions,
      search = op$search,
      mode = .cluster_mode(p)
    ),
    class = c("summary.svyplan_cluster", "list")
  )
}

#' Structured summary of a multi-indicator cluster size
#' @keywords internal
#' @noRd
.summary_multi_cluster <- function(object) {
  object <- .summary_label_rows(object)
  p <- object$params
  op <- object$operational
  dom_cols <- p$domain_cols %||% character(0)
  detail <- object$detail
  fp <- .prec_cluster_from_indicator_fit(object, list(), field = TRUE)
  prec <- .summary_join_prec(detail, fp$detail, dom_cols)

  indicators <- data.frame(
    domain = .summary_domain_label(detail, dom_cols),
    name = as.character(detail$name),
    target_scale = detail$.target_scale,
    target = detail$.target,
    achieved = .summary_achieved(detail$.target_scale, prec),
    cv = prec$.cv,
    cv_target = detail$.cv_target,
    needs = detail$.n,
    binding = if (length(dom_cols) > 0L) {
      field <- op$domains$.binding[match(
        .summary_domain_label(detail, dom_cols),
        .summary_domain_label(op$domains, dom_cols)
      )]
      as.character(detail$name) == as.character(field)
    } else {
      as.character(detail$name) == as.character(op$binding)
    }
  )

  indicators <- cbind(indicators,
                      .summary_cluster_row_assumptions(object, detail, dom_cols))

  domains <- NULL
  domain_values <- NULL
  design <- NULL
  cost <- NULL
  if (length(dom_cols) > 0L) {
    od <- op$domains
    stage_cols <- names(object$n)
    domain_values <- .summary_domain_values(od, dom_cols)
    domains <- data.frame(
      domain = .summary_domain_label(od, dom_cols),
      od[stage_cols],
      total_n = od$.total_n,
      cv = od$.cv,
      cost = od$.cost,
      binding = od$.binding,
      check.names = FALSE
    )
  } else {
    design <- .summary_cluster_design(object)
    cost <- .summary_stage_cost(op$n, p$stage_cost, p$fixed_cost)
  }
  worst <- which.max(indicators$cv)
  binding <- op$binding
  binding_cause <- object$binding_cause
  if (length(dom_cols) > 0L) {
    field_worst <- which.max(op$domains$.cv)
    binding <- op$domains$.binding[field_worst]
    floored <- .cluster_domain_floored(object$domains$.total_n,
                                       p[["min_n_domain"]], p[["joint"]])
    binding_cause <- if (floored[field_worst]) "min_n_domain" else "precision"
  }

  structure(
    list(
      kind = "multi",
      stages = object$stages,
      question = .summary_cluster_question(p, several = is.null(p$budget)),
      design = design,
      domains = domains,
      domain_values = domain_values,
      total = c(total_n = op$total_n, cost = op$cost),
      worst = list(name = indicators$name[worst],
                   domain = indicators$domain[worst],
                   cv = indicators$cv[worst]),
      binding = binding,
      binding_cause = binding_cause,
      indicators = indicators,
      cost = cost,
      assumptions = c(
        .summary_constant_assumptions(indicators),
        list(
          eligibility_rate = p[["eligibility_rate"]],
          eligible_per_unit = p[["eligible_per_unit"]],
          prop_method = if ("p" %in% names(object$indicators) &&
                            any(!is.na(object$indicators$p))) {
            p[["prop_method"]]
          },
          stage_cost = p$stage_cost,
          fixed_cost = p$fixed_cost,
          allocation = if (isTRUE(p[["joint"]])) "joint",
          min_n_domain = p[["min_n_domain"]],
          fpc = if (identical(p[["fpc"]], "stage")) "stage"
        )
      ),
      search = object$operational$search,
      mode = .cluster_mode(p)
    ),
    class = c("summary.svyplan_cluster", "list")
  )
}

#' Identify domains held at their min_n_domain floor
#' @keywords internal
#' @noRd
.cluster_domain_floored <- function(total_n, min_n_domain, joint = FALSE) {
  if (is.null(min_n_domain)) {
    return(rep(FALSE, length(total_n)))
  }
  floors <- rep_len(as.numeric(min_n_domain), length(total_n))
  !is.na(floors) & abs(total_n - floors) <= 1e-6 * pmax(1, floors)
}

#' Per-row design quantities of a multi-indicator cluster fit
#'
#' The columns that apply to the design, in the order the single-indicator
#' summary lists them. `alpha` applies when a row states a `moe` or `rmoe`,
#' and `df` when a row carries one.
#' @keywords internal
#' @noRd
.summary_cluster_row_assumptions <- function(object, detail, dom_cols) {
  ind <- object$indicators
  n_b <- as.integer(object$stages) - 1L
  icc_cols <- paste0("icc_", seq_len(n_b))
  vr_cols <- paste0("var_ratio_", seq_len(n_b))
  rr_cols <- c(paste0("resp_rate_", seq_len(n_b)), "resp_rate")
  cols <- c(icc_cols, "unit_relvar", vr_cols, rr_cols)
  if (any(detail$.target_scale %in% c("moe", "rmoe"))) {
    cols <- c(cols, "alpha")
  }
  if ("df" %in% names(ind) && any(!is.na(ind$df))) {
    cols <- c(cols, "df")
  }
  row <- match(.summary_row_key(detail, dom_cols),
               .summary_row_key(ind, dom_cols))
  out <- lapply(cols, function(col) {
    if (col %in% c(rr_cols, vr_cols) && !col %in% names(ind)) {
      return(rep(1, nrow(detail)))
    }
    ind[[col]][row]
  })
  names(out) <- cols
  as.data.frame(out)
}

#' Names of the per-row design quantities a summary table carries
#' @keywords internal
#' @noRd
.summary_row_assumption_cols <- function(ind) {
  c(
    grep("^icc_[0-9]+$", names(ind), value = TRUE),
    intersect("unit_relvar", names(ind)),
    grep("^var_ratio_[0-9]+$", names(ind), value = TRUE),
    grep("^resp_rate_[0-9]+$", names(ind), value = TRUE),
    intersect(c("resp_rate", "alpha", "df"), names(ind))
  )
}

#' Per-row design quantities equal on every row
#' @keywords internal
#' @noRd
.summary_constant_assumptions <- function(ind) {
  out <- list()
  for (col in .summary_row_assumption_cols(ind)) {
    v <- ind[[col]]
    if (!.summary_is_constant(v)) {
      next
    }
    if (grepl("^resp_rate_[0-9]+$", col) && isTRUE(v[1L] == 1)) {
      next
    }
    out[[col]] <- v[1L]
  }
  out
}

#' @keywords internal
#' @noRd
.summary_is_constant <- function(v) {
  length(unique(signif(v, 12))) == 1L
}

## Printing

#' Padded label-value block under an optional title
#' @keywords internal
#' @noRd
.print_summary_block <- function(rows, title = NULL) {
  if (!is.null(title)) {
    cat(title, "\n", sep = "")
  }
  width <- max(nchar(names(rows), type = "width"))
  for (nm in names(rows)) {
    cat(sprintf("  %s%s  %s\n", nm,
                strrep(" ", width - nchar(nm, type = "width")), rows[[nm]]))
  }
  invisible(NULL)
}

#' Table with labels on the left and numbers on the right
#' @keywords internal
#' @noRd
.print_summary_table <- function(df, left = character(0)) {
  cols <- lapply(seq_along(df), function(j) {
    nm <- names(df)[j]
    v <- c(nm, as.character(df[[j]]))
    v[is.na(v)] <- ""
    pad <- max(nchar(v, type = "width")) - nchar(v, type = "width")
    if (nm %in% left) paste0(v, strrep(" ", pad)) else paste0(strrep(" ", pad), v)
  })
  lines <- do.call(paste, c(cols, sep = "  "))
  cat(paste0(" ", sub(" +$", "", lines)), sep = "\n")
  invisible(NULL)
}

#' Assumptions as label: value lines
#' @keywords internal
#' @noRd
.print_summary_assumptions <- function(a) {
  labels <- c(
    alpha = "alpha", deff = "design effect", resp_rate = "response rate",
    eligibility_rate = "eligibility rate",
    eligible_per_unit = "eligible per unit", N = "population size",
    df = "design df", prop_method = "proportion method",
    domain_sampling = "domain sampling", min_n_domain = "min_n_domain",
    stage_cost = "stage cost", fixed_cost = "fixed cost",
    allocation = "allocation", fpc = "finite population correction",
    occasions = "occasions", overlap = "overlap",
    overlap_cor = "overlap correlation", cor_decay = "correlation decay"
  )
  a <- a[!vapply(a, is.null, logical(1L))]
  if (length(a) == 0L) {
    return(invisible(NULL))
  }
  cat("\nAssumptions\n")
  for (nm in names(a)) {
    v <- a[[nm]]
    shown <- if (is.character(v)) {
      paste(v, collapse = ", ")
    } else if (nm %in% c("N", "stage_cost", "fixed_cost") ||
               grepl("^N_[0-9]+$", nm)) {
      shown_v <- format(v, big.mark = "", scientific = FALSE, trim = TRUE)
      if (identical(nm, "N") && !is.null(names(v))) {
        shown_v <- paste(names(v), "=", shown_v)
      }
      paste(shown_v, collapse = ", ")
    } else {
      paste(sprintf("%.4g", v), collapse = ", ")
    }
    lbl <- if (nm %in% names(labels)) {
      labels[[nm]]
    } else if (grepl("^resp_rate_([0-9]+)$", nm)) {
      sprintf("stage %s response rate", sub("^resp_rate_", "", nm))
    } else if (grepl("^N_([0-9]+)$", nm)) {
      sprintf("stage %s population size", sub("^N_", "", nm))
    } else {
      nm
    }
    cat(sprintf("  %s: %s\n", lbl, shown))
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.fmt_summary_size <- function(v) {
  if (isTRUE(all.equal(v, round(v)))) {
    format(round(v), scientific = FALSE)
  } else {
    sprintf("%.2f", v)
  }
}

#' @keywords internal
#' @noRd
.fmt_summary_target <- function(scale, value) {
  ifelse(is.na(scale), "", sprintf("%s %.4g", scale, value))
}

#' Shorten a label to keep a table inside 80 columns
#' @keywords internal
#' @noRd
.summary_label <- function(x, width = 20L) {
  vapply(as.character(x), .abbrev_bethel_id, character(1L), width = width,
         USE.NAMES = FALSE)
}

#' @rdname summary.svyplan_n
#' @export
print.summary.svyplan_n <- function(x, ...) {
  .check_unused_dots(...)
  if (identical(x$kind, "multi")) {
    .print_summary_multi_n(x)
  } else if (identical(x$kind, "change")) {
    .print_summary_change_n(x)
  } else if (identical(x$kind, "pooled")) {
    .print_summary_pooled_n(x)
  } else {
    .print_summary_scalar_n(x)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.print_summary_change_n <- function(x) {
  cat("Sample size for change between two occasions\n\n")
  if (!is.null(x$target)) {
    cat(sprintf("Question: per-occasion n for %s = %.4g\n",
                x$target$scale, x$target$value))
  }
  fmt_pair <- function(v) {
    if (length(v) == 2L) {
      paste(.fmt_summary_size(v[1L]), "then", .fmt_summary_size(v[2L]))
    } else {
      paste(.fmt_summary_size(v[1L]), "per occasion")
    }
  }
  rows <- c(
    Continuous = fmt_pair(x$n$continuous),
    Field = fmt_pair(x$n$field)
  )
  if (!is.null(x$n$respondents)) {
    rows[["Expected respondents"]] <- paste(sprintf("%.1f", x$n$respondents),
                                            collapse = " then ")
  }
  rows[["Effective (n / deff)"]] <- paste(sprintf("%.1f", x$n$effective),
                                          collapse = " then ")
  cat("\n")
  .print_summary_block(rows, "Sizes")

  pr <- x$precision
  labels <- c(se = "SE", moe = "MOE", rmoe = "Relative MOE", cv = "CV")
  prows <- character(0)
  for (nm in names(labels)) {
    if (!is.na(pr[[nm]])) {
      prows[[labels[[nm]]]] <- sprintf("%.4g", pr[[nm]])
    }
  }
  cat("\n")
  .print_summary_block(prows, "Precision at the field size")
  .print_summary_assumptions(x$assumptions)
}

#' @keywords internal
#' @noRd
.print_summary_pooled_n <- function(x) {
  cat(sprintf("Sample size for pooled estimate (%d occasions)\n\n", x$occasions))
  if (!is.null(x$target)) {
    cat(sprintf("Question: per-occasion n for %s = %.4g\n",
                x$target$scale, x$target$value))
  }
  s <- x$sizes
  rows <- c(
    "Continuous (per occasion)" = .fmt_summary_size(s[["continuous"]]),
    "Field (per occasion)" = .fmt_summary_size(s[["field"]]),
    "Total field (all occasions)" = .fmt_summary_size(s[["total_field"]])
  )
  if (!is.na(s["respondents"])) {
    rows[["Expected respondents"]] <- sprintf("%.1f per occasion", s[["respondents"]])
  }
  rows[["Effective (n / deff)"]] <- sprintf("%.1f per occasion", s[["effective"]])
  cat("\n")
  .print_summary_block(rows, "Sizes")
  if (!is.null(x$overlap_line) && nzchar(x$overlap_line)) {
    cat("\n", x$overlap_line, sep = "")
  }

  pr <- x$precision
  labels <- c(se = "SE", moe = "MOE", rmoe = "Relative MOE", cv = "CV")
  prows <- character(0)
  for (nm in names(labels)) {
    if (!is.na(pr[[nm]])) {
      prows[[labels[[nm]]]] <- sprintf("%.4g", pr[[nm]])
    }
  }
  cat("\n")
  .print_summary_block(prows, "Precision at the field size")
  .print_summary_assumptions(x$assumptions)
}

#' @keywords internal
#' @noRd
.print_summary_scalar_n <- function(x) {
  type_label <- x$type
  method_label <- if (!is.null(x$method)) paste0(" (", x$method, ")") else ""
  cat(sprintf("Sample size for %s%s\n\n", type_label, method_label))
  if (!is.null(x$target)) {
    cat(sprintf("Question: the n that gives %s = %.4g for %s\n",
                x$target$scale, x$target$value, x$estimand))
  }
  if (identical(x$binding, "min_cases")) {
    cat(sprintf("binding: min_cases (min_cases = %.4g)\n", x$min_cases))
  }

  s <- x$sizes
  rows <- c(Continuous = .fmt_summary_size(s[["continuous"]]),
            Field = .fmt_summary_size(s[["field"]]))
  if (!is.na(s["respondents"])) {
    rows[["Expected respondents"]] <- sprintf("%.1f", s[["respondents"]])
  }
  if (!is.na(s["eligible"])) {
    rows[["Eligible respondents"]] <- sprintf("%.1f", s[["eligible"]])
  }
  rows[["Effective (n / deff)"]] <- sprintf("%.1f", s[["effective"]])
  cat("\n")
  .print_summary_block(rows, "Sizes")

  pr <- x$precision
  labels <- c(se = "SE", moe = "MOE", rmoe = "Relative MOE", cv = "CV")
  rows <- character(0)
  for (nm in names(labels)) {
    if (!is.na(pr[[nm]])) {
      rows[[labels[[nm]]]] <- sprintf("%.4g", pr[[nm]])
    }
  }
  if (!is.null(x$expected_cases)) {
    rows[["Expected cases"]] <- sprintf("%.1f", x$expected_cases)
  }
  if (!is.null(x$interval)) {
    rows[["Interval"]] <- sprintf("%.4g to %.4g (%s%%)", x$interval[1L],
                                  x$interval[2L], format(100 * x$level))
  }
  cat("\n")
  .print_summary_block(rows, "Precision at the field size")
  .print_summary_assumptions(x$assumptions)
}

#' @keywords internal
#' @noRd
.print_summary_multi_n <- function(x) {
  cat("Multi-indicator sample size\n\n")
  cat("Question: the n that meets every indicator's target\n\n")
  s <- x$sizes
  rows <- c(Continuous = .fmt_summary_size(s[["continuous"]]),
            Field = .fmt_summary_size(s[["field"]]))
  if (!is.na(s["respondents"])) {
    rows[["Expected respondents"]] <- sprintf("%.1f", s[["respondents"]])
  }
  rows[["Binding"]] <- sprintf("%s (%s)", .summary_label(x$binding),
                               x$binding_cause %||% "precision")
  .print_summary_block(rows, "Sizes")

  ind <- x$indicators
  tab <- data.frame(
    Indicator = .summary_label(ind$name),
    Target = .fmt_summary_target(ind$target_scale, ind$target),
    Achieved = ifelse(is.na(ind$achieved), "", sprintf("%.4g", ind$achieved)),
    check.names = FALSE
  )
  if (is.null(x$assumptions$deff)) {
    tab$Deff <- sprintf("%.2f", ind$deff)
  }
  if (is.na(s["respondents"]) &&
      !isTRUE(all.equal(ind$respondents, ind$field))) {
    tab$Resp. <- sprintf("%.1f", ind$respondents)
  }
  tab$Needs <- .whole_count(ind$needs)
  tab$Binding <- ifelse(ind$binding, "*", "")
  if (any(ind$constraint != "precision")) {
    tab$`Set by` <- c(precision = "", min_cases = "cases",
                      min_respondents = "floor")[ind$constraint]
  }

  if (is.null(x$domains)) {
    cat("\nIndicators at the field size\n\n")
    .print_summary_table(tab, left = c("Indicator", "Binding", "Set by"))
  } else {
    d <- x$domains
    cat("\nDomains at the field size\n\n")
    .print_summary_table(
      cbind(.summary_domain_table(x$domain_values),
            data.frame(Field = vapply(d$field, .fmt_summary_size, ""),
                       Binding = .summary_label(d$binding))),
      left = c(names(x$domain_values), "Binding")
    )
    .print_summary_by_domain(tab, ind$domain, d$domain,
                             left = c("Indicator", "Binding", "Set by"))
  }
  .print_summary_assumptions(x$assumptions)
}

#' @rdname summary.svyplan_n
#' @export
print.summary.svyplan_cluster <- function(x, ...) {
  .check_unused_dots(...)
  heading <- sub(" ?[Aa]llocation$", "",
                 .cluster_heading(x$search, multi = identical(x$kind, "multi")))
  cat(trimws(sprintf("%s %d-stage allocation", heading, x$stages)), "\n\n",
      sep = "")
  cat(paste(strwrap(sprintf("Question: %s", x$question), width = 80L,
                    exdent = 2L), collapse = "\n"), "\n", sep = "")
  .print_stage_search(x$search, mode = x$mode)

  if (!is.null(x$design)) {
    d <- x$design
    fmt <- function(v, row) {
      ifelse(row == "cv", sprintf("%.4f", v),
             ifelse(row == "cost", sprintf("%.0f", v),
                    vapply(v, .fmt_summary_size, "")))
    }
    labels <- d$row
    labels[labels == "total_n"] <- "Total n"
    labels[labels == "cv"] <- if (identical(x$kind, "multi")) "Worst CV" else
      "CV"
    labels[labels == "cost"] <- "Cost"
    cat("\n")
    .print_summary_table(
      data.frame(` ` = labels, Field = fmt(d$field, d$row),
                 Continuous = fmt(d$continuous, d$row), check.names = FALSE),
      left = " "
    )
  }

  if (identical(x$kind, "single")) {
    f <- x$field
    rows <- c(
      sprintf("%.1f", f[["respondents"]]),
      sprintf("%.2f", f[["net_take"]]),
      sprintf("%.3f", f[["deff"]]),
      sprintf("%.1f", f[["effective"]])
    )
    names(rows) <- c(
      if (isTRUE(x$eligible)) "Eligible respondents" else
        "Expected respondents",
      "Net take per PSU", "Design effect (no FPC)", "Effective n (no FPC)"
    )
    if (!is.na(f["df"])) {
      rows[["Design df"]] <- format(f[["df"]])
    }
    cat("\n")
    .print_summary_block(rows, "Field design")
  } else {
    .print_summary_multi_cluster(x)
  }

  if (!is.null(x$cost)) {
    cat("\n")
    .print_summary_block(
      stats::setNames(sprintf("%.0f", x$cost), names(x$cost)),
      "Cost by stage"
    )
  }
  .print_summary_assumptions(x$assumptions)
  if (identical(x$kind, "multi")) {
    .print_summary_row_assumptions(x)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.print_summary_multi_cluster <- function(x) {
  ind <- x$indicators
  worst_label <- .summary_label(x$worst$name)
  if (!is.na(x$worst$domain)) {
    worst_label <- sprintf("%s in %s", worst_label,
                           .summary_label(x$worst$domain))
  }
  cat(sprintf("\nworst cv: %.4f (%s)\n", x$worst$cv, worst_label))
  cat(sprintf("binding: %s (%s)\n", .summary_label(x$binding),
              x$binding_cause %||% "precision"))
  short <- !is.na(ind$cv) & ind$cv > ind$cv_target * (1 + 1e-9)
  if (any(short)) {
    who <- .summary_label(ind$name[short])
    if (any(!is.na(ind$domain))) {
      who <- sprintf("%s in %s", who, .summary_label(ind$domain[short]))
    }
    cat(sprintf("short of target: %s\n", paste(who, collapse = ", ")))
  }

  if (!is.null(x$domains)) {
    d <- x$domains
    stage_cols <- setdiff(names(d), c("domain", "total_n", "cv", "cost",
                                      "binding"))
    tab <- d[stage_cols]
    tab$`Total n` <- d$total_n
    tab$CV <- sprintf("%.4f", d$cv)
    tab$Cost <- sprintf("%.0f", d$cost)
    tab$Binding <- .summary_label(d$binding)
    cat("\nDomains at the field design\n\n")
    .print_summary_table(cbind(.summary_domain_table(x$domain_values), tab),
                         left = c(names(x$domain_values), "Binding"))
    cat(sprintf("\n  Total n  %s\n  Cost     %.0f\n",
                format(x$total[["total_n"]]), x$total[["cost"]]))
  }

  tab <- data.frame(
    Indicator = .summary_label(ind$name),
    Target = .fmt_summary_target(ind$target_scale, ind$target),
    Achieved = ifelse(is.na(ind$achieved), "", sprintf("%.4g", ind$achieved)),
    CV = sprintf("%.4f", ind$cv),
    Needs = .whole_count(ind$needs),
    Binding = ifelse(ind$binding, "*", ""),
    check.names = FALSE
  )
  if (all(ind$target_scale == "cv", na.rm = TRUE)) {
    tab$CV <- NULL
  }
  if (is.null(x$domains)) {
    cat("\nIndicators at the field design\n\n")
    .print_summary_table(tab, left = c("Indicator", "Binding"))
  } else {
    .print_summary_by_domain(tab, ind$domain, x$domains$domain,
                             left = c("Indicator", "Binding"))
  }
}

#' Per-row design quantities that vary, as a table wrapped at 80 columns
#' @keywords internal
#' @noRd
.print_summary_row_assumptions <- function(x) {
  ind <- x$indicators
  cols <- .summary_row_assumption_cols(ind)
  cols <- cols[!vapply(ind[cols], .summary_is_constant, logical(1L))]
  if (length(cols) == 0L) {
    return(invisible(NULL))
  }
  key <- data.frame(Indicator = .summary_label(ind$name), check.names = FALSE)
  if (!is.null(x$domains)) {
    key <- cbind(data.frame(Domain = .summary_label(ind$domain)), key)
  }
  vals <- as.data.frame(lapply(ind[cols], function(v) {
    ifelse(is.na(v), "", sprintf("%.4g", v))
  }), check.names = FALSE)
  cat("\nAssumptions by indicator\n")
  .print_summary_table_wrapped(key, vals, left = names(key))
}

#' A table split into blocks of value columns, each within 80 columns
#'
#' Every block repeats the key columns.
#' @keywords internal
#' @noRd
.print_summary_table_wrapped <- function(key, vals, left = character(0),
                                         width = 80L) {
  col_width <- function(df) {
    vapply(seq_along(df), function(j) {
      v <- c(names(df)[j], as.character(df[[j]]))
      v[is.na(v)] <- ""
      max(nchar(v, type = "width"))
    }, numeric(1L))
  }
  used <- 1 + sum(col_width(key)) + 2 * (ncol(key) - 1L)
  w <- col_width(vals) + 2
  block <- integer(0)
  start <- 1L
  for (j in seq_along(w)) {
    if (j > start && used + sum(w[start:j]) > width) {
      block <- c(block, j - 1L)
      start <- j
    }
  }
  ends <- c(block, length(w))
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  for (b in seq_along(ends)) {
    cat("\n")
    .print_summary_table(cbind(key, vals[starts[b]:ends[b]]), left = left)
  }
  invisible(NULL)
}

#' Domain columns of a printed domain table, each cut at 20 characters
#' @keywords internal
#' @noRd
.summary_domain_table <- function(values) {
  out <- as.data.frame(lapply(values, .summary_label), check.names = FALSE)
  names(out) <- names(values)
  out
}

#' One table per domain, in the order of the domain table
#'
#' Rows are grouped by the domain label, which `.summary_domain_label()`
#' makes distinct for distinct domains.
#' @keywords internal
#' @noRd
.print_summary_by_domain <- function(tab, row_domain, domains, left) {
  for (lev in domains) {
    cat(sprintf("\nIndicators in %s\n\n", .summary_label(lev, width = 60L)))
    .print_summary_table(tab[row_domain == lev, , drop = FALSE], left = left)
  }
  invisible(NULL)
}

#' A rate shared by every row, shown only when it differs from 1
#' @keywords internal
#' @noRd
.summary_constant_rate <- function(v) {
  u <- unique(v)
  if (length(u) == 1L && !isTRUE(all.equal(u, 1))) u else NULL
}
