# ==============================================================================
# 01_event_study_utils.R
#
# Shared helpers for every event-study analysis in this package:
#   - build_event_time(): adds event_time + trt_m1..trt_m6 / trt_p0..trt_p6
#   - run_event_study(): fits the primary feols specification and returns
#     tidy coefficients + the pre-trend joint Wald test
#
# Primary specification (Methods / Supplement "Regression Specifications"):
#   Y_it = subject FE + visit_year FE + visit_age FE + sum_tau delta_tau 1{event_time=tau} + e_it
#   reference period: tau = -2 (trt_m1, forced to 0 / omitted)
#   earliest pre-period (tau <= -14) also omitted (collapsed into the fixed effects)
#   event_time floored to nearest 2-year bin, capped to [-14, 12]
#   vcov: heteroskedasticity-robust ("hetero")
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(ggplot2)
  library(patchwork)
})

build_event_time <- function(data, vax_date_col) {
  vax_date <- data[[vax_date_col]]
  et <- floor(as.numeric(data$visit_date - vax_date) / (2 * 365)) * 2
  et <- pmin(et, 12)
  et <- pmax(et, -14)
  data$event_time <- et

  for (i in 1:6) data[[paste0("trt_m", i)]] <- as.integer(data$event_time == (-i * 2))
  for (i in 0:6) data[[paste0("trt_p", i)]] <- as.integer(data$event_time == (i * 2))
  data$trt_m1 <- 0L   # omitted reference category (t = -2)

  data
}

es_regs <- c(paste0("trt_m", 6:2), paste0("trt_p", 0:6))

coef_map_main <- c(
  trt_m6 = "12 years before", trt_m5 = "10 years before", trt_m4 = "8 years before",
  trt_m3 = "6 years before",  trt_m2 = "4 years before",
  trt_p0 = "ZVL vaccination year", trt_p1 = "2 years after", trt_p2 = "4 years after",
  trt_p3 = "6 years after", trt_p4 = "8 years after", trt_p5 = "10 years after",
  trt_p6 = "12 years after"
)

# Abbreviated axis-tick labels, matching the published figure style.
coef_map_short <- c(
  trt_m6 = "12 yrs before", trt_m5 = "10 yrs before", trt_m4 = "8 yrs before",
  trt_m3 = "6 yrs before",  trt_m2 = "4 yrs before", trt_m1 = "2 yrs before",
  trt_p0 = "ZVL vaccination year", trt_p1 = "2 yrs after", trt_p2 = "4 yrs after",
  trt_p3 = "6 yrs after", trt_p4 = "8 yrs after", trt_p5 = "10 yrs after",
  trt_p6 = "12 yrs after"
)

# "Benefit/Harm" arrow annotation at the left edge, matching the published
# style for single-series event-study plots. harm_coded = TRUE for outcomes
# where a higher value is worse (e.g. P(dementia), P(CASI<86)), which flips
# the label positions (Harm on top, Benefit on bottom) -- confirmed against
# the real Figure S26 image, which flips exactly this way for its three
# harm-coded dementia-probability outcomes.
benefit_harm_annotation <- function(d, harm_coded = FALSE) {
  rng <- range(c(d$ci_low, d$ci_high), na.rm = TRUE)
  x0 <- min(d$x) - 0.45
  top_label <- if (harm_coded) "Harm" else "Benefit"
  bottom_label <- if (harm_coded) "Benefit" else "Harm"
  list(
    ggplot2::annotate("segment", x = x0, xend = x0, y = 0, yend = rng[2] * 0.35,
                       arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")), color = "grey40"),
    ggplot2::annotate("segment", x = x0, xend = x0, y = 0, yend = rng[1] * 0.35,
                       arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")), color = "grey40"),
    ggplot2::annotate("text", x = x0, y = rng[2] * 0.5, label = top_label, color = "grey40", size = 3.2),
    ggplot2::annotate("text", x = x0, y = rng[1] * 0.5, label = bottom_label, color = "grey40", size = 3.2)
  )
}

pretrend_wald_p <- function(mod, terms_pre = c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2")) {
  terms_pre <- intersect(terms_pre, names(stats::coef(mod)))
  if (length(terms_pre) < 2) return(NA_real_)
  b <- stats::coef(mod)[terms_pre]
  V <- stats::vcov(mod)[terms_pre, terms_pre, drop = FALSE]
  if (anyNA(b) || anyNA(V)) return(NA_real_)
  Vinv <- tryCatch(solve(V), error = function(e) MASS::ginv(V))
  W <- as.numeric(t(b) %*% Vinv %*% b)
  1 - stats::pchisq(W, df = length(b))
}

run_event_study <- function(data, outcome, fe = "subject + visit_year + visit_age",
                             regs = es_regs, vcov = "hetero") {
  fml <- as.formula(paste0(outcome, " ~ ", paste(regs, collapse = " + "), " | ", fe))
  mod <- fixest::feols(fml, data = data, vcov = vcov)
  list(mod = mod, pretrend_p = pretrend_wald_p(mod))
}

# ------------------------------------------------------------------------
# Difference-of-coefficients table, used for every "Difference in
# event-study treatment effects: by X" Table S in the supplement, and for
# every Figure 4 / subgroup HTE panel.
#
# Two separate models are used for two different purposes:
#
#  (1) Displayed Estimate/CI/P/BH columns: from two separate-sample
#      regressions (group==1, group==0), each with its own fixed-effect
#      absorption:
#        diff = b1 - b0,  se = sqrt(se1^2 + se0^2)
#
#  (2) Pre-trend joint Wald p-value: from a separate, pooled regression on
#      the combined (group 0 + group 1) sample, with the event-time dummies
#      as main effects and interacted with the group indicator:
#        casi_total ~ trt_m6...trt_p6 + group:trt_m6...group:trt_p6 | FE
#      followed by a joint Wald test on the five pre-period interaction
#      coefficients. This is not mathematically equivalent to combining the
#      two separate models' covariance matrices (the fixed-effect
#      absorption and residual-variance pooling differ), so the pooled
#      specification is estimated directly rather than derived from (1).
# ------------------------------------------------------------------------
diff_table <- function(data, split_var, label1 = "Group 1", label0 = "Group 0",
                        outcome = "casi_total", fe = "subject + visit_year + visit_age") {
  d1 <- data %>% dplyr::filter(.data[[split_var]] == 1)
  d0 <- data %>% dplyr::filter(.data[[split_var]] == 0)
  fml <- as.formula(paste0(outcome, " ~ ", paste(es_regs, collapse = " + "), " | ", fe))
  m1 <- fixest::feols(fml, data = d1, vcov = "hetero")
  m0 <- fixest::feols(fml, data = d0, vcov = "hetero")

  b1 <- coef(m1)[es_regs]; b0 <- coef(m0)[es_regs]
  se1 <- se(m1)[es_regs];  se0 <- se(m0)[es_regs]
  d <- b1 - b0
  se_d <- sqrt(se1^2 + se0^2)
  z <- d / se_d
  pval <- 2 * pnorm(-abs(z))
  ci_lo <- d - 1.96 * se_d
  ci_hi <- d + 1.96 * se_d
  bh <- p.adjust(pval, method = "BH")

  pre_terms <- c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2")
  pooled <- dplyr::bind_rows(d1, d0)
  inter_terms <- paste0(es_regs, ":", split_var)
  fml_pooled <- as.formula(paste0(
    outcome, " ~ ", paste(es_regs, collapse = " + "), " + ",
    paste(inter_terms, collapse = " + "), " | ", fe
  ))
  m_pooled <- fixest::feols(fml_pooled, data = pooled, vcov = "hetero")
  inter_pre_terms <- paste0(pre_terms, ":", split_var)
  b_int <- stats::coef(m_pooled)[inter_pre_terms]
  V_int <- stats::vcov(m_pooled)[inter_pre_terms, inter_pre_terms]
  Vinv <- tryCatch(solve(V_int), error = function(e) MASS::ginv(V_int))
  W <- as.numeric(t(b_int) %*% Vinv %*% b_int)
  pretrend_p <- 1 - pchisq(W, df = length(pre_terms))

  tbl <- tibble::tibble(term = es_regs, label = coef_map_main[es_regs],
                         estimate = d, ci_low = ci_lo, ci_high = ci_hi,
                         pvalue = pval, bh = bh)
  list(tbl = tbl, pretrend_p = pretrend_p, n1 = nobs(m1), n0 = nobs(m0),
       m1 = m1, m0 = m0, m_pooled = m_pooled, label1 = label1, label0 = label0)
}

run_panel <- function(data, split_var, label1, label0, name, outdir = "output/tables") {
  res <- diff_table(data, split_var, label1, label0)
  readr::write_csv(res$tbl, sprintf("%s/TableS_HTE_%s.csv", outdir, name))
  cat(sprintf("\n--- %s: %s (N=%d) vs %s (N=%d) ---\n", name, label1, res$n1, label0, res$n0))
  cat(sprintf("Pre-trend joint Wald p (on subgroup DIFFERENCE): %.4f\n", res$pretrend_p))
  print(res$tbl)
  res
}

# ------------------------------------------------------------------------
# Standard event-study coefficient plot (single series, with omitted t=-2
# reference point drawn in at zero).
# ------------------------------------------------------------------------
es_plot <- function(tbl, title = "", ytitle = "Effect on CASI Total Score", benefit_harm = TRUE,
                     xtitle = "Time relative to the year of ZVL vaccination", harm_coded = FALSE,
                     vax0_label = NULL, ylim = NULL) {
  ref_row <- tibble::tibble(term = "trt_m1", label = "2 years before", estimate = 0,
                             ci_low = 0, ci_high = 0, pvalue = NA, bh = NA)
  d <- dplyr::bind_rows(tbl, ref_row)
  ord <- c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2","trt_m1",
           "trt_p0","trt_p1","trt_p2","trt_p3","trt_p4","trt_p5","trt_p6")
  d$term <- factor(d$term, levels = ord)
  d <- d[order(d$term), ]
  d$x <- seq_len(nrow(d))
  tick_labels <- coef_map_short
  if (!is.null(vax0_label)) tick_labels["trt_p0"] <- vax0_label
  p <- ggplot(d, aes(x = x, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 6.5, linetype = "dashed", color = "red")
  if (benefit_harm) p <- p + benefit_harm_annotation(d, harm_coded = harm_coded)
  p +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = d$x, labels = tick_labels[as.character(d$term)]) +
    (if (!is.null(ylim)) scale_y_continuous(limits = ylim) else NULL) +
    labs(x = xtitle, y = ytitle, title = title) +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}

# Two-series (separate-sample HTE) plot, e.g. Figure 4 panels
es_plot_two <- function(mod1, mod0, label1, label0, title = "", ytitle = "Effect on CASI Total Score") {
  get_tbl <- function(mod, grp) {
    tbl <- es_table(mod)
    ref <- tibble::tibble(term = "trt_m1", label = "2 years before", estimate = 0,
                           ci_low = 0, ci_high = 0, pvalue = NA, bh = NA)
    dplyr::bind_rows(tbl, ref) %>% dplyr::mutate(group = grp)
  }
  d <- dplyr::bind_rows(get_tbl(mod1, label1), get_tbl(mod0, label0))
  ord <- c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2","trt_m1",
           "trt_p0","trt_p1","trt_p2","trt_p3","trt_p4","trt_p5","trt_p6")
  d$term <- factor(d$term, levels = ord)
  d <- d[order(d$group, d$term), ]
  xmap <- setNames(seq_along(ord), ord)
  d$x <- xmap[as.character(d$term)]
  dodge <- ifelse(d$group == label1, -0.12, 0.12)
  d$x <- d$x + dodge
  ggplot(d, aes(x = x, y = estimate, color = group, shape = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 6.5, linetype = "dashed", color = "red") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_along(ord),
                        labels = coef_map_short[ord] %>% { .[is.na(.)] <- "2 yrs before"; . }) +
    scale_color_manual(values = c("black", "grey50")) +
    labs(x = "Time relative to the year of ZVL vaccination", y = ytitle, title = title, color = NULL, shape = NULL) +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}

es_table <- function(mod) {
  s <- summary(mod)
  ct <- as.data.frame(s$coeftable)
  ct$term <- rownames(ct)
  ci <- as.data.frame(confint(mod))
  ci$term <- rownames(ci)
  out <- merge(ct, ci, by = "term")
  names(out) <- c("term", "estimate", "se", "tstat", "pvalue", "ci_low", "ci_high")
  out$bh <- p.adjust(out$pvalue, method = "BH")
  out$label <- coef_map_main[out$term]
  # match() returns NA for any es_regs term the model dropped (e.g. no
  # observations in that bin) -- drop those rather than keep an all-NA row.
  out <- out[match(es_regs, out$term), ]
  out <- out[!is.na(out$term), ]
  out[, c("term", "label", "estimate", "ci_low", "ci_high", "pvalue", "bh")]
}

# ------------------------------------------------------------------------
# Multi-panel grid of es_plot()s, matching how the paper displays a set of
# related outcomes as lettered (A), (B), ... sub-panels of one Figure.
# res_list: a list of lists, each with $tbl (an es_table()) and $label.
# ------------------------------------------------------------------------
panel_grid <- function(res_list, ytitle_fun = function(label) paste0("Effect on ", label), ncol = 3,
                        benefit_harm = TRUE, border = FALSE,
                        xtitle_fun = function(label) "Time relative to the year of ZVL vaccination",
                        harm_coded = FALSE, vax0_label_fun = function(label) NULL, ylim = NULL) {
  letters_ <- LETTERS[seq_along(res_list)]
  plots <- lapply(seq_along(res_list), function(i) {
    r <- res_list[[i]]
    p <- es_plot(r$tbl, title = sprintf("(%s) %s", letters_[i], r$label), ytitle = ytitle_fun(r$label),
                 benefit_harm = benefit_harm, xtitle = xtitle_fun(r$label), harm_coded = harm_coded,
                 vax0_label = vax0_label_fun(r$label), ylim = ylim) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
                      axis.title = ggplot2::element_text(size = 8), axis.text = ggplot2::element_text(size = 6))
    if (border) p <- p + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA))
    p
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

# ------------------------------------------------------------------------
# Illustrative ZVL-effect-vs-age-related-decline trajectory (Figure 2 at
# age 75, Figure S6 at age 80 -- same construction, different illustration
# age; see each figure's caption).
# ------------------------------------------------------------------------
build_illustrative_trajectory <- function(panel, vax_illustration_age) {
  primary <- panel %>% dplyr::filter(!is.na(zost1_date), !is.na(casi_total)) %>%
    build_event_time("zost1_date")
  es_res <- run_event_study(primary, "casi_total")
  es_tbl <- es_table(es_res$mod)

  age_fe <- fixest::fixef(es_res$mod)$visit_age
  age_fe <- age_fe[as.numeric(names(age_fe)) %>% dplyr::between(65, 90)]
  age_fe_centered <- age_fe - age_fe[["65"]]
  baseline_65 <- mean(primary$casi_total[round(primary$visit_age) == 65], na.rm = TRUE)
  age_curve <- tibble::tibble(age_bin = as.numeric(names(age_fe_centered)),
                               mean_casi = baseline_65 + unname(age_fe_centered))

  post <- es_tbl %>% dplyr::filter(grepl("^trt_p", term)) %>%
    dplyr::mutate(years_after = as.numeric(sub("trt_p", "", term)) * 2,
                  age = vax_illustration_age + years_after)
  age_fun <- approxfun(age_curve$age_bin, age_curve$mean_casi, rule = 2)
  post <- post %>% dplyr::mutate(age_curve_value = age_fun(age),
                                  illustrative_casi = age_curve_value + estimate,
                                  illustrative_lo = age_curve_value + ci_low,
                                  illustrative_hi = age_curve_value + ci_high)

  without_zvl <- age_curve %>% dplyr::transmute(age = age_bin, casi = mean_casi, series = "Trajectory without ZVL")
  with_zvl <- dplyr::bind_rows(
    age_curve %>% dplyr::filter(age_bin <= vax_illustration_age) %>% dplyr::transmute(age = age_bin, casi = mean_casi),
    post %>% dplyr::transmute(age, casi = illustrative_casi)
  ) %>% dplyr::mutate(series = "Trajectory with ZVL")
  ribbon_df <- post %>% dplyr::transmute(age, ymin = illustrative_lo, ymax = illustrative_hi)
  traj <- dplyr::bind_rows(without_zvl, with_zvl) %>%
    dplyr::mutate(series = factor(series, levels = c("Trajectory with ZVL", "Trajectory without ZVL")))

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = ribbon_df, ggplot2::aes(x = age, ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.3) +
    ggplot2::geom_line(data = traj, ggplot2::aes(x = age, y = casi, linetype = series), color = "black", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = vax_illustration_age, linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = vax_illustration_age + 0.5, y = max(traj$casi), label = "ZVL vaccination",
                       hjust = 0, vjust = 1, size = 3.5) +
    ggplot2::scale_linetype_manual(values = c("Trajectory with ZVL" = "solid", "Trajectory without ZVL" = "dotted")) +
    ggplot2::labs(x = "Age (years)", y = "CASI Total Score", linetype = NULL, title = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  list(plot = p, data = dplyr::bind_rows(without_zvl, with_zvl))
}
