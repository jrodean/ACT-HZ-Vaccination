# ==============================================================================
# 18_figS39_S40_age_groups.R
#
# Figure S39 (clinical age-at-vaccination bands: <75, 75-79, 80-84, 85+) and
# Figure S40 (quartiles of age at vaccination, group medians ~69/76/83/89).
# Nonparametric companion to Figure S38's continuous quadratic model: "a
# single model in which all
# event-time indicators are interacted with indicators for age at ZVL
# vaccination group" -- i.e. event_time x group fully interacted (no shared
# main effect), group indicators absorbed by subject FE since time-invariant.
# Reported at 2, 4, 6, 8 years after vaccination, same as Figure S38, as
# points with 95% CIs faceted by years-after-vaccination.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
baseline_date <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, baseline_date = visit_date, b_visit_age)
age_vax <- panel %>% distinct(subject, .keep_all = TRUE) %>% select(subject, zost1_date) %>%
  left_join(baseline_date, by = "subject") %>%
  mutate(age_at_zost1 = b_visit_age + as.numeric(zost1_date - baseline_date) / 365.25) %>%
  filter(!is.na(age_at_zost1)) %>% select(subject, age_at_zost1)

primary0 <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date") %>%
  left_join(age_vax, by = "subject")

taus <- c(trt_p1 = 2, trt_p2 = 4, trt_p3 = 6, trt_p4 = 8)

run_group_model <- function(primary, group_labels, fig_label, outname, panels_shown = c(2, 4, 6, 8)) {
  # safe level names (display labels can contain characters like "<" that
  # break formula parsing); mapped back to display labels for plotting.
  # Explicitly construct each term x group dummy as its own 0/1 column
  # (rather than relying on formula-level factor expansion, which fixest
  # does not handle the same way base R model.matrix does) -- fully
  # interacted, no shared main effect, matching "all event-time indicators
  # are interacted with indicators for age at ZVL vaccination group".
  safe_levels <- paste0("grp", 1:4)
  d <- primary %>% filter(!is.na(age_group))

  int_terms <- character(0)
  for (g in 1:4) {
    for (term in es_regs) {
      colname <- paste0(term, "_", safe_levels[g])
      d[[colname]] <- d[[term]] * as.integer(d$age_group == g)
      int_terms <- c(int_terms, colname)
    }
  }
  fml <- as.formula(paste0("casi_total ~ ", paste(int_terms, collapse = " + "), " | subject + visit_year + visit_age"))
  m <- fixest::feols(fml, data = d, vcov = "hetero")
  cat(sprintf("%s: N rows=%d, N subjects=%d\n", outname, nobs(m), n_distinct(d$subject)))

  ct <- as.data.frame(summary(m)$coeftable); ct$term <- rownames(ct)
  ci <- as.data.frame(confint(m)); ci$term <- rownames(ci)
  res <- merge(ct, ci, by = "term")
  names(res) <- c("term", "estimate", "se", "tstat", "pvalue", "ci_low", "ci_high")

  label_map <- setNames(group_labels, safe_levels)
  rows <- list()
  for (tau_term in names(taus)) {
    for (g in safe_levels) {
      key <- paste0(tau_term, "_", g)
      hit <- res[res$term == key, ]
      if (nrow(hit) == 0) next
      rows[[length(rows) + 1]] <- tibble::tibble(
        years_after = taus[tau_term], age_group = unname(label_map[g]),
        estimate = hit$estimate, ci_low = hit$ci_low, ci_high = hit$ci_high
      )
    }
  }
  pred <- bind_rows(rows) %>%
    filter(years_after %in% panels_shown) %>%
    mutate(panel_label = factor(paste(years_after, "yrs after ZVL"), levels = paste(panels_shown, "yrs after ZVL")),
           # age_group ascending bottom-to-top on the y-axis (<75/Q1 at bottom,
           # 85+/Q4 at top) -- matches the published figure; ggplot's discrete
           # y-axis plots the FIRST factor level at the bottom, so levels are
           # kept in ascending order (not reversed)
           age_group = factor(age_group, levels = group_labels))

  dir.create("output/tables", showWarnings = FALSE); dir.create("output/figures", showWarnings = FALSE)
  readr::write_csv(pred, sprintf("output/tables/%s_data.csv", outname))

  # Benefit/Harm arrow annotation, placed above x=0 within each facet panel,
  # matching the house style used elsewhere in the paper's single-series plots
  xrange <- range(c(pred$ci_low, pred$ci_high), na.rm = TRUE)
  xpad <- diff(xrange) * 0.12
  arrow_y <- length(group_labels) + 0.65
  harm_x <- -xpad * 1.6; benefit_x <- xpad * 1.6

  p <- ggplot(pred, aes(x = estimate, y = age_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, color = "black", linewidth = 0.5) +
    geom_point(size = 2.2, color = "black") +
    annotate("segment", x = 0, xend = harm_x, y = arrow_y, yend = arrow_y,
             arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = "grey30", linewidth = 0.4) +
    annotate("segment", x = 0, xend = benefit_x, y = arrow_y, yend = arrow_y,
             arrow = arrow(length = unit(0.1, "cm"), type = "closed"), color = "grey30", linewidth = 0.4) +
    annotate("text", x = harm_x, y = arrow_y, label = "Harm", hjust = 1.05, vjust = 0.5, fontface = "italic", size = 4.2) +
    annotate("text", x = benefit_x, y = arrow_y, label = "Benefit", hjust = -0.05, vjust = 0.5, fontface = "italic", size = 4.2) +
    coord_cartesian(ylim = c(0.5, length(group_labels) + 1.1), clip = "off") +
    facet_wrap(~panel_label, nrow = 1) +
    labs(x = "Effect on CASI Total Score", y = "Age at ZVL vaccination (years)") +
    theme_bw(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black", size = 16),
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 13, color = "black"),
      plot.margin = margin(18, 8, 8, 8)
    )
  ggsave(sprintf("output/figures/%s.png", outname), p, width = 3.4 * length(panels_shown), height = 5, dpi = 600)
  cat(sprintf("saved output/figures/%s.png\n", outname))
  invisible(pred)
}

cat("################ Figure S39: clinical age bands (<75, 75-79, 80-84, 85+) ################\n")
primary_clin <- primary0 %>%
  mutate(age_group = case_when(
    age_at_zost1 < 75 ~ 1, age_at_zost1 < 80 ~ 2, age_at_zost1 < 85 ~ 3, TRUE ~ 4
  ))
run_group_model(primary_clin, c("<75", "75-79", "80-84", "85+"),
                "Figure S39", "FigureS39_age_bands", panels_shown = c(2, 4, 6, 8))

cat("\n################ Figure S40: quartiles of age at vaccination ################\n")
q <- quantile((age_vax %>% distinct(subject, .keep_all = TRUE))$age_at_zost1, c(.25, .5, .75))
primary_quart <- primary0 %>%
  mutate(age_group = case_when(
    age_at_zost1 < q[1] ~ 1, age_at_zost1 < q[2] ~ 2, age_at_zost1 < q[3] ~ 3, TRUE ~ 4
  ))
qlabs <- c("Q1 (median age 69)", "Q2 (median age 76)", "Q3 (median age 83)", "Q4 (median age 89)")
# Figure S40 caption explicitly says "four, six, and eight years after
# vaccination" -- no 2-year panel, unlike S39
run_group_model(primary_quart, qlabs, "Figure S40", "FigureS40_age_quartiles", panels_shown = c(2, 4, 6, 8))
