# ==============================================================================
# 04_robustness_checks.R
#
#   - Extended to all observed event-time periods (no +/-14/12 capping)
#   - Negative control: other vaccines (pneumococcal, Td) on CASI Total Score
#   - Negative control outcomes: top-5 chronic conditions (CCW), top-5
#     hospitalization categories (IP), using the primary ZVL event-study spec
#   - Falsification test: pseudo-vaccination dates (pre-2006, ZVL-eventual
#     recipients only)
#   - Attrition check: sample restricted to participants observed from -6 to
#     +6 years relative to vaccination
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# simple plot for broom::tidy()-style tables (term/estimate/conf.low/conf.high)
tidy_plot <- function(tbl, title, xlabs = NULL, ytitle = "Effect on CASI Total Score") {
  d <- tbl %>% mutate(x = row_number())
  if (!is.null(xlabs)) d$lab <- xlabs else d$lab <- d$term
  ggplot(d, aes(x = x, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, color = "grey40") +
    geom_point(size = 2) +
    scale_x_continuous(breaks = d$x, labels = d$lab) +
    labs(x = "Time relative to ZVL vaccination", y = ytitle, title = title) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

primary0 <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total))

# ------------------------------------------------------------------------
# (1) Figure S4: extended to all observed periods. Matches the original
# Stata specification (analysis_Xi_1225.do, "Robustness: full periods"):
# event time capped at -20/+14 years (not the primary's -14/+12), standard
# errors clustered by vaccination year (not heteroskedasticity-robust).
# ------------------------------------------------------------------------
build_event_time_extended <- function(data, vax_date_col) {
  vax_date <- data[[vax_date_col]]
  et <- floor(as.numeric(data$visit_date - vax_date) / (2 * 365)) * 2
  et <- pmin(et, 14)
  et <- pmax(et, -20)
  data$event_time <- et
  for (i in 1:10) data[[paste0("trt_m", i)]] <- as.integer(data$event_time == (-i * 2))
  for (i in 0:7) data[[paste0("trt_p", i)]] <- as.integer(data$event_time == (i * 2))
  data$trt_m1 <- 0L
  data$vax_year <- as.integer(format(vax_date, "%Y"))
  data
}

ext_data <- build_event_time_extended(primary0, "zost1_date")
regs_ext <- c(paste0("trt_m", 10:2), paste0("trt_p", 0:7))
mod_ext <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(regs_ext, collapse=" + "),
                                            " | subject + visit_year + visit_age")),
                          data = ext_data, cluster = ~vax_year)
tbl_ext <- broom::tidy(mod_ext, conf.int = TRUE)
readr::write_csv(tbl_ext, "output/tables/TableS4_extended_all_periods.csv")
cat(sprintf("=== Figure S4: extended to all observed periods (N=%d) ===\n", nobs(mod_ext)))
print(tbl_ext, n = Inf)
p_ext <- tidy_plot(tbl_ext, "Figure S4. Extended to all observed event-time periods")
ggsave("output/figures/FigureS4_extended_all_periods.png", p_ext, width = 10, height = 5.5, dpi = 200)

# ------------------------------------------------------------------------
# (2) Negative control: other vaccines (pneumococcal, Td)
#
# Sample: all ever-recipients of the given vaccine among the ZVL-eventually-
# vaccinated cohort, with no calendar-year restriction.
# ------------------------------------------------------------------------
run_other_vax <- function(date_col, label) {
  d <- panel %>%
    filter(!is.na(.data[[date_col]]), !is.na(casi_total)) %>%
    build_event_time(date_col)
  res <- run_event_study(d, "casi_total")
  tbl <- es_table(res$mod)
  readr::write_csv(tbl, sprintf("output/tables/TableS_negctrl_vax_%s.csv", label))
  cat(sprintf("\n=== Negative control vaccine: %s (N=%d) ===\n", label, nobs(res$mod)))
  print(tbl)
  list(tbl = tbl, pretrend_p = res$pretrend_p)
}
res_pneu <- run_other_vax("pneu1_date", "pneumococcal")
res_td   <- run_other_vax("td1_date", "td")

p_othervax <- panel_grid(
  list(list(tbl = res_pneu$tbl, label = "Pneumococcal vaccine"),
       list(tbl = res_td$tbl, label = "Tetanus-diphtheria (Td) vaccine")),
  ncol = 2
)
ggsave("output/figures/FigureS7_other_vaccines.png", p_othervax, width = 10, height = 5, dpi = 200)

# ------------------------------------------------------------------------
# (3) Negative control outcomes: the fixed outcome list used in the paper.
#   Figure S9, panels A-E (chronic conditions): hypertension, cataract,
#     rheumatoid arthritis/osteoarthritis, anemia, hyperlipidemia
#   Figure S10, panels A-E (hospitalizations): osteoarthritis, pneumonia,
#     congestive heart failure, acute cerebrovascular disease, cardiac
#     dysrhythmias
#   Figure S9, panel F (Charlson Comorbidity Index, excl. dementia):
#     cci_excludedementia = max(cci - dementia, 0), min-max normalized;
#     constructed in 00_data_prep.R from the raw cci_* event codes using
#     the standard 17-item Charlson/Deyo weights.
# ------------------------------------------------------------------------
primary <- primary0 %>% build_event_time("zost1_date")

negctrl_ccw_vars <- c("ccw_htn", "ccw_cata", "ccw_raoa", "ccw_anem", "ccw_hlip")
negctrl_ip_vars  <- c("ip_oarth", "ip_pneum", "ip_conhf", "ip_acbvd", "ip_dysrh")

# Negative control outcomes (chronic conditions, hospitalizations) match the
# original specification (negativecontroloutcomes.R): event time capped at
# [-12, +12] (not the primary's [-14, +12], with tail periods binned into
# the -12/+12 endpoints), and no requirement that CASI Total Score be
# non-missing, since CASI is not the outcome for these negative controls.
negctrl_base <- panel %>% filter(!is.na(zost1_date)) %>%
  mutate(
    event_time_raw = floor(as.numeric(visit_date - zost1_date) / (2 * 365)) * 2,
    event_time = pmin(pmax(event_time_raw, -12), 12)
  )
for (i in 1:6) negctrl_base[[paste0("trt_m", i)]] <- as.integer(negctrl_base$event_time == (-i * 2))
for (i in 0:6) negctrl_base[[paste0("trt_p", i)]] <- as.integer(negctrl_base$event_time == (i * 2))
negctrl_base$trt_m1 <- 0L

negctrl_outcome <- function(outcome) {
  d <- negctrl_base %>% filter(!is.na(.data[[outcome]]))
  res <- run_event_study(d, outcome)
  tbl <- es_table(res$mod)
  readr::write_csv(tbl, sprintf("output/tables/TableS_negctrl_outcome_%s.csv", outcome))
  list(outcome = outcome, n = nobs(res$mod), pretrend_p = res$pretrend_p, tbl = tbl)
}
negctrl_ccw_labs <- c(ccw_htn = "Hypertension", ccw_cata = "Cataract", ccw_raoa = "Rheumatoid arthritis/OA",
                       ccw_anem = "Anemia", ccw_hlip = "Hyperlipidemia")
negctrl_ip_labs  <- c(ip_oarth = "Osteoarthritis", ip_pneum = "Pneumonia", ip_conhf = "Congestive heart failure",
                       ip_acbvd = "Acute cerebrovascular disease", ip_dysrh = "Cardiac dysrhythmias")
negctrl_ccw <- lapply(negctrl_ccw_vars, negctrl_outcome)
negctrl_ip  <- lapply(negctrl_ip_vars, negctrl_outcome)
for (i in seq_along(negctrl_ccw)) negctrl_ccw[[i]]$label <- unname(negctrl_ccw_labs[negctrl_ccw_vars[i]])
for (i in seq_along(negctrl_ip))  negctrl_ip[[i]]$label  <- unname(negctrl_ip_labs[negctrl_ip_vars[i]])
cat("\n=== Negative control chronic-condition / hospitalization outcomes: trt_p6 (12y) estimates ===\n")
for (r in c(negctrl_ccw, negctrl_ip)) {
  p6 <- r$tbl %>% filter(term == "trt_p6")
  cat(sprintf("%-15s N=%-6d  12y est=%.4f (p=%.3f)\n", r$outcome, r$n, p6$estimate, p6$pvalue))
}

p_negctrl_ccw <- panel_grid(negctrl_ccw, ytitle_fun = function(l) "Effect on outcome", ncol = 3)
ggsave("output/figures/FigureS9_negctrl_chronic_conditions.png", p_negctrl_ccw, width = 14, height = 6, dpi = 200)
p_negctrl_ip <- panel_grid(negctrl_ip, ytitle_fun = function(l) "Effect on outcome", ncol = 3)
ggsave("output/figures/FigureS10_negctrl_hospitalizations.png", p_negctrl_ip, width = 14, height = 6, dpi = 200)

# ------------------------------------------------------------------------
# (3b) Figure S9, panel F: Charlson Comorbidity Index, excluding dementia --
# negative control (a sixth panel of the chronic-conditions negative-control
# figure, not a separate figure). Constructed directly in 00_data_prep.R
# from the raw cci_* event codes using the standard 17-item Charlson/Deyo
# weights; cci_excludedementia = max(cci - dementia, 0), min-max normalized.
# ------------------------------------------------------------------------
negctrl_cci <- negctrl_outcome("cci_excludedementia")
negctrl_cci$label <- "Charlson Comorbidity Index (excl. dementia)"
cat(sprintf("\ncci_excludedementia  N=%-6d  12y est=%.4f (p=%.3f)  pretrend p=%.4f\n",
            negctrl_cci$n, (negctrl_cci$tbl %>% filter(term=="trt_p6"))$estimate,
            (negctrl_cci$tbl %>% filter(term=="trt_p6"))$pvalue, negctrl_cci$pretrend_p))
p_negctrl_cci <- es_plot(negctrl_cci$tbl, title = "",
                          ytitle = "Effect on normalized CCI (excl. dementia)")
ggsave("output/figures/FigureS9f_negctrl_cci.png", p_negctrl_cci, width = 8, height = 5.5, dpi = 200)
# note for future runs: output/REPLICATION_STATUS.md and output/VALIDATION_TABLE.md
# are hand-maintained audit docs, not generated by this pipeline -- `rm -rf output`
# before a clean re-run deletes them too; recreate/restore from version control.

# ------------------------------------------------------------------------
# (4) Falsification test: pseudo-vaccination dates (pre-2006 period only)
# ------------------------------------------------------------------------
set.seed(20261015)
pseudo_base <- panel %>%
  filter(!is.na(zost1_date)) %>%
  distinct(subject) %>%
  left_join(panel %>% group_by(subject) %>% summarise(baseline_date = min(visit_date), .groups="drop"),
            by = "subject") %>%
  mutate(
    cutoff = as.Date("2006-01-01"),
    span_days = as.numeric(cutoff - baseline_date),
    span_days = pmax(span_days, 1),
    pseudo_date = baseline_date + as.integer(runif(n()) * span_days)
  )

pseudo_dat <- panel %>%
  inner_join(pseudo_base %>% select(subject, pseudo_date), by = "subject") %>%
  filter(visit_date < as.Date("2006-01-01"), !is.na(casi_total)) %>%
  build_event_time("pseudo_date")

res_pseudo <- run_event_study(pseudo_dat, "casi_total")
tbl_pseudo <- es_table(res_pseudo$mod)
readr::write_csv(tbl_pseudo, "output/tables/TableS_falsification_pseudo_dates.csv")
cat(sprintf("\n=== Falsification test (pseudo-vaccination dates), N=%d ===\n", nobs(res_pseudo$mod)))
print(tbl_pseudo)
cat(sprintf("Pretrend p: %.3f\n", res_pseudo$pretrend_p))
p_pseudo <- es_plot(tbl_pseudo, title = "")
ggsave("output/figures/FigureS11_falsification_pseudo_dates.png", p_pseudo, width = 8, height = 5.5, dpi = 200)

# ------------------------------------------------------------------------
# (5) Attrition check: observed from -6 to +6 years, endpoints binned there
# ------------------------------------------------------------------------
primary_6 <- primary0 %>%
  mutate(
    event_time = floor(as.numeric(visit_date - zost1_date) / (2*365)) * 2,
    event_time = pmin(event_time, 6),
    event_time = pmax(event_time, -6)
  )
elig_6 <- primary_6 %>%
  group_by(subject) %>%
  summarise(has_m6 = any(event_time <= -6), has_p6 = any(event_time >= 6), .groups = "drop") %>%
  filter(has_m6, has_p6) %>% pull(subject)

d6 <- primary_6 %>% filter(subject %in% elig_6) %>%
  mutate(
    trt_m1 = 0L,
    trt_m2 = as.integer(event_time == -2), trt_m3 = as.integer(event_time == -4),
    trt_p0 = as.integer(event_time == 0), trt_p1 = as.integer(event_time == 2),
    trt_p2 = as.integer(event_time == 4), trt_p3 = as.integer(event_time == 6)
  )
regs_6 <- c("trt_m3","trt_m2","trt_p0","trt_p1","trt_p2","trt_p3")
mod6 <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(regs_6, collapse=" + "),
                                         " | subject + visit_year + visit_age")),
                       data = d6, vcov = "hetero")
tbl6 <- broom::tidy(mod6, conf.int = TRUE)
readr::write_csv(tbl6, "output/tables/TableS_attrition_restricted.csv")
cat(sprintf("\n=== Attrition-restricted (-6 to +6y), N individuals=%d ===\n", length(elig_6)))
print(tbl6, n = Inf)
p_attr <- tidy_plot(tbl6, "Figure S. Restricted to individuals observed in all event-time periods")
ggsave("output/figures/FigureS12_attrition_restricted.png", p_attr, width = 8, height = 5.5, dpi = 200)

saveRDS(list(extended = tbl_ext, pneu = res_pneu, td = res_td,
             negctrl_ccw = negctrl_ccw, negctrl_ip = negctrl_ip, negctrl_cci = negctrl_cci,
             falsification = tbl_pseudo, attrition = tbl6),
        "output/tables/robustness_results.rds")
cat("\nAll robustness-check result tables saved to output/tables/\n")
