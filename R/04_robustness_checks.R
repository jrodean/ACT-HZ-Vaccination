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
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(haven) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# simple plot for broom::tidy()-style tables (term/estimate/conf.low/conf.high),
# in the package's house style (theme_classic, no gridlines, Benefit/Harm
# arrow, red dashed vaccination-year line). term names follow the trt_m{n}/
# trt_p{n} convention, so labels and the omitted trt_m1 reference point (not
# estimated, plotted at zero) are derived from the term names directly.
term_label <- function(term) {
  n <- as.integer(sub("trt_[mp]", "", term))
  dplyr::case_when(term == "trt_p0" ~ "ZVL vaccination year",
                    grepl("^trt_m", term) ~ paste(2 * n, "yrs before"),
                    TRUE ~ paste(2 * n, "yrs after"))
}
tidy_plot <- function(tbl, ytitle = "Effect on CASI Total Score") {
  ref_row <- tibble::tibble(term = "trt_m1", estimate = 0, conf.low = 0, conf.high = 0)
  d <- dplyr::bind_rows(tbl, ref_row)
  ord <- c(paste0("trt_m", 20:1), paste0("trt_p", 0:20))
  d <- d[order(match(d$term, ord)), ]
  d$x <- seq_len(nrow(d))
  d$lab <- term_label(d$term)
  vax_x <- d$x[d$term == "trt_p0"]
  p <- ggplot(d, aes(x = x, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = vax_x, linetype = "dashed", color = "red") +
    benefit_harm_annotation(d %>% dplyr::rename(ci_low = conf.low, ci_high = conf.high)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, linetype = "dashed", color = "grey50") +
    geom_point(size = 2) +
    scale_x_continuous(breaks = d$x, labels = d$lab) +
    labs(x = "Time relative to the year of ZVL vaccination", y = ytitle, title = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank())
  p
}

primary0 <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total))

vax_frozen <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  mutate(visit_date = as.Date(visit_date, origin = "1960-01-01"))

# ------------------------------------------------------------------------
# (1) Figure S4: extended to all observed periods. Ported from the
# original Stata (analysis_Xi_1225.do, "Robustness: full periods"):
# loads the frozen event_all_120125.dta directly, 2-year bins using
# 365-day years (not 365.25), SEs clustered by vaccination year. The
# dummies are built from the RAW (uncapped) event_time -- the original's
# "replace event_time=14/-20" comes AFTER trt_m1..trt_m11/trt_p0..trt_p7
# are already generated and only caps the (unused) displayed variable,
# so observations beyond -20/+14 fall into the same omitted reference
# bucket as -2 yrs (trt_m1, forced to 0) rather than into trt_m10/trt_p7.
# ------------------------------------------------------------------------
s4_frozen <- vax_frozen %>%
  mutate(zost1_date = as.Date(zost1_date, origin = "1960-01-01")) %>%
  filter(!is.na(casi_total), !is.na(zost1_date))

build_event_time_extended <- function(data, vax_date_col) {
  vax_date <- data[[vax_date_col]]
  et <- floor(as.numeric(data$visit_date - vax_date) / (2 * 365)) * 2
  data$event_time <- et
  for (i in 1:10) data[[paste0("trt_m", i)]] <- as.integer(et == (-i * 2))
  for (i in 0:7) data[[paste0("trt_p", i)]] <- as.integer(et == (i * 2))
  data$trt_m1 <- 0L
  data$vax_year <- as.integer(format(vax_date, "%Y"))
  data
}

ext_data <- build_event_time_extended(s4_frozen, "zost1_date")
regs_ext <- c(paste0("trt_m", 10:2), paste0("trt_p", 0:7))
mod_ext <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(regs_ext, collapse=" + "),
                                            " | subject + visit_year + visit_age")),
                          data = ext_data, cluster = ~vax_year)
tbl_ext <- broom::tidy(mod_ext, conf.int = TRUE)
readr::write_csv(tbl_ext, "output/tables/TableS4_extended_all_periods.csv")
cat(sprintf("=== Figure S4: extended to all observed periods (N=%d) ===\n", nobs(mod_ext)))
print(tbl_ext, n = Inf)
p_ext <- tidy_plot(tbl_ext)
ggsave("output/figures/FigureS4_extended_all_periods.png", p_ext, width = 10, height = 5.5, dpi = 200)

# ------------------------------------------------------------------------
# (2) Negative control: other vaccines (pneumococcal, Td) -- ported from
# Xi/v0326/othervaccines.do exactly: no zost1_date restriction (the
# original explicitly comments that filter OUT), vaccination year in
# [2006,2017], event_time capped at +12 but DROPPED (not capped) below
# -14, on the original's own frozen event_all_120125.dta.
# ------------------------------------------------------------------------
run_other_vax <- function(date_col, label) {
  d <- vax_frozen %>%
    mutate(vdate = as.Date(.data[[date_col]], origin = "1960-01-01")) %>%
    filter(!is.na(vdate)) %>%
    mutate(vyear = as.integer(format(vdate, "%Y"))) %>%
    filter(vyear >= 2006, vyear < 2018) %>%
    mutate(event_time = floor(as.numeric(visit_date - vdate) / (2*365)) * 2,
           event_time = pmin(event_time, 12)) %>%
    filter(event_time >= -14)
  for (i in 1:6) d[[paste0("trt_m", i)]] <- as.integer(d$event_time == (-i * 2))
  for (i in 0:6) d[[paste0("trt_p", i)]] <- as.integer(d$event_time == (i * 2))
  d$trt_m1 <- 0L
  res <- run_event_study(d, "casi_total")
  tbl <- es_table(res$mod)
  readr::write_csv(tbl, sprintf("output/tables/TableS_negctrl_vax_%s.csv", label))
  cat(sprintf("\n=== Negative control vaccine: %s (N participants=%d, N visits=%d) ===\n",
              label, length(unique(d$subject)), nobs(res$mod)))
  print(tbl)
  list(tbl = tbl, pretrend_p = res$pretrend_p)
}
res_pneu <- run_other_vax("vaxpneu_date", "pneumococcal")
res_td   <- run_other_vax("vaxtd_date", "td")

othervax_xtitle <- c("Pneumococcal vaccine" = "Time relative to the year of pneumococcal vaccination",
                      "Tetanus-diphtheria (Td) vaccine" = "Time relative to the year of Td vaccination")
othervax_vax0 <- c("Pneumococcal vaccine" = "Pneumococcal vaccination year",
                    "Tetanus-diphtheria (Td) vaccine" = "Td vaccination year")
p_othervax <- panel_grid(
  list(list(tbl = res_pneu$tbl, label = "Pneumococcal vaccine"),
       list(tbl = res_td$tbl, label = "Tetanus-diphtheria (Td) vaccine")),
  ytitle_fun = function(l) "Effect on CASI Total Score",
  xtitle_fun = function(l) othervax_xtitle[[l]],
  vax0_label_fun = function(l) othervax_vax0[[l]],
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

# order matches the live Supplement's Figure S9 panel lettering: (A) Cataract, (B) Hypertension, (C) RA/OA, (D) Anemia, (E) Hyperlipidemia
negctrl_ccw_vars <- c("ccw_cata", "ccw_htn", "ccw_raoa", "ccw_anem", "ccw_hlip")
negctrl_ip_vars  <- c("ip_oarth", "ip_pneum", "ip_conhf", "ip_acbvd", "ip_dysrh")

# Negative control outcomes (chronic conditions, hospitalizations) use event
# time capped at [-12, +12] (not the primary's [-14, +12], with tail periods
# binned into the -12/+12 endpoints), and no requirement that CASI Total
# Score be non-missing, since CASI is not the outcome for these negative
# controls. Matches the original negativecontroloutcomes.R specification:
# i(event_time, ref=-2, keep=seq(-12,12,2)) | subject+visit_year+visit_age,
# vcov="hetero", event time binned with a 365.25-day year. The source data
# is NOT event_all_120125.dta alone -- the original builds its working
# event_all from event_all2.dta (visit years through 2022) with
# event_all_120125.dta's visit_year>=2023 rows appended on top, then drops
# visits after each subject's end of follow-up (days_to_withdraw/
# days_to_death/days_to_enrend vs. days since baseline_date) -- both of
# which the package was previously missing (it used event_all_120125.dta
# alone, uncensored, which inflated later post-vaccination estimates,
# worst for ccw_hlip). d_ccw_*->ccw_* rename (not the raw non-binary
# ccw_* code) and ip_* non-missing->1 are unchanged from before.
negctrl_cols <- c("subject", "visit_date", "visit_year", "visit_age", "zost1_date",
                   "casi_total", "baseline_date", "cci", "dementia",
                   "days_to_withdraw", "days_to_death", "days_to_enrend",
                   paste0("d_", negctrl_ccw_vars), negctrl_ip_vars)
e2 <- haven::read_dta("data/raw/event_all2.dta") %>% select(all_of(negctrl_cols)) %>%
  mutate(zost1_date = as.Date(as.numeric(zost1_date), origin = "1960-01-01"),
         baseline_date = as.Date(as.numeric(baseline_date), origin = "1960-01-01"))
e1_new <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  filter(visit_year >= 2023) %>% select(all_of(negctrl_cols)) %>%
  mutate(across(all_of(c("zost1_date", "baseline_date", negctrl_ip_vars)),
                ~ as.Date(as.numeric(.x), origin = "1960-01-01")))
event_all_full <- bind_rows(e2, e1_new)
cci_ed_min <- min(pmax(event_all_full$cci - event_all_full$dementia, 0), na.rm = TRUE)
cci_ed_max <- max(pmax(event_all_full$cci - event_all_full$dementia, 0), na.rm = TRUE)
event_all_nc <- event_all_full %>%
  mutate(days_to_end = pmin(days_to_withdraw, days_to_death, days_to_enrend, na.rm = TRUE),
         days_to_end = ifelse(is.infinite(days_to_end), NA, days_to_end),
         days_since_baseline = as.integer(visit_date - baseline_date)) %>%
  filter(is.na(days_to_end) | days_since_baseline <= days_to_end) %>%
  filter(!is.na(zost1_date)) %>%
  mutate(
    event_time_raw = floor(as.numeric(visit_date - zost1_date) / (2 * 365.25)) * 2,
    event_time = pmin(pmax(event_time_raw, -12), 12)
  ) %>%
  mutate(across(all_of(paste0("d_", negctrl_ccw_vars)), ~ .x, .names = "{sub('^d_', '', .col)}")) %>%
  mutate(across(all_of(negctrl_ip_vars), ~ as.integer(!is.na(.x))))
negctrl_base <- event_all_nc
for (i in 1:6) negctrl_base[[paste0("trt_m", i)]] <- as.integer(negctrl_base$event_time == (-i * 2))
for (i in 0:6) negctrl_base[[paste0("trt_p", i)]] <- as.integer(negctrl_base$event_time == (i * 2))
negctrl_base$trt_m1 <- 0L

negctrl_outcome <- function(outcome, base = negctrl_base, file_suffix = outcome) {
  d <- base %>% filter(!is.na(.data[[outcome]]))
  res <- run_event_study(d, outcome)
  tbl <- es_table(res$mod)
  readr::write_csv(tbl, sprintf("output/tables/TableS_negctrl_outcome_%s.csv", file_suffix))
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

# Figure S9 panel E (Hyperlipidemia) in the live Supplement is the CASI
# Total Score>=90-restricted version, not the full-sample ccw_hlip above --
# confirmed by direct point-for-point match against the real published
# panel (full-sample estimates run ~1.5-2x too negative in the post-
# vaccination tail; the CASI>=90-restricted ones match closely). Swap it
# in for the panel-grid entry; the full-sample version is still written to
# its own CSV above for the record.
eligible_casi90 <- event_all_nc %>%
  filter(!is.na(subject), !is.na(casi_total)) %>%
  group_by(subject) %>%
  summarise(min_casi = min(casi_total, na.rm = TRUE), .groups = "drop") %>%
  filter(min_casi >= 90) %>%
  pull(subject)
negctrl_hlip_casi90 <- negctrl_outcome("ccw_hlip", base = negctrl_base %>% filter(subject %in% eligible_casi90),
                                        file_suffix = "ccw_hlip_casi90")
negctrl_hlip_casi90$label <- "Hyperlipidemia"
negctrl_ccw[[which(negctrl_ccw_vars == "ccw_hlip")]] <- negctrl_hlip_casi90
cat("\n=== Negative control chronic-condition / hospitalization outcomes: trt_p6 (12y) estimates ===\n")
for (r in c(negctrl_ccw, negctrl_ip)) {
  p6 <- r$tbl %>% filter(term == "trt_p6")
  cat(sprintf("%-15s N=%-6d  12y est=%.4f (p=%.3f)\n", r$outcome, r$n, p6$estimate, p6$pvalue))
}

p_negctrl_ccw <- panel_grid(negctrl_ccw, ytitle_fun = function(l) "Percentage points", ncol = 3,
                             benefit_harm = FALSE, border = TRUE, ylim = c(-0.2, 0.2))
ggsave("output/figures/FigureS9_negctrl_chronic_conditions.png", p_negctrl_ccw, width = 14, height = 6, dpi = 200)
p_negctrl_ip <- panel_grid(negctrl_ip, ytitle_fun = function(l) "Percentage points", ncol = 3,
                            benefit_harm = FALSE, border = TRUE, ylim = c(-0.2, 0.2))
ggsave("output/figures/FigureS10_negctrl_hospitalizations.png", p_negctrl_ip, width = 14, height = 6, dpi = 200)

# ------------------------------------------------------------------------
# (3b) Figure S9, panel F: Charlson Comorbidity Index, excluding dementia --
# negative control (a sixth panel of the chronic-conditions negative-control
# figure, not a separate figure). Ported from negativecontroloutcomes.R:
# cci_excludedementia = max(cci - dementia, 0), min-max normalized over the
# SAME combined event_all2.dta+event_all_120125.dta(2023+) source as (3)
# above (cci_ed_min/max computed pre-filtering, same as the original), not
# the frozen file alone.
# ------------------------------------------------------------------------
negctrl_base_cci <- event_all_nc %>%
  mutate(cci_excludedementia = (pmax(cci - dementia, 0) - cci_ed_min) / (cci_ed_max - cci_ed_min))
for (i in 1:6) negctrl_base_cci[[paste0("trt_m", i)]] <- as.integer(negctrl_base_cci$event_time == (-i * 2))
for (i in 0:6) negctrl_base_cci[[paste0("trt_p", i)]] <- as.integer(negctrl_base_cci$event_time == (i * 2))
negctrl_base_cci$trt_m1 <- 0L
negctrl_cci <- negctrl_outcome("cci_excludedementia", base = negctrl_base_cci)
negctrl_cci$label <- "Charlson Comorbidity Index (excl. dementia)"
cat(sprintf("\ncci_excludedementia  N=%-6d  12y est=%.4f (p=%.3f)  pretrend p=%.4f\n",
            negctrl_cci$n, (negctrl_cci$tbl %>% filter(term=="trt_p6"))$estimate,
            (negctrl_cci$tbl %>% filter(term=="trt_p6"))$pvalue, negctrl_cci$pretrend_p))
p_negctrl_cci <- es_plot(negctrl_cci$tbl, title = "",
                          ytitle = "Effect on normalized CCI (excl. dementia)", benefit_harm = FALSE) +
  theme(panel.border = element_rect(color = "black", fill = NA))
ggsave("output/figures/FigureS9f_negctrl_cci.png", p_negctrl_cci, width = 8, height = 5.5, dpi = 200)

# ------------------------------------------------------------------------
# (3c) Standalone copy of the CASI>=90-restricted hyperlipidemia panel
# (already computed above and used as Figure S9 panel E), saved on its own
# axis as a companion figure.
# ------------------------------------------------------------------------
cat(sprintf("\nccw_hlip (CASI>=90 restricted, = Figure S9 panel E)  N=%-6d  12y est=%.4f (p=%.3f)\n",
            negctrl_hlip_casi90$n,
            (negctrl_hlip_casi90$tbl %>% filter(term == "trt_p6"))$estimate,
            (negctrl_hlip_casi90$tbl %>% filter(term == "trt_p6"))$pvalue))
p_negctrl_hlip_casi90 <- es_plot(negctrl_hlip_casi90$tbl, title = "",
                                  ytitle = "Percentage points", benefit_harm = FALSE) +
  theme(panel.border = element_rect(color = "black", fill = NA))
ggsave("output/figures/FigureS9g_negctrl_hlip_casi90.png", p_negctrl_hlip_casi90,
       width = 8, height = 5.5, dpi = 200)

# ------------------------------------------------------------------------
# (4) Falsification test: pseudo-vaccination dates (pre-2006 period only).
# Ported exactly from the original Xi/v0926/reviewer1.R, PART D ("random
# pseudo-vaccination-date falsification", the Fig S7/S11 analysis): seed 1,
# one subject's own earliest visit date as b0, subjects restricted to
# b0 < 2006-01-01, one runif(1) draw per subject via a sorted group_by.
# ------------------------------------------------------------------------
set.seed(1)
pseudo_dat <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  group_by(subject) %>%
  mutate(b0 = suppressWarnings(min(visit_date, na.rm = TRUE))) %>%
  filter(b0 < as.Date("2006-01-01")) %>%
  mutate(pseudo_date = b0 + runif(1) * as.numeric(as.Date("2006-01-01") - b0)) %>%
  ungroup() %>%
  filter(visit_date < as.Date("2006-01-01")) %>%
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
# (5) Figure S12 -- "balanced panel [-6,+6]" block, ported literally from
# analysis_Xi_1225.do's "Controlling for Compositional Changes over time"
# section: sample=1 if max_et>6 & min_et<-6 (strict), with max_et/min_et
# computed from event_all_120125.dta itself (same working data, not a
# fuller external visit history), reghdfe ... vce(robust).
# ------------------------------------------------------------------------
raw_frozen <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  mutate(zost1_date = as.Date(zost1_date, origin = "1960-01-01"),
         visit_date  = as.Date(visit_date,  origin = "1960-01-01"))
primary_6 <- raw_frozen %>% filter(!is.na(zost1_date)) %>%
  mutate(event_time = floor(as.numeric(visit_date - zost1_date) / (2*365)) * 2)
elig_6 <- primary_6 %>%
  group_by(subject) %>%
  summarise(min_et = min(event_time), max_et = max(event_time), .groups = "drop") %>%
  filter(max_et > 6, min_et < -6) %>% pull(subject)

d6 <- primary_6 %>% filter(subject %in% elig_6) %>%
  mutate(
    trt_m1 = 0L,
    trt_m2 = as.integer(event_time == -4), trt_m3 = as.integer(event_time == -6),
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
p_attr <- tidy_plot(tbl6)
ggsave("output/figures/FigureS12_attrition_restricted.png", p_attr, width = 8, height = 5.5, dpi = 200)

saveRDS(list(extended = tbl_ext, pneu = res_pneu, td = res_td,
             negctrl_ccw = negctrl_ccw, negctrl_ip = negctrl_ip, negctrl_cci = negctrl_cci,
             falsification = tbl_pseudo, attrition = tbl6),
        "output/tables/robustness_results.rds")
cat("\nAll robustness-check result tables saved to output/tables/\n")
