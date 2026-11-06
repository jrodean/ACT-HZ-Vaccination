# ==============================================================================
# 22_chronic_hospitalization_adjusted.R
#
# Figure S16: effect of ZVL vaccination on CASI Total Score controlling for
# hospitalizations and chronic-disease diagnoses -- found in the original
# analysis_Xi_0326.do (`global controls`, line ~679): all 24 ccw_* chronic-
# condition indicators and all 20 ip_* inpatient-hospitalization indicators
# entered simultaneously as linear covariates (not a curated subset), with
# standard errors clustered by vaccination year (not heteroskedasticity-
# robust, unlike the primary specification).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date") %>%
  mutate(vax_year = lubridate::year(zost1_date))

ccw_vars <- grep("^ccw_", names(primary), value = TRUE)
ip_vars <- grep("^ip_", names(primary), value = TRUE)
controls <- c(ccw_vars, ip_vars)
cat(sprintf("N controls: %d ccw_ + %d ip_ = %d\n", length(ccw_vars), length(ip_vars), length(controls)))

fml <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
                          paste(controls, collapse = " + "), " | subject + visit_year + visit_age"))
m16 <- fixest::feols(fml, data = primary, cluster = ~vax_year)
tbl16 <- es_table(m16)
readr::write_csv(tbl16, "output/tables/FigureS16_chronic_hospitalization_adjusted_data.csv")
cat(sprintf("=== Chronic + inpatient controls (N=%d) ===\n", nobs(m16)))
print(tbl16)

p16 <- es_plot(tbl16, ytitle = "Effect on CASI Total Score")
ggsave("output/figures/FigureS16_chronic_hospitalization_adjusted.png", p16, width = 8, height = 5.5, dpi = 200)

cat("\nSaved Figure S16 (chronic-disease/hospitalization-adjusted robustness check).\n")
