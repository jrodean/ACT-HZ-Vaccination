# ==============================================================================
# 22_chronic_hospitalization_adjusted.R
#
# Figure S16. The effect of ZVL vaccination on CASI Total Score when
# controlling for hospitalization events and new chronic disease diagnoses.
#
# Matches equations (3)-(4) in the Supplement's Technical Details
# ("Regression Specifications"):
#   (3) Adjustment for chronic disease onset: adds absorbing (ever-on)
#       indicators for first diagnosis of diabetes, hypertension,
#       cardiovascular disease, COPD, depression, and cancer.
#   (4) Adjustment for acute medical events: further adds time-varying
#       indicators for serious inpatient hospitalizations (stroke,
#       myocardial infarction, pneumonia, septicemia, gastrointestinal
#       hemorrhage), active only in the year of occurrence.
#
# "Cardiovascular disease" is operationalized as ischemic heart disease
# (ccw_ihd), the standard CCW claims-based proxy; the exact original choice
# among the available cardiovascular condition flags (ccw_ihd, ccw_myoc,
# ccw_hfai) was not confirmed.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

chronic_vars <- c("ccw_diab", "ccw_htn", "ccw_ihd", "ccw_copd", "ccw_depr", "ccw_canc")
acute_vars   <- c("ip_acbvd", "ip_myoca", "ip_pneum", "ip_septi", "ip_ghemo")

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date")

missing_chronic <- setdiff(chronic_vars, names(primary))
missing_acute <- setdiff(acute_vars, names(primary))
if (length(missing_chronic) || length(missing_acute)) {
  stop("Missing covariate columns: ", paste(c(missing_chronic, missing_acute), collapse = ", "))
}

fml_chronic <- as.formula(paste0(
  "casi_total ~ ", paste(es_regs, collapse = " + "), " + ", paste(chronic_vars, collapse = " + "),
  " | subject + visit_year + visit_age"
))
m_chronic <- fixest::feols(fml_chronic, data = primary, vcov = "hetero")
tbl_chronic <- es_table(m_chronic)
readr::write_csv(tbl_chronic, "output/tables/TableS16_chronic_adjusted.csv")

fml_acute <- as.formula(paste0(
  "casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
  paste(c(chronic_vars, acute_vars), collapse = " + "),
  " | subject + visit_year + visit_age"
))
m_acute <- fixest::feols(fml_acute, data = primary, vcov = "hetero")
tbl_acute <- es_table(m_acute)
readr::write_csv(tbl_acute, "output/tables/TableS16_chronic_acute_adjusted.csv")

cat(sprintf("=== Eq. (3) chronic-disease-adjusted (N=%d) ===\n", nobs(m_chronic)))
print(tbl_chronic)
cat(sprintf("\n=== Eq. (4) chronic+acute-adjusted (N=%d) ===\n", nobs(m_acute)))
print(tbl_acute)

res_primary <- run_event_study(primary, "casi_total")
tbl_primary <- es_table(res_primary$mod)

combined <- bind_rows(
  tbl_primary %>% mutate(spec = "Primary (unadjusted)"),
  tbl_chronic %>% mutate(spec = "+ Chronic disease onset (eq. 3)"),
  tbl_acute %>% mutate(spec = "+ Chronic disease + acute events (eq. 4)")
) %>%
  mutate(spec = factor(spec, levels = c("Primary (unadjusted)", "+ Chronic disease onset (eq. 3)",
                                         "+ Chronic disease + acute events (eq. 4)")))

ref_rows <- combined %>% distinct(spec) %>%
  mutate(term = "trt_m1", label = "2 years before", estimate = 0, ci_low = 0, ci_high = 0)
plot_df <- bind_rows(combined, ref_rows)
ord <- c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2","trt_m1",
         "trt_p0","trt_p1","trt_p2","trt_p3","trt_p4","trt_p5","trt_p6")
plot_df$term <- factor(plot_df$term, levels = ord)
plot_df <- plot_df %>% arrange(spec, term) %>% group_by(spec) %>% mutate(x = as.integer(term)) %>% ungroup()

p16 <- ggplot(plot_df, aes(x = x, y = estimate, color = spec, shape = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = seq_along(ord), labels = coef_map_short[ord] %>% { .[is.na(.)] <- "2 yrs before"; . }) +
  labs(x = "Time relative to the year of ZVL vaccination", y = "Effect on CASI Total Score", color = NULL, shape = NULL,
       title = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave("output/figures/FigureS16_chronic_hospitalization_adjusted.png", p16, width = 9, height = 6, dpi = 200)

cat("\nSaved Figure S16 (chronic-disease/hospitalization-adjusted robustness check).\n")
