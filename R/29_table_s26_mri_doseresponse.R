# ==============================================================================
# 29_table_s26_mri_doseresponse.R
#
# Table S26 / Figure S29: structural brain MRI, within-ZVL dose-response.
# Ported from Xi/v0926/mri_specification_ladder.R exactly -- among ZVL-
# vaccinated participants with a research-protocol MRI scan, association
# between years since vaccination and 7 structural outcomes (WM
# hyperintensity volume from two pipelines, total brain volume, cortical
# WM volume total/left/right, chronic infarct), under two covariate specs.
# Raw MRI data (not part of this package's other raw extracts) copied from
# Dropbox/Research/ACT study .../2026-09 ACT data for proposal 398 - MRI
# (from Pascal)/ into data/raw/.
# ==============================================================================
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(readr); library(lubridate)
  library(splines); library(ggplot2)
})

amvols <- read_csv("data/raw/ACT_398_mri_amvols_2026_09_02.csv", show_col_types = FALSE) %>%
  mutate(Scan_Date = mdy(Scan_Date))
cdes <- read_csv("data/raw/ACT_398_mri_cdes_2026_09_02.csv", show_col_types = FALSE) %>%
  mutate(Scan_Date = mdy(Scan_Date))

raw <- read_dta("data/raw/event_all_120125.dta")
raw$zost1_date <- as.Date(ifelse(raw$zost1_date == 0, NA, raw$zost1_date), origin = "1960-01-01")
raw$visit_date <- as.Date(ifelse(raw$visit_date == 0, NA, raw$visit_date), origin = "1960-01-01")

vax_group <- raw %>%
  group_by(subject) %>%
  summarise(
    ever_zvl = any(!is.na(zost1_date)),
    ever_rzv = any(!is.na(vaxsh_year1)),
    zost1_date = suppressWarnings(min(zost1_date, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    hz_group = case_when(ever_zvl ~ "ZVL", ever_rzv & !ever_zvl ~ "RZV", TRUE ~ "Never"),
    zost1_date = as.Date(ifelse(is.finite(zost1_date), zost1_date, NA), origin = "1970-01-01")
  ) %>%
  filter(hz_group == "ZVL")

age_anchor <- raw %>%
  filter(!is.na(visit_date), !is.na(visit_age)) %>%
  group_by(subject) %>%
  slice_min(visit_date, n = 1, with_ties = FALSE) %>%
  transmute(subject, anchor_date = visit_date, anchor_age = as.numeric(visit_age))

health_covariates <- raw %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  summarise(
    diabetes = suppressWarnings(as.numeric(first(na.omit(diabetes)))),
    hypertension = suppressWarnings(as.numeric(first(na.omit(hypertension)))),
    stroke = suppressWarnings(as.numeric(first(na.omit(stroke)))),
    bmi = suppressWarnings(as.numeric(first(na.omit(bmi)))),
    smoke_current = suppressWarnings(as.numeric(first(na.omit(smoke_current)))),
    education = suppressWarnings(as.numeric(first(na.omit(education)))),
    female = as.integer(first(na.omit(as.numeric(gender))) == 2),
    white = as.integer(first(na.omit(as.numeric(race_c))) == 7),
    .groups = "drop"
  )

build_analysis_set <- function(df) {
  df %>%
    inner_join(vax_group, by = c("SUBJECT" = "subject")) %>%
    left_join(age_anchor, by = c("SUBJECT" = "subject")) %>%
    left_join(health_covariates, by = c("SUBJECT" = "subject")) %>%
    mutate(
      age_at_scan = anchor_age + as.numeric(Scan_Date - anchor_date) / 365.25,
      scan_year = year(Scan_Date),
      yrs_since_vax = as.numeric(Scan_Date - zost1_date) / 365.25
    ) %>%
    filter(Scan_Type == "R", yrs_since_vax >= 0)
}

am <- build_analysis_set(amvols) %>%
  mutate(
    log_wmh = log(WMH_vol_cm3 + 0.1),
    log_wm_hypointensities = log(vol_WM_hypointensities + 0.1),
    log_cortwm_total = log(vol_Total_cortical_wm_vol + 0.1),
    log_cortwm_left = log(vol_Left_hemi_cortical_wm_vol + 0.1),
    log_cortwm_right = log(vol_Right_hemi_cortical_wm_vol + 0.1)
  )

cd <- build_analysis_set(cdes) %>%
  filter(CDE_Status == "Complete") %>%
  mutate(chronic_infarct = as.numeric(chronic_infarct_present))

# vol_Brain_Seg_Vol_omit_Vent_Surf is identical to vol_Brain_Seg_Vol_omit_Vent (r=1.00);
# ventricle-inclusive vol_Brain_Seg_Vol is a worse atrophy proxy and dropped; Left/Right
# WM & non-WM hypointensities are constant zero in this extract; vol_Non_WM_hypointensities
# is 97.7% zero; acute_infarct has only 2 events (near-complete separation) and is dropped.

COVARIATES <- c("diabetes", "hypertension", "stroke", "smoke_current", "bmi", "education", "female", "white")
SPECS <- list(
  "Age only"       = "ns(age_at_scan, 3)",
  "+ Vascular/SES" = paste("ns(age_at_scan, 3) + ns(scan_year, 3)", paste(COVARIATES, collapse = " + "), sep = " + ")
)
# "+ Vascular/SES" is cumulative. Scan year alone is not reported separately: years
# since vaccination is highly collinear with scan year here (r=0.85), since vaccination
# occurred in a narrow window (2007-2017) while scanning concentrated later (2014-2024).

run_spec_ladder <- function(data, outcome_var, outcome_label, family = gaussian()) {
  needed_cols <- c(outcome_var, "yrs_since_vax", "age_at_scan", COVARIATES)
  complete_case_data <- data %>% filter(complete.cases(across(all_of(needed_cols))))
  lapply(names(SPECS), function(spec_name) {
    formula <- as.formula(paste(outcome_var, "~ yrs_since_vax +", SPECS[[spec_name]]))
    model <- if (family$family == "binomial") {
      glm(formula, data = complete_case_data, family = binomial)
    } else {
      lm(formula, data = complete_case_data)
    }
    coefs <- summary(model)$coefficients
    tibble::tibble(
      outcome = outcome_label, spec = spec_name,
      n_scans = nobs(model), n_subjects = n_distinct(complete_case_data$SUBJECT),
      estimate = coefs["yrs_since_vax", "Estimate"], se = coefs["yrs_since_vax", "Std. Error"],
      p = coefs["yrs_since_vax", ncol(coefs)]
    )
  }) %>% bind_rows()
}

OUTCOME_ORDER <- c(
  "WM hypointensities (log)", "WMH volume (log)", "Total brain volume (cm³)",
  "Total cortical WM volume (log)", "Left cortical WM volume (log)",
  "Right cortical WM volume (log)", "Chronic infarct (log-odds)"
)

results <- bind_rows(
  run_spec_ladder(am, "log_wm_hypointensities", "WM hypointensities (log)"),
  run_spec_ladder(am, "log_wmh", "WMH volume (log)"),
  run_spec_ladder(am, "vol_Brain_Seg_Vol_omit_Vent", "Total brain volume (cm³)"),
  run_spec_ladder(am, "log_cortwm_total", "Total cortical WM volume (log)"),
  run_spec_ladder(am, "log_cortwm_left", "Left cortical WM volume (log)"),
  run_spec_ladder(am, "log_cortwm_right", "Right cortical WM volume (log)"),
  run_spec_ladder(cd, "chronic_infarct", "Chronic infarct (log-odds)", family = binomial())
) %>%
  mutate(
    conf.low = estimate - 1.96 * se, conf.high = estimate + 1.96 * se,
    outcome = factor(outcome, levels = OUTCOME_ORDER),
    spec = factor(spec, levels = rev(names(SPECS)))
  ) %>%
  arrange(outcome, spec)

write_csv(results, "output/tables/TableS26_mri_doseresponse.csv")
cat("\nTable S26 (MRI dose-response):\n"); print(results, n = Inf)

# Robustness check cited in the table notes: clustering SEs by participant for WMH
# volume (age-only) widens p from .064 to .095.
wmh_cc <- am %>% filter(complete.cases(across(all_of(c("log_wmh", "yrs_since_vax", "age_at_scan")))))
wmh_model <- lm(log_wmh ~ yrs_since_vax + ns(age_at_scan, 3), data = wmh_cc)
if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
  clustered_p <- lmtest::coeftest(wmh_model, vcov = sandwich::vcovCL(wmh_model, cluster = wmh_cc$SUBJECT))["yrs_since_vax", 4]
  cat(sprintf("WMH volume, age-only: model-based p = %.3f, participant-clustered p = %.3f\n",
              summary(wmh_model)$coefficients["yrs_since_vax", 4], clustered_p))
}

# ------------------------------------------------------------------------
# Figure S29: forest plot, one facet per outcome, house style.
# ------------------------------------------------------------------------
BENEFIT_IS_POSITIVE <- c(
  "WM hypointensities (log)" = FALSE, "WMH volume (log)" = FALSE,
  "Total brain volume (cm³)" = TRUE, "Total cortical WM volume (log)" = TRUE,
  "Left cortical WM volume (log)" = TRUE, "Right cortical WM volume (log)" = TRUE,
  "Chronic infarct (log-odds)" = FALSE
)

panel_range <- results %>%
  group_by(outcome) %>%
  summarise(lo = min(conf.low), hi = max(conf.high), .groups = "drop") %>%
  mutate(span = hi - lo, gap = span * 0.015, arrow_len = span * 0.18, text_pad = span * 0.02,
         benefit_positive = BENEFIT_IS_POSITIVE[as.character(outcome)])

annotation_arrows <- bind_rows(
  panel_range %>% transmute(outcome, x = -gap, xend = -(gap + arrow_len), y = 2.78),
  panel_range %>% transmute(outcome, x = gap, xend = gap + arrow_len, y = 2.78)
)
annotation_labels <- bind_rows(
  panel_range %>% transmute(outcome, x = -(gap + arrow_len + text_pad), y = 2.78,
                             label = ifelse(benefit_positive, "Harm", "Benefit"), hjust = 1),
  panel_range %>% transmute(outcome, x = gap + arrow_len + text_pad, y = 2.78,
                             label = ifelse(benefit_positive, "Benefit", "Harm"), hjust = 0)
)

p29 <- ggplot(results, aes(x = estimate, y = spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, colour = "grey60") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.55, colour = "grey45") +
  geom_point(size = 2.1, colour = "black") +
  geom_segment(data = annotation_arrows, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.09, "cm"), type = "closed"), colour = "grey30", linewidth = 0.5) +
  geom_text(data = annotation_labels, aes(x = x, y = y, label = label, hjust = hjust), inherit.aes = FALSE,
            size = 2.6, colour = "grey30", fontface = "italic") +
  facet_wrap(~outcome, scales = "free_x", nrow = 2) +
  coord_cartesian(ylim = c(0.4, 3.05), clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.28)) +
  labs(x = "Estimated change per additional year since ZVL vaccination (95% CI)", y = "Covariate adjustment") +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "grey95", colour = "black"),
        strip.text = element_text(size = 8.5, face = "bold"),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        axis.text.y = element_text(size = 9), axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8), plot.margin = margin(5.5, 10, 5.5, 5.5))
ggsave("output/figures/FigureS29_mri_doseresponse.png", p29, width = 10, height = 4.8, dpi = 300, bg = "white")
cat("Saved output/figures/FigureS29_mri_doseresponse.png\n")
