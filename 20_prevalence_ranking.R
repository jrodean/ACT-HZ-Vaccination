# ==============================================================================
# 20_prevalence_ranking.R
#
# Figure S8. Prevalence ranking of chronic condition diagnoses and primary
# diagnoses for inpatient hospitalizations, in the primary analysis cohort.
# The top five conditions in each category are used as negative-control
# outcomes in Figures S9 and S10.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total))

ccw_labels <- c(
  ccw_cata = "Cataract", ccw_htn = "Hypertension", ccw_raoa = "Rheumatoid arthritis/osteoarthritis",
  ccw_anem = "Anemia", ccw_hlip = "Hyperlipidemia", ccw_afib = "Atrial fibrillation",
  ccw_asth = "Asthma", ccw_bphp = "Benign prostatic hyperplasia", ccw_canb = "Cancer (bladder)",
  ccw_canc = "Cancer (colorectal)", ccw_cane = "Cancer (endometrial)", ccw_canl = "Cancer (lung)",
  ccw_canp = "Cancer (prostate)", ccw_ckd = "Chronic kidney disease", ccw_copd = "COPD",
  ccw_depr = "Depression", ccw_diab = "Diabetes", ccw_glau = "Glaucoma",
  ccw_hfai = "Heart failure", ccw_hipf = "Hip fracture", ccw_hthy = "Hypothyroidism",
  ccw_ihd = "Ischemic heart disease", ccw_myoc = "Myocardial infarction", ccw_oste = "Osteoporosis",
  ccw_stro = "Stroke"
)
ip_labels <- c(
  ip_oarth = "Osteoarthritis", ip_pneum = "Pneumonia", ip_conhf = "Congestive heart failure",
  ip_acbvd = "Acute cerebrovascular disease", ip_dysrh = "Cardiac dysrhythmias"
)

ccw_vars <- intersect(names(ccw_labels), names(primary))
ip_vars  <- intersect(names(ip_labels), names(primary))

prevalence_table <- function(vars, labels) {
  tibble::tibble(
    var = vars,
    label = unname(labels[vars]),
    prevalence = vapply(vars, function(v) mean(primary[[v]] == 1, na.rm = TRUE), numeric(1)),
    n_cases = vapply(vars, function(v) sum(primary[[v]] == 1, na.rm = TRUE), numeric(1))
  ) %>% arrange(desc(prevalence))
}

ccw_prev <- prevalence_table(ccw_vars, ccw_labels)
ip_prev  <- prevalence_table(ip_vars, ip_labels)

readr::write_csv(ccw_prev, "output/tables/FigureS8_ccw_prevalence_ranking_data.csv")
readr::write_csv(ip_prev, "output/tables/FigureS8_ip_prevalence_ranking_data.csv")

cat("=== Top 5 chronic conditions by prevalence (negative controls used in Figure S9) ===\n")
print(head(ccw_prev, 5))
cat("\n=== Top 5 hospitalization categories by prevalence (negative controls used in Figure S10) ===\n")
print(head(ip_prev, 5))

plot_prevalence <- function(df, subtitle) {
  df <- df %>% mutate(label = factor(label, levels = rev(label)))
  ggplot(df, aes(x = label, y = prevalence)) +
    geom_col(fill = "grey40") +
    coord_flip() +
    labs(x = NULL, y = "Prevalence in primary analysis cohort", subtitle = subtitle) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size = 11)
}

p_ccw <- plot_prevalence(ccw_prev, "Chronic condition diagnoses")
p_ip  <- plot_prevalence(ip_prev, "Primary diagnoses for inpatient hospitalizations")
p_combined <- patchwork::wrap_plots(p_ccw, p_ip, ncol = 1)
ggsave("output/figures/FigureS8_prevalence_ranking.png", p_combined, width = 9, height = 10, dpi = 200)

cat("\nSaved Figure S8 (prevalence ranking) tables and figure.\n")

# Cross-check against the actual outcome list used in Figure S9
# (04_robustness_checks.R / negctrl_ccw_vars): hypertension, cataract,
# rheumatoid arthritis/osteoarthritis, anemia, hyperlipidemia. This
# reconstruction's top-5-by-prevalence ranking agrees on 4/5 (hypertension,
# cataract, RA/OA, hyperlipidemia) but ranks glaucoma above anemia, while
# Figure S9 uses anemia. The hospitalization ranking (Figure S10) matches
# exactly on all 5. This discrepancy does not affect Figures S9/S10
# themselves (their outcome list is fixed independently of this script),
# only this descriptive justification figure; the original prevalence
# computation underlying the published Figure S8 was not located to
# confirm whether it used a different cohort or a documented clinical
# exclusion of glaucoma.
cat("\nNOTE: reconstructed CCW ranking disagrees with the actual Figure S9",
    "outcome list at rank 5 (glaucoma here vs. anemia in Figure S9) --",
    "see comment above this line.\n")
