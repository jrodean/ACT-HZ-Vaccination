# ==============================================================================
# 20_prevalence_ranking.R
#
# Figure S8. Prevalence ranking of chronic condition diagnoses and primary
# diagnoses for inpatient hospitalizations, in the primary analysis cohort.
# The top five conditions in each category are used as negative-control
# outcomes in Figures S9 and S10.
#
# Ported from the original NCO_mostpopular10.ipynb: prevalence is the share
# of SUBJECTS ever flagged positive (collapsing each subject's visit-level
# indicators to a single "ever" value via max()), not the share of visit-
# rows flagged positive -- the two differ because a condition diagnosed
# late in follow-up contributes fewer positive rows per subject than one
# diagnosed early, which shifts the ranking at the margin (this is what
# previously made this reconstruction disagree with Figure S9's outcome
# list at rank 5: glaucoma vs. anemia). Condition labels also match the
# original notebook's terminology.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total))

ccw_labels <- c(
  ccw_hthy = "Acquired Hypothyroidism", ccw_myoc = "Acute Myocardial Infarction",
  ccw_alzd = "Alzheimer's Disease and Related Disorders or Senile Dementia",
  ccw_anem = "Anemia", ccw_asth = "Asthma", ccw_afib = "Atrial Fibrillation",
  ccw_bphp = "Benign Prostatic Hyperplasia", ccw_canb = "Cancer, Female/Male Breast",
  ccw_canc = "Cancer, Colorectal", ccw_cane = "Cancer, Endometrial", ccw_canl = "Cancer, Lung",
  ccw_canp = "Cancer, Prostate", ccw_cata = "Cataract", ccw_ckd = "Chronic Kidney Disease",
  ccw_copd = "Chronic Obstructive Pulmonary Disease and Bronchiectasis", ccw_depr = "Depression",
  ccw_diab = "Diabetes", ccw_glau = "Glaucoma", ccw_hfai = "Heart Failure",
  ccw_hipf = "Hip/Pelvic Fracture", ccw_hlip = "Hyperlipidemia", ccw_htn = "Hypertension",
  ccw_ihd = "Ischemic Heart Disease", ccw_oste = "Osteoporosis",
  ccw_raoa = "Rheumatoid Arthritis / Osteoarthritis (RA/OA)", ccw_stro = "Stroke / Transient Ischemic Attack"
)
ip_labels <- c(
  ip_myoca = "Acute Myocardial Infarction", ip_ather = "Coronary Atherosclerosis and Other Heart Disease",
  ip_chest = "Nonspecific Chest Pain", ip_dysrh = "Cardiac Dysrhythmias",
  ip_conhf = "Congestive Heart Failure (Nonhypertensive)", ip_acbvd = "Acute Cerebrovascular Disease",
  ip_pneum = "Pneumonia", ip_copdb = "COPD and Bronchiectasis",
  ip_apneu = "Aspiration Pneumonitis (Food/Vomitus)",
  ip_rfail = "Respiratory Failure / Insufficiency / Arrest (Adult)",
  ip_ghemo = "Gastrointestinal Hemorrhage", ip_arenf = "Acute and Unspecified Renal Failure",
  ip_uti = "Urinary Tract Infections", ip_septi = "Septicemia (Except in Labor)",
  ip_oarth = "Osteoarthritis", ip_hipfx = "Fracture of Neck of Femur (Hip)",
  ip_dcomp = "Complication of Device, Implant, or Graft", ip_synco = "Syncope",
  ip_rehab = "Rehabilitation / Prosthesis Fitting / Device Adjustment",
  ip_fluid = "Fluid and Electrolyte Disorders"
)

ccw_vars <- intersect(names(ccw_labels), names(primary))
ip_vars  <- intersect(names(ip_labels), names(primary))

# Subject-level "ever" prevalence: collapse each subject's visit-level 0/1
# indicator to its max (1 if positive at any visit), then average across
# subjects -- matches subject_ever = df.groupby("subject")[cols].max() in
# the original notebook.
subject_ever <- primary %>%
  group_by(subject) %>%
  summarise(across(all_of(c(ccw_vars, ip_vars)), ~ as.integer(any(.x == 1, na.rm = TRUE))), .groups = "drop")

prevalence_table <- function(vars, labels) {
  tibble::tibble(
    var = vars,
    label = unname(labels[vars]),
    prevalence = vapply(vars, function(v) mean(subject_ever[[v]], na.rm = TRUE), numeric(1)),
    n_cases = vapply(vars, function(v) sum(subject_ever[[v]], na.rm = TRUE), numeric(1))
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
    labs(x = NULL, y = "Share of subjects ever diagnosed/hospitalized", subtitle = subtitle) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size = 11)
}

p_ccw <- plot_prevalence(ccw_prev, "Chronic condition diagnoses")
p_ip  <- plot_prevalence(ip_prev, "Primary diagnoses for inpatient hospitalizations")
p_combined <- patchwork::wrap_plots(p_ccw, p_ip, ncol = 1)
ggsave("output/figures/FigureS8_prevalence_ranking.png", p_combined, width = 9, height = 10, dpi = 200)

cat("\nSaved Figure S8 (prevalence ranking) tables and figure.\n")

# Under this subject-level "ever" prevalence measure, the top 5 chronic
# conditions (hypertension, cataract, RA/OA, hyperlipidemia, anemia) match
# the outcome list used in Figure S9 exactly. The top 5 hospitalization
# categories do not fully match Figure S10's fixed list: congestive heart
# failure ranks 11th here (osteoarthritis, cardiac dysrhythmias, acute
# cerebrovascular disease, pneumonia, and myocardial infarction rank
# higher).
