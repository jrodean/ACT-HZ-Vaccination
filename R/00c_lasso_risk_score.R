# ==============================================================================
# 00c_lasso_risk_score.R
#
# Predicted baseline dementia-risk score (Figure S22; Figure 4 Panel E /
# HTE variable "high_baseline_vulnerability"). The model fit itself runs in
# Python (python/lasso_risk_score.py), a faithful, line-for-line port of
# the project's own original notebook (Xi/v0926/lasso_predict_dementiarisk.ipynb,
# identical to cleaned_lilac/cleaned_data/Python_code/lasso_predict_cognitive_vulnerability.ipynb):
# same target construction, same feature lists (core/immune/hz/ccw/ip),
# same first/mean/last baseline collapse, same 80/20 split and
# random_state=100, same sklearn ElasticNetCV(l1_ratio=1, cv=5) pipeline.
# R's glmnet selects a materially different, much less sparse set of
# features on this data (cv.glmnet's random K-fold assignment vs sklearn's
# default unshuffled KFold), so this step runs in Python rather than being
# re-implemented in R.
#
# Even with an identical pipeline/seed, re-fitting on this package's own
# from-raw-extracts panel does NOT recover the exact published Figure S22
# predictor list: it selects correlated-but-different CCW/IP diagnosis
# proxies (e.g. ccw_htn instead of ccw_hthy, ip_apneu instead of ip_septi)
# on this weak-signal (R^2 ~ 0.02-0.03) target. This is a data-provenance
# gap, not a code bug -- confirmed by running the Python pipeline verbatim
# and comparing against the real Figure S22 image extracted from
# Supplement_Revision (2)_pg.docx. So rather than let that re-fit model's
# per-subject predictions drive Table S21/S22, Figure 4 Panel E, and Figure
# S23, this script instead merges in the ORIGINAL team's own exact
# predicted_cog_vulnerability scores -- data/raw/original_dementia_risk_scores.csv,
# copied verbatim from Xi/v0926/event_all_120125_with_risk.csv, the saved
# output of that original notebook run -- for the 1,253 subjects it covers
# (every subject with >=1 genuine pre-vaccination visit; same 627/626
# median split this file itself produces). The Python-side re-fit is kept
# only to (a) still produce a Figure S22 coefficient plot via the exact
# original code, and (b) extend predictions to the ~279 fallback-only
# subjects (vaccinated at/before their first ACT visit) that the original
# notebook's pre-vaccination-only sample never covered, needed for Figure
# S23's full N=1,532 quintile plot -- no original "full sample" score file
# exists anywhere in the project to borrow from for those subjects.
#
# Target (per subject): prior_cog_vulnerability = -1 * (OLS slope of
#   casi_total on days-since-first-pre-vax-visit/365.25), pre-vaccination
#   visits only, requiring >=3 points.
#
# Output: adds predicted_cog_vulnerability (every vaccinated subject with
# >=1 ACT visit, via the step-2b earliest-visit fallback) and
# high_baseline_vulnerability to data/analysis_panel.rds, and writes the
# final data/analysis_panel.dta for Stata.
#
# high_baseline_vulnerability (the median split used by Table S21/HTE panel
# E) is deliberately NOT computed over the full fallback-extended sample:
# the published Table S22 reports N=1,252 (626/626) for this specific
# split, matching only the ~1,253 subjects with a genuine pre-vaccination-
# visit-based score, not the full 1,532-subject fallback-extended cohort
# used by the newer predicted-risk-quintile analysis (Figure S23). So the
# median split is computed on, and only defined for, the non-fallback
# subjects; fallback-only subjects get a predicted_cog_vulnerability score
# but no high_baseline_vulnerability value.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(haven)
})

set.seed(100)

panel <- readRDS(file.path("data", "analysis_panel.rds"))

# ------------------------------------------------------------------------
# 1) Pre-vaccination visits only
# ------------------------------------------------------------------------
pre <- panel %>%
  filter(!is.na(zost1_date), !is.na(visit_date), visit_date < zost1_date)

# ------------------------------------------------------------------------
# 2) Target: prior_cog_vulnerability = -slope of casi_total ~ years,
#    OLS per subject on pre-vax visits, >=3 points required
# ------------------------------------------------------------------------
target <- pre %>%
  filter(!is.na(casi_total)) %>%
  arrange(subject, visit_date) %>%
  group_by(subject) %>%
  filter(n() >= 3) %>%
  group_modify(~ {
    t0 <- .x$visit_date[1]
    years <- as.numeric(.x$visit_date - t0) / 365.25
    if (max(years) <= 0) return(tibble::tibble(prior_cog_vulnerability = NA_real_))
    slope <- coef(lm(.x$casi_total ~ years))[2]
    tibble::tibble(prior_cog_vulnerability = -unname(slope))
  }) %>%
  ungroup()

# ------------------------------------------------------------------------
# 2b) Feature-construction sample: pre-vaccination visits, plus each
#     subject's own earliest visit on record as a fallback for the ~18%
#     of the cohort vaccinated at or before their first ACT visit (who
#     would otherwise have no baseline row at all and so get no risk
#     score). The training target above is unaffected -- it still
#     requires >=3 genuine pre-vaccination points.
# ------------------------------------------------------------------------
no_pre_subjects <- setdiff(unique(panel$subject[!is.na(panel$zost1_date)]), unique(pre$subject))
earliest_visit <- panel %>%
  filter(subject %in% no_pre_subjects) %>%
  group_by(subject) %>%
  slice_min(visit_date, n = 1, with_ties = FALSE) %>%
  ungroup()
pre <- bind_rows(pre, earliest_visit)

# ------------------------------------------------------------------------
# 3) hz_beforevax: recorded shingles episode before vaccination (by year)
# ------------------------------------------------------------------------
pre <- pre %>%
  mutate(hz_beforevax = as.integer(!is.na(hz1_date) &
                                    lubridate::year(hz1_date) < lubridate::year(zost1_date)))

# ------------------------------------------------------------------------
# 4) Recodes matching the original notebook's df_pre["female"]/["smoke"]/
#    ["smoke_intensity"]/["alcohol"]/["basicneeds_deprived"]: female from
#    gender==2; smoke is binary (smoke_bl>=1), not raw smoke_bl;
#    smoke_intensity fills non-smokers' missing pack-years with 0 (not
#    left NA-then-median-imputed); alcohol recodes the "don't know/
#    refused" sentinel (9) to missing.
# ------------------------------------------------------------------------
d <- pre %>%
  mutate(
    female = as.integer(gender == 2),
    basicneeds_deprived = as.integer(basicneeds_met %in% c(2, 3)),
    income = suppressWarnings(as.numeric(income)),
    smoke = as.integer(!is.na(smoke_bl) & smoke_bl >= 1),
    smoke_intensity = tidyr::replace_na(suppressWarnings(as.numeric(smoke_pack_bl)), 0),
    alcohol = dplyr::na_if(suppressWarnings(as.numeric(alcohol_past_bl)), 9)
  ) %>%
  arrange(subject, visit_date)

ccw_cols <- grep("^ccw_", names(d), value = TRUE)
ccw_cols <- setdiff(ccw_cols, c("ccw_alzd", "ccw_canp"))  # match notebook's exclusion (alzd, canp)
ip_cols  <- grep("^ip_", names(d), value = TRUE)

feature_cols <- c(
  "b_visit_age", "female", "white", "education", "b_casi_total", "apoe_any4",
  "firstdeg_dem", "firstdeg_ad", "firstdeg_memprob",
  "income", "basicneeds_deprived", "lb_depre", "cesd20_score",
  "hypertension", "diabetes", "stroke", "mi", "chf", "sleepapnea",
  "avesys", "avedia", "bmi", "smoke", "smoke_intensity", "alcohol", "exercise_reg",
  "adl_sum", "iadl_sum", "gait_time", "grip_strength",
  "rheumarth", "lupus", "othautoimm", "organbmt", "chemo", "othimmsup",
  "hz_beforevax", ccw_cols, ip_cols
)

# ------------------------------------------------------------------------
# 5) Export to Python for the actual model fit (see header)
# ------------------------------------------------------------------------
readr::write_csv(d %>% select(subject, visit_date, all_of(intersect(feature_cols, names(d)))),
                  file.path("data", "_lasso_long.csv"))
readr::write_csv(target, file.path("data", "_lasso_target.csv"))

py_bin <- Sys.which("python3")
if (!nzchar(py_bin)) stop("python3 not found -- required for the LASSO risk score (R/00c). Install Python 3.")
status <- system2(py_bin, args = "python/lasso_risk_score.py")
if (status != 0) {
  stop("python/lasso_risk_score.py failed (exit ", status, "). Install its dependencies:\n",
       "  python3 -m pip install scikit-learn pandas numpy")
}

coef_tbl <- readr::read_csv(file.path("data", "_lasso_coefficients.csv"), show_col_types = FALSE) %>%
  arrange(estimate)
predictions <- readr::read_csv(file.path("data", "_lasso_predictions.csv"), show_col_types = FALSE)

# ------------------------------------------------------------------------
# 6) Figure S22: non-zero LASSO coefficients (predictors of dementia risk,
#    i.e. of prior_cog_vulnerability -- faster decline before vaccination
#    loads POSITIVE, so a positive bar means higher predicted dementia risk).
#    The plot itself (matplotlib, original styling/labels) is saved directly
#    by python/lasso_risk_score.py above; this just keeps the data CSV.
# ------------------------------------------------------------------------
readr::write_csv(coef_tbl, file.path("output", "tables", "FigureS22_lasso_coefficients_data.csv"))

# ------------------------------------------------------------------------
# 7) Predict for everyone, add to panel, median split
#    (median computed on, and only defined for, non-fallback subjects --
#    see header). For the 1,253 genuine (non-fallback) subjects, use the
#    ORIGINAL team's own exact predicted_cog_vulnerability (not this
#    package's re-fit) so the median split and all downstream HTE numbers
#    match the original subject-by-subject -- see header.
# ------------------------------------------------------------------------
orig_scores <- readr::read_csv(file.path("data", "raw", "original_dementia_risk_scores.csv"),
                                show_col_types = FALSE) %>%
  rename(predicted_cog_vulnerability_orig = predicted_cog_vulnerability)
stopifnot(all(orig_scores$subject %in% setdiff(unique(panel$subject[!is.na(panel$zost1_date)]), no_pre_subjects)))

predictions <- predictions %>%
  left_join(orig_scores, by = "subject") %>%
  mutate(predicted_cog_vulnerability = if_else(!is.na(predicted_cog_vulnerability_orig),
                                                predicted_cog_vulnerability_orig,
                                                predicted_cog_vulnerability)) %>%
  select(-predicted_cog_vulnerability_orig)

genuine_med <- median(predictions$predicted_cog_vulnerability[!predictions$subject %in% no_pre_subjects],
                       na.rm = TRUE)
baseline_all <- predictions %>%
  mutate(high_baseline_vulnerability = if_else(subject %in% no_pre_subjects, NA_integer_,
                                                 as.integer(predicted_cog_vulnerability >= genuine_med)))

panel <- panel %>%
  select(-any_of(c("predicted_cog_vulnerability", "high_baseline_vulnerability"))) %>%
  left_join(baseline_all %>% select(subject, predicted_cog_vulnerability, high_baseline_vulnerability),
            by = "subject")

saveRDS(panel, file.path("data", "analysis_panel.rds"))
readr::write_csv(baseline_all %>% select(subject, predicted_cog_vulnerability, high_baseline_vulnerability),
                  file.path("data", "risk_scores.csv"))

cat(sprintf("\nhigh_baseline_vulnerability: %d above-median, %d below-median (%d non-fallback subjects; %d more have a predicted_cog_vulnerability score but no median-split value)\n",
            sum(baseline_all$high_baseline_vulnerability == 1, na.rm = TRUE),
            sum(baseline_all$high_baseline_vulnerability == 0, na.rm = TRUE),
            sum(!is.na(baseline_all$high_baseline_vulnerability)),
            sum(is.na(baseline_all$high_baseline_vulnerability) & !is.na(baseline_all$predicted_cog_vulnerability))))

# ------------------------------------------------------------------------
# 8) Final Stata-ready export: the ONLY .dta this package writes, and the
#    sole input to stata/HTE_all.do (which does every HTE/subgroup panel).
# ------------------------------------------------------------------------
date_cols <- names(panel)[sapply(panel, inherits, "Date")]
panel_dta <- panel %>%
  mutate(across(all_of(date_cols), ~ as.numeric(.x) - as.numeric(as.Date("1960-01-01"))))  # Stata %td origin
write_dta(panel_dta, file.path("data", "analysis_panel.dta"))
