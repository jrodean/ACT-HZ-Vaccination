# ==============================================================================
# 00c_lasso_risk_score.R
#
# Predicted baseline dementia-risk score (Figure 4 Panel E / HTE variable
# "high_baseline_vulnerability"). Ports the project's existing
# `lasso_predict_dementiarisk.ipynb` (Python/sklearn ElasticNetCV) to R
# (glmnet), matching its target definition, feature list, and baseline
# collapse exactly. The underlying optimizer differs across languages
# (sklearn's coordinate descent + its own CV-fold RNG vs. glmnet's), so the
# selected lambda / exact coefficients will not be bit-identical to the
# original notebook, but the method (pure LASSO, i.e. ElasticNetCV with
# l1_ratio=1 <-> cv.glmnet with alpha=1, 5-fold CV) is the same.
#
# Target (per subject): prior_cog_vulnerability = -1 * (OLS slope of
#   casi_total on days-since-first-pre-vax-visit/365.25), fit on PRE-
#   vaccination visits only, requiring >=3 points.
#
# Output: adds predicted_cog_vulnerability and high_baseline_vulnerability
# (median split) to data/analysis_panel.rds, and writes the final
# data/analysis_panel.dta for the Stata HTE pipeline.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmnet)
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
# 3) hz_beforevax: recorded shingles episode before vaccination (by year)
# ------------------------------------------------------------------------
pre <- pre %>%
  mutate(hz_beforevax = as.integer(!is.na(hz1_date) &
                                    lubridate::year(hz1_date) < lubridate::year(zost1_date)))

# ------------------------------------------------------------------------
# 4) Baseline collapse: first / mean / last pre-vax value per subject
#    (exactly matching the notebook's first_vars/mean_vars/last_vars)
# ------------------------------------------------------------------------
first_vars <- c("b_visit_age", "gender", "white", "education",
                "b_casi_total", "apoe_any4",
                "firstdeg_dem", "firstdeg_ad", "firstdeg_memprob")

# basicneeds_deprived and smoke_intensity are MEAN (not last) in the
# notebook; smoke/alcohol are recoded versions of smoke_bl/alcohol_past_bl,
# not the raw fields -- matches the notebook's df_pre["smoke"],
# df_pre["smoke_intensity"], df_pre["alcohol"], df_pre["basicneeds_deprived"]
mean_vars <- c("income", "basicneeds_deprived", "cesd20_score", "avesys", "avedia", "bmi",
               "adl_sum", "iadl_sum", "gait_time", "grip_strength",
               "smoke_intensity")

ccw_cols <- grep("^ccw_", names(pre), value = TRUE)
ccw_cols <- setdiff(ccw_cols, c("ccw_alzd", "ccw_canp"))  # match notebook's exclusion (alzd, canp)
ip_cols  <- grep("^ip_", names(pre), value = TRUE)

# lb_depre (lookback depression diagnosis) was present in the panel all
# along (person_core keeps every lb_* column) but missing from this feature
# list -- added to match the notebook's last_vars exactly.
last_vars <- c("hypertension", "diabetes", "stroke", "mi", "chf", "sleepapnea",
               "lb_depre", "smoke", "alcohol", "exercise_reg", "hz_beforevax",
               "lupus", "othautoimm", "organbmt", "chemo", "othimmsup",
               ccw_cols, ip_cols)

d <- pre %>%
  mutate(
    basicneeds_deprived = as.integer(basicneeds_met %in% c(2, 3)),
    income = suppressWarnings(as.numeric(income)),
    # recoded versions matching the notebook's df_pre["smoke"]/["smoke_intensity"]/["alcohol"]:
    # smoke = binary (smoke_bl>=1), not raw smoke_bl; smoke_intensity fills
    # non-smokers' missing pack-years with 0 (not left NA-then-median-imputed);
    # alcohol recodes the "don't know/refused" sentinel (9) to missing
    smoke = as.integer(!is.na(smoke_bl) & smoke_bl >= 1),
    smoke_intensity = tidyr::replace_na(suppressWarnings(as.numeric(smoke_pack_bl)), 0),
    alcohol = dplyr::na_if(suppressWarnings(as.numeric(alcohol_past_bl)), 9)
  ) %>%
  arrange(subject, visit_date)

num_or_na <- function(x) suppressWarnings(as.numeric(x))

base_first <- d %>%
  group_by(subject) %>%
  summarise(across(all_of(intersect(first_vars, names(d))), ~ dplyr::first(.x[!is.na(.x)], default = NA)),
            .groups = "drop")

base_mean <- d %>%
  group_by(subject) %>%
  summarise(across(all_of(intersect(mean_vars, names(d))), ~ mean(num_or_na(.x), na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(across(-subject, ~ ifelse(is.nan(.x), NA_real_, .x)))

base_last <- d %>%
  group_by(subject) %>%
  summarise(across(all_of(intersect(last_vars, names(d))),
                    ~ dplyr::last(.x[!is.na(.x)], default = NA)),
            .groups = "drop")

baseline_all <- base_first %>%
  left_join(base_mean, by = "subject") %>%
  left_join(base_last, by = "subject") %>%
  left_join(target, by = "subject") %>%
  mutate(female = as.integer(gender == 2)) %>%
  select(-gender)

# ------------------------------------------------------------------------
# 5) Feature list (matches the notebook's core + immune + hz + ccw + ip)
# ------------------------------------------------------------------------
core_features <- c("b_visit_age", "female", "white", "education",
                    "b_casi_total", "apoe_any4",
                    "firstdeg_dem", "firstdeg_ad", "firstdeg_memprob",
                    "income", "basicneeds_deprived",
                    "lb_depre", "cesd20_score",
                    "hypertension", "diabetes", "stroke", "mi", "chf", "sleepapnea",
                    "avesys", "avedia", "bmi",
                    "smoke", "smoke_intensity", "alcohol", "exercise_reg",
                    "adl_sum", "iadl_sum", "gait_time", "grip_strength")
immune_features <- c("rheumarth", "lupus", "othautoimm", "organbmt", "chemo", "othimmsup")
hz_features <- "hz_beforevax"

feature_list <- intersect(c(core_features, immune_features, hz_features, ccw_cols, ip_cols),
                           names(baseline_all))
# drop all-missing features
miss_rate <- sapply(baseline_all[feature_list], function(x) mean(is.na(x)))
feature_list <- feature_list[miss_rate < 1]

cat(sprintf("N features used: %d\n", length(feature_list)))

# ------------------------------------------------------------------------
# 6) Train/test split (80/20), median imputation (fit on train), LASSO
#    via glmnet (alpha = 1, pure LASSO; cv.glmnet picks lambda by 5-fold CV)
# ------------------------------------------------------------------------
X_all <- baseline_all %>% select(all_of(feature_list)) %>% mutate(across(everything(), as.numeric))
y_all <- baseline_all$prior_cog_vulnerability

obs <- !is.na(y_all)
X_obs <- X_all[obs, ]
y_obs <- y_all[obs]

n_obs <- nrow(X_obs)
set.seed(100)
test_idx <- sample.int(n_obs, size = round(0.2 * n_obs))
X_tr <- X_obs[-test_idx, ]; y_tr <- y_obs[-test_idx]
X_te <- X_obs[test_idx, ];  y_te <- y_obs[test_idx]

fill_vals <- sapply(X_tr, median, na.rm = TRUE)
fill_vals[is.na(fill_vals) | is.nan(fill_vals)] <- 0  # columns all-NA in the train split
impute <- function(df) {
  for (nm in names(df)) df[[nm]][is.na(df[[nm]])] <- fill_vals[nm]
  df
}
X_tr_i <- impute(X_tr); X_te_i <- impute(X_te); X_all_i <- impute(X_all)

tr_mat <- as.matrix(X_tr_i)
center <- colMeans(tr_mat)
scale_ <- apply(tr_mat, 2, sd)
scale_[scale_ == 0 | is.na(scale_)] <- 1  # zero-variance columns: center only, don't divide by 0

X_tr_m  <- scale(tr_mat, center = center, scale = scale_)
X_te_m  <- scale(as.matrix(X_te_i),  center = center, scale = scale_)
X_all_m <- scale(as.matrix(X_all_i), center = center, scale = scale_)

cvfit <- cv.glmnet(X_tr_m, y_tr, alpha = 1, nfolds = 5)
yhat_te <- predict(cvfit, newx = X_te_m, s = "lambda.min")
rmse <- sqrt(mean((y_te - yhat_te)^2))
r2 <- 1 - sum((y_te - yhat_te)^2) / sum((y_te - mean(y_te))^2)
cat(sprintf("TEST RMSE: %.4f\n", rmse))
cat(sprintf("TEST R2: %.4f\n", r2))

coefs <- as.matrix(coef(cvfit, s = "lambda.min"))
nonzero <- coefs[coefs[, 1] != 0, , drop = FALSE]
nonzero <- nonzero[order(-abs(nonzero[, 1])), , drop = FALSE]
cat("\nNon-zero LASSO coefficients (standardized X), largest |coef| first:\n")
print(nonzero)

# ------------------------------------------------------------------------
# Figure S22: non-zero LASSO coefficients (predictors of dementia risk,
# i.e. of prior_cog_vulnerability -- faster decline before vaccination
# loads POSITIVE, so a positive bar means higher predicted dementia risk).
# ------------------------------------------------------------------------
coef_tbl <- tibble::tibble(term = rownames(nonzero), estimate = nonzero[, 1]) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::arrange(estimate)
readr::write_csv(coef_tbl, file.path("output", "tables", "FigureS22_lasso_coefficients_data.csv"))

suppressPackageStartupMessages(library(ggplot2))
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
p22 <- ggplot(coef_tbl, aes(x = estimate, y = reorder(term, estimate), fill = estimate > 0)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"), guide = "none") +
  labs(x = "Standardized LASSO coefficient (higher = higher predicted dementia risk)",
       y = NULL, title = NULL) +
  theme_minimal(base_size = 10)
ggsave(file.path("output", "figures", "FigureS22_lasso_coefficients.png"), p22,
       width = 8, height = max(4, 0.25 * nrow(coef_tbl)), dpi = 200)
cat("\nSaved Figure S22 (LASSO dementia-risk predictors).\n")

# ------------------------------------------------------------------------
# 7) Predict for everyone, add to panel, median split
# ------------------------------------------------------------------------
baseline_all$predicted_cog_vulnerability <- as.numeric(predict(cvfit, newx = X_all_m, s = "lambda.min"))
med <- median(baseline_all$predicted_cog_vulnerability, na.rm = TRUE)
baseline_all <- baseline_all %>%
  mutate(high_baseline_vulnerability = as.integer(predicted_cog_vulnerability >= med))

panel <- panel %>%
  select(-any_of(c("predicted_cog_vulnerability", "high_baseline_vulnerability"))) %>%
  left_join(baseline_all %>% select(subject, predicted_cog_vulnerability, high_baseline_vulnerability),
            by = "subject")

saveRDS(panel, file.path("data", "analysis_panel.rds"))
readr::write_csv(baseline_all %>% select(subject, predicted_cog_vulnerability, high_baseline_vulnerability),
                  file.path("data", "risk_scores.csv"))

cat(sprintf("\nhigh_baseline_vulnerability: %d above-median, %d below-median (of %d with a score)\n",
            sum(baseline_all$high_baseline_vulnerability == 1, na.rm = TRUE),
            sum(baseline_all$high_baseline_vulnerability == 0, na.rm = TRUE),
            sum(!is.na(baseline_all$predicted_cog_vulnerability))))

# ------------------------------------------------------------------------
# 8) Final Stata-ready export: the ONLY .dta this package writes, and the
#    sole input to stata/HTE_all.do (which does every HTE/subgroup panel).
# ------------------------------------------------------------------------
date_cols <- names(panel)[sapply(panel, inherits, "Date")]
panel_dta <- panel %>%
  mutate(across(all_of(date_cols), ~ as.numeric(.x) - as.numeric(as.Date("1960-01-01"))))  # Stata %td origin
write_dta(panel_dta, file.path("data", "analysis_panel.dta"))
cat("Wrote data/analysis_panel.dta for the Stata HTE pipeline.\n")
