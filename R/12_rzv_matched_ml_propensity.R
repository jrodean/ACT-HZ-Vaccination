# ==============================================================================
# 12_rzv_matched_ml_propensity.R
#
# Extends 11_rzv_matched_richer_covariates.R's 6-covariate greedy Euclidean
# match, which still had a pretrend violation and poor balance from naive
# Euclidean distance on a sparse, zero-inflated comorbidity count:
#   (1) richer covariates, including PRIOR VACCINATION BEHAVIOR (flu/
#       pneumococcal/Td uptake before baseline) as a direct proxy for
#       preventive-care-seeking -- the more plausible selection mechanism;
#   (2) propensity scores via random forest (ranger) instead of assuming a
#       linear/Euclidean covariate-treatment relationship;
#   (3) matching on logit(propensity) with a caliper (0.2 SD, Austin 2011),
#       dropping unmatched treated units instead of forcing 100% match;
#   (4) SMD balance, before/after, for every covariate.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ranger) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
event <- read_csv("data/raw/ACT_296_event_table_2025_12_01.csv", show_col_types = FALSE)
names(event) <- tolower(names(event))

baseline_date <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, baseline_date = visit_date)

cci_events <- event %>% filter(stringr::str_starts(event, "cci_")) %>% distinct(event) %>% pull(event)
cci_onset <- event %>% filter(event %in% cci_events) %>%
  distinct(subject, event, days_to_event) %>% inner_join(baseline_date, by = "subject") %>%
  mutate(onset_date = baseline_date + days_to_event) %>%
  group_by(subject, event) %>% summarise(onset_date = min(onset_date), .groups = "drop")
baseline_cci <- cci_onset %>% inner_join(baseline_date, by = "subject") %>%
  filter(onset_date <= baseline_date) %>%
  group_by(subject) %>% summarise(cci_count = n_distinct(event), .groups = "drop")

raw_health_vars <- c("adl_sum","iadl_sum","income","basicneeds_met","diabetes","chf",
                      "hypertension","stroke","sleepapnea","smoke_bl","alcohol_past_bl",
                      "avesys","avedia","bmi","gait_time","grip_strength")

baseline <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, b_visit_age, education, b_casi_total, exercise_reg,
         all_of(raw_health_vars), pneu1_date, td1_date, flu1_date,
         zost1_date, shin1_date) %>%
  left_join(baseline_cci, by = "subject") %>%
  left_join(baseline_date, by = "subject") %>%
  mutate(
    cci_count = tidyr::replace_na(cci_count, 0),
    ever_flu_pre  = as.integer(!is.na(flu1_date)  & flu1_date  <= baseline_date),
    ever_pneu_pre = as.integer(!is.na(pneu1_date) & pneu1_date <= baseline_date),
    ever_td_pre   = as.integer(!is.na(td1_date)   & td1_date   <= baseline_date)
  )

# median-impute sparser health covariates (<10% missing each) with a
# missing-indicator flag, rather than listwise-deleting on 16+ variables
# (cesd20_score, at 46% missing, is dropped entirely -- too sparse to trust)
for (v in raw_health_vars) {
  flagcol <- paste0(v, "_miss")
  baseline[[flagcol]] <- as.integer(is.na(baseline[[v]]))
  med <- stats::median(baseline[[v]], na.rm = TRUE)
  baseline[[v]][is.na(baseline[[v]])] <- med
}

match_vars <- c("b_visit_age","education","b_casi_total","cci_count","exercise_reg",
                 raw_health_vars, paste0(raw_health_vars, "_miss"),
                 "ever_flu_pre","ever_pneu_pre","ever_td_pre")
match_vars <- match_vars[sapply(baseline[match_vars], function(x) length(unique(x[!is.na(x)])) > 1)]
cat(sprintf("Matching on %d covariates: %s\n", length(match_vars), paste(match_vars, collapse=", ")))

smd <- function(df, var, treat_col = "is_treated") {
  t <- df[[var]][df[[treat_col]] == 1]; c <- df[[var]][df[[treat_col]] == 0]
  (mean(t, na.rm=TRUE) - mean(c, na.rm=TRUE)) / sqrt((var(t, na.rm=TRUE) + var(c, na.rm=TRUE)) / 2)
}

match_and_run <- function(vax_date_col, label, seed = 42) {
  core_vars <- c("b_visit_age","education","b_casi_total","exercise_reg")
  cov_all <- baseline %>% filter(if_all(all_of(core_vars), ~ !is.na(.)))
  treated_ids <- cov_all %>% filter(!is.na(.data[[vax_date_col]])) %>% pull(subject)
  never_pool  <- cov_all %>% filter(is.na(zost1_date), is.na(shin1_date)) %>% pull(subject)

  d <- cov_all %>% filter(subject %in% c(treated_ids, never_pool)) %>%
    mutate(is_treated = as.integer(subject %in% treated_ids))
  cat(sprintf("\n=== %s: propensity model, N treated=%d, N pool=%d ===\n", label, length(treated_ids), length(never_pool)))

  rf <- ranger(as.formula(paste("is_treated ~", paste(match_vars, collapse=" + "))),
               data = d %>% mutate(is_treated = factor(is_treated)),
               probability = TRUE, num.trees = 500, seed = seed)
  ps <- predict(rf, data = d)$predictions[, "1"]
  ps <- pmin(pmax(ps, 1e-4), 1 - 1e-4)
  d$logit_ps <- log(ps / (1 - ps))

  caliper <- 0.2 * sd(d$logit_ps)
  zt <- d %>% filter(is_treated == 1)
  zc <- d %>% filter(is_treated == 0)
  set.seed(seed)
  zt <- zt[sample(nrow(zt)), ]
  avail <- zc
  matches <- vector("list", nrow(zt))
  for (i in seq_len(nrow(zt))) {
    if (nrow(avail) == 0) break
    dist <- abs(avail$logit_ps - zt$logit_ps[i])
    j <- which.min(dist)
    if (dist[j] <= caliper) {
      matches[[i]] <- tibble::tibble(treated_subject = zt$subject[i], control_subject = avail$subject[j], dist = dist[j])
      avail <- avail[-j, ]
    }
  }
  match_tbl <- dplyr::bind_rows(matches)
  cat(sprintf("Caliper=%.3f (logit PS). Matched %d / %d treated (%.1f%%)\n",
              caliper, nrow(match_tbl), nrow(zt), 100*nrow(match_tbl)/nrow(zt)))

  matched_long <- d %>% filter(subject %in% c(match_tbl$treated_subject, match_tbl$control_subject))
  cat("\nStandardized mean differences, before vs after matching:\n")
  smd_before <- sapply(match_vars, function(v) smd(d, v))
  smd_after  <- sapply(match_vars, function(v) smd(matched_long, v))
  bal_df <- data.frame(var = match_vars, smd_before = round(smd_before, 3), smd_after = round(smd_after, 3))
  print(bal_df, row.names = FALSE)
  cat(sprintf("Max |SMD| after matching: %.3f  (rule of thumb: <0.1 is good balance)\n", max(abs(smd_after), na.rm=TRUE)))

  pseudo_dates <- baseline %>% filter(subject %in% match_tbl$treated_subject) %>%
    select(treated_subject = subject, pseudo_date = all_of(vax_date_col)) %>%
    inner_join(match_tbl, by = "treated_subject") %>% select(control_subject, pseudo_date)

  treated_d <- panel %>% filter(subject %in% match_tbl$treated_subject, !is.na(casi_total)) %>%
    build_event_time(vax_date_col) %>% mutate(is_treated = 1L)
  control_panel <- panel %>% filter(subject %in% pseudo_dates$control_subject, !is.na(casi_total)) %>%
    left_join(pseudo_dates, by = c("subject" = "control_subject"))
  control_d <- build_event_time(control_panel, "pseudo_date") %>% mutate(is_treated = 0L)
  pooled <- bind_rows(treated_d, control_d)

  inter_terms <- paste0(es_regs, ":is_treated")
  fml <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse="+"), " + ", paste(inter_terms, collapse="+"),
                            " | subject+visit_year+visit_age"))
  m_did <- fixest::feols(fml, data = pooled, vcov = "hetero")

  pre_inter <- intersect(paste0(c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2"), ":is_treated"), names(coef(m_did)))
  b <- coef(m_did)[pre_inter]; V <- vcov(m_did)[pre_inter, pre_inter, drop=FALSE]
  Vinv <- tryCatch(solve(V), error=function(e) MASS::ginv(V))
  W <- as.numeric(t(b) %*% Vinv %*% b)
  pretrend_p <- 1 - pchisq(W, df=length(b))

  ct <- as.data.frame(summary(m_did)$coeftable); ct$term <- rownames(ct)
  cat(sprintf("\nJoint pretrend Wald p (ML-propensity matched DiD): %.4f\n", pretrend_p))
  print(ct[ct$term %in% inter_terms, c("term","Estimate","Std. Error","Pr(>|t|)")], row.names = FALSE)

  list(m_did = m_did, pretrend_p = pretrend_p, match_tbl = match_tbl)
}

cat("################################################################\n")
cat("RZV: ML-propensity matched to never-treated\n")
cat("################################################################\n")
res_rzv <- match_and_run("shin1_date", "RZV")

cat("\n\n################################################################\n")
cat("ZVL: ML-propensity matched to never-treated\n")
cat("################################################################\n")
res_zvl <- match_and_run("zost1_date", "ZVL")

cat(sprintf("\n\n=== SUMMARY ===\nRZV pretrend p: %.4f\nZVL pretrend p: %.4f\n", res_rzv$pretrend_p, res_zvl$pretrend_p))
