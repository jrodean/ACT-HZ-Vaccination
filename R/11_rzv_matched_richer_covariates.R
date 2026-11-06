# ==============================================================================
# 11_rzv_matched_richer_covariates.R
#
# 10_rzv_matched_covid_check.R's matched DiD design, matching only on
# baseline age, education, and CASI Total Score, shows a pretrend
# violation (future vaccine recipients already diverge from never-treated
# peers years before vaccination) -- that covariate set doesn't identify a
# valid counterfactual.
#
# Enriches matching with overall health burden: Charlson Comorbidity Index
# (CCI) components from the raw event table (cci_* codes), plus baseline
# ADL limitations (frailty) and regular-exercise status, under the
# hypothesis that selection into vaccination operates through general
# health and healthcare engagement, not only cognition/age/education.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
event <- read_csv("data/raw/ACT_296_event_table_2025_12_01.csv", show_col_types = FALSE)
names(event) <- tolower(names(event))

baseline_date <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, baseline_date = visit_date)

cci_events <- event %>% filter(stringr::str_starts(event, "cci_")) %>% distinct(event) %>% pull(event)
cat("CCI components found:", paste(cci_events, collapse=", "), "\n")

cci_onset <- event %>%
  filter(event %in% cci_events) %>%
  distinct(subject, event, days_to_event) %>%
  inner_join(baseline_date, by = "subject") %>%
  mutate(onset_date = baseline_date + days_to_event) %>%
  group_by(subject, event) %>% summarise(onset_date = min(onset_date), .groups = "drop")

baseline_cci <- cci_onset %>%
  inner_join(baseline_date, by = "subject") %>%
  filter(onset_date <= baseline_date) %>%
  group_by(subject) %>% summarise(cci_count = n_distinct(event), .groups = "drop")

baseline <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, b_visit_age, education, b_casi_total, adl_sum, exercise_reg,
         zost1_date, shin1_date) %>%
  left_join(baseline_cci, by = "subject") %>%
  mutate(cci_count = tidyr::replace_na(cci_count, 0))

cat(sprintf("\nN subjects with baseline_cci info: %d / %d (rest assumed 0 comorbidities)\n",
            sum(!is.na(baseline$cci_count)), nrow(baseline)))
cat("CCI count distribution:\n"); print(table(baseline$cci_count))

match_vars <- c("b_visit_age", "education", "b_casi_total", "cci_count", "adl_sum", "exercise_reg")

match_and_run <- function(vax_date_col, label, seed = 42) {
  cov_all <- baseline %>% filter(if_all(all_of(match_vars), ~ !is.na(.)))

  treated_ids <- cov_all %>% filter(!is.na(.data[[vax_date_col]])) %>% pull(subject)
  never_pool  <- cov_all %>% filter(is.na(zost1_date), is.na(shin1_date)) %>% pull(subject)

  z <- cov_all %>% filter(subject %in% c(treated_ids, never_pool)) %>%
    mutate(across(all_of(match_vars), ~ as.numeric(scale(.))))

  zt <- z %>% filter(subject %in% treated_ids)
  zc <- z %>% filter(subject %in% never_pool)

  set.seed(seed)
  zt <- zt[sample(nrow(zt)), ]
  avail <- zc
  matches <- vector("list", nrow(zt))
  for (i in seq_len(nrow(zt))) {
    if (nrow(avail) == 0) break
    d2 <- rowSums(sapply(match_vars, function(v) (avail[[v]] - zt[[v]][i])^2))
    j <- which.min(d2)
    matches[[i]] <- tibble::tibble(treated_subject = zt$subject[i], control_subject = avail$subject[j], dist = sqrt(d2[j]))
    avail <- avail[-j, ]
  }
  match_tbl <- dplyr::bind_rows(matches)

  cat(sprintf("\n=== %s matching (richer covariates: %s) ===\n", label, paste(match_vars, collapse=", ")))
  cat(sprintf("Treated N=%d, matched N=%d (%.1f%%), never-treated pool N=%d\n",
              nrow(zt), nrow(match_tbl), 100*nrow(match_tbl)/nrow(zt), nrow(zc)))
  cat(sprintf("Mean match distance (standardized, %d dims): %.3f  (median %.3f)\n",
              length(match_vars), mean(match_tbl$dist), median(match_tbl$dist)))

  # covariate balance check
  bal <- baseline %>% filter(subject %in% match_tbl$treated_subject) %>% select(all_of(match_vars)) %>%
    summarise(across(everything(), mean, na.rm=TRUE)) %>% mutate(group="treated")
  bal2 <- baseline %>% filter(subject %in% match_tbl$control_subject) %>% select(all_of(match_vars)) %>%
    summarise(across(everything(), mean, na.rm=TRUE)) %>% mutate(group="matched_control")
  cat("\nCovariate balance (means):\n"); print(bind_rows(bal, bal2))

  pseudo_dates <- baseline %>% filter(subject %in% match_tbl$treated_subject) %>%
    select(treated_subject = subject, pseudo_date = all_of(vax_date_col)) %>%
    inner_join(match_tbl, by = "treated_subject") %>%
    select(control_subject, pseudo_date)

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
  cat(sprintf("\nJoint pretrend Wald p (matched DiD, richer covariates): %.4f\n", pretrend_p))
  print(ct[ct$term %in% inter_terms, c("term","Estimate","Std. Error","Pr(>|t|)")], row.names = FALSE)

  list(m_did = m_did, pretrend_p = pretrend_p, match_tbl = match_tbl)
}

cat("################################################################\n")
cat("RZV matched to never-treated, richer covariates\n")
cat("################################################################\n")
res_rzv <- match_and_run("shin1_date", "RZV")

cat("\n\n################################################################\n")
cat("ZVL matched to never-treated, richer covariates\n")
cat("################################################################\n")
res_zvl <- match_and_run("zost1_date", "ZVL")

cat(sprintf("\n\n=== SUMMARY: pretrend p, richer-covariate matched DiD ===\nRZV: %.4f\nZVL: %.4f\n(compare to age/education/CASI-only matching, which showed p<.001 for both)\n",
            res_rzv$pretrend_p, res_zvl$pretrend_p))
