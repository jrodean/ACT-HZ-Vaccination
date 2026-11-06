# ==============================================================================
# 13_sequential_matching_search.R
#
# Searches over (covariate set x matching method) combinations for a
# never-treated comparison group with genuinely parallel pre-trends. A
# pre-trend p > 0.05 alone is not sufficient (a small, badly-balanced match
# can pass just from losing power), so a step counts as a genuine pass only
# if all three hold: pretrend_p > 0.10, max |SMD| after matching < 0.10, and
# match_rate > 0.70 (not a tiny, cherry-picked remnant).
#
# Covariates are added in a fixed, theory-driven order ("healthy vaccinee"
# mechanisms first: preventive-care-seeking, then comorbidity/frailty, then
# individual conditions, then lifestyle, then SES), crossed with 4 matching
# methods per step: euclidean/mahalanobis NN (no caliper), and ps_logit/
# ps_rf propensity-score NN (caliper 0.2 SD(logit), rf via ranger).
#
# Every attempt is logged, not just the first pass, so the output is a
# transparent record of the whole search rather than one cherry-picked run.
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
  filter(onset_date <= baseline_date) %>% group_by(subject) %>%
  summarise(cci_count = n_distinct(event), .groups = "drop")

raw_health_vars <- c("adl_sum","iadl_sum","income","basicneeds_met","diabetes","chf",
                      "hypertension","stroke","sleepapnea","smoke_bl","alcohol_past_bl",
                      "avesys","avedia","bmi","gait_time","grip_strength")

baseline <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, b_visit_age, education, b_casi_total, exercise_reg,
         all_of(raw_health_vars), pneu1_date, td1_date, flu1_date, zost1_date, shin1_date) %>%
  left_join(baseline_cci, by = "subject") %>%
  left_join(baseline_date, by = "subject") %>%
  mutate(
    cci_count = tidyr::replace_na(cci_count, 0),
    ever_flu_pre  = as.integer(!is.na(flu1_date)  & flu1_date  <= baseline_date),
    ever_pneu_pre = as.integer(!is.na(pneu1_date) & pneu1_date <= baseline_date),
    ever_td_pre   = as.integer(!is.na(td1_date)   & td1_date   <= baseline_date)
  )
for (v in raw_health_vars) {
  baseline[[paste0(v, "_miss")]] <- as.integer(is.na(baseline[[v]]))
  baseline[[v]][is.na(baseline[[v]])] <- stats::median(baseline[[v]], na.rm = TRUE)
}

core_vars <- c("b_visit_age", "education", "b_casi_total")
# theory-driven order: preventive-care-seeking proxies first, then general
# health burden/frailty, then individual conditions, then lifestyle, then SES
candidate_order <- c("ever_flu_pre", "ever_pneu_pre", "ever_td_pre",
                      "cci_count", "adl_sum", "iadl_sum", "exercise_reg",
                      "diabetes", "chf", "hypertension", "stroke", "sleepapnea",
                      "smoke_bl", "alcohol_past_bl", "avesys", "avedia",
                      "bmi", "gait_time", "grip_strength",
                      "income", "basicneeds_met")

smd <- function(df, var, treat_col = "is_treated") {
  t <- df[[var]][df[[treat_col]] == 1]; c <- df[[var]][df[[treat_col]] == 0]
  sdp <- sqrt((stats::var(t, na.rm=TRUE) + stats::var(c, na.rm=TRUE)) / 2)
  if (is.na(sdp) || sdp == 0) return(NA_real_)
  (mean(t, na.rm=TRUE) - mean(c, na.rm=TRUE)) / sdp
}

run_match <- function(d, match_vars, method, seed = 42) {
  zt <- d %>% filter(is_treated == 1)
  zc <- d %>% filter(is_treated == 0)

  if (method == "euclidean") {
    z <- d %>% mutate(across(all_of(match_vars), ~ as.numeric(scale(.))))
    zt_s <- z %>% filter(is_treated == 1); zc_s <- z %>% filter(is_treated == 0)
    dist_fn <- function(i, avail) rowSums(sapply(match_vars, function(v) (avail[[v]] - zt_s[[v]][i])^2))
    caliper <- Inf
  } else if (method == "mahalanobis") {
    X <- as.matrix(d[, match_vars]); S <- stats::cov(X)
    Sinv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
    zt_s <- zt; zc_s <- zc
    dist_fn <- function(i, avail) {
      Xt <- as.numeric(zt[match_vars][i, ]); Xa <- as.matrix(avail[match_vars])
      diffs <- sweep(Xa, 2, Xt, "-")
      rowSums((diffs %*% Sinv) * diffs)
    }
    caliper <- Inf
  } else if (method %in% c("ps_logit", "ps_rf")) {
    if (method == "ps_logit") {
      m <- stats::glm(as.formula(paste("is_treated ~", paste(match_vars, collapse=" + "))),
                       data = d, family = binomial())
      ps <- stats::predict(m, type = "response")
    } else {
      rf <- ranger(as.formula(paste("is_treated ~", paste(match_vars, collapse=" + "))),
                   data = d %>% mutate(is_treated = factor(is_treated)),
                   probability = TRUE, num.trees = 300, seed = seed)
      ps <- stats::predict(rf, data = d)$predictions[, "1"]
    }
    ps <- pmin(pmax(ps, 1e-4), 1 - 1e-4)
    d$logit_ps <- log(ps / (1 - ps))
    zt_s <- d %>% filter(is_treated == 1); zc_s <- d %>% filter(is_treated == 0)
    caliper <- 0.2 * stats::sd(d$logit_ps)
    dist_fn <- function(i, avail) abs(avail$logit_ps - zt_s$logit_ps[i])
  }

  set.seed(seed)
  ord <- sample(nrow(zt_s))
  avail <- zc_s
  matches <- vector("list", nrow(zt_s))
  for (idx in seq_along(ord)) {
    i <- ord[idx]
    if (nrow(avail) == 0) break
    dd <- dist_fn(i, avail)
    j <- which.min(dd)
    if (dd[j] <= caliper) {
      matches[[idx]] <- tibble::tibble(treated_subject = zt_s$subject[i], control_subject = avail$subject[j])
      avail <- avail[-j, ]
    }
  }
  dplyr::bind_rows(matches)
}

pretrend_test <- function(panel, match_tbl, vax_date_col) {
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
  m <- tryCatch(fixest::feols(fml, data = pooled, vcov = "hetero", notes = FALSE), error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  pre_inter <- intersect(paste0(c("trt_m6","trt_m5","trt_m4","trt_m3","trt_m2"), ":is_treated"), names(coef(m)))
  if (length(pre_inter) < 2) return(NA_real_)
  b <- coef(m)[pre_inter]; V <- vcov(m)[pre_inter, pre_inter, drop=FALSE]
  Vinv <- tryCatch(solve(V), error=function(e) MASS::ginv(V))
  W <- as.numeric(t(b) %*% Vinv %*% b)
  1 - pchisq(W, df=length(b))
}

search <- function(vax_date_col, label) {
  cov_all <- baseline %>% filter(if_all(all_of(core_vars), ~ !is.na(.)))
  treated_ids <- cov_all %>% filter(!is.na(.data[[vax_date_col]])) %>% pull(subject)
  never_pool  <- cov_all %>% filter(is.na(zost1_date), is.na(shin1_date)) %>% pull(subject)
  d_full <- cov_all %>% filter(subject %in% c(treated_ids, never_pool)) %>%
    mutate(is_treated = as.integer(subject %in% treated_ids))

  methods <- c("euclidean", "mahalanobis", "ps_logit", "ps_rf")
  log <- list()
  active_vars <- core_vars
  found_pass <- FALSE

  for (step in 0:length(candidate_order)) {
    if (step > 0) active_vars <- c(core_vars, candidate_order[1:step])
    d <- d_full
    for (method in methods) {
      match_tbl <- tryCatch(run_match(d, active_vars, method), error = function(e) tibble::tibble())
      if (nrow(match_tbl) == 0) next
      match_rate <- nrow(match_tbl) / length(treated_ids)
      matched_long <- d %>% filter(subject %in% c(match_tbl$treated_subject, match_tbl$control_subject))
      smds <- sapply(active_vars, function(v) smd(matched_long, v))
      max_smd <- max(abs(smds), na.rm = TRUE)
      p <- pretrend_test(panel, match_tbl, vax_date_col)

      passed <- !is.na(p) && p > 0.10 && !is.na(max_smd) && max_smd < 0.10 && match_rate > 0.70
      log[[length(log) + 1]] <- tibble::tibble(
        group = label, step = step, n_vars = length(active_vars), method = method,
        last_var_added = if (step == 0) "core only" else candidate_order[step],
        n_matched = nrow(match_tbl), match_rate = round(match_rate, 3),
        max_smd = round(max_smd, 3), pretrend_p = round(p, 4), passed = passed
      )
      cat(sprintf("[%s] step=%2d (+%-16s) method=%-11s N=%4d match_rate=%.2f max|SMD|=%.3f pretrend_p=%.4f %s\n",
                  label, step, if (step==0) "-" else candidate_order[step], method,
                  nrow(match_tbl), match_rate, max_smd, p, if (passed) "<-- PASS (all 3 criteria)" else ""))
      if (passed) { found_pass <- TRUE; break }
    }
    if (found_pass) break
  }
  dplyr::bind_rows(log)
}

cat("################################################################\n")
cat("Sequential search: RZV vs never-treated\n")
cat("################################################################\n")
log_rzv <- search("shin1_date", "RZV")

cat("\n\n################################################################\n")
cat("Sequential search: ZVL vs never-treated\n")
cat("################################################################\n")
log_zvl <- search("zost1_date", "ZVL")

full_log <- bind_rows(log_rzv, log_zvl)
readr::write_csv(full_log, "output/tables/sequential_matching_search_log.csv")
cat("\n\n=== Any passes? ===\n")
print(full_log %>% filter(passed))
cat(sprintf("\nTotal specifications tried: %d\n", nrow(full_log)))
cat("Full log saved to output/tables/sequential_matching_search_log.csv\n")
