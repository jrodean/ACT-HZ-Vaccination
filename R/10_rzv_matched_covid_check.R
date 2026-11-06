# ==============================================================================
# 10_rzv_matched_covid_check.R
#
# In the RZV-only event-study, subject+year+age FEs can only be identified
# from within-RZV-group calendar-year variation. For people vaccinated
# 2018-2019, 2-6-years-after visits fall almost entirely in 2020-2024, so
# "post-treatment" and "COVID era" are nearly collinear and the model
# cannot separate an RZV effect from a general COVID-era trend using RZV
# data alone.
#
# Fix: match each RZV recipient to a never-HZ-vaccinated ACT participant on
# baseline age, education, and CASI Total Score (same covariates as the
# manuscript's matched attrition check), assign the control the same
# calendar date as the RZV partner's real shin1_date (pseudo-vaccination
# date), then re-estimate as a pooled DiD: casi_total ~ es_regs +
# es_regs:is_treated | subject+year+age. The calendar-year FE is now
# identified from treated and matched controls in the same real years, so a
# COVID-era shock common to both is absorbed rather than attributed to RZV.
#
# The original reviewer-response analysis addresses this same concern more
# simply, by restricting to pre-2020 RZV vaccinees / pre-2020 visits rather
# than a matched DiD (a direct COVID-era x event-time term is collinear
# with the calendar-year fixed effect and isn't estimable) -- see R/09's C1
# section. This script is a complementary, more rigorous alternative, not a
# substitute for it.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")

baseline <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, b_visit_age, education, b_casi_total, zost1_date, shin1_date)

match_and_run <- function(vax_date_col, label, seed = 42) {
  treated_ids <- baseline %>%
    filter(!is.na(.data[[vax_date_col]]), !is.na(b_visit_age), !is.na(education), !is.na(b_casi_total)) %>%
    pull(subject)

  never_pool <- baseline %>%
    filter(is.na(zost1_date), is.na(shin1_date), !is.na(b_visit_age), !is.na(education), !is.na(b_casi_total)) %>%
    pull(subject)

  cov <- baseline %>% filter(subject %in% c(treated_ids, never_pool)) %>%
    select(subject, b_visit_age, education, b_casi_total)
  z <- cov %>% mutate(across(c(b_visit_age, education, b_casi_total), ~ as.numeric(scale(.))))

  zt <- z %>% filter(subject %in% treated_ids)
  zc <- z %>% filter(subject %in% never_pool)

  set.seed(seed)
  zt <- zt[sample(nrow(zt)), ]  # randomize match order to avoid systematic order bias
  avail <- zc
  matches <- vector("list", nrow(zt))
  for (i in seq_len(nrow(zt))) {
    if (nrow(avail) == 0) break
    d2 <- (avail$b_visit_age - zt$b_visit_age[i])^2 +
          (avail$education  - zt$education[i])^2 +
          (avail$b_casi_total - zt$b_casi_total[i])^2
    j <- which.min(d2)
    matches[[i]] <- tibble::tibble(treated_subject = zt$subject[i], control_subject = avail$subject[j], dist = sqrt(d2[j]))
    avail <- avail[-j, ]
  }
  match_tbl <- dplyr::bind_rows(matches)

  cat(sprintf("\n=== %s matching ===\n", label))
  cat(sprintf("Treated N=%d, matched N=%d (%.1f%%), never-treated pool N=%d\n",
              nrow(zt), nrow(match_tbl), 100*nrow(match_tbl)/nrow(zt), nrow(zc)))
  cat(sprintf("Mean match distance (standardized): %.3f  (median %.3f)\n", mean(match_tbl$dist), median(match_tbl$dist)))

  pseudo_dates <- baseline %>% filter(subject %in% match_tbl$treated_subject) %>%
    select(treated_subject = subject, pseudo_date = .data[[vax_date_col]]) %>%
    inner_join(match_tbl, by = "treated_subject") %>%
    select(control_subject, pseudo_date)

  treated_d <- panel %>% filter(subject %in% match_tbl$treated_subject, !is.na(casi_total)) %>%
    build_event_time(vax_date_col) %>% mutate(is_treated = 1L)

  control_panel <- panel %>% filter(subject %in% pseudo_dates$control_subject, !is.na(casi_total)) %>%
    left_join(pseudo_dates, by = c("subject" = "control_subject"))
  control_panel$pseudo_vax_date <- control_panel$pseudo_date
  control_d <- build_event_time(control_panel, "pseudo_vax_date") %>% mutate(is_treated = 0L)

  pooled <- bind_rows(treated_d, control_d %>% select(-pseudo_date))

  # (1) naive within-treated-only model (what we had before -- for comparison)
  m_naive <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs, collapse="+"), " | subject+visit_year+visit_age")),
                            data = treated_d, vcov = "hetero")

  # (2) matched DiD: interaction terms are the "true" treatment effect once calendar
  #     trend is pinned down by matched never-treated controls in the same years
  inter_terms <- paste0(es_regs, ":is_treated")
  fml <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse="+"), " + ", paste(inter_terms, collapse="+"),
                            " | subject+visit_year+visit_age"))
  m_did <- fixest::feols(fml, data = pooled, vcov = "hetero")

  cat("\n-- naive within-group (no matched control) --\n")
  print(round(coef(m_naive), 3))
  cat("\n-- matched DiD interaction terms (RZV/ZVL effect net of calendar trend) --\n")
  inter_present <- intersect(inter_terms, names(coef(m_did)))
  print(round(coef(m_did)[inter_present], 3))
  cat("\n-- matched DiD p-values --\n")
  ct <- as.data.frame(summary(m_did)$coeftable); ct$term <- rownames(ct)
  print(ct[ct$term %in% inter_present, c("term","Estimate","Std. Error","Pr(>|t|)")], row.names = FALSE)

  list(m_naive = m_naive, m_did = m_did, match_tbl = match_tbl, pooled = pooled, treated_d = treated_d)
}

cat("################################################################\n")
cat("RZV matched to never-treated controls (pseudo-date = partner's shin1_date)\n")
cat("################################################################\n")
res_rzv <- match_and_run("shin1_date", "RZV")

cat("\n\n################################################################\n")
cat("ZVL matched to never-treated controls (pseudo-date = partner's zost1_date)\n")
cat("################################################################\n")
res_zvl <- match_and_run("zost1_date", "ZVL")

saveRDS(list(rzv = res_rzv, zvl = res_zvl), "output/tables/rzv_zvl_matched_did.rds")
cat("\nSaved output/tables/rzv_zvl_matched_did.rds\n")
