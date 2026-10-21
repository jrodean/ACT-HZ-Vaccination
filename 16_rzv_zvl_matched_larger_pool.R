# ==============================================================================
# 16_rzv_zvl_matched_larger_pool.R
#
# The matched design in 15_rzv_zvl_direct_matched.R, comparing RZV-only
# (2018-2019) against the ZVL-only pool (N=850), has poor enrollment-year
# balance (SMD=1.70, mean gap 6.9 years): there are not enough late-
# enrolling, ZVL-exclusive people, since most people who enrolled late and
# wanted HZ vaccination after roughly 2016 received RZV once it existed.
#
# This script expands the ZVL candidate pool to all ZVL recipients
# (N=1,531, including those who later also received RZV), which has many
# more late enrollees. A matched dual-vaccinated ZVL candidate's follow-up
# is censored at their own subsequent RZV dose (if any), so their post-ZVL
# coefficients are not contaminated by a later RZV dose.
#
# Matching variables: enrollment year (stratified first, as in script 15),
# enrollment age (age at first ACT visit) and vaccination age (age at the
# actual ZVL/RZV dose) as explicit, separate dimensions, plus education and
# baseline CASI.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
baseline_date <- panel %>% filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, baseline_date = visit_date)

baseline <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, enroll_year = visit_year, b_visit_age, education, b_casi_total,
         zost1_date, shin1_date, shin2_date) %>%
  left_join(baseline_date, by = "subject") %>%
  mutate(
    age_at_zost1 = b_visit_age + as.numeric(zost1_date - baseline_date) / 365.25,
    age_at_shin1 = b_visit_age + as.numeric(shin1_date - baseline_date) / 365.25
  )

rzv_only_1819_ids <- baseline %>% filter(!is.na(shin1_date), is.na(zost1_date),
                                          as.numeric(format(shin1_date, "%Y")) %in% 2018:2019) %>% pull(subject)
zvl_pool_ids <- baseline %>% filter(!is.na(zost1_date), !is.na(age_at_zost1)) %>% pull(subject)
cat(sprintf("ZVL candidate pool (ALL ZVL recipients, any RZV status): N=%d\n", length(zvl_pool_ids)))
cat(sprintf("RZV-only (2018-2019) targets: N=%d\n", length(rzv_only_1819_ids)))

cov <- baseline %>% filter(subject %in% c(zvl_pool_ids, rzv_only_1819_ids),
                            !is.na(b_visit_age), !is.na(education), !is.na(b_casi_total))
cov$vax_age <- ifelse(cov$subject %in% rzv_only_1819_ids, cov$age_at_shin1, cov$age_at_zost1)

z <- cov %>% mutate(across(c(b_visit_age, vax_age, education, b_casi_total), ~ as.numeric(scale(.)))) %>%
  rename_with(~ paste0(., "_z"), c(b_visit_age, vax_age, education, b_casi_total))

treated <- z %>% filter(subject %in% rzv_only_1819_ids)
pool    <- z %>% filter(subject %in% zvl_pool_ids)

set.seed(42)
treated <- treated[sample(nrow(treated)), ]
avail <- pool
matches <- vector("list", nrow(treated))
for (i in seq_len(nrow(treated))) {
  if (nrow(avail) == 0) break
  w <- 0
  repeat {
    cand <- avail %>% filter(abs(enroll_year - treated$enroll_year[i]) <= w)
    if (nrow(cand) > 0 || w >= 15) break
    w <- w + 1
  }
  if (nrow(cand) == 0) next
  d2 <- (cand$b_visit_age_z - treated$b_visit_age_z[i])^2 +
        (cand$vax_age_z     - treated$vax_age_z[i])^2 +
        (cand$education_z   - treated$education_z[i])^2 +
        (cand$b_casi_total_z - treated$b_casi_total_z[i])^2
  j <- which.min(d2)
  matches[[i]] <- tibble::tibble(rzv_subject = treated$subject[i], zvl_subject = cand$subject[j], enroll_gap = w)
  avail <- avail %>% filter(subject != cand$subject[j])
}
match_tbl <- dplyr::bind_rows(matches)
cat(sprintf("\nMatched %d / %d RZV-only(18-19) recipients (%.1f%%)\n",
            nrow(match_tbl), nrow(treated), 100*nrow(match_tbl)/nrow(treated)))
cat("Enrollment-year gap distribution:\n"); print(table(match_tbl$enroll_gap))
cat(sprintf("Mean |enrollment year gap|: %.2f  (vs. 6.90 with the ZVL-only-restricted pool)\n", mean(match_tbl$enroll_gap)))

n_dual_matched <- sum(match_tbl$zvl_subject %in% (baseline %>% filter(!is.na(shin1_date)) %>% pull(subject)))
cat(sprintf("Matched ZVL partners who are ALSO dual-vaccinated (will be censored at their own RZV date): %d / %d (%.1f%%)\n",
            n_dual_matched, nrow(match_tbl), 100*n_dual_matched/nrow(match_tbl)))

bal <- cov %>% filter(subject %in% c(match_tbl$rzv_subject, match_tbl$zvl_subject)) %>%
  mutate(grp = if_else(subject %in% match_tbl$rzv_subject, "RZV", "ZVL"))
smd <- function(v) {
  t <- bal[[v]][bal$grp=="RZV"]; c <- bal[[v]][bal$grp=="ZVL"]
  (mean(t,na.rm=TRUE)-mean(c,na.rm=TRUE)) / sqrt((var(t,na.rm=TRUE)+var(c,na.rm=TRUE))/2)
}
cat("\nCovariate balance (SMD) after matching:\n")
for (v in c("enroll_year","b_visit_age","vax_age","education","b_casi_total")) {
  cat(sprintf("%-14s SMD=%.3f  (RZV mean=%.2f, ZVL mean=%.2f)\n", v, smd(v),
              mean(bal[[v]][bal$grp=="RZV"],na.rm=TRUE), mean(bal[[v]][bal$grp=="ZVL"],na.rm=TRUE)))
}

# build event-time data, CENSORING any matched ZVL partner's follow-up at
# their own subsequent RZV dose (if they have one) to avoid contamination
# (panel already has shin1_date for every row, no join needed)
zvl_d <- panel %>% filter(subject %in% match_tbl$zvl_subject, !is.na(casi_total),
                           is.na(shin1_date) | visit_date < shin1_date) %>%
  build_event_time("zost1_date")
rzv_d <- panel %>% filter(subject %in% match_tbl$rzv_subject, !is.na(casi_total)) %>% build_event_time("shin1_date")

res_zvl <- run_event_study(zvl_d, "casi_total")
res_rzv <- run_event_study(rzv_d, "casi_total")
cat(sprintf("\nZVL (matched, censored at own RZV) N=%d rows=%d, pretrend p=%.4f\n", n_distinct(zvl_d$subject), nrow(zvl_d), res_zvl$pretrend_p))
cat(sprintf("RZV (matched) N=%d rows=%d, pretrend p=%.4f\n", n_distinct(rzv_d$subject), nrow(rzv_d), res_rzv$pretrend_p))
cat("ZVL coefs:\n"); print(round(coef(res_zvl$mod),3))
cat("RZV coefs:\n"); print(round(coef(res_rzv$mod),3))

zvl_d$is_rzv <- 0; rzv_d$is_rzv <- 1
pooled <- bind_rows(zvl_d, rzv_d)
inter_terms <- paste0(es_regs, ":is_rzv")
fml <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse=" + "), " + ", paste(inter_terms, collapse=" + "), " | subject + visit_year + visit_age"))
m <- fixest::feols(fml, data = pooled, vcov = "hetero", notes = FALSE)
post_terms <- paste0("trt_p",0:6); inter_post <- intersect(paste0(post_terms,":is_rzv"), names(coef(m)))
b <- coef(m)[inter_post]; V <- vcov(m)[inter_post,inter_post,drop=FALSE]
Vinv <- tryCatch(solve(V), error=function(e) MASS::ginv(V))
W <- as.numeric(t(b) %*% Vinv %*% b); p_diff <- 1-pchisq(W, df=length(b))
cat(sprintf("\nJoint Wald test, RZV vs ZVL (larger pool, censored, matched on enroll year+age+vax age): p=%.4f (df=%d)\n", p_diff, length(b)))
ct <- as.data.frame(summary(m)$coeftable); ct$term <- rownames(ct)
print(ct[ct$term %in% inter_post, c("term","Estimate","Std. Error","Pr(>|t|)")], row.names = FALSE)

p <- es_plot_two(res_rzv$mod, res_zvl$mod, "RZV (2018-2019)", "ZVL (matched, censored at own RZV)",
                  title = NULL,
                  ytitle = "Effect on CASI Total Score")
ggsave("output/figures/RZV_ZVL_matched_larger_pool.png", p, width = 10, height = 6, dpi = 220)
cat("saved plot\n")
