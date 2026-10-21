# ==============================================================================
# 15_rzv_zvl_direct_matched.R
#
# The vaccine-exclusive, calendar-restricted comparison (ZVL-only 2014-2017,
# N=139 vs RZV-only 2018-2019, N=418) has a small, noisy ZVL side because
# it is filtered to a 4-year window. This script instead directly matches
# each RZV-only (2018-2019) recipient to the best available ZVL-only
# recipient from the full ZVL-only pool (2006-2017, N=850), stratifying on
# enrollment year first (as in 14_enrollment_matched_pretrend.R), then
# refining on age, education, and baseline CASI — a treated-vs-treated
# match rather than treated-vs-never-treated, retaining more of the
# ZVL-only pool's useful variation than calendar filtering alone.
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
  left_join(baseline_date, by = "subject")

zvl_only_ids <- baseline %>% filter(!is.na(zost1_date), is.na(shin1_date)) %>% pull(subject)
rzv_only_1819_ids <- baseline %>% filter(!is.na(shin1_date), is.na(zost1_date),
                                          as.numeric(format(shin1_date, "%Y")) %in% 2018:2019) %>% pull(subject)
cat(sprintf("ZVL-only pool (any year, 2006-2017): N=%d\n", length(zvl_only_ids)))
cat(sprintf("RZV-only (2018-2019): N=%d\n", length(rzv_only_1819_ids)))

cov <- baseline %>% filter(subject %in% c(zvl_only_ids, rzv_only_1819_ids), !is.na(b_visit_age), !is.na(education), !is.na(b_casi_total))
treated <- cov %>% filter(subject %in% rzv_only_1819_ids)   # RZV is the "treated" we're finding matches FOR
pool    <- cov %>% filter(subject %in% zvl_only_ids)

z <- cov %>% mutate(across(c(b_visit_age, education, b_casi_total), ~ as.numeric(scale(.)))) %>%
  rename_with(~ paste0(., "_z"), c(b_visit_age, education, b_casi_total))
treated <- treated %>% left_join(z %>% select(subject, ends_with("_z")), by = "subject")
pool    <- pool    %>% left_join(z %>% select(subject, ends_with("_z")), by = "subject")

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
        (cand$education_z  - treated$education_z[i])^2 +
        (cand$b_casi_total_z - treated$b_casi_total_z[i])^2
  j <- which.min(d2)
  matches[[i]] <- tibble::tibble(rzv_subject = treated$subject[i], zvl_subject = cand$subject[j], enroll_gap = w)
  avail <- avail %>% filter(subject != cand$subject[j])
}
match_tbl <- dplyr::bind_rows(matches)
cat(sprintf("\nMatched %d / %d RZV-only(18-19) recipients to a ZVL-only partner (%.1f%%)\n",
            nrow(match_tbl), nrow(treated), 100*nrow(match_tbl)/nrow(treated)))
cat("Enrollment-year gap distribution (0 = exact year match):\n")
print(table(match_tbl$enroll_gap))
cat(sprintf("Mean |enrollment year gap|: %.2f\n", mean(match_tbl$enroll_gap)))

# covariate balance on the matched sets
bal <- baseline %>% filter(subject %in% c(match_tbl$rzv_subject, match_tbl$zvl_subject)) %>%
  mutate(grp = if_else(subject %in% match_tbl$rzv_subject, "RZV", "ZVL"))
smd <- function(v) {
  t <- bal[[v]][bal$grp=="RZV"]; c <- bal[[v]][bal$grp=="ZVL"]
  (mean(t,na.rm=TRUE)-mean(c,na.rm=TRUE)) / sqrt((var(t,na.rm=TRUE)+var(c,na.rm=TRUE))/2)
}
cat("\nCovariate balance (SMD) after matching:\n")
for (v in c("enroll_year","b_visit_age","education","b_casi_total")) {
  cat(sprintf("%-14s SMD=%.3f  (RZV mean=%.2f, ZVL mean=%.2f)\n", v, smd(v),
              mean(bal[[v]][bal$grp=="RZV"],na.rm=TRUE), mean(bal[[v]][bal$grp=="ZVL"],na.rm=TRUE)))
}

# run the event studies and formal difference test on the matched sets
zvl_d <- panel %>% filter(subject %in% match_tbl$zvl_subject, !is.na(casi_total)) %>% build_event_time("zost1_date")
rzv_d <- panel %>% filter(subject %in% match_tbl$rzv_subject, !is.na(casi_total)) %>% build_event_time("shin1_date")

res_zvl <- run_event_study(zvl_d, "casi_total")
res_rzv <- run_event_study(rzv_d, "casi_total")
cat(sprintf("\nZVL (matched) N=%d, pretrend p=%.4f\n", n_distinct(zvl_d$subject), res_zvl$pretrend_p))
cat(sprintf("RZV (matched) N=%d, pretrend p=%.4f\n", n_distinct(rzv_d$subject), res_rzv$pretrend_p))

zvl_d$is_rzv <- 0; rzv_d$is_rzv <- 1
pooled <- bind_rows(zvl_d, rzv_d)
inter_terms <- paste0(es_regs, ":is_rzv")
fml <- as.formula(paste0("casi_total ~ ", paste(es_regs, collapse=" + "), " + ", paste(inter_terms, collapse=" + "), " | subject + visit_year + visit_age"))
m <- fixest::feols(fml, data = pooled, vcov = "hetero", notes = FALSE)
post_terms <- paste0("trt_p",0:6); inter_post <- intersect(paste0(post_terms,":is_rzv"), names(coef(m)))
b <- coef(m)[inter_post]; V <- vcov(m)[inter_post,inter_post,drop=FALSE]
Vinv <- tryCatch(solve(V), error=function(e) MASS::ginv(V))
W <- as.numeric(t(b) %*% Vinv %*% b); p_diff <- 1-pchisq(W, df=length(b))
cat(sprintf("\nJoint Wald test, RZV vs ZVL (matched, direct): p=%.4f (df=%d)\n", p_diff, length(b)))
ct <- as.data.frame(summary(m)$coeftable); ct$term <- rownames(ct)
print(ct[ct$term %in% inter_post, c("term","Estimate","Std. Error","Pr(>|t|)")], row.names = FALSE)

p <- es_plot_two(res_rzv$mod, res_zvl$mod, "RZV (2018-2019, matched)", "ZVL (enrollment-year-matched)",
                  title = NULL,
                  ytitle = "Effect on CASI Total Score")
ggsave("output/figures/RZV_ZVL_direct_matched.png", p, width = 10, height = 6, dpi = 220)
cat("saved plot\n")
