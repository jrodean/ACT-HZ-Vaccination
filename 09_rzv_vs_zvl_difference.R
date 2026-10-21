# ==============================================================================
# 09_rzv_vs_zvl_difference.R
#
# Formal pooled-interaction test of whether the RZV event-study effect on
# CASI Total Score differs from the ZVL effect, rather than comparing
# confidence intervals visually (Figures S18/S19). Three versions:
#   (a) full sample, each group anchored to its own vaccine's first-dose date
#   (b) RZV restricted to vaccinations before 2020 (ZVL is already all pre-2018)
#   (c) "comparable conditions": both groups narrowed to the years nearest
#       the 2017/2018 ZVL->RZV transition (ZVL 2014-2017, RZV 2018-2019) --
#       the closest feasible approximation of calendar-matched conditions,
#       since ZVL and RZV never actually co-occur in this data (ZVL 2006-2017,
#       RZV 2018-2024, zero calendar overlap).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")

build_group <- function(panel, date_col, is_rzv, year_min = NULL, year_max = NULL) {
  d <- panel %>% filter(!is.na(.data[[date_col]]), !is.na(casi_total))
  if (!is.null(year_min)) d <- d %>% filter(as.numeric(format(.data[[date_col]], "%Y")) >= year_min)
  if (!is.null(year_max)) d <- d %>% filter(as.numeric(format(.data[[date_col]], "%Y")) <= year_max)
  d <- build_event_time(d, date_col)
  d$is_rzv <- is_rzv
  d
}

rzv_diff_test <- function(zvl_d, rzv_d, label) {
  pooled <- bind_rows(zvl_d, rzv_d)
  inter_terms <- paste0(es_regs, ":is_rzv")
  fml <- as.formula(paste0(
    "casi_total ~ ", paste(es_regs, collapse = " + "), " + ",
    paste(inter_terms, collapse = " + "), " | subject + visit_year + visit_age"
  ))
  m <- fixest::feols(fml, data = pooled, vcov = "hetero")

  post_terms <- paste0("trt_p", 0:6)
  inter_post_wanted <- paste0(post_terms, ":is_rzv")
  inter_post <- intersect(inter_post_wanted, names(stats::coef(m)))
  dropped <- setdiff(inter_post_wanted, names(stats::coef(m)))
  if (length(dropped) > 0) cat(sprintf("NOTE: dropped (not estimable -- insufficient RZV follow-up at this horizon): %s\n", paste(dropped, collapse=", ")))
  b <- stats::coef(m)[inter_post]
  V <- stats::vcov(m)[inter_post, inter_post, drop = FALSE]
  Vinv <- tryCatch(solve(V), error = function(e) MASS::ginv(V))
  W <- as.numeric(t(b) %*% Vinv %*% b)
  p_diff <- 1 - pchisq(W, df = length(b))

  # also report the individual interaction coefficients (RZV minus ZVL, per event-time point)
  ct <- as.data.frame(summary(m)$coeftable)
  ct$term <- rownames(ct)
  inter_rows <- ct[ct$term %in% inter_post_wanted, c("term","Estimate","Std. Error","Pr(>|t|)")]
  names(inter_rows) <- c("term","diff_estimate","se","pvalue")

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("N (ZVL)=%d  N (RZV)=%d\n", nobs(fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs,collapse="+")," | subject+visit_year+visit_age")), data=zvl_d, vcov="hetero")),
                                           nobs(fixest::feols(as.formula(paste0("casi_total ~ ", paste(es_regs,collapse="+")," | subject+visit_year+visit_age")), data=rzv_d, vcov="hetero"))))
  cat(sprintf("Joint Wald test, RZV vs ZVL POST-period coefficients differ: p = %.4f\n", p_diff))
  print(inter_rows, row.names = FALSE)
  list(p_diff = p_diff, inter = inter_rows, n_zvl = nrow(distinct(zvl_d, subject)), n_rzv = nrow(distinct(rzv_d, subject)))
}

cat("################################################################\n")
cat("(a) FULL SAMPLE -- each group anchored to its own vaccine's first dose\n")
cat("################################################################\n")
zvl_full <- build_group(panel, "zost1_date", is_rzv = 0)
rzv_full <- build_group(panel, "shin1_date", is_rzv = 1)
res_full <- rzv_diff_test(zvl_full, rzv_full, "Full sample (ZVL 2006-2017 vs RZV 2018-2024)")

cat("\n################################################################\n")
cat("(b) RZV RESTRICTED TO PRE-2020 (removes COVID-vaccine-era RZV doses)\n")
cat("    ZVL unrestricted (already entirely 2006-2017, unaffected by this cut)\n")
cat("################################################################\n")
rzv_pre2020 <- build_group(panel, "shin1_date", is_rzv = 1, year_max = 2019)
res_pre2020 <- rzv_diff_test(zvl_full, rzv_pre2020, "ZVL (2006-2017) vs RZV restricted to 2018-2019")

cat("\n################################################################\n")
cat("(c) 'COMPARABLE CONDITIONS' -- both groups narrowed to years nearest\n")
cat("    the 2017/2018 ZVL->RZV transition (closest feasible calendar match;\n")
cat("    true overlap is impossible, ZVL and RZV never co-occur in this data)\n")
cat("################################################################\n")
zvl_narrow <- build_group(panel, "zost1_date", is_rzv = 0, year_min = 2014, year_max = 2017)
rzv_narrow <- build_group(panel, "shin1_date", is_rzv = 1, year_min = 2018, year_max = 2019)
res_narrow <- rzv_diff_test(zvl_narrow, rzv_narrow, "ZVL restricted to 2014-2017 vs RZV restricted to 2018-2019")

cat("\n\n================ SUMMARY ================\n")
cat(sprintf("(a) Full sample:            p = %.4f  (N_ZVL=%d, N_RZV=%d)\n", res_full$p_diff, res_full$n_zvl, res_full$n_rzv))
cat(sprintf("(b) RZV pre-2020 only:      p = %.4f  (N_ZVL=%d, N_RZV=%d)\n", res_pre2020$p_diff, res_pre2020$n_zvl, res_pre2020$n_rzv))
cat(sprintf("(c) Comparable conditions:  p = %.4f  (N_ZVL=%d, N_RZV=%d)\n", res_narrow$p_diff, res_narrow$n_zvl, res_narrow$n_rzv))
