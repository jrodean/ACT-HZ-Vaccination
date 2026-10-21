# ==============================================================================
# 03_domains_and_dementia_outcomes.R
#
# Figure 3 + Tables S8-S16 (9 CASI domain-specific event studies) +
# Table S23-S25 / Figure S26 (CASI<86 / CDR global / DSM-IV dementia) +
# Figure S27/S28 (CDR subdomains, >=1 "dementia" and >0 "MCI or dementia"
# thresholds).
#
# Specification identical to the primary analysis (01_event_study_utils.R),
# only the outcome variable changes.
#
# Dementia-diagnosis outcome coding matches the original Stata analysis
# (analysis_Xi_1225.do): CDR- and DSM-IV-based diagnosis are generated with
# the Stata idiom `var > threshold & var < 8`. Under Stata's missing-value
# ordering (missing sorts as larger than any finite value), this evaluates
# to FALSE whenever the underlying variable is missing, so a missing CDR/
# DSM-IV code is coded as 0 ("no diagnosis") rather than excluded from the
# regression. CDR- and DSM-IV-based outcomes are therefore estimated on the
# full vaccinated cohort (N = 1,532 individuals), not restricted to visits
# with a completed CDR/DSM-IV evaluation. CASI Total Score < 86 has no
# missing-value subtlety and is restricted to visits with a non-missing CASI
# Total Score, matching the primary analysis sample (N = 1,531 individuals).
#
# CDR scores take half-point values (e.g., 0.5 = questionable/MCI), so ">0"
# and ">=1" are distinct thresholds. "Any Stage Dementia Based on CDR"
# (cdr_total_any, used for Table S24) uses ">0" per the original Stata
# code. The CDR-subdomain outcomes report both thresholds for comparison:
# ">=1" (dementia-level impairment, Figure S27) and ">0" (MCI or dementia,
# Figure S28).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(purrr); library(ggplot2); library(patchwork) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")

cdr_domains <- c("memory", "orient", "judge", "comm", "home", "pers")

primary <- panel %>%
  filter(!is.na(zost1_date)) %>%
  build_event_time("zost1_date") %>%
  mutate(
    dsmivd = as.integer(!is.na(dsmiv_dx) & dsmiv_dx %in% 1:6),
    casi85 = as.integer(casi_total < 86),
    cdr_total_any = as.integer(!is.na(cdr_total) & cdr_total > 0 & cdr_total < 8),
    across(all_of(paste0("cdr_", cdr_domains)),
           ~ as.integer(!is.na(.) & . >= 1 & . < 8), .names = "{.col}_any"),
    across(all_of(paste0("cdr_", cdr_domains)),
           ~ as.integer(!is.na(.) & . > 0 & . < 8), .names = "{.col}_0")
  )

# Table numbers match the current Supplement "Figures and Tables" captions:
# Table S8 Short-Term Memory, S9 Orientation, S10 Long-Term Memory,
# S11 Mental Manipulation, S12 Language, S13 Verbal Fluency, S14 Attention,
# S15 Visual Construction, S16 Abstraction and Judgment.
domain_outcomes <- c(casi_stm = "Short-Term Memory (Table S8)", casi_ori = "Orientation (Table S9)",
                      casi_ltm = "Long-Term Memory (Table S10)", casi_mmc = "Mental Manipulation (Table S11)",
                      casi_lang = "Language (Table S12)", casi_flu = "Verbal Fluency (Table S13)",
                      casi_att = "Attention (Table S14)", casi_draw = "Visual Construction (Table S15)",
                      casi_ataj = "Abstraction and Judgment (Table S16)")

# Table S23 (CASI<86), Table S24 (CDR global dementia), Table S25 (DSM-IV
# dementia); combined as Figure S26 ("Effect of ZVL vaccination on dementia
# diagnoses").
dementia_outcomes <- c(casi85 = "CASI Total Score < 86 (Table S23)",
                        cdr_total_any = "CDR global dementia (Table S24)",
                        dsmivd = "DSM-IV dementia diagnosis (Table S25)")

cdr_sub_any <- setNames(paste0("cdr_", cdr_domains, "_any"),
                         c("Memory","Orientation","Judgment and Problem Solving",
                           "Community Affairs","Home and Hobbies","Personal Care"))
cdr_sub_0 <- setNames(paste0("cdr_", cdr_domains, "_0"),
                       c("Memory","Orientation","Judgment and Problem Solving",
                         "Community Affairs","Home and Hobbies","Personal Care"))

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

run_and_save <- function(outcome, label, data, file_prefix) {
  res <- run_event_study(data, outcome)
  tbl <- es_table(res$mod)
  readr::write_csv(tbl, sprintf("output/tables/%s_%s.csv", file_prefix, outcome))
  list(outcome = outcome, label = label, n = nobs(res$mod),
       pretrend_p = res$pretrend_p, tbl = tbl)
}
# panel_grid() is defined in R/01_event_study_utils.R

cat("=== Figure 3 / Tables S8-S16 (CASI domain subscores) ===\n")
domain_table_num <- c(casi_stm = "S8", casi_ori = "S9", casi_ltm = "S10", casi_mmc = "S11",
                       casi_lang = "S12", casi_flu = "S13", casi_att = "S14", casi_draw = "S15",
                       casi_ataj = "S16")
domain_res <- imap(domain_outcomes, function(label, outcome) {
  run_and_save(outcome, label, primary %>% filter(!is.na(.data[[outcome]])),
               sprintf("Table%s_domain", domain_table_num[[outcome]]))
})

for (nm in names(domain_res)) {
  r <- domain_res[[nm]]
  p2 <- r$tbl %>% filter(term == "trt_p1")
  p12 <- r$tbl %>% filter(term == "trt_p6")
  cat(sprintf("%-25s 2y=%.3f (p=%.3f)  12y=%.3f (p=%.3f)  pretrend p=%.3f  N=%d\n",
              r$label, p2$estimate, p2$pvalue, p12$estimate, p12$pvalue, r$pretrend_p, r$n))
}
cat("paper: Orientation 2y=0.203(p=.015) 12y=0.970(p=.001); Short-Term Memory 2y=0.285(p=.001) 12y=0.914(p=.002)\n")

p_fig3 <- panel_grid(domain_res, ytitle_fun = function(l) "Effect on CASI subscore", ncol = 3)
ggsave("output/figures/Figure3_domains.png", p_fig3, width = 14, height = 12, dpi = 200)

cat("\n=== Figure S26 / Tables S23-S25: dementia-diagnosis outcomes ===\n")
dem_table_num <- c(casi85 = "S23", cdr_total_any = "S24", dsmivd = "S25")
dem_res <- imap(dementia_outcomes, function(label, outcome) {
  run_and_save(outcome, label, primary %>% filter(!is.na(.data[[outcome]])),
               sprintf("Table%s_dementia", dem_table_num[[outcome]]))
})
for (nm in names(dem_res)) {
  r <- dem_res[[nm]]
  p12 <- r$tbl %>% filter(term == "trt_p6")
  cat(sprintf("%-40s 12y after: %.3f pp (p=%.3f)  N=%d\n", r$label, 100*p12$estimate, p12$pvalue, r$n))
}
cat("paper (abs. pp at max follow-up): CASI<86 -11.9 (p=.020); CDR global -10.4 (p=.025); DSM-IV -8.8 (p=.048)\n")
cat("(ZVL effect is protective/negative on these harm-coded outcomes.)\n")

p_figS26 <- panel_grid(dem_res, ytitle_fun = function(l) "Probability", ncol = 3)
ggsave("output/figures/FigureS26_dementia_outcomes.png", p_figS26, width = 14, height = 5, dpi = 200)

cat("\n=== Figure S27: CDR subdomains indicating dementia (>=1) ===\n")
cdr1_res <- imap(cdr_sub_any, function(outcome, label) {
  run_and_save(outcome, label, primary %>% filter(!is.na(.data[[outcome]])), "TableS27_cdrsub_dementia")
})
p_figS27 <- panel_grid(cdr1_res, ytitle_fun = function(l) "Probability", ncol = 3)
ggsave("output/figures/FigureS27_cdr_subdomains_dementia.png", p_figS27, width = 14, height = 8, dpi = 200)

cat("\n=== Figure S28: CDR subdomains indicating MCI or dementia (>0) ===\n")
cdr0_res <- imap(cdr_sub_0, function(outcome, label) {
  run_and_save(outcome, label, primary %>% filter(!is.na(.data[[outcome]])), "TableS28_cdrsub_mci")
})
p_figS28 <- panel_grid(cdr0_res, ytitle_fun = function(l) "Probability", ncol = 3)
ggsave("output/figures/FigureS28_cdr_subdomains_mci.png", p_figS28, width = 14, height = 8, dpi = 200)

saveRDS(list(domain = domain_res, dementia = dem_res, cdr_any = cdr1_res, cdr_0 = cdr0_res),
        "output/tables/domains_dementia_results.rds")
cat("\nSaved all domain/dementia/CDR-subdomain result tables AND figures to output/\n")
