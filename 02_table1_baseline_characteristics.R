# ==============================================================================
# 02_table1_baseline_characteristics.R
#
# Baseline characteristics of the primary ZVL event-study analysis cohort
# ("Table S. Baseline characteristics of the primary analysis cohort" /
#  manuscript "Sample characteristics" paragraph).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  build_event_time("zost1_date")

baseline <- primary %>%
  distinct(subject, .keep_all = TRUE) %>%
  transmute(
    subject,
    b_visit_age, gender, white, education,
    b_casi_total, b_casi_att, b_casi_mmc, b_casi_ori, b_casi_ltm, b_casi_stm,
    b_casi_lang, b_casi_draw, b_casi_flu, b_casi_ataj,
    college, adl_sum
  )

n_ind <- n_distinct(primary$subject)
n_tests <- sum(!is.na(primary$casi_total))

fmt_cont <- function(x) sprintf("%.1f (SD %.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_pct  <- function(x) sprintf("%.1f%%", 100 * mean(x, na.rm = TRUE))

tbl <- tibble::tibble(
  Characteristic = c(
    "N individuals", "N CASI cognitive tests",
    "Age at baseline, years", "Female, %", "White, %", "Years of education",
    "CASI Total Score at baseline",
    "CASI Attention", "CASI Mental Manipulation", "CASI Orientation",
    "CASI Long-term Memory", "CASI Short-term Memory", "CASI Language",
    "CASI Visual Construction (Drawing)", "CASI Verbal Fluency",
    "CASI Abstraction and Judgment"
  ),
  Value = c(
    format(n_ind, big.mark = ","),
    format(n_tests, big.mark = ","),
    fmt_cont(baseline$b_visit_age),
    fmt_pct(baseline$gender == 2),
    fmt_pct(baseline$white == 1),
    fmt_cont(baseline$education),
    sprintf("%.2f", mean(baseline$b_casi_total, na.rm = TRUE)),
    fmt_cont(baseline$b_casi_att), fmt_cont(baseline$b_casi_mmc),
    fmt_cont(baseline$b_casi_ori), fmt_cont(baseline$b_casi_ltm),
    fmt_cont(baseline$b_casi_stm), fmt_cont(baseline$b_casi_lang),
    fmt_cont(baseline$b_casi_draw), fmt_cont(baseline$b_casi_flu),
    fmt_cont(baseline$b_casi_ataj)
  )
)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(tbl, "output/tables/Table1_baseline_characteristics.csv")
print(tbl, n = Inf)

cat("\n--- paper text comparison ---\n")
cat("paper: N=1,532; 10,962 tests; age 71.8 (SD 5.4); 58% female; 88.3% White;\n")
cat("       education 15.9 (SD 3.1); baseline CASI 94.98\n")
