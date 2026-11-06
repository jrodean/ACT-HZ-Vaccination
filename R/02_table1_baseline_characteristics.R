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

# Cohort = any ZVL recipient (zost1_date non-missing); NOT additionally
# restricted to rows with non-missing casi_total -- verified against the
# frozen event_all_120125.dta that this exact rule gives 1,532 individuals /
# 10,962 visit-rows, matching the paper's "1,532 individuals with 10,962
# CASI cognitive tests" exactly (see 00b_validate_data_prep.R).
primary <- panel %>%
  filter(!is.na(zost1_date)) %>%
  build_event_time("zost1_date")

vaxzo_age_frozen <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  distinct(subject, .keep_all = TRUE) %>% select(subject, vaxzo_age)

baseline <- primary %>%
  distinct(subject, .keep_all = TRUE) %>%
  left_join(vaxzo_age_frozen, by = "subject") %>%
  rename(vax_age = vaxzo_age) %>%
  transmute(
    subject,
    b_visit_age, b_visit_year, vax_age, gender, white, education, educat4, race_c,
    b_casi_total, b_casi_att, b_casi_mmc, b_casi_ori, b_casi_ltm, b_casi_stm,
    b_casi_lang, b_casi_draw, b_casi_flu, b_casi_ataj,
    college, adl_sum
  )

n_ind <- n_distinct(primary$subject)
n_tests <- nrow(primary)  # "CASI tests" = ZVL-cohort visit-rows, per paper text (see 00b)

fmt_cont <- function(x) sprintf("%.2f (SD %.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_pct  <- function(x) sprintf("%.1f%%", 100 * mean(x, na.rm = TRUE))
yr_bin <- cut(baseline$b_visit_year, c(1994, 2000, 2005, 2010, 2015, 2020),
              labels = c("1995-2000", "2001-2005", "2006-2010", "2011-2015", "2016-2020"))

tbl <- tibble::tibble(
  Characteristic = c(
    "N individuals", "N CASI cognitive tests",
    "Female", "1995-2000", "2001-2005", "2006-2010", "2011-2015", "2016-2020",
    "<70 years (baseline)", "70-74 years (baseline)", "75-79 years (baseline)", "80+ years (baseline)",
    "<70 years (vaccination)", "70-79 years (vaccination)", "80-89 years (vaccination)", "90+ years (vaccination)",
    "Less than high school", "Completed high school", "Some college", "Completed college",
    "American Indian or Alaska Native", "Asian", "Black or African American", "White",
    "Some other race or More than one race", "Unknown or Not Reported",
    "CASI Total Score", "CASI Attention", "CASI Mental Manipulation", "CASI Orientation",
    "CASI Long-term Memory", "CASI Short-term Memory", "CASI Language",
    "CASI Visual Construction (Drawing)", "CASI Verbal Fluency", "CASI Abstraction and Judgment"
  ),
  Value = c(
    format(n_ind, big.mark = ","),
    format(n_tests, big.mark = ","),
    fmt_pct(baseline$gender == 2),
    fmt_pct(yr_bin == "1995-2000"), fmt_pct(yr_bin == "2001-2005"), fmt_pct(yr_bin == "2006-2010"),
    fmt_pct(yr_bin == "2011-2015"), fmt_pct(yr_bin == "2016-2020"),
    fmt_pct(baseline$b_visit_age < 70), fmt_pct(baseline$b_visit_age >= 70 & baseline$b_visit_age < 75),
    fmt_pct(baseline$b_visit_age >= 75 & baseline$b_visit_age < 80), fmt_pct(baseline$b_visit_age >= 80),
    fmt_pct(baseline$vax_age < 70), fmt_pct(baseline$vax_age >= 70 & baseline$vax_age < 80),
    fmt_pct(baseline$vax_age >= 80 & baseline$vax_age < 90), fmt_pct(baseline$vax_age >= 90),
    fmt_pct(baseline$educat4 == 0), fmt_pct(baseline$educat4 == 1),
    fmt_pct(baseline$educat4 == 2), fmt_pct(baseline$educat4 == 3),
    fmt_pct(baseline$race_c == 1), fmt_pct(baseline$race_c == 2), fmt_pct(baseline$race_c == 3),
    fmt_pct(baseline$race_c == 7), fmt_pct(baseline$race_c %in% c(6, 8)), fmt_pct(baseline$race_c == 9),
    fmt_cont(baseline$b_casi_total), fmt_cont(baseline$b_casi_att), fmt_cont(baseline$b_casi_mmc),
    fmt_cont(baseline$b_casi_ori), fmt_cont(baseline$b_casi_ltm), fmt_cont(baseline$b_casi_stm),
    fmt_cont(baseline$b_casi_lang), fmt_cont(baseline$b_casi_draw), fmt_cont(baseline$b_casi_flu),
    fmt_cont(baseline$b_casi_ataj)
  )
)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(tbl, "output/tables/Table1_baseline_characteristics.csv")
readr::write_csv(tbl, "output/tables/TableS2_baseline_characteristics.csv")
print(tbl, n = Inf)

cat("\n--- paper text comparison ---\n")
cat("paper: N=1,532; 10,962 tests; age 71.8 (SD 5.4); 58% female; 88.3% White;\n")
cat("       education 15.9 (SD 3.1); baseline CASI 94.98\n")
