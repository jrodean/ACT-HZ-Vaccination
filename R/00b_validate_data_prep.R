# ==============================================================================
# 00b_validate_data_prep.R
#
# Sanity-checks the panel built in 00_data_prep.R against figures reported
# in the manuscript text (Sample characteristics paragraph).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })

panel <- readRDS(file.path("data", "analysis_panel.rds"))

cat("=== Cohort-level checks (paper text, Sample characteristics) ===\n")

n_total <- n_distinct(panel$subject)
cat(sprintf("Total ACT participants in linked cohort: %d  (paper: 5,659)\n", n_total))

hz_status <- panel %>%
  distinct(subject, zost1_date, shin1_date) %>%
  mutate(
    any_zvl = !is.na(zost1_date),
    any_rzv = !is.na(shin1_date),
    grp = case_when(
      any_zvl & any_rzv  ~ "both",
      any_zvl & !any_rzv ~ "zvl_only",
      !any_zvl & any_rzv ~ "rzv_only",
      TRUE               ~ "never"
    )
  )

tab <- hz_status %>% count(grp) %>% mutate(pct = round(100 * n / sum(n), 1))
print(tab)
cat("paper: never 57.0% (3,225); zvl_only 15.0% (850); both 12.1% (682); rzv_only 16.0% (902)\n")
cat(sprintf("any HZ vaccine: %d (%.1f%%)  (paper: 2,434, 43.0%%)\n",
            sum(hz_status$any_zvl | hz_status$any_rzv),
            100 * mean(hz_status$any_zvl | hz_status$any_rzv)))

cat("\n=== Primary ZVL event-study cohort ===\n")
# Eligibility rule used by every downstream script: any ZVL recipient
# (zost1_date non-missing), full stop -- NOT additionally restricted to
# rows with non-missing casi_total. Verified directly against the frozen
# event_all_120125.dta used for the submitted manuscript: filtering there
# to !is.na(zost1_date) alone gives exactly 1,532 individuals and 10,962
# visit-rows, matching the paper's "1,532 individuals with 10,962 CASI
# cognitive tests" exactly (873 of those 10,962 rows have casi_total
# missing -- "CASI tests" counts ZVL-cohort visits, not non-missing CASI
# scores). Individual regressions still drop missing-DV rows internally,
# so this only affects cohort-size/N reporting, not estimated effects.
primary <- panel %>% filter(!is.na(zost1_date))

cat(sprintf("N individuals: %d  (paper: 1,532)\n", n_distinct(primary$subject)))
cat(sprintf("N CASI tests in this cohort: %d  (paper: 10,962)\n", nrow(primary)))

baseline_primary <- primary %>% distinct(subject, .keep_all = TRUE)

cat(sprintf("Mean baseline age: %.1f (SD %.1f)  (paper: 71.8, SD 5.4)\n",
            mean(baseline_primary$b_visit_age, na.rm = TRUE),
            sd(baseline_primary$b_visit_age, na.rm = TRUE)))
cat(sprintf("%% female: %.0f%%  (paper: 58%%)\n",
            100 * mean(baseline_primary$gender == 2, na.rm = TRUE)))
cat(sprintf("%% White: %.1f%%  (paper: 88.3%%)\n",
            100 * mean(baseline_primary$white, na.rm = TRUE)))
cat(sprintf("Mean years education: %.1f (SD %.1f)  (paper: 15.9, SD 3.1)\n",
            mean(baseline_primary$education, na.rm = TRUE),
            sd(baseline_primary$education, na.rm = TRUE)))
cat(sprintf("Mean baseline CASI total: %.2f  (paper: 94.98)\n",
            mean(baseline_primary$b_casi_total, na.rm = TRUE)))
