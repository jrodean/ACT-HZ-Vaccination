# ==============================================================================
# 14_enrollment_matched_pretrend.R
#
# RZV recipients skew heavily toward late ACT enrollment (63% enrolled
# 2010 or later), while the never-treated pool skews heavily toward early
# enrollment (68% enrolled 1994-1996); only about 252 never-treated people
# enrolled 2010-2020 at all. This common-support constraint is a plausible
# source of the pretrend violations seen in the matched designs in scripts
# 10-12, which match on demographics and health but not on enrollment era.
#
# This script: (1) restricts/stratifies matching on enrollment year first
# (nearest available, with an expanding window per subject if needed),
# refining with age, education, and baseline CASI within that stratum;
# (2) assigns each matched never-treated control a pseudo-vaccination date
# equal to their matched RZV (or ZVL) partner's real vaccination date;
# (3) fits the standard single-group event-study model on the matched
# control group alone, using their pseudo-dates, and reports its own joint
# Wald pre-trend p-value (the same test used everywhere else in this
# package), rather than a pooled treated-and-control interaction model.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")

baseline <- panel %>%
  filter(days_since_baseline == 0) %>%
  distinct(subject, .keep_all = TRUE) %>%
  select(subject, enroll_year = visit_year, b_visit_age, education, b_casi_total,
         zost1_date, shin1_date)

match_enrollment_first <- function(vax_date_col, label, max_window = 10, seed = 42) {
  cov <- baseline %>% filter(!is.na(b_visit_age), !is.na(education), !is.na(b_casi_total))
  treated <- cov %>% filter(!is.na(.data[[vax_date_col]]))
  never   <- cov %>% filter(is.na(zost1_date), is.na(shin1_date))

  cat(sprintf("\n=== %s: enrollment-year-stratified matching ===\n", label))
  cat(sprintf("Treated N=%d, never-treated pool N=%d\n", nrow(treated), nrow(never)))

  z_second <- treated %>% select(subject, b_visit_age, education, b_casi_total) %>%
    bind_rows(never %>% select(subject, b_visit_age, education, b_casi_total)) %>%
    mutate(across(c(b_visit_age, education, b_casi_total), ~ as.numeric(scale(.))))
  z_second <- z_second %>% rename_with(~ paste0(., "_z"), -subject)

  treated <- treated %>% left_join(z_second, by = "subject")
  never   <- never   %>% left_join(z_second, by = "subject")

  set.seed(seed)
  treated <- treated[sample(nrow(treated)), ]
  avail <- never
  matches <- vector("list", nrow(treated))
  windows_used <- integer(nrow(treated))

  for (i in seq_len(nrow(treated))) {
    if (nrow(avail) == 0) break
    w <- 0
    repeat {
      cand <- avail %>% filter(abs(enroll_year - treated$enroll_year[i]) <= w)
      if (nrow(cand) > 0 || w >= max_window) break
      w <- w + 1
    }
    if (nrow(cand) == 0) next
    d2 <- (cand$b_visit_age_z - treated$b_visit_age_z[i])^2 +
          (cand$education_z  - treated$education_z[i])^2 +
          (cand$b_casi_total_z - treated$b_casi_total_z[i])^2
    j_local <- which.min(d2)
    chosen <- cand$subject[j_local]
    matches[[i]] <- tibble::tibble(treated_subject = treated$subject[i], control_subject = chosen,
                                    enroll_gap = w, treated_enroll_year = treated$enroll_year[i],
                                    control_enroll_year = cand$enroll_year[j_local])
    windows_used[i] <- w
    avail <- avail %>% filter(subject != chosen)
  }
  match_tbl <- dplyr::bind_rows(matches)

  cat(sprintf("Matched %d / %d treated (%.1f%%)\n", nrow(match_tbl), nrow(treated), 100*nrow(match_tbl)/nrow(treated)))
  cat("Distribution of enrollment-year gap actually used (0 = exact year match):\n")
  print(table(match_tbl$enroll_gap))
  cat(sprintf("Mean |enrollment year gap|: %.2f years (vs. unrestricted matching, which ignored this entirely)\n",
              mean(match_tbl$enroll_gap)))

  # ---- simple pretrend test: fit the matched CONTROL group alone, using
  # their pseudo-dates, with the SAME single-group pretrend test used for
  # the primary analysis everywhere else in this package ----
  pseudo_dates <- baseline %>% filter(subject %in% match_tbl$treated_subject) %>%
    select(treated_subject = subject, pseudo_date = all_of(vax_date_col)) %>%
    inner_join(match_tbl, by = "treated_subject") %>% select(control_subject, pseudo_date)

  control_panel <- panel %>% filter(subject %in% pseudo_dates$control_subject, !is.na(casi_total)) %>%
    left_join(pseudo_dates, by = c("subject" = "control_subject"))
  control_d <- build_event_time(control_panel, "pseudo_date")

  res_control <- run_event_study(control_d, "casi_total")
  cat(sprintf("\nMatched never-treated control group (N=%d subjects, %d rows), own pretrend test\n",
              n_distinct(control_d$subject), nrow(control_d)))
  cat(sprintf("Pretrend joint Wald p (does the matched control show a spurious trend around its pseudo-date?): %.4f\n",
              res_control$pretrend_p))
  print(round(coef(res_control$mod), 3))

  # also report the treated group's own pretrend, for direct comparison
  treated_d <- panel %>% filter(subject %in% match_tbl$treated_subject, !is.na(casi_total)) %>%
    build_event_time(vax_date_col)
  res_treated <- run_event_study(treated_d, "casi_total")
  cat(sprintf("\nFor comparison, matched TREATED group's own pretrend p: %.4f\n", res_treated$pretrend_p))

  list(match_tbl = match_tbl, res_control = res_control, res_treated = res_treated)
}

cat("################################################################\n")
cat("RZV, enrollment-year-stratified matching\n")
cat("################################################################\n")
res_rzv <- match_enrollment_first("shin1_date", "RZV")

cat("\n\n################################################################\n")
cat("ZVL, enrollment-year-stratified matching\n")
cat("################################################################\n")
res_zvl <- match_enrollment_first("zost1_date", "ZVL")

cat("\n\n=== SUMMARY ===\n")
cat(sprintf("RZV: matched control pretrend p = %.4f  (treated group's own pretrend p = %.4f)\n",
            res_rzv$res_control$pretrend_p, res_rzv$res_treated$pretrend_p))
cat(sprintf("ZVL: matched control pretrend p = %.4f  (treated group's own pretrend p = %.4f)\n",
            res_zvl$res_control$pretrend_p, res_zvl$res_treated$pretrend_p))
