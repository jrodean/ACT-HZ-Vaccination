# ==============================================================================
# 27_finalize_table_s22.R
#
# Appends the four subgroup splits built only in stata/HTE_all.do (frailty,
# prior other vaccination, prior herpes simplex, income) to Table S22,
# using the sample sizes that script logs to TableS22_stata_sample_sizes.csv.
# Runs after the Stata step in run_all.R, since that CSV does not exist
# until stata/HTE_all.do has run.
# ==============================================================================
suppressPackageStartupMessages(library(dplyr))

r_part <- readr::read_csv("output/tables/TableS22_hte_sample_sizes.csv", show_col_types = FALSE)
stata_path <- "output/tables/TableS22_stata_sample_sizes.csv"

if (file.exists(stata_path)) {
  stata_part <- readr::read_csv(stata_path, show_col_types = FALSE) %>%
    mutate(Panel = paste0(Panel, ifelse(Level == 1, " -- yes", " -- no"))) %>%
    select(Panel, N_individuals, N_person_visits)
  tbl_s22_full <- bind_rows(r_part, stata_part)
  readr::write_csv(tbl_s22_full, "output/tables/TableS22_hte_sample_sizes.csv")
  cat(sprintf("Table S22 finalized: %d rows (%d from R, %d from Stata).\n",
              nrow(tbl_s22_full), nrow(r_part), nrow(stata_part)))
} else {
  cat("TableS22_stata_sample_sizes.csv not found -- run stata/HTE_all.do first. Table S22 left as the R-only partial version.\n")
}
