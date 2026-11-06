# ==============================================================================
# 06_figure1_plot.R
#
# Coefficient plot + companion table for Figure 1 (main CASI event study).
# All of Figure 4 (panels A-E) and the 7 subgroup-split HTE figures are
# produced in Stata -- see stata/HTE_all.do -- since that is the style
# actually used for every published HTE panel (confirmed against the real
# embedded Figure 4A image).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(haven) })
source("R/01_event_study_utils.R")

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

panel <- haven::read_dta("data/raw/event_all_120125.dta") %>%
  mutate(zost1_date = as.Date(zost1_date, origin = "1960-01-01"),
         visit_date  = as.Date(visit_date,  origin = "1960-01-01"))

primary <- panel %>%
  filter(!is.na(zost1_date)) %>%
  build_event_time("zost1_date")

res1 <- run_event_study(primary, "casi_total")
tbl1 <- es_table(res1$mod)
p1 <- es_plot(tbl1, title = NULL)
ggsave("output/figures/Figure1_main_event_study.png", p1, width = 8, height = 5.5, dpi = 300)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(tbl1, "output/tables/TableS_main_casi_total.csv")
cat(sprintf("Main CASI pretrend joint Wald p = %.4f (paper: 0.452)\n", res1$pretrend_p))

cat("Saved Figure 1 plot + companion table to output/\n")
