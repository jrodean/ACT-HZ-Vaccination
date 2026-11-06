# ==============================================================================
# 08_sunab_estimator.R
#
# Figure S20. Sun & Abraham (2021) interaction-weighted event-study
# estimator, as a robustness check against heterogeneous treatment-effect
# dynamics under staggered ZVL adoption timing, via fixest::sunab() --
# cohort = calendar year of first ZVL dose; subjects vaccinated after 2016
# excluded so 2016 is the last-treated/control cohort; event-time window
# -12 to +8 years (asymmetric -- confirmed against the published figure's
# own axis range, which runs "12 yrs before" to "8 yrs after").
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(fixest); library(ggplot2) })
source("R/01_event_study_utils.R")

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

primary <- panel %>%
  filter(!is.na(zost1_date), !is.na(casi_total)) %>%
  mutate(
    vax_year = lubridate::year(zost1_date),
    event_time = floor(as.numeric(visit_date - zost1_date) / (2 * 365)) * 2,
    rel = pmin(pmax(event_time, -12), 8)
  ) %>%
  filter(vax_year <= 2016)

m_sa <- fixest::feols(casi_total ~ sunab(vax_year, rel, ref.p = -2) | subject + visit_year + visit_age,
                       data = primary, vcov = "hetero")

ct <- broom::tidy(m_sa, conf.int = TRUE)
ct$rel <- suppressWarnings(as.numeric(sub(".*rel::(-?[0-9]+).*", "\\1", ct$term)))
ct <- ct[!is.na(ct$rel), ]
sat <- ct %>% transmute(rel, estimate, ci_low = conf.low, ci_high = conf.high, pvalue = p.value) %>%
  arrange(rel)

readr::write_csv(sat, "output/tables/TableS20_sunab_estimator.csv")
cat("=== Sun-Abraham (interaction-weighted) event study ===\n")
print(sat, n = Inf)

agg <- tryCatch(fixest::aggregate(m_sa, agg = "att"), error = function(e) NULL)
if (!is.null(agg)) cat(sprintf("\nOverall ATT (post-period average): %.3f (SE %.3f)\n", agg[1,1], agg[1,2]))

# Same plain point + dashed-error-bar house style as every other event-study
# panel (confirmed against the published figure directly: no connecting
# line, "Effect on CASI Total Score" with no estimator suffix).
sat_tbl <- sat %>%
  mutate(term = ifelse(rel < 0, paste0("trt_m", -rel / 2), paste0("trt_p", rel / 2))) %>%
  select(term, estimate, ci_low, ci_high, pvalue)
p <- es_plot(sat_tbl, ytitle = "Effect on CASI Total Score")
ggsave("output/figures/FigureS20_sunab_estimator.png", p, width = 8, height = 5.5, dpi = 300)

cat("\nSaved Sun-Abraham table and figure.\n")
