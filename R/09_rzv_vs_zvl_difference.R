# ==============================================================================
# 09_rzv_vs_zvl_difference.R
#
# RZV event study and a formal ZVL-vs-RZV difference test, ported from the
# original reviewer-response analysis (reviewer3.R PART C, "C1"/"C4"):
#
#   C1: RZV event study, calendar-year event time (not the primary ZVL
#       spec's day-based time), window [-8,+8], three subsets -- full RZV
#       cohort, vaccinated 2018-2019 only (pre-COVID-vaccine-era doses),
#       and visits before 2020 only. (The original's "C2" COVID-era x
#       event-time term is not estimable -- collinear with the calendar-
#       year fixed effect -- so the original uses these subsets instead;
#       see also R/10 for a matched-control alternative to the same
#       COVID-era concern.)
#   C4: ZVL vs RZV head-to-head in one pooled model, each subject's own
#       vaccine's first-dose year as its event-time origin, narrow
#       [-4,+4] window (2-year bins), with an RZV interaction term giving
#       the formal RZV-minus-ZVL difference at each horizon.
#
# Both use calendar-year (not day-level) event time and FE = subject +
# visit_year + visit_age, vcov = "hetero", matching the original exactly.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })

panel <- readRDS("data/analysis_panel.rds")

draw <- panel %>%
  mutate(zvl_year = ifelse(!is.na(zost1_date), lubridate::year(zost1_date), NA_integer_),
         rzv_year = ifelse(!is.na(shin1_date), lubridate::year(shin1_date), NA_integer_),
         hz_group = case_when(!is.na(zost1_date) ~ "ZVL",
                               !is.na(shin1_date) & is.na(zost1_date) ~ "RZV",
                               TRUE ~ "Never"))

cget <- function(m, k) if (k %in% names(coef(m)))
  sprintf("%.2f (p=%.3f)", coef(m)[k], 2 * pnorm(-abs(coef(m)[k] / sqrt(vcov(m)[k, k])))) else "--"

# ------------------------------------------------------------------------
# C1: RZV event study, full cohort vs pre-2020 vaccinees vs pre-2020 visits
# ------------------------------------------------------------------------
cat("################################################################\n")
cat("C1: RZV event study -- full cohort vs pre-2020 vaccinees vs pre-2020 visits\n")
cat("################################################################\n")
dr <- draw %>% filter(hz_group == "RZV", !is.na(rzv_year)) %>%
  mutate(et = pmin(pmax(2 * floor((visit_year - rzv_year) / 2), -8), 8))
for (i in 1:4) dr[[paste0("r_m", i)]] <- as.integer(dr$et == -2 * i)
for (i in 0:4) dr[[paste0("r_p", i)]] <- as.integer(dr$et == 2 * i)
dr$r_m1 <- 0L
rr <- paste(c(paste0("r_m", 4:2), paste0("r_p", 0:4)), collapse = " + ")
runr <- function(d) fixest::feols(as.formula(paste0("casi_total ~ ", rr, " | subject + visit_year + visit_age")),
                                   filter(d, !is.na(casi_total)), vcov = "hetero", fixef.iter = 1e5)

c1 <- tibble::tibble(
  subset = c("full RZV cohort", "vaccinated 2018-2019 only", "visits before 2020 only"),
  n_ind = c(n_distinct(dr$subject),
            n_distinct(dr$subject[dr$rzv_year <= 2019]),
            n_distinct(dr$subject[dr$visit_year <= 2019])),
  eff_p1 = c(cget(runr(dr), "r_p1"),
             cget(runr(filter(dr, rzv_year <= 2019)), "r_p1"),
             cget(runr(filter(dr, visit_year <= 2019)), "r_p1")),
  eff_p2 = c(cget(runr(dr), "r_p2"),
             cget(runr(filter(dr, rzv_year <= 2019)), "r_p2"),
             cget(runr(filter(dr, visit_year <= 2019)), "r_p2")))
readr::write_csv(c1, "output/tables/TableS9_rzv_covid_subsets.csv")
print(as.data.frame(c1), row.names = FALSE)

# ------------------------------------------------------------------------
# C4: ZVL vs RZV head-to-head, one pooled model, +/-4y window
# ------------------------------------------------------------------------
cat("\n################################################################\n")
cat("C4: ZVL vs RZV head-to-head, pooled model with an RZV interaction\n")
cat("################################################################\n")
hh <- draw %>% filter(hz_group %in% c("ZVL", "RZV")) %>%
  mutate(vy = dplyr::coalesce(zvl_year, rzv_year),
         et = pmin(pmax(2 * floor((visit_year - vy) / 2), -4), 4),
         rzv = as.integer(hz_group == "RZV"))
for (i in 1:2) hh[[paste0("h_m", i)]] <- as.integer(hh$et == -2 * i)
for (i in 0:2) hh[[paste0("h_p", i)]] <- as.integer(hh$et == 2 * i)
hh$h_m1 <- 0L
hk <- c("h_m2", "h_p0", "h_p1", "h_p2")
mh <- fixest::feols(as.formula(paste0("casi_total ~ ", paste(hk, collapse = " + "), " + ",
  paste(sprintf("%s:rzv", hk), collapse = " + "), " | subject + visit_year + visit_age")),
  filter(hh, !is.na(casi_total)), vcov = "hetero", fixef.iter = 1e5)
c4 <- broom::tidy(mh, conf.int = TRUE) %>%
  filter(term %in% c(hk, paste0(hk, ":rzv"))) %>%
  transmute(term, estimate = round(estimate, 2), ci = sprintf("[%.2f, %.2f]", conf.low, conf.high),
            p = signif(p.value, 2))
readr::write_csv(c4, "output/tables/TableS9_zvl_vs_rzv_headtohead.csv")
print(as.data.frame(c4), row.names = FALSE)
cat("\nThe *:rzv rows are the formal RZV-minus-ZVL difference at each horizon.\n")
