# ==============================================================================
# 24_descriptive_figures_S1_S42.R
#
# Figure S1. Herpes zoster vaccination uptake in ACT participants over time.
# Figure S42. Observed follow-up across event times relative to ZVL
#   vaccination.
#
# Both are purely descriptive (no regression). Figure S42 counts, for each
# event-time point, how many of the primary analysis cohort are observed at
# least that far before (for pre-periods) or after (for post-periods) their
# own vaccination date.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

panel <- readRDS("data/analysis_panel.rds")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# Figure S1: HZ vaccination uptake over time. Stacked-area share of the
# linked cohort vaccinated against HZ (either vaccine) vs. unvaccinated,
# by calendar year, matching the published figure's style exactly
# (stacked-area %, not a cumulative-count line plot): red dashed line
# marks 2018 (RZV introduction / ZVL phase-out).
# ------------------------------------------------------------------------
uptake <- panel %>%
  distinct(subject, .keep_all = TRUE) %>%
  transmute(subject, zost_year = as.integer(format(zost1_date, "%Y")),
            shin_year = as.integer(format(shin1_date, "%Y")))

years <- 2006:2024
share_uptake <- tibble::tibble(year = years) %>%
  rowwise() %>%
  mutate(
    vaccinated_n = sum(uptake$zost_year <= year | uptake$shin_year <= year, na.rm = TRUE),
    total_n = nrow(uptake),
    share_vaccinated = vaccinated_n / total_n
  ) %>%
  ungroup() %>%
  transmute(year, Vaccinated = share_vaccinated, Unvaccinated = 1 - share_vaccinated) %>%
  tidyr::pivot_longer(c(Unvaccinated, Vaccinated), names_to = "status", values_to = "share") %>%
  mutate(status = factor(status, levels = c("Unvaccinated", "Vaccinated")))

readr::write_csv(share_uptake, "output/tables/FigureS1_vaccination_uptake_data.csv")

p1 <- ggplot(share_uptake, aes(x = year, y = share, fill = status)) +
  geom_area() +
  geom_vline(xintercept = 2018, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(2006, 2024, 2), expand = c(0, 0)) +
  scale_fill_manual(values = c(Unvaccinated = "#F8766D", Vaccinated = "#56B4E9")) +
  labs(x = "Year", y = "Shares", fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/figures/FigureS1_vaccination_uptake.png", p1, width = 8, height = 5.5, dpi = 200)
cat("Saved Figure S1 (vaccination uptake over time).\n")

# ------------------------------------------------------------------------
# Figure S42: observed follow-up across event times relative to ZVL
# vaccination, for the primary analysis cohort.
# ------------------------------------------------------------------------
primary0 <- panel %>% filter(!is.na(zost1_date), !is.na(casi_total))

et <- primary0 %>%
  mutate(event_time = floor(as.numeric(visit_date - zost1_date) / (2 * 365)) * 2,
         event_time = pmax(pmin(event_time, 14), -20)) %>%
  distinct(subject, event_time)

subj_range <- et %>% group_by(subject) %>% summarise(min_et = min(event_time), max_et = max(event_time), .groups = "drop")

event_grid <- seq(-20, 14, by = 2)
followup <- tibble::tibble(event_time = event_grid) %>%
  rowwise() %>%
  mutate(
    n_reach = if (event_time < 0) sum(subj_range$min_et <= event_time)
              else sum(subj_range$max_et >= event_time)
  ) %>%
  ungroup()

readr::write_csv(followup, "output/tables/TableS42_followup.csv")
cat("\n=== Figure S42: observed follow-up by event time ===\n")
print(followup, n = Inf)

followup$lab <- ifelse(followup$event_time == 0, "ZVL vaccination year",
                        ifelse(followup$event_time < 0, paste(-followup$event_time, "yrs before"),
                               paste(followup$event_time, "yrs after")))
p42 <- ggplot(followup, aes(x = event_time, y = n_reach)) +
  geom_col(fill = "#4C6E8C") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = followup$event_time, labels = followup$lab) +
  labs(x = "Time relative to the year of ZVL vaccination", y = "Number of individuals observed") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/figures/FigureS42_followup.png", p42, width = 8, height = 5, dpi = 200)

cat("\nSaved Figure S42 (observed follow-up across event times).\n")
