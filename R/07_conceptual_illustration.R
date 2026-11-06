# ==============================================================================
# 07_conceptual_illustration.R
#
# Figure S41. Illustrating event study design with individual vaccination
# adoption time. Uses no ACT data -- three hypothetical individuals
# vaccinated at different times, illustrating how the event-study
# transforms raw trajectories (Panel A) by removing the common age-related
# trend (Panel B), aligning on the vaccination event (Panel C), and
# removing each individual's own baseline level (Panel D).
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ---- Parameters ----
calendar_T <- 12
time <- 0:(calendar_T - 1)
units <- c("Unit A","Unit B","Unit C")
treat_time <- c(4,6,8)               # vaccination times
base_levels <- c(20,30,40)           # spacing for A-C

cols_units <- c("Unit A"="#1f77b4", "Unit B"="#ff7f0e", "Unit C"="#2ca02c")

theme_clean <- theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    plot.title = element_text(face="bold", size=18, hjust=0.5)
  )

# ---- Panel A: Original Time Series (downward trend, slope flattens) ----
dfA <- tibble(
  unit = rep(units, each = calendar_T),
  t    = rep(time, times = length(units)),
  base = rep(base_levels, each = calendar_T),
  tau  = rep(treat_time, each = calendar_T)
) |>
  mutate(
    y = base - 2.5 * pmin(t, tau) + (-2.5 + 2) * pmax(0, t - tau)
  )

vline_df_A <- dfA |> distinct(unit, tau, base) |>
  mutate(y_tau = base - 2.5 * tau)

label_df_A <- vline_df_A |>
  transmute(
    unit,
    x = tau + 0.10,
    y = y_tau - 1.0,
    label = "Vaccination Event"
  )

gA <- ggplot(dfA, aes(t, y, group = unit, color = unit)) +
  geom_line(linewidth = 1) +
  geom_text(
    data = label_df_A,
    aes(x = x, y = y, label = label, color = unit),
    hjust = 0.5, vjust = 1.2, size = 6, lineheight = 0.95, show.legend = FALSE
  ) +
  scale_color_manual(
    values = cols_units,
    labels = c("Unit A" = "Individual 1", "Unit B" = "Individual 2", "Unit C" = "Individual 3")
  ) +
  scale_x_continuous(breaks = seq(0, calendar_T - 1, by = 2)) +
  coord_cartesian(ylim = c(-5, 40)) +
  labs(title = "Panel A: Individual Cognition Trajectory", x = "Calendar Time", y = "Outcome") +
  theme_clean +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1)))

# ---- Panel B: Remove Age-Related Change (flat pre, slope +0.3 post) ----
dfB <- tibble(
  unit = rep(units, each = calendar_T),
  t    = rep(time, times = length(units)),
  base = rep(base_levels, each = calendar_T),
  tau  = rep(treat_time, each = calendar_T)
) |>
  mutate(
    y = base - 10 + 2 * pmax(0, t - tau)
  )

vline_df_B <- dfB |> distinct(unit, tau, base) |>
  mutate(y_tau = base - 10)

label_df_B <- vline_df_B |>
  transmute(
    unit,
    x = tau + 0.10,
    y = y_tau - 0.8,
    label = "Vaccination Event"
  )

gB <- ggplot(dfB, aes(t, y, group = unit, color = unit)) +
  geom_line(linewidth = 1) +
  geom_text(
    data = label_df_B,
    aes(x = x, y = y, label = label, color = unit),
    hjust = 0.5, vjust = 1, size = 6, lineheight = 0.95, show.legend = FALSE
  ) +
  scale_color_manual(
    values = cols_units,
    labels = c("Unit A" = "Individual 1", "Unit B" = "Individual 2", "Unit C" = "Individual 3")
  ) +
  scale_x_continuous(breaks = seq(0, calendar_T - 1, by = 2)) +
  coord_cartesian(ylim = c(-5, 40)) +
  labs(title = "Panel B: Remove Age-Related Change", x = "Calendar Time", y = NULL) +
  theme_clean +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1)))

# ---- Panel C: Align at Event (solid pre, dashed forward counterfactual) ----
dfC <- dfB |> mutate(event_time = t - tau)

pre_params <- dfC %>%
  filter(event_time < 0) %>%
  group_by(unit) %>%
  summarise(
    mx = mean(event_time),
    my = mean(y),
    vx = var(event_time),
    cxy = cov(event_time, y),
    .groups = "drop"
  ) %>%
  mutate(
    slope = ifelse(is.na(vx) | vx == 0, 0, cxy / vx),
    intercept = my - slope * mx
  )

max_ev <- max(dfC$event_time)
cf_df <- pre_params %>%
  group_by(unit) %>%
  do({
    et <- seq(0, max_ev, by = 1)
    tibble(event_time = et,
           y = .$intercept + .$slope * et)
  }) %>%
  ungroup()

gC <- ggplot(dfC, aes(event_time, y, group = unit, color = unit)) +
  geom_line(linewidth = 1) +
  geom_line(data = cf_df,
            aes(event_time, y, color = unit, group = unit),
            linewidth = 1, linetype = "22", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "22", color = "red", linewidth = 0.6) +
  scale_color_manual(
    values = cols_units,
    labels = c("Unit A" = "Individual 1",
               "Unit B" = "Individual 2",
               "Unit C" = "Individual 3")
  ) +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-5, 40)) +
  labs(title = "Panel C: Align at Vaccination Event", x = "Event Time", y = "Outcome") +
  annotate("text", x = -0.4, y = -3.5, label = "Normalized\nVaccination Event",
           color = "red", angle = 0, hjust = 1, vjust = 0, size = 5) +
  theme_clean +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.8)))

# ---- Panel D: Relative to own baseline (different linetypes for overlap) ----
end_t <- 10
dfD <- dfC %>%
  group_by(unit) %>%
  mutate(y0 = y[event_time == 0][1],
         y_rel = y - y0 + 5) %>%
  ungroup()

gD <- ggplot(dfD, aes(event_time, y_rel, group = unit, color = unit)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "22", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = 5, linetype = "22", color = "grey40", linewidth = 0.7) +
  scale_color_manual(
    values = cols_units,
    labels = c("Unit A" = "Individual 1",
               "Unit B" = "Individual 2",
               "Unit C" = "Individual 3")
  ) +
  coord_cartesian(xlim = c(-10, end_t), ylim = c(-5, 40)) +
  labs(title = "Panel D: Remove Individual Baseline Level",
       x = "Event Time", y = NULL) +
  annotate("text", x = -0.4, y = -3.5, label = "Normalized\nVaccination Event",
           color = "red", angle = 0, hjust = 1, vjust = 0, size = 5) +
  theme_clean +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1,0),
    legend.background = element_rect(fill = "white", color =  "grey80"),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.8)))

# ---- Combine & Save ----
final_plot <- (gA + gB) / (gC + gD)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/FigureS41_conceptual_illustration.pdf", final_plot,
       width = 15, height = 11, units = "in")
ggsave("output/figures/FigureS41_conceptual_illustration.png", final_plot,
       width = 15, height = 11, units = "in", dpi = 300)

message("Saved: FigureS_conceptual_illustration.pdf and .png")
