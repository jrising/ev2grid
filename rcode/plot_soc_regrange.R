setwd("~/research/ev2grid/ev2grid")

library(tidyverse)
library(lubridate)

# Ensure plots directory exists
PLOTS_DIR <- "plots_r"
if (!dir.exists(PLOTS_DIR)) {
  dir.create(PLOTS_DIR, recursive = TRUE)
}

# Constants from the simulation
SOC_MIN <- 0.3
SOC_MAX <- 0.95
VEHICLE_CAPACITY <- 75.7  # kWh
NUM_VEHICLES <- 4
REGNEUTRAL <- 0.5

# Read the regulation range data
df <- read_csv("results/bytime_regrange.csv")
df <- df %>% mutate(datetime = as_datetime(datetime))

# Convert regrange_kw to regrange_soc (SOC units)
# regrange_kw is in kW, convert to fraction change of total plugged capacity
df <- df %>% mutate(
  regrange_soc = ifelse(
    vehicles_plugged > 0,
    (REGNEUTRAL / 2) * regrange_kw / (vehicles_plugged * VEHICLE_CAPACITY),
    0
  )
)

# Calculate regulation band bounds (the range vehicle could operate within)
df <- df %>% mutate(
  reg_upper = pmin(soc_plugged + regrange_soc, SOC_MAX),
  reg_lower = pmax(soc_plugged - regrange_soc, SOC_MIN)
)

# For visualization, when no regulation is offered, set bounds to SOC level
df <- df %>% mutate(
  reg_upper = ifelse(regrange_soc == 0, soc_plugged, reg_upper),
  reg_lower = ifelse(regrange_soc == 0, soc_plugged, reg_lower)
)

# Add plot time (middle of hour)
df <- df %>% mutate(plot_time = datetime + minutes(30))
df_day <- df %>% filter(date(datetime) == as_date("2023-07-17"))

plot2 <- ggplot(df_day, aes(x = plot_time)) +
  # Regulation range ribbons
  geom_ribbon(aes(ymin = reg_lower, ymax = reg_upper, fill = Approach), alpha = 0.25) +
  # SOC lines
  geom_line(aes(y = soc_plugged, color = Approach), size = 1) +
  # Driving period markers
  geom_vline(xintercept = as_datetime(c("2023-07-17 09:00:00", "2023-07-17 17:00:00"))) +
  # SOC limits
  geom_hline(yintercept = c(SOC_MIN, SOC_MAX), linetype = "dotted", color = "red", alpha = 0.5) +
  scale_x_datetime(date_labels = "%H:%M", expand = c(0, 0),
                   limits = c(as_datetime("2023-07-17 00:00:00"), as_datetime("2023-07-18 00:00:00"))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.3, 0.5, 0.7, 0.95, 1)) +
  labs(x = NULL, y = "State of Charge (SOC)", fill = "Regulation band:", color = "Approach:") +
  theme_bw()

print(plot2)

ggsave(filename = file.path(PLOTS_DIR, "soc_regrange_comparison.png"), plot = plot2, width = 10, height = 6, dpi = 300)
