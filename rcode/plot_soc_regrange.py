import pandas as pd
import numpy as np
from plotnine import *
from datetime import datetime, time
from pathlib import Path

# Ensure plots directory exists
PLOTS_DIR = Path("plots_python")
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Constants from the Julia simulation (customer.jl)
SOC_MIN = 0.3
SOC_MAX = 0.95
VEHICLE_CAPACITY = 75.7  # kWh
NUM_VEHICLES = 4

# Read the regulation range data
df = pd.read_csv("results/bytime_regrange.csv")
df['datetime'] = pd.to_datetime(df['datetime'])

# Convert regrange_kwh to regrange_soc (SOC units)
# regrange_kwh is in kWh, convert to fraction of total plugged capacity
df['regrange_soc'] = np.where(
    df['vehicles_plugged'] > 0,
    df['regrange_kwh'] / (df['vehicles_plugged'] * VEHICLE_CAPACITY),
    0
)

# Calculate regulation band bounds (the range vehicle could operate within)
# Upper bound: min(soc_plugged + regrange_soc, SOC_MAX)
# Lower bound: max(soc_plugged - regrange_soc, SOC_MIN)
df['reg_upper'] = np.minimum(df['soc_plugged'] + df['regrange_soc'], SOC_MAX)
df['reg_lower'] = np.maximum(df['soc_plugged'] - df['regrange_soc'], SOC_MIN)

# For visualization, when no regulation is offered, set bounds to SOC level
df.loc[df['regrange_soc'] == 0, 'reg_upper'] = df.loc[df['regrange_soc'] == 0, 'soc_plugged']
df.loc[df['regrange_soc'] == 0, 'reg_lower'] = df.loc[df['regrange_soc'] == 0, 'soc_plugged']

# Add plot time (middle of hour)
df['plot_time'] = df['datetime'] + pd.Timedelta(minutes=30)
df_day = df[df['datetime'].dt.date == pd.to_datetime("2023-07-17").date()].copy()

plot2 = (ggplot(df_day, aes(x='plot_time')) +
         # Regulation range ribbons
         geom_ribbon(aes(ymin='reg_lower', ymax='reg_upper', fill='Approach'), 
                     alpha=0.25) +
         # SOC lines
         geom_line(aes(y='soc_plugged', color='Approach'), size=1) +
         # Driving period markers
         geom_vline(xintercept=[datetime(2023, 7, 17, 9, 0, 0),
                                datetime(2023, 7, 17, 17, 0, 0)],
                    linetype='dashed', color='gray') +
         # SOC limits
         geom_hline(yintercept=[SOC_MIN, SOC_MAX], linetype='dotted', 
                    color='red', alpha=0.5) +
         scale_x_datetime(date_labels='%H:%M',
                         limits=[datetime(2023, 7, 17, 0, 0, 0),
                                datetime(2023, 7, 18, 0, 0, 0)]) +
         scale_y_continuous(limits=[0, 1], breaks=[0, 0.3, 0.5, 0.7, 0.95, 1]) +
         labs(x=None, y='State of Charge (SOC)',
              title='Comparison of SOC and Regulation Range Across Strategies',
              fill='Regulation Band', color='SOC') +
         theme_bw() +
         theme(legend_position='bottom'))

print(plot2)
plot2.save(filename=str(PLOTS_DIR / "soc_regrange_comparison.png"), 
           width=12, height=6, dpi=300)

