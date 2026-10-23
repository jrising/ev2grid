import pandas as pd
import numpy as np
from plotnine import *
from datetime import datetime
from pathlib import Path

# Ensure plots directory exists
PLOTS_DIR = Path("plots_python")
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Constants from the Julia simulation (customer.jl)
SOC_MIN = 0.3
SOC_MAX = 0.95
VEHICLE_CAPACITY = 75.7  # kWh
REGNEUTRAL = 0.5
MIN_PLUGGED_FOR_BAND = 0.5  # vehicles; hide regulation band below this


def _drop_tail_by_approach(df, dropfromend):
    if dropfromend <= 0:
        return df.copy()

    group_columns = [
        column for column in ['Scenario', 'Approach']
        if column in df.columns
    ]
    if not group_columns:
        return df.iloc[:-dropfromend].copy()

    groups = []
    for _, group in df.groupby(group_columns, sort=False, observed=True):
        groups.append(group.iloc[:-dropfromend])
    return pd.concat(groups, ignore_index=True)


def _prepare_standard_plot_data(df, dropfromend=12):
    # Python equivalent of src/plotting.jl's plot_standard data preparation.
    plot_df = _drop_tail_by_approach(df, dropfromend)
    # vehicles_plugged is a continuous count that decays toward (but never
    # exactly reaches) zero while driving, so a `> 0` test would keep showing
    # soc_plugged. Weight by the fleet split instead, as fill_soc_matrix does.
    group_columns = [
        column for column in ['Scenario', 'Approach']
        if column in plot_df.columns
    ]
    if group_columns:
        fleet_size = plot_df.groupby(
            group_columns, observed=True
        )['vehicles_plugged'].transform('max')
    else:
        fleet_size = plot_df['vehicles_plugged'].max()
    plot_df['soc'] = (
        (fleet_size - plot_df['vehicles_plugged']) * plot_df['soc_driving'] +
        plot_df['vehicles_plugged'] * plot_df['soc_plugged']
    ) / fleet_size

    return plot_df


def plot_standard(df, dropfromend=12, background_layers=None):
    plot_df = _prepare_standard_plot_data(df, dropfromend)
    group_columns = [
        column for column in ['Scenario', 'Approach']
        if column in plot_df.columns
    ]
    if group_columns:
        plot_df['LineGroup'] = plot_df[group_columns].astype(str).agg(
            ': '.join, axis=1
        )
    else:
        plot_df['LineGroup'] = 'SOC'

    line_aesthetics = {
        'x': 'datetime',
        'y': 'soc',
        'group': 'LineGroup',
    }
    if 'Approach' in plot_df.columns:
        line_aesthetics['color'] = 'Approach'

    if 'Scenario' in plot_df.columns:
        line_aesthetics['linetype'] = 'Scenario'
    line_mapping = aes(**line_aesthetics)

    pp = ggplot()
    if background_layers:
        for make_layer in background_layers:
            pp += make_layer(plot_df)

    pp += geom_line(data=plot_df, mapping=line_mapping, size=1)
    if 'Scenario' in plot_df.columns:
        pp += scale_linetype_manual(values={
            'Arbitrage only': 'dashed',
            'With regulation services': 'solid',
        })
    pp += labs(x=None, y='State of Charge (SOC)',
               color='Approach:', linetype='Scenario:')
    pp += theme_bw()

    return pp


def add_standard_annotations(pp):
    return (pp +
            geom_vline(xintercept=[datetime(2023, 7, 17, 9, 0, 0),
                                   datetime(2023, 7, 17, 17, 0, 0)]) +
            geom_hline(yintercept=[SOC_MIN, SOC_MAX], linetype='dotted',
                       color='red', alpha=0.5) +
            scale_x_datetime(name='', date_labels='%H:%M', expand=(0, 0),
                             limits=[datetime(2023, 7, 17, 0, 0, 0),
                                     datetime(2023, 7, 18, 0, 0, 0)]) +
            scale_y_continuous(limits=[0, 1],
                               breaks=[0, 0.3, 0.5, 0.7, 0.95, 1]))


# Read the regulation range data
regulation_df = pd.read_csv("results/bytime_regrange.csv")
regulation_df['datetime'] = pd.to_datetime(regulation_df['datetime'])
regulation_df['Scenario'] = 'With regulation services'

# Convert regrange_kw to regrange_soc (SOC units)
# regrange_kw is in kW, convert to fraction change of total plugged capacity
regulation_df['regrange_soc'] = np.where(
    regulation_df['vehicles_plugged'] > 0,
    (REGNEUTRAL / 2) * regulation_df['regrange_kw'] /
    (regulation_df['vehicles_plugged'] * VEHICLE_CAPACITY),
    0
)

# Calculate regulation band bounds (the range vehicle could operate within)
# Upper bound: min(soc_plugged + regrange_soc, SOC_MAX)
# Lower bound: max(soc_plugged - regrange_soc, SOC_MIN)
regulation_df['reg_upper'] = np.minimum(
    regulation_df['soc_plugged'] + regulation_df['regrange_soc'], SOC_MAX
)
regulation_df['reg_lower'] = np.maximum(
    regulation_df['soc_plugged'] - regulation_df['regrange_soc'], SOC_MIN
)

# For visualization, when no regulation is offered, set bounds to SOC level
no_regulation_offer = regulation_df['regrange_soc'] == 0
regulation_df.loc[no_regulation_offer, 'reg_upper'] = regulation_df.loc[
    no_regulation_offer, 'soc_plugged'
]
regulation_df.loc[no_regulation_offer, 'reg_lower'] = regulation_df.loc[
    no_regulation_offer, 'soc_plugged'
]

# Read the no-regulation data and combine both cases on one set of axes.
no_regulation_df = pd.read_csv("results/bytime.csv")
no_regulation_df['datetime'] = pd.to_datetime(no_regulation_df['datetime'])
no_regulation_df['Scenario'] = 'Arbitrage only'

scenario_order = [
    'Arbitrage only',
    'With regulation services',
]
comparison_df = pd.concat(
    [no_regulation_df, regulation_df], ignore_index=True, sort=False
)
comparison_df['Scenario'] = pd.Categorical(
    comparison_df['Scenario'], categories=scenario_order, ordered=True
)


def regulation_ribbon(plot_df):
    ribbon_df = plot_df[
        plot_df['Scenario'] == 'With regulation services'
    ].copy()
    # Only show the band while vehicles are plugged in. Split into contiguous
    # plugged-in segments so the ribbon doesn't bridge the driving hours.
    plugged = ribbon_df['vehicles_plugged'] >= MIN_PLUGGED_FOR_BAND
    segment = (plugged != plugged.groupby(
        [ribbon_df['Approach']], observed=True
    ).shift(fill_value=False)).groupby(
        ribbon_df['Approach'], observed=True
    ).cumsum()
    ribbon_df['band_group'] = (ribbon_df['Approach'].astype(str) + '_' +
                               segment.astype(str))
    ribbon_df = ribbon_df[plugged]
    return geom_ribbon(
        data=ribbon_df,
        mapping=aes(x='datetime', ymin='reg_lower', ymax='reg_upper',
                    fill='Approach', group='band_group'),
        alpha=0.25,
        inherit_aes=False
    )


# Plot 1: arbitrage only (no regulation services)
arbitrage_df = comparison_df[
    comparison_df['Scenario'] == 'Arbitrage only'
].copy()
arbitrage_df['Scenario'] = arbitrage_df['Scenario'].cat.remove_unused_categories()

arbitrage_plot = add_standard_annotations(
    plot_standard(arbitrage_df) + guides(linetype=False)
)

arbitrage_plot.show()

arbitrage_plot.save(filename=str(PLOTS_DIR / "soc_arbitrage.png"),
                    width=10, height=6, dpi=300)

# Plot 2: arbitrage plus regulation services scenarios
comparison_plot = add_standard_annotations(
    plot_standard(comparison_df, background_layers=[regulation_ribbon]) +
    labs(fill='Regulation band:')
)

comparison_plot.show()

comparison_plot.save(filename=str(PLOTS_DIR / "soc_comparison.png"),
                     width=10, height=6, dpi=300)
