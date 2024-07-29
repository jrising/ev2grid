################################################################################
# Filename: cusomter.jl
# Description: Information about the fleet of cars available for charging.
#
# Notes:
# - This file contains functions/scripts for <specific task/purpose>.
# - <Any special notes or instructions>.
#
# Usage:
# To use this <function/script>, call <function_name>(<params>) or
# run the script with `julia <your_filename_here>.jl`.
#
################################################################################


# Inputs from current state
vehicles = 4.
vehicle_capacity = 75.7 # kWh / vehicle
driving_energy_perhour = 329 * 42 / 8 / 1e3 # Wh/mile * miles/day * day/8 hr

energy_9am = vehicles * 329 * 42 / 1e3 + 0.3 * vehicle_capacity * vehicles # Requirement at 9am: Wh/mile * miles/day
enerfrac_9am = energy_9am / (vehicle_capacity * vehicles)

"""
Specify the energy scheduled to be needed at any given future time (ignoring unscheduled activity)
"""
function enerfrac_scheduled(dt::DateTime)
    date_part = Dates.Date(dt)
    dt_9am = DateTime(date_part, Dates.Time(9, 0, 0))

    if dt_9am - periodstep(1) ≤ dt < dt_9am
        return enerfrac_9am
    else
        return enerfrac_min
    end
end

# Preferences
weight_portion_below = 1.
