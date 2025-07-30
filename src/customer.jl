################################################################################
# Filename: cusomter.jl
# Information about the fleet of cars available for charging.
################################################################################

using Dates

# Number of vehicles (must be a float)
vehicles = 4.

# Vehicle battery capacity in kilowatt-hours (kWh) per vehicle
vehicle_capacity = 75.7

# Energy consumed by each vehicle per hour (kWh)
driving_energy_perhour = 329 * 42 / 8 / 1e3 # Wh/mile * miles/day * day/8 hr

# Fraction of energy required at all times
soc_min = 0.3

# Maximum fraction of energy allowed
soc_max = 0.95

# Total energy requirement at 9am in kWh
energy_9am = vehicles * 329 * 42 / 1e3 + 0.3 * vehicle_capacity * vehicles # Requirement at 9am: Wh/mile * miles/day
# Fraction of energy required at 9am relative to total capacity
soc_9am = energy_9am / (vehicle_capacity * vehicles)

# Actions
max_charging_kw = 6.6 # Maximum charging rate in kW for Level 2
fracpower_min = -max_charging_kw / vehicle_capacity # discharge in terms of fraction of energy
fracpower_max = max_charging_kw / vehicle_capacity # charging in terms of fraction of energy
efficiency = 0.95 # EFF

"""
    soc_scheduled(dt::DateTime, drive_starts_time::Time) -> Float64

Calculate the energy fraction required at any given future time based on a schedule,
ignoring unscheduled activity.

# Arguments
- `dt::DateTime`: The datetime for which the energy fraction is being calculated.
- `drive_starts_time::Time`: The time of day when driving beings

# Returns
- `Float64`: The fraction of energy required.
"""
function soc_scheduled(dt::DateTime, drive_starts_time::Time)
    date_part = Dates.Date(dt)
    dt_9am = DateTime(date_part, drive_starts_time)

    if dt_9am - periodstep(1) ≤ dt < dt_9am
        return soc_9am
    else
        return soc_min
    end
end

"""
    smart_config()

Returns FF (== PP) so that charging discretely moves between state spaces
"""
function smart_config(FF_proposed)
    max_action = max_charging_kw / vehicle_capacity
    if max_action > soc_max - soc_min
        return FF_proposed
    end

    ## max_action = A * (soc_max - soc_min) / (FF - 2)
    FF_options = floor(FF_proposed / 2):ceil(FF_proposed * 2)
    aas = [(FF - 2) * max_action / (soc_max - soc_min) for FF in FF_options]
    soln = findmin(abs.(aas - round.(aas)))

    return Int64(FF_options[soln[2]])
end
