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
enerfrac_min = 0.3

# Maximum fraction of energy allowed
enerfrac_max = 1.0

# Total energy requirement at 9am in kWh
energy_9am = vehicles * 329 * 42 / 1e3 + 0.3 * vehicle_capacity * vehicles # Requirement at 9am: Wh/mile * miles/day
# Fraction of energy required at 9am relative to total capacity
enerfrac_9am = energy_9am / (vehicle_capacity * vehicles)


"""
    enerfrac_scheduled(dt::DateTime) -> Float64

Calculate the energy fraction required at any given future time based on a schedule,
ignoring unscheduled activity.

# Arguments
- `dt::DateTime`: The datetime for which the energy fraction is being calculated.

# Returns
- `Float64`: The fraction of energy required.
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

# Penalty for having a portion of vehicles below the minimum energy threshold
weight_portion_below = 1.
