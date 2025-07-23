using Dates, StatsBase

include("bizutils.jl")

## Unexpected changes in vehicles
prob_event = 0.01 # per hour, so 1 per 4 days
prob_event_vehicles = [0.5, .25, .125, .125]
prob_event_return = 0.25
prob_delayed_return = 0.1

"""
Simulates one time step where all available vehicles are moved to the driving state.
Updates energy fractions and the number of available vehicles, and ensures no `NaN` values in the output.
Adjust state to reflect the end of the period starting at dt1.
"""
function simustep_alldrive(vehicles_plugged_1::Float64, vehicles_avail_1::Float64, soc_avail_1::Float64, soc_driving_1::Float64)
    ## Move all cars to driving mode
    soc_avail_2 = soc_avail_1
    if vehicles_plugged_1 < vehicles && vehicles - vehicles_plugged_1 + vehicles_avail_1 > 0
        soc_driving_afterdrive = (soc_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
        soc_driving_2 = (soc_driving_afterdrive * (vehicles - vehicles_plugged_1) + soc_avail_1 * vehicles_avail_1) / (vehicles - vehicles_plugged_1 + vehicles_avail_1)
    else
        soc_driving_2 = soc_avail_1
    end
    vehicles_avail_2 = 0.

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(soc_avail_2)
    @assert !isnan(soc_driving_2)

    return vehicles_avail_2, soc_avail_2, soc_driving_2
end

"""
Simulates one time step where all vehicles are plugged in.
Adjusts energy fractions based on the change in vehicle states and returns the new state variables.
"""
function simustep_allplug(vehicles_plugged_1::Float64, vehicles_avail_1::Float64, soc_avail_1::Float64, soc_driving_1::Float64)
    ## Move all cars to plugged in
    if vehicles_avail_1 + vehicles - vehicles_plugged_1 > 0
        soc_avail_2 = (soc_avail_1 * vehicles_avail_1 + soc_driving_1 * (vehicles - vehicles_plugged_1)) / (vehicles_avail_1 + vehicles - vehicles_plugged_1)
    else
        soc_avail_2 = soc_avail_1
    end
    soc_driving_2 = soc_driving_1
    vehicles_avail_2 = vehicles_avail_1 + vehicles - vehicles_plugged_1

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(soc_avail_2)
    @assert !isnan(soc_driving_2)

    return vehicles_avail_2, soc_avail_2, soc_driving_2
end

"""
Simulates the response of vehicles to an event, where a specified
number of vehicles are needed and moved to the driving state.

Updates the state variables accordingly and ensures output values are not `NaN`.
"""
function simustep_event(vehicles_needed::Float64, vehicles_plugged_1::Float64, vehicles_avail_1::Float64, soc_avail_1::Float64, soc_driving_1::Float64)
    ## Move vehicles_needed cars to driving mode
    vehicles_sent = min(vehicles_avail_1, vehicles_needed)
    vehicles_avail_2 = vehicles_avail_1 - vehicles_sent
    soc_avail_2 = soc_avail_1
    if vehicles_plugged_1 < vehicles && vehicles - vehicles_plugged_1 + vehicles_sent > 0
        soc_driving_afterdrive = (soc_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
        soc_driving_2 = (soc_driving_afterdrive * (vehicles - vehicles_plugged_1) + soc_avail_1 * vehicles_sent) / (vehicles - vehicles_plugged_1 + vehicles_sent)
    else
        soc_driving_2 = soc_avail_1
    end

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(soc_avail_2)
    @assert !isnan(soc_driving_2)

    return vehicles_avail_2, soc_avail_2, soc_driving_2
end

"""
Performs a base simulation step without any events, adjusting energy fractions based on current driving status.
Maintains the current number of available vehicles.
"""
function simustep_base(vehicles_plugged_1::Float64, vehicles_avail_1::Float64, soc_avail_1::Float64, soc_driving_1::Float64)
    soc_avail_2 = soc_avail_1
    if vehicles_plugged_1 < vehicles
        soc_driving_2 = (soc_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
        if soc_driving_2 < 0
            soc_driving_2 = 0.
        end
    else
        soc_driving_2 = soc_driving_1
    end
    vehicles_avail_2 = vehicles_avail_1

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(soc_avail_2)
    @assert !isnan(soc_driving_2)

    return vehicles_avail_2, soc_avail_2, soc_driving_2
end

"""
Returns the appropriate simulation step function for the given `dt1` based on a deterministic schedule.
This function determines whether vehicles should drive or plug based on time of day.
"""
function get_simustep_deterministic(dt1::DateTime, drive_starts_time, park_starts_time)
    if drive_starts_time == park_starts_time # can happen under stochastic times
        return simustep_base
    end

    date_part = Dates.Date(dt1)
    dt_9am = DateTime(date_part, drive_starts_time)
    dt_5pm = DateTime(date_part, park_starts_time)

    if dt_9am - periodstep(1) ≤ dt1 < dt_5pm - periodstep(1)
        return simustep_alldrive
    elseif dt_5pm - periodstep(1) ≤ dt1 < dt_5pm
        return simustep_allplug
    else
        return simustep_base
    end
end

"""
Returns a simulation step function for the given `dt1` with stochastic elements, introducing possible events like vehicular needs and delayed returns.
The function chooses between `simustep_alldrive`, `simustep_allplug`, or `simustep_event` probabilistically.
"""
function get_simustep_stochastic(dt1::DateTime, drive_starts_time, park_starts_time)
    date_part = Dates.Date(dt1)
    dt_drive_start = DateTime(date_part, drive_starts_time)
    dt_park_start = DateTime(date_part, park_starts_time)

    rand_delayed_return = rand()
    rand_event_return = rand()
    rand_event = rand()

    if dt_park_start - periodstep(1) ≤ dt1 < dt_park_start && rand_delayed_return < prob_delayed_return
        push!(event_log, (time = dt1, event = :delayed_return, vehicles_affected = vehicles)) ## if there is a delayed return right before cars are about to park, simulate a base step
        return simustep_base
    end
    if rand_event_return < prob_event_return ## this is like a common event that makes everyone stay out longer or plug in at the same time
        if dt_drive_start - periodstep(1) ≤ dt1 < dt_park_start - periodstep(1)
            push!(event_log, (time = dt1, event = :alldrive, vehicles_affected = vehicles))
            return simustep_alldrive
        else
            push!(event_log, (time = dt1, event = :allplug, vehicles_affected = vehicles))
            return simustep_allplug
        end
    end
    if rand_event < prob_event
        vehicles_needed = sample(1:vehicles, Weights(prob_event_vehicles)) ## probability of an emergency requiring an uncertain number of vehicles
        push!(event_log, (time = dt1, event = :emergency, vehicles_needed = vehicles_needed))
        return (vehicles_plugged_1::Float64, vehicles_avail_1::Float64, soc_avail_1::Float64, soc_driving_1::Float64) -> simustep_event(vehicles_needed, vehicles_plugged_1, vehicles_avail_1, soc_avail_1, soc_driving_1)
    end

    return get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
end

"""
Returns a random draw of Normal(starts_time, hours_sd), rounded to the closest timestep
timestep is number of hours per timestep (may be less than 1)
"""
function get_stochastic_times(count::Int64, time::Time, hours_sd::Float64, timestep::Float64)
    delta = hours_sd * randn(count)
    timediff = periodstep.(round.(Int, delta / timestep))
    time + timediff
end
