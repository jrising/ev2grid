using Dates, StatsBase

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
function simustep_alldrive(vehicles_plugged_1::Float64, vehicles_avail_1::Float64, enerfrac_avail_1::Float64, enerfrac_driving_1::Float64)
    ## Move all cars to driving mode
    enerfrac_avail_2 = enerfrac_avail_1
    if vehicles_plugged_1 < vehicles && vehicles - vehicles_plugged_1 + vehicles_avail_1 > 0
        enerfrac_driving_afterdrive = (enerfrac_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
        enerfrac_driving_2 = (enerfrac_driving_afterdrive * (vehicles - vehicles_plugged_1) + enerfrac_avail_1 * vehicles_avail_1) / (vehicles - vehicles_plugged_1 + vehicles_avail_1)
    else
        enerfrac_driving_2 = enerfrac_avail_1
    end
    vehicles_avail_2 = 0.

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(enerfrac_avail_2)
    @assert !isnan(enerfrac_driving_2)

    return vehicles_avail_2, enerfrac_avail_2, enerfrac_driving_2
end

"""
Simulates one time step where all vehicles are plugged in.
Adjusts energy fractions based on the change in vehicle states and returns the new state variables.
"""
function simustep_allplug(vehicles_plugged_1::Float64, vehicles_avail_1::Float64, enerfrac_avail_1::Float64, enerfrac_driving_1::Float64)
    ## Move all cars to plugged in
    if vehicles_avail_1 + vehicles - vehicles_plugged_1 > 0
        enerfrac_avail_2 = (enerfrac_avail_1 * vehicles_avail_1 + enerfrac_driving_1 * (vehicles - vehicles_plugged_1)) / (vehicles_avail_1 + vehicles - vehicles_plugged_1)
    else
        enerfrac_avail_2 = enerfrac_avail_1
    end
    enerfrac_driving_2 = enerfrac_driving_1
    vehicles_avail_2 = vehicles_avail_1 + vehicles - vehicles_plugged_1

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(enerfrac_avail_2)
    @assert !isnan(enerfrac_driving_2)

    return vehicles_avail_2, enerfrac_avail_2, enerfrac_driving_2
end

"""
Simulates the response of vehicles to an event, where a specified
number of vehicles are needed and moved to the driving state.

Updates the state variables accordingly and ensures output values are not `NaN`.
"""
function simustep_event(vehicles_needed::Float64, vehicles_plugged_1::Float64, vehicles_avail_1::Float64, enerfrac_avail_1::Float64, enerfrac_driving_1::Float64)
    ## Move vehicles_needed cars to driving mode
    vehicles_sent = min(vehicles_avail_1, vehicles_needed)
    vehicles_avail_2 = vehicles_avail_1 - vehicles_sent
    enerfrac_avail_2 = enerfrac_avail_1
    if vehicles_plugged_1 < vehicles && vehicles - vehicles_plugged_1 + vehicles_sent > 0
        enerfrac_driving_afterdrive = (enerfrac_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
        enerfrac_driving_2 = (enerfrac_driving_afterdrive * (vehicles - vehicles_plugged_1) + enerfrac_avail_1 * vehicles_sent) / (vehicles - vehicles_plugged_1 + vehicles_sent)
    else
        enerfrac_driving_2 = enerfrac_avail_1
    end

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(enerfrac_avail_2)
    @assert !isnan(enerfrac_driving_2)

    return vehicles_avail_2, enerfrac_avail_2, enerfrac_driving_2
end

"""
Performs a base simulation step without any events, adjusting energy fractions based on current driving status.
Maintains the current number of available vehicles.
"""
function simustep_base(vehicles_plugged_1::Float64, vehicles_avail_1::Float64, enerfrac_avail_1::Float64, enerfrac_driving_1::Float64)
    enerfrac_avail_2 = enerfrac_avail_1
    if vehicles_plugged_1 < vehicles
        enerfrac_driving_2 = (enerfrac_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
        if enerfrac_driving_2 < 0
            enerfrac_driving_2 = 0.
        end
    else
        enerfrac_driving_2 = enerfrac_driving_1
    end
    vehicles_avail_2 = vehicles_avail_1

    @assert !isnan(vehicles_avail_2)
    @assert !isnan(enerfrac_avail_2)
    @assert !isnan(enerfrac_driving_2)

    return vehicles_avail_2, enerfrac_avail_2, enerfrac_driving_2
end

"""
Returns the appropriate simulation step function for the given `dt1` based on a deterministic schedule.
This function determines whether vehicles should drive or plug based on time of day.
"""
function get_simustep_deterministic(dt1::DateTime)
    date_part = Dates.Date(dt1)
    dt_9am = DateTime(date_part, Dates.Time(9, 0, 0))
    dt_5pm = DateTime(date_part, Dates.Time(17, 0, 0))

    if dt_9am - periodstep(1) ≤ dt1 < dt_9am
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
function get_simustep_stochastic(dt1::DateTime)
    date_part = Dates.Date(dt1)
    dt_9am = DateTime(date_part, Dates.Time(9, 0, 0))
    dt_5pm = DateTime(date_part, Dates.Time(17, 0, 0))

    if dt_5pm - periodstep(1) ≤ dt1 < dt_5pm && rand() < prob_delayed_return
        return simustep_base
    end
    if rand() < prob_event_return
        if dt_9am - periodstep(1) ≤ dt1 < dt_5pm - periodstep(1)
            return simustep_alldrive
        else
            return simustep_allplug
        end
    end
    if rand() < prob_event
        vehicles_needed = sample(1:vehicles, Weights(prob_event_vehicles))
        return (vehicles_plugged_1::Float64, vehicles_avail_1::Float64, enerfrac_avail_1::Float64, enerfrac_driving_1::Float64) -> simustep_event(vehicles_needed, vehicles_plugged_1, vehicles_avail_1, enerfrac_avail_1, enerfrac_driving_1)
    end

    return get_simustep_deterministic(dt1)
end