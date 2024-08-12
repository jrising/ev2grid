using Dates

"""
Simulates the energy fraction of vehicles without executing any strategic adjustments over `SS` timesteps. This function models the evolution of energy fractions for parked and driving vehicles, beginning from an initial date and state, allows optional stochastic adjustments, and outputs the simulation results over time in the form of a DataFrame.

# Arguments
- `dt0::DateTime`: The starting date and time for the simulation.
- `vehicles_plugged_1::Float64`: The initial number of vehicles that are plugged in.
- `enerfrac_plugged_1::Float64`: The initial energy fraction of plugged vehicles.
- `enerfrac_driving_1::Float64`: The initial energy fraction of vehicles that are driving.
- `stochastic::Bool=false`: Optional argument to enable stochastic simulation of energy needs.

# Returns
- `DataFrame`: A DataFrame containing columns for each timestep, including datetime, needed energy fraction, fraction and energy below threshold, plugged vehicle energy, and driving vehicle energy.
"""
function simu_inactive(dt0::DateTime, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64, stochastic::Bool=false)
    # datetime, enerfrac_needed, vehicles_plugged_1, portion_below, enerfrac_toadd_below, enerfrac_avail_1, enerfrac_driving_1
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64}[]

    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

    push!(rows, (dt0 + periodstep(1), enerfrac_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], enerfrac_plugged_1, enerfrac_driving_1))

    for tt in 1:(SS-1)
        println(tt)
        enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

        if stochastic
            simustep = get_simustep_stochastic(dt0 + periodstep(tt))
        else
            simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        end

        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], enerfrac_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt + 1), enerfrac_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], enerfrac_plugged_2, enerfrac_driving_2))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :enerfrac_needed, :vehicles_plugged, :portion_below, :enerfrac_below, :enerfrac_plugged, :enerfrac_driving])
end

"""
Simulates the energy fraction of vehicles while applying a given strategy over several timesteps. The strategy is provided as an array, guiding energy adjustments at each timestep. Outputs the simulation results in the form of a DataFrame, capturing both current states and decisions made.

# Arguments
- `dt0::DateTime`: The starting date and time for the simulation.
- `strat::AbstractArray{Int}`: An array representing the strategic decisions over time.
- `vehicles_plugged_1::Float64`: The initial number of vehicles that are plugged in.
- `enerfrac_plugged_1::Float64`: The initial energy fraction of plugged vehicles.
- `enerfrac_driving_1::Float64`: The initial energy fraction of vehicles that are driving.
- `stochastic::Bool=false`: Optional argument to enable stochastic simulation of energy needs.

# Returns
- `DataFrame`: A DataFrame containing columns for each timestep including datetime, needed energy fraction, fractions and energy below and above threshold, plugged and driving vehicle energy, energy fraction change (`denerfrac`), and state representation.
"""
function simu_strat(dt0::DateTime, strat::AbstractArray{Int}, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64, stochastic::Bool=false)
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}}[]

    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

    push!(rows, (dt0, enerfrac_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], vehicle_split[3], enerfrac_plugged_1, enerfrac_driving_1, missing, missing))

    for tt in 1:(size(strat)[1]-1)
        denerfrac = make_actions(enerfrac_plugged_1)

        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1))
        ppbase = strat[tt, statebase...]
        ppceil1 = strat[tt, makeindex1(statebase, stateceil1)]
        ppceil2 = strat[tt, makeindex2(statebase, stateceil2)]
        ppceil3 = strat[tt, makeindex3(statebase, stateceil3)]
        denerfrac_pp = ((probbase1 + probbase2 + probbase3) .* denerfrac[ppbase] + (1 .- probbase1) .* denerfrac[ppceil1] + (1 .- probbase2) .* denerfrac[ppceil2] + (1 .- probbase3) .* denerfrac[ppceil3]) / 3;

        enerfrac_plugged_2 = enerfrac_plugged_1 + denerfrac_pp

        enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(enerfrac_plugged_2, enerfrac_needed)

        if stochastic
            simustep = get_simustep_stochastic(dt0 + periodstep(tt))
        else
            simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        end
        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], enerfrac_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt), enerfrac_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], vehicle_split[3], enerfrac_plugged_2, enerfrac_driving_2, denerfrac_pp, statebase))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :enerfrac_needed, :vehicles_plugged, :portion_below, :enerfrac_below, :enerfrac_above, :enerfrac_plugged, :enerfrac_driving, :denerfrac, :state])
end

