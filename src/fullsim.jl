using Dates

"""
Simulates the energy fraction of vehicles without executing any strategic adjustments over `SS` timesteps. This function models the evolution of energy fractions for parked and driving vehicles, beginning from an initial date and state, allows optional stochastic adjustments, and outputs the simulation results over time in the form of a DataFrame.

# Arguments
- `dt0::DateTime`: The starting date and time for the simulation.
- `vehicles_plugged_1::Float64`: The initial number of vehicles that are plugged in.
- `soc_plugged_1::Float64`: The initial energy fraction of plugged vehicles.
- `soc_driving_1::Float64`: The initial energy fraction of vehicles that are driving.
- `stochastic::Bool=false`: Optional argument to enable stochastic simulation of energy needs.

# Returns
- `DataFrame`: A DataFrame containing columns for each timestep, including datetime, needed energy fraction, fraction and energy below threshold, plugged vehicle energy, and driving vehicle energy.
"""
function simu_inactive(dt0::DateTime, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, stochastic::Bool=false)
    # datetime, soc_needed, vehicles_plugged_1, portion_below, soc_toadd_below, soc_avail_1, soc_driving_1
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64}[]

    soc_needed = soc_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(soc_plugged_1, soc_needed)

    push!(rows, (dt0 + periodstep(1), soc_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], soc_plugged_1, soc_driving_1))

    for tt in 1:(SS-1)
        println(tt)
        soc_needed = soc_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(soc_plugged_1, soc_needed)

        if stochastic
            simustep = get_simustep_stochastic(dt0 + periodstep(tt))
        else
            simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        end

        vehicles_plugged_2, soc_plugged_2, soc_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], soc_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt + 1), soc_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], soc_plugged_2, soc_driving_2))
        vehicles_plugged_1, soc_plugged_1, soc_driving_1 = vehicles_plugged_2, soc_plugged_2, soc_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :soc_needed, :vehicles_plugged, :portion_below, :soc_below, :soc_plugged, :soc_driving])
end

"""
Simulates the energy fraction of vehicles while applying a given strategy over several timesteps. The strategy is provided as an array, guiding energy adjustments at each timestep. Outputs the simulation results in the form of a DataFrame, capturing both current states and decisions made.

# Arguments
- `dt0::DateTime`: The starting date and time for the simulation.
- `strat::AbstractArray{Int}`: An array representing the strategic decisions over time.
- `vehicles_plugged_1::Float64`: The initial number of vehicles that are plugged in.
- `soc_plugged_1::Float64`: The initial energy fraction of plugged vehicles.
- `soc_driving_1::Float64`: The initial energy fraction of vehicles that are driving.
- `stochastic::Bool=false`: Optional argument to enable stochastic simulation of energy needs.

# Returns
- `DataFrame`: A DataFrame containing columns for each timestep including datetime, needed energy fraction, fractions and energy below and above threshold, plugged and driving vehicle energy, energy fraction change (`dsoc`), and state representation.
"""
function simu_strat(dt0::DateTime, strat::AbstractArray{Int}, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, stochastic::Bool=false)
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}}[]

    soc_needed = soc_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(soc_plugged_1, soc_needed)

    push!(rows, (dt0, soc_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], vehicle_split[3], soc_plugged_1, soc_driving_1, missing, missing))

    for tt in 1:(size(strat)[1]-1)
        dsoc = make_actions(soc_plugged_1)

        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((vehicles_plugged_1, soc_plugged_1, soc_driving_1))
        ppbase = strat[tt, statebase...]
        ppceil1 = strat[tt, makeindex1(statebase, stateceil1)]
        ppceil2 = strat[tt, makeindex2(statebase, stateceil2)]
        ppceil3 = strat[tt, makeindex3(statebase, stateceil3)]
        dsoc_pp = ((probbase1 + probbase2 + probbase3) .* dsoc[ppbase] + (1 .- probbase1) .* dsoc[ppceil1] + (1 .- probbase2) .* dsoc[ppceil2] + (1 .- probbase3) .* dsoc[ppceil3]) / 3;

        soc_plugged_2 = soc_plugged_1 + dsoc_pp

        soc_needed = soc_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(soc_plugged_2, soc_needed)

        if stochastic
            simustep = get_simustep_stochastic(dt0 + periodstep(tt))
        else
            simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        end
        vehicles_plugged_2, soc_plugged_2, soc_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], soc_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt), soc_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], vehicle_split[3], soc_plugged_2, soc_driving_2, dsoc_pp, statebase))
        vehicles_plugged_1, soc_plugged_1, soc_driving_1 = vehicles_plugged_2, soc_plugged_2, soc_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :soc_needed, :vehicles_plugged, :portion_below, :soc_below, :soc_above, :soc_plugged, :soc_driving, :dsoc, :state])
end

