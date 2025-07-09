using Dates

"""
Simulates the energy fraction of vehicles while applying a given changes in energy over several timesteps. The strategy is provided as an array, guiding energy adjustments at each timestep. Outputs the simulation results in the form of a DataFrame, capturing both current states and decisions made.

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
function fullsimulate(dt0::DateTime, get_dsoc::Function, get_regrange::Function, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, drive_starts_time::Time, park_starts_time::Time, stochastic::Bool=false)
    global event_log = [] ## clear out event_log before simulating
    ## vehicles_plugged_1, soc_plugged_1, soc_driving_1 = 0., 0.5, 0.5 # debugging
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}, Float64, Float64, Float64, Float64}[]

    soc_needed = soc_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(soc_plugged_1, soc_needed)

    for tt in 1:(SS-1)
        dsoc = get_dsoc(tt, (vehicles_plugged_1, soc_plugged_1, soc_driving_1))
        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((vehicles_plugged_1, soc_plugged_1, soc_driving_1))

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)

        valuep = value_power_action(price, dsoc, vehicle_split[1], vehicles_plugged_1)
        valuepns = value_power_newstate(price, vehicle_split[1], soc_needed - vehicle_split[2], vehicles_plugged_1)
        valuee = value_energy(vehicle_split[1], vehicle_split[3], soc_needed, vehicles_plugged_1)

        pricedfrow = pricedf[pricedf.datetime .== dt1, :] # only works if timestep is whole hours
        regprice = pricedfrow.predpe[1]
        valuer = regprice * get_regrange(tt)

        push!(rows, (dt0 + periodstep(tt), soc_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], vehicle_split[3], soc_plugged_1, soc_driving_1, dsoc, statebase, valuep, valuepns, valuee, valuer))

        ## Apply action
        soc_plugged_2 = soc_plugged_1 + dsoc

        ## Apply simulation
        if stochastic
            simustep = get_simustep_stochastic(dt1, drive_starts_time, park_starts_time)
        else
            simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
        end
        soc_needed = soc_scheduled(dt1)
        vehicle_split = split_below(soc_plugged_2, soc_needed)
        vehicles_plugged_2, soc_plugged_2, soc_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), soc_plugged_2, soc_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        vehicles_plugged_1, soc_plugged_1, soc_driving_1 = vehicles_plugged_2, soc_plugged_2, soc_driving_2
        vehicle_split = split_below(soc_plugged_2, soc_needed)
    end

    push!(rows, (dt0 + periodstep(SS), soc_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], vehicle_split[3], soc_plugged_1, soc_driving_1, missing, missing, 0., 0., 0., 0.))

    df = DataFrame(rows)
    rename!(df, [:datetime, :soc_needed, :vehicles_plugged, :portion_below, :soc_below, :soc_above, :soc_plugged, :soc_driving, :dsoc, :state, :valuep, :valuepns, :valuee, :valuer])
end

function fullsimulate(dt0::DateTime, strat::AbstractArray{Int}, regrange::Vector{Float64}, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, drive_starts_time::Time, park_starts_time::Time, stochastic::Bool=false)
    function get_dsoc(tt, state)
        # Determine action for this period
        soc_range = [0.; range(soc_min, soc_max, FF-1)];

        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate(state)
        dsoc_base = make_actions(soc_range[statebase[2]], soc_range)
        dsoc_ceil = make_actions(soc_range[stateceil2], soc_range)

        ppbase = strat[tt, statebase...]
        ppceil1 = strat[tt, makeindex1(statebase, stateceil1)]
        ppceil2 = strat[tt, makeindex2(statebase, stateceil2)]
        ppceil3 = strat[tt, makeindex3(statebase, stateceil3)]
        dsoc_pp = ((probbase1 + probbase2 + probbase3) .* dsoc_base[ppbase] + (1 .- probbase1) .* dsoc_ceil[ppceil1] + (1 .- probbase2) .* dsoc_ceil[ppceil2] + (1 .- probbase3) .* dsoc_ceil[ppceil3]) / 3;

        dsoc_pp
    end
    fullsimulate(dt0, get_dsoc, tt -> regrange[tt], vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time, stochastic)
end

function fullsimulate(dt0::DateTime, dsoc::Vector{Float64}, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, drive_starts_time::Time, park_starts_time::Time, stochastic::Bool=false)
    fullsimulate(dt0, (tt, state) -> dsoc[tt], (tt) -> 0., vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time, stochastic)
end

function fullsimulate_modify(dt0::DateTime, strat::AbstractArray{Int}, changes::Dict{Int64, Float64}, regrange::Vector{Float64}, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, stochastic::Bool=false)
    function get_dsoc(tt, state)
        if tt ∈ keys(changes)
            return changes[tt]
        end
        # Determine action for this period
        soc_range = [0.; range(soc_min, soc_max, FF-1)];
        dsoc = make_actions(state[2], soc_range)

        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate(state)

        ppbase = strat[tt, statebase...]
        ppceil1 = strat[tt, makeindex1(statebase, stateceil1)]
        ppceil2 = strat[tt, makeindex2(statebase, stateceil2)]
        ppceil3 = strat[tt, makeindex3(statebase, stateceil3)]
        dsoc_pp = ((probbase1 + probbase2 + probbase3) .* dsoc[ppbase] + (1 .- probbase1) .* dsoc[ppceil1] + (1 .- probbase2) .* dsoc[ppceil2] + (1 .- probbase3) .* dsoc[ppceil3]) / 3;

        dsoc_pp
    end
    fullsimulate(dt0, get_dsoc, tt -> regrange[tt], vehicles_plugged_1, soc_plugged_1, soc_driving_1, stochastic)
end

function fullsimulate_with_events(dt0::DateTime, get_dsoc::Function, get_regrange::Function, vehicles_plugged_1::Float64, soc_plugged_1::Float64, soc_driving_1::Float64, drive_starts_time::Time, park_starts_time::Time)
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}, Float64, Float64, Float64, Float64}[]

    soc_needed = soc_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(soc_plugged_1, soc_needed)

    for tt in 1:(SS-1)
        dsoc = get_dsoc(tt, (vehicles_plugged_1, soc_plugged_1, soc_driving_1))
        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((vehicles_plugged_1, soc_plugged_1, soc_driving_1))

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = value_power_action(price, dsoc, vehicles_plugged_1)
        valuepns = value_power_newstate(price, vehicle_split[1], soc_needed - vehicle_split[2], vehicles_plugged_1)
        valuee = value_energy(vehicle_split[1], vehicle_split[3], soc_needed, vehicles_plugged_1)

        pricedfrow = pricedf[pricedf.datetime .== dt1, :]
        regprice = pricedfrow.predpe[1]
        valuer = regprice * get_regrange(tt)

        push!(rows, (dt1, soc_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], vehicle_split[3], soc_plugged_1, soc_driving_1, dsoc, statebase, valuep, valuepns, valuee, valuer))

        soc_plugged_2 = soc_plugged_1 + dsoc
        ev = findfirst(e -> e.time == dt1, event_log)
        ev = ev !== nothing ? event_log[ev] : nothing

        if ev !== nothing
            if ev.event == :allplug
                simustep = simustep_allplug
            elseif ev.event == :alldrive 
                simustep = simustep_alldrive
            elseif ev.event == :emergency
                simustep = (v_p, v_a, s_a, s_d) -> simustep_event(ev.vehicles_needed, v_p, v_a, s_a, s_d)
            elseif ev.event == :delayed_return
                park_starts_time += Minute(60)
                simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
            else
                simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
            end
        else
            simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
        end

        soc_needed = soc_scheduled(dt1)
        vehicle_split = split_below(soc_plugged_2, soc_needed)
        vehicles_plugged_2, soc_plugged_2, soc_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), soc_plugged_2, soc_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        vehicles_plugged_1, soc_plugged_1, soc_driving_1 = vehicles_plugged_2, soc_plugged_2, soc_driving_2
        vehicle_split = split_below(soc_plugged_2, soc_needed)
    end

    push!(rows, (dt0 + periodstep(SS), soc_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], vehicle_split[3], soc_plugged_1, soc_driving_1, missing, missing, 0., 0., 0., 0.))

    df = DataFrame(rows)
    rename!(df, [:datetime, :soc_needed, :vehicles_plugged, :portion_below, :soc_below, :soc_above, :soc_plugged, :soc_driving, :dsoc, :state, :valuep, :valuepns, :valuee, :valuer])
end

