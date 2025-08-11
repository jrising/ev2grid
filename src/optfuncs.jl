include("bizutils.jl")
include("value.jl")
include("retail.jl")
include("optutils.jl")
include("simulate.jl")

function optimize(dt0::DateTime, SS::Int, drive_starts_time::Time, park_starts_time::Time, hours_sd::Float64=0.)
    global event_log = [] ## clear out event_log before optimizing
    ## vehicles_plugged, soc_plugged, soc_driving
    strat = zeros(Int64, SS-1, EE, FF, FF);
    VVall = zeros(Float64, SS, EE, FF, FF);

    # Construct dimensions
    soc_range = [0.; range(soc_min, soc_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    dsoc_FF = [make_actions(soc_plugged, soc_range) for soc_plugged=soc_range];
    dsoc = [dsoc_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, soc_driving=soc_range];

    soc0_byaction = repeat(reshape(soc_range, 1, 1, FF, 1), PP, EE, 1, FF);
    soc1_byaction = soc0_byaction .+ dsoc;

    # STEP 1: Calculate V[S] under every scenario
    soc_needed = soc_scheduled(dt0 + periodstep(SS), drive_starts_time)
    vehicle_split = split_below.(soc_range, soc_needed)
    value_energy_bysoc = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff=1:FF]
    VV2 = repeat(reshape(value_energy_bysoc, EE, FF, 1), 1, 1, FF)
    VVall[end, :, :, :] = VV2

    # Take draws of drive and start times, if that's stochastic
    if hours_sd > 0
        drive_times = get_stochastic_times(mcdraws, drive_starts_time, hours_sd, timestep)
        park_times = get_stochastic_times(mcdraws, park_starts_time, hours_sd, timestep)
    end

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = [value_power_action(price, dsoc[pp, ee, ff1, ff2], vehicle_split[ff1][1], vehicles_plugged_range[ee]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        soc_needed = soc_scheduled(dt1, drive_starts_time)
        vehicle_split = split_below.(soc_range, soc_needed);
        valuepns = [value_power_newstate(price, vehicle_split[ff12][1], soc_needed - vehicle_split[ff12][2], vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF];
        vehicle_split = split_below.(soc_range, soc_needed);
        valuee = [value_energy(vehicle_split[ff12][1], vehicle_split[ff12][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF];

        ff12_byaction = discrete_roundbelow.(soc1_byaction, soc_min, soc_max, FF);
        valuepns_byaction = [valuepns[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];
        valuee_byaction = [valuee[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        VV1byactsummc = zeros(Float64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
            elseif hours_sd > 0
                simustep = get_simustep_deterministic(dt1, drive_times[mc], park_times[mc])
            else
                simustep = get_simustep_stochastic(dt1, drive_starts_time, park_starts_time)
            end
            # Note: We impose costs from soc-below vehicles, but do not adjust state because it pushes up plugged-in soc every period
            statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]),
                                               vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][3], soc_range[ff2]),
                                      vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

            state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2);

            VV1byactthismc = combinebyact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3)
            VV1byactsummc += VV1byactthismc;
        end

        VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction;
        VV1byact[isnan.(VV1byact)] .= -Inf

        bestact = dropdims(argmax(VV1byact, dims=1), dims=1);
        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);
        if any(strat[tt, :, :, :] .== 0)
            break
        end

        VV2 = VV1byact[bestact]
        VV2[isnan.(VV2)] .= -Inf
        VVall[tt, :, :, :] = VV2
    end

    return strat, VVall
end

function optimize_regrange_probstate(dt0::DateTime, probstate::Array{Float64, 4},  drive_starts_time::Time, park_starts_time::Time,)
    """
    Optimize the cost using Bellman optimization for a stochastic process.

    Arguments:
    - SS: Planning horizon (timesteps).

    Returns:
    - strat: SSxEE matrix representing the optimal strategy.

    """

    strat = zeros(Int64, SS-1, EE, FF, FF);
    optregrange = zeros(Float64, SS-1);

    # Construct dimensions
    soc_range = [0.; range(soc_min, soc_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    dsoc_FF = [make_actions(soc_plugged, soc_range) for soc_plugged=soc_range];
    dsoc = [dsoc_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, soc_driving=soc_range];
    energy_dsoc_byact = [vehicles_plugged_range[ee] * vehicle_capacity * dsoc[pp, ee, ff1, ff2] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF]

    soc0_byaction = repeat(reshape(soc_range, 1, 1, FF, 1), PP, EE, 1, FF);
    soc1_byaction = soc0_byaction .+ dsoc;

    # STEP 1: Calculate V[S] under every scenario
    soc_needed = soc_scheduled(dt0 + periodstep(SS), drive_starts_time)
    vehicle_split = split_below.(soc_range, soc_needed)
    value_energy_bysoc = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff=1:FF]
    VV2 = repeat(reshape(value_energy_bysoc, EE, FF, 1), 1, 1, FF)

    # Determine the energy available for each state
    # Assumes same soc_needed as midnight
    energy_minallow = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * 0.3 for ee=1:EE, ff=1:FF]
    energy_maxallow = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * 0.95 for ee=1:EE, ff=1:FF]
    energy_bystate = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * vehicle_split[ff][3] for ee=1:EE, ff=1:FF]
    energy_maxregrange_bystate = min.(energy_maxallow .- energy_bystate, energy_bystate - energy_minallow)

    regrange_range = collect(range(0., maximum(energy_maxregrange_bystate), RR))

    probfail = zeros(Float64, EE, FF, FF); # Sum over periods looking forward

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = [value_power_action(price, dsoc[pp, ee, ff1, ff2], vehicle_split[ff1][1], vehicles_plugged_range[ee]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF]

        soc_needed = soc_scheduled(dt1, drive_starts_time)
        vehicle_split = split_below.(soc_range, soc_needed)
        valuepns = [value_power_newstate(price, vehicle_split[ff12][1], soc_needed - vehicle_split[ff12][2], vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF]
        valuee = [value_energy.(vehicle_split[ff12][1], vehicle_split[ff12][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF];

        ff12_byaction = discrete_roundbelow.(soc1_byaction, soc_min, soc_max, FF);
        valuepns_byaction = [valuepns[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];
        valuee_byaction = [valuee[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        pricedfrow = pricedf[pricedf.datetime .== dt1, :] # only works if timestep is whole hours
        regprice = pricedfrow.predpe[1] # XXX: Later use uncertainty

        bestregrange = -1.
        beststrat = zeros(Int64, EE, FF, FF);
        bestVV2 = zeros(Float64, EE, FF, FF);
        bestprobfail = zeros(Float64, EE, FF, FF);
        bestvalue = -Inf

        for regrange in regrange_range
            VV1byactsummc = zeros(Float64, PP, EE, FF, FF);
            probfailsummc = zeros(Float64, PP, EE, FF, FF);

            regrange_good = (energy_bystate .- regrange .> energy_minallow) .& (energy_bystate .+ regrange .< energy_maxallow)
            regrange_fail_bystate = 1. .- repeat(regrange_good, 1, 1, FF);
            regrange_fail_byact = repeat(reshape(regrange_fail_bystate, 1, EE, FF, FF), PP);
            # Also disallow actions that would overextend our total charge range
            regrange_good_byact = [(energy_bystate[ee, ff1] .- regrange + energy_dsoc_byact[pp, ee, ff1, ff2] .> energy_minallow[ee, ff1]) .& (energy_bystate[ee, ff1] .+ regrange + energy_dsoc_byact[pp, ee, ff1, ff2] .< energy_maxallow[ee, ff1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF]

            for mc in 1:mcdraws
                if mcdraws == 1
                    simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
                else
                    simustep = get_simustep_stochastic(dt1, drive_starts_time, park_starts_time)
                end
                # Note: We impose costs from soc-below vehicles, but do not adjust state because it pushes up plugged-in soc every period
                statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]),
                                                   vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][3], soc_range[ff2]),
                                          vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

                state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2);
                VV1byactthismc = combinebyact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
                VV1byactsummc += VV1byactthismc;

                if regrange > 0
                    ## Calculate probability of falling outside of allowed range
                    probfailthismc = combinebyact(probfail .+ regrange_fail_bystate, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3) .+ (1 .- regrange_good_byact);
                else
                    # no additional probfail
                    probfailthismc = combinebyact(probfail, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
                end

                probfailsummc += probfailthismc
            end

            VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction;
            VV1byact[isnan.(VV1byact)] .= -Inf

            probfailbyact = probfailsummc / mcdraws

            bestact = dropdims(argmax(VV1byact - probfail_penalty * probfailbyact, dims=1), dims=1);

            thisVV2 = VV1byact[bestact];
            thisprobfail = probfailbyact[bestact]
            thisvalue = sum(probstate[tt, :, :, :] .* (thisVV2 .+ regprice * regrange * (1 .- thisprobfail)))
            if thisvalue > bestvalue
                bestvalue = thisvalue
                beststrat = Base.Fix2(getindex, 1).(bestact);
                bestVV2 = thisVV2;
                bestprobfail = thisprobfail
                bestregrange = regrange
            end
        end

        optregrange[tt] = bestregrange
        strat[tt, :, :, :] .= beststrat
        VV2 = bestVV2
        probfail = bestprobfail
    end

    return strat, probfail, optregrange
end

function optimize_regrange_probstate_outer_loop(dt0, soc_plugged_1, soc_driving_1, vehicles_plugged_1, drive_starts_time, park_starts_time)
    global probstate = zeros(SS-1, EE, FF, FF);
    df = fullsimulate(dt0, zeros(SS-1), vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
    for ii in 1:(nrow(df) - 1)
        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((df.vehicles_plugged[ii], df.soc_plugged[ii], df.soc_driving[ii]))
        probstate[ii, statebase...] = (probbase1 + probbase2 + probbase3) / 3
        probstate[ii, makeindex1(statebase, stateceil1)] += (1 - probbase1) / 3
        probstate[ii, makeindex2(statebase, stateceil2)] += (1 - probbase2) / 3
        probstate[ii, makeindex3(statebase, stateceil3)] += (1 - probbase3) / 3
    end

    probstate = probstate / 2 .+ repeat(ones(1, EE, FF, FF) / (EE * FF * FF), SS-1) / 2;

    for ll in 1:10
        println("Loop $(ll)")
        strat, probfail, optregrange = optimize_regrange_probstate(dt0, probstate, drive_starts_time, park_starts_time);

        dfall = nothing
        for ii in 1:mcdraws
            local df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
            energy = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) .* df.soc_above
            energy_minallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.3
            energy_maxallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.95
            df[!, :regrange_avail] = min.(energy_maxallow - energy, energy - energy_minallow)

            if dfall == nothing
                dfall = df
            else
                append!(dfall, df)
            end
        end

        global probstate = zeros(SS-1, EE, FF, FF);
        for tt in 1:SS-1
            dt1 = dt0 + periodstep(tt - 1)

            for state in unique(dfall.state[dfall.datetime .== dt1])
                probstate[tt, state...] = sum(map(rowstate -> rowstate == state, dfall.state[dfall.datetime .== dt1])) / sum(dfall.datetime .== dt1)
            end
        end
    end
    return probstate
end

