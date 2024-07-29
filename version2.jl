## Solution concept:
# Return both standard V(state) and expected profit to aggregator (state)
# Impose different levels of regulation-available range and solve within these
# Choose whichever provides the maximum expected profit in the known initial state
# Current version assumes available always except 8am - 6pm; future versions could impose more complicated timeseries of reg-available

include("version1.jl")

energy_min = 0.
energy_max = vehicle_capacity * vehicles

probfail_penalty = 10.
probfail_limit = 0.01

"""
Optimize the cost using Bellman optimization for a stochastic process.

Arguments:
- SS: Planning horizon (timesteps).

Returns:
- strat: SSxEE matrix representing the optimal strategy.

"""
function optimize(dt0::DateTime, SS::Int, regrange::Float64)
    strat = zeros(Int64, SS-1, EE, FF, FF);

    # Construct dimensions
    enerfrac_range = [0.; range(enerfrac_min, enerfrac_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    denerfrac_FF = [make_actions(enerfrac_plugged) for enerfrac_plugged=enerfrac_range];
    denerfrac = [denerfrac_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, enerfrac_driving=enerfrac_range];

    enerfrac0_byaction = repeat(reshape(enerfrac_range, 1, 1, FF, 1), PP, EE, 1, FF);

    # STEP 1: Calculate V[S] under every scenario
    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(SS))
    vehicle_split = split_below.(enerfrac_range, enerfrac_needed)
    value_energy_byenerfrac = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], enerfrac_needed) for ff=1:FF]
    VV2 = repeat(reshape(value_energy_byenerfrac, 1, FF), EE, 1, FF)

    # Determine the energy available for each state
    # Assumes same enerfrac_needed as midnight
    energy_minallow = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * 0.3 for ee=1:EE, ff=1:FF]
    energy_maxallow = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * 0.95 for ee=1:EE, ff=1:FF]
    energy_bystate = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * vehicle_split[ff][3] for ee=1:EE, ff=1:FF]
    regrange_good = (energy_bystate .- regrange .> energy_minallow) .& (energy_bystate .+ regrange .< energy_maxallow)
    regrange_fail_bystate = 1. .- repeat(regrange_good, 1, 1, FF);
    regrange_fail_byact = repeat(reshape(regrange_fail_bystate, 1, EE, FF, FF), PP);

    probfail = zeros(Float64, EE, FF, FF); # Sum over periods looking forward

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        enerfrac1_byaction = enerfrac0_byaction .+ denerfrac;

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = value_power_action.(price, denerfrac);

        enerfrac_needed = enerfrac_scheduled(dt1)
        vehicle_split = split_below.(enerfrac_range, enerfrac_needed)
        valuepns = [value_power_newstate.(price, vehicle_split[ff12][1], enerfrac_needed - vehicle_split[ff12][2]) for ff12 in 1:FF]
        valuee = [value_energy.(vehicle_split[ff12][1], vehicle_split[ff12][3], enerfrac_needed) for ff12 in 1:FF];

        ff12_byaction = discrete_roundbelow.(enerfrac1_byaction, enerfrac_min, enerfrac_max, FF);
        valuepns_byaction = valuepns[ff12_byaction];
        valuee_byaction = valuee[ff12_byaction];

        # Check if we are on the market
        date_part = Dates.Date(dt1)
        if date_part != Dates.Date(dt0)
            dt_8am = DateTime(date_part, Dates.Time(8, 0, 0))
            dt_6pm = DateTime(date_part, Dates.Time(18, 0, 0))
            if dt1 < dt_8am || dt1 ≥ dt_6pm
                onmarket = true
            else
                onmarket = false
            end
        else
            onmarket = false
        end

        VV1byactsummc = zeros(Float64, PP, EE, FF, FF);
        probfailsummc = zeros(Float64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1)
            else
                simustep = get_simustep_stochastic(dt1)
            end
            state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 =
                simustate2(simustep, vehicles_plugged_range, vehicle_split, ff12_byaction, enerfrac1_byaction, enerfrac_range);
            VV1byactthismc = calcVV1byact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
            VV1byactsummc += VV1byactthismc;

            if onmarket
                ## Calculate probability of falling outside of allowed range
                probfailthismc = calcVV1byact(probfail .+ regrange_fail_bystate, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3) .+ regrange_fail_byact;
            else
                # no additional probfail
                probfailthismc = calcVV1byact(probfail, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
            end
            probfailsummc += probfailthismc
        end

        VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction * hourly_valuee;
        VV1byact[isnan.(VV1byact)] .= -Inf

        probfailbyact = probfailsummc / mcdraws

        bestact = dropdims(argmax(VV1byact - probfail_penalty * probfailbyact, dims=1), dims=1);
        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);

        VV2 = VV1byact[bestact]
        VV2[isnan.(VV2)] .= -Inf

        probfail = probfailbyact[bestact]
    end

    return strat, probfail
end

dt0 = DateTime("2024-07-15T12:00:00")
vehicles_plugged_1 = 4.
enerfrac_plugged_1 = 0.5

##for regrange in range(energy_min, energy_max * (.90 - .3) / 2, length=10)

regrange = (energy_max * (.90 - .3) / 2) / 2
strat, probfail = optimize(dt0, SS, regrange);

mcdraws = 20

dfall = nothing
for ii in 1:mcdraws
    df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, 0., true)
    energy = vehicle_capacity * df.vehicles_plugged * (1 - df.portion_below) * df.enerfrac_above
    energy_minallow = vehicle_capacity * df.vehicles_plugged * (1 - df.portion_below) * 0.3
    energy_maxallow = vehicle_capacity * df.vehicles_plugged * (1 - df.portion_below) * 0.95
    df[!, :regrange_avail] = min.(energy_maxallow - energy, energy - energy_minallow)

    if dfall == nothing
        dfall = df
    else
        extend!(dfall, df)
    end
end

## TODO: Calculate available regrange under probfail_limit
dfall %>% group_by(datetime) %>% summarize(regrange=quantile(regrange, 1 - probfail_limit))
## TODO: replace regrange by vector, and cue onmarket by > 0
## TODO: random gradient search to find optimum
