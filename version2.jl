## Solution concept:
# Return both standard V(state) and expected profit to aggregator (state)
# Impose different levels of regulation-available range and solve within these
# Choose whichever provides the maximum expected profit in the known initial state
# Current version assumes available always except 8am - 6pm; future versions could impose more complicated timeseries of reg-available

using CSV
include("version1.jl")

energy_min = 0.
energy_max = vehicle_capacity * vehicles

RR = 5 # number of possible regrange values

probfail_penalty = 10.
probfail_limit = 0.01

pricedf = CSV.read("predprice.csv", DataFrame)
pricedf[!, :datetime] = DateTime.(replace.(pricedf.datetime, " " => "T"))

"""
Optimize the cost using Bellman optimization for a stochastic process.

Arguments:
- SS: Planning horizon (timesteps).

Returns:
- strat: SSxEE matrix representing the optimal strategy.

"""
function optimize(dt0::DateTime, probstate::Array{Float64, 4})
    strat = zeros(Int64, SS-1, EE, FF, FF);
    optregrange = zeros(Float64, SS-1);

    # Construct dimensions
    enerfrac_range = [0.; range(enerfrac_min, enerfrac_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    denerfrac_FF = [make_actions(enerfrac_plugged) for enerfrac_plugged=enerfrac_range];
    denerfrac = [denerfrac_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, enerfrac_driving=enerfrac_range];
    energy_denerfrac_byact = [vehicles_plugged_range[ee] * vehicle_capacity * denerfrac[pp] for pp=1:PP, ee=1:EE]

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
    energy_maxregrange_bystate = min.(energy_maxallow .- energy_bystate, energy_bystate - energy_minallow)

    regrange_range = collect(range(0., maximum(energy_maxregrange_bystate), RR))

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
            regrange_good_byact = [(energy_bystate[ee, ff1] .- regrange + energy_denerfrac_byact[pp, ee] .> energy_minallow[ee, ff1]) .& (energy_bystate[ee, ff1] .+ regrange + energy_denerfrac_byact[pp, ee] .< energy_maxallow[ee, ff1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF]

            for mc in 1:mcdraws
                if mcdraws == 1
                    simustep = get_simustep_deterministic(dt1)
                else
                    simustep = get_simustep_stochastic(dt1)
                end
                # Note: We impose costs from enerfrac-below vehicles, but do not adjust state because it pushes up plugged-in enerfrac every period
                statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]), enerfrac1_byaction[pp, ee, ff1, ff2], enerfrac_range[ff2]), vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

                state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2)
                VV1byactthismc = combinebyact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
                VV1byactsummc += VV1byactthismc;

                if regrange > 0
                    ## Calculate probability of falling outside of allowed range
                    probfailthismc = combinebyact(probfail .+ regrange_fail_bystate, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3) .+ regrange_fail_byact .+ (1 .- regrange_good_byact);
                else
                    # no additional probfail
                    probfailthismc = combinebyact(probfail, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
                end

                probfailsummc += probfailthismc
            end

            VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction * hourly_valuee;
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

dt0 = DateTime("2023-07-15T12:00:00")
vehicles_plugged_1 = 4.
enerfrac_plugged_1 = 0.5

probstate = zeros(SS-1, EE, FF, FF);
df = simu_inactive(dt0, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)
for ii in 1:nrow(df)
    statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((df.vehicles_plugged[ii], df.enerfrac_plugged[ii], df.enerfrac_driving[ii]))
    probstate[statebase] = (probbase1 + probbase2 + probbase3) / 3
    probstate[makeindex1.(state2base, state2ceil1)] = (1 .- probbase1) / 3
    probstate[makeindex2.(state2base, state2ceil2)] = (1 .- probbase2) / 3
    probstate[makeindex3.(state2base, state2ceil3)] = (1 .- probbase3) / 3
end

probstate = probstate / 2 .+ repeat(ones(1, EE, FF, FF) / (EE * FF * FF), SS-1) / 2;

mcdraws = 1
strat, probfail, optregrange = optimize(dt0, probstate);
df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, 0.)
df[!, :optregrange] = optregrange
pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label="")
subdf = df[df.optregrange .> 0, :]
plot!(pp, subdf.datetime, subdf.enerfrac_plugged + subdf.optregrange / (vehicles * vehicle_capacity), seriestype=:point)
plot!(pp, subdf.datetime, subdf.enerfrac_plugged - subdf.optregrange / (vehicles * vehicle_capacity), seriestype=:point)
pp

mcdraws = 20

dfall = nothing
for ii in 1:mcdraws
    df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, 0., true)
    energy = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) .* df.enerfrac_above
    energy_minallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.3
    energy_maxallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.95
    df[!, :regrange_avail] = min.(energy_maxallow - energy, energy - energy_minallow)

    if dfall == nothing
        dfall = df
    else
        append!(dfall, df)
    end
end

probstate = zeros(SS-1, EE, FF, FF);
for tt in 1:SS-1
    dt1 = dt0 + periodstep(tt)

    for state in unique(dfall.state[dfall.datetime .== dt1])
        probstate[tt, state...] = sum(map(rowstate -> rowstate == state, dfall.state[dfall.datetime .== dt1])) / sum(dfall.datetime .== dt1)
    end
end
