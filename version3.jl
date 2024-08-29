## Solution concept:
# Return both standard V(state) and expected profit to aggregator (state)
# Impose different levels of regulation-available range and solve within these
# Choose whichever provides the maximum expected profit in the known initial state
# Current version assumes available always except 8am - 6pm; future versions could impose more complicated timeseries of reg-available

using Random

include("src/bizutils.jl")
include("src/customer.jl")
include("src/simulate.jl")
include("src/retail.jl")
include("src/config.jl")
include("src/value.jl")
include("src/optutils.jl")
include("src/fullsim.jl")
include("src/plotting.jl")

energy_min = 0.
energy_max = vehicle_capacity * vehicles

probfail_penalty = 100.
portion_below_penalty = 1000.

"""
Optimize the cost using Bellman optimization for a stochastic process.

Arguments:
- SS: Planning horizon (timesteps).

Returns:
- strat: SSxEE matrix representing the optimal strategy.

"""
function optimize(dt0::DateTime, regrange::Vector{Float64})
    strat = zeros(Int64, SS-1, EE, FF, FF);

    # Construct dimensions
    soc_range = [0.; range(soc_min, soc_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    dsoc_FF = [make_actions(soc_plugged) for soc_plugged=soc_range];
    dsoc = [dsoc_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, soc_driving=soc_range];
    energy_dsoc_byact = [vehicles_plugged_range[ee] * vehicle_capacity * dsoc[pp] for pp=1:PP, ee=1:EE]

    soc0_byaction = repeat(reshape(soc_range, 1, 1, FF, 1), PP, EE, 1, FF);
    soc1_byaction = soc0_byaction .+ dsoc;

    # STEP 1: Calculate V[S] under every scenario
    soc_needed = soc_scheduled(dt0 + periodstep(SS));
    vehicle_split = split_below.(soc_range, soc_needed);
    value_energy_bysoc = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], soc_needed) for ff=1:FF];
    VV2 = repeat(reshape(value_energy_bysoc, 1, FF), EE, 1, FF);

    # Determine the energy available for each state
    # Assumes same soc_needed as midnight
    energy_minallow = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * 0.3 for ee=1:EE, ff=1:FF];
    energy_maxallow = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * 0.95 for ee=1:EE, ff=1:FF];
    energy_bystate = [vehicle_capacity * vehicles_plugged_range[ee] * (1. - vehicle_split[ff][1]) * vehicle_split[ff][3] for ee=1:EE, ff=1:FF];
    energy_maxregrange_bystate = min.(energy_maxallow .- energy_bystate, energy_bystate - energy_minallow);

    probfail = zeros(Float64, EE, FF, FF); # Sum over periods looking forward

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        # println(tt)

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = value_power_action.(price, dsoc);

        soc_needed = soc_scheduled(dt1)
        vehicle_split = split_below.(soc_range, soc_needed)
        valuepns = [value_power_newstate.(price, vehicle_split[ff12][1], soc_needed - vehicle_split[ff12][2]) for ff12 in 1:FF]
        valuee = [value_energy.(vehicle_split[ff12][1], vehicle_split[ff12][3], soc_needed) for ff12 in 1:FF];

        ff12_byaction = discrete_roundbelow.(soc1_byaction, soc_min, soc_max, FF);
        valuepns_byaction = valuepns[ff12_byaction];
        valuee_byaction = valuee[ff12_byaction];
        portion_below_byaction = [vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        pricedfrow = pricedf[pricedf.datetime .== dt1, :] # only works if timestep is whole hours
        regprice = pricedfrow.predpe[1] # XXX: Later use uncertainty

        regrange_good = (energy_bystate .- regrange[tt] .> energy_minallow) .& (energy_bystate .+ regrange[tt] .< energy_maxallow)
        regrange_fail_bystate = 1. .- repeat(regrange_good, 1, 1, FF);
        regrange_fail_byact = repeat(reshape(regrange_fail_bystate, 1, EE, FF, FF), PP);
        # Also disallow actions that would overextend our total charge range
        regrange_good_byact = [(energy_bystate[ee, ff1] .- regrange[tt] + energy_dsoc_byact[pp, ee] .> energy_minallow[ee, ff1]) .& (energy_bystate[ee, ff1] .+ regrange[tt] + energy_dsoc_byact[pp, ee] .< energy_maxallow[ee, ff1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        VV1byactsummc = zeros(Float64, PP, EE, FF, FF);
        probfailsummc = zeros(Float64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1)
            else
                simustep = get_simustep_stochastic(dt1)
            end
            # Note: We impose costs from soc-below vehicles, but do not adjust state because it pushes up plugged-in soc every period
            statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]), soc1_byaction[pp, ee, ff1, ff2], soc_range[ff2]), vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

            state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2);
            VV1byactthismc = combinebyact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
            VV1byactsummc += VV1byactthismc;

            if regrange[tt] > 1 # can have very small values that result in failures
                ## Calculate probability of falling outside of allowed range
                probfailthismc = combinebyact(probfail .+ regrange_fail_bystate, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3) .+ (1 .- regrange_good_byact);
            else
                # no additional probfail
                probfailthismc = combinebyact(probfail, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3);
            end

            probfailsummc += probfailthismc
        end

        probfailbyact = probfailsummc / mcdraws;
        VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction + regprice * regrange[tt] * (1 .- probfailbyact);
        VV1byact[isnan.(VV1byact)] .= -Inf;

        bestact = dropdims(argmax(VV1byact - probfail_penalty * sqrt.(probfailbyact) - portion_below_penalty * sqrt.(portion_below_byaction), dims=1), dims=1);

        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);
        VV2 = VV1byact[bestact];
        probfail = probfailbyact[bestact];
    end

    return strat, probfail
end

dt0 = DateTime("2023-07-15T12:00:00")
vehicles_plugged_1 = 0.
soc_plugged_1 = 0.5
soc_driving_1 = 0.5
mcdraws = 20

function calcvalue(df, regrange)
    df[!, :kw] = (df.dsoc .* vehicle_capacity .* df.vehicles_plugged .* (1 .- df.portion_below) / timestep) .+
        (max_charging_kw * df.portion_below .* df.vehicles_plugged) .+ regrange
    sum(df.valuep .+ df.valuepns .+ df.valuee .+ df.valuer .- portion_below_penalty * sqrt.(df.portion_below)) -
        get_demand_cost(dt0, collect(skipmissing(df.kw)), timestep)
end

function calcmaxrange!(df)
    df[!, :energy_soc] = vehicle_capacity .* df.vehicles_plugged .* (1 .- df.portion_below) .* df.soc_above
    df[!, :energy_minallow] = vehicle_capacity .* df.vehicles_plugged .* (1 .- df.portion_below) * 0.3
    df[!, :energy_maxallow] = vehicle_capacity .* df.vehicles_plugged .* (1 .- df.portion_below) * 0.95
    df[!, :regrange_avail] = min.(df.energy_maxallow .- df.energy_soc, df.energy_soc .- df.energy_minallow)
    df[!, :regrange_maxrange] = ((df.energy_maxallow .- df.energy_minallow) / 2) .* (1. .- disallows)
end

df = fullsimulate(dt0, zeros(SS-1), vehicles_plugged_1, 0.8, 0.8) # Make sure all can drive
# Don't allow regrange when vehicles return, because no time to recharge if very low energy
# Don't allow in period that vehicles leave, because they aren't there by end of period
disallows = [false; df.vehicles_plugged[2:end] - df.vehicles_plugged[1:end-1] .> 1.0] .| [df.vehicles_plugged[2:end] - df.vehicles_plugged[1:end-1] .< -1.0; false];

calcmaxrange!(df)

regrange = df.regrange_avail .* (1. .- disallows)

df = fullsimulate(dt0, ones(Int, SS-1, EE, FF, FF), regrange, vehicles_plugged_1, 0.8, 0.8)
calcmaxrange!(df)

valuebest = calcvalue(df, regrange)
dfbest = df
regrangebest = regrange
stratbest = ones(Int, SS-1, EE, FF, FF);
lastimprove = 0
for ll in 1:10000
    println("Loop $(ll)")
    Random.seed!(20240826);
    strat, probfail = optimize(dt0, regrange);

    df = fullsimulate(dt0, strat, regrange, vehicles_plugged_1, soc_plugged_1, soc_driving_1, mcdraws > 1);
    value = calcvalue(df, regrange)

    if value > valuebest
        println("Improved.")
        valuebest = value
        regrangebest = regrange
        stratbest = strat
        lastimprove = ll

        calcmaxrange!(df)

        dfbest = df
    elseif ll >= 100 && (ll - lastimprove > 10)
        break
    end

    regrange = rand.(Truncated.(Normal.(min.(dfbest.regrange_maxrange, regrangebest), exp(-ll / 100)),
                                0., dfbest.regrange_maxrange))
end

df = fullsimulate(dt0, stratbest, regrangebest, vehicles_plugged_1, soc_plugged_1, 0.)
df[!, :optregrange] = [0.; regrangebest[1:end-1]]

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=df.optregrange / (vehicles * vehicle_capacity), label="Regulation Range")
plot!(size=(700,400))
# savefig("optregrange.pdf")
