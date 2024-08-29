using Dates, DataFrames
using StatsBase, Random
using Plots
using HolidayCalendars, RDates
using ArgCheck

# Assume this is the day before, and we can observe current storage
# If there are any cars that need recharging to 30% level, handle that outside this system

## Order of events:
## In state1 at start of timestep
## Apply action throughout timestep
## Construct below-split and perform simulation
## Transition to state2 at end of timestep

include("src/bizutils.jl")
include("src/customer.jl")
include("src/simulate.jl")
include("src/retail.jl")
include("src/config.jl")
include("src/value.jl")
include("src/optutils.jl")
include("src/fullsim.jl")
include("src/plotting.jl")

"""
Optimize the cost using Bellman optimization for a stochastic process.

Arguments:
- SS: Planning horizon (timesteps).

Returns:
- strat: SSxEE matrix representing the optimal strategy.

"""
function optimize(dt0::DateTime, SS::Int)
    strat = zeros(Int64, SS-1, EE, FF, FF);

    # Construct dimensions
    soc_range = [0.; range(soc_min, soc_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    dsoc_FF = [make_actions(soc_plugged) for soc_plugged=soc_range];
    dsoc = [dsoc_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, soc_driving=soc_range];

    soc0_byaction = repeat(reshape(soc_range, 1, 1, FF, 1), PP, EE, 1, FF);
    soc1_byaction = soc0_byaction .+ dsoc;

    # STEP 1: Calculate V[S] under every scenario
    soc_needed = soc_scheduled(dt0 + periodstep(SS))
    vehicle_split = split_below.(soc_range, soc_needed)
    value_energy_bysoc = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff=1:FF]
    VV2 = repeat(reshape(value_energy_bysoc, EE, FF, 1), 1, 1, FF)

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = [value_power_action(price, dsoc[pp, ee, ff1, ff2], vehicles_plugged_range[ee]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF]

        soc_needed = soc_scheduled(dt1)
        vehicle_split = split_below.(soc_range, soc_needed)
        valuepns = [value_power_newstate(price, vehicle_split[ff12][1], soc_needed - vehicle_split[ff12][2], vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF]
        valuee = [value_energy.(vehicle_split[ff12][1], vehicle_split[ff12][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF];

        ff12_byaction = discrete_roundbelow.(soc1_byaction, soc_min, soc_max, FF);
        valuepns_byaction = [valuepns[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];
        valuee_byaction = [valuee[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        VV1byactsummc = zeros(Float64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1)
            else
                simustep = get_simustep_stochastic(dt1)
            end
            # Note: We impose costs from soc-below vehicles, but do not adjust state because it pushes up plugged-in soc every period
            statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]), soc1_byaction[pp, ee, ff1, ff2], soc_range[ff2]), vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

            state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2)

            VV1byactthismc = combinebyact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3)
            VV1byactsummc += VV1byactthismc;
        end

        VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction;
        VV1byact[isnan.(VV1byact)] .= -Inf

        bestact = dropdims(argmax(VV1byact, dims=1), dims=1);
        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);

        VV2 = VV1byact[bestact]
        VV2[isnan.(VV2)] .= -Inf
    end

    return strat
end

dt0 = DateTime("2023-07-17T12:00:00")
mcdraws = 1
@time strat = optimize(dt0, SS);

df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5)
plot_standard(df)
plot!(size=(700,400))
savefig("version1-det.pdf")

mcdraws = 20
@time strat = optimize(dt0, SS);

df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5)
plot_standard(df)
plot!(size=(700,400))
savefig("version1-sto.pdf")

## How do the value parameters affect the final level?
mcdraws = 1
results = DataFrame(weight_portion_above=Float64[], weight_portion_below=Float64[], ratio_exponent=Float64[], socend=Float64[])
for wp_above in range(0, .2, 5)
    for wp_below in range(0, .2, 5)
        for re in [.25, .5, 1, 2]
            global weight_portion_above = wp_above
            global weight_portion_below = wp_below
            global ratio_exponent = re
            strat = optimize(dt0, SS);
            df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5)
            push!(results, [weight_portion_above, weight_portion_below, ratio_exponent, df.soc_plugged[end]])
        end
    end
end
