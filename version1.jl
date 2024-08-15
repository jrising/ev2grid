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

# General configuration

timestep = 1. # 1 hour
SS = 36 # project for 1.5 days
# Noon to following midnight
mcdraws = 1 # 1 for deterministic

# Actions
fracpower_min = -50. / vehicle_capacity # discharge in terms of fraction of energy
fracpower_max = 50. / vehicle_capacity # charging in terms of fraction of energy
PP = 8 # discretized power choices (excluding no-change action)
efficiency = 0.95 # EFF

# States
#   vehicles_plugged: The number of plugged in cars
#   soc_plugged: The fraction of energy available for plugged-in cars
#   soc_driving: The fraction of energy available for driving cars
# Number of states for each dimension:
EE = 5 # 0 - 4 cars
FF = 11 # For both soc_plugged and soc_driving, 0 - 1

include("src/value.jl")
include("src/optutils.jl")

## Checks on configuration parameters

@argcheck mcdraws > 0

"""
Return matrix of changes in energy, across actions.
"""
function make_actions(soc0::Float64)
    fracpower = range(fracpower_min, fracpower_max, length=PP)

    dsocs = timestep * fracpower
    soc1s = max.(0., min.(1., soc0 .+ dsocs))

    [0; soc1s .- soc0]
end

## make_actions(0.5)

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

    # STEP 1: Calculate V[S] under every scenario
    soc_needed = soc_scheduled(dt0 + periodstep(SS))
    vehicle_split = split_below.(soc_range, soc_needed)
    value_energy_bysoc = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], soc_needed) for ff=1:FF]
    VV2 = repeat(reshape(value_energy_bysoc, 1, FF), EE, 1, FF)

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        # if dt0 + periodstep(tt) < DateTime(Dates.Date(dt0 + periodstep(tt)), Dates.Time(9, 0, 0))
        #     break
        # end

        soc1_byaction = soc0_byaction .+ dsoc;

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

        VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction * hourly_valuee;
        VV1byact[isnan.(VV1byact)] .= -Inf

        # if all(.!isfinite.(VV1byact[:, 1, 2, 2]))
        #     println("Cannot support .3")
        #     break
        # end

        bestact = dropdims(argmax(VV1byact, dims=1), dims=1);
        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);

        VV2 = VV1byact[bestact]
        VV2[isnan.(VV2)] .= -Inf
    end

    return strat
end

dt0 = DateTime("2023-07-15T12:00:00")
@time strat = optimize(dt0, SS);

# 0.945021 seconds (14.10 M allocations: 393.177 MiB, 3.75% gc time)

include("src/fullsim.jl")
include("src/plotting.jl")

df = fullsimulate(dt0, strat, zeros(SS-1), 4., 0.5, 0.5)
plot_standard(df)
