using Dates, DataFrames
using StatsBase, Random, Distributions
using Roots
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

include("utils.jl")
include("customer.jl")
include("simulate.jl")
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
#   enerfrac_plugged: The fraction of energy available for plugged-in cars
#   enerfrac_driving: The fraction of energy available for driving cars
# Number of states for each dimension:
EE = 5 # 0 - 4 cars
FF = 11 # For both enerfrac_plugged and enerfrac_driving, 0 - 1

## Unexpected changes in vehicles
prob_event = 0.01 # per hour, so 1 per 4 days
prob_event_vehicles = [0.5, .25, .125, .125]
prob_event_return = 0.25
prob_delayed_return = 0.1

enerfrac_min = 0.3
enerfrac_max = 1.0

include("value.jl")
include("optutils.jl")

## Checks on configuration parameters

@argcheck mcdraws > 0

"""
Return matrix of changes in energy, across actions.
"""
function make_actions(enerfrac0::Float64)
    fracpower = range(fracpower_min, fracpower_max, length=PP)

    denerfracs = timestep * fracpower
    enerfrac1s = max.(0., min.(1., enerfrac0 .+ denerfracs))

    [0; enerfrac1s .- enerfrac0]
end

## make_actions(0.5)

"""
For a given enerfrac_plugged, determine the portion of the plugged-in fleet below the prescribed level
Returns:
 - portion of vehicles below the level
 - enerfrac for low-power portion
 - enerfrac for high-power portion
"""
function split_below(enerfrac_plugged::Float64, enerfrac_needed::Float64)
    # Limits come from truncated(Normal(mu, 0.05), lower=0., upper=1.)
    if enerfrac_plugged ≤ 0.03989422804014327
        return 1.0, 0.0, enerfrac_needed
    elseif enerfrac_plugged ≥ 0.9601057719598567
        return 0.0, enerfrac_needed, 1.0
    end
    ## Use normal distribution, so I can calculate means of truncated
    # Need to get a truncated normal with the desired mean
    f(mu) = mean(truncated(Normal(mu, 0.05), lower=0., upper=1.)) - enerfrac_plugged
    mu = find_zero(f, (0., 1.))
    dist = truncated(Normal(mu, 0.05), lower=0., upper=1.)
    portion_below = cdf(dist, enerfrac_needed)
    portion_below, mean(truncated(dist, upper=enerfrac_needed)), mean(truncated(dist, lower=enerfrac_needed))
end

split_below(0.5, 0.5)
split_below(0.5, 0.3)
split_below(0.7, 0.3)
split_below(0.3, 0.3)
split_below(0.1, 0.3)

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

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        # if dt0 + periodstep(tt) < DateTime(Dates.Date(dt0 + periodstep(tt)), Dates.Time(9, 0, 0))
        #     break
        # end

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

        VV1byactsummc = zeros(Float64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1)
            else
                simustep = get_simustep_stochastic(dt1)
            end
            state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 =
                simustate2(simustep, vehicles_plugged_range, vehicle_split, ff12_byaction, enerfrac1_byaction, enerfrac_range)
            VV1byactthismc = calcVV1byact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3)
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

dt0 = DateTime("2024-07-15T12:00:00")
@time strat = optimize(dt0, SS);

# 0.945021 seconds (14.10 M allocations: 393.177 MiB, 3.75% gc time)

include("fullsim.jl")
