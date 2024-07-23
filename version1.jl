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

# General configuration

timestep = 1. # 1 hour
SS = 36 # project for 1.5 days
# Noon to following midnight
mcdraws = 1 # 1 for deterministic

hourly_valuee = .1

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

energy_min = 0.
energy_max = vehicle_capacity * vehicles

## Unexpected changes in vehicles
prob_event = 0.01 # per hour, so 1 per 4 days
prob_event_vehicles = [0.5, .25, .125, .125]
prob_event_return = 0.25
prob_delayed_return = 0.1

enerfrac_min = 0.3
enerfrac_max = 1.0

## Checks on configuration parameters

@argcheck mcdraws > 0

## Additional calculations

function is_peak(hourstart::DateTime)
    # Define the start and end times for on-peak hours
    peak_start = Time(12, 0)
    peak_end = Time(20, 0)

    # Check if the date is a weekend or a weekday holiday
    weekday = dayofweek(hourstart)
    if weekday in [6, 7]  # Saturday (6) or Sunday (7)
        return false
    end

    # List of holidays
    cal = calendar(CALENDARS, "NEW YORK")
    if is_holiday(cal, Date(hourstart))
        return false
    end

    # Check the time
    if Time(hourstart) ≥ peak_start && Time(hourstart) < peak_end
        return true
    else
        return false
    end
end

function get_retail_price(hourstart::DateTime)
    # Define the months for summer and winter seasons
    summer_months = 6:9  # June through September
    winter_months = [10:12; 1:5]  # October through May

    if month(hourstart) in summer_months
        return is_peak(hourstart) ? 0.1473 : 0.07242
    else
        return is_peak(hourstart) ? 0.1720 : 0.08356
    end
end

if false
    df = DataFrame(hourstart=DateTime("2024-01-01T00:00:00"):Hour(1):DateTime("2024-12-31T00:00:00"))
    df[!, :price] = get_retail_price.(df.hourstart)

    pp = plot(df.hourstart, df.price, seriestype=:steppost, label="")
    plot!(pp, size=(1000, 400))
    savefig("retailprice.png")
end

"""
Translate a continuous value to a discrete state.

Arguments:
- xx: Continuous value to be discretized.
- xmin: Minimum value of the continuous space.
- xmax: Maximum value of the continuous space.
- num: Number of discrete states.

Returns:
- ii: Discrete state corresponding to the continuous value.
"""
discrete_float(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min((num - 1) * (xx - xmin) / (xmax - xmin) + 1, num))
discrete_round(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min(round(Int, (num - 1) * (xx - xmin) / (xmax - xmin)) + 1, num))
discrete_floatbelow(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min((num - 2) * (xx - xmin) / (xmax - xmin) + 2, num))
discrete_roundbelow(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min(round(Int, (num - 2) * (xx - xmin) / (xmax - xmin)) + 2, num))

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


## Value of having the energy at a given level
function value_energy(portion_below::Float64, enerfrac_avail::Float64, enerfrac_needed::Float64)
    if enerfrac_avail < enerfrac_needed
        return -Inf # Only happens if action made it so
    elseif enerfrac_avail > 1.
        return (1. - portion_below) - weight_portion_below * portion_below
    end

    return (1. - portion_below) * sqrt((enerfrac_avail - enerfrac_needed) / (enerfrac_max - enerfrac_needed)) - weight_portion_below * portion_below
end

value_energy(1., .3, .3)
value_energy(0., 0.5, 0.3)
value_energy(0., 0.5, 0.48)
value_energy(0., 0.5, 0.3)
value_energy(0., 0.6, 0.48)

function value_power_action(price::Float64, denerfrac::Float64)
    denergy = denerfrac * vehicle_capacity * vehicles # charging - discharging in terms of kWh

    if denergy > 0
        -price * denergy / efficiency # cost of energy
    else
        price * denergy # payment for energy
    end
end

function value_power_newstate(price::Float64, portion_below::Float64, enerfrac_toadd_below::Float64)
    denergy_below = enerfrac_toadd_below * vehicle_capacity * vehicles * portion_below
    return -price * denergy_below / efficiency
end

"""
For a given enerfrac_plugged, determine the portion of the plugged-in fleet below the prescribed level
Returns:
 - portion of vehicles below the level
 - enerfrac for low-power portion
 - enerfrac for high-power portion
"""
function split_below(enerfrac_plugged::Float64, enerfrac_needed::Float64)
    if enerfrac_plugged ≤ 0.0
        return 1.0, 0.0, enerfrac_needed
    elseif enerfrac_plugged ≥ 1.0
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


function adjust_below(tup::Tuple{Float64, Float64, Float64}, enerfrac_below::Float64, vehicles_below::Float64)
    vehicles_plugged = tup[1] + vehicles_below
    @assert vehicles_plugged ≤ vehicles + 1e-8

    if vehicles_plugged > 0
        enerfrac_plugged = (enerfrac_below * vehicles_below + tup[2] * tup[1]) / vehicles_plugged
    else
        enerfrac_plugged = tup[2]
    end
    (vehicles_plugged, enerfrac_plugged, tup[3])
end

function asstate_float(tup::Tuple{Float64, Float64, Float64})
    (discrete_float(tup[1], 0., vehicles, EE), discrete_floatbelow(tup[2], enerfrac_min, enerfrac_max, FF), discrete_floatbelow(tup[3], enerfrac_min, enerfrac_max, FF))
end

function asstate_round(tup::Tuple{Float64, Float64, Float64})
    (discrete_round(tup[1], 0., vehicles, EE), discrete_roundbelow(tup[2], enerfrac_min, enerfrac_max, FF), discrete_roundbelow(tup[3], enerfrac_min, enerfrac_max, FF))
end

basestate(tup::Tuple{Float64, Float64, Float64}) = floor.(Int, tup)
ceilstate(tup::Tuple{Float64, Float64, Float64}, kk::Int) = ceil(Int, tup[kk])

function getstatedim(tup::Tuple{Float64, Float64, Float64}, kk::Int)
    tup[kk]
end


makeindex1(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2ceil, state2base[2], state2base[3])
makeindex2(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2base[1], state2ceil, state2base[3])
makeindex3(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2base[1], state2base[2], state2ceil)

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
    vehicles_plugged_range = range(0., vehicles, EE);

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

        VV1byactsummc = zeros(Int64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1)
            else
                simustep = get_simustep_stochastic(dt1)
            end
            # Note: We impose costs from enerfrac-below vehicles, but do not adjust state because it pushes up plugged-in enerfrac every period
            statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]), enerfrac1_byaction[pp, ee, ff1, ff2], enerfrac_range[ff2]), vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

            state2 = asstate_float.(statevar2);
            state2base = basestate.(state2);
            state2ceil1 = ceilstate.(state2, 1);
            state2ceil2 = ceilstate.(state2, 2);
            state2ceil3 = ceilstate.(state2, 3);

            probbase1 = (state2ceil1 - getstatedim.(state2, 1));
            probbase2 = (state2ceil2 - getstatedim.(state2, 2));
            probbase3 = (state2ceil3 - getstatedim.(state2, 3));

            ## Now index into VV2 three times and combine
            VV1base = VV2[CartesianIndex.(state2base)];
            VV1ceil1 = VV2[makeindex1.(state2base, state2ceil1)]; # makeindex1 much faster than anonymous functions
            VV1ceil2 = VV2[makeindex2.(state2base, state2ceil2)];
            VV1ceil3 = VV2[makeindex3.(state2base, state2ceil3)];

            VV1byactthismc = ((probbase1 + probbase2 + probbase3) .* VV1base + (1 .- probbase1) .* VV1ceil1 + (1 .- probbase2) .* VV1ceil2 + (1 .- probbase3) .* VV1ceil3) / 3;
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

"""
Simulate without a strategy for SS time steps.
"""
function simu_inactive(dt0::DateTime, SS::Int, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    # datetime, enerfrac_needed, vehicles_plugged_1, portion_below, enerfrac_toadd_below, enerfrac_avail_1, enerfrac_driving_1
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64}[]

    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

    push!(rows, (dt0 + periodstep(1), enerfrac_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], enerfrac_plugged_1, enerfrac_driving_1))

    for tt in 1:(SS-1)
        println(tt)
        enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

        simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], enerfrac_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt + 1), enerfrac_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], enerfrac_plugged_2, enerfrac_driving_2))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :enerfrac_needed, :vehicles_plugged, :portion_below, :enerfrac_below, :enerfrac_plugged, :enerfrac_driving])
end

vehicles_plugged_1 = 4.
enerfrac_plugged_1 = 0.5
enerfrac_driving_1 = 0.5
df = simu_inactive(dt0, SS, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)

pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label="")

"""
Plot a strategy over time.
"""
function simu_strat(dt0::DateTime, strat::AbstractArray{Int}, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}}[]

    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

    push!(rows, (dt0, enerfrac_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], enerfrac_plugged_1, enerfrac_driving_1, missing, missing))

    for tt in 1:(size(strat)[1]-1)
        denerfrac = make_actions(enerfrac_plugged_1)

        index = asstate_round((vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1))
        pp = strat[tt, index...]

        enerfrac_plugged_2 = enerfrac_plugged_1 + denerfrac[pp]

        enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(enerfrac_plugged_2, enerfrac_needed)

        if mcdraws == 1
            simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        else
            simustep = get_simustep_stochastic(dt0 + periodstep(tt))
        end
        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], enerfrac_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt), enerfrac_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], enerfrac_plugged_2, enerfrac_driving_2, denerfrac[pp], index))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :enerfrac_needed, :vehicles_plugged, :portion_below, :enerfrac_toadd_below, :enerfrac_plugged, :enerfrac_driving, :denerfrac, :state])
end

vehicles_plugged_1 = 4.

pp = nothing
for enerfrac_plugged_1 in range(enerfrac_min, enerfrac_max, FF-1)
    enerfrac_driving_1 = enerfrac_plugged_1
    df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)

    if pp == nothing
        pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1, legend=false)
    else
        plot!(pp, df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1, legend=false)
    end
end

pp
