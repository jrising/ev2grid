using Plots
using Dates, DataFrames
using HolidayCalendars, RDates

# Assume this is the day before, and we can observe current storage
# If there are any cars that need recharging to 30% level, handle that outside this system

## Order of events:
## In state1 at start of timestep
## Apply action throughout timestep
## Transition to state2 at end of timestep

# Inputs from current state
vehicles = 4.
vehicle_capacity = 75.7 # kWh
driving_energy_perhour = 329 * 42 / 8 / 1e3 # Wh/mile * miles/day * day/8 hr

# General configuration

timestep = 1. # 1 hour
SS = 36 # project for 1.5 days
# Noon to following midnight

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

energy_max_expected = vehicle_capacity * vehicles * ones(SS) # Battery capacity for 4 Ford Mach-Es
energy_max_expected[1:5] .= 0 # Unavailable from noon to 5pm
energy_max_expected[13 .+ (9:16)] .= 0 # Unavailable 8 hours per day: 9am - 5pm
energy_min_expected = energy_max_expected * .3 # 30% reserved for emergencies
energy_9am = vehicles * 329 * 42 / 1e3 + maximum(energy_min_expected) # Requirement at 9am: Wh/mile * miles/day

energy_min = 0.
energy_max = vehicle_capacity * vehicles

enerfrac_9am = energy_9am / (vehicle_capacity * vehicles)
enerfrac_min = 0.3
enerfrac_max = 1.0

## Additional calculations
function periodstep(steps::Int)
    Hour(round(Int, timestep * steps)) + Minute((timestep * steps - round(Int, timestep * steps)) * 60)
end

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

df = DataFrame(hourstart=DateTime("2024-01-01T00:00:00"):Hour(1):DateTime("2024-12-31T00:00:00"))
df[!, :price] = get_retail_price.(df.hourstart)

pp = plot(df.hourstart, df.price, seriestype=:steppost, label="")
plot!(pp, size=(1000, 400))
savefig("retailprice.png")

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
function value_energy(dt::DateTime, enerfrac::Float64)
    if enerfrac < enerfrac_min
        return -Inf
    elseif enerfrac > 1.
        return 1.
    end

    date_part = Dates.Date(dt)
    dt_9am = DateTime(date_part, Dates.Time(9, 0, 0))

    if dt_9am - periodstep(1) ≤ dt < dt_9am
        if enerfrac < enerfrac_9am
            return -Inf
        else
            return sqrt((enerfrac - enerfrac_9am) / (enerfrac_max - enerfrac_9am))
        end
    else
        return sqrt((enerfrac - enerfrac_min) / (enerfrac_max - enerfrac_min))
    end
end

value_energy(DateTime("2024-07-15T07:34:56"), 0.)
value_energy(DateTime("2024-07-15T07:34:56"), 0.4)
value_energy(DateTime("2024-07-15T08:34:56"), 0.4)
value_energy(DateTime("2024-07-15T09:34:56"), 0.4)
value_energy(DateTime("2024-07-15T08:34:56"), 0.5)

function value_power(dt::DateTime, denerfrac::Float64)
    denergy = denerfrac * vehicle_capacity # charging - discharging in terms of kWh
    if denergy > 0
        -get_retail_price(dt) * denergy / efficiency # cost of energy
    else
        get_retail_price(dt) * denergy # payment for energy
    end
end

"""
Adjust state to reflect the end of the period starting at dt1.
"""
function simustep(dt1::DateTime, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    date_part = Dates.Date(dt1)
    dt_9am = DateTime(date_part, Dates.Time(9, 0, 0))
    dt_5pm = DateTime(date_part, Dates.Time(17, 0, 0))

    if dt_9am - periodstep(1) ≤ dt1 < dt_9am
        ## Move all cars to driving mode
        enerfrac_plugged_2 = enerfrac_plugged_1
        if vehicles_plugged_1 < vehicles
            enerfrac_driving_afterdrive = (enerfrac_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
            enerfrac_driving_2 = (enerfrac_driving_afterdrive * (vehicles - vehicles_plugged_1) + enerfrac_plugged_1 * vehicles_plugged_1) / vehicles
        else
            enerfrac_driving_2 = enerfrac_plugged_1
        end
        vehicles_plugged_2 = 0.
    elseif dt_5pm - periodstep(1) ≤ dt1 < dt_5pm
        ## Move all cars to plugged in
        enerfrac_plugged_2 = (enerfrac_plugged_1 * vehicles_plugged_1 + enerfrac_driving_1 * (vehicles - vehicles_plugged_1)) / vehicles
        enerfrac_driving_2 = enerfrac_driving_1
        vehicles_plugged_2 = vehicles
    else
        enerfrac_plugged_2 = enerfrac_plugged_1
        if vehicles_plugged_1 < vehicles
            enerfrac_driving_2 = (enerfrac_driving_1 * (vehicles - vehicles_plugged_1) * vehicle_capacity - driving_energy_perhour * timestep) / ((vehicles - vehicles_plugged_1) * vehicle_capacity)
            if enerfrac_driving_2 < 0
                enerfrac_driving_2 = 0.
            end
        else
            enerfrac_driving_2 = enerfrac_driving_1
        end
        vehicles_plugged_2 = vehicles_plugged_1
    end

    return vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
end

## Test simustep
simustep(DateTime("2024-07-15T07:34:56"), 4., 1., 0.)
simustep(DateTime("2024-07-15T08:34:56"), 4., 1., 0.)
simustep(DateTime("2024-07-15T09:34:56"), 0., 1., 1.)
simustep(DateTime("2024-07-15T17:34:56"), 0., 1., .9)

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

# function getstatedim(tup::Tuple{Int, Int, Int}, kk::Int)
#     tup[kk]
# end

function makeindex1(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int)
    CartesianIndex(state2ceil, state2base[2], state2base[3])
end

function makeindex2(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int)
    CartesianIndex(state2base[1], state2ceil, state2base[3])
end

function makeindex3(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int)
    CartesianIndex(state2base[1], state2base[2], state2ceil)
end


"""
Optimize the cost using Bellman optimization for a stochastic process.

Arguments:
- SS: Planning horizon (timesteps).

Returns:
- strat: SSxEE matrix representing the optimal strategy.

"""
function optimize(dt0::DateTime, SS::Int)
    strat = zeros(Int64, SS-1, EE, FF, FF);

    enerfrac_range = [0.; range(enerfrac_min, enerfrac_max, FF-1)];
    enerfrac0 = repeat(enerfrac_range', EE, 1, FF);

    # Construct dimensions
    vehicles_plugged_range = range(0., vehicles, EE);

    # Construct exogenous change levels
    denerfrac_FF = [make_actions(enerfrac_plugged) for enerfrac_plugged=enerfrac_range];
    denerfrac = [denerfrac_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, enerfrac_driving=enerfrac_range];

    enerfrac0_byaction = repeat(reshape(enerfrac0, 1, EE, FF, FF), PP, 1, 1, 1);

    # STEP 1: Calculate V[S] under every scenario
    VV2 = value_energy.(dt0 + periodstep(SS), enerfrac0);

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        # if dt0 + periodstep(tt) < DateTime(Dates.Date(dt0 + periodstep(tt)), Dates.Time(9, 0, 0))
        #     break
        # end

        enerfrac1_byaction = enerfrac0_byaction .+ denerfrac;

        valuep = value_power.(dt0 + periodstep(tt), denerfrac);
        valuee = value_energy.(dt0 + periodstep(tt), enerfrac1_byaction);

        statevar2 = [simustep(dt0 + periodstep(tt), vehicles_plugged_range[ee], enerfrac1_byaction[pp, ee, ff1, ff2], enerfrac_range[ff2]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

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

        VV1byact = ((probbase1 + probbase2 + probbase3) .* VV1base + (1 .- probbase1) .* VV1ceil1 + (1 .- probbase2) .* VV1ceil2 + (1 .- probbase3) .* VV1ceil3) / 3 + valuep + valuee * hourly_valuee;
        VV1byact[isnan.(VV1byact)] .= -Inf

        bestact = dropdims(argmax(VV1byact, dims=1), dims=1);
        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);

        VV2 = VV1byact[bestact]
        VV2[isnan.(VV2)] .= -Inf
    end

    return strat
end

dt0 = DateTime("2024-07-15T12:00:00")
strat = optimize(dt0, SS)

"""
Simulate without a strategy for SS time steps.
"""
function simu_inactive(dt1::DateTime, SS::Int, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    rows = Tuple{DateTime, Float64, Float64, Float64}[]
    push!(rows, (dt1 + periodstep(1), vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1))
    for tt in 1:(SS-1)
        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = simustep(dt1 + periodstep(tt), vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)
        push!(rows, (dt1 + periodstep(tt + 1), vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :vehicles_plugged, :enerfrac_plugged, :enerfrac_driving])
end

vehicles_plugged_1 = 4.
enerfrac_plugged_1 = 0.5
enerfrac_driving_1 = 0.5
df = simu_inactive(dt0, SS, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)

pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label="")

"""
Plot a strategy over time.
"""
function simu_strat(dt1::DateTime, strat::AbstractArray{Int}, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    rows = Tuple{DateTime, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}}[]
    push!(rows, (dt1, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1, missing, missing))

    for tt in 1:(size(strat)[1]-1)
        denerfrac = make_actions(enerfrac_plugged_1)

        index = asstate_round((vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1))
        pp = strat[tt, index...]

        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = simustep(dt1 + periodstep(tt), vehicles_plugged_1, enerfrac_plugged_1 + denerfrac[pp], enerfrac_driving_1)
        push!(rows, (dt1 + periodstep(tt), vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2, denerfrac[pp], index))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :vehicles_plugged, :enerfrac_plugged, :enerfrac_driving, :denerfrac, :state])
end

vehicles_plugged_1 = 4.

pp = nothing
for enerfrac_plugged_1 in range(enerfrac_min, enerfrac_max, FF-1)
    enerfrac_driving_1 = enerfrac_plugged_1
    df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)

    if pp == nothing
        pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1)
    else
        plot!(pp, df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1)
    end
end

pp
