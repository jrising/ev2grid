include("../src/calc_benefits.jl")

SS = 36
global mcdraws = 1
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

soc_plugged_1 = 0.5
soc_driving_1 = 0.5

RR = 5 # number of possible regrange values
probfail_penalty = 10.

benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, drive_starts_time, park_starts_time)

vehicles_plugged_1 = vehicles_plugged_scheduled(dt0 + periodstep(1), drive_starts_time, park_starts_time)
df1 = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
df1[!, :Approach] .= "Optimized"

df2 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, soc_min, soc_max, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
df2[!, :Approach] .= "Rule of Thumb"

df3 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
df3[!, :Approach] .= "Baseline"

CSV.write("../results/bytime.csv", [df1; df2; df3])

alldf = []
for drive_starts_hour in 1:23
    println(drive_starts_hour)
    drive_starts_time = Time(drive_starts_hour, 0, 0)
    park_starts_time = Time((drive_starts_hour + 8) % 24, 0, 0)

    benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, drive_starts_time, park_starts_time)

    vehicles_plugged_1 = vehicles_plugged_scheduled(dt0 + periodstep(1), drive_starts_time, park_starts_time)
    df1 = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    df1[!, :Approach] .= "Optimized"
    df1[!, :start_hour] .= drive_starts_hour

    df2 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, soc_min, soc_max, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    df2[!, :Approach] .= "Rule of Thumb"
    df2[!, :start_hour] .= drive_starts_hour

    df3 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    df3[!, :Approach] .= "Baseline"
    df3[!, :start_hour] .= drive_starts_hour

    push!(alldf, df1)
    push!(alldf, df2)
    push!(alldf, df3)
end

CSV.write("../results/bytime-xstart.csv", vcat(alldf...))


## TODO: Create same plots for regrange benefits


vehicles_plugged_1 = vehicles_plugged_scheduled(dt0 + periodstep(1), drive_starts_time, park_starts_time)

probstate = optimize_regrange_probstate_outer_loop(dt0, soc_plugged_1, soc_driving_1, vehicles_plugged_1, drive_starts_time, park_starts_time)
strat, probfail, optregrange = optimize_regrange_probstate(dt0, probstate, drive_starts_time, park_starts_time);
print(optregrange)

df1 = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
df1[!, :Approach] .= "Optimized"
df1[!, :regrange_kw] = [0.; optregrange]  # Add regrange column

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)
df2 = fullsimulate(dt0, dsoc_func, regrange_func, vehicles_plugged_1,  soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
df2[!, :Approach] .= "Rule of Thumb"
df2[!, :regrange_kw] = [regrange_func(tt) for tt in 1:SS]

df3 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
df3[!, :Approach] .= "Baseline"
df3[!, :regrange_kw] = zeros(SS)  # Baseline offers no regulation

CSV.write("results/bytime_regrange.csv", [df1; df2; df3])

alldf = []
for drive_starts_hour in 1:23
    println(drive_starts_hour)
    drive_starts_time = Time(drive_starts_hour, 0, 0)
    park_starts_time = Time((drive_starts_hour + 8) % 24, 0, 0)

    vehicles_plugged_1 = vehicles_plugged_scheduled(dt0 + periodstep(1), drive_starts_time, park_starts_time)

    probstate = optimize_regrange_probstate_outer_loop(dt0, soc_plugged_1, soc_driving_1, vehicles_plugged_1, drive_starts_time, park_starts_time)
    strat, probfail, optregrange = optimize_regrange_probstate(dt0, probstate, drive_starts_time, park_starts_time);

    df1 = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
    df1[!, :Approach] .= "Optimized"
    df1[!, :start_hour] .= drive_starts_hour
    df1[!, :regrange_kw] = [0.; optregrange]  # Add regrange column

    dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)
    df2 = fullsimulate(dt0, dsoc_func, regrange_func, vehicles_plugged_1,  soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
    df2[!, :Approach] .= "Rule of Thumb"
    df2[!, :start_hour] .= drive_starts_hour
    df2[!, :regrange_kw] = [regrange_func(tt) for tt in 1:SS]  # Add regrange column

    df3 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
    df3[!, :Approach] .= "Baseline"
    df3[!, :start_hour] .= drive_starts_hour
    df3[!, :regrange_kw] = zeros(SS)  # Baseline offers no regulation

    push!(alldf, df1)
    push!(alldf, df2)
    push!(alldf, df3)
end

CSV.write("results/bytime-xstart_regrange.csv", vcat(alldf...))

