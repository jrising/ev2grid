include("../src/calc_benefits.jl")

SS = 36
global mcdraws = 1
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

soc_plugged_1 = 0.5
soc_driving_1 = 0.5

benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, drive_starts_time, park_starts_time)

vehicles_plugged_1 = vehicles_plugged_scheduled(dt0, drive_starts_time, park_starts_time)
df1 = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
df1[!, :Approach] .= "Optimized"

df2 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, soc_min, soc_max, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
df2[!, :Approach] .= "Rule of Thumb"

df3 = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
df3[!, :Approach] .= "Baseline"

CSV.write("../results/bytime.csv", [df1; df2; df3])

