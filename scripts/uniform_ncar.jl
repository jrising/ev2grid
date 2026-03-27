include("../src/calc_benefits.jl")

SS = 36
global mcdraws = 1
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Drive start time
park_starts_time = Time(17, 0, 0)  # Park start time
RR = 5 # number of possible regrange values
probfail_penalty = 10.

alldf = []
soc_range = [0.; range(soc_min, soc_max, FF-1)]
for soc_plugged_1 in soc_range
    println("NEW SOC LEVEL")
    println(soc_plugged_1)
    soc_driving_1 = soc_plugged_1

    vehicles_plugged_1 = vehicles_plugged_scheduled(dt0 + periodstep(1), drive_starts_time, park_starts_time)

    probstate = optimize_regrange_probstate_outer_loop(dt0, soc_plugged_1, soc_driving_1, vehicles_plugged_1, drive_starts_time, park_starts_time)
    strat, probfail, optregrange = optimize_regrange_probstate(dt0, probstate, drive_starts_time, park_starts_time);

    df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
    df[!, :optregrange] = [0.; optregrange]
    df[!, :soc_1] .= soc_plugged_1

    push!(alldf, df)
end

CSV.write("../results/bytime-xsoc.csv", vcat(alldf...))
