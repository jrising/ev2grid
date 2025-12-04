include("../src/calc_benefits.jl")

SS = 36
global mcdraws = 1
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

soc_plugged_1 = 0.5
soc_driving_1 = 0.5

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


## TODO: Continue from here

benefits_list_optimized = []
    benefits_list_rot = []
    benefits_list_rot_baseline = []


    for mc in mcdraws
        ## first simulate stochastic events and calculate
        df = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
        benefits_stoch = sum(df[!, "valuep"])
        push!(benefits_list_optimized, benefits_stoch)
        # pass events from event_log into rule of thumb simulation
        events = event_log
        df = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, soc_min, soc_max, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time; events=events)
        benefits_rot = sum(df[!, "valuep"])
        push!(benefits_list_rot, benefits_rot)
        df = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time; events=events)
        benefits_rot_baseline = sum(df[!, "valuep"])
        push!(benefits_list_rot_baseline, benefits_rot_baseline)

    end
    mean_benefits_optimized = mean(benefits_list_optimized)
    mean_benefits_rot = mean(benefits_list_rot)
    mean_benefits_rot_baseline = mean(benefits_list_rot_baseline)





# # Then, run a bunch of simulations to calculate mean benefits under optimal strategy and rule of thumb with events
mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline = run_rule_of_thumb_stochastic_events_simulation(dt0, strat, mcdraws, drive_starts_time, park_starts_time)

println("Stochastic Optimized Benefits mean: ", mean_benefits_optimized)
println("Stochastic Rule of Thumb Benefits mean: ", mean_benefits_rot)
println("Baseline Rule of Thumb Benefits mean: ", mean_benefits_rot_baseline)
