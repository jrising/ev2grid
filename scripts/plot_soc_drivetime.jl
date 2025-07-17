include("../src/calc_benefits.jl")

SS = 36
global mcdraws = 1
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

# Plot of time (midnight to midnight) against drive time with colors as SOC 
benefit_matrix = fill(NaN, 24, 24)  # rows: drive duration (0 to 23 hours), columns: drive start hour (0 to 23)
test_start_times = [Time(h, 0, 0) for h in 0:23]
number_vehicles = EE - 1

global soc_matrix_optimized = fill(NaN, 24, 36)
global soc_matrix_rot = fill(NaN, 24, 36)
global soc_matrix_stoch = fill(NaN, 24, 36)
global soc_matrix_stoch_rot = fill(NaN, 24, 36)
global soc_matrix_stoch_baseline = fill(NaN, 24, 36)

## edit this to include optimized, rule of thumb, stochastic optimized, stochastic rule of thumb, and baseline rule of thumb strategies

for drive_start in test_start_times 
    park_time = Time(mod(hour(drive_start) + 8, 24), 0, 0) ## Always drive for 8 hours
    ## optimized strategy 
    global mcdraws = 1
    @time strat, VV = optimize(dt0, SS, drive_start, park_time);

    df_optimized = fullsimulate(dt0, strat, zeros(SS-1), 4., 0.5, 0.5, drive_start, park_time)

    # rule of thumb strategy
    df_rot = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_start), (tt) -> 0., 4., 0.5, 0.5, drive_start, park_time)

    ## stochastic optimized strategy
    global mcdraws = 100
    @time strat, VV = optimize(dt0, SS, drive_start, park_time);
    
    optimized_df = nothing 
    rule_of_thumb_df = nothing
    baseline_rule_of_thumb_df = nothing

    for i in 1:mcdraws
        ## first simulate stochastic events and calculate
        df_stoch = fullsimulate(dt0, strat, zeros(SS-1), 4., 0.5, 0.5, drive_start, park_time)
        # pass events from event_log into rule of thumb simulation
        events = event_log
        df_rot_stoch = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_start), (tt) -> 0., 4., 0.5, 0.5, drive_start, park_time; events=events)
        df_baseline_stoch = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_start), (tt) -> 0., 4., 0.5, 0.5, drive_start, park_time; events=events)
        ## save dataframes
        if i == 1 
            optimized_df = df_stoch
            rule_of_thumb_df = df_rot_stoch
            baseline_rule_of_thumb_df = df_baseline_stoch
        else
            for col in [:soc_plugged, :vehicles_plugged, :soc_driving] 
                optimized_df[!, col] .+= df_stoch[!, col]
                rule_of_thumb_df[!, col] .+= df_rot_stoch[!, col]
                baseline_rule_of_thumb_df[!, col] .+= df_baseline_stoch[!, col]
            end
        end 

    end

    for col in [:soc_plugged, :vehicles_plugged, :soc_driving]
        optimized_df[!, col] ./= mcdraws
        rule_of_thumb_df[!, col] ./= mcdraws
        baseline_rule_of_thumb_df[!, col] ./= mcdraws
    end

    
    global soc_matrix_optimized = fill_soc_matrix(soc_matrix_optimized, df_optimized, drive_start)
    global soc_matrix_rot = fill_soc_matrix(soc_matrix_rot, df_rot, drive_start)
    global soc_matrix_stoch = fill_soc_matrix(soc_matrix_stoch, optimized_df, drive_start)
    global soc_matrix_stoch_rot = fill_soc_matrix(soc_matrix_stoch_rot, rule_of_thumb_df, drive_start)
    global soc_matrix_stoch_baseline = fill_soc_matrix(soc_matrix_stoch_baseline, baseline_rule_of_thumb_df, drive_start)

end 

using Serialization

# # write it out
open("soc_matrix_optimized.bin", "w") do io
    serialize(io, soc_matrix_optimized)
end

open("soc_matrix_rot.bin", "w") do io
    serialize(io, soc_matrix_rot)
end

open("soc_matrix_stoch.bin", "w") do io
    serialize(io, soc_matrix_stoch)
end

open("soc_matrix_stoch_rot.bin", "w") do io
    serialize(io, soc_matrix_stoch_rot)
end

open("soc_matrix_stoch_baseline.bin", "w") do io
    serialize(io, soc_matrix_stoch_baseline)
end

# io = open("soc_matrix_optimized.bin", "r")
# soc_matrix_optimized = deserialize(io)
# close(io)

# io = open("soc_matrix_rot.bin", "r")
# soc_matrix_rot = deserialize(io)
# close(io)

# io = open("soc_matrix_stoch.bin", "r")
# soc_matrix_stoch = deserialize(io)
# close(io)

# io = open("soc_matrix_stoch_rot.bin", "r")
# soc_matrix_stoch_rot = deserialize(io)
# close(io)

# io = open("soc_matrix_stoch_baseline.bin", "r")
# soc_matrix_stoch_baseline = deserialize(io)
# close(io)

plot_soc_heatmap(soc_matrix_optimized, "Optimized Strategy SOC Heatmap")
plot_soc_heatmap(soc_matrix_rot, "Rule of Thumb SOC Heatmap")
plot_soc_heatmap(soc_matrix_stoch, "Stochastic Strategy SOC Heatmap")
plot_soc_heatmap(soc_matrix_stoch_rot, "Stochastic Rule of Thumb SOC Heatmap")
plot_soc_heatmap(soc_matrix_stoch_baseline, "Stochastic Baseline Rule of Thumb SOC Heatmap")

