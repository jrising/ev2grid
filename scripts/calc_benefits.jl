include("../src/calc_benefits.jl")

SS = 36
mcdraws = 100
dt0 = DateTime("2023-07-17T12:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time


# First, run the stochastic simulation to develop strategy for charging with stochastic events
benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)

# # Then, run a bunch of simulations to calculate mean benefits under optimal strategy and rule of thumb with events
mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline = run_rule_of_thumb_stochastic_events_simulation(dt0, strat, mcdraws, drive_starts_time, park_starts_time)

println("Stochastic Optimized Benefits mean: ", mean_benefits_optimized)
println("Stochastic Rule of Thumb Benefits mean: ", mean_benefits_rot)
println("Baseline Rule of Thumb Benefits mean: ", mean_benefits_rot_baseline)


## create heat map of value for parking and starting driving at different times
test_start_times = [Time(h, 0, 0) for h in 0:23]
test_park_times = [Time(h, 0, 0) for h in 0:23]

benefits_dict = Dict{String, Dict{String, Tuple{Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}}}}()

for drive_starts_time in test_start_times
    benefits_dict[string(drive_starts_time)] = Dict()
    for park_starts_time in test_park_times
        if park_starts_time <= drive_starts_time
            benefits_dict[string(drive_starts_time)][string(park_starts_time)] = (missing, missing, missing, missing, missing)
            continue
        end
        rule_benefits = run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
        optimized_benefits = run_optimized_simulation(dt0, SS, drive_starts_time, park_starts_time)
        benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)
        mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline = run_rule_of_thumb_stochastic_events_simulation(dt0, strat, mcdraws, drive_starts_time, park_starts_time)


        benefits_dict[string(drive_starts_time)][string(park_starts_time)] = (rule_benefits, optimized_benefits, mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline)
    end
end
using Serialization

# write it out
open("benefits_dict.bin", "w") do io
    serialize(io, benefits_dict)
end

# # later, to read it back:
# using Serialization

# io = open("benefits_dict.bin", "r")
# benefits_dict_restored = deserialize(io)
# close(io)
# Export to LaTeX
# export_to_latex_benefits_table(benefits_dict)

# Call the plotting function using precomputed benefits and save the heatmap
plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times, benefits_dict_restored, "Rule of Thumb", 1)
savefig("benefits_rot.png")

plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times, benefits_dict_restored, "Optimized", 2)
savefig("benefits_optimized.png")

plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times, benefits_dict_restored, "Stochastic Optimized", 3)
savefig("benefits_optimized_stoch.png")

plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times, benefits_dict_restored, "Stochastic Rule of Thumb", 4)
savefig("benefits_optimized_stoch_rot.png")

plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times, benefits_dict_restored, "Stochastic Baseline Rule of Thumb", 5)
savefig("benefits_optimized_stoch_rot_baseline.png")
