using DataFrames
using Printf
using Plots

include("bizutils.jl")
include("customer.jl")
include("simulate.jl")
include("retail.jl")
include("config.jl")
include("value.jl")
include("optutils.jl")
include("fullsim.jl")
include("plotting.jl")
include("thumbrule.jl")
include("optfuncs.jl")


function run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
    df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function rectangle(w, h, x, y)
        return Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
end

function run_optimized_simulation(dt0, SS, drive_starts_time, park_starts_time)
    mcdraws = 1
    @time strat, VV = optimize(dt0, SS, drive_starts_time, park_starts_time);

    df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5, drive_starts_time, park_starts_time)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)
    mcdraws = 100
    @time strat, VV = optimize(dt0, SS, drive_starts_time, park_starts_time);

    df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5, drive_starts_time, park_starts_time, true)
    benefits = sum(df[!, "valuep"]) ## conduct a single run for a simulation with stochastic events

    # open("event_log_debug.txt", "w") do io ## to help with debugging print out event_log
    #     for ev in event_log
    #         println(io, ev)
    #     end
    # end

    return benefits, strat
end


function export_to_latex_benefits_table(benefits_dict)
    filename = "results/benefits_table.tex"
    open(filename, "w") do io
        write(io, "\\documentclass{article}\n")
        write(io, "\\usepackage{multirow}\n")
        write(io, "\\usepackage{booktabs}\n")
        write(io, "\\begin{document}\n")
        write(io, "\\begin{table}[h!]\n")
        write(io, "\\centering\n")

        num_drive_starts = length(keys(benefits_dict)) ## Number of different start times

        write(io, "\\begin{tabular}{|c|" * repeat("c|", num_drive_starts + 1) * "}\n")
        write(io, "\\hline\n")
        write(io, "\\multicolumn{1}{|c|}{} & \\multicolumn{$(num_drive_starts)}{|c|}{\\textbf{Drive Start Time}} \\\\ \\hline\n")
        write(io, "\\textbf{Park Start Time} ")
        for drive_start in keys(benefits_dict)
            write(io, @sprintf("& \\textbf{%s} ", drive_start))
        end
        write(io, " \\\\ \\hline\n")

        first_drive_start = first(keys(benefits_dict))
        park_times = sort(collect(keys(benefits_dict[first_drive_start])))

        for park_start in park_times
            park_start_str = string(park_start)
            write(io, @sprintf("%s", park_start_str))

            for drive_start in keys(benefits_dict)
                if park_starts_time <= drive_starts_time
                    continue
                end

                if haskey(benefits_dict[drive_start], park_start_str)
                    benefits = benefits_dict[drive_start][park_start_str]
                    write(io, @sprintf(" & (%.2f, %.2f, %.2f)", benefits[1], benefits[2], benefits[3]))
                else
                    write(io, " & -")  # Placeholder for missing values
                end
            end
            write(io, " \\\\ \n")
        end

        write(io, "\\hline\n")
        write(io, "\\end{tabular}\n")
        write(io, "\\caption{Each cell contains three values: (1) Rule-of-thumb benefits, (2) Optimized benefits, and (3) Stochastic optimized benefits, calculated for different combinations of drive start and park start times.}\n")
        write(io, "\\end{table}\n")
        write(io, "\\end{document}\n")
    end
end

function run_rule_of_thumb_stochastic_events_simulation(dt0, strat, mcdraws, drive_starts_time, park_starts_time)
    # global event_log = [] ## to help with debugging read in previously printed out event_log
    # open("event_log_debug.txt", "r") do io
    #     for line in eachline(io)
    #         push!(event_log, eval(Meta.parse(line)))
    #     end
    # end
    benefits_list_optimized = []
    benefits_list_rot = []
    benefits_list_rot_baseline = []

    for mc in mcdraws
        ## first simulate stochastic events and calculate
        df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5, drive_starts_time, park_starts_time)
        benefits_stoch = sum(df[!, "valuep"])
        push!(benefits_list_optimized, benefits_stoch)
        # pass events from event_log into rule of thumb simulation
        events = event_log
        df = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time; events=events)
        benefits_rot = sum(df[!, "valuep"])
        push!(benefits_list_rot, benefits_rot)
        df = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time; events=events)
        benefits_rot_baseline = sum(df[!, "valuep"])
        push!(benefits_list_rot_baseline, benefits_rot_baseline)

    end
    mean_benefits_optimized = mean(benefits_list_optimized)
    mean_benefits_rot = mean(benefits_list_rot)
    mean_benefits_rot_baseline = mean(benefits_list_rot_baseline)

    return mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline
end


SS = 36
mcdraws = 100
dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time


# First, run the stochastic simulation to develop strategy for charging with stochastic events
# benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)

# # # Then, run a bunch of simulations to calculate mean benefits under optimal strategy and rule of thumb with events
# mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline = run_rule_of_thumb_stochastic_events_simulation(dt0, strat, mcdraws, drive_starts_time, park_starts_time)

# println("Stochastic Optimized Benefits mean: ", mean_benefits_optimized)
# println("Stochastic Rule of Thumb Benefits mean: ", mean_benefits_rot)
# println("Baseline Rule of Thumb Benefits mean: ", mean_benefits_rot_baseline)


## create heat map of value for parking and starting driving at different times
test_start_times = [Time(h, 0, 0) for h in 0:23]
test_park_times = [Time(h, 0, 0) for h in 0:23]

# benefits_dict = Dict{String, Dict{String, Tuple{Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}, Union{Float64, Missing}}}}()

# for drive_starts_time in test_start_times
#     benefits_dict[string(drive_starts_time)] = Dict()
#     for park_starts_time in test_park_times
#         if park_starts_time <= drive_starts_time
#             benefits_dict[string(drive_starts_time)][string(park_starts_time)] = (missing, missing, missing, missing, missing)
#             continue 
#         end
#         rule_benefits = run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
#         optimized_benefits = run_optimized_simulation(dt0, SS, drive_starts_time, park_starts_time)
#         benefits_stoch, strat = run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)
#         mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline = run_rule_of_thumb_stochastic_events_simulation(dt0, strat, mcdraws, drive_starts_time, park_starts_time)


#         benefits_dict[string(drive_starts_time)][string(park_starts_time)] = (rule_benefits, optimized_benefits, mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline)
#     end
# end
# using Serialization

# # write it out
# open("benefits_dict.bin", "w") do io
#     serialize(io, benefits_dict)
# end

# # later, to read it back:
using Serialization

io = open("benefits_dict.bin", "r")
benefits_dict_restored = deserialize(io)
close(io)
# Export to LaTeX
# export_to_latex_benefits_table(benefits_dict)


function plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times, benefits_dict, title, index)
    # Create a 24x24 matrix for benefits, initialized with NaN
    benefit_matrix = fill(NaN, 24, 24)  # rows: drive duration (0 to 23 hours), columns: drive start hour (0 to 23)

    for drive_start in test_start_times
        for park_start in test_park_times
            if park_starts_time <= drive_starts_time
                continue
            end

            rule_benefits_tuple = benefits_dict[string(drive_start)][string(park_start)]
            benefit = rule_benefits_tuple[index]
            if !ismissing(benefit)
                benefit_matrix[hour(park_start)+1, hour(drive_start)+1,] = benefit
            end
        end
    end

    # Create a heatmap with drive start time on the x-axis and driving duration on the y-axis
    heatmap(0:23, 0:23, benefit_matrix,
            xlabel = "Drive Start Time (Hour)",
            ylabel = "Drive Park Time (Hour)",
            title = title,
            colorbar_title = "Benefit")
    vline!([12, 20], linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
    hline!([12, 20], linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
    # Annotate "peak pricing" at the corresponding locations
    # annotate!(12, 23.5, text("Peak Pricing", :black, :center, 10, rotation=90))
    # annotate!(0.5, 20, text("Peak Pricing", :black, :left, 10))
end


## July 7 plot of time (midnight to midnight) against drive time with colors as SOC 
# benefit_matrix = fill(NaN, 24, 24)  # rows: drive duration (0 to 23 hours), columns: drive start hour (0 to 23)
test_start_times = [Time(h, 0, 0) for h in 0:23]
number_vehicles = EE - 1

# soc_matrix = fill(NaN, 24, 36)

# for drive_start in test_start_times 
#     park_time = Time(mod(hour(drive_start) + 8, 24), 0, 0) ## Always drive for 8 hours
#     @time strat, VV = optimize(dt0, SS, drive_start, park_time);

#     df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5, drive_start, park_time)
#     df[!, "number_vehicles"] .= number_vehicles
#     soc = ((df[!, "number_vehicles"] .- df[!, "vehicles_plugged"]).*df[!, "soc_driving"] .+ df[!,"vehicles_plugged"].*df[!,"soc_plugged"]) ./ number_vehicles
#     time = df[!, "datetime"]
#     for (i, t) in enumerate(time)
#         soc_matrix[hour(drive_start)+1, i] = soc[i]
#     end
# end 

using Serialization

# write it out
# open("soc_matrix.bin", "w") do io
#     serialize(io, soc_matrix)
# end

io = open("soc_matrix.bin", "r")
soc_matrix = deserialize(io)
close(io)

# heatmap(0:35, 0:23, soc_matrix,

heatmap(0:23, 0:23, soc_matrix[:, 1:end-12],
        xlabel = "Time (Midnight to Midnight)",
        ylabel = "Drive Start Time (Hour)",
        title = "Optimized Strategy",
        colorbar_title = "SOC",
        aspect_ratio = 1,
        c=:viridis, 
        xlims=(0,23), 
        ylims=(0,23))
# vline!([8], ymin=0, ymax=1, linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
# hline!([1], xmin=0, xmax=8, linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
for p in 0:23
    if p < 16
        plot!(rectangle(8, 1, p, p), fill=false, label=nothing, linewidth=2, color=:black, fillcolor=nothing)
    else 
        plot!(rectangle(24-p, 1, p, p), fill=false, label=nothing, linewidth=2, color=:black, fillcolor=nothing)
        plot!(rectangle(p - 16, 1, 0, p), fill=false, label=nothing, linewidth=2, color=:black, fillcolor=nothing)
    end
end
vline!([12, 20], linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
hline!([12, 20], linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)

savefig("soc_heatmap.png")


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