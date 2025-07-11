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

