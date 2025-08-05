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
    vehicles_plugged_1 = 4.
    df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, 0.3, 0.95), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function rectangle(w, h, x, y)
        return Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
end

function run_optimized_simulation(dt0, SS, drive_starts_time, park_starts_time)
    vehicles_plugged_1 = 4.
    global mcdraws = 1
    @time strat, VV = optimize(dt0, SS, drive_starts_time, park_starts_time);

    df = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function run_optimized_regrange_simulation(dt0, SS, drive_starts_time, park_starts_time)
    vehicles_plugged_1 = 4.

    global probstate = zeros(SS-1, EE, FF, FF);
    df = fullsimulate(dt0, zeros(SS-1), vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
    for ii in 1:(nrow(df) - 1)
        statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((df.vehicles_plugged[ii], df.soc_plugged[ii], df.soc_driving[ii]))
        probstate[ii, statebase...] = (probbase1 + probbase2 + probbase3) / 3
        probstate[ii, makeindex1(statebase, stateceil1)] += (1 - probbase1) / 3
        probstate[ii, makeindex2(statebase, stateceil2)] += (1 - probbase2) / 3
        probstate[ii, makeindex3(statebase, stateceil3)] += (1 - probbase3) / 3
    end

    global probstate = probstate / 2 .+ repeat(ones(1, EE, FF, FF) / (EE * FF * FF), SS-1) / 2;

    for ll in 1:10
        println("Loop $(ll)")
        strat, probfail, optregrange = optimize_regrange(dt0, probstate, drive_starts_time, park_starts_time);

        dfall = nothing
        for ii in 1:mcdraws
            local df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
            energy = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) .* df.soc_above
            energy_minallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.3
            energy_maxallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.95
            df[!, :regrange_avail] = min.(energy_maxallow - energy, energy - energy_minallow)

            if dfall == nothing
                dfall = df
            else
                append!(dfall, df)
            end
        end

        global probstate = zeros(SS-1, EE, FF, FF);
        for tt in 1:SS-1
            dt1 = dt0 + periodstep(tt - 1)

            for state in unique(dfall.state[dfall.datetime .== dt1])
                probstate[tt, state...] = sum(map(rowstate -> rowstate == state, dfall.state[dfall.datetime .== dt1])) / sum(dfall.datetime .== dt1)
            end
        end
    end

    strat, probfail, optregrange = optimize_regrange(dt0, probstate, drive_starts_time, park_starts_time);

    df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
    df[!, :optregrange] = [0.; optregrange]
    benefits = sum(df[!, "valuep"]) + sum(df[!, "valuer"])
    return benefits
end 

function run_rule_of_thumb_regrange_simulation(dt0, SS, drive_starts_time, park_starts_time)
    vehicles_plugged_1 = 4.
    global mcdraws = 1
    df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state), (tt, state) -> thumbrule_regrange(tt, state), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    benefits = sum(df[!, "valuep"])
    return benefits
end 

function run_optimized_stochastic_simulation(dt0, SS, drive_starts_time, park_starts_time)
    vehicles_plugged_1 = 4.
    global mcdraws = 100
    @time strat, VV = optimize(dt0, SS, drive_starts_time, park_starts_time);

    df = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time, true)
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

    vehicles_plugged_1 = 4.

    for mc in mcdraws
        ## first simulate stochastic events and calculate
        df = fullsimulate(dt0, strat, zeros(SS-1), vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
        benefits_stoch = sum(df[!, "valuep"])
        push!(benefits_list_optimized, benefits_stoch)
        # pass events from event_log into rule of thumb simulation
        events = event_log
        df = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, 0.3, 0.95), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time; events=events)
        benefits_rot = sum(df[!, "valuep"])
        push!(benefits_list_rot, benefits_rot)
        df = fullsimulate_with_events(dt0, (tt, state) -> get_dsoc_thumbrule_baseline(tt, state), (tt) -> 0., vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time; events=events)
        benefits_rot_baseline = sum(df[!, "valuep"])
        push!(benefits_list_rot_baseline, benefits_rot_baseline)

    end
    mean_benefits_optimized = mean(benefits_list_optimized)
    mean_benefits_rot = mean(benefits_list_rot)
    mean_benefits_rot_baseline = mean(benefits_list_rot_baseline)

    return mean_benefits_optimized, mean_benefits_rot, mean_benefits_rot_baseline
end


function fill_soc_matrix(soc_matrix, df, drive_start)
    df[!, "number_vehicles"] .= number_vehicles
    soc = ((df[!, "number_vehicles"] .- df[!, "vehicles_plugged"]).*df[!, "soc_driving"] .+ df[!,"vehicles_plugged"].*df[!,"soc_plugged"]) ./ number_vehicles
    time = df[!, "datetime"]
    # print(df)
    for (i, t) in enumerate(time)
        soc_matrix[hour(drive_start)+1, i] = soc[i]
    end
    # print(soc_matrix)
    return soc_matrix
end

function plot_soc_heatmap(soc_matrix, label)

    heatmap(0:23, 0:23, soc_matrix[:, 1:end-12],
        xlabel = "Time (Midnight to Midnight)",
        ylabel = "Drive Start Time (Hour)",
        title = label,
        colorbar_title = "SOC",
        aspect_ratio = 1,
        c=:viridis, 
        xlims=(0,23), 
        ylims=(0,23))

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

    savefig(label * "_soc_heatmap.png")


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


function plot_rule_of_thumb_benefits_difference(dt0, test_start_times, test_park_times, benefits_dict, title, index)
    # Create a 24x24 matrix for benefits, initialized with NaN
    benefit_matrix = fill(NaN, 24, 24)  # rows: drive duration (0 to 23 hours), columns: drive start hour (0 to 23)

    for drive_start in test_start_times
        for park_start in test_park_times
            if park_starts_time <= drive_starts_time
                continue
            end

            rule_benefits_tuple = benefits_dict[string(drive_start)][string(park_start)]
            benefit = rule_benefits_tuple[index] - rule_benefits_tuple[5] ## difference between rule of thumb and baseline rule of thumb
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
            colorbar_title = "Difference in benefit")
    vline!([12, 20], linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
    hline!([12, 20], linestyle = :dash, color = :red, linewidth = 1.5, label = nothing)
    # Annotate "peak pricing" at the corresponding locations
    # annotate!(12, 23.5, text("Peak Pricing", :black, :center, 10, rotation=90))
    # annotate!(0.5, 20, text("Peak Pricing", :black, :left, 10))
end
