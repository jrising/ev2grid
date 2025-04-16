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
    benefits = sum(df[!, "valuep"])

    return benefits
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

function run_rule_of_thumb_stochastic_events_simulation(dt0, drive_starts_time, park_starts_time)
    # global variable to track random events from stochastic simulation
    global event_log = []

    for ev in event_log
        # Only consider events that occur within 1 hour (3600 seconds) of the current simulation time tt
        if ev.event == :delayed_return
            # For a delayed return, postpone park time by an hour
            park_starts_time += Minute(60)
        elseif ev.event == :emergency
            # For emergency events, adjust drive start time based on how many vehicles are needed.
            frac = ev.vehicles_needed / vehicles  # "vehicles" assumed to be the total available
            ## adjust start time for these vehicles
            drive_starts_time = ev.time
            ## returns for these vehicles are increasingly probable with time
            delay = Minute(round(60 * frac))
            drive_starts_time += delay
        end
    end

    
    df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time, true)
    
    benefits = sum(df[!, "valuep"])
    
    return benefits, event_log
end


SS = 36
mcdraws = 100
dt0 = DateTime("2023-07-17T12:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time


benefits_det = run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
benefits, event_log = run_rule_of_thumb_stochastic_events_simulation(dt0, drive_starts_time, park_starts_time)

test_start_times = [Time(h, 0, 0) for h in 0:23]
test_park_times = [Time(h, 0, 0) for h in 0:23]

# benefits_dict = Dict{String, Dict{String, Tuple{Float64, Float64, Float64}}}()

# for drive_starts_time in test_start_times
#     benefits_dict[string(drive_starts_time)] = Dict()

#     for park_starts_time in test_park_times
#         rule_benefits = run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
#         optimized_benefits = run_optimized_simulation(dt0, SS, drive_starts_time, park_starts_time)
#         stochastic_benefits = run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)

#         benefits_dict[string(drive_starts_time)][string(park_starts_time)] = (rule_benefits, optimized_benefits, stochastic_benefits)
#     end
# end

# # Export to LaTeX
# export_to_latex_benefits_table(benefits_dict)


function plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times)
    # Create a 24x24 matrix for benefits
    benefit_matrix = zeros(Float64, 24, 24)  # rows: drive duration (0 to 23 hours), columns: drive start hour (0 to 23)
    
    for drive_start in test_start_times
        for park_start in test_park_times
            # Convert drive_start and park_start to DateTime using dt0's date
            drive_dt = DateTime(year(dt0), month(dt0), day(dt0), hour(drive_start), minute(drive_start), second(drive_start))
            park_dt = DateTime(year(dt0), month(dt0), day(dt0), hour(park_start), minute(park_start), second(park_start))
            
            # If park_dt is earlier than drive_dt, assume it's on the next day
            if park_dt < drive_dt
                park_dt += Dates.Day(1)
            end
            
            # Compute driving duration in hours
            duration = Int((park_dt - drive_dt) / Hour(1))
            
            # Compute the rule-of-thumb benefit for this combination
            benefit = run_rule_of_thumb_simulation(dt0, drive_start, park_start)
            
            # Map the benefit to the matrix: column index = drive start hour, row index = duration
            benefit_matrix[duration+1, hour(drive_start)+1] = benefit
        end
    end
    
    # Create a heatmap with drive start time on the x-axis and driving duration on the y-axis
    heatmap(0:23, 0:23, benefit_matrix,
            xlabel = "Drive Start Time (Hour)",
            ylabel = "Drive Duration (Hours)",
            title = "Rule-of-Thumb Benefits",
            colorbar_title = "Benefit")
end

# Call the plotting function and save the heatmap
plot_rule_of_thumb_benefits(dt0, test_start_times, test_park_times)
savefig("benefits.png")