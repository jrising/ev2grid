using DataFrames
using Printf

include("src/bizutils.jl")
include("src/customer.jl")
include("src/simulate.jl")
include("src/retail.jl")
include("src/config.jl")
include("src/value.jl")
include("src/optutils.jl")
include("src/fullsim.jl")
include("src/plotting.jl")
include("src/thumbrule.jl")
include("src/optfuncs.jl")


function run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
    df = fullsimulate(dt0, (tt, state) -> get_dsoc(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time)
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



SS = 36
mcdraws = 100
dt0 = DateTime("2023-07-17T12:00:00")

test_start_times = [Time(7,0,0), Time(8,0,0), Time(9,0,0)]
test_park_times = [Time(16,0,0), Time(17,0,0), Time(18,0,0)]

benefits_dict = Dict{String, Dict{String, Tuple{Float64, Float64, Float64}}}()

for drive_starts_time in test_start_times
    benefits_dict[string(drive_starts_time)] = Dict()

    for park_starts_time in test_park_times
        rule_benefits = run_rule_of_thumb_simulation(dt0, drive_starts_time, park_starts_time)
        optimized_benefits = run_optimized_simulation(dt0, SS, drive_starts_time, park_starts_time)
        stochastic_benefits = run_optimized_stochastic_simulation(dt0, SS, mcdraws, drive_starts_time, park_starts_time)

        benefits_dict[string(drive_starts_time)][string(park_starts_time)] = (rule_benefits, optimized_benefits, stochastic_benefits)
    end
end

# Export to LaTeX
export_to_latex_benefits_table(benefits_dict)