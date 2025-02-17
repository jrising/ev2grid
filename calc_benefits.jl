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

function get_dsoc(tt, state, drive_starts_time)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)

    dt_drive = DateTime(Dates.Date(dt1)+ Dates.Day(1), drive_starts_time)
    time_available = (dt_drive - dt1).value / 3600 / 1000 # in hours

    # Calculate floor SOC needed to ramp up to 80% by 9am

    charge_needed = max(0,0.80 - soc_plugged_1)
    time_required = charge_needed / (fracpower_max * timestep)
    SOC_floor = max(0.3, soc_plugged_1 - time_required * fracpower_max) ## Charge can't go below 0.3
    
    # calculate if a price change occurs before we need to start driving 
    # take dt1 and drive_starts_time and test each hour in between for if it switches peak to non peak 
    # alternative here is to hard code that the pricing switches at 8 pm and 12 pm 

    price_switch = false  # Default: No switch found
    current_peak_status = is_peak(dt1)

    for t in 1:time_available
        dt_check = dt1 + Hour(t)
        if is_peak(dt_check) != current_peak_status  # Detects a switch
            price_switch = true 
            break  # Exit loop when first switch is found
        end
    end
         
    ## don't I need some logic to say that if they have to charge depending on the amount of driving they do? 


    if is_peak(dt1) && SOC_floor <= soc_plugged_1 && price_switch
        soc_goal = SOC_floor # Discharge if peak price, above floor, and price change coming
    elseif !is_peak(dt1)
        soc_goal = 0.95 # Charge if off peak
    elseif SOC_floor > soc_plugged_1
        soc_goal = 0.80 # Charge if below floor
    end


    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end

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
    
    df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5, drive_starts_time, park_starts_time)
    benefits = sum(df[!, "valuep"])

    return benefits
end


function export_to_latex_benefits_table(benefits_dict)
    filename = "benefits_table.tex"
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

function optimize(dt0::DateTime, SS::Int, drive_starts_time::Time, park_starts_time::Time)
    strat = zeros(Int64, SS-1, EE, FF, FF);

    # Construct dimensions
    soc_range = [0.; range(soc_min, soc_max, FF-1)];
    vehicles_plugged_range = collect(range(0., vehicles, EE));

    # Construct exogenous change levels
    dsoc_FF = [make_actions(soc_plugged, soc_range) for soc_plugged=soc_range];
    dsoc = [dsoc_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, soc_driving=soc_range];

    soc0_byaction = repeat(reshape(soc_range, 1, 1, FF, 1), PP, EE, 1, FF);
    soc1_byaction = soc0_byaction .+ dsoc;

    # STEP 1: Calculate V[S] under every scenario
    soc_needed = soc_scheduled(dt0 + periodstep(SS))
    vehicle_split = split_below.(soc_range, soc_needed)
    value_energy_bysoc = [value_energy(vehicle_split[ff][1], vehicle_split[ff][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff=1:FF]
    VV2 = repeat(reshape(value_energy_bysoc, EE, FF, 1), 1, 1, FF)

    # STEP 2: Determine optimal action for t = S-1 and back
    for tt in (SS-1):-1:1
        println(tt)

        dt1 = dt0 + periodstep(tt)
        price = get_retail_price(dt1)
        valuep = [value_power_action(price, dsoc[pp, ee, ff1, ff2], vehicles_plugged_range[ee]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        soc_needed = soc_scheduled(dt1)
        vehicle_split = split_below.(soc_range, soc_needed);
        valuepns = [value_power_newstate(price, vehicle_split[ff12][1], soc_needed - vehicle_split[ff12][2], vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF];
        valuee = [value_energy(vehicle_split[ff12][1], vehicle_split[ff12][3], soc_needed, vehicles_plugged_range[ee]) for ee=1:EE, ff12=1:FF];

        ff12_byaction = discrete_roundbelow.(soc1_byaction, soc_min, soc_max, FF);
        valuepns_byaction = [valuepns[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];
        valuee_byaction = [valuee[ee, ff12_byaction[pp, ee, ff1, ff2]] for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

        VV1byactsummc = zeros(Float64, PP, EE, FF, FF);

        for mc in 1:mcdraws
            if mcdraws == 1
                simustep = get_simustep_deterministic(dt1, drive_starts_time, park_starts_time)
            else
                simustep = get_simustep_stochastic(dt1, drive_starts_time, park_starts_time)
            end
            # Note: We impose costs from soc-below vehicles, but do not adjust state because it pushes up plugged-in soc every period
            statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]),
                                               vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][3], soc_range[ff2]),
                                      vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

            state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2);

            VV1byactthismc = combinebyact(VV2, state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3)
            VV1byactsummc += VV1byactthismc;
        end

        VV1byact = VV1byactsummc / mcdraws + valuep + valuepns_byaction + valuee_byaction;
        VV1byact[isnan.(VV1byact)] .= -Inf

        bestact = dropdims(argmax(VV1byact, dims=1), dims=1);
        strat[tt, :, :, :] .= Base.Fix2(getindex, 1).(bestact);
        if any(strat[tt, :, :, :] .== 0)
            break
        end

        VV2 = VV1byact[bestact]
        VV2[isnan.(VV2)] .= -Inf
    end

    return strat, VV2
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