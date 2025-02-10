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

function get_dsoc(tt, state)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)

    # Calculate floor SOC needed to ramp up to 80% by 9am
    dt_drive = DateTime(Dates.Date(dt1)+ Dates.Day(1), Dates.Time(9, 0, 0))
    time_available = (dt_drive - dt1).value / 3600 / 1000 # in hours

    charge_needed = max(0,0.95 - soc_plugged_1)
    time_required = charge_needed / (fracpower_max * timestep)
    SOC_floor = max(0,soc_plugged_1 - time_required * fracpower_max)
    
    # Set goal SOC

    soc_goal = 0.95 # Set a default to charge to 0.95 if no other conditions are met

    if is_peak(dt1) && SOC_floor <= soc_plugged_1 
        soc_goal = SOC_floor # Discharge if peak price and above floor
    elseif !is_peak(dt1) && soc_plugged_1 < 0.95
        soc_goal = 0.95 # Charge if off peak and below 95%
    elseif SOC_floor > soc_plugged_1
        soc_goal = 0.95 # Charge if below floor
    end

    # Check if time available is less than time required 
    if time_available < time_required
        soc_goal = min(soc_goal, soc_plugged_1 + time_available * fracpower_max)
    end

    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end

function run_rule_of_thumb_simulation()
    dt0 = DateTime("2023-07-17T12:00:00")
    df = fullsimulate(dt0, (tt, state) -> get_dsoc(tt, state), (tt) -> 0., 0., 0.5, 0.5)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function run_optimized_simulation(SS)
    dt0 = DateTime("2023-07-17T12:00:00")
    mcdraws = 1
    @time strat, VV = optimize(dt0, SS);


    df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function run_optimized_stochastic_simulation(SS, mcdraws)
    dt0 = DateTime("2023-07-17T12:00:00")
    mcdraws = 100
    @time strat, VV = optimize(dt0, SS);
    
    df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5)
    benefits = sum(df[!, "valuep"])

    return benefits
end

function export_to_latex(benefits_dict)
    filename = "benefits_table.tex"
    open(filename, "w") do io
        write(io, "\\begin{table}[h!]\n")
        write(io, "\\centering\n")
        write(io, "\\begin{tabular}{|c|c|}\n")
        write(io, "\\hline\n")
        write(io, "Simulation Type & Benefits \\\\ \n")
        write(io, "\\hline\n")
        for (sim_type, benefits) in benefits_dict
            write(io, @sprintf("%s & \$%.2f\$ \\\\ \n", sim_type, benefits))
        end
        write(io, "\\hline\n")
        write(io, "\\end{tabular}\n")
        write(io, "\\caption{Benefits for different simulation types}\n")
        write(io, "\\end{table}\n")
    end
end

function optimize(dt0::DateTime, SS::Int)
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
                simustep = get_simustep_deterministic(dt1)
            else
                simustep = get_simustep_stochastic(dt1)
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


SS = 144  # Example value, you can change this as needed
mcdraws = 100  # Example value for stochastic draws
dt0 = DateTime("2023-07-17T12:00:00")

benefits_dict = Dict{String, Float64}()

# Run simulations and collect benefits
benefits_dict["Rule of Thumb"] = run_rule_of_thumb_simulation()
benefits_dict["Optimized"] = run_optimized_simulation(SS)
benefits_dict["Optimized with Stochastic Draws"] = run_optimized_stochastic_simulation(SS, mcdraws)

# Export benefits to LaTeX
export_to_latex(benefits_dict)