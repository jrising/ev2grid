# import Pkg; Pkg.add("DataFrames")
# import Pkg; Pkg.add("Roots")
# import Pkg; Pkg.add("Distributions")
# import Pkg; Pkg.add("StatsBase")
# import Pkg; Pkg.add("HolidayCalendars")
# import Pkg; Pkg.add("RDates")
# import Pkg; Pkg.add("ArgCheck")
# import Pkg; Pkg.add("CSV")
# import Pkg; Pkg.add("Plots")
using DataFrames

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

    # Check if it is a peak period 


    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end

dt0 = DateTime("2023-07-17T12:00:00")
df = fullsimulate(dt0, (tt, state) -> get_dsoc(tt, state), (tt) -> 0., 0., 0.5, 0.5)
benefits = sum(df[!, "valuep"])
print(benefits)

plot_standard(df)
plot!(size=(700,400))
savefig("version1-rot2.pdf")