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

dt0 = DateTime("2023-07-17T12:00:00")
drive_starts_time = Dates.Time(7, 0, 0)
park_starts_time = Dates.Time(17, 0, 0)

df = fullsimulate(dt0, (tt, state) -> get_dsoc(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time)
benefits = sum(df[!, "valuep"])
print(benefits)

plot_standard(df)
plot!(size=(700,400))
savefig("version1-rot2.pdf")