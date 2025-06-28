using Dates
function get_dsoc_thumbrule1(tt, state, drive_starts_time)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)
    ## find when the users start to drive
    if drive_starts_time < Dates.Time(dt1)
        dt_drive = DateTime(Dates.Date(dt1)+ Dates.Day(1), drive_starts_time) ## assumes you start driving the next day
    else
        dt_drive = DateTime(Dates.Date(dt1), drive_starts_time) ## assumes you start driving later that day
    end 

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
         
    soc_goal = 0.8 ## if nothing else set to 80% 

    if is_peak(dt1) && SOC_floor <= soc_plugged_1 && price_switch
        soc_goal = SOC_floor # Discharge if peak price, above floor, and price change coming
    elseif !is_peak(dt1)
        soc_goal = 0.95 # Charge if off peak
    elseif SOC_floor > soc_plugged_1
        soc_goal = 0.80 # Charge if below floor
    end


    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end

function get_dsoc_thumbrule_baseline(tt,state,drive_starts_time)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)
    soc_goal = 0.8
    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end 