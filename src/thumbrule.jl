using Dates
function get_dsoc_thumbrule1(tt, state, drive_starts_time, floor, ceiling)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)
    ## find when the users start to drive
    if drive_starts_time < Dates.Time(dt1)
        dt_drive = DateTime(Dates.Date(dt1)+ Dates.Day(1), drive_starts_time) ## assumes you start driving the next day
    else
        dt_drive = DateTime(Dates.Date(dt1), drive_starts_time) ## assumes you start driving later that day
    end 

    time_available = (dt_drive - dt1).value / 3600 / 1000 # in hours

    # Calculate floor SOC needed to ramp up to 80% by 9am (should incorporate the time available)

    charge_needed = max(0,0.80 - soc_plugged_1)
    time_required = charge_needed / (fracpower_max * timestep)
    SOC_floor = max(floor, soc_plugged_1 - time_required * fracpower_max) ## Charge can't go below 0.3
    
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
        soc_goal = ceiling # Charge if off peak
    elseif SOC_floor > soc_plugged_1
        soc_goal = 0.80 # Charge if below floor
    end

    ## should I add this?
    # if time_available < time_required # if not enought time, set goal to 0.8 to ramp up to drive time
    #     soc_goal = 0.8
    # end 



    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end

function get_dsoc_thumbrule_baseline(tt,state)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)
    soc_goal = 0.8
    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end 

function calculate_no_arbitrage_with_ramp(tt, state, drive_starts_time, floor, ceiling)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)
    
    if drive_starts_time < Dates.Time(dt1)
        dt_drive = DateTime(Dates.Date(dt1)+ Dates.Day(1), drive_starts_time) ## assumes you start driving the next day
    else
        dt_drive = DateTime(Dates.Date(dt1), drive_starts_time) ## assumes you start driving later that day
    end 

    # Calculate floor SOC needed to ramp up to 80% by drive time

    charge_needed = max(0,0.80 - soc_plugged_1)
    time_required = charge_needed / (fracpower_max * timestep)
    SOC_floor = soc_plugged_1 - time_required * fracpower_max ## Calculate ramp to get to 0.8 by drive time 
    
    time_available = (dt_drive - dt1).value / 3600 / 1000 # in hours


    if time_available < time_required # if not enought time, set goal to 0.8
        soc_goal = 0.8
    else
        soc_goal = (ceiling + floor) / 2
    end 
         
    soc_goal = max(soc_goal, SOC_floor) ## if the floor is above the goal, set the goal to the floor

    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end 

function calculate_middle_charging_band_rot(tt, state)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)

    middle_of_band = fracpower_max*timestep
    soc_goal = 0.95 - middle_of_band
    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)

end

function thumbrule_regrange(drive_starts_time, park_starts_time)
    vehicles_plugged_1 = 4.
    
    # ## allow energy arbitrage so long as there is space in the middle of your charging rate band to fluctuate as opposed to being fixed at 0.625 when at level 3 charging. 
    # ## 1. Calculate how much we can charge in 1 period. 
    # fracpower_max
    ## 2. Get the range left for ROT1 as 0.3 + maxcharge to 0.95 - maxcharge .
    range_left = (0.95 - fracpower_max)*timestep * vehicles_plugged_1 * vehicle_capacity - (0.3 + fracpower_min)*timestep * vehicles_plugged_1 * vehicle_capacity 

    regrange_func = (tt) -> begin
        current_time = Dates.Time(tt % 24, 0, 0)
        if drive_starts_time <= park_starts_time
            # Don't park on the next day
            if current_time >= drive_starts_time && current_time <= park_starts_time
                return 0.0  # During drive hours
            else
                return regrange_value  # Outside drive hours
            end
        else
            # Midnight crossing case
            if current_time >= drive_starts_time || current_time <= park_starts_time
                return 0.0  # During drive hours
            else
                return regrange_value  # Outside drive hours
            end
        end
    end

    ## 3. If that’s a negative range, then just do ROT1 with a 0.625 target.
    ## If it's a positive range, allow arbitrage within the range and regrange + frac_power_max up to 0.95 - frac_power_min to 0.3 

    ## need to adjust plan for regrange_value conditional on being plugged in (outside of drive_starts_time and park_starts_time)
    if range_left < 0 
        regrange_value = (0.95 - 0.3) * timestep * vehicles_plugged_1 * vehicle_capacity
        df = fullsimulate(dt0, (tt, state) -> calculate_no_arbitrage_with_ramp(tt, state, drive_starts_time, 0.3, 0.95), (tt) -> regrange_value, vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)

    else 
        ceiling = 0.95 - fracpower_max 
        floor = 0.3 + fracpower_min 
        regrange_value = (fracpower_max- fracpower_min) * timestep * vehicles_plugged_1 * vehicle_capacity 
        df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, floor, ceiling), regrange_func, vehicles_plugged_1,  0.5, 0.5, drive_starts_time, park_starts_time)

    end


    ## 4. Offer as regrange min(0.95 - SOC, SOC - 0.3). 


    benefits = sum(df[!, "valuep"]) + sum(df[!, "valuer"])
    return benefits
end 