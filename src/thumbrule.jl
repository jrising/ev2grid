using Dates

function get_dsoc_thumbrule1(tt, state, drive_starts_time, floor, ceiling, drive_time_charge_level)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)

    time_available = hours_to_drive(dt1, drive_starts_time)

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

    soc_goal = drive_time_charge_level ## if nothing else set to 80%

    if is_peak(dt1) && floor <= soc_plugged_1 && price_switch
        soc_goal = floor # Discharge if peak price, above floor, and price change coming
    elseif !is_peak(dt1)
        soc_goal = ceiling # Charge if off peak
    elseif floor > soc_plugged_1
        soc_goal = drive_time_charge_level # Charge if below floor
    end

    ## add safe_charging as a fallback
    safe_charging(dt1, drive_starts_time, soc_goal, floor, drive_time_charge_level)
end

function get_dsoc_thumbrule_baseline(tt, state, drive_time_charge_level)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)
    soc_goal = drive_time_charge_level
    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end


function safe_charging(dt, drive_starts_time, soc_goal, soc_floor, drive_time_charge_level)
    soc_floor = safe_soc_floor(dt, drive_starts_time, soc_floor, drive_time_charge_level)
    soc_goal = max(soc_goal, soc_floor) ## if the floor is above the goal, set the goal to the floor

    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)

end

function thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)
    vehicles_plugged_1 = vehicles_plugged_scheduled(dt0, drive_starts_time, park_starts_time)

    # ## allow energy arbitrage so long as there is space in the middle of your charging rate band to fluctuate as opposed to being fixed at 0.625 when at level 3 charging.
    # ## 1. Calculate how much we can charge in 1 period.
    # fracpower_max
    ## 2. Get the range left for ROT1 as 0.3 + maxcharge to 0.95 - maxcharge .
    ceiling = soc_max - fracpower_max
    floor = soc_min + fracpower_min

    ## consider the case where the regrange overlaps the edges of the battery
    if ceiling < drive_time_charge_level
        ## figure out how many hours of buffer is needed
        n = (drive_time_charge_level - ceiling) / fracpower_max
        buffer = ceil(n)
    else
        buffer = 0
    end

    regrange_func = (tt) -> begin
        dt1 = dt0 + periodstep(tt)
        current_time = Time(dt1)
        if drive_starts_time <= park_starts_time
            # Don't park on the next day
            if current_time >= drive_starts_time && current_time <= park_starts_time - Hour(buffer)
                return 0.0  # During drive hours
            else
                return regrange_value  # Outside drive hours
            end
        else
            # Midnight crossing case
            if current_time >= drive_starts_time || current_time <= park_starts_time - Hour(buffer)
                return 0.0  # During drive hours
            else
                return regrange_value  # Outside drive hours
            end
        end
    end

    ## 3. If that’s a negative range, then just do ROT1 with a 0.625 target.
    ## If it's a positive range, allow arbitrage within the range and regrange + frac_power_max up to 0.95 - frac_power_min to 0.3

    regrange_value = min(soc_max - soc_min, fracpower_max - fracpower_min) * timestep * vehicles_plugged_1 * vehicle_capacity

    if ceiling < floor
        ## no arbitrage case, focus on reg services
        soc_goal = (soc_max- soc_min) / 2
        df = fullsimulate(dt0, (tt, state) -> safe_charging(tt, state, drive_starts_time, soc_goal, floor, drive_time_charge_level), (tt) -> regrange_value, vehicles_plugged_1, 0.5, 0.5, drive_starts_time, park_starts_time)
    else
        ## include arbitrage between ceiling and floor
        df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time, floor, ceiling, drive_time_charge_level), regrange_func, vehicles_plugged_1,  0.5, 0.5, drive_starts_time, park_starts_time)
    end


    ## 4. Offer as regrange min(0.95 - SOC, SOC - 0.3).

    benefits = sum(df[!, "valuep"]) + sum(df[!, "valuer"])
    return benefits
end

function hours_to_drive(dt, drive_starts_time)
    if drive_starts_time < Time(dt)
        dt_drive = DateTime(Date(dt) + Day(1), drive_starts_time) ## assumes you start driving the next day
    else
        dt_drive = DateTime(Date(dt), drive_starts_time) ## assumes you start driving later that day
    end

    return (dt_drive - dt).value / Hour(1)

function safe_soc_floor(dt, drive_starts_time, soc_floor, drive_time_charge_level)
    timesteps_available = hours_to_drive(dt, drive_starts_time) / timestep
    return max(soc_floor, drive_time_charge_level - fracpower_max * (timesteps_available - 1))
end
