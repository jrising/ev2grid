using Dates
using HolidayCalendars, RDates

"""
    is_peak(hourstart::DateTime) -> Bool

Determine whether a given datetime falls within peak energy cost hours.

# Arguments
- `hourstart::DateTime`: The datetime to check, given as the start of an hour.

# Returns
- Bool: `true` if the datetime is within on-peak hours on a non-holiday weekday, otherwise `false`.

# Examples
```julia
is_peak(DateTime("2024-07-29T13:00:00")) # returns true
is_peak(DateTime("2024-07-29T09:00:00")) # returns false
"""
function is_peak(hourstart::DateTime)
    # Define the start and end times for on-peak hours
    peak_start = Time(12, 0)
    peak_end = Time(20, 0)

    # Check if the date is a weekend or a weekday holiday
    weekday = dayofweek(hourstart)
    if weekday in [6, 7]  # Saturday (6) or Sunday (7)
        return false
    end

    # List of holidays
    cal = calendar(CALENDARS, "NEW YORK")
    if is_holiday(cal, Date(hourstart))
        return false
    end

    # Check the time
    if Time(hourstart) ≥ peak_start && Time(hourstart) < peak_end
        return true
    else
        return false
    end
end

"""
    get_retail_price(hourstart::DateTime) -> Float64

Get the retail electricity price for a given datetime.

# Arguments
 - hourstart::DateTime: The datetime to check.

# Returns
 - Float64: The electricity price per kWh based on season and peak hours.

# Examples
get_retail_price(DateTime("2024-07-29T13:00:00")) # returns 0.1473
get_retail_price(DateTime("2024-01-01T02:00:00")) # returns 0.08356
"""
function get_retail_price(hourstart::DateTime)
    # Define the months for summer and winter seasons
    summer_months = 6:9  # June through September
    winter_months = [10:12; 1:5]  # October through May

    if month(hourstart) in summer_months
        return is_peak(hourstart) ? 0.1473 : 0.07242
    else
        return is_peak(hourstart) ? 0.1720 : 0.08356
    end
end

"""
    get_retail_price(dt0::DateTime, kwbytimestep::Vector{Float64}, timestep::Float64) -> Float64

Estimate the demand pricing cost for Delmarva Power, given a vector a kW values all within a given month.

# Arguments
 - dt0::DateTime: The time of the first kW timestep value.
 - kwbytimestep::Vector{Float64}: A vector of kW used during each timestep.
 - timestep::Float64: The length of a timestep, in hours (typically a value of 1. or less).

# Returns
 - Float64: The demand pricing cost, according to this data.
"""
function get_demand_cost(dt0::DateTime, kwbytimestep::Vector{Float64}, timestep::Float64)
    max60_peak = 0.
    max60_nonpeak = 0.
    timestepsper60 = round(Int, 1. / timestep)
    for tt in 1:(length(kwbytimestep) - timestepsper60 + 1)
        max60 = maximum(kwbytimestep[tt:(tt + timestepsper60 - 1)])
        if is_peak(dt0 + Hour(round(Int, timestep * (tt - 1))))
            max60_peak = max(max60, max60_peak)
        else
            max60_nonpeak = max(max60, max60_nonpeak)
        end
    end

    18.73 * max(round(max60_peak), round(max60_nonpeak / 3))
end
