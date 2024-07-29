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
