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

if false
    df = DataFrame(hourstart=DateTime("2024-01-01T00:00:00"):Hour(1):DateTime("2024-12-31T00:00:00"))
    df[!, :price] = get_retail_price.(df.hourstart)

    pp = plot(df.hourstart, df.price, seriestype=:steppost, label="")
    plot!(pp, size=(1000, 400))
    savefig("retailprice.png")
end
