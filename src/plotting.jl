using Plots

function plot_standard(df::DataFrame, dropfromend=12)
    df = df[1:nrow(df) - dropfromend, :]
    df.soc_plugged[(df.vehicles_plugged .== 0) .& [false; df.vehicles_plugged[1:end-1] .== 0]] .= NaN
    df.soc_driving[(df.vehicles_plugged .== vehicles) .& [false; df.vehicles_plugged[1:end-1] .== vehicles]] .= NaN
    pp = plot(df.datetime, df.vehicles_plugged ./ vehicles, seriestype=:line, label="Vehicles plugged-in")
    plot!(pp, df.datetime, df.soc_plugged, seriestype=:line, label="Plugged-in SOC")
    plot!(pp, df.datetime, df.soc_driving, seriestype=:line, label="Driving SOC")
    pp
end
