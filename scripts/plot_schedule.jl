using Dates, Plots, DataFrames

include("../src/bizutils.jl")
include("../src/customer.jl")

timestep = 1. # 1 hour
df = DataFrame(dt=DateTime("2023-06-01T00:00:00"):Hour(1):DateTime("2023-06-03T00:00:00"))
df[!, :soc] = soc_scheduled.(df.dt, Dates.Time(9, 0, 0))

pp = plot(df.dt, df.soc, seriestype=:steppost, label="")
plot!(pp, size=(1000, 400))
savefig("schedule.png")

include("../src/config.jl")
include("../src/optutils.jl")
include("../src/retail.jl")
include("../src/value.jl")
include("../src/simulate.jl")
include("../src/fullsim.jl")

SS = nrow(df) - 1
df2 = fullsimulate(df.dt[1] - periodstep(1), zeros(SS-1), 4., 0.7, 0.7, false)

pp = plot(df2.datetime, df2.vehicles_plugged / 4, seriestype=:steppost, label="Vehicles Plugged-In")
plot!(pp, df2.datetime, df2.soc_needed, seriestype=:steppost, label="Energy Needed")
plot!(pp, df2.datetime, (df2.soc_plugged .* df2.vehicles_plugged .+ df2.soc_driving .* (vehicles .- df2.vehicles_plugged)) / vehicles, seriestype=:steppost, label="Energy of Plugged-In")
plot!(pp, size=(1000, 400))
savefig("schedule.png")

baseline = copy(df2)

pp = plot(baseline.datetime, baseline.vehicles_plugged / 4, seriestype=:steppost, label="Vehicles Plugged-In")

for mc in 1:100
    local df2 = fullsimulate(df.dt[1] - periodstep(1), zeros(SS-1), 4., 0.7, 0.7, true)
    if all(df2.vehicles_plugged .== baseline.vehicles_plugged)
        continue
    end

    plot!(pp, df2.datetime, df2.vehicles_plugged / 4, seriestype=:steppost, la=.25, legend=false)
end

plot!(pp, baseline.datetime, baseline.vehicles_plugged / 4, seriestype=:steppost, linewidth=5, label="Vehicles Plugged-In")
pp
