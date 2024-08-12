using Dates, Plots, DataFrames

include("utils.jl")
include("src/customer.jl")

df = DataFrame(dt=DateTime("2024-06-01T00:00:00"):Hour(1):DateTime("2024-06-03T00:00:00"))
df[!, :enerfrac] = enerfrac_scheduled.(df.dt)

pp = plot(df.dt, df.enerfrac, seriestype=:steppost, label="")
plot!(pp, size=(1000, 400))
savefig("schedule.png")

include("simulate.jl")
include("fullsim.jl")

SS = nrow(df) - 1
df2 = simu_inactive(df.dt[1] - periodstep(1), 4., 0.7, 0.7, false)

pp = plot(df2.datetime, df2.vehicles_plugged / 4, seriestype=:steppost, label="Vehicles Plugged-In")
plot!(pp, df2.datetime, df2.enerfrac_needed, seriestype=:steppost, label="Energy Needed")
plot!(pp, df2.datetime, (df2.enerfrac_plugged .* df2.vehicles_plugged .+ df2.enerfrac_driving .* (vehicles .- df2.vehicles_plugged)) / vehicles, seriestype=:steppost, label="Energy of Plugged-In")
plot!(pp, size=(1000, 400))
savefig("schedule.png")

baseline = copy(df2)

## Unexpected changes in vehicles
prob_event = 0.01 # per hour, so 1 per 4 days
prob_event_vehicles = [0.5, .25, .125, .125]
prob_event_return = 0.25
prob_delayed_return = 0.1

pp = plot(baseline.datetime, baseline.vehicles_plugged / 4, seriestype=:steppost, label="Vehicles Plugged-In")

for mc in 1:100
    df2 = simu_inactive(df.dt[1] - periodstep(1), 4., 0.7, 0.7, true)
    if all(df2.vehicles_plugged .== baseline.vehicles_plugged)
        continue
    end

    println("OK!")
    plot!(pp, df2.datetime, df2.vehicles_plugged / 4, seriestype=:steppost, la=.25, legend=false)
end

plot!(pp, baseline.datetime, baseline.vehicles_plugged / 4, seriestype=:steppost, linewidth=5, label="Vehicles Plugged-In")
pp
