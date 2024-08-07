using DataFrames
using Plots

SS = 36 # project for 1.5 days
timestep = 1. # 1 hour
include("utils.jl")
include("simulate.jl")
include("src/customer.jl")
include("fullsim.jl")

dt0 = DateTime("2024-07-15T12:00:00")
vehicles_plugged_1 = 4.
enerfrac_plugged_1 = 0.5
enerfrac_driving_1 = 0.5
df = simu_inactive(dt0, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)

pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label="")

## TODO: Need a strat for this

vehicles_plugged_1 = 4.

pp = nothing
for enerfrac_plugged_1 in range(enerfrac_min, enerfrac_max, FF-1)
    local enerfrac_driving_1 = enerfrac_plugged_1
    local df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)
    global pp

    if pp == nothing
        pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1, legend=false)
    else
        plot!(pp, df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1, legend=false)
    end
end

pp
