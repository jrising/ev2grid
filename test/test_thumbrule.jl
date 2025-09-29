using Test
using Plots
using DataFrames

include("../src/bizutils.jl")
include("../src/customer.jl")
include("../src/simulate.jl")
include("../src/retail.jl")
include("../src/config.jl")
include("../src/value.jl")
include("../src/optutils.jl")
include("../src/fullsim.jl")
include("../src/thumbrule.jl")

dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)

[regrange_func(tt) for tt in 1:36]



df = fullsimulate(dt0, dsoc_func, regrange_func, 4., 0.5, 0.5, drive_starts_time, park_starts_time)

# Collect regrange values for each timestep
regrange_vals = [regrange_func(tt) for tt in 1:size(df, 1)] /  4 / vehicle_capacity

regrange_upper = min.(soc_max * ones(length(df.soc_plugged)), df.soc_plugged .+ 0.5 .* regrange_vals)
regrange_upper[1] = df.soc_plugged[1]  # Ensure the first value is equal to soc_plugged at the start
regrange_lower = regrange_upper .- regrange_vals

# Plot soc_driving and soc_plugged
plot(df.datetime, df.soc_driving, label="SOC Driving", lw=2)
plot!(df.datetime, df.soc_plugged, label="SOC Plugged", lw=2)

# Plot regrange bounds as dashed lines
plot!(df.datetime, regrange_upper, label="Regrange Upper", lw=2, linestyle=:dash)
plot!(df.datetime, regrange_lower, label="Regrange Lower", lw=2, linestyle=:dash)

xlabel!("Time")
ylabel!("Value")
title!("SOC and Regrange Bounds Over Simulation")
savefig("plot_thumbrule_1.png")

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)

print([regrange_func(tt) for tt in 1:36])

df = fullsimulate(dt0, dsoc_func, regrange_func, 4., 0.9, 0.9, drive_starts_time, park_starts_time)

# Collect regrange values for each timestep
regrange_vals = [regrange_func(tt) for tt in 1:size(df, 1)] /  4 / vehicle_capacity

regrange_upper = min.(soc_max * ones(length(df.soc_plugged)), df.soc_plugged .+ 0.5 .* regrange_vals)
regrange_upper[1] = df.soc_plugged[1]  # Ensure the first value is equal to soc_plugged at the start
regrange_lower = regrange_upper .- regrange_vals

# Plot soc_driving and soc_plugged
plot(df.datetime, df.soc_driving, label="SOC Driving", lw=2)
plot!(df.datetime, df.soc_plugged, label="SOC Plugged", lw=2)

# Plot regrange bounds as dashed lines
plot!(df.datetime, regrange_upper, label="Regrange Upper", lw=2, linestyle=:dash)
plot!(df.datetime, regrange_lower, label="Regrange Lower", lw=2, linestyle=:dash)

xlabel!("Time")
ylabel!("Value")
title!("SOC and Regrange Bounds Over Simulation")
savefig("plot_thumbrule_2.png")

print([regrange_func(tt) for tt in 1:36])

max_charging_kw = 50 # higher charging rate to illustrate thumbrule
fracpower_min = -max_charging_kw / vehicle_capacity # discharge in terms of fraction of energy
fracpower_max = max_charging_kw / vehicle_capacity # charging in terms of fraction of energy

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)

df = fullsimulate(dt0, dsoc_func, regrange_func, 4., 0.7, 0.7, drive_starts_time, park_starts_time)

# Collect regrange values for each timestep
regrange_vals = [regrange_func(tt) for tt in 1:size(df, 1)] /  4 / vehicle_capacity

regrange_upper = min.(soc_max * ones(length(df.soc_plugged)), df.soc_plugged .+ 0.5 .* regrange_vals)
regrange_upper[1] = df.soc_plugged[1]  # Ensure the first value is equal to soc_plugged at the start
regrange_lower = regrange_upper .- regrange_vals

# Plot soc_driving and soc_plugged
plot(df.datetime, df.soc_driving, label="SOC Driving", lw=2)
plot!(df.datetime, df.soc_plugged, label="SOC Plugged", lw=2)

# Plot regrange bounds as dashed lines
plot!(df.datetime, regrange_upper, label="Regrange Upper", lw=2, linestyle=:dash)
plot!(df.datetime, regrange_lower, label="Regrange Lower", lw=2, linestyle=:dash)

xlabel!("Time")
ylabel!("Value")
title!("SOC and Regrange Bounds Over Simulation")
savefig("plot_thumbrule_3.png")
