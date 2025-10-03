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
include("../src/plotting.jl")

dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)

regrange_vals = [regrange_func(tt) for tt in 1:SS] /  vehicles / vehicle_capacity

df = fullsimulate(dt0, dsoc_func, regrange_func, 4., 0.5, 0.5, drive_starts_time, park_starts_time)

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=regrange_vals, label="Reg. Range")

df = fullsimulate(dt0, dsoc_func, regrange_func, 4., 0.9, 0.9, drive_starts_time, park_starts_time)

# Collect regrange values for each timestep
regrange_vals = [regrange_func(tt) for tt in 1:SS] /  vehicles / vehicle_capacity

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=regrange_vals, label="Reg. Range")

max_charging_kw = 50 # higher charging rate to illustrate thumbrule
fracpower_min = -max_charging_kw / vehicle_capacity # discharge in terms of fraction of energy
fracpower_max = max_charging_kw / vehicle_capacity # charging in terms of fraction of energy

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)

df = fullsimulate(dt0, dsoc_func, regrange_func, 4., 0.7, 0.7, drive_starts_time, park_starts_time)

# Collect regrange values for each timestep
regrange_vals = [regrange_func(tt) for tt in 1:SS] /  vehicles / vehicle_capacity

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=regrange_vals, label="Reg. Range")
