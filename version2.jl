## Solution concept:
# Return both standard V(state) and expected profit to aggregator (state)
# Impose different levels of regulation-available range and solve within these
# Choose whichever provides the maximum expected profit in the known initial state
# Current version assumes available always except 8am - 6pm; future versions could impose more complicated timeseries of reg-available

include("src/bizutils.jl")
include("src/customer.jl")
include("src/simulate.jl")
include("src/retail.jl")
include("src/config.jl")
include("src/value.jl")
include("src/optutils.jl")
include("src/fullsim.jl")
include("src/plotting.jl")
include("src/optfuncs.jl")

energy_min = 0.
energy_max = vehicle_capacity * vehicles

RR = 5 # number of possible regrange values

probfail_penalty = 10.

dt0 = DateTime("2023-07-17T00:00:00")
soc_plugged_1 = 0.5
soc_driving_1 = 0.5
mcdraws = 1

drive_starts_time = Dates.Time(9, 0, 0)
park_starts_time = Dates.Time(17, 0, 0)

vehicles_plugged_1 = find_starting_vehicles_plugged(dt0, drive_starts_time, park_starts_time)


probstate = optimize_regrange_probstate_outer_loop(dt0, soc_plugged_1, soc_driving_1, vehicles_plugged_1, drive_starts_time, park_starts_time)
strat, probfail, optregrange = optimize_regrange_probstate(dt0, probstate, drive_starts_time, park_starts_time);

df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
df[!, :optregrange] = [0.; optregrange]

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=df.optregrange / (vehicles * vehicle_capacity), label="Reg. Range")
plot!(size=(700,400))
savefig("optregrange-v1.pdf")
