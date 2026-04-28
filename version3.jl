## Solution concept:
# Return both standard V(state) and expected profit to aggregator (state)
# Impose different levels of regulation-available range and solve within these
# Choose whichever provides the maximum expected profit in the known initial state
# Current version assumes available always except 8am - 6pm; future versions could impose more complicated timeseries of reg-available

using Random

include("src/bizutils.jl")
include("src/customer.jl")
include("src/simulate.jl")
include("src/retail.jl")
include("src/config.jl")
include("src/value.jl")
include("src/optutils.jl")
include("src/fullsim.jl")
include("src/plotting.jl")

energy_min = 0.
energy_max = vehicle_capacity * vehicles

probfail_penalty = 100.
portion_below_penalty = 1000.

dt0 = DateTime("2023-07-17T12:00:00")
vehicles_plugged_1 = 0.
soc_plugged_1 = 0.5
soc_driving_1 = 0.5
mcdraws = 1

stratbest, regrangebest = optimize_regrange_given_outer_loop(dt0, soc_plugged_1, soc_driving_1, vehicles_plugged_1, drive_starts_time, park_starts_time)

df = fullsimulate(dt0, stratbest, regrangebest, vehicles_plugged_1, soc_plugged_1, 0.)
df[!, :optregrange] = [0.; regrangebest[1:end-1]]

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=(regneutral / 2) * df.optregrange / (vehicles * vehicle_capacity), label="Regulation Range")
plot!(size=(700,400))
savefig("optregrange-v2.pdf")
