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
vehicles_plugged_1 = 4.
soc_plugged_1 = 0.5
soc_driving_1 = 0.5
mcdraws = 1

drive_starts_time = Dates.Time(9, 0, 0)
park_starts_time = Dates.Time(17, 0, 0)

probstate = zeros(SS-1, EE, FF, FF);
df = fullsimulate(dt0, zeros(SS-1), vehicles_plugged_1, soc_plugged_1, soc_driving_1, drive_starts_time, park_starts_time)
for ii in 1:(nrow(df) - 1)
    statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3 = breakstate((df.vehicles_plugged[ii], df.soc_plugged[ii], df.soc_driving[ii]))
    probstate[ii, statebase...] = (probbase1 + probbase2 + probbase3) / 3
    probstate[ii, makeindex1(statebase, stateceil1)] += (1 - probbase1) / 3
    probstate[ii, makeindex2(statebase, stateceil2)] += (1 - probbase2) / 3
    probstate[ii, makeindex3(statebase, stateceil3)] += (1 - probbase3) / 3
end

probstate = probstate / 2 .+ repeat(ones(1, EE, FF, FF) / (EE * FF * FF), SS-1) / 2;

for ll in 1:10
    println("Loop $(ll)")
    strat, probfail, optregrange = optimize_regrange(dt0, probstate, drive_starts_time, park_starts_time);

    dfall = nothing
    for ii in 1:mcdraws
        local df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
        energy = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) .* df.soc_above
        energy_minallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.3
        energy_maxallow = vehicle_capacity * df.vehicles_plugged .* (1 .- df.portion_below) * 0.95
        df[!, :regrange_avail] = min.(energy_maxallow - energy, energy - energy_minallow)

        if dfall == nothing
            dfall = df
        else
            append!(dfall, df)
        end
    end

    global probstate = zeros(SS-1, EE, FF, FF);
    for tt in 1:SS-1
        dt1 = dt0 + periodstep(tt - 1)

        for state in unique(dfall.state[dfall.datetime .== dt1])
            probstate[tt, state...] = sum(map(rowstate -> rowstate == state, dfall.state[dfall.datetime .== dt1])) / sum(dfall.datetime .== dt1)
        end
    end
end

strat, probfail, optregrange = optimize_regrange(dt0, probstate, drive_starts_time, park_starts_time);

df = fullsimulate(dt0, strat, optregrange, vehicles_plugged_1, soc_plugged_1, 0., drive_starts_time, park_starts_time)
df[!, :optregrange] = [0.; optregrange]

pp = plot_standard(df)
plot!(pp, df.datetime, df.soc_plugged, ribbon=df.optregrange / (vehicles * vehicle_capacity), label="Reg. Range")
plot!(size=(700,400))
savefig("optregrange-v1.pdf")
