using DataFrames

include("src/bizutils.jl")
include("src/customer.jl")
include("src/simulate.jl")
include("src/retail.jl")
include("src/config.jl")
include("src/value.jl")
include("src/optutils.jl")
include("src/fullsim.jl")
include("src/plotting.jl")

function get_dsoc(tt, state)
    vehicles_plugged_1, soc_plugged_1, soc_driving_1 = state
    dt1 = dt0 + periodstep(tt)

    if is_peak(dt1)
        soc_goal = 0.8
    else
        soc_goal = 0.95
    end

    return max(min(soc_goal - soc_plugged_1, timestep * fracpower_max), timestep * fracpower_min)
end

dt0 = DateTime("2023-07-17T12:00:00")

df = fullsimulate(dt0, get_dsoc, (tt) -> 0., 0., 0.5, 0.5)
benefits = sum(df[!, "valuep"])
plot_standard(df)
plot!(size=(700,400))
savefig("version1-rot.pdf")
