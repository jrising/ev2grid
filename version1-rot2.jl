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
include("src/thumbrule.jl")


dt0 = DateTime("2023-07-17T12:00:00")
drive_starts_time = Dates.Time(7, 0, 0)
park_starts_time = Dates.Time(17, 0, 0)

df = fullsimulate(dt0, (tt, state) -> get_dsoc_thumbrule1(tt, state, drive_starts_time), (tt) -> 0., 0., 0.5, 0.5, drive_starts_time, park_starts_time)
benefits = sum(df[!, "valuep"])
print(benefits)

plot_standard(df)
plot!(size=(700,400))
savefig("version1-rot2.pdf")