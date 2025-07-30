using Plots
using Dates

include("../src/customer.jl")
include("../src/config.jl")
include("../src/optfuncs.jl")
include("../src/fullsim.jl")

dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Dates.Time(9, 0, 0)
park_starts_time = Dates.Time(17, 0, 0)

@time strat, VVall = optimize(dt0, SS, drive_starts_time, park_starts_time);

df = fullsimulate(dt0, strat, zeros(SS-1), 4., 0.5, 0.5, drive_starts_time, park_starts_time)
benefits = sum(df.valuep[1:24] + df.valuepns[1:24])

results = DataFrame(drive_starts_hour=Int64[], park_starts_hour=Int64[], value=Union{Float64, Missing}[])
for drive_starts_hour in 1:22
    for park_starts_hour in 1:23
        if park_starts_hour ≤ drive_starts_hour
            push!(results, [drive_starts_hour, park_starts_hour, missing])
            continue
        end

        drive_starts_time = Dates.Time(drive_starts_hour, 0, 0)
        park_starts_time = Dates.Time(park_starts_hour % 24, 0, 0)

        soc_dispersion = 0.05

        @time strat, VVall = optimize(dt0, SS, drive_starts_time, park_starts_time);

        soc_dispersion = 0.0

        df = fullsimulate(dt0, strat, zeros(SS-1), (drive_starts_hour == 1) ? 0. : 4., 0.5, 0.5, drive_starts_time, park_starts_time)
        benefits = sum(df.valuep[1:24] + df.valuepns[1:24])

        push!(results, [drive_starts_hour, park_starts_hour, benefits])
    end
end

function dataframe_to_matrix(df, xcol, ycol, vcol)
    xs = sort(unique(df[!, xcol]))
    ys = sort(unique(df[!, ycol]))
    matrix = [df[findall((df[!, xcol] .== x) .& (df[!, ycol] .== y)), vcol][1] for y in ys, x in xs]
    return matrix
end

results2 = dataframe_to_matrix(results, :drive_starts_hour, :park_starts_hour, :value)
heatmap(1:23, 1:24, results2, xlabel="Drive time", ylabel="Park time")
plot!(size=(600,400))
savefig("../version1-mv1.png")
