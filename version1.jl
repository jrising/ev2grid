using Dates, DataFrames
using StatsBase, Random
using Plots
using HolidayCalendars, RDates
using ArgCheck

# Assume this is the day before, and we can observe current storage
# If there are any cars that need recharging to 30% level, handle that outside this system

## Order of events:
## In state1 at start of timestep
## Apply action throughout timestep
## Construct below-split and perform simulation
## Transition to state2 at end of timestep

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

soc_dispersion = 0.05

dt0 = DateTime("2023-07-17T00:00:00")
mcdraws = 1
drive_starts_time = Dates.Time(9, 0, 0)
park_starts_time = Dates.Time(17, 0, 0)

@time strat, VVall = optimize(dt0, SS, drive_starts_time, park_starts_time);

hourly_ticks = collect(DateTime("2023-07-17T00:00:00"):Hour(1):DateTime("2023-07-18T11:00:00"))
six_hour_ticks = hourly_ticks[1:6:length(hourly_ticks)]
hour_labels = string.(hour.(six_hour_ticks))
Plots.heatmap(transpose(max.(VVall[:, 5, :, 5], -10)); xticks = (1:6:length(hourly_ticks), hour_labels))

Plots.heatmap(transpose(max.(strat[:, 5, :, 5], -10)); xticks = (1:6:length(hourly_ticks), hour_labels))

df = fullsimulate(dt0, strat, zeros(SS-1), 4., 0.5, 0.5, drive_starts_time, park_starts_time)
benefits = sum(df[!, "valuep"])
plot_standard(df)
plot!(size=(700,400))
savefig("version1-det.pdf")

mcdraws = 100
@time strat, VV = optimize(dt0, SS, drive_starts_time, park_starts_time);

Plots.heatmap(transpose(max.(VV[:, 5, :, 5], -10)); xticks = (1:6:length(hourly_ticks), hour_labels))

df = fullsimulate(dt0, strat, zeros(SS-1), 4., 0.5, 0.5, drive_starts_time, park_starts_time)
benefits_sto = sum(df[!, "valuep"])
plot_standard(df)
plot!(size=(700,400))
savefig("version1-sto.pdf")

## How do the value parameters affect the penultimate level?
mcdraws = 1
results = DataFrame(weight_portion_above=Float64[], weight_portion_below=Float64[], ratio_exponent=Float64[], socend=Float64[])
for wp_above in range(0, 1., 5)
    for wp_below in range(0, 0.5, 5)
        for re in [.25, .5, 1, 2]
            global weight_portion_above = wp_above
            global weight_portion_below = wp_below
            global ratio_exponent = re
            strat, VV = optimize(dt0, SS, drive_starts_time, park_starts_time);
            df = fullsimulate(dt0, strat, zeros(SS-1), 0., 0.5, 0.5, drive_starts_time, park_starts_time)
            push!(results, [weight_portion_above, weight_portion_below, ratio_exponent, sum(df.soc_plugged .* df.vehicles_plugged) / sum(df.vehicles_plugged)])
        end
    end
end

function dataframe_to_matrix(df, xcol, ycol, vcol)
    xs = sort(unique(df[!, xcol]))
    ys = sort(unique(df[!, ycol]))
    matrix = [df[findall((df[!, xcol] .== x) .& (df[!, ycol] .== y)), vcol][1] for y in ys, x in xs]
    return matrix
end

value_matrix = dataframe_to_matrix(results[results.ratio_exponent .== 0.5, :], :weight_portion_above, :weight_portion_below, :socend)
heatmap(range(0, 1., 5), range(0, 0.5, 5), value_matrix, xlabel="Weight of above portion", ylabel="Weight of below portion")
plot!(size=(600,400))
savefig("version1-ves.png")

value_matrix = dataframe_to_matrix(results[results.ratio_exponent .== 1., :], :weight_portion_above, :weight_portion_below, :socend)
heatmap(range(0, 1., 5), range(0, 0.5, 5), value_matrix, xlabel="Weight of above portion", ylabel="Weight of below portion")
savefig("version1-ves2.png")
