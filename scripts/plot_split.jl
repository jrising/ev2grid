using DataFrames
using Plots

include("src/bizutils.jl")

pp = plot([0.3, 0.3], [0, 1], seriestype=:path, color=:red, linewidth=2, label="Vertical Line")

df = DataFrame(soc_plugged=Float64[], portion_below=Float64[], soc_below=Float64[], soc_above=Float64[])
for soc_plugged in range(0., 1., 11)
    if soc_plugged >= 0.3 && soc_plugged <= 0.5
        f(mu) = mean(truncated(Normal(mu, 0.05), lower=0., upper=1.)) - soc_plugged
        mu = find_zero(f, (0., 1.))
        dist = truncated(Normal(mu, 0.05), lower=0., upper=1.)
        xx = range(0., 1., 1000)
        yy = pdf.(dist, xx)
        if soc_plugged == 0.3
            plot!(pp, xx, yy / maximum(yy), seriestype=:line, color=:gray, label="SOC Distribution")
        else
            plot!(pp, xx, yy / maximum(yy), seriestype=:line, color=:gray, label=false)
        end
    end

    portion_below, soc_below, soc_above = split_below(soc_plugged, 0.3)
    push!(df, [soc_plugged, portion_below, soc_below, soc_above])
end

plot!(pp, df.soc_plugged, df.portion_below, seriestype=:line, label="Portion Below")
plot!(pp, df.soc_plugged, df.soc_below, seriestype=:line, label="SOC Below")
plot!(pp, df.soc_plugged, df.soc_above, seriestype=:line, label="SOC Above")

pp
