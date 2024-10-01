using DataFrames
using Plots

include("../src/retail.jl")

df = DataFrame(hourstart=DateTime("2024-01-01T00:00:00"):Hour(1):DateTime("2024-12-31T00:00:00"))
df[!, :price] = get_retail_price.(df.hourstart)

pp = plot(df.hourstart, df.price, seriestype=:steppost, label="")
plot!(pp, size=(1000, 400))
savefig("retailprice.png")
