"""
Simulate without a strategy for SS time steps.
"""
function simu_inactive(dt0::DateTime, SS::Int, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    # datetime, enerfrac_needed, vehicles_plugged_1, portion_below, enerfrac_toadd_below, enerfrac_avail_1, enerfrac_driving_1
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64}[]

    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

    push!(rows, (dt0 + periodstep(1), enerfrac_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], enerfrac_plugged_1, enerfrac_driving_1))

    for tt in 1:(SS-1)
        println(tt)
        enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

        simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], enerfrac_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt + 1), enerfrac_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], enerfrac_plugged_2, enerfrac_driving_2))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :enerfrac_needed, :vehicles_plugged, :portion_below, :enerfrac_below, :enerfrac_plugged, :enerfrac_driving])
end

vehicles_plugged_1 = 4.
enerfrac_plugged_1 = 0.5
enerfrac_driving_1 = 0.5
df = simu_inactive(dt0, SS, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)

pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label="")

"""
Plot a strategy over time.
"""
function simu_strat(dt0::DateTime, strat::AbstractArray{Int}, vehicles_plugged_1::Float64, enerfrac_plugged_1::Float64, enerfrac_driving_1::Float64)
    rows = Tuple{DateTime, Float64, Float64, Float64, Float64, Float64, Float64, Union{Missing, Float64}, Union{Missing, Tuple{Int, Int, Int}}}[]

    enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(1))
    vehicle_split = split_below(enerfrac_plugged_1, enerfrac_needed)

    push!(rows, (dt0, enerfrac_needed, vehicles_plugged_1, vehicle_split[1], vehicle_split[2], enerfrac_plugged_1, enerfrac_driving_1, missing, missing))

    for tt in 1:(size(strat)[1]-1)
        denerfrac = make_actions(enerfrac_plugged_1)

        index = asstate_round((vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1))
        pp = strat[tt, index...]

        enerfrac_plugged_2 = enerfrac_plugged_1 + denerfrac[pp]

        enerfrac_needed = enerfrac_scheduled(dt0 + periodstep(tt))
        vehicle_split = split_below(enerfrac_plugged_2, enerfrac_needed)

        if mcdraws == 1
            simustep = get_simustep_deterministic(dt0 + periodstep(tt))
        else
            simustep = get_simustep_stochastic(dt0 + periodstep(tt))
        end
        vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2 = adjust_below(simustep(vehicles_plugged_1, vehicles_plugged_1 * (1. - vehicle_split[1]), vehicle_split[3], enerfrac_driving_1), vehicle_split[2], vehicles_plugged_1 * vehicle_split[1])
        push!(rows, (dt0 + periodstep(tt), enerfrac_needed, vehicles_plugged_2, vehicle_split[1], vehicle_split[2], enerfrac_plugged_2, enerfrac_driving_2, denerfrac[pp], index))
        vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1 = vehicles_plugged_2, enerfrac_plugged_2, enerfrac_driving_2
    end

    df = DataFrame(rows)
    rename!(df, [:datetime, :enerfrac_needed, :vehicles_plugged, :portion_below, :enerfrac_toadd_below, :enerfrac_plugged, :enerfrac_driving, :denerfrac, :state])
end

vehicles_plugged_1 = 4.

pp = nothing
for enerfrac_plugged_1 in range(enerfrac_min, enerfrac_max, FF-1)
    local enerfrac_driving_1 = enerfrac_plugged_1
    local df = simu_strat(dt0, strat, vehicles_plugged_1, enerfrac_plugged_1, enerfrac_driving_1)
    global pp

    if pp == nothing
        pp = plot(df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1, legend=false)
    else
        plot!(pp, df.datetime, (df.enerfrac_plugged .* df.vehicles_plugged + df.enerfrac_driving .* (vehicles .- df.vehicles_plugged)) / vehicles, seriestype=:line, label=enerfrac_plugged_1, legend=false)
    end
end

pp
