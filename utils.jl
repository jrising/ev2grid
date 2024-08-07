using Roots, Distributions

timestep = 1. # 1 hour

function periodstep(steps::Int)
    Hour(round(Int, timestep * steps)) + Minute((timestep * steps - round(Int, timestep * steps)) * 60)
end

"""
For a given enerfrac_plugged, determine the portion of the plugged-in fleet below the prescribed level
Returns:
 - portion of vehicles below the level
 - enerfrac for low-power portion
 - enerfrac for high-power portion
"""
function split_below(enerfrac_plugged::Float64, enerfrac_needed::Float64)
    # Limits come from truncated(Normal(mu, 0.05), lower=0., upper=1.)
    if enerfrac_plugged ≤ 0.03989422804014327
        return 1.0, 0.0, enerfrac_needed
    elseif enerfrac_plugged ≥ 0.9601057719598567
        return 0.0, enerfrac_needed, 1.0
    end
    ## Use normal distribution, so I can calculate means of truncated
    # Need to get a truncated normal with the desired mean
    f(mu) = mean(truncated(Normal(mu, 0.05), lower=0., upper=1.)) - enerfrac_plugged
    mu = find_zero(f, (0., 1.))
    dist = truncated(Normal(mu, 0.05), lower=0., upper=1.)
    portion_below = cdf(dist, enerfrac_needed)
    portion_below, mean(truncated(dist, upper=enerfrac_needed)), mean(truncated(dist, lower=enerfrac_needed))
end

split_below(0.5, 0.5)
split_below(0.5, 0.3)
split_below(0.7, 0.3)
split_below(0.3, 0.3)
split_below(0.1, 0.3)

function adjust_below(tup::Tuple{Float64, Float64, Float64}, enerfrac_below::Float64, vehicles_below::Float64)
    vehicles_plugged = tup[1] + vehicles_below
    @assert vehicles_plugged ≤ vehicles + 1e-8

    if vehicles_plugged > 0
        enerfrac_plugged = (enerfrac_below * vehicles_below + tup[2] * tup[1]) / vehicles_plugged
    else
        enerfrac_plugged = tup[2]
    end
    (vehicles_plugged, enerfrac_plugged, tup[3])
end
