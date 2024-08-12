using Roots, Distributions

"""
Convert a given number of steps into a period expressed in hours and minutes.

# Arguments
- `steps::Int`: The number of timestep increments to convert.

# Returns
- `Period`: A period of time represented in hours and minutes.
"""
function periodstep(steps::Int)
    Hour(round(Int, timestep * steps)) + Minute((timestep * steps - round(Int, timestep * steps)) * 60)
end

"""
For a SOC for plugged-in vehicles (`enerfrac_plugged`), determine the portion of the fleet that is below a
prescribed level (`enerfrac_needed`). Assumes the fleet follows a truncated normal distribution.

# Arguments
- `enerfrac_plugged::Float64`: The current fraction of the fleet that is plugged in.
- `enerfrac_needed::Float64`: The prescribed energy fraction required.

# Returns
- `Tuple{Float64, Float64, Float64}`:
  1. Portion of vehicles below the required energy fraction.
  2. Energy fraction for the low-power portion of the fleet.
  3. Energy fraction for the high-power portion of the fleet.
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

"""
Adjusts the fraction of vehicles and energy plugged based on the vehicles and energy considered below a certain level.

# Arguments
- `tup::Tuple{Float64, Float64, Float64}`: A tuple containing:
  1. The number of vehicles plugged-in.
  2. The state-of-charge for plugged-in vehicles.
  3. The state-of-charge for driving vehicles.
- `enerfrac_below::Float64`: The lower threshold for state-of-charge
- `vehicles_below::Float64`: The number of vehicles below the given threshold.

# Returns
- `Tuple{Float64, Float64, Float64}`:
  1. Adjusted number of vehicles plugged-in.
  2. Adjusted state-of-charge for plugged-in vehicles.
  3. Original state-of-charge for driving vehicles.
"""
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
