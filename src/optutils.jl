## These functions assume the existence of global variables vehicles, EE, soc_min, soc_max, FF

"""
Translate a continuous value to a discrete state.

Arguments:
- xx: Continuous value to be discretized.
- xmin: Minimum value of the continuous space.
- xmax: Maximum value of the continuous space.
- num: Number of discrete states.

Returns:
- ii: Discrete state corresponding to the continuous value, ensuring it's within the allowable range.
"""
discrete_float(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min((num - 1) * (xx - xmin) / (xmax - xmin) + 1, num))

"""
Translate a continuous value to the nearest discrete state using rounding.

Arguments:
- xx: Continuous value to be discretized.
- xmin: Minimum value of the continuous space.
- xmax: Maximum value of the continuous space.
- num: Number of discrete states.

Returns:
- ii: Discrete state rounded nearest to the continuous value, constrained to within bounds.
"""
discrete_round(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min(round(Int, (num - 1) * (xx - xmin) / (xmax - xmin)) + 1, num))

"""
Translate a continuous value to a discrete state, biased towards a lower state.

Arguments:
- xx: Continuous value to be discretized.
- xmin: Minimum value of the continuous space.
- xmax: Maximum value of the continuous space.
- num: Number of discrete states.

Returns:
- ii: Discrete state corresponding to the continuous value, with an additional "below xmin" state.
"""
discrete_floatbelow(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min((num - 2) * (xx - xmin) / (xmax - xmin) + 2, num))

"""
Translate a continuous value to the nearest lower discrete state using rounding.

Arguments:
- xx: Continuous value to be discretized.
- xmin: Minimum value of the continuous space.
- xmax: Maximum value of the continuous space.
- num: Number of discrete states.

Returns:
- ii: Discrete state rounded nearest to the continuous value, with an additional "below xmin" state.
"""
discrete_roundbelow(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min(round(Int, (num - 2) * (xx - xmin) / (xmax - xmin)) + 2, num))

"""
Convert a tuple of continuous values into discrete states. This
assumes the first state is number of vehicles, the second is the
state-of-charge for plugged-in vehicles, and the third is the
state-of-charge for driving vehicles.

Arguments:
- tup: Tuple of three Float64 values representing continuous states.

Returns:
- Tuple of discrete states for each element in the input tuple.
"""
function asstate_float(tup::Tuple{Float64, Float64, Float64})
    (discrete_float(tup[1], 0., vehicles, EE), discrete_floatbelow(tup[2], soc_min, soc_max, FF), discrete_floatbelow(tup[3], soc_min, soc_max, FF))
end

"""
Floor the values in a tuple to convert them to integers.

Arguments:
- tup: Tuple of three Float64 values.

Returns:
- Tuple of floored integers.
"""
basestate(tup::Tuple{Float64, Float64, Float64}) = floor.(Int, tup)

"""
Ceil one element of a tuple to turn it into an integer.

Arguments:
- tup: Tuple of three Float64 values.
- kk: Index of the element to ceil.

Returns:
- Ceiled integer of the specified element.
"""
ceilstate(tup::Tuple{Float64, Float64, Float64}, kk::Int) = ceil(Int, tup[kk])

"""
Retrieve the dimension of a tuple.

Arguments:
- tup: Tuple of three Float64 values.
- kk: Index of the dimension to retrieve.

Returns:
- The value at the specified index in the tuple.
"""
function getstatedim(tup::Tuple{Float64, Float64, Float64}, kk::Int)
    tup[kk]
end

"""
Creates a Cartesian index from a tuple of base states and a ceiling state for the first dimension.

Arguments:
- state2base: Tuple of three Int64 values representing base states.
- state2ceil: Ceil value of the first dimension.

Returns:
- CartesianIndex for the given state.
"""
makeindex1(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2ceil, state2base[2], state2base[3])

"""
Creates a Cartesian index from a tuple of base states and a ceiling state for the second dimension.

Arguments:
- state2base: Tuple of three Int64 values representing base states.
- state2ceil: Ceil value of the second dimension.

Returns:
- CartesianIndex for the given state.
"""
makeindex2(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2base[1], state2ceil, state2base[3])

"""
Creates a Cartesian index from a tuple of base states and a ceiling state for the third dimension.

Arguments:
- state2base: Tuple of three Int64 values representing base states.
- state2ceil: Ceil value of the third dimension.

Returns:
- CartesianIndex for the given state.
"""
makeindex3(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2base[1], state2base[2], state2ceil)

"""
Break a state down into its base and ceiling components with probabilities.

Arguments:
- tup: Tuple of three Float64 values representing a state.

Returns:
- Seven-element tuple consisting of base states, ceiling states, and transition probabilities.
"""

function breakstate(tup::Tuple{Float64, Float64, Float64})
    state = asstate_float(tup);
    statebase = basestate(state);
    stateceil1 = ceilstate(state, 1);
    stateceil2 = ceilstate(state, 2);
    stateceil3 = ceilstate(state, 3);

    probbase1 = (stateceil1 - getstatedim(state, 1));
    probbase2 = (stateceil2 - getstatedim(state, 2));
    probbase3 = (stateceil3 - getstatedim(state, 3));

    return statebase, stateceil1, probbase1, stateceil2, probbase2, stateceil3, probbase3
end

"""
Break multiple states down into base and ceiling components with probabilities.

Arguments:
- tuparr: Array of three-element Float64 tuples representing states.

Returns:
- Seven-element tuple containing arrays of base states, ceiling states, and transition probabilities.
"""
function breakstate(tuparr::Array{Tuple{Float64, Float64, Float64}, 4})
    state2 = asstate_float.(tuparr);
    state2base = basestate.(state2);
    state2ceil1 = ceilstate.(state2, 1);
    state2ceil2 = ceilstate.(state2, 2);
    state2ceil3 = ceilstate.(state2, 3);

    probbase1 = (state2ceil1 - getstatedim.(state2, 1));
    probbase2 = (state2ceil2 - getstatedim.(state2, 2));
    probbase3 = (state2ceil3 - getstatedim.(state2, 3));

    return state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3
end

"""
Combine multiple action values based on probability and state indices.

Arguments:
- VV2: 3D Array of Float64 values representing action values.
- state2base: Array of three-element Int64 tuples representing base states.
- state2ceil1, state2ceil2, state2ceil3: Arrays of Int64 values for ceiling states corresponding to each dimension.
- probbase1, probbase2, probbase3: Arrays of Float64 values representing probabilities for each dimension.

Returns:
- 3D array of combined action values considering transitions and probabilities.
"""
function combinebyact(VV2::Array{Float64, 3}, state2base::Array{Tuple{Int64, Int64, Int64}, 4},
                      state2ceil1::Array{Int64, 4}, probbase1::Array{Float64, 4},
                      state2ceil2::Array{Int64, 4}, probbase2::Array{Float64, 4},
                      state2ceil3::Array{Int64, 4}, probbase3::Array{Float64, 4})
    ## Index into VV2 three times and combine
    VV1base = VV2[CartesianIndex.(state2base)];
    VV1ceil1 = VV2[makeindex1.(state2base, state2ceil1)]; # makeindex1 much faster than anonymous functions
    VV1ceil2 = VV2[makeindex2.(state2base, state2ceil2)];
    VV1ceil3 = VV2[makeindex3.(state2base, state2ceil3)];

    VV1byactthismc = ((probbase1 + probbase2 + probbase3) .* VV1base + (1 .- probbase1) .* VV1ceil1 + (1 .- probbase2) .* VV1ceil2 + (1 .- probbase3) .* VV1ceil3) / 3;
end
