## These functions assume the existence of global variables vehicles, EE, enerfrac_min, enerfrac_max, FF

"""
Translate a continuous value to a discrete state.

Arguments:
- xx: Continuous value to be discretized.
- xmin: Minimum value of the continuous space.
- xmax: Maximum value of the continuous space.
- num: Number of discrete states.

Returns:
- ii: Discrete state corresponding to the continuous value.
"""
discrete_float(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min((num - 1) * (xx - xmin) / (xmax - xmin) + 1, num))
discrete_round(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min(round(Int, (num - 1) * (xx - xmin) / (xmax - xmin)) + 1, num))
discrete_floatbelow(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min((num - 2) * (xx - xmin) / (xmax - xmin) + 2, num))
discrete_roundbelow(xx::Float64, xmin::Float64, xmax::Float64, num::Int) = max(1, min(round(Int, (num - 2) * (xx - xmin) / (xmax - xmin)) + 2, num))

function asstate_float(tup::Tuple{Float64, Float64, Float64})
    (discrete_float(tup[1], 0., vehicles, EE), discrete_floatbelow(tup[2], enerfrac_min, enerfrac_max, FF), discrete_floatbelow(tup[3], enerfrac_min, enerfrac_max, FF))
end

basestate(tup::Tuple{Float64, Float64, Float64}) = floor.(Int, tup)
ceilstate(tup::Tuple{Float64, Float64, Float64}, kk::Int) = ceil(Int, tup[kk])

function getstatedim(tup::Tuple{Float64, Float64, Float64}, kk::Int)
    tup[kk]
end

makeindex1(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2ceil, state2base[2], state2base[3])
makeindex2(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2base[1], state2ceil, state2base[3])
makeindex3(state2base::Tuple{Int64, Int64, Int64}, state2ceil::Int) = CartesianIndex(state2base[1], state2base[2], state2ceil)

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

function simustate2(simustep::Function, vehicles_plugged_range::Vector{Float64}, vehicle_split::Vector{Tuple{Float64, Float64, Float64}},
                    ff12_byaction::Array{Int64, 4}, enerfrac1_byaction::Array{Float64, 4}, enerfrac_range::Vector{Float64})
    # Note: We impose costs from enerfrac-below vehicles, but do not adjust state because it pushes up plugged-in enerfrac every period
    statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]), enerfrac1_byaction[pp, ee, ff1, ff2], enerfrac_range[ff2]), vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];

    state2 = asstate_float.(statevar2);
    state2base = basestate.(state2);
    state2ceil1 = ceilstate.(state2, 1);
    state2ceil2 = ceilstate.(state2, 2);
    state2ceil3 = ceilstate.(state2, 3);

    probbase1 = (state2ceil1 - getstatedim.(state2, 1));
    probbase2 = (state2ceil2 - getstatedim.(state2, 2));
    probbase3 = (state2ceil3 - getstatedim.(state2, 3));

    return state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3
end

function calcVV1byact(VV2::Array{Float64, 3}, state2base::Array{Tuple{Int64, Int64, Int64}, 4}, state2ceil1::Array{Int64, 4}, probbase1::Array{Float64, 4},
                      state2ceil2::Array{Int64, 4}, probbase2::Array{Float64, 4}, state2ceil3::Array{Int64, 4}, probbase3::Array{Float64, 4})
    ## Index into VV2 three times and combine
    VV1base = VV2[CartesianIndex.(state2base)];
    VV1ceil1 = VV2[makeindex1.(state2base, state2ceil1)]; # makeindex1 much faster than anonymous functions
    VV1ceil2 = VV2[makeindex2.(state2base, state2ceil2)];
    VV1ceil3 = VV2[makeindex3.(state2base, state2ceil3)];

    VV1byactthismc = ((probbase1 + probbase2 + probbase3) .* VV1base + (1 .- probbase1) .* VV1ceil1 + (1 .- probbase2) .* VV1ceil2 + (1 .- probbase3) .* VV1ceil3) / 3;
end
