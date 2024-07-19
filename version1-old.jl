"""
Return matrix of changes in energy, across actions.
"""
function make_action_Fplugged(E_avail_max::Float64, FF::Int, P_min::Float64, P_max::Float64, PP::Int)
    power = range(P_min, P_max, length=PP)
    frac = range(0, 1, length=FF)
    energy = frac * E_avail_max

    denergy_matrix = -timestep * repeat(power, 1, FF)
    energy0_matrix = repeat(energy', PP, 1)
    energy1_matrix = max.(0., min.(E_avail_max, energy0_matrix + denergy_matrix))

    [zeros(FF)'; energy1_matrix - energy0_matrix]
end

## E_avail_max = E_max_expected[12] - E_min_expected[12]
## make_action_Fplugged(E_avail_max, FF, P_min, P_max, PP)

"""
Constructs a mapping between before-states and after-states based on given parameters.

Arguments:
- EE: Number of discretized states of the energy.
- P_min: Greatest discharge of energy.
- P_max: Greatest charing of energy.
- PP: Number of possible changes in energy.

Returns:
- state1: State before actions take effect.
- state2a: Possible resulting states (first option).
- state2b: Possible resulting states (second option).
- xprob2: Probability matrix of transitions to state2b.
"""
function make_actions(E_avail_max::Float64, FF::Int, P_min::Float64, P_max::Float64, PP::Int)
    # Initialize matrices to store states and probabilities
    state1 = repeat((1:FF)', 1+PP, 1) # as indices

    state2a = collect(1:FF)' # as indices
    state2b = collect(1:FF)' # as indices
    xprob2 = repeat([1.0], 1, FF) # probability of going to state2b

    frac = range(0, 1, length=FF)
    energy = frac * E_avail_max
    actions_energy = make_action_Fplugged(E_avail_max, FF, P_min, P_max, PP)

    for ii in 2:size(actions_energy)[1]
        newfrac = (energy + actions_energy[ii, :]) / E_avail_max
        state_continuous = discrete.(newfrac, 0., 1., FF, false)
        state2a = [state2a; ceil.(state_continuous)']
        state2b = [state2b; floor.(state_continuous)']
        xprob2 = [xprob2; 1.0 .- (state_continuous - floor.(state_continuous))']
    end

    return state1, state2a, state2b, xprob2
end
## state1, state2a, state2b, xprob2 = make_actions(E_avail_max, FF, P_min, P_max, PP)


function optimize(dt0::DateTime, SS::Int)
    strat = zeros(Int64, SS-1, EE)

    # STEP 1: Calculate V[S] under every scenario
    VV2 = repeat(value.(dt0 + periodstep(SS), range(enerfrac_min, enerfrac_max, length=EE)), 1, FF, FF)

    # STEP 2: Determine optimal action for t = S-1 and back
        [make_action_Fplugged(enerfrac_plugged) for enerfrac_plugged=range(0., 1., FF), enerfrac_driving=range(0., 1., FF), vehicles_plugged=range(0., vehicles, EE)]

    for tt in (SS-1):-1:1
        println(tt)
        VV1 = zeros(C, F)

        for Fi in 1:F
            sums1 = zeros(size(state2a))
            sums2 = zeros(size(state2b))

            for ii in 1:N
                P_f2 = simustep(P_f[Fi], process, fuel_sigma)
                Fi2 = discrete(P_f2, P_f_min, P_f_max, F)

                sums1 += VV2[(Fi2-1)*C .+ state2a]
                sums2 += VV2[(Fi2-1)*C .+ state2b]
            end

            ecost_later1 = reshape(sums1, size(state2a)) / N
            ecost_later2 = reshape(sums2, size(state2b)) / N
            later = (xprob2 .* ecost_later2 .+ (1 .- xprob2) .* ecost_later1)
            ecost_now = ecost(S_C[state1], Q_const, q_c, q_f, P_f[Fi], alpha_wind, alpha_solar, declining, e_of_scale)
            values = ecost_now + exp(-discount) * later

            indices = argmin(values, dims=1)
            strat[tt, :, Fi] .= [ind[1].I[1] for ind in eachcol(indices)]
            VV1[:, Fi] .= vec(values[indices])
        end

        VV2 = copy(VV1)
    end

    return strat
end

"""
Return matrix of changes in energy, across actions.
"""
function make_action_Fplugged(enerfrac0::Float64)
    fracpower = range(fracpower_min, fracpower_max, length=PP)

    denerfracs = timestep * fracpower
    enerfrac1s = max.(0., min.(1., enerfrac0 .+ denerfracs))

    [0; enerfrac1s .- enerfrac0]
end

## make_action_Fplugged(0.5)
