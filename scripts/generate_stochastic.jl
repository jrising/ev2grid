using Dates, DataFrames, SparseArrays

include("../src/customer.jl")
include("../src/bizutils.jl")
include("../src/simulate.jl")
include("../src/optutils.jl")

timestep = 1. # 1 hour
SS = 36 # project for 1.5 days
EE = 5 # 0 - 4 cars
FF = 11 # For both soc_plugged and soc_driving, 0 - 1
PP = 8 # discretized power choices (excluding no-change action)

dt0 = DateTime("2023-07-15T12:00:00")

soc_range = [0.; range(soc_min, soc_max, FF-1)];
vehicles_plugged_range = collect(range(0., vehicles, EE));
soc0_byaction = repeat(reshape(soc_range, 1, 1, FF, 1), PP, EE, 1, FF);

dsoc_FF = [make_actions(soc_plugged) for soc_plugged=soc_range];
dsoc = [dsoc_FF[ff][pp] for pp=1:PP, vehicles_plugged=vehicles_plugged_range, ff=1:FF, soc_driving=soc_range];

rowarr = zeros(PP, EE, FF, FF);
colarr = zeros(EE, FF, FF);
rowcar2lin = LinearIndices(rowarr);
colcar2lin = LinearIndices(colarr);

for tt in (SS-1):-1:1
    transprob = spzeros(PP * EE * FF * FF, EE * FF * FF);

    soc1_byaction = soc0_byaction .+ dsoc;

    dt1 = dt0 + periodstep(tt)
    soc_needed = soc_scheduled(dt1)
    vehicle_split = split_below.(soc_range, soc_needed)
    ff12_byaction = discrete_roundbelow.(soc1_byaction, soc_min, soc_max, FF);

    for mc in 1:1000
        simustep = get_simustep_stochastic(dt1)

        statevar2 = [adjust_below(simustep(vehicles_plugged_range[ee], vehicles_plugged_range[ee] * (1. - vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]), soc1_byaction[pp, ee, ff1, ff2], soc_range[ff2]), vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][2], vehicles_plugged_range[ee] * vehicle_split[ff12_byaction[pp, ee, ff1, ff2]][1]) for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF];
        state2base, state2ceil1, probbase1, state2ceil2, probbase2, state2ceil3, probbase3 = breakstate(statevar2)

        for pp=1:PP, ee=1:EE, ff1=1:FF, ff2=1:FF
            rowindex = rowcar2lin[pp, ee, ff1, ff2]
            state2baseii = state2base[pp, ee, ff1, ff2]
            colindex = colcar2lin[state2baseii...]

            transprob[rowindex, colindex] += (probbase1[pp, ee, ff1, ff2] + probbase2[pp, ee, ff1, ff2] + probbase3[pp, ee, ff1, ff2]) / 3
            colindex1b = colcar2lin[makeindex1(state2baseii, state2ceil1[pp, ee, ff1, ff2])]
            transprob[rowindex, colindex1b] += (1. - probbase1[pp, ee, ff1, ff2]) / 3
            colindex2b = colcar2lin[makeindex2(state2baseii, state2ceil2[pp, ee, ff1, ff2])]
            transprob[rowindex, colindex2b] += (1. - probbase2[pp, ee, ff1, ff2]) / 3
            colindex3b = colcar2lin[makeindex3(state2baseii, state2ceil3[pp, ee, ff1, ff2])]
            transprob[rowindex, colindex3b] += (1. - probbase3[pp, ee, ff1, ff2]) / 3
        end
    end

    transprob /= 1000.

    ## TODO: Store matrix here
end
