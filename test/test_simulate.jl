timestep = 1. # 1 hour
include("utils.jl")
include("simulate.jl")

## Test simustep
@assert get_simustep_deterministic(DateTime("2024-07-15T07:34:56"))(4., 4., 1., 0.) == (4.0, 1.0, 0.0)
@assert get_simustep_deterministic(DateTime("2024-07-15T08:34:56"))(4., 4., 1., 0.) == (0.0, 1.0, 1.0)
@assert get_simustep_deterministic(DateTime("2024-07-15T09:34:56"))(0., 0., 1., 1.)[3] < 1.0
@assert get_simustep_deterministic(DateTime("2024-07-15T17:34:56"))(0., 0., 1., .9)[3] < 0.9
