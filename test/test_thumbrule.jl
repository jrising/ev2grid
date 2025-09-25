using Test

include("../src/bizutils.jl")
include("../src/customer.jl")
include("../src/simulate.jl")
include("../src/retail.jl")
include("../src/config.jl")
include("../src/value.jl")
include("../src/optutils.jl")
include("../src/fullsim.jl")
include("../src/thumbrule.jl")

dt0 = DateTime("2023-07-17T00:00:00")
drive_starts_time = Time(9, 0, 0)  # Example drive start time
park_starts_time = Time(17, 0, 0)  # Example park start time

dsoc_func, regrange_func = thumbrule_regrange(dt0, drive_starts_time, park_starts_time, drive_time_charge_level)

[regrange_func(tt) for tt in 1:36]
