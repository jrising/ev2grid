using Test

@testset verbose=true "Optimizer tests" begin
    include("test_retail.jl")
    include("test_bizutils.jl")
    include("test_simulate.jl")
end

