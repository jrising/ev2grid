using Test

include("utils.jl")
include("src/customer.jl")

# Simple unit tests for the enerfrac_scheduled function
@testset "Customer configuration Tests" begin

    # Test case 1: Before 9 AM
    dt_early = DateTime(2024, 7, 31, 7, 0, 0)
    @test enerfrac_scheduled(dt_early) == enerfrac_min

    # Test case 2: 9 AM window
    dt_at_9am = DateTime(2024, 7, 31, 8, 0, 0)
    @test enerfrac_scheduled(dt_at_9am) == enerfrac_9am

    # Test case 3: Just after 9 AM
    dt_late = DateTime(2024, 7, 31, 9, 1, 0)
    @test enerfrac_scheduled(dt_late) == enerfrac_min

    # Test case 4: Within the scheduled period (yesterday)
    dt_scheduled = DateTime(2024, 7, 30, 8, 0, 0)
    @test enerfrac_scheduled(dt_scheduled) == enerfrac_9am

end