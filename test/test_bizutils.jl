using Test

timestep = 1.
include("../src/bizutils.jl")

# Test for periodstep function
@testset "Testing periodstep" begin
    # Test with 0 steps
    @test periodstep(0) == Hour(0) + Minute(0)

    # Test with 1 step (assuming timestep is defined as 1 hour)
    global timestep = 1.0   # Example value, adjust as necessary for your use case
    @test periodstep(1) == Hour(1) + Minute(0)

    # Test with fractional steps
    global timestep = 0.5   # 30 minutes
    @test periodstep(1) == Hour(0) + Minute(30)
    @test periodstep(2) == Hour(1) + Minute(0)
end

# Test for split_below function
@testset "Testing split_below" begin
    # Edge cases
    @test split_below(0.0, 0.3) == (1.0, 0.0, 0.3)
    @test split_below(1.0, 0.3) == (0.0, 0.3, 1.0)

    # General case
    portion_below, low_mean, high_mean = split_below(0.5, 0.3)
    @test portion_below < 1e-3
    @test low_mean ≈ 0.29 atol=0.01
    @test high_mean ≈ 0.5 atol=0.01

    # Additional test with different values
    portion_below, low_mean, high_mean = split_below(0.3, 0.3)
    @test portion_below ≈ 0.5
    @test low_mean ≈ 0.26 atol=0.01
    @test high_mean ≈ 0.34 atol=0.01
end

# Test for adjust_below function
@testset "Testing adjust_below" begin
    global vehicles = 4.

    # Initial condition
    tup = (2.0, 0.3, 0.5)
    vehicles_below = 2.0

    # Test with valid input
    vehicles_plugged, soc_plugged, _ = adjust_below(tup, 0.3, vehicles_below)
    @test vehicles_plugged ≈ 4.0
    @test soc_plugged ≈ 0.3

    # Test when no vehicles are below
    tup = (2.0, 0.7, 0.5)
    vehicles_below = 0.0
    vehicles_plugged, soc_plugged, _ = adjust_below(tup, 0.3, vehicles_below)
    @test vehicles_plugged ≈ 2.0       # Remains the same
    @test soc_plugged ≈ 0.7       # Remains the same
end

