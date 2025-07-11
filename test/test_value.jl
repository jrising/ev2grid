using Test

# Constants for the tests
const weight_portion_below = 1.0
const soc_max = 1.0
const vehicle_capacity = 75.7
const vehicles = 4.
const efficiency = 0.95

# Unit tests for value_energy function
@testset "Test value_energy function" begin
    @test value_energy(1.0, 0.3, 0.3, vehicles) ≈ 0.95
    @test value_energy(0.0, 0.5, 0.3, vehicles) ≈ 0.7071067811865476 - weight_portion_below * 0.0
    @test value_energy(0.0, 0.5, 0.48, vehicles) ≈ (1.0 - 0.0) * sqrt( (0.5 - 0.48) / (soc_max - 0.48) ) - weight_portion_below * 0.0
    @test value_energy(0.0, 0.5, 0.3, vehicles) ≈ 0.7071067811865476 - weight_portion_below * 0.0
    @test value_energy(0.0, 0.6, 0.48, vehicles) ≈ (1.0 - 0.0) * sqrt( (0.6 - 0.48) / (soc_max - 0.48) ) - weight_portion_below * 0.0
    @test value_energy(0.0, 0.2, 0.3, vehicles) === -Inf  # Testing impossible state
end

# Unit tests for value_power_action function
@testset "Test value_power_action function" begin
    @test value_power_action(5.0, 0.1, 0, vehicles) ≈ -5.0 * (0.1 * vehicle_capacity * vehicles) / efficiency
    @test value_power_action(5.0, -0.1, 0, vehicles) ≈ 5.0 * (-0.1 * vehicle_capacity * vehicles)
end

# Unit tests for value_power_newstate function
@testset "Test value_power_newstate function" begin
    @test value_power_newstate(5.0, 0.5, 0.1, vehicles) ≈ -5.0 * (0.1 * vehicle_capacity * vehicles * 0.5) / efficiency
end

println("All tests passed!")
