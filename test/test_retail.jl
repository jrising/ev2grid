using Test

include("../src/retail.jl")

@testset "Peak hours and retail price tests" begin

    # Test is_peak function
    @test is_peak(DateTime("2024-07-29T13:00:00")) == true
    @test is_peak(DateTime("2024-07-29T09:00:00")) == false
    @test is_peak(DateTime("2024-01-01T15:00:00")) == false # Assuming it’s a holiday
    @test is_peak(DateTime("2024-01-02T15:00:00")) == true

    # Test get_retail_price function
    @test get_retail_price(DateTime("2024-07-29T13:00:00")) == 0.1473
    @test get_retail_price(DateTime("2024-07-29T09:00:00")) == 0.07242
    @test get_retail_price(DateTime("2024-01-01T02:00:00")) == 0.08356  # Assuming it's a holiday
    @test get_retail_price(DateTime("2024-01-02T15:00:00")) == 0.1720

end

# Unit test
@testset "get_demand_cost" begin
    dt0 = DateTime(2024, 8, 1, 0, 0) # Starting at midnight
    kwbytimestep = [1.0, 2.5, 3.0, 0.5, 4.0, 3.5, 2.0, 6.0] # Power usage in kW
    timestep = 1.0 / 4 # 15 minutes

    result = get_demand_cost(dt0, kwbytimestep, timestep)
    expected_result = 18.73 * (6 / 3) # The expected demand cost based on the mock data

    @test result ≈ expected_result
end
