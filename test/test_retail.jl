using Test

include("src/retail.jl")

function test_is_peak_and_get_retail_price()
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

