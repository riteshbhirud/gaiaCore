using Test
using GaiaDB
using DataFrames
using Statistics

const BASE_URL = get(ENV, "GAIADB_URL", "http://localhost:3000")

@testset "GaiaDB.jl" begin

    @testset "list_locations" begin
        df = list_locations(BASE_URL)
        @test df isa DataFrame
        @test nrow(df) >= 20
        @test "location_id" in names(df)
        @test "city" in names(df)
        @test "latitude" in names(df)
        @test "longitude" in names(df)
    end

    @testset "list_locations with filters" begin
        df = list_locations(BASE_URL; city = "FRESNO")
        @test df isa DataFrame
        @test nrow(df) > 0
        @test all(df.city .== "FRESNO")

        df2 = list_locations(BASE_URL; limit = 3)
        @test nrow(df2) == 3
    end

    @testset "get_data_sources" begin
        df = get_data_sources(BASE_URL)
        @test df isa DataFrame
    end

    @testset "get_exposures" begin
        df = get_exposures(BASE_URL)
        @test df isa DataFrame
    end

    @testset "initialize_agents_from_sdoh" begin
        agents = initialize_agents_from_sdoh(BASE_URL, 500)
        @test length(agents) == 500
        @test all(a -> a isa EpiAgent, agents)

        # Every agent has valid fields
        @test all(a -> 1 <= a.income_quantile <= 5, agents)
        @test all(a -> 0.0 <= a.transit_dependency <= 1.0, agents)
        @test all(a -> a.home_location_id != a.work_location_id, agents)

        # Transit dependency is NOT uniform — mean should be near the
        # pop-density-weighted tract mean (~17%), not 50%
        td = [a.transit_dependency for a in agents]
        @test mean(td) < 0.35  # well below 0.5 (uniform midpoint)
        @test mean(td) > 0.05  # but not near zero

        # Income quintiles are populated across all 5
        for q in 1:5
            @test count(a -> a.income_quantile == q, agents) > 0
        end
    end

end
