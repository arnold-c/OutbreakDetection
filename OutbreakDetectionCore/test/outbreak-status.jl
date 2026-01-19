@testset "outbreak-status.jl" begin
    using OutbreakDetectionCore
    using StatsBase
    using StaticArrays
    using StructArrays

    # Helper function to create mock SEIRRun
    function create_mock_seir_run(incidence::Vector{Int64})
        n = length(incidence)
        states = [SVector{5, Int64}(0, 0, 0, 0, 0) for _ in 1:n]
        reff = zeros(Float64, n)
        return OutbreakDetectionCore.SEIRRun(
            states = states, incidence = incidence, Reff = reff
        )
    end

    @testset "Single simulation - single outbreak" begin
        # Create incidence with one clear outbreak
        inc = [
            repeat([1], 9)...,      # Days 1-9: baseline
            repeat([12], 51)...,    # Days 10-60: outbreak (51 days, 612 cases)
            repeat([1], 39)...,     # Days 61-99: baseline
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        # Test return types and lengths
        @test outbreak_status_vecs isa Vector{Vector{Bool}}
        @test outbreak_bounds_vecs isa Vector{Matrix{Int64}}
        @test length(outbreak_status_vecs) == 1
        @test length(outbreak_bounds_vecs) == 1

        # Test outbreak status vector
        status = outbreak_status_vecs[1]
        @test length(status) == 99
        @test all(status[1:9] .== false)      # Before outbreak
        @test all(status[10:60] .== true)    # During outbreak
        @test all(status[61:99] .== false)    # After outbreak

        # Test outbreak bounds
        bounds = outbreak_bounds_vecs[1]
        @test size(bounds, 1) == 1        # One outbreak
        @test size(bounds, 2) == 4        # [start, end, duration, size]
        @test bounds[1, 1] == 10          # Start day
        @test bounds[1, 2] == 60          # End day
        @test bounds[1, 3] == 51          # Duration
        @test bounds[1, 4] == 612         # Total cases (12 * 51)
    end

    @testset "Single simulation - multiple outbreaks" begin
        # Create incidence with two outbreaks
        inc = [
            repeat([1], 9)...,      # Days 1-9: baseline
            repeat([12], 51)...,    # Days 10-60: outbreak 1
            repeat([1], 19)...,     # Days 61-79: baseline
            repeat([8], 81)...,     # Days 80-160: outbreak 2
            repeat([1], 39)...,     # Days 161-199: baseline
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        # Test two outbreaks detected
        @test size(bounds, 1) == 2

        # Test first outbreak
        @test all(status[10:60] .== true)
        @test bounds[1, 1] == 10
        @test bounds[1, 2] == 60

        # Test second outbreak
        @test all(status[80:160] .== true)
        @test bounds[2, 1] == 80
        @test bounds[2, 2] == 160

        # Test baseline periods
        @test all(status[1:9] .== false)
        @test all(status[61:79] .== false)
        @test all(status[161:199] .== false)
    end

    @testset "Single simulation - no outbreaks" begin
        # All incidence below threshold
        inc = repeat([2], 100)

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        # No outbreaks detected
        @test all(status .== false)
        @test size(bounds, 1) == 0  # Empty matrix
    end

    @testset "Single simulation - period above threshold but too short" begin
        # Period above threshold but doesn't meet minimum duration
        inc = [
            repeat([1], 9)...,      # Days 1-9: baseline
            repeat([15], 20)...,    # Days 10-29: above threshold but only 20 days (< 30)
            repeat([1], 70)...,     # Days 30-99: baseline
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        # Period not classified as outbreak (too short)
        @test all(status .== false)
        @test size(bounds, 1) == 0
    end

    @testset "Single simulation - period above threshold but too small" begin
        # Period meets duration but not size requirement
        inc = [
            repeat([1], 9)...,      # Days 1-9: baseline
            repeat([8], 40)...,     # Days 10-49: 40 days, 320 cases (< 500)
            repeat([1], 50)...,     # Days 50-99: baseline
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        # Period not classified as outbreak (too small)
        @test all(status .== false)
        @test size(bounds, 1) == 0
    end

    @testset "Multiple simulations" begin
        # Create three different simulations
        inc1 = [
            repeat([1], 9)...,
            repeat([12], 51)...,    # Outbreak
            repeat([1], 39)...,
        ]

        inc2 = [
            repeat([2], 19)...,
            repeat([20], 41)...,    # Outbreak
            repeat([2], 39)...,
        ]

        inc3 = repeat([1], 99)      # No outbreak

        seir_results = StructVector(
            [
                create_mock_seir_run(inc1),
                create_mock_seir_run(inc2),
                create_mock_seir_run(inc3),
            ]
        )

        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        # Test correct number of simulations
        @test length(outbreak_status_vecs) == 3
        @test length(outbreak_bounds_vecs) == 3

        # Test first simulation
        @test all(outbreak_status_vecs[1][10:60] .== true)
        @test size(outbreak_bounds_vecs[1], 1) == 1

        # Test second simulation
        @test all(outbreak_status_vecs[2][20:60] .== true)
        @test size(outbreak_bounds_vecs[2], 1) == 1

        # Test third simulation (no outbreak)
        @test all(outbreak_status_vecs[3] .== false)
        @test size(outbreak_bounds_vecs[3], 1) == 0
    end

    @testset "Edge case - outbreak at start" begin
        # Outbreak starts at day 1
        inc = [
            repeat([15], 40)...,    # Days 1-40: outbreak
            repeat([1], 59)...,     # Days 41-99: baseline
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        @test all(status[1:40] .== true)
        @test bounds[1, 1] == 1
        @test bounds[1, 2] == 40
    end

    @testset "Edge case - outbreak at end" begin
        # Outbreak ends at last day
        inc = [
            repeat([1], 59)...,     # Days 1-59: baseline
            repeat([15], 40)...,    # Days 60-99: outbreak
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        @test all(status[60:99] .== true)
        @test bounds[1, 1] == 60
        @test bounds[1, 2] == 99
    end

    @testset "Edge case - entire period is outbreak" begin
        # All days are outbreak
        inc = repeat([15], 99)

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        @test all(status .== true)
        @test size(bounds, 1) == 1
        @test bounds[1, 1] == 1
        @test bounds[1, 2] == 99
    end

    @testset "Different outbreak specifications" begin
        inc = [
            repeat([1], 9)...,
            repeat([12], 51)...,
            repeat([1], 39)...,
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])

        # Test with stricter threshold
        outbreak_spec_strict = OutbreakSpecification(15, 30, 500)
        status_strict, bounds_strict = create_outbreak_status_vecs(
            seir_results, outbreak_spec_strict
        )
        @test all(status_strict[1] .== false)  # 12 < 15, no outbreak

        # Test with looser duration requirement
        outbreak_spec_loose_dur = OutbreakSpecification(5, 10, 500)
        status_loose, bounds_loose = create_outbreak_status_vecs(
            seir_results, outbreak_spec_loose_dur
        )
        @test all(status_loose[1][10:60] .== true)  # Still outbreak

        # Test with looser size requirement
        outbreak_spec_loose_size = OutbreakSpecification(5, 30, 100)
        status_loose_size, bounds_loose_size = create_outbreak_status_vecs(
            seir_results, outbreak_spec_loose_size
        )
        @test all(status_loose_size[1][10:60] .== true)  # Still outbreak
    end

    @testset "Outbreak bounds matrix structure" begin
        inc = [
            repeat([1], 9)...,
            repeat([12], 51)...,
            repeat([1], 19)...,
            repeat([8], 81)...,
            repeat([1], 39)...,
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        _, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        bounds = outbreak_bounds_vecs[1]

        # Verify matrix has correct structure
        @test size(bounds, 2) == 4  # [start, end, duration, size]

        # Verify duration calculation
        for i in 1:size(bounds, 1)
            @test bounds[i, 3] == bounds[i, 2] - bounds[i, 1] + 1
        end

        # Verify bounds are in chronological order
        if size(bounds, 1) > 1
            for i in 1:(size(bounds, 1) - 1)
                @test bounds[i, 2] < bounds[i + 1, 1]  # No overlap
            end
        end
    end

    @testset "Consistency between status and bounds" begin
        inc = [
            repeat([1], 9)...,
            repeat([12], 51)...,
            repeat([1], 19)...,
            repeat([8], 81)...,
            repeat([1], 39)...,
        ]

        seir_results = StructVector([create_mock_seir_run(inc)])
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        status = outbreak_status_vecs[1]
        bounds = outbreak_bounds_vecs[1]

        # Verify that status vector matches bounds
        for i in 1:size(bounds, 1)
            start_idx = bounds[i, 1]
            end_idx = bounds[i, 2]
            @test all(status[start_idx:end_idx] .== true)
        end

        # Verify that non-outbreak periods are false
        outbreak_days = Set{Int}()
        for i in 1:size(bounds, 1)
            for day in bounds[i, 1]:bounds[i, 2]
                push!(outbreak_days, day)
            end
        end

        for day in 1:length(status)
            if day ∉ outbreak_days
                @test status[day] == false
            end
        end
    end
end
