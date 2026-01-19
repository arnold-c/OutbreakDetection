@testset "outbreak-classification.jl - OutbreakSpecification method" begin
    using OutbreakDetectionCore
    using StatsBase

    # Test classify_all_outbreaks! with OutbreakSpecification
    @testset "classify_all_outbreaks! with OutbreakSpecification" begin
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        inc_vec = [
            repeat([1], 9)...,
            repeat([12], 51)..., # Outbreak: 51 days, 612 cases
            repeat([1], 19)...,
            repeat([15], 11)..., # NOT outbreak: 11 days < 30
            repeat([1], 9)...,
            repeat([8], 81)..., # Outbreak: 81 days, 648 cases
            repeat([1], 19)...,
        ]

        abovethreshold_vec = vec(inc_vec .>= outbreak_spec.outbreak_threshold)
        abovethresholdrle = StatsBase.rle(abovethreshold_vec)
        all_outbreak_bounds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )

        outbreak_status = zeros(Bool, length(inc_vec))
        classify_all_outbreaks!(
            outbreak_status, all_outbreak_bounds, inc_vec, outbreak_spec
        )

        # Check that first period (10-60) is classified as outbreak
        @test all(outbreak_status[10:60] .== true)

        # Check that second period (80-90) is NOT classified as outbreak (too short)
        @test all(outbreak_status[80:90] .== false)

        # Check that third period (100-180) is classified as outbreak
        @test all(outbreak_status[100:180] .== true)

        # Check that non-outbreak periods are false
        @test all(outbreak_status[1:9] .== false)
        @test all(outbreak_status[61:79] .== false)
    end

    @testset "classify_all_outbreaks! - different specifications" begin
        inc_vec = [
            repeat([1], 9)...,
            repeat([12], 51)...,
            repeat([1], 39)...,
        ]

        # Test with stricter threshold
        outbreak_spec_strict = OutbreakSpecification(15, 30, 500)
        abovethreshold_vec = vec(inc_vec .>= outbreak_spec_strict.outbreak_threshold)
        abovethresholdrle = StatsBase.rle(abovethreshold_vec)
        all_outbreak_bounds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )
        outbreak_status = zeros(Bool, length(inc_vec))
        classify_all_outbreaks!(
            outbreak_status, all_outbreak_bounds, inc_vec, outbreak_spec_strict
        )
        @test all(outbreak_status .== false)  # 12 < 15, no outbreak

        # Test with looser duration requirement
        outbreak_spec_loose = OutbreakSpecification(5, 10, 500)
        abovethreshold_vec = vec(inc_vec .>= outbreak_spec_loose.outbreak_threshold)
        abovethresholdrle = StatsBase.rle(abovethreshold_vec)
        all_outbreak_bounds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )
        outbreak_status = zeros(Bool, length(inc_vec))
        classify_all_outbreaks!(
            outbreak_status, all_outbreak_bounds, inc_vec, outbreak_spec_loose
        )
        @test all(outbreak_status[10:60] .== true)  # Still outbreak
    end

    @testset "classify_all_outbreaks! - edge cases" begin
        # Test with no periods above threshold
        inc_vec = repeat([1], 100)
        outbreak_spec = OutbreakSpecification(5, 30, 500)

        abovethreshold_vec = vec(inc_vec .>= outbreak_spec.outbreak_threshold)
        abovethresholdrle = StatsBase.rle(abovethreshold_vec)
        all_outbreak_bounds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )
        outbreak_status = zeros(Bool, length(inc_vec))
        classify_all_outbreaks!(
            outbreak_status, all_outbreak_bounds, inc_vec, outbreak_spec
        )
        @test all(outbreak_status .== false)

        # Test with entire period above threshold
        inc_vec = repeat([15], 100)
        abovethreshold_vec = vec(inc_vec .>= outbreak_spec.outbreak_threshold)
        abovethresholdrle = StatsBase.rle(abovethreshold_vec)
        all_outbreak_bounds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 5
        )
        outbreak_status = zeros(Bool, length(inc_vec))
        classify_all_outbreaks!(
            outbreak_status, all_outbreak_bounds, inc_vec, outbreak_spec
        )
        @test all(outbreak_status .== true)
    end
end
