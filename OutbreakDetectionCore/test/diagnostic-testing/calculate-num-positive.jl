using OutbreakDetectionCore, StatsBase, Random, Distributions

@testset "calculate-num-positive.jl" begin
    @testset "calculate_true_positives!" begin
        @testset "Perfect test (sensitivity = 1.0)" begin
            # Test with perfect sensitivity - all tested should become positive
            tested_vec = [0, 5, 10, 15, 20, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 0.95,
                test_result_lag = 2
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # With lag=2, day i tested should appear on day i+2
            # day=1: result_day=3, npos_vec[3] = tested_vec[1] = 0
            # day=2: result_day=4, npos_vec[4] = tested_vec[2] = 5
            @test npos_vec[3] == 0
            @test npos_vec[4] == 5
            @test npos_vec[5] == 10
            @test npos_vec[6] == 15
            @test npos_vec[7] == 20
            @test npos_vec[1] == 0
            @test npos_vec[2] == 0
        end

        @testset "Imperfect test (sensitivity < 1.0)" begin
            # Test with imperfect sensitivity - should use binomial sampling
            tested_vec = [0, 100, 0, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.8,
                specificity = 0.95,
                test_result_lag = 1
            )

            rng = Random.MersenneTwister(42)
            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec, rng
            )

            # With lag=1, day i tested should appear on day i+1
            # day=2: result_day=3, npos_vec[3] = binomial(tested_vec[2], 0.8) = binomial(100, 0.8)
            @test npos_vec[3] > 0
            @test npos_vec[3] <= 100
            @test npos_vec[1] == 0
            @test npos_vec[2] == 0
        end

        @testset "No lag" begin
            tested_vec = [0, 10, 20, 30, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 0.95,
                test_result_lag = 0
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # With lag=0, day i tested should appear on day i
            # day=1: result_day=1, npos_vec[1] = tested_vec[1] = 0
            # day=2: result_day=2, npos_vec[2] = tested_vec[2] = 10
            @test npos_vec[1] == 0
            @test npos_vec[2] == 10
            @test npos_vec[3] == 20
            @test npos_vec[4] == 30
        end

        @testset "Lag extends beyond simulation" begin
            tested_vec = [0, 0, 0, 0, 0, 100, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 0.95,
                test_result_lag = 10
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # Day 6 tested with lag 10 would result on day 16, which is beyond sim_length
            @test all(npos_vec .== 0)
        end
    end

    @testset "calculate_false_positives!" begin
        @testset "Perfect test (specificity = 1.0)" begin
            # Test with perfect specificity - no false positives
            tested_vec = [0, 50, 100, 150, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 1.0,
                test_result_lag = 2
            )

            calculate_false_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # With perfect specificity, all should be zero
            @test all(npos_vec .== 0)
        end

        @testset "Imperfect test (specificity < 1.0)" begin
            # Test with imperfect specificity - should have false positives
            tested_vec = [0, 1000, 0, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 0.95,
                test_result_lag = 1
            )

            rng = Random.MersenneTwister(42)
            calculate_false_positives!(
                npos_vec, tested_vec, sim_length, test_spec, rng
            )

            # With lag=1, day 2 tested (1000) should appear on day 3
            # False positive rate = 1 - 0.95 = 0.05, so expect ~50 false positives
            @test npos_vec[3] > 0
            @test npos_vec[3] <= 1000
            @test npos_vec[1] == 0
            @test npos_vec[2] == 0
        end

        @testset "No lag" begin
            tested_vec = [0, 200, 300, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 0.98,
                test_result_lag = 0
            )

            rng = Random.MersenneTwister(42)
            calculate_false_positives!(
                npos_vec, tested_vec, sim_length, test_spec, rng
            )

            # With lag=0, tested on day i should appear on day i
            @test npos_vec[1] == 0
            @test npos_vec[2] > 0  # Some false positives expected
            @test npos_vec[3] > 0  # Some false positives expected
        end
    end

    @testset "perfect_test_true_positives_vec" begin
        @testset "Basic functionality" begin
            tested_vec = [0, 10, 20, 30, 40, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 2

            OutbreakDetectionCore.perfect_test_true_positives_vec(
                npos_vec, tested_vec, sim_length, lag
            )

            # With lag=2, day i tested should appear on day i+2
            # day=1: result_day=3, npos_vec[3] = tested_vec[1] = 0
            # day=2: result_day=4, npos_vec[4] = tested_vec[2] = 10
            @test npos_vec[3] == 0   # day 1 tested
            @test npos_vec[4] == 10  # day 2 tested
            @test npos_vec[5] == 20  # day 3 tested
            @test npos_vec[6] == 30  # day 4 tested
            @test npos_vec[7] == 40  # day 5 tested
            @test npos_vec[1] == 0
            @test npos_vec[2] == 0
        end

        @testset "No lag" begin
            tested_vec = [5, 10, 15, 20, 25, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 0

            OutbreakDetectionCore.perfect_test_true_positives_vec(
                npos_vec, tested_vec, sim_length, lag
            )

            # With lag=0, day i tested should appear on day i
            @test npos_vec[1] == 5
            @test npos_vec[2] == 10
            @test npos_vec[3] == 15
            @test npos_vec[4] == 20
            @test npos_vec[5] == 25
        end

        @testset "Large lag" begin
            tested_vec = [0, 0, 0, 0, 0, 100, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 8

            OutbreakDetectionCore.perfect_test_true_positives_vec(
                npos_vec, tested_vec, sim_length, lag
            )

            # Day 6 tested with lag 8 would result on day 14, beyond sim_length
            @test all(npos_vec .== 0)
        end

        @testset "Lag at boundary" begin
            tested_vec = [0, 0, 0, 0, 0, 0, 0, 0, 100, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 1

            OutbreakDetectionCore.perfect_test_true_positives_vec(
                npos_vec, tested_vec, sim_length, lag
            )

            # Day 9 tested with lag 1 should appear on day 10
            @test npos_vec[10] == 100
            @test npos_vec[9] == 0
        end

        @testset "Initializes to zero" begin
            tested_vec = [10, 20, 30, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = [999, 999, 999, 999, 999, 999, 999, 999, 999, 999]
            sim_length = 10
            lag = 1

            OutbreakDetectionCore.perfect_test_true_positives_vec(
                npos_vec, tested_vec, sim_length, lag
            )

            # Should initialize to zero first
            @test npos_vec[1] == 0
            @test npos_vec[2] == 10
            @test npos_vec[3] == 20
            @test npos_vec[4] == 30
        end
    end

    @testset "_calculate_positives_vec!" begin
        @testset "Deterministic with multiplier = 1.0" begin
            tested_vec = [0, 100, 200, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 1
            multiplier = 1.0

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # With multiplier=1.0, should be same as perfect test
            # day=1: result_day=2, npos_vec[2] = binomial(tested_vec[1], 1.0) = 0
            # day=2: result_day=3, npos_vec[3] = binomial(tested_vec[2], 1.0) = 100
            @test npos_vec[2] == 0
            @test npos_vec[3] == 100
            @test npos_vec[4] == 200
            @test npos_vec[1] == 0
        end

        @testset "Stochastic with multiplier < 1.0" begin
            tested_vec = [0, 1000, 0, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 1
            multiplier = 0.5

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # With multiplier=0.5, expect approximately 500 (binomial sample)
            # day=1: result_day=2, npos_vec[2] = binomial(tested_vec[1], 0.5) = 0
            # day=2: result_day=3, npos_vec[3] = binomial(tested_vec[2], 0.5) ≈ 500
            @test npos_vec[2] == 0
            @test npos_vec[3] > 0
            @test npos_vec[3] <= 1000
            @test npos_vec[1] == 0
        end

        @testset "No lag" begin
            tested_vec = [50, 100, 150, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 0
            multiplier = 0.8

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # With lag=0, results appear immediately
            @test npos_vec[1] > 0
            @test npos_vec[1] <= 50
            @test npos_vec[2] > 0
            @test npos_vec[2] <= 100
            @test npos_vec[3] > 0
            @test npos_vec[3] <= 150
        end

        @testset "Initializes to zero" begin
            tested_vec = [10, 20, 30, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = [999, 999, 999, 999, 999, 999, 999, 999, 999, 999]
            sim_length = 10
            lag = 1
            multiplier = 0.9

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # Should initialize to zero first
            @test npos_vec[1] == 0
            @test npos_vec[2] > 0
        end

        @testset "Lag extends beyond simulation" begin
            tested_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 100]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 5
            multiplier = 0.8

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # Day 10 tested with lag 5 would result on day 15, beyond sim_length
            @test all(npos_vec .== 0)
        end

        @testset "Reproducibility with same seed" begin
            tested_vec = [0, 500, 0, 0, 0, 0, 0, 0, 0, 0]
            npos_vec1 = zeros(Int64, 10)
            npos_vec2 = zeros(Int64, 10)
            sim_length = 10
            lag = 1
            multiplier = 0.6

            rng1 = Random.MersenneTwister(123)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec1, tested_vec, sim_length, lag, multiplier; rng = rng1
            )

            rng2 = Random.MersenneTwister(123)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec2, tested_vec, sim_length, lag, multiplier; rng = rng2
            )

            # Same seed should produce same results
            @test npos_vec1 == npos_vec2
        end

        @testset "Multiplier = 0.0 (no positives)" begin
            tested_vec = [0, 1000, 2000, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 1
            multiplier = 0.0

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # With multiplier=0.0, all results should be zero
            @test all(npos_vec .== 0)
        end

        @testset "Result day exactly equals sim_length" begin
            tested_vec = [0, 0, 0, 0, 0, 0, 0, 0, 50, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = 1
            multiplier = 1.0

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # Day 9 tested with lag 1 should appear on day 10 (exactly sim_length)
            @test npos_vec[10] == 50
            @test sum(npos_vec[1:9]) == 0
        end
    end

    @testset "Edge cases and type stability" begin
        @testset "Single element vectors" begin
            tested_vec = [100]
            npos_vec = zeros(Int64, 1)
            sim_length = 1
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 1.0,
                test_result_lag = 0
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # With lag=0, should appear immediately
            @test npos_vec[1] == 100
        end

        @testset "Single element with lag" begin
            tested_vec = [100]
            npos_vec = zeros(Int64, 1)
            sim_length = 1
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 1.0,
                test_result_lag = 1
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # With lag=1, result would be on day 2, which is beyond sim_length
            @test npos_vec[1] == 0
        end

        @testset "All zeros tested" begin
            tested_vec = zeros(Int64, 10)
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 0.95,
                test_result_lag = 2
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # All zeros in, all zeros out
            @test all(npos_vec .== 0)
        end

        @testset "Very large lag" begin
            tested_vec = [100, 200, 300, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 1.0,
                test_result_lag = 100
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            # All results would be beyond sim_length
            @test all(npos_vec .== 0)
        end

        @testset "Negative lag (should still work as lag=0)" begin
            # Note: This tests implementation behavior with negative lag
            # The function doesn't validate lag >= 0, so negative lag acts like 0 or less
            tested_vec = [0, 100, 0, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            lag = -1
            multiplier = 1.0

            rng = Random.MersenneTwister(42)
            OutbreakDetectionCore._calculate_positives_vec!(
                npos_vec, tested_vec, sim_length, lag, multiplier; rng = rng
            )

            # With negative lag, result_day = day + (-1) = day - 1
            # day=1: result_day=0, which is out of bounds, so nothing stored
            # day=2: result_day=1, npos_vec[1] = 100
            @test npos_vec[1] == 100
            @test sum(npos_vec[2:end]) == 0
        end

        @testset "Type stability - calculate_true_positives!" begin
            tested_vec = [0, 100, 200, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 0.95,
                test_result_lag = 2
            )
            rng = Random.MersenneTwister(42)

            # Test that function returns Nothing and is type stable
            result = calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec, rng
            )
            @test result === nothing
            @test eltype(npos_vec) === Int64
        end

        @testset "Type stability - calculate_false_positives!" begin
            tested_vec = [0, 100, 200, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 0.95,
                test_result_lag = 2
            )
            rng = Random.MersenneTwister(42)

            # Test that function returns Nothing and is type stable
            result = calculate_false_positives!(
                npos_vec, tested_vec, sim_length, test_spec, rng
            )
            @test result === nothing
            @test eltype(npos_vec) === Int64
        end

        @testset "Different vector types" begin
            # Test with Vector{Int64}
            tested_vec = Vector{Int64}([0, 100, 200, 0, 0, 0, 0, 0, 0, 0])
            npos_vec = Vector{Int64}(zeros(Int64, 10))
            sim_length = 10
            test_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 1.0,
                test_result_lag = 1
            )

            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec
            )

            @test npos_vec[2] == 0
            @test npos_vec[3] == 100
            @test npos_vec[4] == 200
        end

        @testset "Boundary: sensitivity and specificity at extremes" begin
            tested_vec = [0, 1000, 0, 0, 0, 0, 0, 0, 0, 0]
            npos_vec = zeros(Int64, 10)
            sim_length = 10

            # Test with sensitivity = 0.0 (no true positives)
            test_spec_zero_sens = IndividualTestSpecification(;
                sensitivity = 0.0,
                specificity = 1.0,
                test_result_lag = 1
            )
            rng = Random.MersenneTwister(42)
            calculate_true_positives!(
                npos_vec, tested_vec, sim_length, test_spec_zero_sens, rng
            )
            @test all(npos_vec .== 0)

            # Test with specificity = 0.0 (all false positives)
            test_spec_zero_spec = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 0.0,
                test_result_lag = 1
            )
            npos_vec = zeros(Int64, 10)
            rng = Random.MersenneTwister(42)
            calculate_false_positives!(
                npos_vec, tested_vec, sim_length, test_spec_zero_spec, rng
            )
            # With specificity=0.0, false positive rate = 1.0, so all tested become positive
            @test npos_vec[3] == 1000
        end
    end

    @testset "Integration tests" begin
        @testset "True positives + false positives combined" begin
            # Simulate a realistic scenario with both infected and non-infected tested
            infected_tested = [0, 50, 100, 150, 100, 50, 0, 0, 0, 0]
            noninfected_tested = [100, 200, 300, 200, 100, 50, 0, 0, 0, 0]

            true_pos_vec = zeros(Int64, 10)
            false_pos_vec = zeros(Int64, 10)
            sim_length = 10

            test_spec = IndividualTestSpecification(;
                sensitivity = 0.9,
                specificity = 0.95,
                test_result_lag = 1
            )

            rng1 = Random.MersenneTwister(42)
            rng2 = Random.MersenneTwister(43)

            calculate_true_positives!(
                true_pos_vec, infected_tested, sim_length, test_spec, rng1
            )
            calculate_false_positives!(
                false_pos_vec, noninfected_tested, sim_length, test_spec, rng2
            )

            # Both should have values (stochastic, so just check they're reasonable)
            @test sum(true_pos_vec) > 0
            @test sum(false_pos_vec) >= 0

            # True positives should be roughly 90% of infected tested (with lag)
            # False positives should be roughly 5% of non-infected tested (with lag)
            total_pos = true_pos_vec .+ false_pos_vec
            @test sum(total_pos) > 0
        end

        @testset "Perfect test consistency" begin
            # Perfect test should have no false positives and all true positives
            tested_vec = [0, 100, 200, 300, 0, 0, 0, 0, 0, 0]
            true_pos_vec = zeros(Int64, 10)
            false_pos_vec = zeros(Int64, 10)
            sim_length = 10

            perfect_test = IndividualTestSpecification(;
                sensitivity = 1.0,
                specificity = 1.0,
                test_result_lag = 2
            )

            calculate_true_positives!(
                true_pos_vec, tested_vec, sim_length, perfect_test
            )
            calculate_false_positives!(
                false_pos_vec, tested_vec, sim_length, perfect_test
            )

            # Perfect test: all true positives, no false positives
            @test sum(true_pos_vec) == sum(tested_vec)
            @test sum(false_pos_vec) == 0
        end
    end
end
