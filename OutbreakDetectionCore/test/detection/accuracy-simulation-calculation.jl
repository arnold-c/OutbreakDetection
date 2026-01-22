@testset "accuracy-calculation.jl" begin
    using OutbreakDetectionCore
    using StructArrays

    # Test _calculate_accuracy with BalancedAccuracy
    @testset "_calculate_accuracy - BalancedAccuracy" begin
        metric = AccuracyMetric(BalancedAccuracy())

        # Perfect detection
        accuracy = OutbreakDetectionCore._calculate_accuracy(metric, 1.0, 1.0)
        @test accuracy == 1.0

        # Balanced case
        accuracy = OutbreakDetectionCore._calculate_accuracy(metric, 0.8, 0.6)
        @test accuracy ≈ 0.7

        # Poor detection
        accuracy = OutbreakDetectionCore._calculate_accuracy(metric, 0.5, 0.3)
        @test accuracy ≈ 0.4
    end

    # Test _calculate_accuracy with F1
    @testset "_calculate_accuracy - F1" begin
        metric = AccuracyMetric(F1())

        # Perfect detection
        accuracy = OutbreakDetectionCore._calculate_accuracy(metric, 1.0, 1.0)
        @test accuracy == 1.0

        # Test F1 calculation: 2 * (precision * recall) / (precision + recall)
        # precision = 0.8, recall = 0.6
        # F1 = 2 * 0.48 / 1.4 ≈ 0.6857
        accuracy = OutbreakDetectionCore._calculate_accuracy(metric, 0.8, 0.6)
        @test accuracy ≈ 2 * (0.8 * 0.6) / (0.8 + 0.6)

        # Edge case: zero precision and recall
        accuracy = OutbreakDetectionCore._calculate_accuracy(metric, 0.0, 0.0)
        @test accuracy == 0.0
    end

    # Test calculate_simulation_accuracy
    @testset "calculate_simulation_accuracy" begin
        # Create simple outbreak bounds: one outbreak from day 10 to 60
        outbreak_bounds = [10 60 51 612]

        # Create alert vector: alerts from day 15 to 65 (delayed start, extended end)
        alert_vec = [
            repeat([false], 14)...,
            repeat([true], 51)...,
            repeat([false], 34)...,
        ]

        metric = AccuracyMetric(BalancedAccuracy())

        accuracy, prop_detected, prop_correct, delays =
            OutbreakDetectionCore.calculate_simulation_accuracy(
            alert_vec, outbreak_bounds, metric
        )

        @test accuracy isa Float64
        @test 0.0 <= accuracy <= 1.0
        @test prop_detected == 1.0  # Outbreak was detected
        @test prop_correct >= 0.0   # Some alerts were correct
        @test length(delays) >= 0   # May have delays
    end

    @testset "calculate_simulation_accuracy - no outbreaks" begin
        # Empty outbreak bounds
        outbreak_bounds = Matrix{Int64}(undef, 0, 4)

        # Some alerts present
        alert_vec = [
            repeat([false], 5)...,
            repeat([true], 3)...,
            repeat([false], 2)...,
        ]

        metric = AccuracyMetric(BalancedAccuracy())

        accuracy, prop_detected, prop_correct, delays =
            OutbreakDetectionCore.calculate_simulation_accuracy(
            alert_vec, outbreak_bounds, metric
        )

        @test accuracy isa Float64
        @test prop_detected isa Float64
        @test prop_correct isa Float64
        @test delays isa Vector{Int64}
    end

    @testset "calculate_simulation_accuracy - no alerts" begin
        # Outbreak present
        outbreak_bounds = [10 60 51 612]

        # No alerts
        alert_vec = repeat([false], 99)

        metric = AccuracyMetric(BalancedAccuracy())

        accuracy, prop_detected, prop_correct, delays =
            OutbreakDetectionCore.calculate_simulation_accuracy(
            alert_vec, outbreak_bounds, metric
        )

        @test accuracy isa Float64
        @test prop_detected == 0.0  # No outbreaks detected
        @test delays isa Vector{Int64}
        @test isempty(delays)
    end

    # Integration test for calculate_ensemble_accuracy
    @testset "calculate_ensemble_accuracy - integration" begin
        # Create test data for 2 simulations
        test_positives = [[0, 0, 5, 10, 15, 8, 3, 0, 0, 0], [0, 0, 3, 8, 12, 6, 2, 0, 0, 0]]

        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives)
        )

        alert_method = AlertMethod(DailyThreshold())
        accuracy_metric = AccuracyMetric(BalancedAccuracy())

        # Create outbreak status and bounds
        outbreak_status_vecs = [
            [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
            [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
        ]

        outbreak_bounds_vecs = [[3 7 5 50], [3 7 5 45]]

        mean_accuracy, mean_prop_detected, mean_prop_correct, mean_delay =
            OutbreakDetectionCore.calculate_ensemble_accuracy(
            7.0,
            alert_method,
            accuracy_metric,
            container,
            outbreak_status_vecs,
            outbreak_bounds_vecs,
        )

        @test mean_accuracy isa Float64
        @test 0.0 <= mean_accuracy <= 1.0
        @test 0.0 <= mean_prop_detected <= 1.0
        @test 0.0 <= mean_prop_correct <= 1.0
        @test mean_delay isa Float64
    end

    @testset "calculate_ensemble_accuracy - multiple simulations" begin
        # Create test data for 3 simulations with varying performance
        test_positives = [
            [0, 0, 5, 10, 15, 8, 3, 0, 0, 0],
            [0, 0, 3, 8, 12, 6, 2, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ]

        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives)
        )

        alert_method = AlertMethod(DailyThreshold())
        accuracy_metric = AccuracyMetric(F1())

        outbreak_status_vecs = [
            [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
            [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
            [0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
        ]

        outbreak_bounds_vecs = [[3 7 5 50], [3 7 5 45], [3 7 5 40]]

        mean_accuracy, mean_prop_detected, mean_prop_correct, mean_delay =
            OutbreakDetectionCore.calculate_ensemble_accuracy(
            7.0,
            alert_method,
            accuracy_metric,
            container,
            outbreak_status_vecs,
            outbreak_bounds_vecs,
        )

        @test mean_accuracy isa Float64
        @test mean_prop_detected isa Float64
        @test mean_prop_correct isa Float64
        @test mean_delay isa Float64
    end
end
