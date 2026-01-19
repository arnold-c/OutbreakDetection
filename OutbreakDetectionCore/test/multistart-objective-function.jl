@testset "multistart-objective-function.jl" begin
    using OutbreakDetectionCore
    using StatsBase
    using StructArrays

    # Test generate_alerts with DailyThreshold
    @testset "generate_alerts - DailyThreshold" begin
        test_positives = [0, 0, 5, 10, 15, 8, 3, 0, 0, 0]
        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives)
        )

        # Test with threshold = 7
        alerts = OutbreakDetectionCore.generate_alerts(container, 7.0)

        @test alerts isa AbstractVector{Bool}
        @test length(alerts) == 10
        @test alerts == [false, false, false, true, true, true, false, false, false, false]

        # Test with threshold = 12
        alerts_high = OutbreakDetectionCore.generate_alerts(container, 12.0)
        @test alerts_high == [
            false,
            false,
            false,
            false,
            true,
            false,
            false,
            false,
            false,
            false,
        ]

        # Test with threshold = 20 (no alerts)
        alerts_none = OutbreakDetectionCore.generate_alerts(container, 20.0)
        @test all(.!alerts_none)
    end

    # Test generate_alerts with MovingAverage
    @testset "generate_alerts - MovingAverage" begin
        test_positives_ma = [0.0, 0.0, 2.5, 7.5, 10.0, 11.0, 9.0, 5.5, 1.5, 0.0]
        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives_ma)
        )

        # Test with threshold = 8.0
        alerts = OutbreakDetectionCore.generate_alerts(container, 8.0)

        @test alerts isa AbstractVector{Bool}
        @test length(alerts) == 10
        @test alerts == [
            false,
            false,
            false,
            false,
            true,
            true,
            true,
            false,
            false,
            false,
        ]
    end

    # Test generate_alerts with DailyThresholdMovingAverage
    @testset "generate_alerts - DailyThresholdMovingAverage" begin
        test_positives_daily = [0, 0, 5, 10, 15, 8, 3, 0, 0, 0]
        test_positives_ma = [0.0, 0.0, 2.5, 7.5, 10.0, 11.0, 9.0, 5.5, 1.5, 0.0]

        container = TestPositiveContainer(
            DualAlertTestPositiveContainer(
                test_positives = test_positives_daily,
                movingavg_test_positives = test_positives_ma,
            ),
        )

        # Test with threshold = 9.0
        # Should trigger on daily: indices 4,5 (10,15)
        # Should trigger on MA: indices 5,6 (11.0, 9.0)
        # Combined (OR): indices 4,5,6
        alerts = OutbreakDetectionCore.generate_alerts(container, 9.0)

        @test alerts isa AbstractVector{Bool}
        @test length(alerts) == 10
        # Daily: 10>=9 (idx 4), 15>=9 (idx 5), 8<9 (idx 6)
        # MA: 10.0>=9 (idx 5), 11.0>=9 (idx 6), 9.0>=9 (idx 7)
        # Combined OR: indices 4,5,6,7
        @test alerts == [
            false,
            false,
            false,
            true,
            true,
            true,
            true,
            false,
            false,
            false,
        ]
    end

    # Test calculate_accuracy with BalancedAccuracy
    @testset "calculate_accuracy - BalancedAccuracy" begin
        metric = AccuracyMetric(BalancedAccuracy())

        # Perfect detection
        accuracy = OutbreakDetectionCore.calculate_accuracy(metric, 1.0, 1.0)
        @test accuracy == 1.0

        # Balanced case
        accuracy = OutbreakDetectionCore.calculate_accuracy(metric, 0.8, 0.6)
        @test accuracy ≈ 0.7

        # Poor detection
        accuracy = OutbreakDetectionCore.calculate_accuracy(metric, 0.5, 0.3)
        @test accuracy ≈ 0.4
    end

    # Test calculate_accuracy with F1
    @testset "calculate_accuracy - F1" begin
        metric = AccuracyMetric(F1())

        # Perfect detection
        accuracy = OutbreakDetectionCore.calculate_accuracy(metric, 1.0, 1.0)
        @test accuracy == 1.0

        # Test F1 calculation: 2 * (precision * recall) / (precision + recall)
        # precision = 0.8, recall = 0.6
        # F1 = 2 * 0.48 / 1.4 ≈ 0.6857
        accuracy = OutbreakDetectionCore.calculate_accuracy(metric, 0.8, 0.6)
        @test accuracy ≈ 2 * (0.8 * 0.6) / (0.8 + 0.6)

        # Edge case: zero precision and recall
        accuracy = OutbreakDetectionCore.calculate_accuracy(metric, 0.0, 0.0)
        @test accuracy == 0.0
    end

    # Test calculate_simulation_accuracy
    @testset "calculate_simulation_accuracy" begin
        # Create simple outbreak bounds: one outbreak from day 10 to 60
        outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
            lower_bounds = [10],
            upper_bounds = [60],
            duration = [51],
            num_infections_during_bounds = [612]
        )

        # Create alert vector: alerts from day 15 to 65 (delayed start, extended end)
        alert_vec = [
            repeat([false], 14)...,
            repeat([true], 51)...,
            repeat([false], 34)...,
        ]

        metric = AccuracyMetric(BalancedAccuracy())

        accuracy = OutbreakDetectionCore.calculate_simulation_accuracy(
            alert_vec, outbreak_thresholds, metric
        )

        @test accuracy isa Float64
        @test 0.0 <= accuracy <= 1.0
    end

    # Integration test for calculate_ensemble_accuracy
    @testset "calculate_ensemble_accuracy - integration" begin
        # Create test data for 2 simulations
        test_positives_vec = [[0, 0, 5, 10, 15, 8, 3, 0, 0, 0], [0, 0, 3, 8, 12, 6, 2, 0, 0, 0]]

        alert_method = AlertMethod(DailyThreshold())
        container = OutbreakDetectionCore.create_test_positive_container(
            alert_method,
            test_positives_vec
        )

        accuracy_metric = AccuracyMetric(BalancedAccuracy())

        # Create outbreak thresholds
        outbreak_thresholds = StructArrays.StructVector(
            [
                OutbreakDetectionCore.OutbreakThresholds(
                    lower_bounds = [3],
                    upper_bounds = [7],
                    duration = [5],
                    num_infections_during_bounds = [50]
                ),
                OutbreakDetectionCore.OutbreakThresholds(
                    lower_bounds = [3],
                    upper_bounds = [7],
                    duration = [5],
                    num_infections_during_bounds = [45]
                ),
            ]
        )

        mean_accuracy = OutbreakDetectionCore.calculate_ensemble_accuracy(
            7.0,
            accuracy_metric,
            container,
            outbreak_thresholds,
        )

        @test mean_accuracy isa Float64
        @test 0.0 <= mean_accuracy <= 1.0
    end
end
