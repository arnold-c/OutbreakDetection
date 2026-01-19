@testset "threshold-optimization.jl" begin
    using OutbreakDetectionCore
    using StatsBase
    using StructArrays

    # Test data setup
    @testset "Test data setup" begin
        # Create simple test positive data
        test_positives_daily = [
            [0, 0, 5, 10, 15, 8, 3, 0, 0, 0],
            [0, 0, 3, 8, 12, 6, 2, 0, 0, 0],
        ]

        test_positives_ma = [
            [0.0, 0.0, 2.5, 7.5, 10.0, 11.0, 9.0, 5.5, 1.5, 0.0],
            [0.0, 0.0, 1.5, 5.5, 7.5, 8.5, 6.5, 4.0, 1.0, 0.0],
        ]

        @test length(test_positives_daily) == 2
        @test length(test_positives_ma) == 2
    end

    # Test generate_alerts with DailyThreshold
    @testset "generate_alerts - DailyThreshold" begin
        test_positives = [[0, 0, 5, 10, 15, 8, 3, 0, 0, 0]]
        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives)
        )
        alert_method = AlertMethod(DailyThreshold())

        # Test with threshold = 7
        alerts = OutbreakDetectionCore.generate_alerts(7.0, alert_method, container, 1)

        @test alerts isa AbstractVector{Bool}
        @test length(alerts) == 10
        @test alerts == [false, false, false, true, true, true, false, false, false, false]

        # Test with threshold = 12
        alerts_high = OutbreakDetectionCore.generate_alerts(
            12.0, alert_method, container, 1
        )
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
        alerts_none = OutbreakDetectionCore.generate_alerts(
            20.0, alert_method, container, 1
        )
        @test all(.!alerts_none)
    end

    # Test generate_alerts with MovingAverage
    @testset "generate_alerts - MovingAverage" begin
        test_positives_ma = [[0.0, 0.0, 2.5, 7.5, 10.0, 11.0, 9.0, 5.5, 1.5, 0.0]]
        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives_ma)
        )
        alert_method = AlertMethod(MovingAverage(window = 3))

        # Test with threshold = 8.0
        alerts = OutbreakDetectionCore.generate_alerts(8.0, alert_method, container, 1)

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
        test_positives_daily = [[0, 0, 5, 10, 15, 8, 3, 0, 0, 0]]
        test_positives_ma = [[0.0, 0.0, 2.5, 7.5, 10.0, 11.0, 9.0, 5.5, 1.5, 0.0]]

        container = TestPositiveContainer(
            DualAlertTestPositiveContainer(
                test_positives = test_positives_daily,
                movingavg_test_positives = test_positives_ma,
            ),
        )
        alert_method = AlertMethod(DailyThresholdMovingAverage(window = 3))

        # Test with threshold = 9.0
        # Should trigger on daily: indices 4,5 (10,15)
        # Should trigger on MA: indices 5,6 (11.0, 9.0)
        # Combined (OR): indices 4,5,6
        alerts = OutbreakDetectionCore.generate_alerts(9.0, alert_method, container, 1)

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

    # Test create_outbreak_status_vecs
    @testset "create_outbreak_status_vecs" begin
        using StaticArrays

        # Create mock SEIR results
        inc1 = [
            repeat([1], 9)...,
            repeat([12], 51)..., # Outbreak
            repeat([1], 39)...,
        ]
        inc2 = [
            repeat([2], 19)...,
            repeat([20], 41)..., # Outbreak
            repeat([2], 39)...,
        ]

        # Create proper SEIRRun objects
        states1 = [SVector{5, Int64}(0, 0, 0, 0, 0) for _ in 1:99]
        states2 = [SVector{5, Int64}(0, 0, 0, 0, 0) for _ in 1:99]
        reff1 = zeros(Float64, 99)
        reff2 = zeros(Float64, 99)

        seir_results = StructVector(
            [
                OutbreakDetectionCore.SEIRRun(
                    states = states1, incidence = inc1, Reff = reff1
                ),
                OutbreakDetectionCore.SEIRRun(
                    states = states2, incidence = inc2, Reff = reff2
                ),
            ]
        )

        outbreak_spec = OutbreakSpecification(5, 30, 500)

        outbreak_status_vecs, outbreak_bounds_vecs = create_outbreak_status_vecs(
            seir_results, outbreak_spec
        )

        @test length(outbreak_status_vecs) == 2
        @test length(outbreak_bounds_vecs) == 2

        # Check first simulation
        @test length(outbreak_status_vecs[1]) == 99
        @test all(outbreak_status_vecs[1][10:60] .== true)  # Outbreak period
        @test all(outbreak_status_vecs[1][1:9] .== false)    # Before outbreak

        # Check outbreak bounds for first simulation
        @test size(outbreak_bounds_vecs[1], 2) == 4  # [start, end, duration, size]
        @test outbreak_bounds_vecs[1][1, 1] == 10    # Start
        @test outbreak_bounds_vecs[1][1, 2] == 60    # End

        # Check second simulation
        @test length(outbreak_status_vecs[2]) == 99
        @test all(outbreak_status_vecs[2][20:60] .== true)  # Outbreak period
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
end
