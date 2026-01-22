@testset "alert-generation.jl" begin
    using OutbreakDetectionCore

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
        # Daily: 10>=9 (idx 4), 15>=9 (idx 5), 8<9 (idx 6)
        # MA: 10.0>=9 (idx 5), 11.0>=9 (idx 6), 9.0>=9 (idx 7)
        # Combined OR: indices 4,5,6,7
        alerts = OutbreakDetectionCore.generate_alerts(container, 9.0)

        @test alerts isa AbstractVector{Bool}
        @test length(alerts) == 10
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

    @testset "generate_alerts - edge cases" begin
        # Test with all zeros
        test_positives = [0, 0, 0, 0, 0]
        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives)
        )

        alerts = OutbreakDetectionCore.generate_alerts(container, 1.0)
        @test all(.!alerts)

        # Test with all above threshold
        test_positives = [10, 10, 10, 10, 10]
        container = TestPositiveContainer(
            SingleAlertTestPositiveContainer(test_positive_results = test_positives)
        )

        alerts = OutbreakDetectionCore.generate_alerts(container, 5.0)
        @test all(alerts)
    end
end
