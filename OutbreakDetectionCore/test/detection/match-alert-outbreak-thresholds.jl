@testset "match-alert-outbreak-bounds.jl" begin
    using OutbreakDetectionCore

    @testset "MatchedThresholds: Perfect overlap" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10, 50],
            upper_bounds = [30, 70],
            duration = [21, 21],
            num_infections_during_bounds = [100, 200]
        )

        alerts = Thresholds(
            lower_bounds = [15, 55],
            upper_bounds = [25, 65],
            duration = [11, 11]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test matched.outbreak_indices_with_alerts == [1, 2]
        @test matched.alert_indices_per_outbreak == [[1], [2]]
        @test matched.n_outbreaks == 2
        @test matched.n_alerts == 2
        @test calculate_sensitivity(matched) == 1.0
        @test calculate_ppv(matched) == 1.0
        @test get_matched_alert_indices(matched) == [1, 2]
    end

    @testset "MatchedThresholds: Partial detection" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10, 30, 50],
            upper_bounds = [20, 40, 60],
            duration = [11, 11, 11],
            num_infections_during_bounds = [100, 150, 200]
        )

        alerts = Thresholds(
            lower_bounds = [12, 35],
            upper_bounds = [18, 38],
            duration = [7, 4]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test matched.outbreak_indices_with_alerts == [1, 2]
        @test matched.alert_indices_per_outbreak == [[1], [2]]
        @test matched.n_outbreaks == 3
        @test matched.n_alerts == 2
        @test calculate_sensitivity(matched) ≈ 2 / 3
        @test calculate_ppv(matched) == 1.0
    end

    @testset "MatchedThresholds: False alarms" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10, 30],
            upper_bounds = [20, 40],
            duration = [11, 11],
            num_infections_during_bounds = [100, 150]
        )

        alerts = Thresholds(
            lower_bounds = [12, 35, 50, 70],
            upper_bounds = [18, 38, 55, 75],
            duration = [7, 4, 6, 6]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test matched.outbreak_indices_with_alerts == [1, 2]
        @test matched.alert_indices_per_outbreak == [[1], [2]]
        @test matched.n_outbreaks == 2
        @test matched.n_alerts == 4
        @test calculate_sensitivity(matched) == 1.0
        @test calculate_ppv(matched) == 0.5  # 2 correct out of 4 total
    end

    @testset "MatchedThresholds: First outbreak only rule" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10, 30, 50],
            upper_bounds = [20, 40, 60],
            duration = [11, 11, 11],
            num_infections_during_bounds = [100, 150, 200]
        )

        # Alert spans all three outbreaks
        alerts = Thresholds(
            lower_bounds = [15],
            upper_bounds = [55],
            duration = [41]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        # Alert should only match first outbreak
        @test matched.outbreak_indices_with_alerts == [1]
        @test matched.alert_indices_per_outbreak == [[1]]
        @test matched.n_outbreaks == 3
        @test matched.n_alerts == 1
        @test calculate_sensitivity(matched) ≈ 1 / 3  # Only 1 outbreak detected
        @test calculate_ppv(matched) == 1.0  # 1 alert, 1 correct
    end

    @testset "MatchedThresholds: Multiple alerts per outbreak" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10],
            upper_bounds = [50],
            duration = [41],
            num_infections_during_bounds = [300]
        )

        alerts = Thresholds(
            lower_bounds = [15, 25, 35],
            upper_bounds = [20, 30, 40],
            duration = [6, 6, 6]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test matched.outbreak_indices_with_alerts == [1]
        @test matched.alert_indices_per_outbreak == [[1, 2, 3]]
        @test matched.n_outbreaks == 1
        @test matched.n_alerts == 3
        @test calculate_sensitivity(matched) == 1.0
        @test calculate_ppv(matched) == 1.0
        @test get_matched_alert_indices(matched) == [1, 2, 3]
    end

    @testset "MatchedThresholds: No matches" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10, 30],
            upper_bounds = [20, 40],
            duration = [11, 11],
            num_infections_during_bounds = [100, 150]
        )

        alerts = Thresholds(
            lower_bounds = [50, 70],
            upper_bounds = [55, 75],
            duration = [6, 6]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test isempty(matched.outbreak_indices_with_alerts)
        @test isempty(matched.alert_indices_per_outbreak)
        @test matched.n_outbreaks == 2
        @test matched.n_alerts == 2
        @test calculate_sensitivity(matched) == 0.0
        @test calculate_ppv(matched) == 0.0
        @test isempty(get_matched_alert_indices(matched))
    end

    @testset "MatchedThresholds: Empty inputs" begin
        # No outbreaks
        outbreaks_empty = OutbreakThresholds(
            lower_bounds = Int64[],
            upper_bounds = Int64[],
            duration = Int64[],
            num_infections_during_bounds = Int64[]
        )

        alerts = Thresholds(
            lower_bounds = [10, 20],
            upper_bounds = [15, 25],
            duration = [6, 6]
        )

        matched1 = match_outbreak_detection_bounds(outbreaks_empty, alerts)
        @test isempty(matched1.outbreak_indices_with_alerts)
        @test matched1.n_outbreaks == 0
        @test matched1.n_alerts == 2
        @test isnan(calculate_sensitivity(matched1))
        @test calculate_ppv(matched1) == 0.0

        # No alerts
        outbreaks = OutbreakThresholds(
            lower_bounds = [10],
            upper_bounds = [20],
            duration = [11],
            num_infections_during_bounds = [100]
        )

        alerts_empty = Thresholds(
            lower_bounds = Int64[],
            upper_bounds = Int64[],
            duration = Int64[]
        )

        matched2 = match_outbreak_detection_bounds(outbreaks, alerts_empty)
        @test isempty(matched2.outbreak_indices_with_alerts)
        @test matched2.n_outbreaks == 1
        @test matched2.n_alerts == 0
        @test calculate_sensitivity(matched2) == 0.0
        @test isnan(calculate_ppv(matched2))
    end

    @testset "MatchedThresholds: Alert starts before outbreak" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [20],
            upper_bounds = [40],
            duration = [21],
            num_infections_during_bounds = [150]
        )

        alerts = Thresholds(
            lower_bounds = [15],
            upper_bounds = [25],
            duration = [11]
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test matched.outbreak_indices_with_alerts == [1]
        @test matched.alert_indices_per_outbreak == [[1]]
        @test calculate_sensitivity(matched) == 1.0
        @test calculate_ppv(matched) == 1.0
    end

    @testset "MatchedThresholds: Boundary conditions" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [10],
            upper_bounds = [20],
            duration = [11],
            num_infections_during_bounds = [100]
        )

        # Alert starts exactly at outbreak start
        alerts1 = Thresholds(
            lower_bounds = [10],
            upper_bounds = [15],
            duration = [6]
        )

        matched1 = match_outbreak_detection_bounds(outbreaks, alerts1)
        @test matched1.outbreak_indices_with_alerts == [1]
        @test calculate_sensitivity(matched1) == 1.0

        # Alert ends exactly at outbreak end
        alerts2 = Thresholds(
            lower_bounds = [15],
            upper_bounds = [20],
            duration = [6]
        )

        matched2 = match_outbreak_detection_bounds(outbreaks, alerts2)
        @test matched2.outbreak_indices_with_alerts == [1]
        @test calculate_sensitivity(matched2) == 1.0

        # Alert starts exactly at outbreak end + 1 (no overlap)
        alerts3 = Thresholds(
            lower_bounds = [21],
            upper_bounds = [25],
            duration = [5]
        )

        matched3 = match_outbreak_detection_bounds(outbreaks, alerts3)
        @test isempty(matched3.outbreak_indices_with_alerts)
        @test calculate_sensitivity(matched3) == 0.0
    end

    @testset "MatchedThresholds: Sparse efficiency scenario" begin
        # Many outbreaks, few detected (sparse case)
        outbreaks = OutbreakThresholds(
            lower_bounds = collect(10:10:1000),  # 100 outbreaks
            upper_bounds = collect(19:10:1009),
            duration = fill(10, 100),
            num_infections_during_bounds = fill(100, 100)
        )

        # Only 5 alerts, matching first 5 outbreaks
        alerts = Thresholds(
            lower_bounds = [12, 22, 32, 42, 52],
            upper_bounds = [18, 28, 38, 48, 58],
            duration = fill(7, 5)
        )

        matched = match_outbreak_detection_bounds(outbreaks, alerts)

        @test matched.outbreak_indices_with_alerts == [1, 2, 3, 4, 5]
        @test length(matched.alert_indices_per_outbreak) == 5
        @test matched.n_outbreaks == 100
        @test matched.n_alerts == 5
        @test calculate_sensitivity(matched) == 0.05  # 5/100
        @test calculate_ppv(matched) == 1.0  # 5/5
    end

    @testset "MatchedThresholds: Validation tests" begin
        # Test length mismatch assertion
        @test_throws AssertionError MatchedThresholds(
            outbreak_indices_with_alerts = [1, 2],
            alert_indices_per_outbreak = [[1]],  # Length mismatch
            n_matched_outbreaks = 2,
            n_matched_alerts = 1,
            n_outbreaks = 5,
            n_alerts = 10
        )

        # Test negative n_outbreaks
        @test_throws AssertionError MatchedThresholds(
            outbreak_indices_with_alerts = Int64[],
            alert_indices_per_outbreak = Vector{Int64}[],
            n_matched_outbreaks = 0,
            n_matched_alerts = 0,
            n_outbreaks = -1,
            n_alerts = 10
        )

        # Test outbreak index out of range
        @test_throws AssertionError MatchedThresholds(
            outbreak_indices_with_alerts = [1, 6],  # 6 > n_outbreaks
            alert_indices_per_outbreak = [[1], [2]],
            n_matched_outbreaks = 2,
            n_matched_alerts = 2,
            n_outbreaks = 5,
            n_alerts = 10
        )

        # Test alert index out of range
        @test_throws AssertionError MatchedThresholds(
            outbreak_indices_with_alerts = [1],
            alert_indices_per_outbreak = [[1, 11]],  # 11 > n_alerts
            n_matched_outbreaks = 1,
            n_matched_alerts = 2,
            n_outbreaks = 5,
            n_alerts = 10
        )

        # Test n_matched_outbreaks mismatch
        @test_throws AssertionError MatchedThresholds(
            outbreak_indices_with_alerts = [1, 2],
            alert_indices_per_outbreak = [[1], [2]],
            n_matched_outbreaks = 3,  # Should be 2
            n_matched_alerts = 2,
            n_outbreaks = 5,
            n_alerts = 10
        )

        # Test n_matched_alerts mismatch
        @test_throws AssertionError MatchedThresholds(
            outbreak_indices_with_alerts = [1, 2],
            alert_indices_per_outbreak = [[1], [2, 3]],
            n_matched_outbreaks = 2,
            n_matched_alerts = 5,  # Should be 3
            n_outbreaks = 5,
            n_alerts = 10
        )
    end

    @testset "Helper functions: get_matched_alert_indices" begin
        matched = MatchedThresholds(
            outbreak_indices_with_alerts = [1, 3, 5],
            alert_indices_per_outbreak = [[2, 5], [7], [9, 10, 11]],
            n_matched_outbreaks = 3,
            n_matched_alerts = 6,
            n_outbreaks = 10,
            n_alerts = 20
        )

        indices = get_matched_alert_indices(matched)
        @test indices == [2, 5, 7, 9, 10, 11]
        @test issorted(indices)
        @test allunique(indices)
    end

    @testset "AlertFilteringStrategy: AllAlerts (default behavior)" begin
        # Alert starts before outbreak but extends into it
        outbreaks = OutbreakThresholds(
            lower_bounds = [20],
            upper_bounds = [40],
            duration = [21],
            num_infections_during_bounds = [150]
        )

        alerts = Thresholds(
            lower_bounds = [15],
            upper_bounds = [25],
            duration = [11]
        )

        # Default (AllAlerts) should match
        matched_default = match_outbreak_detection_bounds(outbreaks, alerts)
        @test matched_default.outbreak_indices_with_alerts == [1]
        @test calculate_sensitivity(matched_default) == 1.0
        @test calculate_ppv(matched_default) == 1.0

        # Explicit AllAlerts should match
        matched_all = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(AllAlerts())
        )
        @test matched_all.outbreak_indices_with_alerts == [1]
        @test calculate_sensitivity(matched_all) == 1.0
        @test calculate_ppv(matched_all) == 1.0
    end

    @testset "AlertFilteringStrategy: PostOutbreakStartAlerts excludes pre-outbreak alerts" begin
        # Alert starts before outbreak but extends into it
        outbreaks = OutbreakThresholds(
            lower_bounds = [20],
            upper_bounds = [40],
            duration = [21],
            num_infections_during_bounds = [150]
        )

        alerts = Thresholds(
            lower_bounds = [15],
            upper_bounds = [25],
            duration = [11]
        )

        # PostOutbreakStartAlerts should NOT match (alert starts at 15 < outbreak start 20)
        matched_filtered = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
        )
        @test isempty(matched_filtered.outbreak_indices_with_alerts)
        @test calculate_sensitivity(matched_filtered) == 0.0
        @test calculate_ppv(matched_filtered) == 0.0
    end

    @testset "AlertFilteringStrategy: PostOutbreakStartAlerts includes post-outbreak alerts" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [20],
            upper_bounds = [40],
            duration = [21],
            num_infections_during_bounds = [150]
        )

        # Alert starts after outbreak start
        alerts = Thresholds(
            lower_bounds = [25],
            upper_bounds = [35],
            duration = [11]
        )

        # PostOutbreakStartAlerts should match (alert starts at 25 >= outbreak start 20)
        matched_filtered = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
        )
        @test matched_filtered.outbreak_indices_with_alerts == [1]
        @test calculate_sensitivity(matched_filtered) == 1.0
        @test calculate_ppv(matched_filtered) == 1.0
    end

    @testset "AlertFilteringStrategy: PostOutbreakStartAlerts boundary - alert at outbreak start" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [20],
            upper_bounds = [40],
            duration = [21],
            num_infections_during_bounds = [150]
        )

        # Alert starts exactly at outbreak start
        alerts = Thresholds(
            lower_bounds = [20],
            upper_bounds = [30],
            duration = [11]
        )

        # PostOutbreakStartAlerts should match (alert starts at 20 == outbreak start 20)
        matched_filtered = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
        )
        @test matched_filtered.outbreak_indices_with_alerts == [1]
        @test calculate_sensitivity(matched_filtered) == 1.0
        @test calculate_ppv(matched_filtered) == 1.0
    end

    @testset "AlertFilteringStrategy: Mixed alerts with PostOutbreakStartAlerts" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [20, 60],
            upper_bounds = [40, 80],
            duration = [21, 21],
            num_infections_during_bounds = [150, 200]
        )

        # Mix of pre-outbreak and post-outbreak alerts
        alerts = Thresholds(
            lower_bounds = [15, 25, 55, 65],  # 15 < 20, 25 >= 20, 55 < 60, 65 >= 60
            upper_bounds = [25, 35, 65, 75],
            duration = [11, 11, 11, 11]
        )

        # AllAlerts: all 4 alerts match (2 per outbreak)
        matched_all = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(AllAlerts())
        )
        @test matched_all.outbreak_indices_with_alerts == [1, 2]
        @test matched_all.alert_indices_per_outbreak == [[1, 2], [3, 4]]
        @test matched_all.n_matched_alerts == 4
        @test calculate_sensitivity(matched_all) == 1.0
        @test calculate_ppv(matched_all) == 1.0

        # PostOutbreakStartAlerts: only alerts 2 and 4 match
        matched_filtered = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
        )
        @test matched_filtered.outbreak_indices_with_alerts == [1, 2]
        @test matched_filtered.alert_indices_per_outbreak == [[2], [4]]
        @test matched_filtered.n_matched_alerts == 2
        @test calculate_sensitivity(matched_filtered) == 1.0 # all outbreaks detected
        @test calculate_ppv(matched_filtered) == 0.5  # 2 matched out of 4 total alerts
    end

    @testset "AlertFilteringStrategy: PostOutbreakStartAlerts with missed outbreaks" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [20, 60],
            upper_bounds = [40, 80],
            duration = [21, 21],
            num_infections_during_bounds = [150, 200]
        )

        # Only pre-outbreak alerts for first outbreak, post-outbreak for second
        alerts = Thresholds(
            lower_bounds = [15, 65],  # 15 < 20 (pre-outbreak), 65 >= 60 (post-outbreak)
            upper_bounds = [25, 75],
            duration = [11, 11]
        )

        # AllAlerts: both outbreaks detected
        matched_all = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(AllAlerts())
        )
        @test matched_all.outbreak_indices_with_alerts == [1, 2]
        @test calculate_sensitivity(matched_all) == 1.0

        # PostOutbreakStartAlerts: only second outbreak detected
        matched_filtered = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
        )
        @test matched_filtered.outbreak_indices_with_alerts == [2]
        @test matched_filtered.alert_indices_per_outbreak == [[2]]
        @test calculate_sensitivity(matched_filtered) == 0.5  # 1 out of 2 outbreaks
        @test calculate_ppv(matched_filtered) == 0.5  # 1 matched out of 2 total alerts
    end

    @testset "AlertFilteringStrategy: PostOutbreakStartAlerts with no valid alerts" begin
        outbreaks = OutbreakThresholds(
            lower_bounds = [20, 60],
            upper_bounds = [40, 80],
            duration = [21, 21],
            num_infections_during_bounds = [150, 200]
        )

        # All alerts start before their respective outbreaks
        alerts = Thresholds(
            lower_bounds = [15, 55],
            upper_bounds = [25, 65],
            duration = [11, 11]
        )

        # PostOutbreakStartAlerts: no outbreaks detected
        matched_filtered = match_outbreak_detection_bounds(
            outbreaks, alerts, AlertFilteringStrategy(PostOutbreakStartAlerts())
        )
        @test isempty(matched_filtered.outbreak_indices_with_alerts)
        @test calculate_sensitivity(matched_filtered) == 0.0
        @test calculate_ppv(matched_filtered) == 0.0
    end
end
