@testset "detection-metric-functions.jl" begin
    using OutbreakDetectionCore
    using NaNMath

    # Helper function to create test data
    function create_test_thresholds()
        outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
            lower_bounds = [10, 30, 50],
            upper_bounds = [20, 40, 60],
            duration = [11, 11, 11],
            num_infections_during_bounds = [100, 150, 200]
        )

        alert_thresholds = OutbreakDetectionCore.Thresholds(
            lower_bounds = [15, 55],
            upper_bounds = [25, 65],
            duration = [11, 11]
        )

        return outbreak_thresholds, alert_thresholds
    end

    @testset "calculate_detection_delay" begin
        @testset "basic detection delay calculation" begin
            outbreak_thresholds, alert_thresholds = create_test_thresholds()

            # Outbreak 1 (idx 1) matched to alert 1: delay = 15 - 10 = 5
            # Outbreak 3 (idx 3) matched to alert 2: delay = 55 - 50 = 5
            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1, 3],
                alert_indices_per_outbreak = [[1], [2]],
                n_matched_outbreaks = 2,
                n_matched_alerts = 2,
                n_outbreaks = 3,
                n_alerts = 2
            )

            delays = OutbreakDetectionCore.calculate_detection_delay(
                matched, outbreak_thresholds, alert_thresholds
            )

            @test delays == [5, 5]
        end

        @testset "varying delays" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30],
                upper_bounds = [20, 40],
                duration = [11, 11],
                num_infections_during_bounds = [100, 150]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [12, 35],
                upper_bounds = [22, 45],
                duration = [11, 11]
            )

            # Outbreak 1: delay = 12 - 10 = 2
            # Outbreak 2: delay = 35 - 30 = 5
            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1, 2],
                alert_indices_per_outbreak = [[1], [2]],
                n_matched_outbreaks = 2,
                n_matched_alerts = 2,
                n_outbreaks = 2,
                n_alerts = 2
            )

            delays = OutbreakDetectionCore.calculate_detection_delay(
                matched, outbreak_thresholds, alert_thresholds
            )

            @test delays == [2, 5]
        end

        @testset "no matched outbreaks" begin
            outbreak_thresholds, alert_thresholds = create_test_thresholds()

            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = Int64[],
                alert_indices_per_outbreak = Vector{Int64}[],
                n_matched_outbreaks = 0,
                n_matched_alerts = 0,
                n_outbreaks = 3,
                n_alerts = 2
            )

            delays = OutbreakDetectionCore.calculate_detection_delay(
                matched, outbreak_thresholds, alert_thresholds
            )

            @test isempty(delays)
        end

        @testset "alert before outbreak start (negative delay)" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [20],
                upper_bounds = [30],
                duration = [11],
                num_infections_during_bounds = [100]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [15],
                upper_bounds = [25],
                duration = [11]
            )

            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1],
                alert_indices_per_outbreak = [[1]],
                n_matched_outbreaks = 1,
                n_matched_alerts = 1,
                n_outbreaks = 1,
                n_alerts = 1
            )

            delays = OutbreakDetectionCore.calculate_detection_delay(
                matched, outbreak_thresholds, alert_thresholds
            )

            @test delays == [-5]  # Alert started 5 days before outbreak
        end

        @testset "multiple alerts per outbreak (uses first)" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10],
                upper_bounds = [30],
                duration = [21],
                num_infections_during_bounds = [200]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [15, 20, 25],
                upper_bounds = [18, 23, 28],
                duration = [4, 4, 4]
            )

            # Outbreak matched to alerts 1, 2, 3 - should use first (alert 1)
            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1],
                alert_indices_per_outbreak = [[1, 2, 3]],
                n_matched_outbreaks = 1,
                n_matched_alerts = 3,
                n_outbreaks = 1,
                n_alerts = 3
            )

            delays = OutbreakDetectionCore.calculate_detection_delay(
                matched, outbreak_thresholds, alert_thresholds
            )

            @test delays == [5]  # Uses first alert: 15 - 10 = 5
        end
    end

    @testset "calculate_unavoidable_cases" begin
        @testset "matched outbreaks only" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30],
                upper_bounds = [20, 40],
                duration = [11, 11],
                num_infections_during_bounds = [100, 150]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [15, 35],
                upper_bounds = [25, 45],
                duration = [11, 11]
            )

            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1, 2],
                alert_indices_per_outbreak = [[1], [2]],
                n_matched_outbreaks = 2,
                n_matched_alerts = 2,
                n_outbreaks = 2,
                n_alerts = 2
            )

            # Create incidence vector
            incidence_vec = zeros(Int64, 50)
            incidence_vec[10:14] .= 10  # 5 days * 10 = 50 cases before alert 1
            incidence_vec[30:34] .= 20  # 5 days * 20 = 100 cases before alert 2

            unavoidable = OutbreakDetectionCore.calculate_unavoidable_cases(
                matched, outbreak_thresholds, alert_thresholds, incidence_vec
            )

            @test unavoidable == [50, 100]
        end

        @testset "includes unmatched outbreaks" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30, 50],
                upper_bounds = [20, 40, 60],
                duration = [11, 11, 11],
                num_infections_during_bounds = [100, 150, 200]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [15],
                upper_bounds = [25],
                duration = [11]
            )

            # Only outbreak 1 is matched; outbreaks 2 and 3 are unmatched
            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1],
                alert_indices_per_outbreak = [[1]],
                n_matched_outbreaks = 1,
                n_matched_alerts = 1,
                n_outbreaks = 3,
                n_alerts = 1
            )

            # Create incidence vector
            incidence_vec = zeros(Int64, 70)
            incidence_vec[10:14] .= 10  # Outbreak 1: 50 cases before alert
            incidence_vec[15:20] .= 5   # Outbreak 1: after alert (not counted)
            incidence_vec[30:40] .= 15  # Outbreak 2: 11 * 15 = 165 cases (all unavoidable)
            incidence_vec[50:60] .= 20  # Outbreak 3: 11 * 20 = 220 cases (all unavoidable)

            unavoidable = OutbreakDetectionCore.calculate_unavoidable_cases(
                matched, outbreak_thresholds, alert_thresholds, incidence_vec
            )

            @test unavoidable == [50, 165, 220]
        end

        @testset "alert at outbreak start (zero unavoidable)" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10],
                upper_bounds = [20],
                duration = [11],
                num_infections_during_bounds = [100]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [10],
                upper_bounds = [20],
                duration = [11]
            )

            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = [1],
                alert_indices_per_outbreak = [[1]],
                n_matched_outbreaks = 1,
                n_matched_alerts = 1,
                n_outbreaks = 1,
                n_alerts = 1
            )

            incidence_vec = zeros(Int64, 30)
            incidence_vec[10:20] .= 10

            unavoidable = OutbreakDetectionCore.calculate_unavoidable_cases(
                matched, outbreak_thresholds, alert_thresholds, incidence_vec
            )

            @test unavoidable == [0]
        end

        @testset "no outbreaks" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = Int64[],
                upper_bounds = Int64[],
                duration = Int64[],
                num_infections_during_bounds = Int64[]
            )

            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [10],
                upper_bounds = [20],
                duration = [11]
            )

            matched = OutbreakDetectionCore.MatchedThresholds(
                outbreak_indices_with_alerts = Int64[],
                alert_indices_per_outbreak = Vector{Int64}[],
                n_matched_outbreaks = 0,
                n_matched_alerts = 0,
                n_outbreaks = 0,
                n_alerts = 1
            )

            incidence_vec = zeros(Int64, 30)

            unavoidable = OutbreakDetectionCore.calculate_unavoidable_cases(
                matched, outbreak_thresholds, alert_thresholds, incidence_vec
            )

            @test isempty(unavoidable)
        end
    end

    @testset "calculate_alert_duration" begin
        @testset "uniform durations" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [5, 15, 25],
                upper_bounds = [10, 20, 30],
                duration = [6, 6, 6]
            )

            durations = OutbreakDetectionCore.calculate_alert_duration(
                alert_thresholds
            )

            @test durations == [6, 6, 6]
        end

        @testset "varying durations" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [5, 15, 25],
                upper_bounds = [10, 22, 35],
                duration = [6, 8, 11]
            )

            durations = OutbreakDetectionCore.calculate_alert_duration(
                alert_thresholds
            )

            @test durations == [6, 8, 11]
        end

        @testset "no alerts" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = Int64[],
                upper_bounds = Int64[],
                duration = Int64[]
            )

            durations = OutbreakDetectionCore.calculate_alert_duration(
                alert_thresholds
            )

            @test isempty(durations)
        end

        @testset "single alert" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [10],
                upper_bounds = [25],
                duration = [16]
            )

            durations = OutbreakDetectionCore.calculate_alert_duration(
                alert_thresholds
            )

            @test durations == [16]
        end
    end

    @testset "calculate_outbreak_duration" begin
        @testset "uniform durations" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30, 50],
                upper_bounds = [20, 40, 60],
                duration = [11, 11, 11],
                num_infections_during_bounds = [100, 150, 200]
            )

            durations = OutbreakDetectionCore.calculate_outbreak_duration(
                outbreak_thresholds
            )

            @test durations == [11, 11, 11]
        end

        @testset "varying durations" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30, 50],
                upper_bounds = [20, 45, 65],
                duration = [11, 16, 16],
                num_infections_during_bounds = [100, 200, 150]
            )

            durations = OutbreakDetectionCore.calculate_outbreak_duration(
                outbreak_thresholds
            )

            @test durations == [11, 16, 16]
        end

        @testset "no outbreaks" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = Int64[],
                upper_bounds = Int64[],
                duration = Int64[],
                num_infections_during_bounds = Int64[]
            )

            durations = OutbreakDetectionCore.calculate_outbreak_duration(
                outbreak_thresholds
            )

            @test isempty(durations)
        end

        @testset "single outbreak" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10],
                upper_bounds = [30],
                duration = [21],
                num_infections_during_bounds = [250]
            )

            durations = OutbreakDetectionCore.calculate_outbreak_duration(
                outbreak_thresholds
            )

            @test durations == [21]
        end
    end

    @testset "calculate_proportion_timeseries_in_alert" begin
        @testset "basic proportion calculation" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [5, 15],
                upper_bounds = [10, 20],
                duration = [6, 6]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_alert(
                alert_thresholds, 100
            )

            @test proportion ≈ 0.12  # (6 + 6) / 100
        end

        @testset "varying alert durations" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [10, 30, 50],
                upper_bounds = [20, 35, 60],
                duration = [11, 6, 11]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_alert(
                alert_thresholds, 200
            )

            @test proportion ≈ 0.14  # (11 + 6 + 11) / 200
        end

        @testset "no alerts" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = Int64[],
                upper_bounds = Int64[],
                duration = Int64[]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_alert(
                alert_thresholds, 100
            )

            @test proportion ≈ 0.0
        end

        @testset "zero timeseries length" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [10],
                upper_bounds = [20],
                duration = [11]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_alert(
                alert_thresholds, 0
            )

            @test isnan(proportion)
        end

        @testset "entire timeseries in alert" begin
            alert_thresholds = OutbreakDetectionCore.Thresholds(
                lower_bounds = [1],
                upper_bounds = [100],
                duration = [100]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_alert(
                alert_thresholds, 100
            )

            @test proportion ≈ 1.0
        end
    end

    @testset "calculate_proportion_timeseries_in_outbreak" begin
        @testset "basic proportion calculation" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30],
                upper_bounds = [20, 45],
                duration = [11, 16],
                num_infections_during_bounds = [100, 200]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_outbreak(
                outbreak_thresholds, 100
            )

            @test proportion ≈ 0.27  # (11 + 16) / 100
        end

        @testset "varying outbreak durations" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10, 30, 60],
                upper_bounds = [25, 40, 75],
                duration = [16, 11, 16],
                num_infections_during_bounds = [150, 100, 200]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_outbreak(
                outbreak_thresholds, 200
            )

            @test proportion ≈ 0.215  # (16 + 11 + 16) / 200
        end

        @testset "no outbreaks" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = Int64[],
                upper_bounds = Int64[],
                duration = Int64[],
                num_infections_during_bounds = Int64[]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_outbreak(
                outbreak_thresholds, 100
            )

            @test proportion ≈ 0.0
        end

        @testset "zero timeseries length" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [10],
                upper_bounds = [30],
                duration = [21],
                num_infections_during_bounds = [200]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_outbreak(
                outbreak_thresholds, 0
            )

            @test isnan(proportion)
        end

        @testset "entire timeseries in outbreak" begin
            outbreak_thresholds = OutbreakDetectionCore.OutbreakThresholds(
                lower_bounds = [1],
                upper_bounds = [100],
                duration = [100],
                num_infections_during_bounds = [1000]
            )

            proportion = OutbreakDetectionCore.calculate_proportion_timeseries_in_outbreak(
                outbreak_thresholds, 100
            )

            @test proportion ≈ 1.0
        end
    end
end
