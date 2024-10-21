using OutbreakDetectionUtils, StatsBase

@testset "diag-testing-functions.jl" begin
    @testset "Moving average" begin
        daily_testpositives = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        movingavg_testpositives = calculate_movingavg(
            daily_testpositives,
            5,
        )

        @test isequal(
            movingavg_testpositives,
            [
                StatsBase.mean([1]),
                StatsBase.mean([1, 2]),
                StatsBase.mean([1, 2, 3]),
                StatsBase.mean([1, 2, 3, 4]),
                StatsBase.mean([1, 2, 3, 4, 5]),
                StatsBase.mean([2, 3, 4, 5, 6]),
                StatsBase.mean([3, 4, 5, 6, 7]),
                StatsBase.mean([4, 5, 6, 7, 8]),
                StatsBase.mean([5, 6, 7, 8, 9]),
                StatsBase.mean([6, 7, 8, 9, 10]),
            ],
        )

        @test length(movingavg_testpositives) == length(daily_testpositives)

        @test begin
            daily_testpositives = zeros(Int64, 10, 2)
            daily_testpositives[:, 1] .= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            calculate_movingavg!(
                @view(daily_testpositives[:, 2]),
                @view(daily_testpositives[:, 1]),
                5,
            )

            isequal(
                daily_testpositives[:, 2],
                Int64.(
                    round.([
                        StatsBase.mean([1]),
                        StatsBase.mean([1, 2]),
                        StatsBase.mean([1, 2, 3]),
                        StatsBase.mean([1, 2, 3, 4]),
                        StatsBase.mean([1, 2, 3, 4, 5]),
                        StatsBase.mean([2, 3, 4, 5, 6]),
                        StatsBase.mean([3, 4, 5, 6, 7]),
                        StatsBase.mean([4, 5, 6, 7, 8]),
                        StatsBase.mean([5, 6, 7, 8, 9]),
                        StatsBase.mean([6, 7, 8, 9, 10]),
                    ])
                ),
            )
        end
    end

    @testset "Detect outbreak" begin
        incvec = [1, 3, 10, 15, 20, 3, 1]
        avgvec = calculate_movingavg(incvec, 3)
        threshold = 5

        @test begin
            outbreakvec = detectoutbreak(incvec, threshold)
            outbreakvec == [0, 0, 1, 1, 1, 0, 0]
        end

        @test begin
            outbreakvec = detectoutbreak(avgvec, threshold)
            outbreakvec == [0, 0, 0, 1, 1, 1, 1]
        end

        @test begin
            outbreakvec = detectoutbreak(incvec, avgvec, threshold)
            outbreakvec == [0, 0, 1, 1, 1, 1, 1]
        end
    end

    @testset "Matched bounds" begin
        @test begin
            outbreakbounds = [
                10 60 51 600
                100 180 81 700
                300 340 41 800
                380 410 31 900
                500 540 41 1000
                600 660 61 1100
            ]
            detectionbounds = [
                5 15 11
                17 40 24
                50 80 31
                90 105 16
                110 160 51
                390 420 31
                495 550 56
                590 595 6
            ]

            isequal(
                OutbreakDetectionUtils.match_outbreak_detection_bounds(
                    outbreakbounds, detectionbounds
                ),
                (
                    [
                        10 60 5 15 600
                        10 60 17 40 600
                        10 60 50 80 600
                        100 180 90 105 700
                        100 180 110 160 700
                        380 410 390 420 900
                        500 540 495 550 1000
                    ],
                    [51, 81, 41, 31, 41, 61],
                    [11, 24, 31, 16, 51, 31, 56, 6],
                    [600, 700, 800, 900, 1000, 1100],
                    [3, 2, 0, 1, 1, 0],
                ),
            )
        end

        @test begin
            matched_bounds = [
                10 60 5 15 600
                10 60 17 40 600
                10 60 50 80 600
                100 180 90 105 700
                100 180 110 160 700
                380 410 390 420 900
                500 540 495 550 1000
            ]
            filtered_bounds = [
                10 60 5 15 600
                100 180 90 105 700
                380 410 390 420 900
                500 540 495 550 1000
            ]

            isequal(
                filter_first_matched_bounds(matched_bounds), filtered_bounds
            )
        end
    end

    @testset "Cases before and after alert" begin
        @test begin
            incarr = [
                repeat([1], 9)..., # First alert triggered here
                repeat([12], 51)..., # Outbreak here
                repeat([1], 19)...,
                repeat([15], 11)..., # NOT outbreak here (i = 80 to 90)
                repeat([1], 9)...,
                repeat([8], 81)..., # Outbreak here
                repeat([1], 199)...,
                repeat([30], 31)..., # Outbreak here; third alert triggered here
                repeat([1], 89)...,  # Fourth alert triggered here
                repeat([25], 41)..., # Outbreak here
            ]

            matched_bounds = [
                10 60 5 15 612
                100 180 90 105 648
                380 410 390 420 930
                500 540 495 550 1025
            ]

            delay_vec = [-5, -10, 10, -5]
            isequal(
                calculate_cases_before_after_alert(
                    incarr, matched_bounds, delay_vec
                ),
                (
                    [0, 0, 300, 0],
                    [0.0, 0.0, 300 / 930, 0.0],
                    [612, 648, 630, 1025],
                    [1.0, 1.0, 630 / 930, 1.0],
                ),
            )
        end
    end

    @testset "Outbreak detection characteristics" begin
        outbreakbounds = [
            10 60 51 12*51
            100 180 81 8*81
            380 410 31 30*31
            500 540 41 25*41
            600 660 61 1100
        ]

        detectionbounds = [
            5 15 11
            17 40 24
            50 80 31
            90 105 16
            110 160 51
            390 420 31
            495 550 56
            590 595 6
        ]

        outbreak_detection_chars = calculate_outbreak_detection_characteristics(
            outbreakbounds, detectionbounds
        )

        accuracy = mean([4 / 5, 7 / 8])
        matched_bounds = [
            10 60 5 15 12*51
            10 60 17 40 12*51
            10 60 50 80 12*51
            100 180 90 105 8*81
            100 180 110 160 8*81
            380 410 390 420 30*31
            500 540 495 550 25*41
        ]
        noutbreaks = 5
        nalerts = 8
        outbreak_duration_vec = [51, 81, 31, 41, 61]
        alert_duration_vec = [11, 24, 31, 16, 51, 31, 56, 6]
        detected_outbreak_size = [
            12 * 51, 8 * 81, 30 * 31, 25 * 41
        ]
        missed_outbreak_size = [1100]
        n_true_outbreaks_detected = 4
        n_missed_outbreaks = 1
        n_correct_alerts = 7
        n_false_alerts = 1
        alertsperoutbreak = [3, 2, 1, 1, 0]
        periodsumvec = [
            12 * 51,
            8 * 81,
            30 * 31,
            25 * 41,
            1100,
        ]
        perc_true_outbreaks_detected = 4 / 5
        perc_true_outbreaks_missed = 1 / 5
        falsealert_trueoutbreak_prop = 1 / 5
        correctalert_trueoutbreak_prop = 7 / 5
        trueoutbreak_alerts_prop = 5 / 8
        outbreaksmissed_alerts_prop = 1 / 8
        perc_alerts_false = 1 / 8
        perc_alerts_correct = 7 / 8

        expected_outbreak_detection_chars = (;
            accuracy,
            matched_bounds,
            noutbreaks,
            nalerts,
            outbreak_duration_vec,
            alert_duration_vec,
            detected_outbreak_size,
            missed_outbreak_size,
            n_true_outbreaks_detected,
            n_missed_outbreaks,
            n_correct_alerts,
            n_false_alerts,
            alertsperoutbreak,
            periodsumvec,
            perc_true_outbreaks_detected,
            perc_true_outbreaks_missed,
            falsealert_trueoutbreak_prop,
            correctalert_trueoutbreak_prop,
            trueoutbreak_alerts_prop,
            outbreaksmissed_alerts_prop,
            perc_alerts_false,
            perc_alerts_correct,
        )

        for variable in propertynames(outbreak_detection_chars)
            @testset "$(string(variable))" begin
                actual_value = getproperty(outbreak_detection_chars, variable)
                expected_value = getproperty(
                    expected_outbreak_detection_chars, variable
                )
                @test isequal(actual_value, expected_value)
            end
        end

        @test isequal(
            outbreak_detection_chars,
            expected_outbreak_detection_chars,
        )

        @test begin
            infectious_tested_vec = [
                repeat([1], 9)...,
                repeat([12], 51)..., # Outbreak here (i = 10 to 60)
                repeat([1], 19)...,
                repeat([15], 11)..., # NOT outbreak here (i = 80 to 90)
                repeat([1], 9)...,
                repeat([8], 81)..., # Outbreak here (i = 100 to 180)
                repeat([1], 199)...,
                repeat([30], 31)..., # Outbreak here (i = 380 to 410)
                repeat([1], 89)...,
                repeat([25], 41)..., # Outbreak here (i = 500 to 540)
            ]

            noise_tested_vec = [
                repeat([5], 9)...,
                repeat([4], 51)..., # Outbreak here
                repeat([6], 39)...,
                repeat([5], 81)..., # Outbreak here
                repeat([4], 199)...,
                repeat([6], 31)..., # Outbreak here
                repeat([5], 89)...,
                repeat([4], 41)..., # Outbreak here
            ]

            outbreakbounds = [
                10 60 612
                100 180 648
                380 410 930
                500 540 1025
            ]

            isequal(
                calculate_n_outbreak_tests(
                    infectious_tested_vec, noise_tested_vec, outbreakbounds
                ),
                (12 * 51 + 8 * 81 + 30 * 31 + 25 * 41) +
                (4 * 51 + 5 * 81 + 6 * 31 + 4 * 41),
            )
        end
    end
end
