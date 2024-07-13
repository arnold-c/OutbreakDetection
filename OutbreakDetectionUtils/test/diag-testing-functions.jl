@testset "diag-testing-functions.jl" begin
    @testset "Moving average" begin
        using OutbreakDetectionUtils, Statistics

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
        using OutbreakDetectionUtils

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
        using OutbreakDetectionUtils

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

        @test isequal(
            filter_first_matched_bounds(matched_bounds), filtered_bounds
        )
    end

    @testset "Cases before and after alert" begin
        using OutbreakDetectionUtils
        @test begin
            incarr = [
                repeat([1], 9)...,
                repeat([12], 51)...,
                repeat([1], 39)...,
                repeat([8], 81)...,
                repeat([1], 199)...,
                repeat([30], 31)...,
                repeat([1], 89)...,
                repeat([25], 41)...,
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
                    [612, 648, 630, 1025],
                    [1.0, 1.0, 630 / 930, 1.0],
                ),
            )
        end
    end

    @testset "Outbreak detection characteristics" begin
        using OutbreakDetectionUtils
        @test begin
            outbreakbounds = [
                2 4 500 100
                10 60 600 10
                100 180 700 2000
                300 340 800 20
                380 410 900 400
                500 540 1000 5000
            ]
            detectionbounds = [
                5 15
                17 40
                50 80
                90 105
                110 160
                390 420
                495 550
                590 595
            ]

            isequal(
                calculate_outbreak_detection_characteristics(
                    outbreakbounds, detectionbounds
                ),
                (
                    matched_bounds = [
                        10 60 5 15 600
                        10 60 17 40 600
                        10 60 50 80 600
                        100 180 90 105 700
                        100 180 110 160 700
                        380 410 390 420 900
                        500 540 495 550 1000
                    ],
                    noutbreaks = 6,
                    nalerts = 8,
                    detected_outbreak_size = [600, 700, 900, 1000],
                    missed_outbreak_size = [500, 800],
                    n_true_outbreaks_detected = 4,
                    n_missed_outbreaks = 2,
                    n_correct_alerts = 7,
                    n_false_alerts = 1,
                    alertsperoutbreak = [0, 3, 2, 0, 1, 1],
                    periodsumvec = [500, 600, 700, 800, 900, 1000],
                    perc_true_outbreaks_detected = 4 / 6,
                    perc_true_outbreaks_missed = 2 / 6,
                    falsealert_trueoutbreak_prop = 1 / 6,
                    correctalert_trueoutbreak_prop = 7 / 6,
                    trueoutbreak_alerts_prop = 6 / 8,
                    outbreaksmissed_alerts_prop = 2 / 8,
                    perc_alerts_false = 1 / 8,
                    perc_alerts_correct = 7 / 8,
                    # detectiondelays = [-5, -10, 10, -5],
                ),
            )
        end

        @test begin
            infectious_tested_vec = [
                repeat([1], 9)...,
                repeat([12], 51)..., # Outbreak here
                repeat([1], 39)...,
                repeat([8], 81)..., # Outbreak here
                repeat([1], 199)...,
                repeat([30], 31)..., # Outbreak here
                repeat([1], 89)...,
                repeat([25], 41)..., # Outbreak here
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
                calculate_n_outbreak_tests(infectious_tested_vec, noise_tested_vec, outbreakbounds),
                (12*51 + 8*81 + 30*31 + 25*41) + (4*51 + 5*81 + 6*31 + 4*41),
            )
        end
    end
end
