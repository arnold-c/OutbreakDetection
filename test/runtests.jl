using DrWatson, Test
@quickactivate "OutbreakDetection"

using OutbreakDetection

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Detection tests" begin
    @testset "Moving average" begin
        @test begin
            daily_testpositives = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            movingavg_testpositives = calculate_movingavg(
                daily_testpositives,
                5
            )

            isequal(
                movingavg_testpositives,
                [
                    mean([1]),
                    mean([1, 2]),
                    mean([1, 2, 3]),
                    mean([1, 2, 3, 4]),
                    mean([1, 2, 3, 4, 5]),
                    mean([2, 3, 4, 5, 6]),
                    mean([3, 4, 5, 6, 7]),
                    mean([4, 5, 6, 7, 8]),
                    mean([5, 6, 7, 8, 9]),
                    mean([6, 7, 8, 9, 10]),
                ],
            )
        end
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
                        mean([1]),
                        mean([1, 2]),
                        mean([1, 2, 3]),
                        mean([1, 2, 3, 4]),
                        mean([1, 2, 3, 4, 5]),
                        mean([2, 3, 4, 5, 6]),
                        mean([3, 4, 5, 6, 7]),
                        mean([4, 5, 6, 7, 8]),
                        mean([5, 6, 7, 8, 9]),
                        mean([6, 7, 8, 9, 10]),
                    ])
                ),
            )
        end
    end
    @testset "Delay detection" begin
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
                calculate_cases_after_alert(incarr, matched_bounds, delay_vec),
                (
                    [612, 648, 630, 1025],
                    [1.0, 1.0, 630 / 930, 1.0],
                ),
            )
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
