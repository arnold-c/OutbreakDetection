using DrWatson, Test
@quickactivate "OutbreakDetection"

using OutbreakDetection

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "OutbreakDetection tests" begin
    @test 1 == 1
end

@testset "Detection tests" begin
    @testset "Delay detection" begin
        @test begin
            outbreakbounds = [
                2 4 500
                10 60 600
                100 180 700
                300 340 800
                380 410 900
                500 540 1000
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
                    unmatched_bounds = [
                        2 4 500
                        300 340 800
                    ],
                    noutbreaks = 6,
                    ndetectoutbreaks = 8,
                    n_true_outbreaks_detected = 4,
                    n_missed_outbreaks = 2,
                    n_correct_alerts = 7,
                    n_false_alerts = 1,
                    alertsperoutbreak = [0, 3, 2, 0, 1, 1],
                    perc_true_outbreaks_detected = 4 / 6,
                    perc_true_outbreaks_missed = 2 / 6,
                    falsealert_trueoutbreak_prop = 1 / 6,
                    correctalert_trueoutbreak_prop = 7 / 6,
                    trueoutbreak_alerts_prop = 6 / 8,
                    outbreaksmissed_alerts_prop = 2 / 8,
                    perc_alerts_false = 1 / 8,
                    perc_alerts_correct = 7 / 8,
                    detectiondelays = [-5, -10, 10, -5],
                ),
            )
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
