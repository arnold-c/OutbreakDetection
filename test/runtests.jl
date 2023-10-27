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
                2 4
                10 60
                100 180
                300 340
                380 410
                500 540
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
                    [-5, -10, 10, -5],
                    [
                        10 60 5 15
                        10 60 17 40
                        10 60 50 80
                        100 180 90 105
                        100 180 110 160
                        380 410 390 420
                        500 540 495 550
                    ],
                    2,
                    1,
                    [0, 3, 2, 0, 1, 1],
                ),
            )
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
