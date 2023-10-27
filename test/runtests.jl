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
            calculate_outbreak_detection_delay(
                outbreakbounds, detectionbounds
            ) == (
                [Inf, -5.0, -10.0, Inf, 10.0, -5.0],
                [2.0 4.0 Inf Inf
                    10.0 60.0 5.0 15.0
                    100.0 180.0 90.0 105.0
                    300.0 340.0 Inf Inf
                    380.0 410.0 390.0 420.0
                    500.0 540.0 495.0 550.0],
                2.0,
            )
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
