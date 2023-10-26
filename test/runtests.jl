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
            outbreakbounds = [2 4; 10 60; 100 180; 300 340; 380 400]
            detectionbounds = [5 15; 12 40; 50 80; 90 105; 110 160; 390 420]
            calculate_outbreak_detection_delay(
                outbreakbounds, detectionbounds
            ) ==
            (
                [-10, 2, 10, -10, 10],
                [2 4 -10 -10
                    10 60 12 40
                    100 180 110 160
                    300 340 -10 -10
                    380 400 390 420])
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits = 3), " minutes")
