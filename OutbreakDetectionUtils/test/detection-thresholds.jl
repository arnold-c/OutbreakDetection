@testset "detection-thresholds.jl" begin
    using OutbreakDetectionUtils, StatsBase

    outbreakthreshold = 5
    minoutbreakdur = 30
    minoutbreaksize = 500

    inc_vec = [
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
    inc_rle = StatsBase.rle(inc_vec .> outbreakthreshold)
    outbreak_thresholds = calculate_outbreak_thresholds(
        inc_rle; ncols = 5
    )

    @testset "Outbreak thresholds" begin
        @test begin
            isequal(
                outbreak_thresholds,
                [
                    10 60 0 0 0
                    80 90 0 0 0
                    100 180 0 0 0
                    380 410 0 0 0
                    500 540 0 0 0
                ],
            )
        end
    end

    @testset "Outbreak duration" begin
        @test isequal(
            OutbreakDetectionUtils.calculate_outbreak_duration(
                outbreak_thresholds[1, :]
            ),
            51,
        )

        @test isequal(
            OutbreakDetectionUtils.calculate_outbreak_duration(
                outbreak_thresholds[2, 1],
                outbreak_thresholds[2, 2],
            ),
            11,
        )

        OutbreakDetectionUtils.calculate_outbreak_duration!(
            outbreak_thresholds
        )
        @test isequal(
            outbreak_thresholds[:, 3],
            [51, 11, 81, 31, 41],
        )
    end

    @testset "Outbreak size" begin
        @test isequal(
            OutbreakDetectionUtils.calculate_outbreak_size(
                inc_vec,
                outbreak_thresholds[1, 1],
                outbreak_thresholds[1, 2],
            ),
            12 * 51,
        )
    end

    @testset "Classifying Outbreaks" begin
        @testset "Single Outbreak" begin
            @test isequal(
                OutbreakDetectionUtils.classify_outbreak(
                    outbreak_thresholds[1, 1],
                    outbreak_thresholds[1, 2],
                    minoutbreakdur,
                    OutbreakDetectionUtils.calculate_outbreak_size(
                        inc_vec,
                        outbreak_thresholds[1, 1],
                        outbreak_thresholds[1, 2],
                    ),
                    minoutbreaksize,
                ),
                12 * 51 > minoutbreaksize,
            )
            @test isequal(
                OutbreakDetectionUtils.classify_outbreak(
                    outbreak_thresholds[1, 3],
                    minoutbreakdur,
                    12 * 51,
                    minoutbreaksize,
                ),
                12 * 51 > minoutbreaksize,
            )
        end

        @testset "Multiple Outbreaks" begin
            outbreakstatus_vec = zeros(Float64, length(inc_vec))
            classify_all_outbreaks!(
                outbreakstatus_vec,
                outbreak_thresholds,
                inc_vec,
                minoutbreakdur,
                minoutbreaksize,
            )

            @test isequal(
                outbreakstatus_vec,
                [
                    repeat([0], 9)...,
                    repeat([1], 51)..., # Outbreak here (i = 10 to 60)
                    repeat([0], 19)...,
                    repeat([0], 11)..., # NOT outbreak here (i = 80 to 90)
                    repeat([0], 9)...,
                    repeat([1], 81)..., # Outbreak here (i = 100 to 180)
                    repeat([0], 199)...,
                    repeat([1], 31)..., # Outbreak here (i = 380 to 410)
                    repeat([0], 89)...,
                    repeat([1], 41)..., # Outbreak here (i = 500 to 540)
                ],
            )

            @test isequal(
                outbreak_thresholds,
                [
                    10 60 51 12*51 1
                    80 90 11 15*11 0
                    100 180 81 8*81 1
                    380 410 31 30*31 1
                    500 540 41 25*41 1
                ],
            )
        end
    end

    @testset "Filtering only outbreaks thresholds" begin
        @test isequal(
            OutbreakDetectionUtils.filter_only_outbreaks(
                outbreak_thresholds
            ),
            [
                10 60 51 12*51 1
                100 180 81 8*81 1
                380 410 31 30*31 1
                500 540 41 25*41 1
            ],
        )
    end
end
