using Test
using OutbreakDetection
using OutbreakDetectionCore

@testset "sort_test_specifications" begin
    @testset "sorts by test_result_lag descending, then by (sensitivity, specificity) ascending" begin
        # Create test specifications with different combinations
        specs = [
            IndividualTestSpecification(sensitivity = 0.9, specificity = 0.95, test_result_lag = 1),
            IndividualTestSpecification(sensitivity = 0.8, specificity = 0.9, test_result_lag = 2),
            IndividualTestSpecification(sensitivity = 0.85, specificity = 0.92, test_result_lag = 2),
            IndividualTestSpecification(sensitivity = 0.95, specificity = 0.98, test_result_lag = 1),
        ]

        sorted_specs = sort_test_specifications(specs)

        # Expected order:
        # 1. lag=2, sens=0.8, spec=0.90
        # 2. lag=2, sens=0.85, spec=0.92
        # 3. lag=1, sens=0.9, spec=0.95
        # 4. lag=1, sens=0.95, spec=0.98

        @test sorted_specs[1].test_result_lag == 2
        @test sorted_specs[1].sensitivity == 0.8
        @test sorted_specs[1].specificity == 0.9

        @test sorted_specs[2].test_result_lag == 2
        @test sorted_specs[2].sensitivity == 0.85
        @test sorted_specs[2].specificity == 0.92

        @test sorted_specs[3].test_result_lag == 1
        @test sorted_specs[3].sensitivity == 0.9
        @test sorted_specs[3].specificity == 0.95

        @test sorted_specs[4].test_result_lag == 1
        @test sorted_specs[4].sensitivity == 0.95
        @test sorted_specs[4].specificity == 0.98
    end

    @testset "handles single element" begin
        specs = [IndividualTestSpecification(sensitivity = 0.9, specificity = 0.95, test_result_lag = 1)]
        sorted_specs = sort_test_specifications(specs)

        @test length(sorted_specs) == 1
        @test sorted_specs[1] == specs[1]
    end

    @testset "handles empty collection" begin
        specs = IndividualTestSpecification[]
        sorted_specs = sort_test_specifications(specs)

        @test isempty(sorted_specs)
    end

    @testset "sorts by sensitivity when specificity is equal" begin
        specs = [
            IndividualTestSpecification(sensitivity = 0.9, specificity = 0.95, test_result_lag = 1),
            IndividualTestSpecification(sensitivity = 0.8, specificity = 0.95, test_result_lag = 1),
        ]

        sorted_specs = sort_test_specifications(specs)

        # Should be sorted by sensitivity ascending when lag and specificity are equal
        @test sorted_specs[1].sensitivity == 0.8
        @test sorted_specs[2].sensitivity == 0.9
    end

    @testset "sorts by specificity when sensitivity is equal" begin
        specs = [
            IndividualTestSpecification(sensitivity = 0.9, specificity = 0.98, test_result_lag = 1),
            IndividualTestSpecification(sensitivity = 0.9, specificity = 0.95, test_result_lag = 1),
        ]

        sorted_specs = sort_test_specifications(specs)

        # Should be sorted by specificity ascending when lag and sensitivity are equal
        @test sorted_specs[1].specificity == 0.95
        @test sorted_specs[2].specificity == 0.98
    end
end
