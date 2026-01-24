@testset "struct-version-hash.jl" begin
    using OutbreakDetectionCore
    using Test

    # Define test structs for testing hash behavior
    struct TestStruct1
        field1::Int
        field2::String
    end

    struct TestStruct2
        field1::Int
        field2::String
    end

    struct TestStruct3
        field1::Int
        field2::String
        field3::Float64
    end

    struct TestStruct4
        field2::String
        field1::Int
    end

    struct TestStruct5
        field1::Float64
        field2::String
    end

    @testset "struct_definition_hash" begin
        @testset "returns UInt64 hash" begin
            hash_val = OutbreakDetectionCore.struct_definition_hash(TestStruct1)

            @test hash_val isa UInt64
        end

        @testset "generates consistent hashes for same struct" begin
            hash1 = OutbreakDetectionCore.struct_definition_hash(TestStruct1)
            hash2 = OutbreakDetectionCore.struct_definition_hash(TestStruct1)

            @test hash1 == hash2
        end

        @testset "generates different hashes for different structs with same fields" begin
            # TestStruct1 and TestStruct2 have identical fields but different names
            # The hash should include the struct name, so they must be different
            hash1 = OutbreakDetectionCore.struct_definition_hash(TestStruct1)
            hash2 = OutbreakDetectionCore.struct_definition_hash(TestStruct2)

            @test hash1 isa UInt64
            @test hash2 isa UInt64
            @test hash1 != hash2  # Must be different due to different struct names
        end

        @testset "generates different hashes when fields are added" begin
            hash1 = OutbreakDetectionCore.struct_definition_hash(TestStruct1)
            hash3 = OutbreakDetectionCore.struct_definition_hash(TestStruct3)

            @test hash1 != hash3
        end

        @testset "generates different hashes when field order changes" begin
            hash1 = OutbreakDetectionCore.struct_definition_hash(TestStruct1)
            hash4 = OutbreakDetectionCore.struct_definition_hash(TestStruct4)

            @test hash1 != hash4
        end

        @testset "generates different hashes when field types change" begin
            hash1 = OutbreakDetectionCore.struct_definition_hash(TestStruct1)
            hash5 = OutbreakDetectionCore.struct_definition_hash(TestStruct5)

            @test hash1 != hash5
        end

        @testset "works with OptimizationScenario" begin
            hash_val = OutbreakDetectionCore.struct_definition_hash(
                OptimizationScenario
            )

            @test hash_val isa UInt64
            @test hash_val != 0
        end

        @testset "works with OptimizationResult" begin
            hash_val = OutbreakDetectionCore.struct_definition_hash(OptimizationResult)

            @test hash_val isa UInt64
            @test hash_val != 0
        end

        @testset "OptimizationScenario and OptimizationResult have different hashes" begin
            scenario_hash = OutbreakDetectionCore.struct_definition_hash(
                OptimizationScenario
            )
            result_hash = OutbreakDetectionCore.struct_definition_hash(
                OptimizationResult
            )

            @test scenario_hash != result_hash
        end
    end

    @testset "get_optimization_struct_hashes" begin
        @testset "returns NamedTuple with correct keys" begin
            hashes = get_optimization_struct_hashes()

            @test hashes isa NamedTuple
            @test haskey(hashes, :scenario_hash)
            @test haskey(hashes, :result_hash)
            @test length(keys(hashes)) == 2
        end

        @testset "returns UInt64 values" begin
            hashes = get_optimization_struct_hashes()

            @test hashes.scenario_hash isa UInt64
            @test hashes.result_hash isa UInt64
        end

        @testset "returns non-zero hashes" begin
            hashes = get_optimization_struct_hashes()

            @test hashes.scenario_hash != 0
            @test hashes.result_hash != 0
        end

        @testset "scenario and result hashes are different" begin
            hashes = get_optimization_struct_hashes()

            @test hashes.scenario_hash != hashes.result_hash
        end

        @testset "generates consistent hashes across calls" begin
            hashes1 = get_optimization_struct_hashes()
            hashes2 = get_optimization_struct_hashes()

            @test hashes1.scenario_hash == hashes2.scenario_hash
            @test hashes1.result_hash == hashes2.result_hash
            @test hashes1 == hashes2
        end

        @testset "hashes match individual struct_definition_hash calls" begin
            hashes = get_optimization_struct_hashes()

            scenario_hash = OutbreakDetectionCore.struct_definition_hash(
                OptimizationScenario
            )
            result_hash = OutbreakDetectionCore.struct_definition_hash(
                OptimizationResult
            )

            @test hashes.scenario_hash == scenario_hash
            @test hashes.result_hash == result_hash
        end

        @testset "can be compared for equality" begin
            hashes1 = get_optimization_struct_hashes()
            hashes2 = get_optimization_struct_hashes()
            hashes3 = (
                scenario_hash = hashes1.scenario_hash + 1,
                result_hash = hashes1.result_hash,
            )

            @test hashes1 == hashes2
            @test hashes1 != hashes3
        end

        @testset "can be used in dictionary operations" begin
            hashes = get_optimization_struct_hashes()

            # Test that it can be used as a dictionary value
            dict = Dict("current" => hashes)
            @test dict["current"] == hashes

            # Test comparison with stored value
            stored_hashes = dict["current"]
            current_hashes = get_optimization_struct_hashes()
            @test stored_hashes == current_hashes
        end
    end
end
