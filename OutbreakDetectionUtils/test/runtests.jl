using OutbreakDetectionUtils
using Test
using Aqua
using JET

@testset "OutbreakDetectionUtils.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OutbreakDetectionUtils; ambiguities = false)
        @testset "Ambiguities" begin
            Aqua.test_ambiguities(OutbreakDetectionUtils)
        end
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(OutbreakDetectionUtils; target_defined_modules = true)
    end
    include("diag-testing-functions.jl")
end
