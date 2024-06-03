using OutbreakDetectionUtils
using Test
using Aqua
using JET

@testset "OutbreakDetectionUtils.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OutbreakDetectionUtils)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(OutbreakDetectionUtils; target_defined_modules = true)
    end
    # Write your tests here.
end
