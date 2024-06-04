using ODPlots
using Test
using Aqua
using JET

@testset "ODPlots.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ODPlots)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(ODPlots; target_defined_modules = true)
    end
    # Write your tests here.
end
