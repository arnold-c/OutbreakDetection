using OutbreakDetectionCore
using Test
using Aqua
using JET

@testset "OutbreakDetectionCore.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OutbreakDetectionCore; ambiguities = false)
        @testset "Ambiguities" begin
            Aqua.test_ambiguities(OutbreakDetectionCore)
        end
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(OutbreakDetectionCore; target_defined_modules = true)
    end
    include("SEIR-model.jl")
    include("cleaning-functions.jl")
    include("diag-testing-functions.jl")
    include("detection-thresholds.jl")
    include("ensemble-functions.jl")
    include("noise-functions.jl")
    include("collect-thresholds-vec_functions.jl")
    include("optimal-threshold-functions.jl")
    include("optimization-wrapper.jl")
    include("outbreak-status.jl")
    include("outbreak-classification-outbreakspec.jl")
    include("alert-generation.jl")
    include("accuracy-calculation.jl")
    include("match-alert-outbreak-bounds.jl")
    include("detection-metric-functions.jl")
    include("missing-results.jl")
end
