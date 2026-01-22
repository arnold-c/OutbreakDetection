using OutbreakDetectionCore
using Test
using Aqua
using JET

@testset "OutbreakDetectionCore.jl" begin
    # @testset "Code quality (Aqua.jl)" begin
    #     Aqua.test_all(OutbreakDetectionCore; ambiguities = false)
    #     @testset "Ambiguities" begin
    #         Aqua.test_ambiguities(OutbreakDetectionCore)
    #     end
    # end
    # @testset "Code linting (JET.jl)" begin
    #     JET.test_package(OutbreakDetectionCore; target_defined_modules = true)
    # end

    # Simulation tests
    include("simulation/seir-model.jl")
    # include("simulation/simulate-ensemble-seir-results.jl")  # Empty test file

    # Detection tests
    # include("detection/outbreak-thresholds-calculation.jl")  # API changed significantly
    # include("detection/alert-generation.jl")  # Not tested yet
    # include("detection/accuracy-simulation-calculation.jl")  # Not tested yet
    # include("detection/match-alert-outbreak-thresholds.jl")  # Not tested yet
    # include("detection/detection-metric-functions.jl")  # Not tested yet

    # Diagnostic testing tests
    # include("diagnostic-testing/calculate-num-positive.jl")  # Many functions missing

    # Noise tests
    # include("noise/noise-generation.jl")  # API changed significantly

    # Optimal thresholds tests
    # include("optimal-thresholds/optimal-metrics-calculation.jl")  # Not tested yet

    # Threshold optimization tests
    # include("threshold-optimization/missing-results.jl")  # Not tested yet
    # include("threshold-optimization/multistart-objective-function.jl")  # Not tested yet
    # include("threshold-optimization/optimization-wrapper.jl")  # Not tested yet
    # include("threshold-optimization/threshold-optimization.jl")  # Not tested yet

    # Other tests that need to be categorized
    # include("collect-thresholds-vec_functions.jl")
    # include("outbreak-status.jl")
    # include("outbreak-classification-outbreakspec.jl")
    # include("new-type-system.jl")
end
