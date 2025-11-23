module OutbreakDetectionCore

# Package imports
using Bootstrap: Bootstrap
using Bumper: Bumper, @no_escape, @alloc
using Chain: Chain
using DataFrames: DataFrames
using DataFramesMeta: DataFramesMeta
using Dates: Dates
using Distributions: Distributions
using DrWatson: DrWatson, @dict, @tagsave
using FLoops: FLoops
using FreqTables: FreqTables
using JLD2: JLD2
using LabelledArrays: SLVector, SLArray
using LightSumTypes: LightSumTypes, @sum_type
using LinearAlgebra: LinearAlgebra
using MultistartOptimization: MultistartOptimization
using NaNMath: NaNMath
using NLopt: NLopt
using ProgressMeter: ProgressMeter, Progress, next!
using Random: Random
using REPL.TerminalMenus: RadioMenu, request
using StaticArrays: StaticArrays
using Statistics: Statistics
using StatsBase: StatsBase
using StructArrays: StructArrays, StructVector
using StyledStrings
using Tables: Tables
using Try: Try
using UnPack: UnPack
using XLSX: XLSX

# Utilities
include("./utilities/DrWatson-helpers.jl")

# Types
include("./types/time-parameters.jl")
include("./types/state-parameters.jl")
include("./types/dynamics-parameters.jl")
include("./types/simulation-results.jl")
include("./types/ensemble-specifications.jl")
include("./types/outbreak-specifications.jl")
include("./types/detection-specifications.jl")
include("./types/test-specifications.jl")
include("./types/noise-specifications.jl")
include("./types/scenario-specifications.jl")
include("./types/scenario-parameters.jl")
include("./types/optimization-types.jl")
include("./types/optimization-results.jl")

# Simulation (needed before constants)
include("./simulation/transmission-functions.jl")

# Constants
include("./constants/test-constants.jl")

# Simulation (continued)
include("./simulation/endemic-equilibrium.jl")
include("./simulation/seir-model.jl")
include("./simulation/ensemble-simulation.jl")
include("./simulation/ensemble-outbreak-detection.jl")

# Noise
include("./noise/noise-mean-incidence.jl")
include("./noise/noise-dynamics-parameters.jl")
include("./noise/noise-recreation.jl")
include("./noise/noise-parameters-optimization.jl")
include("./noise/noise-generation.jl")

# Utilities (data processing)
include("./utilities/cleaning-functions.jl")
include("./utilities/collect-thresholds-vec_functions.jl")
include("./utilities/create-combinations.jl")

# Detection
include("./detection/outbreak-thresholds.jl")
include("./detection/outbreak-classification.jl")
include("./detection/moving-average.jl")
include("./detection/alert-detection.jl")
include("./detection/detection-characteristics.jl")

# Diagnostic Testing
include("./diagnostic-testing/calculate-num-tested.jl")
include("./diagnostic-testing/calculate-num-positive.jl")
include("./diagnostic-testing/calculate-test-positivity.jl")
include("./diagnostic-testing/create-test-arrays.jl")

# Optimal Thresholds
include("./optimal-thresholds/threshold-calculation.jl")
include("./optimal-thresholds/threshold-summaries.jl")
include("./optimal-thresholds/results-dataframes.jl")
include("./optimal-thresholds/results-export.jl")

# Threshold Optimization
include("./optimization-functions/common/objective-function.jl")
include("./optimization-functions/threshold-optimization/optimization-setup.jl")
include("./optimization-functions/threshold-optimization/optimization-wrapper.jl")

# Scenario Optimization
include("./optimization-functions/scenario-optimization/scenario-creation.jl")
include("./optimization-functions/scenario-optimization/scenario-confirmation.jl")
include("./optimization-functions/scenario-optimization/results-loading.jl")
include("./optimization-functions/scenario-optimization/results-reshaping.jl")
include("./optimization-functions/scenario-optimization/optimization-wrapper.jl")

end
