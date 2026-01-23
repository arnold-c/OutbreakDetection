module OutbreakDetectionCore

# Package imports

using AutoHashEquals: AutoHashEquals
using BangBang: BangBang
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
using LightSumTypes: LightSumTypes, @sumtype
using LinearAlgebra: LinearAlgebra
using MultistartOptimization: MultistartOptimization
using NaNMath: NaNMath
using NLopt: NLopt
using OhMyThreads: OhMyThreads
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

# Generic Utilities
include("./utilities/DrWatson-helpers.jl")

# Types
include("./types/time-parameters.jl")
include("./types/state-parameters.jl")
include("./types/dynamics-parameters.jl")
include("./types/simulation-results.jl")
include("./types/noise-specifications.jl")
include("./types/ensemble-specifications.jl")
include("./types/outbreak-specifications.jl")
include("./types/alert-methods.jl")
include("./types/outbreak-detection-specifications.jl")
include("./types/thresholds.jl")
include("./types/test-specifications.jl")
include("./types/accuracy-metrics.jl")
include("./types/optimization-methods.jl")
include("./types/optimization-parameters.jl")
include("./types/optimization-scenario.jl")
include("./types/optimization-results.jl")
include("./types/optimization-tracker.jl")
include("./types/test-positive-containers.jl")

# Simulation (needed before constants)
include("./simulation/transmission-functions.jl")
include("./simulation/simulate-ensemble-seir-results.jl")

# Constants
include("./constants/test-constants.jl")

# Simulation (continued)
include("./simulation/endemic-equilibrium.jl")
include("./simulation/seir-model.jl")
include("./simulation/validate-ensemble-outbreaks.jl")

# Noise
include("./noise/noise-dynamics-parameters.jl")
include("./noise/noise-generation.jl")
include("./noise/noise-mean-incidence.jl")
include("./noise/noise-parameters-optimization.jl")
include("./noise/noise-recreation.jl")

# Utilities (data processing)
include("./utilities/calculate-moving-average.jl")
include("./utilities/create-ensemble-specification.jl")
include("./utilities/group-structvectors.jl")
include("./utilities/calculate-mean-incidence.jl")
include("./utilities/vaccination-distribution-sample.jl")
include("./utilities/test-descriptions.jl")
# TODO: update to work with current structvector implementation
include("./utilities/results-reshaping.jl")

# Diagnostic Testing
include("./diagnostic-testing/calculate-num-tested.jl")
include("./diagnostic-testing/calculate-num-positive.jl")
include("./diagnostic-testing/create-test-positive-vectors.jl")
include("./diagnostic-testing/create-test-positive-container.jl")

# Detection
include("./detection/accuracy-ensemble-calculation.jl")
include("./detection/accuracy-simulation-calculation.jl")
include("./detection/alert-generation.jl")
include("./detection/classify-outbreaks.jl")
include("./detection/detection-metric-functions.jl")
include("./detection/match-alert-outbreak-thresholds.jl")
include("./detection/outbreak-thresholds-calculation.jl")
include("./detection/threshold-bounds-calculation.jl")

# Threshold Optimization
include("./threshold-optimization/checkpoint-loading.jl")
include("./threshold-optimization/checkpoint-save.jl")
include("./threshold-optimization/confirm-proceed-with-optimization.jl")
include("./threshold-optimization/evaluate-missing-optimizations.jl")
include("./threshold-optimization/missing-results.jl")
include("./threshold-optimization/multistart-objective-function.jl")
include("./threshold-optimization/optimization-wrapper.jl")
include("./threshold-optimization/result-loading.jl")
include("./threshold-optimization/results-retrieval.jl")
include("./threshold-optimization/scenario-creation.jl")
include("./threshold-optimization/threshold-optimization.jl")

# Optimal Thresholds
include("./optimal-thresholds/optimal-metrics-calculation.jl")
# TODO: Update summaries below to be able to remove these files/adapt the underlying functions
include("./optimal-thresholds/results-dataframes.jl")
include("./optimal-thresholds/threshold-calculation.jl")
include("./optimal-thresholds/threshold-summaries.jl")

end
