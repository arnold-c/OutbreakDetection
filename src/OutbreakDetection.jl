module OutbreakDetection

using DrWatson: DrWatson
using StatsBase: StatsBase
using UnPack: @unpack
using DataFrames: DataFrames
using Chain: @chain
using NaNMath: NaNMath
using FLoops: FLoops
using OutbreakDetectionCore: OutbreakDetectionCore
using Match: Match
using StructArrays: StructVector, StructArray
using LaTeXStrings
using GLMakie
using ColorSchemes


include("plotting-helpers.jl")
include("sort_test_specifications.jl")
include("schematic-plot/simulation-setup.jl")
include("schematic-plot/plot.jl")

include("line_plots.jl")

@static if false
    # Include manuscript scripts for LSP support
    # These scripts are not loaded at runtime but help the LSP
    # recognize plotting functions and provide autocomplete
    include("../scripts/optimal-thresholds_optims.jl")
    include("../scripts/manuscript/optimal-thresholds_plots.jl")
    include("../scripts/manuscript/optimal-thresholds_checks.jl")
    include("../scripts/manuscript/optimal-thresholds_tables.jl")
    include("../scripts/manuscript/supplemental_plots.jl")
    include("../scripts/manuscript/supplemental_tables.jl")
    include("../scripts/manuscript/plotting-setup.jl")
    include("../scripts/manuscript/schematic-plot.jl")
end

end
