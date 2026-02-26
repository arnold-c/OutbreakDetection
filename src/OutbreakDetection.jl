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
using StaticArrays: StaticArrays
using Dates: Dates
using LightSumTypes: LightSumTypes
using StyledStrings
using Printf

include("./utilities.jl")
include("./plotting-helpers.jl")
include("./sort_test_specifications.jl")
include("./schematic-plot/simulation-setup.jl")
include("./schematic-plot/plot.jl")

include("./test-descriptions.jl")
include("./compare-optimal-solution-results.jl")
include("./optimal-thresholds_df-utilities.jl")
include("./threshold-summaries.jl")
include("./optimal-thresholds_wide-df.jl")

include("./line_plots.jl")

@static if false
    # Include manuscript scripts for LSP support
    # These scripts are not loaded at runtime but help the LSP
    # recognize plotting functions and provide autocomplete
    include("../scripts/plotting-setup.jl")
    include("../scripts/optimal-thresholds_optims.jl")
    include("../scripts/optimal-thresholds_plots.jl")
    include("../scripts/optimal-thresholds_comparisons.jl")
    include("../scripts/optimal-thresholds_supplement-plots.jl")
    include("../scripts/optimal-thresholds_supplement-tables.jl")
    include("../scripts/schematic-plot.jl")
end

end
