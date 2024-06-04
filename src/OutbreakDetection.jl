module OutbreakDetection

using GLMakie
using ColorSchemes
using UnPack
using DataFrames
using Chain
using NaNMath: NaNMath
using FLoops

include("plotting-helpers.jl")

include("single-sim_plots.jl")
export single_seir_plot, single_seir_statespace_plot, single_seir_beta_plot

include("quantile-plots.jl")
export create_seir_quantiles_plot

include("ensemble-sim_single-scenario_plots.jl")

include("threshold_comparison_plots.jl")
# export single_seir_plot

@static if false
    include("../scripts/single-sim_plots.jl")
    include("../scripts/ensemble-diag-testing_scenarios_plots.jl")
end

end
