#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DataFrames
using CategoricalArrays
using Match: Match
using CSV: CSV
using Chain: Chain

using OutbreakDetectionCore:
    create_optimal_thresholds_df, create_optimal_threshold_summary_df

import OutbreakDetectionCore: create_optimal_threshold_summary_df

using OutbreakDetection: line_plot

#%%
tablesdir(args...) = projectdir("tables", args...)
appendix_plotdir(args...) = plotsdir("supplemental", args...)
appendix_tabledir(args...) = tablesdir("supplemental", args...)
include(scriptsdir("plotting-setup.jl"))

#%%
include(scriptsdir("dynamics-constants.jl"))

#%%
include(scriptsdir("schematic-plot.jl"))

#%%
include(scriptsdir("optimal-thresholds_loading.jl"));

#%%
include(scriptsdir("optimal-thresholds_plots.jl"));

#%%
include(scriptsdir("supplemental_tables.jl"))

#%%
include(scriptsdir("supplemental_plots.jl"))

#%%
include(scriptsdir("optimal-thresholds_checks.jl"))
