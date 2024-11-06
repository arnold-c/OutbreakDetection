#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DataFrames
using CategoricalArrays
using Match: Match
using CSV: CSV
using Chain: Chain

using OutbreakDetectionUtils:
    create_optimal_thresholds_df, create_optimal_threshold_summary_df

import OutbreakDetectionUtils: create_optimal_threshold_summary_df

#%%
manuscriptdir(args...) = DrWatson.projectdir("manuscript", args...)
manuscript_files(args...) = manuscriptdir("manuscript_files", args...)
manuscript_plotdir(args...) = manuscript_files("plots", args...)
manuscript_tabledir(args...) = manuscript_files("tables", args...)
appendix_files(args...) = manuscriptdir("supplemental-appendix_files", args...)

#%%
include(manuscriptdir("optimal-thresholds_loading.jl"));

#%%
include(manuscriptdir("optimal-thresholds_tables.jl"));

#%%
include(manuscriptdir("optimal-thresholds_plots.jl"));

#%%
# include(manuscriptdir("optimal-thresholds_checks.jl"))
