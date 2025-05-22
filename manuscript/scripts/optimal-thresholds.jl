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

using OutbreakDetection: line_plot

#%%
manuscriptdir(args...) = DrWatson.projectdir("manuscript", args...)
manuscript_scripts(args...) = manuscriptdir("scripts", args...)
manuscript_files(args...) = manuscriptdir("manuscript_files", args...)
manuscript_plotdir(args...) = manuscript_files("plots", args...)
manuscript_tabledir(args...) = manuscript_files("tables", args...)
appendix_files(args...) = manuscriptdir("supplemental_files", args...)
appendix_plotdir(args...) = appendix_files("plots", args...)
appendix_tabledir(args...) = appendix_files("tables", args...)
include(manuscript_scripts("plotting-setup.jl"))

#%%
include(manuscript_scripts("schematic-plot.jl"))

#%%
include(manuscript_scripts("optimal-thresholds_loading.jl"));

#%%
include(manuscript_scripts("optimal-thresholds_plots.jl"));

#%%
include(manuscript_scripts("supplemental_tables.jl"))

#%%
include(manuscript_scripts("supplemental_plots.jl"))

#%%
include(manuscript_scripts("optimal-thresholds_checks.jl"))
