module OutbreakDetection

using DrWatson
using GLMakie
using ColorSchemes
using UnPack
using DataFrames
using Chain
using NaNMath: NaNMath
using FLoops
using OutbreakDetectionUtils

include("plotting-helpers.jl")
export ACCURACY_COLOR, DAILY_SENSITIVITY_COLOR, DAILY_SPECIFICITY_COLOR,
    DAILY_PPV_COLOR, DAILY_NPV_COLOR, DETECTION_DELAY_COLOR,
    N_ALERTS_PER_OUTBREAK_COLOR, N_FALSE_ALERTS_COLOR, N_ALERTS_COLOR,
    N_OUTBREAKS_COLOR, N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR,
    PERC_OUTBREAKS_MISSED_COLOR, PERC_ALERTS_CORRECT_COLOR,
    PERC_ALERTS_FALSE_COLOR

include("line_plots.jl")
export line_plot,
    collect_OptimalThresholdCharacteristics

@static if false
    # Include manuscript scripts for LSP support
    # These scripts are not loaded at runtime but help the LSP
    # recognize plotting functions and provide autocomplete
    include("../scripts/manuscript/optimal-thresholds.jl")
    include("../scripts/manuscript/optimal-thresholds_loading.jl")
    include("../scripts/manuscript/optimal-thresholds_plots.jl")
    include("../scripts/manuscript/optimal-thresholds_checks.jl")
    include("../scripts/manuscript/optimal-thresholds_tables.jl")
    include("../scripts/manuscript/supplemental_plots.jl")
    include("../scripts/manuscript/supplemental_tables.jl")
    include("../scripts/manuscript/plotting-setup.jl")
    include("../scripts/manuscript/schematic-plot.jl")
end

end
