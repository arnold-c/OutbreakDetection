seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
seir_state_labels = ["S", "E", "I", "R", "N"]

ACCURACY_COLOR = "#004643"
DAILY_SENSITIVITY_COLOR = "#2A3965"
DAILY_SPECIFICITY_COLOR = "#C31D60"
DAILY_PPV_COLOR = "#22D37D"
DAILY_NPV_COLOR = "#6f366bff"
DETECTION_DELAY_COLOR = "#AE560A"
N_ALERTS_PER_OUTBREAK_COLOR = "#86B1A3"
N_FALSE_ALERTS_COLOR = "#D06778"
N_ALERTS_COLOR = "#00857E"
N_OUTBREAKS_COLOR = "#F4A157"
N_MISSED_OUTBREAKS_COLOR = "#5E5C6C"
PERC_OUTBREAKS_DETECTED_COLOR = "#F0780F"
PERC_OUTBREAKS_MISSED_COLOR = "#3A3842"
PERC_ALERTS_CORRECT_COLOR = "#004643"
PERC_ALERTS_FALSE_COLOR = "#852938"

function time_function(t; annual=true)
    return annual == true ? t / 365.0 : t
end

function annual_label(; annual=true)
    return annual == true ? "Time (years)" : "Time (days)"
end

function lims_check(ax, kwargs)
    if haskey(kwargs, :xlims)
        xlims!(ax, kwargs[:xlims])
    end

    if haskey(kwargs, :ylims)
        ylims!(ax, kwargs[:ylims])
    end
end

