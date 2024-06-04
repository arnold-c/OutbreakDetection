seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
seir_state_labels = ["S", "E", "I", "R", "N"]

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

