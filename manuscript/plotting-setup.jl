using CairoMakie

function theme_adjustments()
    return Theme(;
        fontsize = 16,
        Axis = (;
            xlabelsize = 20,
            ylabelsize = 20,
            xlabelfont = :bold,
            ylabelfont = :bold,
        ),
        Colorbar = (;
            labelsize = 20,
            labelfont = :bold,
        ),
    )
end
custom_theme = merge(theme_adjustments(), theme_minimal())

set_theme!(
    custom_theme;
    fontsize = 16,
    linewidth = 4,
)

update_theme!(; size = (1300, 800))
CairoMakie.activate!(; type = "svg", px_per_unit = 1.0)
