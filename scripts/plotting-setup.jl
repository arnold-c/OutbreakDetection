using CairoMakie

function theme_adjustments()
    return Theme(;
        fontsize = 24,
        Axis = (;
            xlabelsize = 28,
            ylabelsize = 28,
            xlabelfont = :bold,
            ylabelfont = :bold,
        ),
        Colorbar = (;
            labelsize = 24,
            labelfont = :bold,
        ),
    )
end
custom_theme = merge(theme_adjustments(), theme_minimal())

set_theme!(
    custom_theme;
    fontsize = 34,
    linewidth = 12,
)

update_theme!(; size = (1300, 800))
CairoMakie.activate!()

alpha = 0.4
nbanks = 1
facet_fontsize = 50
legendsize = 38
xlabelsize = 38
ylabelsize = 38
legend_rowsize = Makie.Relative(0.05)
xlabel_rowsize = Makie.Relative(0.01)
show_x_facet_label = true
show_y_facet_label = false
markersize = 35
