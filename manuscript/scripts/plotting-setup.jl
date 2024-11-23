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
    # fontsize = 16,
    linewidth = 6,
)

update_theme!(; size = (1300, 800))
CairoMakie.activate!(; type = "svg", px_per_unit = 1.0)

alpha = 0.4
nbanks = 1
facet_fontsize = 28
legendsize = 24
xlabelsize = 28
ylabelsize = 26
legend_rowsize = Makie.Relative(0.05)
xlabel_rowsize = Makie.Relative(0.01);
