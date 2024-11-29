#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetection: line_plot

include(projectdir("manuscript", "scripts", "plotting-setup.jl"))
include(projectdir("manuscript", "scripts", "optimal-thresholds_loading.jl"))

#%%
poisson_optimal_threshold_characteristics = filter(
    optimal_thresholds ->
        get_noise_description(optimal_thresholds.noise_specification[1]) ==
        "poisson",
    optimal_threshold_characteristics,
);

dynamical_optimal_threshold_characteristics = filter(
    optimal_thresholds ->
        get_noise_description(optimal_thresholds.noise_specification[1]) ==
        "dynamical, in-phase",
    optimal_threshold_characteristics,
);

#%%
unique(
    map(
        char -> get_noise_description(char[1].noise_specification),
        poisson_optimal_threshold_characteristics,
    ),
)

#%%
poisson_accuracy_plot = line_plot(
    permutedims(poisson_optimal_threshold_characteristics);
    outcome = :accuracy,
    ylabel = "Detection Accuracy",
    alpha = alpha,
    facet_fontsize = 28,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.5, 1.0),
    force = true,
    plotdirpath = DrWatson.projectdir("manuscript"),
    plotname = "line_accuracy_plot",
    save_plot = false,
    clinical_hline = clinical_hline,
    nbanks = 1,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
)

save(
    plotsdir("optimal-thresholds_accuracy-plot_poisson.svg"),
    poisson_accuracy_plot;
    size = (1300, 400),
)

#%%
dynamical_accuracy_plot = line_plot(
    permutedims(dynamical_optimal_threshold_characteristics);
    outcome = :accuracy,
    ylabel = "Detection Accuracy",
    alpha = alpha,
    facet_fontsize = 28,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.5, 1.0),
    force = true,
    plotdirpath = DrWatson.projectdir("manuscript"),
    plotname = "line_accuracy_plot",
    save_plot = false,
    clinical_hline = clinical_hline,
    nbanks = 1,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
)

save(
    plotsdir("optimal-thresholds_accuracy-plot_dynamical.svg"),
    dynamical_accuracy_plot;
    size = (1300, 400),
)
