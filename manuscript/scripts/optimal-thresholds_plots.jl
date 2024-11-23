#%%
accuracy_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :accuracy,
    ylabel = "Outbreak Detection\nAccuracy",
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
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

save(
    manuscript_plotdir("optimal-thresholds_accuracy-plot.svg"),
    accuracy_plot,
)

#%%
delays_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :detectiondelays,
    ylabel = "Detection Delays\n(Days)",
    alpha = alpha,
    hlines = (0.0),
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (-100, 100),
    force = true,
    plotdirpath = DrWatson.projectdir("manuscript"),
    plotname = "line_delays_plot",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

save(
    manuscript_plotdir("optimal-thresholds_delays-plot.svg"),
    delays_plot,
)

#%%
prop_alert_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :proportion_timeseries_in_alert,
    ylabel = "Proportion of Time Series\nIn Alert",
    alpha = alpha,
    hlines = (0.0),
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 0.35),
    force = true,
    plotdirpath = DrWatson.projectdir("manuscript"),
    plotname = "line_prop_alert_plot",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

save(
    manuscript_plotdir("optimal-thresholds_prop-alert-plot.svg"),
    prop_alert_plot,
)

#%%
unavoidable_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :unavoidable_cases,
    ylabel = "Unavoidable Cases",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 2.4e4),
    force = true,
    plotdirpath = DrWatson.projectdir("manuscript"),
    plotname = "line_unavoidable_plot",
    save_plot = true,
    clinical_hline = clinical_hline,
    cases_scaling = gha_2022_scale_population_per_annum,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

save(
    manuscript_plotdir("optimal-thresholds_unavoidable-plot.svg"),
    unavoidable_plot,
)
