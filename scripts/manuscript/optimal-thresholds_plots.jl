#%%
threshold_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :alert_threshold,
    ylabel = "Alert Threshold",
    alpha = alpha,
    facet_fontsize = 28,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0.0, 18.0),
    force = true,
    plotdirpath = manuscript_plotdir(),
    plotname = "optimal-thresholds_alert-threshold-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
accuracy_plot = line_plot(
    optimal_threshold_characteristics;
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
    plotdirpath = manuscript_plotdir(),
    plotname = "optimal-thresholds_accuracy-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
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
    plotdirpath = manuscript_plotdir(),
    plotname = "optimal-thresholds_delays-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

#%%
prop_alert_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :proportion_timeseries_in_alert,
    ylabel = "Proportion of Time\nIn Alert",
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
    plotdirpath = manuscript_plotdir(),
    plotname = "optimal-thresholds_prop-alert-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)
