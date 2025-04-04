#%%
prop_outbreak_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :proportion_timeseries_in_outbreak,
    ylabel = "Proportion of Time\nSeries In Outbreak",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    force = true,
    save_plot = false,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

save(
    appendix_plotdir("optimal-thresholds_prop-outbreak-plot.svg"),
    prop_outbreak_plot,
)

#%%
alert_duration_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :alert_duration_vec,
    ylabel = "Alert Duration\n(Days)",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 170),
    force = true,
    save_plot = false,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

save(
    appendix_plotdir("optimal-thresholds_alert-duration-plot.svg"),
    alert_duration_plot,
)

#%%
nalerts_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :nalerts,
    ylabel = "Number of Alerts",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 250),
    force = true,
    save_plot = false,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

save(
    appendix_plotdir("optimal-thresholds_n-alerts-plot.svg"),
    nalerts_plot,
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
    plotdirpath = appendix_plotdir(),
    plotname = "optimal-thresholds_unavoidable-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    cases_scaling = gha_2022_scale_population_per_annum,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

#%%
f1_threshold_plot = line_plot(
    f1_optimal_threshold_characteristics;
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
    plotdirpath = appendix_plotdir(),
    plotname = "f1_optimal-thresholds_alert-threshold-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
f1_accuracy_plot = line_plot(
    f1_optimal_threshold_characteristics;
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
    plotdirpath = appendix_plotdir(),
    plotname = "f1_optimal-thresholds_accuracy-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
f1_delays_plot = line_plot(
    f1_optimal_threshold_characteristics;
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
    plotdirpath = appendix_plotdir(),
    plotname = "f1_optimal-thresholds_delays-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

#%%
f1_prop_alert_plot = line_plot(
    f1_optimal_threshold_characteristics;
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
    plotdirpath = appendix_plotdir(),
    plotname = "f1_optimal-thresholds_prop-alert-plot",
    plotformat = "svg",
    save_plot = true,
    clinical_hline = clinical_hline,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)
