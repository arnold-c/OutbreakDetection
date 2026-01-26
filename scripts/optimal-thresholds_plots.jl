#%%
using DrWatson: DrWatson
using OutbreakDetectionCore: OutbreakDetectionCore
using OutbreakDetection
using Try: Try

include(
    DrWatson.scriptsdir("plotting-setup.jl")
)

# Make sure these values are present in the optimization script
alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage(7))
accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.BalancedAccuracy())
threshold_bounds = (; lower = 0.0, upper = 20.0)
plotdirpath = DrWatson.plotsdir()

#%%
optimized_filedir = OutbreakDetectionCore.outdir("ensemble", "threshold-optimization")
optimized_threshold_results = Try.unwrap(
    OutbreakDetectionCore.load_previous_optimization_results_structvector(
        optimized_filedir,
        "threshold-optimization.jld2",
        joinpath(optimized_filedir, "checkpoints"),
    )
)

#%%
threshold_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :alert_threshold,
    ylabel = "Alert Threshold",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0.0, 20.0),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_alert-threshold-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
accuracy_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :accuracies,
    ylabel = "Detection Accuracy",
    alpha = alpha,
    facet_fontsize = 28,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0.5, 1.0),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_accuracy-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
delays_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :detection_delays,
    ylabel = "Detection Delays\n(Days)",
    alpha = alpha,
    hlines = (0.0),
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (-20, 100),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_delays-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

#%%
prop_alert_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :proportion_timeseries_in_alert,
    ylabel = "Proportion of Time\nIn Alert",
    alpha = alpha,
    hlines = (0.0),
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 0.15),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_prop-alert-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
)

#%%
unavoidable_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :unavoidable_cases,
    ylabel = "Unavoidable Cases",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 1.0e3),
    force = true,
    # TODO: previously had a scaling factor for pop size
    # Look into and confirm not needed
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_unavoidable-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)
