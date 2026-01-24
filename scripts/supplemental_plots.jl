#%%
using DrWatson: DrWatson
using OutbreakDetectionCore: OutbreakDetectionCore
using OutbreakDetection
using Try: Try

include(
    DrWatson.scriptsdir("plotting-setup.jl")
)

# Make sure these values are present in the optimization script
alert_method = OutbreakDetectionCore.AlertMethod(OutbreakDetectionCore.MovingAverage())
accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.BalancedAccuracy())
threshold_bounds = (; lower = 0.0, upper = 20.0)
plotdirpath = DrWatson.plotsdir("appendix")

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
prop_outbreak_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :proportion_timeseries_in_outbreak,
    ylabel = "Proportion of Time\nSeries In Outbreak",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0.0, 0.3),
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_prop-outbreak-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
alert_duration_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :alert_durations,
    ylabel = "Alert Duration\n(Days)",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 130),
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_alert-duration-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
nalerts_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :n_alerts,
    ylabel = "Number of Alerts",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 450),
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_n-alerts-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
noutbreaks_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :n_outbreaks,
    ylabel = "Number of outbreaks",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 60),
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_n-outbreaks-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
prop_alerts_correct_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :proportion_alerts_correct,
    ylabel = "Proportion of\nAlerts Correct",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 1.0),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_prop-alerts-correct-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
prop_outbreaks_detected_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :proportion_outbreaks_detected,
    ylabel = "Proportion of\nOutbreaks Detected",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 1.0),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_prop-outbreaks-detected-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
# F1 metrics for supplement
accuracy_metric = OutbreakDetectionCore.AccuracyMetric(OutbreakDetectionCore.F1())

#%%
f1_accuracy_plot = line_plot(
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
    plotname = "f1_optimal-thresholds_accuracy-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
f1_delays_plot = line_plot(
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
    plotname = "f1_optimal-thresholds_delays-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
f1_prop_alert_plot = line_plot(
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
    plotname = "f1_optimal-thresholds_prop-alert-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
f1_prop_alerts_correct_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :proportion_alerts_correct,
    ylabel = "Proportion of\nAlerts Correct",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 1.0),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_prop-alerts-correct-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)

#%%
f1_prop_outbreaks_detected_plot = line_plot(
    optimized_threshold_results;
    alert_method = alert_method,
    accuracy_metric = accuracy_metric,
    threshold_bounds = threshold_bounds,
    outcome = :proportion_outbreaks_detected,
    ylabel = "Proportion of\nOutbreaks Detected",
    alpha = alpha,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    show_x_facet_label = show_x_facet_label,
    show_y_facet_label = show_y_facet_label,
    ylims = (0, 1.0),
    force = true,
    plotdirpath = plotdirpath,
    plotname = "optimal-thresholds_prop-outbreaks-detected-plot",
    plotformat = "svg",
    save_plot = true,
    nbanks = nbanks,
    legend_rowsize = legend_rowsize,
    xlabel_rowsize = xlabel_rowsize,
    size = (1300, 800),
)
