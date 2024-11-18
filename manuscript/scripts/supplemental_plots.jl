#%%
prop_outbreak_plot = line_plot(
    optimal_threshold_characteristics;
    outcome = :proportion_timeseries_in_outbreak,
    ylabel = "Proportion of Time\nSeries In Outbreak",
    alpha = alpha,
    facet_fontsize = 20,
    legendsize = 20,
    xlabelsize = 22,
    ylabelsize = 22,
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
    facet_fontsize = 20,
    legendsize = 20,
    xlabelsize = 22,
    ylabelsize = 22,
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
    facet_fontsize = 20,
    legendsize = 20,
    xlabelsize = 22,
    ylabelsize = 22,
    show_x_facet_label = true,
    show_y_facet_label = false,
    ylims = (0, 350),
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
