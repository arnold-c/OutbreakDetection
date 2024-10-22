#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils:
    IndividualTestSpecification, DynamicalNoiseSpecification
using StatsBase: mean

include(srcdir("makie-plotting-setup.jl"))

include(projectdir("manuscript", "optimal-thresholds-loading.jl"));

if false
    include("optimal-thresholds-loading.jl")
end

#%%
mean_elisa_0d_delays = map(
    test_chars -> mean(vcat(getproperty(test_chars, :detectiondelays)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(1.0, 1.0, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_90_8x_dynamical_delays = map(
    test_chars -> mean(vcat(getproperty(test_chars, :detectiondelays)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.9, 0.9, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_85_8x_dynamical_delays = map(
    test_chars -> mean(vcat(getproperty(test_chars, :detectiondelays)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.85, 0.85, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mapreduce(
    vcat,
    (
        ("ELISA", mean_elisa_0d_delays),
        ("90%", mean_rdt_90_8x_dynamical_delays),
        ("85%", mean_rdt_85_8x_dynamical_delays),
    ),
) do (label, mean_delays_vec)
    Dict(label => round.(extrema(mean_delays_vec); digits = 1))
end

#%%
mean_elisa_0d_alert_duration = map(
    test_chars -> mean(vcat(getproperty(test_chars, :alert_duration_vec)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(1.0, 1.0, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_90_8x_dynamical_alert_duration = map(
    test_chars -> mean(vcat(getproperty(test_chars, :alert_duration_vec)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.9, 0.9, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mean_rdt_85_8x_dynamical_alert_duration = map(
    test_chars -> mean(vcat(getproperty(test_chars, :alert_duration_vec)...)),
    filter(
        chars ->
            chars.individual_test_specification ==
            IndividualTestSpecification(0.85, 0.85, 0),
        dynamical_noise_optimal_solutions,
    ).outbreak_threshold_chars,
)

mapreduce(
    vcat,
    (
        ("ELISA", mean_elisa_0d_alert_duration),
        ("90%", mean_rdt_90_8x_dynamical_alert_duration),
        ("85%", mean_rdt_85_8x_dynamical_alert_duration),
    ),
) do (label, mean_alert_duration_vec)
    Dict(label => round.(extrema(mean_alert_duration_vec); digits = 1))
end

#%%
# outbreak_proportion_line_plot = line_plot(
#     optimal_threshold_characteristics;
#     outcome = :proportion_timeseries_in_outbreak,
#     ylabel = "Proportion of Time Series\nIn Outbreak",
#     plotdirpath = baseplotdirpath,
#     facet_fontsize = 18,
#     labelsize = 20,
#     show_x_facet_label = true,
#     show_y_facet_label = false,
#     ylims = (0.0, 0.25),
#     force = true,
#     save_plot = false,
# clinical_hline = clinical_hline,
# )
#
# #%%

# #%%
# nalerts_line_plot = line_plot(
#     optimal_threshold_characteristics;
#     outcome = :nalerts,
#     ylabel = "Number of Alerts",
#     plotdirpath = baseplotdirpath,
#     facet_fontsize = 18,
#     labelsize = 20,
#     show_x_facet_label = true,
#     show_y_facet_label = false,
#     ylims = (0, 350),
#     force = true,
#     save_plot = false,
# clinical_hline = clinical_hline,
# )
#
# #%%
# nalerts_per_outbreak_line_plot = line_plot(
#     optimal_threshold_characteristics;
#     outcome = :n_alerts_per_outbreak,
#     ylabel = "Number of Alerts per Outbreak",
#     plotdirpath = baseplotdirpath,
#     facet_fontsize = 18,
#     labelsize = 20,
#     show_x_facet_label = true,
#     show_y_facet_label = false,
#     ylims = (0, 9),
#     force = true,
#     save_plot = false,
# clinical_hline = clinical_hline,
# )
#
# #%%
# alert_outbreak_proportion_line_plot = line_plot(
#     optimal_threshold_characteristics;
#     outcome = :alert_outbreak_timeseries_prop_diff,
#     ylabel = "Proportion of Time Series\nIn Alert - Outbreak",
#     plotdirpath = baseplotdirpath,
#     hlines = (0.0),
#     facet_fontsize = 18,
#     labelsize = 20,
#     show_x_facet_label = true,
#     show_y_facet_label = false,
#     ylims = (-0.15, 0.20),
#     force = true,
#     save_plot = false,
# clinical_hline = clinical_hline,
# )
