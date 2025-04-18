#%%
using DrWatson
@quickactivate "OutbreakDetection"

using OutbreakDetectionUtils
using OutbreakDetection
using FLoops
using JLD2
using CSV
using DataFrames

include(srcdir("makie-plotting-setup.jl"))

include(srcdir("ensemble-parameters.jl"))
include(projectdir("manuscript", "scripts", "plotting-setup.jl"))

if false
    include("../src/ensemble-parameters.jl")
end

#%%
model_types_vec = [("seasonal-infectivity-import", "tau-leaping")]

#%%
N_vec = [500_000]
nsims_vec = [100]
init_states_prop_dict = [
    Dict(:s_prop => 0.05, :e_prop => 0.00, :i_prop => 0.00, :r_prop => 0.95)
]

ensemble_state_p_vec = create_combinations_vec(
    StateParameters, (N_vec, init_states_prop_dict)
)

#%%
tmin_vec = [0.0]
tstep_vec = [1.0]
tmax_vec = [365.0 * 100]

time_p_vec = vec(
    map(
        Iterators.product(tmin_vec, tstep_vec, tmax_vec)
    ) do (tmin, tstep, tmax)
        SimTimeParameters(; tmin = tmin, tmax = tmax, tstep = tstep)
    end,
)

#%%
beta_force_vec = [0.2]
annual_births_per_k_vec = [27]
seed = 1234

#%%
latent_per_days_vec = [LATENT_PER_DAYS]
dur_inf_days_vec = [DUR_INF_DAYS]
R_0_vec = [16.0]
sigma_vec = 1 ./ latent_per_days_vec
gamma_vec = 1 ./ dur_inf_days_vec
vaccination_coverage_vec = [0.8]

#%%
ensemble_spec_vec = create_ensemble_spec_combinations(
    beta_force_vec,
    [cos],
    sigma_vec,
    gamma_vec,
    annual_births_per_k_vec,
    R_0_vec,
    vaccination_coverage_vec,
    N_vec,
    init_states_prop_dict,
    model_types_vec,
    time_p_vec,
    nsims_vec,
)

#%%
outbreak_threshold_vec = [5]
min_outbreak_dur_vec = [30]
min_outbreak_size_vec = [500]

outbreak_spec_vec = create_combinations_vec(
    OutbreakSpecification,
    (outbreak_threshold_vec, min_outbreak_dur_vec, min_outbreak_size_vec),
)

#%%
# TODO: How to make the mean of the noise a proportion of the mean of incidence
# when incidence depends on the dynamics parameters?
# Could pass variables to ensemble function and calculate each simulations and
# scenario's noise mean, but that would break implementation using NoiseSpecification
# struct currently
poisson_noise_mean_scaling_vec = [
    1.0,
    collect(2.0:2.0:8.0)...,
]

poisson_noise_spec_vec = create_combinations_vec(
    PoissonNoiseSpecification,
    (["poisson"], poisson_noise_mean_scaling_vec),
)

dynamical_noise_R0 = [5.0]
dynamical_noise_latent_period = [7]
dynamical_noise_duration_infection = [14]
dynamical_noise_correlation = [
    "in-phase"
    # "out-of-phase",
    # "none",
]
dynamical_noise_mean_scaling_vec = [0.15]
dynamical_noise_vaccination_coverage_vec = [
    0.8538, 0.7383, 0.5088, 0.2789, 0.0492
]
dynamical_noise_spec_vec = create_combinations_vec(
    DynamicalNoiseSpecification,
    (
        ["dynamical"],
        dynamical_noise_R0,
        dynamical_noise_latent_period,
        dynamical_noise_duration_infection,
        dynamical_noise_correlation,
        dynamical_noise_mean_scaling_vec,
        dynamical_noise_vaccination_coverage_vec,
    ),
)

noise_spec_vec = vcat(
    poisson_noise_spec_vec,
    dynamical_noise_spec_vec,
)

#%%
# This will get updated in optimization so just use placeholder
alertthreshold_vec = 5.0
moveavglag_vec = [7]
perc_clinic_vec = [1.0]
perc_clinic_test_vec = collect(0.1:0.1:1.0)
alert_method_vec = ["movingavg"]

outbreak_detection_spec_vec = create_combinations_vec(
    OutbreakDetectionSpecification,
    (
        alertthreshold_vec,
        moveavglag_vec,
        perc_clinic_vec,
        perc_clinic_test_vec,
        alert_method_vec,
    ),
)

#%%
test_spec_vec = [
    # PROTOTYPE_RDT_TEST_SPECS...,
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    # CLINICAL_TEST_SPECS...,
    # IndividualTestSpecification(0.98, 0.98, 0),
    # IndividualTestSpecification(0.98, 0.98, 3),
    # IndividualTestSpecification(0.98, 0.98, 7),
    # IndividualTestSpecification(0.98, 0.98, 14),
    IndividualTestSpecification(1.0, 1.0, 0),
    # IndividualTestSpecification(1.0, 1.0, 3),
    # IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

#%%
accuracy_functions = [
    arithmetic_mean,
    calculate_f_beta_score,
]

#%%
optim_df = OutbreakDetectionUtils.run_scenario_optimizations(
    ensemble_spec_vec,
    outbreak_spec_vec,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec,
    MSO,
    accuracy_functions;
    executor = ThreadedEx(),
    force = true,
    save_df = true,
    return_df = true,
)

#%%
population_df = CSV.read(
    datadir("input-populations.csv"),
    DataFrame; delim = ',',
    header = true,
)

gha_2022_pop = only(
    population_df[population_df.ISO3_code .== "GHA", "2022"]
)
gha_2022_scale_population =
    gha_2022_pop / ensemble_state_specification.init_states.N

gha_2022_scale_population_per_annum = gha_2022_scale_population / 100

#%%
all_visit_clinic_optims_mean = subset(
    optim_df, :accuracy_function => f -> f .== arithmetic_mean
)
all_visit_clinic_optims_mean_thresholds_chars = reshape_optim_df_to_matrix(
    all_visit_clinic_optims_mean
);

#%%
imp_85_test_100pc_tested_threshold_chars = filter(
    c ->
        c.individual_test_specification ==
        IndividualTestSpecification(0.85, 0.85, 0) &&
            c.percent_clinic_tested == 1.0,
    all_visit_clinic_optims_mean_thresholds_chars[2, 4],
)[1];

imp_90_test_100pc_tested_threshold_chars = filter(
    c ->
        c.individual_test_specification ==
        IndividualTestSpecification(0.90, 0.90, 0) &&
            c.percent_clinic_tested == 1.0,
    all_visit_clinic_optims_mean_thresholds_chars[2, 4],
)[1];

#%%
mean(
    imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.proportion_timeseries_in_alert,
)
mean(
    imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.proportion_timeseries_in_alert,
)

#%%
mean(
    imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.noutbreaks
)
mean(
    imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.noutbreaks
)

#%%
mean(imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.nalerts)
mean(imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.nalerts)

#%%
mean(
    imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.n_false_alerts,
)
mean(
    imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.n_false_alerts,
)

#%%
mean(
    imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.n_missed_outbreaks,
)
mean(
    imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.n_missed_outbreaks,
)

#%%
mean(
    mapreduce(
        mean,
        vcat,
        imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.detectiondelays,
    ),
)
mean(
    mapreduce(
        mean,
        vcat,
        imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.detectiondelays,
    ),
)

#%%
mean(
    imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.perc_true_outbreaks_detected,
)
mean(
    imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.perc_true_outbreaks_detected,
)

#%%
mean(
    imp_85_test_100pc_tested_threshold_chars.outbreak_threshold_chars.perc_alerts_correct,
)
mean(
    imp_90_test_100pc_tested_threshold_chars.outbreak_threshold_chars.perc_alerts_correct,
)

#%%
imp_85_test_100pc_tested_threshold_chars.accuracy
imp_90_test_100pc_tested_threshold_chars.accuracy

#%%
for accuracy_function in unique(optim_df.accuracy_function)
    accuracy_function_str = String(nameof(accuracy_function))

    all_visit_clinic_optims = subset(
        optim_df,
        :outbreak_detection_spec =>
            x ->
                getproperty.(x, :percent_visit_clinic) .== 1.0,
        :accuracy_function => y -> y .== accuracy_function,
    )

    all_visit_clinic_optims_threshold_chars = reshape_optim_df_to_matrix(
        all_visit_clinic_optims
    )

    line_plot(
        all_visit_clinic_optims_threshold_chars;
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
        save_plot = true,
        clinical_hline = false,
        colors = OutbreakDetection.lineplot_colors,
        plotdirpath = DrWatson.plotsdir(
            "ensemble", "scenario-optimizations", "perc_visit_clinic_1.0",
            accuracy_function_str,
        ),
        plotname = "line_threshold_plot",
        plotformat = "png",
        size = (1300, 800),
    )

    line_plot(
        all_visit_clinic_optims_threshold_chars;
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
        save_plot = true,
        clinical_hline = false,
        colors = OutbreakDetection.lineplot_colors,
        plotdirpath = DrWatson.plotsdir(
            "ensemble", "scenario-optimizations", "perc_visit_clinic_1.0",
            accuracy_function_str,
        ),
        plotname = "line_accuracy_plot",
        plotformat = "png",
        size = (1300, 800),
    )

    line_plot(
        all_visit_clinic_optims_threshold_chars;
        outcome = :f1_score,
        ylabel = "F1 Score",
        alpha = alpha,
        facet_fontsize = 28,
        legendsize = legendsize,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        show_x_facet_label = true,
        show_y_facet_label = false,
        ylims = (0.5, 1.0),
        force = true,
        save_plot = true,
        clinical_hline = false,
        colors = OutbreakDetection.lineplot_colors,
        plotdirpath = DrWatson.plotsdir(
            "ensemble", "scenario-optimizations", "perc_visit_clinic_1.0",
            accuracy_function_str,
        ),
        plotname = "line_f1-score_plot",
        plotformat = "png",
        size = (1300, 800),
    )

    line_plot(
        all_visit_clinic_optims_threshold_chars;
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
        ylims = (-50, 100),
        force = true,
        save_plot = true,
        clinical_hline = false,
        colors = OutbreakDetection.lineplot_colors,
        plotdirpath = DrWatson.plotsdir(
            "ensemble", "scenario-optimizations", "perc_visit_clinic_1.0",
            accuracy_function_str,
        ),
        plotname = "line_delays_plot",
        plotformat = "png",
        size = (1300, 800),
    )

    line_plot(
        all_visit_clinic_optims_threshold_chars;
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
        ylims = (0, 0.5),
        force = true,
        save_plot = true,
        clinical_hline = false,
        colors = OutbreakDetection.lineplot_colors,
        plotdirpath = DrWatson.plotsdir(
            "ensemble", "scenario-optimizations", "perc_visit_clinic_1.0",
            accuracy_function_str,
        ),
        plotname = "line_prop_alert_plot",
        plotformat = "png",
        size = (1300, 800),
    )

    line_plot(
        all_visit_clinic_optims_threshold_chars;
        outcome = :unavoidable_cases,
        ylabel = "Unavoidable Cases",
        alpha = alpha,
        facet_fontsize = facet_fontsize,
        legendsize = legendsize,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        show_x_facet_label = true,
        show_y_facet_label = false,
        ylims = (0, 3.0e4),
        force = true,
        save_plot = true,
        clinical_hline = false,
        colors = OutbreakDetection.lineplot_colors,
        cases_scaling = gha_2022_scale_population_per_annum,
        plotdirpath = DrWatson.plotsdir(
            "ensemble", "scenario-optimizations", "perc_visit_clinic_1.0",
            accuracy_function_str,
        ),
        plotname = "line_unavoidable_plot",
        plotformat = "png",
        size = (1300, 800),
    )
end
