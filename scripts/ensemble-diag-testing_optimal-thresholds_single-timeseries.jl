#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using CSV: CSV
using DataFrames
using DataFramesMeta
using Statistics

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
optimal_threshold_test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    CLINICAL_TEST_SPECS...,
    IndividualTestSpecification(1.0, 1.0, 0),
]

optimal_threshold_alertthreshold_vec = collect(1:1:30)

R_0_vec = collect(8.0:4.0:20.0)

ensemble_dynamics_spec_vec = create_combinations_vec(
    DynamicsParameters,
    (
        [ensemble_state_specification.init_states.N],
        [27],
        [0.2],
        [SIGMA],
        [GAMMA],
        R_0_vec,
        [0.8],
    ),
)

ensemble_spec_vec = create_combinations_vec(
    EnsembleSpecification,
    (
        [ensemble_model_type],
        [ensemble_state_specification],
        ensemble_dynamics_spec_vec,
        [ensemble_time_specification],
        [ensemble_nsims],
    ),
)

sim_number = 1

#%%
for (ensemble_noise_specification, ensemble_specification) in
    Iterators.product(ensemble_noise_specification_vec, ensemble_spec_vec)
    @info "Creating plots and tables for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification))"
    println("==============================================")

    optimal_threshold_comparison_params = (
        alertthreshold_vec = optimal_threshold_alertthreshold_vec,
        ensemble_specification = ensemble_specification,
        noise_specification = ensemble_noise_specification,
        outbreak_specification = ensemble_outbreak_specification,
        moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
        percent_visit_clinic = ensemble_percent_visit_clinic,
    )

    optimal_thresholds_vec = calculate_OptimalThresholdCharacteristics(
        ensemble_percent_clinic_tested_vec,
        optimal_threshold_test_spec_vec,
        optimal_threshold_comparison_params,
    )

    ensemble_single_scenario_inc_file = get_ensemble_file(
        ensemble_specification, ensemble_outbreak_specification
    )

    incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]

    noisearr, poisson_noise_prop = create_noise_arr(
        ensemble_noise_specification,
        incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )

    noise_specification_path = getdirpath(ensemble_noise_specification)
    noise_specification_filename = replace(
        noise_specification_path,
        "/" => "_",
    )

    baseplotdirpath = joinpath(
        plotsdir("ensemble/optimal-thresholds"),
        "R0_$(ensemble_specification.dynamics_parameters.R_0)",
        noise_specification_path,
        "single-scenario",
    )

    mkpath(baseplotdirpath)

    unique_percent_clinic_tested = unique(
        optimal_thresholds_vec.percent_clinic_tested
    )

    for percent_clinic_tested in unique_percent_clinic_tested
        optimal_thresholds_chars = filter(
            optimal_thresholds ->
                optimal_thresholds.percent_clinic_tested ==
                percent_clinic_tested ||
                    (
                        optimal_thresholds.percent_clinic_tested ==
                        1.0 &&
                        optimal_thresholds.individual_test_specification ==
                        CLINICAL_CASE_TEST_SPEC
                    ),
            optimal_thresholds_vec,
        )

        for optimal_thresholds_char in optimal_thresholds_chars
            ind_test_spec =
                optimal_thresholds_char.individual_test_specification

            detection_specification = OutbreakDetectionSpecification(
                optimal_thresholds_char.alert_threshold,
                ensemble_moving_avg_detection_lag,
                ensemble_percent_visit_clinic,
                percent_clinic_tested,
            )

            testarr = create_testing_arrs(
                incarr,
                noisearr,
                detection_specification,
                ind_test_spec
            )

            plot = incidence_testing_plot(
                incarr,
                noisearr,
                testarr,
                detection_specification,
                ensemble_time_specification;
                sim = sim_number,
                plotitle = "% Clinic Tested: $(percent_clinic_tested), Sens: $(ind_test_spec.sensitivity), Spec: $(ind_test_spec.specificity), Lag: $(ind_test_spec.test_result_lag)",
            )

            save(
                joinpath(
                    baseplotdirpath,
                    "single-scenario_clinic-tested-$(percent_clinic_tested)_sens-$(ind_test_spec.sensitivity)_spec-$(ind_test_spec.specificity)_lag-$(ind_test_spec.test_result_lag).png",
                ),
                plot;
                size = (2200, 1600),
            )
        end
    end

    @info "Timeseries $(sim_number) saved for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification))"
    println()
    println("==============================================")
    println()
end
