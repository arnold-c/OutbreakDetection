using StructArrays

function calculate_OptimalThresholdCharacteristics(
    percent_clinic_tested_vec,
    ind_test_spec_vec,
    base_parameters
)
    optimal_thresholds_vec = Vector{OptimalThresholdCharacteristics}(
        undef, length(percent_clinic_tested_vec) * length(ind_test_spec_vec)
    )

    @showprogress for (i, (percent_clinic_tested, ind_test_spec)) in enumerate(
        Iterators.product(percent_clinic_tested_vec, ind_test_spec_vec)
    )
        optimal_thresholds_vec[i] = calculate_optimal_threshold(
            percent_clinic_tested,
            ind_test_spec,
            base_parameters
        )
    end

    return StructArray(optimal_thresholds_vec)
end

function calculate_optimal_threshold(
    percent_clinic_tested,
    individual_test_specification,
    base_parameters
)
    @unpack detectthreshold_vec,
    ensemble_specification,
    noise_specification,
    outbreak_specification,
    moving_avg_detection_lag,
    test_result_lag,
    percent_visit_clinic = base_parameters

    accuracy_array = zeros(Float64, length(detectthreshold_vec))

    ensemble_scenario_spec_vec = map(
        threshold -> ScenarioSpecification(
            ensemble_specification,
            outbreak_specification,
            noise_specification,
            OutbreakDetectionSpecification(
                threshold,
                moving_avg_detection_lag,
                percent_visit_clinic,
                percent_clinic_tested,
                test_result_lag,
            ),
            individual_test_specification,
        ),
        detectthreshold_vec,
    )

    for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        ensemble_chars_file = get_ensemble_file(ensemble_scenario_spec)

        accuracy_array[i] = median(ensemble_chars_file["OT_chars"].accuracy)
    end

    @views optimal_accuracy = NaNMath.maximum(accuracy_array)
    optimal_threshold_index = findfirst(==(optimal_accuracy), accuracy_array)
    optimal_threshold = detectthreshold_vec[optimal_threshold_index]

    return OptimalThresholdCharacteristics(
        ensemble_scenario_spec_vec[optimal_threshold_index],
        optimal_threshold,
        optimal_accuracy,
    )
end
