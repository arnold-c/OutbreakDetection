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

    ensemble_scenario_spec_vec = Vector{ScenarioSpecification}(
        undef,
        length(detectthreshold_vec)
    )
    ensemble_chars_vec = Vector(
        undef, length(ensemble_scenario_spec_vec)
    )

    outbreak_detect_spec_vec = map(
        threshold -> OutbreakDetectionSpecification(
            threshold,
            moving_avg_detection_lag,
            percent_visit_clinic,
            percent_clinic_tested,
            test_result_lag,
        ),
        detectthreshold_vec,
    )

    ensemble_scenario_spec_vec = map(
        outbreak_detect_spec -> ScenarioSpecification(
            ensemble_specification,
            outbreak_specification,
            noise_specification,
            outbreak_detect_spec,
            individual_test_specification,
        ),
        outbreak_detect_spec_vec,
    )

    accuracy_array = zeros(Int64, 2, length(detectthreshold_vec))

    @floop for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        ensemble_chars_file = get_ensemble_file(ensemble_scenario_spec)

        ensemble_chars_vec[i] = (
            OT_chars = ensemble_chars_file["OT_chars"],
            outbreak_detect_spec = ensemble_scenario_spec.outbreak_detection_specification,
            ind_test_spec = ensemble_scenario_spec.individual_test_specification,
        )

        accuracy_array[1, i] = float(
            ensemble_scenario_spec.outbreak_detection_specification.detection_threshold,
        )
        accuracy_array[2, i] = median(ensemble_chars_file["OT_chars"].accuracy)
    end

    optimal_accuracy = NaNMath.maximum(accuracy_array[2, :])
    optimal_threshold = accuracy_array[
        1, findfirst(==(optimal_accuracy), accuracy_array[2, :])
    ]
    return (threshold = optimal_threshold, accuracy = optimal_accuracy)
end
