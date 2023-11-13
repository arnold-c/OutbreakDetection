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

        if optimal_thresholds_vec[i].detection_threshold !=
            optimal_thresholds_vec[i].scenario_specification.outbreak_detection_specification.detection_threshold
            @error "Warning. The optimal threshold is not the same as the detection threshold for the i = $i scenario specification."
        end
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
        OT_chars = get_ensemble_file(ensemble_scenario_spec)["OT_chars"]
        accuracy = median(OT_chars.accuracy)
        detection_threshold =
            ensemble_scenario_spec.outbreak_detection_specification.detection_threshold

        if i == 1
            optimal_accuracy = accuracy
            optimal_threshold = detection_threshold
            optimal_OT_chars = OT_chars
            continue
        end
        if !isnan(accuracy) && accuracy > optimal_accuracy
            optimal_accuracy = accuracy
            optimal_threshold = detection_threshold
            optimal_OT_chars = OT_chars
        end
    end

    return OptimalThresholdCharacteristics(
        optimal_OT_chars,
        individual_test_specification,
        percent_clinic_tested,
        optimal_threshold,
        optimal_accuracy,
    )
end
