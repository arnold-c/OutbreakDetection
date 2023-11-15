using StructArrays

function calculate_OptimalThresholdCharacteristics(
    percent_clinic_tested_vec,
    ind_test_spec_vec,
    base_parameters
)
    clinical_case_test_spec = IndividualTestSpecification(1.0, 0.0)

    non_clinical_case_test_spec_vec = filter(
        spec -> spec != clinical_case_test_spec,
        ind_test_spec_vec
    )
    non_clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics
    }(
        undef,
        length(percent_clinic_tested_vec) *
        length(non_clinical_case_test_spec_vec),
    )

    @showprogress for (i, (percent_clinic_tested, ind_test_spec)) in enumerate(
        Iterators.product(
            percent_clinic_tested_vec, non_clinical_case_test_spec_vec
        ),
    )
        non_clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            percent_clinic_tested,
            ind_test_spec,
            base_parameters
        )
    end

    clinical_case_optimal_thresholds_vec = calculate_optimal_threshold(
        1.0,
        clinical_case_test_spec,
        base_parameters
    )

    return StructArray(
        vcat(
            non_clinical_case_optimal_thresholds_vec,
            clinical_case_optimal_thresholds_vec,
        ),
    )
end

function calculate_optimal_threshold(
    percent_clinic_tested,
    individual_test_specification,
    base_parameters
)
    @unpack alertthreshold_vec,
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
        alertthreshold_vec,
    )

    optimal_accuracy = 0.0
    optimal_threshold = 0
    optimal_OT_chars = 0

    for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        OT_chars = get_ensemble_file(ensemble_scenario_spec)["OT_chars"]
        accuracy = median(OT_chars.accuracy)
        alert_threshold =
            ensemble_scenario_spec.outbreak_detection_specification.alert_threshold

        if i == 1
            optimal_accuracy = accuracy
            optimal_threshold = alert_threshold
            optimal_OT_chars = OT_chars
            continue
        end
        if !isnan(accuracy) && accuracy > optimal_accuracy
            optimal_accuracy = accuracy
            optimal_threshold = alert_threshold
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
