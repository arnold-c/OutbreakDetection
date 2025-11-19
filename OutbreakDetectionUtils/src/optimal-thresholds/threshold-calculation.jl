using UnPack: UnPack
using StatsBase: StatsBase
using ProgressMeter: ProgressMeter
using StructArrays: StructArrays

export calculate_optimal_threshold, calculate_OptimalThresholdCharacteristics

"""
    calculate_OptimalThresholdCharacteristics(percent_clinic_tested_vec, 
                                              ind_test_spec_vec, 
                                              base_parameters)

Calculate optimal threshold characteristics for multiple test specifications.
"""
function calculate_OptimalThresholdCharacteristics(
        percent_clinic_tested_vec,
        ind_test_spec_vec,
        base_parameters,
    )
    non_clinical_case_test_spec_vec = filter(
        spec -> !(spec in CLINICAL_TEST_SPECS),
        ind_test_spec_vec,
    )

    non_clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics,
    }(
        undef,
        length(percent_clinic_tested_vec) *
            length(non_clinical_case_test_spec_vec),
    )

    ProgressMeter.@showprogress for (
            i, (percent_clinic_tested, ind_test_spec),
        ) in enumerate(
            Iterators.product(
                percent_clinic_tested_vec, non_clinical_case_test_spec_vec
            ),
        )
        non_clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            percent_clinic_tested,
            ind_test_spec,
            base_parameters,
        )
    end

    clinical_case_test_spec_vec = filter(
        spec -> spec in CLINICAL_TEST_SPECS,
        ind_test_spec_vec,
    )

    if length(clinical_case_test_spec_vec) == 0
        return StructArrays.StructArray(
            non_clinical_case_optimal_thresholds_vec
        )
    end

    clinical_case_optimal_thresholds_vec = Vector{
        OptimalThresholdCharacteristics,
    }(
        undef,
        length(clinical_case_test_spec_vec),
    )

    for (i, ind_test_spec) in enumerate(clinical_case_test_spec_vec)
        clinical_case_optimal_thresholds_vec[i] = calculate_optimal_threshold(
            1.0,
            ind_test_spec,
            base_parameters,
        )
    end

    return StructArrays.StructArray(
        vcat(
            non_clinical_case_optimal_thresholds_vec,
            clinical_case_optimal_thresholds_vec,
        ),
    )
end

"""
    calculate_optimal_threshold(percent_clinic_tested, 
                                individual_test_specification, 
                                base_parameters)

Calculate optimal alert threshold for a given test specification.
"""
function calculate_optimal_threshold(
        percent_clinic_tested,
        individual_test_specification,
        base_parameters,
    )
    UnPack.@unpack alertthreshold_vec,
        ensemble_specification,
        noise_specification,
        outbreak_specification,
        moving_avg_detection_lag,
        percent_visit_clinic,
        alertmethod = base_parameters

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
                alertmethod,
            ),
            individual_test_specification,
        ),
        alertthreshold_vec,
    )

    optimal_accuracy = 0.0
    optimal_threshold = 0
    optimal_OT_chars = 0

    for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        scenario_chars_file = get_ensemble_file(ensemble_scenario_spec)
        OT_chars = scenario_chars_file["OT_chars"]
        accuracy = StatsBase.median(OT_chars.accuracy)
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
        noise_specification,
        percent_clinic_tested,
        optimal_threshold,
        optimal_accuracy,
    )
end
