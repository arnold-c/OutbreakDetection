export collect_threshold_char_vec

function collect_threshold_char_vec(percent_clinic_tested, base_parameters)
    UnPack.@unpack test_spec_vec,
        alertthreshold_vec,
        ensemble_specification,
        noise_specification,
        outbreak_specification,
        moving_avg_detection_lag,
        percent_visit_clinic,
        alertmethod = base_parameters

    non_clinical_case_test_spec_vec = filter(
        spec -> spec != CLINICAL_CASE_TEST_SPEC,
        test_spec_vec,
    )

    ensemble_scenario_spec_vec = Vector{ScenarioSpecification}(
        undef,
        length(non_clinical_case_test_spec_vec) * length(alertthreshold_vec) +
            length(alertthreshold_vec),
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
            alertmethod,
        ),
        alertthreshold_vec,
    )
    clinical_case_outbreak_detect_spec_vec = map(
        threshold -> OutbreakDetectionSpecification(
            threshold,
            moving_avg_detection_lag,
            percent_visit_clinic,
            1.0,
            alertmethod,
        ),
        alertthreshold_vec,
    )

    for (i, (ind_test_spec, outbreak_detect_spec)) in enumerate(
            Iterators.product(
                non_clinical_case_test_spec_vec, outbreak_detect_spec_vec
            ),
        )
        ensemble_scenario_spec = ScenarioSpecification(
            ensemble_specification,
            outbreak_specification,
            noise_specification,
            outbreak_detect_spec,
            ind_test_spec,
        )

        ensemble_scenario_spec_vec[i] = ensemble_scenario_spec
    end

    ensemble_scenario_spec_vec[(end - length(clinical_case_outbreak_detect_spec_vec) + 1):end] .= create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_specification],
            [outbreak_specification],
            [noise_specification],
            clinical_case_outbreak_detect_spec_vec,
            # TODO: update this to calculate for all detection thresholds
            [CLINICAL_CASE_TEST_SPEC],
        ),
    )

    FLoops.@floop for (i, ensemble_scenario_spec) in
        pairs(ensemble_scenario_spec_vec)
        ensemble_chars_file = get_ensemble_file(ensemble_scenario_spec)

        ensemble_chars_vec[i] = (
            OT_chars = ensemble_chars_file["OT_chars"],
            outbreak_detect_spec = ensemble_scenario_spec.outbreak_detection_specification,
            ind_test_spec = ensemble_scenario_spec.individual_test_specification,
            noise_specification = noise_specification,
        )
    end

    return sort!(
        ensemble_chars_vec;
        by = x -> (
            x.outbreak_detect_spec.alert_threshold,
            x.ind_test_spec.specificity,
            x.ind_test_spec.test_result_lag,
        ),
    )
end
