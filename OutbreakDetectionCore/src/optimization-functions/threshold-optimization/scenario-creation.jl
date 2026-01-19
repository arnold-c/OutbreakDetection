export create_scenarios_structvector

function create_scenarios_structvector(specification_vecs::ScenarioSpecificationVecs)
    UnPack.@unpack ensemble_specification_vec,
        noise_level_vec,
        noise_type_description_vec,
        test_specification_vec,
        percent_tested_vec,
        alert_method_vec,
        accuracy_metric_vec,
        threshold_bounds_vec,
        outbreak_specification_vec = specification_vecs

    combinations = Iterators.product(
        ensemble_specification_vec,
        noise_level_vec,
        noise_type_description_vec,
        test_specification_vec,
        percent_tested_vec,
        alert_method_vec,
        accuracy_metric_vec,
        threshold_bounds_vec,
        outbreak_specification_vec,
    )
    n_combinations = length(combinations)

    scenarios_vec = Vector{OptimizationScenario}(undef, n_combinations)

    for (
            i, (
                ensemble_spec,
                noise_level,
                noise_type_description,
                test_spec,
                percent_tested,
                alert_method,
                accuracy_metric,
                threshold_bounds,
                outbreak_spec,
            ),
        ) in enumerate(combinations)

        scenarios_vec[i] = OptimizationScenario(;
            ensemble_specification = ensemble_spec,
            noise_level = noise_level,
            noise_type_description = noise_type_description,
            test_specification = test_spec,
            percent_tested = percent_tested,
            alert_method = alert_method,
            accuracy_metric = accuracy_metric,
            threshold_bounds = threshold_bounds,
            outbreak_specification = outbreak_spec,
        )
    end

    return StructVector(scenarios_vec)
end
