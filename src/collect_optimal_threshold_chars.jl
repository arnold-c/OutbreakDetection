function collect_OptimalThresholdCharacteristics(
        noise_spec_vec,
        ensemble_percent_clinic_tested_vec,
        optimal_threshold_test_spec_vec,
        optimal_threshold_core_params;
        clinical_hline = false,
    )
    noise_descriptions = OutbreakDetectionUtils.get_noise_description.(noise_spec_vec)
    unique_noise_descriptions = unique(noise_descriptions)

    shape_noise_specifications =
        map(unique_noise_descriptions) do noise_description
        filter(
            noise_spec ->
            noise_description == OutbreakDetectionUtils.get_noise_description(noise_spec),
            noise_spec_vec,
        )
    end

    @assert length(vcat(shape_noise_specifications...)) ==
        length(unique_noise_descriptions) *
        length(shape_noise_specifications[1])

    if clinical_hline
        optimal_threshold_test_spec_vec = vcat(
            optimal_threshold_test_spec_vec, CLINICAL_CASE_TEST_SPEC
        )
    end

    optimal_thresholds_vecs = Array{
        StructArray{OutbreakDetectionUtils.OptimalThresholdCharacteristics},
    }(
        undef,
        length(unique_noise_descriptions),
        length(shape_noise_specifications[1]),
    )

    for (i, noise_description) in pairs(unique_noise_descriptions)
        label_noise_description = Match.@match noise_description begin
            "poisson" => "Poisson Noise"
            "dynamical, in-phase" => "Dynamical Noise: In-Phase"
            _ => "Other Noise"
        end

        for (j, noise_spec) in pairs(shape_noise_specifications[i])
            optimal_threshold_comparison_params = (
                noise_specification = noise_spec,
                optimal_threshold_core_params...,
            )

            optimal_thresholds_vecs[i, j] = OutbreakDetectionUtils.calculate_OptimalThresholdCharacteristics(
                ensemble_percent_clinic_tested_vec,
                optimal_threshold_test_spec_vec,
                optimal_threshold_comparison_params,
            )
        end
    end
    return optimal_thresholds_vecs
end
