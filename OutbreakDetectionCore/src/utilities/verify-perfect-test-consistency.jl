export verify_perfect_test_consistency

function verify_perfect_test_consistency(results)
    perfect_test = IndividualTestSpecification(1.0, 1.0, 0)
    perfect_results = filter(r -> r.test_specification == perfect_test, results)

    if isempty(perfect_results)
        @warn "No perfect test scenarios found."
        return
    end

    unique_percents = unique(perfect_results.percent_tested)
    all_consistent = true

    for pct in unique_percents
        subset = filter(r -> r.percent_tested == pct, perfect_results)
        if isempty(subset)
            continue
        end

        first_res = subset[1]

        for i in 2:length(subset)
            res = subset[i]

            # Check optimal threshold
            if !isapprox(res.optimal_threshold, first_res.optimal_threshold; atol = 1.0e-6)
                @error "Mismatch in optimal_threshold for percent_tested=$pct" first =
                    first_res.optimal_threshold current = res.optimal_threshold first_noise_level =
                    first_res.noise_level noise_level = res.noise_level noise_type = res.noise_type_description
                all_consistent = false
            end

            # Check accuracies (mean)
            if !isapprox(StatsBase.mean(res.accuracies), StatsBase.mean(first_res.accuracies); atol = 1.0e-6)
                @error "Mismatch in mean accuracy for percent_tested=$pct" first =
                    StatsBase.mean(first_res.accuracies) current = StatsBase.mean(res.accuracies) first_noise_level =
                    first_res.noise_level noise_level = res.noise_level noise_type = res.noise_type_description
                all_consistent = false
            end
        end
    end

    return if all_consistent
        @info "Verification passed: Perfect test scenarios produce identical metrics across noise types/levels."
    else
        @error "Verification failed."
    end
end
