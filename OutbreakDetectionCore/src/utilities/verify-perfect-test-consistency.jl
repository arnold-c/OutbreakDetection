export verify_perfect_test_consistency

function verify_perfect_test_consistency(results)
    perfect_test = IndividualTestSpecification(1.0, 1.0, 0)
    perfect_results = filter(r -> r.test_specification == perfect_test, results)

    if isempty(perfect_results)
        @warn "No perfect test scenarios found."
        return
    end

    # Group by all parameters that legitimately affect results
    # (excluding noise_level and noise_type_description which should NOT affect perfect test results)
    unique_groups = unique(
        [
            (
                    r.percent_tested,
                    r.accuracy_metric,
                    r.alert_method,
                    r.threshold_bounds,
                )
                for r in perfect_results
        ]
    )

    all_consistent = true

    for (pct, acc_metric, alert_method, threshold_bounds) in unique_groups
        subset = filter(
            r ->
            r.percent_tested == pct &&
                r.accuracy_metric == acc_metric &&
                r.alert_method == alert_method &&
                r.threshold_bounds == threshold_bounds,
            perfect_results
        )

        if isempty(subset)
            continue
        end

        first_res = subset[1]

        for i in eachindex(subset)
            if i == 1
                # skip the first as would be a identity comparison
                continue
            end
            res = subset[i]

            # Check optimal threshold
            if !isapprox(res.optimal_threshold, first_res.optimal_threshold; atol = 1.0e-6)
                @error "Mismatch in optimal_threshold" percent_tested = pct accuracy_metric = acc_metric alert_method = alert_method first =
                    first_res.optimal_threshold current = res.optimal_threshold first_noise_level =
                    first_res.noise_level first_noise_type = first_res.noise_type_description noise_level = res.noise_level noise_type = res.noise_type_description
                all_consistent = false
            end

            # Check accuracies (mean)
            if !isapprox(StatsBase.mean(res.accuracies), StatsBase.mean(first_res.accuracies); atol = 1.0e-6)
                @error "Mismatch in mean accuracy" percent_tested = pct accuracy_metric = acc_metric alert_method = alert_method first =
                    StatsBase.mean(first_res.accuracies) current = StatsBase.mean(res.accuracies) first_noise_level =
                    first_res.noise_level first_noise_type = first_res.noise_type_description noise_level = res.noise_level noise_type = res.noise_type_description
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
