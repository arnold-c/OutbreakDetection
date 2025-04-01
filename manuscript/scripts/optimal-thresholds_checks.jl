#%%
function compare_optimal_solution_extrema(
    optimal_solutions_vec,
    characteristic::Symbol,
    test_spec_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.85, 0.85, 0),
    ],
)
    @assert length(optimal_solutions_vec) == length(test_spec_vec)

    return mapreduce(
        ((optimal_solutions, test_spec),) -> extract_optimal_solution_extrema(
            optimal_solutions,
            characteristic,
            test_spec,
        ),
        vcat,
        zip(optimal_solutions_vec, test_spec_vec),
    )
end

function extract_optimal_solution_extrema(
    optimal_solutions,
    characteristic::Symbol,
    test_spec;
    digits = 1,
)
    test_extrema =
        round.(
            extrema(
                map(
                    test_chars ->
                        mean(vcat(getproperty(test_chars, characteristic)...)),
                    extract_test_optimal_solutions(
                        optimal_solutions, test_spec
                    ),
                ),
            );
            digits = digits,
        )
    return eval(:(Dict($test_spec => $test_extrema)))
end

function extract_test_optimal_solutions(optimal_solutions, test_spec)
    return filter(
        chars ->
            chars.individual_test_specification == test_spec,
        optimal_solutions,
    ).outbreak_threshold_chars
end

#%%
compare_optimal_solution_extrema(
    vcat(
        [perfect_test_optimal_solutions],
        repeat([dynamical_noise_rdt_optimal_solutions], 2),
    ), :detectiondelays
)

compare_optimal_solution_extrema(
    vcat(
        [perfect_test_optimal_solutions],
        repeat([dynamical_noise_rdt_optimal_solutions], 2),
    ), :alert_duration_vec
)
