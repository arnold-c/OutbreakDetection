using StyledStrings

#%%
function compare_optimal_solution_extrema(
    noise_optimal_solutions_vec,
    common_optimal_solutions_vec,
    characteristic_vec::AbstractVector{T},
    test_spec_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.85, 0.85, 0),
    ];
    digits = 1,
) where {T<:Symbol}
    num_rdt_tests =
        length(test_spec_vec) -
        sum(contains.(get_test_description.(test_spec_vec), "Perfect test"))

    for n in noise_optimal_solutions_vec
        println(styled"{red: $(get_noise_magnitude(n[1].noise_specification))}")

        compare_optimal_solution_extrema(
            vcat(
                common_optimal_solutions_vec,
                repeat([n], num_rdt_tests),
            ),
            characteristic_vec,
            test_spec_vec;
            digits = digits,
        )
    end
end

function compare_optimal_solution_extrema(
    optimal_solutions_vec,
    characteristic_vec::AbstractVector{T},
    test_spec_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.85, 0.85, 0),
    ];
    digits = 1,
) where {T<:Symbol}
    for c in characteristic_vec
        Base.display(
            compare_optimal_solution_extrema(
                optimal_solutions_vec,
                c,
                test_spec_vec;
                digits = digits,
            ),
        )
        println()
    end
end

function compare_optimal_solution_extrema(
    optimal_solutions_vec,
    characteristic::Symbol,
    test_spec_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.85, 0.85, 0),
    ];
    digits = 1,
)
    @assert length(optimal_solutions_vec) == length(test_spec_vec)

    println("Characteristic: $characteristic =======================")
    return mapreduce(
        ((optimal_solutions, test_spec),) -> extract_optimal_solution_extrema(
            optimal_solutions,
            characteristic,
            test_spec;
            digits = digits,
        ),
        merge,
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
    return Dict(test_spec => test_extrema)
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
    [
        poisson_noise_rdt_optimal_solutions,
        dynamical_noise_rdt_optimal_solutions,
    ],
    [perfect_test_optimal_solutions],
    [:detectiondelays, :alert_duration_vec, :accuracy];
    digits = 2,
)
