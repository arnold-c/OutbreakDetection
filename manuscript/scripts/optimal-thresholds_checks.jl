using StyledStrings
using StructArrays

#%%
function compare_optimal_solution_mean_extrema(
    common_optimal_solutions_vec,
    noise_optimal_solutions_vec,
    characteristic_vec::AbstractVector{T},
    test_spec_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.85, 0.85, 0),
    ];
    digits = 1,
    color = :red,
) where {T<:Symbol}
    num_perfect_test = sum(
        contains.(get_test_description.(test_spec_vec), "Perfect test")
    )

    num_rdt_tests = length(test_spec_vec) - num_perfect_test

    for n in noise_optimal_solutions_vec
        println(
            styled"{$color: $(get_noise_magnitude(n[1].noise_specification))}"
        )

        optimal_solutions_vec = vcat(
            repeat(common_optimal_solutions_vec, num_perfect_test),
            repeat([n], num_rdt_tests),
        )

        @assert optimal_solutions_vec[1] ==
            optimal_solutions_vec[num_perfect_test]

        @assert optimal_solutions_vec[num_perfect_test + 1] ==
            optimal_solutions_vec[end]

        compare_optimal_solution_mean_extrema(
            optimal_solutions_vec,
            characteristic_vec,
            test_spec_vec;
            digits = digits,
        )
    end
end

#%%
function compare_optimal_solution_mean_extrema(
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
        sol = compare_optimal_solution_mean_extrema(
            optimal_solutions_vec,
            c,
            test_spec_vec;
            digits = digits,
        )
        Base.display(sol)
        println()
    end
end

#%%
function compare_optimal_solution_mean_extrema(
    optimal_solutions_vec::AbstractVector{S},
    characteristic::Symbol,
    test_spec_vec::AbstractVector{T} = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.85, 0.85, 0),
    ];
    digits = 1,
) where {S<:StructArrays.StructVector,T<:IndividualTestSpecification}
    @assert length(optimal_solutions_vec) == length(test_spec_vec) "\nlength(optimal_solutions_vec) = $(length(optimal_solutions_vec))\nlength(test_spec_vec) = $(length(test_spec_vec)) "

    println("Characteristic: $characteristic =======================")
    return mapreduce(
        ((optimal_solutions, test_spec),) ->
            compare_optimal_solution_mean_extrema(
                optimal_solutions,
                characteristic,
                test_spec;
                digits = digits,
            ),
        merge,
        zip(optimal_solutions_vec, test_spec_vec),
    )
end

#%%
function compare_optimal_solution_mean_extrema(
    optimal_solutions::StructArrays.StructVector{S},
    characteristic::Symbol,
    test_spec::T;
    digits = 1,
) where {S<:OptimalThresholdCharacteristics,T<:IndividualTestSpecification}
    test_vals = map(
        test_chars -> mean(vcat(getproperty(test_chars, characteristic)...)),
        extract_test_optimal_solutions(optimal_solutions, test_spec),
    )

    test_extrema = round.(extrema(test_vals); digits = digits)
    test_mean = round.(mean(test_vals); digits = digits)
    test_median = round.(median(test_vals); digits = digits)
    return Dict(
        test_spec =>
            (; extrema = test_extrema, mean = test_mean, median = test_median),
    )
end

function extract_test_optimal_solutions(optimal_solutions, test_spec)
    test_optimal_solutions = filter(
        chars ->
            chars.individual_test_specification == test_spec,
        optimal_solutions,
    )
    return test_optimal_solutions.outbreak_threshold_chars
end

#%%
all_perfect_test_optimal_solutions =
    subset(
        optim_df,
        :test_spec => ByRow(
            in([IndividualTestSpecification(1.0, 1.0, 14),
                IndividualTestSpecification(1.0, 1.0, 0),
            ]),
        ),
        :noise_spec => ByRow(t -> isa(t, PoissonNoiseSpecification)),
    ) |>
    vec ∘ reshape_optim_df_to_matrix;

all_static_noise_rdt_test_optimal_solutions =
    subset(
        optim_df,
        :test_spec => ByRow(
            in([
                IndividualTestSpecification(0.9, 0.9, 0),
                IndividualTestSpecification(0.85, 0.85, 0),
            ]),
        ),
        :noise_spec => ByRow(t -> isa(t, PoissonNoiseSpecification)),
    ) |>
    vec ∘ reshape_optim_df_to_matrix;

all_dynamical_noise_rdt_test_optimal_solutions =
    subset(
        optim_df,
        :test_spec => ByRow(
            in([
                IndividualTestSpecification(0.9, 0.9, 0),
                IndividualTestSpecification(0.85, 0.85, 0),
            ]),
        ),
        :noise_spec => ByRow(t -> isa(t, DynamicalNoiseSpecification)),
    ) |>
    vec ∘ reshape_optim_df_to_matrix;

#%%
for (i, color) in zip(
    eachindex(all_perfect_test_optimal_solutions),
    [:red, :blue, :green, :magenta, :cyan],
)
    compare_optimal_solution_mean_extrema(
        [all_perfect_test_optimal_solutions[i]],
        [
            all_static_noise_rdt_test_optimal_solutions[i],
            all_dynamical_noise_rdt_test_optimal_solutions[i],
        ],
        [
            # :detectiondelays,
            # :alert_duration_vec,
            :accuracy
        ],
        [
            IndividualTestSpecification(1.0, 1.0, 14),
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
            IndividualTestSpecification(0.85, 0.85, 0),
        ];
        digits = 2,
        color = color,
    )
end
