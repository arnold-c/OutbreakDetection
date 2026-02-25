export compare_optimal_solution_mean_extrema

function compare_optimal_solution_mean_extrema(
        common_optimal_solutions_vec,
        noise_optimal_solutions_vec,
        characteristic_vec::AbstractVector{T},
        test_spec_vec = [
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.9, 0.9, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.85, 0.85, 0),
        ];
        digits = 1,
        color = :red,
    ) where {T <: Symbol}
    num_perfect_test = sum(
        contains.(OutbreakDetectionCore.get_test_description.(test_spec_vec), "Perfect test")
    )

    num_rdt_tests = length(test_spec_vec) - num_perfect_test

    for n in noise_optimal_solutions_vec
        @assert !isempty(n) "Noise-specific optimization results cannot be empty"
        println(
            styled"{$color: Noise level: $(round(n[1].noise_level; digits = digits))}"
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
    return
end

function compare_optimal_solution_mean_extrema(
        optimal_solutions_vec,
        characteristic_vec::AbstractVector{T},
        test_spec_vec = [
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.9, 0.9, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.85, 0.85, 0),
        ];
        digits = 1,
    ) where {T <: Symbol}
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
    return
end

function compare_optimal_solution_mean_extrema(
        optimal_solutions_vec::AbstractVector{<:StructVector{OutbreakDetectionCore.OptimizationResult}},
        characteristic::Symbol,
        test_spec_vec::AbstractVector{T} = [
            OutbreakDetectionCore.IndividualTestSpecification(1.0, 1.0, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.9, 0.9, 0),
            OutbreakDetectionCore.IndividualTestSpecification(0.85, 0.85, 0),
        ];
        digits = 1,
    ) where {T <: OutbreakDetectionCore.IndividualTestSpecification}
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
        optimal_solutions::StructVector{OutbreakDetectionCore.OptimizationResult},
        characteristic::Symbol,
        test_spec::T;
        digits = 1,
    ) where {T <: OutbreakDetectionCore.IndividualTestSpecification}
    test_vals = map(
        result -> extract_characteristic_mean(result, characteristic),
        extract_test_optimal_solutions(optimal_solutions, test_spec),
    )

    @assert !isempty(test_vals) "No optimization results found for $test_spec"

    test_extrema = round.(extrema(test_vals); digits = digits)
    test_mean = round.(StatsBase.mean(test_vals); digits = digits)
    test_median = round.(StatsBase.median(test_vals); digits = digits)
    return Dict(
        test_spec =>
            (; extrema = test_extrema, mean = test_mean, median = test_median),
    )
end

function extract_test_optimal_solutions(optimal_solutions, test_spec)
    return filter(
        chars ->
        chars.test_specification == test_spec,
        optimal_solutions,
    )
end

function extract_characteristic_mean(
        result::OutbreakDetectionCore.OptimizationResult,
        characteristic::Symbol,
    )
    vals = getproperty(result, characteristic)
    return StatsBase.mean(flatten_characteristic_values(vals))
end

flatten_characteristic_values(vals::AbstractVector{<:Number}) = vals

function flatten_characteristic_values(vals::AbstractVector{<:AbstractVector})
    if isempty(vals)
        return Int[]
    end
    return reduce(vcat, vals)
end
