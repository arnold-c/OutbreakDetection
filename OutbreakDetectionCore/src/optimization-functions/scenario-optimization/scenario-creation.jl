export calculate_combination_to_run, create_scenario_grid, _create_noise_for_level

"""
    calculate_combination_to_run(ensemble_specifications, outbreak_specifications, 
                                 noise_specifications, outbreak_detection_specifications, 
                                 individual_test_specifications, optim_method, 
                                 accuracy_functions; scenario_parameter_symbols=...)

Calculate all combinations of scenarios to run for optimization.
"""
function calculate_combination_to_run(
        ensemble_specifications,
        outbreak_specifications,
        noise_specifications,
        outbreak_detection_specifications,
        individual_test_specifications,
        optim_method::TMethod = MSO,
        accuracy_functions = [arithmetic_mean, calculate_f_beta_score];
        scenario_parameter_symbols = [
            :ensemble_spec,
            :outbreak_spec,
            :noise_spec,
            :outbreak_detection_spec,
            :test_spec,
            :optimization_method,
            :accuracy_function,
        ],
    ) where {TMethod <: Type{<:OptimizationMethods}}
    @assert mapreduce(
        f -> in(f, [arithmetic_mean, calculate_f_beta_score]),
        +,
        unique(accuracy_functions),
    ) == length(unique(accuracy_functions))

    return DataFrames.DataFrame(
        Iterators.product(
            ensemble_specifications,
            outbreak_specifications,
            noise_specifications,
            outbreak_detection_specifications,
            individual_test_specifications,
            [optim_method],
            accuracy_functions,
        ),
        scenario_parameter_symbols,
    )
end

"""
    create_scenario_grid(; kwargs...)

Create a grid of scenario parameters for optimization.

Generates all combinations of provided parameter vectors.

# Keyword Arguments
- `R_0_values::Vector{Float64}`: R_0 values to test
- `noise_levels::Vector{String}`: Noise levels ("low", "medium", "high")
- `test_sensitivities::Vector{Float64}`: Test sensitivities to test
- `percent_tested_values::Vector{Float64}`: Proportions tested to test
- `outbreak_thresholds::Vector{Float64}`: Outbreak thresholds to test
- `base_target::TargetDiseaseDynamicsParameters`: Base target disease parameters
- `base_noise::DynamicalNoiseSpecification`: Base noise specification
- `base_test::IndividualTestSpecification`: Base test specification

# Returns
- `Vector{ScenarioParameters}`: All parameter combinations

# Examples
```julia
scenarios = create_scenario_grid(
    R_0_values = [12.0, 16.0, 20.0],
    noise_levels = ["low", "medium", "high"],
    test_sensitivities = [0.8, 0.9, 1.0],
    percent_tested_values = [0.3, 0.5, 0.7],
    outbreak_thresholds = [0.01, 0.02, 0.05],
    base_target = target_params,
    base_noise = noise_spec,
    base_test = test_spec
)

# Results in 3 × 3 × 3 × 3 × 3 = 243 scenarios
length(scenarios)  # 243
```
"""
function create_scenario_grid(;
        R_0_values::Vector{Float64},
        noise_levels::Vector{String},
        test_sensitivities::Vector{Float64},
        percent_tested_values::Vector{Float64},
        outbreak_thresholds::Vector{Float64},
        base_target::TargetDiseaseDynamicsParameters,
        base_noise::DynamicalNoiseSpecification,
        base_test::IndividualTestSpecification,
    )
    scenarios = ScenarioParameters[]

    for R_0 in R_0_values
        target = TargetDiseaseDynamicsParameters(
            R_0 = R_0,
            latent_period_days = base_target.latent_period_days,
            infectious_duration_days = base_target.infectious_duration_days,
            beta_force = base_target.beta_force,
            seasonality = base_target.seasonality,
            life_expectancy_years = base_target.life_expectancy_years,
            population_N = base_target.population_N,
        )

        for noise_level in noise_levels
            # Map noise level to parameters
            noise = _create_noise_for_level(noise_level, base_noise)

            for sensitivity in test_sensitivities
                test = IndividualTestSpecification(
                    sensitivity = sensitivity,
                    specificity = base_test.specificity,
                    test_result_lag = base_test.test_result_lag,
                )

                for percent_tested in percent_tested_values
                    for threshold in outbreak_thresholds
                        scenario = ScenarioParameters(
                            target_dynamics = target,
                            noise_spec = noise,
                            test_spec = test,
                            percent_clinic_tested = percent_tested,
                            outbreak_threshold = threshold,
                        )
                        push!(scenarios, scenario)
                    end
                end
            end
        end
    end

    return scenarios
end

"""
    _create_noise_for_level(level, base_noise)

Internal function to create noise specification for a given level.

# Arguments
- `level::String`: Noise level ("low", "medium", "high")
- `base_noise::DynamicalNoiseSpecification`: Base noise specification

# Returns
- `DynamicalNoiseSpecification`: Noise specification for the level
"""
function _create_noise_for_level(level::String, base_noise::DynamicalNoiseSpecification)
    # Map noise levels to vaccination coverage ranges
    # These will be optimized to match target noise levels
    vaccination_bounds = if level == "low"
        [0.0, 0.5]  # Lower vaccination = higher noise
    elseif level == "medium"
        [0.3, 0.7]
    elseif level == "high"
        [0.5, 0.9]  # Higher vaccination = lower noise
    else
        error("Unknown noise level: $level")
    end

    return DynamicalNoiseSpecification(
        R_0 = base_noise.R_0,
        latent_period = base_noise.latent_period,
        duration_infection = base_noise.duration_infection,
        correlation = base_noise.correlation,
        poisson_component = base_noise.poisson_component,
        vaccination_bounds = vaccination_bounds,
    )
end
