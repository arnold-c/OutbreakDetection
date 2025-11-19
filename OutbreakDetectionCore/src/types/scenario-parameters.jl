export ScenarioParameters, ThresholdOptimizationScenario

using AutoHashEquals: AutoHashEquals

"""
    ScenarioParameters

Base parameters for a scenario with automatic hashing and equality.

Uses AutoHashEquals for efficient comparison and caching.

# Fields
- `target_dynamics::TargetDiseaseDynamicsParameters`: Target disease parameters
- `noise_spec::NoiseSpecification`: Noise specification (DynamicalNoiseSpecification or PoissonNoise)
- `test_spec::IndividualTestSpecification`: Testing specification
- `percent_clinic_tested::Float64`: Proportion tested at clinic (0-1)
- `outbreak_threshold::Float64`: Threshold for outbreak detection

# Examples
```julia
target = TargetDiseaseDynamicsParameters(
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.2
)

noise = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 1.0
)

test = IndividualTestSpecification(
    sensitivity = 1.0,
    specificity = 0.8,
    test_result_lag = 0
)

scenario = ScenarioParameters(
    target_dynamics = target,
    noise_spec = noise,
    test_spec = test,
    percent_clinic_tested = 0.5,
    outbreak_threshold = 0.01
)

# Automatic hashing and equality
scenario1 == scenario2  # Uses AutoHashEquals
hash(scenario)  # Automatic hash function
```
"""
AutoHashEquals.@auto_hash_equals struct ScenarioParameters
    target_dynamics::TargetDiseaseDynamicsParameters
    noise_spec::NoiseSpecification
    test_spec::IndividualTestSpecification
    percent_clinic_tested::Float64
    outbreak_threshold::Float64

    function ScenarioParameters(
            target_dynamics,
            noise_spec,
            test_spec,
            percent_clinic_tested,
            outbreak_threshold,
        )
        @assert 0.0 <= percent_clinic_tested <= 1.0 "Percent tested must be in [0, 1]"
        @assert outbreak_threshold > 0.0 "Outbreak threshold must be positive"

        return new(
            target_dynamics,
            noise_spec,
            test_spec,
            percent_clinic_tested,
            outbreak_threshold,
        )
    end
end

"""
    ThresholdOptimizationScenario

Complete scenario specification for threshold optimization.

Combines scenario parameters with ensemble and time specifications.

# Fields
- `scenario_params::ScenarioParameters`: Core scenario parameters
- `ensemble_spec::EnsembleSpecification`: Ensemble simulation parameters
- `time_params::SimTimeParameters`: Time parameters
- `nsims::Int64`: Number of simulations
- `seed::Int64`: Random seed for reproducibility

# Examples
```julia
scenario = ThresholdOptimizationScenario(
    scenario_params = scenario_params,
    ensemble_spec = ensemble_spec,
    time_params = time_params,
    nsims = 1000,
    seed = 1234
)
```
"""
AutoHashEquals.@auto_hash_equals struct ThresholdOptimizationScenario
    scenario_params::ScenarioParameters
    ensemble_spec::EnsembleSpecification
    time_params::SimTimeParameters
    nsims::Int64
    seed::Int64

    function ThresholdOptimizationScenario(
            scenario_params, ensemble_spec, time_params, nsims, seed
        )
        @assert nsims > 0 "Number of simulations must be positive"
        @assert seed >= 0 "Seed must be non-negative"

        return new(scenario_params, ensemble_spec, time_params, nsims, seed)
    end
end
