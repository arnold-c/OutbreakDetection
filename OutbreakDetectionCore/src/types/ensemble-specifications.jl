export EnsembleSpecification

"""
    EnsembleSpecification

Complete specification for ensemble disease simulation parameters and configuration.

This struct defines all parameters needed to run an ensemble of SEIR disease simulations,
including initial population states, time parameters, disease dynamics for both emergent
and null scenarios, noise characteristics, and output configuration. It serves as the
primary configuration object for ensemble-based epidemiological modeling and analysis.

The `dirpath` field is automatically constructed from the component specifications to
create a unique filesystem path for storing ensemble simulation results, incorporating
all relevant parameter values in a hierarchical directory structure.

# Fields

  - `label::String`: The disease name
  - `state_parameters::StateParameters`: Initial population state configuration (S, E, I, R compartments)
  - `time_parameters::SimTimeParameters`: Simulation time configuration including burnin, duration, and timestep
  - `dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for outbreak scenarios
  - `null_dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for null/baseline scenarios
  - `dynamical_noise_params::DynamicalNoiseParameters`: Specification for dynamical noise in disease parameters
  - `nsims::Int64`: Number of simulation replicates in the ensemble
  - `dirpath::String`: Automatically generated directory path for storing ensemble results

# Constructor

    EnsembleSpecification(; label, state_parameters, time_parameters, dynamics_parameter_specification,
                         null_dynamics_parameter_specification, dynamical_noise_params, nsims, dirpath)

# Example

```julia
# Create ensemble specification with all components
ensemble_spec = EnsembleSpecification(;
    state_parameters = StateParameters(;
        N = 500_000, s_prop = 0.9, e_prop = 0.05, i_prop = 0.05
    ),
    time_parameters = SimTimeParameters(;
        burnin = Year(10), tmax = 20.0, tstep = 0.1
    ),
    dynamics_parameter_specification = emergent_dynamics,
    null_dynamics_parameter_specification = null_dynamics,
    dynamical_noise_params = noise_spec,
    nsims = 1000,
    dirpath = "/path/to/results",
)

# Access ensemble configuration
println("Running \$(ensemble_spec.nsims) simulations")
println("Results stored in: \$(ensemble_spec.dirpath)")
```

# See Also

  - [`StateParameters`](@ref): Initial population state configuration
  - [`SimTimeParameters`](@ref): Simulation time parameters
  - [`DynamicsParameterSpecification`](@ref): Disease dynamics specifications
  - [`DynamicalNoiseParameters`](@ref): Dynamical noise parameters
  - [`simulate_ensemble_seir_results`](@ref): Function that uses EnsembleSpecification for simulation
"""
AutoHashEquals.@auto_hash_equals Base.@kwdef struct EnsembleSpecification
    label::String
    state_parameters::StateParameters
    time_parameters::SimTimeParameters
    dynamics_parameter_specification::DynamicsParameterSpecification
    dynamical_noise_params::DynamicalNoiseParameters
    nsims::Int64
    dirpath::String
end

"""
    EnsembleSpecification(label, state_parameters, time_parameters, dynamics_parameter_specification,
                         null_dynamics_parameter_specification, dynamical_noise_params, nsims)

Construct an EnsembleSpecification with automatic directory path generation.

This constructor creates a complete ensemble specification by combining all the individual
component specifications and automatically generating a hierarchical directory path that
uniquely identifies this ensemble configuration. The directory path incorporates all
relevant parameter values to ensure unique storage locations for different ensemble runs.

The generated directory structure follows this hierarchy:

  - Base: "ensemble/seasonal-infectivity-import/tau-leaping"
  - Population parameters: N, initial recovered proportion
  - Simulation parameters: number of simulations, time parameters
  - Disease dynamics: R₀, latent period, infectious period for both emergent and null scenarios
  - Noise parameters: R₀ noise, latent period noise, infectious period noise
  - Vaccination parameters: burnin and post-burnin coverage ranges
  - Birth rate and transmission parameters

# Arguments

  - `label::String`: The disease name
  - `state_parameters::StateParameters`: Initial population state configuration
  - `time_parameters::SimTimeParameters`: Simulation time configuration
  - `dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for outbreak scenarios
  - `null_dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for baseline scenarios
  - `dynamical_noise_params::DynamicalNoiseParameters`: Dynamical noise specification
  - `nsims::Int64`: Number of simulation replicates

# Returns

  - `EnsembleSpecification`: Complete ensemble specification with auto-generated directory path

# Example

```julia
# Create ensemble specification with automatic path generation
ensemble_spec = EnsembleSpecification(
    label,
    state_params,
    time_params,
    emergent_dynamics,
    null_dynamics,
    noise_spec,
    1000,  # number of simulations
)

# The dirpath is automatically generated based on all parameters
println(ensemble_spec.dirpath)
# Output: "ensemble/seasonal-infectivity-import/tau-leaping/N_500000/r_0.95/nsims_1000/..."
```
"""
function EnsembleSpecification(
        label::String,
        state_parameters::StateParameters,
        time_parameters::SimTimeParameters,
        dynamics_parameter_specification::DynamicsParameterSpecification,
        dynamical_noise_params::DynamicalNoiseParameters,
        nsims::Int64,
    )
    dirpath = outdir(
        "ensemble",
        "seasonal-infectivity-import",
        "tau-leaping",
        "$label",
        "N_$(state_parameters.init_states.N)",
        "r_$(state_parameters.init_state_props.r_prop)",
        "nsims_$(nsims)",
        "R0_$(dynamics_parameter_specification.R_0)",
        "latent_period_$(round(1 / dynamics_parameter_specification.sigma; digits = 2))",
        "infectious_period_$(round(1 / dynamics_parameter_specification.gamma; digits = 2))",
        "noise_R0_$(dynamical_noise_params.R_0)",
        "noise_latent_period_$(round(dynamical_noise_params.latent_period; digits = 2))",
        "noise_infectious_period_$(round(dynamical_noise_params.infectious_duration; digits = 2))",
        "min_vaccination_coverage_$(dynamics_parameter_specification.min_vaccination_coverage)",
        "max_vaccination_coverage_$(dynamics_parameter_specification.max_vaccination_coverage)",
        "births_per_k_$(dynamics_parameter_specification.annual_births_per_k)",
        "beta_force_$(dynamics_parameter_specification.beta_force)",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)",
    )

    return EnsembleSpecification(
        label,
        state_parameters,
        time_parameters,
        dynamics_parameter_specification,
        dynamical_noise_params,
        nsims,
        dirpath,
    )
end
