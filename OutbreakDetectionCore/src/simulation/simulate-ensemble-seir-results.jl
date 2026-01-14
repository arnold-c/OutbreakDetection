export generate_single_ensemble

"""
    generate_single_ensemble(ensemble_spec; seed = 1234)

Generate a single ensemble of SEIR simulations.

This function creates an ensemble of stochastic SEIR model runs, simulating
dynamics scenario (with specified vaccination coverage).
Each simulation in the ensemble uses the same seasonal transmission pattern
but different random seeds to capture stochastic variation.

# Arguments
- `ensemble_spec::EnsembleSpecification`: Specification containing all
  parameters for the ensemble simulation, including state parameters, time
  parameters, dynamics specifications, and the number of simulations to run.
- `seed::Int64 = 1234`: Base random seed for reproducibility. Each simulation
  in the ensemble uses `seed + (sim - 1)` to ensure different but reproducible
  random trajectories.

# Returns
- `StructVector{SEIRRun}`: A struct vector containing `nsims` individual SEIR simulation
  results.
"""
function generate_single_ensemble(
        ensemble_spec::EnsembleSpecification;
        seed::Int64 = 1234
    )
    UnPack.@unpack state_parameters,
        dynamics_parameter_specification,
        time_parameters,
        nsims = ensemble_spec

    UnPack.@unpack tlength = time_parameters

    init_states_sv = StaticArrays.SVector(state_parameters.init_states)

    # Get concrete type to avoid abstract element types
    seir_results = Vector{SEIRRun}(undef, nsims)

    beta_vec = zeros(Float64, tlength)
    calculate_beta_amp!(
        beta_vec,
        dynamics_parameter_specification,
        time_parameters
    )

    for sim in eachindex(seir_results)
        run_seed = seed + (sim - 1)

        local seir_res = _create_run_seir_results(
            dynamics_parameter_specification,
            init_states_sv,
            beta_vec,
            time_parameters,
            run_seed
        )
        seir_results[sim] = seir_res
    end

    return StructArrays.StructVector(seir_results)
end

"""
    _create_run_seir_results(dynamics_parameter_specification, init_states,
                             beta_vec, time_parameters, run_seed = 1234)

Internal helper function to create a single SEIR simulation run.

This function constructs a `DynamicsParameters` object from the specification
and runs the SEIR model with the provided initial conditions, transmission
pattern, and time parameters.

# Arguments
- `dynamics_parameter_specification::DynamicsParameterSpecification`:
  Specification containing disease dynamics parameters including transmission
  rates, recovery rates, and vaccination coverage ranges.
- `init_states::StaticArrays.SVector{5, Int64}`: Initial state vector with
  five compartments [S, E, I, R, V] representing susceptible, exposed,
  infectious, recovered, and vaccinated populations.
- `beta_vec::Vector{Float64}`: Pre-calculated seasonal transmission rate
  pattern over the simulation time period.
- `time_parameters::SimTimeParameters`: Time-related parameters including
  simulation length, time step, and burnin period.
- `run_seed::Int64 = 1234`: Random seed for this specific simulation run.

# Returns
- `SEIRRun`: A structure containing the simulation results including state
  trajectories, incidence time series, and effective reproduction number
  (Reff) time series.

# Notes
This is an internal function (indicated by the leading underscore) used by
`generate_single_ensemble` to avoid code duplication when creating both
emergent and null scenario simulations.
"""
function _create_run_seir_results(
        dynamics_parameter_specification::DynamicsParameterSpecification,
        init_states::StaticArrays.SVector{5, Int64},
        beta_vec::Vector{Float64},
        time_parameters::SimTimeParameters,
        run_seed = 1234
    )
    dynp = DynamicsParameters(
        dynamics_parameter_specification;
        seed = run_seed
    )

    return seir_mod(
        init_states,
        dynp,
        beta_vec,
        time_parameters;
        seed = run_seed,
    )
end
