export EnsembleSpecification

"""
    EnsembleSpecification

Specification for an ensemble of simulations.

# Fields
- `modeltypes::Tuple`: Types of models to run
- `state_parameters::StateParameters`: Initial state parameters
- `dynamics_parameters::DynamicsParameters`: Disease dynamics parameters
- `time_parameters::SimTimeParameters`: Time parameters
- `nsims::Integer`: Number of simulations
- `dirpath::AbstractString`: Directory path for output

# Constructor
    EnsembleSpecification(modeltypes, state_parameters, dynamics_parameters,
                         time_parameters, nsims)

Creates an `EnsembleSpecification` with automatically generated directory path.
"""
struct EnsembleSpecification{
        T1 <: Tuple,
        T2 <: StateParameters,
        T3 <: DynamicsParameters,
        T4 <: SimTimeParameters,
        T5 <: Integer,
        T6 <: AbstractString,
    }
    modeltypes::T1
    state_parameters::T2
    dynamics_parameters::T3
    time_parameters::T4
    nsims::T5
    dirpath::T6
end

function EnsembleSpecification(
        modeltypes::Tuple,
        state_parameters::StateParameters,
        dynamics_parameters::DynamicsParameters,
        time_parameters::SimTimeParameters,
        nsims::Int64,
    )
    dirpath = outdir(
        "ensemble",
        modeltypes...,
        "N_$(state_parameters.init_states.N)",
        "r_$(state_parameters.init_state_props.r_prop)",
        "nsims_$(nsims)",
        "R0_$(dynamics_parameters.R_0)",
        "latent_period_$(round(1 / dynamics_parameters.sigma; digits = 2))",
        "infectious_period_$(round(1 / dynamics_parameters.gamma; digits = 2))",
        "vaccination_coverage_$(dynamics_parameters.vaccination_coverage)",
        "births_per_k_$(dynamics_parameters.annual_births_per_k)",
        "beta_force_$(dynamics_parameters.beta_force)",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)",
    )

    return EnsembleSpecification(
        modeltypes,
        state_parameters,
        dynamics_parameters,
        time_parameters,
        nsims,
        dirpath,
    )
end
