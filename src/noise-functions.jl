# module NoiseFunctions
#
# export create_noise_arr, create_noise_arr!

# include("ensemble-functions.jl")
# using .EnsembleFunctions
using UnPack
using Match

function create_noise_arr(
    noise_specification::DynamicalNoiseSpecification;
    ensemble_specification::EnsembleSpecification,
    seed = 1234,
    kwargs...,
)
    seed *= 10

    @unpack state_parameters, dynamics_parameters, time_parameters, nsims =
        ensemble_specification
    @unpack tlength = time_parameters

    noise_beta_force = @match noise_specification.correlation begin
        "none" => 0.0
        _ => dynamics_parameters.beta_force
    end

    noise_seasonality = @match noise_specification.correlation begin
        "out-of-phase" => @match dynamics_parameters.seasonality begin
            cos => sin
            sin => cos
        end
        _ => dynamics_parameters.seasonality
    end

    noise_dynamics_parameters = DynamicsParameters(
        dynamics_parameters.beta_mean,
        noise_beta_force,
        noise_seasonality,
        dynamics_parameters.sigma,
        dynamics_parameters.gamma,
        dynamics_parameters.mu,
        dynamics_parameters.annual_births_per_k,
        calculate_import_rate(
            dynamics_parameters.mu,
            noise_specification.R_0,
            state_parameters.init_states.N,
        ),
        noise_specification.R_0,
        dynamics_parameters.vaccination_coverage,
    )

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states),2}(
        undef,
        tlength,
        nsims
    )

    ensemble_inc_vecs = Array{typeof(SVector(0)),2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_beta_arr = zeros(Float64, tlength)

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            noise_dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    return convert_svec_to_matrix(ensemble_inc_vecs)
end

function create_noise_arr(
    noise_specification::PoissonNoiseSpecification;
    incarr,
    seed = 1234,
    kwargs...,
)
    noise_arr = zeros(Int64, size(incarr, 1), size(incarr, 3))

    create_poisson_noise_arr!(
        noise_arr, incarr, noise_specification.noise_mean_scaling; seed = seed
    )

    return noise_arr
end

function create_poisson_noise_arr!(
    noise_arr, incarr, noise_mean_scaling; seed = 1234
)
    Random.seed!(seed)
    @inbounds for sim in axes(incarr, 3)
        @views noise_arr[:, sim] .= rand(
            Poisson(noise_mean_scaling * mean(incarr[:, 1, sim])),
            size(incarr, 1),
        )
    end
end

# end
