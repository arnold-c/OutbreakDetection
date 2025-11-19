using UnPack: UnPack
using StaticArrays: StaticArrays

export setup_optimization

"""
    setup_optimization(ensemble_param_dict)

Set up ensemble simulation for optimization.

Returns ensemble_inc_arr and ensemble_thresholds_vec.
"""
function setup_optimization(ensemble_param_dict)
    UnPack.@unpack ensemble_spec,
        seed,
        outbreak_spec = ensemble_param_dict

    UnPack.@unpack state_parameters,
        dynamics_parameters, time_parameters,
        nsims =
        ensemble_spec

    UnPack.@unpack tstep, tlength, trange = time_parameters

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states), 2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_inc_vecs = Array{typeof(StaticArrays.SVector(0)), 2}(
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
            dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_spec
    )

    return ensemble_inc_arr, ensemble_thresholds_vec
end
