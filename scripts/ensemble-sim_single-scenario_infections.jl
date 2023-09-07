#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack

include("ensemble-sim.jl")

#%%
ensemble_spec = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    500_000,
    Dict(
        :s_prop => 0.1,
        :e_prop => 0.01,
        :i_prop => 0.01,
        :r_prop => 0.88,
    ),
    1_000,
    10,
    0.2,
    ensemble_time_p,
)

ensemble_sol = get_ensemble_file(
    "sol", ensemble_spec
)

ensemble_quants = get_ensemble_file(
    "95", ensemble_spec
)

@unpack ensemble_seir_arr,
ensemble_jump_arr, ensemble_change_arr, ensemble_dynamics_p,
ensemble_param_dict = ensemble_sol
@unpack ensemble_seir_summary, caption = ensemble_quants
@unpack time_p = ensemble_param_dict

#%%
outbreak_spec = OutbreakSpecification(5, 30, 500)

inc_infec_arr = create_inc_infec_arr(ensemble_jump_arr, outbreak_spec)
