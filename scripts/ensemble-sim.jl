#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("ensemble-functions.jl"))

#%%
N_vec = convert.(Int64, [5e5])
nsims_vec = [1_000]
init_states_prop_map = [
    Dict(
        :s_init_prop => 0.1,
        :e_init_prop => 0.01,
        :i_init_prop => 0.01,
        :r_init_prop => 0.88,
    ),
]
tstep_vec = [1.0]
tmax_vec = [365.0 * 100]
beta_force_vec = collect(0.0:0.1:0.4)
births_per_k_min = 5
births_per_k_max = 20
births_per_k_step = 5
births_per_k_vec = collect(births_per_k_min:births_per_k_step:births_per_k_max)
seed = 1234

latent_per_days = 8
dur_inf_days = 5
R_0 = 10.0
sigma = 1 / latent_per_days
gamma = 1 / dur_inf_days
transmission_p = EnsembleTransmissionParameters(R_0, sigma, gamma)
tmin = 0.0
tmax = 365.0 * 100
tstep = 1.0
ensemble_time_p = EnsembleTimeParameters(tmin, tmax, tstep)

base_param_dict = @dict(
    N = N_vec,
    init_states_prop = init_states_prop_map,
    transmission_p = transmission_p,
    time_p = time_p,
    nsims = nsims_vec,
    tstep = tstep_vec,
    tmax = tmax_vec,
    beta_force = beta_force_vec,
    births_per_k = births_per_k_vec,
    seed = seed,
)

sol_param_dict = dict_list(
    base_param_dict
)

#%%
prog = Progress(length(sol_param_dict))
run_ensemble_jump_prob(sol_param_dict; prog = prog)

#%%
quantile_ints = [95, 80]

summ_param_dict = @chain base_param_dict begin
    deepcopy(_)
    push!(_, :quantiles => quantile_ints)
    dict_list(_)
end;

#%%
prog = Progress(length(summ_param_dict))
summarize_ensemble_jump_prob(summ_param_dict; prog = prog)

#%%
ensemble_spec = EnsembleSpecification(
    ("seasonal-infectivity-import", "tau-leaping"),
    500_000,
    0.88,
    1_000,
    20,
    0.2,
    ensemble_time_p,
)

ensemble_sol = get_ensemble_file(
    "sol", ensemble_spec
)

ensemble_quants = get_ensemble_file(
    "95", ensemble_spec
)

@unpack ensemble_seir_arr, ensemble_jump_arr, ensemble_change_arr = ensemble_sol
@unpack ensemble_seir_summary, caption, param_dict = ensemble_quants
