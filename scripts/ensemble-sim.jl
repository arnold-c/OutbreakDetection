#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

#%%
N_vec = convert.(Int64, [5e5])
nsims_vec = [1_000]
init_states_prop_dict = [
    Dict(
        :s_prop => 0.1,
        :e_prop => 0.01,
        :i_prop => 0.01,
        :r_prop => 0.88,
    )
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

ensemble_base_dynamics_p = DynamicsParameters(;
    sigma = sigma,
    gamma = gamma,
    R_0 = R_0
)

ensemble_time_p = SimTimeParameters(;
    tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0
)

base_param_dict = @dict(
    N = N_vec,
    init_states_prop = init_states_prop_dict,
    base_dynamics_p = ensemble_base_dynamics_p,
    time_p = ensemble_time_p,
    nsims = nsims_vec,
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
