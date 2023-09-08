#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using Chain

include("../src/OutbreakDetection.jl")
using .OutbreakDetection

#%%
model_types_vec = [("seasonal-infectivity-import", "tau-leaping")]

#%%
N_vec = [500_000]
nsims_vec = [1_000]
init_states_prop_dict = [
    Dict(:s_prop => 0.1, :e_prop => 0.01, :i_prop => 0.01, :r_prop => 0.88)
]

ensemble_state_p_vec = create_combinations_vec(
    StateParameters, (N_vec, init_states_prop_dict)
)

#%%
tmin_vec = [0.0]
tstep_vec = [1.0]
tmax_vec = [365.0 * 100]

time_p_vec = vec(
    map(
        Iterators.product(tmin_vec, tstep_vec, tmax_vec)
    ) do (tmin, tstep, tmax)
        SimTimeParameters(; tmin = tmin, tmax = tmax, tstep = tstep)
    end,
)

#%%
beta_force_vec = collect(0.0:0.1:0.4)
annual_births_per_k_min = 5
annual_births_per_k_max = 20
annual_births_per_k_step = 5
annual_births_per_k_vec = collect(
    annual_births_per_k_min:annual_births_per_k_step:annual_births_per_k_max
)
seed = 1234

#%%
latent_per_days_vec = [8]
dur_inf_days_vec = [5]
R_0_vec = [10.0]
sigma_vec = 1 ./ latent_per_days_vec
gamma_vec = 1 ./ dur_inf_days_vec

#%%
ensemble_spec_vec = create_ensemble_spec_combinations(
    beta_force_vec,
    sigma_vec,
    gamma_vec,
    annual_births_per_k_vec,
    R_0_vec,
    N_vec,
    init_states_prop_dict,
    model_types_vec,
    time_p_vec,
    nsims_vec,
)

#%%
base_param_dict = @dict(
    ensemble_spec = ensemble_spec_vec,
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

quantile_param_dict = @chain base_param_dict begin
    deepcopy(_)
    push!(_, :quantiles => quantile_ints)
    dict_list(_)
end;

#%%
prog = Progress(length(quantile_param_dict))
summarize_ensemble_jump_prob(quantile_param_dict; prog = prog)
