#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter

includet(srcdir("Julia/DrWatson-helpers.jl"))
includet(funsdir("ensemble-functions.jl"))

#%%
N_vec = convert.(Int64, [5e5])
nsims_vec = [1000]
u₀_prop_map = [
    Dict(:s => 0.1, :e => 0.01, :i => 0.01, :r => 0.88)
]
dt_vec = [1.0]
tmax_vec = [365.0 * 100]
beta_force_vec = collect(0.0:0.1:0.4)
μ_min = 5
μ_max = 20
μ_step = 5.0
n_μs = length(μ_min:μ_step:μ_max)
μ_vec = zeros(Float64, n_μs)
μ_vec = convert.(Int64, collect(μ_min:μ_step:μ_max))
seed = 1234

latent_per_days = 8
dur_inf_days = 5
R₀ = 10.0
sigma = 1 / latent_per_days
gamma = 1 / dur_inf_days
transmission_p = EnsembleTransmissionParameters(R₀, sigma, gamma)
tmin = 0.0
tmax = 365.0 * 100
tstep = 1.0
ensemble_time_p = EnsembleTimeParameters(tmin, tmax, tstep)

base_param_dict = @dict(
    N = N_vec,
    u₀_prop = u₀_prop_map,
    transmission_p = transmission_p,
    time_p = time_p,
    nsims = nsims_vec,
    dt = dt_vec,
    tmax = tmax_vec,
    beta_force = beta_force_vec,
    births_per_k = μ_vec,
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
    500000,
    0.88,
    1000,
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
