#%%
using DrWatson
@quickactivate "OutbreakDetection"

using Random
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
β_force_vec = collect(0.0:0.1:0.4)
μ_min = 5
μ_max = 20
μ_step = 5.0
n_μs = length(μ_min:μ_step:μ_max)
μ_vec = zeros(Float64, n_μs)
μ_vec = convert.(Int64, collect(μ_min:μ_step:μ_max))

Random.seed!(1234)
base_param_dict = @dict(
    N = N_vec,
    u₀_prop = u₀_prop_map,
    nsims = nsims_vec,
    dt = dt_vec,
    tmax = tmax_vec,
    β_force = β_force_vec,
    births_per_k = μ_vec,
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

