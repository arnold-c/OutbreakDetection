"""
This is a simulation of an SIR model with a time-dependent infectivity rate.
It uses the JumpProcesses.jl package and the Coevolve method.
It uses VariableRateJumps to simulate the infectivity rate of an individual being a function of time of the year i.e. the infectivity rate peaks once a year.
All jumps are manually defined using JumpProcesses, not using the ModelingToolkit or Catalyst interfaces.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics, DifferentialEquations, ModelingToolkit
using BenchmarkTools, JLD2, Random, ProgressMeter

CairoMakie.activate!()
set_aog_theme!()

#%%
# Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/transmission-functions.jl"))
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
N = 4e5
s = 0.9
i = 0.1
r = 1.0 - (s + i)
u₀ = convert.(Int64, [N * s, N * i, N * r, N])
δt = 0.5
tlower = 0.0
tmax = 365.0 * 100
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

dur_inf = 8
R₀ = 10.0
γ = 1 / dur_inf
μ = 1 / (62 * 365)
β = calculate_beta(R₀, γ, μ, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. A_scale scales the amplitude of cosine function to 1/A_scale. The maximum infectivity is β
A_scale = 2.0
ε = (1.06 * μ * (R₀ - 1)) / sqrt(sum(u₀)) # Commuter imports - see p210 Keeling & Rohani
p = (β, γ, μ, A_scale, ε)

#%%
birth_rate(u, p, t) = p[3] * u[4]  # μ*N
function birth_affect!(integrator)
    integrator.u[1] += 1  # S -> S + 1
    integrator.u[4] += 1  # N -> N + 1
    return nothing
end
birth_jump = ConstantRateJump(birth_rate, birth_affect!)

recov_rate(u, p, t) = p[2] * u[2]         # γ*I
function recov_affect!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
    return nothing
end
recov_jump = ConstantRateJump(recov_rate, recov_affect!)

S_death_rate(u, p, t) = p[3] * u[1]  # μ*S
function S_death_affect!(integrator)
    integrator.u[1] -= 1  # S -> S - 1
    integrator.u[4] -= 1  # N -> N - 1
    return nothing
end
S_death_jump = ConstantRateJump(S_death_rate, S_death_affect!)

I_death_rate(u, p, t) = p[3] * u[2]  # μ*I
function I_death_affect!(integrator)
    integrator.u[2] -= 1  # I -> I - 1
    integrator.u[4] -= 1  # N -> N - 1
    return nothing
end
I_death_jump = ConstantRateJump(I_death_rate, I_death_affect!)

R_death_rate(u, p, t) = p[3] * u[3]  # μ*R
function R_death_affect!(integrator)
    integrator.u[3] -= 1  # R -> R - 1
    integrator.u[4] -= 1  # N -> N - 1
    return nothing
end

R_death_jump = ConstantRateJump(R_death_rate, R_death_affect!)

import_rate(u, p, t) = p[5] * u[4] / R₀   # ε*N/R₀
function import_affect!(integrator)
    integrator.u[2] += 1    # I -> I + 1
    if integrator.u[1] - 1 < 0.0  # S -> S - 1
        integrator.u[4] += 1    # N -> N + 1
    else
        integrator.u[1] - 1
    end
    return nothing
end
import_jump = ConstantRateJump(import_rate, import_affect!)

# Place infection at the end as it is a VariableRateJump, which is ordered after ConstantRateJumps in the dependency graph.
# Amplitude = 0.5 * (cos(2pi * t / 365) + 1)) * 1/scale + (scale - 1)/scale
function infec_rate(u, p, t)
    # β*S*I * 0.5*cos(2πt/365) + 0.5
    Amplitude = 0.5 * (cos(2π * t / 365) + 1) * 1 / p[4] + (p[4] - 1) / p[4]
    return p[1] * u[1] * u[2] * Amplitude
end

# Lower bound on infection rate
lrate(u, p, t) = 0.0
# Upper bound on infection rate
urate(u, p, t) = p[1] * u[1] * u[2]  # β*S*I
# Interval that the infection rate must be within the bounds for
rateinterval(u, p, t) = 1 / (2 * urate(u, p, t))

function infec_affect!(integrator)
    integrator.u[1] -= 1     # S -> S - 1
    integrator.u[2] += 1     # I -> I + 1
    return nothing
end
infec_jump = VariableRateJump(
    infec_rate, infec_affect!; lrate, urate, rateinterval
)

jumps = JumpSet(
    birth_jump, recov_jump, S_death_jump, I_death_jump, R_death_jump,
    import_jump,
    infec_jump,
)

# ConstantRateJumps are ordered before VariableRateJumps in the dependency graph, otherwise within the same ordering presented to the JumpProblem, i.e., birth, recovery ...
# Each row of the dependency graph is a list of rates that must be recalculated when the row's jump is executed, as it affects the underlying states each of the rates depends on.
dep_graph = [
    [1, 3, 6, 7],       # Birth, S death, import, infection
    [2, 4, 5, 7],       # Recovery, I death, R death, infection
    [3, 1, 6, 7],       # S death, birth, import, infection
    [4, 1, 2, 6, 7],    # I death, birth, recovery, import, infection
    [5, 1, 6],          # R death, birth, import
    [6, 1, 2, 3, 4, 7], # Import, birth, recovery, S death, I death, infection
    [7, 2, 3, 4],    # Infection, recovery, S death, I death
]

#%%
season_infec_prob = DiscreteProblem(u₀, tspan, p)
# Bounded VariableJumpRate problems require the Coevolve() algorithm
season_infec_jump_prob = JumpProblem(
    season_infec_prob, Coevolve(), jumps; dep_graph
)
season_infec_sol = solve(season_infec_jump_prob, SSAStepper())

season_infec_sol_df = create_sir_df(season_infec_sol, [:S, :I, :R, :N])

#%%
sircolors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

draw_sir_plot(season_infec_sol_df; colors = sircolors, annual = true)

#%%
# Visualize seasonal changes in the infection rate
t = 0:0.01:720
y = Vector{Float64}(undef, length(t))
y_scale = similar(y)
@. y = 0.5 * (cos(2pi * t / 365) + 1)
@. y_scale =
    (0.5 * (cos(2pi * t / 365) + 1)) * 1 / A_scale + (A_scale - 1) / A_scale
fig = lines(t, y; color = "black", linewidth = 2)
lines!(t, y_scale; color = "green", linewidth = 2)
vlines!(0:365:720; color = (:red, 0.5), linestyle = :dash, linewidth = 2)
hlines!([0.0, 0.5, 1.0])
fig

#%%
function run_ensemble_jump_prob(param_dict)
    @unpack N, u₀_prop, nsims, δt, tmax = param_dict
    @unpack s, i, r = u₀_prop

    u₀ = convert.(Int64, [s * N, i * N, r * N, N])
    u0_dict = Dict()
    for (k, v) in zip([:S, :I, :R, :N], u₀)
        u0_dict[k] = v
    end

    tspan = (0.0, tmax)

    β = calculate_beta(R₀, γ, μ, 1, N)
    ε = (1.06 * μ * (R₀ - 1)) / sqrt(sum(u₀)) # Commuter imports - see p210 Keeling & Rohani
    p = (β, γ, μ, A_scale, ε)

    remade_ensemble_prob = remake(seasonal_infec_ensemble_prob; u0 = u₀, p = p)

    ensemble_sol = solve(
        remade_ensemble_prob,
        SSAStepper(),
        EnsembleThreads();
        trajectories = nsims,
        saveat = δt,
    )

    ensemble_array = create_sir_all_sims_array(ensemble_sol, nsims)

    return @strdict ensemble_array u0_dict
end

function run_ensemble_summary(param_dict)
    @unpack N, u₀_prop, nsims, quantiles, δt, tmax = param_dict
    @unpack s, i, r = u₀_prop

    sim_name = savename(
        "jump_sol",
        param_dict,
        "jld2";
        allowedtypes = (Symbol, Dict, String, Real),
        accesses = [:N, :u₀_prop, :nsims, :tmax, :δt],
        expand = ["u₀_prop"],
        sort = false,
    )
    sim_path = joinpath(
        datadir(
            "seasonal-infectivity-import",
            "jump",
            "N_$N",
            "r_$r",
            "nsims_$nsims",
            "tmax_$tmax",
            "deltat_$δt",
        ),
        sim_name,
    )

    sol_data = load(sim_path)
    @unpack ensemble_array, u0_dict = sol_data
    S = u0_dict[:S]
    I = u0_dict[:I]
    R = u0_dict[:R]

    qlow = round(0.5 - quantiles / 200; digits = 3)
    qhigh = round(0.5 + quantiles / 200; digits = 3)

    qs = [qlow, 0.5, qhigh]

    ensemble_summary = create_sir_all_sim_quantiles(
        ensemble_array; quantiles = qs
    )

    caption = "nsims = $nsims, N = $N, S = $S, I = $I, R = $R, δt = $δt, quantile int = $quantiles"

    return @strdict ensemble_summary caption u0_dict
end

#%%
# mutable struct DrWatsonEnsembleSummary
#     u₀_params::Dict
#     time_params::Dict
#     quantiles::Int
#     ensemble_summary::Array
# end

# test = DrWatsonEnsembleSummary(
#     Dict(:N => 1e3, :s => 0.9, :i => 0.1, :r => 0.0),
#     Dict(:δt => 1.0, :tmax => 720.0),
#     95,
#     fill(nothing, (1, 1))
# )
# DrWatson.allaccess(c::DrWatsonEnsembleSummary) = (:u₀_params, :time_params, :quantiles)
# DrWatson.default_prefix(c::DrWatsonEnsembleSummary) = "jump_summ"
# DrWatson.default_allowed(c::DrWatsonEnsembleSummary) = (Dict, Dict, Int, Array)
# DrWatson.default_expand(c::DrWatsonEnsembleSummary) = ["u₀_params", "time_params"]

# savename(test)

#%%
seasonal_infec_jump_prob = JumpProblem(
    season_infec_prob, Coevolve(), jumps; dep_graph = dep_graph,
    save_positions = (false, false),
)

seasonal_infec_ensemble_prob = EnsembleProblem(seasonal_infec_jump_prob)

#%%
N_vec = convert.(Int64, [1e3, 1e4, 5e5])
nsims_vec = [10, 100, 1000]
u₀_prop_map = [
    Dict(:s => 0.9, :i => 0.1, :r => 0.0),
    Dict(:s => 0.1, :i => 0.1, :r => 0.8)
]
δt_vec = [0.5, 1.0]
tmax_vec = [365.0 * 100]

Random.seed!(1234)
base_param_dict = @dict(
    N = N_vec,
    u₀_prop = u₀_prop_map,
    nsims = nsims_vec,
    δt = δt_vec,
    tmax = tmax_vec,
)

sol_param_dict = dict_list(
    base_param_dict
)

#%%
map(
    p -> @produce_or_load(
        run_ensemble_jump_prob,
        p,
        datadir(
            "seasonal-infectivity-import",
            "jump",
            "N_$(p[:N])",
            "r_$(p[:u₀_prop][:r])",
            "nsims_$(p[:nsims])",
            "tmax_$(p[:tmax])",
            "deltat_$(p[:δt])",
        ),
        prefix = "jump_sol";
        filename = savename(
            p;
            allowedtypes = (Symbol, Dict, String, Real),
            accesses = [:N, :u₀_prop, :nsims, :tmax, :δt],
            expand = ["u₀_prop"],
            sort = false,
        ),
        loadfile = false
    ),
    sol_param_dict,
)

#%%
quantile_ints = [95, 90, 80, 50]

summ_param_dict = @chain base_param_dict begin
    deepcopy(_)
    push!(_, :quantiles => quantile_ints)
    dict_list(_)
end

#%%
map(
    p -> @produce_or_load(
        run_ensemble_summary,
        p,
        datadir(
            "seasonal-infectivity-import",
            "jump",
            "N_$(p[:N])",
            "r_$(p[:u₀_prop][:r])",
            "nsims_$(p[:nsims])",
            "tmax_$(p[:tmax])",
            "deltat_$(p[:δt])",
        ),
        prefix = "jump_quants";
        filename = savename(
            p;
            allowedtypes = (Symbol, Dict, String, Real),
            accesses = [:N, :u₀_prop, :nsims, :tmax, :δt, :quantiles],
            expand = ["u₀_prop"],
            sort = false,
        ),
        loadfile = false
    ),
    summ_param_dict,
)

#%%
δt = 1.0
sim_files = []
quantile_files = []
for (root, dirs, files) in walkdir(
    datadir(
        "seasonal-infectivity-import", "jump", "N_10000", "r_0.0"
    ),
)
    for (i, file) in enumerate(files)
        if occursin("jump_quants", file)
            push!(quantile_files, joinpath(root, file))
        end
        if occursin("jump_sol", file)
            push!(sim_files, joinpath(root, file))
        end
    end
end

sim_data = nothing
for (i, file) in enumerate(sim_files)
    if occursin(r"jump_sol.*.nsims=1000_.*.δt=1.0", file)
        sim_data = load(sim_files[i])
    end
end

@unpack ensemble_array = sim_data

summ_data = nothing
for (i, file) in enumerate(quantile_files)
    if occursin(r"jump_quants.*.nsims=1000_.*.δt=1.0.*.quantiles=95", file)
        summ_data = load(quantile_files[i])
    end
end

@unpack ensemble_summary, caption = summ_data

create_sir_quantiles_plot(
    ensemble_summary; annual = true, caption = caption
)

#%%
f = Figure()
ax = Axis(f[1, 1])

for sim in 1:size(ensemble_array, 3)
    lines!(
        ax,
        1:size(ensemble_array, 2),
        ensemble_array[2, :, sim];
        color = "red",
        alpha = 0.01,
    )
end

hlines!(ax, 5, color = "black", linestyle = :dash, label = "5")

f

#%%
above_5 = zeros(size(ensemble_array, 3), size(ensemble_array, 2))

for j in 1:size(ensemble_array, 2)
    for k in 1:size(ensemble_array, 3)
        above_5[k, j] = ensemble_array[2, j, k] > 5
    end
end

plot(above_5[1, :])

test = zeros(length(above_5[1, :]), 2)
test[:, 1] = above_5[1, :]

test[:, 2] = reduce(vcat, map(
    len -> repeat([len], len),
    rle(test[:, 1])[2]
))

scatter(1:size(test, 1), test[:, 2]; color = test[:, 1])