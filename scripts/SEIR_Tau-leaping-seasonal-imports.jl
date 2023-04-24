"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics, ColorSchemes, Colors
using DifferentialEquations, ModelingToolkit
using BenchmarkTools, JLD2, Random, ProgressMeter, StatsBase, Distributions
using IterTools, FLoops

CairoMakie.activate!()
set_aog_theme!()

#%%
# Revise will keep track of file change_arr and reload functions as necessary
includet(srcdir("Julia/transmission-functions.jl"))
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
N = 5e5
s = 0.1
e = 0.01
i = 0.01
r = 1.0 - (s + e + i)
u₀ = convert.(Int64, [N * s, N * e, N * i, N * r, N])
τ = 1.0
tlower = 0.0
tmax = 365.0 * 100
trange = tlower:τ:tmax
tlength = length(trange)

latent_per = 8
dur_inf = 5
R₀ = 10.0
σ = 1 / latent_per
γ = 1 / dur_inf
μ = 1 / (62.5 * 365)
# β₀ is the average transmission rate
β₀ = calculate_beta(R₀, γ, μ, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. β₁ scales the amplitude of cosine function
β₁ = 0.2
ε = (1.06 * μ * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
p = (β₀, β₁, σ, γ, μ, ε, R₀)

Random.seed!(1234)

#%%
function calculate_beta_amp(β_mean, β_force, t)
    return β_mean * (1 + β_force * cos(2pi * t / 365))
end

function seir_mod(u, p, trange; retβarr = false, type = "stoch")
    tlength = length(trange)
    dt = step(trange)

    state_arr = zeros(Float64, size(u, 1), tlength)

    change_arr = zeros(Float64, size(u, 1), tlength)

    jump_arr = zeros(Float64, 9, tlength)

    if retβarr == true
        beta_arr = zeros(Float64, tlength)
        seir_mod!(
            state_arr, change_arr, jump_arr, beta_arr, u, p, trange; dt = dt,
            type = type,
        )
        return state_arr, change_arr, jump_arr, beta_arr
    else
        seir_mod!(
            state_arr, change_arr, jump_arr, u, p, trange; dt = dt, type = type
        )
        return state_arr, change_arr, jump_arr
    end
end

function seir_mod_loop!(
    state_arr, change_arr, jump_arr, j, p, t, dt; type = type
)
    # Unpack the state variables for easier use
    S = state_arr[1, j - 1]
    E = state_arr[2, j - 1]
    I = state_arr[3, j - 1]
    R = state_arr[4, j - 1]
    N = state_arr[5, j - 1]

    # Unpack the parameters for easier use
    β_mean, β_force, σ, γ, μ, ε, R₀ = p
    β_t = calculate_beta_amp(β_mean, β_force, t)

    # Calculate the rates of each event
    infec_rate = β_t * S * I        # 1
    latent_rate = σ * E             # 2
    recov_rate = γ * I              # 3
    birth_rate = μ * N              # 4
    S_death_rate = μ * S            # 5
    E_death_rate = μ * E            # 6
    I_death_rate = μ * I            # 7
    R_death_rate = μ * R            # 8
    import_rate = ε * N / R₀        # 9

    rates = [
        infec_rate, latent_rate, recov_rate, birth_rate, S_death_rate,
        E_death_rate, I_death_rate, R_death_rate, import_rate]

    # Calculate the number of jumps for each event
    if type == "stoch"
        jump_arr[:, j] = map(
            r -> rand(Poisson(r * dt)),
            rates,
        )
    elseif type == "det"
        jump_arr[:, j] = map(
            r -> r * dt,
            rates,
        )
    else
        return ("Type must be stoch or det")
    end

    # Calculate the change in each state
    change_arr[1, j] =
        jump_arr[4, j] - jump_arr[1, j] - jump_arr[5, j] - jump_arr[9, j]
    change_arr[2, j] =
        jump_arr[1, j] - jump_arr[2, j] - jump_arr[6, j] + jump_arr[9, j]
    change_arr[3, j] = jump_arr[2, j] - jump_arr[3, j] - jump_arr[7, j]
    change_arr[4, j] = jump_arr[3, j] - jump_arr[8, j]
    change_arr[5, j] = sum(change_arr[1:4, j])

    # Check that the change in each state does not result in a negative state, 
    # and if it is, set the change to the negative of the current state
    for state in 1:size(state_arr, 1)
        if state_arr[state, j - 1] + change_arr[state, j] < 0
            change_arr[state, j] = -state_arr[state, j - 1]
        end
    end

    @. state_arr[:, j] = state_arr[:, j - 1] + change_arr[:, j]

    return nothing
end

function seir_mod!(
    state_arr, change_arr, jump_arr, u, p, trange; dt, type = "stoch"
)
    S0, I0, R0, N0 = u

    for (j, t) in pairs(trange)
        if j == 1
            state_arr[:, j] = u
            continue
        end

        seir_mod_loop!(
            state_arr, change_arr, jump_arr, j, p, t, dt; type = type
        )
    end

    return nothing
end

function seir_mod!(
    state_arr, change_arr, jump_arr, beta_arr, u, p, trange; dt, type = "stoch"
)
    S0, I0, R0, N0 = u
    β_mean, β_force, σ, γ, μ, ε, R₀ = p

    for (j, t) in pairs(trange)
        β_t = calculate_beta_amp(β_mean, β_force, t)
        beta_arr[j] = β_t

        if j == 1
            state_arr[:, j] = u

            continue
        end

        seir_mod_loop!(
            state_arr, change_arr, jump_arr, j, p, t, dt; type = type
        )
    end

    return nothing
end

#%%
seir_array, change_array, jump_array = seir_mod(
    u₀, p, trange; retβarr = true, type = "stoch"
);

seir_df = create_sir_df(seir_array, trange, [:S, :E, :I, :R, :N])

seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
state_labels = ["S", "E", "I", "R", "N"]
draw_sir_plot(
    seir_df;
    annual = true,
    colors = seircolors,
    labels = state_labels,
    # ylims = (0, 1000),
    # xlims = (90, 95),
)

#%%
@chain DataFrame(Tables.table(jump_array')) begin
    hcat(trange, _)
    rename!([
        "time", "Infect", "Latent", "Recov", "Birth", "S_death", "E_death",
        "I_death", "R_death",
        "Import",
    ])
    stack(_, Not("time"); variable_name = :Jump, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)",
        :Number;
        color = :Jump,
        layout = :Jump,
    ) *
    visual(Lines; linewidth = 1)
    draw(;
        facet = (; linkyaxes = :none), axis = (; limits = (nothing, nothing))
    )
end

#%%
change_labels = ["dS", "dE", "dI", "dR", "dN"]
@chain DataFrame(Tables.table(change_array')) begin
    hcat(trange, _)
    rename!([
        "time", change_labels...
    ])
    stack(_, Not("time"); variable_name = :Change, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)", :Number;
        color = :Change => sorter(change_labels...), layout = :Change,
    ) *
    visual(Lines; linewidth = 1)
    draw(;
        facet = (; linkyaxes = :none),
        palettes = (; color = seircolors),
        axis = (; limits = ((0, 100), nothing)),
    )
end

#%%
@chain DataFrame(Tables.table(seir_array')) begin
    hcat(trange, _)
    rename!(["time", state_labels...])
    data(_) *
    mapping(:I, :S; color = :time) *
    visual(Lines)
    draw
end

#%%
@chain DataFrame(Tables.table(β_arr)) begin
    hcat(trange, _)
    rename!([:time, :β_t])
    stack(_, Not("time"); variable_name = :β, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)", :Number;
    ) *
    visual(Lines; linewidth = 1)
    draw(; facet = (; linkyaxes = :none), axis = (; limits = ((0, 3), nothing)))
end

#%%
################################################################################
########################## Bifurcation Diagrams ################################
################################################################################

μ_min = 10
μ_max = 100
μ_step = 2.0
n_μs = length(μ_min:μ_step:μ_max)
μ_vec = zeros(Float64, n_μs)
μ_vec .= collect(μ_min:μ_step:μ_max) ./ (1000 * 365)

ε_vec = (1.06 .* μ_vec .* (R₀ - 1)) ./ sqrt(N)

bifurc_μ_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_μs);
bifurc_μ_change_arr = zeros(Float64, size(u₀, 1), tlength, n_μs);
bifurc_μ_jump_arr = zeros(Float64, 9, tlength, n_μs);

prog = Progress(n_μs)
@floop for (k, μ_run) in pairs(μ_vec)
    ε = ε_vec[k]

    seir = @view bifurc_μ_seir_arr[:, :, k]
    change = @view bifurc_μ_change_arr[:, :, k]
    jump = @view bifurc_μ_jump_arr[:, :, k]

    params = (β₀, β₁, σ, γ, μ_run, ε, R₀)
    seir_mod!(seir, change, jump,
        u₀, params, trange; dt = τ, type = "det",
    )
    next!(prog)
end

years = (40 * 365):365:(tlength - 365)
bifurc_μ_annual_summary = zeros(Float64, 5, length(years), n_μs);

for μ in eachindex(μ_vec), state in eachindex(state_labels),
    (year, day) in pairs(years)

    bifurc_μ_annual_summary[state, year, μ] = maximum(
        bifurc_μ_seir_arr[state, day:(day + 364), μ]
    )
end

bifurc_μ_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_μ_seir_arr[2, (40 * 365):(40 * 365 + 364), 10]

#%%
bifurc_μ_fig = Figure()
bifurc_μ_ax = Axis(
    bifurc_μ_fig[1, 1]; xlabel = "μ (per 1000, per annum)", ylabel = "Max. I"
)

for year in eachindex(years)
    scatter!(
        bifurc_μ_ax,
        μ_min:μ_step:μ_max,
        bifurc_μ_annual_summary[2, year, :];
        markersize = 4,
        color = :black,
    )
end

bifurc_μ_fig

#%%
β₁_min = 0.0
β₁_max = 1.0
β₁_step = 0.02
n_β₁s = length(β₁_min:β₁_step:β₁_max)
β₁_vec = zeros(Float64, n_β₁s)
β₁_vec .= collect(β₁_min:β₁_step:β₁_max)

bifurc_β₁_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_β₁s);
bifurc_β₁_change_arr = zeros(Float64, size(u₀, 1), tlength, n_β₁s);
bifurc_β₁_jump_arr = zeros(Float64, 9, tlength, n_β₁s);

μ = 50 / (1000 * 365)
ε = 1.06 * μ * (R₀ - 1) / sqrt(N)

@showprogress for (k, β₁) in pairs(β₁_vec)
    seir = @view bifurc_β₁_seir_arr[:, :, k]
    change = @view bifurc_β₁_change_arr[:, :, k]
    jump = @view bifurc_β₁_jump_arr[:, :, k]

    p = (β₀, β₁, σ, γ, μ, ε, R₀)

    seir_mod!(seir, change, jump, u₀, p, trange; dt = τ, type = "det")
end

years = (40 * 365):365:(tlength - 365)
bifurc_β₁_annual_summary = zeros(Float64, 5, length(years), n_β₁s)

@floop for (k, β₁) in pairs(β₁_vec), (year, day) in pairs(years),
    state in eachindex(state_labels)
    bifurc_β₁_annual_summary[state, year, k] = maximum(
        bifurc_β₁_seir_arr[state, day:(day + 364), k]
    )
end

bifurc_β₁_seir_arr[2, (40 * 365):(40 * 365 + 364), 1] ==
bifurc_β₁_seir_arr[2, (40 * 365):(40 * 365 + 364), 11]

#%%
bifurc_β₁_fig = Figure()
bifurc_β₁_ax = Axis(
    bifurc_β₁_fig[1, 1]; xlabel = "β₁ (seasonality)", ylabel = "Max. I"
)

for year in eachindex(years)
    scatter!(
        bifurc_β₁_ax,
        β₁_min:β₁_step:β₁_max,
        bifurc_β₁_annual_summary[2, year, :];
        markersize = 4,
        color = :black,
    )
end

bifurc_β₁_fig

#%%
bifurc_μ_β₁_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_μs, n_β₁s);
bifurc_μ_β₁_change_arr = zeros(Float64, size(u₀, 1), tlength, n_μs, n_β₁s);
bifurc_μ_β₁_jump_arr = zeros(Float64, 9, tlength, n_μs, n_β₁s);

prog = Progress(length(μ_vec) * length(β₁_vec))
@floop for (β₁_pair, μ_pair) in
           IterTools.product(pairs(β₁_vec), pairs(μ_vec))
    k = μ_pair[1]
    μ = μ_pair[2]
    l = β₁_pair[1]
    β₁ = β₁_pair[2]
    ε = ε_vec[k]

    seir = @view bifurc_μ_β₁_seir_arr[:, :, k, l]
    change = @view bifurc_μ_β₁_change_arr[:, :, k, l]
    jump = @view bifurc_μ_β₁_jump_arr[:, :, k, l]

    params = (β₀, β₁, σ, γ, μ, ε, R₀)

    seir_mod!(seir, change, jump, u₀, params, trange; dt = τ, type = "det")
    next!(prog)
end

#%%
years = (40 * 365):365:(tlength - 365)

bifurc_μ_β₁_annual_summary = zeros(
    Float64, size(u₀, 1), length(years), n_μs, n_β₁s
);
prog = Progress(length(years) * 5 * length(β₁_vec))
@floop for (β₁, state, years) in
           IterTools.product(eachindex(β₁_vec), eachindex(u₀), pairs(years))
    year = years[1]
    day = years[2]

    bifurc_μ_β₁_annual_summary[state, year, :, β₁] = maximum(
        bifurc_μ_β₁_seir_arr[state, day:(day + 364), :, β₁]; dims = 1
    )
    next!(prog)
end

bifurc_μ_β₁_cycle_summary = zeros(Float64, n_β₁s, n_μs, 5);
prog = Progress(5 * length(β₁_vec) * length(μ_vec))
@floop for (β₁, μ, state) in
           IterTools.product(eachindex(β₁_vec), eachindex(μ_vec), 1:5)
    bifurc_μ_β₁_cycle_summary[β₁, μ, state] = length(
        Set(round.(bifurc_μ_β₁_annual_summary[state, :, μ, β₁]))
    )
    next!(prog)
end

#%%
# Because Julia is column-major, we need to transpose the heatmap as it reads the data one column at a time, and in our original matrix, each column represents a different β₁ value.
bifurc_μ_β₁_fig, bifurc_μ_β₁_ax, bifurc_μ_β₁_hm = heatmap(
    μ_vec .* (1000 * 365), β₁_vec, bifurc_μ_β₁_cycle_summary[:, :, 2]'
)
Colorbar(bifurc_μ_β₁_fig[:, end + 1], bifurc_μ_β₁_hm)

bifurc_μ_β₁_ax.xlabel = "μ (per 1000, per annum)"
bifurc_μ_β₁_ax.ylabel = "β₁ (seasonality)"

bifurc_μ_β₁_fig

#%%
################################################################################
########################## Ensemble Analysis ###################################
################################################################################
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
function run_ensemble_jump_prob(param_dict)
    @unpack N, u₀_prop, nsims, dt, tmax, β_force, births_per_k = param_dict
    @unpack s, e, i, r = u₀_prop

    u₀ = convert.(Int64, [s * N, e * N, i * N, r * N, N])
    u0_dict = Dict()
    for (k, v) in zip([:S, :E, :I, :R, :N], u₀)
        u0_dict[k] = v
    end

    tspan = (0.0, tmax)

    μ = births_per_k / (1000 * 365)

    β₀ = calculate_beta(R₀, γ, μ, 1, N)
    ε = (1.06 * μ * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
    p = (β₀, β_force, σ, γ, μ, ε, R₀)

    ensemble_seir_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
    ensemble_change_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
    ensemble_jump_arr = zeros(Int64, 9, tlength, nsims)

    @floop for k in 1:nsims
        @views seir = ensemble_seir_arr[:, :, k]
        @views change = ensemble_change_arr[:, :, k]
        @views jump = ensemble_jump_arr[:, :, k]

        seir_mod!(seir, change, jump, u₀, p, trange; dt = τ)
    end

    return @strdict ensemble_seir_arr ensemble_change_arr ensemble_jump_arr u0_dict param_dict
end

#%%
prog = Progress(length(sol_param_dict))
for p in sol_param_dict
    @produce_or_load(
        run_ensemble_jump_prob,
        p,
        datadir(
            "seasonal-infectivity-import",
            "tau-leaping",
            "N_$(p[:N])",
            "r_$(p[:u₀_prop][:r])",
            "nsims_$(p[:nsims])",
            "births_per_k_$(p[:births_per_k])",
            "beta_force_$(p[:β_force])",
            "tmax_$(p[:tmax])",
            "deltat_$(p[:dt])",
        );
        prefix = "SEIR_tau_sol",
        filename = savename(
            p;
            allowedtypes = (Symbol, Dict, String, Real),
            accesses = [
                :N, :u₀_prop, :nsims, :tmax, :dt, :births_per_k, :β_force
            ],
            expand = ["u₀_prop"],
            sort = false,
        ),
        loadfile = false
    )
    next!(prog)
end

#%%
function run_ensemble_summary(param_dict)
    @unpack N, u₀_prop, nsims, dt, tmax, β_force, births_per_k, quantiles =
        param_dict
    @unpack s, e, i, r = u₀_prop

    sim_name = savename(
        "SEIR_tau_sol",
        param_dict,
        "jld2";
        allowedtypes = (Symbol, Dict, String, Real),
        accesses = [:N, :u₀_prop, :nsims, :tmax, :dt, :births_per_k, :β_force],
        expand = ["u₀_prop"],
        sort = false,
    )
    sim_path = joinpath(
        datadir(
            "seasonal-infectivity-import",
            "tau-leaping",
            "N_$N",
            "r_$r",
            "nsims_$nsims",
            "births_per_k_$births_per_k",
            "beta_force_$β_force",
            "tmax_$tmax",
            "deltat_$dt",
        ),
        sim_name,
    )

    sol_data = load(sim_path)
    @unpack ensemble_seir_arr, u0_dict = sol_data
    S = u0_dict[:S]
    E = u0_dict[:E]
    I = u0_dict[:I]
    R = u0_dict[:R]

    qlow = round(0.5 - quantiles / 200; digits = 3)
    qhigh = round(0.5 + quantiles / 200; digits = 3)

    qs = [qlow, 0.5, qhigh]

    ensemble_seir_summary = create_sir_all_sim_quantiles(
        ensemble_seir_arr; quantiles = qs
    )

    caption = "nsims = $nsims, N = $N, S = $S, I = $I, R = $R, β_force = $β_force,\nbirths per k/annum = $births_per_k dt = $dt, quantile int = $quantiles"

    return @strdict ensemble_seir_summary caption u0_dict param_dict
end

#%%
quantile_ints = [95, 80]

summ_param_dict = @chain base_param_dict begin
    deepcopy(_)
    push!(_, :quantiles => quantile_ints)
    dict_list(_)
end;

#%%
prog = Progress(length(summ_param_dict))
for p in summ_param_dict
    @produce_or_load(
        run_ensemble_summary,
        p,
        datadir(
            "seasonal-infectivity-import",
            "tau-leaping",
            "N_$(p[:N])",
            "r_$(p[:u₀_prop][:r])",
            "nsims_$(p[:nsims])",
            "births_per_k_$(p[:births_per_k])",
            "beta_force_$(p[:β_force])",
            "tmax_$(p[:tmax])",
            "deltat_$(p[:dt])",
        ),
        prefix = "SEIR_tau_quants";
        filename = savename(
            p;
            allowedtypes = (Symbol, Dict, String, Real),
            accesses = [
                :N, :u₀_prop, :nsims, :tmax, :dt, :births_per_k, :β_force,
                :quantiles,
            ],
            expand = ["u₀_prop"],
            sort = false,
        ),
        loadfile = false
    )
    next!(prog)
end

#%%
sim_files = []
quantile_files = []
for (root, dirs, files) in walkdir(
    datadir(
        "seasonal-infectivity-import", "tau-leaping", "N_500000", "r_0.88",
        "nsims_1000",
    ),
)
    for (i, file) in enumerate(files)
        if occursin("SEIR_tau_quants", file)
            push!(quantile_files, joinpath(root, file))
        end
        if occursin("SEIR_tau_sol", file)
            push!(sim_files, joinpath(root, file))
        end
    end
end

sim_data = nothing
for (i, file) in enumerate(sim_files)
    if occursin(
        r"SEIR_tau_sol.*.nsims=1000_.*.births_per_k=20.*.β_force=0.2", file
    )
        sim_data = load(sim_files[i])
    end
end

@unpack ensemble_seir_arr, ensemble_jump_arr, ensemble_change_arr = sim_data

summ_data = nothing
for (i, file) in enumerate(quantile_files)
    if occursin(
        r"SEIR_tau_quant.*.nsims=1000_.*.births_per_k=20.*.β_force=0.2.*.quantiles=95.jld2",
        file,
    )
        summ_data = load(quantile_files[i])
    end
end

@unpack ensemble_seir_summary, caption, param_dict = summ_data

create_sir_quantiles_plot(
    ensemble_seir_summary; labels = state_labels, colors = seircolors,
    annual = true, caption = caption, δt = param_dict[:dt], xlims = (80, 100),
    ylims = (0, 1000),
)

#%%
################################################################################
########################## Above-Below Analysis ################################
################################################################################
inc_infec_arr = zeros(
    Int64, 4, size(ensemble_jump_arr, 2), size(ensemble_jump_arr, 3)
)

prog = Progress(size(ensemble_jump_arr, 3))
@floop for sim in 1:size(ensemble_jump_arr, 3)
    # Copy new infections to array
    inc_infec_arr[1, :, sim] = @view(ensemble_jump_arr[1, :, sim])
    # Calculate if new infection is above or below threshold
    inc_infec_arr[2, :, sim] = @view(inc_infec_arr[1, :, sim]) .> 5

    # Calculate the total number of infections above threshold in a consecutive string of days
    ## Calculate the number of consecutive days of infection above or below threshold
    above5rle = rle(@view(inc_infec_arr[2, :, sim]))

    ## Calculate upper and lower indices of consecutive days of infection
    above5accum = accumulate(+, above5rle[2])
    above5uppers = above5accum[findall(==(1), above5rle[1])]
    above5lowers = filter(
        x -> x <= maximum(above5uppers),
        above5accum[findall(==(0), above5rle[1])] .+ 1,
    )

    for (lower, upper) in zip(above5lowers, above5uppers)
        # Calculate number of infections between lower and upper indices
        period_sum = sum(@view(inc_infec_arr[1, lower:upper, sim]))
        inc_infec_arr[3, lower:upper, sim] .= period_sum

        # Determine if there is an outbreak between lower and upper indices
        if upper - lower >= 30 && period_sum >= 500
            inc_infec_arr[4, lower:upper, sim] .= 1
        end
    end

    next!(prog)
end

#%%
above5fig = Figure()
above5ax_prev = Axis(above5fig[1, 1]; ylabel = "Prevalence")
above5ax_inc = Axis(above5fig[2, 1]; ylabel = "Incidence")
above5ax_periodsum = Axis(
    above5fig[3, 1]; xlabel = "Time (years)", ylabel = "Period Sum"
)

linkxaxes!(above5ax_prev, above5ax_inc, above5ax_periodsum)

times = collect(0:param_dict[:dt]:tmax) ./ 365

lines!(above5ax_prev, times, ensemble_seir_arr[2, :, 1])
lines!(above5ax_inc, times, inc_infec_arr[1, :, 1])
outbreak_fig = barplot!(
    above5ax_periodsum,
    times,
    inc_infec_arr[3, :, 1];
    color = inc_infec_arr[4, :, 1],
    colormap = [:blue, :red],
)

map(hidexdecorations!, [above5ax_prev, above5ax_inc])

map(
    ax -> xlims!(ax, (92, 94)),
    [above5ax_prev, above5ax_inc, above5ax_periodsum],
)
ylims!(above5ax_periodsum, (0, 10000))
ylims!(above5ax_inc, (0, 300))

axislegend(
    above5ax_periodsum,
    [PolyElement(; color = col) for col in [:blue, :red]],
    ["Not Outbreak", "Outbreak"],
)

above5fig

#%%
################################################################################
########################## Background Noise ####################################
################################################################################
background_ode!(du, u, p, t) = (du .= 0.0)
background_noise!(du, u, p, t) = (du .= 1.0)

sde_condition(u, t, integrator) = true
function sde_affect!(integrator)
    if integrator.u[1] < 0.0
        integrator.u[1] = -integrator.u[1]
    end
end

sde_cb = DiscreteCallback(
    sde_condition, sde_affect!; save_positions = (false, false)
)

#%%
# Noise should be incidence, not prevalence
noise_u₀ = round.([0.5 * maximum(inc_infec_arr[1, 200:end, 1])])
tspan = (tlower, tmax)
noise_prob = SDEProblem(background_ode!, background_noise!, noise_u₀, tspan, p)
noise_sol = solve(
    noise_prob, SRIW1(); callback = sde_cb, dt = param_dict[:dt],
    adaptive = false,
)
noise_df = rename(DataFrame(noise_sol), [:time, :noise])

lines(noise_df[:, :time], noise_df[:, :noise])

#%%
noise_arr = zeros(
    Float64, 3, size(ensemble_jump_arr, 2), size(ensemble_jump_arr, 3)
)
@floop for sim in 1:size(ensemble_jump_arr, 3)
    noise_prob =
        noise_prob = SDEProblem(
            background_ode!,
            background_noise!,
            [0.5 * maximum(ensemble_seir_arr[2, 200:end, sim])],
            tspan,
            p,
        )

    noise_sol = solve(
        noise_prob, SRIW1(); callback = sde_cb, dt = param_dict[:dt],
        adaptive = false,
    )

    noise_arr[1, :, sim] = noise_sol[1, :]

    for day in 2:size(noise_arr, 2)
        noise_arr[2, day, sim] =
            noise_arr[1, day, sim] - noise_arr[1, day - 1, sim]
    end
end

#%% 
noise_fig = Figure()
noise_ax = Axis(
    noise_fig[1, 1]; xlabel = "Time (years)", ylabel = "Noise Prevalence"
)

for sim in eachrow(noise_arr)
    lines!(noise_ax, times, sim; color = (:red, 0.05))
end

noise_fig