"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using WGLMakie, AlgebraOfGraphics, ColorSchemes, Colors
using DifferentialEquations, ModelingToolkit
using BenchmarkTools, JLD2, Random, ProgressMeter, StatsBase, Distributions
using IterTools, FLoops, FreqTables, ThreadsX, ProtoStructs

WGLMakie.activate!()
#= CairoMakie.activate!(type = "pdf") =#
set_aog_theme!()
# Set depending on size of screen
update_theme!(; resolution = (1300, 900))

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
seir_array, change_array, jump_array, β_arr = seir_mod(
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

    local p = (β₀, β₁, σ, γ, μ, ε, R₀)

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
        global sim_data = load(sim_files[i])
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
function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    outbreakuppers = outbreakaccum[findall(==(1), outbreakrle[1])]
    outbreaklowers = filter(
        x -> x <= maximum(outbreakuppers),
        outbreakaccum[findall(==(0), outbreakrle[1])] .+ 1,
    )

    return (outbreaklowers, outbreakuppers)
end

#%%
outbreakthreshold = 5
minoutbreakdur = 30
minoutbreaksize = 500

function create_inc_infec_arr!(
    incarr, ensemblejumparr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    prog = Progress(size(ensemblejumparr, 3))
    @floop for sim in axes(ensemblejumparr, 3)
        # Copy new infections to array
        incarr[:, 1, sim] = @view(ensemblejumparr[1, :, sim])
        # Calculate if new infection is above or below threshold
        incarr[:, 2, sim] =
            @view(incarr[:, 1, sim]) .>= outbreakthreshold

        # Calculate the total number of infections above threshold in a consecutive string of days
        # Calculate the number of consecutive days of infection above or below threshold
        above5rle = rle(@view(incarr[:, 2, sim]))

        ## Calculate upper and lower indices of consecutive days of infection
        above5lowers, above5uppers = calculate_outbreak_thresholds(above5rle)

        for (lower, upper) in zip(above5lowers, above5uppers)
            # Calculate number of infections between lower and upper indices
            period_sum = sum(@view(incarr[lower:upper, 1, sim]))
            incarr[lower:upper, 3, sim] .= period_sum

            # Determine if there is an outbreak between lower and upper indices
            if upper - lower >= minoutbreakdur && period_sum >= minoutbreaksize
                incarr[lower:upper, 4, sim] .= 1
            end
        end

        next!(prog)
    end
end

function create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)
    incarr = zeros(
        Int64, size(ensemble_jump_arr, 2), 4, size(ensemble_jump_arr, 3)
    )

    create_inc_infec_arr!(
        incarr,
        ensemble_jump_arr,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize,
    )

    return incarr
end

inc_infec_arr = create_inc_infec_arr(
    ensemble_jump_arr, outbreakthreshold, minoutbreakdur, minoutbreaksize
)

#%%
outbreakcols = [ColorSchemes.magma[i] for i in (200, 20)]

above5fig = Figure()
above5ax_prev = Axis(above5fig[1, 1]; ylabel = "Prevalence")
above5ax_inc = Axis(above5fig[2, 1]; ylabel = "Incidence")
above5ax_periodsum = Axis(
    above5fig[3, 1]; xlabel = "Time (years)", ylabel = "Period Sum"
)

linkxaxes!(above5ax_prev, above5ax_inc, above5ax_periodsum)

times = collect(0:param_dict[:dt]:tmax) ./ 365

lines!(above5ax_prev, times, ensemble_seir_arr[2, :, 1])
lines!(above5ax_inc, times, inc_infec_arr[:, 1, 1])
outbreak_fig = barplot!(
    above5ax_periodsum,
    times,
    inc_infec_arr[:, 3, 1];
    color = inc_infec_arr[:, 4, 1],
    colormap = outbreakcols,
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
    [PolyElement(; color = col) for col in outbreakcols],
    ["Not Outbreak", "Outbreak"],
)

above5fig

#%%
################################################################################
########################## Background Noise ####################################
################################################################################
background_ode!(du, u, p, t) = (du .= 0.0)
background_noise!(du, u, p, t) = (du .= 0.1)

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
noise_u₀ = [10.0]
tspan = (tlower, tmax)
noise_prob = SDEProblem(background_ode!, background_noise!, noise_u₀, tspan, p)
noise_sol = solve(
    noise_prob, SRIW1(); callback = sde_cb, dt = param_dict[:dt],
    adaptive = false,
)
noise_df = rename(DataFrame(noise_sol), [:time, :noise])

lines(noise_df[:, :time], noise_df[:, :noise])

#%%
# Set noise arr to 3D array (even though not necessary), so it has the same 
# dimensions as the other arrays
noise_arr = zeros(
    Float64, size(ensemble_jump_arr, 2), 1, size(ensemble_jump_arr, 3)
)
@floop for sim in axes(ensemble_jump_arr, 3)
    noise_prob = SDEProblem(
        background_ode!,
        background_noise!,
        noise_u₀,
        tspan,
        p
    )

    noise_sol = solve(
        noise_prob, SRIW1(); callback = sde_cb, dt = param_dict[:dt],
        adaptive = false,
    )

    noise_arr[:, 1, sim] = @view(noise_sol[1, :])
    # Set first noise incidence to 0 as no new noise in the first time step
    noise_arr[1, 1, sim] = 0.0
end

#%% 
# noise_fig = Figure()
# noise_ax = Axis(
#     noise_fig[1, 1]; xlabel = "Time (years)", ylabel = "Noise Incidence"
# )

# for sim in axes(noise_arr, 3)
#     lines!(noise_ax, times, noise_arr[:, 1, sim]; color = (:red, 0.1))
# end

# noise_fig

#%%
################################################################################
############################### Testing ########################################
################################################################################
testlag = 3
moveavglag = 7
detectthreshold = 10
perc_clinic = 0.3
perc_clinic_test = 0.8
perc_tested = perc_clinic * perc_clinic_test
testsens = 0.9
testspec = 0.9

testing_arr = zeros(Int64, tlength, 8, size(inc_infec_arr, 3));
post_odds_arr = zeros(Float64, tlength, 2, size(inc_infec_arr, 3));

function calculate_tested!(outarr, outarr_ind, inarr, perc_tested, sim)
    @. outarr[:, outarr_ind, sim] = round(@view(inarr[:, 1, sim]) * perc_tested)
end

function calculate_pos!(
    npos_vec,
    tested_vec,
    ntested,
    lag,
    sens,
    spec;
    noise = false
)
    if noise
        for day in eachindex(tested_vec)
            if day + lag <= ntested
                npos_vec[day + lag] = Int64(
                    round(tested_vec[day] * (1.0 - spec))
                )
            end
        end
    else
        for day in eachindex(tested_vec)
            if day + lag <= ntested
                npos_vec[day + lag] = Int64(
                    round(tested_vec[day] * sens)
                )
            end
        end
    end

    return nothing
end

function calculate_pos(
    tested_vec,
    lag,
    sens,
    spec;
    noise = false,
)
    ntested = length(tested_vec)
    npos = zeros(Int64, ntested)

    calculate_pos!(
        npos,
        tested_vec,
        ntested,
        lag,
        sens,
        spec;
        noise = noise,
    )

    return npos
end

#%%
function calculate_movingavg!(invec, outvec, testlag, avglag; Float = true)
    if Float
        avgfunc =
            (invec, day, avglag) -> mean(@view(invec[(day - avglag + 1):day]))
    else
        avgfunc =
            (invec, day, avglag) -> Int64(round(
                mean(@view(invec[(day - avglag + 1):day]))
            ))
    end
    for day in eachindex(invec)
        if day >= testlag + avglag + 1
            outvec[day] = avgfunc(invec, day, avglag)
        end
    end
end

function calculate_movingavg(invec, testlag, avglag)
    outvec = zeros(Float64, size(invec, 1), 1)

    calculate_movingavg!(invec, outvec, testlag, avglag)

    return outvec
end

function detectoutbreak!(outbreakvec, incvec, avgvec, threshold, avglag)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function detectoutbreak(incvec, avgvec, threshold, avglag)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold, avglag)

    return outbreak
end

#%%
function create_testing_arr!(
    testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
    testspec,
    detectthreshold, moveavglag,
)
    ntested = size(testarr, 1)

    # prog = Progress(size(incarr, 3))
    @floop for sim in 1:size(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(testarr, 1, incarr, perc_tested, sim)

        # Number of noise individuals tested
        calculate_tested!(testarr, 2, noisearr, perc_tested, sim)

        # Number of test positive INFECTED individuals
        calculate_pos!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            ntested,
            testlag,
            testsens,
            testspec;
            noise = false,
        )

        # Number of test positive NOISE individuals
        calculate_pos!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            ntested,
            testlag,
            testsens,
            testspec;
            noise = true,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            testlag, moveavglag;
            Float = false,
        )

        # TOTAL Test positive individuals trigger outbreak response 
        detectoutbreak!(
            @view(testarr[:, 7, sim]),
            @view(testarr[:, 5, sim]),
            @view(testarr[:, 6, sim]),
            detectthreshold, moveavglag,
        )

        # # Posterior prob of infectious / total test positive
        @. @view(posoddsarr[:, 1, sim]) =
            @view(testarr[:, 3, sim]) / @view(testarr[:, 5, sim])
        calculate_movingavg!(
            @view(posoddsarr[:, 1, sim]),
            @view(posoddsarr[:, 2, sim]),
            testlag, moveavglag,
        )

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 8, sim] =
            @view(testarr[:, 7, sim]) == @view(incarr[:, 4, sim])

        # next!(prog)
    end

    return nothing
end

function create_testing_arr(
    incarr, noisearr, perc_tested, testlag, testsens, testspec, detectthreshold,
    moveavglag,
)
    testarr = zeros(Int64, size(incarr, 1), 6, size(incarr, 3))
    posoddsarr = zeros(Float64, size(incarr, 1), 2, size(incarr, 3))

    create_testing_arr!(
        testarr, incarr, noisearr, posoddsarr, perc_tested, testlag, testsens,
        testspec,
        detectthreshold, moveavglag,
    )

    return testarr
end

#%%
# @benchmark create_testing_arr!(
#     $testing_arr,
#     $inc_infec_arr,
#     $noise_arr,
#     $post_odds_arr,
#     $perc_tested,
#     $testlag,
#     $testsens,
#     $testspec,
#     $detectthreshold,
#     $moveavglag,
# )

create_testing_arr!(
    testing_arr,
    inc_infec_arr,
    noise_arr,
    post_odds_arr,
    perc_tested,
    testlag,
    testsens,
    testspec,
    detectthreshold,
    moveavglag,
)

#%%
inc_test_fig = Figure()
inc_test_ax1 = Axis(inc_test_fig[1, 1]; ylabel = "Incidence")
inc_test_ax2 = Axis(inc_test_fig[2, 1]; ylabel = "Test Positive")
inc_test_ax3 = Axis(
    inc_test_fig[3, 1];
    xlabel = "Time (years)",
    ylabel = "7d Avg Test Positive"
)

lines!(
    inc_test_ax1, times, inc_infec_arr[:, 1, 1];
    color = inc_infec_arr[:, 4, 1],
    colormap = outbreakcols,
)
lines!(
    inc_test_ax2, times, testing_arr[:, 3, 1];
    color = testing_arr[:, 7, 1],
    colormap = outbreakcols,
)
lines!(
    inc_test_ax3, times, testing_arr[:, 6, 1];
    color = testing_arr[:, 7, 1],
    colormap = outbreakcols,
)

linkxaxes!(inc_test_ax1, inc_test_ax2, inc_test_ax3)

map(hidexdecorations!, [inc_test_ax1, inc_test_ax2])

# map(
#     ax -> xlims!(ax, (1750 / 365, 1850 / 365)),
#     [inc_test_ax1, inc_test_ax2, inc_test_ax3],
# )
map(ax -> ylims!(ax, (0, 20)), [inc_test_ax1, inc_test_ax2, inc_test_ax3])

hlines!(
    inc_test_ax1, 5;
    color = :black,
    linestyle = :dash,
    linewidth = 2,
)

map(
    ax -> hlines!(
        ax,
        detectthreshold;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    ),
    [inc_test_ax2, inc_test_ax3],
)

Legend(
    inc_test_fig[1, 2],
    [PolyElement(; color = col) for col in outbreakcols],
    ["Not Outbreak", "Outbreak"],
    "True\nOutbreak Status",
)

Legend(
    inc_test_fig[2:3, 2],
    [PolyElement(; color = col) for col in outbreakcols],
    ["Not Outbreak", "Outbreak"],
    "Detected\nOutbreak Status",
)

inc_test_fig

#%%
testing_fig = Figure()
plot_test_sims = 4
plot_test_coords = 1:(plot_test_sims ÷ 2)
for (sim, ax) in
    zip(1:plot_test_sims, IterTools.product(plot_test_coords, plot_test_coords))
    row = ax[1]
    col = ax[2]

    fig_ax = Symbol("testing_ax_" * string(sim))

    @eval $(fig_ax) = Axis(
        testing_fig[$row, $col]; xlabel = "Time (years)",
        ylabel = "Tested"
    )

    @eval lines!(
        $(fig_ax), times, testing_arr[:, 1, $sim];
        color = :red,
        label = "Infectious",
    )
    @eval lines!(
        $(fig_ax), times, testing_arr[:, 2, $sim];
        color = :blue,
        label = "Noise",
    )
    @eval lines!(
        $(fig_ax), times, testing_arr[:, 5, $sim];
        color = :black,
        label = "Total Positive",
    )
end

linkxaxes!(testing_ax_1, testing_ax_2)
linkxaxes!(testing_ax_3, testing_ax_4)

linkyaxes!(testing_ax_1, testing_ax_3)
linkyaxes!(testing_ax_2, testing_ax_4)

map(hidexdecorations!, [testing_ax_1, testing_ax_3])
map(hideydecorations!, [testing_ax_3, testing_ax_4])

Legend(
    testing_fig[3, :],
    testing_ax_1,
    "Type of Individual";
    orientation = :horizontal,
)

testing_fig

#%%
@proto struct OutbreakThresholdChars{A,B,C,D}
    crosstab::A
    tp::B
    tn::B
    fp::B
    fn::B
    sens::C
    spec::C
    ppv::C
    npv::C
    noutbreaks::B
    ndetectoutbreaks::B
    outbreakbounds::D
    detectoutbreakbounds::D
end

#%%
function calculate_ot_characterstics(test_arr, infec_arr, ind)
    crosstab = freqtable(testing_arr[:, 5, ind], inc_infec_arr[:, 4, ind])

    tp = crosstab[2, 2]
    tn = crosstab[1, 1]
    fp = crosstab[2, 1]
    fn = crosstab[1, 2]

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return crosstab, tp, tn, fp, fn, sens, spec, ppv, npv
end

function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

#%%
OT_chars = ThreadsX.map(
    axes(inc_infec_arr, 3)
) do sim
    outbreakrle = rle(@view(inc_infec_arr[:, 4, sim]))
    detectrle = rle(@view(testing_arr[:, 7, sim]))

    OutbreakThresholdChars(
        calculate_ot_characterstics(testing_arr, inc_infec_arr, sim)...,
        calculate_noutbreaks(outbreakrle),
        calculate_noutbreaks(detectrle),
        reduce(hcat, collect(calculate_outbreak_thresholds(outbreakrle))),
        reduce(hcat, collect(calculate_outbreak_thresholds(detectrle))),
    )
end

#%%
OT_chars[1].crosstab
OT_chars[1].outbreakbounds
OT_chars[1].detectoutbreakbounds
OT_chars[1].noutbreaks
OT_chars[1].ndetectoutbreaks

# Note that an outbreak isn't detected continously!
testing_arr[80:100, :, 1]

#%%
otchars_vec = zeros(Float64, length(OT_chars), 6);

@floop for sim in eachindex(OT_chars)
    otchars_vec[sim, 1] = OT_chars[sim].sens
    otchars_vec[sim, 2] = OT_chars[sim].spec
    otchars_vec[sim, 3] = OT_chars[sim].ppv
    otchars_vec[sim, 4] = OT_chars[sim].npv
    otchars_vec[sim, 5] = OT_chars[sim].noutbreaks
    otchars_vec[sim, 6] = OT_chars[sim].ndetectoutbreaks
end

#%%
outbreak_dist_fig = Figure()
outbreak_dist_ax = Axis(
    outbreak_dist_fig[1, 1]; xlabel = "Proportion of Time Series with Outbreak"
)

hist!(
    outbreak_dist_ax,
    vec(sum(@view(inc_infec_arr[:, 4, :]); dims = 1)) ./ size(inc_infec_arr, 1);
    bins = 0.0:0.01:0.7,
    color = (:blue, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "True Outbreaks",
)

hist!(
    outbreak_dist_ax,
    vec(sum(@view(testing_arr[:, 7, :]); dims = 1)) ./ size(testing_arr, 1);
    bins = 0.0:0.01:0.7,
    color = (:red, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "Tested Outbreaks",
)

Legend(outbreak_dist_fig[1, 2], outbreak_dist_ax, "Outbreak Proportion")

outbreak_dist_fig

#%%
noutbreaks_fig = Figure()
noutbreaks_ax = Axis(noutbreaks_fig[1, 1]; xlabel = "Number of Outbreaks")

hist!(
    noutbreaks_ax,
    @view(otchars_vec[:, 5]);
    bins = 0.0:10.0:450.0,
    color = (:blue, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "True Outbreaks",
)

hist!(
    noutbreaks_ax,
    @view(otchars_vec[:, 6]);
    bins = 0.0:10.0:450.0,
    color = (:red, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    normalization = :pdf,
    label = "Tested Outbreaks",
)

Legend(noutbreaks_fig[1, 2], noutbreaks_ax, "# Outbreaks")

noutbreaks_fig

#%%
sens_spec_fig = Figure()
sens_spec_ax = Axis(sens_spec_fig[1, 1]; xticks = 0.0:0.1:1.0)

hist!(
    sens_spec_ax,
    @view(otchars_vec[:, 1]);
    bins = 0.0:0.01:1.01,
    color = (:blue, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "Sensitivity",
    normalization = :pdf,
)

hist!(
    sens_spec_ax,
    @view(otchars_vec[:, 2]);
    bins = 0.0:0.01:1.01,
    color = (:red, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "Specificity",
    normalization = :pdf,
)

vlines!(
    sens_spec_ax,
    [mean(@view(otchars_vec[:, i])) for i in 1:2];
    color = :black,
    linestyle = :dash,
    linewidth = 4,
)

Legend(sens_spec_fig[1, 2], sens_spec_ax, "Characteristic")

sens_spec_fig

#%%
ppv_npv_fig = Figure()
ppv_npv_ax = Axis(ppv_npv_fig[1, 1]; xticks = 0.0:0.1:1.0)

hist!(
    ppv_npv_ax,
    @view(otchars_vec[:, 3]);
    bins = 0.0:0.01:1.01,
    color = (:green, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "PPV",
    normalization = :pdf,
)

hist!(
    ppv_npv_ax,
    @view(otchars_vec[:, 4]);
    bins = 0.0:0.01:1.01,
    color = (:purple, 0.5),
    strokecolor = :black,
    strokewidth = 1,
    label = "NPV",
    normalization = :pdf,
)

vlines!(
    ppv_npv_ax,
    [mean(@view(otchars_vec[:, i])) for i in 3:4];
    color = :black,
    linestyle = :dash,
    linewidth = 4,
)

Legend(ppv_npv_fig[1, 2], ppv_npv_ax, "Characteristic")

ppv_npv_fig

#%%
# TODO: Calculate threshold chars over range of test lags and detection thresholds. Use constant R0 and ind test chars for the moment.

@proto struct MeanOutbreakThresholdChars{A,B}
    testlag::A
    incthreshold::A
    tp::A
    tn::A
    fp::A
    fn::A
    sens::B
    spec::B
    ppv::B
    npv::B
    noutbreaks::A
    ndetectoutbreaks::A
end

testlag_vec = 0:1:2
detectthreshold_vec = [1, collect(5:5:10)...]

function create_OTchars_struct(incarr, testarr)
    ThreadsX.map(
        axes(incarr, 3)
    ) do sim
        outbreakrle = rle(@view(incarr[:, 4, sim]))
        return detectrle = rle(@view(testarr[:, 5, sim]))

        OutbreakThresholdChars(
            calculate_ot_characterstics(testarr, incarr, sim)...,
            calculate_noutbreaks(outbreakrle),
            calculate_noutbreaks(detectrle),
            reduce(hcat, collect(calculate_outbreak_thresholds(outbreakrle))),
            reduce(hcat, collect(calculate_outbreak_thresholds(detectrle))),
        )
    end
end

function calculate_mean_ot_chars(
    testlag,
    detectthreshold;
    ensemblejumparr = ensemble_jump_arr,
    outbreakthreshold = outbreakthreshold,
    minoutbreakdur = minoutbreakdur,
    minoutbreaksize = minoutbreaksize,
    noisearr = noise_arr,
    perctested = perc_tested,
    testsens = testsens,
    testspec = testspec,
    moveavglag = moveavglag,
)
    @info "Creating Incidence Array"
    incarr = create_inc_infec_arr(
        ensemble_jump_arr,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize
    )

    @info "Creating Testing Array"
    testing_arr = zeros(
        Int64, size(incarr, 1), 6, size(incarr, 3)
    )

    create_testing_arr!(
        testing_arr,
        incarr,
        noise_arr,
        perc_tested,
        testlag,
        testsens,
        testspec,
        detectthreshold,
        moveavglag,
    )

    @info "Calculating OT characteristics"
    OT_chars = create_OTchars_struct(incarr, testing_arr)

    return OT_chars
end

#%%
test = calculate_mean_ot_chars(7, 10)

#%%
for (testlag, detectthreshold) in
    IterTools.product(testlag_vec, detectthreshold_vec)
end
