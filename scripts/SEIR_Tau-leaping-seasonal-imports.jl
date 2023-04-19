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
μ = 100 / (1000 * 365)
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

function sir_mod(u, p, trange; retβamp = false, type = "stoch")
    tlength = length(trange)
    dt = step(trange)

    state_arr = zeros(Float64, size(u, 1), tlength)

    change_arr = zeros(Float64, size(u, 1), tlength)

    jump_arr = zeros(Float64, 9, tlength)

    if retβamp == true
        beta_arr = zeros(Float64, 1, tlength)
        sir_mod!(
            state_arr, change_arr, jump_arr, beta_arr, u, p, trange; dt = dt,
            type = type,
        )
        return state_arr, change_arr, jump_arr, beta_arr
    else
        sir_mod!(
            state_arr, change_arr, jump_arr, u, p, trange; dt = dt, type = type
        )
        return state_arr, change_arr, jump_arr
    end
end

function sir_mod_loop!(
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

function sir_mod!(
    state_arr, change_arr, jump_arr, u, p, trange; dt, type = "stoch"
)
    S0, I0, R0, N0 = u

    for (j, t) in pairs(trange)
        if j == 1
            state_arr[:, j] = u
            continue
        end

        sir_mod_loop!(state_arr, change_arr, jump_arr, j, p, t, dt; type = type)
    end

    return nothing
end

function sir_mod!(
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

        sir_mod_loop!(state_arr, change_arr, jump_arr, j, p, t, dt; type = type)
    end

    return nothing
end

#%%
sir_array, change_array, jump_array, β_arr = sir_mod(
    u₀, p, trange; retβamp = true, type = "stoch"
)
sir_df = create_sir_df(sir_array, trange, [:S, :E, :I, :R, :N])

seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
state_labels = ["S", "E", "I", "R", "N"]
draw_sir_plot(
    sir_df;
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
@chain DataFrame(Tables.table(sir_array')) begin
    hcat(trange, _)
    rename!(["time", state_labels...])
    data(_) *
    mapping(:I, :S; color = :time) *
    visual(Lines)
    draw
end

#%%
@chain DataFrame(Tables.table(β_arr')) begin
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
nsims = 100
ensemble_seir_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
ensemble_change_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
ensemble_jump_arr = zeros(Int64, 9, tlength, nsims)

for k in 1:nsims
    @views seir = ensemble_seir_arr[:, :, k]
    @views change = ensemble_change_arr[:, :, k]
    @views jump = ensemble_jump_arr[:, :, k]

    sir_mod!(seir, change, jump, u₀, p, trange; dt = τ)
end

quantile_ints = [95, 90, 80, 50]

ensemble_summary = zeros(Float64, 3, tlength, length(u₀), length(quantile_ints))

for q in eachindex(quantile_ints)
    qlow = round(0.5 - quantile_ints[q] / 200; digits = 3)
    qhigh = round(0.5 + quantile_ints[q] / 200; digits = 3)
    quantiles = [qlow, 0.5, qhigh]

    @views quantile_array = ensemble_summary[:, :, :, q]
    create_sir_all_sim_quantiles!(
        ensemble_seir_arr, quantile_array; quantiles = quantiles
    )
end

create_sir_quantiles_plot(
    ensemble_summary[:, :, :, 1]; annual = true, δt = τ, colors = seircolors,
    labels = state_labels,
)

#%%
μ_min = 10
μ_max = 100
μ_step = 1.0
n_μs = length(μ_min:μ_step:μ_max)
μ_vec = zeros(Float64, n_μs)
μ_vec .= collect(μ_min:μ_step:μ_max) ./ (1000 * 365)

bifurc_μ_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_μs);
bifurc_μ_change_arr = zeros(Float64, size(u₀, 1), tlength, n_μs);
bifurc_μ_jump_arr = zeros(Float64, 9, tlength, n_μs);

@showprogress for (k, μ_run) in pairs(μ_vec)
    seir = @view bifurc_μ_seir_arr[:, :, k]
    change = @view bifurc_μ_change_arr[:, :, k]
    jump = @view bifurc_μ_jump_arr[:, :, k]

    params = (β₀, β₁, σ, γ, μ_run, ε, R₀)
    sir_mod!(seir, change, jump,
        u₀, params, trange; dt = τ, type = "det"
        )
end

years = (40 * 365):365:(tlength - 365)
bifurc_μ_annual_summary = zeros(Float64, length(years), n_μs, 5)
for (i, year) in pairs(years), state in 1:5, μ in eachindex(μ_vec)
    bifurc_μ_annual_summary[i, μ, state] = maximum(bifurc_μ_seir_arr[state, year:(year + 364), μ])
end

bifurc_μ_seir_arr[2, 40*365:(40*365 + 364), 1] == bifurc_μ_seir_arr[2, 40*365:(40*365 + 364), 91]

#%%
bifurc_μ_fig = Figure()
bifurc_μ_ax = Axis(bifurc_μ_fig[1, 1], xlabel = "μ (per 1000, per annum)", ylabel = "Max. I")

for year in eachindex(years)
    scatter!(
        bifurc_μ_ax,
        μ_min:μ_step:μ_max,
        bifurc_μ_annual_summary[year, :, 2];
        markersize = 4,
        color = :black
    )
end

bifurc_μ_fig

#%%
β₁_min = 0.0
β₁_max = 1.0
β₁_step = 0.01
n_β₁s = length(β₁_min:β₁_step:β₁_max)
β₁_vec = zeros(Float64, n_β₁s)
β₁_vec .= collect(β₁_min:β₁_step:β₁_max)

bifurc_β₁_seir_arr = zeros(Float64, size(u₀, 1), tlength, n_β₁s);
bifurc_β₁_change_arr = zeros(Float64, size(u₀, 1), tlength, n_β₁s);
bifurc_β₁_jump_arr = zeros(Float64, 9, tlength, n_β₁s);

μ = 20 / (1000 * 365)

@showprogress for (k, β₁) in pairs(β₁_vec)
    seir = @view bifurc_β₁_seir_arr[:, :, k]
    change = @view bifurc_β₁_change_arr[:, :, k]
    jump = @view bifurc_β₁_jump_arr[:, :, k]

    p = (β₀, β₁, σ, γ, μ, ε, R₀)

    sir_mod!(seir, change, jump, u₀, p, trange; dt = τ, type = "det")
end

years = (40 * 365):365:(tlength - 365)
bifurc_β₁_annual_summary = zeros(Float64, length(years), n_β₁s, 5)
for (i, year) in pairs(years), state in 1:5, (k, β₁) in pairs(β₁_vec)
    bifurc_β₁_annual_summary[i, k, state] = maximum(bifurc_β₁_seir_arr[state, year:(year + 364), k])
end

bifurc_β₁_seir_arr[2, 40*365:(40*365 + 364), 1] == bifurc_β₁_seir_arr[2, 40*365:(40*365 + 364), 91]

#%%
bifurc_β₁_fig = Figure()
bifurc_β₁_ax = Axis(bifurc_β₁_fig[1, 1], xlabel = "β₁ (seasonality)", ylabel = "Max. I")

for year in eachindex(years)
    scatter!(
        bifurc_β₁_ax,
        β₁_min:β₁_step:β₁_max,
        bifurc_β₁_annual_summary[year, :, 2];
        markersize = 4,
    )
end

bifurc_β₁_fig