"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics
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
i = 0.01
r = 1.0 - (s + i)
u₀ = convert.(Int64, [N * s, N * i, N * r, N])
τ = 1.0
tlower = 0.0
tmax = 365.0 * 100
trange = tlower:τ:tmax
tlength = length(trange)

dur_inf = 8
R₀ = 10.0
γ = 1 / dur_inf
μ = 1 / (62.5 * 365)
# β₀ is the average transmission rate
β₀ = calculate_beta(R₀, γ, μ, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. β₁ scales the amplitude of cosine function
β₁ = 0.2
ε = (1.06 * μ * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
p = (β₀, β₁, γ, μ, ε, R₀)

Random.seed!(1234)

#%%
function calculate_beta_amp(β_mean, β_force, t)
    return β_mean * (1 + β_force * cos(2pi * t / 365))
end

function sir_mod(u, p, trange; retβamp = false, type = "stoch")
    tlength = length(trange)
    dt = step(trange)

    state_arr = zeros(Int64, size(u, 1), tlength)

    change_arr = zeros(Int64, size(u, 1), tlength)

    jump_arr = zeros(Int64, 7, tlength)

    if retβamp == true
        beta_arr = zeros(Float64, 1, tlength)
        sir_mod!(
            state_arr, change_arr, jump_arr, beta_arr, u, p, trange; dt = dt, type = type
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
    state_arr, change_arr, jump_arr, β_t, j, dt; type = type
)
    # Unpack the state variables for easier use
    S = state_arr[1, j - 1]
    I = state_arr[2, j - 1]
    R = state_arr[3, j - 1]
    N = state_arr[4, j - 1]

    # Calculate the rates of each event
    infec_rate = β_t * S * I
    recov_rate = γ * I
    birth_rate = μ * N
    S_death_rate = μ * S
    I_death_rate = μ * I
    R_death_rate = μ * R
    import_rate = ε * N / R₀

    # Calculate the number of jumps for each event
    if type == "stoch"
        jump_arr[:, j] = map(
            r -> rand(Poisson(r)),
            [infec_rate, recov_rate, birth_rate, S_death_rate, I_death_rate,
                R_death_rate, import_rate],
        )
    elseif type == "det"
        jump_arr[:, j] = map(
            r -> round(r * dt),
            [infec_rate, recov_rate, birth_rate, S_death_rate, I_death_rate,
                R_death_rate, import_rate],
        )
    else
        return("Type must be stoch or det")
    end

    # Calculate the change in each state
    change_arr[1, j] = jump_arr[3, j] - jump_arr[1, j] - jump_arr[4, j] - jump_arr[7, j]
    change_arr[2, j] = jump_arr[1, j] - jump_arr[2, j] - jump_arr[5, j] + jump_arr[7, j]
    change_arr[3, j] = jump_arr[2, j] - jump_arr[6, j]
    change_arr[4, j] = sum(change_arr[1:3, j])

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

function sir_mod!(state_arr, change_arr, jump_arr, u, p, trange; dt, type = "stoch")
    S0, I0, R0, N0 = u
    β_mean, β_force, γ, μ, ε, R₀ = p

    for (j, t) in pairs(trange)
        if j == 1
            state_arr[:, j] = u
            continue
        end

        β_t = calculate_beta_amp(β_mean, β_force, t)

        sir_mod_loop!(state_arr, change_arr, jump_arr, β_t, j, dt; type = type)
    end

    return nothing
end

function sir_mod!(state_arr, change_arr, jump_arr, beta_arr, u, p, trange; dt, type = "stoch")
    S0, I0, R0, N0 = u
    β_mean, β_force, γ, μ, ε, R₀ = p

    for (j, t) in pairs(trange)
        β_t = calculate_beta_amp(β_mean, β_force, t)
        beta_arr[j] = β_t

        if j == 1
            state_arr[:, j] = u

            continue
        end

        sir_mod_loop!(state_arr, change_arr, jump_arr, β_t, j, dt; type = type)
    end

    return nothing
end

#%%
sir_array, change_array, jump_array, β_arr = sir_mod(
    u₀, p, trange; retβamp = true, type = "det"
)
sir_df = create_sir_df(sir_array, trange, [:S, :I, :R, :N])

sircolors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]
state_labels = ["S", "I", "R", "N"]
draw_sir_plot(sir_df; annual = true, labels = state_labels)

#%%
@chain DataFrame(Tables.table(jump_array')) begin
    hcat(trange, _)
    rename!([
        "time", "Infect", "Recov", "Birth", "S_death", "I_death", "R_death",
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
change_labels = ["dS", "dI", "dR", "dN"]
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
        palettes = (; color = sircolors),
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
ensemble_sir_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
ensemble_change_arr = zeros(Int64, size(u₀, 1), tlength, nsims)
ensemble_jump_arr = zeros(Int64, 7, tlength, nsims)

for k in 1:nsims
    @views sir = ensemble_sir_arr[:, :, k]
    @views change = ensemble_change_arr[:, :, k]
    @views jump = ensemble_jump_arr[:, :, k]

    sir_mod!(sir, change, jump, u₀, p, trange; dt = τ)
end

quantile_ints = [95, 90, 80, 50]

ensemble_summary = zeros(Float64, 3, tlength, length(u₀), length(quantile_ints))

for q in eachindex(quantile_ints)
    qlow = round(0.5 - quantile_ints[q] / 200; digits = 3)
    qhigh = round(0.5 + quantile_ints[q] / 200; digits = 3)
    quantiles = [qlow, 0.5, qhigh]

    @views quantile_array = ensemble_summary[:, :, :, q]
    create_sir_all_sim_quantiles!(
        ensemble_sir_arr, quantile_array; quantiles = quantiles
    )
end

create_sir_quantiles_plot(
    ensemble_summary[:, :, :, 1]; annual = true, δt = τ
)

#%%
μ_min = 10
μ_max = 60
μ_step = 0.1
n_μs = length(μ_min:μ_step:μ_max)
μ_vec = zeros(Float64, n_μs)
μ_vec .= collect(μ_min:μ_step:μ_max) ./ (1000 * 365)

bifurc_sir_arr = zeros(Int64, size(u₀, 1), tlength, n_μs);
bifurc_change_arr = zeros(Int64, size(u₀, 1), tlength, n_μs);
bifurc_jump_arr = zeros(Int64, 7, tlength, n_μs);

@showprogress for (k, μ) in pairs(μ_vec)
    @views sir = bifurc_sir_arr[:, :, k]
    @views change = bifurc_change_arr[:, :, k]
    @views jump = bifurc_jump_arr[:, :, k]

    p = (β₀, β₁, γ, μ, ε, R₀)

    sir_mod!(sir, change, jump, u₀, p, trange; dt = τ)
end

bifurc_summary = DataFrame(
    Dict(
        :μ => μ_min:μ_step:μ_max,
        :S => [median(bifurc_sir_arr[1, end, :, :]; dims = 2)...],
        :I => [median(bifurc_sir_arr[2, end, :, :]; dims = 2)...],
        :R => [median(bifurc_sir_arr[3, end, :, :]; dims = 2)...],
        :N => [median(bifurc_sir_arr[4, end, :, :]; dims = 2)...],
    ),
)

@chain bifurc_summary begin
    data(_) *
    mapping(:μ, :I) *
    visual(Scatter)
    draw
end

@chain bifurc_summary begin
    data(_) *
    mapping(:μ, :N) *
    visual(Scatter)
    draw
end

#%%
bifurc_fig = Figure()
bifurc_ax = Axis(bifurc_fig[1, 1])

for sim in 1:bifurc_sims
    scatter!(
        bifurc_ax,
        μ_min:μ_step:μ_max,
        bifurc_sir_arr[2, end, :, sim];
        colorrange = (1, 10),
        color = :black,
    )
end

bifurc_fig