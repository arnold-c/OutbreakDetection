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
β = calculate_beta(R₀, γ, μ, 1, N)
# Adjust the scale of the seasonal variation in infectivity i.e. A_scale scales the amplitude of cosine function to 1/A_scale. The maximum infectivity is β
βmax_scale = 2.0
βmin_scale = 0.5
ε = (1.06 * μ * (R₀ - 1)) / sqrt(N) # Commuter imports - see p210 Keeling & Rohani
p = (β, γ, μ, ε, R₀, βmax_scale, βmin_scale)

Random.seed!(1234)

#%%
function calculate_beta_amp(Amin, Amax, t)
    return 0.5 * (Amax - Amin) * cos(2pi * t / 365) + ((Amax - Amin) / 2 + Amin)
end

function sir_mod(u, p, trange; retβamp = false)
    tlength = length(trange)
    dt = step(trange)

    state_arr = zeros(Int64, size(u, 1), tlength)

    change_arr = zeros(Int64, size(u, 1), tlength)

    jump_arr = zeros(Int64, 7, tlength)

    if retβamp == true
        beta_amps = zeros(Float64, 2, tlength)
        sir_mod!(
            state_arr, change_arr, jump_arr, beta_amps, u, p, trange; dt = dt
        )
        return state_arr, change_arr, jump_arr, beta_amps
    else
        sir_mod!(state_arr, change_arr, jump_arr, u, p, trange; dt = dt)
        return state_arr, change_arr, jump_arr
    end
end

function sir_mod_loop!(state_arr, change_arr, jump_arr, β_scaled, j, dt)
    S = state_arr[1, j - 1]
    I = state_arr[2, j - 1]
    R = state_arr[3, j - 1]
    N = state_arr[4, j - 1]

    infec_rate = β_scaled * S * I
    infect_num = rand(Poisson(infec_rate * dt))

    recov_rate = γ * I
    recov_num = rand(Poisson(recov_rate * dt))

    birth_rate = μ * N
    birth_num = rand(Poisson(birth_rate * dt))

    S_death_rate = μ * S
    S_death_num = rand(Poisson(S_death_rate * dt))

    I_death_rate = μ * I
    I_death_num = rand(Poisson(I_death_rate * dt))

    R_death_rate = μ * R
    R_death_num = rand(Poisson(R_death_rate * dt))

    import_rate = ε * N / R₀
    import_num = rand(Poisson(import_rate * dt))

    change_arr[1, j] = birth_num - infect_num - S_death_num - import_num
    change_arr[2, j] = infect_num - recov_num - I_death_num + import_num
    change_arr[3, j] = recov_num - R_death_num
    change_arr[4, j] = sum(change_arr[1:3, j])

    for state in 1:size(state_arr, 1)
        if state_arr[state, j - 1] + change_arr[state, j] < 0
            change_arr[state, j] = -state_arr[state, j - 1]
        end
    end

    @. state_arr[:, j] = state_arr[:, j - 1] + change_arr[:, j]

    jump_arr[:, j] = [
        infect_num,
        recov_num,
        birth_num,
        S_death_num,
        I_death_num,
        R_death_num,
        import_num,
    ]
end

function sir_mod!(state_arr, change_arr, jump_arr, u, p, trange; dt)
    S0, I0, R0, N0 = u
    β, γ, μ, ε, R₀, βmax_scale, βmin_scale = p

    for (j, t) in pairs(trange)
        if j == 1
            state_arr[:, j] = u
            continue
        end

        β_amp = calculate_beta_amp(βmax_scale, βmin_scale, t)
        β_scaled = β * β_amp

        sir_mod_loop!(state_arr, change_arr, jump_arr, β_scaled, j, dt)
    end

    return nothing
end

function sir_mod!(state_arr, change_arr, jump_arr, beta_amps_arr, u, p, trange; dt)
    S0, I0, R0, N0 = u
    β, γ, μ, ε, R₀, βmax_scale, βmin_scale = p

    for (j, t) in pairs(trange)

        if j == 1
            state_arr[:, j] = u
            continue
        end

        β_amp = calculate_beta_amp(βmax_scale, βmin_scale, t)
        β_scaled = β * β_amp
        
        sir_mod_loop!(state_arr, change_arr, jump_arr, β_scaled, j, dt)

        beta_amps_arr[1, j] = β_amp
        beta_amps_arr[2, j] = β_scaled
    end

    return nothing
end

#%%
sir_array, change_array, jump_array, β_amps = sir_mod(u₀, p, trange; retβamp = true)
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
        facet = (; linkyaxes = :none), axis = (; limits = ((0, 1), nothing))
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
@chain DataFrame(Tables.table(β_amps')) begin
    hcat(trange, _)
    rename!([:time, :β_amp, :β_scaled])
    stack(_, Not("time"); variable_name = :β, value_name = :Number)
    data(_) *
    mapping(
        :time => (t -> t / 365) => "Time (years)", :Number;
        color = :β => sorter(["β_amp", "β_scaled"]), layout = :β,
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