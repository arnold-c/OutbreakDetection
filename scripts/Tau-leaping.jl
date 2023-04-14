"""
This is a simulation of an SIR model that uses Tau-leaping.
All jumps are manually defined.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics, DifferentialEquations, ModelingToolkit
using BenchmarkTools, JLD2, Random, ProgressMeter, StatsBase, Distributions

CairoMakie.activate!()
set_aog_theme!()

#%%
# Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/transmission-functions.jl"))
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
N = 5e5
s = 0.1
i = 0.01
r = 1.0 - (s + i)
u₀ = convert.(Int64, [N * s, N * i, N * r, N])
δt = 1.0
tlower = 0.0
tmax = 365.0 * 100
tspan = (tlower, tmax)
trange = tlower:δt:tmax
tlength = length(trange)

dur_inf = 8
R₀ = 10.0
γ = 1 / dur_inf
μ = 1 / (62 * 365)
β = calculate_beta(R₀, γ, μ, 1, N)
ε = (1.06 * μ * (R₀ - 1)) / sqrt(sum(u₀)) # Commuter imports - see p210 Keeling & Rohani
p = (β, γ, μ, ε, R₀, δt)

Random.seed!(1234)

#%%
function sir_mod(u, p, tlength)
    S0, I0, R0, N0 = u
    β, γ, μ, ε, R₀, δt = p

    states = zeros(4, tlength)
    changes = similar(states)
    states[:, 1] = [S0, I0, R0, N0]

    jumps = zeros(7, tlength)

    for i in 2:tlength
        S = states[1, i - 1]
        I = states[2, i - 1]
        R = states[3, i - 1]
        N = states[4, i - 1]

        infec_rate = β * S * I
        infect_num = rand(Poisson(infec_rate * δt))

        recov_rate = γ * I
        recov_num = rand(Poisson(recov_rate * δt))

        birth_rate = μ * N
        birth_num = rand(Poisson(birth_rate * δt))

        S_death_rate = μ * S
        S_death_num = rand(Poisson(S_death_rate * δt))

        I_death_rate = μ * I
        I_death_num = rand(Poisson(I_death_rate * δt))

        R_death_rate = μ * R
        R_death_num = rand(Poisson(R_death_rate * δt))

        import_rate = ε * N / R₀
        import_num = rand(Poisson(import_rate * δt))

        changes[1, i] = birth_num - infect_num - S_death_num - import_num
        changes[2, i] = infect_num - recov_num - I_death_num + import_num
        changes[3, i] = recov_num - R_death_num
        changes[4, i] = sum(changes[1:3, i])

        @. states[:, i] = states[:, i - 1] + changes[:, i]

        jumps[:, i] = [
            infect_num,
            recov_num,
            birth_num,
            S_death_num,
            I_death_num,
            R_death_num,
            import_num
        ]
    end

    return states, changes, jumps
end

#%%
sir_array, change_array, jump_array = sir_mod(u₀, p, tlength)
sir_df = create_sir_df(sir_array, trange, [:S, :I, :R, :N])


#%%
draw_sir_plot(sir_df; annual = true, labels = ["S", "I", "R", "N"])

#%%
@chain DataFrame(Tables.table(jump_array')) begin
    hcat(trange, _)
    rename!([
        "time", "Infect", "Recov", "Birth", "S_death", "I_death", "R_death"
    ])
    stack(_, Not("time"); variable_name = :Jump, value_name = :Number)
    data(_) *
    mapping(:time => "Time (days)", :Number; color = :Jump) *
    visual(Lines; linewidth = 4)
    draw
end

#%%
@chain DataFrame(Tables.table(change_array')) begin
    hcat(trange, _)
    rename!([
        "time", "dS", "dI", "dR", "dN"
    ])
    stack(_, Not("time"); variable_name = :Change, value_name = :Number)
    data(_) *
    mapping(:time => "Time (days)", :Number; color = :Change) *
    visual(Lines; linewidth = 4)
    draw
end