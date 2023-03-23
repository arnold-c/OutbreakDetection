"""
This script just performs a simple SIR model using an ODE solver.
"""

using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

u₀ = [999.0, 1.0, 0.0]
tspan = (0.0, 250.0)

function sir!(du, u, p, t)
    (β, γ, μ) = p
    (S, I, R) = u
    N = S + I + R
    du[1] = -β * S * I + μ * (N - S)
    du[2] = β * S * I - γ * I - μ * I
    du[3] = γ * I - μ * R
    return nothing
end

R₀ = 2
γ = 1 / 8
μ = 1 / (62 * 365)
β = calculate_beta(R₀, γ, μ, 1, sum(u₀))

p = [β, γ, μ]

ode_prob = ODEProblem(sir!, u₀, tspan, p)
ode_sol = solve(ode_prob; saveat = 0.1)

ode_sol_df = create_sir_df(ode_sol)

colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(ode_sol_df; colors = colors)
