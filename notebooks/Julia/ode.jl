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
    (β, γ) = p
    (S, I, R) = u
    du[1] = -β * S * I
    du[2] = β * S * I - γ * I
    du[3] = γ * I
    nothing
end


R₀ = 2
γ = 1/8
β = R₀ * γ / sum(u₀)

p = [β, γ]

prob = ODEProblem(sir!, u₀, tspan, p)
sol = solve(prob)


sol_df = @chain Tables.table(sol') begin
    DataFrame(_)
    rename!(:Column1 => :S, :Column2 => :I, :Column3 => :R)
    @rtransform! :N = :S + :I + :R
    hcat(_, DataFrame(time = sol.t))
    stack(_, [:S, :I, :R, :N], variable_name = :State, value_name = :Number)
end

colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]
sir_plot = data(sol_df) *
    mapping(
        :time => "Time (days)", :Number,
        color = :State => sorter("S", "I", "R", "N")
    ) *
    visual(Lines, linewidth = 4)

draw(
    sir_plot;
    palettes = (; color = colors),
    # axis = (; limits = ((0.0, 2.0), nothing))
    )