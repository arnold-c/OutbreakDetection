"""
This is a basic Gillespie simulation of an SIR model.
It uses the JumpProcesses.jl package and the Direct Gillespie method.
All jumps are constant rate jumps and are manually defined using JumpProcesses, not using the ModelingToolkit or Catalyst interfaces.
"""

using DrWatson
@quickactivate "OutbreakDetection"

using JumpProcesses, Statistics, DataFrames, DataFramesMeta, LinearAlgebra
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

u₀ = [999, 1, 0]
tspan = (0.0, 250.0)

function calculateR0(
        β::AbstractFloat,
        γ::AbstractFloat,
        μ::AbstractFloat,
        𝐂::Matrix,
        pop_matrix::Vector
    )
    size(𝐂, 1) == size(𝐂, 2) ? nothing : error("𝐂 must be square")
    size(𝐂, 1) == size(pop_matrix, 1) ? nothing : error("𝐂 and pop_matrix must have the same number of rows")

    𝚩 = β * 𝐂
    
    𝐅 = 𝚩 .* pop_matrix
    𝐕 = Diagonal(repeat([γ + μ], size(𝐂, 1)))

    𝐅𝐕⁻¹ = 𝐅 * inv(𝐕)
    eigenvals, eigenvectors = eigen(𝐅𝐕⁻¹)
    
    R₀ = maximum(real(eigenvals))
    
    return R₀
end

calculateR0(0.00025, 1/8, 0.0, ones(1, 1), [500])

function calculate_beta(
        R₀::AbstractFloat,
        γ::AbstractFloat,
        μ::AbstractFloat,
        𝐂::Matrix,
        pop_matrix::Vector
    )
    size(𝐂, 1) == size(𝐂, 2) ? nothing : error("𝐂 must be square")
    size(𝐂, 1) == size(pop_matrix, 1) ? nothing : error("𝐂 and pop_matrix must have the same number of rows")

    𝐅 = 𝐂 .* pop_matrix
    𝐕 = Diagonal(repeat([γ + μ], size(𝐂, 1)))

    𝐅𝐕⁻¹ = 𝐅 * inv(𝐕)
    eigenvals, eigenvectors = eigen(𝐅𝐕⁻¹)
    β = R₀ / maximum(real(eigenvals))
    
    return β
end

calculate_beta(2.0, 1/8, 0.0, ones(1, 1), [1000])

R₀ = 2.0
γ = 1/8
μ = 0.0

β = calculate_beta(R₀, γ, μ, ones(1, 1), [sum(u₀)])
p = (β, γ)

infec_rate(u, p, t) = p[1] * u[1] * u[2]  # β*S*I
function infec_affect!(integrator)
    integrator.u[1] -= 1         # S -> S - 1
    integrator.u[2] += 1         # I -> I + 1
    nothing
end
infec_jump = ConstantRateJump(infec_rate, infec_affect!)

recov_rate(u, p, t) = p[2] * u[2]         # ν*I
function recov_affect!(integrator)
    integrator.u[2] -= 1        # I -> I - 1
    integrator.u[3] += 1        # R -> R + 1
    nothing
end
recov_jump = ConstantRateJump(recov_rate, recov_affect!)


prob = DiscreteProblem(u₀, tspan, p)

jump_prob = JumpProblem(prob, Direct(), infec_jump, recov_jump)
sol = solve(jump_prob, SSAStepper())

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