"""
This is a deterministic SIR model using the ModelingToolkit interface.
The final section of code calculates the R₀ value, as well as the β parameter, for the model, using the Next Generation Matrix method.
"""
#%%
using DrWatson
@quickactivate "OutbreakDetection"

using DifferentialEquations, ModelingToolkit, Statistics
using DataFrames, DataFramesMeta, LinearAlgebra, Symbolics
using CairoMakie, AlgebraOfGraphics

CairoMakie.activate!()
set_aog_theme!()

#%%
includet(srcdir("Julia/transmission-functions.jl")) # Revise will keep track of file changes and reload functions as necessary
includet(srcdir("Julia/plotting-functions.jl"))
includet(srcdir("Julia/cleaning-functions.jl"))

#%%
δt = 0.005
tlower = 0.0
tmax = 250.0
tspan = (tlower, tmax)
tlength = length(tlower:δt:tmax)

#%%
@parameters β γ μ
@variables t S(t) I(t) R(t)
D = Differential(t)

#%%
R₀ = 10.0
p = Dict(γ => 1 / 8, μ => 1 / (62 * 365))
u₀ = Dict(S => 999.0, I => 1.0, R => 0.0)
push!(p, β => calculate_beta(R₀, p[γ], p[μ], 1, sum(values(u₀))))

#%%
eqs = [D(S) ~ -β * S * I + μ * (I + R), D(I) ~ β * S * I - γ * I - μ * I,
    D(R) ~ γ * I - μ * R]

#%%
@named de = ODESystem(eqs, t, [S, I, R], [β, γ, μ]; tspan = tspan)
de_simple = structural_simplify(de)
push!(
    p,
    β => calculate_beta(;
        ode = de_simple, R₀ = R₀, param = p, C = [1.0], pop_matrix = [sum(values(u₀))]
    ),
)

#%%
ode_prob = ODEProblem(de_simple, u₀, tspan, p)
ode_sol = solve(ode_prob, Tsit5(); saveat = δt)
ode_sol_df = create_sir_df(ode_sol)

#%%
colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"]

create_sir_plot(ode_sol_df; colors = colors)

#%%
function calculateR0(
    ; ode::A, nic::B, nac::B, param::Dict{Num,C}, S⁺::C
) where {A<:ODESystem,B<:Int,C<:AbstractFloat}
    Jac = calculate_jacobian(ode)[(nac + 1):(nac + nic * nac),
        (nac + 1):(nac + nic * nac)]

    F = substitute(Jac, Dict(γ => 0.0, μ => 0.0))
    V = substitute(Jac, Dict(β => 0.0))
    FV⁻¹ = F * -inv(V)
    all_eigenvals =
        convert.(
            Float64, Symbolics.value.(eigvals(eigen(substitute(FV⁻¹, Dict(S => S⁺, p...)))))
        )
    R0 = maximum(real(all_eigenvals))

    return R0
end

#%%
calculateR0(; ode = de_simple, nic = 1, nac = 1, param = p, S⁺ = 1000.0)


#%%
function calculate_beta(
    ; ode::ODESystem, nic::Int, nac::Int, R₀::T, param::Dict{Num,T}, C::Array{T},
    pop_matrix::Array{T},
) where {T<:AbstractFloat}
    size(C, 1) == size(C, 2) ? nothing : error("𝐂 must be square")
    if size(C, 1) == size(pop_matrix, 1)
        nothing
    else
        error("C and pop_matrix must have the same number of rows")
    end

    Jac = calculate_jacobian(ode)[(nac + 1):(nac + nic * nac),
        (nac + 1):(nac + nic * nac)]

    F = C .* pop_matrix
    # F = substitute(Jac, Dict(γ => 0.0, μ => 0.0))
    V = substitute(Jac, Dict(β => 0.0))
    FV⁻¹ = F * -inv(V)
    eigenvals =
        convert.(Float64, Symbolics.value.(eigvals(eigen(substitute(FV⁻¹, Dict(p...))))))
    beta = R0 / maximum(real(Symbolics.value.(eigenvals)))

    return beta
end

#%%
nic = 1 # number of infective compartments
nac = 1 # number of age classes

#%%
calculate_beta(;
    ode = de_simple, nic = 1, nac = 1, R₀ = R₀, param = p, C = [1.0],
    pop_matrix = [sum(values(u₀))],
)