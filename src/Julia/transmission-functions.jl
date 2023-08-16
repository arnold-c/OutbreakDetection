using DrWatson
@quickactivate "OutbreakDetection"

using DataFrames, DataFramesMeta, LinearAlgebra
using ModelingToolkit, DifferentialEquations
using FLoops


"""
    calculate_beta(R_0, gamma, mu, C, pop_matrix)

Calculate the value beta for a given set of parameters and contact matrix.

```jldoctest
julia> calculate_beta(2.0, 1 / 8, 0.0, ones(1, 1), [1000])
0.00025
```

"""
# TODO: Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.
function calculate_beta(
    R_0::T, gamma::T, mu::T, C::Array{T}, pop_matrix::Array{T}
) where {T<:AbstractFloat}
    size(C, 1) == size(C, 2) ? nothing : error("C must be square")
    if size(C, 1) == size(pop_matrix, 1)
        nothing
    else
        error("C and pop_matrix must have the same number of rows")
    end

    F = C .* pop_matrix
    V = Diagonal(repeat([gamma + mu], size(C, 1)))

    FV⁻¹ = F * inv(V)
    eigenvals, eigenvectors = eigen(FV⁻¹)
    beta = R_0 / maximum(real(eigenvals))

    return beta
end

function calculate_beta(R_0, gamma, mu, C, pop_matrix)
    return calculate_beta(
        convert(Float64, R_0),
        convert(Float64, gamma),
        convert(Float64, mu),
        convert(Array{Float64}, [C]),
        convert(Array{Float64}, [pop_matrix]),
    )
end

function calculate_beta(
    ode::S, nic::T, nac::T, R_0::U, param::Dict{Num,U}, C::Array{U},
    pop_matrix::Array{U},
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    size(C, 1) == size(C, 2) ? nothing : error("C must be square")
    if size(C, 1) == size(pop_matrix, 1)
        nothing
    else
        error("C and pop_matrix must have the same number of rows")
    end

    Jac = calculate_jacobian(ode)[(nac + 1):(nac + nic * nac),
        (nac + 1):(nac + nic * nac)]

    F = C .* pop_matrix
    # F = substitute(Jac, Dict(gamma => 0.0, mu => 0.0))
    V = substitute(Jac, Dict(beta => 0.0))
    FV⁻¹ = F * -inv(V)
    eigenvals =
        convert.(
            Float64,
            Symbolics.value.(eigvals(eigen(substitute(FV⁻¹, Dict(param...))))),
        )
    beta = R_0 / maximum(real(Symbolics.value.(eigenvals)))

    return beta
end

function calculate_beta(
    ode::S, nic::T, nac::T, R_0::U, param::Vector{Pair{Num,U}}, C::Array{U},
    pop_matrix::Array{U},
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    return calculate_beta(
        ode, nic, nac, R_0, Dict(param), C, pop_matrix
    )
end

function calculate_beta(
    ode::S, nic::T, nac::T, R_0::U, param::Dict{Num,U}, C::Array{U},
    pop_matrix::Array{T},
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    return calculate_beta(
        ode, nic, nac, R_0, Dict(param), C, convert.(Float64, pop_matrix)
    )
end

"""
    calculateR0(beta, gamma, mu, C, pop_matrix)

Calculate the basic reproduction number R_0 for a given set of parameters and contact matrix.

```jldoctest
julia> calculateR0(0.00025, 1 / 8, 0.0, ones(1, 1), [1000])
2.0
```

* * *

**TODO** Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.

* * *
"""
function calculateR0(
    beta::T, gamma::T, mu::T, C::Array{T}, pop_matrix::Array{T}
) where {T<:AbstractFloat}
    size(C, 1) == size(C, 2) ? nothing : error("C must be square")
    if size(C, 1) == size(pop_matrix, 1)
        nothing
    else
        error("C and pop_matrix must have the same number of rows")
    end

    B = beta * C

    F = B .* pop_matrix
    V = Diagonal(repeat([gamma + mu], size(C, 1)))

    FV⁻¹ = F * inv(V)
    eigenvals, eigenvectors = eigen(FV⁻¹)

    R_0 = maximum(real(eigenvals))

    return R_0
end

function calculateR0(beta, gamma, mu, C, pop_matrix)
    return calculateR0(
        convert(Float64, beta),
        convert(Float64, gamma),
        convert(Float64, mu),
        convert(Array{Float64}, [C]),
        convert(Array{Float64}, [pop_matrix]),
    )
end

function calculateR0(
    ode::S, nic::T, nac::T, param::Dict{Num,U}, S⁺::U
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    Jac = calculate_jacobian(ode)[(nac + 1):(nac + nic * nac),
        (nac + 1):(nac + nic * nac)]

    F = substitute(Jac, Dict(gamma => 0.0, mu => 0.0))
    V = substitute(Jac, Dict(beta => 0.0))
    FV⁻¹ = F * -inv(V)
    all_eigenvals =
        convert.(
            Float64,
            Symbolics.value.(
                eigvals(eigen(substitute(FV⁻¹, Dict(S => S⁺, param...))))
            ),
        )
    R0 = maximum(real(all_eigenvals))

    return R0
end

function calculateR0(
    ode::S, nic::T, nac::T, param::Vector{Pair{Num,U}}, S⁺::U
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    return ngmR0(
        ode, nic, nac, Dict(param), S⁺
    )
end
