"""
    calculate_beta(R₀, γ, μ, C, pop_matrix)

Calculate the value β for a given set of parameters and contact matrix.

```jldoctest
julia> calculate_beta(2.0, 1 / 8, 0.0, ones(1, 1), [1000])
0.00025
```

* * *

**TODO** Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.

* * *
"""
function calculate_beta(
    R₀::T, γ::T, μ::T, C::Array{T}, pop_matrix::Array{T}
) where {T<:AbstractFloat}
    size(C, 1) == size(C, 2) ? nothing : error("C must be square")
    if size(C, 1) == size(pop_matrix, 1)
        nothing
    else
        error("C and pop_matrix must have the same number of rows")
    end

    F = C .* pop_matrix
    V = Diagonal(repeat([γ + μ], size(C, 1)))

    FV⁻¹ = F * inv(V)
    eigenvals, eigenvectors = eigen(FV⁻¹)
    β = R₀ / maximum(real(eigenvals))

    return β
end

function calculate_beta(R₀, γ, μ, C, pop_matrix)
    return calculate_beta(
        convert(Float64, R₀),
        convert(Float64, γ),
        convert(Float64, μ),
        convert(Array{Float64}, [C]),
        convert(Array{Float64}, [pop_matrix]),
    )
end

function calculate_beta(
    ode::S, nic::T, nac::T, R₀::U, param::Dict{Num,U}, C::Array{U},
    pop_matrix::Array{U}
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
    # F = substitute(Jac, Dict(γ => 0.0, μ => 0.0))
    V = substitute(Jac, Dict(β => 0.0))
    FV⁻¹ = F * -inv(V)
    eigenvals =
        convert.(
            Float64, Symbolics.value.(eigvals(eigen(substitute(FV⁻¹, Dict(param...)))))
        )
    beta = R₀ / maximum(real(Symbolics.value.(eigenvals)))

    return beta
end

function calculate_beta(
    ode::S, nic::T, nac::T, R₀::U, param::Vector{Pair{Num,U}}, C::Array{U},
    pop_matrix::Array{U},
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    return calculate_beta(
        ode, nic, nac, R₀, Dict(param), C, pop_matrix
    )
end

"""
    calculateR0(β, γ, μ, C, pop_matrix)

Calculate the basic reproduction number R₀ for a given set of parameters and contact matrix.

```jldoctest
julia> calculateR0(0.00025, 1 / 8, 0.0, ones(1, 1), [1000])
2.0
```

* * *

**TODO** Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.

* * *
"""
function calculateR0(
    β::T, γ::T, μ::T, C::Array{T}, pop_matrix::Array{T}
) where {T<:AbstractFloat}
    size(C, 1) == size(C, 2) ? nothing : error("C must be square")
    if size(C, 1) == size(pop_matrix, 1)
        nothing
    else
        error("C and pop_matrix must have the same number of rows")
    end

    B = β * C

    F = B .* pop_matrix
    V = Diagonal(repeat([γ + μ], size(C, 1)))

    FV⁻¹ = F * inv(V)
    eigenvals, eigenvectors = eigen(FV⁻¹)

    R₀ = maximum(real(eigenvals))

    return R₀
end

function calculateR0(β, γ, μ, C, pop_matrix)
    return calculateR0(
        convert(Float64, β),
        convert(Float64, γ),
        convert(Float64, μ),
        convert(Array{Float64}, [C]),
        convert(Array{Float64}, [pop_matrix]),
    )
end

function calculateR0(
    ode::S, nic::T, nac::T, param::Dict{Num,U}, S⁺::U
) where {S<:ODESystem,T<:Int,U<:AbstractFloat}
    Jac = calculate_jacobian(ode)[(nac + 1):(nac + nic * nac),
        (nac + 1):(nac + nic * nac)]

    F = substitute(Jac, Dict(γ => 0.0, μ => 0.0))
    V = substitute(Jac, Dict(β => 0.0))
    FV⁻¹ = F * -inv(V)
    all_eigenvals =
        convert.(
            Float64,
            Symbolics.value.(eigvals(eigen(substitute(FV⁻¹, Dict(S => S⁺, param...))))),
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