"""
    calculate_beta(R₀, γ, μ, 𝐂, pop_matrix)

Calculate the value β for a given set of parameters and contact matrix.

```jldoctest
julia> calculate_beta(2.0, 1/8, 0.0, ones(1, 1), [1000])
0.00025
```

---

**TODO** Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.

---

"""
function calculate_beta(
        R₀::T,
        γ::T,
        μ::T,
        𝐂::Array{T},
        pop_matrix::Array{T}
    ) where {T<:AbstractFloat}
    size(𝐂, 1) == size(𝐂, 2) ? nothing : error("𝐂 must be square")
    size(𝐂, 1) == size(pop_matrix, 1) ? nothing : error("𝐂 and pop_matrix must have the same number of rows")

    𝐅 = 𝐂 .* pop_matrix
    𝐕 = Diagonal(repeat([γ + μ], size(𝐂, 1)))

    𝐅𝐕⁻¹ = 𝐅 * inv(𝐕)
    eigenvals, eigenvectors = eigen(𝐅𝐕⁻¹)
    β = R₀ / maximum(real(eigenvals))
    
    return β
end

function calculate_beta(R₀, γ, μ, 𝐂, pop_matrix)
    calculate_beta(
        convert(Float64, R₀),
        convert(Float64, γ),
        convert(Float64, μ),
        convert(Array{Float64}, [𝐂]),
        convert(Array{Float64}, [pop_matrix])
    )
end


"""
    calculateR0(β, γ, μ, 𝐂, pop_matrix)

Calculate the basic reproduction number R₀ for a given set of parameters and contact matrix.

```jldoctest
julia> calculateR0(0.00025, 1/8, 0.0, ones(1, 1), [1000])
2.0
```

---

**TODO** Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.

---

"""
function calculateR0(
        β::T,
        γ::T,
        μ::T,
        𝐂::Array{T},
        pop_matrix::Array{T}
    ) where {T<:AbstractFloat}
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

function calculateR0(β, γ, μ, 𝐂, pop_matrix)
    calculateR0(
        convert(Float64, β),
        convert(Float64, γ),
        convert(Float64, μ),
        convert(Array{Float64}, [𝐂]),
        convert(Array{Float64}, [pop_matrix])
    )
end