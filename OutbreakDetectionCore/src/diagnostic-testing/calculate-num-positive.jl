export calculate_positives_vec!

"""
    calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64
    )

Calculate test-positive individuals accounting for test result lag.

Applies the test performance multiplier (sensitivity or 1-specificity) to the number
tested, with results appearing after the specified lag period.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of individuals tested
- `sim_length`: Length of simulation
- `lag`: Test result lag in days
- `tested_multiplier`: Test performance multiplier (sensitivity or 1-specificity)
"""
function calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64;
        rng::Random.AbstractRNG = Random.default_rng()
    )
    # Initialize to zero
    fill!(npos_vec, 0)

    @inbounds for day in eachindex(npos_vec)
        result_day = day + lag
        if result_day <= sim_length
            npos_vec[result_day] = rand(rng, Distributions.Binomial(tested_vec[day], tested_multiplier))
        end
    end

    return nothing
end
