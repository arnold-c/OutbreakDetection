export calculate_true_positives_vec!,
    calculate_false_positives_vec!,
    _calculate_positives_vec!

"""
    calculate_true_positives(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        test_specification::IndividualTestSpecification,
        rng::Random.AbstractRNG = Random.default_rng()
    )

Calculate true positive test results from infected individuals.

This function applies the test sensitivity to the number of infected individuals
tested, accounting for test result lag. For perfect tests (sensitivity = 1.0),
all tested infected individuals will test positive. For imperfect tests, the
number of positives is drawn from a binomial distribution.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of infected individuals tested
- `sim_length`: Length of simulation
- `test_specification`: Test specifications including sensitivity and lag
- `rng`: Random number generator (defaults to `Random.default_rng()`)

# Details
The function:
1. Checks if the test has perfect sensitivity (1.0)
2. For perfect tests, uses `perfect_test_true_positives_vec` for efficiency
3. For imperfect tests, calls `calculate_positives_vec!` with test sensitivity
4. Accounts for test result lag in both cases
5. Returns nothing (modifies `npos_vec` in-place)

# See Also
- [`calculate_false_positives`](@ref): Calculate false positive results
- [`calculate_positives_vec!`](@ref): Core calculation function
"""
function calculate_true_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        test_specification::IndividualTestSpecification;
        rng::Random.AbstractRNG = Random.default_rng()
    )

    if test_specification.sensitivity == 1.0
        perfect_test_true_positives_vec(
            npos_vec,
            tested_vec,
            sim_length,
            test_specification.test_result_lag
        )
        return nothing
    end

    _calculate_positives_vec!(
        npos_vec,
        tested_vec,
        sim_length,
        test_specification.test_result_lag,
        test_specification.sensitivity;
        rng = rng
    )
    return nothing
end

"""
    perfect_test_true_positives_vec(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64
    )

Calculate true positive test results for a perfect test (100% sensitivity).

For a perfect test, all tested infected individuals will test positive. This
function is an optimized version that avoids the binomial sampling used for
imperfect tests, since the result is deterministic.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of infected individuals tested
- `sim_length`: Length of simulation
- `lag`: Test result lag in days

# Details
The function:
1. Initializes the output vector to zero
2. For each day, shifts tested counts by the lag period
3. Stores the shifted counts in the output vector
4. Returns nothing (modifies `npos_vec` in-place)

# Performance Notes
Uses `@inbounds` for optimized array access since indices are guaranteed valid.

# See Also
- [`calculate_true_positives`](@ref): Public interface for true positive calculation
- [`calculate_positives_vec!`](@ref): General calculation for imperfect tests
"""
function perfect_test_true_positives_vec(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
    )
    # Initialize to zero
    fill!(npos_vec, 0)

    return @inbounds for day in eachindex(npos_vec)
        result_day = day + lag
        if result_day <= sim_length
            npos_vec[result_day] = tested_vec[day]
        end
    end
end

"""
    calculate_false_positives(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        test_specification::IndividualTestSpecification,
        rng::Random.AbstractRNG = Random.default_rng()
    )

Calculate false positive test results from non-infected individuals.

This function applies the false positive rate (1 - specificity) to the number
of non-infected individuals tested, accounting for test result lag. For perfect
tests (specificity = 1.0), there are no false positives. For imperfect tests,
the number of false positives is drawn from a binomial distribution.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of non-infected individuals tested
- `sim_length`: Length of simulation
- `test_specification`: Test specifications including specificity and lag
- `rng`: Random number generator (defaults to `Random.default_rng()`)

# Details
The function:
1. Checks if the test has perfect specificity (1.0)
2. For perfect tests, initializes output to zero (no false positives)
3. For imperfect tests, calls `calculate_positives_vec!` with false positive rate
4. Accounts for test result lag in the calculation
5. Returns nothing (modifies `npos_vec` in-place)

# See Also
- [`calculate_true_positives`](@ref): Calculate true positive results
- [`calculate_positives_vec!`](@ref): Core calculation function
"""
function calculate_false_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        test_specification::IndividualTestSpecification;
        rng::Random.AbstractRNG = Random.default_rng()
    )

    if test_specification.specificity == 1.0
        fill!(npos_vec, 0)
        return nothing
    end

    _calculate_positives_vec!(
        npos_vec,
        tested_vec,
        sim_length,
        test_specification.test_result_lag,
        1.0 - test_specification.specificity;
        rng = rng
    )
    return nothing
end


"""
    _calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64;
        rng::Random.AbstractRNG = Random.default_rng()
    )

Calculate test-positive individuals accounting for test result lag.

Applies the test performance multiplier (sensitivity or 1-specificity) to the number
tested, with results appearing after the specified lag period. Uses binomial sampling
to account for stochastic variation in test outcomes.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of individuals tested
- `sim_length`: Length of simulation
- `lag`: Test result lag in days
- `tested_multiplier`: Test performance multiplier (sensitivity or 1-specificity)

# Keyword Arguments
- `rng`: Random number generator (defaults to `Random.default_rng()`)

# Details
The function:
1. Initializes the output vector to zero
2. For each day, shifts tested counts by the lag period
3. Applies binomial sampling with the given multiplier to account for stochastic variation
4. Stores results in the output vector
5. Returns nothing (modifies `npos_vec` in-place)

# Performance Notes
Uses `@inbounds` for optimized array access since indices are guaranteed valid.

# See Also
- [`calculate_true_positives`](@ref): Calculate true positive results
- [`calculate_false_positives`](@ref): Calculate false positive results
"""
function _calculate_positives_vec!(
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
            npos_vec[result_day] = Int64(round(tested_vec[day] * tested_multiplier))
        end
    end

    return nothing
end
